/**
 * The chamber at rest: steady temperature and vapor profiles between a cold plate (z = 0) and a
 * warm, alcohol-soaked top (z = H), in still air.
 *
 * Vapor evaporates at the top, diffuses down and condenses on the plate. Its partial pressure
 * falls roughly linearly with height while the saturation pressure falls exponentially with
 * temperature, so the vapor ends up supersaturated in between. Nothing here places the
 * sensitive layer; it is wherever the solved supersaturation exceeds the ion threshold.
 *
 * Model, in one dimension:
 *  - Vapor: Stefan diffusion through stagnant air. With u = -ln(1 - x) for mole fraction x, the
 *    molar flux is N = -c D du/dz, so u is linear in the "resistance" coordinate
 *    zeta = integral dz / (c D) wherever nothing condenses.
 *  - Fog limit: supersaturation cannot exceed the homogeneous-nucleation threshold S_hom(T).
 *    Past it the vapor makes its own droplets, which remove the excess. So u solves an obstacle
 *    problem: u'' >= 0, u <= u_max, and u'' = 0 off the obstacle. The solution is the lower
 *    convex hull of the obstacle with the two boundary values pinned.
 *  - Heat: conduction through air plus the latent heat released where vapor condenses.
 *  - Mist: the condensation sink is carried by droplets that nucleate, grow and fall. Their
 *    nucleation rate is solved from the sink, top down.
 */
import { airConductivity } from './air'
import { ATM, K_B, R_GAS } from './constants'
import { growthCoefficient, stokesCoefficient } from './droplet'
import {
  IPA,
  IPA_MOLECULE_KG,
  liquidDensity,
  molecularVolume,
  saturationPressure,
  surfaceTension,
  vaporDiffusivity,
  vaporizationEnthalpy,
} from './ipa'
import { homogeneousRate, ionNucleationRate, type Liquid, thresholdSaturation } from './nucleation'

/** Homogeneous nucleation counts as "on" at one droplet per cm^3 per second. */
export const HOMOGENEOUS_ONSET = 1e6 // 1/(m^3 s)
/** An ion counts as nucleating when its rate reaches one per second. */
export const ION_ONSET = 1 // 1/s

export function ipaLiquid(tempK: number): Liquid {
  return {
    sigma: surfaceTension(tempK),
    molecularVolume: molecularVolume(tempK),
    molecularMass: IPA_MOLECULE_KG,
    dielectric: IPA.dielectric,
    saturationPa: saturationPressure(tempK),
  }
}

// Thresholds are smooth in T; tabulate once on a 1 K grid.
const TH_T0 = 180
const TH_N = 161
let thresholds: { ion: Float64Array; hom: Float64Array } | undefined
function thresholdTables() {
  if (thresholds) return thresholds
  const ion = new Float64Array(TH_N)
  const hom = new Float64Array(TH_N)
  for (let i = 0; i < TH_N; i++) {
    const t = TH_T0 + i
    const liquid = ipaLiquid(t)
    ion[i] = thresholdSaturation((s) => ionNucleationRate(liquid, t, s), ION_ONSET)
    hom[i] = thresholdSaturation((s) => homogeneousRate(liquid, t, s), HOMOGENEOUS_ONSET)
  }
  thresholds = { ion, hom }
  return thresholds
}
function lerpTable(table: Float64Array, tempK: number): number {
  const x = Math.min(Math.max(tempK - TH_T0, 0), TH_N - 1.000001)
  const i = Math.floor(x)
  return table[i] + (x - i) * (table[i + 1] - table[i])
}
/** Supersaturation at which an ion nucleates within about a second. */
export const ionThreshold = (tempK: number) => lerpTable(thresholdTables().ion, tempK)
/** Supersaturation at which the vapor fogs on its own. */
export const fogThreshold = (tempK: number) => lerpTable(thresholdTables().hom, tempK)

export interface ChamberConfig {
  /** Distance from plate to top, m. */
  height: number
  plateK: number
  topK: number
  pressurePa: number
  /** Grid intervals. */
  cells: number
}

export const DEFAULT_CHAMBER: ChamberConfig = {
  height: 0.12,
  plateK: 233.15, // -40 C: a plate on dry ice, with losses
  topK: 293.15,
  pressurePa: ATM,
  cells: 480,
}

export interface ChamberProfile {
  config: ChamberConfig
  /** Node heights, m (cells + 1 values). */
  z: Float64Array
  tempK: Float64Array
  /** Vapor partial pressure, Pa. */
  vaporPa: Float64Array
  saturation: Float64Array
  ionThreshold: Float64Array
  fogThreshold: Float64Array
  /** Vapor removed by condensation in the volume, mol/(m^3 s). */
  sink: Float64Array
  /** Mist droplets nucleated, 1/(m^3 s). */
  mistRate: Float64Array
  /** Downward vapor flux leaving the top and arriving at the plate, mol/(m^2 s). */
  fluxTop: number
  fluxPlate: number
  /** Bottom and top of the layer where S >= the ion threshold, m. Both 0 if there is none. */
  sensitiveBottom: number
  sensitiveTop: number
  /** Fraction of the condensation sink that the solved mist population accounts for. */
  mistClosure: number
  iterations: number
}

/** Lower convex hull of (x_i, y_i) with x increasing; returns y on the hull at every x_i. */
export function lowerConvexHull(x: Float64Array, y: Float64Array): Float64Array {
  const hull: number[] = []
  for (let i = 0; i < x.length; i++) {
    while (hull.length >= 2) {
      const a = hull[hull.length - 2]
      const b = hull[hull.length - 1]
      // Drop b if it lies on or above the chord a -> i.
      const cross = (x[b] - x[a]) * (y[i] - y[a]) - (y[b] - y[a]) * (x[i] - x[a])
      if (cross <= 0) hull.pop()
      else break
    }
    hull.push(i)
  }
  const out = new Float64Array(x.length)
  for (let h = 0; h + 1 < hull.length; h++) {
    const a = hull[h]
    const b = hull[h + 1]
    for (let i = a; i <= b; i++) out[i] = y[a] + ((y[b] - y[a]) * (x[i] - x[a])) / (x[b] - x[a])
  }
  return out
}

/** Solve a tridiagonal system in place (Thomas algorithm). */
function solveTridiagonal(
  lower: Float64Array,
  diag: Float64Array,
  upper: Float64Array,
  rhs: Float64Array,
) {
  const n = diag.length
  for (let i = 1; i < n; i++) {
    const m = lower[i] / diag[i - 1]
    diag[i] -= m * upper[i - 1]
    rhs[i] -= m * rhs[i - 1]
  }
  rhs[n - 1] /= diag[n - 1]
  for (let i = n - 2; i >= 0; i--) rhs[i] = (rhs[i] - upper[i] * rhs[i + 1]) / diag[i]
  return rhs
}

export function solveChamber(config: ChamberConfig = DEFAULT_CHAMBER): ChamberProfile {
  const n = config.cells + 1
  const dz = config.height / config.cells
  const p = config.pressurePa
  const z = Float64Array.from({ length: n }, (_, i) => i * dz)
  const temp = Float64Array.from(
    z,
    (h) => config.plateK + ((config.topK - config.plateK) * h) / config.height,
  )
  const u = new Float64Array(n)
  const sink = new Float64Array(n)
  const zeta = new Float64Array(n)
  const obstacle = new Float64Array(n)
  let iterations = 0

  const solveVapor = () => {
    // Resistance coordinate and the fog-limit obstacle for the current temperatures.
    for (let i = 0; i < n; i++) {
      const sat = saturationPressure(temp[i])
      const limit = i === 0 || i === n - 1 ? sat : Math.min(fogThreshold(temp[i]) * sat, 0.9 * p)
      obstacle[i] = -Math.log(1 - limit / p)
      if (i > 0) {
        const tm = 0.5 * (temp[i] + temp[i - 1])
        zeta[i] = zeta[i - 1] + dz / ((p / (R_GAS * tm)) * vaporDiffusivity(tm, p))
      }
    }
    u.set(lowerConvexHull(zeta, obstacle))
    // Sink = divergence of the flux N = -du/dzeta (negative: vapor moves down).
    for (let i = 1; i < n - 1; i++) {
      const above = -(u[i + 1] - u[i]) / (zeta[i + 1] - zeta[i])
      const below = -(u[i] - u[i - 1]) / (zeta[i] - zeta[i - 1])
      // Vapor condenses only where the profile rides the fog limit; elsewhere the flux is
      // constant and any difference is rounding.
      const onLimit = obstacle[i] - u[i] <= 1e-12 * obstacle[i]
      sink[i] = onLimit ? Math.max(0, -(above - below) / dz) : 0
    }
  }

  const solveHeat = () => {
    const lower = new Float64Array(n)
    const diag = new Float64Array(n)
    const upper = new Float64Array(n)
    const rhs = new Float64Array(n)
    diag[0] = 1
    rhs[0] = config.plateK
    diag[n - 1] = 1
    rhs[n - 1] = config.topK
    for (let i = 1; i < n - 1; i++) {
      const kUp = airConductivity(0.5 * (temp[i] + temp[i + 1]))
      const kDown = airConductivity(0.5 * (temp[i] + temp[i - 1]))
      lower[i] = kDown
      upper[i] = kUp
      diag[i] = -(kUp + kDown)
      rhs[i] = -vaporizationEnthalpy(temp[i]) * sink[i] * dz * dz
    }
    return solveTridiagonal(lower, diag, upper, rhs)
  }

  for (iterations = 1; iterations <= 200; iterations++) {
    solveVapor()
    const next = solveHeat()
    let change = 0
    for (let i = 0; i < n; i++) {
      change = Math.max(change, Math.abs(next[i] - temp[i]))
      temp[i] += 0.5 * (next[i] - temp[i])
    }
    if (change < 1e-5) break
  }
  solveVapor()

  const vaporPa = new Float64Array(n)
  const saturation = new Float64Array(n)
  const ionTh = new Float64Array(n)
  const fogTh = new Float64Array(n)
  for (let i = 0; i < n; i++) {
    vaporPa[i] = p * (1 - Math.exp(-u[i]))
    saturation[i] = vaporPa[i] / saturationPressure(temp[i])
    ionTh[i] = ionThreshold(temp[i])
    fogTh[i] = fogThreshold(temp[i])
  }

  let sensitiveBottom = 0
  let sensitiveTop = 0
  for (let i = 1; i < n - 1; i++) {
    if (saturation[i] >= ionTh[i]) {
      if (sensitiveTop === 0) sensitiveBottom = z[i]
      sensitiveTop = z[i]
    }
  }

  const { mistRate, closure } = solveMist(z, temp, saturation, sink, p)

  return {
    config,
    z,
    tempK: temp,
    vaporPa,
    saturation,
    ionThreshold: ionTh,
    fogThreshold: fogTh,
    sink,
    mistRate,
    fluxTop: (u[n - 1] - u[n - 2]) / (zeta[n - 1] - zeta[n - 2]),
    fluxPlate: (u[1] - u[0]) / (zeta[1] - zeta[0]),
    sensitiveBottom,
    sensitiveTop,
    mistClosure: closure,
    iterations,
  }
}

/**
 * Which droplet population carries the sink? A droplet born at height z0 falls at c r^2 while
 * growing as r dr/dt = g (S - 1), so d(r^4)/dz = -4 g (S - 1) / c and
 *   r^4(z) = 4 [Gamma(z0) - Gamma(z)],  Gamma(z) = integral_0^z g (S - 1) / c dz.
 * The liquid a droplet gains crossing a cell is then known, and the nucleation rate in each cell
 * follows from the sink there minus what droplets from above already take: a triangular system
 * solved from the top cell down, with rates clamped at zero.
 */
function solveMist(
  z: Float64Array,
  temp: Float64Array,
  saturation: Float64Array,
  sink: Float64Array,
  p: number,
) {
  const n = z.length
  const dz = z[1] - z[0]
  const gamma = new Float64Array(n)
  for (let i = 1; i < n; i++) {
    const tm = 0.5 * (temp[i] + temp[i - 1])
    const sm = 0.5 * (saturation[i] + saturation[i - 1])
    gamma[i] =
      gamma[i - 1] + (dz * growthCoefficient(tm, p) * Math.max(sm - 1, 0)) / stokesCoefficient(tm)
  }
  // Cell i spans z[i] - dz/2 .. z[i] + dz/2; Gamma at its faces by interpolation.
  const face = (i: number) => 0.5 * (gamma[Math.max(i, 0)] + gamma[Math.min(i + 1, n - 1)]) // top face of cell i
  const mass = (r4: number, tempK: number) =>
    (4 / 3) * Math.PI * liquidDensity(tempK) * Math.max(r4, 0) ** 0.75
  const mistRate = new Float64Array(n)
  // uptake[i]: kg/(m^2 s) absorbed in cell i by droplets born above it.
  const uptake = new Float64Array(n)
  let demand = 0
  let met = 0
  for (let j = n - 2; j >= 1; j--) {
    const want = sink[j] * IPA.molarMass * dz // kg/(m^2 s) to remove in this cell
    demand += want
    const own = mass(4 * (gamma[j] - face(j - 1)), temp[j]) // gain of a droplet born mid-cell
    const rate = own > 0 ? Math.max(0, (want - uptake[j]) / own) : 0 // droplets/(m^2 s) born here
    mistRate[j] = rate / dz
    met += Math.min(want, uptake[j] + rate * own)
    if (rate > 0) {
      for (let i = j - 1; i >= 1; i--) {
        uptake[i] +=
          rate *
          (mass(4 * (gamma[j] - face(i - 1)), temp[i]) - mass(4 * (gamma[j] - face(i)), temp[i]))
      }
    }
  }
  return { mistRate, closure: demand > 0 ? met / demand : 1 }
}

/** Linear interpolation of a profile array at height zq. */
export function sampleProfile(profile: ChamberProfile, values: Float64Array, zq: number): number {
  const dz = profile.z[1] - profile.z[0]
  const x = Math.min(Math.max(zq / dz, 0), profile.z.length - 1.000001)
  const i = Math.floor(x)
  return values[i] + (x - i) * (values[i + 1] - values[i])
}

/** Molecules of vapor per m^3 at a node. */
export function vaporDensityAt(profile: ChamberProfile, i: number): number {
  return profile.vaporPa[i] / (K_B * profile.tempK[i])
}
