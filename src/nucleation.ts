/**
 * When does vapor condense? Classical nucleation theory for a pure vapor (homogeneous) and
 * Thomson's theory for a droplet forming around an ion.
 *
 * A liquid embryo of radius r in vapor at supersaturation S has free energy
 *   dG(r) = -(4 pi r^3 / 3 v) kT ln S + 4 pi r^2 sigma [+ charge term for an ion]
 * The surface term makes small embryos unstable: they must climb a barrier before growth pays.
 * An ion lowers the barrier because polarizable liquid around a charge stores less field energy.
 *
 * Both theories are known to be approximate (they treat a cluster of tens of molecules as bulk
 * liquid). They give thresholds within tens of percent of experiment, which is what sets where
 * tracks can form. The functions take the liquid's properties so a test can run them for water.
 */
import { E_CHARGE, EPS0, K_B } from './constants'

export interface Liquid {
  /** Surface tension, N/m. */
  sigma: number
  /** Volume per molecule in the liquid, m^3. */
  molecularVolume: number
  /** Mass per molecule, kg. */
  molecularMass: number
  /** Relative permittivity of the liquid. */
  dielectric: number
  /** Saturation vapor pressure at this temperature, Pa. */
  saturationPa: number
}

/** Radius of the droplet in unstable equilibrium with vapor at supersaturation S (Kelvin), m. */
export function kelvinRadius(liquid: Liquid, tempK: number, s: number): number {
  return (2 * liquid.sigma * liquid.molecularVolume) / (K_B * tempK * Math.log(s))
}

/** Homogeneous nucleation rate, droplets per m^3 per second (Becker-Doering-Zeldovich form). */
export function homogeneousRate(liquid: Liquid, tempK: number, s: number): number {
  if (s <= 1) return 0
  const kt = K_B * tempK
  const { sigma, molecularVolume: v, molecularMass: m } = liquid
  const barrier = (16 * Math.PI * sigma ** 3 * v * v) / (3 * (kt * Math.log(s)) ** 2)
  const n1 = (s * liquid.saturationPa) / kt
  return Math.sqrt((2 * sigma) / (Math.PI * m)) * v * n1 * n1 * Math.exp(-barrier / kt)
}

/** Coefficient of the charge term: dG_charge(r) = chargeCoeff / r (plus a constant). */
function chargeCoeff(liquid: Liquid): number {
  return ((E_CHARGE * E_CHARGE) / (8 * Math.PI * EPS0)) * (1 - 1 / liquid.dielectric)
}

export interface IonBarrier {
  /** Height of the barrier between the stable ion cluster and the critical droplet, J. 0 if none. */
  height: number
  /** Radius of the critical droplet, m (0 when there is no barrier). */
  criticalRadius: number
  /** Curvature of dG at the critical radius, J/m^2 (negative). */
  curvature: number
}

/** Thomson barrier for a singly charged ion. */
export function ionBarrier(liquid: Liquid, tempK: number, s: number): IonBarrier {
  const none = { height: 0, criticalRadius: 0, curvature: 0 }
  if (s <= 1) return { height: Number.POSITIVE_INFINITY, criticalRadius: 0, curvature: 0 }
  const kt = K_B * tempK
  const a = (4 * Math.PI * kt * Math.log(s)) / liquid.molecularVolume // dG = -(a/3) r^3 + ...
  const b = 8 * Math.PI * liquid.sigma
  const c = chargeCoeff(liquid)
  // dG'(r) = -a r^2 + b r - c / r^2. Multiply by r^2: g(r) = -a r^4 + b r^3 - c.
  const g = (r: number) => -a * r ** 4 + b * r ** 3 - c
  const rPeak = (3 * b) / (4 * a)
  if (g(rPeak) <= 0) return none // dG falls monotonically: every ion grows
  const bisect = (lo: number, hi: number) => {
    let l = lo
    let h = hi
    const rising = g(lo) < g(hi)
    for (let i = 0; i < 200; i++) {
      const mid = 0.5 * (l + h)
      if (g(mid) < 0 === rising) l = mid
      else h = mid
    }
    return 0.5 * (l + h)
  }
  const rStable = bisect(1e-12, rPeak)
  const rCritical = bisect(rPeak, b / a) // g(b/a) = -c < 0
  const dG = (r: number) => (-a / 3) * r ** 3 + (b / 2) * r * r + c / r
  return {
    height: dG(rCritical) - dG(rStable),
    criticalRadius: rCritical,
    curvature: -2 * a * rCritical + b + (2 * c) / rCritical ** 3,
  }
}

/** Supersaturation above which an ion grows with no barrier at all (Thomson's limit). */
export function ionBarrierlessSaturation(liquid: Liquid, tempK: number): number {
  // The barrier vanishes where the two roots of dG' merge: r^3 = c / (2 pi sigma) and
  // kT ln S = (v / 4 pi r^2) (8 pi sigma r - c / r^2) there, which reduces to 3 sigma v / (2 r).
  const r = Math.cbrt(chargeCoeff(liquid) / (2 * Math.PI * liquid.sigma))
  return Math.exp((3 * liquid.sigma * liquid.molecularVolume) / (2 * r * K_B * tempK))
}

/**
 * Rate at which one ion becomes a growing droplet, per second: molecules striking the critical
 * droplet, times the Zeldovich factor, times the Boltzmann factor of the barrier.
 */
export function ionNucleationRate(liquid: Liquid, tempK: number, s: number): number {
  if (s <= 1) return 0
  const kt = K_B * tempK
  const pv = s * liquid.saturationPa
  const flux = pv / Math.sqrt(2 * Math.PI * liquid.molecularMass * kt) // molecules / m^2 / s
  const barrier = ionBarrier(liquid, tempK, s)
  if (barrier.height === 0) {
    // No barrier: growth is limited only by arrival of molecules on a ~1 nm cluster.
    const r = Math.cbrt(chargeCoeff(liquid) / (2 * Math.PI * liquid.sigma))
    return flux * 4 * Math.PI * r * r
  }
  const r = barrier.criticalRadius
  const area = 4 * Math.PI * r * r
  const drdn = liquid.molecularVolume / area
  const zeldovich = Math.sqrt((-barrier.curvature * drdn * drdn) / (2 * Math.PI * kt))
  return flux * area * zeldovich * Math.exp(-barrier.height / kt)
}

/** Smallest S at which `rate(S)` reaches `target`, by bisection on an increasing function. */
export function thresholdSaturation(rate: (s: number) => number, target: number): number {
  let lo = 1.0001
  let hi = 60
  if (rate(hi) < target) return hi
  for (let i = 0; i < 80; i++) {
    const mid = 0.5 * (lo + hi)
    if (rate(mid) < target) lo = mid
    else hi = mid
  }
  return 0.5 * (lo + hi)
}
