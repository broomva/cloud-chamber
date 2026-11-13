/**
 * The running chamber: radiation arrives, ions nucleate droplets, droplets grow and fall.
 *
 * Droplets live in flat arrays so the renderer can upload them directly. One entry is a
 * "super-droplet" standing for `weight` real droplets that share a position and radius. Sparse
 * tracks (muons, electrons) use weight 1, one entry per real droplet; an alpha track has millions
 * of ions per metre and would not fit, so its entries stand for a few dozen droplets each.
 */
import { airDensity, airMeanFreePath, airViscosity } from './air'
import {
  type ChamberConfig,
  type ChamberProfile,
  DEFAULT_CHAMBER,
  ipaLiquid,
  solveChamber,
} from './chamber'
import { E_CHARGE, K_B } from './constants'
import { growthCoefficient, stokesCoefficient } from './droplet'
import { vaporDiffusivity } from './ipa'
import { ionNucleationRate } from './nucleation'
import { Rng } from './rng'
import {
  ambientElectrons,
  type Box,
  cosmicMuons,
  enterBox,
  radonAlphas,
  type Source,
} from './sources'
import {
  DEFAULT_TRANSPORT,
  type Particle,
  type Segment,
  type Species,
  type TransportOptions,
  transport,
} from './transport'

/** Mobility of small ions in air near 1 atm, m^2/(V s). Positive and negative ions differ by ~30%. */
export const ION_MOBILITY = 1.6e-4
/** How long an ion stays available to nucleate before it is lost, s. */
export const ION_LIFETIME = 1
/** Radius of a fresh ion column, m. Sets when neighbouring droplets start to share vapor. */
export const TRACK_RADIUS = 1e-4
/**
 * Offset in the downdraft of a thin column of droplets: a column of length L and radius R that
 * carries weight f per unit length sinks at (f / 4 pi eta) (ln(L / R) + COLUMN_OFFSET) in slow
 * viscous flow. The value is the mean, over droplets filling a cylinder uniformly, of the flow
 * all the others induce; a test recomputes it by direct summation. (A rigid rod has 0.84.)
 */
export const COLUMN_OFFSET = 0.45
/** Euler's constant, in Lamb's drag law. */
const EULER_GAMMA = 0.5772156649

/**
 * How much slower a droplet grows because its neighbours on the track draw on the same vapor.
 * Treating the track as a line sink in quasi-steady radial diffusion, growth is divided by
 * 1 + 2 n r L, with n droplets per metre of radius r and L = ln(R_far / R_track) growing as the
 * depleted zone spreads: L = (1/2) ln(1 + 4 D t / R_track^2). Heat and vapor are taken to spread
 * at the same rate, so the latent-heat factor in the growth law cancels.
 */
export function vaporShare(
  perMetre: number,
  radius: number,
  diffusivity: number,
  age: number,
): number {
  if (perMetre <= 0) return 1
  const reach = 0.5 * Math.log(1 + (4 * diffusivity * age) / (TRACK_RADIUS * TRACK_RADIUS))
  return 1 + 2 * perMetre * radius * reach
}

/**
 * The factor K in the downdraft of a track, v = (f / 4 pi eta) K, for weight per length f.
 * Slow viscous flow gives K = ln(L / R) + COLUMN_OFFSET (`slender`). That logarithm counts air
 * dragged along out to a distance of order L, but inertia stops the flow following the track
 * beyond about nu / v. Lamb's solution of Oseen's equations for a long cylinder gives
 * K = 1/2 - gamma - ln(v R / 4 nu), and v itself is `unit` K, so K solves
 *   K + ln K = 1/2 - gamma - ln(unit R / 4 nu).
 * The left side rises with K, so there is one root; the smaller of it and `slender` applies.
 * `unit` is the speed at K = 1, m/s, and `kinematic` the kinematic viscosity, m^2/s.
 */
export function downdraftFactor(slender: number, unit: number, kinematic: number): number {
  if (slender <= 0 || unit <= 0) return 0
  const target = 0.5 - EULER_GAMMA - Math.log((unit * TRACK_RADIUS) / (4 * kinematic))
  if (slender + Math.log(slender) <= target) return slender
  // Newton on g(K) = K + ln K - target, which is increasing and concave: from a point at or
  // below the root it climbs to it without overshooting.
  let k = Math.min(slender, Math.exp(Math.min(target, 0)))
  for (let i = 0; i < 6; i++) k -= (k + Math.log(k) - target) / (1 + 1 / k)
  return Math.min(k, slender)
}

/** Most super-droplets created per metre of track. */
export const MAX_LINE_DENSITY = 2.5e5

export interface EngineOptions {
  seed: number
  chamber: ChamberConfig
  box: { halfX: number; halfY: number }
  /** Capacity of the droplet arrays. */
  capacity: number
  /** Background mist from the fog zone. */
  mist: boolean
  /** Use the built-in background sources (muons, radon, gamma electrons). */
  background: boolean
  /** Let a track's droplets sink together (off only to show what the term does). */
  collectiveSettling: boolean
  transport: Partial<TransportOptions>
}

export const DEFAULT_ENGINE: EngineOptions = {
  seed: 1,
  chamber: DEFAULT_CHAMBER,
  box: { halfX: 0.1, halfY: 0.1 },
  capacity: 400_000,
  mist: true,
  background: true,
  collectiveSettling: true,
  transport: {},
}

export interface EventRecord {
  time: number
  source: string
  species: Species
  kineticMeV: number
  /** Super-droplets created by this event. */
  droplets: number
}

/** Kind tags stored per droplet, for colouring and for tests. */
export const KIND = { mist: 0, alpha: 1, electron: 2, muon: 3 } as const

export class Engine {
  readonly profile: ChamberProfile
  readonly box: Box
  readonly sources: Source[] = []
  readonly rng: Rng
  readonly options: EngineOptions
  readonly transportOptions: TransportOptions

  time = 0
  count = 0
  /** Droplets lost because the arrays were full. A healthy run keeps this at 0. */
  dropped = 0
  /** Droplets that reached the plate. */
  landed = 0
  /** Droplets created so far (real droplets, not array entries). */
  created = 0
  /** Ion pairs made in the gas so far, whether or not they went on to nucleate. */
  ionPairs = 0
  readonly events: EventRecord[] = []

  // Per-droplet state.
  readonly pos: Float32Array // x, y, z (m)
  readonly radius: Float32Array // m
  readonly weight: Float32Array // real droplets represented
  readonly lineDensity: Float32Array // real droplets per metre of the parent track
  readonly birth: Float32Array // time the droplet appears (after the ion's nucleation delay)
  readonly settle: Float32Array // ln(L / R) + COLUMN_OFFSET of the parent track, 0 for a lone droplet
  readonly kind: Uint8Array

  // Per-node tables derived from the profile.
  private readonly growth: Float64Array // g (S - 1), m^2/s
  private readonly stokes: Float64Array // 1/(m s)
  private readonly meanFreePath: Float64Array
  private readonly ionRate: Float64Array // 1/s
  private readonly ionDiffusivity: Float64Array
  private readonly vaporDiff: Float64Array
  private readonly kinematic: Float64Array // kinematic viscosity of the gas, m^2/s
  private readonly dz: number
  private readonly mistTotal: number // droplets per second over the whole box
  private readonly mistCdf: Float64Array
  private mistWeight = 1

  constructor(options: Partial<EngineOptions> = {}) {
    this.options = { ...DEFAULT_ENGINE, ...options }
    const o = this.options
    this.rng = new Rng(o.seed)
    this.profile = solveChamber(o.chamber)
    this.box = { halfX: o.box.halfX, halfY: o.box.halfY, height: o.chamber.height }
    this.transportOptions = {
      ...DEFAULT_TRANSPORT,
      halfX: o.box.halfX,
      halfY: o.box.halfY,
      ...o.transport,
    }
    const n = this.profile.z.length
    this.dz = this.profile.z[1] - this.profile.z[0]
    this.growth = new Float64Array(n)
    this.stokes = new Float64Array(n)
    this.meanFreePath = new Float64Array(n)
    this.ionRate = new Float64Array(n)
    this.ionDiffusivity = new Float64Array(n)
    this.vaporDiff = new Float64Array(n)
    this.kinematic = new Float64Array(n)
    this.mistCdf = new Float64Array(n)
    const p = o.chamber.pressurePa
    let mist = 0
    for (let i = 0; i < n; i++) {
      const t = this.profile.tempK[i]
      const s = this.profile.saturation[i]
      this.growth[i] = growthCoefficient(t, p) * (s - 1)
      this.stokes[i] = stokesCoefficient(t)
      this.meanFreePath[i] = airMeanFreePath(t, p)
      this.ionRate[i] = ionNucleationRate(ipaLiquid(t), t, s)
      this.ionDiffusivity[i] = (ION_MOBILITY * K_B * t) / E_CHARGE
      this.vaporDiff[i] = vaporDiffusivity(t, p)
      this.kinematic[i] = airViscosity(t) / (airDensity(t, p) * 1e3)
      mist += this.profile.mistRate[i] * this.dz
      this.mistCdf[i] = mist
    }
    this.mistTotal = mist * 4 * o.box.halfX * o.box.halfY

    this.pos = new Float32Array(3 * o.capacity)
    this.radius = new Float32Array(o.capacity)
    this.weight = new Float32Array(o.capacity)
    this.lineDensity = new Float32Array(o.capacity)
    this.birth = new Float32Array(o.capacity)
    this.settle = new Float32Array(o.capacity)
    this.kind = new Uint8Array(o.capacity)

    if (o.background) {
      const mid = 0.5 * (o.chamber.plateK + o.chamber.topK)
      this.sources.push(
        cosmicMuons(this.box),
        radonAlphas(this.box),
        ambientElectrons(this.box, mid, p),
      )
    }
    // Keep the mist within a third of the arrays by letting one entry stand for several droplets.
    const steadyMist = this.mistTotal * 4 // droplets alive if each lasts ~4 s
    this.mistWeight = Math.max(1, Math.ceil(steadyMist / (o.capacity / 3)))
  }

  /** Chance that an ion made at height z seeds a droplet before it is lost. */
  nucleationChance(z: number): number {
    const { i, f } = this.node(z)
    return (
      1 - Math.exp(-(this.ionRate[i] + f * (this.ionRate[i + 1] - this.ionRate[i])) * ION_LIFETIME)
    )
  }

  private node(z: number): { i: number; f: number } {
    const x = Math.min(Math.max(z / this.dz, 0), this.profile.z.length - 1.000001)
    const i = Math.floor(x)
    return { i, f: x - i }
  }

  /**
   * Add one super-droplet that starts growing at time `birth`. `lineDensity` is the number of
   * real droplets per metre sharing its vapor (0 for a droplet on its own).
   */
  add(
    x: number,
    y: number,
    z: number,
    weight: number,
    lineDensity: number,
    birth: number,
    kind: number,
  ): boolean {
    if (this.count >= this.options.capacity) {
      this.dropped += weight
      return false
    }
    const k = this.count++
    this.created += weight
    this.pos[3 * k] = x
    this.pos[3 * k + 1] = y
    this.pos[3 * k + 2] = z
    this.radius[k] = 0
    this.weight[k] = weight
    this.lineDensity[k] = lineDensity
    this.birth[k] = birth
    this.settle[k] = 0
    this.kind[k] = kind
    return true
  }

  /** Turn one track segment into droplets-to-be. Returns the number of entries added. */
  private deposit(seg: Segment): number {
    this.ionPairs += seg.ionPairs
    const ions = 2 * seg.ionPairs // both ions of a pair can seed a droplet
    if (ions <= 0) return 0
    const mid = seg.start[2] + 0.5 * seg.dir[2] * seg.length
    const { i, f } = this.node(mid)
    const rate = this.ionRate[i] + f * (this.ionRate[i + 1] - this.ionRate[i])
    // The threshold is not sharp: an ion nucleates if its turn comes before it is lost.
    const chance = 1 - Math.exp(-rate * ION_LIFETIME)
    if (chance < 1e-6) return 0
    const droplets = ions * chance
    const perMetre = seg.length > 0 ? droplets / seg.length : droplets / TRACK_RADIUS
    const weight = Math.max(1, Math.round(perMetre / MAX_LINE_DENSITY))
    const expected = droplets / weight
    let entries = Math.floor(expected)
    if (this.rng.next() < expected - entries) entries++
    const diffusivity = this.ionDiffusivity[i]
    const kind = KIND[seg.species]
    let added = 0
    for (let e = 0; e < entries; e++) {
      const along = this.rng.next() * seg.length
      // The ion wanders while it waits to nucleate. The wait is exponential, given that it
      // ended within the ion's lifetime.
      const delay = -Math.log(1 - this.rng.next() * chance) / rate
      const spread = Math.sqrt(2 * diffusivity * delay)
      const x = seg.start[0] + seg.dir[0] * along + spread * this.rng.normal()
      const y = seg.start[1] + seg.dir[1] * along + spread * this.rng.normal()
      const z = seg.start[2] + seg.dir[2] * along + spread * this.rng.normal()
      if (
        z <= 0 ||
        z >= this.box.height ||
        Math.abs(x) > this.box.halfX ||
        Math.abs(y) > this.box.halfY
      )
        continue
      if (this.add(x, y, z, weight, perMetre, this.time + delay, kind)) added++
    }
    return added
  }

  /** Send one particle through the chamber now. */
  fire(particle: Particle, source = 'manual'): EventRecord {
    let droplets = 0
    const first = this.count
    const result = transport(
      particle,
      this.profile,
      this.rng,
      (seg) => {
        droplets += this.deposit(seg)
      },
      this.transportOptions,
    )
    // The droplets of one track fall together (see advance); that needs the track's length.
    for (let k = first; k < this.count; k++) {
      const width = Math.max(TRACK_RADIUS, 1 / this.lineDensity[k])
      this.settle[k] = Math.max(
        0,
        Math.log(Math.max(result.pathLength, width) / width) + COLUMN_OFFSET,
      )
    }
    const record = {
      time: this.time,
      source,
      species: particle.species,
      kineticMeV: particle.kinetic,
      droplets,
    }
    this.events.push(record)
    if (this.events.length > 200) this.events.shift()
    return record
  }

  /** Advance the chamber by dt seconds. */
  step(dt: number): void {
    for (const source of this.sources) {
      const n = this.rng.poisson(source.rate * dt)
      for (let k = 0; k < n; k++) {
        const particle = enterBox(source.sample(this.rng), this.box)
        if (particle) this.fire(particle, source.id)
      }
    }
    if (this.options.mist && this.mistTotal > 0) this.spawnMist(dt)
    this.advance(dt)
    this.time += dt
  }

  private spawnMist(dt: number): void {
    const n = this.rng.poisson((this.mistTotal * dt) / this.mistWeight)
    const cdf = this.mistCdf
    const total = cdf[cdf.length - 1]
    for (let k = 0; k < n; k++) {
      const target = this.rng.next() * total
      let lo = 0
      let hi = cdf.length - 1
      while (hi - lo > 1) {
        const mid = (lo + hi) >> 1
        if (cdf[mid] < target) lo = mid
        else hi = mid
      }
      const z = (hi - this.rng.next()) * this.dz
      const x = (2 * this.rng.next() - 1) * this.box.halfX
      const y = (2 * this.rng.next() - 1) * this.box.halfY
      this.add(x, y, z, this.mistWeight, 0, this.time + this.rng.next() * dt, KIND.mist)
    }
  }

  /** Grow and drop every droplet; remove those that land or evaporate. */
  private advance(dt: number): void {
    const { pos, radius, weight, lineDensity, birth, settle, kind } = this
    const end = this.time + dt
    let write = 0
    for (let k = 0; k < this.count; k++) {
      let z = pos[3 * k + 2]
      let r = radius[k]
      const age = end - birth[k]
      let keep = true
      if (age > 0) {
        const h = Math.min(dt, age)
        const { i, f } = this.node(z)
        const growth = this.growth[i] + f * (this.growth[i + 1] - this.growth[i])
        const r2 = r * r + (2 * growth * h) / vaporShare(lineDensity[k], r, this.vaporDiff[i], age)
        if (r2 <= 0) keep = false
        else {
          r = Math.sqrt(r2)
          const kn =
            (this.meanFreePath[i] + f * (this.meanFreePath[i + 1] - this.meanFreePath[i])) / r
          const slip = 1 + kn * (1.257 + 0.4 * Math.exp(-1.1 / kn))
          const stokes = this.stokes[i] + f * (this.stokes[i + 1] - this.stokes[i])
          // A droplet's own fall, plus the downdraft of its whole track: the droplets weigh on
          // the same air and the column sinks together. Weight per length f = n m g gives
          // v = (f / 4 pi eta) K = 1.5 c n r^3 K, with c the Stokes coefficient.
          const unit = 1.5 * stokes * lineDensity[k] * r * r * r
          const draft = this.options.collectiveSettling
            ? unit * downdraftFactor(settle[k], unit, this.kinematic[i])
            : 0
          z -= (stokes * r * r * slip + draft) * h
          if (z <= 0) {
            keep = false
            this.landed += weight[k]
          }
        }
      }
      if (keep) {
        if (write !== k) {
          pos[3 * write] = pos[3 * k]
          pos[3 * write + 1] = pos[3 * k + 1]
          weight[write] = weight[k]
          lineDensity[write] = lineDensity[k]
          birth[write] = birth[k]
          settle[write] = settle[k]
          kind[write] = kind[k]
        }
        pos[3 * write + 2] = z
        radius[write] = r
        write++
      }
    }
    this.count = write
  }

  /** Real droplets currently represented, by kind. */
  census(): Record<keyof typeof KIND, number> {
    const out = { mist: 0, alpha: 0, electron: 0, muon: 0 }
    const names = ['mist', 'alpha', 'electron', 'muon'] as const
    for (let k = 0; k < this.count; k++)
      if (this.radius[k] > 0) out[names[this.kind[k]]] += this.weight[k]
    return out
  }
}
