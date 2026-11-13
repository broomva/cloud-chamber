/**
 * Follow one charged particle through the chamber gas and report where it leaves ions.
 *
 * Between collisions the particle loses energy continuously to soft collisions (restricted
 * stopping power) and that energy becomes ion pairs. Hard collisions are simulated one by one:
 * elastic scattering off nuclei bends the path, and knock-on electrons (delta rays) above the
 * cut are launched and followed as particles themselves.
 */
import { AIR, airDensity } from './air'
import { type ChamberProfile, sampleProfile } from './chamber'
import { EV, MASS_MEV, ME_C2 } from './constants'
import type { Rng } from './rng'
import { deflect, gasScattering, labDeflection, type Projectile, sampleMu } from './scattering'
import {
  alphaStopping,
  betheStopping,
  electronDeltaRate,
  electronStopping,
  heavyDeltaRate,
  kinematics,
  maxEnergyTransfer,
  sampleElectronDelta,
  sampleHeavyDelta,
} from './stopping'

export type Species = 'alpha' | 'electron' | 'muon'

export const PROJECTILES: Record<Species, Projectile> = {
  alpha: { z: 2, massMeV: MASS_MEV.alpha, isElectron: false },
  electron: { z: 1, massMeV: MASS_MEV.electron, isElectron: true },
  muon: { z: 1, massMeV: MASS_MEV.muon, isElectron: false },
}

export interface Particle {
  species: Species
  /** Kinetic energy, MeV. */
  kinetic: number
  /** Position, m. The plate is z = 0. */
  pos: [number, number, number]
  /** Unit direction. */
  dir: [number, number, number]
}

/** A straight piece of track and the ion pairs made along it. */
export interface Segment {
  species: Species
  start: [number, number, number]
  dir: [number, number, number]
  /** Length, m (0 for an energy dump at a point). */
  length: number
  ionPairs: number
}

export interface TransportOptions {
  /** Knock-on electrons above this energy are followed as particles, MeV. */
  cutMeV: number
  /** Below this kinetic energy a particle stops and leaves its energy where it is, MeV. */
  absorbMeV: Record<Species, number>
  /** Longest step, m. */
  maxStep: number
  /** Largest fraction of its energy a particle may lose in one step. */
  maxLossFraction: number
  /** Half-widths of the box in x and y, m. */
  halfX: number
  halfY: number
  /** Follow knock-on electrons (off only for tests that isolate the primary). */
  followDeltas: boolean
}

export const DEFAULT_TRANSPORT: TransportOptions = {
  cutMeV: 1e-3,
  absorbMeV: { alpha: 5e-3, electron: 1e-3, muon: 0.1 },
  maxStep: 5e-4,
  maxLossFraction: 0.05,
  halfX: 0.1,
  halfY: 0.1,
  followDeltas: true,
}

export interface TransportResult {
  /** Path length travelled by the primary, m. */
  pathLength: number
  /** Where the primary stopped or left the box. */
  end: [number, number, number]
  /** True if the primary stopped inside the chamber. */
  stopped: boolean
  /** Energy turned into ion pairs along all tracks (primary and knock-ons), MeV. */
  depositedMeV: number
  /** Energy handed to recoiling nuclei, MeV. */
  nuclearMeV: number
  elasticCollisions: number
  deltaRays: number
  steps: number
}

function restrictedStopping(species: Species, kinetic: number, cut: number): number {
  if (species === 'alpha') return alphaStopping(kinetic, cut)
  if (species === 'electron') return electronStopping(kinetic, cut)
  return betheStopping(1, MASS_MEV.muon, kinetic, cut)
}

function deltaRate(species: Species, kinetic: number, cut: number): number {
  if (species === 'electron') return electronDeltaRate(kinetic, cut)
  const p = PROJECTILES[species]
  return heavyDeltaRate(p.z, p.massMeV, kinetic, cut)
}

/** Mean energy per ion pair for the species, MeV. */
export function wValue(species: Species): number {
  return (species === 'alpha' ? AIR.wAlphaEv : AIR.wElectronEv) * EV
}

function randomAzimuth(rng: Rng): number {
  return 2 * Math.PI * rng.next()
}

/**
 * Transport a particle and every knock-on electron it makes. `emit` receives each track segment
 * as it is produced.
 */
export function transport(
  primary: Particle,
  profile: ChamberProfile,
  rng: Rng,
  emit: (segment: Segment) => void,
  options: TransportOptions = DEFAULT_TRANSPORT,
): TransportResult {
  const height = profile.config.height
  const pressure = profile.config.pressurePa
  const result: TransportResult = {
    pathLength: 0,
    end: [...primary.pos],
    stopped: false,
    depositedMeV: 0,
    nuclearMeV: 0,
    elasticCollisions: 0,
    deltaRays: 0,
    steps: 0,
  }
  const stack: { particle: Particle; primary: boolean }[] = [
    { particle: { ...primary, pos: [...primary.pos], dir: [...primary.dir] }, primary: true },
  ]

  while (stack.length > 0) {
    const item = stack.pop()
    if (!item) break
    const { particle, primary: isPrimary } = item
    const { species } = particle
    const proj = PROJECTILES[species]
    const w = wValue(species)
    const absorb = options.absorbMeV[species]
    let { kinetic } = particle
    let [x, y, z] = particle.pos
    let dir = particle.dir
    let inside = true
    // Distances to the next hard collisions, in g/cm^2 so they survive density changes.
    let toElastic = rng.exp(1)
    let toDelta = rng.exp(1)

    while (inside && kinetic > absorb) {
      result.steps++
      const rho = airDensity(sampleProfile(profile, profile.tempK, z), pressure) // g/cm^3
      const stop = restrictedStopping(species, kinetic, options.cutMeV) // MeV cm^2/g
      const scatter = gasScattering(proj, kinetic)
      const dRate = deltaRate(species, kinetic, options.cutMeV) // per g/cm^2

      // Step length in g/cm^2: the nearest of the next collision, the loss cap and the step cap.
      const elasticLen = toElastic / scatter.perGram
      const deltaLen = dRate > 0 ? toDelta / dRate : Infinity
      const lossLen = (options.maxLossFraction * kinetic) / stop
      const capLen = options.maxStep * 100 * rho
      let len = Math.min(elasticLen, deltaLen, lossLen, capLen)
      const event = len === elasticLen ? 'elastic' : len === deltaLen ? 'delta' : 'none'

      // Clip to the chamber walls.
      let stepM = len / rho / 100
      let exit = false
      const clip = (pos: number, d: number, lo: number, hi: number) => {
        if (d > 0 && pos + d * stepM > hi) {
          stepM = (hi - pos) / d
          exit = true
        } else if (d < 0 && pos + d * stepM < lo) {
          stepM = (lo - pos) / d
          exit = true
        }
      }
      clip(x, dir[0], -options.halfX, options.halfX)
      clip(y, dir[1], -options.halfY, options.halfY)
      clip(z, dir[2], 0, height)
      if (exit) len = stepM * 100 * rho

      const loss = Math.min(stop * len, kinetic)
      emit({ species, start: [x, y, z], dir, length: stepM, ionPairs: loss / w })
      result.depositedMeV += loss
      kinetic -= loss
      x += dir[0] * stepM
      y += dir[1] * stepM
      z += dir[2] * stepM
      if (isPrimary) result.pathLength += stepM
      toElastic -= len * scatter.perGram
      toDelta -= len * dRate
      if (exit) {
        inside = false
        break
      }
      if (kinetic <= absorb) break

      if (event === 'elastic') {
        result.elasticCollisions++
        toElastic = rng.exp(1)
        let pick = rng.next()
        let el = scatter.elements[scatter.elements.length - 1]
        for (const candidate of scatter.elements) {
          if (pick < candidate.share) {
            el = candidate
            break
          }
          pick -= candidate.share
        }
        const lab = labDeflection(sampleMu(el.eta, rng.next()), el.massRatio)
        dir = deflect(dir, lab.cos, randomAzimuth(rng))
        const recoil = lab.energyFraction * kinetic
        if (recoil > 0) {
          // The struck nucleus stops within micrometres; its energy ionizes on the spot.
          kinetic -= recoil
          result.nuclearMeV += recoil
          result.depositedMeV += recoil
          emit({ species, start: [x, y, z], dir, length: 0, ionPairs: recoil / w })
        }
      } else if (event === 'delta') {
        result.deltaRays++
        toDelta = rng.exp(1)
        const phi = randomAzimuth(rng)
        let tDelta: number
        let cosDelta: number
        if (species === 'electron') {
          tDelta = sampleElectronDelta(kinetic, options.cutMeV, () => rng.next())
          // Two-body kinematics for both electrons.
          cosDelta = Math.sqrt((tDelta * (kinetic + 2 * ME_C2)) / (kinetic * (tDelta + 2 * ME_C2)))
          const rest = kinetic - tDelta
          const cosPrimary = Math.sqrt(
            (rest * (kinetic + 2 * ME_C2)) / (kinetic * (rest + 2 * ME_C2)),
          )
          const deltaDir = deflect(dir, cosDelta, phi)
          dir = deflect(dir, cosPrimary, phi + Math.PI)
          kinetic = rest
          if (options.followDeltas) {
            stack.push({
              particle: { species: 'electron', kinetic: tDelta, pos: [x, y, z], dir: deltaDir },
              primary: false,
            })
          }
        } else {
          tDelta = sampleHeavyDelta(kinetic, proj.massMeV, options.cutMeV, rng.next(), () =>
            rng.next(),
          )
          const tMax = maxEnergyTransfer(kinetic, proj.massMeV)
          const pe = Math.sqrt(tDelta * (tDelta + 2 * ME_C2))
          const pMax = Math.sqrt(tMax * (tMax + 2 * ME_C2))
          cosDelta = Math.min(1, ((tDelta / pe) * pMax) / tMax)
          kinetic -= tDelta
          if (options.followDeltas) {
            stack.push({
              particle: {
                species: 'electron',
                kinetic: tDelta,
                pos: [x, y, z],
                dir: deflect(dir, cosDelta, phi),
              },
              primary: false,
            })
          }
        }
        if (!options.followDeltas) {
          result.depositedMeV += tDelta
          emit({
            species: 'electron',
            start: [x, y, z],
            dir,
            length: 0,
            ionPairs: tDelta / wValue('electron'),
          })
        }
      }
    }

    if (inside && kinetic > 0) {
      // Too slow to follow: what is left ionizes where the particle stops.
      emit({ species, start: [x, y, z], dir, length: 0, ionPairs: kinetic / w })
      result.depositedMeV += kinetic
    }
    if (isPrimary) {
      result.end = [x, y, z]
      result.stopped = inside
    }
  }
  return result
}

/** Speed as a fraction of c, for reporting. */
export function betaOf(species: Species, kineticMeV: number): number {
  return Math.sqrt(kinematics(kineticMeV, PROJECTILES[species].massMeV).beta2)
}
