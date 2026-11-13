/**
 * Where the radiation comes from. Each source says how often it fires and samples a particle.
 *
 * Rates are set from stated physical quantities (a muon intensity, a radon concentration, a
 * gamma dose rate), not chosen for looks. Nuclear data: alpha energies from the evaluated decay
 * data for Am-241, Rn-222, Po-218 and Po-214; beta end points for Sr-90 and Y-90.
 */
import { airDensity } from './air'
import { ALPHA_FS, ME_C2 } from './constants'
import type { Rng } from './rng'
import type { Particle } from './transport'

export interface Box {
  halfX: number
  halfY: number
  height: number
}

export interface Source {
  id: string
  label: string
  /** Events per second for this chamber. */
  rate: number
  sample(rng: Rng): Particle
}

function isotropic(rng: Rng): [number, number, number] {
  const cos = 2 * rng.next() - 1
  const sin = Math.sqrt(1 - cos * cos)
  const phi = 2 * Math.PI * rng.next()
  return [sin * Math.cos(phi), sin * Math.sin(phi), cos]
}

/** Vertical intensity of muons at sea level, per m^2 per s per sr (PDG, cosmic rays review). */
export const MUON_VERTICAL_INTENSITY = 70

/**
 * Cosmic-ray muons. Intensity falls as cos^2 of the zenith angle, so the flux through a
 * horizontal surface is (pi/2) I_v. Muons are started on the plane of the chamber top, over an
 * area widened by `margin` so that slanted tracks entering through the walls are included; the
 * few beyond that margin (about 1% of the flux at the default) are left out.
 *
 * Energies come from a rough log-normal around a few GeV. Every value is near minimum
 * ionization, so the tracks do not depend on the shape of this spectrum.
 */
export function cosmicMuons(box: Box, margin = 3 * box.height): Source {
  const lx = 2 * (box.halfX + margin)
  const ly = 2 * (box.halfY + margin)
  return {
    id: 'muon',
    label: 'Cosmic muons',
    rate: (Math.PI / 2) * MUON_VERTICAL_INTENSITY * lx * ly,
    sample(rng) {
      // Flux through a horizontal plane ~ cos^3 sin: cos = u^(1/4).
      const cos = rng.next() ** 0.25
      const sin = Math.sqrt(1 - cos * cos)
      const phi = 2 * Math.PI * rng.next()
      return {
        species: 'muon',
        kinetic: 3000 * Math.exp(0.9 * rng.normal()),
        pos: [(rng.next() - 0.5) * lx, (rng.next() - 0.5) * ly, box.height],
        dir: [sin * Math.cos(phi), sin * Math.sin(phi), -cos],
      }
    },
  }
}

/** Advance a particle that starts outside the box to where it enters, or return null on a miss. */
export function enterBox(p: Particle, box: Box): Particle | null {
  const lo = [-box.halfX, -box.halfY, 0]
  const hi = [box.halfX, box.halfY, box.height]
  let tNear = 0
  let tFar = Infinity
  for (let i = 0; i < 3; i++) {
    if (Math.abs(p.dir[i]) < 1e-12) {
      if (p.pos[i] < lo[i] || p.pos[i] > hi[i]) return null
      continue
    }
    const a = (lo[i] - p.pos[i]) / p.dir[i]
    const b = (hi[i] - p.pos[i]) / p.dir[i]
    tNear = Math.max(tNear, Math.min(a, b))
    tFar = Math.min(tFar, Math.max(a, b))
  }
  if (tNear >= tFar) return null
  const t = tNear + 1e-9
  return { ...p, pos: [p.pos[0] + p.dir[0] * t, p.pos[1] + p.dir[1] * t, p.pos[2] + p.dir[2] * t] }
}

/**
 * Radon in the chamber air. Rn-222 and its short-lived daughters Po-218 and Po-214 each emit an
 * alpha. `concentration` is Bq/m^3 of radon; `equilibrium` is the fraction of the daughters that
 * stay airborne instead of plating out (about 0.4 indoors).
 */
export function radonAlphas(box: Box, concentration = 40, equilibrium = 0.4): Source {
  const volume = 4 * box.halfX * box.halfY * box.height
  const lines = [
    { energy: 5.4895, weight: 1 }, // Rn-222
    { energy: 6.0024, weight: equilibrium }, // Po-218
    { energy: 7.6868, weight: equilibrium }, // Po-214
  ]
  const total = lines.reduce((s, l) => s + l.weight, 0)
  return {
    id: 'radon',
    label: 'Radon in the air',
    rate: concentration * volume * total,
    sample(rng) {
      let pick = rng.next() * total
      let energy = lines[0].energy
      for (const line of lines) {
        if (pick < line.weight) {
          energy = line.energy
          break
        }
        pick -= line.weight
      }
      return {
        species: 'alpha',
        kinetic: energy,
        pos: [
          (2 * rng.next() - 1) * box.halfX,
          (2 * rng.next() - 1) * box.halfY,
          rng.next() * box.height,
        ],
        dir: isotropic(rng),
      }
    },
  }
}

/** A point alpha source (Am-241) at `pos`, emitting into a cone around `axis`. */
export function americiumSource(
  pos: [number, number, number],
  axis: [number, number, number],
  rate = 4,
  halfAngle = 0.5,
): Source {
  const norm = Math.hypot(...axis)
  const a: [number, number, number] = [axis[0] / norm, axis[1] / norm, axis[2] / norm]
  // Any unit vector perpendicular to the axis, and the third to complete the frame.
  const ref: [number, number, number] = Math.abs(a[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0]
  const u: [number, number, number] = [
    a[1] * ref[2] - a[2] * ref[1],
    a[2] * ref[0] - a[0] * ref[2],
    a[0] * ref[1] - a[1] * ref[0],
  ]
  const un = Math.hypot(...u)
  const e1: [number, number, number] = [u[0] / un, u[1] / un, u[2] / un]
  const e2: [number, number, number] = [
    a[1] * e1[2] - a[2] * e1[1],
    a[2] * e1[0] - a[0] * e1[2],
    a[0] * e1[1] - a[1] * e1[0],
  ]
  return {
    id: 'americium',
    label: 'Am-241 alpha source',
    rate,
    sample(rng) {
      const cos = 1 - rng.next() * (1 - Math.cos(halfAngle))
      const sin = Math.sqrt(1 - cos * cos)
      const phi = 2 * Math.PI * rng.next()
      const c = sin * Math.cos(phi)
      const s = sin * Math.sin(phi)
      return {
        species: 'alpha',
        // 5.4856 MeV (85%) and 5.4428 MeV (13%); the remaining 2% of weaker lines go to the second.
        kinetic: rng.next() < 0.848 ? 5.4856 : 5.4428,
        pos: [...pos],
        dir: [
          a[0] * cos + e1[0] * c + e2[0] * s,
          a[1] * cos + e1[1] * c + e2[1] * s,
          a[2] * cos + e1[2] * c + e2[2] * s,
        ],
      }
    },
  }
}

/**
 * Beta-minus spectrum of a first-forbidden unique decay:
 *   N(T) ~ F(Z, E) p E (Q - T)^2 (p^2 + q^2)
 * with the non-relativistic Fermi function F = x / (1 - e^-x), x = 2 pi alpha Z E / p, for the
 * pull of the daughter nucleus (charge Z) on the electron. Sampled by rejection.
 */
export function sampleBeta(endpointMeV: number, zDaughter: number, rng: Rng): number {
  const shape = (t: number) => {
    const e = t + ME_C2
    const p2 = t * (t + 2 * ME_C2)
    const p = Math.sqrt(p2)
    const q = endpointMeV - t
    const x = (2 * Math.PI * ALPHA_FS * zDaughter * e) / p
    return (x / (1 - Math.exp(-x))) * p * e * q * q * (p2 + q * q)
  }
  let peak = 0
  for (let i = 1; i < 400; i++) peak = Math.max(peak, shape((i / 400) * endpointMeV))
  for (let i = 0; i < 10000; i++) {
    const t = rng.next() * endpointMeV
    if (rng.next() * peak * 1.02 <= shape(t)) return t
  }
  return endpointMeV / 3
}

/** A point Sr-90 source in equilibrium with Y-90: two betas per decay, end points 0.5459 and 2.2785 MeV (ENSDF). */
export function strontiumSource(pos: [number, number, number], rate = 6): Source {
  return {
    id: 'strontium',
    label: 'Sr-90 beta source',
    rate,
    sample(rng) {
      return {
        species: 'electron',
        kinetic: rng.next() < 0.5 ? sampleBeta(0.5459, 39, rng) : sampleBeta(2.2785, 40, rng),
        pos: [...pos],
        dir: isotropic(rng),
      }
    },
  }
}

/** Klein-Nishina: sample the kinetic energy and direction cosine of a Compton electron. */
export function sampleCompton(photonMeV: number, rng: Rng): { kinetic: number; cos: number } {
  const k = photonMeV / ME_C2
  const epsMin = 1 / (1 + 2 * k)
  for (let i = 0; i < 10000; i++) {
    // eps = scattered / incident photon energy; d sigma/d eps ~ (1/eps + eps)(1 - eps sin^2 / (1 + eps^2)).
    const eps = epsMin ** rng.next() // 1/eps part by inversion; eps part via rejection weight
    const t = (1 - eps) / (k * eps) // 1 - cos(theta_photon)
    const sin2 = t * (2 - t)
    const weight = ((1 + eps * eps) / 2) * (1 - (eps * sin2) / (1 + eps * eps))
    if (rng.next() <= weight) {
      const kinetic = photonMeV * (1 - eps)
      // Electron angle from momentum conservation along the photon direction.
      const pe = Math.sqrt(kinetic * (kinetic + 2 * ME_C2))
      const cosPhoton = 1 - t
      return { kinetic, cos: Math.min(1, (photonMeV - photonMeV * eps * cosPhoton) / pe) }
    }
  }
  return { kinetic: photonMeV / 2, cos: 1 }
}

/** Energy of the K-40 gamma line, the strongest single line in the natural background, MeV. */
export const K40_GAMMA = 1.4608

/**
 * Electrons knocked loose in the gas by background gamma rays. The rate follows from an ambient
 * absorbed dose rate in air (nGy/h): energy absorbed in the chamber gas per second, divided by
 * the mean Compton-electron energy. All photons are given the K-40 line energy. Electrons that
 * enter from the walls are not included.
 */
export function ambientElectrons(
  box: Box,
  gasTempK: number,
  pressurePa: number,
  doseNGyPerHour = 60,
): Source {
  const volume = 4 * box.halfX * box.halfY * box.height
  const massKg = airDensity(gasTempK, pressurePa) * 1e3 * volume
  const mevPerSecond = (doseNGyPerHour * 1e-9 * massKg) / 3600 / 1.602176634e-13
  // Mean Compton-electron energy for the line, from the sampler itself (fixed seed).
  let mean = 0
  const probe = { s: 12345 }
  const lcg = () => {
    probe.s = (Math.imul(probe.s, 1664525) + 1013904223) >>> 0
    return (probe.s + 0.5) / 4294967296
  }
  const fake = { next: lcg } as unknown as Rng
  for (let i = 0; i < 4000; i++) mean += sampleCompton(K40_GAMMA, fake).kinetic / 4000
  return {
    id: 'ambient',
    label: 'Background gamma electrons',
    rate: mevPerSecond / mean,
    sample(rng) {
      const { kinetic, cos } = sampleCompton(K40_GAMMA, rng)
      const photon = isotropic(rng)
      const sin = Math.sqrt(1 - cos * cos)
      const phi = 2 * Math.PI * rng.next()
      // Build the electron direction around the photon direction.
      const ref: [number, number, number] = Math.abs(photon[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0]
      const u = [
        photon[1] * ref[2] - photon[2] * ref[1],
        photon[2] * ref[0] - photon[0] * ref[2],
        photon[0] * ref[1] - photon[1] * ref[0],
      ]
      const un = Math.hypot(u[0], u[1], u[2])
      const e1 = [u[0] / un, u[1] / un, u[2] / un]
      const e2 = [
        photon[1] * e1[2] - photon[2] * e1[1],
        photon[2] * e1[0] - photon[0] * e1[2],
        photon[0] * e1[1] - photon[1] * e1[0],
      ]
      const c = sin * Math.cos(phi)
      const s = sin * Math.sin(phi)
      return {
        species: 'electron',
        kinetic,
        pos: [
          (2 * rng.next() - 1) * box.halfX,
          (2 * rng.next() - 1) * box.halfY,
          rng.next() * box.height,
        ],
        dir: [
          photon[0] * cos + e1[0] * c + e2[0] * s,
          photon[1] * cos + e1[1] * c + e2[1] * s,
          photon[2] * cos + e1[2] * c + e2[2] * s,
        ],
      }
    },
  }
}
