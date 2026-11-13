/**
 * Elastic (Coulomb) scattering off the nuclei of the gas, one collision at a time.
 *
 * Each collision follows the screened Rutherford cross-section
 *   d sigma / d Omega = A^2 / (1 - cos theta + 2 eta)^2
 * with Moliere's screening angle: eta = chi_a^2 / 4,
 *   chi_a^2 = chi_0^2 (1.13 + 3.76 (alpha z Z / beta)^2),  chi_0 = hbar / (p a_TF).
 * In mu = (1 - cos theta) / 2 the total is sigma = pi A^2 / (eta (1 + eta)) and the angle can
 * be sampled exactly: mu = eta xi / (1 + eta - xi).
 *
 * In air a fast particle makes tens to hundreds of these collisions per centimetre, so they are
 * all simulated; nothing is replaced by a multiple-scattering formula. The curl of a slow
 * electron's track and the rare kink in an alpha's come from this.
 */
import { AIR } from './air'
import { A_BOHR_CM, ALPHA_FS, AMU_MEV, HBAR_C_MEV_CM, ME_C2, R_E_CM } from './constants'
import { type Kinematics, kinematics } from './stopping'

export interface Projectile {
  z: number
  massMeV: number
  /** Electrons also scatter off atomic electrons, which makes Z^2 into Z (Z + 1). */
  isElectron: boolean
}

export interface ElementScattering {
  /** Cross-section per atom, cm^2. */
  sigma: number
  /** Screening parameter. */
  eta: number
  /** Projectile mass over target mass (0 for light projectiles: the nucleus does not recoil). */
  massRatio: number
}

/** Cross-section and screening for one element. */
export function elementScattering(
  proj: Projectile,
  kin: Kinematics,
  zTarget: number,
  aTarget: number,
): ElementScattering {
  const beta = Math.sqrt(kin.beta2)
  // Heavy projectiles: work in the centre of mass (non-relativistic, which alphas here are).
  const heavy = proj.massMeV > 1000
  const massRatio = heavy ? proj.massMeV / (aTarget * AMU_MEV) : 0
  const reduce = 1 / (1 + massRatio) // reduced mass / projectile mass
  const p = kin.p * reduce
  const zz = proj.isElectron ? Math.sqrt(zTarget * (zTarget + 1)) : proj.z * zTarget
  const amp = (zz * R_E_CM * ME_C2) / (p * beta)
  const screenLength =
    (0.8853 * A_BOHR_CM) /
    (proj.z > 1 ? Math.sqrt(proj.z ** (2 / 3) + zTarget ** (2 / 3)) : Math.cbrt(zTarget))
  const chi0 = HBAR_C_MEV_CM / (p * screenLength)
  const born = (ALPHA_FS * proj.z * zTarget) / beta
  const eta = (chi0 * chi0 * (1.13 + 3.76 * born * born)) / 4
  return { sigma: (Math.PI * amp * amp) / (eta * (1 + eta)), eta, massRatio }
}

export interface GasScattering {
  /** Collisions per cm at unit mass density (multiply by g/cm^3), i.e. cm^2/g. */
  perGram: number
  /** Per-element share of the collisions and its screening. */
  elements: { share: number; eta: number; massRatio: number }[]
}

/** Elastic scattering of a projectile in air at the given kinetic energy. */
export function gasScattering(proj: Projectile, kineticMeV: number): GasScattering {
  const kin = kinematics(kineticMeV, proj.massMeV)
  let perGram = 0
  const parts = AIR.elements.map((el) => {
    const s = elementScattering(proj, kin, el.z, el.a)
    const rate = (s.sigma * el.massFraction * 6.02214076e23) / el.a
    perGram += rate
    return { rate, eta: s.eta, massRatio: s.massRatio }
  })
  return {
    perGram,
    elements: parts.map((x) => ({ share: x.rate / perGram, eta: x.eta, massRatio: x.massRatio })),
  }
}

/** Sample mu = (1 - cos theta) / 2 for one collision. `xi` is uniform in (0, 1). */
export function sampleMu(eta: number, xi: number): number {
  return (eta * xi) / (1 + eta - xi)
}

/**
 * Turn a centre-of-mass deflection into what the projectile does in the lab: its lab cosine and
 * the fraction of its kinetic energy handed to the recoiling nucleus.
 */
export function labDeflection(
  mu: number,
  massRatio: number,
): { cos: number; energyFraction: number } {
  const cosCm = 1 - 2 * mu
  if (massRatio === 0) return { cos: cosCm, energyFraction: 0 }
  const sinCm = Math.sqrt(Math.max(0, 1 - cosCm * cosCm))
  const thetaLab = Math.atan2(sinCm, cosCm + massRatio)
  return {
    cos: Math.cos(thetaLab),
    energyFraction: (4 * massRatio * mu) / ((1 + massRatio) * (1 + massRatio)),
  }
}

/** Rotate unit vector d by polar angle (given as its cosine) and azimuth phi. */
export function deflect(
  d: [number, number, number],
  cosTheta: number,
  phi: number,
): [number, number, number] {
  const sinTheta = Math.sqrt(Math.max(0, 1 - cosTheta * cosTheta))
  const [x, y, z] = d
  const cp = Math.cos(phi)
  const sp = Math.sin(phi)
  const t = Math.sqrt(x * x + y * y)
  if (t < 1e-10) {
    // Exactly along z: any perpendicular frame will do.
    return [sinTheta * cp, sinTheta * sp, z > 0 ? cosTheta : -cosTheta]
  }
  const nx = x * cosTheta + (sinTheta * (x * z * cp - y * sp)) / t
  const ny = y * cosTheta + (sinTheta * (y * z * cp + x * sp)) / t
  const nz = z * cosTheta - t * sinTheta * cp
  const norm = Math.hypot(nx, ny, nz)
  return [nx / norm, ny / norm, nz / norm]
}
