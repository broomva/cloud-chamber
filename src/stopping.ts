/**
 * How fast charged particles lose energy to the electrons of the gas.
 *
 * Mass stopping powers are in MeV cm^2/g, energies in MeV. Formulas follow the PDG review
 * "Passage of particles through matter" (heavy particles) and ICRU Report 37 (electrons).
 */
import { AIR } from './air'
import { EV, K_BETHE, MASS_MEV, ME_C2 } from './constants'
import { ASTAR_ELECTRONIC, ASTAR_ENERGY } from './data/astar-air'

const I_MEV = AIR.meanExcitationEv * EV

export interface Kinematics {
  gamma: number
  beta2: number
  /** beta * gamma. */
  bg: number
  /** Momentum, MeV/c. */
  p: number
}

export function kinematics(kineticMeV: number, massMeV: number): Kinematics {
  const gamma = 1 + kineticMeV / massMeV
  const bg2 = gamma * gamma - 1
  return { gamma, beta2: bg2 / (gamma * gamma), bg: Math.sqrt(bg2), p: massMeV * Math.sqrt(bg2) }
}

/** Largest energy a particle of mass M can give a free electron in one collision, MeV. */
export function maxEnergyTransfer(kineticMeV: number, massMeV: number): number {
  const { gamma, bg } = kinematics(kineticMeV, massMeV)
  const r = ME_C2 / massMeV
  return (2 * ME_C2 * bg * bg) / (1 + 2 * gamma * r + r * r)
}

/** Sternheimer density-effect correction delta(beta gamma). Zero in air below beta gamma ~ 55. */
export function densityEffect(bg: number): number {
  const { a, k, x0, x1, cbar, delta0 } = AIR.sternheimer
  const x = Math.log10(bg)
  if (x >= x1) return 2 * Math.LN10 * x - cbar
  if (x >= x0) return 2 * Math.LN10 * x - cbar + a * (x1 - x) ** k
  return delta0 * 10 ** (2 * (x - x0))
}

/**
 * Bethe stopping number L for a heavy particle, with losses above `cutMeV` left out (restricted
 * loss). Pass cutMeV = Infinity for the full mean loss.
 */
function betheNumber(kineticMeV: number, massMeV: number, cutMeV: number): number {
  const { beta2, bg } = kinematics(kineticMeV, massMeV)
  const tMax = maxEnergyTransfer(kineticMeV, massMeV)
  const tUp = Math.min(cutMeV, tMax)
  return (
    0.5 * Math.log((2 * ME_C2 * bg * bg * tUp) / (I_MEV * I_MEV)) -
    (beta2 / 2) * (1 + tUp / tMax) -
    densityEffect(bg) / 2
  )
}

/**
 * Mean electronic stopping power of a heavy charged particle from the Bethe formula.
 * `cutMeV` restricts the loss to collisions that transfer less than the cut. Valid while the
 * particle is much faster than the atomic electrons (it reproduces the PDG muon table and NIST
 * PSTAR above a few MeV; see the tests). It has no shell or charge-exchange corrections, so it
 * is not used for slow alphas.
 */
export function betheStopping(
  z: number,
  massMeV: number,
  kineticMeV: number,
  cutMeV = Infinity,
): number {
  if (kineticMeV <= 0) return 0
  const { beta2 } = kinematics(kineticMeV, massMeV)
  const number = betheNumber(kineticMeV, massMeV, cutMeV)
  return Math.max((K_BETHE * z * z * AIR.zOverA * number) / beta2, 0)
}

const LOG_E = ASTAR_ENERGY.map(Math.log)
const LOG_S = ASTAR_ELECTRONIC.map(Math.log)
export const ALPHA_TABLE_MIN = ASTAR_ENERGY[0]
export const ALPHA_TABLE_MAX = ASTAR_ENERGY[ASTAR_ENERGY.length - 1]

/**
 * Electronic stopping power of an alpha particle in air: the measured ICRU 49 table (NIST
 * ASTAR), interpolated log-log. A slow alpha keeps capturing and losing electrons, and no short
 * formula gives its Bragg peak, so this one quantity is taken from measurement. With a finite
 * cut, the share of the loss carried by delta rays above the cut is removed using the ratio of
 * restricted to full Bethe stopping numbers.
 */
export function alphaStopping(kineticMeV: number, cutMeV = Infinity): number {
  if (kineticMeV <= 0) return 0
  const le = Math.log(Math.min(Math.max(kineticMeV, ALPHA_TABLE_MIN), ALPHA_TABLE_MAX))
  let lo = 0
  let hi = LOG_E.length - 1
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1
    if (LOG_E[mid] <= le) lo = mid
    else hi = mid
  }
  const f = (le - LOG_E[lo]) / (LOG_E[hi] - LOG_E[lo])
  let s = Math.exp(LOG_S[lo] + f * (LOG_S[hi] - LOG_S[lo]))
  // Below the table the loss is proportional to velocity.
  if (kineticMeV < ALPHA_TABLE_MIN) s *= Math.sqrt(kineticMeV / ALPHA_TABLE_MIN)
  const mass = MASS_MEV.alpha
  if (cutMeV < maxEnergyTransfer(kineticMeV, mass)) {
    const full = betheNumber(kineticMeV, mass, Infinity)
    const restricted = betheNumber(kineticMeV, mass, cutMeV)
    if (full > 0 && restricted > 0) s *= restricted / full
  }
  return s
}

/** Delta rays above `cutMeV` produced per g/cm^2 by a heavy particle (PDG eq. 34.8, spin 0). */
export function heavyDeltaRate(
  z: number,
  massMeV: number,
  kineticMeV: number,
  cutMeV: number,
): number {
  const tMax = maxEnergyTransfer(kineticMeV, massMeV)
  if (cutMeV >= tMax) return 0
  const { beta2 } = kinematics(kineticMeV, massMeV)
  return (
    ((K_BETHE / 2) *
      z *
      z *
      AIR.zOverA *
      (1 / cutMeV - 1 / tMax - (beta2 / tMax) * Math.log(tMax / cutMeV))) /
    beta2
  )
}

/** Sample a delta-ray energy above the cut from d sigma/dT ~ (1 - beta^2 T/Tmax) / T^2. */
export function sampleHeavyDelta(
  kineticMeV: number,
  massMeV: number,
  cutMeV: number,
  u1: number,
  u2: () => number,
): number {
  const tMax = maxEnergyTransfer(kineticMeV, massMeV)
  const { beta2 } = kinematics(kineticMeV, massMeV)
  let u = u1
  for (let i = 0; i < 100; i++) {
    const t = 1 / (1 / cutMeV - u * (1 / cutMeV - 1 / tMax)) // 1/T^2 by inversion
    if (u2() <= 1 - (beta2 * t) / tMax) return t
    u = u2()
  }
  return cutMeV
}

/** Moller cross-section shape: d sigma/d eps for energy fraction eps given to the slower electron. */
function mollerShape(eps: number, tau: number): number {
  const g = tau / (tau + 1)
  return (
    1 / (eps * eps) +
    1 / ((1 - eps) * (1 - eps)) +
    g * g -
    (2 * tau + 1) / ((tau + 1) * (tau + 1)) / (eps * (1 - eps))
  )
}

/** Prefactor of the Moller cross-section per g/cm^2 per unit eps: (K/2)(Z/A) / (beta^2 tau). */
function mollerPrefactor(kineticMeV: number): number {
  const { beta2 } = kinematics(kineticMeV, ME_C2)
  return ((K_BETHE / 2) * AIR.zOverA) / (beta2 * (kineticMeV / ME_C2))
}

/** Composite Simpson integral of f over [a, b] in log-spaced panels. */
function integrateLog(f: (x: number) => number, a: number, b: number, panels = 200): number {
  if (b <= a) return 0
  const la = Math.log(a)
  const h = (Math.log(b) - la) / panels
  let sum = 0
  for (let i = 0; i <= panels; i++) {
    const x = Math.exp(la + i * h)
    const w = i === 0 || i === panels ? 1 : i % 2 === 1 ? 4 : 2
    sum += w * f(x) * x
  }
  return (sum * h) / 3
}

/** Delta rays above `cutMeV` produced per g/cm^2 by an electron (Moller). */
export function electronDeltaRate(kineticMeV: number, cutMeV: number): number {
  const epsCut = cutMeV / kineticMeV
  if (epsCut >= 0.5) return 0
  const tau = kineticMeV / ME_C2
  return mollerPrefactor(kineticMeV) * integrateLog((e) => mollerShape(e, tau), epsCut, 0.5)
}

/** Sample the energy of a Moller delta ray above the cut, by rejection under 1/eps^2. */
export function sampleElectronDelta(
  kineticMeV: number,
  cutMeV: number,
  rand: () => number,
): number {
  const epsCut = cutMeV / kineticMeV
  const tau = kineticMeV / ME_C2
  // The shape times eps^2 is largest at one of the two ends of [epsCut, 1/2].
  const bound = Math.max(mollerShape(epsCut, tau) * epsCut * epsCut, mollerShape(0.5, tau) * 0.25)
  for (let i = 0; i < 1000; i++) {
    const eps = 1 / (1 / epsCut - rand() * (1 / epsCut - 2))
    if (rand() * bound <= mollerShape(eps, tau) * eps * eps) return eps * kineticMeV
  }
  return cutMeV
}

/**
 * Collision stopping power of an electron in air (ICRU 37, Bethe with Moller exchange). With a
 * finite `cutMeV`, the energy carried off by delta rays above the cut is subtracted.
 */
export function electronStopping(kineticMeV: number, cutMeV = Infinity): number {
  if (kineticMeV <= 0) return 0
  const tau = kineticMeV / ME_C2
  const { beta2, bg } = kinematics(kineticMeV, ME_C2)
  const f = (1 - beta2) * (1 + (tau * tau) / 8 - (2 * tau + 1) * Math.LN2)
  const number =
    Math.log((kineticMeV * kineticMeV) / (I_MEV * I_MEV)) +
    Math.log(1 + tau / 2) +
    f -
    densityEffect(bg)
  const total = ((K_BETHE / 2) * AIR.zOverA * number) / beta2
  const epsCut = cutMeV / kineticMeV
  if (epsCut >= 0.5) return Math.max(total, 0)
  const carried =
    mollerPrefactor(kineticMeV) *
    kineticMeV *
    integrateLog((e) => e * mollerShape(e, tau), epsCut, 0.5)
  return Math.max(total - carried, 0)
}
