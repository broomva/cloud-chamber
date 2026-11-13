/**
 * One droplet: how fast it grows from supersaturated vapor and how fast it falls.
 */
import { airConductivity, airMeanFreePath, airViscosity } from './air'
import { G_ACCEL, R_GAS } from './constants'
import {
  IPA,
  liquidDensity,
  saturationPressure,
  vaporDiffusivity,
  vaporizationEnthalpy,
} from './ipa'

/**
 * Growth coefficient g in r dr/dt = g (S - 1), m^2/s (Mason's equation). Two resistances in
 * series: vapor must diffuse to the droplet, and the latent heat it releases must conduct away.
 */
export function growthCoefficient(tempK: number, pressurePa: number): number {
  const rhoL = liquidDensity(tempK)
  const rv = R_GAS / IPA.molarMass // specific gas constant of the vapor, J/(kg K)
  const latent = vaporizationEnthalpy(tempK) / IPA.molarMass // J/kg
  const heat = ((latent / (rv * tempK) - 1) * latent * rhoL) / (airConductivity(tempK) * tempK)
  const vapor =
    (rhoL * rv * tempK) / (vaporDiffusivity(tempK, pressurePa) * saturationPressure(tempK))
  return 1 / (heat + vapor)
}

/** Terminal fall speed, m/s: Stokes drag with the Cunningham slip correction. */
export function fallSpeed(radius: number, tempK: number, pressurePa: number): number {
  const kn = airMeanFreePath(tempK, pressurePa) / radius
  const slip = 1 + kn * (1.257 + 0.4 * Math.exp(-1.1 / kn))
  return ((2 * liquidDensity(tempK) * G_ACCEL * radius * radius) / (9 * airViscosity(tempK))) * slip
}

/** The Stokes coefficient c in v = c r^2 (no slip), 1/(m s). */
export function stokesCoefficient(tempK: number): number {
  return (2 * liquidDensity(tempK) * G_ACCEL) / (9 * airViscosity(tempK))
}
