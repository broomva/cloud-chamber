/**
 * Dry air as a stopping medium and as a carrier gas.
 *
 * Stopping constants are from the PDG "Atomic and nuclear properties" page for air (dry, 1 atm),
 * 2024 edition: https://pdg.lbl.gov/2024/AtomicNuclearProperties/HTML/air_dry_1_atm.html
 * and its muon table (Sternheimer coefficients):
 * https://pdg.lbl.gov/2024/AtomicNuclearProperties/MUE/muE_air_dry_1_atm.txt
 */
import { N_A, R_GAS } from './constants'

export const AIR = {
  zOverA: 0.49919, // mol/g
  meanExcitationEv: 85.7,
  radiationLength: 36.62, // g/cm^2
  molarMass: 0.0289647, // kg/mol
  /** Sternheimer density-effect coefficients. */
  sternheimer: { a: 0.1091, k: 3.3994, x0: 1.7418, x1: 4.2759, cbar: 10.5961, delta0: 0 },
  /** Composition by mass fraction (PDG). */
  elements: [
    { symbol: 'N', z: 7, a: 14.007, massFraction: 0.755267 },
    { symbol: 'O', z: 8, a: 15.999, massFraction: 0.231781 },
    { symbol: 'Ar', z: 18, a: 39.948, massFraction: 0.012827 },
    { symbol: 'C', z: 6, a: 12.0107, massFraction: 0.000124 },
  ],
  /** Mean energy per ion pair, eV (ICRU Report 31): electrons and muons, and alpha particles. */
  wElectronEv: 33.97,
  wAlphaEv: 35.08,
} as const

/** Ideal-gas density of dry air in g/cm^3. */
export function airDensity(tempK: number, pressurePa: number): number {
  return ((pressurePa * AIR.molarMass) / (R_GAS * tempK)) * 1e-3
}

/** Atoms per cm^3 of each element at the given mass density (g/cm^3). */
export function atomDensities(rho: number): number[] {
  return AIR.elements.map((e) => (rho * e.massFraction * N_A) / e.a)
}

/** Dynamic viscosity of air, Pa s. Sutherland's law (White, Viscous Fluid Flow). */
export function airViscosity(tempK: number): number {
  return (1.716e-5 * (tempK / 273.15) ** 1.5 * (273.15 + 110.4)) / (tempK + 110.4)
}

/** Thermal conductivity of air, W/(m K). Sutherland's law (White, Viscous Fluid Flow). */
export function airConductivity(tempK: number): number {
  return (0.0241 * (tempK / 273.15) ** 1.5 * (273.15 + 194)) / (tempK + 194)
}

/** Mean free path of air molecules, m, from the viscosity (kinetic theory). */
export function airMeanFreePath(tempK: number, pressurePa: number): number {
  const rho = (pressurePa * AIR.molarMass) / (R_GAS * tempK)
  return (airViscosity(tempK) / rho) * Math.sqrt((Math.PI * AIR.molarMass) / (2 * R_GAS * tempK))
}
