/**
 * Isopropanol (2-propanol), the working vapor.
 *
 * Sources, all read 2026-11-12:
 *  - NIST Chemistry WebBook, C67630 phase-change data: normal boiling point 355.5 K, critical
 *    temperature 508.3 K, and the enthalpy-of-vaporization fit (Majer and Svoboda, 298-380 K).
 *    https://webbook.nist.gov/cgi/cbook.cgi?ID=C67630&Mask=4
 *  - PubChem CID 3776 (HSDB): surface tension at 10, 25, 50 and 75 C; density 0.78509 g/cm^3 at
 *    25 C with a correction of 0.00082 per C.
 *  - Dielectric constant 18.23 at 25 C (Wikipedia data page, citing the CRC Handbook).
 *  - Diffusion in air: Fuller, Schettler and Giddings, Ind. Eng. Chem. 58 (1966) 18.
 */
import { ATM, K_B, N_A, R_GAS } from './constants'

export const IPA = {
  molarMass: 0.060096, // kg/mol
  boilingK: 355.5,
  criticalK: 508.3,
  dielectric: 18.23,
} as const

/** Mass of one molecule, kg. */
export const IPA_MOLECULE_KG = IPA.molarMass / N_A

/** Enthalpy of vaporization, J/mol. NIST fit, extrapolated below its 298 K lower limit. */
export function vaporizationEnthalpy(tempK: number): number {
  const tr = tempK / IPA.criticalK
  return 53380 * Math.exp(0.708 * tr) * (1 - tr) ** 0.6538
}

// ln(p_sat / 1 atm) on a 0.5 K grid, by integrating Clausius-Clapeyron (ideal vapor, negligible
// liquid volume) down and up from the normal boiling point.
const GRID_T0 = 150
const GRID_STEP = 0.5
const GRID_N = 521 // up to 410 K
const LN_P = (() => {
  const integrand = (t: number) => vaporizationEnthalpy(t) / (R_GAS * t * t)
  const out = new Float64Array(GRID_N)
  const anchor = (IPA.boilingK - GRID_T0) / GRID_STEP
  const simpson = (a: number, b: number) =>
    ((b - a) / 6) * (integrand(a) + 4 * integrand((a + b) / 2) + integrand(b))
  const at = (i: number) => GRID_T0 + i * GRID_STEP
  const i0 = Math.floor(anchor)
  out[i0] = -simpson(at(i0), IPA.boilingK)
  for (let i = i0 - 1; i >= 0; i--) out[i] = out[i + 1] - simpson(at(i), at(i + 1))
  for (let i = i0 + 1; i < GRID_N; i++) out[i] = out[i - 1] + simpson(at(i - 1), at(i))
  return out
})()

/** Saturation vapor pressure over a flat liquid surface, Pa. */
export function saturationPressure(tempK: number): number {
  const x = Math.min(Math.max((tempK - GRID_T0) / GRID_STEP, 0), GRID_N - 1.000001)
  const i = Math.floor(x)
  const f = x - i
  return ATM * Math.exp(LN_P[i] + f * (LN_P[i + 1] - LN_P[i]))
}

/** Liquid density, kg/m^3. */
export function liquidDensity(tempK: number): number {
  return 785.09 - 0.82 * (tempK - 298.15)
}

/** Volume of one molecule in the liquid, m^3. */
export function molecularVolume(tempK: number): number {
  return IPA_MOLECULE_KG / liquidDensity(tempK)
}

/** Surface tension, N/m. A line through the four HSDB points (22.90 - 0.0789 t mN/m, t in C). */
export function surfaceTension(tempK: number): number {
  return (22.9 - 0.0789 * (tempK - 273.15)) * 1e-3
}

/** Diffusion coefficient of isopropanol vapor in air, m^2/s (Fuller-Schettler-Giddings). */
export function vaporDiffusivity(tempK: number, pressurePa: number): number {
  const volumeAir = 20.1
  const volumeIpa = 3 * 16.5 + 8 * 1.98 + 5.48 // C3H8O atomic diffusion volumes
  const cm2PerS =
    (1e-3 * tempK ** 1.75 * Math.sqrt(1 / 60.096 + 1 / 28.9647)) /
    ((pressurePa / ATM) * (Math.cbrt(volumeAir) + Math.cbrt(volumeIpa)) ** 2)
  return cm2PerS * 1e-4
}

/** Vapor molecules per m^3 at the given partial pressure. */
export function vaporNumberDensity(partialPa: number, tempK: number): number {
  return partialPa / (K_B * tempK)
}
