/**
 * Physical constants. SI units unless a name says otherwise; particle energies are in MeV.
 * Values are CODATA 2018 (the SI-2019 defining constants are exact).
 */
export const K_B = 1.380649e-23 // J/K
export const N_A = 6.02214076e23 // 1/mol
export const R_GAS = K_B * N_A // J/(mol K)
export const E_CHARGE = 1.602176634e-19 // C
export const EPS0 = 8.8541878128e-12 // F/m
export const G_ACCEL = 9.80665 // m/s^2
export const ATM = 101325 // Pa

export const ME_C2 = 0.51099895 // MeV, electron rest energy
export const R_E_CM = 2.8179403262e-13 // cm, classical electron radius
export const ALPHA_FS = 7.2973525693e-3 // fine-structure constant
export const HBAR_C_MEV_CM = 197.3269804e-13 // MeV cm
export const A_BOHR_CM = 5.29177210903e-9 // cm
export const AMU_MEV = 931.49410242 // MeV per atomic mass unit

/** 4 pi N_A r_e^2 m_e c^2, the Bethe prefactor, in MeV cm^2 / mol (0.307075). */
export const K_BETHE = 4 * Math.PI * N_A * R_E_CM * R_E_CM * ME_C2

export const MASS_MEV = {
  electron: ME_C2,
  muon: 105.6583755,
  proton: 938.27208816,
  alpha: 3727.3794066,
} as const

export const EV = 1e-6 // MeV per eV
