import { expect, test } from 'bun:test'
import { airConductivity, airDensity, airMeanFreePath, airViscosity } from '../src/air'
import { ATM } from '../src/constants'
import { fallSpeed, growthCoefficient, stokesCoefficient } from '../src/droplet'
import {
  liquidDensity,
  saturationPressure,
  surfaceTension,
  vaporDiffusivity,
  vaporizationEnthalpy,
} from '../src/ipa'
import { rel } from './reference'

const MMHG = 133.322
const C = 273.15

test('vapor pressure against the CRC table: within 8% from 2.4 to 82.5 C, and 8 to 12% low at -26.1 C', () => {
  // The model is Clausius-Clapeyron integrated from the boiling point with NIST's enthalpy fit.
  // None of these points is an input. The fit is outside its range below 25 C and the error
  // grows as it cools: the coldest point is pinned from both sides as a known bias, and the
  // chamber's plate is colder still, where nothing here checks it.
  const table: [number, number][] = [
    [2.4, 10],
    [23.8, 40],
    [39.5, 100],
    [67.8, 400],
    [82.5, 760],
  ]
  for (const [tc, mmHg] of table)
    expect(rel(saturationPressure(tc + C) / MMHG, mmHg)).toBeLessThan(0.08)
  const cold = saturationPressure(-26.1 + C) / MMHG - 1
  expect(cold).toBeLessThan(-0.08)
  expect(cold).toBeGreaterThan(-0.12)
})

test('vapor pressure at room temperature matches HSDB (45.4 mmHg at 25 C; 4.4 kPa at 20 C)', () => {
  expect(rel(saturationPressure(25 + C) / MMHG, 45.4)).toBeLessThan(0.03)
  expect(rel(saturationPressure(20 + C), 4400)).toBeLessThan(0.04)
})

test('enthalpy of vaporization is 45.39 kJ/mol at 25 C (NIST, HSDB)', () => {
  expect(rel(vaporizationEnthalpy(298.15), 45390)).toBeLessThan(0.001)
})

test('surface tension passes through the four HSDB points', () => {
  const points: [number, number][] = [
    [10, 22.11],
    [25, 20.93],
    [50, 18.96],
    [75, 16.98],
  ]
  for (const [tc, mN] of points) expect(rel(surfaceTension(tc + C) * 1e3, mN)).toBeLessThan(0.002)
})

test('liquid density is 0.78509 g/cm^3 at 25 C and rises on cooling', () => {
  expect(liquidDensity(298.15)).toBeCloseTo(785.09, 6)
  expect(liquidDensity(233.15)).toBeGreaterThan(liquidDensity(298.15))
})

test('air density at 20 C and 1 atm is the PDG value, 1.205e-3 g/cm^3', () => {
  expect(rel(airDensity(293.15, ATM), 1.205e-3)).toBeLessThan(0.002)
})

test('diffusivity scales as T^1.75 / P (Fuller-Schettler-Giddings)', () => {
  expect(vaporDiffusivity(300, ATM) / vaporDiffusivity(250, ATM)).toBeCloseTo(
    (300 / 250) ** 1.75,
    10,
  )
  expect(vaporDiffusivity(300, ATM / 2)).toBeCloseTo(2 * vaporDiffusivity(300, ATM), 12)
})

test('Sutherland laws give positive, rising transport coefficients; mean free path near 65 nm at 20 C', () => {
  expect(airViscosity(293.15)).toBeGreaterThan(airViscosity(233.15))
  expect(airConductivity(293.15)).toBeGreaterThan(airConductivity(233.15))
  const mfp = airMeanFreePath(293.15, ATM)
  expect(mfp).toBeGreaterThan(60e-9)
  expect(mfp).toBeLessThan(70e-9)
})

test('fall speed is Stokes for large droplets and exceeds it for small ones (slip)', () => {
  const t = 243.15
  expect(rel(fallSpeed(50e-6, t, ATM), stokesCoefficient(t) * 50e-6 ** 2)).toBeLessThan(0.003)
  expect(fallSpeed(0.2e-6, t, ATM)).toBeGreaterThan(1.3 * stokesCoefficient(t) * 0.2e-6 ** 2)
})

test('droplet growth is slower than vapor diffusion alone would allow (latent heat must conduct away)', () => {
  const t = 243.15
  const diffusionOnly =
    (vaporDiffusivity(t, ATM) * saturationPressure(t) * 0.060096) /
    (8.314462618 * t * liquidDensity(t))
  const g = growthCoefficient(t, ATM)
  expect(g).toBeLessThan(diffusionOnly)
  expect(g).toBeGreaterThan(0.2 * diffusionOnly)
})
