import { expect, test } from 'bun:test'
import { AIR, airConductivity } from '../src/air'
import { DEFAULT_CHAMBER, lowerConvexHull, solveChamber } from '../src/chamber'
import { R_GAS } from '../src/constants'
import { IPA, vaporizationEnthalpy } from '../src/ipa'
import { rel } from './reference'

const profile = solveChamber()
const n = profile.z.length
const dz = profile.z[1] - profile.z[0]

test('lower convex hull: a tent under a bump', () => {
  const x = Float64Array.from([0, 1, 2, 3, 4])
  const y = Float64Array.from([0, 5, 1, 5, 0])
  expect([...lowerConvexHull(x, y)]).toEqual([0, 0, 0, 0, 0])
  const y2 = Float64Array.from([0, -2, 5, -2, 0])
  expect([...lowerConvexHull(x, y2)]).toEqual([0, -2, -2, -2, 0])
})

test('vapor is exactly saturated at the plate and at the top, supersaturated in between', () => {
  expect(profile.saturation[0]).toBeCloseTo(1, 9)
  expect(profile.saturation[n - 1]).toBeCloseTo(1, 9)
  for (let i = 1; i < n - 1; i++) expect(profile.saturation[i]).toBeGreaterThan(1)
})

test('supersaturation never exceeds the fog limit', () => {
  for (let i = 0; i < n; i++)
    expect(profile.saturation[i]).toBeLessThanOrEqual(profile.fogThreshold[i] * (1 + 1e-9))
})

test('a sensitive layer emerges above a -40 C plate: centimetres deep, starting millimetres above it', () => {
  expect(profile.sensitiveBottom).toBeGreaterThan(0.002)
  expect(profile.sensitiveBottom).toBeLessThan(0.015)
  expect(profile.sensitiveTop - profile.sensitiveBottom).toBeGreaterThan(0.03)
  expect(profile.sensitiveTop).toBeLessThan(0.75 * DEFAULT_CHAMBER.height)
})

test('a -20 C plate is not cold enough: no sensitive layer, no fog', () => {
  const warm = solveChamber({ ...DEFAULT_CHAMBER, plateK: 253.15 })
  expect(warm.sensitiveTop).toBe(0)
  expect(warm.sink.every((q) => q === 0)).toBe(true)
  const cold = solveChamber({ ...DEFAULT_CHAMBER, plateK: 213.15 })
  expect(cold.sensitiveTop - cold.sensitiveBottom).toBeGreaterThan(
    profile.sensitiveTop - profile.sensitiveBottom,
  )
})

test('with the top at 20 C, the warmest plate that still gives tracks is between -22 and -21 C', () => {
  expect(solveChamber({ ...DEFAULT_CHAMBER, plateK: 251.15 }).sensitiveTop).toBeGreaterThan(0)
  expect(solveChamber({ ...DEFAULT_CHAMBER, plateK: 252.15 }).sensitiveTop).toBe(0)
})

test('the plate limit depends on the top: with the top at 35 C a -20 C plate does give a layer', () => {
  const hot = solveChamber({ ...DEFAULT_CHAMBER, plateK: 253.15, topK: 308.15 })
  expect(hot.sensitiveTop - hot.sensitiveBottom).toBeGreaterThan(0.02)
  expect(solveChamber({ ...DEFAULT_CHAMBER, plateK: 253.15, topK: 293.15 }).sensitiveTop).toBe(0)
})

test('vapor is conserved: what leaves the top either condenses in the fog zone or reaches the plate', () => {
  let condensed = 0
  for (let i = 0; i < n; i++) condensed += profile.sink[i] * dz
  expect(profile.fluxTop).toBeGreaterThan(profile.fluxPlate)
  expect(rel(profile.fluxPlate + condensed, profile.fluxTop)).toBeLessThan(1e-6)
})

test('heat is conserved: conduction into the plate is conduction from the top plus latent heat', () => {
  const flux = (i: number) =>
    (airConductivity(0.5 * (profile.tempK[i] + profile.tempK[i + 1])) *
      (profile.tempK[i + 1] - profile.tempK[i])) /
    dz
  let latent = 0
  for (let i = 0; i < n; i++)
    latent += vaporizationEnthalpy(profile.tempK[i]) * profile.sink[i] * dz
  expect(latent).toBeGreaterThan(0)
  expect(rel(flux(0), flux(n - 2) + latent)).toBeLessThan(1e-4)
})

test('the mist population accounts for the whole condensation sink', () => {
  expect(profile.mistClosure).toBeGreaterThan(0.99)
  expect(profile.mistRate.every((r) => r >= 0)).toBe(true)
})

test('the gas is stably stratified: density falls with height despite the heavy vapor on top', () => {
  // Cold air at the bottom is what keeps the chamber free of convection, so 1D diffusion applies.
  const density = (i: number) => {
    const x = profile.vaporPa[i] / profile.config.pressurePa
    return (
      (profile.config.pressurePa * ((1 - x) * AIR.molarMass + x * IPA.molarMass)) /
      (R_GAS * profile.tempK[i])
    )
  }
  for (let i = 1; i < n; i++) expect(density(i)).toBeLessThan(density(i - 1))
})

test('the solve converges and temperature rises monotonically from plate to top', () => {
  expect(profile.iterations).toBeLessThan(100)
  for (let i = 1; i < n; i++) expect(profile.tempK[i]).toBeGreaterThan(profile.tempK[i - 1])
})
