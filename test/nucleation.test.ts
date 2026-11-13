import { expect, test } from 'bun:test'
import { fogThreshold, ionThreshold, ipaLiquid } from '../src/chamber'
import { K_B, N_A } from '../src/constants'
import {
  homogeneousRate,
  ionBarrier,
  ionBarrierlessSaturation,
  ionNucleationRate,
  kelvinRadius,
  type Liquid,
  thresholdSaturation,
} from '../src/nucleation'
import { rel } from './reference'

// Water at 0 C: surface tension 75.65 mN/m (IAPWS), density 999.84 kg/m^3, permittivity 87.9,
// vapor pressure 611 Pa.
const WATER: Liquid = {
  sigma: 0.07565,
  molecularVolume: 0.018015 / N_A / 999.84,
  molecularMass: 0.018015 / N_A,
  dielectric: 87.9,
  saturationPa: 611,
}

test("classical ion nucleation brackets Wilson's water threshold: the rate criterion is low, the barrierless limit high", () => {
  // C.T.R. Wilson (1897-99) found negative ions start to nucleate water near fourfold
  // supersaturation, S = 4.2 (a textbook value; no source for it is stored in this repo).
  // The simulation uses the rate criterion (one nucleation per ion per second). For water that
  // gives about 3.2, a quarter below Wilson. Thomson's barrierless limit gives about 4.3.
  // So the thresholds this package computes for isopropanol carry an error of that order.
  const rate = thresholdSaturation((s) => ionNucleationRate(WATER, 273.15, s), 1)
  const barrierless = ionBarrierlessSaturation(WATER, 273.15)
  expect(rate).toBeLessThan(4.2)
  expect(barrierless).toBeGreaterThan(4.2)
  expect(rate / 4.2).toBeGreaterThan(0.7)
  expect(rate / 4.2).toBeLessThan(0.82)
  expect(rel(barrierless, 4.2)).toBeLessThan(0.05)
})

test('the sensitive layer exists only because the rate criterion puts the ion threshold below the fog threshold', () => {
  // With the barrierless limit as the ion threshold there would be no window at 243 K: it lies
  // above the fog threshold. That real isopropanol chambers show tracks is the evidence that the
  // ordering used here is the right one; the size of the window is a model result.
  const t = 243
  expect(ionThreshold(t)).toBeLessThan(fogThreshold(t))
  expect(ionBarrierlessSaturation(ipaLiquid(t), t)).toBeGreaterThan(fogThreshold(t))
})

test('the ion barrier vanishes exactly at the barrierless limit', () => {
  const liquid = ipaLiquid(243.15)
  const limit = ionBarrierlessSaturation(liquid, 243.15)
  expect(ionBarrier(liquid, 243.15, limit * 1.001).height).toBe(0)
  expect(ionBarrier(liquid, 243.15, limit * 0.999).height).toBeGreaterThan(0)
  expect(ionBarrier(liquid, 243.15, limit * 0.999).height / (K_B * 243.15)).toBeLessThan(0.01)
})

test('the ion barrier is lower than the homogeneous barrier at the same supersaturation', () => {
  for (const t of [223.15, 243.15, 263.15]) {
    const liquid = ipaLiquid(t)
    const s = 3
    const homogeneous =
      (16 * Math.PI * liquid.sigma ** 3 * liquid.molecularVolume ** 2) /
      (3 * (K_B * t * Math.log(s)) ** 2)
    expect(ionBarrier(liquid, t, s).height).toBeLessThan(homogeneous)
  }
})

test('the Kelvin radius is where the homogeneous free energy peaks', () => {
  const t = 243.15
  const liquid = ipaLiquid(t)
  const s = 3.5
  const dG = (r: number) =>
    (-(4 * Math.PI * r ** 3) / (3 * liquid.molecularVolume)) * K_B * t * Math.log(s) +
    4 * Math.PI * r * r * liquid.sigma
  const r = kelvinRadius(liquid, t, s)
  expect(dG(r)).toBeGreaterThan(dG(r * 0.99))
  expect(dG(r)).toBeGreaterThan(dG(r * 1.01))
})

test('rates rise steeply with supersaturation', () => {
  const t = 243.15
  const liquid = ipaLiquid(t)
  expect(homogeneousRate(liquid, t, 4.2)).toBeGreaterThan(100 * homogeneousRate(liquid, t, 3.8))
  expect(ionNucleationRate(liquid, t, 3.6)).toBeGreaterThan(100 * ionNucleationRate(liquid, t, 3.2))
  expect(homogeneousRate(liquid, t, 1)).toBe(0)
})

test('there is a window where ions nucleate but the vapor does not fog, at every chamber temperature', () => {
  // Without this window a cloud chamber could not show tracks at all.
  let previousIon = Infinity
  for (let t = 205; t <= 290; t += 5) {
    const ion = ionThreshold(t)
    const fog = fogThreshold(t)
    expect(ion).toBeGreaterThan(1.5)
    expect(fog).toBeGreaterThan(ion)
    expect(ion).toBeLessThan(previousIon) // colder vapor needs more supersaturation
    previousIon = ion
  }
})

test('pinned model values: isopropanol thresholds at 243 K', () => {
  // Not a comparison with measurement: these pin the model so a change to a material constant
  // or to the theory shows up as a failing test.
  expect(ionThreshold(243)).toBeCloseTo(3.34, 1)
  expect(fogThreshold(243)).toBeCloseTo(3.98, 1)
})
