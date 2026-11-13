import { expect, test } from 'bun:test'
import { MASS_MEV, ME_C2 } from '../src/constants'
import { ASTAR_ELECTRONIC, ASTAR_ENERGY, ASTAR_NUCLEAR } from '../src/data/astar-air'
import { Rng } from '../src/rng'
import {
  alphaStopping,
  betheStopping,
  densityEffect,
  electronDeltaRate,
  electronStopping,
  heavyDeltaRate,
  maxEnergyTransfer,
  sampleElectronDelta,
  sampleHeavyDelta,
} from '../src/stopping'
import { loadTable, rel } from './reference'

const estar = await loadTable('estar-air.txt')
const pstar = await loadTable('pstar-air.txt')
const astar = await loadTable('astar-air.txt')
const muon = (await loadTable('pdg-muon-air.txt')).filter((row) => row.length >= 12)

test('electron collision stopping power matches NIST ESTAR: 0.1% up to 3 MeV, 0.4% up to 100 MeV', () => {
  // Above a few MeV the density effect sets in. ESTAR computes it directly; this code uses the
  // Sternheimer fit, which is where the larger tolerance comes from.
  const low = estar.filter((r) => r[0] >= 0.01 && r[0] <= 3)
  const high = estar.filter((r) => r[0] > 3 && r[0] <= 100)
  expect(low.length).toBeGreaterThan(30)
  expect(high.length).toBeGreaterThan(15)
  for (const [energy, collision] of low)
    expect(rel(electronStopping(energy), collision)).toBeLessThan(0.001)
  for (const [energy, collision] of high)
    expect(rel(electronStopping(energy), collision)).toBeLessThan(0.004)
})

test('muon ionization loss matches the PDG table from 10 MeV to 1 GeV within 0.3%', () => {
  const rows = muon.filter((r) => r[0] >= 10 && r[0] <= 1000)
  expect(rows.length).toBeGreaterThan(20)
  for (const row of rows)
    expect(rel(betheStopping(1, MASS_MEV.muon, row[0]), row[2])).toBeLessThan(0.003)
})

test('minimum ionization in air is 1.815 MeV cm^2/g (PDG), near 260 MeV for a muon', () => {
  let best = Infinity
  let at = 0
  for (let t = 100; t <= 600; t += 1) {
    const s = betheStopping(1, MASS_MEV.muon, t)
    if (s < best) {
      best = s
      at = t
    }
  }
  expect(rel(best, 1.815)).toBeLessThan(0.002)
  expect(Math.abs(at - 257.2)).toBeLessThan(15)
})

test('Bethe for protons matches NIST PSTAR within 1% from 10 MeV up; the shell-correction gap opens below', () => {
  const rows = pstar.filter((r) => r[0] >= 10 && r[0] <= 1000)
  expect(rows.length).toBeGreaterThan(25)
  for (const [energy, electronic] of rows)
    expect(rel(betheStopping(1, MASS_MEV.proton, energy), electronic)).toBeLessThan(0.01)
  // Bethe has no shell corrections, so it runs high for slow protons: 1 to 2% at 5 MeV, more at
  // 1 MeV. This is a systematic, pinned from both sides, and the reason slow alphas use the table.
  const at = (energy: number) => {
    const found = pstar.find((r) => r[0] === energy)
    if (!found) throw new Error('no PSTAR row')
    return betheStopping(1, MASS_MEV.proton, energy) / found[1] - 1
  }
  expect(at(5)).toBeGreaterThan(0.01)
  expect(at(5)).toBeLessThan(0.02)
  expect(at(1)).toBeGreaterThan(at(5))
})

test('Bethe for alphas approaches NIST ASTAR at high energy and fails near the Bragg peak', () => {
  for (const [energy, electronic] of astar.filter((r) => r[0] >= 20 && r[0] <= 1000)) {
    expect(rel(betheStopping(2, MASS_MEV.alpha, energy), electronic)).toBeLessThan(0.02)
  }
  const peak = astar.find((r) => r[0] === 0.7)
  expect(peak).toBeDefined()
  if (peak) expect(rel(betheStopping(2, MASS_MEV.alpha, 0.7), peak[1])).toBeGreaterThan(0.1)
})

test('the shipped alpha table is the fetched NIST ASTAR table, row for row', () => {
  expect(ASTAR_ENERGY.length).toBe(astar.length)
  for (let i = 0; i < astar.length; i++) {
    expect(ASTAR_ENERGY[i]).toBe(astar[i][0])
    expect(ASTAR_ELECTRONIC[i]).toBe(astar[i][1])
    expect(ASTAR_NUCLEAR[i]).toBe(astar[i][2])
  }
})

test('alpha stopping interpolates the table and peaks near 0.7 MeV (the Bragg peak)', () => {
  expect(alphaStopping(5.5)).toBeCloseTo(711.7, 1)
  expect(alphaStopping(0.7)).toBeGreaterThan(alphaStopping(0.3))
  expect(alphaStopping(0.7)).toBeGreaterThan(alphaStopping(1.5))
  // Between two rows it stays between their values.
  const mid = alphaStopping(5.25)
  expect(mid).toBeLessThan(alphaStopping(5.0))
  expect(mid).toBeGreaterThan(alphaStopping(5.5))
})

test('the density effect is zero in air below beta gamma = 55 and grows above', () => {
  expect(densityEffect(10)).toBe(0)
  expect(densityEffect(1000)).toBeGreaterThan(densityEffect(100))
})

test('maximum energy transfer: 2 m c^2 beta^2 gamma^2 for a heavy particle', () => {
  // 5.5 MeV alpha: beta^2 = 0.00295, so about 3.0 keV.
  expect(
    rel(maxEnergyTransfer(5.5, MASS_MEV.alpha), 2 * ME_C2 * ((1 + 5.5 / MASS_MEV.alpha) ** 2 - 1)),
  ).toBeLessThan(0.001)
})

test('restricted loss plus the energy carried by delta rays is the full Bethe loss (muon)', () => {
  const t = 1000
  const cut = 0.001
  const rng = new Rng(11)
  const rate = heavyDeltaRate(1, MASS_MEV.muon, t, cut)
  let mean = 0
  const n = 400000
  for (let i = 0; i < n; i++)
    mean += sampleHeavyDelta(t, MASS_MEV.muon, cut, rng.next(), () => rng.next()) / n
  // The 1/T^2 spectrum is heavy-tailed, so the sampled mean converges slowly: allow 4%.
  const total = betheStopping(1, MASS_MEV.muon, t, cut) + rate * mean
  expect(rel(total, betheStopping(1, MASS_MEV.muon, t))).toBeLessThan(0.04)
})

test('restricted loss plus Moller delta rays is the full electron loss', () => {
  const t = 1
  const cut = 0.001
  const rng = new Rng(12)
  const rate = electronDeltaRate(t, cut)
  let mean = 0
  const n = 200000
  for (let i = 0; i < n; i++) mean += sampleElectronDelta(t, cut, () => rng.next()) / n
  expect(rel(electronStopping(t, cut) + rate * mean, electronStopping(t))).toBeLessThan(0.02)
})

test('a 1 keV cut leaves no delta rays for alphas below ~1.8 MeV and many for fast muons', () => {
  expect(heavyDeltaRate(2, MASS_MEV.alpha, 1.5, 0.001)).toBe(0)
  expect(heavyDeltaRate(2, MASS_MEV.alpha, 5.5, 0.001)).toBeGreaterThan(0)
  expect(heavyDeltaRate(1, MASS_MEV.muon, 3000, 0.001)).toBeGreaterThan(50)
})
