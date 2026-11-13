import { expect, test } from 'bun:test'
import { ME_C2 } from '../src/constants'
import { Rng } from '../src/rng'
import {
  ambientElectrons,
  americiumSource,
  cosmicMuons,
  enterBox,
  K40_GAMMA,
  MUON_VERTICAL_INTENSITY,
  radonAlphas,
  sampleBeta,
  sampleCompton,
} from '../src/sources'
import { loadTable, rel } from './reference'

const box = { halfX: 0.1, halfY: 0.1, height: 0.12 }

test('muons cross the chamber top at (pi/2) I_v per unit area, about 0.66 per cm^2 per minute', () => {
  const source = cosmicMuons(box)
  const rng = new Rng(1)
  const n = 400000
  let throughTop = 0
  let cosSum = 0
  for (let i = 0; i < n; i++) {
    const p = source.sample(rng)
    if (Math.abs(p.pos[0]) <= box.halfX && Math.abs(p.pos[1]) <= box.halfY) {
      throughTop++
      cosSum += -p.dir[2]
    }
  }
  const perSecond = (source.rate * throughTop) / n
  const expected = (Math.PI / 2) * MUON_VERTICAL_INTENSITY * 0.04
  expect(rel(perSecond, expected)).toBeLessThan(0.03)
  expect((expected / 400) * 60).toBeCloseTo(0.66, 1)
  // Flux through a horizontal plane goes as cos^3 sin, whose mean cosine is 4/5.
  expect(cosSum / throughTop).toBeCloseTo(0.8, 2)
})

test('enterBox moves an outside particle to the wall it crosses and rejects misses', () => {
  const hit = enterBox(
    { species: 'muon', kinetic: 1, pos: [0.2, 0, 0.12], dir: [-Math.SQRT1_2, 0, -Math.SQRT1_2] },
    box,
  )
  expect(hit).not.toBeNull()
  if (hit) {
    expect(hit.pos[0]).toBeCloseTo(0.1, 6)
    expect(hit.pos[2]).toBeCloseTo(0.02, 6)
  }
  // The same ray started further out passes under the plate.
  expect(
    enterBox(
      { species: 'muon', kinetic: 1, pos: [0.3, 0, 0.12], dir: [-Math.SQRT1_2, 0, -Math.SQRT1_2] },
      box,
    ),
  ).toBeNull()
  expect(
    enterBox({ species: 'muon', kinetic: 1, pos: [0.5, 0, 0.12], dir: [0, 0, -1] }, box),
  ).toBeNull()
})

test('radon: 40 Bq/m^3 in a 4.8 litre chamber gives 0.19 Rn decays/s and 1.8 alphas per decay', () => {
  const source = radonAlphas(box, 40, 0.4)
  expect(source.rate).toBeCloseTo(40 * 0.0048 * 1.8, 9)
  const rng = new Rng(2)
  const energies = new Set<number>()
  for (let i = 0; i < 500; i++) energies.add(source.sample(rng).kinetic)
  expect([...energies].sort()).toEqual([5.4895, 6.0024, 7.6868])
})

test('the Am-241 source emits inside its cone', () => {
  const source = americiumSource([0, 0, 0.03], [1, 0, 0], 4, 0.3)
  const rng = new Rng(3)
  for (let i = 0; i < 500; i++) {
    const p = source.sample(rng)
    expect(Math.hypot(...p.dir)).toBeCloseTo(1, 12)
    expect(p.dir[0]).toBeGreaterThanOrEqual(Math.cos(0.3) - 1e-12)
    expect([5.4856, 5.4428]).toContain(p.kinetic)
  }
})

test('beta spectra: mean energies are 3 to 6% above the evaluated ones (IAEA/ENSDF), a known bias', async () => {
  // The spectrum shape uses a non-relativistic Fermi function, which under-weights slow
  // electrons. The bias is systematic, so it is pinned from both sides: a correct relativistic
  // Fermi function would bring it under 3% and should then replace this band.
  const rows = await loadTable('iaea-beta.txt')
  expect(rows.length).toBe(2)
  for (const [zDaughter, endpointKeV, meanKeV] of rows) {
    const rng = new Rng(4)
    const n = 60000
    let mean = 0
    for (let i = 0; i < n; i++) mean += sampleBeta(endpointKeV / 1000, zDaughter, rng) / n
    const excess = (mean * 1000) / meanKeV - 1
    expect(excess).toBeGreaterThan(0.03)
    expect(excess).toBeLessThan(0.06)
  }
})

test('Compton electrons: never above the Compton edge, mean energy equal to a direct Klein-Nishina integral', () => {
  const k = K40_GAMMA / ME_C2
  const edge = (K40_GAMMA * 2 * k) / (1 + 2 * k)
  const rng = new Rng(5)
  const n = 150000
  let mean = 0
  let max = 0
  for (let i = 0; i < n; i++) {
    const s = sampleCompton(K40_GAMMA, rng)
    mean += s.kinetic / n
    max = Math.max(max, s.kinetic)
    expect(s.cos).toBeGreaterThanOrEqual(0) // the electron always goes forward
  }
  expect(max).toBeLessThanOrEqual(edge)
  expect(max).toBeGreaterThan(0.99 * edge)
  // Independent quadrature of the Klein-Nishina cross-section over the photon angle.
  let num = 0
  let den = 0
  for (let i = 0; i < 20000; i++) {
    const cos = -1 + (2 * (i + 0.5)) / 20000
    const ratio = 1 / (1 + k * (1 - cos)) // scattered / incident photon energy
    const weight = ratio * ratio * (ratio + 1 / ratio - (1 - cos * cos))
    num += weight * K40_GAMMA * (1 - ratio)
    den += weight
  }
  expect(rel(mean, num / den)).toBeLessThan(0.01)
})

test('background electron rate follows from the dose rate: 60 nGy/h in this chamber is about one per second', () => {
  const source = ambientElectrons(box, 263.15, 101325, 60)
  expect(source.rate).toBeGreaterThan(0.5)
  expect(source.rate).toBeLessThan(2)
  expect(ambientElectrons(box, 263.15, 101325, 120).rate).toBeCloseTo(2 * source.rate, 9)
})
