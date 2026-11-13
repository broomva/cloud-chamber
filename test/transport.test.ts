import { expect, test } from 'bun:test'
import { AIR, airDensity } from '../src/air'
import { solveChamber } from '../src/chamber'
import { ATM, MASS_MEV, ME_C2 } from '../src/constants'
import { Rng } from '../src/rng'
import { deflect, gasScattering, sampleMu } from '../src/scattering'
import { betheStopping, heavyDeltaRate, kinematics } from '../src/stopping'
import { DEFAULT_TRANSPORT, PROJECTILES, type Segment, transport } from '../src/transport'
import { loadTable, rel } from './reference'

// A large box of uniform air at 20 C and 1 atm, the conditions of the NIST and PDG tables.
const uniform = solveChamber({
  height: 40,
  plateK: 293.15,
  topK: 293.15,
  pressurePa: ATM,
  cells: 10,
})
const wide = { ...DEFAULT_TRANSPORT, halfX: 20, halfY: 20 }
const rho = airDensity(293.15, ATM)
const astar = await loadTable('astar-air.txt')
const estar = await loadTable('estar-air.txt')
const row = (table: number[][], energy: number) => {
  const found = table.find((r) => r[0] === energy)
  if (!found) throw new Error(`no row for ${energy} MeV`)
  return found
}

function fireMany(
  species: 'alpha' | 'electron' | 'muon',
  kinetic: number,
  n: number,
  seed: number,
  options = wide,
) {
  const rng = new Rng(seed)
  const out = { displacement: 0, path: 0, deposited: 0, nuclear: 0, deltas: 0, spread: 0 }
  const ds: number[] = []
  for (let i = 0; i < n; i++) {
    const r = transport(
      { species, kinetic, pos: [0, 0, 20], dir: [1, 0, 0] },
      uniform,
      rng,
      () => {},
      options,
    )
    const d = Math.hypot(r.end[0], r.end[1], r.end[2] - 20)
    ds.push(d)
    out.displacement += d / n
    out.path += r.pathLength / n
    out.deposited += r.depositedMeV / n
    out.nuclear += r.nuclearMeV / n
    out.deltas += r.deltaRays / n
  }
  out.spread = Math.sqrt(ds.reduce((s, d) => s + (d - out.displacement) ** 2, 0) / n)
  return out
}

test('alpha range: simulated 5.5 MeV and 2 MeV alphas stop at the NIST ASTAR projected range (1.5%)', () => {
  for (const energy of [5.5, 2]) {
    const sim = fireMany('alpha', energy, 120, 7)
    const projected = row(astar, energy)[5] / rho / 100 // g/cm^2 -> m
    expect(rel(sim.displacement, projected)).toBeLessThan(0.015)
  }
})

test('alpha range straggling is about 1% and the path is nearly straight', () => {
  const sim = fireMany('alpha', 5.5, 200, 8)
  // Measured straggling for 5 MeV alphas in air is about 1% of the range. The band is wide
  // because 200 tracks fix a standard deviation only to about 5%, and a rare nuclear kink in one
  // track can add more.
  expect(sim.spread / sim.displacement).toBeGreaterThan(0.004)
  expect(sim.spread / sim.displacement).toBeLessThan(0.03)
  expect(sim.displacement / sim.path).toBeGreaterThan(0.99) // ASTAR detour factor: 0.9925
})

test('energy is conserved: a stopped alpha turns all 5.5 MeV into ionization', () => {
  const sim = fireMany('alpha', 5.5, 20, 9)
  expect(sim.deposited).toBeCloseTo(5.5, 9)
  // A few keV goes to recoiling nuclei: nuclear stopping, which comes out of the collisions.
  expect(sim.nuclear).toBeGreaterThan(0.002)
  expect(sim.nuclear).toBeLessThan(0.03)
})

test('electron path length matches the NIST ESTAR CSDA range at 100 keV and 20 keV (3%)', () => {
  for (const energy of [0.1, 0.02]) {
    const sim = fireMany('electron', energy, 250, 3)
    const csda = row(estar, energy)[4] / rho / 100
    expect(rel(sim.path, csda)).toBeLessThan(0.03)
    // And unlike an alpha it wanders: it ends well short of its path length.
    expect(sim.displacement).toBeLessThan(0.75 * sim.path)
  }
})

test('with no delta rays split off, a muon deposits exactly the Bethe loss along its path', () => {
  const rng = new Rng(4)
  const r = transport(
    { species: 'muon', kinetic: 1000, pos: [-19.9, 0, 20], dir: [1, 0, 0] },
    uniform,
    rng,
    () => {},
    { ...wide, cutMeV: 1e9 },
  )
  expect(r.deltaRays).toBe(0)
  const perGram = r.depositedMeV / (r.pathLength * 100 * rho)
  expect(rel(perGram, betheStopping(1, MASS_MEV.muon, 1000))).toBeLessThan(0.002)
})

test('delta rays appear at the rate the cross-section predicts', () => {
  const rng = new Rng(5)
  let deltas = 0
  let grams = 0
  for (let i = 0; i < 40; i++) {
    const r = transport(
      { species: 'muon', kinetic: 1000, pos: [-19.9, 0, 20], dir: [1, 0, 0] },
      uniform,
      rng,
      () => {},
      { ...wide, followDeltas: false },
    )
    deltas += r.deltaRays
    grams += r.pathLength * 100 * rho
  }
  const expected = heavyDeltaRate(1, MASS_MEV.muon, 1000, wide.cutMeV) * grams
  expect(Math.abs(deltas - expected)).toBeLessThan(4 * Math.sqrt(expected))
})

test('a minimum-ionizing muon leaves 60 to 75 ion pairs per cm in air at 1 atm', () => {
  // dE/dx_min * density / W = 1.815e6 eV cm^2/g * 1.205e-3 g/cm^3 / 33.97 eV = 64 per cm.
  const rng = new Rng(6)
  let pairs = 0
  let metres = 0
  for (let i = 0; i < 20; i++) {
    const r = transport(
      { species: 'muon', kinetic: 300, pos: [-19.9, 0, 20], dir: [1, 0, 0] },
      uniform,
      rng,
      (s: Segment) => {
        pairs += s.ionPairs
      },
      wide,
    )
    metres += r.pathLength
  }
  const perCm = pairs / (metres * 100)
  expect(perCm).toBeGreaterThan(55)
  expect(perCm).toBeLessThan(75)
}, 20000)

/** The rotation as first written: near the z axis it replaced the direction instead of turning it. */
function buggyDeflect(
  d: [number, number, number],
  cosTheta: number,
  phi: number,
): [number, number, number] {
  if (Math.abs(d[2]) > 0.99999) {
    const sinTheta = Math.sqrt(Math.max(0, 1 - cosTheta * cosTheta))
    return [sinTheta * Math.cos(phi), sinTheta * Math.sin(phi), d[2] > 0 ? cosTheta : -cosTheta]
  }
  return deflect(d, cosTheta, phi)
}

/** Width of the angular distribution after 1% of a radiation length, over the Highland value. */
function scatteringWidthOverHighland(rotate: typeof deflect): number {
  // 100 MeV electrons through 1% of a radiation length of air, about 3100 collisions each.
  const kinetic = 100
  const thickness = 0.01 * AIR.radiationLength
  const gas = gasScattering(PROJECTILES.electron, kinetic)
  const kin = kinematics(kinetic, ME_C2)
  const rng = new Rng(5)
  const angles: number[] = []
  for (let i = 0; i < 2500; i++) {
    let d: [number, number, number] = [0, 0, 1]
    const collisions = rng.poisson(gas.perGram * thickness)
    for (let k = 0; k < collisions; k++) {
      let pick = rng.next()
      let el = gas.elements[gas.elements.length - 1]
      for (const candidate of gas.elements) {
        if (pick < candidate.share) {
          el = candidate
          break
        }
        pick -= candidate.share
      }
      d = rotate(d, 1 - 2 * sampleMu(el.eta, rng.next()), 2 * Math.PI * rng.next())
    }
    angles.push(Math.atan2(d[0], d[2]))
  }
  angles.sort((a, b) => a - b)
  const core = angles.slice(Math.floor(0.01 * angles.length), Math.ceil(0.99 * angles.length))
  const rms = Math.sqrt(core.reduce((s, v) => s + v * v, 0) / core.length)
  // The rms of a Gaussian cut to its central 98% is 0.9346 of its sigma.
  const theta0 = rms / 0.9346
  const highland =
    (13.6 / (Math.sqrt(kin.beta2) * kin.p)) *
    Math.sqrt(0.01) *
    (1 + 0.038 * Math.log(0.01 / kin.beta2))
  return theta0 / highland
}

test('multiple scattering: many single collisions add up to the Highland width (PDG) within its 11% accuracy', () => {
  expect(Math.abs(scatteringWidthOverHighland(deflect) - 1)).toBeLessThan(0.11)
})

test('this check is what caught the near-axis rotation bug: with the old rotation the width is a quarter low', () => {
  const ratio = scatteringWidthOverHighland(buggyDeflect)
  expect(ratio).toBeLessThan(0.8)
  expect(ratio).toBeGreaterThan(0.6)
})

test('deflect keeps unit length and turns by exactly the requested angle, including near the z axis', () => {
  const rng = new Rng(1)
  for (const start of [
    [0, 0, 1],
    [1e-6, 0, 1],
    [0.001, -0.002, 1],
    [1, 2, 3],
    [0, 0, -1],
  ] as [number, number, number][]) {
    const norm = Math.hypot(...start)
    const d: [number, number, number] = [start[0] / norm, start[1] / norm, start[2] / norm]
    for (let i = 0; i < 50; i++) {
      const cos = 2 * rng.next() - 1
      const out = deflect(d, cos, 2 * Math.PI * rng.next())
      expect(Math.hypot(...out)).toBeCloseTo(1, 12)
      expect(out[0] * d[0] + out[1] * d[1] + out[2] * d[2]).toBeCloseTo(cos, 9)
    }
  }
})

test('a particle stops at the wall it reaches', () => {
  const small = solveChamber({
    height: 0.1,
    plateK: 293.15,
    topK: 293.15,
    pressurePa: ATM,
    cells: 10,
  })
  const rng = new Rng(2)
  const r = transport(
    { species: 'muon', kinetic: 3000, pos: [0, 0, 0.05], dir: [0.6, 0, 0.8] },
    small,
    rng,
    () => {},
  )
  expect(r.stopped).toBe(false)
  expect(r.end[2]).toBeCloseTo(0.1, 9)
  expect(r.end[0]).toBeCloseTo(0.0375, 4)
})
