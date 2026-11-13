/**
 * Prints the measured numbers quoted in README.md, so each one can be regenerated:
 *
 *   bun run tools/report.ts
 */
import { join } from 'node:path'
import { AIR, airDensity, airViscosity } from '../src/air'
import {
  DEFAULT_CHAMBER,
  fogThreshold,
  ionThreshold,
  ipaLiquid,
  sampleProfile,
  solveChamber,
} from '../src/chamber'
import { ATM, MASS_MEV, ME_C2, N_A } from '../src/constants'
import { Engine, KIND, TRACK_RADIUS } from '../src/engine'
import { ionBarrierlessSaturation, ionNucleationRate, thresholdSaturation } from '../src/nucleation'
import { topViewScale } from '../src/render'
import { Rng } from '../src/rng'
import { deflect, gasScattering, sampleMu } from '../src/scattering'
import { betheStopping, kinematics } from '../src/stopping'
import { DEFAULT_TRANSPORT, PROJECTILES, transport } from '../src/transport'

async function table(name: string): Promise<number[][]> {
  const text = await Bun.file(join(import.meta.dir, '..', 'test', 'reference', name)).text()
  return text
    .split('\n')
    .filter((l) => l.trim() && !l.startsWith('#'))
    .map((l) => l.trim().split(/\s+/).map(Number))
}
const row = (rows: number[][], energy: number) => rows.find((r) => r[0] === energy) ?? []
const astar = await table('astar-air.txt')
const estar = await table('estar-air.txt')
const line = (label: string, value: string) => console.log(`${label.padEnd(58)} ${value}`)

const uniform = solveChamber({
  height: 40,
  plateK: 293.15,
  topK: 293.15,
  pressurePa: ATM,
  cells: 10,
})
const wide = { ...DEFAULT_TRANSPORT, halfX: 20, halfY: 20 }
const rho = airDensity(293.15, ATM)

console.log('Chamber (12 cm tall, plate -40 C, top 20 C)')
const p = solveChamber(DEFAULT_CHAMBER)
line(
  'sensitive layer, mm above the plate',
  `${(p.sensitiveBottom * 1e3).toFixed(1)} to ${(p.sensitiveTop * 1e3).toFixed(1)}`,
)
let mist = 0
let fogLo = Infinity
let fogHi = 0
for (let i = 0; i < p.z.length; i++) {
  mist += p.mistRate[i] * (p.z[1] - p.z[0])
  if (p.sink[i] > 0) {
    fogLo = Math.min(fogLo, p.z[i])
    fogHi = Math.max(fogHi, p.z[i])
  }
}
line(
  'fog zone (vapor at the fog limit), mm',
  `${(fogLo * 1e3).toFixed(0)} to ${(fogHi * 1e3).toFixed(0)}`,
)
line('mist droplets nucleated per cm^2 of plate per second', (mist * 1e-4).toFixed(0))
line(
  'vapor leaving the top that condenses before the plate',
  `${((1 - p.fluxPlate / p.fluxTop) * 100).toFixed(0)}%`,
)
line(
  'ion / fog threshold at 243 K',
  `${ionThreshold(243).toFixed(2)} / ${fogThreshold(243).toFixed(2)}`,
)
line(
  'window between them, % of the fog threshold',
  ((1 - ionThreshold(243) / fogThreshold(243)) * 100).toFixed(0),
)
let lo = 233.15
let hi = 253.15
for (let i = 0; i < 30; i++) {
  const mid = (lo + hi) / 2
  if (solveChamber({ ...DEFAULT_CHAMBER, plateK: mid }).sensitiveTop > 0) lo = mid
  else hi = mid
}
line('warmest plate that still gives a sensitive layer, top 20 C', (lo - 273.15).toFixed(1))
for (const topC of [10, 35]) {
  let cold = 213.15
  let warm = 273.15
  for (let i = 0; i < 30; i++) {
    const mid = (cold + warm) / 2
    if (solveChamber({ ...DEFAULT_CHAMBER, plateK: mid, topK: topC + 273.15 }).sensitiveTop > 0)
      cold = mid
    else warm = mid
  }
  line(
    `warmest plate that still gives a sensitive layer, top ${topC} C`,
    (cold - 273.15).toFixed(1),
  )
}

console.log('\nTransport in air at 20 C, 1 atm')
for (const energy of [5.5, 2]) {
  const rng = new Rng(7)
  let d = 0
  const n = 120
  for (let i = 0; i < n; i++) {
    const r = transport(
      { species: 'alpha', kinetic: energy, pos: [0, 0, 20], dir: [1, 0, 0] },
      uniform,
      rng,
      () => {},
      wide,
    )
    d += Math.hypot(r.end[0], r.end[1], r.end[2] - 20) / n
  }
  line(
    `alpha ${energy} MeV: simulated range / NIST ASTAR projected, cm`,
    `${(d * 100).toFixed(3)} / ${(row(astar, energy)[5] / rho).toFixed(3)}`,
  )
}
for (const energy of [0.1, 0.02]) {
  const rng = new Rng(3)
  let path = 0
  const n = 250
  for (let i = 0; i < n; i++)
    path +=
      transport(
        { species: 'electron', kinetic: energy, pos: [0, 0, 20], dir: [1, 0, 0] },
        uniform,
        rng,
        () => {},
        wide,
      ).pathLength / n
  line(
    `electron ${energy * 1000} keV: simulated path / NIST ESTAR CSDA, cm`,
    `${(path * 100).toFixed(2)} / ${(row(estar, energy)[4] / rho).toFixed(2)}`,
  )
}
{
  const gas = gasScattering(PROJECTILES.electron, 100)
  const kin = kinematics(100, ME_C2)
  const rng = new Rng(5)
  const angles: number[] = []
  for (let i = 0; i < 2500; i++) {
    let d: [number, number, number] = [0, 0, 1]
    const collisions = rng.poisson(gas.perGram * 0.01 * AIR.radiationLength)
    for (let k = 0; k < collisions; k++) {
      let pick = rng.next()
      let el = gas.elements[gas.elements.length - 1]
      for (const c of gas.elements) {
        if (pick < c.share) {
          el = c
          break
        }
        pick -= c.share
      }
      d = deflect(d, 1 - 2 * sampleMu(el.eta, rng.next()), 2 * Math.PI * rng.next())
    }
    angles.push(Math.atan2(d[0], d[2]))
  }
  angles.sort((a, b) => a - b)
  const core = angles.slice(25, 2475)
  const theta0 = Math.sqrt(core.reduce((s, v) => s + v * v, 0) / core.length) / 0.9346
  const highland =
    (13.6 / (Math.sqrt(kin.beta2) * kin.p)) * 0.1 * (1 + 0.038 * Math.log(0.01 / kin.beta2))
  line(
    'scattering width, simulated / Highland (PDG), mrad',
    `${(theta0 * 1e3).toFixed(2)} / ${(highland * 1e3).toFixed(2)}`,
  )
}

console.log('\nTracks in the chamber, fired level at 3 cm')
const quiet = { background: false, mist: false }
for (const [species, kinetic, x0] of [
  ['alpha', 5.5, -0.05],
  ['muon', 3000, -0.0999],
] as const) {
  // Seconds until 90% of the real droplets have landed: the median of five seeds, because one
  // hard knock-on electron high in the chamber can hold a single run up.
  const landed = (collectiveSettling: boolean) => {
    const times: number[] = []
    for (let seed = 1; seed <= 5; seed++) {
      const e = new Engine({ seed, ...quiet, collectiveSettling })
      e.fire({ species, kinetic, pos: [x0, 0, 0.03], dir: [1, 0, 0] })
      let t = 60
      for (let i = 1; i <= 3600; i++) {
        e.step(1 / 60)
        if (e.landed >= 0.9 * e.created) {
          t = i / 60
          break
        }
      }
      times.push(t)
    }
    times.sort((a, b) => a - b)
    return times[2].toFixed(1)
  }
  // Largest Reynolds number of the sinking column, on the track diameter.
  const reynolds = () => {
    const e = new Engine({ seed: 5, ...quiet })
    e.fire({ species, kinetic, pos: [x0, 0, 0.03], dir: [1, 0, 0] })
    let peak = 0
    const before = new Float32Array(e.options.capacity)
    for (let i = 0; i < 180 && e.count > 0; i++) {
      for (let k = 0; k < e.count; k++) before[k] = e.pos[3 * k + 2]
      const n = e.count
      e.step(1 / 60)
      if (e.count !== n) continue // entries were removed and the arrays compacted: skip this step
      const speeds: number[] = []
      for (let k = 0; k < n; k++)
        if (e.radius[k] > 0) speeds.push((before[k] - e.pos[3 * k + 2]) * 60)
      speeds.sort((a, b) => a - b)
      const typical = speeds[Math.floor(speeds.length / 2)] ?? 0
      const t = sampleProfile(e.profile, e.profile.tempK, 0.03)
      peak = Math.max(
        peak,
        (typical * 2 * TRACK_RADIUS) / (airViscosity(t) / (airDensity(t, ATM) * 1e3)),
      )
    }
    return peak.toFixed(1)
  }
  const e = new Engine({ seed: 5, ...quiet })
  e.fire({ species, kinetic, pos: [x0, 0, 0.03], dir: [1, 0, 0] })
  let real = 0
  for (let k = 0; k < e.count; k++) real += e.weight[k]
  const entries = e.count
  for (let i = 0; i < 12; i++) e.step(1 / 60)
  let lo = Infinity
  let hi = -Infinity
  for (let k = 0; k < e.count; k++) {
    lo = Math.min(lo, e.pos[3 * k])
    hi = Math.max(hi, e.pos[3 * k])
  }
  for (let i = 0; i < 6; i++) e.step(1 / 60)
  let radius = 0
  let n = 0
  for (let k = 0; k < e.count; k++)
    if (e.kind[k] === KIND[species] && e.radius[k] > 0) {
      radius += e.radius[k]
      n++
    }
  line(
    `${species}: real droplets / array entries / droplets per entry`,
    `${Math.round(real)} / ${entries} / ${(real / entries).toFixed(1)}`,
  )
  line(`${species}: track length along x at 0.2 s, cm`, ((hi - lo) * 100).toFixed(2))
  line(`${species}: mean droplet radius at 0.3 s, micrometres`, ((radius / n) * 1e6).toFixed(1))
  line(
    `${species}: seconds until 90% of droplets land, with / without downdraft`,
    `${landed(true)} / ${landed(false)}`,
  )
  line(
    `${species}: peak Reynolds number of the column (median droplet, track diameter)`,
    reynolds(),
  )
}

console.log('\nHow tracks fade as the plate warms (5.5 MeV alpha at 3 cm, top at 20 C)')
for (const plateC of [-40, -25, -22, -21, -20, -18, -15]) {
  const e = new Engine({
    seed: 5,
    ...quiet,
    chamber: { ...DEFAULT_CHAMBER, plateK: plateC + 273.15 },
  })
  e.fire({ species: 'alpha', kinetic: 5.5, pos: [-0.05, 0, 0.03], dir: [1, 0, 0] })
  let real = 0
  for (let k = 0; k < e.count; k++) real += e.weight[k]
  line(
    `plate ${plateC} C: share of ions that seed a droplet`,
    `${((100 * real) / (2 * e.ionPairs)).toFixed(real / (2 * e.ionPairs) < 0.01 ? 3 : 1)}%`,
  )
}

console.log('\nIon nucleation theory against water (273.15 K)')
const water = {
  sigma: 0.07565,
  molecularVolume: 0.018015 / N_A / 999.84,
  molecularMass: 0.018015 / N_A,
  dielectric: 87.9,
  saturationPa: 611,
}
line(
  'threshold by the rate criterion used here (1 per ion per s)',
  thresholdSaturation((x) => ionNucleationRate(water, 273.15, x), 1).toFixed(2),
)
line("Thomson's barrierless limit", ionBarrierlessSaturation(water, 273.15).toFixed(2))
line("Wilson's measured value for negative ions (textbook)", '4.2')
line(
  'isopropanol at 243 K: rate criterion / barrierless / fog',
  `${ionThreshold(243).toFixed(2)} / ${ionBarrierlessSaturation(ipaLiquid(243), 243).toFixed(2)} / ${fogThreshold(243).toFixed(2)}`,
)

console.log('\nHow much the "one nucleation per ion per second" convention matters')
function plateLimit(onset: number): number {
  let cold = 233.15
  let warm = 258.15
  for (let i = 0; i < 24; i++) {
    const mid = (cold + warm) / 2
    const q = solveChamber({ ...DEFAULT_CHAMBER, plateK: mid })
    let sensitive = false
    for (let k = 1; k < q.z.length - 1 && !sensitive; k += 4) {
      const t = q.tempK[k]
      const liquid = ipaLiquid(t)
      sensitive =
        q.saturation[k] >= thresholdSaturation((x) => ionNucleationRate(liquid, t, x), onset)
    }
    if (sensitive) cold = mid
    else warm = mid
  }
  return cold - 273.15
}
for (const onset of [1, 100])
  line(
    `warmest plate with a sensitive layer, onset ${onset} per ion per s, C`,
    plateLimit(onset).toFixed(1),
  )

console.log('\nOther numbers quoted in the README')
{
  const row07 = row(astar, 0.7)
  line(
    'Bethe formula above NIST ASTAR at 0.7 MeV (alpha Bragg peak)',
    `${((betheStopping(2, MASS_MEV.alpha, 0.7) / row07[1] - 1) * 100).toFixed(0)}%`,
  )
  line(
    'width of a 3-pixel sprite in the 720-pixel top view, mm',
    ((3 / topViewScale(720, 0.1)) * 1000).toFixed(2),
  )
  const e = new Engine({ seed: 5, ...quiet })
  e.fire({ species: 'muon', kinetic: 3000, pos: [-0.0999, 0, 0.03], dir: [1, 0, 0] })
  const off: number[] = []
  for (let k = 0; k < e.count; k++) off.push(Math.hypot(e.pos[3 * k + 1], e.pos[3 * k + 2] - 0.03))
  off.sort((a, b) => a - b)
  line(
    'muon droplets: median distance from the track line, micrometres',
    (off[Math.floor(off.length / 2)] * 1e6).toFixed(0),
  )
  for (const energy of [5.4895, 7.6868]) {
    const a = new Engine({ seed: 5, ...quiet })
    a.fire({ species: 'alpha', kinetic: energy, pos: [-0.09, 0, 0.03], dir: [1, 0, 0] })
    let lo = Infinity
    let hi = -Infinity
    for (let k = 0; k < a.count; k++) {
      lo = Math.min(lo, a.pos[3 * k])
      hi = Math.max(hi, a.pos[3 * k])
    }
    line(`radon-chain alpha ${energy} MeV: track length at 3 cm, cm`, ((hi - lo) * 100).toFixed(1))
  }
}
