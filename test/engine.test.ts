import { expect, test } from 'bun:test'
import { airDensity, airViscosity } from '../src/air'
import { DEFAULT_CHAMBER, sampleProfile } from '../src/chamber'
import { ATM, G_ACCEL, MASS_MEV } from '../src/constants'
import { growthCoefficient } from '../src/droplet'
import { COLUMN_OFFSET, downdraftFactor, Engine, KIND, TRACK_RADIUS } from '../src/engine'
import { liquidDensity, vaporDiffusivity } from '../src/ipa'
import { Rng } from '../src/rng'
import { cosmicMuons } from '../src/sources'
import { betheStopping } from '../src/stopping'
import { loadTable, rel } from './reference'

const quiet = { background: false, mist: false } as const

function summary(e: Engine, kind: number) {
  let n = 0
  let r = 0
  let z = 0
  let xMin = Infinity
  let xMax = -Infinity
  for (let k = 0; k < e.count; k++) {
    if (e.kind[k] !== kind || e.radius[k] <= 0) continue
    n++
    r += e.radius[k]
    z += e.pos[3 * k + 2]
    xMin = Math.min(xMin, e.pos[3 * k])
    xMax = Math.max(xMax, e.pos[3 * k])
  }
  return { n, radius: r / n, z: z / n, extent: xMax - xMin }
}

test('the same seed gives the same chamber, droplet for droplet', () => {
  const run = () => {
    const e = new Engine({ seed: 42 })
    for (let i = 0; i < 90; i++) e.step(1 / 60)
    let hash = 0
    for (let k = 0; k < e.count; k++)
      hash = (hash * 31 + Math.round(e.pos[3 * k] * 1e7) + Math.round(e.radius[k] * 1e12)) | 0
    return { count: e.count, hash, events: e.events.length }
  }
  const a = run()
  expect(a.count).toBeGreaterThan(1000)
  expect(run()).toEqual(a)
  const other = new Engine({ seed: 43 })
  for (let i = 0; i < 90; i++) other.step(1 / 60)
  expect(other.count).not.toBe(a.count)
})

test('an alpha track in the chamber is as long as the NIST range in gas of the local density', async () => {
  // Fired level at 3 cm, where the gas is about -20 C and 16% denser than at room temperature.
  const astar = await loadTable('astar-air.txt')
  const range55 = astar.find((r) => r[0] === 5.5)
  if (!range55) throw new Error('no ASTAR row')
  // One track can be a few percent off (range straggling, a kink from a nuclear collision), so
  // the comparison uses the median of seven.
  const extents: number[] = []
  let rho = 0
  for (let seed = 1; seed <= 7; seed++) {
    const e = new Engine({ seed, ...quiet })
    e.fire({ species: 'alpha', kinetic: 5.5, pos: [-0.05, 0, 0.03], dir: [1, 0, 0] })
    for (let i = 0; i < 12; i++) e.step(1 / 60)
    extents.push(summary(e, KIND.alpha).extent)
    rho = airDensity(sampleProfile(e.profile, e.profile.tempK, 0.03), ATM)
  }
  extents.sort((x, y) => x - y)
  expect(rel(extents[3], range55[5] / rho / 100)).toBeLessThan(0.025)
  expect(extents[3]).toBeLessThan(0.0416) // shorter than in room-temperature air
})

test('a muon track is thin, sinks, and is gone in a few seconds', () => {
  // 20 cm of track at roughly 100-150 droplets per cm (two ions per pair, cold dense gas). A hard
  // knock-on electron can add thousands of droplets to one track, so the count is judged by its
  // median over several seeds, not by a single run.
  const counts: number[] = []
  for (let seed = 1; seed <= 9; seed++) {
    const probe = new Engine({ seed, ...quiet })
    counts.push(
      probe.fire({ species: 'muon', kinetic: 3000, pos: [-0.0999, 0, 0.03], dir: [1, 0, 0] })
        .droplets,
    )
  }
  counts.sort((x, y) => x - y)
  expect(counts[4]).toBeGreaterThan(1800)
  expect(counts[4]).toBeLessThan(3500)
  const e = new Engine({ seed: 6, ...quiet })
  const record = e.fire({ species: 'muon', kinetic: 3000, pos: [-0.0999, 0, 0.03], dir: [1, 0, 0] })
  for (let i = 0; i < 30; i++) e.step(1 / 60)
  const early = summary(e, KIND.muon)
  expect(early.z).toBeLessThan(0.029)
  expect(early.z).toBeGreaterThan(0.015)
  expect(early.radius).toBeGreaterThan(5e-6)
  expect(early.radius).toBeLessThan(20e-6)
  for (let i = 0; i < 150; i++) e.step(1 / 60)
  expect(summary(e, KIND.muon).n).toBe(0) // all landed within 3 s
  expect(e.landed).toBeGreaterThan(0.9 * record.droplets)
  expect(e.dropped).toBe(0)
})

test('track downdraft: the engine matches a direct sum of the flow each droplet induces (3%), which also fixes the offset', () => {
  // Two droplets with the same history; only the first is told it belongs to a 4 cm track. The
  // difference in how far they fall in one step is the engine's downdraft. Compare it with the
  // mean fall speed induced on each droplet by all the others (Oseen tensor), for droplets
  // scattered in a thin horizontal cylinder. Both sides are slow viscous flow in unbounded gas;
  // the plate below and inertia are not in this check.
  const perMetre = 1e4
  const z0 = 0.03
  for (const length of [0.04, 0.2]) {
    const e = new Engine({ seed: 1, ...quiet })
    e.add(-0.01, 0, z0, 1, perMetre, 0, KIND.alpha)
    e.add(0.01, 0, z0, 1, perMetre, 0, KIND.alpha)
    e.radius[0] = 1e-6
    e.radius[1] = 1e-6
    e.settle[0] = Math.log(length / TRACK_RADIUS) + COLUMN_OFFSET
    // One step long enough that the difference in height is far above the 2 nm resolution of
    // the single-precision positions. The engine uses the radius at the end of the step for
    // the whole step, and so does the comparison below.
    const dt = 0.01
    e.step(dt)
    const draft = (e.pos[5] - e.pos[2]) / dt
    expect(draft).toBeGreaterThan(0)

    const count = 2500
    const rng = new Rng(3)
    const pts = new Float64Array(3 * count)
    for (let i = 0; i < count; i++) {
      const rr = TRACK_RADIUS * Math.sqrt(rng.next())
      const phi = 2 * Math.PI * rng.next()
      pts[3 * i] = (rng.next() - 0.5) * length
      pts[3 * i + 1] = rr * Math.cos(phi)
      pts[3 * i + 2] = rr * Math.sin(phi)
    }
    let induced = 0 // mean fall speed in units of F / (8 pi eta), F the weight of one droplet
    for (let i = 0; i < count; i++)
      for (let j = 0; j < count; j++) {
        if (i === j) continue
        const dz = pts[3 * i + 2] - pts[3 * j + 2]
        const d = Math.hypot(pts[3 * i] - pts[3 * j], pts[3 * i + 1] - pts[3 * j + 1], dz)
        induced += (1 + (dz * dz) / (d * d)) / d / count
      }
    // The offset the sum implies: induced = 2 (N / L) (ln(L / R) + offset).
    const offset = induced / ((2 * count) / length) - Math.log(length / TRACK_RADIUS)
    expect(Math.abs(offset - COLUMN_OFFSET)).toBeLessThan(0.1)
    const temp = sampleProfile(e.profile, e.profile.tempK, z0)
    const weight = (4 / 3) * Math.PI * e.radius[0] ** 3 * liquidDensity(temp) * G_ACCEL
    // Scale the sum from this many droplets to the engine's line density.
    const direct =
      (weight / (8 * Math.PI * airViscosity(temp))) * induced * ((perMetre * length) / count)
    expect(rel(draft, direct)).toBeLessThan(0.03)
  }
})

test("track downdraft: inertia caps it, and the capped factor is the root of Lamb's equation at every speed", () => {
  const nu = 1.1e-5
  const slender = Math.log(0.2 / TRACK_RADIUS) + COLUMN_OFFSET
  // Slow: the viscous law applies unchanged.
  expect(downdraftFactor(slender, 1e-7, nu)).toBe(slender)
  expect(downdraftFactor(0, 1e-3, nu)).toBe(0)
  // Faster: K solves K + ln K = 1/2 - gamma - ln(unit R / 4 nu). Solve that here by bisection,
  // independently of the engine's Newton iteration, from light columns to very heavy ones.
  let previousSpeed = 0
  for (const unit of [1e-5, 1e-4, 1e-3, 3e-3, 0.01, 0.035, 0.065, 0.2, 1]) {
    const target = 0.5 - 0.5772156649 - Math.log((unit * TRACK_RADIUS) / (4 * nu))
    let lo = 1e-9
    let hi = 50
    for (let i = 0; i < 100; i++) {
      const mid = 0.5 * (lo + hi)
      if (mid + Math.log(mid) < target) lo = mid
      else hi = mid
    }
    const k = downdraftFactor(slender, unit, nu)
    expect(rel(k, Math.min(lo, slender))).toBeLessThan(1e-6)
    expect(k).toBeGreaterThan(0)
    // A heavier column never sinks more slowly.
    expect(unit * k).toBeGreaterThan(previousSpeed)
    previousSpeed = unit * k
  }
  // Inputs seen on a default-chamber alpha track, where the first version of the solver returned 0.
  expect(downdraftFactor(6.7152352333, 0.1010818365, 0.000011560139)).toBeCloseTo(1.23347, 4)
})

test('the downdraft is what clears an alpha track quickly, and it speeds a muon track up too', () => {
  // Seconds until 90% of the droplets (real droplets, by weight) have landed.
  const landing = (species: 'alpha' | 'muon', collectiveSettling: boolean) => {
    const e = new Engine({ seed: 6, ...quiet, collectiveSettling })
    const start = species === 'alpha' ? -0.05 : -0.0999
    e.fire({
      species,
      kinetic: species === 'alpha' ? 5.5 : 3000,
      pos: [start, 0, 0.03],
      dir: [1, 0, 0],
    })
    for (let i = 1; i <= 3600; i++) {
      e.step(1 / 60)
      if (e.landed >= 0.9 * e.created) return i / 60
    }
    return Infinity
  }
  const alpha = landing('alpha', true)
  expect(alpha).toBeGreaterThan(0.3)
  expect(alpha).toBeLessThan(5)
  expect(landing('alpha', false)).toBeGreaterThan(5 * alpha)
  expect(landing('muon', false)).toBeGreaterThan(landing('muon', true))
})

test('a lone droplet grows as r^2 = 2 g (S - 1) t', () => {
  const e = new Engine({ seed: 7, ...quiet })
  const z0 = 0.03
  e.add(0, 0, z0, 1, 0, 0, KIND.mist)
  for (let i = 0; i < 100; i++) e.step(1e-3)
  expect(e.count).toBe(1)
  const t = sampleProfile(e.profile, e.profile.tempK, z0)
  const s = sampleProfile(e.profile, e.profile.saturation, z0)
  // It falls a fraction of a millimetre in this time, so the conditions are those at z0.
  expect(rel(e.radius[0], Math.sqrt(2 * growthCoefficient(t, ATM) * (s - 1) * 0.1))).toBeLessThan(
    0.01,
  )
  expect(z0 - e.pos[2]).toBeGreaterThan(0)
  expect(z0 - e.pos[2]).toBeLessThan(5e-4)
})

test('shared vapor: droplets in the engine grow as a direct solve of radial diffusion says (4%)', () => {
  // The engine divides a droplet's growth by 1 + 2 n r L. Here a droplet is put in the engine
  // with n neighbours per metre and its radius is compared with an explicit finite-difference
  // solution of the supersaturation field around a line of such droplets. Both sides take the
  // droplets to fill a core of the same radius and heat to spread like vapor; what is checked is
  // the quasi-steady logarithm, for a dense alpha-like column and a sparse muon-like one.
  const z0 = 0.03
  for (const [perMetre, seconds] of [
    [8e6, 0.5],
    [1.6e4, 0.15],
  ]) {
    const e = new Engine({ seed: 1, ...quiet })
    e.add(0, 0, z0, 1, perMetre, 0, KIND.alpha)
    const steps = Math.round(seconds / 1e-3)
    for (let i = 0; i < steps; i++) e.step(1e-3)
    expect(e.count).toBe(1)

    const temp = sampleProfile(e.profile, e.profile.tempK, z0)
    const excess = sampleProfile(e.profile, e.profile.saturation, z0) - 1
    const g = growthCoefficient(temp, ATM)
    const diff = vaporDiffusivity(temp, ATM)
    const cell = 2e-5
    const cells = 1200
    const core = Math.round(TRACK_RADIUS / cell)
    const perVolume = perMetre / (Math.PI * TRACK_RADIUS * TRACK_RADIUS)
    const s = new Float64Array(cells).fill(excess)
    const next = new Float64Array(cells)
    const h = (0.2 * cell * cell) / diff
    let rd = 0
    for (let t = 0; t < seconds; t += h) {
      let inCore = 0
      for (let i = 0; i < core; i++) inCore += s[i] * (2 * i + 1)
      inCore /= core * core
      for (let i = 0; i < cells - 1; i++) {
        const up = (s[i + 1] - s[i]) / cell
        const down = i > 0 ? (s[i] - s[i - 1]) / cell : 0
        const sink = i < core ? 4 * Math.PI * diff * perVolume * rd * s[i] : 0
        next[i] = s[i] + h * ((diff * ((i + 1) * up - i * down)) / ((i + 0.5) * cell) - sink)
      }
      next[cells - 1] = excess
      s.set(next)
      rd = Math.sqrt(Math.max(0, rd * rd + 2 * g * inCore * h))
    }
    expect(rel(e.radius[0], rd)).toBeLessThan(0.04)
    if (perMetre > 1e6) expect(e.radius[0]).toBeLessThan(0.25 * Math.sqrt(2 * g * excess * seconds)) // starved
  }
})

test('alpha droplets are starved of vapor: smaller than muon droplets of the same age', () => {
  const a = new Engine({ seed: 8, ...quiet })
  a.fire({ species: 'alpha', kinetic: 5.5, pos: [-0.05, 0, 0.03], dir: [1, 0, 0] })
  const m = new Engine({ seed: 8, ...quiet })
  m.fire({ species: 'muon', kinetic: 3000, pos: [-0.0999, 0, 0.03], dir: [1, 0, 0] })
  for (let i = 0; i < 18; i++) {
    a.step(1 / 60)
    m.step(1 / 60)
  }
  const alpha = summary(a, KIND.alpha)
  const muon = summary(m, KIND.muon)
  expect(alpha.radius).toBeLessThan(0.6 * muon.radius)
  // But there are so many of them that the track scatters more light per unit length.
  let alphaLight = 0
  for (let k = 0; k < a.count; k++) alphaLight += a.weight[k] * a.radius[k] ** 2
  let muonLight = 0
  for (let k = 0; k < m.count; k++) muonLight += m.weight[k] * m.radius[k] ** 2
  expect(alphaLight / alpha.extent).toBeGreaterThan(3 * (muonLight / muon.extent))
})

test('the edge of the sensitive layer is soft: tracks fade out as the plate warms, they do not switch off', () => {
  // The layer is defined by one nucleation per ion per second. Just outside it an ion still has
  // a small chance, so an alpha over a plate slightly too warm leaves a sparse, faint track.
  const fraction = (plateC: number) => {
    const e = new Engine({
      seed: 5,
      ...quiet,
      chamber: { ...DEFAULT_CHAMBER, plateK: plateC + 273.15 },
    })
    e.fire({ species: 'alpha', kinetic: 5.5, pos: [-0.05, 0, 0.03], dir: [1, 0, 0] })
    let real = 0
    for (let k = 0; k < e.count; k++) real += e.weight[k]
    return real / (2 * e.ionPairs)
  }
  expect(fraction(-40)).toBeGreaterThan(0.95)
  const at20 = fraction(-20)
  expect(at20).toBeGreaterThan(0)
  expect(at20).toBeLessThan(0.05)
  expect(fraction(-15)).toBeLessThan(0.001)
  const cold = new Engine({ seed: 5, ...quiet })
  expect(cold.nucleationChance(0.03)).toBeGreaterThan(0.999)
  expect(cold.nucleationChance(0.1)).toBeLessThan(1e-6)
})

test('nucleation delays are spread over the ion lifetime, not piled up at its end', () => {
  const e = new Engine({ seed: 5, ...quiet, chamber: { ...DEFAULT_CHAMBER, plateK: 252.15 } })
  e.fire({ species: 'alpha', kinetic: 5.5, pos: [-0.05, 0, 0.03], dir: [1, 0, 0] })
  expect(e.count).toBeGreaterThan(200)
  let late = 0
  let early = 0
  for (let k = 0; k < e.count; k++) {
    expect(e.birth[k]).toBeLessThanOrEqual(1)
    if (e.birth[k] > 0.99) late++
    if (e.birth[k] < 0.5) early++
  }
  expect(late / e.count).toBeLessThan(0.03)
  expect(early / e.count).toBeGreaterThan(0.4)
  expect(early / e.count).toBeLessThan(0.6)
})

test('droplets are counted as droplets and conserved: all that are created land', () => {
  const e = new Engine({ seed: 4, ...quiet })
  e.fire({ species: 'alpha', kinetic: 5.5, pos: [-0.05, 0, 0.03], dir: [1, 0, 0] })
  let weights = 0
  for (let k = 0; k < e.count; k++) weights += e.weight[k]
  // `created` counts real droplets: about two per ion pair here, many per array entry.
  expect(e.created).toBe(weights)
  expect(rel(e.created, 2 * e.ionPairs)).toBeLessThan(0.02)
  expect(e.created).toBeGreaterThan(5 * e.count)
  for (let i = 0; i < 3600 && e.count > 0; i++) e.step(1 / 60)
  expect(e.count).toBe(0)
  expect(e.landed).toBe(e.created)
  expect(e.dropped).toBe(0)
})

test('no droplets form outside the sensitive layer', () => {
  const e = new Engine({ seed: 9, ...quiet })
  const high = e.profile.sensitiveTop + 0.02
  const record = e.fire({ species: 'muon', kinetic: 3000, pos: [-0.0999, 0, high], dir: [1, 0, 0] })
  expect(record.droplets).toBe(0)
  expect(e.ionPairs).toBeGreaterThan(500) // the ions are there; the vapor is not supersaturated enough
})

test('cosmic muons ionize the chamber gas at the rate their flux and energy loss imply', () => {
  // Omnidirectional muon fluence rate is (2 pi / 3) I_v. Times ion pairs per metre of track, that
  // is ion pairs per m^3 per s. Fast knock-on electrons leave a chamber this small, so the
  // simulated rate must sit between the estimates from the restricted and the full energy loss.
  const e = new Engine({ seed: 2, ...quiet, capacity: 16 })
  e.sources.push(cosmicMuons(e.box))
  const seconds = 60
  for (let i = 0; i < seconds * 10; i++) e.step(0.1)
  const volume = 4 * e.box.halfX * e.box.halfY * e.box.height
  let rho = 0
  for (let i = 0; i < e.profile.z.length; i++)
    rho += airDensity(e.profile.tempK[i], ATM) / e.profile.z.length
  const fluence = ((2 * Math.PI) / 3) * 70
  const pairsPerMetre = (cut: number) =>
    (betheStopping(1, MASS_MEV.muon, 3000, cut) * rho * 100) / 33.97e-6
  const rate = e.ionPairs / volume / seconds
  expect(rate).toBeGreaterThan(0.95 * fluence * pairsPerMetre(0.01))
  expect(rate).toBeLessThan(1.05 * fluence * pairsPerMetre(Infinity))
  // About one ion pair per cm^3 per second, the order of the measured sea-level cosmic ionization.
  expect(rate * 1e-6).toBeGreaterThan(0.5)
  expect(rate * 1e-6).toBeLessThan(2)
})

test('the mist settles to a steady population and the arrays never overflow', () => {
  const e = new Engine({ seed: 10, background: false })
  for (let i = 0; i < 480; i++) e.step(1 / 60)
  const at8 = e.census().mist
  for (let i = 0; i < 240; i++) e.step(1 / 60)
  const at12 = e.census().mist
  expect(at8).toBeGreaterThan(5000)
  expect(rel(at12, at8)).toBeLessThan(0.1)
  expect(e.dropped).toBe(0)
})
