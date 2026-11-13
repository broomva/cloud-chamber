/**
 * Browser side of the render gate. Runs a fixed scenario, draws it, reads the pixels back and
 * measures the picture. Bundled by tools/render.ts and driven from there.
 */
import { airDensity } from '../src/air'
import { DEFAULT_CHAMBER, sampleProfile } from '../src/chamber'
import { Engine } from '../src/engine'
import { type Camera, DEFAULT_DRAW, Renderer, topViewScale } from '../src/render'
import { americiumSource, strontiumSource } from '../src/sources'

export interface Measurement {
  /** Pixels brighter than the threshold. */
  lit: number
  /** Brightest red, green or blue value anywhere in the frame, 0..255. 0 means exactly black. */
  brightest: number
  /** Extent of the lit pixels along the track's principal axis, with the sprite width removed, m. */
  lengthM: number
  /** RMS distance of lit pixels from the best straight line through them, m. */
  residualM: number
  /** Mean and standard deviation of luma over the frame, 0..1. */
  mean: number
  std: number
  /** Gas density where the test particle was fired, g/cm^3. */
  gasDensity: number
  /** Super-droplets alive when the frame was drawn. */
  droplets: number
  /** Fraction of pixels that changed between this frame and one a second later (scene only). */
  changed: number
}

declare global {
  interface Window {
    runScenario: (name: string, size: number) => Measurement
  }
}

const TOP: Camera = { mode: 'top', yaw: 0, pitch: 0, distance: 1, target: [0, 0, 0] }
const ORBIT: Camera = {
  mode: 'orbit',
  yaw: -1.9,
  pitch: 0.42,
  distance: 0.42,
  target: [0, 0, 0.025],
}
const FIRE_HEIGHT = 0.03
const THRESHOLD = 24

function measure(renderer: Renderer, pixelsPerMetre: number, sprite: number) {
  const { width, height, data } = renderer.readPixels()
  const xs: number[] = []
  const ys: number[] = []
  let sum = 0
  let sum2 = 0
  let brightest = 0
  for (let y = 0; y < height; y++)
    for (let x = 0; x < width; x++) {
      const i = 4 * (y * width + x)
      const luma = 0.2126 * data[i] + 0.7152 * data[i + 1] + 0.0722 * data[i + 2]
      brightest = Math.max(brightest, data[i], data[i + 1], data[i + 2])
      sum += luma / 255
      sum2 += (luma / 255) ** 2
      if (luma > THRESHOLD) {
        xs.push(x)
        ys.push(y)
      }
    }
  const count = width * height
  const mean = sum / count
  const out = {
    lit: xs.length,
    brightest,
    lengthM: 0,
    residualM: 0,
    mean,
    std: Math.sqrt(Math.max(0, sum2 / count - mean * mean)),
  }
  if (xs.length < 2) return out
  // Principal axis of the lit pixels.
  const mx = xs.reduce((s, v) => s + v, 0) / xs.length
  const my = ys.reduce((s, v) => s + v, 0) / ys.length
  let sxx = 0
  let sxy = 0
  let syy = 0
  for (let k = 0; k < xs.length; k++) {
    sxx += (xs[k] - mx) ** 2
    sxy += (xs[k] - mx) * (ys[k] - my)
    syy += (ys[k] - my) ** 2
  }
  const angle = 0.5 * Math.atan2(2 * sxy, sxx - syy)
  const ax = Math.cos(angle)
  const ay = Math.sin(angle)
  let lo = Infinity
  let hi = -Infinity
  let perp = 0
  for (let k = 0; k < xs.length; k++) {
    const along = (xs[k] - mx) * ax + (ys[k] - my) * ay
    const across = -(xs[k] - mx) * ay + (ys[k] - my) * ax
    lo = Math.min(lo, along)
    hi = Math.max(hi, along)
    perp += across * across
  }
  // A sprite lights about half its diameter beyond the droplet at each end of the track.
  out.lengthM = Math.max(0, hi - lo + 1 - sprite / 2) / pixelsPerMetre
  out.residualM = Math.sqrt(perp / xs.length) / pixelsPerMetre
  return out
}

function canvasOf(size: number): HTMLCanvasElement {
  const canvas = document.getElementById('c') as HTMLCanvasElement
  canvas.width = size
  canvas.height = size
  return canvas
}

let renderer: Renderer | undefined

window.runScenario = (name, size) => {
  const canvas = canvasOf(size)
  renderer ??= new Renderer(canvas)
  const quiet = { background: false, mist: false }
  const measuring = { ...DEFAULT_DRAW, scenery: false, exposure: 6, sprite: 3 }
  const frame = 1 / 60
  let engine: Engine
  let camera = TOP
  let options = measuring
  let changed = 0

  if (name === 'empty') {
    engine = new Engine({ seed: 1, ...quiet })
    for (let i = 0; i < 30; i++) engine.step(frame)
  } else if (name === 'alpha' || name === 'alpha-warm-plate') {
    // -15 C is well past the soft edge of the sensitive layer (it fades out near -21 C). The
    // warm frame is read after the ion lifetime, when every droplet that will form has formed.
    const chamber = name === 'alpha' ? DEFAULT_CHAMBER : { ...DEFAULT_CHAMBER, plateK: 258.15 }
    engine = new Engine({ seed: 5, ...quiet, chamber })
    engine.fire({
      species: 'alpha',
      kinetic: 5.5,
      pos: [-0.05, 0.01, FIRE_HEIGHT],
      dir: [Math.cos(0.3), Math.sin(0.3), 0],
    })
    const frames = name === 'alpha' ? 12 : 66
    for (let i = 0; i < frames; i++) engine.step(frame)
  } else if (name === 'muon') {
    engine = new Engine({ seed: 6, ...quiet })
    const norm = Math.hypot(1, 0.3)
    engine.fire({
      species: 'muon',
      kinetic: 3000,
      pos: [-0.0999, -0.03, FIRE_HEIGHT],
      dir: [1 / norm, 0.3 / norm, 0],
    })
    for (let i = 0; i < 30; i++) engine.step(frame)
  } else if (name === 'electron') {
    engine = new Engine({ seed: 7, ...quiet })
    engine.fire({ species: 'electron', kinetic: 0.03, pos: [0, 0, FIRE_HEIGHT], dir: [1, 0, 0] })
    for (let i = 0; i < 30; i++) engine.step(frame)
  } else if (name === 'scene') {
    engine = new Engine({ seed: 11 })
    engine.sources.push(
      americiumSource([-0.098, 0.03, 0.03], [1, -0.2, 0], 2, 0.35),
      strontiumSource([0.098, -0.04, 0.03], 3),
    )
    for (let i = 0; i < 360; i++) engine.step(frame)
    camera = ORBIT
    options = {
      ...DEFAULT_DRAW,
      exposure: Number(new URLSearchParams(location.search).get('exposure') ?? 0.12),
      sprite: Number(new URLSearchParams(location.search).get('sprite') ?? 3),
    }
    renderer.draw(engine, camera, options)
    const before = renderer.readPixels().data
    for (let i = 0; i < 60; i++) engine.step(frame)
    renderer.draw(engine, camera, options)
    const after = renderer.readPixels().data
    let differ = 0
    for (let i = 0; i < before.length; i += 4)
      if (Math.abs(before[i + 1] - after[i + 1]) > 8) differ++
    changed = differ / (before.length / 4)
  } else {
    throw new Error(`unknown scenario: ${name}`)
  }

  renderer.draw(engine, camera, options)
  const m = measure(renderer, topViewScale(size, engine.box.halfX), options.sprite)
  return {
    ...m,
    gasDensity: airDensity(
      sampleProfile(engine.profile, engine.profile.tempK, FIRE_HEIGHT),
      engine.profile.config.pressurePa,
    ),
    droplets: engine.count,
    changed,
  }
}
