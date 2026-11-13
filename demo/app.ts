/**
 * The interactive page: a running chamber, a few controls, and a plot of the solved profiles so
 * the sensitive layer can be seen moving when the plate temperature changes.
 */
import { DEFAULT_CHAMBER } from '../src/chamber'
import { Engine } from '../src/engine'
import { type Camera, DEFAULT_DRAW, Renderer } from '../src/render'
import { americiumSource, strontiumSource } from '../src/sources'
import type { Particle } from '../src/transport'

const $ = <T extends HTMLElement>(id: string) => document.getElementById(id) as T
const canvas = $<HTMLCanvasElement>('view')
const plot = $<HTMLCanvasElement>('plot')
const renderer = new Renderer(canvas)
const camera: Camera = {
  mode: 'orbit',
  yaw: -1.9,
  pitch: 0.42,
  distance: 0.42,
  target: [0, 0, 0.025],
}

const state = {
  plateC: -40,
  topC: 20,
  background: true,
  mist: true,
  americium: false,
  strontium: false,
  tint: false,
  paused: false,
}
let engine: Engine
let seed = 1

function rebuild() {
  engine = new Engine({
    seed: seed++,
    chamber: { ...DEFAULT_CHAMBER, plateK: state.plateC + 273.15, topK: state.topC + 273.15 },
    background: state.background,
    mist: state.mist,
  })
  const z = fireHeight()
  if (state.americium)
    engine.sources.push(americiumSource([-0.098, 0.03, z], [1, -0.2, 0], 0.5, 0.35))
  if (state.strontium) engine.sources.push(strontiumSource([0.098, -0.04, z], 2))
  drawPlot()
  const p = engine.profile
  const depth = p.sensitiveTop - p.sensitiveBottom
  let best = 0
  for (let i = 0; i < p.z.length; i++) best = Math.max(best, engine.nucleationChance(p.z[i]))
  $('layer').textContent =
    depth > 0
      ? `Sensitive layer: ${(p.sensitiveBottom * 1000).toFixed(0)} to ${(p.sensitiveTop * 1000).toFixed(0)} mm above the plate`
      : best > 0.001
        ? `No sensitive layer. At the best height ${(best * 100).toFixed(best < 0.1 ? 1 : 0)}% of ions still seed a droplet, so tracks are faint`
        : 'No sensitive layer: the plate is not cold enough for ions to seed droplets'
}

/** Middle of the sensitive layer, or 2 cm up if there is none. */
function fireHeight(): number {
  const p = engine.profile
  return p.sensitiveTop > 0 ? 0.5 * (p.sensitiveBottom + Math.min(p.sensitiveTop, 0.04)) : 0.02
}

function fire(species: Particle['species']) {
  const z = fireHeight()
  const angle = (Math.random() - 0.5) * 0.8
  const dir: [number, number, number] = [Math.cos(angle), Math.sin(angle), 0]
  if (species === 'alpha')
    engine.fire({ species, kinetic: 5.4856, pos: [-0.06, (Math.random() - 0.5) * 0.1, z], dir })
  else if (species === 'electron')
    engine.fire({
      species,
      kinetic: 0.05 + Math.random() * 0.5,
      pos: [-0.0999, (Math.random() - 0.5) * 0.1, z],
      dir,
    })
  else engine.fire({ species, kinetic: 3000, pos: [-0.0999, (Math.random() - 0.5) * 0.1, z], dir })
}

function drawPlot() {
  const ctx = plot.getContext('2d')
  if (!ctx) return
  const dpr = Math.min(window.devicePixelRatio || 1, 2)
  const w = plot.clientWidth
  const h = plot.clientHeight
  plot.width = w * dpr
  plot.height = h * dpr
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0)
  const css = getComputedStyle(document.documentElement)
  const fg = css.getPropertyValue('--foreground').trim() || '#e8eaed'
  const muted = css.getPropertyValue('--muted-foreground').trim() || '#8b929c'
  const p = engine.profile
  const sMax = 6
  const left = 34
  const bottom = h - 22
  const x = (s: number) => left + ((w - left - 8) * (s - 1)) / (sMax - 1)
  const y = (z: number) => bottom - ((bottom - 8) * z) / p.config.height
  ctx.clearRect(0, 0, w, h)
  if (p.sensitiveTop > 0) {
    ctx.fillStyle = 'rgba(120, 190, 255, 0.16)'
    ctx.fillRect(left, y(p.sensitiveTop), w - left - 8, y(p.sensitiveBottom) - y(p.sensitiveTop))
  }
  ctx.strokeStyle = muted
  ctx.lineWidth = 1
  ctx.beginPath()
  ctx.moveTo(left, 8)
  ctx.lineTo(left, bottom)
  ctx.lineTo(w - 8, bottom)
  ctx.stroke()
  ctx.fillStyle = muted
  ctx.font = '10px system-ui, sans-serif'
  ctx.textAlign = 'center'
  for (let s = 1; s <= sMax; s++) ctx.fillText(String(s), x(s), bottom + 12)
  ctx.textAlign = 'right'
  for (const cm of [0, 4, 8, 12])
    if (cm / 100 <= p.config.height) ctx.fillText(`${cm}`, left - 4, y(cm / 100) + 3)
  const curve = (values: Float64Array, colour: string, dash: number[]) => {
    ctx.strokeStyle = colour
    ctx.setLineDash(dash)
    ctx.lineWidth = dash.length ? 1 : 1.8
    ctx.beginPath()
    for (let i = 0; i < p.z.length; i++) {
      const px = x(Math.min(values[i], sMax))
      if (i === 0) ctx.moveTo(px, y(p.z[i]))
      else ctx.lineTo(px, y(p.z[i]))
    }
    ctx.stroke()
    ctx.setLineDash([])
  }
  curve(p.fogThreshold, '#e0a050', [4, 3])
  curve(p.ionThreshold, '#78beff', [4, 3])
  curve(p.saturation, fg, [])
}

// Controls.
function bindRange(id: string, key: 'plateC' | 'topC', label: string) {
  const input = $<HTMLInputElement>(id)
  const show = () => {
    $(`${id}-value`).textContent = `${label} ${input.value} °C`
  }
  input.addEventListener('input', () => {
    state[key] = Number(input.value)
    show()
    rebuild()
  })
  show()
}
function bindCheck(
  id: string,
  key: 'background' | 'mist' | 'americium' | 'strontium' | 'tint',
  needsRebuild = true,
) {
  const input = $<HTMLInputElement>(id)
  input.addEventListener('change', () => {
    state[key] = input.checked
    if (needsRebuild) rebuild()
  })
}
bindRange('plate', 'plateC', 'Plate')
bindRange('top', 'topC', 'Top')
bindCheck('background', 'background')
bindCheck('mist', 'mist')
bindCheck('americium', 'americium')
bindCheck('strontium', 'strontium')
bindCheck('tint', 'tint', false)
$('fire-alpha').addEventListener('click', () => fire('alpha'))
$('fire-electron').addEventListener('click', () => fire('electron'))
$('fire-muon').addEventListener('click', () => fire('muon'))
$('pause').addEventListener('click', () => {
  state.paused = !state.paused
  $('pause').textContent = state.paused ? 'Resume' : 'Pause'
})

// Orbit with the pointer, zoom with the wheel.
let dragging = false
let lastX = 0
let lastY = 0
canvas.addEventListener('pointerdown', (e) => {
  dragging = true
  lastX = e.clientX
  lastY = e.clientY
  canvas.setPointerCapture(e.pointerId)
})
canvas.addEventListener('pointerup', () => {
  dragging = false
})
canvas.addEventListener('pointermove', (e) => {
  if (!dragging) return
  camera.yaw -= (e.clientX - lastX) * 0.006
  camera.pitch = Math.min(1.5, Math.max(0.03, camera.pitch + (e.clientY - lastY) * 0.006))
  lastX = e.clientX
  lastY = e.clientY
})
canvas.addEventListener(
  'wheel',
  (e) => {
    e.preventDefault()
    camera.distance = Math.min(1.2, Math.max(0.12, camera.distance * Math.exp(e.deltaY * 0.001)))
  },
  { passive: false },
)

const NAMES: Record<string, string> = { alpha: 'alpha', electron: 'electron', muon: 'muon' }
let last = performance.now()
let shown = 0
function frame(now: number) {
  const dpr = Math.min(window.devicePixelRatio || 1, 2)
  const w = Math.round(canvas.clientWidth * dpr)
  const h = Math.round(canvas.clientHeight * dpr)
  if (canvas.width !== w || canvas.height !== h) {
    canvas.width = w
    canvas.height = h
  }
  const elapsed = Math.min((now - last) / 1000, 0.05)
  last = now
  if (!state.paused) {
    const steps = Math.max(1, Math.ceil(elapsed / (1 / 60)))
    for (let i = 0; i < steps; i++) engine.step(elapsed / steps)
  }
  renderer.draw(engine, camera, {
    ...DEFAULT_DRAW,
    exposure: 0.09,
    sprite: 3 * dpr,
    colourByParticle: state.tint,
  })
  if (now - shown > 250) {
    shown = now
    const recent = engine.events
      .filter((ev) => ev.droplets > 0)
      .slice(-4)
      .reverse()
    $('events').textContent = recent.length
      ? recent
          .map(
            (ev) =>
              `${NAMES[ev.species]} ${ev.kineticMeV >= 100 ? `${(ev.kineticMeV / 1000).toFixed(1)} GeV` : `${ev.kineticMeV.toFixed(2)} MeV`}`,
          )
          .join(' · ')
      : 'no tracks yet'
    $('count').textContent = `${engine.count.toLocaleString()} droplets in flight`
  }
  requestAnimationFrame(frame)
}

rebuild()
requestAnimationFrame(frame)
