/**
 * Headless render gate. Runs fixed scenarios in real WebGL2 (Chromium via playwright-core,
 * SwiftShader by default), measures the rendered frames and exits non-zero unless the pictures
 * say what the physics says:
 *
 *   empty             no radiation, no mist                 -> a black frame
 *   alpha             one 5.5 MeV alpha, seen from above    -> a track as long as the NIST ASTAR
 *                                                              range in gas of the local density
 *   alpha-warm-plate  the same alpha over a -15 C plate,    -> a black frame (no sensitive layer)
 *                     read 1.1 s later
 *   muon              one muon across the chamber           -> a straight line, wall to wall
 *   electron          one 30 keV electron                   -> a short track, shorter than its
 *                                                              NIST ESTAR path length
 *   scene             everything on, perspective view       -> a frame with structure that moves
 *
 * Any console error or warning also fails the gate. PNGs and report.json land in --out
 * (default renders/, gitignored).
 *
 *   bun run tools/render.ts [--size 720] [--out renders]
 *   CHAMBER_GPU=1 uses the machine's GPU instead of SwiftShader.
 */
import { mkdir, writeFile } from 'node:fs/promises'
import { join, resolve } from 'node:path'
import { chromium } from 'playwright-core'
import type { Measurement } from './render-entry'

function arg(name: string, fallback: string): string {
  const i = process.argv.indexOf(`--${name}`)
  return (i >= 0 ? process.argv[i + 1] : undefined) ?? fallback
}
const size = Number(arg('size', '720'))
const outDir = resolve(import.meta.dir, '..', arg('out', 'renders'))

async function reference(file: string, energy: number, column: number): Promise<number> {
  const text = await Bun.file(join(import.meta.dir, '..', 'test', 'reference', file)).text()
  for (const line of text.split('\n')) {
    const cells = line.trim().split(/\s+/).map(Number)
    if (!line.startsWith('#') && cells[0] === energy) return cells[column]
  }
  throw new Error(`${file}: no row for ${energy} MeV`)
}
const alphaRange = await reference('astar-air.txt', 5.5, 5) // projected range, g/cm^2
const electronRange = await reference('estar-air.txt', 0.03, 4) // CSDA range, g/cm^2

const build = await Bun.build({
  entrypoints: [join(import.meta.dir, 'render-entry.ts')],
  target: 'browser',
  format: 'esm',
})
if (!build.success) {
  for (const log of build.logs) console.error(log)
  process.exit(1)
}
const [artifact] = build.outputs
if (!artifact) throw new Error('bundle produced no output')
const page = `<!doctype html><html><head><meta charset="utf-8"><style>html,body{margin:0;background:#000}canvas{display:block}</style></head>
<body><canvas id="c"></canvas><script type="module">${await artifact.text()}</script></body></html>`

// Loopback only.
const server = Bun.serve({
  hostname: '127.0.0.1',
  port: 0,
  fetch: () => new Response(page, { headers: { 'content-type': 'text/html' } }),
})
const gpu = process.env.CHAMBER_GPU === '1'
const browser = await chromium.launch({
  headless: true,
  args: gpu
    ? ['--ignore-gpu-blocklist', '--enable-gpu']
    : ['--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--ignore-gpu-blocklist'],
})

type Check = (m: Measurement) => string[]
const within = (label: string, value: number, expected: number, tolerance: number) =>
  Math.abs(value - expected) <= tolerance * expected
    ? []
    : [
        `${label}: ${value.toPrecision(4)} vs expected ${expected.toPrecision(4)} (tolerance ${(tolerance * 100).toFixed(0)}%)`,
      ]

const scenarios: { name: string; check: Check }[] = [
  {
    name: 'empty',
    check: (m) =>
      m.brightest === 0
        ? []
        : [`an empty chamber is not black: brightest channel value ${m.brightest}`],
  },
  {
    name: 'alpha',
    check: (m) => [
      ...within('alpha track length (m)', m.lengthM, alphaRange / m.gasDensity / 100, 0.04),
      ...(m.residualM < 6e-4
        ? []
        : [`alpha track is not straight: rms ${(m.residualM * 1e3).toFixed(2)} mm`]),
    ],
  },
  {
    name: 'alpha-warm-plate',
    check: (m) =>
      m.brightest === 0
        ? []
        : [`the frame over a -15 C plate is not black: brightest channel value ${m.brightest}`],
  },
  {
    name: 'muon',
    check: (m) => [
      ...within('muon track length (m)', m.lengthM, 0.2 * Math.hypot(1, 0.3), 0.03),
      ...(m.residualM < 6e-4
        ? []
        : [`muon track is not straight: rms ${(m.residualM * 1e3).toFixed(2)} mm`]),
    ],
  },
  {
    name: 'electron',
    check: (m) => [
      ...(m.lit > 20 ? [] : [`electron track barely visible: ${m.lit} lit pixels`]),
      ...(m.lengthM < electronRange / m.gasDensity / 100
        ? []
        : [`electron track longer than its path length: ${m.lengthM.toPrecision(3)} m`]),
    ],
  },
  {
    name: 'scene',
    check: (m) => [
      ...(m.std > 0.02 ? [] : [`flat scene: luma std ${m.std.toFixed(4)}`]),
      ...(m.mean < 0.5 ? [] : [`washed-out scene: mean luma ${m.mean.toFixed(3)}`]),
      ...(m.changed > 0.002
        ? []
        : [`scene does not move: ${(m.changed * 100).toFixed(2)}% of pixels changed in a second`]),
      ...(m.droplets > 5000 ? [] : [`too few droplets in the scene: ${m.droplets}`]),
    ],
  },
]

const report: { name: string; measurement: Measurement; problems: string[]; ms: number }[] = []
let failed = false
try {
  const tab = await browser.newPage({
    viewport: { width: size, height: size },
    deviceScaleFactor: 1,
  })
  const messages: string[] = []
  tab.on('console', (m) => {
    if (m.type() !== 'error' && m.type() !== 'warning') return
    if (m.text().includes('GPU stall due to ReadPixels')) return // the readback itself
    messages.push(`${m.type()}: ${m.text()}`)
  })
  tab.on('pageerror', (e) => messages.push(`pageerror: ${e.message}`))
  await tab.goto(`http://127.0.0.1:${server.port}/`)
  await tab
    .waitForFunction(() => typeof window.runScenario === 'function')
    .catch((error) => {
      throw new Error(
        `page never became ready: ${messages.join(' | ') || (error as Error).message}`,
      )
    })
  await mkdir(outDir, { recursive: true })
  for (const scenario of scenarios) {
    messages.length = 0
    const started = performance.now()
    let measurement: Measurement = {
      lit: 0,
      brightest: 0,
      lengthM: 0,
      residualM: 0,
      mean: 0,
      std: 0,
      gasDensity: 1,
      droplets: 0,
      changed: 0,
    }
    const problems: string[] = []
    try {
      measurement = await tab.evaluate(([n, s]) => window.runScenario(n, s), [
        scenario.name,
        size,
      ] as const)
      await tab.locator('#c').screenshot({ path: join(outDir, `${scenario.name}.png`) })
      problems.push(...scenario.check(measurement))
    } catch (error) {
      problems.push(`scenario failed: ${(error as Error).message.split('\n')[0]}`)
    }
    problems.push(...messages)
    if (problems.length) failed = true
    const ms = Math.round(performance.now() - started)
    report.push({ name: scenario.name, measurement, problems, ms })
    console.log(
      `${scenario.name.padEnd(17)} lit=${String(measurement.lit).padStart(6)} max=${String(measurement.brightest).padStart(3)} length=${(measurement.lengthM * 100).toFixed(2).padStart(6)}cm rms=${(measurement.residualM * 1e3).toFixed(2)}mm std=${measurement.std.toFixed(3)} droplets=${String(measurement.droplets).padStart(6)} ${String(ms).padStart(5)}ms  ${problems.length ? `FAIL ${problems.join(' | ')}` : 'ok'}`,
    )
  }
} finally {
  await browser.close()
  server.stop(true)
}
await writeFile(join(outDir, 'report.json'), `${JSON.stringify(report, null, 2)}\n`)
console.log(
  `${report.length} scenarios, ${report.filter((r) => r.problems.length).length} failing -> ${outDir}`,
)
console.log(
  `reference: alpha 5.5 MeV projected range ${alphaRange} g/cm^2 (NIST ASTAR); electron 30 keV CSDA range ${electronRange} g/cm^2 (NIST ESTAR)`,
)
process.exit(failed ? 1 : 0)
