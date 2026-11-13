/**
 * Opens the built demo page in headless Chromium and uses it the way a reader would: lets it
 * run, fires each particle, drags the view, and warms the plate until the tracks stop. Fails on
 * any console error, on a black canvas, or if the page's own readout disagrees with the solver.
 *
 *   bun run tools/check-demo.ts [--page <built html>] [--out renders]
 */
import { mkdir } from 'node:fs/promises'
import { join, resolve } from 'node:path'
import { chromium } from 'playwright-core'
import { DEFAULT_CHAMBER, solveChamber } from '../src/chamber'

function arg(name: string, fallback: string): string {
  const i = process.argv.indexOf(`--${name}`)
  return (i >= 0 ? process.argv[i + 1] : undefined) ?? fallback
}
const pagePath = resolve(arg('page', join(import.meta.dir, '..', 'public', 'index.html')))
const outDir = resolve(import.meta.dir, '..', arg('out', 'renders'))
const html = await Bun.file(pagePath).text()
const server = Bun.serve({
  hostname: '127.0.0.1',
  port: 0,
  fetch: () => new Response(html, { headers: { 'content-type': 'text/html' } }),
})
const browser = await chromium.launch({
  headless: true,
  args: ['--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--ignore-gpu-blocklist'],
})
const problems: string[] = []
try {
  const tab = await browser.newPage({
    viewport: { width: 1100, height: 900 },
    deviceScaleFactor: 1,
  })
  tab.on('console', (m) => {
    if (
      (m.type() === 'error' || m.type() === 'warning') &&
      !m.text().includes('GPU stall due to ReadPixels')
    )
      problems.push(`${m.type()}: ${m.text()}`)
  })
  tab.on('pageerror', (e) => problems.push(`pageerror: ${e.message}`))
  await tab.goto(`http://127.0.0.1:${server.port}/`)
  await mkdir(outDir, { recursive: true })

  const litFraction = () =>
    tab.evaluate(() => {
      const canvas = document.getElementById('view') as HTMLCanvasElement
      const gl = canvas.getContext('webgl2')
      if (!gl) return -1
      const px = new Uint8Array(canvas.width * canvas.height * 4)
      gl.readPixels(0, 0, canvas.width, canvas.height, gl.RGBA, gl.UNSIGNED_BYTE, px)
      let lit = 0
      for (let i = 0; i < px.length; i += 4) if (px[i + 1] > 40) lit++
      return lit / (px.length / 4)
    })

  // 1. It runs: after a few seconds there are droplets on screen and the readouts are filled in.
  await tab.waitForTimeout(4000)
  const running = await litFraction()
  if (running < 0.002)
    problems.push(`canvas is nearly black after 4 s: ${(running * 100).toFixed(2)}% lit`)
  const expected = solveChamber(DEFAULT_CHAMBER)
  const label = (await tab.locator('#layer').textContent()) ?? ''
  const want = `${(expected.sensitiveBottom * 1000).toFixed(0)} to ${(expected.sensitiveTop * 1000).toFixed(0)} mm`
  if (!label.includes(want)) problems.push(`layer readout "${label}" does not contain "${want}"`)
  const count = (await tab.locator('#count').textContent()) ?? ''
  if (!/\d/.test(count)) problems.push(`droplet count not shown: "${count}"`)

  // 2. Buttons fire particles: with background and mist off, the canvas goes dark, then a fired
  // alpha lights it.
  await tab.locator('#background').uncheck()
  await tab.locator('#mist').uncheck()
  await tab.waitForTimeout(300)
  const dark = await litFraction()
  await tab.locator('#fire-alpha').click()
  await tab.locator('#fire-muon').click()
  await tab.locator('#fire-electron').click()
  await tab.waitForTimeout(250)
  const fired = await litFraction()
  if (!(fired > dark + 0.0002)) problems.push(`firing particles lit nothing: ${dark} -> ${fired}`)
  await tab.locator('#view').screenshot({ path: join(outDir, 'demo-fired.png') })

  // 3. Dragging orbits the camera without errors.
  const box = await tab.locator('#view').boundingBox()
  if (box) {
    await tab.mouse.move(box.x + 200, box.y + 200)
    await tab.mouse.down()
    await tab.mouse.move(box.x + 320, box.y + 240, { steps: 5 })
    await tab.mouse.up()
  }

  // 4. Warm the plate to -15 C: the page must say there is no sensitive layer, and a fired alpha
  // must leave nothing.
  await tab.locator('#plate').fill('-15')
  await tab.waitForTimeout(200)
  const warm = (await tab.locator('#layer').textContent()) ?? ''
  if (!warm.startsWith('No sensitive layer: the plate is not cold enough'))
    problems.push(`at -15 C the page says: "${warm}"`)
  await tab.locator('#fire-alpha').click()
  await tab.waitForTimeout(1500) // past the ion lifetime, so any droplet that will form has formed
  const warmLit = await litFraction()
  if (warmLit > dark + 0.0002)
    problems.push(`a track appeared over a -15 C plate: ${(warmLit * 100).toFixed(3)}% lit`)

  // 5. Back to the default with everything on, for a picture of the whole page.
  await tab.locator('#plate').fill('-40')
  await tab.locator('#background').check()
  await tab.locator('#mist').check()
  await tab.locator('#americium').check()
  await tab.waitForTimeout(5000)
  await tab.screenshot({ path: join(outDir, 'demo-page.png'), fullPage: true })
  console.log(
    `lit: running ${(running * 100).toFixed(2)}%, dark ${(dark * 100).toFixed(3)}%, fired ${(fired * 100).toFixed(3)}%, warm plate ${(warmLit * 100).toFixed(3)}%`,
  )
} finally {
  await browser.close()
  server.stop(true)
}
for (const p of problems) console.error(`FAIL ${p}`)
console.log(problems.length ? `${problems.length} problem(s)` : 'demo page ok')
process.exit(problems.length ? 1 : 0)
