/**
 * Builds the interactive page into one self-contained HTML file.
 *
 *   bun run tools/build-demo.ts [--out <path>]   default: public/index.html
 */
import { join, resolve } from 'node:path'

const i = process.argv.indexOf('--out')
const out = resolve(
  i >= 0 ? (process.argv[i + 1] ?? '') : join(import.meta.dir, '..', 'public', 'index.html'),
)

const build = await Bun.build({
  entrypoints: [join(import.meta.dir, '..', 'demo', 'app.ts')],
  target: 'browser',
  format: 'esm',
  minify: true,
})
if (!build.success) {
  for (const log of build.logs) console.error(log)
  process.exit(1)
}
const [artifact] = build.outputs
if (!artifact) throw new Error('bundle produced no output')
const script = (await artifact.text()).replaceAll('</script', '<\\/script')
const page = await Bun.file(join(import.meta.dir, '..', 'demo', 'page.html')).text()
if (!page.includes('/*BUNDLE*/')) throw new Error('demo/page.html has no /*BUNDLE*/ marker')
const html = page.replace('/*BUNDLE*/', () => script)
await Bun.write(out, html)
console.log(`${out} (${(html.length / 1024).toFixed(0)} KB)`)
