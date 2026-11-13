/** Reference tables committed under test/reference, each with its fetch URL and date in the header. */
import { join } from 'node:path'

export async function loadTable(name: string): Promise<number[][]> {
  const text = await Bun.file(join(import.meta.dir, 'reference', name)).text()
  return text
    .split('\n')
    .filter((line) => line.trim() && !line.startsWith('#'))
    .map((line) => line.trim().split(/\s+/).map(Number))
}

/** Relative difference |a - b| / |b|. */
export const rel = (a: number, b: number) => Math.abs(a - b) / Math.abs(b)
