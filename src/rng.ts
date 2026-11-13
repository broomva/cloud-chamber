/** Seeded random numbers, so a run is reproducible from its seed. */

// CHAMBER_SEED_SHIFT moves every seed by a constant. `bun run test:seeds` uses it to check that
// no test passes only on the seed it was written with.
const SEED_SHIFT = Number(
  (globalThis as { process?: { env?: Record<string, string | undefined> } }).process?.env
    ?.CHAMBER_SEED_SHIFT ?? 0,
)

export class Rng {
  private s: number
  constructor(seed: number) {
    this.s = (seed + SEED_SHIFT) >>> 0 || 0x9e3779b9
  }

  /** Uniform in (0, 1): mulberry32, never exactly 0 so a logarithm of it is finite. */
  next(): number {
    this.s = (this.s + 0x6d2b79f5) >>> 0
    let t = this.s
    t = Math.imul(t ^ (t >>> 15), t | 1)
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61)
    return (((t ^ (t >>> 14)) >>> 0) + 0.5) / 4294967296
  }

  /** Exponential with the given mean. */
  exp(mean: number): number {
    return -mean * Math.log(this.next())
  }

  /** Standard normal (Box-Muller). */
  normal(): number {
    return Math.sqrt(-2 * Math.log(this.next())) * Math.cos(2 * Math.PI * this.next())
  }

  /** Poisson count. Knuth's product for small means, a rounded normal above 30. */
  poisson(mean: number): number {
    if (mean <= 0) return 0
    if (mean > 30) return Math.max(0, Math.round(mean + Math.sqrt(mean) * this.normal()))
    const limit = Math.exp(-mean)
    let n = 0
    let p = this.next()
    while (p > limit) {
      n++
      p *= this.next()
    }
    return n
  }
}
