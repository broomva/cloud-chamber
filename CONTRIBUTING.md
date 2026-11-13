# Contributing

Issues and pull requests are welcome, especially corrections to the physics.

## Setup

```bash
bun install            # Bun 1.3.14 (CI pins it; the committed demo page is byte-identical only on it)
bun run typecheck && bun run lint && bun run test
bunx playwright-core install chromium-headless-shell   # once, for the next two
bun run render         # headless WebGL2 gate
bun run demo && bun run check-demo
```

## The one rule

Every claim the simulation makes is either **checked against a published reference** or
**labelled as a model prediction** in the README. A change to the physics comes with a test that
would fail without it:

- If a reference exists (NIST, PDG, IAEA, CRC/HSDB), compare against it and cite it in the test.
- If none exists, test conservation laws and consistency, and say in the README that the result
  is a prediction.
- Known biases are pinned from both sides: a test fails if the model gets better or worse
  without anyone noticing.

Monte Carlo tests must hold under `bun run test:seeds`, which reruns the suite with eight seed
shifts. A test that only passes at seed 1 is not a test.

## Pull requests

- Run the checks above. CI runs the same ones.
- If you change `demo/`, `src/render.ts` or `tools/build-demo.ts`, run `bun run demo` and commit
  the regenerated `public/index.html`. CI rebuilds it and fails on any difference.
- Say what the change does to the numbers printed by `bun run report`.

Be kind in issues and reviews. Disagreement about physics is welcome; disrespect is not.
