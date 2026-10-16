# Notes for Claude Code

## What this project is for
The research question is how much of a crease pattern's folding history its geometry forces.
The owner's main interest is *educational analysis of practical crease patterns*: helping
beginners read a pattern top-down the way experienced folders do. That means spotting the
earliest structural creases, the base or design paradigm, and gadgets such as grafts, point
splits and level shifters.

The current code answers the bottom-up half: folds and the forced order. `docs/ROADMAP.md`
describes the top-down half. Prefer general, principled theory over collections of special
cases.

## Working conventions
- The owner is a mathematician working in formal methods, with high standards of proof, and
  prefers minimal, precise answers.
- State every claim exactly. Separate what is proved from what was computed, and say which
  search restrictions a computed claim depends on.
- Verify every number that goes into the paper by running code.
- The paper must read as a standalone document. Do not describe how the method or the paper
  used to be.

## Code map (src/fold_history.py)
- **Directions.** Directions are reals in [0, 2pi), compared with `eq`/`EPS`. `Pattern`
  clusters crease lines within `TOL` and snaps them to multiples of 2pi/N when some
  N <= `MAXN` (64) fits.
- **Local histories.**
  - `maps` gives the folding maps of the sectors at a vertex; `collision_free` is condition (C).
  - `histories` enumerates histories using visible creases only.
  - `lazy_histories` adds hidden rays and collapse steps; `ctarget` sets the number of
    collapses.
  - `Pattern.options(p, k, c)` returns the options at vertex p with k hidden pairs and c
    collapses.
- **Global model.** `Model` builds the z3 problem.
  - Times are integer `t(key)`. Keys are crease ids for visible creases and `((p, r), 1|2)`
    for hidden ones.
  - Each vertex chooses one option through a `sel[p, i]` Boolean.
  - `path` traces cancelled segments: bends at earlier creases, `MAXBEND` = 2, joins at
    vertices.
- **High-degree vertices.** A vertex with more than `LAZY_T` (2000) histories becomes a
  `'LAZY'` option. `Model.symbolic` encodes the wedge and parity conditions, as bit-vectors
  when N is a power of two. `refine` and `first_failure` add collision conflicts lazily.
- **Search.** `solve` does core-guided expansion: vertices in an unsat core climb the
  `ladder` (more hidden pairs, then collapses). It then minimizes `cost` = hidden rays
  + `W_COLLAPSE` * collapse cost.
- **Analysis.** `analyse` finds the forced-together folds and the forced order; `Checker`
  re-validates lazy vertices on every query. `report` writes the figures.

## Known warts
- `Pattern.candidates` contains a leftover no-op (`or True`). The snapping branch used when
  no N fits is convoluted and has not been tested in isolation.
- The symbolic encoding with real arithmetic (no N) is very slow in z3; it is untested on
  real cases.
- `KINDS` is a global counter used for diagnostics.
- Search restrictions (the paper, Section 8.5): cancelled segments start at visible vertices
  and bend at most twice. Passing through a visible vertex uses only a sufficient condition.
  At most 4 hidden pairs and at most one collapse per vertex. No layer order is computed.

## Quick checks
- `python3 tools/verify_examples.py` takes about a minute. It must report 15 folds,
  30 relations and 19 Hasse edges for the running pattern, and 6 folds and 24 linear
  extensions for the body pattern.
- `python3 src/fold_history.py --triage examples/komatsu_gentle_dragon.fold` must report
  24 hidden rays.
