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
    collapses. Hidden-only steps use bisector images as fold lines (a search restriction).
  - `Pattern.options(p, k, c)` returns the options at vertex p with k hidden rays and c
    collapses.
- **Global model.** `Model` builds the z3 problem.
  - Times are integer `t(key)`. Keys are crease ids for visible creases and `((p, r), 1|2)`
    for hidden ones.
  - Each vertex chooses one option through a `sel[p, i]` Boolean.
  - `path` traces cancelled segments: bends at earlier creases, `MAXBEND` = 2, joins at
    vertices. Given a model, it returns only the branch the model takes; `geometry` uses this
    to pin every segment's path.
- **High-degree vertices.** A vertex with more than `LAZY_T` (2000) histories becomes a
  `'LAZY'` option. `Model.symbolic` encodes the wedge and parity conditions, as bit-vectors
  when N is a power of two. `refine` and `first_failure` add collision conflicts lazily.
- **Search.** `solve` does core-guided expansion over the levels (c, k), c <= 1, k <= 4, in
  lexicographic order: vertices in an unsat core get their next level. It then minimizes
  `cost` = hidden rays + `W_COLLAPSE` * collapse cost over the options given so far, and
  reports `lower_bound` (sum of first levels) and `certified_minimal` (cost == lower bound).
- **Analysis.** `analyse` finds the forced-together folds and the forced order with the
  geometry fixed: `Checker` fixes each vertex's option and (by default) pins every cancelled
  segment's path, and re-validates lazy vertices on every query. `report` writes the figures
  (matplotlib fallback for the Hasse diagram when Graphviz is missing).
- **Exact local search (src/exact_local.py).** Decides whether one vertex has a local history
  with k hidden rays in arbitrary directions: fixed cyclic order + step sequence makes (L),
  flatness and (C) linear (z3 over reals and integers), DFS over steps with pruning.
- **Lookup tables (src/lookup.py, branch lookup-tables).** `build n k DIR` enumerates history types
  (symbolic visible directions, no hidden-only steps) into `DIR/deg{n}_k{k}.json`; `build n k DIR
  --scaffold` (branch scaffold-lines) the types with a scaffold step, i.e. a hidden-only step whose
  line later steps pin (`eliminate`), into `DIR/deg{n}_k{k}_s.json`. `Table` and `options_for`
  replace `Pattern.options` at runtime (`--table=DIR`); `--rigid --cells=DIR` filters by the exact
  cells of src/cells.py; `--dl` runs the global search with z3's difference-logic arithmetic.
  See docs/lookup_tables.md.
- **Tools.** `tools/common.py`: `resolve(path, geometry='found'|'shortest')`, `geometries`
  (all minimal geometries, paths included).

## Known warts
- The snapping branch of `Pattern.__init__` used when no N fits is convoluted and has not
  been tested in isolation.
- The symbolic encoding with real arithmetic (no N) is very slow in z3; it is untested on
  real cases.
- `KINDS` is a global counter used for diagnostics.
- Search restrictions (the paper, Section 8.5): cancelled segments start at visible vertices
  and bend at most twice. Passing through a visible vertex uses only a sufficient condition.
  Hidden-only steps use bisector images. At most 4 hidden rays (pairwise distinct directions)
  and at most one collapse per vertex. Minimality is over the options given during the search
  unless `certified_minimal`. No layer order is computed.

## Running things
- No Python on the owner's machine has the dependencies; use a venv built from the system
  Python (`/usr/bin/python3 -m venv`; the pyenv builds lack `bz2`, which networkx needs).
- The owner asked to keep CPU use low: run one solver process at a time, under `nice -n 19`.

## Quick checks
- `python3 tools/verify_examples.py` takes about a minute. It must report 15 folds,
  30 relations and 19 Hasse edges for the running pattern, and 6 folds and 24 linear
  extensions for the body pattern.
- `python3 src/fold_history.py --triage examples/komatsu_gentle_dragon.fold` must report
  24 hidden rays.
