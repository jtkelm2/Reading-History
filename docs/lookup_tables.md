# Lookup tables of local histories

The solver needs, at every interior vertex, the list of its local histories. It can search for them
at runtime (`Pattern.options`, the method of `paper/reading_history.tex`), or look them up in tables
computed once per vertex degree (`src/lookup.py`). This document describes the tables, the exact
realizability cells that filter them by what paper can do, how both were checked, and what they
give on the corpus.

```
python3 src/fold_history.py --table=tables pattern.fold                                  # tables
python3 src/fold_history.py --table=tables --rigid=reverse --cells=cells --dl pattern.fold   # physical histories only
```

- `--rigid=reverse` keeps only histories realizable with rigid half-turns and reverse folds;
  `--rigid=halfturn` keeps half-turns only.
- `--cells=cells` takes that realizability from the precomputed cells instead of deciding it at
  runtime.
- `--dl` runs the global search with z3's dense difference-logic arithmetic.

## Terms

| term here and in the code | in the encyclopedia and `paper/flatvertex.tex` | meaning |
|---|---|---|
| visible ray | crease | a crease of the pattern at the vertex |
| hidden ray | scaffold crease | a ray that a history makes and later erases |
| scaffold step | free line | a step whose rays are all new hidden rays, so no existing ray fixes its line |
| physical | physical | realizable rigidly, or with reverse folds, in the layered model of `paper/flatvertex.tex` |

## History types

For a degree n and a number k of hidden rays, `python3 src/lookup.py build n k tables` enumerates
every *history type* valid at some flat vertex of degree n. A history type consists of:
- a cyclic order of the n visible and k hidden rays;
- a sequence of single-line steps, each making rays or erasing hidden ones.

The visible directions are symbolic (a₀ = 0 < a₁ < … < a_{n−1}). The enumeration is a depth-first
search over steps that keeps a prefix exactly when z3 finds some flat vertex and some hidden
directions satisfying it: (P), (L) through Theorem 11, flatness of every state, (C), and the cyclic
order.

Restrictions:
1. at most k hidden rays, pairwise distinct and distinct from the visible rays;
2. no collapses;
3. degree n ≤ 6, so vertices of higher degree get no options;
4. no free hidden directions (next section).

Rotations are normalised: visible ray 0 is in the first step that contains a visible ray. The
runtime tries every rotation. Mirror images are separate types.

### Anchored types and scaffold types

In an *anchored* step, some ray (a visible ray, or a hidden ray being erased) fixes the fold line,
and the step's new hidden rays are preimages of that line. When every step is anchored, every
hidden direction is an integer combination of the visible directions mod π. Every condition (L)
and every flatness condition is then a congruence c·a ≡ 0 (mod π). These congruences are the
type's *index*. Anchored types are in `tables/deg{n}_k{k}.json`.

A *scaffold step* makes only hidden rays. `lookup.py build n k tables --scaffold` enumerates the
types with at least one scaffold step into `tables/deg{n}_k{k}_s.json`. The direction (mod π) of a
scaffold step's first ray is a new variable λ_j, and its other rays follow from it. Every condition
is then a congruence c·a + d·λ ≡ 0 (mod π) with integers c, d. `eliminate` removes the λ_j by
unimodular integer row operations:
- **Pivots.** For each line j in turn, Euclid's algorithm on column j leaves at most one row with a
  nonzero entry there, its *pivot row*, which has no entries in the earlier columns.
- **Unpinned lines are dropped.** If some column has no pivot, that line can move continuously
  without violating any congruence: a free hidden direction. Such types are counted as `unpinned`
  and dropped.
- **Index.** Otherwise pivot row j reads d_jj λ_j ≡ r_j(a, λ_{j+1}, …) (mod π), which has exactly
  |d_jj| solutions for every right-hand side. So the vertices at which some lines satisfy all
  congruences are exactly those satisfying the rows without lines; these rows are the index.
- **Branches.** A branch chooses u_j ∈ [0, |d_jj|) for each line, giving
  λ_j = (r_j + u_j π)/d_jj. Here r_j substitutes the later lines unreduced, so that λ_j is affine in
  a and the branch labels do not jump where a line crosses 0 ≡ π.

### Runtime

For each rotation of the vertex, the runtime evaluates the index congruences, keeps the types whose
congruences all hold, and replays each kept type at the vertex. For a scaffold type it replays every
branch. The replay checks the strict inequalities and (C) and computes the hidden directions. There
is no search over steps. Lookups are computed per vertex and per depth on demand, and cached per
vertex shape.

## Exact realizability cells

`src/cells.py` computes, for every type, a decision tree over the vertices: `cells/` holds one file
per table file. Each leaf has a value: 2 (rigidly realizable), 1 (realizable with reverse folds
only) or 0 (neither). Realizability is in the layered model of `paper/flatvertex.tex` (`src/rigid.py`).

**Statement.** Write a vertex as b = (b_1, …, b_{n−1}), its directions in units of π. Then the set
of vertices at which a type is valid and realizable is a finite union of relatively open polyhedra,
each cut out by sign conditions on finitely many affine forms in b with rational coefficients.

**Why.** Replaying a type and deciding its realizability (`src/rigid.py`, through `src/geom.py`)
consults the vertex only through signs of affine forms:
- every normalisation mod 2π is an integer certified by two signs;
- equality, betweenness, sorting and midpoints reduce to signs.

Within a region where every consulted sign is fixed, the computation therefore returns the same
outcome. Splitting a region by the feasible signs of each undecided form gives a finite decision
tree.

**How.**
1. **Pieces.** A type's index congruences and Kawasaki's condition cut out one affine subspace for
   each choice of integer values. That subspace is clipped to the order constraints and kept if it
   keeps its dimension, which is 1 for every degree-4 type and at most 2 for degree 6. Scaffold
   types have one tree per (piece, branch).
2. **Regions.** Regions are relatively open convex polygons, stored by the vertices of their
   closure in exact integer arithmetic. A form's sign on a region is read off those vertices, and
   an undecided form splits the region by exact clipping.
3. **Symbolic execution.** The replay and `rigid.skeleton` run on exact affine forms. Only the
   geometry is symbolic: given the skeletons, the choice of labelings, up sides and senses is
   combinatorial.

Scaffold lines are affine in b with denominators dividing the product of the pivots. A scaffold
table's forms therefore have numerators over `den` = 4 · lcm(products of pivots), stored in the
cells meta (between 8 and 240).

## Builds

One core, z3 4.16, `nice -n 19`. Every build is complete.

| degree | k | anchored types | scaffold types (unpinned, dropped) | table build | cell leaves | cells build |
|---|---|---|---|---|---|---|
| 2 | 0–2 | 1 / 0 / 0 | – / 0 / 1 (– / 0 / 2) | < 1 s | 1 + 2 | < 1 s |
| 4 | 0 | 1 | – | 0.0 s | 3 | 0.0 s |
| 4 | 1 | 2 | 0 | 0.1 s | 8 | 0.0 s |
| 4 | 2 | 12 | 12 (12) | 1.4 + 1.3 s | 114 + 90 | 0.4 + 0.4 s |
| 4 | 3 | 82 | 156 (66) | 19.7 + 41.9 s | 1,310 + 1,860 | 6.1 + 8.7 s |
| 4 | 4 | 556 | 2,754 (1,619) | 303.8 + 1,582.0 s | 12,938 + 60,984 | 72.0 + 294.9 s |
| 6 | 0 | 4 | – | 0.1 s | 68 | 0.4 s |
| 6 | 1 | 54 | 0 | 4.5 s | 1,422 | 11.9 s |
| 6 | 2 | 688 | 282 (152) | 144.5 + 62.0 s | 32,058 + 12,164 | 341.2 + 125.4 s |
| 6 | 3 | 9,476 | not built | 4,528.4 s | 675,096 | 9,249.7 s |

Where a cell holds two numbers, they are for the anchored and scaffold files. Degree-6 scaffold
types with k = 3 are not built; at the anchored cost they would take several hours.

## Checks

| check | what it compares | result |
|---|---|---|
| `tools/validate_table.py tables PATTERNS` | anchored types against the runtime local search at the 66 vertex shapes of degree ≤ 6 in the corpus and examples, every depth | 299 of 299 comparisons equal |
| `tools/validate_scaffold.py tables PATTERNS` | scaffold types against the local search's scaffold histories whose lines are pinned (`tools/pinning.py`), 63 shapes, 257 shape × depth comparisons | all 2,259 search histories are in the tables (5,400); every table history is pinned by the numeric measure |
| `tools/verify_cells.py CELLS TABLE PATTERNS` | each cell tree against the numeric computation, at every leaf's witness vertex and every rotation of every corpus vertex shape | 0 disagreements; anchored tables: 723,018 witnesses, 4,362,610 pairs |
| `tools/verify_flatvertex.py` | every computed claim of `paper/flatvertex.tex` against `src/rigid.py` | all hold |
| `tools/grid_search.py N RAYS --reverse` | independent breadth-first search over layered states, fold lines on a grid | finds the 5-step scaffold histories that the tables give at (22.5°, 90°, 157.5°, 90°) and (67.5°, 90°, 112.5°, 90°) |
| `tools/verify_examples.py --table=tables` | the paper's numbers for its two worked examples, from the tables | equal |

## Corpus results

The 48 patterns of `tools/curated_corpus.txt`, at 170 s each on one core. The runs are in
`results/lookup/`, made with `tools/run_corpus.sh tools/curated_corpus.txt OUT 170 1 "FLAGS"`.

| run | options | resolved | local fail | global unsat | timeout | time of resolved: median / max |
|---|---|---|---|---|---|---|
| `search` | local search (no tables) | 41 | 0 | 0 | 7 | 0.1 / 59.1 s |
| `tables_anchored` | `--table` (anchored types only) | 29 | 18 | 0 | 1 | 0.2 / 94.5 s |
| `tables_dl` | `--table --dl` | 30 | 17 | 0 | 1 | 0.1 / 49.9 s |
| `rigid_anchored` | `--rigid=halfturn --cells` (anchored types only) | 6 | 41 | 1 | 0 | 0.2 / 0.3 s |
| `physical_anchored` | `--rigid=reverse --cells` (anchored types only) | 23 | 22 | 2 | 1 | 0.3 / 4.6 s |
| `physical` | `--rigid=reverse --cells` | 28 | 17 | 0 | 3 | 0.2 / 41.5 s |
| `physical_dl` | `--rigid=reverse --cells --dl` | 30 | 17 | 0 | 1 | 0.2 / 52.9 s |

The table is printed by `python3 tools/lookup_summary.py results/lookup/*.jsonl`.

What these runs show:
- **The physical filter.** Many cheapest local histories are not physically realizable. Requiring
  reverse folds to be available raises costs on some patterns. With anchored types only, it loses
  six of them:
  - Dutch Rabbit, Squirrel, Kirby and kitsune have vertices of shape (22.5°, 90°, 157.5°, 90°) or
    (67.5°, 90°, 112.5°, 90°), with no physical anchored history;
  - Color-change Pinwheel and Marnie have local options everywhere, but no consistent global
    resolution.

  With half-turns only, few patterns resolve, so reverse folds are needed almost wherever there is
  a cancellation.
- **Scaffold types.** All six patterns resolve with physical scaffold histories. The two vertex
  shapes above are reached by 5-step histories that begin with a free line; for example, at
  {0, 1, 5, 12} × 22.5° the steps are {2, 10}, {5, 15}, {12}/✕2, {0}/✕10, {1}/✕15. Each such line
  is pinned by a later step. 054 Hexagon, which no anchored run resolves, resolves too.
- **The price of being physical.** With all types and `--dl`, the filtered and unfiltered runs
  resolve the same 30 patterns. The filter raises the number of hidden rays on 8 of them and
  lowers it on none: 048 (18 to 19), 070 (11 to 12), 073 (22 to 24), 114 (16 to 17),
  141 (17 to 24), 292 (5 to 7), 303 (6 to 7) and 315 (19 to 37).
- **`--dl`.** z3's difference-logic arithmetic is what lets kitsune and Marnie finish within the
  budget.

**Realizability survey** (`tools/rigid_survey.py tables [--reverse] PATTERNS`;
`results/lookup/survey_*.jsonl`). This covers the anchored types at the 66 vertex shapes of
degree ≤ 6 in the corpus and examples:

| hidden rays k | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| histories in the tables | 54 | 111 | 608 | 3,418 | 655 |
| realizable with half-turns only | 49 | 0 | 10 | 0 | 0 |
| realizable with reverse folds | 54 | 76 | 301 | 866 | 130 |

The reason is visible at v\* of the running pattern. Before its third step, the two faces that must
move together are, in every reachable stack, either the middle two layers or the top and bottom
layers: an inside or an outside reverse fold. A half-turn moves a block of layers above or below
the staying ones, so it can do neither.

## Limits

- **Degree.** No tables beyond degree 6, so vertices of degree ≥ 8 get no options.
- **Generic vertices.** Every type needs at least one condition on the angles. Counting argument:
  4 + k rays must be made or erased and each of at most 2 + k steps fixes one line. That leaves at
  least two conditions, one of which flatness can make redundant. Computed: every degree-4 type
  lives on a line or a point of the square of flat vertices. A degree-4 vertex in general position
  therefore has no option in the tables, although it folds rigidly by a collapse of all four
  creases. Collapse steps are not in the tables.
- **Minimality.** Costs are minima over the options given during core-guided expansion; see
  `certified_minimal` in each run.

The encyclopedia (`docs/encyclopedia/`) displays the degree-4 tables and cells stratum by stratum.
