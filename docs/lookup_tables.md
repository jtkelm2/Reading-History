# Lookup tables of local histories: first build

Branch `lookup-tables`. To roll back: `git checkout overnight` (the state before this work,
commit `93ad7f3`). Without `--table`, the solver behaves exactly as on `overnight`.

## What the tables contain

For each degree n and number k of hidden rays, `src/lookup.py build n k` enumerates every
*history type* that is valid at some flat vertex of degree n:

- a cyclic order of the n visible rays and the k hidden rays, and
- a sequence of single-line steps, each making visible or hidden rays or erasing hidden ones.

The visible directions are symbolic (a₀ = 0 < a₁ < … < a_{n−1}). The enumeration is a
depth-first search over steps that keeps a prefix exactly when z3 finds some flat vertex, and some
hidden directions, satisfying it: (P), (L) through Theorem 11, flatness of every state, (C), and
the cyclic order. It is exhaustive within these restrictions:

1. **No free hidden directions.** A step whose rays are all new hidden rays is excluded. Every
   other step contains a visible ray being made or a hidden ray being erased, whose image fixes
   the fold line, so every hidden direction is an integer combination of visible directions
   modulo π.
2. **Depth.** At most k hidden rays per vertex, with k as in the table of builds below; hidden
   directions pairwise distinct and distinct from the visible ones (as in the search).
3. **Degree.** n ≤ 6. Vertices of higher degree have no entries.
4. **No collapses, and (C) rather than rigid realizability** (see below).

By rotation symmetry the first step contains visible ray 0; the runtime tries every rotation.

Each type records the congruences c·a ≡ 0 (mod π) that it needs, with hidden directions
substituted. Since hidden directions are a + jπ for integer combinations a, every condition (L)
and every flatness condition reduces to such a congruence in the visible directions. These are
the equation part of the type's cell; the rest of the cell (the cyclic order and (C)) is linear
inequalities.

## Runtime

For a vertex, for each rotation of its rays: evaluate the finite list of congruence vectors at the
vertex, take the types whose required congruences all hold, and verify each by replaying its
steps at the vertex, which checks the inequalities and (C). No search over steps happens at
runtime. Entries are cached per vertex shape. A vertex the tables do not cover gets no options,
and the solver reports `local fail` with that vertex's degree.

`python3 src/fold_history.py --table=tables pattern.fold` uses the tables;
`tools/run_corpus.sh LIST OUT SECONDS 1 --table=tables` runs the corpus with them.

## Why rigid realizability is not applied yet

Rigid realizability (flatvertex, Definition 2.2) is not a property of a history type alone. It
depends on the layer orders, which depend on the cyclic order of the images of the rays, which
changes inside a congruence cell. Applying it in the tables needs cells refined by the order of
image points (still linear conditions, but many more of them), or a per-vertex check over the
choices of labelling and up-side at every step, which is a search at runtime. Deferred.

## Tables built

| Degree | Hidden rays k | Types | Build time | Complete |
|---|---|---|---|---|
| 2 | 0 / 1 / 2 | 1 / 0 / 0 | < 1 s | yes |
| 4 | 0 / 1 / 2 / 3 / 4 | 1 / 2 / 12 / 82 / 556 | 0.0 / 0.1 / 1.4 / 19.7 / 303.8 s | yes |
| 6 | 0 / 1 / 2 / 3 | 4 / 54 / 688 / 9476 | 0.1 / 4.5 / 144.5 / 4528.4 s | yes |

Degree 2 has no types with hidden rays: its first step uses both visible rays, so no later step
has an anchor. All builds: z3 4.16, one core, `nice -n 19`.

## Validation

`tools/validate_table.py tables` compares, at every vertex shape of degree ≤ 6 in the corpus and the
paper's examples (66 shapes, each once up to rotation and reflection) and every depth built, the
tables' histories with those of the local search, restricted to histories with no hidden-only
step. Result: 299 of 299 comparisons equal (with the degree-6 table to k = 3). `tools/verify_examples.py --table=tables` reproduces the
paper's numbers for the running and body patterns.

## First results on the corpus

`results/lookup_first_build/`: `corpus_search.jsonl` (local search, as on `overnight`),
`corpus_tables.jsonl` (`--table=tables`), `comparison.md`. 170 s per pattern, one core.

- **Resolved:** search 41, tables 28. Every pattern the tables resolve, the search resolves too.
- **Cost on the 28 both resolve:** equal on 23; higher with the tables on 5 (016: 15 vs 13,
  048: 25 vs 19, 058: 17 vs 16, 070: 13 vs 9, 114: 18 vs 16).
- **Time on those 28:** search 106.9 s in total, tables 39.8 s. The local phase with the tables
  takes 0.04–0.16 s per pattern; the rest is the global z3 problem.
- **Why the tables fail (20 patterns):** 13 have a vertex of degree > 6 (10 of degree 8 first,
  and one each of degree 12, 16, 20); 6 have a vertex of degree ≤ 6 with no history within the
  tables (degree 4: 054, 084, 325; degree 6: 015, 093, 315); 304 times out in both.
- **Why the tables do worse** (the search's optimum inspected vertex by vertex): the optima of
  048, 058, 114 and 093 use histories with a hidden-only step, which the tables exclude, and those
  of 048 and 093 also a degree-6 vertex with 3 hidden rays (table depth 2); the optima of 054,
  084, 325 and 015 use collapses. The optimum of 016 uses only histories the tables contain: its
  higher cost with the tables comes from core-guided expansion taking a different path, and
  neither cost is certified minimal. 070 was not diagnosed (the search timed out in the
  diagnostic run, after 13 s in the corpus run).

## With the degree-6 table to k = 3

`results/lookup_first_build/corpus_tables_k3.jsonl`; `comparison.md` is regenerated from it.

- **Resolved:** 29. New: Marnie (315), cost 19 (not certified) in 12.6 s; the search times out
  on it. Dutch Rabbit (048) now reaches the search's cost, 19, in 2.0 s (search: 25.9 s).
- **Cost on the 28 resolved by both:** equal on 24; higher with the tables on 4 (016, 058,
  070, 114).
- **Time on those 28:** search 106.9 s, tables 113.7 s. The difference is 070, which takes
  94.5 s with the larger degree-6 table (5.7 s without k = 3): more options make the global
  core-guided problem slower.

## Addendum: degree 8 (k ≤ 1), outside the first build

`tables_deg8/` (degree 8: k = 0, 33 types, 2.9 s; k = 1, 1244 types, 249 s; both complete),
used together with `tables/`; `results/lookup_first_build/corpus_tables_deg8.jsonl`.

- **Resolved:** 34, five more than without degree 8: Pig (045), 4 Stripes 1 (091), Horned Owl
  (096), Fox (113), Sus (186). Their costs equal the search's except the Horned Owl's (21 vs 16).
- **Still failing:** 7 patterns with vertices of degree 12, 14, 16 or 20 (including the Gentle
  Dragon, degree 12); Butterfly Base (018), whose degree-8 vertex has no history with at most
  one hidden ray; the 5 patterns of degree ≤ 6 above; 304 (timeout in both).

## Rigid realizability of the tables' histories

`src/rigid.py` implements flatvertex.tex, Definitions 2.1–2.2. Checks against that note: H₁ is
rigidly realizable (the note expected this but had not checked it); H₂ and H₃ are not; after steps
1–3 of H₂ all 64 choice combinations survive and give exactly the eight layer orders in the proof
of Theorem 2.3. Erasing steps also pass positive tests (fold in half and unfold; fold twice, unfold
both).

Survey (`tools/rigid_survey.py tables …`, `results/lookup_first_build/rigid_survey.jsonl`): at the
66 vertex shapes of degree ≤ 6 in the corpus and examples, the histories the tables give, and how
many are rigidly realizable:

| hidden rays k | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| histories | 54 | 111 | 608 | 3418 | 655 |
| rigidly realizable | 49 | 0 | 10 | 0 | 0 |

22 of the 66 shapes have no rigidly realizable history within the tables, among them v\* of the
running pattern and p± of the body pattern; with `--rigid`, both examples end in `local fail`.

The reason, worked out at v\*: before its third step ({225° made, hidden 0° erased}), in every
reachable stack the two faces that must move together are either the middle two layers or the top
and bottom layers. That is an inside or an outside reverse fold. A rigid half-turn moves a block
of layers above or below the staying ones, so it cannot do either. Cancellation in the tables is
mostly of this kind, so rigidity as defined in flatvertex rejects it. Using rigidity as a filter
needs a layered model of reverse folds first.

## Reverse folds

`src/rigid.py` now gives every moving piece a sense (the side it is lifted over); each moving run
keeps the common sense or is split at one spine into halves with opposite senses. With all senses
equal this is flatvertex's rigid half-turn. The definitions, lemmas and computed examples are in
flatvertex.tex, Section 3 (`tools/verify_flatvertex.py` asserts every computed claim there).

Survey with reverse folds (`tools/rigid_survey.py tables --reverse …`,
`results/lookup_first_build/rigid_survey_reverse.jsonl`; 723 s):

| hidden rays k | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| histories in the tables | 54 | 111 | 608 | 3418 | 655 |
| realizable, half-turns only | 49 | 0 | 10 | 0 | 0 |
| realizable with reverse folds | 54 | 76 | 301 | 866 | 130 |

Shapes with no realizable history within the tables: 22 with half-turns only, 2 with reverse folds:
(22.5, 90, 157.5, 90), the Horned Owl's staircase vertex, and (67.5, 90, 112.5, 90). For 5 shapes the
least number of hidden rays rises. With `--table=tables --rigid=reverse`, the running and body
patterns resolve with the same folds and forced order as without the filter (15 folds, 30 relations,
19 Hasse edges; 6 folds, 9 relations, 24 linear extensions); the check takes 86 s and 26 s, because
realizability is decided at runtime by a search over labelings, up sides and senses.

## Exact realizability cells (refining the cells for rigidity)

`src/cells.py`, `cells/deg{n}_k{k}.json` (one decision tree per table type), checked by
`tools/verify_cells.py`.

**Statement.** For a history type t of degree n, write the vertex as b = (b_1, …, b_{n−1}), its
directions in units of π (b_0 = 0), with 0 < b_1 < … < b_{n−1} < 2 and Kawasaki's condition
b_1 − b_2 + … + b_{n−1} = 1. The set of vertices at which t is valid and realizable (rigidly, or
with reverse folds) is a finite union of relatively open polyhedra, each cut out by sign conditions
on finitely many affine forms in b with rational coefficients.

**Why.** Replaying t and deciding its realizability (src/rigid.py, through src/geom.py) consults
the vertex only through signs of affine forms: every normalisation mod 2π is an integer q chosen
at a point of the current region and certified by the two signs of x − q·2π and x − (q+1)·2π;
equality, betweenness, sorting and midpoints reduce to such signs. Within a region where every
consulted sign is fixed, the computation asks the same questions, gets the same answers and
returns the same outcome (induction over the questions). Splitting a region into the feasible
signs of an undecided form therefore yields a finite decision tree whose leaves partition the
vertices into polyhedral cells with constant outcome.

**How it is computed.**
1. *Pieces.* A type's required congruences c·b ≡ 0 (mod 1), with Kawasaki's condition, cut out an
   affine subspace for each choice of the integers c·b = m on a basis of them. Each piece is that
   subspace clipped to the order constraints, kept only if it has the subspace's dimension. The
   dimension is 1 for every degree-4 type, and 2 for 10,120 of the 10,222 degree-6 types (1 for the other 102).
2. *Regions* are relatively open convex polygons (or segments, points), stored by the vertices of
   their closure, in exact integer arithmetic. An affine form f is > 0 on a region iff f ≥ 0 at all
   vertices and not 0 at all of them; = 0 iff 0 at all vertices; otherwise all three signs occur and
   the region is split by exact clipping. The witness point used for normalisations is the centroid.
3. *Symbolic execution.* The replay and the skeletons of src/rigid.py run on exact affine forms;
   an undecided sign raises a split, and the computation restarts in each part. Results computed
   in a region (replays, skeletons, normalisations) are reused in all its sub-regions. Only the
   geometry is symbolic: given the skeletons, the search over labelings, up sides and senses is
   combinatorial, so each leaf records both models (2: rigidly realizable; 1: only with reverse
   folds; 0: neither).

**Checks.** For every leaf a witness vertex is stored; at it the numeric computation must give the
leaf's value and the tree lookup must reach it. At every rotation of every corpus vertex shape, the
tree lookup must equal the numeric computation for every type. Result for all tables: 723,018 witnesses and
4,362,610 (vertex rotation, type) pairs, no disagreement.

**Cost** (one core, `nice -n 19`):

| table | types | seconds | leaves (cells) | sign queries |
|---|---|---|---|---|
| degree 4, k = 0 … 4 | 653 | 78.5 | 14,373 | 1.17 M |
| degree 6, k = 0 | 4 | 0.4 | 68 | 7 k |
| degree 6, k = 1 | 54 | 11.9 | 1,422 | 0.20 M |
| degree 6, k = 2 | 688 | 341.2 | 32,058 | 6.0 M |
| degree 6, k = 3 | 9,476 | 9,249.7 | 675,096 | 160.2 M |

In all, 10,875 types in 2.7 hours. The trees are small, 3.1 MB for all tables (2.75 MB for degree 6,
k = 3), because realizability is constant on most pieces and such trees collapse to one value per
piece: 6,922 of the 9,476 degree-6, k = 3 types, 392 of 688 at k = 2, 430 of 556 for degree 4,
k = 4. The slowest type took 8.0 s. The witnesses, kept only for the
checks, are about 60 times larger.

How the cost came down (all measured):
- one degree-4, k = 2 type: 145 s with a z3 oracle evaluating the witness through z3; 41 s reading
  the witness point once per region;
- all 12 degree-4, k = 2 types: 7.1 s with skeletons (each step's geometry once, the choices
  combinatorially; numerically the skeletons are 77 times faster than the reference on the
  corpus histories, with identical results); 0.4 s with the polygon oracle;
- 20 degree-6, k = 2 types: 275.7 s with three z3 calls per new form; 132.2 s with one; 91.0 s
  with results reused within a region (the trees of these three versions are identical); 13.1 s
  with the polygon oracle, whose trees are rooted at the pieces and were checked with
  tools/verify_cells.py instead.

### The corpus with the physical filter from the cells

`results/lookup_first_build/corpus_cells_reverse.jsonl` and `corpus_cells_halfturn.jsonl`
(`--table=tables --rigid=reverse|halfturn --cells=cells`, 170 s per pattern). No history was
decided by search: every realizability came from a tree.

| | resolved | total time on the 23 resolved with and without the filter | Pygmy Jerboa (070) |
|---|---|---|---|
| tables, no filter | 29 | 107.2 s | 94.5 s |
| tables, realizable with reverse folds | 23 | 15.2 s | 4.6 s |
| tables, rigidly realizable (half-turns only) | 6 | | |

- On the 23 patterns resolved both ways, the cost is the same for 19 and higher with the filter
  for 4 (016: 18 vs 15, 070: 16 vs 13, 292: 7 vs 5, 303: 7 vs 6): the cheapest local histories
  there were not physically realizable.
- The filter loses 6 patterns. Dutch Rabbit, Squirrel, Kirby and kitsune have a degree-4 vertex with
  no realizable history within the tables; Color-change Pinwheel and Marnie keep local options at
  every vertex, but no consistent global resolution.
- The time is the global solve: removing non-physical histories shrinks the problem, most for the
  hard patterns (070: 94.5 s to 4.6 s, 016: 4.9 s to 2.5 s, 302: 3.9 s to 1.6 s). The trees add
  0.09 s of loading per pattern (0.14 s instead of 0.05 s), which is why the median rises
  slightly (0.30 s instead of 0.20 s).
- With half-turns only, 6 patterns resolve: reverse folds are needed almost everywhere a
  cancellation is.

## The patterns lost to the physical filter

Of the 29 patterns the tables resolve, 6 are lost with `--rigid=reverse`:
- **Local:** Dutch Rabbit, Squirrel, Kirby and kitsune each have vertices of one of two shapes,
  (22.5°, 90°, 157.5°, 90°) and (67.5°, 90°, 112.5°, 90°) (sectors α, 90°, 180° − α, 90°), with no
  realizable history within the tables (degree 4, at most 4 hidden rays, every step anchored).
- **Global:** Color-change Pinwheel and Marnie keep local options at every vertex but end in an
  unsatisfiable core (8 and 10 vertices). They need the same feature, pinned scaffold lines:
  `docs/report_pinwheel_marnie.md`.

**The first missing feature: scaffold lines.** `tools/grid_search.py` searches breadth-first over
layered states of one vertex, with every fold line on the 22.5° grid, steps as in flatvertex (any
even set of preimages of a line; visible creases may also be erased), and all labelings, up sides
and senses. Both shapes are reached in 5 steps:

| shape (rays in units of 22.5°) | history: rays of each step (made / erased) |
|---|---|
| {0, 1, 5, 12} | {2, 10} / –; {5, 15} / – (reverse fold); {12} / {2}; {0} / {10} (reverse fold); {1} / {15} |
| {0, 3, 7, 12} | {1, 9} / –; {6, 12} / – (reverse fold); {3} / {9}; {7} / {1} (reverse fold); {0} / {6} |

Both are realizable with reverse folds and not rigidly. Each begins with a step that makes only
hidden rays: a fold along a line through the vertex that carries no final crease and is erased
later. The tables exclude such steps, because they required every step to contain a visible ray or
an erased hidden ray (so that every hidden direction is fixed when it is made). The line is not
free: moving it by ±0.001 rad, with every later hidden ray recomputed from its anchor, breaks (L)
at step 3, where one of its rays is erased together with a visible crease. Its direction is fixed,
but by a later step.

## Scaffold lines in the tables

A *scaffold step* makes only hidden rays. `lookup.py build n k DIR --scaffold` enumerates the
types with at least one scaffold step into `DIR/deg{n}_k{k}_s.json`; the anchored types stay in
`deg{n}_k{k}.json`, unchanged (rebuilding them with the new code reproduces the committed files).

**Enumeration.** As for anchored types, with two changes:
- A scaffold step is allowed. The direction (mod π) of its first ray becomes a new variable λ_j, and
  its other rays are filled in from it by Theorem 11.
- Visible ray 0 is in the first step that contains a visible ray (the rotation normalisation; the
  first step may now contain none).

Every condition (L) and every flatness condition is then a congruence c·a + d·λ ≡ 0 (mod π) with
integer c, d.

**Elimination (`eliminate`).** Unimodular integer row operations leave equivalent congruences.
For each line j in turn, Euclid's algorithm on column j leaves at most one row with a nonzero entry
there; that pivot row is set aside. Pivot row j has no entries in the columns before j. If some
column has no pivot, the lines are not pinned. That line can move continuously, with the earlier
lines following it through their pivot rows, without violating any congruence, and the strict
inequalities leave a continuum of valid lines. Such types are counted (`unpinned`) and dropped, so there are still no free hidden
directions. Otherwise the pivot rows are triangular in the lines, and each reads
d_jj λ_j ≡ r_j(a, λ_{j+1}, …) (mod π), which has exactly |d_jj| solutions mod π for every
right-hand side. Hence the vertices for which some lines satisfy all the congruences are exactly
those satisfying the rows left without lines. These rows are the type's index congruences, as for
anchored types. A *branch* is a choice of u_j ∈ [0, |d_jj|) for each line:
λ_j = (r_j + u_j π)/d_jj, where r_j substitutes the later lines unreduced. These are affine in a,
so a branch is a continuous choice of lines; reducing a line mod π before substituting it would
make the labels jump where the line crosses 0 ≡ π. Rays use λ_j mod π.

**Runtime.** For a scaffold type whose index congruences hold, every branch is replayed; the replay
places a scaffold step's first ray at λ_j or λ_j + π, whichever lies between its neighbours.

**Exact cells.** Per (piece, branch), as for anchored types. The lines are affine in b with
denominators dividing the product of the pivots, so the forms of a scaffold table have numerators
over 4 · lcm(products of pivots), stored as `den` in the cells meta (8 to 240 here).

| Degree | k | Types | Unpinned (dropped) | Build | Cells: leaves | Cells: build | den |
|---|---|---|---|---|---|---|---|
| 2 | 1 / 2 | 0 / 1 | 0 / 2 | < 1 s | 0 / 2 | < 1 s | 4 / 8 |
| 4 | 1 / 2 / 3 / 4 | 0 / 12 / 156 / 2754 | 0 / 12 / 66 / 1619 | 0.0 / 1.3 / 41.9 / 1582.0 s | 0 / 90 / 1860 / 60984 | 0.0 / 0.4 / 8.7 / 294.9 s | 4 / 8 / 24 / 240 |
| 6 | 1 / 2 | 0 / 282 | 0 / 152 | 0.0 / 62.0 s | 0 / 12164 | 0.0 / 125.4 s | 4 / 8 |

Degree 6 with k = 3 is not built. The anchored degree-6, k = 3 table took 4528 s and its cells
9250 s; the scaffold version is expected to cost several times more.

**Checks.**
- `tools/verify_cells.py` on every scaffold cells file: every leaf witness, and every (corpus vertex
  rotation, type, branch), agree with the numeric computation. Largest: degree 4, k = 4, with 60984
  witnesses and 1537088 corpus checks; 0 disagree.
- `tools/validate_scaffold.py tables` over the corpus (63 vertex shapes, 257 shape × depth
  comparisons) compares with the solver's own local search, whose scaffold lines are bisector
  images.
  - Every scaffold history of the search whose lines are pinned (`tools/pinning.py`) is in the
    tables: 2259 of them, against 5400 in the tables.
  - Every scaffold history of the tables is pinned by the numeric measure, independently of the
    elimination: 0 exceptions.
- At the two shapes of Dutch Rabbit, Squirrel, Kirby and kitsune, the tables give 2 histories each
  realizable with reverse folds; one at each is the 5-step history of the grid search above.
- `tools/verify_examples.py`, with and without `--table=tables`, reproduces the paper's numbers.

### The corpus with scaffold lines

`results/scaffold_lines/corpus_cells_reverse_dl.jsonl` and `corpus_cells_reverse.jsonl`
(`--table=tables --rigid=reverse --cells=cells`, with and without `--dl`, 170 s per pattern). Every
realizability came from a tree; none was decided by search.

| | resolved | timeouts | total time on the 28 resolved both ways |
|---|---|---|---|
| tables without the filter (reference) | 29 | 1 | |
| physical filter, no scaffold lines | 23 | 1 | |
| physical filter, scaffold lines | 28 | 3 (kitsune, Marnie, 304) | 86.0 s |
| physical filter, scaffold lines, `--dl` | 30 | 1 (304) | 18.2 s |

With `--dl`, all six patterns lost to the physical filter resolve:

| pattern | hidden rays (lower bound) | time | unfiltered tables |
|---|---|---|---|
| Dutch Rabbit | 19 (16) | 1.7 s | 19 |
| Squirrel | 24 (22) | 1.5 s | 22 |
| Kirby | 17 (16) | 1.8 s | 18 |
| kitsune | 24 (17) | 52.9 s | 17 |
| Color-change Pinwheel | 16 (12) | 0.2 s | 17 |
| Marnie | 37 (10) | 15.6 s | 19 |

None of these six is certified minimal: each cost is a minimum over the options the search gave,
and Marnie's varies between runs.
- **054 Hexagon** resolves (8 hidden rays). Without scaffold lines it was a local failure even with
  no filter.
- **Physical resolutions can be cheaper.** 016 needs 13 hidden rays (18 with the filter and no
  scaffold lines, 15 with no filter) and 070 needs 12 (16 and 13).
- **They can also be dearer**, where the unfiltered resolution used histories that are not
  realizable: Squirrel, kitsune and Marnie.
- **Time.** The 23 patterns resolved before take 12.7 s in total, against 15.2 s before scaffold
  lines. Lookups are computed per vertex and per depth on demand, so the extra types cost little.
  The median is 0.2 s.
- **`--dl`** is what lets kitsune and Marnie finish: under the default arithmetic both time out at
  150 s.
