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
