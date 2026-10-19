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
| 6 | 0 / 1 / 2 | 4 / 54 / 688 | 0.1 / 4.5 / 144.5 s | yes |

Degree 2 has no types with hidden rays: its first step uses both visible rays, so no later step
has an anchor. All builds: z3 4.16, one core, `nice -n 19`.

## Validation

`tools/validate_table.py tables` compares, at every vertex shape of degree ≤ 6 in the corpus and the
paper's examples (74 shapes, each once up to rotation and reflection) and every depth built, the
tables' histories with those of the local search, restricted to histories with no hidden-only
step. Result: 265 of 265 comparisons equal. `tools/verify_examples.py --table=tables` reproduces the
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
