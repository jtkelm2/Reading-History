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
