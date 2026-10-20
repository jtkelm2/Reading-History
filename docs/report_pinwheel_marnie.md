# Color-change Pinwheel and Marnie: what they need

Color-change Pinwheel and Marnie are the two patterns lost to the physical filter
(`--table=tables --rigid=reverse --cells=cells`) by a global conflict: they have local options at
every vertex, but no consistent global resolution.

## Summary

Both patterns need the same feature as Dutch Rabbit, Squirrel, Kirby and kitsune: **pinned scaffold
lines**. A scaffold line is a hidden-only step, a fold along a line that carries no final crease and
whose rays are erased later. It is pinned when its direction is fixed by the later steps. With
scaffold lines added, both patterns resolve with physical local histories everywhere. They need no
deeper histories than the tables already hold (degree 4: k ≤ 4, degree 6: k ≤ 3), and no unpinned
scaffold lines.

## Method

`tools/scaffold_witness.py` builds the global model from the physical table options (every depth
the tables cover). With `--scaffolds` it also adds the scaffold histories of the solver's own local
search that pass two checks:
- they are realizable with reverse folds;
- every scaffold line in them is pinned.

The local search takes scaffold lines from bisector images, so these are a subset of all pinned
scaffold lines.

`tools/pinning.py` measures pinning as follows:
- Give each scaffold step a line parameter.
- Express every later condition, (L) and flatness, as a function of the parameters. These
  functions are piecewise linear.
- free_dims = (number of scaffold lines) − (rank of the Jacobian).

The Jacobian is computed numerically. Every reported free direction is confirmed by a finite step
along it that leaves every condition satisfied.

All claims below are computed with z3, within the global model's own restrictions: cancelled
segments bend at most twice, and so on.

## Color-change Pinwheel (058)

- **Tables alone.** The model is unsatisfiable even with every option at every depth. The core is 8
  degree-4 vertices: 7, 8, 11, 12, 15, 17, 20, 24.
- **With pinned scaffold lines.** The solver's search finds 16 hidden rays, the same cost as the
  old search without the physical filter. Scaffold lines appear at exactly four vertices: 8, 15, 17
  and 24.
- **The shape.** All four vertices have the same shape, (22.5°, 45°, 157.5°, 135°) up to rotation.
  At vertex 8 (creases at 67.5°, 90°, 135°, 292.5°) the history is:
  1. fold the scaffold line 0°/180°, which bisects the 135° sector;
  2. fold the two creases bounding that sector, 67.5° and 292.5°, together;
  3. make 135°, erasing 0°;
  4. make 90°, erasing 180°.

  It is realizable with reverse folds, not rigidly.
- **How the line is pinned.** Step 2 pins the scaffold line. Folding the creases a and b together
  over it requires 2λ ≡ a + b (mod π), so λ is their bisector or its perpendicular.

## Marnie (315)

- **Tables alone.** The adjacent degree-4 vertices 10 and 12 cannot both choose. With every option
  of the whole pattern in the model and only these two required to choose, the model is
  unsatisfiable, with core {10, 12}. z3's default arithmetic and its difference-logic solver agree.
- **With pinned scaffold lines.**
  - The pair 10, 12 becomes satisfiable.
  - The whole pattern resolves. The solver's search finds 47 hidden rays, with scaffold lines at 8
    vertices: 10, 19, 20, 22, 23, 32, 34 and 38. Vertices 23 and 38 have degree 6. The result is
    not certified minimal; the lower bound is 10.
  - A second, independent witness has 62 hidden rays and scaffold lines at 11 vertices. Every
    assertion of the model was checked against it: `results/lookup_first_build/scaffold_witness_315.json`.
- **The cost of being physical.** The unfiltered tables resolve Marnie with 19 hidden rays, but 8
  of its local histories are not realizable even with reverse folds.

## Side findings

- **The search's scaffold histories are mostly not pinned.** At Pinwheel's vertex 8, 135 of the
  160 physical scaffold histories with k ≤ 4 have one free dimension. In those, the bisector rule
  chose a line that the history leaves free, so they are free hidden directions in disguise.
  Neither pattern needs them.
- **Solver arithmetic.** On the full Marnie problem, z3's default arithmetic returned unknown after
  360 s. Its dense difference-logic solver (`smt.arith.solver=3`) decided it. The constraints are
  difference constraints on integer times together with Boolean choices. The Bellman–Ford variant
  (`=1`) returned unknown on a satisfiable instance ("incomplete"), so only `=3` is a candidate.
- **Memory.** Models with thousands of options at once used enough memory to push about 6 GB of
  other processes to swap. All solver runs are now capped at 6 GB of address space.

## Recommendations

1. **Implement pinned scaffold lines in the tables (degree ≤ 6, current depths).** This one feature
   accounts for all six lost patterns.
2. **Keep excluding unpinned scaffold lines.** They are free hidden directions, and none of the six
   patterns needs them.
3. **Use the dense difference-logic solver** (`smt.arith.solver=3`) in the global solve, after a
   corpus regression check (same results, timings compared). Marnie is not practical without it.
4. **Expect higher costs.** Physical resolutions use more hidden rays than unfiltered ones: Marnie
   needs 47 against 19, and the minimum is unknown.

## Reproduce

```
python3 tools/scaffold_witness.py PATTERN --witness                  # tables alone
python3 tools/scaffold_witness.py PATTERN --scaffolds --solve        # the solver's search
python3 tools/scaffold_witness.py PATTERN --scaffolds --witness --out=W.json
python3 tools/scaffold_witness.py MARNIE --need=10,12 [--scaffolds]  # the pair
```
