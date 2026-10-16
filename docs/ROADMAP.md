# Open problems and next steps

## 1. Top-down structural reading, for teaching (the main goal)
Experienced folders read a crease pattern from the top down: the design paradigm, the
earliest structural creases, which regions become which flaps, and where gadgets sit. Much
of this is inverse design. For uniaxial bases Lang's map from crease pattern to base to tree
(T = G(F(S))) is explicit, and the active paths are the primary structural creases. Steps:
- Classify the design paradigm: uniaxial/tree, box pleating, 22.5 degrees, tessellation.
- Recover the tree and the circle or river packing where the paradigm allows; this is
  classical geometry, with no satisfiability solving.
- Recognize gadgets (grafts, point splits, level shifters, colour changes, sinks) as local
  templates, or as deviations from the base structure.
- Use the forced order as a supporting layer, for the earliest creases and the cancelled
  precreases. Its top layers can probably be approximated cheaply as long lines that pass
  straight through every vertex they meet.
- Generators (TreeMaker, box-pleating tools) give training or test data with ground-truth
  structure, if a learned component is wanted.

## 2. Question 25 of the paper
Does every flat vertex have a local history with cancellations, using single-line steps only?
- A single vertex is the problem of folding a closed chain. The open-chain (segment) case is
  solved by crimps and end folds (Arkin et al.).
- Lead: let a hidden ray act as a temporary second switch point for each segment-style
  simple fold, then erase it. The staircase vertex of the paper has exactly this shape.
- Smallest test: the vertices of the Bern-Hayes crossover. Search with at most 4 hidden
  pairs found no history.
- Literature still to read: Akitaya, Mitani, Kanamori, Fukui, "Unfolding simple folds from
  crease patterns" (OSME 2014); Demaine-O'Rourke, chapters 12-13; Lang, *Origami Design
  Secrets*.

## 3. Search engineering
- **Fish Base:** shrink each unsat core to a minimal one before expanding. Successive cores
  had 4, 19, 6 and 16 vertices, and each round was slower than the last.
- **Butterfly Base, Marnie, Dutch Rabbit, Sheep:** time-outs not yet diagnosed. These are
  good tests of the bend limit and the pass-through approximation.
- **Cancelled segments without a bend limit:** trace the path implied by the current model
  (it reflects at earlier creases and passes later ones) and refute it if it does not close.
- **Exact pass-through:** add the passing segment's hidden rays to the history at the vertex.

## 4. From folds to moves
Two folds made by one move both map onto the move's fold line under the folding map at that
moment. This is a global test for merging folds (for example mirror pairs) into moves.

## 5. Physical validity and colour
Couple the method with Flat-Folder's layer-order solver, so that intermediate states can be
checked and mountain/valley labels predicted.

## 6. Validation against ground truth
Check that known folding sequences (crane, Yakko, kabuto, published designer diagrams) are
linear extensions of the forced order, and that the recovered cancelled segments match the
known precreases.

## 7. Related work to cite
- Akitaya (2013): the extended crease pattern (one edge per fold) and a per-vertex complexity
  of degree minus two. These overlap with our folds and our excess-line cost.
- Akitaya, Demaine, Ku, "Simple folding is really hard": their definition of a simple fold
  forbids erasing creases.
