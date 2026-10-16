# fold-history

Reading the folding history of a flat-folded crease pattern from local histories at its vertices.

Given a crease pattern in FOLD format, `src/fold_history.py` finds the creases that must have been
made together (*folds*), the order between folds that every consistent history forces, and creases
that must once have existed and been erased (*cancelled segments*). It draws the folds on the
pattern and the Hasse diagram of the forced order. The theory and worked examples are in
`paper/reading_history.pdf`.

Everything computed is a necessary condition for foldability by the moves allowed. The order of
layers is not computed, so a resolution is evidence for a folding history, not a certificate.

## Layout

```
src/fold_history.py      the solver and analysis (single file)
paper/                   reading_history.tex (self-contained) and its two figure PNGs; `make`
examples/                the paper's two worked patterns, the Gentle Dragon, a crossover gadget
tools/                   scripts that reproduce the paper's numbers and tables, and run the corpus
results/                 corpus results (48 patterns) and the dragon run log
docs/ROADMAP.md          open problems and next steps
legacy/                  earlier versions of the script and of the paper, kept for reference
```

## Setup

```
pip install -r requirements.txt     # z3-solver, networkx, matplotlib
# Graphviz (the `neato` binary) is needed for the Hasse diagram PNG.
# A LaTeX installation with TikZ is needed to build the paper.
```

## Usage

```
python3 src/fold_history.py examples/example2_body.fold
python3 src/fold_history.py --out=out/dragon examples/komatsu_gentle_dragon.fold   # about 2 minutes
python3 src/fold_history.py --triage --budget=150 pattern.fold                      # solve only, JSON line
```

Options: `--out=PREFIX` (default: the input's basename), `--budget=SECONDS` for the solve,
`--tol=RADIANS` to fix the direction tolerance (otherwise tried in the order 1e-6, 1e-4, 5e-4,
2e-3 until every vertex is flat), `--triage` to skip the analysis and the figures.

The full run writes `PREFIX_folds.png` (folds coloured by depth, cancelled segments dashed),
`PREFIX_hasse.png` and `.dot`, and `PREFIX_folds.txt` (the creases in each fold and the Hasse
edges), and prints the cancelled segments and a JSON summary.

## Reproducing the paper

```
python3 tools/verify_examples.py                                   # running and body patterns
python3 tools/enumerate_resolutions.py examples/komatsu_gentle_dragon.fold
python3 tools/paper_tables.py                                      # LaTeX fragments into build/
tools/fetch_flatfolder_corpus.sh && tools/run_corpus.sh            # the 48-pattern corpus
cd paper && make
```

`paper/reading_history.tex` is self-contained: the fragments from `tools/paper_tables.py` were
pasted into it, so regenerate and paste them if the results change.

## Status

On the curated corpus (`tools/curated_corpus.txt`, from the Flat-Folder collection), 42 of 48
patterns are resolved; see `results/corpus_results.md`. Unresolved within a few minutes on
one core: Butterfly Base, Fish Base, Dutch Rabbit, Sheep, Assigned Solvable (15k vertices),
Marnie. See `docs/ROADMAP.md`.

## Credits

The corpus and two example files come from Flat-Folder by Jason S. Ku
(https://github.com/origamimagiro/flat-folder, MIT licence); the designs belong to their authors.
