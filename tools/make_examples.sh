#!/bin/sh
# Run the full analysis on every pattern in examples/ and write its figures to out/examples/:
# NAME_folds.png, NAME_hasse.png, NAME_hasse.dot and NAME_folds.txt (see README.md).
# Usage: tools/make_examples.sh [extra solver flags, e.g. "--table=tables --rigid=reverse --cells=cells --dl"]
set -e
cd "$(dirname "$0")/.."
mkdir -p out/examples
for f in examples/*.fold; do
  name=$(basename "$f" .fold)
  echo "== $name"
  python3 src/fold_history.py --out="out/examples/$name" $1 "$f" | tail -1
done
ls out/examples
