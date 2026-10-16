#!/bin/sh
# Triage every pattern named in a list (default: tools/curated_corpus.txt) and append one JSON
# line per pattern to an output file.  Usage: tools/run_corpus.sh [list] [out.jsonl] [seconds]
# Patterns are looked up in data/flat-folder/examples/instagram (run fetch_flatfolder_corpus.sh).
cd "$(dirname "$0")/.."
LIST=${1:-tools/curated_corpus.txt}; OUT=${2:-results/corpus_run.jsonl}; T=${3:-170}
while read f; do
  [ -z "$f" ] && continue
  res=$(timeout "$T" python3 src/fold_history.py --triage --budget=$((T - 20)) "data/flat-folder/examples/instagram/$f" 2>/dev/null | tail -1)
  [ -z "$res" ] && res="{\"file\": \"$f\", \"status\": \"timeout\"}"
  echo "$res" >> "$OUT"; echo "$res" | cut -c1-160
done < "$LIST"
