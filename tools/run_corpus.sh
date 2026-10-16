#!/bin/sh
# Triage every pattern named in a list (default: tools/curated_corpus.txt) and write one JSON line
# per pattern, in list order, to an output file.  Each pattern runs single-threaded; JOBS patterns
# run at a time.  Usage: tools/run_corpus.sh [list] [out.jsonl] [seconds] [jobs]
# Patterns are looked up in data/flat-folder/examples/instagram (run fetch_flatfolder_corpus.sh).
cd "$(dirname "$0")/.."
LIST=${1:-tools/curated_corpus.txt}; OUT=${2:-results/corpus_run.jsonl}; T=${3:-170}; JOBS=${4:-1}
TMP=$(mktemp -d)
# xargs appends "NNN file" to the fixed arguments: $1 = seconds, $2 = tmp dir, $3 = line number, $4 = file
grep -v '^$' "$LIST" | nl -nrz -w3 | xargs -P "$JOBS" -L 1 sh -c '
  res=$(timeout "$1" python3 src/fold_history.py --triage --budget=$(($1 - 20)) "data/flat-folder/examples/instagram/$4" 2>/dev/null | tail -1)
  [ -z "$res" ] && res="{\"file\": \"$4\", \"status\": \"timeout\"}"
  echo "$res" > "$2/$3"; echo "$res" | cut -c1-160' sh "$T" "$TMP"
cat "$TMP"/* > "$OUT"; rm -r "$TMP"
