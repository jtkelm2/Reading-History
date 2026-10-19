#!/bin/sh
# Triage every pattern named in a list (default: tools/curated_corpus.txt) and write one JSON line
# per pattern, in list order, to an output file.  Each pattern runs single-threaded; JOBS patterns
# run at a time.  Usage: tools/run_corpus.sh [list] [out.jsonl] [seconds] [jobs] [extra solver flags]
# Patterns are looked up in data/flat-folder/examples/instagram (run fetch_flatfolder_corpus.sh).
cd "$(dirname "$0")/.."
LIST=${1:-tools/curated_corpus.txt}; OUT=${2:-results/corpus_run.jsonl}; T=${3:-170}; JOBS=${4:-1}; EXTRA=${5:-}
TMP=$(mktemp -d)
# xargs appends "NNN file" to the fixed arguments: $1 = seconds, $2 = tmp dir, $3 = extra flags,
# $4 = line number, $5 = file
grep -v '^$' "$LIST" | nl -nrz -w3 | xargs -P "$JOBS" -L 1 sh -c '
  res=$(timeout "$1" python3 src/fold_history.py --triage --budget=$(($1 - 20)) $3 "data/flat-folder/examples/instagram/$5" 2>/dev/null | tail -1)
  [ -z "$res" ] && res="{\"file\": \"$5\", \"status\": \"timeout\"}"
  echo "$res" > "$2/$4"; echo "$res" | cut -c1-160' sh "$T" "$TMP" "$EXTRA"
cat "$TMP"/* > "$OUT"; rm -r "$TMP"
