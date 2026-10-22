# Results

All runs use the 48 patterns of `tools/curated_corpus.txt` (fetched by
`tools/fetch_flatfolder_corpus.sh`), 170 s per pattern on one core, one JSON line per pattern as
printed by `src/fold_history.py --triage`. Fields: `status` (`ok`, `local fail`, `global unsat`,
`timeout`), `hidden` (hidden rays), `cost`, `lower_bound`, `certified_minimal`, `t` (seconds),
and for table runs `local_s` and `table` (lookup statistics).

| file | run |
|---|---|
| `corpus_results.jsonl`, `corpus_results.md` | the corpus table of `paper/reading_history.tex` (computed with an earlier version of the local search; see `docs/ROADMAP.md`, §3) |
| `lookup/search.jsonl` | local search, current code |
| `lookup/tables_anchored.jsonl` | `--table=tables`, anchored types only |
| `lookup/tables_dl.jsonl` | `--table=tables --dl`, all types |
| `lookup/rigid_anchored.jsonl` | `--table=tables --rigid=halfturn --cells=cells`, anchored types only |
| `lookup/physical_anchored.jsonl` | `--table=tables --rigid=reverse --cells=cells`, anchored types only |
| `lookup/physical.jsonl` | `--table=tables --rigid=reverse --cells=cells` |
| `lookup/physical_dl.jsonl` | `--table=tables --rigid=reverse --cells=cells --dl` |
| `lookup/survey_rigid.jsonl`, `lookup/survey_reverse.jsonl` | `tools/rigid_survey.py`: how many of the anchored tables' histories are realizable, per vertex shape |

"Anchored types only" means the tables before scaffold types were added (`tables/*_s.json`
absent). `python3 tools/lookup_summary.py results/lookup/*.jsonl` summarises the runs, and
`python3 tools/compare_runs.py A.jsonl B.jsonl` compares two runs pattern by pattern.
