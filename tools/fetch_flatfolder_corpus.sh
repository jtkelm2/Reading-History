#!/bin/sh
# Fetch the example collection distributed with Flat-Folder (MIT licence, J. S. Ku)
# into data/flat-folder/examples.  Needs git with network access.
set -e
cd "$(dirname "$0")/.."
mkdir -p data
if [ ! -d data/flat-folder ]; then
  git clone --depth 1 --filter=blob:none --sparse https://github.com/origamimagiro/flat-folder.git data/flat-folder
  (cd data/flat-folder && git sparse-checkout set examples)
fi
ls data/flat-folder/examples
