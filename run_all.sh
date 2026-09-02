#!/usr/bin/env bash
set -e
cd "$(dirname "$0")"
for f in theorems/T0*.py; do echo "== $f"; python3 "$f" 2>&1 | grep -v -e Plink -e warnings.warn; done
echo "== tests"; python3 -m pytest -q tests
