#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

pixi install --locked
pixi run test
pixi run lint
pixi run format-check
pixi run python tests/create_dryrun_fixture.py --build grch38 --input-mode fastq --output tmp/codex/dryrun-grch38
env TMPDIR=tmp/codex XDG_CACHE_HOME=tmp/codex/cache \
  pixi run python -m snakemake --dry-run --cores 1 --quiet all \
  --configfile tmp/codex/dryrun-grch38/config.yaml
pixi run python tests/create_dryrun_fixture.py --build grch37 --output tmp/codex/dryrun-grch37
env TMPDIR=tmp/codex XDG_CACHE_HOME=tmp/codex/cache \
  pixi run python -m snakemake --dry-run --cores 1 --quiet all \
  --configfile tmp/codex/dryrun-grch37/config.yaml

rm -rf tmp/codex/dryrun-grch37 tmp/codex/dryrun-grch38 tmp/codex/pytest
echo "Organoid pipeline unit checks and GRCh37/GRCh38 complete DAG dry runs passed."
