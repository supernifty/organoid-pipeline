#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

pixi install --locked
pixi run test
pixi run lint
pixi run format-check
pixi run python tests/create_dryrun_fixture.py --output tmp/codex/dryrun
env TMPDIR=tmp/codex XDG_CACHE_HOME=tmp/codex/cache \
  .pixi/envs/default/bin/snakemake --dry-run --cores 1 --quiet all \
  --configfile tmp/codex/dryrun/config.yaml

rm -rf tmp/codex/dryrun tmp/codex/pytest
echo "Organoid pipeline unit checks and complete DAG dry run passed."
