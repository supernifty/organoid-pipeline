#!/usr/bin/env bash
set -euo pipefail

REFERENCE="${1:-resources/reference/genome.fa}"
OUTPUT="${2:-resources/reference/regions.bed.gz}"
FAI="${REFERENCE}.fai"

if [ ! -f "$FAI" ]; then
    echo "ERROR: FASTA index not found: $FAI" >&2
    echo "Run: samtools faidx $REFERENCE" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT")"

awk 'BEGIN { OFS="\t" } { print $1, 0, $2 }' "$FAI" | bgzip -c > "$OUTPUT"
tabix -f -p bed "$OUTPUT"

echo "Wrote genome-wide regions:"
echo "  $OUTPUT"
echo "  $OUTPUT.tbi"
