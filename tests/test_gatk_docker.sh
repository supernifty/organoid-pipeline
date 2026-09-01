#!/bin/bash
set -euo pipefail

echo "Checking GATK Docker CLI..."

if ! command -v docker >/dev/null 2>&1; then
    echo "WARNING: docker not installed - skipping GATK Docker CLI check"
    exit 0
fi

if ! docker info >/dev/null 2>&1; then
    echo "WARNING: docker daemon not accessible - skipping GATK Docker CLI check"
    exit 0
fi

help_output="$(docker run --rm --platform linux/amd64 broadinstitute/gatk:4.4.0.0 gatk Mutect2 --help 2>&1)"

printf '%s\n' "$help_output" | grep -q -- "--germline-resource" || {
    echo "ERROR: Mutect2 help did not include --germline-resource"
    exit 1
}

printf '%s\n' "$help_output" | grep -q -- "--panel-of-normals" || {
    echo "ERROR: Mutect2 help did not include --panel-of-normals"
    exit 1
}

printf '%s\n' "$help_output" | grep -q -- "--interval-padding" || {
    echo "ERROR: Mutect2 help did not include --interval-padding"
    exit 1
}

printf '%s\n' "$help_output" | grep -q -- "--max-mnp-distance" || {
    echo "ERROR: Mutect2 help did not include --max-mnp-distance"
    exit 1
}

check_gatk_flag() {
    local tool="$1"
    local flag="$2"
    local output
    output="$(docker run --rm --platform linux/amd64 broadinstitute/gatk:4.4.0.0 gatk "$tool" --help 2>&1)"
    printf '%s\n' "$output" | grep -q -- "$flag" || {
        echo "ERROR: $tool help did not include $flag"
        exit 1
    }
}

check_gatk_flag HaplotypeCaller --emit-ref-confidence
check_gatk_flag HaplotypeCaller --interval-padding
check_gatk_flag GenotypeGVCFs --include-non-variant-sites
check_gatk_flag SelectVariants --exclude-non-variants
check_gatk_flag SelectVariants --select-type-to-include
check_gatk_flag VariantFiltration --filter-expression
check_gatk_flag VariantFiltration --filter-name

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST_DIR="$PIPELINE_DIR/tmp/codex/gatk-germline-filter"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
trap 'rm -rf "$TEST_DIR"' EXIT

cat > "$TEST_DIR/input.vcf" <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##INFO=<ID=QD,Number=1,Type=Float,Description="Quality by depth">
##INFO=<ID=SOR,Number=1,Type=Float,Description="Strand odds ratio">
##INFO=<ID=FS,Number=1,Type=Float,Description="Fisher strand">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Mapping quality rank sum">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Read position rank sum">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	9	.	A	AT	30	PASS	QD=2;FS=200;ReadPosRankSum=-20
1	10	.	A	C	30	PASS	QD=2;SOR=3;FS=60;MQ=40;MQRankSum=-12.5;ReadPosRankSum=-8
1	11	.	A	C	30	PASS	QD=1.9;SOR=3;FS=60;MQ=40;MQRankSum=-12.5;ReadPosRankSum=-8
1	12	.	A	C	29	PASS	QD=2;SOR=3;FS=60;MQ=40;MQRankSum=-12.5;ReadPosRankSum=-8
1	13	.	A	C	30	PASS	QD=2;SOR=3.1;FS=60;MQ=40;MQRankSum=-12.5;ReadPosRankSum=-8
1	14	.	A	C	30	PASS	QD=2;SOR=3;FS=60.1;MQ=40;MQRankSum=-12.5;ReadPosRankSum=-8
1	15	.	A	C	30	PASS	QD=2;SOR=3;FS=60;MQ=39.9;MQRankSum=-12.5;ReadPosRankSum=-8
1	16	.	A	C	30	PASS	QD=2;SOR=3;FS=60;MQ=40;MQRankSum=-12.6;ReadPosRankSum=-8
1	17	.	A	C	30	PASS	QD=2;SOR=3;FS=60;MQ=40;MQRankSum=-12.5;ReadPosRankSum=-8.1
1	20	.	A	AT	30	PASS	QD=2;FS=200;ReadPosRankSum=-20
1	21	.	A	AT	30	PASS	QD=1.9;FS=200;ReadPosRankSum=-20
1	22	.	A	AT	29	PASS	QD=2;FS=200;ReadPosRankSum=-20
1	23	.	A	AT	30	PASS	QD=2;FS=200.1;ReadPosRankSum=-20
1	24	.	A	AT	30	PASS	QD=2;FS=200;ReadPosRankSum=-20.1
1	25	.	A	C,AT	30	PASS	QD=2;FS=201;ReadPosRankSum=-20
1	30	.	A	<NON_REF>	.	PASS	.
1	31	.	A	.	.	PASS	.
EOF

run_gatk() {
    docker run --rm --platform linux/amd64 -v "$TEST_DIR:/data" -w /data \
        broadinstitute/gatk:4.4.0.0 gatk "$@"
}

run_gatk SelectVariants -V /data/input.vcf --exclude-non-variants true \
    --select-type-to-include SNP -O /data/snps.vcf.gz
run_gatk SelectVariants -V /data/input.vcf --exclude-non-variants true \
    --select-type-to-include INDEL --select-type-to-include MIXED -O /data/indels.vcf.gz
run_gatk VariantFiltration -V /data/snps.vcf.gz \
    --filter-name SNP_QD --filter-expression 'QD < 2.0' \
    --filter-name SNP_QUAL --filter-expression 'QUAL < 30.0' \
    --filter-name SNP_SOR --filter-expression 'SOR > 3.0' \
    --filter-name SNP_FS --filter-expression 'FS > 60.0' \
    --filter-name SNP_MQ --filter-expression 'MQ < 40.0' \
    --filter-name SNP_MQRankSum --filter-expression 'MQRankSum < -12.5' \
    --filter-name SNP_ReadPosRankSum --filter-expression 'ReadPosRankSum < -8.0' \
    -O /data/snps.filtered.vcf.gz
run_gatk VariantFiltration -V /data/indels.vcf.gz \
    --filter-name INDEL_QD --filter-expression 'QD < 2.0' \
    --filter-name INDEL_QUAL --filter-expression 'QUAL < 30.0' \
    --filter-name INDEL_FS --filter-expression 'FS > 200.0' \
    --filter-name INDEL_ReadPosRankSum --filter-expression 'ReadPosRankSum < -20.0' \
    -O /data/indels.filtered.vcf.gz
run_gatk MergeVcfs -I /data/snps.filtered.vcf.gz -I /data/indels.filtered.vcf.gz \
    -O /data/filtered.vcf.gz

test -s "$TEST_DIR/filtered.vcf.gz.tbi"
observed="$(gzip -dc "$TEST_DIR/filtered.vcf.gz" | awk '!/^#/ {print $2 "\t" $7}')"
expected="$(cat <<'EOF'
9	PASS
10	PASS
11	SNP_QD
12	SNP_QUAL
13	SNP_SOR
14	SNP_FS
15	SNP_MQ
16	SNP_MQRankSum
17	SNP_ReadPosRankSum
20	PASS
21	INDEL_QD
22	INDEL_QUAL
23	INDEL_FS
24	INDEL_ReadPosRankSum
25	INDEL_FS
EOF
)"
if [ "$observed" != "$expected" ]; then
    echo "ERROR: unexpected germline hard-filter output" >&2
    diff -u <(printf '%s\n' "$expected") <(printf '%s\n' "$observed") || true
    exit 1
fi

echo "GATK Docker CLI and germline hard-filter checks passed."
