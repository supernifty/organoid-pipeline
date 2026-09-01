# Optional SEQC2 benchmark

SEQC2 HCC1395/HCC1395BL is a technical caller benchmark, not biological validation for monoclonal organoids. The v1.2 high-confidence region spans about 2.48 Gb; raw WGS alignments are large. Nothing here downloads data during ordinary installation, tests, or CI.

1. Copy `config.example.yaml` to an ignored local configuration.
2. Choose specific HCC1395 and HCC1395BL WGS replicates from NCBI SRA project `SRP162370` or the NCBI SEQC2 repository.
3. Record exact URLs, sizes, and SHA-256 checksums for both alignments, truth, high-confidence territory, and reference. Preview a selected resource with `pixi run python scripts/fetch_verified.py --url URL --output benchmarks/seqc2/inputs/FILE --sha256 SHA256`; add `--execute` only after review. The command checks remote size, workspace capacity, and the final digest.
4. Measure source depth with mosdepth over the supplied high-confidence regions.
5. Run `scripts/downsample_alignment.py` with each configured seed. The script uses samtools template-name hashing, preserving read pairs, and records achieved rather than assumed depth.
6. Preview storage and every command with `pixi run python scripts/run_seqc2_benchmark.py --config benchmarks/seqc2/config.local.yaml`. Add `--execute` only after review. The runner generates isolated manifests and work directories for tumour 6× against baseline 6×, 15×, and 30× for every seed, then runs exact-allele metrics.
7. Restrict normalized calls and truth to `High-Confidence_Regions_v1.2.bed`, then run `workflow/scripts/benchmark_metrics.py` for exact-allele precision, recall, F1, and false positives per callable gigabase. Compare Mutect2, Strelka2, intersection, and stringent catalogs; aggregate seeds and SBS96 outputs.

Truth resources are published in the [NCBI SEQC2 somatic working-group repository](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/). The full benchmark remains explicit because selecting replicates and downloading hundreds of gigabytes requires operator review.
