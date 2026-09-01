#!/usr/bin/env python3
import argparse
import gzip
import heapq
import re
from contextlib import ExitStack


CONTIG_RE = re.compile(r"^##contig=<ID=([^,>]+)")


def open_text(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode, encoding="utf-8")


def is_passable(filter_value, allowed_filters):
    if filter_value in ("", ".", "PASS"):
        return True
    return all(part in allowed_filters for part in filter_value.split(";"))


def read_header(handle, path):
    headers = []
    contigs = []
    for line in handle:
        if not line.startswith("#"):
            return headers, contigs, line
        headers.append(line)
        match = CONTIG_RE.match(line)
        if match:
            contigs.append(match.group(1))
    return headers, contigs, None


def data_lines(handle, first_line):
    if first_line is not None:
        yield first_line
    yield from handle


def parse_record(line, path, line_number, contig_rank):
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 7:
        raise ValueError(f"{path}:{line_number}: expected at least 7 VCF columns")
    contig = fields[0]
    if contig not in contig_rank:
        raise ValueError(
            f"{path}:{line_number}: contig {contig!r} is absent from the Mutect2 VCF dictionary"
        )
    try:
        position = int(fields[1])
    except ValueError as error:
        raise ValueError(f"{path}:{line_number}: invalid VCF position {fields[1]!r}") from error
    locus = (contig_rank[contig], position)
    allele = (contig, fields[1], fields[3], fields[4])
    return locus, allele, fields[6]


def grouped_records(handle, first_line, path, contig_rank):
    current_locus = None
    group = []
    for line_number, line in enumerate(data_lines(handle, first_line), start=1):
        locus, allele, filter_value = parse_record(line, path, line_number, contig_rank)
        if current_locus is not None and locus < current_locus:
            raise ValueError(f"{path}: VCF records are not coordinate sorted")
        if current_locus is not None and locus != current_locus:
            yield current_locus, group
            group = []
        current_locus = locus
        group.append((allele, filter_value, line))
    if current_locus is not None:
        yield current_locus, group


def grouped_strelka_records(path, contig_rank, allowed_filters):
    with open_text(path, "r") as handle:
        _, _, first_line = read_header(handle, path)
        for locus, records in grouped_records(handle, first_line, path, contig_rank):
            alleles = {}
            for allele, filter_value, _ in records:
                passable = is_passable(filter_value, allowed_filters)
                alleles[allele] = alleles.get(allele, False) or passable
            yield locus, alleles


def merged_strelka_groups(paths, contig_rank, allowed_filters):
    iterators = [iter(grouped_strelka_records(path, contig_rank, allowed_filters)) for path in paths]
    heap = []
    for index, iterator in enumerate(iterators):
        try:
            locus, alleles = next(iterator)
        except StopIteration:
            continue
        heapq.heappush(heap, (locus, index, alleles, iterator))

    while heap:
        locus = heap[0][0]
        merged = {}
        while heap and heap[0][0] == locus:
            _, index, alleles, iterator = heapq.heappop(heap)
            for allele, passable in alleles.items():
                merged[allele] = merged.get(allele, False) or passable
            try:
                next_locus, next_alleles = next(iterator)
            except StopIteration:
                continue
            heapq.heappush(heap, (next_locus, index, next_alleles, iterator))
        yield locus, merged


def main():
    parser = argparse.ArgumentParser(description="Intersect normalized Mutect2 and Strelka VCFs.")
    parser.add_argument("--mutect2-vcf", required=True)
    parser.add_argument("--strelka-vcf", nargs="+", required=True)
    parser.add_argument("--allowed-filters", nargs="*", default=[])
    parser.add_argument("--output-vcf", required=True)
    args = parser.parse_args()

    allowed_filters = set(args.allowed_filters)
    total_mutect2 = 0
    total_strelka = 0
    written = 0

    with ExitStack() as stack:
        mutect2_handle = stack.enter_context(open_text(args.mutect2_vcf, "r"))
        output_handle = stack.enter_context(open(args.output_vcf, "w", encoding="utf-8"))
        headers, contigs, first_line = read_header(mutect2_handle, args.mutect2_vcf)
        if not contigs:
            raise ValueError(f"{args.mutect2_vcf}: VCF header contains no ##contig dictionary")
        if len(contigs) != len(set(contigs)):
            raise ValueError(f"{args.mutect2_vcf}: VCF header contains duplicate contig IDs")
        contig_rank = {contig: index for index, contig in enumerate(contigs)}
        output_handle.writelines(headers)

        strelka_groups = iter(
            merged_strelka_groups(args.strelka_vcf, contig_rank, allowed_filters)
        )
        try:
            strelka_locus, strelka_alleles = next(strelka_groups)
        except StopIteration:
            strelka_locus, strelka_alleles = None, {}

        for mutect2_locus, records in grouped_records(
            mutect2_handle, first_line, args.mutect2_vcf, contig_rank
        ):
            while strelka_locus is not None and strelka_locus < mutect2_locus:
                total_strelka += len(strelka_alleles)
                try:
                    strelka_locus, strelka_alleles = next(strelka_groups)
                except StopIteration:
                    strelka_locus, strelka_alleles = None, {}

            matches = strelka_alleles if strelka_locus == mutect2_locus else {}
            for allele, filter_value, line in records:
                total_mutect2 += 1
                if matches.get(allele, False) and is_passable(filter_value, allowed_filters):
                    output_handle.write(line)
                    written += 1

        while strelka_locus is not None:
            total_strelka += len(strelka_alleles)
            try:
                strelka_locus, strelka_alleles = next(strelka_groups)
            except StopIteration:
                strelka_locus = None

    print(
        f"Intersected {written} variants from {total_mutect2} Mutect2 records "
        f"against {total_strelka} Strelka alleles"
    )


if __name__ == "__main__":
    main()
