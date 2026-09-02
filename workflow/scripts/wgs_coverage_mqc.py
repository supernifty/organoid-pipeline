#!/usr/bin/env python3
"""Convert mosdepth region output to a MultiQC custom-content table."""

from __future__ import annotations

import argparse
import gzip
from collections import OrderedDict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--regions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--role", choices=("baseline", "organoid"), default="organoid")
    parser.add_argument("--expected-depth", type=float)
    parser.add_argument("--warning-fraction", type=float, default=0.5)
    parser.add_argument("--thresholds", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    totals: OrderedDict[str, list[float]] = OrderedDict()
    with gzip.open(args.regions, "rt") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"Malformed mosdepth regions row at line {line_number}")
            bases = int(fields[2]) - int(fields[1])
            total, count = totals.setdefault(fields[0], [0.0, 0])
            totals[fields[0]] = [total + float(fields[-1]) * bases, count + bases]
    if not totals:
        raise ValueError("mosdepth produced no WGS callable regions")
    threshold_totals = OrderedDict()
    if args.thresholds:
        with gzip.open(args.thresholds, "rt") as handle:
            threshold_columns = None
            data_rows = 0
            for line_number, line in enumerate(handle, 1):
                fields = line.rstrip("\n").split("\t")
                if (
                    line_number == 1
                    and len(fields) >= 4
                    and fields[0].lstrip("#").lower() == "chrom"
                    and fields[1:4] == ["start", "end", "region"]
                ):
                    threshold_columns = len(fields) - 4
                    continue
                if len(fields) < 5:
                    raise ValueError(
                        f"Malformed mosdepth thresholds row at line {line_number}: "
                        "expected at least 5 tab-separated fields"
                    )
                if threshold_columns is None:
                    threshold_columns = len(fields) - 4
                elif len(fields) - 4 != threshold_columns:
                    raise ValueError(
                        f"Malformed mosdepth thresholds row at line {line_number}: "
                        f"expected {threshold_columns} threshold columns"
                    )
                try:
                    region_bases = int(fields[2]) - int(fields[1])
                    values = [int(value) for value in fields[4:]]
                except ValueError as error:
                    raise ValueError(
                        f"Malformed mosdepth thresholds row at line {line_number}: "
                        "coordinates and threshold counts must be integers"
                    ) from error
                if region_bases < 0:
                    raise ValueError(
                        f"Malformed mosdepth thresholds row at line {line_number}: end precedes start"
                    )
                data_rows += 1
                for index, value in enumerate(values, 1):
                    covered, bases = threshold_totals.get(index, (0, 0))
                    threshold_totals[index] = (covered + value, bases + region_bases)
            if data_rows == 0:
                raise ValueError("mosdepth thresholds output contains no data rows")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as out:
        out.write("# id: 'wgs_per_contig_coverage'\n")
        out.write("# section_name: 'WGS callable coverage by contig'\n")
        out.write("# description: 'Mean depth across canonical callable territory.'\n")
        out.write("# format: 'tsv'\n")
        out.write("# plot_type: 'table'\n")
        threshold_names = ("pct_1x", "pct_3x", "pct_5x", "pct_10x", "pct_20x", "pct_50x")
        out.write(
            "sample\trole\tmean_autosomal_depth\texpected_depth\tdepth_status\t"
            + "\t".join(threshold_names[: len(threshold_totals)])
            + ("\t" if threshold_totals else "")
            + "\t".join(totals)
            + "\n"
        )
        means = [str(total / count) for total, count in totals.values()]
        autosomal = [
            value
            for contig, value in totals.items()
            if contig.removeprefix("chr").isdigit() and 1 <= int(contig.removeprefix("chr")) <= 22
        ]
        total_bases = sum(count for _, count in autosomal)
        mean_autosomal = sum(total for total, _ in autosomal) / total_bases if total_bases else 0.0
        status = "not_configured"
        if args.expected_depth is not None:
            status = "LOW" if mean_autosomal < args.expected_depth * args.warning_fraction else "OK"
        percentages = [
            str(100 * covered / bases if bases else 0.0)
            for covered, bases in threshold_totals.values()
        ]
        out.write(
            f"{args.sample}\t{args.role}\t{mean_autosomal}\t{args.expected_depth if args.expected_depth is not None else '.'}\t{status}\t"
            + "\t".join(percentages)
            + ("\t" if percentages else "")
            + "\t".join(means)
            + "\n"
        )


if __name__ == "__main__":
    main()
