#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter a mutational signature definition TSV to a configured subset."
    )
    parser.add_argument("--definition", required=True, help="Input signature definition TSV.")
    parser.add_argument("--output", required=True, help="Filtered output TSV.")
    parser.add_argument(
        "--signatures",
        required=True,
        nargs="+",
        help="Signature names to keep, e.g. SBS1 SBS5.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    wanted = args.signatures
    wanted_set = set(wanted)
    seen = {}

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(args.definition, "r", encoding="utf-8", newline="") as fin:
        reader = csv.reader(fin, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Definition file is empty: {args.definition}") from exc

        if not header or header[0] != "Sig":
            raise ValueError(
                f"Expected first column to be 'Sig' in definition file: {args.definition}"
            )

        selected_rows = []
        for row in reader:
            if not row:
                continue
            signature_name = row[0]
            if signature_name in wanted_set:
                selected_rows.append(row)
                seen[signature_name] = True

    missing = [name for name in wanted if name not in seen]
    if missing:
        missing_str = ", ".join(missing)
        raise ValueError(
            f"Configured signatures missing from {args.definition}: {missing_str}"
        )

    selected_rows.sort(key=lambda row: wanted.index(row[0]))

    with open(output_path, "w", encoding="utf-8", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(selected_rows)


if __name__ == "__main__":
    main()
