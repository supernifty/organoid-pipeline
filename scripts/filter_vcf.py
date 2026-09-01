#!/usr/bin/env python3
import argparse
import gzip


def open_text(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode, encoding="utf-8")


def parse_info(info_field):
    info = {}
    if info_field in ("", "."):
        return info
    for entry in info_field.split(";"):
        if "=" in entry:
            key, value = entry.split("=", 1)
            info[key] = value
        else:
            info[entry] = True
    return info


def parse_number(value):
    if value in ("", "."):
        raise ValueError("missing numeric value")
    return float(value.split(",")[0])


def sample_depth(format_map):
    try:
        return parse_number(format_map.get("DP", "."))
    except ValueError:
        ad = format_map.get("AD", ".")
        if ad in ("", "."):
            raise ValueError("missing sample depth")
        values = [parse_number(value) for value in ad.split(",") if value not in ("", ".")]
        if not values:
            raise ValueError("missing sample depth")
        return sum(values)


def main():
    parser = argparse.ArgumentParser(description="Filter VCF by AF and depth.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sample")
    parser.add_argument("--tumour-sample")
    parser.add_argument("--normal-sample")
    parser.add_argument("--af", type=float, required=True)
    parser.add_argument("--dp", type=float)
    parser.add_argument("--tumour-dp", type=float)
    parser.add_argument("--normal-dp", type=float)
    parser.add_argument("--dp-field", default="DP")
    parser.add_argument("--info-af", action="store_true")
    parser.add_argument("--pass-only", action="store_true")
    args = parser.parse_args()

    if args.tumour_sample or args.normal_sample or args.tumour_dp is not None or args.normal_dp is not None:
        if not all(
            value is not None
            for value in (args.tumour_sample, args.normal_sample, args.tumour_dp, args.normal_dp)
        ):
            raise ValueError(
                "Sample-aware mode requires --tumour-sample, --normal-sample, "
                "--tumour-dp, and --normal-dp"
            )
        sample_aware = True
    else:
        if args.sample is None or args.dp is None:
            raise ValueError("Legacy mode requires --sample and --dp")
        sample_aware = False

    sample_index = None
    tumour_index = None
    normal_index = None

    with open_text(args.input, "r") as fin, open(args.output, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                if line.startswith("#CHROM"):
                    columns = line.rstrip("\n").split("\t")
                    samples = columns[9:]
                    if sample_aware:
                        for sample in (args.tumour_sample, args.normal_sample):
                            if sample not in samples:
                                raise ValueError(f"Sample {sample} not found in {args.input}")
                        tumour_index = samples.index(args.tumour_sample)
                        normal_index = samples.index(args.normal_sample)
                    else:
                        if args.sample not in samples:
                            raise ValueError(f"Sample {args.sample} not found in {args.input}")
                        sample_index = samples.index(args.sample)
                continue

            if (sample_aware and (tumour_index is None or normal_index is None)) or (
                not sample_aware and sample_index is None
            ):
                raise ValueError(f"Encountered records before header in {args.input}")

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            if args.pass_only and fields[6] not in ("PASS", "."):
                continue

            info = parse_info(fields[7])
            format_keys = fields[8].split(":")
            if sample_aware:
                tumour_map = dict(zip(format_keys, fields[9 + tumour_index].split(":")))
                normal_map = dict(zip(format_keys, fields[9 + normal_index].split(":")))

                try:
                    tumour_depth = sample_depth(tumour_map)
                    normal_depth = sample_depth(normal_map)
                except ValueError:
                    continue

                try:
                    if args.info_af:
                        af = parse_number(info.get("AF", "."))
                    else:
                        af = parse_number(tumour_map.get("AF", "."))
                except ValueError:
                    continue

                if tumour_depth < args.tumour_dp or normal_depth < args.normal_dp or af < args.af:
                    continue
            else:
                sample_values = fields[9 + sample_index].split(":")
                format_map = dict(zip(format_keys, sample_values))

                try:
                    dp = parse_number(info.get(args.dp_field, "."))
                except ValueError:
                    try:
                        dp = parse_number(format_map.get(args.dp_field, "."))
                    except ValueError:
                        continue

                try:
                    if args.info_af:
                        af = parse_number(info.get("AF", "."))
                    else:
                        af = parse_number(format_map.get("AF", "."))
                except ValueError:
                    continue

                if dp < args.dp or af < args.af:
                    continue

            fout.write(line)


if __name__ == "__main__":
    main()
