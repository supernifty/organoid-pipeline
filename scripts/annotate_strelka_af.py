#!/usr/bin/env python3
import argparse
import gzip


def open_text(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode, encoding="utf-8")


def parse_number_list(value):
    if value in ("", "."):
        return []
    numbers = []
    for part in value.split(","):
        if part in ("", "."):
            continue
        numbers.append(float(part))
    return numbers


def upsert_info(info_field, key, value):
    if info_field in ("", "."):
        return f"{key}={value}"

    entries = info_field.split(";")
    replaced = False
    for index, entry in enumerate(entries):
        if entry.startswith(f"{key}="):
            entries[index] = f"{key}={value}"
            replaced = True
            break
    if not replaced:
        entries.append(f"{key}={value}")
    return ";".join(entries)


def calculate_snv_af(ref, alt, format_map):
    if "AD" in format_map and format_map["AD"] not in ("", "."):
        counts = parse_number_list(format_map["AD"])
        if len(counts) >= 2:
            ref_count = counts[0]
            alt_count = counts[1]
            total = ref_count + alt_count
            return 0.0 if total == 0 else alt_count / total

    ref_key = f"{ref}U"
    alt_key = f"{alt}U"
    ref_count = sum(parse_number_list(format_map.get(ref_key, ".")))
    alt_count = sum(parse_number_list(format_map.get(alt_key, ".")))
    total = ref_count + alt_count
    return 0.0 if total == 0 else alt_count / total


def calculate_indel_af(format_map):
    tar = sum(parse_number_list(format_map.get("TAR", ".")))
    tir = sum(parse_number_list(format_map.get("TIR", ".")))
    total = tar + tir
    return 0.0 if total == 0 else tir / total


def calculate_combined_depth(format_keys, sample_fields):
    if "DP" not in format_keys:
        return None
    dp_index = format_keys.index("DP")
    depth = 0.0
    found = False
    for sample_field in sample_fields:
        values = sample_field.split(":")
        if dp_index >= len(values):
            continue
        parsed = parse_number_list(values[dp_index])
        if not parsed:
            continue
        depth += parsed[0]
        found = True
    return depth if found else None


def main():
    parser = argparse.ArgumentParser(description="Annotate Strelka somatic VCF with INFO/AF.")
    parser.add_argument("--mode", choices=["snv", "indel"], required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    af_header_present = False
    dp_header_present = False
    sample_index = None

    with open_text(args.input, "r") as fin, open(args.output, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("##"):
                if line.startswith("##INFO=<ID=AF,"):
                    af_header_present = True
                if line.startswith("##INFO=<ID=DP,"):
                    dp_header_present = True
                fout.write(line)
                continue

            if line.startswith("#CHROM"):
                if not af_header_present:
                    fout.write('##INFO=<ID=AF,Number=1,Type=Float,Description="Calculated allele frequency">\n')
                if args.mode == "indel" and not dp_header_present:
                    fout.write('##INFO=<ID=DP,Number=1,Type=Float,Description="Combined sample depth">\n')
                fout.write(line)
                columns = line.rstrip("\n").split("\t")
                samples = columns[9:]
                if args.sample not in samples:
                    raise ValueError(f"Sample {args.sample} not found in {args.input}")
                sample_index = samples.index(args.sample)
                continue

            if sample_index is None:
                raise ValueError(f"Encountered records before header in {args.input}")

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                fout.write(line)
                continue

            ref = fields[3]
            alt = fields[4].split(",")[0]
            format_keys = fields[8].split(":")
            sample_values = fields[9 + sample_index].split(":")
            format_map = dict(zip(format_keys, sample_values))
            if args.mode == "indel":
                af = calculate_indel_af(format_map)
                depth = calculate_combined_depth(format_keys, fields[9:])
                if depth is not None:
                    fields[7] = upsert_info(fields[7], "DP", f"{depth:.6g}")
            else:
                af = calculate_snv_af(ref, alt, format_map)
            fields[7] = upsert_info(fields[7], "AF", f"{af:.6g}")
            fout.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    main()
