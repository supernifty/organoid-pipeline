#!/usr/bin/env python3
"""
Migrate samples.yaml from old format to new format.

Old format:
  samples:
    S1: ['in/S1_R1.fastq.gz', 'in/S1_R2.fastq.gz']  # list
  tumours:
    T1: S1

New format:
  samples:
    S1:
      fastq_1: in/S1_R1.fastq.gz
      fastq_2: in/S1_R2.fastq.gz
  tumours:
    T1: S1
"""

import sys
import yaml
import argparse

def migrate_samples(old_samples):
    """Convert old list format to new dict format."""
    new_samples = {}
    for sample, fastqs in old_samples.items():
        if isinstance(fastqs, list) and len(fastqs) >= 2:
            new_samples[sample] = {
                "fastq_1": fastqs[0],
                "fastq_2": fastqs[1]
            }
        elif isinstance(fastqs, dict):
            # Already in new format
            new_samples[sample] = fastqs
        else:
            print(f"Warning: Unexpected format for sample {sample}", file=sys.stderr)
    return new_samples

def migrate(input_path, output_path):
    with open(input_path) as f:
        old_data = yaml.safe_load(f)

    new_data = {
        "samples": migrate_samples(old_data.get("samples", {}))
    }

    if "tumours" in old_data:
        new_data["tumours"] = old_data["tumours"]

    with open(output_path, "w") as f:
        yaml.dump(new_data, f, default_flow_style=False, sort_keys=False)

    print(f"Migrated {len(new_data['samples'])} samples to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Migrate samples.yaml format")
    parser.add_argument("input", help="Input samples.yaml (old format)")
    parser.add_argument("output", nargs="?", help="Output samples.yaml (new format)")
    args = parser.parse_args()

    output = args.output or "config/samples.yaml"
    migrate(args.input, output)