#!/usr/bin/env python3
"""
calc_fastq_stats.py
───────────────────
Reads a FASTQ file and computes per-read statistics:
  • GC content (%)
  • Read length (bp)
  • Mean Phred quality score

Output is a CSV with columns:
  read_id, gc_content, read_length, mean_quality

Usage:
    python calc_fastq_stats.py --input reads.fastq.gz --output qc_results.csv
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate per-read QC statistics from a FASTQ file."
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Path to input FASTQ file (.fastq or .fastq.gz)",
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Path to output CSV file",
    )
    return parser.parse_args()


def open_fastq(path: str):
    """Return a file handle, transparently handling gzip."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def compute_stats(record) -> dict:
    """Compute GC%, length, and mean quality for a single SeqRecord."""
    seq = record.seq
    qualities = record.letter_annotations["phred_quality"]

    gc = gc_fraction(seq) * 100.0
    length = len(seq)
    mean_qual = float(np.mean(qualities)) if qualities else 0.0

    return {
        "read_id": record.id,
        "gc_content": round(gc, 2),
        "read_length": length,
        "mean_quality": round(mean_qual, 2),
    }


def main():
    args = parse_args()
    input_path = args.input
    output_path = args.output

    # Ensure parent directory exists
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    fieldnames = ["read_id", "gc_content", "read_length", "mean_quality"]
    read_count = 0

    with open_fastq(input_path) as fh_in, open(output_path, "w", newline="") as fh_out:
        writer = csv.DictWriter(fh_out, fieldnames=fieldnames)
        writer.writeheader()

        for record in SeqIO.parse(fh_in, "fastq"):
            stats = compute_stats(record)
            writer.writerow(stats)
            read_count += 1

            if read_count % 50_000 == 0:
                print(f"  Processed {read_count:,} reads …", file=sys.stderr)

    print(f"✓ Finished — {read_count:,} reads written to {output_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
