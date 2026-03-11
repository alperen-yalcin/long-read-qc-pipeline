#!/usr/bin/env python3
"""
plot_qc_metrics.py
──────────────────
Reads the per-read QC CSV produced by calc_fastq_stats.py and:
  1. Generates distribution plots for GC content, read length, and mean quality.
  2. Prints summary statistics (mean, median, std, min, max) to stdout.
  3. Saves a combined figure to a PNG file.

Usage:
    python plot_qc_metrics.py --input qc_results.csv --output qc_distributions.png
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for headless servers

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


# ── Aesthetic defaults ───────────────────────────────────────────────────
sns.set_theme(style="whitegrid", context="talk", palette="muted")
FIGSIZE = (18, 5)
DPI = 200


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Visualise per-read QC metrics from a CSV file."
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Path to input CSV (output of calc_fastq_stats.py)",
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Path to output PNG figure",
    )
    return parser.parse_args()


def print_summary(df: pd.DataFrame) -> None:
    """Print key summary statistics for the three QC metrics."""
    metrics = {
        "GC Content (%)": "gc_content",
        "Read Length (bp)": "read_length",
        "Mean Quality (Phred)": "mean_quality",
    }

    print("\n" + "=" * 60)
    print("  SUMMARY STATISTICS")
    print("=" * 60)

    for label, col in metrics.items():
        series = df[col]
        print(f"\n  ▸ {label}")
        print(f"      Mean   : {series.mean():.2f}")
        print(f"      Median : {series.median():.2f}")
        print(f"      Std    : {series.std():.2f}")
        print(f"      Min    : {series.min():.2f}")
        print(f"      Max    : {series.max():.2f}")

    print("\n" + "=" * 60 + "\n")


def plot_distributions(df: pd.DataFrame, output_path: str) -> None:
    """Create a 1×3 figure with histograms + KDE for each metric."""
    fig, axes = plt.subplots(1, 3, figsize=FIGSIZE)
    fig.suptitle("Long-Read QC — Metric Distributions", fontsize=16, fontweight="bold", y=1.02)

    # ── 1. GC Content ────────────────────────────────────────────────────
    ax = axes[0]
    sns.histplot(df["gc_content"], bins=50, kde=True, color="#4c72b0", ax=ax)
    ax.set_xlabel("GC Content (%)")
    ax.set_ylabel("Number of Reads")
    ax.set_title("GC Content Distribution")
    ax.axvline(df["gc_content"].median(), color="red", ls="--", lw=1.2, label=f'Median {df["gc_content"].median():.1f}%')
    ax.legend(fontsize=9)

    # ── 2. Read Length ───────────────────────────────────────────────────
    ax = axes[1]
    sns.histplot(df["read_length"], bins=50, kde=True, color="#55a868", ax=ax)
    ax.set_xlabel("Read Length (bp)")
    ax.set_ylabel("Number of Reads")
    ax.set_title("Read Length Distribution")
    ax.axvline(df["read_length"].median(), color="red", ls="--", lw=1.2, label=f'Median {df["read_length"].median():,.0f} bp')
    ax.legend(fontsize=9)

    # ── 3. Mean Quality Score ────────────────────────────────────────────
    ax = axes[2]
    sns.histplot(df["mean_quality"], bins=50, kde=True, color="#c44e52", ax=ax)
    ax.set_xlabel("Mean Phred Quality Score")
    ax.set_ylabel("Number of Reads")
    ax.set_title("Mean Quality Distribution")
    ax.axvline(df["mean_quality"].median(), color="blue", ls="--", lw=1.2, label=f'Median Q{df["mean_quality"].median():.1f}')
    ax.legend(fontsize=9)

    # ── Save ─────────────────────────────────────────────────────────────
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"✓ Figure saved to {output_path}")


def main():
    args = parse_args()
    df = pd.read_csv(args.input)

    print(f"Loaded {len(df):,} reads from {args.input}")

    print_summary(df)
    plot_distributions(df, args.output)


if __name__ == "__main__":
    main()
