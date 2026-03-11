# ═══════════════════════════════════════════════════════════════════════
#  Snakefile — Long-Read Sequencing QC Pipeline
# ═══════════════════════════════════════════════════════════════════════
#  Runs NanoPlot for general QC, then a custom per-read stats script,
#  and finally generates distribution plots.
#
#  Usage:
#      snakemake --cores <N>
# ═══════════════════════════════════════════════════════════════════════

configfile: "config.yml"

SAMPLE = config["sample"]          # e.g. "data/sample"
FASTQ  = SAMPLE + ".fastq.gz"


# ── Target rule ──────────────────────────────────────────────────────
rule all:
    input:
        "results/nanoplot/NanoStats.txt",
        "results/qc_results.csv",
        "results/qc_distributions.png",


# ── Rule 1: NanoPlot QC ─────────────────────────────────────────────
rule nanoplot_qc:
    """Run NanoPlot to generate comprehensive long-read QC reports."""
    input:
        fastq=FASTQ,
    output:
        stats="results/nanoplot/NanoStats.txt",
    params:
        outdir="results/nanoplot",
    threads: 4
    shell:
        """
        NanoPlot \
            --fastq {input.fastq} \
            --outdir {params.outdir} \
            --threads {threads} \
            --loglevel WARNING
        """


# ── Rule 2: Custom per-read statistics ──────────────────────────────
rule calc_stats:
    """Compute GC%, length, and mean quality for every read."""
    input:
        fastq=FASTQ,
    output:
        csv="results/qc_results.csv",
    shell:
        """
        python scripts/calc_fastq_stats.py \
            --input {input.fastq} \
            --output {output.csv}
        """


# ── Rule 3: Visualisation ───────────────────────────────────────────
rule plot_metrics:
    """Generate distribution plots and print summary statistics."""
    input:
        csv="results/qc_results.csv",
    output:
        png="results/qc_distributions.png",
    shell:
        """
        python scripts/plot_qc_metrics.py \
            --input {input.csv} \
            --output {output.png}
        """
