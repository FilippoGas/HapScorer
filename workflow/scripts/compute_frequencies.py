# workflow/scripts/compute_frequencies.py

"""
Computes the populational frequency of each haplotype.

This script processes the genotype matrix and the translated haplotype metadata.
It unpivots the genotype matrix, extracts individual alleles, and aggregates
their occurrences to compute the absolute and relative percentage frequency
of each unique transcript haplotype within the cohort.

Supports dual execution modes:
1. Native Snakemake workflow integration via the `script:` directive.
2. Standalone Command Line Interface (CLI) execution via `argparse` for manual debugging.

Expected Snakemake Inputs (or CLI Arguments):
    - **genotypes** (`--genotypes`): CSV file mapping sample genotypes to haplotype identifiers.
    - **haplotypes_aa** (`--haplotypes-aa`): CSV file containing the translated haplotypes metadata.

Expected Snakemake Outputs (or CLI Arguments):
    - **haplotypes_annotated** (`--out-haplotypes`): Annotated CSV file updated with the `freq_perc` metric.
"""

import argparse
import logging
import os
import sys
import warnings
from pathlib import Path
from typing import TYPE_CHECKING, Union

from utils import init_snakemake_log, setup_logger, setup_warnings

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake = Snakemake()

logger = setup_logger(__name__)
setup_warnings()


def parse_args() -> argparse.Namespace:
    """
    Parses command-line arguments for the frequency computation workflow.

    Returns:
        args (:obj:`argparse.Namespace`): Parsed execution parameters.
    """
    parser = argparse.ArgumentParser(description="Computes haplotype frequencies using Polars.")
    parser.add_argument("--genotypes", required=True, help="Path to input genotypes CSV.")
    parser.add_argument("--haplotypes-aa", required=True, help="Path to input haplotypes AA CSV.")
    parser.add_argument("--out-haplotypes", required=True, help="Path to output finalized annotated haplotypes CSV.")
    parser.add_argument("--threads", type=int, default=1, help="Maximum number of threads for Polars.")
    parser.add_argument("--verbose", action="store_true", help="Enable detailed DEBUG logging output.")
    return parser.parse_args()


def run_pipeline(
    genotypes_path: Union[str, Path],
    haplotypes_path: Union[str, Path],
    out_path: Union[str, Path],
    threads: int,
) -> None:
    """
    Executes the frequency computation pipeline using Polars.

    Args:
        genotypes_path (:obj:`str` or :obj:`pathlib.Path`): Path to the genotype matrix CSV.
        haplotypes_path (:obj:`str` or :obj:`pathlib.Path`): Path to the translated haplotypes CSV.
        out_path (:obj:`str` or :obj:`pathlib.Path`): Destination path for the annotated haplotypes CSV file.
        threads (:obj:`int`): Number of CPU worker threads to spawn for multiprocessing execution of `polars` python module.
    """
    genotypes_path = Path(genotypes_path)
    haplotypes_path = Path(haplotypes_path)
    out_path = Path(out_path)

    # Set the thread limit at the operating system level
    os.environ["POLARS_MAX_THREADS"] = str(threads)
    import polars as pl

    lf_geno = pl.scan_csv(genotypes_path)
    lf_haplo = pl.scan_csv(haplotypes_path)
    logger.debug(f"Scanning input LazyFrames: '{genotypes_path}' and '{haplotypes_path}'")

    # Extract sample columns (all except the transcript identifier)
    sample_cols = [c for c in lf_geno.collect_schema().names() if c != "ensembl_transcript_id"]
    n_samples = len(sample_cols)
    total_alleles = n_samples * 2

    if total_alleles == 0:
        raise ValueError("No samples found in the genotypes matrix.")

    logger.info(
        f"Initiating haplotype frequency computation across {n_samples} samples ({total_alleles} total alleles) "
        f"using Polars with {threads} worker thread(s)..."
    )

    # Unpivot the genotype matrix to a long format
    melted = lf_geno.unpivot(
        index="ensembl_transcript_id", on=sample_cols, variable_name="sample", value_name="genotype"
    )

    # Split the "hap1/hap2" format directly into two new columns
    melted = melted.with_columns(
        [
            pl.col("genotype").str.split_exact("/", 1).struct.field("field_0").alias("hap1"),
            pl.col("genotype").str.split_exact("/", 1).struct.field("field_1").alias("hap2"),
        ]
    )

    # Reconstruct the full haplotype_transcript_id for both alleles
    hap1_ids = melted.select(
        pl.concat_str([pl.col("ensembl_transcript_id"), pl.col("hap1")], separator="-").alias("haplotype_transcript_id")
    )
    hap2_ids = melted.select(
        pl.concat_str([pl.col("ensembl_transcript_id"), pl.col("hap2")], separator="-").alias("haplotype_transcript_id")
    )

    # Concatenate alleles and count frequencies
    all_ids = pl.concat([hap1_ids, hap2_ids])
    counts = all_ids.group_by("haplotype_transcript_id").len(name="freq_count")

    # Join counts back to the main haplotype lazyframe
    lf_haplo = lf_haplo.join(counts, on="haplotype_transcript_id", how="left")

    # Fill missing frequencies with 0 and compute percentage
    lf_haplo = lf_haplo.with_columns(
        (pl.col("freq_count").fill_null(0) / total_alleles * 100).round(3).alias("freq_perc")
    ).drop("freq_count")

    # Materialize the data by executing the entire computation plan
    df_final = lf_haplo.collect()

    zero_freq_df = df_final.filter(pl.col("freq_perc") == 0.0)
    zero_freq_count = len(zero_freq_df)
    if zero_freq_count > 0:
        zero_ids = zero_freq_df["haplotype_transcript_id"].to_list()

        # Detailed breakdown available only in DEBUG mode to prevent log pollution
        logger.debug(
            f"Complete list of {zero_freq_count} haplotype record(s) with 0.0% populational frequency: {', '.join(zero_ids)}"
        )

        # Dynamic warning formatting based on element threshold
        if zero_freq_count <= 10:
            ids_fmt = ", ".join(zero_ids)
        else:
            ids_fmt = f"{', '.join(zero_ids[:5])}, and {zero_freq_count - 5} more"

        warnings.warn(
            f"Detected {zero_freq_count} haplotype record(s) with 0.0% populational frequency "
            f"absent from sample genotype calls in '{genotypes_path}' ({ids_fmt}).",
            UserWarning,
        )

    logger.info(
        f"Materialized Polars computation graph. Computed populational frequencies for {len(df_final)} unique transcript haplotype records."
    )

    # Relocate 'freq_perc' before 'hgnc_symbol' to match R output formatting
    cols = df_final.columns
    cols.remove("freq_perc")
    if "hgnc_symbol" in cols:
        idx = cols.index("hgnc_symbol")
        cols.insert(idx, "freq_perc")
        df_final = df_final.select(cols)

    # Uppercase Boolean formatting for parity with R
    bool_cols = ["start_gain", "multiple_starts_in_5", "stop_loss_3_translated"]
    valid_cols = [c for c in bool_cols if c in df_final.columns]

    if valid_cols:
        df_final = df_final.with_columns([pl.col(c).cast(pl.String).str.to_uppercase() for c in valid_cols])

    # Export
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df_final.write_csv(out_path, null_value="NA")

    logger.info(f"Successfully exported annotated haplotype dataset with frequency metrics to '{out_path}'.")
    logger.info("Populational haplotype frequency computation workflow completed successfully.")


if __name__ == "__main__":
    is_snakemake = False

    if "snakemake" in globals():
        is_snakemake = True
        init_snakemake_log(snakemake, base_dir=Path.cwd())
    else:
        args = parse_args()

    genotypes_path = snakemake.input["genotypes"] if is_snakemake else args.genotypes
    haplotypes_path = snakemake.input["haplotypes"] if is_snakemake else args.haplotypes_aa
    out_path = snakemake.output["haplotypes_annotated"] if is_snakemake else args.out_haplotypes
    threads = snakemake.threads if is_snakemake else args.threads

    if not is_snakemake and args.verbose:
        logger.setLevel(logging.DEBUG)

    run_pipeline(
        genotypes_path=genotypes_path,
        haplotypes_path=haplotypes_path,
        out_path=out_path,
        threads=threads,
    )
