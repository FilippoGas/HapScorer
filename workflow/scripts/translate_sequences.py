# workflow/scripts/translate_sequences.py

"""
Translates mutated coding nucleotide sequences (CDS) into amino acid sequences.

This script processes intermediate datasets containing haplotype-specific
nucleotide sequences ('nn_sequence'). It translates these sequences utilizing
standard codon tables and subsequently truncates the resulting amino acid
sequences at the first encountered premature or canonical stop codon, mirroring
the legacy R pipeline's behavior.

Supports dual execution modes:
1. Native Snakemake workflow integration via the `script:` directive.
2. Standalone Command Line Interface (CLI) execution via `argparse` for manual debugging.

Expected Snakemake Inputs (or CLI Arguments):
    - **mutated_cds_csv** (`--mutated-cds`): Path to the dataset containing mutated nucleotide sequences.

Expected Snakemake Outputs (or CLI Arguments):
    - **mutated_aa_csv** (`--out-aa`): Path for the generated dataset with translated amino acid sequences.

Optional Parameters:
    - **threads** (`--threads`): Number of CPU cores to utilize for multiprocessing (default: 1).
"""

import argparse
import multiprocessing as mp
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Union

import pandas as pd
from Bio.Data import CodonTable
from Bio.Seq import translate

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake = Snakemake()

# =============================================================================
# CONSTANTS
# =============================================================================

# Canonical and non-canonical initiation codons recognized under NCBI Standard Genetic Code (Table 1)
START_CODONS = set(CodonTable.unambiguous_dna_by_id[1].start_codons)

# =============================================================================
# HELPER FUNCTIONS & VARIANT PARSING
# =============================================================================


def parse_args() -> argparse.Namespace:
    """
    Parses command-line arguments for the sequence translation workflow.

    Returns:
        args (:obj:`argparse.Namespace`): An object containing the parsed execution parameters.
    """
    parser = argparse.ArgumentParser(description="Translates nucleotide sequences to amino acids.")
    parser.add_argument(
        "--mutated-cds",
        required=True,
        help="Path to input mutated nucleotide sequences CSV.",
    )
    parser.add_argument(
        "--out-aa",
        required=True,
        help="Output path for the translated amino acid sequences CSV.",
    )
    parser.add_argument("--threads", type=int, default=1, help="Number of CPU threads to utilize.")
    return parser.parse_args()


def translate_sequence(nn_seq: Union[str, float]) -> str:
    """
    Translates a single nucleotide sequence string into an amino acid sequence.

    Uses Biopython's module-level translate to natively return a string.
    Truncates the translation product at the first encountered stop codon via `to_stop=True`.
    Safely handles invalid, missing (NA), or excessively short sequences by
    returning an empty string. Trims sequences to a strict multiple of 3.

    Args:
        nn_seq (:obj:`str` or :obj:`float`): The input nucleotide sequence string, or float(NaN) if missing.

    Returns:
        translated_seq (:obj:`str`): The truncated amino acid sequence, or an empty string if the input is invalid.
    """
    if pd.isna(nn_seq) or not isinstance(nn_seq, str) or len(nn_seq) < 3:
        return ""

    # Ensure the sequence length is a strict multiple of 3
    valid_length = len(nn_seq) - (len(nn_seq) % 3)
    trimmed_seq = nn_seq[:valid_length]

    try:
        # Returns a string directly and stops at the first in-frame stop codon
        aa_seq = translate(trimmed_seq, to_stop=True)

        # Enforce translation of initial canonical and alternative start codons to Methionine ('M')
        # to maintain functional parity with R Biostrings default behavior
        if aa_seq and trimmed_seq[:3] in START_CODONS:
            aa_seq = "M" + aa_seq[1:]

        return aa_seq
    except Exception:
        # Fallback for unexpected nucleotide anomalies
        return ""


def run_pipeline(
    mutated_cds_path: Union[str, Path],
    out_aa_path: Union[str, Path],
    threads: int,
) -> None:
    """
    Executes the comprehensive sequence translation pipeline.

    Loads the haplotype sequence dataset, parallelizes the translation process
    across designated worker threads using `multiprocessing.Pool`, appends the
    resulting amino acid sequences, and exports the finalized CSV dataset.

    Args:
        mutated_cds_path (:obj:`str` or :obj:`pathlib.Path`): Path to the input mutated nucleotide sequences CSV file.
        out_aa_path (:obj:`str` or :obj:`pathlib.Path`): Output file path for the translated amino acid sequences CSV file.
        threads (:obj:`int`): Number of CPU worker threads to spawn for multiprocessing execution.
    """
    mutated_cds_path = Path(mutated_cds_path)
    out_aa_path = Path(out_aa_path)

    df_haplo = pd.read_csv(mutated_cds_path)

    if "nn_sequence" not in df_haplo.columns:
        raise ValueError("The input CSV dataset lacks the required 'nn_sequence' column.")

    # Utilize multiprocessing pool for scalable and efficient sequence translation
    with mp.Pool(threads) as pool:
        aa_sequences = pool.map(translate_sequence, df_haplo["nn_sequence"])

    # Append the derived translations back into the primary dataframe
    df_haplo["aa_sequence"] = aa_sequences

    # Uppercase Boolean formatting for parity with R
    bool_cols = ["start_gain", "multiple_starts_in_5", "stop_loss_3_translated"]
    for col in bool_cols:
        df_haplo[col] = df_haplo[col].map({True: "TRUE", False: "FALSE"})

    # Ensure the output directory hierarchy exists prior to export
    out_aa_path.parent.mkdir(parents=True, exist_ok=True)
    df_haplo.to_csv(out_aa_path, index=False, na_rep="NA")


# =============================================================================
# MAIN LOGIC
# =============================================================================

if __name__ == "__main__":
    is_snakemake = False

    if "snakemake" in globals():
        is_snakemake = True
        if hasattr(snakemake, "log") and snakemake.log:
            log_f = open(snakemake.log[0], "w", buffering=1)
            sys.stdout = log_f
            sys.stderr = log_f
    else:
        args = parse_args()

    mutated_cds_path = snakemake.input["mutated_cds_csv"] if is_snakemake else args.mutated_cds
    out_aa_path = snakemake.output["mutated_aa_csv"] if is_snakemake else args.out_aa
    threads = snakemake.threads if is_snakemake else args.threads

    run_pipeline(
        mutated_cds_path=mutated_cds_path,
        out_aa_path=out_aa_path,
        threads=threads,
    )
