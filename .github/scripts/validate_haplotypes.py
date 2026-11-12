# .github/scripts/validate_haplotypes.py

"""
Validates the structural integrity, schema correctness, and biological sanity of
generated haplotype, genotype, and sequence datasets across evaluated chromosomes.
"""

import re
import sys
from pathlib import Path

import pandas as pd

# Parse command-line argument for base directory (default: current directory)
base_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("")

# Dynamically discover generated chromosome output files
haplo_dir = base_dir / "results/haplotype_sequences/haplotypes"
haplo_files = sorted(
    haplo_dir.glob("*_haplotypes_annotations.csv"),
    key=lambda p: [int(text) if text.isdigit() else text for text in re.split(r"(\d+)", p.name)],
)
chromosomes = [f.name.replace("_haplotypes_annotations.csv", "") for f in haplo_files]

assert len(chromosomes) > 0, f"No annotated haplotype files found in {haplo_dir}"

print(f"=== INITIATING HAPLOTYPE, GENOTYPE & SEQUENCE BIOLOGICAL VALIDATION ON: {base_dir} ===")

required_anno_cols = [
    "ensembl_gene_id",
    "ensembl_transcript_id",
    "haplotype_transcript_id",
    "haplotype_gene_id",
    "variants",
    "rsid",
    "variant_types",
    "hgnc_symbol",
    "strand",
    "start_gain",
    "multiple_starts_in_5",
    "stop_loss_3_translated",
    "freq_perc",
]

required_seq_cols = [
    "haplotype_transcript_id",
    "ensembl_transcript_id",
    "nn_sequence",
    "aa_sequence",
]

for chrom in chromosomes:
    ann_file = base_dir / f"results/haplotype_sequences/haplotypes/{chrom}_haplotypes_annotations.csv"
    seq_file = base_dir / f"results/haplotype_sequences/haplotypes/{chrom}_haplotypes_sequences.csv"
    geno_file = base_dir / f"results/haplotype_sequences/genotypes/{chrom}_genotypes.csv"

    # File Integrity Verification (Existence and Non-Zero Size)
    assert ann_file.exists(), f"Haplotype annotations output missing for chromosome {chrom}"
    assert seq_file.exists(), f"Haplotype sequences output missing for chromosome {chrom}"
    assert geno_file.exists(), f"Genotype output missing for chromosome {chrom}"

    assert ann_file.stat().st_size > 100, f"Annotations file for {chrom} is empty or corrupted"
    assert seq_file.stat().st_size > 100, f"Sequences file for {chrom} is empty or corrupted"
    assert geno_file.stat().st_size > 100, f"Genotype file for {chrom} is empty or corrupted"

    # Data Loading
    df_anno = pd.read_csv(ann_file)
    df_seq = pd.read_csv(seq_file)
    df_geno = pd.read_csv(geno_file)

    # Schema and Structure Validation
    assert set(required_anno_cols).issubset(df_anno.columns), f"Missing required columns in {ann_file.name}"
    assert set(required_seq_cols).issubset(df_seq.columns), f"Missing required columns in {seq_file.name}"
    assert "ensembl_transcript_id" in df_geno.columns, f"Missing transcript ID column in {geno_file.name}"
    assert len(df_geno.columns) > 1, f"No sample columns detected in {geno_file.name}"
    assert len(df_anno) > 0, f"No haplotype records found in {ann_file.name}"
    assert len(df_seq) > 0, f"No sequence records found in {seq_file.name}"
    assert len(df_geno) > 0, f"No genotype records found in {geno_file.name}"
    assert len(df_anno) == len(df_seq), f"Record count mismatch between {ann_file.name} and {seq_file.name}"

    # Biological & Haplotype Integrity Verification
    # Wild-Type Baseline Integrity (.0)
    wt_haplos = df_anno[df_anno["haplotype_gene_id"].astype(str).str.endswith(".0")]
    assert len(wt_haplos["ensembl_transcript_id"].unique()) == len(
        df_anno["ensembl_transcript_id"].unique()
    ), f"Missing wild-type (.0) entry for some transcripts in {ann_file.name}"

    # Verify Wild-Type Biological Flags
    for _, wt_row in wt_haplos.iterrows():
        assert (
            str(wt_row["variants"]).lower() == "wt"
        ), f"WT haplotype has non-WT variant string in {wt_row['haplotype_transcript_id']}"
        assert (
            str(wt_row["start_gain"]).upper() == "FALSE"
        ), f"WT haplotype marked with start_gain in {wt_row['haplotype_transcript_id']}"
        assert (
            str(wt_row["multiple_starts_in_5"]).upper() == "FALSE"
        ), f"WT haplotype marked with multiple_starts_in_5 in {wt_row['haplotype_transcript_id']}"
        assert (
            str(wt_row["stop_loss_3_translated"]).upper() == "FALSE"
        ), f"WT haplotype marked with stop_loss_3_translated in {wt_row['haplotype_transcript_id']}"

    # Cohort Populational Frequency Sanity (Sum per transcript ~ 100%)
    freq_sums = df_anno.groupby("ensembl_transcript_id")["freq_perc"].sum()
    invalid_freq_txs = freq_sums[abs(freq_sums - 100.0) > 0.1]
    assert (
        invalid_freq_txs.empty
    ), f"Transcript haplotype frequency sum mismatch (expected 100%): {invalid_freq_txs.to_dict()}"

    # Biological Sequence Integrity & Translation Consistency
    df_combined = pd.merge(
        df_anno,
        df_seq[["haplotype_transcript_id", "nn_sequence", "aa_sequence"]],
        on="haplotype_transcript_id",
    )

    wt_m_starts = 0
    total_wt_aa = 0

    for _, row in df_combined.iterrows():
        nn_seq = str(row["nn_sequence"]) if pd.notna(row["nn_sequence"]) else ""
        aa_seq = str(row["aa_sequence"]) if pd.notna(row["aa_sequence"]) else ""
        is_wt = str(row["haplotype_gene_id"]).endswith(".0")

        if nn_seq:
            # Check IUPAC DNA base compliance
            assert set(nn_seq.upper()).issubset(
                set("ACGTN")
            ), f"Invalid DNA characters detected in {row['haplotype_transcript_id']}: {nn_seq}"

        if aa_seq:
            # Check nucleotide-to-amino acid length proportionality
            assert (
                len(nn_seq) >= len(aa_seq) * 3
            ), f"Nucleotide sequence length is insufficient for translated amino acid sequence in {row['haplotype_transcript_id']}"

            # Track Methionine initiation rate for Wild-Type (.0) haplotypes
            if is_wt:
                total_wt_aa += 1
                if aa_seq.startswith("M"):
                    wt_m_starts += 1

    # Verify that > 90% of Wild-Type protein sequences initiate with Methionine (M)
    if total_wt_aa > 0:
        m_rate = wt_m_starts / total_wt_aa
        assert (
            m_rate >= 0.90
        ), f"Wild-type translation start check failed for {chrom}: only {m_rate:.2%} of transcripts initiate with Methionine (M) (expected >= 90%)"

    # Referencing Integrity: All genotype alleles must be declared in haplotypes dataset
    sample_cols = [c for c in df_geno.columns if c != "ensembl_transcript_id"]
    used_haplo_ids = set()
    for col in sample_cols:
        for genotype in df_geno[col].dropna():
            h1, h2 = str(genotype).split("/")
            used_haplo_ids.add(h1)
            used_haplo_ids.add(h2)

    defined_haplo_ids = set(df_anno["haplotype_gene_id"])
    missing_ids = used_haplo_ids - defined_haplo_ids
    assert not missing_ids, f"Undeclared haplotype identifiers found in {geno_file.name}: {missing_ids}"

    # Transcript Alignment Verification
    assert set(df_geno["ensembl_transcript_id"]).issubset(
        set(df_anno["ensembl_transcript_id"])
    ), f"Transcript ID mismatch between genotypes and haplotypes in {chrom}"

    print(
        f"[{chrom}] Successfully validated: {len(df_anno)} haplotypes across {len(df_geno)} transcripts (All biological invariants verified)."
    )

print("=== ALL HAPLOTYPE, GENOTYPE & SEQUENCE BIOLOGICAL VALIDATION TESTS COMPLETED SUCCESSFULLY! ===")
