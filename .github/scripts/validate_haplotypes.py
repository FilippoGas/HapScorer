# .github/scripts/validate_haplotypes.py

"""
Validates the structural integrity, schema correctness, and biological sanity of
generated haplotype and genotype datasets across evaluated chromosomes.
"""

import sys
from pathlib import Path

import pandas as pd

# Parse command-line argument for base directory (default: root directory)
base_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("")
chromosomes = ["chr3", "chr6", "chr8", "chr21", "chr22"]

print(f"=== INITIATING HAPLOTYPE & GENOTYPE VALIDATION ON: {base_dir} ===")

required_haplo_cols = [
    "ensembl_gene_id",
    "ensembl_transcript_id",
    "haplotype_transcript_id",
    "haplotype_gene_id",
    "variants",
    "rsid",
    "variant_types",
    "hgnc_symbol",
    "strand",
]

for chrom in chromosomes:
    haplo_file = base_dir / f"results/haplotype_sequences/haplotypes/{chrom}_haplotypes.csv"
    geno_file = base_dir / f"results/haplotype_sequences/genotypes/{chrom}_genotypes.csv"

    # File Integrity Verification (Existence and Non-Zero Size)
    assert haplo_file.exists(), f"Haplotype output missing for chromosome {chrom}"
    assert geno_file.exists(), f"Genotype output missing for chromosome {chrom}"

    assert haplo_file.stat().st_size > 100, f"Haplotype file for {chrom} is empty or corrupted"
    assert geno_file.stat().st_size > 100, f"Genotype file for {chrom} is empty or corrupted"

    # Data Loading
    df_haplo = pd.read_csv(haplo_file)
    df_geno = pd.read_csv(geno_file)

    # Column Schema Validation
    assert set(required_haplo_cols).issubset(df_haplo.columns), f"Invalid header format in {haplo_file.name}"
    assert "ensembl_transcript_id" in df_geno.columns, f"Missing transcript ID column in {geno_file.name}"
    assert len(df_geno.columns) > 1, f"No sample columns detected in {geno_file.name}"

    # Data Volume Verification
    assert len(df_haplo) > 0, f"No haplotype records found in {haplo_file.name}"
    assert len(df_geno) > 0, f"No genotype records found in {geno_file.name}"

    # Biological & Haplotype Integrity Verification
    # Wild-Type (.0) mandatory presence per transcript
    wt_haplos = df_haplo[df_haplo["haplotype_gene_id"].astype(str).str.endswith(".0")]
    assert len(wt_haplos["ensembl_transcript_id"].unique()) == len(
        df_haplo["ensembl_transcript_id"].unique()
    ), f"Missing wild-type (.0) entry for some transcripts in {haplo_file.name}"

    # Referencing Integrity: All genotype alleles must be declared in haplotypes.csv
    sample_cols = [c for c in df_geno.columns if c != "ensembl_transcript_id"]
    used_haplo_ids = set()
    for col in sample_cols:
        for genotype in df_geno[col].dropna():
            h1, h2 = str(genotype).split("/")
            used_haplo_ids.add(h1)
            used_haplo_ids.add(h2)

    defined_haplo_ids = set(df_haplo["haplotype_gene_id"])
    missing_ids = used_haplo_ids - defined_haplo_ids
    assert not missing_ids, f"Undeclared haplotype identifiers found in {geno_file.name}: {missing_ids}"

    # Transcript Alignment Verification
    assert set(df_geno["ensembl_transcript_id"]).issubset(
        set(df_haplo["ensembl_transcript_id"])
    ), f"Transcript ID mismatch between genotypes and haplotypes in {chrom}"

    print(f"[{chrom}] Successfully validated: {len(df_haplo)} haplotypes across {len(df_geno)} transcripts.")

print("=== ALL HAPLOTYPE AND GENOTYPE VALIDATION TESTS COMPLETED SUCCESSFULLY! ===")
