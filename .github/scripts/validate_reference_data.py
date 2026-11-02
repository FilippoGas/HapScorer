# .github/scripts/validate_reference_data.py

import sys
from pathlib import Path

import pandas as pd

# Parse command-line argument for base directory (default: root directory)
base_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("")

print(f"=== INITIATING VALIDATION ON: {base_dir} ===")

# Construct dynamic file paths
file_transcripts = base_dir / "resources/reference_data/protein_coding_transcripts.csv"
file_cds = base_dir / "resources/reference_data/wt_cds.csv"
file_aa = base_dir / "resources/reference_data/wt_aa.csv"

# File Integrity Verification (Existence and Non-Zero Size)
assert file_transcripts.exists(), "Transcript file missing"
assert file_cds.exists(), "CDS file missing"
assert file_aa.exists(), "Amino acid file missing"

assert file_transcripts.stat().st_size > 1000, "Transcript file is empty or corrupted"
assert file_cds.stat().st_size > 1000, "CDS file is empty or corrupted"
assert file_aa.stat().st_size > 1000, "Amino acid file is empty or corrupted"

# Data Loading
df_tx = pd.read_csv(file_transcripts)
df_cds = pd.read_csv(file_cds)
df_aa = pd.read_csv(file_aa)

# Column Structure Validation
cols_tx = ["ensembl_transcript_id", "hgnc_symbol", "chromosome_name", "transcript_biotype"]
cols_seq = ["ensembl_transcript_id"]

assert set(cols_tx).issubset(df_tx.columns), "Invalid header format in transcript metadata file"
assert set(cols_seq).issubset(df_cds.columns), "Invalid header format in CDS sequence file"
assert set(cols_seq).issubset(df_aa.columns), "Invalid header format in amino acid sequence file"

# Data Volume Verification (expected human protein-coding transcript count > 10,000)
print(f"Total transcripts retrieved: {len(df_tx)}")
assert len(df_tx) > 10000, "Insufficient transcript count retrieved"

# Biological Sanity Checks
# A) Positive control: verify presence of canonical marker gene (TP53)
has_tp53 = "TP53" in df_tx["hgnc_symbol"].values
assert has_tp53, "Positive control failed: TP53 gene missing from dataset"

# B) Accession ID validation: verify Ensembl transcript identifier syntax (ENST prefix)
valid_ids = df_tx["ensembl_transcript_id"].astype(str).str.startswith("ENST").all()
assert valid_ids, "Invalid transcript identifiers detected (missing mandatory 'ENST' prefix)"

# C) Nucleotide sequence integrity: verify IUPAC character compliance in CDS (A, C, G, T, N)
cds_cols = [c for c in ["coding", "cds", "sequence"] if c in df_cds.columns]
assert len(cds_cols), "No valid CDS sequence column found"
cds_col = cds_cols[0]
sample_cds = df_cds[cds_col].dropna().head(100)
valid_dna = sample_cds.astype(str).str.contains(r"^[ACGTNacgtn]+$", regex=True).all()
assert valid_dna, "Non-IUPAC/invalid characters detected within CDS sequences"

# D) Protein translation integrity: verify canonical start codon methionine ('M') frequency
aa_cols = [c for c in ["peptide", "protein_sequence", "sequence"] if c in df_aa.columns]
assert len(aa_cols), "No valid amino acid sequence column found"
aa_col = aa_cols[0]
sample_aa = df_aa[aa_col].dropna().head(100)
start_with_m = (sample_aa.astype(str).str.startswith("M")).mean() > 0.8
assert (
    start_with_m
), "Translation start site check failed: less than 80% of protein sequences initiate with Methionine (M)"

print("=== ALL VALIDATION TESTS COMPLETED SUCCESSFULLY! ===")
