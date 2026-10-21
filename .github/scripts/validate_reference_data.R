# .github/scripts/validate_reference_data.R

library(tidyverse)

# Parse command-line argument for base directory (default: root directory)
args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) > 0) args[1] else ""

message(paste("=== INITIATING VALIDATION ON:", base_dir, "==="))

# Construct dynamic file paths
file_transcripts <- file.path(base_dir, "resources/reference_data/protein_coding_transcripts.csv")
file_cds         <- file.path(base_dir, "resources/reference_data/wt_cds.csv")
file_aa          <- file.path(base_dir, "resources/reference_data/wt_aa.csv")

# 1. File Integrity Verification (Existence and Non-Zero Size)
stopifnot("Transcript file missing" = file.exists(file_transcripts))
stopifnot("CDS file missing"        = file.exists(file_cds))
stopifnot("Amino acid file missing" = file.exists(file_aa))

stopifnot("Transcript file is empty or corrupted" = file.info(file_transcripts)$size > 1000)
stopifnot("CDS file is empty or corrupted"        = file.info(file_cds)$size > 1000)
stopifnot("Amino acid file is empty or corrupted" = file.info(file_aa)$size > 1000)

# 2. Data Loading
df_tx  <- read_csv(file_transcripts, show_col_types = FALSE)
df_cds <- read_csv(file_cds, show_col_types = FALSE)
df_aa  <- read_csv(file_aa, show_col_types = FALSE)

# 3. Column Structure Validation
cols_tx  <- c("ensembl_transcript_id", "hgnc_symbol", "chromosome_name", "transcript_biotype")
cols_seq <- c("ensembl_transcript_id")

stopifnot("Invalid header format in transcript metadata file" = all(cols_tx %in% colnames(df_tx)))
stopifnot("Invalid header format in CDS sequence file"        = all(cols_seq %in% colnames(df_cds)))
stopifnot("Invalid header format in amino acid sequence file" = all(cols_seq %in% colnames(df_aa)))

# 4. Data Volume Verification (expected human protein-coding transcript count > 10,000)
message(paste("Total transcripts retrieved:", nrow(df_tx)))
stopifnot("Insufficient transcript count retrieved" = nrow(df_tx) > 10000)

# 5. Biological Sanity Checks
# A) Positive control: verify presence of canonical marker gene (TP53)
has_tp53 <- "TP53" %in% df_tx$hgnc_symbol
stopifnot("Positive control failed: TP53 gene missing from dataset" = has_tp53)

# B) Accession ID validation: verify Ensembl transcript identifier syntax (ENST prefix)
valid_ids <- all(str_starts(df_tx$ensembl_transcript_id, "ENST"))
stopifnot("Invalid transcript identifiers detected (missing mandatory 'ENST' prefix)" = valid_ids)

# C) Nucleotide sequence integrity: verify IUPAC character compliance in CDS (A, C, G, T, N)
cds_col <- intersect(colnames(df_cds), c("coding", "cds", "sequence"))[1]
sample_cds <- head(df_cds[[cds_col]], 100)
valid_dna  <- all(str_detect(sample_cds, "^[ACGTNacgtn]+$"))
stopifnot("Non-IUPAC/invalid characters detected within CDS sequences" = valid_dna)

# D) Protein translation integrity: verify canonical start codon methionine ('M') frequency
aa_col <- intersect(colnames(df_aa), c("peptide", "protein_sequence", "sequence"))[1]
sample_aa <- head(df_aa[[aa_col]], 100)
start_with_m <- mean(str_starts(sample_aa, "M")) > 0.8
stopifnot("Translation start site check failed: less than 80% of protein sequences initiate with Methionine (M)" = start_with_m)

message("=== ALL VALIDATION TESTS COMPLETED SUCCESSFULLY! ===")
