# workflow/scripts/sequences_download.R

# Filippo Gastaldello - 29/04/24

# Use biomaRt to download annotations and sequences of all human protein coding
# transcripts (non MT) together with transcript annotations

# Redirect standard output and messages to Snakemake log file
log_file <- file(snakemake@log[[1]], open = "at")
sink(log_file)
sink(log_file, type = "message")

library(tidyverse)
library(biomaRt)

options(timeout = 300)
httr::set_config(httr::config(connecttimeout = 60, timeout = 300))

# Retrieve execution parameters from Snakemake params
p_biomart     <- snakemake@params[["biomart"]]
p_dataset     <- snakemake@params[["dataset"]]
p_host        <- snakemake@params[["host"]]
p_biotype     <- snakemake@params[["biotype"]]
p_filters     <- snakemake@params[["filters"]]
p_batch_size  <- as.numeric(snakemake@params[["batch_size"]])
p_sleep_delay <- as.numeric(snakemake@params[["sleep_delay"]])

cat(sprintf("Connecting to Ensembl BioMart...\n  Host: %s\n  BioMart: %s\n  Dataset: %s\n",
            p_host, p_biomart, p_dataset))
ensembl <- useMart(biomart = p_biomart,
                   dataset = p_dataset,
                   host    = p_host)
cat("Successfully established connection to Ensembl BioMart.\n\n")

# get ensembl ID of all protein coding transcripts
protein_coding_transcripts_featurepage <- getBM(attributes = c("ensembl_transcript_id",
                                                               "hgnc_symbol",
                                                               "transcript_tsl",
                                                               "chromosome_name",
                                                               "start_position",
                                                               "end_position",
                                                               "strand",
                                                               "transcript_biotype",
                                                               "uniprot_isoform"),
                                                filters = "transcript_biotype",
                                                values = p_biotype,
                                                mart = ensembl,
                                                useCache=FALSE)

# remove:
# transcript with tsl>3
# non protein transcripts (retained introns and non coding exons only composed by UTRs)
# transcripts non aligned to chromosomes
# mitochondrial transcripts
# Apply declarative filters safely without dynamic code execution
protein_coding_transcripts <- protein_coding_transcripts_featurepage

if (!is.null(p_filters) && is.list(p_filters) && length(p_filters) > 0) {

  # 1. Filter by allowed Transcript Support Levels (TSL)
  if (!is.null(p_filters$allowed_tsl) && length(p_filters$allowed_tsl) > 0) {
    allowed_tsl_vec <- unlist(p_filters$allowed_tsl)
    cat(sprintf("Filtering transcripts by allowed TSL: %s\n", paste(allowed_tsl_vec, collapse = ", ")))
    protein_coding_transcripts <- protein_coding_transcripts %>%
      filter(transcript_tsl %in% allowed_tsl_vec)
  }

  # 2. Exclude specified chromosomes
  if (!is.null(p_filters$excluded_chromosomes) && length(p_filters$excluded_chromosomes) > 0) {
    excluded_chrom_vec <- unlist(p_filters$excluded_chromosomes)
    cat(sprintf("Excluding chromosomes: %s\n", paste(excluded_chrom_vec, collapse = ", ")))
    protein_coding_transcripts <- protein_coding_transcripts %>%
      filter(!(chromosome_name %in% excluded_chrom_vec))
  }

  # 3. Filter chromosome names using regular expression matching
  if (!is.null(p_filters$chromosome_regex) && nchar(trimws(p_filters$chromosome_regex)) > 0) {
    chrom_regex <- as.character(p_filters$chromosome_regex)
    cat(sprintf("Filtering chromosome names with regex pattern: '%s'\n\n", chrom_regex))
    protein_coding_transcripts <- protein_coding_transcripts %>%
      filter(str_detect(chromosome_name, regex(chrom_regex)))
  }

} else {
  cat("No filtering rules provided. Retaining all retrieved transcripts.\n")
}

# get coordinate information for CDS of the previously retrieved transcripts
protein_coding_transcripts_structurepage <- getBM(attributes = c("ensembl_gene_id",
                                                                 "ensembl_transcript_id",
                                                                 "ensembl_peptide_id",
                                                                 "transcript_start",
                                                                 "transcript_end",
                                                                 "transcript_length",
                                                                 "5_utr_start",
                                                                 "5_utr_end",
                                                                 "3_utr_start",
                                                                 "3_utr_end",
                                                                 "exon_chrom_start",
                                                                 "exon_chrom_end"),
                                                  filters = "ensembl_transcript_id",
                                                  values = protein_coding_transcripts$ensembl_transcript_id,
                                                  mart = ensembl)

# merge annotation together
protein_coding_transcripts <- protein_coding_transcripts_structurepage %>% left_join(protein_coding_transcripts_featurepage, by = "ensembl_transcript_id")


# prepare dataframe to store sequences both with, and without UTRs
sequences_cds <- data.frame(coding = character(),
                        ensembl_transcript_id = character())
sequences_3utr <- data.frame('3utr' = character(),
                        ensembl_transcript_id = character())
sequences_5utr <- data.frame('5utr' = character(),
                        ensembl_transcript_id = character())
sequences_aa <- data.frame(peptide = character(),
                        ensembl_transcript_id = character())
# get sequence for every transcript (cycle at batches of size 'p_batch_size' transcripts to avoid time out)
start <- 1; end <- p_batch_size; done <- FALSE

# time tracking
start_time <- Sys.time()

while (!done) {

  # check if end index exceeded dataframe length
  if (end >= length(unique(protein_coding_transcripts$ensembl_transcript_id))) {
    end <- length(unique(protein_coding_transcripts$ensembl_transcript_id))
    done <- TRUE
  }

  # concat new sequences to existing dataframe
  sequences_cds <- bind_rows(sequences_cds, getSequence(id = unique(protein_coding_transcripts$ensembl_transcript_id)[start:end],
                                                    type = "ensembl_transcript_id",
                                                    seqType = "coding",
                                                    mart = ensembl))
  sequences_5utr <- bind_rows(sequences_5utr, getSequence(id = unique(protein_coding_transcripts$ensembl_transcript_id)[start:end],
                                                    type = "ensembl_transcript_id",
                                                    seqType = "5utr",
                                                    mart = ensembl))
  sequences_3utr <- bind_rows(sequences_3utr, getSequence(id = unique(protein_coding_transcripts$ensembl_transcript_id)[start:end],
                                                     type = "ensembl_transcript_id",
                                                     seqType = "3utr",
                                                     mart = ensembl))
  sequences_aa <- bind_rows(sequences_aa, getSequence(id = unique(protein_coding_transcripts$ensembl_transcript_id)[start:end],
                                                type = "ensembl_transcript_id",
                                                seqType = "peptide",
                                                mart = ensembl))

  # print progress
  step_time <- Sys.time()
  duration <- difftime(step_time, start_time, units = "sec")
  cat(paste("Downloaded", end, "of", length(unique(protein_coding_transcripts$ensembl_transcript_id)), "sequences in", round(seconds_to_period(duration), 2)," \n"))

  start <- end+1; end <- end+p_batch_size
  Sys.sleep(p_sleep_delay)
}

# print total elapsed time
end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "sec")
cat(paste("Total elapsed time:", round(seconds_to_period(duration), 2), "\n"))

sequences_3utr$X3utr <- NULL
sequences_5utr$X5utr <- NULL

# Merge 5utr cds and 3utr
sequences <- sequences_5utr %>% left_join(sequences_cds, by = 'ensembl_transcript_id')
sequences <- sequences %>% left_join(sequences_3utr, by = 'ensembl_transcript_id')

# Remove transcripts with unavailable sequences
sequences <- sequences %>% filter(!`5utr`=="Sequence unavailable")
sequences <- sequences %>% filter(!`3utr`=="Sequence unavailable")
sequences <- sequences %>% filter(!coding=="Sequence unavailable")

# merge transcript annotations with sequences
sequences_aa <- sequences_aa %>% dplyr::rename(sequence = peptide)


# Discard transcripts without sequence
sequences_aa <- sequences_aa %>% filter(ensembl_transcript_id %in% sequences$ensembl_transcript_id)
protein_coding_transcripts <- protein_coding_transcripts %>% filter(ensembl_transcript_id %in% sequences$ensembl_transcript_id)

# save sequences and annotations to files to path specified in snakemake "download_wt_sequences" rule.
write_csv(sequences, file = snakemake@output[["sequences"]])
write_csv(sequences_aa, file = snakemake@output[["sequences_aa"]])
write_csv(protein_coding_transcripts, file = snakemake@output[["annotations"]])

# Close logging sinks
sink(type = "message")
sink()
