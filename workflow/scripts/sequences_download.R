# workflow/scripts/sequences_download.R

# Filippo Gastaldello - 29/04/24

# Use biomaRt to download annotations and sequences of all human protein coding
# transcripts (non MT) together with transcript annotations

library(tidyverse)
library(biomaRt)

options(timeout = 300)
httr::set_config(httr::config(connecttimeout = 60, timeout = 300))

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 113)

ensembl_filters <- listFilters(ensembl)
ensembl_attributes <- listAttributes(ensembl)

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
                                                values = "protein_coding",
                                                mart = ensembl,
                                                useCache=FALSE)

# remove:
# transcript with tsl>3
# non protein transcripts (retained introns and non coding exons only composed by UTRs)
# transcripts non aligned to chromosomes
# mitochondrial transcripts
protein_coding_transcripts <- protein_coding_transcripts_featurepage %>% filter(transcript_tsl=="tsl1" |
                                                                                transcript_tsl=="tsl2" |
                                                                                transcript_tsl=="tsl3",
                                                                                !chromosome_name == "MT",
                                                                                str_detect(chromosome_name, regex("^\\d{1,2}$")))

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
# get sequence for every transcript (cycle at batches of size 'stepsize' transcripts to avoid time out)
stepsize <- 10000; start <- 1; end <- stepsize; done <- FALSE

# time tracking
start_time <- Sys.time()

while (!done) {

  # check if end index exceeded dataframe length
  if (end > length(unique(protein_coding_transcripts$ensembl_transcript_id))) {
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

  start <- end+1; end <- end+stepsize
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
