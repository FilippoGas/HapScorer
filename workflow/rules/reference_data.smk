# workflow/rules/reference_data.smk

# Fallback default settings for reference resources
res_cfg = config.setdefault("resources", {})
res_cfg.setdefault("ref_data_dir", "resources/reference_data")

REF_DIR = Path(res_cfg["ref_data_dir"])

# --- Fallback default settings for reference data module ---
ref_cfg = config.setdefault("reference_data", {})
ENSEMBL_BIOMART = ref_cfg.setdefault("ensembl_biomart", "ENSEMBL_MART_ENSEMBL")
ENSEMBL_DATASET = ref_cfg.setdefault("ensembl_dataset", "hsapiens_gene_ensembl")
ENSEMBL_HOST = ref_cfg.setdefault("ensembl_host", "https://oct2024.archive.ensembl.org")

# --- Rule-specific settings for download_wt_sequences ---
dl_seq_cfg = ref_cfg.setdefault("download_wt_sequences", {})
TRANSCRIPT_BIOTYPE = dl_seq_cfg.setdefault("transcript_biotype", "protein_coding")
FILTERS = dl_seq_cfg.setdefault(
    "filters",
    {
        "allowed_tsl": ["tsl1", "tsl2", "tsl3"],
        "excluded_chromosomes": ["MT"],
        "chromosome_regex": r"^\d{1,2}$",
    },
)
BATCH_SIZE = dl_seq_cfg.setdefault("batch_size", 1000)
SLEEP_DELAY = dl_seq_cfg.setdefault("sleep_delay", 2)


rule download_wt_sequences:
    """
    Download wild-type protein-coding sequence annotations, CDS, and AA sequences from Ensembl.

    Executes the biomaRt retrieval script to generate reference transcript tables.

    :output sequences: CSV containing combined 5'UTR, CDS, and 3'UTR nucleotide sequences.
    :output sequences_aa: CSV containing translated wild-type amino acid sequences.
    :output annotations: CSV containing protein-coding transcript genomic annotations.

    """
    output:
        sequences=REF_DIR / "wt_cds.csv",
        sequences_aa=REF_DIR / "wt_aa.csv",
        annotations=REF_DIR / "protein_coding_transcripts.csv",
    log:
        "logs/reference_data/download_wt_sequences.log",
    benchmark:
        "benchmarks/reference_data/download_wt_sequences.tsv"
    retries: 3
    conda:
        "../envs/reference_data.yaml"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=60,
    params:
        host=ENSEMBL_HOST,
        biomart=ENSEMBL_BIOMART,
        dataset=ENSEMBL_DATASET,
        biotype=TRANSCRIPT_BIOTYPE,
        filters=FILTERS,
        batch_size=BATCH_SIZE,
        sleep_delay=SLEEP_DELAY,
    message:
        "Downloading Ensembl wild-type sequences and annotations via biomaRt"
    script:
        "../scripts/sequences_download.R"
