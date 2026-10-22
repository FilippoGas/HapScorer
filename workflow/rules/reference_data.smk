# workflow/rules/reference_data.smk

# Fallback default settings for reference resources
res_cfg = config.setdefault("resources", {})
res_cfg.setdefault("ref_data_dir", "resources/reference_data")

REF_DIR = Path(res_cfg["ref_data_dir"])


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
    message:
        "Downloading Ensembl wild-type sequences and annotations via biomaRt"
    script:
        "../scripts/sequences_download.R"
