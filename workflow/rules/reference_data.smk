# workflow/rules/reference_data.smk

# Fallback default settings for reference resources
res_cfg = config.setdefault("resources", {})
REF_DIR = Path(res_cfg.setdefault("ref_data_dir", "resources/reference_data"))

# --- Fallback default settings for reference data module ---
ref_cfg = config.setdefault("reference_data", {})
ENSEMBL_RELEASE = ref_cfg.setdefault("release", "113")

# --- Rule-specific settings for parse_ensembl_bulk_sequences ---
dl_seq_cfg = ref_cfg.setdefault("parse_ensembl_bulk_sequences", {})
TRANSCRIPT_BIOTYPE = dl_seq_cfg.setdefault("transcript_biotype", "protein_coding")
FILTERS = dl_seq_cfg.setdefault(
    "filters",
    {
        "allowed_tsl": ["tsl1", "tsl2", "tsl3"],
        "excluded_chromosomes": ["MT"],
        "chromosome_regex": r"^\d{1,2}$",
    },
)
STOP_CODONS = dl_seq_cfg.setdefault("stop_codons", ["TAA", "TAG", "TGA"])


rule download_ensembl_bulk_data:
    """
    Download genomic feature annotations (GTF) and sequence datasets (FASTA)
    from Ensembl FTP servers.

    Executes a Python script to reliably download the data, verify network
    transfer integrity using Ensembl's legacy UNIX `sum` checksums, and
    subsequently generate cryptographic SHA-256 hashes for long-term
    data provenance and storage integrity.

    :output gtf: Compressed GTF file containing genomic feature annotations.
    :output cds: Compressed FASTA file containing coding sequences.
    :output pep: Compressed FASTA file containing amino acid (peptide) sequences.
    :output cdna: Compressed FASTA file containing complementary DNA (transcript) sequences.
    :output uniprot: Compressed TSV file containing UniProt cross-reference mappings.

    """
    output:
        gtf=REF_DIR / f"Homo_sapiens.GRCh38.{ENSEMBL_RELEASE}.gtf.gz",
        cds=REF_DIR / "Homo_sapiens.GRCh38.cds.all.fa.gz",
        pep=REF_DIR / "Homo_sapiens.GRCh38.pep.all.fa.gz",
        cdna=REF_DIR / "Homo_sapiens.GRCh38.cdna.all.fa.gz",
        uniprot=REF_DIR / f"Homo_sapiens.GRCh38.{ENSEMBL_RELEASE}.uniprot.tsv.gz",
    log:
        "logs/reference_data/download_ensembl_bulk_data.log",
    benchmark:
        "benchmarks/reference_data/download_ensembl_bulk_data.tsv"
    retries: 3
    conda:
        "../envs/reference_data.yaml"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=60,
    params:
        release=ENSEMBL_RELEASE,
    message:
        "Downloading and verifying Ensembl bulk reference files (Release {params.release})"
    script:
        "../scripts/download_ensembl_bulk.py"


rule parse_ensembl_bulk_sequences:
    """
    Extract and filter wild-type protein-coding sequence annotations, CDS,
    UTRs, and amino acid sequences from local Ensembl bulk reference files.

    Parses pre-downloaded GTF and FASTA datasets to generate reference
    transcript tables offline, replacing remote API queries.

    :input gtf: Compressed GTF file containing genomic feature annotations.
    :input cds: Compressed FASTA file containing coding sequences.
    :input pep: Compressed FASTA file containing peptide sequences.
    :input cdna: Compressed FASTA file containing transcript sequences.

    :output sequences: CSV containing combined 5'UTR, CDS, and 3'UTR nucleotide sequences.
    :output sequences_aa: CSV containing translated wild-type amino acid sequences.
    :output annotations: CSV containing protein-coding transcript genomic annotations.

    """
    input:
        gtf=REF_DIR / f"Homo_sapiens.GRCh38.{ENSEMBL_RELEASE}.gtf.gz",
        cds=REF_DIR / "Homo_sapiens.GRCh38.cds.all.fa.gz",
        pep=REF_DIR / "Homo_sapiens.GRCh38.pep.all.fa.gz",
        cdna=REF_DIR / "Homo_sapiens.GRCh38.cdna.all.fa.gz",
        uniprot=REF_DIR / f"Homo_sapiens.GRCh38.{ENSEMBL_RELEASE}.uniprot.tsv.gz",
    output:
        sequences=REF_DIR / "wt_cds.csv",
        sequences_aa=REF_DIR / "wt_aa.csv",
        annotations=REF_DIR / "protein_coding_transcripts.csv",
    log:
        "logs/reference_data/parse_ensembl_bulk_sequences.log",
    benchmark:
        "benchmarks/reference_data/parse_ensembl_bulk_sequences.tsv"
    conda:
        "../envs/reference_data.yaml"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=60,
    params:
        biotype=TRANSCRIPT_BIOTYPE,
        filters=FILTERS,
        stop_codons=STOP_CODONS,
    message:
        "Extracting Ensembl wild-type sequences and annotations via pre-downloaded resources"
    script:
        "../scripts/parse_ensembl_bulk.py"
