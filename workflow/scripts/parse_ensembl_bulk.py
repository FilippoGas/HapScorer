# workflow/scripts/parse_ensembl_bulk.py

"""
Parses Ensembl GTF and FASTA files to extract coding sequences, calculate
Untranslated Regions (UTRs), and generate BioMart-formatted annotations.

This script is designed to be executed as a node within a Snakemake workflow.
It filters transcripts based on predefined biotypes and custom rules (e.g., Transcript
Support Level, chromosome exclusions), maps structural UTR coordinates to their
containing exons, and integrates UniProt cross-references to produce a comprehensive
transcriptomic annotation dataset.

Expected Snakemake Inputs:
    - gtf: Ensembl GTF annotation file.
    - cdna: FASTA file containing cDNA sequences.
    - cds: FASTA file containing Coding DNA Sequences.
    - pep: FASTA file containing peptide (amino acid) sequences.
    - uniprot: TSV file for UniProt isoform mappings.

Expected Snakemake Outputs:
    - sequences: CSV file with extracted 5' UTR, coding, and 3' UTR sequences.
    - sequences_aa: CSV file with formatted amino acid sequences (including stop codons).
    - annotations: CSV file containing exon-level genomic coordinates, transcript metadata,
      and UniProt cross-references formatted for BioMart compatibility.
"""

import gzip
import re
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Set, Union

import pandas as pd
from Bio import SeqIO

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake = Snakemake()

# Redirect standard output and standard error to the Snakemake log file
log_file = open(snakemake.log[0], "w", buffering=1)
sys.stdout = log_file
sys.stderr = log_file


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================


def extract_gtf_attribute(attribute_series: pd.Series, key: str) -> pd.Series:
    """
    Extracts the value associated with a specified key from a GTF attribute column.

    This function utilizes regular expression matching to parse the 9th column
    of a standard GTF format file, capturing the string value enclosed in
    double quotes for a given attribute key.

    Args:
        attribute_series (pd.Series): A pandas Series containing the raw GTF
            attribute strings (e.g., 'gene_id "ENSG000..."; transcript_id "ENST000...";').
        key (str): The specific attribute key to extract (e.g., "transcript_id" or "gene_name").

    Returns:
        pd.Series: A pandas Series containing the extracted values as strings.
        Missing keys will result in NaN values.
    """
    # Pattern explanation: matches `key "value";` and captures `value`
    pattern = rf'{key} "([^"]+)"'
    return attribute_series.str.extract(pattern, expand=False)


def load_fasta_dict(fasta_path: Union[str, Path], id_regex: str, group: int = 1) -> Dict[str, str]:
    """
    Parses a FASTA file and constructs a dictionary mapping sequence identifiers
    to their corresponding nucleotide or amino acid sequences.

    This function automatically detects and handles gzip-compressed files. It extracts
    a specific sequence identifier from the FASTA description header using a provided
    regular expression, allowing for robust ID matching regardless of the database
    header format.

    Args:
        fasta_path (str or pathlib.Path): The file path to the FASTA file.
            Can be uncompressed (.fa, .fasta) or gzip-compressed (.gz).
        id_regex (str): A regular expression pattern used to capture the target
            identifier from the FASTA record description.
        group (int, optional): The regular expression capture group to use as
            the dictionary key. Defaults to 1.

    Returns:
        dict: A dictionary where keys are the extracted sequence identifiers (str)
        and values are the corresponding biological sequences (str).
    """
    fasta_dict = {}

    # It detects whether the file is compressed and uses the appropriate opening mode for text ("rt")
    is_gz = str(fasta_path).endswith(".gz")
    open_func = gzip.open if is_gz else open

    with open_func(fasta_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            match = re.search(id_regex, record.description)
            if match:
                seq_id = match.group(group)
                fasta_dict[seq_id] = str(record.seq)

    return fasta_dict


def map_utr(
    gtf: pd.DataFrame, feature_name: str, start_col: str, end_col: str, final_tx_ids: Set[str], exons_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Extracts specified Untranslated Region (UTR) features from a GTF dataset
    and maps them strictly to their containing exon coordinates.

    This function filters the GTF DataFrame for a specific UTR feature, standardizes
    the Ensembl transcript IDs by stripping version numbers, and retains only
    the transcripts present in a validated subset. It then performs an inner join
    with exon spatial data to identify and retain only those UTR boundaries that
    fall correctly within the genomic coordinates of their parent exon.

    Args:
        gtf (pd.DataFrame): The parsed GTF data containing feature annotations.
        feature_name (str): The UTR feature type to extract as defined in the GTF
            (e.g., "five_prime_utr" or "three_prime_utr").
        start_col (str): The target column name to assign to the UTR start coordinate.
        end_col (str): The target column name to assign to the UTR end coordinate.
        final_tx_ids (set or list): A collection of validated Ensembl transcript IDs
            to retain during filtering.
        exons_df (pd.DataFrame): A DataFrame containing the standardized transcript IDs
            and their associated exon genomic boundaries ('exon_chrom_start', 'exon_chrom_end').

    Returns:
        pd.DataFrame: A mapped DataFrame containing the transcript ID, the containing
        exon boundaries, and the specific UTR coordinates.
    """
    utr = gtf[gtf["feature"] == feature_name].copy()
    utr["ensembl_transcript_id"] = extract_gtf_attribute(utr["attribute"], "transcript_id").str.replace(
        r"\.\d+$", "", regex=True
    )
    utr = utr[utr["ensembl_transcript_id"].isin(final_tx_ids)][["ensembl_transcript_id", "start", "end"]].rename(
        columns={"start": start_col, "end": end_col}
    )

    # Merge with exons to find the containing exon
    merged = pd.merge(utr, exons_df, on="ensembl_transcript_id")
    # Keep only where UTR is within the exon coordinates
    contained = merged[
        (merged[start_col] >= merged["exon_chrom_start"]) & (merged[end_col] <= merged["exon_chrom_end"])
    ]
    return contained[["ensembl_transcript_id", "exon_chrom_start", "exon_chrom_end", start_col, end_col]]


# =============================================================================
# MAIN LOGIC
# =============================================================================


def main():
    """
    Primary execution function for Snakemake integration.
    """
    # Retrieve inputs from Snakemake input
    input_gtf = snakemake.input["gtf"]
    input_cdna = snakemake.input["cdna"]
    input_cds = snakemake.input["cds"]
    input_pep = snakemake.input["pep"]
    input_uniprot = snakemake.input["uniprot"]

    # Retrieve outputs from Snakemake output
    output_sequences = snakemake.output["sequences"]
    output_sequences_aa = snakemake.output["sequences_aa"]
    output_annotations = snakemake.output["annotations"]

    # Retrieve execution parameters params
    p_biotype = snakemake.params["biotype"]
    p_filters = snakemake.params["filters"]

    print("Parsing local Ensembl GTF and FASTA files...\n")

    # Load GTF Annotations
    print("Loading GTF annotations...")
    gtf_cols = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    # Compression is inferred automatically if the file ends in .gz
    gtf_data = pd.read_csv(input_gtf, sep="\t", comment="#", names=gtf_cols, dtype={"seqname": str})

    print("Extracting transcript features...")
    # Filter for transcripts
    tx_gtf = gtf_data[gtf_data["feature"] == "transcript"].copy()

    # Extract attributes needed for transcripts
    tx_gtf["transcript_id"] = extract_gtf_attribute(tx_gtf["attribute"], "transcript_id")
    tx_gtf["gene_id"] = extract_gtf_attribute(tx_gtf["attribute"], "gene_id")
    tx_gtf["gene_name"] = extract_gtf_attribute(tx_gtf["attribute"], "gene_name")
    tx_gtf["transcript_support_level"] = extract_gtf_attribute(tx_gtf["attribute"], "transcript_support_level")
    tx_gtf["transcript_biotype"] = extract_gtf_attribute(tx_gtf["attribute"], "transcript_biotype")

    # Filter by biotype
    tx_gtf: pd.DataFrame = tx_gtf[tx_gtf["transcript_biotype"] == p_biotype].copy()

    # Format transcript dataframe
    tx_gtf["ensembl_transcript_id"] = tx_gtf["transcript_id"].str.replace(r"\.\d+$", "", regex=True)
    tx_gtf["ensembl_gene_id"] = tx_gtf["gene_id"].str.replace(r"\.\d+$", "", regex=True)
    tx_gtf["strand"] = tx_gtf["strand"].map({"-": -1, "+": 1})
    # Clean TSL (e.g., "1 (assigned to previous version)" -> "1")
    tx_gtf["transcript_tsl"] = tx_gtf["transcript_support_level"].str.extract(r"(\d+)", expand=False)

    protein_coding_transcripts = tx_gtf[
        [
            "ensembl_transcript_id",
            "ensembl_gene_id",
            "gene_name",
            "transcript_tsl",
            "seqname",
            "start",
            "end",
            "strand",
            "transcript_biotype",
        ]
    ].rename(
        columns={
            "gene_name": "hgnc_symbol",
            "seqname": "chromosome_name",
            "start": "transcript_start",
            "end": "transcript_end",
        }
    )

    # Apply declarative filters
    if p_filters and isinstance(p_filters, dict):
        # Filter by TSL
        allowed_tsl = p_filters.get("allowed_tsl")
        if allowed_tsl:
            allowed_tsl_vec = [str(t).replace("tsl", "") for t in allowed_tsl]
            print(f"Filtering transcripts by allowed TSL: {', '.join(allowed_tsl_vec)}")
            protein_coding_transcripts = protein_coding_transcripts[
                protein_coding_transcripts["transcript_tsl"].isin(allowed_tsl_vec)
            ]

        # Exclude chromosomes
        excluded_chroms = p_filters.get("excluded_chromosomes")
        if excluded_chroms:
            print(f"Excluding chromosomes: {', '.join(excluded_chroms)}")
            protein_coding_transcripts = protein_coding_transcripts[
                ~protein_coding_transcripts["chromosome_name"].isin(excluded_chroms)
            ]

        # Regex chromosome filtering
        chrom_regex = p_filters.get("chromosome_regex")
        if chrom_regex and str(chrom_regex).strip():
            print(f"Filtering chromosome names with regex pattern: '{chrom_regex}'\n")
            protein_coding_transcripts = protein_coding_transcripts[
                protein_coding_transcripts["chromosome_name"].str.contains(chrom_regex, regex=True)
            ]
    else:
        print("No filtering rules provided. Retaining all retrieved transcripts.\n")

    # Load FASTA sequences
    print("Loading FASTA sequences...")

    cdna_seqs = load_fasta_dict(input_cdna, r"^(ENST\d+)")
    cds_seqs = load_fasta_dict(input_cds, r"^(ENST\d+)")
    # For peptides, header looks like ENSP... transcript:ENST...
    pep_seqs = load_fasta_dict(input_pep, r"transcript:(ENST\d+)", group=1)

    # Intersection of valid transcripts
    filtered_tx_ids = set(protein_coding_transcripts["ensembl_transcript_id"])
    valid_tx = filtered_tx_ids.intersection(cdna_seqs.keys(), cds_seqs.keys(), pep_seqs.keys())
    valid_tx = sorted(list(valid_tx))

    # Calculate UTRs and format sequences
    print("Extracting 5'UTR and 3'UTR regions...")
    seq_records = []
    aa_records = []

    for tx in valid_tx:
        cdna = cdna_seqs[tx]
        cds = cds_seqs[tx]
        pep = pep_seqs[tx]

        # Find CDS location within cDNA
        cds_start_idx = cdna.find(cds)
        if cds_start_idx != -1:
            utr5 = cdna[:cds_start_idx]
            utr3 = cdna[cds_start_idx + len(cds) :]

            if utr5 and utr3 and cds:
                seq_records.append({"ensembl_transcript_id": tx, "5utr": utr5, "coding": cds, "3utr": utr3})

                # Handle Stop Codons for peptides
                last_codon = cds[-3:]
                final_pep = pep + "*" if last_codon in ["TAA", "TAG", "TGA"] else pep
                aa_records.append({"sequence": final_pep, "ensembl_transcript_id": tx})

    sequences = pd.DataFrame(seq_records)
    sequences_aa = pd.DataFrame(aa_records)

    final_tx_ids = set(sequences["ensembl_transcript_id"])

    # Extraction of Structural Coordinates
    print("Expanding output to exon-level (BioMart format)...")

    # Genes
    genes_gtf = gtf_data[gtf_data["feature"] == "gene"].copy()
    genes_gtf["ensembl_gene_id"] = extract_gtf_attribute(genes_gtf["attribute"], "gene_id").str.replace(
        r"\.\d+$", "", regex=True
    )
    genes_df = (
        genes_gtf[["ensembl_gene_id", "start", "end"]]
        .rename(columns={"start": "start_position", "end": "end_position"})
        .drop_duplicates()
    )

    # Exons
    exons_gtf = gtf_data[gtf_data["feature"] == "exon"].copy()
    exons_gtf["ensembl_transcript_id"] = extract_gtf_attribute(exons_gtf["attribute"], "transcript_id").str.replace(
        r"\.\d+$", "", regex=True
    )
    exons_df = exons_gtf[exons_gtf["ensembl_transcript_id"].isin(final_tx_ids)][
        ["ensembl_transcript_id", "start", "end"]
    ].rename(columns={"start": "exon_chrom_start", "end": "exon_chrom_end"})

    # Transcript lengths
    exons_df["exon_len"] = exons_df["exon_chrom_end"] - exons_df["exon_chrom_start"] + 1
    tx_lengths = (
        exons_df.groupby("ensembl_transcript_id")["exon_len"]
        .sum()
        .reset_index()
        .rename(columns={"exon_len": "transcript_length"})
    )
    exons_df = exons_df.drop(columns=["exon_len"])

    # Peptides
    cds_gtf = gtf_data[gtf_data["feature"] == "CDS"].copy()
    cds_gtf["ensembl_transcript_id"] = extract_gtf_attribute(cds_gtf["attribute"], "transcript_id").str.replace(
        r"\.\d+$", "", regex=True
    )
    cds_gtf["ensembl_peptide_id"] = extract_gtf_attribute(cds_gtf["attribute"], "protein_id").str.replace(
        r"\.\d+$", "", regex=True
    )
    pep_df = cds_gtf[cds_gtf["ensembl_transcript_id"].isin(final_tx_ids)][
        ["ensembl_transcript_id", "ensembl_peptide_id"]
    ].drop_duplicates()

    # UTR Coordinates Mapping
    utr5_mapped = map_utr(gtf_data, "five_prime_utr", "5_utr_start", "5_utr_end", final_tx_ids, exons_df)
    utr3_mapped = map_utr(gtf_data, "three_prime_utr", "3_utr_start", "3_utr_end", final_tx_ids, exons_df)

    # Combine Exons with UTRs
    exons_with_utr = pd.merge(
        exons_df, utr5_mapped, on=["ensembl_transcript_id", "exon_chrom_start", "exon_chrom_end"], how="left"
    )
    exons_with_utr = pd.merge(
        exons_with_utr, utr3_mapped, on=["ensembl_transcript_id", "exon_chrom_start", "exon_chrom_end"], how="left"
    ).drop_duplicates()

    # Load UniProt mappings
    print("Loading UniProt mappings...")
    uniprot_cols = [
        "gene_stable_id",
        "transcript_stable_id",
        "protein_stable_id",
        "xref",
        "db_name",
        "info_type",
        "source_identity",
        "xref_identity",
        "linkage_type",
    ]
    uniprot_df = pd.read_csv(input_uniprot, sep="\t", names=uniprot_cols, header=0, dtype=str)

    uniprot_df = uniprot_df[uniprot_df["db_name"] == "Uniprot_isoform"].copy()
    uniprot_df["ensembl_transcript_id"] = uniprot_df["transcript_stable_id"].str.replace(r"\.\d+$", "", regex=True)
    uniprot_df = (
        uniprot_df[uniprot_df["ensembl_transcript_id"].isin(final_tx_ids)][["ensembl_transcript_id", "xref"]]
        .rename(columns={"xref": "uniprot_isoform"})
        .dropna()
    )

    # Aggregate multiple isoforms with semicolon
    uniprot_df = (
        uniprot_df.groupby("ensembl_transcript_id")["uniprot_isoform"]
        .apply(lambda x: ";".join(x.unique()))
        .reset_index()
    )

    # Generate Final Annotations
    print("Generating Final Annotations...")
    # Filter to final transcripts
    pct = protein_coding_transcripts[protein_coding_transcripts["ensembl_transcript_id"].isin(final_tx_ids)].copy()

    # Joins
    pct = pd.merge(pct, genes_df, on="ensembl_gene_id", how="left")
    pct = pd.merge(pct, tx_lengths, on="ensembl_transcript_id", how="left")
    pct = pd.merge(pct, pep_df, on="ensembl_transcript_id", how="left")
    pct = pd.merge(pct, exons_with_utr, on="ensembl_transcript_id", how="left")
    pct = pd.merge(pct, uniprot_df, on="ensembl_transcript_id", how="left")

    # Format TSL and missing Uniprot
    pct["transcript_tsl"] = pct["transcript_tsl"].apply(lambda x: f"tsl{x}" if pd.notna(x) else pd.NA)
    pct["uniprot_isoform"] = pct["uniprot_isoform"].fillna("")

    # Enforce Exact BioMart Column Order
    final_cols = [
        "ensembl_gene_id",
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
        "exon_chrom_end",
        "hgnc_symbol",
        "transcript_tsl",
        "chromosome_name",
        "start_position",
        "end_position",
        "strand",
        "transcript_biotype",
        "uniprot_isoform",
    ]

    # Convert floats back to Int64 (to avoid .0 in output for missing UTR coordinates)
    float_cols = ["5_utr_start", "5_utr_end", "3_utr_start", "3_utr_end"]
    for col in float_cols:
        if col in pct.columns:
            pct[col] = pct[col].astype("Int64")

    protein_coding_transcripts = pct[final_cols]

    # Save Results
    print("Saving results to CSV...")
    sequences.to_csv(output_sequences, index=False)
    sequences_aa.to_csv(output_sequences_aa, index=False)
    protein_coding_transcripts.to_csv(output_annotations, index=False, na_rep="NA")

    print("Done.\n")


if __name__ == "__main__":
    main()
