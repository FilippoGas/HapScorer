# workflow/scripts/parse_ensembl_bulk.py

"""
Parses Ensembl GTF and FASTA files to extract coding sequences, calculate
Untranslated Regions (UTRs), and generate BioMart-formatted annotations.

Filters transcripts based on predefined biotypes and custom rules (e.g., Transcript
Support Level, chromosome exclusions), maps structural UTR coordinates to their
containing exons, and integrates UniProt cross-references to produce a comprehensive
transcriptomic annotation dataset.

Supports dual execution modes:
1. Native Snakemake workflow integration via the `script:` directive.
2. Standalone Command Line Interface (CLI) execution via `argparse` for manual debugging.

Expected Snakemake Inputs (or CLI Arguments):
    - **gtf** (`--gtf`): Ensembl GTF annotation file.
    - **cdna** (`--cdna`): FASTA file containing cDNA sequences.
    - **cds** (`--cds`): FASTA file containing Coding DNA Sequences.
    - **pep** (`--pep`): FASTA file containing peptide (amino acid) sequences.
    - **uniprot** (`--uniprot`): TSV file for UniProt isoform mappings.

Expected Snakemake Outputs (or CLI Arguments):
    - **sequences** (`--out-sequences`): CSV file with extracted 5' UTR, coding, and 3' UTR sequences.
    - **sequences_aa** (`--out-sequences-aa`): CSV file with formatted amino acid sequences.
    - **annotations** (`--out-annotations`): CSV file containing exon-level genomic coordinates.

Optional Parameters:
    - **biotype** (`--biotype`): Target transcript biotype filter (default: 'protein_coding').
    - **allowed_tsl** (`--allowed-tsl`): List of allowed TSL values (e.g. 1 2).
    - **excluded_chromosomes** (`--excluded-chromosomes`): List of chromosomes to exclude.
    - **chromosome_regex** (`--chromosome-regex`): Regex pattern to filter chromosome names.
    - **ncbi_genetic_code** (`--ncbi-genetic-code`): NCBI Genetic Code Table ID used to derive translation termination codons (default: 1).
"""

import argparse
import gzip
import logging
import re
import sys
import warnings
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set, Union

import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable
from utils import init_snakemake_log, setup_logger, setup_warnings

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake = Snakemake()

logger = setup_logger(__name__)
setup_warnings()


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================


def parse_args() -> argparse.Namespace:
    """
    Parses command-line arguments for the Ensembl bulk parsing workflow.

    Returns:
        args (:obj:`argparse.Namespace`): An object containing the parsed execution parameters.
    """
    parser = argparse.ArgumentParser(description="Parses Ensembl GTF and FASTA files for transcript annotations.")
    parser.add_argument("--gtf", required=True, help="Path to Ensembl GTF file.")
    parser.add_argument("--cdna", required=True, help="Path to cDNA FASTA file.")
    parser.add_argument("--cds", required=True, help="Path to CDS FASTA file.")
    parser.add_argument("--pep", required=True, help="Path to peptide FASTA file.")
    parser.add_argument("--uniprot", required=True, help="Path to UniProt mappings TSV file.")
    parser.add_argument("--out-sequences", required=True, help="Output CSV path for extracted UTR/CDS sequences.")
    parser.add_argument("--out-sequences-aa", required=True, help="Output CSV path for amino acid sequences.")
    parser.add_argument("--out-annotations", required=True, help="Output CSV path for BioMart-like annotations.")
    parser.add_argument("--biotype", default="protein_coding", help="Transcript biotype filter.")
    parser.add_argument("--allowed-tsl", nargs="*", help="Allowed TSL values (e.g., 1 2).")
    parser.add_argument("--excluded-chromosomes", nargs="*", help="Chromosomes to exclude.")
    parser.add_argument("--chromosome-regex", default="", help="Regex pattern to filter chromosomes.")
    parser.add_argument("--ncbi-genetic-code", type=int, default=1, help="NCBI Genetic Code Table ID (default: 1).")
    parser.add_argument("--verbose", action="store_true", help="Enable detailed DEBUG logging output.")
    return parser.parse_args()


def extract_gtf_attribute(attribute_series: pd.Series, key: str) -> pd.Series:
    """
    Extracts the value associated with a specified key from a GTF attribute column.

    This function utilizes regular expression matching to parse the 9th column
    of a standard GTF format file, capturing the string value enclosed in
    double quotes for a given attribute key.

    Args:
        attribute_series (:obj:`pandas.Series`): A pandas Series containing raw GTF
            attribute strings (e.g., 'gene_id "ENSG000..."; transcript_id "ENST000...";').
        key (:obj:`str`): The specific attribute key to extract (e.g., "transcript_id" or "gene_name").

    Returns:
        extracted_series (:obj:`pandas.Series`): A pandas Series containing the extracted values as strings.
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
        fasta_path (:obj:`str` or :obj:`pathlib.Path`): The file path to the FASTA file.
            Can be uncompressed (.fa, .fasta) or gzip-compressed (.gz).
        id_regex (:obj:`str`): A regular expression pattern used to capture the target
            identifier from the FASTA record description.
        group (:obj:`int`, optional): The regular expression capture group to use as
            the dictionary key. Defaults to 1.

    Returns:
        fasta_dict (:obj:`dict` of :obj:`str`): A dictionary where keys are the extracted sequence identifiers
            and values are the corresponding biological sequences.
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
        gtf (:obj:`pandas.DataFrame`): The parsed GTF data containing feature annotations.
        feature_name (:obj:`str`): The UTR feature type to extract as defined in the GTF
            (e.g., "five_prime_utr" or "three_prime_utr").
        start_col (:obj:`str`): The target column name to assign to the UTR start coordinate.
        end_col (:obj:`str`): The target column name to assign to the UTR end coordinate.
        final_tx_ids (:obj:`set` of :obj:`str`): A collection of validated Ensembl transcript IDs
            to retain during filtering.
        exons_df (:obj:`pandas.DataFrame`): A DataFrame containing the standardized transcript IDs
            and their associated exon genomic boundaries ('exon_chrom_start', 'exon_chrom_end').

    Returns:
        mapped_utr (:obj:`pandas.DataFrame`): A mapped DataFrame containing the transcript ID, the containing
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


def get_two_base_map(table_id: int = 1) -> Dict[str, str]:
    """
    Derives 2-base partial codon mappings (4-fold degenerate site prefixes) dynamically.

    Evaluates the specified NCBI Genetic Code Table to identify all 2-nucleotide
    prefix combinations where every possible 3rd nucleotide (A, C, G, T) translates to
    the exact same amino acid residue. This enables unambiguous translation of partial
    terminal codons in incomplete or frame-trimmed coding DNA sequences (CDS).

    Args:
        table_id (:obj:`int`, optional): NCBI Genetic Code Table ID used for codon
            translation mapping (default: 1).

    Returns:
        two_base_map (:obj:`dict` of :obj:`str`): A dictionary mapping 2-base nucleotide
            prefix strings (e.g., 'CC') to their corresponding single-letter amino acid
            symbol (e.g., 'P').
    """
    table = CodonTable.unambiguous_dna_by_id[table_id]
    forward_table = table.forward_table
    bases = ["A", "C", "G", "T"]
    two_base_map = {}

    for b1 in bases:
        for b2 in bases:
            prefix = b1 + b2
            # Extract amino acid residues for all 4 possible 3rd positions
            residues = {forward_table.get(prefix + b3) for b3 in bases}
            # If all 4 codons map to the exact same valid residue (excluding stop codons)
            if len(residues) == 1 and None not in residues:
                two_base_map[prefix] = residues.pop()

    return two_base_map


# =============================================================================
# MAIN LOGIC
# =============================================================================


def run_pipeline(
    input_gtf: Union[str, Path],
    input_cdna: Union[str, Path],
    input_cds: Union[str, Path],
    input_pep: Union[str, Path],
    input_uniprot: Union[str, Path],
    output_sequences: Union[str, Path],
    output_sequences_aa: Union[str, Path],
    output_annotations: Union[str, Path],
    p_biotype: str = "protein_coding",
    p_filters: Optional[Dict[str, Any]] = None,
    p_ncbi_genetic_code: int = 1,
) -> None:
    """
    Executes the complete Ensembl GTF/FASTA parsing and annotation pipeline.

    Filters transcripts based on biotype and custom filtering criteria, computes UTR coordinates
    mapped to parent exons, extracts nucleotide and peptide sequences, and exports standardized CSV files.

    Args:
        input_gtf (:obj:`str` or :obj:`pathlib.Path`): Path to Ensembl GTF annotation file.
        input_cdna (:obj:`str` or :obj:`pathlib.Path`): Path to cDNA FASTA file.
        input_cds (:obj:`str` or :obj:`pathlib.Path`): Path to CDS FASTA file.
        input_pep (:obj:`str` or :obj:`pathlib.Path`): Path to peptide FASTA file.
        input_uniprot (:obj:`str` or :obj:`pathlib.Path`): Path to UniProt mappings TSV file.
        output_sequences (:obj:`str` or :obj:`pathlib.Path`): Output CSV path for extracted nucleotide sequences.
        output_sequences_aa (:obj:`str` or :obj:`pathlib.Path`): Output CSV path for amino acid sequences.
        output_annotations (:obj:`str` or :obj:`pathlib.Path`): Output CSV path for BioMart-formatted annotations.
        p_biotype (:obj:`str`): Target transcript biotype filter.
        p_filters (:obj:`dict`, optional): Dictionary containing transcript filtering rules.
        p_ncbi_genetic_code (:obj:`int`, optional): NCBI Genetic Code Table ID used to derive stop codons (default: 1).
    """
    # Derive stop codons dynamically from NCBI Genetic Code Table
    stop_codons_set = set(CodonTable.unambiguous_dna_by_id[p_ncbi_genetic_code].stop_codons)

    logger.info(
        f"Initiating bulk parsing of Ensembl reference data using NCBI Genetic Code Table {p_ncbi_genetic_code}..."
    )

    # Load GTF Annotations
    logger.info(f"Loading GTF annotation dataset from '{input_gtf}'...")
    gtf_cols = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    # Compression is inferred automatically if the file ends in .gz
    gtf_data = pd.read_csv(input_gtf, sep="\t", comment="#", names=gtf_cols, dtype={"seqname": str})

    logger.info(f"Extracting transcript features filtered by biotype '{p_biotype}'...")
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
            logger.info(f"Filtering transcripts by allowed TSL: {', '.join(allowed_tsl_vec)}")
            protein_coding_transcripts = protein_coding_transcripts[
                protein_coding_transcripts["transcript_tsl"].isin(allowed_tsl_vec)
            ]

        # Exclude chromosomes
        excluded_chroms = p_filters.get("excluded_chromosomes")
        if excluded_chroms:
            logger.info(f"Excluding chromosomes: {', '.join(excluded_chroms)}")
            protein_coding_transcripts = protein_coding_transcripts[
                ~protein_coding_transcripts["chromosome_name"].isin(excluded_chroms)
            ]

        # Regex chromosome filtering
        chrom_regex = p_filters.get("chromosome_regex")
        if chrom_regex and str(chrom_regex).strip():
            logger.info(f"Filtering chromosome names with regex pattern: '{chrom_regex}'\n")
            protein_coding_transcripts = protein_coding_transcripts[
                protein_coding_transcripts["chromosome_name"].str.contains(chrom_regex, regex=True)
            ]
    else:
        logger.info("No filtering rules provided. Retaining all retrieved transcripts.\n")

    # Load FASTA sequences
    logger.info("Loading FASTA sequence datasets (cDNA, CDS, Peptide)...")

    cdna_seqs = load_fasta_dict(input_cdna, r"^(ENST\d+)")
    cds_seqs = load_fasta_dict(input_cds, r"^(ENST\d+)")
    # For peptides, header looks like ENSP... transcript:ENST...
    pep_seqs = load_fasta_dict(input_pep, r"transcript:(ENST\d+)", group=1)

    # Intersection of valid transcripts
    filtered_tx_ids = set(protein_coding_transcripts["ensembl_transcript_id"])
    valid_tx = filtered_tx_ids.intersection(cdna_seqs.keys(), cds_seqs.keys(), pep_seqs.keys())
    valid_tx = sorted(list(valid_tx))

    missing_fasta_txs = filtered_tx_ids - set(valid_tx)
    if missing_fasta_txs:
        sample_missing = sorted(missing_fasta_txs)
        sample_str = ", ".join(sample_missing[:5]) + (
            f" and {len(sample_missing) - 5} more" if len(sample_missing) > 5 else ""
        )
        warnings.warn(
            f"Detected {len(missing_fasta_txs)} transcript(s) in GTF annotation dataset lacking complete "
            f"sequence representation across cDNA, CDS, or peptide FASTA files (e.g., {sample_str}). "
            f"These will be omitted from downstream datasets.",
            UserWarning,
        )

    logger.info(
        f"Identified {len(valid_tx)} valid protein-coding transcripts with complete FASTA sequence concordance."
    )

    # Calculate UTRs and format sequences
    logger.info("Extracting 5' UTR and 3' UTR boundaries and formatting wild-type amino acid sequences...")
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

                # Handle partial amino acids and Stop Codons for peptides
                partial_aa = get_two_base_map(p_ncbi_genetic_code).get(cds[-2:], "") if len(cds) % 3 == 2 else ""
                stop_symbol = "*" if len(cds) >= 3 and len(cds) % 3 == 0 and cds[-3:] in stop_codons_set else ""
                final_pep = pep + partial_aa + stop_symbol
                aa_records.append({"sequence": final_pep, "ensembl_transcript_id": tx})
        else:
            logger.debug(f"Transcript '{tx}' CDS sequence could not be aligned within its cDNA sequence. Skipping.")

    sequences = pd.DataFrame(seq_records)
    sequences_aa = pd.DataFrame(aa_records)

    final_tx_ids = set(sequences["ensembl_transcript_id"])

    # Extraction of Structural Coordinates
    logger.info("Expanding structural coordinates to exon-level boundaries (BioMart format)...")

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
    logger.info(f"Loading UniProt isoform mappings from '{input_uniprot}'...")
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
    logger.debug(f"Cross-referenced {len(uniprot_df)} target transcript IDs with valid UniProt isoform mappings.")

    # Generate Final Annotations
    logger.info("Generating consolidated transcript annotations DataFrame...")
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
    logger.info("Saving results to CSV...")
    sequences.to_csv(output_sequences, index=False)
    sequences_aa.to_csv(output_sequences_aa, index=False)
    protein_coding_transcripts.to_csv(output_annotations, index=False, na_rep="NA")
    logger.info(
        f"Successfully exported extracted sequences ({len(sequences)} records) to '{output_sequences}', "
        f"wild-type amino acid sequences ({len(sequences_aa)} records) to '{output_sequences_aa}', "
        f"and exon-level annotations ({len(protein_coding_transcripts)} records) to '{output_annotations}'."
    )
    logger.info("Ensembl bulk reference data parsing pipeline completed successfully.")


if __name__ == "__main__":
    is_snakemake = False

    if "snakemake" in globals():
        is_snakemake = True
        init_snakemake_log(snakemake, base_dir=Path.cwd())
    else:
        args = parse_args()

    input_gtf = snakemake.input["gtf"] if is_snakemake else args.gtf
    input_cdna = snakemake.input["cdna"] if is_snakemake else args.cdna
    input_cds = snakemake.input["cds"] if is_snakemake else args.cds
    input_pep = snakemake.input["pep"] if is_snakemake else args.pep
    input_uniprot = snakemake.input["uniprot"] if is_snakemake else args.uniprot

    output_sequences = snakemake.output["sequences"] if is_snakemake else args.out_sequences
    output_sequences_aa = snakemake.output["sequences_aa"] if is_snakemake else args.out_sequences_aa
    output_annotations = snakemake.output["annotations"] if is_snakemake else args.out_annotations

    p_biotype = snakemake.params.get("biotype", "protein_coding") if is_snakemake else args.biotype
    p_filters = (
        snakemake.params.get("filters", {})
        if is_snakemake
        else {
            "allowed_tsl": args.allowed_tsl,
            "excluded_chromosomes": args.excluded_chromosomes,
            "chromosome_regex": args.chromosome_regex,
        }
    )
    p_ncbi_genetic_code = snakemake.params.get("ncbi_genetic_code", 1) if is_snakemake else args.ncbi_genetic_code

    if not is_snakemake and args.verbose:
        logger.setLevel(logging.DEBUG)

    run_pipeline(
        input_gtf=input_gtf,
        input_cdna=input_cdna,
        input_cds=input_cds,
        input_pep=input_pep,
        input_uniprot=input_uniprot,
        output_sequences=output_sequences,
        output_sequences_aa=output_sequences_aa,
        output_annotations=output_annotations,
        p_biotype=p_biotype,
        p_filters=p_filters,
        p_ncbi_genetic_code=p_ncbi_genetic_code,
    )
