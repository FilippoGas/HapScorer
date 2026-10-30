# workflow/scripts/extract_haplotypes_and_genotypes.py

"""
Extracts haplotypes and genotypes from phased VCF files.

This script parses phased Variant Call Format (VCF) files alongside transcript
annotations to identify unique haplotypic combinations within exonic regions.
It generates comprehensive datasets detailing unique haplotypes per transcript
and their corresponding genotype assignments across all provided samples.

Supports dual execution modes:
1. Native Snakemake workflow integration via the `script:` directive.
2. Standalone Command Line Interface (CLI) execution via `argparse` for manual debugging.

Expected Snakemake Inputs (or CLI Arguments):
    - phased_vcf (--vcf): Path to the phased, bgzipped, and indexed VCF file.
    - annotations (--annotations): Path to the CSV file containing transcript annotations.

Expected Snakemake Outputs (or CLI Arguments):
    - haplotypes_csv (--out-haplotypes): Output CSV path for unique haplotypes metadata.
    - genotypes_csv (--out-genotypes): Output CSV path for sample genotype assignments.

Optional Parameters:
    - threads (--threads): Number of CPU cores to utilize for multiprocessing (default: 1).
"""

import argparse
import multiprocessing as mp
import sys
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, Union

import pandas as pd
import pysam

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake = Snakemake()

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================


def parse_args() -> argparse.Namespace:
    """
    Parses command-line arguments for the haplotype extraction workflow.

    Returns:
        argparse.Namespace: An object containing the parsed execution parameters.
    """
    parser = argparse.ArgumentParser(description="Extracts haplotypes and genotypes from a phased VCF.")
    parser.add_argument("--vcf", required=True, help="Path to the phased VCF file.")
    parser.add_argument("--annotations", required=True, help="CSV file containing transcript annotations.")
    parser.add_argument("--out-haplotypes", required=True, help="Output path for the haplotypes CSV.")
    parser.add_argument("--out-genotypes", required=True, help="Output path for the genotypes CSV.")
    parser.add_argument("--threads", type=int, default=1, help="Number of CPU cores to utilize.")
    return parser.parse_args()


def process_transcript(
    transcript_id: str, annotations_df: pd.DataFrame, vcf_path: Union[str, Path], samples: List[str]
) -> Optional[Tuple[List[Dict[str, Any]], Dict[str, str]]]:
    """
    Processes a single transcript to extract variants, identify unique haplotypes,
    and construct the sample genotype matrix.

    This function interrogates the VCF strictly within the specified exonic boundaries
    of the transcript. It evaluates the mutational landscape across all provided samples,
    aggregating unique allelic combinations into distinct haplotypes.

    Args:
        transcript_id (str): The Ensembl transcript identifier to process.
        annotations_df (pd.DataFrame): DataFrame containing exonic coordinates.
        vcf_path (str or pathlib.Path): File system path to the indexed phased VCF.
        samples (List[str]): List of sample identifiers present in the VCF header.

    Returns:
        Optional[Tuple[List[Dict[str, Any]], Dict[str, str]]]: A tuple containing the
            haplotype metadata records and the transcript-specific genotype matrix.
            Returns None if the transcript lacks exonic annotations.
    """
    vcf_path = Path(vcf_path)

    tx_anno = annotations_df[annotations_df["ensembl_transcript_id"] == transcript_id]
    if tx_anno.empty:
        return None

    gene_id = tx_anno["ensembl_gene_id"].iloc[0]
    chrom = str(tx_anno["chromosome_name"].iloc[0])
    wt_id = f"{gene_id}.0"

    transcript_genotypes = {"ensembl_transcript_id": transcript_id}
    haplotypes_data = []

    for sample in samples:
        transcript_genotypes[sample] = f"{wt_id}/{wt_id}"

    vcf = pysam.VariantFile(vcf_path)
    sample_variants = {sample: {0: [], 1: []} for sample in samples}
    has_variants = False

    for _, row in tx_anno.iterrows():
        start = int(row["exon_chrom_start"]) - 1
        end = int(row["exon_chrom_end"])

        try:
            for record in vcf.fetch(chrom, start, end):
                has_variants = True
                info_ann = record.info.get("ANN", [""])[0]
                variant_type = info_ann.split("|")[1] if "|" in info_ann else ""
                nuc_change = info_ann.split("|")[9] if len(info_ann.split("|")) > 9 else ""

                variant_string = f"{chrom}:{record.pos}.{record.ref}>{record.alts[0]}"
                rsid = record.id if record.id else ""

                for sample in samples:
                    call = record.samples[sample]
                    gt = call["GT"]

                    if gt[0] not in (0, None):
                        sample_variants[sample][0].append(
                            {"var_str": variant_string, "rsid": rsid, "vtype": variant_type, "nchange": nuc_change}
                        )
                    if gt[1] not in (0, None):
                        sample_variants[sample][1].append(
                            {"var_str": variant_string, "rsid": rsid, "vtype": variant_type, "nchange": nuc_change}
                        )
        except ValueError:
            continue

    vcf.close()

    if not has_variants:
        haplotypes_data.append(
            {
                "ensembl_gene_id": gene_id,
                "ensembl_transcript_id": transcript_id,
                "haplotype_gene_id": wt_id,
                "variants": "wt",
                "rsid": "",
                "variant_types": "",
                "nucleotide_change": "",
            }
        )
        return haplotypes_data, transcript_genotypes

    unique_variant_combos = {}
    haplotype_counter = 1

    for sample in samples:
        hap_ids = []
        for allele in [0, 1]:
            vars_list = sample_variants[sample][allele]

            if not vars_list:
                hap_ids.append(wt_id)
                continue

            var_signature = ",".join(v["var_str"] for v in vars_list)

            if var_signature not in unique_variant_combos:
                hap_id = f"{gene_id}.{haplotype_counter}"
                unique_variant_combos[var_signature] = {
                    "ensembl_gene_id": gene_id,
                    "ensembl_transcript_id": transcript_id,
                    "haplotype_gene_id": hap_id,
                    "variants": var_signature,
                    "rsid": ",".join(v["rsid"] for v in vars_list),
                    "variant_types": ",".join(v["vtype"] for v in vars_list),
                    "nucleotide_change": ",".join(v["nchange"] for v in vars_list),
                }
                haplotype_counter += 1
                hap_ids.append(hap_id)
            else:
                hap_ids.append(unique_variant_combos[var_signature]["haplotype_gene_id"])

        transcript_genotypes[sample] = f"{hap_ids[0]}/{hap_ids[1]}"

    haplotypes_data.append(
        {
            "ensembl_gene_id": gene_id,
            "ensembl_transcript_id": transcript_id,
            "haplotype_gene_id": wt_id,
            "variants": "wt",
            "rsid": "",
            "variant_types": "",
            "nucleotide_change": "",
        }
    )
    haplotypes_data.extend(list(unique_variant_combos.values()))

    return haplotypes_data, transcript_genotypes


def run_pipeline(
    vcf_path: Union[str, Path],
    annotations_path: Union[str, Path],
    out_haplotypes: Union[str, Path],
    out_genotypes: Union[str, Path],
    threads: int,
) -> None:
    """
    Executes the full haplotype extraction and sample genotype matrix generation workflow.

    This function coordinates the extraction of sample metadata directly from the VCF
    header, filters transcript annotations for the targeted chromosome, and parallelizes
    transcript processing across multiple CPU worker processes. It aggregates the resulting
    haplotypic variant signatures and genotype assignments, incorporates transcript
    strand orientation, and exports the final standardized datasets to CSV files.

    Args:
        vcf_path (str or pathlib.Path): Path to the phased, bgzipped, and indexed VCF file.
        annotations_path (str or pathlib.Path): Path to the CSV file containing exon-level annotations.
        out_haplotypes (str or pathlib.Path): Output file path for the generated haplotypes metadata CSV.
        out_genotypes (str or pathlib.Path): Output file path for the sample genotype matrix CSV.
        threads (int): Number of parallel CPU worker processes to spawn.

    Returns:
        None
    """
    vcf_path = Path(vcf_path)
    annotations_path = Path(annotations_path)
    out_haplotypes = Path(out_haplotypes)
    out_genotypes = Path(out_genotypes)

    with pysam.VariantFile(vcf_path) as vcf:
        samples = list(vcf.header.samples)

    anno_df = pd.read_csv(annotations_path)
    chrom_from_vcf = vcf_path.name.split("chr")[-1].split(".")[0]
    anno_df = anno_df[anno_df["chromosome_name"].astype(str) == chrom_from_vcf]
    transcripts = anno_df["ensembl_transcript_id"].unique()

    pool = mp.Pool(threads)
    worker_func = partial(process_transcript, annotations_df=anno_df, vcf_path=vcf_path, samples=samples)
    results = pool.map(worker_func, transcripts)
    pool.close()
    pool.join()

    all_haplotypes = []
    all_genotypes = []

    for res in results:
        if res:
            all_haplotypes.extend(res[0])
            all_genotypes.append(res[1])

    df_haplotypes = pd.DataFrame(all_haplotypes)
    df_haplotypes["haplotype_transcript_id"] = (
        df_haplotypes["ensembl_transcript_id"] + "-" + df_haplotypes["haplotype_gene_id"]
    )
    meta_info = anno_df[["ensembl_transcript_id", "hgnc_symbol", "strand"]].drop_duplicates()
    df_haplotypes = pd.merge(df_haplotypes, meta_info, on="ensembl_transcript_id", how="left")

    df_haplotypes = df_haplotypes[
        [
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
    ]

    df_haplotypes.to_csv(out_haplotypes, index=False, na_rep="NA")
    pd.DataFrame(all_genotypes).to_csv(out_genotypes, index=False, na_rep="NA")


# =============================================================================
# MAIN LOGIC
# =============================================================================

if __name__ == "__main__":
    is_snakemake = False

    if "snakemake" in globals():
        is_snakemake = True
        if hasattr(snakemake, "log") and snakemake.log:
            log_f = open(snakemake.log[0], "w", buffering=1)
            sys.stdout = log_f
            sys.stderr = log_f
    else:
        args = parse_args()

    vcf_path = snakemake.input["phased_vcf"] if is_snakemake else args.vcf
    annotations_path = snakemake.input["annotations"] if is_snakemake else args.annotations
    out_haplotypes = snakemake.output["haplotypes_csv"] if is_snakemake else args.out_haplotypes
    out_genotypes = snakemake.output["genotypes_csv"] if is_snakemake else args.out_genotypes
    threads = snakemake.threads if is_snakemake else args.threads

    run_pipeline(
        vcf_path=vcf_path,
        annotations_path=annotations_path,
        out_haplotypes=out_haplotypes,
        out_genotypes=out_genotypes,
        threads=threads,
    )
