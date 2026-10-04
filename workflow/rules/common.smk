# workflow/rules/common.smk

import pandas as pd

# ------------------------------------------------------------------------------
# Sample Metadata Table Loading & Schema Validation
# ------------------------------------------------------------------------------

# Load raw sample table from path defined in config.yaml
samples_df = pd.read_csv(config["samples"], sep="\t")

# Validate table structure against schema
validate(samples_df, schema="../schemas/samples.schema.yaml")

# Normalize var_type column to lower case for case-insensitive indexing
samples_df["var_type"] = samples_df["var_type"].str.lower()

# Set composite index after successful validation
SAMPLES = samples_df.set_index(["sample", "var_type"], drop=False)


# ------------------------------------------------------------------------------
# Dynamic Rule Input & Parameter Functions
# ------------------------------------------------------------------------------


def get_raw_vcf(wildcards):
    """
    Retrieve the raw VCF file path for a specified sample and variant type.

    Queries the global sample metadata registry (`SAMPLES` DataFrame) to extract
    the file system path corresponding to a given sample identifier and variant
    category (e.g., SNP or INDEL) resolved from Snakemake wildcards.

    Args:
        wildcards (snakemake.io.Wildcards): Snakemake wildcards object containing
            the attributes `sample` (sample name/ID) and `var_type` (variant
            classification).

    Returns:
        str: Relative or absolute file path to the target input Variant Call
        Format (VCF) file.

    Raises:
        KeyError: If the `(sample, var_type)` composite index key is not present
            within the `SAMPLES` DataFrame.
    """
    return SAMPLES.loc[(wildcards.sample, wildcards.var_type.lower()), "vcf_path"]


def get_filter_expression(wildcards):
    """Dynamically construct a bcftools filter query string from workflow parameters.

    Evaluates global workflow configuration settings (`config`) to assemble a
    boolean filtering query compatible with `bcftools filter -i`. The constructed
    expression selectively includes high-confidence variants flagged as "PASS"
    and/or specific Variant Quality Score Recalibration (VQSR) sensitivity tranches
    tailored to the variant class (SNP vs. INDEL).

    Args:
        wildcards (snakemake.io.Wildcards): Snakemake wildcards object containing
            the attribute `var_type` (e.g., 'snp', 'snps', 'indel', or 'indels').

    Returns:
        str: A formatted logical query string combining PASS criteria and exact
        VQSR tranche filter conditions (e.g.,
        `FILTER="PASS" || FILTER="VQSRTrancheSNP99.90to99.95" || FILTER="VQSRTrancheSNP99.95to100.00"`).
    """
    var_type = "indel" if wildcards.var_type.lower() in ["indel", "indels"] else "snp"
    prefix = "INDEL" if var_type == "indel" else "SNP"

    tranches = config["variant_processing"]["filter_variants"]["tranches"].get(
        var_type, []
    )
    include_pass = config["variant_processing"]["filter_variants"].get(
        "include_pass", True
    )

    conditions = []
    if include_pass:
        conditions.append('FILTER="PASS"')

    for tranche in tranches:
        conditions.append(f'FILTER="VQSRTranche{prefix}{tranche}"')

    return " || ".join(conditions)
