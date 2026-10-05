# workflow/rules/variant_processing.smk


rule filter_variants:
    """
    Filter SNPs or INDELs for a given sample based on config-defined VQSR tranches.

    Applies ``bcftools filter`` to isolate high-confidence variants using logical
    expressions generated dynamically from workflow configuration settings.

    :input vcf: Input VCF file retrieved via ``get_raw_vcf`` from sample metadata.
    :output vcf: Filtered compressed VCF file.
    :output csi: Index file for the filtered VCF.
    """
    input:
        vcf=get_raw_vcf,
    output:
        vcf=temp("results/filtered/{sample}_{var_type}_filtered.vcf.gz"),
        csi=temp("results/filtered/{sample}_{var_type}_filtered.vcf.gz.csi"),
    log:
        "logs/filter_variants/{sample}_{var_type}.log",
    benchmark:
        "benchmarks/filter_variants/{sample}_{var_type}.tsv"
    conda:
        "../envs/variant_processing.yaml"
    threads: 2
    resources:
        mem_mb=2000,
    params:
        filter_expr=get_filter_expression,
    message:
        "Filtering {wildcards.var_type} variants for {wildcards.sample} with filter '{params.filter_expr}'"
    shell:
        """
        bcftools filter --threads {threads} -i '{params.filter_expr}' {input.vcf} -Oz -o {output.vcf} >{log} 2>&1
        bcftools index --threads {threads} {output.vcf} >>{log} 2>&1
        """


rule concat_sample_variants:
    """
    Concatenate filtered SNP and INDEL VCF files into a single per-sample VCF.

    Uses ``bcftools concat`` with the ``-a`` flag to resolve overlapping positions
    and produce a unified per-sample variant call set.

    :input snps: Filtered SNP VCF file.
    :input indels: Filtered INDEL VCF file.
    :input snps_csi: CSI index for filtered SNPs.
    :input indels_csi: CSI index for filtered INDELs.
    :output vcf: Combined per-sample compressed VCF file.
    :output csi: Index file for the combined per-sample VCF.
    """
    input:
        snps="results/filtered/{sample}_snp_filtered.vcf.gz",
        indels="results/filtered/{sample}_indel_filtered.vcf.gz",
        snps_csi="results/filtered/{sample}_snp_filtered.vcf.gz.csi",
        indels_csi="results/filtered/{sample}_indel_filtered.vcf.gz.csi",
    output:
        vcf="results/concat/{sample}.vcf.gz",
        csi="results/concat/{sample}.vcf.gz.csi",
    log:
        "logs/concat_sample_variants/{sample}.log",
    benchmark:
        "benchmarks/concat_sample_variants/{sample}.tsv"
    conda:
        "../envs/variant_processing.yaml"
    threads: 2
    resources:
        mem_mb=2000,
    message:
        "Concatenating SNP and INDEL variants for sample {wildcards.sample}"
    shell:
        """
        bcftools concat --threads {threads} -a {input.snps} {input.indels} -Oz -o {output.vcf} >{log} 2>&1
        bcftools index --threads {threads} {output.vcf} >>{log} 2>&1
        """


rule merge_and_sort_samples:
    """
    Merge concatenated VCF files across all cohort samples into a multi-sample VCF,
    sort coordinates, and build the index required for chromosome splitting.

    Executes ``bcftools merge`` piped into ``bcftools sort``, followed by indexing.

    :input vcfs: List of per-sample concatenated VCF files.
    :input csis: List of index files for each per-sample VCF.
    :output vcf: Final unified multi-sample VCF file (``all_samples.vcf.gz``).
    :output csi: Index file for the merged multi-sample VCF.
    """
    input:
        vcfs=expand("results/concat/{sample}.vcf.gz", sample=SAMPLES["sample"].unique()),
        csis=expand(
            "results/concat/{sample}.vcf.gz.csi", sample=SAMPLES["sample"].unique()
        ),
    output:
        vcf=protected("results/merged/all_samples.vcf.gz"),
        csi=protected("results/merged/all_samples.vcf.gz.csi"),
    log:
        "logs/merge_and_sort_samples/all_samples.log",
    benchmark:
        "benchmarks/merge_and_sort_samples/all_samples.tsv"
    conda:
        "../envs/variant_processing.yaml"
    threads: 4
    resources:
        mem_mb=4000,
    message:
        "Merging, sorting, and indexing cohort VCF across all samples"
    shell:
        """
        bcftools merge --threads {threads} --force-single {input.vcfs} -Ou \
            | bcftools sort -T {resources.tmpdir}/ -Oz -o {output.vcf} >{log} 2>&1
        bcftools index --threads {threads} {output.vcf} >>{log} 2>&1
        """
