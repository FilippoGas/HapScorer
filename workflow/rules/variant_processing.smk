# workflow/rules/variant_processing.smk


rule concat_sample_variants:
    """
    Concatenate pre-filtered SNP and INDEL VCF files into a single per-sample VCF.

    Uses ``bcftools concat`` with the ``-a`` flag to resolve overlapping positions
    and produce a unified per-sample variant call set.

    :input snps: Pre-filtered SNP VCF file retrieved from sample metadata.
    :input indels: Pre-filtered INDEL VCF file retrieved from sample metadata.
    :output vcf: Combined per-sample compressed VCF file.
    :output csi: Index file for the combined per-sample VCF.
    """
    input:
        snps=lambda wildcards: SAMPLES.loc[(wildcards.sample, "snp"), "vcf_path"],
        indels=lambda wildcards: SAMPLES.loc[(wildcards.sample, "indel"), "vcf_path"],
    output:
        vcf="results/concat/{sample}.vcf.gz",
        csi="results/concat/{sample}.vcf.gz.csi",
    log:
        "logs/concat_sample_variants/{sample}.log",
    benchmark:
        "benchmarks/concat_sample_variants/{sample}.tsv"
    shadow:
        "minimal"
    conda:
        "../envs/variant_processing.yaml"
    threads: 2
    resources:
        mem_mb=2000,
    message:
        "Concatenating SNP and INDEL variants for sample {wildcards.sample}"
    shell:
        """
        bcftools index --threads {threads} {input.snps} -o {input.snps}.csi >{log} 2>&1
        bcftools index --threads {threads} {input.indels} -o {input.indels}.csi >>{log} 2>&1
        bcftools concat --threads {threads} -a {input.snps} {input.indels} -Oz -o {output.vcf} >>{log} 2>&1
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
