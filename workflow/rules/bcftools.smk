# workflow/rules/bfctools.smk


wildcard_constraints:
    var_type="SNP|INDEL",


rule filter_variants:
    """
    Filter SNPs or INDELs for a given sample based on a VQSR quality threshold.

    Applies ``bcftools filter`` to isolate high-quality variants.

    :input vcf: Raw compressed VCF file (``data/{sample}_{var_type}.vcf.gz``).
    :output vcf: Filtered VCF file.
    :output csi: Index file for the filtered VCF.
    :param threshold: Minimum VQSR quality score required to pass filtering.
    """
    input:
        vcf="data/{sample}_{var_type}.vcf.gz",
    output:
        vcf=temp("results/filtered/{sample}_{var_type}_filtered.vcf.gz"),
        csi=temp("results/filtered/{sample}_{var_type}_filtered.vcf.gz.csi"),
    log:
        "logs/filter_variants/{sample}_{var_type}.log",
    benchmark:
        "benchmarks/filter_variants/{sample}_{var_type}.tsv"
    conda:
        "../envs/bcftools.yaml"
    threads: 2
    resources:
        mem_mb=4000,
        runtime=60,
    params:
        threshold=config["bcftools"]["vqsr_threshold"],
    message:
        "Filtering {wildcards.var_type} variants for {wildcards.sample} (VQSR > {params.threshold})"
    shell:
        """
        bcftools filter --threads {threads} -i 'VQSR > {params.threshold}' {input.vcf} -Oz -o {output.vcf} >{log} 2>&1
        bcftools index --threads {threads} {output.vcf} >>{log} 2>&1
        """


rule concat_sample_variants:
    """
    Concatenate filtered SNP and INDEL VCF files into a single per-sample VCF.

    Uses ``bcftools concat`` with the ``-a`` flag to allow overlap resolution.

    :input snps: Filtered SNP VCF file.
    :input indels: Filtered INDEL VCF file.
    :input snps_csi: CSI index for filtered SNPs.
    :input indels_csi: CSI index for filtered INDELs.
    :output vcf: Combined per-sample VCF file.
    :output csi: Index file for the combined VCF.
    """
    input:
        snps="results/filtered/{sample}_SNP_filtered.vcf.gz",
        indels="results/filtered/{sample}_INDEL_filtered.vcf.gz",
        snps_csi="results/filtered/{sample}_SNP_filtered.vcf.gz.csi",
        indels_csi="results/filtered/{sample}_INDEL_filtered.vcf.gz.csi",
    output:
        vcf="results/concat/{sample}.vcf.gz",
        csi="results/concat/{sample}.vcf.gz.csi",
    log:
        "logs/concat_sample_variants/{sample}.log",
    benchmark:
        "benchmarks/concat_sample_variants/{sample}.tsv"
    conda:
        "../envs/bcftools.yaml"
    threads: 2
    resources:
        mem_mb=4000,
        runtime=60,
    message:
        "Concatenating SNP and INDEL variants for sample {wildcards.sample}"
    shell:
        """
        bcftools concat --threads {threads} -a {input.snps} {input.indels} -Oz -o {output.vcf} >{log} 2>&1
        bcftools index --threads {threads} {output.vcf} >>{log} 2>&1
        """


rule merge_all_samples:
    """
    Merge concatenated VCF files across all cohort samples into a multi-sample VCF.

    Executes ``bcftools merge`` on all sample-level files defined in the workflow context.

    :input vcfs: List of per-sample concatenated VCF files.
    :input csis: List of index files for each per-sample VCF.
    :output vcf: Final unified multi-sample VCF file (``all_samples.vcf.gz``).
    """
    input:
        vcfs=expand("results/concat/{sample}.vcf.gz", sample=SAMPLES),
        csis=expand("results/concat/{sample}.vcf.gz.csi", sample=SAMPLES),
    output:
        vcf=protected("results/merged/all_samples.vcf.gz"),
    log:
        "logs/merge_all_samples/all_samples.log",
    benchmark:
        "benchmarks/merge_all_samples/all_samples.tsv"
    conda:
        "../envs/bcftools.yaml"
    threads: 4
    resources:
        mem_mb=8000,
        runtime=120,
    message:
        "Merging per-sample VCFs into cohort VCF across all samples"
    shell:
        """
        bcftools merge {input.vcfs} -Oz -o {output.vcf} >{log} 2>&1
        """
