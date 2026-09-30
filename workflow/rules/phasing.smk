# workflow/rules/phasing.smk

# Autosomal chromosomes (1–22)
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)]


rule split_vcf_by_chr:
    """
    Extract a single target chromosome from the merged cohort VCF and build its index.

    Applies ``bcftools view -r`` to isolate variants for a specific chromosome,
    followed by ``bcftools index`` to generate spatial indexing required by SHAPEIT5.

    :input vcf: Merged multi-sample VCF file (``results/merged/all_samples.vcf.gz``).
    :input csi: CSI index for the merged multi-sample VCF.
    :output vcf: Per-chromosome VCF file.
    :output csi: CSI index file for the per-chromosome VCF.
    """
    input:
        vcf="results/merged/all_samples.vcf.gz",
        csi="results/merged/all_samples.vcf.gz.csi",
    output:
        vcf=temp("results/phasing/split/all_samples_{chr}.vcf.gz"),
        csi=temp("results/phasing/split/all_samples_{chr}.vcf.gz.csi"),
    log:
        "logs/split_vcf_by_chr/{chr}.log",
    benchmark:
        "benchmarks/split_vcf_by_chr/{chr}.tsv"
    conda:
        "../envs/phasing.yaml"
    message:
        "Extracting and indexing chromosome {wildcards.chr} from merged cohort VCF"
    shell:
        """
        bcftools view --threads {threads} -r {wildcards.chr} {input.vcf} -Oz -o {output.vcf} >{log} 2>&1
        bcftools index --threads {threads} {output.vcf} >>{log} 2>&1
        """
