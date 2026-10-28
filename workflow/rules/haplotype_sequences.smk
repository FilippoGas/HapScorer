# workflow/rules/haplotype_sequences.smk


rule extract_haplotypes_and_genotypes:
    """
    Extract unique haplotype combinations and generate sample genotype matrices.

    Parses phased VCF files and transcript annotations using pysam to extract
    exonic variants, identify unique haplotypic combinations for each transcript,
    and construct sample genotype mappings across all individual samples.

    :input phased_vcf: Phased VCF file containing genomic variants for a specific chromosome.
    :input annotations: CSV file containing exon-level transcript structural annotations.
    :output haplotypes_csv: Output CSV file detailing unique transcript haplotypes and variant signatures.
    :output genotypes_csv: Output CSV matrix mapping sample genotypes to haplotype identifiers.

    """
    input:
        phased_vcf="results/phasing/phased/{chr}.phased.vcf.gz",
        annotations="resources/reference_data/protein_coding_transcripts.csv",
    output:
        haplotypes_csv="results/haplotype_sequences/haplotypes/{chr}_haplotypes.csv",
        genotypes_csv="results/haplotype_sequences/genotypes/{chr}_genotypes.csv",
    log:
        "logs/haplotype_sequences/extract_haplotypes_and_genotypes/{chr}.log",
    benchmark:
        "benchmarks/haplotype_sequences/extract_haplotypes_and_genotypes/{chr}.tsv"
    conda:
        "../envs/haplotype_sequences.yaml"
    threads: 1
    resources:
        mem_mb=7500,
        runtime=30,
    message:
        "Extracting haplotypes and genotypes for chromosome {wildcards.chr}"
    script:
        "../scripts/extract_haplotypes_and_genotypes.py"


# rule extract_haplotypes:
#     input:
#         phased_vcf="results/phasing/phased/{chr}.phased.vcf.gz",
#         annotations="resources/reference_data/protein_coding_transcripts.csv",
#         samples="results/phasing/phased/samples.txt",
#     output:
#         haplotypes="results/haplotypes/{chr}_haplotypes.csv",
#     conda:
#         "../envs/Renv.yaml"
#     threads: 4
#     script:
#         "../scripts/extract_haplotype.R"

# rule genotypes:
#     input:
#         phased_vcf="results/phasing/phased/{chr}.phased.vcf.gz",
#         annotations="resources/reference_data/protein_coding_transcripts.csv",
#         samples="results/phasing/phased/samples.txt",
#         haplotypes="results/haplotypes/{chr}_haplotypes.csv",
#     output:
#         "results/genotypes/{chr}_genotypes.csv",
#     conda:
#         "../envs/Renv.yaml"
#     threads: 4
#     script:
#         "../scripts/genotype.R"
