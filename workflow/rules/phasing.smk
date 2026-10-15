# workflow/rules/phasing.smk

# Injects fallback defaults if not defined in config.yaml
resources_cfg = config.setdefault("resources", {})
resources_cfg.setdefault("genetic_maps_tar", "resources/genetic_maps.b38.tar.gz")
resources_cfg.setdefault("genetic_maps_dir", "resources/maps")

phasing_cfg = config.setdefault("phasing", {})
phase_common_variants_cfg = phasing_cfg.setdefault("phase_common_variants", {})

# Autosomal chromosomes (1–22)
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)]

FILTER_MAF = phase_common_variants_cfg.setdefault("filter_maf", 0.001)

MAPS_TAR = Path(config["resources"]["genetic_maps_tar"])
MAPS_DIR = Path(config["resources"]["genetic_maps_dir"])


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
        vcf="results/variant_processing/merged/all_samples.vcf.gz",
        csi="results/variant_processing/merged/all_samples.vcf.gz.csi",
    output:
        vcf=temp("results/phasing/split/all_samples_{chr}.vcf.gz"),
        csi=temp("results/phasing/split/all_samples_{chr}.vcf.gz.csi"),
    log:
        "logs/phasing/split_vcf_by_chr/{chr}.log",
    benchmark:
        "benchmarks/phasing/split_vcf_by_chr/{chr}.tsv"
    conda:
        "../envs/phasing.yaml"
    threads: 2
    resources:
        mem_mb=2000,
    message:
        "Extracting and indexing chromosome {wildcards.chr} from merged cohort VCF"
    shell:
        """
        bcftools view --threads {threads} -r {wildcards.chr} {input.vcf} -Oz -o {output.vcf} >{log} 2>&1
        bcftools index --threads {threads} {output.vcf} >>{log} 2>&1
        """


rule unpack_genetic_map:
    """
    Extract the per-chromosome genetic map file from the bundled tarball resource.

    Executes ``tar`` to decompress and extract a specific chromosome's genetic map
    file directly into the configured resources directory.

    :input tar: Bundled tarball containing genetic map files for all chromosomes.
    :output gmap: Extracted per-chromosome genetic map file (``{chr}.b38.gmap.gz``).
    """
    input:
        tar=MAPS_TAR,
    output:
        gmap=MAPS_DIR / "{chr}.b38.gmap.gz",
    log:
        "logs/phasing/unpack_genetic_map/{chr}.log",
    benchmark:
        "benchmarks/phasing/unpack_genetic_map/{chr}.tsv"
    conda:
        "../envs/phasing.yaml"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=15,
    params:
        maps_dir=MAPS_DIR,
    message:
        "Extracting genetic map for chromosome {wildcards.chr}"
    shell:
        """
        mkdir -p {params.maps_dir}
        tar -xzf {input.tar} -C {params.maps_dir} {wildcards.chr}.b38.gmap.gz >{log} 2>&1
        """


rule phase_common_variants:
    """
    Perform common variant haplotype phasing on a per-chromosome basis using SHAPEIT5.

    Executes ``SHAPEIT5_phase_common`` on an individual chromosome VCF using
    a chromosome-specific genetic map, followed by indexing of the output VCF.

    :input vcf: Unphased per-chromosome VCF file.
    :input csi: CSI index for the unphased per-chromosome VCF.
    :input gmap: Recombination genetic map file for the target chromosome.
    :output bcf: Phased multi-sample BCF file for the target chromosome.
    :output csi: CSI index file for the phased BCF.
    :param region: Chromosome region target identifier passed to SHAPEIT5.
    :param filter_maf: Minor allele frequency threshold for filtering target variants.

    """
    input:
        vcf="results/phasing/split/all_samples_{chr}.vcf.gz",
        csi="results/phasing/split/all_samples_{chr}.vcf.gz.csi",
        gmap=MAPS_DIR / "{chr}.b38.gmap.gz",
    output:
        bcf=protected("results/phasing/phased/phased_{chr}.bcf"),
        csi=protected("results/phasing/phased/phased_{chr}.bcf.csi"),
    log:
        "logs/phasing/phase_common_variants/{chr}.log",
    benchmark:
        "benchmarks/phasing/phase_common_variants/{chr}.tsv"
    conda:
        "../envs/phasing.yaml"
    threads: 4
    resources:
        mem_mb=8000,
        runtime=120,
    params:
        region="{chr}",
        filter_maf=FILTER_MAF,
    message:
        "Running SHAPEIT5_phase_common for chromosome {wildcards.chr}"
    shell:
        """
        SHAPEIT5_phase_common \
            --input {input.vcf} \
            --map {input.gmap} \
            --region {params.region} \
            --filter-maf {params.filter_maf} \
            --output {output.bcf} \
            --thread {threads} >{log} 2>&1
        """
