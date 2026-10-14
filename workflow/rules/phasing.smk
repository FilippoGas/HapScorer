# workflow/rules/phasing.smk

# Injects fallback defaults if not defined in config.yaml
resources_cfg = config.setdefault("resources", {})
resources_cfg.setdefault("genetic_maps_tar", "resources/genetic_maps.b38.tar.gz")
resources_cfg.setdefault("genetic_maps_dir", "resources/maps")

# Autosomal chromosomes (1–22)
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)]

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
        "logs/unpack_genetic_map/{chr}.log",
    benchmark:
        "benchmarks/unpack_genetic_map/{chr}.tsv"
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
