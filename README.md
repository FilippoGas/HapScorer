<!-- README.md -->

# Transcript Haplotype & Sequence Processing Pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.27.1-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/FilippoGas/HapScorer/workflows/Tests/badge.svg?branch=main)](https://github.com/FilippoGas/HapScorer/actions?query=branch%3Amain+workflow%3ATests)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
<!-- [![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/FilippoGas/HapScorer) -->

This Snakemake workflow implements an end-to-end, reproducible computational framework to extract phased transcript haplotypes and construct mutated protein sequences. These outputs are specifically designed to evaluate the functional effects of genetic variants using Protein Language Models (PLMs), expanding upon the framework introduced in **HapScoreDB** [[1]](#references).

By evaluating co-occurring coding variants in their complete sequence context, the pipeline models intramolecular interactions and provides the foundational data to predict protein-level fitness across population-scale datasets or custom genomic case studies.

## Table of Contents
- [Workflow Architecture](#workflow-architecture)
- [Pipeline Deliverables](#pipeline-deliverables)
- [Deployment & Usage](#deployment--usage)
- [Configuration](#configuration)
- [Continuous Integration](#continuous-integration)
- [Authors](#authors)
- [Acknowledgments](#acknowledgments)
- [References](#references)

## Workflow Architecture

The pipeline automates the processing of cohort genomic variants into translated, haplotype-resolved protein-coding sequences and metadata through five core modules:

1. **Variant Processing & Merging:** Concatenates per-sample SNPs and INDELs using `bcftools` to build a unified multi-sample VCF.
2. **Chromosomal Phasing:** Splits variants per chromosome and phases common variants with **SHAPEIT5**.
3. **Reference Data Management:** Downloads and parses Ensembl bulk reference files (GTF, CDS, Pep, cDNA) with automatic checksum validation and exponential backoff retry mechanisms.
4. **Haplotype & Sequence Reconstruction:** Identifies unique exonic variant combinations per transcript, applies 5' UTR start-gains and 3' UTR stop-loss extensions, and translates nucleotide sequences into protein sequences using Biopython according to NCBI genetic code tables.
5. **Frequency Computation & Deliverable Separation:** Aggregates cohort allele occurrences using high-performance `polars` multithreading and exports disaggregated lightweight metadata and dedicated sequence mapping tables.

## Pipeline Deliverables

Outputs are generated per chromosome and stored in the `results/haplotype_sequences/` directory. The data is deliberately disaggregated to separate metadata from heavy sequence strings:

| Output File Path | Description |
| :--- | :--- |
| `haplotypes/{chr}_haplotypes_annotations.csv` | Lightweight protected CSV containing transcript/gene IDs, variant signatures, biological flags (`start_gain`, `multiple_starts_in_5`, `stop_loss_3_translated`), and populational frequency percentages (`freq_perc`). |
| `haplotypes/{chr}_haplotypes_sequences.csv` | Dedicated sequence mapping CSV connecting `haplotype_transcript_id` and `ensembl_transcript_id` to nucleotide (`nn_sequence`) and amino acid (`aa_sequence`) sequences. |
| `genotypes/{chr}_genotypes.csv` | Matrix mapping cohort sample identifiers to phased haplotype allele pairs per transcript (e.g., `ENSG...0/ENSG...1`). |

## Deployment & Usage

### Prerequisites & Requirements
- **Snakemake:** `≥ 9.27.1`
- **Environment Management:** Conda / Mamba is required for software isolation (`--sdm conda`).
- **Core Dependencies:** Python `≥ 3.14.8`.

Clone the repository and navigate into the workflow root directory:
```bash
git clone https://github.com/FilippoGas/HapScorer.git
cd HapScorer
```

### 1. Dry Run

Verify pipeline rule definitions, configuration, and dependency graphs prior to execution:

```bash
snakemake --dry-run
```

### 2. Test Suite Execution

Validate the installation and test the execution graph using pre-packaged, miniaturized test datasets:

```bash
snakemake --cores 4 --sdm conda --directory .test
```

### 3. Production Execution

Execute the complete pipeline locally across all cohort chromosomes using isolated Conda environments:

```bash
snakemake --cores 16 --sdm conda
```

## Configuration

All workflow variables, input file paths, and resource limits are controlled via the `config/` directory.

### Sample Table ([`config/samples.tsv`](config/samples.tsv))

Input VCF paths and variant classifications must be declared in a tab-separated format:

| Column | Type | Description |
| --- | --- | --- |
| `sample` | String | Unique sample identifier (alphanumeric, `_`, `-`). |
| `var_type` | String | Variant classification (`snp` or `indel`, case-insensitive). |
| `vcf_path` | String | Path to raw per-sample VCF file. |

### Global Parameters ([`config/config.yaml`](config/config.yaml))

Key execution controls include:

* `ncbi_genetic_code`: NCBI Genetic Code Table ID used for translation (default: `1`).
* `reference_data.release`: Ensembl release version (default: `"113"`).
* `reference_data.download_ensembl_dataset.max_retries`: Maximum network retry attempts (default: `3`).
* `phasing.phase_common_variants.filter_maf`: MAF threshold for SHAPEIT5 phasing (default: `0.01`).

For further configuration options, see the dedicated Wiki configuration page [here](). <!-- TODO -->

## Continuous Integration

The biological sanity and schema integrity of the generated deliverables can be strictly validated offline. This is integrated into GitHub Actions but can be run locally via:

```bash
python .github/scripts/validate_haplotypes.py <path_to_working_directory>
```

This validation script verifies IUPAC DNA base compliance, protein translation Methionine initiation rates ($\ge 90\%$), cohort frequency sums ($\sim 100\%$), and cross-table referential integrity.

## Authors

* **[Alex Callegaro](https://github.com/ale-0x)** (*Lead Developer & Maintainer - Snakemake Architecture, Python Re-engineering, Optimization & Documentation*) - University of Trento - [ORCID]()
* **[Filippo Gastaldello](https://github.com/FilippoGas)** (*Code Review & Supervision*) - CIBIO, University of Trento & COSBI - [ORCID]()

## Acknowledgments

This workflow builds upon and re-engineers the original R prototype developed by **[Filippo Gastaldello](https://github.com/FilippoGas)** (CIBIO, University of Trento & COSBI) for the *HapScoreDB* study [[1]](#references). Special thanks to him for providing the initial methodological reference and for reviewing the Pull Requests throughout the development of this repository.

## References

> **[1]** Mazza, F., Gastaldello, F., Dalfovo, D., Lattanzi, G., Romanel, A. “HapScoreDB: a database of protein language model functional scores for haplotype-resolved protein sequences”. In: _Nucleic Acids Research_ 54.D1 (Jan. 6, 2026), pp. D1087–D1097. DOI:
[10.1093/nar/gkaf1184](https://doi.org/10.1093/nar/gkaf1184)

> **[2]** Mölder F., Jablonski K. P., Letcher B., Hall M. B., van Dyken P. C., Tomkins-Tinch C. H., Sochat V., Forster J., Vieira F. G., Meesters C., Lee S., Twardziok S. O., Kanitz A., VanCampen J., Malladi V., Wilm A., Holtgrewe M., Rahmann S., Nahnsen S., Köster J. “Sustainable data analysis with Snakemake”. In: _F1000Research_ 10, 33 (Jan. 18, 2021).
DOI: [10.12688/f1000research.29032.3](https://doi.org/10.12688/f1000research.29032.3)
