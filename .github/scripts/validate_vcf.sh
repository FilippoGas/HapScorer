#!/usr/bin/env bash
set -euo pipefail

# Ensure mandatory arguments are provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <result_vcf> <expected_vcf> [analysis_prefix]" >&2
    exit 1
fi

RESULT_VCF="$1"
EXPECTED_VCF="$2"
PREFIX="${3:-vcf_comparison}"
ISEC_DIR="isec_${PREFIX}"

echo "[$(date +'%Y-%m-%d %H:%M:%S')] INFO: Initiating variant concordance evaluation for ${PREFIX}..."

# Verify input file existence
if [[ ! -f "$RESULT_VCF" ]]; then
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] ERROR: Result VCF not found at ${RESULT_VCF}" >&2
    exit 1
fi

if [[ ! -f "$EXPECTED_VCF" ]]; then
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] ERROR: Expected reference VCF not found at ${EXPECTED_VCF}" >&2
    exit 1
fi

# Perform VCF intersection
# isec output: 0000.vcf (records unique to RESULT_VCF), 0001.vcf (records unique to EXPECTED_VCF)
bcftools isec -p "$ISEC_DIR" "$RESULT_VCF" "$EXPECTED_VCF"

DIFF_RESULT=$(bcftools view -H "${ISEC_DIR}/0000.vcf" | wc -l)
DIFF_EXPECTED=$(bcftools view -H "${ISEC_DIR}/0001.vcf" | wc -l)

# Assert zero dataset-exclusive variants
if [ "$DIFF_RESULT" -ne 0 ] || [ "$DIFF_EXPECTED" -ne 0 ]; then
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] ERROR: Concordance validation failed for ${PREFIX}." >&2
    echo "Dataset-exclusive variants detected -> Result: ${DIFF_RESULT} | Expected: ${DIFF_EXPECTED}" >&2
    exit 1
fi

rm -rf "$ISEC_DIR"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] INFO: Concordance validated successfully. No exclusive variants detected for ${PREFIX}."
