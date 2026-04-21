#!/bin/bash
set -euo pipefail

#
# Strip a single gnomAD v4.1 joint VCF down to the INFO fields we need,
# rename _joint -> legacy (v4.0-style) names, bgzip + tabix the output.
#
# Usage:
#   strip_gnomad41_one.sh INPUT_VCF OUTPUT_VCF [RENAME_FILE]
#
# If RENAME_FILE is omitted, one is written next to OUTPUT_VCF.
#
# Requires bcftools + tabix on PATH.
#

INPUT_VCF="${1:?Usage: $0 INPUT_VCF OUTPUT_VCF [RENAME_FILE]}"
OUTPUT_VCF="${2:?Usage: $0 INPUT_VCF OUTPUT_VCF [RENAME_FILE]}"
RENAME_FILE="${3:-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=gnomad41_fields.sh
source "${SCRIPT_DIR}/gnomad41_fields.sh"

OUTPUT_DIR="$(dirname "$OUTPUT_VCF")"
mkdir -p "$OUTPUT_DIR"

if [ -z "$RENAME_FILE" ]; then
    RENAME_FILE="${OUTPUT_DIR}/rename_annots.txt"
    write_rename_file "$RENAME_FILE"
fi

KEEP_EXPR="$(build_keep_expr)"

echo "Processing $INPUT_VCF -> $OUTPUT_VCF"
echo "  keep expression: $KEEP_EXPR"
echo "  rename file:     $RENAME_FILE"

# Pipeline:
#   1. bcftools annotate --remove  : drop every INFO field except KEEP_FIELDS
#   2. bcftools annotate --rename-annots : apply _joint -> legacy renames
#      (split across two invocations because bcftools can't remove+rename in one pass
#       when the rename targets the kept fields)
#   3. sed : rewrite Number=A -> Number=1 in INFO header lines.
#      gnomAD VCFs are already split to one ALT per row, so A-typed fields are
#      effectively scalar; downstream tools read the renamed fields as Number=1.
#   4. bcftools view -O z : bgzip the output
bcftools annotate \
    --remove "${KEEP_EXPR}" \
    "${INPUT_VCF}" \
    | bcftools annotate \
        --rename-annots "${RENAME_FILE}" \
    | sed 's/,Number=A,/,Number=1,/' \
    | bcftools view -O z -o "${OUTPUT_VCF}"

tabix -f -p vcf "${OUTPUT_VCF}"

echo "Done: ${OUTPUT_VCF}"
