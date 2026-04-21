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

# Pass 1 — pure bcftools pipeline, streamed with BCF (-O u) between stages.
# Writes the bgzipped intermediate with header attributes untouched (AC/AF/AN
# stay declared Number=A here, which is VCF-spec-correct for reserved names).
#
#   1. bcftools norm --multiallelics -any : split any multi-allelic rows first so
#      every later stage operates on biallelic data. gnomAD v4.x is usually
#      pre-split, but this is a cheap safety net.
#   2. bcftools annotate --remove : drop every INFO field except KEEP_FIELDS
#   3. bcftools annotate --rename-annots : rename _joint -> legacy v4.0 names,
#      writing the final bgzipped intermediate directly via -O z.
TMP_VCF="${OUTPUT_VCF}.tmp.vcf.gz"
TMP_HDR="${OUTPUT_VCF}.tmp.hdr"
trap 'rm -f "$TMP_VCF" "$TMP_HDR"' EXIT

bcftools norm \
    --multiallelics -any \
    "${INPUT_VCF}" \
    -O u \
    | bcftools annotate \
        --remove "${KEEP_EXPR}" \
        -O u \
    | bcftools annotate \
        --rename-annots "${RENAME_FILE}" \
        -O z \
        -o "${TMP_VCF}"

# Pass 2 — rewrite Number=A -> Number=1 in the header for reserved names.
# Downstream parsers expect Number=1 to match v4.0 output. This has to happen
# after the last bcftools stage reads anything: newer bcftools (>=1.23) treats
# Number!=A on reserved names as a sanity failure when reading such a header.
# bcftools reheader swaps the header without running those sanity checks.
bcftools view -h "${TMP_VCF}" | sed 's/,Number=A,/,Number=1,/' > "${TMP_HDR}"
bcftools reheader -h "${TMP_HDR}" -o "${OUTPUT_VCF}" "${TMP_VCF}"

tabix -f -p vcf "${OUTPUT_VCF}"

echo "Done: ${OUTPUT_VCF}"
