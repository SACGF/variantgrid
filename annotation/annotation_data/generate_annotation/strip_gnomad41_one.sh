#!/bin/bash
set -euo pipefail

#
# Strip a single gnomAD v4.1 joint VCF down to the INFO fields we need,
# rename _joint -> legacy (v4.0-style) names, rename contigs (e.g.
# NC_000001.11 -> 1) using a bcftools chrom map, bgzip + tabix the output.
#
# Usage:
#   strip_gnomad41_one.sh INPUT_VCF OUTPUT_VCF CHROM_MAP
#
# RENAME_FILE for the _joint -> v4.0 name TSV is auto-generated next to
# OUTPUT_VCF; override by exporting RENAME_FILE in the environment.
#
# Requires bcftools + tabix on PATH.
#

INPUT_VCF="${1:?Usage: $0 INPUT_VCF OUTPUT_VCF CHROM_MAP}"
OUTPUT_VCF="${2:?Usage: $0 INPUT_VCF OUTPUT_VCF CHROM_MAP}"
CHROM_MAP="${3:?Usage: $0 INPUT_VCF OUTPUT_VCF CHROM_MAP}"
RENAME_FILE="${RENAME_FILE:-}"
NON_PAR_HEADER_FILE="${NON_PAR_HEADER_FILE:-}"

if [ ! -f "$CHROM_MAP" ]; then
    echo "CHROM_MAP not found: $CHROM_MAP" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=gnomad41_fields.sh
source "${SCRIPT_DIR}/gnomad41_fields.sh"

OUTPUT_DIR="$(dirname "$OUTPUT_VCF")"
mkdir -p "$OUTPUT_DIR"

if [ -z "$RENAME_FILE" ]; then
    RENAME_FILE="${OUTPUT_DIR}/rename_annots.txt"
    write_rename_file "$RENAME_FILE"
fi

if [ -z "$NON_PAR_HEADER_FILE" ]; then
    NON_PAR_HEADER_FILE="${OUTPUT_DIR}/non_par.hdr"
    write_non_par_header_file "$NON_PAR_HEADER_FILE"
fi

KEEP_EXPR="$(build_keep_expr)"

# Detect sex-chrom inputs from the gnomAD source filename
# (gnomad.joint.v4.1.sites.chrX.vcf.bgz). Only chrX/chrY need the awk
# stamping pass; autosomes get only the non_par header declaration (free,
# via --header-lines on an existing streaming stage).
IS_SEX_CHROM=0
if [[ "$(basename "$INPUT_VCF")" =~ \.chr([XY])\. ]]; then
    IS_SEX_CHROM=1
fi

echo "Processing $INPUT_VCF -> $OUTPUT_VCF"
echo "  keep expression:  $KEEP_EXPR"
echo "  rename file:      $RENAME_FILE"
echo "  non_par header:   $NON_PAR_HEADER_FILE"
echo "  chrom map:        $CHROM_MAP"
echo "  sex chrom:        $IS_SEX_CHROM"

# Pass 1 — pure bcftools pipeline, streamed with BCF (-O u) between stages.
# Writes the bgzipped intermediate with header attributes untouched (AC/AF/AN
# stay declared Number=A here, which is VCF-spec-correct for reserved names).
#
#   1. bcftools norm --multiallelics -any : split any multi-allelic rows first so
#      every later stage operates on biallelic data. gnomAD v4.x is usually
#      pre-split, but this is a cheap safety net.
#   2. bcftools annotate --remove : drop every INFO field except KEEP_FIELDS
#   3. bcftools annotate --rename-annots + --rename-chrs + --header-lines :
#      rename _joint -> legacy v4.0 names AND rename contigs (e.g.
#      NC_000001.11 -> 1) AND inject the synthetic non_par INFO header
#      declaration. Putting --header-lines on this existing streaming stage
#      means autosomes pay zero extra cost while still advertising non_par
#      in their header (so bcftools concat downstream — which takes its
#      header from the first input file — keeps the flag for chrX/chrY).
#
# For chrX/chrY only, an additional decompress -> awk -> bgzip pass
# (stamp_non_par_from_witness) sets the non_par flag from witness FAF
# fields and strips the witnesses. v4.1 dropped the top-level non_par flag
# v4.0 had; gnomAD only populates faf99_joint_XX/_XY on non-PAR records,
# so their presence is the marker.
TMP_VCF="${OUTPUT_VCF}.tmp.vcf.gz"
TMP_HDR="${OUTPUT_VCF}.tmp.hdr"
trap 'rm -f "$TMP_VCF" "$TMP_HDR"' EXIT

if [ "$IS_SEX_CHROM" = "1" ]; then
    # chrX/chrY: stage 3 emits uncompressed VCF, awk stamps non_par,
    # bgzip writes TMP_VCF. One decomp/recomp per sex chrom (2 of 24).
    bcftools norm \
        --multiallelics -any \
        "${INPUT_VCF}" \
        -O u \
        | bcftools annotate \
            --remove "${KEEP_EXPR}" \
            -O u \
        | bcftools annotate \
            --rename-annots "${RENAME_FILE}" \
            --rename-chrs "${CHROM_MAP}" \
            --header-lines "${NON_PAR_HEADER_FILE}" \
            -O v \
        | stamp_non_par_from_witness \
        | bgzip -c > "${TMP_VCF}"
else
    # Autosomes: stage 3 writes bgzipped TMP_VCF directly. The non_par
    # INFO declaration lands in the header via --header-lines but no
    # autosomal record will ever have the flag set.
    bcftools norm \
        --multiallelics -any \
        "${INPUT_VCF}" \
        -O u \
        | bcftools annotate \
            --remove "${KEEP_EXPR}" \
            -O u \
        | bcftools annotate \
            --rename-annots "${RENAME_FILE}" \
            --rename-chrs "${CHROM_MAP}" \
            --header-lines "${NON_PAR_HEADER_FILE}" \
            -O z \
            -o "${TMP_VCF}"
fi

# Pass 2 — rewrite Number=A -> Number=1 in the header for reserved names.
# Downstream parsers expect Number=1 to match v4.0 output. This has to happen
# after the last bcftools stage reads anything: newer bcftools (>=1.23) treats
# Number!=A on reserved names as a sanity failure when reading such a header.
# bcftools reheader swaps the header without running those sanity checks.
bcftools view -h "${TMP_VCF}" | sed 's/,Number=A,/,Number=1,/' > "${TMP_HDR}"
bcftools reheader -h "${TMP_HDR}" -o "${OUTPUT_VCF}" "${TMP_VCF}"

tabix -f -p vcf "${OUTPUT_VCF}"

echo "Done: ${OUTPUT_VCF}"
