#!/bin/bash

# Liftover the denovo-db hg19/GRCh37 VCF to GRCh38 using bcftools +liftover.
# Source: https://denovo-db.gs.washington.edu/denovo-db/Download.jsp
# denovo-db ships hg19 (UCSC chr-prefixed) only — we rename to Ensembl naming
# (1, 2, ..., MT) before lift to match VariantGrid's GRCh37/GRCh38 fastas + chain.
#
# See: claude/issue_1531_denovodb_plan.md

set -euo pipefail

if [[ -z "${ANNOTATION_DIR:-}" ]]; then
  echo "Please set 'ANNOTATION_DIR' to point at the annotation data root" >&2
  exit 1
fi

INPUT_GRCH37=${ANNOTATION_DIR}/VEP/annotation_data/GRCh37/denovo-db.variants.v.1.6.1.vcf.gz
OUTPUT_GRCH38=${ANNOTATION_DIR}/VEP/annotation_data/GRCh38/denovo-db.variants.v.1.6.1.GRCh38.vcf.gz
REJECTED_GRCH38=${ANNOTATION_DIR}/VEP/annotation_data/GRCh38/denovo-db.variants.v.1.6.1.GRCh38.rejected.vcf.gz

SRC_FASTA=${ANNOTATION_DIR}/fasta/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
TGT_FASTA=${ANNOTATION_DIR}/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
CHAIN=${ANNOTATION_DIR}/liftover/GRCh37_to_GRCh38.chain.gz

if [[ ! -e "${INPUT_GRCH37}" ]]; then
  echo "Input not found: ${INPUT_GRCH37}" >&2
  echo "Download denovo-db hg19 VCFs (SSC + non-SSC) from" >&2
  echo "  https://denovo-db.gs.washington.edu/denovo-db/Download.jsp" >&2
  echo "and merge into ${INPUT_GRCH37}" >&2
  exit 1
fi

for f in "${SRC_FASTA}" "${TGT_FASTA}" "${CHAIN}"; do
  if [[ ! -e "${f}" ]]; then
    echo "Missing reference file: ${f}" >&2
    exit 1
  fi
done

if ! bcftools plugin -l 2>/dev/null | grep -q '^liftover$'; then
  echo "bcftools +liftover plugin not found. Need bcftools >= 1.18 with plugins." >&2
  echo "If plugins are installed but not on the default search path, set e.g.:" >&2
  echo "  export BCFTOOLS_PLUGINS=/usr/share/bcftools/plugins" >&2
  exit 1
fi

# Strip 'chr' prefix from every contig in the input VCF (chrM -> MT specifically;
# everything else just drops the prefix). Built from the input header so unplaced
# contigs like chrGL000209.1 -> GL000209.1 are handled too. Ensembl naming is
# what our GRCh37 / GRCh38 fastas + chain use.
CHR_RENAME=$(mktemp)
trap 'rm -f "${CHR_RENAME}"' EXIT
bcftools view -h "${INPUT_GRCH37}" \
  | awk '/^##contig=<ID=/ { match($0, /ID=[^,>]+/); id = substr($0, RSTART+3, RLENGTH-3); print id }' \
  | awk 'BEGIN{OFS="\t"} { new = $1; if (new == "chrM") new = "MT"; else sub(/^chr/, "", new); print $1, new }' \
  > "${CHR_RENAME}"

mkdir -p "$(dirname "${OUTPUT_GRCH38}")"

echo "Renaming chromosomes (chr1 -> 1, chrM -> MT) and lifting GRCh37 -> GRCh38..."
# denovo-db has a handful of records whose REF doesn't match the GRCh37 fasta
# (e.g. 1:1720573 REF=AAGA vs. reference AGAA). bcftools +liftover errors out on
# the first such mismatch, killing the whole run, so drop them up front with
# `bcftools norm --check-ref x` (excludes bad sites; count is reported to stderr).
bcftools annotate --rename-chrs "${CHR_RENAME}" "${INPUT_GRCH37}" -Ou \
  | bcftools norm --check-ref x --fasta-ref "${SRC_FASTA}" -Ou \
  | bcftools sort -Ou \
  | bcftools +liftover -Ou -- \
      -s "${SRC_FASTA}" \
      -f "${TGT_FASTA}" \
      -c "${CHAIN}" \
      --reject "${REJECTED_GRCH38}" \
      --reject-type z \
      --write-src \
  | bcftools sort -Oz -o "${OUTPUT_GRCH38}"

bcftools index -t "${OUTPUT_GRCH38}"
[[ -e "${REJECTED_GRCH38}" ]] && bcftools index -t "${REJECTED_GRCH38}" || true

echo
echo "===== Liftover summary ====="
KEPT=$(bcftools view -H "${OUTPUT_GRCH38}" | wc -l)
LOST=0
[[ -e "${REJECTED_GRCH38}" ]] && LOST=$(bcftools view -H "${REJECTED_GRCH38}" | wc -l)
echo "kept:     ${KEPT}"
echo "rejected: ${LOST}"
echo "Output:   ${OUTPUT_GRCH38}"
echo "Rejected: ${REJECTED_GRCH38}"
