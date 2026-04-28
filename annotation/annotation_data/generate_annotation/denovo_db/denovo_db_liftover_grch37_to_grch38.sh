#!/bin/bash

# Liftover the denovo-db hg19/GRCh37 VCF to GRCh38 using bcftools +liftover.
# Source: https://denovo-db.gs.washington.edu/denovo-db/Download.jsp
# denovo-db ships hg19 (UCSC chr-prefixed) only — we rename to Ensembl naming
# (1, 2, ..., MT) before lift to match VariantGrid's GRCh37/GRCh38 fastas + chain.
#
# See: claude/issue_1531_denovodb_plan.md

set -euo pipefail

if [[ -z "${VEP_ANNOTATION_DIR:-}" ]]; then
  echo "Please set 'VEP_ANNOTATION_DIR' to point at the annotation data root" >&2
  exit 1
fi

INPUT_GRCH37=${VEP_ANNOTATION_DIR}/annotation_data/GRCh37/denovo-db.variants.v.1.6.1.vcf.gz
OUTPUT_GRCH38=${VEP_ANNOTATION_DIR}/annotation_data/GRCh38/denovo-db.variants.v.1.6.1.GRCh38.vcf.gz
REJECTED_GRCH38=${VEP_ANNOTATION_DIR}/annotation_data/GRCh38/denovo-db.variants.v.1.6.1.GRCh38.rejected.vcf.gz

SRC_FASTA=${VEP_ANNOTATION_DIR}/annotation_data/fasta/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
TGT_FASTA=${VEP_ANNOTATION_DIR}/annotation_data/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
CHAIN=${VEP_ANNOTATION_DIR}/annotation_data/liftover/GRCh37_to_GRCh38.chain.gz

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
  exit 1
fi

# chr1 -> 1, ..., chrM -> MT  (Ensembl naming used by our GRCh37 / GRCh38 fastas + chain)
CHR_RENAME=$(mktemp)
trap 'rm -f "${CHR_RENAME}"' EXIT
{
  for n in $(seq 1 22) X Y; do
    echo -e "chr${n}\t${n}"
  done
  echo -e "chrM\tMT"
} > "${CHR_RENAME}"

mkdir -p "$(dirname "${OUTPUT_GRCH38}")"

echo "Renaming chromosomes (chr1 -> 1, chrM -> MT) and lifting GRCh37 -> GRCh38..."
bcftools annotate --rename-chrs "${CHR_RENAME}" "${INPUT_GRCH37}" -Ou \
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
if [[ "${LOST}" -gt 0 ]]; then
  echo "rejection reasons:"
  bcftools view "${REJECTED_GRCH38}" \
    | grep -oP 'REJECT_REASON=[^;[:space:]]*' \
    | sort | uniq -c | sort -rn
fi
echo "Output:   ${OUTPUT_GRCH38}"
