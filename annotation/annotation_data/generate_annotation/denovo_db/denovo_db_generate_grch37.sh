#!/bin/bash

# Build the GRCh37 denovo-db VCF used by VEP --custom annotation.
#
# Why this exists: the upstream denovo-db VCF lacks StudyName / PubmedID /
# PrimaryPhenotype as INFO fields — those only appear in the TSV. We rebuild
# the VCF from the TSV(s) so we can carry the clinically useful fields through
# to VEP.
#
# Inputs (download from https://denovo-db.gs.washington.edu/denovo-db/Download.jsp):
#   ${VEP_ANNOTATION_DIR}/annotation_data/GRCh37/denovo-db.ssc-samples.variants.v.1.6.1.tsv.gz
#   ${VEP_ANNOTATION_DIR}/annotation_data/GRCh37/denovo-db.non-ssc-samples.variants.v.1.6.1.tsv.gz
#
# Output:
#   ${VEP_ANNOTATION_DIR}/annotation_data/GRCh37/denovo-db.variants.v.1.6.1.vcf.gz

set -euo pipefail

if [[ -z "${VEP_ANNOTATION_DIR:-}" ]]; then
  echo "Please set 'VEP_ANNOTATION_DIR' to point at the annotation data root" >&2
  exit 1
fi

SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")
GRCH37_DIR=${VEP_ANNOTATION_DIR}/annotation_data/GRCh37
SSC_TSV=${GRCH37_DIR}/denovo-db.ssc-samples.variants.v.1.6.1.tsv.gz
NONSSC_TSV=${GRCH37_DIR}/denovo-db.non-ssc-samples.variants.v.1.6.1.tsv.gz
OUTPUT=${GRCH37_DIR}/denovo-db.variants.v.1.6.1.vcf.gz

for f in "${SSC_TSV}" "${NONSSC_TSV}"; do
  if [[ ! -e "${f}" ]]; then
    echo "Missing input TSV: ${f}" >&2
    echo "Download from https://denovo-db.gs.washington.edu/denovo-db/Download.jsp" >&2
    exit 1
  fi
done

echo "Converting denovo-db TSV(s) to VCF and sorting..."
python3 "${SCRIPT_DIR}/denovo_db_tsv_to_vcf.py" "${SSC_TSV}" "${NONSSC_TSV}" \
  | bcftools sort -Oz -o "${OUTPUT}"

bcftools index -t "${OUTPUT}"
echo "Wrote ${OUTPUT} ($(bcftools view -H "${OUTPUT}" | wc -l) variants)"
