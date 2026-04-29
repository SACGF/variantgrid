#!/bin/bash

# Build the GRCh37 denovo-db VCF used by VEP --custom annotation.
#
# Why this exists: the upstream denovo-db VCF lacks StudyName / PubmedID /
# PrimaryPhenotype as INFO fields — those only appear in the TSV. We rebuild
# the VCF from the TSV(s) so we can carry the clinically useful fields through
# to VEP.
#
# Inputs (download from https://denovo-db.gs.washington.edu/denovo-db/Download.jsp):
#   ${VEP_ANNOTATION_DIR}/annotation_data/GRCh37/denovo-db.ssc-samples.variants.tsv.gz
#   ${VEP_ANNOTATION_DIR}/annotation_data/GRCh37/denovo-db.non-ssc-samples.variants.tsv.gz
#
# Output:
#   ${VEP_ANNOTATION_DIR}/annotation_data/GRCh37/denovo-db.variants.v.<VERSION>.vcf.gz
# (VERSION is read from the '##version=denovo-db.v.<VERSION>' line at the top of each TSV.)

set -euo pipefail

if [[ -z "${VEP_ANNOTATION_DIR:-}" ]]; then
  echo "Please set 'VEP_ANNOTATION_DIR' to point at the annotation data root" >&2
  exit 1
fi

SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")
GRCH37_DIR=${VEP_ANNOTATION_DIR}/annotation_data/GRCh37
SSC_TSV=${GRCH37_DIR}/denovo-db.ssc-samples.variants.tsv.gz
NONSSC_TSV=${GRCH37_DIR}/denovo-db.non-ssc-samples.variants.tsv.gz

for f in "${SSC_TSV}" "${NONSSC_TSV}"; do
  if [[ ! -e "${f}" ]]; then
    echo "Missing input TSV: ${f}" >&2
    echo "Download from https://denovo-db.gs.washington.edu/denovo-db/Download.jsp" >&2
    exit 1
  fi
done

# Pull the version from the '##version=denovo-db.v.<VERSION>' line at the top of the TSV.
# Use process substitution + a single `read` so we don't SIGPIPE zcat (which would
# trip `set -o pipefail`).
read -r FIRST_LINE < <(zcat "${SSC_TSV}")
VERSION=$(echo "${FIRST_LINE}" | sed -n 's/^##version=denovo-db\.v\.\(\S\+\).*$/\1/p')
if [[ -z "${VERSION}" ]]; then
  echo "Could not parse '##version=denovo-db.v.<VERSION>' from first line of ${SSC_TSV}" >&2
  exit 1
fi
echo "Detected denovo-db version: ${VERSION}"
OUTPUT=${GRCH37_DIR}/denovo-db.variants.v.${VERSION}.vcf.gz

echo "Converting denovo-db TSV(s) to VCF and sorting..."
python3 "${SCRIPT_DIR}/denovo_db_tsv_to_vcf.py" "${SSC_TSV}" "${NONSSC_TSV}" \
  | bcftools sort -Oz -o "${OUTPUT}"

bcftools index -t "${OUTPUT}"
echo "Wrote ${OUTPUT} ($(bcftools view -H "${OUTPUT}" | wc -l) variants)"
