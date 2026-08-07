#!/bin/bash

# GRCh38 common/uncommon pre-process VCF as the INTERSECTION of gnomAD 4.0 and 4.1 AF>5
# so a VCF's 'common' partition is valid (skippable) when filtering rare under either version.
# @see https://github.com/SACGF/variantgrid/issues/1582

set -e

if [[ -z ${VEP_ANNOTATION_DIR} ]]; then
  echo "Please set 'VEP_ANNOTATION_DIR' variable to point where gnomAD files are" >&1
  exit 1
fi

GRCH38_DIR=${VEP_ANNOTATION_DIR}/annotation_data/GRCh38
OUTPUT=${GRCH38_DIR}/gnomad4.0_and_4.1_GRCh38_af_greater_than_5.intersection.stripped.vcf.gz

# -n=2 keeps sites present in both files; -w1 writes the records from the first (4.0) file
bcftools isec -n=2 -w1 --output-type z --output ${OUTPUT} \
  ${GRCH38_DIR}/gnomad4.0_GRCh38_af_greater_than_5.stripped.vcf.gz \
  ${GRCH38_DIR}/gnomad4.1_GRCh38_af_greater_than_5.stripped.vcf.gz
tabix -f ${OUTPUT}
