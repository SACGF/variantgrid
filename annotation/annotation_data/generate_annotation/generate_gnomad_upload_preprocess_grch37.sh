#!/bin/bash

# this is for pre-processing VCFs to split common/uncommon variants
# @see https://github.com/SACGF/variantgrid/issues/521#issuecomment-1042527037

set -e

if [[ -z ${VEP_ANNOTATION_DIR} ]]; then
  echo "Please set 'VEP_ANNOTATION_DIR' variable to point where gnomAD files are" >&1
  exit 1
fi

VARIANTGRID_DIR=$(dirname $(dirname $(dirname $(dirname ${BASH_SOURCE[0]}))))
CHROM_MAPPING_DIR=${VARIANTGRID_DIR}/snpdb/genome

GNOMAD_VERSION_GRCH37=2.1.1


GRCH37_DIR=${VEP_ANNOTATION_DIR}/annotation_data/GRCh37
CHROM_MAPPING_37=${CHROM_MAPPING_DIR}/chrom_mapping_GRCh37.map
OUTPUT_NAME_GRCH37=gnomad${GNOMAD_VERSION_GRCH37}_GRCh37_af_greater_than_5
OUTPUT_AF_GRCH37=${GRCH37_DIR}/${OUTPUT_NAME_GRCH37}.vcf.gz
OUTPUT_STRIPPED_GRCH37=${GRCH37_DIR}/${OUTPUT_NAME_GRCH37}.stripped.vcf.gz
if [[ ! -e ${OUTPUT_STRIPPED_GRCH37} ]]; then
  echo "Generating ${OUTPUT_STRIPPED_GRCH37}..."
  if [[ ! -e ${OUTPUT_AF_GRCH37} ]]; then
    echo "Extracting AF>0.05..."
    bcftools view -i "AF>0.05" ${GRCH37_DIR}/gnomad${GNOMAD_VERSION_GRCH37}_GRCh37_combined_af.vcf.bgz --output-type z --output-file ${OUTPUT_AF_GRCH37}
    tabix ${OUTPUT_AF_GRCH37}
  fi
  echo "Stripping INFO"
  bcftools annotate --rename-chrs ${CHROM_MAPPING_37} -x ^INFO/AF --output-type z --output ${OUTPUT_STRIPPED_GRCH37} ${OUTPUT_AF_GRCH37}
  tabix ${OUTPUT_STRIPPED_GRCH37}
fi
