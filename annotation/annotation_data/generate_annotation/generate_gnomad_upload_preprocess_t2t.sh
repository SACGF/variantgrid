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

GENOME_BUILD=T2T-CHM13v2.0
GNOMAD_VERSION_T2T=4.1.t2t_liftover

ANNOTATION_DATA_DIR=${VEP_ANNOTATION_DIR}/annotation_data/T2T-CHM13v2.0
CHROM_MAPPING=${CHROM_MAPPING_DIR}/chrom_mapping_${GENOME_BUILD}.map
GNOMAD_VCF=${ANNOTATION_DATA_DIR}/gnomad4.1.t2t_liftover_T2T-CHM13v2.0_combined_af.vcf.bgz
OUTPUT_NAME=gnomad${GNOMAD_VERSION_T2T}_${GENOME_BUILD}_af_greater_than_5
OUTPUT_AF=${ANNOTATION_DATA_DIR}/${OUTPUT_NAME}.vcf.gz
OUTPUT_STRIPPED=${GRCH38_DIR}/${OUTPUT_NAME}.stripped.vcf.gz
if [[ ! -e ${OUTPUT_STRIPPED} ]]; then
  echo "Generating ${OUTPUT_STRIPPED}..."
  if [[ ! -e ${OUTPUT_AF} ]]; then
    echo "Extracting AF>0.05..."
    bcftools view -i "AF>0.05" ${GNOMAD_VCF} --output-type z --output-file ${OUTPUT_AF}
    tabix ${OUTPUT_AF}
  fi
  echo "Stripping info"
  bcftools annotate --rename-chrs ${CHROM_MAPPING} -x ^INFO/AF --output-type z --output ${OUTPUT_STRIPPED} ${OUTPUT_AF}
  tabix ${OUTPUT_STRIPPED}
fi
