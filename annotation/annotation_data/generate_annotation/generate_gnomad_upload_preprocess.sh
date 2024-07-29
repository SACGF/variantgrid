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
GNOMAD_VERSION_GRCH38=4.0



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


GRCH38_DIR=${VEP_ANNOTATION_DIR}/annotation_data/GRCh38
CHROM_MAPPING_38=${CHROM_MAPPING_DIR}/chrom_mapping_GRCh38.map
OUTPUT_NAME_GRCH38=gnomad${GNOMAD_VERSION_GRCH38}_GRCh38_af_greater_than_5
OUTPUT_AF_GRCH38=${GRCH38_DIR}/${OUTPUT_NAME_GRCH38}.vcf.gz
OUTPUT_STRIPPED_GRCH38=${GRCH38_DIR}/${OUTPUT_NAME_GRCH38}.stripped.vcf.gz
if [[ ! -e ${OUTPUT_STRIPPED_GRCH38} ]]; then
  echo "Generating ${OUTPUT_STRIPPED_GRCH38}..."
  if [[ ! -e ${OUTPUT_AF_GRCH38} ]]; then
    echo "Extracting AF>0.05..."
    bcftools view -i "AF>0.05" ${GRCH38_DIR}/gnomad${GNOMAD_VERSION_GRCH38}_GRCh38_combined_af.vcf.bgz --output-type z --output-file ${OUTPUT_AF_GRCH38}
    tabix ${OUTPUT_AF_GRCH38}
  fi
  echo "Stripping info"
  bcftools annotate --rename-chrs ${CHROM_MAPPING_38} -x ^INFO/AF --output-type z --output ${OUTPUT_STRIPPED_GRCH38} ${OUTPUT_AF_GRCH38}
  tabix ${OUTPUT_STRIPPED_GRCH38}
fi