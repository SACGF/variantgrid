#!/bin/bash

# this is for pre-processing VCFs to split common/uncommon variants
# @see https://github.com/SACGF/variantgrid/issues/521#issuecomment-1042527037

set -e

if [[ -z ${VEP_ANNOTATION_DIR} ]]; then
  echo "Please set 'VEP_ANNOTATION_DIR' variable to point where gnomAD files are" >&1
  exit 1
fi

GENOME_BUILD=T2T-CHM13v2.0
# 4.1.t2t_liftover_ref_swaps is gnomAD's liftover plus the sites it drops (gnomad_t2t_ref_swaps.sh, #2788)
GNOMAD_VERSION_T2T=${GNOMAD_VERSION_T2T:-4.1.t2t_liftover_ref_swaps}

ANNOTATION_DATA_DIR=${VEP_ANNOTATION_DIR}/annotation_data/T2T-CHM13v2.0
GNOMAD_VCF=${ANNOTATION_DATA_DIR}/gnomad${GNOMAD_VERSION_T2T}_T2T-CHM13v2.0_combined_af.vcf.bgz
OUTPUT_NAME=gnomad${GNOMAD_VERSION_T2T}_${GENOME_BUILD}_af_greater_than_5
OUTPUT_AF=${ANNOTATION_DATA_DIR}/${OUTPUT_NAME}.vcf.gz
OUTPUT_STRIPPED=${ANNOTATION_DATA_DIR}/${OUTPUT_NAME}.stripped.vcf.gz
if [[ ! -e ${OUTPUT_STRIPPED} ]]; then
  echo "Generating ${OUTPUT_STRIPPED}..."
  if [[ ! -e ${OUTPUT_AF} ]]; then
    echo "Extracting AF>0.05..."
    bcftools view -i "AF>0.05" ${GNOMAD_VCF} --output-type z --output-file ${OUTPUT_AF}
    tabix ${OUTPUT_AF}
  fi
  # As text: the combined file has an INFO field its header doesn't declare (non_par), which
  # bcftools annotate refuses to write. Its contigs are already accessions, so there's nothing to rename
  echo "Stripping info"
  zcat ${OUTPUT_AF} \
    | awk -F'\t' 'BEGIN {OFS = "\t"} /^#/ {print; next}
                  {af = "."; n = split($8, info, ";"); for (i = 1; i <= n; i++) if (info[i] ~ /^AF=/) af = substr(info[i], 4)
                   $8 = "AF=" af; print}' \
    | bgzip > ${OUTPUT_STRIPPED}.tmp
  mv ${OUTPUT_STRIPPED}.tmp ${OUTPUT_STRIPPED}
  tabix ${OUTPUT_STRIPPED}
fi
