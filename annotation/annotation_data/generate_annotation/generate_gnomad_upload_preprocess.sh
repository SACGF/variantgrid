#!/bin/bash

# this is for pre-processing VCFs to split common/uncommon variants
# @see https://github.com/SACGF/variantgrid/issues/521#issuecomment-1042527037

set -e

if [[ -z ${VEP_ANNOTATION_DIR} ]]; then
  echo "Please set 'VEP_ANNOTATION_DIR' variable to point where gnomAD files are" >&1
  exit 1
fi

GRCH37_DIR=${VEP_ANNOTATION_DIR}/annotation_data/GRCh37
GRCH38_DIR=${VEP_ANNOTATION_DIR}/annotation_data/GRCh38

bcftools view -i "AF>0.05" ${GRCH37_DIR}/gnomad2.1.1_GRCh37_combined_af.vcf.bgz --output-type z --output-file gnomad_GRCh37_af_greater_than_5.vcf.bgz
tabix gnomad_GRCh37_af_greater_than_5.vcf.bgz
bcftools annotate -x ^INFO/AF --rename-chrs chr_contig_GRCh37.map --output-type z --output gnomad_GRCh37_af_greater_than_5.contigs.vcf.bgz gnomad_GRCh37_af_greater_than_5.vcf.bgz
tabix gnomad_GRCh37_af_greater_than_5.contigs.vcf.bgz

bcftools view -i "AF>0.05" ${GRCH38_DIR}/gnomad3.1_GRCh38_merged.vcf.bgz --output-type z --output-file gnomad_GRCh38_af_greater_than_5.vcf.bgz
tabix gnomad_GRCh38_af_greater_than_5.vcf.bgz
bcftools annotate -x ^INFO/AF --rename-chrs chr_contig_GRCh38.map --output-type z --output gnomad_GRCh38_af_greater_than_5.contigs.vcf.bgz gnomad_GRCh38_af_greater_than_5.vcf.bgz
tabix gnomad_GRCh38_af_greater_than_5.contigs.vcf.bgz
