#!/bin/bash

export PATH=${PATH}:/hpcfs/groups/phoenix-hpc-sacgf/tools/tabix-0.2.6:/hpcfs/groups/phoenix-hpc-sacgf/tools/bcftools/current/bcftools

# THIS_DIR=$(realpath "$(dirname "${BASH_SOURCE[0]}")")
THIS_DIR=/hpcfs/groups/phoenix-hpc-sacgf/reference/hg38/Misce/gnomAD4/sv
cd ${THIS_DIR}

# Structural variants
SV_COLUMNS=INFO/SVLEN,INFO/SVTYPE,INFO/END
COLS=INFO/AC,INFO/AN,INFO/AF
OTHER_COUNTS=INFO/N_HOMREF,INFO/N_HET,INFO/N_HOMALT,INFO/POPMAX_AF,INFO/PAR
SUBPOPS=INFO/afr_AF,INFO/amr_AF,INFO/asj_AF,INFO/eas_AF,INFO/fin_AF,INFO/mid_AF,INFO/nfe_AF,INFO/oth_AF,INFO/sas_AF

KEEP_COLUMNS=${SV_COLUMNS},${COLS},${OTHER_COUNTS},${SUBPOPS}
MAPPING_DIR=$(dirname ${THIS_DIR})
CHROM_MAPPING_FILE=${MAPPING_DIR}/chrom_mapping_GRCh38.map
MERGE_VCF=gnomad.v4.0.sv.merged.vcf

# gnomad v4
merge_args=()
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
  GNOMAD_VCF=gnomad.v4.0.sv.chr${chrom}.vcf.gz
  if [ ! -e ${GNOMAD_VCF} ]; then
    wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.0/genome_sv/${GNOMAD_VCF} &
    wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.0/genome_sv/${GNOMAD_VCF}.tbi &
  fi
done

echo "Waiting for downloads to finish..."
wait
echo "Downloads done"

for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
  GNOMAD_VCF=gnomad.v4.0.sv.chr${chrom}.vcf.gz
  OUTPUT_VCF=gnomad.v4.0.sv.chr${chrom}.converted.vcf.gz
  echo "Going from ${GNOMAD_VCF} -> ${OUTPUT_VCF}"

  # Dont' normalize as is mostly "N" refs
  bcftools annotate --exclude 'AC=0' --remove "^${KEEP_COLUMNS}" ${GNOMAD_VCF} -o ${OUTPUT_VCF} &
  merge_args+=(${OUTPUT_VCF})
done

echo "Waiting for annotate jobs to finish..."
wait

bcftools concat --output-type v --output ${MERGE_VCF} ${merge_args[@]};
bgzip ${MERGE_VCF}
tabix -p vcf ${MERGE_VCF}.gz
