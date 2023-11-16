#!/bin/bash

THIS_DIR=$(realpath "$(dirname "${BASH_SOURCE[0]}")")

# Structural variants
SV_COLUMNS=INFO/SVLEN,INFO/SVTYPE,INFO/END
COLS=INFO/AC,INFO/AN,INFO/AF
OTHER_COUNTS=INFO/N_HOMREF,INFO/N_HET,INFO/N_HOMALT
SUBPOPS=INFO/afr_AF,INFO/amr_AF,INFO/asj_AF,INFO/eas_AF,INFO/fin_AF,INFO/mid_AF,INFO/nfe_AF,INFO/oth_AF,INFO/sas_AF

KEEP_COLUMNS=${SV_COLUMNS},${COLS},${OTHER_COUNTS},${SUBPOPS}
CHROM_MAPPING_FILE=${THIS_DIR}/../../../snpdb/genome/chrom_mapping_GRCh38.map
GENOME_FASTA=/data/annotation/fasta/GCF_000001405.40_GRCh38.p14_genomic.fna.gz

# gnomad v4
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
  GNOMAD_VCF=gnomad.v4.0.sv.chr${chrom}.vcf.gz
  wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.0/genome_sv/${GNOMAD_VCF}
  wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.0/genome_sv/${GNOMAD_VCF}.tbi

  OUTPUT_VCF=
  # bcftools annotate --exclude 'AC=0' --remove '^{KEEP_COLUMNS}' --rename-chrs={CHROM_MAPPING_FILE} | vt normalize - -r ${GENOME_FASTA} -o + | vt uniq + -o ${OUTPUT_VCF}

done


# OTHER_INFOS = ["AC_popmax", "AN_popmax", "AF_popmax", "popmax", "nhomalt", "nhomalt_popmax", "nonpar"]
# GNOMAD_SUB_POPS = ["afr", "amr", "asj", "eas", "fin", "mid", "nfe", "oth", "sas"]  # Will get AF for each

# These have been removed in v4 - "AC_popmax", "AN_popmax", "AF_popmax"
# nonpar is now "par"


