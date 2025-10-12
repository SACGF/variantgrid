#!/bin/sh

# BAMs are located eg:
# /tau/data/clinical_hg38/idt_exome/Exome_21_110_211126_NB501008_0158_AH3LJHBGXL/1_BAM/26_GU_TRIO_PID_AFF_AE_2131009254B_recal_reads.hg38.bam
# So only need to go in 4
# find /tau/data/clinical_hg38 -maxdepth 4 -name "*.bam"

OUTPUT_DIR=$1
FLOWCELL_DIRS=${OUTPUT_DIR}/find_flowcells.txt

while read -r line
do
	ls -d -1 ${line}/1_BAM/*.bam 2> /dev/null || true
done < "${FLOWCELL_DIRS}"
