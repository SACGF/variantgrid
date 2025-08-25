#!/bin/sh

# BAMs are located eg:
# /tau/data/clinical_hg38/idt_exome/Exome_21_110_211126_NB501008_0158_AH3LJHBGXL/1_BAM/26_GU_TRIO_PID_AFF_AE_2131009254B_recal_reads.hg38.bam
# So only need to go in 4
find /tau/data/clinical_hg38 -maxdepth 4 -name "*.bam"

