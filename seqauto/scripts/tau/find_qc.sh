#!/bin/sh

OUTPUT_DIR=$1
FLOWCELL_DIRS=${OUTPUT_DIR}/find_flowcells.txt

while read -r line
do
	ls -d -1 ${line}/4_QC/exec_stats/*qc_summary.txt ${line}/4_QC/bam_stats/samples/*.per_gene_coverage.tsv.gz ${line}/0_goi/*.txt 2> /dev/null || true
done < "${FLOWCELL_DIRS}"
