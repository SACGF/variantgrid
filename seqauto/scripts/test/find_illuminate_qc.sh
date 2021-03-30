#!/bin/sh


OUTPUT_DIR=$1
FLOWCELL_DIRS=${OUTPUT_DIR}/find_flowcells.txt

while read -r line
do
	ls -d -1 ${line}/4_QC/sequencing_stats/Illuminate/illuminate_report.txt 2> /dev/null || true
done < "${FLOWCELL_DIRS}"
