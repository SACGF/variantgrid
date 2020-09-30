#!/bin/sh


OUTPUT_DIR=$1
FLOWCELL_DIRS=${OUTPUT_DIR}/find_flowcells.txt
BASECALL_SUB_DIR=Data/Intensities/BaseCalls

while read -r line
do
	# Files could be in subdirectories based on projects etc, hence we still use find.
	BASECALL_DIR=${line}/${BASECALL_SUB_DIR}
	if [ -e ${BASECALL_DIR} ]; then
		find ${BASECALL_DIR} -name "*.fastq.gz"	
	fi
done < "${FLOWCELL_DIRS}"



