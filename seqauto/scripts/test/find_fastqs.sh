#!/bin/sh


OUTPUT_DIR=$1
FLOWCELL_DIRS=${OUTPUT_DIR}/find_flowcells.txt
FASTQ_DIR=Data/Intensities/Basecalls

while read -r line
do
	BASECALL_DIR=${line}/Data/Intensities/BaseCalls
	if [ -e ${BASECALL_DIR} ]; then
		find ${BASECALL_DIR} -name "*.fastq.gz"	
	fi
done < "${FLOWCELL_DIRS}"



