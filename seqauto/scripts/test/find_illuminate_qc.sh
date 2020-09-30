#!/bin/sh


OUTPUT_DIR=$1
FLOWCELL_DIRS=${OUTPUT_DIR}/find_flowcells.txt
FASTQ_DIR=Data/Intensities/Basecalls

while read -r line
do
	ILLUMINATE_DIR=${line}/Illuminate
	ls -d -1 ${ILLUMINATE_DIR}/* 2> /dev/null || true
done < "${FLOWCELL_DIRS}"



