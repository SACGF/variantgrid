#!/bin/sh

OUTPUT_DIR=$1
FASTQS=${OUTPUT_DIR}/find_fastqs.txt

while read -r line
do
	FASTQ_DIR=$(dirname ${line})
	FASTQC=${FASTQ_DIR}/FastQC/$(basename ${line} .fastq.gz)_fastqc/fastqc_data.txt
	if [ -e ${FASTQC} ]; then
		echo ${FASTQC}
	fi
done < "${FASTQS}"



