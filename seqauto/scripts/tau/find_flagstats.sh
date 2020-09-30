#!/bin/sh

OUTPUT_DIR=$1
BAMS=${OUTPUT_DIR}/find_bams.txt

while read -r line
do
	BAM_DIR=$(dirname ${line})
	FLAGSTATS=${BAM_DIR}/flagstats/$(basename ${line} .bam).flagstat.txt
	if [ -e ${FLAGSTATS} ]; then
		echo ${FLAGSTATS}
	fi
done < "${BAMS}"



