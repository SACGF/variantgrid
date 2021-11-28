#!/bin/bash
# Based on: https://groups.google.com/a/soe.ucsc.edu/forum/#!topic/genome/n3As-4DRyLY

GENOME=$1
if [[ $GENOME != *.chrom.sizes ]]; then
	echo "Expected 1st argument to be xx.chrom.sizes" >&2
	exit 1
fi

OUTPUT_BIGWIG=$2
if [[ $OUTPUT_BIGWIG != *.bw ]]; then
        echo "2nd argument (output file) must end in .bw" >&2
        exit 1
fi

if [[ -e ${OUTPUT_BIGWIG} ]]; then
	echo "Output '${OUTPUT}' already exists!" >&2
	exit 1;
fi

TEMP_DIR="temp_$(basename ${OUTPUT_BIGWIG} .bw)"
mkdir ${TEMP_DIR}

for i in "${@:3}"
do
	echo "i=${i}"
	if [[ $i != *.wigFix.gz ]]; then
        	echo "'${i} must end in .wigFix.gz" >&2
        	exit 1
	fi

	NAME=$(basename $i .gz)
	echo $NAME
	gzcat $i | grep -v ^track > ${TEMP_DIR}/${NAME}.temp
	./wigToBigWig -clip -fixedSummaries -keepAllChromosomes ${TEMP_DIR}/${NAME}.temp ${GENOME} ${TEMP_DIR}/${NAME}.bw
done

./bigWigCat ${OUTPUT_BIGWIG} ${TEMP_DIR}/chr*.wigFix.bw
# rm -rf ${TEMP_DIR}
