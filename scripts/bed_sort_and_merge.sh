#!/bin/bash

BED_FILE=$1

if [[ -z ${BED_FILE} ]]; then
	echo "Usage $(basename $0): bed_file" >&2
	exit 1;
fi

if [ ! -e ${BED_FILE} ]; then
	echo "No such bed file: '${BED_FILE}'" >&2
	exit 1;
fi


VARIANT_GRID_HOME="$(dirname $0)/.."
cd ${VARIANT_GRID_HOME}

bedtools sort -i ${BED_FILE} | bedtools merge -i stdin


