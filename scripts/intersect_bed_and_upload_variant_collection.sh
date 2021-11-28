#!/bin/bash

BED_FILE=$1
VARIANT_COLLECTION_ID=$2

echo "${0} ${BED_FILE} ${VARIANT_COLLECTION_ID}" >&2

if [[ -z ${BED_FILE} || -z ${VARIANT_COLLECTION_ID} ]]; then
	echo "Usage $(basename $0): bed_file variant_collection_id" >&2
	exit 1;
fi

if [ ! -e ${BED_FILE} ]; then
	echo "No such bed file: '${BED_FILE}'" >&2
	exit 1;
fi


VARIANT_GRID_HOME="$(dirname $0)/.."
cd ${VARIANT_GRID_HOME}

bedtools intersect -wa -u -header -sorted -a stdin -b ${BED_FILE} | python3 manage.py stdin_to_variant_collection ${VARIANT_COLLECTION_ID}
