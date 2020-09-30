#!/bin/bash

GENERATE_NEW="true"
GENERATE_NEW="false"

VARIANT_GRID_HOME="$(dirname $0)/.."
cd ${VARIANT_GRID_HOME}

if [[ ${GENERATE_NEW} = "true" ]]; then
	NUM_RECORDS=1000
	RANDOM_FILE_NAME=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 32 | head -n 1).random_${NUM_RECORDS}.vcf
	VCF_FILE=/tmp/${RANDOM_FILE_NAME}
	echo "Generating random ${VCF_FILE}"
	python3 manage.py random_vcf --num ${NUM_RECORDS} > ${VCF_FILE}
else
	VCF_FILE=${VARIANT_GRID_HOME}/seqauto/test_data/aligned/cancer_truSight/151120_AHISEQTEST/combined.vcf
fi

NAME=151120_HISEQTEST_combined

echo "Importing ${VCF_FILE}"
python3 manage.py import_vcf ${VCF_FILE} --name=${NAME} --user=dlawrence

