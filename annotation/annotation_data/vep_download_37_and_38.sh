#!/bin/bash

BASE_DIR=$(dirname ${BASH_SOURCE[0]})
ANNOTATION_DIR=/data/annotation

mkdir -p ${ANNOTATION_DIR}
cd ${ANNOTATION_DIR}

${BASE_DIR}/vep_all_builds_download.sh

${BASE_DIR}/vep_grch37_data_download.sh
${BASE_DIR}/vep_grch37_data_download_generated.sh

${BASE_DIR}/vep_grch38_data_download.sh
${BASE_DIR}/vep_grch38_data_download_generated.sh
