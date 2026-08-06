#!/bin/bash

# Public data downloaded straight from source. Everything else (plugin/custom data files) comes
# from "python3 manage.py download_annotation_data", which works the file list out from settings.

BASE_DIR=$(dirname ${BASH_SOURCE[0]})
ANNOTATION_DIR=/data/annotation

mkdir -p ${ANNOTATION_DIR}
cd ${ANNOTATION_DIR}

${BASE_DIR}/vep_all_builds_download.sh

${BASE_DIR}/vep_grch37_data_download.sh

${BASE_DIR}/vep_grch38_data_download.sh
