#!/bin/bash

CDOT_VERSION=0.2.14
THIS_DIR=$(realpath "$(dirname "${BASH_SOURCE[0]}")")
VG_DIR=${THIS_DIR}/../..
DOWNLOAD_DIR=/tmp

echo "Downloading data in ${DOWNLOAD_DIR}"
cd ${DOWNLOAD_DIR}

wget \
  https://github.com/SACGF/cdot/releases/download/v${CDOT_VERSION}/cdot-${CDOT_VERSION}.ensembl.grch37.json.gz \
  https://github.com/SACGF/cdot/releases/download/v${CDOT_VERSION}/cdot-${CDOT_VERSION}.ensembl.grch38.json.gz \
  https://github.com/SACGF/cdot/releases/download/v${CDOT_VERSION}/cdot-${CDOT_VERSION}.refseq.grch37.json.gz \
  https://github.com/SACGF/cdot/releases/download/v${CDOT_VERSION}/cdot-${CDOT_VERSION}.refseq.grch38.json.gz

cd ${VG_DIR}

python3 manage.py import_gene_annotation --annotation-consortium=RefSeq --genome-build=GRCh37 --json-file ${DOWNLOAD_DIR}/cdot-${CDOT_VERSION}.refseq.grch37.json.gz
python3 manage.py import_gene_annotation --annotation-consortium=RefSeq --genome-build=GRCh38 --json-file ${DOWNLOAD_DIR}/cdot-${CDOT_VERSION}.refseq.grch38.json.gz
python3 manage.py import_gene_annotation --annotation-consortium=Ensembl --genome-build=GRCh37 --json-file ${DOWNLOAD_DIR}/cdot-${CDOT_VERSION}.ensembl.grch37.json.gz
python3 manage.py import_gene_annotation --annotation-consortium=Ensembl --genome-build=GRCh38 --json-file ${DOWNLOAD_DIR}/cdot-${CDOT_VERSION}.ensembl.grch38.json.gz

