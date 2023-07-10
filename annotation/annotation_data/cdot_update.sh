#!/bin/bash

CDOT_VERSION=0.2.20
THIS_DIR=$(realpath "$(dirname "${BASH_SOURCE[0]}")")
VG_DIR=${THIS_DIR}/../..
DOWNLOAD_DIR=/tmp

python3 manage.py import_gene_annotation --annotation-consortium=RefSeq --genome-build=GRCh37 --json-file ${DOWNLOAD_DIR}/cdot-${CDOT_VERSION}.refseq.grch37.json.gz
python3 manage.py import_gene_annotation --annotation-consortium=RefSeq --genome-build=GRCh38 --json-file ${DOWNLOAD_DIR}/cdot-${CDOT_VERSION}.refseq.grch38.json.gz
python3 manage.py import_gene_annotation --annotation-consortium=Ensembl --genome-build=GRCh37 --json-file ${DOWNLOAD_DIR}/cdot-${CDOT_VERSION}.ensembl.grch37.json.gz
python3 manage.py import_gene_annotation --annotation-consortium=Ensembl --genome-build=GRCh38 --json-file ${DOWNLOAD_DIR}/cdot-${CDOT_VERSION}.ensembl.grch38.json.gz

