#!/bin/bash

# This needs illuminate installed via Pip (run with env runner)
set -e

FLOWCELL_DIR=$1
ILLUMINATE_DIR=${FLOWCELL_DIR}/Illuminate

mkdir -p ${ILLUMINATE_DIR}
illuminate --quality --index --tile ${FLOWCELL_DIR} > ${ILLUMINATE_DIR}/illuminate_report.txt
