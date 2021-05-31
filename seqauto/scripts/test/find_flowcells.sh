#!/bin/bash

TEST_DATA_DIR=$(dirname ${BASH_SOURCE[0]})/../../test_data
TEST_DATA_DIR=$(cd "${TEST_DATA_DIR}"; pwd) # Absolute path

# RTAComplete.txt is in 4_QC/sequencing_stats so need to drop down 2 dirs (with dirname)
find ${TEST_DATA_DIR}/clinical_hg38 -maxdepth 6 -name "RTAComplete.txt" -exec dirname {} \; | xargs dirname | xargs dirname

