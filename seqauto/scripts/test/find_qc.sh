#!/bin/bash

# This needs to return all QC objects

TEST_DATA_DIR=$(dirname ${BASH_SOURCE[0]})/../../test_data
TEST_DATA_DIR=$(cd "${TEST_DATA_DIR}"; pwd) # Absolute path

find ${TEST_DATA_DIR} -name "*_stats.txt" -o -name "*_qc_summary.txt" -o -name "*.per_gene_coverage.tsv.gz"
