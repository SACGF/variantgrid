#!/bin/sh

# This needs to return all QC objects
TAU_DATA_DIR=/tau/data/clinical_hg38

find ${TAU_DATA_DIR} -name "*_stats.txt" -o -name "*_qc_summary.txt" -o -name "*.per_gene_coverage.tsv.gz"
