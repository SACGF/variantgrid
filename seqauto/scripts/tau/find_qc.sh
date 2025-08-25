#!/bin/sh

# This needs to return all QC objects

# We want to be able to skip entire directories with .variantgrid_skip_flowcell in them - ie not descend any further
find /tau/data/clinical_hg38 -mindepth 1 -maxdepth 1 -type d ! -exec test -e '{}/.variantgrid_skip_flowcell' \; -print0 \
| xargs -0 -I{} find '{}' -name "*_stats.txt" -o -name "*_qc_summary.txt" -o -name "*.per_gene_coverage.tsv.gz" -o  -path "*/0_goi/*.txt"
