#!/bin/sh

# We want to be able to skip entire directories with .variantgrid_skip_flowcell in them - ie not descend any further
find /tau/data/clinical_hg38 -mindepth 1 -maxdepth 1 -type d ! -exec test -e '{}/.variantgrid_skip_flowcell' \; -print0 \
| xargs -0 -I{} find '{}' \( -name \*.vcf -o -name \*.vcf.gz \)
