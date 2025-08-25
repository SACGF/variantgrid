#!/bin/sh

# Only look in idt_* dirs - as that's where new stuff is put
# RTAComplete.txt is in 4_QC/sequencing_stats so need to drop down 2 dirs (with dirname)
# find /tau/data/clinical_hg38 -maxdepth 5 -name "RTAComplete.txt" -exec dirname {} \; | xargs dirname | xargs dirname

# We want to be able to skip entire directories with .variantgrid_skip_flowcell in them - ie not descend any further
find /tau/data/clinical_hg38 -mindepth 1 -maxdepth 1 -type d ! -exec test -e '{}/.variantgrid_skip_flowcell' \; -print0 | xargs -0 -I{} find '{}' -maxdepth 4 -name 'RTAComplete.txt' | xargs -n1 dirname | xargs -n1 dirname | xargs -n1 dirname

