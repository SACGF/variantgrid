#!/bin/sh

# Only look in idt_* dirs - as that's where new stuff is put
# RTAComplete.txt is in 4_QC/sequencing_stats so need to drop down 2 dirs (with dirname)
find /tau/data/clinical_hg38 -maxdepth 5 -name "RTAComplete.txt" -exec dirname {} \; | xargs dirname | xargs dirname

