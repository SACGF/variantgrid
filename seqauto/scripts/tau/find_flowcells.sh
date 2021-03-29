#!/bin/sh

# Only look in idt_* dirs - as that's where new stuff is put
find /tau/data/clinical_hg38 -maxdepth 3 -name "SampleSheet.csv" -exec dirname {} \;

