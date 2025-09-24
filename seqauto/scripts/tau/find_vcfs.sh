#!/bin/sh

#!/bin/sh


OUTPUT_DIR=$1
FLOWCELL_DIRS=${OUTPUT_DIR}/find_flowcells.txt

while read -r line
do
	ls -d -1 ${line}/2_variants/*.vcf.gz ${line}/2_variants/**/*.vcf ${line}/2_variants/**/*.vcf.gz 2> /dev/null || true
done < "${FLOWCELL_DIRS}"
