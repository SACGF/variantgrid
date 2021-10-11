#!/bin/bash

# v81 (points to 75) and earlier at GTFs that don't have transcript versions - just skip them

#82 is first GFF3 for GRCh37
#83 has no data
#84 is 82 again
#86 is 85 again
for release in 82 85 87; do
  filename=Homo_sapiens.GRCh37.${release}.gff3.gz
  if [[ ! -e ${filename} ]]; then
    wget ftp://ftp.ensembl.org/pub/grch37/release-${release}/gff3/homo_sapiens/${filename}
  fi
done