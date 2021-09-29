#!/bin/bash

# 75 is last GRCh37 in main directory
for release in 60 65 70 75; do
  filename=Homo_sapiens.GRCh37.${release}.gtf.gz
  if [[ ! -e ${filename} ]]; then
    wget ftp://ftp.ensembl.org/pub/release-${release}/gtf/homo_sapiens/${filename}
  fi
done

filename=Homo_sapiens.GRCh37.81.gtf.gz
if [[ ! -e ${filename} ]]; then
  wget ftp://ftp.ensembl.org/pub/grch37/release-81/gtf/homo_sapiens/${filename}
fi

#82 is first GFF3 for GRCh37
for release in 82 83 84 85 86 87; do
  filename=Homo_sapiens.GRCh37.${release}.gff3.gz
  if [[ ! -e ${filename} ]]; then
    wget ftp://ftp.ensembl.org/pub/grch37/release-${release}/gff3/homo_sapiens/${filename}
  fi
done