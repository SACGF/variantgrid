#!/bin/bash

filename=Homo_sapiens.GRCh38.78.gtf.gz
if [[ ! -e ${filename} ]]; then
  wget ftp://ftp.ensembl.org/pub/release-78/gtf/homo_sapiens/${filename}
fi

for release in 76 77 78 79 80; do
  filename=Homo_sapiens.GRCh38.${release}.gtf.gz
  if [[ ! -e ${filename} ]]; then
    wget ftp://ftp.ensembl.org/pub/release-${release}/gtf/homo_sapiens/${filename}
  fi
done

#81 is first GFF3 for GRCh38
for release in 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104; do
  filename=Homo_sapiens.GRCh38.${release}.gff3.gz
  if [[ ! -e ${filename} ]]; then
    wget ftp://ftp.ensembl.org/pub/release-${release}/gff3/homo_sapiens/${filename}
  fi
done