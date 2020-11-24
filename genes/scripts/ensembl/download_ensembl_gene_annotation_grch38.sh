#!/bin/bash

ftp://ftp.ensembl.org/pub/release-78/gtf/homo_sapiens/Homo_sapiens.GRCh38.78.gtf.gz

for release in 76 77 78 79 80; do
  wget ftp://ftp.ensembl.org/pub/release-${release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${release}.gtf.gz
done

#81 is first GFF3 for GRCh38
for release in 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101; do
  wget ftp://ftp.ensembl.org/pub/release-${release}/gff3/homo_sapiens/Homo_sapiens.GRCh38.${release}.gff3.gz
done