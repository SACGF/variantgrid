#!/bin/bash

# 75 is last GRCh37 in main directory
for release in 60 65 70 75; do
  wget ftp://ftp.ensembl.org/pub/release-${release}/gtf/homo_sapiens/Homo_sapiens.GRCh37.${release}.gtf.gz
done

wget ftp://ftp.ensembl.org/pub/grch37/release-81/gtf/homo_sapiens/Homo_sapiens.GRCh37.81.gtf.gz

#82 is first GFF3 for GRCh37
for release in 82 83 84 85 86 87; do
  wget ftp://ftp.ensembl.org/pub/grch37/release-${release}/gff3/homo_sapiens/Homo_sapiens.GRCh37.${release}.gff3.gz
done