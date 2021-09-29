#!/bin/bash

# Having troubles with corrupted files downloading via FTP from NCBI via IPv6, http works ok

filename=ref_GRCh37.p5_top_level.gff3.gz
if [[ ! -e ${filename} ]]; then
  wget http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/BUILD.37.3/GFF/${filename}
fi

filename=ref_GRCh37.p9_top_level.gff3.gz
if [[ ! -e ${filename} ]]; then
  wget http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.103/GFF/${filename}
fi

filename=ref_GRCh37.p10_top_level.gff3.gz
if [[ ! -e ${filename} ]]; then
  wget http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.104/GFF/${filename}
fi

filename=ref_GRCh37.p13_top_level.gff3.gz
if [[ ! -e ${filename} ]]; then
  wget http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/GFF/${filename}
fi

# These all have the same name, so rename them based on release ID
for release in 105.20190906 105.20201022; do
  filename=GCF_000001405.25_GRCh37.p13_genomic.${release}.gff.gz
  if [[ ! -e ${filename} ]]; then
    wget http://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/${release}/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz --output-document=${filename}
  fi
done
