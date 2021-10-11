#!/bin/bash

# Having troubles with corrupted files downloading via FTP from NCBI via IPv6, http works ok

filename=ref_GRCh38_top_level.gff3.gz
if [[ ! -e ${filename} ]]; then
  wget http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.106/GFF/${filename}
fi

filename=ref_GRCh38.p2_top_level.gff3.gz
if [[ ! -e ${filename} ]]; then
  wget http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.107/GFF/${filename}
fi

filename=ref_GRCh38.p7_top_level.gff3.gz
if [[ ! -e ${filename} ]]; then
  wget http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.108/GFF/${filename}
fi

filename=ref_GRCh38.p12_top_level.gff3.gz
if [[ ! -e ${filename} ]]; then
  wget http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.109/GFF/${filename}
fi

filename=GCF_000001405.38_GRCh38.p12_genomic.gff.gz
if [[ ! -e ${filename} ]]; then
  wget http://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109/GCF_000001405.38_GRCh38.p12/${filename}
fi
# These all have the same name, so rename them based on release ID

for release in 109.20190607 109.20190905 109.20191205 109.20200228 109.20200522 109.20200815 109.20201120 109.20210226 109.20210514; do
  filename=GCF_000001405.39_GRCh38.p13_genomic.${release}.gff.gz
  if [[ ! -e ${filename} ]]; then
    wget http://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/${release}/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz --output-document=${filename}
  fi
done
