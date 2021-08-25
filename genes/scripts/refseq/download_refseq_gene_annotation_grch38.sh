#!/bin/bash

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.106/GFF/ref_GRCh38_top_level.gff3.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.107/GFF/ref_GRCh38.p2_top_level.gff3.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.108/GFF/ref_GRCh38.p7_top_level.gff3.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.109/GFF/ref_GRCh38.p12_top_level.gff3.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.gff.gz

# These all have the same name, so rename them based on release ID

for release in 109.20190607 109.20190905 109.20191205 109.20200228 109.20200522 109.20200815 109.20201120 109.20210226 109.20210514; do
  wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/${release}/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz --output-document=GCF_000001405.39_GRCh38.p13_genomic.${release}.gff.gz
done
