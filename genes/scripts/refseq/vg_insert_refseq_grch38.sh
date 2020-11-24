#!/bin/bash

VG_DIR=$(dirname $(dirname $(dirname $(dirname ${BASH_SOURCE[0]}))))

# Since we're replacing, insert from oldest -> newest
python3.8 ${VG_DIR}/manage.py import_gene_annotation --genome-build=GRCh38 --replace --annotation-consortium=RefSeq \
  --gff ref_GRCh38_top_level.gff3.gz \
  --genePred ref_GRCh38_top_level.genePred

python3.8 ${VG_DIR}/manage.py import_gene_annotation --genome-build=GRCh38 --replace --annotation-consortium=RefSeq \
  --gff ref_GRCh38.p2_top_level.gff3.gz \
  --genePred ref_GRCh38.p2_top_level.genePred

python3.8 ${VG_DIR}/manage.py import_gene_annotation --genome-build=GRCh38 --replace --annotation-consortium=RefSeq \
  --gff ref_GRCh38.p7_top_level.gff3.gz \
  --genePred ref_GRCh38.p7_top_level.genePred

python3.8 ${VG_DIR}/manage.py import_gene_annotation --genome-build=GRCh38 --replace --annotation-consortium=RefSeq \
  --gff ref_GRCh38.p12_top_level.gff3.gz \
  --genePred ref_GRCh38.p12_top_level.genePred

python3.8 ${VG_DIR}/manage.py import_gene_annotation --genome-build=GRCh38 --replace --annotation-consortium=RefSeq \
  --gff GCF_000001405.38_GRCh38.p12_genomic.gff.gz \
  --genePred GCF_000001405.38_GRCh38.p12_genomic.genePred

python3.8 ${VG_DIR}/manage.py import_gene_annotation --genome-build=GRCh38 --replace --annotation-consortium=RefSeq \
  --gff GCF_000001405.39_GRCh38.p13_genomic.109.20190607.gff.gz \
  --genePred GCF_000001405.39_GRCh38.p13_genomic.109.20190607.genePred

python3.8 ${VG_DIR}/manage.py import_gene_annotation --genome-build=GRCh38 --replace --annotation-consortium=RefSeq \
  --gff GCF_000001405.39_GRCh38.p13_genomic.109.20190905.gff.gz \
  --genePred GCF_000001405.39_GRCh38.p13_genomic.109.20190905.genePred

python3.8 ${VG_DIR}/manage.py import_gene_annotation --genome-build=GRCh38 --replace --annotation-consortium=RefSeq \
  --gff GCF_000001405.39_GRCh38.p13_genomic.109.20191205.gff.gz \
  --genePred GCF_000001405.39_GRCh38.p13_genomic.109.20191205.genePred

python3.8 ${VG_DIR}/manage.py import_gene_annotation --genome-build=GRCh38 --replace --annotation-consortium=RefSeq \
  --gff GCF_000001405.39_GRCh38.p13_genomic.109.20200228.gff.gz \
  --genePred GCF_000001405.39_GRCh38.p13_genomic.109.20200228.genePred

python3.8 ${VG_DIR}/manage.py import_gene_annotation --genome-build=GRCh38 --replace --annotation-consortium=RefSeq \
  --gff GCF_000001405.39_GRCh38.p13_genomic.109.20200522.gff.gz \
  --genePred GCF_000001405.39_GRCh38.p13_genomic.109.20200522.genePred

python3.8 ${VG_DIR}/manage.py import_gene_annotation --genome-build=GRCh38 --replace --annotation-consortium=RefSeq \
  --gff GCF_000001405.39_GRCh38.p13_genomic.109.20200815.gff.gz \
  --genePred GCF_000001405.39_GRCh38.p13_genomic.109.20200815.genePred
