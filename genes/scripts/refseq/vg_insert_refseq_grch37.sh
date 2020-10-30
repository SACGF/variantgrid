#!/bin/bash

VG_DIR=$(dirname $(dirname $(dirname $(dirname ${BASH_SOURCE[0]}))))

# Since we're replacing, insert from oldest -> newest
python3.8 ${VG_DIR}/manage.py import_gene_annotation --genome-build=GRCh37 --replace --annotation-consortium=RefSeq \
  --gff ref_GRCh37.p5_top_level.gff3.gz \
  --genePred ref_GRCh37.p5_top_level.genePred

python3.8 ${VG_DIR}/manage.py import_gene_annotation --genome-build=GRCh37 --replace --annotation-consortium=RefSeq \
  --gff ref_GRCh37.p9_top_level.gff3.gz \
  --genePred ref_GRCh37.p9_top_level.genePred

python3.8 ${VG_DIR}/manage.py import_gene_annotation --genome-build=GRCh37 --replace --annotation-consortium=RefSeq \
  --gff ref_GRCh37.p10_top_level.gff3.gz \
  --genePred ref_GRCh37.p10_top_level.genePred

python3.8 ${VG_DIR}/manage.py import_gene_annotation --genome-build=GRCh37 --replace --annotation-consortium=RefSeq \
  --gff ref_GRCh37.p13_top_level.gff3.gz \
  --genePred ref_GRCh37.p13_top_level.genePred

python3.8 ${VG_DIR}/manage.py import_gene_annotation --genome-build=GRCh37 --replace --annotation-consortium=RefSeq \
  --gff GCF_000001405.25_GRCh37.p13_genomic.105.20190906.gff.gz \
  --genePred GCF_000001405.25_GRCh37.p13_genomic.105.20190906.genePred

python3.8 ${VG_DIR}/manage.py import_gene_annotation --genome-build=GRCh37 --replace --annotation-consortium=RefSeq \
  --gff GCF_000001405.25_GRCh37.p13_genomic.105.20201022.gff.gz \
  --genePred GCF_000001405.25_GRCh37.p13_genomic.105.20201022.genePred
