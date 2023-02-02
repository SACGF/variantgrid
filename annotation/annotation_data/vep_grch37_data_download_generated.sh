#!/bin/bash

# Behind registration wall
# CosmicCodingMuts.normal.grch37.vcf.gz
# CosmicCodingMuts.normal.grch37.vcf.gz.tbi

BASE_URL=https://variantgrid.com/download/annotation/VEP/annotation_data/GRCh37

cd annotation_data/GRCh37

wget ${BASE_URL}/dbNSFP4.3a.grch37.stripped.gz
wget ${BASE_URL}/dbNSFP4.3a.grch37.stripped.gz.tbi

wget ${BASE_URL}/dbscSNV1.1_GRCh37.txt.gz
wget ${BASE_URL}/dbscSNV1.1_GRCh37.txt.gz.tbi

wget ${BASE_URL}/gnomad2.1.1_GRCh37_combined_af.vcf.bgz
wget ${BASE_URL}/gnomad2.1.1_GRCh37_combined_af.vcf.bgz.tbi

wget ${BASE_URL}/hg19.phastCons46way.placental.bw
wget ${BASE_URL}/hg19.phyloP46way.placental.bw

# Preprocessor files
# gnomad_GRCh37_af_greater_than_5.contigs.vcf.bgz


# Behind registration wall
# spliceai_scores.raw.indel.hg19.vcf.gz
# spliceai_scores.raw.indel.hg19.vcf.gz.tbi

# spliceai_scores.raw.snv.hg19.vcf.gz
# spliceai_scores.raw.snv.hg19.vcf.gz.tbi

cd ../..
