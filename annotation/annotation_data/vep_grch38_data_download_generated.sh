#!/bin/bash

# Behind registration wall
# CosmicCodingMuts.normal.grch37.vcf.gz
# CosmicCodingMuts.normal.grch37.vcf.gz.tbi

BASE_URL=https://variantgrid.com/static/annotation_download

cd annotation_data/GRCh38

wget ${BASE_URL}/dbNSFP4.3a.grch38.stripped.gz
wget ${BASE_URL}/dbNSFP4.3a.grch38.stripped.gz.tbi

wget ${BASE_URL}/dbscSNV1.1_GRCh38.txt.gz
wget ${BASE_URL}/dbscSNV1.1_GRCh38.txt.gz.tbi


wget ${BASE_URL}/gnomad2.1.1_GRCh38_combined_af.vcf.bgz
wget ${BASE_URL}/gnomad2.1.1_GRCh38_combined_af.vcf.bgz.tbi

wget ${BASE_URL}/gnomad3.1_GRCh38_merged.vcf.bgz
wget ${BASE_URL}/gnomad3.1_GRCh38_merged.vcf.bgz.tbi

# Preprocessor files
# gnomad_GRCh37_af_greater_than_5.contigs.vcf.bgz

# Behind registration wall

# annotation_data/GRCh38/CosmicCodingMuts_v95_20211101_grch38.normal.vcf.gz",
#             "mastermind": "annotation_data/GRCh38/mastermind_cited_variants_reference-2022.04.02-grch38.vcf.gz",
#             "maxentscan": "annotation_data/all_builds/maxentscan",
#            "repeatmasker": "annotation_data/GRCh38/repeatmasker_hg38.bed.gz",
#            "spliceai_snv": "annotation_data/GRCh38/spliceai_scores.raw.snv.hg38.vcf.gz",
#            "spliceai_indel": "annotation_data/GRCh38/spliceai_scores.raw.indel.hg38.vcf.gz",

cd ../..
