#!/bin/bash

# Behind registration wall
# CosmicCodingMuts.normal.grch37.vcf.gz
# CosmicCodingMuts.normal.grch37.vcf.gz.tbi

BASE_URL=https://variantgrid.com/download/annotation/VEP/annotation_data/GRCh38

cd annotation_data/GRCh38

wget ${BASE_URL}/CosmicCodingMuts_v95_20211101_grch38.normal.vcf.gz
wget ${BASE_URL}/CosmicCodingMuts_v95_20211101_grch38.normal.vcf.gz.tbi
wget ${BASE_URL}/TOPMED_GRCh38_20180418.vcf.gz
wget ${BASE_URL}/TOPMED_GRCh38_20180418.vcf.gz.tbi
wget ${BASE_URL}/UK10K_COHORT.20160215.sites.GRCh38.vcf.gz
wget ${BASE_URL}/UK10K_COHORT.20160215.sites.GRCh38.vcf.gz.tbi
wget ${BASE_URL}/dbNSFP4.3a.grch38.stripped.gz
wget ${BASE_URL}/dbNSFP4.3a.grch38.stripped.gz.tbi
wget ${BASE_URL}/dbscSNV1.1_GRCh38.txt.gz
wget ${BASE_URL}/dbscSNV1.1_GRCh38.txt.gz.tbi
wget ${BASE_URL}/gnomad2.1.1_GRCh38_combined_af.vcf.bgz
wget ${BASE_URL}/gnomad2.1.1_GRCh38_combined_af.vcf.bgz.tbi
wget ${BASE_URL}/gnomad3.1_GRCh38_merged.vcf.bgz
wget ${BASE_URL}/gnomad3.1_GRCh38_merged.vcf.bgz.tbi
wget ${BASE_URL}/gnomad_GRCh38_af_greater_than_5.contigs.vcf.bgz
wget ${BASE_URL}/gnomad_GRCh38_af_greater_than_5.contigs.vcf.bgz.tbi
wget ${BASE_URL}/hg38.phastCons100way.bw
wget ${BASE_URL}/hg38.phastCons30way.bw
wget ${BASE_URL}/hg38.phyloP100way.bw
wget ${BASE_URL}/hg38.phyloP30way.bw
wget ${BASE_URL}/mastermind_cited_variants_reference-2022.04.02-grch38.vcf.gz
wget ${BASE_URL}/mastermind_cited_variants_reference-2022.04.02-grch38.vcf.gz.tbi
wget ${BASE_URL}/repeatmasker_hg38.bed.gz
wget ${BASE_URL}/repeatmasker_hg38.bed.gz.tbi
wget ${BASE_URL}/spliceai_scores.raw.indel.hg38.vcf.gz
wget ${BASE_URL}/spliceai_scores.raw.indel.hg38.vcf.gz.tbi
wget ${BASE_URL}/spliceai_scores.raw.snv.hg38.vcf.gz
wget ${BASE_URL}/spliceai_scores.raw.snv.hg38.vcf.gz.tbi

cd ../..
