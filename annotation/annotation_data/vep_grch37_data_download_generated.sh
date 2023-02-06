#!/bin/bash

# Files commented out are available on the web (vep_grch37_data_download.sh)

BASE_URL=https://variantgrid.com/download/annotation/VEP/annotation_data/GRCh37

cd annotation_data/GRCh37

wget ${BASE_URL}/CosmicCodingMuts_v95_20211101_grch37.normal.vcf.gz
wget ${BASE_URL}/CosmicCodingMuts_v95_20211101_grch37.normal.vcf.gz.tbi
#wget ${BASE_URL}/TOPMED_GRCh37.vcf.gz
#wget ${BASE_URL}/TOPMED_GRCh37.vcf.gz.tbi
#wget ${BASE_URL}/UK10K_COHORT.20160215.sites.vcf.gz
#wget ${BASE_URL}/UK10K_COHORT.20160215.sites.vcf.gz.tbi
wget ${BASE_URL}/dbNSFP4.3a.grch37.stripped.gz
wget ${BASE_URL}/dbNSFP4.3a.grch37.stripped.gz.tbi
wget ${BASE_URL}/dbscSNV1.1_GRCh37.txt.gz
wget ${BASE_URL}/dbscSNV1.1_GRCh37.txt.gz.tbi
wget ${BASE_URL}/gnomad2.1.1_GRCh37_combined_af.vcf.bgz
wget ${BASE_URL}/gnomad2.1.1_GRCh37_combined_af.vcf.bgz.tbi
wget ${BASE_URL}/gnomad_GRCh37_af_greater_than_5.contigs.vcf.bgz
wget ${BASE_URL}/gnomad_GRCh37_af_greater_than_5.contigs.vcf.bgz.tbi
#wget ${BASE_URL}/hg19.100way.phastCons.bw
#wget ${BASE_URL}/hg19.100way.phyloP100way.bw
#wget ${BASE_URL}/hg19.phastCons46way.placental.bw
#wget ${BASE_URL}/hg19.phyloP46way.placental.bw
wget ${BASE_URL}/mastermind_cited_variants_reference-2022.04.02-grch37.vcf.gz
wget ${BASE_URL}/mastermind_cited_variants_reference-2022.04.02-grch37.vcf.gz.tbi
wget ${BASE_URL}/repeatmasker_hg19.bed.gz
wget ${BASE_URL}/repeatmasker_hg19.bed.gz.tbi
wget ${BASE_URL}/spliceai_scores.raw.indel.hg19.vcf.gz
wget ${BASE_URL}/spliceai_scores.raw.indel.hg19.vcf.gz.tbi
wget ${BASE_URL}/spliceai_scores.raw.snv.hg19.vcf.gz
wget ${BASE_URL}/spliceai_scores.raw.snv.hg19.vcf.gz.tbi


cd ../..
