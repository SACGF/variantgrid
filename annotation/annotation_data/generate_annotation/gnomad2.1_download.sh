#!/bin/bash

# gnomad v4.0

# Exomes
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
  wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.${chrom}.vcf.bgz
  wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.${chrom}.vcf.bgz.tbi
done

# Genomes
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
  wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.${chrom}.vcf.bgz
  wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.${chrom}.vcf.bgz.tbi
done

# Structural
wget https://gnomad-public-us-east-1.s3.amazonaws.com/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz
wget https://gnomad-public-us-east-1.s3.amazonaws.com/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz.tbi
