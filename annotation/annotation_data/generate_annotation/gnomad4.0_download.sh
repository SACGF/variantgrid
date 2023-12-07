#!/bin/bash

# gnomad v4.0

# Exomes
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
  wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr${chrom}.vcf.bgz
  wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr${chrom}.vcf.bgz.tbi
done

# Genomes
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
  wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.chr${chrom}.vcf.bgz
  wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.chr${chrom}.vcf.bgz.tbi
done

# Structural
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
  GNOMAD_VCF=gnomad.v4.0.sv.chr${chrom}.vcf.gz
  wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.0/genome_sv/${GNOMAD_VCF}
  wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.0/genome_sv/${GNOMAD_VCF}.tbi
done