#!/bin/bash

# This downloads the files we can get straight from the web

DISK_FREE_CWD=$(df -Ph . | tail -1 | awk '{print $4}')

echo "This will take a lot of space - make sure you're have a lot of disk space!"
echo "You are in $(pwd), on a drive with ${DISK_FREE_CWD} free."


echo "Fasta"
# Need to bgzip see https://asia.ensembl.org/info/docs/tools/vep/script/vep_cache.html#fasta
wget --quiet -O - https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz | gzip -d | bgzip > GCF_000001405.39_GRCh38.p13_genomic.fna.gz
samtools faidx GCF_000001405.39_GRCh38.p13_genomic.fna.gz

echo "VEP Cache"
wget ftp://ftp.ensembl.org/pub/release-97/variation/indexed_vep_cache/homo_sapiens_vep_100_GRCh38.tar.gz

echo "Conservation"
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons30way/hg38.phastCons30way.bw

wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP30way/hg38.phyloP30way.bw

echo "Population databases"
wget ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/TOPMED_GRCh38_20180418.vcf.gz
wget ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/TOPMED_GRCh38_20180418.vcf.gz.tbi

wget ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/UK10K_COHORT.20160215.sites.GRCh38.vcf.gz
wget ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/UK10K_COHORT.20160215.sites.GRCh38.vcf.gz.tbi
