#!/bin/bash

# This downloads the files we can get straight from the web

DISK_FREE_CWD=$(df -Ph . | tail -1 | awk '{print $4}')

echo "This will take a lot of space - make sure you're have a lot of disk space!"
echo "You are in $(pwd), on a drive with ${DISK_FREE_CWD} free."

samtools --version &> /dev/null || { echo >&2 "Samtools not installed!"; exit 1; }

echo "Downloading GRCh38 Fasta"
mkdir -p fasta
cd fasta
# Need to bgzip see https://asia.ensembl.org/info/docs/tools/vep/script/vep_cache.html#fasta
FASTA_FILE=GCF_000001405.39_GRCh38.p13_genomic.fna.gz
if [ ! -e ${FASTA_FILE} ]; then
  wget --quiet -O - https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz | gzip -d | bgzip > ${FASTA_FILE}
  samtools faidx ${FASTA_FILE}
fi
cd ..

mkdir -p annotation_data/GRCh38
cd annotation_data/GRCh38

echo "Conservation"
if [ ! -e hg38.phastCons100way.bw ]; then
  wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw
fi

if [ ! -e hg38.phastCons30way.bw ]; then
  wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons30way/hg38.phastCons30way.bw
fi

if [ ! -e hg38.phyloP100way.bw ]; then
  wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw
fi

if [ ! -e hg38.phyloP30way.bw ]; then
  wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP30way/hg38.phyloP30way.bw
fi

echo "Population databases"
if [ ! -e TOPMED_GRCh38_20180418.vcf.gz ]; then
  wget ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/TOPMED_GRCh38_20180418.vcf.gz
fi

if [ ! -e TOPMED_GRCh38_20180418.vcf.gz.tbi ]; then
  wget ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/TOPMED_GRCh38_20180418.vcf.gz.tbi
fi

if [ ! -e UK10K_COHORT.20160215.sites.GRCh38.vcf.gz ]; then
  wget ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/UK10K_COHORT.20160215.sites.GRCh38.vcf.gz
fi

if [ ! -e UK10K_COHORT.20160215.sites.GRCh38.vcf.gz.tbi ]; then
  wget ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/UK10K_COHORT.20160215.sites.GRCh38.vcf.gz.tbi
fi

cd ../..