#!/bin/bash

DISK_FREE_CWD=$(df -Ph . | tail -1 | awk '{print $4}')
SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")
MERGE_SCRIPT="${SCRIPT_DIR}/merge_wigs_to_bigwig.sh"

echo "This will take a lot of space - make sure you're have a lot of disk space!"
echo "You are in $(pwd), on a drive with ${DISK_FREE_CWD} free."

echo "Fasta"
# Need to bgzip see https://asia.ensembl.org/info/docs/tools/vep/script/vep_cache.html#fasta
wget --quiet -O - https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.fna.gz | gzip -d | bgzip >  GCF_000001405.25_GRCh37.p13_genomic.fna.gz
samtools faidx GCF_000001405.25_GRCh37.p13_genomic.fna.gz

echo "VEP Cache"
wget ftp://ftp.ensembl.org/pub/release-97/variation/indexed_vep_cache/homo_sapiens_vep_100_GRCh37.tar.gz

echo "Tools"
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigCat
chmod a+x wigToBigWig bigWigCat

echo "Conservation"
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons.bw
wget --recursive --no-parent -R "index.html*" ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons46way/placentalMammals/

wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw
wget --recursive --no-parent -R "index.html*" ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP46way/placentalMammals/

echo "Population databases"
wget ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/TOPMED_GRCh37.vcf.gz
wget ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/TOPMED_GRCh37.vcf.gz.tbi

wget ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/UK10K_COHORT.20160215.sites.vcf.gz
wget ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/UK10K_COHORT.20160215.sites.vcf.gz.tbi

# Gnomad is done in own script + processed

echo "Merging phastCons46way wigs -> genome bigWig"
$(MERGE_SCRIPT) hg19.chrom.sizes hg19.46way.phastCons.bw hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons46way/placentalMammals/*.wigFix.gz

echo "Merging phyloP46way wigs -> genome bigWig"
$(MERGE_SCRIPT) hg19.chrom.sizes hg19.phyloP46way.bw hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP46way/placentalMammals/*.wigFix.gz
