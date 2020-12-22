#!/bin/sh

# Run this where you want data to go...

wget https://github.com/brentp/somalier/releases/download/v0.2.12/somalier
#wget https://github.com/brentp/somalier/files/3412453/sites.hg19.vcf.gz
wget https://github.com/brentp/somalier/files/3412454/sites.hg38.nochr.vcf.gz
wget https://github.com/brentp/somalier/files/3412455/sites.GRCh37.vcf.gz
# wget https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz
wget https://raw.githubusercontent.com/brentp/somalier/master/scripts/ancestry-labels-1kg.tsv
wget https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz?download=1 --output-document=1kg.somalier.tar.gz

chmod a+x somalier
tar xfz 1kg.somalier.tar.gz
