# This is just a dumping ground, I need to clean it up...

##################################
# MaxEntScan
##################################

wget http://hollywood.mit.edu/burgelab/maxent/download/fordownload.tar.gz
tar xvfz ~/Downloads/fordownload.tar.gz
mv fordownload/ maxentscan


# Download repeat masker bed file:
# UCSC Table browser. "Repeats" group and "Repeatmasker" track. Select output format as BED.
##################################
# RepeatMasker
##################################

# hg19
zgrep -v "#" repeatmasker.bed.gz  | sort -k1,1 -k2,2n -k3,3n -t$'\t' | bgzip -c > repeatmasker_hg19.bed.gz
tabix -p bed repeatmasker_hg19.bed.gz



##################################
# dbNSFP
##################################

echo dbNSFP

# https://asia.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#dbnsfp
head -n1 dbNSFP4.0b2a_variant.chr1 > h

zgrep -h -v ^#chr dbNSFP4.0b2a_variant.chr* | awk '$8 != "."' | sort --temporary-directory=/media/dlawrence/SpinningIron/tmp/ -k8,8 -k9,9n - | cat h - | bgzip -c > dbNSFP4.0b2a.hg19.gz
tabix -s 8 -b 9 -e 9 dbNSFP4.0b2a.hg19.gz

# Hg38:

zgrep -h -v ^#chr dbNSFP4.0b2a_variant.chr* | sort --temporary-directory=/media/dlawrence/SpinningIron/tmp/ -k1,1 -k2,2n - | cat h - | bgzip -c > dbNSFP4.0b2a.hg38.gz


##################################
# dbSNV
##################################


echo dbSNV

#Hg19:
cat dbscSNV1.1.chr* | grep -v ^chr | cat h - | bgzip -c > dbscSNV1.1_GRCh37.txt.gz
tabix -s 1 -b 2 -e 2 -c c dbscSNV1.1_GRCh37.txt.gz

#h38:
cat dbscSNV1.1.chr* | grep -v ^chr | awk '$5 != "."' | sort --temporary-directory=/media/dlawrence/SpinningIron/tmp -k5,5 -k6,6n | cat h - | bgzip -c > dbscSNV1.1_GRCh38.txt.gz

# Post | awk bit above as email to Ensembl or on github
