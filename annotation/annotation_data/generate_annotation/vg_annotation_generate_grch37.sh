#!/bin/bash

# This creates VariantGrid GRCh37 annotation from scratch - it may take a long time

SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")
MERGE_SCRIPT="${SCRIPT_DIR}/merge_wigs_to_bigwig.sh"

echo "Tools"
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigCat
chmod a+x wigToBigWig bigWigCat


# Gnomad is done in own script + processed

echo "Merging phastCons46way wigs -> genome bigWig"
$(MERGE_SCRIPT) hg19.chrom.sizes hg19.46way.phastCons.bw hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons46way/placentalMammals/*.wigFix.gz

echo "Merging phyloP46way wigs -> genome bigWig"
$(MERGE_SCRIPT) hg19.chrom.sizes hg19.phyloP46way.bw hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP46way/placentalMammals/*.wigFix.gz
