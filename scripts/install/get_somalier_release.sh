#!/bin/sh

# Downloads somalier and the data it needs into the current directory - run it in whatever
# SOMALIER["annotation_base_dir"] points at, then check it with 'manage.py deployment_check'
#
# v0.3.5 fixed the allele order bug (brentp/somalier#163) by taking the sites VCF at relate time, so a
# deployment on it wants SOMALIER["compensate_allele_order"] = False and its extracts rebuilt with
# 'manage.py somalier_existing_vcfs --clear'

set -eu

SOMALIER_VERSION="0.3.5"
SOMALIER_BINARY="somalier-${SOMALIER_VERSION}"

command -v wget > /dev/null || { echo "wget is required" >&2; exit 1; }

# Anything already downloaded is left alone, so a re-run only fetches what's missing. The partial
# name means a download killed half way isn't mistaken for a complete file next time
download() {
    if [ -s "$2" ]; then
        echo "Already have $2 - skipping"
    else
        echo "Downloading $2"
        wget --output-document "$2.partial" "$1"
        mv "$2.partial" "$2"
    fi
}

download "https://github.com/brentp/somalier/releases/download/v${SOMALIER_VERSION}/somalier" "${SOMALIER_BINARY}"
chmod a+x "${SOMALIER_BINARY}"
if ! "./${SOMALIER_BINARY}" 2>&1 | head -1 | grep -q "${SOMALIER_VERSION}"; then
    echo "'${SOMALIER_BINARY}' doesn't run and report version ${SOMALIER_VERSION}" >&2
    exit 1
fi
# settings.SOMALIER["annotation"]["command"] is 'somalier', so an upgrade is re-running this script
ln -sfn "${SOMALIER_BINARY}" somalier

# Sites files - build specific, and the filenames have to match settings.SOMALIER["annotation"]["sites"]
download "https://github.com/brentp/somalier/files/3412454/sites.hg38.nochr.vcf.gz" "sites.hg38.nochr.vcf.gz"
download "https://github.com/brentp/somalier/files/3412455/sites.GRCh37.vcf.gz" "sites.GRCh37.vcf.gz"
download "https://github.com/brentp/somalier/files/9954286/sites.chm13v2.T2T.vcf.gz" "sites.chm13v2.T2T.vcf.gz"
# The chr-prefixed builds, which no deployment uses:
# download "https://github.com/brentp/somalier/files/3412453/sites.hg19.vcf.gz" "sites.hg19.vcf.gz"
# download "https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz" "sites.hg38.vcf.gz"

# Ancestry: the labels, and the 2,504 1kg samples they label
download "https://raw.githubusercontent.com/brentp/somalier/master/scripts/ancestry-labels-1kg.tsv" "ancestry-labels-1kg.tsv"
download "https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz?download=1" "1kg.somalier.tar.gz"
if [ -d 1kg-somalier ]; then
    echo "Already have 1kg-somalier/ - skipping"
else
    tar xfz 1kg.somalier.tar.gz
    chmod -R a+r 1kg-somalier
fi

echo "Installed $(./somalier 2>&1 | head -1) in $(pwd)"
