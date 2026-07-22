#!/bin/bash

# Tested on Ubuntu 20, 22, 24, 26

set -e # Fail on error

# Install packages that vary in name/availability between Ubuntu versions, one at a
# time, warning (but not aborting the deploy) when a package has no candidate. Accepts
# "pkgA||pkgB" to try an alternate name (e.g. a rename across releases) before warning.
apt_install_optional() {
    for spec in "$@"; do
        installed=false
        # shellcheck disable=SC2086
        IFS='|' read -ra names <<< "$spec"
        for pkg in "${names[@]}"; do
            [ -z "$pkg" ] && continue
            if apt-get install -y "$pkg"; then
                installed=true
                break
            fi
        done
        if [ "$installed" = false ]; then
            echo "WARNING: could not install '$spec' - no candidate on this Ubuntu version, skipping" >&2
        fi
    done
}

echo "Installing packages via apt-get"
apt-get update
apt-get install -y git-core python3-dev python3-pip python3-gdal gfortran libopenblas-dev liblapack-dev zlib1g-dev bedtools bcftools samtools libbz2-dev liblzma-dev libssl-dev libcurl4-openssl-dev curl
apt-get install -y libblas-dev libpng-dev libjpeg-dev libpq-dev libxft-dev zlib1g-dev libxml2-dev libxslt1-dev rustc
# These have churned across Ubuntu versions - libatlas-base-dev was dropped in 26 (libopenblas-dev/liblapack-dev
# above provide optimized BLAS/LAPACK), and libfreetype6-dev was renamed to libfreetype-dev. Install soft so a
# future rename/removal warns instead of aborting the deploy.
apt_install_optional "libfreetype-dev||libfreetype6-dev" "libatlas-base-dev"
apt-get install -y nginx redis-server rabbitmq-server postgresql postgresql-contrib libpq-dev
# actually only needed for development
apt-get install -y sassc

# Stuff required to install on shariant server, might just be moving dependencies
apt-get install -y pkg-config

# Installing Python now done elsewhere so you can do in virtual env