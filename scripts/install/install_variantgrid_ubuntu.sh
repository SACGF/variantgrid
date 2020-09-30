#!/bin/bash

set -e # Fail on error

# We can wipe this once we move this to public servers
if [ -z ${GITHUB_USER} ]; then
    echo "Please set the variable GITHUB_USER used to download variantgrid and bioinformatics projects"
    exit 1; 
fi
if [ -z ${GITHUB_PASSWORD} ]; then
    echo "Please set the variable GITHUB_PASSWORD used to download variantgrid and bioinformatics projects"
    exit 1; 
fi

if [ -z ${VG_INSTALL_DIR} ]; then
    echo "Please set the variable VG_INSTALL_DIR - where you want to install variantgrid, eg 'export VG_INSTALL_DIR=/opt/variantgrid'"
    exit 1; 
fi

if [ -z ${SYSTEM_VARIANTGRID_USER} ]; then
    echo "Please set the variable SYSTEM_VARIANTGRID_USER - either your username (if installing a dev environment) or the user VG services will run as."
    exit 1; 
fi


export VG_LOG_DIR=/var/log/variantgrid

apt-get install -y git-core

# Make user
id -u ${SYSTEM_VARIANTGRID_USER} &>/dev/null || useradd ${SYSTEM_VARIANTGRID_USER} --create-home
mkdir -p ${VG_INSTALL_DIR}
chown ${SYSTEM_VARIANTGRID_USER} ${VG_INSTALL_DIR}
chgrp ${SYSTEM_VARIANTGRID_USER} ${VG_INSTALL_DIR}

echo "##################"
echo "# Code / LIBRARIES"
echo "##################"

mkdir -p ${VG_LOG_DIR}
chown ${SYSTEM_VARIANTGRID_USER} ${VG_LOG_DIR}

cd ${VG_INSTALL_DIR}
if [ ! -e ${VG_INSTALL_DIR}/.git ]; then
    # Keep files owned by variantgrid user
    su ${SYSTEM_VARIANTGRID_USER} -c 'git clone https://${GITHUB_USER}:${GITHUB_PASSWORD}@github.com/sacgf/variantgrid.git .'
fi

./scripts/install/ubuntu_18_install_dependencies.sh

echo "Download Natural Language Tool Kit data to help us parse sentences."
su ${SYSTEM_VARIANTGRID_USER} -c 'python3 -m nltk.downloader punkt'
su ${SYSTEM_VARIANTGRID_USER} -c 'python3 -m nltk.downloader averaged_perceptron_tagger'


##########
# Database
##########
# Creates user snpdb and database snpdb
su postgres -c "psql < ${VG_INSTALL_DIR}/dbscripts/pgsql_database_create.sql"

echo "Creating server settings file"
su ${SYSTEM_VARIANTGRID_USER} -c 'cp ${VG_INSTALL_DIR}/variantgrid/settings/example_settings.py ${VG_INSTALL_DIR}/variantgrid/settings/$(hostname).py'

# Perform Django migrations to empty database
# Sometimes the dependencies for migrations fail, which is why we explicitly run sites and auth first.
su ${SYSTEM_VARIANTGRID_USER} -c "python3 ${VG_INSTALL_DIR}/manage.py migrate auth"
su ${SYSTEM_VARIANTGRID_USER} -c "python3 ${VG_INSTALL_DIR}/manage.py migrate sites"

echo "Performing VG migrations" 

su ${SYSTEM_VARIANTGRID_USER} -c "python3 ${VG_INSTALL_DIR}/manage.py migrate"

# if you get the below error you will manually need to delete a row from django_sites
# there might be an entry that's simply "2, example.com, example.com"
#
# Applying snpdb.0002_initial_data...Traceback (most recent call last):
#  File "/usr/local/lib/python3.6/dist-packages/django/db/backends/utils.py", line 85, in _execute
#    return self.cursor.execute(sql, params)
#psycopg2.IntegrityError: duplicate key value violates unique constraint "django_site_pkey"
# DETAIL:  Key (id)=(2) already exists.

su ${SYSTEM_VARIANTGRID_USER} -c "python3 ${VG_INSTALL_DIR}/manage.py collectstatic"

echo
echo "Running createsuperuser"
echo

su ${SYSTEM_VARIANTGRID_USER} -c "python3 ${VG_INSTALL_DIR}/manage.py createsuperuser"

echo 
echo "Please follow instructions at https://bitbucket.org/sacgf/variantgrid/wiki/Annotation%20Setup"
echo
echo "Your system settings file is: ${VG_INSTALL_DIR}/variantgrid/settings/$(hostname).py"