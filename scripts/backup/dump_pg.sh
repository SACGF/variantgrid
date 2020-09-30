#!/usr/bin/env bash
 
# script dumps a database on localhost using pg_dump
# create different env files with variables for script:
#
# env file contains:
# export DATABASE_NAME="snpdb" - database to dump
# export BACKUP_DIR=... - full path to backups main directory

envfile=$1
 
function echo_log {
  echo "$(date +"%Y-%m-%d %H:%M:%S"): $1"
}

homedir=$( dirname "${BASH_SOURCE[0]}")
cd $homedir
 
[ -z "$envfile" ] && echo "USAGE: ./$(basename $0) env_file" && exit 1
 
if [ ! -f ${envfile} ]; then
  echo_log "cannot find env file: ${envfile}"
  exit 1
fi
 
. ${envfile}

[ -z "DATABASE_NAME" ] && echo "$(basename $0): ERROR: variable DATABASE_NAME empty or undefined" && exit 1
[ -z "BACKUP_DIR" ] && echo "$(basename $0): ERROR: variable BACKUP_DIR empty or undefined" && exit 1
if [ ! -d ${BACKUP_DIR} ]; then
  echo_log "cannot find backup main directory: ${BACKUP_DIR}"
  #exit 1
fi

DATE_STAMP=$(date +"%Y%m%d")
BACKUP_FILENAME="${BACKUP_DIR}/$(hostname).${DATABASE_NAME}.${DATE_STAMP}_dump.sql.gz"

if [ -e ${BACKUP_FILENAME} ]; then
  echo_log "File '${BACKUP_FILENAME}' already exists!"
  exit 1
fi

echo_log "executing: pg_dump ${DATABASE_NAME} | gzip > ${BACKUP_FILENAME};"
pg_dump ${DATABASE_NAME} | gzip > ${BACKUP_FILENAME};
if [[ $? -ne 0 ]]; then
  echo_log "$(basename $0): error in pg_dump"
  exit 1
fi

echo_log "DUMP COMPLETED"
