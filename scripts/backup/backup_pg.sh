#!/usr/bin/env bash
 
#script makes full backup of postgresql instance using pg_basebackup
#create different env files with variables for script:
#
#env file contains:
#export BACKUP_MAIN_DIR=... - full path to backups main directory
#export BACKUP_WAL_DIR=... - full path to wal main directory
#export BACKUP_PG_IP=... - IP of postgresql instance to be backuped - 127.0.0.1 for local
 
envfile=$1
 
function echo_log {
  echo "$(date +"%Y-%m-%d %H:%M:%S.%N"): $1"
}

homedir=$( dirname "${BASH_SOURCE[0]}")
cd $homedir
 
[ -z "$envfile" ] && echo "USAGE: ./$(basename $0) env_file" && exit 1
 
if [ ! -f ${envfile} ]; then
  echo_log "cannot find env file: ${envfile}"
  exit 1
fi
 
#homedir=$( dirname "${BASH_SOURCE[0]}")
#cd $homedir
. ./${envfile}
 
[ -z "$BACKUP_MAIN_DIR" ] && echo "$(basename $0): ERROR: variable BACKUP_MAIN_DIR empty or undefined" && exit 1
[ -z "$BACKUP_WAL_DIR" ] && echo "$(basename $0): ERROR: variable BACKUP_WAL_DIR empty or undefined" && exit 1
[ -z "$BACKUP_PG_IP" ] && echo "$(basename $0): ERROR: variable BACKUP_PG_IP empty or undefined" && exit 1
 
backupstamp=$(date +"%Y%m%d%H%M%S")
echo_log "======= pg_basebackup ${BACKUP_PG_IP} - ${backupstamp} ======="
backupmaindir=${BACKUP_MAIN_DIR}
if [ ! -d ${backupmaindir} ]; then
  echo_log "cannot find backup main directory: ${backupmaindir}"
  exit 1
fi
echo_log "backup main dir: ${backupmaindir}"
 
mkdir -p ${backupmaindir}/${backupstamp}
if [[ $? -ne 0 ]]; then
  echo_log "$(basename $0): cannot create backup directory - ${backupmaindir}/${backupstamp}"
  exit 1
fi
 
#move existing wal files to archive
backupwaldir=${BACKUP_WAL_DIR}
if [ ! -d ${backupwaldir} ]; then
  echo_log "cannot find backup wal directory: ${backupwaldir}"
  exit 1
fi
echo_log "backup wal dir: ${backupwaldir}"

mkdir -p ${backupwaldir}/${backupstamp}
if [[ $? -ne 0 ]]; then
  echo_log "$(basename $0): cannot create backup wal directory - ${backupwaldir}/${backupstamp}"
  exit 1
fi

find ${backupwaldir} -maxdepth 1 -type f -exec mv {} ${backupwaldir}/${backupstamp} \;
tar -czvf ${backupwaldir}/${backupstamp}/wal_backup.tar -C ${backupwaldir}/${backupstamp} .
find ${backupwaldir}/${backupstamp} -maxdepth 1 -type f ! -name '*.tar' -delete

#Run full backup process
echo_log "backup dir: ${backupmaindir}/${backupstamp}"
pg_basebackup -h ${BACKUP_PG_IP} -D ${backupmaindir}/${backupstamp} -X stream -U postgres -P -F t
if [[ $? -ne 0 ]]; then
  echo_log "$(basename $0): error in pg_basebackup"
  exit 1
fi
 
 
echo_log "BACKUP COMPLETED"
