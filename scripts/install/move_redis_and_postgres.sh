#!/bin/bash

set -e # Fail on error

# We can wipe this once we move this to public servers
if [ -z ${BITBUCKET_USER} ]; then
    echo "Please set the variable BITBUCKET_USER used to download variantgrid and bioinformatics projects"
    exit 1; 
fi
if [ -z ${BITBUCKET_PASSWORD} ]; then
    echo "Please set the variable BITBUCKET_PASSWORD used to download variantgrid and bioinformatics projects"
    exit 1; 
fi

if [ -z ${VG_INSTALL_DIR} ]; then
    echo "Please set the variable VG_INSTALL_DIR - where you want to install variantgrid, eg 'export VG_INSTALL_DIR=/opt/variantgrid'"
    exit 1; 
fi

if [ -z ${ANNOTATION_BASEDIR} ]; then
    echo "Please set the variable ANNOTATION_BASEDIR - where you want to store the annotations"
    exit 1; 
fi

if [ -z ${SYSTEM_VARIANTGRID_USER} ]; then
    echo "Please set the variable SYSTEM_VARIANTGRID_USER - either your username (if installing a dev environment) or the user VG services will run as."
    exit 1; 
fi

export POSTGRES_DIR=/var/lib/postgresql/9.5
export POSTGRES_MAIN_DIR=${POSTGRES_DIR}/main
export POSTGRES_PID=${POSTGRES_MAIN_DIR}/postmaster.pid
export DATABASE_DIR=${BASE_DIR}/variantgrid_database
export REDIS_OLD_DATABASE_DIR=/var/lib/redis/
export REDIS_NEW_DATABASE_DIR=${BASE_DIR}/redis_database
export REDIS_DB_FILE_NAME=dump.rdb
export REDIS_ORIG_DATABASE_FILE=${REDIS_OLD_DATABASE_DIR}/${REDIS_DB_FILE_NAME}
export REDIS_NEW_DATABASE_FILE=${REDIS_NEW_DATABASE_DIR}/${REDIS_DB_FILE_NAME}
export REDIS_PID=/var/run/redis/redis-server.pid
export REDIS_SERVICE_SCRIPT=/etc/systemd/system/redis.service
export SACGF_DIR=${BASE_DIR}/sacgf
export IVAT_DIR=${SACGF_DIR}/ivat
export REFERENCE_DIR=${SACGF_DIR}/reference
export ANNOTATION_SCRATCH_DIR=${BASE_DIR}/annotation_scratch

function wait_for_file {
    file_name=$1;
    want_file_to_exist=$2; # true/false
    max_tries=$3;
    
    for ((i=0; i<$max_tries ; i++ )); do
        set +e
        test -e ${file_name};
        file_exists=$?
        set -e
        
        echo "${file_name} exist test returned = ${file_exists}";
        if [[ "$want_file_to_exist" = true && $file_exists == 0 ]]; then
            echo "both exist";
            return 0;
        fi;
        if [[ "$want_file_to_exist" = false && $file_exists != 0 ]]; then
            echo "both don't exist";
            return 0;
        fi;
        sleep 1;
    done
    
    if [ "$want_file_to_exist" = true ]; then
        echo "Couldn't find '${file_name}' in ${max_tries} seconds";
    else
        echo "File '${file_name}' didn't disappear in ${max_tries} seconds";
    fi
    return 1;
}


echo "Making sure that Postgres is installed / running"
wait_for_file ${POSTGRES_PID} true 10

if [[ -L "${POSTGRES_MAIN_DIR}" && -d "${POSTGRES_MAIN_DIR}" ]]; then
    echo "Postgres db dir ${POSTGRES_MAIN_DIR} is already a symbolic link";
else
    # The default postgres database is on / so we want it on the bigger drive
    echo "Moving Postgres database from ${POSTGRES_MAIN_DIR} => ${DATABASE_DIR}"

    service postgresql stop

    # Make sure Postgres has stopped
    wait_for_file ${POSTGRES_PID} false 10

    # Move data
    cp -aRv ${POSTGRES_MAIN_DIR} ${DATABASE_DIR}
    mv ${POSTGRES_MAIN_DIR} ${POSTGRES_MAIN_DIR}.OLD
    su postgres -c 'ln -s ${DATABASE_DIR} ${POSTGRES_MAIN_DIR}'

    # Make sure Postgres is running (from installation)
    service postgresql start
    wait_for_file ${POSTGRES_PID} true 10
fi


if [[ ! -z ${MOVE_REDIS_DATABASE} ]]; then
    echo "Making sure that Redis is installed / running"
    wait_for_file ${REDIS_PID} true 10
    
    if [[ -L "${REDIS_ORIG_DATABASE_FILE}" ]]; then
        echo "Redis DB: ${REDIS_ORIG_DATABASE_FILE} is already a symbolic link";
    else
        # The default redis database is on / so we want it on the bigger drive
        echo "Moving Redis database from ${REDIS_ORIG_DATABASE_FILE} => ${REDIS_NEW_DATABASE_DIR}"
    
        service redis stop
    
        # Make sure Redis has stopped
        wait_for_file ${REDIS_PID} false 10
    
        # Modify the service script
        echo "Added on $(date) by $(whoami) via $0" >> ${REDIS_SERVICE_SCRIPT}
        echo "ReadWriteDirectories=-${REDIS_NEW_DATABASE_DIR}" >> ${REDIS_SERVICE_SCRIPT}
    
        # Move data
        mkdir -p REDIS_NEW_DATABASE_DIR
        chown redis ${REDIS_NEW_DATABASE_DIR}
        mv ${REDIS_ORIG_DATABASE_FILE} ${REDIS_NEW_DATABASE_DIR}
        su redis -c 'ln -s ${REDIS_NEW_DATABASE_FILE} ${REDIS_ORIG_DATABASE_FILE}'
    
        # Make sure Redis is running (from installation)
        service redis start
        wait_for_file ${REDIS_PID} true 10
    fi
fi
