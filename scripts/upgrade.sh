#!/bin/bash

VG_DIR=$(dirname ${BASH_SOURCE[0]})
export PYTHONPATH=${PYTHONPATH}:${VG_DIR};

python3 "./scripts/migrator/migrator.py" $1
