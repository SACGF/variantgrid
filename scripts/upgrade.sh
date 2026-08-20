#!/bin/bash

set -e

VG_DIR=$(dirname ${BASH_SOURCE[0]})
export PYTHONPATH=${PYTHONPATH}:${VG_DIR};

# Install now, in case user did manual git pull and it will fail running migrator
"${VG_DIR}/install_requirements.sh"

python3 "./scripts/migrator/migrator.py" $1
