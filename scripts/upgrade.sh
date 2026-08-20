#!/bin/bash

set -e

VG_DIR=$(dirname ${BASH_SOURCE[0]})
export PYTHONPATH=${PYTHONPATH}:${VG_DIR};

# The migrator's first step calls manage.py, which imports settings -> INSTALLED_APPS -> every 3rd
# party package, so the current requirements have to be installed before it starts
"${VG_DIR}/install_requirements.sh"

python3 "./scripts/migrator/migrator.py" $1
