#!/bin/bash
# Deploy upgrader - see `python3 manage.py upgrade --help` for the modes (--quick, --auto-manage, --steps)

set -e

VG_DIR=$(dirname "${BASH_SOURCE[0]}")/..

# Install before Django starts, in case of a manual git pull. A pull from the upgrader installs what it brings
"${VG_DIR}/scripts/install_requirements.sh"

exec python3 "${VG_DIR}/manage.py" upgrade "$@"
