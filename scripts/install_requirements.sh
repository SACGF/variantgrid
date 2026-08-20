#!/bin/bash
# Install pinned dependencies. Assumes the venv is already activated.
# Set VG_INSTALL_REQUIREMENTS=0 to skip, so a broken install can't hold up an upgrade.

VG_DIR=$(dirname "${BASH_SOURCE[0]}")/..

case "${VG_INSTALL_REQUIREMENTS:-1}" in
    0|no|false)
        echo "VG_INSTALL_REQUIREMENTS=${VG_INSTALL_REQUIREMENTS} - skipping requirements install"
        exit 0
        ;;
esac

if command -v uv > /dev/null; then
    echo "Installing requirements with uv"
    uv pip install -r "${VG_DIR}/requirements.txt"
else
    echo "uv not found - installing requirements with pip"
    python3 -m pip install --quiet -r "${VG_DIR}/requirements.txt"
fi

STATUS=$?
if [[ ${STATUS} -ne 0 ]]; then
    echo >&2
    echo "Installing requirements failed. To carry on without it (dependency changes in this" >&2
    echo "upgrade will be missing, so migrate/collectstatic may fail) re-run with:" >&2
    echo "    VG_INSTALL_REQUIREMENTS=0 ./scripts/upgrade.sh" >&2
fi
exit ${STATUS}
