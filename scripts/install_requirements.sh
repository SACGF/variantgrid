#!/bin/bash
# Install pinned dependencies. Assumes the venv is already activated.
# Deployments install requirements.txt only. Developers pass --dev to add the lint/type-check tools
# in requirements-dev.txt.
# Set VG_INSTALL_REQUIREMENTS=0 to skip, so a broken install can't hold up an upgrade.

VG_DIR=$(dirname "${BASH_SOURCE[0]}")/..

case "${VG_INSTALL_REQUIREMENTS:-1}" in
    0|no|false)
        echo "VG_INSTALL_REQUIREMENTS=${VG_INSTALL_REQUIREMENTS} - skipping requirements install"
        exit 0
        ;;
esac

REQUIREMENTS=(-r "${VG_DIR}/requirements.txt")
if [[ "$1" == "--dev" ]]; then
    REQUIREMENTS+=(-r "${VG_DIR}/requirements-dev.txt")
fi

# uv.toml (the 7-day cooldown) needs uv 0.9.25+ - older uv fails to parse it
UV_MIN_VERSION=0.9.25

if command -v uv > /dev/null; then
    UV_VERSION=$(uv --version | awk '{print $2}')
    if [[ $(printf '%s\n' "${UV_MIN_VERSION}" "${UV_VERSION}" | sort -V | head -1) != "${UV_MIN_VERSION}" ]]; then
        echo "uv ${UV_VERSION} is older than ${UV_MIN_VERSION}, which uv.toml needs - upgrade it (uv self update, or pip install -U uv if pip installed it)" >&2
        exit 1
    fi
    echo "Installing requirements with uv"
    uv pip install "${REQUIREMENTS[@]}"
else
    echo "uv not found - installing requirements with pip"
    python3 -m pip install --quiet "${REQUIREMENTS[@]}"
fi

STATUS=$?
if [[ ${STATUS} -ne 0 ]]; then
    echo >&2
    echo "Installing requirements failed. To carry on without it (dependency changes in this" >&2
    echo "upgrade will be missing, so migrate/collectstatic may fail) re-run with:" >&2
    echo "    VG_INSTALL_REQUIREMENTS=0 ./scripts/upgrade.sh" >&2
fi
exit ${STATUS}
