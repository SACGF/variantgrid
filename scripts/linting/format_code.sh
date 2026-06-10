#!/bin/bash

VG_DIR=$(dirname $0)/../..

# Ruff replaces the old autopep8 pass. Safe auto-fixes only: whitespace,
# unused imports, import ordering, modern-syntax. Config: ruff.toml
ruff check "${VG_DIR}" --fix

# Full black-style reformatting is opt-in (large diff) - uncomment to adopt:
# ruff format "${VG_DIR}"
