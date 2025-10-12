#!/bin/bash

HOST="variantgrid@frgeneseq02"

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 script [args...]" >&2
  exit 1
fi

# This relies on there being a .ssh key setup in variantgrid@frgeneseq02
# And it being automatically applied in ~/.ssh/config (on the VM)

# Quote the script path and forward args as-is
ssh -C "$HOST" "bash -s" -- < "$1" "${@:2}"
