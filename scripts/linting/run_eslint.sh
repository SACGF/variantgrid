export VARIANTGRID_DIR=$(dirname $0)/../..

npx eslint ${VARIANTGRID_DIR}/variantgrid/static_files > eslint.txt
