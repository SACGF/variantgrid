# Full report to lint.txt (what `scripts/vg health` reads), then the import contracts and cycles.
# --errors-only [pylint args...]: only the E-class messages, to stdout - what CI runs (agent-maps.yml)
export VARIANTGRID_DIR=$(dirname $0)/../..
export DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE:-variantgrid.settings}
# pylint_model_utils lives next to this script
export PYTHONPATH=$(dirname $0)${PYTHONPATH:+:${PYTHONPATH}}

PACKAGES="analysis annotation classification eventlog flags genes library manual mme ontology pathtests patients pedigree scripts seqauto snpdb sync uicore upload variantgrid variantopedia"
PYLINT_ARGS="--load-plugins pylint_django,pylint_model_utils --rcfile=${VARIANTGRID_DIR}/config/pylint3.rc"

if [ "$1" = "--errors-only" ]; then
    shift
    exec pylint ${PYLINT_ARGS} --errors-only "$@" ${PACKAGES}
fi

pylint ${PYLINT_ARGS} ${PACKAGES} > lint.txt

(cd ${VARIANTGRID_DIR} && uvx --from import-linter lint-imports --config .importlinter)
${VARIANTGRID_DIR}/scripts/vg imports cycles
