export VARIANTGRID_DIR=$(dirname $0)/../..
export DJANGO_SETTINGS_MODULE=variantgrid.settings

pylint --load-plugins pylint_django --rcfile=${VARIANTGRID_DIR}/config/pylint3.rc analysis annotation classification eventlog flags genes library manual ontology pathtests patients pedigree scripts seqauto snpdb sync uicore upload variantgrid variantopedia > lint.txt
