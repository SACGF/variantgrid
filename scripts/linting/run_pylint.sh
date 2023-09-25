export VARIANTGRID_DIR=$(dirname $0)/../..

pylint --rcfile=${VARIANTGRID_DIR}/config/pylint3.rc analysis annotation auth classification eventlog flags genes library manual ontology pathtests patients pedigree scripts seqauto snpdb sync uicore upload variantgrid variantopedia > lint.txt
