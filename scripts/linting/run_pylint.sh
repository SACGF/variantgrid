export VARIANTGRID_DIR=$(dirname $0)/../..

pylint --rcfile=${VARIANTGRID_DIR}/config/pylint3.rc analysis annotation auth eventlog expression genes library pathtests patients pedigree sapath scripts seqauto snpdb upload variantclassification variantdetails variantgrid variantopedia
