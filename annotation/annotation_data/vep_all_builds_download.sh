#!/bin/bash

FILENAME=fordownload.tar.gz

set -e

if [[ ! -e annotation_data/all_builds/maxentscan ]]; then
  mkdir -p annotation_data/all_builds
  cd annotation_data/all_builds
  echo "Downloading MaxEntScan data (build independent)"
  wget http://hollywood.mit.edu/burgelab/maxent/download/${FILENAME}
  tar xvfz ${FILENAME}
  mv fordownload maxentscan
  rm ${FILENAME}
  cd ../..
fi
