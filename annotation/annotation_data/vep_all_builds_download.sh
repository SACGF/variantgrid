#!/bin/bash

echo "Downloading MaxEntScan data (build independent)"

mkdir -p annotation_data/all_builds/maxentscan
if [[ ! -e annotation_data/all_builds/maxentscan ]]; then
  cd annotation_data/all_builds
  wget http://hollywood.mit.edu/burgelab/maxent/download/fordownload.tar.gz
  tar xvfz ~/Downloads/fordownload.tar.gz
  mv fordownload maxentscan
  cd ../..
fi
