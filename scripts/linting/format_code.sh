#!/bin/bash

VG_DIR=$(dirname $0)/../..

autopep8 ${VG_DIR} --recursive --select=W293,W391,E203,E242,E251,E252,E261,E27,E303,W291,W292,W293,W391 --in-place
