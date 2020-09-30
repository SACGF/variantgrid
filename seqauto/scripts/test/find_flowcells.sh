#!/bin/bash

TEST_DATA_DIR=$(dirname ${BASH_SOURCE[0]})/../../test_data
TEST_DATA_DIR=$(cd "${TEST_DATA_DIR}"; pwd) # Absolute path

find ${TEST_DATA_DIR}/unaligned -maxdepth 3 -name "RTAComplete.txt" -exec dirname {} \;

