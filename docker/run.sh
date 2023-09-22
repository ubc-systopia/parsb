#!/bin/bash

# configure, build, and run parsb

export PARSB_ROOT_DIR="/root/parsb"
export CMAKE_C_COMPILER="/usr/bin/clang-17"
export CMAKE_CXX_COMPILER="/usr/bin/clang++-17"

# run
${PARSB_ROOT_DIR}/build/parsb \
    -f ${GRAPH_PATH} \
    -s \
    -o ${OUTPUT_PATH} \
    -p ${PERCENT} \
    -t ${NUM_THREADS}