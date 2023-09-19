#!/bin/bash

# configure, build, and run parsb

export PARSB_ROOT_DIR="~/parsb"
export CMAKE_C_COMPILER="/usr/bin/clang-17"
export CMAKE_CXX_COMPILER="/usr/bin/clang++-17"

# configure
cmake \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER} \
    -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER} \
    -DBSIZE:STRING=${BSIZE} \
    -DBWIDTH:STRING=${BWIDTH} \
    -DTIME:STRING=${TIME} \
    -S${PARSB_ROOT_DIR} \
    -B${PARSB_ROOT_DIR}/build \
    -G Ninja

# build
cmake --build ${PARSB_ROOT_DIR}/build --config Release --target parsb --

# run
${PARSB_ROOT_DIR}/build/parsb \
    -f ${GRAPH_PATH} \
    -s \
    -o ${OUTPUT_PATH} \
    -p ${PERCENT} \
    -t ${NUM_THREADS}