FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get upgrade -y 

ENV DEPS="cmake ninja-build git ca-certificates lsb-release wget software-properties-common gnupg libtbb-dev" 
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y ${DEPS}

# install clang++-17
RUN bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)"

ENV CXX="/usr/bin/clang++-17"

# clone parsb
# clone, configure, and build abseil-cpp
RUN cd ~ \
&& git clone --recurse-submodules https://github.com/ubc-systopia/parsb.git \
&& cd ~/parsb \
&& mkdir install && cd install \
&& git clone https://github.com/abseil/abseil-cpp.git \
&& cd ~/parsb/install/abseil-cpp/ \
&& mkdir build && cd build \
&& cmake .. \
&& cmake --build . --target all

# configure, build parsb
ENV PARSB_ROOT_DIR="~/parsb"
ENV BSIZE=1024
ENV BWIDTH=65536
ENV TIME=0

ENV CMAKE_C_COMPILER="/usr/bin/clang-17"
ENV CMAKE_CXX_COMPILER="/usr/bin/clang++-17"

RUN cd ~ \
&& cd parsb \
&& cmake -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER} \
    -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER} \
    -DBSIZE:STRING=${BSIZE} \
    -DBWIDTH:STRING=${BWIDTH} \
    -DTIME:STRING=${TIME} \
    -S${PARSB_ROOT_DIR} \
    -B${PARSB_ROOT_DIR}/build \
    -G Ninja \
&& cmake --build ${PARSB_ROOT_DIR}/build --config Release --target parsb --