FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get upgrade -y 

# libtbb-dev required for ips4o
ENV DEPS="cmake ninja-build git ca-certificates lsb-release wget software-properties-common gnupg libtbb-dev" 
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y ${DEPS}

# install clang++-17
RUN bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)"

# install omp
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y libomp-17-dev

