# parsb
A parallel implementation of the SlashBurn vertex ordering algorithm using ips4o sort, Afforest, and Spray Block Reductions

# Repository Structure

# Compilation and Build

## Cloning

`git clone --recurse-submodules https://github.com/ubc-systopia/parsb.git`

## Requirements and Dependencies

1. `parsb` was developed and tested using Ubuntu clang version 17.0.0.
2. 
| Name   | Version |
| ------ | ------- |
| CMake  | >= 3.16 |
| OpenMP | >= 5.1  |
| Ninja  | >= 1.10.1  |
| abseil | tmp     |


### Installing abseil
Intructions taken from https://abseil.io/docs/cpp/quickstart-cmake.html


```bash
$ cd <root directory for pasrb>

# clone
mkdir install && cd install
git clone https://github.com/abseil/abseil-cpp.git

# configure
cd abseil-cpp
mkdir build && cd build
cmake ..

# build
cmake --build . --target all
```

## Building Parallel SlashBurn

### Configure
1. ensure submodules are installed
   1. `git submodule update --init --recursive`
2. Compile time params:
   1. `BSIZE`: 
   2. `BWIDTH`: 
   3. `TIME`:

e.g.
```bash

export PARSB_ROOT_DIR=<root directory for parsb>
export BSIZE=1024
export BWIDTH=65536
export TIME=0

cmake \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_C_COMPILER:STRING=/usr/bin/clang \
    -DCMAKE_CXX_COMPILER:STRING=/usr/bin/clang++ \
    -DBSIZE:STRING=${BSIZE} \
    -DBWIDTH:STRING=${BWIDTH} \
    -DTIME:STRING=${TIME} \
    -S${PARSB_ROOT_DIR} \
    -B${PARSB_ROOT_DIR}/build \
    -G Ninja

```

### Build

```bash
cmake --build ${PARSB_ROOT_DIR}/build --config Release --target parsb --
```


# Running Parallel SlashBurn