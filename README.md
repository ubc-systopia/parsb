# parsb
A parallel implementation of the SlashBurn vertex ordering algorithm using ips4o sort, Afforest, and Spray Block Reductions

# Docker
## Requirements
| Name   | Version   |
| ------ | --------- |
| Docker | >= 24.0.6 |
| Python | >= 3.10   |

Use the [atrostan/parsb](https://hub.docker.com/repository/docker/atrostan/parsb/general) Docker image at  to run `parsb`.

1. Create a virtual environment to run the docker container   
(Required packages are listed in `./docker/requirements.txt`)  
   ```bash
   python -m venv ./docker/venv
   source ./docker/venv/bin/activate
   pip install -r requirements.txt
   ```
3. Run `parsb`  
   e.g.
   ```bash
   python docker/run.py \
      --graph-path "./data/graphs/librec-ciaodvd-trust.net" \
      --output-path "./data/graphs/sb"  \
      --p 0.005 \
      --num-threads 8 \
      --block-width 65536 
   ```
   This will:  
     1. pull docker image from [atrostan/parsb](https://hub.docker.com/repository/docker/atrostan/parsb/general) (if not already pulled)
     2. (re)compile `parsb` using the input arguments
     3. run `parsb` on the given edgelist

# Compilation and Build
Follow these instructions if you wish to configure, build, and run `parsb` locally.
## Cloning

`git clone --recurse-submodules https://github.com/ubc-systopia/parsb.git`

## Requirements and Dependencies
`parsb` was developed and tested using Ubuntu clang version 17.0.0.

| Name   | Version   |
| ------ | --------- |
| CMake  | >= 3.16   |
| OpenMP | >= 5.1    |
| Ninja  | >= 1.10.1 |
| abseil | latest    |


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