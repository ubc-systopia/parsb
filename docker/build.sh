export PARSB_ROOT_DIR="/root/parsb"
export CMAKE_C_COMPILER="/usr/bin/clang-17"
export CMAKE_CXX_COMPILER="/usr/bin/clang++-17"

cd ${PARSB_ROOT_DIR}
mkdir install 
cd install 
git clone https://github.com/abseil/abseil-cpp.git 
cd ${PARSB_ROOT_DIR}/install/abseil-cpp/ 
mkdir build 
cd build 
cmake .. 
cmake --build . --target all

# configure, build parsb
export BSIZE=1024
export BWIDTH=65536
export TIME=0

cd ${PARSB_ROOT_DIR}
cmake -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER} \
    -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER} \
    -DBSIZE:STRING=${BSIZE} \
    -DBWIDTH:STRING=${BWIDTH} \
    -DTIME:STRING=${TIME} \
    -S${PARSB_ROOT_DIR} \
    -B${PARSB_ROOT_DIR}/build \
    -G Ninja \

cmake --build ${PARSB_ROOT_DIR}/build --config Release --target parsb --