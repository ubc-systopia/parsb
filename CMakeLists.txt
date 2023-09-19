# Parallel SlashBurn CMake
cmake_minimum_required(VERSION 3.16)
project(par_slashburn)
set(CMAKE_CXX_STANDARD 17)

# fetch fmt-lib, spray
include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG master
  GIT_PROGRESS TRUE)

add_compile_definitions(__TBB_DYNAMIC_LOAD_ENABLED=0)
FetchContent_MakeAvailable(fmt)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -pedantic -march=native -D_REENTRANT=1 -fopenmp")
message(${CMAKE_CXX_COMPILER_ID})

FetchContent_Declare(
  spray
  GIT_REPOSITORY https://github.com/atrostan/spray.git
  GIT_TAG master)

FetchContent_GetProperties(spray)
if(NOT spray_POPULATED)
  FetchContent_Populate(spray)
endif()

message(spray_SOURCE_DIR="${spray_SOURCE_DIR}")

FetchContent_MakeAvailable(spray)

add_library(spray INTERFACE)
target_include_directories(spray INTERFACE ${spray_SOURCE_DIR}/include)

set(PARSB_HEADERS
    include/benchmark.h
    include/bitmap.h
    include/builder.h
    include/command_line.h
    include/generator.h
    include/graph.h
    include/platform_atomics.h
    include/pvector.h
    include/reader.h
    include/timer.h
    include/writer.h
    include/SlashBurn.h
    include/util.h)

set(IPS4O_USE_OPENMP "true")
set(ONEDPL_PAR_BACKEND "openmp")
set(IPS4O_DISABLE_PARALLEL OFF)
set(ips4o_SOURCE_DIR "./ips4o")

add_subdirectory(${ips4o_SOURCE_DIR})
find_package(TBB REQUIRED)

# include abseil
add_subdirectory(install/abseil-cpp)

add_executable(parsb main.cpp src/SlashBurn.cpp ${PARSB_HEADERS})

target_include_directories(parsb PRIVATE include/)

# compile time constants for spray
target_compile_definitions(parsb PUBLIC BSIZE=${BSIZE})
target_compile_definitions(parsb PUBLIC BWIDTH=${BWIDTH})

# flag --  time subroutines for each iteration of parallel slashburn
target_compile_definitions(parsb PUBLIC TIME=${TIME})

target_link_libraries(parsb absl::flat_hash_set absl::flat_hash_map fmt::fmt
                      ips4o spray)
