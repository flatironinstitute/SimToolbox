#! /bin/bash

export USER_LOCAL=$HOME/local
export SYSTEM_LOCAL=/usr/local

#cmake \
  #-D CMAKE_CXX_COMPILER=mpicxx \
  #-D CMAKE_C_COMPILER=mpicc \
  #-D CMAKE_BUILD_TYPE=RelWithDebInfo \
  #-D Eigen3_DIR="${USER_LOCAL}/share/eigen3/cmake" \
  #-D TRNG_INCLUDE_DIR="${USER_LOCAL}/include" \
  #-D TRNG_LIBRARY="${USER_LOCAL}/lib/libtrng4.a" \
  #-D Trilinos_DIR="${USER_LOCAL}/lib/cmake/Trilinos" \
  #-D yaml-cpp_DIR="${USER_LOCAL}/lib64/cmake/yaml-cpp" \
#../

cmake \
  -D CMAKE_CXX_COMPILER=mpicxx \
  -D CMAKE_C_COMPILER=mpicc \
  -D ENABLE_TEST=ON \
  -D CMAKE_BUILD_TYPE=Debug \
  -D CMAKE_CXX_FLAGS="-O0 -march=broadwell -DDEBUG" \
  -D SFTPATH="${HOME}/local/" \
  -D Eigen3_DIR="${SFTPATH}/share/eigen3/cmake" \
  -D TRNG_INCLUDE_DIR="${USER_LOCAL}/include" \
  -D TRNG_LIBRARY="${USER_LOCAL}/lib/libtrng4.a" \
../
#-D yaml-cpp_DIR="${USER_LOCAL}/lib64/cmake/yaml-cpp" \
