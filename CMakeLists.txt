cmake_minimum_required(VERSION 3.10) # same requirement as trilinos
project(SimToolbox)

# global settings
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(SFTPATH $ENV{HOME}/local)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
# this does not change the compiler to use mpi compiler, do cmake -D
# CMAKE_CXX_COMPILER=mpicxx linking executable with mpi
# target_link_libraries(exe PRIVATE OpenMP::OpenMP_CXX)
# target_link_libraries(exe PRIVATE MPI::MPI_CXX)

# Get VTK
find_package(VTK REQUIRED PATHS ${SFTPATH} NO_DEFAULT_PATH)
if(VTK_FOUND)
  # add_definitions("-D USE_VTK")
  message("Found VTK at ${VTK_DIR}")
endif()

find_package(TRNG REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(yaml-cpp REQUIRED)
message("   YAML_CPP_INCLUDE_DIRS = ${YAML_CPP_INCLUDE_DIR}")
message("   YAML_CPP_LIBRARY_DIR = ${YAML_CPP_LIBRARY_DIR}")
message("   YAML_CPP_LIBRARIES = ${YAML_CPP_LIBRARIES}")
set(YAML_CPP_LIBRARIES "${YAML_CPP_LIBRARIES}::${YAML_CPP_LIBRARIES}")


# Get Trilinos as one entity
find_package(Trilinos REQUIRED)
# clean up the too-long list
list(REMOVE_DUPLICATES Trilinos_INCLUDE_DIRS)
list(REMOVE_DUPLICATES Trilinos_TPL_INCLUDE_DIRS)
list(REVERSE Trilinos_LIBRARIES)
list(REMOVE_DUPLICATES Trilinos_LIBRARIES)
list(REVERSE Trilinos_LIBRARIES)
list(REVERSE Trilinos_TPL_LIBRARIES)
list(REMOVE_DUPLICATES Trilinos_TPL_LIBRARIES)
list(REVERSE Trilinos_TPL_LIBRARIES)
# Echo trilinos build info just for fun
message("\nFound Trilinos!  Here are the details: ")
message("   Trilinos_DIR = ${Trilinos_DIR}")
message("   Trilinos_VERSION = ${Trilinos_VERSION}")
message("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
message("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
message("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
message("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
message("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
message("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
message("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
message("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
message("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
message("End of Trilinos details\n")

message(
  STATUS
    "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS"
)
message(STATUS "SFTPATH: ${SFTPATH}")

enable_testing()

add_subdirectory(Util)
add_subdirectory(Trilinos)
add_subdirectory(MPI)
add_subdirectory(Collision)
add_subdirectory(Boundary)
add_subdirectory(Constraint)
add_subdirectory(Sylinder)
