# A CMake build script for the Fortran part of Ipopt.
# It is separate from the main build script to be able to use
# cmake_add_fortran_subdirectory and build Fortran part with a
# MinGW gfortran on Windows while the rest is built with Visual C++.

cmake_minimum_required(VERSION 2.8)
project(IPOPTFORT C Fortran)

include(FortranCInterface)

# Organize output files.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${IPOPTFORT_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${IPOPTFORT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${IPOPTFORT_BINARY_DIR}/lib)

# Set F77_FUNC & F77_FUNC_ in cache as they are used in the parent project.
set(F77_FUNC "${FortranCInterface_GLOBAL_MACRO}" CACHE INTERNAL "")
set(F77_FUNC_ "${FortranCInterface_GLOBAL__MACRO}" CACHE INTERNAL "")

if (FortranCInterface_GLOBAL_CASE STREQUAL "UPPER")
  add_definitions(-DUPPER)
elseif (FortranCInterface_GLOBAL_SUFFIX STREQUAL "_")
  add_definitions(-DAdd_)
elseif (FortranCInterface_GLOBAL_SUFFIX STREQUAL "__")
  add_definitions(-DAdd__)
endif ()

add_definitions(-DMUMPS_ARITH=MUMPS_ARITH_d)
if (WIN32)
  add_definitions(-DWITHOUT_PTHREAD=1)
endif ()

include_directories(
  ../ThirdParty/Mumps/MUMPS/include ../ThirdParty/Mumps/MUMPS/libseq)

set(mumps_src )
foreach (src
  src/tools_common_mod.F
  src/dmumps_comm_buffer.F
  src/dmumps_struc_def.F
  src/mumps_ooc_common.F
  src/mumps_static_mapping.F
  src/dmumps_ooc_buffer.F
  src/dmumps_load.F
  src/dmumps_ooc.F
  src/dmumps_part1.F
  src/dmumps_part2.F
  src/dmumps_part3.F
  src/dmumps_part4.F
  src/dmumps_part5.F
  src/dmumps_part6.F
  src/dmumps_part7.F
  src/dmumps_part8.F
  src/mumps_part9.F
  src/mumps_sol_es.F
  src/mumps_c.c
  src/mumps_common.c
  src/mumps_orderings.c
  src/mumps_io.c
  src/mumps_io_basic.c
  src/mumps_io_thread.c
  src/mumps_io_err.c
  src/mumps_size.c
  libseq/mpi.f
  libseq/mpic.c
  libseq/elapse.c)
  set(mumps_src ${mumps_src} ../ThirdParty/Mumps/MUMPS/${src})
endforeach ()
  
add_library(ipoptfort ${mumps_src} src/Algorithm/LinearSolvers/IpMa28Partition.F)

find_package(Threads)
target_link_libraries(ipoptfort ${CMAKE_THREAD_LIBS_INIT})

if (CMAKE_GNUtoMS)
  target_link_libraries(ipoptfort
    ${CMAKE_CURRENT_BINARY_DIR}/../ThirdParty/lapack/lib/libblas.dll.a)
else ()
  target_link_libraries(ipoptfort blas)
endif ()

if (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  target_link_libraries(ipoptfort gfortran)
endif ()
