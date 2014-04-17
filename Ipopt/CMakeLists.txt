# A CMake build script for the Fortran part of Ipopt.
# It is separate from the main build script to be able to use
# cmake_add_fortran_subdirectory and build Fortran part with a
# MinGW gfortran on Windows while the rest is built with Visual C++.

cmake_minimum_required(VERSION 2.8)
project(IPOPTFORT C)

enable_language(Fortran)
include(FortranCInterface)

# Organize output files.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${IPOPTFORT_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${IPOPTFORT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${IPOPTFORT_BINARY_DIR}/lib)

# Set F77_FUNC & F77_FUNC_ in cache as they are used in the parent project.
set(F77_FUNC "${FortranCInterface_GLOBAL_MACRO}" CACHE INTERNAL "")
set(F77_FUNC_ "${FortranCInterface_GLOBAL__MACRO}" CACHE INTERNAL "")

# Convert the information about Fortran name mangling into the format
# required by MUMPS.
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
  ../BuildTools/headers src/Common
  ../ThirdParty/Mumps/MUMPS/include ../ThirdParty/Mumps/MUMPS/libseq)

add_definitions(-DIPOPT_BUILD -DHAVE_CONFIG_H)

# Prepend directory to each filename passed as varargs and assign the result to var.
function (prepend_dir var dir)
  set(paths )
  foreach (filename ${ARGN})
    set(paths ${paths} ${dir}/${filename})
  endforeach ()
  set(${var} ${paths} PARENT_SCOPE)
endfunction ()

prepend_dir(LinearSolverLoader_src src/contrib/LinearSolverLoader
  LibraryHandler.c LibraryHandler.h
  HSLLoader.c HSLLoader.h
  PardisoLoader.c PardisoLoader.h)

prepend_dir(mumps_src ../ThirdParty/Mumps/MUMPS
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

add_library(ipoptfort ${mumps_src} ${LinearSolverLoader_src}
  src/Algorithm/LinearSolvers/IpMa28Partition.F)
find_library(DL_LIBRARY dl)
if (DL_LIBRARY)
  target_link_libraries(ipoptfort ${DL_LIBRARY})
endif ()

find_package(Threads)
target_link_libraries(ipoptfort ${CMAKE_THREAD_LIBS_INIT})

if (CMAKE_GNUtoMS)
  target_link_libraries(ipoptfort
    ${CMAKE_CURRENT_BINARY_DIR}/../ThirdParty/lapack/lib/libblas.dll.a)
else ()
  target_link_libraries(ipoptfort blas)
  add_ipopt_test(hs071_f hs071_f.f)
endif ()

if (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  # Get the list of library search directories from gfortran.
  if (NOT GFORTRAN_LIB_DIRS)
    exec_program(gfortran ARGS -print-search-dirs OUTPUT_VARIABLE output)
    if (output MATCHES "libraries\\: \\=(.*)")
      string(REGEX MATCHALL "[^:]+" dirs "${CMAKE_MATCH_1}")
      set(GFORTRAN_LIB_DIRS ${dirs} CACHE INTERNAL "")
    endif ()
  endif ()

  find_library(GFORTRAN_LIBRARY NAMES libgfortran.a PATHS ${GFORTRAN_LIB_DIRS})
  if (GFORTRAN_LIBRARY)
    target_link_libraries(ipoptfort ${GFORTRAN_LIBRARY})
  endif ()
endif ()
