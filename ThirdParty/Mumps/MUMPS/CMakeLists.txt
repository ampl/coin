cmake_minimum_required(VERSION 2.8)
project(MUMPS Fortran)

# Work-around for CMake issue 0009220.
if (DEFINED CMAKE_Fortran_COMPILER AND CMAKE_Fortran_COMPILER MATCHES "^$")
  set(CMAKE_Fortran_COMPILER CMAKE_Fortran_COMPILER-NOTFOUND)
endif ()

include(FortranCInterface)

set(F77_FUNC "${FortranCInterface_GLOBAL_MACRO}" CACHE INTERNAL "")
set(F77_FUNC_ "${FortranCInterface_GLOBAL__MACRO}" CACHE INTERNAL "")

if (FortranCInterface_GLOBAL_CASE STREQUAL "UPPER")
  add_definitions(-DUPPER)
elseif (FortranCInterface_GLOBAL_SUFFIX STREQUAL "_")
  add_definitions(-DAdd_)
elseif (FortranCInterface_GLOBAL_SUFFIX STREQUAL "__")
  add_definitions(-DAdd__)
endif ()

include_directories(include libseq)

add_library(mumps
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

find_package(Threads)
target_link_libraries(mumps ${CMAKE_THREAD_LIBS_INIT})

if (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  target_link_libraries(mumps gfortran)
endif ()
