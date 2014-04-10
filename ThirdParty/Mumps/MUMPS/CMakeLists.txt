cmake_minimum_required(VERSION 2.6)
project(MUMPS Fortran)

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
  src/mumps_sol_es.F)
