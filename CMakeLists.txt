# CMake build script for CBC.
# Author: Victor Zverovich

cmake_minimum_required(VERSION 2.6)

project(CBC)

function(add_coin_library_internal name path)
  set(sources )
  foreach (arg ${ARGN})
    set(sources ${sources} ${path}/${arg})
  endforeach ()
  include_directories(${path})
  add_library(${name} ${sources})
  string(TOUPPER ${name} uppercase_name)
  add_definitions(-D${uppercase_name}_BUILD)
endfunction()

function(read_automake filename var)
  file(READ ${filename} text)
  string(REPLACE "\\\n" "" text "${text}")
  string(REGEX MATCHALL "[^\n]+" lines "${text}")
  set(${var} ${lines} PARENT_SCOPE)
endfunction()

# Adds a COIN-OR library located at the specified path.
# Additional arguments specify library dependencies.
# By default the last path component is used as a library name.
# This can be overriden with NAME <libname> in arguments.
function(add_coin_library path)
  get_filename_component(name ${path} NAME)
  set(src_path ${path})
  if (EXISTS ${path}/src)
    set(src_path ${path}/src)
  endif ()
  # Process optional arguments.
  set(dependencies )
  set(expect_name False)
  foreach (arg ${ARGN})
    if (expect_name)
      set(name ${arg})
      set(expect_name False)
    elseif (arg STREQUAL "NAME")
      set(expect_name True)
    else ()
      set(dependencies ${dependencies} ${arg})
    endif ()
  endforeach ()
  # Extract the source list from Makefile.am.
  set(sources )
  read_automake(${src_path}/Makefile.am LINES)
  foreach (line ${LINES})
    if (line MATCHES "^lib${name}_la_SOURCES [+]?= (.*)")
      string(REGEX MATCHALL "[^ \t]+" src "${CMAKE_MATCH_1}")
      set(sources ${sources} ${src})
    endif ()
  endforeach ()
  add_coin_library_internal(${name} ${src_path} ${sources})
  target_link_libraries(${name} ${dependencies})
endfunction()

macro(add_coin_package name)
  add_coin_library_internal(${name} ${name}/src ${ARGN})
  string(TOUPPER ${name} uppercase_name)
  set(COIN_HAS_${uppercase_name} 1)
endmacro()

# Process config.h using autoconf rules.
function(configure file)
  file(STRINGS ${file}.in CONFIG)
  list(LENGTH CONFIG length)
  math(EXPR length "${length} - 1")
  set(CONFIG_OUT )
  foreach (i RANGE ${length})
    list(GET CONFIG ${i} line)
    if (line MATCHES "^#( *)undef (.*)")
      set(space "${CMAKE_MATCH_1}")
      set(var ${CMAKE_MATCH_2})
      if (${var}_C)
        # Variables with suffix "_C" are defined unquoted.
        set(line "#${space}define ${var} ${${var}_C}")
      elseif (NOT DEFINED ${var} OR (var MATCHES "HAVE_.*_H" AND NOT ${var}))
        set(line "/* #${space}undef ${var} */")
      else ()
        if (NOT (var MATCHES ^HAVE OR ${var} EQUAL 1))
          set(value \"${${var}}\")
        elseif (${var})
          set(value 1)
        else ()
          set(value 0)
        endif ()
        set(line "#${space}define ${var} ${value}")
       endif ()
    endif ()
    string(REPLACE "<semi>" ";" line "${line}")
    set(CONFIG_OUT "${CONFIG_OUT}${line}\n")
  endforeach ()
  file(WRITE ${file}
"/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

${CONFIG_OUT}")
endfunction()

include(CheckIncludeFiles)
check_include_files(math.h HAVE_MATH_H)

include(CheckSymbolExists)
include(CheckTypeSize)
check_symbol_exists(isfinite math.h COIN_HAS_ISFINITE)
if (COIN_HAS_ISFINITE)
  set(COIN_C_FINITE_C isfinite)
endif ()
check_symbol_exists(int64_t inttypes.h COIN_HAS_INT64_T)
if (COIN_HAS_INT64_T)
  set(COIN_INT64_T_C int64_t)
  set(COIN_UINT64_T_C uint64_t)
else ()
  check_type_size("long long" LONG_LONG_SIZE)
  if (LONG_LONG_SIZE EQUAL 8)
    set(COIN_INT64_T_C "long long")
    set(COIN_UINT64_T_C "unsigned long long")
  endif ()
endif ()

# Get version numbers and parts of config.h from configure.ac.
file(READ configure.ac LINES)
# Replace semicolons with "<semi>" to avoid CMake messing with them.
string(REPLACE ";" "<semi>" LINES "${LINES}")
# Split into lines keeping newlines to avoid foreach skipping empty ones.
string(REGEX MATCHALL "[^\n]*\n" LINES "${LINES}")
set(ah_command FALSE)
foreach (line "${EXTRA_CONFIG}" ${LINES})
  string(REPLACE ";" "" line "${line}")
  if (ah_command)
    # Do nothing.
  elseif (line MATCHES "AC_INIT\\(([^,]*), *([^,]*), *([^)]*)\\)")
    set(PACKAGE ${CMAKE_MATCH_1})
    set(VERSION ${CMAKE_MATCH_2})
    set(PACKAGE_BUGREPORT ${CMAKE_MATCH_3})
    foreach (var VERSION PACKAGE PACKAGE_BUGREPORT)
      if (${${var}} MATCHES "^\\[(.*)\\]$")
        set(${var} ${CMAKE_MATCH_1})
      endif ()
    endforeach ()
    message(STATUS "Got VERSION=${VERSION} from configure.ac")
  endif ()
endforeach ()
set(PACKAGE_NAME ${PACKAGE})
string(TOLOWER ${PACKAGE} PACKAGE_TARNAME)
set(PACKAGE_URL "")
set(PACKAGE_VERSION ${VERSION})
string(TOUPPER ${PACKAGE} UPPERCASE_PACKAGE)
set(${UPPERCASE_PACKAGE}_VERSION ${VERSION})
set(PACKAGE_STRING "${PACKAGE} ${VERSION}")

add_coin_library(CoinUtils)

add_coin_package(Clp
  ClpConfig.h
  ClpCholeskyBase.cpp ClpCholeskyBase.hpp
  ClpCholeskyDense.cpp ClpCholeskyDense.hpp
  ClpConstraint.cpp ClpConstraint.hpp
  ClpConstraintLinear.cpp ClpConstraintLinear.hpp
  ClpConstraintQuadratic.cpp ClpConstraintQuadratic.hpp
  Clp_C_Interface.cpp Clp_C_Interface.h
  ClpDualRowDantzig.cpp ClpDualRowDantzig.hpp
  ClpDualRowPivot.cpp ClpDualRowPivot.hpp
  ClpDualRowSteepest.cpp ClpDualRowSteepest.hpp
  ClpDummyMatrix.cpp ClpDummyMatrix.hpp
  ClpDynamicExampleMatrix.cpp ClpDynamicExampleMatrix.hpp
  ClpDynamicMatrix.cpp ClpDynamicMatrix.hpp
  ClpEventHandler.cpp ClpEventHandler.hpp
  ClpFactorization.cpp ClpFactorization.hpp
  ClpGubDynamicMatrix.cpp ClpGubDynamicMatrix.hpp
  ClpGubMatrix.cpp ClpGubMatrix.hpp
  ClpHelperFunctions.cpp ClpHelperFunctions.hpp
  ClpInterior.cpp ClpInterior.hpp
  ClpLinearObjective.cpp ClpLinearObjective.hpp
  ClpMatrixBase.cpp ClpMatrixBase.hpp
  ClpMessage.cpp ClpMessage.hpp
  ClpModel.cpp ClpModel.hpp
  ClpNetworkBasis.cpp ClpNetworkBasis.hpp
  ClpNetworkMatrix.cpp ClpNetworkMatrix.hpp
  ClpNonLinearCost.cpp ClpNonLinearCost.hpp
  ClpNode.cpp ClpNode.hpp
  ClpObjective.cpp ClpObjective.hpp
  ClpPackedMatrix.cpp ClpPackedMatrix.hpp
  ClpParameters.hpp
  ClpPlusMinusOneMatrix.cpp ClpPlusMinusOneMatrix.hpp
  ClpPredictorCorrector.cpp ClpPredictorCorrector.hpp
  ClpPdco.cpp ClpPdco.hpp
  ClpPdcoBase.cpp ClpPdcoBase.hpp
  ClpLsqr.cpp ClpLsqr.hpp
  ClpPresolve.cpp ClpPresolve.hpp
  ClpPrimalColumnDantzig.cpp ClpPrimalColumnDantzig.hpp
  ClpPrimalColumnPivot.cpp ClpPrimalColumnPivot.hpp
  ClpPrimalColumnSteepest.cpp ClpPrimalColumnSteepest.hpp
  ClpQuadraticObjective.cpp ClpQuadraticObjective.hpp
  ClpSimplex.cpp ClpSimplex.hpp
  ClpSimplexDual.cpp ClpSimplexDual.hpp
  ClpSimplexNonlinear.cpp ClpSimplexNonlinear.hpp
  ClpSimplexOther.cpp ClpSimplexOther.hpp
  ClpSimplexPrimal.cpp ClpSimplexPrimal.hpp
  ClpSolve.cpp ClpSolve.hpp
  Idiot.cpp Idiot.hpp
  IdiSolve.cpp)
target_link_libraries(Clp CoinUtils)

add_coin_library(Clp/src/OsiClp Osi Clp)

add_coin_library(Osi/src/Osi)

add_coin_library(Cgl)
add_coin_library(Cgl/src/CglClique)
add_coin_library(Cgl/src/CglDuplicateRow)
add_coin_library(Cgl/src/CglFlowCover)
add_coin_library(Cgl/src/CglGomory)
add_coin_library(Cgl/src/CglKnapsackCover)
add_coin_library(Cgl/src/CglLandP)
add_coin_library(Cgl/src/CglMixedIntegerRounding2)
add_coin_library(Cgl/src/CglOddHole)
add_coin_library(Cgl/src/CglPreProcess CglClique CglDuplicateRow)
add_coin_library(Cgl/src/CglProbing)
add_coin_library(Cgl/src/CglRedSplit)
add_coin_library(Cgl/src/CglResidualCapacity)
add_coin_library(Cgl/src/CglTwomir)

add_coin_library(Cbc
  OsiClp CglFlowCover CglGomory CglKnapsackCover
  CglMixedIntegerRounding2 CglPreProcess CglProbing)

add_coin_library(Cbc NAME CbcSolver
  Cbc Cgl CglLandP CglRedSplit CglResidualCapacity CglTwomir)

add_executable(cbc Cbc/src/CoinSolve.cpp)
target_link_libraries(cbc CbcSolver)

add_definitions(-DHAVE_CONFIG_H)
configure(${CBC_SOURCE_DIR}/CoinUtils/src/config.h)
configure(${CBC_SOURCE_DIR}/Osi/src/Osi/config_osi.h)
configure(${CBC_SOURCE_DIR}/Cbc/src/config.h)
