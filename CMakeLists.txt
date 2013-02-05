# CMake build script for CBC.
# Author: Victor Zverovich

cmake_minimum_required(VERSION 2.6)

project(CBC)

# Output binaries to bin/ to avoid name clash between "Cbc" (directory)
# and "cbc" (executable) on case-insensitive filesystems (on Mac OS X).
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "bin/")

include(CheckCXXCompilerFlag)
if (COIN_DISABLE_WARNINGS)
  # Disable warnings.
  if (CMAKE_COMPILER_IS_GNUCXX)
    # Removing -Wextra seems to be the only way to disable warnings
    # "base class 'C' should be explicitly initialized in the copy constructor"
    # in older GCC.
    string(REPLACE "-Wextra" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    foreach(flag -Wno-format-security -Wno-overloaded-virtual
                 -Wno-unused-variable -Wno-unused-but-set-variable
                 -Wno-int-to-pointer-cast -Wno-uninitialized)
      set(var HAVE${flag})
      string(REPLACE "-" "_" var ${var})
      check_cxx_compiler_flag(${flag} ${var})
      if (${var})
        add_definitions(${flag})
      endif ()
    endforeach ()
  elseif (MSVC)
    add_definitions(-D_SCL_SECURE_NO_WARNINGS)
    foreach (flag_var
             CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
             CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
      string(REPLACE "/W3" "/w" ${flag_var} "${${flag_var}}")
    endforeach()
  endif ()
endif ()

# Parses an automake input processing conditional parts such as:
#   if VAR
#   ...
#   else
#   ...
#   endif
function(parse_am filename index keep)
  set(lines )
  while (1)
    list(LENGTH in_lines length)
    if (NOT (index LESS length))
      break()
    endif ()
    list(GET in_lines ${index} line)
    if (line MATCHES "^[ ]*endif")
      break ()
    elseif (line MATCHES "^[ ]*else")
      break ()
    elseif (line MATCHES "^[ ]*if[ \t]*([^ \t]*)")
      set(keep_then FALSE)
      set(keep_else TRUE)
      if (${CMAKE_MATCH_1})
        set(keep_then TRUE)
        set(keep_else FALSE)
      endif ()
      #message("if (${CMAKE_MATCH_1}): ${keep_then}")
      math(EXPR index "${index} + 1")
      parse_am(${filename} ${index} ${keep_then})
      set(index ${out_index})
      set(lines ${lines} ${out_lines})
      list(GET in_lines ${index} line)
      if (line MATCHES "^[ ]*else")
        math(EXPR index "${index} + 1")
        parse_am(${filename} ${index} ${keep_else})
        set(index ${out_index})
        set(lines ${lines} ${out_lines})
        list(GET in_lines ${index} line)
      endif ()
      if (NOT (line MATCHES "^[ ]*endif"))
        message(WARNING "${filename}:${index}: expected 'else' or 'endif'")
      endif ()
    elseif (keep)
      set(lines ${lines} "${line}")
    endif ()
    math(EXPR index "${index} + 1")
  endwhile ()
  set(out_index ${index} PARENT_SCOPE)
  set(out_lines ${lines} PARENT_SCOPE)
endfunction()

# Reads an automake file (Makefile.am).
function(read_am_file filename var)
  message("Parsing ${filename}")
  file(READ ${filename} text)
  string(REPLACE "\\\n" "" text "${text}")
  string(REGEX MATCHALL "[^\n]+" in_lines "${text}")
  parse_am(${filename} 0 TRUE)
  set(${var} ${out_lines} PARENT_SCOPE)
endfunction()

# Extracts the source list from Makefile.am.
function(get_sources path sources_var var)
  set(sources )
  read_am_file(${path}/Makefile.am LINES)
  foreach (line ${LINES})
    if (line MATCHES "^${sources_var} [+]?= (.*)")
      string(REGEX MATCHALL "[^ \t]+" src "${CMAKE_MATCH_1}")
      foreach (s ${src})
        if (EXISTS ${CBC_SOURCE_DIR}/${path}/${s})
          set(sources ${sources} ${path}/${s})
        endif ()
      endforeach ()
    endif ()
  endforeach ()
  set(${var} ${sources} PARENT_SCOPE)
endfunction ()

# Adds a COIN-OR library located at the specified path.
# Additional arguments specify library dependencies.
# By default the last path component is used as a library name.
# This can be overriden with NAME <libname> in arguments.
function(add_coin_library path)
  get_filename_component(name ${path} NAME)
  set(src_path ${path})
  set(test_path ${path}/test)
  if (EXISTS ${CBC_SOURCE_DIR}/${path}/src)
    set(src_path ${path}/src)
  endif ()
  # Process optional arguments.
  set(dependencies )
  set(expect_name False)
  set(expect_test_path False)
  set(notests False)
  foreach (arg ${ARGN})
    if (expect_name)
      set(name ${arg})
      set(expect_name False)
    elseif (expect_test_path)
      set(test_path ${arg})
      set(expect_test_path False)
    elseif (arg STREQUAL "NAME")
      set(expect_name True)
    elseif (arg STREQUAL "TEST_PATH")
      set(expect_test_path True)
    elseif (arg STREQUAL "NOTESTS")
      set(notests True)
    else ()
      set(dependencies ${dependencies} ${arg})
    endif ()
  endforeach ()
  # Add the library.
  get_sources(${src_path} lib${name}_la_SOURCES sources)
  add_library(${name} ${sources})
  include_directories(${src_path})
  string(TOUPPER ${name} uppercase_name)
  set_target_properties(${name}
    PROPERTIES COMPILE_DEFINITIONS ${uppercase_name}_BUILD)
  set(COIN_HAS_${uppercase_name} 1 PARENT_SCOPE)
  target_link_libraries(${name} ${dependencies})
  # Add tests.
  if (NOT notests AND EXISTS ${CBC_SOURCE_DIR}/${test_path})
    read_am_file(${test_path}/Makefile.am LINES)
    foreach (line ${LINES})
      if (line MATCHES "^(.*)_SOURCES [+]?= (.*)")
        set(test ${name}-${CMAKE_MATCH_1})
        set(sources )
        string(REGEX MATCHALL "[^ \t]+" src "${CMAKE_MATCH_2}")
        foreach (s ${src})
          set(sources ${sources} ${test_path}/${s})
        endforeach ()
        add_executable(${test} ${sources})
        set_target_properties(${test}
          PROPERTIES COMPILE_DEFINITIONS ${uppercase_name}_BUILD)
        target_link_libraries(${test} ${name})
        if (test MATCHES "osiUnitTest$")
          target_link_libraries(${test} OsiCommonTests)
        endif ()
        if (test MATCHES "^Cgl-")
          add_test(${test} bin/${test} Data/Sample Data/Netlib)
        else ()
          add_test(${test} bin/${test} -mpsDir=Data/Sample -netlibDir=Data/Netlib)
        endif ()
      endif ()
    endforeach ()
  endif ()
endfunction()

function(add_cgl_library path)
  add_coin_library(${path} Cgl ${ARGN})
  get_filename_component(name ${path} NAME)
  set_target_properties(${name} PROPERTIES COMPILE_DEFINITIONS CGL_BUILD)
endfunction()

# Processes config.h using autoconf rules.
function(configure package)
  # Get version numbers and parts of config.h from configure.ac.
  file(READ ${package}/configure.ac LINES)
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
      set(PACKAGE_NAME ${CMAKE_MATCH_1})
      set(VERSION ${CMAKE_MATCH_2})
      set(PACKAGE_BUGREPORT ${CMAKE_MATCH_3})
      # Remove surrounding '[' and ']'.
      foreach (var VERSION PACKAGE_NAME PACKAGE_BUGREPORT)
        if (${${var}} MATCHES "^\\[(.*)\\]$")
          set(${var} ${CMAKE_MATCH_1})
        endif ()
      endforeach ()
      message(STATUS "Got VERSION=${VERSION} from configure.ac")
    endif ()
  endforeach ()
  string(TOLOWER ${PACKAGE_NAME} PACKAGE)
  set(PACKAGE_TARNAME ${PACKAGE})
  set(PACKAGE_URL "")
  set(PACKAGE_VERSION ${VERSION})
  string(TOUPPER ${PACKAGE} UPPERCASE_PACKAGE)
  set(${UPPERCASE_PACKAGE}_VERSION ${VERSION})
  set(PACKAGE_STRING "${PACKAGE_NAME} ${VERSION}")
  if (VERSION MATCHES "(.*)\\.(.*)\\.(.*)")
    set(${UPPERCASE_PACKAGE}_VERSION_MAJOR ${CMAKE_MATCH_1})
    set(${UPPERCASE_PACKAGE}_VERSION_MINOR ${CMAKE_MATCH_2})
    set(${UPPERCASE_PACKAGE}_VERSION_RELEASE ${CMAKE_MATCH_3})
  endif ()
  set(COIN_${UPPERCASE_PACKAGE}_CHECKLEVEL 0)
  set(COIN_${UPPERCASE_PACKAGE}_VERBOSITY 0)
  # Process config files.
  foreach (config_file ${ARGN})
    file(STRINGS ${package}/${config_file}.in CONFIG)
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
          if (${var} MATCHES "^-?[0-9]+$")
            set(value ${${var}})
          elseif (${var})
            set(value \"${${var}}\")
          endif ()
          set(line "#${space}define ${var} ${value}")
         endif ()
      endif ()
      string(REPLACE "<semi>" ";" line "${line}")
      set(CONFIG_OUT "${CONFIG_OUT}${line}\n")
    endforeach ()
    get_filename_component(name ${config_file} NAME)
    file(WRITE ${package}/${config_file}
"/* ${config_file}.  Generated by configure.  */
${CONFIG_OUT}")
  endforeach ()
endfunction()

# Check C headers.
include(CheckIncludeFiles)
foreach (header bzlib dlfcn endian inttypes memory stdlib stdint
                strings string sys/stat sys/types unistd windows zlib)
  string(TOUPPER ${header} var)
  string(REPLACE "/" "_" var ${var})
  check_include_files(${header}.h HAVE_${var}_H)
endforeach ()
if (HAVE_STDLIB_H)
  set(STDC_HEADERS 1)
endif ()

# Check C++ headers.
include(CheckIncludeFileCXX)
check_include_file_cxx(cfloat HAVE_CFLOAT)
check_include_file_cxx(cmath HAVE_CMATH)
if (NOT HAVE_CMATH)
  check_include_files(math.h HAVE_MATH_H)
endif ()

include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
  #include <math.h>
  int main() { finite(4.2); }" COIN_HAS_FINITE)
include(CheckSymbolExists)
if (COIN_HAS_FINITE)
  set(COIN_C_FINITE_C finite)
else ()
  check_symbol_exists(isfinite math.h COIN_HAS_ISFINITE)
  if (COIN_HAS_ISFINITE)
    set(COIN_C_FINITE_C isfinite)
  endif ()
endif ()
check_cxx_source_compiles("
  #include <cmath>
  int main() { std::isnan(4.2); }" COIN_HAS_STD_ISNAN)
if (COIN_HAS_STD_ISNAN)
  set(COIN_C_ISNAN_C std::isnan)
else ()
  check_cxx_source_compiles("
    #include <math.h>
    int main() { isnan(4.2); }" COIN_HAS_ISNAN)
  if (COIN_HAS_FINITE)
    set(COIN_C_ISNAN_C isnan)
  endif ()
endif ()

set(COINUTILS_LIBS )

find_package(ZLIB)
if (ZLIB_FOUND)
  set(COIN_HAS_ZLIB 1)
  set(COINUTILS_LIBS ${COINUTILS_LIBS} ${ZLIB_LIBRARIES})
endif ()

find_package(BZip2)
if (BZIP2_FOUND)
  set(COIN_HAS_BZLIB 1)
  set(COINUTILS_LIBS ${COINUTILS_LIBS} ${BZIP2_LIBRARIES})
endif ()

include(CheckTypeSize)
check_cxx_source_compiles("
  #include <stdint.h>
  int main() {
    int64_t x;
    uint64_t y;
    intptr_t z;
  }" COIN_HAS_STDINT)
if (COIN_HAS_STDINT)
  set(COIN_INT64_T_C int64_t)
  # Define COIN_UINT64_T as int64_t for compatibility,
  # may be a bug in CoinUtils.
  set(COIN_UINT64_T_C int64_t)
  set(COIN_INTPTR_T_C intptr_t)
else ()
  check_type_size("long long" LONG_LONG_SIZE)
  if (LONG_LONG_SIZE EQUAL 8)
    set(COIN_INT64_T_C "long long")
    set(COIN_UINT64_T_C "unsigned long long")
  endif ()
endif ()

set(CBC_DEFAULT_SOLVER clp)
string(SUBSTRING ${CBC_DEFAULT_SOLVER} 0 1 solver_head)
string(SUBSTRING ${CBC_DEFAULT_SOLVER} 1 -1 solver_tail)
string(TOUPPER ${solver_head} solver_head)
set(OSICBC_DFLT_SOLVER_C Osi${solver_head}${solver_tail}SolverInterface)
set(OSICBC_DFLT_SOLVER_HPP ${OSICBC_DFLT_SOLVER_C}.hpp)
string(TOUPPER ${CBC_DEFAULT_SOLVER} solver_uppercase)
set(OSICBC_DFLT_SOLVER_${solver_uppercase} 1)
set(COIN_HAS_COINDEPEND 1)
set(COIN_HAS_OSITESTS 1)
set(COIN_HAS_SAMPLE 1)
set(COINUTILS_MEMPOOL_MAXPOOLED -1)
set(GLPK_HAS_INTOPT 1)

enable_testing()

include_directories(BuildTools/headers)

if (COIN_HAS_ASL)
  add_definitions(-DCOIN_HAS_ASL)
endif ()

add_coin_library(CoinUtils ${COINUTILS_LIBS})
add_coin_library(Clp CoinUtils)
add_coin_library(Clp/src/OsiClp Osi Clp)
add_coin_library(Osi/src/Osi CoinUtils)
add_coin_library(Osi/src/OsiCommonTest
  NAME OsiCommonTests TEST_PATH Osi/test Osi OsiClp)

add_coin_library(Cgl Osi)
add_cgl_library(Cgl/src/CglClique)
add_cgl_library(Cgl/src/CglDuplicateRow)
add_cgl_library(Cgl/src/CglFlowCover)
add_cgl_library(Cgl/src/CglGomory)
add_cgl_library(Cgl/src/CglKnapsackCover)
add_cgl_library(Cgl/src/CglLandP)
add_cgl_library(Cgl/src/CglMixedIntegerRounding)
add_cgl_library(Cgl/src/CglMixedIntegerRounding2)
add_cgl_library(Cgl/src/CglOddHole)
add_cgl_library(Cgl/src/CglPreProcess CglClique CglDuplicateRow)
add_cgl_library(Cgl/src/CglProbing)
add_cgl_library(Cgl/src/CglRedSplit)
add_cgl_library(Cgl/src/CglResidualCapacity)
add_cgl_library(Cgl/src/CglSimpleRounding)
add_cgl_library(Cgl/src/CglTwomir)
target_link_libraries(Cgl-unitTest
  OsiClp CglClique CglFlowCover CglGomory CglKnapsackCover CglLandP
  CglMixedIntegerRounding CglMixedIntegerRounding2 CglOddHole CglProbing
  CglRedSplit CglResidualCapacity CglSimpleRounding CglTwomir)

add_coin_library(Cbc
  OsiClp CglFlowCover CglGomory CglKnapsackCover
  CglMixedIntegerRounding2 CglPreProcess CglProbing CglTwomir)
add_coin_library(Cbc/src/OsiCbc Osi Cbc)
target_link_libraries(Cbc-osiUnitTest OsiCbc)
target_link_libraries(Cbc-gamsTest CbcSolver)

add_coin_library(Cbc NAME CbcSolver NOTESTS
  Cbc Cgl CglLandP CglRedSplit CglResidualCapacity CglTwomir ${ASL_LIBRARIES})
set_target_properties(CbcSolver PROPERTIES COMPILE_DEFINITIONS CBC_BUILD)

add_executable(cbc-executable Cbc/src/CoinSolve.cpp)
set_target_properties(cbc-executable
  PROPERTIES OUTPUT_NAME cbc IMPORT_PREFIX imp)
target_link_libraries(cbc-executable CbcSolver)
add_test(Cbc-unitTest bin/cbc -dirSample Data/Sample -unitTest)

if (NOT MSVC)
  add_definitions(-DHAVE_CONFIG_H)
endif ()

configure(CoinUtils src/config.h src/config_coinutils.h)
configure(Osi src/Osi/config.h src/Osi/config_osi.h)
configure(Clp src/config.h src/config_clp.h)
configure(Cgl src/config.h src/config_cgl.h)
configure(Cbc src/config.h src/config_cbc.h)
