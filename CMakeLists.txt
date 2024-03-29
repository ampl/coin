# CMake build script for COIN-OR projects.
# Author: Victor Zverovich
# This file is published under the Eclipse Public License (EPL).

cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION ${CMAKE_MINIMUM_REQUIRED_VERSION})

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
message(STATUS "CMake Version: ${CMAKE_VERSION}")
if (BUILD_SHARED_LIBS)
  set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
endif ()

set(ASL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/asl)
if (EXISTS ${ASL_DIR}/CMakeLists.txt)
  set(AMPL_CMAKE_MODULE_DIR ${ASL_DIR}/support/cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${AMPL_CMAKE_MODULE_DIR})
  # Define the ARCH variable (32 or 64)
  include(initArchitecture)
  # Set generic architecture compiler flags
  include(setArchitecture)
  getArchitectureFlags(${ARCH} CCFLAGS LLFLAGS)
  #add_compile_options(${CCFLAGS})

  add_subdirectory(${ASL_DIR})
  set(COIN_HAS_ASL 1)
  set(ASL_LIBRARIES asl)
  add_definitions(-DCOIN_HAS_ASL)
  include_directories($<TARGET_PROPERTY:asl,INTERFACE_INCLUDE_DIRECTORIES>)
endif ()

project(COIN)

# Replaces libraries in CMAKE_${lang}_IMPLICIT_LINK_LIBRARIES whose
# names appear in ARGN with their static counterparts.
# Usage:
#   force_static_libs(CXX stdc++)
function (force_static_libs lang)
  set(result )
  foreach (lib ${CMAKE_${lang}_IMPLICIT_LINK_LIBRARIES})
    list(FIND ARGN ${lib} found)
    if (NOT found EQUAL -1)
      find_library(${lib}_LIB
        ${CMAKE_STATIC_LIBRARY_PREFIX}${lib}${CMAKE_STATIC_LIBRARY_SUFFIX}
        ${CMAKE_${lang}_IMPLICIT_LINK_DIRECTORIES})
      if (${lib}_LIB)
        set(lib ${${lib}_LIB})
      endif ()
    endif ()
    set(result ${result} ${lib})
  endforeach ()
  set(CMAKE_${lang}_IMPLICIT_LINK_LIBRARIES ${result} PARENT_SCOPE)
endfunction ()

# Output binaries to bin to avoid a name clash between the Cbc directory
# and the cbc executable on case-insensitive filesystems such as HFS Plus
# (on OS X).
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# CoinSignal.hpp doesn't compile on Mac with -std=c++11 flag
string(REPLACE "-std=c++11" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
string(REPLACE "-stdlib=libc++" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# Adds flags passed as varargs to CMAKE_CXX_FLAGS.
function(add_cxx_flags)
  foreach (flag ${ARGN})
    set(flags "${flags} ${flag}")
  endforeach ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flags}" PARENT_SCOPE)
endfunction()


include(CheckCXXCompilerFlag)
# Always disable certain warnings
if(MSVC)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_WARNINGS)
  # Link with the static runtime
  add_compile_options(/MT$<$<CONFIG:Debug>:d>)
endif()

if (COIN_DISABLE_WARNINGS)
  # Disable warnings.
  if (MSVC)
    add_definitions(/wd4065 /wd4101 /wd4102 /wd4244 /wd4267 /wd4297 /wd4355 /wd4800  /wd4996 /wd4309)
  else ()
    foreach (flag -Wall -Wextra -pedantic)
      string(REPLACE ${flag} "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
      string(REPLACE ${flag} "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    endforeach ()
    if (CMAKE_COMPILER_IS_GNUCXX)
      add_cxx_flags(
        -Wno-format-security -Wno-overloaded-virtual -Wno-unused-result)
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      add_cxx_flags(
        -Wno-constant-conversion -Wno-constant-logical-operand
        -Wno-format-security -Wno-language-extension-token -Wno-newline-eof
        -Wno-parentheses -Wno-parentheses-equality -Wno-self-assign-field
        -Wno-sizeof-pointer-memaccess -Wno-tautological-compare
        -Wno-unused-comparison)
    endif ()
  endif ()
endif ()

option(ENABLE_CBC_PARALLEL "Enables compilation of the SMP version of Cbc" ON)
find_package(Threads)
set(PTHREADS_W32_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/pthreads-w32)
if (NOT ENABLE_CBC_PARALLEL)
  # Do nothing.
elseif (CMAKE_USE_PTHREADS_INIT)
  set(CBC_THREAD 1)
  set(CBC_PTHREADS_LIB ${CMAKE_THREAD_LIBS_INIT})
elseif (WIN32 AND EXISTS ${PTHREADS_W32_DIR})
  set(sources )
  foreach (src attr.c barrier.c cancel.c cleanup.c condvar.c create.c dll.c
               autostatic.c errno.c exit.c fork.c global.c misc.c mutex.c
               nonportable.c private.c rwlock.c sched.c semaphore.c signal.c
               spin.c sync.c tsd.c)
    set(sources ${sources} ${PTHREADS_W32_DIR}/${src})
  endforeach ()
  include_directories(${PTHREADS_W32_DIR})
  add_definitions(-DPTW32_STATIC_LIB=1)
  add_library(pthreads STATIC ${sources})
  set(CBC_THREAD 1)
  set(CBC_PTHREADS_LIB pthreads)
  set(COIN_INSTALL_TARGETS ${COIN_INSTALL_TARGETS} pthreads)
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
function(coin_get_sources path name var noinst_var)
  set(sources_var lib${name}_la_SOURCES)
  set(${noinst_var} FALSE PARENT_SCOPE)
  set(sources )
  read_am_file(${path}/Makefile.am LINES)
  foreach (line ${LINES})
    if (line MATCHES "^[ \t]*noinst_LTLIBRARIES = ([^ ]*)")
      if (${CMAKE_MATCH_1} STREQUAL "lib${name}.la")
        set(${noinst_var} TRUE PARENT_SCOPE)
      endif ()
    elseif (line MATCHES "^[ \t]*(nodist_)?${sources_var} [+]?= (.*)")
      string(REGEX MATCHALL "[^ \t]+" source_list "${CMAKE_MATCH_2}")
      foreach (s ${source_list})
        set(src ${s})
        # Translate variable reference of the form $(name) in the source name.
        if (src MATCHES "(.*)[$][(]([^)]+)[)](.*)")
          set(src "${CMAKE_MATCH_1}${${CMAKE_MATCH_2}}${CMAKE_MATCH_3}")
        endif ()
        if (NOT EXISTS ${COIN_SOURCE_DIR}/${path}/${src})
          # do nothing
        elseif (src MATCHES "\\.[fF]$")
          message("Skipped ${src}")
        else ()
          set(sources ${sources} ${path}/${src})
        endif ()
      endforeach ()
    elseif (line MATCHES "^includecoin_HEADERS \\+?= (.*)")
      string(REGEX MATCHALL "[^ \t]+" header_list "${CMAKE_MATCH_1}")
      set(headers )
      foreach (header ${header_list})
        set(headers ${headers} ${path}/${header})
      endforeach ()
      install(FILES ${headers} DESTINATION include/coin)
    elseif (line MATCHES "\\$\\(install_sh_DATA\\) +([a-zA-Z_]+\\.hp?p?) +[A-Za-z_$()]+/([A-Za-z_]+.hp?p?)")
      install(FILES ${path}/${CMAKE_MATCH_1} RENAME ${CMAKE_MATCH_2} DESTINATION include/coin)
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
  if (EXISTS ${COIN_SOURCE_DIR}/${path}/src)
    set(src_path ${path}/src)
  endif ()
  # Process optional arguments.
  set(objects )
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
      # Check if this is an object library and add its content to
      # objects if yes.
      set(type )
      if (TARGET ${arg})
        get_target_property(type ${arg} TYPE)
      endif ()
      if (type STREQUAL OBJECT_LIBRARY)
        set(objects ${objects} "$<TARGET_OBJECTS:${arg}>")
      else ()
        set(dependencies ${dependencies} ${arg})
      endif ()
    endif ()
  endforeach ()
  # Add the library.
  coin_get_sources(${src_path} ${name} sources noinst)
  set(libtype )
  if (${noinst})
    set(libtype OBJECT)
  endif ()
  add_library(${name} ${libtype} ${sources} ${objects})
  include_directories(${src_path})
  string(TOUPPER ${name} uppercase_name)
  set_target_properties(${name}
    PROPERTIES COMPILE_DEFINITIONS ${uppercase_name}_BUILD)
  set(COIN_HAS_${uppercase_name} 1 PARENT_SCOPE)
  if (NOT ${noinst})
    target_link_libraries(${name} ${dependencies})
    set(COIN_INSTALL_TARGETS ${COIN_INSTALL_TARGETS} ${name} CACHE INTERNAL "")
  endif ()
  # Add tests.
  if (NOT notests AND EXISTS ${COIN_SOURCE_DIR}/${test_path})
    read_am_file(${test_path}/Makefile.am LINES)
    foreach (line ${LINES})
      if (line MATCHES "^(.*)_SOURCES [+]?= (.*)")
        set(test ${name}-${CMAKE_MATCH_1})
        set(sources )
        string(REGEX MATCHALL "[^ \t]+" src "${CMAKE_MATCH_2}")
        foreach (s ${src})
          set(source ${test_path}/${s})
          if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${source})
            set(sources ${sources} ${source})
          endif ()
        endforeach ()
        if (NOT sources)
          message("Ignoring test without sources: ${test}")
        else ()
          add_executable(${test} ${sources})
          add_to_folder(test ${test})
          set_target_properties(${test}
            PROPERTIES COMPILE_DEFINITIONS ${uppercase_name}_BUILD)
          target_link_libraries(${test} ${name})
          if (test MATCHES "osiUnitTest$")
            target_link_libraries(${test} OsiCommonTests)
          endif ()
          if (test MATCHES "^Cgl-")
            add_test(NAME ${test}
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/bin/${test} Cbc/Data/Sample Cbc/Data/Netlib
              WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
          else ()
            add_test(NAME ${test}
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/bin/${test} -mpsDir=Cbc/Data/Sample -netlibDir=Cbc/Data/Netlib
              WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
          endif ()
        endif ()
      endif ()
    endforeach ()
  endif ()
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
        elseif (var MATCHES "F77_FUNC")
          # F77_FUNC shouldn't have space between macro name and args.
          set(line "#${space}define ${var}${${var}}")
        elseif (NOT DEFINED ${var} OR (var MATCHES "HAVE_.*" AND NOT ${var}))
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
foreach (header cassert cctype cfloat cmath cstdarg
                cstddef cstdio cstdlib cstring ctime)
  string(TOUPPER ${header} var)
  check_include_file_cxx(${header} HAVE_${var})
endforeach ()
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
  elseif (MSVC)
    set(COIN_C_FINITE_C _finite)
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
  if (COIN_HAS_ISNAN)
    set(COIN_C_ISNAN_C isnan)
  elseif (MSVC)
    set(COIN_C_ISNAN_C _isnan)
  endif ()
endif ()

set(COINUTILS_LIBS )

find_package(ZLIB)
if (ZLIB_FOUND)
  set(COIN_HAS_ZLIB 1)
  set(COINUTILS_LIBS ${COINUTILS_LIBS} ${ZLIB_LIBRARIES})
endif ()

if (NOT WIN32)
  find_package(BZip2)
  if (BZIP2_FOUND)
    set(COIN_HAS_BZLIB 1)
    set(COINUTILS_LIBS ${COINUTILS_LIBS} ${BZIP2_LIBRARIES})
  endif ()
endif ()

include(CheckTypeSize)
check_type_size(int* SIZEOF_INT_P)

include(CheckTypeSize)
check_cxx_source_compiles("
  #include <stdint.h>
  int main() {
    int64_t x;
    uint64_t y;
    intptr_t z;
  }" COIN_HAS_STDINT)
if (MSVC)
  set(COIN_INT64_T_C __int64)
  set(COIN_UINT64_T_C __int64)
  if (SIZEOF_INT_P EQUAL 8)
    set(COIN_INTPTR_T_C __int64)
  else ()
    set(COIN_INTPTR_T_C __int32)
  endif ()
elseif (COIN_HAS_STDINT)
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

set(CBC_DEFAULT_SOLVER_NAME Clp)
string(TOLOWER ${CBC_DEFAULT_SOLVER_NAME} CBC_DEFAULT_SOLVER)
set(OSICBC_DFLT_SOLVER_C Osi${CBC_DEFAULT_SOLVER_NAME}SolverInterface)
set(OSICBC_DFLT_SOLVER_HPP ${OSICBC_DFLT_SOLVER_C}.hpp)
string(TOUPPER ${CBC_DEFAULT_SOLVER} solver_uppercase)
set(OSICBC_DFLT_SOLVER_${solver_uppercase} 1)
set(COIN_HAS_COINDEPEND 1)
set(COIN_HAS_OSITESTS 1)
set(COIN_HAS_SAMPLE 1)
set(COINUTILS_MEMPOOL_MAXPOOLED -1)
set(GLPK_HAS_INTOPT 1)

if (COIN_HAS_ASL)
  # Workaround for f90cache which accept only modified name for GNU fortran
  # compilers in the form gfortran-x[.y].
  string(REGEX MATCHALL "[^:]+" search_paths "$ENV{PATH}")
  foreach (search_path ${search_paths})
    file(GLOB paths ${search_path}/gfortran-5*)
    if (paths)
      list(GET paths 0 CMAKE_Fortran_COMPILER)
      message("Detected Fortran compiler: ${CMAKE_Fortran_COMPILER}")
      break ()
    endif ()
  endforeach ()

  #if (MP_SETENV)
    # Create a script that sets up a build environment in Fortran subprojects.
    # It is necessary because setenv.cmd produces bogus errors
    # ERROR: The system was unable to find the specified registry key or value.
    # which cause custom build to fail.
    #set(MP_SETENV "C:\\Program Files (x86)\\Microsoft Visual Studio\\2019\\Community\\Common7\\Tools\\VsDevCmd.bat")
    set(MP_SETENV "")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/run-setenv.bat "call \"${MP_SETENV}\" %* 2> nul")
    set(path ${CMAKE_CURRENT_BINARY_DIR}/run-setenv.bat)
    set(FORTRAN_CMAKE_FLAGS -DCMAKE_GNUtoMS_VCVARS=${path} -DVCVARS32=${path})
  endif ()
#endif ()

enable_testing()

include_directories(Cbc/BuildTools/headers)

add_coin_library(Cbc/CoinUtils ${COINUTILS_LIBS})
add_coin_library(Cbc/Clp CoinUtils)
add_coin_library(Cbc/Clp/src/OsiClp Osi Clp)
add_coin_library(Cbc/Osi/src/Osi CoinUtils)
add_coin_library(Cbc/Osi/src/OsiCommonTest
  NAME OsiCommonTests TEST_PATH Osi/test Osi OsiClp)

foreach (name AllDifferent Clique DuplicateRow FlowCover GMI Gomory
              KnapsackCover LandP LiftAndProject MixedIntegerRounding
              MixedIntegerRounding2 OddHole PreProcess Probing RedSplit
              RedSplit2 ResidualCapacity SimpleRounding Twomir ZeroHalf)
  add_coin_library(Cbc/Cgl/src/Cgl${name})
  set(CGL_LIBS ${CGL_LIBS} Cgl${name})
  set_target_properties(Cgl${name} PROPERTIES COMPILE_DEFINITIONS CGL_BUILD)
  add_to_folder(cbc/cgl Cgl${name})
endforeach ()
add_coin_library(Cbc/Cgl ${CGL_LIBS} OsiClp Osi)
set(LIB_RT "")
if(UNIX AND NOT APPLE)
    set(LIB_RT "rt")
endif()
add_coin_library(Cbc/Cbc Cgl OsiClp ${CBC_PTHREADS_LIB} ${LIB_RT})
add_coin_library(Cbc/Cbc/src/OsiCbc Osi Cbc)
target_link_libraries(Cbc-osiUnitTest OsiCbc)
target_link_libraries(Cbc-gamsTest CbcSolver)
target_link_libraries(Cbc-CInterfaceTest CbcSolver)
add_coin_library(Cbc/Cbc NAME CbcSolver NOTESTS Cbc Cgl ${ASL_LIBRARIES})
set_target_properties(CbcSolver PROPERTIES COMPILE_DEFINITIONS CBC_BUILD)

add_executable(cbc-executable Cbc/Cbc/src/CoinSolve.cpp)
set_target_properties(cbc-executable
  PROPERTIES OUTPUT_NAME cbc IMPORT_PREFIX imp)
target_link_libraries(cbc-executable CbcSolver)

add_to_folder(cbc/cgl Cgl)
add_to_folder(cbc CoinUtils Clp OsiClp Osi Cbc OsiCbc CbcSolver cbc-executable)
add_to_folder(test OsiCommonTests)
add_test(Cbc-unitTest bin/cbc -dirSample Cbc/Data/Sample -unitTest)

add_definitions(-DHAVE_CONFIG_H)

# Define FORTRAN_INTEGER_TYPE for Ipopt.
set(FORTRAN_INTEGER_TYPE_C int)

configure(Cbc/CoinUtils src/config.h src/config_coinutils.h)
configure(Cbc/Osi src/Osi/config.h src/Osi/config_osi.h)
configure(Cbc/Clp src/config.h src/config_clp.h)
configure(Cbc/Cgl src/config.h src/config_cgl.h)
configure(Cbc/Cbc src/config.h src/config_cbc.h)

function(add_ipopt_library name path)
  add_coin_library(Ipopt/Ipopt/src/${path} NAME ${name} ${ARGN})
  set_target_properties(${name} PROPERTIES COMPILE_DEFINITIONS IPOPT_BUILD)
  add_to_folder(ipopt ${name})
endfunction()

function(add_ipopt_test name)
  cmake_parse_arguments(add_ipopt_test "" SOURCE_DIR "" ${ARGN})
  if (add_ipopt_test_SOURCE_DIR)
    set(source_dir ${add_ipopt_test_SOURCE_DIR})
  else ()
    set(source_dir ${COIN_SOURCE_DIR})
  endif ()
  set(sources )
  foreach (arg ${add_ipopt_test_UNPARSED_ARGUMENTS})
    set(sources ${sources} ${source_dir}/Ipopt/Ipopt/examples/${name}/${arg})
  endforeach ()
  add_executable(${name} ${sources})
  target_link_libraries(${name} ipopt)
  set_target_properties(${name}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${COIN_BINARY_DIR}/bin)
  add_to_folder(test/ipopt ${name})
endfunction()


# Remove linking to gcc_s, introduced after cmake 3.8.2 for CXX projects
function(removegcc_s)
  foreach(lang ${ARGN})
   string(REPLACE "gcc_s" "" FILTERED "${CMAKE_${lang}_IMPLICIT_LINK_LIBRARIES}")
    set(CMAKE_${lang}_IMPLICIT_LINK_LIBRARIES ${FILTERED} PARENT_SCOPE)
  endforeach()
endfunction()

if (EXISTS ${COIN_SOURCE_DIR}/Ipopt)
  set(BUILD_IPOPT TRUE)
  if (MSVC)
    # Look up gfortran in the MinGW directory. This sets the MINGW_GFORTRAN
    # cache entry if gfortran is found and should be done before including
    # CMakeAddFortranSubdirectory.
    # Find fortran compiler
    if (SIZEOF_INT_P EQUAL 8)
      set(bits 64)
    else ()
      set(bits 32)
    endif ()
  if(${bits} EQUAL 32)
      set(STANDARD_MINGW_INSTALL C:/mingw-w64/i686-8.1.0-posix-dwarf-rt_v6-rev0/mingw32)
    else()
      set(STANDARD_MINGW_INSTALL C:/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64)
    endif()
    set(CHOCO_DEFAULT_INSTALL C:/ProgramData/chocolatey)
    set(ALTERNATIVE_INSTALL C:/mingw${ARCH})
    set(GFORTRANPATHS  ${STANDARD_MINGW_INSTALL}/bin ${CHOCO_DEFAULT_INSTALL}/bin ${ALTERNATIVE_INSTALL}/bin)
    find_program(MINGW_GFORTRAN gfortran PATHS ${GFORTRANPATHS})
    if (MINGW_GFORTRAN)
        # Create a script that sets up a build environment in Fortran subprojects.
      # It is necessary because setenv.cmd produces bogus errors
      # ERROR: The system was unable to find the specified registry key or value.
      # which cause custom build to fail.
      file(WRITE run-setenv.bat "call \"${MP_SETENV}\" %* 2> nul")
      set(path ${CMAKE_CURRENT_BINARY_DIR}/run-setenv.bat)

	  # New versions of the mingw package for chocolatey do not install make where cmake expects it
	  # therefore override it
      find_program(MINGW_MAKE make PATHS ${GFORTRANPATHS})
      file(TO_NATIVE_PATH "${MINGW_MAKE}" MINGW_MAKE_PATH )
      string(REPLACE "\\" "/" MINGW_MAKE_PATH "${MINGW_MAKE_PATH}")
      set(makeline "-DCMAKE_MAKE_PROGRAM:PATH=${MINGW_MAKE_PATH}")

      set(FORTRAN_CMAKE_FLAGS ${makeline} -DCMAKE_GNUtoMS_VCVARS=${path} -DCMAKE_GNUtoMS=ON -DBUILD_SHARED_LIBS=ON -DCMAKE_SH="CMAKE_SH-NOTFOUND")
      set_property(GLOBAL PROPERTY CMAKE_GNU_TO_MSVC 1)
      set(BUILD_IPOPTFORT TRUE)
      message(STATUS "Found GFortran at ${MINGW_GFORTRAN}, using it for VC/gfortran build")
    else()
      message(WARNING "Could not find gfortran, tried:")
      foreach(p ${GFORTRANPATHS})
        message(WARNING ${p})
      endforeach()
    endif ()
  else ()
    enable_language(Fortran OPTIONAL)
    if (CMAKE_Fortran_COMPILER)
      set(BUILD_IPOPTFORT TRUE)
      # Discover and link with gfortran runtime library statically since
      # -static-libgfortran is not always available.
      # NOTE: The value of CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES used to be:
      #   "gfortran;m;quadmath;m;c"
      # but starting from CMake version 3.9.0 it became:
      #   "gfortran;m;gcc_s;gcc;quadmath;m;gcc_s;gcc;c;gcc_s;gcc"
      # and leads to errors such as "/usr/bin/ld: cannot find -lgcc_s" when
      # we use force_static_libs for Fortran. It works fine up to CMake v3.8.2.
      if (NOT BUILD_SHARED_LIBS)
        force_static_libs(Fortran gfortran quadmath)
      endif ()
      if (UNIX)
        set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES
          ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} gcc)
      endif ()
      removegcc_s("CXX")
    endif ()
  endif ()
endif ()

if (MSVC OR NOT BUILD_IPOPTFORT)
  set(F77_FUNC "(name,NAME) name##_")
  set(F77_FUNC_ "${F77_FUNC}")
endif ()

if (BUILD_IPOPT)
  set(HAVE_LINEARSOLVERLOADER 1)
  set(BUILD_LINEARSOLVERLOADER 1)
  set(BUILD_PARDISOLOADER 1)

  check_symbol_exists(drand48 stdlib.h HAVE_DRAND48)
  check_symbol_exists(rand stdlib.h HAVE_RAND)
  check_symbol_exists(std::rand cstdlib HAVE_STD__RAND)
  check_symbol_exists(snprintf stdio.h HAVE_SNPRINTF)
  check_symbol_exists(_snprintf stdio.h HAVE__SNPRINTF)
  check_symbol_exists(va_copy stdarg.h HAVE_VA_COPY)
  check_symbol_exists(vsnprintf "stdio.h;stdarg.h" HAVE_VSNPRINTF)

  string(SUBSTRING ${CMAKE_SHARED_LIBRARY_SUFFIX} 1 -1 SHAREDLIBEXT)

  include_directories(
    Ipopt/Ipopt/src/contrib/LinearSolverLoader
    ThirdParty/Mumps/MUMPS/include ThirdParty/Mumps/MUMPS/libseq)

  if (BUILD_IPOPTFORT)
    # FORTRAN_CMAKE_FLAGS is used to pass additional CMake arguments to
    # Fortran subprojects.
    include(CMakeAddFortranSubdirectory)
    cmake_add_fortran_subdirectory(Ipopt
      PROJECT IPOPTFORT
      ARCHIVE_DIR Ipopt/lib
      RUNTIME_DIR Ipopt/bin
      LIBRARIES ipoptfort
      CMAKE_COMMAND_LINE ${FORTRAN_CMAKE_FLAGS}
      NO_EXTERNAL_INSTALL)
    set(COIN_HAS_BLAS 1)
    set(COIN_HAS_MUMPS 1)
    set(IPOPTFORT ipoptfort)
  endif ()

  # Adds a command to copy external library ${lib} to the binary directory.
  macro (add_copy_library_command lib)
    get_filename_component(libname ${lib} NAME)
    set(dest ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CMAKE_CFG_INTDIR}/${libname})
    add_custom_command(
      OUTPUT ${dest}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${lib} ${dest})
    set(IPOPT_EXTERNAL_LIBS ${IPOPT_EXTERNAL_LIBS} ${dest})
  endmacro ()

  # Copy external libraries.
  if (MSVC)
    add_copy_library_command(${COIN_BINARY_DIR}/Ipopt/bin/libipoptfort.dll)
  endif ()

  add_ipopt_library(ipoptalg Algorithm)
  add_ipopt_library(linsolvers Algorithm/LinearSolvers)
  add_ipopt_library(common Common)
  add_ipopt_library(cgpenalty contrib/CGPenalty)
  add_ipopt_library(linalg LinAlg)
  add_ipopt_library(tmatrices LinAlg/TMatrices)
  add_ipopt_library(ipopt Interfaces
    ipoptalg linsolvers common cgpenalty linalg tmatrices ${IPOPTFORT})

  if (COIN_HAS_ASL)
    add_ipopt_library(ipoptamplinterface Apps/AmplSolver ${ASL_LIBRARIES})
    target_link_libraries(ipopt ipoptamplinterface)
    
    add_executable(ipopt-executable
      Ipopt/Ipopt/src/Apps/AmplSolver/ampl_ipopt.cpp ${IPOPT_EXTERNAL_LIBS})
    target_link_libraries(ipopt-executable ipopt)
    set_target_properties(ipopt-executable PROPERTIES
      COMPILE_DEFINITIONS IPOPT_BUILD OUTPUT_NAME ipopt IMPORT_PREFIX imp)
    add_to_folder(ipopt ipopt-executable)
    if (NOT WIN32)
      add_test(NAME Ipopt-amplTest
        COMMAND sh -c "bin/ipopt Ipopt/Ipopt/test/mytoy.nl | grep 'EXIT: Optimal Solution Found.'")
    endif ()
  endif ()

  if (SIZEOF_INT_P EQUAL 8)
    set(BIT32FCOMMENT C)
  else ()
    set(BIT64FCOMMENT C)
  endif ()
  configure_file(Ipopt/Ipopt/examples/hs071_f/hs071_f.f.in Ipopt/Ipopt/examples/hs071_f/hs071_f.f @ONLY)
  set(srcdir ${CMAKE_CURRENT_SOURCE_DIR})
  configure_file(Ipopt/Ipopt/test/run_unitTests.in Ipopt/Ipopt/test/run_unitTests @ONLY)

  add_ipopt_test(hs071_cpp hs071_main.cpp hs071_nlp.cpp hs071_nlp.hpp)
  add_ipopt_test(hs071_c hs071_c.c)
  if (NOT WIN32)
    add_test(NAME Ipopt-unitTest
      WORKING_DIRECTORY ${COIN_BINARY_DIR}/bin
      COMMAND sh ../Ipopt/Ipopt/test/run_unitTests)
    add_to_folder(ipopt/test Ipopt-unitTest)
  endif ()
  
  configure(Ipopt Ipopt/src/Common/config.h Ipopt/src/Common/config_ipopt.h)
endif()

function(add_coin_test name dir lib)
  set(src )
  foreach (f ${ARGN})
    set(src ${src} ${dir}/${f})
  endforeach ()
  add_executable(${name} ${src})
  target_link_libraries(${name} ${lib})
  add_test(NAME ${name} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/bin/${name}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${dir})
  add_to_folder(test ${name})
endfunction()

if (BUILD_IPOPT AND EXISTS ${COIN_SOURCE_DIR}/Bonmin)
  add_coin_library(Bonmin/Bonmin/src/Algorithms/Branching NAME bonbranching)
  add_coin_library(Bonmin/Bonmin/src/Algorithms/QuadCuts NAME bonquadcuts)
  add_coin_library(Bonmin/Bonmin/src/Algorithms/OaGenerators NAME bonoagenerators)
  add_coin_library(Bonmin/Bonmin/src/Algorithms NAME bonalgorithms)
  add_coin_library(Bonmin/Bonmin/src/Interfaces/Ipopt NAME ipoptinterface)
  add_coin_library(Bonmin/Bonmin/src/Interfaces NAME bonmininterfaces)
  add_coin_library(Bonmin/Bonmin/src/CbcBonmin/Heuristics NAME bonheuristics)
  if (COIN_HAS_ASL)
    add_coin_library(Bonmin/Bonmin/src/Algorithms/Ampl NAME bonamplsetup)
    add_coin_library(Bonmin/Bonmin/src/Interfaces/Ampl NAME bonamplinterface)
    add_dependencies(bonamplinterface ${ASL_LIBRARIES})
    set(BONMIN_DEPS bonamplsetup bonamplinterface ${ASL_LIBRARIES})

    add_executable(bonmin-executable Bonmin/Bonmin/src/Apps/BonMin.cpp)
    set_target_properties(bonmin-executable PROPERTIES
      COMPILE_DEFINITIONS BONMIN_BUILD OUTPUT_NAME bonmin IMPORT_PREFIX imp)
    target_link_libraries(bonmin-executable bonmin)
  endif ()
  add_coin_library(Bonmin/Bonmin/src/CbcBonmin NAME bonmin
    bonalgorithms bonbranching bonquadcuts bonoagenerators
    bonmininterfaces ipoptinterface ipopt bonheuristics CbcSolver ${BONMIN_DEPS})
  target_include_directories(bonmin PUBLIC Bonmin/Bonmin/src/Interfaces/Ampl)

  add_coin_test(Bonmin-unitTest Bonmin/Bonmin/test bonmin InterfaceTest.cpp)
  add_coin_test(Bonmin-CppExample Bonmin/Bonmin/examples/CppExample
    bonmin MyBonmin.cpp MyTMINLP.cpp MyTMINLP.hpp)

  configure(Bonmin Bonmin/src/Interfaces/config.h Bonmin/src/Interfaces/config_bonmin.h)

  add_to_folder(bonmin bonbranching bonquadcuts bonoagenerators bonalgorithms
                       ipoptinterface bonmininterfaces bonheuristics bonamplsetup
                       bonamplinterface bonmin-executable bonmin)
endif ()

function(add_couenne_library name path)
  add_coin_library(Couenne/Couenne/src/${path} NAME ${name})
  set_target_properties(${name} PROPERTIES COMPILE_DEFINITIONS COUENNE_BUILD)
  add_to_folder(couenne ${name})
endfunction()

if (BUILD_IPOPT AND EXISTS ${COIN_SOURCE_DIR}/Couenne)
  include_directories(
    Couenne/Couenne/src Couenne/Couenne/src/main
    Couenne/Couenne/src/expression/operators
    Couenne/Couenne/src/expression/operators/bounds
    Couenne/Couenne/src/expression/partial
    Couenne/Couenne/src/problem/depGraph)
  add_couenne_library(CouenneUtil util)
  add_couenne_library(CouenneExpression expression)
  add_couenne_library(CouenneStandardize standardize)
  add_couenne_library(CouenneProblem problem)
  add_couenne_library(CouenneBoundTightening bound_tightening)
  add_couenne_library(CouenneTwoImplied bound_tightening/twoImpliedBT)
  add_couenne_library(CouenneConvex convex)
  add_couenne_library(CouenneBranch branch)
  add_couenne_library(CouenneDisjunctive disjunctive)
  add_couenne_library(CouenneInterfaces interfaces)
  add_couenne_library(CouenneHeuristics heuristics)
  add_couenne_library(CouenneCrossConv cut/crossconv)
  add_couenne_library(CouenneSdpCuts cut/sdpcuts)
  if (COIN_HAS_ASL)
    set(COUENNEINTERFACE_FROM_ASL 1)
    add_dependencies(CouenneProblem asl)
  endif ()

  set(DEPENDENCY_LINKING 1)
  add_coin_library(Couenne/Couenne/src NAME Couenne
    CouenneInterfaces CouenneHeuristics CouenneProblem CouenneStandardize
    CouenneExpression CouenneBranch CouenneConvex CouenneBoundTightening
    CouenneTwoImplied CouenneUtil CouenneDisjunctive CouenneCrossConv CouenneSdpCuts
    bonmin)
  if (NOT WIN32)
    target_link_libraries(Couenne lapack)
  endif ()

  add_executable(couenne-executable Couenne/Couenne/src/main/BonCouenne.cpp)
  set_target_properties(couenne-executable PROPERTIES
    COMPILE_DEFINITIONS COUENNE_BUILD OUTPUT_NAME couenne IMPORT_PREFIX imp)
  target_link_libraries(couenne-executable Couenne)

  set(srcdir ${CMAKE_CURRENT_SOURCE_DIR}/Couenne/test)
  configure_file(Couenne/Couenne/test/grandTest.sh.in Couenne/Couenne/test/grandTest.sh @ONLY)
  add_test(NAME Couenne-grandTest
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/Couenne/test.cmake
            $<TARGET_FILE:couenne-executable>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Couenne)

  configure(Couenne Couenne/src/config.h Couenne/src/config_couenne.h)
    add_to_folder(couenne Couenne couenne-executable) 

endif ()

# Fix paths to libstdc++ on OS X.
if (APPLE)
  set(libname libstdc++.6.dylib)
  foreach (solver ipopt bonmin couenne)
    set(target ${solver}-executable)
    if (TARGET ${target})
      add_custom_command(
        TARGET ${target} POST_BUILD
        COMMAND install_name_tool -change /opt/local/lib/libgcc/${libname}
          /usr/lib/${libname} $<TARGET_FILE:${target}>)
    endif ()
  endforeach ()
endif ()

if (COIN_HAS_ASL)
  # Generate a file containing coin project versions.
  set(COIN_VERSIONS_FILE coin-versions.txt)
  add_custom_target(coin-versions cmake -E remove ${COIN_VERSIONS_FILE})
  foreach (solver cbc ipopt bonmin couenne)
    if (TARGET ${solver}-executable)
      add_custom_command(TARGET coin-versions POST_BUILD
        COMMAND $<TARGET_FILE:${solver}-executable> -v -AMPL >> ${COIN_VERSIONS_FILE})
      add_dependencies(coin-versions ${solver}-executable)
    endif ()
  endforeach ()
endif ()

install(TARGETS ${COIN_INSTALL_TARGETS}
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib)
