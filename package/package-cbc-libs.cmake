# Create the archive with the static libraries needed to link a Cbc driver.
# Called by package-solvers.sh next to package-ampl.cmake.
#
#   cmake -DPLATFORM=<platform> -DLIBDIR=<build dir holding the libs> \
#         -DOUTDIR=<dir> -DSOURCE_DIR=<repo> -DVERSIONS_FILE=<coin-versions.txt> \
#         -P package-cbc-libs.cmake
#
# On Windows the libraries are cbcsolver.lib, cbc.lib, ... (plus
# pthreads.lib); elsewhere they are libCbc.a, libCbcSolver.a, ...
cmake_minimum_required(VERSION 3.10)

foreach (var PLATFORM LIBDIR OUTDIR SOURCE_DIR VERSIONS_FILE)
  if (NOT DEFINED ${var})
    message(FATAL_ERROR "package-cbc-libs.cmake: -D${var}=... is required")
  endif ()
  if (NOT var STREQUAL "PLATFORM")
    get_filename_component(${var} "${${var}}" ABSOLUTE)
  endif ()
endforeach ()
if (PLATFORM MATCHES "^macos")
  set(PLATFORM macos64)
endif ()

# Libraries a Cbc driver links against, in link order.
set(cbc_libs CbcSolver Cbc CoinUtils Osi OsiClp Clp Cgl)

set(archive_extension "tar.gz")
if (PLATFORM MATCHES "^mswin")
  set(archive_extension "zip")
endif ()

# Release versions come from the solver changelogs: use the newest one.
set(package_version "")
file(GLOB changelogs "${SOURCE_DIR}/docs/CHANGES.*.md")
foreach (changelog ${changelogs})
  file(READ "${changelog}" contents)
  if (NOT contents MATCHES "## ([0-9]+)[\r\n]")
    message(FATAL_ERROR "Missing dated release heading in ${changelog}")
  endif ()
  if (CMAKE_MATCH_1 STRGREATER package_version)
    set(package_version "${CMAKE_MATCH_1}")
  endif ()
endforeach ()
if (NOT package_version)
  message(FATAL_ERROR "No solver changelog release versions found")
endif ()

# Map each library to its file in LIBDIR and its name in the package.
set(sources "")
set(names "")
foreach (lib ${cbc_libs})
  if (PLATFORM MATCHES "^mswin")
    string(TOLOWER "${lib}" lower)
    list(APPEND sources "${LIBDIR}/${lib}.lib")
    list(APPEND names "${lower}.lib")
  else ()
    list(APPEND sources "${LIBDIR}/lib${lib}.a")
    list(APPEND names "lib${lib}.a")
  endif ()
endforeach ()
if (PLATFORM MATCHES "^mswin")
  list(APPEND sources "${LIBDIR}/pthreads.lib")
  list(APPEND names "pthreads.lib")
endif ()

foreach (path ${sources} "${VERSIONS_FILE}" "${SOURCE_DIR}/Cbc/Cbc/LICENSE")
  if (NOT EXISTS "${path}")
    message(FATAL_ERROR "Missing ${path}")
  endif ()
endforeach ()

file(MAKE_DIRECTORY "${OUTDIR}")
set(stage "${OUTDIR}/.stage-cbc-libs")
file(REMOVE_RECURSE "${stage}")
file(MAKE_DIRECTORY "${stage}/lib" "${stage}/docs")
list(LENGTH sources count)
math(EXPR last "${count} - 1")
foreach (i RANGE ${last})
  list(GET sources ${i} source)
  list(GET names ${i} name)
  configure_file("${source}" "${stage}/lib/${name}" COPYONLY)
endforeach ()
configure_file("${SOURCE_DIR}/Cbc/Cbc/LICENSE" "${stage}/docs/coin-license.txt" COPYONLY)
configure_file("${VERSIONS_FILE}" "${stage}/docs/coin-versions.txt" COPYONLY)

set(archive "${OUTDIR}/cbc-libs.${PLATFORM}.${package_version}.${archive_extension}")
if (archive_extension STREQUAL "zip")
  set(archive_options cf "${archive}" --format=zip)
else ()
  set(archive_options czf "${archive}" --format=gnutar)
endif ()
execute_process(
  COMMAND ${CMAKE_COMMAND} -E env COPYFILE_DISABLE=1
    ${CMAKE_COMMAND} -E tar ${archive_options} .
  WORKING_DIRECTORY "${stage}"
  RESULT_VARIABLE result)
file(REMOVE_RECURSE "${stage}")
if (NOT result EQUAL 0)
  message(FATAL_ERROR "Failed to create ${archive}")
endif ()
message(STATUS "Created ${archive}")
