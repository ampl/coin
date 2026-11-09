# Smoke tests for the packaged solvers, run with CTest so that the results can
# be published as JUnit (Azure Pipelines "Tests" tab).
#
#   cmake -S package/smoke -B <build> -DPLATFORM=<platform> \
#         -DARCHIVE_DIR=<dir with the archives> \
#         -DSOLVERS=cbc,bonmin,couenne,ipopt [-DPREFIX=signed-]
#   ctest --test-dir <build> --output-on-failure --output-junit junit.xml
#
# Each solver archive is unpacked at configure time and becomes one test that
# runs the solver from the unpacked copy. Testing the archives, not the build
# tree, also catches files missing from a package (e.g. libipoptfort.dll).
cmake_minimum_required(VERSION 3.20)
project(solver-smoke-tests NONE)
enable_testing()

foreach (var PLATFORM ARCHIVE_DIR SOLVERS)
  if (NOT DEFINED ${var})
    message(FATAL_ERROR "-D${var}=... is required")
  endif ()
endforeach ()
string(REPLACE "," ";" SOLVERS "${SOLVERS}")
get_filename_component(ARCHIVE_DIR "${ARCHIVE_DIR}" ABSOLUTE)

set(exe "")
set(archive_extension "tar.gz")
if (PLATFORM MATCHES "^mswin")
  set(exe ".exe")
  set(archive_extension "zip")
endif ()

set(model "${CMAKE_CURRENT_LIST_DIR}/../test/camel1.nl")
if (NOT EXISTS "${model}")
  message(FATAL_ERROR "Missing ${model}")
endif ()

foreach (solver ${SOLVERS})
  file(GLOB archives
    "${ARCHIVE_DIR}/${PREFIX}solvers-${solver}.${PLATFORM}.*.${archive_extension}")
  list(LENGTH archives count)
  if (NOT count EQUAL 1)
    message(FATAL_ERROR "Expected one ${PREFIX}solvers-${solver} archive, found: ${archives}")
  endif ()
  set(dir "${CMAKE_BINARY_DIR}/${solver}")
  file(REMOVE_RECURSE "${dir}")
  file(MAKE_DIRECTORY "${dir}")
  # `cmake -E tar` keeps the executable bit, unlike file(ARCHIVE_EXTRACT).
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf "${archives}"
    WORKING_DIRECTORY "${dir}" RESULT_VARIABLE result)
  if (NOT result EQUAL 0)
    message(FATAL_ERROR "Cannot unpack ${archives}")
  endif ()
  # cbc can't read the nonlinear model, so it is only started with -v.
  if (solver STREQUAL "cbc")
    set(args -v -AMPL)
    set(check "")
  else ()
    file(COPY "${model}" DESTINATION "${dir}")
    set(args camel1 -AMPL)
    set(check "${dir}/camel1.sol")
  endif ()
  add_test(NAME ${solver}
    COMMAND ${CMAKE_COMMAND}
      "-DSOLVER=${dir}/${solver}${exe}" "-DARGS=${args}"
      "-DWORKING_DIR=${dir}" "-DSOLUTION=${check}"
      -P "${CMAKE_CURRENT_LIST_DIR}/run-solver.cmake")
  set_tests_properties(${solver} PROPERTIES TIMEOUT 120)
endforeach ()
