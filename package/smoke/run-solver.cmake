# Run one solver for the smoke tests (see CMakeLists.txt in this directory).
#
#   cmake -DSOLVER=<exe> "-DARGS=a;b" -DWORKING_DIR=<dir> \
#         [-DSOLUTION=<.sol file that must report an optimal solution>] -P run-solver.cmake
cmake_minimum_required(VERSION 3.10)

if (SOLUTION)
  file(REMOVE "${SOLUTION}")
endif ()
execute_process(
  COMMAND "${SOLVER}" ${ARGS}
  WORKING_DIRECTORY "${WORKING_DIR}"
  RESULT_VARIABLE result
  OUTPUT_VARIABLE output
  ERROR_VARIABLE output)
message("${output}")
if (NOT result EQUAL 0)
  message(FATAL_ERROR "${SOLVER} failed (${result})")
endif ()

if (SOLUTION)
  if (NOT EXISTS "${SOLUTION}")
    message(FATAL_ERROR "${SOLVER} wrote no ${SOLUTION}")
  endif ()
  file(READ "${SOLUTION}" solution)
  if (NOT solution MATCHES "[Oo]ptimal")
    message(FATAL_ERROR "${SOLVER} did not report an optimal solution:\n${solution}")
  endif ()
endif ()
