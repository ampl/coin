# Builds the AMPL solver-distribution module "coin" from built binaries.
# Run by the "ampl-module" target (see CMakeLists.txt), or by hand:
#
#   cmake -DPLATFORM=linux64 -DBINDIR=build/bin -DOUTDIR=out -DSOURCE_DIR=. \
#         -P package-ampl.cmake
#
# Writes coin-module.<platform>.zip (mswin64) or
# coin-module.<platform>.tar.gz (everything else), containing:
#
#   ampl.<platform>/
#     docs/README.<solver>.txt, docs/CHANGES.<solver>.md   (per solver)
#     docs/coin-license.txt
#     <solver executables>, [libipoptfort.dll on Windows]
#
# README/CHANGES are copied as-is from docs/ in the source tree. metadata.json is
# not written here: it is generated automatically elsewhere.
#
# Optional: -DNAME=<module name> (default coin),
#           -DSOLVERS=<comma-separated list> (default bonmin,couenne,ipopt).
cmake_minimum_required(VERSION 3.10)

foreach (var PLATFORM BINDIR OUTDIR SOURCE_DIR)
  if (NOT DEFINED ${var})
    message(FATAL_ERROR "package-ampl.cmake: -D${var}=... is required")
  endif ()
endforeach ()
if (NOT DEFINED NAME)
  set(NAME coin)
endif ()
if (NOT DEFINED SOLVERS)
  set(SOLVERS bonmin,couenne,ipopt)
endif ()
string(REPLACE "," ";" SOLVERS "${SOLVERS}")

if (PLATFORM STREQUAL "mswin64")
  set(exe ".exe")
  set(archive_ext "zip")
elseif (PLATFORM MATCHES "^(linux64|linuxarm64|macos64|macosarm64)$")
  set(exe "")
  set(archive_ext "tar.gz")
else ()
  message(FATAL_ERROR "Unknown platform '${PLATFORM}'")
endif ()

get_filename_component(OUTDIR "${OUTDIR}" ABSOLUTE)
set(work "${OUTDIR}/.stage")
set(stage "${work}/ampl.${PLATFORM}")
file(REMOVE_RECURSE "${work}")
file(MAKE_DIRECTORY "${stage}/docs")

# Binaries (file(COPY) keeps the executable bit). On Windows the Ipopt-based
# solvers share libipoptfort.dll.
set(binaries "")
foreach (solver ${SOLVERS})
  list(APPEND binaries "${solver}${exe}")
endforeach ()
if (PLATFORM STREQUAL "mswin64")
  list(APPEND binaries libipoptfort.dll)
endif ()
foreach (f ${binaries})
  if (NOT EXISTS "${BINDIR}/${f}")
    message(FATAL_ERROR "Missing ${BINDIR}/${f}")
  endif ()
  file(COPY "${BINDIR}/${f}" DESTINATION "${stage}")
endforeach ()

# Docs. README and CHANGES are required for every solver so a release can't
# ship without them.
foreach (solver ${SOLVERS})
  set(readme "${SOURCE_DIR}/docs/README.${solver}.txt")
  set(changes "${SOURCE_DIR}/docs/CHANGES.${solver}.md")
  foreach (f "${readme}" "${changes}")
    if (NOT EXISTS "${f}")
      message(FATAL_ERROR "Missing ${f}")
    endif ()
  endforeach ()
  file(COPY "${readme}" "${changes}" DESTINATION "${stage}/docs"
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
endforeach ()

if (NOT EXISTS "${SOURCE_DIR}/Cbc/Cbc/LICENSE")
  message(FATAL_ERROR "Missing ${SOURCE_DIR}/Cbc/Cbc/LICENSE")
endif ()
file(COPY "${SOURCE_DIR}/Cbc/Cbc/LICENSE" DESTINATION "${stage}/docs"
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
file(RENAME "${stage}/docs/LICENSE" "${stage}/docs/coin-license.txt")

# Archive: zip for Windows, tar.gz elsewhere.
set(archive "${OUTDIR}/${NAME}-module.${PLATFORM}.${archive_ext}")
file(REMOVE "${archive}")
if (archive_ext STREQUAL "zip")
  set(tar_args cf "${archive}" --format=zip)
else ()
  set(tar_args czf "${archive}")
endif ()
execute_process(
  COMMAND ${CMAKE_COMMAND} -E tar ${tar_args} "ampl.${PLATFORM}"
  WORKING_DIRECTORY "${work}"
  RESULT_VARIABLE result)
if (NOT result EQUAL 0)
  message(FATAL_ERROR "Failed to create ${archive}")
endif ()
file(REMOVE_RECURSE "${work}")
message(STATUS "Created ${archive}")
