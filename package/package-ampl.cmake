# Create individual solver archives from existing binaries.
# Called by the ampl-module target and package-solvers.sh.
cmake_minimum_required(VERSION 3.10)

foreach (var PLATFORM BINDIR OUTDIR SOURCE_DIR VERSIONS_FILE)
  if (NOT DEFINED ${var})
    message(FATAL_ERROR "package-ampl.cmake: -D${var}=... is required")
  endif ()
  if (NOT var STREQUAL "PLATFORM")
    get_filename_component(${var} "${${var}}" ABSOLUTE)
  endif ()
endforeach ()
if (NOT DEFINED SOLVERS)
  set(SOLVERS cbc,bonmin,couenne,ipopt)
endif ()
string(REPLACE "," ";" SOLVERS "${SOLVERS}")
if (NOT DEFINED PACKAGE_DATE)
  string(TIMESTAMP PACKAGE_DATE "%Y%m%d")
endif ()

set(exe "")
if (PLATFORM MATCHES "^mswin")
  set(exe ".exe")
endif ()

# Validate all inputs before writing any release archives.
set(required "${VERSIONS_FILE}" "${SOURCE_DIR}/Cbc/Cbc/LICENSE")
foreach (solver ${SOLVERS})
  list(APPEND required "${BINDIR}/${solver}${exe}"
    "${SOURCE_DIR}/docs/README.${solver}.txt"
    "${SOURCE_DIR}/docs/CHANGES.${solver}.md")
  if (exe AND NOT solver STREQUAL "cbc")
    list(APPEND required "${BINDIR}/libipoptfort.dll")
  endif ()
endforeach ()
foreach (path ${required})
  if (NOT EXISTS "${path}")
    message(FATAL_ERROR "Missing ${path}")
  endif ()
endforeach ()

file(MAKE_DIRECTORY "${OUTDIR}")
foreach (solver ${SOLVERS})
  set(stage "${OUTDIR}/.stage-${solver}")
  file(REMOVE_RECURSE "${stage}")
  file(MAKE_DIRECTORY "${stage}/docs")
  # Preserve executable permissions; distribute documentation as readable files.
  file(COPY "${BINDIR}/${solver}${exe}" DESTINATION "${stage}")
  if (exe AND NOT solver STREQUAL "cbc")
    file(COPY "${BINDIR}/libipoptfort.dll" DESTINATION "${stage}")
  endif ()
  file(COPY "${SOURCE_DIR}/docs/README.${solver}.txt"
    "${SOURCE_DIR}/docs/CHANGES.${solver}.md" DESTINATION "${stage}/docs"
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
  configure_file("${SOURCE_DIR}/Cbc/Cbc/LICENSE" "${stage}/coin-license.txt" COPYONLY)
  configure_file("${VERSIONS_FILE}" "${stage}/coin-versions.txt" COPYONLY)

  set(archive "${OUTDIR}/solvers-${solver}.${PLATFORM}.${PACKAGE_DATE}.tar.gz")
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E env COPYFILE_DISABLE=1
      ${CMAKE_COMMAND} -E tar czf "${archive}" --format=gnutar .
    WORKING_DIRECTORY "${stage}"
    RESULT_VARIABLE result)
  file(REMOVE_RECURSE "${stage}")
  if (NOT result EQUAL 0)
    message(FATAL_ERROR "Failed to create ${archive}")
  endif ()
  message(STATUS "Created ${archive}")
endforeach ()
