# ASL-free builds of the Cbc and CbcSolver libraries.
#
# The top-level build defines COIN_HAS_ASL for every target, which compiles
# the AMPL interface (Cbc_ampl.cpp, CbcLinkedUtils.cpp and parts of
# CbcSolver.cpp) into Cbc and CbcSolver and makes anything linking them also
# need ASL. That is what the cbc executable wants, but not a driver that
# brings its own model reader (such as the MP-based cbc driver), so the
# libraries shipped in the cbc-libs package are built here without it:
#
#   Cbc-noasl, CbcSolver-noasl -> Cbc / CbcSolver (what cbc-libs packages)
#   Cbc, CbcSolver             -> CbcAmpl / CbcSolverAmpl (cbc, bonmin, couenne)

get_directory_property(defs COMPILE_DEFINITIONS)
list(REMOVE_ITEM defs COIN_HAS_ASL)
set_directory_properties(PROPERTIES COMPILE_DEFINITIONS "${defs}")

# Adds <target>-noasl: <target> rebuilt in this directory, linked against the
# -noasl variants of its dependencies and without ASL.
function (add_noasl_library target)
  set(name ${target}-noasl)
  get_target_property(type ${target} TYPE)
  string(REPLACE "_LIBRARY" "" type ${type})
  get_target_property(source_dir ${target} SOURCE_DIR)
  get_target_property(target_sources ${target} SOURCES)
  set(sources )
  foreach (source ${target_sources})
    if (NOT source MATCHES "^\\$<" AND NOT IS_ABSOLUTE "${source}")
      set(source "${source_dir}/${source}")
    endif ()
    list(APPEND sources "${source}")
  endforeach ()
  add_library(${name} ${type} ${sources})

  foreach (property COMPILE_DEFINITIONS COMPILE_OPTIONS INCLUDE_DIRECTORIES)
    get_target_property(value ${target} ${property})
    if (value)
      set_property(TARGET ${name} PROPERTY ${property} "${value}")
    endif ()
  endforeach ()

  get_target_property(deps ${target} LINK_LIBRARIES)
  set(libs )
  foreach (dep ${deps})
    list(FIND ASL_LIBRARIES "${dep}" asl_index)
    if (TARGET ${dep}-noasl)
      list(APPEND libs ${dep}-noasl)
    elseif (asl_index EQUAL -1)
      list(APPEND libs ${dep})
    endif ()
  endforeach ()
  target_link_libraries(${name} ${libs})

  # Take the plain library name and put it where the top-level libraries go
  # (the build root, or <build root>/<config> with multi-config generators),
  # which is where package-solvers.sh and the macOS scripts look for them.
  set_target_properties(${name} PROPERTIES
    OUTPUT_NAME ${target}
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
    FOLDER cbc)
  set_target_properties(${target} PROPERTIES OUTPUT_NAME ${target}Ampl)
endfunction ()

add_noasl_library(Cbc)
add_noasl_library(CbcSolver)

install(TARGETS Cbc-noasl CbcSolver-noasl
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib)

# Linking a program against the ASL-free libraries alone (no asl) fails if
# any ASL code is still in them.
set(test_dir ${COIN_SOURCE_DIR}/Cbc/Cbc/test)
add_executable(Cbc-noasl-CInterfaceTest ${test_dir}/CInterfaceTest.c)
set_target_properties(Cbc-noasl-CInterfaceTest PROPERTIES
  FOLDER test LINKER_LANGUAGE CXX)
target_link_libraries(Cbc-noasl-CInterfaceTest CbcSolver-noasl)
add_test(NAME Cbc-noasl-CInterfaceTest
  COMMAND Cbc-noasl-CInterfaceTest
  WORKING_DIRECTORY ${COIN_SOURCE_DIR}/Cbc/Cbc/test)
