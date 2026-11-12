# Standalone builds of the libraries in the cbc-libs package, for drivers
# that bring their own model reader (such as the MP-based cbc driver).
#
# The top-level build defines COIN_HAS_ASL for every target, which compiles
# the AMPL interface (Cbc_ampl.cpp, CbcLinkedUtils.cpp and parts of
# CbcSolver.cpp) into Cbc and CbcSolver, and builds CoinUtils with gzip/bzip2
# support for compressed MPS files when zlib/bzip2 are found. That is what
# the executables want, but it makes anything linking the libraries also need
# ASL, zlib and bzip2. The libraries that differ are rebuilt here without
# them, and the package takes these:
#
#   Cbc-standalone, CbcSolver-standalone, CoinUtils-standalone
#       -> Cbc / CbcSolver / CoinUtils (what cbc-libs packages)
#   Cbc, CbcSolver, CoinUtils
#       -> CbcAmpl / CbcSolverAmpl / CoinUtilsAmpl (cbc, bonmin, couenne)
#
# Osi, OsiClp, Clp and Cgl use none of these, so the package takes them as
# they are.

# The package's public config headers are generated into the source tree when
# configuring and aren't checked in. Keep a copy in the build directory, where
# package-solvers.sh (and the macOS stage/lipo scripts) find them even when
# the build ran somewhere else (in a container, or another macOS job).
set(config_dir ${PROJECT_BINARY_DIR}/cbc-libs-include)
file(MAKE_DIRECTORY ${config_dir})
foreach (header Cbc/CoinUtils/src/config_coinutils.h Cbc/Osi/src/Osi/config_osi.h
                Cbc/Clp/src/config_clp.h Cbc/Cgl/src/config_cgl.h
                Cbc/Cbc/src/config_cbc.h)
  configure_file(${COIN_SOURCE_DIR}/${header} ${config_dir} COPYONLY)
endforeach ()

get_directory_property(defs COMPILE_DEFINITIONS)
list(REMOVE_ITEM defs COIN_HAS_ASL)
set_directory_properties(PROPERTIES COMPILE_DEFINITIONS "${defs}")

# Dependencies and definitions that standalone libraries drop.
set(excluded_libs ${ASL_LIBRARIES} ${COINUTILS_LIBS})
set(excluded_defs COIN_HAS_ASL COIN_HAS_ZLIB COIN_HAS_BZLIB)

# Adds <target>-standalone: <target> rebuilt in this directory without the
# excluded definitions and dependencies, linked against the -standalone
# variants of its other dependencies.
function (add_standalone_library target)
  set(name ${target}-standalone)
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
      if (property STREQUAL "COMPILE_DEFINITIONS")
        list(REMOVE_ITEM value ${excluded_defs})
      endif ()
      set_property(TARGET ${name} PROPERTY ${property} "${value}")
    endif ()
  endforeach ()

  get_target_property(deps ${target} LINK_LIBRARIES)
  set(libs )
  foreach (dep ${deps})
    list(FIND excluded_libs "${dep}" excluded_index)
    if (TARGET ${dep}-standalone)
      list(APPEND libs ${dep}-standalone)
    elseif (excluded_index EQUAL -1)
      list(APPEND libs ${dep})
    endif ()
  endforeach ()
  if (libs)
    target_link_libraries(${name} ${libs})
  endif ()

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

add_standalone_library(CoinUtils)
add_standalone_library(Cbc)
add_standalone_library(CbcSolver)

install(TARGETS CoinUtils-standalone Cbc-standalone CbcSolver-standalone
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib)

# Link a program against exactly the packaged libraries plus the system
# thread libraries, the way a driver does: it fails to link if any of them
# still needs ASL, zlib or bzip2. The libraries are given as files, not
# targets, so CMake doesn't add the regular libraries' dependencies.
set(packaged_libs )
foreach (lib CbcSolver Cbc Cgl OsiClp Clp Osi CoinUtils)
  if (TARGET ${lib}-standalone)
    set(lib ${lib}-standalone)
  endif ()
  list(APPEND packaged_libs ${lib})
endforeach ()
set(test_dir ${COIN_SOURCE_DIR}/Cbc/Cbc/test)
add_executable(cbc-libs-CInterfaceTest ${test_dir}/CInterfaceTest.c)
set_target_properties(cbc-libs-CInterfaceTest PROPERTIES
  FOLDER test LINKER_LANGUAGE CXX)
add_dependencies(cbc-libs-CInterfaceTest ${packaged_libs})
foreach (lib ${packaged_libs})
  target_link_libraries(cbc-libs-CInterfaceTest $<TARGET_LINKER_FILE:${lib}>)
endforeach ()
target_link_libraries(cbc-libs-CInterfaceTest ${CBC_PTHREADS_LIB} ${LIB_RT})
add_test(NAME cbc-libs-CInterfaceTest
  COMMAND cbc-libs-CInterfaceTest
  WORKING_DIRECTORY ${COIN_SOURCE_DIR}/Cbc/Cbc/test)
