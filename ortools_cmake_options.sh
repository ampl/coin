#!/bin/bash
# CMake options shared by all the OR-Tools builds (sourced by do_build_ortools.sh).
#
# We only need CP-SAT, built as static libraries with all the dependencies
# (abseil, protobuf, re2, zlib, ...) built from source and installed with it:
#  - BUILD_DEPS=ON, INSTALL_BUILD_DEPS=ON: dependencies built and installed
#    together with OR-Tools, so the package is self-contained;
#  - BUILD_SHARED_LIBS=OFF: static libraries, no DLLs/shared libs to distribute;
#  - other solvers off (SCIP, HiGHS, COIN-OR (CBC, CLP), GLPK, PDLP, Gurobi,
#    CPLEX), also MathOpt and FlatZinc: none is used by CP-SAT, and some
#    are already distributed (with other versions) by AMPL.
#    Gurobi cannot be left on: it needs MathOpt. With it off, v9.15 still
#    references the Gurobi functions unconditionally; see
#    patches/ortools-guard-gurobi.patch.
ORTOOLS_CMAKE_OPTIONS="\
 -DBUILD_DEPS=ON \
 -DINSTALL_BUILD_DEPS=ON \
 -DBUILD_SHARED_LIBS=OFF \
 -DBUILD_CXX=ON \
 -DBUILD_SAMPLES=OFF \
 -DBUILD_EXAMPLES=OFF \
 -DBUILD_FLATZINC=OFF \
 -DBUILD_MATH_OPT=OFF \
 -DBUILD_TESTING=OFF \
 -DUSE_SCIP=OFF \
 -DUSE_HIGHS=OFF \
 -DUSE_COINOR=OFF \
 -DUSE_GLPK=OFF \
 -DUSE_PDLP=OFF \
 -DUSE_GUROBI=OFF \
 -DUSE_CPLEX=OFF"
