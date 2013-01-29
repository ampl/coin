# CMake build script for CBC.
# Author: Victor Zverovich

cmake_minimum_required(VERSION 2.6)

project(CBC)

macro(add_coin_library name path)
  set(sources )
  foreach (arg ${ARGN})
    set(sources ${sources} ${path}/${arg})
  endforeach ()
  include_directories(${path})
  add_library(${name} ${sources})
endmacro()

macro(add_coin_library_auto path)
  get_filename_component(name ${path} NAME)
  set(sources )
  file(STRINGS ${path}/Makefile.am LINES)
  set(continued False)
  foreach (line ${LINES})
    if (continued)
      # Do nothing.
    elseif (line MATCHES "^lib${name}_la_SOURCES [+]?= (.*)")
      set(line ${CMAKE_MATCH_1})
    else ()
      set(line )
    endif ()
    if (line)
      if (line MATCHES "^(.*)\\$")
        set(line ${CMAKE_MATCH_1})
        set(continued True)
      else ()
        set(continued False)
      endif ()
      string(REGEX MATCHALL "[^ \t]+" src ${line})
      set(sources ${sources} ${src})
    endif ()
  endforeach ()
  add_coin_library(${name} ${path} ${sources})
  set(dependencies )
  foreach (arg ${ARGN})
    set(dependencies ${dependencies} ${arg})
  endforeach ()
  target_link_libraries(${name} ${dependencies})
endmacro()

macro(add_coin_package name)
  add_coin_library(${name} ${name}/src ${ARGN})
  string(TOUPPER ${name} uppercase_name)
  set(COIN_HAS_${uppercase_name} 1)
  add_definitions(-D${uppercase_name}_BUILD)
endmacro()

# Process config.h using autoconf rules.
macro(configure file)
  file(STRINGS ${file}.in CONFIG)
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
        if (NOT (var MATCHES ^HAVE OR ${var} EQUAL 1))
          set(value \"${${var}}\")
        elseif (${var})
          set(value 1)
        else ()
          set(value 0)
        endif ()
        set(line "#${space}define ${var} ${value}")
       endif ()
    endif ()
    string(REPLACE "<semi>" ";" line "${line}")
    set(CONFIG_OUT "${CONFIG_OUT}${line}\n")
  endforeach ()
  file(WRITE ${file}
"/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

${CONFIG_OUT}")
endmacro()

include(CheckIncludeFiles)
check_include_files(math.h HAVE_MATH_H)

include(CheckSymbolExists)
include(CheckTypeSize)
check_symbol_exists(isfinite math.h COIN_HAS_ISFINITE)
if (COIN_HAS_ISFINITE)
  set(COIN_C_FINITE_C isfinite)
endif ()
check_symbol_exists(int64_t inttypes.h COIN_HAS_INT64_T)
if (COIN_HAS_INT64_T)
  set(COIN_INT64_T_C int64_t)
  set(COIN_UINT64_T_C uint64_t)
else ()
  check_type_size("long long" LONG_LONG_SIZE)
  if (LONG_LONG_SIZE EQUAL 8)
    set(COIN_INT64_T_C "long long")
    set(COIN_UINT64_T_C "unsigned long long")
  endif ()
endif ()

# Get version numbers and parts of config.h from configure.ac.
file(READ configure.ac LINES)
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
    set(PACKAGE ${CMAKE_MATCH_1})
    set(VERSION ${CMAKE_MATCH_2})
    set(PACKAGE_BUGREPORT ${CMAKE_MATCH_3})
    foreach (var VERSION PACKAGE PACKAGE_BUGREPORT)
      if (${${var}} MATCHES "^\\[(.*)\\]$")
        set(${var} ${CMAKE_MATCH_1})
      endif ()
    endforeach ()
    message(STATUS "Got VERSION=${VERSION} from configure.ac")
  endif ()
endforeach ()
set(PACKAGE_NAME ${PACKAGE})
string(TOLOWER ${PACKAGE} PACKAGE_TARNAME)
set(PACKAGE_URL "")
set(PACKAGE_VERSION ${VERSION})
string(TOUPPER ${PACKAGE} UPPERCASE_PACKAGE)
set(${UPPERCASE_PACKAGE}_VERSION ${VERSION})
set(PACKAGE_STRING "${PACKAGE} ${VERSION}")

add_coin_package(CoinUtils
  config_coinutils.h
  CoinUtilsConfig.h
  Coin_C_defines.h
  CoinAlloc.cpp CoinAlloc.hpp
  CoinBuild.cpp CoinBuild.hpp
  CoinDenseVector.cpp CoinDenseVector.hpp
  CoinDistance.hpp
  CoinError.cpp CoinError.hpp
  CoinFactorization.hpp
  CoinFactorization1.cpp
  CoinFactorization2.cpp
  CoinFactorization3.cpp
  CoinFactorization4.cpp
  CoinSimpFactorization.hpp
  CoinSimpFactorization.cpp
  CoinDenseFactorization.hpp
  CoinDenseFactorization.cpp
  CoinOslFactorization.hpp
  CoinOslFactorization.cpp
  CoinOslFactorization2.cpp
  CoinOslFactorization3.cpp
  CoinOslC.h
  CoinFileIO.cpp CoinFileIO.hpp
  CoinFinite.cpp CoinFinite.hpp
  CoinFloatEqual.hpp
  CoinHelperFunctions.hpp
  CoinIndexedVector.cpp CoinIndexedVector.hpp
  CoinLpIO.cpp CoinLpIO.hpp
  CoinMessage.cpp CoinMessage.hpp
  CoinMessageHandler.cpp CoinMessageHandler.hpp
  CoinModel.cpp CoinModel.hpp
  CoinStructuredModel.cpp CoinStructuredModel.hpp
  CoinModelUseful.cpp CoinModelUseful.hpp
  CoinModelUseful2.cpp
  CoinMpsIO.cpp CoinMpsIO.hpp
  CoinPackedMatrix.cpp CoinPackedMatrix.hpp
  CoinPackedVector.cpp CoinPackedVector.hpp
  CoinPackedVectorBase.cpp CoinPackedVectorBase.hpp
  CoinParam.cpp CoinParamUtils.cpp CoinParam.hpp
  CoinPostsolveMatrix.cpp
  CoinPragma.hpp
  CoinPrePostsolveMatrix.cpp
  CoinPresolveDoubleton.cpp CoinPresolveDoubleton.hpp
  CoinPresolveDual.cpp CoinPresolveDual.hpp
  CoinPresolveDupcol.cpp CoinPresolveDupcol.hpp
  CoinPresolveEmpty.cpp CoinPresolveEmpty.hpp
  CoinPresolveFixed.cpp CoinPresolveFixed.hpp
  CoinPresolveForcing.cpp CoinPresolveForcing.hpp
  CoinPresolveHelperFunctions.cpp
  CoinPresolveImpliedFree.cpp CoinPresolveImpliedFree.hpp
  CoinPresolveIsolated.cpp CoinPresolveIsolated.hpp
  CoinPresolveMatrix.cpp CoinPresolveMatrix.hpp
  CoinPresolvePsdebug.cpp CoinPresolvePsdebug.hpp
  CoinPresolveSingleton.cpp CoinPresolveSingleton.hpp
  CoinPresolveSubst.cpp CoinPresolveSubst.hpp
  CoinPresolveTighten.cpp CoinPresolveTighten.hpp
  CoinPresolveTripleton.cpp CoinPresolveTripleton.hpp
  CoinPresolveUseless.cpp CoinPresolveUseless.hpp
  CoinPresolveZeros.cpp CoinPresolveZeros.hpp
  CoinSearchTree.cpp CoinSearchTree.hpp
  CoinShallowPackedVector.cpp CoinShallowPackedVector.hpp
  CoinSignal.hpp
  CoinSmartPtr.hpp
  CoinSnapshot.cpp CoinSnapshot.hpp
  CoinSort.hpp
  CoinTime.hpp
  CoinTypes.hpp
  CoinUtility.hpp
  CoinWarmStart.hpp
  CoinWarmStartBasis.cpp CoinWarmStartBasis.hpp
  CoinWarmStartVector.cpp CoinWarmStartVector.hpp
  CoinWarmStartDual.cpp CoinWarmStartDual.hpp
  CoinWarmStartPrimalDual.cpp CoinWarmStartPrimalDual.hpp)

add_coin_package(Clp
  ClpConfig.h
  ClpCholeskyBase.cpp ClpCholeskyBase.hpp
  ClpCholeskyDense.cpp ClpCholeskyDense.hpp
  ClpConstraint.cpp ClpConstraint.hpp
  ClpConstraintLinear.cpp ClpConstraintLinear.hpp
  ClpConstraintQuadratic.cpp ClpConstraintQuadratic.hpp
  Clp_C_Interface.cpp Clp_C_Interface.h
  ClpDualRowDantzig.cpp ClpDualRowDantzig.hpp
  ClpDualRowPivot.cpp ClpDualRowPivot.hpp
  ClpDualRowSteepest.cpp ClpDualRowSteepest.hpp
  ClpDummyMatrix.cpp ClpDummyMatrix.hpp
  ClpDynamicExampleMatrix.cpp ClpDynamicExampleMatrix.hpp
  ClpDynamicMatrix.cpp ClpDynamicMatrix.hpp
  ClpEventHandler.cpp ClpEventHandler.hpp
  ClpFactorization.cpp ClpFactorization.hpp
  ClpGubDynamicMatrix.cpp ClpGubDynamicMatrix.hpp
  ClpGubMatrix.cpp ClpGubMatrix.hpp
  ClpHelperFunctions.cpp ClpHelperFunctions.hpp
  ClpInterior.cpp ClpInterior.hpp
  ClpLinearObjective.cpp ClpLinearObjective.hpp
  ClpMatrixBase.cpp ClpMatrixBase.hpp
  ClpMessage.cpp ClpMessage.hpp
  ClpModel.cpp ClpModel.hpp
  ClpNetworkBasis.cpp ClpNetworkBasis.hpp
  ClpNetworkMatrix.cpp ClpNetworkMatrix.hpp
  ClpNonLinearCost.cpp ClpNonLinearCost.hpp
  ClpNode.cpp ClpNode.hpp
  ClpObjective.cpp ClpObjective.hpp
  ClpPackedMatrix.cpp ClpPackedMatrix.hpp
  ClpParameters.hpp
  ClpPlusMinusOneMatrix.cpp ClpPlusMinusOneMatrix.hpp
  ClpPredictorCorrector.cpp ClpPredictorCorrector.hpp
  ClpPdco.cpp ClpPdco.hpp
  ClpPdcoBase.cpp ClpPdcoBase.hpp
  ClpLsqr.cpp ClpLsqr.hpp
  ClpPresolve.cpp ClpPresolve.hpp
  ClpPrimalColumnDantzig.cpp ClpPrimalColumnDantzig.hpp
  ClpPrimalColumnPivot.cpp ClpPrimalColumnPivot.hpp
  ClpPrimalColumnSteepest.cpp ClpPrimalColumnSteepest.hpp
  ClpQuadraticObjective.cpp ClpQuadraticObjective.hpp
  ClpSimplex.cpp ClpSimplex.hpp
  ClpSimplexDual.cpp ClpSimplexDual.hpp
  ClpSimplexNonlinear.cpp ClpSimplexNonlinear.hpp
  ClpSimplexOther.cpp ClpSimplexOther.hpp
  ClpSimplexPrimal.cpp ClpSimplexPrimal.hpp
  ClpSolve.cpp ClpSolve.hpp
  Idiot.cpp Idiot.hpp
  IdiSolve.cpp)
target_link_libraries(Clp CoinUtils)

add_coin_library(OsiClp Clp/src/OsiClp
  OsiClpSolverInterface.cpp OsiClpSolverInterface.hpp)
target_link_libraries(OsiClp Osi Clp)

include_directories(Osi/src/Osi)

add_coin_library(Osi Osi/src/Osi
  OsiConfig.h
  OsiAuxInfo.cpp OsiAuxInfo.hpp
  OsiBranchingObject.cpp OsiBranchingObject.hpp
  OsiChooseVariable.cpp OsiChooseVariable.hpp
  OsiColCut.cpp OsiColCut.hpp
  OsiCollections.hpp
  OsiCut.cpp OsiCut.hpp
  OsiCuts.cpp OsiCuts.hpp
  OsiNames.cpp
  OsiPresolve.cpp OsiPresolve.hpp
  OsiRowCut.cpp OsiRowCut.hpp
  OsiRowCutDebugger.cpp OsiRowCutDebugger.hpp
  OsiSolverBranch.cpp OsiSolverBranch.hpp
  OsiSolverInterface.cpp OsiSolverInterface.hpp
  OsiSolverParameters.hpp)

add_coin_package(Cgl
  CglConfig.h
  CglCutGenerator.cpp CglCutGenerator.hpp
  CglMessage.cpp CglMessage.hpp
  CglStored.cpp CglStored.hpp
  CglParam.cpp CglParam.hpp
  CglTreeInfo.cpp CglTreeInfo.hpp)

add_coin_library_auto(Cgl/src/CglClique)
add_coin_library_auto(Cgl/src/CglDuplicateRow)
add_coin_library_auto(Cgl/src/CglFlowCover)
add_coin_library_auto(Cgl/src/CglGomory)
add_coin_library_auto(Cgl/src/CglKnapsackCover)
add_coin_library_auto(Cgl/src/CglLandP)
add_coin_library_auto(Cgl/src/CglMixedIntegerRounding2)
add_coin_library_auto(Cgl/src/CglOddHole)
add_coin_library_auto(Cgl/src/CglPreProcess CglClique CglDuplicateRow)
add_coin_library_auto(Cgl/src/CglProbing)
add_coin_library_auto(Cgl/src/CglRedSplit)
add_coin_library_auto(Cgl/src/CglResidualCapacity)
add_coin_library_auto(Cgl/src/CglTwomir)

add_coin_package(Cbc
  Cbc_C_Interface.cpp Cbc_C_Interface.h
  CbcConfig.h
  CbcBranchActual.hpp
  CbcBranchAllDifferent.cpp CbcBranchAllDifferent.hpp
  CbcBranchBase.hpp
  CbcBranchCut.cpp CbcBranchCut.hpp
  CbcBranchDecision.cpp CbcBranchDecision.hpp
  CbcBranchDefaultDecision.cpp CbcBranchDefaultDecision.hpp
  CbcBranchDynamic.cpp CbcBranchDynamic.hpp
  CbcBranchingObject.cpp CbcBranchingObject.hpp
  CbcBranchLotsize.cpp CbcBranchLotsize.hpp
  CbcBranchToFixLots.cpp CbcBranchToFixLots.hpp
  CbcCompareActual.hpp
  CbcCompareBase.hpp
  CbcCompareDefault.cpp CbcCompareDefault.hpp
  CbcCompareDepth.cpp CbcCompareDepth.hpp
  CbcCompareEstimate.cpp CbcCompareEstimate.hpp
  CbcCompareObjective.cpp CbcCompareObjective.hpp
  CbcConsequence.cpp CbcConsequence.hpp
  CbcClique.cpp CbcClique.hpp
  CbcCompare.hpp
  CbcCountRowCut.cpp CbcCountRowCut.hpp
  CbcCutGenerator.cpp CbcCutGenerator.hpp
  CbcCutModifier.cpp CbcCutModifier.hpp
  CbcCutSubsetModifier.cpp CbcCutSubsetModifier.hpp
  CbcDummyBranchingObject.cpp CbcDummyBranchingObject.hpp
  CbcEventHandler.cpp CbcEventHandler.hpp
  CbcFathom.cpp CbcFathom.hpp
  CbcFathomDynamicProgramming.cpp CbcFathomDynamicProgramming.hpp
  CbcFeasibilityBase.hpp
  CbcFixVariable.cpp CbcFixVariable.hpp
  CbcFullNodeInfo.cpp CbcFullNodeInfo.hpp
  CbcFollowOn.cpp CbcFollowOn.hpp
  CbcGeneral.cpp CbcGeneral.hpp
  CbcGeneralDepth.cpp CbcGeneralDepth.hpp
  CbcHeuristic.cpp CbcHeuristic.hpp
  CbcHeuristicDINS.cpp CbcHeuristicDINS.hpp
  CbcHeuristicDive.cpp CbcHeuristicDive.hpp
  CbcHeuristicDiveCoefficient.cpp CbcHeuristicDiveCoefficient.hpp
  CbcHeuristicDiveFractional.cpp CbcHeuristicDiveFractional.hpp
  CbcHeuristicDiveGuided.cpp CbcHeuristicDiveGuided.hpp
  CbcHeuristicDiveLineSearch.cpp CbcHeuristicDiveLineSearch.hpp
  CbcHeuristicDivePseudoCost.cpp CbcHeuristicDivePseudoCost.hpp
  CbcHeuristicDiveVectorLength.cpp CbcHeuristicDiveVectorLength.hpp
  CbcHeuristicFPump.cpp CbcHeuristicFPump.hpp
  CbcHeuristicGreedy.cpp CbcHeuristicGreedy.hpp
  CbcHeuristicLocal.cpp CbcHeuristicLocal.hpp
  CbcHeuristicPivotAndFix.cpp CbcHeuristicPivotAndFix.hpp
  CbcHeuristicRandRound.cpp CbcHeuristicRandRound.hpp
  CbcHeuristicRENS.cpp CbcHeuristicRENS.hpp
  CbcHeuristicRINS.cpp CbcHeuristicRINS.hpp
  CbcHeuristicVND.cpp CbcHeuristicVND.hpp
  CbcMessage.cpp CbcMessage.hpp
  CbcModel.cpp CbcModel.hpp
  CbcNode.cpp CbcNode.hpp
  CbcNodeInfo.cpp CbcNodeInfo.hpp
  CbcNWay.cpp CbcNWay.hpp
  CbcObject.cpp CbcObject.hpp
  CbcObjectUpdateData.cpp CbcObjectUpdateData.hpp
  CbcPartialNodeInfo.cpp CbcPartialNodeInfo.hpp
  CbcSimpleInteger.cpp CbcSimpleInteger.hpp
  CbcSimpleIntegerDynamicPseudoCost.cpp
  CbcSimpleIntegerDynamicPseudoCost.hpp
  CbcSimpleIntegerPseudoCost.cpp
  CbcSimpleIntegerPseudoCost.hpp
  CbcSOS.cpp CbcSOS.hpp
  CbcStatistics.cpp CbcStatistics.hpp
  CbcStrategy.cpp CbcStrategy.hpp
  CbcSubProblem.cpp CbcSubProblem.hpp
  CbcThread.cpp CbcThread.hpp
  CbcTree.cpp CbcTree.hpp
  CbcTreeLocal.cpp CbcTreeLocal.hpp)
target_link_libraries(Cbc
  OsiClp CglFlowCover CglGomory CglKnapsackCover
  CglMixedIntegerRounding2 CglPreProcess CglProbing)

add_coin_library(CbcSolver Cbc/src
  CbcCbcParam.cpp
  Cbc_ampl.cpp Cbc_ampl.h
  CbcLinked.cpp CbcLinked.hpp CbcLinkedUtils.cpp
  unitTestClp.cpp CbcSolver.cpp
  CbcSolverHeuristics.cpp CbcSolverHeuristics.hpp
  CbcSolverAnalyze.cpp CbcSolverAnalyze.hpp
  CbcSolverExpandKnapsack.cpp CbcSolverExpandKnapsack.hpp)
target_link_libraries(CbcSolver
  Cbc Cgl CglLandP CglRedSplit CglResidualCapacity CglTwomir)

add_executable(cbc Cbc/src/CoinSolve.cpp)
target_link_libraries(cbc CbcSolver)

add_definitions(-DHAVE_CONFIG_H)
configure(${CBC_SOURCE_DIR}/CoinUtils/src/config.h)
configure(${CBC_SOURCE_DIR}/Osi/src/Osi/config_osi.h)
configure(${CBC_SOURCE_DIR}/Cbc/src/config.h)
