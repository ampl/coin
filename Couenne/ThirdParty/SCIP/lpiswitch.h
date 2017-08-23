/* Copyright (C) 2012 Stefan Vigerske and others
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * Author: Stefan Vigerske
 */

#ifndef LPISWITCH_H_
#define LPISWITCH_H_

#include "scip/type_retcode.h"

#ifdef __cplusplus
extern "C" {
#endif

enum SCIP_LPISW_LPSolver
{
   SCIP_LPISW_CLP = 0,
   SCIP_LPISW_CPLEX = 1,
   SCIP_LPISW_GUROBI = 2,
   SCIP_LPISW_MOSEK = 3,
   SCIP_LPISW_NONE = 4,
   SCIP_LPISW_QSOPT = 5,
   SCIP_LPISW_SOPLEX = 6,
   SCIP_LPISW_SOPLEX2 = 7,
   SCIP_LPISW_XPRESS = 8,
   SCIP_LPISW_NSOLVERS = 9
};
typedef enum SCIP_LPISW_LPSolver SCIP_LPISW_LPSOLVER;

static const char* const SCIP_LPISW_SOLVERNAMES[SCIP_LPISW_NSOLVERS] =
{ "clp", "cplex", "gurobi", "mosek", "none", "qsopt", "soplex", "soplex2", "xpress" };

/** gets enum of currently used LP solver */
extern
SCIP_LPISW_LPSOLVER SCIPlpiSwitchGetCurrentSolver();

/** sets static function pointers to LPI functions of a particular solver */
extern
SCIP_RETCODE SCIPlpiSwitchSetSolver(
   SCIP_LPISW_LPSOLVER   lpsolver            /**< LP Solver to use */
   );

/** sets static function pointers to LPI functions of default solver (determined during configure) */
extern
SCIP_RETCODE SCIPlpiSwitchSetDefaultSolver();

#ifdef __cplusplus
}
#endif

#endif /* LPISWITCH_H_ */
