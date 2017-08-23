/* Copyright (C) 2012 Stefan Vigerske and others
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * Author: Stefan Vigerske
 */

#include <stdlib.h>

#include "config_scip.h"
#undef SCIP_VERSION

#include "scip/scip.h"
#include "lpiswitch.h"

/*
 *  declaration of LPI functions for all available LP solvers
 */

#ifdef COIN_HAS_CLP
#define LPISW_RENAME_SUFFIX Clp
#include "lpiswitch_rename.h"
#undef __SCIP_LPI_H__
#include "lpi/lpi.h"
#include "lpiswitch_unrename.h"
#endif

#ifdef COIN_HAS_CPX
#undef LPISW_RENAME_SUFFIX
#define LPISW_RENAME_SUFFIX Cplex
#include "lpiswitch_rename.h"
#undef __SCIP_LPI_H__
#include "lpi/lpi.h"
#include "lpiswitch_unrename.h"
#endif

#ifdef COIN_HAS_GRB
#undef LPISW_RENAME_SUFFIX
#define LPISW_RENAME_SUFFIX Gurobi
#include "lpiswitch_rename.h"
#undef __SCIP_LPI_H__
#include "lpi/lpi.h"
#include "lpiswitch_unrename.h"
#endif

#ifdef COIN_HAS_MSK
#undef LPISW_RENAME_SUFFIX
#define LPISW_RENAME_SUFFIX Mosek
#include "lpiswitch_rename.h"
#undef __SCIP_LPI_H__
#include "lpi/lpi.h"
#include "lpiswitch_unrename.h"
#endif

#ifdef COIN_HAS_QSO
#undef LPISW_RENAME_SUFFIX
#define LPISW_RENAME_SUFFIX Qsopt
#include "lpiswitch_rename.h"
#undef __SCIP_LPI_H__
#include "lpi/lpi.h"
#include "lpiswitch_unrename.h"
#endif

#undef LPISW_RENAME_SUFFIX
#define LPISW_RENAME_SUFFIX None
#include "lpiswitch_rename.h"
#undef __SCIP_LPI_H__
#include "lpi/lpi.h"
#include "lpiswitch_unrename.h"

#ifdef COIN_HAS_SOPLEX
#undef LPISW_RENAME_SUFFIX
#define LPISW_RENAME_SUFFIX Soplex
#include "lpiswitch_rename.h"
#undef __SCIP_LPI_H__
#include "lpi/lpi.h"
#include "lpiswitch_unrename.h"

#undef LPISW_RENAME_SUFFIX
#define LPISW_RENAME_SUFFIX Soplex2
#include "lpiswitch_rename.h"
#undef __SCIP_LPI_H__
#include "lpi/lpi.h"
#include "lpiswitch_unrename.h"
#endif

#ifdef COIN_HAS_XPR
#undef LPISW_RENAME_SUFFIX
#define LPISW_RENAME_SUFFIX Xpress
#include "lpiswitch_rename.h"
#undef __SCIP_LPI_H__
#include "lpi/lpi.h"
#include "lpiswitch_unrename.h"
#endif

/*
 *  macros for function pointer declarations
 */

#define LPISW_DECL_LPIGETSOLVERNAME(x) const char* x (void)

#define LPISW_DECL_LPIGETSOLVERDESC(x) const char* x (void)

#define LPISW_DECL_LPIGETSOLVERPOINTER(x) void* x (SCIP_LPI* lpi)

#define LPISW_DECL_LPICREATE(x) SCIP_RETCODE x ( \
   SCIP_LPI**            lpi,         \
   SCIP_MESSAGEHDLR*     messagehdlr, \
   const char*           name,        \
   SCIP_OBJSEN           objsen       \
   )

#define LPISW_DECL_LPIFREE(x) SCIP_RETCODE x (SCIP_LPI** lpi)

#define LPISW_DECL_LPILOADCOLLP(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,       \
   SCIP_OBJSEN           objsen,    \
   int                   ncols,     \
   const SCIP_Real*      obj,       \
   const SCIP_Real*      lb,        \
   const SCIP_Real*      ub,        \
   char**                colnames,  \
   int                   nrows,     \
   const SCIP_Real*      lhs,       \
   const SCIP_Real*      rhs,       \
   char**                rownames,  \
   int                   nnonz,     \
   const int*            beg,       \
   const int*            ind,       \
   const SCIP_Real*      val        \
   )

#define LPISW_DECL_LPIADDCOLS(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,       \
   int                   ncols,     \
   const SCIP_Real*      obj,       \
   const SCIP_Real*      lb,        \
   const SCIP_Real*      ub,        \
   char**                colnames,  \
   int                   nnonz,     \
   const int*            beg,       \
   const int*            ind,       \
   const SCIP_Real*      val        \
   )

#define LPISW_DECL_LPIDELCOLS(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,       \
   int                   firstcol,  \
   int                   lastcol    \
   )

#define LPISW_DECL_LPIDELCOLSET(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,   \
   int*                  dstat  \
   )

#define LPISW_DECL_LPIADDROWS(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,      \
   int                   nrows,    \
   const SCIP_Real*      lhs,      \
   const SCIP_Real*      rhs,      \
   char**                rownames, \
   int                   nnonz,    \
   const int*            beg,      \
   const int*            ind,      \
   const SCIP_Real*      val       \
   )

#define LPISW_DECL_LPIDELROWS(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,       \
   int                   firstrow,  \
   int                   lastrow    \
   )

#define LPISW_DECL_LPIDELROWSET(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,   \
   int*                  dstat  \
   )

#define LPISW_DECL_LPICLEAR(x) SCIP_RETCODE x ( SCIP_LPI* lpi )

#define LPISW_DECL_LPICHGBOUNDS(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,   \
   int                   ncols, \
   const int*            ind,   \
   const SCIP_Real*      lb,    \
   const SCIP_Real*      ub     \
   )

#define LPISW_DECL_LPICHGSIDES(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,    \
   int                   nrows,  \
   const int*            ind,    \
   const SCIP_Real*      lhs,    \
   const SCIP_Real*      rhs     \
   )

#define LPISW_DECL_LPICHGCOEF(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,    \
   int                   row,    \
   int                   col,    \
   SCIP_Real             newval  \
   )

#define LPISW_DECL_LPICHGOBJSEN(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,    \
   SCIP_OBJSEN           objsen  \
   )

#define LPISW_DECL_LPICHGOBJ(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,    \
   int                   ncols,  \
   int*                  ind,    \
   SCIP_Real*            obj     \
   )

#define LPISW_DECL_LPISCALEROW(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,     \
   int                   row,     \
   SCIP_Real             scaleval \
   )

#define LPISW_DECL_LPISCALECOL(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,      \
   int                   col,      \
   SCIP_Real             scaleval  \
   )

#define LPISW_DECL_LPIGETNROWS(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,    \
   int*                  nrows   \
   )

#define LPISW_DECL_LPIGETNCOLS(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,    \
   int*                  ncols   \
   )

#define LPISW_DECL_LPIGETOBJSEN(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,    \
   SCIP_OBJSEN*          objsen  \
   )

#define LPISW_DECL_LPIGETNNONZ(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,   \
   int*                  nnonz  \
   )

#define LPISW_DECL_LPIGETCOLS(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,      \
   int                   firstcol, \
   int                   lastcol,  \
   SCIP_Real*            lb,       \
   SCIP_Real*            ub,       \
   int*                  nnonz,    \
   int*                  beg,      \
   int*                  ind,      \
   SCIP_Real*            val       \
   )

#define LPISW_DECL_LPIGETROWS(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   int                   firstrow,           \
   int                   lastrow,            \
   SCIP_Real*            lhs,                \
   SCIP_Real*            rhs,                \
   int*                  nnonz,              \
   int*                  beg,                \
   int*                  ind,                \
   SCIP_Real*            val                 \
   )

#define LPISW_DECL_LPIGETCOLNAMES(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   int                   firstcol,           \
   int                   lastcol,            \
   char**                colnames,           \
   char*                 namestorage,        \
   int                   namestoragesize,    \
   int*                  storageleft         \
   )

#define LPISW_DECL_LPIGETROWNAMES(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   int                   firstrow,           \
   int                   lastrow,            \
   char**                rownames,           \
   char*                 namestorage,        \
   int                   namestoragesize,    \
   int*                  storageleft         \
   )

#define LPISW_DECL_LPIGETOBJ(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   int                   firstcol,           \
   int                   lastcol,            \
   SCIP_Real*            vals                \
   )

#define LPISW_DECL_LPIGETBOUNDS(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   int                   firstcol,           \
   int                   lastcol,            \
   SCIP_Real*            lbs,                \
   SCIP_Real*            ubs                 \
   )

#define LPISW_DECL_LPIGETSIDES(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   int                   firstrow,           \
   int                   lastrow,            \
   SCIP_Real*            lhss,               \
   SCIP_Real*            rhss                \
   )

#define LPISW_DECL_LPIGETCOEF(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   int                   row,                \
   int                   col,                \
   SCIP_Real*            val                 \
   )

#define LPISW_DECL_LPISOLVEPRIMAL(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPISOLVEDUAL(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPISOLVEBARRIER(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   SCIP_Bool             crossover           \
   )

#define LPISW_DECL_LPISTARTSTRONGBRANCH(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPIENDSTRONGBRANCH(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPISTRONGBRANCHFRAC(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   int                   col,                \
   SCIP_Real             psol,               \
   int                   itlim,              \
   SCIP_Real*            down,               \
   SCIP_Real*            up,                 \
   SCIP_Bool*            downvalid,          \
   SCIP_Bool*            upvalid,            \
   int*                  iter                \
   )

#define LPISW_DECL_LPISTRONGBRANCHESFRAC(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   int*                  cols,               \
   int                   ncols,              \
   SCIP_Real*            psols,              \
   int                   itlim,              \
   SCIP_Real*            down,               \
   SCIP_Real*            up,                 \
   SCIP_Bool*            downvalid,          \
   SCIP_Bool*            upvalid,            \
   int*                  iter                \
   )

#define LPISW_DECL_LPISTRONGBRANCHINT(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   int                   col,                \
   SCIP_Real             psol,               \
   int                   itlim,              \
   SCIP_Real*            down,               \
   SCIP_Real*            up,                 \
   SCIP_Bool*            downvalid,          \
   SCIP_Bool*            upvalid,            \
   int*                  iter                \
   )

#define LPISW_DECL_LPISTRONGBRANCHESINT(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   int*                  cols,               \
   int                   ncols,              \
   SCIP_Real*            psols,              \
   int                   itlim,              \
   SCIP_Real*            down,               \
   SCIP_Real*            up,                 \
   SCIP_Bool*            downvalid,          \
   SCIP_Bool*            upvalid,            \
   int*                  iter                \
   )

#define LPISW_DECL_LPIWASSOLVED(x) SCIP_Bool x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPIGETSOLFEASIBILITY(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   SCIP_Bool*            primalfeasible,     \
   SCIP_Bool*            dualfeasible        \
   )

#define LPISW_DECL_LPIEXISTSPRIMALRAY(x) SCIP_Bool x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPIHASPRIMALRAY(x) SCIP_Bool x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPIISPRIMALUNBOUNDED(x) SCIP_Bool x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPIISPRIMALINFEASIBLE(x) SCIP_Bool x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPIISPRIMALFEASIBLE(x) SCIP_Bool x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPIEXISTSDUALRAY(x) SCIP_Bool x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPIHASDUALRAY(x) SCIP_Bool x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPIISDUALUNBOUNDED(x) SCIP_Bool x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPIISDUALINFEASIBLE(x) SCIP_Bool x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPIISDUALFEASIBLE(x) SCIP_Bool x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPIISOPTIMAL(x) SCIP_Bool x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPIISSTABLE(x) SCIP_Bool x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPIISOBJLIMEXC(x) SCIP_Bool x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPIISITERLIMEXC(x) SCIP_Bool x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPIISTIMELIMEXC(x) SCIP_Bool x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPIGETINTERNALSTATUS(x) int x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPIIGNOREINSTABILITY(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   SCIP_Bool*            success             \
   )

#define LPISW_DECL_LPIGETOBJVAL(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   SCIP_Real*            objval              \
   )

#define LPISW_DECL_LPIGETSOL(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   SCIP_Real*            objval,             \
   SCIP_Real*            primsol,            \
   SCIP_Real*            dualsol,            \
   SCIP_Real*            activity,           \
   SCIP_Real*            redcost             \
   )

#define LPISW_DECL_LPIGETPRIMALRAY(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   SCIP_Real*            ray                 \
   )

#define LPISW_DECL_LPIGETDUALFARKAS(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   SCIP_Real*            dualfarkas          \
   )

#define LPISW_DECL_LPIGETITERATIONS(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   int*                  iterations          \
   )

#define LPISW_DECL_LPIGETREALSOLQUALITY(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   SCIP_LPSOLQUALITY     qualityindicator,   \
   SCIP_Real*            quality             \
   )

#define LPISW_DECL_LPIGETBASE(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   int*                  cstat,              \
   int*                  rstat               \
   )

#define LPISW_DECL_LPISETBASE(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   int*                  cstat,              \
   int*                  rstat               \
   )

#define LPISW_DECL_LPIGETBASISIND(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   int*                  bind                \
   )

#define LPISW_DECL_LPIGETBINVROW(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   int                   r,                  \
   SCIP_Real*            coef                \
   )

#define LPISW_DECL_LPIGETBINVCOL(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   int                   c,                  \
   SCIP_Real*            coef                \
   )

#define LPISW_DECL_LPIGETBINVAROW(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   int                   r,                  \
   const SCIP_Real*      binvrow,            \
   SCIP_Real*            coef                \
   )

#define LPISW_DECL_LPIGETBINVACOL(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   int                   c,                  \
   SCIP_Real*            coef                \
   )

#define LPISW_DECL_LPIGETSTATE(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   BMS_BLKMEM*           blkmem,             \
   SCIP_LPISTATE**       lpistate            \
   )

#define LPISW_DECL_LPISETSTATE(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   BMS_BLKMEM*           blkmem,             \
   SCIP_LPISTATE*        lpistate            \
   )

#define LPISW_DECL_LPICLEARSTATE(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPIFREESTATE(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   BMS_BLKMEM*           blkmem,             \
   SCIP_LPISTATE**       lpistate            \
   )

#define LPISW_DECL_LPIHASSTATEBASIS(x) SCIP_Bool x ( \
   SCIP_LPI*             lpi,                \
   SCIP_LPISTATE*        lpistate            \
   )

#define LPISW_DECL_LPIREADSTATE(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   const char*           fname               \
   )

#define LPISW_DECL_LPIWRITESTATE(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   const char*           fname               \
   )

#define LPISW_DECL_LPIGETNORMS(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   BMS_BLKMEM*           blkmem,             \
   SCIP_LPINORMS**       lpinorms            \
   )

#define LPISW_DECL_LPISETNORMS(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   BMS_BLKMEM*           blkmem,             \
   SCIP_LPINORMS*        lpinorms            \
   )

#define LPISW_DECL_LPIFREENORMS(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   BMS_BLKMEM*           blkmem,             \
   SCIP_LPINORMS**       lpinorms            \
   )

#define LPISW_DECL_LPIGETINTPAR(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   SCIP_LPPARAM          type,               \
   int*                  ival                \
   )

#define LPISW_DECL_LPISETINTPAR(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   SCIP_LPPARAM          type,               \
   int                   ival                \
   )

#define LPISW_DECL_LPIGETREALPAR(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   SCIP_LPPARAM          type,               \
   SCIP_Real*            dval                \
   )

#define LPISW_DECL_LPISETREALPAR(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   SCIP_LPPARAM          type,               \
   SCIP_Real             dval                \
   )

#define LPISW_DECL_LPIINFINITY(x) SCIP_Real x ( \
   SCIP_LPI*             lpi                 \
   )

#define LPISW_DECL_LPIISINFINITY(x) SCIP_Bool x ( \
   SCIP_LPI*             lpi,                \
   SCIP_Real             val                 \
   )

#define LPISW_DECL_LPIREADLP(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   const char*           fname               \
   )

#define LPISW_DECL_LPIWRITELP(x) SCIP_RETCODE x ( \
   SCIP_LPI*             lpi,                \
   const char*           fname               \
   )

/*
 *  pointer to currently used LPI functions
 */

static LPISW_DECL_LPIGETSOLVERNAME((*lpiGetSolverName)) = NULL;
static LPISW_DECL_LPIGETSOLVERDESC((*lpiGetSolverDesc)) = NULL;
static LPISW_DECL_LPIGETSOLVERPOINTER((*lpiGetSolverPointer)) = NULL;
static LPISW_DECL_LPICREATE((*lpiCreate)) = NULL;
static LPISW_DECL_LPIFREE((*lpiFree)) = NULL;
static LPISW_DECL_LPILOADCOLLP((*lpiLoadColLP)) = NULL;
static LPISW_DECL_LPIADDCOLS((*lpiAddCols)) = NULL;
static LPISW_DECL_LPIDELCOLS((*lpiDelCols)) = NULL;
static LPISW_DECL_LPIDELCOLSET((*lpiDelColset)) = NULL;
static LPISW_DECL_LPIADDROWS((*lpiAddRows)) = NULL;
static LPISW_DECL_LPIDELROWS((*lpiDelRows)) = NULL;
static LPISW_DECL_LPIDELROWSET((*lpiDelRowset)) = NULL;
static LPISW_DECL_LPICLEAR((*lpiClear)) = NULL;
static LPISW_DECL_LPICHGBOUNDS((*lpiChgBounds)) = NULL;
static LPISW_DECL_LPICHGSIDES((*lpiChgSides)) = NULL;
static LPISW_DECL_LPICHGCOEF((*lpiChgCoef)) = NULL;
static LPISW_DECL_LPICHGOBJSEN((*lpiChgObjsen)) = NULL;
static LPISW_DECL_LPICHGOBJ((*lpiChgObj)) = NULL;
static LPISW_DECL_LPISCALEROW((*lpiScaleRow)) = NULL;
static LPISW_DECL_LPISCALECOL((*lpiScaleCol)) = NULL;
static LPISW_DECL_LPIGETNROWS((*lpiGetNRows)) = NULL;
static LPISW_DECL_LPIGETNCOLS((*lpiGetNCols)) = NULL;
static LPISW_DECL_LPIGETOBJSEN((*lpiGetObjsen)) = NULL;
static LPISW_DECL_LPIGETNNONZ((*lpiGetNNonz)) = NULL;
static LPISW_DECL_LPIGETCOLS((*lpiGetCols)) = NULL;
static LPISW_DECL_LPIGETROWS((*lpiGetRows)) = NULL;
static LPISW_DECL_LPIGETCOLNAMES((*lpiGetColNames)) = NULL;
static LPISW_DECL_LPIGETROWNAMES((*lpiGetRowNames)) = NULL;
static LPISW_DECL_LPIGETOBJ((*lpiGetObj)) = NULL;
static LPISW_DECL_LPIGETBOUNDS((*lpiGetBounds)) = NULL;
static LPISW_DECL_LPIGETSIDES((*lpiGetSides)) = NULL;
static LPISW_DECL_LPIGETCOEF((*lpiGetCoef)) = NULL;
static LPISW_DECL_LPISOLVEPRIMAL((*lpiSolvePrimal)) = NULL;
static LPISW_DECL_LPISOLVEDUAL((*lpiSolveDual)) = NULL;
static LPISW_DECL_LPISOLVEBARRIER((*lpiSolveBarrier)) = NULL;
static LPISW_DECL_LPISTARTSTRONGBRANCH((*lpiStartStrongbranch)) = NULL;
static LPISW_DECL_LPIENDSTRONGBRANCH((*lpiEndStrongbranch)) = NULL;
static LPISW_DECL_LPISTRONGBRANCHFRAC((*lpiStrongbranchFrac)) = NULL;
static LPISW_DECL_LPISTRONGBRANCHESFRAC((*lpiStrongbranchesFrac)) = NULL;
static LPISW_DECL_LPISTRONGBRANCHINT((*lpiStrongbranchInt)) = NULL;
static LPISW_DECL_LPISTRONGBRANCHESINT((*lpiStrongbranchesInt)) = NULL;
static LPISW_DECL_LPIWASSOLVED((*lpiWasSolved)) = NULL;
static LPISW_DECL_LPIGETSOLFEASIBILITY((*lpiGetSolFeasibility)) = NULL;
static LPISW_DECL_LPIEXISTSPRIMALRAY((*lpiExistsPrimalRay)) = NULL;
static LPISW_DECL_LPIHASPRIMALRAY((*lpiHasPrimalRay)) = NULL;
static LPISW_DECL_LPIISPRIMALUNBOUNDED((*lpiIsPrimalUnbounded)) = NULL;
static LPISW_DECL_LPIISPRIMALINFEASIBLE((*lpiIsPrimalInfeasible)) = NULL;
static LPISW_DECL_LPIISPRIMALFEASIBLE((*lpiIsPrimalFeasible)) = NULL;
static LPISW_DECL_LPIEXISTSDUALRAY((*lpiExistsDualRay)) = NULL;
static LPISW_DECL_LPIHASDUALRAY((*lpiHasDualRay)) = NULL;
static LPISW_DECL_LPIISDUALUNBOUNDED((*lpiIsDualUnbounded)) = NULL;
static LPISW_DECL_LPIISDUALINFEASIBLE((*lpiIsDualInfeasible)) = NULL;
static LPISW_DECL_LPIISDUALFEASIBLE((*lpiIsDualFeasible)) = NULL;
static LPISW_DECL_LPIISOPTIMAL((*lpiIsOptimal)) = NULL;
static LPISW_DECL_LPIISSTABLE((*lpiIsStable)) = NULL;
static LPISW_DECL_LPIISOBJLIMEXC((*lpiIsObjlimExc)) = NULL;
static LPISW_DECL_LPIISITERLIMEXC((*lpiIsIterlimExc)) = NULL;
static LPISW_DECL_LPIISTIMELIMEXC((*lpiIsTimelimExc)) = NULL;
static LPISW_DECL_LPIGETINTERNALSTATUS((*lpiGetInternalStatus)) = NULL;
static LPISW_DECL_LPIIGNOREINSTABILITY((*lpiIgnoreInstability)) = NULL;
static LPISW_DECL_LPIGETOBJVAL((*lpiGetObjval)) = NULL;
static LPISW_DECL_LPIGETSOL((*lpiGetSol)) = NULL;
static LPISW_DECL_LPIGETPRIMALRAY((*lpiGetPrimalRay)) = NULL;
static LPISW_DECL_LPIGETDUALFARKAS((*lpiGetDualfarkas)) = NULL;
static LPISW_DECL_LPIGETITERATIONS((*lpiGetIterations)) = NULL;
static LPISW_DECL_LPIGETREALSOLQUALITY((*lpiGetRealSolQuality)) = NULL;
static LPISW_DECL_LPIGETBASE((*lpiGetBase)) = NULL;
static LPISW_DECL_LPISETBASE((*lpiSetBase)) = NULL;
static LPISW_DECL_LPIGETBASISIND((*lpiGetBasisInd)) = NULL;
static LPISW_DECL_LPIGETBINVROW((*lpiGetBInvRow)) = NULL;
static LPISW_DECL_LPIGETBINVCOL((*lpiGetBInvCol)) = NULL;
static LPISW_DECL_LPIGETBINVAROW((*lpiGetBInvARow)) = NULL;
static LPISW_DECL_LPIGETBINVACOL((*lpiGetBInvACol)) = NULL;
static LPISW_DECL_LPIGETSTATE((*lpiGetState)) = NULL;
static LPISW_DECL_LPISETSTATE((*lpiSetState)) = NULL;
static LPISW_DECL_LPICLEARSTATE((*lpiClearState)) = NULL;
static LPISW_DECL_LPIFREESTATE((*lpiFreeState)) = NULL;
static LPISW_DECL_LPIHASSTATEBASIS((*lpiHasStateBasis)) = NULL;
static LPISW_DECL_LPIREADSTATE((*lpiReadState)) = NULL;
static LPISW_DECL_LPIWRITESTATE((*lpiWriteState)) = NULL;
static LPISW_DECL_LPIGETNORMS((*lpiGetNorms)) = NULL;
static LPISW_DECL_LPISETNORMS((*lpiSetNorms)) = NULL;
static LPISW_DECL_LPIFREENORMS((*lpiFreeNorms)) = NULL;
static LPISW_DECL_LPIGETINTPAR((*lpiGetIntpar)) = NULL;
static LPISW_DECL_LPISETINTPAR((*lpiSetIntpar)) = NULL;
static LPISW_DECL_LPIGETREALPAR((*lpiGetRealpar)) = NULL;
static LPISW_DECL_LPISETREALPAR((*lpiSetRealpar)) = NULL;
static LPISW_DECL_LPIINFINITY((*lpiInfinity)) = NULL;
static LPISW_DECL_LPIISINFINITY((*lpiIsInfinity)) = NULL;
static LPISW_DECL_LPIREADLP((*lpiReadLP)) = NULL;
static LPISW_DECL_LPIWRITELP((*lpiWriteLP)) = NULL;

static SCIP_LPISW_LPSOLVER currentsolver = SCIP_LPISW_NSOLVERS;

/** macro to set static function pointers to LPI functions of a particular solver
 */
#define LPISW_ASSIGN_FUNCPTRS(x) \
   lpiGetSolverName = SCIPlpiGetSolverName ## x ; \
   lpiGetSolverDesc = SCIPlpiGetSolverDesc ## x ; \
   lpiGetSolverPointer = SCIPlpiGetSolverPointer ## x ; \
   lpiCreate = SCIPlpiCreate ## x ; \
   lpiFree = SCIPlpiFree ## x ; \
   lpiLoadColLP = SCIPlpiLoadColLP ## x ; \
   lpiAddCols = SCIPlpiAddCols ## x ; \
   lpiDelCols = SCIPlpiDelCols ## x ; \
   lpiDelColset = SCIPlpiDelColset ## x ; \
   lpiAddRows = SCIPlpiAddRows ## x ; \
   lpiDelRows = SCIPlpiDelRows ## x ; \
   lpiDelRowset = SCIPlpiDelRowset ## x ; \
   lpiClear = SCIPlpiClear ## x ; \
   lpiChgBounds = SCIPlpiChgBounds ## x ; \
   lpiChgSides = SCIPlpiChgSides ## x ; \
   lpiChgCoef = SCIPlpiChgCoef ## x ; \
   lpiChgObjsen = SCIPlpiChgObjsen ## x ; \
   lpiChgObj = SCIPlpiChgObj ## x ; \
   lpiScaleRow = SCIPlpiScaleRow ## x ; \
   lpiScaleCol = SCIPlpiScaleCol ## x ; \
   lpiGetNRows = SCIPlpiGetNRows ## x ; \
   lpiGetNCols = SCIPlpiGetNCols ## x ; \
   lpiGetObjsen = SCIPlpiGetObjsen ## x ; \
   lpiGetNNonz = SCIPlpiGetNNonz ## x ; \
   lpiGetCols = SCIPlpiGetCols ## x ; \
   lpiGetRows = SCIPlpiGetRows ## x ; \
   lpiGetColNames = SCIPlpiGetColNames ## x ; \
   lpiGetRowNames = SCIPlpiGetRowNames ## x ; \
   lpiGetObj = SCIPlpiGetObj ## x ; \
   lpiGetBounds = SCIPlpiGetBounds ## x ; \
   lpiGetSides = SCIPlpiGetSides ## x ; \
   lpiGetCoef = SCIPlpiGetCoef ## x ; \
   lpiSolvePrimal = SCIPlpiSolvePrimal ## x ; \
   lpiSolveDual = SCIPlpiSolveDual ## x ; \
   lpiSolveBarrier = SCIPlpiSolveBarrier ## x ; \
   lpiStartStrongbranch = SCIPlpiStartStrongbranch ## x ; \
   lpiEndStrongbranch = SCIPlpiEndStrongbranch ## x ; \
   lpiStrongbranchFrac = SCIPlpiStrongbranchFrac ## x ; \
   lpiStrongbranchesFrac = SCIPlpiStrongbranchesFrac ## x ; \
   lpiStrongbranchInt = SCIPlpiStrongbranchInt ## x ; \
   lpiStrongbranchesInt = SCIPlpiStrongbranchesInt ## x ; \
   lpiWasSolved = SCIPlpiWasSolved ## x ; \
   lpiGetSolFeasibility = SCIPlpiGetSolFeasibility ## x ; \
   lpiExistsPrimalRay = SCIPlpiExistsPrimalRay ## x ; \
   lpiHasPrimalRay = SCIPlpiHasPrimalRay ## x ; \
   lpiIsPrimalUnbounded = SCIPlpiIsPrimalUnbounded ## x ; \
   lpiIsPrimalInfeasible = SCIPlpiIsPrimalInfeasible ## x ; \
   lpiIsPrimalFeasible = SCIPlpiIsPrimalFeasible ## x ; \
   lpiExistsDualRay = SCIPlpiExistsDualRay ## x ; \
   lpiHasDualRay = SCIPlpiHasDualRay ## x ; \
   lpiIsDualUnbounded = SCIPlpiIsDualUnbounded ## x ; \
   lpiIsDualInfeasible = SCIPlpiIsDualInfeasible ## x ; \
   lpiIsDualFeasible = SCIPlpiIsDualFeasible ## x ; \
   lpiIsOptimal = SCIPlpiIsOptimal ## x ; \
   lpiIsStable = SCIPlpiIsStable ## x ; \
   lpiIsObjlimExc = SCIPlpiIsObjlimExc ## x ; \
   lpiIsIterlimExc = SCIPlpiIsIterlimExc ## x ; \
   lpiIsTimelimExc = SCIPlpiIsTimelimExc ## x ; \
   lpiGetInternalStatus = SCIPlpiGetInternalStatus ## x ; \
   lpiIgnoreInstability = SCIPlpiIgnoreInstability ## x ; \
   lpiGetObjval = SCIPlpiGetObjval ## x ; \
   lpiGetSol = SCIPlpiGetSol ## x ; \
   lpiGetPrimalRay = SCIPlpiGetPrimalRay ## x ; \
   lpiGetDualfarkas = SCIPlpiGetDualfarkas ## x ; \
   lpiGetIterations = SCIPlpiGetIterations ## x ; \
   lpiGetRealSolQuality = SCIPlpiGetRealSolQuality ## x ; \
   lpiGetBase = SCIPlpiGetBase ## x ; \
   lpiSetBase = SCIPlpiSetBase ## x ; \
   lpiGetBasisInd = SCIPlpiGetBasisInd ## x ; \
   lpiGetBInvRow = SCIPlpiGetBInvRow ## x ; \
   lpiGetBInvCol = SCIPlpiGetBInvCol ## x ; \
   lpiGetBasisInd = SCIPlpiGetBasisInd ## x ; \
   lpiGetBInvARow = SCIPlpiGetBInvARow ## x ; \
   lpiGetBInvACol = SCIPlpiGetBInvACol ## x ; \
   lpiGetState = SCIPlpiGetState ## x ; \
   lpiSetState = SCIPlpiSetState ## x ; \
   lpiClearState = SCIPlpiClearState ## x ; \
   lpiFreeState = SCIPlpiFreeState ## x ; \
   lpiHasStateBasis = SCIPlpiHasStateBasis ## x ; \
   lpiReadState = SCIPlpiReadState ## x ; \
   lpiWriteState = SCIPlpiWriteState ## x ; \
   lpiGetNorms = SCIPlpiGetNorms ## x ; \
   lpiSetNorms = SCIPlpiSetNorms ## x ; \
   lpiFreeNorms = SCIPlpiFreeNorms ## x ; \
   lpiGetIntpar = SCIPlpiGetIntpar ## x ; \
   lpiSetIntpar = SCIPlpiSetIntpar ## x ; \
   lpiGetRealpar = SCIPlpiGetRealpar ## x ; \
   lpiSetRealpar = SCIPlpiSetRealpar ## x ; \
   lpiInfinity = SCIPlpiInfinity ## x ; \
   lpiIsInfinity = SCIPlpiIsInfinity ## x ; \
   lpiReadLP = SCIPlpiReadLP ## x ; \
   lpiWriteLP = SCIPlpiWriteLP ## x ;

/** gets enum of currently used LP solver */
SCIP_LPISW_LPSOLVER SCIPlpiSwitchGetCurrentSolver()
{
   return currentsolver;
}

/** sets static function pointers to LPI functions of a particular solver */
SCIP_RETCODE SCIPlpiSwitchSetSolver(
   SCIP_LPISW_LPSOLVER   lpsolver            /**< LP Solver to use */
   )
{
   currentsolver = lpsolver;
   switch( lpsolver )
   {
      case SCIP_LPISW_CLP:
      {
#ifdef COIN_HAS_CLP
         LPISW_ASSIGN_FUNCPTRS(Clp);
         return SCIP_OKAY;
#else
         SCIPerrorMessage("CLP not available as LP solver.\n");
         return SCIP_PLUGINNOTFOUND;
#endif
      }

      case SCIP_LPISW_CPLEX:
      {
#ifdef COIN_HAS_CPX
         LPISW_ASSIGN_FUNCPTRS(Cplex);
         return SCIP_OKAY;
#else
         SCIPerrorMessage("Cplex not available as LP solver.\n");
         return SCIP_PLUGINNOTFOUND;
#endif
      }

      case SCIP_LPISW_GUROBI:
      {
#ifdef COIN_HAS_GRB
         LPISW_ASSIGN_FUNCPTRS(Gurobi);
         return SCIP_OKAY;
#else
         SCIPerrorMessage("Gurobi not available as LP solver.\n");
         return SCIP_PLUGINNOTFOUND;
#endif
      }

      case SCIP_LPISW_MOSEK:
      {
#ifdef COIN_HAS_MSK
         LPISW_ASSIGN_FUNCPTRS(Mosek);
         return SCIP_OKAY;
#else
         SCIPerrorMessage("Mosek not available as LP solver.\n");
         return SCIP_PLUGINNOTFOUND;
#endif
      }

      case SCIP_LPISW_NONE:
      {
         LPISW_ASSIGN_FUNCPTRS(None);
         return SCIP_OKAY;
      }

      case SCIP_LPISW_QSOPT:
      {
#ifdef COIN_HAS_QSO
         LPISW_ASSIGN_FUNCPTRS(Qsopt);
         return SCIP_OKAY;
#else
         SCIPerrorMessage("Qsopt not available as LP solver.\n");
         return SCIP_PLUGINNOTFOUND;
#endif
      }

      case SCIP_LPISW_SOPLEX:
      {
#ifdef COIN_HAS_SOPLEX
         LPISW_ASSIGN_FUNCPTRS(Soplex);
         return SCIP_OKAY;
#else
         SCIPerrorMessage("SoPlex not available as LP solver.\n");
         return SCIP_PLUGINNOTFOUND;
#endif
      }

      case SCIP_LPISW_SOPLEX2:
      {
#ifdef COIN_HAS_SOPLEX
         LPISW_ASSIGN_FUNCPTRS(Soplex2);
         return SCIP_OKAY;
#else
         SCIPerrorMessage("SoPlex not available as LP solver.\n");
         return SCIP_PLUGINNOTFOUND;
#endif
      }

      case SCIP_LPISW_XPRESS:
      {
#ifdef COIN_HAS_XPR
         LPISW_ASSIGN_FUNCPTRS(Xpress);
         return SCIP_OKAY;
#else
         SCIPerrorMessage("Xpress not available as LP solver.\n");
         return SCIP_PLUGINNOTFOUND;
#endif
      }
      default: ;
   }
   SCIPerrorMessage("Unknown LP solver id %d\n", lpsolver);

   return SCIP_ERROR;
}

/** sets static function pointers to LPI functions of default solver (determined during configure) */
SCIP_RETCODE SCIPlpiSwitchSetDefaultSolver()
{
   SCIP_CALL( SCIPlpiSwitchSetSolver(SCIP_LPISW_LPDEFAULT) );

   return SCIP_OKAY;
}

/*
 * our implementation of the SCIP LPI
 * calls whatever has been set in our function pointers
 */

#undef __SCIP_LPI_H__
#include "lpi/lpi.h"

/** gets name and version of LP solver */
const char* SCIPlpiGetSolverName(
   void
   )
{
   return (*lpiGetSolverName)();
}

/** gets description of LP solver (developer, webpage, ...) */
const char* SCIPlpiGetSolverDesc(
   void
   )
{
   return (*lpiGetSolverDesc)();
}

/** gets pointer for LP solver - use only with great care
 *
 *  The behavior of this function depends on the solver and its use is
 *  therefore only recommended if you really know what you are
 *  doing. In general, it returns a pointer to the LP solver object.
 */
void* SCIPlpiGetSolverPointer(
   SCIP_LPI*             lpi                 /**< pointer to an LP interface structure */
   )
{
   return (*lpiGetSolverPointer)(lpi);
}

/** creates an LP problem object */
SCIP_RETCODE SCIPlpiCreate(
   SCIP_LPI**            lpi,                /**< pointer to an LP interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   )
{
   return (*lpiCreate)(lpi, messagehdlr, name, objsen);
}

/** deletes an LP problem object */
SCIP_RETCODE SCIPlpiFree(
   SCIP_LPI**            lpi                 /**< pointer to an LP interface structure */
   )
{
   return (*lpiFree)(lpi);
}

/** copies LP data with column matrix into LP solver */
SCIP_RETCODE SCIPlpiLoadColLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen,             /**< objective sense */
   int                   ncols,              /**< number of columns */
   const SCIP_Real*      obj,                /**< objective function values of columns */
   const SCIP_Real*      lb,                 /**< lower bounds of columns */
   const SCIP_Real*      ub,                 /**< upper bounds of columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nrows,              /**< number of rows */
   const SCIP_Real*      lhs,                /**< left hand sides of rows */
   const SCIP_Real*      rhs,                /**< right hand sides of rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array */
   const int*            ind,                /**< row indices of constraint matrix entries */
   const SCIP_Real*      val                 /**< values of constraint matrix entries */
   )
{
   return (*lpiLoadColLP)(lpi, objsen, ncols, obj, lb, ub, colnames, nrows, lhs, rhs, rownames, nnonz, beg, ind, val);
}

/** adds columns to the LP */
SCIP_RETCODE SCIPlpiAddCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to be added */
   const SCIP_Real*      obj,                /**< objective function values of new columns */
   const SCIP_Real*      lb,                 /**< lower bounds of new columns */
   const SCIP_Real*      ub,                 /**< upper bounds of new columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   return (*lpiAddCols)(lpi, ncols, obj, lb, ub, colnames, nnonz, beg, ind, val);
}

/** deletes all columns in the given range from LP */
SCIP_RETCODE SCIPlpiDelCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   )
{
   return (*lpiDelCols)(lpi, firstcol, lastcol);
}

/** deletes columns from SCIP_LPI; the new position of a column must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelColset(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of columns
                                              *   input:  1 if column should be deleted, 0 if not
                                              *   output: new position of column, -1 if column was deleted */
   )
{
   return (*lpiDelColset)(lpi, dstat);
}

/** adds rows to the LP */
SCIP_RETCODE SCIPlpiAddRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to be added */
   const SCIP_Real*      lhs,                /**< left hand sides of new rows */
   const SCIP_Real*      rhs,                /**< right hand sides of new rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   return (*lpiAddRows)(lpi, nrows, lhs, rhs, rownames, nnonz, beg, ind, val);
}

/** deletes all rows in the given range from LP */
SCIP_RETCODE SCIPlpiDelRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{
   return (*lpiDelRows)(lpi, firstrow, lastrow);
}

/** deletes rows from SCIP_LPI; the new position of a row must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelRowset(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   )
{
   return (*lpiDelRowset)(lpi, dstat);
}

/** clears the whole LP */
SCIP_RETCODE SCIPlpiClear(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiClear)(lpi);
}

/** changes lower and upper bounds of columns */
SCIP_RETCODE SCIPlpiChgBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   const int*            ind,                /**< column indices */
   const SCIP_Real*      lb,                 /**< values for the new lower bounds */
   const SCIP_Real*      ub                  /**< values for the new upper bounds */
   )
{
   return (*lpiChgBounds)(lpi, ncols, ind, lb, ub);
}

/** changes left and right hand sides of rows */
SCIP_RETCODE SCIPlpiChgSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to change sides for */
   const int*            ind,                /**< row indices */
   const SCIP_Real*      lhs,                /**< new values for left hand sides */
   const SCIP_Real*      rhs                 /**< new values for right hand sides */
   )
{
   return (*lpiChgSides)(lpi, nrows, ind, lhs, rhs);
}

/** changes a single coefficient */
SCIP_RETCODE SCIPlpiChgCoef(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient to change */
   int                   col,                /**< column number of coefficient to change */
   SCIP_Real             newval              /**< new value of coefficient */
   )
{
   return (*lpiChgCoef)(lpi, row, col, newval);
}

/** changes the objective sense */
SCIP_RETCODE SCIPlpiChgObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{
   return (*lpiChgObjsen)(lpi, objsen);
}

/** changes objective values of columns in the LP */
SCIP_RETCODE SCIPlpiChgObj(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change objective value for */
   int*                  ind,                /**< column indices to change objective value for */
   SCIP_Real*            obj                 /**< new objective values for columns */
   )
{
   return (*lpiChgObj)(lpi, ncols, ind, obj);
}

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
SCIP_RETCODE SCIPlpiScaleRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   return (*lpiScaleRow)(lpi, row, scaleval);
}

/** multiplies a column with a non-zero scalar; the objective value is multiplied with the scalar, and the bounds
 *  are divided by the scalar; for negative scalars, the column's bounds are switched
 */
SCIP_RETCODE SCIPlpiScaleCol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   return (*lpiScaleCol)(lpi, col, scaleval);
}

/** gets the number of rows in the LP */
SCIP_RETCODE SCIPlpiGetNRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nrows               /**< pointer to store the number of rows */
   )
{
   return (*lpiGetNRows)(lpi, nrows);
}

/** gets the number of columns in the LP */
SCIP_RETCODE SCIPlpiGetNCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{
   return (*lpiGetNCols)(lpi, ncols);
}

/** gets the objective sense of the LP */
SCIP_RETCODE SCIPlpiGetObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN*          objsen              /**< pointer to store objective sense */
   )
{
   return (*lpiGetObjsen)(lpi, objsen);
}

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiGetNNonz(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{
   return (*lpiGetNNonz)(lpi, nnonz);
}

/** gets columns from LP problem object; the arrays have to be large enough to store all values;
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiGetCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get from LP */
   int                   lastcol,            /**< last column to get from LP */
   SCIP_Real*            lb,                 /**< buffer to store the lower bound vector, or NULL */
   SCIP_Real*            ub,                 /**< buffer to store the upper bound vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each column in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   return (*lpiGetCols)(lpi, firstcol, lastcol, lb, ub, nnonz, beg, ind, val);
}

/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiGetRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get from LP */
   int                   lastrow,            /**< last row to get from LP */
   SCIP_Real*            lhs,                /**< buffer to store left hand side vector, or NULL */
   SCIP_Real*            rhs,                /**< buffer to store right hand side vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each row in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store row indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   return (*lpiGetRows)(lpi, firstrow, lastrow, lhs, rhs, nnonz, beg, ind, val);
}

/** gets column names */
SCIP_RETCODE SCIPlpiGetColNames(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get name from LP */
   int                   lastcol,            /**< last column to get name from LP */
   char**                colnames,           /**< pointers to column names (of size at least lastcol-firstcol+1) */
   char*                 namestorage,        /**< storage for col names */
   int                   namestoragesize,    /**< size of namestorage (if 0, -storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) */
   )
{
   return (*lpiGetColNames)(lpi, firstcol, lastcol, colnames, namestorage, namestoragesize, storageleft);
}

/** gets row names */
SCIP_RETCODE SCIPlpiGetRowNames(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get name from LP */
   int                   lastrow,            /**< last row to get name from LP */
   char**                rownames,           /**< pointers to row names (of size at least lastrow-firstrow+1) */
   char*                 namestorage,        /**< storage for row names */
   int                   namestoragesize,    /**< size of namestorage (if 0, -storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) */
   )
{
   return (*lpiGetRowNames)(lpi, firstrow, lastrow, rownames, namestorage, namestoragesize, storageleft);
}

/** gets objective coefficients from LP problem object */
SCIP_RETCODE SCIPlpiGetObj(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   SCIP_Real*            vals                /**< array to store objective coefficients */
   )
{
   return (*lpiGetObj)(lpi, firstcol, lastcol, vals);
}

/** gets current bounds from LP problem object */
SCIP_RETCODE SCIPlpiGetBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get bounds for */
   int                   lastcol,            /**< last column to get bounds for */
   SCIP_Real*            lbs,                /**< array to store lower bound values, or NULL */
   SCIP_Real*            ubs                 /**< array to store upper bound values, or NULL */
   )
{
   return (*lpiGetBounds)(lpi, firstcol, lastcol, lbs, ubs);
}

/** gets current row sides from LP problem object */
SCIP_RETCODE SCIPlpiGetSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            lhss,               /**< array to store left hand side values, or NULL */
   SCIP_Real*            rhss                /**< array to store right hand side values, or NULL */
   )
{
   return (*lpiGetSides)(lpi, firstrow, lastrow, lhss, rhss);
}

/** gets a single coefficient */
SCIP_RETCODE SCIPlpiGetCoef(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
   return (*lpiGetCoef)(lpi, row, col, val);
}

/** calls primal simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolvePrimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiSolvePrimal)(lpi);
}

/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolveDual(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiSolveDual)(lpi);
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiSolveBarrier(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{
   return (*lpiSolveBarrier)(lpi, crossover);
}

/** start strong branching - call before any strong branching */
SCIP_RETCODE SCIPlpiStartStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiStartStrongbranch)(lpi);
}

/** end strong branching - call after any strong branching */
SCIP_RETCODE SCIPlpiEndStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiEndStrongbranch)(lpi);
}

/** performs strong branching iterations on one @b fractional candidate */
SCIP_RETCODE SCIPlpiStrongbranchFrac(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< fractional current primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   return (*lpiStrongbranchFrac)(lpi, col, psol, itlim, down, up, downvalid, upvalid, iter);
}

/** performs strong branching iterations on given @b fractional candidates */
SCIP_RETCODE SCIPlpiStrongbranchesFrac(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cols,               /**< columns to apply strong branching on */
   int                   ncols,              /**< number of columns */
   SCIP_Real*            psols,              /**< fractional current primal solution values of columns */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   return (*lpiStrongbranchesFrac)(lpi, cols, ncols, psols, itlim, down, up, downvalid, upvalid, iter);
}

/** performs strong branching iterations on one candidate with @b integral value */
SCIP_RETCODE SCIPlpiStrongbranchInt(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< current integral primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   return (*lpiStrongbranchInt)(lpi, col, psol, itlim, down, up, downvalid, upvalid, iter);
}

/** performs strong branching iterations on given candidates with @b integral values */
SCIP_RETCODE SCIPlpiStrongbranchesInt(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cols,               /**< columns to apply strong branching on */
   int                   ncols,              /**< number of columns */
   SCIP_Real*            psols,              /**< current integral primal solution values of columns */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   return (*lpiStrongbranchesInt)(lpi, cols, ncols, psols, itlim, down, up, downvalid, upvalid, iter);
}

/** returns whether a solve method was called after the last modification of the LP */
SCIP_Bool SCIPlpiWasSolved(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiWasSolved)(lpi);
}

/** gets information about primal and dual feasibility of the current LP solution */
SCIP_RETCODE SCIPlpiGetSolFeasibility(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   return (*lpiGetSolFeasibility)(lpi, primalfeasible, dualfeasible);
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiExistsPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiExistsPrimalRay)(lpi);
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiHasPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiHasPrimalRay)(lpi);
}

/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_Bool SCIPlpiIsPrimalUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiIsPrimalUnbounded)(lpi);
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiIsPrimalInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiIsPrimalInfeasible)(lpi);
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiIsPrimalFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiIsPrimalFeasible)(lpi);
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExistsDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiExistsDualRay)(lpi);
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiHasDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiHasDualRay)(lpi);
}

/** returns TRUE iff LP is proven to be dual unbounded */
SCIP_Bool SCIPlpiIsDualUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiIsDualUnbounded)(lpi);
}

/** returns TRUE iff LP is proven to be dual infeasible */
SCIP_Bool SCIPlpiIsDualInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiIsDualInfeasible)(lpi);
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiIsDualFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiIsDualFeasible)(lpi);
}

/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiIsOptimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiIsOptimal)(lpi);
}

/** returns TRUE iff current LP basis is stable */
SCIP_Bool SCIPlpiIsStable(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiIsStable)(lpi);
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiIsObjlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiIsObjlimExc)(lpi);
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiIsIterlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiIsIterlimExc)(lpi);
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiIsTimelimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiIsTimelimExc)(lpi);
}

/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiGetInternalStatus)(lpi);
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiIgnoreInstability(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   return (*lpiIgnoreInstability)(lpi, success);
}

/** gets objective value of solution */
SCIP_RETCODE SCIPlpiGetObjval(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   return (*lpiGetObjval)(lpi, objval);
}

/** gets primal and dual solution vectors for feasible LPs */
SCIP_RETCODE SCIPlpiGetSol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_Real*            primsol,            /**< primal solution vector, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< dual solution vector, may be NULL if not needed */
   SCIP_Real*            activity,           /**< row activity vector, may be NULL if not needed */
   SCIP_Real*            redcost             /**< reduced cost vector, may be NULL if not needed */
   )
{
   return (*lpiGetSol)(lpi, objval, primsol, dualsol, activity, redcost);
}

/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiGetPrimalRay(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            ray                 /**< primal ray */
   )
{
   return (*lpiGetPrimalRay)(lpi, ray);
}

/** gets dual Farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiGetDualfarkas(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            dualfarkas          /**< dual Farkas row multipliers */
   )
{
   return (*lpiGetDualfarkas)(lpi, dualfarkas);
}

/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiGetIterations(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   return (*lpiGetIterations)(lpi, iterations);
}

/** gets information about the quality of an LP solution
 *
 *  Such information is usually only available, if also a (maybe not optimal) solution is available.
 *  The LPI should return SCIP_INVALID for *quality, if the requested quantity is not available.
 */
SCIP_RETCODE SCIPlpiGetRealSolQuality(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPSOLQUALITY     qualityindicator,   /**< indicates which quality should be returned */
   SCIP_Real*            quality             /**< pointer to store quality number */
   )
{
   return (*lpiGetRealSolQuality)(lpi, qualityindicator, quality);
}

/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
SCIP_RETCODE SCIPlpiGetBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   )
{
   return (*lpiGetBase)(lpi, cstat, rstat);
}

/** sets current basis status for columns and rows */
SCIP_RETCODE SCIPlpiSetBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array with column basis status */
   int*                  rstat               /**< array with row basis status */
   )
{
   return (*lpiSetBase)(lpi, cstat, rstat);
}

/** returns the indices of the basic columns and rows */
SCIP_RETCODE SCIPlpiGetBasisInd(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  bind                /**< basic column n gives value n, basic row m gives value -1-m */
   )
{
   return (*lpiGetBasisInd)(lpi, bind);
}

/** get dense row of inverse basis matrix B^-1 */
SCIP_RETCODE SCIPlpiGetBInvRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the row */
   )
{
   return (*lpiGetBInvRow)(lpi, r, coef);
}

/** get dense column of inverse basis matrix B^-1 */
SCIP_RETCODE SCIPlpiGetBInvCol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the LP;
                                              *   you have to call SCIPlpiGetBasisInd() to get the array which links the
                                              *   B^-1 column numbers to the row and column numbers of the LP!
                                              *   c must be between 0 and nrows-1, since the basis has the size
                                              *   nrows * nrows */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the column */
   )
{
   return (*lpiGetBInvCol)(lpi, c, coef);
}

/** get dense row of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiGetBInvARow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   const SCIP_Real*      binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiGetBInvRow(), or NULL */
   SCIP_Real*            coef                /**< vector to return coefficients */
   )
{
   return (*lpiGetBInvARow)(lpi, r, binvrow, coef);
}

/** get dense column of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiGetBInvACol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number */
   SCIP_Real*            coef                /**< vector to return coefficients */
   )
{
   return (*lpiGetBInvACol)(lpi, c, coef);
}

/** stores LPi state (like basis information) into lpistate object */
SCIP_RETCODE SCIPlpiGetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{
   return (*lpiGetState)(lpi, blkmem, lpistate);
}

/** loads LPi state (like basis information) into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiGetState()
 */
SCIP_RETCODE SCIPlpiSetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE*        lpistate            /**< LPi state information (like basis information) */
   )
{
   return (*lpiSetState)(lpi, blkmem, lpistate);
}

/** clears current LPi state (like basis information) of the solver */
SCIP_RETCODE SCIPlpiClearState(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiClearState)(lpi);
}

/** frees LPi state information */
SCIP_RETCODE SCIPlpiFreeState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{
   return (*lpiFreeState)(lpi, blkmem, lpistate);
}

/** checks, whether the given LPi state contains simplex basis information */
SCIP_Bool SCIPlpiHasStateBasis(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LPi state information (like basis information) */
   )
{
   return (*lpiHasStateBasis)(lpi, lpistate);
}

/** reads LPi state (like basis information from a file */
SCIP_RETCODE SCIPlpiReadState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   return (*lpiReadState)(lpi, fname);
}

/** writes LPi state (like basis information) to a file */
SCIP_RETCODE SCIPlpiWriteState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   return (*lpiWriteState)(lpi, fname);
}

/** stores LPi pricing norms into lpinorms object */
SCIP_RETCODE SCIPlpiGetNorms(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LPi pricing norms information */
   )
{
   return (*lpiGetNorms)(lpi, blkmem, lpinorms);
}

/** loads LPi pricing norms into solver; note that the LP might have been extended with additional
 *  columns and rows since the norms were stored with SCIPlpiGetNorms()
 */
SCIP_RETCODE SCIPlpiSetNorms(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS*        lpinorms            /**< LPi pricing norms information */
   )
{
   return (*lpiSetNorms)(lpi, blkmem, lpinorms);
}

/** frees LPi pricing norms information */
SCIP_RETCODE SCIPlpiFreeNorms(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LPi pricing norms information */
   )
{
   return (*lpiFreeNorms)(lpi, blkmem, lpinorms);
}

/** gets integer parameter of LP */
SCIP_RETCODE SCIPlpiGetIntpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int*                  ival                /**< buffer to store the parameter value */
   )
{
   return (*lpiGetIntpar)(lpi, type, ival);
}

/** sets integer parameter of LP */
SCIP_RETCODE SCIPlpiSetIntpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   return (*lpiSetIntpar)(lpi, type, ival);
}

/** gets floating point parameter of LP */
SCIP_RETCODE SCIPlpiGetRealpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   )
{
   return (*lpiGetRealpar)(lpi, type, dval);
}

/** sets floating point parameter of LP */
SCIP_RETCODE SCIPlpiSetRealpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   return (*lpiSetRealpar)(lpi, type, dval);
}

/** returns value treated as infinity in the LP solver */
SCIP_Real SCIPlpiInfinity(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return (*lpiInfinity)(lpi);
}

/** checks if given value is treated as infinity in the LP solver */
SCIP_Bool SCIPlpiIsInfinity(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real             val                 /**< value to be checked for infinity */
   )
{
   return (*lpiIsInfinity)(lpi, val);
}

/** reads LP from a file */
SCIP_RETCODE SCIPlpiReadLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   return (*lpiReadLP)(lpi, fname);
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiWriteLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   return (*lpiWriteLP)(lpi, fname);
}
