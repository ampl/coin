/* $Id: cons_rowcuts.cpp 870 2012-07-16 21:07:11Z stefan $
 *
 * @file   cons_rowcuts.c
 * @brief  constraint handler for rowcuts constraints
 *         enables separation of convexification cuts during SCIP solution procedure
 * @author Pietro Belotti
 * @author Timo Berthold
 * @license This file is licensed under the Eclipse Public License (EPL)
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

/**@file   cons_rowcuts.c
 * @ingroup CONSHDLRS 
 * @brief  constraint handler for rowcuts constraints
 * @author Pietro Belotti
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "CouenneConfig.h"
#ifdef COIN_HAS_SCIP

#include "cons_rowcuts.h"
#include "scip/cons_linear.h"
#include "scip/scip.h"

#include "CouenneProblem.hpp"

/* constraint handler properties */
#define CONSHDLR_NAME          "rowcuts"
#define CONSHDLR_DESC          "adds row cuts generated by Couenne in a Branch-and-Check fashion"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -9999999 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -9999999 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS         0 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

#define DEFAULT_MAXCUTTINGROUNDS     5 /**< how many rounds of cuts should be applied at most? */

using namespace Couenne;

/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   CouenneCutGenerator*  cutgenerator;       /* CouenneCutGenerator for linearization cuts */
   OsiSolverInterface*   milp;               /* Couenne's MILP relaxation of Couenne's MINLP */
   int                   maxcuttingrounds;   /* how many rounds of cuts should be applied at most */
   int                   ncuttingrounds;     /* how many rounds of cuts have been applied already */
};


/*
 * Local methods
 */

/* tries to find violated linearization cuts and adds them to SCIP */
static
SCIP_RETCODE checkRowcuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< rowcuts constraint handler */
   SCIP_RESULT*          result,             /**< pointer to store result of checking process */
   SCIP_Bool             addcons             /**< shall row cuts be added as linear constraints? */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   CouenneProblem* problem;
   OsiCuts cs;
   SCIP_Real* sol;
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* get Couenne problem data */
   problem = conshdlrdata->cutgenerator-> Problem ();

   /* get variable data, create sol */
   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);
   sol = new SCIP_Real [nvars];
   SCIP_CALL( SCIPgetSolVals(scip, NULL, nvars, vars, sol) );

   /* store solution into MILP data structure */
   conshdlrdata -> milp -> setColSolution(sol);

   /* let Couenne generate linearization cuts */
   problem -> domain () -> push (problem -> nVars (), sol, NULL, NULL);
   conshdlrdata->cutgenerator->genRowCuts(*(conshdlrdata->milp), cs, 0, NULL);
   problem -> domain () -> pop  ();

   if( !addcons )
   {
      *result = (cs.sizeRowCuts() == 0) ? SCIP_FEASIBLE : SCIP_INFEASIBLE;
      return SCIP_OKAY;
   }
   
   for( i = 0; i < cs.sizeRowCuts(); i++ ) 
   {
      CoinPackedVector row;

      SCIP_CONS* cons;
      SCIP_VAR** vars;

      char consname[SCIP_MAXSTRLEN];  
      SCIP_Real* vals;
      int* idxs;

      SCIP_Real lhs;
      SCIP_Real rhs;

      int nvals;
      int j;

      /* get the row corresponding to the cut */
      lhs = cs.rowCut(i).lb();
      rhs = cs.rowCut(i).ub();      
      row = cs.rowCut(i).row();

      /* get row data */
      nvals = row.getNumElements();
      idxs = row.getIndices();
      vals = row.getElements(); 
      
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "rowcut", i);

      /* create an empty linear constraint */
      SCIP_CALL_ABORT( SCIPcreateConsLinear(scip, &cons, consname, 0, NULL, NULL, lhs, rhs, 
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      /* get SCIP variable array */
      vars = SCIPgetVars(scip);      

      /* add variables to constraint */
      for( j = 0; j < nvals; j++ )
      {
         SCIP_CALL( SCIPaddCoefLinear(scip, cons, vars[idxs[j]], vals[j]) );
      }
      
      /* add constraint to SCIP */
      SCIP_CALL( SCIPaddCons(scip, cons) );
#if 0
       SCIP_CALL( SCIPprintCons(scip,cons,NULL) );
#endif
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );        

      *result = SCIP_CONSADDED;
   }
   
   /* store cuts to MILP data structure */
   if (cs.sizeRowCuts ()) {
      conshdlrdata -> milp -> applyCuts (cs);
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#define conshdlrCopyRowcuts NULL

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeRowcuts)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   
   /* free constraint handler data */
   SCIPfreeMemory(scip, &conshdlrdata);
   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}

/** initialization method of constraint handler (called after problem was transformed) */
#define consInitRowcuts NULL

/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitRowcuts NULL

/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#define consInitpreRowcuts NULL

/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define consExitpreRowcuts NULL

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolRowcuts NULL

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#define consExitsolRowcuts NULL

/** frees specific constraint data */
#define consDeleteRowcuts NULL

/** transforms constraint data into data belonging to the transformed problem */ 
#define consTransRowcuts NULL

/** LP initialization method of constraint handler */
#define consInitlpRowcuts NULL

/** separation method of constraint handler for LP solutions */
#define consSepalpRowcuts NULL

/** separation method of constraint handler for arbitrary primal solutions */
#define consSepasolRowcuts NULL

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpRowcuts)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 && conss == NULL); /* there should be no constraints */
   assert(result != NULL);


   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   if( conshdlrdata->ncuttingrounds < conshdlrdata->maxcuttingrounds )
   {
      SCIP_CALL( checkRowcuts(scip, conshdlr, result, TRUE) );
      conshdlrdata->ncuttingrounds++;
   }

   /* we want to accept all solutions, even if we added a constraint that cuts them off */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsRowcuts)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 && conss == NULL); /* there should be no constraints */
   assert(result != NULL);


   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->ncuttingrounds < conshdlrdata->maxcuttingrounds )
   {
      SCIP_CALL( checkRowcuts(scip, conshdlr, result, TRUE) );
      conshdlrdata->ncuttingrounds++;
   }

   /* we want to accept all solutions, even if we added a constraint that cuts them off */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckRowcuts)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 && conss == NULL); /* there should be no constraints */
   assert(result != NULL);


   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING && conshdlrdata->ncuttingrounds < conshdlrdata->maxcuttingrounds )
   {
      SCIP_CALL( checkRowcuts(scip, conshdlr, result, TRUE) );
      conshdlrdata->ncuttingrounds++;
   }

   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
#define consPropRowcuts NULL

/** presolving method of constraint handler */
#define consPresolRowcuts NULL

/** propagation conflict resolving method of constraint handler */
#define consRespropRowcuts NULL

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockRowcuts)
{  /*lint --e{715}*/
   assert(cons == NULL);

   return SCIP_OKAY;
}

/** constraint activation notification method of constraint handler */
#define consActiveRowcuts NULL

/** constraint deactivation notification method of constraint handler */
#define consDeactiveRowcuts NULL

/** constraint enabling notification method of constraint handler */
#define consEnableRowcuts NULL

/** constraint disabling notification method of constraint handler */
#define consDisableRowcuts NULL

/** constraint display method of constraint handler */
#define consPrintRowcuts NULL

/** constraint copying method of constraint handler */
#define consCopyRowcuts NULL

/** constraint parsing method of constraint handler */
#define consParseRowcuts NULL

/*
 * constraint specific interface methods
 */

/** creates the handler for rowcuts constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrRowcuts(
   SCIP*                 scip,               /**< SCIP data structure */
   CouenneCutGenerator*  cutgenerator,       /**< CouenneCutGenerator for linearization cuts */
   OsiSolverInterface*   milp                /**< Couenne's MILP relaxation of Couenne's MINLP */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create rowcuts constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   conshdlrdata->cutgenerator = cutgenerator;
   conshdlrdata->milp = milp;
   conshdlrdata->ncuttingrounds = 0;

#if SCIP_VERSION >= 300
   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         SCIP_PROPTIMING_BEFORELP,
         conshdlrCopyRowcuts,
         consFreeRowcuts, consInitRowcuts, consExitRowcuts,
         consInitpreRowcuts, consExitpreRowcuts, consInitsolRowcuts, consExitsolRowcuts,
         consDeleteRowcuts, consTransRowcuts, consInitlpRowcuts,
         consSepalpRowcuts, consSepasolRowcuts, consEnfolpRowcuts, consEnfopsRowcuts, consCheckRowcuts,
         consPropRowcuts, consPresolRowcuts, consRespropRowcuts, consLockRowcuts,
         consActiveRowcuts, consDeactiveRowcuts,
         consEnableRowcuts, consDisableRowcuts, NULL,
         consPrintRowcuts, consCopyRowcuts, consParseRowcuts, NULL, NULL,
         conshdlrdata) );
#elif SCIP_VERSION >= 210
   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, 
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         SCIP_PROPTIMING_BEFORELP,
         conshdlrCopyRowcuts,
         consFreeRowcuts, consInitRowcuts, consExitRowcuts, 
         consInitpreRowcuts, consExitpreRowcuts, consInitsolRowcuts, consExitsolRowcuts,
         consDeleteRowcuts, consTransRowcuts, consInitlpRowcuts,
         consSepalpRowcuts, consSepasolRowcuts, consEnfolpRowcuts, consEnfopsRowcuts, consCheckRowcuts, 
         consPropRowcuts, consPresolRowcuts, consRespropRowcuts, consLockRowcuts,
         consActiveRowcuts, consDeactiveRowcuts, 
         consEnableRowcuts, consDisableRowcuts, NULL,
         consPrintRowcuts, consCopyRowcuts, consParseRowcuts,
         conshdlrdata) );
#else
   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         conshdlrCopyRowcuts,
         consFreeRowcuts, consInitRowcuts, consExitRowcuts,
         consInitpreRowcuts, consExitpreRowcuts, consInitsolRowcuts, consExitsolRowcuts,
         consDeleteRowcuts, consTransRowcuts, consInitlpRowcuts,
         consSepalpRowcuts, consSepasolRowcuts, consEnfolpRowcuts, consEnfopsRowcuts, consCheckRowcuts,
         consPropRowcuts, consPresolRowcuts, consRespropRowcuts, consLockRowcuts,
         consActiveRowcuts, consDeactiveRowcuts,
         consEnableRowcuts, consDisableRowcuts,
         consPrintRowcuts, consCopyRowcuts, consParseRowcuts,
         conshdlrdata) );
#endif

   /* add rowcuts constraint handler parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/"CONSHDLR_NAME"/maxcuttingrounds",
         "how many rounds of cuts should be applied at most?",
         &conshdlrdata->maxcuttingrounds, FALSE, DEFAULT_MAXCUTTINGROUNDS, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

#endif
