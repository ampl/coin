# Copyright (C) 2011-2012 Stefan Vigerske
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: Makefile.am 937 2010-10-02 19:05:12Z stefan $

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST =

lib_LTLIBRARIES = libcoinscip.la libcoinsciplpis.la

AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/scip/src` $(SCIPLIB_CFLAGS)

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)`

libcoinscip_la_LDFLAGS = $(LT_LDFLAGS)
libcoinsciplpis_la_LDFLAGS = $(LT_LDFLAGS)

libcoinscip_la_SOURCES = \
  scip/src/scip/bitencode.c \
  scip/src/scip/branch.c \
  scip/src/scip/buffer.c \
  scip/src/scip/clock.c \
  scip/src/scip/conflict.c \
  scip/src/scip/cons.c \
  scip/src/scip/cutpool.c \
  scip/src/scip/debug.c \
  scip/src/scip/dialog.c \
  scip/src/scip/disp.c \
  scip/src/scip/event.c \
  scip/src/scip/fileio.c \
  scip/src/scip/heur.c \
  scip/src/scip/history.c \
  scip/src/scip/implics.c \
  scip/src/scip/interrupt.c \
  scip/src/scip/intervalarith.c \
  scip/src/scip/lp.c \
  scip/src/scip/mem.c \
  scip/src/scip/message.c \
  scip/src/scip/misc.c \
  scip/src/scip/nlp.c \
  scip/src/scip/nodesel.c \
  scip/src/scip/paramset.c \
  scip/src/scip/presol.c \
  scip/src/scip/presolve.c \
  scip/src/scip/pricer.c \
  scip/src/scip/pricestore.c \
  scip/src/scip/primal.c \
  scip/src/scip/prob.c \
  scip/src/scip/prop.c \
  scip/src/scip/reader.c \
  scip/src/scip/relax.c \
  scip/src/scip/retcode.c \
  scip/src/scip/scip.c \
  scip/src/scip/scipdefplugins.c \
  scip/src/scip/scipgithash.c \
  scip/src/scip/scipshell.c \
  scip/src/scip/sepa.c \
  scip/src/scip/sepastore.c \
  scip/src/scip/set.c \
  scip/src/scip/sol.c \
  scip/src/scip/solve.c \
  scip/src/scip/stat.c \
  scip/src/scip/tree.c \
  scip/src/scip/var.c \
  scip/src/scip/vbc.c \
  scip/src/scip/branch_allfullstrong.c \
  scip/src/scip/branch_cloud.c \
  scip/src/scip/branch_fullstrong.c \
  scip/src/scip/branch_inference.c \
  scip/src/scip/branch_leastinf.c \
  scip/src/scip/branch_mostinf.c \
  scip/src/scip/branch_pscost.c \
  scip/src/scip/branch_random.c \
  scip/src/scip/branch_relpscost.c \
  scip/src/scip/cons_abspower.c \
  scip/src/scip/cons_and.c \
  scip/src/scip/cons_bivariate.c \
  scip/src/scip/cons_bounddisjunction.c \
  scip/src/scip/cons_conjunction.c \
  scip/src/scip/cons_countsols.c \
  scip/src/scip/cons_cumulative.c \
  scip/src/scip/cons_disjunction.c \
  scip/src/scip/cons_indicator.c \
  scip/src/scip/cons_integral.c \
  scip/src/scip/cons_knapsack.c \
  scip/src/scip/cons_linear.c \
  scip/src/scip/cons_linking.c \
  scip/src/scip/cons_logicor.c \
  scip/src/scip/cons_nonlinear.c \
  scip/src/scip/cons_or.c \
  scip/src/scip/cons_orbitope.c \
  scip/src/scip/cons_pseudoboolean.c \
  scip/src/scip/cons_quadratic.c \
  scip/src/scip/cons_setppc.c \
  scip/src/scip/cons_soc.c \
  scip/src/scip/cons_sos1.c \
  scip/src/scip/cons_sos2.c \
  scip/src/scip/cons_superindicator.c \
  scip/src/scip/cons_varbound.c \
  scip/src/scip/cons_xor.c \
  scip/src/scip/dialog_default.c \
  scip/src/scip/disp_default.c \
  scip/src/scip/heur_actconsdiving.c \
  scip/src/scip/heur_clique.c \
  scip/src/scip/heur_coefdiving.c \
  scip/src/scip/heur_crossover.c \
  scip/src/scip/heur_dins.c \
  scip/src/scip/heur_dualval.c \
  scip/src/scip/heur_feaspump.c \
  scip/src/scip/heur_fixandinfer.c \
  scip/src/scip/heur_fracdiving.c \
  scip/src/scip/heur_guideddiving.c \
  scip/src/scip/heur_intdiving.c \
  scip/src/scip/heur_intshifting.c \
  scip/src/scip/heur_linesearchdiving.c \
  scip/src/scip/heur_localbranching.c \
  scip/src/scip/heur_mutation.c \
  scip/src/scip/heur_nlpdiving.c \
  scip/src/scip/heur_objpscostdiving.c \
  scip/src/scip/heur_octane.c \
  scip/src/scip/heur_oneopt.c \
  scip/src/scip/heur_proximity.c \
  scip/src/scip/heur_pscostdiving.c \
  scip/src/scip/heur_randrounding.c \
  scip/src/scip/heur_rens.c \
  scip/src/scip/heur_rins.c \
  scip/src/scip/heur_rootsoldiving.c \
  scip/src/scip/heur_rounding.c \
  scip/src/scip/heur_shiftandpropagate.c \
  scip/src/scip/heur_shifting.c \
  scip/src/scip/heur_simplerounding.c \
  scip/src/scip/heur_subnlp.c \
  scip/src/scip/heur_trivial.c \
  scip/src/scip/heur_trysol.c \
  scip/src/scip/heur_twoopt.c \
  scip/src/scip/heur_undercover.c \
  scip/src/scip/heur_vbounds.c \
  scip/src/scip/heur_veclendiving.c \
  scip/src/scip/heur_zeroobj.c \
  scip/src/scip/heur_zirounding.c \
  scip/src/scip/message_default.c \
  scip/src/scip/nodesel_bfs.c \
  scip/src/scip/nodesel_breadthfirst.c \
  scip/src/scip/nodesel_dfs.c \
  scip/src/scip/nodesel_estimate.c \
  scip/src/scip/nodesel_hybridestim.c \
  scip/src/scip/nodesel_restartdfs.c \
  scip/src/scip/nodesel_uct.c \
  scip/src/scip/presol_boundshift.c \
  scip/src/scip/presol_components.c \
  scip/src/scip/presol_convertinttobin.c \
  scip/src/scip/presol_domcol.c \
  scip/src/scip/presol_dualinfer.c \
  scip/src/scip/presol_gateextraction.c \
  scip/src/scip/presol_implics.c \
  scip/src/scip/presol_inttobinary.c \
  scip/src/scip/presol_trivial.c \
  scip/src/scip/prop_dualfix.c \
  scip/src/scip/prop_genvbounds.c \
  scip/src/scip/prop_obbt.c \
  scip/src/scip/prop_probing.c \
  scip/src/scip/prop_pseudoobj.c \
  scip/src/scip/prop_redcost.c \
  scip/src/scip/prop_rootredcost.c \
  scip/src/scip/prop_vbounds.c \
  scip/src/scip/reader_bnd.c \
  scip/src/scip/reader_ccg.c \
  scip/src/scip/reader_cip.c \
  scip/src/scip/reader_cnf.c \
  scip/src/scip/reader_fix.c \
  scip/src/scip/reader_fzn.c \
  scip/src/scip/reader_gms.c \
  scip/src/scip/reader_lp.c \
  scip/src/scip/reader_mps.c \
  scip/src/scip/reader_opb.c \
  scip/src/scip/reader_osil.c \
  scip/src/scip/reader_pip.c \
  scip/src/scip/reader_pbm.c \
  scip/src/scip/reader_ppm.c \
  scip/src/scip/reader_rlp.c \
  scip/src/scip/reader_sol.c \
  scip/src/scip/reader_wbo.c \
  scip/src/scip/reader_zpl.c \
  scip/src/scip/sepa_cgmip.c \
  scip/src/scip/sepa_clique.c \
  scip/src/scip/sepa_closecuts.c \
  scip/src/scip/sepa_cmir.c \
  scip/src/scip/sepa_flowcover.c \
  scip/src/scip/sepa_gomory.c \
  scip/src/scip/sepa_impliedbounds.c \
  scip/src/scip/sepa_intobj.c \
  scip/src/scip/sepa_mcf.c \
  scip/src/scip/sepa_oddcycle.c \
  scip/src/scip/sepa_rapidlearning.c \
  scip/src/scip/sepa_strongcg.c \
  scip/src/scip/sepa_zerohalf.c \
  scip/src/tclique/tclique_branch.c \
  scip/src/tclique/tclique_coloring.c \
  scip/src/tclique/tclique_graph.c \
  scip/src/blockmemshell/memory.c \
  scip/src/nlpi/expr.c \
  scip/src/nlpi/nlpi.c \
  scip/src/nlpi/nlpioracle.c \
  scip/src/nlpi/exprinterpret_cppad.cpp \
  scip/src/dijkstra/dijkstra.c \
  scip/src/xml/xmlparse.c

# LP solver interfaces
libcoinsciplpis_la_SOURCES = lpiswitch.c lpiswitch_init.cpp lpiswitch_none.c
if COIN_HAS_SOPLEX
  # for static builds, we need bitencode also in the lpis library,
  # because it is listed after the scip lib on the linker line but
  # has functions that only lpi_spx is using
  libcoinsciplpis_la_SOURCES += scip/src/scip/bitencode.c
  libcoinsciplpis_la_SOURCES += lpiswitch_spx.cpp lpiswitch_spx2.cpp
endif
if COIN_HAS_CLP
  libcoinsciplpis_la_SOURCES += lpiswitch_clp.cpp
endif
if COIN_HAS_CPX
  libcoinsciplpis_la_SOURCES += lpiswitch_cpx.c
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CPXINCDIR)`
endif
if COIN_HAS_MSK
  libcoinsciplpis_la_SOURCES += lpiswitch_msk.c
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(MSKINCDIR)`
endif
if COIN_HAS_XPR
  libcoinsciplpis_la_SOURCES += lpiswitch_xprs.c
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(XPRINCDIR)`
endif
if COIN_HAS_GRB
  libcoinsciplpis_la_SOURCES += lpiswitch_grb.c
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(GRBINCDIR)`
endif
if COIN_HAS_QSO
  libcoinsciplpis_la_SOURCES += lpiswitch_qso.c
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(QSOINCDIR)`
endif


# add IPOPT as nonlinear solver, if available
if COIN_HAS_IPOPT
libcoinscip_la_SOURCES += scip/src/nlpi/nlpi_ipopt.cpp
else
libcoinscip_la_SOURCES += scip/src/nlpi/nlpi_ipopt_dummy.c
endif

thirdpartyincludedir = $(includedir)/coin/ThirdParty
thirdpartyinclude_HEADERS = lpiswitch.h

thirdpartyscipincludedir = $(includedir)/coin/ThirdParty/scip
thirdpartyscipinclude_HEADERS = \
  scip/src/scip/bitencode.h \
  scip/src/scip/branch_allfullstrong.h \
  scip/src/scip/branch_cloud.h \
  scip/src/scip/branch_fullstrong.h \
  scip/src/scip/branch_inference.h \
  scip/src/scip/branch_leastinf.h \
  scip/src/scip/branch_mostinf.h \
  scip/src/scip/branch_pscost.h \
  scip/src/scip/branch_random.h \
  scip/src/scip/branch_relpscost.h \
  scip/src/scip/buffer.h \
  scip/src/scip/clock.h \
  scip/src/scip/cons_abspower.h \
  scip/src/scip/cons_and.h \
  scip/src/scip/cons_bivariate.h \
  scip/src/scip/cons_bounddisjunction.h \
  scip/src/scip/cons_conjunction.h \
  scip/src/scip/cons_countsols.h \
  scip/src/scip/cons_cumulative.h \
  scip/src/scip/cons_disjunction.h \
  scip/src/scip/cons_indicator.h \
  scip/src/scip/cons_integral.h \
  scip/src/scip/cons_knapsack.h \
  scip/src/scip/cons_linear.h \
  scip/src/scip/cons_linking.h \
  scip/src/scip/cons_logicor.h \
  scip/src/scip/cons_nonlinear.h \
  scip/src/scip/cons_or.h \
  scip/src/scip/cons_orbitope.h \
  scip/src/scip/cons_pseudoboolean.h \
  scip/src/scip/cons_quadratic.h \
  scip/src/scip/cons_setppc.h \
  scip/src/scip/cons_soc.h \
  scip/src/scip/cons_sos1.h \
  scip/src/scip/cons_sos2.h \
  scip/src/scip/cons_superindicator.h \
  scip/src/scip/cons_varbound.h \
  scip/src/scip/cons_xor.h \
  scip/src/scip/def.h \
  scip/src/scip/dialog_default.h \
  scip/src/scip/disp_default.h \
  scip/src/scip/event.h \
  scip/src/scip/heur_actconsdiving.h \
  scip/src/scip/heur_clique.h \
  scip/src/scip/heur_coefdiving.h \
  scip/src/scip/heur_crossover.h \
  scip/src/scip/heur_dins.h \
  scip/src/scip/heur_dualval.h \
  scip/src/scip/heur_feaspump.h \
  scip/src/scip/heur_fixandinfer.h \
  scip/src/scip/heur_fracdiving.h \
  scip/src/scip/heur_guideddiving.h \
  scip/src/scip/heur_intdiving.h \
  scip/src/scip/heur_intshifting.h \
  scip/src/scip/heur_linesearchdiving.h \
  scip/src/scip/heur_localbranching.h \
  scip/src/scip/heur_mutation.h \
  scip/src/scip/heur_nlpdiving.h \
  scip/src/scip/heur_objpscostdiving.h \
  scip/src/scip/heur_octane.h \
  scip/src/scip/heur_oneopt.h \
  scip/src/scip/heur_proximity.h \
  scip/src/scip/heur_pscostdiving.h \
  scip/src/scip/heur_randrounding.h \
  scip/src/scip/heur_rens.h \
  scip/src/scip/heur_rins.h \
  scip/src/scip/heur_rootsoldiving.h \
  scip/src/scip/heur_rounding.h \
  scip/src/scip/heur_shiftandpropagate.h \
  scip/src/scip/heur_shifting.h \
  scip/src/scip/heur_simplerounding.h \
  scip/src/scip/heur_subnlp.h \
  scip/src/scip/heur_trivial.h \
  scip/src/scip/heur_trysol.h \
  scip/src/scip/heur_twoopt.h \
  scip/src/scip/heur_undercover.h \
  scip/src/scip/heur_vbounds.h \
  scip/src/scip/heur_veclendiving.h \
  scip/src/scip/heur_zeroobj.h \
  scip/src/scip/heur_zirounding.h \
  scip/src/scip/history.h \
  scip/src/scip/implics.h \
  scip/src/scip/intervalarith.h \
  scip/src/scip/message.h \
  scip/src/scip/message_default.h \
  scip/src/scip/misc.h \
  scip/src/scip/nodesel_bfs.h \
  scip/src/scip/nodesel_breadthfirst.h \
  scip/src/scip/nodesel_dfs.h \
  scip/src/scip/nodesel_estimate.h \
  scip/src/scip/nodesel_hybridestim.h \
  scip/src/scip/nodesel_restartdfs.h \
  scip/src/scip/nodesel_uct.h \
  scip/src/scip/presolve.h \
  scip/src/scip/presol_boundshift.h \
  scip/src/scip/presol_components.h \
  scip/src/scip/presol_convertinttobin.h \
  scip/src/scip/presol_domcol.h \
  scip/src/scip/presol_dualinfer.h \
  scip/src/scip/presol_gateextraction.h \
  scip/src/scip/presol_implics.h \
  scip/src/scip/presol_inttobinary.h \
  scip/src/scip/presol_trivial.h \
  scip/src/scip/prop_dualfix.h \
  scip/src/scip/prop_genvbounds.h \
  scip/src/scip/prop_obbt.h \
  scip/src/scip/prop_probing.h \
  scip/src/scip/prop_pseudoobj.h \
  scip/src/scip/prop_redcost.h \
  scip/src/scip/prop_rootredcost.h \
  scip/src/scip/prop_vbounds.h \
  scip/src/scip/pub_branch.h \
  scip/src/scip/pub_conflict.h \
  scip/src/scip/pub_cons.h \
  scip/src/scip/pub_cutpool.h \
  scip/src/scip/pub_dialog.h \
  scip/src/scip/pub_disp.h \
  scip/src/scip/pub_event.h \
  scip/src/scip/pub_fileio.h \
  scip/src/scip/pub_heur.h \
  scip/src/scip/pub_history.h \
  scip/src/scip/pub_implics.h \
  scip/src/scip/pub_lp.h \
  scip/src/scip/pub_message.h \
  scip/src/scip/pub_misc.h \
  scip/src/scip/pub_nlp.h \
  scip/src/scip/pub_nodesel.h \
  scip/src/scip/pub_paramset.h \
  scip/src/scip/pub_presol.h \
  scip/src/scip/pub_pricer.h \
  scip/src/scip/pub_prop.h \
  scip/src/scip/pub_reader.h \
  scip/src/scip/pub_relax.h \
  scip/src/scip/pub_sepa.h \
  scip/src/scip/pub_sol.h \
  scip/src/scip/pub_tree.h \
  scip/src/scip/pub_var.h \
  scip/src/scip/reader_bnd.h \
  scip/src/scip/reader_ccg.h \
  scip/src/scip/reader_cip.h \
  scip/src/scip/reader_cnf.h \
  scip/src/scip/reader_fix.h \
  scip/src/scip/reader_fzn.h \
  scip/src/scip/reader_gms.h \
  scip/src/scip/reader_lp.h \
  scip/src/scip/reader_mps.h \
  scip/src/scip/reader_opb.h \
  scip/src/scip/reader_osil.h \
  scip/src/scip/reader_pip.h \
  scip/src/scip/reader_pbm.h \
  scip/src/scip/reader_ppm.h \
  scip/src/scip/reader_rlp.h \
  scip/src/scip/reader_sol.h \
  scip/src/scip/reader_wbo.h \
  scip/src/scip/reader_zpl.h \
  scip/src/scip/scipdefplugins.h \
  scip/src/scip/scip.h \
  scip/src/scip/scipshell.h \
  scip/src/scip/sepa_cgmip.h \
  scip/src/scip/sepa_clique.h \
  scip/src/scip/sepa_closecuts.h \
  scip/src/scip/sepa_cmir.h \
  scip/src/scip/sepa_flowcover.h \
  scip/src/scip/sepa_gomory.h \
  scip/src/scip/sepa_impliedbounds.h \
  scip/src/scip/sepa_intobj.h \
  scip/src/scip/sepa_mcf.h \
  scip/src/scip/sepa_oddcycle.h \
  scip/src/scip/sepa_rapidlearning.h \
  scip/src/scip/sepa_strongcg.h \
  scip/src/scip/sepa_zerohalf.h \
  scip/src/scip/set.h \
  scip/src/scip/struct_branch.h \
  scip/src/scip/struct_buffer.h \
  scip/src/scip/struct_clock.h \
  scip/src/scip/struct_conflict.h \
  scip/src/scip/struct_cons.h \
  scip/src/scip/struct_cutpool.h \
  scip/src/scip/struct_dialog.h \
  scip/src/scip/struct_disp.h \
  scip/src/scip/struct_event.h \
  scip/src/scip/struct_heur.h \
  scip/src/scip/struct_history.h \
  scip/src/scip/struct_implics.h \
  scip/src/scip/struct_lp.h \
  scip/src/scip/struct_mem.h \
  scip/src/scip/struct_message.h \
  scip/src/scip/struct_misc.h \
  scip/src/scip/struct_nlp.h \
  scip/src/scip/struct_nodesel.h \
  scip/src/scip/struct_paramset.h \
  scip/src/scip/struct_presol.h \
  scip/src/scip/struct_pricer.h \
  scip/src/scip/struct_pricestore.h \
  scip/src/scip/struct_primal.h \
  scip/src/scip/struct_prob.h \
  scip/src/scip/struct_prop.h \
  scip/src/scip/struct_reader.h \
  scip/src/scip/struct_relax.h \
  scip/src/scip/struct_scip.h \
  scip/src/scip/struct_sepa.h \
  scip/src/scip/struct_sepastore.h \
  scip/src/scip/struct_set.h \
  scip/src/scip/struct_sol.h \
  scip/src/scip/struct_stat.h \
  scip/src/scip/struct_tree.h \
  scip/src/scip/struct_var.h \
  scip/src/scip/struct_vbc.h \
  scip/src/scip/tree.h \
  scip/src/scip/type_branch.h \
  scip/src/scip/type_buffer.h \
  scip/src/scip/type_clock.h \
  scip/src/scip/type_conflict.h \
  scip/src/scip/type_cons.h \
  scip/src/scip/type_cutpool.h \
  scip/src/scip/type_dialog.h \
  scip/src/scip/type_disp.h \
  scip/src/scip/type_event.h \
  scip/src/scip/type_heur.h \
  scip/src/scip/type_history.h \
  scip/src/scip/type_implics.h \
  scip/src/scip/type_interrupt.h \
  scip/src/scip/type_lp.h \
  scip/src/scip/type_mem.h \
  scip/src/scip/type_message.h \
  scip/src/scip/type_misc.h \
  scip/src/scip/type_nlp.h \
  scip/src/scip/type_nodesel.h \
  scip/src/scip/type_paramset.h \
  scip/src/scip/type_presol.h \
  scip/src/scip/type_pricer.h \
  scip/src/scip/type_pricestore.h \
  scip/src/scip/type_primal.h \
  scip/src/scip/type_prob.h \
  scip/src/scip/type_prop.h \
  scip/src/scip/type_reader.h \
  scip/src/scip/type_relax.h \
  scip/src/scip/type_result.h \
  scip/src/scip/type_retcode.h \
  scip/src/scip/type_scip.h \
  scip/src/scip/type_sepa.h \
  scip/src/scip/type_sepastore.h \
  scip/src/scip/type_set.h \
  scip/src/scip/type_sol.h \
  scip/src/scip/type_stat.h \
  scip/src/scip/type_timing.h \
  scip/src/scip/type_tree.h \
  scip/src/scip/type_var.h \
  scip/src/scip/type_vbc.h \
  scip/src/scip/var.h \
  scip/src/scip/githash.c

thirdpartybmsincludedir = $(includedir)/coin/ThirdParty/blockmemshell
thirdpartybmsinclude_HEADERS = scip/src/blockmemshell/memory.h

thirdpartynlpiincludedir = $(includedir)/coin/ThirdParty/nlpi
thirdpartynlpiinclude_HEADERS = \
  scip/src/nlpi/exprinterpret.h \
  scip/src/nlpi/nlpi_ipopt.h \
  scip/src/nlpi/pub_expr.h \
  scip/src/nlpi/struct_expr.h \
  scip/src/nlpi/struct_nlpi.h \
  scip/src/nlpi/type_expr.h \
  scip/src/nlpi/type_exprinterpret.h \
  scip/src/nlpi/type_nlpi.h

thirdpartylpiincludedir = $(includedir)/coin/ThirdParty/lpi
thirdpartylpiinclude_HEADERS = \
  scip/src/lpi/lpi.h \
  scip/src/lpi/type_lpi.h

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = coinscip.pc

test:
	@echo "No test available for SCIP."

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include BuildTools/Makemain.inc
