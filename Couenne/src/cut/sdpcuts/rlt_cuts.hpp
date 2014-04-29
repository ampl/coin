/* $Id: rlt_cuts.hpp 508 2011-02-15 21:52:44Z pbelotti $
 *
 * Name:    rlt_cuts.hpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef RLT_CUTS_HPP
#define RLT_CUTS_HPP

#include <CglCutGenerator.hpp>
#include <tracer.hpp>


void rltCutsGen(const double *sol, int n, OsiCuts &cs, double *lb, double *ub, int m, Tracer *tracer);

#endif
