/*
 *
 *  This file is part of MUMPS 5.9.1, released
 *  on Mon Jul 20 09:00:43 UTC 2026
 *
 *
 *  Copyright 1991-2026 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  Mumps Technologies, University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license 
 *  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
 *  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
 *
 */
/* $Id */
#ifndef MUMPS_SOL_OMP_MEM_MANAGER
#define MUMPS_SOL_OMP_MEM_MANAGER
/*
 * FIXME no error handling
 * FIXME no input sanitation
 */
#include <stdlib.h>
#include "mumps_common.h"
#define init_workspaces     F_SYMBOL(sol_omp_mem_manager_init_workspaces,       SOL_OMP_MEM_MANAGER_INIT_WORKSPACES)
#define resize_workspaces   F_SYMBOL(sol_omp_mem_manager_resize_workspaces,     SOL_OMP_MEM_MANAGER_RESIZE_WORKSPACES)
#define release_workspaces  F_SYMBOL(sol_omp_mem_manager_release_workspaces,    SOL_OMP_MEM_MANAGER_RELEASE_WORKSPACES)
#define get_address         F_SYMBOL(sol_omp_mem_manager_get_address,           SOL_OMP_MEM_MANAGER_GET_ADDRESS)
void MUMPS_CALL init_workspaces(MUMPS_INT* n, MUMPS_INT* p, MUMPS_INT* err);
void MUMPS_CALL resize_workspaces(MUMPS_INT* n, MUMPS_INT* err);
void MUMPS_CALL release_workspaces();
void MUMPS_CALL get_address(MUMPS_INT* id, intptr_t* out);
#endif
