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
/*
 * This module allocates a memory `arena` with as many slots as OpenMP threads.
 *
 * The behaviour is as follows
 *  1) sol_driver.F calls init_workspaces
 *  2) sol_lr.F multiple threads call get_address
 *  3) sol_driver.F calls release_workspaces
 *
 * No concurrency mechanism is provided. The memory slot [i] is expected to be
 * accessed only by thread [i].
 *
 */
#include "mumps_sol_omp_memory_manager.h"
#include <stdalign.h>
#include <stdlib.h>
#include <stdio.h>
static uint8_t **arena = NULL;
static size_t nb_slots = 0;
static size_t surface = 0;
MUMPS_INT _init_workspaces(MUMPS_INT* n, MUMPS_INT p);
void MUMPS_CALL init_workspaces(MUMPS_INT* n, MUMPS_INT* p, MUMPS_INT* err) {
    *err = _init_workspaces( n, *p );
}
MUMPS_INT _init_workspaces(MUMPS_INT* n, MUMPS_INT p) {
    int i;
//    printf("[DEV] (beg) _init workspaces arena = %p\n", arena);
    if (arena != NULL) return 1;
    nb_slots = p;
    arena = malloc( nb_slots * sizeof(uint8_t*) );
    if (arena == NULL) return 1;
    for (i = 0; i < nb_slots; i++) {
        arena[i] = NULL;
    }
    MUMPS_INT error = 0;
    resize_workspaces(n, &error);
    return error;
}
void MUMPS_CALL resize_workspaces(MUMPS_INT* n, MUMPS_INT* err) {
    int i;
    surface = *n;
    *err = 0;
    for (i = 0; i < nb_slots; i++) {
        arena[i] = realloc(arena[i], surface * sizeof(uint8_t));
        if (arena[i] == NULL) {
            *err = 1;
            break;
        }
    }
}
void MUMPS_CALL release_workspaces() {
    int i;
//    printf("[DEV] release workspaces\n");
    for (i = 0; i < nb_slots; i++) {
        if (arena[i] == NULL) continue;
        free(arena[i]);
        arena[i] = NULL;
    }
    free(arena);
    arena = NULL;
}
//     FIXME UB id >= nb_slots
//     AJ I don't want to write this module "defensively"
//        though there could be something that is more defensive that 
//        allow to "get workspace of size X" and such a call would HAVE TO 
//        return a valid address (but you'd have to do a acquire/release
//        explicitly, just in case)
//        [ Get_workspace_of_size_X would be very nice for BLR temporaries ]
void MUMPS_CALL get_address(MUMPS_INT* id, intptr_t* out) {
//    printf("[DEV] arena[id = %ld] := %p = %ld (%ld bytes)\n", *id, arena[*id], arena[*id], surface);
    *out = (intptr_t) arena[*id];
}
