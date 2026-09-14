/*
 -----------------------------------------------------------------------
 NEC Vector Host MUMPS Library Interface
 -----------------------------------------------------------------------
 */

#if ! defined VE_METIS_H
#define VE_METIS_H

/* Includes */

#include <stddef.h>
#include <libvhcall.h>
#include <libgen.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <metis.h>

/* Macros */

#define NECVH_MAX_LIBNAME 2048

#define mumps_abort mumps_abort_

void mumps_abort(void);

/* Functions prototype */

int VE_Metis_setdefaultoptions( int * );
int VE_Metis_nodend( int * , int * , int * , int * ,  int * , int * , int * );

#endif
