/*
 -----------------------------------------------------------------------
 NEC Vector Host MUMPS Library Interface
 -----------------------------------------------------------------------
 */

#if ! defined VE_MUMPS_H
#define VE_MUMPS_H

/* Includes */

#include <stddef.h>
#include <libvhcall.h>
#include <libgen.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

/* Macros */

#define NECVH_MAX_LIBNAME 2048

/* MUMPS functions called prototype */

#define mumps_abort mumps_abort_

void mumps_abort(void);

/* Functions prototype */

void VE_Mumps_ana_h( int TOTEL       , bool COMPUTE_PERM , int N         , int64_t IWLEN , int64_t * PE ,
                     int64_t * PFREE , int * LEN         , int * IW      , int * NV      , int * ELEN   ,
                     int * LAST      , int * NCMPA       , int * DEGREE  , int * HEAD    , int * NEXT   ,
                     int * W         , int * PARENT      );

void VE_Mumps_symqamd( int THRESH   , int * NDENSE    , int N       , int TOTEL           , int64_t IWLEN  ,                  
                       int64_t * PE , int64_t * PFREE , int * LEN   , int * IW            , int * NV       ,                  
                       int * ELEN   , int * LAST      , int * NCMPA , int * DEGREE        , int * HEAD     ,                  
                       int * NEXT   , int * W         , int * PERM  , int * LISTVAR_SCHUR , int SIZE_SCHUR ,                  
                       int AGG6     , int * PARENT    );


void VE_Mumps_wrap_ginp94( int N               , int64_t * IPE , int * IW             , int64_t   LIW8  , int * PERM    ,
                           int * SizeOfBlocks  , int KEEP60    , int * LISTVAR_SCHUR  , int SIZE_SCHUR  , int KEEP378   ,
                           int * COLCOUNT      , int * PARENT  , int * PORDER         , int * IWTMP1    , int * IWTMP2  ,
                           int * IWTMP3        , int * IWTMP4  , int * IWTMP5         , int * INFO      );
#endif
