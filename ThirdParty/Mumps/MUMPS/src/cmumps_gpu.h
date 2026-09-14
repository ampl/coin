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
#ifndef CMUMPS_GPU_H
#define CMUMPS_GPU_H
#include "mumps_compat.h"
#include "mumps_common.h"
#if defined(USE_GPU)
# if defined(UPPER) || defined(MUMPS_WIN32)
#  define F_SYM_ARITH_GPU(lower_case,upper_case) CMUMPS_##upper_case
# elif defined(Add_)
#  define F_SYM_ARITH_GPU(lower_case,upper_case) cmumps_##lower_case##_
# elif defined(Add__)
#  define F_SYM_ARITH_GPU(lower_case,upper_case) cmumps_##lower_case##__
# else
#  define F_SYM_ARITH_GPU(lower_case,upper_case) cmumps_##lower_case
# endif
#define CMUMPS_GET_NDEVICES_C \
  F_SYM_ARITH_GPU(get_ndevices_c,GET_NDEVICES_C)
void MUMPS_CALL CMUMPS_GET_NDEVICES_C(MUMPS_INT *NDEVICES);
#define CMUMPS_GET_DEVICE_PROPERTIES_C \
  F_SYM_ARITH_GPU(get_device_properties_c,GET_DEVICE_PROPERTIES_C)
void MUMPS_CALL CMUMPS_GET_DEVICE_PROPERTIES_C(MUMPS_INT *ARCH, MUMPS_INT *HAS_UNIFIED);
#define CMUMPS_INIT_GPU_C \
  F_SYM_ARITH_GPU(init_gpu_c,INIT_GPU_C)
void MUMPS_CALL
CMUMPS_INIT_GPU_C(MUMPS_INT *K415, MUMPS_INT *K417, MUMPS_INT *K418,
                  MUMPS_INT *K419, MUMPS_INT *K422, MUMPS_INT *K430,
                  MUMPS_INT *K436,
                  MUMPS_INT *GPUBITMAP, MUMPS_INT *MAX_MPIs_PER_GPU,
                  MUMPS_INT info[2]
);
#define CMUMPS_END_GPU_C \
  F_SYM_ARITH_GPU(end_gpu_c,END_GPU_C)
MUMPS_INT MUMPS_CALL
CMUMPS_END_GPU_C(MUMPS_INT *K400, MUMPS_INT *K422
                  , CMUMPS_REAL *flops_gpu
                  , CMUMPS_REAL *time_gpu
                );
#define CMUMPS_GEMM_ON_GPU_C \
  F_SYM_ARITH_GPU(gemm_on_gpu_c,GEMM_ON_GPU_C)
MUMPS_INT MUMPS_CALL
CMUMPS_GEMM_ON_GPU_C(MUMPS_INT *GPU_ACTIVATION, MUMPS_INT *M, MUMPS_INT *N, MUMPS_INT *K);
#define CMUMPS_GEMM_GPU_C \
  F_SYM_ARITH_GPU(gemm_gpu_c,GEMM_GPU_C)
void MUMPS_CALL
CMUMPS_GEMM_GPU_C(MUMPS_INT *GPU_OPTION, char* TransA, char* TransB, MUMPS_INT *M, MUMPS_INT *N, MUMPS_INT *K, CMUMPS_COMPLEX *ALPHA, CMUMPS_COMPLEX *A,
    MUMPS_INT *LDA, CMUMPS_COMPLEX *B, MUMPS_INT *LDB, CMUMPS_COMPLEX *BETA, CMUMPS_COMPLEX *C, MUMPS_INT *LDC, mumps_ftnlen l1);
#define CMUMPS_TRSM_ON_GPU_C \
  F_SYM_ARITH_GPU(trsm_on_gpu_c,TRSM_ON_GPU_C)
MUMPS_INT MUMPS_CALL
CMUMPS_TRSM_ON_GPU_C(MUMPS_INT *GPU_ACTIVATION, MUMPS_INT *M, MUMPS_INT *N, char* Side);
#define CMUMPS_TRSM_GPU_C \
  F_SYM_ARITH_GPU(trsm_gpu_c,TRSM_GPU_C)
void MUMPS_CALL
CMUMPS_TRSM_GPU_C(MUMPS_INT *gpu_async, char* Side, char* Uplo, char* TransA, char* Diag, MUMPS_INT *M, MUMPS_INT *N, CMUMPS_COMPLEX *ALPHA, CMUMPS_COMPLEX *A,
    MUMPS_INT *LDA, CMUMPS_COMPLEX *B, MUMPS_INT *LDB, mumps_ftnlen l1);
/* Synchronous pinning */
#define CMUMPS_GPU_PIN_MEMORY_C \
  F_SYM_ARITH_GPU(gpu_pin_memory_c,GPU_PIN_MEMORY_C)
void MUMPS_CALL
CMUMPS_GPU_PIN_MEMORY_C(CMUMPS_COMPLEX *A, MUMPS_INT8 *ADIM);
/* Synchronous unpinning */
#define CMUMPS_GPU_UNPIN_MEMORY_C \
  F_SYM_ARITH_GPU(gpu_unpin_memory_c,GPU_UNPIN_MEMORY_C)
void MUMPS_CALL
CMUMPS_GPU_UNPIN_MEMORY_C(CMUMPS_COMPLEX *A, MUMPS_INT8 *ADIM);
/* used for both GEMMT and SYRKX: */
#  define CMUMPS_GEMMT_ON_GPU_C \
   F_SYM_ARITH_GPU(gemmt_on_gpu_c,GEMMT_ON_GPU_C)
   MUMPS_INT MUMPS_CALL
   CMUMPS_GEMMT_ON_GPU_C(MUMPS_INT *GPU_ACTIVATION, MUMPS_INT *N, MUMPS_INT *K);
#if defined(GEMMT_AVAILABLE_FOR_GPU) || defined(USE_XKBLAS)
/* GEMMT */
#  define CMUMPS_GEMMT_GPU_C \
   F_SYM_ARITH_GPU(gemmt_gpu_c,GEMMT_GPU_C)
   void MUMPS_CALL
   CMUMPS_GEMMT_GPU_C(MUMPS_INT *GPU_ACTIVATION, char* Uplo, char* TransA, char* TransB, MUMPS_INT *N, MUMPS_INT *K, CMUMPS_COMPLEX *ALPHA, CMUMPS_COMPLEX *A, MUMPS_INT *LDA, CMUMPS_COMPLEX *B, MUMPS_INT *LDB, CMUMPS_COMPLEX *BETA, CMUMPS_COMPLEX *C, MUMPS_INT *LDC, mumps_ftnlen l1);
#else
#  define CMUMPS_SYRKX_GPU_C \
   F_SYM_ARITH_GPU(syrkx_gpu_c,SYRKX_GPU_C)
   void MUMPS_CALL
   CMUMPS_SYRKX_GPU_C(MUMPS_INT *GPU_ACTIVATION, char* Uplo, char* Trans, MUMPS_INT *N, MUMPS_INT *K, CMUMPS_COMPLEX *ALPHA, CMUMPS_COMPLEX *A, MUMPS_INT *LDA, CMUMPS_COMPLEX *B, MUMPS_INT *LDB, CMUMPS_COMPLEX *BETA, CMUMPS_COMPLEX *C, MUMPS_INT *LDC, mumps_ftnlen l1);
#  define CMUMPS_SYMSCHUR_ON_GPU_C \
   F_SYM_ARITH_GPU(symschur_on_gpu_c,SYMSCHUR_ON_GPU_C)
   MUMPS_INT MUMPS_CALL
   CMUMPS_SYMSCHUR_ON_GPU_C(MUMPS_INT *GPU_ACTIVATION, MUMPS_INT *N, MUMPS_INT *K, MUMPS_INT *BLOCKING);
#endif
#if defined(USE_XKBLAS)
/* Allocation and free routines */
#  define CMUMPS_GPU_MALLOC_C \
   F_SYM_ARITH_GPU(gpu_malloc_c,GPU_MALLOC_C)
void MUMPS_CALL CMUMPS_GPU_MALLOC_C(MUMPS_INT8 *address, MUMPS_INT8 *size);
#  define CMUMPS_GPU_FREE_C \
   F_SYM_ARITH_GPU(gpu_free_c,GPU_FREE_C)
void MUMPS_CALL CMUMPS_GPU_FREE_C(MUMPS_INT8 *address);
#if defined(USE_XKBLAS_COPYSCALE)
#  define CMUMPS_GPU_COPY2U_SCALEL_C \
   F_SYM_ARITH_GPU(gpu_copy2u_scalel_c,GPU_COPY2U_SCALEL_C)
   void MUMPS_CALL CMUMPS_GPU_COPY2U_SCALEL_C(MUMPS_INT *gpu_async, MUMPS_INT *irowmax, MUMPS_INT *irowmin, MUMPS_INT *copysize, MUMPS_INT *nfront, MUMPS_INT *npiv, MUMPS_INT *liw, MUMPS_INT *IW, MUMPS_INT *offset_iw, MUMPS_INT8 *la, CMUMPS_COMPLEX *a, MUMPS_INT8 *dpos, MUMPS_INT8 *lpos, MUMPS_INT8 *upos, MUMPS_INT *should_copy);
#endif
#endif
/* Asynchronous pinning */
#define CMUMPS_GPU_ASYNC_PIN_MEMORY_C \
  F_SYM_ARITH_GPU(gpu_async_pin_memory_c,GPU_ASYNC_PIN_MEMORY_C)
void MUMPS_CALL
CMUMPS_GPU_ASYNC_PIN_MEMORY_C(CMUMPS_COMPLEX *A, MUMPS_INT8 *ADIM);
/* Wait all */
#define CMUMPS_GPU_WAIT_ENDOFPINNING_C \
  F_SYM_ARITH_GPU(gpu_wait_endofpinning_c,GPU_WAIT_ENDOFPINNING_C)
void MUMPS_CALL
CMUMPS_GPU_WAIT_ENDOFPINNING_C();
/* Asynchronous unpinning */
#define CMUMPS_GPU_ASYNC_UNPIN_MEMORY_C \
  F_SYM_ARITH_GPU(gpu_async_unpin_memory_c,GPU_ASYNC_UNPIN_MEMORY_C)
void MUMPS_CALL
CMUMPS_GPU_ASYNC_UNPIN_MEMORY_C(CMUMPS_COMPLEX *A, MUMPS_INT8 *ADIM);
#define CMUMPS_GPU_MODIFY_ARITH_TRHRESHOLD \
  F_SYM_ARITH_GPU(gpu_modify_arith_trhreshold,GPU_MODIFY_ARITH_TRHRESHOLD)
void MUMPS_CALL
CMUMPS_GPU_MODIFY_ARITH_TRHRESHOLD(MUMPS_INT *OLD, MUMPS_INT *NEW);
#define CMUMPS_GPU_SET_UNDER_L0_STATUS_C \
  F_SYM_ARITH_GPU(gpu_set_under_l0_status_c,GPU_SET_UNDER_L0_STATUS_C)
void MUMPS_CALL
CMUMPS_GPU_SET_UNDER_L0_STATUS_C(MUMPS_INT *VALEUR);
#else
void MUMPS_CALL cmumps_gpu_return();
#endif
#endif /* CMUMPS_GPU_H */
