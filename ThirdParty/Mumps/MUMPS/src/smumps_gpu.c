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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "smumps_gpu.h"
#if defined(USE_GPU)
#if defined(USE_XKBLAS)
#include "xkblas.h"
#if ( __XKBLAS__ >= 1 || ( __XKBLAS__ == 0 && __XKBLAS_MINOR__ >= 4 ) )
/* We no longer rely on cuda_runtime_api.h */
#else
#include "cuda_runtime_api.h"
#endif
#else
#include "cublasXt.h"
#include "cuda_runtime_api.h"
#endif
static int smumps_arith_threshold;
static int smumps_arith_threshold_above, smumps_arith_threshold_under;
static int smumps_gpu_under_l0;
#if ! defined(USE_XKBLAS)
/* for asynchronous pinning */
# include "mumps_register_thread.h"
static cublasXtHandle_t smumps_handle;
#endif
#include <mpi.h>
/* dev_version... ok not to prefix the global symbols below with smumps_
 * since we rarely mix different arithmetics with dev_version */
double smumps_wtime_pin, smumps_wtime_async_pin, smumps_wtime_gemm, smumps_wtime_gemmt, smumps_wtime_trsm, smumps_wtime_unpin,smumps_wtime_copyscale;
double smumps_size_pin, smumps_size_async_pin, smumps_size_async_unpin, smumps_flops_gemm, smumps_flops_gemmt, smumps_flops_trsm, smumps_flops_copyscale;
double smumps_start_async_pin;
/* for L0-OMP: */
double smumps_wtime_l0_gemm, smumps_wtime_l0_gemmt, smumps_wtime_l0_trsm, smumps_wtime_l0_copyscale;
double smumps_flops_l0_gemm, smumps_flops_l0_gemmt, smumps_flops_l0_trsm, smumps_flops_l0_copyscale;
#if defined (USE_XKBLAS)
/* Internal functions used to translate are defined into XKBlas library.*/
#else
/* Internal functions to translate character to cublasXt types */
static inline cublasOperation_t op_blas2cublasxt( const char* trans )
{
  switch (trans[0]) {
    case 'n':
    case 'N': return CUBLAS_OP_N;
    case 't':
    case 'T': return CUBLAS_OP_T;
    case 'c':
    case 'C': return CUBLAS_OP_C;
    default:
      return -1;
  }
}
static inline cublasSideMode_t sidemode_blas2cublasxt( const char* side )
{
  switch (side[0]) {
    case 'l':
    case 'L': return CUBLAS_SIDE_LEFT;
    case 'r':
    case 'R': return CUBLAS_SIDE_RIGHT;
    default:
      return -1;
  }
}
static inline cublasFillMode_t fillmode_blas2cublasxt( const char* uplo )
{
  switch (uplo[0]) {
    case 'l':
    case 'L': return CUBLAS_FILL_MODE_LOWER;
    case 'u':
    case 'U': return CUBLAS_FILL_MODE_UPPER;
    default:
      printf("Illegal uplo setting:: %c\n", *uplo);
      return -1;
  }
}
static inline cublasDiagType_t diagtype_blas2cublasxt( const char* diag )
{
  switch (diag[0]) {
    case 'n':
    case 'N': return CUBLAS_DIAG_NON_UNIT;
    case 'u':
    case 'U': return CUBLAS_DIAG_UNIT;
    default:
      return -1;
  }
}
#endif
static inline int isEnoughWorkGemm(int m, int n, int k)
{
  return ((2.0 * (double)m * (double)n * (double)k) / ((double)m * (double)n + (double)m * (double)k + (double)n * (double)k) >= (double)smumps_arith_threshold );
}/* isEnoughWorkGemm */
static inline int isEnoughWorkTrsm(int m, int n, char* side)
{
  double comp;
  if (*side == 'L' || *side == 'l')
    comp = (double)m * (double)n * (double)m / ((double)m * (double)m / 2.0 + 2.0*(double)m * (double)n);
  else
    comp = (double)m * (double)n * (double)n / ((double)n * (double)n / 2.0 + 2.0*(double)m * (double)n);
  return (comp >= (double)smumps_arith_threshold );
}/* isEnoughWorkTrsm */
static inline double flopsTrsm(int m, int n, char* side)
{
  double flops;
  if (*side == 'L' || *side == 'l')
    flops = (double)m * (double)n * (double)m;
  else
    flops = (double)m * (double)n * (double)n;
  return flops;
}/* flopsTrsm */
static inline int isEnoughWorkGemmT(int n, int k)
{
  return (((double)n * (double)n * (double)k) / ((double)n * (double)n + (double)n * (double)k + (double)n * (double)k )) >= (double)smumps_arith_threshold;
}/* isEnoughWorkGemmT */
static inline void add_value_to_counter(double *counter_above, double *counter_under, double value)
{
  if (smumps_gpu_under_l0 == 1)
  {
#if defined(_OPENMP)
#pragma omp atomic
#endif
    *counter_under += value;
  } else {
    /* if above L0-OMP, no needs to use an atomic. */
    *counter_above = *counter_above + value;
  }
}/* add_value_to_counter */
/* Functions */
void MUMPS_CALL SMUMPS_GET_NDEVICES_C(MUMPS_INT *NDEVICES)
{
  int error, ndevices_int;
  ndevices_int=0;
#if ( __XKBLAS__ >= 1 || ( __XKBLAS__ == 0 && __XKBLAS_MINOR__ >= 4 ) )
/* No CUDA dependency for recent XKBLAS versions */
  error=xkblas_get_device_count(&ndevices_int);
#else
  error=cudaGetDeviceCount(&ndevices_int);
#endif
  *NDEVICES=(MUMPS_INT)ndevices_int;
}
void MUMPS_CALL SMUMPS_GET_DEVICE_PROPERTIES_C(MUMPS_INT *ARCH, MUMPS_INT *HAS_UNIFIED)
{
/* On exit: *ARCH provides the class of architecture, e.g. A100-V100 / gracehopper / MI300 / unknown
            *HAS_UNIFIED indicates if unified memory is supported by software (bit 0) and/or by hardware (bit 1) */
#if defined(USE_XKBLAS)
#if ( __XKBLAS__ >= 1 || ( __XKBLAS__ == 0 && __XKBLAS_MINOR__ >= 6 ) )
  struct xkblas_dev_prop_t prop;
#endif
#endif
  /* Initialize to values if xkblas_get_device_properties not supported:*/
  *ARCH=-65536;
  *HAS_UNIFIED =-65536; /* this value has bits 0 and 1 equal to 0 */
#if defined(USE_XKBLAS)
#if ( __XKBLAS__ >= 1 || ( __XKBLAS__ == 0 && __XKBLAS_MINOR__ >= 6 ) )
  xkblas_get_device_properties(&prop);
  *ARCH    = prop.arch;
   /* prop.arch = xxx
    *           = yyy
    */
  *HAS_UNIFIED=prop.has_unified;
   /* bit 0: 0 = no software support, 1 = software support
      bit 1: 0 = no hardware support, 1 = hardware support
      Example: value = 3 = 011b on GraceHopper.
               value = 1 = 001b on A100, V100. */
#endif
#endif
  return;
}
#if defined(USE_XKBLAS)
/* SMUMPS_INIT_GPU_C: XKBLAS  version. */
void MUMPS_CALL
SMUMPS_INIT_GPU_C(MUMPS_INT *K415, MUMPS_INT *K417, MUMPS_INT *K418,
                  MUMPS_INT *K419, MUMPS_INT *K422, MUMPS_INT *K430,
                  MUMPS_INT *K436,
                  MUMPS_INT *GPUBITMAP, MUMPS_INT *MAX_MPIs_PER_GPU,
                  MUMPS_INT info[2]
)
/*
 * *K430 is inout. It is not modified except for recent xkblas versions
 * where it is set to 2 when *K436=1 (unified mode). In this case,
 * xkblas_malloc_unified and xkblas_free_unified will
 * be called for most allocations/deallocations
 * concerned by GPU computations.
 *
 * *K436 should NEVER be set to 1 if *K437==0
 *
 * If *K436==1, then the environment variable XKBLAS_UNIFIED is set to 1.
 *
 * *K436 should not be set to 1 if *K437=0, this is tested in MUMPS before calling xkblas_init
 *
 *
 * *K436 should NEVER be set to 1 if *K437==0
 *
 * If *K436==1, then the environment variable XKBLAS_UNIFIED is set to 1.
 *
 * *K436 should not be set to 1 if *K437=0, this is tested in MUMPS before calling xkblas_init
 *
 */
{
  char param[32];
  int ret;
  int tile_size = 1024;
#if defined(MUMPS_WIN32) || defined(__MINGW32__)
  int ierr;
#endif
  smumps_wtime_pin = 0.0;
  smumps_wtime_async_pin = 0.0;
  smumps_wtime_gemm = 0.0;
  smumps_wtime_trsm = 0.0;
  smumps_wtime_copyscale = 0.0;
  smumps_wtime_gemmt = 0.0;
  smumps_wtime_unpin = 0.0;
  smumps_flops_gemm = 0.0;
  smumps_flops_gemmt = 0.0;
  smumps_flops_trsm = 0.0;
  smumps_flops_copyscale = 0.0;
  smumps_size_pin = 0.0;
  smumps_size_async_pin = 0.0;
  smumps_size_async_unpin = 0.0;
  /* for L0-OMP: */
  smumps_wtime_l0_gemm = 0.0;
  smumps_wtime_l0_trsm = 0.0;
  smumps_wtime_l0_copyscale = 0.0;
  smumps_wtime_l0_gemmt = 0.0;
  smumps_flops_l0_gemm = 0.0;
  smumps_flops_l0_gemmt = 0.0;
  smumps_flops_l0_trsm = 0.0;
  smumps_flops_l0_copyscale = 0.0;
  /* parameters preparation */
#if defined(MUMPS_WIN32) || defined(__MINGW32__)
    sprintf(param, "XKBLAS_GPUSET=%d",*GPUBITMAP);
    ierr = _putenv(param);
#else
  sprintf(param, "%d", *GPUBITMAP);
  setenv("XKBLAS_GPUSET",param,1);
#endif
  if (*K417 > 0) /* Note: always the case when SMUMPS_GPU_INIT_C is called */
  {
#if defined(MUMPS_WIN32) || defined(__MINGW32__)
    sprintf(param, "XKBLAS_NGPUS=%d",*K417);
    ierr = _putenv(param);
#else
    sprintf(param, "%d", *K417);
    setenv("XKBLAS_NGPUS",param,1);
#endif
#if ( __XKBLAS__ >= 1 || ( __XKBLAS__ == 0 && __XKBLAS_MINOR__ >= 6 ) )
    if (*K436==1){
      /* Allocation mode: say if xkblas_malloc_unified should be used for S workspace */
#if defined(MUMPS_WIN32) || defined(__MINGW32__)
      sprintf(param, "XKBLAS_UNIFIED=%d",1);
#else
      sprintf(param, "%d", 1);
#endif
      *K430 = 2;
    }
    else {
#if defined(MUMPS_WIN32) || defined(__MINGW32__)
      sprintf(param, "XKBLAS_UNIFIED=%d",0);
#else
      sprintf(param, "%d", 0);
#endif
    }
#if defined(MUMPS_WIN32) || defined(__MINGW32__)
    ierr = _putenv(param);
#else
    setenv("XKBLAS_UNIFIED",param,1);
#endif
#endif
  }
#if defined(MUMPS_WIN32) || defined(__MINGW32__)
  sprintf(param, "XKBLAS_CACHE_LIMIT=%d",90);
#else
  sprintf(param, "%d", 90);
#endif
  if (*MAX_MPIs_PER_GPU > 1)
  {
#if defined(MUMPS_WIN32) || defined(__MINGW32__)
    sprintf(param, "XKBLAS_CACHE_LIMIT=%d",90/ *MAX_MPIs_PER_GPU);
#else
    sprintf(param, "%d", 90 / *MAX_MPIs_PER_GPU);
#endif
  }
#if defined(MUMPS_WIN32) || defined(__MINGW32__)
  ierr = _putenv(param);
#else
  setenv("XKBLAS_CACHE_LIMIT",param,1);
#endif
  if (*K419 > 0)
  { 
    tile_size = *K419;
  }
  xkblas_set_param( tile_size, sizeof(SMUMPS_COMPLEX) );
  ret = xkblas_init();
  if (ret != 0)
  {
    info[0]=-105;
    info[1]=ret;
    return;
  }
  smumps_arith_threshold = 0;
    /* Threshold obtained from *K418 */
    if ( *K418 > 0)
    {
      smumps_arith_threshold = *K418;
    }
  smumps_arith_threshold_above = smumps_arith_threshold;
  smumps_arith_threshold_under = *K415;
  smumps_gpu_under_l0 = 0; /* L0OMP status */
  return;
}/* SMUMPS_INT_GPU_C (using XKblas) */
#else
/* SMUMPS_INIT_GPU_C: CUBLASXT version.
   *K430 and *K436 are ignored in this version */
void MUMPS_CALL
SMUMPS_INIT_GPU_C(MUMPS_INT *K415, MUMPS_INT *K417, MUMPS_INT *K418,
                  MUMPS_INT *K419, MUMPS_INT *K422, MUMPS_INT *K430,
                  MUMPS_INT *K436,
                  MUMPS_INT *GPUBITMAP, MUMPS_INT *MAX_MPIs_PER_GPU,
                  MUMPS_INT info[2]
)
{
  int i, j, ret;
/* the number of devices per MPI is bounded by 32 in MUMPS i
   on might want otherwise to allocate devices of size *K417 */
  int devices[32];
  smumps_wtime_pin = 0.0;
  smumps_wtime_async_pin = 0.0;
  smumps_wtime_gemm = 0.0;
  smumps_wtime_trsm = 0.0;
  smumps_wtime_gemmt = 0.0;
  smumps_wtime_unpin = 0.0;
  smumps_flops_gemm = 0.0;
  smumps_flops_gemmt = 0.0;
  smumps_flops_trsm = 0.0;
  smumps_size_pin = 0.0;
  smumps_size_async_pin = 0.0;
  smumps_size_async_unpin = 0.0;
  /* for L0-OMP: */
  smumps_wtime_l0_gemm = 0.0;
  smumps_wtime_l0_trsm = 0.0;
  smumps_wtime_l0_gemmt = 0.0;
  smumps_flops_l0_gemm = 0.0;
  smumps_flops_l0_gemmt = 0.0;
  smumps_flops_l0_trsm = 0.0;
  ret=cublasXtCreate(&smumps_handle);
  if (ret != CUBLAS_STATUS_SUCCESS)
  {
    info[0]=-102;
    info[1]=ret;
    return;
  }
  j=0;
  for (i=0; (i<32) &&  (j < *K417); ++i)
  {
    if ( *GPUBITMAP & (1 << i) )
      {
        devices[j] = i;
        j++;
      }
  }
  if (j != *K417)
  {
    info[0]=-799;
    info[1]=*GPUBITMAP; /* list of devices in *K416 */
    /* info[1]=j; */ /* number of devices in *K416 */
    return;
  }
  ret= cublasXtDeviceSelect(smumps_handle, *K417, devices);
  if(ret != CUBLAS_STATUS_SUCCESS)
  {
    info[0]=-103;
    info[1]=ret;
    return;
  }
  if ( *K419 > 0) 
  {
    ret=cublasXtSetBlockDim(smumps_handle, *K419);
    if (ret != CUBLAS_STATUS_SUCCESS)
    {
      info[0]=-104;
      info[1]=ret;
      return ;
    }
  }
  smumps_arith_threshold = 0;
    /* Threshold obtained from *K418 */
    if ( *K418 > 0) 
    {
      smumps_arith_threshold = *K418;
    }
  smumps_arith_threshold_above = smumps_arith_threshold;
  smumps_arith_threshold_under = *K415; /* FIXME: use a KEEP */
  smumps_gpu_under_l0 = 0; /* L0OMP status */
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
  if ( *K422 == 2)
  {
    mumps_init_request_pool(1);
  }
#endif
  return;
}/* SMUMPS_INIT_GPU_C (using cublasXt) */
#endif
MUMPS_INT MUMPS_CALL
SMUMPS_END_GPU_C( MUMPS_INT *K400, MUMPS_INT *K422
                  , SMUMPS_REAL *flops_gpu
                  , SMUMPS_REAL *time_gpu
                )
{
  *flops_gpu=smumps_flops_gemm+smumps_flops_gemmt+smumps_flops_trsm+smumps_flops_copyscale;
  *time_gpu=smumps_wtime_gemm+smumps_wtime_gemmt+smumps_wtime_trsm+smumps_wtime_copyscale;
  if (*K400 > 0)
  {
    /* Take into account L0-threads; time is averaged on all threads */
    *flops_gpu=*flops_gpu+smumps_flops_l0_gemm+smumps_flops_l0_gemmt+smumps_flops_l0_trsm+smumps_flops_l0_copyscale;
    *time_gpu=*time_gpu+(smumps_wtime_l0_gemm+smumps_wtime_l0_gemmt+smumps_wtime_l0_trsm+smumps_wtime_l0_trsm) / *K400;
  }
#if defined(USE_XKBLAS)
  xkblas_finalize();
#else
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
  if ( *K422 == 2)
  {
    mumps_finalize_request_pool();
  }
#endif
  if(cublasXtDestroy(smumps_handle) != CUBLAS_STATUS_SUCCESS)
  {
    printf("In MUMPS_END_GPU_C, cublasXtDestroy failed\n");
    return -1;
  }
#endif
  return 0;
}
MUMPS_INT MUMPS_CALL
SMUMPS_GEMM_ON_GPU_C(MUMPS_INT *GPU_ACTIVATION, MUMPS_INT *M, MUMPS_INT *N, MUMPS_INT *K)
{
  int isOnGPU = 0;
  if (*GPU_ACTIVATION > 0 && *M > 0 && *N > 0 && *K > 0)
  {
    if (isEnoughWorkGemm(*M, *N, *K))
    {
      isOnGPU = 1;
    }
  }
  return isOnGPU;
}/* smumps_gemm_on_gpu_C */
MUMPS_INT MUMPS_CALL
SMUMPS_GEMMT_ON_GPU_C(MUMPS_INT *GPU_ACTIVATION, MUMPS_INT *N, MUMPS_INT *K)
{
  int gemmt_on_gpu = 0;
  if (*GPU_ACTIVATION > 0 && *N > 0 && *K > 0)
  {
    if (isEnoughWorkGemmT(*N, *K))
    {
      gemmt_on_gpu = 1;
    }
  }
  return gemmt_on_gpu;
}/* smumps_gemmt_on_gpu_c */
#if ! defined(GEMMT_AVAILABLE_FOR_GPU) && ! defined(USE_XKBLAS)
MUMPS_INT MUMPS_CALL
SMUMPS_SYMSCHUR_ON_GPU_C(MUMPS_INT *GPU_ACTIVATION, MUMPS_INT *N, MUMPS_INT *K, MUMPS_INT *BLOCKING)
{
  int symschur_on_gpu = 0;
  int blsize;
  if (*GPU_ACTIVATION > 0 && *N > 0 && *K > 0)
  {
    double split = 5.0;
    double d_n = *N;
    int max_blocking = *BLOCKING * 2;
    /* round above to avoid one tiny block */
    blsize = (*N+split-1)/split;
    if (blsize > max_blocking)
    {
      /* split more */
      split = ( *N + max_blocking - 1 )/ max_blocking;
      blsize = (*N+split-1)/split;
    }
    double d_k = *K;
    double arith_intensity2;
    arith_intensity2 = (2.0*d_n*blsize*d_k)/(d_n*blsize+d_n*d_k+blsize*d_k);
    if (arith_intensity2 >= smumps_arith_threshold)
    {
      symschur_on_gpu = 1;
      *BLOCKING = blsize;
    }
  }
  return symschur_on_gpu;
}/* smumps_need_pinning_C */
#endif
#if defined(USE_XKBLAS)
void MUMPS_CALL
SMUMPS_GEMM_GPU_C(MUMPS_INT *gpu_async, char* TransA, char* TransB, MUMPS_INT *M, MUMPS_INT *N, MUMPS_INT *K, SMUMPS_COMPLEX *ALPHA, SMUMPS_COMPLEX *A,
                  MUMPS_INT *LDA, SMUMPS_COMPLEX *B, MUMPS_INT *LDB, SMUMPS_COMPLEX *BETA, SMUMPS_COMPLEX *C, MUMPS_INT *LDC, mumps_ftnlen l1)
{
  double start, elapsed;
  double flops;
  start = MPI_Wtime();
  xkblas_sgemm_async(xkblas_blas2cblas_trans(TransA), xkblas_blas2cblas_trans(TransB), *M, *N, *K,
      ALPHA, A, *LDA,
             B, *LDB,
      BETA,  C, *LDC
  );
  xkblas_memory_coherent_async(0,0,*M, *N, C, *LDC, sizeof(SMUMPS_COMPLEX));
  if (*gpu_async == 0)
    {
      xkblas_sync();
      xkblas_memory_invalidate_caches();
    }
  else if (*gpu_async ==2)
    {
    xkblas_sync();
    /* do not invalidate caches */
    }
  elapsed = MPI_Wtime() - start;
  flops = 2*((double)*M*((double)*N)*((double)*K));
  add_value_to_counter(&smumps_wtime_gemm, &smumps_wtime_l0_gemm, elapsed);
  add_value_to_counter(&smumps_flops_gemm, &smumps_flops_l0_gemm, flops);
}/* smumps_gemm_gpu_C (with XKblas) */
#else
void MUMPS_CALL
SMUMPS_GEMM_GPU_C(MUMPS_INT *gpu_async, char* TransA, char* TransB, MUMPS_INT *M, MUMPS_INT *N, MUMPS_INT *K, SMUMPS_COMPLEX *ALPHA, SMUMPS_COMPLEX *A, 
                  MUMPS_INT *LDA, SMUMPS_COMPLEX *B, MUMPS_INT *LDB, SMUMPS_COMPLEX *BETA, SMUMPS_COMPLEX *C, MUMPS_INT *LDC, mumps_ftnlen l1)
{
  cublasStatus_t cublasStatus;
  double start, elapsed;
  double flops;
  start = MPI_Wtime();
  cublasStatus=cublasXtSgemm(smumps_handle, op_blas2cublasxt(TransA),  op_blas2cublasxt(TransB), *M,*N, *K, ALPHA, A, *LDA, B, *LDB, BETA, C, *LDC);
    /* FIXME: decide if the following error notification should be "public" or what to do with it */
  if(cublasStatus != CUBLAS_STATUS_SUCCESS)
  {
    printf("In SMUMPS_GEMM_GPU_C, cublasXtSgemm failed with status %i\n",cublasStatus);
  }
  elapsed = MPI_Wtime() - start;
  flops = 2*((double)*M*((double)*N)*((double)*K));
  add_value_to_counter(&smumps_wtime_gemm, &smumps_wtime_l0_gemm, elapsed);
  add_value_to_counter(&smumps_flops_gemm, &smumps_flops_l0_gemm, flops);
}/* smumps_gemm_gpu_C (with cublasXt) */
#endif
#if defined(USE_XKBLAS)||defined(GEMMT_AVAILABLE_FOR_GPU)
void MUMPS_CALL
SMUMPS_GEMMT_GPU_C(MUMPS_INT *gpu_async, char* Uplo, char* TransA, char* TransB, MUMPS_INT *N, MUMPS_INT *K, SMUMPS_COMPLEX *ALPHA, SMUMPS_COMPLEX *A, MUMPS_INT *LDA, SMUMPS_COMPLEX *B, MUMPS_INT *LDB, SMUMPS_COMPLEX *BETA, SMUMPS_COMPLEX *C, MUMPS_INT *LDC, mumps_ftnlen l1)
{
  double start, elapsed;
  double flops;
  start = MPI_Wtime();
  xkblas_sgemmt_async(xkblas_blas2cblas_fill(Uplo), xkblas_blas2cblas_trans(TransA), xkblas_blas2cblas_trans(TransB), *N, *K,
    ALPHA, A, *LDA,
           B, *LDB,
    BETA,  C, *LDC
  );
  xkblas_memory_coherent_async(xkblas_blas2cblas_fill(Uplo), 0, *N, *N, C, *LDC, sizeof(SMUMPS_COMPLEX));
  xkblas_sync();
  xkblas_memory_invalidate_caches();
  elapsed = MPI_Wtime() - start;
  flops = ((double)*N*((double)*N)*((double)*K));
  add_value_to_counter(&smumps_wtime_gemmt, &smumps_wtime_l0_gemmt, elapsed);
  add_value_to_counter(&smumps_flops_gemmt, &smumps_flops_l0_gemmt, flops);
}/* smumps_gemmt_gpu_C (with XKblas) */
#else
void MUMPS_CALL
SMUMPS_SYRKX_GPU_C(MUMPS_INT *gpu_async, char* Uplo, char* Trans, MUMPS_INT *N, MUMPS_INT *K, SMUMPS_COMPLEX *ALPHA, SMUMPS_COMPLEX *A, MUMPS_INT *LDA, SMUMPS_COMPLEX *B, MUMPS_INT *LDB, SMUMPS_COMPLEX *BETA, SMUMPS_COMPLEX *C, MUMPS_INT *LDC, mumps_ftnlen l1)
{
  cublasStatus_t cublasStatus;
  double start, elapsed;
  double flops;
  start = MPI_Wtime();
  cublasStatus=cublasXtSsyrkx(smumps_handle, fillmode_blas2cublasxt(Uplo), op_blas2cublasxt(Trans), *N, *K, ALPHA, A, *LDA, B, *LDB, BETA, C, *LDC);
  /* flops counted as gemmt (syrkx is a particular case of gemmt kernel) */
  elapsed = MPI_Wtime() - start;
  flops = ((double)*N*((double)*N)*((double)*K));
  add_value_to_counter(&smumps_wtime_gemmt, &smumps_wtime_l0_gemmt, elapsed);
  add_value_to_counter(&smumps_flops_gemmt, &smumps_flops_l0_gemmt, flops);
}
#endif
MUMPS_INT MUMPS_CALL
SMUMPS_TRSM_ON_GPU_C(MUMPS_INT *GPU_ACTIVATION, MUMPS_INT *M, MUMPS_INT *N, char* Side)
{
  int isOnGPU = 0;
  if (*GPU_ACTIVATION > 0 && *M > 0 && *N > 0)
  {
    if (isEnoughWorkTrsm(*M, *N, Side))
    {
      isOnGPU = 1;
    }
  }
  return isOnGPU;
}/* smumps_trsm_on_gpu_C */
#if defined(USE_XKBLAS)
void MUMPS_CALL SMUMPS_GPU_MALLOC_C(MUMPS_INT8 *address, MUMPS_INT8 *size)
{
  void * ptr;
  int ierr;
#if ( __XKBLAS__ >= 1 || ( __XKBLAS__ == 0 && __XKBLAS_MINOR__ >= 6 ) )
  ierr=xkblas_malloc_unified(&ptr,*size);
  if (ierr != 0) ptr=0;
#else
  ptr=malloc(*size);
#endif
  *address=(MUMPS_INT8)(ptr);
}
void MUMPS_CALL SMUMPS_GPU_FREE_C(MUMPS_INT8 *address)
{
#if ( __XKBLAS__ >= 1 || ( __XKBLAS__ == 0 && __XKBLAS_MINOR__ >= 6 ) )
  xkblas_free_unified(address);
#else
  free(address);
#endif
}
#if defined(USE_XKBLAS_COPYSCALE)
void MUMPS_CALL SMUMPS_GPU_COPY2U_SCALEL_C(MUMPS_INT *gpu_async, MUMPS_INT *irowmax, MUMPS_INT *irowmin, MUMPS_INT *copysize, MUMPS_INT *nfront, MUMPS_INT *npiv, MUMPS_INT *liw, MUMPS_INT *IW, MUMPS_INT *offset_iw, MUMPS_INT8 *la, SMUMPS_COMPLEX *A, MUMPS_INT8 *dpos, MUMPS_INT8 *lpos, MUMPS_INT8 *upos, MUMPS_INT *should_copy)
/* Warning: offset_iw = 1 means no offset, upos, lpos, apos also correspond to Fortran array indexing
 * FIXME: A[lpos-1] should be
 */
{
  int M,N;
  double start, elapsed;
  double flops;
  start = MPI_Wtime();
   M = (*irowmax-*irowmin+1);
   N = *npiv; 
xkblas_scopyscale_async( (size_t)M, (size_t)N, (*should_copy), &IW[*offset_iw-1],
                         &A[*dpos-1], (size_t)(*nfront),
                         &A[*lpos-1], (size_t)(*nfront),
                         &A[*upos-1], (size_t)(*nfront)
);
  /* L will be needed on CPU */
  xkblas_memory_coherent_async(0, 0, N, M, &A[*lpos-1], *nfront, sizeof(SMUMPS_COMPLEX));
  /* U: to be checked in which cases CPU really needs it back */
  xkblas_memory_coherent_async(0, 0, M, N, &A[*upos-1], *nfront, sizeof(SMUMPS_COMPLEX));
  /* (Note: one could possibly asynchronously invalidate cache for &A[*dpos-1]
   * to free memory earlier in xkblas) */
  if (*gpu_async == 0)
  {
    xkblas_sync();
    xkblas_memory_invalidate_caches();
  }
  else if (*gpu_async == 2)
  {
    xkblas_sync();
    /* do not invalidate caches, L and U
     * will be reused on GPU for GEMMT */
  }
  elapsed = MPI_Wtime() - start;
  flops = M*N;
  add_value_to_counter(&smumps_wtime_copyscale, &smumps_wtime_l0_copyscale, elapsed);
  add_value_to_counter(&smumps_flops_copyscale, &smumps_flops_l0_copyscale, flops);
} /* SMUMPS_GPU_COPY2U_SCALEL_C */
#endif /* USE_XKBLAS_COPYSCALE */
void MUMPS_CALL
SMUMPS_TRSM_GPU_C(MUMPS_INT *gpu_async, char* Side, char* Uplo, char* TransA, char* Diag, MUMPS_INT *M, MUMPS_INT *N,
                  SMUMPS_COMPLEX *ALPHA, SMUMPS_COMPLEX *A, MUMPS_INT *LDA, SMUMPS_COMPLEX *B, MUMPS_INT *LDB, mumps_ftnlen l1)
{
  double start, elapsed;
  double flops;
  start = MPI_Wtime();
  xkblas_strsm_async(xkblas_blas2cblas_side(Side), xkblas_blas2cblas_fill(Uplo),
                     xkblas_blas2cblas_trans(TransA), xkblas_blas2cblas_diag(Diag),
                     *M, *N,
                     ALPHA, A, *LDA,
                            B, *LDB
  );
  /* gpu_async = 3: do not bring back L21 on CPU because of copy-scale */
  if (*gpu_async != 3)
  {
    xkblas_memory_coherent_async(0, 0, *M, *N, B, *LDB, sizeof(SMUMPS_COMPLEX));
    if (*gpu_async == 0)
    {
      xkblas_sync();
      xkblas_memory_invalidate_caches();
    }
    else if (*gpu_async == 2)
    {
      xkblas_sync();
      /* do not invalidate caches, resulting matrix B could be reused on GPU */
    }
  }
  elapsed = MPI_Wtime() - start;
  flops = flopsTrsm(*M, *N, Side);
  add_value_to_counter(&smumps_wtime_trsm, &smumps_wtime_l0_trsm, elapsed);
  add_value_to_counter(&smumps_flops_trsm, &smumps_flops_l0_trsm, flops);
}/* smumps_trsm_gpu_C (with XKblas) */
#else
void MUMPS_CALL
SMUMPS_TRSM_GPU_C(MUMPS_INT *gpu_async, char* Side, char* Uplo, char* TransA, char* Diag, MUMPS_INT *M, MUMPS_INT *N, SMUMPS_COMPLEX *ALPHA, SMUMPS_COMPLEX *A, 
                  MUMPS_INT *LDA, SMUMPS_COMPLEX *B, MUMPS_INT *LDB, mumps_ftnlen l1)
{
  double start, elapsed;
  double flops;
  start = MPI_Wtime();
  cublasStatus_t cublasStatus;
  cublasStatus=cublasXtStrsm(smumps_handle, sidemode_blas2cublasxt(Side), fillmode_blas2cublasxt(Uplo), op_blas2cublasxt(TransA), diagtype_blas2cublasxt(Diag), *M, *N, ALPHA, A, *LDA, B, *LDB);
  /* FIXME: decide if the following error notification should be "public" or what to do with it */
  if(cublasStatus != CUBLAS_STATUS_SUCCESS)
  {
    printf("In SMUMPS_TRSM_GPU_C, cublasXtStrsm failed with status %i\n",cublasStatus);
  }
  elapsed = MPI_Wtime() - start;
  flops = flopsTrsm(*M, *N, Side);
  add_value_to_counter(&smumps_wtime_trsm, &smumps_wtime_l0_trsm, elapsed);
  add_value_to_counter(&smumps_flops_trsm, &smumps_flops_l0_trsm, flops);
}/* smumps_trsm_gpu_C (with cublasXt) */
#endif
/* With respect to asynchronous pinning/unpinning, it is applied to CublasXt
 * using mumps_register_thread module while if using XKBLAS,
 * asynchronous pinning and unpinning use XKblas API. */
#if defined(USE_XKBLAS)
void MUMPS_CALL
SMUMPS_GPU_ASYNC_PIN_MEMORY_C(SMUMPS_COMPLEX *A, MUMPS_INT8 *ADIM)
{
  uint64_t handle;
  double flops;
  smumps_start_async_pin = MPI_Wtime();
#if defined(DEBUG_GPU)
  printf(" <= GPU_ASYNC_PIN_MEMORY_C: %p of size %li\n",A, *ADIM*sizeof(SMUMPS_COMPLEX));
#endif
  handle = xkblas_register_memory_async( A, *ADIM*sizeof(SMUMPS_COMPLEX) );
  flops = *ADIM*sizeof(SMUMPS_COMPLEX);
#if defined(_OPENMP)
#pragma omp atomic
#endif
  smumps_size_async_pin = smumps_size_async_pin + flops;
#if defined(DEBUG_GPU)
  printf(" <= GPU_ASYNC_PIN_MEMORY_C: %p of size %li, handle= %li\n",A, *ADIM*sizeof(SMUMPS_COMPLEX), handle);
#endif
}/* smumps_gpu_async_pin_memory_C */
void MUMPS_CALL
SMUMPS_GPU_WAIT_ENDOFPINNING_C()
{
  double start, stop, elapsed;
  start = MPI_Wtime();
#if defined(DEBUG_GPU)
  printf(" => GPU_WAIT_ENDOFPINNING_C\n");
#endif
  xkblas_register_memory_waitall();
  stop = MPI_Wtime();
  elapsed = stop - start;
#if defined(_OPENMP)
#pragma omp atomic
#endif
  smumps_wtime_async_pin = smumps_wtime_async_pin + elapsed;
  elapsed = stop - smumps_start_async_pin;
#if defined(_OPENMP)
#pragma omp atomic
#endif
  smumps_wtime_pin = smumps_wtime_pin + elapsed;
#if defined(DEBUG_GPU)
  printf(" <= GPU_WAIT_ENDOFPINNING_C=%lf \n",elapsed);
#endif
}/* smumps_gpu_wait_endofpinning_C */
void MUMPS_CALL
SMUMPS_GPU_ASYNC_UNPIN_MEMORY_C(SMUMPS_COMPLEX *A, MUMPS_INT8 *ADIM)
{
  uint64_t handle;
#if defined(DEBUG_GPU)
  printf(" => GPU_ASYNC_UNPIN_MEMORY_C: %p of size %li \n",A, *ADIM*sizeof(SMUMPS_COMPLEX));
#endif
  handle = xkblas_unregister_memory_async((void *)A, *ADIM*sizeof(SMUMPS_COMPLEX));
  double flops = *ADIM*sizeof(SMUMPS_COMPLEX);
#if defined(_OPENMP)
#pragma omp atomic
#endif
  smumps_size_async_unpin = smumps_size_async_unpin + flops;
#if defined(DEBUG_GPU)
  printf(" <= GPU_ASYNC_UNPIN_MEMORY_C: %p of size %li, handle= %li\n",A, *ADIM*sizeof(SMUMPS_COMPLEX), handle);
#endif
}/* smumps_gpu_async_unpin_memory_c */
#else
void MUMPS_CALL
SMUMPS_GPU_ASYNC_PIN_MEMORY_C(SMUMPS_COMPLEX *A, MUMPS_INT8 *ADIM)
{
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
  smumps_start_async_pin = MPI_Wtime();
  uint64_t handle;
#if defined(DEBUG_GPU)
  printf(" => GPU_ASYNC_PIN_MEMORY_C: %p of size %li \n",A, *ADIM*sizeof(SMUMPS_COMPLEX));
#endif
  handle = mumps_host_register_async((void *)A, *ADIM*sizeof(SMUMPS_COMPLEX));
  double flops = *ADIM*sizeof(SMUMPS_COMPLEX);
#if defined(_OPENMP)
#pragma omp atomic
#endif
  smumps_size_async_pin = smumps_size_async_pin + flops;
#endif
#if defined(DEBUG_GPU)
  printf(" <= GPU_ASYNC_PIN_MEMORY_C: %p of size %li, handle= %li\n",A, *ADIM*sizeof(SMUMPS_COMPLEX), handle);
#endif
}/* smumps_gpu_async_pin_memory_c */
void MUMPS_CALL
SMUMPS_GPU_ASYNC_UNPIN_MEMORY_C(SMUMPS_COMPLEX *A, MUMPS_INT8 *ADIM)
{
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
  uint64_t handle;
#if defined(DEBUG_GPU)
  printf(" => GPU_ASYNC_UNPIN_MEMORY_C: %p of size %li \n",A, *ADIM*sizeof(SMUMPS_COMPLEX));
#endif
  handle = mumps_host_unregister_async((void *)A, *ADIM*sizeof(SMUMPS_COMPLEX));
  double flops = *ADIM*sizeof(SMUMPS_COMPLEX);
#if defined(_OPENMP)
#pragma omp atomic
#endif
  smumps_size_async_unpin = smumps_size_async_unpin + flops;
#endif
#if defined(DEBUG_GPU)
  printf(" <= GPU_ASYNC_UNPIN_MEMORY_C: %p of size %li, handle= %li\n",A, *ADIM*sizeof(SMUMPS_COMPLEX), handle);
#endif
}/* smumps_gpu_async_unpin_memory_c */
void MUMPS_CALL
SMUMPS_GPU_WAIT_ENDOFPINNING_C()
{
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
  double start, stop, elapsed;
  start = MPI_Wtime();
#if defined(DEBUG_GPU)
  printf(" => GPU_WAIT_ENDOFPINNING_C\n");
#endif
  mumps_request_wait_all();
  stop = MPI_Wtime();
  elapsed = stop - start;
#if defined(_OPENMP)
#pragma omp atomic
#endif
  smumps_wtime_async_pin = smumps_wtime_async_pin + elapsed;
  elapsed = stop - smumps_start_async_pin;
#if defined(_OPENMP)
#pragma omp atomic
#endif
  smumps_wtime_pin = smumps_wtime_pin + elapsed;
#endif
#if defined(DEBUG_GPU)
  printf(" <= GPU_WAIT_ENDOFPINNING_C\n");
#endif
}/* smumps_gpu_wait_endofpinning_C */
#endif
void MUMPS_CALL
SMUMPS_GPU_PIN_MEMORY_C(SMUMPS_COMPLEX *A, MUMPS_INT8 *ADIM)
{
#if defined(DEBUG_GPU)
  printf("=> GPU_PIN_MEMORY_C of %p of size %li\n",A, *ADIM*sizeof(SMUMPS_COMPLEX));
#endif
  int error;
  double start, elapsed;
  double flops;
  start = MPI_Wtime();
#if defined(USE_XKBLAS)
  error = xkblas_register_memory( A, *ADIM*sizeof(SMUMPS_COMPLEX) );
  if (error != 0)
  {
    /* FIXME: error management */
    printf("error registering matrix of size %li with code %i\n", *ADIM, error);
  }
#else
  error = cudaHostRegister(A, *ADIM*sizeof(SMUMPS_COMPLEX), cudaHostRegisterPortable);
  if (error != cudaSuccess)
  {
    /* FIXME: error management */
    printf("error registering matrix of size %li with code %i\n", *ADIM, error);
  }
#endif
  elapsed = MPI_Wtime() - start;
#if defined(_OPENMP)
#pragma omp atomic
#endif
  smumps_wtime_pin = smumps_wtime_pin + elapsed;
  flops = *ADIM*sizeof(SMUMPS_COMPLEX);
#if defined(_OPENMP)
#pragma omp atomic
#endif
  smumps_size_pin = smumps_size_pin + flops;
#if defined(DEBUG_GPU)
  printf("<= GPU_PIN_MEMORY_C of %p of size %li, error =%i\n",A, *ADIM*sizeof(SMUMPS_COMPLEX), error);
#endif
}/* smumps_gpu_pin_memory_c */
void MUMPS_CALL
SMUMPS_GPU_MODIFY_ARITH_TRHRESHOLD (MUMPS_INT *OLD, MUMPS_INT *NEW)
{
  /* used to temporarily modify and restore the ARITH_TRHRESHOLD 
     in: NEW value of smumps_arith_threshold
     out: OLD current value of smumps_arith_threshold
  */
  *OLD = smumps_arith_threshold;
  smumps_arith_threshold = *NEW;
}
void MUMPS_CALL
SMUMPS_GPU_SET_UNDER_L0_STATUS_C(MUMPS_INT *VALEUR)
{
  /* If smumps_gpu_under_l0 is 1 it means we are under L0-OMP.
   * If it is 0 it means we are above L0-OMP or that there is no L0-OMP layer.
   */
  smumps_gpu_under_l0 = *VALEUR;
  if (smumps_gpu_under_l0 == 0)
  {
    smumps_arith_threshold = smumps_arith_threshold_above;
  }
  else
  {
    smumps_arith_threshold = smumps_arith_threshold_under;
  }
}
void MUMPS_CALL
SMUMPS_GPU_UNPIN_MEMORY_C(SMUMPS_COMPLEX *A, MUMPS_INT8 *ADIM)
{
#if defined(DEBUG_GPU)
  printf("=> GPU_UNPIN_MEMORY_C of %p of size %li\n",A, *ADIM*sizeof(SMUMPS_COMPLEX));
#endif
  int error;
  double start, elapsed;
  start = MPI_Wtime();
#if defined(USE_XKBLAS)
  error = xkblas_unregister_memory(A, *ADIM*sizeof(SMUMPS_COMPLEX));
  if (error != 0)
  {
    printf("error unregistering matrix with code %i\n", error);
  }
#else
  error = cudaHostUnregister(A);
  if (error != cudaSuccess)
  {
    printf("error unregistering matrix with code %i\n", error);
  }
#endif
  elapsed = MPI_Wtime() - start;
#if defined(_OPENMP)
#pragma omp atomic
  smumps_wtime_unpin = smumps_wtime_unpin + elapsed;
#endif
#if defined(DEBUG_GPU)
  printf("<= GPU_UNPIN_MEMORY_C of %p of size %li, error=%d\n",A, *ADIM*sizeof(SMUMPS_COMPLEX),error);
#endif
}/* smumps_gpu_unpin_memory_c */
#else
void MUMPS_CALL smumps_gpu_return()
{
  /*
    GPU feature will be available in the future
  */
}
#endif
