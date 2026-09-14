C
C  This file is part of MUMPS 5.9.1, released
C  on Mon Jul 20 09:00:43 UTC 2026
C
C
C  Copyright 1991-2026 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
C  Mumps Technologies, University of Bordeaux.
C
C  This version of MUMPS is provided to you free of charge. It is
C  released under the CeCILL-C license 
C  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
C  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
C
      SUBROUTINE MUMPS_PRINT_IF_DEFINED(MPG)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MPG
      IF (MPG.LE.0) RETURN
      write(MPG,*) "================================================="
#if defined(MPI_TO_K_OMP)
#if defined(MTKO_VERBOSE)
      write(MPG, *) 
     & "MUMPS compiled with option -DMTKO_VERBOSE"
#endif
#endif
#if defined(AFFINITY_VERBOSE)
      write(MPG, *) "MUMPS compiled with option -DAFFINITY_VERBOSE"
#endif
#if defined(AVOID_MPI_IN_PLACE)
      WRITE(MPG,*)
     &  " MUMPS compiled with option -DAVOID_MPI_IN_PLACE"
#endif
#if defined(NO_64BIT_ADDRESSES_IN_C)
      WRITE(MPG,*)
     &  " MUMPS compiled with option -DNO_64BIT_ADDRESSES_IN_C"
#endif
#if defined(LARGEMATRICES)
      write(MPG,*) "MUMPS compiled with option -DLARGEMATRICES"
#endif
#if defined(GEMMT_AVAILABLE)
      write(MPG, *) "MUMPS compiled with option -DGEMMT_AVAILABLE"
#endif
#if defined(DETERMINISTIC_PARALLEL_GRAPH)
      write(MPG,*) "MUMPS compiled with option"
     &     ," -DDETERMINISTIC_PARALLEL_GRAPH"
#endif
#if defined(metis)
      write(MPG,*) "MUMPS compiled with option -Dmetis"
#endif
#if defined(metis4)
      write(MPG,*) "MUMPS compiled with option -Dmetis4"
#endif
#if defined(MPI_TO_K_OMP)
      write(MPG,*) "MUMPS version includes compiled code ",
     &             "for MPI_TO_K_OMP"
#  if defined(MTKO_NO_MASK_FUSION)
      write(MPG, *) "MUMPS compiled with option -DMTKO_NO_MASK_FUSION"
#  endif
#  if defined(USE_LIBHWLOC)
      write(MPG, *) "MUMPS compiled with option -DUSE_LIBHWLOC"
#  endif
#  if defined(USE_SCHEDAFFINITY)
      write(MPG, *) "MUMPS compiled with option -DUSE_SCHEDAFFINITY"
#  endif
#  if defined(MTKO_ADVANCED)
      write(MPG, *) "MUMPS compiled with option -DMTKO_ADVANCED"
#    if !defined(USE_LIBHWLOC)&&!defined(USE_SCHEDAFFINITY)
      write(MPG, *) "MUMPS should not be compiled with -DMTKO_ADVANCED",
     &      " without either -DUSE_LIBHWLOC or -DUSE_SCHEDAFFINITY"
#    endif
#  endif
#  if defined(MTKO_ADVANCED_VERBOSE)
#  if ! defined(MTKO_ADVANCED)
      write(MPG, *) "MUMPS should not be compiled with ",
     &              "-DMTKO_ADVANCED_VERBOSE without -DMTKO_ADVANCED"
      CALL MUMPS_ABORT()
#  else
      write(MPG, *) "MUMPS compiled with option -DMTKO_ADVANCED_VERBOSE"
#  endif
#  endif
#endif
#if defined(MUMPS_ALLOC_FROM_C)
      WRITE(MPG,*) "MUMPS compiled with option -DMUMPS_ALLOC_FROM_C"
#endif
#if defined(MUMPS_NOF2003)
      write(MPG,*) "MUMPS compiled with option -DMUMPS_NOF2003"
#endif
#if defined(MUMPS_SCOTCHIMPORTOMPTHREADS)
      write(MPG,*)
     &"MUMPS compiled with option -DMUMPS_SCOTCHIMPORTTHREADS"
#endif
#if defined(MUMPS_WINLARGEFILES)
      write(MPG,*) "MUMPS compiled with option -DMUMPS_WINLARGEFILES"
#endif
#if defined(MUMPS_USE_BLAS2)
      write(MPG,*) "MUMPS compiled with option -DMUMPS_USE_BLAS2"
#endif
#if defined(NOSCALAPACK)
      write(MPG,*) "MUMPS compiled with option -DNOSCALAPACK"
#endif
#if defined(OLD_OOC_NOPANEL)
      write(MPG,*) "MUMPS compiled with option -DOLD_OOC_NOPANEL"
#endif
#if defined(parmetis)
      write(MPG,*) "MUMPS compiled with option -Dparmetis"
#endif
#if defined(parmetis3)
      write(MPG,*) "MUMPS compiled with option -Dparmetis3"
#endif
#if defined(pord)
      write(MPG,*) "MUMPS compiled with option -Dpord"
#endif
#if defined(PRINT_BACKTRACE_ON_ABORT)
      write(MPG, *)
     &  "MUMPS compiled with option -DPRINT_BACKTRACE_ON_ABORT"
#endif
#if defined(ptscotch)
      write(MPG,*) "MUMPS compiled with option -Dptscotch"
#endif
#if defined(scotch)
      write(MPG,*) "MUMPS compiled with option -Dscotch"
#endif
#if defined(NOAGG1)
      write(MPG,*) "MUMPS compiled with option -DNOAGG1"
#endif
#if defined(NOAGG2)
      write(MPG,*) "MUMPS compiled with option -DNOAGG2"
#endif
#if defined(NOAGG3)
      write(MPG,*) "MUMPS compiled with option -DNOAGG3"
#endif
#if defined(NOAGG4)
      write(MPG,*) "MUMPS compiled with option -DNOAGG4"
#endif
#if defined(NOAGG5)
      write(MPG,*) "MUMPS compiled with option -DNOAGG5"
#endif
#if defined(NOAMALGTOFATHER)
      write(MPG,*) "MUMPS compiled with -DNOAMALGTOFATHER"
#endif
#if defined(NO_FDM_DESCBAND)
      write(MPG,*) "MUMPS compiled with -DNO_FDM_DESCBAND"
#endif
#if defined(NO_FDM_MAPROW)
      write(MPG,*) "MUMPS compiled with -DNO_FDM_MAPROW"
#endif
#if defined(NO_SAVE_RESTORE)
      write(MPG,*) "MUMPS compiled with -DNO_SAVE_RESTORE"
#endif
#if defined(NO_SPLIT_IN_BLRGROUPING)
      write(MPG,*) "MUMPS compiled with -DNO_SPLIT_IN_BLRGROUPING"
#endif
#if defined(NODYNAMICCB)
      write(MPG,*) "MUMPS compiled with option -DNODYNAMICCB"
#endif
#if defined(USE_GPU)
      write(MPG, *) "MUMPS compiled with option -DUSE_GPU"
#  if defined(USE_XKBLAS)
      write(MPG, *) "MUMPS compiled with option -DUSE_XKBLAS"
#  endif
#  if defined(USE_XKBLAS_COPYSCALE)
      write(MPG, *) "MUMPS compiled with option -DUSE_XKBLAS_COPYSCALE"
#  endif
#  if defined(GEMMT_AVAILABLE_FOR_GPU)
      write(MPG, *) "MUMPS compiled with option ",
     &              "-DGEMMT_AVAILABLE_FOR_GPU"
#  endif
#endif
#if defined(USE_OLD_SCALING)
      write(MPG,*) "MUMPS compiled with option -DUSE_OLD_SCALING"
#endif
#if defined(VHOFFLOAD)
      write(MPG,*) "MUMPS compiled with -DVHOFFLOAD"
#endif
#if defined(WORKAROUNDINTELILP64MPI2INTEGER)
      write(MPG,*) "MUMPS compiled with option"
     &     ," -DWORKAROUNDINTELILP64MPI2INTEGER"
#endif
#if defined(WORKAROUNDILP64MPICUSTOMREDUCE)
      write(MPG,*) "MUMPS compiled with option"
     &     ," -DWORKAROUNDILP64MPICUSTOMREDUCE"
#endif
#if defined(WORKAROUNDINTELILP64OPENMPLIMITATION)
      write(MPG,*) "MUMPS compiled with option"
     &     ," -DWORKAROUNDINTELILP64OPENMPLIMITATION"
#endif
#if defined(ZERO_TRIANGLE)
      write(MPG,*) "MUMPS compiled with option -DZERO_TRIANGLE"
#endif
#if defined(SINGLEINDOUBLE)
      write(MPG,*) "MUMPS compiled with option -DSINGLEINDOUBLE"
#endif
#if defined(BLR_NOOPENMP)
      write(MPG,*) "MUMPS compiled with option -DBLR_NOOPENMP"
#endif
      write(MPG,*) 
     &   "This MUMPS version includes code for mixed precision BLR"
#if ! defined(NOMIXHPC)
     &   ," including memory accessors during solve"
#endif
#if defined(NOMIXHPC)
      write(MPG,*) "MUMPS compiled with NOMIXHPC"
#endif
#if defined(STAT_ES_SOLVE)
      write(MPG,*) "MUMPS compiled with option -DSTAT_ES_SOLVE"
#endif
#if defined(USE_AVX512_VBMI)
      write(MPG,*) "MUMPS compiled with option -DUSE_AVX512_VBMI"
#endif
#if defined(nrm2_returns_double)
      write(MPG,*) "MUMPS compiled with option -Dnrm2_returns_double"
#endif
      write(MPG,*) "================================================="
      RETURN
      END SUBROUTINE MUMPS_PRINT_IF_DEFINED
