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
      MODULE MUMPS_LR_COMMON
      IMPLICIT NONE
      INTEGER, PARAMETER :: NPREC_STOC_MAX = 10
      INTEGER, PARAMETER :: NPREC_CALC_MAX = 3
C     Global variables
C     Number of precisions used at runtime
      INTEGER, SAVE :: NPREC_STOC
C     All the precisions that could be used within MUMPS
      INTEGER, SAVE :: PREC_DBL, PREC_56, PREC_48, PREC_40,
     &                 PREC_SGL, PREC_24, PREC_16,
C     The storage precision of off-diagonal full-rank
C     blocks when they are uniformally stored i.e.
C     QUANTIZE == .FALSE. Set through KEEP(457). 
     &                 PREC_FR, 
C     The uniform storage precision of diagonal full-rank blocks.
C     Set through KEEP(457) to CLOSE_PREC, else remains PREC_FAC
     &                 PREC_DIAG_FR,
C     PRECISIONS(CLOSE_PREC) < CNTL(7) < PRECISIONS(CLOSE_PREC+1)
     &                 CLOSE_PREC,
C     The precision at which the factorization ran
     &                 PREC_FAC
C     Whether or not the off-diagonal blocks are
C     quantized based on their norms w.r.t. the
C     global matrix.  Set by KEEP(457).
      LOGICAL, SAVE :: QUANTIZE
      DOUBLE PRECISION, SAVE :: PRECISIONS(NPREC_STOC_MAX+1)
      INTEGER, SAVE :: COSTS_PREC(NPREC_STOC_MAX)
      INTEGER, SAVE :: NPREC_CALC
      INTEGER, SAVE :: CALC2STOC(NPREC_CALC_MAX+1)
      CONTAINS
      SUBROUTINE COMPUTE_BLR_VCS(K472, IBCKSZ, MAXSIZE, NASS,
     &           NFRONT, K35)
C     Purpose:
C     compute IBCKSZ the target size for BLR grouping
C
C     Parameters
        INTEGER, INTENT(IN) ::  MAXSIZE, NASS, K472
        INTEGER, INTENT(IN) ::  NFRONT, K35
        INTEGER, INTENT(OUT) :: IBCKSZ
        IF (K472.EQ.1) THEN
#if defined(__ve__)
          IF (NASS.LE.1000) THEN
            IBCKSZ = 256
          ELSEIF (NASS.GT.1000.AND.NASS.LE.5000) THEN
            IBCKSZ = 384
          ELSEIF (NASS.GT.5000.AND.NASS.LE.10000) THEN
            IBCKSZ = 512
          ELSEIF (NASS.GT.10000.AND.NASS.LE.100000) THEN
            IBCKSZ = 512
          ELSE
            IBCKSZ = 768
          ENDIF
#else
          IF (NASS.LE.1000) THEN
            IBCKSZ = 128
          ELSEIF (NASS.GT.1000.AND.NASS.LE.5000) THEN
            IBCKSZ = 256
          ELSEIF (NASS.GT.5000.AND.NASS.LE.10000) THEN
            IBCKSZ = 384
          ELSEIF (NASS.GT.10000.AND.NASS.LE.100000) THEN
            IBCKSZ = 512
          ELSE
            IBCKSZ = 768
          ENDIF
          IF (NFRONT.GT.20*NASS.AND.NFRONT.GT.100000) THEN
            IBCKSZ = max(IBCKSZ,min(NASS,768))
          ENDIF
#endif
          IBCKSZ = min(IBCKSZ,MAXSIZE)
        ELSE
          IBCKSZ = MAXSIZE
        ENDIF
C
      RETURN
      END SUBROUTINE COMPUTE_BLR_VCS
      SUBROUTINE MUMPS_UPD_TREE(NV, NSTEPS, N, FIRST, LPTR, RPTR, F,
     &     VLIST, FILS, FRERE_STEPS, STEP, DAD_STEPS, NE_STEPS, NA, LNA,
     &     PVS, K38, STEP_SCALAPACK_ROOT)
C
C     Purpose:
C     =======
C
C     Update the tree according to order computing during
C     grouping procedure of separators.
C
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: N, NV, NSTEPS, LNA, F, VLIST(NV)
C     INTEGER, congiguous -> assumed-shape arrays: 
C               FILS(N), FRERE_STEPS(NSTEPS),
C               DAD_STEPS(NSTEPS), STEP(N), NE_STEPS(N), NA(LNA)
      INTEGER                :: FILS(:), FRERE_STEPS(:),
     &     DAD_STEPS(:), STEP(:), NE_STEPS(:), NA(:)
      INTEGER, INTENT(INOUT) :: PVS(NSTEPS), LPTR, RPTR
C     Note that K38 can be called with either KEEP(38) or KEEP(20)
      INTEGER, INTENT(INOUT) :: K38
      INTEGER, INTENT(IN)    :: STEP_SCALAPACK_ROOT
      LOGICAL :: FIRST
      INTEGER :: PV, NODE, I
      PV        = VLIST(1)
      NODE      = ABS(STEP(PV))
      PVS(NODE) = PV
C     If node is the first (oldest) son, then update
C     fils(fils(...fils(dad(node))))
      IF(FIRST) THEN
         I = DAD_STEPS(NODE)
         DO WHILE(FILS(I).GT.0)
            I = FILS(I)
         END DO
         FILS(I) = -PV
      END IF
      IF(FRERE_STEPS(NODE) .GT. 0) THEN
C     If node has a younger brother, update frere_steps(node) This
C     relies upon the assumption that the younger brother has been
C     treated before node
         FRERE_STEPS(NODE) = PVS(ABS(STEP(FRERE_STEPS(NODE))))
      ELSE IF(FRERE_STEPS(NODE) .LT. 0) THEN
C     If node is the youngest brother, update frere_steps(node) to make
C     it point to the father
         FRERE_STEPS(NODE) = -PVS(ABS(STEP(DAD_STEPS(NODE))))
      END IF
      IF(DAD_STEPS(NODE) .EQ. 0) THEN
C     If node is a root, update NA
         NA(RPTR) = PV
         RPTR     = RPTR -1
      ELSE
C     Update dad_steps(node)
         DAD_STEPS(NODE) = PVS(ABS(STEP(DAD_STEPS(NODE))))
      END IF
C     If node is a leaf, update NA
      IF(NE_STEPS(NODE) .EQ. 0) THEN
         NA(LPTR) = PV
         LPTR     = LPTR -1
      END IF
C     Update fils
      STEP(VLIST(1)) = ABS(STEP(VLIST(1)))
      IF (STEP(VLIST(1)).EQ.STEP_SCALAPACK_ROOT) THEN
       K38 = VLIST(1)
      ENDIF
      DO I=1, NV-1
         IF(STEP(VLIST(I+1)).GT.0) STEP(VLIST(I+1)) = -STEP(VLIST(I+1))
         FILS(VLIST(I)) = VLIST(I+1)
      END DO
      FILS(VLIST(NV)) = F
      RETURN
      END SUBROUTINE MUMPS_UPD_TREE
      END MODULE MUMPS_LR_COMMON
      SUBROUTINE MUMPS_LRMIXED_CAST_ZZ(A, B, nrows, ncols, ldA, 
     &                                       ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: A(ldA,ncols), B(ldB,ncols)
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          B(I,J) = A(I,J)
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_ZZ
      SUBROUTINE MUMPS_LRMIXED_CAST_DD(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: A(ldA,ncols), B(ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          B(I,J) = A(I,J)
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_DD
      SUBROUTINE MUMPS_LRMIXED_CAST_CC(A, B, nrows, ncols, ldA,
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX :: A(ldA, ncols), B(ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          B(I,J) = A(I,J)
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_CC
      SUBROUTINE MUMPS_LRMIXED_CAST_SS(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        REAL :: A(ldA, ncols), B(ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          B(I,J) = A(I,J)
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_SS
      SUBROUTINE MUMPS_LRMIXED_CAST_ZC(A, B, nrows, ncols, ldA,
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: A(ldA,ncols)
        COMPLEX :: B(ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          B(I,J) =  cmplx(A(I,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_ZC
      SUBROUTINE MUMPS_LRMIXED_CAST_CZ(A, B, nrows, ncols, ldA,
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX :: A(ldA,ncols)
        COMPLEX(kind=8) :: B(ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          B(I,J) =  cmplx(A(I,J),kind=kind(B))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_CZ
      SUBROUTINE MUMPS_LRMIXED_CAST_DS(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: A(ldA,ncols)
        REAL :: B(ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          B(I,J) =  real(A(I,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_DS
      SUBROUTINE MUMPS_LRMIXED_CAST_SD(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        REAL :: A(ldA,ncols)
        DOUBLE PRECISION :: B(ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          B(I,J) =  dble(A(I,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_SD
#if !defined(NOMIXHPC)
C{
      ! {double,float} -> {fp16,...,fp56}
      SUBROUTINE MUMPS_LRMIXED_CAST_S16(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        implicit none
        INTEGER :: nrows, ncols, ldA, ldB
        REAL :: A(ldA,ncols)
        INTEGER*1 :: B(2*ldB,ncols)
        call mumps_float_fp16s(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_S16
      SUBROUTINE MUMPS_LRMIXED_CAST_D16(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        implicit none
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: A(ldA,ncols)
        INTEGER*1 :: B(2*ldB,ncols)
        call mumps_double_fp16s(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_D16
      SUBROUTINE MUMPS_LRMIXED_CAST_S24(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        implicit none
        INTEGER :: nrows, ncols, ldA, ldB
        REAL :: A(ldA,ncols)
        INTEGER*1 :: B(3*ldB,ncols)
        call mumps_float_fp24s(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_S24
      SUBROUTINE MUMPS_LRMIXED_CAST_D24(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        implicit none
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: A(ldA,ncols)
        INTEGER*1 :: B(3*ldB,ncols)
        call mumps_double_fp24s(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_D24
      SUBROUTINE MUMPS_LRMIXED_CAST_D40(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        implicit none
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: A(ldA,ncols)
        INTEGER*1 :: B(5*ldB,ncols)
        call mumps_double_fp40(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_D40
      SUBROUTINE MUMPS_LRMIXED_CAST_D48(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        implicit none
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: A(ldA,ncols)
        INTEGER*1 :: B(6*ldB,ncols)
        call mumps_double_fp48(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_D48
      SUBROUTINE MUMPS_LRMIXED_CAST_D56(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        implicit none
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: A(ldA,ncols)
        INTEGER*1 :: B(7*ldB,ncols)
        call mumps_double_fp56(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_D56
      ! {fp16,...,fp56} -> {double,float}
      SUBROUTINE MUMPS_LRMIXED_CAST_16S(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        implicit none
        INTEGER :: nrows, ncols, ldA, ldB
        INTEGER*1 :: A(2*ldA,ncols)
        REAL :: B(ldB,ncols)
        call mumps_fp16s_float(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_16S
      SUBROUTINE MUMPS_LRMIXED_CAST_16D(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        implicit none
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: B(ldB,ncols)
        INTEGER*1 :: A(2*ldA,ncols)
        call mumps_fp16s_double(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_16D
      SUBROUTINE MUMPS_LRMIXED_CAST_24S(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        implicit none
        INTEGER :: nrows, ncols, ldA, ldB
        REAL :: B(ldB,ncols)
        INTEGER*1 :: A(3*ldA,ncols)
        call mumps_fp24s_float(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_24S
      SUBROUTINE MUMPS_LRMIXED_CAST_24D(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        implicit none
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: B(ldB,ncols)
        INTEGER*1 :: A(3*ldA,ncols)
        call mumps_fp24s_double(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_24D
      SUBROUTINE MUMPS_LRMIXED_CAST_40D(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        implicit none
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: B(ldB,ncols)
        INTEGER*1 :: A(5*ldA,ncols)
        call mumps_fp40_double(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_40D
      SUBROUTINE MUMPS_LRMIXED_CAST_48D(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        implicit none
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: B(ldB,ncols)
        INTEGER*1 :: A(6*ldA,ncols)
        call mumps_fp48_double(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_48D
      SUBROUTINE MUMPS_LRMIXED_CAST_56D(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        implicit none
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: B(ldB,ncols)
        INTEGER*1 :: A(7*ldA,ncols)
        call mumps_fp56_double(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_56D
      ! complex domain {fp16,...,fp56} -> {double,float}
      SUBROUTINE MUMPS_LRMIXED_CAST_16C(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        implicit none
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX :: B(ldB,ncols)
        INTEGER*1 :: A(2*2*ldA,ncols)
        call mumps_fp16c_cfloat(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_16C
      SUBROUTINE MUMPS_LRMIXED_CAST_16Z(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: B(ldB,ncols)
        INTEGER*1 :: A(2*2*ldA,ncols)
        call mumps_fp16c_zdouble(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_16Z
      SUBROUTINE MUMPS_LRMIXED_CAST_24C(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX :: B(ldB,ncols)
        INTEGER*1 :: A(2*3*ldA,ncols)
        call mumps_fp24c_cfloat(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_24C
      SUBROUTINE MUMPS_LRMIXED_CAST_24Z(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: B(ldB,ncols)
        INTEGER*1 :: A(2*3*ldA,ncols)
        call mumps_fp24c_zdouble(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_24Z
      SUBROUTINE MUMPS_LRMIXED_CAST_40Z(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: B(ldB,ncols)
        INTEGER*1 :: A(2*5*ldA,ncols)
        call mumps_fp40z_zdouble(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_40Z
      SUBROUTINE MUMPS_LRMIXED_CAST_48Z(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: B(ldB,ncols)
        INTEGER*1 :: A(2*6*ldA,ncols)
        call mumps_fp48z_zdouble(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_48Z
      SUBROUTINE MUMPS_LRMIXED_CAST_56Z(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: B(ldB,ncols)
        INTEGER*1 :: A(2*7*ldA,ncols)
        call mumps_fp56z_zdouble(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_56Z
      ! complex domain {double,float} -> {fp16,...,fp56}
      SUBROUTINE MUMPS_LRMIXED_CAST_C16(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        implicit none
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX :: A(ldA,ncols)
        INTEGER*1 :: B(2*2*ldB,ncols)
        call mumps_cfloat_fp16c(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_C16
      SUBROUTINE MUMPS_LRMIXED_CAST_Z16(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: A(ldA,ncols)
        INTEGER*1 :: B(2*2*ldB,ncols)
        call mumps_zdouble_fp16c(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_Z16
      SUBROUTINE MUMPS_LRMIXED_CAST_C24(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX :: A(ldA,ncols)
        INTEGER*1 :: B(2*3*ldB,ncols)
        call mumps_cfloat_fp24c(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_C24
      SUBROUTINE MUMPS_LRMIXED_CAST_Z24(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: A(ldA,ncols)
        INTEGER*1 :: B(2*3*ldB,ncols)
        call mumps_zdouble_fp24c(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_Z24
      SUBROUTINE MUMPS_LRMIXED_CAST_Z40(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: A(ldA,ncols)
        INTEGER*1 :: B(2*5*ldB,ncols)
        call mumps_zdouble_fp40z(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_Z40
      SUBROUTINE MUMPS_LRMIXED_CAST_Z48(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: A(ldA,ncols)
        INTEGER*1 :: B(2*6*ldB,ncols)
        call mumps_zdouble_fp48z(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_Z48
      SUBROUTINE MUMPS_LRMIXED_CAST_Z56(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: A(ldA,ncols)
        INTEGER*1 :: B(2*7*ldB,ncols)
        call mumps_zdouble_fp56z(A,B,nrows,ncols,ldA,ldB)
      END SUBROUTINE MUMPS_LRMIXED_CAST_Z56
C} MIXHPC
#else
C{ ! MIXHPC
      SUBROUTINE MUMPS_LRMIXED_CAST_S16(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        REAL :: A(ldA,ncols)
        INTEGER*1 :: B(2*ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_32_16(A(I,J), B(2*I-1,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_S16
      SUBROUTINE MUMPS_LRMIXED_CAST_16S(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        REAL :: B(ldB,ncols)
        INTEGER*1 :: A(2*ldA,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_16_32(A(2*I-1,J),B(I,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_16S
      SUBROUTINE MUMPS_LRMIXED_CAST_D16(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: A(ldA,ncols)
        INTEGER*1 :: B(2*ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_32_16(real(A(I,J)), B(2*I-1,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_D16
      SUBROUTINE MUMPS_LRMIXED_CAST_16D(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: B(ldB,ncols)
        REAL :: TMP
        INTEGER*1 :: A(2*ldA,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_16_32(A(2*I-1,J),TMP)
          B(I,J) = DBLE(TMP)
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_16D
      SUBROUTINE MUMPS_LRMIXED_CAST_S24(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        REAL :: A(ldA,ncols)
        INTEGER*1 :: B(3*ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_32_24(A(I,J), B(3*I-2,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_S24
      SUBROUTINE MUMPS_LRMIXED_CAST_24S(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        REAL :: B(ldB,ncols)
        INTEGER*1 :: A(3*ldA,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_24_32(A(3*I-2,J),B(I,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_24S
      SUBROUTINE MUMPS_LRMIXED_CAST_D24(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: A(ldA,ncols)
        INTEGER*1 :: B(3*ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_32_24(real(A(I,J)), B(3*I-2,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_D24
      SUBROUTINE MUMPS_LRMIXED_CAST_24D(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: B(ldB,ncols)
        INTEGER*1 :: A(3*ldA,ncols)
C       Local variables
        REAL :: TMP
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_24_32(A(3*I-2,J),TMP)
          B(I,J)=dble(TMP)
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_24D
      SUBROUTINE MUMPS_LRMIXED_CAST_D40(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: A(ldA,ncols)
        INTEGER*1 :: B(5*ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_64_40(A(I,J), B(5*I-4,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_D40
      SUBROUTINE MUMPS_LRMIXED_CAST_40D(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: B(ldB,ncols)
        INTEGER*1 :: A(5*ldA,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_40_64(A(5*I-4,J),B(I,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_40D
      SUBROUTINE MUMPS_LRMIXED_CAST_D48(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: A(ldA,ncols)
        INTEGER*1 :: B(6*ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_64_48(A(I,J), B(6*I-5,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_D48
      SUBROUTINE MUMPS_LRMIXED_CAST_48D(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: B(ldB,ncols)
        INTEGER*1 :: A(6*ldA,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_48_64(A(6*I-5,J),B(I,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_48D
      SUBROUTINE MUMPS_LRMIXED_CAST_D56(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: A(ldA,ncols)
        INTEGER*1 :: B(7*ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_64_56(A(I,J), B(7*I-6,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_D56
      SUBROUTINE MUMPS_LRMIXED_CAST_56D(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        DOUBLE PRECISION :: B(ldB,ncols)
        INTEGER*1 :: A(7*ldA,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_56_64(A(7*I-6,J),B(I,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_56D
      SUBROUTINE MUMPS_LRMIXED_CAST_C16(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX :: A(ldA,ncols)
        INTEGER*1 :: B(2*2*ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_32_16(real(A(I,J)), B(4*I-3,J))
          CALL CONVERT_32_16(aimag(A(I,J)), B(4*I-1,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_C16
      SUBROUTINE MUMPS_LRMIXED_CAST_16C(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX :: B(ldB,ncols)
        INTEGER*1 :: A(2*2*ldA,ncols)
C       Local variables
        REAL :: TMP1, TMP2
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_16_32(A(4*I-3,J),TMP1)
          CALL CONVERT_16_32(A(4*I-1,J),TMP2)
          B(i,j)=cmplx( TMP1, TMP2 )
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_16C
      SUBROUTINE MUMPS_LRMIXED_CAST_Z16(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: A(ldA,ncols)
        INTEGER*1 :: B(2*2*ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_32_16(real(real(A(I,J))), B(4*I-3,J))
          CALL CONVERT_32_16(real(aimag(A(I,J))), B(4*I-1,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_Z16
      SUBROUTINE MUMPS_LRMIXED_CAST_16Z(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: B(ldB,ncols)
        INTEGER*1 :: A(2*2*ldA,ncols)
C       Local variables
        REAL :: TMP1, TMP2
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_16_32(A(4*I-3,J),TMP1)
          CALL CONVERT_16_32(A(4*I-1,J),TMP2)
          B(I,J)=cmplx(TMP1, TMP2, kind=kind(B))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_16Z
      SUBROUTINE MUMPS_LRMIXED_CAST_C24(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX :: A(ldA,ncols)
        INTEGER*1 :: B(2*3*ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_32_24(real(A(I,J)), B(6*I-5,J))
          CALL CONVERT_32_24(aimag(A(I,J)), B(6*I-2,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_C24
      SUBROUTINE MUMPS_LRMIXED_CAST_24C(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX :: B(ldB,ncols)
        INTEGER*1 :: A(2*3*ldA,ncols)
C       Local variables
        REAL :: TMP1, TMP2
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_24_32(A(6*I-5,J),TMP1)
          CALL CONVERT_24_32(A(6*I-2,J),TMP2)
          B(I,J)=cmplx(TMP1, TMP2)
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_24C
      SUBROUTINE MUMPS_LRMIXED_CAST_Z24(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: A(ldA,ncols)
        INTEGER*1 :: B(2*3*ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_32_24(real(real(A(I,J))), B(6*I-5,J))
          CALL CONVERT_32_24(real(aimag(A(I,J))), B(6*I-2,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_Z24
      SUBROUTINE MUMPS_LRMIXED_CAST_24Z(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: B(ldB,ncols)
        INTEGER*1 :: A(2*3*ldA,ncols)
C       Local variables
        REAL :: TMP1, TMP2
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_24_32(A(6*I-5,J),TMP1)
          CALL CONVERT_24_32(A(6*I-2,J),TMP2)
          B(I,J)=cmplx(TMP1, TMP2, kind=kind(B))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_24Z
      SUBROUTINE MUMPS_LRMIXED_CAST_Z40(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: A(ldA,ncols)
        INTEGER*1 :: B(2*5*ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_64_40(real(A(I,J)), B(10*I-9,J))
          CALL CONVERT_64_40(aimag(A(I,J)), B(10*I-4,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_Z40
      SUBROUTINE MUMPS_LRMIXED_CAST_40Z(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: B(ldB,ncols)
        INTEGER*1 :: A(2*5*ldA,ncols)
C       Local variables
        DOUBLE PRECISION :: TMP1, TMP2
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_40_64(A(10*I-9,J),TMP1)
          CALL CONVERT_40_64(A(10*I-4,J),TMP2)
          B(I,J)=cmplx(TMP1, TMP2, kind=kind(B))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_40Z
      SUBROUTINE MUMPS_LRMIXED_CAST_Z48(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: A(ldA,ncols)
        INTEGER*1 :: B(2*6*ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_64_48(real(A(I,J)), B(12*I-11,J))
          CALL CONVERT_64_48(aimag(A(I,J)), B(12*I-5,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_Z48
      SUBROUTINE MUMPS_LRMIXED_CAST_48Z(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: B(ldB,ncols)
        INTEGER*1 :: A(2*6*ldA,ncols)
C       Local variables
        DOUBLE PRECISION :: TMP1, TMP2
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_48_64(A(12*I-11,J),TMP1)
          CALL CONVERT_48_64(A(12*I-5,J),TMP2)
          B(I,J)=cmplx(TMP1, TMP2, kind(B))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_48Z
      SUBROUTINE MUMPS_LRMIXED_CAST_Z56(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: A(ldA,ncols)
        INTEGER*1 :: B(2*7*ldB,ncols)
C       Local variables
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_64_56(real(A(I,J)), B(14*I-13,J))
          CALL CONVERT_64_56(aimag(A(I,J)), B(14*I-6,J))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_Z56
      SUBROUTINE MUMPS_LRMIXED_CAST_56Z(A, B, nrows, ncols, ldA, 
     &                                      ldB)
        INTEGER :: nrows, ncols, ldA, ldB
        COMPLEX(kind=8) :: B(ldB,ncols)
        INTEGER*1 :: A(2*7*ldA,ncols)
C       Local variables
        DOUBLE PRECISION :: TMP1, TMP2
        INTEGER :: I, J
        DO J=1,ncols
        DO I=1,nrows
          CALL CONVERT_56_64(A(14*I-13,J),TMP1)
          CALL CONVERT_56_64(A(14*I-6,J),TMP2)
          B(I,J)=cmplx(TMP1, TMP2, kind=kind(B))
        ENDDO
        ENDDO
      END SUBROUTINE MUMPS_LRMIXED_CAST_56Z
C} ! MIXHPC
#endif
       SUBROUTINE CONVERT_32_16(A, B)
         INTEGER*4 A
         INTEGER*4 TMP
         INTEGER*1 B(2) 
         TMP = 0
         CALL MVBITS(A, 16, 8, TMP, 0)
         B(1)=TMP
         CALL MVBITS(A, 24, 8, TMP, 0)
         B(2)=TMP
       END SUBROUTINE CONVERT_32_16
       SUBROUTINE CONVERT_16_32(A, B)
         INTEGER*4 B
         INTEGER*4 TMP
         INTEGER*1 A(3) 
         B = 32768
         TMP=A(1)
         CALL MVBITS(TMP, 0, 8, B, 16)
         TMP=A(2)
         CALL MVBITS(TMP, 0, 8, B, 24)
       END SUBROUTINE CONVERT_16_32
       SUBROUTINE CONVERT_64_40(A, B)
         INTEGER*8 A
         INTEGER*8 TMP
         INTEGER*1 B(5) 
         TMP = 0
         CALL MVBITS(A, 24, 8, TMP, 0)
         B(1)=TMP
         CALL MVBITS(A, 32, 8, TMP, 0)
         B(2)=TMP
         CALL MVBITS(A, 40, 8, TMP, 0)
         B(3)=TMP
         CALL MVBITS(A, 48, 8, TMP, 0)
         B(4)=TMP
         CALL MVBITS(A, 56, 8, TMP, 0)
         B(5)=TMP
       END SUBROUTINE CONVERT_64_40
       SUBROUTINE CONVERT_40_64(A, B)
         INTEGER*8 B
         INTEGER*8 TMP
         INTEGER*1 A(5) 
         B = 8388608_8
         TMP=A(1)
         CALL MVBITS(TMP, 0, 8, B, 24)
         TMP=A(2)
         CALL MVBITS(TMP, 0, 8, B, 32)
         TMP=A(3)
         CALL MVBITS(TMP, 0, 8, B, 40)
         TMP=A(4)
         CALL MVBITS(TMP, 0, 8, B, 48)
         TMP=A(5)
         CALL MVBITS(TMP, 0, 8, B, 56)
       END SUBROUTINE CONVERT_40_64
       SUBROUTINE CONVERT_32_24(A, B)
         INTEGER*4 A
         INTEGER*4 TMP
         INTEGER*1 B(3) 
         TMP = 0
         CALL MVBITS(A, 8, 8, TMP, 0)
         B(1)=TMP
         CALL MVBITS(A, 16, 8, TMP, 0)
         B(2)=TMP
         CALL MVBITS(A, 24, 8, TMP, 0)
         B(3)=TMP
       END SUBROUTINE CONVERT_32_24
       SUBROUTINE CONVERT_24_32(A, B)
         INTEGER*4 B
         INTEGER*4 TMP
         INTEGER*1 A(3) 
         B = 128
         TMP=A(1)
         CALL MVBITS(TMP, 0, 8, B, 8)
         TMP=A(2)
         CALL MVBITS(TMP, 0, 8, B, 16)
         TMP=A(3)
         CALL MVBITS(TMP, 0, 8, B, 24)
       END SUBROUTINE CONVERT_24_32
       SUBROUTINE CONVERT_64_48(A, B)
         INTEGER*8 A
         INTEGER*8 TMP
         INTEGER*1 B(6) 
         TMP = 0
         CALL MVBITS(A, 16, 8, TMP, 0)
         B(1)=TMP
         CALL MVBITS(A, 24, 8, TMP, 0)
         B(2)=TMP
         CALL MVBITS(A, 32, 8, TMP, 0)
         B(3)=TMP
         CALL MVBITS(A, 40, 8, TMP, 0)
         B(4)=TMP
         CALL MVBITS(A, 48, 8, TMP, 0)
         B(5)=TMP
         CALL MVBITS(A, 56, 8, TMP, 0)
         B(6)=TMP
       END SUBROUTINE CONVERT_64_48
       SUBROUTINE CONVERT_48_64(A, B)
         INTEGER*8 B
         INTEGER*8 TMP
         INTEGER*1 A(6) 
         B = 32768_8
         TMP=A(1)
         CALL MVBITS(TMP, 0, 8, B, 16)
         TMP=A(2)
         CALL MVBITS(TMP, 0, 8, B, 24)
         TMP=A(3)
         CALL MVBITS(TMP, 0, 8, B, 32)
         TMP=A(4)
         CALL MVBITS(TMP, 0, 8, B, 40)
         TMP=A(5)
         CALL MVBITS(TMP, 0, 8, B, 48)
         TMP=A(6)
         CALL MVBITS(TMP, 0, 8, B, 56)
       END SUBROUTINE CONVERT_48_64
       SUBROUTINE CONVERT_64_56(A, B)
         INTEGER*8 A
         INTEGER*8 TMP
         INTEGER*1 B(7) 
         TMP = 0
         CALL MVBITS(A, 8, 8, TMP, 0)
         B(1)=TMP
         CALL MVBITS(A, 16, 8, TMP, 0)
         B(2)=TMP
         CALL MVBITS(A, 24, 8, TMP, 0)
         B(3)=TMP
         CALL MVBITS(A, 32, 8, TMP, 0)
         B(4)=TMP
         CALL MVBITS(A, 40, 8, TMP, 0)
         B(5)=TMP
         CALL MVBITS(A, 48, 8, TMP, 0)
         B(6)=TMP
         CALL MVBITS(A, 56, 8, TMP, 0)
         B(7)=TMP
       END SUBROUTINE CONVERT_64_56
       SUBROUTINE CONVERT_56_64(A, B)
         INTEGER*8 B
         INTEGER*8 TMP
         INTEGER*1 A(7) 
         B = 128
         TMP=A(1)
         CALL MVBITS(TMP, 0, 8, B, 8)
         TMP=A(2)
         CALL MVBITS(TMP, 0, 8, B, 16)
         TMP=A(3)
         CALL MVBITS(TMP, 0, 8, B, 24)
         TMP=A(4)
         CALL MVBITS(TMP, 0, 8, B, 32)
         TMP=A(5)
         CALL MVBITS(TMP, 0, 8, B, 40)
         TMP=A(6)
         CALL MVBITS(TMP, 0, 8, B, 48)
         TMP=A(7)
         CALL MVBITS(TMP, 0, 8, B, 56)
       END SUBROUTINE CONVERT_56_64
       SUBROUTINE MUMPS_LRMIXED_PLUS_ZZ(X, Y, nrows, ncols, ldX, 
     &                                       ldY)
         INTEGER :: nrows, ncols, ldX, ldY
         COMPLEX(kind=8) :: X(ldX,ncols), Y(ldY,ncols)
C        Local variables
         INTEGER :: I, J
         DO J=1,ncols
         DO I=1,nrows
           Y(I,J) =  Y(I,J) + X(I,J)
         ENDDO
         ENDDO
       END SUBROUTINE MUMPS_LRMIXED_PLUS_ZZ
       SUBROUTINE MUMPS_LRMIXED_PLUS_DD(X, Y, nrows, ncols, ldX, 
     &                                       ldY)
         INTEGER :: nrows, ncols, ldX, ldY
         DOUBLE PRECISION :: X(ldX,ncols), Y(ldY,ncols)
C        Local variables
         INTEGER :: I, J
         DO J=1,ncols
         DO I=1,nrows
           Y(I,J) =  Y(I,J) + X(I,J)
         ENDDO
         ENDDO
       END SUBROUTINE MUMPS_LRMIXED_PLUS_DD
       SUBROUTINE MUMPS_LRMIXED_PLUS_ZC(X, Y, nrows, ncols, ldX, 
     &                                       ldY)
         INTEGER :: nrows, ncols, ldX, ldY
         COMPLEX(kind=8) :: X(ldX,ncols)
         COMPLEX :: Y(ldY,ncols)
C        Local variables
         INTEGER :: I, J
         DO J=1,ncols
         DO I=1,nrows
           Y(I,J) =  Y(I,J) + cmplx(X(I,J))
         ENDDO
         ENDDO
       END SUBROUTINE MUMPS_LRMIXED_PLUS_ZC
       SUBROUTINE MUMPS_LRMIXED_PLUS_DS(X, Y, nrows, ncols, ldX, 
     &                                       ldY)
         INTEGER :: nrows, ncols, ldX, ldY
         DOUBLE PRECISION :: X(ldX,ncols)
         REAL :: Y(ldY,ncols)
C        Local variables
         INTEGER :: I, J
         DO J=1,ncols
         DO I=1,nrows
           Y(I,J) =  Y(I,J) + real(X(I,J))
         ENDDO
         ENDDO
       END SUBROUTINE MUMPS_LRMIXED_PLUS_DS
       SUBROUTINE MUMPS_LRMIXED_PLUS_CZ(X, Y, nrows, ncols, ldX, 
     &                                       ldY)
         INTEGER :: nrows, ncols, ldX, ldY
         COMPLEX :: X(ldX,ncols)
         COMPLEX(kind=8) :: Y(ldY,ncols)
C        Local variables
         INTEGER :: I, J
         DO J=1,ncols
         DO I=1,nrows
           Y(I,J) =  Y(I,J) + cmplx(X(I,J), kind=kind(Y))
         ENDDO
         ENDDO
       END SUBROUTINE MUMPS_LRMIXED_PLUS_CZ
       SUBROUTINE MUMPS_LRMIXED_PLUS_SD(X, Y, nrows, ncols, ldX, 
     &                                       ldY)
         INTEGER :: nrows, ncols, ldX, ldY
         REAL :: X(ldX,ncols)
         DOUBLE PRECISION :: Y(ldY,ncols)
C        Local variables
         INTEGER :: I, J
         DO J=1,ncols
         DO I=1,nrows
           Y(I,J) =  Y(I,J) + dble(X(I,J))
         ENDDO
         ENDDO
       END SUBROUTINE MUMPS_LRMIXED_PLUS_SD
       SUBROUTINE MUMPS_LRMIXED_PLUS_CC(X, Y, nrows, ncols, ldX, 
     &                                       ldY)
         INTEGER :: nrows, ncols, ldX, ldY
         COMPLEX :: X(ldX,ncols), Y(ldY,ncols)
C        Local variables
         INTEGER :: I, J
         DO J=1,ncols
         DO I=1,nrows
           Y(I,J) =  Y(I,J) + X(I,J)
         ENDDO
         ENDDO
       END SUBROUTINE MUMPS_LRMIXED_PLUS_CC
       SUBROUTINE MUMPS_LRMIXED_PLUS_SS(X, Y, nrows, ncols, ldX, 
     &                                       ldY)
         INTEGER :: nrows, ncols, ldX, ldY
         REAL :: X(ldX,ncols), Y(ldY,ncols)
C        Local variables
         INTEGER :: I, J
         DO J=1,ncols
         DO I=1,nrows
           Y(I,J) =  Y(I,J) + X(I,J)
         ENDDO
         ENDDO
       END SUBROUTINE MUMPS_LRMIXED_PLUS_SS
       SUBROUTINE MUMPS_LRMIXED_PACK(X, NROWS, NCOLS, IROW_BEG,
     &                               NROWS_TO_PACK, NBYTES_PER_ELT,
     &                         BUF, LBUF, POSITION, COMM, IERR)
         INTEGER :: NBYTES_PER_ELT, NROWS, NCOLS, IROW_BEG,
     &              NROWS_TO_PACK
         INTEGER*1 :: X(NROWS*NCOLS*NBYTES_PER_ELT)
         INTEGER, intent(in)  :: COMM, LBUF 
         INTEGER, intent(inout) :: POSITION, IERR
         INTEGER, intent(inout) :: BUF(*)
c        Local variables         
         INTEGER :: J
         INCLUDE 'mpif.h'
         DO J=1,NCOLS
          CALL MPI_PACK( 
     &      X(1+(IROW_BEG-1+(J-1)*NROWS)*NBYTES_PER_ELT),
     &      NROWS_TO_PACK*NBYTES_PER_ELT, MPI_BYTE,
     &      BUF(1), LBUF, POSITION, COMM, IERR)
         ENDDO
       END SUBROUTINE MUMPS_LRMIXED_PACK
       SUBROUTINE MUMPS_LRMIXED_UNPACK(X, NROWS, NCOLS,
     &                                 NBYTES_PER_ELT,
     &                         BUF, LBUF, POSITION, COMM, IERR)
         INTEGER :: NBYTES_PER_ELT, NROWS, NCOLS
         INTEGER*1 :: X(NROWS*NCOLS*NBYTES_PER_ELT)
         INTEGER, intent(in)  :: COMM, LBUF 
         INTEGER, intent(inout) :: POSITION, IERR
         INTEGER, intent(in) :: BUF(*)
         INTEGER :: J
         INCLUDE 'mpif.h'
         DO J=1,NCOLS
           CALL MPI_UNPACK(BUF, LBUF, POSITION,
     &      X(1+(J-1)*NROWS*NBYTES_PER_ELT), NROWS*NBYTES_PER_ELT,
     &     MPI_BYTE, COMM, IERR)
         ENDDO
       END SUBROUTINE MUMPS_LRMIXED_UNPACK
       SUBROUTINE MUMPS_LRMIXED_DGEMM(TRANSA, TRANSB, M, N, K, ALPHA,
     &               A, LDA, B, LDB, BETA, C, LDC)
         INTEGER :: M, N, K, LDA, LDB, LDC
         DOUBLE PRECISION :: A(LDA,*), B(LDB,*), C(LDC,*)
         DOUBLE PRECISION :: ALPHA, BETA
         CHARACTER(len=1) :: TRANSA, TRANSB
         CALL dgemm(TRANSA, TRANSB, M, N, K, ALPHA, 
     &             A, LDA, B, LDB, BETA, C, LDC)
       END SUBROUTINE MUMPS_LRMIXED_DGEMM
CAJ2ALL why are alpha and beta in double precision everywhere ?
CAJ2ALL blas defines alpha/beta in the same precision as the arrays
       SUBROUTINE MUMPS_LRMIXED_SGEMM(TRANSA, TRANSB, M, N, K, ALPHA,
     &               A, LDA, B, LDB, BETA, C, LDC)
         INTEGER :: M, N, K, LDA, LDB, LDC
         REAL :: A(LDA,*), B(LDB,*), C(LDC,*)
         DOUBLE PRECISION :: ALPHA, BETA
         CHARACTER(len=1) :: TRANSA, TRANSB
         CALL sgemm(TRANSA, TRANSB, M, N, K, real(ALPHA),
     &              A, LDA, B, LDB, real(BETA), C, LDC)
       END SUBROUTINE MUMPS_LRMIXED_SGEMM
       SUBROUTINE MUMPS_LRMIXED_ZGEMM(TRANSA, TRANSB, M, N, K, ALPHA,
     &               A, LDA, B, LDB, BETA, C, LDC)
         INTEGER :: M, N, K, LDA, LDB, LDC
         COMPLEX(kind=8) :: A(LDA,*), B(LDB,*), C(LDC,*)
         DOUBLE PRECISION :: ALPHA, BETA
         CHARACTER(len=1) :: TRANSA, TRANSB
         CALL zgemm(TRANSA,TRANSB, M, N, K, cmplx(ALPHA,kind=kind(A)),
     &              A, LDA, B, LDB, cmplx(BETA,kind=kind(A)), C, LDC)
       END SUBROUTINE MUMPS_LRMIXED_ZGEMM
       SUBROUTINE MUMPS_LRMIXED_CGEMM(TRANSA, TRANSB, M, N, K, ALPHA,
     &               A, LDA, B, LDB, BETA, C, LDC)
         INTEGER :: M, N, K, LDA, LDB, LDC
         COMPLEX :: A(LDA,*), B(LDB,*), C(LDC,*)
         DOUBLE PRECISION :: ALPHA, BETA
         CHARACTER(len=1) :: TRANSA, TRANSB
         CALL cgemm(TRANSA, TRANSB, M, N, K, cmplx(ALPHA),
     &              A, LDA, B, LDB, cmplx(BETA), C, LDC)
       END SUBROUTINE MUMPS_LRMIXED_CGEMM
