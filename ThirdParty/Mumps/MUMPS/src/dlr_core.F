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
C
C Note: the last routine of this file, xMUMPS_TRUNCATED_RRQR is derived from
C the LAPACK package, for which BSD 3-clause license applies 
C (see header of the routine).
      MODULE DMUMPS_LR_CORE
      USE MUMPS_LR_COMMON
      USE DMUMPS_LR_TYPE
      USE MUMPS_LR_STATS
      USE DMUMPS_LR_DATA_M
!$    USE OMP_LIB
      IMPLICIT NONE
      CONTAINS
C *******************************!
C ALLOCATION!
C *******************************!
      SUBROUTINE INIT_LRB(LRB_OUT,K,M,N,ISLR)
C This routine simply initializes a LR block but does NOT allocate it
C (allocation occurs somewhere else)
        TYPE(LRB_TYPE), INTENT(OUT) :: LRB_OUT
        INTEGER,INTENT(IN) :: K,M,N
        LOGICAL,INTENT(IN) :: ISLR
        LRB_OUT%M = M
        LRB_OUT%N = N
        LRB_OUT%K = K
        LRB_OUT%ISLR = ISLR
        LRB_OUT%KW_STOC = 0
        NULLIFY(LRB_OUT%Q)
        NULLIFY(LRB_OUT%R)
      END SUBROUTINE INIT_LRB
      SUBROUTINE FRB_GET_PREC(LRB, PREC, LRMIXED_ACTIVATED)
        IMPLICIT NONE
        TYPE(LRB_TYPE), intent(in) :: LRB
        INTEGER, intent(out) :: PREC
        LOGICAL, intent(in) :: LRMIXED_ACTIVATED
        PREC = PREC_FAC
#if !defined(NOMIXHPC)
#else
        RETURN
#endif
        IF (LRB%M*LRB%N.LE.0) RETURN
        IF (.not.LRMIXED_ACTIVATED) THEN
          PREC = PREC_FAC
        ELSE IF (.not.QUANTIZE) THEN
          PREC = PREC_FR
        ELSE IF (ASSOCIATED(LRB%K_STOC)) THEN
        !ELSE IF (ASSOCIATED(LRB%Q_STOC)) THEN
          DO PREC=1,NPREC_STOC
            IF (LRB%K_STOC(PREC).GT.0) THEN
            !IF (associated(LRB%Q_STOC(PREC)%VAL)) THEN
              EXIT
            END IF
          END DO  
          IF (PREC.EQ.NPREC_STOC+1) PREC = PREC_FAC
        ENDIF
      END SUBROUTINE FRB_GET_PREC
C
C
      SUBROUTINE IS_FRONT_BLR_CANDIDATE(INODE, NIV, NFRONT, NASS, 
     &                    BLRON, K473, 
     &                    K489,
     &                    K490, K491, K492, K20, K60, IDAD, K38,
     &                    K123, LRSTATUS, K280, LRGROUPS)
        INTEGER,INTENT(IN) :: INODE, NFRONT, NASS, BLRON, K123,
     &                        K473, 
     &                        K489, K490,
     &                        K491, K492, NIV, K20, K60, IDAD, K38
        INTEGER,INTENT(OUT):: LRSTATUS 
        INTEGER, INTENT(IN):: K280
        INTEGER, INTENT(IN), OPTIONAL :: LRGROUPS(K280)
C       LRSTATUS : see description in mumps_header.h
C         0=FR, 
C         1=LowRank CB only
C         2=LowRank factors/panels only
C         3=LowRank CB+factor/panel
C         4=FR factor/panel
C         5=FR CB+factor/panel
C
C     Local variables
        LOGICAL :: COMPRESS_PANEL, COMPRESS_CB, COULD_COMPRESS_CB
        LRSTATUS = 0
C       Type 3 node is not BLR
        IF (NIV.EQ.3) THEN
         goto 100
        ENDIF
        COMPRESS_PANEL = .FALSE.
        COMPRESS_CB = .FALSE.
        IF (BLRON.NE.0) THEN
            IF ( (K492.LT.0.and.INODE.EQ.abs(K492))
     &           .or.
     &           (K492.GT.0.and.K491.LE.NFRONT.and.K490.LE.NASS)) THEN
              COMPRESS_PANEL = .TRUE.
C             Compression for NASS =1 is useless
              IF (NASS.LE.1) THEN
                COMPRESS_PANEL =.FALSE. 
              ENDIF
              IF (present(LRGROUPS)) THEN
C              If the lr-group of node INODE is negative, then the current
C              front was supposed to have a small separator; even if it is
C              bigger than expected (e.g., because of some delayed pivots), 
C              it probably is more robust NOT to treat this front as LR
C              (but some tests could be done on this matter)            
C              Also note that if the front was selected with KEEP(492) but
C              its separator is too small, it is NOT treated as LR          
               IF (LRGROUPS (INODE) .LT. 0) COMPRESS_PANEL = .FALSE.
              ENDIF
            ENDIF
            IF (  (K489.GT.0.AND.(K489.NE.2.OR.NIV.EQ.2))
     &            .and. 
     &            (  (K492.LT.0.and.INODE.EQ.abs(K492))
     &               .or.
     &               (K492.GT.0.AND.NFRONT-NASS.GT.K491)
     &            ) ) 
     &         THEN
C             We want to select nodes that have a large CB
              COMPRESS_CB = .TRUE.
            ENDIF
        ENDIF
        COULD_COMPRESS_CB = COMPRESS_CB
        IF (.NOT.COMPRESS_PANEL) COMPRESS_CB=.FALSE.
        IF (COMPRESS_PANEL.OR.COMPRESS_CB) THEN
          IF (COMPRESS_CB.AND.(.NOT.COMPRESS_PANEL)) THEN
C           Full rank panel factorization and CB compressed
            LRSTATUS = 1 
          ELSE IF (COMPRESS_PANEL.AND.(.NOT.COMPRESS_CB)) THEN
C           Low rank panel factorization and CB not compressed
            LRSTATUS = 2 
          ELSE
C           Low rank panel factorization and CB compressed
            LRSTATUS = 3 
          ENDIF
        ELSE
C        Full rank factorization and CB not compressed
         LRSTATUS = 0
        ENDIF
C
C       Schur complement cannot be BLR for now
C
        IF ( INODE .EQ. K20 .AND. K60 .NE. 0 ) THEN
          LRSTATUS = 0
        ENDIF
C
C       Do not compress CB of children of root
C
        IF ( IDAD .EQ. K38 .AND. K38 .NE.0 ) THEN
            COMPRESS_CB = .FALSE.
            IF (LRSTATUS.EQ.3) THEN
               LRSTATUS = 2
            ELSE IF (LRSTATUS.EQ.1) THEN
               LRSTATUS = 0
            ! no other cases
            ENDIF
        ENDIF
100     continue
        ! K473=2 : all full-rank becomes BFR, BLR stays BLR
        IF (K473.EQ.2) THEN
            IF (LRSTATUS.EQ.0) LRSTATUS = 4
            IF (LRSTATUS.EQ.1) LRSTATUS = 5
        ! K473=3 : all fronts are BFR
        ELSE IF (K473.EQ.3) THEN
            IF (LRSTATUS.EQ.0.OR.LRSTATUS.EQ.2) LRSTATUS = 4
            IF (LRSTATUS.EQ.1.OR.LRSTATUS.EQ.3) LRSTATUS = 5
        ENDIF
C    
C       Adjust LRSTATUS:
        IF (LRSTATUS.EQ.5) THEN 
          LRSTATUS = 4
        ENDIF
      RETURN
      END SUBROUTINE IS_FRONT_BLR_CANDIDATE
      SUBROUTINE ALLOC_LRB(LRB_OUT,K,M,N,ISLR,IFLAG,IERROR,KEEP8)
        TYPE(LRB_TYPE), INTENT(OUT) :: LRB_OUT
        INTEGER,INTENT(IN) :: K,M,N
        INTEGER,INTENT(INOUT) :: IFLAG, IERROR
        LOGICAL,INTENT(IN) :: ISLR
        INTEGER(8) :: KEEP8(150)
        INTEGER :: MEM, allocok
        DOUBLE PRECISION :: ZERO
        PARAMETER (ZERO = 0.0D0)
C this routine takes care of both
C   - the initialization of the LR block structure/pointers        
C   - the allocation of the LR block        
        LRB_OUT%M = M
        LRB_OUT%N = N
        LRB_OUT%K = K
        LRB_OUT%ISLR = ISLR
#if !defined(NOMIXHPC)
        LRB_OUT%kw_stoc = 0
#endif
        IF ((M.EQ.0).OR.(N.EQ.0)) THEN 
         nullify(LRB_OUT%Q)
         nullify(LRB_OUT%R)
         RETURN
        ENDIF
        IF (ISLR) THEN
          IF (K.EQ.0) THEN
            nullify(LRB_OUT%Q)
            nullify(LRB_OUT%R)
          ELSE
            allocate(LRB_OUT%Q(M,K),LRB_OUT%R(K,N),stat=allocok)
            IF (allocok > 0) THEN
              IFLAG  = -13
              IERROR = K*(M+N)
              RETURN
            ENDIF
          ENDIF
        ELSE
          nullify(LRB_OUT%R)
          allocate(LRB_OUT%Q(M,N),stat=allocok)
          IF (allocok > 0) THEN
            IFLAG  = -13
            IERROR = M*N
            RETURN
          ENDIF
        ENDIF
C       Update memory consumption
        IF (ISLR) THEN
          MEM = M*K + N*K
        ELSE
          MEM = M*N
        ENDIF
C       ATOMIC needed even above L0-threads because several
C       threads may call ALLOC_LRB simultaneously on a front
        CALL MUMPS_DM_FAC_UPD_DYN_MEMCNTS(int(MEM,8),
     &        .TRUE., KEEP8, IFLAG, IERROR, .TRUE., .TRUE.)
      RETURN
      END SUBROUTINE ALLOC_LRB
      SUBROUTINE ALLOC_LRB_FROM_ACC(ACC_LRB, LRB_OUT, K, M, N, LorU,
     &                          IFLAG, IERROR, KEEP8)
C       This routine allocates a LRB block, initializing it with the
C       TRANSPOSE of the accumulator ACC_LRB. 
C       (We transpose it because the accumulator is a LR representation
C       of A^T, not A)        
C       Is called in the context of KEEP(480)>=5 LUA strategy.
        TYPE(LRB_TYPE), INTENT(IN) :: ACC_LRB
        TYPE(LRB_TYPE), INTENT(OUT) :: LRB_OUT
        INTEGER,INTENT(IN) :: K, M, N, LorU
        INTEGER,INTENT(INOUT) :: IFLAG, IERROR
        INTEGER(8) :: KEEP8(150)
        INTEGER :: I
C       We want to compress U^T, therefore no need to tranpose the
C       accumulator, as it is already a LR representation of U^T        
        IF (LorU.EQ.1) THEN
          CALL ALLOC_LRB(LRB_OUT,K,M,N,.TRUE.,IFLAG,IERROR,KEEP8)
          IF (IFLAG.LT.0) RETURN
          DO I=1,K
            LRB_OUT%Q(1:M,I) = ACC_LRB%Q(1:M,I)
            LRB_OUT%R(I,1:N) = -ACC_LRB%R(I,1:N)
          ENDDO
        ELSE
          CALL ALLOC_LRB(LRB_OUT,K,N,M,.TRUE.,IFLAG,IERROR,KEEP8)
          IF (IFLAG.LT.0) RETURN
          DO I=1,K
            LRB_OUT%Q(1:N,I) = ACC_LRB%R(I,1:N)
            LRB_OUT%R(I,1:M) = -ACC_LRB%Q(1:M,I)
          ENDDO
        ENDIF
      END SUBROUTINE ALLOC_LRB_FROM_ACC
      SUBROUTINE REGROUPING2(CUT, NPARTSASS, NASS,
     &                   NPARTSCB, NCB, IBCKSZ, ONLYCB, K472, 
     &                   NFRONT, KEEP)
C
C Because the grouping is done during analysis, we observe small groups
C in the CB of the fronts, which not very good for low-rank and not good
C at all for performance. We thus want to regroup those small groups in
C order to always have an acceptable group size.
C
C this is a basic strategy where we merge neighbors group together until
C the resulting group is large enough.
C
C  Note that NCB might be 0 to indicate that only NASS entries are 
C  considered for regrouping and thus 
C  NFRONT might be different from NCB + NASS !
C  NFRONT only needed for COMPUTE_BLR_VCS
        INTEGER, INTENT(IN) :: IBCKSZ, NASS, NCB, NFRONT, KEEP(500)
        INTEGER, INTENT(INOUT) :: NPARTSCB, NPARTSASS
        INTEGER, POINTER, DIMENSION(:) :: CUT
        INTEGER, POINTER, DIMENSION(:) :: NEW_CUT
        INTEGER :: I, INEW, MINSIZE, NEW_NPARTSASS, allocok
        LOGICAL :: ONLYCB, TRACE
        INTEGER, INTENT(IN) :: K472
        INTEGER :: IBCKSZ2,IFLAG,IERROR
        ALLOCATE(NEW_CUT(max(NPARTSASS,1)+NPARTSCB+1),stat=allocok)
        IF (allocok > 0) THEN
           IFLAG  = -13
           IERROR = max(NPARTSASS,1)+NPARTSCB+1
           write(*,*) 'Allocation problem in BLR routine REGROUPING2:',
     &          ' not enough memory? memory requested = ' , IERROR
           RETURN
        ENDIF
        CALL COMPUTE_BLR_VCS(K472, IBCKSZ2, IBCKSZ, NASS, 
     &                       NFRONT, KEEP(35))
        MINSIZE = int(IBCKSZ2 / 3)
        NEW_NPARTSASS = max(NPARTSASS,1)
        IF (.NOT. ONLYCB) THEN
           NEW_CUT(1) = 1
           INEW = 2
           I = 2
C          Large enough (>MINSIZE) groups are 
C          preserved. 
C          Small groups are regrouped 
C          until reaching a large enough group. 
           DO WHILE (I .LE. NPARTSASS + 1)
              NEW_CUT(INEW) = CUT(I)
              TRACE = .FALSE.
              IF (NEW_CUT(INEW) - NEW_CUT(INEW-1) .GT. MINSIZE) THEN
                 INEW = INEW + 1
                 TRACE = .TRUE.
              ENDIF
              I = I + 1
           END DO
           IF (TRACE) THEN
              INEW = INEW - 1 
           ELSE
C             If last group is small and is not unique
C             then regroup it with previous group
              IF (INEW .NE. 2) THEN
                 NEW_CUT(INEW-1) = NEW_CUT(INEW)
                 INEW = INEW - 1
              ENDIF
           ENDIF
           NEW_NPARTSASS = INEW - 1
        ENDIF
        IF (ONLYCB) THEN
           DO I=1,max(NPARTSASS,1)+1
              NEW_CUT(I) = CUT(I)
           ENDDO
        ENDIF
        IF (NCB .EQ. 0) GO TO 50
        INEW = NEW_NPARTSASS+2
        I = max(NPARTSASS,1) + 2
        DO WHILE (I .LE. max(NPARTSASS,1) + NPARTSCB + 1)
              NEW_CUT(INEW) = CUT(I)
              TRACE = .FALSE.
              IF (NEW_CUT(INEW) - NEW_CUT(INEW-1) .GT. MINSIZE) THEN
                 INEW = INEW + 1
                 TRACE = .TRUE.
              ENDIF
              I = I + 1
        END DO
        IF (TRACE) THEN
           INEW = INEW - 1 
        ELSE
           IF (INEW .NE.  NEW_NPARTSASS+2) THEN
           NEW_CUT(INEW-1) = NEW_CUT(INEW)
              INEW = INEW - 1
           ENDIF
        ENDIF
        NPARTSCB = INEW - 1 - NEW_NPARTSASS
 50     CONTINUE       
        NPARTSASS = NEW_NPARTSASS
        DEALLOCATE(CUT)
        ALLOCATE(CUT(NPARTSASS+NPARTSCB+1),stat=allocok)
        IF (allocok > 0) THEN
           IFLAG  = -13
           IERROR = NPARTSASS+NPARTSCB+1
           write(*,*) 'Allocation problem in BLR routine REGROUPING2:',
     &          ' not enough memory? memory requested = ' , IERROR
           RETURN
        ENDIF
        DO I=1,NPARTSASS+NPARTSCB+1
           CUT(I) = NEW_CUT(I)
        ENDDO
        DEALLOCATE(NEW_CUT)
      END SUBROUTINE REGROUPING2
      SUBROUTINE UPD_MRY_LU_LRGAIN( BLR_PANEL, NBBLOCKS
     &           , LRMIXED_ACTIVATED
     &           )
C Updates the memory gain associated with a given BLR panel
        INTEGER,INTENT(IN) :: NBBLOCKS
        TYPE(LRB_TYPE),  INTENT(IN) :: BLR_PANEL(:)
        LOGICAL, INTENT(IN) :: LRMIXED_ACTIVATED
        DOUBLE PRECISION :: MRY_mixed
        INTEGER :: PREC, PREC_FRB
        DOUBLE PRECISION :: PREC_ENTRIES
        DOUBLE PRECISION :: MRY
        INTEGER :: I
C
        MRY    = 0.0D0
        MRY_mixed = 0.0D0
        DO I = 1, NBBLOCKS
          IF (BLR_PANEL(I)%ISLR) THEN
            MRY = MRY + dble(BLR_PANEL(I)%M*BLR_PANEL(I)%N 
     &            - BLR_PANEL(I)%K*(BLR_PANEL(I)%M + BLR_PANEL(I)%N)) 
            IF (LRMIXED_ACTIVATED) THEN
              MRY_mixed=MRY_mixed + dble(BLR_PANEL(I)%M*BLR_PANEL(I)%N)
              DO PREC=1,NPREC_STOC
                MRY_mixed=MRY_mixed - dble(BLR_PANEL(I)%K_STOC(PREC))
     &          *dble(COSTS_PREC(PREC))/dble(COSTS_PREC(1))
     &          *dble(BLR_PANEL(I)%M + BLR_PANEL(I)%N)
              END DO
            ENDIF
          ELSE ! islr
            MRY = MRY + 0.0d0
            IF (LRMIXED_ACTIVATED) THEN
                CALL FRB_GET_PREC(BLR_PANEL(I), PREC_FRB, 
     &              LRMIXED_ACTIVATED)
                PREC_ENTRIES = dble(BLR_PANEL(I)%M*BLR_PANEL(I)%N)
                MRY_mixed=MRY_mixed + PREC_ENTRIES
                MRY_mixed=MRY_mixed
     &             - dble(COSTS_PREC(PREC_FRB))/dble(COSTS_PREC(1))
     &             * PREC_ENTRIES
            ENDIF
          ENDIF ! islr
        ENDDO
!$OMP   ATOMIC UPDATE
        MRY_LU_LRGAIN  = MRY_LU_LRGAIN  + MRY
!$OMP   END ATOMIC
        IF (LRMIXED_ACTIVATED) THEN
!$OMP     ATOMIC UPDATE
          MRY_LU_LRGAIN_mixed  = MRY_LU_LRGAIN_mixed  + MRY_mixed
!$OMP     END ATOMIC
        ENDIF
      RETURN
      END SUBROUTINE UPD_MRY_LU_LRGAIN
      SUBROUTINE DMUMPS_LRTRSM(A, LA, POSELT_LOCAL, NFRONT, LDA, LRB, 
     &                NIV, SYM, LorU, IW, OFFSET_IW) 
C     -----------
C     Parameters
C     -----------
      INTEGER(8), intent(in)  :: LA
      INTEGER, intent(in)     :: NFRONT, NIV, SYM, LorU, LDA
      INTEGER(8), intent(in)  :: POSELT_LOCAL
      DOUBLE PRECISION, intent(inout)  :: A(LA)
      TYPE(LRB_TYPE), intent(inout)   :: LRB
      INTEGER, OPTIONAL:: OFFSET_IW
      INTEGER, OPTIONAL :: IW(*)
C     -----------
C     Local variables
C     -----------
      INTEGER(8) :: DPOS, POSPV1, POSPV2, OFFDAG
      INTEGER    :: M, N, I, J
      DOUBLE PRECISION, POINTER  :: LR_BLOCK_PTR(:,:)
      DOUBLE PRECISION :: ONE, MONE, ZERO
      DOUBLE PRECISION :: MULT1, MULT2, A11, DETPIV, A22, A12
      PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      N  = LRB%N
      IF (LRB%ISLR) THEN
        M  = LRB%K
        LR_BLOCK_PTR => LRB%R
      ELSE
        M  = LRB%M
        LR_BLOCK_PTR => LRB%Q
      END IF
      IF (M.NE.0) THEN
C           Why is it Right, Lower, Tranpose? 
C           Because A is stored by rows 
C           but BLR_L is stored by columns          
        IF (SYM.EQ.0.AND.LorU.EQ.0) THEN
          CALL dtrsm('R', 'L', 'T', 'N', M, N, ONE,
     &                 A(POSELT_LOCAL), NFRONT,
     &                 LR_BLOCK_PTR(1,1), M)
        ELSE
          CALL dtrsm('R', 'U', 'N', 'U', M, N, ONE,
     &                 A(POSELT_LOCAL), LDA,
     &                 LR_BLOCK_PTR(1,1), M)
          IF (LorU.EQ.0) THEN
C           Now apply D scaling              
            IF (.NOT.present(OFFSET_IW)) THEN
              write(*,*) 'Internal error in ',
     &             'DMUMPS_LRTRSM' 
              CALL MUMPS_ABORT()
            ENDIF
            DPOS = POSELT_LOCAL
            I = 1
            DO
            IF(I .GT. N) EXIT
            IF(IW(OFFSET_IW+I-1) .GT. 0) THEN
C             1x1 pivot
              A11 = ONE/A(DPOS)
              CALL dscal(M, A11, LR_BLOCK_PTR(1,I), 1)
              DPOS = DPOS + int(LDA + 1,8)
              I = I+1
            ELSE
C             2x2 pivot
              POSPV1 = DPOS
              POSPV2 = DPOS+ int(LDA + 1,8)
              OFFDAG = POSPV1+1_8
              A11 = A(POSPV1)
              A22 = A(POSPV2)
              A12 = A(OFFDAG)
              DETPIV = A11*A22 - A12**2
              A22 = A11/DETPIV
              A11 = A(POSPV2)/DETPIV
              A12 = -A12/DETPIV
              DO J = 1,M
                MULT1 = A11*LR_BLOCK_PTR(J,I)+A12*LR_BLOCK_PTR(J,I+1)
                MULT2 = A12*LR_BLOCK_PTR(J,I)+A22*LR_BLOCK_PTR(J,I+1)
                LR_BLOCK_PTR(J,I)   = MULT1
                LR_BLOCK_PTR(J,I+1) = MULT2
              ENDDO
              DPOS = POSPV2 + int(LDA + 1,8)
              I = I+2
            ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDIF
      CALL UPD_FLOP_TRSM(LRB%M, LRB%N, LRB%K, LRB%ISLR, LorU)
      END SUBROUTINE DMUMPS_LRTRSM
      SUBROUTINE DMUMPS_LRGEMM_SCALING(LRB, SCALED, A, LA, DIAG, 
     &          LD_DIAG, IW2, POSELTT, NFRONT, BLOCK, MAXI_CLUSTER) 
C This routine does the scaling (for the symmetric case) before 
C computing the LR product (done in DMUMPS_LRGEMM4)        
        TYPE(LRB_TYPE),INTENT(IN) :: LRB
        INTEGER(8), intent(in)  :: LA
        DOUBLE PRECISION, intent(inout)  :: A(LA)
        DOUBLE PRECISION, intent(inout), DIMENSION(:,:)  :: SCALED
        INTEGER,INTENT(IN) :: LD_DIAG, NFRONT, IW2(*)
        INTEGER(8), INTENT(IN) :: POSELTT
        DOUBLE PRECISION, INTENT(IN), OPTIONAL :: DIAG(*)
C       POSELTT = position of first entry in block to update with
C       gemm
        INTEGER, INTENT(IN) :: MAXI_CLUSTER
C       used as a workspace
C       Only column 1 of BLOCK is used as temporary space in routine
        DOUBLE PRECISION, intent(inout)  :: BLOCK(MAXI_CLUSTER)
        INTEGER :: J, NROWS
        DOUBLE PRECISION :: PIV1, PIV2, OFFDIAG
        IF (LRB%ISLR) THEN
            NROWS = LRB%K
        ELSE 
            NROWS = LRB%M
        ENDIF
        J = 1
        DO WHILE (J <= LRB%N)
C           1x1 pivot              
            IF (IW2(J) > 0) THEN
                SCALED(1:NROWS,J) = DIAG(1+LD_DIAG*(J-1)+J-1) 
     &           * SCALED(1:NROWS,J)
                J = J+1
            ELSE !2x2 pivot
                PIV1    = DIAG(1+LD_DIAG*(J-1)+J-1)
                PIV2    = DIAG(1+LD_DIAG*J+J)
C               off diagonal value of 2by2 pivot is stored 
C               in upper part                
                OFFDIAG = DIAG(1+LD_DIAG*(J-1)+J)
C               save J-th column of the block, before scaling
                BLOCK(1:NROWS)    = SCALED(1:NROWS,J)
                SCALED(1:NROWS,J) = PIV1 * SCALED(1:NROWS,J)
     &            + OFFDIAG * SCALED(1:NROWS,J+1)
                SCALED(1:NROWS,J+1) = OFFDIAG * BLOCK(1:NROWS)
     &            + PIV2 * SCALED(1:NROWS,J+1)
                 J=J+2
            ENDIF
        END DO
      END SUBROUTINE DMUMPS_LRGEMM_SCALING
      SUBROUTINE DMUMPS_LRGEMM4(ALPHA,
     &           LRB1, LRB2, BETA,
     &           A, LA, POSELTT, NFRONT, SYM, 
     &           IFLAG, IERROR,
     &           MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT,
     &           RANK, BUILDQ, 
     &           LUA_ACTIVATED, 
C Start of OPTIONAL arguments        
     &           LorU,
     &           LRB3, MAXI_RANK, 
     &           MAXI_CLUSTER,
     &           DIAG, LD_DIAG, IW2, BLOCK
     &           )
C
C
C Objective:
C ========= 
C    Compute product update of trailing matrix A, with 
C    LRB1 and LRB2 two blocks of LRB_TYPE structure:
C      A = BETA*A + ALPHA*LRB1*LRB2^T 
C
C 
        TYPE(LRB_TYPE),INTENT(IN) :: LRB1,LRB2
        INTEGER(8), intent(in)  :: LA
        DOUBLE PRECISION, intent(inout)  :: A(LA)
        INTEGER,INTENT(IN) :: NFRONT, SYM, TOL_OPT
        INTEGER,INTENT(INOUT) :: IFLAG, IERROR
C       MAXI_CLUSTER, MAXI_RANK: dimensions of the LRB3 allocated space
C       LRB3%Q is of size MAXI_CLUSTER x MAXI_RANK        
C       LRB3%R is of size MAXI_RANK x MAXI_CLUSTER        
C       The first LRB3%K columns of LRB3%Q/rows of LRB3%R are taken by
C       the accumulator. The rest are available for the new LRB.        
        INTEGER(8), INTENT(IN) :: POSELTT
        DOUBLE PRECISION, INTENT(IN), OPTIONAL :: DIAG(*)
        INTEGER,INTENT(IN), OPTIONAL :: LD_DIAG, IW2(*)
C       POSELTT = position of first entry in block to update with the
C       gemm
C       MIDBLK_COMPRESS = whether to compress the middle product
C       Y1*X2 in A = A + X1*(Y1*X2)*Y2 (strategy suggested by CA)
        INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT
C       Threshold used for RRQR if MIDBLK_COMPRESS >= 1 
        DOUBLE PRECISION, intent(in) :: TOLEPS
        DOUBLE PRECISION :: ALPHA,BETA
        LOGICAL, INTENT(OUT) :: BUILDQ
C       used as a workspace
C       Only column 1 of BLOCK is used as temporary
C       space in this routine + we declare BLOCK as
C       a contiguous array.
C       BLOCK could be declared of size MAXI_CLUSTER
        DOUBLE PRECISION, intent(inout), OPTIONAL  :: BLOCK(*)
        INTEGER, INTENT(IN), OPTIONAL :: LorU
        LOGICAL, INTENT(IN) :: LUA_ACTIVATED
C       MAXI_CLUSTER is an optional parameter required either for LUA
C       (to describe dimension of LRB array) or in the symmetric case
C       (to describe dimension of BLOCK array)
        INTEGER, INTENT(IN), OPTIONAL :: MAXI_CLUSTER
        INTEGER, INTENT(IN), OPTIONAL :: MAXI_RANK
        TYPE(LRB_TYPE), INTENT(INOUT), OPTIONAL :: LRB3
        DOUBLE PRECISION, POINTER, DIMENSION(:,:) :: XY_YZ
        DOUBLE PRECISION, ALLOCATABLE, TARGET, DIMENSION(:,:) :: XQ, R_Y
        DOUBLE PRECISION, POINTER, DIMENSION(:,:) :: X, Y, Y1, Y2, Z
        CHARACTER(len=1) :: SIDE, TRANSY
        INTEGER :: K_XY, K_YZ, LDY, LDY1, LDY2, K_Y
        INTEGER :: LDXY_YZ, SAVE_K
        INTEGER :: I, J, RANK, MAXRANK, INFO, LWORK
        DOUBLE PRECISION,    ALLOCATABLE :: RWORK_RRQR(:)
        DOUBLE PRECISION, ALLOCATABLE :: WORK_RRQR(:), TAU_RRQR(:), 
     &                          Y_RRQR(:,:)
        INTEGER, ALLOCATABLE :: JPVT_RRQR(:)
        LOGICAL, PARAMETER :: LRMIXED_ACTIVATED = .FALSE.
        INTEGER :: K_STOC(NPREC_STOC)
        INTEGER :: K453
        LOGICAL, PARAMETER :: MIXED_ADMISSIBILITY = .FALSE.
        INTEGER :: allocok, MREQ
        DOUBLE PRECISION :: ONE, MONE, ZERO
        PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
        PARAMETER (ZERO=0.0D0)
        IF (LRB1%M.EQ.0) THEN
          RETURN
        ENDIF
        IF (LRB2%M.EQ.0) THEN
        ENDIF
        RANK = 0
        BUILDQ = .FALSE.
        IF (LRB1%ISLR.AND.LRB2%ISLR) THEN
            IF ((LRB1%K.EQ.0).OR.(LRB2%K.EQ.0)) THEN
              GOTO 1200
            ENDIF
                allocate(Y(LRB1%K,LRB2%K),stat=allocok)
                IF (allocok > 0) THEN
                    MREQ = LRB1%K*LRB2%K
                    GOTO 1570
                ENDIF
            X    => LRB1%Q
            K_Y  =  LRB1%N
            IF (SYM .EQ. 0) THEN
                Y1  => LRB1%R
            ELSE
C               Symmetric case: column scaling of R1 is done
                allocate(Y1(LRB1%K,LRB1%N),stat=allocok)
                IF (allocok > 0) THEN
                  MREQ = LRB1%K*LRB1%N
                  GOTO 1570
                ENDIF
                DO J=1,LRB1%N
                    DO I=1,LRB1%K
                        Y1(I,J) = LRB1%R(I,J)
                    ENDDO
                ENDDO
                CALL DMUMPS_LRGEMM_SCALING(LRB1, Y1, A, LA, DIAG,
     &                 LD_DIAG, IW2, POSELTT, NFRONT, BLOCK, 
     &                 MAXI_CLUSTER) 
            ENDIF
            LDY1 =  LRB1%K
            Z    => LRB2%Q
            Y2   => LRB2%R
            LDY2 =  LRB2%K
             CALL dgemm('N', 'T', LRB1%K, LRB2%K, K_Y, ONE,
     &            Y1(1,1), LDY1, Y2(1,1), LDY2, ZERO, Y(1,1), LRB1%K )
            IF (MIDBLK_COMPRESS.GE.1) THEN 
C               Compress middle product Y
                LWORK = LRB2%K*(LRB2%K+1)
C               WORK_RRQR(1:LWORK) used for RRQR and dorgqr
                allocate(Y_RRQR(LRB1%K,LRB2%K),
     &               WORK_RRQR(LWORK), RWORK_RRQR(2*LRB2%K), 
     &               TAU_RRQR(MIN(LRB1%K,LRB2%K)),
     &               JPVT_RRQR(LRB2%K),stat=allocok)
                IF (allocok > 0) THEN
                  MREQ = LRB1%K*LRB2%K + LWORK + 2*LRB2%K +
     &                   MIN(LRB1%K,LRB2%K) + LRB2%K
                  GOTO 1570
                ENDIF
                DO J=1,LRB2%K
                    DO I=1,LRB1%K
                        Y_RRQR(I,J) = Y(I,J)
                    ENDDO
                ENDDO
                MAXRANK = MIN(LRB1%K, LRB2%K)-1
                MAXRANK = max (1, int((MAXRANK*KPERCENT/100)))
                JPVT_RRQR = 0
                CALL DMUMPS_TRUNCATED_RRQR(LRB1%K, LRB2%K, Y_RRQR(1,1),
     &               LRB1%K, JPVT_RRQR, TAU_RRQR, WORK_RRQR,
     &               LRB2%K, RWORK_RRQR, TOLEPS, TOL_OPT, RANK, 
     &               MAXRANK, INFO,
     &               K_STOC, LRMIXED_ACTIVATED, K453, 
     &               MIXED_ADMISSIBILITY,
     &               BUILDQ)
                IF (RANK.GT.MAXRANK) THEN 
                    deallocate(Y_RRQR, WORK_RRQR, RWORK_RRQR, TAU_RRQR,
     &                     JPVT_RRQR)
                    BUILDQ = .FALSE.
                ELSE
                    BUILDQ = .TRUE.
                ENDIF
                IF (BUILDQ) THEN
C                 Successfully compressed middle block
                  IF (RANK.EQ.0) THEN
                    deallocate(Y_RRQR, WORK_RRQR, RWORK_RRQR, TAU_RRQR, 
     &                         JPVT_RRQR) 
                      deallocate(Y)
                    nullify(Y)
C                   GOTO 1580 not ok because BUILDQ .EQV. true
C                   would try to free XQ and R_Y that are not allocated
C                   in that case. So we free Y1 now if it was allocated. 
                    IF (SYM .NE. 0) deallocate(Y1)
                    GOTO 1200
                  ELSE
                    allocate(XQ(LRB1%M,RANK), R_Y(RANK,LRB2%K),
     &                       stat=allocok)
                    IF (allocok > 0) THEN
                      MREQ = LRB1%M*RANK + RANK*LRB2%K
                      GOTO 1570
                    ENDIF
                    DO J=1, LRB2%K
                       R_Y(1:MIN(RANK,J),JPVT_RRQR(J)) =
     &                   Y_RRQR(1:MIN(RANK,J),J)
                       IF(J.LT.RANK) R_Y(MIN(RANK,J)+1:
     &                   RANK,JPVT_RRQR(J))= ZERO
                    END DO
C                   LWORK=LRB2%K*(LRB2%K+1), with LRB2%K>RANK
C                   large enough for dorgqr
C                   dorgqr+dgemm is faster than dormqr (less flops and
C                   according to quick test, less time)                    
                    CALL dorgqr 
     &                  (LRB1%K, RANK, RANK, Y_RRQR(1,1),
     &                  LRB1%K, TAU_RRQR(1),  
     &                  WORK_RRQR(1), LWORK, INFO )
                    CALL dgemm('N', 'N', LRB1%M, RANK, LRB1%K, ONE,
     &                    X(1,1), LRB1%M, Y_RRQR(1,1), LRB1%K, ZERO, 
     &                    XQ(1,1), LRB1%M)
                    deallocate(Y_RRQR, WORK_RRQR, RWORK_RRQR, TAU_RRQR, 
     &                         JPVT_RRQR)   
                    nullify(X)
                    X      => XQ
                    K_XY   =  RANK
                      deallocate(Y)
                    nullify(Y)
                    Y      => R_Y
                    LDY    =  RANK
                    K_YZ   =  LRB2%K
                    TRANSY =  'N'
                    SIDE   = 'R'
                  ENDIF
                ENDIF
            ENDIF
            IF (.NOT.BUILDQ) THEN 
C               Did not try or succeed to compress middle product Y
                  LDY    = LRB1%K
                K_XY   = LRB1%K
                K_YZ   = LRB2%K
                TRANSY = 'N'
C  FIXME Actually the formula to decide which side to multiply to is a
C  bit more complex, although tests suggest it doesn't make a big
C  difference (see thesis for true formula)
                IF (LRB1%K .GE. LRB2%K) THEN
                    SIDE = 'L'
                ELSE 
                    SIDE = 'R'
                ENDIF
            ENDIF
        ENDIF
        IF (LRB1%ISLR.AND.(.NOT.LRB2%ISLR)) THEN 
            IF (LRB1%K.EQ.0) THEN
              GOTO 1200
            ENDIF
            SIDE   =  'R'
            K_XY   =  LRB1%K
            TRANSY =  'N' 
            Z      => LRB2%Q 
            X   => LRB1%Q
            LDY =  LRB1%K
            IF (SYM .EQ. 0) THEN
              Y   => LRB1%R
            ELSE
C             Symmetric case: column scaling of R1 is done
              allocate(Y(LRB1%K,LRB1%N),stat=allocok)
              IF (allocok > 0) THEN
                MREQ = LRB1%K*LRB1%N
                GOTO 1570
              ENDIF
              DO J=1,LRB1%N
                  DO I=1,LRB1%K
                      Y(I,J) = LRB1%R(I,J)
                  ENDDO
              ENDDO
              CALL DMUMPS_LRGEMM_SCALING(LRB1, Y, A, LA, DIAG,
     &               LD_DIAG, IW2, POSELTT, NFRONT, BLOCK, 
     &               MAXI_CLUSTER) 
            ENDIF
            K_YZ = LRB2%N
        ENDIF
        IF ((.NOT.LRB1%ISLR).AND.LRB2%ISLR) THEN 
            IF (LRB2%K.EQ.0) THEN
              GOTO 1200
            ENDIF
            SIDE   =  'L'
            K_YZ   =  LRB2%K
            X      => LRB1%Q 
            TRANSY =  'T' 
            K_XY = LRB1%N
            IF (SYM .EQ. 0) THEN
                Y  => LRB2%R
            ELSE
C               Symmetric case: column scaling of R2 is done
                allocate(Y(LRB2%K,LRB2%N),stat=allocok)
                IF (allocok > 0) THEN
                  MREQ = LRB2%K*LRB2%N
                  GOTO 1570
                ENDIF
                DO J=1,LRB2%N
                    DO I=1,LRB2%K
                        Y(I,J) = LRB2%R(I,J)
                    ENDDO
                ENDDO
                CALL DMUMPS_LRGEMM_SCALING(LRB2, Y, A, LA, DIAG,
     &                 LD_DIAG, IW2, POSELTT, NFRONT, BLOCK, 
     &                 MAXI_CLUSTER) 
            ENDIF
            LDY =  LRB2%K
            Z   => LRB2%Q
        ENDIF
        IF ((.NOT.LRB1%ISLR).AND.(.NOT.LRB2%ISLR)) THEN 
            IF (SYM .EQ. 0) THEN
                X => LRB1%Q
            ELSE
C               Symmetric case: column scaling of Q1 done
                allocate(X(LRB1%M,LRB1%N),stat=allocok)
                IF (allocok > 0) THEN
                  MREQ = LRB1%M*LRB1%N
                  GOTO 1570
                ENDIF
C               X = LRB1%Q segfaults  
                DO J=1,LRB1%N
                    DO I=1,LRB1%M
                        X(I,J) = LRB1%Q(I,J)
                    ENDDO
                ENDDO
                CALL DMUMPS_LRGEMM_SCALING(LRB1, X, A, LA, DIAG,
     &                   LD_DIAG, IW2, POSELTT, NFRONT, BLOCK, 
     &                   MAXI_CLUSTER) 
            ENDIF
            SIDE   =  'N'
            Z      => LRB2%Q
            K_XY = LRB1%N  
        ENDIF
        IF (LUA_ACTIVATED) THEN
          SAVE_K = LRB3%K
          IF (SIDE == 'L') THEN
            LRB3%K = LRB3%K+K_YZ
          ELSEIF (SIDE == 'R') THEN 
            LRB3%K = LRB3%K+K_XY
          ENDIF
        ENDIF
        IF (SIDE == 'L') THEN ! LEFT: XY_YZ = X*Y; A = XY_YZ*Z
            IF (.NOT.LUA_ACTIVATED) THEN
C             No LUA: regular LRGEMM
              allocate(XY_YZ(LRB1%M,K_YZ),stat=allocok)
              IF (allocok > 0) THEN
                MREQ = LRB1%M*K_YZ
                GOTO 1570
              ENDIF
              LDXY_YZ = LRB1%M
            ELSE
C             LUA activated: result kept as LR block LRB3
               IF (SAVE_K+K_YZ.GT.MAXI_RANK) THEN
                write(*,*) 'Internal error in DMUMPS_LRGEMM4 1a',
     &           'K_ACC+K_CUR>K_MAX:',SAVE_K,K_YZ,MAXI_RANK
                CALL MUMPS_ABORT()
              ENDIF
              IF (LRB3%M.NE.LRB1%M) THEN
                write(*,*) 'Internal error in DMUMPS_LRGEMM4 1b',
     &           'LRB1%M =/= LRB3%M',LRB1%M,LRB3%M
                CALL MUMPS_ABORT()
              ENDIF
              XY_YZ => LRB3%Q(1:LRB1%M,SAVE_K+1:SAVE_K+K_YZ)
              LDXY_YZ = MAXI_CLUSTER
C             Copy Z in LRB3%R 
              DO I=1,K_YZ
                LRB3%R(SAVE_K+I,1:LRB2%M) = Z(1:LRB2%M,I)
              ENDDO
            ENDIF
            CALL dgemm('N', TRANSY, LRB1%M, K_YZ, K_XY, ONE,
     &             X(1,1), LRB1%M, Y(1,1), LDY, ZERO, XY_YZ(1,1), 
     &             LDXY_YZ)
            IF (.NOT.LUA_ACTIVATED) THEN
              CALL dgemm('N', 'T', LRB1%M, LRB2%M, K_YZ, ALPHA,
     &               XY_YZ(1,1), LRB1%M, Z(1,1), LRB2%M, BETA, 
     &               A(POSELTT), NFRONT)
              deallocate(XY_YZ)
            ENDIF
        ELSEIF (SIDE == 'R') THEN ! RIGHT: XY_YZ = Y*Z; A = X*XY_YZ
            IF (.NOT.LUA_ACTIVATED) THEN
C             No LUA: regular LRGEMM
              allocate(XY_YZ(K_XY,LRB2%M),stat=allocok)
              IF (allocok > 0) THEN
                MREQ = K_XY*LRB2%M
                GOTO 1570
              ENDIF
              LDXY_YZ = K_XY
            ELSE
C             LUA activated: result kept as LR block LRB3
               IF (SAVE_K+K_XY.GT.MAXI_RANK) THEN
                write(*,*) 'Internal error in DMUMPS_LRGEMM4 2a',
     &           'K_ACC+K_CUR>K_MAX:',SAVE_K,K_XY,MAXI_RANK
                CALL MUMPS_ABORT()
              ENDIF
              IF (LRB3%N.NE.LRB2%M) THEN
                write(*,*) 'Internal error in DMUMPS_LRGEMM4 2b',
     &           'LRB2%M =/= LRB3%N',LRB2%M,LRB3%N
                CALL MUMPS_ABORT()
              ENDIF
              XY_YZ => LRB3%R(SAVE_K+1:SAVE_K+K_XY,1:LRB2%M)
              LDXY_YZ = MAXI_RANK
C             Copy X in LRB3%Q 
              DO I=1,K_XY
                LRB3%Q(1:LRB1%M,SAVE_K+I) = X(1:LRB1%M,I)
              ENDDO
            ENDIF
            CALL dgemm(TRANSY, 'T', K_XY, LRB2%M, K_YZ, ONE,
     &             Y(1,1), LDY, Z(1,1), LRB2%M, ZERO, XY_YZ(1,1), 
     &             LDXY_YZ)
            IF (.NOT.LUA_ACTIVATED) THEN
              CALL dgemm('N', 'N', LRB1%M, LRB2%M, K_XY, ALPHA,
     &               X(1,1), LRB1%M, XY_YZ(1,1), K_XY, BETA, A(POSELTT),
     &               NFRONT)
              deallocate(XY_YZ)
            ENDIF
        ELSE ! SIDE == 'N' : NONE; A = X*Z
          CALL dgemm('N', 'T', LRB1%M, LRB2%M, K_XY, ALPHA,
     &             X(1,1), LRB1%M, Z(1,1), LRB2%M, BETA, A(POSELTT),
     &             NFRONT)
        ENDIF
        GOTO 1580
 1570 CONTINUE        
C       Alloc NOT ok!!        
        IFLAG  = -13
        IERROR = MREQ
        RETURN
 1580 CONTINUE       
C       Alloc ok!!        
        IF ((.NOT.LRB1%ISLR).AND.(.NOT.LRB2%ISLR)) THEN 
            IF (SYM .NE. 0) deallocate(X)
        ELSEIF ((.NOT.LRB1%ISLR).AND.LRB2%ISLR) THEN 
            IF (SYM .NE. 0) deallocate(Y)
        ELSEIF (LRB1%ISLR.AND.(.NOT.LRB2%ISLR)) THEN 
            IF (SYM .NE. 0) deallocate(Y)
        ELSE
            IF (SYM .NE. 0) deallocate(Y1)
            IF ((MIDBLK_COMPRESS.GE.1).AND.BUILDQ) THEN
                deallocate(XQ)
                deallocate(R_Y)
            ELSE
                  deallocate(Y)
            ENDIF
        ENDIF
 1200 CONTINUE       
      END SUBROUTINE DMUMPS_LRGEMM4
      SUBROUTINE DMUMPS_DECOMPRESS_ACC(ACC_LRB, MAXI_CLUSTER, 
     &           MAXI_RANK, A, LA, POSELTT, NFRONT, NIV, LorU, 
     &           COUNT_FLOPS)
C
C LUA related routine
C This subroutine decompresses the LR accumulator and sums it to A
C 
        TYPE(LRB_TYPE),INTENT(INOUT) :: ACC_LRB
        INTEGER(8), intent(in)  :: LA
        DOUBLE PRECISION, intent(inout)  :: A(LA)
        INTEGER,INTENT(IN) :: NFRONT, NIV
C       MAXI_CLUSTER, MAXI_RANK are the dimensions of the ACC_LRB
C       allocated space.     
C       ACC_LRB%Q is of size MAXI_CLUSTER x MAXI_RANK        
C       ACC_LRB%R is of size MAXI_RANK x MAXI_CLUSTER        
C       The first ACC_LRB%K columns of ACC_LRB%Q/rows of ACC_LRB%R are
C       taken by the accumulator. The rest are available for the new
C       LRB.        
        INTEGER,INTENT(IN) :: MAXI_CLUSTER, MAXI_RANK
C       POSELTT = position of first entry in block to update with the
C       gemm
        INTEGER(8), INTENT(IN) :: POSELTT
C       Whether to count the flops in the total
        LOGICAL, OPTIONAL :: COUNT_FLOPS
        LOGICAL  :: COUNT_FLOPS_LOC
        INTEGER  :: LorU
        DOUBLE PRECISION :: ONE, MONE, ZERO
        PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
        PARAMETER (ZERO=0.0D0)
        IF (present(COUNT_FLOPS)) THEN
          COUNT_FLOPS_LOC=COUNT_FLOPS
        ELSE
          COUNT_FLOPS_LOC=.TRUE.
        ENDIF
        CALL dgemm('N', 'N', ACC_LRB%M, ACC_LRB%N, ACC_LRB%K,
     &         MONE, ACC_LRB%Q(1,1), MAXI_CLUSTER, ACC_LRB%R(1,1), 
     &         MAXI_RANK, ONE, A(POSELTT), NFRONT)
        ACC_LRB%K = 0
      END SUBROUTINE DMUMPS_DECOMPRESS_ACC
      SUBROUTINE DMUMPS_COMPRESS_FR_UPDATES(ACC_LRB, MAXI_CLUSTER, 
     &          MAXI_RANK, A, LA, POSELTT, NFRONT, NIV,
     &          TOLEPS, TOL_OPT, KPERCENT, BUILDQ, LorU, CB_COMPRESS)
C 
C  LUA related routine
C  
        TYPE(LRB_TYPE),INTENT(INOUT) :: ACC_LRB
        INTEGER(8), intent(in)  :: LA
        DOUBLE PRECISION, intent(inout)  :: A(LA)
        INTEGER,INTENT(IN) :: NFRONT, NIV, LorU, TOL_OPT
C       MAXI_CLUSTER, MAXI_RANK are the dimensions of the ACC_LRB
C       allocated space.     
C       ACC_LRB%Q is of size MAXI_CLUSTER x MAXI_RANK        
C       ACC_LRB%R is of size MAXI_RANK x MAXI_CLUSTER        
C       The first ACC_LRB%K columns of ACC_LRB%Q/rows of ACC_LRB%R are
C       taken by the accumulator. 
        INTEGER,INTENT(IN) :: MAXI_CLUSTER, MAXI_RANK, KPERCENT
        INTEGER(8), INTENT(IN) :: POSELTT
C       POSELTT = position of first entry in block to update with the
C       gemm
C       Threshold used for RRQR
        DOUBLE PRECISION, intent(in) :: TOLEPS
C       BUILDQ = whether we successfully compressed the FR updates        
        LOGICAL, INTENT(OUT) :: BUILDQ
        LOGICAL, INTENT(IN) :: CB_COMPRESS
        DOUBLE PRECISION,    ALLOCATABLE :: RWORK_RRQR(:)
        DOUBLE PRECISION, ALLOCATABLE :: WORK_RRQR(:), TAU_RRQR(:)
        INTEGER, ALLOCATABLE :: JPVT_RRQR(:)
        LOGICAL, PARAMETER :: LRMIXED_ACTIVATED = .FALSE.
        INTEGER :: K_STOC(NPREC_STOC)
        INTEGER :: K453
        LOGICAL, PARAMETER :: MIXED_ADMISSIBILITY = .FALSE.
        INTEGER :: INFO, RANK, MAXRANK, LWORK
        INTEGER :: I, J, M, N
        INTEGER :: allocok, MREQ
        DOUBLE PRECISION :: ONE, MONE, ZERO
        PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
        PARAMETER (ZERO=0.0D0)
              M = ACC_LRB%M
              N = ACC_LRB%N
              MAXRANK = floor(dble(M*N)/dble(M+N))
              MAXRANK = max (1, int((MAXRANK*KPERCENT/100)))
              LWORK = N*(N+1)
              allocate(WORK_RRQR(LWORK), RWORK_RRQR(2*N), 
     &             TAU_RRQR(N),
     &             JPVT_RRQR(N), stat=allocok)
              IF (allocok > 0) THEN
                MREQ = LWORK +4 *N
                GOTO 100
              ENDIF
C             Note the minus sign (because we actually subtract to A in
C             DECOMPRESS_ACC)              
C             Also note we compress A^T as the accumulator is a LR
C             representation of A^T. Cf ALLOC_LRB_FROM_ACC where we
C             retranspose the accumulator to have a LR representation of
C             A to keep
                DO I=1,N
                  ACC_LRB%Q(1:M,I)= 
     &            - A(POSELTT+int(I-1,8)*int(NFRONT,8) :
     &            POSELTT+int(I-1,8)*int(NFRONT,8) + int(M-1,8) )
                END DO  
              JPVT_RRQR = 0
              CALL DMUMPS_TRUNCATED_RRQR(M, N, ACC_LRB%Q(1,1),
     &               MAXI_CLUSTER, JPVT_RRQR(1), TAU_RRQR(1), 
     &               WORK_RRQR(1),
     &               N, RWORK_RRQR(1), TOLEPS, TOL_OPT, 
     &               RANK, MAXRANK, INFO,
     &               K_STOC, LRMIXED_ACTIVATED, K453, 
     &               MIXED_ADMISSIBILITY,
     &               BUILDQ)
              IF (BUILDQ) THEN
C               Successfully compressed ACC_LRB%Q
                DO J=1, N
                   ACC_LRB%R(1:MIN(RANK,J),JPVT_RRQR(J)) =
     &               ACC_LRB%Q(1:MIN(RANK,J),J)
                   IF(J.LT.RANK) ACC_LRB%R(MIN(RANK,J)+1:
     &               RANK,JPVT_RRQR(J))= ZERO
                END DO
                CALL dorgqr 
     &              (M, RANK, RANK, ACC_LRB%Q(1,1),
     &              MAXI_CLUSTER, TAU_RRQR(1),  
     &              WORK_RRQR(1), LWORK, INFO )
C               Set FR array to zero since it is now in the LR array                
                DO I=1,N
                  A( POSELTT+int(I-1,8)*int(NFRONT,8) :
     &              POSELTT+int(I-1,8)*int(NFRONT,8)+int(M-1,8) ) = ZERO
                END DO  
                ACC_LRB%K = RANK
                CALL UPD_FLOP_COMPRESS(ACC_LRB%M, ACC_LRB%N, ACC_LRB%K, 
     &                       ACC_LRB%ISLR, CB_COMPRESS=CB_COMPRESS)
              ELSE
C               Just to count flops
                ACC_LRB%K = RANK
                ACC_LRB%ISLR = .FALSE.
                CALL UPD_FLOP_COMPRESS(ACC_LRB%M, ACC_LRB%N, ACC_LRB%K, 
     &                       ACC_LRB%ISLR, CB_COMPRESS=CB_COMPRESS)
                ACC_LRB%ISLR = .TRUE.
                ACC_LRB%K = 0
              ENDIF
              deallocate(JPVT_RRQR, TAU_RRQR, WORK_RRQR, RWORK_RRQR)
              RETURN
 100        CONTINUE        
C     Alloc NOT ok!!
            write(*,*) 'Allocation problem in BLR routine 
     &        DMUMPS_COMPRESS_FR_UPDATES: ',
     &        'not enough memory? memory requested = ' , MREQ
            CALL MUMPS_ABORT()
            RETURN
      END SUBROUTINE DMUMPS_COMPRESS_FR_UPDATES
      SUBROUTINE DMUMPS_RECOMPRESS_ACC(ACC_LRB, MAXI_CLUSTER, 
     &          MAXI_RANK, A, LA, POSELTT, NFRONT, NIV,
     &          MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT_RMB, 
     &          KPERCENT_LUA, NEW_ACC_RANK)
C 
C  LUA related routine
C  This subroutine recompresses the LR accumulator ACC_LRB
C   1) compute Q1,R1,Q2,R2 such that LRB%Q = Q1*R1 and (LRB%R)^T = Q2*R2        
C   2) perform LRB product Q1*(R1*R2^T)*Q2^T and store result in ACC_LRB        
C  
        TYPE(LRB_TYPE),INTENT(INOUT) :: ACC_LRB
        INTEGER(8), intent(in)  :: LA
        DOUBLE PRECISION, intent(inout)  :: A(LA)
        INTEGER,INTENT(IN) :: NFRONT, NIV, TOL_OPT
        INTEGER :: IFLAG, IERROR
C       MAXI_CLUSTER, MAXI_RANK are the dimensions of the ACC_LRB
C       allocated space.     
C       ACC_LRB%Q is of size MAXI_CLUSTER x MAXI_RANK        
C       ACC_LRB%R is of size MAXI_RANK x MAXI_CLUSTER        
C       The first ACC_LRB%K columns of ACC_LRB%Q/rows of ACC_LRB%R are
C       taken by the accumulator. 
        INTEGER,INTENT(IN) :: MAXI_CLUSTER, MAXI_RANK, KPERCENT_LUA
C       NEW_ACC_RANK is the size of the non-compressed part of the
C       accumulator. One needs to orthogonalize this part against
C       the rest of the accumulator. 
        INTEGER,INTENT(INOUT) :: NEW_ACC_RANK
        INTEGER(8), INTENT(IN) :: POSELTT
C       POSELTT = position of first entry in block to update with the
C       gemm
C       MIDBLK_COMPRESS = whether to compress the middle product
C       Y1*X2 in A = A + X1*(Y1*X2)*Y2 (strategy suggested by CA)
        INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT_RMB
C       Threshold used for RRQR if MIDBLK_COMPRESS >= 2 
        DOUBLE PRECISION, intent(in) :: TOLEPS
        DOUBLE PRECISION,    ALLOCATABLE:: RWORK_RRQR(:)
        DOUBLE PRECISION, ALLOCATABLE:: WORK_RRQR(:), TAU_RRQR(:)
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:),TARGET:: Q1, R1,
     &                                                Q2, R2
        INTEGER, ALLOCATABLE :: JPVT_RRQR(:)
        TYPE(LRB_TYPE)       :: LRB1, LRB2
        INTEGER :: INFO, RANK1, RANK2, RANK, MAXRANK, LWORK
        LOGICAL :: BUILDQ, BUILDQ1, BUILDQ2, SKIP1, SKIP2
        LOGICAL, PARAMETER :: LRMIXED_ACTIVATED = .FALSE.
        INTEGER :: K_STOC(NPREC_STOC)
        INTEGER :: K453
        LOGICAL, PARAMETER :: MIXED_ADMISSIBILITY = .FALSE.
        INTEGER :: I, J, M, N, K
        INTEGER :: allocok, MREQ
        DOUBLE PRECISION :: ONE, MONE, ZERO
        PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
        PARAMETER (ZERO=0.0D0)
              SKIP1 = .FALSE.
              SKIP2 = .FALSE.
              SKIP1 = .TRUE.
 1500 CONTINUE              
              M = ACC_LRB%M
              N = ACC_LRB%N
              K = ACC_LRB%K
              MAXRANK = K-1
              MAXRANK = max (1, int((MAXRANK*KPERCENT_LUA/100)))
              LWORK = K*(K+1)
              IF (.FALSE.) THEN
                CALL DMUMPS_RECOMPRESS_ACC_V2(ACC_LRB,
     &                 MAXI_CLUSTER, MAXI_RANK, A, LA, POSELTT,
     &                 NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS,
     &                 TOL_OPT, KPERCENT_RMB, KPERCENT_LUA, 
     &                 NEW_ACC_RANK)
                K = ACC_LRB%K
                MAXRANK = K-1
                MAXRANK = max (1, int((MAXRANK*KPERCENT_LUA/100)))
                LWORK = K*(K+1)
                SKIP1 = .TRUE.
                SKIP2 = K.EQ.0
              ENDIF
              IF (SKIP1.AND.SKIP2) GOTO 1600
              allocate(Q1(M,K), Q2(N,K), 
     &             WORK_RRQR(LWORK), 
     &             RWORK_RRQR(2*K), 
     &             TAU_RRQR(K),
     &             JPVT_RRQR(K), stat=allocok)
              IF (allocok > 0) THEN
                MREQ = LWORK + M*N + N*K+ 4 * K
                GOTO 100
              ENDIF
              IF (SKIP1) THEN
                BUILDQ1 = .FALSE.
              ELSE
                DO J=1,K
                  DO I=1,M
                    Q1(I,J) = ACC_LRB%Q(I,J)
                  ENDDO
                ENDDO
                JPVT_RRQR = 0
                CALL DMUMPS_TRUNCATED_RRQR(M, K, Q1(1,1),
     &               M, JPVT_RRQR(1), TAU_RRQR(1), WORK_RRQR(1),
     &               K, RWORK_RRQR(1), TOLEPS, TOL_OPT, RANK1, 
     &               MAXRANK, INFO,
     &               K_STOC, LRMIXED_ACTIVATED, K453, 
     &               MIXED_ADMISSIBILITY,
     &               BUILDQ1)
              ENDIF
              IF (BUILDQ1) THEN
C                 Successfully compressed Q1
                  allocate(R1(RANK1,K), stat=allocok)
                  IF (allocok > 0) THEN
                     MREQ = RANK1*K
                     GOTO 100
                  ENDIF
                  DO J=1, K
                     R1(1:MIN(RANK1,J),JPVT_RRQR(J)) =
     &                 Q1(1:MIN(RANK1,J),J)
                     IF(J.LT.RANK1) R1(MIN(RANK1,J)+1:
     &                 RANK1,JPVT_RRQR(J))= ZERO
                  END DO
                  CALL dorgqr 
     &                (M, RANK1, RANK1, Q1(1,1),
     &                M, TAU_RRQR(1),  
     &                WORK_RRQR(1), LWORK, INFO )
              ENDIF
              IF (SKIP2) THEN
                BUILDQ2 = .FALSE.
              ELSE
                DO J=1,K
                  DO I=1,N
                    Q2(I,J) = ACC_LRB%R(J,I)
                  ENDDO
                ENDDO
                JPVT_RRQR = 0
                CALL DMUMPS_TRUNCATED_RRQR(N, K, Q2(1,1),
     &               N, JPVT_RRQR(1), TAU_RRQR(1), WORK_RRQR(1),
     &               K, RWORK_RRQR(1), TOLEPS, TOL_OPT, 
     &               RANK2, MAXRANK, INFO,
     &               K_STOC, LRMIXED_ACTIVATED, K453, 
     &               MIXED_ADMISSIBILITY,
     &               BUILDQ2)
              ENDIF
              IF (BUILDQ2) THEN
C                 Successfully compressed Q2
                  allocate(R2(RANK2,K), stat=allocok)
                  IF (allocok > 0) THEN
                     MREQ = RANK2*K
                     GOTO 100
                  ENDIF
                  DO J=1, K
                     R2(1:MIN(RANK2,J),JPVT_RRQR(J)) =
     &                 Q2(1:MIN(RANK2,J),J)
                     IF(J.LT.RANK2) R2(MIN(RANK2,J)+1:
     &                 RANK2,JPVT_RRQR(J))= ZERO
                  END DO
                  CALL dorgqr 
     &                (N, RANK2, RANK2, Q2(1,1),
     &                N, TAU_RRQR(1),  
     &                WORK_RRQR(1), LWORK, INFO )
              ENDIF
              CALL INIT_LRB(LRB1,RANK1,M,K,BUILDQ1)
              CALL INIT_LRB(LRB2,RANK2,N,K,BUILDQ2)
              IF (BUILDQ1.OR.BUILDQ2) THEN
                IF (BUILDQ1) THEN
                  LRB1%R => R1
                ELSE
                  DO J=1,K
                    DO I=1,M
                      Q1(I,J) = ACC_LRB%Q(I,J)
                    ENDDO
                  ENDDO
                ENDIF
                LRB1%Q => Q1
                IF (BUILDQ2) THEN
                  LRB2%R => R2
                ELSE
                  DO J=1,K
                    DO I=1,N
                      Q2(I,J) = ACC_LRB%R(J,I)
                    ENDDO
                  ENDDO
                ENDIF
                LRB2%Q => Q2
                ACC_LRB%K = 0
                CALL DMUMPS_LRGEMM4(MONE, LRB1, LRB2, ONE,
     &             A, LA, POSELTT, NFRONT, 0, IFLAG, IERROR,
     &             MIDBLK_COMPRESS-1, TOLEPS, TOL_OPT, 
     &             KPERCENT_RMB, 
     &             RANK, BUILDQ, .TRUE., LRB3=ACC_LRB,
     &             MAXI_RANK=MAXI_RANK, MAXI_CLUSTER=MAXI_CLUSTER)
                IF (IFLAG.LT.0) GOTO 100
                CALL UPD_FLOP_UPDATE(LRB1%M, LRB1%N, LRB1%K, LRB1%ISLR,
     &                  LRB2%M, LRB2%N, LRB2%K, LRB2%ISLR,
     &                  MIDBLK_COMPRESS-1, RANK, BUILDQ,
     &                  .TRUE., .FALSE., REC_ACC=.TRUE.) 
              ENDIF
              IF (.NOT. SKIP1)
     &          CALL UPD_FLOP_COMPRESS(LRB1%M, LRB1%N, LRB1%K, 
     &                    LRB1%ISLR, REC_ACC=.TRUE.)
              IF (.NOT. SKIP2)
     &          CALL UPD_FLOP_COMPRESS(LRB2%M, LRB2%N, LRB2%K, 
     &                    LRB2%ISLR, REC_ACC=.TRUE.)
              deallocate(Q1,Q2)
              IF (BUILDQ1) deallocate(R1)
              IF (BUILDQ2) deallocate(R2)
              deallocate(JPVT_RRQR, TAU_RRQR, WORK_RRQR, RWORK_RRQR)
              IF (SKIP1.AND.(RANK2.GT.0)) THEN
                SKIP1 = .FALSE.
                SKIP2 = .TRUE.
                GOTO 1500
              ENDIF
 1600         CONTINUE              
              NEW_ACC_RANK = 0
              RETURN
 100          CONTINUE
C     Alloc NOT ok!!
              write(*,*) 'Allocation problem in BLR routine 
     &          DMUMPS_RECOMPRESS_ACC: ',
     &          'not enough memory? memory requested = ' , MREQ        
              CALL MUMPS_ABORT()             
              RETURN
      END SUBROUTINE DMUMPS_RECOMPRESS_ACC
      RECURSIVE SUBROUTINE DMUMPS_RECOMPRESS_ACC_NARYTREE(
     &          ACC_LRB, MAXI_CLUSTER, MAXI_RANK, A, LA, POSELTT, KEEP8,
     &          NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT, 
     &          KPERCENT_RMB,
     &          KPERCENT_LUA, K478, RANK_LIST, POS_LIST, NB_NODES,
     &          LEVEL, ACC_TMP)
C 
C  LUA related routine
C  This subroutine recompresses the LR accumulator ACC_LRB using a n-ary
C  tree as recompression graph
        TYPE(LRB_TYPE),TARGET,INTENT(INOUT) :: ACC_LRB
        TYPE(LRB_TYPE),TARGET,INTENT(INOUT),OPTIONAL :: ACC_TMP
        INTEGER(8), intent(in)  :: LA
        DOUBLE PRECISION, intent(inout)  :: A(LA)
        INTEGER,INTENT(IN) :: NFRONT, NIV, TOL_OPT
C       MAXI_CLUSTER, MAXI_RANK are the dimensions of the ACC_LRB
C       allocated space.     
C       ACC_LRB%Q is of size MAXI_CLUSTER x MAXI_RANK        
C       ACC_LRB%R is of size MAXI_RANK x MAXI_CLUSTER        
C       The first ACC_LRB%K columns of ACC_LRB%Q/rows of ACC_LRB%R are
C       taken by the accumulator. 
        INTEGER,INTENT(IN) :: MAXI_CLUSTER, MAXI_RANK, KPERCENT_LUA
C       POSELTT = position of first entry in block to update with the
C       gemm
        INTEGER(8), INTENT(IN) :: POSELTT
        INTEGER(8) :: KEEP8(150)
C       MIDBLK_COMPRESS = whether to compress the middle product
C       Y1*X2 in A = A + X1*(Y1*X2)*Y2 (strategy suggested by CA)
        INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT_RMB
C       Threshold used for RRQR if MIDBLK_COMPRESS >= 2 
        DOUBLE PRECISION, intent(in) :: TOLEPS
C       KEEP(478)=-n => use a n-ary tree 
        INTEGER,INTENT(IN) :: K478, NB_NODES, LEVEL
        INTEGER,INTENT(INOUT) :: RANK_LIST(NB_NODES), POS_LIST(NB_NODES)
        TYPE(LRB_TYPE)       :: LRB, ACC_NEW
        TYPE(LRB_TYPE), POINTER :: LRB_PTR
        LOGICAL :: RESORT
        INTEGER :: I, J, M, N, L, NODE_RANK, NARY, IOFF, IMAX, CURPOS
        INTEGER :: NB_NODES_NEW, KTOT, NEW_ACC_RANK
        INTEGER, ALLOCATABLE :: RANK_LIST_NEW(:), POS_LIST_NEW(:)
        DOUBLE PRECISION :: ONE, MONE, ZERO
        PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
        PARAMETER (ZERO=0.0D0)
        INTEGER :: allocok
        RESORT = .FALSE.
        M = ACC_LRB%M
        N = ACC_LRB%N
        NARY = -K478
        IOFF = 0
C       Number of nodes on the next level of the tree        
        NB_NODES_NEW = NB_NODES/NARY
        IF (NB_NODES_NEW*NARY.NE.NB_NODES) THEN
          NB_NODES_NEW = NB_NODES_NEW + 1
        ENDIF
C       Ranks and positions of the nodes on the next level of the tree       
        ALLOCATE(RANK_LIST_NEW(NB_NODES_NEW),POS_LIST_NEW(NB_NODES_NEW),
     &       stat=allocok)
        IF (allocok > 0) THEN
           write(*,*) 'Allocation error of RANK_LIST_NEW/POS_LIST_NEW ',
     &      'in DMUMPS_RECOMPRESS_ACC_NARYTREE'
           call MUMPS_ABORT()
        ENDIF
        DO J=1,NB_NODES_NEW
C         For each node of the next level, we first merge its children
C         (i.e., make sure they are in adjacent memory area) and then
C         recompress them with DMUMPS_RECOMPRESS_ACC.       
          NODE_RANK = RANK_LIST(IOFF+1)
          CURPOS = POS_LIST(IOFF+1)
          IMAX = MIN(NARY,NB_NODES-IOFF)
          IF (IMAX.GE.2) THEN
C           More than one child, we merge and recompress them            
            DO I=2,IMAX
              IF (POS_LIST(IOFF+I).NE.CURPOS+NODE_RANK) THEN
C               Merge (memory copy necessary)
                DO L=0,RANK_LIST(IOFF+I)-1
                  ACC_LRB%Q(1:M,CURPOS+NODE_RANK+L) = 
     &               ACC_LRB%Q(1:M,POS_LIST(IOFF+I)+L)
                  ACC_LRB%R(CURPOS+NODE_RANK+L,1:N) = 
     &               ACC_LRB%R(POS_LIST(IOFF+I)+L,1:N)
                ENDDO
                POS_LIST(IOFF+I) = CURPOS+NODE_RANK
              ENDIF
              NODE_RANK = NODE_RANK+RANK_LIST(IOFF+I)
            ENDDO
            CALL INIT_LRB(LRB,NODE_RANK,M,N,.TRUE.)
            IF (.NOT.RESORT.OR.LEVEL.EQ.0) THEN
              LRB%Q => ACC_LRB%Q(1:M,CURPOS:CURPOS+NODE_RANK)
              LRB%R => ACC_LRB%R(CURPOS:CURPOS+NODE_RANK,1:N)
            ELSE
              LRB%Q => ACC_TMP%Q(1:M,CURPOS:CURPOS+NODE_RANK)
              LRB%R => ACC_TMP%R(CURPOS:CURPOS+NODE_RANK,1:N)
            ENDIF
            NEW_ACC_RANK = NODE_RANK-RANK_LIST(IOFF+1)
            IF (NEW_ACC_RANK.GT.0) THEN
              CALL DMUMPS_RECOMPRESS_ACC(LRB,
     &             MAXI_CLUSTER, MAXI_RANK, A, LA, POSELTT,
     &             NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, 
     &             TOL_OPT, KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
            ENDIF
            RANK_LIST_NEW(J) = LRB%K
            POS_LIST_NEW(J) = CURPOS
          ELSE
C           Only one child (the last of the current level), we leave it
C           alone            
            RANK_LIST_NEW(J) = NODE_RANK
            POS_LIST_NEW(J) = CURPOS
          ENDIF
          IOFF = IOFF+IMAX
        ENDDO
        IF (NB_NODES_NEW.GT.1) THEN
          IF (RESORT) THEN
C           Allocate new accumulator to sort            
            KTOT = SUM(RANK_LIST_NEW)
            CALL INIT_LRB(ACC_NEW,KTOT,M,N,.TRUE.)
            ALLOCATE(ACC_NEW%Q(MAXI_CLUSTER,MAXI_RANK),
     &           ACC_NEW%R(MAXI_RANK,MAXI_CLUSTER), stat=allocok)
            IF (allocok > 0) THEN
               write(*,*) 'Allocation error of ACC_NEW%Q/ACC_NEW%R ',
     &              'in DMUMPS_RECOMPRESS_ACC_NARYTREE'
               call MUMPS_ABORT()
            ENDIF
C           Sort each node by rank
            CALL MUMPS_SORT_INT(NB_NODES_NEW, RANK_LIST_NEW,
     &                                         POS_LIST_NEW)
C           Copy to new accumulator
            CURPOS = 1
            IF (LEVEL.EQ.0) THEN
              LRB_PTR => ACC_LRB
            ELSE
              LRB_PTR => ACC_TMP
            ENDIF
            DO J=1,NB_NODES_NEW
              DO L=0,RANK_LIST_NEW(J)-1
                ACC_NEW%Q(1:M,CURPOS+L) = 
     &               LRB_PTR%Q(1:M,POS_LIST_NEW(J)+L)
                ACC_NEW%R(CURPOS+L,1:N) = 
     &               LRB_PTR%R(POS_LIST_NEW(J)+L,1:N)
              ENDDO
              POS_LIST_NEW(J) = CURPOS
              CURPOS = CURPOS + RANK_LIST_NEW(J)
            ENDDO
            IF (LEVEL.GT.0) THEN
              CALL DEALLOC_LRB(ACC_TMP, KEEP8, 4)
            ENDIF
C           Recusive call          
            CALL DMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELTT, KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               RANK_LIST_NEW, POS_LIST_NEW, NB_NODES_NEW,
     &               LEVEL+1, ACC_NEW)
          ELSE
C           Recusive call          
            CALL DMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELTT, KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               RANK_LIST_NEW, POS_LIST_NEW, NB_NODES_NEW, LEVEL+1)
          ENDIF
        ELSE
          IF (POS_LIST_NEW(1).NE.1) THEN
            write(*,*) 'Internal error in ',
     &      'DMUMPS_RECOMPRESS_ACC_NARYTREE', POS_LIST_NEW(1) 
          ENDIF
          ACC_LRB%K = RANK_LIST_NEW(1)
          IF (RESORT.AND.LEVEL.GT.0) THEN
             DO L=1,ACC_LRB%K
                DO I=1,M
                   ACC_LRB%Q(I,L) = ACC_TMP%Q(I,L)
                ENDDO
                DO I=1,N
                   ACC_LRB%R(L,I) = ACC_TMP%R(L,I)
                ENDDO
            ENDDO
            CALL DEALLOC_LRB(ACC_TMP, KEEP8, 4)
          ENDIF
        ENDIF
        DEALLOCATE(RANK_LIST_NEW, POS_LIST_NEW)
      END SUBROUTINE DMUMPS_RECOMPRESS_ACC_NARYTREE
      SUBROUTINE DMUMPS_RECOMPRESS_ACC_V2(ACC_LRB, MAXI_CLUSTER,
     &          MAXI_RANK, A, LA, POSELTT, NFRONT, NIV,
     &          MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT_RMB, 
     &          KPERCENT_LUA, NEW_ACC_RANK)
C 
C  LUA related routine
C  This subroutine recompresses the LR accumulator ACC_LRB
C  In this version (_V2) we take advantage of the fact part of ACC_LRB%Q
C  is already orthogonal: it's of the form [Q1,Q2]*[R1;R2], where Q2 has
C  NEW_ACC_RANK colums.
C   1) Gram-Schmidt:     orthogonalize Q2 w.r.t. Q1       
C   2) Truncated RRQR:   recompress Q2: Q2 = Q2b*R2b       
C   3) Update ACC_LRB%R: ACC_LRB%R = [R1+(Q1^T*Q2)*R2;R2b*R2] 
        TYPE(LRB_TYPE),INTENT(INOUT) :: ACC_LRB
        INTEGER(8), intent(in)  :: LA
        DOUBLE PRECISION, intent(inout)  :: A(LA)
        INTEGER,INTENT(IN) :: NFRONT, NIV, TOL_OPT
C       MAXI_CLUSTER, MAXI_RANK are the dimensions of the ACC_LRB
C       allocated space.     
C       ACC_LRB%Q is of size MAXI_CLUSTER x MAXI_RANK        
C       ACC_LRB%R is of size MAXI_RANK x MAXI_CLUSTER        
C       The first ACC_LRB%K columns of ACC_LRB%Q/rows of ACC_LRB%R are
C       taken by the accumulator. 
        INTEGER,INTENT(IN) :: MAXI_CLUSTER, MAXI_RANK, KPERCENT_LUA
C       NEW_ACC_RANK is the size of the non-compressed part of the
C       accumulator. One needs to orthogonalize this part against
C       the rest of the accumulator. 
        INTEGER,INTENT(INOUT) :: NEW_ACC_RANK
        INTEGER(8), INTENT(IN) :: POSELTT
C       POSELTT = position of first entry in block to update with the
C       gemm
C       MIDBLK_COMPRESS = whether to compress the middle product
C       Y1*X2 in A = A + X1*(Y1*X2)*Y2 (strategy suggested by CA)
        INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT_RMB
C       Threshold used for RRQR if MIDBLK_COMPRESS >= 2 
        DOUBLE PRECISION, intent(in) :: TOLEPS
        DOUBLE PRECISION,    ALLOCATABLE :: RWORK_RRQR(:)
        DOUBLE PRECISION, ALLOCATABLE :: WORK_RRQR(:), TAU_RRQR(:)
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:), TARGET ::
     &                                Q1, R1, Q2, PROJ
        INTEGER, ALLOCATABLE :: JPVT_RRQR(:)
        INTEGER :: INFO, RANK1, MAXRANK, LWORK
        LOGICAL :: BUILDQ1
        LOGICAL, PARAMETER :: LRMIXED_ACTIVATED = .FALSE.
        INTEGER :: K_STOC(NPREC_STOC)
        INTEGER :: K453
        LOGICAL, PARAMETER :: MIXED_ADMISSIBILITY = .FALSE.
        INTEGER :: I, J, M, N, K, K1
        INTEGER :: allocok, MREQ
        DOUBLE PRECISION :: ONE, MONE, ZERO
        PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
        PARAMETER (ZERO=0.0D0)
              M = ACC_LRB%M
              N = ACC_LRB%N
              K = NEW_ACC_RANK
              K1 = ACC_LRB%K - K
              MAXRANK = K-1
              MAXRANK = max (1, int((MAXRANK*KPERCENT_LUA/100)))
              LWORK = K*(K+1)
              allocate(Q1(M,K), PROJ(K1, K),
     &             WORK_RRQR(LWORK), RWORK_RRQR(2*K), 
     &             TAU_RRQR(K),
     &             JPVT_RRQR(K), stat=allocok)
              IF (allocok > 0) THEN
                MREQ = M * K + K1 * K + LWORK + 4 * K
                GOTO 100
              ENDIF
              DO J=1,K
                  DO I=1,M
                      Q1(I,J) = ACC_LRB%Q(I,J+K1)
                  ENDDO
              ENDDO
C             1/ Orthogonalize Q1 (new columns) w.r.t the rest of the
C             accumulator              
              CALL dgemm('T', 'N', K1, K, M, ONE, ACC_LRB%Q(1,1), 
     &                  MAXI_CLUSTER, Q1(1,1), M, ZERO, PROJ(1,1), K1)
              CALL dgemm('N', 'N', M, K, K1, MONE, ACC_LRB%Q(1,1),
     &                  MAXI_CLUSTER, PROJ(1,1), K1, ONE, Q1(1,1), M)
C             2/ Recompress Q1 
              JPVT_RRQR = 0
              CALL DMUMPS_TRUNCATED_RRQR(M, K, Q1(1,1),
     &               M, JPVT_RRQR(1), TAU_RRQR(1), WORK_RRQR(1),
     &               K, RWORK_RRQR(1), TOLEPS, TOL_OPT, 
     &               RANK1, MAXRANK, INFO,
     &               K_STOC, LRMIXED_ACTIVATED, K453, 
     &               MIXED_ADMISSIBILITY,
     &               BUILDQ1)
              IF (BUILDQ1) THEN
C               Successfully compressed Q1
                allocate(Q2(N,K), stat=allocok)
                IF (allocok > 0) THEN
                   MREQ = N*K
                   GOTO 100
                ENDIF
                DO J=1,K
                    DO I=1,N
                        Q2(I,J) = ACC_LRB%R(J+K1,I)
                    ENDDO
                ENDDO
                CALL dgemm('N', 'T', K1, N, K, ONE, PROJ(1,1), K1,
     &                  Q2(1,1), N, ONE, ACC_LRB%R(1,1), MAXI_RANK)
                IF (RANK1.GT.0) THEN
                   allocate(R1(RANK1,K), stat=allocok)
                   IF (allocok > 0) THEN
                      MREQ = RANK1*K
                      GOTO 100
                   ENDIF
                  DO J=1, K
                     R1(1:MIN(RANK1,J),JPVT_RRQR(J)) =
     &                 Q1(1:MIN(RANK1,J),J)
                     IF(J.LT.RANK1) R1(MIN(RANK1,J)+1:
     &                 RANK1,JPVT_RRQR(J))= ZERO
                  END DO
                  CALL dorgqr 
     &                (M, RANK1, RANK1, Q1(1,1),
     &                M, TAU_RRQR(1),  
     &                WORK_RRQR(1), LWORK, INFO )
                  DO J=1,K
                    DO I=1,M
                        ACC_LRB%Q(I,J+K1) = Q1(I,J)
                    ENDDO
                  ENDDO
                  CALL dgemm('N', 'T', RANK1, N, K, ONE, R1(1,1), RANK1,
     &                   Q2(1,1), N, ZERO, ACC_LRB%R(K1+1,1), MAXI_RANK)
                  deallocate(R1)
                ENDIF
                deallocate(Q2)
                ACC_LRB%K = K1 + RANK1
              ENDIF
              deallocate(PROJ)
              deallocate(Q1, JPVT_RRQR, TAU_RRQR, WORK_RRQR, RWORK_RRQR)
              RETURN
 100          CONTINUE
C     Alloc NOT ok!!
              write(*,*) 'Allocation problem in BLR routine 
     &          DMUMPS_RECOMPRESS_ACC_V2: ',
     &          'not enough memory? memory requested = ' , MREQ        
              CALL MUMPS_ABORT()
              RETURN
      END SUBROUTINE DMUMPS_RECOMPRESS_ACC_V2
      SUBROUTINE MAX_CLUSTER(CUT,CUT_SIZE,MAXI_CLUSTER)
        INTEGER, intent(in) :: CUT_SIZE
        INTEGER, intent(out) :: MAXI_CLUSTER
        INTEGER, DIMENSION(:), intent(in) :: CUT
        INTEGER :: I
        MAXI_CLUSTER = 0
        DO I = 1, CUT_SIZE
          IF (CUT(I+1) - CUT(I) .GE. MAXI_CLUSTER) THEN
            MAXI_CLUSTER = CUT(I+1) - CUT(I)
          END IF
        END DO
      END SUBROUTINE MAX_CLUSTER
      SUBROUTINE DMUMPS_GET_LUA_ORDER(NB_BLOCKS, ORDER, RANK, IWHANDLER,
     &                       SYM, FS_OR_CB, I, J, FRFR_UPDATES, 
     &                       LBANDSLAVE_IN, K474, BLR_U_COL)
C       This subroutines computes the order in which the BLR updates of
C       a given block A(X,Y) will be accumulated.
C       If called from FS update routine (FS_OR_CB = 0), then
C         - if J = 0 (L panel), then X = I and Y = NB_BLOCKS+1
C         - if J = 1 (U panel), then X = NB_BLOCKS+1 and Y = I
C       If called from CB update routine (FS_OR_CB = 1), then 
C         X = I and Y = J.
C       Currently it orders the ranks in increasing order (exact or
C       estimate if KEEP(481) = 1).
C     -----------
C     Parameters
C     -----------
        INTEGER, INTENT(IN) :: NB_BLOCKS, IWHANDLER, SYM, FS_OR_CB, I, J
        INTEGER, INTENT(OUT) :: ORDER(NB_BLOCKS), RANK(NB_BLOCKS), 
     &                          FRFR_UPDATES
        LOGICAL, OPTIONAL, INTENT(IN) :: LBANDSLAVE_IN
        INTEGER, OPTIONAL, INTENT(IN) :: K474
        TYPE(LRB_TYPE), POINTER, OPTIONAL :: BLR_U_COL(:)
C     -----------
C     Local variables
C     -----------
        INTEGER :: K, IND_L, IND_U
        LOGICAL :: LBANDSLAVE
        TYPE(LRB_TYPE), POINTER :: BLR_L(:), BLR_U(:)
        IF (PRESENT(LBANDSLAVE_IN)) THEN
          LBANDSLAVE = LBANDSLAVE_IN
        ELSE
          LBANDSLAVE = .FALSE.
        ENDIF
        IF ((SYM.NE.0).AND.(FS_OR_CB.EQ.0).AND.(J.NE.0)) THEN
C         In LDLT, the FS panel should always be L
         write(6,*)
     &     "Internal error in DMUMPS_GET_LUA_ORDER",
     &     'SYM, FS_OR_CB, J = ',SYM,FS_OR_CB,J
          CALL MUMPS_ABORT()
        ENDIF
        FRFR_UPDATES = 0
        DO K = 1, NB_BLOCKS
C         Initialize ORDER array
          ORDER(K) = K
C         Compute IND_L and IND_U (i.e., X and Y shifted by NB_BLOCKS)                
          IF (FS_OR_CB.EQ.0) THEN ! FS
            IF (J.EQ.0) THEN ! L panel
              IND_L = NB_BLOCKS+I-K
              IND_U = NB_BLOCKS+1-K
            ELSE ! U panel
              IND_L = NB_BLOCKS+1-K
              IND_U = NB_BLOCKS+I-K
            ENDIF
          ELSE ! CB
            IND_L = I-K
            IND_U = J-K
          ENDIF
          IF (LBANDSLAVE) THEN
            IND_L = I
            IF (K474.GE.2) THEN
              IND_U = K
            ENDIF
          ENDIF
C         Retrieve L and U panels (make BLR_U point to BLR_L if LDLT) 
          CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &         IWHANDLER,
     &         0, ! L Panel
     &         K, BLR_L)
          IF (SYM.EQ.0) THEN
            IF (LBANDSLAVE.AND.K474.GE.2) THEN
              BLR_U => BLR_U_COL
            ELSE
              CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &         IWHANDLER,
     &         1, ! L Panel
     &         K, BLR_U)
            ENDIF
          ELSE
            BLR_U => BLR_L
          ENDIF
C         Extract rank of each product (exact or estimate if KEEP(481)=1)          
          IF (BLR_L(IND_L)%ISLR) THEN
            IF (BLR_U(IND_U)%ISLR) THEN
              RANK(K) = min(BLR_L(IND_L)%K, BLR_U(IND_U)%K)
            ELSE
              RANK(K) = BLR_L(IND_L)%K
            ENDIF
          ELSE
            IF (BLR_U(IND_U)%ISLR) THEN
              RANK(K) = BLR_U(IND_U)%K
            ELSE
C              RANK(K) = BLR_U(IND_U)%M
              RANK(K) = -1
              FRFR_UPDATES = FRFR_UPDATES + 1
            ENDIF
          ENDIF
        ENDDO
C       Sort the ranks in increasing order          
        CALL MUMPS_SORT_INT(NB_BLOCKS, RANK, ORDER)
      END SUBROUTINE DMUMPS_GET_LUA_ORDER
      SUBROUTINE DMUMPS_BLR_ASM_NIV1 (A, LA, POSEL1, NFRONT, NASS1, 
     &      IWHANDLER, SON_IW, LIW, LSTK, NELIM, K1, K2, SYM,
     &      KEEP, KEEP8, OPASSW)
C
C Purpose
C =======
C     
C       Called by a level 1 master assembling the contribution
C       block of a level 1 son that has been BLR-compressed        
C
C
C Parameters
C ==========
C
        INTEGER(8) :: LA, POSEL1
        INTEGER :: LIW, NFRONT, NASS1, LSTK, NELIM, K1, K2, IWHANDLER
        DOUBLE PRECISION :: A(LA)
C       INTEGER :: SON_IW(LIW)
        INTEGER :: SON_IW(:) ! contiguity information lost but no copy
        INTEGER :: KEEP(500)
        INTEGER(8) :: KEEP8(150)
        INTEGER :: SYM
        DOUBLE PRECISION, INTENT(INOUT) :: OPASSW
C
C Local variables
C ===============
C
        DOUBLE PRECISION, ALLOCATABLE :: SON_A(:)
        INTEGER(8) :: APOS, SON_APOS, IACHK, JJ2, NFRONT8
        INTEGER :: KK, KK1, allocok, SON_LA
        INTEGER :: IFLAG, IERROR
        TYPE(LRB_TYPE), POINTER :: CB_LRB(:,:), LRB
        INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_DYNAMIC
        INTEGER :: NB_INCB, NB_INASM, NB_BLR, I, J, M, N, II, NPIV,
     &        IBIS, IBIS_END, FIRST_ROW, LAST_ROW, FIRST_COL, LAST_COL,
     &        SON_LDA
        DOUBLE PRECISION :: PROMOTE_COST
        DOUBLE PRECISION :: ONE, ZERO
        PARAMETER (ONE = 1.0D0)
        PARAMETER (ZERO = 0.0D0)
        IFLAG  = 0
        IERROR = 0
        CALL DMUMPS_BLR_RETRIEVE_BEGSBLR_DYN(IWHANDLER,
     &                       BEGS_BLR_DYNAMIC)
        CALL DMUMPS_BLR_RETRIEVE_CB_LRB(IWHANDLER, CB_LRB)
        NB_BLR   = size(BEGS_BLR_DYNAMIC)-1
        NB_INCB  = size(CB_LRB,1)
        NB_INASM = NB_BLR - NB_INCB
        NPIV = BEGS_BLR_DYNAMIC(NB_INASM+1)-1
        NFRONT8 = int(NFRONT,8)
        IF (SYM.EQ.0) THEN
          IBIS_END = NB_INCB*NB_INCB
        ELSE
          IBIS_END = NB_INCB*(NB_INCB+1)/2
        ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL         
!$OMP DO PRIVATE(IBIS, I, J, M, N, SON_LA, SON_LDA, FIRST_ROW,
!$OMP&         LAST_ROW, FIRST_COL, LAST_COL, LRB, SON_A, II, KK,
!$OMP&         APOS, IACHK, KK1, JJ2, PROMOTE_COST, allocok, SON_APOS)
#endif
        DO IBIS = 1,IBIS_END
C         Determining I,J from IBIS
          IF (SYM.EQ.0) THEN
            I = (IBIS-1)/NB_INCB+1
            J = IBIS - (I-1)*NB_INCB
          ELSE
            I = ceiling((1.0D0+sqrt(1.0D0+8.0D0*dble(IBIS)))/2.0D0)-1
            J = IBIS - I*(I-1)/2
          ENDIF
          I = I+NB_INASM
          J = J+NB_INASM
C  TMJY:  -NPIV-NELIM
          IF (I.EQ.NB_INASM+1) THEN
C           first CB block, add NELIM because FIRST_ROW starts at NELIM+1
            FIRST_ROW = BEGS_BLR_DYNAMIC(I)-NPIV+NELIM
          ELSE
            FIRST_ROW = BEGS_BLR_DYNAMIC(I)-NPIV
          ENDIF
          LAST_ROW  = BEGS_BLR_DYNAMIC(I+1)-1-NPIV
          M=LAST_ROW-FIRST_ROW+1
          FIRST_COL = BEGS_BLR_DYNAMIC(J)-NPIV
          LAST_COL  = BEGS_BLR_DYNAMIC(J+1)-1-NPIV
          N = BEGS_BLR_DYNAMIC(J+1)-BEGS_BLR_DYNAMIC(J)
          SON_APOS  = 1_8
          SON_LA    = M*N
          SON_LDA   = N
          LRB => CB_LRB(I-NB_INASM,J-NB_INASM)
          IF (LRB%ISLR.AND.LRB%K.EQ.0) THEN
C           No need to perform extend-add              
            CALL DEALLOC_LRB(LRB, KEEP8, KEEP(34))
            NULLIFY(LRB)
            CYCLE
          ENDIF
          allocate(SON_A(SON_LA),stat=allocok)
          IF (allocok.GT.0) THEN
            write(*,*) 'Not enough memory in DMUMPS_BLR_ASM_NIV1',
     &     ", Memory requested = ", SON_LA
            CALL MUMPS_ABORT()
          ENDIF
C         decompress block        
          IF (LRB%ISLR) THEN
            IF (KEEP(452).GT.0) THEN
                CALL DMUMPS_DECONVERT_LRB(LRB, IFLAG, IERROR, KEEP8)
                IF (IFLAG.LT.0) THEN
                    write(*,*) 'Not enough memory in ',
     &              "DMUMPS_BLR_ASM_NIV1, Memory requested = ", IERROR
                    CALL MUMPS_ABORT()
                ENDIF
            ENDIF
            CALL dgemm('T', 'T', N, M, LRB%K, ONE, LRB%R(1,1), LRB%K,
     &                 LRB%Q(1,1), M, ZERO, SON_A(SON_APOS), SON_LDA)
            PROMOTE_COST = 2.0D0*M*N*LRB%K
            CALL UPD_FLOP_DECOMPRESS(PROMOTE_COST, .TRUE.)
          ELSE
            IF (I.EQ.J.AND.SYM.NE.0) THEN
C             Diag block and LDLT, copy only lower half              
             IF (J-NB_INASM.EQ.1.AND.NELIM.GT.0) THEN
C             The first diagonal block is rectangular !!
C             with NELIM more cols than rows
              DO II=1,M
                DO KK=1,II+NELIM
                  SON_A(SON_APOS+int(II-1,8)*int(SON_LDA,8) + 
     &                           int(KK-1,8))
     &                 = LRB%Q(II,KK)
                ENDDO
              ENDDO  
             ELSE
              DO II=1,M
                DO KK=1,II
                  SON_A(SON_APOS+int(II-1,8)*int(SON_LDA,8) + 
     &                           int(KK-1,8))
     &                 = LRB%Q(II,KK)
                ENDDO
              ENDDO  
             ENDIF
            ELSE
              DO II=1,M
                DO KK=1,N
                  SON_A(SON_APOS+int(II-1,8)*int(SON_LDA,8) + 
     &                           int(KK-1,8))
     &                 = LRB%Q(II,KK)
                ENDDO
              ENDDO  
            ENDIF
          ENDIF
C         Deallocate block
          CALL DEALLOC_LRB(LRB, KEEP8, KEEP(34))
          NULLIFY(LRB)
C         extend add in father                           
          IF (SYM.NE.0.AND.J-NB_INASM.EQ.1.AND.NELIM.GT.0) THEN
C           Case of LDLT with NELIM: first-block column is treated
C           differently as the NELIM are assembled at the end of the
C           father
              DO KK = FIRST_ROW, LAST_ROW
                IACHK = 1_8 + int(KK-FIRST_ROW,8)*int(SON_LDA,8)
                IF (SON_IW(KK+K1-1).LE.NASS1) THEN
C                 Fully summed row of the father => permute destination in
C                 father, symmetric swap to be done 
C                 First NELIM columns
                  APOS = POSEL1 + int(SON_IW(KK+K1-1),8) - 1_8
                  DO KK1 = FIRST_COL, FIRST_COL+NELIM-1
                   JJ2 = APOS + int(SON_IW(K1+KK1-1)-1,8)*NFRONT8
                   A(JJ2) = A(JJ2) + SON_A(IACHK + int(KK1-FIRST_COL,8))
                  ENDDO
C                 Remaining columns
                  APOS = POSEL1 + int(SON_IW(KK+K1-1)-1,8)*NFRONT8
C                 DO KK1 = FIRST_COL+NELIM, LAST_COL
C                 In case I=J and first block, one may have
C                 LAST_COL > KK, but only lower triangular part
C                 should be assembled. We use min(LAST_COL,KK)
C                 below index to cover this case.
                  DO KK1 = FIRST_COL+NELIM, min(LAST_COL,KK)
                   JJ2 = APOS + int(SON_IW(K1+KK1-1),8) - 1_8
                   A(JJ2) = A(JJ2) + SON_A(IACHK + int(KK1-FIRST_COL,8))
                  ENDDO
                ELSE
                  APOS = POSEL1 + int(SON_IW(KK+K1-1)-1,8)*NFRONT8
                  DO KK1 = FIRST_COL, min(LAST_COL,KK)
                   JJ2 = APOS + int(SON_IW(K1+KK1-1),8) - 1_8
                   A(JJ2) = A(JJ2) + SON_A(IACHK + int(KK1-FIRST_COL,8))
                  ENDDO
                ENDIF
              ENDDO
          ELSE
C           Case of LDLT without NELIM or LU: everything is simpler
            DO KK = FIRST_ROW, LAST_ROW
              APOS = POSEL1 + int(SON_IW(KK+K1-1)-1,8)*NFRONT8
              IACHK = 1_8 + int(KK-FIRST_ROW,8)*int(SON_LDA,8)
              IF (I.EQ.J.AND.SYM.NE.0) THEN
C               LDLT diag block: assemble only lower half
                DO KK1 = FIRST_COL, KK
                  JJ2 = APOS + int(SON_IW(K1+KK1-1),8) - 1_8
                  A(JJ2) = A(JJ2) + SON_A(IACHK + int(KK1-FIRST_COL,8))
                ENDDO
              ELSE
                DO KK1 = FIRST_COL, LAST_COL
                  JJ2 = APOS + int(SON_IW(K1+KK1-1),8) - 1_8
                  A(JJ2) = A(JJ2) + SON_A(IACHK + int(KK1-FIRST_COL,8))
                ENDDO
              ENDIF
            ENDDO
          ENDIF
C         Deallocate SON_A
          DEALLOCATE(SON_A)
        ENDDO
#if ! defined(BLR_NOOPENMP)
!$OMP END DO
!$OMP END PARALLEL         
#endif
        CALL DMUMPS_BLR_FREE_CB_LRB(IWHANDLER,
C          Only CB_LRB structure is left to deallocate
     &             .TRUE., KEEP8, KEEP(34))
        IF ((KEEP(486).EQ.3).OR.KEEP(486).EQ.0) THEN
C        Case of FR solve: the BLR structure could not be freed
C        in DMUMPS_END_FACTO_SLAVE and should be freed here
C        Not reachable in case of error: set INFO1 to 0
         CALL DMUMPS_BLR_END_FRONT(IWHANDLER, 0, KEEP8, KEEP(34),
     &                             MTK405=KEEP(405))
        ENDIF
      END SUBROUTINE DMUMPS_BLR_ASM_NIV1
C{
      SUBROUTINE DMUMPS_LRMIXED_CONVERT_PANEL(BLR_PANEL, NB_BLR, IP,
     &    K34, KEEP8, IFLAG, IERROR)
C
C Parameters
C ==========
C
        TYPE(LRB_TYPE), POINTER :: BLR_PANEL(:)
        INTEGER :: NB_BLR, IP, K34
        INTEGER(8) KEEP8(150)
        INTEGER :: IFLAG, IERROR
C
C Local variables
C ===============
C
        INTEGER :: I
#if ! defined(BLR_NOOPENMP)
!$OMP DO SCHEDULE(DYNAMIC)
#endif
        DO I = 1, NB_BLR-IP
          CALL DMUMPS_LRMIXED_CONVERT_LRB(BLR_PANEL(I), K34, KEEP8, 
     &                                    IFLAG, IERROR)
        ENDDO
#if ! defined(BLR_NOOPENMP)
!$OMP END DO NOWAIT
#endif
      END SUBROUTINE DMUMPS_LRMIXED_CONVERT_PANEL
      SUBROUTINE DMUMPS_LRMIXED_CONVERT_LRB(LRB, K34, KEEP8, IFLAG,
     &                                      IERROR)
C       this routine assumes LRB is stored in uniform PREC_FAC precision
C       this is always true as of writing this comment, but MAY NOT BE !
C       especially if mixed compression is introduced before schur
C       complement update to speed it up
C
C Parameters
C ==========
C
        TYPE(LRB_TYPE) :: LRB
        INTEGER  :: K34 
        INTEGER(8) KEEP8(150)
        INTEGER :: IFLAG, IERROR
C
C Local variables
C ===============
C       K34 is the size of an "INTEGER", in bytes (KEEP(34))
C       COSTS_PREC(i) is the number of bytes required for storing a 
C       coefficient in precision i
        INTEGER :: PREC, NB_INT_Q, NB_INT_R, N, M, K, Kprec, K_CUMUL
        INTEGER :: MEM, allocok
#if !defined(NOMIXHPC)
        INTEGER :: PREC_FRB
#endif
        MEM = 0
        M = LRB%M
        N = LRB%N
        IF (LRB%ISLR) THEN
          if (lrb%K.le.0) return
          allocate(LRB%Q_STOC(NPREC_STOC),
     &             LRB%R_STOC(NPREC_STOC), stat=allocok)
          IF (allocok.GT.0) THEN
C           IFLAG/IERROR are shared (it would make
C           sense to have them private and copy them
C           to the shared versions of IFLAG/IERROR
C           before the end of the // region.
!$OMP ATOMIC WRITE
            IFLAG = -13
!$OMP ATOMIC WRITE
            IERROR = 2*NPREC_STOC
            GOTO 100
          ENDIF
          DO PREC=1,NPREC_STOC
            NULLIFY(LRB%Q_STOC(PREC)%VAL)
            NULLIFY(LRB%R_STOC(PREC)%VAL)
          ENDDO
          K = LRB%K
          K_CUMUL = 1
          DO PREC=1,NPREC_STOC
            Kprec = LRB%K_STOC(PREC)
C           Skip empty blocks
            IF (KPREC.EQ.0) THEN
              nullify(LRB%Q_STOC(PREC)%val)
              nullify(LRB%R_STOC(PREC)%val)
            ELSE
C             keep(34) is the size of an "INTEGER", in bytes
C             allocate n coefs dans %val: allocate(Q_STOC(i)%val((n*cost+K34-1)/K34))
              NB_INT_Q = (M*Kprec*COSTS_PREC(PREC) + K34-1) / K34
C             ou bien, équivalent à:
C             NB_INT=CEILING( FLOAT(M*K*COSTS_PREC(PREC))/FLOAT(K34) )
              NB_INT_R = (N*Kprec*COSTS_PREC(PREC) + K34-1) / K34
              allocate(LRB%Q_STOC(PREC)%val(NB_INT_Q),
     &                 LRB%R_STOC(PREC)%val(NB_INT_R),
     &                 stat=allocok)
              IF (allocok .GT. 0) THEN
!$OMP ATOMIC WRITE
                IFLAG = -13
!$OMP ATOMIC WRITE
                IERROR = NB_INT_Q+NB_INT_R
                GOTO 100
              ENDIF
              CALL DMUMPS_LRMIXED_CAST(LRB%Q(1,K_CUMUL),
     &          LRB%Q_STOC(PREC)%val(1), M, Kprec, M, M, 
     &          PREC_FAC, PREC)
              CALL DMUMPS_LRMIXED_CAST(LRB%R(K_CUMUL,1),
     &          LRB%R_STOC(PREC)%val(1), Kprec, N, K, Kprec, 
     &          PREC_FAC, PREC)
C             Update memory consumption
              MEM = MEM 
     &            + 1 + (NB_INT_R*K34 -1) / COSTS_PREC(PREC_FAC)
     &            + 1 + (NB_INT_Q*K34 -1) / COSTS_PREC(PREC_FAC)
            ENDIF
            K_CUMUL = K_CUMUL + Kprec
          ENDDO
#if !defined(NOMIXHPC)
          LRB%KW_STOC = DOT_PRODUCT(
     &        LRB%K_STOC(2:NPREC_STOC),
     &        COSTS_PREC(2:NPREC_STOC))
#endif
C         Deallocate UNIFORM part and keep MIXED 
          CALL DEALLOC_LRB(LRB, KEEP8, K34, .FALSE., .TRUE.)
#if !defined(NOMIXHPC)
       ELSE 
         ! LRMIXED_ACTIVATED could be passed but testing association 
         ! simply works here
         CALL FRB_GET_PREC(LRB, PREC_FRB, associated(LRB%K_STOC))
         IF (PREC_FRB.eq.PREC_FAC.and.PREC_FAC.eq.1) RETURN
         allocate(LRB%Q_STOC(NPREC_STOC), stat=allocok)
         IF (allocok.GT.0) THEN
!$OMP ATOMIC WRITE
           IFLAG = -13
!$OMP ATOMIC WRITE
           IERROR = NPREC_STOC
           GOTO 100
         ENDIF
         DO PREC=1,NPREC_STOC
           NULLIFY(LRB%Q_STOC(PREC)%VAL)
         ENDDO
         PREC = PREC_FRB
         NB_INT_Q = (M*N*COSTS_PREC(PREC) + K34-1) / K34
         allocate(LRB%Q_STOC(PREC)%val(NB_INT_Q), stat=allocok)
         IF (allocok.GT.0) THEN
!$OMP ATOMIC WRITE
           IFLAG = -13
!$OMP ATOMIC WRITE
           IERROR = -NB_INT_Q
           GOTO 100
         ENDIF
         MEM = 1 + ( NB_INT_Q*K34 - 1 )/COSTS_PREC(PREC_FAC)
         CALL DMUMPS_LRMIXED_CAST(LRB%Q(1,1),
     &     LRB%Q_stoc(PREC)%VAL(1),       M, N, M, M, PREC_FAC, PREC)
#endif
C        Deallocate UNIFORM part and keep MIXED 
         CALL DEALLOC_LRB(LRB, KEEP8, K34, .FALSE., .TRUE.)
        ENDIF
        CALL MUMPS_DM_FAC_UPD_DYN_MEMCNTS(int(MEM,8),
     &        .TRUE., KEEP8, IFLAG, IERROR, .TRUE., .TRUE.)
C       Normal return
        RETURN
 100  CONTINUE
C     Return with error: not enough memory to perform
C     conversion from UNIFORM to MIXED precision.
C     Deallocate what was allocated in MIXED and keep UNIFORM
      CALL DEALLOC_LRB(LRB, KEEP8, K34, .TRUE., .FALSE.)
      RETURN
      END SUBROUTINE DMUMPS_LRMIXED_CONVERT_LRB
      SUBROUTINE DMUMPS_DECONVERT_LRB(LRB, IFLAG, IERROR, KEEP8)
C Parameters
C ==========
        TYPE(LRB_TYPE) :: LRB
        INTEGER(8) KEEP8(150)
        INTEGER, INTENT(INOUT) :: IFLAG, IERROR
C Local variables
C ===============
        INTEGER :: PREC, Kprec, K_CUMUL
        TYPE(LRB_TYPE) :: LRB2
        IF (LRB%ISLR .AND. LRB%K.GT.0) THEN
          CALL ALLOC_LRB(LRB2,LRB%K,LRB%M,LRB%N,.true.,IFLAG,
     &     IERROR,KEEP8)
          LRB%Q=>LRB2%Q
          LRB%R=>LRB2%R
          K_CUMUL=1
          DO PREC=1,NPREC_STOC
            Kprec = LRB%K_STOC(PREC)
            IF (KPREC.NE.0) THEN
              CALL DMUMPS_LRMIXED_CAST(LRB%Q_STOC(PREC)%val(1),
     &            LRB%Q(1,K_CUMUL), LRB%M, Kprec, LRB%M, LRB%M, 
     &               PREC, PREC_FAC)
              CALL DMUMPS_LRMIXED_CAST(LRB%R_STOC(PREC)%val(1),
     &            LRB%R(K_CUMUL,1), Kprec, LRB%N, Kprec, LRB%K, 
     &               PREC, PREC_FAC)
            ENDIF
            K_CUMUL = K_CUMUL + Kprec
          END DO
        ELSE
          IF (associated(LRB%Q_STOC)) THEN
            CALL ALLOC_LRB(LRB2,LRB%K,LRB%M,LRB%N,.true.,IFLAG,
     &       IERROR,KEEP8)
            LRB%Q=>LRB2%Q
            DO PREC=1,NPREC_STOC
              IF(ASSOCIATED(LRB%Q_STOC(prec)%VAL)) THEN
                CALL DMUMPS_LRMIXED_CAST(LRB%Q_STOC(PREC)%val(1),
     &              LRB%Q(1,1),LRB%M, LRB%N, LRB%M, LRB%M,
     &                 PREC, PREC_FAC)
                EXIT
              END IF
            END DO
          ENDIF
        END IF
      END SUBROUTINE DMUMPS_DECONVERT_LRB
#if !defined(NOMIXHPC)
C{
      SUBROUTINE DMUMPS_DECONVERT_LRB_R_W(LRB, R_WORKSPACE, LD, N, 
     &         START_PREC, END_PREC,
     &         START_PREC_K, END_PREC_K,
     &         IFLAG, IERROR, KEEP8)
        use iso_c_binding
C Parameters
C ==========
        TYPE(LRB_TYPE), intent(in) :: LRB
        INTEGER(8) KEEP8(150)
        INTEGER, INTENT(IN) :: LD, N
        INTEGER, INTENT(IN) :: START_PREC, END_PREC
        INTEGER, INTENT(IN) :: START_PREC_K, END_PREC_K
        INTEGER, INTENT(INOUT) :: IFLAG, IERROR
        DOUBLE PRECISION, INTENT(INOUT) :: R_WORKSPACE(:,:) !LD,N
C Local variables
C ===============
        INTEGER :: PREC, Kprec, K_CUMUL, Kbeg
        INTEGER :: Krem
        integer(kind=1), pointer, dimension(:) :: stoc_ptr
        integer :: ptr_shape(1)
        IF (LRB%ISLR .AND. LRB%K.GT.0) THEN
          K_CUMUL=1
          DO PREC=START_PREC,END_PREC
            Kprec = LRB%K_STOC(PREC)
            IF (KPREC.eq.0) cycle
            Kbeg = 1
            Krem = Kprec
            if (prec.eq.start_prec) then
              Kbeg = start_prec_k
              Krem = krem - start_prec_k + 1
            end if
            if (prec.eq.end_prec) then
              Krem = krem - (kprec - end_prec_k)
            end if
            ! FIXME GCC flags stop_ptr as creating temporary arrays
            ptr_shape(1) = costs_prec(prec)*Kprec*LRB%N
            call c_f_pointer(c_loc(LRB%R_stoc(prec)%val),
     &         stoc_ptr, ptr_shape)
            CALL DMUMPS_LRMIXED_CAST(
     &          stoc_ptr((kbeg-1)*costs_prec(prec)+1),
     &          R_WORKSPACE(K_CUMUL,1), Krem, LRB%N, Kprec, LD, 
     &             PREC, 1)
            K_CUMUL = K_CUMUL + Krem
          END DO
        END IF
      END SUBROUTINE DMUMPS_DECONVERT_LRB_R_W
      SUBROUTINE DMUMPS_DECONVERT_LRB_Q_W(LRB, Q_WORKSPACE, LD, N,
     &         START_PREC, END_PREC,
     &         START_PREC_K, END_PREC_K,
     &         IFLAG, IERROR, KEEP8)
        use iso_c_binding
C Parameters
C ==========
        TYPE(LRB_TYPE), intent(in) :: LRB
        INTEGER(8) KEEP8(150)
        INTEGER, INTENT(IN) :: LD, N
        INTEGER, INTENT(IN) :: START_PREC, END_PREC
        ! FIXME unused
        INTEGER, INTENT(IN) :: START_PREC_K, END_PREC_K
        INTEGER, INTENT(INOUT) :: IFLAG, IERROR
        DOUBLE PRECISION, INTENT(INOUT) :: Q_WORKSPACE(:,:) !LD,N)
C Local variables
C ===============
        INTEGER :: PREC, Kprec, K_CUMUL
        integer :: krem, kbeg
        integer(kind=1), pointer, dimension(:) :: stoc_ptr
        integer :: ptr_shape(1)
        IF (LRB%ISLR .AND. LRB%K.GT.0) THEN
          K_CUMUL=1
          DO PREC=start_prec,end_prec
            Kprec = LRB%K_STOC(PREC)
            if (kprec.eq.0) cycle
            Kbeg = 1
            Krem = Kprec
            if (prec.eq.start_prec) then
              Kbeg = start_prec_k
              Krem = krem - start_prec_k + 1
            end if
            if (prec.eq.end_prec) then
              Krem = krem - (kprec - end_prec_k)
            end if
            ptr_shape(1) = costs_prec(prec)*Kprec*LRB%M
            call c_f_pointer(c_loc(LRB%q_stoc(prec)%val),
     &         stoc_ptr, ptr_shape)
            CALL DMUMPS_LRMIXED_CAST(
     &          stoc_ptr(LRB%M*(kbeg-1)*costs_prec(prec)+1),
     &          Q_WORKSPACE(1,K_CUMUL), LRB%M, Krem, LRB%M, LD, 
     &             PREC, 1)
            K_CUMUL = K_CUMUL + Krem
          END DO
        END IF
      END SUBROUTINE DMUMPS_DECONVERT_LRB_Q_W
C} MIXHPC
#endif
C} LRMIXED
      END MODULE DMUMPS_LR_CORE
C --------------------------------------------------------------------
      SUBROUTINE DMUMPS_LRMIXED_CAST(A, B, nrows, ncols, ldA, ldB,
     &                               precA, precB)
        USE DMUMPS_LR_TYPE
C       This routine casts a generic integer matrix A of precision precA 
C       into generic integer matrix B of precision precB   
        INTEGER :: nrows, ncols, ldA, ldB, precA, precB
        INTEGER :: A(*), B(*) 
        IF (precA.EQ.prec_dbl .AND. precB.EQ.prec_dbl) THEN
          CALL MUMPS_LRMIXED_CAST_DD(A(1), B(1), nrows, ncols,
     &              ldA, ldB)
        ELSEIF (precA.EQ.prec_dbl .AND. precB.EQ.prec_sgl) THEN
          CALL MUMPS_LRMIXED_CAST_DS(A(1), B(1), nrows, ncols,
     &              ldA, ldB)
        ELSEIF (precA.EQ.prec_sgl .AND. precB.EQ.prec_dbl) THEN
          CALL MUMPS_LRMIXED_CAST_SD(A(1), B(1), nrows, ncols,
     &              ldA, ldB)
        ELSEIF (precA.EQ.prec_sgl .AND. precB.EQ.prec_sgl) THEN
          CALL MUMPS_LRMIXED_CAST_SS(A(1), B(1), nrows, ncols,
     &              ldA, ldB)
        ELSEIF (precA.EQ.prec_sgl .AND. precB.EQ.prec_16) THEN
          CALL MUMPS_LRMIXED_CAST_S16(A(1), B(1), nrows, ncols,
     &              ldA, ldB)
        ELSEIF (precA.EQ.prec_16 .AND. precB.EQ.prec_sgl) THEN
          CALL MUMPS_LRMIXED_CAST_16S(A(1), B(1), nrows, ncols,
     &              ldA, ldB)
        ELSEIF (precA.EQ.prec_dbl .AND. precB.EQ.prec_16) THEN
          CALL MUMPS_LRMIXED_CAST_D16(A(1), B(1), nrows, ncols,
     &              ldA, ldB)
        ELSEIF (precA.EQ.prec_16 .AND. precB.EQ.prec_dbl) THEN
          CALL MUMPS_LRMIXED_CAST_16D(A(1), B(1), nrows, ncols,
     &              ldA, ldB)
        ELSEIF (precA.EQ.prec_dbl .AND. precB.EQ.prec_24) THEN
          CALL MUMPS_LRMIXED_CAST_D24(A(1), B(1), nrows, ncols,
     &              ldA, ldB)
        ELSEIF (precA.EQ.prec_24 .AND. precB.EQ.prec_dbl) THEN
          CALL MUMPS_LRMIXED_CAST_24D(A(1), B(1), nrows, ncols,
     &              ldA, ldB)
        ELSEIF (precA.EQ.prec_sgl .AND. precB.EQ.prec_24) THEN
          CALL MUMPS_LRMIXED_CAST_S24(A(1), B(1), nrows, ncols,
     &              ldA, ldB)
        ELSEIF (precA.EQ.prec_24 .AND. precB.EQ.prec_sgl) THEN
          CALL MUMPS_LRMIXED_CAST_24S(A(1), B(1), nrows, ncols,
     &              ldA, ldB)
        ELSEIF (precA.EQ.prec_dbl .AND. precB.EQ.prec_40) THEN
          CALL MUMPS_LRMIXED_CAST_D40(A(1), B(1), nrows, ncols,
     &              ldA, ldB)
        ELSEIF (precA.EQ.prec_40 .AND. precB.EQ.prec_dbl) THEN
          CALL MUMPS_LRMIXED_CAST_40D(A(1), B(1), nrows, ncols,
     &              ldA, ldB)
        ELSEIF (precA.EQ.prec_dbl .AND. precB.EQ.prec_48) THEN
          CALL MUMPS_LRMIXED_CAST_D48(A(1), B(1), nrows, ncols,
     &              ldA, ldB)
        ELSEIF (precA.EQ.prec_48 .AND. precB.EQ.prec_dbl) THEN
          CALL MUMPS_LRMIXED_CAST_48D(A(1), B(1), nrows, ncols,
     &              ldA, ldB)
C Binary56:
        ELSEIF (precA.EQ.prec_dbl .AND. precB.EQ.prec_56) THEN
          CALL MUMPS_LRMIXED_CAST_D56(A(1), B(1), nrows, ncols,
     &              ldA, ldB)
        ELSEIF (precA.EQ.prec_56 .AND. precB.EQ.prec_dbl) THEN
          CALL MUMPS_LRMIXED_CAST_56D(A(1), B(1), nrows, ncols,
     &              ldA, ldB)
        ELSE 
          write(*,*) 'conversion impossible: format non pris en charge'
          write(*,*) 'precA', precA, 'precB', precB
          CALL MUMPS_ABORT()
        ENDIF
      END SUBROUTINE DMUMPS_LRMIXED_CAST
      SUBROUTINE DMUMPS_LRMIXED_PLUS(X, Y, nrows, ncols, ldX, ldY,  
     &       precX, precY)
        USE DMUMPS_LR_TYPE
C       This routine performs Y <- X+Y where X and Y are generic
C       integer matrices in precision precX and precY, respectively
        INTEGER :: nrows, ncols, ldX, ldY, precX, precY
        INTEGER :: X(*), Y(*)
        IF (precX .EQ. prec_dbl .AND. precY .EQ. prec_dbl) THEN
          CALL MUMPS_LRMIXED_PLUS_DD(X, Y, nrows, ncols, ldX, ldY)
        ELSEIF (precX .EQ. prec_dbl .AND. precY .EQ. prec_sgl) THEN
          CALL MUMPS_LRMIXED_PLUS_DS(X, Y, nrows, ncols, ldX, ldY)
        ELSEIF (precX .EQ. prec_sgl .AND. precY .EQ. prec_dbl) THEN
          CALL MUMPS_LRMIXED_PLUS_SD(X, Y, nrows, ncols, ldX, ldY)
        ELSEIF (precX .EQ. prec_sgl .AND. precY .EQ. prec_sgl) THEN
          CALL MUMPS_LRMIXED_PLUS_SS(X, Y, nrows, ncols, ldX, ldY)
        ELSE 
          write(*,*) 'conversion impossible: format non pris en charge'
          write(*,*) 'precX', precX, 'precY', precY
          CALL MUMPS_ABORT()
        ENDIF
      END SUBROUTINE DMUMPS_LRMIXED_PLUS
      SUBROUTINE DMUMPS_TRUNCATED_RRQR( M, N, A, LDA, JPVT, TAU, WORK,
     &     LDW, RWORK, TOLEPS, TOL_OPT, RANK, MAXRANK, INFO,
     &     K_STOC, LRMIXED_ACTIVATED, K453, MIXED_ADMISSIBILITY,
     &     ISLR)
C     This routine computes a Rank-Revealing QR factorization of a dense
C     matrix A. The factorization is truncated when the absolute value of
C     a diagonal coefficient of the R factor becomes smaller than a
C     prescribed threshold TOLEPS. The resulting partial Q and R factors
C     provide a rank-k approximation of the input matrix A with accuracy
C     TOLEPS.
C     
C     This routine is obtained by merging the LAPACK
C     (http://www.netlib.org/lapack/) CGEQP3 and CLAQPS routines and by
C     applying a minor modification to the outer factorization loop in
C     order to stop computations as soon as possible when the required
C     accuracy is reached.
C
C     Copyright (c) 1992-2017 The University of Tennessee and The 
C     University of Tennessee Research Foundation.  All rights reserved.
C     Copyright (c) 2000-2017 The University of California Berkeley. 
C     All rights reserved.
C     Copyright (c) 2006-2017 The University of Colorado Denver.  
C     All rights reserved.
C
C     Redistribution and use in source and binary forms, with or without
C     modification, are permitted provided that the following conditions
C     are met:
C
C      - Redistributions of source code must retain the above copyright
C        notice, this list of conditions and the following disclaimer.
C
C      - Redistributions in binary form must reproduce the above 
C        copyright notice, this list of conditions and the following 
C        disclaimer listed in this license in the documentation and/or 
C        other materials provided with the distribution.
C
C      - Neither the name of the copyright holders nor the names of its
C        contributors may be used to endorse or promote products derived from
C        this software without specific prior written permission.
C
C      The copyright holders provide no reassurances that the source code
C      provided does not infringe any patent, copyright, or any other
C      intellectual property rights of third parties.  The copyright holders
C      disclaim any liability to any recipient for claims brought against
C      recipient by any third party for infringement of that parties
C      intellectual property rights.
C
C      THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C      "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C      LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C      A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C      OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C      SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C      LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C      DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C      THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C      (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C      OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C
      USE DMUMPS_LR_TYPE
      IMPLICIT NONE
C
      INTEGER            ::  INFO, LDA, LDW, M, N, RANK, MAXRANK
C     TOL_OPT controls the tolerance option used      
C       >0 => use 2-norm (||.||_X = ||.||_2)
C       <0 => use Frobenius-norm (||.||_X = ||.||_F)
C     Furthermore, depending on abs(TOL_OPT):      
C       1 => absolute: ||B_{I(k+1:end),J(k+1:end)}||_X <= TOLEPS     
C       2 => relative to 2-norm of the compressed block: 
C        ||B_{I(k+1:end),J(k+1:end)}||_X <= TOLEPS*||B_{I,J}||_2
C       3 => relative to the max of the 2-norms of the row and column diagonal blocks 
C        ||B_{I(k+1:end),J{k+1:end}}||_X <= TOLEPS*max(||B_{I,I}||_2,||B_{J,J}||_2)
C       4 => relative to the sqrt of product of the 2-norms of the row and column diagonal blocks 
C        ||B_{I(k+1:end),J{k+1:end}}||_X <= TOLEPS*sqrt(||B_{I,I}||_2*||B_{J,J}||_2)
      INTEGER            ::  TOL_OPT
      DOUBLE PRECISION               ::  TOLEPS
      INTEGER            ::  JPVT(*)
      DOUBLE PRECISION               ::  RWORK(*)
      DOUBLE PRECISION            ::  A(LDA,*), TAU(*)
      DOUBLE PRECISION            ::  WORK(LDW,*)
      LOGICAL            ::  ISLR
      DOUBLE PRECISION               ::  TOLEPS_EFF, TRUNC_ERR
      INTEGER, PARAMETER ::  INB=1, INBMIN=2
      INTEGER            :: J, JB, MINMN, NB
      INTEGER            :: OFFSET, ITEMP
      INTEGER            :: LSTICC, PVT, K, RK
      DOUBLE PRECISION               :: TEMP, TEMP2, TOL3Z
      DOUBLE PRECISION            :: AKK
      LOGICAL INADMISSIBLE
      INTEGER           :: K_STOC(NPREC_STOC)
      LOGICAL           :: LRMIXED_ACTIVATED
      INTEGER           :: K453
      LOGICAL           :: MIXED_ADMISSIBILITY
      INTEGER PREC
      INTEGER CURRENT_PREC
#if !defined(NOMIXHPC)
      INTEGER MIN_CURRENT_PREC
#endif
      INTEGER COST_RK ! =sum(K_STOC*COSTS_PREC)
      DOUBLE PRECISION, PARAMETER    :: RZERO=0.0D+0, RONE=1.0D+0
      DOUBLE PRECISION :: ZERO
      DOUBLE PRECISION :: ONE
      PARAMETER          ( ONE = 1.0D+0 )
      PARAMETER          ( ZERO = 0.0D+0 ) 
      DOUBLE PRECISION               :: dlamch
      INTEGER            :: ilaenv, idamax
      EXTERNAL           :: idamax, dlamch
      EXTERNAL           dgeqrf, dormqr, xerbla
      EXTERNAL           ilaenv
      EXTERNAL           dgemm, dgemv, dlarfg, dswap
#if defined(nrm2_returns_double)
      DOUBLE PRECISION, EXTERNAL :: dnrm2
      DOUBLE PRECISION :: TMP_D
#else
      DOUBLE PRECISION, EXTERNAL :: dnrm2
#endif
      INFO = 0
      ISLR = .FALSE.
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( LDW.LT.N ) THEN
            INFO = -8
         END IF
      END IF
      IF( INFO.NE.0 ) THEN
         WRITE(*,999) -INFO
         RETURN
      END IF
      MINMN = MIN(M,N)
      IF( MINMN.EQ.0 ) THEN
         RANK = 0
         RETURN
      END IF
      NB = ilaenv( INB, 'CGEQRF', ' ', M, N, -1, -1 )
      IF (LRMIXED_ACTIVATED) THEN
        CURRENT_PREC=PREC_FAC !#1
        DO PREC=1,NPREC_STOC
          K_STOC(PREC)=0
        ENDDO
        COST_RK=0
      ENDIF
      SELECT CASE(abs(TOL_OPT))
      CASE(1)
        TOLEPS_EFF = TOLEPS
      CASE(2)
C      TOLEPS_EFF will be computed at step K=1 below        
      CASE DEFAULT
        write(*,*) 'Internal error in DMUMPS_TRUNCATED_RRQR: TOL_OPT =',
     &        TOL_OPT
        CALL MUMPS_ABORT()
      END SELECT
      TOLEPS_EFF = TOLEPS
C
C     Avoid pointers (and TARGET attribute on RWORK/WORK)
C     because of implicit interface. An implicit interface
C     is needed to avoid intermediate array copies
C     VN1  => RWORK(1:N)
C     VN2  => RWORK(N+1:2*N)
C     AUXV => WORK(1:LDW,1:1)
C     F    => WORK(1:LDW,2:NB+1)
C     LDF  =  LDW
*     Initialize partial column norms. The first N elements of work
*     store the exact column norms.
      DO J = 1, N
C        VN1( J ) = dnrm2( M, A( 1, J ), 1 )
#if defined(nrm2_returns_double)
         TMP_D      = dnrm2( M, A( 1, J ), 1 )
         RWORK( J ) = TMP_D
#else
         RWORK( J ) = dnrm2( M, A( 1, J ), 1 )
#endif
C        VN2( J ) = VN1( J )
         RWORK( N + J ) = RWORK( J )
         JPVT(J) = J
      END DO
      IF (TOL_OPT.LT.0) THEN
C       Compute TRUNC_ERR for first step              
C       TRUNC_ERR = dnrm2( N, VN1( 1 ), 1 )
#if defined(nrm2_returns_double)
        TMP_D     = dnrm2( N, RWORK( 1 ), 1 )
        TRUNC_ERR = TMP_D
#else
        TRUNC_ERR = dnrm2( N, RWORK( 1 ), 1 )
#endif
      ENDIF
      OFFSET = 0
      TOL3Z  = SQRT(dlamch('Epsilon'))
#if !defined(NOMIXHPC)
      MIN_CURRENT_PREC = NPREC_STOC
#endif
      DO 
         JB     = MIN(NB,MINMN-OFFSET)
         LSTICC = 0
         K      = 0
         DO 
            IF(K.EQ.JB) EXIT
            K   = K+1
            RK  = OFFSET+K
C           PVT = ( RK-1 ) + IDAMAX( N-RK+1, VN1( RK ), 1 )
            PVT = ( RK-1 ) + idamax( N-RK+1, RWORK( RK ), 1 )
            IF (RK.EQ.1) THEN 
C             IF (abs(TOL_OPT).EQ.2) TOLEPS_EFF = VN1(PVT)*TOLEPS
              IF (abs(TOL_OPT).EQ.2) TOLEPS_EFF = RWORK(PVT)*TOLEPS
            ENDIF
            IF (TOL_OPT.GT.0) THEN
C             TRUNC_ERR = VN1(PVT)
              TRUNC_ERR = RWORK(PVT)
C           ELSE
C             TRUNC_ERR has been already computed at previous step
            ENDIF
            IF(TRUNC_ERR.LT.TOLEPS_EFF) THEN
              RANK = RK-1
              ISLR = .TRUE.  
              RETURN
            ENDIF
            IF (LRMIXED_ACTIVATED) THEN
              DO WHILE(TRUNC_ERR.LT.
     &               TOLEPS_EFF/PRECISIONS(CURRENT_PREC+1)/K453)
                CURRENT_PREC = CURRENT_PREC+1
              ENDDO
#if !defined(NOMIXHPC)
              MIN_CURRENT_PREC = MIN(MIN_CURRENT_PREC, CURRENT_PREC)
#endif
              K_STOC(CURRENT_PREC) = K_STOC(CURRENT_PREC) + 1
              COST_RK = COST_RK + COSTS_PREC(CURRENT_PREC)
C             CURRENT_PREC is now the right precision format for the current column, 
C             because
C             TOLEPS_EFF/u(CURRENT_PREC+1) <= TRUNC_ERR < TOLEPS_EFF/u(CURRENT_PREC)
            ENDIF
            IF (.NOT.LRMIXED_ACTIVATED.OR..NOT.MIXED_ADMISSIBILITY) THEN
              INADMISSIBLE = (RK.GT.MAXRANK)
            ELSE
#if !defined(NOMIXHPC)
              INADMISSIBLE = 
     &           (COST_RK.GT.MAXRANK*COSTS_PREC(PREC_FAC))
              IF (QUANTIZE) THEN
                INADMISSIBLE = INADMISSIBLE.OR. 
     &             ((M+N)*COST_RK.GT.M*N*COSTS_PREC(MIN_CURRENT_PREC))
              ENDIF
#else
                INADMISSIBLE = 
     &             (COST_RK.GT.MAXRANK*COSTS_PREC(PREC_FAC))
#endif
            ENDIF
            IF (INADMISSIBLE) THEN
               RANK = RK
               INFO = RK
               ISLR = .FALSE.
               RETURN
            END IF
            IF( PVT.NE.RK ) THEN
               CALL dswap( M, A( 1, PVT ), 1, A( 1, RK ), 1 )
c              CALL dswap( K-1, F( PVT-OFFSET, 1 ), LDF,
c    &              F( K, 1 ), LDF )
               CALL dswap( K-1, WORK( PVT-OFFSET, 2 ), LDW,
     &              WORK( K, 2 ), LDW )
               ITEMP     = JPVT(PVT)
               JPVT(PVT) = JPVT(RK)
               JPVT(RK)  = ITEMP
C              VN1(PVT)  = VN1(RK)
C              VN2(PVT)  = VN2(RK)
               RWORK(PVT)    = RWORK(RK)
               RWORK(N+PVT)  = RWORK(N+RK)
            END IF
*     Apply previous Householder reflectors to column K:
*     A(RK:M,RK) := A(RK:M,RK) - A(RK:M,OFFSET+1:RK-1)*F(K,1:K-1)**H.
            IF( K.GT.1 ) THEN
               CALL dgemv( 'No transpose', M-RK+1, K-1, -ONE,
C    &              A(RK,OFFSET+1), LDA, F(K,1), LDF,
     &              A(RK,OFFSET+1), LDA, WORK(K,2), LDW,
     &              ONE, A(RK,RK), 1 )
            END IF
*     Generate elementary reflector H(k).
            IF( RK.LT.M ) THEN
               CALL dlarfg( M-RK+1, A(RK,RK), A(RK+1,RK), 1, TAU(RK) )
            ELSE
               CALL dlarfg( 1, A(RK,RK), A(RK,RK), 1, TAU(RK) )
            END IF
            AKK      = A(RK,RK)
            A(RK,RK) = ONE
*     Compute Kth column of F:
*     F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)**H*A(RK:M,K).
            IF( RK.LT.N ) THEN
               CALL dgemv( 'Transpose', M-RK+1, N-RK, TAU(RK),
     &              A(RK,RK+1), LDA, A(RK,RK), 1, ZERO,
C    &              F( K+1, K ), 1 )
     &              WORK( K+1, K+1 ), 1 )
            END IF
*     Padding F(1:K,K) with zeros.
            DO J = 1, K
C              F( J, K ) = ZERO
               WORK( J, K+1 ) = ZERO
            END DO
*     Incremental updating of F:
*     F(1:N,K) := F(1:N-OFFSET,K) - 
*             tau(RK)*F(1:N,1:K-1)*A(RK:M,OFFSET+1:RK-1)**H*A(RK:M,RK).
            IF( K.GT.1 ) THEN
               CALL dgemv( 'Transpose', M-RK+1, K-1, -TAU(RK),
     &              A(RK,OFFSET+1), LDA, A(RK,RK), 1, ZERO,
     &              WORK(1,1), 1 )
C    &              AUXV(1,1), 1 )
               CALL dgemv( 'No transpose', N-OFFSET, K-1, ONE,
     &              WORK(1,2), LDW, WORK(1,1), 1, ONE, WORK(1,K+1), 1 )
C    &              F(1,1), LDF, AUXV(1,1), 1, ONE, F(1,K), 1 )
            END IF
*     Update the current row of A:
*     A(RK,RK+1:N) := A(RK,RK+1:N) - A(RK,OFFSET+1:RK)*F(K+1:N,1:K)**H.
            IF( RK.LT.N ) THEN
C              CALL dgemv( 'No Transpose', N-RK, K, -ONE, F( K+1, 1 ), 
               CALL dgemv( 'No Transpose', N-RK, K, -ONE, WORK( K+1,2 ),
     &              LDW,
     &              A( RK, OFFSET+1 ), LDA, ONE, A( RK, RK+1 ), LDA )
            END IF
*     Update partial column norms.
*     
            IF( RK.LT.MINMN ) THEN
               DO J = RK + 1, N
C                 IF( VN1( J ).NE.RZERO ) THEN
                  IF( RWORK( J ).NE.RZERO ) THEN
*     
*     NOTE: The following 4 lines follow from the analysis in
*     Lapack Working Note 176.
*
C                    TEMP = ABS( A( RK, J ) ) / VN1( J )
                     TEMP = ABS( A( RK, J ) ) / RWORK( J )
                     TEMP = MAX( RZERO, ( RONE+TEMP )*( RONE-TEMP ) )
C                    TEMP2 = TEMP*( VN1( J ) / VN2( J ) )**2
                     TEMP2 = TEMP*( RWORK( J ) / RWORK( N+J ) )**2
                     IF( TEMP2 .LE. TOL3Z ) THEN
C                       VN2( J ) = dble( LSTICC )
                        RWORK( N+J ) = dble( LSTICC )
                        LSTICC = J
                     ELSE
C                       VN1( J ) = VN1( J )*SQRT( TEMP )
                        RWORK( J ) = RWORK( J )*SQRT( TEMP )
                     END IF
                  END IF
               END DO
            END IF
            A( RK, RK ) = AKK
            IF (LSTICC.NE.0) EXIT
            IF (TOL_OPT.LT.0) THEN
C             Compute TRUNC_ERR for next step              
C             TRUNC_ERR = dnrm2( N-RK, VN1( RK+1 ), 1 )
#if defined(nrm2_returns_double)
              TMP_D     = dnrm2( N-RK, RWORK( RK+1 ), 1 )
              TRUNC_ERR = TMP_D
#else
              TRUNC_ERR = dnrm2( N-RK, RWORK( RK+1 ), 1 )
#endif
            ENDIF
         END DO
*     Apply the block reflector to the rest of the matrix:
*     A(RK+1:M,RK+1:N) := A(RK+1:M,RK+1:N) -
*     A(RK+1:M,OFFSET+1:RK)*F(K+1:N-OFFSET,1:K)**H.
         IF( RK.LT.MIN(N,M) ) THEN
            CALL dgemm( 'No transpose', 'Transpose', M-RK,
     &           N-RK, K, -ONE, A(RK+1,OFFSET+1), LDA,
C    &           F(K+1,1), LDF, ONE, A(RK+1,RK+1), LDA )
     &           WORK(K+1,2), LDW, ONE, A(RK+1,RK+1), LDA )
         END IF
*     Recomputation of difficult columns.
         DO WHILE( LSTICC.GT.0 ) 
C           ITEMP = NINT( VN2( LSTICC ) )
            ITEMP = NINT( RWORK( N + LSTICC ) )
C           VN1( LSTICC ) = dnrm2( M-RK, A( RK+1, LSTICC ), 1 )
#if defined(nrm2_returns_double)
            TMP_D           = dnrm2( M-RK, A( RK+1, LSTICC ), 1 )
            RWORK( LSTICC ) = TMP_D
#else
            RWORK( LSTICC ) = dnrm2( M-RK, A( RK+1, LSTICC ), 1 )
#endif
*     
*     NOTE: The computation of RWORK( LSTICC ) relies on the fact that 
*     SNRM2 does not fail on vectors with norm below the value of
*     SQRT(DLAMCH('S')) 
*     
C           VN2( LSTICC ) = VN1( LSTICC )
            RWORK( N + LSTICC ) = RWORK( LSTICC )
            LSTICC = ITEMP
         END DO
         IF(RK.GE.MINMN) EXIT
         OFFSET = RK
         IF (TOL_OPT.LT.0) THEN
C          Compute TRUNC_ERR for next step              
C          TRUNC_ERR = dnrm2( N-RK, VN1( RK+1 ), 1 )
#if defined(nrm2_returns_double)
           TMP_D     = dnrm2( N-RK, RWORK( RK+1 ), 1 )
           TRUNC_ERR = TMP_D
#else
           TRUNC_ERR = dnrm2( N-RK, RWORK( RK+1 ), 1 )
#endif
         ENDIF
      END DO
      RANK = RK
      IF (.NOT.LRMIXED_ACTIVATED.OR..NOT.MIXED_ADMISSIBILITY) THEN
        ISLR = .NOT.(RK.GT.MAXRANK)
      ELSE
#if !defined(NOMIXHPC)
        ISLR = .NOT.(COST_RK.GT.MAXRANK*COSTS_PREC(CLOSE_PREC))
#else
        ISLR = .NOT.(COST_RK.GT.MAXRANK*COSTS_PREC(PREC_FAC))
#endif
      ENDIF
      RETURN
 999  FORMAT ('On entry to DMUMPS_TRUNCATED_RRQR, parameter number',
     &            I2,' had an illegal value')
      END SUBROUTINE DMUMPS_TRUNCATED_RRQR
       SUBROUTINE DMUMPS_LRMIXED_STOC2CALC_LRB(LRB, K34, K451, KEEP8,
     &    KEEP_stoc)
        USE DMUMPS_LR_TYPE
C
C Parameters
C ==========
C
        TYPE(LRB_TYPE) :: LRB
        INTEGER  :: K34 
        INTEGER(8) KEEP8(150)
        INTEGER :: K451
        LOGICAL :: KEEP_stoc
C Local variables
C ===============
        INTEGER :: N, M, K, Kprec_calc, Kprec_stoc
        INTEGER :: offset_Q_calc, NB_INT_Q_calc, offset_R_calc, 
     &  NB_INT_R_calc, new_ints_Qcalc, new_ints_Rcalc
        INTEGER :: Pcalc, Pstoc, Pcalc_stoc, PREC, cost_pcalc
        IF (LRB%ISLR .AND. LRB%K.GT.0) THEN
          IF(associated(LRB%Q_calc))THEN
            write(*,*) 'DMUMPS_LRMIXED_STOC2CALC_LRB: It seems that the 
     & conversion has already been done previously'
            CALL MUMPS_ABORT()
          ELSEIF (K451.EQ.0) THEN
C           KEEP 451 0 => precision formats for storage and 
C           computations are identical: no need to allocate anything
            ALLOCATE(LRB%Q_CALC(NPREC_CALC,NPREC_CALC))
            ALLOCATE(LRB%R_CALC(NPREC_CALC,NPREC_CALC))
            DO PREC=1,NPREC_CALC
              LRB%Q_CALC(PREC,PREC)%VAL=>LRB%Q_STOC(PREC)%VAL
              LRB%R_CALC(PREC,PREC)%VAL=>LRB%R_STOC(PREC)%VAL
            ENDDO
            LRB%K_CALC=>LRB%K_STOC
          ELSE  ! All cases of K451.NE.0 !
            IF (.NOT.ASSOCIATED(LRB%Q_STOC)) THEN
              write(*,*) 'Internal error in LRMIXED_STOC2CALC_LRB'
              CALL MUMPS_ABORT()
            ENDIF
C           KEEP 451 1 => precision formats for storage and computations
C           are different: effective conversions and copies needed
            ALLOCATE(LRB%Q_CALC(NPREC_CALC,NPREC_CALC))
            ALLOCATE(LRB%R_CALC(NPREC_CALC,NPREC_CALC))
            ALLOCATE(LRB%K_CALC(NPREC_CALC))
            M = LRB%M
            N = LRB%N
            K = LRB%K 
            DO PCALC=1,NPREC_CALC
              KPREC_CALC=SUM(
     &          LRB%K_STOC(CALC2STOC(PCALC):CALC2STOC(PCALC+1)-1))
              LRB%K_CALC(PCALC)=KPREC_CALC
              PCALC_STOC=CALC2STOC(PCALC)
              COST_PCALC=COSTS_PREC(PCALC_STOC)
              OFFSET_Q_CALC=0
              OFFSET_R_CALC=0
              NEW_INTS_QCALC=0
              NEW_INTS_RCALC=0
              IF(KPREC_CALC.GT.0) THEN
                NB_INT_Q_CALC = (M*KPREC_CALC
     &             *COST_PCALC + K34-1) / K34
                NB_INT_R_CALC = (N*KPREC_CALC
     &             *COST_PCALC + K34-1) / K34
                ALLOCATE(LRB%Q_CALC(PCALC,PCALC)%VAL(NB_INT_Q_CALC))
                ALLOCATE(LRB%R_CALC(PCALC,PCALC)%VAL(NB_INT_R_CALC))
                DO PSTOC=PCALC_STOC, CALC2STOC(PCALC+1)-1
                  KPREC_STOC = LRB%K_STOC(PSTOC)
                  IF(KPREC_STOC.LE.0) THEN
                    NEW_INTS_QCALC=0
                    NEW_INTS_RCALC=0
                  ELSE
                    NEW_INTS_QCALC=(M*KPREC_STOC*COSTS_PREC(PCALC_STOC) 
     &                 + K34-1) / K34
                    NEW_INTS_RCALC=KPREC_STOC*COSTS_PREC(PCALC_STOC)/K34
                    CALL DMUMPS_LRMIXED_CAST(LRB%Q_STOC(PSTOC)%VAL(1),
     &                LRB%Q_CALC(PCALC,PCALC)%VAL(1+OFFSET_Q_CALC), 
     &                M, KPREC_STOC, M, M, PSTOC, PCALC_STOC)
                    CALL DMUMPS_LRMIXED_CAST(LRB%R_STOC(PSTOC)%VAL(1),
     &                LRB%R_CALC(PCALC,PCALC)%VAL(1+OFFSET_R_CALC), 
     &                KPREC_STOC, N, KPREC_STOC, KPREC_CALC, 
     &                PSTOC, PCALC_STOC)
                  ENDIF
                  OFFSET_Q_CALC=OFFSET_Q_CALC+NEW_INTS_QCALC
                  OFFSET_R_CALC=OFFSET_R_CALC+NEW_INTS_RCALC
                ENDDO
            ENDIF
          ENDDO
          ENDIF
        ENDIF
      END SUBROUTINE DMUMPS_LRMIXED_STOC2CALC_LRB
#if defined(NOMIXHPC)
      SUBROUTINE DMUMPS_LRMIXED_LRxFR(LorU, LRB, B, LDB, C, LDC, NCOLS,
     &                   ALPHA, K34, K451, KEEP8, RHS, RB, QRB)
        USE DMUMPS_LR_TYPE
C       This routine performs C <- alpha*LRB*B + C, where LRB is a
C       mixed precision LR block
        TYPE(LRB_TYPE) :: LRB
        INTEGER ::  B(*), C(*)
        INTEGER ::  NCOLS, LDB, LDC
        DOUBLE PRECISION :: ALPHA
        INTEGER  :: K34, K451
        INTEGER(8) KEEP8(150)
C        INTEGER :: IFLAG, IERROR
C        LOGICAL :: KEEP_STOC
C       arrays of pointers:
C       RHS(p)%val is the copy of B in precision p
C       QRB(p)%val is the temporary array used to store LRB%Q_calc(p)*LRB%R_calc(p)*B
        type(LRMIXED_TYPE) :: RHS(*), QRB(*)
C       RB is the temporary array used to store LRB%R_calc(p)*B
        INTEGER :: RB(*)
C       LorU: 0 = L, 1 = U        
        INTEGER :: LorU
C       Local variables
        DOUBLE PRECISION :: GAMMA
        INTEGER :: M, N, K, KP, PREC, PCALC, PCALC_NEXT, I
        INTEGER :: PRECS_NONEMPTY(NPREC_CALC), NB_NONEMPTY
C       MorN = M if LorU=L, N if LorU=U (it is the size of QRB)
        INTEGER :: MorN
        M = LRB%M
        N = LRB%N
        K = LRB%K
        IF (M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0) THEN
          RETURN
        ENDIF
C       If needed, switch from STOC to CALC representation
        CALL DMUMPS_LRMIXED_STOC2CALC_LRB(LRB, K34, K451, KEEP8, .TRUE.)
C       We need to perform the operation 
C           C <- alpha*Q*(R*B) + C       if forward
C           C <- alpha*R^T*(Q^T*B) + C   if backward
C       Example: if forward and NPREC_CALC=3, we need to perform these operations:
C         QRB{3} <- alpha*Q{3} * (R{3}*B{3})
C         QRB{2} <- cast(QRB{3}, prec=2)
C         QRB{2} <- QRB{2} + alpha*Q{2} * (R{2}*B{2})
C         C <- C + cast(QRB{2}, prec=1)
C         C <- C + alpha*Q{1} * (R{1}*B{1})
C       Actually some of these steps may be skipped if precision format 1, 2 or 3 is empty.
C       We select all 'non-empty' precisions (except the working 
C       precision, which is a special case). We store them in the array 
C       PRECS_NONEMPTY of size NB_NONEMPTY:
        NB_NONEMPTY = 0
        DO I = 2, NPREC_CALC
          IF (LRB%K_CALC(I) .GT. 0) THEN
            NB_NONEMPTY = NB_NONEMPTY + 1
            PRECS_NONEMPTY(NB_NONEMPTY) = I
          ENDIF
        ENDDO
C       Loop on precisions >=2, in descending order. We can skip all 
C       empty blocks thanks to the array PRECS_NONEMPTY
C       Note that if NPREC_CALC=2, this loop has only 1 iteration, 
C       or even 0 iterations if LRB%Q_CALC(2) is empty
        DO I = NB_NONEMPTY, 1, -1
            Pcalc = PRECS_NONEMPTY(I)
            KP = LRB%K_CALC(Pcalc)
            PREC = CALC2STOC(Pcalc)
            IF(I.EQ.NB_NONEMPTY) THEN
C               For the first product QRB is empty -> gamma=0
                GAMMA=0.0D0
            ELSE 
                GAMMA=1.0D0
            ENDIF
            IF (LorU .EQ. 0) THEN
C             RB=R{p}*B{p}
C             QRB{P}=ALPHA*Q{P}*RB + GAMMA*QRB{P}:
              MorN = M
              CALL DMUMPS_LRMIXED_GEMM('N', 'N', KP, NCOLS, N, 1.0D0,
     &            LRB%R_CALC(PCALC,PCALC)%VAL(1), KP,
     &            RHS(PCALC)%VAL(1), N, 0.0D0,
     &            RB, KP, PREC)
              CALL DMUMPS_LRMIXED_GEMM('N', 'N', M, NCOLS, KP,
     &            ALPHA, LRB%Q_CALC(PCALC,PCALC)%VAL(1), M, 
     &            RB, KP, GAMMA,
     &            QRB(PCALC)%VAL(1), M, PREC)
            ELSE
C             RB=QT{p}*B{p}
C             QRB{P}=ALPHA*RT{P}*RB + GAMMA*QRB{P}:
              MorN = N
              CALL DMUMPS_LRMIXED_GEMM('T', 'N', KP, NCOLS, M, 1.0D0,
     &            LRB%Q_CALC(PCALC,PCALC)%VAL(1), M,
     &            RHS(PCALC)%VAL(1), M, 0.0D0,
     &            RB, KP, PREC)
              CALL DMUMPS_LRMIXED_GEMM('T', 'N', N, NCOLS, KP,
     &            ALPHA, LRB%R_CALC(PCALC,PCALC)%VAL(1), KP, 
     &            RB, KP, GAMMA,
     &            QRB(PCALC)%VAL(1), N, PREC)
            ENDIF
C           C=beta*C+cast(QRB{2}, prec=1) if this is the last precision (p=2, or more generally i=1),
C           or QRB{p-1}=cast(QRB{p}, prec=p-1) if p>2
            IF(I .EQ. 1) THEN
C             This is the last precision of the loop: we add the result to C
              CALL DMUMPS_LRMIXED_PLUS(QRB(PCALC)%VAL(1), 
     &            C, MorN, NCOLS, MorN, LDC, CALC2STOC(PCALC), 
     &            CALC2STOC(1))
            ELSE
C             Not the last precision: cast current QRB to next precision
              PCALC_NEXT = PRECS_NONEMPTY(I-1)
              CALL DMUMPS_LRMIXED_CAST(QRB(PCALC)%VAL(1), 
     &          QRB(PCALC_NEXT)%VAL(1), MorN, NCOLS, 
     &          MorN, MorN, CALC2STOC(PCALC), CALC2STOC(PCALC_NEXT))
            ENDIF
        ENDDO
C       PCALC=1, the working precision, is a special case and is handled
C       separately:  RB=R{1}*B{1}; C=alpha*Q{1}*RB+C
        PCALC = 1
        PREC  = CALC2STOC(PCALC)
        KP    = LRB%K_CALC(PCALC)
        IF(KP.GT.0) THEN
            IF (LorU .EQ. 0) THEN
              CALL DMUMPS_LRMIXED_GEMM('N', 'N', KP, NCOLS, N, 1.0D0,
     &          LRB%R_CALC(PCALC,PCALC)%VAL(1), KP,
     &          B, LDB, 0.0D0, RB(1), KP, PREC)
              CALL DMUMPS_LRMIXED_GEMM('N', 'N', M, NCOLS, KP,
     &          ALPHA, LRB%Q_CALC(PCALC,PCALC)%VAL(1), M, 
     &          RB(1), KP, 1.0D0, C, LDC, PREC)
            ELSE
              CALL DMUMPS_LRMIXED_GEMM('T', 'N', KP, NCOLS, M, 1.0D0,
     &          LRB%Q_CALC(PCALC,PCALC)%VAL(1), M,
     &          B, LDB, 0.0D0, RB(1), KP, PREC)
              CALL DMUMPS_LRMIXED_GEMM('T', 'N', N, NCOLS, KP,
     &          ALPHA, LRB%R_CALC(PCALC,PCALC)%VAL(1), KP, 
     &          RB(1), KP, 1.0D0, C, LDC, PREC)
            ENDIF
        ENDIF
C       In case of STOC representation, deallocate CALC now
C       TODO Other strategies are possible, for example keep CALC after
C       1st time and deallocate STOC        
C       AJ such a strategy is useful in a compute-bound context where
C       a low compte precision can be leveraged
        IF(LRB%ISLR .AND. K451.NE.0 ) THEN
          DO PCALC=1,NPREC_CALC
            IF(LRB%K_CALC(PCALC).GT.0) THEN
              DEALLOCATE(LRB%Q_CALC(PCALC,PCALC)%VAL)
              DEALLOCATE(LRB%R_CALC(PCALC,PCALC)%VAL)
              NULLIFY(LRB%Q_CALC(PCALC,PCALC)%VAL)
              NULLIFY(LRB%R_CALC(PCALC,PCALC)%VAL)
            ENDIF
          END DO
          DEALLOCATE(LRB%K_CALC)
          NULLIFY(LRB%K_CALC)
        ENDIF
        DEALLOCATE(LRB%Q_CALC)
        DEALLOCATE(LRB%R_CALC)
        NULLIFY(LRB%Q_CALC)
        NULLIFY(LRB%R_CALC)
        RETURN
      END SUBROUTINE DMUMPS_LRMIXED_LRxFR
C ^^ .not.MIXHPC
#else 
C vv MIXHPC
      SUBROUTINE DMUMPS_LRMIXED_LRxFR(LorU,  
     &                  LRB, B, LDB, C, LDC, NCOLS,
     &                  ALPHA, K34, K451, KEEP8, RHS, RB, QRB,
     &                  thread_workspace, tws_size) 
        use iso_c_binding
        USE DMUMPS_LR_TYPE
C       This routine performs C <- alpha*LRB*B + C, where LRB is a
C       mixed precision LR block
        implicit none
C       LorU: 0 = L, 1 = U        
        INTEGER, intent(in):: LorU
        TYPE(LRB_TYPE), intent(inout) :: LRB
        INTEGER, target, intent(inout) ::  B(*)
        INTEGER, target, intent(inout) :: C(*)
        INTEGER, intent(in) ::  NCOLS, LDB, LDC
        DOUBLE PRECISION, intent(in) :: ALPHA
        INTEGER, intent(in) :: K34, K451
        INTEGER(8), intent(in) :: KEEP8(150)
C        INTEGER :: IFLAG, IERROR
C        LOGICAL :: KEEP_STOC
C       arrays of pointers:
C       RHS(p)%val is the copy of B in precision p
C       QRB(p)%val is the temporary array used to store LRB%Q_calc(p)*LRB%R_calc(p)*B
        type(LRMIXED_TYPE), intent(inout) :: RHS(*), QRB(*)
C       RB is the temporary array used to store LRB%R_calc(p)*B
        INTEGER, intent(inout) :: RB(*)
C       the caller should allocate this array which is supposed enough to store
C       any Q_stoc/R_stoc from its precision to the next working precision
        integer, intent(in) :: tws_size
        DOUBLE PRECISION, target, intent(inout) :: 
     &          thread_workspace(tws_size)
!        type(c_ptr), intent(inout) :: 
!     &    thread_ptr
C       Local variables
        DOUBLE PRECISION :: GAMMA
        INTEGER :: M, N, K, PREC, PCALC
        integer :: pcalc_stoc, kprec_stoc, pstoc
C       MorN = M if LorU=L, N if LorU=U (it is the size of QRB)
        INTEGER :: MorN
        INTEGER(kind=1), pointer :: R_ptr(:,:)
        INTEGER(kind=1), pointer :: Q_ptr(:,:)
        INTEGER(kind=1), pointer :: Bptr(:), Cptr(:)
        INTEGER :: LD_B, LD_C
        logical :: gamma_touched
        integer :: last_prec_stoc, last_prec_calc
        integer :: ptr1_shape(1), ptr2_shape(2)
        M = LRB%M
        N = LRB%N
        K = LRB%K
        IF (M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0) THEN
          RETURN
        ENDIF
!        call c_f_pointer(thread_ptr, thread_workspace, (/ tws_size /))
C       If needed, switch from STOC to CALC representation
!        CALL DMUMPS_LRMIXED_STOC2CALC_LRB(LRB, K34, K451, KEEP8, .TRUE.)
C       We need to perform the operation 
C           C <- alpha*Q*(R*B) + C       if forward
C           C <- alpha*R^T*(Q^T*B) + C   if backward
C       Example: if forward and NPREC_CALC=3, we need to perform these operations:
C         QRB{3} <- alpha*Q{3} * (R{3}*B{3})
C         QRB{2} <- cast(QRB{3}, prec=2)
C         QRB{2} <- QRB{2} + alpha*Q{2} * (R{2}*B{2})
C         C <- C + cast(QRB{2}, prec=1)
C         C <- C + alpha*Q{1} * (R{1}*B{1})
C       Actually some of these steps may be skipped if precision format 1, 2 or 3 is empty.
        gamma_touched = .false. 
        IF (LorU .EQ. 0) THEN
          MorN = M
        else
          MorN = N
        end if
        last_prec_stoc = -1
        last_prec_calc = -1
        do PSTOC=NPREC_STOC,1,-1
          KPREC_STOC = LRB%K_STOC(PSTOC)
          if (kprec_stoc.eq.0) cycle
          pcalc = -1
          pcalc_stoc = -1
          do prec=NPREC_CALC,1,-1
            if (calc2stoc(prec).le.pstoc) then
              pcalc = prec
              exit 
            end if
          end do
          pcalc_stoc = calc2stoc(pcalc)
          if (pcalc_stoc.ne.1.and..not.gamma_touched) then
            gamma_touched = .true.
            gamma = 0.0D0
          else 
            GAMMA = 1.0D0
          end if
          if (pcalc_stoc.eq.1) then
             LD_B = LDB
             ptr1_shape(1) = LD_B*NCOLS*COSTS_PREC(1)
             call c_f_pointer(c_loc(B),Bptr,ptr1_shape)
             LD_C = LDC
             ptr1_shape(1) = LD_C*NCOLS*COSTS_PREC(1)
             call c_f_pointer(c_loc(C),Cptr,ptr1_shape)
             ! the last prec. used is less accurate than the working one
             if (last_prec_stoc.gt.pcalc_stoc) then
               CALL DMUMPS_LRMIXED_PLUS(QRB(last_prec_calc)%VAL(1), 
     &             Cptr(1), MorN, NCOLS, MorN, LD_C, last_prec_stoc, 1)
             end if
          else
             if (LorU.eq.0) then
              LD_B = N
             else
              LD_B = M
             endif
             ptr1_shape(1) = LD_B*NCOLS*COSTS_PREC(PCALC_STOC)
             call c_f_pointer(c_loc(RHS(pcalc)%val),Bptr,ptr1_shape)
             LD_C = MorN
             ptr1_shape(1) = LD_C*NCOLS*COSTS_PREC(PCALC_STOC)
             call c_f_pointer(c_loc(QRB(pcalc)%val),Cptr,ptr1_shape)
             ! the last precision used is less accurate than the one
             ! we are about to use
             if (last_prec_stoc.gt.pcalc_stoc) then
C              cast previousQRB to current precision
               CALL DMUMPS_LRMIXED_CAST(QRB(last_prec_calc)%VAL(1), 
     &           Cptr(1), MorN, NCOLS, 
     &           LD_C, LD_C, last_prec_stoc, pcalc_stoc)
             end if
          endif
          !TODO factorize
          IF (LorU .EQ. 0) THEN
            ptr2_shape(1) = Kprec_stoc*COSTS_PREC(PCALC_STOC)
            ptr2_shape(2) = N
            if (K451.eq.0.or.pstoc.eq.pcalc_stoc) then
              call c_f_pointer(c_loc(lrb%R_stoc(PSTOC)%val),
     &                         R_ptr, ptr2_shape)
            else
              call c_f_pointer(c_loc(thread_workspace),
     &                         R_ptr, ptr2_shape)
              CALL DMUMPS_LRMIXED_CAST(LRB%R_STOC(PSTOC)%VAL(1),
     &       R_ptr(1,1), !LRB%R_CALC(PCALC,PCALC)%VAL(1+OFFSET_R_CALC), 
     &       KPREC_STOC, N, KPREC_STOC, KPREC_STOC, 
     &          PSTOC, PCALC_STOC)
            end if
C           RB=R{p}*B{p}
            CALL DMUMPS_LRMIXED_GEMM('N', 'N', KPREC_STOC, NCOLS, N, 
     &          1.0D0, R_ptr(1,1), KPREC_STOC,
     &          Bptr(1), LD_B, 0.0D0,
     &          RB, KPREC_STOC, PCALC_STOC)
            ptr2_shape(1) = M*COSTS_PREC(PCALC_STOC)
            ptr2_shape(2) = Kprec_stoc
            if (K451.eq.0.or.pstoc.eq.pcalc_stoc) then
              call c_f_pointer(c_loc(lrb%Q_stoc(PSTOC)%val),
     &                          Q_ptr, ptr2_shape)
            else
              call c_f_pointer(c_loc(thread_workspace),
     &                         Q_ptr, ptr2_shape)
              CALL DMUMPS_LRMIXED_CAST(LRB%Q_STOC(PSTOC)%VAL(1),
     &       q_ptr(1,1), !LRB%Q_CALC(PCALC,PCALC)%VAL(1+OFFSET_Q_CALC), 
     &       M, KPREC_STOC, M, M, PSTOC, PCALC_STOC)
            end if
C           QRB{P}=ALPHA*Q{P}*RB + GAMMA*QRB{P}:
            CALL DMUMPS_LRMIXED_GEMM('N', 'N', M, NCOLS, KPREC_STOC,
     &          ALPHA, q_ptr(1,1), M,
     &          RB, KPREC_STOC, GAMMA,
     &          Cptr(1), LD_C, PCALC_STOC)
          ELSE
            ptr2_shape(1) = M*COSTS_PREC(PCALC_STOC)
            ptr2_shape(2) = Kprec_stoc
            if (K451.eq.0.or.pstoc.eq.pcalc_stoc) then
              call c_f_pointer(c_loc(lrb%Q_stoc(PSTOC)%val),
     &                          Q_ptr, ptr2_shape)
            else
              call c_f_pointer(c_loc(thread_workspace),
     &                          Q_ptr, ptr2_shape)
              CALL DMUMPS_LRMIXED_CAST(LRB%Q_STOC(PSTOC)%VAL(1),
     &      q_ptr(1,1), !LRB%Q_CALC(PCALC,PCALC)%VAL(1+OFFSET_Q_CALC), 
     &      M, KPREC_STOC, M, M, PSTOC, PCALC_STOC)
            end if
C           RB=QT{p}*B{p}
            CALL DMUMPS_LRMIXED_GEMM('T', 'N', KPREC_STOC, NCOLS, M, 
     &          1.0D0, Q_ptr(1,1), M, Bptr(1), LD_B, 
     &          0.0D0, RB, KPREC_STOC, PCALC_STOC)
            ptr2_shape(1) = Kprec_stoc*COSTS_PREC(PCALC_STOC)
            ptr2_shape(2) = N
            if (K451.eq.0.or.pstoc.eq.pcalc_stoc) then
              call c_f_pointer(c_loc(lrb%R_stoc(PSTOC)%val(1)),
     &                         R_ptr, ptr2_shape)
            else
              call c_f_pointer(c_loc(thread_workspace),
     &                         R_ptr, ptr2_shape)
              CALL DMUMPS_LRMIXED_CAST(LRB%R_STOC(PSTOC)%VAL(1),
     &       R_ptr(1,1), !LRB%R_CALC(PCALC,PCALC)%VAL(1+OFFSET_R_CALC), 
     &       KPREC_STOC, N, KPREC_STOC, KPREC_STOC, 
     &       PSTOC, PCALC_STOC)
            end if
C           QRB{P}=ALPHA*RT{P}*RB + GAMMA*QRB{P}:
            CALL DMUMPS_LRMIXED_GEMM('T', 'N', N, NCOLS, KPREC_STOC,
     &          ALPHA, R_ptr(1,1), KPREC_STOC, RB, KPREC_STOC, 
     &          GAMMA, Cptr(1), LD_C, PCALC_STOC)
          ENDIF
          last_prec_calc = pcalc
          last_prec_stoc = pcalc_stoc
        end do
        if (last_prec_stoc.gt.1) then
           ptr1_shape(1) = LDC*NCOLS*COSTS_PREC(1)
           call c_f_pointer(c_loc(C),Cptr,ptr1_shape)
           CALL DMUMPS_LRMIXED_PLUS(QRB(last_prec_calc)%VAL(1), 
     &         Cptr(1), MorN, NCOLS, MorN, LDC, last_prec_stoc,
     &         CALC2STOC(1))
         end if
      END SUBROUTINE DMUMPS_LRMIXED_LRxFR
      SUBROUTINE DMUMPS_LRMIXED_FRxFR(TRANSA, TRANSB,  
     &                  M, N, K, ALPHA, 
     &                  FRB, iA, jA, LDA, 
     &                  B, LDB, C, LDC,
     &                  RHS, QB, thread_workspace, tws_size)
        use iso_c_binding
        USE DMUMPS_LR_TYPE ! for lrb_type, lrmixed_type
        USE DMUMPS_LR_CORE ! for frb_get_prec
C       This routine performs 
C       C <- alpha*FRB(iA:,jA:)*B(ja:,1:) + C, where FRB is a
C       uniform some-precision FR block
C       B and C are assumed to be working precision blocks
        implicit none
C       LorU: 0 = L, 1 = U        
        CHARACTER, intent(in):: TRANSA, TRANSB
        INTEGER, intent(in) :: M, N, K
        DOUBLE PRECISION, intent(in) :: ALPHA
        TYPE(LRB_TYPE), intent(inout) :: FRB
        INTEGER, intent(in) :: iA, jA, LDA
        INTEGER, target, intent(inout) ::  B(*)
        INTEGER, target, intent(inout) :: C(*)
        INTEGER, intent(in) ::  LDB, LDC
C       arrays of pointers:
C       RHS(p)%val is the copy of B in precision p
C       QB(p)%val is the temporary array used to store FRB%Q_calc(p)*B
        type(LRMIXED_TYPE), intent(inout) :: RHS(*)
        INTEGER, TARGET, INTENT(IN) :: QB(*)
C       the caller should allocate this array which 
C       is supposed enough to store
C       any Q_stoc/R_stoc from its precision to the next working 
C       precision
        integer, intent(in) :: tws_size
        DOUBLE PRECISION, target, intent(inout) :: 
     &          thread_workspace(tws_size)
!        type(c_ptr), intent(inout) :: 
!     &    thread_ptr
C       Local variables
        INTEGER :: PREC, PREC_CALC, PCALC_STOC
        INTEGER(kind=1), pointer :: Aptr(:), Bptr(:)
        DOUBLE PRECISION, pointer :: Cptr(:)
        integer :: i, j, ld_c, ld_b
C       mixed blas are not "blas specification"-compliant so the
C       declaration of DALPHA is needed
        DOUBLE PRECISION :: DALPHA 
        DOUBLE PRECISION :: gamma
        integer :: ptr_shape(1)
        IF (M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0) THEN
          RETURN
        ENDIF
        DALPHA = DBLE(ALPHA)
        CALL FRB_GET_PREC(FRB, PREC, .true.)
        do prec_calc=NPREC_CALC,1,-1
          if (calc2stoc(prec_calc).le.prec) then
            exit 
          end if
        end do
        prec_calc = min(nprec_calc, prec_calc)
        pcalc_stoc = calc2stoc(prec_calc)
        ptr_shape(1) = COSTS_PREC(PREC)*FRB%M*FRB%N
        if (prec.le.1) then
          call c_f_pointer(c_loc(FRB%Q), Aptr, ptr_shape)
        else if (calc2stoc(prec_calc).eq.prec) then
          call c_f_pointer(c_loc(FRB%Q_stoc(prec)%val),
     &                     Aptr, ptr_shape)
        else
          ! TODO casting A(1:i,:) is unnecessary
          ptr_shape(1) = COSTS_PREC(1)*FRB%M*FRB%N
          call c_f_pointer(c_loc(thread_workspace),
     &      Aptr, ptr_shape)
          CALL DMUMPS_LRMIXED_CAST(
     &      FRB%Q_stoc(PREC)%VAL(1),
     &      Aptr(1), FRB%M, FRB%N,
     &      LDA, LDA,
     &      PREC, pcalc_stoc)
        ENDIF
        if (prec_calc.le.1) then
          LD_B = LDB
          ptr_shape(1) = LD_B*COSTS_PREC(PREC)*N
          call c_f_pointer(c_loc(B), Bptr, ptr_shape)
          LD_C = LDC
          ptr_shape(1) = LD_C*COSTS_PREC(PREC)*N
          call c_f_pointer(c_loc(C), Cptr, ptr_shape)
          gamma = 1.0d0
        else
          LD_B = K
          ptr_shape(1) = LD_B*COSTS_PREC(PCALC_STOC) * N
          call c_f_pointer(c_loc(RHS(prec_calc)%val), Bptr, ptr_shape)
          LD_C = M
          ptr_shape(1) = LD_C*COSTS_PREC(PCALC_STOC) * N
          call c_f_pointer(c_loc(QB), Cptr, ptr_shape)
          gamma = 0.0d0
        ENDIF
        i = (ia-1) * costs_prec(pcalc_stoc) + 1
        j = (ja-1) * costs_prec(pcalc_stoc) + 1
        CALL DMUMPS_LRMIXED_GEMM(TRANSA, TRANSB,
     &      M, N, K,
     &      DALPHA, Aptr(i-1+LDA*COSTS_PREC(PCALC_STOC)*(j-1)+1), LDA, 
     &             Bptr(1), LD_B, 
     &      GAMMA, Cptr(1), LD_C, PCALC_STOC)
        if (PCALC_STOC.GT.1) then
           CALL DMUMPS_LRMIXED_PLUS(
     &         Cptr(1), C, M, N, LD_C, LDC, pcalc_stoc, CALC2STOC(1))
        endif 
      END SUBROUTINE DMUMPS_LRMIXED_FRxFR
#endif
      SUBROUTINE DMUMPS_LRMIXED_GEMM(TRANSA, TRANSB, M, N, K, ALPHA,
     &                   A, LDA, B, LDB, BETA, C, LDC, prec)
        USE DMUMPS_LR_TYPE
        INTEGER :: A(*), B(*), C(*)
        INTEGER :: M, N, K, LDA, LDB, LDC, prec
        DOUBLE PRECISION :: ALPHA, BETA
        CHARACTER(len=1) :: TRANSA, TRANSB
        IF (prec.EQ.prec_dbl) THEN
            CALL MUMPS_LRMIXED_DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, 
     &           A, LDA, B, LDB, BETA, C, LDC)
        ELSEIF (prec.EQ.prec_sgl) THEN
            CALL MUMPS_LRMIXED_SGEMM(TRANSA, TRANSB, M, N, K, ALPHA, 
     &           A, LDA, B, LDB, BETA, C, LDC)
        ELSE 
          write(*,*) 'Internal error in MUMPS_LRMIXED_GEMM, prec=',prec
          CALL MUMPS_ABORT()
        ENDIF
      END SUBROUTINE DMUMPS_LRMIXED_GEMM
