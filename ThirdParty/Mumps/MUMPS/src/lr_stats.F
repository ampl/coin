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
      MODULE MUMPS_LR_STATS
      USE MUMPS_LR_COMMON, only: NPREC_STOC_MAX,
     &                           NPREC_CALC, NPREC_STOC, 
     &                           COSTS_PREC, PRECISIONS, CALC2STOC,
     &                           PREC_DBL, PREC_SGL, 
     &                           PREC_56, PREC_48,
     &                           PREC_40, PREC_24,
     &                           PREC_16,
     &                           PREC_FR, PREC_DIAG_FR, CLOSE_PREC,
     &                           PREC_FAC,
     &                           QUANTIZE
      IMPLICIT NONE
C --------------------------------------------------------
C  ALL the variables for the statistics
C --------------------------------------------------------
C memory stuff
      DOUBLE PRECISION :: MRY_CB_FR, 
     &                    MRY_CB_LRGAIN,
     &                    MRY_LU_FR,
     &                    MRY_LU_LRGAIN,
     &                    GLOBAL_MRY_LPRO_COMPR,
     &                    GLOBAL_MRY_LTOT_COMPR
     &                  , MRY_CB_LRGAIN_mixed,
     &                    MRY_LU_LRGAIN_mixed,
     &                    GLOBAL_MRY_LPRO_COMPR_mixed,
     &                    GLOBAL_MRY_LTOT_COMPR_mixed
                   ! per_prec = off-diagonal LR blocks in LR fronts
     &                  , ENTRIES_PER_PREC(NPREC_STOC_MAX),
                   ! per_prec_fr = off-diagonal FR blocks in LR fronts
     &                    ENTRIES_PER_PREC_FR(NPREC_STOC_MAX),
                   ! per_prec_diag = FR diagonal blocks in LR fronts
     &                    ENTRIES_PER_PREC_DIAG(NPREC_STOC_MAX),
                          ! per_prec_cb = idk
     &                    ENTRIES_PER_PREC_CB(NPREC_STOC_MAX)
      INTEGER :: CNT_NODES
C flop stuff
      DOUBLE PRECISION :: FLOP_LRGAIN,
     &                    FLOP_FACTO_FR,
     &                    FLOP_FACTO_LR,
     &                    FLOP_PANEL,
     &                    FLOP_TRSM,
     &                    FLOP_TRSM_FR,
     &                    FLOP_TRSM_LR,
     &                    FLOP_UPDATE_FR,
     &                    FLOP_UPDATE_LR,
     &                    FLOP_UPDATE_LRLR1,
     &                    FLOP_UPDATE_LRLR2,
     &                    FLOP_UPDATE_LRLR3,
     &                    FLOP_UPDATE_FRLR,
     &                    FLOP_UPDATE_FRFR
     &                  , FLOP_LRGAIN_mixed,
     &                    FLOP_FACTO_LR_mixed,
     &                    FLOP_TRSM_LR_mixed,
     &                    FLOP_UPDATE_LR_mixed,
     &                    FLOP_UPDATE_LRLR1_mixed,
     &                    FLOP_UPDATE_LRLR2_mixed,
     &                    FLOP_UPDATE_LRLR3_mixed,
     &                    FLOP_UPDATE_FRLR_mixed,
     &                    FLOP_UPDATE_FRFR_mixed
      DOUBLE PRECISION :: FLOP_COMPRESS,
     &                    FLOP_CB_COMPRESS,
     &                    FLOP_COMPRESS_mixed,
     &                    FLOP_MIDBLK_COMPRESS,
     &                    FLOP_FRSWAP_COMPRESS,
     &                    FLOP_ACCUM_COMPRESS,
     &                    FLOP_DECOMPRESS,
     &                    FLOP_CB_DECOMPRESS,
     &                    FLOP_FRFRONTS
      DOUBLE PRECISION :: FACTOR_PROCESSED_FRACTION
      INTEGER(KIND=8)  :: FACTOR_SIZE
      DOUBLE PRECISION :: TOTAL_FLOP
C timing stuff
      DOUBLE PRECISION :: TIME_UPDATE
      DOUBLE PRECISION :: TIME_UPDATE_LRLR1
      DOUBLE PRECISION :: TIME_UPDATE_LRLR2
      DOUBLE PRECISION :: TIME_UPDATE_LRLR3
      DOUBLE PRECISION :: TIME_UPDATE_FRLR
      DOUBLE PRECISION :: TIME_UPDATE_FRFR
      DOUBLE PRECISION :: TIME_COMPRESS
      DOUBLE PRECISION :: TIME_MIDBLK_COMPRESS
      DOUBLE PRECISION :: TIME_FRSWAP_COMPRESS
      DOUBLE PRECISION :: TIME_CB_COMPRESS
      DOUBLE PRECISION :: TIME_LR_MODULE
      DOUBLE PRECISION :: TIME_UPD_NELIM
      DOUBLE PRECISION :: TIME_LRTRSM
      DOUBLE PRECISION :: TIME_FRTRSM
      DOUBLE PRECISION :: TIME_PANEL
      DOUBLE PRECISION :: TIME_FAC_I
      DOUBLE PRECISION :: TIME_FAC_MQ
      DOUBLE PRECISION :: TIME_FAC_SQ
      DOUBLE PRECISION :: TIME_FRFRONTS
      DOUBLE PRECISION :: TIME_DIAGCOPY
      DOUBLE PRECISION :: TIME_DECOMP
      DOUBLE PRECISION :: TIME_DECOMP_UCFS
      DOUBLE PRECISION :: TIME_LRASM_NIV1
      DOUBLE PRECISION :: TIME_LRASM_LOCASM2
      DOUBLE PRECISION :: TIME_LRASM_MAPLIG1
      DOUBLE PRECISION :: TIME_LRASM_CONTRIB2
      DOUBLE PRECISION :: TIME_FRASM_LOCASM2
      DOUBLE PRECISION :: TIME_FRASM_MAPLIG1
      DOUBLE PRECISION :: TIME_FRASM_CONTRIB2
C analysis timers      
      DOUBLE PRECISION :: TIME_LRANA_LRGROUPING
      DOUBLE PRECISION :: TIME_LRANA_SEPGROUPING
      DOUBLE PRECISION :: TIME_LRANA_GETHALO
      DOUBLE PRECISION :: TIME_LRANA_KWAY
      DOUBLE PRECISION :: TIME_LRANA_GNEW
      DOUBLE PRECISION :: AVG_FLOP_FACTO_LR
      DOUBLE PRECISION :: MIN_FLOP_FACTO_LR
      DOUBLE PRECISION :: MAX_FLOP_FACTO_LR
      INTEGER :: TOTAL_NBLOCKS_ASS, TOTAL_NBLOCKS_CB
      INTEGER :: MIN_BLOCKSIZE_ASS, MAX_BLOCKSIZE_ASS
      INTEGER :: MIN_BLOCKSIZE_CB, MAX_BLOCKSIZE_CB
      DOUBLE PRECISION :: AVG_BLOCKSIZE_ASS, AVG_BLOCKSIZE_CB
      CONTAINS
      SUBROUTINE COLLECT_BLOCKSIZES(CUT,NPARTSASS,NPARTSCB)
        INTEGER, INTENT(IN) :: NPARTSASS, NPARTSCB
        INTEGER, POINTER, DIMENSION(:) :: CUT
        INTEGER :: LOC_MIN_ASS, LOC_MIN_CB, LOC_MAX_ASS, LOC_MAX_CB,
     &             LOC_TOT_ASS, LOC_TOT_CB
        DOUBLE PRECISION :: LOC_AVG_ASS, LOC_AVG_CB 
        INTEGER :: I
        LOC_TOT_ASS = 0
        LOC_TOT_CB = 0
        LOC_AVG_ASS = 0.D0
        LOC_AVG_CB = 0.D0
        LOC_MIN_ASS = 100000
        LOC_MIN_CB = 100000
        LOC_MAX_ASS = 0
        LOC_MAX_CB = 0
        DO I = 1,NPARTSASS
          LOC_AVG_ASS = ( LOC_TOT_ASS * LOC_AVG_ASS
     &                        + CUT(I+1) - CUT(I) )
     &                        / (LOC_TOT_ASS + 1)
          LOC_TOT_ASS = LOC_TOT_ASS + 1
          IF (CUT(I+1) - CUT(I) .LE. LOC_MIN_ASS) THEN
            LOC_MIN_ASS = CUT(I+1) - CUT(I)
          END IF
          IF (CUT(I+1) - CUT(I) .GE. LOC_MAX_ASS) THEN
            LOC_MAX_ASS = CUT(I+1) - CUT(I)
          END IF
        END DO
        DO I = NPARTSASS+1,NPARTSASS+NPARTSCB
          LOC_AVG_CB = ( LOC_TOT_CB * LOC_AVG_CB
     &                        + CUT(I+1) - CUT(I) )
     &                        / (LOC_TOT_CB + 1)
          LOC_TOT_CB = LOC_TOT_CB + 1
          IF (CUT(I+1) - CUT(I) .LE. LOC_MIN_CB) THEN
            LOC_MIN_CB = CUT(I+1) - CUT(I)
          END IF
          IF (CUT(I+1) - CUT(I) .GE. LOC_MAX_CB) THEN
            LOC_MAX_CB = CUT(I+1) - CUT(I)
          END IF
        END DO
        AVG_BLOCKSIZE_ASS = (TOTAL_NBLOCKS_ASS*AVG_BLOCKSIZE_ASS
     &     + LOC_TOT_ASS*LOC_AVG_ASS) / (TOTAL_NBLOCKS_ASS+LOC_TOT_ASS)
        AVG_BLOCKSIZE_CB = (TOTAL_NBLOCKS_CB*AVG_BLOCKSIZE_CB
     &     + LOC_TOT_CB*LOC_AVG_CB) / (TOTAL_NBLOCKS_CB+LOC_TOT_CB)
        TOTAL_NBLOCKS_ASS = TOTAL_NBLOCKS_ASS + LOC_TOT_ASS
        TOTAL_NBLOCKS_CB = TOTAL_NBLOCKS_CB + LOC_TOT_CB
        MIN_BLOCKSIZE_ASS = min(MIN_BLOCKSIZE_ASS,LOC_MIN_ASS)
        MIN_BLOCKSIZE_CB = min(MIN_BLOCKSIZE_CB,LOC_MIN_CB)
        MAX_BLOCKSIZE_ASS = max(MAX_BLOCKSIZE_ASS,LOC_MAX_ASS)
        MAX_BLOCKSIZE_CB = max(MAX_BLOCKSIZE_CB,LOC_MAX_CB)
C        open(unit=20,access='append',file='block_stats.txt')
C        write(20,*) LOC_TOT_ASS, LOC_AVG_ASS,
C     &              LOC_MIN_ASS, LOC_MAX_ASS,
C     &              LOC_TOT_CB, LOC_AVG_CB,
C     &              LOC_MIN_CB, LOC_MAX_CB
C        close(20)
      END SUBROUTINE COLLECT_BLOCKSIZES
      SUBROUTINE UPD_FLOP_DECOMPRESS(F, CB)
          DOUBLE PRECISION, INTENT(IN) :: F
          LOGICAL, INTENT(IN) :: CB
!$OMP     ATOMIC UPDATE
          FLOP_DECOMPRESS = FLOP_DECOMPRESS + F
!$OMP     END ATOMIC
          IF (CB) THEN
!$OMP       ATOMIC UPDATE
            FLOP_CB_DECOMPRESS = FLOP_CB_DECOMPRESS + F
!$OMP       END ATOMIC
          ENDIF
          RETURN
      END SUBROUTINE UPD_FLOP_DECOMPRESS
      SUBROUTINE UPD_FLOP_COMPRESS(LRBM,LRBN,LRBK,ISLR, REC_ACC, 
     &         CB_COMPRESS, FRSWAP)
        INTEGER :: LRBM, LRBN, LRBK
        LOGICAL :: ISLR
        INTEGER(8) :: M,N,K
        DOUBLE PRECISION :: HR_COST,BUILDQ_COST,
     &  HR_AND_BUILDQ_COST
        LOGICAL, OPTIONAL :: REC_ACC, CB_COMPRESS, FRSWAP
        M = int(LRBM,8)
        N = int(LRBN,8)
        K = int(LRBK,8)
        HR_COST =  dble(K*K*K/3_8 + 4_8*K*M*N - (2_8*M+N)*K*K)
        IF (ISLR) THEN 
          BUILDQ_COST = dble(2_8*K*K*M - K*K*K)
        ELSE 
          BUILDQ_COST = 0.0d0
        END IF
        HR_AND_BUILDQ_COST = HR_COST + BUILDQ_COST
!$OMP   ATOMIC UPDATE
        FLOP_COMPRESS = FLOP_COMPRESS + HR_AND_BUILDQ_COST
!$OMP   END ATOMIC
        IF (present(REC_ACC)) THEN
          IF (REC_ACC) THEN
!$OMP       ATOMIC UPDATE
            FLOP_ACCUM_COMPRESS = FLOP_ACCUM_COMPRESS +
     &                            HR_AND_BUILDQ_COST
!$OMP       END ATOMIC
          ENDIF
        ENDIF
        IF (present(CB_COMPRESS)) THEN
          IF (CB_COMPRESS) THEN
!$OMP       ATOMIC UPDATE
            FLOP_CB_COMPRESS = FLOP_CB_COMPRESS +
     &                         HR_AND_BUILDQ_COST 
!$OMP       END ATOMIC
          ENDIF
        ENDIF
        IF (present(FRSWAP)) THEN
          IF (FRSWAP) THEN
!$OMP       ATOMIC UPDATE
            FLOP_FRSWAP_COMPRESS = FLOP_FRSWAP_COMPRESS +
     &                             HR_AND_BUILDQ_COST
!$OMP       END ATOMIC
          ENDIF
        ENDIF
      RETURN
      END SUBROUTINE UPD_FLOP_COMPRESS
      SUBROUTINE UPD_FLOP_TRSM(M, N, K, ISLR, LorU)
          INTEGER,INTENT(IN) :: M, N, K, LorU
          LOGICAL,INTENT(IN) :: ISLR
          DOUBLE PRECISION :: LR_COST, FR_COST, LR_GAIN
          IF (LorU.EQ.0) THEN ! L
            FR_COST = dble(M*N*N)
            IF (ISLR) THEN
              LR_COST = dble(K*N*N)
            ELSE
              LR_COST = FR_COST
            ENDIF
          ELSE ! U
            FR_COST = dble(M-1)*dble(N*N)
            IF (ISLR) THEN
              LR_COST = dble(N-1)*dble(N*K)
            ELSE
              LR_COST = FR_COST
            ENDIF
          ENDIF
          LR_GAIN = FR_COST - LR_COST
!$OMP     ATOMIC UPDATE
          FLOP_LRGAIN  = FLOP_LRGAIN + LR_GAIN
!$OMP     END ATOMIC
      RETURN
      END SUBROUTINE UPD_FLOP_TRSM
      SUBROUTINE UPD_FLOP_UPDATE(LRB1M, LRB1N, LRB1K, LRB1ISLR, 
     &      LRB2M, LRB2N, LRB2K, LRB2ISLR, 
     &      MIDBLK_COMPRESS, RANK_IN, BUILDQ,
     &      IS_SYMDIAG, LUA_ACTIVATED, REC_ACC)
        INTEGER, INTENT(IN) :: LRB1M, LRB1N, LRB1K, LRB2M, LRB2N, LRB2K
        LOGICAL, INTENT(IN) :: LRB1ISLR, LRB2ISLR
C       MIDBLK_COMPRESS = whether we tried compressing the middle
C       product in CMUMPS_LRGEMM3  
C       RANK = rank revealed when MIDBLK_COMPRESS >= 1
C       BUILDQ = if RANK is small enough, we also build Q in
C       CMUMPS_LRGEMM3 to exploit low-rankness
C       IS_SYMDIAG = whether it is a diagonal block AND we are in
C                    the symmetric case 
C       REC_ACC = whether this routine is called from
C       CMUMPS_RECOMPRESS_ACCUMULATOR        
        LOGICAL, INTENT(IN) :: BUILDQ, IS_SYMDIAG, LUA_ACTIVATED
        INTEGER, INTENT(IN) :: RANK_IN, MIDBLK_COMPRESS
        LOGICAL, INTENT(IN), OPTIONAL :: REC_ACC
        DOUBLE PRECISION :: COST_FR, COST_LR, COST_LRLR1, COST_LRLR2,
     &                      COST_LRLR3, COST_FRLR, COST_FRFR, 
     &                      COST_COMPRESS, COST_LR_AND_COMPRESS, LR_GAIN
        DOUBLE PRECISION :: M1,N1,K1,M2,N2,K2,RANK
        LOGICAL :: REC_ACC_LOC
        M1 = dble(LRB1M)
        N1 = dble(LRB1N)
        K1 = dble(LRB1K)
        M2 = dble(LRB2M)
        N2 = dble(LRB2N)
        K2 = dble(LRB2K)
        RANK = dble(RANK_IN)
        COST_LRLR1 = 0.0D0
        COST_LRLR2 = 0.0D0
        COST_LRLR3 = 0.0D0
        COST_FRLR = 0.0D0
        COST_FRFR = 0.0D0
        COST_COMPRESS = 0.0D0
        IF (present(REC_ACC)) THEN
          REC_ACC_LOC = REC_ACC
        ELSE
          REC_ACC_LOC = .FALSE.
        ENDIF
        IF ((.NOT.LRB1ISLR).AND.(.NOT.LRB2ISLR)) THEN
          COST_FRFR = 2.0D0*M1*M2*N1
          COST_LR = 2.0D0*M1*M2*N1
          COST_FR = 2.0D0*M1*M2*N1
        ELSEIF (LRB1ISLR.AND.(.NOT.LRB2ISLR)) THEN
          COST_FRLR = 2.0D0*K1*M2*N1
          COST_LRLR3 = 2.0D0*M1*M2*K1
          COST_LR = COST_FRLR + COST_LRLR3
          COST_FR = 2.0D0*M1*M2*N1
        ELSEIF ((.NOT.LRB1ISLR).AND.LRB2ISLR) THEN
          COST_FRLR = 2.0D0*M1*K2*N1 
          COST_LRLR3 = 2.0D0*M1*M2*K2
          COST_LR = COST_FRLR + COST_LRLR3
          COST_FR = 2.0D0*M1*M2*N1
        ELSE
          IF (MIDBLK_COMPRESS.GE.1) THEN
            COST_COMPRESS =  RANK*RANK*RANK/3.0D0 + 
     &                       4.0D0*RANK*K1*K2 - 
     &                       (2.0D0*K1+K2)*RANK*RANK
            IF (BUILDQ) THEN
              COST_COMPRESS = COST_COMPRESS + 4.0D0*RANK*RANK*K1 
     &                                      - RANK*RANK*RANK
            ENDIF
          ENDIF
          COST_LRLR1 = 2.0D0*K1*K2*N1
          IF ((MIDBLK_COMPRESS.GE.1).AND.BUILDQ) THEN
            COST_LRLR2 = 2.0D0*K1*M1*RANK + 2.0D0*K2*M2*RANK
            COST_LRLR3 = 2.0D0*M1*M2*RANK
          ELSE
            IF (K1 .GE. K2) THEN
              COST_LRLR2 = 2.0D0*K1*M1*K2
              COST_LRLR3 = 2.0D0*M1*M2*K2
            ELSE
              COST_LRLR2 = 2.0D0*K1*M2*K2
              COST_LRLR3 = 2.0D0*M1*M2*K1
            ENDIF
          ENDIF
          COST_LR = COST_LRLR1 + COST_LRLR2 + COST_LRLR3
          COST_FR = 2.0D0*M1*M2*N1
        ENDIF
        IF (IS_SYMDIAG) THEN
C         We do not use the fact that diagonal blocks are triangular 
C         in LDLT. But the total flops are computed assuming we do.            
C         If we want the sum of each part to be equal to the total,
C         we have to cheat:            
          COST_FR = COST_FR/2.0D0
          COST_LRLR3 = COST_LRLR3/2.0D0
          COST_FRFR = COST_FRFR/2.0D0
          COST_LR = COST_LR - COST_LRLR3 - COST_FRFR
        ENDIF
        IF (LUA_ACTIVATED) THEN
C         LUA is activated: the outer update (LR block decompression)
C         is not done
          COST_LR = COST_LR - COST_LRLR3
          COST_LRLR3 = 0.0D0
          IF (REC_ACC_LOC) THEN
C           Special case: this a not a regular LR product, it is
C           called from CMUMPS_RECOMPRESS_ACCUMULATOR 
            COST_LR_AND_COMPRESS = COST_LR + COST_COMPRESS
!$OMP       ATOMIC UPDATE
            FLOP_COMPRESS  = FLOP_COMPRESS + COST_LR_AND_COMPRESS
!$OMP       END ATOMIC
          ENDIF
        ENDIF
        IF (.NOT.REC_ACC_LOC) THEN
!$OMP     ATOMIC UPDATE
          FLOP_COMPRESS  = FLOP_COMPRESS  + COST_COMPRESS
!$OMP     END ATOMIC
          LR_GAIN = COST_FR - COST_LR
!$OMP     ATOMIC UPDATE
          FLOP_LRGAIN = FLOP_LRGAIN + LR_GAIN
!$OMP     END ATOMIC
        ENDIF
      END SUBROUTINE UPD_FLOP_UPDATE
      SUBROUTINE UPD_FLOP_UPDATE_LRLR3(M, N, K, NIV)
        INTEGER,INTENT(IN) :: M, N, K, NIV
        DOUBLE PRECISION :: FLOP_COST
        FLOP_COST = 2.0D0*dble(M)*dble(N)*dble(K)
!$OMP   ATOMIC UPDATE
        FLOP_LRGAIN = FLOP_LRGAIN - FLOP_COST
!$OMP   END ATOMIC
        RETURN
      END SUBROUTINE UPD_FLOP_UPDATE_LRLR3
      SUBROUTINE UPD_FLOP_ROOT(KEEP50, NFRONT, NPIV,
     &           NPROW, NPCOL, MYID)
C
C  Purpose:
C  =======
C
C     Update FLOP_FRFRONTS with flops for the factorization of 
C     root node.
C     To be called only by processes participating to the root node
C     (root%yes = .true.)
C
C  Arguments:
C  =========
        INTEGER, intent(in) :: KEEP50, NFRONT, NPIV,
     &           NPROW, NPCOL, MYID
C  Local variables
C  ===============
        DOUBLE PRECISION :: COST, COST_PER_PROC
        INTEGER, PARAMETER :: LEVEL3 = 3
        CALL MUMPS_GET_FLOPS_COST(NFRONT, NPIV, NFRONT, KEEP50, LEVEL3, 
     &                            COST)
        COST_PER_PROC = dble(int( COST,8) / int(NPROW * NPCOL,8))
C      Simpler for root to forget a few flops and have same 
C      value on all procs.
C      Proc with MYID = 0 gets more flops in case division by NPROW and
C      NPCOL is not exact.
C      IF ( MYID .EQ. 0) THEN
C       COST_PER_PROC = COST - (COST_PER_PROC * dble(NPROW * NPCOL-1))
C      ENDIF
!$OMP   ATOMIC UPDATE
        FLOP_FRFRONTS = FLOP_FRFRONTS + COST_PER_PROC
!$OMP   END ATOMIC
        RETURN
      END SUBROUTINE UPD_FLOP_ROOT
C  Stats front-by-front no longer available
      SUBROUTINE INIT_STATS_GLOBAL()
        INTEGER::PREC
        DO PREC=1,NPREC_STOC_MAX
          ENTRIES_PER_PREC(PREC)        = 0.D0
          ENTRIES_PER_PREC_FR(PREC)     = 0.D0
          ENTRIES_PER_PREC_DIAG(PREC)   = 0.D0
          ENTRIES_PER_PREC_CB(PREC)     = 0.D0
        ENDDO
        MRY_LU_LRGAIN_mixed = 0.D0
        MRY_CB_LRGAIN_mixed = 0.D0
        MRY_LU_FR = 0.D0
        MRY_LU_LRGAIN = 0.D0
        MRY_CB_FR = 0.D0
        MRY_CB_LRGAIN = 0.D0
        FLOP_FACTO_FR = 0.D0 
        FLOP_FACTO_LR = 0.D0
        FLOP_LRGAIN = 0.D0
        FLOP_CB_COMPRESS = 0.D0
        FLOP_CB_DECOMPRESS = 0.D0
        FLOP_DECOMPRESS = 0.D0
        FLOP_UPDATE_FR = 0.D0 
        FLOP_UPDATE_LR = 0.D0
        FLOP_UPDATE_LRLR1 = 0.D0
        FLOP_UPDATE_LRLR2 = 0.D0
        FLOP_UPDATE_LRLR3 = 0.D0
        FLOP_UPDATE_FRLR = 0.D0
        FLOP_UPDATE_FRFR = 0.D0
        FLOP_MIDBLK_COMPRESS = 0.D0
        FLOP_TRSM_FR = 0.D0 
        FLOP_TRSM_LR = 0.D0
        FLOP_COMPRESS = 0.D0
        FLOP_ACCUM_COMPRESS = 0.D0
        FLOP_FRSWAP_COMPRESS = 0.D0
        FLOP_PANEL = 0.D0
        FLOP_TRSM = 0.D0
        FLOP_FRFRONTS = 0.D0
        TOTAL_NBLOCKS_ASS = 0
        TOTAL_NBLOCKS_CB = 0
        AVG_BLOCKSIZE_ASS = 0.D0
        AVG_BLOCKSIZE_CB = 0.D0
        MIN_BLOCKSIZE_ASS = huge(1)
        MAX_BLOCKSIZE_ASS = 0
        MIN_BLOCKSIZE_CB = huge(1)
        MAX_BLOCKSIZE_CB = 0
        CNT_NODES = 0
        TIME_UPDATE = 0.D0 
        TIME_MIDBLK_COMPRESS = 0.D0 
        TIME_UPDATE_LRLR1 = 0.D0 
        TIME_UPDATE_LRLR2 = 0.D0 
        TIME_UPDATE_LRLR3 = 0.D0 
        TIME_UPDATE_FRLR = 0.D0 
        TIME_UPDATE_FRFR = 0.D0 
        TIME_COMPRESS = 0.D0 
        TIME_CB_COMPRESS = 0.D0 
        TIME_LR_MODULE = 0.D0 
        TIME_UPD_NELIM = 0.D0 
        TIME_LRTRSM = 0.D0 
        TIME_FRTRSM = 0.D0 
        TIME_PANEL = 0.D0 
        TIME_FAC_I = 0.D0 
        TIME_FAC_MQ = 0.D0 
        TIME_FAC_SQ = 0.D0 
        TIME_FRFRONTS = 0.D0 
        TIME_DIAGCOPY = 0.D0 
        TIME_FRSWAP_COMPRESS = 0.D0 
        TIME_DECOMP = 0.D0 
        TIME_DECOMP_UCFS = 0.D0 
        TIME_LRASM_NIV1 = 0.D0 
        TIME_LRASM_LOCASM2 = 0.D0 
        TIME_LRASM_MAPLIG1 = 0.D0 
        TIME_LRASM_CONTRIB2 = 0.D0 
        TIME_FRASM_LOCASM2 = 0.D0 
        TIME_FRASM_MAPLIG1 = 0.D0 
        TIME_FRASM_CONTRIB2 = 0.D0 
      END SUBROUTINE INIT_STATS_GLOBAL
      SUBROUTINE UPD_MRY_LU_FR(NASS, NCB, SYM, NELIM)
        INTEGER,INTENT(IN) :: NASS, NCB, SYM, NELIM
        DOUBLE PRECISION :: MRY
        INTEGER :: NPIV
        NPIV = NASS - NELIM
        IF (SYM .GT. 0) THEN
           MRY = dble(NPIV)*(dble(NPIV)+1.D0)/2.D0 
     &         + dble(NPIV)*dble(NCB+NELIM)
        ELSE
           MRY = dble(NPIV)*dble(NPIV) 
     &         + 2.0D0*dble(NPIV)*dble(NCB+NELIM)
        END IF
!$OMP   ATOMIC UPDATE
        MRY_LU_FR  = MRY_LU_FR + MRY
!$OMP   END ATOMIC
      RETURN
      END SUBROUTINE UPD_MRY_LU_FR
      SUBROUTINE UPD_MRY_CB_FR(NROWS, NCOLS, SYM)
        INTEGER,INTENT(IN) :: NROWS, NCOLS, SYM
        DOUBLE PRECISION :: MRY
        IF (SYM.EQ.0) THEN
          MRY = dble(NCOLS)*dble(NROWS)
        ELSE
          MRY = dble(NCOLS-NROWS)*dble(NROWS) +
     &                dble(NROWS)*dble(NROWS+1)/2.D0
        ENDIF
!$OMP   ATOMIC UPDATE
        MRY_CB_FR = MRY_CB_FR + MRY
!$OMP   END ATOMIC
        RETURN
      END SUBROUTINE UPD_MRY_CB_FR
      SUBROUTINE UPD_MRY_CB_LRGAIN(M, N, K
     &            , LRMIXED_ACTIVATED, K_STOC, NPREC_STOC, COSTS_PREC
     &                            )
        INTEGER, INTENT(IN) :: M, N, K 
        LOGICAL, INTENT(IN) :: LRMIXED_ACTIVATED
        INTEGER, INTENT(IN) :: NPREC_STOC
C       INTEGER, INTENT(IN) :: K_STOC(NPREC_STOC)
C       K_STOC is declared assumed shape array to allow passing a
C       pointer array without the risk of temporary array creation
        INTEGER, INTENT(IN) :: K_STOC(:)
        INTEGER, INTENT(IN) :: COSTS_PREC(NPREC_STOC_MAX)
C       Local variables
        DOUBLE PRECISION :: LRGAIND
        DOUBLE PRECISION :: LRGAIND_mixed
        INTEGER :: PREC
        LRGAIND = dble(M*N-(M+N)*K)
!$OMP   ATOMIC UPDATE
        MRY_CB_LRGAIN = MRY_CB_LRGAIN + LRGAIND
!$OMP   END ATOMIC
        IF(LRMIXED_ACTIVATED) THEN
            LRGAIND_mixed = dble(M*N)
            DO PREC=1,NPREC_STOC
               LRGAIND_mixed = LRGAIND_mixed - dble(K_STOC(PREC))
     &          *dble(COSTS_PREC(PREC))/dble(COSTS_PREC(1))
     &          *dble(M + N)
            END DO
!$OMP   ATOMIC UPDATE
            MRY_CB_LRGAIN_mixed = MRY_CB_LRGAIN_mixed + LRGAIND_mixed
!$OMP   END ATOMIC
        ENDIF
        RETURN
      END SUBROUTINE UPD_MRY_CB_LRGAIN
      SUBROUTINE UPD_FLOP_FACTO_FR( NFRONT, NASS, NPIV, SYM, NIV)
          INTEGER,INTENT(IN) :: NFRONT, SYM, NASS, NPIV, NIV
          DOUBLE PRECISION   :: FLOP
          CALL MUMPS_GET_FLOPS_COST(NFRONT, NPIV, NASS, 
     &                              SYM, NIV, FLOP)
!$OMP     ATOMIC UPDATE
          FLOP_FACTO_FR = FLOP_FACTO_FR + FLOP
!$OMP     END ATOMIC
      END SUBROUTINE UPD_FLOP_FACTO_FR
      SUBROUTINE STATS_COMPUTE_FLOP_SLAVE_TYPE2( NROW1, NCOL1,
     &                NASS1, KEEP50, INODE)
C         This routines computes the FR flops performed by a slave on
C         its band of a type 2 front
          INTEGER,INTENT(IN) :: NROW1, NCOL1, KEEP50, NASS1, INODE
          DOUBLE PRECISION   :: NROW2, NCOL2, NASS2
          DOUBLE PRECISION   :: FLOP
          NROW2 = dble(NROW1)
          NCOL2 = dble(NCOL1)
          NASS2 = dble(NASS1)
          IF (KEEP50.EQ.0) THEN
            FLOP = NROW2*NASS2*NASS2         ! trsm
     &              + 2.0D0*NROW2*NASS2*(NCOL2-NASS2) ! gemm
          ELSE
             FLOP =
C              trsm
     &            NROW2*NASS2*NASS2
C              gemm diag
     &          + NROW2*NASS2*NROW2
C              gemm subdiag
     &          + 2.0D0*NROW2*NASS2*(NCOL2-NASS2-NROW2)
          ENDIF
!$OMP     ATOMIC UPDATE
          FLOP_FACTO_FR = FLOP_FACTO_FR + FLOP
!$OMP     END ATOMIC
      RETURN
      END SUBROUTINE STATS_COMPUTE_FLOP_SLAVE_TYPE2
      SUBROUTINE UPD_FLOP_FRFRONTS(NFRONT, NPIV, NASS, SYM, 
     &                                        NIV)
          INTEGER, INTENT(IN) :: NFRONT, NPIV, NASS, SYM, NIV
          DOUBLE PRECISION    :: FLOP_FAC
          CALL MUMPS_GET_FLOPS_COST(NFRONT, NPIV, NASS, 
     &                              SYM, NIV, FLOP_FAC)
!$OMP     ATOMIC UPDATE
          FLOP_FRFRONTS = FLOP_FRFRONTS + FLOP_FAC
!$OMP     END ATOMIC
      RETURN
      END SUBROUTINE UPD_FLOP_FRFRONTS
      SUBROUTINE UPD_FLOP_FRFRONT_SLAVE(NROW1, NCOL1, NASS1,
     &                                        KEEP50, INODE)
          INTEGER,INTENT(IN) :: NROW1, NCOL1, KEEP50, NASS1, INODE
          DOUBLE PRECISION   :: NROW2, NCOL2, NASS2
          DOUBLE PRECISION   :: FLOP
          NROW2 = dble(NROW1)
          NCOL2 = dble(NCOL1)
          NASS2 = dble(NASS1)
          IF (KEEP50.EQ.0) THEN
            FLOP = NROW2*NASS2*NASS2         ! trsm
     &              + 2.0D0*NROW2*NASS2*(NCOL2-NASS2) ! gemm
          ELSE
             FLOP =
C                   trsm
     &            NROW2*NASS2*NASS2
C                   gemm diag
     &          + NROW2*NASS2*NROW2
C                   gemm subdiag
     &          + 2.0D0*NROW2*NASS2*(NCOL2-NASS2-NROW2)
          ENDIF
!$OMP     ATOMIC UPDATE
          FLOP_FRFRONTS = FLOP_FRFRONTS + FLOP
!$OMP     END ATOMIC
      RETURN
      END SUBROUTINE UPD_FLOP_FRFRONT_SLAVE
      SUBROUTINE COMPUTE_GLOBAL_GAINS(NB_ENTRIES_FACTOR, 
     &                FLOP_NUMBER, NB_ENTRIES_FACTOR_withLR,
     &                PROKG, MPG)
        INTEGER(8), INTENT(IN) :: NB_ENTRIES_FACTOR   !=KEEP8(110)
        INTEGER, INTENT(IN)    :: MPG
        LOGICAL, INTENT(IN)    :: PROKG
        DOUBLE PRECISION, INTENT(IN)  :: FLOP_NUMBER   !=RINFOG(3)
        INTEGER(8), INTENT(OUT) :: 
     &                  NB_ENTRIES_FACTOR_withLR !=KEEP8(49)
C       Local variables:
C        DOUBLE PRECISION :: LR_PEAK, FR_PEAK
        IF (NB_ENTRIES_FACTOR < 0) THEN
         IF (PROKG.AND.MPG.GT.0) THEN
          WRITE(MPG,*) "NEGATIVE NUMBER OF ENTRIES IN FACTOR"
          WRITE(MPG,*) "===> OVERFLOW ?"
         END IF
        END IF
        IF (MRY_LU_FR .EQ. 0) THEN
           GLOBAL_MRY_LPRO_COMPR = 100.0D0
        ELSE
           GLOBAL_MRY_LPRO_COMPR = 100.0D0 *
     &                             MRY_LU_LRGAIN/MRY_LU_FR
        ENDIF
        IF (MRY_LU_FR .EQ. 0) THEN
           GLOBAL_MRY_LPRO_COMPR_mixed = 100.0D0
        ELSE
           GLOBAL_MRY_LPRO_COMPR_mixed = 100.0D0 *
     &                             MRY_LU_LRGAIN_mixed/MRY_LU_FR
        ENDIF
C if only the root was processed with LR we need to set the values by
C hand to avoid problems
        IF (MRY_CB_FR .EQ. 0) THEN 
          MRY_CB_FR = 100.0D0
        END IF
        NB_ENTRIES_FACTOR_withLR = NB_ENTRIES_FACTOR -
     &                             int(MRY_LU_LRGAIN,8)
        IF (NB_ENTRIES_FACTOR.EQ.0) THEN
          FACTOR_PROCESSED_FRACTION = 100.0D0
          GLOBAL_MRY_LTOT_COMPR = 100.0D0
          GLOBAL_MRY_LTOT_COMPR_mixed = 100.0D0
        ELSE
          FACTOR_PROCESSED_FRACTION = 100.0D0 *
     &                            MRY_LU_FR/dble(NB_ENTRIES_FACTOR)
          GLOBAL_MRY_LTOT_COMPR = 
     &            100.0D0*MRY_LU_LRGAIN/dble(NB_ENTRIES_FACTOR)
          GLOBAL_MRY_LTOT_COMPR_mixed = 
     &            100.0D0*MRY_LU_LRGAIN_mixed/dble(NB_ENTRIES_FACTOR)
        ENDIF
        TOTAL_FLOP = FLOP_NUMBER
        FLOP_FACTO_LR = FLOP_FACTO_FR - FLOP_LRGAIN + FLOP_COMPRESS 
     &                                              + FLOP_DECOMPRESS
        RETURN
      END SUBROUTINE COMPUTE_GLOBAL_GAINS
      SUBROUTINE SAVEandWRITE_GAINS(LOCAL, K489,
     &         LR_DKEEPSHIFT, LR_TABSIZE, LR_TAB, LR_EPSILON,
     &         N, ICNTL36,
     &         DEPTH, BCKSZ, NASSMIN, NFRONTMIN, SYM, K486,
     &         K450, K452, NPREC_STOC,
     &         NBTHREADS, K472, K475, K478, K480, K481, K483, K484,
     &         K117, K8110, K849, 
     &         NBTREENODES, NPROCS, MPG, PROKG)
        INTEGER, INTENT(IN) :: LOCAL,K489,DEPTH, N,
     &      ICNTL36, BCKSZ,NASSMIN,
     &      NFRONTMIN, K486, NBTREENODES, MPG, 
     &      NBTHREADS, K472, K475, K478, K480, K481, K483, K484, K117,
     &      SYM, NPROCS
C       LR_TAB is a double precision array that corrsponds
C       to a section of the real or complex array DKEEP
C       LR_TAB(xx-LR_DKEEPSHIFT) <-> DKEEP(xx)
        INTEGER, INTENT(IN) :: LR_TABSIZE, LR_DKEEPSHIFT
        DOUBLE PRECISION, INTENT(INOUT) :: LR_TAB(LR_TABSIZE)
C       LR_EPSILON=CNTL(7) or DKEEP(8): just to be able to print it
        DOUBLE PRECISION, INTENT(IN)    :: LR_EPSILON
        INTEGER, INTENT(IN) :: K450, K452, NPREC_STOC
        INTEGER(8), INTENT(IN) :: K8110, K849
        LOGICAL, INTENT(IN) :: PROKG
C  if 0 : only the global memory gain
C  if 1 : the gains for each front and the global memory gain
C       Divide timers corresponding to multithreaded computations by
C       KEEP(249) (number of threads)
        TIME_UPDATE_LRLR1    = TIME_UPDATE_LRLR1/dble(NBTHREADS)
        TIME_UPDATE_LRLR2    = TIME_UPDATE_LRLR2/dble(NBTHREADS)
        TIME_UPDATE_LRLR3    = TIME_UPDATE_LRLR3/dble(NBTHREADS)
        TIME_UPDATE_FRLR     = TIME_UPDATE_FRLR/dble(NBTHREADS)
        TIME_UPDATE_FRFR     = TIME_UPDATE_FRFR/dble(NBTHREADS)
        TIME_MIDBLK_COMPRESS = TIME_MIDBLK_COMPRESS/dble(NBTHREADS)
        IF (PROKG) THEN
        WRITE(MPG,'(/A,A)') 
     & '-------------- Beginning of BLR statistics -------------------',
     & '--------------'
        WRITE(MPG,'(A,I2)') 
     & ' ICNTL(36) BLR variant                            = ', ICNTL36
        WRITE(MPG,'(A,ES8.1)')  
     & ' CNTL(7)   Dropping parameter controlling accuracy = ',
     &                          LR_EPSILON
C        WRITE(MPG,'(A)') 
C     &          ' -----------------------------------'
        WRITE(MPG,'(A)') 
     &          ' Statistics after BLR factorization :'
C        WRITE(MPG,'(A)') 
C     &          ' -----------------------------------'
        WRITE(MPG,'(A,I8)') 
     &    '     Number of BLR fronts                     = ',
     &                          CNT_NODES
C       endif test prokg
        ENDIF  
        IF (PROKG) WRITE(MPG,'(A,F8.1,A)')
     &    '     Fraction of factors in BLR fronts        =',
     &                FACTOR_PROCESSED_FRACTION,'% ' 
        IF (PROKG) THEN
C         Factor : here we are concerned by the number of entries
C         and not the storage that also depends on the factor 
C         factors need be saved in compressed form i.e. KEEP(486)>=2
          WRITE(MPG,'(A)') 
     &  '     Statistics on the number of entries in factors :'
          WRITE(MPG,'(A,ES10.3,A)')
     &  '     INFOG(29) Theoretical nb of entries in factors      ='
     &     ,real(K8110),' (100.0%)'
          WRITE(MPG,'(A,ES10.3,A,F5.1,A)') 
     &  '     INFOG(35) Effective nb of entries  (% of INFOG(29)) ='
     &     ,real(K849),' ('
     &     ,real(100)*(real(K849)/real(max(K8110,1_8)))
     &     ,'%)'
        ENDIF
        IF (PROKG) WRITE(MPG,'(A)')
     &  '     Statistics on operation counts (OPC):'
C       To avoid division by zero        
        TOTAL_FLOP = max(TOTAL_FLOP,epsilon(1.0D0))
        LR_TAB(55-LR_DKEEPSHIFT)=real(TOTAL_FLOP)
        LR_TAB(60-LR_DKEEPSHIFT)=100.0D0
        LR_TAB(56-LR_DKEEPSHIFT)=FLOP_FACTO_LR+FLOP_FRFRONTS
        LR_TAB(61-LR_DKEEPSHIFT)=100.0D0*(FLOP_FACTO_LR+FLOP_FRFRONTS)/
     &                                    TOTAL_FLOP
        IF (PROKG) THEN 
        WRITE(MPG,'(A,ES10.3,A,F5.1,A)') 
     &  '     RINFOG(3) Total theoretical operations counts       ='
     &     ,TOTAL_FLOP,' (',100.0D0*TOTAL_FLOP/TOTAL_FLOP,'%)'
        WRITE(MPG,'(A,ES10.3,A,F5.1,A)') 
     &  '     RINFOG(14) Total effective OPC     (% of RINFOG(3)) ='
     &     ,FLOP_FACTO_LR+FLOP_FRFRONTS,' ('
     &,100.0D0*(FLOP_FACTO_LR+FLOP_FRFRONTS)/TOTAL_FLOP
     &,'%)'
        ENDIF
      IF (PROKG) WRITE(MPG,'(A,A)') 
     & '-------------- End of BLR statistics -------------------------',
     & '--------------'
      RETURN
      END SUBROUTINE SAVEandWRITE_GAINS
      END MODULE MUMPS_LR_STATS
