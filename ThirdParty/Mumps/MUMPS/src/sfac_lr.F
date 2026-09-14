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
      MODULE SMUMPS_FAC_LR
      USE SMUMPS_LR_TYPE
      USE SMUMPS_LR_CORE
      USE MUMPS_LR_COMMON, only: NPREC_STOC_MAX
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE SMUMPS_BLR_UPDATE_TRAILING_LDLT(
     &        A, LA, POSELT, IFLAG, IERROR, NFRONT,
     &        BEGS_BLR, NB_BLR, CURRENT_BLR, BLR_L,
     &        NELIM, IW2, BLOCK,
     &        MAXI_CLUSTER, NPIV, NIV, 
     &        MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT)
!$    USE OMP_LIB
C      -----------
C     Parameters
C     -----------
      INTEGER(8), intent(in)       :: LA
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_BLR, CURRENT_BLR,
     &   NELIM, MAXI_CLUSTER, NPIV, NIV, TOL_OPT
      INTEGER, intent(inout)         :: IFLAG, IERROR
      REAL, intent(inout)    :: A(LA)
C     TYPE(LRB_TYPE),intent(in) :: BLR_L(NB_BLR-CURRENT_BLR)
      TYPE(LRB_TYPE),intent(in) :: BLR_L(:)
C     BLOCK is of size (MAXI_CLUSTER, (OMP_NUM+1)*MAXI_CLUSTER)
C     and contiguos in memory:
      REAL, INTENT(INOUT) :: BLOCK(MAXI_CLUSTER,*)
      INTEGER, intent(in) :: IW2(*)
      INTEGER, DIMENSION(:) :: BEGS_BLR
C     MIDBLK_COMPRESS = whether to compress the middle product
C     Y1*X2 in A = A + X1*(Y1*X2)*Y2 (strategy suggested by CA)
      INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT
C     Threshold used for RRQR if MIDBLK_COMPRESS = TRUE 
      REAL,intent(in) :: TOLEPS
C     -----------
C     Local variables
C     -----------
      INTEGER :: I, NB_BLOCKS_PANEL, J, MID_RANK
      LOGICAL :: BUILDQ
      INTEGER :: OMP_NUM
      INTEGER :: IBIS
#if ! defined(BLR_NOOPENMP)
      INTEGER :: CHUNK
#endif
      INTEGER(8) :: POSELTT, POSELTD
      REAL :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
      PARAMETER (ZERO=0.0E0)
      NB_BLOCKS_PANEL = NB_BLR-CURRENT_BLR
C Position of first diagonal value for applying scaling D from LDLT
      POSELTD = POSELT + int(NFRONT,8) * int(BEGS_BLR(CURRENT_BLR)-1,8)
     &          + int(BEGS_BLR(CURRENT_BLR) - 1,8)
      OMP_NUM = 0 ! shared
C Now update the rest of the trailing submatrix using panel BLR blocks     
#if ! defined(BLR_NOOPENMP)
C Merge the two loops into a single loop for performance
      CHUNK = 1
!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
!$OMP& PRIVATE(I, J, POSELTT, OMP_NUM,
!$OMP&         MID_RANK, BUILDQ)
#endif
C     traverse panel
      DO IBIS = 1, (NB_BLOCKS_PANEL*(NB_BLOCKS_PANEL+1)/2) 
        IF (IFLAG.LT.0) CYCLE
C       Determining I,J from IBIS
C       Notice that IBIS = I(I-1)/2 + J - 1      
        I = CEILING((1.0D0+SQRT(1.0D0+8.0D0*dble(IBIS)))/2.0D0)-1
        J = IBIS - I*(I-1)/2
C#else
C      DO I = 1, NB_BLOCKS_PANEL !traverse panel 
C         DO J = 1, I
C#endif
#if ! defined(BLR_NOOPENMP)         
        OMP_NUM = 0
!$      OMP_NUM = OMP_GET_THREAD_NUM() 
#endif
C Position of first value of the block to be updated with a GEMM in the
C trailing submatrix
            POSELTT = POSELT + int(NFRONT,8) *
     &                int(BEGS_BLR(CURRENT_BLR+I)-1,8)
     &           + int(BEGS_BLR(CURRENT_BLR+J) - 1, 8)
C           Use first column of each MAXI_CLUSTER*MAXI_CLUSTER block
            CALL SMUMPS_LRGEMM4(MONE,
     &            BLR_L(J), BLR_L(I), ONE, A, LA, 
     &            POSELTT, NFRONT, 1, IFLAG, IERROR, 
     &            MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT,
     &            MID_RANK, BUILDQ,
     &            .FALSE., MAXI_CLUSTER=MAXI_CLUSTER,
     &            DIAG=A(POSELTD), LD_DIAG=NFRONT, IW2=IW2,
     &            BLOCK=BLOCK(1:MAXI_CLUSTER,OMP_NUM*MAXI_CLUSTER+1))
            IF (IFLAG.LT.0) CYCLE
            CALL UPD_FLOP_UPDATE(BLR_L(J)%M, BLR_L(J)%N, BLR_L(J)%K, 
     &          BLR_L(J)%ISLR, BLR_L(I)%M, BLR_L(I)%N, BLR_L(I)%K,
     &          BLR_L(I)%ISLR, MIDBLK_COMPRESS, MID_RANK, BUILDQ,
     &          (I.EQ.J), .FALSE.)
         ENDDO
#if ! defined(BLR_NOOPENMP) 
!$OMP END DO
#endif
      END SUBROUTINE SMUMPS_BLR_UPDATE_TRAILING_LDLT
      SUBROUTINE SMUMPS_BLR_SLV_UPD_TRAIL_LDLT(A, LA, POSELT, 
     &        IFLAG, IERROR, NCOL, NROW, JBEG_BLOCK,
     &        A_BLOCFACTO, LA_BLOCFACTO, LD_BLOCFACTO, 
     &        BEGS_BLR_LM, NB_BLR_LM, BLR_LM, ISHIFT_LM,
     &        BEGS_BLR_LS, NB_BLR_LS, BLR_LS, ISHIFT_LS,
     &        CURRENT_BLR_LM, CURRENT_BLR_LS,
     &        IW2, BLOCK,
     &        MAXI_CLUSTER,
     &        MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT
     &        )
!$    USE OMP_LIB      
C     
C     Purpose
C     -------
C     Update of the trailing submatrix in case of the band of a slave
C     of size NCOL * NROW.
C     JBEG_BLOCK.EQ.1 indicates that it is the first right looking 
C     update for a given panel factorization. 
C     Step 2 LS D LS^T for last diagonal block is performed 
C     only when JBEG_BLOCK.EQ.1
C     
C     First block to be updated is at position (1, 1+ISHIFT_LM+NPIV+NELIM)
C     Diagonal block of the band of L starts at position NCOL-ROW
C     Diagonal block sent by master is stored in A_BLOCFACTO
C     array of size (NPIV,NPIV+NELIM), passed in LA_BLOCFACTO,
C     with leading dimension LD_BLOCFACTO
C     BLR_LM corresponds to the block LM (NPIV columns) 
C       sent by the master (in LR form)
C       and is associated to the L block of fully summed variables.
C     BLR_LS corresponds to the block LS of NPIV columns of the band 
C       of the slave that have just been factored/updated and compressed
C       by slave
C     NCOL number of column of the band and leading dimension since
C       band is stored by row.
C     NELIM holds the number of variables non eliminated so far (in this front)
C     NPIV holds the number of pivots eliminated in the current BLR block
C       (therefore NPIV+NELIM corresponds to the size of current LR block
C       that the master has tried to factor).
C     ISHIFT_LM, ISHIFT_LS and used to find respectively the position in  
C     the full-rank band stored in A(POSELT) of the beginning of the
C     block associated to LM and to LS blocks
C
C     ALGO:
C      1/ LS D LM^T
C      2/ LS D LS^T
C
C 
C     -----------
C     Parameters
C     -----------
      INTEGER(8), intent(in)  :: LA, LA_BLOCFACTO
      REAL, intent(inout)  :: A(LA)
      REAL, intent(in)  :: A_BLOCFACTO(LA_BLOCFACTO)
      INTEGER(8), intent(in)  :: POSELT 
      INTEGER, intent(inout)    :: IFLAG, IERROR
      INTEGER, intent(in)     :: NCOL, NROW, IW2(*), TOL_OPT,
     &                           MAXI_CLUSTER, LD_BLOCFACTO,
     &                           JBEG_BLOCK
      INTEGER, intent(in)     :: NB_BLR_LM, NB_BLR_LS, 
     &                           ISHIFT_LM, ISHIFT_LS, 
     &                           CURRENT_BLR_LM, CURRENT_BLR_LS
C     BLOCK is of size MAXI_CLUSTER, NOMP*MAXICLUSTER, although
C     only one column of BLOCK is used per OMP in SMUMPS_LRGEMM2
      REAL, INTENT(INOUT) :: BLOCK(MAXI_CLUSTER,*)
      INTEGER, DIMENSION(:) :: BEGS_BLR_LM, BEGS_BLR_LS
      TYPE(LRB_TYPE),intent(in) :: BLR_LM(NB_BLR_LM-CURRENT_BLR_LM),
     &                             BLR_LS(NB_BLR_LS-CURRENT_BLR_LS)
C     MIDBLK_COMPRESS = whether to compress the middle product
C     Y1*X2 in A = A + X1*(Y1*X2)*Y2 (strategy suggested by CA)
      INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT
C     Threshold used for RRQR if MIDBLK_COMPRESS = TRUE 
      REAL,intent(in) :: TOLEPS
C     -----------
C     Local variables
C     -----------
      INTEGER :: I, NB_BLOCKS_PANEL_LM, NB_BLOCKS_PANEL_LS, J, MID_RANK
      LOGICAL :: BUILDQ
      INTEGER :: OMP_NUM
      INTEGER :: IBIS
#if ! defined(BLR_NOOPENMP)
      INTEGER :: CHUNK
#endif
      INTEGER(8) :: POSELTT
      REAL :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
      PARAMETER (ZERO=0.0E0)
      NB_BLOCKS_PANEL_LM = NB_BLR_LM-CURRENT_BLR_LM
      NB_BLOCKS_PANEL_LS = NB_BLR_LS-CURRENT_BLR_LS
C Now update the rest of the trailing submatrix using panel BLR blocks     
C
C OMP_NUM initialization in case of compilation
C with BLR_NOOPENMP and without OpenMP
      OMP_NUM = 0
#if ! defined(BLR_NOOPENMP)
      CHUNK = 1
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
!$OMP& PRIVATE(I, J, POSELTT, OMP_NUM, MID_RANK, BUILDQ) 
#endif
C     traverse panel
      DO IBIS = 1, (NB_BLOCKS_PANEL_LS*NB_BLOCKS_PANEL_LM) 
        IF (IFLAG.LT.0) CYCLE
C       Determining I,J from IBIS
        I = (IBIS-1)/NB_BLOCKS_PANEL_LM+1
        J = IBIS - (I-1)*NB_BLOCKS_PANEL_LM
C      DO I = 1, NB_BLOCKS_PANEL_LS !traverse panels of L from slave
CC        1/ First step : perform  LS D LM^T 
C         DO J = 1, NB_BLOCKS_PANEL_LM ! traverse panels of L from master
#if ! defined(BLR_NOOPENMP)
        OMP_NUM = 0 ! in case we enter this block w/o openmp (private)
!$      OMP_NUM = OMP_GET_THREAD_NUM() 
#endif
C Position of first value of the block to be updated with a GEMM in the
C trailing submatrix
            POSELTT = POSELT 
     &           + int(NCOL,8) * 
     &             int((BEGS_BLR_LS(CURRENT_BLR_LS+I)+ISHIFT_LS-1),8)
     &           + int((BEGS_BLR_LM(CURRENT_BLR_LM+J)+ISHIFT_LM-1),8)
C            write(6,*) "I=", I, "J=",J, 
C     &         "POSELTT-POSELT=", POSELTT - POSELT
            CALL SMUMPS_LRGEMM4(MONE,
     &            BLR_LM(J), BLR_LS(I), ONE, A, LA, 
     &            POSELTT, NCOL, ! beg of bloc to update and LDA
     &            1, IFLAG, IERROR, 
     &            MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT,
     &            MID_RANK, BUILDQ, 
     &            .FALSE., MAXI_CLUSTER=MAXI_CLUSTER,
     &            DIAG=A_BLOCFACTO, LD_DIAG=LD_BLOCFACTO, IW2=IW2,
     &            BLOCK=BLOCK(1:MAXI_CLUSTER,OMP_NUM*MAXI_CLUSTER+1))
            IF (IFLAG.LT.0) CYCLE
            CALL UPD_FLOP_UPDATE(BLR_LM(J)%M, BLR_LM(J)%N, BLR_LM(J)%K, 
     &          BLR_LM(J)%ISLR, BLR_LS(I)%M, BLR_LS(I)%N, BLR_LS(I)%K,
     &          BLR_LS(I)%ISLR, MIDBLK_COMPRESS, MID_RANK, BUILDQ,
     &           .FALSE., .FALSE.)
         ENDDO
#if ! defined(BLR_NOOPENMP)
!$OMP END DO
         IF (IFLAG.LT.0) RETURN
C        Step 2/ only performed when processing first block
         IF (JBEG_BLOCK.NE.1) RETURN
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
!$OMP& PRIVATE(I, J, POSELTT, MID_RANK, OMP_NUM, BUILDQ) 
#endif
C        2/ LS D LS^T for last diagonal block of size (NROW columns)
C           of the band
         DO IBIS = 1, (NB_BLOCKS_PANEL_LS*(NB_BLOCKS_PANEL_LS+1)/2)
          IF (IFLAG.LT.0) CYCLE
C         Determining I,J from IBIS
C         Notice that IBIS = I(I-1)/2 + J - 1      
          I = CEILING((1.0D0+SQRT(1.0D0+8.0D0*dble(IBIS)))/2.0D0)-1
          J = IBIS - I*(I-1)/2
C        DO J = 1, I
#if ! defined(BLR_NOOPENMP)
          OMP_NUM = 0
!$        OMP_NUM = OMP_GET_THREAD_NUM() 
#endif
          POSELTT = POSELT 
     &        + int(NCOL,8) * 
     &          int((BEGS_BLR_LS(CURRENT_BLR_LS+I)+ISHIFT_LS-1),8)
     &        + int((NCOL-NROW+(BEGS_BLR_LS(CURRENT_BLR_LS+J)-1)),8)
C          write(6,*) "DIAG I=", I, 
C     &         "POSELTT-POSELT=", POSELTT - POSELT
          CALL SMUMPS_LRGEMM4(MONE,
     &            BLR_LS(J),BLR_LS(I), ONE, A, LA, 
     &            POSELTT, NCOL, ! beg of bloc to update and LDA
     &            1, IFLAG, IERROR, 
     &            MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT,
     &            MID_RANK, BUILDQ,
     &            .FALSE., MAXI_CLUSTER=MAXI_CLUSTER,
     &            DIAG=A_BLOCFACTO, LD_DIAG=LD_BLOCFACTO, IW2=IW2,
     &            BLOCK=BLOCK(1:MAXI_CLUSTER,OMP_NUM*MAXI_CLUSTER+1))
            IF (IFLAG.LT.0) CYCLE
            CALL UPD_FLOP_UPDATE(BLR_LS(J)%M, BLR_LS(J)%N, BLR_LS(J)%K, 
     &          BLR_LS(J)%ISLR, BLR_LS(I)%M, BLR_LS(I)%N, BLR_LS(I)%K,
     &          BLR_LS(I)%ISLR, MIDBLK_COMPRESS, MID_RANK, BUILDQ, 
     &            (I.EQ.J), .FALSE.)
C         ENDDO
      ENDDO
#if ! defined(BLR_NOOPENMP)
!$OMP END DO
#endif
      RETURN
      END SUBROUTINE SMUMPS_BLR_SLV_UPD_TRAIL_LDLT
      SUBROUTINE SMUMPS_BLR_UPD_NELIM_VAR_U(
     &        A, LA, POSELT, IFLAG, IERROR, NFRONT,
     &        BEGS_BLR, CURRENT_BLR, BLR_U, NB_BLR, 
     &        FIRST_BLOCK, IBEG_BLR, NPIV, NELIM)
!$    USE OMP_LIB
C     
C     Update of the NELIM variables in the U factors of current BLR
C     panel, in case of restricted pivoting. 
C     FIXME We should parallelize the loop in this routine.     
C     -----------
C     Parameters
C     -----------
      INTEGER(8), intent(in)       :: LA
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_BLR, CURRENT_BLR,
     &                                IBEG_BLR, NPIV, NELIM, FIRST_BLOCK
      INTEGER, intent(inout)         :: IFLAG, IERROR
      REAL, TARGET, intent(inout) :: A(LA)
C     TYPE(LRB_TYPE),TARGET,intent(in) :: BLR_U(NB_BLR-CURRENT_BLR)
      TYPE(LRB_TYPE),TARGET,intent(in) :: BLR_U(:)
      INTEGER, DIMENSION(:) :: BEGS_BLR
C     -----------
C     Local variables
C     -----------
      TYPE(LRB_TYPE), POINTER :: LRB
      INTEGER :: IP
      INTEGER :: allocok
      INTEGER(8) :: LPOS, UPOS
      REAL, ALLOCATABLE,DIMENSION(:,:) :: TEMP_BLOCK
      REAL :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
      PARAMETER (ZERO=0.0E0)
      IF (NELIM.NE.0) THEN 
C       Rows corresponding to non eliminated pivots have not been
C       updated in fac_sq (outside of current blr panel) so we do the
C       updates using BLR_U
        LPOS = POSELT + int(NFRONT,8)*int(NPIV,8) + int(IBEG_BLR-1,8)
#if ! defined(BLR_NOOPENMP)
!$OMP DO PRIVATE(LRB, UPOS)
#endif
        DO IP = FIRST_BLOCK, NB_BLR 
          IF (IFLAG.LT.0) CYCLE
          LRB => BLR_U(IP-CURRENT_BLR)
          UPOS = POSELT + int(NFRONT,8)*int(NPIV,8) 
     &                  + int(BEGS_BLR(IP)-1,8)
          IF (LRB%ISLR) THEN
             IF (LRB%K.GT.0) THEN
               allocate(TEMP_BLOCK( LRB%K, NELIM ), stat=allocok )
               IF (allocok .GT. 0) THEN
                 IFLAG  = -13
                 IERROR = NELIM * LRB%K
                 GOTO 100
               ENDIF
               CALL sgemm('N', 'N', LRB%K, NELIM, LRB%N, ONE,
     &              LRB%R(1,1), LRB%K, A(LPOS), NFRONT,
     &              ZERO, TEMP_BLOCK, LRB%K) 
               CALL sgemm('N', 'N', LRB%M, NELIM, LRB%K, MONE,
     &              LRB%Q(1,1), LRB%M, TEMP_BLOCK, LRB%K,
     &              ONE, A(UPOS), NFRONT) 
               deallocate(TEMP_BLOCK)
             ENDIF
          ELSE
            CALL sgemm('N', 'N', LRB%M, NELIM, LRB%N, MONE,
     &              LRB%Q(1,1), LRB%M, A(LPOS), NFRONT,
     &              ONE, A(UPOS), NFRONT) 
          ENDIF
 100    CONTINUE
        ENDDO
#if ! defined(BLR_NOOPENMP)
!$OMP ENDDO
#endif
      ENDIF
      END SUBROUTINE SMUMPS_BLR_UPD_NELIM_VAR_U
      SUBROUTINE SMUMPS_BLR_UPD_NELIM_VAR_L(
     &        A_U, LA_U, UPOS, A_L, LA_L, LPOS, IFLAG, IERROR, LDU, LDL,
     &        BEGS_BLR_L, CURRENT_BLR, BLR_L, NB_BLR_L, 
     &        FIRST_BLOCK, NELIM, UTRANS)
!$    USE OMP_LIB
C     
C     Update of the NELIM variables in the L factors of current BLR
C     panel. 
C     FIXME We should call this routine for all instances of
C           SMUMPS_BLR_UPDATE_TRAILING.      
C     -----------
C     Parameters
C     -----------
      INTEGER(8), intent(in)       :: LA_U, LA_L
      INTEGER(8), intent(in)       :: UPOS, LPOS
      INTEGER, intent(in)          :: LDU, LDL, NB_BLR_L, CURRENT_BLR,
     &                                NELIM,  FIRST_BLOCK
      CHARACTER(len=1),INTENT(IN)  :: UTRANS
      INTEGER, intent(inout)         :: IFLAG, IERROR
      REAL, TARGET, intent(inout) :: A_L(LA_L), A_U(LA_U)
C     TYPE(LRB_TYPE),intent(in) :: BLR_L(NB_BLR_L-CURRENT_BLR)
      TYPE(LRB_TYPE),intent(in) :: BLR_L(:)
      INTEGER                          :: BEGS_BLR_L(:)
C     -----------
C     Local variables
C     -----------
      INTEGER :: I, NB_BLOCKS_PANEL_L, KL, ML, NL
      INTEGER :: allocok
      INTEGER(8) :: IPOS
      REAL, ALLOCATABLE,DIMENSION(:,:) :: TEMP_BLOCK
      REAL :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
      PARAMETER (ZERO=0.0E0)
      NB_BLOCKS_PANEL_L = NB_BLR_L-CURRENT_BLR
      IF (NELIM.NE.0) THEN 
C       Columns corresponding to non eliminated pivots have not been
C       updated in fac_sq (outside of current blr panel) so we do the
C       updates using BLR_L 
#if ! defined(BLR_NOOPENMP)
!$OMP DO PRIVATE(KL, ML, NL, IPOS)
#endif
        DO I = FIRST_BLOCK-CURRENT_BLR, NB_BLOCKS_PANEL_L
          IF (IFLAG.LT.0) CYCLE
          KL = BLR_L(I)%K 
          ML = BLR_L(I)%M 
          NL = BLR_L(I)%N 
          IPOS = LPOS + int(LDL,8) * 
     &        int(BEGS_BLR_L(CURRENT_BLR+I)-BEGS_BLR_L(CURRENT_BLR+1),8)
          IF (BLR_L(I)%ISLR) THEN
             IF (KL.GT.0) THEN
               allocate(TEMP_BLOCK( NELIM, KL ), stat=allocok )
               IF (allocok .GT. 0) THEN
                 IFLAG  = -13
                 IERROR = NELIM * KL
                 write(*,*) 'Allocation problem in BLR routine 
     &         SMUMPS_BLR_UPD_NELIM_VAR_L: ',
     &         'not enough memory? memory requested = ', IERROR
                 GOTO 100
               ENDIF
               CALL sgemm(UTRANS , 'T' , NELIM, KL, NL , ONE ,
     &              A_U(UPOS) , LDU , BLR_L(I)%R(1,1) , KL ,
     &              ZERO , TEMP_BLOCK , NELIM) 
               CALL sgemm('N' , 'T' , NELIM , ML , KL , MONE ,
     &              TEMP_BLOCK , NELIM , BLR_L(I)%Q(1,1) , ML ,
     &              ONE , A_L(IPOS) , LDL) 
               deallocate(TEMP_BLOCK)
             ENDIF
          ELSE
            CALL sgemm(UTRANS , 'T' , NELIM, ML, NL , MONE ,
     &              A_U(UPOS) , LDU , BLR_L(I)%Q(1,1) , ML ,
     &              ONE , A_L(IPOS) , LDL) 
          ENDIF
 100      CONTINUE
        ENDDO
#if ! defined(BLR_NOOPENMP)
!$OMP ENDDO
#endif
      ENDIF
      END SUBROUTINE SMUMPS_BLR_UPD_NELIM_VAR_L
      SUBROUTINE SMUMPS_BLR_UPDATE_TRAILING(
     &        A, LA, POSELT, IFLAG, IERROR, NFRONT,
     &        BEGS_BLR_L, BEGS_BLR_U, CURRENT_BLR, BLR_L, NB_BLR_L, 
     &        BLR_U, NB_BLR_U, NELIM, LBANDSLAVE, ISHIFT, NIV, SYM,
     &        MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT)
!$    USE OMP_LIB
C     
C     Update of trailing submatrix of the band of a slave
C     Case of LU factorization and case of LDLT update
C     of a slave at the reception of a message BLFAC_SLAVE
C     sent by a previous slave
C
C          
C
C     -----------
C     Parameters
C     -----------
      INTEGER(8), intent(in)       :: LA
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_BLR_L, NB_BLR_U, 
     &                                CURRENT_BLR,
     &                                NELIM, NIV, SYM, TOL_OPT
      INTEGER, intent(inout)       :: IFLAG, IERROR
      LOGICAL, intent(in)          :: LBANDSLAVE
C  ISHIFT used to shift columns of L band of a slave 
C      only used when  LBANDSLAVE is true and correspond to the first column
C      to be updated
C     it should be stored in BEGS_BLR_U but since the ISHIFT first columns
C     of U are not stored, it needs a special processing
      INTEGER, intent(in)          :: ISHIFT
      REAL, intent(inout) :: A(LA)
C     TYPE(LRB_TYPE),TARGET,intent(in) :: BLR_U(NB_BLR_U-CURRENT_BLR)
C     TYPE(LRB_TYPE),TARGET,intent(in) :: BLR_L(NB_BLR_L-CURRENT_BLR)
      TYPE(LRB_TYPE),TARGET,intent(in) :: BLR_U(:)
      TYPE(LRB_TYPE),TARGET,intent(in) :: BLR_L(:)
      INTEGER :: BEGS_BLR_L(:), BEGS_BLR_U(:)
C     MIDBLK_COMPRESS = whether to compress the middle product
C     Y1*X2 in A = A + X1*(Y1*X2)*Y2 (strategy suggested by CA)
      INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT
C     Threshold used for RRQR if MIDBLK_COMPRESS = TRUE 
      REAL,intent(in) :: TOLEPS
C     -----------
C     Local variables
C     -----------
      INTEGER :: I, NB_BLOCKS_PANEL_L, NB_BLOCKS_PANEL_U, 
     &           KL, ML, NL, J, IS, MID_RANK
      INTEGER :: allocok
      LOGICAL :: BUILDQ
      INTEGER :: IBIS
#if ! defined(BLR_NOOPENMP)
      INTEGER :: CHUNK
#endif
      INTEGER(8) :: POSELT_INCB, POSELT_TOP
      REAL, ALLOCATABLE,DIMENSION(:,:) :: TEMP_BLOCK
      REAL :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
      PARAMETER (ZERO=0.0E0)
      NB_BLOCKS_PANEL_L = NB_BLR_L-CURRENT_BLR
      NB_BLOCKS_PANEL_U = NB_BLR_U-CURRENT_BLR
      IF (LBANDSLAVE) THEN
       IS = ISHIFT
      ELSE
       IS = 0
      ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP SINGLE
#endif
C     NELIM=0 when called from fac_process_blocfacto
      IF (NELIM.NE.0) THEN 
C     Columns corresponding to non eliminated pivots have not been
C     updated in fac_sq (outside of current blr panel) so we do the
C     updates using BLR_L 
         DO I = 1, NB_BLOCKS_PANEL_L
C this should be put in a single routine for this purpose .... 
            KL = BLR_L(I)%K 
            ML = BLR_L(I)%M 
            NL = BLR_L(I)%N 
            IF (BLR_L(I)%ISLR) THEN
               IF (KL.GT.0) THEN
C
CFIXME allocate once of size (NELIM,MAXICLUSTER) ouside loop I
               allocate(TEMP_BLOCK( NELIM, KL ), stat=allocok )
               IF (allocok .GT. 0) THEN
                 IFLAG  = -13
                 IERROR = NELIM * KL
                 GOTO 100
               ENDIF
               POSELT_TOP  = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR_U(CURRENT_BLR)-1),8)
     &           + int(BEGS_BLR_U(CURRENT_BLR+1) + IS - NELIM - 1,8)
               POSELT_INCB = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR_L(CURRENT_BLR+I)-1),8)
     &           + int(BEGS_BLR_U(CURRENT_BLR+1)+IS-NELIM-1,8)
               CALL sgemm('N' , 'T' , NELIM, KL, NL , ONE ,
     &                A(POSELT_TOP) , NFRONT , BLR_L(I)%R(1,1) , KL ,
     &                ZERO , TEMP_BLOCK , NELIM) 
               CALL sgemm('N' , 'T' , NELIM , ML , KL , MONE ,
     &                TEMP_BLOCK , NELIM , BLR_L(I)%Q(1,1) , ML ,
     &                ONE , A(POSELT_INCB) , NFRONT) 
               deallocate(TEMP_BLOCK)
               ENDIF
            ELSE
              POSELT_TOP  = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR_L(CURRENT_BLR)-1),8)
     &           + int(BEGS_BLR_U(CURRENT_BLR+1)+IS-NELIM-1,8)
              POSELT_INCB = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR_L(CURRENT_BLR+I)-1),8)
     &           + int(BEGS_BLR_U(CURRENT_BLR+1) + IS - NELIM - 1, 8)
               CALL sgemm('N' , 'T' , NELIM, ML, NL , MONE ,
     &                A(POSELT_TOP) , NFRONT , BLR_L(I)%Q(1,1) , ML ,
     &                ONE , A(POSELT_INCB) , NFRONT) 
            ENDIF
         ENDDO
      ENDIF
 100  CONTINUE
#if ! defined(BLR_NOOPENMP)
!$OMP END SINGLE
#endif
      IF (IFLAG.LT.0) GOTO 200
C Now update the rest of the trailing submatrix using BLR blocks     
#if ! defined(BLR_NOOPENMP)
C Merge the two loops into a single loop for performance
      CHUNK = 1
!$OMP DO SCHEDULE(DYNAMIC,CHUNK) 
!$OMP& PRIVATE(I, J, POSELT_INCB, MID_RANK, BUILDQ)
#endif
      DO IBIS = 1, (NB_BLOCKS_PANEL_L*NB_BLOCKS_PANEL_U) !traverse panel 
        IF (IFLAG.LT.0) CYCLE
C       Determining I, J from IBIS
C       I: number of block in L/row direction
C       J: number of block in column direction (of row panel I)
C         J: nb with respect to the NB_BLR_U/NBLRB_PACKETS received
        I = (IBIS-1)/NB_BLOCKS_PANEL_U+1
        J = IBIS - (I-1)*NB_BLOCKS_PANEL_U
        POSELT_INCB = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR_L(CURRENT_BLR+I)-1),8)
     &           + int(BEGS_BLR_U(CURRENT_BLR+J) +IS - 1,8)
         CALL SMUMPS_LRGEMM4(MONE, BLR_U(J),
     &           BLR_L(I), ONE, A, LA, POSELT_INCB,
     &           NFRONT, 0, IFLAG, IERROR, 
     &           MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &           KPERCENT, MID_RANK, BUILDQ, .FALSE.)
         IF (IFLAG.LT.0) CYCLE
            CALL UPD_FLOP_UPDATE(BLR_U(J)%M, BLR_U(J)%N, BLR_U(J)%K, 
     &          BLR_U(J)%ISLR, BLR_L(I)%M, BLR_L(I)%N, BLR_L(I)%K,
     &          BLR_L(I)%ISLR, MIDBLK_COMPRESS, MID_RANK, BUILDQ,
     &           .FALSE., .FALSE.)
       ENDDO
#if ! defined(BLR_NOOPENMP)
!$OMP END DO
#endif
 200  CONTINUE
      END SUBROUTINE SMUMPS_BLR_UPDATE_TRAILING
      SUBROUTINE SMUMPS_BLR_UPD_PANEL_LEFT_LDLT(
     &        A, LA, POSELT, NFRONT, IWHANDLER,
     &        BEGS_BLR, CURRENT_BLR, NB_BLR, NPARTSASS,
     &        NELIM, IW2, BLOCK, ACC_LUA,
     &        MAXI_CLUSTER, MAXI_RANK, NIV, IFLAG, IERROR,
     &        MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT_RMB,
     &        K480, K479, K478, KPERCENT_LUA, KPERCENT,
     &        KEEP8,
     &        FIRST_BLOCK
     &        )
!$    USE OMP_LIB
C     
C     Left-looking update of the next BLR panel. 
C     -----------
C     Parameters
C     -----------
      INTEGER(8), intent(in)       :: LA
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_BLR, NPARTSASS,
     &                                CURRENT_BLR, IWHANDLER, TOL_OPT,
     &                                NELIM, NIV, K480, K479, K478,
     &                                MAXI_CLUSTER, MAXI_RANK,
     &                                KPERCENT_LUA, KPERCENT
      REAL, intent(inout)    :: A(LA)
      INTEGER, intent(in) :: IW2(*)
      REAL :: BLOCK(MAXI_CLUSTER,*)
      TYPE(LRB_TYPE), POINTER :: ACC_LUA(:)
      INTEGER(8) :: KEEP8(150)
      INTEGER, DIMENSION(:) :: BEGS_BLR
C     MIDBLK_COMPRESS = whether to compress the middle product
C     Y1*X2 in A = A + X1*(Y1*X2)*Y2 (strategy suggested by CA)
      INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT_RMB
C     Threshold used for RRQR if MIDBLK_COMPRESS = TRUE 
      REAL,intent(in) :: TOLEPS
      INTEGER,intent(inout) :: IFLAG, IERROR
C     When LRTRSM is activated, KEEP(475) >= 2, and KEEP(480) >=5,
C     FIRST_BLOCK indicates the first block that can be compressed 
C     before the trsm: FIRST_BLOCK = NPARTSASS-CURRENT_BLR 
C     if KEEP(475) = 2 and FIRST_BLOCK = 1 if KEEP(475) = 3 
      INTEGER,OPTIONAL,intent(in) :: FIRST_BLOCK
C     -----------
C     Local variables
C     -----------
      TYPE(LRB_TYPE), POINTER :: BLR_L(:), NEXT_BLR_L(:)
      TYPE(LRB_TYPE), POINTER :: ACC_LRB
      INTEGER :: OLD_ACC_RANK, MAX_ACC_RANK, NEW_ACC_RANK, FRFR_UPDATES,
     &           I, II, J, JJ, NB_BLOCKS_PANEL, IND_U, IND_L, K_MAX,
     &           MAXRANK, NB_DEC, FR_RANK
      INTEGER :: MID_RANK, allocok
      INTEGER :: J_ORDER(CURRENT_BLR), J_RANK(CURRENT_BLR)
      INTEGER, ALLOCATABLE :: POS_LIST(:), RANK_LIST(:)
      LOGICAL :: BUILDQ, COMPRESSED_FR
      INTEGER :: OFFSET_IW
      INTEGER :: OMP_NUM
#if ! defined(BLR_NOOPENMP)
      INTEGER :: CHUNK
#endif
      INTEGER(8) :: POSELT_INCB, POSELTD
      REAL :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
      PARAMETER (ZERO=0.0E0)
      NB_BLOCKS_PANEL = NB_BLR-CURRENT_BLR
      ACC_LRB => ACC_LUA(1)
      IF (K480.GE.5) THEN
        IF (NB_BLOCKS_PANEL.GT.1) THEN
          CALL SMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &         IWHANDLER,
     &         0, ! L Panel
     &         CURRENT_BLR+1, NEXT_BLR_L)
        ENDIF
        IF (.not.(present(FIRST_BLOCK))) THEN
          write(*,*) "Internal error in 
     &      SMUMPS_BLR_UPD_PANEL_LEFT_LDLT: KEEP(480)=",K480,
     &      ">= 5, but FIRST_BLOCK argument is missing"
          CALL MUMPS_ABORT()
        ENDIF
      ENDIF
      OMP_NUM = 0
#if ! defined(BLR_NOOPENMP)
      CHUNK = 1
!$OMP DO SCHEDULE(DYNAMIC,CHUNK) 
!$OMP& PRIVATE(I, J, JJ, POSELT_INCB, MID_RANK, BUILDQ, K_MAX,
!$OMP&         BLR_L, OMP_NUM, J_ORDER, J_RANK,
!$OMP&         IND_U, IND_L, ACC_LRB, POSELTD, NB_DEC,
!$OMP&         MAX_ACC_RANK, OLD_ACC_RANK, NEW_ACC_RANK,
!$OMP&         FRFR_UPDATES, COMPRESSED_FR, FR_RANK, II, OFFSET_IW)
#endif
      DO I = 1, NB_BLOCKS_PANEL
#if ! defined(BLR_NOOPENMP)
        IF (IFLAG.LT.0) CYCLE
        OMP_NUM = 0
!$      OMP_NUM = OMP_GET_THREAD_NUM() 
        ACC_LRB => ACC_LUA(OMP_NUM+1)
#endif
        POSELT_INCB = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR(CURRENT_BLR+I)-1),8)
     &           + int(BEGS_BLR(CURRENT_BLR+1)-1,8)
C       Initialize actual size of ACC_LRB
        ACC_LRB%N = BEGS_BLR(CURRENT_BLR+I+1)-BEGS_BLR(CURRENT_BLR+I)
        ACC_LRB%M = BEGS_BLR(CURRENT_BLR+2)-BEGS_BLR(CURRENT_BLR+1)
        MAX_ACC_RANK = 0
        NEW_ACC_RANK = 0
        COMPRESSED_FR = .FALSE.
        IF (K480.EQ.2) THEN
C         Pure left-looking facto, no LUA done 
          DO J = 1, CURRENT_BLR
            J_ORDER(J) = J 
          ENDDO
        ELSE
C         LUA activated 
          CALL SMUMPS_GET_LUA_ORDER(CURRENT_BLR, J_ORDER, J_RANK,
     &                              IWHANDLER, 
     &                              1, 0, I, 0, ! SYM, FS, I, LorU 
     &                              FRFR_UPDATES)
        ENDIF
        FR_RANK = 0
        IF ((K480.GE.5).AND.(I.NE.1)) THEN
          IF (I.GT.FIRST_BLOCK) THEN
C           This block is outside the FS x FS zone, we can
C           compress it early to avoid its decompression.                  
            IF (FRFR_UPDATES.EQ.0) THEN
C             No FR update, we compress A and initialize the
C             accumulator with it
              CALL SMUMPS_COMPRESS_FR_UPDATES(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB,
     &               NFRONT, NIV, TOLEPS, TOL_OPT, KPERCENT, 
     &               COMPRESSED_FR, 0, .FALSE.)
              MAX_ACC_RANK = ACC_LRB%K
              NEW_ACC_RANK = ACC_LRB%K
              FR_RANK = ACC_LRB%K
            ENDIF
          ENDIF
        ENDIF
        NB_DEC = FRFR_UPDATES
        DO JJ = 1, CURRENT_BLR
          J = J_ORDER(JJ)
          K_MAX = J_RANK(JJ)
C         Position of first diagonal value for applying scaling D from LDLT
          POSELTD = POSELT + int(NFRONT,8) * int(BEGS_BLR(J)-1,8)
     &          + int(BEGS_BLR(J) - 1,8)
          OFFSET_IW = BEGS_BLR(J)
          IND_L = CURRENT_BLR+I-J
          IND_U = CURRENT_BLR+1-J
          CALL SMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &           IWHANDLER,
     &           0, ! L Panel
     &           J, BLR_L)
          IF (BLR_L(IND_L)%M.EQ.0) THEN
            CYCLE
          ENDIF
          IF (K480.GE.3) THEN
            IF (ACC_LRB%K+K_MAX.GT.MAXI_RANK) THEN
C             FIXME We should still recompress when we're forced to
C             decompress                
              NB_DEC = JJ-1
              CALL SMUMPS_DECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &              MAXI_RANK, A, LA, POSELT_INCB, NFRONT, NIV, 0)
              COMPRESSED_FR = .FALSE.
              MAX_ACC_RANK = 0
            ENDIF
            OLD_ACC_RANK = ACC_LRB%K
          ENDIF
          CALL SMUMPS_LRGEMM4(MONE,
     &              BLR_L(IND_U), BLR_L(IND_L), ONE,
     &              A, LA, POSELT_INCB,
     &              NFRONT, 1, IFLAG, IERROR, 
     &              MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &              KPERCENT_RMB, MID_RANK, BUILDQ, 
     &              (K480.GE.3), LorU=0, 
     &              LRB3=ACC_LRB, MAXI_RANK=MAXI_RANK,
     &              MAXI_CLUSTER=MAXI_CLUSTER,
     &              DIAG=A(POSELTD), LD_DIAG=NFRONT, 
     &              IW2=IW2(OFFSET_IW),
     &              BLOCK=BLOCK(1:MAXI_CLUSTER,OMP_NUM*MAXI_CLUSTER+1))
          IF (IFLAG.LT.0) GOTO 100
            CALL UPD_FLOP_UPDATE(BLR_L(IND_U)%M, BLR_L(IND_U)%N, 
     &        BLR_L(IND_U)%K, BLR_L(IND_U)%ISLR, BLR_L(IND_L)%M, 
     &        BLR_L(IND_L)%N, BLR_L(IND_L)%K, BLR_L(IND_L)%ISLR,
     &        MIDBLK_COMPRESS, MID_RANK, BUILDQ, (I.EQ.1), (K480.GE.3))
C         Update K_RANK when K481=1 (later used by n-ary tree
C         recompression)
          IF ((MIDBLK_COMPRESS.GE.1).AND.BUILDQ) THEN
            J_RANK(JJ) = MID_RANK
          ENDIF
          IF (K480.GE.3) THEN
            NEW_ACC_RANK = NEW_ACC_RANK + ACC_LRB%K - OLD_ACC_RANK
            MAX_ACC_RANK = MAX(MAX_ACC_RANK, ACC_LRB%K - OLD_ACC_RANK)
            IF (K480.EQ.4) THEN
C             LUA with recompression
              IF ((K478.GT.0).AND.((ACC_LRB%K-MAX_ACC_RANK).GE.K478)) 
     &          THEN
C               Recompress accumulator  
                IF (ACC_LRB%K.GT.0) THEN
                  CALL SMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &                 MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB,
     &                 NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS,
     &                 TOL_OPT,
     &                 KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
                  MAX_ACC_RANK = ACC_LRB%K
                ENDIF
              ENDIF
            ENDIF
            IF ((K480.GE.5).AND.(I.NE.1)) THEN
              IF (I.GT.FIRST_BLOCK) THEN
C               This block is outside the FS x FS zone, we can
C               compress it early to avoid its decompression.                  
                IF (JJ.EQ.FRFR_UPDATES) THEN
C                 Last FR update, we compress FR part and initialize the
C                 accumulator with it
                  CALL SMUMPS_COMPRESS_FR_UPDATES(ACC_LRB,
     &                MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB,
     &                NFRONT, NIV, TOLEPS, TOL_OPT, KPERCENT, 
     &                COMPRESSED_FR, 0, .FALSE.)
                  MAX_ACC_RANK = ACC_LRB%K
                  NEW_ACC_RANK = ACC_LRB%K
                  IF (COMPRESSED_FR) THEN
                    J_RANK(JJ) = ACC_LRB%K 
                    NB_DEC = FRFR_UPDATES-1 
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        IF (K480.GE.3) THEN
          IF ((K480.GE.5)) THEN
C FIXME Implement strategy 5/6 with comb tree recompression            
            IF (COMPRESSED_FR.OR.(K480.GE.6)) THEN  
C             Recompress accumulator
C             Even if NOT COMPRESSED_FR, we can still recompress the
C             accumulator to reduce the decompress flops: strategy
C             K480=6
              IF (ACC_LRB%K.GT.0) THEN
                IF (K478.EQ.-1) THEN
C                 K478=-1: star graph recompression  
                  IF (CURRENT_BLR-FRFR_UPDATES.GT.1) THEN
                    CALL SMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
                  ENDIF
                ELSEIF (K478.LE.-2) THEN
C                 K478=-n, n>=2: n-ary tree recompression  
                  IF (FRFR_UPDATES.GT.0) THEN
C                   We use J_RANK(FRFR_UPDATES) to store the rank of the
C                   compressed FR updates
                    allocate(POS_LIST(CURRENT_BLR-NB_DEC),stat=allocok)
                    IF (allocok .GT. 0) THEN
                       IFLAG  = -13
                       IERROR = CURRENT_BLR-NB_DEC
                       write(*,*) 'Allocation problem in BLR routine ',
     &                      'SMUMPS_BLR_UPD_PANEL_LEFT_LDLT: ',
     &                      'not enough memory? memory requested = ',
     &                      IERROR
                       GOTO 100
                    ENDIF
                    POS_LIST(1) = 1
                    DO II = 1,CURRENT_BLR-NB_DEC-1
                      POS_LIST(II+1)=POS_LIST(II)+J_RANK(NB_DEC+II) 
                    ENDDO
                    CALL SMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB, KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               J_RANK(NB_DEC+1:CURRENT_BLR), POS_LIST,
     &               CURRENT_BLR-NB_DEC, 0)
                  ELSE
C                   We allocate a new RANK_LIST array to store FR_RANK
C                   (rank of the compressed FR updates) followed by
C                   J_RANK
                     allocate(POS_LIST(CURRENT_BLR+1),stat=allocok)
                     IF (allocok .GT. 0) THEN
                        IFLAG  = -13
                        IERROR = CURRENT_BLR+1
                        write(*,*) 'Allocation problem in BLR routine ',
     &                       'SMUMPS_BLR_UPD_PANEL_LEFT_LDLT: ',
     &                       'not enough memory? memory requested = ',
     &                       IERROR
                        GOTO 100
                     ENDIF
                    POS_LIST(1) = 1
                    POS_LIST(2) = 1 + FR_RANK
                    DO II = 2,CURRENT_BLR
                      POS_LIST(II+1)=POS_LIST(II)+J_RANK(II-1) 
                    ENDDO
                    allocate(RANK_LIST(CURRENT_BLR+1),stat=allocok)
                    IF (allocok .GT. 0) THEN
                       IFLAG  = -13
                       IERROR = CURRENT_BLR+1
                       write(*,*) 'Allocation problem in BLR routine ',
     &                      'SMUMPS_BLR_UPD_PANEL_LEFT_LDLT: ',
     &                      'not enough memory? memory requested = ',
     &                      IERROR
                       GOTO 100
                    ENDIF
                    RANK_LIST(1) = FR_RANK
                    DO II = 2,CURRENT_BLR+1
                      RANK_LIST(II) = J_RANK(II-1)
                    ENDDO
                    CALL SMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB, KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               RANK_LIST, POS_LIST,
     &               CURRENT_BLR+1, 0)
                    deallocate(RANK_LIST)
                  ENDIF
                  deallocate(POS_LIST)
                ENDIF
              ENDIF
            ENDIF
            MAXRANK = floor(real(ACC_LRB%M*ACC_LRB%N)/real(ACC_LRB%M+
     &                                                     ACC_LRB%N))
            IF (COMPRESSED_FR.AND.(ACC_LRB%K.LE.MAXRANK)) THEN  
C             It's better to keep this block in LR form
              CALL ALLOC_LRB_FROM_ACC(ACC_LRB, NEXT_BLR_L(I-1), 
     &                       ACC_LRB%K, ACC_LRB%M, ACC_LRB%N, 0,
     &                       IFLAG, IERROR, KEEP8)
              IF (IFLAG.LT.0) CYCLE
              ACC_LRB%K = 0
            ELSE
C FIXME       Differentiate the case where we did not succeed in
C             compressing the FR part (in which case we can decompress
C             the LR part and try compressing the sum FR+LR) and the
C             case where the rank of the accumulator exceeds MAXRANK (in
C             which case there's no need to try compressing again)
C             It's better to keep this block in FR form
              IF (I.NE.1) NEXT_BLR_L(I-1)%ISLR=.FALSE.
              CALL SMUMPS_DECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &              MAXI_RANK, A, LA, POSELT_INCB, NFRONT, NIV, 0)
            ENDIF
          ELSE ! KEEP(480) < 5
C           Finished updating block, in K480=3,4 strategies we
C           decompress the accumulator.
            IF ((K480.EQ.4).AND.(K478.EQ.-1).AND.(ACC_LRB%K.GT.0)) THEN
C             K478=-1: star graph recompression  
              IF (CURRENT_BLR-FRFR_UPDATES.GT.1) THEN
                CALL SMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
              ENDIF
            ELSEIF ((K480.EQ.4).AND.(K478.LE.-2).AND.(ACC_LRB%K.GT.0)) 
     &        THEN
C             K478=-n, n>=2: n-ary tree recompression  
               allocate(POS_LIST(CURRENT_BLR-NB_DEC),stat=allocok)
               IF (allocok .GT. 0) THEN
                  IFLAG  = -13
                  IERROR = CURRENT_BLR-NB_DEC
                  GOTO 100
               ENDIF
              POS_LIST(1) = 1
              DO II = 1,CURRENT_BLR-NB_DEC-1
                POS_LIST(II+1)=POS_LIST(II)+J_RANK(NB_DEC+II) 
              ENDDO
              CALL SMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB, KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               J_RANK(NB_DEC+1:CURRENT_BLR), POS_LIST,
     &               CURRENT_BLR-NB_DEC, 0)
              deallocate(POS_LIST)
            ENDIF
            CALL SMUMPS_DECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &            MAXI_RANK, A, LA, POSELT_INCB, NFRONT, NIV, 0)
          ENDIF                               
        ENDIF
 100    CONTINUE        
      ENDDO                                 
#if ! defined(BLR_NOOPENMP)
!$OMP END DO                                
#endif
      END SUBROUTINE SMUMPS_BLR_UPD_PANEL_LEFT_LDLT
      SUBROUTINE SMUMPS_BLR_UPD_PANEL_LEFT(
     &        A, LA, POSELT, NFRONT, IWHANDLER, LorU,
     &        BEGS_BLR, BEGS_BLR_U, CURRENT_BLR, ACC_LUA,
     &        NB_BLR, NPARTSASS, NELIM, NIV, SYM, 
     &        LBANDSLAVE, IFLAG, IERROR, ISHIFT,
     &        MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT_RMB,
     &        K480, K479, K478, KPERCENT_LUA, KPERCENT,
     &        MAXI_CLUSTER, MAXI_RANK, 
     &        K474, FSorCB, BLR_U_COL, KEEP8,
     &        FIRST_BLOCK, BEG_I_IN, END_I_IN)
!$    USE OMP_LIB
C     
C     Left-looking update of the next BLR panel. 
C     Rather than duplicate this routine's code twice to update the L
C     part and the U part, we only update the L part when LorU = 0 and
C     the U part when LorU = 1 or other.  Thus this routine is called
C     twice by fac_front_LU_type1.  
C     -----------
C     Parameters
C     -----------
      INTEGER(8), intent(in)       :: LA
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_BLR, NPARTSASS,
     &                                CURRENT_BLR, IWHANDLER, LorU,
     &                                NELIM, NIV, SYM, K480, K479, K478,
     &                                MAXI_CLUSTER, MAXI_RANK,
     &                                KPERCENT_LUA, KPERCENT, ISHIFT,
     &                                K474, FSorCB
      LOGICAL, intent(in)          :: LBANDSLAVE
      REAL, TARGET, intent(inout) :: A(LA)
      TYPE(LRB_TYPE), POINTER :: ACC_LUA(:), BLR_U_COL(:)
      INTEGER(8) :: KEEP8(150)
      INTEGER, DIMENSION(:) :: BEGS_BLR, BEGS_BLR_U
C     MIDBLK_COMPRESS = whether to compress the middle product
C     Y1*X2 in A = A + X1*(Y1*X2)*Y2 (strategy suggested by CA)
      INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT_RMB, TOL_OPT
C     Threshold used for RRQR if MIDBLK_COMPRESS = TRUE 
      REAL,intent(in) :: TOLEPS
      INTEGER,intent(inout) :: IFLAG, IERROR
C     When LRTRSM is activated, KEEP(475) >= 2, and KEEP(480) >=5,
C     FIRST_BLOCK indicates the first block that can be compressed 
C     before the trsm: FIRST_BLOCK = NPARTSASS-CURRENT_BLR
C     if KEEP(475) = 2 and FIRST_BLOCK = 1 if KEEP(475) = 3 
      INTEGER,OPTIONAL,intent(in) :: FIRST_BLOCK
      INTEGER,OPTIONAL,intent(in) :: BEG_I_IN, END_I_IN
C     -----------
C     Local variables
C     -----------
      TYPE(LRB_TYPE), POINTER :: BLR_U(:), BLR_L(:), NEXT_BLR(:)
      TYPE(LRB_TYPE), POINTER :: ACC_LRB
      INTEGER :: OLD_ACC_RANK, MAX_ACC_RANK, NEW_ACC_RANK, FRFR_UPDATES,
     &           NB_DEC, FR_RANK, MAXRANK, BEG_I, END_I
      INTEGER :: I,II,J,JJ, NB_BLOCKS_PANEL, IND_U, IND_L, K_MAX
      INTEGER :: MID_RANK, allocok
      INTEGER :: J_ORDER(CURRENT_BLR), J_RANK(CURRENT_BLR)
      INTEGER, ALLOCATABLE :: POS_LIST(:), RANK_LIST(:)
      LOGICAL :: BUILDQ, COMPRESSED_FR
#if ! defined(BLR_NOOPENMP)
      INTEGER :: OMP_NUM
      INTEGER :: CHUNK
#endif
      INTEGER(8) :: POSELT_INCB
      REAL :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
      PARAMETER (ZERO=0.0E0)
      IF (NIV.EQ.2.AND.LorU.EQ.0) THEN
        IF (LBANDSLAVE) THEN
          NB_BLOCKS_PANEL = NB_BLR
        ELSE
          NB_BLOCKS_PANEL = NPARTSASS-CURRENT_BLR
        ENDIF
      ELSE
        NB_BLOCKS_PANEL = NB_BLR-CURRENT_BLR
      ENDIF
      ACC_LRB => ACC_LUA(1)
      IF (K480.GE.5) THEN
        IF (NB_BLOCKS_PANEL.GT.1) THEN
          CALL SMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &         IWHANDLER,
     &         LorU, ! L Panel
     &         CURRENT_BLR+1, NEXT_BLR)
        ENDIF
        IF (.not.(present(FIRST_BLOCK))) THEN
          write(*,*) "Internal error in 
     &      SMUMPS_BLR_UPD_PANEL_LEFT: KEEP(480)=",K480,
     &      ">=5, but FIRST_BLOCK argument is missing"
          CALL MUMPS_ABORT()
        ENDIF
      ENDIF
C To avoid updating the diagonal block twice in LU,
C we only do the update in the L panel case
      IF (LorU.EQ.0) THEN ! L panel update
          BEG_I = 1
      ELSE ! U panel update
          BEG_I = 2
      ENDIF
      END_I = NB_BLOCKS_PANEL
      IF (K474.EQ.3) THEN
        IF(present(BEG_I_IN)) THEN
          BEG_I = BEG_I_IN - CURRENT_BLR
        ENDIF
        IF(present(END_I_IN)) THEN
          END_I = END_I_IN - CURRENT_BLR
        ENDIF
      ENDIF
#if ! defined(BLR_NOOPENMP)
      CHUNK = 1
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
!$OMP& PRIVATE(I, J, JJ, POSELT_INCB, MID_RANK, BUILDQ,
!$OMP&         BLR_U, BLR_L, J_ORDER, J_RANK, K_MAX,
!$OMP&         IND_U, IND_L, OMP_NUM, ACC_LRB,
!$OMP&         MAX_ACC_RANK, OLD_ACC_RANK, NEW_ACC_RANK,
!$OMP&         FRFR_UPDATES, FR_RANK, COMPRESSED_FR)
#endif
      DO I = BEG_I, END_I
        IF (IFLAG.LT.0) CYCLE
#if ! defined(BLR_NOOPENMP)
        OMP_NUM = 0
!$      OMP_NUM = OMP_GET_THREAD_NUM() 
        ACC_LRB => ACC_LUA(OMP_NUM+1)
#endif
        IF (LorU.EQ.0) THEN ! L panel update
          IF (LBANDSLAVE) THEN
            POSELT_INCB = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR(I+1)-1),8)
     &           + int(BEGS_BLR_U(2)+ISHIFT-1,8)
C           Initialize actual size of ACC_LRB
            ACC_LRB%N = BEGS_BLR(I+2)-BEGS_BLR(I+1)
            ACC_LRB%M = BEGS_BLR_U(3)-BEGS_BLR_U(2)
            IF (K474.GE.2) THEN
C             LUA-NIV2 V2 strategy
              BLR_U => BLR_U_COL
            ENDIF
          ELSE
            POSELT_INCB = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR(CURRENT_BLR+I)-1),8)
     &           + int(BEGS_BLR(CURRENT_BLR+1)-1,8)
C           Initialize actual size of ACC_LRB
            ACC_LRB%N = BEGS_BLR(CURRENT_BLR+I+1)
     &                 -BEGS_BLR(CURRENT_BLR+I)
            ACC_LRB%M = BEGS_BLR(CURRENT_BLR+2)-BEGS_BLR(CURRENT_BLR+1)
          ENDIF
        ELSE ! U panel update
          POSELT_INCB = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR(CURRENT_BLR+1)-1),8)
     &           + int(BEGS_BLR(CURRENT_BLR+I)-1,8)
C         Initialize actual size of ACC_LRB
          ACC_LRB%N = BEGS_BLR(CURRENT_BLR+2)-BEGS_BLR(CURRENT_BLR+1)
          ACC_LRB%M = BEGS_BLR(CURRENT_BLR+I+1)-BEGS_BLR(CURRENT_BLR+I)
        ENDIF
        MAX_ACC_RANK = 0
        NEW_ACC_RANK = 0
        COMPRESSED_FR = .FALSE.
        IF (K480.EQ.2) THEN
C         Pure left-looking facto, no LUA done 
          DO J = 1, CURRENT_BLR
            J_ORDER(J) = J 
          ENDDO
        ELSE
C         LUA activated 
          CALL SMUMPS_GET_LUA_ORDER(CURRENT_BLR, J_ORDER, J_RANK, 
     &                              IWHANDLER, 
     &                              0, 0, I, LorU, ! SYM, FS, I, LorU 
     &                              FRFR_UPDATES, 
     &                              LBANDSLAVE, K474, BLR_U_COL)
        ENDIF
        FR_RANK = 0
        IF ((K480.GE.5).AND.(I.NE.1)) THEN
          IF (I.GT.FIRST_BLOCK) THEN
C           This block is outside the FRTRSM zone, we can
C           compress it early to avoid its decompression.                  
            IF (FRFR_UPDATES.EQ.0) THEN
C             No FR update, we compress A and initialize the
C             accumulator with it
              CALL SMUMPS_COMPRESS_FR_UPDATES(ACC_LRB,
     &                MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB,
     &                NFRONT, NIV, TOLEPS, TOL_OPT, KPERCENT, 
     &                COMPRESSED_FR, LorU, .FALSE.)
              MAX_ACC_RANK = ACC_LRB%K
              NEW_ACC_RANK = ACC_LRB%K
              FR_RANK = ACC_LRB%K
            ENDIF
          ENDIF
        ENDIF
        NB_DEC = FRFR_UPDATES
        DO JJ = 1, CURRENT_BLR
          J = J_ORDER(JJ)
          K_MAX = J_RANK(JJ)
          IF (LorU.EQ.0) THEN ! L panel update
            IF (LBANDSLAVE) THEN
              IND_L = I
              IF (K474.LT.2) THEN
                IND_U = CURRENT_BLR+1-J
              ELSE
                IND_U = J
              ENDIF
            ELSE
              IND_L = CURRENT_BLR+I-J
              IND_U = CURRENT_BLR+1-J
            ENDIF
          ELSE ! U panel update
            IND_L = CURRENT_BLR+1-J
            IND_U = CURRENT_BLR+I-J
          ENDIF
          CALL SMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &           IWHANDLER,
     &           0, ! L Panel
     &           J, BLR_L)
          IF (BLR_L(IND_L)%M.EQ.0) THEN
            CYCLE
          ENDIF
          IF (.NOT.LBANDSLAVE.OR.K474.LT.2) THEN
            CALL SMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &           IWHANDLER,
     &           1, ! U Panel
     &           J, BLR_U)
          ENDIF
        IF (K480.GE.3) THEN
          IF (ACC_LRB%K+K_MAX.GT.MAXI_RANK) THEN
C           FIXME We should still recompress when we're forced to
C           decompress                
            NB_DEC = JJ-1
            CALL SMUMPS_DECOMPRESS_ACC(ACC_LRB, MAXI_CLUSTER, 
     &            MAXI_RANK, A, LA, POSELT_INCB, NFRONT, NIV, LorU)
            COMPRESSED_FR = .FALSE.
            MAX_ACC_RANK = 0
          ENDIF
          OLD_ACC_RANK = ACC_LRB%K
        ENDIF
          CALL SMUMPS_LRGEMM4(MONE,
     &            BLR_U(IND_U), BLR_L(IND_L), ONE,
     &            A, LA, POSELT_INCB,
     &            NFRONT, 0, IFLAG, IERROR, 
     &            MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &            KPERCENT_RMB, MID_RANK, BUILDQ,
     &            (K480.GE.3), LorU=LorU, 
     &            LRB3=ACC_LRB, MAXI_RANK=MAXI_RANK,
     &            MAXI_CLUSTER=MAXI_CLUSTER
     &            )
        IF (IFLAG.LT.0) GOTO 100
        CALL UPD_FLOP_UPDATE(BLR_U(IND_U)%M, BLR_U(IND_U)%N, 
     &        BLR_U(IND_U)%K, BLR_U(IND_U)%ISLR, BLR_L(IND_L)%M, 
     &        BLR_L(IND_L)%N, BLR_L(IND_L)%K, BLR_L(IND_L)%ISLR,
     &            MIDBLK_COMPRESS, MID_RANK, BUILDQ, 
     &            .FALSE., (K480.GE.3))
C         Update K_RANK when K481=1 (later used by n-ary tree
C         recompression)
          IF ((MIDBLK_COMPRESS.GE.1).AND.BUILDQ) THEN
            J_RANK(JJ) = MID_RANK
          ENDIF
        IF (K480.GE.3) THEN
          NEW_ACC_RANK = NEW_ACC_RANK + ACC_LRB%K - OLD_ACC_RANK
          MAX_ACC_RANK = MAX(MAX_ACC_RANK, ACC_LRB%K - OLD_ACC_RANK)
          IF (K480.EQ.4) THEN
C           LUA with recompression
              IF ((K478.GT.0).AND.((ACC_LRB%K-MAX_ACC_RANK).GE.K478)) 
     &          THEN
C             Recompress accumulator  
              CALL SMUMPS_RECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &            MAXI_RANK, A, LA, POSELT_INCB, NFRONT, NIV,
     &            MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT_RMB,
     &            KPERCENT_LUA, NEW_ACC_RANK)
              MAX_ACC_RANK = ACC_LRB%K
            ENDIF
          ENDIF
        ENDIF
        IF ((K480.GE.5).AND.(I.NE.1)) THEN
          IF (I.GT.FIRST_BLOCK) THEN
C           This block is outside the FRTRSM zone, we can
C           compress it early to avoid its decompression.                  
            IF (JJ.EQ.FRFR_UPDATES) THEN
C             Last FR update, we compress FR part and initialize the
C             accumulator with it
              CALL SMUMPS_COMPRESS_FR_UPDATES(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB,
     &               NFRONT, NIV, TOLEPS, TOL_OPT, KPERCENT,
     &               COMPRESSED_FR, LorU, .FALSE.)
              MAX_ACC_RANK = ACC_LRB%K
              NEW_ACC_RANK = ACC_LRB%K
              IF (COMPRESSED_FR) THEN
                J_RANK(JJ) = ACC_LRB%K 
                NB_DEC = FRFR_UPDATES-1 
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        ENDDO
        IF (K480.GE.3) THEN
          IF ((K480.GE.5)) THEN
C FIXME Implement strategy 5/6 with comb tree recompression            
            IF (COMPRESSED_FR.OR.(K480.GE.6)) THEN  
C             Recompress accumulator
C             Even if NOT COMPRESSED_FR, we can still recompress the
C             accumulator to reduce the decompress flops: strategy
C             K480=6
              IF (ACC_LRB%K.GT.0) THEN
                IF (K478.EQ.-1) THEN
C                 K478=-1: star graph recompression  
                  IF (CURRENT_BLR-FRFR_UPDATES.GT.1) THEN
                    CALL SMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
                  ENDIF
                ELSEIF (K478.LE.-2) THEN
C                 K478=-n, n>=2: n-ary tree recompression  
                  IF (FRFR_UPDATES.GT.0) THEN
C                   We use J_RANK(FRFR_UPDATES) to store the rank of the
C                   compressed FR updates
                    allocate(POS_LIST(CURRENT_BLR-NB_DEC),stat=allocok)
                    IF (allocok .GT. 0) THEN
                       IFLAG  = -13
                       IERROR = CURRENT_BLR-NB_DEC
                       write(*,*) 'Allocation problem in BLR routine ',
     &                      'SMUMPS_BLR_UPD_PANEL_LEFT: ',
     &                      'not enough memory? memory requested = ',
     &                      IERROR
                       GOTO 100
                    ENDIF
                    POS_LIST(1) = 1
                    DO II = 1,CURRENT_BLR-NB_DEC-1
                      POS_LIST(II+1)=POS_LIST(II)+J_RANK(NB_DEC+II) 
                    ENDDO
                    CALL SMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB, KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               J_RANK(NB_DEC+1:CURRENT_BLR), POS_LIST,
     &               CURRENT_BLR-NB_DEC, 0)
                  ELSE
C                   We allocate a new RANK_LIST array to store FR_RANK
C                   (rank of the compressed FR updates) followed by
C                   J_RANK
                     allocate(POS_LIST(CURRENT_BLR+1),stat=allocok)
                     IF (allocok .GT. 0) THEN
                        IFLAG  = -13
                        IERROR = CURRENT_BLR+1
                        write(*,*) 'Allocation problem in BLR routine ',
     &                       'SMUMPS_BLR_UPD_PANEL_LEFT: ',
     &                       'not enough memory? memory requested = ',
     &                       IERROR
                        GOTO 100
                     ENDIF
                    POS_LIST(1) = 1
                    POS_LIST(2) = 1 + FR_RANK
                    DO II = 2,CURRENT_BLR
                      POS_LIST(II+1)=POS_LIST(II)+J_RANK(II-1) 
                    ENDDO
                    allocate(RANK_LIST(CURRENT_BLR+1),stat=allocok)
                    IF (allocok .GT. 0) THEN
                       IFLAG  = -13
                       IERROR = CURRENT_BLR+1
                       write(*,*) 'Allocation problem in BLR routine ',
     &                      'SMUMPS_BLR_UPD_PANEL_LEFT: ',
     &                      'not enough memory? memory requested = ',
     &                      IERROR
                       GOTO 100
                    ENDIF
                    RANK_LIST(1) = FR_RANK
                    DO II = 2,CURRENT_BLR+1
                      RANK_LIST(II) = J_RANK(II-1)
                    ENDDO
                    CALL SMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB, KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               RANK_LIST, POS_LIST,
     &               CURRENT_BLR+1, 0)
                    deallocate(RANK_LIST)
                  ENDIF
                  deallocate(POS_LIST)
                ENDIF
              ENDIF
            ENDIF
            MAXRANK = FLOOR(REAL(ACC_LRB%M*ACC_LRB%N)/REAL(ACC_LRB%M+
     &                                                     ACC_LRB%N))
            IF (COMPRESSED_FR.AND.(ACC_LRB%K.LE.MAXRANK)) THEN  
C             It's better to keep this block in LR form
              CALL ALLOC_LRB_FROM_ACC(ACC_LRB, NEXT_BLR(I-1), 
     &                    ACC_LRB%K, ACC_LRB%M, ACC_LRB%N, LorU,
     &                    IFLAG, IERROR, KEEP8)
              IF (IFLAG.LT.0) CYCLE
              ACC_LRB%K = 0
            ELSE
C FIXME       Differentiate the case where we did not succeed in
C             compressing the FR part (in which case we can decompress
C             the LR part and try compressing the sum FR+LR) and the
C             case where the rank of the accumulator exceeds MAXRANK (in
C             which case there's no need to try compressing again)
C             It's better to keep this block in FR form
              IF (I.NE.1) NEXT_BLR(I-1)%ISLR=.FALSE.
              CALL SMUMPS_DECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &              MAXI_RANK, A, LA, POSELT_INCB, NFRONT, NIV, LorU)
            ENDIF
          ELSE ! KEEP(480) < 5
C           Finished updating block, in K480=3,4 strategies we
C           decompress the accumulator.
            IF ((K480.EQ.4).AND.(K478.EQ.-1).AND.(ACC_LRB%K.GT.0)) THEN
C             K478=-1: star graph recompression  
              IF (CURRENT_BLR-FRFR_UPDATES.GT.1) THEN
                CALL SMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
              ENDIF
            ELSEIF ((K480.EQ.4).AND.(K478.LE.-2).AND.(ACC_LRB%K.GT.0)) 
     &      THEN
C             K478=-n, n>=2: n-ary tree recompression  
              allocate(POS_LIST(CURRENT_BLR-NB_DEC),stat=allocok)
              IF (allocok .GT. 0) THEN
                 IFLAG  = -13
                 IERROR = CURRENT_BLR-NB_DEC
                 GOTO 100
              ENDIF
              POS_LIST(1) = 1
              DO II = 1,CURRENT_BLR-NB_DEC-1
              POS_LIST(II+1)=POS_LIST(II)+J_RANK(NB_DEC+II) 
              ENDDO
              CALL SMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB, KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               J_RANK(NB_DEC+1:CURRENT_BLR), POS_LIST,
     &               CURRENT_BLR-NB_DEC, 0)
              deallocate(POS_LIST)
            ENDIF
            CALL SMUMPS_DECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &            MAXI_RANK, A, LA, POSELT_INCB, NFRONT, NIV, LorU)
          ENDIF
        ENDIF
 100    CONTINUE        
      ENDDO
#if ! defined(BLR_NOOPENMP)
!$OMP END DO NOWAIT
#endif
      END SUBROUTINE SMUMPS_BLR_UPD_PANEL_LEFT
      SUBROUTINE SMUMPS_BLR_UPD_CB_LEFT_LDLT(A, LA, POSELT, NFRONT,
     &        BEGS_BLR, BEGS_BLR_DYN, NB_INCB, NB_INASM, NASS,      
     &        IWHANDLER, 
     &        IW2, BLOCK, ACC_LUA,
     &        MAXI_CLUSTER, MAXI_RANK, NIV, IFLAG, IERROR,
     &        MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT_RMB,
     &        K480, K479, K478, KPERCENT_LUA, KPERCENT, KEEP8)
!$    USE OMP_LIB 
C    This subroutine does a left-looking BLR update of the contribution
C    block.         
C     -----------
C     Parameters
C     -----------
      INTEGER(8), intent(in)       :: LA
      REAL, intent(inout)       :: A(LA)
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_INCB, NB_INASM
      INTEGER, INTENT(IN)          :: NIV, IWHANDLER, MAXI_CLUSTER,
     &                                MAXI_RANK, K480, K479, K478, NASS,
     &                                KPERCENT_LUA, KPERCENT
      INTEGER, intent(inout)         :: IFLAG, IERROR
      INTEGER(8) :: KEEP8(150)
C     BEGS_BLR is BEGS_BLR_STATIC, used for updating the CB
C     In LDLT we also need BEGS_LR_DYNAMIC to apply scaling      
      INTEGER, DIMENSION(:) :: BEGS_BLR
      INTEGER, DIMENSION(:) :: BEGS_BLR_DYN
C     BLOCK of size (MAXI_CLUSTER, NOMP*MAXI_CLUSTER)
      REAL, INTENT(INOUT) :: BLOCK(MAXI_CLUSTER,*)
      INTEGER, intent(in) :: IW2(*)
      TYPE(LRB_TYPE), POINTER :: ACC_LUA(:)
C     MIDBLK_COMPRESS = whether to compress the middle product
C     Y1*X2 in A = A + X1*(Y1*X2)*Y2 (strategy suggested by CA)
      INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT_RMB, TOL_OPT
C     Threshold used for RRQR if MIDBLK_COMPRESS = TRUE 
      REAL,intent(in) :: TOLEPS
C     -----------
C     Local variables
C     -----------
      INTEGER :: M, N, allocok
      INTEGER :: I, II, J, K, KK, IND_L, IND_U, K_MAX, IBIS,
     &           K_ORDER(NB_INASM), K_RANK(NB_INASM), NB_DEC
      INTEGER, ALLOCATABLE :: POS_LIST(:), RANK_LIST(:)
      INTEGER(8) :: POSELT_BLOCK, POSELTD
      INTEGER :: MID_RANK, FRFR_UPDATES, MAXRANK, FR_RANK
      LOGICAL :: BUILDQ, COMPRESSED_FR
      TYPE(LRB_TYPE), POINTER :: BLR_L(:)
      TYPE(LRB_TYPE), POINTER :: ACC_LRB
      INTEGER :: OLD_ACC_RANK, MAX_ACC_RANK, NEW_ACC_RANK
      INTEGER :: OFFSET_IW
      INTEGER :: OMP_NUM
#if ! defined(BLR_NOOPENMP)
      INTEGER :: CHUNK
#endif
      REAL :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
      PARAMETER (ZERO=0.0E0)
      ACC_LRB => ACC_LUA(1)
      OMP_NUM = 0
C     Loop over all the potential BLR blocks in the CB
#if ! defined(BLR_NOOPENMP)
      CHUNK = 1
!$OMP DO SCHEDULE(DYNAMIC,CHUNK) 
!$OMP& PRIVATE(I, J, K, KK, POSELT_BLOCK, MID_RANK, BUILDQ,
!$OMP&         BLR_L, IND_U, IND_L, M, N, K_ORDER, K_RANK,
!$OMP&         K_MAX, OMP_NUM, ACC_LRB, POSELTD,
!$OMP&         MAX_ACC_RANK, OLD_ACC_RANK, NEW_ACC_RANK, 
!$OMP&         FRFR_UPDATES, FR_RANK, NB_DEC, II)
#endif
      DO IBIS = 1,NB_INCB*(NB_INCB+1)/2
        IF (IFLAG.LT.0) CYCLE     
C       Determining I,J from IBIS
C       Notice that IBIS = I(I-1)/2 + J - 1      
        I = CEILING((1.0D0+SQRT(1.0D0+8.0D0*dble(IBIS)))/2.0D0)-1
        J = IBIS - I*(I-1)/2
        I = I+NB_INASM
        J = J+NB_INASM
C       DO I = NB_INASM + 1 , NB_INASM + NB_INCB
#if ! defined(BLR_NOOPENMP)
        OMP_NUM = 0
!$      OMP_NUM = OMP_GET_THREAD_NUM() 
        ACC_LRB => ACC_LUA(OMP_NUM+1)
#endif
C         DO J = NB_INASM + 1 , I
        MAX_ACC_RANK = 0
        NEW_ACC_RANK = 0
        M = BEGS_BLR(I+1)-BEGS_BLR(I)
        N = BEGS_BLR(J+1)-BEGS_BLR(J)
        POSELT_BLOCK = POSELT + int(NFRONT,8)*int(BEGS_BLR(I)-1,8) + 
     &           int(BEGS_BLR(J)-1,8)
C       M and N are inversed because ACC_LRB is a LR approximation of
C       A^T, not A 
        ACC_LRB%M = N
        ACC_LRB%N = M
        IF (K480.EQ.2) THEN
C         Pure left-looking facto, no LUA done 
          DO K = 1, NB_INASM
            K_ORDER(K) = K 
          ENDDO
        ELSE
C         LUA activated 
          CALL SMUMPS_GET_LUA_ORDER(NB_INASM, K_ORDER, K_RANK, 
     &                             IWHANDLER, 
     &                             1, 1, I, J, ! SYM, CB, I, J
     &                             FRFR_UPDATES)
        ENDIF
        FR_RANK = 0
        IF ((K480.GE.5).AND.(I.NE.J)) THEN
          IF (FRFR_UPDATES.EQ.0) THEN
C           No FR update, we compress A and initialize the
C           accumulator with it
            CALL SMUMPS_COMPRESS_FR_UPDATES(ACC_LRB,
     &                MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &                NFRONT, NIV, TOLEPS, TOL_OPT, KPERCENT, 
     &                COMPRESSED_FR, 0, .TRUE.)
            FR_RANK = ACC_LRB%K
            MAX_ACC_RANK = ACC_LRB%K
            NEW_ACC_RANK = ACC_LRB%K
          ENDIF
        ENDIF
        NB_DEC = FRFR_UPDATES
        DO KK = 1, NB_INASM
          K = K_ORDER(KK)
          K_MAX = K_RANK(KK)
C         Position of first diagonal value for applying scaling D from LDLT
          POSELTD = POSELT + int(NFRONT,8) * int(BEGS_BLR_DYN(K)-1,8)
     &         + int(BEGS_BLR_DYN(K) - 1,8)
          OFFSET_IW = BEGS_BLR_DYN(K)
          IND_L = I-K
          IND_U = J-K
          CALL SMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &             IWHANDLER,
     &             0, ! L Panel
     &             K, BLR_L)
          IF (BLR_L(IND_L)%M.EQ.0) THEN
            CYCLE
          ENDIF
          IF (K480.GE.3) THEN ! K480 >= 3: LUA
            IF (ACC_LRB%K+K_MAX.GT.MAXI_RANK) THEN
C             FIXME We should still recompress when we're forced to
C             decompress                
              NB_DEC = KK-1
              CALL SMUMPS_DECOMPRESS_ACC(ACC_LRB, 
     &              MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK, 
     &              NFRONT, NIV, 2)
              COMPRESSED_FR = .FALSE.
              MAX_ACC_RANK = 0
            ENDIF
            OLD_ACC_RANK = ACC_LRB%K
          ENDIF
          CALL SMUMPS_LRGEMM4(MONE,
     &              BLR_L(IND_U), BLR_L(IND_L), ONE,
     &              A, LA, POSELT_BLOCK,
     &              NFRONT, 1, IFLAG, IERROR, 
     &              MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &              KPERCENT_RMB, MID_RANK, BUILDQ,
     &              (K480.GE.3), LorU=2, 
     &              LRB3=ACC_LRB, MAXI_RANK=MAXI_RANK,
     &              MAXI_CLUSTER=MAXI_CLUSTER,
     &              DIAG=A(POSELTD), LD_DIAG=NFRONT, 
     &              IW2=IW2(OFFSET_IW),
     &              BLOCK=BLOCK(1:MAXI_CLUSTER,OMP_NUM*MAXI_CLUSTER+1))
          IF (IFLAG.LT.0) GOTO 100
            CALL UPD_FLOP_UPDATE(BLR_L(IND_U)%M, BLR_L(IND_U)%N, 
     &        BLR_L(IND_U)%K, BLR_L(IND_U)%ISLR, BLR_L(IND_L)%M, 
     &        BLR_L(IND_L)%N, BLR_L(IND_L)%K, BLR_L(IND_L)%ISLR,
     &           MIDBLK_COMPRESS, MID_RANK, BUILDQ, 
     &           (I.EQ.J), (K480.GE.3))
C         Update K_RANK when K481=1 (later used by n-ary tree
C         recompression)
          IF ((MIDBLK_COMPRESS.GE.1).AND.BUILDQ) THEN
            K_RANK(KK) = MID_RANK
          ENDIF
          IF (K480.GE.3) THEN
            NEW_ACC_RANK = NEW_ACC_RANK + ACC_LRB%K - OLD_ACC_RANK
            MAX_ACC_RANK = MAX(MAX_ACC_RANK, ACC_LRB%K - OLD_ACC_RANK)
            IF (K480.EQ.4) THEN
C             LUA with recompression
              IF ((K478.GT.0).AND.((ACC_LRB%K-MAX_ACC_RANK).GE.K478)) 
     &          THEN
C               Recompress accumulator  
                IF (ACC_LRB%K.GT.0) THEN
                  CALL SMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS,
     &               TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
                  MAX_ACC_RANK = ACC_LRB%K
                ENDIF
              ENDIF
            ENDIF
            IF ((K480.GE.5).AND.(I.NE.J)) THEN
              IF (KK.EQ.FRFR_UPDATES) THEN
C               Last FR update, we compress FR part and initialize the
C               accumulator with it
                CALL SMUMPS_COMPRESS_FR_UPDATES(ACC_LRB,
     &                   MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &                   NFRONT, NIV, TOLEPS, TOL_OPT, KPERCENT, 
     &                   COMPRESSED_FR, 0, .TRUE.)
                IF (COMPRESSED_FR) THEN
                  K_RANK(KK) = ACC_LRB%K 
                  NB_DEC = FRFR_UPDATES-1 
                ENDIF
                MAX_ACC_RANK = ACC_LRB%K
                NEW_ACC_RANK = ACC_LRB%K
              ENDIF
            ENDIF
          ENDIF
        END DO
        IF (K480.GE.3) THEN
          IF ((K480.GE.5)) THEN
C FIXME Implement strategy 5/6 with comb tree recompression            
            IF (COMPRESSED_FR.OR.(K480.GE.6)) THEN  
C             Recompress accumulator
C             Even if NOT COMPRESSED_FR, we can still recompress the
C             accumulator to reduce the decompress flops: strategy
C             K480=6
              IF (ACC_LRB%K.GT.0) THEN
                IF (K478.EQ.-1) THEN
C                 K478=-1: star graph recompression  
                  IF (NB_INASM-FRFR_UPDATES.GT.1) THEN               
                    CALL SMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS,
     &               TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
                  ENDIF
                ELSEIF (K478.LE.-2) THEN
C                 K478=-n, n>=2: n-ary tree recompression  
                  IF (FRFR_UPDATES.GT.0) THEN
C                   We use K_RANK(FRFR_UPDATES) to store the rank of the
C                   compressed FR updates
                    allocate(POS_LIST(NB_INASM-NB_DEC),stat=allocok)
                    IF (allocok .GT. 0) THEN
                       IFLAG  = -13
                       IERROR = NB_INASM-NB_DEC
                       write(*,*) 'Allocation problem in BLR routine ',
     &                      'SMUMPS_BLR_UPD_CB_LEFT_LDLT: ',
     &                      'not enough memory? memory requested = ',
     &                      IERROR
                       GOTO 100
                    ENDIF
                    POS_LIST(1) = 1
                    DO II = 1,NB_INASM-NB_DEC-1
                      POS_LIST(II+1)=POS_LIST(II)+K_RANK(NB_DEC+II) 
                    ENDDO
                    CALL SMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               K_RANK(NB_DEC+1:NB_INASM), POS_LIST,
     &               NB_INASM-NB_DEC, 0)
                  ELSE
C                   We allocate a new RANK_LIST array to store FR_RANK
C                   (rank of the compressed FR updates) followed by
C                   K_RANK
                    allocate(POS_LIST(NB_INASM+1),stat=allocok)
                    IF (allocok .GT. 0) THEN
                       IFLAG  = -13
                       IERROR = NB_INASM+1
                       write(*,*) 'Allocation problem in BLR routine ',
     &                      'SMUMPS_BLR_UPD_CB_LEFT_LDLT: ',
     &                      'not enough memory? memory requested = ',
     &                      IERROR
                       GOTO 100
                    ENDIF
                    POS_LIST(1) = 1
                    POS_LIST(2) = 1 + FR_RANK
                    DO II = 2,NB_INASM
                      POS_LIST(II+1)=POS_LIST(II)+K_RANK(II-1) 
                    ENDDO
                    allocate(RANK_LIST(NB_INASM+1),stat=allocok)
                    IF (allocok .GT. 0) THEN
                       IFLAG  = -13
                       IERROR = NB_INASM+1
                       write(*,*) 'Allocation problem in BLR routine ',
     &                      'SMUMPS_BLR_UPD_CB_LEFT_LDLT: ',
     &                      'not enough memory? memory requested = ',
     &                      IERROR
                       GOTO 100
                    ENDIF
                    RANK_LIST(1) = FR_RANK
                    DO II = 2,NB_INASM+1
                      RANK_LIST(II) = K_RANK(II-1)
                    ENDDO
                    CALL SMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               RANK_LIST, POS_LIST,
     &               NB_INASM+1, 0)
                    deallocate(RANK_LIST)
                  ENDIF
                  deallocate(POS_LIST)
                ENDIF
              ENDIF
            ENDIF
            MAXRANK = FLOOR(REAL(ACC_LRB%M*ACC_LRB%N)/REAL(ACC_LRB%M+
     &                                                     ACC_LRB%N))
            IF (COMPRESSED_FR.AND.(ACC_LRB%K.LE.MAXRANK)) THEN  
C             It's better to keep this block in LR form
C             We decompress anyway (until we have a BLR solve) but do
C             not count the cost of doing so             
C FIXME       Probably better (for timing purposes) to isolate this
C             decompression outside the main PARALLEL loop              
C FIXME       Even better would be to have a CB BLR panel and decompress
C             it at the end like the L/U panels                 
              CALL SMUMPS_DECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &              MAXI_RANK, A, LA, POSELT_BLOCK, NFRONT, NIV, 2,
     &              COUNT_FLOPS=.FALSE.)
            ELSE
C             It's better to keep this block in FR form
              CALL SMUMPS_DECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &              MAXI_RANK, A, LA, POSELT_BLOCK, NFRONT, NIV, 2)
            ENDIF
          ELSE ! KEEP(480) < 5
C           Finished updating block, in K480=3,4 strategies we
C           decompress the accumulator.
            IF ((K480.EQ.4).AND.(K478.EQ.-1).AND.(ACC_LRB%K.GT.0)) THEN
C             K478=-1: star graph recompression  
              IF (NB_INASM-FRFR_UPDATES.GT.1) THEN               
                CALL SMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
              ENDIF
            ELSEIF ((K480.EQ.4).AND.(K478.LE.-2).AND.(ACC_LRB%K.GT.0)) 
     &           THEN
C             K478=-n, n>=2: n-ary tree recompression  
              allocate(POS_LIST(NB_INASM-NB_DEC),stat=allocok)
              IF (allocok .GT. 0) THEN
                 IFLAG  = -13
                 IERROR = NB_INASM-NB_DEC
                 GOTO 100
              ENDIF
              POS_LIST(1) = 1
              DO II = 1,NB_INASM-NB_DEC-1
                POS_LIST(II+1)=POS_LIST(II)+K_RANK(NB_DEC+II) 
              ENDDO
              CALL SMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &               KEEP8, NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, 
     &               TOL_OPT, KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               K_RANK(NB_DEC+1:NB_INASM), POS_LIST,
     &               NB_INASM-NB_DEC, 0)
              deallocate(POS_LIST)
            ENDIF
            CALL SMUMPS_DECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &            MAXI_RANK, A, LA, POSELT_BLOCK, NFRONT, NIV, 2)
          ENDIF                               
        ENDIF
 100    CONTINUE     
      END DO
#if ! defined(BLR_NOOPENMP)
!$OMP END DO 
#endif
      END SUBROUTINE SMUMPS_BLR_UPD_CB_LEFT_LDLT
      SUBROUTINE SMUMPS_BLR_UPD_CB_LEFT(A, LA, POSELT, NFRONT,
     &        BEGS_BLR, BEGS_BLR_U, NB_ROWS, NB_INCB, NB_INASM, NASS,
     &        IWHANDLER, NIV, LBANDSLAVE, IFLAG, IERROR,
     &        MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT_RMB,
     &        ACC_LUA, K480, K479, K478, KPERCENT_LUA,
     &        K452,
     &        KPERCENT,
     &        MAXI_CLUSTER, MAXI_RANK,
     &        K474, FSorCB, BLR_U_COL, COMPRESS_CB, CB_LRB, KEEP8)
!$    USE OMP_LIB
C    This subroutine does a left-looking BLR update of the contribution
C    block.         
C     -----------
C     Parameters
C     -----------
      INTEGER(8), intent(in)       :: LA
      REAL, intent(inout)       :: A(LA)
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_ROWS, NB_INCB, NB_INASM
      INTEGER, INTENT(IN)          :: NIV, IWHANDLER, MAXI_CLUSTER,
     &                                MAXI_RANK, KPERCENT_LUA, KPERCENT
      INTEGER, INTENT(IN)          :: K480, K479, K478, NASS, K474,
     &                                FSorCB
      INTEGER, intent(inout)         :: IFLAG, IERROR
      INTEGER, DIMENSION(:) :: BEGS_BLR, BEGS_BLR_U
#if defined(MUMPS_NOF2003)
      TYPE(LRB_TYPE), POINTER :: CB_LRB(:,:)
#else
      TYPE(LRB_TYPE), POINTER, intent(inout) :: CB_LRB(:,:)
#endif
      TYPE(LRB_TYPE), POINTER :: ACC_LUA(:), BLR_U_COL(:)
      INTEGER(8) :: KEEP8(150)
C     MIDBLK_COMPRESS = whether to compress the middle product
C     Y1*X2 in A = A + X1*(Y1*X2)*Y2 (strategy suggested by CA)
      INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT_RMB, TOL_OPT
C     Threshold used for RRQR if MIDBLK_COMPRESS = TRUE 
      REAL,intent(in) :: TOLEPS
      LOGICAL, intent(in) :: LBANDSLAVE, COMPRESS_CB
      INTEGER, intent(in) :: K452
C     -----------
C     Local variables
C     -----------
      INTEGER :: M, N, allocok
      INTEGER :: I, II, J, K, KK, IND_L, IND_U, IBIS,
     &           K_ORDER(NB_INASM), K_RANK(NB_INASM)
      INTEGER, ALLOCATABLE :: POS_LIST(:), RANK_LIST(:)
      INTEGER(8) :: POSELT_BLOCK
      INTEGER :: MID_RANK, K_MAX, FRFR_UPDATES, NB_DEC
      LOGICAL :: BUILDQ, COMPRESSED_FR
      TYPE(LRB_TYPE), POINTER :: BLR_U(:), BLR_L(:)
      TYPE(LRB_TYPE), POINTER :: ACC_LRB, LRB
      INTEGER :: OLD_ACC_RANK, MAX_ACC_RANK, NEW_ACC_RANK, MAXRANK,
     &           FR_RANK
#if ! defined(BLR_NOOPENMP)
      INTEGER :: OMP_NUM
      INTEGER :: CHUNK
#endif
      REAL :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
      PARAMETER (ZERO=0.0E0)
      ACC_LRB => ACC_LUA(1)
C     Loop over all the potential BLR blocks in the CB
#if ! defined(BLR_NOOPENMP)
      CHUNK = 1
!$OMP DO SCHEDULE(DYNAMIC,CHUNK) 
!$OMP& PRIVATE(I, J, K, KK, POSELT_BLOCK, MID_RANK, BUILDQ,
!$OMP&         BLR_U, BLR_L, IND_U, IND_L, M, N, 
!$OMP&         ACC_LRB, OMP_NUM, K_MAX, K_ORDER, K_RANK,
!$OMP&         MAX_ACC_RANK, OLD_ACC_RANK, NEW_ACC_RANK,
!$OMP&         FRFR_UPDATES, LRB)
#endif
      DO IBIS = 1,NB_ROWS*NB_INCB
        IF (IFLAG.LT.0) CYCLE     
C       Determining I,J from IBIS
        I = (IBIS-1)/NB_INCB+1
        J = IBIS - (I-1)*NB_INCB
        IF (.NOT.LBANDSLAVE) THEN
          I = I+NB_INASM
        ENDIF
        J = J+NB_INASM
C       DO I = NB_INASM + 1 , NB_INASM + NB_INCB
#if ! defined(BLR_NOOPENMP)
        OMP_NUM=0 ! in case compilation is w/o openmp
!$      OMP_NUM = OMP_GET_THREAD_NUM() 
        ACC_LRB => ACC_LUA(OMP_NUM+1)
#endif
C       DO J = NB_INASM + 1 , NB_INASM + NB_INCB
        MAX_ACC_RANK = 0
        NEW_ACC_RANK = 0
        IF (LBANDSLAVE) THEN
          M = BEGS_BLR(I+2)-BEGS_BLR(I+1)
          IF (K474.EQ.1) THEN
C           LUA-NIV2 V1 strategy
            POSELT_BLOCK = POSELT + int(NFRONT,8)*int(BEGS_BLR(I+1)-1,8)
     &           +int(NASS,8) + int(BEGS_BLR_U(J-NB_INASM+1)-1,8)
            N = BEGS_BLR_U(J-NB_INASM+2)-BEGS_BLR_U(J-NB_INASM+1)
          ELSEIF (K474.GE.2) THEN
C           LUA-NIV2 V2 strategy
            BLR_U => BLR_U_COL
            POSELT_BLOCK = POSELT + int(NFRONT,8)*int(BEGS_BLR(I+1)-1,8)
     &          + int(NASS-1,8)
            N = BEGS_BLR_U(3)-BEGS_BLR_U(2)
          ELSE
C     LUA-NIV2 V0 shoudl not do LUA on slave
            write(*,*) 'Internal error in SMUMPS_BLR_UPD_CB_LEFT',
     &        LBANDSLAVE,K474       
            CALL MUMPS_ABORT()
          ENDIF
        ELSE
          M = BEGS_BLR(I+1)-BEGS_BLR(I)
          POSELT_BLOCK = POSELT + int(NFRONT,8)*int(BEGS_BLR(I)-1,8) + 
     &           int(BEGS_BLR_U(J)-1,8)
          N = BEGS_BLR_U(J+1)-BEGS_BLR_U(J)
        ENDIF
C       M and N are inversed because ACC_LRB is a LR approximation of
C       A^T, not A 
        ACC_LRB%M = N
        ACC_LRB%N = M
        IF (K480.EQ.2) THEN
C         Pure left-looking facto, no LUA done 
          DO K = 1, NB_INASM
            K_ORDER(K) = K 
          ENDDO
        ELSE
C         LUA activated 
          CALL SMUMPS_GET_LUA_ORDER(NB_INASM, K_ORDER, K_RANK,
     &                             IWHANDLER, 
     &                             0, 1, I, J, ! SYM, CB, I, J
     &                             FRFR_UPDATES,
     &                             LBANDSLAVE, K474, BLR_U_COL)
        ENDIF
        COMPRESSED_FR = .FALSE.
        FR_RANK = 0
        DO KK = 1, NB_INASM
          IF ((K480.GE.5.OR.COMPRESS_CB).AND.I.NE.J) THEN
            IF (KK-1.EQ.FRFR_UPDATES) THEN
C             Last FR update, we compress FR part and initialize the
C             accumulator with it
C             Note that if FRFR_UPDATES = NB_INASM (i.e. all updates
C             are FR-FR) then this routine is never called and the
C             corresponding CB block is not compressed               
              CALL SMUMPS_COMPRESS_FR_UPDATES(ACC_LRB,
     &                 MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &                 NFRONT, NIV, TOLEPS, TOL_OPT, KPERCENT, 
     &                 COMPRESSED_FR, 0, .TRUE.)
              IF (COMPRESSED_FR) THEN
                K_RANK(KK) = ACC_LRB%K 
                NB_DEC = FRFR_UPDATES-1 
              ENDIF
              MAX_ACC_RANK = ACC_LRB%K
              NEW_ACC_RANK = ACC_LRB%K
              FR_RANK = ACC_LRB%K
            ENDIF
          ENDIF
          K = K_ORDER(KK)
          K_MAX = K_RANK(KK)
          IF (LBANDSLAVE) THEN
            IND_L = I
            IF (K474.LT.2) THEN
              IND_U = J-K
            ELSE
              IND_U = K
            ENDIF
          ELSE
            IND_L = I-K
            IND_U = J-K
          ENDIF
          CALL SMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &             IWHANDLER,
     &             0, ! L Panel
     &             K, BLR_L)
          IF (BLR_L(IND_L)%M.EQ.0) THEN
            CYCLE
          ENDIF
          IF (.NOT.LBANDSLAVE.OR.K474.LT.2) THEN
            CALL SMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &             IWHANDLER,
     &             1, ! U Panel
     &             K, BLR_U)
          ENDIF
          IF (K480.GE.3) THEN
            IF (ACC_LRB%K+K_MAX.GT.MAXI_RANK) THEN
C             Since we have been forced to decompress the accumulator back to
C             FR by lack of space we set COMPRESSED_FR back to FALSE. 
              COMPRESSED_FR = .FALSE.
              NB_DEC = KK-1
              CALL SMUMPS_DECOMPRESS_ACC(ACC_LRB, 
     &              MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK, 
     &              NFRONT, NIV, 2)
              MAX_ACC_RANK = 0
            ENDIF
            OLD_ACC_RANK = ACC_LRB%K
          ENDIF
            CALL SMUMPS_LRGEMM4(MONE,
     &              BLR_U(IND_U), BLR_L(IND_L), ONE,
     &              A, LA, POSELT_BLOCK,
     &              NFRONT, 0, IFLAG, IERROR,
     &              MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &              KPERCENT_RMB, MID_RANK, BUILDQ,
     &              (K480.GE.3), LorU=2,
     &              LRB3=ACC_LRB, MAXI_RANK=MAXI_RANK,
     &              MAXI_CLUSTER=MAXI_CLUSTER)
          IF (IFLAG.LT.0) GOTO 100
            CALL UPD_FLOP_UPDATE(BLR_U(IND_U)%M, BLR_U(IND_U)%N, 
     &        BLR_U(IND_U)%K, BLR_U(IND_U)%ISLR, BLR_L(IND_L)%M, 
     &        BLR_L(IND_L)%N, BLR_L(IND_L)%K, BLR_L(IND_L)%ISLR,
     &           MIDBLK_COMPRESS, MID_RANK, BUILDQ, 
     &           .FALSE., (K480.GE.3))
C         Update K_RANK when K481=1 (later used by n-ary tree
C         recompression)
          IF ((MIDBLK_COMPRESS.GE.1).AND.BUILDQ) THEN
            K_RANK(KK) = MID_RANK
          ENDIF
          IF (K480.GE.3) THEN
            NEW_ACC_RANK = NEW_ACC_RANK + ACC_LRB%K - OLD_ACC_RANK
            MAX_ACC_RANK = MAX(MAX_ACC_RANK, ACC_LRB%K - OLD_ACC_RANK)
            IF (K480.EQ.4) THEN
C             LUA with recompression
              IF ((K478.GT.0).AND.((ACC_LRB%K-MAX_ACC_RANK).GE.K478)) 
     &          THEN
C               Recompress accumulator  
                CALL SMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &                 MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &                 NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS,
     &                 TOL_OPT,
     &                 KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
                MAX_ACC_RANK = ACC_LRB%K
              ENDIF
            ENDIF
          ENDIF
        END DO
        IF (K480.GE.3) THEN
          IF (K480.GE.5.OR.COMPRESS_CB) THEN
            IF (K480.GE.5.AND.(COMPRESSED_FR.OR.K480.GE.6)) THEN  
C             Recompress accumulator
C             Even if NOT COMPRESSED_FR, we can still recompress the
C             accumulator to reduce the decompress flops: strategy
C             K480=6
              IF (ACC_LRB%K.GT.0) THEN
                IF (K478.EQ.-1) THEN
C                 K478=-1: star graph recompression  
                  IF (NB_INASM-FRFR_UPDATES.GT.1) THEN               
                    CALL SMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
                  ENDIF
                ELSEIF (K478.LE.-2) THEN
C                 K478=-n, n>=2: n-ary tree recompression  
                  IF (FRFR_UPDATES.GT.0) THEN
C                   We use K_RANK(FRFR_UPDATES) to store the rank of the
C                   compressed FR updates
                    allocate(POS_LIST(NB_INASM-NB_DEC),stat=allocok)
                    IF (allocok .GT. 0) THEN
                       IFLAG  = -13
                       IERROR = NB_INASM-NB_DEC
                       GOTO 100
                    ENDIF
                    POS_LIST(1) = 1
                    DO II = 1,NB_INASM-NB_DEC-1
                      POS_LIST(II+1)=POS_LIST(II)+K_RANK(NB_DEC+II) 
                    ENDDO
                    CALL SMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               K_RANK(NB_DEC+1:NB_INASM), POS_LIST,
     &               NB_INASM-NB_DEC, 0)
                  ELSE
C                   We allocate a new RANK_LIST array to store FR_RANK
C                   (rank of the compressed FR updates) followed by
C                   K_RANK
                    allocate(POS_LIST(NB_INASM+1),stat=allocok)
                    IF (allocok .GT. 0) THEN
                       IFLAG  = -13
                       IERROR = NB_INASM+1
                       GOTO 100
                    ENDIF
                    POS_LIST(1) = 1
                    POS_LIST(2) = 1 + FR_RANK
                    DO II = 2,NB_INASM
                      POS_LIST(II+1)=POS_LIST(II)+K_RANK(II-1) 
                    ENDDO
                    allocate(RANK_LIST(NB_INASM+1),stat=allocok)
                    IF (allocok .GT. 0) THEN
                       IFLAG  = -13
                       IERROR = NB_INASM+1
                       GOTO 100
                    ENDIF
                    RANK_LIST(1) = FR_RANK
                    DO II = 2,NB_INASM+1
                      RANK_LIST(II) = K_RANK(II-1)
                    ENDDO
                    CALL SMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               RANK_LIST, POS_LIST,
     &               NB_INASM+1, 0)
                    deallocate(RANK_LIST)
                  ENDIF
                  deallocate(POS_LIST)
                ENDIF
              ENDIF
            ENDIF
            MAXRANK = FLOOR(REAL(ACC_LRB%M*ACC_LRB%N)/REAL(ACC_LRB%M+
     &                                                     ACC_LRB%N))
            IF (COMPRESSED_FR.AND.(ACC_LRB%K.LE.MAXRANK)) THEN
C             It's better to keep this block in LR form
              LRB => CB_LRB(I-NB_INASM,J-NB_INASM)
              CALL ALLOC_LRB_FROM_ACC(ACC_LRB, LRB,
     &                       ACC_LRB%K, ACC_LRB%M, ACC_LRB%N, 0,
     &                       IFLAG, IERROR, KEEP8)
              CALL UPD_MRY_CB_LRGAIN(LRB%M, LRB%N, LRB%K
     &                          , K452 .GT. 0, LRB%K_STOC, 
     &                    NPREC_STOC, COSTS_PREC
     &                              )
              ACC_LRB%K = 0
              IF (IFLAG.LT.0) GOTO 100
            ELSE
C             It's better to keep this block in FR form
              CALL SMUMPS_DECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &              MAXI_RANK, A, LA, POSELT_BLOCK, NFRONT, NIV, 2)
              LRB => CB_LRB(I-NB_INASM,J-NB_INASM)
              CALL ALLOC_LRB(LRB, ACC_LRB%K, ACC_LRB%N, ACC_LRB%M, 
     &                 .FALSE., IFLAG, IERROR, KEEP8)
              IF (IFLAG.LT.0) GOTO 100
              DO II=1,ACC_LRB%N
                LRB%Q(II,1:ACC_LRB%M) =
     &          A( POSELT_BLOCK+int((II-1),8)*int(NFRONT,8) :
     &            POSELT_BLOCK+int((II-1),8)*int(NFRONT,8)
     &                        +int(ACC_LRB%M-1,8) )
              END DO  
            ENDIF
          ELSE ! KEEP(480) < 5
C           Finished updating block, in K480=3,4 strategies we
C           decompress the accumulator.
            IF ((K480.EQ.4).AND.(K478.EQ.-1).AND.(ACC_LRB%K.GT.0)) THEN
C             K478=-1: star graph recompression  
              IF (NB_INASM-FRFR_UPDATES.GT.1) THEN               
                CALL SMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
              ENDIF
            ELSEIF ((K480.EQ.4).AND.(K478.LE.-2).AND.(ACC_LRB%K.GT.0)) 
     &        THEN
C             IF (IBIS.EQ.2.AND.ACC_LRB%K.EQ.53) THEN
C             K478=-n, n>=2: n-ary tree recompression  
              allocate(POS_LIST(NB_INASM-NB_DEC),stat=allocok)
              IF (allocok .GT. 0) THEN
                 IFLAG  = -13
                 IERROR = NB_INASM-NB_DEC
                 GOTO 100
              ENDIF
              POS_LIST(1) = 1
              DO II = 1,NB_INASM-NB_DEC-1
                POS_LIST(II+1)=POS_LIST(II)+K_RANK(NB_DEC+II) 
              ENDDO
              CALL SMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &               KEEP8,NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, 
     &               TOL_OPT, KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               K_RANK(NB_DEC+1:NB_INASM), POS_LIST,
     &               NB_INASM-NB_DEC, 0)
              deallocate(POS_LIST)
            ENDIF
            CALL SMUMPS_DECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &            MAXI_RANK, A, LA, POSELT_BLOCK, NFRONT, NIV, 2)
          ENDIF
        ENDIF                               
 100    CONTINUE       
      END DO
#if ! defined(BLR_NOOPENMP)
!$OMP END DO 
#endif
#if ! defined(BLR_NOOPENMP)
!$OMP MASTER
#endif
      IF (COMPRESS_CB) THEN
        CALL UPD_MRY_CB_FR(NFRONT-NASS, NFRONT-NASS, 0)
      ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP END MASTER
#endif
      END SUBROUTINE SMUMPS_BLR_UPD_CB_LEFT
      SUBROUTINE SMUMPS_DECOMPRESS_PANEL(A, LA, POSELT, LDA11,
     &        LDA21, COPY_DENSE_BLOCKS,
     &        BEGS_BLR_DIAG, BEGS_BLR_FIRST_OFFDIAG,
     &        NB_BLR, BLR_PANEL, CURRENT_BLR, DIR, DECOMP_TIMER,
     &        BEG_I_IN, END_I_IN, ONLY_NELIM_IN
     &        )
!$    USE OMP_LIB 
C     ------------
C     Purpose
C     -------
C     We promote a panel of (NB_BLR-CURRENT_BLR-1) off diagonal BLR blocks stored
C     in BLR_PANEL in the corresponding band of the front located at A(POSELT).
C      IF called by a slave  (COPY_DENSE_BLOCKS.EQV.TRUE.) then
C      also non LR blocks are copied back in band A (POSELT ..)
C     BEGS_BLR_DIAG: Position in front of the beginning of the diagonal block with respect to NFRONT
C     BEGS_BLR_FIRST_OFFDIAG : Position in front of the beginning of first off_diagonal block
C     -----------
C     Parameters
C     -----------
      INTEGER(8), intent(in)       :: LA
      REAL, intent(inout)       :: A(LA)
      INTEGER(8), intent(in)       :: POSELT 
      LOGICAL, intent(in)          :: COPY_DENSE_BLOCKS  
      INTEGER, intent(in)          :: NB_BLR, CURRENT_BLR
      INTEGER, intent(in)          :: BEGS_BLR_DIAG, 
     &                                BEGS_BLR_FIRST_OFFDIAG
C     TYPE(LRB_TYPE), intent(inout) :: BLR_PANEL(NB_BLR-CURRENT_BLR)
      TYPE(LRB_TYPE), intent(inout) :: BLR_PANEL(:)
C     V = Vertical panel, H = Horizontal panel
      CHARACTER(len=1) :: DIR
C     LDA11: leading dimension for L for rows <= NPIV
C     LDA21: leading dimension for L for rows > NPIV
C     In case of call from solve, leading dimensions are LDA11=NFRONT
C     and LDA21=NPIV in LU, and LDA11=LDA21=NPIV in LDLT.
C     During factorization, all leading dimensions are NFRONT.
      INTEGER, intent(in) :: LDA11, LDA21
      INTEGER, intent(in) :: DECOMP_TIMER
C     BEG_I and END_I are optional parameters which control which blocks
C     of the panel are decompressed; if not given the entire panel is
C     decompressed.
C     ONLY_NELIM_IN is also optional; it is used in the context of the UCFS
C     factorization; if given, only the ONLY_NELIM last rows (U) / columns
C     (L) of each block are decompressed.
      INTEGER,OPTIONAL,intent(in) :: BEG_I_IN, END_I_IN, ONLY_NELIM_IN
C     -----------
C     Local variables
C     -----------
      INTEGER :: IP, M, N, BIP, BIP_START, BEG_I, END_I, ONLY_NELIM
#if ! defined(BLR_NOOPENMP)
      INTEGER :: CHUNK
#endif
      INTEGER :: K, I
      DOUBLE PRECISION :: PROMOTE_COST
      INTEGER(8) :: POSELT_BLOCK, LD_BLK_IN_FRONT
      REAL :: ONE, ALPHA, ZERO
      PARAMETER (ONE = 1.0E0, ALPHA=-1.0E0)
      PARAMETER (ZERO = 0.0E0)
      IF(present(BEG_I_IN)) THEN
        BEG_I = BEG_I_IN
      ELSE
        BEG_I = CURRENT_BLR+1
      ENDIF
      IF(present(END_I_IN)) THEN
        END_I = END_I_IN
      ELSE
        END_I = NB_BLR
      ENDIF
      IF(present(ONLY_NELIM_IN)) THEN
        ONLY_NELIM = ONLY_NELIM_IN
C       For now we only count flops in UCFS factorization case        
      ELSE
        ONLY_NELIM = 0
      ENDIF
      LD_BLK_IN_FRONT = int(LDA11,8)
      ! pointer to first off diag block
      BIP_START       = BEGS_BLR_FIRST_OFFDIAG
      IF (BEG_I .NE. CURRENT_BLR+1) THEN
        DO I = 1, BEG_I - CURRENT_BLR - 1
          BIP_START  = BIP_START +  BLR_PANEL(I)%M
        ENDDO
      ENDIF
#if defined(BLR_NOOPENMP)
C     BIP is incremented in routine
      BIP             = BIP_START
#endif
#if ! defined(BLR_NOOPENMP)
C     BIP is recomputed in routine from BIP_START
      CHUNK   = 1
C     Note: all local variables (BIP, BIP_START,
C           POSELT_BLOCK,...) are private
!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
#endif
C     Loop over all the off diagonal BLR blocks in the panel
      DO IP = BEG_I, END_I
C       With multithreading, to increase BIP correctly for each thread,
C       we need to recompute the BIP depending on the loop index.
C       Warning: in case of multithreading we have no guarantee that
C       each thread processes IP values in increasing order, so we
C       must start the loop below at I=1 to compute BIP.
#if ! defined(BLR_NOOPENMP)
        BIP = BIP_START
        DO I = BEG_I, IP-1 
          BIP  = BIP +  BLR_PANEL(I-CURRENT_BLR)%M
        ENDDO
#endif
        IF (DIR .eq. 'V') THEN
           IF (BIP .LE. LDA21) THEN
C            Leading dimension LDA11
             POSELT_BLOCK = POSELT + int(LDA11,8)*int(BIP-1,8) + 
     &              int(BEGS_BLR_DIAG - 1,8)
           ELSE
C            Leading dimension LDA11 until NPIV, LDA21 afterwards
C            Move to first position after pivot block
             POSELT_BLOCK = POSELT + int(LDA11,8)*int(LDA21,8)+
     &              int(BEGS_BLR_DIAG - 1,8)
C            Shift with (BIP-1-NPIV)
             POSELT_BLOCK = POSELT_BLOCK +
     &                      int(LDA21,8)*int(BIP-1-LDA21,8)
             LD_BLK_IN_FRONT=int(LDA21,8)
           ENDIF
        ELSE ! DIR = 'H'
         POSELT_BLOCK = POSELT + int(LDA11,8)*int(BEGS_BLR_DIAG-1,8)
     &              + int(BIP-1,8)
        ENDIF
        M = BLR_PANEL(IP-CURRENT_BLR)%M
        N = BLR_PANEL(IP-CURRENT_BLR)%N
        IF(present(ONLY_NELIM_IN)) THEN
          ONLY_NELIM = ONLY_NELIM_IN
        ELSE
          ONLY_NELIM = N
        ENDIF
        K = BLR_PANEL(IP-CURRENT_BLR)%K
        IF (BLR_PANEL(IP-CURRENT_BLR)%ISLR) THEN
          IF (K.EQ.0) THEN
            IF (DIR .eq. 'V') THEN
              DO I = 1, M
                IF (BIP+I-1.GT.LDA21) THEN
C                 this part of the block is > NPIV, leading dimension of A 
C                 should change
                  LD_BLK_IN_FRONT = int(LDA21,8)
                ENDIF
                A(POSELT_BLOCK+int(I-1,8)*LD_BLK_IN_FRONT :
     &                   POSELT_BLOCK+int(I-1,8)*LD_BLK_IN_FRONT
     &                          + int(N-1,8)) = ZERO
              ENDDO
            ELSE ! DIR = 'H'
              DO I = N-ONLY_NELIM+1, N
              A(POSELT_BLOCK+int(I-1,8)*int(LDA11,8):
     &          POSELT_BLOCK+int(I-1,8)*int(LDA11,8) + int(M-1,8)) 
     &              = ZERO
              ENDDO
            ENDIF
            GOTO 1800
          ENDIF
          IF (DIR .eq. 'V') THEN
C            LD_BLK_IN_FRONT: LDA11 or LDA21
            IF (DIR .eq.'V' .AND. BIP .LE. LDA21
     &                 .AND. BIP + M - 1 .GT. LDA21) THEN
C             Block a cheval          
              CALL sgemm('T', 'T', N, LDA21-BIP+1, K, ONE ,
     &            BLR_PANEL(IP-CURRENT_BLR)%R(1,1) , K, 
     &            BLR_PANEL(IP-CURRENT_BLR)%Q(1,1) , M, 
     &            ZERO, A(POSELT_BLOCK), int(LD_BLK_IN_FRONT))
              CALL sgemm('T', 'T', N, BIP+M-LDA21-1, K, ONE ,
     &            BLR_PANEL(IP-CURRENT_BLR)%R(1,1) , K, 
     &            BLR_PANEL(IP-CURRENT_BLR)%Q(LDA21-BIP+2,1) , M, 
     &            ZERO, A(POSELT_BLOCK+int(LDA21-BIP,8)*int(LDA11,8)),
     &            LDA21)
            ELSE
              CALL sgemm('T', 'T', N, M, K, ONE ,
     &            BLR_PANEL(IP-CURRENT_BLR)%R(1,1) , K, 
     &            BLR_PANEL(IP-CURRENT_BLR)%Q(1,1) , M, 
     &            ZERO, A(POSELT_BLOCK), int(LD_BLK_IN_FRONT))
            ENDIF
          ELSE ! DIR = 'H'
             CALL sgemm('N', 'N', M, ONLY_NELIM, K, ONE,
     &          BLR_PANEL(IP-CURRENT_BLR)%Q(1,1), M, 
     &          BLR_PANEL(IP-CURRENT_BLR)%R(1,N-ONLY_NELIM+1), K, ZERO,
     &          A(POSELT_BLOCK+int(N-ONLY_NELIM,8)*int(LDA11,8)), LDA11)
          ENDIF
          PROMOTE_COST = 2.0D0*M*K*ONLY_NELIM
          IF(present(ONLY_NELIM_IN)) THEN
            CALL UPD_FLOP_DECOMPRESS(PROMOTE_COST, .FALSE.)
          ENDIF
        ELSE  IF (COPY_DENSE_BLOCKS) THEN
          IF (DIR .eq. 'V') THEN
            DO I = 1, M
              IF (BIP+I-1.GT.LDA21) THEN
C               this part of the block is > NPIV, leading dimension of A 
C               should change
                LD_BLK_IN_FRONT = int(LDA21,8)
              ENDIF
                A(POSELT_BLOCK+int(I-1,8)*LD_BLK_IN_FRONT :
     &            POSELT_BLOCK+int(I-1,8)*LD_BLK_IN_FRONT
     &                                           + int(N-1,8)) 
     &          = BLR_PANEL(IP-CURRENT_BLR)%Q(I,1:N)
            ENDDO
          ELSE ! DIR = 'H'
            DO I = N-ONLY_NELIM+1, N
              A(POSELT_BLOCK+int(I-1,8)*int(LDA11,8):
     &           POSELT_BLOCK+int(I-1,8)*int(LDA11,8) + int(M-1,8))
     &        = BLR_PANEL(IP-CURRENT_BLR)%Q(1:M,I)
            ENDDO
          ENDIF
        ENDIF
 1800   CONTINUE
#if defined(BLR_NOOPENMP)
        BIP  = BIP +  BLR_PANEL(IP-CURRENT_BLR)%M
#endif
      ENDDO
#if ! defined(BLR_NOOPENMP)
!$OMP END DO
#endif
      END SUBROUTINE SMUMPS_DECOMPRESS_PANEL
      SUBROUTINE SMUMPS_COMPRESS_CB(A, LA, POSELT, LDA,
C
     &        BEGS_BLR, BEGS_BLR_U, NB_ROWS, NB_COLS, NB_INASM,
     &        NROWS, NCOLS, INODE,   
     &        IWHANDLER, SYM, NIV, IFLAG, IERROR,
     &        TOLEPS, TOL_OPT, KPERCENT, K489, CB_LRB,
     &        WORK, TAU, JPVT, LWORK, RWORK, BLOCK,
     &        MAXI_CLUSTER, KEEP8,
C //PIV :
     &        NFS4FATHER, NPIV, NVSCHUR_K253, KEEP,
C             Optional
     &        M_ARRAY,
C Optional parameter Type 1 nodes or LDLT T2 nodes
     &        NELIM, 
C Optional parameter for LDLT T2 nodes
     &        NBROWSinF
     &        )
!$    USE OMP_LIB
C    This subroutine compresses the CB block.
C
C    In case of LDLT Type 2 matrices and KEEP(219).NE.0
C    then the maximum of the absolute values
C    per column is computed prior to compression and 
C    is saved in BLR structure to be sent together with 
C    the compressed CB.
C    The max is computed only on columns corresponding
C    to the fully summed columns of the father 
C    NFS4FATHER that should include at this point
C    all the non eliminated columns of the son son
C    that will be sent to the father.
C    (see SMUMPS_COMPUTE_MAXPERCOL to compute the max)
C    For LDLT T2 node, NBROWSinF is needed to skip
C    rows mapped on father
C    
C
C    NROWS: note that for Level 1 nodes, NRROWS
C           does not include the non-eliminated variables of the node.
C           (thus for T1 nodes: NFRONT=NROWS+NPIV+NELIM)
C     -----------
C     Parameters
C     -----------
      INTEGER(8), intent(in)       :: LA
      REAL, intent(inout)       :: A(LA)
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: LDA, NB_ROWS, NB_COLS, NB_INASM
      INTEGER, INTENT(IN)          :: NIV, IWHANDLER, MAXI_CLUSTER, 
     &                                KPERCENT, TOL_OPT, LWORK
      INTEGER, INTENT(IN)          :: K489, NROWS, NCOLS, INODE, SYM
      INTEGER, intent(inout)         :: IFLAG, IERROR
C     TYPE(LRB_TYPE), TARGET, intent(inout) :: CB_LRB(NB_ROWS,NB_COLS)
C     Contiguity info not critical, avoid copy when arg. is pointer
C     by declaring an assumed-shape array
      TYPE(LRB_TYPE), TARGET, intent(inout) :: CB_LRB(:,:)
      INTEGER, DIMENSION(:) :: BEGS_BLR, BEGS_BLR_U
      REAL, TARGET, DIMENSION(:) :: RWORK
      REAL, TARGET, DIMENSION(:,:) :: BLOCK
      REAL, TARGET, DIMENSION(:) :: WORK, TAU
      INTEGER, TARGET, DIMENSION(:) :: JPVT
      INTEGER(8) :: KEEP8(150)
      REAL,intent(in) :: TOLEPS
C //PIV
C     NPIV number of pivots eliminated in front
C     NVSCHUR_K253 not accessed in case of LU, holds
C                  nb of rows that need not be considered
C                  for //PIV (i.e. Schur+fwd in facto var)
      INTEGER, INTENT(in) :: NFS4FATHER, NPIV, NVSCHUR_K253, KEEP(500)
      REAL, OPTIONAL :: M_ARRAY(max(NFS4FATHER,1))
      INTEGER, intent(in), OPTIONAL :: NELIM
      INTEGER, intent(in), OPTIONAL :: NBROWSinF
C     -----------
C     Local variables
C     -----------
      INTEGER :: M, N, INFO
      INTEGER :: I, J, IBIS, IBIS_END, RANK, MAXRANK, II, JJ
      INTEGER(8) :: POSELT_BLOCK
      LOGICAL :: ISLR
      TYPE(LRB_TYPE), POINTER :: LRB
      INTEGER :: OMP_NUM
C //PIV
      INTEGER(8) :: POSA, ASIZE
      INTEGER    :: NROWS_CM
#if ! defined(BLR_NOOPENMP)
      INTEGER :: CHUNK
#endif
      INTEGER, POINTER :: K_STOC(:)
      REAL, POINTER, DIMENSION(:) :: RWORK_THR
      REAL, POINTER, DIMENSION(:,:) :: BLOCK_THR
      REAL, POINTER, DIMENSION(:) :: WORK_THR, TAU_THR
      INTEGER, POINTER, DIMENSION(:) :: JPVT_THR
      REAL :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
      PARAMETER (ZERO=0.0E0)
      NULLIFY(K_STOC)
#if ! defined(BLR_NOOPENMP)
!$OMP MASTER
#endif
C    FIXME FWD in facto  ?
C        see SMUMPS_BUF_SEND_CONTRIB_TYPE2
C       KEEP(253) last rows should not be compressed ?
C       why  KEEP253_LOC and should be suppressed 
C       only for last slave ?
C
      IF ( (KEEP(219).NE.0).AND.(KEEP(50).EQ.2).AND.
     &     (NFS4FATHER.GT.0) ) THEN
C     {
C //PIV
C      NROWS_CM = nb of rows on which MAX is to be computed
C      NFS4FATHER includes the NELIM columns
       IF (NIV.EQ.1) THEN
C       NROWS does not include NELIM 
        NROWS_CM  = NROWS  - (NFS4FATHER-NELIM)
       ELSE
        NROWS_CM  = NROWS  - NBROWSinF
       ENDIF
       IF (NROWS_CM-NVSCHUR_K253.GT.0)  THEN
         IF (NIV.EQ.1) THEN
          POSA     = POSELT
     &             + int(LDA,8)*int(NPIV+NFS4FATHER,8)
     &             + int(NPIV,8)
          ASIZE    = int(LDA,8)*int(LDA,8)
     &             - int(LDA,8)*int(NPIV+NFS4FATHER,8)
     &             - int(NPIV,8) 
         ELSE
          POSA     = POSELT 
     &             + int(LDA,8)*int(NBROWSinF,8)
     &             + int(NPIV,8) 
          ASIZE    = int(NROWS,8)*int(LDA,8)
     &             - int(LDA,8)*int(NBROWSinF,8)
     &             - int(NPIV,8)
         ENDIF
C        Compute M_ARRAY
         CALL SMUMPS_COMPUTE_MAXPERCOL (
     &      A(POSA), ASIZE, LDA, 
     &      NROWS_CM-NVSCHUR_K253,
     &      M_ARRAY(1), NFS4FATHER, .FALSE., 
     &      -9999)
       ELSE 
C         Nothing to compute just set M_ARRAY to 0
          DO I=1, NFS4FATHER
            M_ARRAY(I) = ZERO
          ENDDO
       ENDIF
C     }
      ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP END MASTER
!$OMP BARRIER
#endif
C     Loop over all the potential BLR blocks in the CB
      IF (SYM.EQ.0.OR.NIV.EQ.2) THEN
C       Note that in symmetric NIV2, we are allocating the full
C       rectangular band for simplicity but blocks on the upper
C       triangular part at the end of the band are skipped        
        IBIS_END = NB_ROWS*NB_COLS
      ELSE
        IBIS_END = NB_ROWS*(NB_COLS+1)/2
      ENDIF
#if ! defined(BLR_NOOPENMP)
      CHUNK = 1
!$OMP DO SCHEDULE(DYNAMIC,CHUNK) 
!$OMP& PRIVATE(I, J, POSELT_BLOCK, M, N, OMP_NUM, INFO, RANK,
!$OMP&         MAXRANK, ISLR, II, JJ, LRB)
#endif
      DO IBIS = 1,IBIS_END
        IF (IFLAG.LT.0) CYCLE     
        OMP_NUM = 0 ! private
#if ! defined(BLR_NOOPENMP)
!$      OMP_NUM = OMP_GET_THREAD_NUM()
#endif
        BLOCK_THR => BLOCK(1:MAXI_CLUSTER,OMP_NUM*MAXI_CLUSTER+1:
     &                  (OMP_NUM+1)*MAXI_CLUSTER) 
        JPVT_THR  => JPVT(OMP_NUM*MAXI_CLUSTER+1:
     &                 (OMP_NUM+1)*MAXI_CLUSTER) 
        TAU_THR   => TAU(OMP_NUM*MAXI_CLUSTER+1:
     &                 (OMP_NUM+1)*MAXI_CLUSTER) 
        WORK_THR  => WORK(OMP_NUM*LWORK+1:
     &                 (OMP_NUM+1)*LWORK) 
        RWORK_THR => RWORK(OMP_NUM*2*MAXI_CLUSTER+1:
     &                  (OMP_NUM+1)*2*MAXI_CLUSTER) 
C       Determining I,J from IBIS
        IF (SYM.EQ.0.OR.NIV.EQ.2) THEN
          I = (IBIS-1)/NB_COLS+1
          J = IBIS - (I-1)*NB_COLS
        ELSE
          I = CEILING((1.0D0+SQRT(1.0D0+8.0D0*dble(IBIS)))/2.0D0)-1
          J = IBIS - I*(I-1)/2
        ENDIF
        IF (NIV.EQ.1) THEN
          I = I+NB_INASM
          J = J+NB_INASM
        ELSE
C         BEGS_BLR_COL          
          J = J+NB_INASM
          IF (SYM.NE.0) THEN
C           Skip block on upper triangular part      
            IF (BEGS_BLR_U(J).GE.BEGS_BLR(I+2)+NCOLS-NROWS-1+
     &          BEGS_BLR_U(NB_INASM+1)) THEN 
              CYCLE
            ENDIF
          ENDIF
        ENDIF
        IF (NIV.EQ.1) THEN
          M = BEGS_BLR(I+1)-BEGS_BLR(I)
          POSELT_BLOCK = POSELT + int(LDA,8)*int(BEGS_BLR(I)-1,8) + 
     &           int(BEGS_BLR_U(J)-1,8)
          IF (I .EQ. NB_INASM+1 .AND. present(NELIM)) THEN
            POSELT_BLOCK = POSELT_BLOCK + int(NELIM,8)*int(LDA,8)
            M = M - NELIM
          ENDIF
          N = BEGS_BLR_U(J+1)-BEGS_BLR_U(J)
        ELSE
C         BEGS_BLR_{L,LS} are of size NB_BLR_{L,U,LS}+2 rather than +1 in T2
C         slaves. However, BEGS_BLR_COL is of size NB_BLR_COL+1.
          M = BEGS_BLR(I+2)-BEGS_BLR(I+1)
          POSELT_BLOCK = POSELT + int(LDA,8)*int(BEGS_BLR(I+1)-1,8)
     &           + int(BEGS_BLR_U(J)-1,8)
          IF (SYM.EQ.0) THEN
            N = BEGS_BLR_U(J+1)-BEGS_BLR_U(J)
          ELSE
C           In LDLT T2, crop diagonal blocks to updated part            
C           This is a bit tricky: the update has been done following the
C           row clustering but the compression is done following the
C           column one; therefore the diagonal blocks defined by
C           BEGS_BLR_COL may not be entirely updated and need to be
C           cropped according to BEGS_BLR_LS. 
            N = min(BEGS_BLR_U(J+1), BEGS_BLR(I+2) + NCOLS - NROWS -1
     &              + BEGS_BLR_U(NB_INASM+1)) - BEGS_BLR_U(J) 
          ENDIF
        ENDIF
        JPVT_THR(1:MAXI_CLUSTER) = 0
        IF (NIV.EQ.1) THEN
          LRB => CB_LRB(I-NB_INASM,J-NB_INASM)
        ELSE
C         BEGS_BLR_COL          
          LRB => CB_LRB(I,J-NB_INASM)
        ENDIF
        IF (K489.EQ.3) THEN
C           Save CB blocks in the structure but don't compress for now
            MAXRANK = 1
            RANK = MAXRANK+1
            INFO = 0
            ISLR = .FALSE.
            GOTO 3800
        ENDIF
C       Copy current block to a temp block
        DO II=1,M
          BLOCK_THR(II,1:N)=
     &    A( POSELT_BLOCK+int(II-1,8)*int(LDA,8) :
     &    POSELT_BLOCK+int(II-1,8)*int(LDA,8)+int(N-1,8) )
        ENDDO  
C       Compute Householder vectors and R matrix
        MAXRANK = floor(real(M*N)/real(M+N))
C       Equal to MAXRANK we gain in flops
C       with constant memory, further reducing MAXRANK then
C       QR might stop earlier on FR blocks 
C       so that is will be more efficient
C       old: MAXRANK = max(1, MAXRANK/2)
C       Compress factors controlled by KPERCENT 
        MAXRANK = max (1, int((MAXRANK*KPERCENT/100)))
        allocate(K_STOC(NPREC_STOC))
        K_STOC = 0
        CALL SMUMPS_TRUNCATED_RRQR( M, N,
     &       BLOCK_THR(1,1),
     &       MAXI_CLUSTER, JPVT_THR(1), 
     &       TAU_THR(1), 
     &       WORK_THR(1), N, 
     &       RWORK_THR(1), 
     &       TOLEPS, TOL_OPT, RANK, MAXRANK, INFO,
     &       K_STOC(1), (KEEP(452).GT.0), KEEP(453), KEEP(454).GE.1,
     &       ISLR)
 3800 CONTINUE
        IF (INFO < 0) THEN
           WRITE(*,*) " PROBLEM IN ARGUMENT NUMBER ",INFO,
     &                 " OF TRUNCATED_RRQR WHILE COMPRESSING A CB BLOCK"
           write(6,*) " M, N=", M, N
           CALL MUMPS_ABORT()
        END IF
        CALL ALLOC_LRB(LRB, RANK, M, N, ISLR, IFLAG, IERROR, KEEP8)
        IF (IFLAG.LT.0) CYCLE
        IF (ISLR) THEN 
!         block I,J is LR 
          IF (associated(LRB%K_STOC)) THEN
            deallocate(LRB%K_STOC)
            nullify(LRB%K_STOC)
          ENDIF
          if (associated(K_STOC)) then
              LRB%K_STOC=>K_STOC
          end if
          nullify(K_STOC)
           IF (RANK .GT. 0) THEN 
C compute R from the upper triangular original R computed
C by mygeqp3: copy r into the blr block data structure column-by-column
C in the order given by jpvt
               DO JJ=1,N
                  DO II=1,MIN(RANK,JJ)
                     LRB%R(II,JPVT_THR(JJ)) = BLOCK_THR(II,JJ)
                  ENDDO 
                  IF(JJ.LT.RANK) LRB%R(MIN(RANK,JJ)+1:RANK,JPVT_THR(JJ))
     &                 = ZERO
               ENDDO
C compute Q from H and store it in %Q field
               CALL sorgqr 
     &           (M, RANK, RANK,
     &           BLOCK_THR(1,1), 
     &           MAXI_CLUSTER, TAU_THR(1), 
     &           WORK_THR(1), LWORK, INFO )
               DO II=1,RANK
                 DO JJ= 1, M
                  LRB%Q(JJ,II) = BLOCK_THR(JJ,II)
                 ENDDO
               END DO
               IF (INFO < 0) THEN
                 WRITE(*,*) " PROBLEM IN ARGUMENT NUMBER ",INFO,
     &                     " OF CUNGQR WHILE COMPRESSING A CB BLOCK"
                 CALL MUMPS_ABORT()
               END IF
               IF(KEEP(452).GT.0) THEN
                  CALL SMUMPS_LRMIXED_CONVERT_LRB(LRB, KEEP(34), KEEP8, 
     &                                    IFLAG, IERROR)
                  IF (IFLAG .LT. 0) CYCLE
               ENDIF
               IF (K489.NE.3) THEN
                 CALL UPD_FLOP_COMPRESS(LRB%M, LRB%N, LRB%K, 
     &                    LRB%ISLR, CB_COMPRESS=.TRUE.)
               ENDIF
          END IF
          CALL UPD_MRY_CB_LRGAIN(LRB%M, LRB%N, LRB%K
     &                  , (KEEP(452).GT.0), LRB%K_STOC, 
     &                    NPREC_STOC, COSTS_PREC
     &                           )
        ELSE !block I,J is not LR
          IF (associated(K_STOC)) THEN
            LRB%K_STOC => K_STOC
          ENDIF
          nullify(K_STOC)
C          Copy current block to a temp block
           DO II=1,M
             LRB%Q(II,1:N) =
     &       A( POSELT_BLOCK+int((II-1),8)*int(LDA,8) :
     &         POSELT_BLOCK+int((II-1),8)*int(LDA,8)
     &                     +int(N-1,8) )
           END DO  
           IF (K489.NE.3) THEN
               CALL UPD_FLOP_COMPRESS(LRB%M, LRB%N, LRB%K, 
     &                    LRB%ISLR, CB_COMPRESS=.TRUE.)
           ENDIF
C          back to zero to ease operations
           LRB%K = -1
        END IF
      END DO
#if ! defined(BLR_NOOPENMP)
!$OMP END DO 
#endif
#if ! defined(BLR_NOOPENMP)
!$OMP MASTER
#endif
      CALL UPD_MRY_CB_FR(NROWS, NCOLS, SYM)
#if ! defined(BLR_NOOPENMP)
!$OMP END MASTER
#endif
      END SUBROUTINE SMUMPS_COMPRESS_CB
      SUBROUTINE SMUMPS_COMPRESS_PANEL(
     &        A, LA, POSELT, IFLAG, IERROR, NFRONT,
     &        BEGS_BLR, NB_BLR, TOLEPS, TOL_OPT, K458, K473, BLR_PANEL,
     &        CURRENT_BLR,
     &        DIR, WORK, TAU, JPVT, 
     &        LWORK, RWORK, BLOCK,
     &        MAXI_CLUSTER, NELIM, 
     &        LBANDSLAVE, NPIV, ISHIFT, NIV, KPERCENT, 
     &        KEEP8,
     &        LRMIXED_ACTIVATED, K453, MIXED_ADMISSIBILITY,
     &        K480, BEG_I_IN, END_I_IN, FRSWAP
C     &        ,DIAG_NORM
     &        )
!$    USE OMP_LIB
C     -----------
C     Parameters
C     -----------
      INTEGER(8), intent(in)       :: LA
      REAL, intent(inout)       :: A(LA)
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_BLR, CURRENT_BLR, NIV
      INTEGER, intent(inout)          :: IFLAG, IERROR
c     TYPE(LRB_TYPE), intent(inout) :: BLR_PANEL(NB_BLR-CURRENT_BLR)
C     Contiguity information not so much needed=>assumed shape array
C     to allow calling it with a pointer / assumed-shape argument
      TYPE(LRB_TYPE), TARGET, intent(inout) :: BLR_PANEL(:)
      REAL, TARGET, DIMENSION(:) :: RWORK
      REAL, TARGET, DIMENSION(:,:) :: BLOCK
      REAL, TARGET, DIMENSION(:) :: WORK, TAU
      INTEGER, TARGET, DIMENSION(:) :: JPVT
      INTEGER :: BEGS_BLR(:)
      INTEGER(8) :: KEEP8(150)
      LOGICAL, intent(in)          :: LRMIXED_ACTIVATED
      INTEGER, intent(in)          :: K453
      LOGICAL, intent(in)          :: MIXED_ADMISSIBILITY
      INTEGER, OPTIONAL, intent(in) :: K480
      INTEGER,OPTIONAL,intent(in) :: BEG_I_IN, END_I_IN
      LOGICAL, OPTIONAL, intent(in) :: FRSWAP
C  ISHIFT used to shift columns of L band of a slave (only used
C         when  LBANDSLAVE is true)
C  it should stored in BEGS_BLR_U but since the ISHIFT first columns
C  of U are not stored, it need a special processing
      INTEGER, intent(in)          :: NPIV, ISHIFT, KPERCENT,
     &                                K458, K473, TOL_OPT
      LOGICAL, intent(in)          :: LBANDSLAVE
      INTEGER                      :: MAXI_CLUSTER, LWORK, NELIM
      REAL,intent(in)              :: TOLEPS
C     V = Vertical panel, H = Horizontal panel
      CHARACTER(len=1) :: DIR
C     -----------
C     Local variables
C     -----------
      INTEGER :: IP, NB_BLOCKS_PANEL, M, N, RANK, MAXRANK
      INTEGER :: INFO, I, J, K, IS, BEG_I, END_I
      INTEGER(8) :: POSELT_BLOCK
      LOGICAL :: ISLR
      REAL :: ONE, ALPHA, ZERO
      PARAMETER (ONE = 1.0E0, ALPHA=-1.0E0)
      PARAMETER (ZERO = 0.0D0)
      TYPE(LRB_TYPE), POINTER :: LRB
      INTEGER :: OMP_NUM
      REAL, POINTER, DIMENSION(:) :: RWORK_THR
      REAL, POINTER, DIMENSION(:,:) :: BLOCK_THR
      REAL, POINTER, DIMENSION(:) :: WORK_THR, TAU_THR
      INTEGER, POINTER, DIMENSION(:) :: JPVT_THR
#if ! defined(BLR_NOOPENMP)
      INTEGER :: CHUNK
#endif
      INTEGER, POINTER :: K_STOC(:)
      NULLIFY(K_STOC)
      IF(present(BEG_I_IN)) THEN
        BEG_I = BEG_I_IN
      ELSE
        BEG_I = CURRENT_BLR+1
      ENDIF
      IF(present(END_I_IN)) THEN
        END_I = END_I_IN
      ELSE
        END_I = NB_BLR
      ENDIF
      IF (LBANDSLAVE) THEN
       IS = ISHIFT
      ELSE
       IS=0
      ENDIF
      IF (DIR .eq. 'V') THEN
         IF (LBANDSLAVE) THEN
          N = NPIV
         ELSE
          N = BEGS_BLR(CURRENT_BLR+1)-BEGS_BLR(CURRENT_BLR)-NELIM
         ENDIF
      ELSE IF (DIR .eq. 'H') THEN
        N = BEGS_BLR(CURRENT_BLR+1)-BEGS_BLR(CURRENT_BLR)-NELIM
      ELSE
         WRITE(*,*) " WRONG ARGUMENT IN SMUMPS_COMPRESS_PANEL "
         CALL MUMPS_ABORT()
      END IF
        IF (N.LE.0) THEN
         RETURN
        ENDIF
      NB_BLOCKS_PANEL = NB_BLR-CURRENT_BLR
C     Loop over all the BLR blocks in the panel      
      OMP_NUM = 0 ! in case of compilation with -DBLR_NOOPENMP
#if ! defined(BLR_NOOPENMP)
      CHUNK = 1
!$OMP DO PRIVATE(INFO, POSELT_BLOCK, RANK, MAXRANK, I, J, OMP_NUM, LRB)
!$OMP&   SCHEDULE(DYNAMIC,CHUNK)
#endif
      DO IP = BEG_I, END_I
C     Compute a rank-revealing QR on IP-th block.
C     This is done in-place (the full-rank version of the block is not
C     needed anymore => th promoted version will replace it later. 
        IF (IFLAG.LT.0) CYCLE
#if ! defined(BLR_NOOPENMP)
        OMP_NUM = 0 ! private
!$      OMP_NUM = OMP_GET_THREAD_NUM()
#endif
        BLOCK_THR => BLOCK(1:MAXI_CLUSTER,OMP_NUM*MAXI_CLUSTER+1:
     &                  (OMP_NUM+1)*MAXI_CLUSTER) 
        JPVT_THR  => JPVT(OMP_NUM*MAXI_CLUSTER+1:
     &                 (OMP_NUM+1)*MAXI_CLUSTER) 
        TAU_THR   => TAU(OMP_NUM*MAXI_CLUSTER+1:
     &                 (OMP_NUM+1)*MAXI_CLUSTER) 
        WORK_THR  => WORK(OMP_NUM*LWORK+1:
     &                 (OMP_NUM+1)*LWORK) 
        RWORK_THR => RWORK(OMP_NUM*2*MAXI_CLUSTER+1:
     &                  (OMP_NUM+1)*2*MAXI_CLUSTER) 
        LRB => BLR_PANEL(IP-CURRENT_BLR)
        RANK = 0
        M = BEGS_BLR(IP+1)-BEGS_BLR(IP)
        IF (N.EQ.0.OR.M.EQ.0) THEN
         CYCLE
        ENDIF
        IF (DIR .eq. 'V') THEN
          POSELT_BLOCK = POSELT + 
     &              int(NFRONT,8) * int(BEGS_BLR(IP)-1,8) + 
     &              int(BEGS_BLR(CURRENT_BLR) + IS - 1,8)
        ELSE ! DIR = H
          POSELT_BLOCK = POSELT +
     &              int(NFRONT,8)*int(BEGS_BLR(CURRENT_BLR)-1,8) + 
     &              int( BEGS_BLR(IP) - 1,8)
        ENDIF
        IF (present(K480)) then 
        IF (K480.GE.5) THEN
C         In LUA KEEP(480)>=5 strategy, we already have some blocks in
C         LR form (coming from the recompressed accumulators that we
C         never decompressed).            
          IF (LRB%ISLR) THEN
C           This block is already in LR form, skip it            
             IF (M.NE.LRB%M) THEN
              write(*,*) 'Internal error in SMUMPS_COMPRESS_PANEL',
     &                    ' M size inconsistency',M,
     &                    LRB%M
              CALL MUMPS_ABORT()
            ENDIF
            IF (N.NE.LRB%N) THEN
              write(*,*) 'Internal error in SMUMPS_COMPRESS_PANEL',
     &                    ' N size inconsistency',N,
     &                    LRB%N
              CALL MUMPS_ABORT()
            ENDIF
            MAXRANK = floor(real(M*N)/real(M+N))
            IF (LRB%K.GT.MAXRANK) THEN
              write(*,*) 'Internal error in SMUMPS_COMPRESS_PANEL',
     &                    ' MAXRANK inconsistency',MAXRANK,
     &                    LRB%K
              CALL MUMPS_ABORT()
            ENDIF
C            write(*,*) 'SKIPPING compress',IP-CURRENT_BLR,
C     &      BLR_PANEL(IP-CURRENT_BLR)%K
            GOTO 3000
          ENDIF
        ENDIF
        ENDIF
        JPVT_THR(1:MAXI_CLUSTER) = 0
        IF (K473.EQ.1 .OR. IP .LT. BEG_I+K458) THEN
C           BFR facto, skip the compress            
            MAXRANK = 1
            RANK = MAXRANK+1
            INFO = 0
            ISLR = .FALSE.
            GOTO 3800
        ENDIF
C     Copy current block to a temp block
        IF (DIR .eq. 'V') THEN
            DO I=1,M
C               FIXME: BLOCK_THR is contiguous but this is not exploited
C               FIXME: access BLOCK_THR with wrong stride
                BLOCK_THR(I,1:N)=
     &          A( POSELT_BLOCK+int(I-1,8)*int(NFRONT,8) :
     &          POSELT_BLOCK+int(I-1,8)*int(NFRONT,8)+int(N-1,8) )
            END DO  
        ELSE ! DIR = H
            DO I=1,N
C               FIXME: BLOCK_THR is contiguous but this is not exploited
                BLOCK_THR(1:M,I)=
     &          A( POSELT_BLOCK+int(I-1,8)*int(NFRONT,8) :
     &          POSELT_BLOCK+int(I-1,8)*int(NFRONT,8)+int(M-1,8) )
            END DO  
        END IF
C     Compute Householder vectors and R matrix
        MAXRANK = floor(real(M*N)/real(M+N))
C     Equal to MAXRANK we gain in flops
C     with constant memory, further reducing MAXRANK then
C     QR might stop earlier on FR blocks 
C     so that is will be more efficient
C     old: MAXRANK = max(1, MAXRANK/2)
C     Compress factors controlled by KPERCENT 
        MAXRANK = max (1, int((MAXRANK*KPERCENT/100)))
        allocate(K_STOC(NPREC_STOC))
        K_STOC = 0 
C       Implicit interface but pointers
C       known to point to contiguous data
        CALL SMUMPS_TRUNCATED_RRQR( M, N,
     &       BLOCK_THR(1,1),
     &       MAXI_CLUSTER, JPVT_THR(1), 
     &       TAU_THR(1), 
     &       WORK_THR(1), N, 
     &       RWORK_THR(1), 
     &       TOLEPS, TOL_OPT, RANK, MAXRANK, INFO,
     &       K_STOC(1), LRMIXED_ACTIVATED, K453, MIXED_ADMISSIBILITY,
     &       ISLR)
 3800 CONTINUE
        IF (INFO < 0) THEN
           WRITE(*,*) " PROBLEM IN ARGUMENT NUMBER ",INFO,
     &                 " OF TRUNCATED_RRQR WHILE COMPRESSING A BLOCK "
           CALL MUMPS_ABORT()
        END IF
        CALL ALLOC_LRB(LRB, RANK,
     &                 M, N, ISLR, IFLAG, IERROR, KEEP8)
        IF (IFLAG.LT.0) CYCLE
        IF ((M.EQ.0).OR.(N.EQ.0)) THEN
          IF (associated(K_STOC)) THEN
            DEALLOCATE(K_STOC)
            NULLIFY(K_STOC)
          ENDIF
          GOTO 3000
        ENDIF
        IF (ISLR) THEN 
!         block IP is LR 
          IF (associated(LRB%K_STOC)) THEN
          ENDIF
          IF (associated(K_STOC)) THEN
            LRB%K_STOC=>K_STOC
          ENDIF
          nullify(K_STOC)
           IF (RANK .EQ. 0) THEN !RANK = 0
           ELSE !LOW RANK BUT NOT ZERO
C compute R from the upper triangular original R computed
C by mygeqp3: copy r into the blr block data structure column-by-column
C in the order given by jpvt
               DO J=1,N
                 DO K=1,min(RANK,J)
                   LRB%R(K, JPVT_THR(J)) = BLOCK_THR(K,J)
                 ENDDO
                 IF(J.LT.RANK) THEN
                   LRB%R(J+1:RANK,JPVT_THR(J)) = ZERO
                 ENDIF
               ENDDO
C compute Q from the Householder vectors and store it in %Q field
               CALL sorgqr 
     &           (M, RANK, RANK,
     &           BLOCK_THR(1,1), 
     &           MAXI_CLUSTER, TAU_THR(1), 
     &           WORK_THR(1), LWORK, INFO )
               DO I=1,RANK
                 DO K=1,M
                   LRB%Q(K,I) = BLOCK_THR(K,I)
                 ENDDO
               END DO
               IF (INFO < 0) THEN
                 WRITE(*,*) " PROBLEM IN ARGUMENT NUMBER ",INFO,
     &                     " OF CUNGQR WHILE COMPRESSING A BLOCK "
                 CALL MUMPS_ABORT()
               END IF
               IF (present(FRSWAP)) THEN
                 CALL UPD_FLOP_COMPRESS(LRB%M, LRB%N, LRB%K, LRB%ISLR, 
     &                                  FRSWAP=FRSWAP)
               ELSE
                 CALL UPD_FLOP_COMPRESS(LRB%M, LRB%N, LRB%K, LRB%ISLR)
               ENDIF
          END IF
        ELSE !A is not LR
#if defined(NOMIXHPC)
          IF (associated(K_STOC)) THEN
            deallocate(K_STOC)
            nullify(K_STOC)
          ENDIF
#else
          IF (associated(LRB%K_STOC)) THEN
            deallocate(LRB%K_STOC)
            nullify(LRB%K_STOC)
          ENDIF
          IF (associated(K_STOC)) THEN
            LRB%K_STOC => K_STOC
          ENDIF
          nullify(K_STOC)
#endif
C          Copy current block to a temp block
           IF (DIR .eq. 'V') THEN
               DO I=1,M
                   LRB%Q(I,1:N) =
     &             A( POSELT_BLOCK+int((I-1),8)*int(NFRONT,8) :
     &               POSELT_BLOCK+int((I-1),8)*int(NFRONT,8)
     &                           +int(N-1,8) )
               END DO  
           ELSE ! DIR = H
               DO I=1,N
                   LRB%Q(1:M,I) =
     &             A( POSELT_BLOCK+int((I-1),8)*int(NFRONT,8) :
     &               POSELT_BLOCK+int((I-1),8)*int(NFRONT,8)
     &                           +int(M-1,8) )
               END DO  
           END IF
           IF (K473.EQ.0) THEN
             IF (present(FRSWAP)) THEN
               CALL UPD_FLOP_COMPRESS(LRB%M, LRB%N, LRB%K, LRB%ISLR, 
     &                                  FRSWAP=FRSWAP)
             ELSE
               CALL UPD_FLOP_COMPRESS(LRB%M, LRB%N, LRB%K, LRB%ISLR)
             ENDIF
           ENDIF
C back to zero to ease operations
           LRB%K = -1
        END IF
 3000   CONTINUE
      END DO ! end loop on all the blr block of current panel
#if ! defined(BLR_NOOPENMP)
C Adding NOWAIT to have better looking traces
C and most importantly to remove synchro between L and U      
!$OMP END DO NOWAIT
#endif
      RETURN
      END SUBROUTINE SMUMPS_COMPRESS_PANEL
      SUBROUTINE SMUMPS_BLR_PANEL_LRTRSM(
     &                A,
     &                LA, POSELT, NFRONT,
     &                IBEG_BLOCK, NB_BLR,
     &                BLR_LorU,
     &                CURRENT_BLR, FIRST_BLOCK, LAST_BLOCK,
     &                NIV, SYM, LorU, LBANDSLAVE,
     &                IW, OFFSET_IW, NASS) 
!$    USE OMP_LIB
C     -----------
C     Parameters
C     -----------
      INTEGER(8), intent(in)  :: LA
      INTEGER, intent(in)     :: NFRONT, NB_BLR, CURRENT_BLR,
     &                           NIV, SYM, LorU
      LOGICAL, intent(in)     :: LBANDSLAVE
      INTEGER(8), intent(in)  :: POSELT 
      INTEGER, intent(in)     :: IBEG_BLOCK, FIRST_BLOCK, LAST_BLOCK
      INTEGER, OPTIONAL, intent(in)     :: NASS
      REAL, intent(inout)  :: A(LA)
C     TYPE(LRB_TYPE), intent(inout)   :: BLR_LorU(NB_BLR-CURRENT_BLR)
C     Contiguity information not critical and pointer passed
C     => use assumed shape array declaration
      TYPE(LRB_TYPE), intent(inout)   :: BLR_LorU(:)
      INTEGER, OPTIONAL :: OFFSET_IW
      INTEGER, OPTIONAL :: IW(*)
C     -----------
C     Local variables
C     -----------
      INTEGER(8) :: POSELT_LOCAL
      INTEGER    :: IP, LDA
#if ! defined(BLR_NOOPENMP)
      INTEGER :: CHUNK
#endif
      REAL :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
      PARAMETER (ZERO=0.0E0)
      LDA = NFRONT
C     Size of the block really eliminated
      IF (LorU.EQ.0.AND.SYM.NE.0.AND.NIV.EQ.2
     &             .AND.(.NOT.LBANDSLAVE)) THEN
        IF (present(NASS)) THEN
          LDA = NASS
       ELSE
          write(*,*) 'Internal error in SMUMPS_BLR_PANEL_LRTRSM'
          CALL MUMPS_ABORT()
        ENDIF
      ENDIF
      IF (LBANDSLAVE) THEN
        POSELT_LOCAL = POSELT
      ELSE
        POSELT_LOCAL = POSELT + 
     &      int(IBEG_BLOCK-1,8)*int(LDA,8) + int(IBEG_BLOCK - 1,8)
      ENDIF
C     Loop on BLR blocks to apply TRSM  
#if ! defined(BLR_NOOPENMP)
      CHUNK = 1
!$OMP DO 
!$OMP& SCHEDULE(DYNAMIC,CHUNK)
#endif  
      DO IP = FIRST_BLOCK, LAST_BLOCK
        CALL SMUMPS_LRTRSM(A, LA, POSELT_LOCAL, NFRONT, LDA,
     &             BLR_LorU(IP-CURRENT_BLR), NIV, SYM, LorU,
     &             IW, OFFSET_IW) 
      END DO
#if ! defined(BLR_NOOPENMP)
C Adding NOWAIT to have better looking traces      
!$OMP END DO NOWAIT
#endif          
      END SUBROUTINE SMUMPS_BLR_PANEL_LRTRSM
      END MODULE SMUMPS_FAC_LR
