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
      SUBROUTINE MUMPS_AB_FREE_LMAT ( LMAT, K147 )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      IMPLICIT NONE
      TYPE(LMATRIX_T)     :: LMAT
      INTEGER, INTENT(IN) :: K147
      INTEGER :: J
      IF (associated(LMAT%COL)) THEN
        IF (LMAT%NBCOL_LOC.GT.0) THEN
         DO J = 1,LMAT%NBCOL_LOC, K147
           IF (associated(LMAT%COL(J)%IRN)) THEN
             DEALLOCATE(LMAT%COL(J)%IRN)
             NULLIFY(LMAT%COL(J)%IRN)
           ENDIF
         ENDDO
        ENDIF
        DEALLOCATE(LMAT%COL)
        NULLIFY(LMAT%COL)
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_AB_FREE_LMAT
      SUBROUTINE MUMPS_AB_FREE_GCOMP ( GCOMP, MEMCNT )
      USE MUMPS_ANA_BLK_M, ONLY : COMPACT_GRAPH_T
      IMPLICIT NONE
      TYPE(COMPACT_GRAPH_T)      :: GCOMP
      INTEGER(8), INTENT(INOUT)  :: MEMCNT
      IF (associated(GCOMP%IPE)) 
     &        MEMCNT = MEMCNT - int(size(GCOMP%IPE),8)
      IF (associated(GCOMP%ADJ)) 
     &        MEMCNT = MEMCNT - GCOMP%SIZEADJALLOCATED - 5
      IF (associated(GCOMP%IPE)) THEN
          DEALLOCATE(GCOMP%IPE)
          NULLIFY(GCOMP%IPE)
      ENDIF
      IF (associated(GCOMP%ADJ)) THEN
          DEALLOCATE(GCOMP%ADJ)
          NULLIFY(GCOMP%ADJ)
      ENDIF
      GCOMP%NZG   = 0_8
      GCOMP%FIRST = 0
      GCOMP%LAST  = 0
      RETURN
      END SUBROUTINE MUMPS_AB_FREE_GCOMP
      SUBROUTINE MUMPS_AB_COMPUTE_SIZEOFBLOCK (
     &          NBLK, NDOF, BLKPTR, BLKVAR,
     &          SIZEOFBLOCKS, DOF2BLOCK )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NBLK, NDOF
      INTEGER, INTENT(IN) :: BLKPTR(NBLK+1), BLKVAR(NDOF)
      INTEGER, INTENT(OUT):: SIZEOFBLOCKS(NBLK), DOF2BLOCK(NDOF)
      INTEGER :: IB, I, IDOF
      DO IB=1, NBLK
        SIZEOFBLOCKS(IB)= BLKPTR(IB+1)-BLKPTR(IB)
        DO I=BLKPTR(IB), BLKPTR(IB+1)-1
          IDOF = BLKVAR(I)
          DOF2BLOCK(IDOF) = IB
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_AB_COMPUTE_SIZEOFBLOCK 
      SUBROUTINE MUMPS_AB_COORD_TO_LMAT (  MYID, 
     &     NBLK, NDOF, NNZ, IRN, JCN,
     &     DOF2BLOCK, 
     &     IFLAG, IERROR, LP, LPOK,
     &     LMAT, OFFDIAG, KEEP
#if ! defined(NOANA_BLKAUTO)
     &     , OPTION, MAPCOL, SIZE_MAPCOL
#endif
     &     )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: MYID, NBLK, NDOF
      INTEGER(8), INTENT(IN) :: NNZ
      INTEGER, INTENT(IN) :: IRN(max(1_8,NNZ)), JCN(max(1_8,NNZ))
      INTEGER, INTENT(IN) :: DOF2BLOCK(NDOF)
      INTEGER             :: LP, IFLAG, IERROR
      LOGICAL, INTENT(IN) :: LPOK
      INTEGER(8), INTENT(OUT) :: OFFDIAG
      INTEGER, INTENT(INOUT) :: KEEP(500)
      TYPE(LMATRIX_T)     :: LMAT
#if ! defined(NOANA_BLKAUTO)
      INTEGER, INTENT(IN)  :: OPTION, SIZE_MAPCOL 
      INTEGER, INTENT(OUT) :: MAPCOL(SIZE_MAPCOL)
#endif
      LOGICAL :: LU_WITH_SYM_STRUCT, NO_DUPPLICATES
      INTEGER, ALLOCATABLE, DIMENSION(:) :: FLAG
      INTEGER :: allocok, IERROR_LOC
      INTEGER :: I, J, JJ, JJB, IIB, IB, JB, PT
      INTEGER(8) :: I8, NB_8, PT_CUR, SIZE_TO_ALLOC_8
      INTEGER, POINTER, DIMENSION(:) :: BLOCK_PTR
      LU_WITH_SYM_STRUCT = (KEEP(202).EQ.1)
      NO_DUPPLICATES     = (KEEP(203).EQ.1)
      LMAT%NBCOL = NBLK
      LMAT%NZL   = 0_8
      LMAT%NBCOL_LOC = NBLK
      LMAT%FIRST     = 1
      LMAT%CONTIGUOUS_COL = .FALSE.
      ALLOCATE(LMAT%COL(NBLK),
     &         FLAG(NBLK), STAT=allocok)
      IF (allocok.NE.0) THEN
           IFLAG = -7
           CALL MUMPS_SET_IERROR(2_8*int(NBLK,8), IERROR)
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate of LMAT%COL"
           END IF
           RETURN
      ENDIF
      DO IB=1,LMAT%NBCOL_LOC
       LMAT%COL(IB)%NBINCOL = 0
      ENDDO
      DO IB=1,LMAT%NBCOL
       FLAG(IB)             = 0
      ENDDO
#if ! defined(NOANA_BLKAUTO)
      IF (OPTION.EQ.3) THEN 
        DO I=1, SIZE_MAPCOL
           MAPCOL(I) = 0
        ENDDO
      ENDIF
#endif
      IERROR_LOC = 0 
      OFFDIAG    = 0
      DO I8=1, NNZ
         I = IRN(I8)
         J = JCN(I8)
         IF ( (I.GT.NDOF).OR.(J.GT.NDOF).OR.(I.LT.1)
     &                     .OR.(J.LT.1)) THEN
           IERROR_LOC = IERROR_LOC + 1
         ELSE
          IB  = DOF2BLOCK(I)
          JB  = DOF2BLOCK(J)
          JJB = min(IB,JB)
#if ! defined(NOANA_BLKAUTO)
          IF (OPTION.EQ.3.AND.IB.EQ.JB) MAPCOL(IB)=MYID+1
#endif
          IF (IB.NE.JB) THEN
           OFFDIAG = OFFDIAG + 1
           IF (LU_WITH_SYM_STRUCT.AND.IB.LT.JB) CYCLE
           LMAT%NZL = LMAT%NZL+1_8
           LMAT%COL(JJB)%NBINCOL =  LMAT%COL(JJB)%NBINCOL + 1
           ENDIF
         ENDIF
      ENDDO
      IF (IERROR_LOC.GE.1) THEN
         IF (mod(IFLAG,2) .EQ. 0) THEN 
           IFLAG = IFLAG+1
           IERROR = IERROR_LOC
         ENDIF
      ENDIF
      DO JB=1,LMAT%NBCOL_LOC, KEEP(147)
       SIZE_TO_ALLOC_8 = 0_8
       DO JJ=JB, min(JB+KEEP(147)-1,LMAT%NBCOL_LOC)
         SIZE_TO_ALLOC_8 = SIZE_TO_ALLOC_8 + 
     &                     int(LMAT%COL(JJ)%NBINCOL, 8)
       ENDDO
       ALLOCATE(LMAT%COL(JB)%IRN(max(SIZE_TO_ALLOC_8,1_8)),
     &          STAT=allocok)
       IF (allocok.NE.0) THEN
          IFLAG  = -7
          CALL MUMPS_SET_IERROR(max(SIZE_TO_ALLOC_8,1_8), IERROR)
          IF ( LPOK ) THEN
             WRITE(LP, *) " ERROR allocate block of ", 
     &               KEEP(147), " columns in LMAT%COL",
     &               " of size ", max(SIZE_TO_ALLOC_8, 1_8)
          END IF
          RETURN
       ENDIF
       IF (KEEP(147).GT.1) THEN
         BLOCK_PTR => LMAT%COL(JB)%IRN
         NB_8    = int(LMAT%COL(JB)%NBINCOL, 8) 
         PT_CUR = 1_8 + NB_8 
         DO JJ=JB+1, min(JB+KEEP(147)-1,LMAT%NBCOL_LOC)
           NB_8 =  int(LMAT%COL(JJ)%NBINCOL,8)
           IF (NB_8.GT.0) THEN
            LMAT%COL(JJ)%IRN => BLOCK_PTR(PT_CUR: PT_CUR+NB_8-1_8)
            PT_CUR = PT_CUR + NB_8
           ENDIF
         ENDDO
       ENDIF
      ENDDO
      DO I8=1, NNZ
         I = IRN(I8)
         J = JCN(I8)
         IF ( (I.LE.NDOF).AND.(J.LE.NDOF).AND.(I.GE.1)
     &                     .AND.(J.GE.1)) THEN
          IB  = DOF2BLOCK(I)
          JB  = DOF2BLOCK(J)
          IF (LU_WITH_SYM_STRUCT.AND.IB.LT.JB) CYCLE
          JJB = min(IB,JB)
          IIB = max(IB,JB)
          IF (IIB.NE.JJB) THEN
             PT        = FLAG(JJB)+1
             FLAG(JJB) = PT
             LMAT%COL(JJB)%IRN(PT) = IIB
          ENDIF
         ENDIF
      ENDDO
      IF (  ( NDOF.NE.NBLK ) 
     &     .OR. 
     &      ( KEEP(50).NE.0 .AND. .NOT.NO_DUPPLICATES )
     &     .OR.
     &      ( KEEP(50).EQ.0.AND.
     &       .NOT.(LU_WITH_SYM_STRUCT.AND.NO_DUPPLICATES) ) 
     &     )
     &   THEN
        CALL MUMPS_AB_LOCALCLEAN_LMAT ( MYID,
     &     LMAT, KEEP(147), 
     &     FLAG(1), LMAT%NBCOL, IFLAG, IERROR, LP, LPOK
     &   )
      ENDIF
      IF (allocated(FLAG)) DEALLOCATE(FLAG)
      RETURN
      END SUBROUTINE MUMPS_AB_COORD_TO_LMAT
      SUBROUTINE MUMPS_AB_LOCALCLEAN_LMAT (  MYID,
     &     LMAT, K147, FLAG, NBCOL, IFLAG, IERROR, LP, LPOK
     &     )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MYID, LP
      LOGICAL, INTENT(IN) :: LPOK
      INTEGER, INTENT(IN)  :: NBCOL, K147
      INTEGER, INTENT(OUT) :: FLAG(NBCOL)
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
      TYPE(LMATRIX_T), INTENT(INOUT)  :: LMAT
      INTEGER, POINTER, DIMENSION(:)     :: PTCLEAN
      INTEGER :: allocok, IB, JB, JJ, LMAT_FIRST, NBCOL_LOC
      INTEGER(8) :: NB_8, PT_CUR, SIZE_TO_ALLOC_8
      LMAT_FIRST = LMAT%FIRST
      NBCOL_LOC  = LMAT%NBCOL_LOC
      DO JB=1, NBCOL
        FLAG(JB) = 0
      ENDDO
      LMAT%NZL = 0_8
      DO JB=1, NBCOL_LOC, K147
       SIZE_TO_ALLOC_8 = 0_8
       DO JJ=JB, min(JB+K147-1,NBCOL_LOC)
         IF ( LMAT%COL(JJ)%NBINCOL.EQ.0) CYCLE
         DO IB=1,  LMAT%COL(JJ)%NBINCOL
          IF (FLAG(LMAT%COL(JJ)%IRN(IB)).EQ.JJ) THEN
           LMAT%COL(JJ)%IRN(IB)=0
          ELSE
           SIZE_TO_ALLOC_8 = SIZE_TO_ALLOC_8 + 1_8
           LMAT%NZL = LMAT%NZL+1_8
           FLAG(LMAT%COL(JJ)%IRN(IB)) = JJ
          ENDIF
         ENDDO
       ENDDO
       IF (SIZE_TO_ALLOC_8.GT.0) THEN
         ALLOCATE(PTCLEAN(SIZE_TO_ALLOC_8), STAT=allocok)
         IF (allocok.NE.0) THEN
           IFLAG     = -7
           CALL MUMPS_SET_IERROR(SIZE_TO_ALLOC_8, IERROR)
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate PTCLEAN of size", 
     &                     IERROR
           END IF
           RETURN
         ENDIF
         PT_CUR = 1_8
         DO JJ=JB, min(JB+K147-1,NBCOL_LOC)
           NB_8   = 0_8
           DO IB=1,  LMAT%COL(JJ)%NBINCOL
            IF (LMAT%COL(JJ)%IRN(IB).NE.0) THEN
             PTCLEAN(PT_CUR+NB_8)=LMAT%COL(JJ)%IRN(IB)
             NB_8 = NB_8 + 1_8
            ENDIF
           ENDDO
           LMAT%COL(JJ)%NBINCOL = int(NB_8)
           IF (JJ.GT.JB) 
     &     LMAT%COL(JJ)%IRN     => PTCLEAN(PT_CUR: PT_CUR+NB_8-1_8)
           PT_CUR = PT_CUR + NB_8
         ENDDO
         deallocate(LMAT%COL(JB)%IRN)
         LMAT%COL(JB)%IRN => PTCLEAN
         NULLIFY(PTCLEAN)
       ELSE
         if (associated(LMAT%COL(JB)%IRN)) 
     &           deallocate(LMAT%COL(JB)%IRN)
         NULLIFY(LMAT%COL(JB)%IRN)
       ENDIF
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_AB_LOCALCLEAN_LMAT
      SUBROUTINE MUMPS_AB_CLEANLMAT_TO_LUMAT( 
     &     LMAT, LUMAT, K147, INFO, ICNTL )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      IMPLICIT NONE
      TYPE(LMATRIX_T)        :: LMAT, LUMAT
      INTEGER, INTENT(IN)    :: K147
      INTEGER, INTENT(IN)    :: ICNTL(60)
      INTEGER, INTENT(INOUT) :: INFO(80)
      INTEGER    :: IB, JJ, IIB, JB, allocok, LP, MPG, NB
      INTEGER(8) :: NB_8, PT_CUR, SIZE_TO_ALLOC_8
      INTEGER, POINTER, DIMENSION(:) :: BLOCK_PTR
      LOGICAL LPOK, PROKG
      LP    = ICNTL( 1 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      MPG   = ICNTL( 3 )
      PROKG = ( MPG .GT. 0 .and. (ICNTL(4).GE.2) ) 
      LUMAT%NBCOL  = LMAT%NBCOL
      LUMAT%NZL    = 2_8*LMAT%NZL
      LUMAT%FIRST          = LMAT%FIRST
      LUMAT%NBCOL_LOC      = LMAT%NBCOL_LOC
      LUMAT%CONTIGUOUS_COL = LMAT%CONTIGUOUS_COL
      ALLOCATE( LUMAT%COL(LUMAT%NBCOL_LOC),STAT=allocok)
      IF (allocok.NE.0) THEN
           INFO( 1 ) = -7
           INFO( 2 ) = LUMAT%NBCOL_LOC
           IF ( LPOK ) THEN
             WRITE(LP, *) " ERROR allocating LUMAT%COL "
           END IF
           RETURN
      ENDIF
      DO JB=1,  LMAT%NBCOL_LOC
         LUMAT%COL(JB)%NBINCOL = LMAT%COL(JB)%NBINCOL
      ENDDO
      DO JB=1,  LMAT%NBCOL_LOC
        DO IB=1, LMAT%COL(JB)%NBINCOL
          IIB=LMAT%COL(JB)%IRN(IB)
          LUMAT%COL(IIB)%NBINCOL = LUMAT%COL(IIB)%NBINCOL + 1
        ENDDO
      ENDDO
      DO JB=1,  LUMAT%NBCOL_LOC, K147
        SIZE_TO_ALLOC_8 = 0_8
        DO JJ=JB, min(JB+K147-1,LUMAT%NBCOL_LOC)
         SIZE_TO_ALLOC_8 = SIZE_TO_ALLOC_8 + 
     &                     int(LUMAT%COL(JJ)%NBINCOL, 8)
        ENDDO
        ALLOCATE(LUMAT%COL(JB)%IRN(SIZE_TO_ALLOC_8), STAT=allocok)
        IF (allocok.NE.0) THEN
           INFO(1)  = -7
           CALL MUMPS_SET_IERROR(SIZE_TO_ALLOC_8, INFO(2))
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate block of ", 
     &                K147, " columns in LUMAT%COL",
     &                " of size ", SIZE_TO_ALLOC_8
           END IF
           RETURN
        ENDIF
        IF (K147.GT.1) THEN
          BLOCK_PTR => LUMAT%COL(JB)%IRN
          NB_8    = int(LUMAT%COL(JB)%NBINCOL, 8) 
          PT_CUR = 1_8 + NB_8 
          DO JJ=JB+1, min(JB+K147-1,LUMAT%NBCOL_LOC)
            NB_8 =  int(LUMAT%COL(JJ)%NBINCOL,8)
            IF (NB_8.GT.0) THEN
             LUMAT%COL(JJ)%IRN => BLOCK_PTR(PT_CUR: PT_CUR+NB_8-1_8)
             PT_CUR = PT_CUR + NB_8
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      DO JB=1, LMAT%NBCOL_LOC
        LUMAT%COL(JB)%NBINCOL = 0
      ENDDO
      DO JB=1, LMAT%NBCOL_LOC
        DO IB=1, LMAT%COL(JB)%NBINCOL
          IIB=LMAT%COL(JB)%IRN(IB)
          NB = LUMAT%COL(JB)%NBINCOL+1
          LUMAT%COL(JB)%NBINCOL = NB
          LUMAT%COL(JB)%IRN(NB)= IIB
          NB = LUMAT%COL(IIB)%NBINCOL+1
          LUMAT%COL(IIB)%NBINCOL = NB
          LUMAT%COL(IIB)%IRN(NB)= JB
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_AB_CLEANLMAT_TO_LUMAT
      SUBROUTINE MUMPS_AB_LMAT_TO_CLEAN_G( MYID, UNFOLD, 
     &     READY_FOR_ANA_F,
     &     LMAT, GCOMP, INFO, ICNTL
     &     , MEMCNT )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T, COMPACT_GRAPH_T
#if defined(DETERMINISTIC_PARALLEL_GRAPH)
      USE MUMPS_TOOLS_COMMON_M, ONLY : MUMPS_MERGESORT,
     &                                 MUMPS_MERGESWAP1
#endif
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: MYID
      LOGICAL, INTENT(IN)    :: UNFOLD, READY_FOR_ANA_F
      TYPE(LMATRIX_T)        :: LMAT
      TYPE(COMPACT_GRAPH_T)  :: GCOMP
      INTEGER, INTENT(IN)    :: ICNTL(60)
      INTEGER, INTENT(INOUT) :: INFO(80)
      INTEGER(8), INTENT(INOUT) :: MEMCNT
      INTEGER    :: IB, IIB, JJB, allocok, LP, MPG 
      INTEGER    :: JFIRST, NG_LOCAL
      INTEGER(8) :: JPOS, SIZEGCOMPALLOCATED
      INTEGER(8), ALLOCATABLE, DIMENSION(:) :: IQ
#if defined(DETERMINISTIC_PARALLEL_GRAPH)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: WORK
      INTEGER(8) :: IFIRST, ILAST
      INTEGER :: L
#endif
      LOGICAL LPOK, PROKG
      LP    = ICNTL( 1 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      MPG   = ICNTL( 3 )
      PROKG = ( MPG .GT. 0 .and. (ICNTL(4).GE.2) ) 
      GCOMP%NG    = LMAT%NBCOL
      GCOMP%FIRST = LMAT%FIRST
      GCOMP%LAST  = LMAT%FIRST+LMAT%NBCOL_LOC-1
      JFIRST     = GCOMP%FIRST
      NG_LOCAL   = LMAT%NBCOL_LOC
      IF (GCOMP%FIRST.GT.GCOMP%LAST) THEN 
       GCOMP%NZG   = 0_8
       RETURN
      ENDIF
      IF (UNFOLD) THEN 
       GCOMP%NZG   = 2_8*LMAT%NZL
       SIZEGCOMPALLOCATED = GCOMP%NZG + int(GCOMP%NG,8)+1_8
      ELSE IF (READY_FOR_ANA_F) THEN
       GCOMP%NZG = LMAT%NZL
       SIZEGCOMPALLOCATED = GCOMP%NZG + int(GCOMP%NG,8)+1_8
      ELSE
       GCOMP%NZG          = LMAT%NZL
       SIZEGCOMPALLOCATED = GCOMP%NZG
      ENDIF
      GCOMP%SIZEADJALLOCATED= SIZEGCOMPALLOCATED
      ALLOCATE( GCOMP%ADJ(SIZEGCOMPALLOCATED), 
     &          GCOMP%IPE(NG_LOCAL+1),
     &          IQ(NG_LOCAL),STAT=allocok)
      IF (allocok.NE.0) THEN
           INFO( 1 ) = -7
           CALL MUMPS_SET_IERROR(
     &        GCOMP%NZG + 3_8*int(GCOMP%NG,8)+1_8, INFO(2))
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocating graph in",
     &                     " MUMPS_AB_LMAT_TO_CLEAN_G"
           END IF
           RETURN
      ENDIF
      DO JJB=1, NG_LOCAL
         IQ(JJB)=0_8
      ENDDO
      IF (UNFOLD) THEN 
       DO JJB=1, NG_LOCAL
        DO IB=1, LMAT%COL(JJB)%NBINCOL
          IIB=LMAT%COL(JJB)%IRN(IB)
          IQ(JJB)=IQ(JJB)+1
          IQ(IIB)=IQ(IIB)+1
        ENDDO
       ENDDO
      ELSE
       DO JJB=1, NG_LOCAL
         IQ(JJB) = LMAT%COL(JJB)%NBINCOL
       ENDDO
      ENDIF
      GCOMP%IPE(1) = 1_8
      DO JJB=1, NG_LOCAL
        GCOMP%IPE(JJB+1) = GCOMP%IPE(JJB)+IQ(JJB)
      ENDDO
      IF (UNFOLD) THEN 
       DO JJB=1, GCOMP%NG
        IQ(JJB)= GCOMP%IPE(JJB)
       ENDDO
       DO JJB=1, NG_LOCAL
        DO IB=1, LMAT%COL(JJB)%NBINCOL
          IIB=LMAT%COL(JJB)%IRN(IB)
          GCOMP%ADJ(IQ(IIB))= JJB
          IQ(IIB)           = IQ(IIB)+1_8
          GCOMP%ADJ(IQ(JJB))= IIB
          IQ(JJB)           = IQ(JJB)+1_8
        ENDDO
       ENDDO
      ELSE
       DO JJB=1, NG_LOCAL
        JPOS =  GCOMP%IPE(JJB)
        DO IB=1, LMAT%COL(JJB)%NBINCOL
          IIB=LMAT%COL(JJB)%IRN(IB)
          GCOMP%ADJ(JPOS)= IIB
          JPOS           = JPOS+1_8
        ENDDO
       ENDDO
      ENDIF
      DEALLOCATE(IQ)
#if defined(DETERMINISTIC_PARALLEL_GRAPH)
      IF (.NOT.READY_FOR_ANA_F) THEN
        ALLOCATE(WORK(0:GCOMP%NG),stat=allocok)
        IF (allocok.NE.0) THEN
          INFO( 1 ) = -7
          INFO( 2 ) = GCOMP%NG+1
          IF ( LPOK ) THEN
             WRITE(LP, *) " ERROR allocating WORK in",
     &                    " MUMPS_AB_LMAT_TO_CLEAN_G"
          END IF
          RETURN
        ENDIF
        DO JJB=1, NG_LOCAL
          IFIRST = GCOMP%IPE(JJB)
          ILAST= GCOMP%IPE(JJB+1)-1
          L = int(ILAST-IFIRST+1)
          IF ( L .GE. 2 ) THEN 
            IF (L .GE. GCOMP%NG ) THEN
              WRITE(*,*) " Internal error in MUMPS_AB_LMAT_TO_CLEAN_G",
     &        L, GCOMP%NG
              CALL MUMPS_ABORT()
            ENDIF
            CALL MUMPS_MERGESORT( L,
     &      GCOMP%ADJ(IFIRST:ILAST), WORK(0:L+1) )
            CALL MUMPS_MERGESWAP1( L,
     &      WORK(0:L+1), GCOMP%ADJ(IFIRST:ILAST) )
          ENDIF
        ENDDO
        DEALLOCATE(WORK)
      ENDIF
#endif
      MEMCNT = MEMCNT + GCOMP%SIZEADJALLOCATED 
     &          + (GCOMP%LAST-GCOMP%FIRST+1)+1 
     &          + 5
      RETURN
      END SUBROUTINE MUMPS_AB_LMAT_TO_CLEAN_G
      SUBROUTINE MUMPS_AB_COL_DISTRIBUTION ( OPTION,
     &     INFO, ICNTL, COMM, NBLK, MYID, NPROCS,
     &     LMAT, MAPCOL, CONTIGUOUS_COL )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR
      INTEGER, INTENT(IN) :: OPTION, NBLK
      INTEGER, INTENT(IN) :: ICNTL(60), COMM, MYID, NPROCS
      INTEGER             :: INFO(80)
      TYPE(LMATRIX_T)     :: LMAT
      INTEGER, INTENT(OUT):: MAPCOL(NBLK)
      LOGICAL, INTENT(OUT):: CONTIGUOUS_COL
      INTEGER    :: LP, SIZE_NZROW, I
      LOGICAL    :: LPOK
      INTEGER(8) :: NZL, NNZ
      INTEGER, DIMENSION(:), ALLOCATABLE :: NZ_ROW
      LP  = ICNTL( 1 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      IF (OPTION.EQ.1) THEN
       NNZ        = -9999 
       SIZE_NZROW = 1
      ELSE
       NZL  = LMAT%NZL
       SIZE_NZROW = NBLK
      ENDIF
      ALLOCATE(NZ_ROW(SIZE_NZROW), STAT=IERR)
      IF (IERR.NE.0) THEN
           INFO(1)  = -7
           INFO(2)  = SIZE_NZROW
           IF ( LPOK ) THEN
              WRITE(LP, *) 
     &     " ERROR allocate in MUMPS_AB_COL_DISTRIBUTION ", INFO(2)
           END IF
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF (INFO(1).LT.0) GOTO 500
      IF (OPTION.NE.1) THEN
        DO I = 1, LMAT%NBCOL_LOC
         MAPCOL(I) = LMAT%COL(I)%NBINCOL
        ENDDO
        CALL MUMPS_BIGALLREDUCE(.FALSE.,
     &        MAPCOL(1), NZ_ROW(1), NBLK, 
     &        MPI_INTEGER, MPI_SUM, COMM, IERR)
        CALL MPI_ALLREDUCE(NZL, NNZ, 1, 
     &        MPI_INTEGER8, MPI_SUM, COMM, IERR)
      ENDIF
      CALL MUMPS_AB_COMPUTE_MAPCOL (OPTION, INFO, ICNTL, MYID,
     &   NNZ, NZ_ROW(1), SIZE_NZROW, NBLK, NPROCS, MAPCOL(1),
     &   CONTIGUOUS_COL)
 500  CONTINUE
      IF (allocated(NZ_ROW)) DEALLOCATE(NZ_ROW)
      RETURN
      END SUBROUTINE MUMPS_AB_COL_DISTRIBUTION
      SUBROUTINE MUMPS_AB_COMPUTE_MAPCOL (OPTION, INFO, ICNTL, 
     &    MYID, NNZ, NZ_ROW, SIZE_NZROW, NBLK,  NPROCS, MAPCOL,
     &    CONTIGUOUS_COL )
#if ! defined(NOANA_BLKAUTO)
#endif
      INTEGER, INTENT(IN) :: OPTION, MYID, SIZE_NZROW, NBLK
      INTEGER, INTENT(IN) :: ICNTL(60), NPROCS
      INTEGER             :: INFO(80)
      INTEGER(8)          :: NNZ
      INTEGER, INTENT(IN) :: NZ_ROW(SIZE_NZROW)
      INTEGER, INTENT(OUT):: MAPCOL(NBLK)
      LOGICAL, INTENT(OUT) :: CONTIGUOUS_COL
      INTEGER    :: I, J, P, F, LP, IERR
      LOGICAL    :: LPOK
      INTEGER(8) :: SHARE, T
      INTEGER, DIMENSION(:), ALLOCATABLE :: FIRST
      LP  = ICNTL( 1 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      ALLOCATE(FIRST(NPROCS+1), STAT=IERR)
      IF (IERR.NE.0) THEN
           INFO(1)  = -7
           INFO(2)  = NPROCS+1
           IF ( LPOK ) THEN
              WRITE(LP, *) 
     &     " ERROR allocate in MUMPS_AB_COL_DISTRIBUTION ", INFO(2)
           END IF
           GOTO 500
      ENDIF
      DO I=1,NPROCS+1
       FIRST(I) = 0
      ENDDO
      CONTIGUOUS_COL = .TRUE.
      IF (OPTION.EQ.1) THEN
       SHARE = int(NBLK/NPROCS,8)
       DO I=1, NPROCS
          FIRST(I) = (I-1)*int(SHARE)+1
       END DO
       FIRST(NPROCS+1)=NBLK+1
      ELSE
       SHARE = (NNZ-1_8)/int(NPROCS,8) + 1_8
         P = 0
         T = 0_8
         F = 1
         DO I=1, NBLK
            T = T+int(NZ_ROW(I),8)
            IF (
     &           (T .GE. SHARE) .OR.
     &           ((NBLK-I).EQ.(NPROCS-P-1)) .OR.
     &           (I.EQ.NBLK)
     &           ) THEN
               P             = P+1
               IF(P.EQ.NPROCS) THEN
                  FIRST(P) = F
                  EXIT
               ELSE
                  FIRST(P) = F
                  F             = I+1
                  T             = 0_8
               END IF
            END IF
            IF ((I.EQ.NBLK).AND.(P.NE.NPROCS)) THEN 
             DO J=P,NPROCS
              FIRST(J) = FIRST(P)
             ENDDO
            ENDIF
         END DO
         FIRST(NPROCS+1) = NBLK+1
      ENDIF
      DO I=1,NPROCS
        DO J=FIRST(I), FIRST(I+1)-1
          MAPCOL(J) = I-1
        ENDDO
      ENDDO
      IF (allocated(FIRST))  DEALLOCATE(FIRST)
 500  CONTINUE
      RETURN
      END SUBROUTINE MUMPS_AB_COMPUTE_MAPCOL
      SUBROUTINE MUMPS_AB_CHECK_FINALIZE_MAPCOL (
     &            MYID, COMM, NPROCS, NPROCS_PARAORD,
     &            MAPCOL, SIZE_MAPCOL, 
     &            MAPCOL_VALIDATED, CONTIGUOUS_COL,
     &            LP, LPOK, INFO )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER, INTENT(IN)    :: MYID, COMM, NPROCS, SIZE_MAPCOL
      INTEGER, INTENT(IN)    :: LP, NPROCS_PARAORD
      INTEGER, INTENT(INOUT) :: MAPCOL(SIZE_MAPCOL), INFO(2)
      LOGICAL, INTENT(OUT)   :: MAPCOL_VALIDATED, CONTIGUOUS_COL
      LOGICAL, INTENT(IN)    :: LPOK
      INTEGER :: I, IP, IERR, IFIRST, REMAP(0:1)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: VISITED
#if defined(AVOID_MPI_IN_PLACE)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAPCOL_loc
#endif
      INTEGER :: allocok
      MAPCOL_VALIDATED=.FALSE.
      ALLOCATE( 
#if defined(AVOID_MPI_IN_PLACE)
     &         MAPCOL_loc(SIZE_MAPCOL), 
#endif
     &         VISITED(NPROCS), 
     &         STAT=allocok) 
      IF (allocok.NE.0) THEN
             INFO(1)  = -7
             INFO(2)  =  NPROCS
#if defined(AVOID_MPI_IN_PLACE)
     &                   + SIZE_MAPCOL
#endif
             IF ( LPOK ) THEN
                WRITE(LP, *) 
     &       " ERROR allocate in MUMPS_AB_CHECK_FINALIZE_MAPCOL",
     &       " of size ", INFO(2)
             END IF
      ENDIF
      CALL MUMPS_PROPINFO(  LP, INFO(1),
     &       COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
#if defined(AVOID_MPI_IN_PLACE)
      MAPCOL_loc = MAPCOL
#endif
#if defined(AVOID_MPI_IN_PLACE)
      CALL MUMPS_BIGALLREDUCE(.FALSE., MAPCOL_loc, 
#else
      CALL MUMPS_BIGALLREDUCE(.TRUE., MPI_IN_PLACE,
#endif
     &        MAPCOL(1), SIZE_MAPCOL,
     &        MPI_INTEGER, MPI_MAX, COMM, IERR)
      MAPCOL_VALIDATED=.TRUE.
      IF (MAPCOL(1).EQ.0) THEN
       DO I=2, SIZE_MAPCOL
         IF (MAPCOL(I).NE.0) THEN 
              MAPCOL(1)=MAPCOL(I)
              EXIT
         ENDIF
       ENDDO
      ENDIF
      IF (MAPCOL(1).EQ.0) THEN
        MAPCOL_VALIDATED=.FALSE.
        GOTO 500
      ENDIF
      DO I=2, SIZE_MAPCOL
        IF (MAPCOL(I).EQ.0) MAPCOL(I) = MAPCOL(I-1)
      ENDDO
      CONTIGUOUS_COL = .TRUE.
      VISITED        = 0
      IFIRST         = 1
      IP             = MAPCOL(1)
      DO I=2, SIZE_MAPCOL
       IF (VISITED(MAPCOL(I)).NE.0) THEN
         CONTIGUOUS_COL=.FALSE.
         EXIT
       ENDIF
       IF (MAPCOL(I).NE.IP) THEN 
            VISITED(IP) = I-IFIRST
            IP          = MAPCOL(I)
            IFIRST      = I
       ENDIF
       IF (I.EQ.SIZE_MAPCOL) THEN
             VISITED(IP) = SIZE_MAPCOL-IFIRST+1
       ENDIF
      ENDDO
      IF (.NOT.CONTIGUOUS_COL) then 
      ENDIF
 500  CONTINUE
      IF ( MAPCOL_VALIDATED ) THEN 
       DO I = 1, SIZE_MAPCOL
         MAPCOL(I) = MAPCOL(I) - 1
       ENDDO
       IF (NPROCS_PARAORD.LT.NPROCS) THEN
        IF (.NOT.CONTIGUOUS_COL) THEN
          REMAP(0)= 0
          REMAP(1)= NPROCS_PARAORD-1
          DO I=1, SIZE_MAPCOL
            IP = MAPCOL(I)
            IF (IP.GE.NPROCS_PARAORD) THEN
             MAPCOL(I) = REMAP (mod(IP,2))
            ENDIF
          ENDDO
        ELSE
          IF (MAPCOL(1).GE.NPROCS_PARAORD) THEN
            DO I=2, SIZE_MAPCOL
              IF (MAPCOL(I).LT.NPROCS_PARAORD) THEN
                MAPCOL(1)=MAPCOL(I)
                EXIT
              ENDIF
            ENDDO
            IF (MAPCOL(1).GE.NPROCS_PARAORD) THEN
             MAPCOL(1)= 0
            ENDIF
          ENDIF
          DO I=2, SIZE_MAPCOL
            IP = MAPCOL(I)
            IF (IP.GE.NPROCS_PARAORD) THEN
              MAPCOL(I) = MAPCOL(I-1)
            ENDIF
          ENDDO
         ENDIF
       ENDIF
      ENDIF
#if defined(AVOID_MPI_IN_PLACE)
      IF (allocated(MAPCOL_loc)) DEALLOCATE(MAPCOL_loc)
#endif
      IF (allocated(VISITED))    DEALLOCATE(VISITED)
      RETURN
      END SUBROUTINE MUMPS_AB_CHECK_FINALIZE_MAPCOL
      SUBROUTINE MUMPS_AB_BUILD_DCLEAN_LUMATRIX ( 
     &     MAPCOLonLUMAT, MAPCOL_IN_NSTEPS, CONTIGUOUS_COL,
     &     INFO, ICNTL, KEEP, COMM, MYID, NBLK, NPROCS,
     &     LMAT, MAPCOL, SIZEMAPCOL,
     &     STEP, SIZESTEP, 
     &     LUMAT)
      USE MUMPS_ANA_BLK_M 
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      LOGICAL, INTENT(IN) :: MAPCOLonLUMAT, MAPCOL_IN_NSTEPS
      LOGICAL, INTENT(INOUT) :: CONTIGUOUS_COL
      INTEGER, INTENT(IN) :: MYID, NPROCS, NBLK, SIZEMAPCOL
      INTEGER, INTENT(IN) :: ICNTL(60), COMM, KEEP(500)
      INTEGER, INTENT(IN) :: SIZESTEP
      INTEGER, INTENT(IN) :: STEP(SIZESTEP)
      INTEGER, INTENT(INOUT) ::  INFO(80)
      TYPE(LMATRIX_T), INTENT(INOUT)  :: LMAT
      INTEGER, INTENT(INOUT)          :: MAPCOL(SIZEMAPCOL)
      TYPE(LMATRIX_T), INTENT(OUT)    :: LUMAT
      INTEGER ::      NBLKloc, IERR, JB, IB, LP, NB, I,
     &                NBRECORDS, JJ
      INTEGER(8)  ::  NNZ, NZ_locMAX8, NSEND8, NLOCAL8
      INTEGER(8) :: NB_8, PT_CUR, SIZE_TO_ALLOC_8
      INTEGER, POINTER, DIMENSION(:) :: BLOCK_PTR
      LOGICAL :: LPOK
      INTEGER, ALLOCATABLE, DIMENSION(:) :: WT, WNBINCOL
      INTEGER OPTION
      PARAMETER (OPTION=2)
      NBLKloc = LMAT%NBCOL
      IF (NBLKloc.NE.NBLK) THEN
       write(6,*) "Internal error in MUMPS_AB_BUILD_DCLEAN_LUMATRIX ",
     &       "NBCOL, NBCOL_LOC, NBLK=", LMAT%NBCOL_LOC, NBLKloc, NBLK
      ENDIF
      LP  = ICNTL( 1 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      ALLOCATE(WT(NBLK),WNBINCOL(NBLK),STAT=IERR)
      IF (IERR.NE.0) THEN
           INFO(1) = -7
           INFO(2) = 2*LMAT%NBCOL_LOC
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate of LUMAT%COL; WT"
           END IF
      ENDIF
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
       IF ( INFO(1) .LT. 0 ) GOTO 500
      WT = 0
      DO JB=1, LMAT%NBCOL_LOC
       WT(JB) = LMAT%COL(JB)%NBINCOL
      ENDDO
      DO JB=1,LMAT%NBCOL_LOC
       IF ( LMAT%COL(JB)%NBINCOL.EQ.0) CYCLE
       DO IB=1,  LMAT%COL(JB)%NBINCOL
        I = LMAT%COL(JB)%IRN(IB)
        WT(I)= WT(I)+1
       ENDDO
      ENDDO
      CALL MUMPS_BIGALLREDUCE(.FALSE.,
     &        WT(1), WNBINCOL(1), NBLK, 
     &        MPI_INTEGER, MPI_SUM, COMM, IERR)
      IF (allocated(WT)) DEALLOCATE(WT)
      IF (MAPCOLonLUMAT) THEN
       NNZ = 0_8
       DO I=1, NBLK
        NNZ=NNZ+int(WNBINCOL(I),8)
       ENDDO
       CALL  MUMPS_AB_COMPUTE_MAPCOL (OPTION, INFO, ICNTL,
     &     MYID, NNZ, WNBINCOL(1), NBLK, 
     &     NBLK, NPROCS, MAPCOL(1), CONTIGUOUS_COL)
       CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
       IF ( INFO(1) .LT. 0 ) GOTO 500
      ENDIF
      LUMAT%NBCOL          = NBLK
      LUMAT%NZL            = 0_8
      LUMAT%FIRST          = 1
      LUMAT%NBCOL_LOC      = NBLK
      LUMAT%CONTIGUOUS_COL = CONTIGUOUS_COL
      IF (CONTIGUOUS_COL) THEN
        LUMAT%NBCOL_LOC = 0
        LUMAT%FIRST     = -9991
        DO JB=1,NBLK
         IF (MAPCOL(JB).EQ.MYID) THEN
           IF ( LUMAT%FIRST .LT.0 ) LUMAT%FIRST=JB
           LUMAT%NBCOL_LOC = LUMAT%NBCOL_LOC + 1
         ENDIF
        ENDDO
        IF (LUMAT%FIRST.LT.0) THEN
          LUMAT%FIRST     = 1
        ENDIF
      ENDIF
      IF (LUMAT%NBCOL_LOC.GT.0) THEN
       ALLOCATE(LUMAT%COL(LUMAT%NBCOL_LOC), STAT=IERR)
       IF (IERR.NE.0) THEN
           INFO(1) = -7
           INFO(2) = LUMAT%NBCOL_LOC
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate of LUMAT%COL; WT"
           END IF
       ENDIF
      ENDIF
      IF ( INFO(1) .GE. 0 ) THEN
       DO JB=1,NBLK
        NB =  WNBINCOL(JB)
        IF (MAPCOL_IN_NSTEPS) THEN
         IF (MAPCOL(abs(STEP(JB))).EQ.MYID) THEN
          LUMAT%NZL       = LUMAT%NZL + int(NB,8)
         ELSE
          NB = 0
         ENDIF
         LUMAT%COL(JB)%NBINCOL = NB
        ELSE
         IF (LUMAT%NBCOL.EQ.LUMAT%NBCOL_LOC) 
     &      LUMAT%COL(JB)%NBINCOL = 0
         IF (MAPCOL(JB).EQ.MYID) THEN
          LUMAT%NZL       = LUMAT%NZL + int(NB,8)
          LUMAT%COL(JB-LUMAT%FIRST+1)%NBINCOL = NB
         ENDIF
        ENDIF
       ENDDO
       DO JB=1,LUMAT%NBCOL_LOC, KEEP(147)
        SIZE_TO_ALLOC_8 = 0_8
        DO JJ=JB, min(JB+KEEP(147)-1,LUMAT%NBCOL_LOC)
         SIZE_TO_ALLOC_8 = SIZE_TO_ALLOC_8 + LUMAT%COL(JJ)%NBINCOL
        ENDDO
        IF (SIZE_TO_ALLOC_8.GT.0) THEN
         ALLOCATE(LUMAT%COL(JB)%IRN(SIZE_TO_ALLOC_8), 
     &            STAT=IERR)
         IF (IERR.NE.0) THEN
           INFO(1)  = -7
           CALL MUMPS_SET_IERROR(SIZE_TO_ALLOC_8, INFO(2))
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate of LUMAT%COL"
           END IF
           EXIT
         ENDIF
        ENDIF
        IF (KEEP(147).GT.1.AND.SIZE_TO_ALLOC_8.GT.0) THEN
         BLOCK_PTR => LUMAT%COL(JB)%IRN
         NB_8    = int(LUMAT%COL(JB)%NBINCOL, 8) 
         PT_CUR = 1_8 + NB_8 
         DO JJ=JB+1, min(JB+KEEP(147)-1,LUMAT%NBCOL_LOC)
           NB_8 =  int(LUMAT%COL(JJ)%NBINCOL,8)
           IF (NB_8.GT.0) THEN
            LUMAT%COL(JJ)%IRN => 
     &                      BLOCK_PTR(PT_CUR: PT_CUR+NB_8-1_8)
            PT_CUR = PT_CUR + NB_8
           ENDIF
         ENDDO
        ENDIF
       ENDDO
      ENDIF
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      IF (allocated(WNBINCOL)) DEALLOCATE(WNBINCOL)
      CALL MPI_ALLREDUCE(LUMAT%NZL, NZ_locMAX8, 1, MPI_INTEGER8,
     &                   MPI_MAX, COMM, IERR)
      NBRECORDS = KEEP(39)
      IF (NZ_locMAX8 .LT. int(NBRECORDS,8)) THEN
            NBRECORDS = int(NZ_locMAX8)
      ENDIF
      CALL MUMPS_AB_DIST_LMAT_TO_LUMAT ( 
     &  .TRUE.,   
     &  MAPCOL_IN_NSTEPS,  
     &  INFO, ICNTL, COMM, MYID, NBLK, NPROCS,
     &  LMAT, MAPCOL, SIZEMAPCOL, STEP, SIZESTEP, 
     &  LUMAT, NBRECORDS, NSEND8, NLOCAL8
     &  )
      CALL MUMPS_AB_FREE_LMAT(LMAT, KEEP(147))
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      ALLOCATE(WT(LUMAT%NBCOL), STAT=IERR)
      IF (IERR.NE.0) THEN
           INFO(1) = -7
           INFO(2) = 2*NBLK
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate of LUMAT%COL; WT"
           END IF
           GOTO 400
      ENDIF
      CALL MUMPS_AB_LOCALCLEAN_LMAT ( MYID,
     &     LUMAT, KEEP(147), 
     &     WT(1), LUMAT%NBCOL, INFO(1), INFO(2), LP, LPOK
     & )
      DEALLOCATE(WT)
 400  CONTINUE
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      GOTO 600
  500 CONTINUE
      IF (allocated(WT)) DEALLOCATE(WT)
      IF (allocated(WNBINCOL)) DEALLOCATE(WNBINCOL)
  600 CONTINUE
      RETURN
      END SUBROUTINE MUMPS_AB_BUILD_DCLEAN_LUMATRIX
      SUBROUTINE MUMPS_INI_REDIST_LUMAT (
     &  INFO, ICNTL, KEEP, COMM, MYID, NBLK, 
     &  NPROCS, MAPCOL_IN_NSTEPS, 
     &  LUMAT, PROCNODE_STEPS, NSTEPS, 
     &  MAPCOL, SIZE_MAPCOL,
     &  LUMAT_REMAP, NBRECORDS, STEP,
     &  NPROCS_PARAORD, PARAORD_to_idCOMM, 
     &  SIZE_PARAORD_to_idCOMM
     &  )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER  :: IERR, MASTER
      PARAMETER (MASTER=0)
      INTEGER, INTENT(IN) :: MYID, NBLK, NSTEPS, KEEP(500)
      INTEGER, INTENT(IN)   :: ICNTL(60), COMM, NPROCS
      INTEGER               :: INFO(80)
      INTEGER, INTENT(IN) :: PROCNODE_STEPS(NSTEPS)
      TYPE(LMATRIX_T), INTENT(IN)  :: LUMAT
      INTEGER, INTENT(IN)          :: STEP(NBLK)
      TYPE(LMATRIX_T), INTENT(INOUT) :: LUMAT_REMAP
      INTEGER, INTENT(OUT) :: NBRECORDS
      LOGICAL, INTENT(IN)  :: MAPCOL_IN_NSTEPS
      INTEGER, INTENT(IN)  :: SIZE_MAPCOL
      INTEGER, INTENT(OUT) :: MAPCOL(SIZE_MAPCOL)
      INTEGER, INTENT(IN)  :: SIZE_PARAORD_to_idCOMM
      INTEGER, INTENT(IN)   :: NPROCS_PARAORD, 
     &       PARAORD_to_idCOMM(SIZE_PARAORD_to_idCOMM)
      INTEGER :: I, LP, MP, ISTEP, JB, JJ, NB
      INTEGER(8) :: NB_8, PT_CUR, SIZE_TO_ALLOC_8, NNZ
      INTEGER, POINTER, DIMENSION(:) :: BLOCK_PTR
      LOGICAL  :: CONTIGUOUS_COL
      LOGICAL :: LPOK
      INTEGER, ALLOCATABLE, DIMENSION(:) :: WT, WNBINCOL
      INTEGER MUMPS_PROCNODE
      INTEGER(8) :: NZ_locMAX8
      INTEGER OPTION
      LP  = ICNTL( 1 )
      MP  = ICNTL( 2 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      CONTIGUOUS_COL = .FALSE.
      ALLOCATE(WT(NBLK), WNBINCOL(NBLK), STAT=IERR)
      IF (IERR.NE.0) THEN
        INFO(1) = -7
        INFO(2) = 2*NBLK
        IF ( LPOK ) THEN
           WRITE(LP, *) " ERROR allocate WT"
        END IF
      ENDIF
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      WT = 0
      DO JB=1, LUMAT%NBCOL_LOC
       WT(JB+LUMAT%FIRST-1) = LUMAT%COL(JB)%NBINCOL
      ENDDO
      CALL MUMPS_BIGALLREDUCE(.FALSE., WT(1), WNBINCOL(1),
     &     NBLK, MPI_INTEGER, MPI_SUM, COMM, IERR)
      IF (allocated(WT)) DEALLOCATE(WT)
      IF (MAPCOL_IN_NSTEPS) THEN
        CONTIGUOUS_COL = .FALSE.
        IF (MYID.EQ.MASTER) THEN
         IF (KEEP(381).GT.1) THEN
           DO ISTEP=1, NSTEPS
            MAPCOL(ISTEP) = mod(ISTEP,NPROCS)
           ENDDO
         ELSE
           DO ISTEP=1, NSTEPS
            MAPCOL(ISTEP) = 
     &       MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),KEEP(199)) 
           ENDDO
         ENDIF
        ENDIF
        CALL MPI_BCAST( MAPCOL(1), NSTEPS, MPI_INTEGER,
     &       MASTER, COMM, IERR )
        CALL MPI_BCAST( STEP(1), NBLK, MPI_INTEGER,
     &       MASTER, COMM, IERR )
      ELSE
       NNZ = 0_8
       DO I=1, NBLK
        NNZ=NNZ+int(WNBINCOL(I),8)
       ENDDO
       IF (KEEP(244).NE.1) THEN
         IF (NNZ.LT.1000*NPROCS_PARAORD) THEN
           OPTION = 1 
         ELSE
           OPTION = 2 
         ENDIF
       ELSE
         OPTION = 2   
       ENDIF
       CALL  MUMPS_AB_COMPUTE_MAPCOL (OPTION, INFO, ICNTL,
     &     MYID, NNZ, WNBINCOL(1), NBLK, 
     &     NBLK, NPROCS_PARAORD, MAPCOL(1), CONTIGUOUS_COL)
       CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
       IF ( INFO(1) .LT. 0 ) GOTO 500
      ENDIF
      IF (NPROCS_PARAORD.LT.NPROCS)
     &  THEN
         DO I=1, SIZE_MAPCOL
           MAPCOL(I) = PARAORD_to_idCOMM(MAPCOL(I)+1)
         ENDDO
      ENDIF
      LUMAT_REMAP%NZL    = 0_8  
      LUMAT_REMAP%NBCOL  = NBLK
      LUMAT_REMAP%CONTIGUOUS_COL     = CONTIGUOUS_COL
      IF (CONTIGUOUS_COL) THEN
        LUMAT_REMAP%NBCOL_LOC = 0
        LUMAT_REMAP%FIRST     = -9991
        DO JB=1,NBLK
         IF (MAPCOL(JB).EQ.MYID) THEN
           IF ( LUMAT_REMAP%FIRST .LT.0 ) LUMAT_REMAP%FIRST=JB
           LUMAT_REMAP%NBCOL_LOC = LUMAT_REMAP%NBCOL_LOC + 1
         ENDIF
        ENDDO
        IF (LUMAT_REMAP%FIRST.LT.0) THEN
          LUMAT_REMAP%FIRST     = 1
        ENDIF
      ELSE
        LUMAT_REMAP%NBCOL_LOC          = NBLK
        LUMAT_REMAP%FIRST              = 1
      ENDIF
      IF (LUMAT_REMAP%NBCOL_LOC.GT.0) THEN
        ALLOCATE(LUMAT_REMAP%COL(LUMAT_REMAP%NBCOL_LOC), STAT=IERR)
        IF (IERR.NE.0) THEN
             INFO(1) = -7
             INFO(2) = LUMAT_REMAP%NBCOL_LOC
             IF ( LPOK ) THEN
                WRITE(LP, *) " ERROR allocate of LUMAT_REMAP%COL"
             END IF
        ENDIF
      ENDIF
      IF ( INFO(1) .GE. 0 ) THEN
       DO JB=1,NBLK
        NB =  WNBINCOL(JB)
        IF (MAPCOL_IN_NSTEPS) THEN
         IF (MAPCOL(abs(STEP(JB))).EQ.MYID) THEN
          LUMAT_REMAP%NZL       = LUMAT_REMAP%NZL + int(NB,8)
         ELSE
          NB = 0
         ENDIF
         LUMAT_REMAP%COL(JB)%NBINCOL = NB
        ELSE
         IF (LUMAT_REMAP%NBCOL.EQ.LUMAT_REMAP%NBCOL_LOC) 
     &      LUMAT_REMAP%COL(JB)%NBINCOL = 0
         IF (MAPCOL(JB).EQ.MYID) THEN
          LUMAT_REMAP%NZL       = LUMAT_REMAP%NZL + int(NB,8)
          LUMAT_REMAP%COL(JB-LUMAT_REMAP%FIRST+1)%NBINCOL = NB
         ENDIF
        ENDIF
       ENDDO
       DO JB=1,LUMAT_REMAP%NBCOL_LOC, KEEP(147)
        SIZE_TO_ALLOC_8 = 0_8
        DO JJ=JB, min(JB+KEEP(147)-1,LUMAT_REMAP%NBCOL_LOC)
         SIZE_TO_ALLOC_8 = SIZE_TO_ALLOC_8 + LUMAT_REMAP%COL(JJ)%NBINCOL
        ENDDO
        IF (SIZE_TO_ALLOC_8.GT.0) THEN
         ALLOCATE(LUMAT_REMAP%COL(JB)%IRN(SIZE_TO_ALLOC_8), 
     &            STAT=IERR)
         IF (IERR.NE.0) THEN
           INFO(1)  = -7
           CALL MUMPS_SET_IERROR(SIZE_TO_ALLOC_8, INFO(2))
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate of LUMAT_REMAP%COL"
           END IF
           EXIT
         ENDIF
        ENDIF
        IF (KEEP(147).GT.1.AND.SIZE_TO_ALLOC_8.GT.0) THEN
         BLOCK_PTR => LUMAT_REMAP%COL(JB)%IRN
         NB_8    = int(LUMAT_REMAP%COL(JB)%NBINCOL, 8) 
         PT_CUR = 1_8 + NB_8 
         DO JJ=JB+1, min(JB+KEEP(147)-1,LUMAT_REMAP%NBCOL_LOC)
           NB_8 =  int(LUMAT_REMAP%COL(JJ)%NBINCOL,8)
           IF (NB_8.GT.0) THEN
            LUMAT_REMAP%COL(JJ)%IRN => 
     &                      BLOCK_PTR(PT_CUR: PT_CUR+NB_8-1_8)
            PT_CUR = PT_CUR + NB_8
           ENDIF
         ENDDO
        ENDIF
       ENDDO
      ENDIF
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      IF (allocated(WNBINCOL)) DEALLOCATE(WNBINCOL)
       CALL MPI_ALLREDUCE(LUMAT_REMAP%NZL, NZ_locMAX8, 1, MPI_INTEGER8,
     &                   MPI_MAX, COMM, IERR)
        NBRECORDS = KEEP(39)
      IF (NZ_locMAX8 .LT. int(NBRECORDS,8)) THEN
            NBRECORDS = int(NZ_locMAX8)
      ENDIF
      GOTO 600
  500 CONTINUE
      IF (allocated(WT)) DEALLOCATE(WT)
      IF (allocated(WNBINCOL)) DEALLOCATE(WNBINCOL)
  600 CONTINUE
      RETURN
      END SUBROUTINE MUMPS_INI_REDIST_LUMAT
      SUBROUTINE MUMPS_AB_DCOORD_TO_DCOMPG ( 
     &     MYID, NPROCS, COMM, 
     &     NBLK, NDOF, NNZ, 
     &     IRN, JCN, DOF2BLOCK, 
     &     ICNTL, INFO, KEEP, 
     &     LUMAT, LUMAT_AVAIL, GCOMP, READY_FOR_ANA_F
#if ! defined(NOANA_BLKAUTO)
     &     , BLKFORMAT_DETECTED
#endif
     &     )
#if ! defined(NOANA_BLKAUTO)
#endif
      USE MUMPS_ANA_BLK_M, ONLY: LMATRIX_T, COMPACT_GRAPH_T
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, MASTER
      PARAMETER( MASTER = 0 )
      INTEGER, INTENT(IN)    :: MYID, NPROCS, NDOF
      INTEGER, INTENT(INOUT) :: NBLK
      INTEGER(8), INTENT(IN) :: NNZ
      INTEGER, INTENT(IN)    :: IRN(max(1_8,NNZ)), JCN(max(1_8,NNZ))
      LOGICAL, INTENT(IN)    :: READY_FOR_ANA_F
      INTEGER, INTENT(INOUT) :: DOF2BLOCK(NDOF)
      INTEGER, INTENT(IN)    :: ICNTL(60), COMM
      INTEGER, INTENT(INOUT) :: KEEP(500), INFO(80)
#if ! defined(NOANA_BLKAUTO)
      LOGICAL, INTENT(INOUT) :: BLKFORMAT_DETECTED
#endif
      TYPE(COMPACT_GRAPH_T)  :: GCOMP
      TYPE(LMATRIX_T)        :: LUMAT  
      LOGICAL, INTENT(OUT)   :: LUMAT_AVAIL
      TYPE(LMATRIX_T)        :: LMAT
      INTEGER :: IDUMMY_ARRAY(1)
      INTEGER :: allocok, LP, MPG
      LOGICAL :: LPOK, PROKG
      INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: MAPCOL
#if ! defined(NOANA_BLKAUTO)
      INTEGER, DIMENSION(:), POINTER             :: MAPCOL_PTR
      INTEGER :: SIZE_MAPCOL, K244
      INTEGER, TARGET :: IDUMMY_TARGET(1)
      LOGICAL :: MAPCOL_VALIDATED, REMAP_CONTIG
#endif
      LOGICAL :: MAPCOLonLUMAT, MAPCOL_IN_NSTEPS, CONTIGUOUS_COL,
     &           MAPCOLallocated
      INTEGER(8) :: MEMCNT, IDUMMY8
      INTEGER OPTION
      LP  = ICNTL( 1 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      MPG = ICNTL( 3 )
      PROKG = ( MPG .GT. 0 .and. MYID .eq. MASTER )
      IF (KEEP(13).EQ.2) THEN 
         OPTION= 3
      ELSE
        OPTION =KEEP(337) 
      ENDIF
#if ! defined(NOANA_BLKAUTO)
      IF (KEEP(13).EQ.2) THEN 
        BLKFORMAT_DETECTED = .FALSE.
      ENDIF
#endif
      MEMCNT = 0_8
      MAPCOLonLUMAT        = .FALSE.
      MAPCOL_IN_NSTEPS     = .FALSE.
      CONTIGUOUS_COL       = .TRUE.
      MAPCOLallocated      = .FALSE.
#if ! defined(NOANA_BLKAUTO)
      REMAP_CONTIG         = .FALSE.
      MAPCOL_VALIDATED     = .FALSE.
#endif
      LUMAT_AVAIL          = .FALSE.
#if ! defined(NOANA_BLKAUTO)
      IF (KEEP(13).EQ.2.AND.OPTION.EQ.3) THEN
        SIZE_MAPCOL = NBLK
        ALLOCATE(MAPCOL(SIZE_MAPCOL), STAT=allocok) 
        IF (allocok.NE.0) THEN
             INFO(1)  = -7
             INFO(2)  = NBLK
             IF ( LPOK ) THEN
                WRITE(LP, *) " ERROR allocate MAPCOL of size", 
     &           INFO(2)
           END IF
        ENDIF
        CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
        IF ( INFO(1) .LT. 0 ) GOTO 500
        MAPCOL_PTR => MAPCOL(1:SIZE_MAPCOL)
        MAPCOLallocated = .TRUE.
      ELSE
         SIZE_MAPCOL = 1
         MAPCOL_PTR => IDUMMY_TARGET
      ENDIF
#endif
      IF (KEEP(14).EQ.0.OR.KEEP(13).NE.2) THEN
         CALL MPI_BCAST( DOF2BLOCK, NDOF, MPI_INTEGER, MASTER, 
     &     COMM, IERR )
      ENDIF
#if ! defined(NOANA_BLKAUTO)
#endif
      CALL MUMPS_AB_COORD_TO_LMAT (  MYID, 
     &          NBLK, NDOF, NNZ, IRN, JCN,
     &          DOF2BLOCK, 
     &          INFO(1), INFO(2), LP, LPOK, 
     &          LMAT, IDUMMY8, KEEP
#if ! defined(NOANA_BLKAUTO)
     &          ,OPTION, MAPCOL_PTR(1), SIZE_MAPCOL 
#endif
     &  )
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
#if ! defined(NOANA_BLKAUTO)
      IF (OPTION.EQ.3) THEN
       CALL MUMPS_AB_CHECK_FINALIZE_MAPCOL ( 
     &     MYID, COMM, NPROCS, NPROCS,
     &     MAPCOL_PTR(1), SIZE_MAPCOL, 
     &     MAPCOL_VALIDATED, CONTIGUOUS_COL,
     &     LP, LPOK, INFO(1) )
         IF (.NOT.MAPCOL_VALIDATED) THEN
           IF (KEEP(337).EQ.3) THEN
             OPTION = 1
           ELSE
             OPTION = KEEP(337)
           ENDIF
         ENDIF
       ENDIF
      IF ( INFO(1) .LT. 0 ) GOTO 500
#endif
      IF (.NOT.MAPCOLallocated) THEN
        ALLOCATE(MAPCOL(NBLK), STAT=allocok) 
        IF (allocok.NE.0) THEN
             INFO(1)  = -7
             INFO(2)  = NBLK
             IF ( LPOK ) THEN
                WRITE(LP, *) " ERROR allocate MAPCOL of size", 
     &           INFO(2)
           END IF
        ENDIF
        CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
        IF ( INFO(1) .LT. 0 ) GOTO 500
        MAPCOLallocated = .TRUE.
      ENDIF
      IF ( MAPCOLallocated
#if ! defined(NOANA_BLKAUTO)
     &    .AND. (.NOT.MAPCOL_VALIDATED) 
#endif
     &   ) THEN 
        CALL MUMPS_AB_COL_DISTRIBUTION (OPTION, 
     &     INFO, ICNTL, COMM, NBLK, MYID, NPROCS,
     &     LMAT, MAPCOL, CONTIGUOUS_COL)
        CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
        IF ( INFO(1) .LT. 0 ) GOTO 500
      ENDIF
      CALL MUMPS_AB_BUILD_DCLEAN_LUMATRIX (
     &    MAPCOLonLUMAT, MAPCOL_IN_NSTEPS, CONTIGUOUS_COL,
     &    INFO, ICNTL, KEEP, COMM,  MYID, NBLK, NPROCS,
     &    LMAT, MAPCOL, NBLK, 
     &    IDUMMY_ARRAY(1), 
     &    1, 
     &    LUMAT)
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      LUMAT_AVAIL = .TRUE.
#if ! defined(NOANA_BLKAUTO)
      IF (KEEP(13).EQ.2) THEN
         IF (.NOT.LUMAT%CONTIGUOUS_COL) REMAP_CONTIG=.TRUE.
         K244= KEEP(244)
         CALL MUMPS_AB_DLUMAT_TO_BLK_DLUMAT ( 
     &        INFO, ICNTL, KEEP, COMM, MYID, 
     &        NPROCS, NDOF, LUMAT, LP, LPOK,
     &        DOF2BLOCK(1), NBLK, BLKFORMAT_DETECTED,
     &        MAPCOL_PTR(1), SIZE_MAPCOL,
     &        REMAP_CONTIG, K244,
     &        NPROCS, IDUMMY_ARRAY(1), 1
     &        )
         KEEP(244)=K244
         CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &        COMM, MYID )
         IF ( INFO(1) .LT. 0 ) GOTO 500
      ENDIF
#endif
       IF (allocated(MAPCOL))  DEALLOCATE(MAPCOL)
       CALL MUMPS_AB_LMAT_TO_CLEAN_G ( MYID, .FALSE., 
     &         READY_FOR_ANA_F,
     &         LUMAT, GCOMP, INFO, ICNTL, MEMCNT
     &      )
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      IF (KEEP(494).EQ.0.OR.KEEP(487).NE.1) THEN
       IF (LUMAT_AVAIL) THEN
         CALL MUMPS_AB_FREE_LMAT(LUMAT, KEEP(147))
         LUMAT_AVAIL=.FALSE.
       ENDIF
      ENDIF
      GOTO 600
  500 CONTINUE
      IF (allocated(MAPCOL))  DEALLOCATE(MAPCOL)
      CALL MUMPS_AB_FREE_LMAT(LMAT, KEEP(147))
      CALL MUMPS_AB_FREE_LMAT(LUMAT, KEEP(147))
      LUMAT_AVAIL=.FALSE.
  600 CONTINUE
      RETURN
      END SUBROUTINE MUMPS_AB_DCOORD_TO_DCOMPG
      SUBROUTINE MUMPS_AB_DCOORD_TO_DLUMAT (
     &     MYID, NPROCS, COMM, 
     &     NPROCS_PARAORD, PARAORD_to_idCOMM,
     &     NBLK, NDOF, NNZ, 
     &     IRN, JCN, DOF2BLOCK, 
     &     ICNTL, INFO, KEEP, KEEP8,
     &     LUMAT, LUMAT_AVAIL
#if ! defined(NOANA_BLKAUTO)
     &     , BLKFORMAT_DETECTED, CENTRALIZED_MAT
#endif
     &     )
#if ! defined(NOANA_BLKAUTO)
#endif
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, MASTER
      PARAMETER( MASTER = 0 )
      INTEGER, INTENT(IN)    :: MYID, NPROCS, NDOF
      INTEGER, INTENT(INOUT) :: NBLK
      INTEGER(8), INTENT(IN) :: NNZ
      INTEGER, INTENT(IN)    :: NPROCS_PARAORD, 
     &                          PARAORD_to_idCOMM(NPROCS_PARAORD)
      INTEGER, INTENT(IN)    :: IRN(max(1_8,NNZ)), JCN(max(1_8,NNZ))
      INTEGER, INTENT(INOUT) :: DOF2BLOCK(NDOF)
      INTEGER, INTENT(IN)    :: ICNTL(60), COMM
      INTEGER, INTENT(INOUT) :: KEEP(500), INFO(80)
      INTEGER(8), INTENT(INOUT) :: KEEP8(150)
#if ! defined(NOANA_BLKAUTO)
      LOGICAL, INTENT(IN)       :: CENTRALIZED_MAT
      LOGICAL, INTENT(INOUT)    :: BLKFORMAT_DETECTED
#endif
      TYPE(LMATRIX_T)        :: LUMAT  
      LOGICAL, INTENT(OUT)   :: LUMAT_AVAIL
      TYPE(LMATRIX_T)        :: LMAT
      INTEGER :: IDUMMY_ARRAY(1)
      INTEGER :: allocok, LP, MPG, I
      LOGICAL :: LPOK, PROKG
      INTEGER(8) :: OFFDIAG
      INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: MAPCOL
#if ! defined(NOANA_BLKAUTO)
      INTEGER, DIMENSION(:), POINTER             :: MAPCOL_PTR
      INTEGER :: SIZE_MAPCOL, IP, IPP, K244
      INTEGER, TARGET :: IDUMMY_TARGET(1)
      LOGICAL :: MAPCOL_VALIDATED 
      INTEGER :: NBCOLonMYID(0:NPROCS-1)
      INTEGER :: FIRSTCOLonMYID(0:NPROCS-1)
      LOGICAL :: REMAP_CONTIG
#endif
      LOGICAL :: MAPCOLonLUMAT, MAPCOL_IN_NSTEPS, CONTIGUOUS_COL,
     &           MAPCOLallocated
      INTEGER OPTION
      LP  = ICNTL( 1 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      MPG = ICNTL( 3 )
      PROKG = ( MPG .GT. 0 .and. MYID .eq. MASTER )
      IF (KEEP(13).EQ.2) THEN 
         OPTION= 3
      ELSE
        OPTION =KEEP(337) 
        IF (KEEP(72).EQ.1 .and. NBLK.LE.10) THEN
         OPTION =1
        ENDIF
      ENDIF
      MAPCOLonLUMAT        = .FALSE.
      MAPCOL_IN_NSTEPS     = .FALSE.
      CONTIGUOUS_COL       = .TRUE.
      MAPCOLallocated      = .FALSE.
#if ! defined(NOANA_BLKAUTO)
      REMAP_CONTIG         = .FALSE.
      MAPCOL_VALIDATED     = .FALSE.
#endif
      LUMAT_AVAIL          = .FALSE.
#if ! defined(NOANA_BLKAUTO)
      IF (KEEP(13).EQ.2.AND.OPTION.EQ.3) THEN
        SIZE_MAPCOL = NBLK
        ALLOCATE(MAPCOL(SIZE_MAPCOL), STAT=allocok) 
        IF (allocok.NE.0) THEN
             INFO(1)  = -7
             INFO(2)  = NBLK
             IF ( LPOK ) THEN
                WRITE(LP, *) " ERROR allocate MAPCOL of size", 
     &           INFO(2)
           END IF
        ENDIF
        CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
        IF ( INFO(1) .LT. 0 ) GOTO 500
        MAPCOL_PTR => MAPCOL(1:SIZE_MAPCOL)
        MAPCOLallocated = .TRUE.
      ELSE
         SIZE_MAPCOL = 1
         MAPCOL_PTR => IDUMMY_TARGET
      ENDIF
#endif
      IF (KEEP(14).EQ.0.OR.KEEP(13).NE.2) THEN
         CALL MPI_BCAST( DOF2BLOCK, NDOF, MPI_INTEGER, MASTER, 
     &     COMM, IERR )
      ENDIF
      CALL MUMPS_AB_COORD_TO_LMAT (  MYID, 
     &          NBLK, NDOF, NNZ, IRN, JCN,
     &          DOF2BLOCK, 
     &          INFO(1), INFO(2), LP, LPOK, 
     &          LMAT, OFFDIAG,  KEEP
#if ! defined(NOANA_BLKAUTO)
     &          ,OPTION, MAPCOL_PTR(1), SIZE_MAPCOL 
#endif
     &          )
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
#if ! defined(NOANA_BLKAUTO)
      IF (OPTION.EQ.3) THEN
       IF (NPROCS_PARAORD.LT.NPROCS) THEN
         REMAP_CONTIG=.TRUE.
       ENDIF
       IF (CENTRALIZED_MAT) THEN
        MAPCOL_VALIDATED = .TRUE.
        DO I=1, SIZE_MAPCOL
          MAPCOL_PTR(I) = MASTER 
        ENDDO
       ELSE
         CALL MUMPS_AB_CHECK_FINALIZE_MAPCOL ( 
     &     MYID, COMM, NPROCS,  NPROCS, 
     &     MAPCOL_PTR(1), SIZE_MAPCOL,  
     &     MAPCOL_VALIDATED,  CONTIGUOUS_COL,
     &     LP, LPOK, INFO(1) )
       ENDIF
       IF (.NOT.MAPCOL_VALIDATED) THEN
         IF (KEEP(337).EQ.3) THEN
           OPTION = 1
         ELSE
           OPTION = KEEP(337)
         ENDIF
       ELSE
         IF (CENTRALIZED_MAT) REMAP_CONTIG = .TRUE.
       ENDIF
      ENDIF
      IF ( INFO(1) .LT. 0 ) GOTO 500
#endif
      CALL MPI_ALLREDUCE (OFFDIAG, KEEP8(127), 1, MPI_INTEGER8,
     &     MPI_SUM, COMM, IERR)
      KEEP8(127) = KEEP8(127)+3_8*int(NDOF,8)
      KEEP8(126) = KEEP8(127)-2_8*int(NDOF,8)
      IF (.NOT.MAPCOLallocated) THEN
        ALLOCATE(MAPCOL(NBLK), STAT=allocok) 
        IF (allocok.NE.0) THEN
             INFO(1)  = -7
             INFO(2)  = NBLK
             IF ( LPOK ) THEN
                WRITE(LP, *) " ERROR allocate MAPCOL of size", 
     &           INFO(2)
             END IF
        ENDIF
        CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &       COMM, MYID )
        IF ( INFO(1) .LT. 0 ) GOTO 500
        MAPCOLallocated = .TRUE.
      ENDIF
      IF (.NOT.MAPCOLonLUMAT
#if ! defined(NOANA_BLKAUTO)
     &    .AND. (.NOT.MAPCOL_VALIDATED) 
#endif
     &    ) THEN 
        CALL MUMPS_AB_COL_DISTRIBUTION (OPTION, 
     &     INFO, ICNTL, COMM, NBLK, MYID, NPROCS_PARAORD, 
     &     LMAT, MAPCOL, CONTIGUOUS_COL)
        CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
        IF ( INFO(1) .LT. 0 ) GOTO 500
      ENDIF
      IF (NPROCS_PARAORD.LE.NPROCS
#if ! defined(NOANA_BLKAUTO)
     &    .AND..NOT.MAPCOL_VALIDATED
#endif
     &   ) THEN
         DO I=1, NBLK
           MAPCOL(I) = PARAORD_to_idCOMM(MAPCOL(I)+1)
         ENDDO
      ENDIF
#if ! defined(NOANA_BLKAUTO)
      IF (KEEP(13).EQ.2.AND.OPTION.EQ.3
     &    .AND..NOT.CENTRALIZED_MAT.AND.CONTIGUOUS_COL) THEN
       IF (NPROCS_PARAORD.GT.NBLK) THEN
           INFO(1)=-1981
       ELSE
         NBCOLonMYID    = 0
         FIRSTCOLonMYID = 0
         DO I= 1, NBLK
          IP = MAPCOL(I)
          NBCOLonMYID(IP) = NBCOLonMYID(IP)+1
          IF (FIRSTCOLonMYID(IP).EQ.0) THEN
             FIRSTCOLonMYID(IP)=I 
          ENDIF
         ENDDO
         DO IP=0, NPROCS_PARAORD-1
           IF (NBCOLonMYID(IP).EQ.0) THEN
             DO IPP=0, NPROCS_PARAORD-1
               IF (NBCOLonMYID(IPP).GT.1) THEN
                 MAPCOL(IP)      = FIRSTCOLonMYID(IPP)
                 FIRSTCOLonMYID(IPP) =  FIRSTCOLonMYID(IPP)+1
                 NBCOLonMYID(IPP) = NBCOLonMYID(IPP) -1
                 EXIT
               ENDIF
             ENDDO
           ENDIF
         ENDDO
       ENDIF
       CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
        IF ( INFO(1) .LT. 0 ) GOTO 500
      ENDIF
#endif
      CALL MUMPS_AB_BUILD_DCLEAN_LUMATRIX (
     &    MAPCOLonLUMAT, MAPCOL_IN_NSTEPS, CONTIGUOUS_COL,
     &    INFO, ICNTL, KEEP, COMM,  MYID, NBLK, NPROCS,
     &    LMAT, MAPCOL, NBLK, 
     &    IDUMMY_ARRAY, 1, 
     &    LUMAT)
      LUMAT_AVAIL = .TRUE.
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
#if ! defined(NOANA_BLKAUTO)
      IF (KEEP(13).EQ.2) THEN
         IF (.NOT.LUMAT%CONTIGUOUS_COL) REMAP_CONTIG=.TRUE.
         K244= KEEP(244)
         CALL MUMPS_AB_DLUMAT_TO_BLK_DLUMAT ( 
     &        INFO, ICNTL, KEEP, COMM, MYID, 
     &        NPROCS, NDOF, LUMAT, LP, LPOK,
     &        DOF2BLOCK(1), NBLK, BLKFORMAT_DETECTED, 
     &        MAPCOL_PTR(1), SIZE_MAPCOL, 
     &        REMAP_CONTIG, K244, 
     &        NPROCS_PARAORD, PARAORD_to_idCOMM, 
     &        NPROCS_PARAORD
     &        )
         KEEP(244)=K244
      ENDIF
#endif
      IF (allocated(MAPCOL))  DEALLOCATE(MAPCOL)
      GOTO 600
  500 CONTINUE
      IF (allocated(MAPCOL))  DEALLOCATE(MAPCOL)
      CALL MUMPS_AB_FREE_LMAT(LMAT, KEEP(147))
      CALL MUMPS_AB_FREE_LMAT(LUMAT, KEEP(147))
  600 CONTINUE
      IF (allocated(MAPCOL))  DEALLOCATE(MAPCOL)
      RETURN
      END SUBROUTINE MUMPS_AB_DCOORD_TO_DLUMAT
      SUBROUTINE MUMPS_AB_DCOORD_TO_DTREE_LUMAT ( 
     &     MYID, NPROCS, COMM, 
     &     NBLK, NDOF, NNZ, 
     &     IRN, JCN, 
     &     PROCNODE_STEPS, NSTEPS, STEP,
     &     ICNTL, INFO, KEEP, 
     &     MAPCOL, LUMAT)
      USE MUMPS_ANA_BLK_M, ONLY: LMATRIX_T
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, MASTER
      PARAMETER( MASTER = 0 )
      INTEGER, INTENT(IN)    :: MYID, NPROCS, NBLK, NDOF, NSTEPS
      INTEGER(8), INTENT(IN) :: NNZ
      INTEGER, INTENT(IN)    :: IRN(max(1_8,NNZ)), JCN(max(1_8,NNZ))
      INTEGER, INTENT(IN)    :: ICNTL(60), COMM
      INTEGER, INTENT(IN)    :: PROCNODE_STEPS(NSTEPS)
      INTEGER, INTENT(IN)    :: STEP(NBLK)
      INTEGER, INTENT(INOUT) :: KEEP(500), INFO(80)
      INTEGER, INTENT(OUT)   :: MAPCOL(NSTEPS)
      TYPE(LMATRIX_T)        :: LUMAT  
      INTEGER, DIMENSION(:), allocatable:: DOF2BLOCK
      TYPE(LMATRIX_T)        :: LMAT
      INTEGER :: allocok, LP
      LOGICAL :: LPOK
      INTEGER :: IDOF, ISTEP
      LOGICAL :: MAPCOL_IN_NSTEPS, MAPCOLonLUMAT, CONTIGUOUS_COL
      INTEGER(8) :: IDUMMY8
#if ! defined(NOANA_BLKAUTO)
      INTEGER  :: IDUMMY_ARRAY(1)
#endif
      INTEGER OPTION
      PARAMETER (OPTION=2)
      INTEGER MUMPS_PROCNODE
      LP  = ICNTL( 1 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      MAPCOLonLUMAT    = .FALSE.
      MAPCOL_IN_NSTEPS = .TRUE.
      IF (KEEP(14).EQ.1) THEN
         CALL MUMPS_ABORT()
      ENDIF
      allocate(DOF2BLOCK(NDOF), STAT=allocok)
      IF (allocok.NE.0) THEN
           INFO( 1 ) = -7
           INFO( 2 ) = NDOF
           IF ( LPOK ) WRITE(LP, 150) ' DOF2BLOCK'
      ENDIF
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      DO IDOF=1, NDOF
       DOF2BLOCK(IDOF) = IDOF
      ENDDO
      CALL MUMPS_AB_COORD_TO_LMAT (  MYID, 
     &          NBLK, NDOF, NNZ, IRN, JCN,
     &          DOF2BLOCK, 
     &          INFO(1), INFO(2), LP, LPOK, 
     &          LMAT, IDUMMY8, KEEP
#if ! defined(NOANA_BLKAUTO)
     &          , 0, IDUMMY_ARRAY, 1
#endif
     &          )
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      IF (allocated(DOF2BLOCK))  DEALLOCATE(DOF2BLOCK)
      IF (MYID.EQ.MASTER) THEN
       DO ISTEP=1, NSTEPS
        MAPCOL(ISTEP) = 
     &              MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),KEEP(199))
       ENDDO
      ENDIF
      CALL MPI_BCAST( MAPCOL(1), NSTEPS, MPI_INTEGER,
     &     MASTER, COMM, IERR )
      CALL MPI_BCAST( STEP(1), NBLK, MPI_INTEGER,
     &     MASTER, COMM, IERR )
      CONTIGUOUS_COL = .FALSE.
      CALL MUMPS_AB_BUILD_DCLEAN_LUMATRIX(
     &    MAPCOLonLUMAT, MAPCOL_IN_NSTEPS, CONTIGUOUS_COL,
     &    INFO, ICNTL, KEEP, COMM,  MYID, NBLK, NPROCS,
     &    LMAT, MAPCOL, NSTEPS,
     &    STEP, NBLK, LUMAT)
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      GOTO 600
  500 CONTINUE
      IF (allocated(DOF2BLOCK))  DEALLOCATE(DOF2BLOCK)
      CALL MUMPS_AB_FREE_LMAT(LMAT, KEEP(147))
      CALL MUMPS_AB_FREE_LMAT(LUMAT, KEEP(147))
  600 CONTINUE
      RETURN
 150  FORMAT(
     & /' ** FAILURE IN MUMPS_AB_DCOORD_TO_DTREE_LUMAT, ', 
     &  ' DYNAMIC ALLOCATION OF ',
     &     A30)
      END SUBROUTINE MUMPS_AB_DCOORD_TO_DTREE_LUMAT
      SUBROUTINE MUMPS_AB_DIST_LMAT_TO_LUMAT ( 
     &  UNFOLD,
     &  MAPCOL_IN_NSTEPS,
     &  INFO, ICNTL, COMM, MYID, NBLK, SLAVEF,
     &  LMAT, MAPCOL, SIZE_MAPCOL, STEP, SIZE_STEP, 
     &  LUMAT, NBRECORDS, NSEND8, NLOCAL8
     &  )
#if ! defined(NOANA_BLKAUTO)
#endif
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER  :: IERR, MASTER, MSGSOU
      PARAMETER (MASTER=0)
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      LOGICAL, INTENT(IN) :: UNFOLD, MAPCOL_IN_NSTEPS
      INTEGER, INTENT(IN) :: MYID, SLAVEF, NBLK 
      INTEGER, INTENT(IN) :: SIZE_MAPCOL, SIZE_STEP 
      INTEGER, INTENT(IN) :: ICNTL(60), COMM, NBRECORDS
      INTEGER             :: INFO(80)
      TYPE(LMATRIX_T), INTENT(IN)  :: LMAT
      INTEGER, INTENT(IN)          :: MAPCOL(SIZE_MAPCOL)
      INTEGER, INTENT(IN)          :: STEP(SIZE_STEP)
      TYPE(LMATRIX_T), INTENT(INOUT) :: LUMAT
      INTEGER(8), INTENT(OUT) :: NSEND8, NLOCAL8
      INTEGER :: LP, MP, allocok
      INTEGER :: IB, JB, I, II, ISEND, JSEND, ITOSEND
      LOGICAL :: LPOK 
      INTEGER :: NBTOSEND
      INTEGER :: LMAT_FIRST, LUMAT_FIRST
      INTEGER END_MSG_2_RECV
      INTEGER KPROBE, FREQPROBE
      INTEGER                            :: SIZE_PTLOC
      INTEGER, ALLOCATABLE, DIMENSION(:) :: PTLOC
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: BUFI
      INTEGER, ALLOCATABLE, DIMENSION(:) :: BUFRECI
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IACT, IREQI
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: SEND_ACTIVE
      INTEGER  :: DEST
      LOGICAL  :: FLAG
      LP  = ICNTL( 1 )
      MP  = ICNTL( 2 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      IF (UNFOLD) THEN
        NBTOSEND = 2
      ELSE
        NBTOSEND = 1
      ENDIF
      NSEND8  = 0_8
      NLOCAL8 = 0_8
      END_MSG_2_RECV = SLAVEF-1
      ALLOCATE( IACT(SLAVEF), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating IACT in matrix distribution'
        END IF
        INFO(1) = -7
        INFO(2) = SLAVEF
        GOTO 20
      END IF
      ALLOCATE( IREQI(SLAVEF), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating IREQI in matrix distribution'
        END IF
        INFO(1) = -7
        INFO(2) = SLAVEF
        GOTO 20
      END IF
      ALLOCATE( SEND_ACTIVE(SLAVEF), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating SEND_ACTIVE in matrix distribution'
        END IF
        INFO(1) = -7
        INFO(2) = SLAVEF
        GOTO 20
      END IF
      ALLOCATE( BUFI( NBRECORDS * 2 + 1, 2, SLAVEF ), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating int buffer for matrix distribution'
        END IF
        INFO(1) = -7
        INFO(2) = ( NBRECORDS * 2 + 1 ) * SLAVEF * 2
        GOTO 20
      END IF
      ALLOCATE( BUFRECI( NBRECORDS * 2 + 1 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &    '** Error allocating int recv buffer for matrix distribution'
        END IF
        INFO(1) = -7
        INFO(2) = NBRECORDS * 2 + 1
        GOTO 20
      END IF
      SIZE_PTLOC =  max(LUMAT%NBCOL_LOC,1)
      ALLOCATE( PTLOC( SIZE_PTLOC ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &    '** Error allocating int recv buffer for matrix distribution'
        END IF
        INFO(1) = -7
        INFO(2) = NBLK
        GOTO 20
      END IF
 20   CONTINUE
      CALL MUMPS_PROPINFO( ICNTL, INFO, COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 100
      DO I = 1, SLAVEF
        BUFI( 1, 1, I ) = 0
      END DO
      DO I = 1, SLAVEF
        BUFI( 1, 2, I ) = 0
      END DO
      DO I = 1, SLAVEF
        SEND_ACTIVE( I ) = .FALSE.
        IACT( I ) = 1
      END DO
      DO I = 1, SIZE_PTLOC
        PTLOC(I) = 0
      END DO
      KPROBE = 0
      FREQPROBE = max(1,NBRECORDS/10)
      IF (SLAVEF .EQ. 1) FREQPROBE = huge(FREQPROBE)
      LUMAT_FIRST = LUMAT%FIRST
      LMAT_FIRST  = LMAT%FIRST
      DO JB=1, LMAT%NBCOL_LOC
       IF ( LMAT%COL(JB)%NBINCOL.EQ.0) CYCLE
       DO II=1,  LMAT%COL(JB)%NBINCOL
        KPROBE = KPROBE + 1
        IF ( KPROBE .eq. FREQPROBE ) THEN
          KPROBE = 0
          CALL MPI_IPROBE( MPI_ANY_SOURCE, LMATDIST, COMM,
     &                     FLAG, STATUS, IERR )
          IF ( FLAG ) THEN
            MSGSOU = STATUS( MPI_SOURCE )
            CALL MPI_RECV( BUFRECI(1), NBRECORDS * 2 + 1, 
     &                 MPI_INTEGER,
     &                 MSGSOU, LMATDIST, COMM, STATUS, IERR )
            CALL MUMPS_AB_LMAT_TREAT_RECV_BUF(
     &             MYID, BUFRECI(1), NBRECORDS, LUMAT, 
     &             PTLOC(1), SIZE_PTLOC, END_MSG_2_RECV
     &             )
          END IF
        END IF
        IB = LMAT%COL(JB)%IRN(II)
        DO ITOSEND=1,NBTOSEND
         IF (ITOSEND.EQ.1) THEN
          IF (MAPCOL_IN_NSTEPS) THEN
           DEST  = MAPCOL(abs(STEP(JB+LMAT_FIRST-1)))
          ELSE
           DEST  = MAPCOL(JB+LMAT_FIRST-1)
          ENDIF
          ISEND = IB
          JSEND = JB+LMAT_FIRST-1
         ELSE
          IF (MAPCOL_IN_NSTEPS) THEN
           DEST  = MAPCOL(abs(STEP(IB)))
          ELSE
           DEST  = MAPCOL(IB)
          ENDIF
          ISEND = JB+LMAT_FIRST-1
          JSEND = IB
         ENDIF
         IF (DEST.EQ.MYID) THEN
          LUMAT%COL(JSEND-LUMAT_FIRST+1)
     &                     %IRN(1+PTLOC(JSEND-LUMAT_FIRST+1))= ISEND
          PTLOC(JSEND-LUMAT_FIRST+1)    = PTLOC(JSEND-LUMAT_FIRST+1) + 1
          NLOCAL8 = NLOCAL8 + 1_8
         ELSE
          NSEND8 = NSEND8 + 1_8
          CALL MUMPS_AB_LMAT_FILL_BUFFER( 
     &     DEST, ISEND, JSEND, 
     &     BUFI, BUFRECI, PTLOC, SIZE_PTLOC,
     &     NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, 
     &     SEND_ACTIVE, LUMAT, END_MSG_2_RECV
     &      )
         ENDIF
        ENDDO
       ENDDO
      ENDDO
      DEST = -3
      CALL MUMPS_AB_LMAT_FILL_BUFFER(DEST, ISEND, JSEND, 
     &     BUFI, BUFRECI, PTLOC, SIZE_PTLOC,
     &     NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, 
     &     SEND_ACTIVE, LUMAT, END_MSG_2_RECV
     &      )
      DO WHILE ( END_MSG_2_RECV .NE. 0 )
        CALL MPI_RECV( BUFRECI(1), NBRECORDS * 2 + 1, MPI_INTEGER,
     &                 MPI_ANY_SOURCE, LMATDIST, COMM, STATUS, IERR )
        CALL MUMPS_AB_LMAT_TREAT_RECV_BUF(
     &             MYID, BUFRECI(1), NBRECORDS, LUMAT, 
     &             PTLOC(1), SIZE_PTLOC, END_MSG_2_RECV
     &             )
      END DO
      DO I = 1, SLAVEF
        IF ( SEND_ACTIVE( I ) ) THEN
          CALL MPI_WAIT( IREQI( I ), STATUS, IERR )
        END IF
      END DO
 100  CONTINUE
      IF (ALLOCATED(PTLOC))   DEALLOCATE( PTLOC )
      IF (ALLOCATED(BUFI))    DEALLOCATE( BUFI )
      IF (ALLOCATED(BUFRECI)) DEALLOCATE( BUFRECI )
      IF (ALLOCATED(IACT))    DEALLOCATE( IACT )
      IF (ALLOCATED(IREQI))   DEALLOCATE( IREQI )
      IF (ALLOCATED(SEND_ACTIVE)) DEALLOCATE( SEND_ACTIVE )
      RETURN
      END SUBROUTINE MUMPS_AB_DIST_LMAT_TO_LUMAT
      SUBROUTINE MUMPS_AB_LMAT_TREAT_RECV_BUF (
     &   MYID, BUFI, NBRECORDS, LUMAT,  
     &   PTLOC, SIZE_PTLOC, END_MSG_2_RECV
     & )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER, INTENT(IN)   :: MYID, NBRECORDS, SIZE_PTLOC
      INTEGER, INTENT(IN)   :: BUFI( NBRECORDS * 2 + 1 ) 
      INTEGER, INTENT(INOUT):: END_MSG_2_RECV, PTLOC(SIZE_PTLOC)
      TYPE(LMATRIX_T), INTENT(INOUT) :: LUMAT
      INTEGER :: IREC, NB_REC, IB, JB, LUMAT_FIRST
      NB_REC = BUFI( 1 )
      IF ( NB_REC .LE. 0 ) THEN
        END_MSG_2_RECV = END_MSG_2_RECV - 1
        NB_REC = - NB_REC
      END IF
      IF ( NB_REC .eq. 0 ) RETURN
      LUMAT_FIRST = LUMAT%FIRST
      DO IREC = 1, NB_REC
       IB = BUFI( IREC * 2 )
       JB = BUFI( IREC * 2 + 1 )
       JB = JB - LUMAT_FIRST + 1
       LUMAT%COL(JB)%IRN(1+PTLOC(JB))= IB
       PTLOC(JB) = PTLOC(JB) + 1
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_AB_LMAT_TREAT_RECV_BUF
      SUBROUTINE MUMPS_AB_LMAT_FILL_BUFFER (
     &     DEST, ISEND, JSEND, 
     &     BUFI, BUFRECI, PTLOC, SIZE_PTLOC,
     &     NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, 
     &     SEND_ACTIVE, LUMAT, END_MSG_2_RECV
     &      )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER, INTENT(IN)  :: DEST, ISEND, JSEND, SLAVEF, COMM, MYID, 
     &                        NBRECORDS, SIZE_PTLOC
      INTEGER, INTENT(INOUT) :: END_MSG_2_RECV, PTLOC(SIZE_PTLOC)
      TYPE(LMATRIX_T), INTENT(INOUT)  :: LUMAT
      LOGICAL, INTENT(INOUT) ::  SEND_ACTIVE(SLAVEF)
      INTEGER, INTENT(INOUT) ::  IREQI(SLAVEF), IACT(SLAVEF)
      INTEGER, INTENT(INOUT) ::  BUFI( NBRECORDS * 2 + 1, 2, SLAVEF )
      INTEGER, INTENT(INOUT) ::  BUFRECI( NBRECORDS * 2 + 1)
      INTEGER :: IBEG, IEND, ISLAVE, TAILLE_SEND_I, IREQ, MSGSOU,
     &           NBREC, IERR
      LOGICAL :: FLAG
      IF ( DEST .eq. -3 ) THEN
        IBEG = 1
        IEND = SLAVEF
      ELSE
        IBEG = DEST + 1
        IEND = DEST + 1
      END IF
      DO ISLAVE = IBEG, IEND
        NBREC = BUFI(1,IACT(ISLAVE),ISLAVE)
        IF ( DEST .eq. -3 ) THEN
          BUFI(1,IACT(ISLAVE),ISLAVE) = - NBREC
        END IF
        IF ( DEST .eq. -3 .or. NBREC + 1 > NBRECORDS ) THEN
          DO WHILE ( SEND_ACTIVE( ISLAVE ) )
            CALL MPI_TEST( IREQI( ISLAVE ), FLAG, STATUS, IERR )
            IF ( .NOT. FLAG ) THEN
                CALL MPI_IPROBE( MPI_ANY_SOURCE, LMATDIST, COMM,
     &                           FLAG, STATUS, IERR )
                IF ( FLAG ) THEN
                  MSGSOU = STATUS(MPI_SOURCE)
                  CALL MPI_RECV( BUFRECI(1), 2*NBRECORDS+1,
     &                  MPI_INTEGER, MSGSOU, LMATDIST, COMM,
     &                  STATUS, IERR )
                  CALL MUMPS_AB_LMAT_TREAT_RECV_BUF(
     &             MYID, BUFRECI, NBRECORDS, LUMAT, 
     &             PTLOC(1), SIZE_PTLOC, END_MSG_2_RECV
     &             )
                END IF
            ELSE
              SEND_ACTIVE( ISLAVE ) = .FALSE.
            END IF
          END DO
          IF ( ISLAVE - 1 .ne. MYID ) THEN
            TAILLE_SEND_I = NBREC * 2 + 1
            CALL MPI_ISEND( BUFI(1, IACT(ISLAVE), ISLAVE ),
     &          TAILLE_SEND_I,
     &          MPI_INTEGER, ISLAVE - 1, LMATDIST, COMM,
     &          IREQI( ISLAVE ), IERR )
            SEND_ACTIVE( ISLAVE ) = .TRUE.
          ELSE
            IF (NBREC.NE.0) THEN 
              write(*,*) " Internal error in ", 
     &                    " MUMPS_AB_LMAT_FILL_BUFFER "
              CALL MUMPS_ABORT()
            ENDIF
          END IF
          IACT( ISLAVE ) = 3 - IACT( ISLAVE )
          BUFI( 1, IACT( ISLAVE ), ISLAVE ) = 0
        END IF
        IF ( DEST .ne. -3 ) THEN
          IREQ = BUFI(1,IACT(ISLAVE),ISLAVE) + 1
          BUFI(1,IACT(ISLAVE),ISLAVE) = IREQ
          BUFI(IREQ*2,IACT(ISLAVE),ISLAVE)  = ISEND
          BUFI(IREQ*2+1,IACT(ISLAVE),ISLAVE) = JSEND
        END IF
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_AB_LMAT_FILL_BUFFER
      SUBROUTINE MUMPS_AB_GATHER_GRAPH (
     &    ICNTL, KEEP, COMM, MYID, NPROCS, INFO, 
     &    GCOMP_DIST, GCOMP)
      USE MUMPS_ANA_BLK_M, ONLY : COMPACT_GRAPH_T
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, MASTER
      PARAMETER( MASTER = 0 )
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      TYPE(COMPACT_GRAPH_T), INTENT(IN)  :: GCOMP_DIST
      INTEGER, INTENT(IN)    :: MYID, NPROCS,  ICNTL(60), COMM,
     &                          KEEP(500)
      INTEGER, INTENT(INOUT) ::  INFO(80)
      TYPE(COMPACT_GRAPH_T) :: GCOMP
       INTEGER    :: NG, allocok, LP, MPG, I, J, K, PI
       INTEGER    :: INDX, NB_BLOCK_SENT, MAX_NBBLOCK_loc, NRECV,
     &               BLOCKSIZE, SIZE_SENT, NB_BLOCKS, NBBLOCK_loc 
       INTEGER    :: JFIRST_GDIST, JLAST_GDIST, NG_LOCAL_GDIST, 
     &                NBNONEMPTY_PI
       INTEGER(4) :: IOVFLO
       INTEGER(8) :: NZG, NZG_CENT, I8, IBEG8, IEND8, 
     &               SIZEGCOMPALLOCATED
       LOGICAL    :: LPOK, PROKG
       INTEGER(8), ALLOCATABLE, DIMENSION(:) :: IQ
      INTEGER, ALLOCATABLE :: REQPTR(:)
      INTEGER(8), ALLOCATABLE :: GPTR(:), GPTR_cp(:)
      LP  = ICNTL( 1 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      MPG = ICNTL( 3 )
      PROKG = ( MPG .GT. 0 .and. MYID .eq. MASTER )
      PROKG = (PROKG.AND.(ICNTL(4).GE.2))
      IOVFLO = huge(IOVFLO)  
      BLOCKSIZE = int(max(100000_8,int(IOVFLO,8)/200_8))
      NZG =  GCOMP_DIST%NZG
      NG  =  GCOMP_DIST%NG
      JFIRST_GDIST    = GCOMP_DIST%FIRST
      JLAST_GDIST     = GCOMP_DIST%LAST
      NG_LOCAL_GDIST  = JLAST_GDIST - JFIRST_GDIST + 1
      CALL MPI_REDUCE( NZG, NZG_CENT, 1, MPI_INTEGER8, 
     &     MPI_SUM, MASTER, COMM, IERR )
      IF (MYID.EQ.MASTER) THEN
        GCOMP%NZG    = NZG_CENT
        GCOMP%NG     = NG
        GCOMP%FIRST  = 1
        GCOMP%LAST   = NG
        SIZEGCOMPALLOCATED     = NZG_CENT+int(NG,8)+1_8
        GCOMP%SIZEADJALLOCATED = SIZEGCOMPALLOCATED
        ALLOCATE( GCOMP%ADJ(SIZEGCOMPALLOCATED),
     &          GCOMP%IPE(NG+1), 
     &          GPTR( NPROCS ),
     &          GPTR_cp( NPROCS ),
     &          REQPTR( NPROCS-1 ),
     &          IQ(NG+1),STAT=allocok)
        IF (allocok.NE.0) THEN
           INFO( 1 ) = -7
           CALL MUMPS_SET_IERROR(
     &        NZG_CENT + 3_8*int(NG,8)+3_8+3_8*int(NPROCS,8)-1_8, 
     &        INFO(2))
           IF ( LPOK )
     &     WRITE(LP, *) " ERROR allocating graph in",
     &                  " MUMPS_AB_GATHER_GRAPH"
        ENDIF
      ELSE
        ALLOCATE( IQ(NG_LOCAL_GDIST), STAT=allocok)
        IF (allocok.NE.0) THEN
           INFO( 1 ) = -7
           INFO( 2 ) = NG+1
           IF ( LPOK )
     &     WRITE(LP, *) " ERROR allocating pointers",
     &                  " MUMPS_AB_GATHER_GRAPH"
        END IF
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF (INFO(1).LT.0) GOTO 500
      DO I=1,NG_LOCAL_GDIST
       IQ(I) = int(GCOMP_DIST%IPE(I+1)-GCOMP_DIST%IPE(I))
      ENDDO
      IF (MYID.EQ.MASTER) THEN
       DO J=1, NG
        GCOMP%IPE(J) = 0
       ENDDO
       IF (NG_LOCAL_GDIST.GT.0) THEN
        DO J=1, NG_LOCAL_GDIST
          GCOMP%IPE(J+JFIRST_GDIST-1) = IQ(J)
        ENDDO
       ENDIF
       DO PI = 1, NPROCS - 1
         CALL MPI_RECV( NBNONEMPTY_PI, 1, 
     &           MPI_INTEGER, PI,
     &           GATHERG_NB, COMM, STATUS, IERR )
         IF (NBNONEMPTY_PI.GT.0) THEN
           CALL MPI_RECV( J, 1, 
     &           MPI_INTEGER, PI,
     &           GATHERG_FIRST, COMM, STATUS, IERR )
           CALL MPI_RECV( GCOMP%IPE(J), NBNONEMPTY_PI, 
     &           MPI_INTEGER8, PI,
     &           GATHERG_IPE, COMM, STATUS, IERR )
         ENDIF
       ENDDO
      ELSE
        CALL MPI_SEND( NG_LOCAL_GDIST, 1, MPI_INTEGER, MASTER,
     &       GATHERG_NB, COMM, IERR )
        IF (NG_LOCAL_GDIST.GT.0) THEN
          CALL MPI_SEND( JFIRST_GDIST, 1, MPI_INTEGER, MASTER,
     &       GATHERG_FIRST, COMM, IERR )
          CALL MPI_SEND( IQ(1), NG_LOCAL_GDIST, 
     &       MPI_INTEGER8, MASTER,
     &       GATHERG_IPE, COMM, IERR )
        ENDIF
      ENDIF
      IF (MYID.EQ.MASTER) THEN
       IQ(1) = 1_8
       DO I=1,NG
         IQ(I+1) = IQ(I) + GCOMP%IPE(I)
         GCOMP%IPE(I) = IQ(I)
       ENDDO
       GCOMP%IPE(NG+1) = IQ(NG+1)
       DEALLOCATE(IQ)
      ELSE
       DEALLOCATE(IQ)
      ENDIF
      IF (MYID.EQ.MASTER) THEN
        NB_BLOCK_SENT = 0
        MAX_NBBLOCK_loc  = 0
        DO I = 1, NPROCS - 1
            CALL MPI_RECV( GPTR( I+1 ), 1, 
     &           MPI_INTEGER8, I,
     &           GATHERG_NZG, COMM, STATUS, IERR )
         NBBLOCK_loc = ceiling(dble(GPTR(I+1))/dble(BLOCKSIZE))
         MAX_NBBLOCK_loc = max(MAX_NBBLOCK_loc, NBBLOCK_loc)
         NB_BLOCK_SENT = NB_BLOCK_SENT + NBBLOCK_loc
        ENDDO
        GPTR( 1 ) = NZG + 1_8
         DO I = 2, NPROCS
            GPTR( I ) = GPTR( I ) + GPTR( I-1 )
         END DO
      ELSE
        CALL MPI_SEND( NZG, 1, MPI_INTEGER8, MASTER,
     &        GATHERG_NZG, COMM, IERR )
      ENDIF
      IF (MYID.EQ.MASTER) THEN
        DO I=1, NPROCS
         GPTR_cp(I) = GPTR(I) 
        ENDDO
        IF (NZG.GT.0_8) THEN
         DO I8=1, NZG
          GCOMP%ADJ(I8) = GCOMP_DIST%ADJ(I8)
         ENDDO
        ENDIF
        NB_BLOCKS = 0
        DO K = 1, MAX_NBBLOCK_loc
         NRECV = 0
         DO I = 1, NPROCS - 1
            IBEG8     = GPTR_cp( I )
            IF ( IBEG8 .LT. GPTR(I+1))  THEN
              NRECV = NRECV + 1
              IEND8 = min(IBEG8+int(BLOCKSIZE,8)-1_8, 
     &                    GPTR(I+1)-1_8)
              GPTR_cp( I ) = IEND8 + 1_8
              SIZE_SENT   = int(IEND8 -  IBEG8 + 1_8)
              NB_BLOCKS   = NB_BLOCKS + 1
              CALL MPI_IRECV( GCOMP%ADJ(IBEG8), SIZE_SENT, 
     &           MPI_INTEGER,
     &           I, GATHERG_ADJ, COMM, REQPTR(I), IERR )
            ELSE
             REQPTR( I ) = MPI_REQUEST_NULL
            ENDIF
         END DO
         DO I = 1, NRECV
             CALL MPI_WAITANY
     &           ( NPROCS-1, REQPTR, INDX, 
     &           STATUS, IERR )
         ENDDO
      END DO
        DEALLOCATE( REQPTR )
        DEALLOCATE( GPTR )
        DEALLOCATE( GPTR_cp )
      ELSE
        IF (NZG.EQ.0) GOTO 600
        DO I8=1_8, NZG, int(BLOCKSIZE,8)
         SIZE_SENT = BLOCKSIZE
         IF (NZG-I8+1_8.LT.int(BLOCKSIZE,8)) THEN
              SIZE_SENT = int(NZG-I8+1_8)
         ENDIF
         CALL MPI_SEND( 
     &           GCOMP_DIST%ADJ(I8), SIZE_SENT,
     &           MPI_INTEGER, MASTER,
     &           GATHERG_ADJ, COMM, IERR )
        ENDDO
      ENDIF
       GOTO 600
 500  CONTINUE
      IF (MYID.EQ.MASTER) THEN
       IF (associated(GCOMP%ADJ)) THEN
          DEALLOCATE(GCOMP%ADJ)
          nullify(GCOMP%ADJ)
       ENDIF
       IF (associated(GCOMP%IPE)) THEN
        DEALLOCATE(GCOMP%IPE)
        nullify(GCOMP%IPE)
       ENDIF
      ENDIF
 600  CONTINUE
      IF (allocated(IQ)) DEALLOCATE(IQ)
      RETURN
      END SUBROUTINE MUMPS_AB_GATHER_GRAPH
#if ! defined(NOANA_BLKAUTO)
      SUBROUTINE MUMPS_REDIST_LUMAT_TOCONTLUMAT (
     &       INFO, ICNTL, KEEP, COMM, MYID, NBLK,
     &       NPROCS, LUMAT, K147, 
     &       NPROCS_PARAORD, PARAORD_to_idCOMM, 
     &       SIZE_PARAORD_to_idCOMM
     &     )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T, COL_LMATRIX_T
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER, INTENT(IN)   :: MYID, NBLK, KEEP(500)
      INTEGER, INTENT(IN)   :: ICNTL(60), COMM, NPROCS
      INTEGER               :: INFO(80)
      TYPE(LMATRIX_T), INTENT(INOUT) :: LUMAT
      INTEGER, INTENT(IN)   :: K147
      INTEGER, INTENT(IN)  :: SIZE_PARAORD_to_idCOMM
      INTEGER, INTENT(IN)   :: NPROCS_PARAORD,
     &       PARAORD_to_idCOMM(SIZE_PARAORD_to_idCOMM)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAPCOL_CONT
      INTEGER :: IDUMMY_ARRAY(1), NBRECORDS, allocok
      TYPE(LMATRIX_T) :: LUMAT_REMAP
      INTEGER(8) :: NSEND8, NLOCAL8
      INTEGER*8  :: swap8
      INTEGER    :: swap
      LOGICAL    :: swapl
      TYPE(COL_LMATRIX_T), DIMENSION(:), POINTER :: PTR_COL
      ALLOCATE(MAPCOL_CONT(NBLK), stat=allocok)
      IF (allocok .ne.0) then
           INFO(1)  = -7
           INFO(2) =  NBLK
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &  COMM, MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
      CALL MUMPS_INI_REDIST_LUMAT ( 
     &     INFO, ICNTL, KEEP, COMM, MYID, NBLK,
     &     NPROCS, .FALSE.,  
     &     LUMAT, IDUMMY_ARRAY, 1, 
     &     MAPCOL_CONT, NBLK,
     &     LUMAT_REMAP, NBRECORDS, IDUMMY_ARRAY(1), 
     &     NPROCS_PARAORD, PARAORD_to_idCOMM, 
     &     SIZE_PARAORD_to_idCOMM
     &     )
      IF ( INFO(1).LT.0 ) GOTO 500
      CALL MUMPS_AB_DIST_LMAT_TO_LUMAT ( 
     &  .FALSE., 
     &  .FALSE.,  
     &  INFO, ICNTL, COMM, MYID, NBLK, NPROCS,
     &  LUMAT, MAPCOL_CONT, NBLK, IDUMMY_ARRAY, 1,
     &  LUMAT_REMAP, NBRECORDS, NSEND8, NLOCAL8
     & )
      swap                = LUMAT_REMAP%NBCOL 
      LUMAT_REMAP%NBCOL   = LUMAT%NBCOL
      LUMAT%NBCOL         = swap
      swap                  = LUMAT_REMAP%NBCOL_LOC 
      LUMAT_REMAP%NBCOL_LOC = LUMAT%NBCOL_LOC
      LUMAT%NBCOL_LOC       = swap
      swap                = LUMAT_REMAP%FIRST
      LUMAT_REMAP%FIRST   = LUMAT%FIRST
      LUMAT%FIRST         = swap
      swap8               = LUMAT_REMAP%NZL
      LUMAT_REMAP%NZL     = LUMAT%NZL
      LUMAT%NZL           = swap8
      swapl                     = LUMAT_REMAP%CONTIGUOUS_COL
      LUMAT_REMAP%CONTIGUOUS_COL= LUMAT%CONTIGUOUS_COL
      LUMAT%CONTIGUOUS_COL      = swapl
      PTR_COL              => LUMAT_REMAP%COL
      LUMAT_REMAP%COL      => LUMAT%COL
      LUMAT%COL            => PTR_COL
      CALL MUMPS_AB_FREE_LMAT(LUMAT_REMAP,K147)
 500  CONTINUE
      RETURN
      END SUBROUTINE MUMPS_REDIST_LUMAT_TOCONTLUMAT
      SUBROUTINE MUMPS_AB_DLUMAT_TO_BLK_DLUMAT (
     &        INFO, ICNTL, KEEP, COMM, MYID,
     &        NPROCS, NDOF, LUMAT, LP, LPOK, 
     &        DOF2BLOCK, NBLK, BLKFORMAT_DETECTED,
     &        MAPCOL, SIZE_MAPCOL,
     &        REMAP_CONTIG,  K244,
     &        NPROCS_PARAORD, PARAORD_to_idCOMM, 
     &        SIZE_PARAORD_to_idCOMM
     &        )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, MASTER
      PARAMETER( MASTER = 0 )
      INTEGER, INTENT(IN)   :: MYID, KEEP(500)
      INTEGER, INTENT(INOUT):: K244
      INTEGER, INTENT(IN)   :: ICNTL(60), COMM, NPROCS
      INTEGER               :: INFO(80)
      INTEGER, INTENT(IN)   :: NDOF, SIZE_MAPCOL
      INTEGER, INTENT(IN)            :: LP
      LOGICAL, INTENT(IN)            :: LPOK, REMAP_CONTIG
      INTEGER, INTENT(INOUT)         :: MAPCOL(SIZE_MAPCOL)
      INTEGER, INTENT(OUT)           :: NBLK
      TYPE(LMATRIX_T), INTENT(INOUT) :: LUMAT  
      INTEGER, INTENT(OUT)           :: DOF2BLOCK(NDOF)
      LOGICAL, INTENT(INOUT)         :: BLKFORMAT_DETECTED
      INTEGER, INTENT(IN)            :: NPROCS_PARAORD, 
     &                                  SIZE_PARAORD_to_idCOMM
      INTEGER, INTENT(IN) :: 
     &         PARAORD_to_idCOMM(SIZE_PARAORD_to_idCOMM)
      INTEGER :: NBLK_LOC, SHIFT, NBLK_I, SHIFT_I
      INTEGER :: allocok, I
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NBLK_ARRAY_IN, 
     &                                      NBLK_ARRAY_OUT, 
     &                                      TMP_DOF2BLOCK
      ALLOCATE(NBLK_ARRAY_IN(0:NPROCS-1), 
     &         NBLK_ARRAY_OUT(0:NPROCS-1), 
     &         TMP_DOF2BLOCK(NDOF),
     &         STAT=allocok)
      IF (allocok.NE.0) THEN
           INFO(1) = -7
           INFO(2) = 2*NPROCS + NDOF
           GOTO 500
      ENDIF
           CALL MUMPS_AB_LOCALCOMPUTE_BLKFORMAT ( 
     &        MYID, NDOF, LUMAT, KEEP(10), 
     &        INFO(1), INFO(2),  LP, LPOK,
     &        TMP_DOF2BLOCK, NBLK_LOC,
     &        MAPCOL, SIZE_MAPCOL
     &           )
      NBLK_ARRAY_IN       = 0
      NBLK_ARRAY_IN(MYID) = NBLK_LOC
      CALL MPI_ALLREDUCE(NBLK_ARRAY_IN(0), NBLK_ARRAY_OUT(0), NPROCS, 
     &        MPI_INTEGER, MPI_SUM, COMM, IERR)
      NBLK = 0
      DO I=0, NPROCS-1
       NBLK = NBLK + NBLK_ARRAY_OUT(I)
      ENDDO
      IF (NBLK.GT.
     &  int((real(KEEP(440))/real(100))*real(NDOF))
     &   ) THEN
         BLKFORMAT_DETECTED = .FALSE.
         NBLK                = NDOF
         DO I = 1, NDOF
           DOF2BLOCK(I)=I
         ENDDO
         DEALLOCATE(TMP_DOF2BLOCK)
      ELSE
        BLKFORMAT_DETECTED = .TRUE.
        SHIFT = 0
        IF (MYID.GT.0) THEN
         DO I=0, MYID-1 
           SHIFT = SHIFT + NBLK_ARRAY_OUT(I)
         ENDDO
        ENDIF
        DEALLOCATE(NBLK_ARRAY_IN, NBLK_ARRAY_OUT)
        DOF2BLOCK = TMP_DOF2BLOCK
        IF (SHIFT.GT.0) THEN
          DO I=LUMAT%FIRST,LUMAT%FIRST+LUMAT%NBCOL_LOC-1
               IF (TMP_DOF2BLOCK(I).EQ.0) CYCLE
               TMP_DOF2BLOCK(I) = TMP_DOF2BLOCK(I) + SHIFT
          ENDDO
        ENDIF
        CALL MPI_ALLREDUCE(TMP_DOF2BLOCK(1), DOF2BLOCK(1), NDOF,
     &      MPI_INTEGER, MPI_SUM, COMM, IERR)
         DEALLOCATE(TMP_DOF2BLOCK)
        IF (LUMAT%CONTIGUOUS_COL) THEN
         SHIFT_I = SHIFT+1
         NBLK_I  = NBLK_LOC
        ELSE
         SHIFT_I = 1
         NBLK_I  = NBLK
        ENDIF
        CALL MUMPS_AB_LOC_COMPRESSLUMAT (
     &        MYID, NDOF, LUMAT, KEEP(147), 
     &        INFO(1), INFO(2), LP, LPOK,
     &        DOF2BLOCK, NBLK, 
     &        SHIFT_I, NBLK_I
     &        )
      ENDIF
      IF (REMAP_CONTIG) THEN
        IF ((KEEP(244).NE.1) .AND. (NBLK.LT.2*NPROCS_PARAORD)) THEN
          K244= 1  
          IF (.NOT.LUMAT%CONTIGUOUS_COL) THEN
            CALL MUMPS_REDIST_LUMAT_TOCONTLUMAT (
     &       INFO, ICNTL, KEEP, COMM, MYID, NBLK,
     &       NPROCS, LUMAT, KEEP(147), 
     &       NPROCS, PARAORD_to_idCOMM, 
     &       SIZE_PARAORD_to_idCOMM
     &     )
          ENDIF
        ELSE
          CALL MUMPS_REDIST_LUMAT_TOCONTLUMAT (
     &       INFO, ICNTL, KEEP, COMM, MYID, NBLK,
     &       NPROCS, LUMAT, KEEP(147), 
     &       NPROCS_PARAORD, PARAORD_to_idCOMM, 
     &       SIZE_PARAORD_to_idCOMM
     &     )
         ENDIF
       ENDIF
  500 CONTINUE
      RETURN
      END SUBROUTINE MUMPS_AB_DLUMAT_TO_BLK_DLUMAT 
      SUBROUTINE MUMPS_AB_LOC_COMPRESSLUMAT (
     &        MYID, N, LUMAT, K147, 
     &        IFLAG, IERROR, LP, LPOK,
     &        DOF2BLOCK, NBLK, 
     &        FIRSTcomp, NBLK_LOCcomp
     &        )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T, COL_LMATRIX_T
      IMPLICIT NONE
      INTEGER, INTENT(IN)            :: MYID, N, K147, LP
      LOGICAL, INTENT(IN)            :: LPOK
      TYPE(LMATRIX_T), INTENT(INOUT) :: LUMAT
      INTEGER :: IFLAG, IERROR
      INTEGER, INTENT(IN)  :: DOF2BLOCK(N), NBLK 
      INTEGER, INTENT(IN)  :: FIRSTcomp, NBLK_LOCcomp
      TYPE(LMATRIX_T), TARGET                    :: LUMATcomp
      TYPE(COL_LMATRIX_T), DIMENSION(:), POINTER :: PTR_COL
      INTEGER, ALLOCATABLE, DIMENSION(:)         :: MAP
      INTEGER, POINTER, DIMENSION(:)             :: BLOCK_PTR
      INTEGER    :: allocok
      INTEGER    :: I, J, JJ, IB, JB, NB, swap, JDIAG
      INTEGER(8) :: NB_8, PT_CUR, SIZE_TO_ALLOC_8
      INTEGER*8  :: NZ8, swap8
      LOGICAL    :: swapl
      nullify(BLOCK_PTR)
      LUMATcomp%NBCOL          = NBLK
      LUMATcomp%NBCOL_LOC      = NBLK_LOCcomp
      LUMATcomp%FIRST          = FIRSTcomp
      LUMATcomp%CONTIGUOUS_COL = LUMAT%CONTIGUOUS_COL
      ALLOCATE(LUMATcomp%COL(NBLK_LOCcomp),
     &         MAP(NBLK), STAT=allocok)
      IF (allocok.NE.0) THEN
           IFLAG = -7
           CALL MUMPS_SET_IERROR(int(NBLK_LOCcomp,8)+int(N,8), IERROR)
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate of LUMATcomp%COL, MAP "
           END IF
           RETURN
      ENDIF
      DO IB=1,LUMATcomp%NBCOL_LOC
       LUMATcomp%COL(IB)%NBINCOL = 0
      ENDDO
      DO IB=1,NBLK
       MAP(IB)             = 0
      ENDDO
      DO J=1, LUMAT%NBCOL_LOC
        JDIAG =  J + LUMAT%FIRST - 1
        JB    = DOF2BLOCK(JDIAG)
        IF (  LUMATcomp%COL(JB-FIRSTcomp+1)%NBINCOL .GT.0 ) CYCLE
        NB=0
        DO I=1,  LUMAT%COL(J)%NBINCOL
         IB=DOF2BLOCK(LUMAT%COL(J)%IRN(I))
         IF (MAP(IB) .EQ. J) CYCLE
         MAP(IB) = J
         IF (IB.EQ.JB) CYCLE
         NB      = NB+1
        ENDDO
        LUMATcomp%COL(JB-FIRSTcomp+1)%NBINCOL = NB
      ENDDO
      DO JB=1,  LUMATcomp%NBCOL_LOC, K147
        SIZE_TO_ALLOC_8 = 0_8
        DO JJ=JB, min(JB+K147-1,LUMATcomp%NBCOL_LOC)
         SIZE_TO_ALLOC_8 = SIZE_TO_ALLOC_8 + 
     &                     int(LUMATcomp%COL(JJ)%NBINCOL, 8)
        ENDDO
        ALLOCATE(LUMATcomp%COL(JB)%IRN(SIZE_TO_ALLOC_8),STAT=allocok)
        IF (allocok.NE.0) THEN
           IFLAG  = -7
           CALL MUMPS_SET_IERROR(SIZE_TO_ALLOC_8, IERROR)
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate block of ", 
     &                K147, " columns in LUMATcomp%COL",
     &                " of size ", SIZE_TO_ALLOC_8
           END IF
           GOTO 500
        ENDIF
        IF (K147.GT.1) THEN
          BLOCK_PTR => LUMATcomp%COL(JB)%IRN
          NB_8    = int(LUMATcomp%COL(JB)%NBINCOL, 8) 
          PT_CUR = 1_8 + NB_8 
          DO JJ=JB+1, min(JB+K147-1,LUMATcomp%NBCOL_LOC)
            NB_8 =  int(LUMATcomp%COL(JJ)%NBINCOL,8)
            IF (NB_8.GT.0) THEN
             LUMATcomp%COL(JJ)%IRN => 
     &                     BLOCK_PTR(PT_CUR:PT_CUR+NB_8-1_8)
             PT_CUR = PT_CUR + NB_8
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      DO JB=1, LUMATcomp%NBCOL_LOC
        LUMATcomp%COL(JB)%NBINCOL = 0
      ENDDO
      DO IB=1,NBLK
       MAP(IB)             = 0
      ENDDO
      NZ8 = 0_8
      DO J=1, LUMAT%NBCOL_LOC
        JDIAG =  J + LUMAT%FIRST - 1
        JB = DOF2BLOCK(JDIAG)
        IF (  LUMATcomp%COL(JB-FIRSTcomp+1)%NBINCOL .GT.0 ) CYCLE
        NB = 0
        DO I=1,  LUMAT%COL(J)%NBINCOL
         IB=DOF2BLOCK(LUMAT%COL(J)%IRN(I))
         IF (MAP(IB) .EQ. J) CYCLE
         MAP(IB) = J
         IF (JB.EQ.IB) CYCLE
         NB      = NB+1
         LUMATcomp%COL(JB-FIRSTcomp+1)%IRN(NB) = IB
        ENDDO
        LUMATcomp%COL(JB-FIRSTcomp+1)%NBINCOL = NB
        NZ8  = NZ8 +int(NB,8)
      ENDDO
      LUMATcomp%NZL = NZ8
      swap                = LUMATcomp%NBCOL 
      LUMATcomp%NBCOL     = LUMAT%NBCOL
      LUMAT%NBCOL         = swap
      swap                = LUMATcomp%NBCOL_LOC 
      LUMATcomp%NBCOL_LOC = LUMAT%NBCOL_LOC
      LUMAT%NBCOL_LOC     = swap
      swap                = LUMATcomp%FIRST
      LUMATcomp%FIRST     = LUMAT%FIRST
      LUMAT%FIRST         = swap
      swap8               = LUMATcomp%NZL
      LUMATcomp%NZL       = LUMAT%NZL
      LUMAT%NZL           = swap8
      swapl                        = LUMATcomp%CONTIGUOUS_COL
      LUMATcomp%CONTIGUOUS_COL     = LUMAT%CONTIGUOUS_COL
      LUMAT%CONTIGUOUS_COL         = swapl
      PTR_COL              => LUMATcomp%COL
      LUMATcomp%COL        => LUMAT%COL
      LUMAT%COL            => PTR_COL
 500  CONTINUE
      IF (allocated(MAP)) DEALLOCATE(MAP)
      CALL MUMPS_AB_FREE_LMAT (LUMATcomp, K147)
      RETURN
      END SUBROUTINE MUMPS_AB_LOC_COMPRESSLUMAT
      SUBROUTINE MUMPS_AB_LOCALCOMPUTE_BLKFORMAT (
     &        MYID, N, LUMAT, KEEP10,  
     &        IFLAG, IERROR, LP, LPOK,
     &        DOF2BLOCK, NBLK, 
     &        MAPCOL, SIZE_MAPCOL
     &        )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      USE MUMPS_TOOLS_COMMON_M, ONLY : MUMPS_MERGESORTK8, 
     &                                 MUMPS_MERGESWAP48
      IMPLICIT NONE
      INTEGER, INTENT(IN)        :: MYID, N, LP, KEEP10
      LOGICAL, INTENT(IN)        :: LPOK 
      TYPE(LMATRIX_T), INTENT(IN):: LUMAT
      INTEGER, INTENT(IN)        :: SIZE_MAPCOL
      INTEGER, INTENT(IN)        :: MAPCOL(SIZE_MAPCOL)
      INTEGER :: IFLAG, IERROR
      INTEGER, INTENT(OUT) :: DOF2BLOCK(N), NBLK
      INTEGER, ALLOCATABLE, DIMENSION(:)    :: MARK, JDIAG
      INTEGER(8), ALLOCATABLE, DIMENSION(:) :: HASH
      INTEGER    :: allocok, I, II, J, JJ, JJJ, JD, NBCOL_LOC, 
     &              JSUIV, FIRST
      INTEGER(8) :: H
      LOGICAL    :: IDENTICAL
      NBCOL_LOC      = LUMAT%NBCOL_LOC
      FIRST          = LUMAT%FIRST
      NBLK      = 0
      DOF2BLOCK = 0
      IF (NBCOL_LOC.EQ.0) GOTO 600
      ALLOCATE(HASH(NBCOL_LOC), JDIAG(NBCOL_LOC), 
     &         MARK(0:max(N,NBCOL_LOC+1)), 
     &         STAT=allocok)
      IF (allocok.NE.0) THEN
           IFLAG = -7
           CALL MUMPS_SET_IERROR(
     &      int(max(N+1,NBCOL_LOC+2),8)*int(KEEP10,8)
     &      + int(NBCOL_LOC,8) +int(N,8), IERROR)
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate",
     &        " during MUMPS_AB_LOCALCOMPUTE_BLKFORMAT", IERROR
           END IF
           GOTO 600
      ENDIF
      FIRST    =  LUMAT%FIRST
      DO J=1, NBCOL_LOC
        JD      = J + FIRST - 1
        JDIAG(J) = JD
        IF (.NOT.LUMAT%CONTIGUOUS_COL) THEN
         IF (MAPCOL(J).NE.MYID) THEN
           HASH(J) = 0_8
           CYCLE
         ENDIF
        ENDIF
        H = int(JD,8)
        DO I=1,  LUMAT%COL(J)%NBINCOL
         II = LUMAT%COL(J)%IRN(I)
         H  = H + int(II,8)
        ENDDO
        HASH(J)  = H
      ENDDO
        CALL MUMPS_MERGESORTK8(NBCOL_LOC, HASH(1:NBCOL_LOC), 
     &         MARK(0:NBCOL_LOC+1) )
        CALL MUMPS_MERGESWAP48(NBCOL_LOC, MARK(0:NBCOL_LOC+1), 
     &          JDIAG(1:NBCOL_LOC), HASH(1:NBCOL_LOC))
      MARK      = 0
      DO J=1, NBCOL_LOC
        JJ       = JDIAG(J)   
        IF (.NOT.LUMAT%CONTIGUOUS_COL) THEN
           IF (MAPCOL(JJ).NE.MYID) CYCLE
        ENDIF
        MARK(JJ) = J
        DO  I=1,  LUMAT%COL(JJ-FIRST+1)%NBINCOL
          II = LUMAT%COL(JJ-FIRST+1)%IRN(I)
          MARK(II) = J
        ENDDO
        IF (DOF2BLOCK(JJ).EQ.0) THEN
          NBLK          = NBLK + 1
          DOF2BLOCK(JJ) = NBLK
          DO JSUIV=J+1, NBCOL_LOC
            JJJ  = JDIAG(JSUIV)
            IF (HASH(J).NE.HASH(JSUIV)) EXIT
            IF ( 
     &        (LUMAT%COL(JJ-FIRST+1)%NBINCOL .NE. 
     &         LUMAT%COL(JJJ-FIRST+1)%NBINCOL) 
     &         ) CYCLE
            IF (DOF2BLOCK(JJJ).EQ.0) THEN
              IDENTICAL = .TRUE.
              DO  I=1,  LUMAT%COL(JJJ-FIRST+1)%NBINCOL
                II = LUMAT%COL(JJJ-FIRST+1)%IRN(I)
                IF (MARK(II).NE.J) THEN
                    IDENTICAL = .FALSE.
                    EXIT
                ENDIF
              ENDDO
              IF (IDENTICAL) DOF2BLOCK(JJJ) = NBLK
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      IF (ALLOCATED(HASH))  DEALLOCATE(HASH)
      IF (ALLOCATED(JDIAG)) DEALLOCATE(JDIAG)
      IF (ALLOCATED(MARK))  DEALLOCATE(MARK)
  600 CONTINUE
      RETURN
      END SUBROUTINE MUMPS_AB_LOCALCOMPUTE_BLKFORMAT
      SUBROUTINE MUMPS_AB_FINALIZE_BLKFORMAT (
     &         N, DOF2BLOCK, NBLK, 
     &         BLKPTR, BLKVAR, SIZEOFBLOCKS
     &         )
       IMPLICIT NONE
      INTEGER, INTENT(IN)        :: N
      INTEGER, INTENT(IN) :: DOF2BLOCK(N), NBLK
      INTEGER, INTENT(OUT) :: BLKPTR(NBLK+1), BLKVAR(N), 
     &                        SIZEOFBLOCKS(NBLK)
      INTEGER I, IB
      DO IB=1,NBLK
        SIZEOFBLOCKS(IB)= 0
      ENDDO
      DO I=1,N
        IB = DOF2BLOCK(I)
        SIZEOFBLOCKS(IB) = SIZEOFBLOCKS(IB)+1
      ENDDO
      BLKPTR(1) = SIZEOFBLOCKS(1) + 1
      DO IB=2,NBLK
       BLKPTR(IB) = BLKPTR(IB-1)+SIZEOFBLOCKS(IB)
      ENDDO
      DO I=1,N
        IB                 = DOF2BLOCK(I)
        BLKPTR(IB)         = BLKPTR(IB)-1
        BLKVAR(BLKPTR(IB)) = I
      ENDDO
      BLKPTR(NBLK+1)= N+1
      RETURN
      END SUBROUTINE MUMPS_AB_FINALIZE_BLKFORMAT
#endif
