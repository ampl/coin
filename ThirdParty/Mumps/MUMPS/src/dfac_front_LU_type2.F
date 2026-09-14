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
      MODULE DMUMPS_FAC2_LU_M
      CONTAINS
      SUBROUTINE DMUMPS_FAC2_LU( COMM_LOAD, ASS_IRECV, 
     &           N, INODE, FPERE, IW, LIW, A, LA,
     &           UU, NOFFW, NPVW, NBTINYW,
     &           DET_EXPW, DET_MANTW, DET_SIGNW,
     &             COMM, MYID, BUFR, LBUFR,LBUFR_BYTES,NBFIN,LEAF,
     &             IFLAG, IERROR, IPOOL,LPOOL,
     &             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &             LRLUS, COMP,
     &             PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP,
     &             PIMASTER, PAMASTER,
     &             NSTK_S,PERM,PROCNODE_STEPS, root, roota,
     &             OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &             FILS, DAD, PTRARW, PTRAIW,
     &             PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &             LPTRAR, NELT, FRTPTR, FRTELT, SEUIL,
     &             ISTEP_TO_INIV2, TAB_POS_IN_PERE, AVOID_DELAYED,
     &             DKEEP,PIVNUL_LIST_STRUCT
     &               , LRGROUPS
     &             )
!$    USE OMP_LIB
      USE DMUMPS_FAC_FRONT_AUX_M
      USE DMUMPS_FAC_FRONT_TYPE2_AUX_M
      USE DMUMPS_OOC
      USE MUMPS_BUF_COMMON, ONLY : MUMPS_BUF_TEST
      USE DMUMPS_FAC_LR
      USE DMUMPS_LR_CORE
      USE DMUMPS_LR_TYPE
      USE MUMPS_LR_STATS
      USE DMUMPS_ANA_LR, ONLY : GET_CUT
      USE DMUMPS_LR_DATA_M
!$    USE OMP_LIB
      USE MUMPS_INTR_TYPES, ONLY : MUMPS_ROOT_STRUC
      USE DMUMPS_INTR_TYPES, ONLY : DMUMPS_ROOT_STRUC
      USE MUMPS_PIVNUL_MOD, ONLY : PIVNUL_LIST_STRUCT_T
      USE MUMPS_TOOLS_COMMON_M, ONLY : MUMPS_LR_COMPRESS_PANEL,
     &                                 MUMPS_LR_FORCE_BFR
      IMPLICIT NONE
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N, INODE, FPERE, LIW
      INTEGER, intent(inout) :: NOFFW, NPVW, NBTINYW
      INTEGER, intent(inout) :: DET_EXPW, DET_SIGNW
      DOUBLE PRECISION, intent(inout) :: DET_MANTW
      INTEGER(8) :: LA
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      DOUBLE PRECISION UU, SEUIL
      TYPE (MUMPS_ROOT_STRUC) :: root
      TYPE (DMUMPS_ROOT_STRUC) :: roota
      INTEGER COMM, MYID, LBUFR, LBUFR_BYTES
      INTEGER LPTRAR, NELT
      INTEGER ICNTL(60), KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER NBFIN, SLAVEF,
     &        IFLAG, IERROR, LEAF, LPOOL
      INTEGER(8) :: POSFAC, IPTRLU, LRLU, LRLUS
      INTEGER IWPOS, IWPOSCB, COMP 
      INTEGER FRTPTR( N + 1 ), FRTELT( NELT )
      INTEGER BUFR( LBUFR ), IPOOL(LPOOL),
     &        ITLOC(N+KEEP(253)), FILS(N), DAD( KEEP(28) ),
     &        ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER(8), INTENT(IN) :: PTRARW(LPTRAR), PTRAIW(LPTRAR)
      INTEGER(8), INTENT(IN) :: PTR8ARR(KEEP(193))
      INTEGER, INTENT(IN) :: NINCOLARR(KEEP(194))
      INTEGER, INTENT(IN) :: NINROWARR(KEEP(195))
      INTEGER, INTENT(IN) :: PTRDEBARR(KEEP(196))
      DOUBLE PRECISION :: RHS_MUMPS(KEEP8(85))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28)),
     &        STEP(N), PIMASTER(KEEP(28)),
     &        NSTK_S(KEEP(28)), PERM(N),
     &        PROCNODE_STEPS(KEEP(28))
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION DBLARR(KEEP8(26))
      INTEGER INTARR(KEEP8(27))
      LOGICAL AVOID_DELAYED
      TYPE(PIVNUL_LIST_STRUCT_T)     :: PIVNUL_LIST_STRUCT
      DOUBLE PRECISION DKEEP(230)
      INTEGER :: LRGROUPS(KEEP(280))
      INTEGER INOPV, IFINB, NFRONT, NPIV, IBEG_BLOCK, IEND_BLOCK
      INTEGER :: IBEG_BLOCK_FOR_IPIV
      INTEGER NASS, NBKJIB_ORIG, XSIZE
      INTEGER NBLR_ORIG, IBEG_BLR, IEND_BLR
      INTEGER Inextpiv
      LOGICAL LASTPANEL 
      INTEGER(8) :: POSELT
      INTEGER IOLDPS, allocok, K263, J
      INTEGER idummy 
      DOUBLE PRECISION    UUTEMP
      LOGICAL STATICMODE
      DOUBLE PRECISION SEUIL_LOC
      INTEGER , ALLOCATABLE, DIMENSION ( : ) :: IPIV
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, LNextPiv2beWritten, 
     &        UNextPiv2beWritten, IFLAG_OOC,
     &        PP_FIRST2SWAP_L, PP_FIRST2SWAP_U,
     &        PP_LastPIVRPTRFilled_L,
     &        PP_LastPIVRPTRFilled_U
      TYPE(IO_BLOCK) :: MonBloc 
      LOGICAL LAST_CALL
      INTEGER CURRENT_BLR, NELIM, NELIM_LOC
      LOGICAL LR_ACTIVATED, COMPRESS_PANEL, FORCE_BFR
      LOGICAL OOCWRITE_COMPATIBLE_WITH_BLR,
     &        OOC_EFFECTIVE_ON_FRONT, 
     &        OOC_EFF_AND_WRITE_BYPANEL
      INTEGER :: IROW_L, NVSCHUR, NSLAVES
      INTEGER :: PIVOT_OPTION, LAST_COL, FIRST_COL
      INTEGER :: PARPIV_T1
      INTEGER FIRST_BLOCK, LAST_BLOCK
      INTEGER :: INFO_TMP(2)
      INTEGER :: MAXI_RANK
      INTEGER HF, NPARTSASS, NPARTSCB, NB_BLR, END_I, NB_BLR_TMP
      INTEGER MAXI_CLUSTER, LWORK
      TYPE(LRB_TYPE), DIMENSION(1), TARGET  :: BLR_DUMMY
      INTEGER, POINTER, DIMENSION(:)        :: PTDummy
      TYPE(LRB_TYPE), POINTER, DIMENSION(:) :: ACC_LUA
      INTEGER, POINTER, DIMENSION(:)        :: BEGS_BLR
      TYPE(LRB_TYPE), POINTER, DIMENSION(:) :: BLR_L, BLR_U, BLR_SEND
      DOUBLE PRECISION, POINTER, DIMENSION(:)        :: DIAG
      TYPE(LRB_TYPE), POINTER, DIMENSION(:) :: BLR_PANEL
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_TMP, BEGS_BLR_STATIC
      INTEGER :: DIAGSIZ_STA, DIAGSIZ_DYN, DPOS, LorU, I, IP, MEM,
     &           MEM_TOT
      INTEGER(8) :: POSELT_DIAG
      CHARACTER(len=1) :: DIR
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:), TAU(:)
      INTEGER, ALLOCATABLE :: JPVT(:)
      DOUBLE PRECISION, ALLOCATABLE :: RWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: BLOCK(:,:)
      INTEGER :: OMP_NUM
      INTEGER(8) :: UPOS, LPOS
      INTEGER :: MY_NUM
#if defined(USE_GPU)
      LOGICAL :: run_on_gpu_slave, run_on_gpu_master
#endif
      LOGICAL :: POSTPROCESS_BLR, LR_ACTIV_BlocFACTO, LR_ACTIV_POSTPROC,
     &           MAXI_CLUSTER_UPDATED
      INTEGER :: END_I_IN_LOC
      LOGICAL :: FRSWAP_LOC
      INTEGER :: K473_LOC
      LOGICAL :: SWAP_OCCURRED
      INCLUDE 'mumps_headers.h'
      NULLIFY(BLR_L,BLR_U) 
      NULLIFY(PTDummy)
      NULLIFY(ACC_LUA)
      NULLIFY(BEGS_BLR)
      NULLIFY(BLR_L, BLR_U, BLR_SEND)
      NULLIFY(DIAG)
      NULLIFY(BLR_PANEL)
      NULLIFY( BEGS_BLR_TMP, BEGS_BLR_STATIC)
#if defined(USE_GPU)
      run_on_gpu_slave=.FALSE.
      run_on_gpu_master=.FALSE.
#endif
      IF (KEEP(206).GE.1) THEN
        Inextpiv = 1   
      ELSE 
        Inextpiv = 0   
      ENDIF
      SWAP_OCCURRED = .FALSE.
      idummy  = 0
      IOLDPS = PTLUST_S(STEP( INODE ))
      POSELT = PTRAST(STEP( INODE ))
      XSIZE  = KEEP(IXSZ)
      NFRONT = IW(IOLDPS+XSIZE)
      NASS   = abs(IW(IOLDPS+2+XSIZE))
      IW(IOLDPS+3+XSIZE) =  -99999
      LR_ACTIVATED   = (IW(IOLDPS+XXLR).GT.0)
      FORCE_BFR      = MUMPS_LR_FORCE_BFR(IW(IOLDPS+XXLR))
      COMPRESS_PANEL = MUMPS_LR_COMPRESS_PANEL(IW(IOLDPS+XXLR))
      COMPRESS_PANEL = FORCE_BFR.OR.COMPRESS_PANEL
      POSTPROCESS_BLR = 
     &     (KEEP(117).EQ.1.AND. IW(IOLDPS+XXG).NE.NotCandidateGPU)
      LR_ACTIV_BlocFACTO =  (LR_ACTIVATED.AND..NOT.POSTPROCESS_BLR)
      IF (
     &         (KEEP(486).EQ.2) 
     &       ) THEN
        LR_ACTIV_POSTPROC  = (LR_ACTIVATED .AND.POSTPROCESS_BLR)
      ELSE
        POSTPROCESS_BLR    = .FALSE.
        LR_ACTIV_POSTPROC  = .FALSE.
        LR_ACTIV_BlocFACTO =  LR_ACTIVATED
      ENDIF
      OOCWRITE_COMPATIBLE_WITH_BLR = 
     &          ( .NOT.LR_ACTIVATED.OR.(.NOT.COMPRESS_PANEL).OR.
     &            (KEEP(486).NE.2) 
     &          )
      OOC_EFFECTIVE_ON_FRONT= ((KEEP(201).EQ.1).AND. 
     &                         OOCWRITE_COMPATIBLE_WITH_BLR)
      PARPIV_T1 = 0
      INOPV = 0
      SEUIL_LOC = SEUIL
      IF(KEEP(97) .EQ. 0) THEN
         STATICMODE = .FALSE.
      ELSE
         STATICMODE = .TRUE.
      ENDIF
      IF (AVOID_DELAYED) THEN
         STATICMODE = .TRUE.
         UUTEMP=UU
         SEUIL_LOC = max(SEUIL,epsilon(SEUIL))
      ELSE
         UUTEMP=UU
      ENDIF
      FORCE_BFR = MUMPS_LR_FORCE_BFR( IW(IOLDPS+XXLR) )
      IF (FORCE_BFR) THEN
        K473_LOC = 1
      ELSE
        K473_LOC = KEEP(473)
      ENDIF
      CALL DMUMPS_SET_INNERBLOCKSIZE( NBKJIB_ORIG, NASS, KEEP)
      IF (.NOT. LR_ACTIV_BlocFACTO) THEN 
          NBLR_ORIG     = KEEP(420)
#if defined(USE_GPU)
          IF ( IW(IOLDPS+XXG).NE.NotCandidateGPU ) THEN
           IF (KEEP(417) .GE. 1) THEN
             NSLAVES   = IW(IOLDPS+5+XSIZE)
             CALL DMUMPS_RUNONGPU_BLOCKSIZEST2 (
     &         NFRONT, NASS, NSLAVES, KEEP, 
     &         run_on_gpu_slave, run_on_gpu_master,
     &         NBLR_ORIG, NBKJIB_ORIG)
           ENDIF
          ENDIF
#endif
      ELSE
          NBLR_ORIG  = -9999 
      ENDIF
      IF ((KEEP(114).EQ.1) .AND. 
     &    (KEEP(116).GT.0) .AND. ((NFRONT-NASS-KEEP(253)).GT.0) 
     &   ) THEN
         NSLAVES = IW(IOLDPS+5+XSIZE)
         IROW_L = IOLDPS+6+XSIZE+NSLAVES+NASS
         CALL DMUMPS_GET_SIZE_SCHUR_IN_FRONT ( 
     &     N, 
     &     NFRONT-NASS-KEEP(253), 
     &     KEEP(116), 
     &     IW(IROW_L), PERM, 
     &     NVSCHUR )
      ELSE
         NVSCHUR = 0
      ENDIF
      IF (LR_ACTIVATED) THEN
         K263 = 1   
      ELSE
         K263 = KEEP(263)
         IF (K263 .NE. 0 .AND. NASS/NBLR_ORIG < 4) THEN
           IF ( NBLR_ORIG .GT. NBKJIB_ORIG * 4 ) THEN
             NBLR_ORIG = max(NBKJIB_ORIG, (NASS+3)/4)
           ELSE
             K263 = 0
           ENDIF
         ENDIF
      ENDIF
#if defined(USE_GPU)
      IF (run_on_gpu_master) K263 = 1
#endif
      PIVOT_OPTION = KEEP(468)
      IF ( UUTEMP == 0.0D0 .AND.
     &    .NOT.(
     &      OOC_EFFECTIVE_ON_FRONT
     &         )
     &   ) THEN
          IF (K263.EQ.1.AND.(.NOT.LR_ACTIV_BlocFACTO)) THEN
            PIVOT_OPTION = 0
          ENDIF
      ENDIF
      IEND_BLOCK  = 0
      IEND_BLR    = 0
      IF (.NOT. LR_ACTIV_BlocFACTO) THEN
        CURRENT_BLR = -99999
      ELSE
        CURRENT_BLR = 0
      ENDIF
      ALLOCATE( IPIV( NASS ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        WRITE(*,*) MYID,' : DMUMPS_FAC2_LU :failed to allocate ',
     &  NASS, ' integers'
        IFLAG  = -13
        IERROR =NASS
        GO TO 500
      END IF
      CALL MUMPS_GETI8(LAFAC,IW(IOLDPS+XXR))
      LIWFAC    = IW(IOLDPS+XXI)
      IF ( OOC_EFFECTIVE_ON_FRONT ) THEN
          LNextPiv2beWritten = 1 
          UNextPiv2beWritten = 1 
          PP_FIRST2SWAP_L = LNextPiv2beWritten 
          PP_FIRST2SWAP_U = UNextPiv2beWritten 
          MonBloc%LastPanelWritten_L = 0
          MonBloc%LastPanelWritten_U = 0
          MonBloc%INODE    = INODE
          MonBloc%MASTER   = .TRUE.
          MonBloc%Typenode = 2
          MonBloc%NROW     = NASS
          MonBloc%NCOL     = NFRONT
          MonBloc%NFS      = NASS
          MonBloc%Last     = .FALSE.   
          MonBloc%LastPiv  = -68877    
          NULLIFY(MonBloc%INDICES)   
      ENDIF
      IF (LR_ACTIV_blocFACTO) THEN
             PIVOT_OPTION = 4
             IF (KEEP(475).EQ.1) THEN
                 PIVOT_OPTION = 3
             ELSEIF (KEEP(475).EQ.2) THEN
                 PIVOT_OPTION = 3
             ELSEIF (KEEP(475).EQ.3.AND.K263.EQ.1) THEN
               IF (UUTEMP == 0.0D0) THEN
                 PIVOT_OPTION = 0
               ELSE
                 PIVOT_OPTION = 1
               ENDIF
             ENDIF
      ENDIF
      IF (LR_ACTIVATED) THEN
             CNT_NODES = CNT_NODES + 1 
      ENDIF
      HF = 6 + IW(IOLDPS+5+XSIZE)+XSIZE
      OOC_EFF_AND_WRITE_BYPANEL  = ( (PIVOT_OPTION.GE.3) .AND.
     &                                     OOC_EFFECTIVE_ON_FRONT )
      IF (LR_ACTIVATED) THEN
         CALL GET_CUT(IW(IOLDPS+HF:IOLDPS+HF+NFRONT-1), NASS,
     &        NFRONT-NASS, LRGROUPS, NPARTSCB, 
     &        NPARTSASS, BEGS_BLR)
         CALL REGROUPING2(BEGS_BLR, NPARTSASS, NASS, NPARTSCB,
     &        NFRONT-NASS, KEEP(488), .FALSE., KEEP(472), 
     &        NFRONT, KEEP(1))     
         NB_BLR = NPARTSASS + NPARTSCB
         call MAX_CLUSTER(BEGS_BLR,NB_BLR,MAXI_CLUSTER)
         MAXI_RANK = KEEP(479)*MAXI_CLUSTER
         LWORK = MAXI_CLUSTER*(MAXI_CLUSTER+1)
         OMP_NUM = 1
#if ! defined(BLR_NOOPENMP)
!$       OMP_NUM = OMP_GET_MAX_THREADS()
#endif
         IF (LR_ACTIV_BlocFACTO) THEN
           ALLOCATE(BLOCK(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &             RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &             TAU(MAXI_CLUSTER*OMP_NUM),
     &             JPVT(MAXI_CLUSTER*OMP_NUM), 
     &             WORK(LWORK*OMP_NUM),stat=allocok)
           IF (allocok > 0) THEN
             IFLAG  = -13
             IERROR = OMP_NUM*(LWORK + MAXI_CLUSTER*(MAXI_CLUSTER+4))
             GOTO 480
           ENDIF
           ALLOCATE(ACC_LUA(OMP_NUM),stat=allocok)
           IF (allocok > 0) THEN
            IFLAG  = -13 
            IERROR = OMP_NUM
            GOTO 480
           ENDIF
           IF (KEEP(480).GE.3) THEN
             DO MY_NUM=1,OMP_NUM
               CALL ALLOC_LRB(ACC_LUA(MY_NUM), MAXI_RANK,
     &                      MAXI_CLUSTER, MAXI_CLUSTER, .TRUE., 
     &                      IFLAG, IERROR, KEEP8)
               IF (IFLAG.LT.0) GOTO 500
               ACC_LUA(MY_NUM)%K = 0
             ENDDO
           ENDIF
         ENDIF
      ENDIF
      IF (LR_ACTIVATED.AND.
     &    (KEEP(480).NE.0
     &       .OR.
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &    )
     &   ) THEN
        INFO_TMP(1) = IFLAG
        INFO_TMP(2) = IERROR
        CALL DMUMPS_BLR_INIT_FRONT(IW(IOLDPS+XXF), INFO_TMP)
        IFLAG  = INFO_TMP(1) 
        IERROR = INFO_TMP(2) 
        IF (IFLAG.LT.0) GOTO 500
        CALL DMUMPS_BLR_SAVE_INIT(IW(IOLDPS+XXF), 
     &              .FALSE., 
     &              .TRUE., 
     &              .FALSE., 
     &              NPARTSASS, 
     &              BEGS_BLR, PTDummy, 
     &              huge(NPARTSASS),  
     &              INFO_TMP)
        IFLAG  = INFO_TMP(1) 
        IERROR = INFO_TMP(2) 
        IF (IFLAG.LT.0) GOTO 500
      ENDIF
      LASTPANEL = .FALSE.
      DO WHILE (IEND_BLR < NASS ) 
        MAXI_CLUSTER_UPDATED = .FALSE.
        IBEG_BLR = IW(IOLDPS+1+KEEP(IXSZ)) + 1 
        IF (.NOT. LR_ACTIV_BlocFACTO) THEN
          IEND_BLR = min(IEND_BLR + NBLR_ORIG, NASS)
        ELSE
          CURRENT_BLR = CURRENT_BLR + 1
          IEND_BLR = BEGS_BLR(CURRENT_BLR+1)-1
          BEGS_BLR( CURRENT_BLR ) = IBEG_BLR
          IF ( IEND_BLR - IBEG_BLR + 1 .GT. MAXI_CLUSTER ) THEN
            MAXI_CLUSTER = IEND_BLR - IBEG_BLR + 1
            LWORK = MAXI_CLUSTER*(MAXI_CLUSTER+1)
            MAXI_CLUSTER_UPDATED = .TRUE.
           ENDIF
        ENDIF
        IF (LR_ACTIV_BlocFACTO) THEN
          IF ( MAXI_CLUSTER_UPDATED ) THEN
            DEALLOCATE(BLOCK, WORK, RWORK, TAU, JPVT)
            ALLOCATE(BLOCK(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &             RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &             TAU(MAXI_CLUSTER*OMP_NUM),
     &             JPVT(MAXI_CLUSTER*OMP_NUM), 
     &             WORK(LWORK*OMP_NUM),stat=allocok)
            IF (allocok > 0) THEN
              IFLAG  = -13
              IERROR = OMP_NUM*(LWORK + MAXI_CLUSTER*(MAXI_CLUSTER+4))
              GOTO 480
            ENDIF
            IF (KEEP(480).GE.3) THEN
              DO MY_NUM=1,OMP_NUM
                CALL DEALLOC_LRB(ACC_LUA(MY_NUM), KEEP8, KEEP(34))
                CALL ALLOC_LRB(ACC_LUA(MY_NUM), MAXI_RANK,
     &                         MAXI_CLUSTER, MAXI_CLUSTER, .TRUE.,
     &                         IFLAG, IERROR, KEEP8)
                IF (IFLAG.LT.0) GOTO 500
                ACC_LUA(MY_NUM)%K = 0
              ENDDO
            ENDIF
          ENDIF
        ENDIF
        DO WHILE (IEND_BLOCK < IEND_BLR ) 
          IBEG_BLOCK = IW(IOLDPS+1+KEEP(IXSZ)) + 1
          IF (KEEP(405).EQ.0) THEN
            KEEP(425)=max(KEEP(425),IEND_BLOCK-IBEG_BLOCK)
          ELSE
!$OMP       ATOMIC UPDATE
            KEEP(425)=max(KEEP(425),IEND_BLOCK-IBEG_BLOCK)
!$OMP       END ATOMIC
          ENDIF
          IEND_BLOCK = min(IEND_BLOCK + NBKJIB_ORIG, IEND_BLR)
  50      CONTINUE  
            IF (K263.EQ.0) THEN
              IBEG_BLOCK_FOR_IPIV = IBEG_BLOCK
            ELSE
              IBEG_BLOCK_FOR_IPIV = IBEG_BLR
            ENDIF
            CALL DMUMPS_FAC_I(NFRONT,NASS,NASS,
     &      IBEG_BLOCK_FOR_IPIV,IEND_BLOCK,N,INODE,
     &      IW,LIW,A,LA,INOPV,NOFFW,NBTINYW,
     &      DET_EXPW, DET_MANTW, DET_SIGNW,
     &      IFLAG,IERROR,IOLDPS,POSELT,UU,SEUIL_LOC,KEEP,KEEP8,
     &      DKEEP(1), PIVNUL_LIST_STRUCT, SWAP_OCCURRED,
     &      PP_FIRST2SWAP_L,  MonBloc%LastPanelWritten_L,
     &      PP_LastPIVRPTRFilled_L,
     &      PP_FIRST2SWAP_U,  MonBloc%LastPanelWritten_U,
     &      PP_LastPIVRPTRFilled_U,
     &      PIVOT_OPTION, LR_ACTIV_BlocFACTO, IEND_BLR,
     &      Inextpiv, OOC_EFFECTIVE_ON_FRONT,
     &      NVSCHUR, PARPIV_T1,
     &      TIPIV=IPIV 
     &      )
            IF (IFLAG.LT.0) GOTO 500  
          IF (INOPV.EQ.1) THEN
              IF (STATICMODE) THEN
                INOPV = -1
                GOTO 50
              ENDIF
              LASTPANEL = .TRUE.
          ELSE IF (INOPV .LE. 0) THEN 
            INOPV = 0
            IF (PIVOT_OPTION.GE.3) THEN
              LAST_COL = NFRONT
            ELSEIF (PIVOT_OPTION.EQ.2) THEN
              LAST_COL = NASS
            ELSE
              LAST_COL = IEND_BLR
            ENDIF
            CALL DMUMPS_FAC_MQ(IBEG_BLOCK, IEND_BLOCK,
     &              NFRONT, NASS, IW(IOLDPS+1+XSIZE),
     &              LAST_COL, A, LA, POSELT, IFINB,
     &              LR_ACTIV_BlocFACTO)
            IW(IOLDPS+1+XSIZE) = IW(IOLDPS+1+XSIZE) + 1
            NPVW = NPVW + 1
            IF (IFINB.EQ.0) THEN
              GOTO 50 
            ELSE IF (IFINB .EQ. -1) THEN
              LASTPANEL = .TRUE.
            ENDIF
          ENDIF
          NPIV = IW(IOLDPS+1+XSIZE)
          IF (K263.EQ.0) THEN
            NELIM = IEND_BLR - NPIV
            CALL DMUMPS_SEND_FACTORED_PANEL( COMM_LOAD, ASS_IRECV, 
     &           N, INODE, FPERE, IW, LIW, IOLDPS, POSELT, A, LA,
     &           NFRONT, IBEG_BLOCK, NPIV, IPIV,NASS,LASTPANEL,idummy, 
     &           COMM, MYID, BUFR, LBUFR, LBUFR_BYTES,NBFIN,LEAF,
     &           IFLAG, IERROR, IPOOL,LPOOL,
     &           SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &           LRLUS, COMP, PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP,
     &           PIMASTER, PAMASTER, NSTK_S,PERM,PROCNODE_STEPS,
     &           root, roota, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &           FILS, DAD, PTRARW, PTRAIW,
     &           PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &           INTARR,DBLARR, ICNTL,KEEP,KEEP8,
     &           DKEEP,ND,FRERE, LPTRAR, NELT, FRTPTR, FRTELT, 
     &           ISTEP_TO_INIV2, TAB_POS_IN_PERE, NELIM
     &           , .FALSE.  
     &           , LR_ACTIV_POSTPROC
     &           , NPARTSASS, CURRENT_BLR
     &           , BLR_DUMMY, LRGROUPS
     &           )
          END IF
          IF ( IFLAG .LT. 0 ) GOTO 500
          IF ( OOC_EFF_AND_WRITE_BYPANEL ) THEN
            MonBloc%LastPiv= IW(IOLDPS+1+XSIZE)
            STRAT          = STRAT_TRY_WRITE
            LAST_CALL      = .FALSE.
            CALL DMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT, TYPEF_U,
     &           A(POSELT), LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
            IF (IFLAG_OOC .LT. 0 ) THEN
              IFLAG = IFLAG_OOC
              GOTO 500
            ENDIF
          ENDIF
          NPIV       =  IW(IOLDPS+1+XSIZE)
          IF ( IEND_BLR .GT. IEND_BLOCK ) THEN
              CALL MUMPS_BUF_TEST()
              IF (PIVOT_OPTION.GE.3) THEN
                LAST_COL = NFRONT
              ELSEIF (PIVOT_OPTION.EQ.2) THEN
                LAST_COL = NASS
              ELSE 
                LAST_COL = IEND_BLR
              ENDIF
              CALL DMUMPS_FAC_SQ(IBEG_BLOCK, IEND_BLOCK,
     &            NPIV, NFRONT, IEND_BLR, LAST_COL,
     &            A, LA, POSELT, 
     &            -77777, 
     &            .TRUE., .FALSE., .TRUE.,
     &            .FALSE.,  
     &            LR_ACTIV_BlocFACTO)         
          ENDIF
          CALL MUMPS_BUF_TEST()
        END DO 
        NPIV   = IW(IOLDPS+1+XSIZE)
        IF (LR_ACTIV_BlocFACTO) THEN
          ALLOCATE(BLR_U(NB_BLR-CURRENT_BLR),stat=allocok)
          IF (allocok > 0) THEN
             IFLAG  = -13
             IERROR = NB_BLR-CURRENT_BLR
             GOTO 500
          ENDIF 
          ALLOCATE(BLR_L(NPARTSASS-CURRENT_BLR),stat=allocok)
          IF (allocok > 0) THEN
             IFLAG  = -13
             IERROR = NPARTSASS-CURRENT_BLR
             GOTO 500
          ENDIF 
          NELIM = IEND_BLR - NPIV
          IF (NELIM .EQ. IEND_BLR - IBEG_BLR + 1) THEN
            IF (KEEP(480).GE.2
     &       .OR.
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &          ) THEN
              DO J=1,NB_BLR-CURRENT_BLR
                 BLR_U(J)%M=0
                 BLR_U(J)%N=0
                 BLR_U(J)%K=0
                 BLR_U(J)%ISLR=.FALSE.
                 NULLIFY(BLR_U(J)%Q)
                 NULLIFY(BLR_U(J)%R)
              ENDDO
              CALL DMUMPS_BLR_SAVE_PANEL_LORU (
     &              IW(IOLDPS+XXF),
     &              1, 
     &        CURRENT_BLR, BLR_U, 0) 
              DO J=1,NPARTSASS-CURRENT_BLR
                 BLR_L(J)%M=0
                 BLR_L(J)%N=0
                 BLR_L(J)%K=0
                 BLR_L(J)%ISLR=.FALSE.
                 NULLIFY(BLR_L(J)%Q)
                 NULLIFY(BLR_L(J)%R)
              ENDDO
              CALL DMUMPS_BLR_SAVE_PANEL_LORU (
     &              IW(IOLDPS+XXF),
     &              0, 
     &        CURRENT_BLR, BLR_L, 0) 
              NULLIFY(BLR_L)
              NULLIFY(BLR_U)
            ENDIF
            GOTO 101
          ENDIF
            END_I=NB_BLR
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL PRIVATE(FIRST_BLOCK,LAST_BLOCK)
#endif
          CALL DMUMPS_COMPRESS_PANEL(A, LA, POSELT, IFLAG, IERROR, 
     &       NFRONT,
     &       BEGS_BLR, NB_BLR, DKEEP(8), KEEP(466),
     &       KEEP(458), K473_LOC, BLR_U, 
     &       CURRENT_BLR,
     &       'H', WORK, TAU, JPVT, LWORK, RWORK,
     &       BLOCK, MAXI_CLUSTER, NELIM, 
     &       .FALSE., 0, 0, 2, KEEP(483), KEEP8,
     &       (KEEP(450).GT.0), KEEP(453), KEEP(454).EQ.2,
     &       END_I_IN=END_I
     &        )
          IF (IFLAG.LT.0) GOTO 300
          IF ((KEEP(480).NE.0.AND.NB_BLR.GT.CURRENT_BLR)
     &       .OR.
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &       ) THEN
            IF (KEEP(480).LT.5) THEN
              CALL DMUMPS_BLR_SAVE_PANEL_LORU(
     &            IW(IOLDPS+XXF),
     &            1, 
     &        CURRENT_BLR, BLR_U, 0) 
            ENDIF
          ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif          
          IF (PIVOT_OPTION.LT.3) THEN
            IF (PIVOT_OPTION.LT.2) THEN
              FIRST_BLOCK = CURRENT_BLR+1
            ELSE
              FIRST_BLOCK = NPARTSASS+1
            ENDIF
              LAST_BLOCK=NB_BLR
            CALL DMUMPS_BLR_PANEL_LRTRSM(A, LA, POSELT, NFRONT,
     &           IBEG_BLR, 
     &           NB_BLR, BLR_U, CURRENT_BLR,
     &           FIRST_BLOCK, LAST_BLOCK, 2, 0, 1,
     &           .FALSE.)
          ENDIF
 300      CONTINUE         
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif          
        ENDIF
 101    CONTINUE       
        IF (LR_ACTIV_blocFACTO.OR. (K263.NE.0.AND.PIVOT_OPTION.GE.3)) 
     &    THEN
          NELIM = IEND_BLR - NPIV
          BLR_SEND=>BLR_DUMMY
          IF (associated(BLR_U)) THEN
            BLR_SEND=>BLR_U
          ENDIF
          CALL DMUMPS_SEND_FACTORED_PANEL( COMM_LOAD, ASS_IRECV, 
     &         N, INODE, FPERE, IW, LIW, IOLDPS, POSELT, A, LA, NFRONT,
     &         IBEG_BLR, NPIV, IPIV, NASS,LASTPANEL, idummy, 
     &         COMM, MYID, BUFR, LBUFR, LBUFR_BYTES,NBFIN,LEAF,
     &         IFLAG, IERROR, IPOOL,LPOOL,
     &         SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &         LRLUS, COMP, PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP,
     &         PIMASTER, PAMASTER, NSTK_S,PERM,PROCNODE_STEPS,
     &         root, roota, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &         FILS, DAD, PTRARW, PTRAIW,
     &         PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &         INTARR,DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND,FRERE,
     &         LPTRAR, NELT, FRTPTR, FRTELT, 
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE, NELIM
     &         , LR_ACTIV_BlocFACTO, LR_ACTIV_POSTPROC
     &         , NPARTSASS, CURRENT_BLR
     &         , BLR_SEND, LRGROUPS
     &         )
        ENDIF
        IF (.NOT.LR_ACTIV_BlocFACTO) THEN
          LAST_COL = NFRONT
          IF (PIVOT_OPTION.EQ.2) THEN
            FIRST_COL = NASS
          ELSE
            FIRST_COL = NPIV
          ENDIF
          IF (IEND_BLR.LT.NASS .OR. PIVOT_OPTION.LT.3) THEN
#if defined(USE_GPU)
            IF ( run_on_gpu_master ) THEN
              IF ( IW(IOLDPS+XXG) .EQ. PinningOnTheWay ) THEN
                CALL DMUMPS_GPU_WAIT_ENDOFPINNING_C()
                IW( IOLDPS + XXG ) = MemPinned
              ELSE IF ( KEEP(422) .EQ. 1 .AND. 
     &                  IW(IOLDPS+XXG) .EQ. MemNotPinned ) THEN
                CALL DMUMPS_GPU_PIN_MEMORY_C(A(POSELT),
     &          int(NASS,8)*int(NFRONT,8))
                IW( IOLDPS + XXG ) = MemPinned
              ENDIF
              CALL DMUMPS_FAC_SQ_GPU(INODE,
     &        KEEP, IBEG_BLR, IEND_BLR, NPIV, NFRONT, NASS,
     &        LAST_COL, A, LA, POSELT, FIRST_COL, .TRUE.,
     &        (PIVOT_OPTION.LT.3), .TRUE.)
            ELSE
#endif
            CALL DMUMPS_FAC_SQ(IBEG_BLR, IEND_BLR,
     &           NPIV, NFRONT, NASS, LAST_COL,
     &           A, LA, POSELT, FIRST_COL, .TRUE., (PIVOT_OPTION.LT.3),
     &           .TRUE., (KEEP(377).EQ.1), 
     &           LR_ACTIV_BlocFACTO)
#if defined(USE_GPU)
            ENDIF
#endif
          ENDIF
          IF (K263.NE.0 .AND. PIVOT_OPTION.LT.3) THEN
            NELIM = IEND_BLR - NPIV
            BLR_SEND=>BLR_DUMMY
            IF (associated(BLR_U)) THEN
              BLR_SEND=>BLR_U
            ENDIF
            CALL DMUMPS_SEND_FACTORED_PANEL( COMM_LOAD, ASS_IRECV, 
     &           N, INODE, FPERE, IW, LIW, IOLDPS, POSELT, A, LA,
     &           NFRONT, IBEG_BLR, NPIV, IPIV, NASS,LASTPANEL, idummy, 
     &           COMM, MYID, BUFR, LBUFR, LBUFR_BYTES,NBFIN,LEAF,
     &           IFLAG, IERROR, IPOOL,LPOOL,
     &           SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &           LRLUS, COMP, PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP,
     &           PIMASTER, PAMASTER, NSTK_S,PERM,PROCNODE_STEPS,
     &           root, roota, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &           FILS, DAD, PTRARW, PTRAIW,
     &           PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &           INTARR,DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND,FRERE,
     &           LPTRAR, NELT, FRTPTR, FRTELT, 
     &           ISTEP_TO_INIV2, TAB_POS_IN_PERE, NELIM
     &           , LR_ACTIV_BlocFACTO, LR_ACTIV_POSTPROC
     &           , NPARTSASS, CURRENT_BLR
     &           , BLR_SEND, LRGROUPS
     &           )
          ENDIF
        ELSE
         NELIM = IEND_BLR - NPIV
         IF (NELIM .EQ. IEND_BLR - IBEG_BLR + 1) THEN
            IF (KEEP(480).GE.2) THEN
              IF (IEND_BLR.LT.NASS) THEN
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL
#endif
                CALL DMUMPS_BLR_UPD_PANEL_LEFT(A, LA, POSELT,
     &          NFRONT, IW(IOLDPS+XXF), 0, 
     &          BEGS_BLR, BEGS_BLR, CURRENT_BLR, ACC_LUA,
     &          NB_BLR, NPARTSASS, NELIM,
     &          2, 0, 
     &          .FALSE., IFLAG, IERROR, 0,
     &          KEEP(481), DKEEP(11), KEEP(466), KEEP(477), 
     &          KEEP(480), KEEP(479), KEEP(478), KEEP(476), 
     &          KEEP(483), MAXI_CLUSTER, MAXI_RANK, 
     &          KEEP(474), 0, BLR_U, KEEP8
     &          )
                IF (IFLAG.LT.0) GOTO 600
                CALL DMUMPS_BLR_UPD_PANEL_LEFT(A, LA, POSELT,
     &          NFRONT, IW(IOLDPS+XXF), 1, 
     &          BEGS_BLR, BEGS_BLR, CURRENT_BLR, ACC_LUA,
     &          NB_BLR, NPARTSASS, NELIM,
     &          2, 0, 
     &          .FALSE., IFLAG, IERROR, 0,
     &          KEEP(481), DKEEP(11), KEEP(466), KEEP(477), 
     &          KEEP(480), KEEP(479), KEEP(478), KEEP(476), 
     &          KEEP(483), MAXI_CLUSTER, MAXI_RANK,
     &          KEEP(474), 0, BLR_U, KEEP8,
     &          END_I_IN=END_I
     &          )
 600            CONTINUE                
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif
                IF (IFLAG.LT.0) GOTO 500
              ENDIF
            ENDIF
            IF (KEEP(486).EQ.3) THEN
              IF (KEEP(480).EQ.0) THEN
                DEALLOCATE(BLR_U,BLR_L)
                NULLIFY(BLR_L)
                NULLIFY(BLR_U)
              ENDIF
            ENDIF
            GOTO 100
         ENDIF
         IF (KEEP(475).EQ.0) THEN
           IF (IEND_BLR.LT.NFRONT) THEN
              CALL DMUMPS_FAC_SQ(IBEG_BLR, IEND_BLR,
     &            NPIV, NFRONT, NASS, 
     &            -77777, 
     &            A, LA, POSELT,
     &            -77777, 
     &            .TRUE., .FALSE., .FALSE.,
     &            .FALSE.,  
     &            LR_ACTIVATED)
           ENDIF
          ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL PRIVATE(UPOS,LPOS,FIRST_BLOCK,LAST_BLOCK)
#endif
         CALL DMUMPS_COMPRESS_PANEL(A, LA, POSELT, IFLAG, IERROR, 
     &        NFRONT,
     &        BEGS_BLR, NPARTSASS, DKEEP(8), KEEP(466),
     &        KEEP(458), K473_LOC, 
     &        BLR_L,
     &        CURRENT_BLR, 'V', WORK, TAU, JPVT, LWORK, RWORK,
     &        BLOCK, MAXI_CLUSTER, NELIM,
     &        .FALSE., 0, 0,
     &        2, KEEP(483), KEEP8
     &        , (KEEP(450).GT.0), KEEP(453), KEEP(454).EQ.2
     &        )
#if ! defined(BLR_NOOPENMP)
!$OMP MASTER
#endif
         IF ((KEEP(480).NE.0.AND.NB_BLR.GT.CURRENT_BLR)
     &       .OR.
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &       ) THEN
            IF (KEEP(480).LT.5) THEN
              CALL DMUMPS_BLR_SAVE_PANEL_LORU (
     &             IW(IOLDPS+XXF),
     &             0, 
     &        CURRENT_BLR, BLR_L, 0) 
            ENDIF
          ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP END MASTER
!$OMP BARRIER
#endif          
          IF (IFLAG.LT.0) GOTO 400
          IF (KEEP(475).GT.0) THEN
            CALL DMUMPS_BLR_PANEL_LRTRSM(A, LA, POSELT, NFRONT,
     &                IBEG_BLR,
     &                NPARTSASS, BLR_L, CURRENT_BLR, CURRENT_BLR+1, 
     &                NPARTSASS, 2, 0, 0, .FALSE.) 
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif          
          ENDIF
          IF (KEEP(480).GE.2) THEN
            UPOS = POSELT+int(BEGS_BLR(CURRENT_BLR)-1,8)*int(NFRONT,8)
     &                   +int(BEGS_BLR(CURRENT_BLR+1)-NELIM-1,8)
            LPOS = POSELT+int(BEGS_BLR(CURRENT_BLR+1)-1,8)*int(NFRONT,8)
     &                   +int(BEGS_BLR(CURRENT_BLR+1)-NELIM-1,8)
            CALL DMUMPS_BLR_UPD_NELIM_VAR_L(A, LA, UPOS, A, LA, LPOS,
     &        IFLAG, IERROR, NFRONT, NFRONT,
     &        BEGS_BLR, CURRENT_BLR, BLR_L, NPARTSASS,
     &        CURRENT_BLR+1, NELIM, 'N') 
            IF (IFLAG.LT.0) GOTO 444
            IF (IEND_BLR.LT.NASS) THEN
              CALL DMUMPS_BLR_UPD_PANEL_LEFT(A, LA, POSELT,
     &          NFRONT, IW(IOLDPS+XXF), 0, 
     &          BEGS_BLR, BEGS_BLR, CURRENT_BLR, ACC_LUA,
     &          NB_BLR, NPARTSASS, NELIM,
     &          2, 0, 
     &          .FALSE., IFLAG, IERROR, 0,
     &          KEEP(481), DKEEP(11), KEEP(466), KEEP(477), 
     &          KEEP(480), KEEP(479), KEEP(478), KEEP(476), 
     &          KEEP(483), MAXI_CLUSTER, MAXI_RANK,
     &          KEEP(474), 0, BLR_U, KEEP8
     &          )
              IF (IFLAG.LT.0) GOTO 442
              CALL DMUMPS_BLR_UPD_PANEL_LEFT(A, LA, POSELT,
     &          NFRONT, IW(IOLDPS+XXF), 1, 
     &          BEGS_BLR, BEGS_BLR, CURRENT_BLR, ACC_LUA,
     &          NB_BLR, NPARTSASS, NELIM,
     &          2, 0, 
     &          .FALSE., IFLAG, IERROR, 0,
     &          KEEP(481), DKEEP(11), KEEP(466), KEEP(477), 
     &          KEEP(480), KEEP(479), KEEP(478), KEEP(476), 
     &          KEEP(483), MAXI_CLUSTER, MAXI_RANK,
     &          KEEP(474), 0, BLR_U, KEEP8,
     &          END_I_IN=END_I
     &          )
 442          CONTINUE
            ENDIF
 444        CONTINUE
          ELSE
            CALL DMUMPS_BLR_UPDATE_TRAILING(A, LA, POSELT, 
     &        IFLAG, IERROR, NFRONT,
     &        BEGS_BLR, BEGS_BLR, CURRENT_BLR, BLR_L, NPARTSASS,
     &        BLR_U, NB_BLR, NELIM, .FALSE., 0,
     &        2, 0, 
     &        KEEP(481), DKEEP(11), KEEP(466), KEEP(477) 
     &        )
          ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif
          IF (IFLAG.LT.0) GOTO 400
          IF (KEEP(475).GT.0) THEN
            FIRST_BLOCK = CURRENT_BLR+1
            IF (KEEP(486).EQ.2.AND.UU.EQ.0) THEN
              LAST_BLOCK = CURRENT_BLR
            ELSE
              LAST_BLOCK = NPARTSASS
            ENDIF
            CALL DMUMPS_DECOMPRESS_PANEL(A, LA, POSELT, NFRONT,
     &      NFRONT, .TRUE.,  
     &      BEGS_BLR(CURRENT_BLR),
     &      BEGS_BLR(CURRENT_BLR+1), NPARTSASS, BLR_L, CURRENT_BLR, 'V',
     &      1,
     &      BEG_I_IN=FIRST_BLOCK, END_I_IN=LAST_BLOCK)
#if ! defined(BLR_NOOPENMP)
#endif          
          ENDIF
         IF (KEEP(475).GE.2) THEN
           IF (KEEP(475).EQ.2) THEN
             FIRST_BLOCK = NPARTSASS+1
           ELSE
             FIRST_BLOCK = CURRENT_BLR+1
           ENDIF
           IF (KEEP(486).NE.2) THEN
             LAST_BLOCK = END_I
           ELSEIF(UU.GT.0) THEN
             LAST_BLOCK = NPARTSASS
           ELSE
             LAST_BLOCK = CURRENT_BLR
           ENDIF
           CALL DMUMPS_DECOMPRESS_PANEL(A, LA, POSELT, NFRONT,
     &       NFRONT, .TRUE.,   
     &       BEGS_BLR(CURRENT_BLR),
     &       BEGS_BLR(CURRENT_BLR+1), NB_BLR, BLR_U, CURRENT_BLR, 'H',
     &       1,
     &       BEG_I_IN=FIRST_BLOCK, END_I_IN=LAST_BLOCK)
         ENDIF
 400      CONTINUE
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif          
          IF (IFLAG.LT.0) GOTO 500
          IF (KEEP(486).EQ.3) THEN
            IF (KEEP(480).EQ.0.OR.NB_BLR.EQ.CURRENT_BLR) THEN
              CALL DEALLOC_BLR_PANEL(BLR_U, NB_BLR-CURRENT_BLR,
     &                      KEEP8, KEEP(34))
              CALL DEALLOC_BLR_PANEL(BLR_L, NPARTSASS-CURRENT_BLR, 
     &                      KEEP8, KEEP(34))
              DEALLOCATE(BLR_U,BLR_L)
            ENDIF
          ENDIF
          NULLIFY(BLR_L)
          NULLIFY(BLR_U)
        ENDIF
        IF ( OOC_EFFECTIVE_ON_FRONT ) THEN
          STRAT            = STRAT_TRY_WRITE
          MonBloc%LastPiv  = NPIV
          LAST_CALL= .FALSE.
          CALL DMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT, TYPEF_BOTH_LU, 
     &           A(POSELT), LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
          IF (IFLAG_OOC .LT. 0 ) THEN
            IFLAG = IFLAG_OOC
            GOTO 500
          ENDIF
        ENDIF
 100    CONTINUE
      END DO 
      IF (LR_ACTIVATED) THEN
        NPIV     = IW(IOLDPS+1+XSIZE)
        IBEG_BLR = NPIV + 1
        IF (LR_ACTIV_POSTPROC) THEN
          BEGS_BLR( NPARTSASS  + 1 )  = IBEG_BLR
          NELIM   = NASS - IW(IOLDPS+1+XSIZE)
        ELSE
          BEGS_BLR( CURRENT_BLR + 1 ) = IBEG_BLR
        ENDIF
        IF ( (KEEP(486).EQ.2) 
     &       ) THEN
          CALL DMUMPS_BLR_RETRIEVE_BEGSBLR_STA(IW(IOLDPS+XXF),
     &                      BEGS_BLR_STATIC)
          NPIV = IW(IOLDPS+1+XSIZE)
          IF (LR_ACTIV_POSTPROC) THEN
             DO IP=1,NPARTSASS
               IF ( BEGS_BLR_STATIC(IP).GE. NPIV ) THEN
                 BEGS_BLR(IP)        = NPIV
               ENDIF
             ENDDO
             MAXI_CLUSTER = max (MAXI_CLUSTER,  
     &                           NASS - BEGS_BLR(NPARTSASS) + 1)
          ENDIF
          IF ( (UU.GT.0 .AND. SWAP_OCCURRED)
     &                       .OR.(LR_ACTIV_POSTPROC)) THEN
            allocate(BEGS_BLR_TMP(NB_BLR+1),stat=allocok)
            IF (allocok > 0) THEN
              IFLAG  = -13
              IERROR = NB_BLR+1
              GOTO 500
            ENDIF
            DO IP=1,NB_BLR+1
               BEGS_BLR_TMP(IP) = BEGS_BLR_STATIC(IP)
            ENDDO
          ENDIF
        ENDIF
        IF (LR_ACTIV_POSTPROC) THEN
          LWORK = MAXI_CLUSTER*(MAXI_CLUSTER+1)
          ALLOCATE(BLOCK(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &              RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &              TAU(MAXI_CLUSTER*OMP_NUM),
     &              JPVT(MAXI_CLUSTER*OMP_NUM), 
     &              WORK(LWORK*OMP_NUM),
     &              stat=allocok)
          IF (allocok > 0) THEN
            IFLAG  = -13
            IERROR = OMP_NUM*(LWORK + MAXI_CLUSTER*(MAXI_CLUSTER+4))
            GOTO 480
          ENDIF
          DO IP=1,NPARTSASS
              DO LorU=0,1
                IF (LorU.EQ.0) THEN
                  NB_BLR_TMP = NPARTSASS
                ELSE
                  NB_BLR_TMP = NB_BLR
                ENDIF
                ALLOCATE(BLR_PANEL(NB_BLR_TMP-IP),stat=allocok)
                IF (allocok > 0) THEN
                  IFLAG  = -13
                  IERROR = NB_BLR_TMP-IP
                  GOTO 480
                ENDIF
                IF (BEGS_BLR(IP+1)-BEGS_BLR(IP) .EQ.0) THEN
                 DO J=1, NB_BLR_TMP-IP
                   BLR_PANEL(J)%M=0
                   BLR_PANEL(J)%N=0
                   BLR_PANEL(J)%K=0
                   BLR_PANEL(J)%ISLR=.FALSE.
                   NULLIFY(BLR_PANEL(J)%Q)
                   NULLIFY(BLR_PANEL(J)%R)
                 ENDDO
                ENDIF
                CALL DMUMPS_BLR_SAVE_PANEL_LORU (
     &              IW(IOLDPS+XXF),
     &              LorU, 
     &              IP, BLR_PANEL, 0) 
              ENDDO
          ENDDO
        ENDIF
        IF (
     &         (KEEP(486).EQ.2) 
     &       ) THEN
          MEM_TOT = 0
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL
!$OMP& PRIVATE(IP, LorU, DIR, NELIM_LOC, BLR_PANEL, 
!$OMP&         END_I_IN_LOC, FRSWAP_LOC)
#endif
#if ! defined(BLR_NOOPENMP)
!$OMP DO PRIVATE(DIAG, DIAGSIZ_STA, DIAGSIZ_DYN, DPOS, POSELT_DIAG, MEM,
!$OMP&           allocok)
!$OMP&   REDUCTION(+:MEM_TOT)
#endif
          DO IP=1,NPARTSASS
            IF (IFLAG.LT.0) CYCLE
            DIAGSIZ_DYN = BEGS_BLR(IP+1)-BEGS_BLR(IP)
            DIAGSIZ_STA = BEGS_BLR_STATIC(IP+1)-BEGS_BLR(IP)
            MEM = DIAGSIZ_DYN*(2*DIAGSIZ_STA-DIAGSIZ_DYN)
            IF (MEM.EQ.0) CYCLE 
            MEM_TOT = MEM_TOT + MEM
            ALLOCATE(DIAG(MEM), stat=allocok)
            IF (allocok > 0) THEN
              IFLAG  = -13
              IERROR = MEM
              CYCLE
            ENDIF 
            DPOS = 1
            POSELT_DIAG = POSELT + int(BEGS_BLR(IP)-1,8)*int(NFRONT,8)
     &                             + int(BEGS_BLR(IP)-1,8)
            DO I=1,DIAGSIZ_STA
              IF (I.LE.DIAGSIZ_DYN) THEN
                DIAG(DPOS:DPOS+DIAGSIZ_STA-1) =
     &                 A(POSELT_DIAG:POSELT_DIAG+int(DIAGSIZ_STA-1,8))
                DPOS = DPOS + DIAGSIZ_STA
              ELSE
                DIAG(DPOS:DPOS+DIAGSIZ_DYN-1) =
     &                 A(POSELT_DIAG:POSELT_DIAG+int(DIAGSIZ_DYN-1,8))
                DPOS = DPOS + DIAGSIZ_DYN
              ENDIF
              POSELT_DIAG = POSELT_DIAG + int(NFRONT,8)
            ENDDO
            CALL DMUMPS_BLR_SAVE_DIAG_BLOCK(
     &            IW(IOLDPS+XXF),
     &            IP, DIAG, KEEP(34))
          ENDDO
#if ! defined(BLR_NOOPENMP)
!$OMP ENDDO
!$OMP SINGLE
#endif
          ENTRIES_PER_PREC_DIAG(PREC_DIAG_FR) =
     &       ENTRIES_PER_PREC_DIAG(PREC_DIAG_FR) + MEM_TOT
          MEM_TOT = 1 + 
     &      (MEM_TOT*COSTS_PREC(PREC_DIAG_FR)-1)/COSTS_PREC(PREC_FAC)
          CALL MUMPS_DM_FAC_UPD_DYN_MEMCNTS(int(MEM_TOT,8),
     &          .FALSE., KEEP8, IFLAG, IERROR, .TRUE., .TRUE.)
#if ! defined(BLR_NOOPENMP)
!$OMP END SINGLE
#endif
          IF (IFLAG.LT.0) GOTO 460
          IF ((UU.GT.0.AND.SWAP_OCCURRED).OR.(LR_ACTIV_POSTPROC)) THEN
            DO IP=1,NPARTSASS
              NELIM_LOC = BEGS_BLR_TMP(IP+1)-BEGS_BLR(IP+1)
              DO LorU=0,1
                CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &            IW(IOLDPS+XXF), LorU, IP, BLR_PANEL)
#if ! defined(BLR_NOOPENMP)
!$OMP SINGLE
#endif
                IF (LR_ACTIV_BlocFACTO) THEN
                    CALL DEALLOC_BLR_PANEL(BLR_PANEL, NPARTSASS-IP,
     &                                 KEEP8, KEEP(34))
                ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP END SINGLE
#endif
                IF (LorU.EQ.0) THEN 
                  DIR = 'V'
                ELSE
                  DIR = 'H'
                ENDIF
                IF (LR_ACTIV_POSTPROC) THEN
                  IF (LorU.EQ.0) THEN
                    END_I_IN_LOC = NPARTSASS
                  ELSE
                    END_I_IN_LOC = NB_BLR
                  ENDIF
                  FRSWAP_LOC   = .FALSE.
                ELSE
                  END_I_IN_LOC = NPARTSASS
                  FRSWAP_LOC   = .TRUE.
                ENDIF
                CALL DMUMPS_COMPRESS_PANEL(A, LA, POSELT, IFLAG,
     &            IERROR, NFRONT, BEGS_BLR_TMP,
     &            NB_BLR, DKEEP(8), KEEP(466), KEEP(458), K473_LOC,
     &            BLR_PANEL, IP,
     &            DIR, WORK, TAU, JPVT, LWORK, RWORK,
     &            BLOCK, MAXI_CLUSTER, NELIM_LOC,
     &            .FALSE., 0, 0,
     &            2, KEEP(483), KEEP8,
     &          (KEEP(450).GT.0), KEEP(453), KEEP(454).GE.1,
     &          END_I_IN=END_I_IN_LOC, FRSWAP=FRSWAP_LOC
     &          )
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif
                IF (IFLAG.LT.0) GOTO 440
              ENDDO
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
!$OMP SINGLE
#endif
              BEGS_BLR_TMP(IP+1) = BEGS_BLR(IP+1)
#if ! defined(BLR_NOOPENMP)
!$OMP END SINGLE
#endif
            ENDDO
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif
 440        CONTINUE
          ENDIF 
 460      CONTINUE          
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif
          IF ( (UU.GT.0 .AND. SWAP_OCCURRED)
     &                       .OR.(LR_ACTIV_POSTPROC)) THEN
            deallocate(BEGS_BLR_TMP)
          ENDIF
        ENDIF  
        IF (IFLAG.LT.0) GOTO 500
        IF (
     &         (KEEP(486).EQ.2) 
     &       ) THEN
            CALL DMUMPS_BLR_SAVE_BEGS_BLR_DYN(IW(IOLDPS+XXF),
     &        BEGS_BLR)
        ENDIF
        IF (
     &        ((KEEP(486).EQ.2)
     &       ) .AND. (KEEP(450).GT.0)) THEN
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL PRIVATE(IP, BLR_PANEL)
#endif
          DO IP=1,NPARTSASS
                IF ( BEGS_BLR(IP).EQ. BEGS_BLR(IP+1)) CYCLE
                CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &             IW(IOLDPS+XXF), 0, IP, BLR_PANEL)
                CALL DMUMPS_LRMIXED_CONVERT_PANEL(BLR_PANEL, NPARTSASS,
     &                                          IP, KEEP(34), KEEP8, 
     &                                          IFLAG, IERROR)
                CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &             IW(IOLDPS+XXF), 1, IP, BLR_PANEL)
                CALL DMUMPS_LRMIXED_CONVERT_PANEL(BLR_PANEL, NB_BLR,
     &                                          IP, KEEP(34), KEEP8,
     &                                          IFLAG, IERROR)
          ENDDO
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif
          IF (IFLAG.LT.0) GOTO 480
        ENDIF
        NELIM = NASS-IW(IOLDPS+1+XSIZE)
        IF (IFLAG.GE.0) THEN
            CALL UPD_MRY_LU_FR(NASS, NFRONT-NASS, 0, NELIM)
            DO IP=1,NPARTSASS
              IF ( BEGS_BLR(IP).EQ. BEGS_BLR(IP+1)) CYCLE
              CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &               IW(IOLDPS+XXF), 0, IP, BLR_PANEL)
              CALL UPD_MRY_LU_LRGAIN(BLR_PANEL, NPARTSASS-IP
     &                , KEEP(450) .GT. 0
     &               )
              CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &               IW(IOLDPS+XXF), 1, IP, BLR_PANEL)
              CALL UPD_MRY_LU_LRGAIN(BLR_PANEL, NB_BLR-IP
     &                , KEEP(450) .GT. 0
     &                )
            ENDDO
            CALL UPD_FLOP_FACTO_FR(NFRONT, NASS, NASS-NELIM, 0, 2)
        ENDIF
      ENDIF
#if defined(USE_GPU)
        IF ( run_on_gpu_master 
     &       .AND. IW(IOLDPS+XXG) .EQ. MemPinned ) THEN
          IF ( KEEP(422) .EQ. 2 ) THEN
            CALL DMUMPS_GPU_ASYNC_UNPIN_MEMORY_C(A(POSELT),
     &           int(NASS,8)*int(NFRONT,8))
            IW(IOLDPS+XXG) = MemNotPinned
          ELSE  
            CALL DMUMPS_GPU_UNPIN_MEMORY_C(A(POSELT),
     &           int(NASS,8)*int(NFRONT,8))
            IW(IOLDPS+XXG) = MemNotPinned
          ENDIF
        ENDIF
#endif
      IF (KEEP(486).NE.0) THEN
        IF (.NOT.LR_ACTIV_BlocFACTO) THEN
          CALL UPD_FLOP_FRFRONTS(NFRONT, NPIV, NASS, 0, 2)
        ENDIF
      ENDIF
      IF (IFLAG.LT.0) GOTO 500
      IF ( OOC_EFFECTIVE_ON_FRONT ) THEN
          STRAT        = STRAT_WRITE_MAX   
          MonBloc%Last = .TRUE.
          MonBloc%LastPiv  = IW(IOLDPS+1+XSIZE)
          LAST_CALL    = .TRUE.
          CALL DMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT, TYPEF_BOTH_LU,
     &           A(POSELT), LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC, LAST_CALL )
          IF (IFLAG_OOC .LT. 0 ) THEN
            IFLAG = IFLAG_OOC
            GOTO 500
          ENDIF
          CALL DMUMPS_OOC_PP_TRYRELEASE_SPACE (IWPOS, 
     &      IOLDPS, IW, LIW, MonBloc , NFRONT, KEEP)
      ENDIF
      GOTO 700
 480  CONTINUE
 500  CONTINUE
      CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
 700  CONTINUE
      IF (LR_ACTIVATED) THEN
        IF (allocated(RWORK))  DEALLOCATE(RWORK)
        IF (allocated(WORK)) DEALLOCATE(WORK)
        IF (allocated(TAU)) DEALLOCATE(TAU)
        IF (allocated(JPVT)) DEALLOCATE(JPVT)
        IF (allocated(BLOCK)) DEALLOCATE(BLOCK)
        IF (associated(ACC_LUA)) THEN
         IF (KEEP(480).GE.3) THEN
           DO MY_NUM=1,OMP_NUM
             CALL DEALLOC_LRB(ACC_LUA(MY_NUM), KEEP8, KEEP(34))
           ENDDO
         ENDIF
         DEALLOCATE(ACC_LUA)
        ENDIF
        IF (associated(BEGS_BLR)) THEN
           DEALLOCATE(BEGS_BLR)
           NULLIFY(BEGS_BLR)
        ENDIF
      ENDIF
      IF (LR_ACTIVATED.AND.KEEP(480).NE.0) THEN
        IF (.NOT.
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &     ) THEN
          CALL DMUMPS_BLR_FREE_ALL_PANELS(IW(IOLDPS+XXF), 2, 
     &                    KEEP8, KEEP(34))
        ENDIF
      ENDIF
      IF ( LR_ACTIVATED .AND. KEEP(486).EQ. 2 .AND.
     &     KEEP(251) .EQ. 2) THEN
          CALL DMUMPS_BLR_FREE_ALL_PANELS(IW(IOLDPS+XXF), 0, 
     &                       KEEP8, KEEP(34))
      ENDIF
      IF (LR_ACTIVATED) THEN
        IF (.NOT.
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &     ) THEN
          CALL DMUMPS_BLR_END_FRONT(IW(IOLDPS+XXF), IFLAG, KEEP8,
     &                              KEEP(34))
        ENDIF
      ENDIF
      DEALLOCATE( IPIV )
      RETURN
      END SUBROUTINE DMUMPS_FAC2_LU
      END MODULE DMUMPS_FAC2_LU_M
