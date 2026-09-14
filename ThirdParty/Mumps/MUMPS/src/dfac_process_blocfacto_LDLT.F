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
      RECURSIVE SUBROUTINE DMUMPS_PROCESS_SYM_BLOCFACTO( 
     &   COMM_LOAD, ASS_IRECV,
     &   BUFR, LBUFR,
     &   LBUFR_BYTES, PROCNODE_STEPS, MSGSOU,
     &   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     &   A, LA, PTRIST, PTRAST, NSTK_S, PERM,
     &   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     &   MYID, COMM, IFLAG, IERROR, NBFIN,
     &
     &    PTLUST_S, PTRFAC, root, roota, OPASSW, OPELIW,
     &    ITLOC, RHS_MUMPS, FILS, DAD,
     &    PTRARW, PTRAIW,
     &    PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &    INTARR, DBLARR,
     &    ICNTL,KEEP,KEEP8,DKEEP,IPOOL, LPOOL, LEAF, ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &               , LRGROUPS
     &    )
      USE DMUMPS_OOC, ONLY : IO_BLOCK
      USE MUMPS_OOC_COMMON, ONLY : TYPEF_L, STRAT_TRY_WRITE
      USE MUMPS_TOOLS_COMMON_M, ONLY : MUMPS_LR_COMPRESS_CB,
     &                                 MUMPS_LR_COMPRESS_PANEL,
     &                                 MUMPS_LR_FORCE_BFR
      USE MUMPS_LOAD
      USE MUMPS_BUF_COMMON, ONLY : MUMPS_BUF_SEND_1INT
      USE DMUMPS_BUF, ONLY : DMUMPS_BUF_SEND_BLFAC_SLAVE
      USE DMUMPS_LR_CORE
      USE DMUMPS_LR_TYPE
      USE MUMPS_LR_STATS
      USE DMUMPS_FAC_LR
      USE DMUMPS_ANA_LR, ONLY : GET_CUT
      USE DMUMPS_LR_DATA_M
      USE MUMPS_INTR_TYPES, ONLY : MUMPS_ROOT_STRUC
      USE DMUMPS_INTR_TYPES, ONLY : DMUMPS_ROOT_STRUC
      USE DMUMPS_DYNAMIC_MEMORY_M, ONLY : DMUMPS_DM_SET_DYNPTR,
     &    DMUMPS_DM_ALLOC_S_WK, DMUMPS_DM_FREE_S_WK
      USE DMUMPS_FAC_FRONT_AUX_M, ONLY : DMUMPS_GET_SIZE_SCHUR_IN_FRONT
!$    USE OMP_LIB
      IMPLICIT NONE
      INCLUDE 'mumps_headers.h'
      TYPE (MUMPS_ROOT_STRUC) :: root
      TYPE (DMUMPS_ROOT_STRUC) :: roota
      INTEGER ICNTL( 60 ), KEEP( 500 )
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION    DKEEP(230)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER N, SLAVEF, IWPOS, IWPOSCB, LIW
      INTEGER(8) IPTRLU, LRLU, LRLUS, LA, POSFAC
      INTEGER COMP
      INTEGER IFLAG, IERROR, NBFIN, MSGSOU
      INTEGER PROCNODE_STEPS(KEEP(28)), PTRIST(KEEP(28)),
     &        NSTK_S(KEEP(28))
      INTEGER(8) PTRAST(KEEP(28)), PTRFAC(KEEP(28)), PAMASTER(KEEP(28))
      INTEGER PERM(N), STEP(N), 
     & PIMASTER(KEEP(28))
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER, intent(in) :: LRGROUPS(KEEP(280))
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER COMM, MYID
      INTEGER PTLUST_S(KEEP(28)),
     &        ITLOC(N+KEEP(253)), FILS(N), DAD(KEEP(28)), ND(KEEP(28))
      DOUBLE PRECISION :: RHS_MUMPS(KEEP8(85))
      INTEGER(8), INTENT(IN) :: PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER(8), INTENT(IN) :: PTR8ARR(KEEP(193))
      INTEGER, INTENT(IN) :: NINCOLARR(KEEP(194))
      INTEGER, INTENT(IN) :: NINROWARR(KEEP(195))
      INTEGER, INTENT(IN) :: PTRDEBARR(KEEP(196))
      INTEGER FRERE_STEPS(KEEP(28))
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION FLOP1
      INTEGER INTARR( KEEP8(27) )
      DOUBLE PRECISION DBLARR( KEEP8(26) )
      INTEGER LEAF, LPOOL 
      INTEGER IPOOL( LPOOL )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER PIVI
      INTEGER (8) POSPV1,POSPV2,OFFDAG,LPOS1
      INTEGER J2
      DOUBLE PRECISION MULT1,MULT2, A11, DETPIV, A22, A12
      INTEGER :: NFS4FATHER, NVSCHUR_K253, NSLAVES_L, IROW_L
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: M_ARRAY
      INTEGER NBROWSinF
      INTEGER :: BLFAC_NBCOLS_ALREADY_SENT, BLFAC_NBLRB_ALREADY_SENT
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER LP
      INTEGER INODE, POSITION, NPIV, IERR
      INTEGER NEWCOL_RECV, JBEG_BLOCK, NCOL_GEMM_FR, 
     &        SHIFT_LPOS, SHIFT_UPOS
      INTEGER :: IFLAG_OOC
      INTEGER(8) :: POSBLOCFACTO
      INTEGER :: LD_BLOCFACTO 
      INTEGER(8) :: LA_BLOCFACTO 
      INTEGER(8) :: LA_PTR 
      INTEGER(8) :: POSELT
      DOUBLE PRECISION, DIMENSION(:), POINTER :: A_PTR
      INTEGER IOLDPS, NASS1, NROW1, NCOL1, NPIV1
      INTEGER NSLAV1, HS, ISW, DEST
      INTEGER ICT11
      INTEGER(8) LPOS, LPOS2, DPOS, UPOS
      INTEGER (8) IPOS, KPOS
      INTEGER I, IPIV, FPERE, NSLAVES_TOT,
     &        NSLAVES_FOLLOW, NB_BLOC_FAC
      INTEGER IPOSK, JPOSK, NPIVSENT, Block, IROW, BLSIZE
      INTEGER allocok, TO_UPDATE_CPT_END
      INTEGER(8) :: LUIP21K
      DOUBLE PRECISION, DIMENSION(:), POINTER :: UIP21K
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DYN_BLOCFACTO
      INTEGER, DIMENSION(:), ALLOCATABLE :: LIST_SLAVES_FOLLOW
      INTEGER, DIMENSION(:), ALLOCATABLE :: DYN_PIVINFO
      LOGICAL LASTPANEL
      LOGICAL LASTBL_INPANEL, LASTBL_INLASTPANEL
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      DOUBLE PRECISION ONE,ALPHA
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      INTEGER LIWFAC, STRAT, NextPivDummy
      TYPE(IO_BLOCK) :: MonBloc
      LOGICAL LAST_CALL
      INTEGER J
      LOGICAL COUNTER_WAS_HUGE
      INTEGER TO_UPDATE_CPT_RECUR
      INTEGER :: LR_ACTIVATED_INT
      LOGICAL :: LR_ACTIV_BlocFACTO, LR_ACTIVATED, COMPRESS_CB, 
     &           COMPRESS_PANEL, FORCE_BFR
      LOGICAL :: LR_ACTIV_POSTPROC
      LOGICAL :: DYNAMIC_ALLOC
      LOGICAL OOCWRITE_COMPATIBLE_WITH_BLR
      INTEGER :: CURRENT_BLR, INFO_TMP(2), ISHIFT
      INTEGER :: NELIM, NB_BLR_LM, NB_BLR_LS,  
     &           MAXI_CLUSTER_LM, MAXI_CLUSTER_LS, MAXI_CLUSTER, 
     &           NPARTSASS, NPARTSCB, NPARTSASS_COL, 
     &           NB_BLR_COL, MAXI_CLUSTER_COL
       INTEGER :: NPARTSASS_MASTER, IPANEL, IP, NB_ACCESSES_INIT, 
     &            NB_ACCESSES_LEFT_INIT
      TYPE (LRB_TYPE), DIMENSION(:), ALLOCATABLE :: BLR_LM 
      TYPE (LRB_TYPE), DIMENSION(:), POINTER     :: BLR_LS
      TYPE(LRB_TYPE), POINTER, DIMENSION(:,:) :: CB_LRB
      INTEGER                                 :: SIZE_BEGS_BLR_COL
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_LM, BEGS_BLR_LS, 
     &                                  BEGS_BLR_COL
      LOGICAL KEEP_BEGS_BLR_LS, KEEP_BEGS_BLR_COL, KEEP_BLR_LS
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: WORK, TAU
      INTEGER, ALLOCATABLE, DIMENSION(:) :: JPVT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: BLOCKLR
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: RWORK
      INTEGER :: OMP_NUM, LWORK
      INTEGER :: II,JJ
#if defined(USE_GPU)
      INTEGER, EXTERNAL :: DMUMPS_TRSM_ON_GPU_C, DMUMPS_GEMM_ON_GPU_C,
     &                     DMUMPS_GEMMT_ON_GPU_C
      INTEGER :: run_trsm_on_gpu, run_gemm_on_gpu, run_gemmt_on_gpu,
     &           async_gemm
      INTEGER(8) :: TMP_SIZE_PIN, TMP_SIZE_UNPIN
#endif
      LOGICAL :: NOTHING_WAS_SENT
      INTEGER :: KEEP430_LOC
      INTEGER :: NB, IB, IBEG, IEND
!$    INTEGER :: NOMP
!$    LOGICAL :: OMP_FLAG
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      OMP_NUM = 1
#if ! defined(BLR_NOOPENMP)
!$    OMP_NUM = OMP_GET_MAX_THREADS()
#endif
      KEEP(174)=KEEP(174)+1
      KEEP(175)=max(KEEP(174),KEEP(175))
      LP = ICNTL(1)
      IF (ICNTL(4) .LE. 0) LP = -1
      POSITION = 0
      TO_UPDATE_CPT_END = -654321
      NULLIFY(UIP21K)
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NPIV, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, FPERE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      LASTPANEL = (NPIV.LE.0)
      IF (LASTPANEL) THEN 
         NPIV = -NPIV
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NB_BLOC_FAC, 1,
     &                 MPI_INTEGER, COMM, IERR )
      ENDIF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NEWCOL_RECV, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, JBEG_BLOCK, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NELIM, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, 
     &                 NPARTSASS_MASTER, 1,
     &                 MPI_INTEGER, COMM, IERR )
       NPARTSASS_COL = NPARTSASS_MASTER
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, IPANEL,
     &                 1, MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, LR_ACTIVATED_INT, 1,
     &                 MPI_INTEGER, COMM, IERR )
      LR_ACTIV_BlocFACTO    = (LR_ACTIVATED_INT.EQ.1)
      IF (JBEG_BLOCK.EQ.1) THEN
        NCOL_GEMM_FR  = NEWCOL_RECV - NPIV
        SHIFT_LPOS = NPIV
        SHIFT_UPOS = NPIV
      ELSE
        SHIFT_LPOS = JBEG_BLOCK - 1
        IF (LR_ACTIV_BlocFACTO) THEN
          NCOL_GEMM_FR  = -99993
          SHIFT_UPOS = -99994
        ELSE
          NCOL_GEMM_FR  = NEWCOL_RECV
          SHIFT_UPOS = 0
        ENDIF
      ENDIF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NSLAVES_TOT, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, 
     &               SIZE_BEGS_BLR_COL, 1,
     &               MPI_INTEGER, COMM, IERR )
      KEEP_BEGS_BLR_COL =.FALSE.
      NULLIFY(BEGS_BLR_COL)
      IF (SIZE_BEGS_BLR_COL.GT.0) THEN
          ALLOCATE(BEGS_BLR_COL(SIZE_BEGS_BLR_COL),
     &        stat=allocok)
          IF (allocok > 0 ) THEN
             IFLAG = -13
             IERROR = SIZE_BEGS_BLR_COL
             IF (LP > 0) WRITE(LP,*) MYID,
     &       ': ERROR allocation during DMUMPS_PROCESS_BLOCFACTO_LDLT'
             GOTO 700
          ENDIF
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, 
     &                 BEGS_BLR_COL, SIZE_BEGS_BLR_COL,
     &                 MPI_INTEGER, COMM, IERR )
      ENDIF
      KEEP_BEGS_BLR_LS  =.FALSE. 
      NULLIFY(BEGS_BLR_LS)
      KEEP_BLR_LS       =.FALSE.
      NULLIFY(BLR_LS)
      NULLIFY(BEGS_BLR_LM)
      IF ( LR_ACTIV_BlocFACTO ) THEN
        LA_BLOCFACTO = int(NPIV,8) * int(NPIV+NELIM,8)
        LD_BLOCFACTO = max(NPIV+NELIM,1)
      ELSE
        LA_BLOCFACTO = int(NPIV,8) * int(NEWCOL_RECV,8)
        LD_BLOCFACTO = max(NEWCOL_RECV,1)
      ENDIF
      IF (LR_ACTIV_BlocFACTO) THEN
        DYNAMIC_ALLOC = .TRUE.
      ELSE
        DYNAMIC_ALLOC = .FALSE.
      ENDIF
      IF ( .NOT. DYNAMIC_ALLOC ) THEN
        IF ( NPIV .EQ. 0 ) THEN
          IPIV = 1  
          POSBLOCFACTO = 1_8 
        ELSE 
          CALL DMUMPS_GET_SIZE_NEEDED(
     &       NPIV, LA_BLOCFACTO, .FALSE.,
     &       KEEP(1), KEEP8(1), 
     &       N, IW, LIW, A, LA,
     &       LRLU, IPTRLU,
     &       IWPOS, IWPOSCB, PTRIST, PTRAST,
     &       STEP, PIMASTER, PAMASTER, LRLUS,
     &       KEEP(IXSZ),COMP,DKEEP(97),
     &       MYID, SLAVEF, PROCNODE_STEPS, DAD,
     &       IFLAG, IERROR)
          IF (IFLAG.LT.0) GOTO 700
          LRLU  = LRLU - LA_BLOCFACTO
          LRLUS = LRLUS - LA_BLOCFACTO
          KEEP8(69) = KEEP8(69) + LA_BLOCFACTO
          KEEP8(67) = min(LRLUS, KEEP8(67))
          KEEP8(68) = max(KEEP8(69), KEEP8(68))
          POSBLOCFACTO = POSFAC
          POSFAC = POSFAC + LA_BLOCFACTO
          IPIV = IWPOS
          IWPOS = IWPOS + NPIV
          CALL MUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &               LA-LRLUS,0_8,LA_BLOCFACTO,KEEP,KEEP8,LRLUS)
        ENDIF
      ELSE 
        ALLOCATE(DYN_PIVINFO(max(1,NPIV)),
     &           DYN_BLOCFACTO(max(1_8,LA_BLOCFACTO)),
     &           stat=allocok)
        IF (allocok.GT.0) THEN
           IF (LP > 0 ) WRITE(LP,*) MYID,
     &     ": ALLOCATION FAILURE FOR DYN_PIVINFO and DYN_BLOCFACTO IN ",
     &     "DMUMPS_PROCESS_SYM_BLOCFACTO"
           IFLAG = -13
           CALL MUMPS_SET_IERROR(max(1_8,LA_BLOCFACTO), IERROR)
           GOTO 700
        ENDIF
        KEEP8(130)=KEEP8(130)+max(1_8,LA_BLOCFACTO)
        KEEP8(131)=max(KEEP8(130),KEEP8(131))
        KEEP8(73) = KEEP8(73) + max(1_8,LA_BLOCFACTO)
        KEEP8(69) = KEEP8(69) + max(1_8,LA_BLOCFACTO)
        KEEP8(74) = max(KEEP8(74), KEEP8(73))
        KEEP8(68) = max(KEEP8(68), KEEP8(69))
        POSBLOCFACTO = 1_8
        IPIV = 1
      ENDIF
      IF (NPIV.GT.0) THEN
         IF (DYNAMIC_ALLOC) THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 DYN_PIVINFO, NPIV,
     &                 MPI_INTEGER, COMM, IERR )
         ELSE
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 IW( IPIV ), NPIV,
     &                 MPI_INTEGER, COMM, IERR )
         ENDIF
        IF (DYNAMIC_ALLOC) THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                DYN_BLOCFACTO, int(LA_BLOCFACTO),
     &                 MPI_DOUBLE_PRECISION,
     &                 COMM, IERR )
        ELSE
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 A(POSBLOCFACTO), int(LA_BLOCFACTO),
     &                 MPI_DOUBLE_PRECISION,
     &                 COMM, IERR )
        ENDIF
        IF ( LR_ACTIV_BlocFACTO ) THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 NB_BLR_LM, 1, MPI_INTEGER,
     &                 COMM, IERR )
          ALLOCATE(BLR_LM(max(NB_BLR_LM,1)), stat=allocok)
          IF ( allocok .GT. 0 ) THEN
             IF (LP > 0 ) WRITE(LP,*) MYID,
     &            ": ALLOCATION FAILURE FOR BLR_LM IN ",
     &            "DMUMPS_PROCESS_SYM_BLOCFACTO"
             IFLAG = -13
             IERROR = max(NB_BLR_LM,1)
             GOTO 700
          END IF
          ALLOCATE(BEGS_BLR_LM(NB_BLR_LM+2), stat=allocok)
          IF ( allocok .GT. 0 ) THEN
             IF (LP > 0 ) WRITE(LP,*) MYID,
     &            ": ALLOCATION FAILURE FOR BEGS_BLR_LM IN ",
     &            "DMUMPS_PROCESS_SYM_BLOCFACTO"
             IFLAG = -13
             IERROR = NB_BLR_LM+2
             GOTO 700
          END IF
          CALL DMUMPS_MPI_UNPACK_LR_PARTIAL(
     &          BUFR, LBUFR, LBUFR_BYTES, POSITION, NPIV, NELIM, 
     &          'V', BLR_LM, NB_BLR_LM, JBEG_BLOCK,
     &          BEGS_BLR_LM(1), KEEP8, COMM, IERR, IFLAG, IERROR)
          IF (IFLAG.LT.0) GOTO 700
        ENDIF 
      ENDIF 
      IF (PTRIST(STEP( INODE )) .EQ. 0) THEN
          CALL DMUMPS_TREAT_DESCBAND( INODE, COMM_LOAD, ASS_IRECV,
     &      BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &      IWPOS, IWPOSCB, IPTRLU,
     &      LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &      PTLUST_S, PTRFAC,
     &      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &      IFLAG, IERROR, COMM,
     &      PERM, IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &      root, roota, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &      FILS, DAD, PTRARW, PTRAIW,
     &      PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &      INTARR, DBLARR, ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &      LPTRAR, NELT, FRTPTR, FRTELT, 
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.
     &               , LRGROUPS
     &        )
          IF ( IFLAG .LT. 0 ) GOTO 550
      ENDIF
      IF ( IW( PTRIST(STEP(INODE)) + 3 + KEEP(IXSZ)) .EQ. 0 ) THEN
       DO WHILE ( IW(PTRIST(STEP(INODE)) + XXNBPR) .NE. 0)
        BLOCKING  = .TRUE.
        SET_IRECV = .FALSE.
        MESSAGE_RECEIVED = .FALSE.
        CALL DMUMPS_TRY_RECVTREAT( COMM_LOAD,
     &    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    MPI_ANY_SOURCE, CONTRIB_TYPE2,
     &    STATUS,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM,
     &    PERM, IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &    root, roota, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, DAD, PTRARW, PTRAIW,
     &    PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &    INTARR, DBLARR, ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &               , LRGROUPS
     &      )
        IF ( IFLAG .LT. 0 ) GOTO 550
      END  DO
      ENDIF
        SET_IRECV = .TRUE.
        BLOCKING  = .FALSE.
        MESSAGE_RECEIVED = .TRUE.
        CALL DMUMPS_TRY_RECVTREAT( COMM_LOAD, ASS_IRECV,
     &    BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    MPI_ANY_SOURCE, MPI_ANY_TAG, 
     &    STATUS,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM,
     &    PERM, IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &    root, roota, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, DAD, PTRARW, PTRAIW,
     &    PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &    INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &               , LRGROUPS
     &    )
      IOLDPS = PTRIST(STEP(INODE))
      CALL DMUMPS_DM_SET_DYNPTR( IW(IOLDPS+XXS), A, LA,
     &     PTRAST(STEP(INODE)), IW(IOLDPS+XXD), IW(IOLDPS+XXR),
     &     A_PTR, POSELT, LA_PTR )
      NASS1  = IW( IOLDPS + 1 + KEEP(IXSZ))
      LR_ACTIVATED   = (IW(IOLDPS+XXLR).GT.0)
      FORCE_BFR = MUMPS_LR_FORCE_BFR(IW(IOLDPS+XXLR))
      COMPRESS_PANEL = MUMPS_LR_COMPRESS_PANEL(IW(IOLDPS+XXLR))
      COMPRESS_PANEL = FORCE_BFR.OR.COMPRESS_PANEL
      IF (
     &         (KEEP(486).EQ.2) 
     &       ) THEN
        LR_ACTIV_POSTPROC =
     &         COMPRESS_PANEL.AND.(.NOT.LR_ACTIV_BlocFACTO)
      ELSE
        LR_ACTIV_POSTPROC = .FALSE.
      ENDIF
      OOCWRITE_COMPATIBLE_WITH_BLR = 
     &          ( (.NOT.COMPRESS_PANEL).OR.
     &            (KEEP(486).NE.2) 
     &          )
      IF ( NASS1 < 0 ) THEN
        NASS1 = -NASS1
        IW( IOLDPS + 1 + KEEP(IXSZ)) = NASS1
        IF (KEEP(55) .EQ. 0) THEN 
          CALL DMUMPS_ASM_SLAVE_ARROWHEADS(INODE, STEP(INODE),
     &       N, IW, LIW,
     &       IOLDPS, A_PTR(POSELT), LA_PTR, 1_8, KEEP, KEEP8, ITLOC,
     &       FILS,
     &       PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &       INTARR, DBLARR, KEEP8(27), KEEP8(26), RHS_MUMPS, 
     &       LRGROUPS)
        ELSE
          CALL DMUMPS_ASM_SLAVE_ELEMENTS(INODE, N, NELT, IW, LIW,
     &       IOLDPS, A_PTR(POSELT), LA_PTR, 1_8, KEEP, KEEP8, ITLOC,
     &       FILS, PTRAIW,
     &       PTRARW, INTARR, DBLARR, KEEP8(27), KEEP8(26),
     &       FRTPTR, FRTELT, RHS_MUMPS, LRGROUPS)
        ENDIF
      ENDIF
      NROW1              = IW( IOLDPS + 2 +KEEP(IXSZ))
      NCOL1  = IW( IOLDPS + 3 +KEEP(IXSZ)) + IW( IOLDPS + KEEP(IXSZ))
      IF (JBEG_BLOCK.EQ.1) THEN
        NPIV1  = IW( IOLDPS + 3 +KEEP(IXSZ))
      ELSE
        NPIV1  = IW( IOLDPS + 3 +KEEP(IXSZ)) - NPIV
      ENDIF
      LASTBL_INPANEL     = JBEG_BLOCK+NEWCOL_RECV.GT.NASS1-NPIV1
      LASTBL_INLASTPANEL = LASTPANEL .AND. LASTBL_INPANEL
      NSLAV1 = IW( IOLDPS + 5 + KEEP(IXSZ))
      NSLAVES_FOLLOW = NSLAV1 - XTRA_SLAVES_SYM
      IF (SIZE_BEGS_BLR_COL.GT.0) THEN
         DO I=1, SIZE_BEGS_BLR_COL
           IF (BEGS_BLR_COL(I).GT.NCOL1+1) THEN
             BEGS_BLR_COL(I)=NCOL1+1
           ENDIF
         ENDDO
      ENDIF
      HS     = 6 + NSLAV1 + KEEP(IXSZ)
      IF ( LASTBL_INLASTPANEL ) THEN
        TO_UPDATE_CPT_END = ( NSLAVES_TOT - NSLAVES_FOLLOW - 1 ) * 
     &                       NB_BLOC_FAC
      END IF
#if defined(USE_GPU)
      run_trsm_on_gpu =0
      run_gemm_on_gpu =0
      run_gemmt_on_gpu=0
#endif
      IF (NPIV.GT.0) THEN
#if defined(USE_GPU)
        IF ( .NOT. LR_ACTIV_BlocFACTO .OR. KEEP(475) .EQ. 0) THEN
          run_trsm_on_gpu=DMUMPS_TRSM_ON_GPU_C(
     &                  KEEP(417), NPIV, NROW1, 'L')
        ENDIF
        IF ( .NOT. LR_ACTIV_BlocFACTO ) THEN
          IF (NCOL_GEMM_FR.GT.0) THEN
            run_gemm_on_gpu=DMUMPS_GEMM_ON_GPU_C(
     &                  KEEP(417), NCOL_GEMM_FR, NROW1, NPIV )
          ENDIF
          IF (JBEG_BLOCK.EQ.1) THEN
            run_gemmt_on_gpu=DMUMPS_GEMMT_ON_GPU_C(
     &                  KEEP(417), NROW1, NPIV )
          ENDIF
        ENDIF
        IF ( KEEP(422).EQ.2 ) THEN
          IF (IW(IOLDPS+XXG) .EQ. PinningOnTheWay) THEN
            CALL DMUMPS_GPU_WAIT_ENDOFPINNING_C()
            IW(IOLDPS+XXG) = MemPinned
          ENDIF
        ELSE IF (KEEP(422) .EQ.1 ) THEN
          IF ( IW(IOLDPS+XXG) .NE. MemPinned ) THEN
            IF ( run_trsm_on_gpu.EQ.1 .OR. run_gemm_on_gpu.EQ.1
     &         .OR. run_gemmt_on_gpu.EQ.1 ) THEN
              CALL DMUMPS_GPU_PIN_MEMORY_C(A_PTR(POSELT), LA_PTR)
              IW(IOLDPS+XXG) = MemPinned
              IW(IOLDPS+XXG) = MemPinned
            ENDIF
          ENDIF
        ENDIF
#endif
        IF (JBEG_BLOCK.EQ.1) THEN
          ICT11 = IOLDPS+HS+NROW1+NPIV1 - 1
          DO I = 1, NPIV
            IF (DYNAMIC_ALLOC) THEN
              PIVI = abs(DYN_PIVINFO(I))
            ELSE
              PIVI = abs(IW(IPIV+I-1))
            ENDIF
            IF (PIVI.EQ.I) CYCLE
            ISW = IW(ICT11+I)
            IW(ICT11+I) = IW(ICT11+PIVI)
            IW(ICT11+PIVI) = ISW
            IPOS = POSELT + int(NPIV1 + I - 1,8)
            KPOS = POSELT + int(NPIV1 + PIVI - 1,8)
            CALL dswap(NROW1, A_PTR(IPOS), NCOL1, A_PTR(KPOS), NCOL1)
          ENDDO
        IF (LR_ACTIV_BlocFACTO) THEN
          LUIP21K = 1_8
        ELSE
          LUIP21K=int(NPIV,8)*int(NROW1,8)
        ENDIF
        KEEP430_LOC=min(KEEP(430),1)
#if defined(USE_GPU)
#if defined(USE_XKBLAS)
        IF (run_gemmt_on_gpu.EQ.1) THEN
          IF (KEEP(430).EQ.2) THEN
            KEEP430_LOC=2
          ENDIF
        ENDIF
#endif
#endif
        CALL DMUMPS_DM_ALLOC_S_WK( UIP21K, LUIP21K, allocok,
     &                             KEEP430_LOC, KEEP(35) )
        IF ( allocok .GT. 0 ) THEN
          IF (LP > 0 ) WRITE(LP,*) MYID,
     &": ALLOCATION FAILURE FOR UIP21K IN DMUMPS_PROCESS_SYM_BLOCFACTO"
          IFLAG = -13
          IERROR = NPIV * NROW1
          GOTO 700
        END IF
        KEEP8(130)=KEEP8(130)+LUIP21K
        KEEP8(131)=max(KEEP8(130),KEEP8(131))
        KEEP8(73) = KEEP8(73) + LUIP21K
        KEEP8(69) = KEEP8(69) + LUIP21K
        KEEP8(74) = max(KEEP8(74), KEEP8(73))
        KEEP8(68) = max(KEEP8(68), KEEP8(69))
        IF (.NOT.LR_ACTIV_BlocFACTO) THEN
#         if defined(USE_GPU)
          IF ( KEEP(422).NE.0 .AND. KEEP(423) .EQ.1 .AND.
     &        ( run_gemmt_on_gpu.EQ.1 .OR. run_gemm_on_gpu.EQ.1 ) ) THEN
            TMP_SIZE_PIN = 0_8
            IF (run_gemm_on_gpu .EQ. 1) THEN
              TMP_SIZE_PIN=TMP_SIZE_PIN + LA_BLOCFACTO
              IF (DYNAMIC_ALLOC) THEN
                CALL DMUMPS_GPU_PIN_MEMORY_C(
     &               DYN_BLOCFACTO, LA_BLOCFACTO)
              ELSE
                CALL DMUMPS_GPU_PIN_MEMORY_C(
     &               A(POSBLOCFACTO), LA_BLOCFACTO)
              ENDIF
            ENDIF
            IF (run_gemmt_on_gpu .EQ. 1.AND.(JBEG_BLOCK.EQ.1)) THEN
              TMP_SIZE_PIN=TMP_SIZE_PIN + int(NPIV,8) * int(NROW1,8)
              CALL DMUMPS_GPU_PIN_MEMORY_C(
     &             UIP21K, int(NPIV,8)*int(NROW1,8))
            ENDIF
          ENDIF
#endif
        ENDIF
        ENDIF 
        IF ( NSLAVES_FOLLOW .NE. 0 .and. NPIV .NE. 0 ) THEN
          ALLOCATE( LIST_SLAVES_FOLLOW ( NSLAVES_FOLLOW ),
     &            stat = allocok )
          IF ( allocok .GT. 0 ) THEN
            IF (LP > 0 ) WRITE(LP,*) MYID,
     &": ALLOCATION FAILURE FOR LIST_SLAVES_FOLLOW
     & IN DMUMPS_PROCESS_SYM_BLOCFACTO"
            IFLAG = -13
            IERROR = NSLAVES_FOLLOW
            GOTO 700
          END IF
          LIST_SLAVES_FOLLOW(1:NSLAVES_FOLLOW)=
     &    IW(IOLDPS+6+XTRA_SLAVES_SYM+KEEP(IXSZ):
     &     IOLDPS+5+XTRA_SLAVES_SYM+KEEP(IXSZ)+NSLAVES_FOLLOW)
        END IF
       IF ( (JBEG_BLOCK.EQ.1) .AND.
     &       ((.NOT. LR_ACTIV_BlocFACTO).OR.KEEP(475).EQ.0)
     &     ) THEN
#if defined(USE_GPU)
          IF (run_trsm_on_gpu.EQ.1) THEN
            IF (DYNAMIC_ALLOC) THEN
              CALL DMUMPS_TRSM_GPU_C(0, 
     &               'L','U','T','U', NPIV, NROW1, ONE,
     &               DYN_BLOCFACTO, LD_BLOCFACTO,
     &               A_PTR(POSELT+int(NPIV1,8)), NCOL1 )
            ELSE
              CALL DMUMPS_TRSM_GPU_C(0, 
     &             'L','U','T','U', NPIV, NROW1, ONE,
     &             A( POSBLOCFACTO ), LD_BLOCFACTO,
     &             A_PTR(POSELT+int(NPIV1,8)), NCOL1)
            ENDIF
          ELSE
#endif
            IF (DYNAMIC_ALLOC) THEN
              CALL dtrsm( 'L', 'U', 'T', 'U', NPIV, NROW1, ONE,
     &             DYN_BLOCFACTO, LD_BLOCFACTO,
     &             A_PTR(POSELT+int(NPIV1,8)), NCOL1)
            ELSE
              CALL dtrsm( 'L', 'U', 'T', 'U', NPIV, NROW1, ONE,
     &             A( POSBLOCFACTO ), LD_BLOCFACTO,
     &             A_PTR(POSELT+int(NPIV1,8)), NCOL1)
            ENDIF
#if defined(USE_GPU)
          ENDIF
#endif
        ENDIF
        IF (JBEG_BLOCK.EQ.1) THEN
        IF (.NOT.LR_ACTIV_BlocFACTO.OR.KEEP(475).EQ.0) THEN
#if defined(USE_GPU)
#endif
         NB = KEEP(360)
!$       NOMP = OMP_GET_MAX_THREADS()
!$       OMP_FLAG = (NOMP.GT.1.AND. (int(NROW1/NB).GE.NOMP))
!$OMP    PARALLEL DO 
!$OMP&   PRIVATE (IB, II, IBEG, IEND, I, J, UPOS, LPOS, DPOS, 
!$OMP&            PIVI, A11, A12, A22, POSPV1, POSPV2, 
!$OMP&            OFFDAG, DETPIV, LPOS1, MULT1, MULT2
!$OMP&            )
!$OMP&   SCHEDULE(DYNAMIC,1) IF (OMP_FLAG)
         DO IB=1, NROW1, NB
          IBEG = IB
          IEND = min(IB+NB-1, NROW1)
          IF (.NOT.LR_ACTIV_BlocFACTO) THEN
            LPOS = POSELT + int(NPIV1,8)
            UPOS = 1_8
            LPOS = LPOS + int((IBEG-1),8)*int(NCOL1,8)
            UPOS = UPOS + int((IBEG-1),8)*int(NPIV,8)
            DO II = IBEG, IEND
             DO J = 0, NPIV-1
               UIP21K( UPOS+J ) = A_PTR(LPOS+J)
             ENDDO
             LPOS = LPOS + int(NCOL1,8)
             UPOS = UPOS + int(NPIV,8)
            END DO
          ENDIF
          LPOS = POSELT + int(NPIV1,8)
          LPOS = LPOS + int((IBEG-1),8)*int(NCOL1,8)
          IF (DYNAMIC_ALLOC) THEN
            DPOS = 1_8
          ELSE
            DPOS = POSBLOCFACTO
          ENDIF
          I = 1
          DO
            IF(I .GT. NPIV) EXIT
            IF (DYNAMIC_ALLOC) THEN
              PIVI = DYN_PIVINFO(I)
            ELSE
              PIVI = IW(IPIV+I-1)
            ENDIF
            IF(PIVI .GT. 0) THEN
              IF (DYNAMIC_ALLOC) THEN
                A11 = ONE/DYN_BLOCFACTO(DPOS)
              ELSE
                A11 = ONE/A(DPOS)
              ENDIF
              CALL dscal( IEND-IBEG+1, A11, A_PTR(LPOS), NCOL1 )
              LPOS = LPOS + 1_8
              DPOS = DPOS + int(LD_BLOCFACTO + 1,8)
              I = I+1
            ELSE
              POSPV1 = DPOS
              POSPV2 = DPOS+ int(LD_BLOCFACTO + 1,8)
              OFFDAG = POSPV1+1_8
              IF (DYNAMIC_ALLOC) THEN
                A11 = DYN_BLOCFACTO(POSPV1)
                A22 = DYN_BLOCFACTO(POSPV2)
                A12 = DYN_BLOCFACTO(OFFDAG)
                DETPIV = A11*A22 - A12**2
                A22 = A11/DETPIV
                A11 = DYN_BLOCFACTO(POSPV2)/DETPIV
                A12 = -A12/DETPIV
              ELSE
                A11 = A(POSPV1)
                A22 = A(POSPV2)
                A12 = A(OFFDAG)
                DETPIV = A11*A22 - A12**2
                A22 = A11/DETPIV
                A11 = A(POSPV2)/DETPIV
                A12 = -A12/DETPIV
              ENDIF
              LPOS1 = LPOS 
              DO J2 = 1, IEND-IBEG+1
                 MULT1 = A11*A_PTR(LPOS1)+A12*A_PTR(LPOS1+1_8)
                 MULT2 = A12*A_PTR(LPOS1)+A22*A_PTR(LPOS1+1_8)
                 A_PTR(LPOS1) = MULT1
                 A_PTR(LPOS1+1_8) = MULT2
                 LPOS1 = LPOS1 + int(NCOL1,8)
              ENDDO
              LPOS = LPOS + 2_8
              DPOS = POSPV2 + int(LD_BLOCFACTO + 1,8)
              I = I+2
            ENDIF
          ENDDO
         ENDDO
!$OMP    END PARALLEL DO
        ENDIF
        ENDIF 
      ENDIF 
      COMPRESS_CB = .FALSE.
      IF ( LR_ACTIVATED ) THEN
        COMPRESS_CB    = MUMPS_LR_COMPRESS_CB(IW(IOLDPS+XXLR))
      ENDIF
       IF (NROW1.GT.0) THEN
         IF (NPIV.GT.0.AND.NROW1.LE.0) THEN
            CALL MUMPS_ABORT()  
         ENDIF
       IF (LR_ACTIVATED) THEN
        IF (NPIV1.NE.0.OR.JBEG_BLOCK.NE.1) THEN
          IF (LR_ACTIV_BlocFACTO) THEN
           CALL DMUMPS_BLR_RETRIEVE_BEGS_BLR_L (IW(IOLDPS+XXF), 
     &                  BEGS_BLR_LS)
           KEEP_BEGS_BLR_LS = .TRUE.  
           NB_BLR_LS = size(BEGS_BLR_LS) - 2
           NPARTSCB  = NB_BLR_LS
          ENDIF
        ELSE
          CALL GET_CUT(IW(IOLDPS+HS:IOLDPS+HS+NROW1-1), 0,
     &                    NROW1, LRGROUPS, NPARTSCB, 
     &                    NPARTSASS, BEGS_BLR_LS)
              CALL REGROUPING2(BEGS_BLR_LS, NPARTSASS, NASS1, NPARTSCB,
     &                        NROW1-0, KEEP(488), .TRUE., KEEP(472),
     &                        NCOL1, KEEP(1))
             NB_BLR_LS = NPARTSCB
        ENDIF
        IF (LR_ACTIV_BlocFACTO) THEN
          IF (NPIV.GT.0) THEN
           call MAX_CLUSTER(BEGS_BLR_LM(2:NB_BLR_LM+2),NB_BLR_LM,
     &                            MAXI_CLUSTER_LM)
          ELSE
           MAXI_CLUSTER_LM = 0
          ENDIF
          call MAX_CLUSTER(BEGS_BLR_LS,NB_BLR_LS+1,MAXI_CLUSTER_LS)
          MAXI_CLUSTER=max(MAXI_CLUSTER_LS,MAXI_CLUSTER_LM,NPIV)
        ENDIF
        IF (NPIV1.EQ.0.AND.(JBEG_BLOCK.EQ.1))  THEN
          INFO_TMP(1) = IFLAG
          INFO_TMP(2) = IERROR
          NB_ACCESSES_INIT = 1
          IF ( (KEEP(486).EQ.2) 
     &       ) THEN
            NB_ACCESSES_INIT = huge(NPARTSASS_MASTER)
          END IF
          INFO_TMP(1) = IFLAG
          INFO_TMP(2) = IERROR
          IF (IFLAG.LT.0) GOTO 700
          CALL DMUMPS_BLR_SAVE_INIT(IW(IOLDPS+XXF), 
     &              .TRUE., 
     &              .TRUE., 
     &              .TRUE., 
     &              NPARTSASS_COL, 
     &              BEGS_BLR_LS, 
     &              BEGS_BLR_COL,  
     &              NB_ACCESSES_INIT, 
     &              INFO_TMP)
          IF (associated(BEGS_BLR_COL)) THEN 
            DEALLOCATE(BEGS_BLR_COL)
            NULLIFY(BEGS_BLR_COL)
          ENDIF
          IFLAG  = INFO_TMP(1) 
          IERROR = INFO_TMP(2) 
          IF (IFLAG.LT.0) GOTO 700
        ENDIF
        IF (LR_ACTIV_BlocFACTO) THEN
         LWORK = MAXI_CLUSTER*(MAXI_CLUSTER+1)
         ALLOCATE(BLOCKLR(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &            RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &            TAU(MAXI_CLUSTER*OMP_NUM),
     &            JPVT(MAXI_CLUSTER*OMP_NUM), 
     &            WORK(LWORK*OMP_NUM),
     &            stat=allocok)
         IF (allocok > 0 ) THEN
           IFLAG  = -13
           IERROR = OMP_NUM*(LWORK + MAXI_CLUSTER*(MAXI_CLUSTER+4))
           GOTO 700
         ENDIF
         CURRENT_BLR = 1
         IF (JBEG_BLOCK.EQ.1.AND.NPIV.GT.0) THEN
           CURRENT_BLR = 1
           ALLOCATE(BLR_LS(NB_BLR_LS), stat=allocok)
           IF (allocok > 0 ) THEN
             IFLAG  = -13
             IERROR = NB_BLR_LS
             GOTO 700
           ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL
#endif
           CALL DMUMPS_COMPRESS_PANEL_I_NOOPT
     &        (A_PTR(POSELT), LA_PTR, 1_8,
     &        IFLAG, IERROR, NCOL1,
     &        BEGS_BLR_LS(1), size(BEGS_BLR_LS), NB_BLR_LS+1,
     &        DKEEP(8), KEEP(466), 0, 
     &        KEEP(473), BLR_LS(1), 
     &        CURRENT_BLR, 'V', WORK, TAU, JPVT, LWORK, RWORK,
     &        BLOCKLR, MAXI_CLUSTER, NELIM, 
     &        .TRUE.,  
     &        NPIV, NPIV1,
     &        2, KEEP(483), KEEP8, 
     &        (KEEP(450).GT.0), KEEP(453), KEEP(454).EQ.2,
     &        OMP_NUM)
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif
           IF (IFLAG.LT.0) GOTO 300
           IF (KEEP(475).GE.1) THEN
             IF (DYNAMIC_ALLOC) THEN
               CALL DMUMPS_BLR_PANEL_LRTRSM(
     &               DYN_BLOCFACTO, LA_BLOCFACTO, 1_8,
     &               LD_BLOCFACTO, -6666, 
     &               NB_BLR_LS+1, 
     &               BLR_LS, CURRENT_BLR, CURRENT_BLR+1, NB_BLR_LS+1, 
     &               2, 1, 0, 
     &               .TRUE., 
     &               DYN_PIVINFO, OFFSET_IW=1)
             ELSE
             CALL DMUMPS_BLR_PANEL_LRTRSM(A, LA, POSBLOCFACTO, 
     &               LD_BLOCFACTO, -6666, 
     &               NB_BLR_LS+1, 
     &               BLR_LS, CURRENT_BLR, CURRENT_BLR+1, NB_BLR_LS+1, 
     &               2, 1, 0, 
     &               .TRUE., 
     &               IW, OFFSET_IW=IPIV)
             ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif          
             IF (KEEP(486).NE.2) THEN
               CALL DMUMPS_DECOMPRESS_PANEL_I_NOOPT(
     &         A_PTR(POSELT), LA_PTR, 1_8,
     &         NCOL1, NCOL1,
     &         .TRUE.,    
     &         NPIV1+1,   
     &         1,         
     &         NB_BLR_LS+1, BLR_LS(1), CURRENT_BLR, 'V', 1)
             ENDIF
          ENDIF
 300      CONTINUE         
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif          
          IF (IFLAG.LT.0) GOTO 700
         ENDIF 
        ENDIF 
       ENDIF  
       ENDIF   
      IF ( (KEEP(201).eq.1) .AND. 
     &    (OOCWRITE_COMPATIBLE_WITH_BLR .OR. NPIV.EQ.0) 
     &    .AND. (JBEG_BLOCK.EQ.1) 
     &   ) THEN
        MonBloc%INODE = INODE
        MonBloc%MASTER = .FALSE.
        MonBloc%Typenode = 2
        MonBloc%NROW = NROW1  
        MonBloc%NCOL = NCOL1  
        MonBloc%NFS  = NASS1
        MonBloc%LastPiv = NPIV1 + NPIV 
        MonBloc%LastPanelWritten_L = -9999 
        MonBloc%LastPanelWritten_U = -9999 
        NULLIFY(MonBloc%INDICES)
        MonBloc%Last = LASTPANEL 
        STRAT = STRAT_TRY_WRITE 
        NextPivDummy      = -8888 
        LIWFAC = IW(IOLDPS+XXI)
        LAST_CALL = .FALSE.
        CALL DMUMPS_OOC_IO_LU_PANEL_I( STRAT, TYPEF_L,
     &       A_PTR(POSELT),
     &       LA_PTR, MonBloc, NextPivDummy, NextPivDummy,
     &       IW(IOLDPS), LIWFAC, MYID, KEEP8(31), IFLAG_OOC, LAST_CALL)
        IF ( IFLAG_OOC .LT. 0 )THEN
          IFLAG = IFLAG_OOC
          GOTO 700
        ENDIF
      ENDIF
      IF (NPIV.GT.0) THEN
        IF (LR_ACTIV_BlocFACTO) THEN
          IF (JBEG_BLOCK.NE.1) THEN
            CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &           IW(IOLDPS+XXF),
     &           0, 
     &           IPANEL, BLR_LS)
            KEEP_BLR_LS = .TRUE.  
          ENDIF
          IF (NELIM.GT.0.AND.JBEG_BLOCK.EQ.1) THEN
            LPOS2 = POSELT + int(NPIV1,8)
            UPOS = 1_8+int(SHIFT_UPOS,8)
            LPOS  = LPOS2 + int(SHIFT_LPOS,8)
            IF (DYNAMIC_ALLOC) THEN
              CALL DMUMPS_BLR_UPD_NELIM_VAR_L_I(
     &        DYN_BLOCFACTO, LA_BLOCFACTO, UPOS,
     &        A_PTR(POSELT), LA_PTR, LPOS-POSELT+1_8,
     &        IFLAG, IERROR, LD_BLOCFACTO, NCOL1,
     &        BEGS_BLR_LS(1), size(BEGS_BLR_LS),
     &        CURRENT_BLR, BLR_LS(1), NB_BLR_LS+1, 
     &        CURRENT_BLR+1, NELIM, 'N')
            ELSE
              CALL DMUMPS_BLR_UPD_NELIM_VAR_L_I(
     &        A(POSBLOCFACTO), LA_BLOCFACTO, UPOS,
     &        A_PTR(POSELT), LA_PTR, LPOS-POSELT+1_8,
     &        IFLAG, IERROR, LD_BLOCFACTO, NCOL1,
     &        BEGS_BLR_LS(1), size(BEGS_BLR_LS),
     &        CURRENT_BLR, BLR_LS(1), NB_BLR_LS+1, 
     &        CURRENT_BLR+1, NELIM, 'N')
            ENDIF
          ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL
#endif          
          IF (DYNAMIC_ALLOC) THEN
            CALL DMUMPS_BLR_SLV_UPD_TRAIL_LDLT_I(
     &        A_PTR(POSELT), LA_PTR, 1_8, 
     &        IFLAG, IERROR, NCOL1, NROW1, JBEG_BLOCK, 
     &        DYN_BLOCFACTO, LA_BLOCFACTO,
     &        LD_BLOCFACTO, 
     &        BEGS_BLR_LM(1), size(BEGS_BLR_LM), NB_BLR_LM+1,
     &        BLR_LM(1), NPIV1, 
     &        BEGS_BLR_LS(1), size(BEGS_BLR_LS), NB_BLR_LS+1,
     &        BLR_LS(1), 0, 
     &        CURRENT_BLR, CURRENT_BLR,   
     &        DYN_PIVINFO, 
     &        BLOCKLR,
     &        MAXI_CLUSTER, OMP_NUM,
     &        KEEP(481), DKEEP(11), KEEP(466), KEEP(477) 
     &        )
          ELSE
            CALL DMUMPS_BLR_SLV_UPD_TRAIL_LDLT_I(
     &        A_PTR(POSELT), LA_PTR, 1_8, 
     &        IFLAG, IERROR, NCOL1, NROW1, JBEG_BLOCK, 
     &        A(POSBLOCFACTO), LA_BLOCFACTO,
     &        LD_BLOCFACTO, 
     &        BEGS_BLR_LM(1), size(BEGS_BLR_LM), NB_BLR_LM+1,
     &        BLR_LM(1), NPIV1, 
     &        BEGS_BLR_LS(1), size(BEGS_BLR_LS), NB_BLR_LS+1,
     &        BLR_LS(1), 0, 
     &        CURRENT_BLR, CURRENT_BLR,   
     &        IW(IPIV), 
     &        BLOCKLR,
     &        MAXI_CLUSTER, OMP_NUM,
     &        KEEP(481), DKEEP(11), KEEP(466), KEEP(477) 
     &        )
          ENDIF
          IF (IFLAG.LT.0) GOTO 400
 400      CONTINUE          
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif          
          IF (IFLAG.LT.0) GOTO 700
          IF (LASTBL_INPANEL) THEN 
            CALL UPD_MRY_LU_LRGAIN(BLR_LS, NPARTSCB
     &               , KEEP(450) .GT. 0
     &               )
          ENDIF
          CALL DEALLOC_BLR_PANEL(BLR_LM, NB_BLR_LM, KEEP8, KEEP(34))
          DEALLOCATE(BLR_LM)
          IF ( JBEG_BLOCK.EQ.1
     &       )
     &     THEN
            IF ( (KEEP(486).EQ.2) 
     &       ) THEN
                NB_ACCESSES_LEFT_INIT = huge(NB_ACCESSES_LEFT_INIT)
            ELSE
                NB_ACCESSES_LEFT_INIT = NCOL1 - NPIV1 - NROW1
            ENDIF
            CALL DMUMPS_BLR_SAVE_PANEL_LORU(
     &          IW(IOLDPS+XXF),
     &          0,   
     &          IPANEL, BLR_LS, NB_ACCESSES_LEFT_INIT)
             KEEP_BLR_LS = .TRUE.
          ENDIF
        ELSE 
          IF (NPIV .GT. 0 .AND. NCOL_GEMM_FR.GT.0)THEN
          LPOS2 = POSELT + int(NPIV1,8)
          LPOS  = LPOS2 + int(SHIFT_LPOS,8)
#if defined(USE_GPU)
          IF (run_gemm_on_gpu .EQ.1) THEN
            IF (KEEP(427) .LT.0 .AND. run_gemmt_on_gpu .EQ.1) THEN
              async_gemm=1
            ELSE 
              async_gemm=0
            ENDIF
            IF (DYNAMIC_ALLOC) THEN
              UPOS = 1_8+int(SHIFT_UPOS,8)
              CALL DMUMPS_GEMM_GPU_C(async_gemm,
     &             'N','N', NCOL_GEMM_FR, NROW1, NPIV,
     &             ALPHA, DYN_BLOCFACTO(UPOS), NEWCOL_RECV,
     &             A_PTR(LPOS2), NCOL1, ONE, A_PTR(LPOS), NCOL1)
            ELSE
              UPOS = POSBLOCFACTO+int(SHIFT_UPOS,8)
              CALL DMUMPS_GEMM_GPU_C(async_gemm,
     &             'N','N', NCOL_GEMM_FR, NROW1, NPIV,
     &             ALPHA,A(UPOS),NEWCOL_RECV,
     &             A_PTR(LPOS2), NCOL1, ONE, A_PTR(LPOS),NCOL1)
            ENDIF
          ELSE
#endif
            IF (DYNAMIC_ALLOC) THEN
              UPOS = 1_8+int(SHIFT_UPOS,8)
              CALL dgemm('N','N', NCOL_GEMM_FR, NROW1, NPIV,
     &             ALPHA, DYN_BLOCFACTO(UPOS), NEWCOL_RECV,
     &             A_PTR(LPOS2), NCOL1, ONE, A_PTR(LPOS), NCOL1)
            ELSE
              UPOS = POSBLOCFACTO+int(SHIFT_UPOS,8)
              CALL dgemm('N','N', NCOL_GEMM_FR, NROW1, NPIV,
     &             ALPHA,A(UPOS), NEWCOL_RECV,
     &             A_PTR(LPOS2), NCOL1, ONE, A_PTR(LPOS), NCOL1)
            ENDIF
#if defined(USE_GPU)
          ENDIF
#endif
          ENDIF
          IF (JBEG_BLOCK.EQ.1) THEN
          DPOS = POSELT + int(NCOL1 - NROW1,8)
#if defined(USE_GPU)
          IF (run_gemmt_on_gpu.EQ.1) THEN
            LPOS2 = POSELT + int(NPIV1,8)
            UPOS  = 1_8
#if defined(USE_XKBLAS)
            CALL DMUMPS_GEMMT_GPU_C(0, 'U', 'T' , 'N' , NROW1, NPIV,
     &           ALPHA,
     &           UIP21K( UPOS ), NPIV,
     &           A_PTR( LPOS2 ), NCOL1, ONE,
     &           A_PTR( DPOS ), NCOL1 )
#else
            CALL DMUMPS_SYRKX_GPU_C(0, 'U', 'T', NROW1, NPIV,
     &           ALPHA,
     &           UIP21K( UPOS ), NPIV,
     &           A_PTR( LPOS2 ), NCOL1, ONE,
     &           A_PTR( DPOS ), NCOL1 )
#endif
          ELSE
#endif
#if defined(GEMMT_AVAILABLE)
            IF ( KEEP(421).EQ. -1 .OR.
     &        ( KEEP(421) .GT. 0 .AND. NROW1 .GT. KEEP(421) ) ) THEN
              LPOS2 = POSELT + int(NPIV1,8)
              UPOS  = 1_8
              CALL dgemmt( 'U', 'T', 'N', NROW1, NPIV, ALPHA,
     &         UIP21K( UPOS ), NPIV,
     &         A_PTR( LPOS2 ), NCOL1, ONE,
     &         A_PTR( DPOS ), NCOL1 )
            ELSE
#endif
              IF ( NROW1 .GT. KEEP(7) ) THEN
                BLSIZE = KEEP(8)
              ELSE
                BLSIZE = NROW1
              ENDIF
              IF ( NROW1 .GT. 0 ) THEN
                DO IROW = 1, NROW1, BLSIZE
                  Block = min( BLSIZE, NROW1 - IROW + 1 )
                  DPOS  = POSELT + int(NCOL1 - NROW1,8)
     &                + int( IROW - 1, 8 ) * int( NCOL1 + 1, 8 )
                  LPOS2 = POSELT + int(NPIV1,8)
     &                + int( IROW - 1, 8 ) * int( NCOL1, 8 )
                  UPOS  = int( IROW - 1, 8 ) * int(NPIV, 8) + 1_8
                  DO I = 1, Block
                  CALL dgemv( 'T', NPIV, Block-I+1, ALPHA,
     &             A_PTR( LPOS2 + int(I - 1,8) * int(NCOL1,8) ), NCOL1,
     &             UIP21K( UPOS + int(NPIV,8) * int( I - 1, 8 ) ),
     &             1, ONE, A_PTR(DPOS+int(NCOL1+1,8)*int(I-1,8)),NCOL1 )
                  END DO
                  IF ( NROW1-IROW+1-Block .ne. 0 )
     &            CALL dgemm( 'T', 'N', Block, NROW1-IROW+1-Block,
     &            NPIV, ALPHA,
     &            UIP21K( UPOS ), NPIV,
     &            A_PTR( LPOS2 + int(Block,8) * int(NCOL1,8) ), NCOL1,
     &            ONE,
     &            A_PTR( DPOS + int(Block,8) * int(NCOL1,8) ), NCOL1 )
                ENDDO
              ENDIF
#if defined(GEMMT_AVAILABLE)
            ENDIF 
#endif
#if defined(USE_GPU)
          ENDIF
#endif
          ENDIF  
#if defined(USE_GPU)
          IF ( KEEP(422).NE.0 .AND. KEEP(423).EQ.1. AND.
     &        ( run_gemmt_on_gpu.EQ.1 .OR. run_gemm_on_gpu.EQ.1 ) ) THEN
              TMP_SIZE_UNPIN = 0_8
              IF (run_gemm_on_gpu .EQ. 1) THEN
                TMP_SIZE_UNPIN = TMP_SIZE_UNPIN + LA_BLOCFACTO
                IF (DYNAMIC_ALLOC) THEN
                  CALL DMUMPS_GPU_UNPIN_MEMORY_C(
     &                 DYN_BLOCFACTO(1), LA_BLOCFACTO)
                ELSE
                  CALL DMUMPS_GPU_UNPIN_MEMORY_C(
     &                 A(POSBLOCFACTO),LA_BLOCFACTO)
                ENDIF
              ENDIF
              IF (run_gemmt_on_gpu .EQ. 1.AND. JBEG_BLOCK.EQ.1) THEN
                TMP_SIZE_UNPIN = TMP_SIZE_UNPIN +
     &                           int(NPIV,8)*int(NROW1,8)
                CALL DMUMPS_GPU_UNPIN_MEMORY_C(
     &               UIP21K(1),int(NPIV,8)*int(NROW1,8))
              ENDIF
          ENDIF
#endif
        ENDIF 
        IF (LASTBL_INPANEL) THEN
          FLOP1 = dble(NROW1) * dble(NPIV) *
     &           dble( 2 * (NASS1-NPIV1)  - NPIV + NROW1 +1 )
          FLOP1 = -FLOP1
          CALL MUMPS_LOAD_UPDATE( 1, .FALSE., FLOP1, KEEP,KEEP8 )
        ENDIF
      ENDIF 
      IF  (JBEG_BLOCK.EQ.1) THEN
        IW(IOLDPS+KEEP(IXSZ))   = IW(IOLDPS+KEEP(IXSZ)) - NPIV   
        IW(IOLDPS+3+KEEP(IXSZ)) = IW(IOLDPS+3+KEEP(IXSZ)) + NPIV 
      ENDIF
      IF (LASTBL_INLASTPANEL) THEN
        IW(IOLDPS+1+KEEP(IXSZ)) = IW(IOLDPS+3+KEEP(IXSZ))
      ENDIF
      IF ( .NOT. LR_ACTIV_BlocFACTO ) THEN
      IF (DYNAMIC_ALLOC) THEN
        IF (allocated(DYN_PIVINFO)  ) DEALLOCATE(DYN_PIVINFO)
        IF (allocated(DYN_BLOCFACTO)) THEN
           KEEP8(130) = KEEP8(130)-max(1_8,LA_BLOCFACTO)
           DEALLOCATE(DYN_BLOCFACTO)
           KEEP8(69) = KEEP8(69) - max(1_8,LA_BLOCFACTO)
           KEEP8(73) = KEEP8(73) - max(1_8,LA_BLOCFACTO)
        ENDIF
      ELSE
      LRLU  = LRLU + LA_BLOCFACTO
      LRLUS = LRLUS + LA_BLOCFACTO
      KEEP8(69) = KEEP8(69) - LA_BLOCFACTO
      POSFAC = POSFAC - LA_BLOCFACTO
      IWPOS = IWPOS - NPIV
      CALL MUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &                LA-LRLUS,0_8,-LA_BLOCFACTO,KEEP,KEEP8,LRLUS)
      ENDIF 
      ENDIF 
      IF (LR_ACTIV_BlocFACTO) THEN
        IF (allocated(RWORK))  DEALLOCATE(RWORK)
        IF (allocated(WORK)) DEALLOCATE(WORK)
        IF (allocated(TAU)) DEALLOCATE(TAU)
        IF (allocated(JPVT)) DEALLOCATE(JPVT)
        IF (allocated(BLOCKLR)) DEALLOCATE(BLOCKLR)
      ENDIF
      IF ( NSLAVES_FOLLOW .NE. 0 .and. NPIV .NE. 0 
     &     .AND. JBEG_BLOCK.EQ.1 ) THEN
         IPOSK = NPIV1 + 1
         JPOSK = NCOL1 - NROW1 + 1
           NPIVSENT = NPIV
           BLFAC_NBCOLS_ALREADY_SENT = 0
           BLFAC_NBLRB_ALREADY_SENT  = 0  
           IERR = -1
           DO WHILE ( IERR .eq. -1 )
            IF (DYNAMIC_ALLOC) THEN
              CALL DMUMPS_BUF_SEND_BLFAC_SLAVE(
     &                    INODE, NPIVSENT, FPERE, 
     &                    IPOSK, JPOSK,
     &                    UIP21K, LUIP21K, NROW1,
     &                    NSLAVES_FOLLOW,
     &                    LIST_SLAVES_FOLLOW(1),
     &                    COMM, KEEP,
     &           LR_ACTIV_BlocFACTO, BLR_LS, IPANEL, 
     &           BLFAC_NBCOLS_ALREADY_SENT, BLFAC_NBLRB_ALREADY_SENT,
     &           NOTHING_WAS_SENT,
     &           DYN_BLOCFACTO, LA_BLOCFACTO,
     &           1_8, LD_BLOCFACTO,
     &           DYN_PIVINFO, MAXI_CLUSTER,
     &                  IERR, IERROR )
            ELSE
              CALL DMUMPS_BUF_SEND_BLFAC_SLAVE(
     &                    INODE, NPIVSENT, FPERE,
     &                    IPOSK, JPOSK,
     &                    UIP21K, LUIP21K, NROW1,
     &                    NSLAVES_FOLLOW,
     &                    LIST_SLAVES_FOLLOW(1),
     &                    COMM, KEEP,
     &           LR_ACTIV_BlocFACTO, BLR_LS, IPANEL, 
     &           BLFAC_NBCOLS_ALREADY_SENT, BLFAC_NBLRB_ALREADY_SENT,
     &           NOTHING_WAS_SENT,
     &           A, LA, 
     &           POSBLOCFACTO, LD_BLOCFACTO,
     &           IW(IPIV), MAXI_CLUSTER,
     &                  IERR, IERROR )
            ENDIF
            IF (IERR.EQ.-13) THEN
             IFLAG  = IERR
             IF (LP > 0 ) WRITE(LP,*) MYID,
     &     ": ALLOCATION FAILURE within DMUMPS_BUF_SEND_BLFAC_SLAVE",
     &     " during DMUMPS_PROCESS_SYM_BLOCFACTO", IERROR
             GOTO 700
            ENDIF
            IF (IERR .EQ. -1 .AND. NOTHING_WAS_SENT) THEN
              IOLDPS = PTRIST(STEP(INODE))
              IF ( IW(IOLDPS+6+KEEP(IXSZ)) .eq.
     &              huge(IW(IOLDPS+6+KEEP(IXSZ))) ) THEN
                    COUNTER_WAS_HUGE=.TRUE.
                    IW(IOLDPS+6+KEEP(IXSZ)) = 1
              ELSE
                    COUNTER_WAS_HUGE=.FALSE.
              ENDIF
              TO_UPDATE_CPT_RECUR =
     &                      ( NSLAVES_TOT - NSLAVES_FOLLOW - 1 ) *
     &                       (2*NASS1/KEEP(6))
              IW(IOLDPS+6+KEEP(IXSZ)) =
     &             IW(IOLDPS+6+KEEP(IXSZ)) - TO_UPDATE_CPT_RECUR - 10
              BLOCKING = .FALSE.
              SET_IRECV= .TRUE.
              MESSAGE_RECEIVED = .FALSE.
              CALL DMUMPS_TRY_RECVTREAT( COMM_LOAD, ASS_IRECV,
     &         BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &         MPI_ANY_SOURCE, MPI_ANY_TAG,
     &         STATUS, 
     &         BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &         IWPOS, IWPOSCB, IPTRLU,
     &         LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &         PTLUST_S, PTRFAC,
     &         PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &         IFLAG, IERROR, COMM,
     &         PERM, IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &         root, roota, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &         FILS, DAD, PTRARW, PTRAIW,
     &         PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &         INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &         LPTRAR, NELT, FRTPTR, FRTELT, 
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &               , LRGROUPS
     &           )
              IOLDPS = PTRIST(STEP(INODE))
              IW(IOLDPS+6+KEEP(IXSZ)) =
     &             IW(IOLDPS+6+KEEP(IXSZ)) + TO_UPDATE_CPT_RECUR + 10
              IF ( COUNTER_WAS_HUGE .AND.
     &             IW(IOLDPS+6+KEEP(IXSZ)).EQ.1 ) THEN
                IW(IOLDPS+6+KEEP(IXSZ)) = huge(IW(IOLDPS+6+KEEP(IXSZ)))
              ENDIF
              IF ( IFLAG .LT. 0 ) GOTO 550
            END IF
           END DO
           IF ( IERR .eq. -2 ) THEN
             IF (LP > 0 ) THEN 
               WRITE(LP,*) MYID,
     &": FAILURE, SEND BUFFER TOO SMALL DURING
     & DMUMPS_PROCESS_SYM_BLOCFACTO"
             ENDIF
             IFLAG = -17
             IERROR = 5 * KEEP(34) + NPIV * NROW1 * KEEP(35)
             GOTO 700
           END IF
           IF ( IERR .eq. -3 ) THEN
              IF (LP > 0 ) WRITE(LP,*) MYID,
     &": FAILURE, RECV BUFFER TOO SMALL DURING
     & DMUMPS_PROCESS_SYM_BLOCFACTO"
             IFLAG = -20
             IERROR = 5 * KEEP(34) + NPIV * NROW1 * KEEP(35)
             GOTO 700
           END IF
           DEALLOCATE(LIST_SLAVES_FOLLOW)
      END IF
      IF ( LR_ACTIV_BlocFACTO ) THEN
        IF (NPIV.GT.0 
     &     .AND. KEEP(486).EQ.3
     &    ) THEN 
          IOLDPS = PTRIST(STEP(INODE))
          CALL DMUMPS_BLR_DEC_AND_TRYFREE_L(IW(IOLDPS+XXF),IPANEL,
     &                       KEEP8, KEEP(34), NEWCOL_RECV)
        ENDIF 
        IF (DYNAMIC_ALLOC) THEN
          IF (allocated(DYN_PIVINFO))   DEALLOCATE(DYN_PIVINFO)
          IF (allocated(DYN_BLOCFACTO)) THEN
             KEEP8(130) = KEEP8(130)-max(1_8,LA_BLOCFACTO)
             DEALLOCATE(DYN_BLOCFACTO)
             KEEP8(69) = KEEP8(69) - max(1_8,LA_BLOCFACTO)
             KEEP8(73) = KEEP8(73) - max(1_8,LA_BLOCFACTO)
          ENDIF
        ELSE IF (NPIV .GT. 0) THEN
          LRLU  = LRLU + LA_BLOCFACTO
          LRLUS = LRLUS + LA_BLOCFACTO
          KEEP8(69) = KEEP8(69) - LA_BLOCFACTO
          POSFAC = POSFAC - LA_BLOCFACTO
          IWPOS = IWPOS - NPIV
      CALL MUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &             LA-LRLUS,0_8,-LA_BLOCFACTO,KEEP,KEEP8,LRLUS)
        ENDIF 
      ENDIF   
      IF ( NPIV .NE. 0 )  THEN
        IF (associated(UIP21K)) THEN
          CALL DMUMPS_DM_FREE_S_WK( UIP21K, KEEP430_LOC )
          NULLIFY( UIP21K )
          KEEP8(130) = KEEP8(130)-LUIP21K
          KEEP8(69) = KEEP8(69) - LUIP21K
          KEEP8(73) = KEEP8(73) - LUIP21K
        ENDIF
      ENDIF
      IOLDPS = PTRIST(STEP(INODE))
      IF (LR_ACTIVATED ) THEN
        COMPRESS_CB    = MUMPS_LR_COMPRESS_CB(IW(IOLDPS+XXLR))
      ENDIF
      CALL DMUMPS_DM_SET_DYNPTR( IW(IOLDPS+XXS), A, LA,
     &     PTRAST(STEP(INODE)), IW(IOLDPS+XXD), IW(IOLDPS+XXR),
     &     A_PTR, POSELT, LA_PTR )
      IF (LASTBL_INLASTPANEL) THEN
         IF ( KEEP(486) .NE. 0) THEN
           IF (LR_ACTIV_BlocFACTO) THEN
             CALL STATS_COMPUTE_FLOP_SLAVE_TYPE2(NROW1, NCOL1, NASS1,
     &             KEEP(50), INODE)
           ELSE
             CALL UPD_FLOP_FRFRONT_SLAVE(NROW1, NCOL1, NASS1,
     &             KEEP(50), INODE)
           ENDIF
         ENDIF
         IF ( IW(IOLDPS+6+KEEP(IXSZ)).EQ.
     &     huge(IW(IOLDPS+6+KEEP(IXSZ))) ) THEN
           IW(IOLDPS+6+KEEP(IXSZ)) =  1
         ENDIF
         IW(IOLDPS+6+KEEP(IXSZ)) = IW(IOLDPS+6+KEEP(IXSZ))
     &                           - TO_UPDATE_CPT_END 
     &                           - 1 
      END IF
      IF ( IW(IOLDPS+6+KEEP(IXSZ) ) .eq. 0
     &       .and. KEEP(50) .ne. 0 .and. NSLAVES_FOLLOW .eq. 0
     &       .and. NSLAVES_TOT.NE.1 ) THEN
          DEST = MUMPS_PROCNODE( PROCNODE_STEPS(STEP(INODE)),
     &                           KEEP(199) )
          CALL MUMPS_BUF_SEND_1INT( INODE, DEST, END_NIV2_LDLT,
     &                              COMM, KEEP, IERR )
          IF ( IERR .LT. 0 ) THEN
            write(*,*) ' Internal error in PROCESS_SYM_BLOCFACTO.'
            IFLAG = -99
            GOTO 700
          END IF
      ENDIF
      IOLDPS = PTRIST(STEP(INODE))
      IF (IW(IOLDPS+6+KEEP(IXSZ)) .eq. 0 
     &   ) THEN 
          NELIM =  IW( IOLDPS + 4 + KEEP(IXSZ)) - 
     &              IW( IOLDPS + 3  + KEEP(IXSZ))
          NPIV1  = IW( IOLDPS + 3 +KEEP(IXSZ))
          IF (LR_ACTIVATED) THEN
            MAXI_CLUSTER    = 0
            MAXI_CLUSTER_LS = 0
            IF (LR_ACTIV_POSTPROC.OR. COMPRESS_CB) THEN
             IF (.not.associated(BEGS_BLR_LS)) THEN
               CALL DMUMPS_BLR_RETRIEVE_BEGS_BLR_L (IW(IOLDPS+XXF), 
     &                     BEGS_BLR_LS)
               KEEP_BEGS_BLR_LS = .TRUE.
             ENDIF
             NB_BLR_LS        = size(BEGS_BLR_LS) - 2
             call MAX_CLUSTER(BEGS_BLR_LS,NB_BLR_LS+1,
     &            MAXI_CLUSTER_LS)
              CALL DMUMPS_BLR_RETRIEVE_BEGS_BLR_C (IW(IOLDPS+XXF), 
     &                    BEGS_BLR_COL, NPARTSASS_COL )
              KEEP_BEGS_BLR_COL = .TRUE.  
              NB_BLR_COL   = size(BEGS_BLR_COL) - 1
            ENDIF
            IF (LR_ACTIV_POSTPROC) THEN
              NPARTSASS = 1
              NPARTSCB  = NB_BLR_LS
              DO IP=1,NPARTSASS_MASTER
                  IF ( BEGS_BLR_COL(IP).GT. NPIV1 ) THEN
                    BEGS_BLR_COL(IP)   = NPIV1+1
                  ENDIF
              ENDDO
              BEGS_BLR_COL(NPARTSASS_MASTER+1) = NPIV1+1
              call MAX_CLUSTER(BEGS_BLR_COL,NPARTSASS_MASTER,
     &           MAXI_CLUSTER_COL)
              MAXI_CLUSTER = max(MAXI_CLUSTER_COL, MAXI_CLUSTER_LS)
              IF (allocated(BLOCKLR)) DEALLOCATE(BLOCKLR)
              IF (allocated(RWORK))  DEALLOCATE(RWORK)
              IF (allocated(TAU)) DEALLOCATE(TAU)
              IF (allocated(JPVT)) DEALLOCATE(JPVT)
              IF (allocated(WORK)) DEALLOCATE(WORK)
              LWORK = MAXI_CLUSTER*(MAXI_CLUSTER+1)
              ALLOCATE(BLOCKLR(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &                  RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &                  TAU(MAXI_CLUSTER*OMP_NUM),
     &                  JPVT(MAXI_CLUSTER*OMP_NUM), 
     &                  WORK(LWORK*OMP_NUM),
     &                  stat=allocok)
              IF (allocok > 0 ) THEN
                IFLAG  = -13
                IERROR = OMP_NUM*(LWORK + MAXI_CLUSTER*(MAXI_CLUSTER+4))
                GOTO 700
              ENDIF
              DO IP=1,NPARTSASS_MASTER
                IF (BEGS_BLR_COL(IP+1)-BEGS_BLR_COL(IP).EQ.0) THEN
                  NULLIFY(BLR_LS)
                ELSE 
                  ALLOCATE(BLR_LS(NB_BLR_LS),stat=allocok)
                  IF (allocok > 0) THEN
                    IFLAG  = -13
                    IERROR = NB_BLR_LS
                    GOTO 700
                  ENDIF
                  DO J=1, NB_BLR_LS
                    BLR_LS(J)%M=0
                    BLR_LS(J)%N=0
                    BLR_LS(J)%K=0
                    BLR_LS(J)%ISLR=.FALSE.
                    NULLIFY(BLR_LS(J)%Q)
                    NULLIFY(BLR_LS(J)%R)
                  ENDDO
                ENDIF
                CALL DMUMPS_BLR_SAVE_PANEL_LORU (
     &              IW(IOLDPS+XXF),
     &              0, 
     &              IP, BLR_LS, 0) 
              ENDDO
              CURRENT_BLR = 1 
              KEEP_BLR_LS = .TRUE.  
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL
!$OMP& PRIVATE(IP, BLR_LS, ISHIFT)
#endif
              DO IP=1,NPARTSASS_MASTER
                IF (BEGS_BLR_COL(IP+1)-BEGS_BLR_COL(IP).EQ.0) CYCLE
                CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &            IW(IOLDPS+XXF), 0,  
     &            IP, BLR_LS)
                ISHIFT = BEGS_BLR_COL(IP)-1
                CALL DMUMPS_COMPRESS_PANEL_I_NOOPT
     &            (A_PTR(POSELT), LA_PTR, 1_8,
     &            IFLAG, IERROR, NCOL1,
     &            BEGS_BLR_LS(1), size(BEGS_BLR_LS), NB_BLR_LS+1,
     &            DKEEP(8), KEEP(466), 0, 
     &            KEEP(473), BLR_LS(1), 
     &            1, 'V', WORK, TAU, JPVT, LWORK, RWORK,
     &            BLOCKLR, MAXI_CLUSTER, NELIM, 
     &            .TRUE.,  
     &            BEGS_BLR_COL(IP+1)-BEGS_BLR_COL(IP), ISHIFT,
     &            2, KEEP(483), KEEP8,
     &            (KEEP(450).GT.0), KEEP(453), KEEP(454).EQ.2,
     &            OMP_NUM )
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif
                  IF (IFLAG.LT.0) GOTO 440
              ENDDO
 440          CONTINUE
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif          
              IF (IFLAG.GE.0) THEN
                DO IP=1,NPARTSASS_MASTER
                  IF (BEGS_BLR_COL(IP+1)-BEGS_BLR_COL(IP).EQ.0) CYCLE
                  CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &              IW(IOLDPS+XXF), 0,  
     &              IP, BLR_LS)
                  CALL UPD_MRY_LU_LRGAIN(BLR_LS, NB_BLR_LS
     &               , KEEP(450) .GT. 0
     &               )
                ENDDO
              ENDIF
            ENDIF
            IF (COMPRESS_CB) THEN
              allocate(CB_LRB(NB_BLR_LS,NB_BLR_COL-NPARTSASS_COL),
     &                 stat=allocok)
              IF (allocok > 0) THEN
                IFLAG  = -13
                IERROR = NB_BLR_LS*(NB_BLR_COL-NPARTSASS_COL)
                GOTO 700
              ENDIF
              DO II=1,NB_BLR_LS
              DO JJ=1,NB_BLR_COL-NPARTSASS_COL
                CB_LRB(II,JJ)%M=0
                CB_LRB(II,JJ)%N=0
                NULLIFY(CB_LRB(II,JJ)%Q)
                NULLIFY(CB_LRB(II,JJ)%R)
                CB_LRB(II,JJ)%ISLR = .FALSE.
              ENDDO
              ENDDO
              CALL DMUMPS_BLR_SAVE_CB_LRB(IW(IOLDPS+XXF),CB_LRB)
              NFS4FATHER = -9999
              IF ( (KEEP(219).NE.0).AND.(KEEP(50).EQ.2) ) THEN
               CALL DMUMPS_BLR_RETRIEVE_NFS4FATHER ( IW(IOLDPS+XXF),
     &             NFS4FATHER )
                NFS4FATHER = max(NFS4FATHER,0) + NELIM
              ENDIF
              ALLOCATE(M_ARRAY(max(1,NFS4FATHER)), stat=allocok)
              IF ( allocok .GT. 0 ) THEN
                IF (LP > 0 ) WRITE(LP,*) MYID,
     &          ": ALLOCATION FAILURE FOR M_ARRAY ",
     &          "DMUMPS_PROCESS_SYM_BLOCFACTO"
                IFLAG = -13
                IERROR = max(1,NFS4FATHER)
              ENDIF
              IF (.NOT.LR_ACTIV_POSTPROC) THEN
                BEGS_BLR_COL(1+NPARTSASS_COL) = 
     &               BEGS_BLR_COL(1+NPARTSASS_COL) - NELIM
              ENDIF
              CALL MAX_CLUSTER(
     &         BEGS_BLR_COL(max(NPARTSASS_MASTER,1)+1:NB_BLR_COL+1),
     &         NB_BLR_COL-max(NPARTSASS_MASTER,1),MAXI_CLUSTER_COL
     &        )
              MAXI_CLUSTER=max(MAXI_CLUSTER_LS, MAXI_CLUSTER_COL)
              IF (allocated(BLOCKLR)) DEALLOCATE(BLOCKLR)
              IF (allocated(RWORK))  DEALLOCATE(RWORK)
              IF (allocated(TAU)) DEALLOCATE(TAU)
              IF (allocated(JPVT)) DEALLOCATE(JPVT)
              IF (allocated(WORK)) DEALLOCATE(WORK)
              LWORK = MAXI_CLUSTER*(MAXI_CLUSTER+1)
              ALLOCATE(BLOCKLR(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &                  RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &                  TAU(MAXI_CLUSTER*OMP_NUM),
     &                  JPVT(MAXI_CLUSTER*OMP_NUM), 
     &                  WORK(LWORK*OMP_NUM),
     &                  stat=allocok)
              IF (allocok > 0 ) THEN
                IFLAG  = -13
                IERROR = OMP_NUM*(LWORK + MAXI_CLUSTER*(MAXI_CLUSTER+4))
                GOTO 700
              ENDIF
              NBROWSinF    = 0
              NVSCHUR_K253 =  0
              IF ( (KEEP(219).NE.0).AND.(KEEP(50).EQ.2).AND.
     &             NFS4FATHER.GT.0  ) THEN
               CALL DMUMPS_COMPUTE_NBROWSinF (
     &                N, INODE, FPERE, KEEP, 
     &                IOLDPS, HS, 
     &                IW, LIW, 
     &                NROW1, NCOL1, NPIV1,
     &                NELIM, NFS4FATHER,
     &                NBROWSinF
     &                )
               IF ((KEEP(114).EQ.1) .AND. (KEEP(116).GT.0) ) THEN
                  NSLAVES_L = IW(PTRIST(STEP(INODE)) + 5 + KEEP(IXSZ))
                  IROW_L    = PTRIST(STEP(INODE)) + 6 + NSLAVES_L + 
     &                        KEEP(IXSZ)
                  CALL DMUMPS_GET_SIZE_SCHUR_IN_FRONT ( 
     &                 N, 
     &                 NROW1,
     &                 KEEP(116), 
     &                 IW(IROW_L),
     &                 PERM, NVSCHUR_K253 )
               ELSE IF (KEEP(253).NE.0) THEN
                  NSLAVES_L = IW(PTRIST(STEP(INODE)) + 5 + KEEP(IXSZ))
                  IROW_L    = PTRIST(STEP(INODE)) + 6 + NSLAVES_L + 
     &                        KEEP(IXSZ)
                  CALL DMUMPS_GET_SIZE_SCHUR_IN_FRONT ( 
     &                 N, 
     &                 NROW1,
     &                 0,  
     &                 IW(IROW_L),
     &                 PERM, NVSCHUR_K253 )
               ENDIF
              ENDIF
              IF (IFLAG.LT.0) GOTO 700
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL
#endif
              CALL DMUMPS_COMPRESS_CB_I(
     &        A_PTR(POSELT), LA_PTR, 1_8, NCOL1,
     &        BEGS_BLR_LS(1), size(BEGS_BLR_LS),
     &        BEGS_BLR_COL(1), size(BEGS_BLR_COL),
     &        NB_BLR_LS, NB_BLR_COL-NPARTSASS_COL,
     &        NPARTSASS_COL, 
     &        NROW1, NCOL1-NPIV1, INODE,
     &        IW(IOLDPS+XXF), 1, 2, IFLAG, IERROR,
     &        DKEEP(12), KEEP(466), KEEP(484), KEEP(489),
     &        CB_LRB(1,1),
     &        WORK, TAU, JPVT, LWORK, RWORK, BLOCKLR,
     &        MAXI_CLUSTER, KEEP8, OMP_NUM,
     &        NFS4FATHER, NPIV1, NVSCHUR_K253, KEEP(1), 
     &        M_ARRAY
     &        , NELIM, NBROWSinF
     &        )
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif
              IF (IFLAG.LT.0) GOTO 650
              IF ( (KEEP(219).NE.0).AND.(KEEP(50).EQ.2).AND.
     &             NFS4FATHER.GT.0  ) THEN
                 INFO_TMP(1) = IFLAG
                 INFO_TMP(2) = IERROR
                 CALL DMUMPS_BLR_SAVE_M_ARRAY( IW(IOLDPS+XXF),
     &            M_ARRAY, INFO_TMP)
                 IFLAG  = INFO_TMP(1) 
                 IERROR = INFO_TMP(2) 
              ENDIF
              DEALLOCATE(M_ARRAY)
 650          CONTINUE         
            ENDIF 
            IF (IFLAG.LT.0) GOTO 700
          ENDIF  
        CALL DMUMPS_END_FACTO_SLAVE( COMM_LOAD, ASS_IRECV, 
     &    N, INODE, FPERE, 
     &    root, roota,
     &    MYID, COMM,
     &    
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     &    PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP, PIMASTER,
     &    PAMASTER,
     &    NSTK_S, COMP, IFLAG, IERROR, PERM,
     &    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, DAD, PTRARW, PTRAIW,
     &    PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &    INTARR,DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND,FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &               , LRGROUPS
     &      )
      ENDIF 
      GOTO 550
 700  CONTINUE
      CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
 550  CONTINUE
      IF (LR_ACTIVATED) THEN
        IF (allocated(RWORK))  DEALLOCATE(RWORK)
        IF (allocated(WORK)) DEALLOCATE(WORK)
        IF (allocated(TAU)) DEALLOCATE(TAU)
        IF (allocated(JPVT)) DEALLOCATE(JPVT)
        IF (allocated(BLOCKLR)) DEALLOCATE(BLOCKLR)
        IF (.NOT.KEEP_BEGS_BLR_LS) THEN
          IF (associated(BEGS_BLR_LS)) DEALLOCATE(BEGS_BLR_LS)
        ENDIF
        IF (.NOT.KEEP_BLR_LS) THEN
          IF (associated(BLR_LS)) THEN
            CALL DEALLOC_BLR_PANEL(BLR_LS, NB_BLR_LS, KEEP8, KEEP(34))
            DEALLOCATE(BLR_LS)
           ENDIF
        ENDIF
        IF (associated(BEGS_BLR_LM)) DEALLOCATE(BEGS_BLR_LM)
        IF (.NOT.KEEP_BEGS_BLR_COL) THEN
              IF (associated(BEGS_BLR_COL)) THEN 
                DEALLOCATE( BEGS_BLR_COL)
              ENDIF
        ENDIF
      ENDIF
      KEEP(174)=KEEP(174)-1
      RETURN
      END SUBROUTINE DMUMPS_PROCESS_SYM_BLOCFACTO
