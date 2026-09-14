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
      RECURSIVE SUBROUTINE DMUMPS_PROCESS_BLOCFACTO(
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
     &    PTRARW, PTRAIW, PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &    INTARR, DBLARR,
     &    ICNTL, KEEP,KEEP8, DKEEP, 
     &    IPOOL, LPOOL, LEAF, ND, FRERE_STEPS,
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
      USE DMUMPS_LR_CORE
      USE DMUMPS_LR_TYPE
      USE MUMPS_LR_STATS
      USE DMUMPS_FAC_LR
      USE DMUMPS_ANA_LR, ONLY : GET_CUT
      USE DMUMPS_LR_DATA_M
      USE MUMPS_INTR_TYPES, ONLY : MUMPS_ROOT_STRUC
      USE DMUMPS_INTR_TYPES, ONLY : DMUMPS_ROOT_STRUC
      USE DMUMPS_DYNAMIC_MEMORY_M, ONLY : DMUMPS_DM_SET_DYNPTR
!$    USE OMP_LIB
      IMPLICIT NONE
      INCLUDE 'mumps_headers.h'
      TYPE (MUMPS_ROOT_STRUC) :: root
      TYPE (DMUMPS_ROOT_STRUC) :: roota
      INTEGER ICNTL( 60 ), KEEP( 500 )
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION DKEEP(230)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER N, SLAVEF, IWPOS, IWPOSCB, LIW
      INTEGER(8) :: IPTRLU, LRLU, LRLUS, LA
      INTEGER(8) :: POSFAC
      INTEGER COMP
      INTEGER IFLAG, IERROR, NBFIN, MSGSOU
      INTEGER PROCNODE_STEPS(KEEP(28)), PTRIST(KEEP(28)),
     &        NSTK_S(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER PERM(N), STEP(N), 
     & PIMASTER(KEEP(28))
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER, intent(in) :: LRGROUPS(KEEP(280))
      INTEGER COMM, MYID
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
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
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER INODE, POSITION, NPIV, IERR, LP
      INTEGER NCOL_RECV, JBEG_BLOCK, NCOL_GEMM, SHIFT_LPOS, SHIFT_UPOS
      INTEGER SHIFT_BEGS_BLR_U
      INTEGER :: IFLAG_OOC
      INTEGER(8) :: POSBLOCFACTO
      INTEGER :: LD_BLOCFACTO 
      INTEGER(8) :: LA_BLOCFACTO 
      INTEGER(8) :: LA_PTR 
      INTEGER(8) :: POSELT
      DOUBLE PRECISION, DIMENSION(:), POINTER :: A_PTR
      INTEGER IOLDPS, LCONT1, NASS1, NROW1, NCOL1, NPIV1, NPIV_TOT
      INTEGER NSLAV1, HS, ISW, ISHIFT
      INTEGER (8) :: LPOS, UPOS, LPOS2, IPOS, KPOS
      INTEGER ICT11
      INTEGER I, IPIV, FPERE
      LOGICAL LASTPANEL, KEEP_BEGS_BLR_L, KEEP_BEGS_BLR_COL
      LOGICAL LASTBL_INPANEL, LASTBL_INLASTPANEL
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      DOUBLE PRECISION ONE,ALPHA
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      INTEGER LIWFAC, STRAT, NextPivDummy
      TYPE(IO_BLOCK) :: MonBloc
      LOGICAL LAST_CALL
      INTEGER :: INFO_TMP(2)
      INTEGER :: IDUMMY(1)
      INTEGER :: NELIM, NPARTSASS_MASTER, NPARTSASS_MASTER_AUX,
     &           IPANEL, IP,
     &           CURRENT_BLR, 
     &           NB_BLR_L, NB_BLR_U, NB_BLR_COL
      INTEGER :: SIZE_BEGS_BLR_COL
      TYPE (LRB_TYPE), POINTER, DIMENSION(:,:) :: CB_LRB
      TYPE (LRB_TYPE), DIMENSION(:), POINTER :: BLR_U, BLR_L
      LOGICAL :: LR_ACTIV_BlocFACTO, COMPRESS_CB, COMPRESS_PANEL, 
     &           LR_ACTIVATED, FORCE_BFR
      LOGICAL :: LR_ACTIV_POSTPROC
      LOGICAL OOCWRITE_COMPATIBLE_WITH_BLR
      INTEGER :: LR_ACTIVATED_INT
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_L, BEGS_BLR_U,
     & BEGS_BLR_COL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: WORK, TAU
      INTEGER, ALLOCATABLE, DIMENSION(:) :: JPVT
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: RWORK
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: BLOCK
      INTEGER :: OMP_NUM
      INTEGER NPARTSASS, NPARTSCB, MAXI_CLUSTER, LWORK,
     &        MAXI_CLUSTER_L, MAXI_CLUSTER_U, MAXI_CLUSTER_COL
      INTEGER  :: allocok
#if defined(USE_GPU)
      INTEGER, EXTERNAL :: DMUMPS_GEMM_ON_GPU_C, DMUMPS_TRSM_ON_GPU_C
      INTEGER :: run_gemm_on_gpu, run_trsm_on_gpu, async_trsm
#endif
      KEEP_BEGS_BLR_COL = .FALSE.
      KEEP_BEGS_BLR_L   = .FALSE.
      nullify(BEGS_BLR_L)
      NB_BLR_U         = -7654321
      SHIFT_BEGS_BLR_U = 0
      NULLIFY(BEGS_BLR_U)
      NULLIFY(BEGS_BLR_COL)
      MAXI_CLUSTER = 0
      CURRENT_BLR  = 1 
      FPERE    = -1
      POSITION = 0
      OMP_NUM = 1
#if ! defined(BLR_NOOPENMP)
!$           OMP_NUM = OMP_GET_MAX_THREADS()
#endif
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NPIV, 1,
     &                 MPI_INTEGER, COMM, IERR )
      LASTPANEL = (NPIV.LE.0)
      IF (LASTPANEL) THEN 
         NPIV = -NPIV
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, FPERE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      ENDIF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NCOL_RECV, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, JBEG_BLOCK, 1,
     &                 MPI_INTEGER, COMM, IERR )
      IF (JBEG_BLOCK.EQ.1) THEN
        NCOL_GEMM = NCOL_RECV - NPIV
        SHIFT_LPOS = NPIV
        SHIFT_UPOS = NPIV
      ELSE
        NCOL_GEMM = NCOL_RECV
        SHIFT_LPOS = JBEG_BLOCK-1
        SHIFT_UPOS = 0
      ENDIF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NELIM, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, 
     &                 NPARTSASS_MASTER , 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, IPANEL,
     &                 1, MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, LR_ACTIVATED_INT,
     &                 1, MPI_INTEGER, COMM, IERR )
      LR_ACTIV_BlocFACTO    = (LR_ACTIVATED_INT.EQ.1)
      IF ( LR_ACTIV_BlocFACTO ) THEN
       IF (JBEG_BLOCK.NE.1) THEN
          LA_BLOCFACTO = 0_8
        ELSE
          LA_BLOCFACTO = int(NPIV,8) * int(NPIV+NELIM,8)
        ENDIF
      ELSE
        LA_BLOCFACTO = int(NPIV,8) * int(NCOL_RECV,8)
      ENDIF
      CALL DMUMPS_GET_SIZE_NEEDED(
     &      NPIV, LA_BLOCFACTO, .FALSE.,
     &      KEEP(1), KEEP8(1),
     &      N, IW, LIW, A, LA,
     &      LRLU, IPTRLU,
     &      IWPOS, IWPOSCB, PTRIST, PTRAST,
     &      STEP, PIMASTER, PAMASTER, LRLUS,
     &      KEEP(IXSZ),COMP,DKEEP(97),MYID,SLAVEF, PROCNODE_STEPS, 
     &      DAD, IFLAG, IERROR)
      IF (IFLAG.LT.0) GOTO 700
      LRLU  = LRLU - LA_BLOCFACTO
      LRLUS = LRLUS - LA_BLOCFACTO
      KEEP8(67) = min(LRLUS, KEEP8(67))
      KEEP8(69) = KEEP8(69) + LA_BLOCFACTO
      KEEP8(68) = max(KEEP8(69), KEEP8(68))
      POSBLOCFACTO = POSFAC
      POSFAC = POSFAC + LA_BLOCFACTO
      CALL MUMPS_LOAD_MEM_UPDATE(.FALSE., .FALSE.,
     &               LA-LRLUS,0_8,LA_BLOCFACTO,KEEP,KEEP8,LRLUS)
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, 
     &               SIZE_BEGS_BLR_COL, 1,
     &               MPI_INTEGER, COMM, IERR )
      IF (SIZE_BEGS_BLR_COL.GT.0) THEN
          ALLOCATE(BEGS_BLR_COL(SIZE_BEGS_BLR_COL),
     &        stat=allocok)
          IF (allocok > 0 ) THEN
             IFLAG = -13
             IERROR = SIZE_BEGS_BLR_COL
             LP = ICNTL(1)
             IF (ICNTL(4) .LE. 0) LP=-1
             IF (LP > 0) WRITE(LP,*) MYID,
     &            ': ERROR allocation during DMUMPS_PROCESS_BLOCFACTO'
             GOTO 700
          ENDIF
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, 
     &                 BEGS_BLR_COL, SIZE_BEGS_BLR_COL,
     &                 MPI_INTEGER, COMM, IERR )
      ENDIF
      IF ((NPIV .EQ. 0) 
     &     ) THEN
        IPIV=1
      ELSE
        IPIV = IWPOS
        IWPOS = IWPOS + NPIV
        IF (NPIV .GT. 0.AND.JBEG_BLOCK.EQ.1) THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 IW( IPIV ), NPIV,
     &                 MPI_INTEGER, COMM, IERR )
        ENDIF
        IF ( LR_ACTIV_BlocFACTO .AND. JBEG_BLOCK.EQ.1) THEN
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 A(POSBLOCFACTO), NPIV*(NPIV+NELIM),
     &                 MPI_DOUBLE_PRECISION,
     &                 COMM, IERR )
            LD_BLOCFACTO = NPIV+NELIM
        ENDIF
        IF ( LR_ACTIV_BlocFACTO ) THEN
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 NB_BLR_U, 1, MPI_INTEGER,
     &                 COMM, IERR )
            ALLOCATE(BLR_U(max(NB_BLR_U,1)), stat=allocok)
            IF (allocok > 0 ) THEN
               IFLAG = -13
               IERROR = max(NB_BLR_U,1)
               LP = ICNTL(1)
               IF (ICNTL(4) .LE. 0) LP=-1
               IF (LP > 0) WRITE(LP,*) MYID,
     &              ': ERROR allocation during DMUMPS_PROCESS_BLOCFACTO'
               GOTO 700
            ENDIF
            ALLOCATE(BEGS_BLR_U(NB_BLR_U+2), stat=allocok)
            IF (allocok > 0 ) THEN
               IFLAG = -13
               IERROR = NB_BLR_U+2
               LP = ICNTL(1)
               IF (ICNTL(4) .LE. 0) LP=-1
               IF (LP > 0) WRITE(LP,*) MYID,
     &              ': ERROR allocation during DMUMPS_PROCESS_BLOCFACTO'
               GOTO 700
            ENDIF
            CALL DMUMPS_MPI_UNPACK_LR_PARTIAL 
     &                        (BUFR, LBUFR, LBUFR_BYTES, 
     &                             POSITION, NPIV, NELIM, 'H',
     &                             BLR_U(1), NB_BLR_U, 
     &                             JBEG_BLOCK,
     &                             BEGS_BLR_U(1),
     &                             KEEP8, COMM, IERR, IFLAG, IERROR)
            IF (JBEG_BLOCK.NE.1) SHIFT_BEGS_BLR_U = 1
            IF (IFLAG.LT.0) GOTO 700
        ELSE
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 A(POSBLOCFACTO), NPIV*NCOL_RECV,
     &                 MPI_DOUBLE_PRECISION,
     &                 COMM, IERR )
          LD_BLOCFACTO = NCOL_RECV
        ENDIF
      ENDIF 
      IF (PTRIST(STEP( INODE )) .EQ. 0) THEN
          CALL DMUMPS_TREAT_DESCBAND( INODE, COMM_LOAD,
     &    ASS_IRECV,
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
     &    )
          IF ( IFLAG .LT. 0 ) GOTO 550
      ENDIF
      IF ( IW( PTRIST(STEP(INODE)) + 3 +KEEP(IXSZ)) .EQ. 0 ) THEN
       DO WHILE ( IW(PTRIST(STEP(INODE)) + XXNBPR) .NE. 0)
        BLOCKING = .TRUE.
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
      LCONT1 = IW( IOLDPS + KEEP(IXSZ))
      NASS1  = IW( IOLDPS + 1 + KEEP(IXSZ))
      LR_ACTIVATED   = (IW(IOLDPS+XXLR).GT.0)
      FORCE_BFR      = MUMPS_LR_FORCE_BFR(IW(IOLDPS+XXLR))
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
     &          ( (.NOT.COMPRESS_PANEL)
     &            .OR. (KEEP(486).NE.2) 
     &          )
      LASTBL_INPANEL     = JBEG_BLOCK+NCOL_RECV.GT.LCONT1
      LASTBL_INLASTPANEL = LASTPANEL .AND. LASTBL_INPANEL
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
      NROW1  = IW( IOLDPS + 2 +KEEP(IXSZ))
      NPIV1  = IW( IOLDPS + 3 +KEEP(IXSZ))
      NPIV_TOT = NPIV1 + NPIV
      COMPRESS_CB = MUMPS_LR_COMPRESS_CB(IW(IOLDPS+XXLR))
      IF (LR_ACTIVATED) THEN 
        IF (COMPRESS_CB.AND.NPIV.EQ.0) THEN
           COMPRESS_CB = .FALSE.
           IW(IOLDPS+XXLR) = IW(IOLDPS+XXLR) -1
        ENDIF
      ENDIF
      NSLAV1 = IW( IOLDPS + 5 + KEEP(IXSZ))
      HS     = 6 + NSLAV1 + KEEP(IXSZ)
      NCOL1  = LCONT1 + NPIV1
#if defined(USE_GPU)
      run_trsm_on_gpu=0
      run_gemm_on_gpu=0
#endif
      IF (NPIV.GT.0) THEN
#if defined(USE_GPU)
        IF ( .NOT. LR_ACTIV_BlocFACTO .OR. KEEP(475) .EQ. 0) THEN
          run_trsm_on_gpu=DMUMPS_TRSM_ON_GPU_C(
     &                  KEEP(417), NPIV, NROW1, 'L')
        ENDIF
        IF ( .NOT. LR_ACTIV_BlocFACTO ) THEN
          IF (NCOL_GEMM.GT.0) THEN
            run_gemm_on_gpu=DMUMPS_GEMM_ON_GPU_C(
     &                  KEEP(417), NCOL_GEMM, NROW1, NPIV )
          ENDIF
        ENDIF
        IF ( KEEP(422).EQ.2 ) THEN
          IF (IW(IOLDPS+XXG) .EQ. PinningOnTheWay) THEN
            CALL DMUMPS_GPU_WAIT_ENDOFPINNING_C()
            IW(IOLDPS+XXG) = MemPinned
          ENDIF
        ELSE IF (KEEP(422) .EQ.1 ) THEN
          IF ( IW(IOLDPS+XXG) .NE. MemPinned ) THEN
            IF ( run_trsm_on_gpu.EQ.1 .OR. run_gemm_on_gpu.EQ.1 ) THEN
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
            IF (IW(IPIV+I-1).EQ.I) CYCLE
            ISW = IW(ICT11+I)
            IW(ICT11+I) = IW(ICT11+IW(IPIV+I-1))
            IW(ICT11+IW(IPIV+I-1)) = ISW
            IPOS = POSELT + int(NPIV1 + I - 1,8)
            KPOS = POSELT + int(NPIV1 + IW(IPIV+I-1) - 1,8)
            CALL dswap(NROW1, A_PTR(IPOS), NCOL1, A_PTR(KPOS), NCOL1)
          ENDDO
        ENDIF
        LPOS2 = POSELT + int(NPIV1,8)
        LPOS  = LPOS2  + int(SHIFT_LPOS,8)
#if defined(USE_GPU)
        IF ( KEEP(422).NE.0 .AND. KEEP(423) .EQ. 1 .AND.
     &       ( run_gemm_on_gpu.EQ.1 ) ) THEN
             CALL DMUMPS_GPU_PIN_MEMORY_C(A(POSBLOCFACTO),LA_BLOCFACTO)
        ENDIF
#endif
        IF ( (JBEG_BLOCK.EQ.1) .AND.
     &       ((.NOT. LR_ACTIV_BlocFACTO).OR.KEEP(475).EQ.0)
     &     ) THEN
#if defined(USE_GPU)
          IF (run_trsm_on_gpu.EQ.1.and.run_gemm_on_gpu.EQ.1) THEN
            IF (KEEP(427).LT.0) THEN
              async_trsm=1 
            ELSEIF (KEEP(427) .EQ.0) THEN
              async_trsm=0 
            ELSE
              async_trsm=2 
            ENDIF
            CALL DMUMPS_TRSM_GPU_C(async_trsm,
     &             'L','L','N','N',NPIV, NROW1, ONE, 
     &             A(POSBLOCFACTO), LD_BLOCFACTO,
     &             A_PTR(LPOS2), NCOL1)
          ELSE
#endif
            CALL dtrsm('L','L','N','N', NPIV, NROW1, ONE,
     &             A(POSBLOCFACTO), LD_BLOCFACTO,
     &             A_PTR(LPOS2), NCOL1)
#if defined(USE_GPU)
          ENDIF
#endif
        ENDIF
      ENDIF 
      IF (LR_ACTIVATED) THEN
        IF (NPIV.NE.0) THEN
        IF ( (NPIV1.EQ.0).AND.(JBEG_BLOCK.EQ.1) 
     &     ) THEN
          IOLDPS = PTRIST(STEP(INODE))
          CALL GET_CUT(IW(IOLDPS+HS:IOLDPS+HS+NROW1-1), 0,
     &                    NROW1, LRGROUPS, NPARTSCB, 
     &                    NPARTSASS, BEGS_BLR_L)
          CALL REGROUPING2(BEGS_BLR_L, NPARTSASS, NASS1, NPARTSCB,
     &                        NROW1-0, KEEP(488), .TRUE., KEEP(472),
     &                        NCOL1, KEEP(1))
          NB_BLR_L =  NPARTSCB
          IF (IFLAG.LT.0) GOTO 700
          CALL DMUMPS_BLR_SAVE_INIT(IW(IOLDPS+XXF), 
     &           .FALSE.,       
     &           .TRUE.,        
     &           .TRUE.,        
     &           NPARTSASS_MASTER, 
     &           BEGS_BLR_L, 
     &           BEGS_BLR_COL, 
     &           huge(NPARTSASS_MASTER),
     &           INFO_TMP)
          IF (associated(BEGS_BLR_COL)) DEALLOCATE(BEGS_BLR_COL)
          IF (IFLAG.LT.0) GOTO 700
        ELSE           
         IF (LR_ACTIV_BlocFACTO) THEN
          CALL DMUMPS_BLR_RETRIEVE_BEGS_BLR_L (IW(IOLDPS+XXF), 
     &                  BEGS_BLR_L)
          KEEP_BEGS_BLR_L = .TRUE.  
          NB_BLR_L  = size(BEGS_BLR_L) - 2
          NPARTSASS = 1
          NPARTSCB  = NB_BLR_L
         ENDIF 
        ENDIF 
        ENDIF 
      ENDIF 
        IF ( (NPIV .GT. 0)
     &   ) THEN
        IF (LR_ACTIV_BlocFACTO) THEN
        IF ( (JBEG_BLOCK.EQ.1)
     &     ) THEN
        call MAX_CLUSTER(BEGS_BLR_L,NB_BLR_L+1,MAXI_CLUSTER_L)
        call MAX_CLUSTER(BEGS_BLR_U(1+SHIFT_BEGS_BLR_U:NB_BLR_U+2),
     &                   NB_BLR_U+1-SHIFT_BEGS_BLR_U,
     &                   MAXI_CLUSTER_U)
        IF (SHIFT_BEGS_BLR_U.EQ.1) 
     &      MAXI_CLUSTER_U = max(MAXI_CLUSTER_U,NPIV+NELIM)
        IF (LASTBL_INLASTPANEL.AND.COMPRESS_CB) THEN
          MAXI_CLUSTER=max(MAXI_CLUSTER_U+NELIM,MAXI_CLUSTER_L)
        ELSE
          MAXI_CLUSTER=max(MAXI_CLUSTER_U,MAXI_CLUSTER_L)
        ENDIF
        LWORK = MAXI_CLUSTER*(MAXI_CLUSTER+1)
        ALLOCATE(BLOCK(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &       RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &       TAU(MAXI_CLUSTER*OMP_NUM),
     &       JPVT(MAXI_CLUSTER*OMP_NUM), 
     &       WORK(LWORK*OMP_NUM), stat=allocok)
        IF (allocok > 0 ) THEN
           IFLAG = -13
           IERROR = MAXI_CLUSTER*OMP_NUM*MAXI_CLUSTER
     &          + 2*MAXI_CLUSTER*OMP_NUM
     &          + MAXI_CLUSTER*OMP_NUM
     &          + MAXI_CLUSTER*OMP_NUM
     &          + LWORK*OMP_NUM
           LP = ICNTL(1)
           IF (ICNTL(4) .LE. 0) LP=-1
           IF (LP > 0) WRITE(LP,*) MYID,
     &          ': ERROR allocation during DMUMPS_PROCESS_BLOCFACTO'
           GOTO 700
        ENDIF
        CURRENT_BLR=1 
        ALLOCATE(BLR_L(NB_BLR_L), stat=allocok)
        IF (allocok > 0 ) THEN
              IFLAG = -13
              IERROR = NB_BLR_L
              LP = ICNTL(1)
              IF (ICNTL(4) .LE. 0) LP=-1
              IF (LP > 0) WRITE(LP,*) MYID,
     &             ': ERROR allocation during DMUMPS_PROCESS_BLOCFACTO'
              GOTO 700
           ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL
#endif
        CALL DMUMPS_COMPRESS_PANEL_I_NOOPT
     &        (A_PTR(POSELT), LA_PTR, 1_8,
     &        IFLAG, IERROR, NCOL1,
     &        BEGS_BLR_L(1), size(BEGS_BLR_L), NB_BLR_L+1,
     &        DKEEP(8), KEEP(466), 0, 
     &        KEEP(473), BLR_L(1), 
     &        CURRENT_BLR, 'V', WORK, TAU, JPVT, LWORK, RWORK,
     &        BLOCK, MAXI_CLUSTER, NELIM, 
     &        .TRUE.,  
     &        NPIV, NPIV1,
     &        2, KEEP(483), KEEP8,
     &        (KEEP(450).GT.0), KEEP(453), KEEP(454).EQ.2,
     &        OMP_NUM )
#if ! defined(BLR_NOOPENMP)
!$OMP MASTER
#endif
          CALL DMUMPS_BLR_SAVE_PANEL_LORU (
     &         IW(IOLDPS+XXF),
     &         0, 
     &         IPANEL, BLR_L, 0) 
#if ! defined(BLR_NOOPENMP)
!$OMP END MASTER
!$OMP BARRIER
#endif
          IF (IFLAG.LT.0) GOTO 300
          IF ((KEEP(475).GE.1).AND.(JBEG_BLOCK.EQ.1)) THEN
            CALL DMUMPS_BLR_PANEL_LRTRSM(A, LA, POSBLOCFACTO, 
     &              LD_BLOCFACTO, -6666, 
     &              NB_BLR_L+1,  
     &              BLR_L, CURRENT_BLR, CURRENT_BLR+1, NB_BLR_L+1, 
     &              2, 0, 0,  
     &              .TRUE.) 
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif          
            IF (KEEP(486).NE.2) THEN
              CALL DMUMPS_DECOMPRESS_PANEL_I_NOOPT(
     &        A_PTR(POSELT), LA_PTR, 1_8,
     &        NCOL1, NCOL1,
     &        .TRUE.,    
     &        NPIV1+1,   
     &        1,         
     &        NB_BLR_L+1, BLR_L(1), CURRENT_BLR, 'V', 1)
            ENDIF
          ENDIF
 300      CONTINUE
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif          
          IF (IFLAG.LT.0) GOTO 700
      ELSE 
          CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &           IW(IOLDPS+XXF),
     &           0, 
     &           IPANEL, BLR_L)
          CURRENT_BLR=1 
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
      IF ( (NPIV .GT. 0)
     &   ) THEN
        IF (LR_ACTIV_BlocFACTO) THEN
          IF (NELIM.GT.0.AND.JBEG_BLOCK.EQ.1) THEN
            UPOS = 1_8+int(SHIFT_UPOS,8)
            CALL DMUMPS_BLR_UPD_NELIM_VAR_L_I(
     &        A(POSBLOCFACTO), LA_BLOCFACTO, UPOS,
     &        A_PTR(POSELT), LA_PTR, LPOS-POSELT+1_8,
     &        IFLAG, IERROR, LD_BLOCFACTO, NCOL1,
     &        BEGS_BLR_L(1), size(BEGS_BLR_L),
     &        CURRENT_BLR, BLR_L(1), NB_BLR_L+1, 
     &        CURRENT_BLR+1, NELIM, 'N')
          ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL
#endif
         CALL DMUMPS_BLR_UPDATE_TRAILING_I(
     &        A_PTR(POSELT), LA_PTR, 1_8, 
     &        IFLAG, IERROR, NCOL1,
     &        BEGS_BLR_L(1), size(BEGS_BLR_L),
     &        BEGS_BLR_U(1), size(BEGS_BLR_U), CURRENT_BLR, 
     &        BLR_L(1), NB_BLR_L+1, 
     &        BLR_U(1), NB_BLR_U+1,  
     &        0,      
     &        .TRUE., 
     &        NPIV1,  
     &        2, 0, 
     &        KEEP(481), DKEEP(11), KEEP(466), KEEP(477) 
     &        )
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif          
          IF (IFLAG.LT.0) GOTO 700
        ELSE 
#if defined(USE_GPU)
          IF (run_gemm_on_gpu .EQ.1) THEN
            UPOS = POSBLOCFACTO+int(SHIFT_UPOS,8)
            CALL DMUMPS_GEMM_GPU_C(0,'N','N', NCOL_GEMM, NROW1, NPIV,
     &             ALPHA,A(UPOS),NCOL_RECV,
     &             A_PTR(LPOS2), NCOL1, ONE, A_PTR(LPOS), NCOL1)
          ELSE
#endif
            UPOS = POSBLOCFACTO+int(SHIFT_UPOS,8)
            CALL dgemm('N','N', NCOL_GEMM, NROW1, NPIV,
     &             ALPHA,A(UPOS), NCOL_RECV,
     &             A_PTR(LPOS2), NCOL1, ONE, A_PTR(LPOS), NCOL1)
#if defined(USE_GPU)
          ENDIF
#endif
#if defined(USE_GPU)
          IF (KEEP(422).NE.0 .AND. KEEP(423).EQ.1 .AND.
     &        run_gemm_on_gpu.EQ.1) THEN
              CALL DMUMPS_GPU_UNPIN_MEMORY_C(
     &               A(POSBLOCFACTO),LA_BLOCFACTO)
          ENDIF
#endif
        ENDIF 
      ENDIF 
      IF (LASTBL_INPANEL) THEN
        IW(IOLDPS + KEEP(IXSZ))    = IW(IOLDPS + KEEP(IXSZ)) - NPIV
        IW(IOLDPS + 3 + KEEP(IXSZ))= IW(IOLDPS + 3 + KEEP(IXSZ)) + NPIV
        IF (LASTPANEL) THEN
          IW(IOLDPS+1+KEEP(IXSZ) ) = IW(IOLDPS+3+KEEP(IXSZ) )
        ENDIF
      ENDIF
      IF ( .not. LASTBL_INLASTPANEL .AND. 
     &   (IW(IOLDPS+1+KEEP(IXSZ)) .EQ. IW(IOLDPS + 3+KEEP(IXSZ))) ) THEN
        write(*,*) 'Internal ERROR 1 **** IN BLACFACTO '
        CALL MUMPS_ABORT()
      ENDIF
      IF (LR_ACTIV_BlocFACTO) THEN
        IF ((NPIV.GT.0)
     &     ) THEN
          CALL DEALLOC_BLR_PANEL( BLR_U, NB_BLR_U, KEEP8, KEEP(34))
          DEALLOCATE(BLR_U)
           IF (KEEP(486).NE.3) THEN
            IF (KEEP(450).GT.0) THEN
              IF (LASTBL_INPANEL) THEN
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL
#endif
                CALL DMUMPS_LRMIXED_CONVERT_PANEL(BLR_L, NB_BLR_L,
     &                                            0, KEEP(34), KEEP8,
     &                                            IFLAG, IERROR)
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif
              ENDIF
              IF (IFLAG.LT.0) GOTO 700
            ENDIF
            IF (LASTBL_INPANEL) THEN
              CALL UPD_MRY_LU_LRGAIN(BLR_L, NPARTSCB
     &               , KEEP(450) .GT. 0
     &               )
            ENDIF
          ENDIF
      ENDIF 
      ENDIF 
        LRLU  = LRLU + LA_BLOCFACTO
        LRLUS = LRLUS + LA_BLOCFACTO
        KEEP8(69) = KEEP8(69) - LA_BLOCFACTO
        POSFAC = POSFAC - LA_BLOCFACTO
        CALL MUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &             LA-LRLUS,0_8,-LA_BLOCFACTO,KEEP,KEEP8,LRLUS)
        IWPOS = IWPOS - NPIV
      IF (LASTBL_INPANEL) THEN
        FLOP1 = dble( NPIV1*NROW1 ) +
     &        dble(NROW1*NPIV1)*dble(2*NCOL1-NPIV1-1)
     &   -
     &        dble((NPIV1+NPIV)*NROW1 ) -
     &        dble(NROW1*(NPIV1+NPIV))*dble(2*NCOL1-NPIV1-NPIV-1)
        CALL MUMPS_LOAD_UPDATE( 1, .FALSE., FLOP1, KEEP,KEEP8 )
      ENDIF
      IF (LR_ACTIV_BlocFACTO.AND.LASTBL_INPANEL.AND.
     &     (KEEP(486).EQ.3)
     &   ) THEN
          IF (NPIV.NE.0) THEN
           CALL DMUMPS_BLR_FORCE_FREE_PANEL_L(IW(IOLDPS+XXF), IPANEL, 
     &                            KEEP8, KEEP(34)) 
           nullify(BLR_L)
          ENDIF
      ENDIF
      IF (LASTBL_INLASTPANEL) THEN
        IF (LR_ACTIV_POSTPROC) THEN
          CALL DMUMPS_BLR_RETRIEVE_BEGS_BLR_L (IW(IOLDPS+XXF), 
     &                  BEGS_BLR_L)
          KEEP_BEGS_BLR_L = .TRUE.  
          NB_BLR_L  = size(BEGS_BLR_L) - 2
          NPARTSASS = 1
          NPARTSCB  = NB_BLR_L
          call MAX_CLUSTER(BEGS_BLR_L,NB_BLR_L+1,MAXI_CLUSTER_L)
           CALL DMUMPS_BLR_RETRIEVE_BEGS_BLR_C (IW(IOLDPS+XXF), 
     &                  BEGS_BLR_COL, NPARTSASS_MASTER_AUX)
           KEEP_BEGS_BLR_COL = .TRUE.
           NPIV_TOT =  IW(IOLDPS + 3 + KEEP(IXSZ))
           DO IP=1,NPARTSASS_MASTER
               IF ( BEGS_BLR_COL(IP).GT. NPIV_TOT ) THEN
                 BEGS_BLR_COL(IP)   = NPIV_TOT+1
               ENDIF
           ENDDO
           BEGS_BLR_COL(NPARTSASS_MASTER+1) = NPIV_TOT+1
           call MAX_CLUSTER(BEGS_BLR_COL,NPARTSASS_MASTER_AUX,
     &           MAXI_CLUSTER_U)
           MAXI_CLUSTER = max (MAXI_CLUSTER_L, MAXI_CLUSTER_U)
           LWORK = MAXI_CLUSTER*(MAXI_CLUSTER+1)
           ALLOCATE(BLOCK(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &        RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &        TAU(MAXI_CLUSTER*OMP_NUM),
     &        JPVT(MAXI_CLUSTER*OMP_NUM), 
     &        WORK(LWORK*OMP_NUM), stat=allocok)
           IF (allocok > 0 ) THEN
             IFLAG = -13
             IERROR = MAXI_CLUSTER*OMP_NUM*MAXI_CLUSTER
     &         + 2*MAXI_CLUSTER*OMP_NUM
     &         + MAXI_CLUSTER*OMP_NUM
     &         + MAXI_CLUSTER*OMP_NUM
     &         + LWORK*OMP_NUM
             LP = ICNTL(1)
             IF (ICNTL(4) .LE. 0) LP=-1
             IF (LP > 0) WRITE(LP,*) MYID,
     &        ': ERROR allocation during DMUMPS_PROCESS_BLOCFACTO'
             GOTO 700
           ENDIF
           DO IP=1,NPARTSASS_MASTER
              IF (BEGS_BLR_COL(IP+1)-BEGS_BLR_COL(IP).EQ.0) THEN
                NULLIFY(BLR_L)
              ELSE
                ALLOCATE(BLR_L(NB_BLR_L),stat=allocok)
                IF (allocok > 0) THEN
                  IFLAG  = -13
                  IERROR = NB_BLR_L
                  GOTO 700
                ENDIF
                DO I=1, NB_BLR_L
                   BLR_L(I)%M=0
                   BLR_L(I)%N=0
                   BLR_L(I)%K=0
                   BLR_L(I)%ISLR=.FALSE.
                   NULLIFY(BLR_L(I)%Q)
                   NULLIFY(BLR_L(I)%R)
                ENDDO
              ENDIF
              CALL DMUMPS_BLR_SAVE_PANEL_LORU (
     &            IW(IOLDPS+XXF),
     &            0, 
     &            IP, BLR_L, 0) 
           ENDDO
           CURRENT_BLR = 1 
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL
!$OMP& PRIVATE(IP, BLR_L, ISHIFT)
#endif
           DO IP=1,NPARTSASS_MASTER
             IF (BEGS_BLR_COL(IP+1)-BEGS_BLR_COL(IP).EQ.0) CYCLE
             CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &            IW(IOLDPS+XXF), 0,  
     &            IP, BLR_L)
             ISHIFT = BEGS_BLR_COL(IP)-1
             CALL DMUMPS_COMPRESS_PANEL_I_NOOPT
     &             (A_PTR(POSELT), LA_PTR, 1_8,
     &             IFLAG, IERROR, NCOL1,
     &             BEGS_BLR_L(1), size(BEGS_BLR_L), NB_BLR_L+1,
     &             DKEEP(8), KEEP(466), 0, 
     &             KEEP(473), BLR_L(1), 
     &             1, 'V', WORK, TAU, JPVT, LWORK, RWORK,
     &             BLOCK, MAXI_CLUSTER, NELIM, 
     &             .TRUE.,  
     &             BEGS_BLR_COL(IP+1)-BEGS_BLR_COL(IP), ISHIFT,
     &             2, KEEP(483), KEEP8,
     &            (KEEP(450).GT.0), KEEP(453), KEEP(454).EQ.2,
     &            OMP_NUM )
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif
               IF (IFLAG.LT.0) GOTO 440
               IF (
     &               ((KEEP(486).EQ.2)
     &              ) .AND. (KEEP(450).GT.0)) THEN
                   CALL DMUMPS_LRMIXED_CONVERT_PANEL(BLR_L, NB_BLR_L,
     &                                          0, KEEP(34), KEEP8,
     &                                          IFLAG, IERROR)
               ENDIF
           ENDDO
 440       CONTINUE
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif          
           DO IP=1,NPARTSASS_MASTER
             IF (BEGS_BLR_COL(IP+1)-BEGS_BLR_COL(IP).EQ.0) CYCLE
             CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &            IW(IOLDPS+XXF), 0,  
     &            IP, BLR_L)
             CALL UPD_MRY_LU_LRGAIN(BLR_L, NB_BLR_L
     &               , KEEP(450) .GT. 0
     &               )
           ENDDO
        ENDIF
        IF (IFLAG.LT.0) GOTO 700
        IF (KEEP(486).NE.0) THEN
          IF (LR_ACTIV_BlocFACTO) THEN
            CALL STATS_COMPUTE_FLOP_SLAVE_TYPE2(NROW1, NCOL1, NASS1,
     &              KEEP(50), INODE)
          ELSE
            CALL UPD_FLOP_FRFRONT_SLAVE(NROW1, NCOL1, NASS1,
     &              KEEP(50), INODE)
          ENDIF
        ENDIF
        IF (LR_ACTIVATED) THEN
         IF (COMPRESS_CB) THEN
           CALL DMUMPS_BLR_RETRIEVE_BEGS_BLR_C (IW(IOLDPS+XXF), 
     &                  BEGS_BLR_COL, NPARTSASS_MASTER_AUX)
           KEEP_BEGS_BLR_COL = .TRUE.
           IF (.NOT.LR_ACTIV_POSTPROC) THEN
             BEGS_BLR_COL(1+NPARTSASS_MASTER) = 
     &               BEGS_BLR_COL(1+NPARTSASS_MASTER) - NELIM
           ENDIF
           NB_BLR_COL = size(BEGS_BLR_COL) - 1
            call MAX_CLUSTER(BEGS_BLR_L,NB_BLR_L+1,MAXI_CLUSTER_L)
            call MAX_CLUSTER(BEGS_BLR_COL,NB_BLR_COL,MAXI_CLUSTER_COL)
            MAXI_CLUSTER_COL=MAXI_CLUSTER_COL+NELIM
            IF ( (MAXI_CLUSTER.LT.MAXI_CLUSTER_COL).OR.
     &           (MAXI_CLUSTER.LT.MAXI_CLUSTER_L) ) THEN
             MAXI_CLUSTER = max(MAXI_CLUSTER_COL,MAXI_CLUSTER_L) 
             IF (allocated(RWORK))  DEALLOCATE(RWORK)
             IF (allocated(WORK)) DEALLOCATE(WORK)
             IF (allocated(TAU)) DEALLOCATE(TAU)
             IF (allocated(JPVT)) DEALLOCATE(JPVT)
             IF (allocated(BLOCK)) DEALLOCATE(BLOCK)
             LWORK = MAXI_CLUSTER*(MAXI_CLUSTER+1)
             ALLOCATE(BLOCK(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &         RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &         TAU(MAXI_CLUSTER*OMP_NUM),
     &         JPVT(MAXI_CLUSTER*OMP_NUM), 
     &         WORK(LWORK*OMP_NUM), stat=allocok)
             IF (allocok > 0 ) THEN
               IFLAG = -13
               IERROR = MAXI_CLUSTER*OMP_NUM*MAXI_CLUSTER
     &           + 2*MAXI_CLUSTER*OMP_NUM
     &           + MAXI_CLUSTER*OMP_NUM
     &           + MAXI_CLUSTER*OMP_NUM
     &           + LWORK*OMP_NUM
               LP = ICNTL(1)
               IF (ICNTL(4) .LE. 0) LP=-1
               IF (LP > 0) WRITE(LP,*) MYID,
     &          ': ERROR allocation during DMUMPS_PROCESS_BLOCFACTO'
               GOTO 700
             ENDIF
           ENDIF 
           allocate(CB_LRB(NB_BLR_L,NB_BLR_COL-NPARTSASS_MASTER),
     &                 stat=allocok)
           IF (allocok > 0) THEN
             IFLAG  = -13
             IERROR = NB_BLR_L*(NB_BLR_COL-NPARTSASS_MASTER)
             GOTO 700
           ENDIF
           CALL DMUMPS_BLR_SAVE_CB_LRB(IW(IOLDPS+XXF),CB_LRB)
         ENDIF 
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL
#endif
         IF (COMPRESS_CB) THEN
           CALL DMUMPS_COMPRESS_CB_I(
     &      A_PTR(POSELT), LA_PTR, 1_8, NCOL1,
     &      BEGS_BLR_L(1), size(BEGS_BLR_L),
     &      BEGS_BLR_COL(1), size(BEGS_BLR_COL),
     &      NB_BLR_L, NB_BLR_COL-NPARTSASS_MASTER,
     &      NPARTSASS_MASTER, 
     &      NROW1, NCOL1-NPIV1-NPIV, INODE,
     &      IW(IOLDPS+XXF), 0, 2, IFLAG, IERROR,
     &      DKEEP(12), KEEP(466), KEEP(484), KEEP(489),
     &      CB_LRB(1,1),
     &      WORK, TAU, JPVT, LWORK, RWORK, BLOCK,
     &      MAXI_CLUSTER, KEEP8, OMP_NUM,
     &      -9999, -9999, -9999, KEEP(1),
     &      IDUMMY, 0, -9999 )
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif
         ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif
         IF (IFLAG.LT.0) GOTO 700
         IF ( KEEP(251).EQ.2 .AND. KEEP(486).EQ.2 ) THEN
           CALL DMUMPS_BLR_FREE_ALL_PANELS( IW(IOLDPS+XXF),
     &                          0, 
     &                          KEEP8, KEEP(34) )
         ENDIF
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
      IF (COMPRESS_PANEL) THEN
        IF (allocated(RWORK))  DEALLOCATE(RWORK)
        IF (allocated(WORK)) DEALLOCATE(WORK)
        IF (allocated(TAU)) DEALLOCATE(TAU)
        IF (allocated(JPVT)) DEALLOCATE(JPVT)
        IF (allocated(BLOCK)) DEALLOCATE(BLOCK)
        IF (associated(BEGS_BLR_COL)) THEN
           IF (.NOT. KEEP_BEGS_BLR_COL) DEALLOCATE(BEGS_BLR_COL)
        ENDIF
        IF (associated(BEGS_BLR_L)) THEN
            IF (.NOT. KEEP_BEGS_BLR_L) DEALLOCATE(BEGS_BLR_L)
        ENDIF
        IF ((NPIV.GT.0)
     &     ) THEN
          IF (associated(BEGS_BLR_U)) DEALLOCATE(BEGS_BLR_U)
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_PROCESS_BLOCFACTO
      SUBROUTINE DMUMPS_MPI_UNPACK_LR(
     &           BUFR, LBUFR, LBUFR_BYTES, POSITION,
     &                             NPIV, NELIM, DIR, 
     &                             BLR_U, NB_BLOCK_U,
     &                             BEGS_BLR_U, KEEP8,
     &                             COMM, IERR, IFLAG, IERROR)
      USE DMUMPS_LR_CORE, ONLY : LRB_TYPE, ALLOC_LRB
      USE DMUMPS_LR_TYPE
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LBUFR
      INTEGER, INTENT(IN) :: LBUFR_BYTES
      INTEGER, INTENT(IN) :: BUFR(LBUFR)
      INTEGER, INTENT(INOUT) :: POSITION
      INTEGER, INTENT(IN)    :: NB_BLOCK_U, NELIM, NPIV
      CHARACTER(len=1) :: DIR
      INTEGER, INTENT(IN) :: COMM
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
      INTEGER, INTENT(OUT) :: IERR
      TYPE (LRB_TYPE), INTENT(OUT), 
     &          DIMENSION(max(NB_BLOCK_U,1)):: BLR_U
      INTEGER, INTENT(OUT), DIMENSION(NB_BLOCK_U+2)  :: BEGS_BLR_U 
      INTEGER(8) :: KEEP8(150)
      LOGICAL :: ISLR
      INTEGER :: ISLR_INT, I
      INTEGER :: K, M, N
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      IERR = 0
      IF (size(BLR_U) .NE. 
     &     MAX(NB_BLOCK_U,1) ) THEN
        WRITE(*,*) "Internal error 1 in DMUMPS_MPI_UNPACK",
     &             NB_BLOCK_U,size(BLR_U)
        CALL MUMPS_ABORT()
      ENDIF
      BEGS_BLR_U(1) = 1
      BEGS_BLR_U(2) = NPIV+NELIM+1 
      DO I = 1, NB_BLOCK_U
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 ISLR_INT, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 K, 1,
     &                 MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 M, 1,
     &                 MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 N, 1,
     &                 MPI_INTEGER, COMM, IERR )
        BEGS_BLR_U(I+2) = BEGS_BLR_U(I+1) + M
        IF (ISLR_INT .eq. 1) THEN
          ISLR = .TRUE.
        ELSE
          ISLR = .FALSE.
        ENDIF
        CALL ALLOC_LRB( BLR_U(I), K, M, N, ISLR, 
     &             IFLAG, IERROR, KEEP8 )
        IF (IFLAG.LT.0) RETURN
        IF (ISLR) THEN
          IF (K .GT. 0) THEN
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                     BLR_U(I)%Q(1,1), M*K, MPI_DOUBLE_PRECISION,
     &                     COMM, IERR )
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                     BLR_U(I)%R(1,1), N*K, MPI_DOUBLE_PRECISION,
     &                     COMM, IERR)
          ENDIF
        ELSE
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                     BLR_U(I)%Q(1,1), M*N, MPI_DOUBLE_PRECISION,
     &                     COMM, IERR)
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_MPI_UNPACK_LR
      SUBROUTINE DMUMPS_MPI_UNPACK_LR_PARTIAL(
     &           BUFR, LBUFR, LBUFR_BYTES, POSITION,
     &                             NPIV, NELIM, DIR, 
     &                             BLR_U, NB_BLOCK_U, 
     &                             JBEG_BLOCK, 
     &                             BEGS_BLR_U, KEEP8,
     &                             COMM, IERR, IFLAG, IERROR)
      USE DMUMPS_LR_CORE, ONLY : LRB_TYPE, ALLOC_LRB
      USE DMUMPS_LR_TYPE
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LBUFR
      INTEGER, INTENT(IN) :: LBUFR_BYTES
      INTEGER, INTENT(IN) :: BUFR(LBUFR)
      INTEGER, INTENT(INOUT) :: POSITION
      INTEGER, INTENT(IN)    :: NB_BLOCK_U, NELIM, NPIV, JBEG_BLOCK
      CHARACTER(len=1) :: DIR
      INTEGER, INTENT(IN) :: COMM
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
      INTEGER, INTENT(OUT) :: IERR
      TYPE (LRB_TYPE), INTENT(OUT), 
     &          DIMENSION(max(NB_BLOCK_U,1)):: BLR_U
      INTEGER, INTENT(OUT), DIMENSION(NB_BLOCK_U+2)  :: BEGS_BLR_U 
      INTEGER(8) :: KEEP8(150)
      LOGICAL :: ISLR
      INTEGER :: ISLR_INT, I
      INTEGER :: K, M, N
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      IERR = 0
      IF (size(BLR_U) .NE. 
     &     MAX(NB_BLOCK_U,1) ) THEN
        WRITE(*,*) "Internal error 1 in DMUMPS_MPI_UNPACK",
     &             NB_BLOCK_U,size(BLR_U)
        CALL MUMPS_ABORT()
      ENDIF
      BEGS_BLR_U(1) = 1
      BEGS_BLR_U(2) = NPIV+NELIM+1 
      IF (JBEG_BLOCK.NE.1) THEN
        BEGS_BLR_U(2) = JBEG_BLOCK
      ENDIF
      DO I = 1, NB_BLOCK_U
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 ISLR_INT, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 K, 1,
     &                 MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 M, 1,
     &                 MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 N, 1,
     &                 MPI_INTEGER, COMM, IERR )
        BEGS_BLR_U(I+2) = BEGS_BLR_U(I+1) + M
        IF (ISLR_INT .eq. 1) THEN
          ISLR = .TRUE.
        ELSE
          ISLR = .FALSE.
        ENDIF
        CALL ALLOC_LRB( BLR_U(I), K, M, N, ISLR, 
     &             IFLAG, IERROR, KEEP8 )
        IF (IFLAG.LT.0) RETURN
        IF (ISLR) THEN
          IF (K .GT. 0) THEN
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                     BLR_U(I)%Q(1,1), M*K, MPI_DOUBLE_PRECISION,
     &                     COMM, IERR )
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                     BLR_U(I)%R(1,1), N*K, MPI_DOUBLE_PRECISION,
     &                     COMM, IERR)
          ENDIF
        ELSE
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                     BLR_U(I)%Q(1,1), M*N, MPI_DOUBLE_PRECISION,
     &                     COMM, IERR)
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_MPI_UNPACK_LR_PARTIAL
