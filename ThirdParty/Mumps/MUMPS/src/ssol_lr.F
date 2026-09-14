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
      MODULE SMUMPS_SOL_LR
      USE SMUMPS_LR_TYPE
      USE SMUMPS_LR_CORE
      USE MUMPS_LR_STATS
      USE SMUMPS_LR_DATA_M, only: BLR_ARRAY
      USE MUMPS_L0_OMP_M, ONLY : thread_id
      USE ISO_C_BINDING
!$    USE OMP_LIB
      IMPLICIT NONE
      CONTAINS
C     INCOMPATIBILITIES:
C     LRSOLVE is not compatible with old option KEEP(350) = 0
C              and is not compatible with (RHSINTR_BYROWS)
      SUBROUTINE SMUMPS_SOL_FWD_LR_SU 
     &    (INODE, N, IWHDLR, NPIV_GLOBAL, NSLAVES,
     &     IW, IPOS_INIT, LIW, 
     &     LIELL, WCB, LWCB, 
     &     LD_WCBPIV, LD_WCBCB,
     &     PPIV_INIT, PCB_INIT, 
     &     RHSINTR, LRHSINTR, NRHS, 
     &     POSINRHSINTR_FWD, JBDEB, JBFIN, 
     &     MTYPE, KEEP, KEEP8, OOCWRITE_COMPATIBLE_WITH_BLR,
     &     IFLAG, IERROR )
C
C  Purpose:
C  =======
C
C      Compute LY=B or UTY=B of a BLR node 
C      Solve and Update are performed by panel of BLR blocks
C
C      At each iteration (BLR block of the FS block) of the loop :
C       we process a COMPLETE panel (block pivot and OFF diag block) 
C       We first compute a TRSV/TRSM on the pivot block 
C       We then update the COMPLETE RHS with GEMV/GEMM
C
C      NPIV_GLOBAL is needed to know when WCB changes LDA from LD_WCBPIV
C      to LD_WCBCB
C  
       INTEGER, INTENT(IN) :: INODE, N, IWHDLR, NPIV_GLOBAL, NSLAVES
       INTEGER, INTENT(IN) :: MTYPE, LIELL, KEEP(500)
       INTEGER(8), INTENT(IN) :: KEEP8(150)
       INTEGER, INTENT(IN) :: LIW, IPOS_INIT, LRHSINTR
       INTEGER, INTENT(IN) :: IW(LIW), POSINRHSINTR_FWD(N)
       INTEGER(8), INTENT(IN) :: LWCB, PPIV_INIT, PCB_INIT
       INTEGER, INTENT(IN) :: LD_WCBPIV, LD_WCBCB, NRHS, JBDEB, JBFIN
       REAL, INTENT(INOUT) :: WCB(LWCB)
       INTEGER, INTENT(INOUT) :: IFLAG, IERROR
       REAL, INTENT(INOUT) :: RHSINTR(LRHSINTR, NRHS)
       LOGICAL, INTENT(IN)    :: OOCWRITE_COMPATIBLE_WITH_BLR
C      MTYPE
C             1 = Ly=b bloc of factors
C             else  UT y = b
C      LD_WCBPIV: Leading dimension of WCB for pivot block 
C      LD_WCBCB : Leading dimension of WCB for CB block
C                 including NELIM(panel) variables
C      NSLAVES 
C             0  = Type 1 node
C             >0 = Type 2 node
C
C Local variables
        INTEGER :: I, NPARTSASS, NB_BLR , NELIM, LDADIAG,
     &             DIAGSIZ_DYN, DIAGSIZ_STA, IBEG_BLR, IEND_BLR,
     &             LD_CB, NRHS_B, IPOS, KCB
        INTEGER(8) :: PPIV, PCB
        INTEGER    :: LAST_BLR
        REAL, POINTER, DIMENSION(:) :: DIAG
        TYPE(LRB_TYPE), POINTER, DIMENSION(:) :: BLR_PANEL
#if !defined(NOMIXHPC)
        REAL, ALLOCATABLE, TARGET :: WORKSPACE_INCASE(:)
        REAL, POINTER :: WORKSPACE(:)
        INTEGER :: WSMAX
        INTEGER :: CACHE_L2
        integer(c_intptr_t) :: th_ws_int
        type(c_ptr), target :: th_ws_ptr
#endif
        integer :: th_num
C
C    Constants
C    =========
      REAL :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
      PARAMETER (ZERO=0.0E0)
C       Current block size of RHS beeing processed
        NRHS_B = JBFIN-JBDEB+1
C
        IF (MTYPE.EQ.1) THEN
C       Uncompress L panels in A
         IF (associated(BLR_ARRAY(IWHDLR)%PANELS_L))
     &      THEN
            NPARTSASS=size(BLR_ARRAY(IWHDLR)%PANELS_L)
            NB_BLR = size(BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC) -1
         ELSE
           WRITE(6,*) " Internal error in SMUMPS_SOL_FWD_SU_MASTER"
         ENDIF
        ELSE
         IF (associated(BLR_ARRAY(IWHDLR)%PANELS_U))
     &      THEN
            NPARTSASS=size(BLR_ARRAY(IWHDLR)%PANELS_U)
            NB_BLR = size(BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC) -1
         ENDIF
        ENDIF
C       Set LAST_BLR for call to SMUMPS_SOL_FWD_BLR_UPDATE
        IF (NSLAVES.EQ.0 .OR. (KEEP(50).eq.0 .and. MTYPE .NE.1)) THEN
C         Type 1 or unsym transposed type 2:
C         update using until last (NB_BLR) block
          LAST_BLR = NB_BLR
        ELSE
C         Type 2 symmetric: 
          LAST_BLR = NPARTSASS
        ENDIF
C
C Pointeur
C
C     IW(IPOS_INIT+1): starting of index list of the front
      IPOS = IPOS_INIT
      PPIV = PPIV_INIT 
      th_num = 0
      if (KEEP(405).EQ.1) then
          th_num = thread_id - 1
      endif
#if !defined(NOMIXHPC)
      WSMAX = 0
      CACHE_L2 = KEEP(365) / COSTS_PREC(1)
      IF (PREC_DIAG_FR.gt.1) THEN
        DO I=1, NPARTSASS
          WSMAX = MAX(WSMAX,
     &        BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE)
        END DO
        IF (WSMAX.GT.CACHE_L2) THEN
          ALLOCATE(WORKSPACE_INCASE(WSMAX))
          WORKSPACE => WORKSPACE_INCASE
        ELSE
          call mumps_sol_omp_mem_manager_get_address(th_num, th_ws_int)
          th_ws_ptr = transfer(th_ws_int, th_ws_ptr)
          call c_f_pointer(th_ws_ptr, WORKSPACE, (/ CACHE_L2 /))
        ENDIF
      ENDIF
#endif
      DO I=1, NPARTSASS
        IBEG_BLR = BLR_ARRAY(IWHDLR)%BEGS_BLR_DYNAMIC(I)
        IEND_BLR = BLR_ARRAY(IWHDLR)%BEGS_BLR_DYNAMIC(I+1) -1
        DIAGSIZ_DYN = BLR_ARRAY(IWHDLR)%BEGS_BLR_DYNAMIC(I+1) -
     &           IBEG_BLR
        DIAGSIZ_STA = BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC(I+1)  - 
     &           IBEG_BLR 
        IF (KEEP(50).NE.0) THEN
         LDADIAG = DIAGSIZ_DYN
        ELSE
         LDADIAG = DIAGSIZ_STA
        ENDIF 
        IF (IEND_BLR.EQ.NPIV_GLOBAL) THEN
          PCB   = PCB_INIT
        ELSE
          PCB  = PPIV + int(DIAGSIZ_DYN,8)
        ENDIF
C       Skip empty blocks
        IF ( DIAGSIZ_DYN.EQ.0) CYCLE
C       NELIM has to be defined by panel
        NELIM = DIAGSIZ_STA - DIAGSIZ_DYN
C       or also NELIM =    BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC(I+1) -
C     &          BLR_ARRAY(IWHDLR)%BEGS_BLR_DYNAMIC(I+1)
         IF ( MTYPE .EQ. 1 ) THEN
          BLR_PANEL => BLR_ARRAY(IWHDLR)%PANELS_L(I)%LRB_PANEL
         ELSE
C        ---------
C        U^T y = b
C        ---------
C         Since U^T=QR already => nothing special to do
          BLR_PANEL => BLR_ARRAY(IWHDLR)%PANELS_U(I)%LRB_PANEL
         END IF
C       ===========================================================
C       I/ Process Diagonal block : 
C           I.0/ Cast to working precision if necessary
C           I.1/ Triangular compute part of solution
C           I.2/ Part of the update of the RHS concerning the NELIM
C                var at this step
C       ===========================================================
C        I.0
         DIAG => BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_BLOCK 
#if !defined(NOMIXHPC)
         IF (PREC_DIAG_FR.gt.1) THEN
           DIAG(1:BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE) => 
     &       WORKSPACE(1:BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE)
           CALL SMUMPS_LRMIXED_CAST(
     &       BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_STOC%VAL(1), 
     &       DIAG(1),
     &       BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE, 1,
     &       BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE,
     &       BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE,
     &       PREC_DIAG_FR, 1)
         ELSE IF (.not.ASSOCIATED(DIAG)) THEN
           call c_f_pointer(
     &        c_loc(BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_STOC%VAL),
     &        DIAG, (/ BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE/) )
         ENDIF
#endif
C        I.1:
         CALL SMUMPS_SOLVE_FWD_TRSOLVE (DIAG(1), 
     &         int(size(DIAG),8), 1_8,
     &         DIAGSIZ_DYN , LDADIAG, NRHS_B, WCB, LWCB, NPIV_GLOBAL,
     &         PPIV, MTYPE, KEEP)
C       I.2 
         IF (NELIM.GT.0) THEN
C       
           KCB = int(PCB-PPIV_INIT+1)
C          Update NELIM block
           IF (IEND_BLR.EQ.NPIV_GLOBAL) THEN 
C            Switch from PIV to CB part             
             LD_CB = LD_WCBCB
           ELSE
             LD_CB = LD_WCBPIV
           ENDIF
C          LD of the diag block is DIAGSIZ_DYN on the part concerned           
           IF (MTYPE.EQ.1) THEN
             IF (KCB.LE.NPIV_GLOBAL .AND.
     &          KCB+NELIM-1.GT.NPIV_GLOBAL) THEN
               CALL sgemm('T', 'N', NPIV_GLOBAL-KCB+1, NRHS_B,
     &              DIAGSIZ_DYN, MONE,
     &              DIAG(1+DIAGSIZ_DYN*LDADIAG), DIAGSIZ_DYN, 
     &              WCB(PPIV), LD_WCBPIV,
     &              ONE, WCB(PCB), LD_CB) 
               CALL sgemm('T', 'N', KCB+NELIM-NPIV_GLOBAL-1,
     &              NRHS_B, DIAGSIZ_DYN, MONE,
     &              DIAG(1+DIAGSIZ_DYN*LDADIAG +
     &              (NPIV_GLOBAL-KCB+1)*DIAGSIZ_DYN),
     &              DIAGSIZ_DYN, 
     &              WCB(PPIV), LD_WCBPIV,
     &              ONE, WCB(PCB_INIT), LD_WCBCB) 
             ELSE
               CALL sgemm('T', 'N', NELIM, NRHS_B, DIAGSIZ_DYN, MONE,
     &              DIAG(1+DIAGSIZ_DYN*LDADIAG), DIAGSIZ_DYN, 
     &              WCB(PPIV), LD_WCBPIV,
     &              ONE, WCB(PCB), LD_CB) 
             ENDIF
           ELSE
             IF (KCB.LE.NPIV_GLOBAL .AND.
     &          KCB+NELIM-1.GT.NPIV_GLOBAL) THEN
               CALL sgemm('N', 'N', NPIV_GLOBAL-KCB+1,
     &              NRHS_B, DIAGSIZ_DYN, MONE,
     &              DIAG(1+DIAGSIZ_DYN), DIAGSIZ_STA, 
     &              WCB(PPIV), LD_WCBPIV,
     &              ONE, WCB(PCB), LD_CB) 
               CALL sgemm('N', 'N', KCB+NELIM-NPIV_GLOBAL-1,
     &              NRHS_B, DIAGSIZ_DYN, MONE,
     &              DIAG(1+DIAGSIZ_DYN+NPIV_GLOBAL-KCB+1), 
     &              DIAGSIZ_STA, 
     &              WCB(PPIV), LD_WCBPIV,
     &              ONE, WCB(PCB_INIT), LD_WCBCB) 
             ELSE
               CALL sgemm('N', 'N', NELIM, NRHS_B, DIAGSIZ_DYN, MONE,
     &              DIAG(1+DIAGSIZ_DYN), DIAGSIZ_STA, 
     &              WCB(PPIV), LD_WCBPIV,
     &              ONE, WCB(PCB), LD_CB) 
             ENDIF
           ENDIF
         ENDIF
C       ===========================================================
C       II/ Update RHS with computed solution
C        CALL SMUMPS_SOL_FWD_BLR_UPDATE 
C       ===========================================================
C       LAST_BLR=NPARTSASS or NB_BLR (see above,
C       outside parallel region) is the last panel
C       block to use for the update
        CALL SMUMPS_SOL_FWD_BLR_UPDATE (
     &        WCB, LWCB, 1, LD_WCBPIV, PPIV_INIT, 1,
     &        WCB, LWCB, LD_WCBCB, PCB_INIT,
     &        PPIV, !POSDIAG   
     &        NRHS_B, NPIV_GLOBAL,
     &        BLR_PANEL, LAST_BLR, I,
     &        BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC,
     &        KEEP8, KEEP(34), 
#if !defined(NOMIXHPC)
     &        KEEP(148),
     &        KEEP(365),
     &        KEEP(405), ! I need to know if I'm under L0
#endif
     &        KEEP(448), KEEP(450), KEEP(451), 
#if !defined(NOMIXHPC)
     &        KEEP(455),
#endif
     &        .FALSE.,   ! IS_T2_SLAVE
     &        IFLAG, IERROR)
        IF (IFLAG.LT.0) RETURN
C
C       ===========================================================
C       III/ 
C         Perform Y = DIAG * X for symmetric matrices (A = LDL^T)
C         then Reload WCB back to RHSINTR
C       ===========================================================
#if !defined(NOMIXHPC)
         IF (PREC_DIAG_FR.gt.1) THEN
           DIAG(1:BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE) => 
     &       WORKSPACE(1:BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE)
           CALL SMUMPS_LRMIXED_CAST(
     &       BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_STOC%VAL(1), DIAG(1),
     &       BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE, 1,
     &       BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE,
     &       BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE,
     &       PREC_DIAG_FR, 1)
         ELSE IF (.not.ASSOCIATED(DIAG)) THEN
           call c_f_pointer(
     &        c_loc(BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_STOC%VAL),
     &        DIAG, (/ BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE/) )
         ENDIF
#endif
        CALL SMUMPS_SOL_LD_AND_RELOAD_PANEL (
     &  INODE, N, DIAGSIZ_DYN, LIELL, NELIM, NSLAVES,
     &  PPIV,
     &  IW, IPOS, LIW, 
     &  DIAG(1), 
     &  int(size(DIAG),8), 1_8, 
!     &  BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE, 1_8,
     &  WCB, LWCB, LD_WCBPIV, 
     &  RHSINTR, LRHSINTR, NRHS, 
     &  POSINRHSINTR_FWD, JBDEB, JBFIN, 
     &  MTYPE, KEEP, OOCWRITE_COMPATIBLE_WITH_BLR,
     &  .TRUE.  ! ignore KEEP(459) -> do not subdivide into panels
     &  )
        PPIV = PPIV +  int(DIAGSIZ_DYN,8)
        IPOS = IPOS + DIAGSIZ_DYN
C     -- end of loops on panels
      ENDDO
#if !defined(NOMIXHPC)
      IF (ALLOCATED(WORKSPACE_INCASE)) DEALLOCATE(WORKSPACE_INCASE)
#endif
      RETURN
      END SUBROUTINE SMUMPS_SOL_FWD_LR_SU
      SUBROUTINE SMUMPS_SOL_SLAVE_LR_U
     &    (INODE, IWHDLR, NPIV_GLOBAL,
     &     WCB, LWCB, 
     &     LDX, LDY,
     &     PTRX_INIT, PTRY_INIT, 
     &     JBDEB, JBFIN, 
     &     MTYPE, KEEP, KEEP8, IFLAG, IERROR )
C
C  Purpose:
C  =======
C
C      Compute Y <- Y - L X or L^T X (depending on MTYPE) where L is
C      the rectangular matrix of a type 2 slave stored in BLR format.
C
C      At each iteration of the loop, we process a panel
C      and update Y.
C
C      This routine is called for type 2 slaves (fwd and bwd).
C
C  Arguments:
C  =========
C
C      WCB(PTRX_INIT): points to X, leading dimension LDX
C      WCB(PTRY_INIT): points to Y, leading dimension LDY
C      MTYPE         : 1: Y <- Y - L X
C                      0: Y <- Y - L^T X
C      
       INTEGER, INTENT(IN) :: INODE, IWHDLR, NPIV_GLOBAL
       INTEGER, INTENT(IN) :: MTYPE, KEEP(500)
       INTEGER(8), INTENT(IN) :: KEEP8(150)
       INTEGER(8), INTENT(IN) :: LWCB, PTRX_INIT, PTRY_INIT
       INTEGER, INTENT(IN) :: LDX, LDY, JBDEB, JBFIN
       REAL, INTENT(INOUT) :: WCB(LWCB)
       INTEGER, INTENT(INOUT) :: IFLAG, IERROR
C
C  Local variables
C  ===============
C
       INTEGER :: I, NPARTSASS, NB_BLR , NRHS_B
       INTEGER(8) :: PTRX, PTRY
       TYPE(LRB_TYPE), POINTER, DIMENSION(:) :: BLR_PANEL
       INTEGER :: M, PREC, NB_INT_B, POS_B, IBEG_BLR, 
     &            NROWS, K34
       INTEGER :: allocok
       TYPE(LRMIXED_TYPE), ALLOCATABLE :: PT_RHS(:,:)
       INTEGER, ALLOCATABLE, TARGET, DIMENSION(:) :: TEMP_RHS
C
C    Constants
C    =========
      REAL :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
      PARAMETER (ZERO=0.0E0)
C     Current block size of RHS beeing processed
      NRHS_B = JBFIN-JBDEB+1
C
      IF (associated(BLR_ARRAY(IWHDLR)%PANELS_L))
     &   THEN
         NPARTSASS=size(BLR_ARRAY(IWHDLR)%PANELS_L)
         NB_BLR = size(BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC)
         NB_BLR = NB_BLR - 2
      ELSE
        WRITE(6,*) " Internal error 1 in SMUMPS_SOL_SLAVE_LR_U"
        CALL MUMPS_ABORT()
      ENDIF
C     Copy RHS in each precision
      IF (KEEP(448).GT.0 .AND. MTYPE.NE.1 .AND. KEEP(450).GT.0) THEN
        K34 = KEEP(34)
        NROWS = BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC(NB_BLR+2) 
     &          - BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC(2)
C       Count size of all copies of RHS in each prec        
        NB_INT_B = 0
        DO PREC=2,NPREC_CALC
          NB_INT_B = NB_INT_B + ( NROWS * NRHS_B *
     &                        COSTS_PREC(CALC2STOC(PREC)) +K34-1) / K34
        ENDDO
C       Workspace containing all copies         
        allocate(TEMP_RHS(NB_INT_B), stat=allocok)
        if (allocok.gt.0) then
          IFLAG = -13
          IERROR = NB_INT_B
          write(*,*) 'Allocation problem in BLR routine 
     &    SMUMPS_SOL_SLAVE_LR_U',
     &    'not enough memory? memory requested = ', IERROR
          CALL MUMPS_ABORT()
        endif
C       Pointer that points to each block copy in each precision
        allocate(PT_RHS(NPREC_CALC,NB_BLR), stat=allocok)
        if (allocok.gt.0) then
          IFLAG = -13
          IERROR = NPREC_CALC * NB_BLR
          write(*,*) 'Allocation problem in BLR routine 
     &    SMUMPS_SOL_SLAVE_LR_U',
     &    'not enough memory? memory requested = ', IERROR
          CALL MUMPS_ABORT()
        endif
C       Nullify unused pointers (PREC=1)
        DO I=1,NB_BLR
          nullify(PT_RHS(1,I)%val)
        ENDDO
C       Cast RHS
        POS_B = 1
        DO I = 1,NB_BLR
          IBEG_BLR = BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC(I+1)
          M = BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC(I+2) 
     &                  - IBEG_BLR
          PTRX = PTRX_INIT +int(IBEG_BLR-1,8)
          DO PREC=2,NPREC_CALC
            PT_RHS(PREC,I)%val => TEMP_RHS(POS_B:)
            POS_B = POS_B + (M * NRHS_B * 
     &                 COSTS_PREC(CALC2STOC(PREC)) +K34-1) / K34
            CALL SMUMPS_LRMIXED_CAST(WCB(PTRX),
     &          PT_RHS(PREC,I)%val(1), M, 
     &          NRHS_B, LDX, M, 1, CALC2STOC(PREC))
          ENDDO
        ENDDO
#if defined(MUMPS_NOF2003)
      ELSE
C       PT_RHS can not be declared allocatable in
C       SMUMPS_SOL_BWD_BLR_UPDATE, it is declared
C       as an array and so it must be allocated.
        ALLOCATE(PT_RHS(1,1))
#endif
      ENDIF
      PTRX = PTRX_INIT
      PTRY = PTRY_INIT
      DO I = 1, NPARTSASS
        BLR_PANEL => BLR_ARRAY(IWHDLR)%PANELS_L(I)%LRB_PANEL
        IF (associated(BLR_PANEL)) THEN
C         =================================
C         Update RHS with computed solution
C         =================================
          IF (MTYPE.EQ.1) THEN ! Fwd
            CALL SMUMPS_SOL_FWD_BLR_UPDATE (
     &          WCB, LWCB, 1, LDX, -99999_8, 1,
     &          WCB, LWCB, LDY, PTRY,
     &          PTRX,
     &          NRHS_B, NPIV_GLOBAL,
     &          BLR_PANEL, NB_BLR, 0, ! update starts at CURRENT_PANEL+1
     &          BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC(2:NB_BLR+2),
     &          KEEP8, KEEP(34), 
#if !defined(NOMIXHPC)
     &          KEEP(148),
     &          KEEP(365),
     &        KEEP(405), ! I need to know if I'm under L0
#endif
     &          KEEP(448), KEEP(450), KEEP(451), 
#if !defined(NOMIXHPC)
     &          KEEP(455),
#endif
     &          .TRUE., IFLAG, IERROR )
          ELSE !Bwd
            IF (.NOT.allocated(PT_RHS))  THEN
               ALLOCATE(PT_RHS(1,1))
            ENDIF
            CALL SMUMPS_SOL_BWD_BLR_UPDATE (
c    &          WCB, LWCB, 1, LDX, -99999, 1,
c    &          WCB, LWCB, LDY, PTRY,
c    &          PTRX,
     &          WCB, LWCB, 1, LDY, -99999_8, 1,
     &          WCB, LWCB, LDX, PTRX,
     &          PTRY,
     &          NRHS_B, NPIV_GLOBAL,
     &          BLR_PANEL, NB_BLR, 0, ! update starts at CURRENT_PANEL+1
     &          BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC(2:NB_BLR+2),
     &          KEEP8, KEEP(34), 
#if !defined(NOMIXHPC)
     &          KEEP(148),
     &          KEEP(365),
     &        KEEP(405), ! I need to know if I'm under L0
#endif
     &          KEEP(448), KEEP(450), KEEP(451),
#if !defined(NOMIXHPC)
     &          KEEP(455),
#endif
     &          .TRUE.,
     &          PT_RHS, 
     &          IFLAG, IERROR )
          ENDIF
          IF (MTYPE .EQ. 1) THEN
            PTRX = PTRX + BLR_PANEL(1)%N
          ELSE
            PTRY = PTRY + BLR_PANEL(1)%N
          ENDIF
          IF (IFLAG.LT.0) RETURN
        ENDIF
C
C     -- end of loops on panels
      ENDDO  
      IF (allocated(PT_RHS))   deallocate(PT_RHS)
      IF (allocated(TEMP_RHS)) deallocate(TEMP_RHS)
      RETURN
      END SUBROUTINE SMUMPS_SOL_SLAVE_LR_U
      SUBROUTINE GEMM_Q_FWD(m, nrhs_b, k, npiv, 
     &          Q, TMP, ldT,
     &          arraypiv, ldpiv, arraycb, lcb, ldcb, 
     &          ibeg_block, iend_block, is_t2_slave, 
     &          poscb, pospiv, pospivcol, ibeg_tmp)
          implicit none
C         arguments
          integer, intent(in) :: m, nrhs_b, k, npiv
          REAL, dimension(:,:), intent(inout)     :: Q
          REAL, dimension(ldt, *), intent(inout)     :: TMP
          integer(8), intent(in) :: lcb
          integer, intent(in)    :: ldpiv
          REAL, intent(inout) :: arraypiv(ldpiv,*)
          REAL, intent(inout) :: arraycb(lcb)
          integer, intent(in) :: ldt, ldcb
          integer, intent(in) :: ibeg_block, iend_block
          logical, intent(in) :: is_t2_slave
          integer(8), intent(in) :: poscb, pospiv
          integer, intent(in)    :: pospivcol
          integer, intent(in) :: ibeg_tmp
C         local variables
          integer :: posblock
          REAL :: ONE, MONE, ZERO
          PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
          PARAMETER (ZERO=0.0E0)
          IF (IS_T2_SLAVE) THEN
            POSBLOCK = POSCB+int(IBEG_BLOCK-1,8)
            CALL sgemm('N', 'N', M, NRHS_B, K, MONE,
     &        Q(1,1), M, TMP(ibeg_tmp,1), LDT,
     &        ONE, ARRAYCB(POSBLOCK), LDCB) 
          ELSEIF (IBEG_BLOCK.LE.NPIV.AND.
     &            IEND_BLOCK.GT.NPIV) THEN
C           Case of overlap between the FS and non FS
            POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
            CALL sgemm('N', 'N', NPIV-IBEG_BLOCK+1,NRHS_B, K,
     &        MONE, Q(1,1), M, TMP(ibeg_tmp,1), LDT,
     &        ONE, ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV) 
            CALL sgemm('N', 'N', IBEG_BLOCK+M-NPIV-1, 
     &        NRHS_B, K, 
     &        MONE, Q(NPIV-IBEG_BLOCK+2,1), M, 
     &        TMP(ibeg_tmp,1), LDT, ONE, ARRAYCB(POSCB), LDCB) 
          ELSEIF (IBEG_BLOCK.LE.NPIV) THEN
C           Entirely on FS side              
            POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
            CALL sgemm('N', 'N', M, NRHS_B, K, MONE,
     &        Q(1,1), M, TMP(ibeg_tmp,1), LDT,
     &        ONE, ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV) 
          ELSE
C           Entirely on CB side              
            POSBLOCK = POSCB+int(IBEG_BLOCK-1-NPIV,8)
            CALL sgemm('N', 'N', M, NRHS_B, K, MONE,
     &        Q(1,1), M, TMP(ibeg_tmp,1), LDT,
     &        ONE, ARRAYCB(POSBLOCK), LDCB) 
          ENDIF
      RETURN
      END SUBROUTINE GEMM_Q_FWD
      ! TODO it would be better to have these routines to factor some
      ! convoluted calls
C      SUBROUTINE GEMM_LR_FWD(
C                m, nrhs_b, k, npiv, 
C     &          Q, TMP, ldT,
C     &          arraypiv, ldpiv, arraycb, lcb, ldcb, 
C     &          ibeg_block, iend_block, is_t2_slave, 
C     &          poscb, pospiv, pospivcol, ibeg_tmp, prec)
C          INTEGER(8), INTENT(IN) :: LPIV, LCB, POSPIV, POSCB, POSDIAG
C          INTEGER, INTENT(IN)    :: LAST_BLR, NRHS_B,  LDPIV, LDCB,
C     &                              CURRENT_BLR, NPIV, K34, PREC
C          INTEGER, INTENT(IN)    :: LPIVCOL, POSPIVCOL, KINT, KREM, M
C          integer(8), intent(in) :: poscb, pospiv
C          logical, intent(in) :: is_t2_slave
C          REAL, INTENT(INOUT) :: ARRAYPIV(LPIV,LPIVCOL)
C          REAL, INTENT(INOUT) :: ARRAYCB(LCB)
C          INTEGER, intent(in) :: IBEG_BLOCK, IEND_BLOCK
C! local variables
C          REAL :: ONE, MONE, ZERO
C          INTEGER(8) :: POSBLOCK
C#if 0
C          PARAMETER (ONE=(1.0E0,0.0E0), MONE=(-1.0E0,0.0E0))
C          PARAMETER (ZERO=(0.0E0,0.0E0))
C#else
C          PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
C          PARAMETER (ZERO=0.0E0)
C#endif
C          IF (IS_T2_SLAVE) THEN
C            POSBLOCK = POSCB + int(IBEG_BLOCK-1,8)
C            CALL sgemm('N', 'N', M, NRHS_B, KREM, MONE,
C     &            Q_ptr(1,1), M, 
C     &            ARRAYPIV(POSDIAG+KINT-1,POSPIVCOL), LDPIV,
C     &            ONE, ARRAYCB(POSBLOCK), LDCB) 
C          ELSEIF (IBEG_BLOCK.LE.NPIV.AND.IEND_BLOCK.GT.NPIV) THEN
CC           Case of overlap between the FS and non FS
C            POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
C            CALL sgemm('N', 'N', NPIV-IBEG_BLOCK+1, NRHS_B, 
C     &            KREM, MONE,
C     &            Q_ptr(1,1), M,
C     &            ARRAYPIV(POSDIAG+KINT-1,POSPIVCOL), LDPIV,
C     &            ONE, ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV) 
C            CALL sgemm('N', 'N', IBEG_BLOCK+M-NPIV-1, NRHS_B,
C     &            KREM, MONE,
C     &            Q_ptr(NPIV-IBEG_BLOCK+2,1), M,
C     &            ARRAYPIV(POSDIAG+KINT-1,POSPIVCOL),
C     &            LDPIV, ONE, ARRAYCB(POSCB), LDCB) 
C          ELSEIF (IBEG_BLOCK.LE.NPIV) THEN
CC           Entirely on FS side              
C            POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
C            CALL sgemm('N', 'N', M, NRHS_B, KREM, MONE,
C     &            Q_ptr(1,1), M,
C     &            ARRAYPIV(POSDIAG+KINT-1,POSPIVCOL),LDPIV,
C     &            ONE, ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV)
C          ELSE
CC           Entirely on CB side              
C            POSBLOCK = POSCB + int(IBEG_BLOCK-1-NPIV,8)
C            CALL sgemm('N', 'N', M, NRHS_B, KREM, MONE,
C     &            Q_ptr(1,1), M,
C     &            ARRAYPIV(POSDIAG+KINT-1,POSPIVCOL),LDPIV,
C     &            ONE, ARRAYCB(POSBLOCK), LDCB) 
C          ENDIF
C      END SUBROUTINE GEMM_LR_FWD
      SUBROUTINE GEMM_Q_BWD(m, nrhs_b, k, npiv, 
     &          Q, TMP, ldT,
     &          arraypiv, lpiv, ldpiv, arraycb, lcb, ldcb, 
     &          ibeg_block, iend_block, is_t2_slave, 
     &          poscb, pospiv, pospivcol, ibeg_tmp)
          implicit none
C         arguments
          integer, intent(in) :: m, nrhs_b, k, npiv
          REAL, dimension(:, :), intent(inout)     :: Q
          REAL, dimension(ldt, *), intent(inout)     :: TMP
          integer(8), intent(in) :: lcb, lpiv
          REAL, intent(inout) :: arraypiv(lpiv,*)
          REAL, intent(inout) :: arraycb(lcb)
          integer, intent(in) :: ldt, ldcb, ldpiv
          integer, intent(in) :: ibeg_block, iend_block
          logical, intent(in) :: is_t2_slave
          integer(8), intent(in) :: poscb, pospiv
          integer, intent(in)    :: pospivcol
          integer, intent(in) :: ibeg_tmp
C         local variables
          integer(8) :: posblock
          REAL :: ONE, MONE, ZERO
          PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
          PARAMETER (ZERO=0.0E0)
          IF (IS_T2_SLAVE) THEN
            POSBLOCK = POSCB +int(IBEG_BLOCK-1,8)
            CALL sgemm('T', 'N', K, NRHS_B, M, ONE,
     &        Q(1,1), M,  
     &        ARRAYCB(POSBLOCK), LDCB, ZERO, 
     &        TMP(ibeg_tmp,1), ldt)
          ELSEIF (IBEG_BLOCK.LE.NPIV.AND.
     &            IEND_BLOCK.GT.NPIV) THEN
C           Case of overlap between the FS and non FS
            POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
            CALL sgemm('T', 'N', K, NRHS_B, 
     &        NPIV-IBEG_BLOCK+1,
     &        ONE, Q(1,1), M,  
     &        ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV,
     &        ZERO, TMP(ibeg_tmp, 1), ldt)
            CALL sgemm('T', 'N', 
     &        K, NRHS_B, IBEG_BLOCK+M-NPIV-1,
     &        ONE, Q(NPIV-IBEG_BLOCK+2,1), M,  
     &        ARRAYCB(POSCB), LDCB, 
     &        ONE,
     &        TMP(ibeg_tmp,1), ldt)
          ELSEIF (IBEG_BLOCK.LE.NPIV) THEN
C           Entirely on FS side              
            POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
            CALL sgemm('T', 'N', K, NRHS_B, M, ONE,
     &        Q(1,1), M,  
     &        ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV,
     &        ZERO, TMP(ibeg_tmp, 1), ldt)
          ELSE
C           Entirely on CB side              
            POSBLOCK = POSCB+int(IBEG_BLOCK-1-NPIV,8)
            CALL sgemm('T', 'N', K, NRHS_B, M, ONE,
     &        Q(1,1), M,  
     &        ARRAYCB(POSBLOCK), LDCB, ZERO, 
     &        TMP(ibeg_tmp, 1), ldt)
          ENDIF
      RETURN
      END SUBROUTINE GEMM_Q_BWD
#if !defined(NOMIXHPC)
      SUBROUTINE SMUMPS_SOL_FWD_BLR_UPDATE (
     &        ARRAYPIV, LPIV, LPIVCOL, LDPIV, POSPIV, POSPIVCOL,
     &        ARRAYCB, LCB, LDCB, POSCB, 
     &        POSDIAG,
     &        NRHS_B, NPIV,
     &        BLR_PANEL, LAST_BLR,
     &        CURRENT_BLR, BEGS_BLR_STATIC,
     &        KEEP8, K34, K148, K365, K405, K448, K450, K451, K455, 
     &        IS_T2_SLAVE, IFLAG, IERROR )
C
C  -------
C  Purpose
C  -------
C     
C     BLR update (gemm) of ARRAYPIB/CB for all blocks of the current panel.
C     On output, number of columns of current panel 
C
C     -----------
C     Parameters
C     -----------
C
C     ARRAYPIV is the array containing the FS (eliminated) rows of the solution (WCB
C     or RHSINTR depending on K350 in the future). It is of size LPIV
C     and leading dimension LDPIV. POSPIV points to the first eliminated row
C     associated to the current front in ARRAYPIV. 
C     ARRAYCB is the array containing the NFS+NELIM rows of the solution
C     (always WCB). It is of size LCB and leading dimension LDCB. 
C     POSCB points to the first non-eliminated row associated to the current front
C     in ARRAYCB. 
C     POSDIAG points to the beginning of the current panel CURRENT_BLR
C     (in ARRAYPIV).        
C     NPIV is the global NPIV of the front, used for masters only
C     LPIVCOL is either 1 (ARRAYPIV=WCB) or NRHS_B (ARRAYPIV=RHSINTR)        
        implicit none
      INTEGER(8), INTENT(IN) :: LPIV, LCB, POSPIV, POSCB, POSDIAG
      INTEGER, INTENT(IN)    :: LPIVCOL, POSPIVCOL
      REAL, TARGET, INTENT(INOUT) :: ARRAYPIV(LPIV,LPIVCOL)
      REAL, TARGET, INTENT(INOUT) :: ARRAYCB(LCB)
      INTEGER, INTENT(IN)    :: LAST_BLR, NRHS_B,  LDPIV, LDCB,
     &                          CURRENT_BLR, NPIV
      INTEGER, INTENT(IN)    :: K34, K148, K365, K405
      INTEGER, INTENT(IN)    :: K448, K450, K451, K455
C     TYPE(LRB_TYPE), TARGET,INTENT(IN) ::
C    &                           BLR_PANEL(LAST_BLR-CURRENT_BLR)
C     Assumed shape array declaration preferred: contiguity info
C     lost (+size lost in declaration), but no copy-in/copy-out:
      TYPE(LRB_TYPE), TARGET,INTENT(IN) ::
     &                           BLR_PANEL(:)
      INTEGER :: BEGS_BLR_STATIC(:)
      LOGICAL, INTENT(IN) :: IS_T2_SLAVE
      INTEGER(8), INTENT(IN) :: KEEP8(150)
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
C     -----------
C     Local variables
C     -----------
      INTEGER :: I, J1, J2, K, M, N, IBEG_BLOCK, IEND_BLOCK
      INTEGER :: KMAX, WSMAX
      INTEGER(8) :: POSBLOCK
      INTEGER :: PREC, allocok
      TYPE(LRB_TYPE), POINTER :: LRB
      REAL, ALLOCATABLE,DIMENSION(:) :: TEMP_BLOCK
      integer(c_intptr_t) :: th_ws_int
      type(c_ptr) :: th_ws_ptr
      REAL, POINTER :: thread_workspaces(:) 
      integer(c_int), allocatable, target, dimension(:) :: val_workspace
      type(c_ptr) :: val_ptr
      integer :: kcumul, kmagic, krem, kint, kprec
      integer :: cache_L2, start_prec
      REAL, pointer :: Q_ptr(:,:), R_ptr(:,:), Rhs_ptr(:)
      integer(kind=1), pointer :: stoc_ptr(:), a_ptr(:), b_ptr(:)
      integer :: ptr2_shape(2), ptr1_shape(1)
      REAL :: ONE, MONE, ZERO
      ! for cache stuff 
      DOUBLE PRECISION :: total_consumption, current_consumption
      DOUBLE PRECISION :: overall_consumption
      integer :: budget, budget_ws
      integer :: ws_footprint, prec_footprint
      integer :: start_k, current_k, current_prec, used_k
      real :: nu
      integer :: overall_k
      integer :: end_k, end_prec
      PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
      PARAMETER (ZERO=0.0E0)
      INTEGER, ALLOCATABLE, TARGET, DIMENSION(:) :: TEMP_RHS, 
     &  TEMP_RB, TEMP_QRB
      INTEGER :: POS_B, POS_QRB
      TYPE(LRMIXED_TYPE), ALLOCATABLE :: PT_RHS(:), PT_QRB(:)
      INTEGER, POINTER:: PT_RB(:)
      INTEGER :: MMAX
      INTEGER :: NB_INT_B, NB_INT_RB, NB_INT_QRB
      INTEGER :: LD_rhs, LD_R
      INTEGER :: PREC_CALC 
      LOGICAL :: A_CHEVAL
      LOGICAL :: MIXED_STORED
      double precision :: alpha
      integer :: th_num
      CACHE_L2 = K365
      MIXED_STORED = K450.GT.0.OR.K148.EQ.2
      KMAX = -1
      WSMAX = -1
      MMAX = -1
      DO I = CURRENT_BLR+1, LAST_BLR
        IBEG_BLOCK = BEGS_BLR_STATIC(I)
        IEND_BLOCK = BEGS_BLR_STATIC(I+1)-1
        IF (IBEG_BLOCK .EQ. IEND_BLOCK + 1) CYCLE
        LRB => BLR_PANEL(I-CURRENT_BLR)
        KMAX = max(KMAX, LRB%K)
        MMAX = max(MMAX, LRB%M)
        A_CHEVAL =  (.NOT.IS_T2_SLAVE) .AND. (IBEG_BLOCK.LE.NPIV) 
     &        .AND. (IEND_BLOCK.GT.NPIV)
        IF (.NOT.MIXED_STORED) CYCLE
C        CALL COMPUTE_WORKSPACE_SIZE(LRB, WSMAX,
C     &    K34, K148, K448, K450, K455, CACHE_L2,NRHS_B)
      ENDDO
C     Initializing shared variables:
      IF(K448.GT.0.AND.MIXED_STORED.AND.CURRENT_BLR.LT.LAST_BLR) THEN
          N = BLR_PANEL(1)%N
          NB_INT_B = 0
          DO PREC=2,NPREC_CALC
            NB_INT_B = NB_INT_B +
     &        (N*NRHS_B*COSTS_PREC(CALC2STOC(PREC)) +K34-1) / K34
          ENDDO
          allocate(TEMP_RHS(NB_INT_B), stat=allocok)
          if (allocok.gt.0) then
            IFLAG = -13
            IERROR = NB_INT_B
            write(*,*) 'Allocation problem in BLR routine 
     &      SMUMPS_SOL_FWD_BLR_UPDATE ',
     &      'not enough memory? memory requested = ', IERROR
            CALL MUMPS_ABORT()
          endif
          allocate(PT_RHS(NPREC_CALC), stat=allocok)
          if (allocok.gt.0) then
            IFLAG = -13
            IERROR = NPREC_CALC
            write(*,*) 'Allocation problem in BLR routine 
     &      SMUMPS_SOL_FWD_BLR_UPDATE ',
     &      'not enough memory? memory requested = ', IERROR
            CALL MUMPS_ABORT()
          endif
C         Nullify unused pointers (PREC=1)
C         Cast RHS in each precision in TEMP_RHS 
          POS_B=1
          DO PREC=2,NPREC_CALC
            PT_RHS(PREC)%val => TEMP_RHS(POS_B:)
            POS_B = POS_B+
     &        (N*NRHS_B*COSTS_PREC(CALC2STOC(PREC)) +K34-1) / K34
            CALL SMUMPS_LRMIXED_CAST(ARRAYPIV(POSDIAG,POSPIVCOL),
     &        PT_RHS(PREC)%val(1), N, 
     &        NRHS_B, LDPIV, N, 1, CALC2STOC(PREC))
          ENDDO
      ENDIF
      WSMAX = CACHE_L2 / 8
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL PRIVATE(TEMP_BLOCK, 
!$OMP& th_num, th_ws_ptr, th_ws_int,
!$OMP& thread_workspaces, Q_ptr, R_ptr, stoc_ptr,
!$OMP& ptr1_shape, ptr2_shape,
!$OMP& TEMP_QRB, TEMP_RB,
!$OMP& j1, j2,
!$OMP& PT_QRB, PT_RB, POS_QRB, PREC, NB_INT_RB, NB_INT_QRB,
!$OMP& alpha, ld_rhs, LD_R, rhs_ptr, val_ptr, val_workspace,
!$OMP& a_ptr, b_ptr,
!$OMP& allocok
!$OMP& )
#endif 
      th_num = 0
!$    th_num = omp_get_thread_num()
      if (K405.EQ.1) then
        th_num = thread_id - 1
      endif
      IF (K448.GT.0.AND.MIXED_STORED) THEN
C       Allocate RB and QRB workspaces for each thread        
        IF(CURRENT_BLR.LT.LAST_BLR) THEN
          NB_INT_QRB = 0
          DO PREC=2,NPREC_CALC
            NB_INT_QRB = NB_INT_QRB + 
     &         (MMAX*NRHS_B*COSTS_PREC(CALC2STOC(PREC)) +K34-1) / K34
          ENDDO
          NB_INT_RB=(KMAX*NRHS_B*COSTS_PREC(CALC2STOC(1)) + K34-1) / K34
          allocate(TEMP_RB(NB_INT_RB), stat=allocok)
          if (allocok.gt.0) then
!$OMP ATOMIC WRITE
            IFLAG = -13
!$OMP END ATOMIC
            IERROR = NB_INT_RB
            write(*,*) 'Allocation problem in BLR routine 
     &      SMUMPS_SOL_FWD_BLR_UPDATE over temp_rb',
     &      'not enough memory? memory requested = ', IERROR
          endif
          allocate(TEMP_QRB(NB_INT_QRB), stat=allocok)
          if (allocok.gt.0) then
!$OMP ATOMIC WRITE
            IFLAG = -13
!$OMP END ATOMIC
            IERROR = NB_INT_QRB
            write(*,*) 'Allocation problem in BLR routine 
     &      SMUMPS_SOL_FWD_BLR_UPDATE over temp_qrb',
     &      'not enough memory? memory requested = ', IERROR
          endif
          allocate(PT_QRB(NPREC_CALC), stat=allocok)
          if (allocok.gt.0) then
!$OMP ATOMIC WRITE
            IFLAG = -13
!$OMP END ATOMIC
            IERROR = NPREC_CALC
            write(*,*) 'Allocation problem in BLR routine 
     &      SMUMPS_SOL_FWD_BLR_UPDATE over pt_qrb',
     &      'not enough memory? memory requested = ', IERROR,
     &      'stat value', allocok
          endif
          nullify(PT_QRB(1)%val)
          POS_QRB = 1
          DO PREC=2,NPREC_CALC
            PT_QRB(PREC)%VAL => TEMP_QRB(POS_QRB:)
            POS_QRB = POS_QRB + 
     &         (MMAX*NRHS_B*COSTS_PREC(CALC2STOC(PREC)) +K34-1) / K34
          ENDDO
          PT_RB=>TEMP_RB(1:NB_INT_RB)
        ENDIF
      ENDIF
      IF (MIXED_STORED) THEN
        call mumps_sol_omp_mem_manager_get_address(th_num, th_ws_int)
        th_ws_ptr = transfer(th_ws_int, th_ws_ptr)
        call c_f_pointer(th_ws_ptr, thread_workspaces, (/ WSMAX /))
      ENDIF
      IF (KMAX.GT.0) THEN
C       TEMP_BLOCK is allocated even when K450.GT.0 because we might
C       need it when A_CHEVAL is true        
        allocate(TEMP_BLOCK(KMAX*NRHS_B), stat=allocok )
        IF (allocok .GT. 0) THEN
!$OMP ATOMIC WRITE
          IFLAG  = -13
!$OMP END ATOMIC
          IERROR = NRHS_B * KMAX
          write(*,*) 'Allocation problem in BLR routine 
     &    SMUMPS_SOL_FWD_BLR_UPDATE over temp_block: ',
     &    'not enough memory? memory requested = ', IERROR
        ENDIF
      ENDIF
! TODO try to remove this barrier
#if ! defined(BLR_NOOPENMP)
!$OMP DO SCHEDULE(DYNAMIC,1)
!$OMP& PRIVATE(IBEG_BLOCK, IEND_BLOCK, LRB, K, M, N,
!$OMP&         start_prec, end_prec,
!$OMP&         prec, kcumul, krem, kint, kmagic, kprec,
!$OMP&         current_consumption, current_k, current_prec,
!$OMP&         overall_consumption, total_consumption, nu,
!$OMP&         budget, budget_ws,
!$OMP&         overall_k, ws_footprint, prec_footprint,
!$OMP&         used_k, start_k, end_k,
!$OMP&         prec_calc, A_CHEVAL,
!$OMP&         POSBLOCK)
#endif
        DO I = CURRENT_BLR+1, LAST_BLR
          IF (IFLAG.LT.0) CYCLE
          IBEG_BLOCK = BEGS_BLR_STATIC(I)
          IEND_BLOCK = BEGS_BLR_STATIC(I+1)-1
          IF (IBEG_BLOCK .EQ. IEND_BLOCK + 1) CYCLE
          LRB => BLR_PANEL(I-CURRENT_BLR)
          K = LRB%K 
          M = LRB%M
          N = LRB%N 
C          write(*,*) "[DEV] SZ FWD BLR", M, N, K
          IF (LRB%ISLR) THEN
            IF (K.LE.0) CYCLE
C           Check if block is "a cheval": if so, do not use LRMIXED
C           for simplicity (we could try to code this in the future)
            A_CHEVAL =  (.NOT.IS_T2_SLAVE) .AND. (IBEG_BLOCK.LE.NPIV) 
     &            .AND. (IEND_BLOCK.GT.NPIV)
            IF (K448.EQ.0.OR..NOT.MIXED_STORED.OR. A_CHEVAL) THEN
              ! TODO can we allocate/deallocate workspaces just before
              ! computing ?
              ! we could but that feels stupid and overrelient on the
              ! operating system ; we can allocate too big workspaces
              ! beforehand and just use a bytes in those though
              kprec = lrb%k - lrb%k_stoc(1)
              if (.NOT.MIXED_STORED) then
                CALL sgemm('N', 'N', K, NRHS_B, N, ONE,
     &                LRB%R(1,1), K, ARRAYPIV(POSDIAG,POSPIVCOL), 
     &                LDPIV, ZERO, TEMP_BLOCK(1), K) 
              else if (k455.ge.3) then
                kcumul = 1
                if (lrb%k_stoc(1).gt.0) then
                  ! handle working arithmetic 
                  kprec = lrb%k_stoc(1)
                  ptr2_shape(1) = Kprec
                  ptr2_shape(2) = LRB%N
                  call c_f_pointer(c_loc(LRB%R_stoc(1)%val),
     &                             R_ptr, ptr2_shape)
                  CALL sgemm('N', 'N', Kprec, NRHS_B, N, ONE,
     &                  R_ptr(1,1), Kprec, 
     &                  ARRAYPIV(POSDIAG,POSPIVCOL), 
     &                  LDPIV, ZERO, TEMP_BLOCK(kcumul), K)
                end if 
                kprec = lrb%k - lrb%k_stoc(1)
                kcumul = kcumul + lrb%k_stoc(1)
                overall_k = lrb%k_stoc(1)
                total_consumption = lrb%n*( kprec * costs_prec(1) + 
     &                                  lrb%kw_stoc ) / 
     &                              ((4.0/5.0)*dble(cache_L2))
                start_prec = 2
                budget = 0
                budget_ws = 0
                start_k = 1
                current_k = 0
                overall_k = lrb%k_stoc(1)
                nu = 1.0d0
                if (k455.eq.4) nu = 0.2d0
                if (k455.eq.5.and.lrb%N.gt.200) nu = 0.2d0
                current_consumption = 0.0d0
                overall_consumption = 0.0d0
                do prec = 2, nprec_stoc
                  if (lrb%k_stoc(prec).eq.0) cycle
                  krem = lrb%k_stoc(prec)
                  current_prec = 1
                  do
                    ! avoid small blas when cache is barely over filled
                    if (k455.eq.3.and.nu.eq.1.0d0
     &          .and. total_consumption - overall_consumption.gt.1.0D0
     &          .and. total_consumption - overall_consumption.lt.1.5D0
     &              ) then
                      nu = (total_consumption - overall_consumption)/2.0
                      nu = nu + 0.02
                    end if
                    kmagic = max(dble(nu*(4.0/5.0)*CACHE_L2
     &                        - budget - budget_ws),0.0d0) /
     &                       (lrb%N*(costs_prec(prec)+costs_prec(1)))
C                    if (kmagic.lt.krem.and.kmagic.gt.8)
C     &                kmagic = kmagic - mod(kmagic,8)
                    used_k = min(kmagic, krem)
                    ws_footprint = lrb%n * used_k * costs_prec(1)
                    prec_footprint = lrb%n * used_k * costs_prec(prec)
                    budget = budget + prec_footprint
                    budget_ws = budget_ws + ws_footprint
                    current_consumption = current_consumption +
     &                dble(prec_footprint+ws_footprint) / 
     &                ((80.0/100.0)*dble(cache_L2))
                    overall_consumption = overall_consumption +
     &                dble(prec_footprint+ws_footprint) / 
     &                ((80.0/100.0)*dble(cache_L2))
                    current_prec = current_prec + used_k
                    current_k = current_k + used_k
                    krem = krem - used_k
                    overall_k = overall_k + used_k
                    if (krem.ge.1.and.k455.eq.4
     &                  .and.current_consumption.gt.0.3.or.
     &                  krem.ge.1.or.prec.eq.nprec_stoc.or.
     &                  overall_k.eq.k) then
                      R_ptr(1:current_k,1:N) => 
     &                    thread_workspaces(1:(N*current_k))
                      end_prec = prec
                      end_k = lrb%k_stoc(prec) - krem
                      ! if we switched prec but even a row wont fit
                      if (krem.eq.lrb%k_stoc(prec)) then
                        end_prec = prec - 1
                        end_k = lrb%k_stoc(end_prec)
                      end if
                      CALL SMUMPS_DECONVERT_LRB_R_W(LRB, R_ptr, 
     &                    current_k, N,
     &                    start_prec, end_prec, 
     &                    start_k, end_k,
     &                    IFLAG, IERROR, KEEP8)
C                      write(*,*) "MYSTUFF", n, current_k, k
                      CALL sgemm('N', 'N', current_k, NRHS_B, N, ONE,
     &                    R_ptr(1,1), current_k, 
     &                    ARRAYPIV(POSDIAG,POSPIVCOL),  LDPIV, 
     &                    ZERO, TEMP_BLOCK(kcumul), K)
                      kcumul = kcumul + current_k
                      start_prec = prec
                      budget = 0
                      budget_ws = 0
                      start_k = current_prec
                      current_k = 0
                      current_consumption = 0.0d0
                    end if
                    if (krem.lt.1) then 
                      exit
                    end if
                  end do
                end do
              else
                kcumul = lrb%k + 1
                start_prec=nprec_stoc
                ! if the sum of R^T (except working prec.) fit in cache
                if (kprec.gt.0.and.
     &              (k455.eq.0.or.k455.ne.1.and.
     &              N*LRB%kw_stoc +
     &              N*costs_prec(1) * kprec < (4.0/5)*cache_L2)) then
                  kcumul = kcumul - kprec
                  R_ptr(1:kprec,1:N) => 
     &                thread_workspaces(1:(N*Kprec))
                  CALL SMUMPS_DECONVERT_LRB_R_W(LRB, 
     &                R_ptr, kprec, N,
     &                2, NPREC_STOC, 1, LRB%K_STOC(NPREC_STOC),
     &                IFLAG, IERROR, KEEP8)
C                  write(*,*) "MYSTUFF", n, kprec, k
                  CALL sgemm('N', 'N', Kprec, NRHS_B, N, ONE,
     &                R_ptr(1,1), Kprec, 
     &                ARRAYPIV(POSDIAG,POSPIVCOL), 
     &                LDPIV, ZERO, TEMP_BLOCK(kcumul), K) 
                  ! working prec is adressed in following do loop
                  start_prec=1
                endif
                ! else (e.g. if any R^T_prec is too big)
                do prec=start_prec,1,-1
                   kprec = lrb%k_stoc(prec)
                   if (kprec.eq.0) cycle
                   kcumul = kcumul - kprec
                   if (prec.eq.1) then
                     KMAGIC=Kprec
                     ptr2_shape(1) = Kprec
                     ptr2_shape(2) = LRB%N
                     call c_f_pointer(c_loc(LRB%R_stoc(prec)%val),
     &                                R_ptr, ptr2_shape)
                   else
                     ptr1_shape(1) = costs_prec(prec)*Kprec*LRB%N
                     call c_f_pointer(c_loc(LRB%R_stoc(prec)%val),
     &                  stoc_ptr, ptr1_shape)
                     KMAGIC=(80.0/100.0)*(CACHE_L2-costs_prec(1)*N)/
     &                 ( (costs_prec(1)+costs_prec(prec))*N + 
     &                 NRHS_B*costs_prec(1) )
                     KMAGIC = KMAGIC - MODULO(KMAGIC,8)
                   endif
                   if (k455.eq.1) kmagic = kprec
                   do kint=1,Kprec,KMAGIC
                     !Krem may be stupidly small
                     krem = min(KMAGIC,Kprec-Kint+1)
C                     write(*,*) "[MYSTATS] R", 
C     &                  N, Krem, K, KPREC, COSTS_PREC(PREC)
                     if (prec.gt.1) then
                       R_ptr(1:Krem,1:LRB%N) => 
     &                     thread_workspaces(1:(LRB%N*Krem))
C                        write(*,*) "MYSTUFF", n, krem, k
                       CALL SMUMPS_LRMIXED_CAST(
     &                     stoc_ptr(
     &                       costs_prec(prec)*(kint-1)+1),
     &                     R_ptr(1,1), Krem, LRB%N, 
     &                     Kprec, Krem, PREC, 1)
                     endif
                     CALL sgemm('N', 'N', Krem, NRHS_B, N, ONE,
     &                     R_ptr(1,1), Krem, 
     &                     ARRAYPIV(POSDIAG,POSPIVCOL), 
     &                     LDPIV, ZERO, TEMP_BLOCK(kcumul+kint-1), K) 
                   end do
                end do
              end if
              if (.not.MIXED_STORED) then
                CALL GEMM_Q_FWD(M, NRHS_B, K, NPIV, 
     &              LRB%Q, TEMP_BLOCK, K, 
     &              ARRAYPIV, ldpiv, ARRAYCB, lcb, ldcb, 
     &              ibeg_block, iend_block, is_t2_slave, 
     &              poscb, pospiv, pospivcol, 1)
              ! FIXME the splitting of updates here should
              ! be more abstract
              else if (k455.ge.3) then
                kcumul = 1
                if (lrb%k_stoc(1).gt.0) then
                  ! handle working arithmetic 
                  kprec = lrb%k_stoc(1)
                  ptr2_shape(1) = M
                  ptr2_shape(2) = Kprec
                  call c_f_pointer(c_loc(LRB%Q_stoc(1)%val),
     &                             Q_ptr, ptr2_shape)
                  CALL GEMM_Q_FWD(M, NRHS_B, Kprec, NPIV, 
     &                Q_ptr, TEMP_BLOCK, K, 
     &                ARRAYPIV, ldpiv, ARRAYCB, lcb, ldcb, 
     &                ibeg_block, iend_block, is_t2_slave, 
     &                poscb, pospiv, pospivcol, kcumul)
                end if 
                overall_k = lrb%k_stoc(1)
                kprec = lrb%k - lrb%k_stoc(1)
                kcumul = kcumul + lrb%k_stoc(1)
                total_consumption = m*( kprec * costs_prec(1) + 
     &                                  lrb%kw_stoc ) / 
     &                              ((4.0/5.0)*dble(cache_L2))
                start_prec = 2
                budget = 0
                budget_ws = 0
                start_k = 1
                current_k = 0
                nu = 1.0d0
                if (k455.eq.4) nu = 0.2d0
                if (k455.eq.5.and.lrb%M.gt.200) nu = 0.2d0
                current_consumption = 0.0d0
                overall_consumption = 0.0d0
                overall_k = lrb%k_stoc(1)
                do prec = 2, nprec_stoc
                  if (lrb%k_stoc(prec).eq.0) cycle
                  krem = lrb%k_stoc(prec)
                  current_prec = 1
                  do
                    ! avoid small blas when cache is barely over filled
                    if (k455.eq.3.and.nu.eq.1.0d0
     &          .and. total_consumption - overall_consumption.gt.1.0D0
     &          .and. total_consumption - overall_consumption.lt.1.5D0
     &              ) then
                      nu = (total_consumption - overall_consumption)/2.0
                      nu = nu + 0.02
                    end if
                    kmagic = max(dble(nu*(4.0/5.0)*CACHE_L2
     &                        - budget - budget_ws),0.0d0) /
     &                       (lrb%M*(costs_prec(prec)+costs_prec(1)))
                    used_k = min(kmagic, krem)
                    ws_footprint = lrb%M * used_k * costs_prec(1)
                    prec_footprint = lrb%M * used_k * costs_prec(prec)
                    budget = budget + prec_footprint
                    budget_ws = budget_ws + ws_footprint
                    current_consumption = current_consumption +
     &                dble(prec_footprint+ws_footprint) / 
     &                ((80.0/100.0)*dble(cache_L2))
                    overall_consumption = overall_consumption +
     &                dble(prec_footprint+ws_footprint) / 
     &                ((80.0/100.0)*dble(cache_L2))
                    current_prec = current_prec + used_k
                    current_k = current_k + used_k
                    krem = krem - used_k
                    overall_k = overall_k + used_k
                    if (krem.ge.1.and.k455.eq.4
     &                  .and.current_consumption.gt.0.3.or.
     &                  krem.ge.1.or.prec.eq.nprec_stoc.or.
     &                  overall_k.eq.k) then
                      end_prec = prec
                      end_k = lrb%k_stoc(prec) - krem
                      ! if we switched prec but even a row wont fit
                      if (krem.eq.lrb%k_stoc(prec)) then
                        end_prec = prec - 1
                        end_k = lrb%k_stoc(end_prec)
                      end if
                      Q_ptr(1:M,1:current_k) => 
     &                    thread_workspaces(1:(M*current_k))
                      CALL SMUMPS_DECONVERT_LRB_Q_W(LRB, Q_ptr, 
     &                    M, current_k,
     &                    start_prec, end_prec, 
     &                    start_k, end_k,
     &                    IFLAG, IERROR, KEEP8)
C                      write(*,*) "MYSTUFF", m, current_k, k
                      CALL GEMM_Q_FWD(M, NRHS_B, current_K, NPIV, 
     &                    Q_ptr, TEMP_BLOCK, K, 
     &                    ARRAYPIV, ldpiv, ARRAYCB, lcb, ldcb, 
     &                    ibeg_block, iend_block, is_t2_slave, 
     &                    poscb, pospiv, pospivcol, kcumul)
                      kcumul = kcumul + current_k
                      start_prec = prec
                      budget = 0
                      budget_ws = 0
                      start_k = current_prec
                      current_k = 0
                      current_consumption = 0.0d0
                    end if
                    if (krem.lt.1) then 
                      exit
                    end if
                  end do
                end do
              else
                kcumul = lrb%k + 1
                start_prec=nprec_stoc
                ! if the sum of R^T (except working prec.) fit in cache
                kprec = lrb%k - lrb%k_stoc(1)
                if (kprec.gt.0.and.
     &              (k455.eq.0.or.k455.ne.1.and.
     &              M*lrb%kw_stoc +
     &              M*costs_prec(1) * kprec < (4.0/5)*cache_L2)) then
                  kcumul = kcumul - kprec
C                  write(*,*) "[MYSTATS] Q", 
C     &                M, kprec, k, -1, -1
                  Q_ptr(1:M,1:kprec) =>
     &                thread_workspaces(1:(M*Kprec))
                  CALL SMUMPS_DECONVERT_LRB_Q_W(LRB, 
     &                Q_ptr, M, kprec,
     &                2, NPREC_STOC, 1, LRB%K_STOC(NPREC_STOC),
     &                IFLAG, IERROR, KEEP8)
C                  write(*,*) "MYSTUFF", m, kprec, k
                  CALL GEMM_Q_FWD(M, NRHS_B, Kprec, NPIV, 
     &                Q_ptr, TEMP_BLOCK, K, 
     &                ARRAYPIV, ldpiv, ARRAYCB, lcb, ldcb, 
     &                ibeg_block, iend_block, is_t2_slave, 
     &                poscb, pospiv, pospivcol, kcumul)
                  start_prec=1
                endif
                ! else (e.g. if any R^T_prec is too big)
                do prec=start_prec,1,-1
                  kprec = lrb%k_stoc(prec)
                  if (kprec.eq.0) cycle
                  kcumul = kcumul - kprec
                  if (prec.eq.1) then
                    kmagic = kprec
                    ptr2_shape(1) = M
                    ptr2_shape(2) = Kprec
                    call c_f_pointer(c_loc(LRB%Q_stoc(prec)%val),
     &                               Q_ptr, ptr2_shape)
                  else
                    ptr1_shape(1) = costs_prec(prec)*LRB%M*Kprec
                    call c_f_pointer(c_loc(LRB%Q_stoc(prec)%val),
     &                    stoc_ptr, ptr1_shape)
                    KMAGIC=(80.0/100.0)*(CACHE_L2-costs_prec(1)*M)/
     &              ( (costs_prec(1)+costs_prec(prec))*M + 
     &                NRHS_B*costs_prec(1) )
                    KMAGIC = KMAGIC - MODULO(KMAGIC,8)
                  endif
                  if (k455.eq.1) kmagic = kprec
                  do kint=1,Kprec,KMAGIC
                    !Krem may be stupidly small
                    krem = min(KMAGIC,Kprec-Kint+1)
C                    write(*,*) "[MYSTATS] Q", 
C     &                  M, Krem, K, KPREC, COSTS_PREC(PREC)
                    if (prec.gt.1) then
                      Q_ptr(1:M,1:Krem) =>
     &                  thread_workspaces(1:(M*krem))
C                      write(*,*) "MYSTUFF", m, krem, k
                      CALL SMUMPS_LRMIXED_CAST(
     &                  stoc_ptr(m*costs_prec(prec)*
     &                             (kint-1)+1),
     &                  Q_ptr(1,1), M, Krem,
     &                  M, M, PREC, 1)
                    endif
                    CALL GEMM_Q_FWD(M, NRHS_B, Krem, NPIV, 
     &                  Q_ptr, TEMP_BLOCK, K,
     &                  ARRAYPIV, ldpiv, ARRAYCB, lcb, ldcb, 
     &                  ibeg_block, iend_block, is_t2_slave,
     &                  poscb, pospiv, pospivcol, kcumul+kint-1)
                  end do
                end do
              end if
            ELSE
              IF (IS_T2_SLAVE) THEN
                POSBLOCK = POSCB+int(IBEG_BLOCK-1,8)
                CALL SMUMPS_LRMIXED_LRxFR(0, LRB,
     &            ARRAYPIV(POSDIAG,POSPIVCOL), LDPIV,
     &            ARRAYCB(POSBLOCK), LDCB,
     &            NRHS_B, -1.0D0,
     &            K34, K451, KEEP8, PT_RHS, PT_RB(1), PT_QRB,
     &            thread_workspaces, WSMAX)
              ELSEIF (IBEG_BLOCK.LE.NPIV) THEN
C               Entirely on FS side              
                POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
                CALL SMUMPS_LRMIXED_LRxFR(0, LRB,
     &            ARRAYPIV(POSDIAG,POSPIVCOL), LDPIV,
     &            ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV,
     &            NRHS_B, -1.0D0, 
     &            K34, K451, KEEP8, PT_RHS, PT_RB(1), PT_QRB,
     &            thread_workspaces, WSMAX)
              ELSE
C               Entirely on CB side              
                POSBLOCK = POSCB+int(IBEG_BLOCK-1-NPIV,8)
                CALL SMUMPS_LRMIXED_LRxFR(0, LRB,
     &            ARRAYPIV(POSDIAG,POSPIVCOL), LDPIV,
     &            ARRAYCB(POSBLOCK), LDCB, 
     &            NRHS_B, -1.0D0,
     &            K34, K451, KEEP8, PT_RHS, PT_RB(1), PT_QRB,
     &            thread_workspaces, WSMAX)
              ENDIF
            ENDIF
! LRMIXED
          ELSE ! NOT ISLR
            CALL FRB_GET_PREC(LRB, PREC, MIXED_STORED)
            IF (K448.GE.1.AND.K450.NE.0) THEN
              IF (IS_T2_SLAVE) THEN
                POSBLOCK = POSCB + int(IBEG_BLOCK-1,8)
                CALL SMUMPS_LRMIXED_FRxFR('N', 'N',
     &                M, NRHS_B, N, MONE, 
     &                LRB, 1, 1, M,
     &                ARRAYPIV(POSDIAG,POSPIVCOL),LDPIV,
     &                ARRAYCB(POSBLOCK), LDCB,
     &                PT_RHS, TEMP_QRB, thread_workspaces, WSMAX)
              ELSEIF (IBEG_BLOCK.LE.NPIV.AND.IEND_BLOCK.GT.NPIV) THEN
C               Case of overlap between the FS and non FS
                POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
                CALL SMUMPS_LRMIXED_FRxFR('N', 'N',
     &                NPIV-IBEG_BLOCK+1, NRHS_B, N, MONE, 
     &                LRB, 1, 1, M,
     &                ARRAYPIV(POSDIAG,POSPIVCOL),LDPIV,
     &                ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV,
     &                PT_RHS, TEMP_QRB, thread_workspaces, WSMAX)
                CALL SMUMPS_LRMIXED_FRxFR('N', 'N',
     &                M-(NPIV-IBEG_BLOCK+1), NRHS_B, N, MONE, 
     &                LRB, NPIV-IBEG_BLOCK+2, 1, M,
     &                ARRAYPIV(POSDIAG,POSPIVCOL),LDPIV,
     &                ARRAYCB(POSCB), LDCB,
     &                PT_RHS, TEMP_QRB, thread_workspaces, WSMAX)
              ELSEIF (IBEG_BLOCK.LE.NPIV) THEN
C               Entirely on FS side              
                POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
                CALL SMUMPS_LRMIXED_FRxFR('N', 'N',
     &                M, NRHS_B, N, MONE, 
     &                LRB, 1, 1, M,
     &                ARRAYPIV(POSDIAG,POSPIVCOL), LDPIV,
     &                ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV,
     &                PT_RHS, TEMP_QRB, thread_workspaces, WSMAX)
              ELSE
C               Entirely on CB side              
                POSBLOCK = POSCB + int(IBEG_BLOCK-1-NPIV,8)
                CALL SMUMPS_LRMIXED_FRxFR('N', 'N',
     &                M, NRHS_B, N, MONE, 
     &                LRB, 1, 1, M,
     &                ARRAYPIV(POSDIAG,POSPIVCOL),LDPIV,
     &                ARRAYCB(POSBLOCK), LDCB,
     &                PT_RHS, TEMP_QRB, thread_workspaces, WSMAX)
              ENDIF
            ELSE
              if (prec.le.1.and.associated(LRB%Q)) then
                KMAGIC = N
                ptr2_shape(1) = M
                ptr2_shape(2) = KMAGIC
                call c_f_pointer(c_loc(LRB%Q),
     &                           Q_ptr, ptr2_shape)
              else
                KMAGIC = (80.0/100.0)*CACHE_L2/
     &            ((costs_prec(1)+costs_prec(prec))*M)
                KMAGIC=MIN(KMAGIC,N)
                Q_ptr(1:M,1:KMAGIC) =>
     &            thread_workspaces(1:(M*KMAGIC))
                ptr1_shape(1) = costs_prec(prec)*M*N
                call c_f_pointer(c_loc(LRB%Q_stoc(prec)%val),
     &            stoc_ptr, ptr1_shape)
              ENDIF
              DO KINT=1,N,KMAGIC
                krem = min(KMAGIC,N-Kint+1)
                if (prec.gt.1) then
                  CALL SMUMPS_LRMIXED_CAST(
     &              stoc_ptr(costs_prec(prec)*M*(kint-1)+1),
     &              Q_ptr(1,1), M, KREM,
     &              M, M, PREC, 1)
                ENDIF
! LRMIXED
                ! FIXME use GEMM_Q_FWD
                IF (IS_T2_SLAVE) THEN
                  POSBLOCK = POSCB + int(IBEG_BLOCK-1,8)
                  CALL sgemm('N', 'N', M, NRHS_B, KREM, MONE,
     &                  Q_ptr(1,1), M, 
     &                  ARRAYPIV(POSDIAG+KINT-1,POSPIVCOL), LDPIV,
     &                  ONE, ARRAYCB(POSBLOCK), LDCB) 
                ELSEIF (IBEG_BLOCK.LE.NPIV.AND.IEND_BLOCK.GT.NPIV) THEN
C                 Case of overlap between the FS and non FS
                  POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
                  CALL sgemm('N', 'N', NPIV-IBEG_BLOCK+1, NRHS_B, 
     &                  KREM, MONE,
     &                  Q_ptr(1,1), M,
     &                  ARRAYPIV(POSDIAG+KINT-1,POSPIVCOL), LDPIV,
     &                  ONE, ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV) 
                  CALL sgemm('N', 'N', IBEG_BLOCK+M-NPIV-1, NRHS_B,
     &                  KREM, MONE,
     &                  Q_ptr(NPIV-IBEG_BLOCK+2,1), M,
     &                  ARRAYPIV(POSDIAG+KINT-1,POSPIVCOL),
     &                  LDPIV, ONE, ARRAYCB(POSCB), LDCB) 
                ELSEIF (IBEG_BLOCK.LE.NPIV) THEN
C                 Entirely on FS side              
                  POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
                  CALL sgemm('N', 'N', M, NRHS_B, KREM, MONE,
     &                  Q_ptr(1,1), M,
     &                  ARRAYPIV(POSDIAG+KINT-1,POSPIVCOL),LDPIV,
     &                  ONE, ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV)
                ELSE
C                 Entirely on CB side              
                  POSBLOCK = POSCB + int(IBEG_BLOCK-1-NPIV,8)
                  CALL sgemm('N', 'N', M, NRHS_B, KREM, MONE,
     &                  Q_ptr(1,1), M,
     &                  ARRAYPIV(POSDIAG+KINT-1,POSPIVCOL),LDPIV,
     &                  ONE, ARRAYCB(POSBLOCK), LDCB) 
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ENDDO
#if ! defined(BLR_NOOPENMP)
!$OMP END DO
#endif
        IF (allocated(TEMP_BLOCK)) deallocate(TEMP_BLOCK)
        IF (allocated(PT_QRB))   deallocate(PT_QRB)
        IF (allocated(TEMP_RB))  deallocate(TEMP_RB)
        IF (allocated(TEMP_QRB)) deallocate(TEMP_QRB)
        IF (allocated(val_workspace)) deallocate(val_workspace)
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL         
#endif 
        IF (allocated(PT_RHS))   deallocate(PT_RHS)
        IF (allocated(TEMP_RHS)) deallocate(TEMP_RHS)
        RETURN
      END SUBROUTINE SMUMPS_SOL_FWD_BLR_UPDATE
! (MIXHPC)
#else 
! (.not.MIXHPC)
      SUBROUTINE SMUMPS_SOL_FWD_BLR_UPDATE (
     &        ARRAYPIV, LPIV, LPIVCOL, LDPIV, POSPIV, POSPIVCOL,
     &        ARRAYCB, LCB, LDCB, POSCB, 
     &        POSDIAG,
     &        NRHS_B, NPIV,
     &        BLR_PANEL, LAST_BLR,
     &        CURRENT_BLR, BEGS_BLR_STATIC,
     &        KEEP8, K34, K448, K450, K451, IS_T2_SLAVE, IFLAG, IERROR )
C
C  -------
C  Purpose
C  -------
C     
C     BLR update (gemm) of ARRAYPIB/CB for all blocks of the current panel.
C     On output, number of columns of current panel 
C
C     -----------
C     Parameters
C     -----------
C
C     ARRAYPIV is the array containing the FS (eliminated) rows of the solution (WCB
C     or RHSINTR depending on K350 in the future). It is of size LPIV
C     and leading dimension LDPIV. POSPIV points to the first eliminated row
C     associated to the current front in ARRAYPIV. 
C     ARRAYCB is the array containing the NFS+NELIM rows of the solution
C     (always WCB). It is of size LCB and leading dimension LDCB. 
C     POSCB points to the first non-eliminated row associated to the current front
C     in ARRAYCB. 
C     POSDIAG points to the beginning of the current panel CURRENT_BLR
C     (in ARRAYPIV).        
C     NPIV is the global NPIV of the front, used for masters only
C     LPIVCOL is either 1 (ARRAYPIV=WCB) or NRHS_B (ARRAYPIV=RHSINTR)        
      INTEGER(8), INTENT(IN) :: LPIV, LCB, POSPIV, POSCB, POSDIAG
      INTEGER, INTENT(IN)    :: LPIVCOL, POSPIVCOL
      REAL, INTENT(INOUT) :: ARRAYPIV(LPIV,LPIVCOL)
      REAL, INTENT(INOUT) :: ARRAYCB(LCB)
      INTEGER, INTENT(IN)    :: LAST_BLR, NRHS_B,  LDPIV, LDCB,
     &                          CURRENT_BLR, NPIV, K34, K448, K450, K451
C     TYPE(LRB_TYPE), TARGET,INTENT(IN) ::
C    &                           BLR_PANEL(LAST_BLR-CURRENT_BLR)
C     Assumed shape array declaration preferred: contiguity info
C     lost (+size lost in declaration), but no copy-in/copy-out:
      TYPE(LRB_TYPE), TARGET,INTENT(IN) ::
     &                           BLR_PANEL(:)
      INTEGER :: BEGS_BLR_STATIC(:)
      LOGICAL, INTENT(IN) :: IS_T2_SLAVE
      INTEGER(8), INTENT(IN) :: KEEP8(150)
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
C     -----------
C     Local variables
C     -----------
      INTEGER :: I, K, M, N, IBEG_BLOCK, IEND_BLOCK
      INTEGER :: KMAX
      INTEGER :: MMAX
      INTEGER(8) :: POSBLOCK
      INTEGER :: allocok
      TYPE(LRB_TYPE), POINTER :: LRB
      REAL, ALLOCATABLE,DIMENSION(:) :: TEMP_BLOCK
      REAL :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
      PARAMETER (ZERO=0.0E0)
      INTEGER, ALLOCATABLE, TARGET, DIMENSION(:) :: TEMP_RHS, 
     &  TEMP_RB, TEMP_QRB
      INTEGER :: POS_B, POS_QRB
      TYPE(LRMIXED_TYPE), ALLOCATABLE :: PT_RHS(:), PT_QRB(:)
      INTEGER, POINTER:: PT_RB(:)
      INTEGER :: PREC, NB_INT_B, NB_INT_RB, NB_INT_QRB
      LOGICAL :: A_CHEVAL
      KMAX = -1
      MMAX = -1
      DO I = CURRENT_BLR+1, LAST_BLR
        KMAX = max(KMAX, BLR_PANEL(I-CURRENT_BLR)%K)
        MMAX = max(MMAX, BLR_PANEL(I-CURRENT_BLR)%M)
      ENDDO
C     Initializing shared variables:
      IF(K448.GE.1.AND.K450.GT.0.AND. (CURRENT_BLR.LT.LAST_BLR)) THEN
          N = BLR_PANEL(1)%N
          NB_INT_B = 0
          DO PREC=2,NPREC_CALC
            NB_INT_B = NB_INT_B +
     &        (N*NRHS_B*COSTS_PREC(CALC2STOC(PREC)) +K34-1) / K34
          ENDDO
          allocate(TEMP_RHS(NB_INT_B), stat=allocok)
          if (allocok.gt.0) then
            IFLAG = -13
            IERROR = NB_INT_B
            write(*,*) 'Allocation problem in BLR routine 
     &      SMUMPS_SOL_FWD_BLR_UPDATE ',
     &      'not enough memory? memory requested = ', IERROR
            CALL MUMPS_ABORT()
          endif
          allocate(PT_RHS(NPREC_CALC), stat=allocok)
          if (allocok.gt.0) then
            IFLAG = -13
            IERROR = NPREC_CALC
            write(*,*) 'Allocation problem in BLR routine 
     &      SMUMPS_SOL_FWD_BLR_UPDATE ',
     &      'not enough memory? memory requested = ', IERROR
            CALL MUMPS_ABORT()
          endif
C         Cast RHS in each precision in TEMP_RHS 
          POS_B=1
          DO PREC=2,NPREC_CALC
            PT_RHS(PREC)%val => TEMP_RHS(POS_B:)
            POS_B = POS_B+
     &        (N*NRHS_B*COSTS_PREC(CALC2STOC(PREC)) +K34-1) / K34
            CALL SMUMPS_LRMIXED_CAST(ARRAYPIV(POSDIAG,POSPIVCOL),
     &        PT_RHS(PREC)%val(1), N, 
     &        NRHS_B, LDPIV, N, 1, CALC2STOC(PREC))
          ENDDO
      ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL PRIVATE(TEMP_BLOCK, allocok
!$OMP& , TEMP_QRB, TEMP_RB,
!$OMP& PT_QRB, PT_RB, POS_QRB, PREC, NB_INT_RB, NB_INT_QRB
!$OMP& )
#endif 
      IF (K448.GE.1.AND.K450.GT.0) THEN
C       Allocate RB and QRB workspaces for each thread        
        IF(CURRENT_BLR.LT.LAST_BLR) THEN
          NB_INT_QRB = 0
          DO PREC=2,NPREC_CALC
            NB_INT_QRB = NB_INT_QRB + 
     &         (MMAX*NRHS_B*COSTS_PREC(CALC2STOC(PREC)) +K34-1) / K34
          ENDDO
          NB_INT_RB=(KMAX*NRHS_B*COSTS_PREC(CALC2STOC(1)) + K34-1) / K34
          allocate(TEMP_RB(NB_INT_RB), stat=allocok)
          if (allocok.gt.0) then
!$OMP ATOMIC WRITE
            IFLAG = -13
!$OMP END ATOMIC
            IERROR = NB_INT_RB
            write(*,*) 'Allocation problem in BLR routine 
     &      SMUMPS_SOL_FWD_BLR_UPDATE over temp_rb',
     &      'not enough memory? memory requested = ', IERROR
          endif
          allocate(TEMP_QRB(NB_INT_QRB), stat=allocok)
          if (allocok.gt.0) then
!$OMP ATOMIC WRITE
            IFLAG = -13
!$OMP END ATOMIC
            IERROR = NB_INT_QRB
            write(*,*) 'Allocation problem in BLR routine 
     &      SMUMPS_SOL_FWD_BLR_UPDATE over temp_qrb',
     &      'not enough memory? memory requested = ', IERROR
          endif
          allocate(PT_QRB(NPREC_CALC), stat=allocok)
          if (allocok.gt.0) then
!$OMP ATOMIC WRITE
            IFLAG = -13
!$OMP END ATOMIC
            IERROR = NPREC_CALC
            write(*,*) 'Allocation problem in BLR routine 
     &      SMUMPS_SOL_FWD_BLR_UPDATE over pt_qrb',
     &      'not enough memory? memory requested = ', IERROR
          endif
          nullify(PT_QRB(1)%val)
          POS_QRB = 1
          DO PREC=2,NPREC_CALC
            PT_QRB(PREC)%VAL => TEMP_QRB(POS_QRB:)
            POS_QRB = POS_QRB + 
     &         (MMAX*NRHS_B*COSTS_PREC(CALC2STOC(PREC)) +K34-1) / K34
          ENDDO
          PT_RB=>TEMP_RB(1:NB_INT_RB)
        ENDIF
      ENDIF
      IF (KMAX.GT.0) THEN
C       TEMP_BLOCK is allocated even when K450.GT.0 because we might
C       need it when A_CHEVAL is true        
        allocate(TEMP_BLOCK(KMAX*NRHS_B), stat=allocok )
        IF (allocok .GT. 0) THEN
          IFLAG  = -13
          IERROR = NRHS_B * KMAX
          write(*,*) 'Allocation problem in BLR routine 
     &    SMUMPS_SOL_FWD_BLR_UPDATE for TEMP_BLOCK: ',
     &    'not enough memory? memory requested = ', IERROR
        ENDIF
      ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif 
#if ! defined(BLR_NOOPENMP)
!$OMP DO SCHEDULE(DYNAMIC,1)
!$OMP& PRIVATE(IBEG_BLOCK, IEND_BLOCK, LRB, K, M, N,
!$OMP&         A_CHEVAL,
!$OMP&         POSBLOCK)
#endif
        DO I = CURRENT_BLR+1, LAST_BLR
          IF (IFLAG.LT.0) CYCLE
          IBEG_BLOCK = BEGS_BLR_STATIC(I)
          IEND_BLOCK = BEGS_BLR_STATIC(I+1)-1
          IF (IBEG_BLOCK .EQ. IEND_BLOCK + 1) CYCLE
          LRB => BLR_PANEL(I-CURRENT_BLR)
          K = LRB%K 
          M = LRB%M
          N = LRB%N 
          IF (LRB%ISLR) THEN
            IF (K.GT.0) THEN
C             Check if block is "a cheval": if so, do not use LRMIXED
C             for simplicity (we could try to code this in the future)
              A_CHEVAL =  (.NOT.IS_T2_SLAVE) .AND. (IBEG_BLOCK.LE.NPIV) 
     &              .AND. (IEND_BLOCK.GT.NPIV)
              IF (K450.GT.0
     &            .AND.(K448.EQ.0.OR. A_CHEVAL)
     &           ) THEN
                CALL SMUMPS_DECONVERT_LRB(LRB, IFLAG, IERROR, KEEP8)
              ENDIF
              IF (K448.EQ.0.OR.K450.EQ.0.OR. A_CHEVAL) THEN
                CALL sgemm('N', 'N', K, NRHS_B, N, ONE,
     &                LRB%R(1,1), K, ARRAYPIV(POSDIAG,POSPIVCOL), LDPIV,
     &                ZERO, TEMP_BLOCK(1), K) 
                IF (IS_T2_SLAVE) THEN
                  POSBLOCK = POSCB+int(IBEG_BLOCK-1,8)
                  CALL sgemm('N', 'N', M, NRHS_B, K, MONE,
     &              LRB%Q(1,1), M, TEMP_BLOCK(1), K,
     &              ONE, ARRAYCB(POSBLOCK), LDCB) 
                ELSEIF (IBEG_BLOCK.LE.NPIV.AND.IEND_BLOCK.GT.NPIV) THEN
C                 Case of overlap between the FS and non FS
                  POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
                  CALL sgemm('N', 'N', NPIV-IBEG_BLOCK+1, NRHS_B, K,
     &              MONE, LRB%Q(1,1), M, TEMP_BLOCK(1), K,
     &              ONE, ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV) 
                  CALL sgemm('N', 'N', IBEG_BLOCK+M-NPIV-1, NRHS_B, K, 
     &              MONE, LRB%Q(NPIV-IBEG_BLOCK+2,1), M, TEMP_BLOCK(1),
     &              K, ONE, ARRAYCB(POSCB), LDCB) 
                ELSEIF (IBEG_BLOCK.LE.NPIV) THEN
C                 Entirely on FS side              
                  POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
                  CALL sgemm('N', 'N', M, NRHS_B, K, MONE,
     &              LRB%Q(1,1), M, TEMP_BLOCK(1), K,
     &              ONE, ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV) 
                ELSE
C                 Entirely on CB side              
                  POSBLOCK = POSCB+int(IBEG_BLOCK-1-NPIV,8)
                  CALL sgemm('N', 'N', M, NRHS_B, K, MONE,
     &              LRB%Q(1,1), M, TEMP_BLOCK(1), K,
     &              ONE, ARRAYCB(POSBLOCK), LDCB) 
                ENDIF
              ELSE
                IF (IS_T2_SLAVE) THEN
                  POSBLOCK = POSCB+int(IBEG_BLOCK-1,8)
                  CALL SMUMPS_LRMIXED_LRxFR(0, LRB,
     &              ARRAYPIV(POSDIAG,POSPIVCOL), LDPIV,
     &              ARRAYCB(POSBLOCK), LDCB, NRHS_B, -1.0D0,
     &              K34, K451, KEEP8, PT_RHS, PT_RB(1), PT_QRB)
                ELSEIF (IBEG_BLOCK.LE.NPIV) THEN
C                 Entirely on FS side              
                  POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
                  CALL SMUMPS_LRMIXED_LRxFR(0, LRB,
     &              ARRAYPIV(POSDIAG,POSPIVCOL), LDPIV,
     &              ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV,
     &              NRHS_B, -1.0D0, K34, K451, KEEP8, 
     &              PT_RHS, PT_RB(1), PT_QRB)
                ELSE
C                 Entirely on CB side              
                  POSBLOCK = POSCB+int(IBEG_BLOCK-1-NPIV,8)
                  CALL SMUMPS_LRMIXED_LRxFR(0, LRB,
     &              ARRAYPIV(POSDIAG,POSPIVCOL), LDPIV,
     &              ARRAYCB(POSBLOCK), LDCB, NRHS_B, -1.0D0,
     &              K34, K451, KEEP8, PT_RHS, PT_RB(1), PT_QRB)
                ENDIF
              ENDIF
            ENDIF
          ELSE ! NOT ISLR
            IF (IS_T2_SLAVE) THEN
              POSBLOCK = POSCB + int(IBEG_BLOCK-1,8)
              CALL sgemm('N', 'N', M, NRHS_B, N, MONE,
     &              LRB%Q(1,1), M,  ARRAYPIV(POSDIAG,POSPIVCOL), LDPIV,
     &              ONE, ARRAYCB(POSBLOCK), LDCB) 
            ELSEIF (IBEG_BLOCK.LE.NPIV.AND.IEND_BLOCK.GT.NPIV) THEN
C             Case of overlap between the FS and non FS
              POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
              CALL sgemm('N', 'N', NPIV-IBEG_BLOCK+1, NRHS_B, N, MONE,
     &              LRB%Q(1,1), M,  ARRAYPIV(POSDIAG,POSPIVCOL), LDPIV,
     &              ONE, ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV) 
              CALL sgemm('N', 'N', IBEG_BLOCK+M-NPIV-1, NRHS_B, N, MONE,
     &              LRB%Q(NPIV-IBEG_BLOCK+2,1), M,
     &              ARRAYPIV(POSDIAG,POSPIVCOL),
     &              LDPIV, ONE, ARRAYCB(POSCB), LDCB) 
            ELSEIF (IBEG_BLOCK.LE.NPIV) THEN
C             Entirely on FS side              
              POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
              CALL sgemm('N', 'N', M, NRHS_B, N, MONE,
     &              LRB%Q(1,1), M,  ARRAYPIV(POSDIAG,POSPIVCOL), LDPIV,
     &              ONE, ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV)
            ELSE
C             Entirely on CB side              
              POSBLOCK = POSCB + int(IBEG_BLOCK-1-NPIV,8)
              CALL sgemm('N', 'N', M, NRHS_B, N, MONE,
     &              LRB%Q(1,1), M,  ARRAYPIV(POSDIAG,POSPIVCOL), LDPIV,
     &              ONE, ARRAYCB(POSBLOCK), LDCB) 
            ENDIF
          ENDIF
          IF (LRB%ISLR .AND. K450.GT.0
     &        .AND.(K448.EQ.0.OR. A_CHEVAL)) THEN
C           Deallocate UNIFORM part, keep MIXED for future solves            
            CALL DEALLOC_LRB(LRB, KEEP8, K34, .FALSE., .TRUE.)
          ENDIF
        ENDDO
#if ! defined(BLR_NOOPENMP)
!$OMP END DO
#endif
        IF (KMAX.GT.0) THEN
          IF (allocated(TEMP_BLOCK)) deallocate(TEMP_BLOCK)
        ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL         
#endif 
        IF (allocated(PT_QRB))   deallocate(PT_QRB)
        IF (allocated(TEMP_RB))  deallocate(TEMP_RB)
        IF (allocated(TEMP_QRB)) deallocate(TEMP_QRB)
        IF (allocated(PT_RHS))   deallocate(PT_RHS)
        IF (allocated(TEMP_RHS)) deallocate(TEMP_RHS)
        RETURN
      END SUBROUTINE SMUMPS_SOL_FWD_BLR_UPDATE
! (MIXHPC)
#endif 
      SUBROUTINE SMUMPS_SOL_BWD_LR_SU 
     &   ( INODE, IWHDLR, NPIV_GLOBAL, NSLAVES,
     &     LIELL, WCB, LWCB, NRHS_B, PTWCB,
     &     RHSINTR, LRHSINTR, NRHS,
     &     IPOSINRHSINTR, JBDEB, 
     &     MTYPE, KEEP, KEEP8,
     &     IFLAG, IERROR )
C
C  Purpose:
C  =======
C
C      Compute UX=Y (MTYPE=1 and KEEP(50).EQ.0) or LTX=Y of a BLR node 
C      Updates followed by Solves are performed by panel of BLR blocks
C
C      At each iteration (BLR block of the FS block) of the loop :
C       we process a COMPLETE panel (OFF diag block then block pivot) 
C       We first compute a 
C        GEMV/GEMM to update part of solution
C        TRSV/TRSM on the pivot block 
C
C  
       INTEGER, INTENT(IN) :: INODE, IWHDLR, NPIV_GLOBAL, NSLAVES
C      PCB_BEG_GLOBAL : Position of first entry of current front in WCB
       INTEGER, INTENT(IN) :: MTYPE, LIELL, KEEP(500)
       INTEGER(8), INTENT(IN) :: KEEP8(150)
       INTEGER, INTENT(IN) :: IPOSINRHSINTR, JBDEB, LRHSINTR, NRHS
       INTEGER(8), INTENT(IN) :: LWCB, PTWCB
C      NRHS_B size of current block of RHS beeing processed
       INTEGER, INTENT(IN) ::  NRHS_B
       INTEGER, INTENT(INOUT) :: IFLAG, IERROR
C      Note that LDA of WCB is LIELL
       REAL, INTENT(INOUT) :: WCB(LWCB)
       REAL RHSINTR(LRHSINTR,NRHS)
C      MTYPE
C             1 and KEEP(50)=0 UX=Y bloc of factors
C             else  LT X = Y
C
C Local variables
       INTEGER :: I, NPARTSASS, NB_BLR, LAST_BLR, 
     &            NELIM_PANEL, LD_WCB,
     &            DIAGSIZ_DYN, DIAGSIZ_STA, LDADIAG,
     &            IEND_BLR, IBEG_BLR
       INTEGER(8) :: PWCB
       INTEGER    :: IPIV_PANEL
       REAL, POINTER, DIMENSION(:) :: DIAG
       TYPE(LRB_TYPE), POINTER, DIMENSION(:) :: BLR_PANEL
       TYPE(LRMIXED_TYPE), ALLOCATABLE :: PT_RHS(:,:)
       INTEGER, ALLOCATABLE, TARGET, DIMENSION(:) :: TEMP_RHS
       INTEGER(8) :: POSBLOCK, POSCB, POSPIV
       INTEGER :: POSPIVCOL
       INTEGER :: M, PREC, NB_INT_B, POS_B
       INTEGER    :: NFRONT_OFFDIAG
       INTEGER :: allocok, K34
#if !defined(NOMIXHPC)
!       COMPLEX, ALLOCATABLE, DIMENSION(:), TARGET :: WORKSPACE
!       INTEGER :: WSMAX
        REAL, ALLOCATABLE, TARGET :: WORKSPACE_INCASE(:)
        REAL, DIMENSION(:), POINTER :: WORKSPACE
        INTEGER :: WSMAX
        INTEGER :: CACHE_L2
        integer(c_intptr_t) :: th_ws_int
        type(c_ptr) :: th_ws_ptr
#endif
        integer :: th_num
C
C    Constants
C    =========
       REAL :: ONE, MONE, ZERO
       PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
       PARAMETER (ZERO=0.0E0)
C
       IF ((MTYPE.EQ.1).AND.(KEEP(50).EQ.0)) THEN
         IF (associated(BLR_ARRAY(IWHDLR)%PANELS_U))
     &      THEN
            NPARTSASS=size(BLR_ARRAY(IWHDLR)%PANELS_U)
            NB_BLR = size(BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC) -1
         ENDIF
       ELSE
C      Uncompress L panels in A
         IF (associated(BLR_ARRAY(IWHDLR)%PANELS_L))
     &      THEN
            NPARTSASS=size(BLR_ARRAY(IWHDLR)%PANELS_L)
            NB_BLR = size(BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC) -1
         ELSE
           WRITE(6,*) " Internal error in SMUMPS_SOL_FWD_SU_MASTER"
         ENDIF
       ENDIF
C  ------------------------------
C  Variables independent of panel:
C  ------------------------------
C     RHS(IPOSINRHSINTR_NPIV_GLOBAL-1) 
C        last fully summed variable of the front
C     WCB(PTWCB+LIELL-1) = last entry in WCB of front variables
      PWCB        = PTWCB + int(NPIV_GLOBAL,8)
      LD_WCB      = LIELL
C  --------------------------------------
C  Variables updated within loop on panel
C  --------------------------------------
C     IPIV_PANEL 
C       RHSINTR(IPIV_PANEL,:) will be updated to hold the 
C                            position of first pivot in current panel
      IF (KEEP(50).EQ.0 .AND. NSLAVES.GT.0 .AND. MTYPE.NE.1) THEN
C       Type 2 unsymmetric transpose, use L panels
C       so we must stop at block NPARTSASS
        LAST_BLR = NPARTSASS
      ELSE
C       Other cases: update also CB part
        LAST_BLR = NB_BLR
      ENDIF
C     Copy RHS in each precision
      IF (KEEP(448).GE.1.AND.KEEP(450).GT.0) THEN
        POSPIV = int(IPOSINRHSINTR,8)
        POSPIVCOL = JBDEB
        POSCB = PWCB
        K34 = KEEP(34)
C       Front size without first diagonal block, which we don't need to copy        
        NFRONT_OFFDIAG = BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC(LAST_BLR+1) 
     &                 - BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC(2)
C       Count size of all copies of RHS in each prec        
        NB_INT_B = 0
        DO PREC=2,NPREC_CALC
          NB_INT_B = NB_INT_B + ( NFRONT_OFFDIAG * NRHS_B *
     &                        COSTS_PREC(CALC2STOC(PREC)) +K34-1) / K34
        ENDDO
C       Workspace containing all copies         
        allocate(TEMP_RHS(NB_INT_B), stat=allocok)
        if (allocok.gt.0) then
          IFLAG = -13
          IERROR = NB_INT_B
          write(*,*) 'Allocation problem in BLR routine 
     &    SMUMPS_SOL_SLAVE_LR_U',
     &    'not enough memory? memory requested = ', IERROR
          CALL MUMPS_ABORT()
        endif
C       Pointer that points to each block copy in each precision
        allocate(PT_RHS(NPREC_CALC,LAST_BLR), stat=allocok)
        if (allocok.gt.0) then
          IFLAG = -13
          IERROR = NPREC_CALC * NB_BLR
          write(*,*) 'Allocation problem in BLR routine 
     &    SMUMPS_SOL_SLAVE_LR_U',
     &    'not enough memory? memory requested = ', IERROR
          CALL MUMPS_ABORT()
        endif
C       Nullify unused pointers (I=1 and PREC=1)
        DO PREC=2,NPREC_CALC
          nullify(PT_RHS(PREC,1)%val)
        ENDDO
        DO I=1,LAST_BLR
          nullify(PT_RHS(1,I)%val)
        ENDDO
C       Cast RHS
        POS_B = 1
        DO I = NPARTSASS+1,LAST_BLR
          IBEG_BLR = BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC(I)
          M = BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC(I+1) 
     &                  - IBEG_BLR
          POSBLOCK = POSCB +int(IBEG_BLR-1-NPIV_GLOBAL,8)
          DO PREC=2,NPREC_CALC
            PT_RHS(PREC,I)%val => TEMP_RHS(POS_B:)
            POS_B = POS_B + (M * NRHS_B * 
     &                 COSTS_PREC(CALC2STOC(PREC)) +K34-1) / K34
            CALL SMUMPS_LRMIXED_CAST(WCB(POSBLOCK),
     &          PT_RHS(PREC,I)%val(1), M, 
     &          NRHS_B, LD_WCB, M, 1, CALC2STOC(PREC))
          ENDDO
        ENDDO
#if defined(MUMPS_NOF2003)
      ELSE
C       PT_RHS can not be declared allocatable in
C       SMUMPS_SOL_BWD_BLR_UPDATE, it is declared
C       as an array and so it must be allocated.
        ALLOCATE(PT_RHS(1,1))
#endif
      ENDIF
      th_num = 0
      if (KEEP(405).EQ.1) then
          th_num = thread_id - 1
      endif
#if !defined(NOMIXHPC)
      WSMAX = 0
      IF (PREC_DIAG_FR.gt.1) THEN
         DO I=1, NPARTSASS
           WSMAX = MAX(WSMAX, 
     &         BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE)
         END DO
CAJ2ALL FIXME this allocation (and the one done in (f/b)wd_update) is
C             not accounted for in the memory counters
      ENDIF
      CACHE_L2 = KEEP(365) / COSTS_PREC(1)
      IF (PREC_DIAG_FR.GT.1) THEN
        IF (WSMAX.GT.CACHE_L2) THEN
CAJ2ALL FIXME this allocation (and the one done in (f/b)wd_update) is
C               not accounted for in the memory counters
          ALLOCATE(WORKSPACE_INCASE(WSMAX))
          WORKSPACE => WORKSPACE_INCASE
        ELSE
          call mumps_sol_omp_mem_manager_get_address(th_num, th_ws_int)
          th_ws_ptr = transfer(th_ws_int, th_ws_ptr)
          call c_f_pointer(th_ws_ptr, workspace, (/ CACHE_L2 /))
        ENDIF
      ENDIF
#endif
C     Panel loop (start from last panel)
      DO I=NPARTSASS,1,-1
        IBEG_BLR = BLR_ARRAY(IWHDLR)%BEGS_BLR_DYNAMIC(I)
        IEND_BLR =  BLR_ARRAY(IWHDLR)%BEGS_BLR_DYNAMIC(I+1) -1
        DIAGSIZ_DYN = BLR_ARRAY(IWHDLR)%BEGS_BLR_DYNAMIC(I+1) -
     &           IBEG_BLR
        DIAGSIZ_STA = BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC(I+1)  - 
     &           IBEG_BLR 
        IF (KEEP(50).NE.0) THEN
         LDADIAG = DIAGSIZ_DYN
        ELSE
         LDADIAG = DIAGSIZ_STA
        ENDIF 
C       Skip empty blocks
        IF (DIAGSIZ_DYN.EQ.0) GOTO 1000
C       NELIM_PANEL has to be defined by panel
        NELIM_PANEL = DIAGSIZ_STA - DIAGSIZ_DYN
C       or also NELIM_PANEL =    
C     &         BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC(I+1) -
C     &          BLR_ARRAY(IWHDLR)%BEGS_BLR_DYNAMIC(I+1)
         IPIV_PANEL = IPOSINRHSINTR + IBEG_BLR -1
         IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0) THEN
C         ---------
C         U X = Y
C         ---------
C         Since U^T=QR already => nothing special to do
          BLR_PANEL => BLR_ARRAY(IWHDLR)%PANELS_U(I)%LRB_PANEL
         ELSE
C         ---------
C         L^T  X = Y
C         ---------
          BLR_PANEL => BLR_ARRAY(IWHDLR)%PANELS_L(I)%LRB_PANEL
         END IF
C       ===========================================================
C       I/ Update Solution with computed parts of solution
C       ===========================================================
            IF (.NOT.allocated(PT_RHS))  THEN
               ALLOCATE(PT_RHS(1,1))
            ENDIF
        CALL SMUMPS_SOL_BWD_BLR_UPDATE ( 
     &        RHSINTR, int(LRHSINTR,8), NRHS, LRHSINTR, 
     &        int(IPOSINRHSINTR,8), JBDEB,
     &        WCB, LWCB, LD_WCB, PWCB, 
     &        int(IPIV_PANEL,8), ! POSDIAG 
     &        NRHS_B, NPIV_GLOBAL,
     &        BLR_PANEL, LAST_BLR,
     &        I, BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC,
     &        KEEP8, KEEP(34), 
#if !defined(NOMIXHPC)
     &        KEEP(148),
     &        KEEP(365),
     &        KEEP(405), ! I need to know if I'm under L0
#endif
     &        KEEP(448), KEEP(450), KEEP(451),
#if !defined(NOMIXHPC)
     &          KEEP(455),
#endif
     &        .FALSE., 
     &        PT_RHS,
     &        IFLAG, IERROR)
        IF (IFLAG.LT.0) RETURN
C       ===========================================================
C       II/ Process Diagonal block : 
C           I.0/ Cast to working precision if necessary
C           I.1/ Update Solution with solution concerning the 
C                NELIM_PANEL var at this step
C           I.2/ Triangular compute part of solution
C       ===========================================================
C        I.0
         DIAG => BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_BLOCK
#if !defined(NOMIXHPC)
C         write(*,*) "[CSV]", INODE, "BWD", PREC_DIAG_FR, 
C     &     BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE
         IF (PREC_DIAG_FR.gt.1) THEN
           DIAG(1:BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE) => 
     &       WORKSPACE(1:BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE)
           CALL SMUMPS_LRMIXED_CAST(
     &       BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_STOC%VAL(1), DIAG(1),
     &       BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE, 1,
     &       BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE,
     &       BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE,
     &       PREC_DIAG_FR, 1)
         ELSE IF (.not.ASSOCIATED(DIAG)) THEN
           call c_f_pointer(
     &        c_loc(BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_STOC%VAL),
     &        DIAG, (/ BLR_ARRAY(IWHDLR)%DIAG_BLOCKS(I)%DIAG_SIZE/) )
         ENDIF
#endif
C        I.1 
         IF (NELIM_PANEL.GT.0) THEN
           IF (MTYPE.EQ.1.AND.KEEP(50).EQ.0) THEN
             IF (IEND_BLR.EQ.NPIV_GLOBAL) THEN
               CALL sgemm('T', 'N', DIAGSIZ_DYN, NRHS_B, NELIM_PANEL, 
     &              MONE, DIAG(1+DIAGSIZ_DYN), DIAGSIZ_STA,  WCB(PWCB),
     &              LD_WCB, ONE , RHSINTR(IPIV_PANEL,JBDEB),LRHSINTR)
             ELSE
               IF (IEND_BLR+1.LE.NPIV_GLOBAL .AND.
     &            IEND_BLR+NELIM_PANEL.GT.NPIV_GLOBAL) THEN
                 CALL sgemm('T', 'N', DIAGSIZ_DYN, NRHS_B,
     &              NPIV_GLOBAL-IEND_BLR, 
     &              MONE, DIAG(1+DIAGSIZ_DYN), DIAGSIZ_STA, 
     &              RHSINTR(IPIV_PANEL+DIAGSIZ_DYN,JBDEB), LRHSINTR,
     &              ONE, RHSINTR(IPIV_PANEL,JBDEB), LRHSINTR)
                 CALL sgemm('T', 'N', DIAGSIZ_DYN, NRHS_B,
     &              IEND_BLR+NELIM_PANEL-NPIV_GLOBAL, 
     &              MONE, DIAG(1+DIAGSIZ_DYN+NPIV_GLOBAL-IEND_BLR),
     &              DIAGSIZ_STA, 
     &              WCB(PWCB), LD_WCB,
     &              ONE, RHSINTR(IPIV_PANEL,JBDEB), LRHSINTR)
               ELSE
                 CALL sgemm('T', 'N', DIAGSIZ_DYN, NRHS_B, NELIM_PANEL, 
     &              MONE, DIAG(1+DIAGSIZ_DYN), DIAGSIZ_STA, 
     &              RHSINTR(IPIV_PANEL+DIAGSIZ_DYN,JBDEB), LRHSINTR,
     &              ONE, RHSINTR(IPIV_PANEL,JBDEB), LRHSINTR)
               ENDIF
             ENDIF
           ELSE
             IF (IEND_BLR.EQ.NPIV_GLOBAL) THEN
               CALL sgemm('N', 'N', DIAGSIZ_DYN, NRHS_B, NELIM_PANEL, 
     &              MONE, DIAG(1+DIAGSIZ_DYN*LDADIAG), DIAGSIZ_DYN, 
     &              WCB(PWCB), LD_WCB, ONE, 
     &              RHSINTR(IPIV_PANEL,JBDEB), LRHSINTR)
             ELSE
               IF (IEND_BLR+1.LE.NPIV_GLOBAL .AND.
     &            IEND_BLR+NELIM_PANEL.GT.NPIV_GLOBAL) THEN
                 CALL sgemm('N', 'N', DIAGSIZ_DYN, NRHS_B,
     &              NPIV_GLOBAL-IEND_BLR, 
     &              MONE, DIAG(1+DIAGSIZ_DYN*LDADIAG), DIAGSIZ_DYN, 
     &              RHSINTR(IPIV_PANEL+DIAGSIZ_DYN,JBDEB), LRHSINTR,
     &              ONE, RHSINTR(IPIV_PANEL,JBDEB), LRHSINTR)
                 CALL sgemm('N', 'N', DIAGSIZ_DYN, NRHS_B,
     &              IEND_BLR+NELIM_PANEL-NPIV_GLOBAL, 
     &              MONE, DIAG(1+DIAGSIZ_DYN*LDADIAG +
     &                  (NPIV_GLOBAL-IEND_BLR)*DIAGSIZ_DYN),
     &              DIAGSIZ_DYN, 
     &              WCB(PWCB), LD_WCB,
     &              ONE, RHSINTR(IPIV_PANEL,JBDEB), LRHSINTR)
               ELSE
                 CALL sgemm('N', 'N', DIAGSIZ_DYN, NRHS_B, NELIM_PANEL, 
     &              MONE, DIAG(1+DIAGSIZ_DYN*LDADIAG), DIAGSIZ_DYN, 
     &              RHSINTR(IPIV_PANEL+DIAGSIZ_DYN,JBDEB), LRHSINTR,
     &              ONE, RHSINTR(IPIV_PANEL,JBDEB), LRHSINTR)
               ENDIF
             ENDIF
           ENDIF
         ENDIF
         IF (IFLAG.LT.0) RETURN
C        I.2:
         CALL SMUMPS_SOLVE_BWD_LR_TRSOLVE  (
     &         DIAG(1), size(DIAG), DIAGSIZ_DYN, NELIM_PANEL, LIELL,
     &         NRHS_B, WCB, LWCB, 
     &         RHSINTR, LRHSINTR, NRHS,
     &         IPIV_PANEL, JBDEB, 
     &         MTYPE, KEEP )
 1000 CONTINUE         
C       Cast block of RHS that has just been computed         
        IF (KEEP(448).GE.1.AND.KEEP(450).GT.0 .AND. I.GT.1) THEN
          IBEG_BLR = BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC(I)
          M = BLR_ARRAY(IWHDLR)%BEGS_BLR_STATIC(I+1) - IBEG_BLR
          DO PREC=2,NPREC_CALC
            PT_RHS(PREC,I)%val => TEMP_RHS(POS_B:)
            POS_B = POS_B + ( M * NRHS_B * 
     &                 COSTS_PREC(CALC2STOC(PREC)) +K34-1) / K34
            IF (IBEG_BLR.GT.NPIV_GLOBAL) THEN
C             CB side
              POSBLOCK = POSCB +int(IBEG_BLR-1-NPIV_GLOBAL,8)
              CALL SMUMPS_LRMIXED_CAST(WCB(POSBLOCK),
     &          PT_RHS(PREC,I)%val(1), M, 
     &          NRHS_B, LD_WCB, M, 1, CALC2STOC(PREC))
            ELSEIF (IEND_BLR.LE.NPIV_GLOBAL) THEN
C             FS side
              POSBLOCK = POSPIV + int(IBEG_BLR-1,8)
              CALL SMUMPS_LRMIXED_CAST(RHSINTR(POSBLOCK,POSPIVCOL),
     &          PT_RHS(PREC,I)%val(1), M, 
     &          NRHS_B, LRHSINTR, M, 1, CALC2STOC(PREC))
            ELSE
C             A cheval              
C             In current version these blocks are handled in double
C             precision, so no cast is necessary              
            ENDIF
          ENDDO
        ENDIF
C     -- end of loops on panels
      ENDDO  
#if !defined(NOMIXHPC)
      IF (allocated(WORKSPACE_INCASE)) DEALLOCATE(WORKSPACE_INCASE)
#endif
      IF (allocated(PT_RHS))   deallocate(PT_RHS)
      IF (allocated(TEMP_RHS)) deallocate(TEMP_RHS)
C
      RETURN
      END SUBROUTINE SMUMPS_SOL_BWD_LR_SU
#if !defined(NOMIXHPC)
      SUBROUTINE SMUMPS_SOL_BWD_BLR_UPDATE (
     &        ARRAYPIV, LPIV, LPIVCOL, LDPIV, POSPIV, POSPIVCOL,
     &        ARRAYCB, LCB, LDCB, POSCB, 
     &        POSDIAG,
     &        NRHS_B, NPIV,
     &        BLR_PANEL, LAST_BLR, CURRENT_BLR,
     &        BEGS_BLR_STATIC,
     &        KEEP8, K34, K148, K365, K405, K448, K450, K451, K455,
     &        IS_T2_SLAVE,
     &        PT_RHS, 
     &        IFLAG, IERROR)
C
C Inspired from SMUMPS_BLR_UPD_NELIM_VAR_L
C     
C     Update of the solution phase of current panel 
C     -----------
C     Parameters
C     -----------
C
C     ARRAYPIV is the array containing the FS (eliminated) rows of the solution (WCB
C     or RHSINTR depending on K350 in the future). It is of size LPIV
C     and leading dimension LDPIV. POSPIV points to the first eliminated row
C     associated to the current front in ARRAYPIV and is not used for
C     type 2 slaves
C     ARRAYCB is the array containing the NFS+NELIM rows of the solution
C     (always WCB). It is of size LCB and leading dimension LDCB. 
C     POSCB points to the first non-eliminated row associated to the current front
C     in ARRAYCB. 
C     POSDIAG points to the beginning of the current panel CURRENT_BLR
C     (in ARRAYPIV).        
C     NPIV is the global NPIV of the front
C     LPIVCOL is either 1 (ARRAYPIV=WCB) or NRHS (ARRAYPIV=RHSINTR)        
      implicit none
      INTEGER(8), INTENT(IN) :: LPIV, LCB, POSPIV, POSCB, POSDIAG
      INTEGER,INTENT(IN)     :: LPIVCOL, POSPIVCOL
      REAL, INTENT(INOUT) :: ARRAYPIV(LPIV,LPIVCOL)
      REAL, INTENT(INOUT) :: ARRAYCB(LCB)
      INTEGER, INTENT(IN)    :: LAST_BLR, NRHS_B,  LDPIV, LDCB,
     &                          CURRENT_BLR, NPIV
      INTEGER, INTENT(IN)    :: K34, K148, K365, K405
      INTEGER, INTENT(IN)    :: K448, K450, K451, K455
c     TYPE(LRB_TYPE), TARGET,INTENT(IN) ::
c    &                          BLR_PANEL(LAST_BLR-CURRENT_BLR)
      TYPE(LRB_TYPE), TARGET,INTENT(IN) ::
     &                          BLR_PANEL(:)
      INTEGER(8), INTENT(IN) :: KEEP8(150)
      LOGICAL, INTENT(IN)    :: IS_T2_SLAVE
      INTEGER :: BEGS_BLR_STATIC(:)
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
      TYPE(LRMIXED_TYPE), target :: PT_RHS(:,:)
C     -----------
C     Local variables
C     -----------
      INTEGER :: I, K, M, N, IBEG_BLOCK, IEND_BLOCK
      INTEGER :: KMAX, WSMAX
      INTEGER :: J1, J2
      INTEGER(8) :: POSBLOCK
      TYPE(LRB_TYPE), POINTER :: LRB
      REAL, ALLOCATABLE, DIMENSION(:) :: TEMP_BLOCK
      REAL, ALLOCATABLE, DIMENSION(:) :: DEST_ARRAY
      integer(c_intptr_t) :: th_ws_int
      type(c_ptr) :: th_ws_ptr
      REAL, POINTER :: thread_workspaces(:) 
      integer, allocatable, target :: val_workspace(:)
      REAL, pointer :: Q_ptr(:,:), R_ptr(:,:)
      INTEGER :: allocok
      integer :: kcumul, kmagic, krem, kint, prec, kprec
      integer :: start_prec, cache_L2
      integer(kind=1), pointer, dimension(:) :: stoc_ptr
      REAL :: ONE, MONE, ZERO
      ! for cache stuff 
      DOUBLE PRECISION :: total_consumption, current_consumption
      DOUBLE PRECISION :: overall_consumption
      integer :: budget, budget_ws
      integer :: ws_footprint, prec_footprint
      integer :: start_k, current_k, current_prec, used_k
      integer :: ptr2_shape(2), ptr1_shape(1)
      real :: nu
      integer :: overall_k, end_prec, end_k
      TYPE(LRMIXED_TYPE),POINTER :: ptPT_RHS(:)
      PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
      PARAMETER (ZERO=0.0E0)
      INTEGER, ALLOCATABLE, TARGET, DIMENSION(:) :: TEMP_RB, TEMP_QRB
      INTEGER :: POS_QRB
      TYPE(LRMIXED_TYPE), ALLOCATABLE :: PT_QRB(:)
      INTEGER, POINTER:: PT_RB(:)
      INTEGER :: NB_INT_RB, NB_INT_QRB
      INTEGER :: PREC_CALC, PREC_
      LOGICAL :: MIXED_STORED
      LOGICAL :: A_CHEVAL
      INTEGER :: th_num
      CACHE_L2 = K365
      MIXED_STORED = .not.(K450.LE.0.AND.K148.NE.2)
      KMAX = -1
      WSMAX = -1
      DO I = CURRENT_BLR+1, LAST_BLR
        LRB => BLR_PANEL(I-CURRENT_BLR)
        KMAX = max(KMAX, LRB%K)
        IF (.NOT.MIXED_STORED) CYCLE
C        CALL COMPUTE_WORKSPACE_SIZE(LRB, WSMAX,
C     &    K34, K148, K448, K450, K455, CACHE_L2,NRHS_B)
      ENDDO
      IF (CURRENT_BLR.LT.LAST_BLR) THEN 
C       N is constant for all blocks in the panel        
        N = BLR_PANEL(1)%N
      ELSE
        RETURN
      ENDIF
      allocate(DEST_ARRAY(N*NRHS_B),stat=allocok)
      IF (allocok .GT. 0) THEN
        IFLAG = -13
        IERROR = N * NRHS_B
        GOTO 100
      ENDIF
      DEST_ARRAY = ZERO
      WSMAX = CACHE_L2 / 8
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL PRIVATE(TEMP_BLOCK, 
!$OMP& th_num, th_ws_ptr, th_ws_int,
!$OMP& Q_ptr, R_ptr, stoc_ptr, thread_workspaces,
!$OMP& ptr1_shape, ptr2_shape,
!$OMP& start_prec, prec, kcumul, kprec, kmagic, kint, krem,
!$OMP& end_prec, 
!$OMP& j1, j2,
!$OMP& TEMP_QRB, TEMP_RB, 
!$OMP& PT_QRB, PT_RB, POS_QRB, PREC_, PREC_CALC, NB_INT_RB, NB_INT_QRB,
!$OMP& ptpt_rhs, val_workspace,
!$OMP& allocok)
#endif 
      th_num = 0
!$    th_num = omp_get_thread_num()
      if (K405.EQ.1) then
        th_num = thread_id - 1
      endif
      IF (K448.GE.1.AND.MIXED_STORED) THEN
C       Allocate RB and QRB workspaces for each thread        
        IF(CURRENT_BLR.LT.LAST_BLR) THEN
          NB_INT_QRB = 0
          DO PREC=2,NPREC_CALC
            NB_INT_QRB = NB_INT_QRB + 
     &         (N*NRHS_B*COSTS_PREC(CALC2STOC(PREC)) +K34-1) / K34
          ENDDO
          NB_INT_RB=(KMAX*NRHS_B*COSTS_PREC(CALC2STOC(1)) + K34-1) / K34
          allocate(TEMP_RB(NB_INT_RB), stat=allocok)
          if (allocok.gt.0) then
!$OMP ATOMIC WRITE
            IFLAG = -13
!$OMP END ATOMIC
            IERROR = NB_INT_RB
            write(*,*) 'Allocation problem in BLR routine 
     &      SMUMPS_SOL_FWD_BLR_UPDATE over temp_rb',
     &      'not enough memory? memory requested = ', IERROR
          endif
          allocate(TEMP_QRB(NB_INT_QRB), stat=allocok)
          if (allocok.gt.0) then
!$OMP ATOMIC WRITE
            IFLAG = -13
!$OMP END ATOMIC
            IERROR = NB_INT_QRB
            write(*,*) 'Allocation problem in BLR routine 
     &      SMUMPS_SOL_FWD_BLR_UPDATE over temp_qrb',
     &      'not enough memory? memory requested = ', IERROR
          endif
          allocate(PT_QRB(NPREC_CALC), stat=allocok)
          if (allocok.gt.0) then
!$OMP ATOMIC WRITE
            IFLAG = -13
!$OMP END ATOMIC
            IERROR = NPREC_CALC
            write(*,*) 'Allocation problem in BLR routine 
     &      SMUMPS_SOL_FWD_BLR_UPDATE over pt_qrb',
     &      'not enough memory? memory requested = ', IERROR,
     &      'stat value', allocok
          endif
          nullify(PT_QRB(1)%val)
          POS_QRB = 1
          DO PREC=2,NPREC_CALC
            PT_QRB(PREC)%VAL => TEMP_QRB(POS_QRB:)
            POS_QRB = POS_QRB + 
     &         (N*NRHS_B*COSTS_PREC(CALC2STOC(PREC)) +K34-1) / K34
          ENDDO
          PT_RB=>TEMP_RB(1:NB_INT_RB)
        ENDIF
      ENDIF
      IF (MIXED_STORED) THEN
        call mumps_sol_omp_mem_manager_get_address(th_num, th_ws_int)
        th_ws_ptr = transfer(th_ws_int, th_ws_ptr)
        call c_f_pointer(th_ws_ptr, thread_workspaces, (/ WSMAX /))
      ENDIF
      IF (KMAX.GT.0) THEN
C       TEMP_BLOCK is allocated even when K450.GT.0 because we might
C       need it when A_CHEVAL is true        
        allocate(TEMP_BLOCK(KMAX*NRHS_B), stat=allocok )
        IF (allocok .GT. 0) THEN
          IFLAG  = -13
          IERROR = NRHS_B * KMAX
          write(*,*) 'Allocation problem in BLR routine 
     &    SMUMPS_SOL_BWD_BLR_UPDATE: ',
     &    'not enough memory? memory requested = ', IERROR
        ENDIF
      ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP DO SCHEDULE(DYNAMIC,1)
!$OMP& PRIVATE(IBEG_BLOCK, IEND_BLOCK, LRB, K, M, POSBLOCK)
!$OMP& PRIVATE(current_consumption, current_k, current_prec,
!$OMP&         total_consumption, nu,
!$OMP&         overall_consumption, budget, budget_ws,
!$OMP&         overall_k, ws_footprint, prec_footprint,
!$OMP&         used_k, start_k, end_k)
!$OMP& PRIVATE(A_CHEVAL)
!$OMP& REDUCTION(+:DEST_ARRAY)
#endif
        DO I = CURRENT_BLR+1, LAST_BLR
          IF (IFLAG.LT.0) CYCLE
          IBEG_BLOCK = BEGS_BLR_STATIC(I)
          IEND_BLOCK = BEGS_BLR_STATIC(I+1)-1
          LRB => BLR_PANEL(I-CURRENT_BLR)
C Case of U: U^T = Q R
          K = LRB%K 
          M = LRB%M   ! always nb of rows of Q
!          write(*,*) "[DEV] SZ BWD BLR", M, N, K
C
          IF (LRB%ISLR) THEN
            IF (K.LE.0) CYCLE
C             Check if block is "a cheval": if so, do not use LRMIXED
C             for simplicity (we could try to code this in the future)
            A_CHEVAL =  (.NOT.IS_T2_SLAVE) .AND. (IBEG_BLOCK.LE.NPIV) 
     &            .AND. (IEND_BLOCK.GT.NPIV)
            IF (K448.EQ.0.OR..NOT.MIXED_STORED.OR. A_CHEVAL) THEN
              if (.NOT.MIXED_STORED) then
                CALL GEMM_Q_BWD(M, NRHS_B, K, NPIV, 
     &              LRB%Q, TEMP_BLOCK, K, 
     &              ARRAYPIV, lpiv, ldpiv, ARRAYCB, lcb, ldcb, 
     &              ibeg_block, iend_block, is_t2_slave, 
     &              poscb, pospiv, pospivcol, 1)
              else if (k455.ge.3) then
                kcumul = 1
                if (lrb%k_stoc(1).gt.0) then
                  ! handle working arithmetic 
                  kprec = lrb%k_stoc(1)
                  ptr2_shape(1) = M
                  ptr2_shape(2) = Kprec
                  call c_f_pointer(c_loc(LRB%Q_stoc(1)%val),
     &                             Q_ptr, ptr2_shape) !(/M, Kprec/))
                  CALL GEMM_Q_BWD(M, NRHS_B, Kprec, NPIV, 
     &                Q_ptr, TEMP_BLOCK, K, 
     &                ARRAYPIV, lpiv, ldpiv, ARRAYCB, lcb, ldcb, 
     &                ibeg_block, iend_block, is_t2_slave, 
     &                poscb, pospiv, pospivcol, kcumul)
                end if 
                overall_k = lrb%k_stoc(1)
                kprec = lrb%k - lrb%k_stoc(1)
                kcumul = kcumul + lrb%k_stoc(1)
                total_consumption = m*( kprec * costs_prec(1) + 
     &                                  lrb%kw_stoc ) / 
     &                              dble((4.0/5.0)*cache_L2)
C                if (kprec.gt.0) then
C                  Q_ptr(1:M,1:kprec) => 
C     &                thread_workspaces(1:(M*kprec))
C                  CALL SMUMPS_DECONVERT_LRB_Q_W(LRB, Q_ptr, 
C     &                M, kprec,
C     &                2, NPREC_STOC, 1, LRB%K_STOC(NPREC_STOC),
C     &                IFLAG, IERROR, KEEP8)
C                  CALL GEMM_Q_BWD(M, NRHS_B, kprec, NPIV, 
C     &                Q_ptr, TEMP_BLOCK, K, 
C     &                ARRAYPIV, ldpiv, ARRAYCB, lcb, ldcb, 
C     &                ibeg_block, iend_block, is_t2_slave, 
C     &                poscb, pospiv, pospivcol, kcumul)
C                end if
                start_prec = 2
                budget = 0
                budget_ws = 0
                start_k = 1
                current_k = 0
                overall_k = lrb%k_stoc(1)
                nu = 1.0d0
                if (k455.eq.4) nu = 0.2d0
                if (k455.eq.5.and.lrb%M.gt.200) nu = 0.2d0
                current_consumption = 0.0d0
                overall_consumption = 0.0d0
                do prec = 2, nprec_stoc
                  if (lrb%k_stoc(prec).eq.0) cycle
                  krem = lrb%k_stoc(prec)
                  current_prec = 1
                  do
                    ! avoid small blas when cache is barely over filled
                    if (k455.eq.3.and.nu.eq.1.0d0
     &          .and. total_consumption - overall_consumption.gt.1.0D0
     &          .and. total_consumption - overall_consumption.lt.1.5D0
     &              ) then
                      nu = (total_consumption - overall_consumption)/2.0
                      nu = nu + 0.02
                    end if
                    kmagic = max(dble(nu*(4.0/5.0)*CACHE_L2
     &                        - budget - budget_ws),0.0d0) /
     &                       (lrb%M*(costs_prec(prec)+costs_prec(1)))
                    used_k = min(kmagic, krem)
                    ws_footprint = lrb%M * used_k * costs_prec(1)
                    prec_footprint = lrb%M * used_k * costs_prec(prec)
                    budget = budget + prec_footprint
                    budget_ws = budget_ws + ws_footprint
                    current_consumption = current_consumption +
     &                dble(prec_footprint+ws_footprint) / 
     &                ((80.0/100.0)*dble(cache_L2))
                    overall_consumption = overall_consumption +
     &                dble(prec_footprint+ws_footprint) / 
     &                ((80.0/100.0)*dble(cache_L2))
                    current_prec = current_prec + used_k
                    current_k = current_k + used_k
                    krem = krem - used_k
                    overall_k = overall_k + used_k
                    if (krem.ge.1.and.k455.eq.4
     &                  .and.current_consumption.gt.0.3.or.
     &                  krem.ge.1.or.prec.eq.nprec_stoc.or.
     &                  overall_k.eq.k) then
                      Q_ptr(1:M,1:current_k) => 
     &                    thread_workspaces(1:(M*current_k))
                      end_prec = prec
                      end_k = lrb%k_stoc(prec) - krem
                      ! if we switched prec but even a row wont fit
                      if (krem.eq.lrb%k_stoc(prec)) then
                        end_prec = prec - 1
                        end_k = lrb%k_stoc(end_prec)
                      end if
                      CALL SMUMPS_DECONVERT_LRB_Q_W(LRB, Q_ptr, 
     &                    M, current_k,
     &                    start_prec, end_prec, 
     &                    start_k, end_k,
     &                    IFLAG, IERROR, KEEP8)
C                      write(*,*) "MYSTUFF", m, current_k, k
                      CALL GEMM_Q_BWD(M, NRHS_B, current_K, NPIV, 
     &                    Q_ptr, TEMP_BLOCK, K, 
     &                    ARRAYPIV, lpiv, ldpiv, ARRAYCB, lcb, ldcb, 
     &                    ibeg_block, iend_block, is_t2_slave, 
     &                    poscb, pospiv, pospivcol, kcumul)
                      kcumul = kcumul + current_k
                      start_prec = prec
                      budget = 0
                      budget_ws = 0
                      start_k = current_prec
                      current_consumption = 0.0d0
                      current_k = 0
                    end if
                    if (krem.lt.1) then 
                      exit
                    end if
                  end do
                end do
              else
                kcumul = lrb%k + 1
                start_prec=nprec_stoc
                ! if the sum of R^T (except working prec.) fit in cache
                kprec = lrb%k - lrb%k_stoc(1)
                if (kprec.gt.0.and.
     &              (k455.eq.0.or.k455.ne.1.and.
     &              M*LRB%kw_stoc +
     &              M*costs_prec(1) * kprec < (4.0/5)*cache_L2)) then
                  kcumul = kcumul - kprec
                  Q_ptr(1:M,1:kprec) =>
     &                thread_workspaces(1:(M*Kprec))
                  CALL SMUMPS_DECONVERT_LRB_Q_W(LRB, 
     &                Q_ptr, M, kprec,
     &                2, NPREC_STOC, 1, LRB%K_STOC(NPREC_STOC),
     &                IFLAG, IERROR, KEEP8)
C                  write(*,*) "MYSTUFF", m, kprec, k
                  CALL GEMM_Q_BWD(M, NRHS_B, Kprec, NPIV, 
     &                Q_ptr, TEMP_BLOCK, K, 
     &                ARRAYPIV, lpiv, ldpiv, ARRAYCB, lcb, ldcb, 
     &                ibeg_block, iend_block, is_t2_slave, 
     &                poscb, pospiv, pospivcol, kcumul)
                  ! working prec is adressed in following do loop
                  start_prec=1
                endif
                ! else (e.g. if any R^T_prec is too big)
                do prec=start_prec,1,-1
                  kprec = lrb%k_stoc(prec)
                  if (kprec.eq.0) cycle
                  kcumul = kcumul - kprec
                  if (prec.eq.1) then
                    kmagic = kprec
                    ptr2_shape(1) = LRB%M
                    ptr2_shape(2) = Kprec
                    call c_f_pointer(c_loc(LRB%Q_stoc(prec)%val),
     &                       Q_ptr, ptr2_shape)!(/LRB%M, Kprec/))
                  else
                    KMAGIC=(80.0/100.0)*(CACHE_L2-costs_prec(1)*M)/
     &              ( (costs_prec(1)+costs_prec(prec))*M + 
     &                NRHS_B*costs_prec(1) )
                    ptr1_shape(1) = costs_prec(prec)*LRB%M*Kprec
                    call c_f_pointer(c_loc(LRB%Q_stoc(prec)%val),
     &                    stoc_ptr, ptr1_shape)
                    KMAGIC = KMAGIC - MODULO(KMAGIC,8)
                  endif
                  if (k455.eq.1) kmagic = kprec
                  do kint=1,Kprec,KMAGIC
                    !Krem may be stupidly small
                    krem = min(KMAGIC,Kprec-Kint+1)
                    if (prec.gt.1) then
                      Q_ptr(1:M,1:Krem) =>
     &                  thread_workspaces(1:(M*krem))
C                      write(*,*) "MYSTUFF", m, krem, k
                      CALL SMUMPS_LRMIXED_CAST(
     &                  stoc_ptr(m*costs_prec(prec)*
     &                             (kint-1)+1),
     &                  Q_ptr(1,1), M, Krem,
     &                  M, M, PREC, 1)
                    endif
                    CALL GEMM_Q_BWD(M, NRHS_B, Krem, NPIV, 
     &                  Q_ptr, TEMP_BLOCK, K, 
     &                  ARRAYPIV, lpiv, ldpiv, ARRAYCB, lcb, ldcb, 
     &                  ibeg_block, iend_block, is_t2_slave, 
     &                  poscb, pospiv, pospivcol, kcumul+kint-1)
                  end do
                end do
              end if
              if (.not.MIXED_STORED) then
                CALL sgemm('T', 'N', N, NRHS_B, K, MONE,
     &            LRB%R(1,1), K, 
     &            TEMP_BLOCK(1), K, ONE,
     &            DEST_ARRAY(1), N)
              else if (k455.ge.3) then
                kcumul = 1
                if (lrb%k_stoc(1).gt.0) then
                  ! handle working arithmetic 
                  kprec = lrb%k_stoc(1)
                  ptr2_shape(1) = Kprec
                  ptr2_shape(2) = LRB%N
                  call c_f_pointer(c_loc(LRB%R_stoc(1)%val),
     &                             R_ptr, ptr2_shape)
                  CALL sgemm('T', 'N', N, NRHS_B, Kprec, MONE,
     &                  R_ptr(1,1), Kprec, 
     &                  TEMP_BLOCK(kcumul), K, ONE,
     &                  DEST_ARRAY(1), N)
                end if 
                kprec = lrb%k - lrb%k_stoc(1)
                kcumul = kcumul + lrb%k_stoc(1)
                overall_k = lrb%k_stoc(1)
                total_consumption = lrb%n*( kprec * costs_prec(1) + 
     &                                  lrb%kw_stoc ) / 
     &                              dble((4.0/5.0)*cache_L2)
C                if (kprec.gt.0) then
C                  R_ptr(1:kprec,1:N) => 
C     &                thread_workspaces(1:(N*kprec))
C                  CALL SMUMPS_DECONVERT_LRB_R_W(LRB, R_ptr, 
C     &                kprec, N,
C     &                2, NPREC_STOC, 1, LRB%K_STOC(NPREC_STOC),
C     &                IFLAG, IERROR, KEEP8)
C                  CALL sgemm('T', 'N', N, NRHS_B, kprec, MONE,
C     &                R_ptr(1,1), kprec, 
C     &                TEMP_BLOCK(kcumul), K, ONE,
C     &                DEST_ARRAY(1), N)
C                end if
                start_prec = 2
                budget = 0
                budget_ws = 0
                start_k = 1
                current_k = 0
                overall_k = lrb%k_stoc(1)
                nu = 1.0d0
                if (k455.eq.4) nu = 0.2d0
                if (k455.eq.5.and.lrb%N.gt.200) nu = 0.2d0
                current_consumption = 0.0d0
                overall_consumption = 0.0d0
                do prec = 2, nprec_stoc
                  if (lrb%k_stoc(prec).eq.0) cycle
                  krem = lrb%k_stoc(prec)
                  current_prec = 1
                  do
                    ! avoid small blas when cache is barely over filled
                    if (k455.eq.3.and.nu.eq.1.0d0
     &          .and. total_consumption - overall_consumption.gt.1.0D0
     &          .and. total_consumption - overall_consumption.lt.1.5D0
     &              ) then
                      nu = (total_consumption - overall_consumption)/2.0
                      nu = nu + 0.02
                    end if
                    kmagic = max(dble(nu*(4.0/5.0)*CACHE_L2
     &                        - budget - budget_ws),0.0d0) /
     &                       (lrb%N*(costs_prec(prec)+costs_prec(1)))
C                    if (kmagic.lt.krem.and.kmagic.gt.8)
C     &                kmagic = kmagic - mod(kmagic,8)
                    used_k = min(kmagic, krem)
                    ws_footprint = lrb%n * used_k * costs_prec(1)
                    prec_footprint = lrb%n * used_k * costs_prec(prec)
                    budget = budget + prec_footprint
                    budget_ws = budget_ws + ws_footprint
                    current_consumption = current_consumption +
     &                dble(prec_footprint+ws_footprint) / 
     &                ((80.0/100.0)*dble(cache_L2))
                    overall_consumption = overall_consumption +
     &                dble(prec_footprint+ws_footprint) / 
     &                ((80.0/100.0)*dble(cache_L2))
                    current_prec = current_prec + used_k
                    current_k = current_k + used_k
                    krem = krem - used_k
                    overall_k = overall_k + used_k
                    if (krem.ge.1.and.k455.eq.4
     &                  .and.current_consumption.gt.0.3.or.
     &                  krem.ge.1.or.prec.eq.nprec_stoc.or.
     &                  overall_k.eq.k) then
                      R_ptr(1:current_k,1:N) => 
     &                    thread_workspaces(1:(N*current_k))
                      end_prec = prec
                      end_k = lrb%k_stoc(prec) - krem
                      ! if we switched prec but even a row wont fit
                      if (krem.eq.lrb%k_stoc(prec)) then
                        end_prec = prec - 1
                        end_k = lrb%k_stoc(end_prec)
                      end if
                      CALL SMUMPS_DECONVERT_LRB_R_W(LRB, R_ptr, 
     &                    current_k, N,
     &                    start_prec, end_prec, 
     &                    start_k, end_k,
     &                    IFLAG, IERROR, KEEP8)
C                      write(*,*) "MYSTUFF", n, current_k, k
                      CALL sgemm('T', 'N', N, NRHS_B, current_k, MONE,
     &                    R_ptr(1,1), current_k, 
     &                    TEMP_BLOCK(kcumul), K, ONE,
     &                    DEST_ARRAY(1), N)
                      kcumul = kcumul + current_k
                      start_prec = prec
                      budget = 0
                      budget_ws = 0
                      start_k = current_prec
                      current_k = 0
                      current_consumption = 0.0d0
                    end if
                    if (krem.lt.1) then 
                      exit
                    end if
                  end do
                end do
              else
                kcumul = lrb%k + 1
                start_prec=nprec_stoc
                ! if the sum of R^T (except working prec.) fit in cache
                kprec = lrb%k - lrb%k_stoc(1)
                if (kprec.gt.0.and.
     &              (k455.eq.0.or.k455.ne.1.and.
     &              N*LRB%kw_stoc +
     &              N*costs_prec(1) * kprec < (4.0/5)*cache_L2)) then
                  kcumul = kcumul - kprec
                  R_ptr(1:kprec,1:N) => 
     &                thread_workspaces(1:(N*Kprec))
                  CALL SMUMPS_DECONVERT_LRB_R_W(LRB, 
     &               R_ptr, kprec, N,
     &                2, NPREC_STOC, 1, LRB%K_STOC(NPREC_STOC),
     &                IFLAG, IERROR, KEEP8)
                  CALL sgemm('T', 'N', N, NRHS_B, Kprec, MONE,
     &                  R_ptr(1,1), Kprec, 
     &                  TEMP_BLOCK(kcumul), K, ONE,
     &                  DEST_ARRAY(1), N)
                  ! working prec is adressed in following do loop
                  start_prec=1
                endif
                ! else (e.g. if any R^T_prec is too big)
                do prec=start_prec,1,-1
                   kprec = lrb%k_stoc(prec)
                   if (kprec.eq.0) cycle
                   kcumul = kcumul - kprec
                   if (prec.eq.1) then
                     KMAGIC=Kprec
                     ptr2_shape(1) = Kprec
                     ptr2_shape(2) = LRB%N
                     call c_f_pointer(c_loc(LRB%R_stoc(prec)%val),
     &                                R_ptr, ptr2_shape)
                   else
                     ptr1_shape(1) = costs_prec(prec)*Kprec*LRB%N
                     call c_f_pointer(c_loc(LRB%R_stoc(prec)%val),
     &                    stoc_ptr, ptr1_shape)
                     KMAGIC=(80.0/100.0)*(CACHE_L2-costs_prec(1)*N)/
     &               ( (costs_prec(1)+costs_prec(prec))*N + 
     &                 NRHS_B*costs_prec(1) )
                     KMAGIC = KMAGIC - MODULO(KMAGIC,8)
                   endif
                   if (k455.eq.1) kmagic = kprec
                   do kint=1,Kprec,KMAGIC
                     !Krem may be stupidly small
                     krem = min(KMAGIC,Kprec-Kint+1)
                     if (prec.gt.1) then
                       R_ptr(1:Krem,1:LRB%N) => 
     &                     thread_workspaces(1:(LRB%N*krem))
                       CALL SMUMPS_LRMIXED_CAST(
     &                     stoc_ptr(costs_prec(prec)*(kint-1)+1),
     &                     R_ptr(1,1), Krem, LRB%N, 
     &                     Kprec, Krem, PREC, 1)
                     endif
                     CALL sgemm('T', 'N', N, NRHS_B, Krem, MONE,
     &                 R_ptr(1,1), Krem, 
     &                 TEMP_BLOCK(kcumul+kint-1), K, ONE,
     &                 DEST_ARRAY(1), N)
                  end do
                end do
             end if
            ELSE
              ptpt_rhs(1:nprec_calc) => pt_rhs(1:nprec_calc,I)
              IF (IS_T2_SLAVE) THEN
                POSBLOCK = POSCB+int(IBEG_BLOCK-1,8)
                CALL SMUMPS_LRMIXED_LRxFR(1, LRB,
     &            ARRAYCB(POSBLOCK), LDCB,
     &            DEST_ARRAY(1), N, NRHS_B, -1.0D0,
     &            K34, K451, KEEP8, ptPT_RHS(1), PT_RB(1), PT_QRB,
     &            thread_workspaces, wsmax)
              ELSEIF (IBEG_BLOCK.LE.NPIV) THEN
C               Entirely on FS side              
                POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
                CALL SMUMPS_LRMIXED_LRxFR(1, LRB,
     &            ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV,
     &            DEST_ARRAY(1), N,
     &            NRHS_B, -1.0D0, K34, K451, KEEP8, 
     &            ptPT_RHS(1), PT_RB(1), PT_QRB,
     &            thread_workspaces, wsmax)
              ELSE
C               Entirely on CB side              
                POSBLOCK = POSCB+int(IBEG_BLOCK-1-NPIV,8)
                CALL SMUMPS_LRMIXED_LRxFR(1, LRB,
     &            ARRAYCB(POSBLOCK), LDCB,
     &            DEST_ARRAY(1), N, NRHS_B, -1.0D0,
     &            K34, K451, KEEP8, ptPT_RHS(1), PT_RB(1), PT_QRB,
     &            thread_workspaces, wsmax)
              ENDIF
            ENDIF
! LRMIXED
          ELSE ! NOT ISLR
            CALL FRB_GET_PREC(LRB, PREC, MIXED_STORED)
            IF (K448.GE.1.AND.K450.NE.0) THEN
              IF (IS_T2_SLAVE) THEN
                POSBLOCK = POSCB + int(IBEG_BLOCK-1,8)
                CALL SMUMPS_LRMIXED_FRxFR('T', 'N',
     &                N, NRHS_B, M, MONE, 
     &                LRB, 1, 1, M,
     &                ARRAYCB(POSBLOCK), LDCB,
     &                DEST_ARRAY(1), N,
     &                PT_RHS(1,I), TEMP_QRB, thread_workspaces, WSMAX)
              ELSEIF (IBEG_BLOCK.LE.NPIV.AND.IEND_BLOCK.GT.NPIV) THEN
C               Case of overlap between the FS and non FS
                POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
                CALL SMUMPS_LRMIXED_FRxFR('T', 'N',
     &                N, NRHS_B, NPIV-IBEG_BLOCK+1, MONE, 
     &                LRB, 1, 1, M,
     &                ARRAYPIV(POSBLOCK,POSPIVCOL),LDPIV,
     &                DEST_ARRAY(1), N,
     &                PT_RHS(1,I), TEMP_QRB, thread_workspaces, WSMAX)
                CALL SMUMPS_LRMIXED_FRxFR('T', 'N',
     &                N, NRHS_B, M-(NPIV-IBEG_BLOCK+1), MONE, 
     &                LRB, NPIV-IBEG_BLOCK+2, 1, M,
     &                ARRAYCB(POSCB), LDCB,
     &                DEST_ARRAY(1), N,
     &                PT_RHS(1,I), TEMP_QRB, thread_workspaces, WSMAX)
              ELSEIF (IBEG_BLOCK.LE.NPIV) THEN
C               Entirely on FS side              
                POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
                CALL SMUMPS_LRMIXED_FRxFR('T', 'N',
     &                N, NRHS_B, M, MONE, 
     &                LRB, 1, 1, M,
     &                ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV,
     &                DEST_ARRAY(1), N,
     &                PT_RHS(1,I), TEMP_QRB, thread_workspaces, WSMAX)
              ELSE
C               Entirely on CB side              
                POSBLOCK = POSCB + int(IBEG_BLOCK-1-NPIV,8)
                CALL SMUMPS_LRMIXED_FRxFR('T', 'N',
     &                N, NRHS_B, M, MONE, 
     &                LRB, 1, 1, M,
     &                ARRAYCB(POSBLOCK),LDCB,
     &                DEST_ARRAY(1), N,
     &                PT_RHS(1,I), TEMP_QRB, thread_workspaces, WSMAX)
              ENDIF
            ELSE
              if (prec.le.1.and.associated(LRB%Q)) then
                KMAGIC = N
                ptr2_shape(1) = M
                ptr2_shape(2) = KMAGIC
                call c_f_pointer(c_loc(LRB%Q),
     &                           Q_ptr, ptr2_shape)
              else
                KMAGIC = (80.0/100.0)*CACHE_L2/
     &            ((costs_prec(1)+costs_prec(prec))*M)
                KMAGIC=MIN(KMAGIC,N)
                Q_ptr(1:M,1:KMAGIC) =>
     &            thread_workspaces(1:(M*KMAGIC))
                ptr1_shape(1) = costs_prec(prec)*M*N
                call c_f_pointer(c_loc(LRB%Q_stoc(prec)%val),
     &            stoc_ptr, ptr1_shape)
              ENDIF
              DO KINT=1,N,KMAGIC
                krem = min(KMAGIC,N-Kint+1)
                if (prec.gt.1) then
                  CALL SMUMPS_LRMIXED_CAST(
     &              stoc_ptr(costs_prec(prec)*M*(kint-1)+1),
     &              Q_ptr(1,1), M, KREM,
     &              M, M, PREC, 1)
                ENDIF
!               CALL GEMM_Q_BWD(M, NRHS_B, krem, NPIV, 
!     &              Q_ptr, DEST_ARRAY, N, 
!     &              ARRAYPIV, ldpiv, ARRAYCB, lcb, ldcb, 
!     &              ibeg_block, iend_block, is_t2_slave, 
!     &              poscb, pospiv, pospivcol, kint)
! LRMIXED
                IF (IS_T2_SLAVE) THEN
                  POSBLOCK = POSCB+int(IBEG_BLOCK-1,8)
                  CALL sgemm('T', 'N', KREM, NRHS_B, M, MONE,
     &                  Q_ptr(1,1), M, ARRAYCB(POSBLOCK), LDCB,
     &                  ONE, DEST_ARRAY(KINT), N)
                ELSEIF (IBEG_BLOCK.LE.NPIV.AND.IEND_BLOCK.GT.NPIV) THEN
C                 Case of overlap between the FS and non FS
                  POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
                  CALL sgemm('T', 'N', KREM, NRHS_B, NPIV-IBEG_BLOCK+1,
     &                  MONE,Q_ptr(1,1), M, 
     &                  ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV, 
     &                  ONE, DEST_ARRAY(KINT), N)
                  CALL sgemm('T', 'N', KREM, NRHS_B,IBEG_BLOCK+M-NPIV-1,
     &                  MONE,Q_ptr(NPIV-IBEG_BLOCK+2,1),M,
     &                  ARRAYCB(POSCB), LDCB, 
     &                  ONE, DEST_ARRAY(KINT), N)
                ELSEIF (IBEG_BLOCK.LE.NPIV) THEN
C                 Entirely on FS side              
                  POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
                  CALL sgemm('T', 'N', KREM, NRHS_B, M, MONE,
     &                  Q_ptr(1,1), M, 
     &                  ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV,
     &                  ONE, DEST_ARRAY(KINT), N)
                ELSE
C                 Entirely on CB side              
                  POSBLOCK = POSCB+int(IBEG_BLOCK-1-NPIV,8)
                  CALL sgemm('T', 'N', KREM, NRHS_B, M, MONE,
     &                  Q_ptr(1,1), M,  ARRAYCB(POSBLOCK), LDCB,
     &                  ONE, DEST_ARRAY(KINT), N)
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ENDDO
#if ! defined(BLR_NOOPENMP)
!$OMP END DO
#endif
       IF (allocated(TEMP_BLOCK)) deallocate(TEMP_BLOCK)
        IF (allocated(PT_QRB))   deallocate(PT_QRB)
        IF (allocated(TEMP_RB))  deallocate(TEMP_RB)
        IF (allocated(TEMP_QRB)) deallocate(TEMP_QRB)
        IF (allocated(val_workspace)) deallocate(val_workspace)
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif 
        IF (IS_T2_SLAVE) THEN
          DO I=1,NRHS_B
            call saxpy(N, ONE, DEST_ARRAY((I-1)*N+1), 1,
     &                 ARRAYPIV(POSDIAG+(I-1)*LDPIV,POSPIVCOL), 1)
          ENDDO
        ELSE
          DO I=1,NRHS_B
            call saxpy(N, ONE, DEST_ARRAY((I-1)*N+1), 1,
     &                 ARRAYPIV(POSDIAG,POSPIVCOL+I-1), 1)
          ENDDO
        ENDIF
 100    CONTINUE
        IF (allocated(DEST_ARRAY)) DEALLOCATE(DEST_ARRAY)
        RETURN
      END SUBROUTINE SMUMPS_SOL_BWD_BLR_UPDATE
!MIXHPC
#else 
!.not.MIXHPC
      SUBROUTINE SMUMPS_SOL_BWD_BLR_UPDATE (
     &        ARRAYPIV, LPIV, LPIVCOL, LDPIV, POSPIV, POSPIVCOL,
     &        ARRAYCB, LCB, LDCB, POSCB, 
     &        POSDIAG,
     &        NRHS_B, NPIV,
     &        BLR_PANEL, LAST_BLR, CURRENT_BLR,
     &        BEGS_BLR_STATIC,
     &        KEEP8, K34, K448, K450, K451, IS_T2_SLAVE,
     &        PT_RHS, 
     &        IFLAG, IERROR)
C
C Inspired from SMUMPS_BLR_UPD_NELIM_VAR_L
C     
C     Update of the solution phase of current panel 
C     -----------
C     Parameters
C     -----------
C
C     ARRAYPIV is the array containing the FS (eliminated) rows of the solution (WCB
C     or RHSINTR depending on K350 in the future). It is of size LPIV
C     and leading dimension LDPIV. POSPIV points to the first eliminated row
C     associated to the current front in ARRAYPIV and is not used for
C     type 2 slaves
C     ARRAYCB is the array containing the NFS+NELIM rows of the solution
C     (always WCB). It is of size LCB and leading dimension LDCB. 
C     POSCB points to the first non-eliminated row associated to the current front
C     in ARRAYCB. 
C     POSDIAG points to the beginning of the current panel CURRENT_BLR
C     (in ARRAYPIV).        
C     NPIV is the global NPIV of the front
C     LPIVCOL is either 1 (ARRAYPIV=WCB) or NRHS (ARRAYPIV=RHSINTR)        
      INTEGER(8), INTENT(IN) :: LPIV, LCB, POSPIV, POSCB, POSDIAG
      INTEGER,INTENT(IN)     :: LPIVCOL, POSPIVCOL
      REAL, INTENT(INOUT) :: ARRAYPIV(LPIV,LPIVCOL)
      REAL, INTENT(INOUT) :: ARRAYCB(LCB)
      INTEGER, INTENT(IN)    :: LAST_BLR, NRHS_B,  LDPIV, LDCB,
     &                          CURRENT_BLR, NPIV, K34, K448, K450, K451
c     TYPE(LRB_TYPE), TARGET,INTENT(IN) ::
c    &                          BLR_PANEL(LAST_BLR-CURRENT_BLR)
      TYPE(LRB_TYPE), TARGET,INTENT(IN) ::
     &                          BLR_PANEL(:)
      INTEGER(8), INTENT(IN) :: KEEP8(150)
      LOGICAL, INTENT(IN)    :: IS_T2_SLAVE
      INTEGER :: BEGS_BLR_STATIC(:)
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
#if defined(MUMPS_NOF2003)
C     Argument cannot be allocatable.
C     Even if not used, PT_RHS has been allocated of size (1,1) outside
      TYPE(LRMIXED_TYPE) :: PT_RHS(:,:)
#else
      TYPE(LRMIXED_TYPE), ALLOCATABLE :: PT_RHS(:,:)
#endif
C     -----------
C     Local variables
C     -----------
      INTEGER :: I, K, M, N, IBEG_BLOCK, IEND_BLOCK
      INTEGER :: KMAX
      INTEGER(8) :: POSBLOCK
      TYPE(LRB_TYPE), POINTER :: LRB
      REAL, ALLOCATABLE, DIMENSION(:) :: TEMP_BLOCK
      REAL, ALLOCATABLE, DIMENSION(:) :: DEST_ARRAY
      INTEGER :: allocok
      REAL :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0E0, MONE=-1.0E0)
      PARAMETER (ZERO=0.0E0)
      INTEGER, ALLOCATABLE, TARGET, DIMENSION(:) :: TEMP_RB, TEMP_QRB
      INTEGER :: POS_QRB
      TYPE(LRMIXED_TYPE), ALLOCATABLE :: PT_QRB(:)
      INTEGER, POINTER:: PT_RB(:)
      INTEGER :: PREC, NB_INT_RB, NB_INT_QRB
      LOGICAL :: A_CHEVAL
      KMAX = -1
      DO I = CURRENT_BLR+1, LAST_BLR
        KMAX = max(KMAX, BLR_PANEL(I-CURRENT_BLR)%K)
      ENDDO
      IF (CURRENT_BLR.LT.LAST_BLR) THEN 
C       N is constant for all blocks in the panel        
        N = BLR_PANEL(1)%N
      ELSE
        RETURN
      ENDIF
      allocate(DEST_ARRAY(N*NRHS_B),stat=allocok)
      IF (allocok .GT. 0) THEN
        IFLAG = -13
        IERROR = N * NRHS_B
        GOTO 100
      ENDIF
      DEST_ARRAY = ZERO
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL PRIVATE(TEMP_BLOCK, allocok
!$OMP& , TEMP_QRB, TEMP_RB, 
!$OMP& PT_QRB, PT_RB, POS_QRB, PREC, NB_INT_RB, NB_INT_QRB
!$OMP& )
#endif 
      IF (K448.GE.1.AND.K450.GT.0) THEN
C       Allocate RB and QRB workspaces for each thread        
        IF(CURRENT_BLR.LT.LAST_BLR) THEN
          NB_INT_QRB = 0
          DO PREC=2,NPREC_CALC
            NB_INT_QRB = NB_INT_QRB + 
     &         (N*NRHS_B*COSTS_PREC(CALC2STOC(PREC)) +K34-1) / K34
          ENDDO
          NB_INT_RB=(KMAX*NRHS_B*COSTS_PREC(CALC2STOC(1)) + K34-1) / K34
          allocate(TEMP_RB(NB_INT_RB))
          allocate(TEMP_QRB(NB_INT_QRB))
          allocate(PT_QRB(NPREC_CALC))
          nullify(PT_QRB(1)%val)
          POS_QRB = 1
          DO PREC=2,NPREC_CALC
            PT_QRB(PREC)%VAL => TEMP_QRB(POS_QRB:)
            POS_QRB = POS_QRB + 
     &         (N*NRHS_B*COSTS_PREC(CALC2STOC(PREC)) +K34-1) / K34
          ENDDO
          PT_RB=>TEMP_RB(1:NB_INT_RB)
        ENDIF
      ENDIF
      IF (KMAX.GT.0) THEN
C       TEMP_BLOCK is allocated even when K450.GT.0 because we might
C       need it when A_CHEVAL is true        
        allocate(TEMP_BLOCK(KMAX*NRHS_B), stat=allocok )
        IF (allocok .GT. 0) THEN
          IFLAG  = -13
          IERROR = NRHS_B * KMAX
          write(*,*) 'Allocation problem in BLR routine 
     &    SMUMPS_SOL_BWD_BLR_UPDATE: ',
     &    'not enough memory? memory requested = ', IERROR
        ENDIF
      ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP DO SCHEDULE(DYNAMIC,1)
!$OMP& PRIVATE(IBEG_BLOCK, IEND_BLOCK, LRB, K, M, POSBLOCK)
!$OMP& PRIVATE(A_CHEVAL)
!$OMP& REDUCTION(+:DEST_ARRAY)
#endif
        DO I = CURRENT_BLR+1, LAST_BLR
          IF (IFLAG.LT.0) CYCLE
          IBEG_BLOCK = BEGS_BLR_STATIC(I)
          IEND_BLOCK = BEGS_BLR_STATIC(I+1)-1
          LRB => BLR_PANEL(I-CURRENT_BLR)
C Case of U: U^T = Q R
          K = LRB%K 
          M = LRB%M   ! always nb of rows of Q
C
          IF (LRB%ISLR) THEN
            IF (K.GT.0) THEN
C             Check if block is "a cheval": if so, do not use LRMIXED
C             for simplicity (we could try to code this in the future)
              A_CHEVAL =  (.NOT.IS_T2_SLAVE) .AND. (IBEG_BLOCK.LE.NPIV) 
     &              .AND. (IEND_BLOCK.GT.NPIV)
              IF (K450.GT.0.AND.(K448.EQ.0.OR.A_CHEVAL)) THEN
                CALL SMUMPS_DECONVERT_LRB(LRB, IFLAG, IERROR, KEEP8)
              ENDIF
              IF (K448.EQ.0.OR.K450.EQ.0 .OR. A_CHEVAL) THEN
                IF (IS_T2_SLAVE) THEN
                  POSBLOCK = POSCB +int(IBEG_BLOCK-1,8)
                  CALL sgemm('T', 'N', K, NRHS_B, M, ONE,
     &              LRB%Q(1,1), M,  
     &              ARRAYCB(POSBLOCK), LDCB, ZERO, 
     &              TEMP_BLOCK(1), K)
                ELSEIF (IBEG_BLOCK.LE.NPIV.AND.IEND_BLOCK.GT.NPIV) THEN
C                 Case of overlap between the FS and non FS
                  POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
                  CALL sgemm('T', 'N', K, NRHS_B, NPIV-IBEG_BLOCK+1,
     &              ONE, LRB%Q(1,1), M,  
     &              ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV,
     &              ZERO, TEMP_BLOCK(1), K)
                  CALL sgemm('T', 'N', K, NRHS_B, IBEG_BLOCK+M-NPIV-1,
     &              ONE, LRB%Q(NPIV-IBEG_BLOCK+2,1), M,  
     &              ARRAYCB(POSCB), LDCB, 
     &              ONE,  ! to sum to the result of the previous product
     &              TEMP_BLOCK(1), K)
                ELSEIF (IBEG_BLOCK.LE.NPIV) THEN
C                 Entirely on FS side              
                  POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
                  CALL sgemm('T', 'N', K, NRHS_B, M, ONE,
     &              LRB%Q(1,1), M,  
     &              ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV,
     &              ZERO, TEMP_BLOCK(1), K)
                ELSE
C                 Entirely on CB side              
                  POSBLOCK = POSCB+int(IBEG_BLOCK-1-NPIV,8)
                  CALL sgemm('T', 'N', K, NRHS_B, M, ONE,
     &              LRB%Q(1,1), M,  
     &              ARRAYCB(POSBLOCK), LDCB, ZERO, 
     &              TEMP_BLOCK(1), K)
                ENDIF
                CALL sgemm('T', 'N', N, NRHS_B, K, MONE,
     &            LRB%R(1,1), K, 
     &            TEMP_BLOCK(1), K, ONE,
     &            DEST_ARRAY(1), N)
              ELSE
                IF (IS_T2_SLAVE) THEN
                  POSBLOCK = POSCB+int(IBEG_BLOCK-1,8)
                  CALL SMUMPS_LRMIXED_LRxFR(1, LRB,
     &              ARRAYCB(POSBLOCK), LDCB,
     &              DEST_ARRAY(1), N, NRHS_B, -1.0D0,
     &              K34, K451, KEEP8, PT_RHS(1,I), PT_RB(1), PT_QRB)
                ELSEIF (IBEG_BLOCK.LE.NPIV) THEN
C                 Entirely on FS side              
                  POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
                  CALL SMUMPS_LRMIXED_LRxFR(1, LRB,
     &              ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV,
     &              DEST_ARRAY(1), N,
     &              NRHS_B, -1.0D0, K34, K451, KEEP8, 
     &              PT_RHS(1,I), PT_RB(1), PT_QRB)
                ELSE
C                 Entirely on CB side              
                  POSBLOCK = POSCB+int(IBEG_BLOCK-1-NPIV,8)
                  CALL SMUMPS_LRMIXED_LRxFR(1, LRB,
     &              ARRAYCB(POSBLOCK), LDCB,
     &              DEST_ARRAY(1), N, NRHS_B, -1.0D0,
     &              K34, K451, KEEP8, PT_RHS(1,I), PT_RB(1), PT_QRB)
                ENDIF
              ENDIF
            ENDIF ! K GT 0
          ELSE ! NOT ISLR
            IF (IS_T2_SLAVE) THEN
              POSBLOCK = POSCB+int(IBEG_BLOCK-1,8)
              CALL sgemm('T', 'N', N, NRHS_B, M, MONE,
     &              LRB%Q(1,1), M,  ARRAYCB(POSBLOCK), LDCB,
     &              ONE, DEST_ARRAY(1), N)
            ELSEIF (IBEG_BLOCK.LE.NPIV.AND.IEND_BLOCK.GT.NPIV) THEN
C             Case of overlap between the FS and non FS
              POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
              CALL sgemm('T', 'N', N, NRHS_B, NPIV-IBEG_BLOCK+1, MONE,
     &              LRB%Q(1,1), M,  ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV,
     &              ONE, DEST_ARRAY(1), N)
              CALL sgemm('T', 'N', N, NRHS_B, IBEG_BLOCK+M-NPIV-1, MONE,
     &              LRB%Q(NPIV-IBEG_BLOCK+2,1), M,  ARRAYCB(POSCB), 
     &              LDCB, ONE, DEST_ARRAY(1), N)
            ELSEIF (IBEG_BLOCK.LE.NPIV) THEN
C             Entirely on FS side              
              POSBLOCK = POSPIV+int(IBEG_BLOCK-1,8)
              CALL sgemm('T', 'N', N, NRHS_B, M, MONE,
     &              LRB%Q(1,1), M,  ARRAYPIV(POSBLOCK,POSPIVCOL), LDPIV,
     &              ONE, DEST_ARRAY(1), N)
            ELSE
C             Entirely on CB side              
              POSBLOCK = POSCB+int(IBEG_BLOCK-1-NPIV,8)
              CALL sgemm('T', 'N', N, NRHS_B, M, MONE,
     &              LRB%Q(1,1), M,  ARRAYCB(POSBLOCK), LDCB,
     &              ONE, DEST_ARRAY(1), N)
            ENDIF
          ENDIF
          IF (LRB%ISLR .AND. K450.GT.0
     &        .AND.(K448.EQ.0.OR. A_CHEVAL)
     &       ) THEN
C           Deallocate UNIFORM part, keep MIXED for future solves            
            CALL DEALLOC_LRB(LRB, KEEP8, K34, .FALSE., .TRUE.)
          ENDIF
        ENDDO
#if ! defined(BLR_NOOPENMP)
!$OMP END DO
#endif
        IF (KMAX.GT.0) THEN
          IF (allocated(TEMP_BLOCK)) deallocate(TEMP_BLOCK)
        ENDIF
        IF (allocated(PT_QRB))   deallocate(PT_QRB)
        IF (allocated(TEMP_RB))  deallocate(TEMP_RB)
        IF (allocated(TEMP_QRB)) deallocate(TEMP_QRB)
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif 
        IF (IS_T2_SLAVE) THEN
          DO I=1,NRHS_B
            call saxpy(N, ONE, DEST_ARRAY((I-1)*N+1), 1,
     &                 ARRAYPIV(POSDIAG+(I-1)*LDPIV,POSPIVCOL), 1)
          ENDDO
        ELSE
          DO I=1,NRHS_B
            call saxpy(N, ONE, DEST_ARRAY((I-1)*N+1), 1,
     &                 ARRAYPIV(POSDIAG,POSPIVCOL+I-1), 1)
          ENDDO
        ENDIF
 100    CONTINUE
        IF (allocated(DEST_ARRAY)) DEALLOCATE(DEST_ARRAY)
        RETURN
      END SUBROUTINE SMUMPS_SOL_BWD_BLR_UPDATE
!MIXHPC
#endif 
      SUBROUTINE COMPUTE_WORKSPACE_SIZE(LRB, WSMAX,
     &      K34, K148, K448, K450, K455, CACHE, NRHS_B)
         IMPLICIT NONE
         TYPE(LRB_TYPE), intent(in) :: LRB
         INTEGER, intent(inout) :: WSMAX
         INTEGER, intent(in) :: K34, K148, K448, K450, K455
         INTEGER, intent(in) :: CACHE
         INTEGER, intent(in) :: NRHS_B
         ! Local variables
         INTEGER :: PREC, KPREC
         INTEGER :: M, N, K
         INTEGER :: KMAGIC, KREM
         ! for cache stuff 
         DOUBLE PRECISION :: TOTAL_CONSUMPTION, CURRENT_CONSUMPTION
         DOUBLE PRECISION :: OVERALL_CONSUMPTION
         INTEGER :: BUDGET, BUDGET_WS
         INTEGER :: WS_FOOTPRINT, PREC_FOOTPRINT
         INTEGER :: CURRENT_K, USED_K
         REAL :: NU
         INTEGER :: OVERALL_K
         INTEGER :: PREC_CALC, PCALC_STOC
         M = LRB%M
         N = LRB%N
         ! AJ is memory that tight ?
         IF (LRB%ISLR) THEN
           IF(LRB%K.LE.0) RETURN
           K = LRB%K
           KPREC = LRB%K - LRB%K_STOC(1)
           SELECT CASE(K455)
             CASE (:0) ! Aggregate all non-arith precs
               WSMAX = MAX(WSMAX, M * KPREC)
               WSMAX = MAX(WSMAX, N * KPREC)
             CASE (1)  ! Separate all precs
               DO PREC=2,NPREC_STOC
                 WSMAX = MAX(WSMAX, M * LRB%K_STOC(PREC))
                 WSMAX = MAX(WSMAX, N * LRB%K_STOC(PREC))
               END DO
             CASE (2) ! Aggregate if possible, otherwise separate
               IF (N*LRB%KW_STOC + N*COSTS_PREC(1)*KPREC <
     &           (80.0/100.0)*CACHE) THEN
                 WSMAX = MAX(WSMAX, KPREC * N)
               END IF
               IF (M*LRB%KW_STOC + M*COSTS_PREC(1)*KPREC < 
     &           (80.0/100.0)*CACHE) THEN
                 WSMAX = MAX(WSMAX, KPREC * M)
               END IF
               DO PREC=2,NPREC_STOC
                 K = LRB%K_STOC(PREC)
                 IF (K.LE.0) CYCLE
                 KMAGIC=(80.0/100.0)*(CACHE-COSTS_PREC(1)*N)/
     &           ( (COSTS_PREC(1)+COSTS_PREC(PREC))*N + 
     &             NRHS_B*COSTS_PREC(1) )
                 WSMAX = MAX(WSMAX, KMAGIC * N)
                 KMAGIC=(80.0/100.0)*(CACHE-COSTS_PREC(1)*M)/
     &           ( (COSTS_PREC(1)+COSTS_PREC(PREC))*M + 
     &             NRHS_B*COSTS_PREC(1) )
                 WSMAX = MAX(WSMAX, KMAGIC * M)
               END DO
             CASE (3:) ! FIT THE MOST IN CACHE
               KPREC = LRB%K - LRB%K_STOC(1)
               BUDGET = 0
               BUDGET_WS = 0
               NU = 1.0D0
               IF (K455.EQ.4) NU = 0.2D0
               IF (K455.EQ.5.AND.M.GT.200) NU = 0.2D0
               CURRENT_CONSUMPTION = 0.0D0
               OVERALL_CONSUMPTION = 0.0D0
               TOTAL_CONSUMPTION = M*( KPREC * COSTS_PREC(1) + 
     &                                 LRB%KW_STOC ) / 
     &                             ((4.0/5.0)*DBLE(CACHE))
               CURRENT_K = 0
               OVERALL_K = LRB%K_STOC(1)
               DO PREC=2,NPREC_STOC
                 IF (LRB%K_STOC(PREC).EQ.0) CYCLE
                 KREM = LRB%K_STOC(PREC) 
                 DO
                   ! AVOID SMALL BLAS WHEN CACHE IS BARELY OVER FILLED
                   IF (K455.EQ.3.AND.NU.EQ.1.0D0
     &         .AND. TOTAL_CONSUMPTION - OVERALL_CONSUMPTION.GT.1.0D0
     &         .AND. TOTAL_CONSUMPTION - OVERALL_CONSUMPTION.LT.1.5D0
     &             ) THEN
                     NU = (TOTAL_CONSUMPTION - OVERALL_CONSUMPTION)/2.0
                     NU = NU + 0.02
                   END IF
                   KMAGIC = MAX(DBLE(NU*(4.0/5.0)*CACHE
     &                       - BUDGET - BUDGET_WS),0.0D0) /
     &                      (M*(COSTS_PREC(PREC)+COSTS_PREC(1)))
                   USED_K = MIN(KMAGIC, KREM)
                   WS_FOOTPRINT = M * USED_K * COSTS_PREC(1)
                   PREC_FOOTPRINT = M * USED_K * COSTS_PREC(PREC)
                   BUDGET = BUDGET + PREC_FOOTPRINT
                   BUDGET_WS = BUDGET_WS + WS_FOOTPRINT
                   CURRENT_CONSUMPTION = CURRENT_CONSUMPTION +
     &               DBLE(PREC_FOOTPRINT+WS_FOOTPRINT) / 
     &               ((80.0/100.0)*DBLE(CACHE))
                   OVERALL_CONSUMPTION = OVERALL_CONSUMPTION +
     &               DBLE(PREC_FOOTPRINT+WS_FOOTPRINT) / 
     &               ((80.0/100.0)*DBLE(CACHE))
                   KREM = KREM - USED_K
                   OVERALL_K = OVERALL_K + USED_K
                   CURRENT_K = CURRENT_K + USED_K
                   IF (KREM.GE.1.AND.K455.EQ.4
     &                 .AND.CURRENT_CONSUMPTION.GT.0.3.OR.
     &                 KREM.GE.1.OR.PREC.EQ.NPREC_STOC.OR.
     &                 OVERALL_K.EQ.LRB%K) THEN
                      WSMAX = MAX(WSMAX, M*CURRENT_K)
                      BUDGET = 0
                      BUDGET_WS = 0
                      CURRENT_K = 0
                      CURRENT_CONSUMPTION = 0.0D0
                   ENDIF
                   IF (KREM.LT.1) EXIT
                 END DO
               END DO
               IF (N.NE.M) THEN
                 BUDGET = 0
                 BUDGET_WS = 0
                 NU = 1.0D0
                 IF (K455.EQ.4) NU = 0.2D0
                 IF (K455.EQ.5.AND.N.GT.200) NU = 0.2D0
                 CURRENT_CONSUMPTION = 0.0D0
                 OVERALL_CONSUMPTION = 0.0D0
                 TOTAL_CONSUMPTION = N*( KPREC * COSTS_PREC(1) + 
     &                                   LRB%KW_STOC ) / 
     &                               ((4.0/5.0)*DBLE(CACHE))
                 CURRENT_K = 0
                 OVERALL_K = LRB%K_STOC(1)
                 DO PREC=2,NPREC_STOC
                   IF (LRB%K_STOC(PREC).EQ.0) CYCLE
                   KREM = LRB%K_STOC(PREC) 
                   DO
                     ! AVOID SMALL BLAS WHEN CACHE IS BARELY OVER FILLED
                      IF (K455.EQ.3.AND.NU.EQ.1.0D0
     &           .AND. TOTAL_CONSUMPTION - OVERALL_CONSUMPTION.GT.1.0D0
     &           .AND. TOTAL_CONSUMPTION - OVERALL_CONSUMPTION.LT.1.5D0
     &               ) THEN
                       NU = (TOTAL_CONSUMPTION - OVERALL_CONSUMPTION)
     &                      /2.0
                       NU = NU + 0.02
                     END IF
                     KMAGIC = MAX(DBLE(NU*(4.0/5.0)*CACHE
     &                         - BUDGET - BUDGET_WS),0.0D0) /
     &                        (N*(COSTS_PREC(PREC)+COSTS_PREC(1)))
C                     IF (KMAGIC.LT.KREM.AND.KMAGIC.GT.8)
C     &                 KMAGIC = KMAGIC - MOD(KMAGIC,8)
                     USED_K = MIN(KMAGIC, KREM)
                     WS_FOOTPRINT = N * USED_K * COSTS_PREC(1)
                     PREC_FOOTPRINT = N * USED_K * COSTS_PREC(PREC)
                     BUDGET = BUDGET + PREC_FOOTPRINT
                     BUDGET_WS = BUDGET_WS + WS_FOOTPRINT
                     CURRENT_CONSUMPTION = CURRENT_CONSUMPTION +
     &                 DBLE(PREC_FOOTPRINT+WS_FOOTPRINT) / 
     &                 ((80.0/100.0)*DBLE(CACHE))
                     OVERALL_CONSUMPTION = OVERALL_CONSUMPTION +
     &                 DBLE(PREC_FOOTPRINT+WS_FOOTPRINT) / 
     &                 ((80.0/100.0)*DBLE(CACHE))
                     KREM = KREM - USED_K
                     CURRENT_K = CURRENT_K + USED_K
                     OVERALL_K = OVERALL_K + USED_K
                     IF (KREM.GE.1.AND.K455.EQ.4
     &                   .AND.CURRENT_CONSUMPTION.GT.0.3.OR.
     &                   KREM.GE.1.OR.PREC.EQ.NPREC_STOC.OR.
     &                   OVERALL_K.EQ.LRB%K) THEN
                        WSMAX = MAX(WSMAX, CURRENT_K * N)
                        BUDGET = 0
                        BUDGET_WS = 0
                        CURRENT_K = 0
                        CURRENT_CONSUMPTION = 0.0D0
                     ENDIF
                     IF (KREM.LT.1) EXIT
                   END DO
                 END DO 
               END IF
           END SELECT
           IF (K448.GE.1) THEN
             DO PREC=2,NPREC_STOC
               !TODO using PREC_CALC compute available cache_size
               !PREC_CALC = -1
               !DO PREC_=NPREC_CALC,1,-1
               !  IF (CALC2STOC(PREC_).LE.PREC) THEN
               !    PREC_CALC = CALC2STOC(PREC_)
               !    EXIT
               !  END IF
               !END DO 
               WSMAX = MAX(WSMAX, N*LRB%K_STOC(PREC))
               WSMAX = MAX(WSMAX, M*LRB%K_STOC(PREC))
             END DO
           END IF
         ELSE ! NOT ISLR
           CALL FRB_GET_PREC(LRB, PREC, K450.GT.0.OR.K148.EQ.2)
           if (K448.GE.1) THEN
             do prec_calc=NPREC_CALC,1,-1
               if (calc2stoc(prec_calc).le.prec) then
                 exit 
               end if
             end do
             prec_calc = min(nprec_calc, prec_calc)
             pcalc_stoc = calc2stoc(prec_calc)
             IF (PREC.NE.pcalc_stoc) THEN
               WSMAX = MAX(WSMAX,
     &            (M*N*COSTS_PREC(pcalc_stoc)-1)/
     &            COSTS_PREC(1)+1)
             END IF
           ELSE
             IF (PREC.GT.1) THEN
               ! THE N IS MAGIC ACTUALLY
               KMAGIC=(80.0/100.0)*CACHE/
     &           ((COSTS_PREC(1)+COSTS_PREC(PREC))*M)
               KMAGIC=MIN(KMAGIC,N)
               WSMAX = MAX(WSMAX, M * KMAGIC)
             END IF
           END IF
         END IF
      END SUBROUTINE COMPUTE_WORKSPACE_SIZE
      END MODULE SMUMPS_SOL_LR
********************************************************************
C     ------------------------------------
C     Implicit interface needed to keep
C     contiguity information on DIAG array
C     when passing an assumed shape array
C     known to be contiguous.
C     ------------------------------------
      SUBROUTINE SMUMPS_SOLVE_BWD_LR_TRSOLVE (
     &           DIAG, LDIAG, NPIV, NELIM, LIELL,
     &           NRHS_B, W, LWC, 
     &           RHSINTR, LRHSINTR, NRHS,
     &           PPIVINRHSINTR, JBDEB, 
     &           MTYPE, KEEP)
C
C Purpose:
C      Compute triangular solution during bwd step
C        IF ( MTYPE .EQ. 1) THEN
C         U X = Y
C         ELSE
C          L^T X = Y
C          (leading dimension of L is different
C           whether KEEP(50)=0 or not => LDAJ is adapted
C         ENDIF
C  
       INTEGER, INTENT(IN) :: MTYPE, LIELL, NPIV, NELIM, KEEP(500)
C      NRHS_B : size of current block of RHS being processed
       INTEGER, INTENT(IN) :: NRHS_B, LDIAG
       INTEGER, INTENT(IN) :: PPIVINRHSINTR, JBDEB, LRHSINTR, NRHS
       INTEGER(8), INTENT(IN) ::  LWC
       REAL, INTENT(IN) :: DIAG(LDIAG)
       REAL, INTENT(INOUT) :: W(LWC)
       REAL RHSINTR(LRHSINTR,NRHS)
C Local variables
C
       INTEGER :: LDAJ
C     
C     Constants
C     =========
      REAL ONE
      PARAMETER (ONE = 1.0E0)
        IF ( MTYPE .eq. 1 ) THEN
!        U x = y
         LDAJ =  NPIV + NELIM
C        write(6,*) " TRSOLVE:NPIV,NRHS_B,LDAJ=", NPIV,NRHS_B,LDAJ
         CALL strsm('L','L','T','N', NPIV, NRHS_B, ONE, DIAG(1),
     &              LDAJ, RHSINTR(PPIVINRHSINTR,JBDEB), 
     &              LRHSINTR)
        ELSE
         IF ( KEEP(50) .EQ. 0 ) THEN
C          L^t x = y
           LDAJ=NPIV+NELIM
         ELSE
C          L^t x = y
           LDAJ=NPIV
         ENDIF
         CALL strsm('L','U','N','U', NPIV, NRHS_B, ONE, DIAG(1),
     &            LDAJ, RHSINTR(PPIVINRHSINTR,JBDEB), LRHSINTR)
c        DO K=1, NRHS_B
C           TO BE VALIDATED
c           CALL strsv('U','N','U', NPIV, DIAG(1), LDAJ,
c    *           W(PPIV)+(K-1)*LIELL), 1)
c        ENDDO
        END IF ! MTYPE
      RETURN
      END SUBROUTINE SMUMPS_SOLVE_BWD_LR_TRSOLVE
