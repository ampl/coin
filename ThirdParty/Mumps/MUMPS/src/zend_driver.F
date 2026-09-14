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
      SUBROUTINE ZMUMPS_END_DRIVER( id, idintr )
      USE ZMUMPS_STRUC_DEF, ONLY: ZMUMPS_STRUC
      USE ZMUMPS_INTR_TYPES, ONLY: ZMUMPS_INTR_STRUC
      IMPLICIT NONE
C
C     Purpose:
C     =======
C
C     Terminate a MUMPS instance. Free all internal data structure and
C     suppress OOC files on disk, if any.
C
C     Argument:
C     ========
C
      TYPE( ZMUMPS_STRUC )      :: id
      TYPE( ZMUMPS_INTR_STRUC ) :: idintr
C
C     Local declarations
C     ==================
      INTEGER IERR
      INTEGER, PARAMETER :: MASTER = 0
C
C     Executable statements
C     =====================
C
C     First, free all MUMPS internal data except communicators created
C     during a call to MUMPS wit JOB=-1
      CALL ZMUMPS_FREE_DATA_ANAFACSOL( id, idintr )
C
C     Allocated during JOB=-1:
      IF ( id%MYID .NE. MASTER .OR. id%KEEP(46) .ne. 0 ) THEN
C       Note that on some very old platforms, COMM_NODES would have been
C       freed inside BLACS_GRIDEXIT, which may cause problems
C       in the call to MPI_COMM_FREE.
        CALL MPI_COMM_FREE( id%COMM_NODES, IERR )
C       Free communicator related to load messages.
        CALL MPI_COMM_FREE( id%COMM_LOAD, IERR )
      END IF
      CALL MUMPS_DESTROY_ARCH_NODE_COMM( id%KEEP(411), id%KEEP(410),
     &                                   id%KEEP(413) )
C     Nullifying id%SCHUR_CINTERFACE here is not necessary,
C     it is freed systematically each time we exit ZMUMPS_DRIVER
C     and reset each time we enter MUMPS through its C interface.
      NULLIFY(id%SCHUR_CINTERFACE)
C
      RETURN
      END SUBROUTINE ZMUMPS_END_DRIVER
C
      SUBROUTINE ZMUMPS_END_ROOT(roota)
      USE ZMUMPS_INTR_TYPES, ONLY: ZMUMPS_ROOT_STRUC
      IMPLICIT NONE
      TYPE(ZMUMPS_ROOT_STRUC) :: roota
      IF (associated(roota%RHS_CNTR_MASTER_ROOT)) THEN
        DEALLOCATE(roota%RHS_CNTR_MASTER_ROOT)
        NULLIFY(roota%RHS_CNTR_MASTER_ROOT)
      ENDIF
      IF (associated(roota%RHS_ROOT))THEN
        DEALLOCATE(roota%RHS_ROOT)
        NULLIFY(roota%RHS_ROOT)
      ENDIF
      CALL ZMUMPS_RR_FREE_POINTERS(roota)
      RETURN
      END SUBROUTINE ZMUMPS_END_ROOT
C
      SUBROUTINE ZMUMPS_FREE_DATA_ANAFACSOL(id, idintr)
      USE ZMUMPS_STRUC_DEF
      USE ZMUMPS_INTR_TYPES, ONLY: ZMUMPS_INTR_STRUC
      IMPLICIT NONE
C
C   Purpose:
C   =======
C   Free all MUMPS internal data, except communicators built during
C   a JOB=-1 call. Called by ZMUMPS_END_DRIVER and ZMUMPS_ANA_DRIVER.
C   Calls ZMUMPS_FREE_DATA_FACTO, which frees most of the data allocated
C   during factorization and solve, except:
C   - scaling arrays, because they are sometimes allocated at analysis
C   - STEP2NODE, which can be reused when analysis does not change
C   Therefore, scaling arrays and STEP2NODE are freed here.
C
C   Arguments
C   =========
      TYPE( ZMUMPS_STRUC )      :: id
      TYPE( ZMUMPS_INTR_STRUC ) :: idintr
C  Local declarations
C  ==================
      LOGICAL I_AM_SLAVE
      INTEGER, PARAMETER :: MASTER = 0
C
C     Executable statements
C     ---------------------
C
      I_AM_SLAVE = ( id%MYID .ne. MASTER .OR. id%KEEP(46) .NE. 0 )
C
C     First, free data from factoriation and solve:
      CALL ZMUMPS_FREE_DATA_FACTO(id,idintr)
C      -------------------------------------
C      Right-hand-side and solutions are
C      always user data, we do not free them
C      -------------------------------------
      IF (associated(id%MEM_DIST))  THEN
         DEALLOCATE(id%MEM_DIST)
         NULLIFY(id%MEM_DIST)
      ENDIF
C     ---------------------------------
C     Allocated by ZMUMPS, Used by user.
C     ZMUMPS deallocates. User should
C     use them before ZMUMPS_END_DRIVER or
C     copy.
C     ---------------------------------
      IF (associated(id%MAPPING)) THEN
        DEALLOCATE(id%MAPPING)
        NULLIFY(id%MAPPING)
      END IF
C     -------------------------------------
C     Always deallocate scaling arrays
C     if they are associated, except
C     when provided by the user (on master)
C     -------------------------------------
      IF ( id%KEEP(52) .NE. -1 .or. id%MYID .ne. MASTER ) THEN
        IF (associated(id%COLSCA)) THEN
          DEALLOCATE(id%COLSCA)
          NULLIFY(id%COLSCA)
        ENDIF
        IF (associated(id%ROWSCA)) THEN
          DEALLOCATE(id%ROWSCA)
          NULLIFY(id%ROWSCA)
        ENDIF
      END IF
C     Begin PRUN_NODES
C     Info for pruning tree 
      IF (associated(id%Step2node))      THEN
        DEALLOCATE(id%Step2node)
        NULLIFY(id%Step2node)
      ENDIF
C     END PRUN_NODES
c     --------------------- 
C     Allocated during analysis:
      IF (associated(id%NE_STEPS))  THEN
        DEALLOCATE(id%NE_STEPS)
        NULLIFY(id%NE_STEPS)
      ENDIF
C     Allocated during analysis:
      IF (associated(id%ND_STEPS))  THEN
        DEALLOCATE(id%ND_STEPS)
        NULLIFY(id%ND_STEPS)
      ENDIF
C     Allocated during analysis:
      IF (associated(id%FRERE_STEPS))  THEN
        DEALLOCATE(id%FRERE_STEPS)
        NULLIFY(id%FRERE_STEPS)
      ENDIF
C     Allocated during analysis:
      IF (associated(id%DAD_STEPS))  THEN
        DEALLOCATE(id%DAD_STEPS)
        NULLIFY(id%DAD_STEPS)
      ENDIF
C     Allocated during analysis:
      IF (associated(id%SYM_PERM))  THEN
          DEALLOCATE(id%SYM_PERM)
          NULLIFY(id%SYM_PERM)
      ENDIF
      IF (associated(id%STEP))  THEN
          DEALLOCATE(id%STEP)
          NULLIFY(id%STEP)
      ENDIF
      IF (associated(id%LRGROUPS))  THEN
          DEALLOCATE(id%LRGROUPS)
          NULLIFY(id%LRGROUPS)
      ENDIF
      IF (associated(id%FILS))      THEN
        DEALLOCATE(id%FILS)
        NULLIFY(id%FILS)
      ENDIF
      IF (associated(id%UNS_PERM))  THEN
        DEALLOCATE(id%UNS_PERM)
        NULLIFY(id%UNS_PERM)
      ENDIF
CC     Allocated during analysis:
      IF (associated(id%FRTPTR))    THEN
        DEALLOCATE(id%FRTPTR)
        NULLIFY(id%FRTPTR)
      ENDIF
C     Allocated during analysis:
      IF (associated(id%FRTELT))    THEN
        DEALLOCATE(id%FRTELT)
        NULLIFY(id%FRTELT)
      ENDIF
C     Allocated during analysis:
      IF (associated(id%NA))        THEN
        DEALLOCATE(id%NA)
        NULLIFY(id%NA)
      ENDIF
C     Allocated during analysis:
      IF (associated(id%PROCNODE_STEPS)) THEN
        DEALLOCATE(id%PROCNODE_STEPS)
        NULLIFY(id%PROCNODE_STEPS)
      ENDIF
      IF (associated(id%ELTPROC))     THEN
        DEALLOCATE(id%ELTPROC)
        NULLIFY(id%ELTPROC)
      ENDIF
C     Allocated at analysis:
      IF (associated(id%CANDIDATES)) THEN
        DEALLOCATE(id%CANDIDATES)
        NULLIFY(id%CANDIDATES)
      ENDIF
C     Allocated at analysis:
      IF (associated(id%I_AM_CAND)) THEN
        DEALLOCATE(id%I_AM_CAND)
        NULLIFY(id%I_AM_CAND)
      ENDIF
C     Allocated at analysis:
      IF (associated(id%ISTEP_TO_INIV2)) THEN
        DEALLOCATE(id%ISTEP_TO_INIV2)
        NULLIFY(id%ISTEP_TO_INIV2)
      ENDIF
C     Node partitionning (only allocated on slaves)
      IF (I_AM_SLAVE) THEN
C      Allocated at analysis:
       IF (associated(id%TAB_POS_IN_PERE)) THEN
        DEALLOCATE(id%TAB_POS_IN_PERE)
        NULLIFY(id%TAB_POS_IN_PERE)
       ENDIF
      ENDIF
      IF (I_AM_SLAVE) THEN
C      Allocated at analysis:
       IF (associated(id%FUTURE_NIV2)) THEN
        DEALLOCATE(id%FUTURE_NIV2)
        NULLIFY(id%FUTURE_NIV2)
       ENDIF
      ENDIF
C     Allocated at analysis:
      IF(associated(id%DEPTH_FIRST))THEN
        DEALLOCATE(id%DEPTH_FIRST)
        NULLIFY(id%DEPTH_FIRST)
      ENDIF
C     Allocated at analysis:
      IF(associated(id%DEPTH_FIRST_SEQ))THEN
        DEALLOCATE(id%DEPTH_FIRST_SEQ)
        NULLIFY(id%DEPTH_FIRST_SEQ)
      ENDIF
C     Allocated at analysis:
      IF(associated(id%SBTR_ID))THEN
        DEALLOCATE(id%SBTR_ID)
        NULLIFY(id%SBTR_ID)
      ENDIF
C     Allocated and initialized at analysis:
      IF(associated(id%SCHED_DEP))THEN
        DEALLOCATE(id%SCHED_DEP)
        NULLIFY(id%SCHED_DEP)
      ENDIF
C     Allocated and initialized at analysis:
      IF(associated(id%SCHED_SBTR))THEN
        DEALLOCATE(id%SCHED_SBTR)
        NULLIFY(id%SCHED_SBTR)
      ENDIF
C     Allocated and initialized at analysis:
      IF(associated(id%SCHED_GRP))THEN
        DEALLOCATE(id%SCHED_GRP)
        NULLIFY(id%SCHED_GRP)
      ENDIF
C     Allocated and initialized at analysis:
      IF(associated(id%CROIX_MANU))THEN
        DEALLOCATE(id%CROIX_MANU)
        NULLIFY(id%CROIX_MANU)
      ENDIF
C     Allocated during analysis:
      IF (associated(id%PTRAR))     THEN
        DEALLOCATE(id%PTRAR)
        NULLIFY(id%PTRAR)
      ENDIF
C     Allocated at analysis:
      IF (associated(id%MEM_SUBTREE)) THEN
        DEALLOCATE(id%MEM_SUBTREE)
        NULLIFY(id%MEM_SUBTREE)
      ENDIF
C     Allocated at analysis:
      IF (associated(id%MY_ROOT_SBTR)) THEN
        DEALLOCATE(id%MY_ROOT_SBTR)
        NULLIFY(id%MY_ROOT_SBTR)
      ENDIF
C     Allocated at analysis:
      IF (associated(id%MY_FIRST_LEAF)) THEN
        DEALLOCATE(id%MY_FIRST_LEAF)
        NULLIFY(id%MY_FIRST_LEAF)
      ENDIF
C     Allocated at analysis:
      IF (associated(id%MY_NB_LEAF)) THEN
        DEALLOCATE(id%MY_NB_LEAF)
        NULLIFY(id%MY_NB_LEAF)
      ENDIF
C     Allocated at analysis:
      IF (associated(id%COST_TRAV)) THEN
        DEALLOCATE(id%COST_TRAV)
        NULLIFY(id%COST_TRAV)
      ENDIF
C     Allocated at analysis:
      IF (associated(id%CB_SON_SIZE)) THEN
        DEALLOCATE(id%CB_SON_SIZE)
        NULLIFY(id%CB_SON_SIZE)
      ENDIF
C     Allocated at analysis:
      IF (associated(id%SUP_PROC)) THEN
         DEALLOCATE(id%SUP_PROC)
         NULLIFY(id%SUP_PROC)
      ENDIF
!     IF(id%KEEP(486).NE.0) THEN
C     Allocated at analysis:
        IF (associated(id%LRGROUPS)) THEN
           DEALLOCATE(id%LRGROUPS)
           NULLIFY(id%LRGROUPS)
        ENDIF
!     ENDIF
C
C     free data concerned when redoing cheap analysis
      CALL ZMUMPS_FREE_DATA_REDO_ANA( id )
C
C     gridinit performed at analysis:
#if ! defined(NOSCALAPACK)
      IF (idintr%root%gridinit_done) THEN
        IF ( id%KEEP(38).NE.0 .and. idintr%root%yes ) THEN
          CALL blacs_gridexit( idintr%root%CNTXT_BLACS )
          idintr%root%gridinit_done = .FALSE.
        END IF
      END IF
#endif
      RETURN
      END SUBROUTINE ZMUMPS_FREE_DATA_ANAFACSOL
      SUBROUTINE ZMUMPS_FREE_DATA_REDO_ANA ( id )
      USE ZMUMPS_STRUC_DEF
      IMPLICIT NONE
C
C   Purpose:
C   =======
C   Free all MUMPS internal data concerned 
C   when redoing a cheap analysis :
C    - data related to MPI2KOMP allocated during analysis
C    - data related to L0OMP allocated during analysis
C    - data related to building arrowheads because 
C      of EARLYT3ROOTINS that might change when of
C      L0-thread (KEEP(400)
C   Arguments
C   =========
      TYPE( ZMUMPS_STRUC ) :: id
C
C     Executable statements
C     ---------------------
CCN#if defined(MPI_TO_K_OMP)
C     Allocated at analysis:
      IF (associated(id%MTKO_PROCS_MAP)) THEN
        DEALLOCATE(id%MTKO_PROCS_MAP)
        NULLIFY(id%MTKO_PROCS_MAP)
      ENDIF
C     Allocated at analysis:
      IF (associated(id%IPOOL_B_L0_OMP)) THEN
        DEALLOCATE(id%IPOOL_B_L0_OMP)
        NULLIFY(id%IPOOL_B_L0_OMP)
      END IF
      IF (associated(id%IPOOL_A_L0_OMP)) THEN
        DEALLOCATE(id%IPOOL_A_L0_OMP)
        NULLIFY(id%IPOOL_A_L0_OMP)
      END IF
      IF (associated(id%PHYS_L0_OMP)) THEN
        DEALLOCATE(id%PHYS_L0_OMP)
        NULLIFY(id%PHYS_L0_OMP)
      END IF
      IF (associated(id%VIRT_L0_OMP)) THEN
        DEALLOCATE(id%VIRT_L0_OMP)
        NULLIFY(id%VIRT_L0_OMP)
      END IF
      IF (associated(id%VIRT_L0_OMP_MAPPING)) THEN
        DEALLOCATE(id%VIRT_L0_OMP_MAPPING)
        NULLIFY(id%VIRT_L0_OMP_MAPPING)
      END IF
      IF (associated(id%PERM_L0_OMP)) THEN
        DEALLOCATE(id%PERM_L0_OMP)
        NULLIFY(id%PERM_L0_OMP)
      END IF
      IF (associated(id%PTR_LEAFS_L0_OMP)) THEN
        DEALLOCATE(id%PTR_LEAFS_L0_OMP)
        NULLIFY(id%PTR_LEAFS_L0_OMP)
      END IF
C     Allocated at analysis:
      IF (associated(id%I4_L0_OMP)) THEN
        DEALLOCATE(id%I4_L0_OMP)
        NULLIFY(id%I4_L0_OMP)
      END IF
C     Allocated at analysis:
      IF (associated(id%I8_L0_OMP)) THEN
        DEALLOCATE(id%I8_L0_OMP)
        NULLIFY(id%I8_L0_OMP)
      END IF
C     =================================================
C     BEGIN Pointers to original matrix
C     allocated during analysis
C     in format ready for assembly during factorization
C     (arrowheads if assembled format)
C     Allocated during analysis:
C     id%PTRAR is allocated in ana_driver and
C     should not be deallocated here (it does not 
C     change in sze)
      IF (associated(id%PTR8ARR))   THEN
        DEALLOCATE(id%PTR8ARR)
        NULLIFY(id%PTR8ARR)
      ENDIF
C     Allocated during analysis:
      IF (associated(id%NINCOLARR))   THEN
        DEALLOCATE(id%NINCOLARR)
        NULLIFY(id%NINCOLARR)
      ENDIF
C     Allocated during analysis:
      IF (associated(id%NINROWARR))   THEN
        DEALLOCATE(id%NINROWARR)
        NULLIFY(id%NINROWARR)
      ENDIF
C     Allocated during analysis:
      IF (associated(id%PTRDEBARR))   THEN
        DEALLOCATE(id%PTRDEBARR)
        NULLIFY(id%PTRDEBARR)
      ENDIF
C     =================================================
      RETURN
      END SUBROUTINE ZMUMPS_FREE_DATA_REDO_ANA
      SUBROUTINE ZMUMPS_FREE_ID_DATA_MODULES(id_FDM_F_ENCODING,
     &  id_BLRARRAY_ENCODING, KEEP8, K34, K148)
      USE MUMPS_FRONT_DATA_MGT_M, only : MUMPS_FDM_STRUC_TO_MOD,
     &                                   MUMPS_FDM_END
      USE ZMUMPS_LR_DATA_M, only : ZMUMPS_BLR_STRUC_TO_MOD,
     &                             ZMUMPS_BLR_END_MODULE
#if defined(SINGLEINDOUBLE)
      USE CMUMPS_LR_DATA_M, only: CMUMPS_BLR_STRUC_TO_MOD,
     &                             CMUMPS_BLR_END_MODULE
#endif
      IMPLICIT NONE
C
C     Purpose:
C     =======
C
C     Free data from modules kept from one phase to the other
C     and referenced through the main MUMPS structure, id.
C
C     Both id%FDM_F_ENCODING and id%BLRARRAY_ENCODING
C     are concerned.
C
C
C
C     Arguments:
C     =========
C
#     if defined(MUMPS_NOF2003)
      CHARACTER, DIMENSION(:), POINTER :: id_BLRARRAY_ENCODING
      CHARACTER, DIMENSION(:), POINTER :: id_FDM_F_ENCODING
#     else
      CHARACTER, DIMENSION(:), POINTER, intent(inout) ::
     &                                            id_BLRARRAY_ENCODING
      CHARACTER, DIMENSION(:), POINTER, intent(inout) ::
     &                                            id_FDM_F_ENCODING
#     endif
      INTEGER(8), intent(inout) :: KEEP8(150)
      INTEGER, intent(in) :: K34, K148
C
      IF (associated(id_FDM_F_ENCODING)) THEN
C       Allow access to FDM_F data for BLR_END_MODULE
        CALL MUMPS_FDM_STRUC_TO_MOD('F', id_FDM_F_ENCODING)
        IF (associated(id_BLRARRAY_ENCODING)) THEN
C         Pass id_BLRARRAY_ENCODING control to module
C         and terminate BLR module of current instance
#if defined(SINGLEINDOUBLE)
          IF (K148 .NE. 1) THEN
C         to cover both cases KEEP(148)=0, 2
#endif
          CALL ZMUMPS_BLR_STRUC_TO_MOD(id_BLRARRAY_ENCODING)
          CALL ZMUMPS_BLR_END_MODULE(0, KEEP8, K34,
     &                               LRSOLVE_ACT_OPT=.TRUE.)
#if defined(SINGLEINDOUBLE)
          ELSE
            CALL CMUMPS_BLR_STRUC_TO_MOD(id_BLRARRAY_ENCODING)
            CALL CMUMPS_BLR_END_MODULE(0, KEEP8, K34, 
     &                               LRSOLVE_ACT_OPT=.TRUE.)
          ENDIF
#endif
        ENDIF
C       ---------------------------------------
C       FDM data structures are still allocated
C       in the module and should be freed
C       ---------------------------------------
        CALL MUMPS_FDM_END('F')
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_FREE_ID_DATA_MODULES
C
C     -----------------------------------------------------------------
C
      SUBROUTINE ZMUMPS_FREE_DATA_FACTO(id,idintr)
C
C     Purpose:
C     -------
C
C     ZMUMPS_FREE_DATA_FACTO frees data that was allocated during
C     factorization and that can be useful for the solve. Afterwards,
C     data from analysis is kept so that a new factorization phase
C     is possible.
C     
C     Module depencies
C     ----------------
      USE ZMUMPS_STRUC_DEF, ONLY: ZMUMPS_STRUC
      USE ZMUMPS_INTR_TYPES, ONLY: ZMUMPS_INTR_STRUC
      USE ZMUMPS_FACSOL_L0OMP_M, ONLY : ZMUMPS_FREE_L0_OMP_FACTORS
#if defined(SINGLEINDOUBLE)
      USE CMUMPS_FACSOL_L0OMP_M, ONLY : CMUMPS_FREE_L0_OMP_FACTORS
#endif
      USE ZMUMPS_DYNAMIC_MEMORY_M, ONLY : ZMUMPS_DM_FREE_S_WK
#if defined(SINGLEINDOUBLE)
      USE CMUMPS_DYNAMIC_MEMORY_M, ONLY : CMUMPS_DM_FREE_S_WK
#endif
      USE MUMPS_BUF_COMMON, ONLY : 
     &    MUMPS_BUF_DEALL_CB,
     &    MUMPS_BUF_DEALL_SMALL_BUF
      IMPLICIT NONE
C
C     Argument:
C     --------
C
C     id is the main MUMPS structure, giving with idintr access
C     to all internal objects allocated by the package.
C
      TYPE( ZMUMPS_STRUC) :: id
      TYPE( ZMUMPS_INTR_STRUC ) :: idintr
C
C     Local declarations
C     ------------------
      INTEGER :: IERR
      LOGICAL :: I_AM_SLAVE
      INTEGER, PARAMETER :: MASTER = 0
C
C     Interface blocks
C     ----------------
      INTERFACE
C     (explicit needed because of pointer arguments)
      SUBROUTINE ZMUMPS_FREE_ID_DATA_MODULES(id_FDM_F_ENCODING,
     &  id_BLRARRAY_ENCODING, KEEP8, K34, K148)
#     if defined(MUMPS_NOF2003)
      CHARACTER, DIMENSION(:), POINTER :: id_BLRARRAY_ENCODING
      CHARACTER, DIMENSION(:), POINTER :: id_FDM_F_ENCODING
#     else
      CHARACTER, DIMENSION(:), POINTER, intent(inout) ::
     &                                            id_BLRARRAY_ENCODING
      CHARACTER, DIMENSION(:), POINTER, intent(inout) ::
     &                                            id_FDM_F_ENCODING
#     endif
      INTEGER(8), intent(inout) :: KEEP8(150)
      INTEGER, intent(in) :: K34, K148
      END SUBROUTINE ZMUMPS_FREE_ID_DATA_MODULES
      END INTERFACE
C
      I_AM_SLAVE = ( id%MYID .ne. MASTER .OR. id%KEEP(46) .NE. 0 )
C
C     Free OOC-related data
C     ---------------------
C     (this includes suppression of OOC files)
      IF (id%KEEP(201).GT.0 .AND. I_AM_SLAVE) THEN
        CALL ZMUMPS_CLEAN_OOC_DATA(id,IERR)
        IF (IERR < 0) THEN
          id%INFO(1) = -90
          id%INFO(2) = 0
        ENDIF
      END IF
      CALL MUMPS_PROPINFO(id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID)
      IF (associated(id%PTLUST_S)) THEN
        DEALLOCATE(id%PTLUST_S)
        NULLIFY(id%PTLUST_S)
      END IF
      IF (associated(id%PTRFAC)) THEN
        DEALLOCATE(id%PTRFAC)
        NULLIFY(id%PTRFAC)
      END IF
      IF (associated(id%IS)) THEN
        DEALLOCATE(id%IS)
        NULLIFY(id%IS)
      ENDIF
      IF (associated(id%PIVNUL_LIST))  THEN
        DEALLOCATE(id%PIVNUL_LIST)
        NULLIFY(id%PIVNUL_LIST)
      ENDIF
      IF (associated(id%ROWSCA_loc)) THEN
        DEALLOCATE(id%ROWSCA_loc)
        NULLIFY(id%ROWSCA_loc)
      ENDIF
      IF (id%KEEP(50) .EQ. 0) THEN
        IF (associated(id%COLSCA_loc)) THEN
          DEALLOCATE(id%COLSCA_loc)
        ENDIF
      ENDIF
      NULLIFY(id%COLSCA_loc)
C     IPIV is used both for ScaLAPACK and RR
C     Keep it outside ZMUMPS_RR_FREE_POINTERS
      IF (associated(idintr%root%IPIV))    THEN
        DEALLOCATE(idintr%root%IPIV)
        NULLIFY(idintr%root%IPIV)
      ENDIF
      CALL ZMUMPS_END_ROOT(idintr%roota)
      IF (associated(id%SINGULAR_VALUES)) THEN
        DEALLOCATE(id%SINGULAR_VALUES)
        NULLIFY(id%SINGULAR_VALUES)
      ENDIF
C     Free module data from factorization:
      CALL ZMUMPS_FREE_ID_DATA_MODULES(id%FDM_F_ENCODING, ! done
     &  id%BLRARRAY_ENCODING, id%KEEP8(1), id%KEEP(34), id%KEEP(148))
C     ---------------------------
C     Deallocate main workarray S 
C     ---------------------------
      IF (id%KEEP8(24).EQ.0_8) THEN
C       -- deallocate only when not provided/allocated by the user
        IF (associated(id%S)) THEN
          CALL ZMUMPS_DM_FREE_S_WK(id%S, id%KEEP(430))
        ENDIF
#if defined(SINGLEINDOUBLE)
        IF (associated(id%LPS)) THEN
          CALL CMUMPS_DM_FREE_S_WK(id%LPS, id%KEEP(430))
          NULLIFY(id%LPS)
        ENDIF
#endif
C       Reset KEEP(430)=0 since S is free
C       KEEP(430) will be redefined during facto
        id%KEEP(430) = 0
C       Update allocated size of S:
        id%KEEP8(23)=0_8
      ENDIF
      NULLIFY(id%S)
      IF (I_AM_SLAVE) THEN
C       ------------------------
C       Deallocate buffer for
C       contrib-blocks (facto/
C       solve). Note that this
C       will cancel all possible
C       pending requests.
C       ------------------------
        CALL MUMPS_BUF_DEALL_CB( IERR )
C       Deallocate buffer for integers (facto/solve)
        CALL MUMPS_BUF_DEALL_SMALL_BUF( IERR )
      END IF
C
      IF (associated(id%L0_OMP_MAPPING)) THEN
        DEALLOCATE(id%L0_OMP_MAPPING)
        NULLIFY(id%L0_OMP_MAPPING)
      END IF
      IF (associated(idintr%L0_OMP_FACTORS)) THEN
        CALL ZMUMPS_FREE_L0_OMP_FACTORS(idintr%L0_OMP_FACTORS)
      END IF
#if defined(SINGLEINDOUBLE)
        CALL CMUMPS_FREE_L0_OMP_FACTORS(idintr%LPL0_OMP_FACTORS)
#endif
C
C     Data allocated during solve
C     ---------------------------
C
C     (or for some of it, factorization -- forward during factorization)
      IF (associated(id%RHSINTR)) THEN
        DEALLOCATE(id%RHSINTR)
        NULLIFY(id%RHSINTR)
        id%KEEP8(25)=0_8
      ENDIF
      IF (associated(id%GLOB2LOC_RHS)) THEN
        DEALLOCATE(id%GLOB2LOC_RHS)
        NULLIFY(id%GLOB2LOC_RHS)
      ENDIF
      IF (id%GLOB2LOC_SOL_ALLOC) THEN
        DEALLOCATE(id%GLOB2LOC_SOL)
        NULLIFY(id%GLOB2LOC_SOL)
        id%GLOB2LOC_SOL_ALLOC = .FALSE.
      ENDIF
C     Allocated during solve:
C     (even in case of fwd in facto)
      IF (associated(id%IPTR_WORKING)) THEN
        DEALLOCATE(id%IPTR_WORKING)
        NULLIFY(id%IPTR_WORKING)
      END IF
      IF (associated(id%WORKING)) THEN 
        DEALLOCATE(id%WORKING)
        NULLIFY(id%WORKING)
      END IF
      RETURN
      END SUBROUTINE ZMUMPS_FREE_DATA_FACTO
      SUBROUTINE ZMUMPS_FREE_DATA_RHSINTR(id)
C
C     Purpose:
C     -------
C       Free RHSINTR related data that might 
C       have been generated after a forward only step (ICNTL(26)=1)
C     Module depencies
C     ----------------
      USE ZMUMPS_STRUC_DEF
      IMPLICIT NONE
C
C     Argument:
C     --------
C
C     id is the main MUMPS structure, giving with idintr access
C     to all internal objects allocated by the package.
C
      TYPE( ZMUMPS_STRUC) :: id
C
      IF (associated(id%RHSINTR)) THEN
        DEALLOCATE(id%RHSINTR)
        NULLIFY(id%RHSINTR)
        id%KEEP8(25)=0_8
        id%LD_RHSINTR = 0
      ENDIF
      IF (associated(id%GLOB2LOC_RHS)) THEN
        DEALLOCATE(id%GLOB2LOC_RHS)
        NULLIFY(id%GLOB2LOC_RHS)
      ENDIF
      IF (id%GLOB2LOC_SOL_ALLOC) THEN
        DEALLOCATE(id%GLOB2LOC_SOL)
        NULLIFY(id%GLOB2LOC_SOL)
        id%GLOB2LOC_SOL_ALLOC = .FALSE.
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_FREE_DATA_RHSINTR
      SUBROUTINE ZMUMPS_CLEAN_OOC_DATA(id,IERR)
      USE ZMUMPS_STRUC_DEF, ONLY : ZMUMPS_STRUC
      IMPLICIT NONE
      TYPE(ZMUMPS_STRUC) :: id
      INTEGER IERR
      IERR=0
      CALL ZMUMPS_OOC_CLEAN_FILES(id,IERR)
      IF(associated(id%OOC_TOTAL_NB_NODES))THEN
         DEALLOCATE(id%OOC_TOTAL_NB_NODES)
         NULLIFY(id%OOC_TOTAL_NB_NODES)
      ENDIF
      IF(associated(id%OOC_INODE_SEQUENCE))THEN
         DEALLOCATE(id%OOC_INODE_SEQUENCE)
         NULLIFY(id%OOC_INODE_SEQUENCE)
      ENDIF
      IF(associated(id%OOC_SIZE_OF_BLOCK))THEN
         DEALLOCATE(id%OOC_SIZE_OF_BLOCK)
         NULLIFY(id%OOC_SIZE_OF_BLOCK)
      ENDIF
      IF(associated(id%OOC_VADDR))THEN
         DEALLOCATE(id%OOC_VADDR)
         NULLIFY(id%OOC_VADDR)
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_CLEAN_OOC_DATA
      SUBROUTINE ZMUMPS_OOC_CLEAN_FILES(id,IERR)
      USE ZMUMPS_STRUC_DEF
      USE MUMPS_OOC_COMMON, ONLY : ERR_STR_OOC,
     &                             DIM_ERR_STR_OOC,
     &                             FILENAMELENGTH
      IMPLICIT NONE
      EXTERNAL MUMPS_OOC_REMOVE_FILE_C
      TYPE(ZMUMPS_STRUC) :: id
      INTEGER IERR
      INTEGER I,J,I1,K
      CHARACTER(len=1):: TMP_NAME(FILENAMELENGTH)
C     Note that Fortran initializes IERR to 0.
C     The C layer modifies it in case of error.
      IERR=0
      K=1
C     WHEN SAVE/RESTORE IS ON, OOC FILES ASSOCIATED TO A SAVED INSTANCE
C     ARE NOT REMOVED
      IF(.NOT. id%ASSOCIATED_OOC_FILES) THEN
        IF(associated(id%OOC_FILE_NAMES).AND.
     &       associated(id%OOC_FILE_NAME_LENGTH))THEN
           DO I1=1,id%OOC_NB_FILE_TYPE
              DO I=1,id%OOC_NB_FILES(I1)
                 DO J=1,id%OOC_FILE_NAME_LENGTH(K)
                    TMP_NAME(J)=id%OOC_FILE_NAMES(K,J)
                 ENDDO
C       Note that termination character '0' is included
                 CALL MUMPS_OOC_REMOVE_FILE_C(IERR, TMP_NAME(1))
                 IF(IERR.LT.0)THEN
                    IF (id%ICNTL(1).GT.0 .AND. id%ICNTL(4).GE.1)THEN
                       WRITE(id%ICNTL(1),*) id%MYID,': ',
     &                      ERR_STR_OOC(1:DIM_ERR_STR_OOC)
                       RETURN
                    ENDIF
                 ENDIF
                 K=K+1
              ENDDO
           ENDDO
        ENDIF
      ENDIF
      IF(associated(id%OOC_FILE_NAMES))THEN
         DEALLOCATE(id%OOC_FILE_NAMES)
         NULLIFY(id%OOC_FILE_NAMES)
      ENDIF
      IF(associated(id%OOC_FILE_NAME_LENGTH))THEN
         DEALLOCATE(id%OOC_FILE_NAME_LENGTH)
         NULLIFY(id%OOC_FILE_NAME_LENGTH)
      ENDIF      
      IF(associated(id%OOC_NB_FILES))THEN
         DEALLOCATE(id%OOC_NB_FILES)
         NULLIFY(id%OOC_NB_FILES)
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_OOC_CLEAN_FILES
