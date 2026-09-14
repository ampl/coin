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
      SUBROUTINE ZMUMPS_FAC_DRIVER(id,idintr)
      USE MUMPS_BUF_COMMON, ONLY : MUMPS_BUF_INI_MYID, MUMPS_BUF_INIT,
     &    MUMPS_BUF_ALLOC_SMALL_BUF, MUMPS_BUF_DEALL_SMALL_BUF,
     &    MUMPS_BUF_DIST_IRECV_SIZE
      USE MUMPS_LOAD
      USE ZMUMPS_OOC, ONLY : ZMUMPS_OOC_INIT_FACTO,
     &                       ZMUMPS_OOC_END_FACTO
      USE ZMUMPS_STRUC_DEF
      USE ZMUMPS_INTR_TYPES, ONLY: ZMUMPS_INTR_STRUC
#if defined(SINGLEINDOUBLE)
      USE CMUMPS_INTR_TYPES, ONLY: CMUMPS_ROOT_STRUC
#endif
      USE ZMUMPS_FACSOL_L0OMP_M, ONLY: ZMUMPS_INIT_L0_OMP_FACTORS
#if defined(SINGLEINDOUBLE)
      USE CMUMPS_FACSOL_L0OMP_M, ONLY: CMUMPS_FREE_L0_OMP_FACTORS,
     &  CMUMPS_INIT_L0_OMP_FACTORS
#endif
      USE ZMUMPS_DYNAMIC_MEMORY_M, ONLY: ZMUMPS_DM_ALLOC_S_WK,
     &    ZMUMPS_DM_FREE_S_WK
#if defined(SINGLEINDOUBLE)
      USE CMUMPS_DYNAMIC_MEMORY_M, ONLY: CMUMPS_DM_ALLOC_S_WK,
     &    CMUMPS_DM_FREE_S_WK
#endif
      USE MUMPS_LR_STATS
      USE ZMUMPS_LR_DATA_M, only: ZMUMPS_BLR_INIT_MODULE, 
     &                            ZMUMPS_BLR_END_MODULE
     &                          , ZMUMPS_BLR_MOD_TO_STRUC
       USE ZMUMPS_FAC_COMPACT_FACTORS_M, ONLY:
     &     ZMUMPS_TRY_COMPACT_FACTORS
#if defined(SINGLEINDOUBLE)
      USE CMUMPS_LR_DATA_M, only: CMUMPS_BLR_INIT_MODULE, 
     &                            CMUMPS_BLR_END_MODULE
     &                          , CMUMPS_BLR_MOD_TO_STRUC
      USE CMUMPS_FAC_COMPACT_FACTORS_M, ONLY:
     &    CMUMPS_TRY_COMPACT_FACTORS
#endif
      USE ZMUMPS_LR_TYPE, only : NPREC_STOC,
     &              ZMUMPS_INIT_LRMIXED_PRECISIONS
#if defined(SINGLEINDOUBLE)
      USE CMUMPS_LR_TYPE, only : NPREC_STOC,
     &              CMUMPS_INIT_LRMIXED_PRECISIONS
#endif
      USE MUMPS_PIVNUL_MOD, only: PIVNUL_LIST_STRUCT_T
      USE MUMPS_FRONT_DATA_MGT_M
#if ! defined(NO_FDM_DESCBAND)
      USE MUMPS_FAC_DESCBAND_DATA_M
#endif
#if ! defined(NO_FDM_MAPROW)
      USE MUMPS_FAC_MAPROW_DATA_M
#endif
#if defined(MPI_TO_K_OMP)
      USE MUMPS_MPITOOMP_M
#endif
!$    USE OMP_LIB
C     Derived datatype to pass pointers with implicit interfaces
      USE ZMUMPS_FAC_S_IS_POINTERS_M, ONLY : ZMUMPS_S_IS_POINTERS_T
#if defined(SINGLEINDOUBLE)
      USE CMUMPS_FAC_S_IS_POINTERS_M, ONLY : CMUMPS_S_IS_POINTERS_T
#endif
      IMPLICIT NONE
C
C  Purpose
C  =======
C
C  Performs scaling, sorting in arrowhead, then
C  distributes the matrix, and perform
C  factorization.
C
C
      INTERFACE
      SUBROUTINE ZMUMPS_ANORMINF(id, ANORMINF, LSCAL, EFF_SIZE_SCHUR)
      USE ZMUMPS_STRUC_DEF
      TYPE (ZMUMPS_STRUC), TARGET :: id
      DOUBLE PRECISION, INTENT(OUT) :: ANORMINF
      LOGICAL, INTENT(IN) :: LSCAL
      INTEGER, INTENT(IN) :: EFF_SIZE_SCHUR
      END SUBROUTINE ZMUMPS_ANORMINF
#if defined(SINGLEINDOUBLE)
      SUBROUTINE ZMUMPS_LPid_TO_Cid(id, idintr, K34, IFLAG, IERROR)
      USE ZMUMPS_STRUC_DEF,  ONLY: ZMUMPS_STRUC 
      USE ZMUMPS_INTR_TYPES, ONLY: ZMUMPS_INTR_STRUC
      USE CMUMPS_STRUC_DEF, ONLY: CMUMPS_STRUC
      USE CMUMPS_LR_DATA_M, BLR_ARRAY_LP => BLR_ARRAY
      USE ZMUMPS_LR_DATA_M,  BLR_ARRAY_C  => BLR_ARRAY
      USE ZMUMPS_FACSOL_L0OMP_M, ONLY: ZMUMPS_INIT_L0_OMP_FACTORS
      USE CMUMPS_FACSOL_L0OMP_M, ONLY: CMUMPS_FREE_L0_OMP_FACTORS
      use iso_c_binding
      TYPE(ZMUMPS_STRUC), INTENT(INOUT)       :: id
      TYPE (ZMUMPS_INTR_STRUC), intent(inout) :: idintr
      INTEGER :: K34
      INTEGER, intent(inout) :: IFLAG, IERROR
      END SUBROUTINE ZMUMPS_LPid_TO_Cid
#endif
      END INTERFACE
C
C  Parameters
C  ==========
C
      TYPE (ZMUMPS_STRUC), TARGET :: id
      TYPE (ZMUMPS_INTR_STRUC)    :: idintr
C
C  MPI
C  ===
C
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR
      INTEGER, PARAMETER :: MASTER = 0
C
C  Local variables
C  ===============
C
      INCLUDE 'mumps_headers.h'
      INTEGER(8) :: NSEND8, NSEND_TOT8
      INTEGER(8) :: NLOCAL8, NLOCAL_TOT8
      INTEGER(4) :: I4
      INTEGER :: LDPTRAR, NELT_arg, NBRECORDS
      INTEGER :: ITMP, JTMP
      INTEGER :: KEEP464COPY, KEEP465COPY
      INTEGER(8) :: KEEP826_SAVE
      INTEGER(8) :: K67, K68, K70, K74, K75
      INTEGER(8) ITMP8
      INTEGER  MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      INTEGER MP, LP, MPG, allocok
      LOGICAL PROK, PROKG, LSCAL, LPOK, COMPUTE_ANORMINF
C     Reception buffer
      INTEGER    :: ZMUMPS_LBUFR, ZMUMPS_LBUFR_BYTES
      INTEGER(8) :: ZMUMPS_LBUFR_BYTES8 ! for intermediate computation
      INTEGER, ALLOCATABLE, DIMENSION(:) :: BUFR
C     Size of send buffers (in bytes)
      INTEGER    :: ZMUMPS_LBUF, ZMUMPS_LBUF_INT
      INTEGER(8) :: ZMUMPS_LBUF8 ! for intermediate computation
C
      INTEGER PTRIST, PTRWB, MAXELT_SIZE,
     &     ITLOC, IPOOL, LPOOL
      INTEGER IRANK, ID_ROOT
      INTEGER KKKK
      INTEGER(8) :: NZ_locMAX8
      INTEGER(8) MEMORY_MD_ARG
      INTEGER(8) MAXS_BASE8, MAXS_BASE_RELAXED8
      DOUBLE PRECISION CNTL4, AVG_FLOPS
      INTEGER MIN_PERLU, MAXIS_ESTIM
      INTEGER SUM_INFO22_THIS_NODE, MAX_SUM_INFO22_THIS_NODE
C
      TYPE (ZMUMPS_S_IS_POINTERS_T) :: S_IS_POINTERS
#if defined(USE_GPU)
#if defined(SINGLEINDOUBLE)
      REAL(kind(0.E0)) :: LPRINFO9, LPRINFO10
#endif
#endif
#if defined(SINGLEINDOUBLE)
C     Parameters to call CMUMPS_FAC_B
      TYPE (CMUMPS_S_IS_POINTERS_T) :: LPS_IS_POINTERS
      COMPLEX(kind=4), DIMENSION(:), ALLOCATABLE :: LPDBLARR
      COMPLEX(kind=4), DIMENSION(:), ALLOCATABLE :: LPRHS_MUMPS
      TYPE(CMUMPS_ROOT_STRUC) :: lproota 
      REAL(kind(0.E0)) :: LPCNTL1
      REAL(kind(0.E0)) :: LPDK8
#endif
      INTEGER    :: MAXIS
      INTEGER(8) :: MAXS
#if defined(SINGLEINDOUBLE)
      INTEGER(8) :: I8
#endif
      INTEGER    :: ICNTL49_LOC, TMP_INFOG_4, ICNTL36_LOC
C     For S argument to arrowhead routines:
      INTEGER(8) :: MAXS_ARG
      COMPLEX(kind=8), TARGET :: S_DUMMY_ARG(1)
      COMPLEX(kind=8), POINTER, DIMENSION(:) :: S_PTR_ARG
      TYPE(PIVNUL_LIST_STRUCT_T)     :: PIVNUL_LIST_STRUCT
C     Arrowheads
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INTARR
      COMPLEX(kind=8), POINTER, DIMENSION(:) :: DBLARR
C     (pointer to point on used-data in some cases--elt-entry)
      DOUBLE PRECISION TMPTIME
      INTEGER  NOMP
      INTEGER NB_THREADS
      DOUBLE PRECISION TIMEAVG, TIMEMAX,
     &                 FLOPAVG, FLOPMAX
      DOUBLE PRECISION TMPFLOP
#if defined(USE_GPU)
      INTEGER :: GPUBITMAP, MAX_MPIs_PER_GPU, KEEP416_PROVIDED,
     &           KEEP416_MASK, NDEVICES
      DOUBLE PRECISION :: TIME_INIT_GPU
      DOUBLE PRECISION :: TIME_INIT_GPU_MIN
      DOUBLE PRECISION :: TIME_INIT_GPU_MAX
#endif
      INTEGER NPIV_CRITICAL_PATH, EFF_SIZE_SCHUR
      DOUBLE PRECISION TIME, TIMEET
      DOUBLE PRECISION ZERO, ONE, MONE
      PARAMETER( ZERO = 0.0D0, ONE = 1.0D0, MONE = -1.0D0)
#if defined(SINGLEINDOUBLE)
      REAL(kind(0.E0)), PARAMETER :: LPZERO = 0.0e0
      REAL(kind(0.E0)) :: LPSEUIL, LPSEUIL_LDLT_NIV2
#endif
      COMPLEX(kind=8) CZERO
      PARAMETER( CZERO = (0.0D0, 0.0D0) )
      INTEGER PERLU, TOTAL_MBYTES, K231, K232, K233, BLR_STRAT
      INTEGER, PARAMETER :: IDUMMY = -9999
      LOGICAL, PARAMETER :: BDUMMY =.FALSE.
      INTEGER, PARAMETER :: PANEL_TABSIZE = 20
      INTEGER COLOUR, COMM_FOR_SCALING ! For Simultaneous scaling
      INTEGER(8) :: LIWK, LIWK8
      INTEGER(8) :: LWK, LWK_REAL, LWRKR_TH, LWRKC_TH
      INTEGER    :: NOMP_MAX
C     I_AM_SLAVE: used to determine if proc has the role of a slave
C     WK_USER_PROVIDED is set to true when WK_USER is provided by user
      LOGICAL I_AM_SLAVE, PERLU_ON, WK_USER_PROVIDED, EARLYT3ROOTINS
      LOGICAL PRINT_MAXAVG, PRINT_NODEINFO
      DOUBLE PRECISION :: ANORMINF, SEUIL, SEUIL_LDLT_NIV2, Thresh_Seuil
      DOUBLE PRECISION :: CNTL1, CNTL3, CNTL5, CNTL6, EPS
      INTEGER N, LPN_LIST,POSBUF
      INTEGER, DIMENSION (:), ALLOCATABLE :: ITMP2
      INTEGER I,K
      INTEGER(8) :: ITEMP8
      INTEGER    :: PARPIV_T1
      INTEGER FRONTWISE
C temporary variables for collecting stats from all processors
      INTEGER, PARAMETER :: LR_DKEEPSHIFT=49, LR_TABSIZE=18
      DOUBLE PRECISION :: LR_TAB(LR_TABSIZE), LR_EPSILON
      DOUBLE PRECISION :: TMP_MRY_LU_FR
      DOUBLE PRECISION :: TMP_MRY_LU_LRGAIN
      DOUBLE PRECISION :: TMP_MRY_LU_LRGAIN_mixed
      DOUBLE PRECISION :: TMP_MRY_CB_LRGAIN_mixed
      INTEGER :: KEEP451_SAVE, KEEP452_SAVE, KEEP483_SAVE, KEEP484_SAVE
      INTEGER :: KEEP399_SAVE, KEEP20_SAVE
      DOUBLE PRECISION :: TMP_MRY_CB_FR
      DOUBLE PRECISION :: TMP_MRY_CB_LRGAIN
      DOUBLE PRECISION :: TMP_FLOP_LRGAIN
      DOUBLE PRECISION :: TMP_FLOP_TRSM
      DOUBLE PRECISION :: TMP_FLOP_PANEL
      DOUBLE PRECISION :: TMP_FLOP_FRFRONTS
      DOUBLE PRECISION :: TMP_FLOP_TRSM_FR
      DOUBLE PRECISION :: TMP_FLOP_TRSM_LR
      DOUBLE PRECISION :: TMP_FLOP_UPDATE_FR
      DOUBLE PRECISION :: TMP_FLOP_UPDATE_LR
      DOUBLE PRECISION :: TMP_FLOP_UPDATE_LRLR3
      DOUBLE PRECISION :: TMP_FLOP_COMPRESS
      DOUBLE PRECISION :: TMP_FLOP_DECOMPRESS
      DOUBLE PRECISION :: TMP_FLOP_MIDBLK_COMPRESS
      DOUBLE PRECISION :: TMP_FLOP_FRSWAP_COMPRESS
      DOUBLE PRECISION :: TMP_FLOP_ACCUM_COMPRESS
      DOUBLE PRECISION :: TMP_FLOP_CB_COMPRESS
      DOUBLE PRECISION :: TMP_FLOP_CB_DECOMPRESS
      DOUBLE PRECISION :: TMP_FLOP_FACTO_FR
      INTEGER :: TMP_CNT_NODES
      DOUBLE PRECISION :: TMP_TIME_UPDATE
      DOUBLE PRECISION :: TMP_TIME_UPDATE_LRLR1
      DOUBLE PRECISION :: TMP_TIME_UPDATE_LRLR2
      DOUBLE PRECISION :: TMP_TIME_UPDATE_LRLR3
      DOUBLE PRECISION :: TMP_TIME_UPDATE_FRLR
      DOUBLE PRECISION :: TMP_TIME_UPDATE_FRFR
      DOUBLE PRECISION :: TMP_TIME_COMPRESS
      DOUBLE PRECISION :: TMP_TIME_MIDBLK_COMPRESS
      DOUBLE PRECISION :: TMP_TIME_FRSWAP_COMPRESS
      DOUBLE PRECISION :: TMP_TIME_CB_COMPRESS
      DOUBLE PRECISION :: TMP_TIME_PANEL
      DOUBLE PRECISION :: TMP_TIME_FAC_I
      DOUBLE PRECISION :: TMP_TIME_FAC_MQ
      DOUBLE PRECISION :: TMP_TIME_FAC_SQ
      DOUBLE PRECISION :: TMP_TIME_LRTRSM
      DOUBLE PRECISION :: TMP_TIME_FRTRSM
      DOUBLE PRECISION :: TMP_TIME_FRFRONTS
      DOUBLE PRECISION :: TMP_TIME_LR_MODULE
      DOUBLE PRECISION :: TMP_TIME_DIAGCOPY
      DOUBLE PRECISION :: TMP_TIME_DECOMP
      DOUBLE PRECISION :: TMP_TIME_DECOMP_UCFS
      DOUBLE PRECISION :: TMP_TIME_LRASM_NIV1
      DOUBLE PRECISION :: TMP_TIME_LRASM_LOCASM2
      DOUBLE PRECISION :: TMP_TIME_LRASM_MAPLIG1
      DOUBLE PRECISION :: TMP_TIME_LRASM_CONTRIB2
      DOUBLE PRECISION :: TMP_TIME_FRASM_LOCASM2
      DOUBLE PRECISION :: TMP_TIME_FRASM_MAPLIG1
      DOUBLE PRECISION :: TMP_TIME_FRASM_CONTRIB2
C
C  Workspace
C
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWK
      COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE :: WK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE    :: WK_REAL
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: WRKR_TH, 
     &                                    WRKC_TH
      INTEGER(8), DIMENSION(:), ALLOCATABLE :: IWK8
      INTEGER, DIMENSION(:), ALLOCATABLE :: BURP
      INTEGER, DIMENSION(:), ALLOCATABLE :: BUCP
      INTEGER, DIMENSION(:), ALLOCATABLE :: BURS
      INTEGER, DIMENSION(:), ALLOCATABLE :: BUCS
      INTEGER(8) :: BUREGISTRE(12)
      INTEGER(8) :: BUINTSZ, BURESZ, NZ_loc8
      INTEGER    :: BUJOB
      INTEGER BUMAXMN, M, SCMYID, SCNPROCS
      DOUBLE PRECISION    SCONEERR, SCINFERR
C
C  Parameters arising from the structure
C  =====================================
C
*     Control parameters: see description in ZMUMPSID
      DOUBLE PRECISION,DIMENSION(:),POINTER::RINFO, RINFOG
      DOUBLE PRECISION,DIMENSION(:),POINTER::    CNTL
      INTEGER,DIMENSION(:),POINTER   :: INFOG, KEEP
      INTEGER, DIMENSION(:), POINTER :: IRN_loc_PTR, JCN_loc_PTR
      DOUBLE PRECISION, DIMENSION(:), POINTER   :: COLSCA_PTR, 
     &                                 ROWSCA_PTR
      COMPLEX(kind=8), DIMENSION(:), POINTER:: A_loc_PTR
      INTEGER, TARGET :: DUMMYIRN_loc(1), DUMMYJCN_loc(1)
      DOUBLE PRECISION, TARGET    :: DUMMYSCA(1)
      COMPLEX(kind=8), TARGET :: DUMMYA_loc(1)
      INTEGER,DIMENSION(:),POINTER::ICNTL
      EXTERNAL MUMPS_GET_POOL_LENGTH
      INTEGER MUMPS_GET_POOL_LENGTH, SIZESCAL
      INTEGER(8) :: TOTAL_BYTES
C
C  External references
C  ===================
      INTEGER, EXTERNAL :: MUMPS_NUMROC
#if defined(MPI_TO_K_OMP)
      INTEGER, DIMENSION(:), allocatable:: CORES_LIST
      INTEGER:: NCORES
      INTEGER:: ERROR
      INTEGER:: MY_WORKING_ID
      LOGICAL:: I_AM_HEAD
#if defined(MTKO_VERBOSE)
C     SETMASK:
C                   = 2 after restoring original mask
C                   = 1 when  MUMPS_MPITOOMP_SET_MASK called
C                   = 0 otherwise (print original mask)
C                   
      INTEGER:: thread_num, SETMASK
#endif
#endif
      INTEGER:: NWORKING
      LOGICAL:: MEM_EFF_ALLOCATED
      INTEGER   :: TOTAL_MBYTES_UNDER_L0
      INTEGER(8):: TOTAL_BYTES_UNDER_L0
C  Fwd in facto:
      COMPLEX(kind=8), DIMENSION(:), POINTER :: RHS_MUMPS
      LOGICAL :: RHS_MUMPS_ALLOCATED
      LOGICAL :: DBLARR_ALLOCATED
      INTEGER :: NB_ACTIVE_FRONTS_ESTIM
      INTEGER :: NB_FRONTS_F_ESTIM
      INTEGER :: KEEP_486_FOR_PRINT
C
C     --------------------------
C     Pointers used as shortcuts
C     --------------------------
      RINFO=>id%RINFO
      RINFOG=>id%RINFOG
      CNTL=>id%CNTL
      INFOG=>id%INFOG
      KEEP=>id%KEEP
      ICNTL=>id%ICNTL
      IF (id%KEEP8(29) .NE. 0) THEN
        IRN_loc_PTR=>id%IRN_loc
        JCN_loc_PTR=>id%JCN_loc
        A_loc_PTR=>id%A_loc
      ELSE
        IRN_loc_PTR=>DUMMYIRN_loc
        JCN_loc_PTR=>DUMMYJCN_loc
        A_loc_PTR=>DUMMYA_loc
      ENDIF
      NOMP = 1 
      N = id%N
C     TIMINGS: reset to 0
      id%DKEEP(92)=0.0D0
      id%DKEEP(93)=0.0D0
      id%DKEEP(94)=0.0D0
      id%DKEEP(95)=0.0D0
      id%DKEEP(96)=0.0D0
      id%DKEEP(97)=0.0D0
      id%DKEEP(98)=0.0D0
      id%DKEEP(99)=0.0D0
      id%DKEEP(56)=0.0D0
C     Count of MPI messages: reset to 0
      id%KEEP(266)=0
      id%KEEP(267)=0
C     -------------------------------------
C     Depending on the type of parallelism,
C     the master can now (soon) potentially
C     have the role of a slave
C     -------------------------------------
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     &             ( id%MYID .eq. MASTER .AND.
     &               KEEP(46) .eq. 1 ) )
      LIWK  = 0_8
      LIWK8 = 0_8
C     RR related
      id%KEEP(17)  = 0
      id%INFOG(28) = 0
      id%INFOG(50) = 0
      id%INFOG(51) = 0
C     Number of symmetric swaps
      id%KEEP8(80)=0_8
C     Largest increase of internal panel size
      id%KEEP(425) =0
C     Dynamic memory during process_blocfacto, in number of scalar entries
      id%KEEP8(130) = 0_8 ! instantaneous
      id%KEEP8(131) = 0_8 ! max
      id%KEEP8(132) = 0_8 ! max of max
      id%KEEP8(133) = 0_8 ! sum of max
C     Measure recursivity =max number of simultaneous calls to
C     ZMUMPS_FAC_PROCESS_BLOCFACTO_LDLT
      id%KEEP(174) = 0
      id%KEEP(175) = 0
C
      KEEP20_SAVE = KEEP(20)
      PRINT_MAXAVG = .NOT.(id%NSLAVES.EQ.1 .AND. KEEP(46).EQ.1)
C     Print per node informtation only in case ther are several
C     compute nodes (id%KEEP(412): #MPI procs on comupte node)
      PRINT_NODEINFO = PRINT_MAXAVG .AND. id%NPROCS .NE. id%KEEP(412)
C
C     Related to forward in facto functionality (referred to as "Fwd in facto")
      NULLIFY(RHS_MUMPS)
      NULLIFY(DBLARR)
      RHS_MUMPS_ALLOCATED = .FALSE.
      DBLARR_ALLOCATED = .FALSE.
C     -----------------------------------------------------------------------
C     Set WK_USER_PROVIDED to true when workspace WK_USER is provided by user
C     We can accept WK_USER to be provided on only one proc and 
C     different values of WK_USER per processor
C     WK_USER(LWK_USER): only on working processes
      WK_USER_PROVIDED = (id%LWK_USER.NE.0 .AND. I_AM_SLAVE)
      IF (WK_USER_PROVIDED) THEN
          IF (id%LWK_USER.GT.0) THEN
            id%KEEP8(24) = int(id%LWK_USER,8)
          ELSE
            id%KEEP8(24) = -int(id%LWK_USER,8)* 1000000_8 
          ENDIF
      ELSE
          id%KEEP8(24) = 0_8
      ENDIF
C     Compute sum of LWK_USER provided by user 
      CALL MPI_REDUCE ( id%KEEP8(24), id%KEEP8(124), 1, MPI_INTEGER8,
     &                  MPI_SUM, MASTER, id%COMM, IERR )
C
C     KEEP8(26) might be modified
C       (element entry format) 
C       but need be restore for 
C       future factorisation
C       with different scaling option
C 
      KEEP826_SAVE = id%KEEP8(26)
C     In case of loop on factorization with
C     different scaling options, initialize
C     DKEEP(4:5) to 0.
      id%DKEEP(4)=-1.0D0
      id%DKEEP(5)=-1.0D0
C  Mapping information used during solve. In case of several facto+solve
C  it has to be recomputed. In case of several solves with the same
C  facto, it is not recomputed.
      IF (associated(id%IPTR_WORKING)) THEN
        DEALLOCATE(id%IPTR_WORKING)
        NULLIFY(id%IPTR_WORKING)
      END IF
      IF (associated(id%WORKING)) THEN 
        DEALLOCATE(id%WORKING)
        NULLIFY(id%WORKING)
      END IF
C
C  Units for printing
C  MP: diagnostics
C  LP: errors
C
      LP  = ICNTL( 1 )
      MP  = ICNTL( 2 )
      MPG = ICNTL( 3 )
      LPOK    = ((LP.GT.0).AND.(id%ICNTL(4).GE.1))
      PROK    = ((MP.GT.0).AND.(id%ICNTL(4).GE.2))
      PROKG   = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
      PROKG   = (PROKG.AND.(id%ICNTL(4).GE.2))
      IF ( PROK ) WRITE( MP, 130 )
      IF ( PROKG ) WRITE( MPG, 130 )
C
C  Prepare work for out-of-core
C
      IF (id%MYID .EQ. MASTER .AND. KEEP(201) .NE. -1) THEN
C       Note that if KEEP(201)=-1, then we have decided
C       at analysis phase that factors will not be stored
C       (neither in memory nor on disk). In that case,
C       ICNTL(22) is ignored.
C       -- ICNTL(22) must be set before facto phase 
C          (=1 OOC on; =0 OOC off)
C          and cannot be changed for subsequent solve phases.
        KEEP(201)=id%ICNTL(22)
        IF (KEEP(201) .EQ. 1) THEN
#         if defined(OLD_OOC_NOPANEL)
            KEEP(201)=2
#         else
            KEEP(201)=1
#         endif
        ELSE
          id%KEEP(201)=0
        ENDIF
      ENDIF
C
C     ----------------------
C     Broadcast ICNTL(49)
      IF (id%MYID.EQ.MASTER) THEN 
         ICNTL49_LOC=id%ICNTL(49)
C        out of range treated as 0
         IF ( (ICNTL49_LOC.GT.2).or.(ICNTL49_LOC.LT.0) )
     &   ICNTL49_LOC = 0
      ENDIF
      CALL MPI_BCAST( ICNTL49_LOC, 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
C     ----------------------
C
C     Broadcast few other KEEP entries that have been decoded
C     and are defined for facto:
C     ----------------------
      CALL MPI_BCAST( KEEP(12), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(19), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(21), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(201), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(459), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(460), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      IF ( KEEP(459) .GE. PANEL_TABSIZE ) THEN
        IF ( LPOK ) THEN
        WRITE(LP,'(A,I4,A,I3)') " ** WARNING ** KEEP(459)=",KEEP(459),
     &  " too large, resetting to",PANEL_TABSIZE-1
        ENDIF
        KEEP(459) = PANEL_TABSIZE - 1
      ENDIF
      PERLU = KEEP(12)
      IF (id%MYID.EQ.MASTER) THEN
C     {
C       KEEP(50)  case
C       ==============
C
C       KEEP(50)  = 0 : matrix is unsymmetric
C       KEEP(50) /= 0 : matrix is symmetric
C       KEEP(50) = 1 : Ask L L^T on the root. Matrix is PSD.
C       KEEP(50) = 2 : Ask for L U on the root
C       KEEP(50) = 3 ... L D L^T ??
C     
        CNTL1 = id%CNTL(1)
C       ---------------------------------------
C       For symmetric (non general) matrices 
C       set (directly) CNTL1 = 0.0
C       ---------------------------------------
        KEEP(17)=0
C       Automatic choice if CNTL(1)<0
C       For rank-revealing (KEEP(19).GT.0) then
C       set CNTL1=0.1 even if SYM=1 
        IF (CNTL1.LT.ZERO) THEN
C        automatic choice
         IF (KEEP(19).GT.0) THEN
          CNTL1=0.1D0
         ELSE IF (KEEP(50).EQ.1) THEN
          CNTL1=ZERO
         ELSE
          CNTL1=0.01D0
         ENDIF
        ENDIF
        IF ( KEEP(50) .eq. 1 ) THEN
          IF (CNTL1 .ne. ZERO ) THEN
            IF ( PROKG ) THEN
              WRITE(MPG,'(A)')
     & '** Warning : SPD solver called, resetting CNTL(1) to 0.0D0'
            END IF
          END IF
          CNTL1 = ZERO
        END IF
C       CNTL1 threshold value must be between 
C       0.0 and 1.0 (for SYM=0) and 0.5 (for SYM=1,2)
        IF (CNTL1.GT.ONE)   CNTL1=ONE
        IF (CNTL1.LT.ZERO)  CNTL1=ZERO
        IF (KEEP(50).NE.0.AND.CNTL1.GT.0.5D0) THEN
          CNTL1 = 0.5D0
        ENDIF
        PARPIV_T1 = id%KEEP(268)
        IF (PARPIV_T1.EQ.77) THEN
         PARPIV_T1 = 0 
#if defined(USE_GPU)
         IF (KEEP(417).GT.0) THEN
          PARPIV_T1 = -2
         ENDIF
#endif
#if defined(__ve__)
         PARPIV_T1 = -2
#endif
        ENDIF
#if !defined(USE_GPU)
        IF (PARPIV_T1.EQ.-3) THEN
          PARPIV_T1 = 0
        ENDIF
#endif
        IF ((PARPIV_T1.LT.-3).OR.(PARPIV_T1.GT.1)) THEN
C        out of range values
         PARPIV_T1 =0
        ENDIF
C       note that KEEP(50).EQ.1 => CNTL1=0.0
        IF (CNTL1.EQ.0.0D0.OR.(KEEP(50).eq.1)) PARPIV_T1 = 0 
C
        IF (PARPIV_T1.EQ.-2) THEN
         IF (KEEP(19).NE.0) THEN
C         switch off PARPIV_T1 if RR activated
C         but do NOT switch off PARPIV_1 with null pivot detection
          PARPIV_T1 = 0
         ENDIF
        ENDIF
        id%KEEP(269) = PARPIV_T1
C       Save few KEEP related to LRMIXED now,
C       they may be modified later
        KEEP451_SAVE = KEEP(451)
        KEEP452_SAVE = KEEP(452)
        KEEP483_SAVE = KEEP(483) 
        KEEP484_SAVE = KEEP(484) 
C     }
      ENDIF
      CALL MPI_BCAST(CNTL1, 1, MPI_DOUBLE_PRECISION,
     &             MASTER, id%COMM, IERR)
      CALL MPI_BCAST( KEEP(269), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
#if defined(SINGLEINDOUBLE)
C     single precision factorization? Set KEEP(148), etc.
      CALL MUMPS_ICNTL47_SETTINGS( 2, id%ICNTL, id%KEEP,
     &id%KEEP8, id%SIZE_SCHUR, id%MYID, id%COMM, id%INFO, LP, LPOK )
      IF (id%INFO(1).LT.0) GOTO 530
#endif
      IF (id%MYID.EQ.MASTER) THEN
C{
C         OMP parallelization of arrowheads
C         out of range are treated as zero
          IF (KEEP(399).LT.-1) KEEP(399)=-1
          KEEP399_SAVE = KEEP(399)
          IF (KEEP(399).EQ.-1) THEN
            IF ((KEEP(54).EQ.0).AND.(id%NPROCS.GT.1)) THEN
              KEEP(399) = 1
            ELSE
              KEEP(399) = 3
            ENDIF
          ENDIF
#if defined(PCPRET)
C         new multithreaded >=2 algo does not compile on PCPRET
          KEEP(399) = 1
#endif
C         -----------------------------------------------------
C         Decoding of ICNTL(35) for factorization: same as
C         at analysis except that we store a copy of ICNTL(35)
C         in KEEP(486) instead of KEEP(494) and need to check
C         compatibility of KEEP(486) and KEEP(494): If LR was
C         not activated during analysis, it cannot be activated 
C         at factorization.
C         ------------------------------------------------------
          id%KEEP(486) = id%ICNTL(35)
#if defined(SINGLEINDOUBLE)
          IF (id%KEEP(473) .GE. 2) id%KEEP(486) = 1
#endif
          IF (id%KEEP(486).EQ.1) THEN
C         -- Automatic BLR option setting
           id%KEEP(486)= 2
          ENDIF
          IF ( id%KEEP(486).EQ.4) id%KEEP(486)=0
          IF ((id%KEEP(486).LT.0).OR.(id%KEEP(486).GT.4)) THEN
C           Out of range values treated as 0
            id%KEEP(486) = 0
          ENDIF
          IF ((KEEP(486).NE.0).AND.(KEEP(494).EQ.0)) THEN
C           To activate BLR during factorization,
C           ICNTL(35) must have been set at analysis.
            IF (LPOK)  THEN
              WRITE(LP,'(A)') 
     &      " *** Error with BLR setting "
              WRITE(LP,'(A)') " *** BLR was not activated during ",
     &      " analysis but is requested during factorization."
            ENDIF
            id%INFO(1)=-54
            id%INFO(2)=0
            GOTO 105
          ENDIF
C         Save value of KEEP(486) before possibly
C         forcing it to 3 in case of discard factors
          KEEP_486_FOR_PRINT=KEEP(486)
          IF (KEEP(201) .EQ. -1 .AND. KEEP(486) .NE.0) THEN
            KEEP(486) = 3
          ENDIF
          KEEP464COPY  = id%ICNTL(38)
          IF (KEEP464COPY.LT.0.OR.KEEP464COPY.GT.1000) THEN
C          Out of range values treated as 1000
           KEEP464COPY = 1000
          ENDIF
          IF (id%KEEP(461).LT.1) THEN
            id%KEEP(461) = 10
          ENDIF
          IF (id%KEEP(462).LT.1) THEN
            id%KEEP(462) = 10
          ENDIF
          KEEP465COPY = id%ICNTL(39)
          IF (KEEP465COPY.LT.0.OR.(KEEP465COPY.GT.1000)) THEN
C          Out of range values treated as 1000
           KEEP465COPY = 1000
          ENDIF
          ICNTL36_LOC = id%ICNTL(36)
C         Out of range values
          IF (ICNTL36_LOC.NE.0) ICNTL36_LOC= 0
          IF (ICNTL36_LOC.EQ.1.OR.ICNTL36_LOC.EQ.3) THEN
            IF (CNTL1.EQ.ZERO .OR. KEEP(468).LE.1) THEN
              KEEP(475) = 3
            ELSE IF ( (KEEP(269).GT.0).OR. (KEEP(269).EQ.-2)) THEN
              KEEP(475) = 2
            ELSE IF (KEEP(468).EQ.2) THEN
              KEEP(475) = 2
            ELSE
              KEEP(475) = 1
            ENDIF
          ELSE
            KEEP(475) = 0
          ENDIF
          KEEP(481)=0
          IF (ICNTL36_LOC.LT.0 .OR. ICNTL36_LOC.GE.2) THEN
C           Only options 1 and 2 are allowed
            KEEP(475) = 0
          ENDIF
C         K489 is set according to ICNTL(37)
          IF (id%ICNTL(37).EQ.0.OR.id%ICNTL(37).EQ.1) THEN
            KEEP(489) = id%ICNTL(37)
          ELSE
C           Other values treated as zero
            KEEP(489) = 0
          ENDIF
          IF (KEEP(79).GE.1) THEN
C          CompressCB incompatible with type4,5,6 nodes
           KEEP(489)=0
          ENDIF
          IF (id%ICNTL(40).EQ.0.OR.id%ICNTL(40).EQ.1) THEN
            KEEP(450) = id%ICNTL(40)
          ELSE
C           Other values treated as zero
            KEEP(450) = 0
          ENDIF
C         Saves are done earlier to ensure they're always done
C         before restoration (case of errors with "GOTO")
C         KEEP451_SAVE = KEEP(451)
C         KEEP452_SAVE = KEEP(452)
C         KEEP483_SAVE = KEEP(483)
C         KEEP484_SAVE = KEEP(484)
          IF (KEEP(451).EQ.-1) THEN
C           Use all storage precisions by default
            KEEP(451) = 1
#if defined(NOMIXHPC)
#if defined(SINGLEINDOUBLE)
C           KEEP(148)=0 means that we are really in double
C           (in case of single computations, do not reset)
            IF (KEEP(148).EQ.0) THEN
#endif
C             ... except in double for LRMIXED_CALC 
C             (use only double and single)
              IF (KEEP(448) .GT. 0) THEN
                KEEP(451) = 0
              ENDIF
#if defined(SINGLEINDOUBLE)
            ENDIF
#endif
#endif
          ENDIF
          IF (KEEP(452).EQ.-1) THEN
            KEEP(452) = KEEP(450)
          ENDIF
C         id%KEEP(476) \in [1,100] 
          IF ((id%KEEP(476).GT.100).OR.(id%KEEP(476).LT.1)) THEN
              id%KEEP(476)=  50
          ENDIF
C         id%KEEP(477) \in [1,100] 
          IF ((id%KEEP(477).GT.100).OR.(id%KEEP(477).LT.1)) THEN
              id%KEEP(477)=  100
          ENDIF
C         id%KEEP(483) \in [1,100] 
          IF ((id%KEEP(483).GT.100).OR.(id%KEEP(483).LT.1)) THEN
              id%KEEP(483)= 80
          ENDIF
C         id%KEEP(484) \in [1,100] 
          IF ((id%KEEP(484).GT.100).OR.(id%KEEP(484).LT.1)) THEN
              id%KEEP(484)=  80
          ENDIF
C         id%KEEP(480)=0,2,3,4,5,6
          IF ((id%KEEP(480).GT.6).OR.(id%KEEP(480).LT.0)
     &                           .OR.(id%KEEP(480).EQ.1)) THEN
               id%KEEP(480)=0
          ENDIF
C         id%KEEP(473)=0,1,2 or 3
          IF (id%KEEP(473).LT.0.OR.id%KEEP(473).GT.3) THEN
             id%KEEP(473)=0
          ENDIF
C         id%KEEP(474)=0,1,2,3
          IF ((id%KEEP(474).GT.3).OR.(id%KEEP(474).LT.0)) THEN
               id%KEEP(474)=0
          ENDIF
C         id%KEEP(479)>0
          IF (id%KEEP(479).LE.0) THEN
               id%KEEP(479)=1
          ENDIF
          IF (id%KEEP(474).NE.0.AND.id%KEEP(480).EQ.0) THEN
              id%KEEP(474) = 0
          ENDIF
          IF (id%KEEP(478).NE.0.AND.id%KEEP(480).LT.4) THEN
              id%KEEP(478) = 0
          ENDIF
          IF (id%KEEP(480).GE.5 .OR.
     &           (id%KEEP(480).NE.0.AND.id%KEEP(474).EQ.3)) THEN
            IF (id%KEEP(475).LT.2) THEN
C             Reset to 3 if 5 or to 4 if 6
              id%KEEP(480) = id%KEEP(480) - 2
              write(*,*) ' Resetting KEEP(480) to ', id%KEEP(480)
            ENDIF
          ENDIF
 105    CONTINUE
C}
      ENDIF  ! id%MYID .EQ. MASTER
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                     id%COMM, id%MYID )
C
      IF (id%INFO(1).LT.0) GOTO 530
#if defined(SINGLEINDOUBLE)
      IF ( KEEP(148) .EQ. 0 ) THEN
#endif
        EPS = epsilon ( ZERO )
#if defined(SINGLEINDOUBLE)
      ELSE
        EPS = epsilon( LPZERO )
      ENDIF
#endif
      CALL MPI_BCAST( KEEP(281), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(399), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(473), 14, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      IF (KEEP(486).NE.0) THEN
        CALL MPI_BCAST( KEEP(489), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
        CALL MPI_BCAST( KEEP464COPY, 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
        CALL MPI_BCAST( KEEP465COPY, 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
C       BCASTS KEEP(448), 450, 451, 452, 453, 454
        CALL MPI_BCAST( KEEP(448), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
        CALL MPI_BCAST( KEEP(450), 5, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      ELSE
C       When BLR is off during facto, 
C       make sure that adaptive BLR (ICNTL(40) "copied" in KEEP(450)
C       is considered off.
C       Doing so that it can be tested also during solve by all mpi 
C       procs even when BLR is off
        KEEP(450)=0
      ENDIF
      IF (id%MYID.EQ.MASTER) THEN
        IF (KEEP(217).GT.2.OR.KEEP(217).LT.0) THEN
          KEEP(217)=0
        ENDIF
        KEEP(214)=KEEP(217)
        IF (KEEP(214).EQ.0) THEN
          IF (KEEP(201).NE.0) THEN ! OOC or no factors
            KEEP(214)=1
          ELSE
            KEEP(214)=2
          ENDIF
          IF (KEEP(486).EQ.2) THEN
            KEEP(214)=1
          ENDIF
        ENDIF
      ENDIF
      CALL MPI_BCAST( KEEP(214), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      IF (KEEP(201).NE.0) THEN
C       -- Low Level I/O strategy
        CALL MPI_BCAST( KEEP(99), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
        CALL MPI_BCAST( KEEP(205), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
        CALL MPI_BCAST( KEEP(211), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
        CALL MPI_BCAST( KEEP(255), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      ENDIF
C     Fwd in facto: explicitly forbid
C     sparse RHS and A-1 computation
      IF (id%KEEP(252).EQ.1 .AND. id%MYID.EQ.MASTER) THEN
        IF (id%ICNTL(20).EQ.1) THEN ! out-of-range => 0
C         NB: in doc ICNTL(20) only accessed during solve
C         In practice, will have failed earlier if RHS not allocated.
C         Still it looks safer to keep this test.
          id%INFO(1)=-43
          id%INFO(2)=20
          IF (LPOK) WRITE(LP,'(A)')
     &       ' ERROR: Sparse RHS is incompatible with forward',
     &       ' performed during factorization (ICNTL(32)=1)'
        ELSE IF (id%ICNTL(30).NE.0) THEN ! out-of-range => 1
          id%INFO(1)=-43
          id%INFO(2)=30
          IF (LPOK) WRITE(LP,'(A)')
     &       ' ERROR: A-1 functionality incompatible with forward',
     &       ' performed during factorization (ICNTL(32)=1)'
        ELSE IF (id%ICNTL(9) .NE. 1) THEN
          id%INFO(1)=-43
          id%INFO(2)=9
          IF (LPOK) WRITE(LP,'(A)')
     &       ' ERROR: Transpose system (ICNTL(9).NE.0) not ',
     &       ' compatible with forward performed during',
     &       ' factorization (ICNTL(32)=1)'
        ENDIF
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
C
      IF (id%INFO(1).LT.0) GOTO 530
C
C     The memory allowed is given by ICNTL(23) in Mbytes
C     0 means that nothing is provided.
C     Save memory available, ICNTL(23) in KEEP8(4)
C
      IF ( ICNTL(23) .GT. 0 ) THEN
        ITMP = 1
      ELSE
        ITMP = 0
      ENDIF
      CALL MPI_ALLREDUCE( ITMP, JTMP, 1, MPI_INTEGER,
     &                 MPI_SUM, id%COMM, IERR)
      IF ( id%MYID.EQ.MASTER ) THEN
C       Negative values considered 0
        ITMP = max(ICNTL(23),0)
      END IF
      CALL MPI_BCAST( ITMP, 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
C     JTMP: nb of procs with nonzero ICNTL(23)
C     ITMP: value of ICNTL(23) on master
      IF ( ITMP .GT. 0 .AND. JTMP .EQ. 1 ) THEN
C       ICNTL(23)>0 only on master
      ELSE
C       Local values of ICNTL(23) are used, note that
C       they could all be zeros
        ITMP = ICNTL(23)
      ENDIF
C 
      ITMP8 = int(ITMP, 8)
      id%KEEP8(4) = ITMP8 * 1000000_8   ! convert to nb of bytes
C     Compute \sum of memories allowed
      CALL MPI_REDUCE( id%KEEP8(4), ITMP8, 1, MPI_INTEGER8,
     &                  MPI_SUM, MASTER, id%COMM, IERR )
      ITMP8 = ITMP8 / 1000000_8 ! Use to print \sum_{ICNTL(23)}
      IF ( PROKG ) THEN
#if defined(MPI_TO_K_OMP)
        IF (id%KEEP(381).gt.1) THEN
          NWORKING = id%NSLAVES/id%KEEP(381)
        ELSE
#endif
          NWORKING = id%NSLAVES
#if defined(MPI_TO_K_OMP)
        ENDIF
#endif
          CALL MUMPS_SETI8TOI4( id%KEEP8(129), TMP_INFOG_4)
          WRITE( MPG, 172 ) 
     &    NWORKING, id%ICNTL(22), KEEP_486_FOR_PRINT,
     &    KEEP(489),
     &    KEEP(450),
#if defined(SINGLEINDOUBLE)
     &    KEEP(148),
#endif
     &    id%ICNTL(49), 
     &    id%KEEP(19),
     &    KEEP(12), 
     &    id%KEEP8(111), TMP_INFOG_4, KEEP(127), KEEP(28),
     &    id%KEEP8(4)/1000000_8, ITMP8, id%KEEP8(124), CNTL1
          IF (KEEP(252).GT.0) 
     &    WRITE(MPG,173) KEEP(253)
          IF (KEEP(269).NE.0) 
     &    WRITE(MPG,176) KEEP(269)
      ENDIF
      IF (KEEP(201).LE.0) THEN
C       In-core version or no factors
        KEEP(IXSZ)=XSIZE_IC
      ELSE IF (KEEP(201).EQ.2) THEN
C       OOC version, no panels
        KEEP(IXSZ)=XSIZE_OOC_NOPANEL
      ELSE IF (KEEP(201).EQ.1) THEN
C     Panel versions:
        IF (KEEP(50).EQ.0) THEN
          KEEP(IXSZ)=XSIZE_OOC_UNSYM
        ELSE
          KEEP(IXSZ)=XSIZE_OOC_SYM
        ENDIF
      ENDIF
      IF ( KEEP(486) .NE. 0 ) THEN !LR is activated
C       Stats initialization for LR
        CALL INIT_STATS_GLOBAL()
      END IF
C     Memory management: allocate id%S etc. from C or Fortran?
      id%KEEP(430) = 0
#if defined(MUMPS_MALLOC_FROM_C)
      id%KEEP(430) = 1
#endif
#if defined(USE_GPU)
C{
      CALL MPI_BCAST( id%KEEP(415), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
C     BCAST KEEP(417), KEEP(418), KEEP(419)
      CALL MPI_BCAST( id%KEEP(417), 3, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%KEEP(422), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
C     BCAST KEEP(431),KEEP(432),KEEP(433),KEEP(434),KEEP(435)
      CALL MPI_BCAST( id%KEEP(431), 5, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      IF (id%KEEP(417).GT.0) then
#if defined(MPI_TO_K_OMP)
        IF (id%KEEP(381) .GT. 1) THEN
          id%INFO(1)=-800
          id%INFO(2)=17
          IF (LPOK)
     &    WRITE(LP,*)
     &   'ERROR: MPI_TO_K_OMP feature currently incompatible with GPU'
          CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &       id%COMM, id%MYID )
          IF ( id%INFO(1) .LT. 0 ) GOTO 530
        ENDIF
#endif
        id%KEEP(428)=0
        CALL ZMUMPS_GET_NDEVICES_C( NDEVICES )
        IF (id%MYID .EQ. MASTER) THEN
          KEEP416_PROVIDED = id%KEEP(416) ! 0 means automatic
          IF (KEEP416_PROVIDED .EQ. 0) THEN
            KEEP416_MASK = NDEVICES
            KEEP416_MASK = ishftc(1, KEEP416_MASK)-1
          ENDIF
        ENDIF
        CALL MPI_BCAST( KEEP416_PROVIDED, 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
        IF (KEEP416_PROVIDED.EQ.0) THEN
          CALL MPI_BCAST( KEEP416_MASK, 1, MPI_INTEGER,
     &                    MASTER, id%COMM, IERR )
        ELSE IF (KEEP416_PROVIDED > 0) THEN
          KEEP416_MASK = id%KEEP(416)
        ENDIF
        IF ( KEEP416_PROVIDED < 0 ) THEN
C         Special case where KEEP416_PROVIDED encodes
C         the maximum number of MPI per GPU.
          IF ( NDEVICES < KEEP(417) ) THEN
C           NDEVICES not big enough:
            id%INFO(1) = -100
            id%INFO(2) = NDEVICES
          ELSE
C           Set GPUBITMAP and MAX_MPI_PER_GPU
            GPUBITMAP = 2**KEEP(417)-1
            MAX_MPIs_PER_GPU = -KEEP416_PROVIDED
          ENDIF
        ELSE
          CALL MUMPS_GPU_DEF_BITMAP( KEEP416_MASK, id%KEEP(417),
     &                      id%MYID, id%NPROCS, id%COMM,
C    &                              id%KEEP(46),
C    &                      id%MYID_NODES, id%NSLAVES, id%COMM,
     &                      GPUBITMAP, MAX_MPIs_PER_GPU,
     &                      id%ICNTL, id%INFO )
        ENDIF
        CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
        IF ( id%INFO(1) .LT. 0 ) GOTO 530
        IF ( id%KEEP(436).EQ. 1 ! unified requested
     &       .AND. ( id%KEEP(438) .LT. 0 .OR.
     &               id%KEEP(438) .EQ. 0 .OR.
     &               id%KEEP(438) .EQ. 2 )) THEN
            IF ( PROKG ) THEN
              WRITE(MPG,'(A)')
     & '** ERROR : Unified memory requested but not supported'
            END IF
            CALL MUMPS_ABORT()
        ENDIF
        CALL MUMPS_SECDEB(TIME_INIT_GPU)
#if defined(SINGLEINDOUBLE)
        IF (KEEP(148) .EQ. 0) THEN
#endif
        CALL ZMUMPS_INIT_GPU_C(id%KEEP(415), id%KEEP(417), id%KEEP(418),
     &                         id%KEEP(419), id%KEEP(422), id%KEEP(430),
     &                         id%KEEP(436),
     &                         GPUBITMAP, MAX_MPIs_PER_GPU,
     &                         id%INFO )
#if defined(SINGLEINDOUBLE)
        ELSE
        CALL CMUMPS_INIT_GPU_C(id%KEEP(415), id%KEEP(417),id%KEEP(418),
     &                         id%KEEP(419), id%KEEP(422), id%KEEP(430),
     &                         id%KEEP(436),
     &                         GPUBITMAP, MAX_MPIs_PER_GPU,
     &                         id%INFO )
        ENDIF
#endif
        CALL MUMPS_SECFIN(TIME_INIT_GPU)
        CALL MPI_REDUCE(TIME_INIT_GPU,TIME_INIT_GPU_MAX,1,
     &      MPI_DOUBLE_PRECISION,
     &      MPI_MAX,MASTER,id%COMM,IERR)
        CALL MPI_REDUCE(TIME_INIT_GPU,TIME_INIT_GPU_MIN,1,
     &      MPI_DOUBLE_PRECISION,
     &      MPI_MIN,MASTER,id%COMM,IERR) 
        IF ( PROK ) THEN
C          print on all procs
           WRITE(MP,'(A,I16)') 
     &     ' ICNTL(51) Number of GPUs                   =',KEEP(417)
             WRITE(MP,'(A,I16)') 
     &     ' Set of GPU devices to be used              =',KEEP(416)
             WRITE(MP,'(A,I16)') 
     &     ' GPU Arith. Intensity threshold             =',KEEP(418)
            WRITE(MP,'(A,I16)')
     &     ' Memory pinning (0=OFF)                     =',KEEP(422)
           WRITE(MP,'(A,I16)')
     &     ' Max number of MPI per GPU                  =',
     &     MAX_MPIs_PER_GPU
           WRITE(MP,'(A,I16,I4)')
     &     ' Local GPU bit set, id%MYID                 =',GPUBITMAP,
     &                                                     id%MYID
            WRITE(MP,'(A,I16)')
     &     ' MPI setting for larger granularity (0=OFF) =',KEEP(198)
           WRITE(MP,'(A,F12.4)')
     &     ' Local time spent in MUMPS_INIT_GPU_C       =',TIME_INIT_GPU
         ENDIF
         IF (PROKG .AND. .NOT. (PROK .AND. MPG .EQ. MP)) THEN
C          print on master
           WRITE(MPG,'(A,I16)') 
     &     ' ICNTL(51) Offload on GPUs                  =',KEEP(417)
             WRITE(MPG,'(A,I16)') 
     &     ' Set of GPU devices to be used              =',KEEP(416)
             WRITE(MPG,'(A,I16)') 
     &     ' GPU Arith. Intensity threshold             =',KEEP(418)
            WRITE(MPG,'(A,I16)')
     &     ' Memory pinning (0=OFF)                     =',KEEP(422)
           WRITE(MPG,'(A,I16)')
     &     ' Max number of MPI per GPU                  =',
     &     MAX_MPIs_PER_GPU
            WRITE(MPG,'(A,I16)')
     &     ' MPI setting for larger granularity (0=OFF) =',KEEP(198)
         ENDIF
         IF (PROKG) THEN
           IF (id%NSLAVES.GT.1) THEN
C            Print min and max time
             WRITE(MPG,'(A,F12.4)')
     &       ' Max time spent in MUMPS_INIT_GPU_C         =',
     &       TIME_INIT_GPU_MAX
             WRITE(MPG,'(A,F12.4)')
     &       ' Min time spent in MUMPS_INIT_GPU_C         =',
     &       TIME_INIT_GPU_MIN
           ELSE 
C            Only 1 MPI worker 
             WRITE(MPG,'(A,F12.4)')
     &       ' Time spent in MUMPS_INIT_GPU_C             =',
     &       TIME_INIT_GPU_MIN
           ENDIF
         ENDIF
         CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0) GOTO 530
      ENDIF
C}
#endif
C
*     **********************************
*     Begin intializations regarding the
*     computation of the determinant
*     **********************************
      IF (id%MYID.EQ.MASTER) KEEP(258)=ICNTL(33)
      CALL MPI_BCAST(KEEP(258), 1, MPI_INTEGER,
     &               MASTER, id%COMM, IERR)
      IF (KEEP(258) .NE. 0) THEN
        KEEP(259) = 0      ! Initial exponent of the local determinant
        KEEP(260) = 1      ! Number of permutations
        id%DKEEP(6)  = 1.0D0  ! real part of the local determinant
        id%DKEEP(7)  = 0.0D0  ! imaginary part of the local determinant
      ENDIF
*     ********************************
*     End intializations regarding the
*     computation of the determinant
*     ********************************
C
         CALL MUMPS_STOP_ON_USER_REQUEST(id%KEEP,id%KEEP8, id%ICNTL,
     &   id%INFO, id%MYID)
         CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0) GOTO 530
*     **********************
*     Begin of Scaling phase
*     **********************
C
C     SCALING MANAGEMENT
C     * Options 1, 3, 4 centralized only
C  
C     * Options 7, 8  : also works for distributed matrix
C
C     At this point, we have the scaling arrays allocated
C     on the master. They have been allocated on the master
C     inside the main MUMPS driver.
C
      CALL MPI_BCAST(KEEP(52), 1, MPI_INTEGER,
     &               MASTER, id%COMM, IERR)
      LSCAL = ((KEEP(52) .GT. 0) .AND. (KEEP(52) .LE. 8))
      IF (LSCAL) THEN
C
        IF ( id%MYID.EQ.MASTER ) THEN
          CALL MUMPS_SECDEB(TIMEET)
        ENDIF
C       -----------------------
C       Retrieve parameters for
C       simultaneous scaling
C       -----------------------
        IF (KEEP(52) .EQ. 7) THEN 
C       -- Cheap setting of SIMSCALING (it is the default in 4.8.4)
           K231= KEEP(231)
           K232= KEEP(232)
           K233= KEEP(233)
        ELSEIF (KEEP(52) .EQ. 8) THEN
C       -- More expensive setting of SIMSCALING (it was the default in 4.8.1,2,3)
           K231= KEEP(239)
           K232= KEEP(240)
           K233= KEEP(241)
        ENDIF
        CALL MPI_BCAST(id%DKEEP(3),1,MPI_DOUBLE_PRECISION,MASTER,
     &       id%COMM,IERR)
C
        IF ( ((KEEP(52).EQ.7).OR.(KEEP(52).EQ.8)) .AND. 
     &       KEEP(54).NE.0 ) THEN
C{        ------------------------------
C         Scaling for distributed matrix
C         We need to allocate scaling
C         arrays on all processors, not
C         only the master.
C         ------------------------------
           IF ( id%MYID .NE. MASTER ) THEN
              IF ( associated(id%COLSCA))
     &             DEALLOCATE( id%COLSCA )
              IF ( associated(id%ROWSCA))
     &             DEALLOCATE( id%ROWSCA )
            ALLOCATE( id%COLSCA(N), stat=IERR)
            IF (IERR .GT.0) THEN
               id%INFO(1)=-13
               id%INFO(2)=N
            ENDIF
            ALLOCATE( id%ROWSCA(N), stat=IERR)
            IF (IERR .GT.0) THEN
               id%INFO(1)=-13
               id%INFO(2)=N
            ENDIF
         ENDIF
         M = N
         BUMAXMN=M
         IF(N > BUMAXMN) BUMAXMN = N
         LIWK = 4_8*int(BUMAXMN,8)
         ALLOCATE (IWK(LIWK), BURP(M), BUCP(N),
     &            BURS(2* id%NPROCS), BUCS(2* id%NPROCS),
     &            stat=allocok)
         IF (allocok > 0) THEN
            id%INFO(1)=-13
            CALL MUMPS_SET_IERROR(LIWK+int(M,8)+int(N,8)+
     &      4_8*int(id%NPROCS,8) , id%INFO(2))
         ENDIF
C        --- Propagate enventual error
         CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF (id%INFO(1).LT.0) GOTO 517
C        -- estimation of memory and construction of partvecs
         BUJOB = 1
C        -- LWK not used
         LWK_REAL = 1_8
         LWRKR_TH = 1_8
         LWRKC_TH = 1_8
         NOMP_MAX   = 1
         ALLOCATE(WK_REAL(LWK_REAL),
     &            WRKR_TH(LWRKR_TH,NOMP_MAX),
     &            WRKC_TH(LWRKC_TH,NOMP_MAX),
     &            stat=allocok)
         IF (allocok > 0) THEN
            id%INFO(1)=-13
            CALL MUMPS_SET_IERROR(
     &             LWK_REAL+ 
     &             LWRKR_TH*int(NOMP_MAX,8)+
     &             LWRKC_TH*int(NOMP_MAX,8),
     &             id%INFO(2))
         ENDIF
C        --- Propagate enventual error
         CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF (id%INFO(1).LT.0) GOTO 517
         CALL ZMUMPS_SIMSCALEABS(
     &        IRN_loc_PTR(1), JCN_loc_PTR(1), A_loc_PTR(1),
     &        id%KEEP8(29),
     &        M, N,  id%NPROCS, id%MYID, id%COMM,
     &        BURP, BUCP,
     &        BURS, BUCS, BUREGISTRE,
     &        IWK, LIWK,
     &        BUINTSZ, BURESZ, BUJOB,
     &        id%ROWSCA(1), id%COLSCA(1), WK_REAL, LWK_REAL, 
     &        WRKR_TH, LWRKR_TH, 
     &        WRKC_TH, LWRKC_TH,
     &        NOMP_MAX, id%KEEP(50),
     &        K231, K232, K233, 
     &        id%DKEEP(3),
     &        SCONEERR, SCINFERR)
         IF(LIWK < BUINTSZ) THEN
            DEALLOCATE(IWK)
            LIWK = BUINTSZ
            ALLOCATE(IWK(LIWK), stat=allocok)
            IF (allocok > 0) THEN
               id%INFO(1)=-13
               CALL MUMPS_SET_IERROR(LIWK,id%INFO(2))
            ENDIF
         ENDIF
         DEALLOCATE(WK_REAL, WRKR_TH, WRKC_TH)
         LWK_REAL     = BURESZ
C
C        -- Set NOMP_MAX from KEEP(281)
         CALL ZMUMPS_SET_NOMP_MAX(id%KEEP(281), id%KEEP(361), 
     &                            N, NOMP_MAX)
C
         IF (NOMP_MAX.LE.1) THEN
C          temp array per thread not used
           LWRKR_TH = 1
           LWRKC_TH = 1
         ELSE
           LWRKR_TH = N
           IF (id%KEEP(50).NE.0) THEN
C            WRKC_TH not used on symmetric matrices
             LWRKC_TH = 1
           ELSE
             LWRKC_TH = N
           ENDIF
         ENDIF
         ALLOCATE(WK_REAL(LWK_REAL),
     &            WRKR_TH(LWRKR_TH,max(NOMP_MAX,1)),
     &            WRKC_TH(LWRKC_TH,max(NOMP_MAX,1)),
     &            stat=allocok)
         IF (allocok > 0) THEN
            id%INFO(1)=-13
            CALL MUMPS_SET_IERROR(
     &             LWK_REAL+ 
     &             LWRKR_TH*int(max(NOMP_MAX,1),8)+
     &             LWRKC_TH*int(max(NOMP_MAX,1),8),
     &             id%INFO(2))
         ENDIF
C        --- Propagate enventual error
         CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF (id%INFO(1).LT.0) GOTO 517
C        -- estimation of memory and construction of partvecs
         BUJOB = 2
         CALL ZMUMPS_SIMSCALEABS(
     &        IRN_loc_PTR(1), JCN_loc_PTR(1), A_loc_PTR(1),
     &        id%KEEP8(29),
     &        M, N,  id%NPROCS, id%MYID, id%COMM,
     &        BURP, BUCP,
     &        BURS, BUCS, BUREGISTRE,
     &        IWK, LIWK,
     &        BUINTSZ, BURESZ, BUJOB,
     &        id%ROWSCA(1), id%COLSCA(1), WK_REAL, LWK_REAL,
     &        WRKR_TH, LWRKR_TH, 
     &        WRKC_TH, LWRKC_TH,
     &        NOMP_MAX, id%KEEP(50),
     &        K231, K232, K233, 
     &        id%DKEEP(3),
     &        SCONEERR, SCINFERR)
         id%DKEEP(4) = SCONEERR
         id%DKEEP(5) = SCINFERR
CXXXX 
         DEALLOCATE(IWK, BURP,BUCP,BURS, BUCS)
         DEALLOCATE(WK_REAL, WRKR_TH, WRKC_TH)
C}
        ELSE IF ( KEEP(54) .EQ. 0 ) THEN
C{        ------------------
C         Centralized matrix
C         ------------------
          IF ((KEEP(52).EQ.7).OR.(KEEP(52).EQ.8))  THEN
C             -------------------------------
C             Create a communicator of size 1
C             -------------------------------
              IF (id%MYID.EQ.MASTER) THEN
                COLOUR = 0
              ELSE
                COLOUR = MPI_UNDEFINED
              ENDIF
              CALL MPI_COMM_SPLIT( id%COMM, COLOUR, 0,
     &             COMM_FOR_SCALING, IERR )
              IF (id%MYID.EQ.MASTER) THEN
                 M = N
                 BUMAXMN=N
                 IF(N > BUMAXMN) BUMAXMN = N
                 LIWK = 1_8
                 ALLOCATE(IWK(LIWK), BURP(1), BUCP(1),
     &                BURS(1), BUCS(1),
     &                stat=allocok)
                 IF (allocok > 0) THEN
                    id%INFO(1)=-13
                    CALL MUMPS_SET_IERROR( LIWK+4_8, id%INFO(2) )
                    GOTO 400
                 ENDIF
                 LWK_REAL = int(M,8) + int(N,8)  
C
C                -- Set NOMP_MAX from KEEP(281)
                 CALL ZMUMPS_SET_NOMP_MAX(id%KEEP(281), id%KEEP(361), 
     &                                    N, NOMP_MAX)
C
                 IF (NOMP_MAX.LE.1) THEN
C                  temp array per thread not used
                   LWRKR_TH = 1
                   LWRKC_TH = 1
                 ELSE
                   LWRKR_TH = N
                   IF (id%KEEP(50).NE.0) THEN
C                    WRKC_TH not used on symmetric matrices
                     LWRKC_TH = 1
                   ELSE
                     LWRKC_TH = N
                   ENDIF
                 ENDIF
                 ALLOCATE(WK_REAL(LWK_REAL),
     &                    WRKR_TH(LWRKR_TH,max(NOMP_MAX,1)),
     &                    WRKC_TH(LWRKC_TH,max(NOMP_MAX,1)),
     &                    stat=allocok)
                 IF (allocok > 0) THEN
                    id%INFO(1)=-13
                    CALL MUMPS_SET_IERROR(
     &                     LWK_REAL+ 
     &                     LWRKR_TH*int(max(NOMP_MAX,1),8)+
     &                     LWRKC_TH*int(max(NOMP_MAX,1),8),
     &                     id%INFO(2))
                 ENDIF
                 CALL MPI_COMM_RANK(COMM_FOR_SCALING, SCMYID, IERR)
                 CALL MPI_COMM_SIZE(COMM_FOR_SCALING, SCNPROCS, IERR)
                 BUJOB = 1
                 CALL ZMUMPS_SIMSCALEABS(
     &                id%IRN(1), id%JCN(1), id%A(1),
     &                id%KEEP8(28),
     &                M, N,  SCNPROCS, SCMYID, COMM_FOR_SCALING,
     &                BURP, BUCP,
     &                BURS, BUCS, BUREGISTRE,
     &                IWK, LIWK,
     &                BUINTSZ, BURESZ, BUJOB,
     &                id%ROWSCA(1), id%COLSCA(1), WK_REAL, LWK_REAL,
     &                WRKR_TH, LWRKR_TH, 
     &                WRKC_TH, LWRKC_TH,
     &                NOMP_MAX, id%KEEP(50),
     &                K231, K232, K233, 
     &                id%DKEEP(3),
     &                SCONEERR, SCINFERR)
                 IF(LWK_REAL < BURESZ) THEN
                    id%INFO(1) = -136
                    GOTO 400
                 ENDIF
                 BUJOB = 2
                 CALL ZMUMPS_SIMSCALEABS(id%IRN(1),
     &                id%JCN(1), id%A(1),
     &                id%KEEP8(28),
     &                M, N,  SCNPROCS, SCMYID, COMM_FOR_SCALING,
     &                BURP, BUCP,
     &                BURS, BUCS, BUREGISTRE,
     &                IWK, LIWK,
     &                BUINTSZ, BURESZ, BUJOB,
     &                id%ROWSCA(1), id%COLSCA(1), WK_REAL, LWK_REAL,
     &                WRKR_TH, LWRKR_TH, 
     &                WRKC_TH, LWRKC_TH,
     &                NOMP_MAX, id%KEEP(50),
     &                K231, K232, K233, 
     &                id%DKEEP(3),
     &                SCONEERR, SCINFERR)
                 id%DKEEP(4) = SCONEERR
                 id%DKEEP(5) = SCINFERR
 400             CONTINUE
                 IF (allocated(WK_REAL)) DEALLOCATE(WK_REAL)
                 IF (allocated(WRKR_TH)) DEALLOCATE(WRKR_TH)
                 IF (allocated(WRKC_TH)) DEALLOCATE(WRKC_TH)
                 IF (allocated(IWK))     DEALLOCATE(IWK)
                 IF (allocated(BURP))    DEALLOCATE(BURP)
                 IF (allocated(BUCP))    DEALLOCATE(BUCP)
                 IF (allocated(BURS))    DEALLOCATE(BURS)
                 IF (allocated(BUCS))    DEALLOCATE(BUCS)
              ENDIF
C             Centralized matrix: make DKEEP(4:5) available to all processors
              CALL MPI_BCAST( id%DKEEP(4),2,MPI_DOUBLE_PRECISION,
     &                        MASTER, id%COMM, IERR )
              IF (id%MYID.EQ.MASTER) THEN
C               Communicator should only be
C               freed on the master process
                CALL MPI_COMM_FREE(COMM_FOR_SCALING, IERR)
              ENDIF
              CALL MUMPS_PROPINFO(ICNTL(1), id%INFO(1),
     &             id%COMM, id%MYID)
              IF (id%INFO(1).LT.0) GOTO 517
          ELSE IF (id%MYID.EQ.MASTER) THEN
C           -------------------
C           Centralized scaling
C           -------------------
            IF (KEEP(52).GT.0 .AND. KEEP(52).LE.6) THEN
C             ---------------------
C             Allocate temporary
C             workspace for scaling
C             ---------------------
              IF (KEEP(52) .eq. 1 ) THEN
C               No workspace indeed needed
                LWK = 1_8
                LWK_REAL = 1_8
              ELSE IF ( KEEP(52) .eq. 3 ) THEN 
                LWK = 1_8
                LWK_REAL = int(N,8)
              ELSE IF ( KEEP(52) .eq. 4 ) THEN
C               Options 3 or 4
                LWK = 1_8
                LWK_REAL = 2_8*int(N,8)
              END IF
C             Real workarray
              ALLOCATE( WK_REAL( LWK_REAL ), stat = IERR )
              IF ( IERR .GT. 0 ) THEN
                id%INFO(1) = -13
                CALL MUMPS_SET_IERROR( LWK_REAL, id%INFO(2) )
                GOTO 137
              END IF
C             Real/complex workarray
              ALLOCATE( WK( LWK ), stat = IERR )
              IF ( IERR .GT. 0 ) THEN
                id%INFO(1) = -13
                CALL MUMPS_SET_IERROR(LWK, id%INFO(2))
                GOTO 137
              END IF
              CALL ZMUMPS_FAC_A(N, id%KEEP8(28), KEEP(52), id%A(1),
     &             id%IRN(1), id%JCN(1),
     &             id%COLSCA(1), id%ROWSCA(1),
     &             WK, LWK, WK_REAL, LWK_REAL, ICNTL(1), id%INFO(1) )
              DEALLOCATE( WK_REAL )
              DEALLOCATE( WK )
            ENDIF
          ENDIF
C}
        ENDIF ! Scaling distributed matrices or centralized
        IF (KEEP(125).NE.0) THEN
C         ------------------------
C           If we enable the scaling of the |A11 A12| block
C           we set to 1 the scaling corresponding to the Schur
C           complement matrix A22
C         ------------------------
            IF ((KEEP(60).GT.0) .and. (KEEP(116).GT.0)) THEN
C            Schur is active, reset Schur entries to ONE
             IF ( ((KEEP(52).EQ.7).OR.(KEEP(52).EQ.8)) .AND.
     &          KEEP(54).NE.0 ) THEN
C             Scaling available on all procs
              DO I=1, N
               IF (id%SYM_PERM(I).GT.id%N-KEEP(116)) THEN
                 id%COLSCA(I) = ONE
                 id%ROWSCA(I) = ONE
               ENDIF
              ENDDO
             ELSE IF ( id%MYID .EQ. MASTER) THEN
C             Scaling available on master
              DO I=1, N
               IF (id%SYM_PERM(I).GT.id%N-KEEP(116)) THEN
                 id%COLSCA(I) = ONE
                 id%ROWSCA(I) = ONE
               ENDIF
              ENDDO
             ENDIF
            ENDIF
        ENDIF 
        IF (id%MYID.EQ.MASTER) THEN
            CALL MUMPS_SECFIN(TIMEET)
            id%DKEEP(92)=TIMEET
          IF (PROKG) WRITE( MPG, 140 ) TIMEET
C         Print inf-norm after last KEEP(233) iterations of
C         scaling option KEEP(52)=7 or 8 (SimScale)
C
          IF (PROKG.AND.(KEEP(52).EQ.7.OR.KEEP(52).EQ.8)) THEN
            IF (K233+K231+K232.GT.0) THEN
              IF (K232.GT.0) WRITE(MPG, 166) id%DKEEP(4)
            ENDIF
          ENDIF
        ENDIF
      ENDIF ! LSCAL
C
C       scaling might also be provided by the user
        LSCAL = (LSCAL .OR. (KEEP(52) .EQ. -1) .OR. KEEP(52) .EQ. -2)
        IF (LSCAL .AND. KEEP(258).NE.0 .AND. id%MYID .EQ. MASTER) THEN
          DO I = 1, id%N
            CALL ZMUMPS_UPDATEDETER_SCALING(id%ROWSCA(I),
     &           id%DKEEP(6),    ! determinant
     &           KEEP(259))   ! exponent of the determinant
          ENDDO
          IF (KEEP(50) .EQ. 0) THEN ! unsymmetric
            DO I = 1, id%N
              CALL ZMUMPS_UPDATEDETER_SCALING(id%COLSCA(I),
     &           id%DKEEP(6),    ! determinant
     &           KEEP(259))   ! exponent of the determinant
            ENDDO
          ELSE
C           -----------------------------------------
C           In this case COLSCA = ROWSCA
C           Since determinant was initialized to 1,
C           compute square of the current determinant
C           rather than going through COLSCA.
C           -----------------------------------------
            CALL ZMUMPS_DETER_SQUARE(id%DKEEP(6), KEEP(259))
          ENDIF
C         Now we should have taken the
C         inverse of the scaling vectors
          CALL ZMUMPS_DETER_SCALING_INVERSE(id%DKEEP(6), KEEP(259))
        ENDIF
C
C       ********************
C       End of Scaling phase
C       At this point: either (matrix is distributed and KEEP(52)=7 or 8)
C       in which case scaling arrays are allocated on all processors,
C       or scaling arrays are only on the host processor.
C       In case of distributed matrix input, we will free the scaling
C       arrays on procs with MYID .NE. 0 after the all-to-all distribution
C       of the original matrix.
C       ********************
C
 137  CONTINUE
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF (id%INFO(1).LT.0) GOTO 517
C     systematically this array now than waiting for
C     the root node. We rely on the fact that it is
C     allocated or not during the solve phase so if
C     it was allocated in a 1st call to facto and not
C     in a second, we don't want the solve to think
C     it was allocated in the second call.
      IF (associated(idintr%roota%RHS_CNTR_MASTER_ROOT)) THEN
        DEALLOCATE (idintr%roota%RHS_CNTR_MASTER_ROOT)
        NULLIFY (idintr%roota%RHS_CNTR_MASTER_ROOT)
      ENDIF
C     Fwd in facto: check that id%NRHS has not changed
      IF ( id%MYID.EQ.MASTER.AND. KEEP(252).EQ.1 .AND.
     &      id%NRHS .NE. id%KEEP(253) ) THEN
C         Error: NRHS should not have
C         changed since the analysis
          id%INFO(1)=-42
          id%INFO(2)=id%KEEP(253)
      ENDIF
      IF (id%KEEP(252) .EQ. 1) THEN
          IF ( id%MYID.NE.MASTER ) THEN
            id%KEEP(254) = N              ! Leading dimension
            id%KEEP8(85) = int(N,8)*int(id%KEEP(253),8) ! Tot size
            ALLOCATE(RHS_MUMPS(id%KEEP8(85)),stat=IERR)
            IF (IERR > 0) THEN
               id%INFO(1)=-13
               CALL MUMPS_SET_IERROR(id%KEEP8(85), id%INFO(2))
               IF (LPOK)
     &         WRITE(LP,*) 'ERROR while allocating RHS on a slave'
               NULLIFY(RHS_MUMPS)
              RHS_MUMPS_ALLOCATED = .FALSE.
            ELSE
              RHS_MUMPS_ALLOCATED = .TRUE.
            ENDIF
          ELSE 
C           Case of non working master
            id%KEEP(254)=id%LRHS              ! Leading dimension
            id%KEEP8(85)=int(id%LRHS,8)*int(id%KEEP(253)-1,8)+
     &                   int(id%N,8) ! Tot size
            RHS_MUMPS=>id%RHS
            RHS_MUMPS_ALLOCATED = .FALSE.
            IF (LSCAL) THEN
C             Scale before broadcast: apply row
C             scaling (remark that we assume no
C             transpose).
              DO K=1, id%KEEP(253)
                DO I=1, N
                  RHS_MUMPS( id%KEEP(254) * (K-1) + I )
     &          = RHS_MUMPS( id%KEEP(254) * (K-1) + I )
     &          * id%ROWSCA(I)
                ENDDO
              ENDDO
            ENDIF
          ENDIF
      ELSE
          id%KEEP8(85)=1_8
          ALLOCATE(RHS_MUMPS(1),stat=IERR)
          IF (IERR > 0) THEN
             id%INFO(1)=-13
             id%INFO(2)=1
             IF (LPOK)
     &            WRITE(LP,*) 'ERREUR while allocating RHS on a slave'
             NULLIFY(RHS_MUMPS)
             RHS_MUMPS_ALLOCATED = .FALSE.
          ELSE
             RHS_MUMPS_ALLOCATED = .TRUE.
#if defined(SINGLEINDOUBLE)
C            RHS_MUMPS will be accessed: we initialize it to avoid
C            a possible floating point exception when accessing.
             RHS_MUMPS(1)=1.0D0
#endif
          ENDIF
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO(1).lt.0 ) GOTO 516
      IF (KEEP(252) .EQ. 1) THEN
C
C         Broadcast the columns of the right-hand side
C         one by one. Leading dimension is keep(254)=N
C         on procs with MYID > 0 but may be larger on
C         the master processor.
          DO I= 1, id%KEEP(253)
            CALL MPI_BCAST(RHS_MUMPS((I-1)*id%KEEP(254)+1), N,
     &           MPI_DOUBLE_COMPLEX, MASTER,id%COMM,IERR)
          END DO
      ENDIF
      IF (id%MYID.EQ. MASTER) THEN
C       Copy the value of ICNTL(24) and make it
C       available on all working processors.
        KEEP(110)=id%ICNTL(24)
C       KEEP(110) defaults to 0 for out of range values
        IF (KEEP(110).NE.1) KEEP(110)=0
      ENDIF
      CALL MPI_BCAST(KEEP(110), 1, MPI_INTEGER,
     &               MASTER, id%COMM, IERR)
C     -----------------------------------------------
C     Depending on the option used for 
C       -detecting null pivots (ICNTL(24)/KEEP(110))
C         CNTL(3) is used to set DKEEP(1)
C               ( A row is considered as null if ||row|| < DKEEP(1) )
C         CNTL(5) is then used to define if a large 
C                 value is set on the diagonal or if a 1 is set
C                 and other values in the row are reset to zeros.
C       -rank revealing on the Schur (ICNTL(56)/KEEP(19))
C         SEUIL* corresponds to the minimum required 
C                absolute value of pivot.
C         SEUIL_LDLT_NIV2 is used only in the 
C                case of SYM=2 within a niv2 node for which 
C                we have only a partial view of the fully summed rows.
      IF (id%MYID .EQ. MASTER) CNTL3 = id%CNTL(3)
      CALL MPI_BCAST(CNTL3, 1, MPI_DOUBLE_PRECISION,
     &               MASTER, id%COMM, IERR)
      IF (id%MYID .EQ. MASTER) CNTL5 = id%CNTL(5)
      CALL MPI_BCAST(CNTL5, 1, MPI_DOUBLE_PRECISION,
     &               MASTER, id%COMM, IERR)
      IF (id%MYID .EQ. MASTER) CNTL6 = id%CNTL(6)
      CALL MPI_BCAST(CNTL6, 1, MPI_DOUBLE_PRECISION,
     &               MASTER, id%COMM, IERR)
      IF (id%MYID .EQ. MASTER) id%DKEEP(8) = id%CNTL(7)
      CALL MPI_BCAST(id%DKEEP(8), 1, MPI_DOUBLE_PRECISION,
     &               MASTER, id%COMM, IERR)
      id%DKEEP(11) = id%DKEEP(8)/id%KEEP(461)
      id%DKEEP(12) = id%DKEEP(8)/id%KEEP(462)
      IF (KEEP(486).EQ.0) id%DKEEP(8) = ZERO
      COMPUTE_ANORMINF = .FALSE.
      IF ( (KEEP(486) .NE. 0).AND. (id%DKEEP(8).LT.ZERO)) THEN
        COMPUTE_ANORMINF = .TRUE.
      ENDIF
      IF (KEEP(19).NE.0) THEN
C       Rank revealing factorisation
        COMPUTE_ANORMINF = .TRUE.
      ENDIF
      IF (KEEP(110).NE.0) THEN
C       Null pivot detection       
        COMPUTE_ANORMINF = .TRUE.
      ENDIF
      IF (id%DKEEP(8).LT.ZERO) THEN
C       Experimental setting of CNTL(7)
        IF (COMPUTE_ANORMINF) THEN
         EFF_SIZE_SCHUR = 0
         CALL ZMUMPS_ANORMINF(  id , ANORMINF, LSCAL, EFF_SIZE_SCHUR )
C        If no schur ANORMINF fine for other cases
         id%DKEEP(8) = abs(id%DKEEP(8))*ANORMINF
        ELSE
         ANORMINF = ZERO
         id%DKEEP(8) = abs(id%DKEEP(8))
        ENDIF
C       ANORMINF need be recomputed in case of schur
        IF ((KEEP(60).GT.0).AND.KEEP(116).GT.0) ANORMINF=ZERO
      ENDIF
      IF (PROKG) THEN
            IF ( (CNTL(7) < ZERO) .AND. COMPUTE_ANORMINF .AND.
     &          (KEEP(486) .NE. 0) ) THEN
C           Warning : using negative values is an experimental and 
C            non recommended setting.
             WRITE(MPG,'(/A,A/,A/,A,A)') 
     &  ' WARNING in BLR input setting: ',
     &  '          CNTL(7) < 0 is experimental: ',
     &  '          Effective BLR threshold = |CNTL(7| x ||A_pre||, ',
     &  '          where A_pre is the preprocessed matrix as defined',
     &  ' in the users guide '
             WRITE(MPG,'(A,3D16.4/)') 
     &  '          Effective BLR threshold, CNTL(7), ||A_pre|| = ', 
     &      id%DKEEP(8), CNTL(7), ANORMINF
            ENDIF
      ENDIF
      IF (KEEP(486).NE.0) THEN
#if defined(SINGLEINDOUBLE)
        IF (KEEP(148).EQ.0) THEN
#endif
          CALL ZMUMPS_INIT_LRMIXED_PRECISIONS(id%KEEP(148), 
     &          id%KEEP(451), id%KEEP(456), id%KEEP(457), id%DKEEP(8))
#if defined(SINGLEINDOUBLE)
        ELSE
C         id%DKEEP(8): casted to single precision explicitely
          LPDK8 = id%DKEEP(8)
          CALL CMUMPS_INIT_LRMIXED_PRECISIONS(id%KEEP(148),
     &          id%KEEP(451), id%KEEP(456), id%KEEP(457), LPDK8)
        ENDIF
#endif
       ENDIF
C     -------------------------------------------------------
C        We compute ANORMINF, when needed, based on
C        the infinite norm of Rowsca *A*Colsca
C        and make it available on all working processes.
      IF (COMPUTE_ANORMINF) THEN
         EFF_SIZE_SCHUR = 0
         IF (KEEP(60).GT.0) EFF_SIZE_SCHUR = KEEP(116)
         CALL ZMUMPS_ANORMINF(  id , ANORMINF, LSCAL, EFF_SIZE_SCHUR )
      ELSE
         ANORMINF = ZERO
      ENDIF
C
       IF ((KEEP(19).NE.0).OR.(KEEP(110).NE.0)) THEN 
         IF (PROKG) THEN 
           IF (KEEP(19).NE.0) THEN
           WRITE(MPG,'(A,1PD16.4)')
     &    ' CNTL(3) for null pivot rows/singularities  =',CNTL3      
           ELSE
           WRITE(MPG,'(A,1PD16.4)')
     &    ' CNTL(3) for null pivot row detection       =',CNTL3      
           ENDIF
         ENDIF
       ENDIF
       IF (KEEP(19).EQ.0) THEN 
C{       -- RR is off
         SEUIL = ZERO
         id%DKEEP(9) = ZERO
C}
       ELSE
C{       -- RR is on
C
C      CNTL(3) is the threshold used in the following to compute  
C      DKEEP(9) the threshold under which the sing val. are considered
C      as null and from which we start to look for a gap between two
C      sing val.
         IF (CNTL3 .LT. ZERO) THEN
           id%DKEEP(9) = abs(CNTL(3))
         ELSE IF  (CNTL3 .GT. ZERO) THEN
           id%DKEEP(9) = CNTL3*ANORMINF
         ELSE  !  (CNTL(3) .EQ. ZERO) THEN
           CALL MUMPS_NPIV_CRITICAL_PATH(
     &     N, KEEP(28), id%STEP(1), id%FRERE_STEPS(1), id%FILS(1),
     &     id%NA(1), id%LNA, id%NE_STEPS(1), KEEP(127), 
     &     NPIV_CRITICAL_PATH )
           id%DKEEP(9) = sqrt(dble(NPIV_CRITICAL_PATH))*EPS*ANORMINF  
         ENDIF
         IF (PROKG) THEN 
           WRITE(MPG, '(A,I16)')
     &    ' ICNTL(56) rank revealing effective value   =',KEEP(19)   
           WRITE(MPG,'(A,1PD16.4)')
     &    ' ...Threshold for singularities on the root =',id%DKEEP(9)
         ENDIF
C        RR postponing considers that pivot rows with norm smaller 
C        than SEUIL should be postponed. 
C        SEUIL should be bigger than DKEEP(9), this means that 
C        DKEEP(13) should be bigger than 1. 
         Thresh_Seuil = id%DKEEP(13)
         IF (id%DKEEP(13).LT.1)  Thresh_Seuil = 10
         SEUIL = id%DKEEP(9)*Thresh_Seuil
         IF (PROKG) WRITE(MPG,'(A,1PD16.4)')
     &   ' ...Threshold for postponing                =',SEUIL
C}
       ENDIF !end KEEP(19).ne.0
       SEUIL_LDLT_NIV2 = SEUIL
C     -------------------------------
C     -- Null pivot row detection
C     -------------------------------
       IF (KEEP(110).EQ.0) THEN 
C{       -- Null pivot is off
C        Initialize DKEEP(1) to a negative value
C        in order to avoid detection of null pivots
C        (test max(AMAX,RMAX,abs(PIVOT)).LE.PIVNUL
C        in ZMUMPS_FAC_I, where PIVNUL=DKEEP(1))
         id%DKEEP(1) = -1.0D0
         id%DKEEP(2) = ZERO
C}
       ELSE
C{       -- Null pivot detection is on
        IF (KEEP(19).NE.0) THEN
C{       -- RR is on
C      RR postponing considers that pivot rows of norm smaller that SEUIL 
C      should be postponed, but pivot rows smaller than DKEEP(1) are 
C      directly added to null space and thus considered as null pivot rows. 
         IF ((id%DKEEP(10).LE.0).OR.(id%DKEEP(10).GT.1)) THEN
C           DKEEP(10) is out of range, set to the default value 10-1
           id%DKEEP(1) = id%DKEEP(9)*1D-1
         ELSE
           id%DKEEP(1) = id%DKEEP(9)*id%DKEEP(10) 
         ENDIF    
C}
        ELSE
C{       -- RR is off
C        -- only Null pivot detection
C         We keep strategy currently used in MUMPS 4.10.0
          IF (CNTL3 .LT. ZERO) THEN
           id%DKEEP(1)  = abs(CNTL(3))
          ELSE IF  (CNTL3 .GT. ZERO) THEN
            id%DKEEP(1)  = CNTL3*ANORMINF
          ELSE !  (CNTL(3) .EQ. ZERO) THEN
c          id%DKEEP(1)  = NPIV_CRITICAL_PATH*EPS*ANORMINF
           CALL MUMPS_NPIV_CRITICAL_PATH(
     &         N, KEEP(28), id%STEP(1), id%FRERE_STEPS(1), id%FILS(1),
     &         id%NA(1), id%LNA, id%NE_STEPS(1), KEEP(127), 
     &         NPIV_CRITICAL_PATH )
           id%DKEEP(1)  = sqrt(dble(NPIV_CRITICAL_PATH))*EPS*ANORMINF
          ENDIF
C}
        ENDIF ! fin rank revealing        
        IF ((KEEP(110).NE.0).AND.(PROKG)) THEN
           WRITE(MPG, '(A,I16)')
     &    ' ICNTL(24) null pivot rows detection        =',KEEP(110)   
           WRITE(MPG,'(A,1PD16.4)')
     &    ' ...Zero pivot detection threshold          =',id%DKEEP(1)
        ENDIF 
        IF (CNTL5.GT.ZERO) THEN
           id%DKEEP(2) = CNTL5 * ANORMINF
           IF (PROKG) WRITE(MPG,'(A,1PD16.4)')
     &    ' ...Fixation for null pivots                =',id%DKEEP(2)
         ELSE
           IF (PROKG) WRITE(MPG,*) '...Infinite fixation '
           IF (id%KEEP(50).EQ.0) THEN
C            Unsym
           ! the user let us choose a fixation. set in NEGATIVE
           ! to detect during facto when to set row to zero !
#if defined(SINGLEINDOUBLE)
             IF (KEEP(148).EQ.0) THEN
#endif
               id%DKEEP(2) = -max(1.0D10*ANORMINF, 
     &                 sqrt(huge(ANORMINF))/1.0D8)
#if defined(SINGLEINDOUBLE)
             ELSE
C              Use single precision fixation since single
C              precision factorization will be used. Having
C              DKEEP(2) exceeding single precision max exponent
C              causes floating-point exceptions in conversion
C              before ZMUMPS_FAC_B or result in an invalid fixation.
               id%DKEEP(2) = -max(1.0D10*ANORMINF, 
     &                sqrt(huge(LPCNTL1))/1.0D8)
C              (we use LCNTL1 because it is declared single precision)
             ENDIF
#endif
           ELSE
C            Sym
             id%DKEEP(2) = ZERO
           ENDIF
         ENDIF
C}
       ENDIF ! fin null pivot detection.
C     Find id of root node if RR is on 
      IF (KEEP(19).NE.0) THEN
        ID_ROOT =MUMPS_PROCNODE(id%PROCNODE_STEPS(id%STEP(KEEP(20))),
     &                          id%KEEP(199))
        IF ( KEEP( 46 )  .NE. 1 ) THEN
          ID_ROOT = ID_ROOT + 1
        END IF
      ENDIF
      IF (KEEP(19).EQ.0.AND.KEEP(53).NE.0) THEN 
C       ICNTL(56)>0 at analysis and =0 at facto
C       save special root index
        KEEP20_SAVE = KEEP(20)
C       suppress special RR treatment
        KEEP(20)   = 0
      ENDIF
C Second pass:  set parameters for null pivot detection
C Allocate PIVNUL_LIST_STRUCT in case of null pivot detection
C and in case of rank revealing
      KEEP(109) = 0
      LPN_LIST  = 0
      IF(KEEP(110) .EQ. 1) THEN
         LPN_LIST = 100
      ENDIF
      IF (KEEP(19).NE.0 .AND.
     &   (ID_ROOT.EQ.id%MYID .OR. id%MYID.EQ.MASTER)) THEN
         LPN_LIST = 100
      ENDIF
      IF (LPN_LIST.GT.0) THEN
        PIVNUL_LIST_STRUCT%SIZE_PIVNUL_LIST = LPN_LIST
        ALLOCATE( PIVNUL_LIST_STRUCT%PIVNUL_LIST(LPN_LIST),stat = IERR )
        IF ( IERR .GT. 0 ) THEN
           id%INFO(1)=-13
           id%INFO(2)=LPN_LIST
        END IF
        PIVNUL_LIST_STRUCT%PIVNUL_LIST(1:LPN_LIST) = 0
      ENDIF
C end set parameter for null pivot detection
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO(1).lt.0 ) GOTO 516
C   --------------------------------------------------------------
C   STATIC PIVOTING
C     -- Static pivoting only when RR and Null pivot detection OFF
C   --------------------------------------------------------------
      KEEP(97) = 0
      IF ((KEEP(19).EQ.0).AND.(KEEP(110).EQ.0)) THEN
        IF (id%MYID .EQ. MASTER) CNTL4 = id%CNTL(4)
        CALL MPI_BCAST( CNTL4, 1, MPI_DOUBLE_PRECISION,
     &                MASTER, id%COMM, IERR )
C 
        IF ( CNTL4 .GE. ZERO ) THEN
         KEEP(97) = 1
         IF ( CNTL4 .EQ. ZERO ) THEN
C           -- set seuil to sqrt(eps)*||A||
            IF(ANORMINF .EQ. ZERO) THEN
               EFF_SIZE_SCHUR = 0
               IF (KEEP(60).GT.0) EFF_SIZE_SCHUR = KEEP(116)
               CALL ZMUMPS_ANORMINF( id , ANORMINF, LSCAL, 
     &              EFF_SIZE_SCHUR )
            ENDIF
            SEUIL = sqrt(EPS) * ANORMINF
         ELSE
            SEUIL = CNTL4
         ENDIF
         SEUIL_LDLT_NIV2 = SEUIL
        ELSE 
         SEUIL = ZERO
        ENDIF
      ENDIF
C     set number of tiny pivots / 2x2 pivots in types 1 /
C     2x2 pivots in types 2, to zero. This is because the
C     user can call the factorization step several times.
      KEEP(98)  = 0
      KEEP(103) = 0
      KEEP(105) = 0
      MAXS      = 1_8
*
*     Start allocations
*     *****************
*
C
C  The slaves can now perform the factorization
C
C
C  or point to user provided data WK_USER when LWK_USER>0
C  =======================
C
C     Compute BLR_STRAT and a first estimation 
C     of MAXS, the size of id%S
      CALL  ZMUMPS_SET_BLRSTRAT_AND_MAXS_K8 (
     &           MAXS_BASE8, MAXS_BASE_RELAXED8,
     &           BLR_STRAT,
     &           id%KEEP(1), id%KEEP8(1))
C   
      MAXS = MAXS_BASE_RELAXED8
      IF (WK_USER_PROVIDED) THEN
C       -- Set MAXS to size of WK_USER_
        MAXS = id%KEEP8(24)
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF (id%INFO(1) .LT. 0) THEN
        GOTO 516
      ENDIF
C
      id%KEEP8(75) = huge(id%KEEP8(75))
      id%KEEP8(76) = huge(id%KEEP8(76))
      IF (I_AM_SLAVE) THEN
C
        IF (id%KEEP8(4) .NE. 0_8) THEN
C
          IF ( .NOT. WK_USER_PROVIDED ) THEN
C           Set MAXS given BLR_STRAT, KEEP(201) and MAXS_BASE_RELAXED8
            CALL ZMUMPS_MEM_ALLOWED_SET_MAXS (
     &          MAXS,
     &          BLR_STRAT, id%KEEP(201), MAXS_BASE_RELAXED8,
     &          id%KEEP(1), id%KEEP8(1), id%MYID, id%N, id%NELT, 
     &          id%NA(1), id%LNA, id%NSLAVES, 
     &          KEEP464COPY,  KEEP465COPY,
     &          id%INFO(1), id%INFO(2)
     &           , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &           size(id%I8_L0_OMP,2)
     &           )
C           Given  MAXS and  max memory allowed KEEP8(4)
C           compute in KEEP8(75) the number of real/complex 
C           available for dynamic allocations 
            CALL  ZMUMPS_MEM_ALLOWED_SET_K75 ( 
     &           MAXS, id%MYID, 
     &           .FALSE., ! UNDER_L0_OMP
     &           N, id%NELT, id%NA(1), id%LNA, id%NSLAVES,
     &           BLR_STRAT, id%KEEP(201), 
     &           id%KEEP, id%KEEP8, id%INFO(1), id%INFO(2)
     &           , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &           size(id%I8_L0_OMP,2)
     &           )
          ELSE
C            KEEP8(75) dow not include MAXS, since WK_USER is provided
             CALL  ZMUMPS_MEM_ALLOWED_SET_K75 ( 
     &           0_8, id%MYID, 
     &           .FALSE., ! UNDER_L0_OMP
     &           N, id%NELT, id%NA(1), id%LNA, id%NSLAVES,
     &           BLR_STRAT, id%KEEP(201), 
     &           id%KEEP, id%KEEP8, id%INFO(1), id%INFO(2)
     &           , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &           size(id%I8_L0_OMP,2)
     &           )
          ENDIF
          IF (KEEP(400) .GT.0) THEN
C            ------------------------------
C            compute KEEP8(75) under L0_OMP
C            ------------------------------
C            Save KEEP8(75) above L0_OMP to reset KEEP8(75) 
C            when starting FAC_PAR_M
             id%KEEP8(76) = id%KEEP8(75)
             CALL  ZMUMPS_MEM_ALLOWED_SET_K75 ( 
     &             0_8, ! MAXS=0_8 
     &             id%MYID, 
     &             .TRUE., ! UNDER_L0_OMP
     &             id%N, id%NELT, id%NA(1), id%LNA, id%NSLAVES,
     &             BLR_STRAT, id%KEEP(201), 
     &             id%KEEP, id%KEEP8, id%INFO(1), id%INFO(2)
     &             , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &             size(id%I8_L0_OMP,2)
     &             )
C            KEEP8(75) holds the number of entries that 
C            can be allocated underL0.
C            It will be used during ZMUMPS_FAC_L0_OMP to adjust the 
C            the size of MUMPS_TPS_ARR(ITH)%LA 
          ENDIF
        ENDIF ! MEM_ALLOWED
C
      ENDIF ! I_AM_SLAVE THEN
C
      IF (I_AM_SLAVE) THEN
        IF ( (KEEP(400).GT.0) .AND. (KEEP(406).EQ.2) ) THEN
C         Compute KEEP8(77) the peak authorized used by 
C         ZMUMPS_PERFORM_COPIES
          CALL ZMUMPS_L0_COMPUTE_PEAK_ALLOWED(
     &          id%MYID, id%N,
     &          id%NELT, id%NA(1), id%LNA, id%NSLAVES,
     &          BLR_STRAT, id%KEEP(201),
     &          id%KEEP(1), id%KEEP8(1), id%INFO(1), id%INFO(2)
     &          , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &          size(id%I8_L0_OMP,2)
     &           )
        ENDIF
      ENDIF ! I_AM_SLAVE)
C
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF (id%INFO(1) .LT. 0) THEN
        GOTO 516
      ENDIF
      CALL MUMPS_SETI8TOI4(MAXS, id%INFO(39))
      CALL ZMUMPS_AVGMAX_STAT8(PROKG, MPG, MAXS, id%NSLAVES,
     & PRINT_MAXAVG,
     & id%COMM, " Effective size of S     (based on INFO(39))=   ")
C
      IF ( I_AM_SLAVE ) THEN
C       ----------------------------------------
C       Initialize some global variables related
C       to communication buffer management
C       ----------------------------------------
        CALL MUMPS_BUF_INI_MYID(id%MYID_NODES)
        CALL MUMPS_BUF_INIT( id%KEEP( 34 ), id%KEEP(35) )
C       ------------------
C       Dynamic scheduling
C       ------------------
        CALL MUMPS_LOAD_SET_INICOST( dble(id%COST_SUBTREES),
     &        KEEP(64), dble(id%DKEEP(15)), KEEP(375), MAXS )
        MEMORY_MD_ARG = min(int(PERLU,8) * ( MAXS_BASE8 / 100_8 + 1_8 ),
C       Restrict freedom from dynamic scheduler when
C       MEM_ALLOWED=ICNTL(23) is small (case where KEEP8(4)-MAXS_BASE8
C       is negative after call to ZMUMPS_MAX_MEM)
     &                      max(0_8, MAXS-MAXS_BASE8))
        CALL MUMPS_LOAD_INIT( MEMORY_MD_ARG, MAXS, id%KEEP,
     & id%KEEP8, id%INFO, id%ISTEP_TO_INIV2, id%CANDIDATES, id%ND_STEPS,
     & id%FILS, id%FRERE_STEPS, id%DAD_STEPS, id%PROCNODE_STEPS,
     & id%STEP, id%NE_STEPS, id%N, id%MAX_SURF_MASTER, id%SUP_PROC,
     & id%COMM_LOAD, id%COMM_NODES,
     & id%DEPTH_FIRST, id%COST_TRAV, id%DEPTH_FIRST_SEQ, id%SBTR_ID,
     & id%NA, id%NSLAVES, id%FUTURE_NIV2,
     & id%NBSA, id%NBSA_LOCAL, id%MEM_SUBTREE, id%MY_FIRST_LEAF,
     & id%MY_NB_LEAF, id%MY_ROOT_SBTR )
        IF (KEEP(201) .GT. 0) THEN
C          -------------------
C          OOC initializations
C          -------------------
           IF (KEEP(201).EQ.1 !PANEL Version
     &         .AND.KEEP(50).EQ.0 ! Unsymmetric
     &         .AND.KEEP(251).NE.2 ! Store L to disk
     &         ) THEN
              id%OOC_NB_FILE_TYPE=2 ! declared in MUMPS_OOC_COMMON
           ELSE
              id%OOC_NB_FILE_TYPE=1 ! declared in MUMPS_OOC_COMMON
           ENDIF
C          ------------------------------
C          Dimension IO buffer, KEEP(100)
C          ------------------------------
           IF (KEEP(205) .GT. 0) THEN
             KEEP(100) = KEEP(205)
           ELSE
             IF (KEEP(201).EQ.1) THEN ! PANEL version
               ITMP8 = int(id%OOC_NB_FILE_TYPE,8) *
     &               2_8 * int(KEEP(226),8)
             ELSE
               ITMP8 = 2_8 * id%KEEP8(119)
             ENDIF
             ITMP8 = ITMP8 +  int(max(KEEP(12),0),8) *
     &               (ITMP8/100_8+1_8)
C            we want to avoid too large IO buffers.
C            12M corresponds to 100Mbytes given to buffers.
             ITMP8 = min(ITMP8, 12000000_8)
             KEEP(100)=int(ITMP8)
           ENDIF
           IF (KEEP(201).EQ.1) THEN
C            Panel version. Force the use of a buffer. 
             IF ( KEEP(99) < 3 ) THEN
               KEEP(99) = KEEP(99) + 3
             ENDIF
           ENDIF
C          --------------------------
C          Reset KEEP(100) to 0 if no
C          buffer is used for OOC.
C          --------------------------
           IF (KEEP(99) .LT.3) KEEP(100)=0
           IF((dble(KEEP(100))*dble(KEEP(35))/dble(2)).GT.
     &          (dble(1999999999)))THEN
             IF (PROKG) THEN
               WRITE(MPG,*)id%MYID,': Warning: DIM_BUF_IO might be
     &  too big for Filesystem'
             ENDIF
           ENDIF
           ALLOCATE (id%OOC_INODE_SEQUENCE(KEEP(28),
     &          id%OOC_NB_FILE_TYPE),
     &          stat=IERR)
           IF ( IERR .GT. 0 ) THEN
              id%INFO(1) = -13
              id%INFO(2) = id%OOC_NB_FILE_TYPE*KEEP(28)
              NULLIFY(id%OOC_INODE_SEQUENCE)
              GOTO 112
           ENDIF
           ALLOCATE (id%OOC_TOTAL_NB_NODES(id%OOC_NB_FILE_TYPE),
     &          stat=IERR)
           IF ( IERR .GT. 0 ) THEN
              id%INFO(1) = -13
              id%INFO(2) = id%OOC_NB_FILE_TYPE
              NULLIFY(id%OOC_TOTAL_NB_NODES)
              GOTO 112
           ENDIF
           ALLOCATE (id%OOC_SIZE_OF_BLOCK(KEEP(28),
     &          id%OOC_NB_FILE_TYPE),
     &          stat=IERR)
           IF ( IERR .GT. 0 ) THEN
              id%INFO(1) = -13
              id%INFO(2) = id%OOC_NB_FILE_TYPE*KEEP(28)
              NULLIFY(id%OOC_SIZE_OF_BLOCK)
              GOTO 112
           ENDIF
           ALLOCATE (id%OOC_VADDR(KEEP(28),id%OOC_NB_FILE_TYPE),
     &          stat=IERR)
           IF ( IERR .GT. 0 ) THEN
              id%INFO(1) = -13
              id%INFO(2) = id%OOC_NB_FILE_TYPE*KEEP(28)
              NULLIFY(id%OOC_VADDR)
              GOTO 112
           ENDIF
        ENDIF
      ENDIF
 112  CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF (id%INFO(1) < 0) THEN
C       LOAD_END must be done but not OOC_END_FACTO
        GOTO 513
      ENDIF
      IF (I_AM_SLAVE) THEN
        IF (KEEP(201) .GT. 0) THEN
           IF ((KEEP(201).EQ.1).OR.(KEEP(201).EQ.2)) THEN
             CALL ZMUMPS_OOC_INIT_FACTO(id%ICNTL(1), id%ICNTL(4), id%N,
     &       id%NSLAVES, id%MYID, MAXS, id%OOC_NB_FILE_TYPE,
     &       id%KEEP, id%KEEP8, id%STEP, id%PROCNODE_STEPS,
     &       id%OOC_SIZE_OF_BLOCK, id%OOC_VADDR, id%INFO,
     &       id%OOC_TMPDIR, id%OOC_PREFIX, id%OOC_NB_FILES,
     &       id%OOC_INODE_SEQUENCE)
          ELSE
             WRITE(*,*) "Internal error in ZMUMPS_FAC_DRIVER"
             CALL MUMPS_ABORT()
           ENDIF
           IF(id%INFO(1).LT.0)THEN
              GOTO 111
           ENDIF
        ENDIF
C       First increment corresponds to the number of
C       floating-point operations for subtrees allocated
C       to the local processor.
        CALL MUMPS_LOAD_UPDATE(0,.FALSE.,dble(id%COST_SUBTREES),
     &          id%KEEP(1),id%KEEP8(1))
        IF (id%INFO(1).LT.0) GOTO 111
      END IF
C     -----------------------
C     Manage main workarray S
C     -----------------------
      EARLYT3ROOTINS = KEEP(200) .EQ.0
     &        .OR.  ( KEEP(200) .LT. 0 .AND. KEEP(400) .EQ. 0 )
#if defined(SINGLEINDOUBLE)
C     We cannot easily manage EARLYT3ROOTINS with single precision
C     factorization since in this case, S should be single precision
C     and allocated during ZMUMPS_FAC_B
      EARLYT3ROOTINS = .FALSE.
C     Note that this is done independently of KEEP(148), because
C     KEEP(148) is known at facto and decision should be
C     identical to the one from analysis.
C     Also, this is done both in single precision and double precision
C     routines, to keep the value coherent when mixing precisions.
#endif
#if defined (LARGEMATRICES)
      IF ( id%MYID .ne. MASTER ) THEN
#endif
      IF (.NOT.WK_USER_PROVIDED) THEN
        IF ( EARLYT3ROOTINS ) THEN
C         Standard allocation strategy
          CALL ZMUMPS_DM_ALLOC_S_WK(id%S, MAXS, IERR,
     &    KEEP(430), KEEP(35))
          id%KEEP8(23) = MAXS
          IF ( IERR .GT. 0 ) THEN
            id%INFO(1) = -13
            CALL MUMPS_SET_IERROR(MAXS, id%INFO(2))
C           On some platforms (IBM for example), an
C           allocation failure returns a non-null pointer.
C           Therefore we nullify S
            NULLIFY(id%S)
            id%KEEP8(23)=0_8
          ENDIF
        ENDIF
      ELSE
#if defined(SINGLEINDOUBLE)
#endif
        id%S => id%WK_USER(1:id%KEEP8(24))
        id%KEEP8(23) = 0_8
#if defined(SINGLEINDOUBLE)
#endif
      ENDIF
#if defined (LARGEMATRICES)
      END IF
#endif
C
 111  CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO(1).LT.0 ) GOTO 514
C     --------------------------
C     Initialization of modules
C     related to data management
C     --------------------------
      NB_ACTIVE_FRONTS_ESTIM = 3
      NB_THREADS = 1
!$    NB_THREADS = OMP_GET_MAX_THREADS()
C     
      NB_ACTIVE_FRONTS_ESTIM = 3*NB_THREADS
      IF (I_AM_SLAVE) THEN
C
        CALL MUMPS_FDM_INIT('A',NB_ACTIVE_FRONTS_ESTIM, id%INFO)
C
        IF ( (KEEP(486).EQ.2) 
     &     .OR. ((KEEP(489).NE.0).AND.(KEEP(400).GT.1))
     &    ) THEN
C         In case of LRSOLVE or CompressCB, 
C         initialize nb of handlers to nb of BLR
C         nodes estimated at analysis
          NB_FRONTS_F_ESTIM = KEEP(470)
        ELSE
          IF (KEEP(489).NE.0) THEN
C          Compress CB and no L0 OMP (or 1 thread under L0): 
C          NB_ACTIVE_FRONTS_ESTIM is too small, 
C          to limit nb of reallocations make it twice larger
           NB_FRONTS_F_ESTIM = 2*NB_ACTIVE_FRONTS_ESTIM
          ELSE
           NB_FRONTS_F_ESTIM = NB_ACTIVE_FRONTS_ESTIM
          ENDIF
        ENDIF
        CALL MUMPS_FDM_INIT('F',NB_FRONTS_F_ESTIM, id%INFO )
        IF (id%INFO(1) .LT. 0 ) GOTO 114
#if ! defined(NO_FDM_DESCBAND)
C         Storage of DESCBAND information
          CALL MUMPS_FDBD_INIT( NB_ACTIVE_FRONTS_ESTIM, id%INFO )
#endif
#if ! defined(NO_FDM_MAPROW)
C         Storage of MAPROW and ROOT2SON information
          CALL MUMPS_FMRD_INIT( NB_ACTIVE_FRONTS_ESTIM, id%INFO )
#endif
#if defined(SINGLEINDOUBLE)
        IF (KEEP(148).EQ.0) THEN
#endif
          CALL ZMUMPS_BLR_INIT_MODULE( NB_FRONTS_F_ESTIM, id%INFO 
     &     )
#if defined(SINGLEINDOUBLE)
        ELSE IF (KEEP(148).EQ.2) THEN
C           call ZMUMPS_BLR_INIT_MODULE is postponed to 
C           ZMUMPS_LPid_TO_Cid with final size of BLR_ARRAY
          CALL CMUMPS_BLR_INIT_MODULE( NB_FRONTS_F_ESTIM, id%INFO 
     &     )
        ELSE
C         call single precision version of BLR module
          CALL CMUMPS_BLR_INIT_MODULE( NB_FRONTS_F_ESTIM, id%INFO
     &    )
        ENDIF
#endif
 114    CONTINUE
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                       id%COMM, id%MYID )
C     GOTO 500: one of the above module initializations failed
      IF ( id%INFO(1).LT.0 ) GOTO 500
C
C
C  Allocate space for matrix in arrowhead form
C  ===========================================
C
C  CASE 1 : Matrix is assembled
C  CASE 2 : Matrix is elemental
C
      IF ( KEEP(55) .eq. 0 ) THEN
C       ------------------------------------
C       Space has been allocated already for
C       the integer part during analysis
C       Only slaves need the arrowheads.
C       ------------------------------------
        IF ( I_AM_SLAVE .and. id%KEEP8(26) .ne. 0_8 ) THEN
          ALLOCATE( DBLARR( id%KEEP8(26) ),
     &              INTARR( id%KEEP8(27) ), stat = IERR )
        ELSE
          ALLOCATE( DBLARR( 1 ),
     &              INTARR( 1 ), stat =IERR )
        END IF
        IF ( IERR .NE. 0 ) THEN
          IF (LPOK) THEN
            WRITE(LP,*) id%MYID,
     &      ': Allocation error for DBLARR(',id%KEEP8(26),')+INTARR(',
     &      id%KEEP8(27),')'
          ENDIF
          id%INFO(1)=-13
          CALL MUMPS_SET_IERROR( max(id%KEEP8(26),1_8)+
     &                           max(id%KEEP8(27),1_8),
     &                            id%INFO(2) )
          GOTO 100
        END IF
        DBLARR_ALLOCATED = .TRUE.
      ELSE
C        --------------------------------
C        Allocate element variables lists
C        --------------------------------
         IF ( I_AM_SLAVE .and. id%KEEP8(27) .ne. 0_8 ) THEN
           ALLOCATE( INTARR( id%KEEP8(27) ), stat = allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             CALL MUMPS_SET_IERROR(id%KEEP8(27), id%INFO(2))
             GOTO 100
           END IF
         ELSE
C          INTARR also allocated of size 1 on non-working master
           ALLOCATE( INTARR(1),stat=allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             id%INFO(2) = 1
             GOTO 100
           END IF
         END IF
C        ---------------------------------------
C        Allocate DBLARR to hold possibly scaled
C        copies of elemental matrices
C        On a working master (hybrid host) and
C        no scaling, avoid the copy and point
C        directly to user data instead.
C        ---------------------------------------
         IF ( I_AM_SLAVE ) THEN
           IF (      id%MYID_NODES .eq. MASTER
     &       .AND.   KEEP(46)   .eq. 1
     &       .AND.   KEEP(52)   .eq. 0 ) THEN
C            -------------------
C            Pointer association
C            -------------------
             DBLARR => id%A_ELT
           ELSE
C            ----------
C            Allocation
C            ----------
             ALLOCATE( DBLARR( max(id%KEEP8(26),1_8) ), stat=allocok )
             IF ( allocok .GT. 0 ) THEN
               id%INFO(1) = -13
               CALL MUMPS_SET_IERROR(max(id%KEEP8(26),1_8), id%INFO(2))
               NULLIFY(DBLARR)
               GOTO 100
             END IF
             DBLARR_ALLOCATED = .TRUE.
           END IF
         ELSE
           ALLOCATE( DBLARR(1), stat = allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             id%INFO(2) = 1
             NULLIFY(DBLARR)
             GOTO 100
           END IF
#if defined(SINGLEINDOUBLE)
C          In case it needs to be converted
           DBLARR(1)=CZERO
#endif
           DBLARR_ALLOCATED = .TRUE.
         END IF
      END IF
 100  CONTINUE
C     ----------------
C     Check for errors
C     ----------------
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                     id%COMM, id%MYID )
      IF ( id%INFO(1).LT.0 ) GOTO 500
C     ------------------------------------------
C     Prepare stuff for the 2D block-cyclic root
C     ------------------------------------------
      IF ( KEEP(38).NE.0 ) THEN
          ALLOCATE(idintr%root%RG2L(id%N), stat=allocok)
          IF ( allocok .NE. 0 ) THEN
            id%INFO(1)=-13
            id%INFO(2)=id%N
          END IF
        IF ( id%INFO(1) .GE.0 ) THEN
          CALL ZMUMPS_INIT_ROOT_FAC( id%N, id%MYID,
     &    idintr%root, id%FILS(1), id%KEEP(1) )
        ENDIF
      ENDIF
C     ----------------
C     Check for errors
C     ----------------
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                     id%COMM, id%MYID )
      IF ( id%INFO(1).LT.0 ) GOTO 500
C
C       -----------------------------------
C
C       DISTRIBUTION OF THE ORIGINAL MATRIX
C
C       -----------------------------------
C
C     TIMINGS: computed (and printed) on the host
C     Next line: global time for distrib(arrowheads,elts)
C     on the host. Synchronization has been performed.
      IF (id%MYID.EQ.MASTER) CALL MUMPS_SECDEB(TIME)
C     -------------------------------------------
C     S_PTR_ARG / MAXS_ARG will be used for id%S
C     argument to arrowhead/element distribution
C     routines: if id%S is not allocated, we pass
C     S_DUMMY_ARG instead, which is not accessed.
C     -------------------------------------------
      IF (EARLYT3ROOTINS) THEN
            S_PTR_ARG => id%S
            MAXS_ARG = MAXS
      ELSE
            S_PTR_ARG => S_DUMMY_ARG
            MAXS_ARG = 1
      ENDIF
C
      IF ( KEEP( 55 ) .eq. 0 ) THEN
C     {
C       ----------------------------
C       Original matrix is assembled
C       Arrowhead format to be used.
C       ----------------------------
C       KEEP8(26) and KEEP8(27) hold the number of entries for real/integer
C       for the matrix in arrowhead format. They have been set by the
C       analysis phase (ZMUMPS_ANA_F and ZMUMPS_ANA_G)
C
C       ------------------------------------------------------------------
C       Blocking is used for sending arrowhead records (I,J,VAL)
C              buffer(1) is used to store number of bytes already packed
C              buffer(2) number of records already packed
C       KEEP(39) : Number of records (blocking factor)
C       ------------------------------------------------------------------
C
C     ---------------------------------------------
C     In case of parallel root compute minimum
C     size of workspace to receive arrowheads
C     of root node. Will be used to check that
C     MAXS is large enough for arrowheads (case
C     of EARLYT3ROOTINS (KEEP(200)=0); if .NOT.
C     EARLYT3ROOTINS (KEEP(200)=1), root will
C     be assembled into id%S later and size of
C     id%S will be checked later)
C     ---------------------------------------------
      IF (EARLYT3ROOTINS .AND. KEEP(38).NE.0 .AND.
     &    KEEP(60) .EQ.0 .AND. I_AM_SLAVE) THEN
        LWK = int(MUMPS_NUMROC( idintr%root%ROOT_SIZE,
     &             idintr%root%MBLOCK,
     &             idintr%root%MYROW, 0, idintr%root%NPROW ),8)
        LWK = max( 1_8, LWK )
        LWK = LWK*
     &        int(MUMPS_NUMROC( idintr%root%ROOT_SIZE,
     &        idintr%root%NBLOCK,
     &        idintr%root%MYCOL, 0, idintr%root%NPCOL ),8)
        LWK = max( 1_8, LWK )
      ELSE
        LWK = 1_8
      ENDIF
C     MAXS must be at least 1, and in case of
C     parallel root, large enough to receive
C     arrowheads of root.
      IF (MAXS .LT. int(LWK,8)) THEN
           id%INFO(1) = -9
           CALL MUMPS_SET_IERROR(LWK, id%INFO(2))
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
      IF ( id%INFO(1).LT.0 ) GOTO 500
C
C     =====================================================
      IF (KEEP(399).GE.2) THEN
C{    Multihtreaded algorithm taking into account all cases
C     =====================================================
C
C        Timing on master.
         IF (id%MYID.EQ.MASTER) THEN
           CALL MUMPS_SECDEB(TIME)
         END IF
C
         NZ_loc8   = 0_8
         NBRECORDS = KEEP(39)
         SIZESCAL = id%N
C        Set  NZ_loc8, A_loc_PTR, IRN_loc_PTR, JCN_loc_PTR
C        and update NBRECORDS
         IF (KEEP(54).EQ.0) THEN
C         centralized matrix
          IF (id%KEEP8(28) .LT. int(NBRECORDS,8)) THEN
              NBRECORDS = int(id%KEEP8(28))
          ENDIF
          IF (id%MYID.EQ.MASTER) THEN
           NZ_loc8 = id%KEEP8(28)
           A_loc_PTR   => id%A
           IRN_loc_PTR => id%IRN
           JCN_loc_PTR => id%JCN
           IF (LSCAL) THEN
             COLSCA_PTR => id%COLSCA
             ROWSCA_PTR => id%ROWSCA
           ELSE
             COLSCA_PTR  => DUMMYSCA
             ROWSCA_PTR  => DUMMYSCA
             SIZESCAL  = 1
           ENDIF
          ELSE
           A_loc_PTR   => DUMMYA_loc
           IRN_loc_PTR => DUMMYIRN_loc
           JCN_loc_PTR => DUMMYJCN_loc
           COLSCA_PTR  => DUMMYSCA
           ROWSCA_PTR  => DUMMYSCA
           SIZESCAL  = 1
          ENDIF
         ELSE
C         distributed matrix
C         ---------------------------------------------------
C         In order to have possibly IRN_loc/JCN_loc/A_loc
C         of size 0, avoid to pass them inside REDISTRIBUTION
C         and pass id instead
C         NZ_locMAX8 gives as a maximum buffer size (send/recv) used
C          an upper bound to limit buffers on small matrices
C         ---------------------------------------------------
          CALL MPI_ALLREDUCE(id%KEEP8(29), NZ_locMAX8, 1, 
     &              MPI_INTEGER8, MPI_MAX, id%COMM, IERR)
          IF (NZ_locMAX8 .LT. int(NBRECORDS,8)) THEN
              NBRECORDS = int(NZ_locMAX8)
          ENDIF
          NZ_loc8 = id%KEEP8(29)
          LSCAL   = (KEEP(52).EQ.7).OR.(KEEP(52).EQ.8)
C         available on all MPI
          IF (LSCAL) THEN
            COLSCA_PTR => id%COLSCA
            ROWSCA_PTR => id%ROWSCA
          ELSE
           COLSCA_PTR  => DUMMYSCA
           ROWSCA_PTR  => DUMMYSCA
           SIZESCAL  = 1
          ENDIF
         ENDIF
#if ! defined(PCPRET)
          IF (id%KEEP(72).EQ.1) THEN
           NBRECORDS = max(3,NBRECORDS/10)
          ENDIF
         CALL ZMUMPS_FAC_DIST_ARROWHEADS_OMP ( id%N,
     &     NZ_loc8, 
C          replace id by:
     &        A_loc_PTR(1), IRN_loc_PTR(1), JCN_loc_PTR(1),
     &        SIZESCAL, LSCAL, COLSCA_PTR(1), ROWSCA_PTR(1), 
     &     DBLARR(1), id%KEEP8(26), INTARR(1),
     &     id%KEEP8(27),
     &     id%PTR8ARR(1), id%NINCOLARR(1),
     &     id%NINROWARR(1), id%PTRDEBARR(1),
     &     id%FILS(1),
     &     KEEP(1), id%KEEP8(1), id%MYID, id%COMM, NBRECORDS,
     &     S_PTR_ARG(1), MAXS_ARG, idintr%root, idintr%roota,
     &     id%PROCNODE_STEPS(1),
     &     id%NPROCS, id%NSLAVES, id%SYM_PERM(1), id%STEP(1),
     &     id%ICNTL(1), id%INFO(1), NSEND8, NLOCAL8,
     &     id%ISTEP_TO_INIV2(1),
     &     id%CANDIDATES(1,1) 
     &     )
         CALL MPI_BARRIER(id%COMM, IERR)
#else 
#endif
C     
         IF ( KEEP(46) .eq. 0 .AND. id%MYID.eq.MASTER ) THEN
C          ------------------------------
C          The host is not working -> had
C          no data from initial matrix
C          ------------------------------
           NSEND8  = 0_8
           NLOCAL8 = 0_8
         END IF
C        -------------------------
C        Check for possible errors
C        -------------------------
         CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
         IF ( id%INFO( 1 ) .LT. 0 ) GOTO 500
         CALL MPI_REDUCE( NSEND8, NSEND_TOT8, 1, MPI_INTEGER8,
     &                 MPI_SUM, MASTER, id%COMM, IERR )
         CALL MPI_REDUCE( NLOCAL8, NLOCAL_TOT8, 1, MPI_INTEGER8,
     &                 MPI_SUM, MASTER, id%COMM, IERR )
         IF ( PROKG ) THEN
           WRITE(MPG,125) NLOCAL_TOT8, NSEND_TOT8
         END IF
#if defined(LARGEMATRICES)
        write(6,*) " -DLARGEMATRICES: deallocate input matrix !!!"
        IF (associated(IRN_loc_PTR)) THEN
           DEALLOCATE(IRN_loc_PTR)
           NULLIFY(IRN_loc_PTR)
        ENDIF
        IF (associated(JCN_loc_PTR)) THEN
           DEALLOCATE(JCN_loc_PTR)
           NULLIFY(JCN_loc_PTR)
        ENDIF
        IF (associated(A_loc_PTR)) THEN
           DEALLOCATE(A_loc_PTR)
           NULLIFY(A_loc_PTR)
        ENDIF
        IF (KEEP(54).EQ.0) THEN
          IF (id%MYID.EQ.MASTER) THEN
            NULLIFY(id%A)
            NULLIFY(id%IRN)
            NULLIFY(id%JCN)
          ENDIF
        ELSE
          NULLIFY(id%A_LOC)
          NULLIFY(id%IRN_LOC)
          NULLIFY(id%JCN_LOC)
        ENDIF
#endif
C}
      ELSE
C{   =======================================================
      IF ( KEEP(54) .eq. 0 ) THEN
C     {
C       ================================================
C       FIRST CASE : MATRIX IS NOT INITIALLY DISTRIBUTED
C       ================================================
C       A small integer workspace is needed to
C       send the arrowheads.
        IF ( id%MYID .eq. MASTER ) THEN
#if defined(LARGEMATRICES)
          ALLOCATE (WK(LWK),stat=IERR)
          IF ( IERR .GT. 0 ) THEN
            id%INFO(1) = -13
            CALL MUMPS_SET_IERROR(LWK, id%INFO(2))
            write(6,*) ' PB1 ALLOC LARGEMAT'
          ENDIF
#endif
        ENDIF
        CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
        IF ( id%INFO(1).LT.0 ) GOTO 500
        IF ( id%MYID .eq. MASTER ) THEN
C
C         --------------------------------
C         MASTER sends arowheads using the
C         global communicator with ranks
C         also in global communicator
C         IWK is used as temporary
C         workspace of size N.
C         --------------------------------
          NBRECORDS = KEEP(39)
          IF (id%KEEP8(28) .LT. int(NBRECORDS,8)) THEN
            NBRECORDS = int(id%KEEP8(28))
          ENDIF
#if defined(LARGEMATRICES)
#if defined(SINGLEINDOUBLE)
          IF (KEEP(148).NE.0) THEN
         WRITE(*,*) "SINGLEINDOUBLE not compatible with -DLARGEMATRICES"
            CALL MUMPS_ABORT()
          ENDIF
#endif
          CALL ZMUMPS_FACTO_SEND_ARROWHEADS(id%N, id%KEEP8(28), id%A(1),
     &      id%IRN(1), id%JCN(1), id%SYM_PERM(1),
     &      LSCAL, id%COLSCA(1), id%ROWSCA(1),
     &      id%MYID, id%NSLAVES, id%PROCNODE_STEPS(1),
     &      NBRECORDS,
     &      id%COMM, idintr%root, idintr%roota, KEEP,id%KEEP8,
     &      id%FILS(1),
     &
     &      INTARR(1), id%KEEP8(27), DBLARR(1), id%KEEP8(26),
     &      id%PTR8ARR(1), id%NINCOLARR(1),
     &      id%NINROWARR(1), id%PTRDEBARR(1),
     &      id%FRERE_STEPS(1), id%STEP(1), WK(1), LWK,
     &      id%ISTEP_TO_INIV2(1), id%I_AM_CAND(1),
     &      id%CANDIDATES(1,1), id%ICNTL(1), id%INFO(1) )
C         write(6,*) '!!! A,IRN,JCN are freed during factorization '
          DEALLOCATE (id%A)
          NULLIFY(id%A)
          DEALLOCATE (id%IRN)
          NULLIFY (id%IRN)
          DEALLOCATE (id%JCN)
          NULLIFY (id%JCN)
          IF (.NOT.WK_USER_PROVIDED) THEN
            IF (EARLYT3ROOTINS) THEN
C             Standard allocation strategy
              CALL ZMUMPS_DM_ALLOC_S_WK(id%S, MAXS, IERR,
     &              KEEP(430), KEEP(35))
              id%KEEP8(23) = MAXS
              IF ( IERR .GT. 0 ) THEN
                id%INFO(1) = -13
                id%INFO(2) = MAXS
                NULLIFY(id%S)
                id%KEEP8(23)=0_8
                write(6,*) ' PB2 ALLOC LARGEMAT',MAXS
                CALL MUMPS_ABORT()
              ENDIF
            ENDIF
          ELSE
            id%S => id%WK_USER(1:id%KEEP8(24))
          ENDIF
          IF (EARLYT3ROOTINS) THEN
            id%S(MAXS-LWK+1_8:MAXS) = WK(1_8:LWK)
          ENDIF
          DEALLOCATE (WK)
#else
          CALL ZMUMPS_FACTO_SEND_ARROWHEADS(id%N, id%KEEP8(28),
     &    id%A(1), id%IRN(1), id%JCN(1), id%SYM_PERM(1),
     &    LSCAL, id%COLSCA(1), id%ROWSCA(1),   
     &    id%MYID, id%NSLAVES, id%PROCNODE_STEPS(1),
     &    NBRECORDS,
     &    id%COMM, idintr%root, idintr%roota, KEEP(1),id%KEEP8(1),
     &    id%FILS(1),
     &
     &    INTARR(1), id%KEEP8(27), DBLARR(1), id%KEEP8(26),
     &    id%PTR8ARR(1), id%NINCOLARR(1),
     &    id%NINROWARR(1), id%PTRDEBARR(1),
     &    id%FRERE_STEPS(1), id%STEP(1), S_PTR_ARG(1), MAXS_ARG,
     &    id%ISTEP_TO_INIV2(1), id%I_AM_CAND(1),
     &    id%CANDIDATES(1,1), id%ICNTL(1), id%INFO(1) )
#endif
        ELSE
          NBRECORDS = KEEP(39)
          IF (id%KEEP8(28) .LT. int(NBRECORDS,8)) THEN
            NBRECORDS = int(id%KEEP8(28))
          ENDIF
          CALL ZMUMPS_FACTO_RECV_ARROWHD2( id%N,
     &       DBLARR(1), id%KEEP8(26),
     &       INTARR(1), id%KEEP8(27),
     &       id%PTR8ARR(1), id%NINCOLARR(1),
     &       id%NINROWARR(1), id%PTRDEBARR(1),
     &       KEEP( 1 ), id%KEEP8(1), id%FILS(1), id%MYID, id%COMM,
     &       NBRECORDS,
     &
     &       S_PTR_ARG(1), MAXS_ARG,
     &       idintr%root, idintr%roota,
     &       id%PROCNODE_STEPS(1), id%NSLAVES,
     &       id%SYM_PERM(1), id%FRERE_STEPS(1), id%STEP(1),
     &       id%ICNTL(1), id%INFO(1) )
        ENDIF
C     }
      ELSE
C     {
C     =============================================
C     SECOND CASE : MATRIX IS INITIALLY DISTRIBUTED
C     =============================================
C     Timing on master.
      IF (id%MYID.EQ.MASTER) THEN
        CALL MUMPS_SECDEB(TIME)
      END IF
      IF ( I_AM_SLAVE ) THEN
C     {
C       ---------------------------------------------------
C       In order to have possibly IRN_loc/JCN_loc/A_loc
C       of size 0, avoid to pass them inside REDISTRIBUTION
C       and pass id instead
C       NZ_locMAX8 gives as a maximum buffer size (send/recv) used
C        an upper bound to limit buffers on small matrices
C       ---------------------------------------------------
        CALL MPI_ALLREDUCE(id%KEEP8(29), NZ_locMAX8, 1, MPI_INTEGER8,
     &                     MPI_MAX, id%COMM_NODES, IERR)
        NBRECORDS = KEEP(39)
        IF (NZ_locMAX8 .LT. int(NBRECORDS,8)) THEN
            NBRECORDS = int(NZ_locMAX8)
        ENDIF
        CALL ZMUMPS_REDISTRIBUTION( id%N,
     &  id%KEEP8(29),
     &  id,
     &  DBLARR(1), id%KEEP8(26), INTARR(1),
     &  id%KEEP8(27),
     &  id%PTR8ARR(1), id%NINCOLARR(1),
     &  id%NINROWARR(1), id%PTRDEBARR(1),
     &  id%FILS(1),
     &  KEEP(1), id%KEEP8(1), id%MYID_NODES,
     &  id%COMM_NODES, NBRECORDS,
     &  S_PTR_ARG(1), MAXS_ARG, idintr%root, idintr%roota,
     &  id%PROCNODE_STEPS(1),
     &  id%NSLAVES, id%SYM_PERM(1), id%STEP(1),
     &  id%ICNTL(1), id%INFO(1), NSEND8, NLOCAL8,
     &  id%ISTEP_TO_INIV2(1),
     &  id%CANDIDATES(1,1) )
        IF ( ( KEEP(52).EQ.7 ).OR. (KEEP(52).EQ.8) ) THEN
C         -------------------------------------------------
C         In that case, scaling arrays have been allocated
C         on all processors. They were useful for matrix
C         distribution. But we now really only need them
C         on the host. In case of distributed solution, we
C         will have to broadcast either ROWSCA or COLSCA
C         (depending on MTYPE) but this is done later.
C
C         In other words, on exit from the factorization,
C         we want to have scaling arrays available only
C         on the host.
C         -------------------------------------------------
          IF ( id%MYID > 0 ) THEN
            IF (associated(id%ROWSCA)) THEN
              DEALLOCATE(id%ROWSCA)
              NULLIFY(id%ROWSCA)
            ENDIF
            IF (associated(id%COLSCA)) THEN
              DEALLOCATE(id%COLSCA)
              NULLIFY(id%COLSCA)
            ENDIF
          ENDIF
        ENDIF
#if defined(LARGEMATRICES)
C      deallocate id%IRN_loc, id%JCN(loc) to free extra space
C      Note that in this case IRN_loc cannot be used
C      anymore during the solve phase for IR and Error analysis.
        IF (associated(id%IRN_loc)) THEN
           DEALLOCATE(id%IRN_loc)
           NULLIFY(id%IRN_loc)
        ENDIF
        IF (associated(id%JCN_loc)) THEN
           DEALLOCATE(id%JCN_loc)
           NULLIFY(id%JCN_loc)
        ENDIF
        IF (associated(id%A_loc)) THEN
           DEALLOCATE(id%A_loc)
           NULLIFY(id%A_loc)
        ENDIF
        write(6,*) ' Warning :', 
     &        ' id%A_loc, IRN_loc, JCN_loc deallocated !!! '
#endif
        IF (PROK) THEN
          WRITE(MP,120) NLOCAL8, NSEND8
        END IF
C     }
      END IF ! I_AM_SLAVE
C     -------------------------
C     Check for possible errors
C     -------------------------
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO( 1 ) .LT. 0 ) GOTO 500
      IF ( KEEP(46) .eq. 0 .AND. id%MYID.eq.MASTER ) THEN
C       ------------------------------
C       The host is not working -> had
C       no data from initial matrix
C       ------------------------------
        NSEND8  = 0_8
        NLOCAL8 = 0_8
      END IF
C     --------------------------
C     Put into some info/infog ?
C     --------------------------
      CALL MPI_REDUCE( NSEND8, NSEND_TOT8, 1, MPI_INTEGER8,
     &                 MPI_SUM, MASTER, id%COMM, IERR )
      CALL MPI_REDUCE( NLOCAL8, NLOCAL_TOT8, 1, MPI_INTEGER8,
     &                 MPI_SUM, MASTER, id%COMM, IERR )
      IF ( PROKG ) THEN
        WRITE(MPG,125) NLOCAL_TOT8, NSEND_TOT8
      END IF
C     }
      ENDIF ! distributed matrix
C     }
      ENDIF !  "old" multithreaded algorithm 
C     }
      ELSE
C     {
C       -------------------
C       Matrix is elemental,
C       provided on the
C       master only
C       -------------------
        IF ( id%MYID.eq.MASTER)
     &   CALL ZMUMPS_MAXELT_SIZE( id%ELTPTR(1),
     &                        id%NELT,
     &                        MAXELT_SIZE )
C
C         Perform the distribution of the elements.
C         A this point,
C           PTRAIW/PTRARW have been computed.
C           INTARR/DBLARR have been allocated
C           ELTPROC gives the mapping of elements
C
        CALL ZMUMPS_ELT_DISTRIB( id%N, id%NELT, id%KEEP8(30),
     &     id%COMM, id%MYID,
     &     id%NSLAVES, id%PTRAR(1),
     &     id%PTRAR(id%NELT+2),
     &     INTARR(1), DBLARR(1), id%KEEP8(27), id%KEEP8(26),
     &     id%KEEP(1), id%KEEP8(1), MAXELT_SIZE,
     &     id%FRTPTR(1), id%FRTELT(1),
     &     S_PTR_ARG(1), MAXS_ARG, id%FILS(1),
     &     id, idintr%root, idintr%roota )
C       ----------------
C       Broadcast errors
C       ----------------
        CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
        IF ( id%INFO( 1 ) .LT. 0 ) GOTO 500
C     }
      END IF ! Element entry
C     ------------------------
C     Time the redistribution:
C     ------------------------
      IF ( id%MYID.EQ.MASTER) THEN
        CALL MUMPS_SECFIN(TIME)
        id%DKEEP(93) = TIME
        IF (PROKG) WRITE(MPG,160) id%DKEEP(93)
      END IF
C     -------------------------------------
C     Small memory optimizaiton: we can now
C     free RG2L on the non working host,
C     -------------------------------------
      IF (id%KEEP(38) .NE. 0 .AND. .NOT. I_AM_SLAVE) THEN
          IF (associated(idintr%root%RG2L)) THEN
            DEALLOCATE(idintr%root%RG2L)
            NULLIFY   (idintr%root%RG2L)
          ENDIF
      ENDIF
#if defined(MPI_TO_K_OMP)
C{
#if defined(MTKO_VERBOSE)
C{
C     Print only if PROK
      IF (id%KEEP(381).GT.1.AND.PROK) THEN
C      original affinity mask not updated yet
       SETMASK     = 0
       thread_num  = 0
!$OMP  PARALLEL private(thread_num)
!$       thread_num = OMP_GET_THREAD_NUM()
!$OMP    CRITICAL
         CALL MUMPS_PRINT_THREAD_MASK(id%MYID,
     &                       thread_num, SETMASK)
!$OMP    END CRITICAL
!$OMP  END PARALLEL
      ENDIF
      MY_WORKING_ID = -1
C}
#endif
C
C Set affinity mask on running procs if asked
C ===========================================
C
      IF ( I_AM_SLAVE ) THEN
        I_AM_HEAD = .TRUE.
C       will be retset to false on non-working slaves
      ELSE
C       false on non-working master
        I_AM_HEAD = .FALSE.
      ENDIF
      IF (id%KEEP(381).gt.1) THEN
        IF (id%KEEP(382).eq.1) THEN
          write(*,*) '[MPI_TO_K_OMP] Unexpected situation: the mask ',
     &         'should have been unset by mumps_driver at the end.'
        ENDIF
        IF (I_AM_SLAVE) THEN
          CALL MPI_COMM_RANK(id%COMM_NODES, MY_WORKING_ID, IERR)
          CALL MUMPS_MPITOOMP_PRE_WORK_MAPS(id%KEEP(381),
     &            id%NSLAVES,
     &            NCORES,
     &            ERROR)
          IF (ERROR.ne.0) THEN
            id%INFO(1) = -62
            id%INFO(2) = id%MYID
          ENDIF
          allocok = 0
          allocate(CORES_LIST(0:NCORES-1), STAT=allocok)
          IF (allocok .ne. 0) THEN
            id%INFO(1) = -62
            id%INFO(2) = id%MYID
          ENDIF
        ENDIF
C       error check
        CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &             id%COMM, id%MYID)
        IF ( id%INFO( 1 ) .LT. 0) GOTO 500
        IF (I_AM_SLAVE) THEN
          CALL MUMPS_MPITOOMP_GET_WORK_MAPS(id%COMM_NODES,
     &            id%KEEP(381), id%NSLAVES, NCORES,
     &            id%MTKO_PROCS_MAP, CORES_LIST, ERROR)
          IF (ERROR.ne.0) THEN
            id%INFO(1) = -62
            id%INFO(2) = id%MYID
          ENDIF
        ENDIF
C       error check
        CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &             id%COMM, id%MYID)
        IF ( id%INFO( 1 ) .LT. 0) GOTO 500
        IF (I_AM_SLAVE) THEN
          IF (MY_WORKING_ID .EQ.
     &        id%MTKO_PROCS_MAP(MY_WORKING_ID)) THEN
            I_AM_HEAD = .TRUE.
            CALL MUMPS_MPITOOMP_SET_MASK(NCORES, id%KEEP(249),
     &                CORES_LIST, .true., ERROR)
            IF (ERROR.ne.0) THEN
              id%INFO(1) = -61
              id%INFO(2) = id%MYID
            ENDIF
          ELSE
            I_AM_HEAD = .FALSE.
          ENDIF
        ENDIF
C       error check
        CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
        IF ( id%INFO( 1 ) .LT. 0 ) GOTO 500
C       The masks have been successfully set on the workers.
        IF (I_AM_SLAVE) THEN
          IF (MY_WORKING_ID .EQ.
     &      id%MTKO_PROCS_MAP(MY_WORKING_ID)) THEN
            id%KEEP(382) = 1
          ELSE
            id%KEEP(382) = -1
          ENDIF
        ELSE
          id%KEEP(382) = -1
        ENDIF
      ENDIF
#if defined(MTKO_VERBOSE)
C{
C     Print only if PROK
      IF (id%KEEP(382).eq.1.AND.PROK) THEN
C     displays new mask only for working MPIs
        SETMASK     = 1
        thread_num  = 0
!$OMP PARALLEL private(thread_num)
!$      thread_num = OMP_GET_THREAD_NUM()
!$OMP   CRITICAL
        CALL MUMPS_PRINT_THREAD_MASK(id%MYID,
     &                       thread_num, SETMASK)
!$OMP   END CRITICAL
!$OMP END PARALLEL
      ENDIF
C}
#endif
#endif
       IF ( KEEP(400) .GT. 0 .AND. KEEP(369).EQ.0) THEN
C{       Check if number of threads is consistent with
C        the one used during analysis for all procs
C           Note that if  KEEP(369)>0 
C           KEEP(400) was set based on 
C           KEEP(369) and KEEP(381) so that 
C           omp_set_num_threads(KEEP(400)) will be called 
C           explicitly before L0_OMP section and thus
C           KEEP(400) cannot be checked here in this way
        NOMP=1
!$      NOMP = omp_get_max_threads()
#if defined(MPI_TO_K_OMP)
C       NOMP takes K_OMP factor into account
C       (#threads was changed on MPI_TO_K_OMP masters)
C       but test should then only be performed on
C       MPI_TO_K_OMP masters!
        IF (I_AM_HEAD) THEN
#endif
          IF ( NOMP .NE. KEEP(400) ) THEN
            id%INFO(1)=-58
            id%INFO(2)=KEEP(400)
            IF (LPOK) WRITE(LP,'(A,A,I5,A,I5)')
     &" FAILURE DETECTED IN FACTORIZATION: #threads for multithreaded",
     &" tree parallelism changed from",KEEP(400)," at analysis to", 
     &  NOMP
          ENDIF
#if defined(MPI_TO_K_OMP)
        ENDIF
#endif
C}
       ENDIF
C      error check done outside previous if bloc 
C      because KEEP(369) might be 0 on some and nonzero on some proc
       CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
       IF ( id%INFO( 1 ) .LT. 0 ) GOTO 500
C
C     TIMINGS:
C     Next line: elapsed time for factorization
      IF (id%MYID.EQ.MASTER)  CALL MUMPS_SECDEB(TIME)
C
C  Allocate buffers on the workers
C  ===============================
C
      IF ( I_AM_SLAVE )  THEN
C
C  Some buffers are required to pack/unpack data and for
C  receiving MPI messages.
C  For packing/unpacking : the buffer must be large
C  enough to send several messages while receives might not
C  be posted yet.
C  It is assumed that the size of an integer is held in KEEP(34)
C  while the size of a complex is held in KEEP(35).
C  BUFR and LBUFR are declared of type integer, since byte is not
C  a standard datatype.
C  We now use KEEP(43) or KEEP(379) and KEEP(44) or KEEP(380)
C  as estimated at analysis to allocate appropriate buffer sizes
C
C  Receive buffer
C  --------------
        IF (KEEP(486).NE.0) THEN
          ZMUMPS_LBUFR_BYTES8 = int(KEEP( 380 ),8) * int(KEEP(35),8)
        ELSE
          ZMUMPS_LBUFR_BYTES8 = int(KEEP( 44 ),8) * int(KEEP(35),8)
        ENDIF
C       ---------------------------------------
C       Ensure a reasonable minimal buffer size
C       ---------------------------------------
        IF (KEEP(72).NE.1) THEN
C        ensure minimum size for small problems
         ZMUMPS_LBUFR_BYTES8 = max( ZMUMPS_LBUFR_BYTES8,
     &                      200000_8 )
        ENDIF
C
C  If there is pivoting, size of the message might still increase.
C  We use a relaxation (so called PERLU) to increase the estimate.
C
C  Note: PERLU is a global estimate for pivoting. 
C  It may happen that one large contribution block size is increased 
C  by more than that.
C  This is why we use an extra factor 2 relaxation coefficient for 
C  the relaxation of
C  the reception buffer in the case where pivoting is allowed.
C  A more dynamic strategy could be applied: if message to
C  be received is larger than expected, reallocate a larger
C  buffer. (But this won't work with IRECV.)
C  Finally, one may want (as we are currently doing it for 
C  most messages)
C  to cut large messages into a series of smaller ones.
C
        IF (KEEP(48).EQ.5) THEN
          MIN_PERLU = 2
        ELSE
          MIN_PERLU = 0
        ENDIF
C
        IF (KEEP(72).NE.1) THEN
          ZMUMPS_LBUFR_BYTES8 = ZMUMPS_LBUFR_BYTES8
     &        + int( dble(max(PERLU/2,MIN_PERLU))*
     &        dble(ZMUMPS_LBUFR_BYTES8)/100D0, 8)
        ELSE
C         on small pb we want to relax buffers
C         for pivoting 
          ZMUMPS_LBUFR_BYTES8 = ZMUMPS_LBUFR_BYTES8
     &        + int( dble(max(PERLU,MIN_PERLU))*
     &        dble(ZMUMPS_LBUFR_BYTES8)/100D0, 8)
        ENDIF
        ZMUMPS_LBUFR_BYTES8 = min(ZMUMPS_LBUFR_BYTES8,
     &                            int(huge(I4)-100,8))
        ZMUMPS_LBUFR_BYTES  = int( ZMUMPS_LBUFR_BYTES8 )
C       ZMUMPS_LBUFR is the size of the buffer as a number of integers,
C       we round ZMUMPS_LBUFR (size in #integers) above to have at least
C       ZMUMPS_LBUFR_BYTES available in the buffer.
        ZMUMPS_LBUFR = (ZMUMPS_LBUFR_BYTES+KEEP(34)-1)/KEEP(34)
C       Finally, make ZMUMPS_LBUFR_BYTES a multiple of KEEP(34) by setting
C       ZMUMPS_LBUFR_BYTES to the number of bytes that will be allocated
        ZMUMPS_LBUFR_BYTES = ZMUMPS_LBUFR*KEEP(34)
        IF (KEEP(48)==5) THEN
C          Since the buffer is going to be allocated, use
C          it as the constraint for memory/granularity
C          in hybrid scheduler
C
           id%KEEP8(21) = id%KEEP8(22) +
     &        int( dble(max(PERLU/2,MIN_PERLU))*
     &        dble(id%KEEP8(22))/100D0,8)
        ENDIF
C
C  Now estimate the size for the buffer for asynchronous
C  sends of contribution blocks (so called CB). We want to be able to send at
C  least KEEP(213)/100 (two in general) messages at the
C  same time.
C
C   Send buffer
C   -----------
        IF (KEEP(486).NE.0) THEN
         ZMUMPS_LBUF8 = int( dble(KEEP(213)) / 100.0D0 *
     &                      dble(KEEP(379)) * dble(KEEP(35)), 8  )
        ELSE
         ZMUMPS_LBUF8 = int( dble(KEEP(213)) / 100.0D0 *
     &                      dble(KEEP(43)) * dble(KEEP(35)), 8  )
        ENDIF
        IF (KEEP(72).NE.1) THEN
C         ensure minimum size for small problems
          ZMUMPS_LBUF8 = max( ZMUMPS_LBUF8, 200000_8 )
          ZMUMPS_LBUF8 = ZMUMPS_LBUF8
     &                 + int( dble(max(PERLU/2,MIN_PERLU))*
     &                   dble(ZMUMPS_LBUF8)/100D0, 8)
        ELSE
C         for very small pb force extra relaxation 
          ZMUMPS_LBUF8 = ZMUMPS_LBUF8
     &                 + int( dble(max(PERLU,MIN_PERLU))*
     &                   dble(ZMUMPS_LBUF8)/100D0, 8)
        ENDIF
C       Make ZMUMPS_LBUF8 small enough to be stored in a standard integer
        ZMUMPS_LBUF8 = min(ZMUMPS_LBUF8, int(huge(I4)-100,8))
C
C       No reason to have send buffer smaller than receive buffer.
C       This should never occur with the formulas above but just
C       in case:
        ZMUMPS_LBUF8 = max(ZMUMPS_LBUF8, ZMUMPS_LBUFR_BYTES8+3*KEEP(34))
        ZMUMPS_LBUF  = int(ZMUMPS_LBUF8)
        IF(id%KEEP(48).EQ.4)THEN
           ZMUMPS_LBUFR_BYTES=ZMUMPS_LBUFR_BYTES*5
           ZMUMPS_LBUF=ZMUMPS_LBUF*5
        ENDIF
C
C  Estimate size of buffer for small messages 
C  Each node can send ( NSLAVES - 1 ) messages to (NSLAVES-1) nodes
C
C  KEEP(56) is the number of nodes of level II.
C  Messages will be sent for the symmetric case
C  for synchronisation issues.
C
C  We take an upperbound
C
        ZMUMPS_LBUF_INT = ( KEEP(56) + id%NSLAVES * id%NSLAVES ) * 5
     &               * KEEP(34)
        IF ( KEEP( 38 ) .NE. 0 ) THEN
C
C
          KKKK = MUMPS_PROCNODE( id%PROCNODE_STEPS(id%STEP(KEEP(38))),
     &                           id%KEEP(199) )
          IF ( KKKK .EQ. id%MYID_NODES ) THEN
             ZMUMPS_LBUF_INT = ZMUMPS_LBUF_INT + 4 * KEEP(34) *
     &         ( id%NSLAVES + id%NE_STEPS(id%STEP(KEEP(38)))
     &      + min(KEEP(56), id%NE_STEPS(id%STEP(KEEP(38)))) * id%NSLAVES
     &         )
          END IF
        END IF
C       At this point, ZMUMPS_LBUFR_BYTES, ZMUMPS_LBUF
C       and ZMUMPS_LBUF_INT have been computed (all
C       are in numbers of bytes).
        IF ( PROK ) THEN
          WRITE( MP, 9999 ) ZMUMPS_LBUFR_BYTES,
     &                      ZMUMPS_LBUF, ZMUMPS_LBUF_INT
        ELSE 
         IF (PROKG) THEN
          WRITE( MPG, 9999 ) ZMUMPS_LBUFR_BYTES,
     &                      ZMUMPS_LBUF, ZMUMPS_LBUF_INT
         ENDIF
        END IF
 9999   FORMAT( /,' Allocated buffers',/,' ------------------',/,
     &  ' Size of reception buffer in bytes ...... = ', I12,
     &  /,
     &  ' Size of async. emission buffer (bytes).. = ', I12,/,
     &  ' Small emission buffer (bytes) .......... = ', I12)
C       --------------------------
C       Allocate small send buffer
C       required for ZMUMPS_FAC_B
C       --------------------------
        CALL MUMPS_BUF_ALLOC_SMALL_BUF( ZMUMPS_LBUF_INT, IERR )
        IF ( IERR .NE. 0 ) THEN
          id%INFO(1)= -13
C         convert to size in integer  id%INFO(2)= ZMUMPS_LBUF_INT
          id%INFO(2)= (ZMUMPS_LBUF_INT+KEEP(34)-1)/KEEP(34)
          IF (LPOK) THEN
            WRITE(LP,*) id%MYID,
     &     ':Allocation error in MUMPS_BUF_ALLOC_SMALL_BUF'
     &     ,id%INFO(2)
          ENDIF
          GO TO 110
        END IF
C
C       --------------------------------------
C       Allocate reception buffer on all procs
C       This is done now.
C       --------------------------------------
        ALLOCATE( BUFR( ZMUMPS_LBUFR ),stat=IERR )
        IF ( IERR .NE. 0 ) THEN
          id%INFO(1) = -13
          id%INFO(2) = ZMUMPS_LBUFR
          IF (LPOK) THEN
            WRITE(LP,*)
     &     ': Allocation error for BUFR(', ZMUMPS_LBUFR,
     &     ') on MPI process',id%MYID
          ENDIF
          GO TO 110
        END IF
C       -----------------------------------------
C       Estimate MAXIS. IS will be allocated in
C       ZMUMPS_FAC_B. It will contain factors and
C       contribution blocks integer information
C       -----------------------------------------
C       Relax integer workspace based on PERLU
        PERLU          = KEEP( 12 )
        IF (KEEP(201).GT.0) THEN
C         OOC panel or non panel (note that
C         KEEP(15)=KEEP(225) if non panel)
          MAXIS_ESTIM   = KEEP(225)
        ELSE
C         In-core or reals for factors not stored
          MAXIS_ESTIM   = KEEP(15)
        ENDIF
        MAXIS = max( 1, int( min(  int(huge(MAXIS),8),
     &          int(MAXIS_ESTIM,8) + 3_8 * max(int(PERLU,8),10_8) *
     &          ( int(MAXIS_ESTIM,8) / 100_8 + 1_8 )
     &                          )   ! min
     &                     )  ! int
     &             ) !max
C       ----------------------------
C       Allocate PTLUST_S and PTRFAC
C       They will be used to access
C       factors in the solve phase.
C       They are also needed for
C       ZMUMPS_FAC_L0_OMP.
C       ----------------------------
        ALLOCATE( id%PTLUST_S( id%KEEP(28) ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
          id%INFO(1)=-13
          id%INFO(2)=id%KEEP(28)
          IF (LPOK) THEN
            WRITE(LP,*) id%MYID,
     &      ': Allocation error for id%PTLUST_S(', id%KEEP(28),')'
          ENDIF
          NULLIFY(id%PTLUST_S)
          GOTO 110
        END IF
        ALLOCATE( id%PTRFAC( id%KEEP(28) ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
          id%INFO(1)=-13
          id%INFO(2)=id%KEEP(28)
          NULLIFY(id%PTRFAC)
          IF (LPOK) THEN
            WRITE(LP,*) id%MYID,
     &      ': Allocation error for id%PTRFAC(', id%KEEP(28),')'
          ENDIF
          GOTO 110
        END IF
C       -----------------------------
C       Reserve temporary workspace :
C       IPOOL, PTRWB, ITLOC, PTRIST
C       PTRWB will be subdivided again
C       in routine ZMUMPS_FAC_B
C       -----------------------------
        PTRIST = 1
        PTRWB  = PTRIST + id%KEEP(28)
        ITLOC  = PTRWB  + 2 * id%KEEP(28)
C Fwd in facto: ITLOC of size id%N + id%KEEP(253)
        IPOOL  = ITLOC  + id%N + id%KEEP(253)
C
C       --------------------------------
C       NA(1) is an upperbound for LPOOL
C       --------------------------------
C       Structure of the pool:
C     ____________________________________________________
C    | Subtrees   |         | Top nodes           | 1 2 3 |
C     ----------------------------------------------------
        LPOOL = MUMPS_GET_POOL_LENGTH(id%NA(1), id%KEEP(1),id%KEEP8(1))
        LIWK  = IPOOL + LPOOL - 1
        ALLOCATE( IWK(  IPOOL + LPOOL - 1 ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
          id%INFO(1)=-13
          id%INFO(2)=IPOOL + LPOOL - 1
          IF (LPOK) THEN
            WRITE(LP,*) id%MYID,
     &      ': Allocation error for IWK(',IPOOL+LPOOL-1,')'
          ENDIF
          GOTO 110
        END IF
        LIWK8 =  2 * id%KEEP(28)
        ALLOCATE(IWK8( 2 * id%KEEP(28)), stat = IERR)
        IF ( IERR .NE. 0 ) THEN
          id%INFO(1)=-13
          id%INFO(2)=2 * id%KEEP(28)
          IF (LPOK) THEN
            WRITE(LP,*) id%MYID,
     &      ': Allocation error for IWKB(', 2*id%KEEP(28),')'
          ENDIF
          GOTO 110
        END IF
C
C  Return to SPMD
C
      ENDIF
C
 110  CONTINUE
C     ----------------
C     Broadcast errors
C     ----------------
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO( 1 ) .LT. 0 ) GOTO 500
C
      IF ( I_AM_SLAVE )  THEN
C       Store size of receive buffers in ZMUMPS_LBUF module
        CALL MUMPS_BUF_DIST_IRECV_SIZE( ZMUMPS_LBUFR_BYTES )
#if defined(USE_GPU)
#endif
        IF (PROK) THEN
          WRITE( MP, 170 ) MAXS, MAXIS, MAXS_BASE8, KEEP(15),
     &    id%KEEP8(26), id%KEEP8(27), id%KEEP8(11), KEEP(26), KEEP(27)
        ENDIF
      END IF
C     ===============================================================
C     Before calling the main driver, ZMUMPS_FAC_B,
C     some statistics should be initialized to 0,
C     even on the host node because they will be
C     used in REDUCE operations afterwards.
C     --------------------------------------------
C     Size of factors written. It will be set to POSFAC in
C     IC, otherwise we accumulate written factors in it.
      id%KEEP8(31)= 0_8
C     Size of factors under L0 will be returned
C     in id%KEEP8(64), not included in KEEP8(31))
C     Number of entries in factors
      id%KEEP8(10) = 0_8
C     KEEP8(8) will hold the volume of extra copies due to
C              in-place stacking in fac_mem_stack.F
      id%KEEP8(8)=0_8
      id%INFO(9:14)=0
      RINFO(2:3)=ZERO
      IF ( I_AM_SLAVE ) THEN
C       ------------------------------------
C       Call effective factorization routine
C       ------------------------------------
        IF ( KEEP(55) .eq. 0 ) THEN
          LDPTRAR = 1 ! PTRAR no longer used (of size 2)
        ELSE
          LDPTRAR = id%NELT + 1
        END IF
        IF ( id%KEEP(55) .NE. 0 ) THEN
          NELT_arg = id%NELT
        ELSE
C         ------------------------------
C         Use size 1 to avoid complaints
C         when using check bound options
C         ------------------------------
          NELT_arg = 1
        END IF
      ENDIF
      IF (I_AM_SLAVE) THEN
        IF (associated(id%L0_OMP_MAPPING))
     &    DEALLOCATE(id%L0_OMP_MAPPING)
        IF (KEEP(400) .GT. 0) THEN
          id%LL0_OMP_MAPPING = KEEP(28)
        ELSE
          id%LL0_OMP_MAPPING = 1
        ENDIF
        ALLOCATE(id%L0_OMP_MAPPING(id%LL0_OMP_MAPPING), stat=allocok)
        IF ( allocok > 0) THEN
            write(*,*) "Problem allocating L0_OMP_MAPPING",
     &              IERR, KEEP(28)
            GOTO 115
        ENDIF
        IF (KEEP(400) .GT. 0) THEN
          id%LL0_OMP_FACTORS = KEEP(400)
        ELSE
          id%LL0_OMP_FACTORS = 1
        ENDIF
#       if defined(SINGLEINDOUBLE)
        IF (KEEP(148).EQ.0) THEN
#       endif
        ALLOCATE(idintr%L0_OMP_FACTORS(id%LL0_OMP_FACTORS),
     &           stat = allocok)
#       if defined(SINGLEINDOUBLE)
        ELSE
        ALLOCATE(idintr%LPL0_OMP_FACTORS(id%LL0_OMP_FACTORS),
     &           stat = allocok)
        ENDIF
#       endif
        IF (allocok > 0) THEN
           id%INFO(1)=-7
           id%INFO(2)=NB_THREADS
           GOTO 115
        ENDIF
#       if defined(SINGLEINDOUBLE)
        IF (KEEP(148).EQ.0) THEN
#       endif
        CALL ZMUMPS_INIT_L0_OMP_FACTORS(idintr%L0_OMP_FACTORS)
#       if defined(SINGLEINDOUBLE)
        ELSE
        CALL CMUMPS_INIT_L0_OMP_FACTORS(idintr%LPL0_OMP_FACTORS)
        ENDIF
#       endif
      ENDIF
 115  CONTINUE
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO( 1 ) .LT. 0 ) GOTO 500
C     Compute DKEEP(17) 
      AVG_FLOPS    =  RINFOG(1)/(dble(id%NSLAVES))
      id%DKEEP(17) = max ( id%DKEEP(18), AVG_FLOPS/dble(50) )
     &               
      IF (PROK.AND.id%MYID.EQ.MASTER) THEN
        IF (id%NSLAVES.LE.1) THEN
         WRITE(MP,'(/A,A,1PD10.3)')
     &' Start factorization with total', 
     &' estimated flops (RINFOG(1))                         = ', 
     &   RINFOG(1)
        ELSE
         WRITE(MP,'(/A,A,1PD10.3,A,1PD10.3)')
     &' Start factorization with total', 
     &' estimated flops RINFOG(1) / Average per MPI proc    = ', 
     &      RINFOG(1), ' / ', AVG_FLOPS
        ENDIF
      ENDIF
      IF (I_AM_SLAVE) THEN
C       IS/S pointers passed to ZMUMPS_FAC_B with
C       implicit interface through intermediate
C       structure S_IS_POINTERS. IS will be allocated
C       during ZMUMPS_FAC_B.
C       In case of L0OMP, id%IS and id%S are allocated during
C       ZMUMPS_FAC_B, and only after L0OMP nodes are processed,
C       in order to limit the global memory peak. 
#if defined(SINGLEINDOUBLE)
        IF (KEEP(148).EQ.0) THEN
#endif
        S_IS_POINTERS%IW => id%IS; NULLIFY(id%IS)
        S_IS_POINTERS%A  => id%S ; NULLIFY(id%S)
        CALL ZMUMPS_FAC_B(id%N,S_IS_POINTERS,MAXS,MAXIS,id%SYM_PERM(1),
     &  id%NA(1),id%LNA,id%NE_STEPS(1),id%ND_STEPS(1), id%FILS(1),
     &  id%STEP(1),id%FRERE_STEPS(1),id%DAD_STEPS(1),id%CANDIDATES(1,1),
     &  id%ISTEP_TO_INIV2(1),id%TAB_POS_IN_PERE(1,1), id%PTRAR(1),
     &LDPTRAR,id%PTR8ARR(1),id%NINCOLARR(1),id%NINROWARR(1),id%PTRDEBARR
     &  (1), IWK(PTRIST),id%PTLUST_S(1),id%PTRFAC(1),IWK(PTRWB),IWK8,
     &  IWK(ITLOC),RHS_MUMPS(1),IWK(IPOOL),LPOOL,CNTL1,ICNTL(1),
     &  id%INFO(1), RINFO(1),KEEP(1),id%KEEP8(1),id%PROCNODE_STEPS(1),
     &  id%NSLAVES,id%COMM_NODES,id%MYID,id%MYID_NODES,BUFR,ZMUMPS_LBUFR
     &  , ZMUMPS_LBUFR_BYTES, ZMUMPS_LBUF, INTARR(1), DBLARR(1),
     &  idintr%root, idintr%roota, NELT_arg, id%FRTPTR(1), id%FRTELT(1),
     &  id%COMM_LOAD,id%ASS_IRECV,SEUIL,SEUIL_LDLT_NIV2,id%MEM_DIST(0),
     &  id%DKEEP(1), PIVNUL_LIST_STRUCT, id%LRGROUPS(1)
     &  ,id%IPOOL_B_L0_OMP(1),id%LPOOL_B_L0_OMP,
     &  id%IPOOL_A_L0_OMP(1),id%LPOOL_A_L0_OMP,id%L_VIRT_L0_OMP,
     &  id%VIRT_L0_OMP(1), id%VIRT_L0_OMP_MAPPING(1),id%L_PHYS_L0_OMP,
     &  id%PHYS_L0_OMP(1), id%PERM_L0_OMP(1), id%PTR_LEAFS_L0_OMP(1),
     &  id%L0_OMP_MAPPING(1),id%LL0_OMP_MAPPING, id%THREAD_LA,
     &  idintr%L0_OMP_FACTORS(1), id%LL0_OMP_FACTORS,
     &  id%I4_L0_OMP(1,1), size(id%I4_L0_OMP,1), size(id%I4_L0_OMP,2),
     &  id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1), size(id%I8_L0_OMP,2)
#if defined(MPI_TO_K_OMP)
     &     ,I_AM_HEAD
#endif
     &     )
        id%IS => S_IS_POINTERS%IW; NULLIFY(S_IS_POINTERS%IW)
        id%S  => S_IS_POINTERS%A ; NULLIFY(S_IS_POINTERS%A)
#if defined(SINGLEINDOUBLE)
        ELSE
C       Arith-dependent parameters are:
C         S_IS_POINTERS => LPS_IS_POINTERS
C
C         idintr%roota => lproota
C
C
C         RINFO, DKEEP => translated in and out
C                         in-place avoids allocations and is OK since
C                         they are small
C         DBLARR, input only: copied to single precision, freeing
C                             DBLARR early when possible, i.e.
C                             when allocated in this routine.
C                Better than in-place copy to avoid memory peak
C
C
C       LPS_IS_POINTERS is a wrapper around IS and LPS
        LPS_IS_POINTERS%IW => id%IS; NULLIFY(id%IS)
        LPS_IS_POINTERS%A  => id%LPS ; NULLIFY(id%LPS)
C       Compress RINFO and DKEEP from double to single precision
        CALL MUMPS_RCOPY_64TO32_64C_IP( RINFO(1),
     &                                  int(size(RINFO),kind=8) )
        CALL MUMPS_RCOPY_64TO32_64C_IP( id%DKEEP(1),
     &                                  int(size(id%DKEEP),kind=8) )
C       arrowheads or elements
        ALLOCATE(LPDBLARR(max(id%KEEP8(26),1_8)), stat=allocok)
        IF (allocok .GT. 0) THEN
          id%INFO(1) = -13
          CALL MUMPS_SET_IERROR(max(id%KEEP8(26),1_8), id%INFO(2))
          GOTO 107
        ENDIF
C       2_8* because MUMPS_RCOPY_64TO32_64C expects reals
        CALL MUMPS_RCOPY_64TO32_64C(DBLARR(1), 2_8*id%KEEP8(26),
     &                              LPDBLARR(1) )
C       Free DBLARR before memory peak in CMUMPS_FAC_B
        IF (DBLARR_ALLOCATED) THEN
          DEALLOCATE(DBLARR)
          NULLIFY(DBLARR)
          DBLARR_ALLOCATED = .FALSE.
        ENDIF
C       RHS_MUMPS -> LPRHS_MUMPS
        ALLOCATE(LPRHS_MUMPS(id%KEEP8(85)), stat=allocok)
        IF (allocok .GT. 0) THEN
          id%INFO(1) = -13
          CALL MUMPS_SET_IERROR(id%KEEP8(85), id%INFO(2))
          GOTO 107
        ENDIF
C       Copy done column by column, in case of master
C       with leading dimension higher than N
        I8=int(N,8)
        I8=I8+I8
        DO K=1, id%KEEP(253)
          CALL MUMPS_RCOPY_64TO32_64C(RHS_MUMPS(1_8+int(K-1,8)*
     &    int(KEEP(254),8)),
     &    I8, LPRHS_MUMPS(1_8+int(K-1,8)*int(id%KEEP(254),8)))
        ENDDO
        IF (RHS_MUMPS_ALLOCATED) THEN
C         Double precision version of RHS_MUMPS will no longer be used
          DEALLOCATE(RHS_MUMPS)
          NULLIFY(RHS_MUMPS)
          RHS_MUMPS_ALLOCATED=.FALSE.
        ENDIF
C       Nullify pointers in roota:
        CALL CMUMPS_INI_ROOT(lproota)
 107    CONTINUE
C     Exit IF I_AM_SLAVE block to call MUMPS_PROPINFO
      ENDIF ! KEEP(148)
      ENDIF ! I_AM_SLAVE
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO( 1 ) .LT. 0 ) THEN
C       LPDBLARR and LPRHS_MUMPS are allocated and freed only around
C       the call to CMUMPS_FAC_B and not at the end of the
C       ZMUMPS_FAC_DRIVER. So we free them now in case of error.
        IF (allocated(LPDBLARR)) DEALLOCATE(LPDBLARR)
        IF (allocated(LPRHS_MUMPS)) DEALLOCATE(LPRHS_MUMPS)
        GOTO 500
      ENDIF
      IF (I_AM_SLAVE) THEN
      IF (KEEP(148).NE.0) THEN
C       CNTL1 => casted to single implicitly
        LPCNTL1 = CNTL1
        LPSEUIL = SEUIL
        LPSEUIL_LDLT_NIV2 = SEUIL_LDLT_NIV2
      CALL CMUMPS_FAC_B(id%N,LPS_IS_POINTERS,MAXS,MAXIS,id%SYM_PERM(1)
     &  , id%NA(1),id%LNA,id%NE_STEPS(1),id%ND_STEPS(1), id%FILS(1),
     &  id%STEP(1),id%FRERE_STEPS(1),id%DAD_STEPS(1),id%CANDIDATES(1,1),
     &  id%ISTEP_TO_INIV2(1),id%TAB_POS_IN_PERE(1,1), id%PTRAR(1),
     &LDPTRAR,id%PTR8ARR(1),id%NINCOLARR(1),id%NINROWARR(1),id%PTRDEBARR
     &  (1), IWK(PTRIST),id%PTLUST_S(1),id%PTRFAC(1),IWK(PTRWB),IWK8,
     &  IWK(ITLOC),LPRHS_MUMPS(1),IWK(IPOOL),LPOOL,LPCNTL1,ICNTL(1),
     &  id%INFO(1), RINFO(1),KEEP(1),id%KEEP8(1),id%PROCNODE_STEPS(1),
     &  id%NSLAVES,id%COMM_NODES,id%MYID,id%MYID_NODES,BUFR,ZMUMPS_LBUFR
     &  , ZMUMPS_LBUFR_BYTES, ZMUMPS_LBUF, INTARR(1), LPDBLARR(1),
     &idintr%root, lproota, NELT_arg, id%FRTPTR(1), id%FRTELT(1),
     &  id%COMM_LOAD,id%ASS_IRECV,LPSEUIL,LPSEUIL_LDLT_NIV2,
     &id%MEM_DIST(0), id%DKEEP(1), PIVNUL_LIST_STRUCT, id%LRGROUPS(1)
     &  ,id%IPOOL_B_L0_OMP(1),id%LPOOL_B_L0_OMP,
     &  id%IPOOL_A_L0_OMP(1),id%LPOOL_A_L0_OMP,id%L_VIRT_L0_OMP,
     &  id%VIRT_L0_OMP(1), id%VIRT_L0_OMP_MAPPING(1),id%L_PHYS_L0_OMP,
     &  id%PHYS_L0_OMP(1), id%PERM_L0_OMP(1), id%PTR_LEAFS_L0_OMP(1),
     &  id%L0_OMP_MAPPING(1),id%LL0_OMP_MAPPING, id%THREAD_LA,
     &  idintr%LPL0_OMP_FACTORS(1), id%LL0_OMP_FACTORS,
     &  id%I4_L0_OMP(1,1), size(id%I4_L0_OMP,1), size(id%I4_L0_OMP,2),
     &  id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1), size(id%I8_L0_OMP,2)
#if defined(MPI_TO_K_OMP)
     &     ,I_AM_HEAD
#endif
     &     )
C       Free temporary arrrowheads
        DEALLOCATE(LPDBLARR)
C       Free temporary LPRHS_MUMPS
        DEALLOCATE(LPRHS_MUMPS)
C       Extend back RINFO and DKEEP from single to double precision
        CALL MUMPS_RCOPY_32TO64_64C_IP( RINFO(1),
     &                                  int(size(RINFO),kind=8) )
        CALL MUMPS_RCOPY_32TO64_64C_IP( id%DKEEP(1),
     &                                  int(size(id%DKEEP),kind=8) )
        id%IS   => LPS_IS_POINTERS%IW; NULLIFY(LPS_IS_POINTERS%IW)
        id%LPS  => LPS_IS_POINTERS%A ; NULLIFY(LPS_IS_POINTERS%A)
C       lproota => roota
C       RHS_CNTR_MASTER_ROOT: allocated in ZMUMPS_FAC_PAR
        IF (associated(lproota%RHS_CNTR_MASTER_ROOT)) THEN
#if defined(MUMPS_NOF2003)
          I8 = int(size(lproota%RHS_CNTR_MASTER_ROOT),8)
#else
          I8 = size(lproota%RHS_CNTR_MASTER_ROOT,kind=8)
#endif
          IF (I8.EQ.0_8) THEN
            WRITE(*,*)"Internal error I8=0 in fac_driver.F"
            CALL MUMPS_ABORT()
          ENDIF
          ALLOCATE(idintr%roota%RHS_CNTR_MASTER_ROOT(I8),
     &             stat=allocok)
          IF (allocok .GT. 0) THEN
            id%INFO(1)=-13
            CALL MUMPS_SET_IERROR(I8, id%INFO(2))
          ELSE
C           2 reals for 1 complex
            I8 = I8 + I8
            CALL MUMPS_RCOPY_32TO64_64C(
     &         lproota%RHS_CNTR_MASTER_ROOT(1), I8,
     &         idintr%roota%RHS_CNTR_MASTER_ROOT(1)
     &         )
          ENDIF
        ENDIF
C       Free temporary lproota
        CALL CMUMPS_END_ROOT(lproota)
      ENDIF
#endif
C
C     ------------------------------
C     Deallocate temporary workspace
C     ------------------------------
      DEALLOCATE( IWK  )
      DEALLOCATE( IWK8 )
      ENDIF
C     Fwd in facto: free RHS_MUMPS in case it was allocated.
      IF (RHS_MUMPS_ALLOCATED) THEN
        DEALLOCATE(RHS_MUMPS)
        RHS_MUMPS_ALLOCATED = .FALSE.
      ENDIF
      NULLIFY(RHS_MUMPS)
C     ---------------------------------
C     Free some workspace corresponding
C     to the original matrix in
C     arrowhead or elemental format.
C                  -----
C     Note : DBLARR may be a pointer
C     in case of element-entry.
C     ---------------------------------
      IF (allocated( INTARR )) DEALLOCATE( INTARR )
      IF (DBLARR_ALLOCATED) THEN
        DEALLOCATE(DBLARR)
        DBLARR_ALLOCATED=.FALSE.
      ENDIF
      NULLIFY(DBLARR)
C     We also free RG2L now
      IF ( KEEP(38) .NE. 0) THEN
        IF (associated(idintr%root%RG2L)) THEN
          DEALLOCATE(idintr%root%RG2L)
          NULLIFY(idintr%root%RG2L)
        ENDIF
      ENDIF
C
C     Memory statistics
C     -----------------------------------
C     If QR (Keep(19)) is not zero, and if
C     the host does not have the information
C     (ie is not slave), send information
C     computed on the slaves during facto
C     to the host.
C     -----------------------------------
C     Note the KEEP(17), KEEP(143) have been bcasted during fac_par_m
      IF ( KEEP(19) .NE. 0 ) THEN
        IF ( KEEP(46) .NE. 1 ) THEN
C         Host was not working during facto_root
C         Send him the information
          IF ( id%MYID .eq. MASTER ) THEN
            CALL MPI_RECV( KEEP(17), 1, MPI_INTEGER, 1, DEFIC_TAG,
     &                   id%COMM, STATUS, IERR )
            CALL MPI_RECV( KEEP(143), 1, MPI_INTEGER, 1, DEFIC_TAG,
     &                   id%COMM, STATUS, IERR )
          ELSE IF ( id%MYID .EQ. 1 ) THEN
            CALL MPI_SEND( KEEP(17), 1, MPI_INTEGER, 0, DEFIC_TAG,
     &                   id%COMM, IERR )
            CALL MPI_SEND( KEEP(143), 1, MPI_INTEGER, 0, DEFIC_TAG,
     &                   id%COMM, IERR )
          END IF
        END IF
      END IF
C     --------------------------------
C     Deallocate communication buffers
C     They will be reallocated
C     in the solve.
C     --------------------------------
      IF (allocated(BUFR)) DEALLOCATE(BUFR)
      CALL MUMPS_BUF_DEALL_SMALL_BUF( IERR )
C
C     Check for errors.
C     After ZMUMPS_FAC_B every slave is aware of an error. 
C     The call below informs the master, in case it is not
C     included in the computations.
#if defined(SINGLEINDOUBLE)
C     In case of SINGLEINDOUBLE, an error may also have occurred
C     when allocating idintr%roota%RHS_CNTR_MASTER_ROOT
#endif
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
C
      CALL ZMUMPS_EXTRACT_SCHUR_REDRHS(id,idintr)
C     return to user singular values
      IF (id%KEEP(19) .NE.0) THEN
       CALL ZMUMPS_EXTRACT_SINGULAR_VALUES(id,idintr)
      ENDIF
      IF (KEEP(201) .GT. 0) THEN
      END IF
      IF (id%MYID.EQ.MASTER) THEN
        CALL MUMPS_SECFIN(TIME)
        id%DKEEP(94)=TIME
        IF (KEEP(400).GT.0) THEN
C         Facto time above L0_OMP = total time - facto time under L0_OMP
          id%DKEEP(96)=id%DKEEP(94)-id%DKEEP(95)
        ENDIF
      ENDIF
C     Time to process root node:
      CALL MPI_REDUCE( id%DKEEP(99), TMPTIME, 1,
     &                 MPI_DOUBLE_PRECISION,
     &                 MPI_MAX, MASTER, id%COMM, IERR )
      id%DKEEP(99)=TMPTIME
C     =====================================================================
C     COMPUTE MEMORY ALLOCATED BY MUMPS, INFO(16)
C     ---------------------------------------------
      MEM_EFF_ALLOCATED = .TRUE.
      CALL ZMUMPS_MAX_MEM( id%KEEP(1),id%KEEP8(1),
     &     id%MYID, N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &     id%KEEP8(30),
     &     id%NSLAVES, TOTAL_MBYTES, .TRUE., id%KEEP(201),
     &     BLR_STRAT, .TRUE., TOTAL_BYTES,
     &     IDUMMY, BDUMMY, MEM_EFF_ALLOCATED
     &     , .FALSE.   ! UNDER_L0_OMP
     &    , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &     size(id%I8_L0_OMP,2)
     &     )
      IF (KEEP(400) .GT. 0 ) THEN ! L0 activated
        CALL ZMUMPS_MAX_MEM( id%KEEP(1),id%KEEP8(1),
     &     id%MYID, N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &     id%KEEP8(30),
     &     id%NSLAVES, TOTAL_MBYTES_UNDER_L0, .TRUE., id%KEEP(201),
     &     BLR_STRAT, .TRUE., TOTAL_BYTES_UNDER_L0,
     &     IDUMMY, BDUMMY, MEM_EFF_ALLOCATED
     &     , .TRUE.   ! UNDER_L0_OMP
     &     , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &       size(id%I8_L0_OMP,2)
     &        )
        TOTAL_MBYTES = max (TOTAL_MBYTES,TOTAL_MBYTES_UNDER_L0) 
        TOTAL_BYTES  = max (TOTAL_BYTES, TOTAL_BYTES_UNDER_L0)
      ENDIF
       IF (id%KEEP8(24).NE.0) THEN
C         WK_USER is not part of memory allocated by MUMPS
C         and is not counted, id%KEEP8(23) should be zero
          id%INFO(16) = TOTAL_MBYTES
       ELSE
C         Note that even for the case of ICNTL(23)>0 
C         we report here the memory effectively allocated
C         that can be smaller than ICNTL(23) !
          id%INFO(16) = TOTAL_MBYTES
       ENDIF
C       ----------------------------------------------------
C       Centralize memory statistics on the host
C       id%INFOG(18) = size of mem in Mbytes for facto,
C                   for the processor using largest memory
C       id%INFOG(19) = size of mem in Mbytes for facto,
C                   sum over all processors
C       ----------------------------------------------------
        CALL MUMPS_MEM_CENTRALIZE( id%MYID, id%COMM,
     &                           id%INFO(16), id%INFOG(18), IRANK )
        CALL ZMUMPS_PRINT_ALLOCATED_MEM( PROK, PROKG, PRINT_MAXAVG,
     &       MP, MPG, id%INFO(16), id%INFOG(18), id%INFOG(19),
     &       id%NSLAVES, IRANK,
     &       id%KEEP(1) )
C     If WK_USER is provided, this memory excludes WK_USER
      IF (PROK ) THEN
          WRITE(MP,'(A,I12) ')
     & ' ** Eff. min. Space MBYTES for facto             (INFO(16)):',
     &                TOTAL_MBYTES
      ENDIF
C     ========================(INFO(16) RELATED)======================
C     ---------------------------------------
C     COMPUTE EFFECTIVE MEMORY USED INFO(22)
C     ---------------------------------------
      PERLU_ON = .TRUE.
      MEM_EFF_ALLOCATED = .FALSE.
      CALL ZMUMPS_MAX_MEM( id%KEEP(1),id%KEEP8(1),
     &     id%MYID, N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &     id%KEEP8(30),
     &     id%NSLAVES, TOTAL_MBYTES, .TRUE., id%KEEP(201),
     &     BLR_STRAT, PERLU_ON, TOTAL_BYTES,
     &     IDUMMY, BDUMMY, MEM_EFF_ALLOCATED
     &     , .FALSE.   ! UNDER_L0_OMP
     &     , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &       size(id%I8_L0_OMP,2)
     &     )
      IF (KEEP(400) .GT. 0 ) THEN ! L0 activated
        CALL ZMUMPS_MAX_MEM( id%KEEP(1),id%KEEP8(1),
     &    id%MYID, N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &    id%KEEP8(30),
     &    id%NSLAVES, TOTAL_MBYTES_UNDER_L0, .TRUE., id%KEEP(201),
     &    BLR_STRAT, PERLU_ON, TOTAL_BYTES_UNDER_L0,
     &    IDUMMY, BDUMMY, MEM_EFF_ALLOCATED
     &    , .TRUE. ! UNDER_L0_OMP
     &    , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &      size(id%I8_L0_OMP,2)
     &  )
        TOTAL_MBYTES = max (TOTAL_MBYTES,TOTAL_MBYTES_UNDER_L0) 
        TOTAL_BYTES  = max (TOTAL_BYTES, TOTAL_BYTES_UNDER_L0)
      ENDIF
C     -- TOTAL_BYTES and TOTAL_MBYTES includes both static 
C     -- (MAXS) and BLR structures computed as the SUM of the PEAKS 
C     -- (KEEP8(67) + KEEP8(70))
      id%KEEP8(7) = TOTAL_BYTES
C     -- INFO(22) holds the effective space (in Mbytes) used by MUMPS
C     -- (it includes part of WK_USER used if provided by user)
      id%INFO(22) = TOTAL_MBYTES
C     ----------------------------------------------------
C     Centralize memory statistics on the host
C       INFOG(21) = size of effective mem (Mbytes) for facto,
C                   for the processor using largest memory
C       INFOG(22) = size of effective mem (Mbytes) for facto,
C                   sum over all processors
C     ----------------------------------------------------
      CALL MUMPS_MEM_CENTRALIZE( id%MYID, id%COMM,
     &                    id%INFO(22), id%INFOG(21), IRANK )
      IF ( PROKG ) THEN
       IF (PRINT_MAXAVG) THEN
        WRITE( MPG,'(A,I12) ')
     &  ' ** Memory effectively used, max in  Mbytes     (INFOG(21)):',
     &  id%INFOG(21)
       ENDIF
       WRITE( MPG,'(A,I12) ')
     &  ' ** Memory effectively used, total in Mbytes    (INFOG(22)):',
     &  id%INFOG(22)
      END IF
      SUM_INFO22_THIS_NODE=0
      CALL MPI_REDUCE( id%INFO(22), SUM_INFO22_THIS_NODE, 1,
     &                 MPI_INTEGER,
     &                 MPI_SUM, 0, id%KEEP(411), IERR )
      CALL MPI_REDUCE( SUM_INFO22_THIS_NODE, MAX_SUM_INFO22_THIS_NODE,
     &                 1, MPI_INTEGER, MPI_MAX, 0, id%COMM, IERR )
      IF (PROKG .AND. PRINT_NODEINFO) THEN
        WRITE(MPG,'(A,I12)')
     &  ' ** Max. effective space per compute node, in MBytes       :',
     &  MAX_SUM_INFO22_THIS_NODE
      ENDIF
C
      IF (I_AM_SLAVE) THEN
       K67 = id%KEEP8(67)
       K68 = id%KEEP8(68)
       K70 = id%KEEP8(70)
       K74 = id%KEEP8(74)
       K75 = id%KEEP8(75)
      ELSE
       K67 = 0_8
       K68 = 0_8
       K70 = 0_8
       K74 = 0_8
       K75 = 0_8
      ENDIF
C     -- Save the number of entries effectively used
C        in main working array S
      CALL MUMPS_SETI8TOI4(K67,id%INFO(21))
C
      IF (id%NPROCS .GT. 1 .AND. id%KEEP(50) .NE. 0) THEN
        CALL MPI_REDUCE( id%KEEP8(131), id%KEEP8(132), 1, MPI_INTEGER8,
     &                   MPI_MAX, MASTER, id%COMM, IERR )
        CALL MPI_REDUCE( id%KEEP8(131), id%KEEP8(133), 1, MPI_INTEGER8,
     &                   MPI_SUM, MASTER, id%COMM, IERR )
        CALL MPI_REDUCE( id%KEEP(175), id%KEEP(176), 1, MPI_INTEGER,
     &                   MPI_MAX, MASTER, id%COMM, IERR )
      ENDIF
C
      IF (KEEP(400) .GT.0 ) THEN
        IF (.NOT. I_AM_SLAVE) THEN
          id%DKEEP(95) = 0.0D0
          id%DKEEP(16) = 0.0D0
        ENDIF
        IF (id%NPROCS .GT. 1) THEN
C         Compute average and max (across MPI's)
          CALL MPI_REDUCE(id%DKEEP(95), TMPTIME, 1,
     &    MPI_DOUBLE_PRECISION, MPI_SUM, MASTER, id%COMM, IERR)
          IF (id%MYID.EQ.MASTER) TIMEAVG = TMPTIME
          CALL MPI_REDUCE(id%DKEEP(16), TMPFLOP, 1,
     &    MPI_DOUBLE_PRECISION, MPI_SUM, MASTER, id%COMM, IERR)
          IF (id%MYID.EQ.MASTER) FLOPAVG = TMPFLOP
          IF (id%MYID.EQ.MASTER) THEN
            TIMEAVG = TIMEAVG / id%NSLAVES
            FLOPAVG = FLOPAVG / id%NSLAVES
          ENDIF
          CALL MPI_REDUCE(id%DKEEP(95), TIMEMAX, 1, 
     &    MPI_DOUBLE_PRECISION, MPI_MAX, MASTER, id%COMM, IERR)
          CALL MPI_REDUCE(id%DKEEP(16), FLOPMAX, 1, 
     &    MPI_DOUBLE_PRECISION, MPI_MAX, MASTER, id%COMM, IERR)
C         (PROKG may only be true on master)
          IF ( PROKG ) THEN
            WRITE(MPG,190) FLOPAVG, FLOPMAX
            WRITE(MPG,188) TIMEAVG, TIMEMAX
          ENDIF
        ELSE
C          Print DKEEP(95) directly without reduction
           IF ( PROKG ) THEN
             WRITE(MPG,189) id%DKEEP(16)
             WRITE(MPG,187) id%DKEEP(95)
           ENDIF
        ENDIF
      ENDIF
      IF ( PROKG ) THEN
          IF ( ( KEEP(38).NE.0 .OR. KEEP(20).NE.0 ) .AND.
     &        KEEP(60) .EQ. 0 ) THEN
            WRITE(MPG,186) id%DKEEP(99)
          ENDIF
C         Elapsed time for factorization:
          IF (id%INFO(1) .GE.0) THEN
            WRITE(MPG,180) id%DKEEP(94)
          ELSE
            WRITE(MPG,185) id%DKEEP(94)
          ENDIF
      ENDIF
C
C  Sum RINFO(2) : total number of flops for assemblies
C  Sum RINFO(3) : total number of flops for eliminations
C     Initialize RINFO(4) in case BLR was not activated
      RINFO(4) = RINFO(3)
C 
C  Should work even if the master does some work
C
      CALL MPI_REDUCE( RINFO(2), RINFOG(2), 2,
     &                 MPI_DOUBLE_PRECISION,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
C     Reduce needed to dimension small working array
C     on all procs during ZMUMPS_GATHER_SOLUTION
      KEEP(247) = 0
      CALL MPI_REDUCE( KEEP(246), KEEP(247), 1, MPI_INTEGER, 
     &                 MPI_MAX, MASTER, id%COMM, IERR)
C
C     Reduce compression times: get max compression times
      CALL MPI_REDUCE( id%DKEEP(97), id%DKEEP(98), 1,
     &     MPI_DOUBLE_PRECISION,
     &     MPI_MAX, MASTER, id%COMM, IERR)
C
      CALL MPI_REDUCE( RINFO(2), RINFOG(2), 2,
     &                 MPI_DOUBLE_PRECISION,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      CALL MUMPS_REDUCEI8( id%KEEP8(31)+id%KEEP8(64),id%KEEP8(6),
     &                     MPI_SUM, MASTER, id%COMM )
C
      IF (id%MYID.EQ.0) THEN
C      In MegaBytes
       RINFOG(16) = dble(id%KEEP8(6)*int(KEEP(35),8))/dble(1D6)
       IF (KEEP(201).LE.0) THEN
        RINFOG(16) = ZERO
       ENDIF
      ENDIF
      CALL MUMPS_REDUCEI8( id%KEEP8(48),id%KEEP8(148), MPI_SUM,
     &                     MASTER, id%COMM )
      CALL MUMPS_SETI8TOI4(id%KEEP8(148), INFOG(9))
C
      CALL MPI_REDUCE( int(id%INFO(10),8), id%KEEP8(128),
     &                 1, MPI_INTEGER8,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      IF (id%MYID.EQ.MASTER) THEN
        CALL MUMPS_SETI8TOI4(id%KEEP8(128), id%INFOG(10))
      ENDIF
C     Use MPI_MAX for this one to get largest front size
      CALL MPI_ALLREDUCE( id%INFO(11), INFOG(11), 1, MPI_INTEGER,
     &                 MPI_MAX, id%COMM, IERR)
C     make maximum effective frontal size available on all procs
C     for solve phase
C     (Note that INFO(11) includes root size on root master)
      KEEP(133) = INFOG(11)
      CALL MPI_REDUCE( id%INFO(12), INFOG(12), 3, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      CALL MPI_REDUCE( id%INFO(40), INFOG(50), 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      CALL MPI_REDUCE( KEEP(103), INFOG(25), 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      KEEP(229) = INFOG(25)
      CALL MPI_REDUCE( KEEP(105), INFOG(25), 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      KEEP(230) = INFOG(25)
C
      id%INFO(25) = KEEP(98)
      CALL MPI_ALLREDUCE( id%INFO(25), INFOG(25), 1, MPI_INTEGER,
     &                 MPI_SUM, id%COMM, IERR)
C     Extra copies due to in-place stacking
      CALL MUMPS_REDUCEI8( id%KEEP8(8), id%KEEP8(108), MPI_SUM,
     &                     MASTER, id%COMM )
C     Entries in factors
      CALL MUMPS_SETI8TOI4(id%KEEP8(10), id%INFO(27))
      CALL MUMPS_REDUCEI8( id%KEEP8(10),id%KEEP8(110), MPI_SUM,
     &                     MASTER, id%COMM )
      CALL MUMPS_SETI8TOI4(id%KEEP8(110), INFOG(29))
C     Initialize INFO(28)/INFOG(35) in case BLR not activated
      id%INFO(28)  = id%INFO(27)
      INFOG(35)    = INFOG(29)
C     ==============================
C     LOW-RANK
C     ==============================
      IF ( KEEP(486) .NE. 0 ) THEN  !LR is activated
C       Compute and Save local amount of flops in case of BLR
        RINFO(4) = dble(FLOP_FRFRONTS + FLOP_FACTO_FR - FLOP_LRGAIN
     &             + FLOP_COMPRESS  + FLOP_FRFRONTS)
C
C       Compute and Save local number of entries in compressed factors
C 
        ITMP8 =  id%KEEP8(10) - int(MRY_LU_LRGAIN,8)
        CALL MUMPS_SETI8TOI4( ITMP8, id%INFO(28))
C
        CALL MPI_REDUCE( MRY_LU_LRGAIN, TMP_MRY_LU_LRGAIN
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
C       TMP_MRY_LU_LRGAIN_mixed need be initialized anyway
        TMP_MRY_LU_LRGAIN_mixed =  ZERO
        IF (id%KEEP(450).GT.0) THEN
          CALL MPI_REDUCE( MRY_LU_LRGAIN_mixed,
     &                   TMP_MRY_LU_LRGAIN_mixed
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        ENDIF
        TMP_MRY_CB_LRGAIN_mixed = ZERO
        IF (id%KEEP(452).GT.0) THEN
          CALL MPI_REDUCE( MRY_CB_LRGAIN_mixed,
     &                       TMP_MRY_CB_LRGAIN_mixed
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
        ENDIF
        CALL MPI_REDUCE( MRY_LU_FR, TMP_MRY_LU_FR
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( MRY_CB_FR, TMP_MRY_CB_FR
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( MRY_CB_LRGAIN, TMP_MRY_CB_LRGAIN
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( FLOP_LRGAIN, TMP_FLOP_LRGAIN
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( FLOP_TRSM_FR, TMP_FLOP_TRSM_FR
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( FLOP_TRSM_LR, TMP_FLOP_TRSM_LR
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( FLOP_UPDATE_FR, TMP_FLOP_UPDATE_FR
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( FLOP_UPDATE_LR, TMP_FLOP_UPDATE_LR
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( FLOP_FRSWAP_COMPRESS,
     &                   TMP_FLOP_FRSWAP_COMPRESS
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( FLOP_MIDBLK_COMPRESS,
     &                   TMP_FLOP_MIDBLK_COMPRESS
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( FLOP_UPDATE_LRLR3, TMP_FLOP_UPDATE_LRLR3
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE(FLOP_ACCUM_COMPRESS, TMP_FLOP_ACCUM_COMPRESS
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( FLOP_TRSM, TMP_FLOP_TRSM
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( FLOP_PANEL, TMP_FLOP_PANEL
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( FLOP_FRFRONTS, TMP_FLOP_FRFRONTS
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( FLOP_COMPRESS, TMP_FLOP_COMPRESS
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( FLOP_DECOMPRESS, TMP_FLOP_DECOMPRESS
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( FLOP_CB_COMPRESS, TMP_FLOP_CB_COMPRESS
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( FLOP_CB_DECOMPRESS,TMP_FLOP_CB_DECOMPRESS
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( FLOP_FACTO_FR, TMP_FLOP_FACTO_FR
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( CNT_NODES,TMP_CNT_NODES
     &                  , 1, MPI_INTEGER,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        IF (id%NPROCS.GT.1) THEN
          FLOP_FACTO_LR = FLOP_FACTO_FR - FLOP_LRGAIN
     &                      + FLOP_COMPRESS + FLOP_FRFRONTS
          CALL MPI_REDUCE( FLOP_FACTO_LR, AVG_FLOP_FACTO_LR
     &                    , 1, MPI_DOUBLE_PRECISION,
     &                     MPI_SUM, MASTER, id%COMM, IERR)
          IF (id%MYID.EQ.MASTER) THEN
            AVG_FLOP_FACTO_LR = AVG_FLOP_FACTO_LR/id%NPROCS
          ENDIF
          CALL MPI_REDUCE( FLOP_FACTO_LR, MIN_FLOP_FACTO_LR
     &                    , 1, MPI_DOUBLE_PRECISION,
     &                     MPI_MIN, MASTER, id%COMM, IERR)
          CALL MPI_REDUCE( FLOP_FACTO_LR, MAX_FLOP_FACTO_LR
     &                    , 1, MPI_DOUBLE_PRECISION,
     &                     MPI_MAX, MASTER, id%COMM, IERR)
        ENDIF ! NPROCS > 1
        CALL MPI_REDUCE( TIME_UPDATE, TMP_TIME_UPDATE
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_UPDATE_LRLR1, TMP_TIME_UPDATE_LRLR1
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_UPDATE_LRLR2, TMP_TIME_UPDATE_LRLR2
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_UPDATE_LRLR3, TMP_TIME_UPDATE_LRLR3
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_UPDATE_FRLR, TMP_TIME_UPDATE_FRLR
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_UPDATE_FRFR, TMP_TIME_UPDATE_FRFR
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_DIAGCOPY, TMP_TIME_DIAGCOPY
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_COMPRESS,TMP_TIME_COMPRESS
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_MIDBLK_COMPRESS, 
     &                   TMP_TIME_MIDBLK_COMPRESS
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_FRSWAP_COMPRESS,
     &                   TMP_TIME_FRSWAP_COMPRESS
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_CB_COMPRESS, TMP_TIME_CB_COMPRESS
     &                  , 1, MPI_DOUBLE_PRECISION, 
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_DECOMP, TMP_TIME_DECOMP
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_DECOMP_UCFS, TMP_TIME_DECOMP_UCFS
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_LRASM_NIV1, TMP_TIME_LRASM_NIV1
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE(TIME_LRASM_LOCASM2, TMP_TIME_LRASM_LOCASM2
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE(TIME_LRASM_MAPLIG1, TMP_TIME_LRASM_MAPLIG1
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE(TIME_LRASM_CONTRIB2, TMP_TIME_LRASM_CONTRIB2
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE(TIME_FRASM_LOCASM2, TMP_TIME_FRASM_LOCASM2
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE(TIME_FRASM_MAPLIG1, TMP_TIME_FRASM_MAPLIG1
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE(TIME_FRASM_CONTRIB2, TMP_TIME_FRASM_CONTRIB2
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_PANEL, TMP_TIME_PANEL
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_FAC_I, TMP_TIME_FAC_I
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_FAC_MQ, TMP_TIME_FAC_MQ
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_FAC_SQ, TMP_TIME_FAC_SQ
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_LRTRSM, TMP_TIME_LRTRSM
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_FRTRSM, TMP_TIME_FRTRSM
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_FRFRONTS, TMP_TIME_FRFRONTS
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        CALL MPI_REDUCE( TIME_LR_MODULE, TMP_TIME_LR_MODULE
     &                  , 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, MASTER, id%COMM, IERR)
        IF (id%MYID.EQ.MASTER) THEN
          IF (id%NPROCS.GT.1) THEN
C rename the stat variable so that COMPUTE_GLOBAL_GAINS can work for any
C number of procs
            MRY_LU_FR            = TMP_MRY_LU_FR
            MRY_LU_LRGAIN        = TMP_MRY_LU_LRGAIN
            MRY_LU_LRGAIN_mixed  = TMP_MRY_LU_LRGAIN_mixed
            MRY_CB_LRGAIN_mixed  = TMP_MRY_CB_LRGAIN_mixed
            MRY_CB_FR            = TMP_MRY_CB_FR
            MRY_CB_LRGAIN        = TMP_MRY_CB_LRGAIN
            FLOP_LRGAIN          = TMP_FLOP_LRGAIN
            FLOP_PANEL           = TMP_FLOP_PANEL
            FLOP_TRSM            = TMP_FLOP_TRSM
            FLOP_TRSM_FR         = TMP_FLOP_TRSM_FR
            FLOP_TRSM_LR         = TMP_FLOP_TRSM_LR
            FLOP_UPDATE_FR       = TMP_FLOP_UPDATE_FR
            FLOP_UPDATE_LR       = TMP_FLOP_UPDATE_LR
            FLOP_UPDATE_LRLR3    = TMP_FLOP_UPDATE_LRLR3
            FLOP_COMPRESS        = TMP_FLOP_COMPRESS
            FLOP_MIDBLK_COMPRESS = TMP_FLOP_MIDBLK_COMPRESS
            FLOP_FRSWAP_COMPRESS = TMP_FLOP_FRSWAP_COMPRESS
            FLOP_ACCUM_COMPRESS  = TMP_FLOP_ACCUM_COMPRESS
            FLOP_CB_COMPRESS     = TMP_FLOP_CB_COMPRESS
            FLOP_DECOMPRESS      = TMP_FLOP_DECOMPRESS
            FLOP_CB_DECOMPRESS   = TMP_FLOP_CB_DECOMPRESS
            FLOP_FRFRONTS        = TMP_FLOP_FRFRONTS
            FLOP_FACTO_FR        = TMP_FLOP_FACTO_FR
            CNT_NODES            = TMP_CNT_NODES
            TIME_UPDATE          = TMP_TIME_UPDATE         /id%NPROCS
            TIME_UPDATE_LRLR1    = TMP_TIME_UPDATE_LRLR1   /id%NPROCS
            TIME_UPDATE_LRLR2    = TMP_TIME_UPDATE_LRLR2   /id%NPROCS
            TIME_UPDATE_LRLR3    = TMP_TIME_UPDATE_LRLR3   /id%NPROCS
            TIME_UPDATE_FRLR     = TMP_TIME_UPDATE_FRLR    /id%NPROCS
            TIME_UPDATE_FRFR     = TMP_TIME_UPDATE_FRFR    /id%NPROCS
            TIME_COMPRESS        = TMP_TIME_COMPRESS       /id%NPROCS
            TIME_MIDBLK_COMPRESS = TMP_TIME_MIDBLK_COMPRESS/id%NPROCS
            TIME_FRSWAP_COMPRESS = TMP_TIME_FRSWAP_COMPRESS/id%NPROCS
            TIME_DIAGCOPY        = TMP_TIME_DIAGCOPY       /id%NPROCS
            TIME_CB_COMPRESS     = TMP_TIME_CB_COMPRESS    /id%NPROCS
            TIME_PANEL           = TMP_TIME_PANEL          /id%NPROCS
            TIME_FAC_I           = TMP_TIME_FAC_I          /id%NPROCS
            TIME_FAC_MQ          = TMP_TIME_FAC_MQ         /id%NPROCS
            TIME_FAC_SQ          = TMP_TIME_FAC_SQ         /id%NPROCS
            TIME_LRTRSM          = TMP_TIME_LRTRSM         /id%NPROCS
            TIME_FRTRSM          = TMP_TIME_FRTRSM         /id%NPROCS
            TIME_FRFRONTS        = TMP_TIME_FRFRONTS       /id%NPROCS
            TIME_LR_MODULE       = TMP_TIME_LR_MODULE      /id%NPROCS
            TIME_DECOMP          = TMP_TIME_DECOMP         /id%NPROCS
            TIME_DECOMP_UCFS     = TMP_TIME_DECOMP_UCFS    /id%NPROCS
            TIME_LRASM_NIV1      = TMP_TIME_LRASM_NIV1     /id%NPROCS
            TIME_LRASM_LOCASM2   = TMP_TIME_LRASM_LOCASM2  /id%NPROCS
            TIME_LRASM_MAPLIG1   = TMP_TIME_LRASM_MAPLIG1  /id%NPROCS
            TIME_LRASM_CONTRIB2  = TMP_TIME_LRASM_CONTRIB2 /id%NPROCS
            TIME_FRASM_LOCASM2   = TMP_TIME_FRASM_LOCASM2  /id%NPROCS
            TIME_FRASM_MAPLIG1   = TMP_TIME_FRASM_MAPLIG1  /id%NPROCS
            TIME_FRASM_CONTRIB2  = TMP_TIME_FRASM_CONTRIB2 /id%NPROCS
          ENDIF
          CALL COMPUTE_GLOBAL_GAINS(id%KEEP8(110),
     &         id%RINFOG(3),
     &         id%KEEP8(49), PROKG, MPG)
C         Number of entries in factor  INFOG(35) in 
C         compressed form is updated as long as 
C         BLR is activated,  this independently of the 
C         fact that factors are saved in LR.
          CALL MUMPS_SETI8TOI4(id%KEEP8(49),  id%INFOG(35))
          FRONTWISE = 0
C         WRITE gains also compute stats stored in DKEEP array
          DO I=1,LR_TABSIZE
            LR_TAB(I)  = id%DKEEP(I+LR_DKEEPSHIFT)
            LR_EPSILON = id%DKEEP(8)
          ENDDO
          CALL SAVEandWRITE_GAINS(FRONTWISE, KEEP(489),
     &                LR_DKEEPSHIFT, LR_TABSIZE, LR_TAB, LR_EPSILON,
     &                N,  id%ICNTL(36),
     &                KEEP(487), KEEP(488), KEEP(490),
     &                KEEP(491), KEEP(50), KEEP(486), 
     &                KEEP(450), KEEP(452), NPREC_STOC,
     &                KEEP(249)*max(KEEP(381), 1),
     &                KEEP(472), KEEP(475), KEEP(478), 
     &                KEEP(480), KEEP(481),
     &                KEEP(483), KEEP(484), KEEP(117), 
     &                id%KEEP8(110), id%KEEP8(49),
     &                KEEP(28), id%NPROCS, MPG, PROKG)
C           flops when BLR activated
            RINFOG(14) = id%DKEEP(56)
          DO I=1,18
            id%DKEEP(I+LR_DKEEPSHIFT)=dble(LR_TAB(I))
          ENDDO
        ELSE
            RINFOG(14) = 0.0D00
        ENDIF
        IF (id%MYID .eq. MASTER) THEN
          KEEP(399) = KEEP399_SAVE
        ENDIF
        IF (id%MYID .eq. MASTER) THEN
C         Reset KEEP(451,452,483,484) to value on entry        
C         (only on master but it will be broadcasted again at nex call)
          KEEP(451) = KEEP451_SAVE
          KEEP(452) = KEEP452_SAVE
          KEEP(483) = KEEP483_SAVE
          KEEP(484) = KEEP484_SAVE
        ENDIF
      ENDIF
C     ==============================
C     NULL PIVOTS AND RANK-REVEALING
C     ==============================
      IF (KEEP(19).EQ.0.AND.KEEP(53).NE.0) THEN 
C       restore KEEP(20)
        KEEP(20) = KEEP20_SAVE
      ENDIF
      IF(KEEP(110) .EQ. 1) THEN
C        -- make available to users the local number of null pivots detected 
C        -- with ICNTL(24) = 1.
         id%INFO(18) = KEEP(109)
         CALL MPI_ALLREDUCE( KEEP(109), KEEP(112), 1, MPI_INTEGER,
     &        MPI_SUM, id%COMM, IERR)
      ELSE
         id%INFO(18)  = 0
         KEEP(109) = 0
         KEEP(112) = 0
      ENDIF
      IF ( associated( id%PIVNUL_LIST) ) THEN
        DEALLOCATE(id%PIVNUL_LIST)
        NULLIFY(id%PIVNUL_LIST)
      ENDIF
C     set INFOG(28) even in case of error 
      IF (id%MYID.EQ.MASTER) THEN
C       INFOG(28) deficiency resulting from ICNTL(24) and ICNTL(56)
C       save in INFOG(51) deficiency found by ICNTL(24)
        INFOG(51)=KEEP(112)
        INFOG(28)=KEEP(112)
        IF (KEEP(17).GT.0) THEN
          INFOG(28)=KEEP(112)+KEEP(17)
        ENDIF
      ENDIF
C
      IF (id%INFO(1).GE.0) THEN
C{     PIVNUL_LIST not meaningful in case of error
C      (do not allocate)
       IF (id%MYID.EQ.MASTER) THEN
        IF ( INFOG(28) .GT. 0 ) THEN
          ALLOCATE(id%PIVNUL_LIST(INFOG(28)), stat=IERR)
          IF ( IERR .GT. 0 ) THEN
            id%INFO(1)=-13
            id%INFO(2)=INFOG(28)
          END IF
        ENDIF
       ELSE
C       id%PIVNUL_LIST(1:KEEP(109)) used during sol_driver on slaves 
C       to initialize id%RHSINTR
        IF (KEEP(109).GT.0) THEN
          ALLOCATE(id%PIVNUL_LIST(KEEP(109)), stat=IERR)
          IF ( IERR .GT. 0 ) THEN
           id%INFO(1)=-13
           id%INFO(2)=INFOG(28)
          END IF
        ENDIF
       ENDIF
C}
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO(1).lt.0 ) GOTO 490
      IF ( (KEEP(19).NE.0) .AND. (KEEP(143) .NE. KEEP(17)) ) THEN
C
C         Raise a warning (on all MPI processes) since determinant or
C         inertia or null pivot list is not consistent with deficiency
C         computed with ICNTL(46)=1.
C
        id%INFO(1) = id%INFO(1)+16
        id%INFO(2) = KEEP(112)+KEEP(143)
        IF (KEEP(118) .GE. 40) THEN
          IF ( PROKG ) THEN
            WRITE(MPG,'(/A,A/,A,A,I8/,A,A,I8/)')
     &      " WARNING: in the context of rank-revealing,",
     &      " the inertia, determinant and pivnul list",  
     &      " are computed with RR (rank-revealing)-LU,", 
     &      " but the deficiency found by RR-LU: ", 
     &      id%INFO(2),
     &      " is different from the deficiency computed",
     &      " with ICNTL(56)>0: ",  KEEP(112)+KEEP(17)
          ENDIF
        ELSE
          IF ( LP .GT. 0 ) THEN
            WRITE(LP,'(/A,A/,A/)')
     &      " ERROR : in the context of rank-revealing,",
     &      " the inertia, determinant and pivnul list",  
     &      " are not correct because RR LU not called "
          ENDIF
        ENDIF
      ENDIF
C     ========================================
C     We now provide to the host the part of
C     PIVNUL_LIST resulting from the processing
C     of the root node and we update id%INFO(18)
C     on the processor holding the root to
C     include null pivots relative to the root
C     ========================================
      IF ( KEEP(109).GT.0 ) THEN
          DO I=1, KEEP(109) 
            id%PIVNUL_LIST(I)= 
     &             PIVNUL_LIST_STRUCT%PIVNUL_LIST(I)
          ENDDO
      ENDIF
      IF (KEEP(17) .NE. 0) THEN
        IF (id%MYID .EQ. ID_ROOT) THEN
C         Include in id%INFO(18) null pivots resulting
C         from deficiency on the root. In this way,
C         the sum of all id%INFO(18) is equal to INFOG(28).
          id%INFO(18)=id%INFO(18)+KEEP(17)
        ENDIF
        IF (ID_ROOT .EQ. MASTER) THEN
          IF (id%MYID.EQ.MASTER) THEN
C           --------------------------------------------------
C           Null pivots of root have been stored in
C           PIVNUL_LIST_STRUCT%PIVNUL_LIST(
C                         KEEP(109)+1:KEEP(109)+KEEP(17) )
C           Shift them at the end of the list because:
C           * this is what we need to build the null space
C           * we would otherwise overwrite them on the host
C             when gathering null pivots from other processors
C           --------------------------------------------------
            DO I= KEEP(17), 1, -1
              id%PIVNUL_LIST(KEEP(112)+I)= 
     &               PIVNUL_LIST_STRUCT%PIVNUL_LIST(KEEP(109)+I)
            ENDDO
          ENDIF
        ELSE
C         ---------------------------------
C         Null pivots of root must be sent
C         from the processor responsible of
C         the root to the host (or MASTER).
C         ---------------------------------
          IF (id%MYID .EQ. ID_ROOT) THEN
            CALL MPI_SEND(
     &               PIVNUL_LIST_STRUCT%PIVNUL_LIST(KEEP(109)+1), 
     &               KEEP(17), MPI_INTEGER, MASTER, ZERO_PIV,
     &               id%COMM, IERR)
          ELSE IF (id%MYID .EQ. MASTER) THEN
            CALL MPI_RECV(id%PIVNUL_LIST(KEEP(112)+1), KEEP(17),
     &               MPI_INTEGER, ID_ROOT, ZERO_PIV,
     &               id%COMM, STATUS, IERR )
          ENDIF
        ENDIF
      ENDIF
C     ===========================
C     gather zero pivots indices
C     on the host node
C     ===========================
C     In case of non working host, the following code also
C     works considering that KEEP(109) is equal to 0 on
C     the non-working host
      IF(KEEP(110) .EQ. 1) THEN
         ALLOCATE(ITMP2(id%NPROCS),stat = IERR )  ! deallocated in 490
         IF ( IERR .GT. 0 ) THEN
            id%INFO(1)=-13
            id%INFO(2)=id%NPROCS
         END IF
         CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
         IF (id%INFO(1).LT.0) GOTO 490
         CALL MPI_GATHER ( KEEP(109),1, MPI_INTEGER, 
     &        ITMP2(1), 1, MPI_INTEGER, 
     &        MASTER, id%COMM, IERR)
         IF(id%MYID .EQ. MASTER) THEN
            POSBUF = ITMP2(1)+1
C           First null pivot of master is in
C           position 1 of global list
            KEEP(220)=1
            DO I = 1,id%NPROCS-1
              IF (ITMP2(I+1).GT.0) THEN
                CALL MPI_RECV(id%PIVNUL_LIST(POSBUF), ITMP2(I+1), 
     &              MPI_INTEGER,I, 
     &              ZERO_PIV, id%COMM, STATUS, IERR)
C              Send position POSBUF of first null pivot of proc I
C              in global list. Will allow to quickly identify during
C              the solve step if one is concerned by a global position
C              K, 0 <= K <= INFOG(28).
               CALL MPI_SEND(POSBUF, 1, MPI_INTEGER, I, ZERO_PIV,
     &              id%COMM, IERR)
               POSBUF = POSBUF + ITMP2(I+1)
              ENDIF
            ENDDO
         ELSE
           IF (KEEP(109).GT.0) THEN
            CALL MPI_SEND( 
     &           PIVNUL_LIST_STRUCT%PIVNUL_LIST(1), KEEP(109), 
     &           MPI_INTEGER,
     &           MASTER,ZERO_PIV, id%COMM, IERR)
            CALL MPI_RECV( KEEP(220), 1, MPI_INTEGER, MASTER, ZERO_PIV,
     &           id%COMM, STATUS, IERR )
           ENDIF
         ENDIF
      ENDIF
      IF (associated( PIVNUL_LIST_STRUCT%PIVNUL_LIST)) THEN
        DEALLOCATE(PIVNUL_LIST_STRUCT%PIVNUL_LIST)
        NULLIFY(PIVNUL_LIST_STRUCT%PIVNUL_LIST)
      ENDIF
C     =====================================
C     Statistics concerning the determinant
C     =====================================
C
C     1/ on the host better take into account null pivots if scaling:
C
C     Since null pivots are excluded from the computation
C     of the determinant, we also exclude the corresponding
C     scaling entries. Since those entries have already been
C     taken into account before the factorization, we multiply
C     the determinant on the host by the scaling values corresponding
C     to pivots in PIVNUL_LIST.
      IF (id%MYID.EQ.MASTER .AND. LSCAL. AND. KEEP(258).NE.0) THEN
        K = min(KEEP(143), KEEP(17))
        K = max(K, 0)
        DO I = 1, KEEP(112)+ K
c       DO I = 1, id%INFOG(28)  ! all null pivots + singular values
          CALL ZMUMPS_UPDATEDETER_SCALING(
     &                            id%ROWSCA(id%PIVNUL_LIST(I)),
     &                            id%DKEEP(6), KEEP(259))
          CALL ZMUMPS_UPDATEDETER_SCALING(
     &                            id%COLSCA(id%PIVNUL_LIST(I)),
     &                            id%DKEEP(6), KEEP(259))
        ENDDO
      ENDIF
C
C     2/ Swap signs depending on pivoting on each proc
C
      IF (KEEP(258).NE.0) THEN
C       Return the determinant in INFOG(34) and RINFOG(12/13)
        IF (KEEP(260).EQ.-1) THEN ! Local to each processor
          id%DKEEP(6)=-id%DKEEP(6)
          id%DKEEP(7)=-id%DKEEP(7)
        ENDIF
C
C       3/ Perform a reduction
C
        CALL ZMUMPS_DETER_REDUCTION(
     &           id%COMM, id%DKEEP(6), KEEP(259),
     &           RINFOG(12), INFOG(34), id%NPROCS)
C
C       4/ Swap sign if needed
C
        IF (id%KEEP(50).EQ.0 .AND. id%MYID.EQ. MASTER) THEN
C         Modify sign of determinant according
C         to unsymmetric permutation (max-trans
C         of max-weighted matching)
          IF (id%KEEP(23).NE.0) THEN
            CALL ZMUMPS_DETER_SIGN_PERM(
     &           RINFOG(12), id%N,
     &           id%UNS_PERM(1) )
C           Remark that RINFOG(12/13) are modified only
C           on the host but will be broadcast on exit
C           from MUMPS (see ZMUMPS_DRIVER)
          ENDIF
        ENDIF
      ENDIF
 490  IF (allocated(ITMP2)) DEALLOCATE(ITMP2)
C
C     =====================================
C     Statistics relative to min/max pivots
C     =====================================
      CALL MPI_REDUCE( id%DKEEP(19), RINFOG(19), 1, 
     &                 MPI_DOUBLE_PRECISION,
     &                 MPI_MIN, MASTER, id%COMM, IERR )
      CALL MPI_REDUCE( id%DKEEP(20), RINFOG(20), 1, 
     &                 MPI_DOUBLE_PRECISION,
     &                 MPI_MIN, MASTER, id%COMM, IERR )
      CALL MPI_REDUCE( id%DKEEP(21), RINFOG(21), 1, 
     &                 MPI_DOUBLE_PRECISION,
     &                 MPI_MAX, MASTER, id%COMM, IERR )
C     =========================================
C     Centralized number of swaps for pivoting
C     =========================================
      CALL MPI_REDUCE( id%KEEP8(80), ITEMP8, 1, MPI_INTEGER8,
     &                 MPI_SUM, MASTER, id%COMM, IERR )
      IF (id%MYID .EQ. MASTER) THEN
        CALL MUMPS_SETI8TOI4(ITEMP8,id%INFOG(48))
      ENDIF
C     ==========================================
C     Centralized largest increase of panel size
C     ==========================================
      CALL MPI_REDUCE( id%KEEP(425), id%INFOG(49), 1, MPI_INTEGER,
     &                 MPI_MAX, MASTER, id%COMM, IERR )
      IF ( PROKG ) THEN
C{    -----------------------------
C     PRINT STATISTICS  (on master)
C     -----------------------------
          WRITE(MPG,99984) RINFOG(2),RINFOG(3),KEEP(52), 
     &                    id%KEEP8(148),
     &                    id%KEEP8(128), INFOG(11), id%KEEP8(110)
          IF (id%KEEP(50) == 0) THEN
            ! off diag pivots
            WRITE(MPG, 99985) INFOG(12)
          END IF
          IF (id%KEEP(50) .NE. 1) THEN
            ! delayed pivots
            WRITE(MPG, 99982) INFOG(13)
          END IF
          IF (KEEP(97) .NE. 0) THEN
            ! tiny pivots
            WRITE(MPG, '(A,D16.4)') 
     &     ' Effective static pivoting thresh., CNTL(4) =', SEUIL
            WRITE(MPG, 99986) INFOG(25)
          ENDIF
          IF (id%KEEP(50) == 2) THEN
            !number of 2x2 pivots in type 1 nodes
             WRITE(MPG, 99988) KEEP(229)
            !number of 2x2 pivots in type 2 nodes
             WRITE(MPG, 99989) KEEP(230)
          ENDIF
          !Total deficiency
          IF (KEEP(110).NE.0.OR.KEEP(19).NE.0)
     &    WRITE(MPG, 99992) KEEP(17)+KEEP(112) ! = INFOG(28)
          !number of zero pivots
          IF (KEEP(110) .NE.0) THEN
              WRITE(MPG, 99991) INFOG(51)  ! KEEP(112)
          ENDIF
          !Deficiency on root
          IF ( KEEP(19) .ne. 0 )
c         IF ( KEEP(17) .ne. 0 )
     &      WRITE(MPG, 99983) KEEP(17)
          !Smallest pivot with also null pivots in abs value
          WRITE(MPG, 99995) RINFOG(19)
          !Smallest pivot in abs value
          WRITE(MPG, 99993) RINFOG(20)
          !Largest pivot in abs value
          WRITE(MPG, 99994) RINFOG(21)
          !value of ICNTL(12) that was effectively used.
          WRITE(MPG, 99996) INFOG(24)
          ! Memory compress
          WRITE(MPG, 99981) INFOG(14)
          ! Extra copies due to ip stack in unsym case
          ! in core case (or OLD_OOC_PANEL)
          IF (id%KEEP8(108) .GT. 0_8) THEN
            WRITE(MPG, 99980) id%KEEP8(108)
          ENDIF
          IF  ((KEEP(60).NE.0) .AND. INFOG(25).GT.0) THEN
          !  Schur on and tiny pivots set in last level 
          ! before the Schur if KEEP(114)=0
           WRITE(MPG, '(A)') 
     & " ** Warning Static pivoting was necessary"
           WRITE(MPG, '(A)') 
     & " ** to factor interior variables with Schur ON"
          ENDIF
          IF (KEEP(258).NE.0) THEN
            WRITE(MPG,99978) RINFOG(12)
            WRITE(MPG,99979) RINFOG(13)
            WRITE(MPG,99977) INFOG(34)
          ENDIF
C}
      END IF
* ==========================================
*
*  End of Factorization Phase
*
* ==========================================
C
C  Goto 500 is done when
C  LOAD_INIT
C  OOC_INIT_FACTO
C  MUMPS_FDM_INIT
#if ! defined(NO_FDM_DESCBAND)
C  MUMPS_FDBD_INIT
#endif
#if ! defined(NO_FDM_MAPROW)
C  MUMPS_FMRD_INIT
#endif
C  are all called.
C
 500  CONTINUE
C     Redo free INTARR and DBLARR in case an error occurred
C     after allocating them and before freeing them.
      IF (associated(DBLARR)) THEN
        DEALLOCATE(DBLARR)
        NULLIFY(DBLARR)
      ENDIF
      IF (allocated(INTARR)) THEN
        DEALLOCATE(INTARR)
      ENDIF
      IF ( KEEP(38) .NE. 0) THEN
        IF (associated(idintr%root%RG2L)) THEN
          DEALLOCATE(idintr%root%RG2L)
          NULLIFY(idintr%root%RG2L)
        ENDIF
      ENDIF
#if defined(MPI_TO_K_OMP)
      IF (allocated(CORES_LIST)) THEN
        deallocate(CORES_LIST)
      ENDIF
#endif
#if ! defined(NO_FDM_DESCBAND)
      IF (I_AM_SLAVE) THEN
        CALL MUMPS_FDBD_END(id%INFO(1))  ! INFO(1): input only
      ENDIF
#endif
#if ! defined(NO_FDM_MAPROW)
      IF (I_AM_SLAVE) THEN
        CALL MUMPS_FMRD_END(id%INFO(1))  ! INFO(1): input only
      ENDIF
#endif
      IF (I_AM_SLAVE) THEN
C{      Terminate BLR module except if it is still needed for solve
        IF ( 
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &      .AND. id%INFO(1).GE.0
     &     ) THEN
C         Store pointer to BLR_ARRAY in MUMPS structure
C         (requires successful factorization otherwise module is freed)
#if defined(SINGLEINDOUBLE)
          IF (KEEP(148) .EQ. 0) THEN
#endif
          CALL ZMUMPS_BLR_MOD_TO_STRUC(id%BLRARRAY_ENCODING)
#if defined(SINGLEINDOUBLE)
         ELSE IF(KEEP(148).EQ.2) THEN
           CALL ZMUMPS_LPid_TO_Cid
     &             (id, idintr,KEEP(34), id%INFO(1), id%INFO(2))
           CALL ZMUMPS_BLR_MOD_TO_STRUC(id%BLRARRAY_ENCODING)
         ELSE
C          Same id%BLRARRAY_ENCODING for single or double prec BLR
          CALL CMUMPS_BLR_MOD_TO_STRUC(id%BLRARRAY_ENCODING)
         ENDIF
#endif
        ELSE
C         INFO(1) positive or negative
#if defined(SINGLEINDOUBLE)
          IF (KEEP(148) .EQ. 0) THEN
#endif
          CALL ZMUMPS_BLR_END_MODULE(id%INFO(1), id%KEEP8, id%KEEP(34))
#if defined(SINGLEINDOUBLE)
         ELSE IF (KEEP(148).EQ.2) THEN
          CALL ZMUMPS_BLR_END_MODULE(id%INFO(1), id%KEEP8, id%KEEP(34))
!          CALL CMUMPS_BLR_END_MODULE(id%INFO(1), id%KEEP8, id%KEEP(34))
         ELSE
          CALL CMUMPS_BLR_END_MODULE(id%INFO(1), id%KEEP8, id%KEEP(34))
         ENDIF
#endif
        ENDIF
C}
      ENDIF
      CALL MUMPS_PROPINFO(ICNTL(1), id%INFO(1),
     &          id%COMM, id%MYID)
      IF (I_AM_SLAVE) THEN
        CALL MUMPS_FDM_END('A')
C       Terminate BLR module except if it is still needed for solve
        IF ( 
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &      .AND. id%INFO(1).GE.0
     &     ) THEN
           CALL MUMPS_FDM_MOD_TO_STRUC('F', id%FDM_F_ENCODING,
     &        id%INFO(1))
           IF (.NOT. associated(id%FDM_F_ENCODING)) THEN
             WRITE(*,*) "Internal error 2 in ZMUMPS_FAC_DRIVER"
           ENDIF
        ELSE
           CALL MUMPS_FDM_END('F')
        ENDIF
      ENDIF
C
C  Goto 514 is done when an
C  error occurred in MUMPS_FDM_INIT
C  or (after FDM_INIT but before
C  OOC_INIT)
C
 514  CONTINUE
      IF ( I_AM_SLAVE ) THEN
         IF ((KEEP(201).EQ.1).OR.(KEEP(201).EQ.2)) THEN
            CALL ZMUMPS_OOC_END_FACTO(id%KEEP,id%KEEP8,
     &      id%OOC_MAX_NB_NODES_FOR_ZONE,id%OOC_TOTAL_NB_NODES,
     &      id%OOC_FILE_NAMES, id%INFO, id%OOC_FILE_NAME_LENGTH,
     &      id%OOC_NB_FILES, IERR)
            IF (id%ASSOCIATED_OOC_FILES) THEN
              id%ASSOCIATED_OOC_FILES = .FALSE.
            ENDIF
            IF (IERR.LT.0 .AND. id%INFO(1) .GE. 0) id%INFO(1) = IERR
         ENDIF
         IF (WK_USER_PROVIDED) THEN
C     at the end of a phase S is always null when WK_USER provided
            NULLIFY(id%S)
         ELSE IF (KEEP(201).NE.0) THEN
C           ----------------------------------------
C           In OOC or if KEEP(201).EQ.-1 we always
C           free S at end of factorization. As id%S
C           may be unassociated in case of error
C           during or before the allocation of id%S,
C           we only free S when it was associated.
C           ----------------------------------------
            IF (associated(id%S))  THEN
              CALL ZMUMPS_DM_FREE_S_WK(id%S, KEEP(430))
C             Reset KEEP(430)=0 since S will be allocated
C             from Fortran during solve
              KEEP(430) = 0
            ENDIF
            NULLIFY(id%S)   ! in all cases
#if defined(SINGLEINDOUBLE)
            IF (associated(id%LPS)) THEN
              CALL CMUMPS_DM_FREE_S_WK(id%LPS, KEEP(430))
            ENDIF
            NULLIFY(id%LPS)
#endif
            id%KEEP8(23)=0_8
         ELSE ! in core
#if defined(SINGLEINDOUBLE)
            IF (KEEP(148) .EQ.0) THEN
#endif
              CALL ZMUMPS_TRY_COMPACT_FACTORS(ICNTL49_LOC,
     &        WK_USER_PROVIDED, id%S, id%KEEP, id%KEEP8,
     &        id%INFO, id%MYID, id%ICNTL, PROK, MP,
     &        ZMUMPS_LBUFR_BYTES8, ZMUMPS_LBUF8,
     &        LIWK, LIWK8 )
#if defined(SINGLEINDOUBLE)
            ELSE
              CALL CMUMPS_TRY_COMPACT_FACTORS(ICNTL49_LOC,
     &        WK_USER_PROVIDED, id%LPS, id%KEEP, id%KEEP8,
     &        id%INFO, id%MYID, id%ICNTL, PROK, MP,
     &        ZMUMPS_LBUFR_BYTES8, ZMUMPS_LBUF8, 
     &        LIWK, LIWK8 )
            ENDIF
#endif
         ENDIF
      ELSE  ! host not working
         IF (WK_USER_PROVIDED) THEN
C     at the end of a phase S is always freed when WK_USER provided
            NULLIFY(id%S)
         ELSE
            IF (associated(id%S))  DEALLOCATE(id%S)
            NULLIFY(id%S)
#if defined(SINGLEINDOUBLE)
            IF (associated(id%LPS)) DEALLOCATE(id%LPS)
            NULLIFY(id%LPS)
#endif
            id%KEEP8(23)=0_8
         END IF
      END IF
C
C     Goto 513 is done in case of error where LOAD_INIT was
C     called but not the scaling nor OOC_INIT_FACTO.
 513  CONTINUE
      IF ( I_AM_SLAVE ) THEN
         CALL MUMPS_LOAD_END( id%INFO(1), id%NSLAVES, IERR )
         IF (IERR.LT.0 .AND. id%INFO(1) .GE. 0) id%INFO(1) = IERR
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
C
C     Goto 516 is done in case of error when GPU initialiwqtion
C     has been performed and scaling was optionally computed but
C     not LOAD_INIT nor OOC_INIT_FACTO. We can then extract
C     scaling arrays in case of error.
 516  CONTINUE
C     --------------------------------------------
C     We now build id%ROWSCA_loc and id%COLSCA_loc
C     in case of successful factorization, in the
C     numbering associated to the fully summed
C     variables of the frontal matrices.
C     This requires the factorization to be
C     successful because otherwise we do not have
C     the final lists of pivots associated to
C     the fronts, including delayed pivots and
C     symmetric/unsymmetric permutations done
C     during the factorization process.
C     --------------------------------------------
      IF (LSCAL .AND. id%INFO(1).GE.0) THEN
        CALL ZMUMPS_EXTRACT_SCALING(id)
C       occurs during scaling extraction, keep the error.
        IF ( id%INFO(1) .LT. 0 ) THEN
          IF (LPOK) THEN
            WRITE(LP,*) "Facto OK but error during EXTRACT_SCALING",
     &      id%INFO(1:2)
          ENDIF
        ENDIF
      ENDIF
C
C     Goto 517 is done when an error occurs when GPU initialization
C     has been performed but not LOAD_INIT or OOC_INIT_FACTO, e.g.
C     when an error occurred during the scaling.
 517  CONTINUE
      IF (associated( PIVNUL_LIST_STRUCT%PIVNUL_LIST)) THEN
        DEALLOCATE(PIVNUL_LIST_STRUCT%PIVNUL_LIST)
        NULLIFY(PIVNUL_LIST_STRUCT%PIVNUL_LIST)
      ENDIF
#if defined(USE_GPU)
      IF (id%KEEP(417).GT.0) then
#if defined(SINGLEINDOUBLE)
        IF (KEEP(148) .EQ. 0) THEN
#endif
        CALL ZMUMPS_END_GPU_C(id%KEEP(400),id%KEEP(422)
     &  , id%RINFO(9), id%RINFO(10)
     &  )
#if defined(SINGLEINDOUBLE)
        ELSE
        CALL CMUMPS_END_GPU_C(id%KEEP(400),id%KEEP(422)
     &  , LPRINFO9, LPRINFO10
     &  )
        id%RINFO(9) =dble(LPRINFO9)
        id%RINFO(10)=dble(LPRINFO10)
        ENDIF
#endif
C       Sum of flops returned in RINFOG(22)
        CALL MPI_REDUCE( id%RINFO(9), id%RINFOG(22), 1,
     &                 MPI_DOUBLE_PRECISION,
     &                 MPI_SUM, MASTER, id%COMM, IERR )
C       Average GPU time returned in RINFOG(23)
        CALL MPI_REDUCE( id%RINFO(10), id%RINFOG(23), 1,
     &                 MPI_DOUBLE_PRECISION,
     &                 MPI_SUM, MASTER, id%COMM, IERR )
        id%RINFOG(23) = id%RINFOG(23) / id%NSLAVES
        IF (PROKG) THEN
          WRITE(MPG,99975) id%RINFOG(22)
          WRITE(MPG,99976) id%RINFOG(23)
        ENDIF
      ENDIF
#endif
C
C     Goto 530 is done when an error occurs before
C     the calls to GPU_INIT, LOAD_INIT and OOC_INIT_FACTO
 530  CONTINUE
C  Fwd in facto: free RHS_MUMPS in case
C  it was allocated.
      IF (RHS_MUMPS_ALLOCATED) DEALLOCATE(RHS_MUMPS)
      NULLIFY(RHS_MUMPS)
      RHS_MUMPS_ALLOCATED = .FALSE.
C
      id%KEEP8(26) = KEEP826_SAVE
      RETURN
 120  FORMAT(/' Local redistrib: data local/sent           =',I16,I16)
 125  FORMAT(/' Redistrib: total data local/sent           =',I16,I16)
 130  FORMAT(//'****** FACTORIZATION STEP ********'/)
 140  FORMAT(/' Statistics on the scaling phase'
     & /' Elapsed time for scaling                   =',F12.4)
 160  FORMAT(
     & ' Elapsed time to reformat/distribute matrix =',F12.4/)
 166  FORMAT(' Max difference from 1 after scaling the entries',
     &       ' for ONE-NORM (option 7/8)   =',D9.2/)
 170  FORMAT(' STATISTICS PRIOR NUMERICAL FACTORIZATION ...'/
     &        ' Size of internal working array S           =',I16/
     &        ' Size of internal working array IS          =',I16/
     &        ' Minimum (ICNTL(14)=0) size of S            =',I16/
     &        ' Minimum (ICNTL(14)=0) size of IS           =',I16/
     &        ' Real space for original matrix             =',I16/
     &        ' Integer space for original matrix          =',I16/
     &        ' INFO(3) Real space for factors (estimated) =',I16/ 
     &        ' INFO(4) Integer space for factors (estim.) =',I16/
     &        ' Maximum frontal size (estimated)           =',I16)
 172  FORMAT(' GLOBAL STATISTICS PRIOR NUMERICAL FACTORIZATION ...'/
     &        ' Number of working processes                =',I16/
     &        ' ICNTL(22) Out-of-core option               =',I16/
     &        ' ICNTL(35) BLR activation (eff. choice)     =',I16/
     &        ' ICNTL(37) BLR CB compression (eff. choice) =',I16/
     &        ' ICNTL(40) Mixed/adaptive precision BLR     =',I16/
#if defined(SINGLEINDOUBLE)
     &        ' ICNTL(47) Single prec fact double instance =',I16/
#endif
     &        ' ICNTL(49) Compact workarray S (end facto.) =',I16/
     &        ' ICNTL(56) Effective value during facto.    =',I16/
     &        ' ICNTL(14) Memory relaxation                =',I16/
     &        ' INFOG(3) Real space for factors (estimated)=',I16/
     &        ' INFOG(4) Integer space for factors (estim.)=',I16/
     &        ' Maximum frontal size (estimated)           =',I16/
     &        ' Number of nodes in the tree                =',I16/
     &        ' ICNTL(23) Memory allowed (value on host)   =',I16/
     &        '           Sum over all procs               =',I16/
     &        ' Memory provided by user, sum of LWK_USER   =',I16/
     &        ' Effective threshold for pivoting, CNTL(1)  =',D16.4)
 173  FORMAT( ' Perform forward during facto, NRHS         =',I16)
 176  FORMAT( ' Relaxed pivoting (0=OFF)                   =',I16)
 180  FORMAT(/' Elapsed time for factorization                     =',
     &          F12.4)
 185  FORMAT(/' Elapsed time for (failed) factorization            =',
     & F12.4)
 186  FORMAT(/' Elapsed time to process root node                  =',
     & F12.4)
 187  FORMAT( ' Elapsed time under L0                      =',F12.4)
 188  FORMAT( ' Elapsed time under L0 (avg/max across MPI) =',
     &          F12.4,F12.4)
 189  FORMAT(/' Flops under L0 layer                       =',1PD12.3)
 190  FORMAT(/' Flops under L0 layer  (avg/max across MPI) =',
     &          1PD12.3,1PD12.3)
#if defined(USE_GPU)
99975 FORMAT( ' RINFOG(22) Total flops on GPU                      =',
     &          1PD12.3)
99976 FORMAT( ' RINFOG(23) Time spent on GPU (avg per MPI)         =',
     &          F12.4)
#endif
99977 FORMAT( ' INFOG(34)  Determinant (base 2 exponent)           =',
     &          I16)
99978 FORMAT( ' RINFOG(12) Determinant (real part)                 =',
     &          F16.8)
99979 FORMAT( ' RINFOG(12) Determinant (imaginary part)            =',
     &          F16.8)
99980 FORMAT( ' Extra copies due to In-Place stacking              =',
     &          I16)
99981 FORMAT( ' INFOG (14) Number of memory compress               =',
     &          I16)
99982 FORMAT( ' INFOG (13) Number of delayed pivots                =',
     &      I16)
99983 FORMAT( ' Nb of singularities detected by ICNTL(56)          =',
     &      I16)
99991 FORMAT( ' INFOG (51) Nb of null pivots detected by ICNTL(24) =',
     &      I16)
99992 FORMAT( ' INFOG (28) Estimated deficiency                    =',
     &      I16)
99995 FORMAT( ' RINFOG(19) Smallest pivot WITH perturbed pivots    =',
     &     1PD10.3)
99993 FORMAT( ' RINFOG(20) Smallest pivot WITHOUT perturbed pivots =',
     &     1PD10.3)
99994 FORMAT( ' RINFOG(21) Largest pivot in absolute value         =',
     &     1PD10.3)
99996 FORMAT( ' INFOG (24) Effective value of ICNTL(12)            =',
     &     I16)
99984 FORMAT(/'Leaving factorization with ...'/
     &        ' RINFOG (2) Operations in node assembly             =',
     &          1PD10.3/                                   
     &        ' ------ (3) Operations in node elimination          =',
     &          1PD10.3/                                       
     &        ' ICNTL  (8) Scaling effectively used                =',
     &          I16/                                           
     &        ' INFOG  (9) Real space for factors                  =',
     &          I16/                                           
     &        ' INFOG (10) Integer space for factors               =',
     &          I16/                                           
     &        ' INFOG (11) Maximum front size                      =',
     &          I16/                                           
     &        ' INFOG (29) Number of entries in factors            =',
     &          I16)                                           
99985 FORMAT( ' INFOG (12) Number of off diagonal pivots           =',
     &          I16)                                           
99986 FORMAT( ' INFOG (25) Number of tiny pivots(static)           =',
     &          I16)                                           
99988 FORMAT( ' Number of 2x2 pivots in type 1 nodes               =',
     &          I16)
99989 FORMAT( ' Number of 2x2 pivots in type 2 nodes               =',
     &          I16)
C
      END SUBROUTINE ZMUMPS_FAC_DRIVER
C
C    
#if defined(SINGLEINDOUBLE)
      SUBROUTINE ZMUMPS_LPid_TO_Cid(id, idintr, K34, IFLAG, IERROR)
      USE ZMUMPS_STRUC_DEF, ONLY: ZMUMPS_STRUC
      USE ZMUMPS_INTR_TYPES, ONLY: ZMUMPS_INTR_STRUC
      ! USE CMUMPS_INTR_TYPES, ONLY: CMUMPS_ROOT_STRUC
      USE CMUMPS_STRUC_DEF, ONLY: CMUMPS_STRUC
      USE CMUMPS_LR_DATA_M, BLR_ARRAY_LP => BLR_ARRAY
      USE ZMUMPS_LR_DATA_M,  BLR_ARRAY_C  => BLR_ARRAY
      USE ZMUMPS_FACSOL_L0OMP_M, ONLY: ZMUMPS_INIT_L0_OMP_FACTORS
      USE CMUMPS_FACSOL_L0OMP_M, ONLY: CMUMPS_FREE_L0_OMP_FACTORS
      use iso_c_binding
      IMPLICIT NONE
C
C     Purpose:
C     =======
C         Move data from Low Precision (LP) instance to 
C         higher precision instance that will be used
C         during the solve phase.
C     Parameters:
C     ===========
          TYPE(ZMUMPS_STRUC), INTENT(INOUT)       :: id
          TYPE (ZMUMPS_INTR_STRUC), intent(inout) :: idintr
          INTEGER :: K34
          INTEGER, intent(inout) :: IFLAG, IERROR
C
C     Local variables:
c     ================
          INTEGER :: IFLAG_LOC, IERROR_LOC
          INTEGER :: allocok
          INTEGER :: I
          IF (associated(idintr%LPL0_OMP_FACTORS)) THEN
            allocate(idintr%L0_OMP_FACTORS(
     &                   size(idintr%lpl0_omp_factors)))
            CALL ZMUMPS_INIT_L0_OMP_FACTORS(idintr%l0_omp_factors)
            DO I=1, size(idintr%LPL0_OMP_FACTORS)
              IF (associated(idintr%LPL0_OMP_FACTORS(I)%A)) THEN
                idintr%L0_OMP_FACTORS(I)%LA =
     &                idintr%LPL0_OMP_FACTORS(I)%LA
                allocate(idintr%L0_OMP_FACTORS(I)%A(
     &                size(idintr%LPL0_OMP_FACTORS(I)%A) ),
     &                stat = allocok)
                IF (allocok .GT. 0) THEN
                  IFLAG = -13
                  IERROR = size(idintr%LPL0_OMP_FACTORS(I)%A)
                  return
                ENDIF
                idintr%L0_OMP_FACTORS(I)%A =
     &                dcmplx(idintr%LPL0_OMP_FACTORS(I)%A)
              ENDIF
            ENDDO
            CALL CMUMPS_FREE_L0_OMP_FACTORS(idintr%lpl0_omp_factors)
          ENDIF
          IF (ASSOCIATED(BLR_ARRAY_LP)) THEN
           IF ( SIZE(BLR_ARRAY_LP) > 0) THEN
            CALL ZMUMPS_BLR_INIT_MODULE (SIZE(BLR_ARRAY_LP), id%INFO 
     &      )
!$OMP       PARALLEL PRIVATE(IFLAG_LOC, IERROR_LOC)
            IFLAG_LOC = 0
C
C
!$OMP       DO SCHEDULE(DYNAMIC,1)
            DO I=1,SIZE(BLR_ARRAY_LP)
              IF (IFLAG_LOC.EQ.0)
     &        CALL ZMUMPS_blr_struc_t_LP2C (
     &          BLR_ARRAY_LP(i),BLR_ARRAY_C(i),K34,
     &          IFLAG_LOC, IERROR_LOC)
            ENDDO
!$OMP       END DO 
!$OMP       CRITICAL
C           IFLAG > 0 on entry to routine indicates warning. 
C           It should be changed
C           only in case of error (IFLAG_LOC > 0)
C           IFLAG should hold the smallest negative value.
            IF (IFLAG_LOC.LT.0.AND.IFLAG_LOC.LT.IFLAG) THEN
                IFLAG = IFLAG_LOC
                IERROR = IERROR_LOC
            ENDIF
!$OMP       END CRITICAL
!$OMP       END PARALLEL
            CALL CMUMPS_BLR_END_MODULE(id%INFO(1), id%KEEP8, K34)
           ENDIF
          ENDIF
          CALL ZMUMPS_BLR_MOD_TO_STRUC(id%BLRARRAY_ENCODING)
          CALL ZMUMPS_BLR_STRUC_TO_MOD(id%BLRARRAY_ENCODING)
          RETURN
      END SUBROUTINE ZMUMPS_LPid_TO_Cid
C
      SUBROUTINE ZMUMPS_blr_struc_t_LP2C (
     &          blr_struc1, blr_struc2, K34, IFLAG, IERROR)
C
          USE CMUMPS_LR_DATA_M, BLR_STRUC_T1 => BLR_STRUC_T
          USE ZMUMPS_LR_DATA_M, BLR_STRUC_T2 => BLR_STRUC_T
          TYPE(BLR_STRUC_T1) :: blr_struc1
          TYPE(BLR_STRUC_T2) :: blr_struc2
          INTEGER, intent(in) :: K34
          INTEGER, intent(inout) :: IFLAG, IERROR
          INTEGER :: k, i, j
          INTEGER :: allocok
          blr_struc2%IsSYM   = blr_struc1%IsSYM
          blr_struc2%IsT2    = blr_struc1%IsT2
          blr_struc2%IsSLAVE = blr_struc1%IsSLAVE
          nullify(blr_struc2%PANELS_L)
          nullify(blr_struc2%PANELS_U)
          nullify(blr_struc2%CB_LRB)
          nullify(blr_struc2%DIAG_BLOCKS)
          IF (ASSOCIATED(blr_struc1%PANELS_L)) THEN
           IF (SIZE(blr_struc1%PANELS_L) > 0) THEN
            ALLOCATE(blr_struc2%PANELS_L(SIZE(blr_struc1%PANELS_L)),
     &        stat=allocok)
            IF (allocok .GT. 0) THEN
              IFLAG = -13
              IERROR = size(blr_struc1%PANELS_L)
              return
            ENDIF
           ENDIF
          ENDIF
          IF (ASSOCIATED(blr_struc1%PANELS_U)) THEN
           IF (SIZE(blr_struc1%PANELS_U) > 0) THEN
            ALLOCATE(blr_struc2%PANELS_U(SIZE(blr_struc1%PANELS_U)),
     &        stat=allocok)
            IF (allocok .GT. 0) THEN
              IFLAG = -13
              IERROR = size(blr_struc1%PANELS_U)
              return
            ENDIF
           ENDIF
          ENDIF
          IF (ASSOCIATED(blr_struc1%DIAG_BLOCKS)) THEN
           IF (SIZE(blr_struc1%DIAG_BLOCKS) > 0) THEN
            ALLOCATE(blr_struc2%DIAG_BLOCKS(
     &               SIZE(blr_struc1%DIAG_BLOCKS)),
     &        stat=allocok)
            IF (allocok .GT. 0) THEN
              IFLAG = -13
              IERROR = size(blr_struc1%DIAG_BLOCKS)
              return
            ENDIF
           ENDIF
          ENDIF
          IF (ASSOCIATED(blr_struc1%PANELS_L)) THEN
           IF (SIZE(blr_struc1%PANELS_L) > 0) THEN
            DO k=1,SIZE(blr_struc1%PANELS_L)
              IF (IFLAG.GE.0)
     &          CALL ZMUMPS_blr_panel_type_cast(
     &                                 blr_struc1%PANELS_L(k), 
     &                                 blr_struc2%PANELS_L(k), K34,
     &                                 IFLAG, IERROR)
            ENDDO
           ENDIF
          ENDIF
          IF (ASSOCIATED(blr_struc1%PANELS_U)) THEN
           IF (SIZE(blr_struc1%PANELS_U) > 0) THEN
            DO k=1,SIZE(blr_struc1%PANELS_U)
              IF (IFLAG.GE.0)
     &          CALL ZMUMPS_blr_panel_type_cast(
     &                                 blr_struc1%PANELS_U(k), 
     &                                 blr_struc2%PANELS_U(k), K34,
     &                                 IFLAG, IERROR)
            ENDDO
           ENDIF
          ENDIF
          IF (ASSOCIATED(blr_struc1%CB_LRB)) THEN 
            ALLOCATE(blr_struc2%CB_LRB(SIZE(blr_struc1%CB_LRB, DIM=1), 
     &                                 SIZE(blr_struc1%CB_LRB, DIM=2)),
     &        stat=allocok)
            IF (allocok .GT. 0) THEN
              IFLAG = -13
              IERROR = SIZE(blr_struc1%CB_LRB, DIM=1) * 
     &                 SIZE(blr_struc1%CB_LRB, DIM=2)
              return
            ENDIF
            DO i=1,SIZE(blr_struc1%CB_LRB, DIM=1)
              DO j=1,SIZE(blr_struc1%CB_LRB, DIM=2)
                IF (IFLAG.GE.0)
     &            CALL ZMUMPS_lrb_type_cast(
     &                             blr_struc1%CB_LRB(i,j),
     &                             blr_struc2%CB_LRB(i,j), K34,
     &                             IFLAG, IERROR)
              ENDDO
            ENDDO
            DEALLOCATE(blr_struc1%CB_LRB)
            NULLIFY(blr_struc1%CB_LRB)
          ENDIF
          IF (ASSOCIATED(blr_struc1%DIAG_BLOCKS)) THEN
           IF (SIZE(blr_struc1%DIAG_BLOCKS) > 0) THEN
            DO k=1,SIZE(blr_struc1%DIAG_BLOCKS)
              IF (IFLAG.GE.0)
     &          CALL ZMUMPS_diag_block_type_cast (
     &                                  blr_struc1%DIAG_BLOCKS(k),
     &                                  blr_struc2%DIAG_BLOCKS(k), K34,
     &                                  iflag, ierror)
            ENDDO
           ENDIF
          ENDIF
          IF (ASSOCIATED(blr_struc1%PANELS_L)) THEN
           IF (SIZE(blr_struc1%PANELS_L) > 0) THEN
            DEALLOCATE(blr_struc1%PANELS_L)
            NULLIFY(blr_struc1%PANELS_L)
           ENDIF
          ENDIF
          IF (ASSOCIATED(blr_struc1%PANELS_U)) THEN
           IF (SIZE(blr_struc1%PANELS_U) > 0) THEN
            DEALLOCATE(blr_struc1%PANELS_U)
            NULLIFY(blr_struc1%PANELS_U)
           ENDIF
          ENDIF
          IF (ASSOCIATED(blr_struc1%DIAG_BLOCKS)) THEN
           IF (SIZE(blr_struc1%DIAG_BLOCKS) > 0) THEN
            DEALLOCATE(blr_struc1%DIAG_BLOCKS)
            NULLIFY(blr_struc1%DIAG_BLOCKS)
           ENDIF
          ENDIF
          blr_struc2%BEGS_BLR_STATIC  => blr_struc1%BEGS_BLR_STATIC
          blr_struc2%BEGS_BLR_DYNAMIC => blr_struc1%BEGS_BLR_DYNAMIC
          blr_struc2%BEGS_BLR_L       => blr_struc1%BEGS_BLR_L
          blr_struc2%BEGS_BLR_COL     => blr_struc1%BEGS_BLR_COL
          NULLIFY(blr_struc1%BEGS_BLR_STATIC )
          NULLIFY(blr_struc1%BEGS_BLR_DYNAMIC)
          NULLIFY(blr_struc1%BEGS_BLR_L      )
          NULLIFY(blr_struc1%BEGS_BLR_COL    )
          blr_struc2%NB_ACCESSES_INIT = blr_struc1%NB_ACCESSES_INIT
          blr_struc2%NB_PANELS        = blr_struc1%NB_PANELS
          blr_struc2%NFS4FATHER       = blr_struc1%NFS4FATHER
C
          IF (ASSOCIATED(blr_struc1%M_ARRAY)) THEN
           IF ( SIZE(blr_struc1%M_ARRAY) > 0) THEN
            ALLOCATE(blr_struc2%M_ARRAY(SIZE(blr_struc1%M_ARRAY)),
     &        stat=allocok)
            IF (allocok .GT. 0) THEN
              IFLAG = -13
              IERROR = SIZE(blr_struc1%M_ARRAY)
              return
            ENDIF
            DO i=1,SIZE(blr_struc1%M_ARRAY)
                blr_struc2%M_ARRAY(i) = DCMPLX(blr_struc1%M_ARRAY(i))
            ENDDO
            !blr_struc2%M_ARRAY = dble(blr_struc1%M_ARRAY, KIND=8)
            DEALLOCATE(blr_struc1%M_ARRAY)
            NULLIFY(blr_struc1%M_ARRAY)
           ENDIF
          ENDIF
      END SUBROUTINE ZMUMPS_blr_struc_t_LP2C
      SUBROUTINE ZMUMPS_diag_block_type_cast(
     &        diag_block1, diag_block2, K34, IFLAG, IERROR)
          USE ZMUMPS_LR_TYPE
          USE CMUMPS_LR_DATA_M, DIAG_BLOCK_TYPE1 => DIAG_BLOCK_TYPE
          USE ZMUMPS_LR_DATA_M, DIAG_BLOCK_TYPE2 => DIAG_BLOCK_TYPE
          use iso_c_binding
          TYPE(DIAG_BLOCK_TYPE1), intent(inout) :: diag_block1
          TYPE(DIAG_BLOCK_TYPE2), intent(inout) :: diag_block2
          INTEGER, intent(in):: K34
          INTEGER, intent(inout):: IFLAG, IERROR
#if !defined(NOMIXHPC)
          integer :: int_sz
#endif
          INTEGER :: allocok
          nullify(diag_block2%diag_block)
#if !defined(NOMIXHPC)
          nullify(diag_block2%diag_stoc%val)
          if (associated(diag_block1%diag_stoc%val)) then
             diag_block2%diag_size = diag_block1%diag_size
             diag_block2%diag_stoc%val => diag_block1%diag_stoc%val
             nullify(diag_block1%diag_block)
          else if (associated(diag_block1%diag_block)) then
             diag_block2%diag_size = size(diag_block1%diag_block)
             int_sz=1+(costs_prec(PREC_FAC)*diag_block2%diag_size-1)/K34
             if (PREC_DIAG_FR.EQ.PREC_FAC) then
               int_sz = 1 +
     &           (costs_prec(PREC_DIAG_FR)*diag_block2%diag_size-1)/K34
               allocate(diag_block2%diag_stoc%val(int_sz), stat=allocok)
               IF (allocok .GT. 0) THEN
                 IFLAG = -13
                 IERROR = int_sz
                 return
               ENDIF
               CALL ZMUMPS_LRMIXED_CAST(
     &           diag_block1%diag_block(1), 
     &           diag_block2%diag_stoc%val(1),
     &           diag_block2%DIAG_SIZE, 1,
     &           diag_block2%DIAG_SIZE,
     &           diag_block2%DIAG_SIZE,
     &           PREC_FAC, PREC_DIAG_FR)
             else
               allocate(diag_block2%diag_block(diag_block2%diag_size),
     &            stat=allocok)
               IF (allocok .GT. 0) THEN
                 IFLAG = -13
                 IERROR = diag_block2%diag_size
                 return
               ENDIF
               CALL ZMUMPS_LRMIXED_CAST(
     &           diag_block1%diag_block(1), 
     &           diag_block2%diag_block(1),
     &           diag_block2%DIAG_SIZE, 1,
     &           diag_block2%DIAG_SIZE,
     &           diag_block2%DIAG_SIZE,
     &           PREC_FAC, PREC_DIAG_FR)
             endif
             deallocate(diag_block1%diag_block)
          endif
#else
          ! NOMIXHPC defined
          if (associated(diag_block1%diag_block)) then
              allocate(diag_block2%diag_block(
     &          size(diag_block1%diag_block)), stat=allocok)
              IF (allocok .GT. 0) THEN
                IFLAG = -13
                IERROR = size(diag_block1%diag_block)
                return
              ENDIF
              CALL ZMUMPS_LRMIXED_CAST(
     &          diag_block1%diag_block(1), 
     &          diag_block2%diag_block(1),
     &          size(diag_block1%DIAG_BLOCK), 1,
     &          size(diag_block1%DIAG_BLOCK),
     &          size(diag_block2%DIAG_BLOCK),
     &          PREC_FAC, PREC_DIAG_FR)
              deallocate(diag_block1%diag_block)
          endif
#endif
      END SUBROUTINE ZMUMPS_diag_block_type_cast
      SUBROUTINE ZMUMPS_blr_panel_type_cast(
     &    blr_panel1, blr_panel2, K34,
     &    IFLAG, IERROR)
          USE CMUMPS_LR_DATA_M, BLR_PANEL_TYPE1 => BLR_PANEL_TYPE
          USE ZMUMPS_LR_DATA_M, BLR_PANEL_TYPE2 => BLR_PANEL_TYPE
          TYPE(BLR_PANEL_TYPE1),intent(inout) :: blr_panel1
          TYPE(BLR_PANEL_TYPE2),intent(inout) :: blr_panel2
          INTEGER, intent(in):: K34
          INTEGER, intent(inout):: IFLAG, IERROR
          INTEGER :: k, allocok
          blr_panel2%NB_ACCESSES_LEFT = blr_panel1%NB_ACCESSES_LEFT 
          nullify(blr_panel2%lrb_panel)
          IF (ASSOCIATED(blr_panel1%LRB_PANEL)) THEN
            IF (SIZE(blr_panel1%LRB_PANEL).GE.0) THEN
              ALLOCATE(blr_panel2%LRB_PANEL(SIZE(blr_panel1%LRB_PANEL)),
     &           stat= allocok)
              IF (allocok .GT. 0) THEN
                IFLAG = -13
                IERROR = size(blr_panel1%lrb_panel)
                return
              ENDIF
              DO k=1,SIZE(blr_panel1%LRB_PANEL)
                CALL ZMUMPS_lrb_type_cast(blr_panel1%LRB_PANEL(k), 
     &                             blr_panel2%LRB_PANEL(k), K34,
     &                             IFLAG, IERROR)
              ENDDO
            ENDIF
            DEALLOCATE(blr_panel1%LRB_PANEL)
            NULLIFY(blr_panel1%LRB_PANEL)
          ENDIF
      END SUBROUTINE ZMUMPS_blr_panel_type_cast
      SUBROUTINE ZMUMPS_lrb_type_cast(
     &                lrb1, lrb2, K34, IFLAG, IERROR)
          USE CMUMPS_LR_TYPE, LRB_TYPE1 => LRB_TYPE
          USE ZMUMPS_LR_TYPE, LRB_TYPE2 => LRB_TYPE
          TYPE(LRB_TYPE1) :: lrb1
          TYPE(LRB_TYPE2) :: lrb2
          INTEGER, intent(in):: K34
          INTEGER, intent(inout):: IFLAG, IERROR
          INTEGER         :: i
          INTEGER         :: int_sz, allocok
          integer         :: m, n
          lrb2%K    = lrb1%K    
          lrb2%M    = lrb1%M    
          lrb2%N    = lrb1%N    
          lrb2%ISLR = lrb1%ISLR 
          nullify(lrb2%Q_stoc)
          nullify(lrb2%R_stoc)
          nullify(lrb2%K_stoc)
          nullify(lrb2%Q)
          nullify(lrb2%R)
          if (associated(lrb1%Q_stoc).or.
     &        associated(lrb1%Q)) then
            allocate(lrb2%Q_stoc(NPREC_STOC), stat=allocok)
            IF (allocok .GT. 0) THEN
              IFLAG = -13
              IERROR = NPREC_STOC
              return
            ENDIF
            do i = 1, NPREC_STOC
              nullify(lrb2%Q_stoc(i)%val)
            end do
          end if
          if (associated(lrb1%R_stoc).or.
     &        associated(lrb1%R)) then
            allocate(lrb2%R_stoc(NPREC_STOC), stat=allocok)
            IF (allocok .GT. 0) THEN
              IFLAG = -13
              IERROR = NPREC_STOC
              return
            ENDIF
            do i = 1, NPREC_STOC
              nullify(lrb2%R_stoc(i)%val)
            end do
          end if
          if (associated(lrb1%Q_stoc)) then
            do i = 1, NPREC_STOC
              lrb2%Q_stoc(i)%val => lrb1%Q_stoc(i)%val
              nullify(lrb1%Q_stoc(i)%val)
            end do
            deallocate(lrb1%Q_stoc)
          end if
          if (associated(lrb1%R_stoc)) then
            do i = 1, NPREC_STOC
              lrb2%R_stoc(i)%val => lrb1%R_stoc(i)%val
              nullify(lrb1%R_stoc(i)%val)
            end do
            deallocate(lrb1%R_stoc)
          end if
          ! this seems easily possible to just swap pointers here
          if (associated(lrb1%k_stoc)) then
            allocate(lrb2%K_stoc(nprec_stoc), stat=allocok)
            IF (allocok .GT. 0) THEN
              IFLAG = -13
              IERROR = nprec_stoc
              return
            ENDIF
            do i = 1, NPREC_STOC
              lrb2%K_stoc(i) = lrb1%K_stoc(i)
            end do
            deallocate(lrb1%K_stoc)
          endif
          if (associated(lrb1%Q)) then
             if (lrb1%islr) then
                 m = lrb1%m
                 n = lrb1%k
             else
                 m = lrb1%m
                 n = lrb1%n
             endif
             int_sz = 1 + (COSTS_PREC(PREC_FAC)*size(lrb1%Q)-1)/K34
#if !defined(NOMIXHPC)
             allocate(lrb2%Q_stoc(PREC_FAC)%val(int_sz), stat=allocok)
             IF (allocok .GT. 0) THEN
               IFLAG = -13
               IERROR = int_sz
               return
             ENDIF
             CALL ZMUMPS_LRMIXED_CAST(
     &         lrb1%Q(1,1), lrb2%Q_stoc(PREC_FAC)%val(1),
     &         m, n,
     &         lrb1%M, lrb2%M,
     &         PREC_FAC, PREC_FAC)
#else
             if (lrb1%islr) then
               allocate(lrb2%Q_stoc(PREC_FAC)%val(int_sz), stat=allocok)
               IF (allocok .GT. 0) THEN
                 IFLAG = -13
                 IERROR = int_sz
                 return
               ENDIF
               CALL ZMUMPS_LRMIXED_CAST(
     &           lrb1%Q(1,1), lrb2%Q_stoc(PREC_FAC)%val(1),
     &           m, n,
     &           lrb1%M, lrb2%M,
     &           PREC_FAC, PREC_SGL)
             else
               allocate(lrb2%Q( size(lrb1%Q, 1), size(lrb1%Q, 2)),
     &           stat=allocok)
               IF (allocok .GT. 0) THEN
                 IFLAG = -13
                 IERROR = size(lrb1%Q, 1) * size(lrb1%Q, 2)
                 return
               ENDIF
               CALL ZMUMPS_LRMIXED_CAST(
     &           lrb1%Q(1,1), lrb2%Q(1,1),
     &           m, n,
     &           lrb1%M, lrb2%M,
     &           PREC_FAC, 1)
             endif
#endif
             deallocate(lrb1%Q)
          end if
          if (associated(lrb1%R)) then
             int_sz = 1 + (COSTS_PREC(PREC_FAC)*size(lrb1%R)-1)/K34
             allocate(lrb2%R_stoc(PREC_FAC)%val(int_sz), stat=allocok)
             IF (allocok .GT. 0) THEN
               IFLAG = -13
               IERROR = int_sz
               return
             ENDIF
             CALL ZMUMPS_LRMIXED_CAST(
     &         lrb1%R(1,1), lrb2%R_stoc(PREC_FAC)%val(1),
     &         lrb1%K, lrb1%N,
     &         lrb1%K, lrb2%K,
     &         PREC_FAC, PREC_FAC)
             deallocate(lrb1%R)
             if (.not.associated(lrb2%k_stoc)) then
               allocate(lrb2%k_stoc(NPREC_STOC), stat=allocok)
               IF (allocok .GT. 0) THEN
                 IFLAG = -13
                 IERROR = NPREC_STOC
                 return
               ENDIF
             end if
             do i = 1, NPREC_STOC
               lrb2%K_stoc(i) = 0
             end do
             lrb2%K_stoc(PREC_FAC) = lrb2%k
          end if
      END SUBROUTINE ZMUMPS_lrb_type_cast
#endif
C
      SUBROUTINE ZMUMPS_PRINT_ALLOCATED_MEM( PROK, PROKG, PRINT_MAXAVG,
     &       MP, MPG, INFO16, INFOG18, INFOG19, NSLAVES, IRANK, KEEP )
      IMPLICIT NONE
C
C  Purpose:
C  =======
C     Print memory allocated during factorization
C     - called at beginning of factorization in full-rank
C     - called at end of factorization in low-rank (because
C       of dynamic allocations)
C
      LOGICAL, INTENT(IN) :: PROK, PROKG, PRINT_MAXAVG
      INTEGER, INTENT(IN) :: MP, MPG, INFO16, INFOG18, INFOG19
      INTEGER, INTENT(IN) :: IRANK, NSLAVES
      INTEGER, INTENT(IN) :: KEEP(500)
C
      IF ( PROKG ) THEN
        IF (PRINT_MAXAVG) THEN
         WRITE( MPG,'(A,I12) ')
     &   ' ** Memory allocated, max in Mbytes             (INFOG(18)):',
     &   INFOG18
        ENDIF
        WRITE( MPG,'(/A,I12) ')
     &   ' ** Memory allocated, total in Mbytes           (INFOG(19)):',
     &    INFOG19
      END IF
      RETURN
      END SUBROUTINE ZMUMPS_PRINT_ALLOCATED_MEM
      SUBROUTINE ZMUMPS_AVGMAX_STAT8(PROKG, MPG, VAL, NSLAVES,
     &     PRINT_MAXAVG, COMM, MSG)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      LOGICAL, intent(in) :: PROKG
      INTEGER, intent(in) :: MPG
      INTEGER(8), intent(in) :: VAL
      INTEGER, intent(in) :: NSLAVES
      LOGICAL, intent(in) :: PRINT_MAXAVG
      INTEGER, intent(in) :: COMM
      CHARACTER*48 MSG 
C  Local
      INTEGER(8) MAX_VAL
      INTEGER IERR, MASTER
      DOUBLE PRECISION LOC_VAL, AVG_VAL
      PARAMETER(MASTER=0)
C
      CALL MUMPS_REDUCEI8( VAL, MAX_VAL, MPI_MAX, MASTER, COMM)
      LOC_VAL = dble(VAL)/dble(NSLAVES)
      CALL MPI_REDUCE( LOC_VAL, AVG_VAL, 1, MPI_DOUBLE_PRECISION,
     &                 MPI_SUM, MASTER, COMM, IERR )
      IF (PROKG) THEN
        IF (PRINT_MAXAVG) THEN
          WRITE(MPG,100) " Average", MSG, int(AVG_VAL,8)
        ELSE
          WRITE(MPG,110)  MSG, MAX_VAL
        ENDIF
      ENDIF
      RETURN
 100  FORMAT(A8,A48,I18)
 110  FORMAT(A48,I18)
      END SUBROUTINE ZMUMPS_AVGMAX_STAT8
C
C
C     ==================================================================
C
      SUBROUTINE ZMUMPS_EXTRACT_SCALING(id)
      USE ZMUMPS_STRUC_DEF
      IMPLICIT NONE
C
C     Purpose:
C     =======
C
C     Extract distributed scaling arrays from ZMUMPS_EXTRACT_SCALING
C     In case of unsymmetric permutation, ROWSCA and COLSCA correspond
C     to Dr and Dc, in the expression Dr A Q Dc. In other terms, Dc
C     is compatbile with the front column indices, it does not
C     correspond to the column indices of A, meaning that Q is not
C     needed to just extract the scaling values.
C
C
      TYPE(ZMUMPS_STRUC) :: id
      INTEGER, EXTERNAL  :: MUMPS_PROCNODE
C
C     MPI
C     ===
C
      INCLUDE 'mpif.h'
C
C     Local declarations
C     ==================
C
      DOUBLE PRECISION, DIMENSION(:), POINTER :: COLSCA
      DOUBLE PRECISION, DIMENSION(:), POINTER :: ROWSCA
      INTEGER, PARAMETER :: MASTER = 0
C
      INTEGER    :: ISTEP, NPIV, LIELL
      INTEGER    :: IERR_MPI, allocok
      INTEGER    :: ISCA
      INTEGER    :: JROW, JCOL, IPOS, JJ ! access to IS
      INTEGER    :: LIW_PASSED
      INTEGER(8) :: LALLOC
C
C     Free and reallocate distributed scaling arrays :
C     - in symmetric, COLSCA_loc points on ROWSCA_loc.
C     - not allocated if KEEP(89)=0
C
      NULLIFY(ROWSCA)
      NULLIFY(COLSCA)
      IF (associated(id%ROWSCA_loc)) THEN
        DEALLOCATE(id%ROWSCA_loc)
        NULLIFY(id%ROWSCA_loc)
      ENDIF
      IF (associated(id%COLSCA_loc)) THEN
        IF (id%KEEP(50) .EQ. 0) THEN
          DEALLOCATE(id%COLSCA_loc)
        ENDIF
        NULLIFY(id%COLSCA_loc)
      ENDIF
C     id%ROWSCA and id%COLSCA are available on master,
C     allocate ROWSCA and COLSCA of order N on other procs
      IF ( id%MYID .EQ. MASTER ) THEN
        ROWSCA => id%ROWSCA
        COLSCA => id%COLSCA
        IF (.NOT. associated(ROWSCA)) THEN
          WRITE(*,*) "Internal error 1 in ZMUMPS_EXTRACT_SCALING"
          CALL MUMPS_ABORT()
        ENDIF
        IF (.NOT. associated(COLSCA)) THEN
          WRITE(*,*) "Internal error 2 in ZMUMPS_EXTRACT_SCALING"
          CALL MUMPS_ABORT()
        ENDIF
      ELSE
        IF (id%KEEP(50).EQ.0) THEN
          ALLOCATE(ROWSCA(id%N),COLSCA(id%N),stat=allocok)
          LALLOC = int(id%N+id%N,8)
        ELSE
          ALLOCATE(ROWSCA(id%N),stat=allocok)
          COLSCA => ROWSCA
          LALLOC = int(id%N,8)
        ENDIF
        IF (allocok .GT. 0) THEN
          id%INFO(1)=-13
          CALL MUMPS_SET_IERROR(LALLOC,id%INFO(2))
        ENDIF
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
C     Jump to 110 in case of error on ROWSCA or COLSCA
C     on one of the MPI processes.
      IF (id%INFO(1) .LT. 0) GOTO 110
C
      IF ( id%KEEP(89) .GT. 0) THEN
        IF (id%KEEP(50).EQ.0) THEN
          ALLOCATE(id%ROWSCA_loc(id%KEEP(89)),
     &             id%COLSCA_loc(id%KEEP(89)),stat=allocok)
          LALLOC = int(id%KEEP(89),8)*2_8
        ELSE
          ALLOCATE(id%ROWSCA_loc(id%KEEP(89)),stat=allocok)
          id%COLSCA_loc => id%ROWSCA_loc
          LALLOC = int(id%KEEP(89),8)
        ENDIF
        IF (allocok .GT. 0) THEN
          id%INFO(1)=-13
          CALL MUMPS_SET_IERROR(LALLOC,id%INFO(2))
        ENDIF
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
C     Jump to 100 in case of error (we free everything)
      IF (id%INFO(1) .LT. 0) GOTO 100
      CALL MPI_BCAST(ROWSCA(1), id%N, MPI_DOUBLE_PRECISION,
     &               MASTER, id%COMM, IERR_MPI)
      IF (id%KEEP(50) .EQ. 0) THEN
        CALL MPI_BCAST(COLSCA(1), id%N, MPI_DOUBLE_PRECISION,
     &                 MASTER, id%COMM, IERR_MPI)
      ENDIF
      LIW_PASSED = max(id%KEEP(32),1)
      ISCA = 1
      IF ( id%MYID .ne. MASTER .OR.
     &     id%KEEP(46) .eq. 1 ) THEN ! I_AM_SLAVE
        DO ISTEP = 1, id%KEEP(28)
          IF ( id%MYID_NODES.EQ. MUMPS_PROCNODE(
     &         id%PROCNODE_STEPS(ISTEP),
     &         id%KEEP(199) ) ) THEN
              CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, id%KEEP,
     &        NPIV, LIELL, IPOS,
     &        id%IS(1), LIW_PASSED, id%PTLUST_S(1), id%STEP(1), id%N)
            IF ( id%KEEP(50) .EQ. 0 ) THEN
C             Row indices:
              JROW = IPOS + 1
C             Column indices:
              JCOL = IPOS + 1 + LIELL
            ELSE
C             Use row indices because column indices may have
C             been set to negative to flag 2x2 pivots
              JROW = IPOS + 1
            ENDIF
            IF (id%KEEP(50).EQ.0) THEN
              DO JJ = 1, NPIV
                id%ROWSCA_loc(ISCA+JJ-1) = ROWSCA(id%IS(JROW+JJ-1))
                id%COLSCA_loc(ISCA+JJ-1) = COLSCA(id%IS(JCOL+JJ-1))
              ENDDO
            ELSE
              DO JJ = 1, NPIV
                id%ROWSCA_loc(ISCA+JJ-1) = ROWSCA(id%IS(JROW+JJ-1))
              ENDDO
            ENDIF
            ISCA = ISCA + NPIV
          ENDIF
        ENDDO
      ENDIF
C     End of EXTRACT_SCALING, we keep id%ROWSCA_loc and id%COLSCA_loc
C     but free ROWSCA and COLSCA
      GOTO 110
 100  CONTINUE
C     Exit with error, free what was allocated
      IF (associated(id%ROWSCA_loc)) THEN
        DEALLOCATE(id%ROWSCA_loc)
        NULLIFY(id%ROWSCA_loc)
      ENDIF
      IF (associated(id%COLSCA_loc)) THEN
        IF (id%KEEP(50) .EQ. 0) THEN
          DEALLOCATE(id%COLSCA_loc)
        ENDIF
        NULLIFY(id%COLSCA_loc)
      ENDIF
 110  CONTINUE
C     Free local ROWSCA and COLSCA arrays
      IF ( id%MYID .NE. 0) THEN
        IF (associated(ROWSCA)) DEALLOCATE(ROWSCA)
        IF ( id%KEEP(50) .EQ. 0 ) THEN
          IF (associated(COLSCA)) DEALLOCATE(COLSCA)
        ENDIF
      ENDIF
      NULLIFY(ROWSCA)
      NULLIFY(COLSCA)
      RETURN
      END SUBROUTINE ZMUMPS_EXTRACT_SCALING
C
C     ==================================================================
C
      SUBROUTINE ZMUMPS_EXTRACT_SCHUR_REDRHS(id,idintr)
      USE ZMUMPS_STRUC_DEF
      USE ZMUMPS_INTR_TYPES, ONLY: ZMUMPS_INTR_STRUC
      IMPLICIT NONE
C
C  Purpose
C  =======
C
C     Extract the Schur and possibly also the reduced right-hand side
C     (if Fwd in facto) from the processor working on Schur and copy
C     it into the user datastructures id%SCHUR and id%REDRHS on the host.
C     This routine assumes that the integer list of the Schur has not
C     been permuted and still corresponds to LISTVAR_SCHUR.
C
C     If the Schur is centralized, the master of the Schur holds the
C     Schur and possibly also the reduced right-hand side.
C     If the Schur is distribued (already built in user's datastructure),
C     then the master of the Schur may hold the reduced right-hand side,
C     in which case it is available in roota%RHS_CNTR_MASTER_ROOT.
C     
      TYPE (ZMUMPS_STRUC)      :: id
      TYPE (ZMUMPS_INTR_STRUC) :: idintr
C
C  Local variables
C  ===============
C
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mumps_headers.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR
      INTEGER, PARAMETER :: MASTER = 0
      INTEGER :: ID_SCHUR, SIZE_SCHUR, LD_SCHUR, IB, BL4
      INTEGER(4) :: I4 ! 32-bit even in 64-bit version
      INTEGER :: ROW_LENGTH, I
      INTEGER(8) :: SURFSCHUR8, BL8, SHIFT8
      INTEGER(8) :: ISCHUR_SRC, ISCHUR_DEST, ISCHUR_SYM, ISCHUR_UNS
C
C  External functions
C  ==================
C
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
C     Quick return in case factorization did not terminate correctly
      IF (id%INFO(1) .LT. 0) RETURN
C     Quick return if Schur option off
      IF (id%KEEP(60) .EQ. 0) RETURN
C     Get Schur id
      ID_SCHUR =MUMPS_PROCNODE(
     &    id%PROCNODE_STEPS(id%STEP(max(id%KEEP(20),id%KEEP(38)))),
     &    id%KEEP(199))
      IF ( id%KEEP( 46 )  .NE. 1 ) THEN
        ID_SCHUR = ID_SCHUR + 1
      END IF
C     Get size of Schur
      IF (id%MYID.EQ.ID_SCHUR) THEN
        IF (id%KEEP(60).EQ.1) THEN
C         Sequential Schur
          LD_SCHUR =
     &    id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))+2+id%KEEP(IXSZ))
          SIZE_SCHUR = LD_SCHUR - id%KEEP(253)
        ELSE
C         Parallel Schur
          LD_SCHUR   = -999999 ! not used
          SIZE_SCHUR = idintr%root%TOT_ROOT_SIZE
        ENDIF
      ELSE IF (id%MYID .EQ. MASTER) THEN
        SIZE_SCHUR = id%KEEP(116)
        LD_SCHUR = -44444 ! Not used
      ELSE
C       Proc is not concerned with Schur, return
        RETURN
      ENDIF
      SURFSCHUR8 = int(SIZE_SCHUR,8)*int(SIZE_SCHUR,8)
C     =====================================
C     Case of parallel Schur: if REDRHS
C     was requested, obtain it directly
C     from idintr%roota%RHS_CNTR_MASTER_ROOT
C     =====================================
      IF (id%KEEP(60) .GT. 1) THEN
        IF (id%KEEP(221).EQ.1 .AND. id%KEEP(252).GT.0) THEN
          DO I = 1, id%KEEP(253)
            IF (ID_SCHUR.EQ.MASTER) THEN ! Necessarily = id%MYID
              CALL zcopy(SIZE_SCHUR,
     &        idintr%roota%RHS_CNTR_MASTER_ROOT((I-1)*SIZE_SCHUR+1), 1,
     &        id%REDRHS((I-1)*id%LREDRHS+1), 1)
            ELSE
              IF (id%MYID.EQ.ID_SCHUR) THEN
C               Send
                CALL MPI_SEND(
     &          idintr%roota%RHS_CNTR_MASTER_ROOT((I-1)*SIZE_SCHUR+1),
     &          SIZE_SCHUR,
     &          MPI_DOUBLE_COMPLEX,
     &          MASTER, TAG_SCHUR,
     &          id%COMM, IERR )
              ELSE ! MYID.EQ.MASTER
C               Receive
                CALL MPI_RECV( id%REDRHS((I-1)*id%LREDRHS+1),
     &             SIZE_SCHUR,
     &             MPI_DOUBLE_COMPLEX, ID_SCHUR, TAG_SCHUR,
     &             id%COMM, STATUS, IERR )
              ENDIF
            ENDIF
          ENDDO
C         ------------------------------
C         In case of parallel Schur, we
C         free roota%RHS_CNTR_MASTER_ROOT
C         ------------------------------
          IF (id%MYID.EQ.ID_SCHUR) THEN
            DEALLOCATE(idintr%roota%RHS_CNTR_MASTER_ROOT)
            NULLIFY   (idintr%roota%RHS_CNTR_MASTER_ROOT)
          ENDIF
        ENDIF
C       return because this is all we need to do
C       in case of parallel Schur complement
        RETURN
      ENDIF
C     ============================
C     Centralized Schur complement
C     ============================
C     PTRAST has been freed at the moment of calling this
C     routine. Schur is available through
C     PTRFAC(IW( PTLUST_S( STEP(KEEP(20)) ) + 4 +KEEP(IXSZ) ))
      IF (id%KEEP(252).EQ.0) THEN
C       CASE 1 (ORIGINAL CODE):
C       Schur is contiguous on ID_SCHUR
        IF ( ID_SCHUR .EQ. MASTER ) THEN ! Necessarily equals id%MYID
C         ---------------------
C         Copy Schur complement
C         ---------------------
          CALL ZMUMPS_COPYI8SIZE( SURFSCHUR8,
     &      id%S(id%PTRFAC(id%STEP(id%KEEP(20)))),
     &      id%SCHUR(1) )
        ELSE
C         -----------------------------------------
C         The processor responsible of the Schur
C         complement sends it to the host processor
C         Use blocks to avoid too large messages.
C         -----------------------------------------
          BL8=int(huge(I4)/id%KEEP(35)/10,8)
          DO IB=1, int((SURFSCHUR8+BL8-1_8) / BL8)
            SHIFT8 = int(IB-1,8) * BL8                ! Where to send
            BL4    = int(min(BL8,SURFSCHUR8-SHIFT8)) ! Size of block
            IF ( id%MYID .eq. ID_SCHUR ) THEN
C             Send Schur complement
              CALL MPI_SEND( id%S( SHIFT8 +
     &          id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     &                    +4+id%KEEP(IXSZ)))),
     &          BL4,
     &          MPI_DOUBLE_COMPLEX,
     &          MASTER, TAG_SCHUR,
     &          id%COMM, IERR )
            ELSE IF ( id%MYID .eq. MASTER ) THEN
C             Receive Schur complement
              CALL MPI_RECV( id%SCHUR(1_8 + SHIFT8),
     &                     BL4,
     &                     MPI_DOUBLE_COMPLEX, ID_SCHUR, TAG_SCHUR,
     &                     id%COMM, STATUS, IERR )
            END IF
          ENDDO
        END IF
      ELSE
C       CASE 2 (Fwd in facto): Schur is not contiguous on ID_SCHUR,
C       process it row by row.
C
C       2.1: We first centralize Schur complement into id%SCHUR
        ISCHUR_SRC = id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     &               +4+id%KEEP(IXSZ)))
        ISCHUR_DEST= 1_8
        DO I=1, SIZE_SCHUR
          ROW_LENGTH = SIZE_SCHUR
          IF (ID_SCHUR.EQ.MASTER) THEN ! Necessarily = id%MYID
            CALL zcopy(ROW_LENGTH, id%S(ISCHUR_SRC), 1,
     &                 id%SCHUR(ISCHUR_DEST),1)
          ELSE
            IF (id%MYID.EQ.ID_SCHUR) THEN
C             Send
              CALL MPI_SEND( id%S(ISCHUR_SRC), ROW_LENGTH,
     &        MPI_DOUBLE_COMPLEX,
     &        MASTER, TAG_SCHUR,
     &        id%COMM, IERR )
            ELSE
C             Recv
              CALL MPI_RECV( id%SCHUR(ISCHUR_DEST),
     &                   ROW_LENGTH,
     &                   MPI_DOUBLE_COMPLEX, ID_SCHUR, TAG_SCHUR,
     &                   id%COMM, STATUS, IERR )
            ENDIF
          ENDIF
          ISCHUR_SRC = ISCHUR_SRC+int(LD_SCHUR,8)
          ISCHUR_DEST= ISCHUR_DEST+int(SIZE_SCHUR,8)
        ENDDO
C       2.2: Get REDRHS on host
C       2.2.1: Symmetric => REDRHS is available in last KEEP(253)
C              rows of Schur structure on ID_SCHUR
C       2.2.2: Unsymmetric => REDRHS corresponds to last KEEP(253)
C              columns. However it must be transposed.
        IF (id%KEEP(221).EQ.1) THEN ! Implies Fwd in facto
          ISCHUR_SYM = id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     &                    +4+id%KEEP(IXSZ))) + int(SIZE_SCHUR,8) *
     &                    int(LD_SCHUR,8)
          ISCHUR_UNS =
     &                 id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     &                    +4+id%KEEP(IXSZ))) + int(SIZE_SCHUR,8)
          ISCHUR_DEST = 1_8
          DO I = 1, id%KEEP(253)
            IF (ID_SCHUR .EQ. MASTER) THEN ! necessarily = id%MYID
              IF (id%KEEP(50) .EQ. 0) THEN
                CALL zcopy(SIZE_SCHUR, id%S(ISCHUR_UNS), LD_SCHUR,
     &                     id%REDRHS(ISCHUR_DEST), 1)
              ELSE
                CALL zcopy(SIZE_SCHUR, id%S(ISCHUR_SYM), 1,
     &                     id%REDRHS(ISCHUR_DEST), 1)
              ENDIF
            ELSE
              IF (id%MYID .NE. MASTER) THEN
                IF (id%KEEP(50) .EQ. 0) THEN
C                 Use id%S(ISCHUR_SYM) as temporary contig. workspace
C                 of size SIZE_SCHUR. 
                  CALL zcopy(SIZE_SCHUR, id%S(ISCHUR_UNS), LD_SCHUR,
     &            id%S(ISCHUR_SYM), 1)
                ENDIF
                CALL MPI_SEND(id%S(ISCHUR_SYM), SIZE_SCHUR,
     &          MPI_DOUBLE_COMPLEX, MASTER, TAG_SCHUR,
     &          id%COMM, IERR )
              ELSE
                CALL MPI_RECV(id%REDRHS(ISCHUR_DEST),
     &          SIZE_SCHUR, MPI_DOUBLE_COMPLEX, ID_SCHUR, TAG_SCHUR,
     &          id%COMM, STATUS, IERR )
              ENDIF
            ENDIF
            IF (id%KEEP(50).EQ.0) THEN
              ISCHUR_UNS = ISCHUR_UNS + int(LD_SCHUR,8)
            ELSE
              ISCHUR_SYM = ISCHUR_SYM + int(LD_SCHUR,8)
            ENDIF
            ISCHUR_DEST = ISCHUR_DEST + int(id%LREDRHS,8)
          ENDDO
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_EXTRACT_SCHUR_REDRHS
      SUBROUTINE ZMUMPS_EXTRACT_SINGULAR_VALUES(id,idintr)
      USE ZMUMPS_STRUC_DEF
      USE ZMUMPS_INTR_TYPES, ONLY: ZMUMPS_INTR_STRUC
      IMPLICIT NONE
C
C  Purpose
C  =======
C
C     
      TYPE (ZMUMPS_STRUC)      :: id
      TYPE (ZMUMPS_INTR_STRUC) :: idintr
C
C  Local variables
C  ===============
C
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mumps_headers.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR
      INTEGER, PARAMETER :: MASTER = 0
      INTEGER :: ID_ROOT, ALLOCOK
C
C  External functions
C  ==================
C
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
C     Quick return in case factorization did not terminate correctly
      IF (id%INFO(1) .LT. 0) RETURN
C     Quick return if Postponing + rank revealing option off
      IF (id%KEEP(19) .EQ. 0) RETURN
C     Get Root id
      ID_ROOT =MUMPS_PROCNODE(
     &    id%PROCNODE_STEPS(id%STEP(id%KEEP(20))),
     &    id%KEEP(199))
      IF ( id%KEEP( 46 )  .NE. 1 ) THEN
        ID_ROOT = ID_ROOT + 1
      END IF  
C     =================================
C     Singular values are stored in 
C     roota%SINGULAR_VALUES
C     We copy it to id%SINGULAR_VALUES
C     =================================
      IF ((ID_ROOT.EQ.id%MYID).AND.(id%MYID.EQ.MASTER)) THEN
C        write(6,*) " singular_values already on host"
         IF (associated(id%SINGULAR_VALUES))
     &        DEALLOCATE(id%SINGULAR_VALUES)
         NULLIFY(id%SINGULAR_VALUES)
         id%NB_SINGULAR_VALUES=idintr%root%NB_SINGULAR_VALUES
         ALLOCATE(id%SINGULAR_VALUES(id%NB_SINGULAR_VALUES)
     &        , stat=ALLOCOK)
         IF (ALLOCOK.GT.0) THEN
            id%INFO(1)=-13
            id%INFO(2)= id%NB_SINGULAR_VALUES
            RETURN
         END IF
         CALL dcopy(id%NB_SINGULAR_VALUES,
     &        idintr%roota%SINGULAR_VALUES(1), 1,
     &        id%SINGULAR_VALUES(1), 1)
      ELSE
         IF (id%MYID.EQ.ID_ROOT) THEN
C     Send
C           write(6,*) " id%MYID sends singular_values "
            CALL MPI_SEND(
     &           idintr%root%NB_SINGULAR_VALUES,
     &           1,
     &           MPI_INTEGER,
     &           MASTER, TAG_ROOT1,
     &           id%COMM, IERR )
            CALL MPI_SEND(
     &           idintr%roota%SINGULAR_VALUES(1),
     &           idintr%root%NB_SINGULAR_VALUES,
     &           MPI_DOUBLE_PRECISION,
     &           MASTER, TAG_ROOT2,
     &           id%COMM, IERR )
         ELSEIF (id%MYID.EQ.MASTER) THEN
C     Receive
            CALL MPI_RECV( id%NB_SINGULAR_VALUES,
     &           1,
     &           MPI_INTEGER, ID_ROOT, TAG_ROOT1,
     &           id%COMM, STATUS, IERR )
            IF (associated(id%SINGULAR_VALUES))
     &           DEALLOCATE(id%SINGULAR_VALUES)
            NULLIFY(id%SINGULAR_VALUES)
            ALLOCATE(id%SINGULAR_VALUES(id%NB_SINGULAR_VALUES)
     &           , stat=ALLOCOK)
            IF (ALLOCOK.GT.0) THEN
               id%INFO(1)=-13
               id%INFO(2)= id%NB_SINGULAR_VALUES
               RETURN
            END IF
            CALL MPI_RECV( id%SINGULAR_VALUES(1),
     &           id%NB_SINGULAR_VALUES,
     &           MPI_DOUBLE_PRECISION, ID_ROOT, TAG_ROOT2,
     &           id%COMM, STATUS, IERR )
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_EXTRACT_SINGULAR_VALUES
      SUBROUTINE ZMUMPS_SET_NOMP_MAX(KEEP281, KEEP361,
     &                               N, NOMP_MAX)
!$    USE OMP_LIB
C
C  Purpose
C  =======
C     set NOMP_MAX from KEEP(281)
C     on output NOMP_MAX >=0
C
C  Parameters
C  ==========
C
      INTEGER, INTENT(IN)    :: KEEP281, KEEP361, N
      INTEGER, INTENT(OUT)   :: NOMP_MAX
C     
C     Local variables
C   
      INTEGER :: NOMP
C
C     out-of-range entries treated as -1
      NOMP_MAX= max(-1, KEEP281)
      NOMP    = 1
!$    NOMP = omp_get_max_threads()
      IF (NOMP_MAX.EQ.-1) THEN
C     automatic setting
       IF (N.LE.KEEP361) THEN
         NOMP_MAX = 0
         RETURN
       ENDIF
       IF (NOMP.GT.1) THEN
C       conservative because of memory allocation
        NOMP_MAX = min(NOMP, 10)
       ELSE
C       no multithreading and all parallel do suppressed
        NOMP_MAX = 0
       ENDIF
      ELSE
C     NOMP_MAX >=0
C     use provided value
       NOMP_MAX = min(NOMP_MAX, NOMP)
      ENDIF
C
      RETURN
      END SUBROUTINE ZMUMPS_SET_NOMP_MAX
