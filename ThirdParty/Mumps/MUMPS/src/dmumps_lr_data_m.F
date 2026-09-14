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
C     The DMUMPS_LR_DATA_M module is designed to manage low-rank data
C     associated to the fronts in a MUMPS instance.
C     The module was initially designed
C     to handle save low-rank panels in case of LDT type2 slaves 
C     has been revisited and generalized to be able to store BLR panels
C     of types 1 and 2, in both LU and LDLT cases, as well as
C     compressed contribution blocks. 
      MODULE DMUMPS_LR_DATA_M
      USE DMUMPS_LR_TYPE
      IMPLICIT NONE
C     This file will later contain a module
C     to manage Low rank information.
      PRIVATE
      PUBLIC :: DMUMPS_BLR_END_FRONT, DMUMPS_BLR_INIT_MODULE,
     &  DMUMPS_BLR_END_MODULE, DMUMPS_BLR_INIT_FRONT, 
     &  DMUMPS_BLR_SAVE_INIT,
     &  DMUMPS_BLR_SAVE_PANEL_LORU, DMUMPS_BLR_RETRIEVE_BEGS_BLR_L,
     &  DMUMPS_BLR_SAVE_BEGS_BLR_C, DMUMPS_BLR_RETRIEVE_BEGS_BLR_C,
     &  DMUMPS_BLR_DEC_AND_RETRIEVE_L, DMUMPS_BLR_RETRIEVE_PANEL_LORU,
     &  DMUMPS_BLR_DEC_AND_TRYFREE_L, DMUMPS_BLR_TRY_FREE_PANEL, 
     &  DMUMPS_BLR_FORCE_FREE_PANEL_L,
     &  DMUMPS_BLR_FREE_CB_LRB, DMUMPS_BLR_FREE_ALL_PANELS,
     &  DMUMPS_BLR_SAVE_CB_LRB,
     &  DMUMPS_BLR_RETRIEVE_CB_LRB, DMUMPS_BLR_RETRIEVE_BEGSBLR_STA,
     &  DMUMPS_BLR_SAVE_BEGS_BLR_DYN, DMUMPS_BLR_RETRIEVE_BEGSBLR_DYN,
     &  DMUMPS_BLR_RETRIEVE_NB_PANELS, DMUMPS_BLR_EMPTY_PANEL_LORU,
     &  DMUMPS_BLR_SAVE_NFS4FATHER, DMUMPS_BLR_RETRIEVE_NFS4FATHER,
     &  DMUMPS_BLR_SAVE_M_ARRAY, DMUMPS_BLR_RETRIEVE_M_ARRAY,
     &  DMUMPS_BLR_FREE_M_ARRAY
     &  , DMUMPS_BLR_STRUC_TO_MOD, DMUMPS_BLR_MOD_TO_STRUC, BLR_ARRAY
C         BLR_ARRAY public requires its types and subtypes to be public
     &  , BLR_STRUC_T, blr_panel_type, diag_block_type
     &  , DMUMPS_BLR_SAVE_DIAG_BLOCK, DMUMPS_BLR_RETRIEVE_DIAG_BLOCK
#if ! defined(NO_SAVE_RESTORE)
     &  , DMUMPS_SAVE_RESTORE_BLR
#endif      
C
C     --------------------------------------
C     This type is meant to store a panel (or the corresponding
C     part within a front) in BLR format
      TYPE blr_panel_type
         integer                 :: NB_ACCESSES_LEFT
         type(LRB_TYPE), pointer :: LRB_PANEL(:)
C        In case of type 1 node:
C          * both in case of L and U,
C            1+NPARTSCB <= size(LRB_PANEL) <= NPARTSASS+NPARTSCB
C        In case of type 2 master:
C          * U only: 1+NPARTSCB <= size(LRB_PANEL) <= NPARTSASS+NPARTSCB
C          * L only: 1 <= size(LRB_PANEL)  <= NPARTSASS
C        In case of type 2 slave:
C          * U: LRB_PANEL is not allocated
C          * L: size(LRB_PANEL) <= NPARTSASS
      END TYPE blr_panel_type
      TYPE diag_block_type
        DOUBLE PRECISION, POINTER   :: DIAG_BLOCK(:)
#if !defined(NOMIXHPC)
        TYPE(LRMIXED_TYPE) :: DIAG_STOC
        INTEGER :: DIAG_SIZE = 0
#endif
      END TYPE diag_block_type
C     --------------------------------------
C     --------------------------------------
C     This type is meant to store factors of level-1 nodes, of bands,
C     or of level-2 fronts in BLR format. This is simply a list of
C     L and U panels.
C    
C     The size of arrays PANELS_L / PANELS_U is the number of panels
C     and should be NPARTSASS
      TYPE BLR_STRUC_T
C        LOCICALS to indicate if we are in the case of 
C            Symmetric factorization (IsSYm)
C            Type 2 node (IsT2)
C            Slave panel (IsSlave)
         LOGICAL   :: IsSYM, IsT2, IsSLAVE
         TYPE(blr_panel_type), DIMENSION (:), POINTER  :: PANELS_L 
         TYPE(blr_panel_type), DIMENSION (:), POINTER  :: PANELS_U
         TYPE(LRB_TYPE), pointer :: CB_LRB(:,:)
         TYPE(diag_block_type), DIMENSION (:), POINTER :: DIAG_BLOCKS
C        BEGS_BLR_STATIC : 
C          Initial static grouping (after regrouping but before pivoting)
C          of the master fronts
         INTEGER, DIMENSION(:), POINTER :: BEGS_BLR_STATIC
C        BEGS_BLR_DYNAMIC :
C          contains the beginning number of each BLR
C          block at the end of the factorization. It includes the shift
C           due to numerical pivoting
         INTEGER, DIMENSION(:), POINTER :: BEGS_BLR_DYNAMIC 
C     Set of pointers to retreive BLR panels (desiged for Type 2 only (to be checked)) 
C     Rows of slaves
         INTEGER, DIMENSION(:), POINTER :: BEGS_BLR_L
C     columns of slaves (only used for LDLT type 2 slaves)
         INTEGER, DIMENSION(:), POINTER :: BEGS_BLR_COL
C        BEGS_BLR_COL: those BLR pointers are needed when CB block of the 
C                      slave is compressed. In this case only the block
C                      stucture of the CB block is used which starts at
C                      BEGS_BLR_COL(NB_PANELS+1)
C                      
C  NB_ACCESSES_INIT and NB_ACCESSES_LEFT:
C  =====================================
C    =============================================================
C    FOR BLR-LUA schemes BLR structure need be kept in all cases :
C    ===========
C
C       NB_ACCESSES_INIT --> huge or normal NB_ACCESSES_INIT +1 and this
C                         even for "first slave" when NSLAVES_PREC is zero.
C
C       For type2 nodes not so easy to keep all neccessary data 
C       without modifying the pattern of communication.
C       It has been decided to first implement LUA in type 1 cases
C       and validate the potential on this. 
C    =============================================================
C    ===========================
C    Case of TYPE 2 slaves LDLT
C    ===========================
C     Algorithm describing the 
C     possible status of panels and corresponding values of 
C     NB_ACCESSES_INIT and NB_ACCESSES_LEFT
C
C        After DMUMPS_BLR_INIT_FRONT:
C         If BLR PANELS will not be stored 
C               (only storage of BEGS_BLR_L and BEGS_BLR_COL)
C               then NB_ACCESSES_INIT = PANELS_NOTUSED 
C         else NB_ACCESSES_INIT >0 
C         endif
C    
C      NB_ACCESSES_LEFT (of a given panel number)
C      ================
C       Assumption DMUMPS_BLR_INIT_FRONT was called with 
C       and NB_ACCESSES_INIT > 0
C       
C       Each time DMUMPS_BLR_SAVE_PANEL(ipanel) is called to 
C       store a new panel NB_ACCESSES_LEFT <-- NB_ACCESSES_INIT
C       Each time DMUMPS_BLR_DEC_AND_RETRIEVE_L(ipanel) is called
C       NB_ACCESSES_LEFT of ipanel is decremented 
C       (but not decremented during DMUMPS_BLR_RETRIEVE_PANEL_LORU)
C
C       Each time DMUMPS_BLR_DEC_AND_TRYFREE is called
C       NB_ACCESSES_LEFT of ipanel is decremented  and 
C       and DMUMPS_BLR_TRY_FREE_PANEL(ipanel) is called
C       
C       IF when calling DMUMPS_BLR_TRY_FREE_PANEL(ipanel) 
C       NB_ACCESSES_LEFT of ipanel is zero then 
C       LRB_PANEL of ipanel is freed.
         INTEGER :: NB_ACCESSES_INIT ! Number of accesses to each panel
         INTEGER :: NB_PANELS ! Number of panels in FS blocks 
C        NB_PANELS: is the number of BLR panel in FS block on master.
C                   it is invariant of the numerical pivoting since
C                   non-eliminated variables will glue to current BLR block.
C                   However when it is recomputed (it should be suppressed
C                   for clearness, even if it works because position of 
C                   panels before the current one are never used) 
C                   at each BLR step it might increase.
C                   NB_PANELS is needed to find the starting of pointers
C                   in BEG_BLR_COL to retreive BLR structure of CB in case
C                   of compress CB.
C
C // PIV and COMPRESS_CB
C        NFS4FATHER used in sym indef and in case of COMPRESS_CB when
C             the father is of type 2,
C            The current node, referred to as the SON, computes
C            an estimation of the number of FS columns in its FATHER
C            present in the current and so called SON node.
C            NFS4FATHER is computed by the master of a SON node (T1 or T2).
C            The number of uneliminated variables of the SON 
C            is not yet know and thus not included in NFS4FATHER.
C            IF son is of type 2 then
C              NFS4FATHER is sent to the slaves of SON and stored 
C              in BLR structure during DMUMPS_PROCESS_DESC_BANDE
C            ELSE
C              NFS4FATHER is saved in  BLR structure during assembly
C            ENDIF
C            It is sent to its final value including NELIM variables
C            of the the son when M_ARRAY is saved 
C            (DMUMPS_BLR_SAVE_M_ARRAY)
         INTEGER   :: NFS4FATHER
C        M_ARRAY holds the maximum per fully-summed columns
C                computed before COMPRESS_CB
         DOUBLE PRECISION,  DIMENSION(:), POINTER :: M_ARRAY
C                   
      END TYPE BLR_STRUC_T
C     --------------------------------------
C  BLR_ARRAY(IWHANDLER), where IWHANDLER is stored in header 
C  of each front (at position IW(IOLDPS+XXF)), holds the
C  BLR structure stored by panel for the front.
C  BLR_ARRAY should be saved in MUMPS structure:
C     1/ either at the end of the facto phase with a routine 
C       of the module doing it.
C     2/ or storing directly this array in module. Then 
C        BLR_ARRAY should be made accessible at all level
C        of facto
C       
C     Validate with 1 and then do 2/ when making MUMPS thread safe ?
C 
      type(BLR_STRUC_T), POINTER, DIMENSION(:), SAVE :: 
     &                                          BLR_ARRAY => null()
C     IF during DMUMPS_BLR_INIT_MODULE the INITIAL_SIZE is zero
C     (BLR acttivated but no BLR fronts)
C
C     Derived datatype to include a BLR_ARRAY pointer,
      TYPE BLR_ARRAY_T
        type(BLR_STRUC_T), POINTER, DIMENSION(:) :: BLR_ARRAY
      END TYPE BLR_ARRAY_T
C
      INTEGER BLR_ARRAY_FREE, PANELS_NOTUSED, PANELS_FREED,
     &        NB_PANELS_NOTINIT, NFS4FATHER_NOTINIT
      PARAMETER (BLR_ARRAY_FREE=-9999, 
     &           PANELS_NOTUSED=-1111, PANELS_FREED=-2222, 
     &           NB_PANELS_NOTINIT=-3333, 
     &           NFS4FATHER_NOTINIT=-4444 )
#if ! defined(NO_SAVE_RESTORE)
      INCLUDE 'mumps_save_restore_modes.h'
#endif
      CONTAINS
      ! this routine only allocates the BLR_FRONTS
      ! array and should be called upon beginning of the factorization
      ! phase.
      SUBROUTINE DMUMPS_BLR_INIT_MODULE(INITIAL_SIZE, INFO
     &  )
C
C      Purpose:
C      =======
C
C      Module initialization
C
C      Arguments
C      =========
C
       INTEGER, INTENT(IN) :: INITIAL_SIZE
       INTEGER, INTENT(INOUT) :: INFO(2)
C
C      Local variables
C      ===============
C
       INTEGER :: I, IERR
C
       nullify(BLR_ARRAY)
       IF (INITIAL_SIZE.LE.0) RETURN
       ALLOCATE(BLR_ARRAY( INITIAL_SIZE ), stat=IERR)
       IF (IERR > 0 ) THEN
         INFO(1)=-13
         INFO(2)=INITIAL_SIZE
         RETURN
       ENDIF
       DO I=1, INITIAL_SIZE
C
C FIXME: might be elegant to test status for debugging
C        BLR_ARRAY(I)%INITFRONTDONE = .FALSE.
         NULLIFY(BLR_ARRAY(I)%PANELS_L)
         NULLIFY(BLR_ARRAY(I)%PANELS_U)
         NULLIFY(BLR_ARRAY(I)%CB_LRB)
         NULLIFY(BLR_ARRAY(I)%DIAG_BLOCKS)
         NULLIFY(BLR_ARRAY(I)%BEGS_BLR_STATIC)
         NULLIFY(BLR_ARRAY(I)%BEGS_BLR_DYNAMIC)
         BLR_ARRAY(I)%NB_ACCESSES_INIT = BLR_ARRAY_FREE
         BLR_ARRAY(I)%NB_PANELS        = NB_PANELS_NOTINIT
         NULLIFY(BLR_ARRAY(I)%BEGS_BLR_L)
         NULLIFY(BLR_ARRAY(I)%BEGS_BLR_COL)
         BLR_ARRAY(I)%NFS4FATHER       = NFS4FATHER_NOTINIT
         NULLIFY(BLR_ARRAY(I)%M_ARRAY)
       ENDDO
       RETURN
       END SUBROUTINE DMUMPS_BLR_INIT_MODULE
       SUBROUTINE DMUMPS_BLR_END_MODULE(INFO1, KEEP8, K34
     &             , LRSOLVE_ACT_OPT
     &            )
C
C      Purpose:
C      =======
C      Module final termination.
C      Called from DMUMPS_END_DRIVER or from DMUMPS_FAC_DRIVER
C
C      Arguments:
C      =========
C
       INTEGER, INTENT(IN) :: INFO1, K34
       LOGICAL, OPTIONAL, INTENT(IN) :: LRSOLVE_ACT_OPT
       INTEGER(8) :: KEEP8(150)
C
C      Local variables:
C      ===============
C
       INTEGER :: I, ILOOP
       IF (.NOT. associated(BLR_ARRAY)) THEN
C        could be not associated if DMUMPS_BLR_INIT_MODULE called
C        with INITIAL_SIZE==0 (BLR activated and no BLR fronts)
         RETURN
       ENDIF
       DO I=1, size(BLR_ARRAY)
C        ILOOP is a copy of I that is passed to DMUMPS_END_FRONT
C        to avoid I being overwritten       
         ILOOP= I
         IF (associated(BLR_ARRAY(I)%PANELS_L).OR.
     &       associated(BLR_ARRAY(I)%PANELS_U).OR.
     &       associated(BLR_ARRAY(I)%CB_LRB).OR.
     &       associated(BLR_ARRAY(I)%DIAG_BLOCKS)
     &      ) THEN
C
C          There are some BLR data to be freed:
C          factors (L, U, or D) or CB blocks
C
           IF (present(LRSOLVE_ACT_OPT)) THEN
             CALL DMUMPS_BLR_END_FRONT(ILOOP, INFO1, KEEP8, K34
     &                              , LRSOLVE_ACT_OPT
     &                              )
           ELSE
             CALL DMUMPS_BLR_END_FRONT(ILOOP, INFO1, KEEP8, K34 )
           ENDIF
         ENDIF
       ENDDO
       DEALLOCATE(BLR_ARRAY)
       NULLIFY(BLR_ARRAY)
       RETURN
       END SUBROUTINE DMUMPS_BLR_END_MODULE
       SUBROUTINE DMUMPS_BLR_MOD_TO_STRUC(id_BLRARRAY_ENCODING)
C
C  Purpose:
C  =======
C      From BLR_ARRAY in MODule to MUMPS STRUCture:
C      save/encode BLR_ARRAY pointer information (pointers only) in
C      MUMPS main structure.
C
C  Argument:
C  ========
#      if defined(MUMPS_NOF2003)
       CHARACTER, DIMENSION(:), POINTER :: id_BLRARRAY_ENCODING
#      else
       CHARACTER, DIMENSION(:), POINTER, intent(inout) ::
     &                                            id_BLRARRAY_ENCODING
#      endif
C
C  Local variables
C  ===============
C
C      Character array of arbitrary size 1
       CHARACTER :: CHAR_ARRAY(1)
       INTEGER :: CHAR_LENGTH, IERR
       TYPE(BLR_ARRAY_T) :: BLR_ARRAY_VAR
       IF (associated(id_BLRARRAY_ENCODING)) THEN
C     Should be unassociated for this to work
         WRITE(*,*) "Internal error 1 in DMUMPS_BLR_MOD_TO_STRUC"
         CALL MUMPS_ABORT()
       ENDIF
       BLR_ARRAY_VAR%BLR_ARRAY => BLR_ARRAY
       CHAR_LENGTH=size(transfer(BLR_ARRAY_VAR,CHAR_ARRAY))
       ALLOCATE(id_BLRARRAY_ENCODING(CHAR_LENGTH), stat=IERR)
       IF (IERR > 0 ) THEN
          WRITE(*,*) "Allocation error in DMUMPS_BLR_MOD_TO_STRUC"
          CALL MUMPS_ABORT()
       ENDIF
C      Fill with derived datatype
       id_BLRARRAY_ENCODING=transfer(BLR_ARRAY_VAR,CHAR_ARRAY)
C      BLR_ARRAY should now only be accessed again
C      after DMUMPS_BLR_STRUC_TO_MOD is called
       NULLIFY(BLR_ARRAY)
       RETURN
       END SUBROUTINE DMUMPS_BLR_MOD_TO_STRUC
C
       SUBROUTINE DMUMPS_BLR_STRUC_TO_MOD(id_BLRARRAY_ENCODING)
C
C  Purpose:
C  =======
C      From to MUMPS STRUCture to variable BLR_ARRAY in MODule:
C      set/decode BLR_ARRAY pointer information 
C      (pointers only, allocation was done earlier) 
C      by decoding id_BLRARRAY_ENCODING characters
C      from main MUMPS structure.
C
C  Arguments
C  =========
C
#      if defined(MUMPS_NOF2003)
       CHARACTER, DIMENSION(:), POINTER :: id_BLRARRAY_ENCODING
#      else
       CHARACTER, DIMENSION(:), POINTER, intent(inout) ::
     &                                            id_BLRARRAY_ENCODING
#      endif
C
C  Local variables
C  ===============
C
       TYPE (BLR_ARRAY_T) :: BLR_ARRAY_VAR
C
       IF (.NOT.associated(id_BLRARRAY_ENCODING)) THEN
         WRITE(*,*) "Internal error 1 in DMUMPS_BLR_STRUC_TO_MOD"
       ENDIF
       BLR_ARRAY_VAR = transfer(id_BLRARRAY_ENCODING,BLR_ARRAY_VAR)
       BLR_ARRAY => BLR_ARRAY_VAR%BLR_ARRAY
C      Module is now responsible for BLR_ARRAY
       DEALLOCATE(id_BLRARRAY_ENCODING)
       NULLIFY(id_BLRARRAY_ENCODING)
       RETURN
       END SUBROUTINE DMUMPS_BLR_STRUC_TO_MOD
       SUBROUTINE DMUMPS_BLR_INIT_FRONT(IWHANDLER, 
     &                       INFO, MTK405)
       USE MUMPS_FRONT_DATA_MGT_M, ONLY : MUMPS_FDM_START_IDX
!$    USE OMP_LIB
C
C      Purpose:
C      =======
C
C      Initialize a BLR structure for a front.
C
C      Arguments:
C      =========
C
       INTEGER, INTENT(INOUT) :: IWHANDLER, INFO(2)
       INTEGER, INTENT(IN), OPTIONAL    :: MTK405
C      (MTK405 is 1 if DMUMPS_BLR_INIT_FRONT can
C       be called by several threads simultaneously)
C
C      Local variables:
C      ===============
C
       TYPE(BLR_STRUC_T), POINTER, DIMENSION(:) :: BLR_ARRAY_TMP
       INTEGER :: OLD_SIZE, NEW_SIZE
       INTEGER :: I
       INTEGER :: IERR
       LOGICAL :: NEEDS_THREAD_SAFETY
C
       NEEDS_THREAD_SAFETY = .FALSE.
       IF (present(MTK405)) THEN
         IF (MTK405 .EQ. 1 ) THEN
           NEEDS_THREAD_SAFETY = .TRUE.
         ENDIF
       ENDIF
       IF ( NEEDS_THREAD_SAFETY ) THEN
!$OMP CRITICAL(critical_blr_idx)
         CALL MUMPS_FDM_START_IDX('F', 'INITF', IWHANDLER, INFO)
!$OMP END CRITICAL(critical_blr_idx)
       ELSE
         CALL MUMPS_FDM_START_IDX('F', 'INITF', IWHANDLER, INFO)
       ENDIF
       IF (IWHANDLER > size(BLR_ARRAY)) THEN
C        Reallocate in a bigger array
C        Because front data management module increases the size
C        of its internal arrays independently of ours and can be
C        called for other purposes than MAPROW (buffer management
C        of type 2 nodes...), IWHANDLER might (in exceptional cases)
C        be larger than what we may expect. This is why we use a max.
         OLD_SIZE = size(BLR_ARRAY)
         NEW_SIZE = max( (OLD_SIZE * 3) / 2 + 1, IWHANDLER)
C        
         ALLOCATE(BLR_ARRAY_TMP(NEW_SIZE),stat=IERR)
         IF (IERR.GT.0) THEN
           INFO(1)=-13
           INFO(2)=NEW_SIZE
           GOTO 500
         ENDIF
         DO I=1, OLD_SIZE
           BLR_ARRAY_TMP(I)=BLR_ARRAY(I)
         ENDDO
C        Similar to what is done in MUMPS_FMRD_INIT
         DO I=OLD_SIZE+1, NEW_SIZE
           NULLIFY(BLR_ARRAY_TMP(I)%PANELS_L)
           NULLIFY(BLR_ARRAY_TMP(I)%PANELS_U)
           NULLIFY(BLR_ARRAY_TMP(I)%CB_LRB)
           NULLIFY(BLR_ARRAY_TMP(I)%DIAG_BLOCKS)
           NULLIFY(BLR_ARRAY_TMP(I)%BEGS_BLR_STATIC)
           NULLIFY(BLR_ARRAY_TMP(I)%BEGS_BLR_DYNAMIC)
           BLR_ARRAY_TMP(I)%NB_ACCESSES_INIT = BLR_ARRAY_FREE
           BLR_ARRAY_TMP(I)%NB_PANELS        = NB_PANELS_NOTINIT
           NULLIFY(BLR_ARRAY_TMP(I)%BEGS_BLR_L)
           NULLIFY(BLR_ARRAY_TMP(I)%BEGS_BLR_COL)
           BLR_ARRAY_TMP(I)%NFS4FATHER       = NFS4FATHER_NOTINIT
           NULLIFY(BLR_ARRAY_TMP(I)%M_ARRAY)
         ENDDO
         DEALLOCATE(BLR_ARRAY)
         BLR_ARRAY => BLR_ARRAY_TMP
         NULLIFY(BLR_ARRAY_TMP)
  500   CONTINUE
       ENDIF
       RETURN
       END SUBROUTINE DMUMPS_BLR_INIT_FRONT
C
C
       SUBROUTINE DMUMPS_BLR_SAVE_INIT(IWHANDLER, 
     &                       IsSYM, IsT2, IsSLAVE, 
     &                       NB_PANELS,
     &                       BEGS_BLR_L, BEGS_BLR_COL, 
     &                       NB_ACCESSES_INIT, INFO)
C
C      Purpose:
C      =======
C
C      Initial SAVE of internal structures of a front.
C      Note that DMUMPS_BLR_INIT_FRONT should have been called
C                and thus IWHANDLER set
C
C
C      Arguments:
C      =========
C
!     The number of BLR col panels associated
!     to fully summed part of front (npartsass)
!     (not used if NB_ACCESSES_INIT==0)
       LOGICAL, INTENT(IN) :: IsSYM, IsT2, IsSLAVE
       INTEGER, INTENT(IN) :: NB_PANELS, IWHANDLER
       INTEGER, INTENT(INOUT) :: INFO(2)
!     Nb of accesses to each panel for updates
       INTEGER, INTENT(IN) :: NB_ACCESSES_INIT 
       INTEGER, INTENT(IN), DIMENSION(:) :: BEGS_BLR_L
C FIXME next should be optional and then tested as present()
       INTEGER,  DIMENSION(:), POINTER :: BEGS_BLR_COL
C
C      Local variables:
C      ===============
C
       INTEGER :: I
       INTEGER :: IERR
C
       IF (NB_PANELS.EQ.0) THEN
        WRITE(6,*) " Internal error 1 in DMUMPS_BLR_SAVE_INIT ", 
     &    NB_PANELS
       ENDIF
       IF (IWHANDLER .LE.0 ) THEN
        WRITE(6,*) " Internal error 2 in DMUMPS_BLR_SAVE_INIT ", 
     &     IWHANDLER
       ENDIF
C
       IF (associated(BEGS_BLR_COL)) THEN
C       case with BEGS_BLR_COL provided
         ALLOCATE(
     &         BLR_ARRAY(IWHANDLER)%BEGS_BLR_COL(size(BEGS_BLR_COL)),
     &         stat=IERR)
         IF (IERR .GT. 0) THEN
           INFO(1)=-13
           INFO(2)=size(BEGS_BLR_COL)
           RETURN
         ENDIF
       ENDIF
       IF (NB_ACCESSES_INIT.EQ.0) THEN
C       NO NEED TO allocate PANELS_L(NB_PANELS) and PANELS_U
        NULLIFY(BLR_ARRAY(IWHANDLER)%PANELS_L)
        NULLIFY(BLR_ARRAY(IWHANDLER)%PANELS_U)
        NULLIFY(BLR_ARRAY(IWHANDLER)%CB_LRB)
        NULLIFY(BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS)
        ALLOCATE(
     &         BLR_ARRAY(IWHANDLER)%BEGS_BLR_L(size(BEGS_BLR_L)),
     &         BLR_ARRAY(IWHANDLER)%BEGS_BLR_STATIC(size(BEGS_BLR_L)),
     &         BLR_ARRAY(IWHANDLER)%BEGS_BLR_DYNAMIC(size(BEGS_BLR_L)),
     &         stat=IERR)
        IF (IERR .GT. 0) THEN
           INFO(1)=-13
           INFO(2)=3*size(BEGS_BLR_L)
           RETURN
        ENDIF
       ELSE
C        write(6,*) " ***LR_DATA_M:DMUMPS_BLR_SAVE_INIT: NB_PANELS=", 
C     &    NB_PANELS
        IF (IsSYM) THEN
C        LDLT 
         ALLOCATE(BLR_ARRAY(IWHANDLER)%PANELS_L(NB_PANELS),
     &          BLR_ARRAY(IWHANDLER)%BEGS_BLR_L(size(BEGS_BLR_L)),
     &          BLR_ARRAY(IWHANDLER)%BEGS_BLR_STATIC(size(BEGS_BLR_L)),
     &          BLR_ARRAY(IWHANDLER)%BEGS_BLR_DYNAMIC(size(BEGS_BLR_L)),
     &          stat=IERR)
        ELSE 
C        LU  (both PANELS_L and PANELS_U are needed)
         ALLOCATE(BLR_ARRAY(IWHANDLER)%PANELS_L(NB_PANELS),
     &          BLR_ARRAY(IWHANDLER)%PANELS_U(NB_PANELS),
     &          BLR_ARRAY(IWHANDLER)%BEGS_BLR_STATIC(size(BEGS_BLR_L)),
     &          BLR_ARRAY(IWHANDLER)%BEGS_BLR_DYNAMIC(size(BEGS_BLR_L)),
     &          BLR_ARRAY(IWHANDLER)%BEGS_BLR_L(size(BEGS_BLR_L)),
     &          stat=IERR)
        ENDIF
        IF (IERR .GT. 0) THEN
           INFO(1)=-13
           IF (IsSYM) THEN
              INFO(2)=NB_PANELS+3*size(BEGS_BLR_L)
           ELSE
              INFO(2)=NB_PANELS+NB_PANELS+3*size(BEGS_BLR_L)
           ENDIF
           RETURN
        ENDIF
        IF (.NOT.IsSLAVE) THEN
C         Diagonal blocks (needed in all cases except slave band)          
          ALLOCATE(BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS(NB_PANELS),
     &          stat=IERR)
          IF (IERR .GT. 0) THEN
            INFO(1)=-13
            INFO(2)=NB_PANELS
            RETURN
          ENDIF
        ENDIF
        DO I=1,NB_PANELS
           NULLIFY(BLR_ARRAY(IWHANDLER)%PANELS_L(I)%LRB_PANEL)
           IF (.NOT.IsSYM) THEN
C            LU : nullify also PANELS_U
             NULLIFY(BLR_ARRAY(IWHANDLER)%PANELS_U(I)%LRB_PANEL)
           ENDIF
           IF (.NOT.IsSLAVE) THEN
             NULLIFY(BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS(I)%DIAG_BLOCK)
           ENDIF
        ENDDO
       ENDIF
       BLR_ARRAY(IWHANDLER)%IsSYM        = IsSYM
       BLR_ARRAY(IWHANDLER)%IsT2         = IsT2
       BLR_ARRAY(IWHANDLER)%IsSLAVE      = IsSLAVE
       BLR_ARRAY(IWHANDLER)%NB_PANELS    = NB_PANELS
       BLR_ARRAY(IWHANDLER)%BEGS_BLR_L   = BEGS_BLR_L
       BLR_ARRAY(IWHANDLER)%BEGS_BLR_STATIC   = BEGS_BLR_L
       BLR_ARRAY(IWHANDLER)%BEGS_BLR_DYNAMIC   = -999991
       IF (NB_ACCESSES_INIT.EQ.0) THEN
         BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT = PANELS_NOTUSED
       ELSE
         BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT = NB_ACCESSES_INIT
       ENDIF
       IF (associated(BEGS_BLR_COL)) THEN
          DO I=1,size(BEGS_BLR_COL)
             BLR_ARRAY(IWHANDLER)%BEGS_BLR_COL(I) = BEGS_BLR_COL(I)
          ENDDO
       ELSE
         NULLIFY( BLR_ARRAY(IWHANDLER)%BEGS_BLR_COL )
       ENDIF
       RETURN
       END SUBROUTINE DMUMPS_BLR_SAVE_INIT
C
C
C
       SUBROUTINE DMUMPS_BLR_END_FRONT(IWHANDLER, INFO1, KEEP8, K34
     &             , LRSOLVE_ACT_OPT, MTK405 )
       USE MUMPS_FRONT_DATA_MGT_M, ONLY : MUMPS_FDM_END_IDX
       use iso_c_binding
C
C      Purpose:
C      =======
C
C      Front final termination.
C
C      Arguments:
C      =========
C
       INTEGER, INTENT(INOUT) :: IWHANDLER
       INTEGER, INTENT(IN) :: INFO1
       INTEGER(8) :: KEEP8(150)
       INTEGER, INTENT(IN) :: K34
       LOGICAL, OPTIONAL, INTENT(IN) :: LRSOLVE_ACT_OPT
       INTEGER, OPTIONAL, INTENT(IN) :: MTK405
C
C      Local variables:
C      ===============
C
       INTEGER :: IPANEL, JPANEL
       INTEGER(8) :: MEM_FREED
       INTEGER :: IDUMMY, JDUMMY
       TYPE(blr_panel_type), POINTER  :: THEPANEL
       LOGICAL :: LRSOLVE_ACT, NEEDS_THREAD_SAFETY
       TYPE(diag_block_type), POINTER  :: THEBLOCK
C
C
C      ---------------------------
C      Look at optional parameters
C      ---------------------------
       LRSOLVE_ACT = .FALSE.
       IF (present(LRSOLVE_ACT_OPT)) LRSOLVE_ACT = LRSOLVE_ACT_OPT
C      enable calling DMUMPS_BLR_TRY_FREE_PANEL when module was not called
       IF (IWHANDLER.LE.0) THEN 
        RETURN
       ENDIF
       NEEDS_THREAD_SAFETY = .FALSE.
       IF (present(MTK405)) THEN
         IF (MTK405 .EQ. 1 ) THEN
           NEEDS_THREAD_SAFETY = .TRUE.
         ENDIF
       ENDIF
       IF (IWHANDLER .GT. size(BLR_ARRAY)) THEN
C       Not necessarily and error if data management module increases the size
C       of its internal arrays independently of ours 
C           WRITE(*,*) " Internal Error 1 in MUMPS_BLR_END_FRONT ",
C       &              IWHANDLER
C           CALL MUMPS_ABORT()
        RETURN
       END IF
       IF (BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT.EQ.BLR_ARRAY_FREE)
     &      RETURN  ! was called while INIT_FRONT NOT called
                    ! case of wrong "wrong" call to DMUMPS_BLR_END_FRONT
                    ! (IWHANDLER was however set for another purpose!)
       IF (BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT.NE.
     &      PANELS_NOTUSED) THEN
        DO IPANEL = 1, size(BLR_ARRAY(IWHANDLER)%PANELS_L)
           THEPANEL => BLR_ARRAY(IWHANDLER)%PANELS_L(IPANEL)
           IF (associated(THEPANEL%LRB_PANEL)) THEN
             IF (INFO1 .GE. 0 
     &            .AND..NOT.LRSOLVE_ACT
     &           ) THEN
C              All panels should be freed earlier except in case of
C              error or in case of LRSOLVE                
               WRITE(*,*) " Internal Error 2a in MUMPS_BLR_END_FRONT ",
     &          IWHANDLER, "NB_ACCESSES_INIT=", 
     &          BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT, 
     &          "Pointer to panel number ",IPANEL," still associated",
     &          " NB_ACCESSES_LEFT= ",THEPANEL%NB_ACCESSES_LEFT
               CALL MUMPS_ABORT()
             ELSE
               CALL DEALLOC_BLR_PANEL(THEPANEL%LRB_PANEL,
     &          size(THEPANEL%LRB_PANEL), KEEP8, K34)
               THEPANEL%NB_ACCESSES_LEFT = PANELS_FREED
             ENDIF
             DEALLOCATE(THEPANEL%LRB_PANEL)
             NULLIFY(THEPANEL%LRB_PANEL)
           ENDIF
        ENDDO
        IF (associated(BLR_ARRAY(IWHANDLER)%PANELS_L)) THEN
         DEALLOCATE(BLR_ARRAY(IWHANDLER)%PANELS_L)
         NULLIFY(BLR_ARRAY(IWHANDLER)%PANELS_L)
        ENDIF
C
        IF (.NOT.BLR_ARRAY(IWHANDLER)%IsSYM) THEN
C        PANELS_U should also be processed
         DO IPANEL = 1, size(BLR_ARRAY(IWHANDLER)%PANELS_U)
           THEPANEL => BLR_ARRAY(IWHANDLER)%PANELS_U(IPANEL)
           IF (associated(THEPANEL%LRB_PANEL)) THEN
             IF (INFO1 .GE. 0 
     &            .AND..NOT.LRSOLVE_ACT
     &          ) THEN
C            All panels should be freed earlier except in case of error
C            or in case of LRSOLVE                
             WRITE(*,*) " Internal Error 2b in MUMPS_BLR_END_FRONT ",
     &          IWHANDLER, "NB_ACCESSES_INIT=", 
     &          BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT, 
     &          "Pointer to panel number ",IPANEL," still associated"
             CALL MUMPS_ABORT()
             ELSE
               CALL DEALLOC_BLR_PANEL(THEPANEL%LRB_PANEL,
     &          size(THEPANEL%LRB_PANEL), KEEP8, K34)
               THEPANEL%NB_ACCESSES_LEFT = PANELS_FREED 
             ENDIF
             DEALLOCATE(THEPANEL%LRB_PANEL)
             NULLIFY(THEPANEL%LRB_PANEL)
           ENDIF
         ENDDO
         IF (associated(BLR_ARRAY(IWHANDLER)%PANELS_U)) THEN
          DEALLOCATE(BLR_ARRAY(IWHANDLER)%PANELS_U)
          NULLIFY(BLR_ARRAY(IWHANDLER)%PANELS_U)
         ENDIF
        ENDIF   ! case of LU and PANELS_U to process
        IF (.NOT.BLR_ARRAY(IWHANDLER)%IsSLAVE) THEN
C         DIAG_BLOCKS should also be processed
          MEM_FREED = 0_8
          DO IPANEL = 1, size(BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS)
            THEBLOCK => BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS(IPANEL)
            IF (associated(THEBLOCK%DIAG_BLOCK)) THEN
              IF (INFO1 .GE. 0 
     &               .AND..NOT.LRSOLVE_ACT
     &           ) THEN
C               All panels should be freed earlier except in case of error
C               or in case of LRSOLVE                
                WRITE(*,*) " Internal Error 3 in MUMPS_BLR_END_FRONT ",
     &             IWHANDLER, "NB_ACCESSES_INIT=", 
     &             BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT, 
     &             "Pointer to panel number ",IPANEL," still associated"
                CALL MUMPS_ABORT()
              ENDIF
              MEM_FREED = MEM_FREED + int(size(THEBLOCK%DIAG_BLOCK),8)
              DEALLOCATE    (THEBLOCK%DIAG_BLOCK)
              NULLIFY       (THEBLOCK%DIAG_BLOCK)
            ENDIF
C{
#if !defined(NOMIXHPC)
C{
            IF (associated(THEBLOCK%DIAG_STOC%VAL)) THEN
              IF (INFO1 .GE. 0 
     &               .AND..NOT.LRSOLVE_ACT
     &           ) THEN
C               All panels should be freed earlier except in case of error
C               or in case of LRSOLVE                
                WRITE(*,*) " Internal Error 3 in MUMPS_BLR_END_FRONT ",
     &             IWHANDLER, "NB_ACCESSES_INIT=", 
     &             BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT, 
     &             "Pointer to panel number ",IPANEL," still associated"
                CALL MUMPS_ABORT()
              ENDIF
              MEM_FREED = MEM_FREED 
     &           + int(1 + (size(THEBLOCK%DIAG_STOC%VAL)*
     &                      K34 - 1 )/COSTS_PREC(1),8)
              DEALLOCATE(THEBLOCK%DIAG_STOC%VAL)
              NULLIFY   (THEBLOCK%DIAG_STOC%VAL)
            ENDIF
C} MIXHPC
#endif
C} LRMIXED
          ENDDO
          IF ( MEM_FREED .GT. 0_8 ) THEN
            CALL MUMPS_DM_FAC_UPD_DYN_MEMCNTS(-MEM_FREED,
     &        NEEDS_THREAD_SAFETY, KEEP8, 
     &        IDUMMY, JDUMMY,  ! IFLAG, IERROR not used for deallocate
     &        .TRUE., .TRUE.)
          ENDIF
          IF (associated(BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS)) THEN
            DEALLOCATE(BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS)
            NULLIFY(BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS)
          ENDIF
        ENDIF   ! case of DIAG_BLOCKS to process
        IF (.NOT.BLR_ARRAY(IWHANDLER)%IsT2.OR.
     &           BLR_ARRAY(IWHANDLER)%IsSLAVE) THEN
C         CB_LRB should also be processed
          IF (associated(BLR_ARRAY(IWHANDLER)%CB_LRB)) THEN
             IF (INFO1 .GE. 0) THEN
C            In current version, all panels should
C            be freed earlier except in case of error.
             WRITE(*,*) 
     &          " Internal Error 4 in MUMPS_BLR_END_FRONT ",
     &          IWHANDLER, "CB block still associated", 
     &          BLR_ARRAY(IWHANDLER)%IsT2,
     &          BLR_ARRAY(IWHANDLER)%IsSLAVE
             CALL MUMPS_ABORT()
            ELSE
C             Case of error (memory need be freed)
C             (Probably called from DMUMPS_BLR_END_MODULE)
              DO IPANEL = 1, size(BLR_ARRAY(IWHANDLER)%CB_LRB,1)
              DO JPANEL = 1, size(BLR_ARRAY(IWHANDLER)%CB_LRB,2)
                CALL DEALLOC_LRB(
     &                 BLR_ARRAY(IWHANDLER)%CB_LRB(IPANEL,JPANEL),
     &                 KEEP8, K34)
              ENDDO
              ENDDO
              DEALLOCATE(BLR_ARRAY(IWHANDLER)%CB_LRB)
              NULLIFY(BLR_ARRAY(IWHANDLER)%CB_LRB)
            ENDIF
          ENDIF
        ENDIF   ! case of CB_LRB to process
C
       ENDIF  ! PANELS used
C      ------------------------------------
C      BEGS_BLR_STATIC, BEGS_BLR_DYNAMIC,
C      BEGS_BLR_L should always be freed,
C      "if (associated...)" only for safety
C      ------------------------------------
       IF (associated(BLR_ARRAY(IWHANDLER)%BEGS_BLR_STATIC)) THEN
C        Should always be initialized
         DEALLOCATE(BLR_ARRAY(IWHANDLER)%BEGS_BLR_STATIC)
         NULLIFY(BLR_ARRAY(IWHANDLER)%BEGS_BLR_STATIC)
       ENDIF
       IF (associated(BLR_ARRAY(IWHANDLER)%BEGS_BLR_DYNAMIC)) THEN
C        Should always be initialized
         DEALLOCATE(BLR_ARRAY(IWHANDLER)%BEGS_BLR_DYNAMIC)
         NULLIFY(BLR_ARRAY(IWHANDLER)%BEGS_BLR_DYNAMIC)
       ENDIF
       IF (associated(BLR_ARRAY(IWHANDLER)%BEGS_BLR_L)) THEN
C        Should always enter
         DEALLOCATE(BLR_ARRAY(IWHANDLER)%BEGS_BLR_L)
         NULLIFY(BLR_ARRAY(IWHANDLER)%BEGS_BLR_L)
       ENDIF
C      Check if BEGS_BLR_COL is associated
       IF (associated(BLR_ARRAY(IWHANDLER)%BEGS_BLR_COL)) THEN
          DEALLOCATE(BLR_ARRAY(IWHANDLER)%BEGS_BLR_COL)
          NULLIFY(BLR_ARRAY(IWHANDLER)%BEGS_BLR_COL)
       ENDIF
       BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT = BLR_ARRAY_FREE
       BLR_ARRAY(IWHANDLER)%NB_PANELS        = NB_PANELS_NOTINIT
       BLR_ARRAY(IWHANDLER)%NFS4FATHER       = NFS4FATHER_NOTINIT
       IF (associated(BLR_ARRAY(IWHANDLER)%M_ARRAY)) THEN
          DEALLOCATE(BLR_ARRAY(IWHANDLER)%M_ARRAY)
          NULLIFY(BLR_ARRAY(IWHANDLER)%M_ARRAY)
       ENDIF
       IF (NEEDS_THREAD_SAFETY) THEN
!$OMP    CRITICAL(critical_blr_idx)
         CALL MUMPS_FDM_END_IDX('F', 'ENDF', IWHANDLER)
!$OMP    END CRITICAL(critical_blr_idx)
       ELSE
         CALL MUMPS_FDM_END_IDX('F', 'ENDF', IWHANDLER)
       ENDIF
       RETURN
       END SUBROUTINE DMUMPS_BLR_END_FRONT
C
C
C
       SUBROUTINE DMUMPS_BLR_SAVE_PANEL_LORU (
     &    IWHANDLER, LORU, IPANEL, LRB_PANEL, NB_ACCESSES_INIT_IN )
C
C      Purpose:
C      =======
C
C      Save a new L or U panel. DMUMPS_BLR_INIT_FRONT has already been called.
C      We only point to the data passed, without a copy.
C 
C      If NB_ACCESSES_INIT_IN > 0 then
C        THEPANEL%NB_ACCESSES_LEFT = NB_ACCESSES_INIT_IN
C      else
C        THEPANEL%NB_ACCESSES_LEFT = 
C     &           BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT
C      endif
C      thus NB_ACCESSES_INIT_IN <=0 => NB_ACCESSES_INIT_IN not used
C
C      Arguments:
C      =========
C
       type(LRB_TYPE), DIMENSION(:), pointer  :: LRB_PANEL
C      The declaration below would work but the one above avoids the risk
C      of passing local data (e.g., an ALLOCATABLE array) with local scope
C      type(LRB_TYPE), DIMENSION(:), target  :: LRB_PANEL
       INTEGER, INTENT(IN) :: IWHANDLER, IPANEL, NB_ACCESSES_INIT_IN
       INTEGER, INTENT(IN) :: LORU
C      LORU :
C             0 = L Panel
C             1 = U Panel
C             others  : U panel
C
C      Local variables
C      ===============
C
       TYPE(blr_panel_type), POINTER  :: THEPANEL
C
       IF ( IWHANDLER .GT. size(BLR_ARRAY) .OR. IWHANDLER .LE. 0 ) THEN
C        Since DMUMPS_BLR_INIT_FRONT has been called, BLR_ARRAY
C        should be large enough
         WRITE(*,*) "Internal error 1 in DMUMPS_BLR_SAVE_PANEL_LORU"
         CALL MUMPS_ABORT()
       ENDIF
       IF (LORU.EQ.0) THEN
         THEPANEL => BLR_ARRAY(IWHANDLER)%PANELS_L(IPANEL)
       ELSE 
         THEPANEL => BLR_ARRAY(IWHANDLER)%PANELS_U(IPANEL)
       ENDIF
C      Number of futur accesses to this panel
       IF (NB_ACCESSES_INIT_IN.GT.0) THEN
         THEPANEL%NB_ACCESSES_LEFT = NB_ACCESSES_INIT_IN
       ELSE
         THEPANEL%NB_ACCESSES_LEFT = 
     &           BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT
       ENDIF
       THEPANEL%LRB_PANEL => LRB_PANEL
C
       RETURN
       END SUBROUTINE DMUMPS_BLR_SAVE_PANEL_LORU
C
       SUBROUTINE DMUMPS_BLR_SAVE_CB_LRB (
     &    IWHANDLER, CB_LRB )
C
C      Purpose:
C      =======
C
C      Save a new CB LR block.
C
C      Arguments:
C      =========
C
#if defined(MUMPS_NOF2003)
       TYPE(LRB_TYPE), POINTER :: CB_LRB(:,:)
#else       
       TYPE(LRB_TYPE), POINTER, INTENT(IN) :: CB_LRB(:,:)
#endif       
       INTEGER, INTENT(IN) :: IWHANDLER
C
       IF ( IWHANDLER .GT. size(BLR_ARRAY) .OR. IWHANDLER .LE. 0 ) THEN
C        Since DMUMPS_BLR_INIT_FRONT has been called, BLR_ARRAY
C        should be large enough
         WRITE(*,*) "Internal error 1 in DMUMPS_BLR_SAVE_CB_LRB"
         CALL MUMPS_ABORT()
       ENDIF
       BLR_ARRAY(IWHANDLER)%CB_LRB => CB_LRB
C
       RETURN
       END SUBROUTINE DMUMPS_BLR_SAVE_CB_LRB
C
       SUBROUTINE DMUMPS_BLR_SAVE_DIAG_BLOCK (
     &    IWHANDLER, IPANEL, D, KEEP34 )
         use iso_c_binding
C
C      Purpose:
C      =======
C
C      Save a the diagonal block of panel IPANEL 
C
C      Arguments:
C      =========
C
       DOUBLE PRECISION,POINTER :: D(:)
       INTEGER, INTENT(IN) :: IWHANDLER, IPANEL
       INTEGER, INTENT(IN) :: KEEP34
C
       IF ( IWHANDLER .GT. size(BLR_ARRAY) .OR. IWHANDLER .LE. 0 ) THEN
C        Since DMUMPS_BLR_INIT_FRONT has been called, BLR_ARRAY
C        should be large enough
         WRITE(*,*) "Internal error 1 in DMUMPS_BLR_SAVE_DIAG_BLOCK"
         CALL MUMPS_ABORT()
       ENDIF
       IF ( BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT.LT.0) THEN
         !  PANELS_NOTUSED or BLR_ARRAY_FREE
         WRITE(*,*) "Internal error 2 in DMUMPS_BLR_SAVE_DIAG_BLOCK"
         CALL MUMPS_ABORT()
       ENDIF
C{
#if !defined(NOMIXHPC)
C{
       BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS(IPANEL)%DIAG_SIZE = size(D)
       IF (.not.associated(
     &       BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS(IPANEL)%DIAG_BLOCK)) THEN
          if (PREC_DIAG_FR.NE.PREC_FAC) then
            allocate(BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS(IPANEL)%DIAG_STOC
     &         %VAL((COSTS_PREC(PREC_DIAG_FR)*size(D)+KEEP34-1)/KEEP34))
            CALL DMUMPS_LRMIXED_CAST(
     &        D(1),
     &        BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS(IPANEL)%DIAG_STOC%VAL(1),
     &        size(D), 1, size(D), size(D), PREC_FAC, PREC_DIAG_FR)
            deallocate(D)
            nullify(D)
            nullify(BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS(IPANEL)%DIAG_BLOCK)
          else
            BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS(IPANEL)%DIAG_BLOCK => D
            nullify(BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS(IPANEL)
     &                                %DIAG_STOC%VAL)
          end if
       ELSE
C} MIXHPC
#endif
C} LRMIXED
          BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS(IPANEL)%DIAG_BLOCK => D
#if !defined(NOMIXHPC)
          nullify(BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS(IPANEL)
     &                                %DIAG_STOC%VAL)
       ENDIF
#endif
       RETURN
       END SUBROUTINE DMUMPS_BLR_SAVE_DIAG_BLOCK
       SUBROUTINE DMUMPS_BLR_SAVE_BEGS_BLR_C (
     &    IWHANDLER, BEGS_BLR_COL, INFO)
C
C      Purpose:
C      =======
C
C      Save BEGS_BLR_COL. DMUMPS_BLR_INIT_FRONT has already been called.
C
C      Arguments:
C      =========
C
       INTEGER, DIMENSION(:), POINTER :: BEGS_BLR_COL
       INTEGER, INTENT(IN) :: IWHANDLER
       INTEGER, INTENT(INOUT) :: INFO(2)
C
C      Local variables:
C      ===============
C
       INTEGER :: I, IERR
       IF ( IWHANDLER .GT. size(BLR_ARRAY) .OR. IWHANDLER .LE. 0 ) THEN
C        Since DMUMPS_BLR_INIT_FRONT has been called, BLR_ARRAY
C        should be large enough
          WRITE(*,*) "Internal error 1 in DMUMPS_BLR_SAVE_BEGS_BLR_C"
          CALL MUMPS_ABORT()
       ENDIF
       IF ( BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT.LT.0) THEN
         !  PANELS_NOTUSED or BLR_ARRAY_FREE
          WRITE(*,*) "Internal error 2 in DMUMPS_BLR_SAVE_BEGS_BLR_C"
          CALL MUMPS_ABORT()
       ENDIF
       ALLOCATE(BLR_ARRAY(IWHANDLER)%BEGS_BLR_COL(size(BEGS_BLR_COL)),
     &          stat=IERR)
       IF (IERR > 0 ) THEN
         INFO(1)=-13
         INFO(2)=size(BEGS_BLR_COL)
         RETURN
       ENDIF
       DO I=1,size(BEGS_BLR_COL)
          BLR_ARRAY(IWHANDLER)%BEGS_BLR_COL(I) = BEGS_BLR_COL(I)
       ENDDO
C
       RETURN
       END SUBROUTINE DMUMPS_BLR_SAVE_BEGS_BLR_C
C
       SUBROUTINE DMUMPS_BLR_SAVE_BEGS_BLR_DYN (
     &    IWHANDLER, BEGS_BLR_DYNAMIC )
C
C      Purpose:
C      =======
C
C      Save BEGS_BLR_DYNAMIC. DMUMPS_BLR_INIT_FRONT has already been called.
C
C      Arguments:
C      =========
C
       INTEGER, DIMENSION(:), POINTER :: BEGS_BLR_DYNAMIC
       INTEGER, INTENT(IN) :: IWHANDLER
       INTEGER :: I
C
C      Local variables
C      ===============
C
       IF ( IWHANDLER .GT. size(BLR_ARRAY) .OR. IWHANDLER .LE. 0 ) THEN
C        Since DMUMPS_BLR_INIT_FRONT has been called, BLR_ARRAY
C     should be large enough
          WRITE(*,*) "Internal error 1 in DMUMPS_BLR_SAVE_BEGS_BLR_DYN"
          CALL MUMPS_ABORT()
       ENDIF
       IF ( BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT.LT.0) THEN
!  PANELS_NOTUSED or BLR_ARRAY_FREE
          WRITE(*,*) "Internal error 2 in DMUMPS_BLR_SAVE_BEGS_BLR_DYN"
          CALL MUMPS_ABORT()
       ENDIF
       DO I=1,size(BEGS_BLR_DYNAMIC)
          BLR_ARRAY(IWHANDLER)%BEGS_BLR_DYNAMIC(I) = BEGS_BLR_DYNAMIC(I)
       ENDDO
C
       RETURN
       END SUBROUTINE DMUMPS_BLR_SAVE_BEGS_BLR_DYN
C
       SUBROUTINE DMUMPS_BLR_RETRIEVE_BEGS_BLR_L 
     &            ( IWHANDLER, BEGS_BLR_L )
C
C      Purpose:
C      =======
C
C      Retrieve BEGS_BLR_L pointer array associated to IWHANDLER
C
C      Arguments:
C      =========
C
       INTEGER, INTENT(IN) :: IWHANDLER
#if defined(MUMPS_NOF2003)
       INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_L
#else
       INTEGER, POINTER, DIMENSION(:), INTENT(OUT) :: BEGS_BLR_L
#endif
C
       IF ( IWHANDLER .GT. size(BLR_ARRAY) .OR. IWHANDLER .LE. 0 ) THEN
C        Since DMUMPS_BLR_INIT_FRONT has been called, BLR_ARRAY
C        should be large enough
         WRITE(*,*) 
     &       "Internal error 1 in DMUMPS_BLR_RETRIEVE_BEGS_BLR_L"
         CALL MUMPS_ABORT()
       ENDIF
       BEGS_BLR_L => BLR_ARRAY(IWHANDLER)%BEGS_BLR_L
       RETURN
       END SUBROUTINE DMUMPS_BLR_RETRIEVE_BEGS_BLR_L
C
       SUBROUTINE DMUMPS_BLR_RETRIEVE_BEGSBLR_STA
     &            ( IWHANDLER, BEGS_BLR_STATIC )
C
C      Purpose:
C      =======
C
C      Retrieve BEGS_BLR_STATIC pointer array associated to IWHANDLER
C
C      Arguments:
C      =========
C
       INTEGER, INTENT(IN) :: IWHANDLER
#if defined(MUMPS_NOF2003)
       INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_STATIC
#else
       INTEGER, POINTER, DIMENSION(:), INTENT(OUT) :: BEGS_BLR_STATIC
#endif
C
       IF ( IWHANDLER .GT. size(BLR_ARRAY) .OR. IWHANDLER .LE. 0 ) THEN
C        Since DMUMPS_BLR_INIT_FRONT has been called, BLR_ARRAY
C        should be large enough
         WRITE(*,*) 
     &       "Internal error 1 in DMUMPS_BLR_RETRIEVE_BEGSBLR_STA"
         CALL MUMPS_ABORT()
       ENDIF
       BEGS_BLR_STATIC => BLR_ARRAY(IWHANDLER)%BEGS_BLR_STATIC
       RETURN
       END SUBROUTINE DMUMPS_BLR_RETRIEVE_BEGSBLR_STA
       SUBROUTINE DMUMPS_BLR_RETRIEVE_BEGSBLR_DYN
     &            ( IWHANDLER, BEGS_BLR_DYNAMIC )
C
C      Purpose:
C      =======
C
C      Retrieve BEGS_BLR_DYNAMIC pointer array associated to IWHANDLER
C
C      Arguments:
C      =========
C
       INTEGER, INTENT(IN) :: IWHANDLER
#if defined(MUMPS_NOF2003)
       INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_DYNAMIC
#else
       INTEGER, POINTER, DIMENSION(:), INTENT(OUT) :: BEGS_BLR_DYNAMIC
#endif
C
       IF ( IWHANDLER .GT. size(BLR_ARRAY) .OR. IWHANDLER .LE. 0 ) THEN
C        Since DMUMPS_BLR_INIT_FRONT has been called, BLR_ARRAY
C        should be large enough
         WRITE(*,*) 
     &       "Internal error 1 in DMUMPS_BLR_RETRIEVE_BEGSBLR_DYN"
         CALL MUMPS_ABORT()
       ENDIF
       BEGS_BLR_DYNAMIC => BLR_ARRAY(IWHANDLER)%BEGS_BLR_DYNAMIC
       RETURN
       END SUBROUTINE DMUMPS_BLR_RETRIEVE_BEGSBLR_DYN
C
       SUBROUTINE DMUMPS_BLR_RETRIEVE_BEGS_BLR_C
     &            ( IWHANDLER, BEGS_BLR_COL, NB_PANELS )
C
C      Purpose:
C      =======
C
C      Retrieve BEGS_BLR_COL pointer array associated to IWHANDLER
C
C      Arguments:
C      =========
C
       INTEGER, INTENT(IN)  :: IWHANDLER
       INTEGER, INTENT(OUT) :: NB_PANELS
#if defined(MUMPS_NOF2003)
       INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_COL
#else
       INTEGER, POINTER, DIMENSION(:), INTENT(OUT) :: BEGS_BLR_COL
#endif
C
       IF ( IWHANDLER .GT. size(BLR_ARRAY) .OR. IWHANDLER .LE. 0 ) THEN
C        Since DMUMPS_BLR_INIT_FRONT has been called, BLR_ARRAY
C        should be large enough
         WRITE(*,*) 
     &   "Internal error 1 in  DMUMPS_BLR_RETRIEVE_BEGS_BLR_C"
         CALL MUMPS_ABORT()
       ENDIF
       BEGS_BLR_COL => BLR_ARRAY(IWHANDLER)%BEGS_BLR_COL
       NB_PANELS     = BLR_ARRAY(IWHANDLER)%NB_PANELS
       RETURN
       END SUBROUTINE DMUMPS_BLR_RETRIEVE_BEGS_BLR_C
C
C
       SUBROUTINE DMUMPS_BLR_RETRIEVE_NB_PANELS
     &            ( IWHANDLER, NB_PANELS )
C
C      Purpose:
C      =======
C
C      Retrieve NB_PANELS associated to IWHANDLER
C
C      Arguments:
C      =========
C
       INTEGER, INTENT(IN)  :: IWHANDLER
       INTEGER, INTENT(OUT) :: NB_PANELS
       IF ( IWHANDLER .GT. size(BLR_ARRAY) .OR. IWHANDLER .LE. 0 ) THEN
C        Since DMUMPS_BLR_INIT_FRONT has been called, BLR_ARRAY
C        should be large enough
         WRITE(*,*) 
     &       "Internal error 1 in DMUMPS_BLR_RETRIEVE_NB_PANELS"
         CALL MUMPS_ABORT()
       ENDIF
       NB_PANELS = BLR_ARRAY(IWHANDLER)%NB_PANELS
       RETURN
       END SUBROUTINE DMUMPS_BLR_RETRIEVE_NB_PANELS
C
C
       SUBROUTINE DMUMPS_BLR_DEC_AND_RETRIEVE_L(IWHANDLER, IPANEL,
     &                                      BEGS_BLR_L, THELRBPANEL, 
     &                                      NBDEC )
C
C      Purpose:
C      =======
C
C      Retrieve BEGS_BLR_L and panel IPANEL in pointer THEPANEL 
C      for front identified by IWHANDLER.
C      Decrement NB_ACCESSES_LEFT of value NBDEC
C
C      Arguments:
C      =========
C
       INTEGER, INTENT(IN) :: IWHANDLER
       INTEGER, INTENT(IN) :: IPANEL, NBDEC
#if defined(MUMPS_NOF2003)
       INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_L
       TYPE(LRB_TYPE), POINTER, DIMENSION(:) :: THELRBPANEL
#else
       INTEGER, POINTER, DIMENSION(:), INTENT(OUT) :: BEGS_BLR_L
       TYPE(LRB_TYPE), INTENT(OUT), DIMENSION(:), POINTER :: THELRBPANEL
#endif
C
       IF ( IWHANDLER .GT. size(BLR_ARRAY) .OR. IWHANDLER .LE. 0 ) THEN
C        Since DMUMPS_BLR_INIT_FRONT has been called, BLR_ARRAY
C        should be large enough
         WRITE(*,*) "Internal error 1 in DMUMPS_BLR_DEC_AND_RETRIEVE_L",
     &    "IPANEL=", IPANEL
         CALL MUMPS_ABORT()
       ENDIF
       IF ( .NOT. associated(BLR_ARRAY(IWHANDLER)%PANELS_L)) THEN
         WRITE(*,*) "Internal error 2 in DMUMPS_BLR_DEC_AND_RETRIEVE_L",
     &    "IPANEL=", IPANEL
         CALL MUMPS_ABORT()
       ENDIF
       IF ( .NOT.
     & associated(BLR_ARRAY(IWHANDLER)%PANELS_L(IPANEL)%LRB_PANEL) )
     &      THEN
         WRITE(*,*) "Internal error 3 in DMUMPS_BLR_DEC_AND_RETRIEVE_L",
     &    "IPANEL=", IPANEL
         CALL MUMPS_ABORT()
       ENDIF
C
       CALL DMUMPS_BLR_RETRIEVE_BEGS_BLR_L( IWHANDLER, BEGS_BLR_L )
       THELRBPANEL =>
     &        BLR_ARRAY(IWHANDLER)%PANELS_L(IPANEL)%LRB_PANEL
C
       BLR_ARRAY(IWHANDLER)%PANELS_L(IPANEL)%NB_ACCESSES_LEFT = 
     & BLR_ARRAY(IWHANDLER)%PANELS_L(IPANEL)%NB_ACCESSES_LEFT - NBDEC
       RETURN
       END SUBROUTINE DMUMPS_BLR_DEC_AND_RETRIEVE_L
C
       LOGICAL FUNCTION DMUMPS_BLR_EMPTY_PANEL_LORU
     & (IWHANDLER, LorU, IPANEL)
C
C      Purpose:
C      =======
C      Check if panel IPANEL for front identified by
C      IWHANDLER is empty (in which case one should not
C      retrieve it for computations).
C
       INTEGER, INTENT(IN) :: LorU, IPANEL, IWHANDLER
       IF ( IWHANDLER .GT. size(BLR_ARRAY) .OR. IWHANDLER .LE. 0 ) THEN
C        Since DMUMPS_BLR_INIT_FRONT has been called, BLR_ARRAY
C        should be large enough
         WRITE(*,*) 
     &    "Internal error 1 in DMUMPS_BLR_EMPTY_PANEL_LORU, ",
     &    "IWHANDLER=", IWHANDLER
         CALL MUMPS_ABORT()
       ENDIF
C
       IF (LorU.EQ.0) THEN
C        L Panel
         IF ( .NOT. associated(BLR_ARRAY(IWHANDLER)%PANELS_L)) THEN
           WRITE(*,*) 
     &     "Internal error 2 in DMUMPS_BLR_EMPTY_PANEL_LORU, ",
     &     "IWHANDLER=", IWHANDLER
           CALL MUMPS_ABORT()
         ENDIF
         DMUMPS_BLR_EMPTY_PANEL_LORU = .NOT.
     &   associated(BLR_ARRAY(IWHANDLER)%PANELS_L(IPANEL)%LRB_PANEL)
       ELSE
         IF ( .NOT. associated(BLR_ARRAY(IWHANDLER)%PANELS_U)) THEN
           WRITE(*,*) 
     &     "Internal error 3 in DMUMPS_BLR_EMPTY_PANEL_LORU, ",
     &     "IWHANDLER=", IWHANDLER
           CALL MUMPS_ABORT()
         ENDIF
C        U panel
         DMUMPS_BLR_EMPTY_PANEL_LORU = .NOT.
     &   associated(BLR_ARRAY(IWHANDLER)%PANELS_U(IPANEL)%LRB_PANEL)
       ENDIF
       RETURN
       END FUNCTION DMUMPS_BLR_EMPTY_PANEL_LORU
C
C
       SUBROUTINE DMUMPS_BLR_RETRIEVE_PANEL_LORU
     &                          (IWHANDLER, LORU, IPANEL,
     &                           THELRBPANEL)
C
C      Purpose:
C      =======
C
C      Retrieve panel IPANEL in pointer THELRBBPANEL for front identified by IWHANDLER.
C
C      Arguments:
C      =========
C
       INTEGER, INTENT(IN) :: IWHANDLER
       INTEGER, INTENT(IN) :: LORU 
C      LORU :
C             0 = L Panel
C             1 = U Panel
C             others  : U panel
       INTEGER, INTENT(IN) :: IPANEL
#if defined(MUMPS_NOF2003)
       TYPE(LRB_TYPE), POINTER, DIMENSION(:) :: THELRBPANEL
#else
       TYPE(LRB_TYPE), INTENT(OUT), DIMENSION(:), POINTER :: THELRBPANEL
#endif
C
       IF ( IWHANDLER .GT. size(BLR_ARRAY) .OR. IWHANDLER .LE. 0 ) THEN
C        Since DMUMPS_BLR_INIT_FRONT has been called, BLR_ARRAY
C        should be large enough
         WRITE(*,*) 
     &    "Internal error 1 in DMUMPS_BLR_RETRIEVE_PANEL_LORU",
     &    "IWHANDLER=", IWHANDLER
         CALL MUMPS_ABORT()
       ENDIF
C
       IF (LORU.EQ.0) THEN
C       L Panel
          IF ( .NOT. associated(BLR_ARRAY(IWHANDLER)%PANELS_L)) THEN
         WRITE(*,*) 
     &    "Internal error 2 in DMUMPS_BLR_RETRIEVE_PANEL_LORU",
     &    " IWHANDLER=", IWHANDLER
         CALL MUMPS_ABORT()
        ENDIF
        IF ( .NOT.
     &  associated(BLR_ARRAY(IWHANDLER)%PANELS_L(IPANEL)%LRB_PANEL) )
     &       THEN
         WRITE(*,*) 
     &    "Internal error 3 in DMUMPS_BLR_RETRIEVE_PANEL_LORU",
     &    " IPANEL=", IPANEL
         CALL MUMPS_ABORT()
        ENDIF
        THELRBPANEL =>
     &        BLR_ARRAY(IWHANDLER)%PANELS_L(IPANEL)%LRB_PANEL
       ELSE
C       U Panel
          IF ( .NOT. associated(BLR_ARRAY(IWHANDLER)%PANELS_U)) THEN
         WRITE(*,*) 
     &    "Internal error 4 in DMUMPS_BLR_RETRIEVE_PANEL_LORU",
     &    " IWHANDLER=", IWHANDLER
         CALL MUMPS_ABORT()
        ENDIF
        IF ( .NOT.
     &  associated(BLR_ARRAY(IWHANDLER)%PANELS_U(IPANEL)%LRB_PANEL) )
     &       THEN
         WRITE(*,*) 
     &    "Internal error 5 in DMUMPS_BLR_RETRIEVE_PANEL_LORU",
     &    " IPANEL=", IPANEL
         CALL MUMPS_ABORT()
        ENDIF
        THELRBPANEL =>
     &        BLR_ARRAY(IWHANDLER)%PANELS_U(IPANEL)%LRB_PANEL
C       NB_ACCESSES_LEFT Not used in current version
       ENDIF
       RETURN
       END SUBROUTINE DMUMPS_BLR_RETRIEVE_PANEL_LORU
       SUBROUTINE DMUMPS_BLR_RETRIEVE_DIAG_BLOCK
     &                          (IWHANDLER, IPANEL,
     &                           THEBLOCK)
C
C      Purpose:
C      =======
C
C      Retrieve panel IPANEL in pointer THEPANEL for front 
C      identified by IWHANDLER.
C
C      FIXME: if we enable to save a null pointer in 
C             DMUMPS_BLR_SAVE_PANEL_LORU
C             one could enable to retrieve a null pointer in 
C             DMUMPS_BLR_RETRIEVE_DIAG_BLOCK
C             This would enable not to allocate panels
C             empty because of pivoting and to easily test
C             for empty panels (as already done
C             in fac_process_end_facto_slave.F and during solve)
C
C      Arguments:
C      =========
C
       INTEGER, INTENT(IN) :: IWHANDLER
       INTEGER, INTENT(IN) :: IPANEL
#if defined(MUMPS_NOF2003)
       DOUBLE PRECISION, POINTER :: THEBLOCK(:)
#else
       DOUBLE PRECISION, POINTER, INTENT(OUT) :: THEBLOCK(:)
#endif
C
       IF ( IWHANDLER .GT. size(BLR_ARRAY) .OR. IWHANDLER .LE. 0 ) THEN
C        Since DMUMPS_BLR_INIT_FRONT has been called, BLR_ARRAY
C        should be large enough
         WRITE(*,*) 
     &    "Internal error 1 in DMUMPS_BLR_RETRIEVE_DIAG_BLOCK",
     &    "IPANEL=", IPANEL
         CALL MUMPS_ABORT()
       ENDIF
C
       IF ( .NOT. associated(BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS)) THEN
         WRITE(*,*) 
     &    "Internal error 2 in DMUMPS_BLR_RETRIEVE_DIAG_BLOCK",
     &    "IPANEL=", IPANEL
         CALL MUMPS_ABORT()
       ENDIF
       IF ( .NOT.
     & associated(BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS(IPANEL)%DIAG_BLOCK))
     &      THEN
         WRITE(*,*) 
     &    "Internal error 3 in DMUMPS_BLR_RETRIEVE_DIAG_BLOCK",
     &    "IPANEL=", IPANEL
         CALL MUMPS_ABORT()
       ENDIF
       THEBLOCK =>
     &        BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS(IPANEL)%DIAG_BLOCK
       RETURN
       END SUBROUTINE DMUMPS_BLR_RETRIEVE_DIAG_BLOCK
       SUBROUTINE DMUMPS_BLR_RETRIEVE_CB_LRB
     &                          (IWHANDLER, THECB)
C
C      Purpose:
C      =======
C
C      Retrieve CB_LRB_PTR pointer for front identified by IWHANDLER.
C
C      Arguments:
C      =========
C
       INTEGER, INTENT(IN) :: IWHANDLER
#if defined(MUMPS_NOF2003)
       TYPE(LRB_TYPE), POINTER :: THECB(:,:)
#else
       TYPE(LRB_TYPE), POINTER, INTENT(OUT) :: THECB(:,:)
#endif
C
       IF ( IWHANDLER .GT. size(BLR_ARRAY) .OR. IWHANDLER .LE. 0 ) THEN
C        Since DMUMPS_BLR_INIT_FRONT has been called, BLR_ARRAY
C        should be large enough
         WRITE(*,*) "Internal error 1 in DMUMPS_BLR_RETRIEVE_CB_LRB"
         CALL MUMPS_ABORT()
       ENDIF
C
       IF ( .NOT. associated(BLR_ARRAY(IWHANDLER)%CB_LRB)) THEN
         WRITE(*,*) "Internal error 2 in DMUMPS_BLR_RETRIEVE_CB_LRB"
         CALL MUMPS_ABORT()
       ENDIF
       THECB => BLR_ARRAY(IWHANDLER)%CB_LRB
       RETURN
       END SUBROUTINE DMUMPS_BLR_RETRIEVE_CB_LRB
       SUBROUTINE DMUMPS_BLR_SAVE_NFS4FATHER
     &            ( IWHANDLER, NFS4FATHER  )
C
C      Purpose:
C      =======
C
C      Save NFS4FATHER associated to IWHANDLER
C
C      Arguments:
C      =========
C
       INTEGER, INTENT(IN) :: IWHANDLER
       INTEGER, INTENT(IN) :: NFS4FATHER 
       IF ( IWHANDLER .GT. size(BLR_ARRAY) .OR. IWHANDLER .LE. 0 ) THEN
C        Since DMUMPS_BLR_INIT_FRONT has been called, BLR_ARRAY
C        should be large enough
         WRITE(*,*) 
     &       "Internal error 1 in DMUMPS_BLR_RETRIEVE_NFS4FATHER"
         CALL MUMPS_ABORT()
       ENDIF
       BLR_ARRAY(IWHANDLER)%NFS4FATHER = NFS4FATHER
       RETURN
       END SUBROUTINE DMUMPS_BLR_SAVE_NFS4FATHER
       SUBROUTINE DMUMPS_BLR_RETRIEVE_NFS4FATHER
     &            ( IWHANDLER, NFS4FATHER  )
C
C      Purpose:
C      =======
C
C      Retrieve NFS4FATHER associated to IWHANDLER
C
C      Arguments:
C      =========
C
       INTEGER, INTENT(IN)  :: IWHANDLER
       INTEGER, INTENT(OUT) :: NFS4FATHER 
       IF ( IWHANDLER .GT. size(BLR_ARRAY) .OR. IWHANDLER .LE. 0 ) THEN
C        Since DMUMPS_BLR_INIT_FRONT has been called, BLR_ARRAY
C        should be large enough
         WRITE(*,*) 
     &       "Internal error 1 in DMUMPS_BLR_RETRIEVE_NFS4FATHER"
         CALL MUMPS_ABORT()
       ENDIF
       NFS4FATHER = BLR_ARRAY(IWHANDLER)%NFS4FATHER
       RETURN
       END SUBROUTINE DMUMPS_BLR_RETRIEVE_NFS4FATHER
C
       SUBROUTINE DMUMPS_BLR_SAVE_M_ARRAY (
     &    IWHANDLER, M_ARRAY, INFO)
C
C      Purpose:
C      =======
C      Save existing array M_ARRAY (a copy is performed so that
C      the argument can be an array rather than a pointer).
C      (size of which includes non eliminated variables of the front)
C       DMUMPS_BLR_INIT_FRONT has already been called.
C
C      Arguments:
C      =========
C
       DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: M_ARRAY
       INTEGER, INTENT(IN) ::  IWHANDLER
       INTEGER, INTENT(INOUT) :: INFO(2)
C
C      Local variables:
C      ===============
C
       INTEGER :: I, IERR
       IF ( IWHANDLER .GT. size(BLR_ARRAY) .OR. IWHANDLER .LE. 0 ) THEN
C        Since DMUMPS_BLR_INIT_FRONT has been called, BLR_ARRAY
C        should be large enough
          WRITE(*,*) "Internal error 1 in DMUMPS_BLR_SAVE_M_ARRAY"
          CALL MUMPS_ABORT()
       ENDIF
C       IF ( BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT.LT.0) THEN
C         !  PANELS_NOTUSED or BLR_ARRAY_FREE
C          WRITE(*,*) "Internal error 2 in DMUMPS_BLR_SAVE_BEGS_BLR_C", 
C     &        BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT
C          CALL MUMPS_ABORT()
C       ENDIF
       ALLOCATE(BLR_ARRAY(IWHANDLER)%M_ARRAY(size(M_ARRAY)),
     &          stat=IERR)
       IF (IERR > 0 ) THEN
         INFO(1)=-13
         INFO(2)=size(M_ARRAY)
         RETURN
       ENDIF
       DO I=1,size(M_ARRAY)
          BLR_ARRAY(IWHANDLER)%M_ARRAY(I) = M_ARRAY(I)
       ENDDO
C
       BLR_ARRAY(IWHANDLER)%NFS4FATHER = size(M_ARRAY)
C
       RETURN
       END SUBROUTINE DMUMPS_BLR_SAVE_M_ARRAY
       SUBROUTINE DMUMPS_BLR_RETRIEVE_M_ARRAY ( IWHANDLER, M_ARRAY)
C
       IMPLICIT NONE
C
C      Purpose:
C      =======
C
C      Retrieve M_ARRAY 
C
C      Arguments:
C      =========
C
       INTEGER, INTENT(IN) ::  IWHANDLER
#if defined(MUMPS_NOF2003)
C      Intent not allowed
       DOUBLE PRECISION, DIMENSION(:), POINTER :: M_ARRAY
#else
       DOUBLE PRECISION, DIMENSION(:), POINTER, INTENT(OUT) :: M_ARRAY
#endif
       IF ( IWHANDLER .GT. size(BLR_ARRAY) .OR. IWHANDLER .LE. 0 ) THEN
C        Since DMUMPS_BLR_INIT_FRONT has been called, BLR_ARRAY
C        should be large enough
          WRITE(*,*) "Internal error 1 in DMUMPS_BLR_RETRIEVE_M_ARRAY"
          CALL MUMPS_ABORT()
       ENDIF
       M_ARRAY =>  BLR_ARRAY(IWHANDLER)%M_ARRAY
C
       RETURN
       END SUBROUTINE DMUMPS_BLR_RETRIEVE_M_ARRAY
C
       SUBROUTINE DMUMPS_BLR_FREE_M_ARRAY ( IWHANDLER )
C
       IMPLICIT NONE
C
C      Purpose:
C      =======
C
C      Free M_ARRAY and reset NFS4FATHER
C
C      Arguments:
C      =========
C
       INTEGER, INTENT(IN) ::  IWHANDLER
       IF ( IWHANDLER .GT. size(BLR_ARRAY) .OR. IWHANDLER .LE. 0 ) THEN
C        Since DMUMPS_BLR_INIT_FRONT has been called, BLR_ARRAY
C        should be large enough
          WRITE(*,*) "Internal error 1 in DMUMPS_BLR_FREE_M_ARRAY"
          CALL MUMPS_ABORT()
       ENDIF
       IF (associated(BLR_ARRAY(IWHANDLER)%M_ARRAY)) THEN
          DEALLOCATE(BLR_ARRAY(IWHANDLER)%M_ARRAY)
          NULLIFY(BLR_ARRAY(IWHANDLER)%M_ARRAY)
       ENDIF
       BLR_ARRAY(IWHANDLER)%NFS4FATHER       = NFS4FATHER_NOTINIT
C
       RETURN
       END SUBROUTINE DMUMPS_BLR_FREE_M_ARRAY
C
       SUBROUTINE DMUMPS_BLR_DEC_AND_TRYFREE_L( IWHANDLER, IPANEL,
     &                            KEEP8, K34, NBDEC)      
       IMPLICIT NONE
C
C      Purpose:
C      =======
C
C      Decrement the number of times it is used and 
C      try to free a panel after it has
C      been used for the last time.
C      Decrement NB_ACCESSES_LEFT of value NBDEC
C
       INTEGER, INTENT(IN) :: IWHANDLER, IPANEL, K34, NBDEC
       INTEGER(8) :: KEEP8(150)
C
C      Local variable
C      ==============
C
C      enable calling DMUMPS_BLR_TRY_FREE_PANEL when module was not called
       IF (IWHANDLER.LE.0) RETURN
C      
       IF ( BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT.LT.0)
         !  PANELS_NOTUSED or BLR_ARRAY_FREE
     &    RETURN
C      Decrement number of accesses
       BLR_ARRAY(IWHANDLER)%PANELS_L(IPANEL)%NB_ACCESSES_LEFT = 
     & BLR_ARRAY(IWHANDLER)%PANELS_L(IPANEL)%NB_ACCESSES_LEFT - NBDEC
C      Try to free panel
       CALL DMUMPS_BLR_TRY_FREE_PANEL (IWHANDLER, IPANEL,
     &                   KEEP8, K34)
       RETURN
       END SUBROUTINE DMUMPS_BLR_DEC_AND_TRYFREE_L
       SUBROUTINE DMUMPS_BLR_FORCE_FREE_PANEL_L( IWHANDLER, IPANEL, 
     &           KEEP8, K34 )
       IMPLICIT NONE
C
C      Purpose:
C      =======
C
C      Force to free a panel independently 
C      of the nb of times it has been used 
C
       INTEGER, INTENT(IN) :: IWHANDLER, IPANEL
       INTEGER(8) :: KEEP8(150)
       INTEGER, INTENT(IN)    :: K34
C
C      Local variable
C      ==============
C
       TYPE(blr_panel_type), POINTER  :: THEPANEL
C
C      enable calling DMUMPS_BLR_TRY_FREE_PANEL when module was not called
       IF (IWHANDLER.LE.0) RETURN
C      
       IF ( BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT.LT.0)
         !  PANELS_L_NOTUSED or BLR_ARRAY_FREE
     &    RETURN
       THEPANEL => BLR_ARRAY(IWHANDLER)%PANELS_L(IPANEL)
!       IF ( THEPANEL%NB_ACCESSES_LEFT .EQ. 0 ) THEN
         IF (associated(THEPANEL%LRB_PANEL)) THEN
          IF (size(THEPANEL%LRB_PANEL) .GT.0) THEN
           CALL DEALLOC_BLR_PANEL(THEPANEL%LRB_PANEL,
     &        size(THEPANEL%LRB_PANEL), KEEP8, K34)
          ENDIF
          DEALLOCATE(THEPANEL%LRB_PANEL)
          NULLIFY(THEPANEL%LRB_PANEL)
         ENDIF
         THEPANEL%NB_ACCESSES_LEFT = PANELS_FREED 
!       ENDIF
       !  DEALLOCATE(THEPANEL)
       RETURN
C
       END SUBROUTINE DMUMPS_BLR_FORCE_FREE_PANEL_L
       SUBROUTINE DMUMPS_BLR_TRY_FREE_PANEL( IWHANDLER, IPANEL, 
     &           KEEP8, K34 )
       IMPLICIT NONE
C
C      Purpose:
C      =======
C
C      Try to free a panel after it has
C      been used for the last time.
C
       INTEGER, INTENT(IN) :: IWHANDLER, IPANEL
       INTEGER(8) :: KEEP8(150)
       INTEGER, INTENT(IN)    :: K34
C
C      Local variable
C      ==============
C
       TYPE(blr_panel_type), POINTER  :: THEPANEL
C
C      enable calling DMUMPS_BLR_TRY_FREE_PANEL when module was not called
       IF (IWHANDLER.LE.0) RETURN
C      
       IF ( BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT.LT.0.OR.
     &      BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT.eq.huge(IPANEL) )
         !  PANELS_L_NOTUSED or BLR_ARRAY_FREE
     &    RETURN
       THEPANEL => BLR_ARRAY(IWHANDLER)%PANELS_L(IPANEL)
       IF ( THEPANEL%NB_ACCESSES_LEFT .EQ. 0 ) THEN
         IF (associated(THEPANEL%LRB_PANEL)) THEN
          IF (size(THEPANEL%LRB_PANEL) .GT.0) THEN
           CALL DEALLOC_BLR_PANEL(THEPANEL%LRB_PANEL,
     &        size(THEPANEL%LRB_PANEL), KEEP8, K34)
          ENDIF
          DEALLOCATE(THEPANEL%LRB_PANEL)
          NULLIFY(THEPANEL%LRB_PANEL)
         ENDIF
         THEPANEL%NB_ACCESSES_LEFT = PANELS_FREED 
       ENDIF
       !  DEALLOCATE(THEPANEL)
       RETURN
C
       END SUBROUTINE DMUMPS_BLR_TRY_FREE_PANEL
C
       SUBROUTINE DMUMPS_BLR_FREE_CB_LRB ( IWHANDLER, FREE_ONLY_STRUCT,
     &                                     KEEP8, K34 )
       IMPLICIT NONE
C
C      Purpose:
C      =======
C
C      Free CB_LRB
C
       INTEGER, INTENT(IN) :: IWHANDLER, K34
       LOGICAL, INTENT(IN) :: FREE_ONLY_STRUCT
       INTEGER(8) :: KEEP8(150)
C
C      Local variable
C      ==============
C
       TYPE(LRB_TYPE), POINTER  :: CB_LRB(:,:)
       INTEGER :: IPANEL, JPANEL
       TYPE(LRB_TYPE), POINTER  :: THELRB 
C
C            If ASM_ON_THE_FLY is TRUE, only CB_LRB structure is
C            left to deallocate
       IF (BLR_ARRAY(IWHANDLER)%IsT2.AND.
     &      .NOT.BLR_ARRAY(IWHANDLER)%IsSLAVE) THEN
         write(*,*) 'Internal error 1 in DMUMPS_BLR_FREE_CB_LRB'
         CALL MUMPS_ABORT()
       ENDIF
       CB_LRB => BLR_ARRAY(IWHANDLER)%CB_LRB
       IF (.NOT.associated(CB_LRB)) THEN
         write(*,*) 'Internal error 2 in DMUMPS_BLR_FREE_CB_LRB'
         CALL MUMPS_ABORT()
       ENDIF
       IF (.NOT.FREE_ONLY_STRUCT) THEN
         DO IPANEL = 1,size(CB_LRB,1)
         DO JPANEL = 1,size(CB_LRB,2)
           THELRB => CB_LRB(IPANEL,JPANEL)
           IF (associated(THELRB)) THEN
             CALL DEALLOC_LRB(THELRB, KEEP8, K34)
           ENDIF
         ENDDO
         ENDDO
       ENDIF
       DEALLOCATE(BLR_ARRAY(IWHANDLER)%CB_LRB)
       NULLIFY(BLR_ARRAY(IWHANDLER)%CB_LRB)
       RETURN
C
       END SUBROUTINE DMUMPS_BLR_FREE_CB_LRB
C
C
       SUBROUTINE DMUMPS_BLR_FREE_ALL_PANELS ( IWHANDLER, 
     &         LorU, KEEP8, K34)
       IMPLICIT NONE
C
C      Purpose:
C      =======
C
C      Free ALL panels
C
       INTEGER, INTENT(IN) :: IWHANDLER, LorU, K34
       INTEGER(8) :: KEEP8(150)
C
C      Local variable
C      ==============
C
       INTEGER             :: IPANEL
       INTEGER             :: IDUMMY, JDUMMY
       TYPE(blr_panel_type), POINTER  :: THEPANEL
       TYPE(diag_block_type), POINTER  :: THEBLOCK
       INTEGER(8) :: MEM_FREED
C
C      enable calling DMUMPS_BLR_TRY_FREE_PANEL when module was not called
       IF (IWHANDLER.LE.0) RETURN
C      
       IF (BLR_ARRAY(IWHANDLER)%NB_ACCESSES_INIT.EQ.
     &      PANELS_NOTUSED) RETURN
       IF (LorU.EQ.0.OR.LorU.EQ.2) THEN
C{        Free L panels         
C        check if pointer is associated
C        in case of error during allocation then 
C FIXME   PA2ALL not clear to me why BLR_ARRAY(IWHANDLER)%PANELS_L
C         is not associated and of size 1 ! 
C         causing segfault at instruction  THEPANEL => BLR_ARRAY...
         IF (associated(BLR_ARRAY(IWHANDLER)%PANELS_L)) THEN
         DO IPANEL = 1, size(BLR_ARRAY(IWHANDLER)%PANELS_L)
         THEPANEL => BLR_ARRAY(IWHANDLER)%PANELS_L(IPANEL)
C         IF ( THEPANEL%NB_ACCESSES_LEFT .EQ. 0 ) THEN
           IF (associated(THEPANEL%LRB_PANEL)) THEN
            IF (size(THEPANEL%LRB_PANEL) .GT.0) THEN
             CALL DEALLOC_BLR_PANEL(THEPANEL%LRB_PANEL,
     &          size(THEPANEL%LRB_PANEL), KEEP8, K34)
            ENDIF
            DEALLOCATE(THEPANEL%LRB_PANEL)
            NULLIFY(THEPANEL%LRB_PANEL)
           ENDIF
           THEPANEL%NB_ACCESSES_LEFT = PANELS_FREED 
C         ENDIF
         ENDDO
         ENDIF
C}
       ENDIF
       IF (LorU.GE.1.AND..NOT.BLR_ARRAY(IWHANDLER)%IsSYM) THEN
C        Free U panels         
C        check if pointer is associated
C        in case of error during allocation then 
         IF (associated(BLR_ARRAY(IWHANDLER)%PANELS_U)) THEN
         DO IPANEL = 1, size(BLR_ARRAY(IWHANDLER)%PANELS_U)
           THEPANEL => BLR_ARRAY(IWHANDLER)%PANELS_U(IPANEL)
          IF (associated(THEPANEL%LRB_PANEL)) THEN
           IF (size(THEPANEL%LRB_PANEL) .GT.0) THEN
            CALL DEALLOC_BLR_PANEL(THEPANEL%LRB_PANEL,
     &        size(THEPANEL%LRB_PANEL), KEEP8, K34)
           ENDIF
           DEALLOCATE(THEPANEL%LRB_PANEL)
           NULLIFY(THEPANEL%LRB_PANEL)
          ENDIF
          THEPANEL%NB_ACCESSES_LEFT = PANELS_FREED 
         ENDDO
         ENDIF
       ENDIF
       IF (LorU.GE.1.AND..NOT.BLR_ARRAY(IWHANDLER)%IsSLAVE) THEN
C        DIAG_BLOCKS should also be processed only if U panel freed
C        check if pointer is associated
C        in case of error during allocation then 
         IF (associated(BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS)) THEN
           MEM_FREED = 0_8
           DO IPANEL = 1, size(BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS)
             THEBLOCK => BLR_ARRAY(IWHANDLER)%DIAG_BLOCKS(IPANEL)
             IF (associated(THEBLOCK%DIAG_BLOCK)) THEN
               MEM_FREED = MEM_FREED + int(size(THEBLOCK%DIAG_BLOCK),8)
               DEALLOCATE(THEBLOCK%DIAG_BLOCK)
               NULLIFY   (THEBLOCK%DIAG_BLOCK)
             ENDIF
#if !defined(NOMIXHPC)
             IF (associated(THEBLOCK%DIAG_STOC%VAL)) THEN
               MEM_FREED = MEM_FREED 
     &            + int(1 + (size(THEBLOCK%DIAG_STOC%VAL)*
     &                       K34 - 1 )/COSTS_PREC(PREC_DIAG_FR),8)
               DEALLOCATE(THEBLOCK%DIAG_STOC%VAL)
               NULLIFY   (THEBLOCK%DIAG_STOC%VAL)
             ENDIF
#endif
           ENDDO
           IF (MEM_FREED .GT. 0 ) THEN
             CALL MUMPS_DM_FAC_UPD_DYN_MEMCNTS(-MEM_FREED,
     &        .TRUE., KEEP8, 
     &        IDUMMY, JDUMMY, ! IFLAG, IERROR not used for deallocate
     &        .TRUE., .TRUE.)
           ENDIF
         ENDIF
       ENDIF
       RETURN
C
       END SUBROUTINE DMUMPS_BLR_FREE_ALL_PANELS
C
#if ! defined(NO_SAVE_RESTORE)
!     --------------------------------- DMUMPS_SAVE_RESTORE_BLR ----------------------------- !
      SUBROUTINE DMUMPS_SAVE_RESTORE_BLR(id_BLRARRAY_ENCODING
     &                 ,unit,MYID,mode
     &                 ,SIZE_GEST,SIZE_VARIABLES
     &                 ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &                 ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                 ,size_read,size_allocated,size_written
     &                 ,INFO)
      INCLUDE 'mpif.h'
C  =======
C  Purpose
C  =======
C
C     This routine is designed to manage a BLR_STRUC_T structure (save, restore, compute memory)
C     
C  ==========
C  Parameters
C  ==========
C
C     BLR_STRUC           : TYPE (BLR_STRUC_T) : the main structure
C
C     unit                : The unit of the file to be written or read
C
C     mode                : the type of operation to be performed by the routine 
C                           memory_save  = compute the size of the save file and of the structure
C                           save         = save the instace
C                           restore      = restore the instace
C
C     TOTAL_FILE_SIZE     : size of the file to be written or read
C
C     TOTAL_STRUC_SIZE    : size of the structure to be saved or restored
C
C     SIZE_INT            : size of an integer      
C
C     SIZE_ARITH_DEP      : size of an arithmetic dependent variable
C
C     SIZE_LOGICAL        : size of a logical     
C
C     INFO                : copies of of INFO(1) and INFO(2) to allow
C                           save/restore of failed instances
C      
      CHARACTER, DIMENSION(:), POINTER :: id_BLRARRAY_ENCODING
      INTEGER,intent(IN)::unit,MYID
      INTEGER,intent(IN) :: mode
      INTEGER,INTENT(OUT) :: SIZE_GEST
      INTEGER(8),intent(OUT) :: SIZE_VARIABLES
      INTEGER(8),intent(IN) :: TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
      INTEGER,intent(INOUT):: INFO(2)
      INTEGER,intent(IN):: SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
      INTEGER(8),intent(INOUT):: size_read,size_allocated,size_written
      INTEGER:: j1,NbRecords,NbSubRecords,size_array1,dummy,allocok,err
      INTEGER:: SIZE_GEST_BLR_ARRAY,SIZE_GEST_BLR_ARRAY_j1
      INTEGER(8):: SIZE_VARIABLES_BLR_ARRAY,SIZE_VARIABLES_BLR_ARRAY_j1
      INTEGER(4) :: I4
      NbRecords=0
      SIZE_GEST_BLR_ARRAY=0
      SIZE_GEST_BLR_ARRAY_j1=0
      SIZE_VARIABLES_BLR_ARRAY=0_8
      SIZE_VARIABLES_BLR_ARRAY_j1=0_8
      SIZE_GEST=0
      SIZE_VARIABLES=0_8
      if(mode.EQ.memory_save_mode.OR.mode.EQ.save_mode) then
         call DMUMPS_BLR_STRUC_TO_MOD(id_BLRARRAY_ENCODING)
      endif
      if(mode.EQ.memory_save_mode) then
         IF(associated(BLR_ARRAY)) THEN
C     NbRecords_BLR_ARRAY = 1 for the size of the array
            NbRecords=1
            SIZE_GEST=SIZE_INT
            SIZE_VARIABLES=0
            DO j1=1,size(BLR_ARRAY,1)
C     Call save_memory on BLR_ARRAY(j1)
               CALL DMUMPS_SAVE_RESTORE_BLR_STRUC(
     &              BLR_ARRAY(j1)
     &              ,unit,MYID,memory_save_mode
     &              ,SIZE_GEST_BLR_ARRAY_j1
     &              ,SIZE_VARIABLES_BLR_ARRAY_j1
     &              ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &              ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &              ,size_read,size_allocated,size_written
     &              ,INFO)
               SIZE_GEST_BLR_ARRAY=SIZE_GEST_BLR_ARRAY+
     &              SIZE_GEST_BLR_ARRAY_j1
               SIZE_VARIABLES_BLR_ARRAY=SIZE_VARIABLES_BLR_ARRAY+
     &              SIZE_VARIABLES_BLR_ARRAY_j1
            ENDDO
         ELSE
            NbRecords=2
            SIZE_GEST=SIZE_INT*2
            SIZE_VARIABLES=0
         ENDIF
      elseif(mode.EQ.save_mode) then
         IF(associated(BLR_ARRAY)) THEN
            NbRecords=1
            SIZE_GEST=SIZE_INT
            SIZE_VARIABLES=0
            write(unit,iostat=err) size(BLR_ARRAY,1)
            if(err.ne.0) then
               INFO(1) = -72
               CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &              INFO(2))
            endif
            IF ( INFO(1) .LT. 0 ) GOTO 100
            DO j1=1,size(BLR_ARRAY,1)
C     Call save on BLR_ARRAY(j1)
               CALL DMUMPS_SAVE_RESTORE_BLR_STRUC(
     &              BLR_ARRAY(j1)
     &              ,unit,MYID,save_mode
     &              ,SIZE_GEST_BLR_ARRAY_j1
     &              ,SIZE_VARIABLES_BLR_ARRAY_j1
     &              ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &              ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &              ,size_read,size_allocated,size_written
     &              ,INFO)
               IF ( INFO(1) .LT. 0 ) GOTO 100
            ENDDO
         ELSE
            NbRecords=2
            SIZE_GEST=SIZE_INT*2
            SIZE_VARIABLES=0
            write(unit,iostat=err) -999
            if(err.ne.0) then
               INFO(1) = -72
               CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &              INFO(2))
            endif
            IF ( INFO(1) .LT. 0 ) GOTO 100
            write(unit,iostat=err) -999
            if(err.ne.0) then
               INFO(1) = -72
               CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &              INFO(2))
            endif
            IF ( INFO(1) .LT. 0 ) GOTO 100
         ENDIF 
      elseif(mode.EQ.restore_mode) then
         nullify(BLR_ARRAY)
         read(unit,iostat=err) size_array1
         if(err.ne.0) THEN
            INFO(1) = -75
            CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &           ,INFO(2))
         endif
         IF ( INFO(1) .LT. 0 ) GOTO 100
         if(size_array1.EQ.-999) then
            NbRecords=2
            SIZE_GEST=SIZE_INT*2
            SIZE_VARIABLES=0
            read(unit,iostat=err) dummy
            if(err.ne.0) THEN
               INFO(1) = -75
               CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &              ,INFO(2))
            endif
            IF ( INFO(1) .LT. 0 ) GOTO 100
         else
C    NbRecords = 1 for the size of the array
            NbRecords=1
            SIZE_GEST=SIZE_INT
            SIZE_VARIABLES=0
            allocate(BLR_ARRAY(size_array1), stat=allocok)
            if (allocok .GT. 0) THEN
               INFO(1) = -78
               CALL MUMPS_SETI8TOI4(
     &              TOTAL_STRUC_SIZE-size_allocated
     &              ,INFO(2))
            endif
            DO j1=1,size_array1
C     call restore on BLR_ARRAY(j1)
               CALL DMUMPS_SAVE_RESTORE_BLR_STRUC(
     &              BLR_ARRAY(j1)
     &              ,unit,MYID,restore_mode
     &              ,SIZE_GEST_BLR_ARRAY_j1
     &              ,SIZE_VARIABLES_BLR_ARRAY_j1
     &              ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &              ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &              ,size_read,size_allocated,size_written
     &              ,INFO)
               SIZE_GEST_BLR_ARRAY=SIZE_GEST_BLR_ARRAY+
     &              SIZE_GEST_BLR_ARRAY_j1
               SIZE_VARIABLES_BLR_ARRAY=SIZE_VARIABLES_BLR_ARRAY+
     &              SIZE_VARIABLES_BLR_ARRAY_j1
            ENDDO
         endif
      endif
      if(mode.EQ.memory_save_mode) then
C     If the size to write (SIZE_VARIABLES) is greater than 2^31
C     Subrecords are created which need to be taken into account in
C     the file size computation
         NbSubRecords=int(SIZE_VARIABLES/huge(I4))
         IF(NbSubRecords.GT.0) then
            NbRecords=NbRecords+NbSubRecords
         ENDIF
      elseif(mode.EQ.save_mode) then
         size_written=size_written+SIZE_VARIABLES
     &        +int(SIZE_GEST,kind=8)
#if defined(MUMPS_NOF2003)
         size_written=size_written
     &        +int(2*SIZE_INT*NbRecords,kind=8)
#endif
      elseif(mode.EQ.restore_mode) then
         size_allocated=size_allocated+SIZE_VARIABLES
         size_read=size_read+SIZE_VARIABLES
     &        +int(SIZE_GEST,kind=8)
#if defined(MUMPS_NOF2003)
         size_read=size_read
     &        +int(2*SIZE_INT*NbRecords,kind=8)
#endif
      endif
      if(mode.EQ.memory_save_mode) then
         SIZE_VARIABLES=SIZE_VARIABLES+SIZE_VARIABLES_BLR_ARRAY
         SIZE_GEST=SIZE_GEST+SIZE_GEST_BLR_ARRAY
#if defined(MUMPS_NOF2003)
C        If the file is not written with access="stream"  (stream is
C        only in case of Fortran 2003), the record's length is written
C        at the beginning and at the end of each record.
C        This is done using 2 INTEGERs so we use 2*SIZE_INT more space
C        for each record
         SIZE_GEST=SIZE_GEST+2*SIZE_INT*NbRecords
#endif
      endif
      call DMUMPS_BLR_MOD_TO_STRUC(id_BLRARRAY_ENCODING)
 100  continue
      RETURN
      END SUBROUTINE DMUMPS_SAVE_RESTORE_BLR
!     ---------- DMUMPS_SAVE_RESTORE_BLR_STRUC ----------------------- !
      SUBROUTINE DMUMPS_SAVE_RESTORE_BLR_STRUC(BLR_STRUC
     &     ,unit,MYID,mode
     &     ,Local_SIZE_GEST,Local_SIZE_VARIABLES
     &     ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &     ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &     ,size_read,size_allocated,size_written
     &     ,INFO)     
      INCLUDE 'mpif.h'
C  =======
C  Purpose
C  =======
C
C     This routine is designed to manage a BLR_STRUC_T structure
C     (save, restore, compute memory)
C     
C  ==========
C  Parameters
C  ==========
C
C     BLR_STRUC           : TYPE (BLR_STRUC_T) : the main structure
C
C     unit                : The unit of the file to be written or read
C
C     mode                : the type of operation to be performed by the routine 
C                           memory_save  = compute the size of the save file and of the structure
C                           save         = save the instace
C                           restore      = restore the instace
C
C     TOTAL_FILE_SIZE     : size of the file to be written or read
C
C     TOTAL_STRUC_SIZE    : size of the structure to be saved or restored
C
C     SIZE_INT            : size of an integer      
C
C     SIZE_ARITH_DEP      : size of an arithmetic dependent variable
C
C     SIZE_LOGICAL        : size of a logical     
C
C     INFO1/INFO2         : copies of of INFO(1) and INFO(2) to allow
C                           save/restore of failed instances
C
      TYPE(BLR_STRUC_T) :: BLR_STRUC
      INTEGER,intent(IN)::unit,MYID
      INTEGER,intent(IN) :: mode
      INTEGER,INTENT(OUT) :: Local_SIZE_GEST
      INTEGER(8),intent(OUT) :: Local_SIZE_VARIABLES
      INTEGER,intent(INOUT):: INFO(2)
      INTEGER,intent(IN):: SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
      INTEGER(8),intent(IN) :: TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
      INTEGER(8),intent(INOUT):: size_read,size_allocated,size_written
      INTEGER, PARAMETER :: NBVARIABLES_BLR_STRUC_T = 15
      INTEGER, PARAMETER :: B_IsSYM=1
      INTEGER, PARAMETER :: B_IsT2=2
      INTEGER, PARAMETER :: B_IsSLAVE=3
      INTEGER, PARAMETER :: B_PANELS_L=4
      INTEGER, PARAMETER :: B_PANELS_U=5
      INTEGER, PARAMETER :: B_CB_LRB=6
      INTEGER, PARAMETER :: B_DIAG_BLOCKS=7
      INTEGER, PARAMETER :: B_BEGS_BLR_STATIC=8
      INTEGER, PARAMETER :: B_BEGS_BLR_DYNAMIC=9
      INTEGER, PARAMETER :: B_BEGS_BLR_L=10
      INTEGER, PARAMETER :: B_BEGS_BLR_COL=11
      INTEGER, PARAMETER :: B_NB_ACCESSES_INIT=12
      INTEGER, PARAMETER :: B_NB_PANELS=13
      INTEGER, PARAMETER :: B_NFS4FATHER=14
      INTEGER, PARAMETER :: B_M_ARRAY=15
      INTEGER(8),dimension(NBVARIABLES_BLR_STRUC_T)::
     &     SIZE_VARIABLES_BLR_STRUC_T
      INTEGER,dimension(NBVARIABLES_BLR_STRUC_T)::SIZE_GEST_BLR_STRUC_T
      INTEGER,dimension(NBVARIABLES_BLR_STRUC_T)::NbRecords_BLR_STRUC_T
      INTEGER:: size_array1,size_array2,dummy,allocok
      INTEGER:: err,i1,j1,j2,NbSubRecords
#if defined(MUMPS_NOF2003)
      INTEGER:: Local_NbRecords
#endif
      INTEGER::SIZE_GEST_PANELS_L,SIZE_GEST_PANELS_L_j1
      INTEGER(8)::SIZE_VARIABLES_PANELS_L,SIZE_VARIABLES_PANELS_L_j1
      INTEGER::SIZE_GEST_PANELS_U,SIZE_GEST_PANELS_U_j1
      INTEGER(8)::SIZE_VARIABLES_PANELS_U,SIZE_VARIABLES_PANELS_U_j1
      INTEGER::SIZE_GEST_CB_LRB,SIZE_GEST_CB_LRB_j1j2
      INTEGER(8)::SIZE_VARIABLES_CB_LRB,SIZE_VARIABLES_CB_LRB_j1j2
      INTEGER::SIZE_GEST_DIAG_BLOCKS,SIZE_GEST_DIAG_BLOCKS_j1
      INTEGER(8)::SIZE_VARIABLES_DIAG_BLOCKS
      INTEGER(8)::SIZE_VARIABLES_DIAG_BLOCKS_j1
      INTEGER(4)::I4
      SIZE_VARIABLES_BLR_STRUC_T(:)=0_8
      SIZE_GEST_BLR_STRUC_T(:)=0
      NbRecords_BLR_STRUC_T(:)=0
      SIZE_GEST_PANELS_L=0
      SIZE_GEST_PANELS_L_j1=0
      SIZE_VARIABLES_PANELS_L=0_8
      SIZE_VARIABLES_PANELS_L_j1=0_8
      SIZE_GEST_PANELS_U=0
      SIZE_GEST_PANELS_U_j1=0
      SIZE_VARIABLES_PANELS_U=0_8
      SIZE_VARIABLES_PANELS_U_j1=0_8
      SIZE_GEST_CB_LRB=0
      SIZE_GEST_CB_LRB_j1j2=0
      SIZE_VARIABLES_CB_LRB=0_8
      SIZE_VARIABLES_CB_LRB_j1j2=0_8
      SIZE_GEST_DIAG_BLOCKS=0
      SIZE_GEST_DIAG_BLOCKS_j1=0
      SIZE_VARIABLES_DIAG_BLOCKS=0_8
      SIZE_VARIABLES_DIAG_BLOCKS_j1=0_8
C     
C     BEGINNING OF THE MAIN LOOP ON ALL VARIABLES OF THE STRUCTURE
C
      DO i1=1,NBVARIABLES_BLR_STRUC_T 
         SELECT CASE(i1)
         CASE(B_IsSYM)
            NbRecords_BLR_STRUC_T(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES_BLR_STRUC_T(i1)=SIZE_LOGICAL
            elseif(mode.EQ.save_mode) then
               SIZE_VARIABLES_BLR_STRUC_T(i1)=SIZE_LOGICAL
               write(unit,iostat=err) BLR_STRUC%IsSYM
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES_BLR_STRUC_T(i1)=SIZE_LOGICAL
               read(unit,iostat=err) BLR_STRUC%IsSYM
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif 
               IF (INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(B_IsT2)
            NbRecords_BLR_STRUC_T(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES_BLR_STRUC_T(i1)=SIZE_LOGICAL
            elseif(mode.EQ.save_mode) then
               SIZE_VARIABLES_BLR_STRUC_T(i1)=SIZE_LOGICAL
               write(unit,iostat=err) BLR_STRUC%IsT2
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES_BLR_STRUC_T(i1)=SIZE_LOGICAL
               read(unit,iostat=err) BLR_STRUC%IsT2
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif 
               IF (INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(B_IsSLAVE)
            NbRecords_BLR_STRUC_T(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES_BLR_STRUC_T(i1)=SIZE_LOGICAL
            elseif(mode.EQ.save_mode) then
               SIZE_VARIABLES_BLR_STRUC_T(i1)=SIZE_LOGICAL
               write(unit,iostat=err) BLR_STRUC%IsSLAVE
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES_BLR_STRUC_T(i1)=SIZE_LOGICAL
               read(unit,iostat=err) BLR_STRUC%IsSLAVE
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif 
               IF (INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(B_BEGS_BLR_STATIC)
            NbRecords_BLR_STRUC_T(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(BLR_STRUC%BEGS_BLR_STATIC)) THEN
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=
     &                 size(BLR_STRUC%BEGS_BLR_STATIC,1)*SIZE_INT
               ELSE
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(BLR_STRUC%BEGS_BLR_STATIC)) THEN
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=
     &                 size(BLR_STRUC%BEGS_BLR_STATIC,1)*SIZE_INT
                  write(unit,iostat=err)
     &                 size(BLR_STRUC%BEGS_BLR_STATIC,1)
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) BLR_STRUC%BEGS_BLR_STATIC
               ELSE
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(BLR_STRUC%BEGS_BLR_STATIC)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=size_array1*SIZE_INT
                  allocate(BLR_STRUC%BEGS_BLR_STATIC(size_array1)
     &                 , stat=allocok)
                  if (allocok .GT. 0) THEN
                     INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,INFO(2))
                  endif
                  read(unit,iostat=err) BLR_STRUC%BEGS_BLR_STATIC
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(B_BEGS_BLR_DYNAMIC)
            NbRecords_BLR_STRUC_T(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(BLR_STRUC%BEGS_BLR_DYNAMIC)) THEN
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=
     &                 size(BLR_STRUC%BEGS_BLR_DYNAMIC,1)*SIZE_INT
               ELSE
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(BLR_STRUC%BEGS_BLR_DYNAMIC)) THEN
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=
     &                 size(BLR_STRUC%BEGS_BLR_DYNAMIC,1)*SIZE_INT
                  write(unit,iostat=err)
     &                 size(BLR_STRUC%BEGS_BLR_DYNAMIC,1)
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) BLR_STRUC%BEGS_BLR_DYNAMIC
               ELSE
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(BLR_STRUC%BEGS_BLR_DYNAMIC)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=size_array1*SIZE_INT
                  allocate(BLR_STRUC%BEGS_BLR_DYNAMIC(size_array1)
     &                 , stat=allocok)
                  if (allocok .GT. 0) THEN
                     INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,INFO(2))
                  endif
                  read(unit,iostat=err) BLR_STRUC%BEGS_BLR_DYNAMIC
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(B_BEGS_BLR_L)
            NbRecords_BLR_STRUC_T(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(BLR_STRUC%BEGS_BLR_L)) THEN
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=
     &                 size(BLR_STRUC%BEGS_BLR_L,1)*SIZE_INT
               ELSE
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(BLR_STRUC%BEGS_BLR_L)) THEN
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=
     &                 size(BLR_STRUC%BEGS_BLR_L,1)*SIZE_INT
                  write(unit,iostat=err)
     &                 size(BLR_STRUC%BEGS_BLR_L,1)
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) BLR_STRUC%BEGS_BLR_L
               ELSE
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(BLR_STRUC%BEGS_BLR_L)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=size_array1*SIZE_INT
                  allocate(BLR_STRUC%BEGS_BLR_L(size_array1)
     &                 , stat=allocok)
                  if (allocok .GT. 0) THEN
                     INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,INFO(2))
                  endif
                  read(unit,iostat=err) BLR_STRUC%BEGS_BLR_L
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(B_BEGS_BLR_COL)
            NbRecords_BLR_STRUC_T(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(BLR_STRUC%BEGS_BLR_COL)) THEN
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=
     &                 size(BLR_STRUC%BEGS_BLR_COL,1)*SIZE_INT
               ELSE
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(BLR_STRUC%BEGS_BLR_COL)) THEN
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=
     &                 size(BLR_STRUC%BEGS_BLR_COL,1)*SIZE_INT
                  write(unit,iostat=err)
     &                 size(BLR_STRUC%BEGS_BLR_COL,1)
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) BLR_STRUC%BEGS_BLR_COL
               ELSE
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(BLR_STRUC%BEGS_BLR_COL)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=size_array1*SIZE_INT
                  allocate(BLR_STRUC%BEGS_BLR_COL(size_array1)
     &                 , stat=allocok)
                  if (allocok .GT. 0) THEN
                     INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,INFO(2))
                  endif
                  read(unit,iostat=err) BLR_STRUC%BEGS_BLR_COL
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(B_NB_ACCESSES_INIT)
            NbRecords_BLR_STRUC_T(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES_BLR_STRUC_T(i1)=SIZE_INT
            elseif(mode.EQ.save_mode) then
               SIZE_VARIABLES_BLR_STRUC_T(i1)=SIZE_INT
               write(unit,iostat=err) BLR_STRUC%NB_ACCESSES_INIT
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES_BLR_STRUC_T(i1)=SIZE_INT
               read(unit,iostat=err) BLR_STRUC%NB_ACCESSES_INIT
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif 
               IF (INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(B_NB_PANELS)
            NbRecords_BLR_STRUC_T(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES_BLR_STRUC_T(i1)=SIZE_INT
            elseif(mode.EQ.save_mode) then
               SIZE_VARIABLES_BLR_STRUC_T(i1)=SIZE_INT
               write(unit,iostat=err) BLR_STRUC%NB_PANELS
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES_BLR_STRUC_T(i1)=SIZE_INT
               read(unit,iostat=err) BLR_STRUC%NB_PANELS
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif 
               IF (INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(B_PANELS_L)
            if(mode.EQ.memory_save_mode) then
               IF(associated(BLR_STRUC%PANELS_L)) THEN
C    NbRecords_BLR_STRUC_T = 1 for the size of the array
                  NbRecords_BLR_STRUC_T(i1)=1
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  DO j1=1,size(BLR_STRUC%PANELS_L,1)
C     Call save_memory on BLR_STRUC%PANELS_L(j1)
                     CALL DMUMPS_SAVE_RESTORE_BLR_PANEL(
     &                    BLR_STRUC%PANELS_L(j1)
     &                    ,unit,MYID,memory_save_mode
     &                    ,SIZE_GEST_PANELS_L_j1
     &                    ,SIZE_VARIABLES_PANELS_L_j1
     &                    ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &                    ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                    ,size_read,size_allocated,size_written
     &                    ,INFO)
                     SIZE_GEST_PANELS_L=SIZE_GEST_PANELS_L+
     &                    SIZE_GEST_PANELS_L_j1
                     SIZE_VARIABLES_PANELS_L=SIZE_VARIABLES_PANELS_L+
     &                    SIZE_VARIABLES_PANELS_L_j1
                  ENDDO
               ELSE
                  NbRecords_BLR_STRUC_T(i1)=2
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(BLR_STRUC%PANELS_L)) THEN
                  NbRecords_BLR_STRUC_T(i1)=1
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  write(unit,iostat=err)
     &                 size(BLR_STRUC%PANELS_L,1)
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  DO j1=1,size(BLR_STRUC%PANELS_L,1)
C     Call save on BLR_STRUC%PANELS_L(j1)
                     CALL DMUMPS_SAVE_RESTORE_BLR_PANEL(
     &                    BLR_STRUC%PANELS_L(j1)
     &                    ,unit,MYID,save_mode
     &                    ,SIZE_GEST_PANELS_L_j1
     &                    ,SIZE_VARIABLES_PANELS_L_j1
     &                    ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &                    ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                    ,size_read,size_allocated,size_written
     &                    ,INFO)                     
                     IF ( INFO(1) .LT. 0 ) GOTO 100
                  ENDDO
               ELSE
                  NbRecords_BLR_STRUC_T(i1)=2
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
               ENDIF 
            elseif(mode.EQ.restore_mode) then
               nullify(BLR_STRUC%PANELS_L)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  NbRecords_BLR_STRUC_T(i1)=2
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  read(unit,iostat=err) dummy
                  if(err.ne.0) THEN
                     INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                    ,INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
               else
C    NbRecords_BLR_STRUC_T = 1 for the size of the array
                  NbRecords_BLR_STRUC_T(i1)=1
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  allocate(BLR_STRUC%PANELS_L(size_array1)
     &                 , stat=allocok)
                  if (allocok .GT. 0) THEN
                     INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,INFO(2))
                  endif
                  DO j1=1,size_array1
C     call restore on BLR_STRUC%PANELS_L(j1)
                     CALL DMUMPS_SAVE_RESTORE_BLR_PANEL(
     &                    BLR_STRUC%PANELS_L(j1)
     &                    ,unit,MYID,restore_mode
     &                    ,SIZE_GEST_PANELS_L_j1
     &                    ,SIZE_VARIABLES_PANELS_L_j1
     &                    ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &                    ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                    ,size_read,size_allocated,size_written
     &                    ,INFO)
                     SIZE_GEST_PANELS_L=SIZE_GEST_PANELS_L+
     &                    SIZE_GEST_PANELS_L_j1
                     SIZE_VARIABLES_PANELS_L=SIZE_VARIABLES_PANELS_L+
     &                    SIZE_VARIABLES_PANELS_L_j1
                  ENDDO
               endif
            endif
            CASE(B_PANELS_U)
            if(mode.EQ.memory_save_mode) then
               IF(associated(BLR_STRUC%PANELS_U)) THEN
C    NbRecords_BLR_STRUC_T = 1 for the size of the array
                  NbRecords_BLR_STRUC_T(i1)=1
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  DO j1=1,size(BLR_STRUC%PANELS_U,1)
C     Call save_memory on BLR_STRUC%PANELS_U(j1)
                     CALL DMUMPS_SAVE_RESTORE_BLR_PANEL(
     &                    BLR_STRUC%PANELS_U(j1)
     &                    ,unit,MYID,memory_save_mode
     &                    ,SIZE_GEST_PANELS_U_j1
     &                    ,SIZE_VARIABLES_PANELS_U_j1
     &                    ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &                    ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                    ,size_read,size_allocated,size_written
     &                    ,INFO)
                     SIZE_GEST_PANELS_U=SIZE_GEST_PANELS_U+
     &                    SIZE_GEST_PANELS_U_j1
                     SIZE_VARIABLES_PANELS_U=SIZE_VARIABLES_PANELS_U+
     &                    SIZE_VARIABLES_PANELS_U_j1
                  ENDDO
               ELSE
                  NbRecords_BLR_STRUC_T(i1)=2
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(BLR_STRUC%PANELS_U)) THEN
                  NbRecords_BLR_STRUC_T(i1)=1
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  write(unit,iostat=err)
     &                 size(BLR_STRUC%PANELS_U,1)
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  DO j1=1,size(BLR_STRUC%PANELS_U,1)
C     Call save on BLR_STRUC%PANELS_U(j1)
                     CALL DMUMPS_SAVE_RESTORE_BLR_PANEL(
     &                    BLR_STRUC%PANELS_U(j1)
     &                    ,unit,MYID,save_mode
     &                    ,SIZE_GEST_PANELS_U_j1
     &                    ,SIZE_VARIABLES_PANELS_U_j1
     &                    ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &                    ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                    ,size_read,size_allocated,size_written
     &                    ,INFO)
                     IF ( INFO(1) .LT. 0 ) GOTO 100
                  ENDDO
               ELSE
                  NbRecords_BLR_STRUC_T(i1)=2
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
               ENDIF 
            elseif(mode.EQ.restore_mode) then
               nullify(BLR_STRUC%PANELS_U)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  NbRecords_BLR_STRUC_T(i1)=2
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  read(unit,iostat=err) dummy
                  if(err.ne.0) THEN
                     INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                    ,INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
               else
C    NbRecords_BLR_STRUC_T = 1 for the size of the array
                  NbRecords_BLR_STRUC_T(i1)=1
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  allocate(BLR_STRUC%PANELS_U(size_array1)
     &                 , stat=allocok)
                  if (allocok .GT. 0) THEN
                     INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,INFO(2))
                  endif
                  DO j1=1,size_array1
C     call restore on BLR_STRUC%PANELS_U(j1)
                     CALL DMUMPS_SAVE_RESTORE_BLR_PANEL(
     &                    BLR_STRUC%PANELS_U(j1)
     &                    ,unit,MYID,restore_mode
     &                    ,SIZE_GEST_PANELS_U_j1
     &                    ,SIZE_VARIABLES_PANELS_U_j1
     &                    ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &                    ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                    ,size_read,size_allocated,size_written
     &                    ,INFO)
                     SIZE_GEST_PANELS_U=SIZE_GEST_PANELS_U+
     &                    SIZE_GEST_PANELS_U_j1
                     SIZE_VARIABLES_PANELS_U=SIZE_VARIABLES_PANELS_U+
     &                    SIZE_VARIABLES_PANELS_U_j1
                  ENDDO
               endif
            endif
         CASE(B_CB_LRB)
            if(mode.EQ.memory_save_mode) then
               IF(associated(BLR_STRUC%CB_LRB)) THEN
C    NbRecords_BLR_STRUC_T = 1 for the sizes of the array
                  NbRecords_BLR_STRUC_T(i1)=1
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  DO j1=1,size(BLR_STRUC%CB_LRB,1)
                     DO j2=1,size(BLR_STRUC%CB_LRB,2)
C     Call save_memory on BLR_STRUC%CB_LRB(j1,j2)
                        CALL DMUMPS_SAVE_RESTORE_LRB(
     &                       BLR_STRUC%CB_LRB(j1,j2)
     &                       ,unit,MYID,memory_save_mode
     &                       ,SIZE_GEST_CB_LRB_j1j2
     &                       ,SIZE_VARIABLES_CB_LRB_j1j2
     &                       ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &                       ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                       ,size_read,size_allocated,size_written
     &                       ,INFO)
                        SIZE_GEST_CB_LRB=SIZE_GEST_CB_LRB+
     &                       SIZE_GEST_CB_LRB_j1j2
                        SIZE_VARIABLES_CB_LRB=SIZE_VARIABLES_CB_LRB+
     &                       SIZE_VARIABLES_CB_LRB_j1j2
                     ENDDO
                  ENDDO
               ELSE
                  NbRecords_BLR_STRUC_T(i1)=2
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*3
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(BLR_STRUC%CB_LRB)) THEN
                  NbRecords_BLR_STRUC_T(i1)=1
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  write(unit,iostat=err)
     &                 size(BLR_STRUC%CB_LRB,1),size(BLR_STRUC%CB_LRB,2)
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  DO j1=1,size(BLR_STRUC%CB_LRB,1)
                     DO j2=1,size(BLR_STRUC%CB_LRB,2)
C     Call save on BLR_STRUC%PANELS_U(j1,j2)
                        CALL DMUMPS_SAVE_RESTORE_LRB(
     &                       BLR_STRUC%CB_LRB(j1,j2)
     &                       ,unit,MYID,save_mode
     &                       ,SIZE_GEST_CB_LRB_j1j2
     &                       ,SIZE_VARIABLES_CB_LRB_j1j2
     &                       ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &                       ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                       ,size_read,size_allocated,size_written
     &                       ,INFO)
                        IF ( INFO(1) .LT. 0 ) GOTO 100
                     ENDDO
                  ENDDO
               ELSE
                  NbRecords_BLR_STRUC_T(i1)=2
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*3
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
               ENDIF 
            elseif(mode.EQ.restore_mode) then
               nullify(BLR_STRUC%CB_LRB)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  NbRecords_BLR_STRUC_T(i1)=2
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*3
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  read(unit,iostat=err) dummy
                  if(err.ne.0) THEN
                     INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                    ,INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
               else
C     NbRecords_BLR_STRUC_T = 1 for the sizes of the array
                  NbRecords_BLR_STRUC_T(i1)=1
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  allocate(BLR_STRUC%CB_LRB(size_array1,size_array2)
     &                 , stat=allocok)
                  if (allocok .GT. 0) THEN
                     INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,INFO(2))
                  endif
                  DO j1=1,size_array1
                     DO j2=1,size_array2
C     call restore on BLR_STRUC%CB_LRB(j1,j2)
                        CALL DMUMPS_SAVE_RESTORE_LRB(
     &                       BLR_STRUC%CB_LRB(j1,j2)
     &                       ,unit,MYID,restore_mode
     &                       ,SIZE_GEST_CB_LRB_j1j2
     &                       ,SIZE_VARIABLES_CB_LRB_j1j2
     &                       ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &                       ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                       ,size_read,size_allocated,size_written
     &                       ,INFO)
                        SIZE_GEST_CB_LRB=SIZE_GEST_CB_LRB+
     &                       SIZE_GEST_CB_LRB_j1j2
                        SIZE_VARIABLES_CB_LRB=SIZE_VARIABLES_CB_LRB+
     &                       SIZE_VARIABLES_CB_LRB_j1j2
                     ENDDO
                  ENDDO
               endif
            endif
         CASE(B_DIAG_BLOCKS)
            if(mode.EQ.memory_save_mode) then
               IF(associated(BLR_STRUC%DIAG_BLOCKS)) THEN
C     NbRecords_BLR_STRUC_T = 1 for the size of the array
                  NbRecords_BLR_STRUC_T(i1)=1
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  DO j1=1,size(BLR_STRUC%DIAG_BLOCKS,1)
                     CALL DMUMPS_SAVE_RESTORE_DIAG_BLOCK(
     &                    BLR_STRUC%DIAG_BLOCKS(j1)
     &                    ,unit,MYID,memory_save_mode
     &                    ,SIZE_GEST_DIAG_BLOCKS_j1
     &                    ,SIZE_VARIABLES_DIAG_BLOCKS_j1
     &                    ,SIZE_INT, SIZE_ARITH_DEP
     &                    ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                    ,size_read,size_allocated,size_written
     &                    ,INFO)
                     SIZE_GEST_DIAG_BLOCKS=SIZE_GEST_DIAG_BLOCKS+
     &                    SIZE_GEST_DIAG_BLOCKS_j1
                     SIZE_VARIABLES_DIAG_BLOCKS=
     &                    SIZE_VARIABLES_DIAG_BLOCKS+
     &                    SIZE_VARIABLES_DIAG_BLOCKS_j1
                  ENDDO
               ELSE
                  NbRecords_BLR_STRUC_T(i1)=2
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(BLR_STRUC%DIAG_BLOCKS)) THEN
                  NbRecords_BLR_STRUC_T(i1)=1
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  write(unit,iostat=err)
     &                 size(BLR_STRUC%DIAG_BLOCKS,1)
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  DO j1=1,size(BLR_STRUC%DIAG_BLOCKS,1)
                     CALL DMUMPS_SAVE_RESTORE_DIAG_BLOCK(
     &                    BLR_STRUC%DIAG_BLOCKS(j1)
     &                    ,unit,MYID,save_mode
     &                    ,SIZE_GEST_DIAG_BLOCKS_j1
     &                    ,SIZE_VARIABLES_DIAG_BLOCKS_j1
     &                    ,SIZE_INT, SIZE_ARITH_DEP
     &                    ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                    ,size_read,size_allocated,size_written
     &                    ,INFO)
                     IF ( INFO(1) .LT. 0 ) GOTO 100
                  ENDDO
               ELSE
                  NbRecords_BLR_STRUC_T(i1)=2
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
               ENDIF 
            elseif(mode.EQ.restore_mode) then
               nullify(BLR_STRUC%DIAG_BLOCKS)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  NbRecords_BLR_STRUC_T(i1)=2
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  read(unit,iostat=err) dummy
                  if(err.ne.0) THEN
                     INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                    ,INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
               else
C     NbRecords_BLR_STRUC_T = 1 for the size of the array
                  NbRecords_BLR_STRUC_T(i1)=1
                  SIZE_GEST_BLR_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_BLR_STRUC_T(i1)=0
                  allocate(BLR_STRUC%DIAG_BLOCKS(size_array1)
     &                 , stat=allocok)
                  if (allocok .GT. 0) THEN
                     INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,INFO(2))
                  endif
                  DO j1=1,size_array1
                     CALL DMUMPS_SAVE_RESTORE_DIAG_BLOCK(
     &                    BLR_STRUC%DIAG_BLOCKS(j1)
     &                    ,unit,MYID,restore_mode
     &                    ,SIZE_GEST_DIAG_BLOCKS_j1
     &                    ,SIZE_VARIABLES_DIAG_BLOCKS_j1
     &                    ,SIZE_INT, SIZE_ARITH_DEP
     &                    ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                    ,size_read,size_allocated,size_written
     &                    ,INFO)
                     SIZE_GEST_DIAG_BLOCKS=SIZE_GEST_DIAG_BLOCKS+
     &                    SIZE_GEST_DIAG_BLOCKS_j1
                     SIZE_VARIABLES_DIAG_BLOCKS=
     &                    SIZE_VARIABLES_DIAG_BLOCKS+
     &                    SIZE_VARIABLES_DIAG_BLOCKS_j1
                  ENDDO
               endif
            endif
         CASE(B_NFS4FATHER)
C     We should not save/restore NFS4FATHER
            NbRecords_BLR_STRUC_T(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES_BLR_STRUC_T(i1)=SIZE_INT
            elseif(mode.EQ.save_mode) then
               SIZE_VARIABLES_BLR_STRUC_T(i1)=SIZE_INT
               write(unit,iostat=err) BLR_STRUC%NFS4FATHER
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES_BLR_STRUC_T(i1)=SIZE_INT
               read(unit,iostat=err) BLR_STRUC%NFS4FATHER
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif 
               IF (INFO(1) .LT. 0 ) GOTO 100
            endif
            CASE(B_M_ARRAY)
C     Nothing to save nor restore. M_ARRAY stays in facto
               if(mode.EQ.restore_mode) then
                  nullify(BLR_STRUC%M_ARRAY)
               endif
         CASE DEFAULT
         END SELECT
         if(mode.EQ.memory_save_mode) then
C     If the size to write (SIZE_VARIABLES_BLR_STRUC_T(i1)) is greater than 2^31
C     Subrecords are created which need to be taken into account in
C     the file size computation
            NbSubRecords=int(SIZE_VARIABLES_BLR_STRUC_T(i1)/huge(I4))
            IF(NbSubRecords.GT.0) then
               NbRecords_BLR_STRUC_T(i1)=NbRecords_BLR_STRUC_T(i1)
     &              +NbSubRecords
            ENDIF
         elseif(mode.EQ.save_mode) then
            size_written=size_written+SIZE_VARIABLES_BLR_STRUC_T(i1)
     &           +int(SIZE_GEST_BLR_STRUC_T(i1),kind=8)
#if defined(MUMPS_NOF2003)
            size_written=size_written
     &           +int(2*SIZE_INT*NbRecords_BLR_STRUC_T(i1),kind=8)
#endif
         elseif(mode.EQ.restore_mode) then
            size_allocated=size_allocated+
     &           SIZE_VARIABLES_BLR_STRUC_T(i1)
            size_read=size_read+SIZE_VARIABLES_BLR_STRUC_T(i1)
     &           +int(SIZE_GEST_BLR_STRUC_T(i1),kind=8)
#if defined(MUMPS_NOF2003)
            size_read=size_read
     &           +int(2*SIZE_INT*NbRecords_BLR_STRUC_T(i1),kind=8)
#endif
         endif
      ENDDO
      if(mode.EQ.memory_save_mode) then
         Local_SIZE_VARIABLES=sum(SIZE_VARIABLES_BLR_STRUC_T)
     &        +SIZE_VARIABLES_PANELS_L
     &        +SIZE_VARIABLES_PANELS_U
     &        +SIZE_VARIABLES_CB_LRB
     &        +SIZE_VARIABLES_DIAG_BLOCKS
         Local_SIZE_GEST=sum(SIZE_GEST_BLR_STRUC_T)
     &        +SIZE_GEST_PANELS_L
     &        +SIZE_GEST_PANELS_U
     &        +SIZE_GEST_CB_LRB
     &        +SIZE_GEST_DIAG_BLOCKS
#if defined(MUMPS_NOF2003)
         Local_NbRecords=sum(NbRecords_BLR_STRUC_T)
         Local_SIZE_GEST=Local_SIZE_GEST+2*SIZE_INT*Local_NbRecords
#endif
      endif
 100  continue
      RETURN
      END SUBROUTINE DMUMPS_SAVE_RESTORE_BLR_STRUC
!--------------------------------- DMUMPS_SAVE_RESTORE_LRB -----------------------------!
      SUBROUTINE DMUMPS_SAVE_RESTORE_LRB(LRB_T
     &     ,unit,MYID,mode
     &     ,Local_SIZE_GEST,Local_SIZE_VARIABLES
     &     ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &     ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &     ,size_read,size_allocated,size_written
     &     ,INFO)
      INCLUDE 'mpif.h'
C  =======
C  Purpose
C  =======
C
C     This routine is designed to manage a LRB_TYPE structure (save, restore, compute memory)
C     
C  ==========
C  Parameters
C  ==========
C
C     LRB_T               : TYPE (LRB_TYPE) : the main structure
C
C     unit                : The unit of the file to be written or read
C
C     mode                : the type of operation to be performed by the routine 
C                           memory_save  = compute the size of the save file and of the structure
C                           save         = save the instace
C                           restore      = restore the instace
C
C     Local_SIZE_GEST     : The management size of the structure (array sizes, records...)
C
C     Local_SIZE_VARIABLE : The size in memory of the structure
C
C     SIZE_INT            : size of an integer      
C
C     SIZE_ARITH_DEP      : size of an arithmetic dependent variable
C
C     SIZE_LOGICAL        : size of a logical
C
C     INFO1/INFO2         : copies of of INFO(1) and INFO(2) to allow save/restore of failled instaces
C
      TYPE(LRB_TYPE) :: LRB_T
      INTEGER,intent(IN)::unit,MYID
      INTEGER,intent(IN) :: mode
      INTEGER,INTENT(OUT) :: Local_SIZE_GEST
      INTEGER(8),intent(OUT) :: Local_SIZE_VARIABLES
      INTEGER,intent(INOUT):: INFO(2)
      INTEGER,intent(IN):: SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
      INTEGER(8),intent(IN) :: TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
      INTEGER(8),intent(INOUT):: size_read,size_allocated,size_written
C LRB_TYPE fields are Q,R,K,M,N,ISLR:
      INTEGER, PARAMETER :: LRB_Q=1
      INTEGER, PARAMETER :: LRB_R=2
      INTEGER, PARAMETER :: LRB_K=3
      INTEGER, PARAMETER :: LRB_M=4
      INTEGER, PARAMETER :: LRB_N=5
      INTEGER, PARAMETER :: LRB_ISLR=6
      INTEGER, PARAMETER :: LRB_K_CALC=7
      INTEGER, PARAMETER :: LRB_K_STOC=8
      INTEGER, PARAMETER :: LRB_Q_CALC=9
      INTEGER, PARAMETER :: LRB_R_CALC=10
      INTEGER, PARAMETER :: LRB_Q_STOC=11
      INTEGER, PARAMETER :: LRB_R_STOC=12
      INTEGER, PARAMETER :: NBVARIABLES_LRB_TYPE=12
      INTEGER(8),dimension(NBVARIABLES_LRB_TYPE)::
     &     SIZE_VARIABLES_LRB_TYPE
      INTEGER,dimension(NBVARIABLES_LRB_TYPE)::
     &     SIZE_GEST_LRB_TYPE
      INTEGER,dimension(NBVARIABLES_LRB_TYPE)::
     &     NbRecords_LRB_TYPE
      INTEGER:: size_array1,size_array2,dummy,allocok
      INTEGER:: err,i1,NbSubRecords
#if defined(MUMPS_NOF2003)
      INTEGER:: Local_NbRecords
#endif
      INTEGER(4) ::I4
      SIZE_VARIABLES_LRB_TYPE(:)=0_8
      SIZE_GEST_LRB_TYPE(:)=0
      NbRecords_LRB_TYPE(:)=0
      DO i1=1,NBVARIABLES_LRB_TYPE
         SELECT CASE(i1)
         CASE(LRB_Q)
            NbRecords_LRB_TYPE(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(LRB_T%Q)) THEN
                  SIZE_GEST_LRB_TYPE(i1)=SIZE_INT*2
                  SIZE_VARIABLES_LRB_TYPE(i1)=
     &                 size(LRB_T%Q,1)*size(LRB_T%Q,2)
     &                 *     SIZE_ARITH_DEP
               ELSE
                  SIZE_GEST_LRB_TYPE(i1)=SIZE_INT*3
                  SIZE_VARIABLES_LRB_TYPE(i1)=0
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(LRB_T%Q)) THEN
                  SIZE_GEST_LRB_TYPE(i1)=SIZE_INT*2
                  SIZE_VARIABLES_LRB_TYPE(i1)=
     &                 size(LRB_T%Q,1)*size(LRB_T%Q,2)
     &                 *     SIZE_ARITH_DEP
                  write(unit,iostat=err) size(LRB_T%Q,1),size(LRB_T%Q,2)
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 300
                  write(unit,iostat=err) LRB_T%Q
               ELSE
                  SIZE_GEST_LRB_TYPE(i1)=SIZE_INT*3
                  SIZE_VARIABLES_LRB_TYPE(i1)=0
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 300
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 300
            elseif(mode.EQ.restore_mode) then
               nullify(LRB_T%Q)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 300
               if(size_array1.EQ.-999) then
                  SIZE_GEST_LRB_TYPE(i1)=SIZE_INT*3
                  SIZE_VARIABLES_LRB_TYPE(i1)=0
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST_LRB_TYPE(i1)=SIZE_INT*2
                  SIZE_VARIABLES_LRB_TYPE(i1)=
     &                 size_array1*size_array2*SIZE_ARITH_DEP
                  allocate(LRB_T%Q(size_array1,size_array2),
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,INFO(2))
                  endif
                  read(unit,iostat=err) LRB_T%Q
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 300
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 300
            endif
         CASE(LRB_R)
            NbRecords_LRB_TYPE(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(LRB_T%R)) THEN
                  SIZE_GEST_LRB_TYPE(i1)=SIZE_INT*2
                  SIZE_VARIABLES_LRB_TYPE(i1)=
     &                 size(LRB_T%R,1)*size(LRB_T%R,2)
     &                 *     SIZE_ARITH_DEP
               ELSE
                  SIZE_GEST_LRB_TYPE(i1)=SIZE_INT*3
                  SIZE_VARIABLES_LRB_TYPE(i1)=0
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(LRB_T%R)) THEN
                  SIZE_GEST_LRB_TYPE(i1)=SIZE_INT*2
                  SIZE_VARIABLES_LRB_TYPE(i1)=
     &                 size(LRB_T%R,1)*size(LRB_T%R,2)
     &                 *     SIZE_ARITH_DEP
                  write(unit,iostat=err) size(LRB_T%R,1),size(LRB_T%R,2)
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 300
                  write(unit,iostat=err) LRB_T%R
               ELSE
                  SIZE_GEST_LRB_TYPE(i1)=SIZE_INT*3
                  SIZE_VARIABLES_LRB_TYPE(i1)=0
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 300
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 300
            elseif(mode.EQ.restore_mode) then
               nullify(LRB_T%R)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 300
               if(size_array1.EQ.-999) then
                  SIZE_GEST_LRB_TYPE(i1)=SIZE_INT*3
                  SIZE_VARIABLES_LRB_TYPE(i1)=0
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST_LRB_TYPE(i1)=SIZE_INT*2
                  SIZE_VARIABLES_LRB_TYPE(i1)=
     &                 size_array1*size_array2*SIZE_ARITH_DEP
                  allocate(LRB_T%R(size_array1,size_array2),
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,INFO(2))
                  endif
                  read(unit,iostat=err) LRB_T%R
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 300
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 300
            endif
         CASE(LRB_K)
            NbRecords_LRB_TYPE(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES_LRB_TYPE(i1)=SIZE_INT
            elseif(mode.EQ.save_mode) then
               SIZE_VARIABLES_LRB_TYPE(i1)=SIZE_INT
               write(unit,iostat=err) LRB_T%K
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 300
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES_LRB_TYPE(i1)=SIZE_INT
               read(unit,iostat=err) LRB_T%K
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif 
               IF ( INFO(1) .LT. 0 ) GOTO 300
            endif
         CASE(LRB_M)
            NbRecords_LRB_TYPE(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES_LRB_TYPE(i1)=SIZE_INT
            elseif(mode.EQ.save_mode) then
               SIZE_VARIABLES_LRB_TYPE(i1)=SIZE_INT
               write(unit,iostat=err) LRB_T%M
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 300
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES_LRB_TYPE(i1)=SIZE_INT
               read(unit,iostat=err) LRB_T%M
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif 
               IF ( INFO(1) .LT. 0 ) GOTO 300
            endif
         CASE(LRB_N)
            NbRecords_LRB_TYPE(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES_LRB_TYPE(i1)=SIZE_INT
            elseif(mode.EQ.save_mode) then
               SIZE_VARIABLES_LRB_TYPE(i1)=SIZE_INT
               write(unit,iostat=err) LRB_T%N
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 300
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES_LRB_TYPE(i1)=SIZE_INT
               read(unit,iostat=err) LRB_T%N
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif 
               IF ( INFO(1) .LT. 0 ) GOTO 300
            endif
         CASE(LRB_ISLR)
            NbRecords_LRB_TYPE(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES_LRB_TYPE(i1)=SIZE_LOGICAL
            elseif(mode.EQ.save_mode) then
               SIZE_VARIABLES_LRB_TYPE(i1)=SIZE_LOGICAL
               write(unit,iostat=err) LRB_T%ISLR
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 300
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES_LRB_TYPE(i1)=SIZE_LOGICAL
               read(unit,iostat=err) LRB_T%ISLR
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif 
               IF ( INFO(1) .LT. 0 ) GOTO 300
            endif   
         CASE (LRB_K_CALC,LRB_Q_CALC,LRB_R_CALC)
C          Do nothing, we only keep LRB_*_STOC* from 1 phase to another
C          This may change in the future if, for example, we compute
C          CALC versions of the blocks at the first solve and keep
C          them for subsequent solves.
            NbRecords_LRB_TYPE(i1)=0
            SIZE_GEST_LRB_TYPE(i1)=0
            SIZE_VARIABLES_LRB_TYPE(i1)=0
         CASE(LRB_K_STOC)
            NbRecords_LRB_TYPE(i1)=2 ! array size + array
            if(mode.EQ.memory_save_mode) then
               IF(associated(LRB_T%K_STOC)) THEN
                 SIZE_GEST_LRB_TYPE(i1)=SIZE_INT
                 SIZE_VARIABLES_LRB_TYPE(i1)=size(LRB_T%K_STOC)*SIZE_INT
               ELSE
                 SIZE_GEST_LRB_TYPE(i1)=SIZE_INT*2
                 SIZE_VARIABLES_LRB_TYPE(i1)=0
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(LRB_T%K_STOC)) THEN
                  SIZE_GEST_LRB_TYPE(i1)=SIZE_INT
                  SIZE_VARIABLES_LRB_TYPE(i1)=size(LRB_T%K_STOC)*
     &                                        SIZE_INT
                  write(unit,iostat=err) size(LRB_T%K_STOC)
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 300
                  write(unit,iostat=err) LRB_T%K_STOC
               ELSE
                  SIZE_GEST_LRB_TYPE(i1)=SIZE_INT*2
                  SIZE_VARIABLES_LRB_TYPE(i1)=0
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 300
                  write(unit,iostat=err) -999
               ENDIF
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 300
            elseif(mode.EQ.restore_mode) then
               nullify(LRB_T%K_STOC)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 300
               if(size_array1.EQ.-999) then
                  SIZE_GEST_LRB_TYPE(i1)=SIZE_INT*2
                  SIZE_VARIABLES_LRB_TYPE(i1)=0
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST_LRB_TYPE(i1)=SIZE_INT
                  SIZE_VARIABLES_LRB_TYPE(i1)=size_array1*SIZE_INT
                  allocate(LRB_T%K_STOC(size_array1),stat=allocok)
                  if (allocok .GT. 0) THEN
                     INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 300
                  read(unit,iostat=err) LRB_T%K_STOC
               endif
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 300
            endif
         CASE(LRB_Q_STOC)
             CALL DMUMPS_SAVE_RESTORE_Q_OR_R_STOC(LRB_T%Q_STOC)
             IF ( INFO(1) .LT. 0 ) GOTO 300
         CASE(LRB_R_STOC)
             CALL DMUMPS_SAVE_RESTORE_Q_OR_R_STOC(LRB_T%R_STOC)
             IF ( INFO(1) .LT. 0 ) GOTO 300
         CASE DEFAULT
         END SELECT
         if(mode.EQ.memory_save_mode) then
C     If the size to write (SIZE_VARIABLES_BLR_STRUC_T(i1)) is greater than 2^31
C     Subrecords are created which need to be taken into account in
C     the file size computation
            NbSubRecords=int(SIZE_VARIABLES_LRB_TYPE(i1)/huge(I4))
            IF(NbSubRecords.GT.0) then
               NbRecords_LRB_TYPE(i1)=
     &              NbRecords_LRB_TYPE(i1)
     &              +NbSubRecords
            ENDIF
         elseif(mode.EQ.save_mode) then
            size_written=size_written+SIZE_VARIABLES_LRB_TYPE(i1)
     &           +int(SIZE_GEST_LRB_TYPE(i1),kind=8)
#if defined(MUMPS_NOF2003)
            size_written=size_written
     &           +int(2*SIZE_INT*NbRecords_LRB_TYPE(i1),kind=8)
#endif
         elseif(mode.EQ.restore_mode) then
            size_allocated=size_allocated+
     &           SIZE_VARIABLES_LRB_TYPE(i1)
            size_read=size_read+SIZE_VARIABLES_LRB_TYPE(i1)
     &           +int(SIZE_GEST_LRB_TYPE(i1),kind=8)
#if defined(MUMPS_NOF2003)
            size_read=size_read
     &           +int(2*SIZE_INT*NbRecords_LRB_TYPE(i1),kind=8)
#endif
         endif
      ENDDO
      if(mode.EQ.memory_save_mode) then
         Local_SIZE_VARIABLES=sum(SIZE_VARIABLES_LRB_TYPE)
         Local_SIZE_GEST=sum(SIZE_GEST_LRB_TYPE)
#if defined(MUMPS_NOF2003)
         Local_NbRecords=sum(NbRecords_LRB_TYPE)
         Local_SIZE_GEST=Local_SIZE_GEST+2*SIZE_INT*Local_NbRecords
#endif
      endif
 300  continue
      RETURN
      CONTAINS
      SUBROUTINE DMUMPS_SAVE_RESTORE_Q_OR_R_STOC(QorR_STOC)
      IMPLICIT NONE
      type(LRMIXED_TYPE), POINTER :: QorR_STOC(:)
      integer j1
      if(mode.eq.memory_save_mode) then
        IF(associated(QorR_STOC)) THEN
C         NbRecords_BLR_LRB_TYPE = 1 for the size of the array
C         SIZE_VARIABLES_LRB_TYPE = 0 because everything is in the routine
          NbRecords_LRB_TYPE(i1)=1 ! array size
          SIZE_GEST_LRB_TYPE(i1)=SIZE_INT
          SIZE_VARIABLES_LRB_TYPE(i1)=0
C         Loop on the precisions
          DO j1=1,size(QorR_STOC,1)
C           Count two records (size + contents or twice -999)
            NbRecords_LRB_TYPE(i1)=NbRecords_LRB_TYPE(i1)+2
            IF (associated(QorR_STOC(j1)%val)) THEN
C             Count 1 GEST int for size + size(QorR_STOC(j1)%val) ints
              SIZE_GEST_LRB_TYPE(i1)=SIZE_GEST_LRB_TYPE(i1)+SIZE_INT
              SIZE_VARIABLES_LRB_TYPE(i1)= SIZE_VARIABLES_LRB_TYPE(i1)+
     &                            size(QorR_STOC(j1)%val)*SIZE_INT
            ELSE ! QorR_STOC(j)%val not associated
C             Count two ints: -999 then -999 again at next line
              SIZE_GEST_LRB_TYPE(i1)=SIZE_GEST_LRB_TYPE(i1)+2*SIZE_INT
              SIZE_VARIABLES_LRB_TYPE(i1)=SIZE_VARIABLES_LRB_TYPE(i1)+0
            ENDIF
          ENDDO
        ELSE ! QorR_STOC is not associated, count 2 records and 2 ints
          NbRecords_LRB_TYPE(i1)=2
          SIZE_GEST_LRB_TYPE(i1)=SIZE_INT*2
          SIZE_VARIABLES_LRB_TYPE(i1)=0
        ENDIF
      else if (mode.eq.save_mode) then
        IF (associated(QorR_STOC)) THEN
          NbRecords_LRB_TYPE(i1)=1 ! array size
          SIZE_GEST_LRB_TYPE(i1)=SIZE_INT
          SIZE_VARIABLES_LRB_TYPE(i1)=0
          write(unit,iostat=err) size(QorR_STOC)
          if(err.ne.0) then
            INFO(1) = -72
            CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &           INFO(2))
          endif
          IF ( INFO(1) .LT. 0 ) GOTO 300
C         Loop on the precisions
          DO j1=1,size(QorR_STOC,1)
            NbRecords_LRB_TYPE(i1)=NbRecords_LRB_TYPE(i1)+2
            IF (associated(QorR_STOC(j1)%val)) THEN
C             2 records (size+array itself)
              SIZE_GEST_LRB_TYPE(i1)=SIZE_GEST_LRB_TYPE(i1)+SIZE_INT
              SIZE_VARIABLES_LRB_TYPE(i1)=SIZE_VARIABLES_LRB_TYPE(i1)+
     &                       size(QorR_STOC(j1)%val)*SIZE_INT
              write(unit,iostat=err) size(QorR_STOC(j1)%val)
              if(err.ne.0) then
                INFO(1) = -72
                CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &               INFO(2))
              endif
              IF ( INFO(1) .LT. 0 ) GOTO 300
              write(unit,iostat=err) QorR_STOC(j1)%val
              if(err.ne.0) then
                INFO(1) = -72
                CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &               INFO(2))
              endif
              IF ( INFO(1) .LT. 0 ) GOTO 300
            ELSE ! QorR_STOC(j1)%val not associated
C             Count two records and two ints
              SIZE_GEST_LRB_TYPE(i1)=SIZE_GEST_LRB_TYPE(i1)+2*SIZE_INT
              SIZE_VARIABLES_LRB_TYPE(i1)=SIZE_VARIABLES_LRB_TYPE(i1)+0
              write(unit,iostat=err) -999
              if(err.ne.0) then
                INFO(1) = -72
                CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &               INFO(2))
              endif
              write(unit,iostat=err) -999
              if(err.ne.0) then
                INFO(1) = -72
                CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &               INFO(2))
              endif
            ENDIF
          ENDDO
        ELSE ! QorR_STOC is not associated, count 2 records and 2 ints
          NbRecords_LRB_TYPE(i1)=2
          SIZE_GEST_LRB_TYPE(i1)=SIZE_INT*2
          SIZE_VARIABLES_LRB_TYPE(i1)=0
          write(unit,iostat=err) -999
          if(err.ne.0) then
             INFO(1) = -72
             CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &            INFO(2))
          endif
          IF ( INFO(1) .LT. 0 ) GOTO 300
          write(unit,iostat=err) -999
        ENDIF
        if(err.ne.0) then
           INFO(1) = -72
           CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &          INFO(2))
        endif
      else if (mode.eq.restore_mode) then
        nullify(QorR_STOC)
        read(unit,iostat=err) size_array1
        if(err.ne.0) THEN
          INFO(1) = -75
          CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &         ,INFO(2))
        endif
        IF ( INFO(1) .LT. 0 ) GOTO 300
        if (size_array1.eq.-999) then
          NbRecords_LRB_TYPE(i1)=2
          SIZE_GEST_LRB_TYPE(i1)=SIZE_INT*2
          SIZE_VARIABLES_LRB_TYPE(i1)=0
          read(unit,iostat=err) dummy
          if(err.ne.0) THEN
            INFO(1) = -75
            CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &           ,INFO(2))
          endif
        else
          NbRecords_LRB_TYPE(i1)=1 ! array size
          SIZE_GEST_LRB_TYPE(i1)=SIZE_INT
          SIZE_VARIABLES_LRB_TYPE(i1)=0
          allocate(QorR_STOC(size_array1),stat=allocok)
          if (allocok .GT. 0) THEN
            INFO(1) = -78
            CALL MUMPS_SETI8TOI4(
     &      TOTAL_STRUC_SIZE-size_allocated
     &      ,INFO(2))
          endif
          IF ( INFO(1) .LT. 0 ) GOTO 300
C         Loop on the precisions to read QorR_STOC(j1) arrays
          DO j1=1,size(QorR_STOC,1)
            nullify(QorR_STOC(j1)%val)
C           2 records per precision
            NbRecords_LRB_TYPE(i1)=NbRecords_LRB_TYPE(i1)+2
            read(unit,iostat=err) size_array2
            if(err.ne.0) THEN
               INFO(1) = -75
               CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &              ,INFO(2))
            endif
            IF (INFO(1).LT.0) GOTO 300
            IF (size_array2 .eq. -999) THEN
C             Count two ints
              SIZE_GEST_LRB_TYPE(i1)=SIZE_GEST_LRB_TYPE(i1)+2*SIZE_INT
              read(unit,iostat=err) dummy
              if(err.ne.0) THEN
                 INFO(1) = -75
                 CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                ,INFO(2))
              endif
              IF (INFO(1).LT.0) GOTO 300
            ELSE
              ALLOCATE(QorR_STOC(j1)%val(size_array2),stat=allocok)
              if (allocok .GT. 0) THEN
                INFO(1) = -78
                CALL MUMPS_SETI8TOI4(
     &          TOTAL_STRUC_SIZE-size_allocated
     &          ,INFO(2))
              endif
              IF ( INFO(1) .LT. 0 ) GOTO 300
C             Count 1 int + size_array2 ints for the actual array
              SIZE_GEST_LRB_TYPE(i1)=SIZE_GEST_LRB_TYPE(i1)+SIZE_INT
              SIZE_VARIABLES_LRB_TYPE(i1)= SIZE_VARIABLES_LRB_TYPE(i1)+
     &                                size(QorR_STOC(j1)%val)*SIZE_INT
              read(unit,iostat=err) QorR_STOC(j1)%val
              if(err.ne.0) THEN
                 INFO(1) = -75
                 CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                ,INFO(2))
              endif
              IF (INFO(1).LT.0) GOTO 300
            ENDIF
            if(err.ne.0) THEN
              INFO(1) = -75
              CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &             ,INFO(2))
            endif
            IF ( INFO(1) .LT. 0 ) GOTO 300
          ENDDO
        endif
      endif
 300  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_SAVE_RESTORE_Q_OR_R_STOC
      END SUBROUTINE DMUMPS_SAVE_RESTORE_LRB
!--------------------------------- DMUMPS_SAVE_RESTORE_BLR_PANEL -----------------------------!
      SUBROUTINE DMUMPS_SAVE_RESTORE_BLR_PANEL(BLR_PANEL_T
     &     ,unit,MYID,mode
     &     ,Local_SIZE_GEST,Local_SIZE_VARIABLES
     &     ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &     ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &     ,size_read,size_allocated,size_written
     &     ,INFO)
      INCLUDE 'mpif.h'
C  =======
C  Purpose
C  =======
C
C     This routine is designed to manage a blr_panel_type structure
C     (save, restore, compute memory)
C     
C  ==========
C  Parameters
C  ==========
C
C     BLR_PANEL_T         : TYPE (blr_panel_type): the main structure
C
C     unit                : The unit of the file to be written or read
C
C     mode                : the type of operation to be performed by the routine 
C                           memory_save  = compute the size of the save file and of the structure
C                           save         = save the instace
C                           restore      = restore the instace
C
C     Local_SIZE_GEST     : The management size of the structure (array sizes, records...)
C
C     Local_SIZE_VARIABLE : The size in memory of the structure
C
C     SIZE_INT            : size of an integer      
C
C     SIZE_ARITH_DEP      : size of an arithmetic dependent variable
C
C     SIZE_LOGICAL        : size of a logical
C
C     INFO1/INFO2         : copies of of INFO(1) and INFO(2) to allow save/restore of failled instaces
C
      TYPE(blr_panel_type) :: BLR_PANEL_T
      INTEGER,intent(IN)::unit,MYID
      INTEGER,intent(IN) :: mode
      INTEGER,INTENT(OUT) :: Local_SIZE_GEST
      INTEGER(8),intent(OUT) :: Local_SIZE_VARIABLES
      INTEGER,intent(INOUT):: INFO(2)
      INTEGER,intent(IN):: SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
      INTEGER(8),intent(IN) :: TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
      INTEGER(8),intent(INOUT):: size_read,size_allocated,size_written
C     blr_panel_type contains NB_ACCESSES_LEFT and LRB_PANEL:
      INTEGER, PARAMETER :: P_NB_ACCESSES_LEFT=1
      INTEGER, PARAMETER :: P_LRB_PANEL=2
      INTEGER, PARAMETER :: NBVARIABLES_BLR_PANEL_TYPE = 2
      INTEGER(8),dimension(NBVARIABLES_BLR_PANEL_TYPE)::
     &     SIZE_VARIABLES_BLR_PANEL_TYPE
      INTEGER,dimension(NBVARIABLES_BLR_PANEL_TYPE)::
     &     SIZE_GEST_BLR_PANEL_TYPE
      INTEGER,dimension(NBVARIABLES_BLR_PANEL_TYPE)::
     &     NbRecords_BLR_PANEL_TYPE
      INTEGER:: size_array1,dummy,allocok
      INTEGER:: err,i1,j1,NbSubRecords
#if defined(MUMPS_NOF2003)
      INTEGER:: Local_NbRecords
#endif
      INTEGER:: SIZE_GEST_LRB_PANEL_j1,SIZE_GEST_LRB_PANEL
      INTEGER(8)::SIZE_VARIABLES_LRB_PANEL_j1,SIZE_VARIABLES_LRB_PANEL
      INTEGER(4)::I4
      SIZE_VARIABLES_BLR_PANEL_TYPE(:)=0_8
      SIZE_GEST_BLR_PANEL_TYPE(:)=0
      NbRecords_BLR_PANEL_TYPE(:)=0
      SIZE_GEST_LRB_PANEL_j1=0
      SIZE_GEST_LRB_PANEL=0
      SIZE_VARIABLES_LRB_PANEL_j1=0_8
      SIZE_VARIABLES_LRB_PANEL=0_8
      DO i1=1,NBVARIABLES_BLR_PANEL_TYPE
         SELECT CASE(i1)
         CASE(P_NB_ACCESSES_LEFT)
            NbRecords_BLR_PANEL_TYPE(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES_BLR_PANEL_TYPE(i1)=SIZE_INT
            elseif(mode.EQ.save_mode) then
               SIZE_VARIABLES_BLR_PANEL_TYPE(i1)=SIZE_INT
               write(unit,iostat=err) BLR_PANEL_T%NB_ACCESSES_LEFT
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 400
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES_BLR_PANEL_TYPE(i1)=SIZE_INT
               read(unit,iostat=err) BLR_PANEL_T%NB_ACCESSES_LEFT
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif 
               IF ( INFO(1) .LT. 0 ) GOTO 400
            endif
         CASE(P_LRB_PANEL)
            if(mode.EQ.memory_save_mode) then
               IF(associated(BLR_PANEL_T%LRB_PANEL)) THEN
C     NbRecords_BLR_PANEL_TYPE = 1 for the size of the array
                  NbRecords_BLR_PANEL_TYPE(i1)=1
                  SIZE_GEST_BLR_PANEL_TYPE(i1)=SIZE_INT
C     SIZE_VARIABLES_BLR_PANEL_TYPE = 0 because everything is in the routine
                  SIZE_VARIABLES_BLR_PANEL_TYPE(i1)=0
                  DO j1=1,size(BLR_PANEL_T%LRB_PANEL,1)
C     Call save_memory on BLR_PANEL_T%LRB_PANEL(j1)
                     CALL DMUMPS_SAVE_RESTORE_LRB(
     &                    BLR_PANEL_T%LRB_PANEL(j1)
     &                    ,unit,MYID,memory_save_mode
     &                    ,SIZE_GEST_LRB_PANEL_j1
     &                    ,SIZE_VARIABLES_LRB_PANEL_j1
     &                    ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &                    ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                    ,size_read,size_allocated,size_written
     &                    ,INFO)
                     SIZE_GEST_LRB_PANEL=SIZE_GEST_LRB_PANEL+
     &                    SIZE_GEST_LRB_PANEL_j1
                     SIZE_VARIABLES_LRB_PANEL=SIZE_VARIABLES_LRB_PANEL+
     &                    SIZE_VARIABLES_LRB_PANEL_j1
                  ENDDO
               ELSE
                  NbRecords_BLR_PANEL_TYPE(i1)=2
                  SIZE_GEST_BLR_PANEL_TYPE(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_PANEL_TYPE(i1)=0
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(BLR_PANEL_T%LRB_PANEL)) THEN
C     NbRecords_BLR_PANEL_TYPE = 1 for the size of the array
                  NbRecords_BLR_PANEL_TYPE(i1)=1
                  SIZE_GEST_BLR_PANEL_TYPE(i1)=SIZE_INT
C     SIZE_VARIABLES_BLR_PANEL_TYPE = 0 because everything is in the routine
                  SIZE_VARIABLES_BLR_PANEL_TYPE(i1)=0
                  write(unit,iostat=err) size(BLR_PANEL_T%LRB_PANEL,1)
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 400
                  DO j1=1,size(BLR_PANEL_T%LRB_PANEL,1)
C     Call save on BLR_PANEL_T%LRB_PANEL(j1)
                     CALL DMUMPS_SAVE_RESTORE_LRB(
     &                    BLR_PANEL_T%LRB_PANEL(j1)
     &                    ,unit,MYID,save_mode
     &                    ,SIZE_GEST_LRB_PANEL_j1
     &                    ,SIZE_VARIABLES_LRB_PANEL_j1
     &                    ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &                    ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                    ,size_read,size_allocated,size_written
     &                    ,INFO)
                     IF ( INFO(1) .LT. 0 ) GOTO 400
                  ENDDO
               ELSE
                  NbRecords_BLR_PANEL_TYPE(i1)=2
                  SIZE_GEST_BLR_PANEL_TYPE(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_PANEL_TYPE(i1)=0
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 400
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 400
               ENDIF 
            elseif(mode.EQ.restore_mode) then
               nullify(BLR_PANEL_T%LRB_PANEL)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 400
               if(size_array1.EQ.-999) then
                  NbRecords_BLR_PANEL_TYPE(i1)=2
                  SIZE_GEST_BLR_PANEL_TYPE(i1)=SIZE_INT*2
                  SIZE_VARIABLES_BLR_PANEL_TYPE(i1)=0
                  read(unit,iostat=err) dummy
                  if(err.ne.0) THEN
                     INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                    ,INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 400
               else
C     NbRecords_BLR_PANEL_TYPE = 1 for the size of the array                  
                  NbRecords_BLR_PANEL_TYPE(i1)=1
                  SIZE_GEST_BLR_PANEL_TYPE(i1)=SIZE_INT
C     SIZE_VARIABLES_BLR_PANEL_TYPE = 0 because everything is in the routine                  
                  SIZE_VARIABLES_BLR_PANEL_TYPE(i1)=0
                  allocate(BLR_PANEL_T%LRB_PANEL(size_array1)
     &                 , stat=allocok)
                  if (allocok .GT. 0) THEN
                     INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,INFO(2))
                  endif
                  DO j1=1,size_array1
C     call restore on BLR_PANEL_T%LRB_PANEL(j1)
                     CALL DMUMPS_SAVE_RESTORE_LRB(
     &                    BLR_PANEL_T%LRB_PANEL(j1)
     &                    ,unit,MYID,restore_mode
     &                    ,SIZE_GEST_LRB_PANEL_j1
     &                    ,SIZE_VARIABLES_LRB_PANEL_j1
     &                    ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &                    ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                    ,size_read,size_allocated,size_written
     &                    ,INFO)
                     SIZE_GEST_LRB_PANEL=SIZE_GEST_LRB_PANEL+
     &                    SIZE_GEST_LRB_PANEL_j1
                     SIZE_VARIABLES_LRB_PANEL=SIZE_VARIABLES_LRB_PANEL+
     &                    SIZE_VARIABLES_LRB_PANEL_j1
                  ENDDO
               endif
            endif
         CASE DEFAULT
         END SELECT
         if(mode.EQ.memory_save_mode) then
C     If the size to write (SIZE_VARIABLES_BLR_PANEL_TYPE(i1)) is greater than 2^31
C     Subrecords are created which need to be taken into account in
C     the file size computation
            NbSubRecords=int(SIZE_VARIABLES_BLR_PANEL_TYPE(i1)/huge(I4))
            IF(NbSubRecords.GT.0) then
               NbRecords_BLR_PANEL_TYPE(i1)=
     &              NbRecords_BLR_PANEL_TYPE(i1)
     &              +NbSubRecords
            ENDIF
         elseif(mode.EQ.save_mode) then
            size_written=size_written+SIZE_VARIABLES_BLR_PANEL_TYPE(i1)
     &           +int(SIZE_GEST_BLR_PANEL_TYPE(i1),kind=8)
#if defined(MUMPS_NOF2003)
            size_written=size_written
     &           +int(2*SIZE_INT*NbRecords_BLR_PANEL_TYPE(i1),kind=8)
#endif
         elseif(mode.EQ.restore_mode) then
            size_allocated=size_allocated+
     &           SIZE_VARIABLES_BLR_PANEL_TYPE(i1)
            size_read=size_read+SIZE_VARIABLES_BLR_PANEL_TYPE(i1)
     &           +int(SIZE_GEST_BLR_PANEL_TYPE(i1),kind=8)
#if defined(MUMPS_NOF2003)
            size_read=size_read
     &           +int(2*SIZE_INT*NbRecords_BLR_PANEL_TYPE(i1),kind=8)
#endif
         endif
      ENDDO
      if(mode.EQ.memory_save_mode) then
         Local_SIZE_VARIABLES=sum(SIZE_VARIABLES_BLR_PANEL_TYPE)+
     &        SIZE_VARIABLES_LRB_PANEL
         Local_SIZE_GEST=sum(SIZE_GEST_BLR_PANEL_TYPE)+
     &        SIZE_GEST_LRB_PANEL
#if defined(MUMPS_NOF2003)
         Local_NbRecords=sum(NbRecords_BLR_PANEL_TYPE)
         Local_SIZE_GEST=Local_SIZE_GEST+2*SIZE_INT*Local_NbRecords
#endif
      endif
 400  continue
      RETURN
      END SUBROUTINE DMUMPS_SAVE_RESTORE_BLR_PANEL
!--------------------------------- DMUMPS_SAVE_RESTORE_DIAG_BLOCK -----------------------------!
      SUBROUTINE DMUMPS_SAVE_RESTORE_DIAG_BLOCK(DIAG_BLOCK_T
     &     ,unit,MYID,mode
     &     ,Local_SIZE_GEST,Local_SIZE_VARIABLES
     &     ,SIZE_INT, SIZE_ARITH_DEP
     &     ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &     ,size_read,size_allocated,size_written
     &     ,INFO)
      INCLUDE 'mpif.h'
C  =======
C  Purpose
C  =======
C
C     This routine is designed to manage a diag_block_type structure
C     (save, restore, compute memory)
C     
C  ==========
C  Parameters
C  ==========
C
C     DIAG_BLOCK          : TYPE (diag_block_type) : the main structure
C
C     unit                : The unit of the file to be written or read
C
C     mode                : the type of operation to be performed by the routine 
C                           memory_save  = compute the size of the save file and of the structure
C                           save         = save the instace
C                           restore      = restore the instace
C
C     Local_SIZE_GEST     : The management size of the structure (array sizes, records...)
C
C     Local_SIZE_VARIABLE : The size in memory of the structure
C
C     SIZE_INT            : size of an integer      
C
C     SIZE_ARITH_DEP      : size of an arithmetic dependent variable
C
C     INFO1/INFO2         : copies of of INFO(1) and INFO(2) to allow save/restore of failled instaces
C
      TYPE(diag_block_type) :: DIAG_BLOCK_T
      INTEGER,intent(IN)::unit,MYID
      INTEGER,intent(IN) :: mode
      INTEGER,INTENT(OUT) :: Local_SIZE_GEST
      INTEGER(8),intent(OUT) :: Local_SIZE_VARIABLES
      INTEGER,intent(INOUT):: INFO(2)
      INTEGER,intent(IN):: SIZE_INT, SIZE_ARITH_DEP
      INTEGER(8),intent(IN) :: TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
      INTEGER(8),intent(INOUT):: size_read,size_allocated,size_written
C     diag_block_type contains only one field "DIAG_BLOCK"
      INTEGER, PARAMETER :: D_DIAG_BLOCK=1
#if !defined(NOMIXHPC)
      INTEGER, PARAMETER :: D_DIAG_STOC=2
      INTEGER, PARAMETER :: D_DIAG_SIZE=3
      INTEGER, PARAMETER :: NBVARIABLES_DIAG_BLOCK_TYPE = 3
#else
      INTEGER, PARAMETER :: NBVARIABLES_DIAG_BLOCK_TYPE = 1
#endif
      INTEGER(8),dimension(NBVARIABLES_DIAG_BLOCK_TYPE)::
     &     SIZE_VARIABLES_DIAG_BLOCK_TYPE
      INTEGER,dimension(NBVARIABLES_DIAG_BLOCK_TYPE)::
     &     SIZE_GEST_DIAG_BLOCK_TYPE
      INTEGER,dimension(NBVARIABLES_DIAG_BLOCK_TYPE)::
     &     NbRecords_DIAG_BLOCK_TYPE
      INTEGER:: size_array1,dummy,allocok
      INTEGER:: err,i1,NbSubRecords
#if defined(MUMPS_NOF2003)
      INTEGER:: Local_NbRecords
#endif
      INTEGER(4) :: I4
      SIZE_VARIABLES_DIAG_BLOCK_TYPE(:)=0_8
      SIZE_GEST_DIAG_BLOCK_TYPE(:)=0
      NbRecords_DIAG_BLOCK_TYPE(:)=0
      DO i1=1,NBVARIABLES_DIAG_BLOCK_TYPE
         SELECT CASE(i1)
         CASE(D_DIAG_BLOCK)
            NbRecords_DIAG_BLOCK_TYPE(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(DIAG_BLOCK_T%DIAG_BLOCK)) THEN
                  SIZE_GEST_DIAG_BLOCK_TYPE(i1)=SIZE_INT
                  SIZE_VARIABLES_DIAG_BLOCK_TYPE(i1)=
     &                 size(DIAG_BLOCK_T%DIAG_BLOCK,1)
     &                 *     SIZE_ARITH_DEP
               ELSE
                  SIZE_GEST_DIAG_BLOCK_TYPE(i1)=SIZE_INT*2
                  SIZE_VARIABLES_DIAG_BLOCK_TYPE(i1)=0
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(DIAG_BLOCK_T%DIAG_BLOCK)) THEN
                  SIZE_GEST_DIAG_BLOCK_TYPE(i1)=SIZE_INT
                  SIZE_VARIABLES_DIAG_BLOCK_TYPE(i1)=
     &                 size(DIAG_BLOCK_T%DIAG_BLOCK,1)
     &                 *     SIZE_ARITH_DEP
                  write(unit,iostat=err) size(DIAG_BLOCK_T%DIAG_BLOCK,1)
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 200
                  write(unit,iostat=err) DIAG_BLOCK_T%DIAG_BLOCK
               ELSE
                  SIZE_GEST_DIAG_BLOCK_TYPE(i1)=SIZE_INT*2
                  SIZE_VARIABLES_DIAG_BLOCK_TYPE(i1)=0
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 200
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 200
            elseif(mode.EQ.restore_mode) then
               nullify(DIAG_BLOCK_T%DIAG_BLOCK)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 200
               if(size_array1.EQ.-999) then
                  SIZE_GEST_DIAG_BLOCK_TYPE(i1)=SIZE_INT*2
                  SIZE_VARIABLES_DIAG_BLOCK_TYPE(i1)=0
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST_DIAG_BLOCK_TYPE(i1)=SIZE_INT
                  SIZE_VARIABLES_DIAG_BLOCK_TYPE(i1)=
     &                 size_array1*SIZE_ARITH_DEP
                  allocate(DIAG_BLOCK_T%DIAG_BLOCK(size_array1),
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,INFO(2))
                     GOTO 200
                  endif
                  read(unit,iostat=err) DIAG_BLOCK_T%DIAG_BLOCK
               endif
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
                  GOTO 200
               endif
            endif
C{
#if !defined(NOMIXHPC)
C{
         CASE(D_DIAG_STOC)
C          Count two records (size + contents or twice -999)
           NbRecords_DIAG_BLOCK_TYPE(i1)=2
           if(mode.eq.memory_save_mode) then
             SIZE_VARIABLES_DIAG_BLOCK_TYPE(i1)=0
             IF (associated(DIAG_BLOCK_T%DIAG_STOC%val)) THEN
C              Count 1 GEST int for size + size(DIAG_BLOCK_T%DIAG_STOC%val) ints
               SIZE_GEST_DIAG_BLOCK_TYPE(i1) = SIZE_INT
               SIZE_VARIABLES_DIAG_BLOCK_TYPE(i1) = 
     &             size(DIAG_BLOCK_T%DIAG_STOC%val)*SIZE_INT
             ELSE ! DIAG_BLOCK_T(j)%val not associated
C              Count two ints: -999 then -999 again at next line
               SIZE_GEST_DIAG_BLOCK_TYPE(i1) =
     &           SIZE_GEST_DIAG_BLOCK_TYPE(i1) + 2*SIZE_INT
               SIZE_VARIABLES_DIAG_BLOCK_TYPE(i1) = 0
             ENDIF
           else if (mode.eq.save_mode) then
             IF (associated(DIAG_BLOCK_T%DIAG_STOC%val)) THEN
C              2 records (size+array itself)
               SIZE_GEST_DIAG_BLOCK_TYPE(i1) = SIZE_INT 
               SIZE_VARIABLES_DIAG_BLOCK_TYPE(i1) =
     &           size(DIAG_BLOCK_T%DIAG_STOC%val)*SIZE_INT
               write(unit,iostat=err) size(DIAG_BLOCK_T%DIAG_STOC%val)
               if(err.ne.0) then
                 INFO(1) = -72
                 CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 200
               write(unit,iostat=err) DIAG_BLOCK_T%DIAG_STOC%val
               if(err.ne.0) then
                 INFO(1) = -72
                 CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 200
             ELSE ! DIAG_BLOCK_T%DIAG_STOC%val not associated
C              Count two records and two ints
               SIZE_GEST_DIAG_BLOCK_TYPE(i1) = 2*SIZE_INT
               SIZE_VARIABLES_DIAG_BLOCK_TYPE(i1) = 0
               write(unit,iostat=err) -999
               if(err.ne.0) then
                 INFO(1) = -72
                 CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                INFO(2))
               endif
               write(unit,iostat=err) -999
               if(err.ne.0) then
                 INFO(1) = -72
                 CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 200
             ENDIF
             if(err.ne.0) then
                INFO(1) = -72
                CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &               INFO(2))
             endif
           else if (mode.eq.restore_mode) then
             read(unit,iostat=err) size_array1
             if(err.ne.0) THEN
               INFO(1) = -75
               CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &              ,INFO(2))
             endif
             IF ( INFO(1) .LT. 0 ) GOTO 200
             if (size_array1.eq.-999) then
               SIZE_GEST_DIAG_BLOCK_TYPE(i1)=SIZE_INT*2
               SIZE_VARIABLES_DIAG_BLOCK_TYPE(i1)=0
               read(unit,iostat=err) dummy
               if(err.ne.0) THEN
                 INFO(1) = -75
                 CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                ,INFO(2))
               endif
             else
C              Count one int for the size of the array
               SIZE_GEST_DIAG_BLOCK_TYPE(i1)=SIZE_INT
               SIZE_VARIABLES_DIAG_BLOCK_TYPE(i1)=0
               allocate(DIAG_BLOCK_T%DIAG_STOC%VAL(size_array1),
     &                  stat=allocok)
               if (allocok .GT. 0) THEN
                 INFO(1) = -78
                 CALL MUMPS_SETI8TOI4(
     &           TOTAL_STRUC_SIZE-size_allocated
     &           ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 200
C              Count size_array1 ints for the actual array
               SIZE_VARIABLES_DIAG_BLOCK_TYPE(i1) =
     &            SIZE_VARIABLES_DIAG_BLOCK_TYPE(i1) +
     &            size(DIAG_BLOCK_T%DIAG_STOC%val)*SIZE_INT
               read(unit,iostat=err) DIAG_BLOCK_T%DIAG_STOC%val
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF (INFO(1).LT.0) GOTO 200
             endif
           endif
         CASE(D_DIAG_SIZE)
             NbRecords_DIAG_BLOCK_TYPE(i1)=1
             SIZE_VARIABLES_DIAG_BLOCK_TYPE(i1)=SIZE_INT
             if(mode.EQ.memory_save_mode) then
             elseif(mode.EQ.save_mode) then
                write(unit,iostat=err) DIAG_BLOCK_T%DIAG_SIZE
                if(err.ne.0) then
                   INFO(1) = -72
                   CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                  INFO(2))
                endif
             elseif(mode.EQ.restore_mode) then
                read(unit,iostat=err) DIAG_BLOCK_T%DIAG_SIZE
                if(err.ne.0) THEN
                   INFO(1) = -75
                   CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                  ,INFO(2))
                endif 
             endif   
             IF ( INFO(1) .LT. 0 ) GOTO 200
C} MIXHPC
#endif
C} LRMIXED
         CASE DEFAULT
         END SELECT
         if(mode.EQ.memory_save_mode) then
C     If the size to write (SIZE_VARIABLES_BLR_STRUC_T(i1)) is greater than 2^31
C     Subrecords are created which need to be taken into account in
C     the file size computation
            NbSubRecords=int(SIZE_VARIABLES_DIAG_BLOCK_TYPE(i1)/
     &      huge(I4))
            IF(NbSubRecords.GT.0) then
               NbRecords_DIAG_BLOCK_TYPE(i1)=
     &              NbRecords_DIAG_BLOCK_TYPE(i1)
     &              +NbSubRecords
            ENDIF
         elseif(mode.EQ.save_mode) then
            size_written=size_written+SIZE_VARIABLES_DIAG_BLOCK_TYPE(i1)
     &           +int(SIZE_GEST_DIAG_BLOCK_TYPE(i1),kind=8)
#if defined(MUMPS_NOF2003)
            size_written=size_written
     &           +int(2*SIZE_INT*NbRecords_DIAG_BLOCK_TYPE(i1),kind=8)
#endif
         elseif(mode.EQ.restore_mode) then
            size_allocated=size_allocated+
     &           SIZE_VARIABLES_DIAG_BLOCK_TYPE(i1)
            size_read=size_read+SIZE_VARIABLES_DIAG_BLOCK_TYPE(i1)
     &           +int(SIZE_GEST_DIAG_BLOCK_TYPE(i1),kind=8)
#if defined(MUMPS_NOF2003)
            size_read=size_read
     &           +int(2*SIZE_INT*NbRecords_DIAG_BLOCK_TYPE(i1),kind=8)
#endif
         endif
      ENDDO
      if(mode.EQ.memory_save_mode) then
         Local_SIZE_VARIABLES=sum(SIZE_VARIABLES_DIAG_BLOCK_TYPE)
         Local_SIZE_GEST=sum(SIZE_GEST_DIAG_BLOCK_TYPE)
#if defined(MUMPS_NOF2003)
         Local_NbRecords=sum(NbRecords_DIAG_BLOCK_TYPE)
         Local_SIZE_GEST=Local_SIZE_GEST+2*SIZE_INT*Local_NbRecords
#endif
      endif
 200  continue
      RETURN
      END SUBROUTINE DMUMPS_SAVE_RESTORE_DIAG_BLOCK
#endif
      END MODULE DMUMPS_LR_DATA_M
