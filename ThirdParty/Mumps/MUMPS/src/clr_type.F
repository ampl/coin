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
      MODULE CMUMPS_LR_TYPE
      USE MUMPS_LR_COMMON, only: NPREC_STOC_MAX, NPREC_CALC_MAX,
     &    NPREC_STOC, PREC_DBL, PREC_56, PREC_48, PREC_40,
     &    PREC_SGL, PREC_24, PREC_16, 
     &    PREC_FR, PREC_DIAG_FR, CLOSE_PREC,
     &    PREC_FAC, ! TODO add prec_sol and modify sol_lr consequently
     &    QUANTIZE,
     &    PRECISIONS, COSTS_PREC,
     &    NPREC_CALC, CALC2STOC
      IMPLICIT NONE
C  -------------------------------
C   type LRB_TYPE
C  -------------------------------
C ISLR   = false : the block is not LR
C        = true  : the block is LR
C 
C Q      = if ISLR : first k columns of the Q factor from the 
C          RRQR facto of the original block B, size M * K
C          else the M * N block 
C 
C R      = if ISLR : first k rows of the R factor. the permutation P
C          from the RRQR facto is already applied so that
C          R is rectangular of size K * N
C          else not allocated
C 
C K = the rank, computed with truncated RRQR
C M = nb rows
C N = nb columns
C 
C  Note that nullifying pointers at initialization is
C  necessary to be able to test correctly if they
C  need to be deallocated.
C     We could consider using iso_c_binding :: c_int8 integer kind
C     instead, to avoid rounding to K34 in places
      TYPE LRMIXED_TYPE
        INTEGER, POINTER :: val(:) => null()
      END TYPE LRMIXED_TYPE
      TYPE LRB_TYPE
        COMPLEX,POINTER,DIMENSION(:,:) :: Q => null()
        COMPLEX,POINTER,DIMENSION(:,:) :: R => null()
        INTEGER :: K,M,N
        LOGICAL :: ISLR
        INTEGER, POINTER :: K_CALC(:)  => null()
        INTEGER, POINTER :: K_STOC(:)  => null()
        type(LRMIXED_TYPE), POINTER :: Q_CALC(:,:) => null()
        type(LRMIXED_TYPE), POINTER :: R_CALC(:,:) => null()
        type(LRMIXED_TYPE), POINTER :: Q_STOC(:) => null()
        type(LRMIXED_TYPE), POINTER :: R_STOC(:) => null()
        ! this holds
        ! dot_product(k_stoc(2:NPREC_STOC),
        !             costs_prec(2:NPREC_STOC))
        integer :: kw_stoc = 0
      END TYPE LRB_TYPE
      CONTAINS
      SUBROUTINE CMUMPS_INIT_LRMIXED_PRECISIONS(
     &      K148, K451, K456, K457, DK8)
C       This routine initializes the global variables related to
C       LRMIXED.
C  - NPREC_STOC is inherited from facto and describe how many precisions are
C    used at most to store a given tile
C  - NPREC_CALC describes the number of precisions that *could* be used at any
C    point during factorization or solve. Currently most hardware support double
C    and single precision computation.
C  - These variables do not dictate the behaviour of MUMPS, KEEPs (and ICNTLs)
C    dictate how these variables can be interpreted (e.g. K448 to rely on
C    multiple precisions during solve phase)
        INTEGER, intent(in) :: K148, K451, K456, K457
        REAL, intent(in) :: DK8
C ------------ LOCAL VARIABLES --------------------
        INTEGER :: PREC
C       Initialize all global variables to 0, then set those that are
C       effectively used (depending on both KEEPs and MUMPS_ARITH) to
C       1,2,3,...
        PREC_DBL = 0
        PREC_56  = 0
        PREC_48  = 0
        PREC_40  = 0
        PREC_SGL = 0
        PREC_24  = 0
        PREC_16  = 0
        ! The following constants are initialized with 1 which is
        ! assumed to be the working precisions
        PREC_FR  = 1
        PREC_DIAG_FR  = 1
        CLOSE_PREC = 1
        QUANTIZE = .FALSE.
        PREC_FAC = 1
        CALC2STOC = 0
        IF (K451.EQ.0) THEN
C         Use same precisions for storage as for computations
          NPREC_STOC = 1
          PREC_SGL = 1
          NPREC_CALC = 1
          CALC2STOC(1) = PREC_SGL
C         Except when you want only one other arithmetic
          if (k456.gt.0) then
            select case(k456)
              case(3)
                NPREC_STOC = 2
                PREC_24=2
              case(2)
                NPREC_STOC = 2
                PREC_16=2
            end select
          endif
        ELSEIF (K451.EQ.1) THEN
C         Use all multiple of 8 bits          
          NPREC_STOC = 3
          PREC_SGL = 1
          PREC_24  = 2
          PREC_16  = 3
          NPREC_CALC = 1
          CALC2STOC(1) = PREC_SGL
        ENDIF
C       Push all stoc/calc for PREC_DBL which will get PRECISIONS = 0
        IF (K148.EQ.2) THEN
          PREC_DBL = 1
          NPREC_STOC = NPREC_STOC + 1
          IF (PREC_56.GT.0)  PREC_56  = PREC_56 + 1
          IF (PREC_48.GT.0)  PREC_48  = PREC_48 + 1
          IF (PREC_40.GT.0)  PREC_40  = PREC_40 + 1
          IF (PREC_SGL.GT.0) PREC_SGL = PREC_SGL + 1
          IF (PREC_24.GT.0)  PREC_24  = PREC_24 + 1
          IF (PREC_16.GT.0)  PREC_16  = PREC_16 + 1
          DO PREC=NPREC_CALC,1,-1
            IF(CALC2STOC(PREC).EQ.PREC_SGL-1)
     &        CALC2STOC(PREC+1) = PREC_SGL
          END DO
          NPREC_CALC = NPREC_CALC + 1
          CALC2STOC(1) = PREC_DBL
        ENDIF
        DO PREC=1,NPREC_STOC
                IF (PREC.EQ.PREC_DBL) THEN
                        COSTS_PREC(PREC) = 8
                        PRECISIONS(PREC) = 2.0D0**(-53)
                ELSEIF (PREC.EQ.PREC_56) THEN
                        COSTS_PREC(PREC) = 7
                        PRECISIONS(PREC) = 2.0D0**(-45)
                ELSEIF (PREC.EQ.PREC_48) THEN
                        COSTS_PREC(PREC) = 6
                        PRECISIONS(PREC) = 2.0D0**(-37)
                ELSEIF (PREC.EQ.PREC_40) THEN
                        COSTS_PREC(PREC) = 5
                        PRECISIONS(PREC) = 2.0D0**(-29)
                ELSEIF (PREC.EQ.PREC_SGL) THEN
                        COSTS_PREC(PREC) = 4
                        PRECISIONS(PREC) = 2.0D0**(-24)
                ELSEIF (PREC.EQ.PREC_24) THEN
                        COSTS_PREC(PREC) = 3
                        PRECISIONS(PREC) = 2.0D0**(-16)
                ELSEIF (PREC.EQ.PREC_16) THEN
                        COSTS_PREC(PREC) = 2
                        PRECISIONS(PREC) = 2.0D0**(-8)
                ENDIF
        END DO
        CALC2STOC(NPREC_CALC+1) = NPREC_STOC+1
        PRECISIONS(NPREC_STOC+1) = 1.0D0 
        COSTS_PREC(:) = 2 * COSTS_PREC(:)
        IF (K148.EQ.2) THEN
C         This bypasses check 0.0D0 <= \rho_i < prec(1)
C         because singular value \rho_i > 0
          PRECISIONS(1) = 0.0D0
          PREC_FR = 2
          PREC_FAC = 2
          PREC_DIAG_FR = 2
          IF (K456.GT.0) THEN
              PREC_DIAG_FR = K456
          ENDIF
#if defined(NOMIXHPC)
          PREC_DIAG_FR = 1
          PREC_FR = 1
#endif
        ELSE  IF (K456.GT.0) THEN
C         AJ maybe this is too bold of a check because it is specific to
C         Dmumps+mixed (w/o custom!)
C         This bypasses check prec(1) <= \rho_i < prec(2) 
C         because singular value \rho_i > 0
          PRECISIONS(2) = 0.0D0
        ENDIF
        DO PREC=1,NPREC_STOC
           IF (DK8.gt.PRECISIONS(PREC)) CLOSE_PREC = PREC
        ENDDO
        IF (BTEST(K457,0)) PREC_FR = CLOSE_PREC
        IF (BTEST(K457,1)) PREC_DIAG_FR = CLOSE_PREC
        IF (BTEST(K457,2)) QUANTIZE = .TRUE.
      END SUBROUTINE CMUMPS_INIT_LRMIXED_PRECISIONS
      SUBROUTINE DEALLOC_LRB(LRB_OUT, KEEP8, K34
     &                       , KEEP_UNIFORM_ARG, KEEP_MIXED_ARG
     &                       )
        TYPE(LRB_TYPE), INTENT(INOUT) :: LRB_OUT
        INTEGER(8) :: KEEP8(150)
        INTEGER :: K34
        LOGICAL, OPTIONAL :: KEEP_UNIFORM_ARG, KEEP_MIXED_ARG
C     local variables
        INTEGER :: MEM, IDUMMY, JDUMMY
        LOGICAL :: KEEP_UNIFORM, KEEP_MIXED
        INTEGER :: PREC
C       Skip empty blocks        
        IF (LRB_OUT%M.EQ.0) RETURN
        IF (LRB_OUT%N.EQ.0) RETURN
        IF(present(KEEP_UNIFORM_ARG)) THEN
          KEEP_UNIFORM = KEEP_UNIFORM_ARG
        ELSE
          KEEP_UNIFORM = .FALSE.
        ENDIF
        IF(present(KEEP_MIXED_ARG)) THEN
          KEEP_MIXED = KEEP_MIXED_ARG
        ELSE
          KEEP_MIXED = .FALSE.
        ENDIF
C       Update memory consumption and deallocate
        MEM = 0
        IF (.NOT.KEEP_UNIFORM) THEN
          IF (LRB_OUT%ISLR) THEN
            IF (associated(LRB_OUT%Q)) THEN
              MEM = MEM + size(LRB_OUT%Q)
              DEALLOCATE (LRB_OUT%Q)
              NULLIFY(LRB_OUT%Q)
            ENDIF
            IF (associated(LRB_OUT%R)) THEN
              MEM = MEM + size(LRB_OUT%R)
              DEALLOCATE (LRB_OUT%R)
              NULLIFY(LRB_OUT%R)
            ENDIF
          ELSE
            IF (associated(LRB_OUT%Q)) THEN
              MEM = MEM + size(LRB_OUT%Q)
              DEALLOCATE (LRB_OUT%Q)
              NULLIFY(LRB_OUT%Q)
            ENDIF
          ENDIF
        ENDIF
        IF (.NOT.KEEP_MIXED) THEN
          IF (associated(LRB_OUT%Q_STOC)) THEN
            DO PREC=1,NPREC_STOC
              IF(associated(LRB_OUT%Q_STOC(PREC)%val))  THEN
                MEM = MEM + (size(LRB_OUT%Q_STOC(PREC)%val)*K34 
     &          + COSTS_PREC(1)-1) / COSTS_PREC(1)
                DEALLOCATE(LRB_OUT%Q_STOC(PREC)%val)
                NULLIFY(LRB_OUT%Q_STOC(PREC)%val)
              ENDIF
            ENDDO
            DEALLOCATE (LRB_OUT%Q_STOC)
            NULLIFY(LRB_OUT%Q_STOC)
          ENDIF
          IF (associated(LRB_OUT%R_STOC)) THEN
            DO PREC=1,NPREC_STOC
              IF(associated(LRB_OUT%R_STOC(PREC)%val))  THEN
                MEM = MEM + (size(LRB_OUT%R_STOC(PREC)%val)*K34 
     &          + COSTS_PREC(1)-1) / COSTS_PREC(1)
                DEALLOCATE(LRB_OUT%R_STOC(PREC)%val)
                NULLIFY(LRB_OUT%R_STOC(PREC)%val)
              ENDIF
            ENDDO
            DEALLOCATE(LRB_OUT%R_STOC)
            NULLIFY(LRB_OUT%R_STOC)
          ENDIF
          IF (associated(LRB_OUT%K_STOC)) THEN
            DEALLOCATE(LRB_OUT%K_STOC)
            NULLIFY(LRB_OUT%K_STOC)
          ENDIF
        ENDIF
        CALL MUMPS_DM_FAC_UPD_DYN_MEMCNTS(-int(MEM,8),
     &        .TRUE., KEEP8, 
     &        IDUMMY, JDUMMY,    ! not used in case of deallocate
     &        .TRUE., .TRUE.)
      END SUBROUTINE DEALLOC_LRB
      SUBROUTINE DEALLOC_BLR_PANEL(BLR_PANEL, IEND, KEEP8, K34, IBEG_IN)
C       DEALLOCATE all B and R blocks allocated
        INTEGER, INTENT(IN)           :: IEND
        TYPE(LRB_TYPE), INTENT(INOUT) :: BLR_PANEL(:)
        INTEGER(8) :: KEEP8(150)
        INTEGER, INTENT(IN) :: K34
        INTEGER, INTENT(IN), OPTIONAL :: IBEG_IN
C     local variables
        INTEGER :: I, IBEG
        IF (present(IBEG_IN)) THEN
          IBEG = IBEG_IN
        ELSE
          IBEG = 1
        ENDIF
        IF (IEND.GE.IBEG) THEN
          IF (BLR_PANEL(1)%M.NE.0) THEN
            DO I=IBEG, IEND
              CALL DEALLOC_LRB(BLR_PANEL(I), KEEP8, K34)
            ENDDO
          ENDIF
        ENDIF
      END SUBROUTINE DEALLOC_BLR_PANEL
      END MODULE CMUMPS_LR_TYPE
