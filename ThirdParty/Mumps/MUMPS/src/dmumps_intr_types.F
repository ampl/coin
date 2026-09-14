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
      MODULE DMUMPS_INTR_TYPES
      USE MUMPS_INTR_TYPES, ONLY: MUMPS_ROOT_STRUC
#if defined(SINGLEINDOUBLE)
      USE SMUMPS_INTR_TYPES, ONLY : SMUMPS_L0OMPFAC_T,
     &                               SMUMPS_ROOT_STRUC
#endif
      PRIVATE
      PUBLIC :: DMUMPS_ROOT_STRUC,
     &          DMUMPS_L0OMPFAC_T,
     &          DMUMPS_INTR_STRUC,
     &          DMUMPS_INIT_INTR_ENCODING,
     &          DMUMPS_FREE_INTR_ENCODING,
     &          DMUMPS_ENCODE_INTR,
     &          DMUMPS_DECODE_INTR
#if defined(SINGLEINDOUBLE)
     &       , SMUMPS_ROOT_STRUC, SMUMPS_L0OMPFAC_T
#endif
C     DMUMPS_ROOT_STRUC no longer contains INTEGERS
      TYPE DMUMPS_ROOT_STRUC
!       Centralized master of root
        DOUBLE PRECISION, DIMENSION(:), POINTER :: RHS_CNTR_MASTER_ROOT
!       Used to access Schur easily from root structure
        DOUBLE PRECISION, DIMENSION(:), POINTER :: SCHUR_POINTER
!       for null_space preprocessing only:
        DOUBLE PRECISION, DIMENSION(:), POINTER :: QR_TAU
!       Fwd in facto: 
!           case of scalapack root: to store RHS in 2D block cyclic
!           format compatible with root distribution
        DOUBLE PRECISION, DIMENSION(:,:), POINTER :: RHS_ROOT
!       for SVD on root (null_space)
        DOUBLE PRECISION, DIMENSION(:,:), POINTER :: SVD_U, SVD_VT
!       for RR on root (null_space)
        DOUBLE PRECISION, DIMENSION(:), POINTER :: SINGULAR_VALUES
!
      END TYPE DMUMPS_ROOT_STRUC
!     multicore
      TYPE DMUMPS_L0OMPFAC_T
         DOUBLE PRECISION, POINTER, DIMENSION(:) :: A
         INTEGER(8) :: LA
      END TYPE DMUMPS_L0OMPFAC_T
C
C     All MUMPS internal datatypes are in an internal structure:
      TYPE DMUMPS_INTR_STRUC
         TYPE (MUMPS_ROOT_STRUC)  :: root
         TYPE (DMUMPS_ROOT_STRUC) :: roota
         TYPE (DMUMPS_L0OMPFAC_T), DIMENSION(:), POINTER ::
     &                             L0_OMP_FACTORS
#if defined(SINGLEINDOUBLE)
         TYPE (SMUMPS_L0OMPFAC_T), DIMENSION(:), POINTER ::
     &                             LPL0_OMP_FACTORS
#endif
      END TYPE DMUMPS_INTR_STRUC
C     =================================================================
      CONTAINS
C     =================================================================
      SUBROUTINE DMUMPS_INIT_INTR_ENCODING(id_intr_ENCODING)
      IMPLICIT NONE
      CHARACTER, DIMENSION(:), POINTER :: id_intr_ENCODING
C     To be called only before JOB=-1
      NULLIFY(id_intr_ENCODING)
      END SUBROUTINE DMUMPS_INIT_INTR_ENCODING
C     =================================================================
      SUBROUTINE DMUMPS_FREE_INTR_ENCODING(id_intr_ENCODING)
      IMPLICIT NONE
      CHARACTER, DIMENSION(:), POINTER :: id_intr_ENCODING
C     To be called only after JOB=-2
      DEALLOCATE(id_intr_ENCODING)
      NULLIFY(id_intr_ENCODING)
      RETURN
      END SUBROUTINE DMUMPS_FREE_INTR_ENCODING
C     =================================================================
      SUBROUTINE DMUMPS_ENCODE_INTR(id_intr_ENCODING, id_intr)
      IMPLICIT NONE
C
C     Arguments:
C     =========
      CHARACTER, DIMENSION(:), POINTER :: id_intr_ENCODING
      TYPE (DMUMPS_INTR_STRUC)         :: id_intr
C
C     Local variables:
C     ===============
      CHARACTER :: CHAR_ARRAY(1)
      INTEGER :: CHAR_LENGTH, IERR
C
      IF (associated(id_intr_ENCODING)) THEN
C       Should be unassociated on entry
        WRITE(*,*) "Internal error in DMUMPS_ENCODE_INTR:",
     &             " id_intr_ENCODING already allocated"
        CALL MUMPS_ABORT()
      ENDIF
      CHAR_LENGTH=size(transfer(id_intr,CHAR_ARRAY))
      ALLOCATE(id_intr_ENCODING(CHAR_LENGTH), stat=IERR)
      IF (IERR > 0 ) THEN
        WRITE(*,*) "Allocation error in MUMPS_ENCODE_INTR"
        CALL MUMPS_ABORT()
      ENDIF
C     Fill with derived datatype
      id_intr_ENCODING=transfer(id_intr,CHAR_ARRAY)
      RETURN
      END SUBROUTINE DMUMPS_ENCODE_INTR
C     =================================================================
      SUBROUTINE DMUMPS_DECODE_INTR(id_intr_ENCODING, id_intr)
      IMPLICIT NONE
      CHARACTER(len=1), DIMENSION(:), POINTER :: id_intr_ENCODING
      TYPE (DMUMPS_INTR_STRUC) :: id_intr
      IF (.NOT.associated(id_intr_ENCODING)) THEN
        WRITE(*,*) "Internal error 1 in DMUMPS_DECODE_INTR"
        CALL MUMPS_ABORT()
      ENDIf
      id_intr=transfer(id_intr_ENCODING,id_intr)
      DEALLOCATE(id_intr_ENCODING)
      NULLIFY(id_intr_ENCODING)
      RETURN
      END SUBROUTINE DMUMPS_DECODE_INTR
      END MODULE DMUMPS_INTR_TYPES
