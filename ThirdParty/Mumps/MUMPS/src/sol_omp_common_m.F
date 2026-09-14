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
C***********************************************************************
      MODULE MUMPS_SOL_L0OMP_M
C
C     Purpose:
C     =======
C     Manage locks for right-looking updates of RHSINTR unde L0 threads
C
!$    USE OMP_LIB, ONLY: OMP_LOCK_KIND
      INTEGER, PARAMETER :: NB_LOCK_MAX = 18
!$    INTEGER(OMP_LOCK_KIND),
!$   &ALLOCATABLE, DIMENSION(:), SAVE :: LOCK_FOR_SCATTER
C
C
      CONTAINS
C
      SUBROUTINE MUMPS_SOL_L0OMP_LI( K400 )
!$    USE OMP_LIB, ONLY: OMP_INIT_LOCK
      IMPLICIT NONE
C
C     Purpose:
C     =======
C     Initialize locks for forward solution with L0-threads feature
C     (LI suffix: Lock Initialization)
C
C     Argument:
C     ========
C     K400: the number of threads for L0-threads;
C           we use min(K400, NB_LOCK_MAX) locks.
C
      INTEGER, INTENT(IN) :: K400
C
!$    INTEGER :: I
C
C     Executable statements
C     =====================
C
!$    IF (K400 .GT. 0) THEN
!$      ALLOCATE(LOCK_FOR_SCATTER(min(NB_LOCK_MAX,K400)))
!$      DO I = 1, min(NB_LOCK_MAX,K400)
!$        CALL OMP_INIT_LOCK(LOCK_FOR_SCATTER(I))
!$      ENDDO
!$    ENDIF
      RETURN
      END SUBROUTINE MUMPS_SOL_L0OMP_LI
      SUBROUTINE MUMPS_SOL_L0OMP_LD( K400 )
!$    USE OMP_LIB, ONLY : OMP_DESTROY_LOCK
      IMPLICIT NONE
C
C     Purpose:
C     =======
C     Destroy locks for forward solution with L0-threads feature
C     (LD suffix: Lock Destruction)
C
C     Argument:
C     ========
C     K400: the number of threads for L0-threads;
C           we use min(K400, NB_LOCK_MAX) locks.
      INTEGER, INTENT(IN) :: K400
C
!$    INTEGER :: I
C
C     Executable statements
C     =====================
C
!$    IF (allocated(LOCK_FOR_SCATTER)) THEN
!$      IF (K400 .GT. 0) THEN
!$        DO I = 1, min(NB_LOCK_MAX,K400)
!$          CALL OMP_DESTROY_LOCK(LOCK_FOR_SCATTER(I))
!$        ENDDO
!$        DEALLOCATE(LOCK_FOR_SCATTER)
!$      ENDIF
!$    ENDIF
      RETURN
      END SUBROUTINE MUMPS_SOL_L0OMP_LD
      END MODULE MUMPS_SOL_L0OMP_M
