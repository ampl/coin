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
      MODULE ZMUMPS_FAC_COMPACT_FACTORS_M
      PRIVATE
      PUBLIC :: ZMUMPS_TRY_COMPACT_FACTORS
      CONTAINS
      SUBROUTINE ZMUMPS_TRY_COMPACT_FACTORS(ICNTL49_LOC,
     &      WK_USER_PROVIDED, S, KEEP, KEEP8, INFO, MYID, ICNTL,
     &      PROK, MP, ZMUMPS_LBUFR_BYTES8, ZMUMPS_LBUF8, 
     &      LIWK, LIWK8 )
!$    USE OMP_LIB
      USE ZMUMPS_DYNAMIC_MEMORY_M, ONLY : ZMUMPS_DM_FREE_S_WK
C
C Purpose
C =======
C    If no factors stored in S and .NOT.WK_USER_PROVIDED deallocate(S)
C    If ICNTL49_LOC = 1, 2 try to compress S
C       Possible values :
C           0 : nothing is done.
C           1 : compact  S while satisfying the
C               memory constraint that might have been provided 
C               with ICNTL(23) feature.
C           2 : compact S. The memory constraint that might have been 
C               provided with ICNTL(23) feature does not apply 
C
C Parameters
C ==========
      INTEGER :: ICNTL49_LOC, MP, MYID
      COMPLEX(kind=8), POINTER, DIMENSION(:) :: S
      INTEGER :: KEEP(500), INFO(80), ICNTL(60)
      LOGICAL :: PROK, WK_USER_PROVIDED
      INTEGER(8) :: ZMUMPS_LBUFR_BYTES8, ZMUMPS_LBUF8
      INTEGER(8) :: KEEP8(150)
      INTEGER(8), INTENT(IN) ::  LIWK, LIWK8
C
C     Local declarations
C     ==================
C
      LOGICAL    :: Compact_S_Authorized
      INTEGER    :: IERR, NOMP
      COMPLEX(kind=8), DIMENSION(:), POINTER :: TMPS
      INTEGER(8) :: TMPpeak, I8
!$    INTEGER(8) :: CHUNK8
      IF (.NOT.WK_USER_PROVIDED) THEN
C{
        IF (KEEP8(31).EQ.0) THEN
C{
C            No factors stored in S
          IF (associated(S)) THEN
            CALL ZMUMPS_DM_FREE_S_WK(S, KEEP(430))
C           Reset KEEP(430)=0 since next allocations of S
C           will be from Fortran
            KEEP(430)=0
            NULLIFY(S)
            KEEP8(23) = 0
          ENDIF
C}
        ELSE IF (ICNTL49_LOC.NE.0) THEN
C{        Factors stored in S, try to compact S
          TMPpeak = KEEP8(73) + KEEP8(31) 
     &     - (ZMUMPS_LBUFR_BYTES8+ZMUMPS_LBUF8)/int(KEEP(35),8)
     &     - KEEP8(26) 
     &     - ((LIWK+LIWK8*KEEP(10)+KEEP8(27))*int(KEEP(34),8))
     &             /int(KEEP(35),8)
          Compact_S_Authorized = .FALSE.
C         Set Compact_S_Authorized 
          IF (KEEP8(4).GT.0_8) THEN
            IF (TMPpeak.LT.KEEP8(75)) 
     &         Compact_S_Authorized=.TRUE.
          ELSE
               Compact_S_Authorized = .TRUE.
          ENDIF
          IF (ICNTL49_LOC.EQ.1.AND..NOT.Compact_S_Authorized) THEN
C{
              INFO(1) = INFO(1) + 4
C               INFO(2) = 
C           New value of ICNTL(23) (in MBytes:
C              ( KEEP8(4) +  (TMPpeak- KEEP8(75))*KEEP(35) )/1000000
C               + 1 for safety
                INFO(2) = int(
     &             ( 
     &               KEEP8(4) +  
     &              (TMPpeak- KEEP8(75))*int(KEEP(35),8) 
     &             ) / 1000000_8  + 1_8
     &                       )
C             In fact increasing INFO(2) will not help 
C             since increasing ICNTL(23) will also increase 
C             MAXS and thus the peak of memory. 
C             Thus setting ICNTL(23) to INFO(2) might not 
C             enable user to Compact_S.
C             Simplest is to advice to set ICNTL(49)=2 
C             or to switch of ICNTL(23) feature.
              IF (PROK) THEN
                WRITE(MP,'(A,I4,A,I2,A,/A,/A,A)') 
     &              " ** WARNING ** on MPI proc= ", MYID, 
     &              " ICNTL(49)= ", ICNTL49_LOC, 
     &              ", but not enough memory to compact S due to ",
     &              "    memory limitation given by ICNTL(23).",
     &              "    ICNTL(23) should be reset to zero or",
     &              " ICNTL(49) should be set to 2 "
              ENDIF
C}
          ELSE IF (
     &        (ICNTL49_LOC.EQ.1.AND.Compact_S_Authorized) 
     &            .OR.
     &        (ICNTL49_LOC.EQ.2) 
C{
     &           ) THEN
C             Try to compact S of size MAXS
              ALLOCATE(TMPS(KEEP8(31)), stat=IERR)
              IF (IERR .GT. 0 ) THEN
                IF (PROK) THEN
                  WRITE(MP,'(A,I4,A,I3,A)') 
     &              " ** WARNING ** on MPI proc= ", MYID, 
     &              " ICNTL(49)= ", ICNTL49_LOC, 
     &              ", but not enough memory to compact S "
                ENDIF
                INFO(1) = INFO(1) + 4
                GOTO 513
              ENDIF
C
!$            NOMP = OMP_GET_MAX_THREADS()
!$            CHUNK8=int(KEEP(361),8)
!$OMP         PARALLEL DO PRIVATE(I8) SCHEDULE(STATIC, CHUNK8)
!$OMP&        IF ( KEEP8(31) > int(KEEP(361),8) .AND. NOMP .GT. 1)
#if defined(__ve__)
!NEC$ IVDEP
#endif
              DO I8=1_8, KEEP8(31)
                TMPS(I8) = S(I8)
              ENDDO
!$OMP         END PARALLEL DO
              CALL ZMUMPS_DM_FREE_S_WK(S, KEEP(430))
C             Reset KEEP(430)=0 since TMPS is allocated 
C             in Fortran and S=>TMPS should be deallocated
C             in Fortran.
              KEEP(430)=0
              S => TMPS;  NULLIFY(TMPS)
              KEEP8(23) = KEEP8(31)
C}
          ENDIF
C}
        ENDIF
C}
      ENDIF
 513  CONTINUE
      RETURN
      END SUBROUTINE ZMUMPS_TRY_COMPACT_FACTORS
      END MODULE ZMUMPS_FAC_COMPACT_FACTORS_M
