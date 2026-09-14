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
      MODULE MUMPS_TOOLS_COMMON_M
C
C     This module with arithmetic-independant utility routines
C     such as integer sorting algorithms
C
      PRIVATE
      PUBLIC :: MUMPS_MERGESORT, MUMPS_MERGESWAP, MUMPS_MERGESWAP1
#if ! defined(NOANA_BLKAUTO)
      PUBLIC :: MUMPS_MERGESORTK8, MUMPS_MERGESWAP48
#endif
      PUBLIC :: MUMPS_LR_COMPRESS_CB
      PUBLIC :: MUMPS_LR_COMPRESS_PANEL
      PUBLIC :: MUMPS_LR_FORCE_BFR
      CONTAINS
      SUBROUTINE MUMPS_MERGESORT(N, K, L)
      IMPLICIT NONE
C     Plain implementation of the merge-sort algorithm
C     as described in:
C
C     D. E. Knuth "The Art of Computer Programming,"
C     vol.3: Sorting and Searching, Addison-Wesley, 1973
      INTEGER    :: N
      INTEGER    :: K(:), L(0:)
      INTEGER    :: P, Q, S, T
      L(0) = 1
      T = N + 1
      DO  P = 1,N - 1
         IF (K(P) <= K(P+1)) THEN
            L(P) = P + 1
         ELSE
            L(T) = - (P+1)
            T = P
       END IF
      END DO
      L(T) = 0
      L(N) = 0
      IF (L(N+1) == 0) THEN
         RETURN 
      ELSE
         L(N+1) = iabs(L(N+1))
      END IF
 200  CONTINUE
      S = 0
      T = N+1
      P = L(S)
      Q = L(T)
      IF(Q .EQ. 0) RETURN
 300  CONTINUE
      IF(K(P) .GT. K(Q)) GOTO 600 
      L(S) = sign(P,L(S))
      S = P
      P = L(P)
      IF (P .GT. 0) GOTO 300
      L(S) = Q
      S = T
      DO
         T = Q
         Q = L(Q)
         IF (Q .LE. 0) EXIT
      END DO
      GOTO 800
 600  CONTINUE
      L(S) = sign(Q, L(S))
      S = Q
      Q = L(Q)
      IF (Q .GT. 0) GOTO 300
      L(S) = P
      S = T
      DO
         T = P
         P = L(P)
         IF (P .LE. 0) EXIT
      END DO
 800  CONTINUE
      P = -P
      Q = -Q
      IF(Q.EQ.0) THEN
         L(S) = sign(P, L(S))
         L(T) = 0
         GOTO 200
      END IF
      GOTO 300
      END SUBROUTINE MUMPS_MERGESORT
#if ! defined(NOANA_BLKAUTO)
C     =============================================
C     MUMPS_MERGESORTK8: same as MUMPS_MERGESORT above, except that
C     K is declared INTEGE(8) instead of INTEGER.
      SUBROUTINE MUMPS_MERGESORTK8(N, K, L)
      IMPLICIT NONE
      INTEGER    :: N
      INTEGER    :: L(0:)
      INTEGER(8) :: K(:)
      INTEGER    :: P, Q, S, T
      L(0) = 1
      T = N + 1
      DO  P = 1,N - 1
         IF (K(P) <= K(P+1)) THEN
            L(P) = P + 1
         ELSE
            L(T) = - (P+1)
            T = P
       END IF
      END DO
      L(T) = 0
      L(N) = 0
      IF (L(N+1) == 0) THEN
         RETURN 
      ELSE
         L(N+1) = iabs(L(N+1))
      END IF
 200  CONTINUE
      S = 0
      T = N+1
      P = L(S)
      Q = L(T)
      IF(Q .EQ. 0) RETURN
 300  CONTINUE
      IF(K(P) .GT. K(Q)) GOTO 600 
      L(S) = sign(P,L(S))
      S = P
      P = L(P)
      IF (P .GT. 0) GOTO 300
      L(S) = Q
      S = T
      DO
         T = Q
         Q = L(Q)
         IF (Q .LE. 0) EXIT
      END DO
      GOTO 800
 600  CONTINUE
      L(S) = sign(Q, L(S))
      S = Q
      Q = L(Q)
      IF (Q .GT. 0) GOTO 300
      L(S) = P
      S = T
      DO
         T = P
         P = L(P)
         IF (P .LE. 0) EXIT
      END DO
 800  CONTINUE
      P = -P
      Q = -Q
      IF(Q.EQ.0) THEN
         L(S) = sign(P, L(S))
         L(T) = 0
         GOTO 200
      END IF
      GOTO 300
      END SUBROUTINE MUMPS_MERGESORTK8
#endif
      SUBROUTINE MUMPS_MERGESWAP1(N, L, A)
      IMPLICIT NONE
      INTEGER   :: I, LP, ISWAP, N
      INTEGER   :: L(0:), A(:)
      LP = L(0)
      I  = 1
      DO 
         IF ((LP==0).OR.(I>N)) EXIT
         DO 
            IF (LP >= I) EXIT
            LP = L(LP)
         END DO
         ISWAP    = A(LP)
         A(LP)   = A(I)
         A(I)    = ISWAP
         ISWAP    = L(LP)
         L(LP) = L(I)
         L(I)  = LP
         LP = ISWAP 
         I  = I + 1
      ENDDO
      END SUBROUTINE MUMPS_MERGESWAP1
      SUBROUTINE MUMPS_MERGESWAP(N, L, A1, A2)
      IMPLICIT NONE
      INTEGER   :: I, LP, ISWAP, N
      INTEGER   :: L(0:), A1(:), A2(:)
      LP = L(0)
      I  = 1
      DO 
         IF ((LP==0).OR.(I>N)) EXIT
         DO 
            IF (LP >= I) EXIT
            LP = L(LP)
         END DO
         ISWAP    = A1(LP)
         A1(LP)   = A1(I)
         A1(I)    = ISWAP
         ISWAP    = A2(LP)
         A2(LP)   = A2(I)
         A2(I)    = ISWAP
         ISWAP    = L(LP)
         L(LP) = L(I)
         L(I)  = LP
         LP = ISWAP 
         I  = I + 1
      ENDDO
      END SUBROUTINE MUMPS_MERGESWAP
#if ! defined(NOANA_BLKAUTO)
C
C     MUMPS_MERGESWAP48: same as MUMPS_MERGESWAP except that
C     A2 is declared INTEGER(8) instead of INTEGER
      SUBROUTINE MUMPS_MERGESWAP48(N, L, A1, A2)
      IMPLICIT NONE
      INTEGER   :: I, LP, ISWAP, N
      INTEGER   :: L(0:), A1(:)
      INTEGER(8):: A2(:), ISWAP2
      LP = L(0)
      I  = 1
      DO 
         IF ((LP==0).OR.(I>N)) EXIT
         DO 
            IF (LP >= I) EXIT
            LP = L(LP)
         END DO
         ISWAP    = A1(LP)
         A1(LP)   = A1(I)
         A1(I)    = ISWAP
         ISWAP2   = A2(LP)
         A2(LP)   = A2(I)
         A2(I)    = ISWAP2
         ISWAP    = L(LP)
         L(LP) = L(I)
         L(I)  = LP
         LP = ISWAP 
         I  = I + 1
      ENDDO
      END SUBROUTINE MUMPS_MERGESWAP48
#endif
C
      FUNCTION MUMPS_LR_COMPRESS_PANEL(LRSTATUS) RESULT(COMPRESS_PANEL)
        IMPLICIT NONE
        INTEGER, intent(in)  :: LRSTATUS
        LOGICAL :: COMPRESS_PANEL
        COMPRESS_PANEL = LRSTATUS.GE.2.AND.LRSTATUS.LT.4
      END FUNCTION MUMPS_LR_COMPRESS_PANEL
      FUNCTION MUMPS_LR_COMPRESS_CB(LRSTATUS) RESULT(COMPRESS_CB)
        IMPLICIT NONE
        INTEGER, intent(in)  :: LRSTATUS
        LOGICAL :: COMPRESS_CB
        COMPRESS_CB = LRSTATUS.EQ.1.OR.LRSTATUS.EQ.3.OR.LRSTATUS.EQ.5
      END FUNCTION MUMPS_LR_COMPRESS_CB
      FUNCTION MUMPS_LR_FORCE_BFR(LRSTATUS) RESULT(FORCE_BFR)
        IMPLICIT NONE
        INTEGER, intent(in)  :: LRSTATUS
        LOGICAL :: FORCE_BFR
        FORCE_BFR = LRSTATUS.GE.4
      END FUNCTION MUMPS_LR_FORCE_BFR
      END MODULE MUMPS_TOOLS_COMMON_M
