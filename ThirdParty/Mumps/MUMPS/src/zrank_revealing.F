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
      SUBROUTINE ZMUMPS_GET_NS_OPTIONS_FACTO(N,KEEP,ICNTL,MPG)
      IMPLICIT NONE
          INTEGER N, KEEP(500), ICNTL(60), MPG
          KEEP(19)=0
          KEEP(19)=ICNTL(56)
          IF ((KEEP(19).LT.1).OR.(KEEP(19).GE.2)) KEEP(19)=0
          IF ( KEEP(53) .LE. 0 .and.
     &         KEEP(19) .NE. 0 ) THEN
             KEEP(19) = 0
             IF ( MPG .GT. 0 ) THEN
      WRITE( MPG,'(A)') '** Warning:  ICNTL(56) null space option'
      WRITE( MPG,'(A)') '** disabled (incompatibility with analysis)'
             END IF
          END IF
          KEEP(21) = min(ICNTL(57),N)
          KEEP(22) = max(ICNTL(55),0)
          IF ( KEEP(19) .ne. 0 .and. KEEP(60) .ne. 0 ) THEN
             IF ( MPG .GT. 0 ) THEN
      WRITE( MPG,'(A)') '** Warning: ICNTL(56)  null space option'
      WRITE( MPG,'(A)') '** disabled (incompatibility with Schur)'
             END IF
            KEEP(19) = 0
          END IF
          RETURN
      END SUBROUTINE ZMUMPS_GET_NS_OPTIONS_FACTO
      SUBROUTINE ZMUMPS_GET_NS_OPTIONS_SOLVE(ICNTL, KEEP, 
     &                                       NRHS, MPG, INFO)
      IMPLICIT NONE
      INTEGER, intent(in)   :: KEEP(500), NRHS, MPG, ICNTL(60)
      INTEGER, intent(inout):: INFO(80)
      IF (KEEP(19).EQ.0.AND.KEEP(110).EQ.0) THEN
        IF (KEEP(111).NE.0) THEN
         INFO(1) = -37
         INFO(2) = 56  
         IF (KEEP(110).EQ.0) INFO(2) = 24 
          IF(MPG.GT.0) THEN
           WRITE( MPG,'(A)')
     &'** ERROR  : Null space computation requirement'
          WRITE( MPG,'(A)')
     &'** not consistent with factorization options'
         ENDIF
         GOTO 333
        ENDIF
      ENDIF
       IF (ICNTL(9).NE.1) THEN
         IF (KEEP(111).NE.0) THEN
          INFO(1) = -37
          INFO(2) = 9  
          IF (MPG.GT.0) THEN
           WRITE(MPG,'(A)')
     &'** ERROR  ICNTL(25) incompatible with '
           WRITE( MPG,'(A)')
     &'** option transposed system (ICNTL(9).ne.1) '
          ENDIF
         ENDIF
         GOTO 333
       ENDIF
      IF (KEEP(19).EQ.2) THEN
       IF ((KEEP(111).NE.0).AND.(KEEP(50).EQ.0)) THEN
         INFO(1) = -37
         INFO(2) = 0  
         IF (MPG.GT.0) THEN
          WRITE(MPG,'(A)')
     &'** ERROR  ICNTL(25) incompatible with '
          WRITE( MPG,'(A)')
     &'** option RRQR (ICNLT(56)=2) and unsym. matrices '
         ENDIF
        ENDIF
        GOTO 333
      ENDIF
      IF (KEEP(111).eq.-1.AND.NRHS.NE.KEEP(112)+KEEP(17))THEN
            INFO(1)=-32
            INFO(2)=NRHS
            GOTO 333
      ENDIF
      IF (KEEP(111).gt.0 .AND. NRHS .NE. 1) THEN
            INFO(1)=-32
            INFO(2)=NRHS
            GOTO 333
      ENDIF
      IF (KEEP(248) .NE.0.AND.KEEP(111).NE.0) THEN
         IF (MPG.GT.0) THEN
            WRITE(MPG,'(A)')
     &           ' ERROR: ICNTL(20) and ICNTL(30) functionalities ',
     &           ' incompatible with null space'
         ENDIF
         INFO(1) = -37
         IF (KEEP(237).NE.0) THEN
            INFO(2) = 30    
            IF (MPG.GT.0) THEN
               WRITE(MPG,'(A)')
     &           ' ERROR: ICNTL(30) functionality ',
     &              ' incompatible with null space'
            ENDIF
         ELSE
            IF (MPG.GT.0) THEN
               WRITE(MPG,'(A)')
     &              ' ERROR: ICNTL(20) functionality ',
     &              ' incompatible with null space'
            ENDIF
            INFO(2) = 20    
         ENDIF
         GOTO 333
      ENDIF
      IF (( KEEP(111) .LT. -1 ) .OR.
     &     (KEEP(111).GT.KEEP(112)+KEEP(17)) .OR.
     &     (KEEP(111) .EQ.-1 .AND. KEEP(112)+KEEP(17).EQ.0))
     &     THEN
         INFO(1)=-36
         INFO(2)=KEEP(111)
         GOTO 333
      ENDIF
      IF (KEEP(221).NE.0.AND.KEEP(111).NE.0) THEN
         INFO(1)=-37
         INFO(2)=26
         GOTO 333
      ENDIF
 333  CONTINUE
      RETURN
      END SUBROUTINE ZMUMPS_GET_NS_OPTIONS_SOLVE
      SUBROUTINE ZMUMPS_RR_INIT_POINTERS(roota)
      USE ZMUMPS_INTR_TYPES, ONLY: ZMUMPS_ROOT_STRUC
      IMPLICIT NONE
      TYPE (ZMUMPS_ROOT_STRUC) :: roota
      NULLIFY(roota%QR_TAU)
      NULLIFY(roota%SVD_U)
      NULLIFY(roota%SVD_VT)
      NULLIFY(roota%SINGULAR_VALUES)
      RETURN
      END SUBROUTINE ZMUMPS_RR_INIT_POINTERS
      SUBROUTINE ZMUMPS_RR_FREE_POINTERS(roota)
      USE ZMUMPS_INTR_TYPES, ONLY: ZMUMPS_ROOT_STRUC
      IMPLICIT NONE
      TYPE (ZMUMPS_ROOT_STRUC) :: roota
      IF (associated(roota%QR_TAU))  THEN
        DEALLOCATE(roota%QR_TAU)
        NULLIFY(roota%QR_TAU)
      ENDIF
      IF (associated(roota%SVD_U))  THEN
        DEALLOCATE(roota%SVD_U)
        NULLIFY(roota%SVD_U)
      ENDIF
      IF (associated(roota%SVD_VT))  THEN
        DEALLOCATE(roota%SVD_VT)
        NULLIFY(roota%SVD_VT)
      ENDIF
      IF (associated(roota%SINGULAR_VALUES))  THEN
        DEALLOCATE(roota%SINGULAR_VALUES)
        NULLIFY(roota%SINGULAR_VALUES)
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_RR_FREE_POINTERS
      SUBROUTINE ZMUMPS_SEQ_SYMMETRIZE(N,A)
      INTEGER N
      COMPLEX(kind=8) A( N, N )
      INTEGER I,J
      DO I = 2, N
        DO J = 1, I - 1
          A( I, J ) = A( J, I )
        END DO
      END DO
      RETURN
      END SUBROUTINE ZMUMPS_SEQ_SYMMETRIZE
      SUBROUTINE ZMUMPS_UXVSBP(N,PERM,X,RN01)
      INTEGER N,PERM(N),I
      COMPLEX(kind=8) RN01(N),X(N)
      DO I=1,N
        RN01(PERM(I))=X(I)
      ENDDO
      DO I=1,N
        X(I)=RN01(I)
      ENDDO
      RETURN
      END SUBROUTINE ZMUMPS_UXVSBP
      SUBROUTINE ZMUMPS_UXVSFP(N,PERM,X,RN01)
      INTEGER N,PERM(N),I
      COMPLEX(kind=8) RN01(N),X(N)
      DO I=1,N
        RN01(I)=X(PERM(I))
      ENDDO
      DO I=1,N
        X(I)=RN01(I)
      ENDDO
      RETURN
      END SUBROUTINE ZMUMPS_UXVSFP
      SUBROUTINE ZMUMPS_SVD_QR_ESTIM_WK( PHASE, MBLOCK, NBLOCK, 
     &           SIZE_ROOT_ARG,
     &           LOCAL_M, LOCAL_N, ROOT_OWNER, KEEP,KEEP8,
     &           LIWK_RR, LWK_RR )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: PHASE, SIZE_ROOT_ARG
      INTEGER, INTENT(IN) :: MBLOCK, NBLOCK, LOCAL_M, LOCAL_N
      LOGICAL, INTENT(IN) :: ROOT_OWNER
      INTEGER, INTENT(IN) :: KEEP(500)
      INTEGER(8), INTENT(IN) :: KEEP8(150)
      INTEGER, INTENT(OUT):: LIWK_RR
      INTEGER(8), INTENT(OUT):: LWK_RR
      INTEGER SIZE_ROOT
      INTEGER NBPOSPONED_ESTIM
      PARAMETER (NBPOSPONED_ESTIM=2000)
      INTEGER SVD_QR,PAR_ROOT
      SVD_QR = KEEP(19)
      PAR_ROOT = KEEP(38)
      LIWK_RR = 0
      LWK_RR = 0_8
      IF (PAR_ROOT.EQ.0) THEN
         IF(ROOT_OWNER) THEN
           IF (PHASE.EQ.0) THEN
                  SIZE_ROOT=SIZE_ROOT_ARG+NBPOSPONED_ESTIM
           ELSE
                  SIZE_ROOT=SIZE_ROOT_ARG
           ENDIF
           IF(SVD_QR.EQ.1) THEN
              LWK_RR=int(3*SIZE_ROOT+1,8)
           ELSEIF(SVD_QR.EQ.2) THEN
              LWK_RR=int(SIZE_ROOT+1,8)
           END IF
         END IF
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_SVD_QR_ESTIM_WK
      SUBROUTINE ZMUMPS_SEQ_FACTO_ROOT_SVD_QR
     &(NN,A,root,roota,WR03,LWR03,KEEP,KEEP8,INFO,LP,DKEEP,
     & GLOBK109,OPELIW,PIVNUL_LIST,LPIVNUL_LIST,
     & ROW_INDICES)
      USE MUMPS_INTR_TYPES, ONLY : MUMPS_ROOT_STRUC
      USE ZMUMPS_INTR_TYPES, ONLY : ZMUMPS_ROOT_STRUC
      IMPLICIT NONE
      TYPE ( MUMPS_ROOT_STRUC ) :: root
      TYPE ( ZMUMPS_ROOT_STRUC ) :: roota
      INTEGER :: NN,LP,LWR03,LWR03_MINSIZE
      COMPLEX(kind=8) :: A(NN*NN)
      INTEGER :: INFO(2),KEEP(500)
      INTEGER(8) :: KEEP8(150)
      DOUBLE PRECISION :: DKEEP(230)
      DOUBLE PRECISION :: OPELIW
      INTEGER :: GLOBK109
      INTEGER :: LPIVNUL_LIST
      INTEGER :: PIVNUL_LIST(LPIVNUL_LIST)
      INTEGER :: ROW_INDICES(NN)
      COMPLEX(kind=8) :: WR03(LWR03)
      INTEGER LDLT,DEFICIENCY
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RWORK
      INTEGER :: I,LDA,LDU,LDVT,J
      INTEGER :: IERR, LAST_BEFORE_GAP_IND
      INTEGER :: LAST_BEFORE_GAPLIMIT_IND, FIRST_AFTER_MinPiv, 
     &           FIRST_AFTER_GAPLIMIT, START_POINT, END_POINT 
      INTEGER :: ALLOCOK,MAXDEF,MINDEF
      DOUBLE PRECISION :: EPS, ZERO, GAPLIMIT, MaxGap, MaxGap1,
     &        MinPiv, Tol_MaxGap
      PARAMETER(ZERO=0.0D0)
      EPS = epsilon(ZERO)
      LDLT=KEEP(50)
      IF ((KEEP(19) .NE. 1).AND.(KEEP(19) .NE. 2)) THEN
        INFO(1)=-107
        INFO(2)= KEEP(19)
        IF ( LP .GT. 0 ) THEN
            WRITE(LP,*) " *** Option ",KEEP(19),
     &           " for null space no more available."
        ENDIF  
        GOTO 100
      ENDIF
      IF(KEEP(19).EQ.1) THEN
              LWR03_MINSIZE=3*NN+1
      ELSEIF(KEEP(19).EQ.2) THEN
              LWR03_MINSIZE=NN+1
      END IF
      MAXDEF=KEEP(21)
      IF ( MAXDEF .LE. 0 ) THEN
        MAXDEF = NN
      ELSE
        MAXDEF = max(MAXDEF - GLOBK109,0)
      ENDIF
      MINDEF = max(KEEP(22) - GLOBK109,0)
      MINDEF = min(MINDEF,NN)
      MAXDEF = min(MAXDEF,NN)
      IF(KEEP(19).EQ.1) THEN
         OPELIW = OPELIW + dble(26)*dble(NN)*dble(NN)*dble(NN)
      ELSEIF(KEEP(19).EQ.2) THEN
         OPELIW = OPELIW +  dble(4)*dble(NN)*dble(NN)*dble(NN)/dble(3)
      ENDIF
      IF (associated(roota%SINGULAR_VALUES)) 
     &     DEALLOCATE(roota%SINGULAR_VALUES)
      NULLIFY(roota%SINGULAR_VALUES)
      root%NB_SINGULAR_VALUES=NN
      ALLOCATE(roota%SINGULAR_VALUES(NN), stat=ALLOCOK)
      IF (ALLOCOK.GT.0) THEN
         INFO(1)=-13
         INFO(2)=NN
         IF ( LP .GT. 0 )
     &        WRITE(LP,*) 'ZMUMPS_SEQ_FACTO_ROOT_SVD_QR: Problem '
     &        //'allocating SINGULAR_VALUES'
         GOTO 100
      END IF
      IF(KEEP(19).EQ.1) THEN     
         IF(associated(roota%SVD_U)) DEALLOCATE(roota%SVD_U)
         NULLIFY(roota%SVD_U)
         ALLOCATE(roota%SVD_U(NN,NN),stat=ALLOCOK )
         IF(ALLOCOK.GT.0) THEN
            INFO(1)=-13
            INFO(2)=NN*NN
            IF ( LP .GT. 0 )
     &           WRITE(LP,*) 'ZMUMPS_SEQ_FACTO_ROOT_SVD_QR: Problem '
     &           //'allocating SVD_U'
            GOTO 100
         END IF
         IF (associated(roota%SVD_VT)) DEALLOCATE(roota%SVD_VT)
         NULLIFY(roota%SVD_VT)
         ALLOCATE(roota%SVD_VT(NN,NN), stat=ALLOCOK)
         IF (ALLOCOK.GT.0) THEN
            INFO(1)=-13
            INFO(2)=NN*NN
            IF ( LP .GT. 0 )
     &           WRITE(LP,*) 'ZMUMPS_SEQ_FACTO_ROOT_SVD_QR: Problem '
     &           //'allocating SVD_VT'
            GOTO 100
         END IF
         IF (allocated(RWORK)) DEALLOCATE(RWORK)
         ALLOCATE(RWORK(5*NN), stat=ALLOCOK)
         IF (ALLOCOK.GT.0) THEN
            INFO(1)=-13
            INFO(2)=5*NN
            IF ( LP .GT. 0 )
     &           WRITE(LP,*) 
     &           'ZMUMPS_SEQ_FACTO_ROOT_SVD_QR: Problem '
     &           //'allocating RWORK'
            GOTO 100
         END IF
      ELSEIF(KEEP(19).EQ.2) THEN
         IF (associated(roota%QR_TAU)) DEALLOCATE(roota%QR_TAU)
         NULLIFY(roota%QR_TAU)
         ALLOCATE(roota%QR_TAU(NN), stat=ALLOCOK)
         IF (ALLOCOK.GT.0) THEN
            INFO(1)=-13
            INFO(2)=NN
            IF ( LP .GT. 0 )
     &           WRITE(LP,*) 
     &           'ZMUMPS_SEQ_FACTO_ROOT_SVD_QR: Problem '
     &           //'allocating QR_TAU'
            GOTO 100
         END IF
         IF (associated(ROOT%IPIV)) DEALLOCATE(ROOT%IPIV)
         NULLIFY(ROOT%IPIV)
         ALLOCATE(ROOT%IPIV(NN), stat=ALLOCOK)
         IF (ALLOCOK.GT.0) THEN
            INFO(1)=-13
            INFO(2)=NN
            IF ( LP .GT. 0 )
     &           WRITE(LP,*) 
     &           'ZMUMPS_SEQ_FACTO_ROOT_SVD_QR: Problem '
     &           //'allocating IPIV'
            GOTO 100
         END IF
         IF (allocated(RWORK)) DEALLOCATE(RWORK)
         ALLOCATE(RWORK(2*NN), stat=ALLOCOK)
         IF (ALLOCOK.GT.0) THEN
            INFO(1)=-13
            INFO(2)=2*NN
            IF ( LP .GT. 0 )
     &           WRITE(LP,*) 
     &           'ZMUMPS_SEQ_FACTO_ROOT_SVD_QR: Problem '
     &           //'allocating RWORK'
            GOTO 100
         END IF
      ENDIF
      IF (LDLT.NE.0) THEN
         CALL ZMUMPS_SEQ_SYMMETRIZE(NN,A)
      END IF
      LDA=NN
      LDU=NN
      LDVT=NN
      IERR = 0
      IF(KEEP(19).EQ.1) THEN
         CALL zgesvd('A','A',NN,NN,A,LDA,roota%SINGULAR_VALUES(1)
     &        ,roota%SVD_U(1,1)
     &        ,LDU,roota%SVD_VT(1,1),LDVT,WR03,LWR03,RWORK,IERR)
      ENDIF
      IF(IERR.NE.0) THEN
         INFO(1)=-107
         INFO(2)=IERR
         IF (LP.GT.0) THEN
            IF(KEEP(19).EQ.1) THEN
               WRITE(LP,*) ' Problem in zgesvd : IERR = ', IERR
            ELSEIF(KEEP(19).EQ.2) THEN
               WRITE(LP,*) ' Problem in zgeqpf : IERR = ', IERR
            ENDIF
           GOTO 100
         END IF
      ENDIF
      IF(KEEP(19).EQ.2) THEN
         DO I=1,NN
            roota%SINGULAR_VALUES(I)=abs(A(I+NN*(I-1)))
         ENDDO
      ENDIF
      DEFICIENCY=0
      MinPiv = DKEEP(20)
      GAPLIMIT = DKEEP(9)
      IF (roota%SINGULAR_VALUES(NN).GT.MinPiv) THEN
              DEFICIENCY = 0
              GOTO 170
      ENDIF
      IF (roota%SINGULAR_VALUES(1).LE.GAPLIMIT) THEN
           DEFICIENCY = NN
           GOTO 170
      ENDIF
      LAST_BEFORE_GAPLIMIT_IND = 0
      LAST_BEFORE_GAP_IND  = 0
      FIRST_AFTER_MinPiv   = 0
      FIRST_AFTER_GAPLIMIT = 0
      MaxGap  = 0
      MaxGap1 = 0
      Tol_MaxGap = DKEEP(24)
      DO I=NN,1,-1
              IF (FIRST_AFTER_MinPiv.GT.0)  exit   
         IF(roota%SINGULAR_VALUES(I).LE.GAPLIMIT) THEN
                  LAST_BEFORE_GAPLIMIT_IND = I
         ELSE IF  ((FIRST_AFTER_GAPLIMIT.EQ.0).AND.
     &               (roota%SINGULAR_VALUES(I).LE.MinPiv)) THEN
                 FIRST_AFTER_GAPLIMIT = I
         ELSE IF (roota%SINGULAR_VALUES(I).GT.MinPiv) THEN
                 FIRST_AFTER_MinPiv = I
                 IF (FIRST_AFTER_GAPLIMIT.EQ.0) FIRST_AFTER_GAPLIMIT = I
         ENDIF
      ENDDO
      START_POINT = LAST_BEFORE_GAPLIMIT_IND
      IF ((LAST_BEFORE_GAPLIMIT_IND.EQ.0).AND.
     &    (FIRST_AFTER_GAPLIMIT.GT. FIRST_AFTER_MinPiv)) 
     &               START_POINT = FIRST_AFTER_GAPLIMIT
      END_POINT = FIRST_AFTER_MinPiv
      IF (FIRST_AFTER_MinPiv.EQ.0) END_POINT = 1
      DO I=START_POINT,END_POINT+1,-1
          IF (roota%SINGULAR_VALUES(I).EQ.0) THEN
             LAST_BEFORE_GAP_IND = I
          ELSE
              MaxGap1 = roota%SINGULAR_VALUES(I-1)*
     &                                    (1/roota%SINGULAR_VALUES(I))
             IF (MaxGap1.GE. Tol_MaxGap) THEN
                  IF (MaxGap1.GE. DKEEP(25)*MaxGap ) THEN 
                     LAST_BEFORE_GAP_IND = I
                     MaxGap =  MaxGap1
                  ENDIF  
             ENDIF  
          ENDIF
      ENDDO
      IF (MaxGap.EQ.ZERO) THEN        
          IF (LAST_BEFORE_GAPLIMIT_IND.EQ.0) THEN
               DEFICIENCY = 0
          ELSE
               DEFICIENCY = NN - LAST_BEFORE_GAPLIMIT_IND +1  
          ENDIF
      ELSE
          DEFICIENCY = NN - LAST_BEFORE_GAP_IND +1  
      ENDIF 
 170  CONTINUE
      DEFICIENCY=min(DEFICIENCY,MAXDEF)
      DEFICIENCY=max(DEFICIENCY,MINDEF)      
      KEEP(17)=DEFICIENCY
      IF(KEEP(19).EQ.2) THEN
         IF(DEFICIENCY.GT.0) THEN
            CALL ztrtrs('U','N','N',NN-DEFICIENCY,DEFICIENCY,
     &           A,LDA,A(LDA*(NN-DEFICIENCY)+1),LDA,IERR)
            IF ( IERR .NE. 0 ) THEN
               IF (LP.GT.0) 
     &         WRITE(LP,*) ' Internal error in ztrtrs: IERR = ',IERR
               CALL MUMPS_ABORT()
            END IF
         END IF
      ENDIF
      DO J=NN-DEFICIENCY+1, NN
         IF(KEEP(19).EQ.1) THEN
            PIVNUL_LIST(J-NN+DEFICIENCY) = ROW_INDICES(J)
         ELSEIF(KEEP(19).EQ.2) THEN
            PIVNUL_LIST(J-NN+DEFICIENCY) = ROW_INDICES(root%IPIV(J))
         ENDIF
      ENDDO
 100  CONTINUE     
      IF (allocated(RWORK)) DEALLOCATE(RWORK)
      RETURN
      END SUBROUTINE ZMUMPS_SEQ_FACTO_ROOT_SVD_QR
      SUBROUTINE ZMUMPS_SEQ_SOLVE_ROOT_SVD_QR
     &  (NRHS,NN,A,root, roota,
     &   IBEG_ROOT_DEF, IEND_ROOT_DEF,
     &   RHS,KEEP,KEEP8,MTYPE,INFO,LWK8,WK, LP)
      USE MUMPS_INTR_TYPES, ONLY : MUMPS_ROOT_STRUC
      USE ZMUMPS_INTR_TYPES, ONLY : ZMUMPS_ROOT_STRUC
      IMPLICIT NONE
      INTEGER NN, NRHS
      INTEGER(8), INTENT(IN) :: LWK8
      TYPE ( MUMPS_ROOT_STRUC ) :: root
      TYPE ( ZMUMPS_ROOT_STRUC ) :: roota
      COMPLEX(kind=8) A(NN*NN)
      INTEGER IBEG_ROOT_DEF, IEND_ROOT_DEF, MTYPE
      INTEGER INFO(2),KEEP(500)
      INTEGER(8) KEEP8(150)
      COMPLEX(kind=8) RHS(NN,NRHS), WK(LWK8)
      INTEGER LP
      INTEGER :: LWK
      COMPLEX(kind=8),DIMENSION(:,:), allocatable :: TEMP_RHS
      INTEGER :: I,IERR,K
      INTEGER :: LDLT,RRSTRAT,DEFICIENCY,LDA,LDRHS
      INTEGER :: ALLOCOK
      DOUBLE PRECISION, PARAMETER    :: RONE=1.0D+0
      COMPLEX(kind=8) ZERO, ONE, MINUSONE
      PARAMETER( ZERO = (0.0D0,0.0D0), ONE = (1.0D0,0.0D0))
      PARAMETER( MINUSONE=(-1.0D0,0.0D0))
      LDLT      = KEEP(50)
      RRSTRAT   = KEEP(19)
      DEFICIENCY = KEEP(17)
      LDA = NN
      LDRHS = NN
      LWK = int(min(int(huge(LWK),8),LWK8))
      IERR = 0
      IF ((RRSTRAT .NE. 1).AND.(RRSTRAT .NE. 2)) THEN
         WRITE(*,*) " *** Internal error ption ",RRSTRAT,
     &        " for null space no more available."
         CALL MUMPS_ABORT()
      ENDIF
      IF (KEEP(111).EQ.0) THEN
         IF(KEEP(19).EQ.1) THEN
           IF ((MTYPE.EQ.0).OR.(LDLT.GT.0)) THEN
               allocate(TEMP_RHS(NN,NRHS), stat=ALLOCOK)
               IF (ALLOCOK.GT.0) THEN
                 INFO(1)=-13
                 INFO(2)=NN*NRHS
                 IF ( LP .GT. 0 )
     &              WRITE(LP,*) 
     &              'ZMUMPS_SEQ_SOLVE_ROOT_SVD_QR: Problem '
     &              //'allocating TEMP_RHS'
                 RETURN
               END IF
               CALL zgemm('Conjugate transpose','N',NN,NRHS,NN,ONE,
     &              roota%SVD_U(1,1),NN,RHS,
     &              NN,ZERO,TEMP_RHS,NN)
               DO I=1,NN-DEFICIENCY
                 TEMP_RHS( I, 1:NRHS ) = 
     &         cmplx(RONE/roota%SINGULAR_VALUES(I),kind=kind(TEMP_RHS))*
     &                          TEMP_RHS( I, 1:NRHS )
               ENDDO
               DO I=NN-DEFICIENCY +1, NN
                  TEMP_RHS(I, 1:NRHS) = ZERO 
               ENDDO
               CALL zgemm('Conjugate transpose','N',NN,NRHS,NN,ONE,
     &               roota%SVD_VT(1,1),NN,
     &              TEMP_RHS, NN,ZERO,RHS,NN)
               DEALLOCATE(TEMP_RHS)
           ELSEIF(MTYPE.EQ.1) THEN
               allocate(TEMP_RHS(NRHS,NN), stat=ALLOCOK)
               IF (ALLOCOK.GT.0) THEN
                  INFO(1)=-13
                  INFO(2)=NN*NRHS
                  IF ( LP .GT. 0 )
     &              WRITE(LP,*) 
     &              'ZMUMPS_SEQ_SOLVE_ROOT_SVD_QR: Problem '
     &              //'allocating TEMP_RHS'
                  RETURN
               END IF
               CALL zgemm('T','C',NRHS,NN, NN,ONE,RHS, NN, 
     &               roota%SVD_VT(1,1),NN,
     &              ZERO,TEMP_RHS,NRHS)
               DO I=1,NN-DEFICIENCY
                  RHS(I, 1:NRHS ) = 
     &         cmplx(RONE/roota%SINGULAR_VALUES(I),kind=kind(TEMP_RHS))*
     &                          TEMP_RHS(1:NRHS,I ) 
               ENDDO
               DO I=NN-DEFICIENCY +1, NN
                  RHS(I,1:NRHS) = ZERO 
               ENDDO
               CALL zgemm('T','C',NRHS,NN,NN,ONE,RHS, NN, 
     &             roota%SVD_U(1,1),NN,
     &             ZERO,TEMP_RHS,NRHS)
                  DO I=1,NRHS
                  RHS(1:NN,I) =TEMP_RHS(I,1:NN )
                  ENDDO
                 DEALLOCATE(TEMP_RHS)
            ENDIF
         ELSEIF(KEEP(19).EQ.2) THEN
            IF ((MTYPE.EQ.0).OR.(LDLT.GT.0)) THEN
               CALL zunmqr('L','Conjugate transpose',NN,NRHS,NN,
     &              A,LDA,roota%QR_TAU(1),
     &              RHS,LDRHS,WK(1),LWK,IERR)
               IF(IERR.LT.0) THEN
                  WRITE(*,*) 
     &         'Error return from zunmqr in root solve: IERR=', IERR
                  RETURN
               END IF
               CALL ztrtrs('U','N','N',NN-DEFICIENCY,NRHS,A,LDA,
     &              RHS,LDRHS,IERR)
               IF ( IERR .LT. 0 ) THEN
                  WRITE(*,*) 
     &         'Error return from ztrtrs in roor solve: IERR =',IERR
                  RETURN
               END IF
               DO I=1,NRHS
                  RHS( NN - DEFICIENCY + 1: NN, I ) = ZERO  
               ENDDO
               DO I=1,NRHS
                  CALL ZMUMPS_UXVSBP(NN,root%IPIV(1),RHS(1,I),WK)
               ENDDO
           ELSEIF(MTYPE.EQ.1) THEN
               DO I=1,NRHS
                 CALL ZMUMPS_UXVSFP(NN,root%IPIV(1),RHS(1,I),WK)
               ENDDO
               CALL ztrtrs('U','T','N',NN-DEFICIENCY,NRHS,
     &              A,LDA,RHS,LDRHS,IERR)
               IF(IERR.NE.0) THEN
                  WRITE(*,*) 'Error return from trtrs: IERR=', IERR
                  STOP    
               END IF
               allocate(TEMP_RHS(NRHS,NN), stat=ALLOCOK)
               IF (ALLOCOK.GT.0) THEN
                  INFO(1)=-13
                  INFO(2)=NN*NRHS
                  IF ( LP .GT. 0 )
     &              WRITE(LP,*) 
     &              'ZMUMPS_SEQ_SOLVE_ROOT_SVD_QR: Problem '
     &              //'allocating TEMP_RHS'
                  RETURN
               END IF
               DO I=1,NRHS
                  TEMP_RHS(I,1:NN-DEFICIENCY)=RHS(1:NN-DEFICIENCY, I)
               ENDDO
               DO I=NN - DEFICIENCY + 1,NN
                  TEMP_RHS( 1: NRHS, I ) = ZERO  
               ENDDO
                CALL zunmqr( 'R','Conjugate transpose',NRHS,NN,NN,A,LDA,
     &              roota%QR_TAU(1),
     &              TEMP_RHS,NRHS,WK,LWK,IERR)
               IF(IERR.LT.0) THEN
                  WRITE(*,*) 'Error return from zunmqr: IERR=', IERR
                  RETURN
               END IF
                DO I=1,NRHS
                  RHS(1:NN, I)= TEMP_RHS(I,1:NN)
               ENDDO
               DEALLOCATE(TEMP_RHS)
            ENDIF
         ENDIF
      ELSE
         IF(KEEP(19).EQ.1) THEN
           IF ((MTYPE.EQ.0).OR.(LDLT.GT.0)) THEN
               DO I=IBEG_ROOT_DEF,IEND_ROOT_DEF
                  RHS(:,I+1-IBEG_ROOT_DEF) = 
     &                 CONJG(roota%SVD_VT(NN-DEFICIENCY+I,:))
               ENDDO
           ELSEIF(MTYPE.EQ.1) THEN
               DO I=IBEG_ROOT_DEF,IEND_ROOT_DEF
                  RHS(:,I+1-IBEG_ROOT_DEF) = 
     &               CONJG(roota%SVD_U(:,NN-DEFICIENCY+I))
               ENDDO
            ENDIF
         ELSEIF(KEEP(19).EQ.2) THEN
           IF((MTYPE.EQ.0).OR.(LDLT.GT.0)) THEN        
               DO I=IBEG_ROOT_DEF,IEND_ROOT_DEF
                  RHS(NN-DEFICIENCY+I,I-IBEG_ROOT_DEF+1) = MINUSONE
                  DO K=1,NN-DEFICIENCY
                     RHS(K,I-IBEG_ROOT_DEF+1)=
     &                    A(K + LDA*(NN-DEFICIENCY+I-1))
                  ENDDO
               ENDDO
               DO I=1,IEND_ROOT_DEF-IBEG_ROOT_DEF+1
                  CALL ZMUMPS_UXVSBP(NN,root%IPIV(1),RHS(1,I),WK)
               ENDDO 
           ELSEIF(MTYPE.EQ.1) THEN
               WRITE(*,*) 'Computation of a null space basis'
     &              // ' of A is unavailable for unsymetric matrices'
            DO I=IBEG_ROOT_DEF,IEND_ROOT_DEF
               RHS(NN-DEFICIENCY+I,I-IBEG_ROOT_DEF+1) = ONE
            ENDDO
            CALL zunmqr('L','N',NN,NRHS,NN, A,LDA,roota%QR_TAU(1),
     &              RHS,LDRHS,WK(1),LWK,IERR)
            ENDIF
         ENDIf
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_SEQ_SOLVE_ROOT_SVD_QR
