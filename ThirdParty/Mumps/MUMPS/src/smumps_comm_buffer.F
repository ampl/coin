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
      MODULE SMUMPS_BUF
      USE MUMPS_BUF_COMMON, ONLY: BUF_CB, SIZE_RBUF_BYTES,
     &      SIZEofINT, SIZEofREAL, OVHSIZE, BUF_ADJUST, BUF_LOOK,
     &      MUMPS_BUF_SIZE_AVAILABLE
      PRIVATE
      INTEGER, SAVE ::  BUF_LMAX_ARRAY
      REAL, DIMENSION(:), ALLOCATABLE
     &       , SAVE, TARGET :: BUF_MAX_ARRAY
      PUBLIC :: BUF_LMAX_ARRAY, BUF_MAX_ARRAY
      PUBLIC :: SMUMPS_BUF_DEALL_MAX_ARRAY,
     &          SMUMPS_BUF_MAX_ARRAY_MINSIZE
      PUBLIC :: SMUMPS_BUF_SEND_CB,
     &        SMUMPS_BUF_SEND_MASTER2SLAVE,
     &        SMUMPS_BUF_SEND_VCB,
     &        SMUMPS_BUF_SEND_MAITRE2,
     &        SMUMPS_BUF_SEND_CONTRIB_TYPE2,
     &        SMUMPS_BUF_SEND_BLOCFACTO,
     &        SMUMPS_BUF_SEND_BLFAC_SLAVE,
     &        SMUMPS_BUF_SEND_CONTRIB_TYPE3,
     &        SMUMPS_BUF_SEND_BACKVEC, 
     &        SMUMPS_MPI_UNPACK_LRB
      CONTAINS
        SUBROUTINE SMUMPS_BUF_DEALL_MAX_ARRAY()
        IMPLICIT NONE
        IF (allocated( BUF_MAX_ARRAY)) DEALLOCATE( BUF_MAX_ARRAY )
        RETURN
        END SUBROUTINE SMUMPS_BUF_DEALL_MAX_ARRAY
        SUBROUTINE SMUMPS_BUF_MAX_ARRAY_MINSIZE(NFS4FATHER,IERR)
        IMPLICIT NONE
        INTEGER IERR, NFS4FATHER
        IERR = 0
        IF (allocated( BUF_MAX_ARRAY)) THEN
          IF (BUF_LMAX_ARRAY .GE. NFS4FATHER) RETURN
          DEALLOCATE( BUF_MAX_ARRAY )
        ENDIF
        BUF_LMAX_ARRAY=max(1,NFS4FATHER)
        ALLOCATE(BUF_MAX_ARRAY(BUF_LMAX_ARRAY),stat=IERR)
        IF ( IERR .GT. 0 ) THEN
           IERR = -1
           RETURN
        END IF
        RETURN
        END SUBROUTINE SMUMPS_BUF_MAX_ARRAY_MINSIZE
        SUBROUTINE SMUMPS_BUF_SEND_CB( NBROWS_ALREADY_SENT,
     &                                INODE, FPERE, NFRONT, LCONT,
     &                                NASS, NPIV,
     &                                IWROW, IWCOL, A, PACKED_CB,
     &                                DEST, TAG, COMM, KEEP, IERR )
        IMPLICIT NONE
        INTEGER DEST, TAG, COMM, IERR
        INTEGER NBROWS_ALREADY_SENT
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INTEGER INODE, FPERE, NFRONT, LCONT, NASS, NPIV 
        INTEGER IWROW( LCONT ), IWCOL( LCONT )
        REAL A( * )
        LOGICAL PACKED_CB
        INCLUDE 'mpif.h'
        INTEGER :: IERR_MPI
        INTEGER NBROWS_PACKET
        INTEGER POSITION, IREQ, IPOS, I, J1
        INTEGER SIZE1, SIZE2, SIZE_PACK, SIZE_AV, SIZE_AV_REALS
        INTEGER IZERO, IONE
        INTEGER SIZECB
        INTEGER LCONT_SENT
        INTEGER DEST2(1)
        PARAMETER( IZERO = 0, IONE = 1 )
        LOGICAL RECV_BUF_SMALLER_THAN_SEND
        DOUBLE PRECISION TMP
        DEST2(1) = DEST
        IERR = 0
        IF (NBROWS_ALREADY_SENT .EQ. 0) THEN
          CALL MPI_PACK_SIZE( 11 + LCONT + LCONT, MPI_INTEGER,
     &                        COMM, SIZE1,  IERR_MPI)
        ELSE
          CALL MPI_PACK_SIZE( 5, MPI_INTEGER, COMM, SIZE1, IERR_MPI)
        ENDIF
        CALL MUMPS_BUF_SIZE_AVAILABLE( BUF_CB, SIZE_AV )
        IF ( SIZE_AV .LT. SIZE_RBUF_BYTES ) THEN
          RECV_BUF_SMALLER_THAN_SEND = .FALSE.
        ELSE
          SIZE_AV = SIZE_RBUF_BYTES
          RECV_BUF_SMALLER_THAN_SEND = .TRUE.
        ENDIF
        SIZE_AV_REALS = ( SIZE_AV - SIZE1 ) / SIZEofREAL
        IF (SIZE_AV_REALS < 0 ) THEN
          NBROWS_PACKET = 0
        ELSE
          IF (PACKED_CB) THEN
            TMP=2.0D0*dble(NBROWS_ALREADY_SENT)+1.0D0
            NBROWS_PACKET = int(
     &                      ( sqrt( TMP * TMP
     &                        + 8.0D0 * dble(SIZE_AV_REALS)) - TMP )
     &                        / 2.0D0 )
          ELSE
            IF (LCONT.EQ.0) THEN
              NBROWS_PACKET = 0
            ELSE
              NBROWS_PACKET = SIZE_AV_REALS / LCONT
            ENDIF
          ENDIF
        ENDIF
 10     CONTINUE
        NBROWS_PACKET = max(0,
     &            min(NBROWS_PACKET, LCONT - NBROWS_ALREADY_SENT))
        IF (NBROWS_PACKET .EQ. 0 .AND. LCONT .NE. 0) THEN
           IF (RECV_BUF_SMALLER_THAN_SEND) THEN
            IERR = -3
            GOTO 100
         ELSE
            IERR = -1
            GOTO 100
          ENDIF
        ENDIF
        IF (PACKED_CB) THEN
          SIZECB = (NBROWS_ALREADY_SENT*NBROWS_PACKET)+(NBROWS_PACKET
     &             *(NBROWS_PACKET+1))/2
        ELSE
          SIZECB = NBROWS_PACKET * LCONT
        ENDIF
        CALL MPI_PACK_SIZE( SIZECB, MPI_REAL,
     &                    COMM, SIZE2,  IERR_MPI )
        SIZE_PACK = SIZE1 + SIZE2
        IF (SIZE_PACK .GT. SIZE_AV ) THEN
          NBROWS_PACKET = NBROWS_PACKET - 1
          IF (NBROWS_PACKET > 0) THEN
             GOTO 10
          ELSE
             IF (RECV_BUF_SMALLER_THAN_SEND) THEN
               IERR=-3
               GOTO 100
            ELSE
               IERR = -1
               GOTO 100
             ENDIF
          ENDIF
        ENDIF
        IF (NBROWS_PACKET + NBROWS_ALREADY_SENT.NE.LCONT .AND.
     &     SIZE_PACK  .LT. SIZE_RBUF_BYTES / 10
     &    .AND. 
     &    .NOT. RECV_BUF_SMALLER_THAN_SEND)
     &       THEN
            IERR = -1
            GOTO 100
        ENDIF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE_PACK, IERR, 
     &                 IONE , DEST2)
        IF (IERR .EQ. -1 .OR. IERR .EQ. -2) THEN
          NBROWS_PACKET = NBROWS_PACKET - 1
          IF ( NBROWS_PACKET > 0 )  GOTO 10
        ENDIF
        IF ( IERR .LT. 0 ) GOTO 100
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( FPERE, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
        IF (PACKED_CB) THEN
          LCONT_SENT=-LCONT
        ELSE
          LCONT_SENT=LCONT
        ENDIF
        CALL MPI_PACK( LCONT_SENT, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NBROWS_ALREADY_SENT, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NBROWS_PACKET, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
        IF (NBROWS_ALREADY_SENT == 0) THEN
          CALL MPI_PACK( LCONT, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( NASS-NPIV, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( LCONT , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( IZERO, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( IONE,  1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( IZERO, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( IWROW, LCONT, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( IWCOL, LCONT, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
        ENDIF
        IF ( LCONT .NE. 0 ) THEN
          J1 = 1 + NBROWS_ALREADY_SENT * NFRONT
          IF (PACKED_CB) THEN
           DO I = NBROWS_ALREADY_SENT+1,
     &            NBROWS_ALREADY_SENT+NBROWS_PACKET
            CALL MPI_PACK( A( J1 ), I, MPI_REAL,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
             J1 = J1 + NFRONT
           END DO
          ELSE
           DO I = NBROWS_ALREADY_SENT+1,
     &            NBROWS_ALREADY_SENT+NBROWS_PACKET
            CALL MPI_PACK( A( J1 ), LCONT, MPI_REAL,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
             J1 = J1 + NFRONT
           END DO
          ENDIF
        END IF
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                  DEST, TAG, COMM, BUF_CB%CONTENT( IREQ ),
     &                  IERR_MPI )
        IF ( SIZE_PACK .LT. POSITION ) THEN
          WRITE(*,*) 'Error Try_send_cb: SIZE, POSITION=',SIZE_PACK,
     &               POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE_PACK .NE. POSITION )
     &    CALL BUF_ADJUST( BUF_CB, POSITION )
        NBROWS_ALREADY_SENT = NBROWS_ALREADY_SENT + NBROWS_PACKET
        IF (NBROWS_ALREADY_SENT .NE. LCONT ) THEN
          IERR = -1
          RETURN
        ENDIF
 100    CONTINUE
        RETURN
        END SUBROUTINE SMUMPS_BUF_SEND_CB
        SUBROUTINE SMUMPS_BUF_SEND_MASTER2SLAVE( NRHS, INODE, IFATH,
     &             EFF_CB_SIZE, LD_CB, LD_PIV, NPIV, 
     &             JBDEB, JBFIN,
     &             CB, SOL,
     &             DEST, COMM, KEEP, IERR )
        IMPLICIT NONE
        INTEGER NRHS, INODE, IFATH, EFF_CB_SIZE, LD_CB, LD_PIV, NPIV 
        INTEGER DEST, COMM, IERR, JBDEB, JBFIN
        REAL CB( LD_CB*(NRHS-1)+EFF_CB_SIZE )
        REAL SOL( max(1, LD_PIV*(NRHS-1)+NPIV) )
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER SIZE, SIZE1, SIZE2, K
        INTEGER POSITION, IREQ, IPOS
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1) = DEST
        IERR = 0
        CALL MPI_PACK_SIZE( 6, MPI_INTEGER, COMM, SIZE1, IERR )
        CALL MPI_PACK_SIZE( NRHS * (EFF_CB_SIZE + NPIV),
     &                      MPI_REAL, COMM,
     &                      SIZE2, IERR_MPI )
        SIZE = SIZE1 + SIZE2
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
      RETURN
        ENDIF
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( IFATH, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( EFF_CB_SIZE  , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NPIV , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( JBDEB , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( JBFIN , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        DO K = 1, NRHS
               CALL MPI_PACK( CB ( 1 + LD_CB * (K-1) ),
     &                        EFF_CB_SIZE, MPI_REAL,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        END DO
        IF ( NPIV .GT. 0 ) THEN
          DO K=1, NRHS
          CALL MPI_PACK( SOL(1+LD_PIV*(K-1)),
     &                         NPIV, MPI_REAL,
     &                         BUF_CB%CONTENT( IPOS ), SIZE,
     &                         POSITION, COMM, IERR_MPI )
          ENDDO
        END IF
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                  DEST, Master2Slave, COMM,
     &                  BUF_CB%CONTENT( IREQ ), IERR_MPI )
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) 'Try_send_master2slave: SIZE, POSITION = ',
     &               SIZE, POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION ) CALL BUF_ADJUST( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE SMUMPS_BUF_SEND_MASTER2SLAVE
        SUBROUTINE SMUMPS_BUF_SEND_VCB( NRHS_B, NODE1, NODE2, NCB, LDW,
     &             LONG,
     &             IW, W, JBDEB, JBFIN,
     &             RHSINTR, NRHS, LRHSINTR, IPOSINRHSINTR, NPIV,
     &             KEEP,
     &             DEST, TAG, COMM, IERR )
        IMPLICIT NONE
        INTEGER LDW, DEST, TAG, COMM, IERR
        INTEGER NRHS_B, NODE1, NODE2, NCB, LONG, JBDEB, JBFIN
        INTEGER IW( max( 1, LONG ) )
        INTEGER, INTENT(IN) :: LRHSINTR, NRHS, IPOSINRHSINTR, NPIV
        REAL W( max( 1, LDW * NRHS_B ) )
        REAL RHSINTR(LRHSINTR,NRHS)
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INCLUDE 'mpif.h'
        INTEGER :: IERR_MPI
        INTEGER POSITION, IREQ, IPOS
        INTEGER SIZE1, SIZE2, SIZE, K
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1)=DEST
        IERR = 0
        IF ( NODE2 .EQ. 0 ) THEN
         CALL MPI_PACK_SIZE( 4+LONG, MPI_INTEGER, COMM, SIZE1,
     &                       IERR_MPI )
        ELSE
         CALL MPI_PACK_SIZE( 6+LONG, MPI_INTEGER, COMM, SIZE1,
     &                       IERR_MPI )
        END IF
        SIZE2 = 0
        IF ( LONG .GT. 0 ) THEN
          CALL MPI_PACK_SIZE( NRHS_B*LONG, MPI_REAL,
     &                        COMM, SIZE2, IERR_MPI )
        END IF
        SIZE = SIZE1 + SIZE2
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        POSITION = 0
        CALL MPI_PACK( NODE1, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        IF ( NODE2 .NE. 0 ) THEN
          CALL MPI_PACK( NODE2, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( NCB, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        ENDIF
        CALL MPI_PACK( JBDEB, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( JBFIN, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( LONG,  1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        IF ( LONG .GT. 0 ) THEN
          CALL MPI_PACK( IW, LONG, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
          IF (NODE2.EQ.0) THEN
            DO K=1, NRHS_B
              IF (NPIV.GT.0) THEN
              CALL MPI_PACK( RHSINTR(IPOSINRHSINTR,JBDEB+K-1), NPIV,
     &                          MPI_REAL,
     &                          BUF_CB%CONTENT( IPOS ), SIZE,
     &                          POSITION, COMM, IERR_MPI )
              ENDIF
              IF (LONG-NPIV .NE.0) THEN
                CALL MPI_PACK( W(NPIV+1+(K-1)*LDW), LONG-NPIV,
     &                          MPI_REAL,
     &                          BUF_CB%CONTENT( IPOS ), SIZE,
     &                          POSITION, COMM, IERR_MPI )
              ENDIF
            END DO
          ELSE
            DO K=1, NRHS_B
              CALL MPI_PACK( W(1+(K-1)*LDW), LONG, MPI_REAL,
     &                          BUF_CB%CONTENT( IPOS ), SIZE,
     &                          POSITION, COMM, IERR_MPI )
            END DO
          ENDIF
        END IF
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                  DEST, TAG, COMM, BUF_CB%CONTENT( IREQ ),
     &                  IERR_MPI )
        IF ( SIZE .NE. POSITION ) CALL BUF_ADJUST( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE SMUMPS_BUF_SEND_VCB
        SUBROUTINE SMUMPS_BUF_SEND_MAITRE2( NBROWS_ALREADY_SENT,
     &  IPERE, ISON, NROW,
     &  IROW, NCOL, ICOL, VAL, LDA, NELIM, TYPE_SON,
     &  NSLAVES, SLAVES, DEST, COMM, IERR, 
     & 
     &  SLAVEF, KEEP,KEEP8, INIV2, TAB_POS_IN_PERE )
        IMPLICIT NONE
        INTEGER NBROWS_ALREADY_SENT
        INTEGER LDA, NELIM, TYPE_SON
        INTEGER IPERE, ISON, NROW, NCOL, NSLAVES
        INTEGER IROW( NROW )
        INTEGER ICOL( NCOL )
        INTEGER SLAVES( NSLAVES )
        REAL VAL(LDA, *)
        INTEGER IPOS, IREQ, DEST, COMM, IERR
        INTEGER SLAVEF, KEEP(500), INIV2
        INTEGER(8) KEEP8(150)
        INTEGER TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER SIZE1, SIZE2, SIZE3, SIZE_PACK, POSITION, I
        INTEGER NBROWS_PACKET, NCOL_SEND
        INTEGER SIZE_AV
        LOGICAL RECV_BUF_SMALLER_THAN_SEND
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1) = DEST
        IERR = 0
        IF ( NELIM .NE. NROW ) THEN
          WRITE(*,*) 'Error in TRY_SEND_MAITRE2:',NELIM, NROW
          CALL MUMPS_ABORT()
        END IF
        IF (NBROWS_ALREADY_SENT .EQ. 0) THEN
          CALL MPI_PACK_SIZE( NROW+NCOL+7+NSLAVES, MPI_INTEGER,
     &                      COMM, SIZE1, IERR_MPI )
          IF ( TYPE_SON .eq. 2 ) THEN
          CALL MPI_PACK_SIZE( NSLAVES+1, MPI_INTEGER,
     &                          COMM, SIZE3, IERR_MPI )
          ELSE
            SIZE3 = 0
          ENDIF
          SIZE1=SIZE1+SIZE3
        ELSE
          CALL MPI_PACK_SIZE(7, MPI_INTEGER,COMM,SIZE1,IERR_MPI)
        ENDIF
        IF ( KEEP(50).ne.0  .AND. TYPE_SON .eq. 2 ) THEN
          NCOL_SEND = NROW
        ELSE
          NCOL_SEND = NCOL
        ENDIF
        CALL MUMPS_BUF_SIZE_AVAILABLE( BUF_CB, SIZE_AV )
        IF (SIZE_AV .LT. SIZE_RBUF_BYTES) THEN
          RECV_BUF_SMALLER_THAN_SEND = .FALSE.
        ELSE
          RECV_BUF_SMALLER_THAN_SEND = .TRUE.
          SIZE_AV = SIZE_RBUF_BYTES
        ENDIF
        IF (NROW .GT. 0 ) THEN 
         NBROWS_PACKET = (SIZE_AV - SIZE1) / NCOL_SEND / SIZEofREAL
         NBROWS_PACKET = min(NBROWS_PACKET, NROW - NBROWS_ALREADY_SENT)
         NBROWS_PACKET = max(NBROWS_PACKET, 0)
        ELSE
          NBROWS_PACKET =0
        ENDIF
        IF (NBROWS_PACKET .EQ. 0 .AND. NROW .NE. 0) THEN
           IF (RECV_BUF_SMALLER_THAN_SEND) THEN
              IERR=-3
              GOTO 100
           ELSE
              IERR=-1
              GOTO 100
          ENDIF
        ENDIF
 10     CONTINUE
        CALL MPI_PACK_SIZE( NBROWS_PACKET * NCOL_SEND,
     &           MPI_REAL,
     &           COMM, SIZE2, IERR_MPI )
        SIZE_PACK = SIZE1 + SIZE2
        IF (SIZE_PACK .GT. SIZE_AV) THEN
          NBROWS_PACKET = NBROWS_PACKET - 1
          IF ( NBROWS_PACKET .GT. 0 ) THEN
            GOTO 10
          ELSE
            IF (RECV_BUF_SMALLER_THAN_SEND) THEN
                IERR = -3
                GOTO 100
             ELSE
                IERR = -1
                GOTO 100
            ENDIF
          ENDIF
        ENDIF
       IF (NBROWS_PACKET + NBROWS_ALREADY_SENT.NE.NROW .AND.
     &   SIZE_PACK - SIZE1  .LT. ( SIZE_RBUF_BYTES - SIZE1 ) / 10
     &  .AND. 
     &   .NOT. RECV_BUF_SMALLER_THAN_SEND)
     &       THEN
           IERR = -1
           GOTO 100
       ENDIF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE_PACK, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
          GOTO 100
        ENDIF
        POSITION = 0
        CALL MPI_PACK( IPERE, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( ISON,  1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NSLAVES, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NROW, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NCOL, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NBROWS_ALREADY_SENT, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NBROWS_PACKET, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI )
        IF (NBROWS_ALREADY_SENT .EQ. 0) THEN
          IF (NSLAVES.GT.0) THEN
            CALL MPI_PACK( SLAVES, NSLAVES, MPI_INTEGER,
     &                BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                POSITION, COMM, IERR_MPI )
          ENDIF
          CALL MPI_PACK( IROW, NROW, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( ICOL, NCOL, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI )
          IF ( TYPE_SON .eq. 2 ) THEN
            CALL MPI_PACK( TAB_POS_IN_PERE(1,INIV2), NSLAVES+1, 
     &                 MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI )
          ENDIF
        ENDIF
        IF (NBROWS_PACKET.GE.1) THEN
          DO I=NBROWS_ALREADY_SENT+1,
     &                   NBROWS_ALREADY_SENT+NBROWS_PACKET
            CALL MPI_PACK( VAL(1,I), NCOL_SEND, 
     &               MPI_REAL,
     &               BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &               POSITION, COMM, IERR_MPI )
          ENDDO
        ENDIF
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                  DEST, MAITRE2, COMM,
     &                  BUF_CB%CONTENT( IREQ ), IERR_MPI )
        IF ( SIZE_PACK .LT. POSITION ) THEN
          write(*,*) 'Try_send_maitre2, SIZE,POSITION=',
     &                SIZE_PACK,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE_PACK .NE. POSITION )
     &    CALL BUF_ADJUST( BUF_CB, POSITION )
        NBROWS_ALREADY_SENT = NBROWS_ALREADY_SENT + NBROWS_PACKET
        IF ( NBROWS_ALREADY_SENT .NE. NROW ) THEN
          IERR = -1
        ENDIF
 100    CONTINUE
        RETURN
      END SUBROUTINE SMUMPS_BUF_SEND_MAITRE2
      SUBROUTINE SMUMPS_BUF_SEND_CONTRIB_TYPE2(NBROWS_ALREADY_SENT,
     &                 NBCOLS_ALREADY_SENT, NBLRB_ALREADY_SENT,
     &  DESC_IN_LU,
     &  IPERE, NFRONT_PERE, NASS_PERE, NFS4FATHER, NSLAVES_PERE,
     &  ISON, NBROW, LMAP, MAPROW, POS_FIRST_ROW_TO_PDEST, 
     &  IW_CBSON, A_CBSON, LA_CBSON,
     &  ISLAVE, PDEST, PDEST_MASTER, COMM, 
     &  IERR, IERR_SIZEREQ, LP,
     &  
     &  KEEP,KEEP8, STEP, N, SLAVEF,
     &  ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &  PACKED_CB, KEEP253_LOC, NVSCHUR,
     &  SON_NIV, MYID )
        USE SMUMPS_LR_TYPE
        USE SMUMPS_LR_DATA_M
        USE MUMPS_BUF_COMMON
        IMPLICIT NONE
        INTEGER NBROWS_ALREADY_SENT
        INTEGER, INTENT(inout)::  NBCOLS_ALREADY_SENT, 
     &                            NBLRB_ALREADY_SENT
        INTEGER, INTENT (in) :: KEEP253_LOC, NVSCHUR
        INTEGER, INTENT (in) :: SON_NIV
        INTEGER, INTENT(in)  :: POS_FIRST_ROW_TO_PDEST
        INTEGER IPERE, ISON, NBROW, MYID
        INTEGER PDEST, ISLAVE, COMM, IERR
        INTEGER, INTENT (in)  :: LP
        INTEGER, INTENT (out) :: IERR_SIZEREQ
        INTEGER PDEST_MASTER, NASS_PERE, NSLAVES_PERE,
     &       NFRONT_PERE, LMAP
        INTEGER MAPROW( LMAP )
        INTEGER IW_CBSON( * )
        REAL A_CBSON( : )
        INTEGER(8) :: LA_CBSON
        LOGICAL DESC_IN_LU, PACKED_CB
        INTEGER   KEEP(500), N , SLAVEF
        INTEGER(8) KEEP8(150)
        INTEGER   STEP(N), 
     &            ISTEP_TO_INIV2(KEEP(71)), 
     &            TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: IERR_MPI
      INTEGER NFS4FATHER,SIZE3,PS1,NCA,LROW1
      INTEGER(8) :: ASIZE
      LOGICAL COMPUTE_MAX
      REAL, POINTER, DIMENSION(:) :: M_ARRAY
      INTEGER NBROWS_PACKET 
      INTEGER NBLRB_TOTAL
      INTEGER NBLRB_PACKET 
      INTEGER MAX_ROW_LENGTH
      INTEGER LROW, NELIM
      INTEGER(8) :: ITMP8
      INTEGER NPIV, NFRONT, HS
      INTEGER SIZE_PACK, SIZE0, SIZE1, SIZE2, POSITION,I
      INTEGER SIZE_INTEGERS, B, SIZE_REALS, TMPSIZE, ONEorTWO, SIZE_AV
      INTEGER SIZE_NEXT_BLOCK
      INTEGER SIZE_LRMIXED, SIZE4
      INTEGER NBINT, L
      INTEGER(8) :: APOS, SHIFTCB_SON, LDA_SON8
      INTEGER IPOS_IN_SLAVE
      INTEGER STATE_SON
      INTEGER INDICE_PERE, NROW, IPOS, IREQ, NOSLA
      INTEGER IONE, J, THIS_ROW_LENGTH
      INTEGER SIZE_DESC_BANDE, DESC_BANDE_BYTES
      LOGICAL RECV_BUF_SMALLER_THAN_SEND
      LOGICAL NOT_ENOUGH_SPACE
      LOGICAL AVOID_TOO_SMALL_GRANULARITY
      INTEGER PDEST2(1)
      LOGICAL CB_IS_LR
      TYPE(LRB_TYPE), POINTER :: CB_LRB(:,:)
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_ROW, BEGS_BLR_COL,
     &                    BEGS_BLR_STA
      INTEGER :: NB_ROW_SHIFT, NB_COL_SHIFT, NASS_SHIFT, PANEL2SEND,
     &           CURRENT_PANEL_SIZE, NB_BLR_ROWS, NB_BLR_COLS,
     &           CB_IS_LR_INT, NCOL_SHIFT, NROW_SHIFT,
     &           PANEL_BEG_OFFSET
      INTEGER :: NPIV_LR, LNEXT
      REAL    :: K170PER1000
      PARAMETER ( IONE=1 )
      INCLUDE 'mumps_headers.h'
      REAL ZERO
      PARAMETER (ZERO = 0.0E0)
      CB_IS_LR = (IW_CBSON(1+XXLR).EQ.1
     &       .OR. IW_CBSON(1+XXLR).EQ.3)
      NBLRB_PACKET = 0
      NBLRB_TOTAL  = 0
      IF (CB_IS_LR) THEN
        CB_IS_LR_INT = 1
      ELSE
        CB_IS_LR_INT = 0
      ENDIF
      AVOID_TOO_SMALL_GRANULARITY =  .TRUE.
      IF (KEEP(173).EQ.1) AVOID_TOO_SMALL_GRANULARITY = .FALSE.
      COMPUTE_MAX = (KEEP(219) .NE. 0) .AND.
     &              (KEEP(50) .EQ. 2) .AND.
     &              (PDEST.EQ.PDEST_MASTER)
     &       .AND. (NBCOLS_ALREADY_SENT.EQ.0)
     &       .AND. (NBROWS_ALREADY_SENT.EQ.0)
      LNEXT          = 0
      MAX_ROW_LENGTH = 0
      IERR_SIZEREQ   = 0
      B              = 0
      ONEorTWO       = 1
      IF (COMPUTE_MAX) THEN
        CALL SMUMPS_BUF_MAX_ARRAY_MINSIZE(NFS4FATHER,IERR)
        IF (IERR .NE. 0) THEN
          IERR         = -4
          IERR_SIZEREQ = max(1,NFS4FATHER)
          RETURN
        ENDIF
      ENDIF
      PDEST2(1) = PDEST
      IERR   = 0
      LROW   = IW_CBSON( 1 + KEEP(IXSZ))
      NELIM  = IW_CBSON( 2 + KEEP(IXSZ))
      NPIV   = IW_CBSON( 4 + KEEP(IXSZ))
      IF ( NPIV .LT. 0 ) THEN
          NPIV = 0
      END IF
      NROW   = IW_CBSON( 3 + KEEP(IXSZ))
      NFRONT = LROW + NPIV
      HS     = 6 + IW_CBSON( 6 + KEEP(IXSZ)) + KEEP(IXSZ)
      IF (CB_IS_LR.AND.NBROW.GT.0) THEN
        CALL SMUMPS_BLR_RETRIEVE_CB_LRB(IW_CBSON(1+XXF), CB_LRB)
        IF (SON_NIV.EQ.1) THEN
          CALL SMUMPS_BLR_RETRIEVE_BEGSBLR_STA(IW_CBSON(1+XXF),
     &                     BEGS_BLR_ROW)
          CALL SMUMPS_BLR_RETRIEVE_BEGSBLR_DYN(IW_CBSON(1+XXF),
     &                     BEGS_BLR_COL)
          NB_BLR_ROWS = size(BEGS_BLR_ROW) - 1
          CALL SMUMPS_BLR_RETRIEVE_NB_PANELS(IW_CBSON(1+XXF), 
     &                    NB_COL_SHIFT)
          NB_ROW_SHIFT = NB_COL_SHIFT
          NASS_SHIFT = BEGS_BLR_ROW(NB_ROW_SHIFT+1)-1
          NPIV_LR  = BEGS_BLR_COL(NB_COL_SHIFT+1)-1
        ELSE
          NPIV_LR=NPIV
          CALL SMUMPS_BLR_RETRIEVE_BEGSBLR_STA(IW_CBSON(1+XXF),
     &                     BEGS_BLR_STA)
          NB_BLR_ROWS = size(BEGS_BLR_STA) - 2
          BEGS_BLR_ROW => BEGS_BLR_STA(2:NB_BLR_ROWS+2)
          CALL SMUMPS_BLR_RETRIEVE_BEGS_BLR_C(IW_CBSON(1+XXF),
     &                     BEGS_BLR_COL, NB_COL_SHIFT)
          NASS_SHIFT = 0
          NB_ROW_SHIFT = 0
        ENDIF
        PANEL2SEND             = -1
        DO I=NB_ROW_SHIFT+1,NB_BLR_ROWS
          IF (BEGS_BLR_ROW(I+1)-1-NASS_SHIFT
     &         .GT.NBROWS_ALREADY_SENT+POS_FIRST_ROW_TO_PDEST-1) THEN
            PANEL2SEND = I
            EXIT
          ENDIF
        ENDDO
        IF (PANEL2SEND.EQ.-1) THEN
          write(*,*) 'Internal error: PANEL2SEND not found'
          CALL MUMPS_ABORT()
        ENDIF
        CURRENT_PANEL_SIZE = BEGS_BLR_ROW(PANEL2SEND+1)
     &                     - BEGS_BLR_ROW(PANEL2SEND) 
        PANEL_BEG_OFFSET = POS_FIRST_ROW_TO_PDEST + 
     &                     NBROWS_ALREADY_SENT -
     &                     BEGS_BLR_ROW(PANEL2SEND) + NASS_SHIFT
        IF (KEEP(50).EQ.0) THEN
          NB_BLR_COLS = size(BEGS_BLR_COL)  - 1
        ELSEIF (SON_NIV.EQ.1) THEN
          NB_BLR_COLS = PANEL2SEND
        ELSE
          NB_BLR_COLS = -1
          NCOL_SHIFT = NPIV_LR
          NROW_SHIFT = LROW - NROW
          DO I=NB_COL_SHIFT+1,size(BEGS_BLR_COL)-1
             IF (BEGS_BLR_COL(I+1)-NCOL_SHIFT.GT.
     &        (  min (
     &           BEGS_BLR_ROW(PANEL2SEND+1)-POS_FIRST_ROW_TO_PDEST, 
     &           NBROW
     &               )
     &              + NROW_SHIFT + POS_FIRST_ROW_TO_PDEST-1 )
     &          ) THEN
              NB_BLR_COLS = I
              EXIT
            ENDIF
          ENDDO
          IF (NB_BLR_COLS.EQ.-1) THEN
            write(*,*) 'Internal error: NB_BLR_COLS not found'
            CALL MUMPS_ABORT()
          ENDIF
          MAX_ROW_LENGTH = BEGS_BLR_ROW(PANEL2SEND)-1+NROW_SHIFT
     &        + min(NBROW-NBROWS_ALREADY_SENT + PANEL_BEG_OFFSET,
     &          CURRENT_PANEL_SIZE)
        ENDIF
        NBLRB_TOTAL = NB_BLR_COLS - NB_COL_SHIFT
      ENDIF
      STATE_SON = IW_CBSON(1+XXS)
      IF (STATE_SON .EQ. S_NOLCBCONTIG) THEN
               LDA_SON8    = int(LROW,8)
               SHIFTCB_SON = int(NPIV,8)*int(NROW,8)
      ELSE IF (STATE_SON .EQ. S_NOLCLEANED) THEN
               LDA_SON8    = int(LROW,8)
               SHIFTCB_SON = 0_8
      ELSE
               LDA_SON8     = int(NFRONT,8)
               SHIFTCB_SON = int(NPIV,8)
      ENDIF
      CALL MUMPS_BUF_SIZE_AVAILABLE( BUF_CB, SIZE_AV )
      IF (PDEST .EQ. PDEST_MASTER) THEN
        SIZE_DESC_BANDE=0 
      ELSE
        SIZE_DESC_BANDE=(11+SLAVEF+KEEP(127)*2)
        SIZE_DESC_BANDE=SIZE_DESC_BANDE+int(real(KEEP(12))*
     &                  real(SIZE_DESC_BANDE)/100.0E0)
        SIZE_DESC_BANDE=max(SIZE_DESC_BANDE,
     &     11+NSLAVES_PERE+NFRONT_PERE+NFRONT_PERE-NASS_PERE)
      ENDIF
      DESC_BANDE_BYTES=SIZE_DESC_BANDE*SIZEofINT
      IF ( SIZE_AV .LT. SIZE_RBUF_BYTES-DESC_BANDE_BYTES ) THEN
        RECV_BUF_SMALLER_THAN_SEND = .FALSE.
      ELSE
        RECV_BUF_SMALLER_THAN_SEND = .TRUE.
        SIZE_AV = SIZE_RBUF_BYTES-DESC_BANDE_BYTES
      ENDIF
      SIZE1=0
      IF(COMPUTE_MAX) THEN
               CALL MPI_PACK_SIZE(1, MPI_INTEGER,
     &            COMM, SIZE0, IERR_MPI )
               IF(NFS4FATHER .GT. 0) THEN
                CALL MPI_PACK_SIZE( NFS4FATHER, MPI_REAL,
     &             COMM, SIZE1, IERR_MPI )
               ENDIF
               SIZE1 = SIZE1+SIZE0
      ENDIF
      ONEorTWO = 1
      IF (PDEST .EQ.PDEST_MASTER) THEN
        L = 0
      ELSE IF (KEEP(50) .EQ. 0) THEN
        L = LROW
      ELSE
        L = LROW + POS_FIRST_ROW_TO_PDEST-LMAP+NBROWS_ALREADY_SENT-1
        ONEorTWO=ONEorTWO+1
      ENDIF
      NBINT = 6 + L + 1 
      IF (CB_IS_LR.AND.NBROW.GT.0) THEN
        NBINT = NBINT + 4*(NB_BLR_COLS-NB_COL_SHIFT) + 3
      ENDIF
      CALL MPI_PACK_SIZE( NBINT, MPI_INTEGER,
     &                    COMM, TMPSIZE, IERR_MPI )
      SIZE1 = SIZE1 + TMPSIZE
      SIZE_AV = SIZE_AV - SIZE1
      NOT_ENOUGH_SPACE=.FALSE.
      IF (SIZE_AV .LT.0 ) THEN
        NBROWS_PACKET = 0
        NOT_ENOUGH_SPACE=.TRUE.
      ELSE
        IF ( KEEP(50) .EQ. 0 ) THEN
          NBROWS_PACKET =
     &       SIZE_AV / ( ONEorTWO*SIZEofINT+LROW*SIZEofREAL)
        ELSE
          B = 2 * ONEorTWO + 
     &      ( -1 + 2 *  LROW + 2 * POS_FIRST_ROW_TO_PDEST -2*LMAP
     &                                  + 2 * NBROWS_ALREADY_SENT )
     &      * SIZEofREAL / SIZEofINT
          NBROWS_PACKET=int(
     &        (dble(-B)+
     &         sqrt
     &         ( (dble(B)*dble(B))+
     &           dble(4)*dble(2)*dble(SIZE_AV)/dble(SIZEofINT) *
     &           dble(SIZEofREAL)/dble(SIZEofINT)
     &         )
     &        )*
     &        dble(SIZEofINT) / dble(2) / dble(SIZEofREAL))
        ENDIF
      ENDIF
 10   CONTINUE
      SIZE_NEXT_BLOCK = 0
      IF (CB_IS_LR) THEN
        IF ( NBROW .GT. 0) THEN
          NBROWS_PACKET = CURRENT_PANEL_SIZE-PANEL_BEG_OFFSET
        ELSE
          NBROWS_PACKET = 0
        ENDIF
      ENDIF
      NBROWS_PACKET = max( 0, NBROWS_PACKET)
      NBROWS_PACKET = min(NBROW-NBROWS_ALREADY_SENT, NBROWS_PACKET)
      NOT_ENOUGH_SPACE = NOT_ENOUGH_SPACE .OR.
     &                   (NBROWS_PACKET .EQ.0.AND. NBROW.NE.0)
      IF (CB_IS_LR.AND.NBROW.GT.0.AND..NOT.NOT_ENOUGH_SPACE) THEN
        CALL MPI_PACK_SIZE( ONEorTWO* NBROWS_PACKET, MPI_INTEGER,
     &                    COMM, TMPSIZE, IERR_MPI )
        CALL SMUMPS_BLR_GET_SIZEREALS_CB_LRB(
     &       SIZE_AV-TMPSIZE, CB_LRB, 
     &       NB_ROW_SHIFT, PANEL2SEND,
     &       NBLRB_ALREADY_SENT, NBLRB_TOTAL, 
     &         NBLRB_PACKET, SIZE_REALS, SIZE_NEXT_BLOCK
     &         , SIZE_LRMIXED
     &         , (KEEP(452).GT.0) 
     &         , KEEP(173)
     &         )
          NOT_ENOUGH_SPACE = NOT_ENOUGH_SPACE.OR.(NBLRB_PACKET.EQ.0)
      ENDIF
      IF ( (NBROWS_PACKET + NBROWS_ALREADY_SENT.NE.NBROW) .AND.
     &      .NOT.CB_IS_LR
     &   ) THEN
       IF (KEEP(50).EQ.0) THEN
        LNEXT = LROW + 1
       ELSE
       MAX_ROW_LENGTH = LROW + POS_FIRST_ROW_TO_PDEST - LMAP
     &                    + NBROWS_ALREADY_SENT
     &                    + NBROWS_PACKET-1
        LNEXT = MAX_ROW_LENGTH + 1  
       ENDIF
       LNEXT = LNEXT + ONEorTWO 
       CALL MPI_PACK_SIZE( LNEXT, 
     &                     MPI_REAL,
     &                     COMM, SIZE_NEXT_BLOCK, IERR_MPI )
      ENDIF 
      IF (SIZE_NEXT_BLOCK.GT.SIZE_RBUF_BYTES) THEN
         IF (LP.GT.0) THEN
            write(LP, '(A/, A/, A, I5,I10,L2/, A, 2I12/,A,2I10)') 
     &      " Error 2 in  SMUMPS_BUF_SEND_CONTRIB_TYPE2,",
     &      " RECV buffer too small for next block",
     &      "    BUF_MYID, ISON,CB_IS_LR         = ", 
     &           BUF_MYID, ISON,CB_IS_LR,
     &      "    SIZE_NEXT_BLOCK, SIZE_RBUF_BYTES= ",
     &           SIZE_NEXT_BLOCK, SIZE_RBUF_BYTES,
     &      "    LNEXT, MAX_ROW_LENGTH           = ",
     &           LNEXT, MAX_ROW_LENGTH
            write(LP,'(A/,A,5I9/,A,4I12/,A,3I9/,A,2I9)')
     &       " Error 2 additional information :",
     &       "    NBROW, NBROWS_ALREADY_SENT, LROW, B, ONEorTWO    = ",
     &            NBROW, NBROWS_ALREADY_SENT, LROW, B, ONEorTWO, 
     &       "    SIZE_AV, SIZE1, SIZE_RBUF_BYTES, DESC_BANDE_BYTES= ",
     &            SIZE_AV, SIZE1, SIZE_RBUF_BYTES, DESC_BANDE_BYTES,
     &       "    NPIV,NROW,NFRONT           = ", NPIV,NROW,NFRONT,
     &       "    NBROWS_PACKET, NBLRB_PACKET= ",
     &            NBROWS_PACKET, NBLRB_PACKET
         ENDIF
         IERR = -3
         IERR_SIZEREQ = SIZE_NEXT_BLOCK
         RETURN
      ENDIF
      IF (NOT_ENOUGH_SPACE) THEN
         IF (RECV_BUF_SMALLER_THAN_SEND) THEN
          IERR = -3
          IERR_SIZEREQ =  SIZE1
          IF (LP.GT.0) THEN
            write(LP,'(A/,A/,A,I5,I9,L2/,A,2I4)')
     &      " Error 1 in SMUMPS_BUF_SEND_CONTRIB_TYPE2, ",
     &      " RECV buffer too small to send any real/complex data",
     &      "    BUF_MYID,ISON, CB_IS_LR = ",
     &           BUF_MYID,ISON, CB_IS_LR,
     &      "    SIZEofINT, SIZEofREAL   = ",
     &           SIZEofINT, SIZEofREAL
            write(LP,'(A/,A,5I9/,A,4I12/,A,3I9/,A,2I9)')
     &      " Error 1 additional information 1:",
     &      "    NBROW, NBROWS_ALREADY_SENT, LROW, B, ONEorTWO    = ",
     &           NBROW, NBROWS_ALREADY_SENT, LROW, B, ONEorTWO, 
     &      "    SIZE_AV, SIZE1, SIZE_RBUF_BYTES, DESC_BANDE_BYTES= ",
     &           SIZE_AV, SIZE1, SIZE_RBUF_BYTES, DESC_BANDE_BYTES,
     &      "    NPIV, NROW, NFRONT         = ", NPIV, NROW, NFRONT,
     &      "    NBROWS_PACKET, NBLRB_PACKET= ",
     &           NBROWS_PACKET, NBLRB_PACKET
            write(LP,'(A/,A,4I9/,A,3I9/,A,2I9)')
     &      " Error 1 additional information 2:",
     &      "    SON_NIV, NELIM, PDEST, PDEST_MASTER      = ",
     &           SON_NIV, NELIM, PDEST, PDEST_MASTER,
     &      "    POS_FIRST_ROW_TO_PDEST, LMAP, L          = ",
     &           POS_FIRST_ROW_TO_PDEST, LMAP, L, 
     &      "    MAX_ROW_LENGTH, LNEXT                    = ",
     &           MAX_ROW_LENGTH, LNEXT                    
          ENDIF
          GOTO 100
       ELSE
          IERR = -1
          GOTO 100
        ENDIF
      ENDIF
      IF (CB_IS_LR.AND.NBROW.GT.0) THEN
        IF (KEEP(50).EQ.0) THEN
          MAX_ROW_LENGTH = -99999
        ELSEIF (SON_NIV.EQ.1) THEN
          MAX_ROW_LENGTH = LROW+POS_FIRST_ROW_TO_PDEST -LMAP
     &                 + NBROWS_ALREADY_SENT
     &                 + NBROWS_PACKET-1
        ENDIF
      ELSE
        IF (KEEP(50).EQ.0) THEN
          MAX_ROW_LENGTH = -99999
          SIZE_REALS = NBROWS_PACKET * LROW
        ELSE
          SIZE_REALS = ( LROW + POS_FIRST_ROW_TO_PDEST - LMAP 
     &                   + NBROWS_ALREADY_SENT ) *
     &      NBROWS_PACKET + ( NBROWS_PACKET *
     &     ( NBROWS_PACKET - 1) ) / 2
          MAX_ROW_LENGTH = LROW + POS_FIRST_ROW_TO_PDEST - LMAP
     &                    + NBROWS_ALREADY_SENT
     &                    + NBROWS_PACKET-1
        ENDIF
      ENDIF
      SIZE_INTEGERS = ONEorTWO* NBROWS_PACKET
      CALL MPI_PACK_SIZE( SIZE_REALS, MPI_REAL,
     &                    COMM, SIZE2,  IERR_MPI )
      IF (CB_IS_LR .AND. (NBROW.GT.0) .AND. (KEEP(452).GT.0)) THEN
        CALL MPI_PACK_SIZE( SIZE_LRMIXED, MPI_INTEGER,
     &                    COMM, SIZE4,  IERR_MPI )
        SIZE2 = SIZE2 + SIZE4
      ENDIF
      CALL MPI_PACK_SIZE( SIZE_INTEGERS, MPI_INTEGER,
     &                    COMM, SIZE3,  IERR_MPI )
      IF (SIZE2 + SIZE3 .GT. SIZE_AV .AND. .NOT.CB_IS_LR) THEN
         NBROWS_PACKET = NBROWS_PACKET -1
         IF (NBROWS_PACKET .GT. 0 .AND..NOT.CB_IS_LR) THEN
           GOTO 10
         ENDIF
          IF (RECV_BUF_SMALLER_THAN_SEND) THEN
             IF (LP.GT.0) THEN
              WRITE (LP,'(A/,A,A/,A,I5,I9,L2)')
     &        " Error 5 in SMUMPS_BUF_SEND_CONTRIB_TYPE2",
     &        " DO NOT SEND: ",
     &        " Receive buffer minus DESC_BAND would be too small ",
     &        "    BUF_MYID, ISON, CB_IS_LR= ", BUF_MYID, ISON, CB_IS_LR
              WRITE(LP,'(A/,A,3I12/,A,3I12)')
     &        " Error 5 additional information 1:",
     &        "    Message size, size available, SIZE_RBUF_BYTES= ",
     &         SIZE1+SIZE2+SIZE3, SIZE_AV+SIZE1, SIZE_RBUF_BYTES,
     &        "    SIZE1, SIZE2, SIZE3                          = ", 
     &             SIZE1, SIZE2, SIZE3 
              write(LP,'(A/,A,5I9/,A,3I12/,A,3I9/,A,2I9)')
     &        " Error 5 additional information 2:",
     &        "    NBROW, NBROWS_ALREADY_SENT, LROW, B, ONEorTWO= ",
     &             NBROW, NBROWS_ALREADY_SENT, LROW, B, ONEorTWO, 
     &        "    SIZE_AV, SIZE1, DESC_BANDE_BYTES             = ",
     &             SIZE_AV, SIZE1, DESC_BANDE_BYTES,
     &        "    NPIV, NROW, NFRONT                           = ", 
     &             NPIV, NROW, NFRONT,
     &        "    NBROWS_PACKET, NBLRB_PACKET                  = ",
     &             NBROWS_PACKET, NBLRB_PACKET
             ENDIF
             IERR_SIZEREQ = SIZE1+SIZE2+SIZE3
             IERR = -3
             GOTO 100
          ELSE
             IERR = -1
             GOTO 100
           ENDIF
      ENDIF
        SIZE_PACK = SIZE1 + SIZE2 + SIZE3
        IERR_SIZEREQ = SIZE_PACK
        K170PER1000 = real(min(KEEP(170),500))/real(1000)
        IF ( .NOT.CB_IS_LR 
     &    .AND. (NBROWS_PACKET + NBROWS_ALREADY_SENT.NE.NBROW)
     &   .AND. .NOT. RECV_BUF_SMALLER_THAN_SEND 
     &   .AND. ( SIZE_PACK  .LT. 
     &         int(real(SIZE_RBUF_BYTES)*K170PER1000) )
     &   .AND.  
     &       ( int(SIZE_PACK,8) + int(SIZE_NEXT_BLOCK,8) .LE. 
     &      int(SIZE_RBUF_BYTES,8)  )
     &     ) THEN
          IERR = -1
          GOTO 100
        ENDIF
        IF ( CB_IS_LR.AND. 
     &      ( NBROWS_PACKET.NE.0 ).AND.
     &      ( NBLRB_ALREADY_SENT+NBLRB_PACKET.NE. NBLRB_TOTAL )
     &     .AND. ( SIZE_PACK  .LT. 
     &         int(real(SIZE_RBUF_BYTES)*K170PER1000) )
     &    .AND. .NOT. RECV_BUF_SMALLER_THAN_SEND 
     &    .AND. AVOID_TOO_SMALL_GRANULARITY
     &    .AND.  ( 
     &      int(SIZE_PACK,8) + int(SIZE_NEXT_BLOCK,8) .LE. 
     &      int(SIZE_RBUF_BYTES,8) )
     &     ) THEN
            IERR = -1
            GOTO 100
        ENDIF
        IF (SIZE_PACK.GT.SIZE_RBUF_BYTES ) THEN
          IERR         = -3
          IERR_SIZEREQ = SIZE_PACK
          IF (LP.GT.0) THEN
            WRITE(LP,'(A/,A/,A,I5,I10,L2)')
     &      " Error 4 in SMUMPS_BUF_SEND_CONTRIB_TYPE2",
     &      "    Receive buffer will be too small ",
     &      "    BUF_MYID, ISON, CB_IS_LR= ", BUF_MYID, ISON, CB_IS_LR
            WRITE(LP,'(A/,A,3I12/,A,4I12/,A,3I12)')
     &      " Error 4 additional information 1:",
     &      "    Message size, size available, SIZE_RBUF_BYTES= ",
     &           SIZE1+SIZE2+SIZE3, SIZE_AV+SIZE1, SIZE_RBUF_BYTES,
     &      "    SIZE_PACK, SIZE_RBUF_BYTES, ISON, IPERE      = ",
     &           SIZE_PACK, SIZE_RBUF_BYTES, ISON, IPERE, 
     &      "    SIZE1, SIZE2, SIZE3                          = ", 
     &           SIZE1, SIZE2, SIZE3 
            WRITE(LP,'(A/,A,5I9/,A,4I12/,A,3I9/,A,2I9)')
     &      " Error 4 additional information 2:",
     &      "    NBROW, NBROWS_ALREADY_SENT, LROW, B, ONEorTWO    = ",
     &           NBROW, NBROWS_ALREADY_SENT, LROW, B, ONEorTWO, 
     &      "    SIZE_AV, SIZE1, SIZE_RBUF_BYTES, DESC_BANDE_BYTES= ",
     &           SIZE_AV, SIZE1, SIZE_RBUF_BYTES, DESC_BANDE_BYTES,
     &      "    NPIV, NROW, NFRONT                               = ", 
     &           NPIV, NROW, NFRONT,
     &      "    NBROWS_PACKET, NBLRB_PACKET                      = ",
     &           NBROWS_PACKET, NBLRB_PACKET
           ENDIF
           GOTO 100
        ENDIF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE_PACK, IERR, 
     &                 IONE , PDEST2)
        IF (IERR .EQ. -1 .OR. IERR.EQ. -2) THEN
          IF (LP.GT.0) THEN
           WRITE(LP,*) 
     &      "Internal Error 6 in SMUMPS_BUF_SEND_CONTRIB_TYPE2",
     &      " SIZE msg/available= ",SIZE_PACK, SIZE_AV + SIZE1,
     &      " IERR= ",IERR
          ENDIF
          NBROWS_PACKET = NBROWS_PACKET - 1
          IF (NBROWS_PACKET > 0 ) GOTO 10
        ENDIF
        IF ( IERR .LT. 0 ) GOTO 100
        POSITION = 0
        CALL MPI_PACK( IPERE, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI  )
        CALL MPI_PACK( ISON, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI  )
        CALL MPI_PACK( NBROW, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI  )
        IF (KEEP(50)==0) THEN
          CALL MPI_PACK( LROW, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI  )
        ELSE
          IF (CB_IS_LR.AND.
     &        NBLRB_ALREADY_SENT+NBLRB_PACKET .EQ. NBLRB_TOTAL) THEN
            CALL MPI_PACK( -MAX_ROW_LENGTH, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI  )
          ELSE
            CALL MPI_PACK( MAX_ROW_LENGTH, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI  )
          ENDIF
        ENDIF
        CALL MPI_PACK( NBROWS_ALREADY_SENT, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI  )
        CALL MPI_PACK( NBROWS_PACKET, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI  )
        CALL MPI_PACK( CB_IS_LR_INT, 1, MPI_INTEGER,
     &          BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &          POSITION, COMM, IERR_MPI  )
        IF ( PDEST .NE. PDEST_MASTER ) THEN
          IF (KEEP(50)==0) THEN
          CALL MPI_PACK( IW_CBSON( HS + NROW +  NPIV + 1 ), LROW,
     &                 MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI  )
          ELSE
           IF (MAX_ROW_LENGTH > 0) THEN
           CALL MPI_PACK( IW_CBSON( HS + NROW +  NPIV + 1 ),
     &                 MAX_ROW_LENGTH,
     &                 MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI  )
           ENDIF
          ENDIF
        END IF
        DO J=NBROWS_ALREADY_SENT+1,NBROWS_ALREADY_SENT+NBROWS_PACKET
           I = POS_FIRST_ROW_TO_PDEST + J -1
           INDICE_PERE=MAPROW(I)
           CALL MUMPS_BLOC2_GET_ISLAVE(
     &          KEEP,KEEP8, IPERE, STEP, N, SLAVEF,
     &          ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &          NASS_PERE,
     &          NFRONT_PERE - NASS_PERE,
     &          NSLAVES_PERE,
     &          INDICE_PERE,
     &          NOSLA,
     &          IPOS_IN_SLAVE )
           INDICE_PERE = IPOS_IN_SLAVE
           CALL MPI_PACK( INDICE_PERE, 1, MPI_INTEGER,
     &          BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &          POSITION, COMM, IERR_MPI )
        ENDDO
        IF (CB_IS_LR.AND.(NBROW.GT.0)) THEN
          CALL SMUMPS_BLR_PACK_CB_LRB(
     &             CB_LRB, NB_ROW_SHIFT,
     &             NBCOLS_ALREADY_SENT, 
     &             NBLRB_ALREADY_SENT, NBLRB_PACKET, 
     &             PANEL2SEND,
     &             PANEL_BEG_OFFSET+1, PANEL_BEG_OFFSET+NBROWS_PACKET,
     &             BUF_CB%CONTENT(IPOS:),
     &             SIZE_PACK, POSITION, COMM, IERR
     &             , (KEEP(452).GT.0), KEEP(34)
     &             )
          GOTO 200
        ENDIF
        DO J=NBROWS_ALREADY_SENT+1,NBROWS_ALREADY_SENT+NBROWS_PACKET
         I = POS_FIRST_ROW_TO_PDEST + J -1
         IF (KEEP(50).ne.0) THEN
            THIS_ROW_LENGTH = LROW + I - LMAP
         ELSE
            THIS_ROW_LENGTH = LROW
         ENDIF
         IF (DESC_IN_LU) THEN 
            IF ( PACKED_CB ) THEN
             IF (NELIM.EQ.0) THEN
               ITMP8 = int(I,8)
             ELSE
               ITMP8 = int(NELIM+I,8)
             ENDIF
             APOS = ITMP8 * (ITMP8-1_8) / 2_8 + 1_8
            ELSE
             APOS = int(I+NELIM-1, 8) * int(LROW,8) + 1_8
            ENDIF
         ELSE
            IF ( PACKED_CB ) THEN
             IF ( LROW .EQ. NROW )  THEN
               ITMP8 = int(I,8)
               APOS  = ITMP8 * (ITMP8-1_8)/2_8 + 1_8
             ELSE
               ITMP8 = int(I + LROW - NROW,8)
               APOS  = ITMP8 * (ITMP8-1_8)/2_8 + 1_8 -
     &                 int(LROW - NROW, 8) * int(LROW-NROW+1,8) / 2_8
             ENDIF
            ELSE
             APOS = int( I - 1, 8 ) * LDA_SON8 + SHIFTCB_SON + 1_8
            ENDIF
         ENDIF
         CALL MPI_PACK( A_CBSON( APOS ), THIS_ROW_LENGTH,
     &        MPI_REAL,
     &        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &        POSITION, COMM, IERR_MPI )
        ENDDO
 200    CONTINUE
        IF (COMPUTE_MAX) THEN
           CALL MPI_PACK(NFS4FATHER,1,
     &          MPI_INTEGER,
     &          BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &          POSITION, COMM, IERR_MPI )
           IF (NFS4FATHER .GT. 0) THEN
            IF (CB_IS_LR) THEN
              CALL SMUMPS_BLR_RETRIEVE_M_ARRAY (
     &            IW_CBSON(1+XXF), M_ARRAY)
              CALL MPI_PACK(M_ARRAY(1), NFS4FATHER,
     &             MPI_REAL,
     &             BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &             POSITION, COMM, IERR_MPI )
              CALL SMUMPS_BLR_FREE_M_ARRAY ( IW_CBSON(1+XXF) )
            ELSE
              BUF_MAX_ARRAY(1:NFS4FATHER) = ZERO
              IF(MAPROW(NROW) .GT. NASS_PERE) THEN
                 DO PS1=1,NROW
                    IF(MAPROW(PS1).GT.NASS_PERE) EXIT
                 ENDDO
                 IF (DESC_IN_LU) THEN
                   IF (PACKED_CB) THEN
                    APOS = int(NELIM+PS1,8) * int(NELIM+PS1-1,8) /
     &                     2_8 + 1_8
                    NCA  = -44444
                    ASIZE  = int(NROW,8) * int(NROW+1,8)/2_8 -
     &                       int(NELIM+PS1,8) * int(NELIM+PS1-1,8)/2_8
                    LROW1  = PS1 + NELIM
                   ELSE
                    APOS = int(PS1+NELIM-1,8) * int(LROW,8) + 1_8
                    NCA = LROW
                    ASIZE = int(NCA,8) * int(NROW-PS1+1,8)
                    LROW1 = LROW
                   ENDIF
                 ELSE
                    IF (PACKED_CB) THEN
                      IF (NPIV.NE.0) THEN
         WRITE(*,*) "Error in PARPIV/SMUMPS_BUF_SEND_CONTRIB_TYPE2"
         CALL MUMPS_ABORT()
                      ENDIF
                      LROW1=LROW-NROW+PS1
                      ITMP8 = int(PS1 + LROW - NROW,8)
                      APOS = ITMP8 * (ITMP8 - 1_8) / 2_8 + 1_8 -
     &                       int(LROW-NROW,8)*int(LROW-NROW+1,8)/2_8
                      ASIZE = int(LROW,8)*int(LROW+1,8)/2_8 -
     &                       ITMP8*(ITMP8-1_8)/2_8
                      NCA   = -555555
                    ELSE
                      APOS = int(PS1-1,8) * LDA_SON8 + 1_8 + SHIFTCB_SON
                      NCA = int(LDA_SON8)
                      ASIZE = LA_CBSON - APOS + 1_8
                      LROW1=-666666
                    ENDIF
                 ENDIF
                 IF ( NROW-PS1+1-KEEP253_LOC-NVSCHUR .GT. 0 ) THEN
                   CALL SMUMPS_COMPUTE_MAXPERCOL(
     &                A_CBSON(APOS),ASIZE,NCA,
     &                NROW-PS1+1-KEEP253_LOC-NVSCHUR,
     &                BUF_MAX_ARRAY,NFS4FATHER,PACKED_CB,LROW1)
                 ENDIF
              ENDIF
              CALL MPI_PACK(BUF_MAX_ARRAY, NFS4FATHER,
     &             MPI_REAL,
     &             BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &             POSITION, COMM, IERR_MPI )
            ENDIF 
           ENDIF
        ENDIF 
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                  PDEST, CONTRIB_TYPE2, COMM,
     &                  BUF_CB%CONTENT( IREQ ), IERR_MPI )
        IF ( SIZE_PACK.LT. POSITION ) THEN
          WRITE(*,*) ' contniv2: SIZE, POSITION =',SIZE_PACK, POSITION
          WRITE(*,*) ' NBROW, LROW = ', NBROW, LROW
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE_PACK .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_CB, POSITION )
        IF (CB_IS_LR) THEN
          IF (NBLRB_ALREADY_SENT+NBLRB_PACKET.EQ.NBLRB_TOTAL) THEN
           NBROWS_ALREADY_SENT = NBROWS_ALREADY_SENT + NBROWS_PACKET
           NBCOLS_ALREADY_SENT = 0
           NBLRB_ALREADY_SENT  = 0
          ELSE
           NBLRB_ALREADY_SENT = NBLRB_ALREADY_SENT + NBLRB_PACKET
          ENDIF
        ELSE
           NBROWS_ALREADY_SENT=NBROWS_ALREADY_SENT + NBROWS_PACKET
        ENDIF
        IF (NBROWS_ALREADY_SENT .NE. NBROW ) THEN
           IERR = -1
        ENDIF
 100    CONTINUE
        RETURN
      END SUBROUTINE SMUMPS_BUF_SEND_CONTRIB_TYPE2
        SUBROUTINE SMUMPS_BUF_SEND_BLOCFACTO( INODE, NFRONT,
     &             NCOL, NPIV, FPERE, LASTPANEL, IPIV, VAL,
     &             PDEST, NDEST, KEEP, NB_BLOC_FAC,
     &             NSLAVES_TOT, WIDTH, COMM,
     &             NELIM, NPARTSASS, CURRENT_BLR_PANEL, 
     &             LR_ACTIV_BlocFACTO, 
     &             BLR_LorU, SIZE_BEGS_BLR_COL, BEGS_BLR_COL,
     &             NBCOLS_ALREADY_SENT, NBLRB_ALREADY_SENT,
     &             IBEG_PANEL, 
     &             ICNTL, IERR )
      USE SMUMPS_LR_TYPE
      IMPLICIT NONE
        INTEGER, intent(in) :: INODE, NCOL, NPIV, 
     &                         FPERE, NFRONT, NDEST
        INTEGER, intent(in) :: IPIV( NPIV )
        REAL, intent(in) :: VAL( NFRONT, * )
        INTEGER, intent(in) :: PDEST( NDEST ) 
        INTEGER, intent(inout) :: KEEP(500)
        INTEGER, intent(in) :: NB_BLOC_FAC,
     &                         NSLAVES_TOT, COMM, WIDTH
        LOGICAL, intent(in) :: LR_ACTIV_BlocFACTO
        LOGICAL, intent(in) :: LASTPANEL
        INTEGER, intent(in) :: NELIM, NPARTSASS, CURRENT_BLR_PANEL, 
     &                         IBEG_PANEL
        TYPE (LRB_TYPE), DIMENSION(:), intent(in) :: BLR_LorU
        INTEGER, intent(in)                   :: SIZE_BEGS_BLR_COL
        INTEGER, POINTER, DIMENSION(:), intent(in) :: BEGS_BLR_COL
        INTEGER, intent(in) :: ICNTL(60)
        INTEGER, intent(inout) :: IERR
        INTEGER, INTENT(inout)::  NBCOLS_ALREADY_SENT, 
     &                            NBLRB_ALREADY_SENT
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER POSITION, IREQ, IPOS, SIZE1, SIZE3, SIZET,
     &          IDEST, IPOSMSG, I, SIZE_MSG_BYTES
        LOGICAL OVERFLOW
        INTEGER(8) :: SIZE2_8, SIZE3_8, SIZET_8, SIZE_OVERFLOW
        INTEGER NPIVSENT
        INTEGER :: LP
        LOGICAL :: LPOK
        LOGICAL :: RECV_BUF_SMALLER_THAN_SEND, NOT_ENOUGH_SPACE
        INTEGER :: DEST_BLOCFACTO, TAG_BLOCFACTO
        INTEGER :: LR_ACTIVATED_INT
        INTEGER :: NBINT, SIZE_AV, SIZE_AV_ADJUSTED
        INTEGER :: JBEG_BLOCK, NCOL_SEND, NCOL_MAX, 
     &             NCOL_DIAG, NEWCOL_SENT
        INTEGER :: NBLRB_PACKET, SIZE_NEXT_BLOCK
        LOGICAL :: AVOID_TOO_SMALL_GRANULARITY
        INTEGER, PARAMETER :: kmaxcol=3
        REAL    :: K170PER1000
        LP               = ICNTL( 1 )
        LPOK             = ( LP.GT.0 .AND. ICNTL(4).GE.1 )
        IERR             = 0
        OVERFLOW         = .FALSE.
        NOT_ENOUGH_SPACE = .FALSE.
        NBLRB_PACKET     = -9988
        NCOL_DIAG        = -9988
        AVOID_TOO_SMALL_GRANULARITY =  .TRUE.
        IF (KEEP(173).EQ.1) AVOID_TOO_SMALL_GRANULARITY = .FALSE.
        SIZE_OVERFLOW      = 0_8
        JBEG_BLOCK         = NBCOLS_ALREADY_SENT + 1
        NCOL_SEND    = NCOL - JBEG_BLOCK + 1
        NEWCOL_SENT  = NCOL_SEND
        CALL MUMPS_BUF_SIZE_AVAILABLE( BUF_CB, SIZE_AV )
        IF ( SIZE_AV .LT. SIZE_RBUF_BYTES ) THEN
          RECV_BUF_SMALLER_THAN_SEND = .FALSE.
        ELSE
          RECV_BUF_SMALLER_THAN_SEND = .TRUE.
          SIZE_AV = SIZE_RBUF_BYTES
        ENDIF
        IF (
     &         (KEEP(50).NE.0) .OR. 
     &         (KEEP(50).EQ.0.AND.JBEG_BLOCK.EQ.1)
     &     ) THEN
          NBINT = NPIV
        ELSE
          NBINT = 0
        ENDIF
        IF ( LASTPANEL ) THEN
          IF ( KEEP(50) .eq. 0 ) THEN
            NBINT = 9 + NBINT
          ELSE
            NBINT = 11 + NBINT
          END IF
        ELSE
          IF ( KEEP(50) .eq. 0 ) THEN
            NBINT = 8 + NBINT 
          ELSE
            NBINT = 10 + NBINT
          END IF
        END IF
        NBINT = NBINT + 1 
        IF (SIZE_BEGS_BLR_COL .GT.0 ) THEN
            NBINT = NBINT + SIZE_BEGS_BLR_COL
        ENDIF
        CALL MPI_PACK_SIZE( NBINT + ( NDEST - 1 ) * OVHSIZE,
     &                          MPI_INTEGER, COMM, SIZE1, IERR_MPI )
        SIZE2_8 = 0_8
        SIZE_AV_ADJUSTED  = SIZE_AV  
        SIZE_NEXT_BLOCK   = 0   
        IF ( (NPIV.GT.0)
     &     ) THEN
          SIZE_AV_ADJUSTED  = SIZE_AV_ADJUSTED - int(SIZE2_8) - SIZE1
          NOT_ENOUGH_SPACE  = (SIZE_AV_ADJUSTED.LE.0)
          IF (.NOT. LR_ACTIV_BlocFACTO)  THEN
            NCOL_MAX  = (SIZE_AV_ADJUSTED) / (NPIV*SIZEofREAL)
            NCOL_MAX  = max(NCOL_MAX,0)
            NCOL_SEND = min( NCOL_SEND, NCOL_MAX)
            NOT_ENOUGH_SPACE = NOT_ENOUGH_SPACE.OR. 
     &                (NCOL_SEND.EQ.0) .OR. 
     &                ((JBEG_BLOCK.EQ.1).AND.(NCOL_MAX.LT.NPIV))
            IF (JBEG_BLOCK.EQ.1) NCOL_SEND = max(NCOL_SEND, NPIV)
            IF (KEEP(173).EQ.1) THEN
              IF (JBEG_BLOCK.EQ.1) THEN
                NCOL_SEND = min(NCOL_SEND, kmaxcol+NPIV)
              ELSE
                NCOL_SEND = min(NCOL_SEND, kmaxcol)
              ENDIF
            ENDIF
            NOT_ENOUGH_SPACE= NOT_ENOUGH_SPACE.OR.
     &                        (NCOL_SEND .GT. NCOL_MAX)
            SIZE3_8 = int(NPIV,8)*int(NCOL_SEND,8)*int(KEEP(35),8)
            IF (SIZE3_8 .GT. int(huge(SIZE3),8))  THEN
              OVERFLOW      = .TRUE.
              SIZE_OVERFLOW = SIZE3_8
            ELSE
              CALL MPI_PACK_SIZE( NPIV*NCOL_SEND, 
     &                      MPI_REAL,
     &                      COMM, SIZE3, IERR_MPI )
              SIZE2_8 = SIZE2_8 + int(SIZE3,8)
            ENDIF
            NEWCOL_SENT = NCOL_SEND
            IF (NBCOLS_ALREADY_SENT+NEWCOL_SENT.LT.NCOL)  THEN
              CALL MPI_PACK_SIZE( NPIV, 
     &                      MPI_REAL,
     &                      COMM, SIZE_NEXT_BLOCK, IERR_MPI )
            ENDIF
          ELSE 
            NCOL_DIAG = -9995
            IF ((KEEP(50).NE.0).OR.(JBEG_BLOCK.EQ.1)) THEN
              SIZE3_8 = int(NPIV,8)*int(NPIV+NELIM,8)*int(KEEP(35),8)
              IF (SIZE3_8 .GT. int(huge(SIZE3),8)) THEN
                OVERFLOW = .TRUE.
                SIZE_OVERFLOW = SIZE3_8
              ELSE
                CALL MPI_PACK_SIZE( NPIV*(NPIV+NELIM), 
     &                        MPI_REAL,
     &                        COMM, SIZE3, IERR_MPI )
                SIZE2_8 = SIZE2_8+int(SIZE3,8)
                NCOL_SEND = NPIV+NELIM
                SIZE_AV_ADJUSTED =  SIZE_AV_ADJUSTED - int(SIZE2_8)
              ENDIF
            ELSE
              NCOL_SEND = 0
            ENDIF
            NCOL_DIAG = NCOL_SEND
            IF (JBEG_BLOCK.EQ.1) THEN
              NEWCOL_SENT = NCOL_DIAG
            ELSE
              NEWCOL_SENT = 0
            ENDIF
            NOT_ENOUGH_SPACE =   ( NOT_ENOUGH_SPACE.OR.
     &                             (SIZE_AV_ADJUSTED.LE.0) )
            CALL SMUMPS_MPI_PACK_SIZE_LR_PARTIAL ( 0,
     &          BLR_LorU, NBLRB_ALREADY_SENT, 
     &          SIZE_AV_ADJUSTED, KEEP(173),
     &          NBLRB_PACKET, NCOL_SEND, SIZE3_8, 
     &          SIZE_NEXT_BLOCK,
     &          COMM, IERR 
     &           )
             NEWCOL_SENT = NEWCOL_SENT + (NCOL_SEND-NCOL_DIAG)
             NOT_ENOUGH_SPACE= ( NOT_ENOUGH_SPACE.OR.
     &            (NEWCOL_SENT.EQ.0).OR. 
     &            (SIZE3_8.GT.int(SIZE_AV_ADJUSTED,8)) )
             IF (SIZE3_8 .GT. int(huge(SIZE3),8)) THEN
              OVERFLOW = .TRUE.
              SIZE_OVERFLOW = SIZE3_8
             ENDIF
             SIZE2_8 = SIZE2_8+SIZE3_8
          ENDIF
        ENDIF
        IF (SIZE_NEXT_BLOCK.GT.SIZE_RBUF_BYTES) THEN
         IERR = -3
         RETURN
        ENDIF
        IF (NOT_ENOUGH_SPACE
     &      ) THEN
          IF (RECV_BUF_SMALLER_THAN_SEND
     &       ) THEN
            IERR = -3
            RETURN
          ELSE
            IERR = -1
            RETURN
          ENDIF
        ENDIF
        SIZET_8 = int(SIZE1,8) + SIZE2_8
        IF (SIZET_8 .GT. int(huge(SIZE3),8)) THEN
              OVERFLOW      = .TRUE.
              SIZE_OVERFLOW = SIZET_8
        ENDIF
        IF (OVERFLOW) THEN
          IERR=-3
          IF (LPOK) WRITE(LP,*)
     & "Integer overflow message inSMUMPS_BUF_SEND_BLOCFACTO",
     & "SIZE_OVERFLOW,NPIV,NFRONT,NELIM=",
     &  SIZE_OVERFLOW, NPIV, NFRONT, NELIM
          RETURN
        ENDIF
        SIZET = int(SIZET_8)
        IF (SIZET.GT.SIZE_RBUF_BYTES) THEN
           IERR = -3 
            RETURN
        ENDIF
        K170PER1000 = real(min(KEEP(170),500))/real(1000)
        IF ( (NBCOLS_ALREADY_SENT+NEWCOL_SENT.LT.NCOL) 
     &    .AND. .NOT. RECV_BUF_SMALLER_THAN_SEND 
     &    .AND.  ( SIZET .LT.  
     &         int(real(SIZE_RBUF_BYTES)*K170PER1000) )
     &    .AND.  ( 
     &      int(SIZET,8) + int(SIZE_NEXT_BLOCK,8) .LE. 
     &      int(SIZE_RBUF_BYTES,8)  )
     &    .AND. AVOID_TOO_SMALL_GRANULARITY
     &     ) THEN
         IERR = -1 
         RETURN
        ENDIF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZET, IERR, 
     &                 NDEST , PDEST)
        IF ( IERR .LT. 0 ) THEN
          RETURN
        ENDIF
        BUF_CB%ILASTMSG = BUF_CB%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_CB%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_CB%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        SIZE_MSG_BYTES = SIZET - OVHSIZE * (NDEST-1) * KEEP(34)
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                        POSITION, COMM, IERR_MPI )
        NPIVSENT = NPIV
        IF (LASTPANEL) NPIVSENT = -NPIV
        CALL MPI_PACK( NPIVSENT, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                        POSITION, COMM, IERR_MPI )
        IF ( LASTPANEL .OR. KEEP(50).ne.0 ) THEN
          CALL MPI_PACK( FPERE, 1, MPI_INTEGER,
     &                   BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                   POSITION, COMM, IERR_MPI )
        END IF
        IF ( LASTPANEL .AND. KEEP(50) .NE. 0 ) THEN
            CALL MPI_PACK( NB_BLOC_FAC, 1, MPI_INTEGER,
     &                   BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                   POSITION, COMM, IERR_MPI )
        END IF
        CALL MPI_PACK( NEWCOL_SENT, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( JBEG_BLOCK, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NELIM, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NPARTSASS, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( CURRENT_BLR_PANEL, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                        POSITION, COMM, IERR_MPI )
        IF (LR_ACTIV_BlocFACTO) THEN
          LR_ACTIVATED_INT = 1
        ELSE
          LR_ACTIVATED_INT = 0
        ENDIF
        CALL MPI_PACK( LR_ACTIVATED_INT, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                        POSITION, COMM, IERR_MPI )
        IF ( KEEP(50) .ne. 0 ) THEN
          CALL MPI_PACK( NSLAVES_TOT, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                        POSITION, COMM, IERR_MPI )
        ENDIF
        CALL MPI_PACK( SIZE_BEGS_BLR_COL, 1, MPI_INTEGER,
     &                          BUF_CB%CONTENT( IPOSMSG ), 
     &                          SIZE_MSG_BYTES,
     &                          POSITION, COMM, IERR_MPI )
        IF ( SIZE_BEGS_BLR_COL.GT.0 ) THEN
          IF (SIZE_BEGS_BLR_COL.GT.0) THEN
              CALL MPI_PACK( BEGS_BLR_COL, SIZE_BEGS_BLR_COL, 
     &                          MPI_INTEGER,
     &                          BUF_CB%CONTENT( IPOSMSG ), 
     &                          SIZE_MSG_BYTES,
     &                          POSITION, COMM, IERR_MPI )
          ENDIF
        ENDIF
        IF ( (NPIV.GT.0)
     &     ) THEN
          IF (
     &         (KEEP(50).NE.0) .OR.
     &         (KEEP(50).EQ.0.AND.JBEG_BLOCK.EQ.1)
     &        ) THEN
             CALL MPI_PACK( IPIV, NPIV, MPI_INTEGER,
     &                      BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                      POSITION, COMM, IERR_MPI )
          ENDIF
          IF (LR_ACTIV_BlocFACTO) THEN
              IF (KEEP(50).NE.0.OR.JBEG_BLOCK.EQ.1) THEN
                DO I = 1, NPIV
                  CALL MPI_PACK( VAL(1,I), NPIV+NELIM,
     &                    MPI_REAL,
     &                    BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                    POSITION, COMM, IERR_MPI )
                END DO
              ENDIF
              CALL SMUMPS_MPI_PACK_LR_PARTIAL( 
     &        BLR_LorU, NBLRB_ALREADY_SENT, NBLRB_PACKET,
     &        BUF_CB%CONTENT(IPOSMSG:
     &             IPOSMSG+(SIZE_MSG_BYTES+KEEP(34)-1)/KEEP(34)-1),
     &        SIZE_MSG_BYTES, POSITION, COMM, IERR,KEEP(34) ) 
          ELSE
            DO I = 1, NPIV
              CALL MPI_PACK( VAL(JBEG_BLOCK,I), NCOL_SEND,
     &                        MPI_REAL,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                        POSITION, COMM, IERR_MPI )
            END DO
          ENDIF
        ENDIF
        DO IDEST = NDEST, 1, -1
          DEST_BLOCFACTO = PDEST(IDEST)
          IF ( KEEP(50) .EQ. 0) THEN
            TAG_BLOCFACTO = BLOC_FACTO
            KEEP(266)=KEEP(266)+1
            CALL MPI_ISEND( BUF_CB%CONTENT( IPOSMSG ), POSITION, 
     &                MPI_PACKED,
     &                DEST_BLOCFACTO, TAG_BLOCFACTO, COMM,
     &                BUF_CB%CONTENT( IREQ + ( IDEST-1 ) * OVHSIZE ),
     &                IERR_MPI )
          ELSE
            KEEP(266)=KEEP(266)+1
            CALL MPI_ISEND( BUF_CB%CONTENT( IPOSMSG ), POSITION, 
     &                MPI_PACKED,
     &                DEST_BLOCFACTO, BLOC_FACTO_SYM, COMM,
     &                BUF_CB%CONTENT( IREQ + ( IDEST-1 ) * OVHSIZE ),
     &                IERR_MPI )
          END IF
        END DO
        IF (NBCOLS_ALREADY_SENT+NEWCOL_SENT.EQ.NCOL
     &     ) THEN
         NBCOLS_ALREADY_SENT = 0
         NBLRB_ALREADY_SENT  = 0
        ELSE
         NBCOLS_ALREADY_SENT = NBCOLS_ALREADY_SENT + NEWCOL_SENT
         IF (LR_ACTIV_BlocFACTO) THEN
          NBLRB_ALREADY_SENT  = NBLRB_ALREADY_SENT + NBLRB_PACKET
         ENDIF
         IERR = -1
        ENDIF
        IF ( SIZE_MSG_BYTES .LT. POSITION ) THEN
          WRITE(*,*) ' Error sending blocfacto : size < position'
          WRITE(*,*) ' INODE= ', INODE, 
     &               ' Size,position= ',SIZE_MSG_BYTES,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE_MSG_BYTES .NE. POSITION ) 
     &                  CALL BUF_ADJUST( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE SMUMPS_BUF_SEND_BLOCFACTO
        SUBROUTINE SMUMPS_BUF_SEND_BLFAC_SLAVE( INODE,
     &           NPIV, FPERE, IPOSK, JPOSK, UIP21K, LUIP21K, NCOLU,
     &           NDEST, PDEST, COMM, KEEP,
     &           LR_ACTIV_BlocFACTO, 
     &           BLR_LS, IPANEL,
     &           NBCOLS_ALREADY_SENT, NBLRB_ALREADY_SENT,
     &           NOTHING_WAS_SENT,
     &           A , LA, POSBLOCFACTO, LD_BLOCFACTO,
     &           IPIV, MAXI_CLUSTER, IERR, IERROR )
      USE SMUMPS_LR_TYPE
        IMPLICIT NONE
        INTEGER INODE, NCOLU, IPOSK, JPOSK, NPIV, NDEST, FPERE
        INTEGER(8) :: LUIP21K
        REAL UIP21K( : )
        INTEGER PDEST( NDEST ) 
        INTEGER   COMM, IERR, IERROR
        INTEGER, intent(inout)  :: KEEP(500)
        LOGICAL, intent(in)     :: LR_ACTIV_BlocFACTO
        INTEGER, intent(inout)  :: NBCOLS_ALREADY_SENT, 
     &                             NBLRB_ALREADY_SENT
        LOGICAL, intent(out)  :: NOTHING_WAS_SENT
        TYPE (LRB_TYPE), DIMENSION(:), POINTER :: BLR_LS
        INTEGER(8), intent(in)  :: LA, POSBLOCFACTO
        INTEGER, intent(in)     :: LD_BLOCFACTO, IPIV(NPIV), 
     &                             MAXI_CLUSTER, IPANEL
        REAL, intent(inout)  :: A(LA)
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER LR_ACTIVATED_INT
        INTEGER POSITION, IREQ, IPOS, SIZE1, SIZET,
     &          IDEST, IPOSMSG, SSS, SIZE3, SIZE_MSG_BYTES
        INTEGER(8) :: SIZE2_8, SIZE3_8, SIZET_8, SIZE_OVERFLOW
        LOGICAL    :: OVERFLOW, LASTBL_INPANEL
        INTEGER :: JBEG_BLOCK, NCOL_SEND, NCOL_MAX
        INTEGER :: SIZE_AV, SIZE_AV_ADJUSTED
        LOGICAL :: RECV_BUF_SMALLER_THAN_SEND, NOT_ENOUGH_SPACE
        INTEGER :: NBLRB_PACKET,  SIZE_NEXT_BLOCK
        LOGICAL :: AVOID_TOO_SMALL_GRANULARITY
        INTEGER, PARAMETER :: kmaxcol=3
        REAL    :: K170PER1000
        IERR          = 0
        OVERFLOW      = .FALSE.
        SIZE_OVERFLOW = 0_8
        JBEG_BLOCK    = NBCOLS_ALREADY_SENT + 1
        NCOL_SEND     = NCOLU - JBEG_BLOCK + 1
        NBLRB_PACKET  = -9977
        NOTHING_WAS_SENT = .TRUE.
        AVOID_TOO_SMALL_GRANULARITY =  .TRUE.
        IF (KEEP(173).EQ.1) AVOID_TOO_SMALL_GRANULARITY = .FALSE.
        CALL MUMPS_BUF_SIZE_AVAILABLE( BUF_CB, SIZE_AV )
        IF ( SIZE_AV .LT. SIZE_RBUF_BYTES ) THEN
          RECV_BUF_SMALLER_THAN_SEND = .FALSE.
        ELSE
          RECV_BUF_SMALLER_THAN_SEND = .TRUE.
          SIZE_AV = SIZE_RBUF_BYTES
        ENDIF
        CALL MPI_PACK_SIZE( 8 + ( NDEST - 1 ) * OVHSIZE,
     &                      MPI_INTEGER, COMM, SIZE1, IERR_MPI )
        SIZE2_8  = 0_8
        SIZE_AV_ADJUSTED = SIZE_AV - SIZE1
        SIZE_NEXT_BLOCK  = 0   
        NOT_ENOUGH_SPACE  = (SIZE_AV_ADJUSTED.LE.0)
        IF (.NOT. LR_ACTIV_BlocFACTO) THEN
           NCOL_MAX  = (SIZE_AV_ADJUSTED) / (NPIV*SIZEofREAL)
           NCOL_MAX  = max(NCOL_MAX,0)
           NCOL_SEND = min( NCOL_SEND, NCOL_MAX)
           IF (KEEP(173).EQ.1) THEN
                NCOL_SEND = min(NCOL_SEND, kmaxcol)
           ENDIF
           NOT_ENOUGH_SPACE = NOT_ENOUGH_SPACE.OR.(NCOL_SEND.EQ.0) 
           SIZE3_8 = int(NPIV,8)*int(NCOL_SEND,8)
           IF (SIZE3_8 .GT. int(huge(SIZE3),8))  THEN
               OVERFLOW      = .TRUE.
               SIZE_OVERFLOW = SIZE3_8
           ELSE
            CALL MPI_PACK_SIZE( abs(NPIV)*NCOL_SEND, 
     &                       MPI_REAL,
     &                       COMM, SIZE3, IERR_MPI )
            SIZE2_8=SIZE2_8 + int(SIZE3,8)
           ENDIF
            IF (NBCOLS_ALREADY_SENT+NCOL_SEND.LT.NCOLU) THEN
              CALL MPI_PACK_SIZE( NPIV, 
     &                      MPI_REAL,
     &                      COMM, SIZE_NEXT_BLOCK, IERR_MPI )
            ENDIF
        ELSE
            NCOL_SEND    = 0
            NOT_ENOUGH_SPACE =   ( NOT_ENOUGH_SPACE.OR.
     &                             (SIZE_AV_ADJUSTED.LE.0) )            
            CALL SMUMPS_MPI_PACK_SIZE_LR_PARTIAL ( 1,
     &          BLR_LS, NBLRB_ALREADY_SENT, 
     &          SIZE_AV_ADJUSTED, KEEP(173),
     &          NBLRB_PACKET, NCOL_SEND, SIZE3_8, 
     &          SIZE_NEXT_BLOCK,
     &          COMM, IERR 
     &           )
             NOT_ENOUGH_SPACE= ( NOT_ENOUGH_SPACE.OR.
     &            (NCOL_SEND.EQ.0).OR.
     &            (SIZE3_8.GT.int(SIZE_AV_ADJUSTED,8)) )
             IF (SIZE3_8 .GT. int(huge(SIZE3),8)) THEN
              OVERFLOW = .TRUE.
              SIZE_OVERFLOW = SIZE3_8
             ENDIF
             SIZE2_8 = SIZE2_8+SIZE3_8
        ENDIF
        IF (SIZE_NEXT_BLOCK.GT.SIZE_RBUF_BYTES) THEN
         IERR = -3
         RETURN
        ENDIF
        IF (NOT_ENOUGH_SPACE) THEN
          IF (RECV_BUF_SMALLER_THAN_SEND) THEN
            IERR = -3
            RETURN
          ELSE
            IERR = -1
            RETURN
          ENDIF
        ENDIF
        SIZET_8 = int(SIZE1,8) + SIZE2_8
        IF (SIZET_8 .GT. int(huge(SIZE3),8)) THEN
              OVERFLOW      = .TRUE.
              SIZE_OVERFLOW = SIZET_8
        ENDIF
        IF (OVERFLOW) THEN
          IERR=-3
          RETURN
         ENDIF
        SIZET = int(SIZET_8)
        IF (SIZET.GT.SIZE_RBUF_BYTES) THEN
         CALL MPI_PACK_SIZE( 6 ,
     &                      MPI_INTEGER, COMM, SSS, IERR_MPI )
         SIZE2_8 = int(SSS,8)+SIZE2_8
         IF (int(SIZE2_8).GT.SIZE_RBUF_BYTES) THEN
           IERR = -2
           RETURN
         ENDIF
        END IF
        K170PER1000 = real(min(KEEP(170),500))/real(1000)
        IF ((NBCOLS_ALREADY_SENT+NCOL_SEND.LT.NCOLU)
     &    .AND. .NOT. RECV_BUF_SMALLER_THAN_SEND 
     &    .AND. ( SIZET .LT.  
     &         int(real(SIZE_RBUF_BYTES)*K170PER1000) )
     &    .AND.  ( 
     &      int(SIZET,8) + int(SIZE_NEXT_BLOCK,8) .LE. 
     &      int(SIZE_RBUF_BYTES,8) )
     &    .AND. AVOID_TOO_SMALL_GRANULARITY
     &     ) THEN
         IERR = -1 
         RETURN
        ENDIF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZET, IERR, 
     &                 NDEST, PDEST)
        IF ( IERR .LT. 0 ) THEN
      RETURN
        ENDIF
        BUF_CB%ILASTMSG = BUF_CB%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_CB%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_CB%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        SIZE_MSG_BYTES = SIZET - OVHSIZE * (NDEST-1) * KEEP(34)
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &                  BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                  POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( IPOSK, 1, MPI_INTEGER,
     &                  BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                  POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( JPOSK+JBEG_BLOCK-1, 1, MPI_INTEGER,
     &                  BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                  POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NPIV, 1, MPI_INTEGER,
     &                  BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                  POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( FPERE, 1, MPI_INTEGER,
     &                  BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                  POSITION, COMM, IERR_MPI )
        LASTBL_INPANEL = (NBCOLS_ALREADY_SENT+NCOL_SEND.EQ.NCOLU)
        IF (LASTBL_INPANEL) THEN
          CALL MPI_PACK( -NCOL_SEND, 1, MPI_INTEGER,
     &                  BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                  POSITION, COMM, IERR_MPI )
        ELSE
          CALL MPI_PACK( NCOL_SEND, 1, MPI_INTEGER,
     &                  BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                  POSITION, COMM, IERR_MPI )
        ENDIF
        IF (LR_ACTIV_BlocFACTO) THEN
          LR_ACTIVATED_INT = 1
        ELSE
          LR_ACTIVATED_INT = 0
        ENDIF
        CALL MPI_PACK( LR_ACTIVATED_INT, 1, MPI_INTEGER,
     &                  BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                  POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( IPANEL, 1, MPI_INTEGER,
     &                  BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                  POSITION, COMM, IERR_MPI )
        IF (LR_ACTIV_BlocFACTO) THEN
           CALL SMUMPS_MPI_PACKSCALE_LR_PARTIAL( BLR_LS,
     &           NBLRB_ALREADY_SENT, NBLRB_PACKET,
     &           BUF_CB%CONTENT( IPOSMSG:
     &            IPOSMSG+(SIZE_MSG_BYTES+KEEP(34)-1)/KEEP(34)-1 ),
     &           SIZE_MSG_BYTES, POSITION, COMM,
     &           A, LA, POSBLOCFACTO, LD_BLOCFACTO, 
     &           IPIV, NPIV, MAXI_CLUSTER, IERR, IERROR )
            IF (IERR.LT.0) RETURN
        ELSE
        CALL MPI_PACK( UIP21K(1_8+int(JBEG_BLOCK-1,8)*int(NPIV,8)),
     &                  NPIV * NCOL_SEND,
     &                  MPI_REAL,
     &                  BUF_CB%CONTENT( IPOSMSG ), SIZE_MSG_BYTES,
     &                  POSITION, COMM, IERR_MPI )
        ENDIF
        NOTHING_WAS_SENT = .FALSE.
        DO IDEST = 1, NDEST
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOSMSG ), POSITION, MPI_PACKED,
     &                  PDEST(IDEST), BLOC_FACTO_SYM_SLAVE, COMM,
     &                  BUF_CB%CONTENT( IREQ + ( IDEST-1 ) * OVHSIZE ),
     &                  IERR_MPI )
        END DO
        IF ( LASTBL_INPANEL ) THEN
         NBCOLS_ALREADY_SENT = 0
         NBLRB_ALREADY_SENT  = 0
        ELSE
         NBCOLS_ALREADY_SENT = NBCOLS_ALREADY_SENT + NCOL_SEND
         IF (LR_ACTIV_BlocFACTO) THEN
           NBLRB_ALREADY_SENT  = NBLRB_ALREADY_SENT + NBLRB_PACKET
         ENDIF
         IERR = -1
        ENDIF
        IF ( SIZE_MSG_BYTES .LT. POSITION ) THEN
          WRITE(*,*) ' Error sending blfac slave : size < position'
          WRITE(*,*) ' Size,position=',SIZE_MSG_BYTES,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE_MSG_BYTES .NE. POSITION ) 
     &                   CALL BUF_ADJUST( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE SMUMPS_BUF_SEND_BLFAC_SLAVE
        SUBROUTINE SMUMPS_BUF_SEND_CONTRIB_TYPE3( N, ISON,
     &             NBCOL_SON, NBROW_SON, INDCOL_SON, INDROW_SON,
     &             LD_SON, VAL_SON, TAG, SUBSET_ROW, SUBSET_COL,
     &             NSUBSET_ROW, NSUBSET_COL,
     &             NSUPROW, NSUPCOL,
     &             NPROW, NPCOL, MBLOCK, RG2L,
     &             NBLOCK, PDEST, COMM, IERR , 
     &             TAB, TABSIZE, TRANSP, SIZE_PACK,
     &             N_ALREADY_SENT, KEEP, BBPCBP,
     &             NELIM_ROOT, NELIM_ROW, NELIM_COL ) 
        IMPLICIT NONE
        INTEGER N, ISON, NBCOL_SON, NBROW_SON, NSUBSET_ROW, NSUBSET_COL
        INTEGER NPROW, NPCOL, MBLOCK, NBLOCK, LD_SON
        INTEGER BBPCBP
        INTEGER PDEST, TAG, COMM, IERR
        INTEGER INDCOL_SON( NBCOL_SON ), INDROW_SON( NBROW_SON )
        INTEGER SUBSET_ROW( NSUBSET_ROW ), SUBSET_COL( NSUBSET_COL )
        INTEGER, INTENT(IN) :: RG2L(N)
        INTEGER NSUPROW, NSUPCOL
        INTEGER(8), INTENT(IN) :: TABSIZE
        INTEGER SIZE_PACK
        INTEGER KEEP(500)
        REAL VAL_SON( LD_SON, * ), TAB(*)
        LOGICAL TRANSP
        INTEGER N_ALREADY_SENT
        INTEGER :: NELIM_ROOT, NELIM_ROW, NELIM_COL
        INCLUDE 'mpif.h'
        INTEGER :: IERR_MPI
        INTEGER SIZE1, SIZE2, SIZE_AV, POSITION
        INTEGER SIZE_CBP, SIZE_TMP
        INTEGER IREQ, IPOS, ITAB
        INTEGER ISUB, JSUB, I, J 
        INTEGER ILOC_ROOT, JLOC_ROOT
        INTEGER IPOS_ROOT, JPOS_ROOT
        INTEGER IONE
        LOGICAL RECV_BUF_SMALLER_THAN_SEND
        INTEGER PDEST2(1)
        PARAMETER ( IONE=1 )
        INTEGER N_PACKET
        INTEGER NSUBSET_ROW_EFF, NSUBSET_COL_EFF, NSUPCOL_EFF
        PDEST2(1) = PDEST
        IERR = 0
        IF ( NSUBSET_ROW * NSUBSET_COL .NE. 0 ) THEN
          CALL MUMPS_BUF_SIZE_AVAILABLE( BUF_CB, SIZE_AV )
          IF (SIZE_AV .LT. SIZE_RBUF_BYTES) THEN
            RECV_BUF_SMALLER_THAN_SEND = .FALSE.
          ELSE
            RECV_BUF_SMALLER_THAN_SEND = .TRUE.
            SIZE_AV = SIZE_RBUF_BYTES
          ENDIF
          SIZE_AV = min(SIZE_AV, SIZE_RBUF_BYTES)
          CALL MPI_PACK_SIZE(8 + NSUBSET_COL,
     &                      MPI_INTEGER, COMM, SIZE1, IERR_MPI )
          SIZE_CBP = 0
          IF (N_ALREADY_SENT .EQ. 0 .AND.
     &        min(NSUPROW,NSUPCOL) .GT.0) THEN
            CALL MPI_PACK_SIZE(NSUPROW, MPI_INTEGER, COMM,
     &           SIZE_CBP, IERR_MPI )
            CALL MPI_PACK_SIZE(NSUPCOL, MPI_INTEGER, COMM,
     &           SIZE_TMP, IERR_MPI )
            SIZE_CBP = SIZE_CBP + SIZE_TMP
            CALL MPI_PACK_SIZE(NSUPROW*NSUPCOL,
     &           MPI_REAL, COMM,
     &           SIZE_TMP, IERR_MPI )
            SIZE_CBP = SIZE_CBP + SIZE_TMP
            SIZE1 = SIZE1 + SIZE_CBP
          ENDIF
          IF (BBPCBP.EQ.1) THEN
            NSUBSET_COL_EFF = NSUBSET_COL - NSUPCOL
            NSUPCOL_EFF = 0
          ELSE
            NSUBSET_COL_EFF = NSUBSET_COL
            NSUPCOL_EFF = NSUPCOL
          ENDIF
          NSUBSET_ROW_EFF = NSUBSET_ROW - NSUPROW
          N_PACKET =
     &    (SIZE_AV - SIZE1) / (SIZEofINT + NSUBSET_COL_EFF * SIZEofREAL)
 10       CONTINUE
          N_PACKET = min( N_PACKET,
     &                    NSUBSET_ROW_EFF-N_ALREADY_SENT )
          IF (N_PACKET .LE. 0 .AND.
     &        NSUBSET_ROW_EFF-N_ALREADY_SENT.GT.0) THEN
             IF (RECV_BUF_SMALLER_THAN_SEND) THEN
              IERR=-3
              GOTO 100
           ELSE
              IERR = -1
              GOTO 100
            ENDIF
          ENDIF
          CALL MPI_PACK_SIZE( 8 + NSUBSET_COL_EFF + N_PACKET,
     &                      MPI_INTEGER, COMM, SIZE1, IERR_MPI )
          SIZE1 = SIZE1 + SIZE_CBP
          CALL MPI_PACK_SIZE( N_PACKET * NSUBSET_COL_EFF,
     &                      MPI_REAL,
     &                      COMM, SIZE2, IERR_MPI )
          SIZE_PACK = SIZE1 + SIZE2
          IF (SIZE_PACK .GT. SIZE_AV) THEN
            N_PACKET = N_PACKET - 1
            IF ( N_PACKET > 0 ) THEN
              GOTO 10
            ELSE
               IF (RECV_BUF_SMALLER_THAN_SEND) THEN
                IERR = -3
                GOTO 100
             ELSE
                IERR = -1
                GOTO 100
              ENDIF
            ENDIF
          ENDIF
          IF (N_PACKET + N_ALREADY_SENT .NE. NSUBSET_ROW - NSUPROW
     &         .AND.
     &         SIZE_PACK .LT. SIZE_RBUF_BYTES / 10
     &         .AND. .NOT. RECV_BUF_SMALLER_THAN_SEND)
     &         THEN
             IERR = -1
             GOTO 100
          ENDIF
        ELSE 
          N_PACKET = 0
          CALL MPI_PACK_SIZE(8,MPI_INTEGER, COMM, SIZE_PACK, IERR_MPI )
        END IF
        IF ( SIZE_PACK.GT.SIZE_RBUF_BYTES ) THEN
           IERR = -3
           GOTO 100
        ENDIF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE_PACK, IERR, 
     &                 IONE, PDEST2
     &               )
        IF ( IERR .LT. 0 ) GOTO 100
        POSITION = 0
        CALL MPI_PACK( ISON, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NSUBSET_ROW, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NSUPROW, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NSUBSET_COL, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NSUPCOL, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( N_ALREADY_SENT, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( N_PACKET, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( BBPCBP, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR_MPI )
        IF ( NSUBSET_ROW * NSUBSET_COL .NE. 0 ) THEN
          IF (N_ALREADY_SENT .EQ. 0 .AND.
     &          min(NSUPROW, NSUPCOL) .GT. 0) THEN
            DO ISUB = NSUBSET_ROW-NSUPROW+1, NSUBSET_ROW
              I =  SUBSET_ROW( ISUB )
              IF ( I .LE. NELIM_COL ) THEN
                IPOS_ROOT = NELIM_ROOT + I - 1
              ELSE
                IPOS_ROOT = RG2L(INDCOL_SON( I ))
              ENDIF
              ILOC_ROOT = MBLOCK
     &                 * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     &                 + mod( IPOS_ROOT - 1, MBLOCK ) + 1
              CALL MPI_PACK( ILOC_ROOT, 1, MPI_INTEGER,
     &                      BUF_CB%CONTENT( IPOS ),
     &                      SIZE_PACK, POSITION, COMM, IERR_MPI )
            ENDDO
            DO ISUB = NSUBSET_COL-NSUPCOL+1, NSUBSET_COL
               J = SUBSET_COL( ISUB )
               JPOS_ROOT = INDROW_SON( J ) - N
               JLOC_ROOT = NBLOCK
     &                  * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                  + mod( JPOS_ROOT - 1, NBLOCK ) + 1
              CALL MPI_PACK( JLOC_ROOT, 1, MPI_INTEGER,
     &                       BUF_CB%CONTENT( IPOS ),
     &                       SIZE_PACK, POSITION, COMM, IERR_MPI )
            ENDDO
            IF ( TABSIZE.GE.int(NSUPROW,8)*int(NSUPCOL,8) ) THEN
              ITAB = 1
              DO JSUB = NSUBSET_ROW - NSUPROW+1, NSUBSET_ROW
                J = SUBSET_ROW(JSUB)
                DO ISUB = NSUBSET_COL - NSUPCOL+1, NSUBSET_COL
                  I = SUBSET_COL(ISUB)
                  TAB(ITAB) = VAL_SON(J, I)
                  ITAB = ITAB + 1
                ENDDO
              ENDDO
              CALL MPI_PACK(TAB(1), NSUPROW*NSUPCOL,
     &         MPI_REAL, 
     &         BUF_CB%CONTENT( IPOS ),
     &         SIZE_PACK, POSITION, COMM, IERR_MPI )
            ELSE
              DO JSUB = NSUBSET_ROW - NSUPROW+1, NSUBSET_ROW
                J = SUBSET_ROW(JSUB)
                DO ISUB = NSUBSET_COL - NSUPCOL+1, NSUBSET_COL
                  I = SUBSET_COL(ISUB)
                  CALL MPI_PACK(VAL_SON(J,I), 1,
     &            MPI_REAL, 
     &            BUF_CB%CONTENT( IPOS ),
     &            SIZE_PACK, POSITION, COMM, IERR_MPI )
                ENDDO
              ENDDO
            ENDIF
          ENDIF
          IF ( .NOT. TRANSP ) THEN
            DO ISUB = N_ALREADY_SENT+1, N_ALREADY_SENT+N_PACKET
              I         = SUBSET_ROW( ISUB )
              IF ( I .LE. NELIM_ROW ) THEN
                IPOS_ROOT = NELIM_ROOT +  I - 1
              ELSE
                IPOS_ROOT = RG2L( INDROW_SON( I ) )
              ENDIF
              ILOC_ROOT = MBLOCK
     &                 * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     &                 + mod( IPOS_ROOT - 1, MBLOCK ) + 1
              CALL MPI_PACK( ILOC_ROOT, 1, MPI_INTEGER,
     &                      BUF_CB%CONTENT( IPOS ),
     &                      SIZE_PACK, POSITION, COMM, IERR_MPI )
            END DO
            DO JSUB = 1, NSUBSET_COL_EFF - NSUPCOL_EFF
              J         = SUBSET_COL( JSUB )
              IF ( J .LE. NELIM_COL ) THEN
                JPOS_ROOT = NELIM_ROOT + J - 1
              ELSE
                JPOS_ROOT = RG2L( INDCOL_SON( J ) )
              ENDIF
              JLOC_ROOT = NBLOCK
     &                  * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                  + mod( JPOS_ROOT - 1, NBLOCK ) + 1
              CALL MPI_PACK( JLOC_ROOT, 1, MPI_INTEGER,
     &                       BUF_CB%CONTENT( IPOS ),
     &                       SIZE_PACK, POSITION, COMM, IERR_MPI )
            END DO
            DO JSUB = NSUBSET_COL_EFF-NSUPCOL_EFF+1, NSUBSET_COL_EFF
               J = SUBSET_COL( JSUB )
               JPOS_ROOT = INDCOL_SON( J ) - N
               JLOC_ROOT = NBLOCK
     &                  * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                  + mod( JPOS_ROOT - 1, NBLOCK ) + 1
              CALL MPI_PACK( JLOC_ROOT, 1, MPI_INTEGER,
     &                       BUF_CB%CONTENT( IPOS ),
     &                       SIZE_PACK, POSITION, COMM, IERR_MPI )
            ENDDO
          ELSE
            DO JSUB = N_ALREADY_SENT+1, N_ALREADY_SENT+N_PACKET
              J         = SUBSET_ROW( JSUB )
              IF ( J .LE. NELIM_COL ) THEN
                IPOS_ROOT = NELIM_ROOT + J - 1
              ELSE
                IPOS_ROOT = RG2L( INDCOL_SON( J ) )
              ENDIF
              ILOC_ROOT = MBLOCK
     &                 * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     &                 + mod( IPOS_ROOT - 1, MBLOCK ) + 1
              CALL MPI_PACK( ILOC_ROOT, 1, MPI_INTEGER,
     &                       BUF_CB%CONTENT( IPOS ),
     &                       SIZE_PACK, POSITION, COMM, IERR_MPI )
            END DO
            DO ISUB = 1, NSUBSET_COL_EFF - NSUPCOL_EFF
              I         = SUBSET_COL( ISUB )
              IF ( I .LE. NELIM_ROW ) THEN 
                JPOS_ROOT = NELIM_ROOT + I - 1
              ELSE
                JPOS_ROOT = RG2L( INDROW_SON( I ) )
              ENDIF
              JLOC_ROOT = NBLOCK
     &                  * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                  + mod( JPOS_ROOT - 1, NBLOCK ) + 1
              CALL MPI_PACK( JLOC_ROOT, 1, MPI_INTEGER,
     &                      BUF_CB%CONTENT( IPOS ),
     &                      SIZE_PACK, POSITION, COMM, IERR_MPI )
            END DO
            DO ISUB = NSUBSET_COL_EFF - NSUPCOL_EFF + 1, NSUBSET_COL_EFF
              I         = SUBSET_COL( ISUB )
              JPOS_ROOT = INDROW_SON(I) - N
              JLOC_ROOT = NBLOCK
     &                  * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                  + mod( JPOS_ROOT - 1, NBLOCK ) + 1
              CALL MPI_PACK( JLOC_ROOT, 1, MPI_INTEGER,
     &                      BUF_CB%CONTENT( IPOS ),
     &                      SIZE_PACK, POSITION, COMM, IERR_MPI )
            ENDDO
          END IF
          IF ( TABSIZE.GE.int(N_PACKET,8)*int(NSUBSET_COL_EFF,8) ) THEN
            IF ( .NOT. TRANSP ) THEN
              ITAB = 1
              DO ISUB = N_ALREADY_SENT+1,
     &                  N_ALREADY_SENT+N_PACKET
                I         = SUBSET_ROW( ISUB )
                DO JSUB = 1, NSUBSET_COL_EFF
                  J              = SUBSET_COL( JSUB )
                  TAB( ITAB )    = VAL_SON(J,I)
                  ITAB           = ITAB + 1
                END DO
              END DO
              CALL MPI_PACK(TAB(1), NSUBSET_COL_EFF*N_PACKET,
     &         MPI_REAL, 
     &         BUF_CB%CONTENT( IPOS ),
     &         SIZE_PACK, POSITION, COMM, IERR_MPI )
            ELSE
              ITAB = 1
              DO JSUB = N_ALREADY_SENT+1, N_ALREADY_SENT+N_PACKET
                J = SUBSET_ROW( JSUB )
                DO ISUB = 1, NSUBSET_COL_EFF
                  I         = SUBSET_COL( ISUB )
                  TAB( ITAB ) = VAL_SON( J, I )
                  ITAB = ITAB + 1
                END DO
              END DO
              CALL MPI_PACK(TAB(1), NSUBSET_COL_EFF*N_PACKET,
     &         MPI_REAL, 
     &         BUF_CB%CONTENT( IPOS ),
     &         SIZE_PACK, POSITION, COMM, IERR_MPI )
            END IF
          ELSE
            IF ( .NOT. TRANSP ) THEN
              DO ISUB = N_ALREADY_SENT+1, N_ALREADY_SENT+N_PACKET
                I         = SUBSET_ROW( ISUB )
                DO JSUB = 1, NSUBSET_COL_EFF
                  J         = SUBSET_COL( JSUB )
                  CALL MPI_PACK( VAL_SON( J, I ), 1,
     &            MPI_REAL,
     &            BUF_CB%CONTENT( IPOS ),
     &            SIZE_PACK, POSITION, COMM, IERR_MPI )
                END DO
              END DO
            ELSE
              DO JSUB = N_ALREADY_SENT+1, N_ALREADY_SENT+N_PACKET
                J = SUBSET_ROW( JSUB )
                DO ISUB = 1, NSUBSET_COL_EFF
                  I         = SUBSET_COL( ISUB )
                  CALL MPI_PACK( VAL_SON( J, I ), 1,
     &            MPI_REAL,
     &            BUF_CB%CONTENT( IPOS ),
     &            SIZE_PACK, POSITION, COMM, IERR_MPI )
                END DO
              END DO
            END IF
          ENDIF
        END IF
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                PDEST, TAG, COMM, BUF_CB%CONTENT( IREQ ),
     &                IERR_MPI )
        IF ( SIZE_PACK .LT. POSITION ) THEN
          WRITE(*,*) ' Error sending contribution to root:Size<positn'
          WRITE(*,*) ' Size,position=',SIZE_PACK,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE_PACK .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_CB, POSITION )
        N_ALREADY_SENT = N_ALREADY_SENT + N_PACKET
        IF (NSUBSET_ROW * NSUBSET_COL .NE. 0) THEN
          IF ( N_ALREADY_SENT.NE.NSUBSET_ROW_EFF ) IERR = -1
        ENDIF
  100   CONTINUE
        RETURN
        END SUBROUTINE SMUMPS_BUF_SEND_CONTRIB_TYPE3
        SUBROUTINE SMUMPS_BUF_SEND_BACKVEC
     &             ( NRHS, INODE, W, LW, LD_W, DEST, MSGTAG,
     &               JBDEB, JBFIN, KEEP, COMM, IERR )
        IMPLICIT NONE
        INTEGER NRHS, INODE,LW,COMM,IERR,DEST,MSGTAG, LD_W
        INTEGER, intent(in) :: JBDEB, JBFIN
        REAL :: W(LD_W, *)
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INCLUDE 'mpif.h'
        INTEGER :: IERR_MPI
        INTEGER SIZE, SIZE1, SIZE2
        INTEGER POSITION, IREQ, IPOS
        INTEGER IONE, K
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        IERR = 0
        DEST2(1) = DEST
        CALL MPI_PACK_SIZE( 4 , MPI_INTEGER, COMM, SIZE1, IERR_MPI )
        CALL MPI_PACK_SIZE( LW*NRHS, MPI_REAL, COMM,
     &                      SIZE2, IERR_MPI )
        SIZE = SIZE1 + SIZE2
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE, DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( LW   , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( JBDEB   , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( JBFIN   , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        DO K=1, NRHS
        CALL MPI_PACK( W(1,K), LW, MPI_REAL,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        END DO
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                  DEST, MSGTAG, COMM,
     &                  BUF_CB%CONTENT( IREQ ), IERR_MPI )
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) 'Try_update: SIZE, POSITION = ',
     &               SIZE, POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION ) CALL BUF_ADJUST( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE SMUMPS_BUF_SEND_BACKVEC
      SUBROUTINE SMUMPS_BLR_GET_SIZEREALS_CB_LRB(
     &        SIZE_AV, CB_LRB, 
     &        NB_ROW_SHIFT, PANEL2SEND,
     &        NBLRB_ALREADY_SENT, NBLRB_TOTAL,
     &        NBLRB_PACKET, SIZE_REALS, SIZE_NEXT_BLOCK
     &        , SIZE_LRMIXED 
     &        , LRMIXED_ACTIVATED
     &        , KEEP173
     &        )
        USE SMUMPS_LR_TYPE
        IMPLICIT NONE
      TYPE(LRB_TYPE), POINTER :: CB_LRB(:,:)
      INTEGER, INTENT(IN) :: NB_ROW_SHIFT, 
     &                       PANEL2SEND, 
     &                       NBLRB_ALREADY_SENT, NBLRB_TOTAL
      INTEGER, INTENT(in) :: SIZE_AV
      INTEGER, INTENT(out) :: NBLRB_PACKET
      INTEGER, INTENT(out) :: SIZE_REALS, SIZE_NEXT_BLOCK
      LOGICAL, INTENT(IN) :: LRMIXED_ACTIVATED
      INTEGER, INTENT(out) :: SIZE_LRMIXED
      INTEGER, INTENT(IN) :: KEEP173
      INTEGER :: J
      TYPE(LRB_TYPE), POINTER :: LRB
      LOGICAL    :: NOT_ENOUGH_SPACE
      INTEGER(8) :: SIZE_BYTE8
      INTEGER    :: SIZE_REALS_TRY, SIZE_NEXT_TMP
      INTEGER :: SIZE_LRMIXED_TRY, PREC
      INTEGER, PARAMETER :: kmax=2
          SIZE_REALS    = 0   
          SIZE_REALS_TRY= 0   
          SIZE_BYTE8    = 0_8 
          NBLRB_PACKET  = 0   
          SIZE_NEXT_BLOCK = 0 
          SIZE_LRMIXED     = 0 
          SIZE_LRMIXED_TRY = 0 
          DO J=NBLRB_ALREADY_SENT+1, NBLRB_TOTAL
            LRB => CB_LRB(PANEL2SEND-NB_ROW_SHIFT,J)
            SIZE_REALS_TRY = 0
            SIZE_LRMIXED_TRY = 0
            IF (LRB%ISLR) THEN
              IF (LRB%K.GT.0) THEN
                IF (LRMIXED_ACTIVATED) THEN
                  SIZE_LRMIXED_TRY = SIZE_LRMIXED_TRY + NPREC_STOC + 2
                  DO PREC=1,NPREC_STOC
                    IF (LRB%K_STOC(PREC).NE.0) THEN
                      SIZE_LRMIXED_TRY = SIZE_LRMIXED_TRY + 2
                      SIZE_LRMIXED_TRY = SIZE_LRMIXED_TRY
     &                           + size(LRB%R_STOC(PREC)%val)
                      SIZE_LRMIXED_TRY = SIZE_LRMIXED_TRY 
     &                           +  size(LRB%Q_STOC(PREC)%val)
                    ENDIF
                  ENDDO
                ELSE
                  SIZE_REALS_TRY = LRB%K*(LRB%M+LRB%N)
                ENDIF
              ENDIF
            ELSE
              SIZE_REALS_TRY =  LRB%M*LRB%N
            ENDIF
            IF (LRMIXED_ACTIVATED) THEN
              SIZE_NEXT_TMP =  SIZE_LRMIXED_TRY*SIZEofINT
     &        + SIZE_REALS_TRY*SIZEofREAL
              SIZE_BYTE8 =  SIZE_BYTE8 + int(SIZE_NEXT_TMP,8)
            ELSE
              SIZE_NEXT_TMP = SIZE_REALS_TRY*SIZEofREAL
              SIZE_BYTE8    = SIZE_BYTE8 + int(SIZE_NEXT_TMP,8)
            ENDIF
            NOT_ENOUGH_SPACE = (int(SIZE_AV,8).LT.SIZE_BYTE8)
            IF (.NOT.NOT_ENOUGH_SPACE .AND. (KEEP173.EQ.1)) THEN
              IF (J.GE.NBLRB_ALREADY_SENT+1+kmax) 
     &                        NOT_ENOUGH_SPACE=.TRUE.
            ENDIF
            IF (NOT_ENOUGH_SPACE) THEN
               SIZE_NEXT_BLOCK = SIZE_NEXT_TMP
              EXIT
            ELSE
              NBLRB_PACKET = NBLRB_PACKET + 1
              SIZE_REALS   = SIZE_REALS + SIZE_REALS_TRY
              IF (LRMIXED_ACTIVATED) THEN
                SIZE_LRMIXED = SIZE_LRMIXED + SIZE_LRMIXED_TRY
              ENDIF
            ENDIF
          ENDDO
          RETURN
      END SUBROUTINE SMUMPS_BLR_GET_SIZEREALS_CB_LRB
      SUBROUTINE SMUMPS_BLR_PACK_CB_LRB(
     &        CB_LRB, NB_ROW_SHIFT, 
     &        NBCOLS_ALREADY_SENT, 
     &        NBLRB_ALREADY_SENT, NBLRB_PACKET, 
     &        PANEL2SEND, IROW_BEG, IROW_END,
     &        BUF, LBUF, POSITION, COMM, IERR
     &        , LRMIXED_ACTIVATED, K34
     &        )
        USE SMUMPS_LR_TYPE
        IMPLICIT NONE
      TYPE(LRB_TYPE), POINTER :: CB_LRB(:,:)
      INTEGER, INTENT(IN) :: NB_ROW_SHIFT, 
     &                       PANEL2SEND, IROW_BEG, IROW_END, 
     &                       NBLRB_ALREADY_SENT, NBLRB_PACKET
      INTEGER, INTENT(INOUT) :: NBCOLS_ALREADY_SENT
      INTEGER, intent(out) :: IERR
      INTEGER, intent(in)  :: COMM, LBUF  
      INTEGER, intent(inout) :: POSITION
      INTEGER, intent(inout) :: BUF(:) 
      INTEGER, INTENT(IN) :: K34
      LOGICAL, intent(in)  :: LRMIXED_ACTIVATED
      INTEGER :: J, IERR_MPI
      INTEGER :: MAXI_CLUSTER
      INCLUDE 'mpif.h'
          IERR = 0
          CALL MPI_PACK( NBLRB_PACKET, 
     &              1, MPI_INTEGER,
     &              BUF(1), LBUF, POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( NBCOLS_ALREADY_SENT, 
     &              1, MPI_INTEGER,
     &              BUF(1), LBUF, POSITION, COMM, IERR_MPI )
          MAXI_CLUSTER = 1
          DO J=NBLRB_ALREADY_SENT+1, NBLRB_ALREADY_SENT+NBLRB_PACKET
            MAXI_CLUSTER = max(MAXI_CLUSTER, 
     &                         CB_LRB(PANEL2SEND-NB_ROW_SHIFT,J)%N)
          ENDDO
          CALL MPI_PACK( MAXI_CLUSTER, 1, MPI_INTEGER,
     &              BUF(1), LBUF, POSITION, COMM, IERR_MPI )
          DO J=NBLRB_ALREADY_SENT+1, NBLRB_ALREADY_SENT+NBLRB_PACKET
            CALL SMUMPS_MPI_PACK_LRB(
     &               CB_LRB(PANEL2SEND-NB_ROW_SHIFT,J),
     &               IROW_BEG, IROW_END,       
     &               BUF, LBUF, POSITION, COMM, IERR
     &               , LRMIXED_ACTIVATED, K34
     &               )
            NBCOLS_ALREADY_SENT = NBCOLS_ALREADY_SENT+ 
     &                           CB_LRB(PANEL2SEND-NB_ROW_SHIFT,J)%N 
          ENDDO
      END SUBROUTINE SMUMPS_BLR_PACK_CB_LRB
      SUBROUTINE SMUMPS_MPI_PACK_SIZE_LR_PARTIAL
     &            ( CALLED_FROM, 
     &              BLR_LorU, NBLRB_ALREADY_SENT, SIZE_AV,
     &              KEEP173, NBLRB_PACKET, NCOL_SEND, SIZE_OUT8,
     &              SIZE_NEXT_BLOCK,
     &              COMM, IERR )
      USE SMUMPS_LR_TYPE
      INTEGER, intent(in)      :: CALLED_FROM
      INTEGER, intent(in)      :: NBLRB_ALREADY_SENT, SIZE_AV, KEEP173
      INTEGER, intent(out)     :: IERR
      INTEGER(8), intent(out)  :: SIZE_OUT8
      INTEGER, intent(inout)   :: NCOL_SEND
      INTEGER, intent(out)     :: NBLRB_PACKET, SIZE_NEXT_BLOCK
      INTEGER, intent(in)  :: COMM
      TYPE (LRB_TYPE), DIMENSION(:), intent(in) :: BLR_LorU
      INTEGER     :: I, IERR_MPI, SIZE_LOC
      INTEGER(8)  :: SIZE_TRY8, SIZE_AV8
      INTEGER     :: kmax
       INTEGER, PARAMETER :: kmax_blfacslave=2
       INTEGER, PARAMETER :: kmax_blocfacto=2
      INCLUDE 'mpif.h'
        kmax = kmax_blocfacto 
        IF (CALLED_FROM.EQ.1) kmax = kmax_blfacslave
        IERR         = 0
        SIZE_OUT8    = 0
        SIZE_TRY8    = SIZE_OUT8
        NBLRB_PACKET = 0
        SIZE_AV8     = int(SIZE_AV,8)
        SIZE_NEXT_BLOCK = 0
        CALL MPI_PACK_SIZE( 1, MPI_INTEGER, COMM, SIZE_LOC,  IERR_MPI )
        SIZE_TRY8 = SIZE_TRY8 + int(SIZE_LOC,8)
        SIZE_OUT8    = SIZE_TRY8
        IF (SIZE_AV8.LT.SIZE_TRY8) THEN
           RETURN
        ENDIF
        IF (KEEP173.EQ.1) THEN
        ENDIF
        DO I = 1+NBLRB_ALREADY_SENT, size(BLR_LorU)
          CALL SMUMPS_MPI_PACK_SIZE_LRB(BLR_LorU(I), SIZE_LOC, COMM, 
     &                                 IERR )
          SIZE_TRY8 = SIZE_TRY8 + int(SIZE_LOC,8)
          IF (SIZE_AV8.LT.SIZE_TRY8) THEN
             SIZE_NEXT_BLOCK = SIZE_LOC
             EXIT
          ENDIF
          NBLRB_PACKET = NBLRB_PACKET +1
          NCOL_SEND    = NCOL_SEND + BLR_LorU(I)%M
          SIZE_OUT8    = SIZE_TRY8
          IF (KEEP173.EQ.1) THEN
              IF (NBLRB_PACKET.GE.kmax) THEN
                SIZE_NEXT_BLOCK = SIZE_LOC
                EXIT
              ENDIF
          ENDIF
        ENDDO
        RETURN
      END SUBROUTINE SMUMPS_MPI_PACK_SIZE_LR_PARTIAL
      SUBROUTINE SMUMPS_MPI_PACK_SIZE_LRB(LRB, SIZE_OUT, COMM, IERR )
      USE SMUMPS_LR_TYPE
      INTEGER, intent(out) :: SIZE_OUT, IERR
      INTEGER, intent(in)  :: COMM
      TYPE (LRB_TYPE), intent(in) :: LRB
      INTEGER :: SIZE_LOC, IERR_MPI
      INCLUDE 'mpif.h'
      IERR = 0
      SIZE_OUT = 0
      CALL MPI_PACK_SIZE( 4,      
     &          MPI_INTEGER, COMM, SIZE_LOC,  IERR_MPI )
      SIZE_OUT = SIZE_OUT + SIZE_LOC
      IF ( LRB%ISLR ) THEN
        IF (LRB%K .GT. 0) THEN
          CALL MPI_PACK_SIZE( LRB%M * LRB%K,
     &       MPI_REAL, COMM, SIZE_LOC,  IERR_MPI )
          SIZE_OUT = SIZE_OUT + SIZE_LOC
          CALL MPI_PACK_SIZE( LRB%K * LRB%N,
     &       MPI_REAL, COMM, SIZE_LOC,  IERR_MPI )
          SIZE_OUT = SIZE_OUT + SIZE_LOC
        ENDIF
      ELSE
        CALL MPI_PACK_SIZE( LRB%M * LRB%N,
     &       MPI_REAL, COMM, SIZE_LOC,  IERR_MPI )
        SIZE_OUT = SIZE_OUT + SIZE_LOC
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_MPI_PACK_SIZE_LRB
      SUBROUTINE SMUMPS_MPI_PACK_LR( BLR_LorU, BUF, LBUF, POSITION,
     &                              COMM, IERR, K34 )
      USE SMUMPS_LR_TYPE
      INTEGER, intent(out) :: IERR
      INTEGER, intent(in)  :: COMM, LBUF, K34 
      INTEGER, intent(inout) :: POSITION
      INTEGER, intent(inout) :: BUF(:) 
      TYPE (LRB_TYPE), DIMENSION(:), intent(in) :: BLR_LorU
      INTEGER I
      INTEGER :: IERR_MPI
      INCLUDE 'mpif.h'
      IERR = 0
      CALL MPI_PACK( size(BLR_LorU), 1, MPI_INTEGER,
     &       BUF(1), LBUF, POSITION, COMM, IERR_MPI )
      DO I = 1, size(BLR_LorU)
        CALL SMUMPS_MPI_PACK_LRB(BLR_LorU(I), 
     &                           1, BLR_LorU(I)%M,   
     &                           BUF, LBUF, POSITION, 
     &                           COMM, IERR
     &                           , .FALSE., K34
     &                           )
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_MPI_PACK_LR
      SUBROUTINE SMUMPS_MPI_PACK_LR_PARTIAL( 
     &           BLR_LorU, NBLRB_ALREADY_SENT, NBLRB_PACKET,
     &           BUF, LBUF, POSITION,
     &           COMM, IERR, K34 )
      USE SMUMPS_LR_TYPE
      INTEGER, intent(out) :: IERR
      INTEGER, intent(in)  :: NBLRB_ALREADY_SENT, NBLRB_PACKET
      INTEGER, intent(in)  :: COMM, LBUF, K34 
      INTEGER, intent(inout) :: POSITION
      INTEGER, intent(inout) :: BUF(:) 
      TYPE (LRB_TYPE), DIMENSION(:), intent(in) :: BLR_LorU
      INTEGER I, IDEB, IFIN
      INTEGER :: IERR_MPI
      INCLUDE 'mpif.h'
      IERR = 0
      CALL MPI_PACK( NBLRB_PACKET, 1, MPI_INTEGER,
     &       BUF(1), LBUF, POSITION, COMM, IERR_MPI )
      IF (NBLRB_PACKET.EQ.0) RETURN
      IDEB = 1+NBLRB_ALREADY_SENT
      IFIN = IDEB + NBLRB_PACKET-1
      DO I =  IDEB, IFIN
        CALL SMUMPS_MPI_PACK_LRB(BLR_LorU(I), 
     &                           1, BLR_LorU(I)%M,   
     &                           BUF, LBUF, POSITION, 
     &                           COMM, IERR
     &                           , .FALSE., K34
     &                           )
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_MPI_PACK_LR_PARTIAL
      SUBROUTINE SMUMPS_MPI_PACK_LRB( LRB, IROW_BEG, IROW_END, 
     &                              BUF, LBUF, POSITION,
     &                              COMM, IERR 
     &                              , LRMIXED_ACTIVATED, K34
     &                              )
      USE SMUMPS_LR_TYPE
      IMPLICIT NONE
      INTEGER, intent(out) :: IERR
      INTEGER, intent(in)  :: COMM, LBUF  
      INTEGER, intent(inout) :: POSITION
      INTEGER, intent(inout) :: BUF(:) 
      TYPE (LRB_TYPE), intent(in) :: LRB
      INTEGER, intent(in)  :: IROW_BEG, IROW_END
      LOGICAL, intent(in)  :: LRMIXED_ACTIVATED
      INTEGER, intent(in)  :: K34
      INTEGER ISLR_INT
      INTEGER :: IERR_MPI
      INTEGER :: J, NROWS_TO_PACK
      INTEGER :: PREC
      INTEGER :: Q_STOC_SIZE, R_STOC_SIZE
      INCLUDE 'mpif.h'
      IERR = 0
      NROWS_TO_PACK = IROW_END - IROW_BEG + 1
      IF (LRB%ISLR) THEN
        ISLR_INT = 1
      ELSE
        ISLR_INT = 0
      ENDIF
      CALL MPI_PACK( ISLR_INT, 1, MPI_INTEGER,
     &     BUF(1), LBUF, POSITION, COMM, IERR_MPI )
      CALL MPI_PACK( LRB%K,
     &     1, MPI_INTEGER,
     &     BUF(1), LBUF, POSITION, COMM, IERR_MPI )
      CALL MPI_PACK( NROWS_TO_PACK,
     &     1, MPI_INTEGER,
     &     BUF(1), LBUF, POSITION, COMM, IERR_MPI )
      CALL MPI_PACK( LRB%N,
     &     1, MPI_INTEGER,
     &     BUF(1), LBUF, POSITION, COMM, IERR_MPI )
      IF (LRB%ISLR) THEN
        IF (LRB%K .GT. 0) THEN
          IF (LRMIXED_ACTIVATED) THEN
            DO PREC=1,NPREC_STOC
              CALL MPI_PACK( LRB%K_STOC(PREC),
     &            1, MPI_INTEGER,
     &            BUF(1), LBUF, POSITION, COMM, IERR_MPI )
              IF (LRB%K_STOC(PREC).GT.0) THEN
                Q_STOC_SIZE = size(LRB%Q_STOC(PREC)%val)
                CALL MPI_PACK( Q_STOC_SIZE,
     &            1, MPI_INTEGER,
     &            BUF(1), LBUF, POSITION, COMM, IERR_MPI )
                CALL MUMPS_LRMIXED_PACK( LRB%Q_STOC(PREC)%val(1), 
     &            LRB%M, LRB%K_STOC(PREC), IROW_BEG, NROWS_TO_PACK,
     &            COSTS_PREC(PREC), 
     &            BUF(1), LBUF, POSITION, COMM, IERR_MPI )
                DO WHILE ( mod(POSITION,K34) .NE. 0 )
                  POSITION = POSITION + 1
                ENDDO
                R_STOC_SIZE = size(LRB%R_STOC(PREC)%val)
                CALL MPI_PACK( R_STOC_SIZE,
     &            1, MPI_INTEGER,
     &            BUF(1), LBUF, POSITION, COMM, IERR_MPI )
                CALL MPI_PACK( LRB%R_STOC(PREC)%val(1), 
     &            R_STOC_SIZE, MPI_INTEGER,
     &            BUF(1), LBUF, POSITION, COMM, IERR_MPI )
              ENDIF
            ENDDO
          ELSE
            DO J=1,LRB%K
              CALL MPI_PACK( LRB%Q(IROW_BEG,J), 
     &          NROWS_TO_PACK, MPI_REAL,
     &          BUF(1), LBUF, POSITION, COMM, IERR_MPI )
            ENDDO
            CALL MPI_PACK( LRB%R(1,1),
     &        LRB%N*LRB%K, MPI_REAL,
     &        BUF(1), LBUF, POSITION, COMM, IERR_MPI )
          ENDIF
        ENDIF
      ELSE
        DO J=1,LRB%N
          CALL MPI_PACK( LRB%Q(IROW_BEG,J), 
     &       NROWS_TO_PACK, MPI_REAL,
     &       BUF(1), LBUF, POSITION, COMM, IERR_MPI )
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_MPI_PACK_LRB
      SUBROUTINE SMUMPS_MPI_UNPACK_LRB(
     &           BUFR, LBUFR, LBUFR_BYTES, POSITION,
     &                             LRB, KEEP8,
     &                             COMM, IFLAG, IERROR
     &                             , LRMIXED_ACTIVATED, K34
     &                             )
      USE SMUMPS_LR_CORE, ONLY : ALLOC_LRB
      USE SMUMPS_LR_TYPE
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LBUFR
      INTEGER, INTENT(IN) :: LBUFR_BYTES
      INTEGER, INTENT(IN) :: BUFR(LBUFR)
      INTEGER, INTENT(INOUT) :: POSITION
      INTEGER, INTENT(IN) :: COMM
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
      TYPE (LRB_TYPE), INTENT(OUT) :: LRB
      INTEGER(8) :: KEEP8(150)
      LOGICAL, intent(in)  :: LRMIXED_ACTIVATED
      INTEGER, intent(in) :: K34
      LOGICAL :: ISLR
      INTEGER :: ISLR_INT
      INTEGER :: K, M, N
      INTEGER :: IERR_MPI
      INTEGER :: allocok
      INTEGER :: PREC, MEM
      INTEGER :: Q_STOC_SIZE, R_STOC_SIZE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &               ISLR_INT, 1, MPI_INTEGER, COMM, IERR_MPI )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &               K, 1,
     &               MPI_INTEGER, COMM, IERR_MPI )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &               M, 1,
     &               MPI_INTEGER, COMM, IERR_MPI )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &               N, 1,
     &               MPI_INTEGER, COMM, IERR_MPI )
      IF (ISLR_INT .eq. 1) THEN
        ISLR = .TRUE.
      ELSE
        ISLR = .FALSE.
      ENDIF
      IF (LRMIXED_ACTIVATED.AND.ISLR) THEN
        LRB%M = M
        LRB%N = N
        LRB%K = K
        LRB%ISLR = ISLR
      ELSE
        CALL ALLOC_LRB( LRB, K, M, N, ISLR, 
     &           IFLAG, IERROR, KEEP8 )
        IF (IFLAG.LT.0) RETURN
      ENDIF
      IF (ISLR) THEN
        IF (K .GT. 0) THEN
          IF (LRMIXED_ACTIVATED) THEN
            ALLOCATE(LRB%K_STOC(NPREC_STOC))
            ALLOCATE(LRB%Q_STOC(NPREC_STOC))
            ALLOCATE(LRB%R_STOC(NPREC_STOC))
            DO PREC=1,NPREC_STOC
              CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &               LRB%K_STOC(PREC), 1,
     &               MPI_INTEGER, COMM, IERR_MPI )
              IF (LRB%K_STOC(PREC).GT.0) THEN
                CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &               Q_STOC_SIZE, 1,
     &               MPI_INTEGER, COMM, IERR_MPI )
                ALLOCATE(LRB%Q_STOC(PREC)%val(Q_STOC_SIZE),
     &          stat=allocok)
                IF (allocok .GT. 0) THEN
                   WRITE(*,*) "Q..%VAL ALLOC FAILED",Q_STOC_SIZE
                   CALL MUMPS_ABORT()
                ENDIF
                CALL MUMPS_LRMIXED_UNPACK( 
     &               LRB%Q_STOC(PREC)%val(1), LRB%M, LRB%K_STOC(PREC),
     &               COSTS_PREC(PREC),
     &               BUFR, LBUFR_BYTES, POSITION, COMM, IERR_MPI)
                DO WHILE (mod(POSITION,K34) .NE. 0)
                  POSITION=POSITION+1
                ENDDO
                CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &               R_STOC_SIZE, 1,
     &               MPI_INTEGER, COMM, IERR_MPI )
                ALLOCATE(LRB%R_STOC(PREC)%val(R_STOC_SIZE),
     &          stat=allocok)
                IF (allocok .GT. 0) THEN
                   WRITE(*,*) "R..%VAL ALLOC FAILED",R_STOC_SIZE
                   CALL MUMPS_ABORT()
                ENDIF
                CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &               LRB%R_STOC(PREC)%val(1), R_STOC_SIZE,
     &               MPI_INTEGER, COMM, IERR_MPI )
                MEM = (R_STOC_SIZE*K34 + COSTS_PREC(1)-1)/COSTS_PREC(1)
     &              + (Q_STOC_SIZE*K34 + COSTS_PREC(1)-1)/COSTS_PREC(1)
                CALL MUMPS_DM_FAC_UPD_DYN_MEMCNTS(int(MEM,8),
     &              .TRUE., KEEP8, IFLAG, IERROR, .TRUE., .TRUE.)
                IF (IFLAG.LT.0) RETURN
              ELSE
                NULLIFY(LRB%Q_STOC(PREC)%val)
                NULLIFY(LRB%R_STOC(PREC)%val)
              ENDIF
            ENDDO
          ELSE
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                     LRB%Q(1,1), M*K, MPI_REAL,
     &                     COMM, IERR_MPI )
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                     LRB%R(1,1), N*K, MPI_REAL,
     &                     COMM, IERR_MPI )
          ENDIF
        ENDIF
      ELSE
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   LRB%Q(1,1), M*N, MPI_REAL,
     &                   COMM, IERR_MPI )
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_MPI_UNPACK_LRB
      SUBROUTINE SMUMPS_MPI_PACKSCALE_LR_PARTIAL 
     &                  ( BLR, NBLRB_ALREADY_SENT, NBLRB_PACKET,
     &                    BUF, LBUF, POSITION, COMM, 
     &                    A , LA, POSELTD, LD_DIAG,
     &                    IPIV, NPIV, MAXI_CLUSTER, 
     &                    IERR, IERROR )
      USE SMUMPS_LR_TYPE
      INTEGER, intent(out)   :: IERR
      INTEGER, intent(inout) :: IERROR
      INTEGER, intent(in)  :: NBLRB_ALREADY_SENT, NBLRB_PACKET
      INTEGER, intent(in)  :: COMM, LBUF 
      INTEGER, intent(inout) :: POSITION
      INTEGER, intent(inout) :: BUF(:)  
      TYPE  (LRB_TYPE), DIMENSION(:), intent(in) :: BLR
      INTEGER(8), intent(in)  :: LA, POSELTD
      INTEGER, intent(in)     :: LD_DIAG, NPIV 
      INTEGER, intent(in)     :: IPIV(NPIV), MAXI_CLUSTER
      REAL, intent(inout)  :: A(LA)
      INTEGER :: IERR_MPI
      INTEGER I, ISLR_INT, J, ALLOCOK, IDEB, IFIN
      REAL, ALLOCATABLE,DIMENSION(:,:) ::  SCALED
      REAL, ALLOCATABLE,DIMENSION(:) ::  BLOCK
      REAL :: PIV1, PIV2, OFFDIAG
      INCLUDE 'mpif.h'
      IERR   = 0
      IERROR = 0
      CALL MPI_PACK( NBLRB_PACKET, 1, MPI_INTEGER,
     &       BUF(1), LBUF, POSITION, COMM, IERR_MPI )
      IF (NBLRB_PACKET.EQ.0) GOTO 600 
      IDEB = 1+NBLRB_ALREADY_SENT
      IFIN = IDEB + NBLRB_PACKET-1
      allocate(BLOCK(MAXI_CLUSTER), STAT=ALLOCOK )
      IF ( ALLOCOK .GT. 0 ) THEN
             IERR = -13
             IERROR = MAXI_CLUSTER
             GOTO 500
      END IF
      allocate(SCALED(MAXI_CLUSTER,2), STAT=ALLOCOK )
      IF ( ALLOCOK .GT. 0 ) THEN
             IERR = -13
             IERROR = 2*MAXI_CLUSTER
             GOTO 500
      END IF
      DO I = IDEB, IFIN
        IF (BLR(I)%ISLR) THEN
          ISLR_INT = 1
        ELSE
          ISLR_INT = 0
        ENDIF
        CALL MPI_PACK( ISLR_INT, 1, MPI_INTEGER,
     &       BUF(1), LBUF, POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( BLR(I)%K,
     &       1, MPI_INTEGER,
     &       BUF(1), LBUF, POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( BLR(I)%M, 
     &       1, MPI_INTEGER,
     &       BUF(1), LBUF, POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( BLR(I)%N,
     &       1, MPI_INTEGER,
     &       BUF(1), LBUF, POSITION, COMM, IERR_MPI )
        IF (BLR(I)%ISLR) THEN
          IF (BLR(I)%K .GT. 0) THEN
            CALL MPI_PACK( BLR(I)%Q(1,1), BLR(I)%M*BLR(I)%K,
     &       MPI_REAL,
     &       BUF(1), LBUF, POSITION, COMM, IERR_MPI )
            J =1
          DO WHILE (J <= BLR(I)%N)
              IF (IPIV(J) > 0) THEN
                SCALED(1:BLR(I)%K,1) = A(POSELTD+LD_DIAG*(J-1)+J-1) 
     &            * BLR(I)%R(1:BLR(I)%K,J)
                J = J+1
              CALL MPI_PACK( SCALED(1,1), BLR(I)%K,
     &           MPI_REAL,
     &           BUF(1), LBUF, POSITION, COMM, IERR_MPI )
              ELSE 
                PIV1    = A(POSELTD+LD_DIAG*(J-1)+J-1)
                PIV2    = A(POSELTD+LD_DIAG*J+J)
                OFFDIAG = A(POSELTD+LD_DIAG*(J-1)+J)
                BLOCK(1:BLR(I)%K)    = BLR(I)%R(1:BLR(I)%K,J)
                SCALED(1:BLR(I)%K,1) = PIV1 * BLR(I)%R(1:BLR(I)%K,J)
     &            + OFFDIAG * BLR(I)%R(1:BLR(I)%K,J+1)
                CALL MPI_PACK( SCALED(1,1), BLR(I)%K,
     &           MPI_REAL,
     &           BUF(1), LBUF, POSITION, COMM, IERR_MPI )
                SCALED(1:BLR(I)%K,2) = OFFDIAG * BLOCK(1:BLR(I)%K)
     &            + PIV2 * BLR(I)%R(1:BLR(I)%K,J+1)
                 J =J+2
                CALL MPI_PACK( SCALED(1,2), BLR(I)%K,
     &           MPI_REAL,
     &           BUF(1), LBUF, POSITION, COMM, IERR_MPI )
              ENDIF
          END DO
        ENDIF
        ELSE
          J = 1
          DO WHILE (J <= BLR(I)%N)
              IF (IPIV(J) > 0) THEN
                SCALED(1:BLR(I)%M,1) = A(POSELTD+LD_DIAG*(J-1)+J-1) 
     &           * BLR(I)%Q(1:BLR(I)%M,J)
                CALL MPI_PACK( SCALED(1,1), BLR(I)%M,
     &           MPI_REAL,
     &           BUF(1), LBUF, POSITION, COMM, IERR_MPI )
                J = J+1
              ELSE 
                PIV1    = A(POSELTD+LD_DIAG*(J-1)+J-1)
                PIV2    = A(POSELTD+LD_DIAG*J+J)
                OFFDIAG = A(POSELTD+LD_DIAG*(J-1)+J)
                BLOCK(1:BLR(I)%M)    = BLR(I)%Q(1:BLR(I)%M,J)
                SCALED(1:BLR(I)%M,1) = PIV1 * BLR(I)%Q(1:BLR(I)%M,J)
     &            + OFFDIAG * BLR(I)%Q(1:BLR(I)%M,J+1)
                CALL MPI_PACK( SCALED(1,1), BLR(I)%M,
     &           MPI_REAL,
     &           BUF(1), LBUF, POSITION, COMM, IERR_MPI )
                SCALED(1:BLR(I)%M,2) = OFFDIAG * BLOCK(1:BLR(I)%M)
     &            + PIV2 * BLR(I)%Q(1:BLR(I)%M,J+1)
                CALL MPI_PACK( SCALED(1,2), BLR(I)%M,
     &           MPI_REAL,
     &           BUF(1), LBUF, POSITION, COMM, IERR_MPI )
                 J=J+2
              ENDIF
          END DO
        ENDIF
      ENDDO
 500  CONTINUE
      IF (allocated(BLOCK)) deallocate(BLOCK)
      IF (allocated(SCALED)) deallocate(SCALED)
 600  CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_MPI_PACKSCALE_LR_PARTIAL
      END MODULE SMUMPS_BUF
