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
        MODULE MUMPS_BUF_COMMON
        PRIVATE
        PUBLIC :: BUF_ADJUST, BUF_LOOK, MUMPS_BUF_SIZE_AVAILABLE
        PUBLIC :: MUMPS_BUF_INIT, 
     &   MUMPS_BUF_INI_MYID,
     &   MUMPS_BUF_ALLOC_CB ,       MUMPS_BUF_DEALL_CB ,
     &   MUMPS_BUF_ALLOC_SMALL_BUF, MUMPS_BUF_DEALL_SMALL_BUF,
     &   MUMPS_BUF_ALLOC_LOAD_BUFFER,MUMPS_BUF_DEALL_LOAD_BUFFER,
     &   MUMPS_BUF_SEND_1INT,       MUMPS_BUF_SEND_DESC_BANDE,
     &   MUMPS_BUF_SEND_MAPLIG,
     &   MUMPS_BUF_SEND_RTNELIND,
     &   MUMPS_BUF_SEND_ROOT2SLAVE, MUMPS_BUF_SEND_ROOT2SON,
     &   MUMPS_BUF_SEND_UPDATE_LOAD, 
     &   MUMPS_BUF_DIST_IRECV_SIZE,
     &   MUMPS_BUF_BCAST_ARRAY, MUMPS_BUF_ALL_EMPTY,
     &   MUMPS_BUF_BROADCAST, MUMPS_BUF_SEND_NOT_MSTR,
     &   MUMPS_BUF_SEND_FILS,
     &   MUMPS_BUF_TEST
#if defined(USE_GPU)
        PUBLIC :: MUMPS_BUF_MAX_OUTER_T2
#endif
        INTEGER NEXT, REQ, CONTENT, OVHSIZE
        PARAMETER( NEXT = 0, REQ = 1, CONTENT = 2, OVHSIZE = 2 )
        PUBLIC :: OVHSIZE
        INTEGER, SAVE, PUBLIC :: SIZEofINT, SIZEofREAL, BUF_MYID
        TYPE MUMPS_COMM_BUFFER_TYPE
          INTEGER LBUF, HEAD, TAIL,LBUF_INT, ILASTMSG
          INTEGER, DIMENSION(:),POINTER :: CONTENT
        END TYPE MUMPS_COMM_BUFFER_TYPE
        TYPE ( MUMPS_COMM_BUFFER_TYPE ), SAVE, PUBLIC :: BUF_CB
        TYPE ( MUMPS_COMM_BUFFER_TYPE ), SAVE, PUBLIC :: BUF_SMALL
        TYPE ( MUMPS_COMM_BUFFER_TYPE ), SAVE, PUBLIC :: BUF_LOAD
        INTEGER, SAVE, PUBLIC :: SIZE_RBUF_BYTES
      CONTAINS
        SUBROUTINE MUMPS_BUF_INI_MYID( MYID )
        IMPLICIT NONE
        INTEGER MYID
        BUF_MYID  = MYID
        RETURN
        END SUBROUTINE MUMPS_BUF_INI_MYID
        SUBROUTINE MUMPS_BUF_INIT( IntSize, RealSize )
        IMPLICIT NONE
        INTEGER IntSize, RealSize
        SIZEofINT = IntSize
        SIZEofREAL = RealSize
        NULLIFY(BUF_CB  %CONTENT)
        NULLIFY(BUF_SMALL%CONTENT)
        NULLIFY(BUF_LOAD%CONTENT)
        BUF_CB%LBUF     = 0
        BUF_CB%LBUF_INT = 0
        BUF_CB%HEAD     = 1
        BUF_CB%TAIL     = 1
        BUF_CB%ILASTMSG = 1
        BUF_SMALL%LBUF     = 0
        BUF_SMALL%LBUF_INT = 0
        BUF_SMALL%HEAD     = 1
        BUF_SMALL%TAIL     = 1
        BUF_SMALL%ILASTMSG = 1
        BUF_LOAD%LBUF     = 0
        BUF_LOAD%LBUF_INT = 0
        BUF_LOAD%HEAD     = 1
        BUF_LOAD%TAIL     = 1
        BUF_LOAD%ILASTMSG = 1
        RETURN
        END SUBROUTINE MUMPS_BUF_INIT
        SUBROUTINE MUMPS_BUF_ALLOC_CB( SIZE, IERR )
        IMPLICIT NONE
        INTEGER SIZE, IERR
        CALL BUF_ALLOC( BUF_CB, SIZE, IERR )
        RETURN
        END SUBROUTINE MUMPS_BUF_ALLOC_CB
        SUBROUTINE MUMPS_BUF_ALLOC_SMALL_BUF( SIZE, IERR )
        IMPLICIT NONE
        INTEGER SIZE, IERR
        CALL BUF_ALLOC( BUF_SMALL, SIZE, IERR )
        RETURN
        END SUBROUTINE MUMPS_BUF_ALLOC_SMALL_BUF
        SUBROUTINE MUMPS_BUF_ALLOC_LOAD_BUFFER( SIZE, IERR )
        IMPLICIT NONE
        INTEGER SIZE, IERR
        CALL BUF_ALLOC( BUF_LOAD, SIZE, IERR )        
        RETURN
        END SUBROUTINE MUMPS_BUF_ALLOC_LOAD_BUFFER
        SUBROUTINE MUMPS_BUF_DEALL_LOAD_BUFFER( IERR )
        IMPLICIT NONE
        INTEGER IERR
        CALL BUF_DEALL( BUF_LOAD, IERR )
        RETURN
        END SUBROUTINE MUMPS_BUF_DEALL_LOAD_BUFFER
        SUBROUTINE MUMPS_BUF_DEALL_CB( IERR )
        IMPLICIT NONE
        INTEGER IERR
        CALL BUF_DEALL( BUF_CB, IERR )
        RETURN
        END SUBROUTINE MUMPS_BUF_DEALL_CB
        SUBROUTINE MUMPS_BUF_DEALL_SMALL_BUF( IERR )
        IMPLICIT NONE
        INTEGER IERR
        CALL BUF_DEALL( BUF_SMALL, IERR )
        RETURN
        END SUBROUTINE MUMPS_BUF_DEALL_SMALL_BUF
        SUBROUTINE BUF_ALLOC( BUF, SIZE, IERR )
        IMPLICIT NONE
        TYPE ( MUMPS_COMM_BUFFER_TYPE ) :: BUF
        INTEGER SIZE, IERR
        IERR         = 0
        BUF%LBUF     = SIZE
        BUF%LBUF_INT = ( SIZE + SIZEofINT - 1 ) / SIZEofINT
        IF ( associated ( BUF%CONTENT ) ) DEALLOCATE( BUF%CONTENT )
        ALLOCATE( BUF%CONTENT( BUF%LBUF_INT ), stat = IERR )
       IF (IERR .NE. 0) THEN
          NULLIFY( BUF%CONTENT )
          IERR         = -1
          BUF%LBUF     =  0
          BUF%LBUF_INT =  0
        END IF
        BUF%HEAD     = 1
        BUF%TAIL     = 1
        BUF%ILASTMSG = 1
        RETURN
        END SUBROUTINE BUF_ALLOC
        SUBROUTINE BUF_DEALL( BUF, IERR )
        IMPLICIT NONE
        TYPE ( MUMPS_COMM_BUFFER_TYPE ) :: BUF
        INTEGER :: IERR
        INCLUDE 'mpif.h'
        INTEGER :: IERR_MPI
        INTEGER :: STATUS(MPI_STATUS_SIZE)
        LOGICAL :: FLAG
        IF ( .NOT. associated ( BUF%CONTENT ) ) THEN
          BUF%HEAD     = 1
          BUF%LBUF     = 0
          BUF%LBUF_INT = 0
          BUF%TAIL     = 1
          BUF%ILASTMSG = 1
          RETURN
        END IF
        DO WHILE ( BUF%HEAD.NE.0 .AND. BUF%HEAD .NE. BUF%TAIL )
          CALL MPI_TEST(BUF%CONTENT( BUF%HEAD + REQ ), FLAG,
     &                  STATUS, IERR_MPI)
          IF ( .not. FLAG ) THEN
            WRITE(*,*) '** Warning: trying to cancel a request.'
            WRITE(*,*) '** This might be problematic'
            CALL MPI_CANCEL( BUF%CONTENT( BUF%HEAD + REQ ), IERR_MPI )
            CALL MPI_REQUEST_FREE( BUF%CONTENT( BUF%HEAD + REQ ),
     &                             IERR_MPI )
          END IF
          BUF%HEAD = BUF%CONTENT( BUF%HEAD + NEXT )
        END DO
        DEALLOCATE( BUF%CONTENT )
        NULLIFY( BUF%CONTENT )
        BUF%LBUF     = 0
        BUF%LBUF_INT = 0
        BUF%HEAD     = 1
        BUF%TAIL     = 1
        BUF%ILASTMSG = 1
        RETURN
        END SUBROUTINE BUF_DEALL
        SUBROUTINE MUMPS_BUF_SEND_1INT( I, DEST, TAG, COMM,
     &                                   KEEP, IERR )
        IMPLICIT NONE
        INTEGER I
        INTEGER DEST, TAG, COMM, IERR
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INCLUDE 'mpif.h'
        INTEGER :: IERR_MPI
        INTEGER IPOS, IREQ, MSG_SIZE, POSITION
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1)=DEST
        IERR = 0
        CALL MPI_PACK_SIZE( 1, MPI_INTEGER,
     &                      COMM, MSG_SIZE, IERR_MPI )
        CALL BUF_LOOK( BUF_SMALL, IPOS, IREQ, MSG_SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
         write(6,*) ' Internal error in MUMPS_BUF_SEND_1INT',
     &       ' Buf size (bytes)= ',BUF_SMALL%LBUF
      RETURN
        ENDIF
        POSITION=0
        CALL MPI_PACK( I, 1,
     &                 MPI_INTEGER, BUF_SMALL%CONTENT( IPOS ),
     &                 MSG_SIZE,
     &                 POSITION, COMM, IERR_MPI )
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_SMALL%CONTENT(IPOS), MSG_SIZE,
     &                  MPI_PACKED, DEST, TAG, COMM,
     &                  BUF_SMALL%CONTENT( IREQ ), IERR_MPI )
        RETURN
        END SUBROUTINE MUMPS_BUF_SEND_1INT
        SUBROUTINE MUMPS_BUF_ALL_EMPTY(CHECK_COMM_NODES,
     &             CHECK_COMM_LOAD,FLAG)
        LOGICAL, INTENT(IN)  :: CHECK_COMM_NODES, CHECK_COMM_LOAD
        LOGICAL, INTENT(OUT) :: FLAG
        LOGICAL FLAG1, FLAG2, FLAG3
        FLAG = .TRUE.
        IF (CHECK_COMM_NODES) THEN
          CALL MUMPS_BUF_EMPTY( BUF_SMALL, FLAG1 )
          CALL MUMPS_BUF_EMPTY( BUF_CB, FLAG2 )
          FLAG = FLAG .AND. FLAG1 .AND. FLAG2
        ENDIF
        IF ( CHECK_COMM_LOAD ) THEN
          CALL MUMPS_BUF_EMPTY( BUF_LOAD, FLAG3 )
          FLAG = FLAG .AND. FLAG3
        ENDIF
        RETURN
        END SUBROUTINE MUMPS_BUF_ALL_EMPTY
        SUBROUTINE MUMPS_BUF_EMPTY( B, FLAG )
        TYPE ( MUMPS_COMM_BUFFER_TYPE ) :: B
        LOGICAL :: FLAG
        INTEGER SIZE_AVAIL
        CALL MUMPS_BUF_SIZE_AVAILABLE(B, SIZE_AVAIL)
        FLAG = ( B%HEAD == B%TAIL )
        RETURN
        END SUBROUTINE MUMPS_BUF_EMPTY
        SUBROUTINE MUMPS_BUF_FREEREQUESTS( B )
        IMPLICIT NONE
        TYPE ( MUMPS_COMM_BUFFER_TYPE ) :: B
        INCLUDE 'mpif.h'
        INTEGER :: IERR_MPI, CURRENT, LAST_NOT_FREE, LAST_NOT_FREE_TAIL
        INTEGER :: STATUS(MPI_STATUS_SIZE)
        LOGICAL :: FLAG, BROADCAST_NOT_FREE
        IF ( B%HEAD .NE. B%TAIL ) THEN
 10       CONTINUE
          CALL MPI_TEST( B%CONTENT( B%HEAD + REQ ), FLAG, STATUS,
     &                   IERR_MPI )
          IF ( FLAG ) THEN
            B%HEAD = B%CONTENT( B%HEAD + NEXT )
            IF ( B%HEAD .EQ. 0 ) B%HEAD = B%TAIL
            IF ( B%HEAD .NE. B%TAIL ) GOTO 10
          ELSE
            LAST_NOT_FREE = B%HEAD
            CURRENT = B%CONTENT( LAST_NOT_FREE + NEXT )
            LAST_NOT_FREE_TAIL = CURRENT
            BROADCAST_NOT_FREE = B%CONTENT(LAST_NOT_FREE+NEXT).EQ.
     &                              LAST_NOT_FREE+OVHSIZE 
            DO WHILE ( CURRENT .NE. 0 )
              IF (BROADCAST_NOT_FREE) THEN
                FLAG = .FALSE.
              ELSE
                CALL MPI_TEST( B%CONTENT( CURRENT + REQ ), FLAG, STATUS,
     &                       IERR_MPI )
              ENDIF
              IF (FLAG) THEN
                CURRENT = B%CONTENT( CURRENT + NEXT )
                B%CONTENT( LAST_NOT_FREE + NEXT ) = CURRENT
              ELSE
                LAST_NOT_FREE = CURRENT
                CURRENT = B%CONTENT( CURRENT + NEXT )
                IF ( CURRENT .NE. 0 ) THEN
                  LAST_NOT_FREE_TAIL = CURRENT
                ELSE
                  LAST_NOT_FREE_TAIL = B%TAIL
                ENDIF
                BROADCAST_NOT_FREE = B%CONTENT(LAST_NOT_FREE+NEXT).EQ.
     &                               LAST_NOT_FREE+OVHSIZE 
              ENDIF
            ENDDO
            IF ( LAST_NOT_FREE_TAIL .NE. 0 ) THEN
              B%TAIL = LAST_NOT_FREE_TAIL
              B%ILASTMSG = LAST_NOT_FREE
            ELSE IF (B%ILASTMSG .NE. LAST_NOT_FREE) THEN
              WRITE(*,*) "ABORT", B%ILASTMSG, LAST_NOT_FREE
              CALL MUMPS_ABORT()
            ENDIF
          END IF
        END IF
        IF ( B%HEAD .EQ. B%TAIL ) THEN
          B%HEAD = 1
          B%TAIL = 1
          B%ILASTMSG = 1
        END IF
        RETURN
        END SUBROUTINE MUMPS_BUF_FREEREQUESTS
        SUBROUTINE MUMPS_BUF_SIZE_AVAILABLE( B, SIZE_AV, SIZE_AV2 )
        IMPLICIT NONE
        TYPE ( MUMPS_COMM_BUFFER_TYPE ) :: B
        INTEGER, INTENT(OUT) :: SIZE_AV
        INTEGER, OPTIONAL, INTENT(OUT) :: SIZE_AV2
        CALL MUMPS_BUF_FREEREQUESTS( B )
        IF ( B%HEAD .LE. B%TAIL ) THEN
          SIZE_AV = max( B%LBUF_INT - B%TAIL, B%HEAD - 2 )
          IF ( B%LBUF_INT - B%TAIL .GT. B%HEAD - 2 ) THEN
            SIZE_AV  = B%LBUF_INT - B%TAIL
            IF (present(SIZE_AV2)) SIZE_AV2 = 0
          ELSE
            SIZE_AV = B%HEAD - 2
            IF (present(SIZE_AV2)) SIZE_AV2 = B%LBUF_INT - B%TAIL
          ENDIF
        ELSE
          SIZE_AV  = B%HEAD - B%TAIL - 1
          IF (present(SIZE_AV2)) SIZE_AV2 = 0
        END IF
        SIZE_AV = max(0,SIZE_AV - OVHSIZE)
        SIZE_AV = SIZE_AV * SIZEofINT
        IF (present(SIZE_AV2)) THEN
          IF (SIZE_AV2 .NE. 0) THEN
            SIZE_AV = max(0,SIZE_AV2 - OVHSIZE)
            SIZE_AV2 = SIZE_AV2 * SIZEofINT
          ENDIF
        ENDIF
        RETURN
        END SUBROUTINE MUMPS_BUF_SIZE_AVAILABLE
        SUBROUTINE MUMPS_BUF_TEST()
        CALL MUMPS_BUF_FREEREQUESTS(BUF_CB)
        RETURN
        END SUBROUTINE MUMPS_BUF_TEST
        SUBROUTINE BUF_LOOK( B, IPOS, IREQ, MSG_SIZE, IERR, 
     &    NDEST , PDEST )
        IMPLICIT NONE
        TYPE ( MUMPS_COMM_BUFFER_TYPE ) :: B
        INTEGER, INTENT(IN)        :: MSG_SIZE
        INTEGER, INTENT(OUT)       :: IPOS, IREQ, IERR
        INTEGER NDEST
        INTEGER, INTENT(IN)        :: PDEST(max(1,NDEST))
        INCLUDE 'mpif.h'
        INTEGER :: MSG_SIZE_INT
        INTEGER :: IBUF
        IERR = 0
        CALL MUMPS_BUF_FREEREQUESTS( B )
        MSG_SIZE_INT = ( MSG_SIZE + ( SIZEofINT - 1 ) ) / SIZEofINT
        MSG_SIZE_INT = MSG_SIZE_INT + OVHSIZE
        IF ( MSG_SIZE_INT .GT. B%LBUF_INT - 1 ) THEN
          IERR = -2
          IPOS = -1
          IREQ = -1
          RETURN
        END IF
        IF ( B%HEAD .LE. B%TAIL ) THEN
          IF ( MSG_SIZE_INT .LE. B%LBUF_INT - B%TAIL ) THEN
            IBUF = B%TAIL
          ELSE IF ( MSG_SIZE_INT .LE. B%HEAD - 2 ) THEN
            IBUF = 1
          ELSE
            IERR = -1
          END IF
        ELSE
          IF ( MSG_SIZE_INT .LE . B%HEAD - B%TAIL - 1) THEN
            IBUF = B%TAIL
          ELSE
            IERR = -1
          ENDIF
        END IF
        IF (IERR .LT. 0) RETURN
        B%CONTENT( B%ILASTMSG + NEXT ) = IBUF
        B%ILASTMSG = IBUF
        B%TAIL = IBUF + MSG_SIZE_INT
        B%CONTENT( IBUF + NEXT ) = 0
        IPOS = IBUF + CONTENT
        IREQ = IBUF + REQ
        RETURN
        END SUBROUTINE BUF_LOOK
        SUBROUTINE BUF_ADJUST( BUF, SIZE )
        IMPLICIT NONE
        TYPE ( MUMPS_COMM_BUFFER_TYPE ) :: BUF
        INTEGER SIZE
        INTEGER SIZE_INT
        SIZE_INT = ( SIZE + SIZEofINT - 1 ) / SIZEofINT
        SIZE_INT = SIZE_INT + OVHSIZE
        BUF%TAIL = BUF%ILASTMSG + SIZE_INT
        RETURN
        END SUBROUTINE BUF_ADJUST
      SUBROUTINE MUMPS_BUF_SEND_DESC_BANDE(
     &             INODE, NBPROCFILS, NLIG, ILIG, NCOL, ICOL,
     &             NASS, NSLAVES_HDR, LIST_SLAVES,
     &             NSLAVES, XXG_MASTER,
     &             ESTIM_NFS4FATHER_ATSON,
     &             DEST, NFRONT, COMM, KEEP, IERR
     &             , LRSTATUS
     &)
      IMPLICIT NONE
        INTEGER COMM, IERR, NFRONT
        INTEGER, intent(in) :: INODE,NLIG, NCOL, NASS
        INTEGER, intent(in) :: NSLAVES_HDR, NSLAVES, XXG_MASTER
        INTEGER, intent(in) :: ESTIM_NFS4FATHER_ATSON
        INTEGER NBPROCFILS, DEST
        INTEGER ILIG( NLIG )
        INTEGER ICOL( NCOL )
        INTEGER LIST_SLAVES( max(NSLAVES_HDR,1) )
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INTEGER, INTENT(IN) :: LRSTATUS
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER SIZE_INT, SIZE_BYTES, POSITION, IPOS, IREQ
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1) = DEST
        IERR = 0
        SIZE_INT = ( 12 + NLIG + NCOL + NSLAVES_HDR )
        SIZE_BYTES = SIZE_INT * SIZEofINT
        IF (SIZE_INT.GT.SIZE_RBUF_BYTES ) THEN
         IERR = -3
      RETURN
        END IF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE_BYTES, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
      RETURN
        ENDIF
        POSITION = IPOS
        BUF_CB%CONTENT( POSITION ) = SIZE_INT
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = INODE
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NBPROCFILS
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NLIG
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NCOL
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NASS
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NFRONT
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NSLAVES_HDR
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NSLAVES
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = LRSTATUS
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = ESTIM_NFS4FATHER_ATSON
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = XXG_MASTER
        POSITION = POSITION + 1
        IF (NSLAVES_HDR.GT.0) THEN
         BUF_CB%CONTENT( POSITION: POSITION + NSLAVES_HDR - 1 ) = 
     &   LIST_SLAVES( 1: NSLAVES_HDR )
         POSITION = POSITION + NSLAVES_HDR
        ENDIF
        BUF_CB%CONTENT( POSITION:POSITION + NLIG - 1 ) = ILIG
        POSITION = POSITION + NLIG
        BUF_CB%CONTENT( POSITION:POSITION + NCOL - 1 ) = ICOL
        POSITION = POSITION + NCOL
        POSITION = POSITION - IPOS
        IF ( POSITION * SIZEofINT .NE. SIZE_BYTES ) THEN
          WRITE(*,*) 'Error in MUMPS_BUF_SEND_DESC_BANDE :',
     &               ' wrong estimated size'
          CALL MUMPS_ABORT()
        END IF
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), SIZE_BYTES,
     &                  MPI_PACKED,
     &                  DEST, MAITRE_DESC_BANDE, COMM,
     &                  BUF_CB%CONTENT( IREQ ), IERR_MPI )
        RETURN
        END SUBROUTINE MUMPS_BUF_SEND_DESC_BANDE
        SUBROUTINE MUMPS_BUF_SEND_MAPLIG( 
     &                INODE, NFRONT, NASS1, NFS4FATHER,
     &                ISON, MYID, NSLAVES, SLAVES_PERE,
     &                TROW, NCBSON,
     &                COMM, IERR,
     &                DEST, NDEST, SLAVEF, 
     & 
     &                KEEP,KEEP8, STEP, N, 
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &
     &                                  )
        IMPLICIT NONE
      INTEGER INODE, NFRONT, NASS1, NCBSON, NSLAVES, 
     &          NDEST
      INTEGER SLAVEF, MYID, ISON
      INTEGER TROW( NCBSON )
      INTEGER DEST( NDEST )
      INTEGER SLAVES_PERE( NSLAVES )
      INTEGER COMM, IERR
      INTEGER KEEP(500), N
      INTEGER(8) KEEP8(150)
      INTEGER STEP(N), 
     &        ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER SIZE_AV, IDEST, NSEND, SIZE, NFS4FATHER
        INTEGER TROW_SIZE, POSITION, INDX, INIV2
        INTEGER IPOS, IREQ
        INTEGER IONE
        PARAMETER ( IONE=1 )
        IERR = 0
        IF ( NDEST .eq. 1 ) THEN
          IF ( DEST(1).EQ.MYID )  GOTO 500
          SIZE = SIZEofINT * ( 7 + NSLAVES + NCBSON )
          IF ( NSLAVES.GT.0 ) THEN
             SIZE = SIZE + SIZEofINT * ( NSLAVES + 1 )
          ENDIF
          IF (SIZE.GT.SIZE_RBUF_BYTES ) THEN
            IERR = -3
            RETURN
          END IF
          CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE, DEST
     &                 )
          IF (IERR .LT. 0 ) THEN
            RETURN
          ENDIF
              POSITION = IPOS
              BUF_CB%CONTENT( POSITION ) = INODE
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = ISON
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NSLAVES
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NFRONT
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NASS1
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NCBSON
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NFS4FATHER
              POSITION = POSITION + 1
              IF ( NSLAVES.GT.0 ) THEN
                INIV2 = ISTEP_TO_INIV2 ( STEP(INODE) )
                BUF_CB%CONTENT( POSITION: POSITION + NSLAVES )
     &          =  TAB_POS_IN_PERE(1:NSLAVES+1,INIV2)
                POSITION = POSITION + NSLAVES + 1
              ENDIF
              IF ( NSLAVES .NE. 0 ) THEN
                BUF_CB%CONTENT( POSITION: POSITION + NSLAVES - 1 )
     &          = SLAVES_PERE( 1: NSLAVES )
                POSITION = POSITION + NSLAVES
              END IF
              BUF_CB%CONTENT( POSITION:POSITION+NCBSON-1 ) =
     &        TROW( 1: NCBSON )
              POSITION = POSITION + NCBSON
              POSITION = POSITION - IPOS
              IF ( POSITION * SIZEofINT .NE. SIZE ) THEN
                WRITE(*,*) 'Error in MUMPS_BUF_SEND_MAPLIG :',
     &                     ' wrong estimated size'
                CALL MUMPS_ABORT()
              END IF
              KEEP(266)=KEEP(266)+1
              CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), SIZE,
     &                        MPI_PACKED,
     &                        DEST( NDEST ), MAPLIG, COMM,
     &                        BUF_CB%CONTENT( IREQ ),
     &                        IERR_MPI )
        ELSE
          NSEND = 0
          DO IDEST = 1, NDEST
            IF ( DEST( IDEST ) .ne. MYID ) NSEND = NSEND + 1
          END DO
          SIZE = SIZEofINT * 
     &         ( ( OVHSIZE + 7 + NSLAVES )* NSEND + NCBSON )
          IF ( NSLAVES.GT.0 ) THEN
           SIZE = SIZE + SIZEofINT * NSEND*( NSLAVES + 1 )
          ENDIF
          CALL MUMPS_BUF_SIZE_AVAILABLE( BUF_CB, SIZE_AV )
          IF ( SIZE_AV .LT. SIZE ) THEN
            IERR = -1
            RETURN
          END IF
          DO IDEST= 1, NDEST
            CALL MUMPS_BLOC2_GET_SLAVE_INFO( 
     &                KEEP,KEEP8, ISON, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                IDEST, NCBSON, 
     &                NDEST, 
     &                TROW_SIZE, INDX  )
            SIZE = SIZEofINT * ( NSLAVES + TROW_SIZE + 7 )
            IF ( NSLAVES.GT.0 ) THEN
             SIZE = SIZE + SIZEofINT * ( NSLAVES + 1 )
            ENDIF
            IF ( MYID .NE. DEST( IDEST ) ) THEN
               IF (SIZE.GT.SIZE_RBUF_BYTES) THEN
                IERR = -3
                RETURN
              ENDIF
              CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE, IERR,
     &                       IONE, DEST(IDEST) )
              IF ( IERR .LT. 0 )  THEN
                WRITE(*,*) 'Internal error MUMPS_BUF_SEND_MAPLIG',
     &                     'IERR after BUF_LOOK=',IERR
                CALL MUMPS_ABORT()
              END IF
              POSITION = IPOS
              BUF_CB%CONTENT( POSITION ) = INODE
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = ISON
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NSLAVES
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NFRONT
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NASS1
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = TROW_SIZE
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NFS4FATHER
              POSITION = POSITION + 1
              IF ( NSLAVES.GT.0 ) THEN
                INIV2 = ISTEP_TO_INIV2 ( STEP(INODE) )
                BUF_CB%CONTENT( POSITION: POSITION + NSLAVES )
     &          =  TAB_POS_IN_PERE(1:NSLAVES+1,INIV2)
                POSITION = POSITION + NSLAVES + 1
              ENDIF
              IF ( NSLAVES .NE. 0 ) THEN
                BUF_CB%CONTENT( POSITION: POSITION + NSLAVES - 1 )
     &          = SLAVES_PERE( 1: NSLAVES )
                POSITION = POSITION + NSLAVES
              END IF
              BUF_CB%CONTENT( POSITION:POSITION+TROW_SIZE-1 ) =
     &        TROW( INDX: INDX + TROW_SIZE - 1 )
              POSITION = POSITION + TROW_SIZE
              POSITION = POSITION - IPOS
              IF ( POSITION * SIZEofINT .NE. SIZE ) THEN
               WRITE(*,*) ' ERROR 1 in TRY_SEND_MAPLIG:',
     &          'Wrong estimated size'
               CALL MUMPS_ABORT()
              END IF
              KEEP(266)=KEEP(266)+1
              CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), SIZE,
     &                        MPI_PACKED,
     &                        DEST( IDEST ), MAPLIG, COMM,
     &                        BUF_CB%CONTENT( IREQ ),
     &                        IERR_MPI )
            END IF
          END DO
        END IF
 500    CONTINUE
        RETURN
        END SUBROUTINE MUMPS_BUF_SEND_MAPLIG
        SUBROUTINE MUMPS_BUF_SEND_RTNELIND( ISON, NELIM,
     &             NELIM_ROW, NELIM_COL, NSLAVES, SLAVES,
     &             DEST, COMM, KEEP, IERR )
        INTEGER ISON, NELIM
        INTEGER NSLAVES, DEST, COMM, IERR
        INTEGER NELIM_ROW( NELIM ), NELIM_COL( NELIM )
        INTEGER SLAVES( NSLAVES )
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER SIZE, POSITION, IPOS, IREQ
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1) = DEST
        IERR = 0
        SIZE = ( 3 + NSLAVES + 2 * NELIM ) * SIZEofINT
        IF (SIZE.GT.SIZE_RBUF_BYTES) THEN
             IERR = -3
      RETURN
        ENDIF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE, DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        POSITION = IPOS
        BUF_CB%CONTENT( POSITION ) = ISON
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NELIM
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NSLAVES
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION: POSITION + NELIM - 1 ) = NELIM_ROW
        POSITION = POSITION + NELIM
        BUF_CB%CONTENT( POSITION: POSITION + NELIM - 1 ) = NELIM_COL
        POSITION = POSITION + NELIM
        BUF_CB%CONTENT( POSITION: POSITION + NSLAVES - 1 ) = SLAVES
        POSITION = POSITION + NSLAVES
        POSITION = POSITION - IPOS
        IF ( POSITION * SIZEofINT .NE. SIZE ) THEN
          WRITE(*,*) 'Error in MUMPS_BUF_SEND_ROOT_NELIM_INDICES:',
     &               'wrong estimated size'
          CALL MUMPS_ABORT()
        END IF
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), SIZE, 
     &                  MPI_PACKED,
     &                  DEST, ROOT_NELIM_INDICES, COMM,
     &                  BUF_CB%CONTENT( IREQ ), IERR_MPI )
        RETURN
        END SUBROUTINE MUMPS_BUF_SEND_RTNELIND
        SUBROUTINE MUMPS_BUF_SEND_ROOT2SON( ISON, NELIM_ROOT,
     &             DEST, COMM, KEEP, IERR )
        IMPLICIT NONE
        INTEGER ISON, NELIM_ROOT, DEST, COMM, IERR
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER IPOS, IREQ, SIZE
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1)=DEST
        IERR = 0
        SIZE = 2 * SIZEofINT
        CALL BUF_LOOK( BUF_SMALL, IPOS, IREQ, SIZE, IERR,
     &                 IONE, DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
          WRITE(*,*) 'Internal error 1 with small buffers '
          CALL MUMPS_ABORT()
        END IF
        IF ( IERR .LT. 0 ) THEN
          RETURN
        ENDIF
        BUF_SMALL%CONTENT( IPOS )     = ISON
        BUF_SMALL%CONTENT( IPOS + 1 ) = NELIM_ROOT
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_SMALL%CONTENT( IPOS ), SIZE, 
     &                  MPI_PACKED,
     &                  DEST, ROOT_2SON, COMM,
     &                  BUF_SMALL%CONTENT( IREQ ), IERR_MPI )
        RETURN
        END SUBROUTINE MUMPS_BUF_SEND_ROOT2SON
        SUBROUTINE MUMPS_BUF_SEND_ROOT2SLAVE
     &  ( TOT_ROOT_SIZE, TOT_CONT2RECV, DEST, COMM, KEEP, IERR )
        IMPLICIT NONE
        INTEGER TOT_ROOT_SIZE, TOT_CONT2RECV, DEST, COMM, IERR
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER SIZE, IPOS, IREQ
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        IERR = 0
        DEST2(1) = DEST
        SIZE = 2 * SIZEofINT
        CALL BUF_LOOK( BUF_SMALL, IPOS, IREQ, SIZE, IERR,
     &                 IONE, DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
          WRITE(*,*) 'Internal error 2 with small buffers '
          CALL MUMPS_ABORT()
        END IF
        IF ( IERR .LT. 0 ) THEN
      RETURN
        ENDIF
        BUF_SMALL%CONTENT( IPOS     ) = TOT_ROOT_SIZE
        BUF_SMALL%CONTENT( IPOS + 1 ) = TOT_CONT2RECV
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_SMALL%CONTENT( IPOS ), SIZE, 
     &                  MPI_PACKED,
     &                  DEST, ROOT_2SLAVE, COMM,
     &                  BUF_SMALL%CONTENT( IREQ ), IERR_MPI )
        RETURN
        END SUBROUTINE MUMPS_BUF_SEND_ROOT2SLAVE
        SUBROUTINE MUMPS_BUF_SEND_UPDATE_LOAD
     &             ( BDC_SBTR,BDC_MEM,BDC_MD, COMM, NPROCS, LOAD,
     &               MEM,SBTR_CUR,
     &               LU_USAGE,
     &               FUTURE_NIV2,
     &               MYID, KEEP, IERR)
        IMPLICIT NONE
        INTEGER COMM, NPROCS, MYID, IERR
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INTEGER FUTURE_NIV2(NPROCS)
        DOUBLE PRECISION LU_USAGE
        DOUBLE PRECISION LOAD
        DOUBLE PRECISION MEM,SBTR_CUR
        LOGICAL BDC_MEM,BDC_SBTR,BDC_MD
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER POSITION, IREQ, IPOS, SIZE1, SIZE2, SIZE
        INTEGER I, NDEST, IDEST, IPOSMSG, WHAT, NREALS
        INTEGER IZERO
        INTEGER MYID2(1)
        PARAMETER ( IZERO=0 )
        IERR = 0
        MYID2(1) = MYID
        NDEST = NPROCS - 1
        NDEST = 0
        DO I = 1, NPROCS
           IF ( I .NE. MYID + 1 .AND. FUTURE_NIV2(I).NE.0) THEN
              NDEST = NDEST + 1
           ENDIF
        ENDDO
        IF ( NDEST .eq. 0 ) THEN
           RETURN
        ENDIF
        CALL MPI_PACK_SIZE( 1 + (NDEST-1) * OVHSIZE, 
     &                       MPI_INTEGER, COMM,
     &                       SIZE1, IERR_MPI )
        NREALS = 1
        IF (BDC_MEM) THEN
          NREALS = 2
        ENDIf
        IF (BDC_SBTR)THEN
          NREALS = 3
        ENDIF
        IF(BDC_MD)THEN
           NREALS=NREALS+1
        ENDIF
        CALL MPI_PACK_SIZE( NREALS, MPI_DOUBLE_PRECISION,
     &                      COMM, SIZE2, IERR_MPI )
        SIZE = SIZE1 + SIZE2
        CALL BUF_LOOK( BUF_LOAD, IPOS, IREQ, SIZE, IERR, 
     &                  IZERO, MYID2 
     &               )
        IF ( IERR .LT. 0 ) THEN
      RETURN
        ENDIF
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        WHAT = 0  
        POSITION = 0
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     &                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                 POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( LOAD, 1, MPI_DOUBLE_PRECISION,
     &                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                 POSITION, COMM, IERR_MPI )
        IF (BDC_MEM) THEN
          CALL MPI_PACK( MEM, 1, MPI_DOUBLE_PRECISION,
     &                   BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                   POSITION, COMM, IERR_MPI )
        END IF
        IF (BDC_SBTR) THEN
          CALL MPI_PACK( SBTR_CUR, 1, MPI_DOUBLE_PRECISION,
     &                   BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                   POSITION, COMM, IERR_MPI )
        END IF
        IF(BDC_MD)THEN
           CALL MPI_PACK( LU_USAGE, 1, MPI_DOUBLE_PRECISION,
     &          BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &          POSITION, COMM, IERR_MPI )
        ENDIF
        IDEST = 0
        DO I = 0, NPROCS - 1
        IF ( I .NE. MYID .AND. FUTURE_NIV2(I+1) .NE. 0) THEN
            IDEST = IDEST + 1
            KEEP(267)=KEEP(267)+1
            CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     &                      POSITION, MPI_PACKED, I,
     &                      UPDATE_LOAD, COMM,
     &                      BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     &                      IERR_MPI )
          END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in MUMPS_BUF_SEND_UPDATE_LOAD'
          WRITE(*,*) ' Size,position=',SIZE,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE MUMPS_BUF_SEND_UPDATE_LOAD
        SUBROUTINE MUMPS_BUF_BROADCAST
     &             ( WHAT, COMM, NPROCS, 
     &               FUTURE_NIV2,
     &               LOAD, UPD_LOAD,
     &               MYID, KEEP267, IERR)
        IMPLICIT NONE
        INTEGER COMM, NPROCS, MYID, IERR, WHAT
        DOUBLE PRECISION LOAD,UPD_LOAD
        INTEGER, INTENT(INOUT) :: KEEP267
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER POSITION, IREQ, IPOS, SIZE1, SIZE2, SIZE
        INTEGER I, NDEST, IDEST, IPOSMSG, NREALS
        INTEGER IZERO
        INTEGER MYID2(1)
        INTEGER FUTURE_NIV2(NPROCS)
        PARAMETER ( IZERO=0 )
        IERR = 0
        IF (WHAT .NE. 2 .AND. WHAT .NE. 3 .AND.
     &       WHAT.NE.6.AND. WHAT.NE.8 .AND.WHAT.NE.9.AND.
     &       WHAT.NE.17) THEN
          WRITE(*,*)
     &          "Internal error 1 in MUMPS_BUF_BROADCAST",WHAT
        END IF
        MYID2(1) = MYID
        NDEST = NPROCS - 1
        NDEST = 0
        DO I = 1, NPROCS
          IF ( I .NE. MYID + 1 .AND. FUTURE_NIV2(I).NE.0) THEN
            NDEST = NDEST + 1
          ENDIF
        ENDDO
        IF ( NDEST .eq. 0 ) THEN
           RETURN
        ENDIF
        CALL MPI_PACK_SIZE( 1 + (NDEST-1) * OVHSIZE, 
     &                       MPI_INTEGER, COMM,
     &                       SIZE1, IERR_MPI )
        IF((WHAT.NE.17).AND.(WHAT.NE.10))THEN
           NREALS = 1
        ELSE
           NREALS = 2
        ENDIF
        CALL MPI_PACK_SIZE( NREALS, MPI_DOUBLE_PRECISION,
     &                      COMM, SIZE2, IERR_MPI )
        SIZE = SIZE1 + SIZE2
        CALL BUF_LOOK( BUF_LOAD, IPOS, IREQ, SIZE, IERR, 
     &                  IZERO, MYID2 
     &               )
        IF ( IERR .LT. 0 ) THEN
      RETURN
        ENDIF
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     &                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                 POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( LOAD, 1, MPI_DOUBLE_PRECISION,
     &                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                 POSITION, COMM, IERR_MPI )
        IF((WHAT.EQ.17).OR.(WHAT.EQ.10))THEN
           CALL MPI_PACK( UPD_LOAD, 1, MPI_DOUBLE_PRECISION,
     &          BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &          POSITION, COMM, IERR_MPI )
        ENDIF
        IDEST = 0
        DO I = 0, NPROCS - 1
          IF ( I .NE. MYID .AND. FUTURE_NIV2(I+1) .NE. 0) THEN
            IDEST = IDEST + 1
            KEEP267 = KEEP267 + 1
            CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     &                      POSITION, MPI_PACKED, I,
     &                      UPDATE_LOAD, COMM,
     &                      BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     &                      IERR_MPI )
          END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in MUMPS_BUF_BROADCAST'
          WRITE(*,*) ' Size,position=',SIZE,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE MUMPS_BUF_BROADCAST
        SUBROUTINE MUMPS_BUF_SEND_FILS
     &             ( WHAT, COMM, NPROCS,
     &               FATHER_NODE,INODE,NCB,KEEP,
     &               MYID,REMOTE, IERR)
        IMPLICIT NONE
        INTEGER COMM, NPROCS, MYID, IERR, WHAT,REMOTE
        INTEGER FATHER_NODE,INODE
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER POSITION, IREQ, IPOS, SIZE
        INTEGER NDEST, IDEST, IPOSMSG
        INTEGER IZERO,NCB,KEEP(500)
        INTEGER MYID2(1)
        PARAMETER ( IZERO=0 )
        MYID2(1) = MYID
        NDEST = 1
        IF ( NDEST .eq. 0 ) THEN
           RETURN
        ENDIF
        IF((KEEP(81).EQ.2).OR.(KEEP(81).EQ.3))THEN
           CALL MPI_PACK_SIZE( 4 + OVHSIZE, 
     &          MPI_INTEGER, COMM,
     &          SIZE, IERR_MPI )
        ELSE
           CALL MPI_PACK_SIZE( 2, 
     &          MPI_INTEGER, COMM,
     &          SIZE, IERR_MPI )
        ENDIF
        CALL BUF_LOOK( BUF_LOAD, IPOS, IREQ, SIZE, IERR, 
     &                  IZERO, MYID2 
     &               )
        IF ( IERR .LT. 0 ) THEN
      RETURN
        ENDIF
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     &                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                 POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( FATHER_NODE, 1, MPI_INTEGER,
     &                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                 POSITION, COMM, IERR_MPI )
        IF((KEEP(81).EQ.2).OR.(KEEP(81).EQ.3))THEN
           CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &          BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &          POSITION, COMM, IERR_MPI )
           CALL MPI_PACK( NCB, 1, MPI_INTEGER,
     &          BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &          POSITION, COMM, IERR_MPI )
        ENDIF
        IDEST = 1
        KEEP(267)=KEEP(267)+1
        CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     &                 POSITION, MPI_PACKED, REMOTE,
     &                 UPDATE_LOAD, COMM,
     &                 BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     &                 IERR_MPI )
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in MUMPS_BUF_SEND_FILS'
          WRITE(*,*) ' Size,position=',SIZE,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE MUMPS_BUF_SEND_FILS
        SUBROUTINE MUMPS_BUF_SEND_NOT_MSTR( COMM, MYID, NPROCS,
     &  MAX_SURF_MASTER, KEEP, IERR)
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER COMM, MYID, IERR, NPROCS
        DOUBLE PRECISION MAX_SURF_MASTER
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INTEGER :: IERR_MPI
        INTEGER IPOS, IREQ, IDEST, IPOSMSG, POSITION, I
        INTEGER IZERO
        INTEGER MYID2(1)
        PARAMETER ( IZERO=0 )
        INTEGER NDEST, NINTS, NREALS, SIZE, SIZE1, SIZE2
        INTEGER WHAT
        IERR = 0
        MYID2(1) = MYID
        NDEST = NPROCS - 1
        NINTS = 1 + ( NDEST-1 ) * OVHSIZE
        NREALS = 1
        CALL MPI_PACK_SIZE( NINTS, 
     &                       MPI_INTEGER, COMM,
     &                       SIZE1, IERR_MPI )
        CALL MPI_PACK_SIZE( NREALS,
     &                       MPI_DOUBLE_PRECISION, COMM,
     &                       SIZE2, IERR_MPI )
        SIZE=SIZE1+SIZE2
        CALL BUF_LOOK( BUF_LOAD, IPOS, IREQ, SIZE, IERR,
     &       IZERO, MYID2 )
        IF ( IERR .LT. 0 ) THEN
      RETURN
        ENDIF
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        WHAT = 4
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( MAX_SURF_MASTER, 1, MPI_DOUBLE_PRECISION,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR_MPI )
        IDEST = 0
        DO I = 0, NPROCS - 1
           IF ( I .ne. MYID ) THEN
              IDEST = IDEST + 1
              KEEP(267)=KEEP(267)+1
              CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     &             POSITION, MPI_PACKED, I,
     &             UPDATE_LOAD, COMM,
     &             BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     &             IERR_MPI )
           END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in MUMPS_BUF_BCAST_ARRAY'
          WRITE(*,*) ' Size,position=',SIZE,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE MUMPS_BUF_SEND_NOT_MSTR
        SUBROUTINE MUMPS_BUF_BCAST_ARRAY( BDC_MEM,
     &      COMM, MYID, NPROCS,
     &      FUTURE_NIV2,
     &      NSLAVES,
     &      LIST_SLAVES,INODE,
     &      MEM_INCREMENT, FLOPS_INCREMENT,CB_BAND, WHAT,
     &      KEEP,
     &      IERR )
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        LOGICAL BDC_MEM
        INTEGER COMM, MYID, NPROCS, NSLAVES, IERR
        INTEGER FUTURE_NIV2(NPROCS)
        INTEGER LIST_SLAVES(NSLAVES),INODE
        DOUBLE PRECISION MEM_INCREMENT(NSLAVES)
        DOUBLE PRECISION FLOPS_INCREMENT(NSLAVES)
        DOUBLE PRECISION CB_BAND(NSLAVES)
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INTEGER :: IERR_MPI
        INTEGER NDEST, NINTS, NREALS, SIZE1, SIZE2, SIZE
        INTEGER IPOS, IPOSMSG, IREQ, POSITION
        INTEGER I, IDEST, WHAT
        INTEGER IZERO
        INTEGER MYID2(1)
        PARAMETER ( IZERO=0 )
        MYID2(1)=MYID
        IERR = 0
        NDEST = 0
        DO I = 1, NPROCS
          IF ( I .NE. MYID + 1 .AND. FUTURE_NIV2(I).NE.0) THEN
            NDEST = NDEST + 1
          ENDIF
        ENDDO
        IF ( NDEST == 0 ) THEN
           RETURN
        ENDIF
        NINTS = 2 +  NSLAVES + ( NDEST - 1 ) * OVHSIZE + 1
        NREALS = NSLAVES
        IF (BDC_MEM) NREALS = NREALS + NSLAVES
        IF(WHAT.EQ.19) THEN 
           NREALS = NREALS + NSLAVES
        ENDIF
        CALL MPI_PACK_SIZE( NINTS, 
     &                       MPI_INTEGER, COMM,
     &                       SIZE1, IERR_MPI )
        CALL MPI_PACK_SIZE( NREALS, MPI_DOUBLE_PRECISION,
     &       COMM, SIZE2, IERR_MPI )
        SIZE = SIZE1+SIZE2
        CALL BUF_LOOK( BUF_LOAD, IPOS, IREQ, SIZE, IERR,
     &       IZERO, MYID2 )
        IF ( IERR .LT. 0 ) THEN
      RETURN
        ENDIF
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NSLAVES, 1, MPI_INTEGER,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( LIST_SLAVES, NSLAVES, MPI_INTEGER,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( FLOPS_INCREMENT, NSLAVES,
     &      MPI_DOUBLE_PRECISION,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR_MPI )
        IF (BDC_MEM) THEN
          CALL MPI_PACK( MEM_INCREMENT, NSLAVES,
     &      MPI_DOUBLE_PRECISION,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR_MPI )
        END IF
        IF(WHAT.EQ.19)THEN
           CALL MPI_PACK( CB_BAND, NSLAVES,
     &          MPI_DOUBLE_PRECISION,
     &          BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &          POSITION, COMM, IERR_MPI )
        ENDIF
        IDEST = 0
        DO I = 0, NPROCS - 1
        IF ( I .NE. MYID .AND. FUTURE_NIV2(I+1) .NE. 0) THEN
            IDEST = IDEST + 1
            KEEP(267)=KEEP(267)+1
            CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     &                      POSITION, MPI_PACKED, I,
     &                      UPDATE_LOAD, COMM,
     &                      BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     &                      IERR_MPI )
          END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in MUMPS_BUF_BCAST_ARRAY'
          WRITE(*,*) ' Size,position=',SIZE,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE MUMPS_BUF_BCAST_ARRAY
        SUBROUTINE MUMPS_BUF_DIST_IRECV_SIZE
     &             ( MUMPS_LBUFR_BYTES)
        IMPLICIT NONE
        INTEGER MUMPS_LBUFR_BYTES 
        SIZE_RBUF_BYTES = MUMPS_LBUFR_BYTES
        RETURN
      END SUBROUTINE MUMPS_BUF_DIST_IRECV_SIZE
#if defined(USE_GPU)
        SUBROUTINE MUMPS_BUF_MAX_OUTER_T2 ( MAX_OUTER_T2, NBSMB,
     &      NFRONT, NASS, NSLAVES, KEEP, NBmin )
        INTEGER, intent(in)  :: MAX_OUTER_T2, NBSMB
        INTEGER, intent(in)  :: NFRONT, NASS, NSLAVES, KEEP(500)
        INTEGER, intent(out) :: NBmin
        DOUBLE PRECISION     :: LARGEST_MESSAGE
         LARGEST_MESSAGE = dble(min(MAX_OUTER_T2,NASS))*
     &          min (
     &               dble(min(MAX_OUTER_T2,NASS)*NBSMB), 
     &               dble(NFRONT)
     &               )
         LARGEST_MESSAGE = LARGEST_MESSAGE*dble(KEEP(35))
         NBmin = int(LARGEST_MESSAGE/dble(SIZE_RBUF_BYTES))+1
        RETURN
        END SUBROUTINE MUMPS_BUF_MAX_OUTER_T2
#endif
      END MODULE MUMPS_BUF_COMMON
