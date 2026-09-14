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
#if ! defined(NO_SAVE_RESTORE)
      MODULE CMUMPS_SAVE_RESTORE
      USE CMUMPS_INTR_TYPES, ONLY: CMUMPS_INTR_STRUC
      USE CMUMPS_SAVE_RESTORE_FILES
      USE MUMPS_FRONT_DATA_MGT_M
      IMPLICIT NONE
      INCLUDE 'mumps_save_restore_modes.h'
      CONTAINS
      SUBROUTINE CMUMPS_REMOVE_SAVED(id)
      USE CMUMPS_OOC 
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER MASTER
      PARAMETER ( MASTER = 0 )
      TYPE (CMUMPS_STRUC) :: id
      CHARACTER(len=LEN_SAVE_FILE) :: RESTOREFILE, INFOFILE
      INTEGER :: fileunit, ierr, SIZE_INT, SIZE_INT8
      INTEGER(8) :: size_read, TOTAL_FILE_SIZE, TOTAL_STRUC_SIZE
      INTEGER :: READ_OOC_FILE_NAME_LENGTH,READ_SYM,READ_PAR,READ_NPROCS
      CHARACTER(len=LEN_SAVE_FILE) :: READ_OOC_FIRST_FILE_NAME
      CHARACTER :: READ_ARITH
      LOGICAL :: READ_INT_TYPE_64
      CHARACTER(len=23) :: READ_HASH
      LOGICAL :: FORTRAN_VERSION_OK
      LOGICAL :: SAME_OOC
      INTEGER :: ICNTL34, MAX_LENGTH, FLAG_SAME, SUM_FLAG_SAME
      TYPE (CMUMPS_STRUC) :: localid
      ierr = 0
      call CMUMPS_GET_SAVE_FILES(id,RESTOREFILE,INFOFILE)
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      CALL MUMPS_FIND_UNIT(fileunit)
      IF ( fileunit .EQ. -1 ) THEN
         id%INFO(1) = -79
         id%INFO(2) = 2
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      open(UNIT=fileunit,FILE=RESTOREFILE
#if ! defined(MUMPS_NOF2003)
     &     ,ACCESS="stream"
#endif
     &     ,STATUS='old',FORM='unformatted',IOSTAT=ierr)
      IF ( ierr.ne.0 ) THEN
        id%INFO(1) = -76
        id%INFO(2) = id%MYID
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      SIZE_INT = id%KEEP(34)
      SIZE_INT8 = id%KEEP(34)*id%KEEP(10)
      size_read = 0_8
      call MUMPS_READ_HEADER(fileunit,ierr,size_read,SIZE_INT,
     &     SIZE_INT8, TOTAL_FILE_SIZE, TOTAL_STRUC_SIZE,
     &     READ_ARITH, READ_INT_TYPE_64,
     &     READ_OOC_FILE_NAME_LENGTH, READ_OOC_FIRST_FILE_NAME,
     &     READ_HASH,READ_SYM,READ_PAR,READ_NPROCS,
     &     FORTRAN_VERSION_OK)
      close(fileunit)
      if (ierr.ne.0) THEN
         id%INFO(1) = -75
         CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &              ,id%INFO(2))
      elseif (.NOT.FORTRAN_VERSION_OK) THEN
         id%INFO(1) = -73
         id%INFO(2) = 1
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      CALL CMUMPS_CHECK_HEADER(id,.TRUE.,READ_INT_TYPE_64,
     &           READ_HASH, READ_NPROCS,
     &           READ_ARITH, READ_SYM, READ_PAR)
      IF ( id%INFO(1) .LT. 0 ) RETURN
      ICNTL34 = -99998
      IF (id%MYID.EQ.MASTER) THEN
        ICNTL34 = id%ICNTL(34)
      ENDIF
      CALL MPI_BCAST( ICNTL34, 1, MPI_INTEGER, MASTER, id%COMM, ierr )
      CALL MPI_BCAST( READ_SYM, 1, MPI_INTEGER, MASTER, id%COMM, ierr )
      CALL CMUMPS_CHECK_FILE_NAME(id, READ_OOC_FILE_NAME_LENGTH,
     &           READ_OOC_FIRST_FILE_NAME, SAME_OOC)
      CALL MPI_ALLREDUCE(READ_OOC_FILE_NAME_LENGTH,MAX_LENGTH,1,
     &                   MPI_INTEGER,MPI_MAX,id%COMM,ierr)
      IF (MAX_LENGTH.NE.-999) THEN
        FLAG_SAME = 0
        IF (SAME_OOC) THEN
          FLAG_SAME = 1
        ENDIF
        CALL MPI_ALLREDUCE(FLAG_SAME,SUM_FLAG_SAME,1,
     &                   MPI_INTEGER,MPI_SUM,id%COMM,ierr)
        IF (SUM_FLAG_SAME.NE.0) THEN
          IF (ICNTL34 .EQ. 1) THEN
            id%ASSOCIATED_OOC_FILES = .TRUE.
          ELSE
            id%ASSOCIATED_OOC_FILES = .FALSE.
          ENDIF
        ELSE
          IF (ICNTL34 .NE. 1) THEN
            localid%COMM = id%COMM
            localid%INFO(1) = 0
            localid%ICNTL(1) = id%ICNTL(1) 
            localid%MYID = id%MYID
            localid%NPROCS = id%NPROCS
            localid%KEEP(10) = id%KEEP(10)
            localid%SAVE_PREFIX = id%SAVE_PREFIX
            localid%SAVE_DIR = id%SAVE_DIR
            call CMUMPS_RESTORE_OOC(localid)
            IF ( localid%INFO(1) .EQ. 0 ) THEN
              localid%ASSOCIATED_OOC_FILES = .FALSE.
              IF (READ_OOC_FILE_NAME_LENGTH.NE.-999) THEN
                call CMUMPS_OOC_CLEAN_FILES(localid,ierr)
                IF ( ierr.ne.0 ) THEN
                  id%INFO(1) = -90
                  id%INFO(2) = id%MYID
                ENDIF
              ENDIF
            ENDIF
            CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &           id%COMM, id%MYID )
            IF ( id%INFO(1) .LT. 0 ) RETURN 
          ENDIF
        ENDIF
      ENDIF
      call MUMPS_CLEAN_SAVED_DATA(id%MYID,ierr,RESTOREFILE,INFOFILE)
      IF (ierr.eq.-79) THEN
        id%INFO(1) = -79
        id%INFO(2) = 2
      ELSE IF (ierr.ne.0) THEN
        id%INFO(1) = -76
        id%INFO(2) = id%MYID
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      RETURN
      END SUBROUTINE CMUMPS_REMOVE_SAVED
      SUBROUTINE CMUMPS_RESTORE_OOC(localid)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      TYPE (CMUMPS_STRUC) :: localid
      INTEGER::ierr,IN,NBVARIABLES,NBVARIABLES_ROOTC
      CHARACTER(len=LEN_SAVE_FILE):: restore_file_ooc,INFO_FILE
      INTEGER(8),allocatable, dimension(:)::SIZE_VARIABLES
      INTEGER(8),allocatable, dimension(:)::SIZE_VARIABLES_ROOTC
      INTEGER,allocatable, dimension(:)::SIZE_GEST
      INTEGER,allocatable, dimension(:)::SIZE_GEST_ROOTC
      INTEGER:: INFO1,INFO2,INFOG1,INFOG2,allocok
      INTEGER(8) :: TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
      TYPE (CMUMPS_INTR_STRUC) :: localidintr
      NBVARIABLES=194
      NBVARIABLES_ROOTC=20
      allocate(SIZE_VARIABLES(NBVARIABLES),  stat=allocok)
      if (allocok .GT. 0) THEN
       localid%INFO(1)  =-13
       localid%INFO(2) = NBVARIABLES
      endif
      CALL MUMPS_PROPINFO( localid%ICNTL(1), localid%INFO(1),
     &     localid%COMM, localid%MYID )
      IF ( localid%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_VARIABLES_ROOTC(NBVARIABLES_ROOTC),  stat=allocok)
      if (allocok .GT. 0) THEN
       localid%INFO(1)  =-13
       localid%INFO(2) = NBVARIABLES_ROOTC
      endif
      CALL MUMPS_PROPINFO( localid%ICNTL(1), localid%INFO(1),
     &     localid%COMM, localid%MYID )
      IF ( localid%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_GEST(NBVARIABLES),  stat=allocok)
      if (allocok .GT. 0) THEN
       localid%INFO(1)  =-13
       localid%INFO(2) = NBVARIABLES
      endif
      CALL MUMPS_PROPINFO( localid%ICNTL(1), localid%INFO(1),
     &     localid%COMM, localid%MYID )
      IF ( localid%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_GEST_ROOTC(NBVARIABLES_ROOTC), stat=allocok)
      if (allocok .GT. 0) THEN
       localid%INFO(1)  =-13
       localid%INFO(2) = NBVARIABLES_ROOTC
      endif
      CALL MUMPS_PROPINFO( localid%ICNTL(1), localid%INFO(1),
     &     localid%COMM, localid%MYID )
      IF ( localid%INFO(1) .LT. 0 ) RETURN
      SIZE_VARIABLES(:)=0_8
      SIZE_VARIABLES_ROOTC(:)=0_8
      SIZE_GEST(:)=0
      SIZE_GEST_ROOTC(:)=0
      TOTAL_FILE_SIZE=0_8
      TOTAL_STRUC_SIZE=0_8
      INFO1 = -999
      INFO2 = -999
      INFOG1 = -999
      INFOG2 = -999
      CALL CMUMPS_GET_SAVE_FILES(localid,restore_file_ooc,INFO_FILE)
      IF ( localid%INFO(1) .LT. 0 ) RETURN
      CALL MUMPS_FIND_UNIT(IN)
      IF ( IN .EQ. -1 ) THEN
         localid%INFO(1) = -79
         localid%INFO(2) = 2
      ENDIF
      CALL MUMPS_PROPINFO( localid%ICNTL(1), localid%INFO(1),
     &     localid%COMM, localid%MYID )
      IF ( localid%INFO(1) .LT. 0 ) RETURN
      open(UNIT=IN,FILE=restore_file_ooc
#if ! defined(MUMPS_NOF2003)
     &     ,ACCESS="stream"
#endif
     &     ,STATUS='old',form='unformatted',iostat=ierr)
      if(ierr.ne.0) THEN
         localid%INFO(1) = -74
         localid%INFO(2) = localid%MYID
      endif
      CALL MUMPS_PROPINFO( localid%ICNTL(1), localid%INFO(1),
     &     localid%COMM, localid%MYID )
      IF ( localid%INFO(1) .LT. 0 ) RETURN
      CALL CMUMPS_SAVE_RESTORE_STRUCTURE(localid,localidintr,IN
     &     ,restore_ooc_mode
     &     ,NBVARIABLES,SIZE_VARIABLES,SIZE_GEST
     &     ,NBVARIABLES_ROOTC,SIZE_VARIABLES_ROOTC,SIZE_GEST_ROOTC
     &     ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &     ,INFO1,INFO2,INFOG1,INFOG2)
      CLOSE(IN)
      deallocate(SIZE_VARIABLES,SIZE_VARIABLES_ROOTC)
      deallocate(SIZE_GEST,SIZE_GEST_ROOTC)
      RETURN
      END SUBROUTINE CMUMPS_RESTORE_OOC
      SUBROUTINE CMUMPS_COMPUTE_MEMORY_SAVE(id,idintr,
     &     TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      TYPE (CMUMPS_STRUC)      :: id
      TYPE (CMUMPS_INTR_STRUC) :: idintr
      INTEGER::NBVARIABLES,NBVARIABLES_ROOTC
      INTEGER(8),allocatable, dimension(:)::SIZE_VARIABLES
      INTEGER(8),allocatable, dimension(:)::SIZE_VARIABLES_ROOTC
      INTEGER,allocatable, dimension(:)::SIZE_GEST
      INTEGER,allocatable, dimension(:)::SIZE_GEST_ROOTC
      INTEGER :: INFO1,INFO2,INFOG1,INFOG2,allocok
      INTEGER(8) :: TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
      NBVARIABLES=194
      NBVARIABLES_ROOTC=20
      allocate(SIZE_VARIABLES(NBVARIABLES),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_VARIABLES_ROOTC(NBVARIABLES_ROOTC), stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES_ROOTC
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_GEST(NBVARIABLES),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_GEST_ROOTC(NBVARIABLES_ROOTC), stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES_ROOTC
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      SIZE_VARIABLES(:)=0_8
      SIZE_VARIABLES_ROOTC(:)=0_8
      SIZE_GEST(:)=0
      SIZE_GEST_ROOTC(:)=0
      TOTAL_FILE_SIZE=0_8
      TOTAL_STRUC_SIZE=0_8
      INFO1 = -999
      INFO2 = -999
      INFOG1 = -999
      INFOG2 = -999
      CALL CMUMPS_SAVE_RESTORE_STRUCTURE(id,idintr,0,memory_save_mode
     &     ,NBVARIABLES,SIZE_VARIABLES,SIZE_GEST
     &     ,NBVARIABLES_ROOTC,SIZE_VARIABLES_ROOTC,SIZE_GEST_ROOTC
     &     ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &     ,INFO1,INFO2,INFOG1,INFOG2)
      deallocate(SIZE_VARIABLES,SIZE_VARIABLES_ROOTC)
      deallocate(SIZE_GEST,SIZE_GEST_ROOTC)
      RETURN
      END SUBROUTINE CMUMPS_COMPUTE_MEMORY_SAVE
      SUBROUTINE CMUMPS_SAVE(id,idintr)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      TYPE (CMUMPS_STRUC)      :: id
      TYPE (CMUMPS_INTR_STRUC) :: idintr
      INTEGER::ierr,OUT,NBVARIABLES,NBVARIABLES_ROOTC,OUTINFO
      CHARACTER(len=LEN_SAVE_FILE):: SAVE_FILE,INFO_FILE
      LOGICAL:: SAVE_FILE_exist,INFO_FILE_exist
      INTEGER(8),allocatable, dimension(:)::SIZE_VARIABLES
      INTEGER(8),allocatable, dimension(:)::SIZE_VARIABLES_ROOTC
      INTEGER,allocatable, dimension(:)::SIZE_GEST
      INTEGER,allocatable, dimension(:)::SIZE_GEST_ROOTC
      INTEGER:: INFO1,INFO2,INFOG1,INFOG2,MPG
      INTEGER(8) :: TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
      LOGICAL :: PROKG
      INTEGER :: I,J,K,H,allocok
      CHARACTER(len=1) :: TMP_OOC_NAMES(350)
      INFO1 = id%INFO(1)
      INFO2 = id%INFO(2)
      INFOG1 = id%INFO(1)
      INFOG2 = id%INFO(1)
      id%INFO(1)=0
      id%INFO(2)=0
      id%INFOG(1)=0
      id%INFOG(2)=0
      MPG= id%ICNTL(3)
      PROKG   = ( MPG .GT. 0 .and. id%MYID .eq. 0 )
      NBVARIABLES=194
      NBVARIABLES_ROOTC=20
      allocate(SIZE_VARIABLES(NBVARIABLES),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_VARIABLES_ROOTC(NBVARIABLES_ROOTC), stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES_ROOTC
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_GEST(NBVARIABLES),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_GEST_ROOTC(NBVARIABLES_ROOTC),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES_ROOTC
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      SIZE_VARIABLES(:)=0_8
      SIZE_VARIABLES_ROOTC(:)=0_8
      SIZE_GEST(:)=0
      SIZE_GEST_ROOTC(:)=0
      TOTAL_FILE_SIZE=0_8
      TOTAL_STRUC_SIZE=0_8
      TMP_OOC_NAMES(:)="?"
      CALL CMUMPS_SAVE_RESTORE_STRUCTURE(id,idintr,0,memory_save_mode
     &     ,NBVARIABLES,SIZE_VARIABLES,SIZE_GEST
     &     ,NBVARIABLES_ROOTC,SIZE_VARIABLES_ROOTC,SIZE_GEST_ROOTC
     &     ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &     ,INFO1,INFO2,INFOG1,INFOG2)
      CALL CMUMPS_GET_SAVE_FILES(id,SAVE_FILE,INFO_FILE)
      IF ( id%INFO(1) .LT. 0 ) RETURN
      inquire(FILE=SAVE_FILE, EXIST=SAVE_FILE_exist)
      IF(SAVE_FILE_exist) THEN
         id%INFO(1) = -70
         id%INFO(2) = 0
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      CALL MUMPS_FIND_UNIT(OUT)
      IF ( OUT .EQ. -1 ) THEN
         id%INFO(1) = -79
         id%INFO(2) = 2
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      open(UNIT=OUT,FILE=SAVE_FILE
#if ! defined(MUMPS_NOF2003)
     &     ,ACCESS="stream"
#endif
     &     ,STATUS='new',form='unformatted',iostat=ierr)
      if(ierr.ne.0) THEN
         id%INFO(1) = -71
         id%INFO(2) = 0
      endif  
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      inquire(FILE=INFO_FILE, EXIST=INFO_FILE_exist)
      IF(INFO_FILE_exist) THEN
         id%INFO(1) = -70
         id%INFO(2) = 0
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      CALL MUMPS_FIND_UNIT(OUTINFO)
      IF ( OUTINFO .EQ. -1 ) THEN
         id%INFO(1) = -79
         id%INFO(2) = 2
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      open(UNIT=OUTINFO,FILE=INFO_FILE,STATUS='new',iostat=ierr)
      if(ierr.ne.0) THEN
         id%INFO(1) = -71
         id%INFO(2) = 0
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      CALL CMUMPS_SAVE_RESTORE_STRUCTURE(id,idintr,OUT,save_mode
     &     ,NBVARIABLES,SIZE_VARIABLES,SIZE_GEST
     &     ,NBVARIABLES_ROOTC,SIZE_VARIABLES_ROOTC,SIZE_GEST_ROOTC
     &     ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &     ,INFO1,INFO2,INFOG1,INFOG2)
      if (id%INFO(1).EQ.0) then
        id%INFO(1)=INFO1
        id%INFO(2)=INFO2
        id%INFOG(1)=INFOG1
        id%INFOG(2)=INFOG2
        CLOSE(OUT)
        if(id%INFOG(1).NE.0) then
           if (PROKG) THEN
             write(MPG,*) "Warning: "
     &           ,"saved instance has negative INFO(1):"
     &           , id%INFOG(1)
           endif
        endif
        IF(PROKG) THEN
           write(MPG,*) "Save done successfully"
           IF(id%KEEP(201).EQ.1) THEN
              K=1
              write(MPG,*) "The corresponding OOC files are:"
              DO I=1,id%OOC_NB_FILE_TYPE
                 DO J=1,id%OOC_NB_FILES(I)
                    DO H=1,id%OOC_FILE_NAME_LENGTH(K)-2
                       TMP_OOC_NAMES(H)=id%OOC_FILE_NAMES(K,H)
                    ENDDO
                    write(MPG,*)
     &                   TMP_OOC_NAMES(1:id%OOC_FILE_NAME_LENGTH(K)-2)
                    K=K+1
                 ENDDO
              ENDDO
           ENDIF
        ENDIF
        write(OUTINFO,*) "Save done by CMUMPS ",
     &       trim(adjustl(id%VERSION_NUMBER)),
     &       " after JOB=",id%KEEP(40)+456789,
     &       " With SYM, PAR =",id%KEEP(50),id%KEEP(46)
        write(OUTINFO,*) "On ",id%NPROCS," processes"
        if((id%ICNTL(18).EQ.0).AND.(id%ICNTL(5).EQ.0)) then
           write(OUTINFO,*) "with N, NNZ ", id%N, id%NNZ
        elseif((id%ICNTL(18).EQ.1).AND.(id%ICNTL(5).EQ. 0)) then
           write(OUTINFO,*) "with N, NNZ_loc=", id%N, id%NNZ_loc
        elseif((id%ICNTL(18).EQ.0).AND.(id%ICNTL(5).EQ. 1)) then
           write(OUTINFO,*) "with N, NELT=", id%N, id%NELT
        endif
        IF(id%KEEP(10).EQ.1) THEN
           write(OUTINFO,*) "With a default integer size of 64 bits"
        ELSE
           write(OUTINFO,*) "With a default integer size of 32 bits"
        ENDIF
#if defined(MUMPS_NOF2003)
        write(OUTINFO,*) "Using MUMPS_NOF2003"
#endif        
        write(OUTINFO,*) ''
        write(OUTINFO,*) "The corresponding save file is:"
        write(OUTINFO,*) trim(adjustl(SAVE_FILE))
        write(OUTINFO,*) "of size",TOTAL_FILE_SIZE, " Bytes"
        IF(id%KEEP(201).EQ.1) THEN
           write(OUTINFO,*) ''
           write(OUTINFO,*) "The corresponding OOC files are:"
           K=1
           DO I=1,id%OOC_NB_FILE_TYPE
              DO J=1,id%OOC_NB_FILES(I)
                 DO H=1,id%OOC_FILE_NAME_LENGTH(K)-2
                    TMP_OOC_NAMES(H)=id%OOC_FILE_NAMES(K,H)
                 ENDDO
                 write(OUTINFO,*)
     &                TMP_OOC_NAMES(1:id%OOC_FILE_NAME_LENGTH(K)-2)
                 K=K+1
              ENDDO
           ENDDO
        ENDIF
        CLOSE(OUTINFO)
      else
         CLOSE(OUT,STATUS='delete')
         CLOSE(OUTINFO,STATUS='delete')         
      endif
      deallocate(SIZE_VARIABLES,SIZE_VARIABLES_ROOTC)
      deallocate(SIZE_GEST,SIZE_GEST_ROOTC)
      if (id%KEEP(201) .GT. 0) THEN
        id%ASSOCIATED_OOC_FILES=.TRUE.
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_SAVE
      SUBROUTINE CMUMPS_RESTORE(id,idintr)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER::ierr,IN,NBVARIABLES,NBVARIABLES_ROOTC
      CHARACTER(len=LEN_SAVE_FILE):: restore_file,INFO_FILE
      INTEGER(8),allocatable, dimension(:)::SIZE_VARIABLES
      INTEGER(8),allocatable, dimension(:)::SIZE_VARIABLES_ROOTC
      INTEGER,allocatable, dimension(:)::SIZE_GEST
      INTEGER,allocatable, dimension(:)::SIZE_GEST_ROOTC
      INTEGER:: INFO1,INFO2,INFOG1,INFOG2,MPG,MP,JOB
      INTEGER(8) :: TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
      LOGICAL :: PROKG
      INTEGER :: I,J,K,H,allocok
      CHARACTER(len=1) :: TMP_OOC_NAMES(350)
      TYPE (CMUMPS_STRUC)      :: id
      TYPE (CMUMPS_INTR_STRUC) :: idintr
      NBVARIABLES=194
      NBVARIABLES_ROOTC=20
      allocate(SIZE_VARIABLES(NBVARIABLES),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_VARIABLES_ROOTC(NBVARIABLES_ROOTC), stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES_ROOTC
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_GEST(NBVARIABLES),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_GEST_ROOTC(NBVARIABLES_ROOTC), stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES_ROOTC
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      SIZE_VARIABLES(:)=0_8
      SIZE_VARIABLES_ROOTC(:)=0_8
      SIZE_GEST(:)=0
      SIZE_GEST_ROOTC(:)=0
      TOTAL_FILE_SIZE=0_8
      TOTAL_STRUC_SIZE=0_8
      TMP_OOC_NAMES(:)="?"
      INFO1 = -999
      INFO2 = -999
      INFOG1 = -999
      INFOG2 = -999
      CALL CMUMPS_GET_SAVE_FILES(id,restore_file,INFO_FILE)
      IF ( id%INFO(1) .LT. 0 ) RETURN
      CALL MUMPS_FIND_UNIT(IN)
      IF ( IN .EQ. -1 ) THEN
         id%INFO(1) = -79
         id%INFO(2) = 2
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      open(UNIT=IN,FILE=restore_file
#if ! defined(MUMPS_NOF2003)
     &     ,ACCESS="stream"
#endif
     &     ,STATUS='old',form='unformatted',iostat=ierr)
      if(ierr.ne.0) THEN
         id%INFO(1) = -74
         id%INFO(2) = id%MYID
      endif  
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      MP= id%ICNTL(2)
      MPG= id%ICNTL(3)
      CALL CMUMPS_SAVE_RESTORE_STRUCTURE(id,idintr,IN,restore_mode
     &     ,NBVARIABLES,SIZE_VARIABLES,SIZE_GEST
     &     ,NBVARIABLES_ROOTC,SIZE_VARIABLES_ROOTC,SIZE_GEST_ROOTC
     &     ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &     ,INFO1,INFO2,INFOG1,INFOG2)
      PROKG   = ( MPG .GT. 0 .and. id%MYID .eq. 0 )
      if(id%INFO(1).EQ.0) then
         id%INFO(1)=INFO1
         id%INFO(2)=INFO2
         id%INFOG(1)=INFOG1
         id%INFOG(2)=INFOG2
         if(id%INFOG(1).NE.0) then
            write(MPG,*) "Warning: "
     &           ,"restored instance has negative INFOG(1):"
     &           , id%INFOG(1)
         endif
         if(MP.GT.0) then
            JOB=id%KEEP(40)+456789
            write(MP,*) "Restore done successfully"
            write(MP,*) "From file ",trim(adjustl(restore_file))
            if((id%ICNTL(18).EQ.0).AND.(id%ICNTL(5).EQ.0)) then
               write(MP,*) "with JOB, N, NNZ ",JOB, id%N,id%NNZ
            elseif((id%ICNTL(18).EQ.1).AND.(id%ICNTL(5).EQ. 0)) then
               write(MP,*) "with JOB, N, NNZ_loc=", JOB, id%N,
     &              id%NNZ_loc
            elseif((id%ICNTL(18).EQ.0).AND.(id%ICNTL(5).EQ. 1)) then
               write(MP,*) "with JOB, N, NELT=", JOB, id%N, id%NELT
            endif
         endif
         IF(PROKG) THEN
            IF(id%KEEP(201).EQ.1) THEN
               K=1
               write(MPG,*) "The corresponding OOC files are:"
               DO I=1,id%OOC_NB_FILE_TYPE
                  DO J=1,id%OOC_NB_FILES(I)
                     DO H=1,id%OOC_FILE_NAME_LENGTH(K)-2
                        TMP_OOC_NAMES(H)=id%OOC_FILE_NAMES(K,H)
                     ENDDO
                     write(MPG,*)
     &                    TMP_OOC_NAMES(1:id%OOC_FILE_NAME_LENGTH(K)-2)
                     K=K+1
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      else
         idintr%root%gridinit_done=.FALSE.
         id%KEEP(140)=1
      endif
      CLOSE(IN)
      deallocate(SIZE_VARIABLES,SIZE_VARIABLES_ROOTC)
      deallocate(SIZE_GEST,SIZE_GEST_ROOTC)
      if (id%KEEP(201) .GT. 0) THEN
        id%ASSOCIATED_OOC_FILES=.TRUE.
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_RESTORE
      SUBROUTINE CMUMPS_SAVE_RESTORE_STRUCTURE(id,idintr,unit,mode
     &     ,NBVARIABLES,SIZE_VARIABLES,SIZE_GEST
     &     ,NBVARIABLES_ROOTC,SIZE_VARIABLES_ROOTC,SIZE_GEST_ROOTC
     &     ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &     ,INFO1,INFO2,INFOG1,INFOG2)
      USE CMUMPS_FACSOL_L0OMP_M, ONLY : CMUMPS_SAVE_RESTORE_L0FACARRAY
      USE CMUMPS_LR_DATA_M, ONLY: CMUMPS_SAVE_RESTORE_BLR
#if defined(SINGLEINDOUBLE)
#endif
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER, intent(in) ::unit,NBVARIABLES,NBVARIABLES_ROOTC
      INTEGER, intent(in) :: mode
      INTEGER(8),dimension(NBVARIABLES)::SIZE_VARIABLES
      INTEGER(8),dimension(NBVARIABLES_ROOTC)::SIZE_VARIABLES_ROOTC
      INTEGER,dimension(NBVARIABLES)::SIZE_GEST
      INTEGER,dimension(NBVARIABLES_ROOTC)::SIZE_GEST_ROOTC
      INTEGER(8) :: TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
      INTEGER:: INFO1,INFO2,INFOG1,INFOG2
      INTEGER:: j,i1,i2,err,ierr
      CHARACTER :: ARITH,READ_ARITH
      INTEGER(8) :: size_written,gest_size,WRITTEN_STRUC_SIZE
      INTEGER:: SIZE_INT, SIZE_INT8, SIZE_RL_OR_DBL, SIZE_ARITH_DEP
      INTEGER:: SIZE_DOUBLE_PRECISION, SIZE_LOGICAL, SIZE_CHARACTER
      INTEGER:: READ_NPROCS, READ_PAR, READ_SYM
      INTEGER,dimension(NBVARIABLES)::NbRecords
      INTEGER,dimension(NBVARIABLES_ROOTC)::NbRecords_ROOTC
      INTEGER:: size_array1,size_array2,dummy,allocok
      INTEGER(8):: size_array_INT8_1,size_array_INT8_2
      LOGICAL:: INT_TYPE_64, READ_INT_TYPE_64, CALL_SAVE_RESTORE_BLR
      INTEGER:: tot_NbRecords,NbSubRecords
      INTEGER(8):: size_read,size_allocated
      INTEGER(8),dimension(NBVARIABLES)::DIFF_SIZE_ALLOC_READ
      INTEGER(8),dimension(NBVARIABLES_ROOTC)::
     &                             DIFF_SIZE_ALLOC_READ_ROOTC
      INTEGER::READ_OOC_FILE_NAME_LENGTH
      CHARACTER(len=LEN_SAVE_FILE):: READ_OOC_FIRST_FILE_NAME
      INTEGER,dimension(4)::OOC_INDICES
      CHARACTER(len=8) :: date
      CHARACTER(len=10) :: time
      CHARACTER(len=5) :: zone
      INTEGER,dimension(8):: values
      CHARACTER(len=23) :: hash,READ_HASH
      LOGICAL:: BASIC_CHECK
      LOGICAL :: FORTRAN_VERSION_OK
      CHARACTER(len=1) :: TMP_OOC_NAMES(350)
      INTEGER(8)::SIZE_VARIABLES_BLR,SIZE_VARIABLES_FRONT_DATA,
     &            SIZE_VARIABLES_L0FAC
      INTEGER :: SIZE_GEST_ROOTA
      INTEGER(8) :: SIZE_VARIABLES_ROOTA
      INTEGER::SIZE_GEST_BLR,SIZE_GEST_FRONT_DATA,SIZE_GEST_L0FAC
#if defined(SINGLEINDOUBLE)
#endif
      INTEGER :: KEEP410_SAVE, KEEP411_SAVE 
      INTEGER(8) :: KEEP883_SAVE, KEEP884_SAVE 
      INTEGER(4) :: I4
      LOGICAL :: IS_SYMMETRIC
      TYPE (CMUMPS_STRUC) :: id
      TYPE (CMUMPS_INTR_STRUC) :: idintr
      INTEGER, PARAMETER :: S_pad13=194
      INTEGER, PARAMETER :: S_ASSOCIATED_OOC_FILES=193
      INTEGER, PARAMETER :: S_pad12=192
      INTEGER, PARAMETER :: S_NB_SINGULAR_VALUES=191
      INTEGER, PARAMETER :: S_SINGULAR_VALUES=190
      INTEGER, PARAMETER :: S_MTKO_PROCS_MAP=189
      INTEGER, PARAMETER :: S_L0_OMP_MAPPING=188
      INTEGER, PARAMETER :: S_PTR_LEAFS_L0_OMP=187
      INTEGER, PARAMETER :: S_PERM_L0_OMP=186
      INTEGER, PARAMETER :: S_VIRT_L0_OMP_MAPPING=185
      INTEGER, PARAMETER :: S_VIRT_L0_OMP=184
      INTEGER, PARAMETER :: S_PHYS_L0_OMP=183
      INTEGER, PARAMETER :: S_IPOOL_A_L0_OMP=182
      INTEGER, PARAMETER :: S_IPOOL_B_L0_OMP=181
      INTEGER, PARAMETER :: S_I8_L0_OMP=180
      INTEGER, PARAMETER :: S_I4_L0_OMP=179
      INTEGER, PARAMETER :: S_THREAD_LA=178
      INTEGER, PARAMETER :: S_LL0_OMP_FACTORS=177
      INTEGER, PARAMETER :: S_LL0_OMP_MAPPING=176
      INTEGER, PARAMETER :: S_L_VIRT_L0_OMP=175
      INTEGER, PARAMETER :: S_L_PHYS_L0_OMP=174
      INTEGER, PARAMETER :: S_LPOOL_B_L0_OMP=173
      INTEGER, PARAMETER :: S_LPOOL_A_L0_OMP=172
      INTEGER, PARAMETER :: S_BLRARRAY_ENCODING=171
      INTEGER, PARAMETER :: S_FDM_F_ENCODING=170
      INTEGER, PARAMETER :: S_pad11=169
      INTEGER, PARAMETER :: S_NBGRP=168
      INTEGER, PARAMETER :: S_LRGROUPS=167
      INTEGER, PARAMETER :: S_INTR_ENCODING=166
      INTEGER, PARAMETER :: S_WORKING=165
      INTEGER, PARAMETER :: S_IPTR_WORKING=164
      INTEGER, PARAMETER :: S_pad10=163
      INTEGER, PARAMETER :: S_SUP_PROC=162
      INTEGER, PARAMETER :: S_PIVNUL_LIST=161
      INTEGER, PARAMETER :: S_OOC_FILE_NAMES=160
      INTEGER, PARAMETER :: S_OOC_FILE_NAME_LENGTH=159
      INTEGER, PARAMETER :: S_pad9=158
      INTEGER, PARAMETER :: S_OOC_NB_FILE_TYPE=157
      INTEGER, PARAMETER :: S_OOC_NB_FILES=156
      INTEGER, PARAMETER :: S_OOC_TOTAL_NB_NODES=155
      INTEGER, PARAMETER :: S_OOC_VADDR=154
      INTEGER, PARAMETER :: S_OOC_SIZE_OF_BLOCK=153
      INTEGER, PARAMETER :: S_OOC_INODE_SEQUENCE=152
      INTEGER, PARAMETER :: S_OOC_MAX_NB_NODES_FOR_ZONE=151
      INTEGER, PARAMETER :: S_INSTANCE_NUMBER=150
      INTEGER, PARAMETER :: S_CB_SON_SIZE=149
      INTEGER, PARAMETER :: S_DKEEP=148
      INTEGER, PARAMETER :: S_LWK_USER=147
      INTEGER, PARAMETER :: S_NBSA_LOCAL=146
      INTEGER, PARAMETER :: S_WK_USER=145
      INTEGER, PARAMETER :: S_CROIX_MANU=144
      INTEGER, PARAMETER :: S_SCHED_SBTR=143
      INTEGER, PARAMETER :: S_SCHED_GRP=142
      INTEGER, PARAMETER :: S_SCHED_DEP=141
      INTEGER, PARAMETER :: S_SBTR_ID=140
      INTEGER, PARAMETER :: S_DEPTH_FIRST_SEQ=139
      INTEGER, PARAMETER :: S_DEPTH_FIRST=138
      INTEGER, PARAMETER :: S_MY_NB_LEAF=137
      INTEGER, PARAMETER :: S_MY_FIRST_LEAF=136
      INTEGER, PARAMETER :: S_MY_ROOT_SBTR=135
      INTEGER, PARAMETER :: S_COST_TRAV=134
      INTEGER, PARAMETER :: S_MEM_SUBTREE=133
      INTEGER, PARAMETER :: S_RHSINTR=132
      INTEGER, PARAMETER :: S_GLOB2LOC_SOL=131
      INTEGER, PARAMETER :: S_pad8=130
      INTEGER, PARAMETER :: S_GLOB2LOC_SOL_ALLOC=129
      INTEGER, PARAMETER :: S_GLOB2LOC_RHS=128
      INTEGER, PARAMETER :: S_MEM_DIST=127
      INTEGER, PARAMETER :: S_I_AM_CAND=126
      INTEGER, PARAMETER :: S_TAB_POS_IN_PERE=125
      INTEGER, PARAMETER :: S_FUTURE_NIV2=124
      INTEGER, PARAMETER :: S_ISTEP_TO_INIV2=123
      INTEGER, PARAMETER :: S_CANDIDATES=122
      INTEGER, PARAMETER :: S_ELTPROC=121
      INTEGER, PARAMETER :: S_LELTVAR=120
      INTEGER, PARAMETER :: S_NELT_loc=119
      INTEGER, PARAMETER :: S_PROCNODE=118
      INTEGER, PARAMETER :: S_LPS=117
      INTEGER, PARAMETER :: S_S=116
      INTEGER, PARAMETER :: S_PTRFAC=115
      INTEGER, PARAMETER :: S_PTLUST_S=114
      INTEGER, PARAMETER :: S_Step2node=113
      INTEGER, PARAMETER :: S_PROCNODE_STEPS=112
      INTEGER, PARAMETER :: S_NA=111
      INTEGER, PARAMETER :: S_PTRDEBARR=110
      INTEGER, PARAMETER :: S_NINROWARR=109
      INTEGER, PARAMETER :: S_NINCOLARR=108
      INTEGER, PARAMETER :: S_PTR8ARR=107
      INTEGER, PARAMETER :: S_PTRAR=106
      INTEGER, PARAMETER :: S_FRTELT=105
      INTEGER, PARAMETER :: S_FRTPTR=104
      INTEGER, PARAMETER :: S_FILS=103
      INTEGER, PARAMETER :: S_DAD_STEPS=102
      INTEGER, PARAMETER :: S_FRERE_STEPS=101
      INTEGER, PARAMETER :: S_ND_STEPS=100
      INTEGER, PARAMETER :: S_NE_STEPS=99
      INTEGER, PARAMETER :: S_STEP=98
      INTEGER, PARAMETER :: S_NBSA=97
      INTEGER, PARAMETER :: S_LNA=96
      INTEGER, PARAMETER :: S_KEEP=95
      INTEGER, PARAMETER :: S_IS=94
      INTEGER, PARAMETER :: S_ASS_IRECV=93
      INTEGER, PARAMETER :: S_NSLAVES=92
      INTEGER, PARAMETER :: S_NPROCS=91
      INTEGER, PARAMETER :: S_MYID=90
      INTEGER, PARAMETER :: S_COMM_LOAD=89
      INTEGER, PARAMETER :: S_MYID_NODES=88
      INTEGER, PARAMETER :: S_COMM_NODES=87
      INTEGER, PARAMETER :: S_INST_Number=86
      INTEGER, PARAMETER :: S_MAX_SURF_MASTER=85
      INTEGER, PARAMETER :: S_KEEP8=84
      INTEGER, PARAMETER :: S_pad7=83
      INTEGER, PARAMETER :: S_SAVE_PREFIX=82
      INTEGER, PARAMETER :: S_SAVE_DIR=81
      INTEGER, PARAMETER :: S_WRITE_PROBLEM=80
      INTEGER, PARAMETER :: S_OOC_PREFIX=79
      INTEGER, PARAMETER :: S_OOC_TMPDIR=78
      INTEGER, PARAMETER :: S_VERSION_NUMBER=77
      INTEGER, PARAMETER :: S_MAPPING=76
      INTEGER, PARAMETER :: S_LISTVAR_SCHUR=75
      INTEGER, PARAMETER :: S_SCHUR_CINTERFACE=74
      INTEGER, PARAMETER :: S_SCHUR=73
      INTEGER, PARAMETER :: S_SIZE_SCHUR=72
      INTEGER, PARAMETER :: S_SCHUR_LLD=71
      INTEGER, PARAMETER :: S_SCHUR_NLOC=70
      INTEGER, PARAMETER :: S_SCHUR_MLOC=69
      INTEGER, PARAMETER :: S_NBLOCK=68
      INTEGER, PARAMETER :: S_MBLOCK=67
      INTEGER, PARAMETER :: S_NPCOL=66
      INTEGER, PARAMETER :: S_NPROW=65
      INTEGER, PARAMETER :: S_UNS_PERM=64
      INTEGER, PARAMETER :: S_SYM_PERM=63
      INTEGER, PARAMETER :: S_METIS_OPTIONS=62
      INTEGER, PARAMETER :: S_RINFOG=61
      INTEGER, PARAMETER :: S_RINFO=60
      INTEGER, PARAMETER :: S_CNTL=59
      INTEGER, PARAMETER :: S_COST_SUBTREES=58
      INTEGER, PARAMETER :: S_INFOG=57
      INTEGER, PARAMETER :: S_INFO=56
      INTEGER, PARAMETER :: S_ICNTL=55
      INTEGER, PARAMETER :: S_pad6=54
      INTEGER, PARAMETER :: S_LD_RHSINTR=53
      INTEGER, PARAMETER :: S_NSOL_loc=52
      INTEGER, PARAMETER :: S_LSOL_loc=51
      INTEGER, PARAMETER :: S_LREDRHS=50
      INTEGER, PARAMETER :: S_LRHS_loc=49
      INTEGER, PARAMETER :: S_Nloc_RHS=48
      INTEGER, PARAMETER :: S_NZ_RHS=47
      INTEGER, PARAMETER :: S_NRHS=46
      INTEGER, PARAMETER :: S_LRHS=45
      INTEGER, PARAMETER :: S_IRHS_loc=44
      INTEGER, PARAMETER :: S_ISOL_loc=43
      INTEGER, PARAMETER :: S_IRHS_PTR=42
      INTEGER, PARAMETER :: S_IRHS_SPARSE=41
      INTEGER, PARAMETER :: S_RHS_loc=40
      INTEGER, PARAMETER :: S_SOL_loc=39
      INTEGER, PARAMETER :: S_RHS_SPARSE=38
      INTEGER, PARAMETER :: S_REDRHS=37
      INTEGER, PARAMETER :: S_RHS=36
      INTEGER, PARAMETER :: S_BLKVAR=35
      INTEGER, PARAMETER :: S_BLKPTR=34
      INTEGER, PARAMETER :: S_pad5=33
      INTEGER, PARAMETER :: S_NBLK=32
      INTEGER, PARAMETER :: S_PERM_IN=31
      INTEGER, PARAMETER :: S_pad4=30
      INTEGER, PARAMETER :: S_A_ELT=29
      INTEGER, PARAMETER :: S_ELTVAR=28
      INTEGER, PARAMETER :: S_ELTPTR=27
      INTEGER, PARAMETER :: S_pad3=26
      INTEGER, PARAMETER :: S_NELT=25
      INTEGER, PARAMETER :: S_pad2=24
      INTEGER, PARAMETER :: S_A_loc=23
      INTEGER, PARAMETER :: S_JCN_loc=22
      INTEGER, PARAMETER :: S_IRN_loc=21
      INTEGER, PARAMETER :: S_NNZ_loc=20
      INTEGER, PARAMETER :: S_pad1=19
      INTEGER, PARAMETER :: S_NZ_loc=18
      INTEGER, PARAMETER :: S_PIVOTS=17
      INTEGER, PARAMETER :: S_COLIND=16
      INTEGER, PARAMETER :: S_ROWIND=15
      INTEGER, PARAMETER :: S_ROWSCA_loc=14
      INTEGER, PARAMETER :: S_COLSCA_loc=13
      INTEGER, PARAMETER :: S_ROWSCA=12
      INTEGER, PARAMETER :: S_COLSCA=11
      INTEGER, PARAMETER :: S_JCN=10
      INTEGER, PARAMETER :: S_IRN=9
      INTEGER, PARAMETER :: S_A=8
      INTEGER, PARAMETER :: S_NNZ=7
      INTEGER, PARAMETER :: S_NZ=6
      INTEGER, PARAMETER :: S_N=5
      INTEGER, PARAMETER :: S_JOB=4
      INTEGER, PARAMETER :: S_PAR=3
      INTEGER, PARAMETER :: S_SYM=2
      INTEGER, PARAMETER :: S_COMM=1
      INTEGER,  PARAMETER :: R_gridinit_done=20
      INTEGER,  PARAMETER :: R_yes=19
      INTEGER,  PARAMETER :: R_RG2L=18
      INTEGER,  PARAMETER :: R_IPIV=17
      INTEGER,  PARAMETER :: R_NB_SINGULAR_VALUES=16
      INTEGER,  PARAMETER :: R_LPIV=15
      INTEGER,  PARAMETER :: R_CNTXT_BLACS=14
      INTEGER,  PARAMETER :: R_DESCRIPTOR=13
      INTEGER,  PARAMETER :: R_TOT_ROOT_SIZE=12
      INTEGER,  PARAMETER :: R_ROOT_SIZE=11
      INTEGER,  PARAMETER :: R_RHS_NLOC=10
      INTEGER,  PARAMETER :: R_SCHUR_LLD=9
      INTEGER,  PARAMETER :: R_SCHUR_NLOC=8
      INTEGER,  PARAMETER :: R_SCHUR_MLOC=7
      INTEGER,  PARAMETER :: R_MYCOL=6
      INTEGER,  PARAMETER :: R_MYROW=5
      INTEGER,  PARAMETER :: R_NPCOL=4
      INTEGER,  PARAMETER :: R_NPROW=3
      INTEGER,  PARAMETER :: R_NBLOCK=2
      INTEGER,  PARAMETER :: R_MBLOCK=1
      OOC_INDICES=(/156,157,159,160/)
      SIZE_INT = id%KEEP(34)
      SIZE_INT8 = id%KEEP(34)*id%KEEP(10)
      SIZE_RL_OR_DBL = id%KEEP(150)
      SIZE_ARITH_DEP = id%KEEP(149)
      SIZE_DOUBLE_PRECISION = 8
      SIZE_LOGICAL = 4
      SIZE_CHARACTER = 1
      size_written=int(0,kind=8)
      tot_NbRecords=0
      NbRecords(:)=0
      NbRecords_ROOTC(:)=0
      size_read=int(0,kind=8)
      size_allocated=int(0,kind=8)
      DIFF_SIZE_ALLOC_READ(:)=0
      DIFF_SIZE_ALLOC_READ_ROOTC(:)=0
      WRITTEN_STRUC_SIZE=int(0,kind=8)
      TMP_OOC_NAMES(:)="?"
      SIZE_VARIABLES_BLR=0_8
      SIZE_GEST_BLR=0
      SIZE_VARIABLES_FRONT_DATA=0_8
      SIZE_GEST_FRONT_DATA=0
      SIZE_VARIABLES_L0FAC=0
      SIZE_GEST_L0FAC=0
#if defined(SINGLEINDOUBLE)
#endif
      if(mode.EQ.memory_save_mode) then
      elseif(mode.EQ.save_mode) then
         write(unit,iostat=err) "MUMPS"
         if(err.ne.0) THEN
            id%INFO(1) = -72
            CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &           id%INFO(2))
         endif 
         CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0 ) GOTO 100
         size_written=size_written+int(5*SIZE_CHARACTER,kind=8)
#if defined(MUMPS_NOF2003)
         size_written=size_written
     &        +int(2*id%KEEP(34)*1,kind=8)
#endif
         IF(id%MYID.EQ.0) THEN
            call date_and_time(date,time,zone,values)
            hash=trim(date)//trim(time)//trim(zone)
         ENDIF
         CALL MPI_BCAST( hash, 23, MPI_CHARACTER, 0, id%COMM, ierr )
         write(unit,iostat=err) hash
         if(err.ne.0) THEN
            id%INFO(1) = -72
            CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &           id%INFO(2))
         endif 
         CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0 ) GOTO 100
         size_written=size_written+int(23*SIZE_CHARACTER,kind=8)
#if defined(MUMPS_NOF2003)
         size_written=size_written
     &        +int(2*id%KEEP(34)*1,kind=8)
#endif
         write(unit,iostat=err) TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
         if(err.ne.0) THEN
            id%INFO(1) = -72
            CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &           id%INFO(2))
         endif 
         CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0 ) GOTO 100
         size_written=size_written+int(2*SIZE_INT8,kind=8)
#if defined(MUMPS_NOF2003)
         size_written=size_written
     &        +int(2*id%KEEP(34)*1,kind=8)
#endif
         ARITH="CMUMPS"(1:1)
         write(unit,iostat=err) ARITH
         if(err.ne.0) THEN
            id%INFO(1) = -72
            CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &           id%INFO(2))
         endif 
         CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0 ) GOTO 100
         size_written=size_written+int(1,kind=8)
#if defined(MUMPS_NOF2003)
         size_written=size_written
     &        +int(2*id%KEEP(34)*1,kind=8)
#endif
         write(unit,iostat=err) id%SYM,id%PAR,id%NPROCS
         if(err.ne.0) THEN
            id%INFO(1) = -72
            CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &           id%INFO(2))
         endif 
         CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0 ) GOTO 100
         size_written=size_written+int(3*SIZE_INT,kind=8)
#if defined(MUMPS_NOF2003)
         size_written=size_written
     &        +int(2*id%KEEP(34)*1,kind=8)
#endif
         IF(id%KEEP(10).EQ.1) THEN
            INT_TYPE_64=.TRUE.
         ELSE
            INT_TYPE_64=.FALSE.
         ENDIF
         write(unit,iostat=err) INT_TYPE_64
         if(err.ne.0) THEN
            id%INFO(1) = -72
            CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &           id%INFO(2))
         endif 
         CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0 ) GOTO 100
         size_written=size_written+int(SIZE_LOGICAL,kind=8)
#if defined(MUMPS_NOF2003)
         size_written=size_written
     &        +int(2*id%KEEP(34)*1,kind=8)
#endif
         IF(associated(id%OOC_FILE_NAME_LENGTH).AND.
     &        associated(id%OOC_FILE_NAMES)) THEN
            write(unit,iostat=err) id%OOC_FILE_NAME_LENGTH(1)
            if(err.ne.0) THEN
               id%INFO(1) = -72
               CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &              id%INFO(2))
            endif 
            CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &           id%COMM, id%MYID )
            IF ( id%INFO(1) .LT. 0 ) GOTO 100
            size_written=size_written+int(SIZE_INT,kind=8)
#if defined(MUMPS_NOF2003)
            size_written=size_written
     &           +int(2*id%KEEP(34)*1,kind=8)
#endif
            TMP_OOC_NAMES(1:id%OOC_FILE_NAME_LENGTH(1))=
     &           id%OOC_FILE_NAMES(1,1:id%OOC_FILE_NAME_LENGTH(1))
            write(unit,iostat=err)
     &           TMP_OOC_NAMES(1:id%OOC_FILE_NAME_LENGTH(1))
            if(err.ne.0) THEN
               id%INFO(1) = -72
               CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &              id%INFO(2))
            endif 
            CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &           id%COMM, id%MYID )
            IF ( id%INFO(1) .LT. 0 ) GOTO 100
            size_written=size_written
     &           +int(id%OOC_FILE_NAME_LENGTH(1)*SIZE_CHARACTER,kind=8)
#if defined(MUMPS_NOF2003)
            size_written=size_written
     &           +int(2*id%KEEP(34)*1,kind=8)
#endif
         ELSE
            write(unit,iostat=err) -999
            if(err.ne.0) THEN
               id%INFO(1) = -72
               CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &              id%INFO(2))
            endif 
            CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &           id%COMM, id%MYID )
            IF ( id%INFO(1) .LT. 0 ) GOTO 100
            size_written=size_written+int(SIZE_INT,kind=8)
#if defined(MUMPS_NOF2003)
            size_written=size_written
     &           +int(2*id%KEEP(34)*1,kind=8)
#endif
            write(unit,iostat=err) -999
            if(err.ne.0) THEN
               id%INFO(1) = -72
               CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &              id%INFO(2))
            endif 
            CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &           id%COMM, id%MYID )
            IF ( id%INFO(1) .LT. 0 ) GOTO 100
            size_written=size_written+int(SIZE_INT,kind=8)
#if defined(MUMPS_NOF2003)
            size_written=size_written
     &           +int(2*id%KEEP(34)*1,kind=8)
#endif
         ENDIF
      elseif((mode.EQ.restore_mode).OR.
     &        (mode.EQ.restore_ooc_mode)) then
        CALL MUMPS_READ_HEADER(unit,err,size_read,SIZE_INT,SIZE_INT8,
     &       TOTAL_FILE_SIZE, TOTAL_STRUC_SIZE, READ_ARITH,
     &       READ_INT_TYPE_64, READ_OOC_FILE_NAME_LENGTH,
     &        READ_OOC_FIRST_FILE_NAME,READ_HASH,
     &        READ_SYM,READ_PAR,READ_NPROCS,FORTRAN_VERSION_OK)
         if (err.ne.0) THEN
            id%INFO(1) = -75
            CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
         elseif (.NOT.FORTRAN_VERSION_OK) THEN
            id%INFO(1) = -73
            id%INFO(2) = 1
         endif
         CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0 ) GOTO 100
         BASIC_CHECK = .false.
         IF (mode.EQ.restore_ooc_mode) THEN
           BASIC_CHECK = .true.
         ENDIF
         CALL CMUMPS_CHECK_HEADER(id,BASIC_CHECK,READ_INT_TYPE_64,
     &              READ_HASH, READ_NPROCS,
     &              READ_ARITH, READ_SYM, READ_PAR)
         IF (id%INFO(1) .LT. 0) GOTO 100
      elseif(mode.EQ.fake_restore_mode) then
         read(unit,iostat=err) READ_HASH
         if(err.ne.0) GOTO 100
         read(unit,iostat=err) TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
         if(err.ne.0) GOTO 100
         IF ( id%INFO(1) .LT. 0 ) GOTO 100
         GOTO 200
      else
         CALL MUMPS_ABORT()
      endif
      DO j=1,size(OOC_INDICES)
         i1=OOC_INDICES(j)
         SELECT CASE(i1)
         CASE(S_OOC_NB_FILES)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%OOC_NB_FILES)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%OOC_NB_FILES,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%OOC_NB_FILES)) THEN
                  write(unit,iostat=err) size(id%OOC_NB_FILES,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%OOC_NB_FILES
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif((mode.EQ.restore_mode).OR.
     &             (mode.EQ.restore_ooc_mode)) then
               nullify(id%OOC_NB_FILES)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%OOC_NB_FILES(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%OOC_NB_FILES
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_OOC_NB_FILE_TYPE)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) id%OOC_NB_FILE_TYPE
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif((mode.EQ.restore_mode).OR.
     &             (mode.EQ.restore_ooc_mode)) then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%OOC_NB_FILE_TYPE
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_OOC_FILE_NAMES)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%OOC_FILE_NAMES)) THEN
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size(id%OOC_FILE_NAMES,1)
     &                 *size(id%OOC_FILE_NAMES,2)*SIZE_CHARACTER
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%OOC_FILE_NAMES)) THEN
                  write(unit,iostat=err) size(id%OOC_FILE_NAMES,1)
     &                 ,size(id%OOC_FILE_NAMES,2)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%OOC_FILE_NAMES
               ELSE
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif((mode.EQ.restore_mode).OR.
     &             (mode.EQ.restore_ooc_mode)) then
               nullify(id%OOC_FILE_NAMES)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size_array1*size_array2
     &                 *SIZE_CHARACTER
                  allocate(id%OOC_FILE_NAMES(size_array1,size_array2),
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%OOC_FILE_NAMES
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_OOC_FILE_NAME_LENGTH)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%OOC_FILE_NAME_LENGTH)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=
     &                 size(id%OOC_FILE_NAME_LENGTH,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%OOC_FILE_NAME_LENGTH)) THEN
                  write(unit,iostat=err) size(id%OOC_FILE_NAME_LENGTH,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%OOC_FILE_NAME_LENGTH
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif((mode.EQ.restore_mode).OR.
     &             (mode.EQ.restore_ooc_mode)) then
               nullify(id%OOC_FILE_NAME_LENGTH)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%OOC_FILE_NAME_LENGTH(size_array1), 
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%OOC_FILE_NAME_LENGTH
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE DEFAULT
         END SELECT
      ENDDO
      if(mode.EQ.restore_ooc_mode) then
         goto 200
      endif
      DO i1=1,NBVARIABLES 
         SELECT CASE(i1)
         CASE(S_COMM)
            if(mode.EQ.memory_save_mode) then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            elseif(mode.EQ.save_mode) then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            elseif(mode.EQ.restore_mode) then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            endif
         CASE(S_SYM)
           CALL MUMPS_SAVE_INT(id%SYM)
           IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_PAR)
           CALL MUMPS_SAVE_INT(id%PAR)
           IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_JOB)
            if(mode.EQ.memory_save_mode) then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            elseif(mode.EQ.save_mode) then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            elseif(mode.EQ.restore_mode) then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            endif
         CASE(S_N)
           CALL MUMPS_SAVE_INT(id%N)
           IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_ICNTL)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT*size(id%ICNTL,1)
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) id%ICNTL
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT*size(id%ICNTL,1)
               read(unit,iostat=err) id%ICNTL
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_INFO)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT*size(id%INFO,1)
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) INFO1,INFO2
     &              ,id%INFO(3:size(id%INFO,1))
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT*size(id%INFO,1)
               read(unit,iostat=err) INFO1,INFO2
     &              ,id%INFO(3:size(id%INFO,1))
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_INFOG)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT*size(id%INFOG,1)
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) INFOG1,INFOG2
     &              ,id%INFOG(3:size(id%INFOG,1))
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT*size(id%INFOG,1)
               read(unit,iostat=err) INFOG1,INFOG2
     &              ,id%INFOG(3:size(id%INFOG,1))
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_COST_SUBTREES)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=SIZE_RL_OR_DBL
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) id%COST_SUBTREES
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=SIZE_RL_OR_DBL
               read(unit,iostat=err) id%COST_SUBTREES
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_CNTL)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=SIZE_RL_OR_DBL*size(id%CNTL,1)
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) id%CNTL
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=SIZE_RL_OR_DBL*size(id%CNTL,1)
               read(unit,iostat=err) id%CNTL
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_RINFO)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=SIZE_RL_OR_DBL*size(id%RINFO,1)
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) id%RINFO
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=SIZE_RL_OR_DBL*size(id%RINFO,1)
               read(unit,iostat=err) id%RINFO
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_RINFOG)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=SIZE_RL_OR_DBL*size(id%RINFOG,1)
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) id%RINFOG
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=SIZE_RL_OR_DBL*size(id%RINFOG,1)
               read(unit,iostat=err) id%RINFOG
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_KEEP8)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT8*size(id%KEEP8,1)
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) id%KEEP8
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT8*size(id%KEEP8,1)
               KEEP883_SAVE=id%KEEP8(83)
               KEEP884_SAVE=id%KEEP8(84)
               read(unit,iostat=err) id%KEEP8
               id%KEEP8(83)=KEEP883_SAVE
               id%KEEP8(84)=KEEP884_SAVE
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_KEEP)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT*size(id%KEEP,1)
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) id%KEEP
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT*size(id%KEEP,1)
               KEEP410_SAVE = id%KEEP(410)
               KEEP411_SAVE = id%KEEP(411)
               read(unit,iostat=err) id%KEEP
               id%KEEP(410) = KEEP410_SAVE
               id%KEEP(411) = KEEP411_SAVE
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_DKEEP)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=SIZE_RL_OR_DBL*size(id%DKEEP,1)
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) id%DKEEP
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=SIZE_RL_OR_DBL*size(id%DKEEP,1)
               read(unit,iostat=err) id%DKEEP
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_NZ)
           CALL MUMPS_SAVE_INT(id%NZ)
           IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_NNZ)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT8
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) id%NNZ
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT8
               read(unit,iostat=err) id%NNZ
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif  
         CASE(S_A)
         CASE(S_IRN)
         CASE(S_JCN)
         CASE(S_COLSCA)
            IF(id%KEEP(52).NE.-1) THEN
              CALL MUMPS_SAVERSTR_REALARRAY(id%COLSCA)
            ELSE
            ENDIF
            IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_ROWSCA)
            IF(id%KEEP(52).NE.-1) THEN
              CALL MUMPS_SAVERSTR_REALARRAY(id%ROWSCA)
            ELSE
            ENDIF
            IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_COLSCA_loc)
            CALL MUMPS_SAVERSTR_REALARRAY(id%COLSCA_loc)
            IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_ROWSCA_loc)
            IS_SYMMETRIC = .FALSE.
            IF (mode.EQ.memory_save_mode .OR.
     &          mode.EQ.save_mode) THEN
              IS_SYMMETRIC = id%KEEP(50).EQ.1 .OR.
     &                       id%KEEP(50).EQ.2
            ELSEIF (mode.EQ.restore_mode) THEN
              IS_SYMMETRIC = READ_SYM.EQ.1 .OR.
     &                       READ_SYM.EQ.2
            ENDIF
            IF ( IS_SYMMETRIC ) THEN
              IF ( mode.EQ.restore_mode ) THEN
                id%ROWSCA_loc => id%COLSCA_loc
              ENDIF
            ELSE
              CALL MUMPS_SAVERSTR_REALARRAY(id%ROWSCA_loc)
            ENDIF
            IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_NZ_loc)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) id%NZ_loc
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%NZ_loc
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_NNZ_loc)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT8
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) id%NNZ_loc
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT8
               read(unit,iostat=err) id%NNZ_loc
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_IRN_loc)
         CASE(S_JCN_loc)
         CASE(S_A_loc)
         CASE(S_NELT)
           CALL MUMPS_SAVE_INT(id%NELT)
           IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_NBLK)
           CALL MUMPS_SAVE_INT(id%NBLK)
           IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_ELTPTR)
         CASE(S_ELTVAR)
         CASE(S_A_ELT)
         CASE(S_PERM_IN)
         CASE(S_BLKPTR)
         CASE(S_BLKVAR)
         CASE(S_COLIND)
         CASE(S_PIVOTS)
         CASE(S_RHS)
         CASE(S_REDRHS)
         CASE(S_ROWIND)
         CASE(S_RHS_SPARSE)
         CASE(S_SOL_loc)
         CASE(S_RHS_loc)
         CASE(S_IRHS_SPARSE)
         CASE(S_IRHS_PTR)
         CASE(S_ISOL_loc)
         CASE(S_IRHS_loc)
         CASE(S_LRHS)
         CALL MUMPS_SAVE_INT(id%LRHS)
         IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_NRHS)
         CALL MUMPS_SAVE_INT(id%NRHS)
         IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_NZ_RHS)
         CALL MUMPS_SAVE_INT(id%NZ_RHS)
         IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_LRHS_loc)
         CALL MUMPS_SAVE_INT(id%LRHS_loc)
         IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_Nloc_RHS)
         CALL MUMPS_SAVE_INT(id%Nloc_RHS)
         IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_LD_RHSINTR)
         CALL MUMPS_SAVE_INT(id%LD_RHSINTR)
         IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_NSOL_loc)
         CALL MUMPS_SAVE_INT(id%NSOL_loc)
         IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_LSOL_loc)
         CALL MUMPS_SAVE_INT(id%LSOL_loc)
         IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_LREDRHS)
         CALL MUMPS_SAVE_INT(id%LREDRHS)
         IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_SYM_PERM)
            CALL CMUMPS_SAVE_INT_SHPTR_ARRAY(id%SYM_PERM
     &           ) 
            IF (id%INFO(1) .LT. 0 ) GOTO 100
         CASE(S_UNS_PERM)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%UNS_PERM)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%UNS_PERM,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%UNS_PERM)) THEN
                  write(unit,iostat=err) size(id%UNS_PERM,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%UNS_PERM
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%UNS_PERM)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%UNS_PERM(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%UNS_PERM
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_NPROW)
         CALL MUMPS_SAVE_INT(id%NPROW)
         IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_NPCOL)
         CALL MUMPS_SAVE_INT(id%NPCOL)
         IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_MBLOCK)
         CALL MUMPS_SAVE_INT(id%MBLOCK)
         IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_NBLOCK)
         CALL MUMPS_SAVE_INT(id%NBLOCK)
         IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_SCHUR_MLOC)
         CALL MUMPS_SAVE_INT(id%SCHUR_MLOC)
         IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_SCHUR_NLOC)
         CALL MUMPS_SAVE_INT(id%SCHUR_NLOC)
         IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_SCHUR_LLD)
         CALL MUMPS_SAVE_INT(id%SCHUR_LLD)
         IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_SIZE_SCHUR)
         CALL MUMPS_SAVE_INT(id%SIZE_SCHUR)
         IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_SCHUR)
         CASE(S_SCHUR_CINTERFACE)
          CASE(S_LISTVAR_SCHUR)
         CASE(S_MAPPING)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%MAPPING)) THEN
                  SIZE_GEST(i1)=SIZE_INT8
                  SIZE_VARIABLES(i1)=id%KEEP8(28)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT8+SIZE_INT
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%MAPPING)) THEN
                  write(unit,iostat=err) id%KEEP8(28)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%MAPPING
               ELSE
                  write(unit,iostat=err) int(-999,8)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%MAPPING)
               read(unit,iostat=err) size_array_INT8_1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array_INT8_1.EQ.int(-999,8)) then
                  SIZE_GEST(i1)=SIZE_INT+SIZE_INT8
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT8
                  SIZE_VARIABLES(i1)=size_array_INT8_1*SIZE_INT
                  allocate(id%MAPPING(size_array_INT8_1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%MAPPING
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_VERSION_NUMBER)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=len(id%VERSION_NUMBER)
     &              *SIZE_CHARACTER
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) id%VERSION_NUMBER
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=len(id%VERSION_NUMBER)
     &              *SIZE_CHARACTER
               read(unit,iostat=err) id%VERSION_NUMBER
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_OOC_TMPDIR)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=len(id%OOC_TMPDIR)
     &              *SIZE_CHARACTER
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) id%OOC_TMPDIR
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=len(id%OOC_TMPDIR)
     &              *SIZE_CHARACTER
               read(unit,iostat=err) id%OOC_TMPDIR
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_OOC_PREFIX)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=len(id%OOC_PREFIX)
     &              *SIZE_CHARACTER
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) id%OOC_PREFIX
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               NbRecords(i1)=1
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=len(id%OOC_PREFIX)
     &              *SIZE_CHARACTER
               read(unit,iostat=err) id%OOC_PREFIX
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_WRITE_PROBLEM)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=len(id%WRITE_PROBLEM)
     &              *SIZE_CHARACTER
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) id%WRITE_PROBLEM
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=len(id%WRITE_PROBLEM)
     &              *SIZE_CHARACTER
               read(unit,iostat=err) id%WRITE_PROBLEM
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_MAX_SURF_MASTER)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT8
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) id%MAX_SURF_MASTER
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT8
               read(unit,iostat=err) id%MAX_SURF_MASTER
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_INST_Number)
           CALL MUMPS_SAVE_INT(id%INST_Number)
           IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_COMM_NODES)
            if(mode.EQ.memory_save_mode) then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            elseif(mode.EQ.save_mode) then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            elseif(mode.EQ.restore_mode) then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            endif
         CASE(S_MYID_NODES)
           CALL MUMPS_SAVE_INT(id%MYID_NODES)
           IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_COMM_LOAD)
            if(mode.EQ.memory_save_mode) then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            elseif(mode.EQ.save_mode) then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            elseif(mode.EQ.restore_mode) then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            endif
         CASE(S_MYID)
           CALL MUMPS_SAVE_INT(id%MYID)
           IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_NPROCS)
           CALL MUMPS_SAVE_INT(id%NPROCS)
           IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_NSLAVES)
           CALL MUMPS_SAVE_INT(id%NSLAVES)
           IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_ASS_IRECV)
           CALL MUMPS_SAVE_INT(id%ASS_IRECV)
           IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_IS)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%IS)) THEN
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=id%KEEP(32)*SIZE_INT
                  DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT*
     &                 (size(id%IS,1)-id%KEEP(32))
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%IS)) THEN
                  write(unit,iostat=err) size(id%IS,1),id%KEEP(32)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%IS(1:id%KEEP(32))
                  DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT*
     &                 (size(id%IS,1)-id%KEEP(32))
               ELSE
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%IS)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size_array2*SIZE_INT
                  DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT*
     &                 (size_array1-size_array2)
                  allocate(id%IS(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%IS(1:size_array2)
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_LNA)
           CALL MUMPS_SAVE_INT(id%LNA)
           IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_NBSA)
           CALL MUMPS_SAVE_INT(id%NBSA)
           IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_STEP)
            CALL CMUMPS_SAVE_INT_SHPTR_ARRAY(id%STEP
     &           ) 
            IF (id%INFO(1) .LT. 0 ) GOTO 100
         CASE(S_NE_STEPS)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%NE_STEPS)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%NE_STEPS,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%NE_STEPS)) THEN
                  write(unit,iostat=err) size(id%NE_STEPS,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%NE_STEPS
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%NE_STEPS)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%NE_STEPS(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%NE_STEPS
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_ND_STEPS)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%ND_STEPS)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%ND_STEPS,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%ND_STEPS)) THEN
                  write(unit,iostat=err) size(id%ND_STEPS,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%ND_STEPS
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%ND_STEPS)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%ND_STEPS(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%ND_STEPS
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_Step2node)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%Step2node)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%Step2node,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%Step2node)) THEN
                  write(unit,iostat=err) size(id%Step2node,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%Step2node
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%Step2node)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%Step2node(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%Step2node
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_FRERE_STEPS)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%FRERE_STEPS)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%FRERE_STEPS,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%FRERE_STEPS)) THEN
                  write(unit,iostat=err) size(id%FRERE_STEPS,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%FRERE_STEPS
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%FRERE_STEPS)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%FRERE_STEPS(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%FRERE_STEPS
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_DAD_STEPS)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%DAD_STEPS)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%DAD_STEPS,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%DAD_STEPS)) THEN
                  write(unit,iostat=err) size(id%DAD_STEPS,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%DAD_STEPS
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%DAD_STEPS)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%DAD_STEPS(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%DAD_STEPS
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_FILS)
            CALL CMUMPS_SAVE_INT_SHPTR_ARRAY(id%FILS
     &           ) 
            IF ( id%INFO(1) .LT. 0 ) GOTO 100
         CASE(S_PTR8ARR)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%PTR8ARR)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%PTR8ARR,1)*SIZE_INT8
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%PTR8ARR)) THEN
                  write(unit,iostat=err) size(id%PTR8ARR,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%PTR8ARR
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               NbRecords(i1)=2
            elseif(mode.EQ.restore_mode) then
               nullify(id%PTR8ARR)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT8
                  allocate(id%PTR8ARR(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%PTR8ARR
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_NINCOLARR)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               NbRecords(i1)=2
               IF(associated(id%NINCOLARR)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%NINCOLARR,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%NINCOLARR)) THEN
                  write(unit,iostat=err) size(id%NINCOLARR,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%NINCOLARR
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%NINCOLARR)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%NINCOLARR(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%NINCOLARR
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_NINROWARR)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               NbRecords(i1)=2
               IF(associated(id%NINROWARR)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%NINROWARR,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%NINROWARR)) THEN
                  write(unit,iostat=err) size(id%NINROWARR,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%NINROWARR
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%NINROWARR)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%NINROWARR(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%NINROWARR
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_PTRDEBARR)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               NbRecords(i1)=2
               IF(associated(id%PTRDEBARR)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%PTRDEBARR,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%PTRDEBARR)) THEN
                  write(unit,iostat=err) size(id%PTRDEBARR,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%PTRDEBARR
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%PTRDEBARR)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%PTRDEBARR(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%PTRDEBARR
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_PTRAR)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%PTRAR)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%PTRAR,1)*SIZE_INT8
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%PTRAR)) THEN
                  write(unit,iostat=err) size(id%PTRAR,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%PTRAR
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               NbRecords(i1)=2
            elseif(mode.EQ.restore_mode) then
               nullify(id%PTRAR)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT8
                  allocate(id%PTRAR(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%PTRAR
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_FRTPTR)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%FRTPTR)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%FRTPTR,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%FRTPTR)) THEN
                  write(unit,iostat=err) size(id%FRTPTR,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%FRTPTR
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%FRTPTR)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%FRTPTR(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%FRTPTR
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_FRTELT)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%FRTELT)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%FRTELT,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%FRTELT)) THEN
                  write(unit,iostat=err) size(id%FRTELT,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  write(unit,iostat=err) id%FRTELT
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%FRTELT)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%FRTELT(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%FRTELT
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_NA)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%NA)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%NA,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%NA)) THEN
                  write(unit,iostat=err) size(id%NA,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  write(unit,iostat=err) id%NA
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%NA)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%NA(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%NA
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_PROCNODE_STEPS)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               NbRecords(i1)=2
               IF(associated(id%PROCNODE_STEPS)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%PROCNODE_STEPS,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%PROCNODE_STEPS)) THEN
                  write(unit,iostat=err) size(id%PROCNODE_STEPS,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%PROCNODE_STEPS
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%PROCNODE_STEPS)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%PROCNODE_STEPS(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%PROCNODE_STEPS
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_PTLUST_S)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%PTLUST_S)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%PTLUST_S,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%PTLUST_S)) THEN
                  write(unit,iostat=err) size(id%PTLUST_S,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%PTLUST_S
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%PTLUST_S)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%PTLUST_S(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%PTLUST_S
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_PTRFAC)
            NbRecords(i1)=2  
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%PTRFAC)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%PTRFAC,1)*SIZE_INT8
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%PTRFAC)) THEN
                  write(unit,iostat=err) size(id%PTRFAC,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  write(unit,iostat=err) id%PTRFAC
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%PTRFAC)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT8
                  allocate(id%PTRFAC(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%PTRFAC
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_S)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%S)) THEN
                  SIZE_GEST(i1)=SIZE_INT8*2
                  SIZE_VARIABLES(i1)=id%KEEP8(31)*SIZE_ARITH_DEP
                  DIFF_SIZE_ALLOC_READ(i1)=
     &                 SIZE_ARITH_DEP*(id%KEEP8(23)-id%KEEP8(31))
               ELSE
                  SIZE_GEST(i1)=SIZE_INT8*2+SIZE_INT
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%S)) THEN
                  write(unit,iostat=err) id%KEEP8(23),id%KEEP8(31)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%S(1:id%KEEP8(31))
                  DIFF_SIZE_ALLOC_READ(i1)=
     &                 SIZE_ARITH_DEP*(id%KEEP8(23)-id%KEEP8(31))
               ELSE
                  write(unit,iostat=err) int(-999,kind=8)
     &                 ,int(-998,kind=8)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%S)
               read(unit,iostat=err) size_array_INT8_1,size_array_INT8_2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array_INT8_1.EQ.int(-999,kind=8)) then
                  SIZE_GEST(i1)=SIZE_INT8*2+SIZE_INT
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT8*2
                  SIZE_VARIABLES(i1)=size_array_INT8_2*SIZE_ARITH_DEP
                  DIFF_SIZE_ALLOC_READ(i1)=
     &                 SIZE_ARITH_DEP*
     &                 (size_array_INT8_1-size_array_INT8_2)
                  allocate(id%S(1:size_array_INT8_1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%S(1:size_array_INT8_2)
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_LPS)    
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%LPS)) THEN
                  SIZE_GEST(i1)=SIZE_INT8*2
                  SIZE_VARIABLES(i1)=id%KEEP8(31)*SIZE_ARITH_DEP/2
                  DIFF_SIZE_ALLOC_READ(i1)=
     &                 SIZE_ARITH_DEP/2*(id%KEEP8(23)-id%KEEP8(31))
               ELSE
                  SIZE_GEST(i1)=SIZE_INT8*2+SIZE_INT
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%LPS)) THEN
                  write(unit,iostat=err) id%KEEP8(23),id%KEEP8(31)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%LPS(1:id%KEEP8(31))
                  DIFF_SIZE_ALLOC_READ(i1)=
     &                 SIZE_ARITH_DEP/2*(id%KEEP8(23)-id%KEEP8(31))
               ELSE
                  write(unit,iostat=err) int(-999,kind=8)
     &                 ,int(-998,kind=8)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%LPS)
               read(unit,iostat=err) size_array_INT8_1,size_array_INT8_2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array_INT8_1.EQ.int(-999,kind=8)) then
                  SIZE_GEST(i1)=SIZE_INT8*2+SIZE_INT
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT8*2
                  SIZE_VARIABLES(i1)=size_array_INT8_2*SIZE_ARITH_DEP/2
                  DIFF_SIZE_ALLOC_READ(i1)=
     &                 SIZE_ARITH_DEP/2*
     &                 (size_array_INT8_1-size_array_INT8_2)
                  allocate(id%LPS(1:size_array_INT8_1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%LPS(1:size_array_INT8_2)
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_PROCNODE)
         CASE(S_NELT_loc)
           CALL MUMPS_SAVE_INT(id%NELT_loc)
           IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_LELTVAR)
           CALL MUMPS_SAVE_INT(id%LELTVAR)
           IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_ELTPROC)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%ELTPROC)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%ELTPROC,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%ELTPROC)) THEN
                  write(unit,iostat=err) size(id%ELTPROC,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%ELTPROC
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%ELTPROC)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%ELTPROC(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%ELTPROC
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_I4_L0_OMP)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%I4_L0_OMP)) THEN
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size(id%I4_L0_OMP,1)
     &                 *size(id%I4_L0_OMP,2)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%I4_L0_OMP)) THEN
                  write(unit,iostat=err) size(id%I4_L0_OMP,1)
     &                 ,size(id%I4_L0_OMP,2)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%I4_L0_OMP
               ELSE
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%I4_L0_OMP)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size_array1*size_array2*SIZE_INT
                  allocate(id%I4_L0_OMP(size_array1,size_array2)
     &                 , stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%I4_L0_OMP
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_I8_L0_OMP)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%I8_L0_OMP)) THEN
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size(id%I8_L0_OMP,1)
     &                 *size(id%I8_L0_OMP,2)*SIZE_INT8
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%I8_L0_OMP)) THEN
                  write(unit,iostat=err) size(id%I8_L0_OMP,1)
     &                 ,size(id%I8_L0_OMP,2)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%I8_L0_OMP
               ELSE
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%I8_L0_OMP)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size_array1*size_array2*SIZE_INT8
                  allocate(id%I8_L0_OMP(size_array1,size_array2)
     &                 , stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%I8_L0_OMP
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_CANDIDATES)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%CANDIDATES)) THEN
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size(id%CANDIDATES,1)
     &                 *size(id%CANDIDATES,2)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%CANDIDATES)) THEN
                  write(unit,iostat=err) size(id%CANDIDATES,1)
     &                 ,size(id%CANDIDATES,2)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%CANDIDATES
               ELSE
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%CANDIDATES)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size_array1*size_array2*SIZE_INT
                  allocate(id%CANDIDATES(size_array1,size_array2)
     &                 , stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%CANDIDATES
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_ISTEP_TO_INIV2)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%ISTEP_TO_INIV2)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%ISTEP_TO_INIV2,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%ISTEP_TO_INIV2)) THEN
                  write(unit,iostat=err) size(id%ISTEP_TO_INIV2,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%ISTEP_TO_INIV2
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%ISTEP_TO_INIV2)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%ISTEP_TO_INIV2(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%ISTEP_TO_INIV2
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_FUTURE_NIV2)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%FUTURE_NIV2)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%FUTURE_NIV2,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%FUTURE_NIV2)) THEN
                  write(unit,iostat=err) size(id%FUTURE_NIV2,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%FUTURE_NIV2
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%FUTURE_NIV2)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%FUTURE_NIV2(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%FUTURE_NIV2
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_TAB_POS_IN_PERE)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%TAB_POS_IN_PERE)) THEN
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size(id%TAB_POS_IN_PERE,1)
     &                 *size(id%TAB_POS_IN_PERE,2)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%TAB_POS_IN_PERE)) THEN
                  write(unit,iostat=err) size(id%TAB_POS_IN_PERE,1)
     &                 ,size(id%TAB_POS_IN_PERE,2)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%TAB_POS_IN_PERE
               ELSE
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%TAB_POS_IN_PERE)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size_array1*size_array2*SIZE_INT
                  allocate(id%TAB_POS_IN_PERE(size_array1,size_array2)
     &                 , stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%TAB_POS_IN_PERE
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_I_AM_CAND)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%I_AM_CAND)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%I_AM_CAND,1)*SIZE_LOGICAL
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%I_AM_CAND)) THEN
                  write(unit,iostat=err) size(id%I_AM_CAND,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%I_AM_CAND
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%I_AM_CAND)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_LOGICAL
                  allocate(id%I_AM_CAND(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%I_AM_CAND
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_MEM_DIST)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%MEM_DIST)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%MEM_DIST,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%MEM_DIST)) THEN
                  write(unit,iostat=err) size(id%MEM_DIST,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  write(unit,iostat=err) id%MEM_DIST
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%MEM_DIST)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%MEM_DIST(0:size_array1-1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%MEM_DIST
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_GLOB2LOC_RHS)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%GLOB2LOC_RHS)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=
     &                 size(id%GLOB2LOC_RHS,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%GLOB2LOC_RHS)) THEN
                  write(unit,iostat=err) size(id%GLOB2LOC_RHS,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%GLOB2LOC_RHS
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%GLOB2LOC_RHS)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%GLOB2LOC_RHS(size_array1), 
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%GLOB2LOC_RHS
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_GLOB2LOC_SOL_ALLOC)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=SIZE_LOGICAL
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) id%GLOB2LOC_SOL_ALLOC
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=SIZE_LOGICAL
               read(unit,iostat=err) id%GLOB2LOC_SOL_ALLOC
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_GLOB2LOC_SOL)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%GLOB2LOC_SOL)) THEN
                  IF(id%GLOB2LOC_SOL_ALLOC) THEN
                     SIZE_GEST(i1)=SIZE_INT
                     SIZE_VARIABLES(i1)=
     &                    size(id%GLOB2LOC_SOL,1)*SIZE_INT
                  ELSE
                     SIZE_GEST(i1)=SIZE_INT*2
                     SIZE_VARIABLES(i1)=0_8
                  ENDIF
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%GLOB2LOC_SOL)) THEN
                  IF(id%GLOB2LOC_SOL_ALLOC) THEN
                     write(unit,iostat=err) size(id%GLOB2LOC_SOL,1)
                     if(err.ne.0) then
                        id%INFO(1) = -72
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_FILE_SIZE-size_written,
     &                       id%INFO(2))
                     endif
                     CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
                     IF ( id%INFO(1) .LT. 0 ) GOTO 100
                     write(unit,iostat=err) id%GLOB2LOC_SOL
                  ELSE
                     write(unit,iostat=err) size(id%GLOB2LOC_SOL,1)
                     if(err.ne.0) then
                        id%INFO(1) = -72
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_FILE_SIZE-size_written,
     &                       id%INFO(2))
                     endif
                     CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
                     IF ( id%INFO(1) .LT. 0 ) GOTO 100
                     write(unit,iostat=err) -999
                  ENDIF
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%GLOB2LOC_SOL)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  if(id%GLOB2LOC_SOL_ALLOC) then
                     SIZE_GEST(i1)=SIZE_INT
                     SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                     allocate(id%GLOB2LOC_SOL(size_array1), 
     &                    stat=allocok)
                     if (allocok .GT. 0) THEN
                        id%INFO(1) = -78
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_STRUC_SIZE-size_allocated
     &                       ,id%INFO(2))
                     endif
                     read(unit,iostat=err) id%GLOB2LOC_SOL
                  else
                     SIZE_GEST(i1)=SIZE_INT*2
                     SIZE_VARIABLES(i1)=0_8
                     read(unit,iostat=err) dummy
                     id%GLOB2LOC_SOL=>id%GLOB2LOC_RHS
                  endif
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_RHSINTR)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%RHSINTR)) THEN
                  SIZE_GEST(i1)=SIZE_INT8
                  SIZE_VARIABLES(i1)=id%KEEP8(25)*SIZE_ARITH_DEP
               ELSE
                  SIZE_GEST(i1)=SIZE_INT8+SIZE_INT
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%RHSINTR)) THEN
                  write(unit,iostat=err) id%KEEP8(25)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%RHSINTR
               ELSE
                  write(unit,iostat=err) int(-999,8)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%RHSINTR)
               read(unit,iostat=err) size_array_INT8_1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array_INT8_1.EQ.int(-999,8)) then
                  SIZE_GEST(i1)=SIZE_INT8+SIZE_INT
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT8
                  SIZE_VARIABLES(i1)=size_array_INT8_1*SIZE_ARITH_DEP
                  allocate(id%RHSINTR(size_array_INT8_1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%RHSINTR
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_MEM_SUBTREE)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%MEM_SUBTREE)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=
     &                 size(id%MEM_SUBTREE,1)*SIZE_DOUBLE_PRECISION
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%MEM_SUBTREE)) THEN
                  write(unit,iostat=err) size(id%MEM_SUBTREE,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%MEM_SUBTREE
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%MEM_SUBTREE)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_DOUBLE_PRECISION
                  allocate(id%MEM_SUBTREE(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%MEM_SUBTREE
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_COST_TRAV)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%COST_TRAV)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=
     &                 size(id%COST_TRAV,1)*SIZE_DOUBLE_PRECISION
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%COST_TRAV)) THEN
                  write(unit,iostat=err) size(id%COST_TRAV,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%COST_TRAV
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%COST_TRAV)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_DOUBLE_PRECISION
                  allocate(id%COST_TRAV(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%COST_TRAV
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_MY_ROOT_SBTR)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%MY_ROOT_SBTR)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%MY_ROOT_SBTR,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%MY_ROOT_SBTR)) THEN
                  write(unit,iostat=err) size(id%MY_ROOT_SBTR,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%MY_ROOT_SBTR
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%MY_ROOT_SBTR)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%MY_ROOT_SBTR(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%MY_ROOT_SBTR
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_MY_FIRST_LEAF)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%MY_FIRST_LEAF)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%MY_FIRST_LEAF,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%MY_FIRST_LEAF)) THEN
                  write(unit,iostat=err) size(id%MY_FIRST_LEAF,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%MY_FIRST_LEAF
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%MY_FIRST_LEAF)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%MY_FIRST_LEAF(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%MY_FIRST_LEAF
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_MY_NB_LEAF)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%MY_NB_LEAF)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%MY_NB_LEAF,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%MY_NB_LEAF)) THEN
                  write(unit,iostat=err) size(id%MY_NB_LEAF,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%MY_NB_LEAF
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%MY_NB_LEAF)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%MY_NB_LEAF(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%MY_NB_LEAF
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_DEPTH_FIRST)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%DEPTH_FIRST)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%DEPTH_FIRST,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%DEPTH_FIRST)) THEN
                  write(unit,iostat=err) size(id%DEPTH_FIRST,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%DEPTH_FIRST
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%DEPTH_FIRST)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%DEPTH_FIRST(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%DEPTH_FIRST
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_DEPTH_FIRST_SEQ)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%DEPTH_FIRST_SEQ)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%DEPTH_FIRST_SEQ,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%DEPTH_FIRST_SEQ)) THEN
                  write(unit,iostat=err) size(id%DEPTH_FIRST_SEQ,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%DEPTH_FIRST_SEQ
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%DEPTH_FIRST_SEQ)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%DEPTH_FIRST_SEQ(size_array1), 
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%DEPTH_FIRST_SEQ
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_SBTR_ID)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%SBTR_ID)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%SBTR_ID,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%SBTR_ID)) THEN
                  write(unit,iostat=err) size(id%SBTR_ID,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%SBTR_ID
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%SBTR_ID)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%SBTR_ID(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%SBTR_ID
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_SCHED_DEP)
         CASE(S_SCHED_GRP)
         CASE(S_CROIX_MANU)
         CASE(S_WK_USER)
         CASE(S_NBSA_LOCAL)
           CALL MUMPS_SAVE_INT(id%NBSA_LOCAL)
           IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_LWK_USER)
            if(mode.EQ.memory_save_mode) then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            elseif(mode.EQ.save_mode) then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            elseif(mode.EQ.restore_mode) then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            endif
         CASE(S_CB_SON_SIZE)
         CASE(S_INSTANCE_NUMBER)
           CALL MUMPS_SAVE_INT(id%INSTANCE_NUMBER)
           IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_OOC_MAX_NB_NODES_FOR_ZONE)
           CALL MUMPS_SAVE_INT(id%OOC_MAX_NB_NODES_FOR_ZONE)
           IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_OOC_INODE_SEQUENCE)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%OOC_INODE_SEQUENCE)) THEN
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size(id%OOC_INODE_SEQUENCE,1)
     &                 *size(id%OOC_INODE_SEQUENCE,2)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%OOC_INODE_SEQUENCE)) THEN
                  write(unit,iostat=err) size(id%OOC_INODE_SEQUENCE,1)
     &                 ,size(id%OOC_INODE_SEQUENCE,2)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%OOC_INODE_SEQUENCE
               ELSE
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%OOC_INODE_SEQUENCE)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size_array1*size_array2*SIZE_INT
                  allocate(id%OOC_INODE_SEQUENCE(size_array1
     &                 ,size_array2), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%OOC_INODE_SEQUENCE
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_OOC_SIZE_OF_BLOCK)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%OOC_SIZE_OF_BLOCK)) THEN
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size(id%OOC_SIZE_OF_BLOCK,1)
     &                 *size(id%OOC_SIZE_OF_BLOCK,2)*SIZE_INT8
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%OOC_SIZE_OF_BLOCK)) THEN
                  write(unit,iostat=err) size(id%OOC_SIZE_OF_BLOCK,1)
     &                 ,size(id%OOC_SIZE_OF_BLOCK,2)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%OOC_SIZE_OF_BLOCK
               ELSE
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%OOC_SIZE_OF_BLOCK)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size_array1*size_array2*SIZE_INT8
                  allocate(id%OOC_SIZE_OF_BLOCK(size_array1
     &                 ,size_array2), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%OOC_SIZE_OF_BLOCK
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_OOC_VADDR)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%OOC_VADDR)) THEN
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size(id%OOC_VADDR,1)
     &                 *size(id%OOC_VADDR,2)*SIZE_INT8
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%OOC_VADDR)) THEN
                  write(unit,iostat=err) size(id%OOC_VADDR,1)
     &                 ,size(id%OOC_VADDR,2)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%OOC_VADDR
               ELSE
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%OOC_VADDR)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size_array1*size_array2*SIZE_INT8
                  allocate(id%OOC_VADDR(size_array1,size_array2),
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%OOC_VADDR
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_OOC_TOTAL_NB_NODES)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%OOC_TOTAL_NB_NODES)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=
     &                 size(id%OOC_TOTAL_NB_NODES,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%OOC_TOTAL_NB_NODES)) THEN
                  write(unit,iostat=err) size(id%OOC_TOTAL_NB_NODES,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%OOC_TOTAL_NB_NODES
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%OOC_TOTAL_NB_NODES)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%OOC_TOTAL_NB_NODES(size_array1), 
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%OOC_TOTAL_NB_NODES
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_OOC_NB_FILES)
         CASE(S_OOC_NB_FILE_TYPE)
         CASE(S_OOC_FILE_NAMES)
         CASE(S_OOC_FILE_NAME_LENGTH)
         CASE(S_PIVNUL_LIST)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%PIVNUL_LIST)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%PIVNUL_LIST,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%PIVNUL_LIST)) THEN
                  write(unit,iostat=err) size(id%PIVNUL_LIST,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%PIVNUL_LIST
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%PIVNUL_LIST)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%PIVNUL_LIST(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%PIVNUL_LIST
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_SUP_PROC)
         CASE(S_IPTR_WORKING)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%IPTR_WORKING)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%IPTR_WORKING,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%IPTR_WORKING)) THEN
                  write(unit,iostat=err) size(id%IPTR_WORKING,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%IPTR_WORKING
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%IPTR_WORKING)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%IPTR_WORKING(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%IPTR_WORKING
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_WORKING)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%WORKING)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%WORKING,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%WORKING)) THEN
                  write(unit,iostat=err) size(id%WORKING,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%WORKING
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%WORKING)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%WORKING(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%WORKING
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_INTR_ENCODING)
            NbRecords(i1)     =0   
            SIZE_GEST(i1)     =0   
            SIZE_VARIABLES(i1)=0_8 
#if defined(SINGLEINDOUBLE)
#endif
            DO i2=1,NBVARIABLES_ROOTC
            SELECT CASE(i2)
            CASE(R_MBLOCK)
               NbRecords_ROOTC(i2)=1
               if(mode.EQ.memory_save_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
               elseif(mode.EQ.save_mode) then
                  write(unit,iostat=err) idintr%root%MBLOCK
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(mode.EQ.restore_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
                  read(unit,iostat=err) idintr%root%MBLOCK
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOTC(i2)
     &                    ,id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE(R_NBLOCK)
               NbRecords_ROOTC(i2)=1
               if(mode.EQ.memory_save_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
               elseif(mode.EQ.save_mode) then
                  write(unit,iostat=err) idintr%root%NBLOCK
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100        
               elseif(mode.EQ.restore_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
                  read(unit,iostat=err) idintr%root%NBLOCK
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOTC(i2)
     &                    ,id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
            CASE(R_NPROW)
               NbRecords_ROOTC(i2)=1
               if(mode.EQ.memory_save_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
               elseif(mode.EQ.save_mode) then
                  write(unit,iostat=err) idintr%root%NPROW
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100                 
               elseif(mode.EQ.restore_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
                  read(unit,iostat=err) idintr%root%NPROW
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOTC(i2)
     &                    ,id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE(R_NPCOL)
               NbRecords_ROOTC(i2)=1
               if(mode.EQ.memory_save_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
               elseif(mode.EQ.save_mode) then
                  write(unit,iostat=err) idintr%root%NPCOL
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100                
               elseif(mode.EQ.restore_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
                  read(unit,iostat=err) idintr%root%NPCOL
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOTC(i2)
     &                    ,id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE(R_MYROW)
               NbRecords_ROOTC(i2)=1
               if(mode.EQ.memory_save_mode) then
                  NbRecords_ROOTC(i2)=1
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
               elseif(mode.EQ.save_mode) then
                  write(unit,iostat=err) idintr%root%MYROW
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(mode.EQ.restore_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
                  read(unit,iostat=err) idintr%root%MYROW
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOTC(i2),
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE(R_MYCOL)
               NbRecords_ROOTC(i2)=1
               if(mode.EQ.memory_save_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
               elseif(mode.EQ.save_mode) then
                  write(unit,iostat=err) idintr%root%MYCOL
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(mode.EQ.restore_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
                  read(unit,iostat=err) idintr%root%MYCOL
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOTC(i2)
     &                    ,id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE(R_SCHUR_MLOC)
               NbRecords_ROOTC(i2)=1
               if(mode.EQ.memory_save_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
               elseif(mode.EQ.save_mode) then
                  write(unit,iostat=err) idintr%root%SCHUR_MLOC
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(mode.EQ.restore_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
                  read(unit,iostat=err) idintr%root%SCHUR_MLOC
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOTC(i2)
     &                    ,id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE(R_SCHUR_NLOC)
               NbRecords_ROOTC(i2)=1
               if(mode.EQ.memory_save_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
               elseif(mode.EQ.save_mode) then
                  write(unit,iostat=err) idintr%root%SCHUR_NLOC
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(mode.EQ.restore_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
                  read(unit,iostat=err) idintr%root%SCHUR_NLOC
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOTC(i2)
     &                    ,id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE(R_SCHUR_LLD)
               NbRecords_ROOTC(i2)=1
               if(mode.EQ.memory_save_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
               elseif(mode.EQ.save_mode) then
                  write(unit,iostat=err) idintr%root%SCHUR_LLD
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(mode.EQ.restore_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
                  read(unit,iostat=err) idintr%root%SCHUR_LLD
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOTC(i2),
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE(R_RHS_NLOC)
               NbRecords_ROOTC(i2)=1
               if(mode.EQ.memory_save_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
               elseif(mode.EQ.save_mode) then
                  write(unit,iostat=err) idintr%root%RHS_NLOC
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(mode.EQ.restore_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
                  read(unit,iostat=err) idintr%root%RHS_NLOC
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOTC(i2),
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE(R_ROOT_SIZE)
               NbRecords_ROOTC(i2)=1
               if(mode.EQ.memory_save_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
               elseif(mode.EQ.save_mode) then
                  write(unit,iostat=err) idintr%root%ROOT_SIZE
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(mode.EQ.restore_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
                  read(unit,iostat=err) idintr%root%ROOT_SIZE
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOTC(i2),
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE(R_TOT_ROOT_SIZE)
               NbRecords_ROOTC(i2)=1
               if(mode.EQ.memory_save_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
               elseif(mode.EQ.save_mode) then
                  write(unit,iostat=err) idintr%root%TOT_ROOT_SIZE
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(mode.EQ.restore_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
                  read(unit,iostat=err) idintr%root%TOT_ROOT_SIZE
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOTC(i2),
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE(R_DESCRIPTOR)
               NbRecords_ROOTC(i2)=1
               if(mode.EQ.memory_save_mode) then
                 SIZE_VARIABLES_ROOTC(i2)=
     &                 size(idintr%root%DESCRIPTOR,1) * SIZE_INT
               elseif(mode.EQ.save_mode) then
                  write(unit,iostat=err) idintr%root%DESCRIPTOR
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(mode.EQ.restore_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT*
     &                 size(idintr%root%DESCRIPTOR,1)
                  read(unit,iostat=err) idintr%root%DESCRIPTOR
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOTC(i2),
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE(R_CNTXT_BLACS)
               NbRecords_ROOTC(i2)=1
               if(mode.EQ.memory_save_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
               elseif(mode.EQ.save_mode) then
                  write(unit,iostat=err) idintr%root%CNTXT_BLACS
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(mode.EQ.restore_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
                  read(unit,iostat=err) idintr%root%CNTXT_BLACS
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOTC(i2),
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE(R_LPIV)
               NbRecords_ROOTC(i2)=1
               if(mode.EQ.memory_save_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
               elseif(mode.EQ.save_mode) then
                  write(unit,iostat=err) idintr%root%LPIV
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(mode.EQ.restore_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
                  read(unit,iostat=err) idintr%root%LPIV
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOTC(i2),
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE(R_RG2L)
            CASE(R_IPIV)
               NbRecords_ROOTC(i2)=2
               if(mode.EQ.memory_save_mode) then
                  IF(associated(idintr%root%IPIV)) THEN
                     SIZE_GEST_ROOTC(i2)=SIZE_INT
                     SIZE_VARIABLES_ROOTC(i2)=
     &                    size(idintr%root%IPIV,1)*SIZE_INT
                  ELSE
                     SIZE_GEST_ROOTC(i2)=SIZE_INT*2
                     SIZE_VARIABLES_ROOTC(i2)=0_8
                  ENDIF
               elseif(mode.EQ.save_mode) then
                  IF(associated(idintr%root%IPIV)) THEN
                     write(unit,iostat=err) size(idintr%root%IPIV,1)
                     if(err.ne.0) then
                        id%INFO(1) = -72
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_FILE_SIZE-size_written
     &                       ,id%INFO(2))
                     endif
                     CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
                     IF ( id%INFO(1) .LT. 0 ) GOTO 100
                     write(unit,iostat=err) idintr%root%IPIV
                  ELSE
                     write(unit,iostat=err) -999
                     if(err.ne.0) then
                        id%INFO(1) = -72
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_FILE_SIZE-size_written
     &                       ,id%INFO(2))
                     endif
                     CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
                     IF ( id%INFO(1) .LT. 0 ) GOTO 100
                     write(unit,iostat=err) -999
                  ENDIF 
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(mode.EQ.restore_mode) then
                  nullify(idintr%root%IPIV)
                  read(unit,iostat=err) size_array1
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOTC(i2),
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  if(size_array1.EQ.-999) then
                     SIZE_GEST_ROOTC(i2)=SIZE_INT*2
                     SIZE_VARIABLES_ROOTC(i2)=0_8
                     read(unit,iostat=err) dummy
                  else
                     SIZE_GEST_ROOTC(i2)=SIZE_INT
                     SIZE_VARIABLES_ROOTC(i2)=size_array1*SIZE_INT
                     allocate(idintr%root%IPIV(size_array1),
     &                        stat=allocok)
                     if (allocok .GT. 0) THEN
                        id%INFO(1) = -78
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_STRUC_SIZE-size_allocated
     &                       ,id%INFO(2))
                     endif
                     read(unit,iostat=err) idintr%root%IPIV
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOTC(i2),
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE(R_yes)
               NbRecords_ROOTC(i2)=1
               if(mode.EQ.memory_save_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_LOGICAL
               elseif(mode.EQ.save_mode) then
                  write(unit,iostat=err) idintr%root%yes
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(mode.EQ.restore_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_LOGICAL
                  read(unit,iostat=err) idintr%root%yes
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read,
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE(R_gridinit_done)
               NbRecords_ROOTC(i2)=1
               if(mode.EQ.memory_save_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_LOGICAL
               elseif(mode.EQ.save_mode) then
                  write(unit,iostat=err) idintr%root%gridinit_done
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(mode.EQ.restore_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_LOGICAL
                  read(unit,iostat=err) idintr%root%gridinit_done
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read,
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE(R_NB_SINGULAR_VALUES)
               NbRecords_ROOTC(i2)=1
               if(mode.EQ.memory_save_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
               elseif(mode.EQ.save_mode) then
                  write(unit,iostat=err) idintr%root%NB_SINGULAR_VALUES
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(mode.EQ.restore_mode) then
                  SIZE_VARIABLES_ROOTC(i2)=SIZE_INT
                  read(unit,iostat=err) idintr%root%NB_SINGULAR_VALUES
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                    ,id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE DEFAULT
            END SELECT
            if(mode.EQ.memory_save_mode) then
               NbSubRecords=int(SIZE_VARIABLES_ROOTC(i2)/huge(I4))
               IF(NbSubRecords.GT.0) then
                  NbRecords_ROOTC(i2)=NbRecords_ROOTC(i2)+NbSubRecords
               ENDIF
            elseif(mode.EQ.save_mode) then
               size_written=size_written+SIZE_VARIABLES_ROOTC(i2)
     &              +int(SIZE_GEST_ROOTC(i2),kind=8)
#if defined(MUMPS_NOF2003)
               size_written=size_written
     &              +int(2*id%KEEP(34)*NbRecords_ROOTC(i2),kind=8)
#endif
            elseif(mode.EQ.restore_mode) then
               size_allocated=size_allocated+SIZE_VARIABLES_ROOTC(i2)+
     &              DIFF_SIZE_ALLOC_READ_ROOTC(i2)
               size_read=size_read+SIZE_VARIABLES_ROOTC(i2)
     &              +int(SIZE_GEST_ROOTC(i2),kind=8)
#if defined(MUMPS_NOF2003)
               size_read=size_read
     &              +int(2*id%KEEP(34)*NbRecords_ROOTC(i2),kind=8)
#endif
            elseif(mode.EQ.fake_restore_mode) then
            endif
         ENDDO
            CALL CMUMPS_SAVE_RESTORE_L0FACARRAY(
     &                  idintr%L0_OMP_FACTORS
     &                 ,unit,id%MYID,mode
     &                 ,SIZE_GEST_L0FAC,SIZE_VARIABLES_L0FAC
     &                 ,SIZE_INT, SIZE_INT8, SIZE_ARITH_DEP
     &                 ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                 ,size_read,size_allocated,size_written
     &                 ,id%INFO(1))
            CALL CMUMPS_SAVE_RESTORE_ROOTA(
     &                 idintr%roota
     &                 ,unit,id%MYID,mode
     &                 ,SIZE_GEST_ROOTA,SIZE_VARIABLES_ROOTA
     &                 ,SIZE_INT, SIZE_INT8, SIZE_ARITH_DEP
     &                 ,SIZE_RL_OR_DBL
     &                 ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                 ,size_read,size_allocated,size_written
     &                 ,id%INFO)
            CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
            IF ( id%INFO(1) .LT. 0 ) GOTO 100
#if defined(SINGLEINDOUBLE)
#endif
         CASE(S_NBGRP)
           CALL MUMPS_SAVE_INT(id%NBGRP)
           IF (id%INFO(1) .LT. 0) GOTO 100
         CASE(S_LRGROUPS)
            CALL CMUMPS_SAVE_INT_SHPTR_ARRAY(id%LRGROUPS
     &           ) 
            IF (id%INFO(1) .LT. 0 ) GOTO 100
         CASE(S_FDM_F_ENCODING)
            NbRecords(i1)=1
            SIZE_GEST(i1)=SIZE_INT
            SIZE_VARIABLES(i1)=0_8
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%FDM_F_ENCODING)) THEN
                  CALL MUMPS_SAVE_RESTORE_FRONT_DATA(id%FDM_F_ENCODING
     &                 ,unit,id%MYID,memory_save_mode
     &                 ,SIZE_GEST_FRONT_DATA,SIZE_VARIABLES_FRONT_DATA
     &                 ,SIZE_INT,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                 ,size_read,size_allocated,size_written
     &                 ,id%INFO(1))                 
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%FDM_F_ENCODING)) THEN
                  write(unit,iostat=err) size(id%FDM_F_ENCODING,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  CALL MUMPS_SAVE_RESTORE_FRONT_DATA(id%FDM_F_ENCODING
     &                 ,unit,id%MYID,save_mode
     &                 ,SIZE_GEST_FRONT_DATA,SIZE_VARIABLES_FRONT_DATA
     &                 ,SIZE_INT,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                 ,size_read,size_allocated,size_written
     &                 ,id%INFO(1))
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%FDM_F_ENCODING)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.NE.-999) then
                  CALL MUMPS_SAVE_RESTORE_FRONT_DATA(id%FDM_F_ENCODING
     &                 ,unit,id%MYID,restore_mode
     &                 ,SIZE_GEST_FRONT_DATA,SIZE_VARIABLES_FRONT_DATA
     &                 ,SIZE_INT,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                 ,size_read,size_allocated,size_written
     &                 ,id%INFO(1))   
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_BLRARRAY_ENCODING)            
            NbRecords(i1)=1
            SIZE_GEST(i1)=SIZE_INT
            SIZE_VARIABLES(i1)=0_8
            CALL_SAVE_RESTORE_BLR = .FALSE.
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%BLRARRAY_ENCODING)) THEN
                 CALL_SAVE_RESTORE_BLR = .TRUE.
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%BLRARRAY_ENCODING)) THEN
                  write(unit,iostat=err) size(id%BLRARRAY_ENCODING,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  CALL_SAVE_RESTORE_BLR = .TRUE.
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               ENDIF 
            elseif(mode.EQ.restore_mode) then
               nullify(id%BLRARRAY_ENCODING)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.NE.-999) then
                 CALL_SAVE_RESTORE_BLR = .TRUE.
               endif
            endif
            IF (CALL_SAVE_RESTORE_BLR) THEN
#if defined(SINGLEINDOUBLE)
#endif
                 CALL CMUMPS_SAVE_RESTORE_BLR(id%BLRARRAY_ENCODING
     &                 ,unit,id%MYID,mode
     &                 ,SIZE_GEST_BLR,SIZE_VARIABLES_BLR
     &                 ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &                 ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                 ,size_read,size_allocated,size_written
     &                 ,id%INFO(1))
#if defined(SINGLEINDOUBLE)
#endif
            ENDIF
            CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &            id%COMM, id%MYID )
            IF ( id%INFO(1) .LT. 0 ) GOTO 100
         CASE(S_SCHED_SBTR)
         CASE(S_LPOOL_A_L0_OMP)
           CALL MUMPS_SAVE_INT(id%LPOOL_A_L0_OMP)
           IF ( id%INFO(1) .LT. 0 ) GOTO 100
         CASE(S_LPOOL_B_L0_OMP)
           CALL MUMPS_SAVE_INT(id%LPOOL_B_L0_OMP)
           IF ( id%INFO(1) .LT. 0 ) GOTO 100
         CASE(S_L_PHYS_L0_OMP)
           CALL MUMPS_SAVE_INT(id%L_PHYS_L0_OMP)
           IF ( id%INFO(1) .LT. 0 ) GOTO 100
         CASE(S_L_VIRT_L0_OMP)     
           CALL MUMPS_SAVE_INT(id%L_VIRT_L0_OMP)
           IF ( id%INFO(1) .LT. 0 ) GOTO 100
         CASE(S_LL0_OMP_MAPPING)
           CALL MUMPS_SAVE_INT(id%LL0_OMP_MAPPING)
           IF ( id%INFO(1) .LT. 0 ) GOTO 100
         CASE(S_LL0_OMP_FACTORS)
           CALL MUMPS_SAVE_INT(id%LL0_OMP_FACTORS)
           IF ( id%INFO(1) .LT. 0 ) GOTO 100
         CASE(S_THREAD_LA)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT8
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) id%THREAD_LA
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT8
               read(unit,iostat=err) id%THREAD_LA
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_IPOOL_A_L0_OMP)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%IPOOL_A_L0_OMP)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=
     &                 size(id%IPOOL_A_L0_OMP,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%IPOOL_A_L0_OMP)) THEN
                  write(unit,iostat=err) size(id%IPOOL_A_L0_OMP,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%IPOOL_A_L0_OMP
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%IPOOL_A_L0_OMP)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%IPOOL_A_L0_OMP(size_array1), 
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%IPOOL_A_L0_OMP
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_IPOOL_B_L0_OMP)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%IPOOL_B_L0_OMP)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=
     &                 size(id%IPOOL_B_L0_OMP,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%IPOOL_B_L0_OMP)) THEN
                  write(unit,iostat=err) size(id%IPOOL_B_L0_OMP,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%IPOOL_B_L0_OMP
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%IPOOL_B_L0_OMP)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%IPOOL_B_L0_OMP(size_array1), 
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%IPOOL_B_L0_OMP
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_PHYS_L0_OMP)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%PHYS_L0_OMP)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=
     &                 size(id%PHYS_L0_OMP,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%PHYS_L0_OMP)) THEN
                  write(unit,iostat=err) size(id%PHYS_L0_OMP,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%PHYS_L0_OMP
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%PHYS_L0_OMP)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%PHYS_L0_OMP(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%PHYS_L0_OMP
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_VIRT_L0_OMP)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%VIRT_L0_OMP)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%VIRT_L0_OMP,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%VIRT_L0_OMP)) THEN
                  write(unit,iostat=err) size(id%VIRT_L0_OMP,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%VIRT_L0_OMP
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%VIRT_L0_OMP)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%VIRT_L0_OMP(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%VIRT_L0_OMP
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_VIRT_L0_OMP_MAPPING)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%VIRT_L0_OMP_MAPPING)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%VIRT_L0_OMP_MAPPING,1)
     &                              *SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%VIRT_L0_OMP_MAPPING)) THEN
                  write(unit,iostat=err) size(id%VIRT_L0_OMP_MAPPING,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%VIRT_L0_OMP_MAPPING
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%VIRT_L0_OMP_MAPPING)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%VIRT_L0_OMP_MAPPING(size_array1),
     &                     stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%VIRT_L0_OMP_MAPPING
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_PERM_L0_OMP)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%PERM_L0_OMP)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%PERM_L0_OMP,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%PERM_L0_OMP)) THEN
                  write(unit,iostat=err) size(id%PERM_L0_OMP,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%PERM_L0_OMP
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%PERM_L0_OMP)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%PERM_L0_OMP(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%PERM_L0_OMP
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_PTR_LEAFS_L0_OMP)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%PTR_LEAFS_L0_OMP)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=
     &                 size(id%PTR_LEAFS_L0_OMP,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%PTR_LEAFS_L0_OMP)) THEN
                  write(unit,iostat=err) size(id%PTR_LEAFS_L0_OMP,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%PTR_LEAFS_L0_OMP
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%PTR_LEAFS_L0_OMP)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%PTR_LEAFS_L0_OMP(size_array1), 
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%PTR_LEAFS_L0_OMP
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_L0_OMP_MAPPING)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%L0_OMP_MAPPING)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=
     &                 size(id%L0_OMP_MAPPING,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%L0_OMP_MAPPING)) THEN
                  write(unit,iostat=err) size(id%L0_OMP_MAPPING,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%L0_OMP_MAPPING
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%L0_OMP_MAPPING)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%L0_OMP_MAPPING(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%L0_OMP_MAPPING
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_SINGULAR_VALUES)
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%SINGULAR_VALUES)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=
     &                 size(id%SINGULAR_VALUES,1)*SIZE_RL_OR_DBL
              ELSE
                 SIZE_GEST(i1)=SIZE_INT*2
                 SIZE_VARIABLES(i1)=0_8
              ENDIF
           elseif(mode.EQ.save_mode) then
              IF(associated(id%SINGULAR_VALUES)) THEN
                 write(unit,iostat=err) size(id%SINGULAR_VALUES,1)
                 if(err.ne.0) then
                    id%INFO(1) = -72
                    CALL MUMPS_SETI8TOI4(size_written,id%INFO(2))
                 endif
                 CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                id%COMM, id%MYID )
                 IF ( id%INFO(1) .LT. 0 ) GOTO 100
                 write(unit,iostat=err) id%SINGULAR_VALUES
              ELSE
                 write(unit,iostat=err) -999
                 if(err.ne.0) then
                    id%INFO(1) = -72
                    CALL MUMPS_SETI8TOI4(size_written,id%INFO(2))
                 endif
                 CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                id%COMM, id%MYID )
                 IF ( id%INFO(1) .LT. 0 ) GOTO 100
                 write(unit,iostat=err) -999
              ENDIF 
              if(err.ne.0) then
                 id%INFO(1) = -72
                 CALL MUMPS_SETI8TOI4(size_written,id%INFO(2))
              endif
              CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &             id%COMM, id%MYID )
              IF ( id%INFO(1) .LT. 0 ) GOTO 100
           elseif(mode.EQ.restore_mode) then
              nullify(id%SINGULAR_VALUES)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_RL_OR_DBL
                  allocate(id%SINGULAR_VALUES(size_array1),
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%SINGULAR_VALUES
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
           endif
        CASE(S_NB_SINGULAR_VALUES)
           CALL MUMPS_SAVE_INT(id%NB_SINGULAR_VALUES)
           IF ( id%INFO(1) .LT. 0 ) GOTO 100
        CASE(S_ASSOCIATED_OOC_FILES)
            if(mode.EQ.memory_save_mode) then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_LOGICAL
            elseif(mode.EQ.save_mode) then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_LOGICAL
            elseif(mode.EQ.restore_mode) then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_LOGICAL
            endif
        CASE(S_SAVE_DIR)
           NbRecords(i1)=1
           if(mode.EQ.memory_save_mode) then
              SIZE_VARIABLES(i1)=len(id%SAVE_DIR)*SIZE_CHARACTER
           elseif(mode.EQ.save_mode) then
              write(unit,iostat=err) id%SAVE_DIR
              if(err.ne.0) then
                 id%INFO(1) = -72
                 CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                id%INFO(2))
              endif
              CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &             id%COMM, id%MYID )
              IF ( id%INFO(1) .LT. 0 ) GOTO 100
           elseif(mode.EQ.restore_mode) then
              SIZE_VARIABLES(i1)=len(id%SAVE_DIR)*SIZE_CHARACTER
              read(unit,iostat=err) id%SAVE_DIR
              if(err.ne.0) THEN
                 id%INFO(1) = -75
                 CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                ,id%INFO(2))
              endif 
              CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &             id%COMM, id%MYID )
              IF ( id%INFO(1) .LT. 0 ) GOTO 100
           endif
        CASE(S_SAVE_PREFIX)
           NbRecords(i1)=1
           if(mode.EQ.memory_save_mode) then
              SIZE_VARIABLES(i1)=len(id%SAVE_PREFIX)*SIZE_CHARACTER
           elseif(mode.EQ.save_mode) then
              write(unit,iostat=err) id%SAVE_PREFIX
              if(err.ne.0) then
                 id%INFO(1) = -72
                 CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                id%INFO(2))
              endif
              CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &             id%COMM, id%MYID )
              IF ( id%INFO(1) .LT. 0 ) GOTO 100
           elseif(mode.EQ.restore_mode) then
              SIZE_VARIABLES(i1)=len(id%SAVE_PREFIX)*SIZE_CHARACTER
              read(unit,iostat=err) id%SAVE_PREFIX
              if(err.ne.0) THEN
                 id%INFO(1) = -75
                 CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                ,id%INFO(2))
              endif 
              CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &             id%COMM, id%MYID )
              IF ( id%INFO(1) .LT. 0 ) GOTO 100
           endif
        CASE(S_MTKO_PROCS_MAP)
           NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id%MTKO_PROCS_MAP)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=
     &                 size(id%MTKO_PROCS_MAP,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id%MTKO_PROCS_MAP)) THEN
                  write(unit,iostat=err) size(id%MTKO_PROCS_MAP,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%MTKO_PROCS_MAP
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id%MTKO_PROCS_MAP)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%MTKO_PROCS_MAP(0:size_array1-1),
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%MTKO_PROCS_MAP
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE(S_METIS_OPTIONS)
            NbRecords(i1)=1
            if(mode.EQ.memory_save_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT*size(id%METIS_OPTIONS,1)
            elseif(mode.EQ.save_mode) then
               write(unit,iostat=err) id%METIS_OPTIONS
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               SIZE_VARIABLES(i1)=SIZE_INT*size(id%METIS_OPTIONS,1)
               read(unit,iostat=err) id%METIS_OPTIONS
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif      
      CASE(S_pad1,S_pad2,S_pad3,S_pad4,S_pad5,S_pad6,S_pad7,
     &     S_pad8,S_pad9,S_pad10,S_pad11,S_pad12,S_pad13)
      CASE DEFAULT
      END SELECT
      if(mode.EQ.memory_save_mode) then
         NbSubRecords=int(SIZE_VARIABLES(i1)/huge(I4))
         IF(NbSubRecords.GT.0) then
            NbRecords(i1)=NbRecords(i1)+NbSubRecords
         ENDIF
      elseif(mode.EQ.save_mode) then
         size_written=size_written+SIZE_VARIABLES(i1)
     &        +int(SIZE_GEST(i1),kind=8)
#if defined(MUMPS_NOF2003)
         size_written=size_written
     &        +int(2*id%KEEP(34)*NbRecords(i1),kind=8)
#endif
      elseif(mode.EQ.restore_mode) then
         size_allocated=size_allocated+SIZE_VARIABLES(i1)+
     &        DIFF_SIZE_ALLOC_READ(i1)
         size_read=size_read+SIZE_VARIABLES(i1)
     &        +int(SIZE_GEST(i1),kind=8)
#if defined(MUMPS_NOF2003)
         size_read=size_read
     &        +int(2*id%KEEP(34)*NbRecords(i1),kind=8)
#endif
      elseif(mode.EQ.fake_restore_mode) then
      endif
      ENDDO
 200  continue
      if(mode.EQ.memory_save_mode) then
         WRITTEN_STRUC_SIZE=sum(SIZE_VARIABLES)+
     &         sum(SIZE_VARIABLES_ROOTC)+
     &         SIZE_VARIABLES_BLR+SIZE_VARIABLES_FRONT_DATA+
     &         SIZE_VARIABLES_L0FAC+
#if defined(SINGLEINDOUBLE)
#endif
     &         SIZE_VARIABLES_ROOTA
         TOTAL_STRUC_SIZE=WRITTEN_STRUC_SIZE
     &        +sum(DIFF_SIZE_ALLOC_READ)
     &        +sum(DIFF_SIZE_ALLOC_READ_ROOTC)
         gest_size=sum(SIZE_GEST)+sum(SIZE_GEST_ROOTC)
     &        +SIZE_GEST_BLR+SIZE_GEST_FRONT_DATA
     &        +SIZE_GEST_L0FAC
     &        +SIZE_GEST_ROOTA
     &        +int(5*SIZE_CHARACTER,kind=8)
     &        +int(23*SIZE_CHARACTER,kind=8)
     &        +int(2*SIZE_INT8,kind=8)+int(1,kind=8)
     &        +int(3*SIZE_INT,kind=8)  
     &        +int(SIZE_LOGICAL,kind=8)
#if defined(SINGLEINDOUBLE)
#endif
         IF(associated(id%OOC_FILE_NAME_LENGTH).AND.
     &        associated(id%OOC_FILE_NAMES)) THEN
            gest_size=gest_size+int(SIZE_INT,kind=8)
     &           +int(id%OOC_FILE_NAME_LENGTH(1)*SIZE_CHARACTER,kind=8)
         ELSE
            gest_size=gest_size+int(2*SIZE_INT,kind=8)
         ENDIF
#if defined(MUMPS_NOF2003)
         tot_NbRecords=sum(NbRecords)+sum(NbRecords_ROOTC)+8
         gest_size=gest_size+int(2*id%KEEP(34)*tot_NbRecords,kind=8)
#endif
         TOTAL_FILE_SIZE=WRITTEN_STRUC_SIZE+gest_size
      elseif(mode.EQ.save_mode) then
      elseif(mode.EQ.restore_mode) then
#if ! defined(NOSCALAPACK)
         if(idintr%root%gridinit_done) then
            idintr%root%CNTXT_BLACS = id%COMM_NODES
            CALL blacs_gridinit( idintr%root%CNTXT_BLACS, 'R',
     &           idintr%root%NPROW, idintr%root%NPCOL )
            idintr%root%gridinit_done = .TRUE.
            idintr%root%DESCRIPTOR(2) = idintr%root%CNTXT_BLACS
         endif
#endif
      elseif(mode.EQ.fake_restore_mode) then
      elseif(mode.EQ.restore_ooc_mode) then
      endif
 100  continue
      RETURN
      CONTAINS
      SUBROUTINE MUMPS_SAVERSTR_REALARRAY(idREAL)
      IMPLICIT NONE
      REAL, DIMENSION(:), POINTER :: idREAL
      NbRecords(i1)=2
      if(mode.EQ.memory_save_mode) then
         IF(associated(idREAL)) THEN
            SIZE_GEST(i1)=SIZE_INT
            SIZE_VARIABLES(i1)=size(idREAL,1)*SIZE_RL_OR_DBL
         ELSE
            SIZE_GEST(i1)=SIZE_INT*2
            SIZE_VARIABLES(i1)=0_8
         ENDIF
      elseif(mode.EQ.save_mode) then
         IF(associated(idREAL)) THEN
            write(unit,iostat=err) size(idreal,1)
            if(err.ne.0) then
               id%INFO(1) = -72
               CALL MUMPS_SETI8TOI4(
     &              TOTAL_FILE_SIZE-size_written,
     &              id%INFO(2))
            endif
            write(unit,iostat=err) idREAL
         ELSE
            write(unit,iostat=err) -999
            if(err.ne.0) then
               id%INFO(1) = -72
               CALL MUMPS_SETI8TOI4(
     &         TOTAL_FILE_SIZE-size_written,
     &         id%INFO(2))
            endif
            write(unit,iostat=err) -999
         ENDIF 
         if(err.ne.0) then
               id%INFO(1) = -72
               CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &         id%INFO(2))
         endif
         CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0 ) GOTO 100
      elseif(mode.EQ.restore_mode) then
         nullify(idREAL)
         read(unit,iostat=err) size_array1
         if (err.ne.0) THEN
           id%INFO(1) = -75
           CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                    ,id%INFO(2))
         endif
         CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0 ) GOTO 100
           if (size_array1.EQ.-999) then
                     SIZE_GEST(i1)=SIZE_INT*2
                     SIZE_VARIABLES(i1)=0_8
                     read(unit,iostat=err) dummy
           else
                     SIZE_GEST(i1)=SIZE_INT
                     SIZE_VARIABLES(i1)=size_array1*SIZE_RL_OR_DBL
                     allocate(idREAL(size_array1), stat=allocok)
                     if (allocok .GT. 0) THEN
                        id%INFO(1) = -78
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_STRUC_SIZE-size_allocated
     &                       ,id%INFO(2))
                     endif
                     read(unit,iostat=err) idREAL
           endif
           CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                          id%COMM, id%MYID )
           IF ( id%INFO(1) .LT. 0 ) GOTO 100
           if (err.ne.0) THEN
               id%INFO(1) = -75
               CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                    ,id%INFO(2))
           endif 
           CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                          id%COMM, id%MYID )
           IF ( id%INFO(1) .LT. 0 ) GOTO 100
      endif
 100  CONTINUE
      RETURN
      END SUBROUTINE MUMPS_SAVERSTR_REALARRAY
      SUBROUTINE MUMPS_SAVE_INT(idINT)
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: idINT
      NbRecords(i1)=1
      if(mode.EQ.memory_save_mode) then
         SIZE_VARIABLES(i1)=SIZE_INT
      elseif(mode.EQ.save_mode) then
         write(unit,iostat=err) idINT
         if(err.ne.0) then
           id%INFO(1) = -72
           CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &          id%INFO(2))
         endif
         CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0 ) GOTO 100
      elseif(mode.EQ.restore_mode) then
         SIZE_VARIABLES(i1)=SIZE_INT
         read(unit,iostat=err) idINT
         if(err.ne.0) THEN
            id%INFO(1) = -75
            CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &           ,id%INFO(2))
         endif 
         CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0 ) GOTO 100
      endif
 100  CONTINUE
      RETURN
      END SUBROUTINE MUMPS_SAVE_INT
      SUBROUTINE CMUMPS_SAVE_INT_SHPTR_ARRAY(id_INTPTR
     &  )
      IMPLICIT NONE
      INTEGER, POINTER, DIMENSION(:) :: id_INTPTR
            NbRecords(i1)=2
            if(mode.EQ.memory_save_mode) then
               IF(associated(id_INTPTR)
     &            ) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id_INTPTR,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(mode.EQ.save_mode) then
               IF(associated(id_INTPTR)
     &            ) THEN
                  write(unit,iostat=err) size(id_INTPTR,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id_INTPTR
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(mode.EQ.restore_mode) then
               nullify(id_INTPTR)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
                  if(err.ne.0) THEN
                    id%INFO(1) = -75
                    CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                   ,id%INFO(2))
                  endif
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id_INTPTR(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  else
                    read(unit,iostat=err) id_INTPTR
                    if(err.ne.0) THEN
                      id%INFO(1) = -75
                      CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                     ,id%INFO(2))
                    endif
                  endif
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
 100  CONTINUE
      RETURN
      END SUBROUTINE CMUMPS_SAVE_INT_SHPTR_ARRAY
      END SUBROUTINE CMUMPS_SAVE_RESTORE_STRUCTURE
      SUBROUTINE CMUMPS_SAVE_RESTORE_ROOTA(
     &                 roota
     &                 ,unit,MYID,mode
     &                 ,SIZE_GEST_ROOTA,SIZE_VARIABLES_ROOTA
     &                 ,SIZE_INT, SIZE_INT8, SIZE_ARITH_DEP
     &                 ,SIZE_RL_OR_DBL
     &                 ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                 ,size_read,size_allocated,size_written
     &                 ,INFO)
      USE CMUMPS_INTR_TYPES, ONLY : CMUMPS_ROOT_STRUC
      IMPLICIT NONE
      TYPE (CMUMPS_ROOT_STRUC) :: roota
      INTEGER,intent(IN)::unit,MYID
      INTEGER,intent(IN) :: mode
      INTEGER,INTENT(OUT) :: SIZE_GEST_ROOTA
      INTEGER(8),intent(OUT) :: SIZE_VARIABLES_ROOTA
      INTEGER,intent(INOUT):: INFO(2)
      INTEGER,intent(IN):: SIZE_INT, SIZE_INT8, SIZE_ARITH_DEP
      INTEGER,intent(IN):: SIZE_RL_OR_DBL
      INTEGER(8),intent(IN) :: TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
      INTEGER(8),intent(INOUT):: size_read,size_allocated,size_written
      INTEGER :: SIZE_GEST, i3
      INTEGER(8) :: SIZE_VARIABLE
      INTEGER, PARAMETER :: NBVARIABLES_ROOTA=7
      INTEGER,  PARAMETER :: RA_SINGULAR_VALUES=7
      INTEGER,  PARAMETER :: RA_SVD_VT=6
      INTEGER,  PARAMETER :: RA_SVD_U=5
      INTEGER,  PARAMETER :: RA_RHS_ROOT=4
      INTEGER,  PARAMETER :: RA_QR_TAU=3
      INTEGER,  PARAMETER :: RA_SCHUR_POINTER=2
      INTEGER,  PARAMETER :: RA_RHS_CNTR_MASTER_ROOT=1
      SIZE_GEST_ROOTA = 0
      SIZE_VARIABLES_ROOTA = 0_8
      DO i3 = 1, NBVARIABLES_ROOTA
            SIZE_GEST     = 0
            SIZE_VARIABLE = 0_8
            SELECT CASE(i3)
            CASE(RA_QR_TAU)
               CALL CMUMPS_SAVE_RESTORE_ARRAY_C1D(
     &               roota%QR_TAU )
            CASE(RA_SVD_U)
               CALL CMUMPS_SAVE_RESTORE_ARRAY_2D(roota%SVD_U)
            CASE(RA_SVD_VT)
            CASE(RA_SINGULAR_VALUES)
               CALL CMUMPS_SAVE_RESTORE_ARRAY_R1D(
     &               roota%SINGULAR_VALUES)
            CASE(RA_RHS_CNTR_MASTER_ROOT)
               CALL CMUMPS_SAVE_RESTORE_ARRAY_C1D(
     &               roota%RHS_CNTR_MASTER_ROOT)
            CASE(RA_RHS_ROOT)
            CASE(RA_SCHUR_POINTER)
            CASE DEFAULT
            END SELECT
            IF ( INFO(1) .LT. 0 ) GOTO 100
            IF (mode.EQ.memory_save_mode) then
              SIZE_VARIABLES_ROOTA = SIZE_VARIABLES_ROOTA +
     &                               SIZE_VARIABLE
              SIZE_GEST_ROOTA      = SIZE_GEST_ROOTA + SIZE_GEST
            ENDIF
      END DO 
 100  CONTINUE
      RETURN
      CONTAINS
               SUBROUTINE CMUMPS_SAVE_RESTORE_ARRAY_2D(PTRARRAY2D)
               IMPLICIT NONE
               COMPLEX, DIMENSION(:,:), POINTER :: PTRARRAY2D
               INTEGER :: LocalNbRecords, err, allocok, dummy
               INTEGER :: size_array1, size_array2
               INTEGER(8) :: sz
               LocalNbRecords = 2
               if(mode.EQ.memory_save_mode) then
                  IF(associated(PTRARRAY2D)) THEN
                     SIZE_GEST = SIZE_INT*2
                     SIZE_VARIABLE = size(PTRARRAY2D,1)
     &                    *size(PTRARRAY2D,2)*SIZE_ARITH_DEP
                  ELSE
                     SIZE_GEST = SIZE_INT*3
                     SIZE_VARIABLE = 0_8
                  ENDIF
               elseif(mode.EQ.save_mode) then
                 IF(associated(PTRARRAY2D)) THEN
                    write(unit,iostat=err) size(PTRARRAY2D,1)
     &                   ,size(PTRARRAY2D,2)
                 ELSE
                    write(unit,iostat=err) -999,-998
                 ENDIF
                 if(err.ne.0) then
                   INFO(1) = -72
                   CALL MUMPS_SETI8TOI4(
     &                  TOTAL_FILE_SIZE-size_written
     &                  ,INFO(2))
                 else
                   size_written=size_written+2*SIZE_INT
                 endif
                 IF ( INFO(1) .LT. 0 ) GOTO 100
                 IF(associated(PTRARRAY2D)) THEN
                   write(unit,iostat=err) PTRARRAY2D
                   sz= int(size(PTRARRAY2D,1),8) *
     &                 int(size(PTRARRAY2D,2),8) *
     &                 SIZE_ARITH_DEP
                 ELSE
                   write(unit,iostat=err) -999
                   sz=SIZE_INT
                 ENDIF
                 if(err.ne.0) then
                    INFO(1) = -72
                    CALL MUMPS_SETI8TOI4(
     &                    TOTAL_FILE_SIZE-size_written
     &                    ,INFO(2))
                 else
                   size_written=size_written+sz
                 endif
                 IF ( INFO(1) .LT. 0 ) GOTO 100
               elseif(mode.EQ.restore_mode) then
                  nullify(PTRARRAY2D)
                  read(unit,iostat=err) size_array1,size_array2
                  if(err.ne.0) THEN
                     INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                    ,INFO(2))
                  else
                    size_read = size_read+2*SIZE_INT
                    size_allocated = size_allocated + 2*SIZE_INT8
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  if(size_array1.EQ.-999) then
                     read(unit,iostat=err) dummy
                     sz = SIZE_INT
                  else
                     sz = int(size_array1,8)*int(size_array2,8)
     &                               * SIZE_ARITH_DEP
                     allocate(PTRARRAY2D(size_array1,
     &                    size_array2),
     &                    stat=allocok)
                     if (allocok .GT. 0) THEN
                        INFO(1) = -78
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_STRUC_SIZE-size_allocated
     &                       ,INFO(2))
                     endif
                     read(unit,iostat=err) PTRARRAY2D
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  if(err.ne.0) THEN
                     INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                    , INFO(2))
                  else
                    size_read = size_read + sz
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
               endif
#if defined(MUMPS_NOF2003)
               IF (mode.EQ.memory_save_mode) THEN
                 SIZE_GEST = SIZE_GEST+2*SIZE_INT*LocalNbRecords
               ELSE IF (mode.EQ.save_mode) THEN
                 size_written = size_written+2*SIZE_INT*LocalNbRecords
               ELSE IF (mode.EQ.restore_mode) THEN
                 size_read = size_read+2*SIZE_INT*LocalNbRecords
               ENDIF
#endif
  100          CONTINUE
               RETURN
            END SUBROUTINE CMUMPS_SAVE_RESTORE_ARRAY_2D
            SUBROUTINE CMUMPS_SAVE_RESTORE_ARRAY_C1D(PTRARRAY1D)
               IMPLICIT NONE
               COMPLEX, DIMENSION(:), POINTER :: PTRARRAY1D
               INTEGER :: LocalNbRecords, err, allocok, dummy
               INTEGER :: size_array1
               INTEGER(8) :: sz
               LocalNbRecords = 2
               if(mode.EQ.memory_save_mode) then
                  IF(associated(PTRARRAY1D)) THEN
                     SIZE_GEST = SIZE_INT
                     SIZE_VARIABLE = size(PTRARRAY1D)*SIZE_ARITH_DEP
                  ELSE
                     SIZE_GEST = SIZE_INT*2
                     SIZE_VARIABLE = 0_8
                  ENDIF
               elseif(mode.EQ.save_mode) then
                 IF(associated(PTRARRAY1D)) THEN
                    write(unit,iostat=err) size(PTRARRAY1D)
                 ELSE
                    write(unit,iostat=err) -999
                 ENDIF
                 if(err.ne.0) then
                   INFO(1) = -72
                   CALL MUMPS_SETI8TOI4(
     &                  TOTAL_FILE_SIZE-size_written
     &                  ,INFO(2))
                 else
                   size_written=size_written+SIZE_INT
                 endif
                 IF ( INFO(1) .LT. 0 ) GOTO 100
                 IF(associated(PTRARRAY1D)) THEN
                   write(unit,iostat=err) PTRARRAY1D
                   sz= int(size(PTRARRAY1D),8)*
     &                 SIZE_ARITH_DEP
                 ELSE
                   write(unit,iostat=err) -999
                   sz=SIZE_INT
                 ENDIF
                 if(err.ne.0) then
                    INFO(1) = -72
                    CALL MUMPS_SETI8TOI4(
     &                    TOTAL_FILE_SIZE-size_written
     &                    ,INFO(2))
                 else
                   size_written=size_written+sz
                 endif
                 IF ( INFO(1) .LT. 0 ) GOTO 100
               elseif(mode.EQ.restore_mode) then
                  nullify(PTRARRAY1D)
                  read(unit,iostat=err) size_array1
                  if(err.ne.0) THEN
                     INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                    ,INFO(2))
                  else
                    size_read = size_read+SIZE_INT
                    size_allocated = size_allocated + SIZE_INT
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  if(size_array1.EQ.-999) then
                     read(unit,iostat=err) dummy
                     sz = SIZE_INT
                  else
                     sz = int(size_array1,8) * SIZE_ARITH_DEP
                     allocate(PTRARRAY1D(size_array1),
     &                    stat=allocok)
                     if (allocok .GT. 0) THEN
                        INFO(1) = -78
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_STRUC_SIZE-size_allocated
     &                       ,INFO(2))
                     endif
                     read(unit,iostat=err) PTRARRAY1D
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  if(err.ne.0) THEN
                     INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                    , INFO(2))
                  else
                    size_read = size_read + sz
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
               endif
#if defined(MUMPS_NOF2003)
               IF (mode.EQ.memory_save_mode) THEN
                 SIZE_GEST = SIZE_GEST+2*SIZE_INT*LocalNbRecords
               ELSE IF (mode.EQ.save_mode) THEN
                 size_written = size_written+2*SIZE_INT*LocalNbRecords
               ELSE IF (mode.EQ.restore_mode) THEN
                 size_read = size_read+2*SIZE_INT*LocalNbRecords
               ENDIF
#endif
  100          CONTINUE
               RETURN
            END SUBROUTINE CMUMPS_SAVE_RESTORE_ARRAY_C1D
            SUBROUTINE CMUMPS_SAVE_RESTORE_ARRAY_R1D(PTRARRAY1D)
               IMPLICIT NONE
               REAL, DIMENSION(:), POINTER :: PTRARRAY1D
               INTEGER :: LocalNbRecords, err, allocok, dummy
               INTEGER :: size_array1
               INTEGER(8) :: sz
               LocalNbRecords = 2
               if(mode.EQ.memory_save_mode) then
                  IF(associated(PTRARRAY1D)) THEN
                     SIZE_GEST = SIZE_INT
                     SIZE_VARIABLE = size(PTRARRAY1D)*SIZE_RL_OR_DBL
                  ELSE
                     SIZE_GEST = SIZE_INT*2
                     SIZE_VARIABLE = 0_8
                  ENDIF
               elseif(mode.EQ.save_mode) then
                 IF(associated(PTRARRAY1D)) THEN
                    write(unit,iostat=err) size(PTRARRAY1D)
                 ELSE
                    write(unit,iostat=err) -999
                 ENDIF
                 if(err.ne.0) then
                   INFO(1) = -72
                   CALL MUMPS_SETI8TOI4(
     &                  TOTAL_FILE_SIZE-size_written
     &                  ,INFO(2))
                 else
                   size_written=size_written+SIZE_INT
                 endif
                 IF ( INFO(1) .LT. 0 ) GOTO 100
                 IF(associated(PTRARRAY1D)) THEN
                   write(unit,iostat=err) PTRARRAY1D
                   sz= int(size(PTRARRAY1D),8)*
     &                 SIZE_RL_OR_DBL
                 ELSE
                   write(unit,iostat=err) -999
                   sz=SIZE_INT
                 ENDIF
                 if(err.ne.0) then
                    INFO(1) = -72
                    CALL MUMPS_SETI8TOI4(
     &                    TOTAL_FILE_SIZE-size_written
     &                    ,INFO(2))
                 else
                   size_written=size_written+sz
                 endif
                 IF ( INFO(1) .LT. 0 ) GOTO 100
               elseif(mode.EQ.restore_mode) then
                  nullify(PTRARRAY1D)
                  read(unit,iostat=err) size_array1
                  if(err.ne.0) THEN
                     INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                    ,INFO(2))
                  else
                    size_read = size_read+SIZE_INT
                    size_allocated = size_allocated + SIZE_INT
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  if(size_array1.EQ.-999) then
                     read(unit,iostat=err) dummy
                     sz = SIZE_INT
                  else
                     sz = int(size_array1,8) * SIZE_RL_OR_DBL
                     allocate(PTRARRAY1D(size_array1),
     &                    stat=allocok)
                     if (allocok .GT. 0) THEN
                        INFO(1) = -78
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_STRUC_SIZE-size_allocated
     &                       ,INFO(2))
                     endif
                     read(unit,iostat=err) PTRARRAY1D
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  if(err.ne.0) THEN
                     INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                    , INFO(2))
                  else
                    size_read = size_read + sz
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
               endif
#if defined(MUMPS_NOF2003)
               IF (mode.EQ.memory_save_mode) THEN
                 SIZE_GEST = SIZE_GEST+2*SIZE_INT*LocalNbRecords
               ELSE IF (mode.EQ.save_mode) THEN
                 size_written = size_written+2*SIZE_INT*LocalNbRecords
               ELSE IF (mode.EQ.restore_mode) THEN
                 size_read = size_read+2*SIZE_INT*LocalNbRecords
               ENDIF
#endif
  100          CONTINUE
               RETURN
            END SUBROUTINE CMUMPS_SAVE_RESTORE_ARRAY_R1D
      END SUBROUTINE CMUMPS_SAVE_RESTORE_ROOTA
      END MODULE CMUMPS_SAVE_RESTORE
#else
      SUBROUTINE CMUMPS_SAVE_RESTORE_RETURN()
      RETURN
      END SUBROUTINE CMUMPS_SAVE_RESTORE_RETURN
#endif
