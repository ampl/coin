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
      MODULE DMUMPS_CONFIG_FILE
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      INTEGER :: LEN_CONFIG_FILE_NAME
      PARAMETER( LEN_CONFIG_FILE_NAME = 1024 )
      CONTAINS
      SUBROUTINE DMUMPS_CONFIG_FROM_FILE(id)
      INCLUDE 'mpif.h'
      CHARACTER(len=LEN_CONFIG_FILE_NAME):: CONFIG_FILE_NAME
      CHARACTER(len=1096)::LINE   
      CHARACTER:: CHARCOM
      INTEGER :: LINE_size,IERR
      CHARACTER(len=15)::PARAM=' '
      CHARACTER(len=10)::INDICE=' '
      CHARACTER(len=1024)::VAL=' ' 
      INTEGER :: LENVAL
      INTEGER :: unit_config
      INTEGER :: i1,ind_first_par,ind_second_par,ind_equal,line_number
      INTEGER :: LP,ind_end
      LOGICAL :: LPOK, IS_CONFIG_FILE_OPEN
      TYPE(DMUMPS_STRUC) id
      LP    = id%ICNTL(1)
      LPOK  = (LP.GT.0) .AND. (id%ICNTL(4).GE.1)
      IS_CONFIG_FILE_OPEN = .FALSE.
      IF ( id%MYID.eq.0 ) THEN
         CALL DMUMPS_GET_CONFIG_FILE(CONFIG_FILE_NAME,id%INFO)
         IF (id%INFO(1) .LT. 0) GOTO 100
         if (CONFIG_FILE_NAME.EQ."NAME_NOT_INITIALIZED") goto 100
         CALL MUMPS_FIND_UNIT(unit_config)
         IF( unit_config .LT. 0 ) THEN
            id%INFO(1) = -79
            id%INFO(2) = 3
            IF(LPOK) THEN
             write(LP,*) 
     &          " ERROR during configuration file feature:", 
     &          " no free file unit could be found "
            ENDIF
            goto 100
         ENDIF
         OPEN(FILE=trim(adjustl(CONFIG_FILE_NAME)), UNIT=unit_config,
     &        STATUS='OLD',iostat=IERR)
         if(IERR.NE.0) THEN
            id%INFO(1)=-63
            IF(LPOK) THEN
               write(LP,*) 
     &          " ERROR during configuration file feature:", 
     &          " cannot open file: ",
     &          trim(adjustl(CONFIG_FILE_NAME))
            ENDIF
            id%INFO(2)=-1
            goto 100
         endif
         IS_CONFIG_FILE_OPEN = .TRUE.
         line_number=0
         do
 999        read(unit_config,'(a)',iostat=IERR) LINE
            if(IERR .LT. 0 ) exit
            line_number=line_number+1
            if(IERR.GT.0) THEN
               id%INFO(1)=-64
               id%INFO(2)=line_number
               IF(LPOK) THEN
                  write(LP,*) 
     &          " ERROR during configuration file feature:", 
     &          " cannot read line ",line_number,
     &          " of file: ",trim(adjustl(CONFIG_FILE_NAME))
               ENDIF
               goto 100
            endif
            CHARCOM=LINE(1:1)
            if(CHARCOM.eq."#") goto 999
            LINE_size=len_trim(LINE)
            IF (LINE_size.EQ.0) cycle
            ind_first_par=-999
            ind_second_par=-999
            ind_equal=-999
            ind_end=LINE_size
            i1=1
            do while (i1.LE.LINE_size)
               if(LINE(i1:i1).EQ."(") then
                  ind_first_par=i1
               elseif(LINE(i1:i1).EQ.")") then
                  ind_second_par=i1
               elseif(LINE(i1:i1).EQ."=") then
                  ind_equal=i1
               elseif(LINE(i1:i1).EQ."#") then
                  ind_end=i1-1
                  goto 998
               endif
               i1=i1+1
            enddo
 998        continue
            if(ind_equal.eq.-999) then
               id%INFO(1)=-65
               id%INFO(2)=line_number
               IF(LPOK) THEN
                  write(LP,*) 
     &            " ERROR during configuration file feature:", 
     &            " syntax error in ",
     &            trim(adjustl(CONFIG_FILE_NAME))
     &            ," at line ",line_number," cannot find ="
               ENDIF
               goto 100
            endif
            if((ind_first_par.eq.-999).and.
     &           (ind_second_par.eq.-999)) then
               PARAM=trim(adjustl(LINE(1:ind_equal-1)))
               INDICE="NULLINDICE"
               VAL=trim(adjustl(LINE(ind_equal+1:ind_end)))
               LENVAL=len_trim(VAL)
            elseif(ind_first_par.LT.ind_second_par) then
               PARAM=trim(adjustl(LINE(1:ind_first_par-1)))
               INDICE=trim(adjustl(
     &              LINE(ind_first_par+1:ind_second_par-1)))
               VAL=trim(adjustl(LINE(ind_equal+1:ind_end)))
               LENVAL=len_trim(VAL)
            else
               id%INFO(1)=-65
               id%INFO(2)=line_number
               IF(LPOK) THEN
                  write(LP,*) 
     &            " ERROR during configuration file feature:", 
     &            " syntax error in ",
     &            trim(adjustl(CONFIG_FILE_NAME))
     &            ," at line ",line_number," problem with",
     &            " parenthesis"
               ENDIF
               goto 100
            endif
            call UPPERCASE(PARAM)
            CALL MPI_BCAST(PARAM,15,MPI_CHARACTER,0,id%COMM,IERR)
            CALL MPI_BCAST(INDICE,10,MPI_CHARACTER,0,id%COMM,IERR)
            CALL MPI_BCAST(LENVAL,1,MPI_INTEGER,0,id%COMM,IERR)
            CALL MPI_BCAST(VAL,LENVAL,MPI_CHARACTER,0,id%COMM,IERR)
            CALL DMUMPS_config_param(id,PARAM,INDICE,VAL(1:LENVAL),IERR)
            if (IERR.NE.0) THEN
               id%INFO(1)=IERR
               IF (IERR.EQ.-79) THEN
                 id%INFO(2)=4
                 IF(LPOK) THEN
                 write(LP,*) 
     &           " ERROR during configuration file feature:", 
     &           " a free file unit could be found for an",
     &           " auxiliary file"
                 ENDIF
               ELSE
                 id%INFO(2)=line_number
                 IF(LPOK) THEN
                 write(LP,*) 
     &           " ERROR ",IERR," during configuration file feature:", 
     &           " error occurred at line ",line_number," of ",
     &           trim(adjustl(CONFIG_FILE_NAME))
                 ENDIF
               ENDIF
               goto 101
            endif
         enddo
 100     continue         
         PARAM="ENDOFCONFIGFILE"
         CALL MPI_BCAST(PARAM,15,MPI_CHARACTER,0
     &        ,id%COMM,IERR)
 101     continue        
         IF (IS_CONFIG_FILE_OPEN) THEN
           close(unit_config)
         ENDIF
      ELSE
         DO 
            CALL MPI_BCAST(PARAM, 15,MPI_CHARACTER,0,id%COMM,IERR)
            IF(trim(adjustl(PARAM)).EQ."ENDOFCONFIGFILE") goto 200
            CALL MPI_BCAST(INDICE, 10,MPI_CHARACTER,0,id%COMM,IERR)
            CALL MPI_BCAST(LENVAL,1,MPI_INTEGER,0,id%COMM,IERR)
            CALL MPI_BCAST(VAL,LENVAL,MPI_CHARACTER,0,id%COMM,IERR)
            CALL DMUMPS_config_param(id,PARAM,INDICE,VAL(1:LENVAL),IERR)
            if(IERR.NE.0) THEN
               goto 200
            endif
         ENDDO
 200     continue
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      RETURN
      END SUBROUTINE DMUMPS_CONFIG_FROM_FILE
      SUBROUTINE DMUMPS_config_param(id,PARAM,INDICE,VAL,IERR)
      TYPE(DMUMPS_STRUC) id
      CHARACTER(len=15)::PARAM
      CHARACTER(len=10)::INDICE
      CHARACTER(len=*)::VAL
      INTEGER :: IERR
      INTEGER :: IERR_MPI
      INCLUDE 'mpif.h'
      IERR=0
      select case(trim(adjustl(PARAM)))
      case("KEEP")
         call DMUMPS_KEEP_config(id,INDICE,VAL,IERR)
      case("DKEEP")
         call DMUMPS_DKEEP_config(id,INDICE,VAL,IERR)
      case("KEEP8")
         call DMUMPS_KEEP8_config(id,INDICE,VAL,IERR)
      case("ICNTL")       
         call DMUMPS_ICNTL_config(id,INDICE,VAL,IERR)
      case("CNTL")
         call DMUMPS_CNTL_config(id,INDICE,VAL,IERR)
      case("METIS_OPTIONS")
         call DMUMPS_METIS_OPTIONS_config(id,INDICE,VAL,IERR)
      case("WRITE_PROBLEM")
         call DMUMPS_WRITE_PROBLEM_config(id,VAL)
      case("BLKPTR_FILE")
         IF (id%MYID.EQ.0) THEN
           call DMUMPS_READ_BLKPTR(id,VAL,IERR)
         ENDIF
         CALL MPI_BCAST(IERR,1,MPI_INTEGER,0,id%COMM,IERR_MPI)
      case("BLKVAR_FILE")
         IF (id%MYID.EQ.0) THEN
           call DMUMPS_READ_BLKVAR(id,VAL,IERR)
         ENDIF
         CALL MPI_BCAST(IERR,1,MPI_INTEGER,0,id%COMM,IERR_MPI)
      case default
         IERR=-65 
         RETURN
      end select
      RETURN
      END SUBROUTINE DMUMPS_config_param
      SUBROUTINE DMUMPS_READ_BLKPTR(id,VAL,IERR)
      TYPE(DMUMPS_STRUC) :: id
      CHARACTER(len=*)   :: VAL
      INTEGER            :: IERR
      INTEGER :: unit_blk, I
      IF (associated(id%BLKPTR)) THEN
        IERR = -67
        GOTO 100
      ENDIF
      CALL MUMPS_FIND_UNIT(unit_blk)
      IF( unit_blk .LT. 0 ) THEN
         IERR = -79
         GOTO 100
      ENDIF
      OPEN(FILE=trim(adjustl(VAL)), UNIT=unit_blk,
     &        STATUS='OLD',iostat=IERR)
      IF ( IERR .NE. 0 ) THEN
         IERR=-66
         goto 100 
      ENDIF
      READ(UNIT=unit_blk,FMT=*,IOSTAT=IERR) id%NBLK
      IF ( IERR .NE. 0 ) THEN
        IERR = -66
        GOTO 110 
      ENDIF
      ALLOCATE(id%BLKPTR(id%NBLK+1),stat=IERR)
      IF (IERR.GT.0) THEN
        IERR=-67
        GOTO 110 
      ENDIF
      DO I=1, id%NBLK + 1
        READ(UNIT=unit_blk,FMT=*,IOSTAT=IERR) id%BLKPTR(I)
        IF (IERR .NE. 0) THEN
          IERR = -66
          GOTO 120
        ENDIF
      ENDDO
      id%KEEP(145)=1
      RETURN
 120  CONTINUE
      DEALLOCATE(id%BLKPTR)
      NULLIFY(id%BLKPTR)
 110  CONTINUE
      CLOSE(unit_blk)
 100  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_READ_BLKPTR
      SUBROUTINE DMUMPS_READ_BLKVAR(id,VAL,IERR)
      TYPE(DMUMPS_STRUC) :: id
      CHARACTER(len=*)   :: VAL
      INTEGER            :: IERR
      INTEGER :: unit_blk, I
      IF (associated(id%BLKVAR)) THEN
        IERR = -67
        GOTO 100
      ENDIF
      CALL MUMPS_FIND_UNIT(unit_blk)
      IF( unit_blk .LT. 0 ) THEN
         IERR = -79
         GOTO 100
      ENDIF
      OPEN(FILE=trim(adjustl(VAL)), UNIT=unit_blk,
     &        STATUS='OLD',iostat=IERR)
      IF ( IERR .NE. 0 ) THEN
         IERR=-66
         goto 100 
      ENDIF
      ALLOCATE(id%BLKVAR(id%N),stat=IERR)
      IF (IERR.GT.0) THEN
        IERR=-67
        GOTO 110 
      ENDIF
      DO I=1, id%N
        READ(UNIT=unit_blk,FMT=*,IOSTAT=IERR) id%BLKVAR(I)
        IF (IERR .NE. 0) THEN
          IERR = -66
          GOTO 120
        ENDIF
      ENDDO
      id%KEEP(146)=1
      RETURN
 120  CONTINUE
      DEALLOCATE(id%BLKVAR)
      NULLIFY(id%BLKVAR)
 110  CONTINUE
      CLOSE(unit_blk)
 100  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_READ_BLKVAR
      SUBROUTINE DMUMPS_config_exit(id)
      TYPE(DMUMPS_STRUC) :: id
      IF (id%KEEP(145).EQ.1) THEN
        IF ( associated(id%BLKPTR) ) THEN
          DEALLOCATE(id%BLKPTR)
        ELSE
          WRITE(*,*) "Internal error 1 in DMUMPS_config_exit,",
     &    id%KEEP(145)
          CALL MUMPS_ABORT()
        ENDIF
        NULLIFY(id%BLKPTR)
        id%KEEP(145)=0
      ENDIF
      IF (id%KEEP(146).EQ.1) THEN
        IF ( associated(id%BLKVAR) ) THEN
          DEALLOCATE(id%BLKVAR)
        ELSE
          WRITE(*,*) "Internal error 1 in DMUMPS_config_exit,",
     &    id%KEEP(146)
          CALL MUMPS_ABORT()
        ENDIF
        NULLIFY(id%BLKVAR)
        id%KEEP(146)=0
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_config_exit
      SUBROUTINE DMUMPS_KEEP_config(id,INDICE,VAL, IERR)
      TYPE(DMUMPS_STRUC) :: id
      CHARACTER(len=10)  :: INDICE
      CHARACTER(len=*)   :: VAL
      INTEGER            :: IERR
      INTEGER::ind
      call get_indice_from_line(INDICE,ind,IERR)
      IF (IERR .NE. 0) GOTO 100
      call Set_INT_value(id%KEEP(ind),VAL,IERR)
      IF (IERR .NE. 0) GOTO 100
 100  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_KEEP_config
      SUBROUTINE DMUMPS_DKEEP_config(id,INDICE,VAL,IERR)
      TYPE(DMUMPS_STRUC) :: id
      CHARACTER(len=10)  :: INDICE
      CHARACTER(len=*) :: VAL
      INTEGER            :: IERR
      INTEGER::ind
      call get_indice_from_line(INDICE,ind,IERR)
      IF (IERR .NE. 0) GOTO 100
      call Set_REAL_value(id%DKEEP(ind),VAL,IERR)
      IF (IERR .NE. 0) GOTO 100
 100  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_DKEEP_config
      SUBROUTINE DMUMPS_KEEP8_config(id,INDICE,VAL,IERR)
      TYPE(DMUMPS_STRUC) :: id
      CHARACTER(len=10)  :: INDICE
      CHARACTER(len=*)   :: VAL
      INTEGER            :: IERR
      INTEGER::ind
      call get_indice_from_line(INDICE,ind,IERR)
      IF (IERR .NE. 0) GOTO 100
      call Set_INT64_value(id%KEEP8(ind),VAL,IERR)
      IF (IERR .NE. 0) GOTO 100
 100  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_KEEP8_config
      SUBROUTINE DMUMPS_ICNTL_config(id,INDICE,VAL,IERR)
      TYPE(DMUMPS_STRUC) :: id
      CHARACTER(len=10)  :: INDICE
      CHARACTER(len=*)   :: VAL
      INTEGER            :: IERR
      INTEGER::ind
      call get_indice_from_line(INDICE,ind,IERR)
      IF (IERR .NE. 0) GOTO 100
      call Set_INT_value(id%ICNTL(ind),VAL,IERR)
      IF (IERR .NE. 0) GOTO 100
 100  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_ICNTL_config
      SUBROUTINE DMUMPS_CNTL_config(id,INDICE,VAL,IERR)
      TYPE(DMUMPS_STRUC) :: id
      CHARACTER(len=10)  :: INDICE
      CHARACTER(len=*)   :: VAL
      INTEGER            :: IERR
      INTEGER::ind
      call get_indice_from_line(INDICE,ind,IERR)
      IF (IERR .NE. 0) GOTO 100
      call Set_REAL_value(id%CNTL(ind),VAL,IERR)
      IF (IERR .NE. 0) GOTO 100
 100  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_CNTL_config
      SUBROUTINE DMUMPS_WRITE_PROBLEM_config(id,VAL)
      TYPE(DMUMPS_STRUC) id
      CHARACTER(len=*) :: VAL
      id%WRITE_PROBLEM=trim(adjustl(VAL))
      END SUBROUTINE DMUMPS_WRITE_PROBLEM_config
      SUBROUTINE DMUMPS_METIS_OPTIONS_config(id,INDICE,VAL,IERR)
      TYPE(DMUMPS_STRUC) :: id
      INTEGER            :: IERR
      CHARACTER(len=10)  :: INDICE
      CHARACTER(len=*)   :: VAL
      INTEGER::ind
      call get_indice_from_line(INDICE,ind, IERR)
      IF (IERR .NE. 0) GOTO 100
      call Set_INT_value(id%METIS_OPTIONS(ind),VAL, IERR)
      IF (IERR .NE. 0) GOTO 100
 100  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_METIS_OPTIONS_config
      SUBROUTINE DMUMPS_GET_CONFIG_FILE(CONFIG_FILE_NAME,INFO)
      CHARACTER(len=LEN_CONFIG_FILE_NAME),intent(out):: CONFIG_FILE_NAME
      INTEGER::len_config_file_C
      CHARACTER(len=LEN_CONFIG_FILE_NAME):: tmp_configfile
      INTEGER :: INFO(2)
      CONFIG_FILE_NAME=''
      tmp_configfile=''
      call mumps_get_config_file_C(len_config_file_C,tmp_configfile)
      if(tmp_configfile(1:len_config_file_C).EQ."NAME_NOT_INITIALIZED")
     &     then
         CONFIG_FILE_NAME="NAME_NOT_INITIALIZED"
      else if (len_config_file_C .GT. LEN_CONFIG_FILE_NAME) then
        INFO(1)=-63
        INFO(2)=LEN_CONFIG_FILE_NAME
      else
         CONFIG_FILE_NAME=
     &        trim(adjustl(tmp_configfile(1:len_config_file_C)))
      endif
      RETURN
      END SUBROUTINE DMUMPS_GET_CONFIG_FILE
      SUBROUTINE UPPERCASE(STR)
      CHARACTER(LEN=*), INTENT(IN OUT) :: STR
      INTEGER :: I, DEL
      DEL = iachar('a') - iachar('A')
      DO I = 1, len_trim(STR)
         IF (lge(STR(I:I),'a') .AND. lle(STR(I:I),'z')) THEN
            STR(I:I) = achar(iachar(STR(I:I)) - DEL)
         END IF
      END DO
      RETURN
      END SUBROUTINE UPPERCASE
      SUBROUTINE get_indice_from_line(string,ind,IERR)
      INTEGER           :: ind
      CHARACTER(len=10) :: string
      INTEGER           :: IERR
      read (string,*,iostat=IERR) ind
      IF (IERR .NE. 0) THEN
        IERR = -65
      ENDIF
      RETURN
      END SUBROUTINE get_indice_from_line      
      SUBROUTINE SET_INT_VALUE(idind,string_value,IERR)
      INTEGER            :: idind
      CHARACTER(len=*)   :: string_value
      INTEGER            :: IERR
      INTEGER:: VAL
      read (string_value,*,IOSTAT=IERR) VAL
      IF (IERR .NE. 0) THEN
        IERR = -65
      ELSE
        idind=VAL
      ENDIF
      RETURN
      END SUBROUTINE SET_INT_VALUE
      SUBROUTINE SET_INT64_VALUE(idind,string_value,IERR)
      INTEGER(8)         :: idind
      CHARACTER(len=*)   :: string_value
      INTEGER            :: IERR
      INTEGER(8):: val
      read (string_value,*,iostat=IERR) val
      IF (IERR .NE. 0) THEN
        IERR = -65
      ELSE
        idind=VAL
      ENDIF
      RETURN
      END SUBROUTINE SET_INT64_VALUE
      SUBROUTINE SET_REAL_VALUE(idind,string_value,IERR)
      DOUBLE PRECISION               :: idind
      CHARACTER(len=*)   :: string_value
      INTEGER            :: IERR
      DOUBLE PRECISION:: val
      READ (string_value,*,IOSTAT=IERR) val
      IF (IERR .NE. 0) THEN
        IERR = -65
      ELSE
        idind=VAL
      ENDIF
      RETURN
      END SUBROUTINE SET_REAL_VALUE
      END MODULE DMUMPS_CONFIG_FILE
      SUBROUTINE DMUMPS_CONFIG_FILE_RETURN()
      RETURN
      END SUBROUTINE DMUMPS_CONFIG_FILE_RETURN
