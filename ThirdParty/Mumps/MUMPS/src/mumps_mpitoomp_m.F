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
#if defined(MPI_TO_K_OMP)
      MODULE MUMPS_MPITOOMP_M
      IMPLICIT NONE
      TYPE PROCS_MAP_STRUC_T
         INTEGER, POINTER :: PROCS_MAP(:)
      END TYPE PROCS_MAP_STRUC_T
      PUBLIC:: MUMPS_MPITOOMP_CHECK,
     &         MUMPS_MPITOOMP_PRE_WORK_LIST,
     &         MUMPS_MPITOOMP_GET_WORK_LIST,
     &         MUMPS_MPITOOMP_PRE_WORK_MAPS,
     &         MUMPS_MPITOOMP_GET_WORK_MAPS,
     &         MUMPS_MPITOOMP_SET_MASK
      INTEGER:: TAG_CORES_LIST_PART
      PARAMETER ( TAG_CORES_LIST_PART = 50 )
      INTEGER:: TAG_AGGREGATED_MAP
      PARAMETER ( TAG_AGGREGATED_MAP = 60 )
      INTEGER:: MPITOOMP_MASTER
      PARAMETER ( MPITOOMP_MASTER = 0 )
      CONTAINS
      SUBROUTINE MUMPS_MPITOOMP_CHECK(COMM, PAR, K414, MTKO_K, NPROCS,
     &                                ERROR)
#if defined(_OPENMP)
      USE OMP_LIB
#endif
      INCLUDE 'mpif.h'
      INTEGER, intent(in):: COMM
      INTEGER, intent(in):: PAR, K414, MTKO_K, NPROCS
      INTEGER, intent(out):: ERROR
      INTEGER :: NOMP
#if defined(_OPENMP)
      INTEGER:: MAX_NOMP, MIN_NOMP
      INTEGER:: MAX_MODULO, MIN_MODULO
#endif
      INTEGER:: IERR
      INTEGER:: MYID
      ERROR = 0
      call MPI_COMM_RANK(COMM, MYID, IERR)
      call MUMPS_CHECK_OMP_C(NOMP)
      IF (NOMP.EQ.0) THEN
        ERROR = -2
      ELSE
#if !defined(_OPENMP)
      ERROR = -2
#else
      call MPI_ALLREDUCE(mod(K414,MTKO_K), MAX_MODULO, 1, MPI_INTEGER,
     &                MPI_MAX, COMM, IERR)
      call MPI_ALLREDUCE(mod(K414,MTKO_K), MIN_MODULO, 1, MPI_INTEGER,
     &                MPI_MIN, COMM, IERR)
      IF (MAX_MODULO .EQ. MIN_MODULO .AND. MAX_MODULO.EQ.0) THEN
        NOMP = 1
!$      NOMP = OMP_GET_MAX_THREADS()
        call MPI_ALLREDUCE(NOMP, MAX_NOMP, 1, MPI_INTEGER, MPI_MAX,
     &                     COMM, IERR)
        call MPI_ALLREDUCE(NOMP, MIN_NOMP, 1, MPI_INTEGER, MPI_MIN,
     &                     COMM, IERR)
        IF (MAX_NOMP .NE. MIN_NOMP) THEN
          ERROR = -3
        ENDIF
      ELSE
        ERROR = -1
      ENDIF
#endif
      ENDIF
      END SUBROUTINE MUMPS_MPITOOMP_CHECK
      SUBROUTINE MUMPS_MPITOOMP_PRE_WORK_LIST(COMM, MTKO_K, NPROCS,
     &           ONLY_ALLOC, NWORKERS,
     &           ERROR)
      INCLUDE 'mpif.h'
      INTEGER, intent(in):: COMM
      INTEGER, intent(in):: MTKO_K, NPROCS
      LOGICAL, intent(in):: ONLY_ALLOC
      INTEGER, intent(out):: NWORKERS
      INTEGER, intent(out):: ERROR
      INTEGER:: MAX_NWORKERS, MIN_NWORKERS
      INTEGER:: IERR
      ERROR = 0
      NWORKERS = NPROCS/MTKO_K
      IF (ONLY_ALLOC) THEN
        MAX_NWORKERS = NWORKERS
        MIN_NWORKERS = NWORKERS
      ELSE
        call MPI_ALLREDUCE(NWORKERS, MAX_NWORKERS, 1, MPI_INTEGER,
     &                MPI_MAX, COMM, IERR)
        call MPI_ALLREDUCE(NWORKERS, MIN_NWORKERS, 1, MPI_INTEGER,
     &                MPI_MIN, COMM, IERR)
      ENDIF
      IF (MAX_NWORKERS .ne. MIN_NWORKERS) THEN
        ERROR = -2
      ENDIF
      END SUBROUTINE MUMPS_MPITOOMP_PRE_WORK_LIST
      SUBROUTINE MUMPS_MPITOOMP_GET_WORK_LIST(COMM, MTKO_K, NPROCS,
     &                NWORKERS, WORKING_LIST, PROCS_MAP, ERROR)
      INTEGER, intent(in):: COMM
      INTEGER, intent(in):: MTKO_K, NPROCS, NWORKERS
      INTEGER, DIMENSION(0:NWORKERS-1), intent(out):: WORKING_LIST
      INTEGER, DIMENSION(:), POINTER:: PROCS_MAP
      INTEGER, intent(out):: ERROR
      INTEGER:: I
      INTEGER:: total
      call MUMPS_MPITOOMP_FILLIN_PROCS_MAP(COMM, MTKO_K, NPROCS,
     &          PROCS_MAP, ERROR)
      total = 0
      DO I=0, NPROCS-1
        IF (PROCS_MAP(I) .eq. I) THEN
          WORKING_LIST(total) = I
          total = total + 1
        ENDIF
      END DO
      END SUBROUTINE MUMPS_MPITOOMP_GET_WORK_LIST
#if defined(MTKO_ADVANCED)
      SUBROUTINE MUMPS_MPITOOMP_INSERT_ID_SORTED(VALUES_SIZE, VALUES,
     &                 TARGET_VALUE, ID, N, SORTED)
      INTEGER, intent(in):: VALUES_SIZE
      INTEGER, DIMENSION(0:VALUES_SIZE-1), intent(in):: VALUES
      INTEGER, intent(in):: TARGET_VALUE, ID, N
      INTEGER, DIMENSION(N+1), intent(inout):: SORTED
      INTEGER:: I,J
      INTEGER:: VAL, VAL_PREV
      INTEGER:: TMP
      LOGICAL:: ok
      VAL = VALUES(ID)
      SORTED(N+1) = ID
      ok = .false.
      I = N+1
      DO WHILE(.not. ok)
        IF (I.GT.1) THEN
          VAL_PREV = VALUES(SORTED(I-1))
          IF (VAL_PREV .GT. VAL) THEN
            TMP = SORTED(I)
            SORTED(I) = SORTED(I-1)
            SORTED(I-1) = TMP
            I = I-1
          ELSE
            ok = .true.
          ENDIF
        ELSE
          ok = .true.
        END IF
      END DO
#if defined(MTKO_ADVANCED_VERBOSE)
      write(*,*) 'MUMPS_MPITOOMP_INSERT_ID_SORTED :',
     &           SORTED(:)
#endif
      END SUBROUTINE MUMPS_MPITOOMP_INSERT_ID_SORTED
      SUBROUTINE MUMPS_MPITOOMP_STATUS_BY_HOST(NPROCS, HOST_NUMBER,
     &                PROC_COUNT_TAB, HEAD_PROC_LINKED_LIST, PROCS_MAP,
     &                MPI_PROC_HOST_TAB, HOST_HEAD_ID_TAB,
     &                HOST_HEAD_COUNT_TAB, HOST_PROC_COUNT_TAB,
     &                ERROR)
      INTEGER, intent(in):: NPROCS
      INTEGER, intent(in):: HOST_NUMBER
      INTEGER, DIMENSION(0:NPROCS-1), intent(in)::
     &                                PROC_COUNT_TAB
      INTEGER, DIMENSION(0:NPROCS-1), intent(in)::
     &                                HEAD_PROC_LINKED_LIST
      INTEGER, DIMENSION(:), POINTER:: PROCS_MAP
      INTEGER, DIMENSION(0:NPROCS-1), intent(out):: MPI_PROC_HOST_TAB
      INTEGER, DIMENSION(HOST_NUMBER), intent(out):: HOST_HEAD_ID_TAB
      INTEGER, DIMENSION(HOST_NUMBER), intent(out):: HOST_HEAD_COUNT_TAB
      INTEGER, DIMENSION(HOST_NUMBER), intent(out):: HOST_PROC_COUNT_TAB
      INTEGER:: ERROR
      INTEGER:: I, K, KINI
      INTEGER:: HOST_COUNTER
      ERROR = 0
      MPI_PROC_HOST_TAB(:) = -1
      HOST_HEAD_COUNT_TAB(:) = 0
      HOST_PROC_COUNT_TAB(:) = 0
      HOST_COUNTER = 1
      K = NPROCS - 1 
      DO WHILE (K .gt. 0 .AND. HOST_COUNTER .LE. HOST_NUMBER)
        DO WHILE (HEAD_PROC_LINKED_LIST(K) .EQ. -1
     &            .OR. MPI_PROC_HOST_TAB(K) .NE. -1)
          K = K - 1
        END DO
        KINI = K
        DO WHILE (K .NE. -1)
          IF (MPI_PROC_HOST_TAB(K).EQ.-1) THEN
            MPI_PROC_HOST_TAB(K) = HOST_COUNTER
            IF (PROCS_MAP(K).EQ.K) THEN
              HOST_HEAD_ID_TAB(HOST_COUNTER) = K
              HOST_HEAD_COUNT_TAB(HOST_COUNTER) =
     &            HOST_HEAD_COUNT_TAB(HOST_COUNTER) + 1
              HOST_PROC_COUNT_TAB(HOST_COUNTER) =
     &            HOST_PROC_COUNT_TAB(HOST_COUNTER)
     &            + PROC_COUNT_TAB(K)
            ENDIF
          ENDIF
          K = HEAD_PROC_LINKED_LIST(K)
        END DO
        HOST_COUNTER = HOST_COUNTER + 1
        K = KINI - 1
      END DO
        DO I = 0, NPROCS-1
          IF (PROCS_MAP(I).NE.I) THEN
            IF (MPI_PROC_HOST_TAB(I).NE.-1) THEN
              MPI_PROC_HOST_TAB(I) = MPI_PROC_HOST_TAB(PROCS_MAP(I))
#if defined(MTKO_ADVANCED_VERBOSE)
            ELSE
              write(*,*) 'STATUS: for ', I, ' host ',
     &             MPI_PROC_HOST_TAB(I), ' from map ',
     &             MPI_PROC_HOST_TAB(PROCS_MAP(I))
#endif
            ENDIF
          ENDIF
        ENDDO
      END SUBROUTINE MUMPS_MPITOOMP_STATUS_BY_HOST
      SUBROUTINE MUMPS_MPITOOMP_AGGREG_BY_NUMBER(NPROCS, MTKO_K,
     &                HEAD_COUNT, PROC_COUNT, PROC_COUNT_TAB,
     &                HEAD_SORTED_TAB, AGGREGATION_PART_MAP,
     &                FINAL_HEAD_COUNT, NB_PB, HEAD_PB_TAB)
      INTEGER, intent(in):: NPROCS, MTKO_K, PROC_COUNT, HEAD_COUNT
      INTEGER, DIMENSION(0:NPROCS-1), intent(inout):: PROC_COUNT_TAB
      INTEGER, DIMENSION(0:HEAD_COUNT-1), intent(inout)::
     &                                    HEAD_SORTED_TAB
      INTEGER, DIMENSION(0:PROC_COUNT-1), intent(inout)::
     &                                    AGGREGATION_PART_MAP
      INTEGER, intent(out):: FINAL_HEAD_COUNT
      INTEGER, intent(out):: NB_PB
      INTEGER, DIMENSION(0:HEAD_COUNT-1), intent(out):: HEAD_PB_TAB
      INTEGER:: I, FIRST_OK, J, START, K
      INTEGER:: VAL
      INTEGER:: CURR_HEAD, NB, MISSING
      INTEGER:: COMPLETION_HEAD
      LOGICAL:: found
      FINAL_HEAD_COUNT = HEAD_COUNT
      NB_PB = 0
      FIRST_OK = -1
      found = .false.
      I = 0
      DO WHILE (.not. found .AND. I .LT. HEAD_COUNT)
        IF (PROC_COUNT_TAB(HEAD_SORTED_TAB(I)) .LT. MTKO_K)
     &     THEN
          I = I + 1
        ELSE
          found = .true.
        ENDIF
      END DO
#if defined(MTKO_ADVANCED_VERBOSE)
      write(*,*) 'before I ', I, ', START ', START
#endif
      START = 0
      IF (found) THEN
        I = I - 1
        FIRST_OK = I
      ELSE
        I = HEAD_COUNT-1
      ENDIF
      DO WHILE (I .GE. START)
#if defined(MTKO_ADVANCED_VERBOSE)
        write(*,*) 'I ',I, ', START ',START
#endif
        CURR_HEAD = HEAD_SORTED_TAB(I)
        NB = PROC_COUNT_TAB(CURR_HEAD)
#if defined(MTKO_ADVANCED_VERBOSE)
        write(*,*) 'CURR_HEAD ',CURR_HEAD, ', NB ',NB
#endif
        MISSING = MTKO_K - NB
        J = START
        found = .false.
        DO WHILE (.not. found .AND. J .LT. I)
#if defined(MTKO_ADVANCED_VERBOSE)
        write(*,*) 'searching complement, J', J
#endif
          VAL = PROC_COUNT_TAB(HEAD_SORTED_TAB(J))
          IF (VAL .LT. MISSING) THEN
            J = J + 1
          ELSE IF (VAL .EQ. MISSING) THEN
            found = .true.
          ELSE
            J = I
          ENDIF
        END DO
#if defined(MTKO_ADVANCED_VERBOSE)
        write(*,*) 'found ', found, 'NB ', NB, 'MISSING ',
     &            MISSING, 'J ', J
#endif
        IF (found) THEN
          COMPLETION_HEAD = HEAD_SORTED_TAB(J)
#if defined(MTKO_ADVANCED_VERBOSE)
          write(*,*) 'COMPLETION_HEAD = ', COMPLETION_HEAD
#endif
          DO K = 0, PROC_COUNT-1
            IF (AGGREGATION_PART_MAP(K) .EQ. CURR_HEAD)
     &         THEN
              AGGREGATION_PART_MAP(K) = COMPLETION_HEAD
            ENDIF
          END DO
          PROC_COUNT_TAB(COMPLETION_HEAD) = MTKO_K
          PROC_COUNT_TAB(CURR_HEAD) = 0
          HEAD_SORTED_TAB(I) = COMPLETION_HEAD
          HEAD_SORTED_TAB(J) = CURR_HEAD
          DO K = J-1, START, -1
            HEAD_SORTED_TAB(K+1) = HEAD_SORTED_TAB(K)
          END DO
          FINAL_HEAD_COUNT = FINAL_HEAD_COUNT - 1
          FIRST_OK = I
          I = I - 1
        ELSE
          HEAD_PB_TAB(NB_PB) = HEAD_SORTED_TAB(I)
          NB_PB = NB_PB + 1
          DO K = I-1, START, -1
            HEAD_SORTED_TAB(K+1) = HEAD_SORTED_TAB(K)
          END DO
        ENDIF
        START = START + 1
      END DO
#if defined(MTKO_ADVANCED_VERBOSE)
      write(*,*) 'End of aggregate by number: FIRST_OK ', FIRST_OK,
     &           ', NB_PB ', NB_PB
#endif
      END SUBROUTINE MUMPS_MPITOOMP_AGGREG_BY_NUMBER
      SUBROUTINE MUMPS_MPITOOMP_AGGREGATE_PASS2(NPROCS, MTKO_K, NB_PB,
     &                PROC_COUNT, PROC_COUNT_TAB, HEAD_PB_TAB,
     &                AGGREGATION_PART_MAP, FINAL_NB_PB)
      INTEGER, intent(in):: NPROCS, MTKO_K, NB_PB, PROC_COUNT
      INTEGER, DIMENSION(0:NPROCS-1), intent(inout):: PROC_COUNT_TAB
      INTEGER, DIMENSION(0:NB_PB-1), intent(inout):: HEAD_PB_TAB
      INTEGER, DIMENSION(0:PROC_COUNT-1), intent(inout)::
     &                                    AGGREGATION_PART_MAP
      INTEGER, intent(out):: FINAL_NB_PB
      INTEGER:: I, J, K, START
      INTEGER:: CURR_HEAD, CURR_COUNT, COMPLETION_HEAD, COMPLETION_COUNT
      LOGICAL:: found
      FINAL_NB_PB = NB_PB
#if defined(MTKO_ADVANCED_VERBOSE)
      write(*,*) 'HEAD_PB_TAB : ', HEAD_PB_TAB(:)
#endif
      START = 0
      I = FINAL_NB_PB-1
      DO WHILE (I .GT. START)
        CURR_HEAD = HEAD_PB_TAB(I)
        CURR_COUNT = PROC_COUNT_TAB(CURR_HEAD)
        COMPLETION_HEAD = HEAD_PB_TAB(START)
        COMPLETION_COUNT = PROC_COUNT_TAB(START)
        J = 0
        K = 0
        DO WHILE (J .LT. (MTKO_K-CURR_COUNT))
          found = .false.
          DO WHILE(.not. found)
            IF (AGGREGATION_PART_MAP(K).EQ. COMPLETION_HEAD)
     &         THEN
              AGGREGATION_PART_MAP(K) = CURR_HEAD
              PROC_COUNT_TAB(CURR_HEAD) =
     &            PROC_COUNT_TAB(CURR_HEAD) + 1
              PROC_COUNT_TAB(COMPLETION_HEAD) =
     &            PROC_COUNT_TAB(COMPLETION_HEAD) - 1
              IF (PROC_COUNT_TAB(COMPLETION_HEAD).EQ.0) THEN
                START = START + 1
                COMPLETION_HEAD = HEAD_PB_TAB(FINAL_NB_PB-1)
                COMPLETION_COUNT =
     &             PROC_COUNT_TAB(COMPLETION_HEAD)
                K = 0
              ENDIF
              found = .true.
            ELSE
              K = K + 1
            ENDIF
          END DO
          J = J + 1
        END DO
        I = I - 1
      END DO
      FINAL_NB_PB = I + 1 - START
#if defined(MTKO_ADVANCED_VERBOSE)
      write(*,*) 'FINAL_NB_PB ', FINAL_NB_PB
#endif
      END SUBROUTINE MUMPS_MPITOOMP_AGGREGATE_PASS2
#endif
      SUBROUTINE MUMPS_MPITOOMP_FILLIN_PROCS_MAP(COMM, MTKO_K, NPROCS,
     &                PROCS_MAP, ERROR)
      INCLUDE 'mpif.h'
      INTEGER, intent(in):: COMM
      INTEGER, intent(in):: MTKO_K, NPROCS
      INTEGER, DIMENSION(:), POINTER:: PROCS_MAP
      INTEGER, intent(out):: ERROR
      INTEGER:: IERR, I,J
      INTEGER:: MYID
      LOGICAL:: found
      CHARACTER(len=MPI_MAX_PROCESSOR_NAME):: PROC_NAME
      INTEGER:: PROC_NAME_LENGTH
      CHARACTER(len=MPI_MAX_PROCESSOR_NAME), DIMENSION(0:NPROCS-1)::
     &          PROC_NAME_TAB
      INTEGER, DIMENSION(0:NPROCS-1):: PROC_NAME_LENGTH_TAB
      INTEGER, DIMENSION(0:NPROCS-1):: PROC_COUNT_TAB
#if defined(MTKO_ADVANCED)
      INTEGER:: cpu_id, node_id
      INTEGER:: K, HEAD_NUMBER, HOST_NUMBER, S
      INTEGER, DIMENSION(0:NPROCS-1):: PROC_NODE_TAB
      INTEGER, DIMENSION(0:NPROCS-1):: HEAD_PROC_LINKED_LIST
      INTEGER, DIMENSION(0:NPROCS-1):: MPI_PROC_HOST_TAB
      INTEGER:: MYHEAD, THISHOST
      INTEGER:: MYHOST, MYHOST_HEAD, MYHOST_HEAD_COUNT
      INTEGER:: MYHOST_PROC_COUNT
      INTEGER, DIMENSION(:), ALLOCATABLE:: HOST_HEAD_ID_TAB
      INTEGER, DIMENSION(:), ALLOCATABLE:: HOST_HEAD_COUNT_TAB
      INTEGER, DIMENSION(:), ALLOCATABLE:: HOST_PROC_COUNT_TAB
      INTEGER, DIMENSION(:), ALLOCATABLE:: HOST_TAB
      INTEGER, DIMENSION(:), ALLOCATABLE:: MYHOST_PROC_IDS_TAB
      INTEGER, DIMENSION(:), ALLOCATABLE:: AGGREGATION_PART_MAP
      INTEGER, DIMENSION(:), ALLOCATABLE:: MYHOST_HEAD_SORTED
      INTEGER, DIMENSION(:), ALLOCATABLE:: MYHOST_HEAD_PB
      INTEGER:: FIRST_OK, START, NB_PB, MISSING, VAL, NB
      INTEGER:: FINAL_NB_PB, FINAL_HEAD_COUNT, MPIID
      INTEGER:: allocok
      INTEGER:: MPI_STAT(MPI_STATUS_SIZE)
      INTEGER:: HOST_OF_MASTER, HEAD_OF_MASTER, COUNTER
      INTEGER, DIMENSION(:), ALLOCATABLE:: AGGREGATED_MAPS, HOST_START,
     &                                     HOST_I
#endif
      call MPI_COMM_RANK(COMM, MYID, IERR)
      ERROR = 0
      call MPI_GET_PROCESSOR_NAME(PROC_NAME, PROC_NAME_LENGTH, IERR)
#if defined(MTKO_ADVANCED)
      CALL MUMPS_GET_CPU_ID(cpu_id)
      CALL MUMPS_GET_NUMA_NODE_OF_CPU(cpu_id, node_id)
#if defined(MTKO_ADVANCED_VERBOSE)
      write(*,*) '[', MYID,  '] runs on proc ', cpu_id, ' node ',
     &           node_id
#endif
#endif
#if defined(MTKO_ADVANCED)
#endif
      DO I = 0, NPROCS-1
        IF (MYID .eq. I) THEN
          PROC_NAME_LENGTH_TAB(I) = PROC_NAME_LENGTH
        ENDIF
        CALL MPI_BCAST(PROC_NAME_LENGTH_TAB(I), 1, MPI_INTEGER, I, COMM,
     &          IERR)
        IF (MYID .eq. I) THEN
          PROC_NAME_TAB(I) = PROC_NAME
        ENDIF
        CALL MPI_BCAST(PROC_NAME_TAB(I), PROC_NAME_LENGTH_TAB(I),
     &          MPI_CHARACTER, I, COMM, IERR)
#if defined(MTKO_ADVANCED)
        IF (MYID .eq. I) THEN
          PROC_NODE_TAB(I) = node_id
        ENDIF
        CALL MPI_BCAST(PROC_NODE_TAB(I), 1, MPI_INTEGER, I, COMM, IERR)
#endif
      END DO
#if defined(MTKO_ADVANCED_VERBOSE)
      IF (MYID .EQ. 0) THEN
        write(*,*) 'PROC_NODE_TAB : ', PROC_NODE_TAB(:)
      ENDIF
#endif
#if defined(MTKO_ADVANCED)
#endif
      DO I=0, NPROCS-1
        PROC_COUNT_TAB(I) = 0
      END DO
#if defined(MTKO_ADVANCED)
      DO I = 0, NPROCS-1
        HEAD_PROC_LINKED_LIST(I) = -1
      END DO
#endif
      PROCS_MAP(0) = 0
      PROC_COUNT_TAB(0) = 1
#if defined(MTKO_ADVANCED)
      HEAD_NUMBER = 1
      HOST_NUMBER = 1
#endif
      DO I = 1, NPROCS-1
        J = I-1
        found = .false.
        DO WHILE (.not. found .and. J .ge. 0)
          IF (STRINGS_EQUALITY(PROC_NAME_LENGTH_TAB(J), PROC_NAME_TAB(J)
     &        , PROC_NAME_LENGTH_TAB(I), PROC_NAME_TAB(I))) THEN
            found = .true.
          ELSE
            J = J - 1
          ENDIF
        END DO
#if !defined(MTKO_ADVANCED)
        IF (found) THEN
          IF (PROC_COUNT_TAB(PROCS_MAP(J)) .lt. MTKO_K) THEN
            PROCS_MAP(I) = PROCS_MAP(J)
            PROC_COUNT_TAB(PROCS_MAP(J)) = PROC_COUNT_TAB(PROCS_MAP(J))
     &                                   + 1
          ELSE
            PROCS_MAP(I) = I
            PROC_COUNT_TAB(I) = 1
          ENDIF
        ELSE
          PROCS_MAP(I) = I
          PROC_COUNT_TAB(I) = 1
        ENDIF
#else
        IF (.not. found) THEN
#if defined(MTKO_ADVANCED_VERBOSE)
          IF (MYID .EQ. 0) THEN
            write(*,*) 'name not found for proc ', I,
     &          '(', PROC_NAME_TAB(I)(1:PROC_NAME_LENGTH_TAB(I)), ')'
          ENDIF
#endif
          PROCS_MAP(I) = I
          PROC_COUNT_TAB(I) = 1
          HEAD_NUMBER = HEAD_NUMBER + 1
          HOST_NUMBER = HOST_NUMBER + 1
        ELSE
          found = .false.
          K = J
          DO WHILE (.not. found .and. K .ne. -1)
#if defined(MTKO_ADVANCED_VERBOSE)
            IF (MYID .EQ. 0) THEN
              write (*,*) 'PROC_NODE_TAB(', K, ') = ', PROC_NODE_TAB(K),
     &                ', PROC_COUNT_TAB = ', PROC_COUNT_TAB(K)
            ENDIF
#endif
            IF (PROC_NODE_TAB(K) .eq. PROC_NODE_TAB(I)
     &          .and. PROC_COUNT_TAB(PROCS_MAP(K)) .lt. MTKO_K) THEN
              found = .true.
            ELSE
              K = HEAD_PROC_LINKED_LIST(K)
            ENDIF
          END DO
#if defined(MTKO_ADVANCED_VERBOSE)
          IF (MYID .EQ. 0) THEN
            write(*,*) 'for proc ', I, ' K is ', K
          ENDIF
#endif
          IF (K .eq. -1) THEN
            PROCS_MAP(I) = I
            PROC_COUNT_TAB(I) = 1
            HEAD_NUMBER = HEAD_NUMBER + 1
          ELSE
            PROCS_MAP(I) = PROCS_MAP(K)
            PROC_COUNT_TAB(PROCS_MAP(K)) = PROC_COUNT_TAB(PROCS_MAP(K))
     &                                     + 1
          ENDIF
          HEAD_PROC_LINKED_LIST(I) = J
#if defined(MTKO_ADVANCED_VERBOSE)
          IF (MYID.EQ.0) THEN
            write(*,*) 'PROCS_MAP : ', PROCS_MAP(0:NPROCS-1)
            write(*,*) 'HEAD_PROC_LINKED_LIST : ',
     &                 HEAD_PROC_LINKED_LIST(0:NPROCS-1)
          ENDIF
#endif
        ENDIF
#endif
      END DO
#if defined(MTKO_ADVANCED)
#if defined(MTKO_ADVANCED_VERBOSE)
      IF (MYID .EQ. 0) THEN
        write(*,*) 'HEAD_NUMBER ', HEAD_NUMBER, ', NPROCS ', NPROCS,
     &           ', MTKO_K ', MTKO_K, ', HOST_NUMBER ', HOST_NUMBER
        write(*,*) 'PROC_COUNT_TAB ', PROC_COUNT_TAB(:)
        write(*,*) 'HEAD_PROC_LINKED_LIST ', HEAD_PROC_LINKED_LIST(:)
        write(*,*) 'PROCS_MAP ', PROCS_MAP(:)
        write(*,*) 'LBOUND PROCS_MAP ', LBOUND(PROCS_MAP), ' UBOUND ',
     &              UBOUND(PROCS_MAP)
      ENDIF
#endif
      IF (HEAD_NUMBER .gt. NPROCS/MTKO_K) THEN
        allocate(HOST_HEAD_ID_TAB(1:HOST_NUMBER),
     &           HOST_HEAD_COUNT_TAB(1:HOST_NUMBER),
     &           HOST_PROC_COUNT_TAB(1:HOST_NUMBER), STAT=allocok)
        IF (allocok .ne. 0) THEN
          ERROR = -1
          GOTO 100
        ENDIF
        call MUMPS_MPITOOMP_STATUS_BY_HOST(NPROCS, HOST_NUMBER,
     &            PROC_COUNT_TAB, HEAD_PROC_LINKED_LIST, PROCS_MAP,
     &            MPI_PROC_HOST_TAB, HOST_HEAD_ID_TAB,
     &            HOST_HEAD_COUNT_TAB, HOST_PROC_COUNT_TAB, ERROR)
        IF (ERROR .NE. 0) THEN
          ERROR = -2
          GOTO 100
        ENDIF
#if defined(MTKO_ADVANCED_VERBOSE)
        IF (MYID.EQ.0) THEN
          write(*,*) MYID,'HOST_HEAD_COUNT_TAB ', HOST_HEAD_COUNT_TAB(:)
          write(*,*) MYID,'HOST_HEAD_ID_TAB ', HOST_HEAD_ID_TAB(:)
          write(*,*) MYID,'HOST_PROC_COUNT_TAB ', HOST_PROC_COUNT_TAB(:)
        ENDIF
#endif
        MYHOST = MPI_PROC_HOST_TAB(PROCS_MAP(MYID))
        MYHOST_HEAD = HOST_HEAD_ID_TAB(MYHOST)
        MYHOST_HEAD_COUNT = HOST_HEAD_COUNT_TAB(MYHOST)
        MYHOST_PROC_COUNT = HOST_PROC_COUNT_TAB(MYHOST)
#if defined(MTKO_ADVANCED_VERBOSE)
        write(*,*) '[',MYID, '] MYHOST ', MYHOST, ', MYHOST_HEAD ',
     &    MYHOST_HEAD, ', MYHOST_HEAD_COUNT', MYHOST_HEAD_COUNT,
     &    ', MYHOST_PROC_COUNT ', MYHOST_PROC_COUNT
#endif
        IF (modulo(MYHOST_PROC_COUNT,MTKO_K).NE.0) THEN
          ERROR = -2
          GOTO 100
        ENDIF
        HEAD_OF_MASTER = MPITOOMP_MASTER
        IF (MYID .EQ. MYHOST_HEAD) THEN
          allocate(AGGREGATION_PART_MAP(0:MYHOST_PROC_COUNT-1),
     &             STAT=allocok)
          IF (allocok .ne. 0) THEN
            ERROR = -1
            GOTO 100
          ENDIF
          IF (MYHOST_PROC_COUNT .EQ. MTKO_K) THEN
#if defined(MTKO_ADVANCED_VERBOSE)
            write(*,*) 'all MPI processes should be mapped to the',
     &                 ' first group head of this node: ', MYHOST_HEAD
#endif
            AGGREGATION_PART_MAP(:) = MYHOST_HEAD
            PROC_COUNT_TAB(MYHOST_HEAD) = MYHOST_PROC_COUNT
          ELSE
            allocate(MYHOST_PROC_IDS_TAB(0:MYHOST_PROC_COUNT-1),
     &           STAT=allocok)
            IF (allocok .ne. 0) THEN
              ERROR = -1
              GOTO 80
            ENDIF
            J = 0 
            DO I = 0, NPROCS-1
              THISHOST = MPI_PROC_HOST_TAB(PROCS_MAP(I))
              IF (THISHOST .EQ. MYHOST) THEN
                MYHOST_PROC_IDS_TAB(J) = I
                J = J + 1
              ENDIF
            END DO
            DO I = 0, MYHOST_PROC_COUNT-1
              AGGREGATION_PART_MAP(I)= PROCS_MAP(MYHOST_PROC_IDS_TAB(I))
            END DO
#if defined(MTKO_ADVANCED_VERBOSE)
            write(*,*) 'AGGREGATION_PART_MAP ', AGGREGATION_PART_MAP(:)
#endif
            IF (MYHOST_PROC_COUNT/MTKO_K .EQ. MYHOST_HEAD_COUNT) THEN
#if defined(MTKO_ADVANCED_VERBOSE)
              write(*,*) 'nothing to do, this is correct for this node'
              write(*,*) 'the AGGREGATION_PART_MAP is going to be send'
     &                   , ' as it is'
#endif
            ELSE
              allocate(MYHOST_HEAD_SORTED(0:MYHOST_HEAD_COUNT-1),
     &           MYHOST_HEAD_PB(0:MYHOST_HEAD_COUNT-1),
     &           STAT=allocok)
              IF (allocok .ne. 0) THEN
                ERROR = -1
                GOTO 40
              ENDIF
#if defined(MTKO_ADVANCED_VERBOSE)
              write(*,*) 'before sort, values: ', PROC_COUNT_TAB
#endif
              K = 0 
              DO I = 0, MYHOST_PROC_COUNT-1
                MPIID = MYHOST_PROC_IDS_TAB(I)
                IF (PROCS_MAP(MPIID).EQ.MPIID) THEN
                  call MUMPS_MPITOOMP_INSERT_ID_SORTED(NPROCS,
     &                      PROC_COUNT_TAB, MTKO_K, MPIID, K,
     &                      MYHOST_HEAD_SORTED)
                  K = K + 1
                ENDIF
              END DO
              call MUMPS_MPITOOMP_AGGREG_BY_NUMBER(NPROCS, MTKO_K,
     &                  MYHOST_HEAD_COUNT, MYHOST_PROC_COUNT,
     &                  PROC_COUNT_TAB,  MYHOST_HEAD_SORTED,
     &                  AGGREGATION_PART_MAP,
     &                  FINAL_HEAD_COUNT, NB_PB, MYHOST_HEAD_PB)
#if defined(MTKO_ADVANCED_VERBOSE)
              write(*,*) 'FINAL_HEAD_COUNT ', FINAL_HEAD_COUNT
              write(*,*) 'NB_PB', NB_PB, ', FIRST_OK ', FIRST_OK
              write(*,*) 'AGGREGATION_PART_MAP ',
     &                    AGGREGATION_PART_MAP(:)
#endif
              IF (FINAL_HEAD_COUNT .NE. MYHOST_PROC_COUNT/MTKO_K) THEN
#if defined(MTKO_ADVANCED_VERBOSE)
                write(*,*) 'not yet ok MYHOST_HEAD_COUNT ',
     &                     FINAL_HEAD_COUNT
#endif
                call MUMPS_MPITOOMP_AGGREGATE_PASS2(NPROCS, MTKO_K,
     &                    NB_PB, MYHOST_PROC_COUNT, PROC_COUNT_TAB,
     &                    MYHOST_HEAD_PB, AGGREGATION_PART_MAP,
     &                    FINAL_NB_PB)
              ENDIF 
40            IF (allocated(MYHOST_HEAD_SORTED)) THEN
                deallocate(MYHOST_HEAD_SORTED)
              ENDIF
              IF (allocated(MYHOST_HEAD_PB)) THEN
                deallocate(MYHOST_HEAD_PB)
              ENDIF
            ENDIF 
            IF (allocated(MYHOST_PROC_IDS_TAB)) THEN
              deallocate(MYHOST_PROC_IDS_TAB)
            ENDIF
            IF (ERROR .NE. 0) THEN
              GOTO 80
            ENDIF
          ENDIF 
          IF (HOST_NUMBER.gt.1) THEN
            HOST_OF_MASTER = MPI_PROC_HOST_TAB(
     &                          PROCS_MAP(MPITOOMP_MASTER))
            HEAD_OF_MASTER = HOST_HEAD_ID_TAB(HOST_OF_MASTER)
            IF (MYHOST .EQ. HOST_OF_MASTER) THEN
#if defined(MTKO_ADVANCED_VERBOSE)
              write(*,*) 'I am MPI of rank ', MYID, ', the group head ',
     &                   ' of the host of the master '
#endif
              allocate(AGGREGATED_MAPS(0:NPROCS-1),
     &                 HOST_START(HOST_NUMBER),
     &                 HOST_I(HOST_NUMBER), STAT=allocok)
              IF (allocok .ne. 0) THEN
                ERROR = -2
                GOTO 60
              ENDIF
              AGGREGATED_MAPS(0:MYHOST_PROC_COUNT-1) =
     &                        AGGREGATION_PART_MAP(:)
              HOST_START(MYHOST) = 0
              START = MYHOST_PROC_COUNT
              DO K = 1, HOST_NUMBER-1
                call MPI_RECV(AGGREGATED_MAPS(START), NPROCS-START,
     &                  MPI_INTEGER, MPI_ANY_SOURCE, TAG_AGGREGATED_MAP,
     &                  COMM, MPI_STAT, IERR)
                THISHOST = MPI_PROC_HOST_TAB(MPI_STAT(MPI_SOURCE))
                call MPI_GET_COUNT(MPI_STAT, MPI_INTEGER,
     &                  COUNTER, IERR)
                HOST_START(THISHOST) = START
                START = START + COUNTER
#if defined(MTKO_ADVANCED_VERBOSE)
                write(*,*) 'received #', COUNTER,
     &                     ' elements from ', MPI_STAT(MPI_SOURCE),
     &                     ' of host ', THISHOST, ' starting at ',
     &                     HOST_START(THISHOST)
#endif
              END DO
              HOST_I(:) = 0
              DO I = 0, NPROCS-1
                THISHOST = MPI_PROC_HOST_TAB(I)
                PROCS_MAP(I) = AGGREGATED_MAPS( HOST_START(THISHOST)
     &                         + HOST_I(THISHOST) )
                HOST_I(THISHOST) = HOST_I(THISHOST) + 1
              END DO
60            IF (allocated(AGGREGATED_MAPS)) THEN
                deallocate(AGGREGATED_MAPS)
              ENDIF
              IF (allocated(HOST_START)) THEN
                deallocate(HOST_START)
              ENDIF
              IF (allocated(HOST_I)) THEN
                deallocate(HOST_I)
              ENDIF
              IF (ERROR .ne. 0) THEN
                GOTO 80
              ENDIF
            ELSE
                call MPI_SEND(AGGREGATION_PART_MAP(0), MYHOST_PROC_COUNT
     &                  , MPI_INTEGER, HEAD_OF_MASTER,
     &                  TAG_AGGREGATED_MAP, COMM, IERR)
            ENDIF
          ELSE
            DO I = 0, NPROCS-1
              PROCS_MAP(I) = AGGREGATION_PART_MAP(I)
            ENDDO
          ENDIF
80        IF (allocated(AGGREGATION_PART_MAP)) THEN
            deallocate(AGGREGATION_PART_MAP)
          ENDIF
          IF (ERROR .NE. 0) THEN
            GOTO 100
          ENDIF
        ENDIF 
        call MPI_BCAST(PROCS_MAP(0), NPROCS, MPI_INTEGER,
     &                 HEAD_OF_MASTER, COMM, IERR)
100     IF (allocated(HOST_HEAD_ID_TAB)) THEN
          deallocate(HOST_HEAD_ID_TAB)
        ENDIF
        IF (allocated(HOST_HEAD_COUNT_TAB)) THEN
          deallocate(HOST_HEAD_COUNT_TAB)
        ENDIF
        IF (allocated(HOST_PROC_COUNT_TAB)) THEN
          deallocate(HOST_PROC_COUNT_TAB)
        ENDIF
      ENDIF 
#endif
#if defined(MTKO_ADVANCED_VERBOSE)
      IF (MYHOST_HEAD .EQ. MYID) THEN
        total = 0
        DO I = 0, NPROCS - 1
          IF (MPI_PROC_HOST_TAB(I) .EQ. MYHOST) THEN
            IF (PROC_COUNT_TAB(I).ne.0 .and.PROC_COUNT_TAB(I).ne.MTKO_K)
     &      THEN
              write (*,*) "problem for MPI proc ", I, " with ",
     &                    PROC_COUNT_TAB(I)
            ENDIF
            total = total + PROC_COUNT_TAB(I)
          ENDIF
        END DO
        write(*,*) '[', MYID, '] total in host ', MYHOST, ' is ',
     &             total
      ENDIF
      IF (PROCS_MAP(MYID) .EQ. MYID) THEN
        write (*,*) '[', MYID, '] PROCS_MAP ', PROCS_MAP(:)
      ENDIF
#else
#endif
      END SUBROUTINE MUMPS_MPITOOMP_FILLIN_PROCS_MAP
      SUBROUTINE MUMPS_MPITOOMP_PRE_WORK_MAPS(MTKO_K, NPROCS,
     &                NCORES,
     &                ERROR)
#if defined(_OPENMP)
      USE OMP_LIB
#endif
#if !defined(MTKO_NO_MASK_FUSION)
#else
#endif
#if !defined(MTKO_NO_MASK_FUSION)
#else
#endif
      INTEGER, intent(in):: MTKO_K, NPROCS
      INTEGER, intent(out):: NCORES
      INTEGER, intent(out):: ERROR
      INTEGER:: NOMP
      ERROR = 0
      NOMP = 1
!$    NOMP = OMP_GET_MAX_THREADS()
      NCORES = NOMP*MTKO_K
#if !defined(MTKO_NO_MASK_FUSION)
      NCORES = MTKO_K
#endif
      END SUBROUTINE MUMPS_MPITOOMP_PRE_WORK_MAPS
      SUBROUTINE MUMPS_MPITOOMP_GET_WORK_MAPS(COMM, MTKO_K, NPROCS,
     &                NCORES, PROCS_MAP, CORES_LIST, ERROR)
#if defined(_OPENMP)
      USE OMP_LIB
#endif
      include 'mpif.h'
#if !defined(MTKO_NO_MASK_FUSION)
#else
#endif
#if !defined(MTKO_NO_MASK_FUSION)
#else
#endif
      INTEGER, intent(in):: COMM
      INTEGER, intent(in):: MTKO_K, NPROCS, NCORES
      INTEGER, DIMENSION(:), POINTER:: PROCS_MAP
      INTEGER, DIMENSION(0:NCORES-1), intent(out):: CORES_LIST
      INTEGER, intent(out):: ERROR
      INTEGER:: MPI_STAT(MPI_STATUS_SIZE)
      INTEGER:: MYID, IERR
      INTEGER:: nb_msg
#if !defined(MTKO_NO_MASK_FUSION)
      INTEGER:: pid
#else
      INTEGER:: cpu_id,thread_num
#endif
      call MPI_COMM_RANK(COMM, MYID, IERR)
      ERROR = 0
#if defined(MTKO_NO_MASK_FUSION)
      CALL MUMPS_GET_CPU_ID(cpu_id)
      IF (NCORES .gt. MTKO_K) THEN
        thread_num = 0
!$OMP PARALLEL private(cpu_id,thread_num)
        CALL MUMPS_GET_CPU_ID(cpu_id)
!$      thread_num = OMP_GET_THREAD_NUM()
        CORES_LIST(thread_num) = cpu_id
!$OMP END PARALLEL
      ELSE
        CORES_LIST(0) = cpu_id
      ENDIF
#else
      CALL MUMPS_GET_PID(pid)
      CORES_LIST(0) = pid;
#endif
      IF (PROCS_MAP(MYID) .eq. MYID) THEN
        nb_msg = 2
        DO WHILE (nb_msg .le. MTKO_K)
          call MPI_RECV(CORES_LIST((nb_msg-1)*(NCORES/MTKO_K)),
     &            NCORES/MTKO_K, MPI_INTEGER, MPI_ANY_SOURCE,
     &            TAG_CORES_LIST_PART, COMM, MPI_STAT, IERR)
          nb_msg = nb_msg + 1
        END DO
      ELSE
        call MPI_SEND(CORES_LIST(0), NCORES/MTKO_K, MPI_INTEGER,
     &          PROCS_MAP(MYID), TAG_CORES_LIST_PART, COMM, IERR)
      END IF
      END SUBROUTINE MUMPS_MPITOOMP_GET_WORK_MAPS
      SUBROUTINE MUMPS_MPITOOMP_SET_MASK(NCORES, NOMP, MASK, CHECK,
     &                ERROR)
#if defined(_OPENMP)
      USE OMP_LIB
#endif
      INTEGER, intent(in):: NCORES
      INTEGER, intent(in):: NOMP
      INTEGER, DIMENSION(0:NCORES-1), intent(in):: MASK
      LOGICAL, intent(in):: CHECK
      INTEGER, intent(out):: ERROR
#if defined(MTKO_NO_MASK_FUSION)      
      INTEGER:: IS_INCLUDED
#endif      
      LOGICAL:: SET
      ERROR = 0
      IF (CHECK) THEN
        SET = .true.
#if defined(MTKO_NO_MASK_FUSION)
        CALL MUMPS_CHECK_CPU_SET_INCLUSION(NCORES, MASK, IS_INCLUDED)
        IF (IS_INCLUDED .LT. 0) THEN
        ELSE
          IF (IS_INCLUDED .EQ. 1) THEN
            SET = .false.
          ELSE
            SET = .true.
          ENDIF
        ENDIF
#endif
        IF (SET) THEN
          call MUMPS_ASSIGN_THREADS_FROM_LIST(NCORES,
     &              NOMP, MASK, 0, ERROR)
          IF (ERROR .LT. 0) THEN
#if defined(WORKAROUNDINTELILP64OPENMPLIMITATION)
!$          CALL omp_set_num_threads(int(NOMP,4))
#else
!$          CALL omp_set_num_threads(NOMP)
#endif
          ENDIF
        ELSE
#if defined(WORKAROUNDINTELILP64OPENMPLIMITATION)
!$        CALL omp_set_num_threads(int(NCORES,4))
#else
!$        CALL omp_set_num_threads(NCORES)
#endif
        ENDIF
      ELSE
        call MUMPS_ASSIGN_THREADS_FROM_LIST(NCORES,
     &            NOMP, MASK, 1, ERROR)
        IF (ERROR .LT. 0) THEN
#if defined(WORKAROUNDINTELILP64OPENMPLIMITATION)
!$        CALL omp_set_num_threads(int(NOMP,4))
#else
!$        CALL omp_set_num_threads(NOMP)
#endif
        ENDIF
      ENDIF
      END SUBROUTINE MUMPS_MPITOOMP_SET_MASK
      FUNCTION STRINGS_EQUALITY(L1, NAME1, L2,NAME2)
      INTEGER, intent(in):: L1
      CHARACTER(len=L1), intent(in):: NAME1
      INTEGER, intent(in):: L2
      CHARACTER(len=L2), intent(in):: NAME2
      LOGICAL:: STRINGS_EQUALITY
      INTEGER:: I
      LOGICAL:: ok
      STRINGS_EQUALITY = .false.
      IF (L1.eq.L2) THEN
        I = 1
        ok = .true.
        DO WHILE (ok .and. I .le. L1)
          IF (NAME1(I:I) .ne. NAME2(I:I)) THEN
            ok = .false.
          ENDIF
          I = I + 1
        END DO
        IF (ok) THEN
          STRINGS_EQUALITY = .true.
        ENDIF
      ENDIF
      END FUNCTION STRINGS_EQUALITY
      SUBROUTINE MUMPS_MPITOOMP_SLEEPSYNC( KEEP, PROCS_MAP,
     &           MYID_NODES, COMM_NODES, NSLAVES, I_AM_SLAVE )
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER, INTENT(IN) :: KEEP(500), MYID_NODES, COMM_NODES, NSLAVES
      LOGICAL, INTENT(IN) :: I_AM_SLAVE
      INTEGER, POINTER :: PROCS_MAP(:)
      INTEGER :: MPI_STAT(MPI_STATUS_SIZE), IDEST, IDUMMY(1), IERR
      LOGICAL :: FLAG
      IDUMMY(1)=1
      IF (I_AM_SLAVE .AND. KEEP(381) .GT. 1) THEN
        IF ( PROCS_MAP(MYID_NODES) .EQ. MYID_NODES ) THEN
          DO IDEST = 0, NSLAVES-1
            IF ( MYID_NODES .NE. IDEST .AND.
     &         MYID_NODES .EQ. PROCS_MAP(IDEST) ) THEN
              CALL MPI_SEND(IDUMMY, 1, MPI_INTEGER, IDEST,
     &           TAG_WAKEUP, COMM_NODES, IERR )
            ENDIF
          ENDDO
        ELSE
          DO WHILE (.TRUE. )
            CALL MPI_IPROBE(PROCS_MAP(MYID_NODES), TAG_WAKEUP,
     &      COMM_NODES, FLAG, MPI_STAT, IERR )
            IF ( FLAG ) THEN
              CALL MPI_RECV(IDUMMY, 1, MPI_INTEGER,
     &        PROCS_MAP(MYID_NODES), TAG_WAKEUP, COMM_NODES,
     &        MPI_STAT, IERR )
              EXIT
            ENDIF
            CALL MUMPS_USLEEP(1000)
          ENDDO
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_MPITOOMP_SLEEPSYNC
      END MODULE MUMPS_MPITOOMP_M
#else
      MODULE MUMPS_MPITOOMP_M
      CONTAINS
      SUBROUTINE MUMPS_MPITOOMP_M_RETURN()
        RETURN
      END SUBROUTINE MUMPS_MPITOOMP_M_RETURN
      END MODULE MUMPS_MPITOOMP_M
#endif
