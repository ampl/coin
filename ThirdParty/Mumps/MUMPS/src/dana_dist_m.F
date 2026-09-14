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
      SUBROUTINE DMUMPS_ANA_COMPUTE_ESTIMATES ( id, idintr, NBLK )
      USE DMUMPS_STRUC_DEF, ONLY: DMUMPS_STRUC
      USE DMUMPS_INTR_TYPES, ONLY: DMUMPS_INTR_STRUC
      USE MUMPS_ANA_OMP_M, ONLY: MUMPS_ANA_L0_OMP
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER   :: IERR, MASTER
      PARAMETER( MASTER = 0 )
      TYPE(DMUMPS_STRUC), TARGET  :: id
      TYPE(DMUMPS_INTR_STRUC) :: idintr
      INTEGER, INTENT(IN)     :: NBLK
      INTEGER, DIMENSION(:), POINTER :: KEEP, INFO, INFOG
      INTEGER :: allocok
      INTEGER(8), DIMENSION(:), POINTER :: KEEP8
      DOUBLE PRECISION, DIMENSION(:), POINTER :: RINFO
      DOUBLE PRECISION, DIMENSION(:), POINTER :: RINFOG
      INTEGER, DIMENSION(:), POINTER :: ICNTL
      INTEGER IRANK
      INTEGER   :: LP, MP, MPG
      LOGICAL   :: PROK, PROKG, LPOK
      LOGICAL   :: I_AM_SLAVE, PERLU_ON, PRINT_MAXAVG 
      LOGICAL   :: SUM_OF_PEAKS, PRINT_NODEINFO
      INTEGER    :: TOTAL_MBYTES_UNDER_L0
      INTEGER(8) :: TOTAL_BYTES_UNDER_L0
      INTEGER :: NBSTATS_I4, NBSTATS_I8
      PARAMETER (NBSTATS_I4=4, NBSTATS_I8=24)
      INTEGER, ALLOCATABLE, TARGET, DIMENSION(:) :: TNSTK_afterL0
      INTEGER, ALLOCATABLE, DIMENSION(:) :: FLAGGED_LEAVES
      INTEGER(8) :: PEAK_UNDER_L0, PEAK_ABOVE_L0
      INTEGER(8) :: SUM_NRLADU, MAX_NRLADU, MIN_NRLADU,
     &              SUM_NRLADU_if_LR_LU, 
     &              SUM_NRLADULR_UD, SUM_NRLADULR_WC,
     &              SUM_NRLNEC, SUM_NRLNEC_ACTIVE, 
     &              MIN_NRLNEC
      INTEGER ::    SUM_NIRADU,
     &              SUM_NIRADU_OOC, 
     &              SUM_NIRNEC, SUM_NIRNEC_OOC
      INTEGER :: LIPOOL_local
      INTEGER, ALLOCATABLE, DIMENSION(:), TARGET :: IPOOL
      INTEGER                                    :: I, LIPOOL
      INTEGER(4) :: I4
      INTEGER, POINTER, DIMENSION(:) :: NE_STEPSPTR
      INTEGER, POINTER, DIMENSION(:) :: IPOOLPTR 
      LOGICAL    :: BDUMMY
      INTEGER(8) :: K8_33relaxed, K8_34relaxed, K8_35relaxed,
     &              K8_50relaxed
      INTEGER(8) :: ENTRIES_IN_FACTORS_LOC_MASTERS
      INTEGER    :: OOC_STRAT, BLR_STRAT, IDUMMY, ISTEP, NBNODES_BLR
      INTEGER(8) :: TOTAL_BYTES, ITMP8
      INTEGER(8) ::  MAX_FRONT_SURFACE_LOCAL_L0, 
     &               MAX_SIZE_FACTOR_L0,
     &               ENTRIES_IN_FACTORS_UNDER_L0, 
     &               ENTRIES_IN_FACTORS_MASTERS_LO
      INTEGER             :: MAXFR_UNDER_L0
      DOUBLE PRECISION    :: COST_SUBTREES_UNDER_L0, OPSA_UNDER_L0
      INTEGER(8)::  SIZECB_UNDER_L0, SIZECB_UNDER_L0_IF_LRCB
      LOGICAL   :: ABOVE_L0
      INTEGER :: locMYID_NODES
      LOGICAL, POINTER :: locI_AM_CAND(:)
      INTEGER LOCAL_M, LOCAL_N
      INTEGER, EXTERNAL :: MUMPS_NUMROC
      INTEGER SBUF_SEND_FR, SBUF_REC_FR
      INTEGER SBUF_SEND_LR, SBUF_REC_LR
      INTEGER TOTAL_MBYTES
      INTEGER(8) SUM_KEEP811_THIS_NODE, MAX_SUM_KEEP811_THIS_NODE
      INTEGER SUM_INFO15_THIS_NODE, MAX_SUM_INFO15_THIS_NODE
      INTEGER SUM_INFO17_THIS_NODE, MAX_SUM_INFO17_THIS_NODE
      INTEGER(8) SBUF_RECOLD8, MIN_BUF_SIZE8
      LOGICAL UPDATE_BUFFER
      INTEGER MIN_BUF_SIZE, SIZE_DESC_BANDE,
     &        MaxBlocSize_FR, MaxBlocSize_BLR,
     &        MIN_BUF_SIZE_FR, MIN_BUF_SIZE_BLR
      INTEGER(8) MAX_SIZE_FACTOR_TMP, KEEP26_I8_TMP
      KEEP   => id%KEEP
      KEEP8  => id%KEEP8
      INFO   => id%INFO
      RINFO  => id%RINFO
      INFOG  => id%INFOG
      RINFOG => id%RINFOG
      ICNTL  => id%ICNTL
      LP  = id%ICNTL( 1 )
      MP  = id%ICNTL( 2 )
      MPG = id%ICNTL( 3 )
      LPOK  = ((LP.GT.0).AND.(id%ICNTL(4).GE.1))
      PROK  = (( MP  .GT. 0 ).AND.(id%ICNTL(4).GE.2))
      PROKG = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
      PROKG = (PROKG.AND.(id%ICNTL(4).GE.2))
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     &     ( id%MYID .eq. MASTER .AND.
     &     id%KEEP(46) .eq. 1 ) )
      PRINT_MAXAVG = .NOT.(id%NSLAVES.EQ.1 .AND. KEEP(46).EQ.1)
      PRINT_NODEINFO = PRINT_MAXAVG .AND. id%NPROCS .NE. id%KEEP(412)
      IDUMMY = 1
      BDUMMY = .FALSE.
      IF ( I_AM_SLAVE ) THEN
         locI_AM_CAND => id%I_AM_CAND
         locMYID_NODES = id%MYID_NODES
            IF ( idintr%root%yes ) THEN
              LOCAL_M = MUMPS_NUMROC(
     &          id%ND_STEPS(id%STEP(KEEP(38))),
     &          idintr%root%MBLOCK, idintr%root%MYROW, 0,
     &          idintr%root%NPROW )
              LOCAL_M = max(1, LOCAL_M)
              LOCAL_N = MUMPS_NUMROC(
     &          id%ND_STEPS(id%STEP(KEEP(38))),
     &          idintr%root%NBLOCK, idintr%root%MYCOL, 0,
     &          idintr%root%NPCOL )
            ELSE
              LOCAL_M = 0
              LOCAL_N = 0
            END IF
            IF  ( KEEP(60) .EQ. 2 .OR. KEEP(60) .EQ. 3 ) THEN
               id%SCHUR_MLOC=LOCAL_M
               id%SCHUR_NLOC=LOCAL_N
               idintr%root%SCHUR_MLOC=LOCAL_M
               idintr%root%SCHUR_NLOC=LOCAL_N
            ENDIF
            IF ( .NOT. associated(id%CANDIDATES)) THEN
               ALLOCATE(id%CANDIDATES(id%NSLAVES+1,1), stat=allocok)
               IF (allocok .gt.0) THEN
                     IF ( LPOK ) THEN
                        WRITE(LP, 150) 'CANDIDATES'
                     END IF
                     INFO(1)= -7
                     INFO(2)= id%NSLAVES+1
               ENDIF
            ENDIF
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &             id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
      IF (KEEP(400) .GT. 0 ) THEN 
          IF ( I_AM_SLAVE ) THEN
             CALL MUMPS_ANA_L0_OMP(
     &       KEEP(400), id%N, KEEP(28),
     &       KEEP(50), id%NSLAVES, id%DAD_STEPS, id%FRERE_STEPS, 
     &       id%FILS, id%NE_STEPS, id%ND_STEPS, id%STEP, 
     &       id%PROCNODE_STEPS, KEEP, KEEP8, locMYID_NODES, 
     &       id%NA, id%LNA, "DMUMPS"(1:1),
     &       id%LPOOL_B_L0_OMP, id%IPOOL_B_L0_OMP,
     &       id%LPOOL_A_L0_OMP, id%IPOOL_A_L0_OMP,
     &       id%L_VIRT_L0_OMP,id%VIRT_L0_OMP, id%VIRT_L0_OMP_MAPPING,
     &       id%L_PHYS_L0_OMP,id%PHYS_L0_OMP, id%PERM_L0_OMP,
     &       id%PTR_LEAFS_L0_OMP,
     &       id%INFO, id%ICNTL)
             IF (id%INFO(1) .GE. 0) THEN
               ALLOCATE(
     &          id%I4_L0_OMP(NBSTATS_I4, KEEP(400)), 
     &          id%I8_L0_OMP(NBSTATS_I8, KEEP(400)), 
     &          TNSTK_afterL0(KEEP(28)),
     &          stat=allocok)
               IF (allocok .gt.0) THEN
                 IF ( LPOK ) THEN
                     WRITE(LP, 150) 'L0_OMP stats'
                 END IF
                 INFO(1)= -7
                 INFO(2)= NBSTATS_I4* KEEP(400) +
     &                    NBSTATS_I8* KEEP(400)*KEEP(10)
     &                   + KEEP(28)
               ENDIF
             ENDIF
          ELSE
            ALLOCATE(id%I4_L0_OMP(1,1), id%I8_L0_OMP(1,1), stat=allocok)
            IF (allocok.gt.0) THEN
              INFO(1)= -7
              INFO(2)= 2 
            ENDIF
          ENDIF
          CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &             id%COMM, id%MYID )
          IF ( INFO(1).LT.0 ) GOTO 500
          IF ( I_AM_SLAVE ) THEN
               CALL DMUMPS_ANA_DISTM_UNDERL0OMP(
     &           id%LPOOL_B_L0_OMP, id%IPOOL_B_L0_OMP(1),
     &           id%L_VIRT_L0_OMP,
     &           id%VIRT_L0_OMP(1), id%VIRT_L0_OMP_MAPPING(1),
     &           id%L_PHYS_L0_OMP, id%PHYS_L0_OMP(1), id%PERM_L0_OMP(1),
     &           id%PTR_LEAFS_L0_OMP(1),
     &           id%KEEP(1), id%N, id%NE_STEPS(1), id%STEP(1),
     &           id%FRERE_STEPS(1), id%FILS(1), id%DAD_STEPS(1), 
     &           id%ND_STEPS(1), 
     &           locMYID_NODES, id%PROCNODE_STEPS(1),
     &           id%I4_L0_OMP(1,1), NBSTATS_I4, 
     &           id%I8_L0_OMP(1,1), NBSTATS_I8, KEEP(400),
     &           SIZECB_UNDER_L0, SIZECB_UNDER_L0_IF_LRCB,
     &           TNSTK_afterL0, MAXFR_UNDER_L0, 
     &           MAX_FRONT_SURFACE_LOCAL_L0, MAX_SIZE_FACTOR_L0, 
     &           ENTRIES_IN_FACTORS_UNDER_L0,
     &           ENTRIES_IN_FACTORS_MASTERS_LO,
     &           COST_SUBTREES_UNDER_L0, OPSA_UNDER_L0, NBNODES_BLR,
     &           INFO(1), INFO(2)
     &           )
               CALL MPI_ALLREDUCE (NBNODES_BLR, KEEP(470), 1, 
     &            MPI_INTEGER, MPI_SUM, id%COMM_NODES, IERR)
          ENDIF
      ELSE
          IF ( I_AM_SLAVE ) THEN
              id%LPOOL_B_L0_OMP = 1
              id%LPOOL_A_L0_OMP = 1
              id%L_VIRT_L0_OMP = 1
              id%L_PHYS_L0_OMP = 1
              id%THREAD_LA = -1_8
              ALLOCATE ( id%VIRT_L0_OMP ( id%L_VIRT_L0_OMP ), 
     &              id%VIRT_L0_OMP_MAPPING ( id%L_VIRT_L0_OMP ), 
     &              id%PERM_L0_OMP ( id%L_PHYS_L0_OMP ),
     &              id%PTR_LEAFS_L0_OMP ( id%L_PHYS_L0_OMP + 1 ),
     &              id%IPOOL_B_L0_OMP ( id%LPOOL_B_L0_OMP ),
     &              id%IPOOL_A_L0_OMP ( id%LPOOL_A_L0_OMP ),
     &              id%PHYS_L0_OMP( id%L_PHYS_L0_OMP ),
     &              id%I4_L0_OMP(1,1), id%I8_L0_OMP(1,1), stat=allocok)
              IF (allocok .gt.0) THEN
                 IF ( LPOK ) THEN
                    WRITE(LP, 150) 'Allocation error in multicore'
                 END IF
                 INFO(1)= -7
                 INFO(2)= id%L_VIRT_L0_OMP
     &                + id%L_PHYS_L0_OMP
     &                + id%L_PHYS_L0_OMP + 1
     &                + id%LPOOL_B_L0_OMP
     &                + id%LPOOL_A_L0_OMP
     &                + id%L_PHYS_L0_OMP + 1 + KEEP(10)
              ENDIF
          ELSE
            ALLOCATE(id%I4_L0_OMP(1,1), id%I8_L0_OMP(1,1), stat=allocok)
            IF (allocok.gt.0) THEN
              INFO(1)= -7
              INFO(2)= 2 
            ENDIF
          ENDIF
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &      id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
      IF ( I_AM_SLAVE ) THEN
          IF (KEEP(400).GT.0) THEN
              IF (id%NSLAVES .GT.1)  THEN
                 ALLOCATE (FLAGGED_LEAVES(KEEP(28)),
     &                stat=allocok)
                 IF (allocok .gt.0) THEN
                      IF ( LPOK ) THEN
                          WRITE(LP, 150) 'L0_OMP FLAGGED LEAVES'
                      END IF
                      INFO(1)= -7
                      INFO(2)= KEEP(28)
                 ENDIF
               ENDIF
          ENDIF
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &      id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
      IF ( I_AM_SLAVE ) THEN
          IF (KEEP(400).GT.0) THEN
              IF (id%NSLAVES .GT.1)  THEN
                 LIPOOL_local=
     &              id%IPOOL_A_L0_OMP(id%LPOOL_A_L0_OMP)
                 CALL DMUMPS_PREP_ANA_DISTM_ABOVEL0(
     &               id%N, id%NSLAVES, id%COMM_NODES, id%MYID_NODES, 
     &               id%STEP(1), id%DAD_STEPS(1),id%ICNTL,LP,LPOK,
     &               id%INFO, 
     &               id%PHYS_L0_OMP(1), id%L_PHYS_L0_OMP, 
     &               id%IPOOL_A_L0_OMP(1), LIPOOL_local,
     &               id%KEEP, TNSTK_afterL0,
     &               FLAGGED_LEAVES
     &              )
                 IF ( INFO(1).LT.0 ) GOTO 75
                 LIPOOL= 0
                 DO ISTEP=1,KEEP(28)
                  IF (FLAGGED_LEAVES(ISTEP).GT.0) LIPOOL=LIPOOL+1
                 ENDDO
                 ALLOCATE( IPOOL(max(LIPOOL,1)), 
     &                          stat=allocok)
                 IF (allocok .gt.0) THEN
                          IF ( LPOK ) THEN
                             WRITE(LP, 150) 'Allocation IPOOL'
                          END IF
                          INFO(1)= -7
                          INFO(2)=  LIPOOL
                 ENDIF
              ELSE
                 LIPOOL = id%IPOOL_A_L0_OMP(id%LPOOL_A_L0_OMP)
              ENDIF
          ELSE
           LIPOOL = id%NA(1)
          ENDIF
      ENDIF
 75   CONTINUE
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &      id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
      IF ( I_AM_SLAVE ) THEN
          IF (KEEP(400) .GT. 0 ) THEN 
            IF (id%NSLAVES .GT.1)  THEN
              IF (LIPOOL .GT.0) THEN
                I =LIPOOL
                DO ISTEP=1, KEEP(28)
                 IF (FLAGGED_LEAVES(ISTEP).GT.0)  THEN
                   IPOOL(I) = FLAGGED_LEAVES(ISTEP)
                   I=I-1
                 ENDIF
                 IF (I.EQ.0) EXIT 
                ENDDO
              ENDIF
              DEALLOCATE(FLAGGED_LEAVES)
              IPOOLPTR => IPOOL
            ELSE
               IPOOLPTR => id%IPOOL_A_L0_OMP
            ENDIF
            ABOVE_L0    =.TRUE.
            NE_STEPSPTR => TNSTK_afterL0(1:KEEP(28))
          ELSE
             IPOOLPTR => id%NA(3:3+max(LIPOOL,1)-1)
             ABOVE_L0                   =.FALSE.
             SIZECB_UNDER_L0            = 0_8
             SIZECB_UNDER_L0_IF_LRCB    = 0_8
             MAX_FRONT_SURFACE_LOCAL_L0 = 0_8
             MAX_SIZE_FACTOR_L0         = 0_8
             ENTRIES_IN_FACTORS_UNDER_L0= 0_8
             ENTRIES_IN_FACTORS_MASTERS_LO = 0_8
             MAXFR_UNDER_L0             = 0
             COST_SUBTREES_UNDER_L0     = 0.0D0
             OPSA_UNDER_L0              = 0.0D0
             NE_STEPSPTR => id%NE_STEPS
          ENDIF
           KEEP(139) = MAXFR_UNDER_L0
           CALL DMUMPS_ANA_DISTM( locMYID_NODES, id%N, id%STEP(1),
     & id%FRERE_STEPS(1), id%FILS(1), IPOOLPTR(1), LIPOOL, NE_STEPSPTR
     & (1), id%DAD_STEPS(1), id%ND_STEPS(1), id%PROCNODE_STEPS(1),
     & id%NSLAVES, ABOVE_L0,SIZECB_UNDER_L0,SIZECB_UNDER_L0_IF_LRCB,
     & MAXFR_UNDER_L0, MAX_FRONT_SURFACE_LOCAL_L0, MAX_SIZE_FACTOR_L0, 
     & ENTRIES_IN_FACTORS_UNDER_L0, ENTRIES_IN_FACTORS_MASTERS_LO,
     & COST_SUBTREES_UNDER_L0, OPSA_UNDER_L0, KEEP8(53), KEEP8(54),
     & KEEP8(11), KEEP(26), KEEP(15), KEEP8(12),  KEEP8(14),  
     & KEEP8(32), KEEP8(33), KEEP8(34), KEEP8(35), KEEP8(50), 
     & KEEP8(36), KEEP8(47), KEEP8(37), KEEP8(38), KEEP8(39),
     & KEEP8(40), KEEP8(41), KEEP8(42), KEEP8(43), KEEP8(44), KEEP8(45),
     & KEEP8(46), KEEP8(51), KEEP8(52), KEEP(224),KEEP(225),KEEP(27),
     & RINFO(1),id%CNTL(1), KEEP(1), KEEP8(1), LOCAL_M, LOCAL_N,
     & SBUF_RECOLD8, SBUF_SEND_FR, SBUF_REC_FR, SBUF_SEND_LR,
     & SBUF_REC_LR, id%COST_SUBTREES, KEEP(28), locI_AM_CAND(1),
     & max(KEEP(56),1), id%ISTEP_TO_INIV2(1), id%CANDIDATES(1,1), 
     & INFO(1), INFO(2), KEEP8(15),MAX_SIZE_FACTOR_TMP, 
     & KEEP8(9), ENTRIES_IN_FACTORS_LOC_MASTERS,
     & idintr%root%yes, idintr%root%NPROW, idintr%root%NPCOL
     &           )
           IF (ALLOCATED(IPOOL)) DEALLOCATE(IPOOL)
           NULLIFY(NE_STEPSPTR,IPOOLPTR)
           IF (KEEP(400) .GT. 0) THEN 
             DEALLOCATE (TNSTK_afterL0)
             SUM_NIRNEC     = 0
             SUM_NIRADU     = 0
             SUM_NIRADU_OOC = 0
             SUM_NIRNEC_OOC = 0
             DO I=1, KEEP(400)
                SUM_NIRADU = SUM_NIRADU        + id%I4_L0_OMP(1,I)
                SUM_NIRNEC = SUM_NIRNEC        + id%I4_L0_OMP(2,I)
                SUM_NIRADU_OOC = SUM_NIRADU_OOC+ id%I4_L0_OMP(3,I)
                SUM_NIRNEC_OOC = SUM_NIRNEC_OOC+ id%I4_L0_OMP(4,I)
             ENDDO
             KEEP(26)  = KEEP(26)  + SUM_NIRADU
             KEEP(224) = KEEP(224) + SUM_NIRADU_OOC
             KEEP(15) = max(KEEP(15),KEEP(26))
             KEEP(225) = max(KEEP(225),KEEP(224))
             KEEP(137) = SUM_NIRNEC
             KEEP(138) = SUM_NIRNEC_OOC
             SUM_NIRNEC =  int( 
     &              (dble(SUM_NIRNEC)*dble(KEEP(34)))/dble(KEEP(35))
     &                     )
             SUM_NIRNEC_OOC =  int( 
     &              (dble(SUM_NIRNEC_OOC)*dble(KEEP(34)))/dble(KEEP(35))
     &                     )
             MAX_NRLADU          = 0_8
             MIN_NRLADU          = id%I8_L0_OMP(1,1)
             SUM_NRLADU          = 0_8
             SUM_NRLNEC          = 0_8
             MIN_NRLNEC          = huge(MIN_NRLNEC)
             SUM_NRLNEC_ACTIVE   = 0_8
             SUM_NRLADU_if_LR_LU = 0_8
             SUM_NRLADULR_UD     = 0_8
             SUM_NRLADULR_WC     = 0_8
             DO I=1, KEEP(400)
               MIN_NRLADU = min(MIN_NRLADU, id%I8_L0_OMP(1,I))
               MAX_NRLADU = max(MAX_NRLADU, id%I8_L0_OMP(1,I))
               SUM_NRLADU = SUM_NRLADU + id%I8_L0_OMP(1,I)
               SUM_NRLNEC = SUM_NRLNEC + id%I8_L0_OMP(2,I)
               MIN_NRLNEC = min(MIN_NRLNEC, id%I8_L0_OMP(2,I))
               SUM_NRLNEC_ACTIVE = SUM_NRLNEC_ACTIVE +
     &                      id%I8_L0_OMP(3,I)
               SUM_NRLADU_if_LR_LU = SUM_NRLADU_if_LR_LU +
     &                      id%I8_L0_OMP(4,I)
               SUM_NRLADULR_UD = SUM_NRLADULR_UD  +
     &                      id%I8_L0_OMP(9,I)
               SUM_NRLADULR_WC = SUM_NRLADULR_WC  +
     &                      id%I8_L0_OMP(10,I)
             ENDDO
             KEEP8(81) = KEEP8(11) 
             KEEP8(11) = KEEP8(11) + SUM_NRLADU
             KEEP8(82) = KEEP8(32)
             KEEP8(32) = KEEP8(32) + SUM_NRLADU_if_LR_LU
             PEAK_UNDER_L0 = SUM_NRLNEC + MIN_NRLNEC + 
     &        int(
     &        (dble(id%N*KEEP(400))*dble(KEEP(34)))/dble(KEEP(35)),
     &               8)
             PEAK_ABOVE_L0 =  KEEP8(53)+ SUM_NRLADU + 
     &            
     &            max( int(SBUF_SEND_FR,8),  100000_8) +
     &            
     &             int( 
     &              (dble(KEEP(15))*dble(KEEP(34)))/dble(KEEP(35)),
     &               8)
             KEEP8(53) = KEEP8(53)+ SUM_NRLADU
             KEEP8(40) = KEEP8(40)+
     &              SUM_NRLADU_if_LR_LU + SUM_NRLADULR_UD
             KEEP8(41) = KEEP8(41)+ SUM_NRLADULR_UD
             KEEP8(42) = KEEP8(42)+ SUM_NRLADULR_WC 
             KEEP8(43) = KEEP8(43)+
     &              SUM_NRLADU_if_LR_LU + SUM_NRLADULR_UD
             KEEP8(44) = KEEP8(44)+
     &              SUM_NRLADU_if_LR_LU + SUM_NRLADULR_WC
             KEEP8(45) = KEEP8(45)+ SUM_NRLADULR_UD
             KEEP8(46) = KEEP8(46)+ SUM_NRLADULR_WC
             KEEP8(51) = KEEP8(51)+ SUM_NRLADU
             KEEP8(52) = KEEP8(52)+  SUM_NRLADULR_UD
           ELSE
             KEEP(137)=0
             KEEP(138)=0
           ENDIF
            id%DKEEP(15)    = RINFO(1)/1000000.0D0
            IF(ASSOCIATED(locI_AM_CAND)) NULLIFY(locI_AM_CAND)
            id%MAX_SURF_MASTER = KEEP8(15)
            KEEP8(19)=MAX_SIZE_FACTOR_TMP
            KEEP( 29 ) = KEEP(15) + 3* max(KEEP(12),10)
     &           * ( KEEP(15) / 100 + 1)
            INFO( 19 ) = KEEP(225) + 3* max(KEEP(12),10)
     &           * ( KEEP(225) / 100 + 1)
            KEEP8(13)  = KEEP8(12) + int(KEEP(12),8) *
     &           ( KEEP8(12) / 100_8 + 1_8 )
            KEEP8(17)  = KEEP8(14) + int(KEEP(12),8) *
     &           ( KEEP8(14) /100_8 +1_8)
            K8_33relaxed  = KEEP8(33) + int(KEEP(12),8) *
     &           ( KEEP8(33) /100_8 +1_8)
            K8_34relaxed  = KEEP8(34) + int(KEEP(12),8) *
     &           ( KEEP8(34) /100_8 +1_8)
            K8_35relaxed  = KEEP8(35) + int(KEEP(12),8) *
     &           ( KEEP8(35) /100_8 +1_8)
            K8_50relaxed  = KEEP8(50) + int(KEEP(12),8) *
     &           ( KEEP8(50) /100_8 +1_8)
         CALL MUMPS_ALLREDUCEI8 ( SBUF_RECOLD8, KEEP8(22), MPI_MAX,
     &                            id%COMM_NODES )
         IF ( (id%NSLAVES.GT.1)
     &      ) THEN
             SBUF_SEND_FR = max(SBUF_SEND_FR,KEEP(27))
             SBUF_SEND_LR = max(SBUF_SEND_LR,KEEP(27))
             SBUF_REC_FR  = max(SBUF_REC_FR ,KEEP(27))
             SBUF_REC_LR  = max(SBUF_REC_LR ,KEEP(27))
         ENDIF
         CALL MPI_ALLREDUCE (SBUF_REC_FR, KEEP(44), 1, 
     &        MPI_INTEGER, MPI_MAX,
     &        id%COMM_NODES, IERR)
         CALL MPI_ALLREDUCE (SBUF_REC_LR, KEEP(380), 1, 
     &        MPI_INTEGER, MPI_MAX,
     &        id%COMM_NODES, IERR)
         IF (KEEP(48)==5) THEN
           KEEP(43)  = KEEP(44)
           KEEP(379) = KEEP(380)
#if defined(MPI_TO_K_OMP)
           IF ( KEEP(381).GT. 1)  THEN
             IF (id%MTKO_PROCS_MAP(id%MYID_NODES) 
     &              .NE.id%MYID_NODES ) THEN
              KEEP(43)  = SBUF_SEND_FR
              KEEP(379) = SBUF_SEND_LR
             ELSE
             ENDIF
           ENDIF
#endif
         ELSE
           KEEP(43)=SBUF_SEND_FR
           KEEP(379)=SBUF_SEND_LR
         ENDIF
         UPDATE_BUFFER = .TRUE.
#if defined(MPI_TO_K_OMP)
         IF ( KEEP(381).GT. 1)  THEN
             IF (id%MTKO_PROCS_MAP(id%MYID_NODES) 
     &              .NE.id%MYID_NODES ) THEN
               UPDATE_BUFFER = .FALSE.
             ENDIF
          ENDIF
#endif
         MIN_BUF_SIZE8 = KEEP8(22) / int(KEEP(238),8)
         MIN_BUF_SIZE8 = min(MIN_BUF_SIZE8,
     &                       int(huge(I4),8)/int(KEEP(35),8) )
         MIN_BUF_SIZE  = max(int( MIN_BUF_SIZE8 ), KEEP(127))
         SIZE_DESC_BANDE=(11+id%NSLAVES+KEEP(127)*2)
         MaxBlocSize_FR  = min (KEEP(420), KEEP(127))
         MaxBlocSize_FR  = MaxBlocSize_FR*MaxBlocSize_FR
         MaxBlocSize_BLR = min (KEEP(142), KEEP(127))
         MaxBlocSize_BLR = MaxBlocSize_BLR*MaxBlocSize_BLR
         MIN_BUF_SIZE_FR  = MIN_BUF_SIZE
         MIN_BUF_SIZE_BLR = MIN_BUF_SIZE
         MIN_BUF_SIZE_FR = min ( MIN_BUF_SIZE_FR, 
     &       int ( min (
     &              dble(KEEP(44)) * 
     &             (dble(abs(KEEP(180))) / dble(100)) ,
     &             dble(huge(I4))/dble(KEEP(35))
     &                 )  )
     &                          )
         MIN_BUF_SIZE_BLR = min ( MIN_BUF_SIZE_BLR, 
     &       int ( min (
     &              dble(KEEP(44)) * 
     &             (dble(abs(KEEP(181))) / dble(100)) ,
     &             dble(huge(I4))/dble(KEEP(35))
     &                 )  )
     &                          )
         IF (KEEP(50).EQ.0) THEN
           KEEP(43) = max(
     &        min(KEEP(43),MaxBlocSize_FR*max(KEEP(171),3)), 
     &        int(KEEP(43)/KEEP(172)) )
           KEEP(44) = max(
     &        min(KEEP(44), MaxBlocSize_FR*max(KEEP(171),3)), 
     &        int(KEEP(44)/KEEP(172)) )
         ELSE
           KEEP(43) = max(
     &        min(KEEP(43),MaxBlocSize_FR*max(KEEP(171),3)), 
     &        int((KEEP(43)*KEEP(178))/KEEP(172)) )
           KEEP(44) = max(
     &        min(KEEP(44), MaxBlocSize_FR*max(KEEP(171),3)), 
     &        int((KEEP(44)*KEEP(178))/KEEP(172)) )
         ENDIF
         KEEP(379) = max(
     &        min(KEEP(379), MaxBlocSize_BLR*max(KEEP(171),3)), 
     &        int(KEEP(379)/KEEP(172)) )
         KEEP(380) = max(
     &       min(KEEP(380),MaxBlocSize_BLR*max(KEEP(171),3)), 
     &       int(KEEP(380)/KEEP(172)) )
         IF (UPDATE_BUFFER) THEN
           KEEP(43) = max(KEEP(43),MIN_BUF_SIZE_FR)  +
     &         (SIZE_DESC_BANDE * KEEP(34))/KEEP(35)
           KEEP(379)= max(KEEP(379),MIN_BUF_SIZE_BLR)+
     &         (SIZE_DESC_BANDE * KEEP(34))/KEEP(35)
         ENDIF 
         IF ( (KEEP(38).NE.0) .OR. UPDATE_BUFFER) THEN
           KEEP(44) = max(KEEP(44),MIN_BUF_SIZE_FR)  +
     &         (SIZE_DESC_BANDE * KEEP(34))/KEEP(35)
           KEEP(380)= max(KEEP(380),MIN_BUF_SIZE_BLR)+
     &         (SIZE_DESC_BANDE * KEEP(34))/KEEP(35)
         ENDIF
          IF ( int(KEEP(43),8)*int(KEEP(35),8)  .GE.
     &                 int(huge(I4)-100,8)) THEN   
                KEEP(43) = huge(KEEP(43))-100
          ENDIF
          IF ( int(KEEP(44),8)*int(KEEP(35),8)  .GE.
     &                 int(huge(I4)-100,8)) THEN   
                KEEP(44) = huge(KEEP(44))-100
          ENDIF
          IF ( int(KEEP(379),8)*int(KEEP(35),8)  .GE.
     &                 int(huge(I4)-100,8)) THEN   
                KEEP(379) = huge(KEEP(379))-100
          ENDIF
          IF ( int(KEEP(380),8)*int(KEEP(35),8)  .GE.
     &                 int(huge(I4)-100,8)) THEN   
                KEEP(380) = huge(KEEP(380))-100
          ENDIF
            IF ( PROK ) THEN
               WRITE(MP,'(A,I16) ') 
     &              ' Estimated INTEGER space for factors         :',
     &              KEEP(26)
               WRITE(MP,'(A,I16) ') 
     &              ' INFO(3), est. real space to store factors   :',
     &              KEEP8(11)
               WRITE(MP,'(A,I16) ') 
     &              ' Estimated number of entries in factors      :',
     &              KEEP8(9)
               WRITE(MP,'(A,I16) ') 
     &              ' Current value of space relaxation parameter :',
     &              KEEP(12)
               WRITE(MP,'(A,I16) ') 
     &              ' Estimated size of IS (In Core factorization):',
     &              KEEP(29)
               WRITE(MP,'(A,I16) ') 
     &              ' Estimated size of S  (In Core factorization):',
     &              KEEP8(13)
               WRITE(MP,'(A,I16) ') 
     &              ' Estimated size of S  (OOC factorization)    :',
     &              KEEP8(17)
            END IF
      ELSE
         ENTRIES_IN_FACTORS_LOC_MASTERS = 0_8
         KEEP8(13) = 0_8
         KEEP(29) = 0
         KEEP8(17)= 0_8
         INFO(19) = 0
         KEEP8(11) = 0_8
         KEEP8(81) = 0_8
         KEEP8(82) = 0_8
         KEEP(26) = 0
         KEEP(27) = 0
         RINFO(1) = 0.0D0
         K8_33relaxed = 0_8
         K8_34relaxed = 0_8
         K8_35relaxed = 0_8
         K8_50relaxed = 0_8
         IF (KEEP(400) .GT.0) THEN
             SUM_NIRNEC     = 0
             SUM_NIRADU     = 0
             SUM_NIRADU_OOC = 0
             SUM_NIRNEC_OOC = 0
             MAX_NRLADU          = 0_8
             MIN_NRLADU          = 0_8
             SUM_NRLADU          = 0_8
             SUM_NRLNEC          = 0_8
             SUM_NRLNEC_ACTIVE   = 0_8
             SUM_NRLADU_if_LR_LU = 0_8
             SUM_NRLADULR_UD     = 0_8
             SUM_NRLADULR_WC     = 0_8
         ENDIF
      END IF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      CALL MUMPS_ALLREDUCEI8( ENTRIES_IN_FACTORS_LOC_MASTERS, 
     &     KEEP8(109), MPI_SUM, id%COMM)
      CALL MUMPS_ALLREDUCEI8( KEEP8(19), KEEP8(119),
     &     MPI_MAX, id%COMM)
      CALL MPI_ALLREDUCE( KEEP(27), KEEP(127), 1,
     &     MPI_INTEGER, MPI_MAX,
     &     id%COMM, IERR)
      KEEP26_I8_TMP = int(KEEP(26),8)
      CALL MUMPS_ALLREDUCEI8( KEEP26_I8_TMP, 
     &     KEEP8(129), MPI_SUM, id%COMM)
      CALL MUMPS_REDUCEI8( KEEP8(11),
     &   KEEP8(111), MPI_SUM,
     &   MASTER, id%COMM )
      CALL MUMPS_SETI8TOI4( KEEP8(111), INFOG(3) )
      RINFO(5) = dble(KEEP8(32)
     &               *int(KEEP(35),8))/1D6
      CALL MUMPS_REDUCEI8( KEEP8(32),
     &                   ITMP8, MPI_SUM,
     &                   MASTER, id%COMM )
      IF (id%MYID.EQ.MASTER) THEN
       RINFOG(15) = dble(ITMP8*int(KEEP(35),8))/1D6
      ENDIF
      CALL MPI_ALLREDUCE( RINFO(1), RINFOG(1), 1,
     &     MPI_DOUBLE_PRECISION, MPI_SUM,
     &     id%COMM, IERR)
      CALL MUMPS_SETI8TOI4( KEEP8(11), INFO(3) )
      INFO ( 4 ) = KEEP(  26 )
      INFO ( 5 ) = KEEP(  27 )
      INFO ( 7 ) = KEEP(  29 )
      CALL MUMPS_SETI8TOI4( KEEP8(13), INFO(8) )
      CALL MUMPS_SETI8TOI4( KEEP8(17), INFO(20) )
      CALL MUMPS_SETI8TOI4( KEEP8(9), INFO(24) )
      CALL MUMPS_SETI8TOI4( K8_33relaxed, INFO(29) )
      CALL MUMPS_SETI8TOI4( K8_34relaxed, INFO(32) )
      CALL MUMPS_SETI8TOI4( K8_35relaxed, INFO(33) )
      CALL MUMPS_SETI8TOI4( K8_50relaxed, INFO(36) )
      CALL MUMPS_SETI8TOI4( KEEP8(129), INFOG(4) )
      INFOG( 5 ) = KEEP( 127 )
      CALL MUMPS_SETI8TOI4( KEEP8(109), INFOG(20) )
      CALL DMUMPS_DIAG_ANA(id%N, NBLK,
     &     id%MYID, id%COMM, KEEP(1), KEEP8(1),
     &     INFO(1), INFOG(1), RINFO(1), RINFOG(1), ICNTL(1),
     &     id%SIZE_SCHUR )
      IF (PROK) WRITE( MP, 112 )
      IF (PROKG .AND. (MPG.NE.MP)) WRITE( MPG, 112 )
      SUM_KEEP811_THIS_NODE=0_8
      CALL MPI_REDUCE( KEEP8(11), SUM_KEEP811_THIS_NODE, 1,
     &                 MPI_INTEGER8,
     &                 MPI_SUM, 0, id%KEEP(411), IERR )
      CALL MPI_REDUCE( SUM_KEEP811_THIS_NODE, MAX_SUM_KEEP811_THIS_NODE,
     &                 1, MPI_INTEGER8, MPI_MAX, 0, id%COMM, IERR )
      IF (PROKG .AND. PRINT_NODEINFO) THEN
        WRITE(MPG,'(A,I16)')
     &  '    Max. estimated space for factors per compute node       :',
     &  MAX_SUM_KEEP811_THIS_NODE 
      ENDIF
         OOC_STRAT = KEEP(201)
         BLR_STRAT = 0  
         IF (KEEP(201) .NE. -1) OOC_STRAT=0 
         PERLU_ON = .FALSE.     
         CALL DMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &        id%KEEP8(30),
     &        id%NSLAVES, TOTAL_MBYTES, .FALSE.,
     &        OOC_STRAT, BLR_STRAT, PERLU_ON, TOTAL_BYTES, 
     &        IDUMMY, BDUMMY, .FALSE., 
     &        .FALSE. 
     &        , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &        size(id%I8_L0_OMP,2) 
     &        )
            IF (KEEP(400) .GT. 0 ) THEN 
              CALL DMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &        id%KEEP8(30),
     &        id%NSLAVES, TOTAL_MBYTES_UNDER_L0, .FALSE.,
     &        OOC_STRAT, BLR_STRAT, PERLU_ON, TOTAL_BYTES_UNDER_L0, 
     &        IDUMMY, BDUMMY, .FALSE., 
     &        .TRUE. 
     &        , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &        size(id%I8_L0_OMP,2) 
     &        )
              TOTAL_MBYTES = max (TOTAL_MBYTES,TOTAL_MBYTES_UNDER_L0) 
              TOTAL_BYTES  = max (TOTAL_BYTES, TOTAL_BYTES_UNDER_L0)
            ENDIF
         KEEP8(2) = TOTAL_BYTES    
         PERLU_ON  = .TRUE.
         CALL DMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &        id%KEEP8(30),
     &        id%NSLAVES, TOTAL_MBYTES, .FALSE.,
     &        OOC_STRAT, BLR_STRAT, PERLU_ON, TOTAL_BYTES, 
     &        IDUMMY, BDUMMY, .FALSE., 
     &        .FALSE. 
     &        , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &        size(id%I8_L0_OMP,2) 
     &        )
         IF (KEEP(400) .GT. 0 ) THEN 
           CALL DMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &          id%MYID, id%N, id%NELT, id%NA(1), id%LNA, 
     &          id%KEEP8(28),
     &          id%KEEP8(30),
     &          id%NSLAVES, TOTAL_MBYTES_UNDER_L0, .FALSE.,
     &          OOC_STRAT, BLR_STRAT, PERLU_ON, TOTAL_BYTES_UNDER_L0, 
     &          IDUMMY, BDUMMY, .FALSE., 
     &          .TRUE. 
     &          , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &          size(id%I8_L0_OMP,2) 
     &        )
           TOTAL_MBYTES = max (TOTAL_MBYTES,TOTAL_MBYTES_UNDER_L0) 
           TOTAL_BYTES  = max (TOTAL_BYTES, TOTAL_BYTES_UNDER_L0)
         ENDIF
         IF ( PROK ) THEN
            WRITE(MP,'(A,I16) ')
     & ' Estimated space in MBytes for IC factorization   (INFO(15)):',
     &           TOTAL_MBYTES
         END IF
         id%INFO(15) = TOTAL_MBYTES
      CALL MUMPS_MEM_CENTRALIZE( id%MYID, id%COMM,
     &     id%INFO(15), id%INFOG(16), IRANK )
      IF ( PROKG ) THEN
       IF (PRINT_MAXAVG) THEN
         WRITE( MPG,'(A,I16) ')
     & '    Maximum estim. space in Mbytes, IC facto.    (INFOG(16)):',
     &        id%INFOG(16)
       ENDIF
       WRITE(MPG,'(A,I16) ')
     & '    Total space in MBytes, IC factorization      (INFOG(17)):'
     &        ,id%INFOG(17)
      END IF
      SUM_INFO15_THIS_NODE=0
      CALL MPI_REDUCE( INFO(15), SUM_INFO15_THIS_NODE, 1, MPI_INTEGER,
     &                 MPI_SUM, 0, id%KEEP(411), IERR )
      CALL MPI_REDUCE( SUM_INFO15_THIS_NODE, MAX_SUM_INFO15_THIS_NODE,
     &                 1, MPI_INTEGER, MPI_MAX, 0, id%COMM, IERR )
      IF ( PROKG .AND. PRINT_NODEINFO ) THEN
        WRITE(MPG,'(A,I16)')
     &  '    Max. estim. space per compute node, in MBytes, IC fact  :',
     &  MAX_SUM_INFO15_THIS_NODE
      ENDIF
         OOC_STRAT = KEEP(201)
         BLR_STRAT = 0  
#if defined(OLD_OOC_NOPANEL)
         IF (OOC_STRAT .NE. -1) OOC_STRAT=2
#else
         IF (OOC_STRAT .NE. -1) OOC_STRAT=1
#endif
         PERLU_ON = .FALSE.     
         CALL DMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &        id%KEEP8(30),
     &        id%NSLAVES, TOTAL_MBYTES, .FALSE.,
     &        OOC_STRAT, BLR_STRAT, PERLU_ON, TOTAL_BYTES, 
     &        IDUMMY, BDUMMY, .FALSE., 
     &        .FALSE. 
     &        , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &        size(id%I8_L0_OMP,2) 
     &         )
         IF (KEEP(400) .GT. 0 ) THEN 
            CALL DMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &        id%KEEP8(30),
     &        id%NSLAVES, TOTAL_MBYTES_UNDER_L0, .FALSE.,
     &        OOC_STRAT, BLR_STRAT, PERLU_ON, TOTAL_BYTES_UNDER_L0, 
     &        IDUMMY, BDUMMY, .FALSE., 
     &        .TRUE. 
     &        , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &        size(id%I8_L0_OMP,2) 
     &        )
            TOTAL_MBYTES = max (TOTAL_MBYTES,TOTAL_MBYTES_UNDER_L0) 
            TOTAL_BYTES  = max (TOTAL_BYTES, TOTAL_BYTES_UNDER_L0)
         ENDIF
         KEEP8(3) = TOTAL_BYTES
         PERLU_ON  = .TRUE.     
         CALL DMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &        id%KEEP8(30),
     &        id%NSLAVES, TOTAL_MBYTES, .FALSE.,
     &        OOC_STRAT, BLR_STRAT, PERLU_ON, TOTAL_BYTES,
     &        IDUMMY, BDUMMY, .FALSE., 
     &        .FALSE. 
     &        , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &        size(id%I8_L0_OMP,2) 
     &         )
         IF (KEEP(400) .GT. 0 ) THEN 
           CALL DMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &        id%KEEP8(30),
     &        id%NSLAVES, TOTAL_MBYTES_UNDER_L0, .FALSE.,
     &        OOC_STRAT, BLR_STRAT, PERLU_ON, TOTAL_BYTES_UNDER_L0,
     &        IDUMMY, BDUMMY, .FALSE., 
     &        .TRUE. 
     &        , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &        size(id%I8_L0_OMP,2) 
     &        )
           TOTAL_MBYTES = max (TOTAL_MBYTES,TOTAL_MBYTES_UNDER_L0) 
           TOTAL_BYTES  = max (TOTAL_BYTES, TOTAL_BYTES_UNDER_L0)
         ENDIF
         id%INFO(17) = TOTAL_MBYTES
      CALL MUMPS_MEM_CENTRALIZE( id%MYID, id%COMM,
     &     id%INFO(17), id%INFOG(26), IRANK )
      IF ( PROKG  ) THEN
       IF (PRINT_MAXAVG) THEN
         WRITE( MPG,'(A,I16) ')
     & '    Maximum estim. space in Mbytes, OOC facto.   (INFOG(26)):',
     &        id%INFOG(26)
       ENDIF
       WRITE(MPG,'(A,I16) ')
     & '    Total space in MBytes,  OOC factorization    (INFOG(27)):'
     &        ,id%INFOG(27)
      END IF
      SUM_INFO17_THIS_NODE=0
      CALL MPI_REDUCE( INFO(17), SUM_INFO17_THIS_NODE, 1, MPI_INTEGER,
     &                 MPI_SUM, 0, id%KEEP(411), IERR )
      CALL MPI_REDUCE( SUM_INFO17_THIS_NODE, MAX_SUM_INFO17_THIS_NODE,
     &                 1, MPI_INTEGER, MPI_MAX, 0, id%COMM, IERR )
      IF (PROKG .AND. PRINT_NODEINFO) THEN
        WRITE(MPG,'(A,I16)')
     &  '    Max. estim. space per compute node, in MBytes, OOC fact :',
     &  MAX_SUM_INFO17_THIS_NODE
      ENDIF
      IF (KEEP(494).NE.0) THEN
         SUM_OF_PEAKS = .TRUE.
         CALL  DMUMPS_MEM_ESTIM_BLR_ALL( SUM_OF_PEAKS, 
     &        KEEP(1), KEEP8(1),
     &        id%MYID, id%COMM, 
     &        id%N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &        id%KEEP8(30), id%NSLAVES,
     &        id%INFO, id%INFOG, PROK, MP, PROKG, MPG 
     &        , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &        size(id%I8_L0_OMP,2) 
     &        )
      END IF
 500  CONTINUE
      IF (allocated(TNSTK_afterL0)) DEALLOCATE(TNSTK_afterL0)
      IF (allocated(FLAGGED_LEAVES)) DEALLOCATE(FLAGGED_LEAVES)
      IF (INFO(1) .LT. 0) THEN
        IF (associated(id%IPOOL_B_L0_OMP)) THEN
          DEALLOCATE(id%IPOOL_B_L0_OMP)
          NULLIFY(id%IPOOL_B_L0_OMP)
        ENDIF
        IF (associated(id%IPOOL_A_L0_OMP)) THEN
          DEALLOCATE(id%IPOOL_A_L0_OMP)
          NULLIFY(id%IPOOL_A_L0_OMP)
        ENDIF
        IF (associated(id%VIRT_L0_OMP)) THEN
          DEALLOCATE(id%VIRT_L0_OMP)
          NULLIFY(id%VIRT_L0_OMP)
        ENDIF
        IF (associated(id%VIRT_L0_OMP_MAPPING)) THEN
          DEALLOCATE(id%VIRT_L0_OMP_MAPPING)
          NULLIFY(id%VIRT_L0_OMP_MAPPING)
        ENDIF
        IF (associated(id%PERM_L0_OMP)) THEN
          DEALLOCATE(id%PERM_L0_OMP)
          NULLIFY(id%PERM_L0_OMP)
        ENDIF
        IF (associated(id%PTR_LEAFS_L0_OMP)) THEN
          DEALLOCATE(id%PTR_LEAFS_L0_OMP)
          NULLIFY(id%PTR_LEAFS_L0_OMP)
        ENDIF
      ENDIF
      RETURN
 112  FORMAT(/' MEMORY ESTIMATIONS ... '/
     &   ' Estimations with standard Full-Rank (FR) factorization:')
150   FORMAT(
     & /' ** FAILURE DURING DMUMPS_ANA_DRIVER, DYNAMIC ALLOCATION OF',
     &     A30)
      END SUBROUTINE DMUMPS_ANA_COMPUTE_ESTIMATES
      SUBROUTINE DMUMPS_ANA_DISTM(MYID, N, STEP, FRERE, FILS, IPOOL,
     & LIPOOL, NE, DAD, ND, PROCNODE, SLAVEF, ABOVE_L0, SIZECB_UNDER_L0,
     & SIZECB_UNDER_L0_IF_LRCB, MAXFR_UNDER_L0,
     & MAX_FRONT_SURFACE_LOCAL_L0, MAX_SIZE_FACTOR_L0, 
     & ENTRIES_IN_FACTORS_UNDER_L0, ENTRIES_IN_FACTORS_MASTERS_LO,
     & COST_SUBTREES_UNDER_LO, OPSA_UNDER_L0, PEAK_FR, PEAK_FR_OOC, 
     & NRLADU, NIRADU, NIRNEC, NRLNEC, NRLNEC_ACTIVE, 
     & NRLADU_if_LR_LU, NRLNEC_if_LR_LU, NRLNEC_if_LR_LUCB, 
     & NRLNECOOC_if_LR_LUCB, NRLNEC_if_LR_CB, NRLADULR_UD, NRLADULR_WC,
     & NRLNECLR_CB_UD, NRLNECLR_LUCB_UD, NRLNECLR_LUCB_WC,
     & PEAK_LRLU_UD,PEAK_OOC_LRLU_UD,PEAK_OOC_LRLU_WC, PEAK_LRLUCB_UD,
     & PEAK_LRLUCB_WC,PEAK_OOC_LRLUCB_UD, PEAK_OOC_LRLUCB_WC,
     & PEAK_LRCB_UD, PEAK_OOC_LRCB_UD, NIRADU_OOC, NIRNEC_OOC, MAXFR,
     & OPSA, UU, KEEP,KEEP8, LOCAL_M, LOCAL_N, SBUF_RECOLD,
     & SBUF_SEND_FR, SBUF_REC_FR, SBUF_SEND_LR, SBUF_REC_LR, 
     & OPS_SUBTREE, NSTEPS, I_AM_CAND,NMB_PAR2, ISTEP_TO_INIV2,
     & CANDIDATES, IFLAG, IERROR, MAX_FRONT_SURFACE_LOCAL,
     & MAX_SIZE_FACTOR, ENTRIES_IN_FACTORS_LOC,
     & ENTRIES_IN_FACTORS_LOC_MASTERS, ROOT_yes, ROOT_NPROW, ROOT_NPCOL
     &     )
      USE DMUMPS_LR_CORE, ONLY : IS_FRONT_BLR_CANDIDATE
      USE MUMPS_TOOLS_COMMON_M, ONLY : MUMPS_LR_COMPRESS_CB,
     &                                 MUMPS_LR_COMPRESS_PANEL,
     &                                 MUMPS_LR_FORCE_BFR
      IMPLICIT NONE
      LOGICAL, intent(in) :: ROOT_yes
      INTEGER, intent(in) :: ROOT_NPROW, ROOT_NPCOL
      INTEGER, intent(in) :: MYID, N, LIPOOL
      LOGICAL, intent(in) :: ABOVE_L0
      INTEGER, intent(in) :: MAXFR_UNDER_L0
      INTEGER(8), intent(in) :: MAX_FRONT_SURFACE_LOCAL_L0, 
     &                          MAX_SIZE_FACTOR_L0, 
     &                          ENTRIES_IN_FACTORS_UNDER_L0,
     &                          ENTRIES_IN_FACTORS_MASTERS_LO
      DOUBLE PRECISION, intent(in) :: COST_SUBTREES_UNDER_LO, 
     &                                OPSA_UNDER_L0
      INTEGER(8), intent(inout) :: SIZECB_UNDER_L0, 
     &                             SIZECB_UNDER_L0_IF_LRCB
      INTEGER, intent(inout) ::  IFLAG, IERROR
      INTEGER  NIRADU, NIRNEC
      INTEGER(8) NRLADU, NRLNEC, NRLNEC_ACTIVE
      INTEGER(8), intent(out) :: NRLADU_if_LR_LU, 
     &  NRLADULR_UD, NRLADULR_WC,
     &  NRLNEC_if_LR_LU, NRLNEC_if_LR_LUCB, NRLNEC_if_LR_CB,
     &  NRLNECOOC_if_LR_LUCB, NRLNECLR_CB_UD,
     &  NRLNECLR_LUCB_UD, NRLNECLR_LUCB_WC
      INTEGER(8), intent(out)::
     &  PEAK_FR, PEAK_FR_OOC,
     &  PEAK_LRLU_UD, 
     &  PEAK_OOC_LRLU_UD, PEAK_OOC_LRLU_WC,
     &  PEAK_LRLUCB_UD, PEAK_LRLUCB_WC,
     &  PEAK_OOC_LRLUCB_UD, PEAK_OOC_LRLUCB_WC,
     &  PEAK_LRCB_UD, PEAK_OOC_LRCB_UD
      INTEGER(8) NRLADU_CURRENT, NRLADU_ROOT_3
      INTEGER NIRADU_OOC, NIRNEC_OOC
      INTEGER MAXFR, NSTEPS
      INTEGER(8) MAX_FRONT_SURFACE_LOCAL
      INTEGER STEP(N)
      INTEGER FRERE(NSTEPS), FILS(N), IPOOL(max(LIPOOL,1)), NE(NSTEPS),
     &        ND(NSTEPS), PROCNODE(NSTEPS), DAD(NSTEPS)
      DOUBLE PRECISION UU
      INTEGER  SLAVEF, KEEP(500), LOCAL_M, LOCAL_N
      INTEGER(8) KEEP8(150)
      INTEGER(8) ENTRIES_IN_FACTORS_LOC,
     &          ENTRIES_IN_FACTORS_LOC_MASTERS
      INTEGER  SBUF_SEND_FR, SBUF_REC_FR
      INTEGER  SBUF_SEND_LR, SBUF_REC_LR
      INTEGER(8) SBUF_RECOLD
      INTEGER  NMB_PAR2
      INTEGER  ISTEP_TO_INIV2( KEEP(71) )
      LOGICAL  I_AM_CAND(NMB_PAR2)
      INTEGER  CANDIDATES( SLAVEF+1, NMB_PAR2 )
      INTEGER  PHASE
      PARAMETER (PHASE=0)
      DOUBLE PRECISION OPSA
      DOUBLE PRECISION OPSA_LOC 
      INTEGER(8) MAX_SIZE_FACTOR
      DOUBLE PRECISION OPS_SUBTREE
      DOUBLE PRECISION OPS_SBTR_LOC 
      INTEGER, ALLOCATABLE, DIMENSION(:) :: TNSTK, LSTKI 
      INTEGER(8), ALLOCATABLE, DIMENSION(:) :: LSTKR
      INTEGER(8), ALLOCATABLE, DIMENSION(:) :: LSTKR_if_LRCB,
     &                                         LSTKRLR_CB_UD,
     &                                         LSTKRLR_CB_WC
      LOGICAL  OUTER_SENDS_FR  
      INTEGER(8) :: SAVE_SIZECB_UNDER_L0, 
     &              SAVE_SIZECB_UNDER_L0_IF_LRCB
      INTEGER SBUFR_FR, SBUFS_FR
      INTEGER SBUFR_LR, SBUFS_LR
      INTEGER(8) SBUFS_CB, SBUFR_CB
      INTEGER ITOP,NELIM,NFR
      INTEGER(8) ISTKR, LSTK
      INTEGER(8) :: NRLADU_CURRENT_MISSING
      INTEGER(8) :: NRLADU_CURRENT_K60_1
      LOGICAL    :: I_PROCESS_SCHUR_K60_1
      INTEGER(8) :: ISTKR_if_LRCB, ISTKRLR_CB_UD, ISTKRLR_CB_WC,
     &              K464_8, K465_8
      INTEGER    :: LRSTATUS, IDUMMY
      INTEGER    :: NBNODES_BLR
      LOGICAL    :: COMPRESS_PANEL, COMPRESS_CB
      LOGICAL    :: FORCE_BFR
      INTEGER ISTKI,  STKI, ISTKI_OOC
      INTEGER K,NSTK, IFATH
      INTEGER INODE, LEAF, IN
      INTEGER LEVEL, MAXITEMPCB
      INTEGER(8) CURRENT_ACTIVE_MEM, MAXTEMPCB
      LOGICAL UPDATE, UPDATEF, MASTER, MASTERF, INSSARBR
      INTEGER LEVELF, NCB, SIZECBI
      INTEGER(8) NCB8
      INTEGER(8) NFR8, NELIM8
      INTEGER(8) :: SIZEFRNOCBLU
      INTEGER    :: IDUMMY_ARRAY(1)
      INTEGER(8) SIZECB, SIZECBINFR, SIZECB_SLAVE
      INTEGER(8) SIZECB_if_LRCB, SIZECB_SLAVE_if_LRCB
      INTEGER(8) SIZECBLR_SLAVE_UD, SIZECBLR_SLAVE_WC
      INTEGER(8) SIZECBLR_UD, SIZECBLR_WC
      INTEGER(8) SIZECBSLR, NCBS8, 
     &           SIZECBS, SIZECBINFRS
      INTEGER NFRS, NELIMS, NCBS, LEVELS, LRSTATUSS
      LOGICAL COMPRESS_CBS
      INTEGER(8) :: PEAK_DYN_LRLU_UD, PEAK_DYN_LRCB_UD,
     &              PEAK_DYN_LRLUCB_UD, PEAK_DYN_LRLU_WC,
     &              PEAK_DYN_LRLUCB_WC
      INTEGER SIZEHEADER, SIZEHEADER_OOC, XSIZE_OOC
      INTEGER EXTRA_PERM_INFO_OOC
      INTEGER NBROWMAX, NSLAVES_LOC, NSLAVES_PASSED,
     &         NELIMF, NFRF, NCBF,
     &         NBROWMAXF, LKJIB_FR, LKJIB_LR,
     &         NBR, NBCOLFAC
      INTEGER(8) LEV3MAXREC, CBMAXR, CBMAXS
      INTEGER ALLOCOK
      INTEGER PANEL_SIZE
      LOGICAL PACKED_CB
      DOUBLE PRECISION OPS_NODE, OPS_NODE_MASTER, OPS_NODE_SLAVE
      INTEGER(8) ENTRIES_NODE_UPPER_PART, ENTRIES_NODE_LOWER_PART
      INTEGER :: FLAG_L0OMP
      PARAMETER (FLAG_L0OMP=-2014)
      INCLUDE 'mumps_headers.h'
      LOGICAL ROOT_OWNER
      INTEGER(8) LWK_RR
      INTEGER LIWK_RR
      INTEGER IROOT, SIZE_ROOT
      INTEGER WHAT
      INTEGER(8) IDUMMY8
      INTRINSIC min, int
      INTRINSIC dble
      INTEGER DMUMPS_OOC_GET_PANEL_SIZE
      EXTERNAL DMUMPS_OOC_GET_PANEL_SIZE
      INTEGER MUMPS_PROCNODE, MUMPS_TYPENODE
      LOGICAL MUMPS_IN_OR_ROOT_SSARBR
      EXTERNAL MUMPS_MAX_SURFCB_NBROWS
      EXTERNAL MUMPS_PROCNODE, MUMPS_TYPENODE, 
     &         MUMPS_IN_OR_ROOT_SSARBR
      logical :: FORCE_CAND, CONCERNED, UPDATES, STACKCB, MASTERSON
      integer :: IFSON, LEVELSON
      IDUMMY_ARRAY(1) = 0
      IF (KEEP(50).eq.2) THEN
        EXTRA_PERM_INFO_OOC = 1
      ELSE IF (KEEP(50).eq.0) THEN
        EXTRA_PERM_INFO_OOC = 2
      ELSE
        EXTRA_PERM_INFO_OOC = 0
      ENDIF
      PACKED_CB=( KEEP(215).EQ.0 .AND. KEEP(50).NE.0 )
      MAX_FRONT_SURFACE_LOCAL=0_8
      MAX_SIZE_FACTOR=0_8
      ALLOCATE( LSTKR(NSTEPS), TNSTK(NSTEPS), 
     &          LSTKI(NSTEPS) , 
     &          LSTKR_if_LRCB(NSTEPS), LSTKRLR_CB_UD(NSTEPS),
     &          LSTKRLR_CB_WC(NSTEPS),
     &          stat=ALLOCOK)
      if (ALLOCOK .GT. 0) THEN
        IFLAG  =-7
        IERROR = 6*NSTEPS
        RETURN
      endif
      LKJIB_FR = max(KEEP(5),KEEP(6))
      OUTER_SENDS_FR = (KEEP(263).NE.0 .OR. 
     &      KEEP(50).EQ.0. AND. (KEEP(468).LT.3 .OR. UU.EQ.0.0D0))
      IF ( OUTER_SENDS_FR ) THEN
         LKJIB_FR = max(LKJIB_FR, KEEP(420))
      ENDIF
      LKJIB_LR = max(LKJIB_FR,KEEP(142))
      IF (KEEP(198).NE.0.AND.SLAVEF.GT.1) THEN
         LKJIB_FR = min(LKJIB_FR*KEEP(179), KEEP(435))
      ENDIF
      TNSTK = NE
      LEAF = LIPOOL+1 
#if defined(OLD_OOC_NOPANEL)
      XSIZE_OOC=XSIZE_OOC_NOPANEL
#else
      IF (KEEP(50).EQ.0) THEN
              XSIZE_OOC=XSIZE_OOC_UNSYM
      ELSE
              XSIZE_OOC=XSIZE_OOC_SYM
      ENDIF
#endif
      SIZEHEADER_OOC = XSIZE_OOC+6 
      SIZEHEADER = XSIZE_IC + 6  
      ISTKR_if_LRCB   = 0_8
      ISTKRLR_CB_UD   = 0_8
      ISTKRLR_CB_WC   = 0_8
      ISTKR           = 0_8
      ISTKI           = 0
      ISTKI_OOC       = 0
      NBNODES_BLR     = 0
      OPSA_LOC   = 0.0D0
      ENTRIES_IN_FACTORS_LOC = 0_8
      ENTRIES_IN_FACTORS_LOC_MASTERS = 0_8
      OPS_SBTR_LOC = 0.0D0
      NRLADU     = 0_8
      NIRADU     = 0
      NIRADU_OOC = 0
      NRLADU_CURRENT_K60_1  = 0_8
      I_PROCESS_SCHUR_K60_1 = .FALSE.
      NRLADU_CURRENT = 0_8
      NRLADULR_UD = 0_8
      NRLADULR_WC = 0_8
      NRLADU_ROOT_3 = 0_8
      NRLNEC_ACTIVE = 0_8
      IF (ABOVE_L0) THEN
        SAVE_SIZECB_UNDER_L0 = SIZECB_UNDER_L0
        SAVE_SIZECB_UNDER_L0_IF_LRCB = SIZECB_UNDER_L0_IF_LRCB
      ELSE
        SAVE_SIZECB_UNDER_L0 = 0_8
        SAVE_SIZECB_UNDER_L0_IF_LRCB = 0_8
      ENDIF
      PEAK_DYN_LRLU_UD   = SAVE_SIZECB_UNDER_L0
      PEAK_DYN_LRCB_UD   = SAVE_SIZECB_UNDER_L0_IF_LRCB
      PEAK_DYN_LRLUCB_UD = SAVE_SIZECB_UNDER_L0_IF_LRCB
      PEAK_DYN_LRLU_WC   = SAVE_SIZECB_UNDER_L0
      PEAK_DYN_LRLUCB_WC = SAVE_SIZECB_UNDER_L0
      NRLNEC     = 0_8
      NRLADU_if_LR_LU      = 0_8
      NRLNEC_if_LR_LU      = 0_8
      NRLNEC_if_LR_CB      = 0_8
      NRLNEC_if_LR_LUCB    = 0_8
      NRLNECOOC_if_LR_LUCB = 0_8
      NRLNECLR_CB_UD       = 0_8
      NRLNECLR_LUCB_UD     = 0_8
      NRLNECLR_LUCB_WC     = 0_8
      NIRNEC     = 0
      NIRNEC_OOC = 0
      MAXFR      = 0
      PEAK_FR           = 0_8
      PEAK_FR_OOC       = 0_8
      PEAK_LRLU_UD      = 0_8
      PEAK_OOC_LRLU_UD  = 0_8
      PEAK_OOC_LRLU_WC  = 0_8
      PEAK_LRLUCB_UD    = 0_8
      PEAK_LRLUCB_WC    = 0_8
      PEAK_OOC_LRLUCB_UD= 0_8
      PEAK_OOC_LRLUCB_WC= 0_8
      PEAK_LRCB_UD      = 0_8
      PEAK_OOC_LRCB_UD  = 0_8
      ITOP       = 0
      MAXTEMPCB  = 0_8
      MAXITEMPCB = 0
      SBUFS_CB   = 1_8
      SBUFS_FR   = 1
      SBUFS_LR   = 1
      SBUFR_CB   = 1_8
      SBUFR_FR   = 1
      SBUFR_LR   = 1
      IF (KEEP(38) .NE. 0 .AND. KEEP(60).EQ.0) THEN
        INODE  = KEEP(38)
        NRLADU_ROOT_3 = int(LOCAL_M,8)*int(LOCAL_N,8)
        NRLADU = NRLADU_ROOT_3
        NRLNEC_ACTIVE = NRLADU_CURRENT
        MAX_SIZE_FACTOR=max(MAX_SIZE_FACTOR,NRLADU_ROOT_3)
        NRLNEC = NRLADU
        NRLADU_if_LR_LU      = NRLADU_ROOT_3
        NRLNECOOC_if_LR_LUCB = NRLNEC_ACTIVE
        NRLNEC_if_LR_LU   = NRLADU
        NRLNEC_if_LR_CB   = NRLADU
        NRLNEC_if_LR_LUCB = NRLADU
        PEAK_LRLU_UD      = max(PEAK_LRLU_UD,
     &             NRLNEC_if_LR_LU + NRLADULR_UD + SIZECB_UNDER_L0)
        PEAK_OOC_LRLU_UD = 
     &             max(PEAK_OOC_LRLU_UD,
     &             NRLNEC_ACTIVE + NRLADULR_UD)
        PEAK_OOC_LRLU_WC = 
     &             max(PEAK_OOC_LRLU_WC,
     &             NRLNEC_ACTIVE + NRLADULR_WC)
        PEAK_LRLUCB_UD = 
     &             max(PEAK_LRLUCB_UD,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_UD)
        PEAK_LRLUCB_WC = 
     &             max(PEAK_LRLUCB_WC,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_WC)
        PEAK_OOC_LRLUCB_UD = 
     &             max(PEAK_OOC_LRLUCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_UD)
        PEAK_OOC_LRLUCB_WC = 
     &             max(PEAK_OOC_LRLUCB_WC,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_WC)
        PEAK_LRCB_UD =
     &             max(PEAK_LRCB_UD,
     &             NRLNEC_if_LR_CB + NRLNECLR_CB_UD)
        PEAK_OOC_LRCB_UD =
     &             max(PEAK_OOC_LRCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_CB_UD)
        IF (ABOVE_L0) THEN
         PEAK_DYN_LRLU_UD   = max( PEAK_DYN_LRLU_UD,
     &        SIZECB_UNDER_L0 + NRLADULR_UD )
         PEAK_DYN_LRCB_UD   = max( PEAK_DYN_LRCB_UD,
     &        SIZECB_UNDER_L0_IF_LRCB + ISTKR_if_LRCB)
         PEAK_DYN_LRLUCB_UD = max( PEAK_DYN_LRLUCB_UD,
     &        SIZECB_UNDER_L0_IF_LRCB + ISTKR_if_LRCB + NRLADULR_UD)
         PEAK_DYN_LRLU_WC   = max( PEAK_DYN_LRLU_WC,
     &        SIZECB_UNDER_L0 + NRLADULR_WC)
         PEAK_DYN_LRLUCB_WC = max( PEAK_DYN_LRLUCB_WC,
     &        SIZECB_UNDER_L0 + ISTKRLR_CB_WC + NRLADULR_WC)
        ENDIF
        IF (MUMPS_PROCNODE(PROCNODE(STEP(INODE)),KEEP(199))
     &                                       .EQ. MYID) THEN
          NIRADU     = SIZEHEADER+2*(ND(STEP(INODE))+KEEP(253))
          NIRADU_OOC = SIZEHEADER_OOC+2*(ND(STEP(INODE))+KEEP(253))
        ELSE
          NIRADU     = SIZEHEADER
          NIRADU_OOC = SIZEHEADER_OOC
        ENDIF
        NIRNEC     = NIRADU
        NIRNEC_OOC = NIRADU_OOC
      ENDIF
      IF((KEEP(24).eq.0).OR.(KEEP(24).eq.1)) THEN
         FORCE_CAND=.FALSE.           
      ELSE
         FORCE_CAND=(mod(KEEP(24),2).eq.0)
      END IF
 90   CONTINUE
      IF (LEAF.NE.1) THEN
         LEAF = LEAF - 1
         INODE = IPOOL(LEAF)
      ELSE
         IF (LIPOOL.NE.0) THEN
           WRITE(MYID+6,*) ' ERROR 1 in DMUMPS_ANA_DISTM '
           CALL MUMPS_ABORT()
         ELSE
           GOTO 115
         ENDIF
      ENDIF
 95   CONTINUE 
      NFR    = ND(STEP(INODE))+KEEP(253)
      NFR8   = int(NFR,8)
      NSTK   = NE(STEP(INODE))
      NELIM = 0 
        IN = INODE
 100    NELIM = NELIM + 1 
      NELIM8=int(NELIM,8)
        IN = FILS(IN)
        IF (IN .GT. 0 ) GOTO 100
      IFSON = -IN
      IFATH = DAD(STEP(INODE))
      MASTER = MUMPS_PROCNODE(PROCNODE(STEP(INODE)),KEEP(199))
     &           .EQ. MYID
      LEVEL  = MUMPS_TYPENODE(PROCNODE(STEP(INODE)),KEEP(199))
      INSSARBR = MUMPS_IN_OR_ROOT_SSARBR(PROCNODE(STEP(INODE)),
     &        KEEP(199))
      UPDATE=.FALSE.
       if(.NOT.FORCE_CAND) then
         UPDATE = ( (MASTER.AND.(LEVEL.NE.3) ).OR. LEVEL.EQ.2 )
       else
         if(MASTER.and.(LEVEL.ne.3)) then
            UPDATE = .TRUE.
         else if(LEVEL.eq.2) then
            if ( I_AM_CAND(ISTEP_TO_INIV2(STEP(INODE)))) THEN
              UPDATE = .TRUE.
            end if
         end if
       end if
      NCB      = NFR-NELIM
      NCB8     = int(NCB,8)
      SIZECBINFR = NCB8*NCB8
      IF (KEEP(50).EQ.0) THEN
        SIZECB = SIZECBINFR
      ELSE
        IFATH = DAD(STEP(INODE))
        IF ( IFATH.NE.KEEP(38) .AND. PACKED_CB ) THEN
          SIZECB    = (NCB8*(NCB8+1_8))/2_8
        ELSE
          SIZECB    = SIZECBINFR
        ENDIF
      ENDIF
      IDUMMY = -99999
      CALL IS_FRONT_BLR_CANDIDATE (INODE, LEVEL, NFR, NELIM, 
     &                KEEP(494), KEEP(473), 1, KEEP(490),
     &                KEEP(491), KEEP(492),
     &                KEEP(20), KEEP(60), DAD(STEP(INODE)), KEEP(38),
     &                KEEP(123), LRSTATUS, IDUMMY)
      FORCE_BFR      = MUMPS_LR_FORCE_BFR(LRSTATUS)
      COMPRESS_PANEL = MUMPS_LR_COMPRESS_PANEL(LRSTATUS)
      COMPRESS_CB    = MUMPS_LR_COMPRESS_CB(LRSTATUS)
      IF (FORCE_BFR.OR.COMPRESS_PANEL.OR.COMPRESS_CB) 
     &  NBNODES_BLR = NBNODES_BLR+1
      IF (COMPRESS_PANEL) THEN
        K464_8 = int(KEEP(464),8)
      ELSE
        K464_8 = 1000_8
      ENDIF
      IF (COMPRESS_CB) THEN
        K465_8 = int(KEEP(465),8)
        SIZECB_if_LRCB = 0_8
        SIZECBLR_UD    = SIZECB*K465_8/1000_8
        SIZECBLR_WC    = SIZECB
      ELSE
        K465_8 = 1000_8
        SIZECBLR_UD    = 0_8
        SIZECBLR_WC    = 0_8
        SIZECB_if_LRCB = SIZECB 
      ENDIF
      SIZECBI      = 2* NCB  + SIZEHEADER
      IF (LEVEL.NE.2) THEN
        NSLAVES_LOC     = -99999999
        SIZECB_SLAVE = -99999997_8
        SIZECB_SLAVE_if_LRCB = SIZECB_SLAVE
        NBROWMAX        = NCB
      ELSE
        IF (KEEP(48) .EQ. 5) THEN
          WHAT = 5 
          IF (FORCE_CAND) THEN
            NSLAVES_LOC=CANDIDATES(SLAVEF+1,
     &                    ISTEP_TO_INIV2(STEP(INODE)))
          ELSE
            NSLAVES_LOC=SLAVEF-1
          ENDIF
          NSLAVES_PASSED=NSLAVES_LOC
        ELSE
          WHAT = 2 
          NSLAVES_PASSED=SLAVEF
          NSLAVES_LOC   =SLAVEF-1
        ENDIF
         CALL MUMPS_MAX_SURFCB_NBROWS(WHAT, KEEP,KEEP8,
     &     NCB, NFR, NSLAVES_PASSED, NBROWMAX, SIZECB_SLAVE
     &    )
       IF (COMPRESS_CB) THEN
        SIZECB_SLAVE_if_LRCB =  0_8
        SIZECBLR_SLAVE_UD = SIZECB_SLAVE*K465_8/1000_8
        SIZECBLR_SLAVE_WC = SIZECB_SLAVE
       ELSE
        SIZECB_SLAVE_if_LRCB = SIZECB_SLAVE
        SIZECBLR_SLAVE_UD    = 0_8
        SIZECBLR_SLAVE_WC    = 0_8
       ENDIF
      ENDIF
      IF (KEEP(60).GT.1) THEN
         IF (MASTER .AND. INODE.EQ.KEEP(38)) THEN
          NIRADU     = NIRADU+SIZEHEADER+2*(ND(STEP(INODE))+KEEP(253))
          NIRADU_OOC = NIRADU_OOC+SIZEHEADER_OOC+
     &                 2*(ND(STEP(INODE))+KEEP(253))
         ENDIF
      ENDIF
      IF (LEVEL.EQ.3) THEN
         IF ( 
     &     KEEP(60).LE.1
     &      ) THEN
           NRLADU_CURRENT = int(LOCAL_M,8)*int(LOCAL_N,8)
           NRLNEC = max(NRLNEC,NRLADU+ISTKR+
     &                 NRLADU_CURRENT)
           NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_ROOT_3 + 
     &                        NRLADU_CURRENT+ISTKR)
           NRLNEC_if_LR_LU   = max(NRLNEC_if_LR_LU,
     &                         NRLADU_if_LR_LU+ISTKR+
     &                         NRLADU_CURRENT)
           NRLNEC_if_LR_CB   = max(NRLNEC_if_LR_CB  ,
     &                         NRLADU+ISTKR_if_LRCB+
     &                         NRLADU_CURRENT)
           NRLNEC_if_LR_LUCB = max(NRLNEC_if_LR_LUCB,
     &                         NRLADU_if_LR_LU+ISTKR_if_LRCB+
     &                         NRLADU_CURRENT)
           NRLNECOOC_if_LR_LUCB = max(NRLNECOOC_if_LR_LUCB,
     &                            NRLADU_ROOT_3 + 
     &                            NRLADU_CURRENT+ISTKR_if_LRCB)
           PEAK_LRLU_UD     = max(PEAK_LRLU_UD,
     &                        NRLNEC_if_LR_LU + NRLADULR_UD)
           PEAK_OOC_LRLU_UD = 
     &             max(PEAK_OOC_LRLU_UD,
     &             NRLNEC_ACTIVE + NRLADULR_UD)
           PEAK_OOC_LRLU_WC = 
     &             max(PEAK_OOC_LRLU_WC,
     &             NRLNEC_ACTIVE + NRLADULR_WC)
           PEAK_LRLUCB_UD = 
     &             max(PEAK_LRLUCB_UD,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_UD)
           PEAK_LRLUCB_WC = 
     &             max(PEAK_LRLUCB_WC,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_WC)
           PEAK_OOC_LRLUCB_UD = 
     &             max(PEAK_OOC_LRLUCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_UD)
           PEAK_OOC_LRLUCB_WC = 
     &             max(PEAK_OOC_LRLUCB_WC,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_WC)
           PEAK_LRCB_UD =
     &             max(PEAK_LRCB_UD,
     &             NRLNEC_if_LR_CB + NRLNECLR_CB_UD)
           PEAK_OOC_LRCB_UD =
     &             max(PEAK_OOC_LRCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_CB_UD)
           IF (ABOVE_L0) THEN
            PEAK_DYN_LRLU_UD   = max( PEAK_DYN_LRLU_UD,
     &           SIZECB_UNDER_L0 + NRLADULR_UD )
            PEAK_DYN_LRCB_UD   = max( PEAK_DYN_LRCB_UD,
     &           SIZECB_UNDER_L0_IF_LRCB + ISTKR_if_LRCB)
            PEAK_DYN_LRLUCB_UD = max( PEAK_DYN_LRLUCB_UD,
     &           SIZECB_UNDER_L0_IF_LRCB + ISTKR_if_LRCB + NRLADULR_UD)
            PEAK_DYN_LRLU_WC   = max( PEAK_DYN_LRLU_WC,
     &           SIZECB_UNDER_L0 + NRLADULR_WC)
            PEAK_DYN_LRLUCB_WC = max( PEAK_DYN_LRLUCB_WC,
     &           SIZECB_UNDER_L0 + ISTKRLR_CB_WC + NRLADULR_WC)
           ENDIF
         ENDIF
         IF (MASTER) THEN 
            IF (NFR.GT.MAXFR) MAXFR = NFR
         ENDIF
      ENDIF
      IF(KEEP(86).EQ.1)THEN
         IF(MASTER.AND.(.NOT.MUMPS_IN_OR_ROOT_SSARBR(
     &        PROCNODE(STEP(INODE)), KEEP(199)))
     &     )THEN
            IF(LEVEL.EQ.1)THEN
               MAX_FRONT_SURFACE_LOCAL=max(MAX_FRONT_SURFACE_LOCAL,
     &              NFR8*NFR8)
               IF (KEEP(268).NE.0) THEN
                 MAX_FRONT_SURFACE_LOCAL=max(MAX_FRONT_SURFACE_LOCAL,
     &              NFR8*NFR8+NELIM8)
               ENDIF
            ELSEIF(LEVEL.EQ.2)THEN
               IF(KEEP(50).EQ.0)THEN
                 MAX_FRONT_SURFACE_LOCAL=max(MAX_FRONT_SURFACE_LOCAL,
     &                 NFR8*NELIM8)
               ELSE
                 MAX_FRONT_SURFACE_LOCAL=max(MAX_FRONT_SURFACE_LOCAL,
     &                 NELIM8*NELIM8)
                 IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
                  MAX_FRONT_SURFACE_LOCAL=max(MAX_FRONT_SURFACE_LOCAL,
     &                  NELIM8*(NELIM8+1_8))
                 ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      IF (LEVEL.EQ.2) THEN
        IF (MASTER) THEN
          IF (KEEP(50).EQ.0) THEN
             SBUFS_FR = max(SBUFS_FR, NFR*LKJIB_FR+LKJIB_FR+4)
             SBUFS_LR = max(SBUFS_LR, NFR*LKJIB_LR+LKJIB_LR+4)
          ELSE
             SBUFS_FR = max(SBUFS_FR, NELIM*LKJIB_FR+NELIM+6)
             SBUFS_LR = max(SBUFS_LR, NELIM*LKJIB_LR+NELIM+6)
          ENDIF
        ELSEIF (UPDATE) THEN
            if (KEEP(50).EQ.0) THEN
              SBUFR_FR   = max(SBUFR_FR, NFR*LKJIB_FR+LKJIB_FR+4)
              SBUFR_LR   = max(SBUFR_LR, NFR*LKJIB_LR+LKJIB_LR+4)
            else
              SBUFR_FR = max( SBUFR_FR, NELIM*LKJIB_FR+NELIM+6 )
              SBUFR_LR = max( SBUFR_LR, NELIM*LKJIB_LR+NELIM+6 )
              SBUFS_FR = max( SBUFS_FR, NBROWMAX*LKJIB_FR+6 )
              SBUFS_LR = max( SBUFS_LR, NBROWMAX*LKJIB_LR+6 )
              SBUFR_FR = max( SBUFR_FR, NBROWMAX*LKJIB_FR+6 )
              SBUFR_LR = max( SBUFR_LR, NBROWMAX*LKJIB_LR+6 )
            endif
        ENDIF
      ENDIF
      IF ( UPDATE ) THEN
          IF ( (MASTER) .AND. (LEVEL.EQ.1) ) THEN
            NIRADU     = NIRADU + 2*NFR + SIZEHEADER
            NIRADU_OOC = NIRADU_OOC + 2*NFR + SIZEHEADER_OOC
            PANEL_SIZE = DMUMPS_OOC_GET_PANEL_SIZE(
     &      2_8*int(KEEP(226),8), NFR, KEEP(227), KEEP(50))
            NIRADU_OOC = NIRADU_OOC +
     &      EXTRA_PERM_INFO_OOC*(2+NELIM + NELIM/PANEL_SIZE+1)
            IF (KEEP(50).EQ.0) THEN
              IF (KEEP(251) .EQ. 2) THEN
                NRLADU_CURRENT = int(NELIM,8) * int(NFR,8)
                SIZEFRNOCBLU   = int(NFR-NELIM,8)*int(NELIM)
              ELSE
                NRLADU_CURRENT = int(NELIM,8)*int(2*NFR-NELIM,8)
                SIZEFRNOCBLU   = 0_8
              ENDIF
              NRLADU = NRLADU + NRLADU_CURRENT
              MAX_SIZE_FACTOR  = max(MAX_SIZE_FACTOR,NRLADU_CURRENT)
            ELSE
              IF ( KEEP(459) .GT. 1 .AND. KEEP(50).NE.0 ) THEN
                CALL MUMPS_LDLTPANEL_STORAGE( NELIM, KEEP, 
     &          IDUMMY_ARRAY(1), NRLADU_CURRENT )
                SIZEFRNOCBLU   = int(NELIM,8)*int(NELIM,8) -
     &                           NRLADU_CURRENT
                NRLADU_CURRENT = NRLADU_CURRENT +
     &                           int(NELIM,8) * int(NFR-NELIM,8)
              ELSE
                SIZEFRNOCBLU   = 0_8
                NRLADU_CURRENT = int(NELIM,8)*int(NFR,8)
              ENDIF
              NRLADU = NRLADU + NRLADU_CURRENT
              MAX_SIZE_FACTOR=max(MAX_SIZE_FACTOR,NRLADU_CURRENT)
             ENDIF
             IF (INODE.EQ.KEEP(20).AND.(KEEP(60).EQ.1)) THEN
                I_PROCESS_SCHUR_K60_1 = .TRUE.
                NRLADU_CURRENT_K60_1  = NRLADU_CURRENT
             ENDIF
             IF (FORCE_BFR.OR.COMPRESS_PANEL) THEN
               NRLADU_if_LR_LU = NRLADU_if_LR_LU + 0_8
               NRLADU_CURRENT_MISSING = NRLADU_CURRENT
               NRLADULR_UD     = NRLADULR_UD + 
     &                           NRLADU_CURRENT*K464_8/1000_8
               NRLADULR_WC     = NRLADULR_WC + 
     &                           NRLADU_CURRENT
             ELSE
               NRLADU_CURRENT_MISSING = 0_8
               NRLADU_if_LR_LU = NRLADU_if_LR_LU + NRLADU_CURRENT
             ENDIF
            SIZECBI        = 2* NCB  + SIZEHEADER
          ELSEIF (LEVEL.EQ.2) THEN
            IF (MASTER) THEN
              NIRADU     = NIRADU+SIZEHEADER +SLAVEF-1+2*NFR 
              NIRADU_OOC = NIRADU_OOC+SIZEHEADER_OOC +SLAVEF-1+2*NFR 
              IF (KEEP(50).EQ.0) THEN
                NRLADU_CURRENT = int(NFR,8)*int(NELIM,8)
                SIZEFRNOCBLU   = 0_8
                NBCOLFAC       = NFR
              ELSE
                NBCOLFAC       = NELIM
                IF ( KEEP(459) .GT. 1 .AND. KEEP(50) .NE. 0 ) THEN
                  CALL MUMPS_LDLTPANEL_STORAGE( NELIM, KEEP, 
     &            IDUMMY_ARRAY(1), NRLADU_CURRENT )
                  SIZEFRNOCBLU = int(NELIM,8)*int(NELIM,8) -
     &                           NRLADU_CURRENT
                ELSE
                  SIZEFRNOCBLU   = 0_8
                  NRLADU_CURRENT = int(NELIM,8)*int(NELIM,8)
                ENDIF
              ENDIF
              PANEL_SIZE = DMUMPS_OOC_GET_PANEL_SIZE(
     &        2_8*int(KEEP(226),8), NBCOLFAC, KEEP(227), KEEP(50))
              NIRADU_OOC = NIRADU_OOC +
     &        EXTRA_PERM_INFO_OOC*(2+NELIM + NELIM/PANEL_SIZE+1)
              NRLADU = NRLADU + NRLADU_CURRENT
              IF (FORCE_BFR.OR.COMPRESS_PANEL) THEN
               NRLADU_if_LR_LU = NRLADU_if_LR_LU + 0_8
               NRLADU_CURRENT_MISSING = NRLADU_CURRENT
               NRLADULR_UD     = NRLADULR_UD + 
     &                           NRLADU_CURRENT*K464_8/1000_8
               NRLADULR_WC     = NRLADULR_WC + 
     &                           NRLADU_CURRENT
              ELSE
               NRLADU_CURRENT_MISSING = 0_8
               NRLADU_if_LR_LU = NRLADU_if_LR_LU + NRLADU_CURRENT
              ENDIF
              MAX_SIZE_FACTOR=max(MAX_SIZE_FACTOR,NRLADU_CURRENT)
               SIZECB     = 0_8
               SIZECBINFR = 0_8
               SIZECB_if_LRCB = 0_8
               SIZECBLR_UD    = 0_8
               SIZECBLR_WC    = 0_8
               SIZECBI    = NCB + SIZEHEADER + SLAVEF - 1
            ELSE
             SIZECB=SIZECB_SLAVE
             SIZECBINFR = SIZECB
             NIRADU       = NIRADU+4+NELIM+NBROWMAX
             NIRADU_OOC   = NIRADU_OOC+4+NELIM+NBROWMAX
             IF (KEEP(50).EQ.0) THEN
               IF (KEEP(251) .EQ. 2) THEN
                 NRLADU_CURRENT = 0_8
                 SIZEFRNOCBLU   = int(NELIM,8)*int(NBROWMAX,8)
               ELSE
                 NRLADU_CURRENT = int(NELIM,8)*int(NBROWMAX,8)
                 SIZEFRNOCBLU   = 0_8
               ENDIF
             ELSE
               NRLADU_CURRENT   = int(NELIM,8)*int(NCB/NSLAVES_LOC,8)
               SIZEFRNOCBLU     = 0_8
             ENDIF
             NRLADU   = NRLADU + NRLADU_CURRENT
             MAX_SIZE_FACTOR=max(MAX_SIZE_FACTOR,NRLADU_CURRENT)
             IF (FORCE_BFR.OR.COMPRESS_PANEL) THEN
               NRLADU_if_LR_LU = NRLADU_if_LR_LU + 0_8
               NRLADU_CURRENT_MISSING = NRLADU_CURRENT
               NRLADULR_UD     = NRLADULR_UD + 
     &                           NRLADU_CURRENT*K464_8/1000_8
               NRLADULR_WC     = NRLADULR_WC + 
     &                           NRLADU_CURRENT
             ELSE
               NRLADU_CURRENT_MISSING = 0_8
               NRLADU_if_LR_LU = NRLADU_if_LR_LU + NRLADU_CURRENT
             ENDIF
             IF (KEEP(50).EQ.0) THEN
               SIZECBI       = 7 + NBROWMAX + NCB
             ELSE
               SIZECBI       = 8 + NBROWMAX + NCB
             ENDIF
             IF (KEEP(50).NE.0) THEN 
                     SIZECBI=SIZECBI+NSLAVES_LOC+
     &                                  XTRA_SLAVES_SYM
             ELSE
                     SIZECBI=SIZECBI+NSLAVES_LOC+
     &                                  XTRA_SLAVES_UNSYM 
             ENDIF
            ENDIF
         ENDIF
         NIRNEC = max(NIRNEC,
     &             NIRADU+ISTKI+SIZECBI+MAXITEMPCB)
         NIRNEC_OOC = max(NIRNEC_OOC,
     &             NIRADU_OOC+ISTKI_OOC+SIZECBI+MAXITEMPCB +
     &             (XSIZE_OOC-XSIZE_IC) ) 
         CURRENT_ACTIVE_MEM = ISTKR+SIZECBINFR
         IF ( KEEP(50).NE.0 .AND. LEVEL.EQ.1 ) THEN
             SIZEFRNOCBLU = SIZEFRNOCBLU + int(NELIM,8)*int(NCB,8)
         ENDIF
         CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM + SIZEFRNOCBLU
         IF (MASTER .AND.  KEEP(219).NE.0.AND.
     &       KEEP(50).EQ.2.AND.LEVEL.EQ.2) THEN
             CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM + int(NELIM,8)
         ENDIF
         IF (SLAVEF.EQ.1) THEN
           NRLNEC_if_LR_CB   = 
     &             max(NRLNEC_if_LR_CB  ,NRLADU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB)
           NRLNEC_if_LR_LUCB = 
     &             max(NRLNEC_if_LR_LUCB,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB+
     &             NRLADU_CURRENT_MISSING)
           NRLNECOOC_if_LR_LUCB = 
     &             max(NRLNECOOC_if_LR_LUCB,
     &             NRLADU_CURRENT+NRLADU_ROOT_3+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB)
         ELSE
           NRLNEC_if_LR_CB   = 
     &             max(NRLNEC_if_LR_CB  ,NRLADU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB+
     &             MAXTEMPCB)
           NRLNEC_if_LR_LUCB = 
     &             max(NRLNEC_if_LR_LUCB,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB+
     &             NRLADU_CURRENT_MISSING+
     &             MAXTEMPCB)
           NRLNECOOC_if_LR_LUCB = max(NRLNECOOC_if_LR_LUCB,
     &             NRLADU_CURRENT+NRLADU_ROOT_3+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB+
     &             MAXTEMPCB)
         ENDIF
         IF (NSTK .NE. 0 .AND. INSSARBR .AND.
     &     KEEP(234).NE.0 .AND. KEEP(55).EQ.0) THEN
           CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM - LSTKR(ITOP)
         ENDIF
         IF (SLAVEF.EQ.1) THEN
           NRLNEC = max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM)
           NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_CURRENT+
     &             NRLADU_ROOT_3+CURRENT_ACTIVE_MEM)
           NRLNEC_if_LR_LU = 
     &             max(NRLNEC_if_LR_LU,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM+NRLADU_CURRENT_MISSING)
         ELSE
           NRLNEC = max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM+MAXTEMPCB)
           NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_CURRENT+
     &             NRLADU_ROOT_3+CURRENT_ACTIVE_MEM+MAXTEMPCB)
           NRLNEC_if_LR_LU = 
     &             max(NRLNEC_if_LR_LU,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM+MAXTEMPCB+
     &             NRLADU_CURRENT_MISSING)
         ENDIF
         NRLNECLR_LUCB_UD =  max(NRLNECLR_LUCB_UD,
     &                           NRLADULR_UD+ISTKRLR_CB_UD)
         NRLNECLR_LUCB_WC =  max(NRLNECLR_LUCB_WC,
     &                           NRLADULR_WC+ISTKRLR_CB_WC)
         PEAK_LRLU_UD = 
     &             max(PEAK_LRLU_UD,
     &             NRLNEC_if_LR_LU + NRLADULR_UD)
         PEAK_OOC_LRLU_UD = 
     &             max(PEAK_OOC_LRLU_UD,
     &             NRLNEC_ACTIVE + NRLADULR_UD)
         PEAK_OOC_LRLU_WC = 
     &             max(PEAK_OOC_LRLU_WC,
     &             NRLNEC_ACTIVE + NRLADULR_WC)
         PEAK_LRLUCB_UD = 
     &             max(PEAK_LRLUCB_UD,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_UD)
         PEAK_LRLUCB_WC = 
     &             max(PEAK_LRLUCB_WC,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_WC)
         PEAK_OOC_LRLUCB_UD = 
     &             max(PEAK_OOC_LRLUCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_UD)
         PEAK_OOC_LRLUCB_WC = 
     &             max(PEAK_OOC_LRLUCB_WC,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_WC)
         PEAK_LRCB_UD =
     &             max(PEAK_LRCB_UD,
     &             NRLNEC_if_LR_CB + NRLNECLR_CB_UD)
         PEAK_OOC_LRCB_UD =
     &             max(PEAK_OOC_LRCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_CB_UD)
         IF (ABOVE_L0) THEN
          PEAK_DYN_LRLU_UD   = max( PEAK_DYN_LRLU_UD,
     &         SIZECB_UNDER_L0 + NRLADULR_UD )
          PEAK_DYN_LRCB_UD   = max( PEAK_DYN_LRCB_UD,
     &         SIZECB_UNDER_L0_IF_LRCB + ISTKR_if_LRCB)
          PEAK_DYN_LRLUCB_UD = max( PEAK_DYN_LRLUCB_UD,
     &         SIZECB_UNDER_L0_IF_LRCB + ISTKR_if_LRCB + NRLADULR_UD)
          PEAK_DYN_LRLU_WC   = max( PEAK_DYN_LRLU_WC,
     &         SIZECB_UNDER_L0 + NRLADULR_WC)
          PEAK_DYN_LRLUCB_WC = max( PEAK_DYN_LRLUCB_WC,
     &         SIZECB_UNDER_L0 + ISTKRLR_CB_WC + NRLADULR_WC)
         ENDIF
         IF (NFR.GT.MAXFR) MAXFR = NFR
         IF (ABOVE_L0.AND.MASTER.AND.(LEVEL.EQ.1))  THEN
           DO WHILE (IFSON.GT.0) 
             IF (TNSTK(STEP(IFSON)).EQ.FLAG_L0OMP) THEN
               LEVELS = MUMPS_TYPENODE(PROCNODE(STEP(IFSON)),KEEP(199))
               NFRS   = ND(STEP(IFSON))+KEEP(253)
               NELIMS= 0 
               IN = IFSON
               DO WHILE (IN.GT.0)
                IN     = FILS(IN)
                NELIMS = NELIMS + 1
               ENDDO
               NCBS     = NFRS-NELIMS
               NCBS8    = int(NCBS,8)
               SIZECBINFRS = NCBS8*NCBS8
               IF (KEEP(50).EQ.0) THEN
                 SIZECBS = SIZECBINFRS
               ELSE
                 IF ( PACKED_CB ) THEN
                   SIZECBS   = (NCBS8*(NCBS8+1_8))/2_8
                 ELSE
                   SIZECBS   = SIZECBINFRS
                 ENDIF
               ENDIF
               IDUMMY = -99999
               CALL IS_FRONT_BLR_CANDIDATE 
     &               (IFSON, LEVELS, NFRS, NELIMS, 
     &                KEEP(494), KEEP(473), 1, KEEP(490),
     &                KEEP(491), KEEP(492),
     &                KEEP(20), KEEP(60), DAD(STEP(IFSON)), KEEP(38),
     &                KEEP(123), LRSTATUSS, IDUMMY)
               COMPRESS_CBS   = MUMPS_LR_COMPRESS_CB(LRSTATUSS)
               IF (COMPRESS_CBS) THEN
                 K465_8    = int(KEEP(465),8)
                 SIZECBSLR = SIZECBS*K465_8/1000_8
               ELSE
                 SIZECBSLR = SIZECBS
               ENDIF
               SIZECB_UNDER_L0 = SIZECB_UNDER_L0 - SIZECBS
               SIZECB_UNDER_L0_IF_LRCB = SIZECB_UNDER_L0_IF_LRCB 
     &                                   - SIZECBSLR
             ENDIF
             IFSON = FRERE(STEP(IFSON))
           ENDDO
         ENDIF
         IF (NSTK.GT.0) THEN
            DO 70 K=1,NSTK
               LSTK = LSTKR(ITOP)
               ISTKR = ISTKR - LSTK
               IF (K==1 .AND. INSSARBR.AND.KEEP(234).NE.0
     &            .AND.KEEP(55).EQ.0) THEN
               ELSE
                 CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM - LSTK
               ENDIF
               LSTK          = LSTKR_if_LRCB(ITOP)
               ISTKR_if_LRCB = ISTKR_if_LRCB - LSTK
               LSTK          = LSTKRLR_CB_UD(ITOP) 
               ISTKRLR_CB_UD = ISTKRLR_CB_UD - LSTK
               LSTK          = LSTKRLR_CB_WC(ITOP) 
               ISTKRLR_CB_WC = ISTKRLR_CB_WC - LSTK
               STKI = LSTKI( ITOP )
               ISTKI = ISTKI - STKI
               ISTKI_OOC = ISTKI_OOC - STKI - (XSIZE_OOC-XSIZE_IC)
               ITOP = ITOP - 1
               IF (ITOP.LT.0) THEN
                  write(*,*) MYID,
     &            ': ERROR 2 in DMUMPS_ANA_DISTM. ITOP = ',ITOP
                  CALL MUMPS_ABORT()
               ENDIF
 70         CONTINUE
         ENDIF 
      ELSE IF (LEVEL.NE.3) THEN
         DO WHILE (IFSON.GT.0) 
            UPDATES=.FALSE.
            MASTERSON = MUMPS_PROCNODE(PROCNODE(STEP(IFSON)),KEEP(199))
     &                  .EQ.MYID
            LEVELSON  = MUMPS_TYPENODE(PROCNODE(STEP(IFSON)),KEEP(199))
            if(.NOT.FORCE_CAND) then
               UPDATES =((MASTERSON.AND.(LEVELSON.NE.3)).OR. 
     &                   LEVELSON.EQ.2)
            else
               if(MASTERSON.and.(LEVELSON.ne.3)) then
                  UPDATES = .TRUE.
               else if(LEVELSON.eq.2) then
                  if ( I_AM_CAND(ISTEP_TO_INIV2(STEP(IFSON)))) then
                    UPDATES = .TRUE.
                  end if
               end if
            end if
            IF (UPDATES) THEN
              LSTK = LSTKR(ITOP)
              ISTKR = ISTKR - LSTK
              LSTK          = LSTKR_if_LRCB(ITOP)
              ISTKR_if_LRCB = ISTKR_if_LRCB - LSTK
              LSTK          =  LSTKRLR_CB_UD(ITOP) 
              ISTKRLR_CB_UD = ISTKRLR_CB_UD - LSTK
              LSTK          =  LSTKRLR_CB_WC(ITOP) 
              ISTKRLR_CB_WC = ISTKRLR_CB_WC - LSTK
              STKI = LSTKI( ITOP )
              ISTKI = ISTKI - STKI
              ISTKI_OOC = ISTKI_OOC - STKI - (XSIZE_OOC-XSIZE_IC)
              ITOP = ITOP - 1
              IF (ITOP.LT.0) THEN
                write(*,*) MYID,
     &          ': ERROR 2 in DMUMPS_ANA_DISTM. ITOP = ',ITOP
                CALL MUMPS_ABORT()
              ENDIF
            ENDIF
            IFSON = FRERE(STEP(IFSON)) 
         END DO
      ENDIF
      IF (
     &        ( (INODE.NE.KEEP(20)).OR.(KEEP(60).EQ.0) ) 
     &       .AND.
     &        ( (INODE.NE.KEEP(38)).OR.(KEEP(60).LE.1) ) 
     &      )
     &  THEN
            ENTRIES_NODE_LOWER_PART = int(NFR-NELIM,8) * int(NELIM,8)
            IF ( KEEP(50).EQ.0 ) THEN
              ENTRIES_NODE_UPPER_PART = int(NFR,8) * int(NELIM,8)
            ELSE
              ENTRIES_NODE_UPPER_PART =
     &        (int(NELIM,8)*int(NELIM+1,8))/2_8
            ENDIF
            IF (KEEP(50).EQ.2 .AND. LEVEL.EQ.3) THEN
              CALL MUMPS_GET_FLOPS_COST(NFR, 
     &           NELIM, NELIM, 0,
     &           1,OPS_NODE)
            ELSE
              CALL MUMPS_GET_FLOPS_COST(NFR, 
     &           NELIM, NELIM,KEEP(50),
     &           1,OPS_NODE)
            ENDIF
            IF (LEVEL.EQ.2) THEN
              CALL MUMPS_GET_FLOPS_COST(NFR, 
     &           NELIM, NELIM,KEEP(50),
     &           2,OPS_NODE_MASTER)
              OPS_NODE_SLAVE=OPS_NODE-OPS_NODE_MASTER
            ENDIF
      ELSE
           OPS_NODE = 0.0D0
           ENTRIES_NODE_UPPER_PART = 0_8
           ENTRIES_NODE_LOWER_PART = 0_8
      ENDIF
      IF ( MASTER )  THEN
        ENTRIES_IN_FACTORS_LOC_MASTERS = 
     &                     ENTRIES_IN_FACTORS_LOC_MASTERS +
     &                            ENTRIES_NODE_UPPER_PART +
     &                            ENTRIES_NODE_LOWER_PART
      ENDIF
      IF (UPDATE.OR.LEVEL.EQ.3) THEN
         IF ( LEVEL .EQ. 3 ) THEN
            IF (ROOT_yes) THEN
              OPSA_LOC = OPSA_LOC + 
     &                   dble( 
     &                        int(OPS_NODE,8)/
     &                        int(ROOT_NPROW*ROOT_NPCOL,8)
     &                       )
              ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC +
     &                            ENTRIES_NODE_UPPER_PART /
     &                            int(ROOT_NPROW*ROOT_NPCOL,8)
              IF (MASTER) THEN
                ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC +
     &                                 mod(ENTRIES_NODE_UPPER_PART,
     &                                 int(SLAVEF,8))
              ENDIF
            ENDIF
         ELSE IF (MASTER .AND. LEVEL.EQ.2) THEN
            OPSA_LOC = OPSA_LOC + OPS_NODE_MASTER
            ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC +
     &                      ENTRIES_NODE_UPPER_PART +
     &                      mod(ENTRIES_NODE_LOWER_PART,
     &                          int(NSLAVES_LOC,8))
         ELSE IF (MASTER .AND. LEVEL.EQ.1) THEN
            OPSA_LOC = OPSA_LOC + OPS_NODE
            ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC +
     &                               ENTRIES_NODE_UPPER_PART +
     &                               ENTRIES_NODE_LOWER_PART
         ELSE IF (UPDATE) THEN 
            OPSA_LOC = OPSA_LOC + 
     &            dble(OPS_NODE_SLAVE)/dble(NSLAVES_LOC)
            ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC
     &                 + ENTRIES_NODE_LOWER_PART / 
     &                 int(NSLAVES_LOC,8)
         ENDIF
         IF (MUMPS_IN_OR_ROOT_SSARBR(PROCNODE(STEP(INODE)),
     &        KEEP(199)) .OR. NE(STEP(INODE))==0) THEN
           IF (LEVEL == 1) THEN
             OPS_SBTR_LOC = OPS_SBTR_LOC + OPS_NODE
           ELSE
             CALL MUMPS_GET_FLOPS_COST(NFR, 
     &           NELIM, NELIM,KEEP(50),
     &           1,OPS_NODE)
             OPS_SBTR_LOC = OPS_SBTR_LOC + OPS_NODE
           ENDIF
         ENDIF
        ENDIF
      IF (IFATH .EQ. 0) THEN
        IF (LEAF.GT.1) THEN
         GOTO 90 
        ELSE
         GOTO 115 
        ENDIF 
      ELSE
         NFRF = ND(STEP(IFATH))+KEEP(253)
         IF (DAD(STEP(IFATH)).EQ.0) THEN
           NELIMF = NFRF-KEEP(253)
           IF (ABOVE_L0) IN=0
         ELSE
           NELIMF = 0
           IN = IFATH
           DO WHILE (IN.GT.0)
              IN = FILS(IN)
              NELIMF = NELIMF+1
           ENDDO
         ENDIF
         NCBF = NFRF - NELIMF
         LEVELF = MUMPS_TYPENODE(PROCNODE(STEP(IFATH)),KEEP(199))
         MASTERF= MUMPS_PROCNODE(PROCNODE(STEP(IFATH)),
     &                           KEEP(199)).EQ.MYID
         UPDATEF= .FALSE.
         if(.NOT.FORCE_CAND) then
            UPDATEF= ((MASTERF.AND.(LEVELF.NE.3)).OR.LEVELF.EQ.2)
         else
            if(MASTERF.and.(LEVELF.ne.3)) then
               UPDATEF = .TRUE.
            else if (LEVELF.eq.2) then
               if ( I_AM_CAND(ISTEP_TO_INIV2(STEP(IFATH)))) THEN
                 UPDATEF = .TRUE.
               end if
            end if
         end if
         CONCERNED  = UPDATEF .OR. UPDATE
         IF (LEVELF .NE. 2) THEN
           NBROWMAXF = -999999
         ELSE
           IF (KEEP(48) .EQ. 5) THEN
               WHAT = 4
               IF (FORCE_CAND) THEN
                 NSLAVES_LOC=CANDIDATES(SLAVEF+1,
     &               ISTEP_TO_INIV2(STEP(IFATH)))
               ELSE
                 NSLAVES_LOC=SLAVEF-1
               ENDIF
           ELSE
               WHAT = 1 
               NSLAVES_LOC=SLAVEF
           ENDIF
           CALL MUMPS_MAX_SURFCB_NBROWS( WHAT, KEEP, KEEP8,
     &     NCBF, NFRF, NSLAVES_LOC, NBROWMAXF, IDUMMY8 
     &          )
         ENDIF
         IF(LEVEL.EQ.1.AND.UPDATE.AND.
     &      (UPDATEF.OR.LEVELF.EQ.2)
     &      .AND.LEVELF.NE.3) THEN
             IF ( INSSARBR .AND. KEEP(234).NE.0) THEN
               NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_CURRENT+
     &           NRLADU_ROOT_3+CURRENT_ACTIVE_MEM)
               NRLNEC = max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM)
               NRLNEC_if_LR_LU = 
     &             max(NRLNEC_if_LR_LU,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM+
     &             NRLADU_CURRENT_MISSING)
               NRLNEC_if_LR_CB   = 
     &             max(NRLNEC_if_LR_CB  ,NRLADU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB)
               NRLNEC_if_LR_LUCB = 
     &             max(NRLNEC_if_LR_LUCB,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB+
     &             NRLADU_CURRENT_MISSING)
               NRLNECOOC_if_LR_LUCB = max(NRLNECOOC_if_LR_LUCB,
     &             NRLADU_CURRENT+NRLADU_ROOT_3+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB)
             ELSE
               NRLNEC = max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM+SIZECB)
               NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_CURRENT+
     &           NRLADU_ROOT_3+CURRENT_ACTIVE_MEM+SIZECB)
               NRLNEC_if_LR_LU = 
     &             max(NRLNEC_if_LR_LU,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM+SIZECB+
     &             NRLADU_CURRENT_MISSING)
               NRLNEC_if_LR_CB   = 
     &             max(NRLNEC_if_LR_CB  ,NRLADU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ 
     &             ISTKR_if_LRCB+ SIZECB)
               NRLNEC_if_LR_LUCB = 
     &             max(NRLNEC_if_LR_LUCB,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ 
     &             ISTKR_if_LRCB+ SIZECB+
     &             NRLADU_CURRENT_MISSING)
               NRLNECOOC_if_LR_LUCB = max(NRLNECOOC_if_LR_LUCB,
     &             NRLADU_CURRENT+NRLADU_ROOT_3+
     &             CURRENT_ACTIVE_MEM-ISTKR+ 
     &             ISTKR_if_LRCB+ SIZECB)
             ENDIF
             PEAK_LRLU_UD = 
     &             max(PEAK_LRLU_UD,
     &             NRLNEC_if_LR_LU + NRLADULR_UD)
             PEAK_OOC_LRLU_UD = 
     &             max(PEAK_OOC_LRLU_UD,
     &             NRLNEC_ACTIVE + NRLADULR_UD)
             PEAK_OOC_LRLU_WC = 
     &             max(PEAK_OOC_LRLU_WC,
     &             NRLNEC_ACTIVE + NRLADULR_WC)
             PEAK_LRLUCB_UD = 
     &             max(PEAK_LRLUCB_UD,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_UD)
             PEAK_LRLUCB_WC = 
     &             max(PEAK_LRLUCB_WC,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_WC)
             PEAK_OOC_LRLUCB_UD = 
     &             max(PEAK_OOC_LRLUCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_UD)
             PEAK_OOC_LRLUCB_WC = 
     &             max(PEAK_OOC_LRLUCB_WC,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_WC)
             PEAK_LRCB_UD =
     &             max(PEAK_LRCB_UD,
     &             NRLNEC_if_LR_CB + NRLNECLR_CB_UD)
             PEAK_OOC_LRCB_UD =
     &             max(PEAK_OOC_LRCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_CB_UD)
             IF (ABOVE_L0) THEN
              PEAK_DYN_LRLU_UD   = max( PEAK_DYN_LRLU_UD,
     &          SIZECB_UNDER_L0 + NRLADULR_UD )
              PEAK_DYN_LRCB_UD   = max( PEAK_DYN_LRCB_UD,
     &          SIZECB_UNDER_L0_IF_LRCB + ISTKR_if_LRCB)
              PEAK_DYN_LRLUCB_UD = max( PEAK_DYN_LRLUCB_UD,
     &          SIZECB_UNDER_L0_IF_LRCB + ISTKR_if_LRCB + NRLADULR_UD)
              PEAK_DYN_LRLU_WC   = max( PEAK_DYN_LRLU_WC,
     &          SIZECB_UNDER_L0 + NRLADULR_WC)
              PEAK_DYN_LRLUCB_WC = max( PEAK_DYN_LRLUCB_WC,
     &          SIZECB_UNDER_L0 + ISTKRLR_CB_WC + NRLADULR_WC)
             ENDIF
         ENDIF
         IF (UPDATE .AND. LEVEL.EQ.2 .AND. .NOT. MASTER) THEN
             NRLNEC =
     &         max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM+NRLADU_CURRENT)
             NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_CURRENT+
     &         NRLADU_ROOT_3+CURRENT_ACTIVE_MEM)
             IF (.NOT.COMPRESS_PANEL.AND..NOT.FORCE_BFR) THEN
              NRLNEC_if_LR_LU = max(
     &             NRLNEC_if_LR_LU,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM+NRLADU_CURRENT)
              NRLNEC_if_LR_CB   = max(
     &             NRLNEC_if_LR_CB  ,NRLADU +
     &             CURRENT_ACTIVE_MEM-ISTKR +
     &             ISTKR_if_LRCB+NRLADU_CURRENT)
              NRLNEC_if_LR_LUCB = max(
     &             NRLNEC_if_LR_LUCB,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM-ISTKR +
     &             ISTKR_if_LRCB+NRLADU_CURRENT)
              NRLNECOOC_if_LR_LUCB = max(NRLNECOOC_if_LR_LUCB,
     &             NRLADU_CURRENT+
     &             NRLADU_ROOT_3+CURRENT_ACTIVE_MEM-ISTKR+
     &             ISTKR_if_LRCB)
              PEAK_LRLU_UD = 
     &             max(PEAK_LRLU_UD,
     &             NRLNEC_if_LR_LU + NRLADULR_UD)
              PEAK_OOC_LRLU_UD = 
     &             max(PEAK_OOC_LRLU_UD,
     &             NRLNEC_ACTIVE + NRLADULR_UD)
              PEAK_OOC_LRLU_WC = 
     &             max(PEAK_OOC_LRLU_WC,
     &             NRLNEC_ACTIVE + NRLADULR_WC)
              PEAK_LRLUCB_UD = 
     &             max(PEAK_LRLUCB_UD,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_UD)
              PEAK_LRLUCB_WC = 
     &             max(PEAK_LRLUCB_WC,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_WC)
              PEAK_OOC_LRLUCB_UD = 
     &             max(PEAK_OOC_LRLUCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_UD)
              PEAK_OOC_LRLUCB_WC = 
     &             max(PEAK_OOC_LRLUCB_WC,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_WC)
              PEAK_LRCB_UD =
     &             max(PEAK_LRCB_UD,
     &             NRLNEC_if_LR_CB + NRLNECLR_CB_UD)
              PEAK_OOC_LRCB_UD =
     &             max(PEAK_OOC_LRCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_CB_UD)
              IF (ABOVE_L0) THEN
               PEAK_DYN_LRLU_UD   = max( PEAK_DYN_LRLU_UD,
     &              SIZECB_UNDER_L0 + NRLADULR_UD )
               PEAK_DYN_LRCB_UD   = max( PEAK_DYN_LRCB_UD,
     &              SIZECB_UNDER_L0_IF_LRCB + ISTKR_if_LRCB)
               PEAK_DYN_LRLUCB_UD = max( PEAK_DYN_LRLUCB_UD,
     &        SIZECB_UNDER_L0_IF_LRCB + ISTKR_if_LRCB + NRLADULR_UD)
               PEAK_DYN_LRLU_WC   = max( PEAK_DYN_LRLU_WC,
     &              SIZECB_UNDER_L0 + NRLADULR_WC)
               PEAK_DYN_LRLUCB_WC = max( PEAK_DYN_LRLUCB_WC,
     &              SIZECB_UNDER_L0 + ISTKRLR_CB_WC + NRLADULR_WC)
              ENDIF
             ENDIF
         ENDIF
        IF (LEVELF.EQ.3) THEN
          IF (LEVEL.EQ.1) THEN
            LEV3MAXREC = int(min(NCB,LOCAL_M),8) *
     &                   int(min(NCB,LOCAL_N),8)
          ELSE
            LEV3MAXREC = min(SIZECB,
     &                 int(min(NBROWMAX,LOCAL_M),8)
     &                *int(min(NCB,LOCAL_N),8))
          ENDIF
          MAXTEMPCB  = max(MAXTEMPCB, LEV3MAXREC)
          MAXITEMPCB = max(MAXITEMPCB,SIZECBI+SIZEHEADER)
          SBUFR_CB   = max(SBUFR_CB, LEV3MAXREC+int(SIZECBI,8))
          NIRNEC = max(NIRNEC,NIRADU+ISTKI+
     &    min(NCB,LOCAL_M)+ min(NCB,LOCAL_N)+SIZEHEADER)
          NIRNEC_OOC = max(NIRNEC_OOC,NIRADU_OOC+ISTKI_OOC+
     &    min(NCB,LOCAL_M)+ min(NCB,LOCAL_N)+SIZEHEADER)
        ENDIF
        IF (CONCERNED) THEN
         IF (LEVELF.EQ.2) THEN
           IF (UPDATE.AND.(LEVEL.NE.2.OR..NOT.MASTER)) THEN
             IF(MASTERF)THEN
                 NBR = min(NBROWMAXF,NBROWMAX)
             ELSE
                 NBR = min(max(NELIMF,NBROWMAXF),NBROWMAX)
             ENDIF
             IF (KEEP(50).EQ.0) THEN
               CBMAXS = int(NBR,8)*int(NCB,8)
             ELSE
               CBMAXS = int(NBR,8)*int(NCB,8) -
     &                  (int(NBR,8)*int(NBR-1,8))/2_8
             ENDIF
           ELSE
              CBMAXS = 0_8
           END IF
           IF (MASTERF) THEN
             IF (LEVEL.EQ.1) THEN
                IF (.NOT.UPDATE) THEN
                  NBR = min(NELIMF, NCB)
                ELSE
                  NBR = 0
                ENDIF
             ELSE
                NBR = min(NELIMF, NBROWMAX)
             ENDIF
             IF (KEEP(50).EQ.0) THEN
                CBMAXR = int(NBR,8)*NCB8
             ELSE
                CBMAXR = int(NBR,8)*int(min(NCB,NELIMF),8)-
     &                   (int(NBR,8)*int(NBR-1,8))/2_8
                CBMAXR = min(CBMAXR, int(NELIMF,8)*int(NELIMF+1,8)/2_8)
                CBMAXR = min(CBMAXR, SIZECB)
                IF ((LEVEL.EQ.1).AND.(.NOT. PACKED_CB)) THEN
                  CBMAXR = min(CBMAXR,(NCB8*(NCB8+1_8))/2_8)
                ENDIF
             ENDIF
           ELSE IF (UPDATEF) THEN
              NBR = min(NBROWMAXF,NBROWMAX)
              CBMAXR = int(NBR,8) * NCB8
              IF (KEEP(50).NE.0) THEN
                CBMAXR = CBMAXR - (int(NBR,8)*(int(NBR-1,8)))/2_8
              ENDIF
           ELSE
              CBMAXR = 0_8
           ENDIF
         ELSEIF (LEVELF.EQ.3) THEN
           CBMAXR = LEV3MAXREC
           IF (UPDATE.AND. .NOT. (MASTER.AND.LEVEL.EQ.2)) THEN
             CBMAXS = LEV3MAXREC
           ELSE
             CBMAXS = 0_8
           ENDIF
         ELSE
           IF (MASTERF) THEN
             CBMAXS = 0_8
             NBR = min(NFRF,NBROWMAX)
             IF ((LEVEL.EQ.1).AND.UPDATE) THEN
                NBR = 0
             ENDIF
             CBMAXR = int(NBR,8)*int(min(NFRF,NCB),8)
             IF (LEVEL.EQ.2)
     &       CBMAXR = min(CBMAXR, SIZECB_SLAVE)
             IF ( KEEP(50).NE.0 )  THEN
              CBMAXR = min(CBMAXR,(int(NFRF,8)*int(NFRF+1,8))/2_8)
             ELSE
              CBMAXR = min(CBMAXR,int(NFRF,8)*int(NFRF,8))
             ENDIF
           ELSE
             CBMAXR = 0_8
             CBMAXS = SIZECB
           ENDIF
         ENDIF
         IF (UPDATE) THEN
           CBMAXS = min(CBMAXS, SIZECB)
           IF ( .not. ( LEVELF .eq. 1 .AND. UPDATEF ) )THEN
              SBUFS_CB = max(SBUFS_CB, CBMAXS+int(SIZECBI,8))
           ENDIF
         ENDIF
         STACKCB = .FALSE.
         IF (UPDATEF) THEN
          STACKCB = .TRUE.
          SIZECBI     = 2 * NCB + SIZEHEADER
          IF (LEVEL.EQ.1) THEN
             IF (KEEP(50).NE.0.AND.LEVELF.NE.3
     &           .AND.PACKED_CB) THEN
                 SIZECB = (NCB8*(NCB8+1_8))/2_8
             ELSE
                 SIZECB = NCB8*NCB8
             ENDIF
             IF (MASTER) THEN
               IF (MASTERF) THEN
                  SIZECBI     = 2+ XSIZE_IC
               ENDIF
             ELSE IF (LEVELF.EQ.1) THEN
               SIZECB  = min(CBMAXR,SIZECB)
               IF (COMPRESS_CB) THEN
                SIZECBLR_UD = min(SIZECBLR_UD,SIZECB)
                SIZECBLR_WC = min(SIZECBLR_WC,SIZECB)
                SIZECB_if_LRCB = min(SIZECB_if_LRCB,SIZECB)
               ENDIF
               SIZECBI    = 2 * NCB +  9 
               SBUFR_CB   = max(SBUFR_CB, int(SIZECBI,8)+SIZECB)
               SIZECBI    =  2 * NCB + SIZEHEADER     
             ELSE 
               SIZECBI    = 2 * NCB +  9 
               SBUFR_CB   = max(SBUFR_CB, 
     &                      min(SIZECB,CBMAXR) + int(SIZECBI,8))
               MAXTEMPCB  = max(MAXTEMPCB, min(SIZECB,CBMAXR)) 
               SIZECBI    =  2 * NCB + SIZEHEADER 
               MAXITEMPCB = max(MAXITEMPCB, SIZECBI)
               IF ( .NOT. MASTERF ) THEN
                 SIZECBI     = 0
               ELSE
                  SIZECBI = NCB + SIZEHEADER + SLAVEF - 1
               ENDIF
               SIZECB      = 0_8
               SIZECBLR_UD    = 0_8
               SIZECBLR_WC    = 0_8
               SIZECB_if_LRCB = 0_8
             ENDIF
          ELSE 
             SIZECB = SIZECB_SLAVE
             SIZECBLR_UD    = SIZECBLR_SLAVE_UD
             SIZECBLR_WC    = SIZECBLR_SLAVE_WC
             SIZECB_if_LRCB = SIZECB_SLAVE_if_LRCB
             MAXTEMPCB  = max(MAXTEMPCB, min(CBMAXR,SIZECB) )
             MAXITEMPCB = max(MAXITEMPCB,NBROWMAX+NCB+SIZEHEADER)
             IF (.NOT. 
     &        (UPDATE.AND.(.NOT.MASTER).AND.(NSLAVES_LOC.EQ.1))
     &          ) 
     &       SBUFR_CB = max(SBUFR_CB, 
     &            min(CBMAXR,SIZECB) + int(NBROWMAX + NCB + 6,8))
             IF (MASTER) THEN
              SIZECB  = 0_8
              SIZECB_SLAVE_if_LRCB =  0_8
              SIZECBLR_UD          = 0_8
              SIZECBLR_WC          = 0_8
              IF (MASTERF) THEN
                SIZECBI = 2 + XSIZE_IC
              ELSE
                SIZECBI = 0
              ENDIF
             ELSE IF (UPDATE) THEN
              IF (MASTERF) THEN
                SIZECBI      =  NCB + SIZEHEADER + SLAVEF - 1
              ELSE
                SIZECBI = 0
              ENDIF
               IF (KEEP(50).EQ.0) THEN
                 SIZECBI = SIZECBI + NBROWMAX + NFR + 
     &                     SIZEHEADER
               ELSE
                 SIZECBI = SIZECBI + NBROWMAX + NFR +
     &                     SIZEHEADER+ NSLAVES_LOC
               ENDIF
             ELSE
              SIZECB      = 0_8
              IF ( MASTERF ) THEN
                SIZECBI = NCB + SIZEHEADER + SLAVEF - 1
              ELSE
                SIZECBI = 0
              ENDIF
              SIZECB_SLAVE_if_LRCB =  0_8
              SIZECBLR_UD          = 0_8
              SIZECBLR_WC          = 0_8
             ENDIF
          ENDIF
         ELSE
           IF (LEVELF.NE.3) THEN
               STACKCB     = .TRUE.
               SIZECB      = 0_8
               SIZECB_SLAVE_if_LRCB =  0_8
               SIZECBLR_UD          = 0_8
               SIZECBLR_WC          = 0_8
               SIZECBI     = 0
               IF ( (LEVEL.EQ.1) .AND. (LEVELF.NE.1) ) THEN
                  IF (PACKED_CB) THEN 
                      SIZECB  = (NCB8*(NCB8+1_8))/2_8
                  ELSE
                      SIZECB  = NCB8*NCB8
                  ENDIF
                  SIZECBI     = 2 * NCB + SIZEHEADER
               ELSE IF (LEVEL.EQ.2) THEN
                 IF (MASTER) THEN
                   SIZECBI=0
                 ELSE 
                   SIZECB  = SIZECB_SLAVE
                   SIZECBLR_UD    = SIZECBLR_SLAVE_UD
                   SIZECBLR_WC    = SIZECBLR_SLAVE_WC
                   SIZECB_if_LRCB = SIZECB_SLAVE_if_LRCB
                   SIZECBI = NBROWMAX + NFR + SIZEHEADER
                 ENDIF 
               ENDIF
           ENDIF
         ENDIF
         IF (STACKCB) THEN
            IF (FRERE(STEP(INODE)).EQ.0) THEN
                  write(*,*) ' ERROR 3 in DMUMPS_ANA_DISTM'
                  CALL MUMPS_ABORT()
           ENDIF
           ITOP = ITOP + 1
           IF ( ITOP .GT. NSTEPS ) THEN
             WRITE(*,*) 'ERROR 4 in DMUMPS_ANA_DISTM '
           ENDIF
           LSTKI(ITOP) = SIZECBI
           ISTKI=ISTKI + SIZECBI
           ISTKI_OOC = ISTKI_OOC + SIZECBI + (XSIZE_OOC-XSIZE_IC)
           LSTKR(ITOP) = SIZECB
           ISTKR = ISTKR + LSTKR(ITOP)
           NRLNEC = max(NRLNEC,NRLADU+ISTKR+MAXTEMPCB)
           NIRNEC = max(NIRNEC,NIRADU+ISTKI+MAXITEMPCB)
           NIRNEC_OOC = max(NIRNEC_OOC,NIRADU_OOC+ISTKI_OOC+
     &          MAXITEMPCB + 
     &          (XSIZE_OOC-XSIZE_IC) ) 
           NRLNEC_if_LR_LU    =  max(NRLNEC_if_LR_LU,
     &                             NRLADU_if_LR_LU+ISTKR+MAXTEMPCB)
           LSTKR_if_LRCB(ITOP) = SIZECB_if_LRCB
           ISTKR_if_LRCB       = ISTKR_if_LRCB +  LSTKR_if_LRCB(ITOP)
           NRLNEC_if_LR_CB     =  max(NRLNEC_if_LR_CB,
     &                             NRLADU_if_LR_LU+ISTKR_if_LRCB + 
     &                             MAXTEMPCB)
           NRLNEC_if_LR_LUCB  =  max(NRLNEC_if_LR_LUCB,
     &                             NRLADU_if_LR_LU+ISTKR_if_LRCB + 
     &                             MAXTEMPCB)
           LSTKRLR_CB_UD(ITOP) = SIZECBLR_UD
           ISTKRLR_CB_UD       = ISTKRLR_CB_UD + LSTKRLR_CB_UD(ITOP)
           LSTKRLR_CB_WC(ITOP) = SIZECBLR_WC
           ISTKRLR_CB_WC       = ISTKRLR_CB_WC + LSTKRLR_CB_WC(ITOP)
           NRLNECLR_CB_UD    =  max(NRLNECLR_CB_UD, ISTKRLR_CB_UD)
           NRLNECLR_LUCB_UD  =  max(NRLNECLR_LUCB_UD, 
     &                              NRLADULR_UD+ISTKRLR_CB_UD)
           NRLNECLR_LUCB_WC  =  max(NRLNECLR_LUCB_WC, 
     &                              NRLADULR_WC+ISTKRLR_CB_WC)
         ENDIF 
        ENDIF 
         TNSTK(STEP(IFATH)) = TNSTK(STEP(IFATH)) - 1
         IF ( TNSTK(STEP(IFATH)) .EQ. 0 ) THEN
            INODE = IFATH 
            GOTO 95
         ELSE
            GOTO 90
         ENDIF
      ENDIF 
 115  CONTINUE
      IF ( KEEP(53) .NE. 0 ) THEN
        IF ( KEEP(38) .ne. 0 ) THEN
          IROOT = KEEP( 38 )
        ELSE
          IROOT = KEEP( 20 )
        END IF
        ROOT_OWNER = ( MYID .eq.
     &  MUMPS_PROCNODE( PROCNODE(STEP(IROOT)), KEEP(199) ) )
        SIZE_ROOT = ND(STEP(IROOT))+KEEP(253)
        CALL DMUMPS_SVD_QR_ESTIM_WK( PHASE, 
     &       KEEP(51), KEEP(51), SIZE_ROOT,
     &       LOCAL_M, LOCAL_N, ROOT_OWNER, KEEP,KEEP8,
     &       LIWK_RR, LWK_RR )
        IF ( NRLNEC-NRLADU .LT. LWK_RR ) THEN
          NRLNEC = NRLADU + LWK_RR
          NRLNEC_if_LR_LU   = NRLNEC_if_LR_LU + LWK_RR
          NRLNEC_if_LR_CB   = NRLNEC_if_LR_CB + LWK_RR
          NRLNEC_if_LR_LUCB = NRLNEC_if_LR_LUCB + LWK_RR
          NRLNEC_ACTIVE     = NRLNEC_ACTIVE     + LWK_RR
          NRLNECOOC_if_LR_LUCB = NRLNECOOC_if_LR_LUCB  + LWK_RR
          PEAK_LRLU_UD = 
     &             max(PEAK_LRLU_UD,
     &             NRLNEC_if_LR_LU + NRLADULR_UD)
          PEAK_OOC_LRLU_UD = 
     &             max(PEAK_OOC_LRLU_UD,
     &             NRLNEC_ACTIVE + NRLADULR_UD)
          PEAK_OOC_LRLU_WC = 
     &             max(PEAK_OOC_LRLU_WC,
     &             NRLNEC_ACTIVE + NRLADULR_WC)
          PEAK_LRLUCB_UD = 
     &             max(PEAK_LRLUCB_UD,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_UD)
          PEAK_LRLUCB_WC = 
     &             max(PEAK_LRLUCB_WC,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_WC)
          PEAK_OOC_LRLUCB_UD = 
     &             max(PEAK_OOC_LRLUCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_UD)
          PEAK_OOC_LRLUCB_WC = 
     &             max(PEAK_OOC_LRLUCB_WC,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_WC)
          PEAK_LRCB_UD =
     &             max(PEAK_LRCB_UD,
     &             NRLNEC_if_LR_CB + NRLNECLR_CB_UD)
          PEAK_OOC_LRCB_UD =
     &             max(PEAK_OOC_LRCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_CB_UD)
          IF (ABOVE_L0) THEN
           PEAK_DYN_LRLU_UD   = max( PEAK_DYN_LRLU_UD,
     &          SIZECB_UNDER_L0 + NRLADULR_UD )
           PEAK_DYN_LRCB_UD   = max( PEAK_DYN_LRCB_UD,
     &          SIZECB_UNDER_L0_IF_LRCB + ISTKR_if_LRCB)
           PEAK_DYN_LRLUCB_UD = max( PEAK_DYN_LRLUCB_UD,
     &      SIZECB_UNDER_L0_IF_LRCB + ISTKR_if_LRCB + NRLADULR_UD)
           PEAK_DYN_LRLU_WC   = max( PEAK_DYN_LRLU_WC,
     &          SIZECB_UNDER_L0 + NRLADULR_WC)
           PEAK_DYN_LRLUCB_WC = max( PEAK_DYN_LRLUCB_WC,
     &          SIZECB_UNDER_L0 + ISTKRLR_CB_WC + NRLADULR_WC)
          ENDIF
        END IF
        IF ( NIRNEC-NIRADU .LT. LIWK_RR ) THEN
          NIRNEC = NIRADU + LIWK_RR
        END IF
        IF ( NIRNEC_OOC-NIRADU_OOC .LT. LIWK_RR ) THEN
          NIRNEC_OOC = NIRADU_OOC + LIWK_RR
        END IF
      END IF
      NRLNEC = max(NRLNEC, NRLADU+int(KEEP(30),8))
      NRLNEC_ACTIVE = max(NRLNEC_ACTIVE, MAX_SIZE_FACTOR+
     &                    int(KEEP(30),8))
      NRLNEC_if_LR_LU = max(NRLNEC_if_LR_LU,
     &         NRLADU_if_LR_LU + int(KEEP(30),8))
      NRLNEC_if_LR_LUCB = max(NRLNEC_if_LR_LUCB,
     &         NRLADU_if_LR_LU + int(KEEP(30),8))
      NRLNEC_if_LR_CB = max(NRLNEC_if_LR_CB,
     &         NRLADU + int(KEEP(30),8))
      NRLNECOOC_if_LR_LUCB = max(NRLNECOOC_if_LR_LUCB,
     &                    MAX_SIZE_FACTOR+ int(KEEP(30),8))
      PEAK_FR      = SAVE_SIZECB_UNDER_L0 + NRLNEC
      PEAK_FR_OOC  = SAVE_SIZECB_UNDER_L0 + NRLNEC_ACTIVE
      PEAK_LRCB_UD =
     &          max(PEAK_LRCB_UD,
     &          NRLNEC_if_LR_CB + NRLNECLR_CB_UD)
      PEAK_LRLU_UD = 
     &          max(PEAK_LRLU_UD,
     &          NRLNEC_if_LR_LU + NRLADULR_UD)
      PEAK_OOC_LRLU_UD = 
     &          max(PEAK_OOC_LRLU_UD,
     &          NRLNEC_ACTIVE + NRLADULR_UD)
      PEAK_OOC_LRLU_WC = 
     &          max(PEAK_OOC_LRLU_WC,
     &          NRLNEC_ACTIVE + NRLADULR_WC)
      PEAK_LRLUCB_UD = 
     &          max(PEAK_LRLUCB_UD,
     &          NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_UD)
      PEAK_LRLUCB_WC = 
     &          max(PEAK_LRLUCB_WC,
     &          NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_WC)
      PEAK_OOC_LRLUCB_UD = 
     &          max(PEAK_OOC_LRLUCB_UD,
     &          NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_UD)
      PEAK_OOC_LRLUCB_WC = 
     &             max(PEAK_OOC_LRLUCB_WC,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_WC)
      IF (KEEP(60).EQ.1) THEN
        IF (I_PROCESS_SCHUR_K60_1) THEN
          NRLADU          = NRLADU - NRLADU_CURRENT_K60_1
          NRLADU_IF_LR_LU = NRLADU_IF_LR_LU - NRLADU_CURRENT_K60_1
        ENDIF
      ENDIF
      IF (ABOVE_L0) THEN
         PEAK_LRCB_UD =
     &         max(PEAK_LRCB_UD,
     &         NRLNEC_if_LR_CB + SAVE_SIZECB_UNDER_L0_IF_LRCB)
         PEAK_LRCB_UD =
     &         max(PEAK_LRCB_UD,
     &         NRLNEC_if_LR_CB + PEAK_DYN_LRCB_UD)
         PEAK_LRLU_UD = max(PEAK_LRLU_UD, 
     &           NRLNEC_if_LR_LU + SAVE_SIZECB_UNDER_L0)
         PEAK_LRLU_UD = max(PEAK_LRLU_UD, 
     &           NRLNEC_if_LR_LU + PEAK_DYN_LRLU_UD)
         PEAK_OOC_LRLU_UD = 
     &         max(PEAK_OOC_LRLU_UD,
     &         NRLNEC_ACTIVE + SAVE_SIZECB_UNDER_L0)
         PEAK_OOC_LRLU_UD = 
     &         max(PEAK_OOC_LRLU_UD,
     &         NRLNEC_ACTIVE + PEAK_DYN_LRLU_UD)
         PEAK_OOC_LRLU_WC = 
     &         max(PEAK_OOC_LRLU_WC,
     &         NRLNEC_ACTIVE + SAVE_SIZECB_UNDER_L0) 
         PEAK_OOC_LRLU_WC = 
     &         max(PEAK_OOC_LRLU_WC,
     &         NRLNEC_ACTIVE + PEAK_DYN_LRLU_WC)
         PEAK_LRLUCB_UD = 
     &         max(PEAK_LRLUCB_UD,
     &         NRLNEC_if_LR_LUCB + SAVE_SIZECB_UNDER_L0_IF_LRCB)
         PEAK_LRLUCB_UD = 
     &         max(PEAK_LRLUCB_UD,
     &         NRLNEC_if_LR_LUCB + PEAK_DYN_LRLUCB_UD)
         PEAK_LRLUCB_WC = 
     &             max(PEAK_LRLUCB_WC,
     &             NRLNEC_if_LR_LUCB + SAVE_SIZECB_UNDER_L0)
         PEAK_LRLUCB_WC = 
     &         max(PEAK_LRLUCB_WC,
     &         NRLNEC_if_LR_LUCB + PEAK_DYN_LRLUCB_WC)
         PEAK_OOC_LRLUCB_UD = 
     &         max(PEAK_OOC_LRLUCB_UD,
     &         NRLNECOOC_if_LR_LUCB + SAVE_SIZECB_UNDER_L0_IF_LRCB)
         PEAK_OOC_LRLUCB_UD = 
     &         max(PEAK_OOC_LRLUCB_UD,
     &         NRLNECOOC_if_LR_LUCB + PEAK_DYN_LRLUCB_UD)
         PEAK_OOC_LRLUCB_WC = 
     &         max(PEAK_OOC_LRLUCB_WC,
     &         NRLNECOOC_if_LR_LUCB + SAVE_SIZECB_UNDER_L0)
         PEAK_OOC_LRLUCB_WC = 
     &         max(PEAK_OOC_LRLUCB_WC,
     &         NRLNECOOC_if_LR_LUCB + PEAK_DYN_LRLUCB_WC)
      ENDIF
      SBUF_RECOLD = max(SBUFR_CB,
     &        MAXTEMPCB+int(MAXITEMPCB,8)) + 10_8
      SBUF_REC_FR = max(SBUFR_FR, int(min(100000_8,SBUFR_CB))) + 17
      SBUF_REC_LR = max(SBUFR_LR, int(min(100000_8,SBUFR_CB))) + 17
      SBUF_REC_FR = SBUF_REC_FR + 2 * KEEP(127) + SLAVEF - 1 + 7
      SBUF_REC_LR = SBUF_REC_LR + 2 * KEEP(127) + SLAVEF - 1 + 7
      SBUF_SEND_FR = max(SBUFS_FR, int(min(100000_8,SBUFR_CB)))+17
      SBUF_SEND_LR = max(SBUFS_LR, int(min(100000_8,SBUFR_CB)))+17
      IF(KEEP(219).NE.0.AND.KEEP(50) .EQ. 2) THEN
         SBUF_RECOLD = SBUF_RECOLD+int(KEEP(108)+1,8)
         SBUF_REC_FR = SBUF_REC_FR+KEEP(108)+1
         SBUF_REC_LR = SBUF_REC_LR+KEEP(108)+1
         SBUF_SEND_FR = SBUF_SEND_FR+KEEP(108)+1
         SBUF_SEND_LR = SBUF_SEND_LR+KEEP(108)+1
      ENDIF
      IF (SLAVEF.EQ.1) THEN 
         SBUF_RECOLD = 1_8
         SBUF_REC_FR = 1
         SBUF_REC_LR = 1
         SBUF_SEND_FR= 1
         SBUF_SEND_LR= 1
      ENDIF
      DEALLOCATE( LSTKR, TNSTK, LSTKI,
     &     LSTKR_if_LRCB, LSTKRLR_CB_UD, LSTKRLR_CB_WC
     &     )
      IF (ABOVE_L0) THEN 
        KEEP(470) = KEEP(470)+  NBNODES_BLR
      ELSE
        KEEP(470) = NBNODES_BLR
      ENDIF
      IF (.NOT.ABOVE_L0) THEN
       PEAK_FR     = NRLNEC
       PEAK_FR_OOC = NRLNEC_ACTIVE
      ENDIF
      MAXFR = max(MAXFR, MAXFR_UNDER_L0)
      MAX_FRONT_SURFACE_LOCAL = max (MAX_FRONT_SURFACE_LOCAL,
     &                               MAX_FRONT_SURFACE_LOCAL_L0)
      MAX_SIZE_FACTOR         = max (MAX_SIZE_FACTOR, 
     &                               MAX_SIZE_FACTOR_L0)
      ENTRIES_IN_FACTORS_LOC_MASTERS = ENTRIES_IN_FACTORS_LOC_MASTERS +
     &                          ENTRIES_IN_FACTORS_MASTERS_LO
      ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC +
     &                         ENTRIES_IN_FACTORS_UNDER_L0
      OPS_SBTR_LOC = OPS_SBTR_LOC + COST_SUBTREES_UNDER_LO
      OPSA_LOC     = OPSA_LOC + OPSA_UNDER_L0
      OPS_SUBTREE = dble(OPS_SBTR_LOC)
      OPSA        = dble(OPSA_LOC)
      RETURN
      END SUBROUTINE DMUMPS_ANA_DISTM
      SUBROUTINE DMUMPS_ANA_DISTM_UNDERL0OMP( 
     &    LPOOL_B_L0_OMP, IPOOL_B_L0_OMP,
     &    L_VIRT_L0_OMP      , VIRT_L0_OMP  , VIRT_L0_OMP_MAPPING,
     &    L_PHYS_L0_OMP      , PHYS_L0_OMP  , PERM_L0_OMP,
     &    PTR_LEAFS_L0_OMP,
     &    KEEP, N, NE, STEP, FRERE, FILS, DAD, ND, 
     &    MYID, PROCNODE, 
     &    I4_L0, NBSTATS_I4, I8_L0, NBSTATS_I8, NBTHREADS,
     &    SIZECB_UNDER_L0, SIZECB_UNDER_L0_IF_LRCB_UD, 
     &    TNSTK, MAXFR,  
     &    MAX_FRONT_SURFACE_LOCAL, MAX_SIZE_FACTOR,
     &    ENTRIES_IN_FACTORS_LOC, 
     &    ENTRIES_IN_FACTORS_LOC_MASTERS, 
     &    OPS_SBTR_LOC, OPSA_LOC, NBNODES_BLR,
     &    IFLAG, IERROR )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LPOOL_B_L0_OMP,
     &                       L_VIRT_L0_OMP, L_PHYS_L0_OMP
      INTEGER, INTENT(IN) :: IPOOL_B_L0_OMP ( LPOOL_B_L0_OMP )
      INTEGER, INTENT(IN) :: VIRT_L0_OMP         ( L_VIRT_L0_OMP )
      INTEGER, INTENT(IN) :: VIRT_L0_OMP_MAPPING ( L_VIRT_L0_OMP )
      INTEGER, INTENT(IN) :: PHYS_L0_OMP         ( L_PHYS_L0_OMP )
      INTEGER, INTENT(IN) :: PERM_L0_OMP         ( L_PHYS_L0_OMP )
      INTEGER, INTENT(IN) :: PTR_LEAFS_L0_OMP    ( L_PHYS_L0_OMP + 1 )
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(INOUT) :: KEEP(500)
      INTEGER, INTENT(IN) :: NE(KEEP(28))
      INTEGER, INTENT(IN) :: STEP(N)
      INTEGER, INTENT(IN) :: FRERE(KEEP(28))
      INTEGER, INTENT(IN) :: FILS(N)
      INTEGER, INTENT(IN) :: DAD(KEEP(28)), ND(KEEP(28))
      INTEGER, INTENT(IN) :: MYID, PROCNODE(KEEP(28))
      INTEGER, INTENT(IN) :: NBSTATS_I4, NBSTATS_I8, NBTHREADS
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
      INTEGER, INTENT(OUT)   :: TNSTK(KEEP(28))
      INTEGER, INTENT(OUT)   :: I4_L0 (NBSTATS_I4, NBTHREADS)
      INTEGER(8), INTENT(OUT):: I8_L0 (NBSTATS_I8, NBTHREADS)
      INTEGER(8), INTENT(OUT):: ENTRIES_IN_FACTORS_LOC, 
     &             ENTRIES_IN_FACTORS_LOC_MASTERS,
     &             SIZECB_UNDER_L0, SIZECB_UNDER_L0_IF_LRCB_UD
      INTEGER, INTENT(OUT)   :: MAXFR, NBNODES_BLR
      INTEGER(8), INTENT(OUT):: MAX_FRONT_SURFACE_LOCAL, 
     &                          MAX_SIZE_FACTOR
      DOUBLE PRECISION, INTENT(OUT)      :: OPS_SBTR_LOC, OPSA_LOC
      INTEGER, ALLOCATABLE, DIMENSION(:) :: LSTKI
      INTEGER(8), ALLOCATABLE, DIMENSION(:) :: LSTKR
      INTEGER(8), ALLOCATABLE, DIMENSION(:) :: LSTKR_if_LRCB,
     &                                         LSTKRLR_CB_UD,
     &                                         LSTKRLR_CB_WC
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: OPSA_LOC_L0_OMP
      INTEGER :: ITH
      INTEGER :: NSTEPS
      INTEGER :: allocok
      INTEGER(8):: ISTKR, ISTKR_if_LRCB, ISTKRLR_CB_UD, 
     &             ISTKRLR_CB_WC
      INTEGER    :: ISTKI, ISTKI_OOC, ITOP
      NSTEPS    = KEEP(28)
      ALLOCATE( LSTKR(NSTEPS), LSTKI(NSTEPS) , 
     &          LSTKR_if_LRCB(NSTEPS), LSTKRLR_CB_UD(NSTEPS), 
     &          LSTKRLR_CB_WC(NSTEPS),
     &          OPSA_LOC_L0_OMP(NBTHREADS),
     &          
     &           stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IFLAG  =-7
        IERROR =  4*NSTEPS+NBTHREADS
        RETURN
      ENDIF
      TNSTK = NE
      OPSA_LOC_L0_OMP(1:NBTHREADS) = 0.0D0
      OPS_SBTR_LOC                 = 0.0D0
      OPSA_LOC                     = 0.0D0
      I4_L0(1:NBSTATS_I4, 1:NBTHREADS) = 0
      I8_L0(1:NBSTATS_I8, 1:NBTHREADS) = 0_8
      NBNODES_BLR                = 0
      SIZECB_UNDER_L0            = 0_8
      SIZECB_UNDER_L0_IF_LRCB_UD = 0_8
      MAXFR                      = 0
      MAX_FRONT_SURFACE_LOCAL    = 0_8  
      MAX_SIZE_FACTOR            = 0_8   
      ENTRIES_IN_FACTORS_LOC         = 0_8
      ENTRIES_IN_FACTORS_LOC_MASTERS = 0_8
      DO ITH = 1, NBTHREADS
       ISTKI      = 0
       ISTKI_OOC  = 0
       ITOP       = 0
       ISTKR           = 0_8
       ISTKR_if_LRCB   = 0_8
       ISTKRLR_CB_UD   = 0_8
       ISTKRLR_CB_WC   = 0_8
       ISTKR           = 0_8
       CALL DMUMPS_ANA_DISTM_UNDERL0_1THR ( ITH,
     &    LPOOL_B_L0_OMP, IPOOL_B_L0_OMP,
     &    L_VIRT_L0_OMP      , VIRT_L0_OMP  , VIRT_L0_OMP_MAPPING,
     &    L_PHYS_L0_OMP      , PHYS_L0_OMP  , PERM_L0_OMP,
     &    PTR_LEAFS_L0_OMP, KEEP, N, NE, NSTEPS, 
     &    STEP, FRERE, FILS, DAD, ND, MYID, PROCNODE, 
     &    ISTKR, ISTKI, ISTKI_OOC, ISTKR_if_LRCB, ISTKRLR_CB_UD, 
     &    ISTKRLR_CB_WC, ITOP,
     &    LSTKI, LSTKR, LSTKR_if_LRCB, LSTKRLR_CB_UD, LSTKRLR_CB_WC,
     &    I4_L0(1,ITH), I4_L0(2,ITH), I4_L0(3,ITH), I4_L0(4,ITH),
     &    I8_L0(1,ITH), I8_L0(2,ITH), I8_L0(3,ITH), I8_L0(4,ITH), 
     &    I8_L0(5,ITH), I8_L0(6,ITH), I8_L0(7,ITH), I8_L0(8,ITH), 
     &    I8_L0(9,ITH), I8_L0(10,ITH), I8_L0(11,ITH), I8_L0(12,ITH), 
     &    I8_L0(13,ITH), I8_L0(14,ITH), I8_L0(15,ITH), I8_L0(16,ITH), 
     &    I8_L0(17,ITH), I8_L0(18,ITH), I8_L0(19,ITH), I8_L0(20,ITH), 
     &    I8_L0(21,ITH), I8_L0(22,ITH), 
     &    NBNODES_BLR, TNSTK, MAXFR,
     &    MAX_FRONT_SURFACE_LOCAL, MAX_SIZE_FACTOR,
     &    ENTRIES_IN_FACTORS_LOC, ENTRIES_IN_FACTORS_LOC_MASTERS, 
     &    OPS_SBTR_LOC, OPSA_LOC_L0_OMP(ITH), IFLAG, IERROR )
       OPSA_LOC = OPSA_LOC + OPSA_LOC_L0_OMP(ITH)
       I8_L0(23,ITH)     = ISTKR
       SIZECB_UNDER_L0   = SIZECB_UNDER_L0 + ISTKR
       I8_L0(24,ITH)     = ISTKR_if_LRCB + ISTKRLR_CB_UD
       SIZECB_UNDER_L0_IF_LRCB_UD = SIZECB_UNDER_L0_IF_LRCB_UD + 
     &                             ISTKR_if_LRCB + ISTKRLR_CB_UD
      ENDDO
      DEALLOCATE( LSTKR, LSTKI , 
     &          LSTKR_if_LRCB, LSTKRLR_CB_UD, 
     &          LSTKRLR_CB_WC,
     &          OPSA_LOC_L0_OMP)
      RETURN
      END SUBROUTINE DMUMPS_ANA_DISTM_UNDERL0OMP
      SUBROUTINE DMUMPS_ANA_DISTM_UNDERL0_1THR ( ITHREAD,
     &    LPOOL_B_L0_OMP, IPOOL_B_L0_OMP,
     &    L_VIRT_L0_OMP      , VIRT_L0_OMP  , VIRT_L0_OMP_MAPPING,
     &    L_PHYS_L0_OMP      , PHYS_L0_OMP  , PERM_L0_OMP,
     &    PTR_LEAFS_L0_OMP, KEEP, N, NE, NSTEPS, STEP, FRERE, FILS, DAD,
     &    ND, MYID, PROCNODE, ISTKR, ISTKI, ISTKI_OOC, ISTKR_if_LRCB,
     &    ISTKRLR_CB_UD, ISTKRLR_CB_WC, ITOP,
     &    LSTKI, LSTKR, LSTKR_if_LRCB, LSTKRLR_CB_UD, LSTKRLR_CB_WC,
     &    NIRADU, NIRNEC, NIRADU_OOC, NIRNEC_OOC, NRLADU, NRLNEC,
     &    NRLNEC_ACTIVE, NRLADU_if_LR_LU, NRLNEC_if_LR_LU,
     &    NRLNEC_if_LR_LUCB, NRLNECOOC_if_LR_LUCB, NRLNEC_if_LR_CB, 
     &    NRLADULR_UD, NRLADULR_WC, NRLNECLR_CB_UD, 
     &    NRLNECLR_LUCB_UD, NRLNECLR_LUCB_WC,
     &    PEAK_LRLU_UD, PEAK_OOC_LRLU_UD, PEAK_OOC_LRLU_WC,
     &    PEAK_LRLUCB_UD, PEAK_LRLUCB_WC, PEAK_OOC_LRLUCB_UD, 
     &    PEAK_OOC_LRLUCB_WC, PEAK_LRCB_UD, PEAK_OOC_LRCB_UD,
     &    NBNODES_BLR, TNSTK, MAXFR,  
     &    MAX_FRONT_SURFACE_LOCAL, MAX_SIZE_FACTOR,
     &    ENTRIES_IN_FACTORS_LOC, ENTRIES_IN_FACTORS_LOC_MASTERS, 
     &    OPS_SBTR_LOC, OPSA_LOC, IFLAG, IERROR )
      USE DMUMPS_LR_CORE, ONLY : IS_FRONT_BLR_CANDIDATE
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ITHREAD, LPOOL_B_L0_OMP,
     &                       L_VIRT_L0_OMP, L_PHYS_L0_OMP
      INTEGER, INTENT(IN) :: IPOOL_B_L0_OMP ( LPOOL_B_L0_OMP )
      INTEGER, INTENT(IN) :: VIRT_L0_OMP         ( L_VIRT_L0_OMP )
      INTEGER, INTENT(IN) :: VIRT_L0_OMP_MAPPING ( L_VIRT_L0_OMP )
      INTEGER, INTENT(IN) :: PHYS_L0_OMP         ( L_PHYS_L0_OMP )
      INTEGER, INTENT(IN) :: PERM_L0_OMP         ( L_PHYS_L0_OMP )
      INTEGER, INTENT(IN) :: PTR_LEAFS_L0_OMP    ( L_PHYS_L0_OMP + 1 )
      INTEGER, INTENT(IN) :: KEEP(500), N, NSTEPS
      INTEGER, INTENT(IN) :: NE(NSTEPS)
      INTEGER, INTENT(IN) :: STEP(N)
      INTEGER, INTENT(IN) :: FRERE(NSTEPS)
      INTEGER, INTENT(IN) :: FILS(N)
      INTEGER, INTENT(IN) :: DAD(NSTEPS), ND(NSTEPS)
      INTEGER, INTENT(IN)             :: MYID, PROCNODE(NSTEPS)
      DOUBLE PRECISION, INTENT(INOUT) :: OPS_SBTR_LOC
      DOUBLE PRECISION, INTENT(OUT)   :: OPSA_LOC
      INTEGER, INTENT(INOUT)          :: NBNODES_BLR
      INTEGER, INTENT(INOUT)   :: TNSTK(NSTEPS)
      INTEGER, INTENT(INOUT)   :: IFLAG, IERROR
      INTEGER, INTENT(INOUT)   :: MAXFR
      INTEGER(8), INTENT(INOUT):: MAX_FRONT_SURFACE_LOCAL, 
     &                            MAX_SIZE_FACTOR
      INTEGER(8), INTENT(INOUT):: ENTRIES_IN_FACTORS_LOC, 
     &                            ENTRIES_IN_FACTORS_LOC_MASTERS
      INTEGER(8), INTENT(INOUT) :: 
     &              ISTKR, ISTKR_if_LRCB,
     &              ISTKRLR_CB_UD, ISTKRLR_CB_WC
      INTEGER, INTENT(INOUT)    :: ISTKI,  ISTKI_OOC, ITOP
      INTEGER, INTENT(INOUT)    :: LSTKI(NSTEPS)
      INTEGER(8), INTENT(INOUT) :: LSTKR(NSTEPS), 
     &                             LSTKR_if_LRCB(NSTEPS),
     &                             LSTKRLR_CB_UD(NSTEPS),
     &                             LSTKRLR_CB_WC(NSTEPS)
      INTEGER, INTENT(OUT)   :: NIRADU, NIRNEC, NIRADU_OOC, NIRNEC_OOC
      INTEGER(8), INTENT(OUT):: NRLADU, NRLNEC, NRLNEC_ACTIVE, 
     &       NRLADU_if_LR_LU, NRLNEC_if_LR_LU, NRLNEC_if_LR_LUCB, 
     &       NRLNECOOC_if_LR_LUCB, NRLNEC_if_LR_CB, 
     &       NRLADULR_UD, NRLADULR_WC, NRLNECLR_CB_UD, 
     &       NRLNECLR_LUCB_UD, NRLNECLR_LUCB_WC,
     &       PEAK_LRLU_UD, PEAK_OOC_LRLU_UD, PEAK_OOC_LRLU_WC,
     &       PEAK_LRLUCB_UD, PEAK_LRLUCB_WC, PEAK_OOC_LRLUCB_UD, 
     &       PEAK_OOC_LRLUCB_WC, PEAK_LRCB_UD, PEAK_OOC_LRCB_UD
      LOGICAL :: INSSARBR
      INTEGER :: VIRTUAL_TASK, PHYSICAL_TASK
      INTEGER :: INODE, IFATH, I
      INTEGER    :: SIZEHEADER, SIZEHEADER_OOC, XSIZE_OOC
      INTEGER    :: EXTRA_PERM_INFO_OOC
      LOGICAL :: PACKED_CB
      INTEGER(8) :: NRLADU_ROOT_3
      INTEGER :: FLAG_L0OMP
      PARAMETER (FLAG_L0OMP=-2014)
      INCLUDE 'mumps_headers.h'
      IF (KEEP(50).eq.2) THEN
        EXTRA_PERM_INFO_OOC = 1
      ELSE IF (KEEP(50).eq.0) THEN
        EXTRA_PERM_INFO_OOC = 2
      ELSE
        EXTRA_PERM_INFO_OOC = 0
      ENDIF
      PACKED_CB=( KEEP(215).EQ.0 .AND. KEEP(50).NE.0 )
      NRLADU_ROOT_3 = 0_8
#if defined(OLD_OOC_NOPANEL)
      XSIZE_OOC=XSIZE_OOC_NOPANEL
#else
      IF (KEEP(50).EQ.0) THEN
              XSIZE_OOC=XSIZE_OOC_UNSYM
      ELSE
              XSIZE_OOC=XSIZE_OOC_SYM
      ENDIF
#endif
      SIZEHEADER_OOC = XSIZE_OOC+6 
      SIZEHEADER = XSIZE_IC + 6  
      DO VIRTUAL_TASK = 1, L_VIRT_L0_OMP - 1
        IF (VIRT_L0_OMP_MAPPING(VIRTUAL_TASK) .EQ. ITHREAD) THEN
          DO PHYSICAL_TASK=
     &       VIRT_L0_OMP ( VIRTUAL_TASK ),
     &       VIRT_L0_OMP ( VIRTUAL_TASK + 1 ) - 1
            DO I = PTR_LEAFS_L0_OMP ( PERM_L0_OMP ( PHYSICAL_TASK ) ),
     &        PTR_LEAFS_L0_OMP ( PERM_L0_OMP ( PHYSICAL_TASK )+1 ) +1,
     &              - 1
              INODE = IPOOL_B_L0_OMP(I)
              IF (INODE .LE. 0) THEN
                CYCLE
              ENDIF
 10           CONTINUE
              IFATH = DAD(STEP(INODE))
              CALL DMUMPS_PROCESS_NODE_UNDERL0 ()
              IF (IFATH .NE. 0) THEN
                  TNSTK( STEP(IFATH) ) = TNSTK( STEP(IFATH) ) - 1
              ENDIF
              IF ( PHYS_L0_OMP ( PERM_L0_OMP ( PHYSICAL_TASK ) )
     &             .EQ. INODE ) THEN
                TNSTK(STEP(INODE)) = FLAG_L0OMP
              ELSE IF ( TNSTK( STEP(IFATH) ) .EQ. 0 ) THEN
                    INODE = IFATH
                    GOTO 10
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      RETURN
      CONTAINS
      SUBROUTINE DMUMPS_PROCESS_NODE_UNDERL0
      USE MUMPS_TOOLS_COMMON_M, ONLY : MUMPS_LR_COMPRESS_CB,
     &                                 MUMPS_LR_COMPRESS_PANEL,
     &                                 MUMPS_LR_FORCE_BFR
      IMPLICIT NONE
      INTEGER    :: LRSTATUS, IDUMMY
      LOGICAL    :: COMPRESS_PANEL, COMPRESS_CB
      LOGICAL    :: FORCE_BFR
      INTEGER    :: STKI
      INTEGER(8) :: LSTK
      INTEGER   :: K, NFR, NFRF, NELIM, NELIMF, NCB, NSTK,
     &             LEVEL, LEVELF, IN,
     &             MAXITEMPCB, PANEL_SIZE, SIZECBI
      INTEGER(8):: NFR8, NCB8, 
     &             K464_8, K465_8,
     &             CURRENT_ACTIVE_MEM,
     &             ENTRIES_NODE_LOWER_PART, ENTRIES_NODE_UPPER_PART, 
     &             NRLADU_CURRENT, NRLADU_CURRENT_MISSING
      INTEGER(8) :: SIZEFRNOCBLU
      INTEGER    :: IDUMMY_ARRAY(1)
      INTEGER(8):: SIZECB, SIZECBINFR 
      INTEGER(8):: SIZECB_if_LRCB 
      INTEGER(8):: SIZECBLR_UD, SIZECBLR_WC
      LOGICAL          :: MASTER, MASTERF, STACKCB
      DOUBLE PRECISION :: OPS_NODE
      INTRINSIC int
      INTEGER DMUMPS_OOC_GET_PANEL_SIZE
      EXTERNAL DMUMPS_OOC_GET_PANEL_SIZE
      INTEGER MUMPS_PROCNODE, MUMPS_TYPENODE
      LOGICAL MUMPS_IN_OR_ROOT_SSARBR
      EXTERNAL MUMPS_MAX_SURFCB_NBROWS
      EXTERNAL MUMPS_PROCNODE, MUMPS_TYPENODE, 
     &         MUMPS_IN_OR_ROOT_SSARBR
      MAXITEMPCB = 0
      STACKCB    = .TRUE.
      NFR    = ND(STEP(INODE))+KEEP(253)
      NFR8   = int(NFR,8)
      NSTK   = NE(STEP(INODE))
      NELIM = 0 
        IN = INODE
 100    NELIM = NELIM + 1 
        IN = FILS(IN)
        IF (IN .GT. 0 ) GOTO 100
      IFATH = DAD(STEP(INODE))
      MASTER = MUMPS_PROCNODE(PROCNODE(STEP(INODE)),KEEP(199))
     &           .EQ. MYID
      LEVEL  = MUMPS_TYPENODE(PROCNODE(STEP(INODE)),KEEP(199))
      INSSARBR = MUMPS_IN_OR_ROOT_SSARBR(PROCNODE(STEP(INODE)),
     &        KEEP(199))
      IDUMMY_ARRAY(1) = 0
      NCB      = NFR-NELIM
      NCB8     = int(NCB,8)
      SIZECBINFR = NCB8*NCB8
      IF (KEEP(50).EQ.0) THEN
        SIZECB = SIZECBINFR
      ELSE
        IF ( PACKED_CB ) THEN
          SIZECB    = (NCB8*(NCB8+1_8))/2_8
        ELSE
          SIZECB    = SIZECBINFR
        ENDIF
      ENDIF
      IDUMMY = -99999
      CALL IS_FRONT_BLR_CANDIDATE (INODE, LEVEL, NFR, NELIM, 
     &                KEEP(494), KEEP(473), 1, KEEP(490),
     &                KEEP(491), KEEP(492),
     &                KEEP(20), KEEP(60), DAD(STEP(INODE)), KEEP(38),
     &                KEEP(123), LRSTATUS, IDUMMY)
      FORCE_BFR      = MUMPS_LR_FORCE_BFR(LRSTATUS)
      COMPRESS_PANEL = MUMPS_LR_COMPRESS_PANEL(LRSTATUS)
      COMPRESS_CB    = MUMPS_LR_COMPRESS_CB(LRSTATUS)
      IF (FORCE_BFR.OR.COMPRESS_PANEL.OR.COMPRESS_CB)
     &  NBNODES_BLR = NBNODES_BLR+1
      IF (COMPRESS_PANEL) THEN
        K464_8 = int(KEEP(464),8)
      ELSE
        K464_8 = 1000_8
      ENDIF
      IF (COMPRESS_CB) THEN
        K465_8 = int(KEEP(465),8)
        SIZECB_if_LRCB = 0_8
        SIZECBLR_UD    = SIZECB*K465_8/1000_8
        SIZECBLR_WC    = SIZECB
      ELSE
        K465_8 = 1000_8
        SIZECBLR_UD    = 0_8
        SIZECBLR_WC    = 0_8
        SIZECB_if_LRCB = SIZECB 
      ENDIF
      SIZECBI      = 2* NCB  + SIZEHEADER
      IF(KEEP(86).EQ.1)THEN
         IF(MASTER.AND.(.NOT.MUMPS_IN_OR_ROOT_SSARBR(
     &        PROCNODE(STEP(INODE)), KEEP(199)))
     &     )THEN
               MAX_FRONT_SURFACE_LOCAL=max(MAX_FRONT_SURFACE_LOCAL,
     &              NFR8*NFR8)
         ENDIF
      ENDIF
      NIRADU     = NIRADU + 2*NFR + SIZEHEADER
      NIRADU_OOC = NIRADU_OOC + 2*NFR + SIZEHEADER_OOC
      PANEL_SIZE = DMUMPS_OOC_GET_PANEL_SIZE(
     &      2_8*int(KEEP(226),8), NFR, KEEP(227), KEEP(50))
      NIRADU_OOC = NIRADU_OOC +
     &      EXTRA_PERM_INFO_OOC*(2+NELIM + NELIM/PANEL_SIZE+1)
      IF (KEEP(50).EQ.0) THEN
              IF (KEEP(251) .EQ. 2) THEN
                NRLADU_CURRENT = int(NELIM,8) * int(NFR,8)
                SIZEFRNOCBLU   = int(NFR-NELIM,8)*int(NELIM)
              ELSE
                NRLADU_CURRENT = int(NELIM,8)*int(2*NFR-NELIM,8)
                SIZEFRNOCBLU   = 0_8
              ENDIF
              NRLADU = NRLADU + NRLADU_CURRENT
              MAX_SIZE_FACTOR  = max(MAX_SIZE_FACTOR,NRLADU_CURRENT)
      ELSE
              IF ( KEEP(459) .GT. 1 .AND. KEEP(50) .NE. 0 ) THEN
                CALL MUMPS_LDLTPANEL_STORAGE( NELIM, KEEP, 
     &          IDUMMY_ARRAY(1), NRLADU_CURRENT )
                SIZEFRNOCBLU = int(NELIM,8)*int(NELIM,8) -
     &          NRLADU_CURRENT
                NRLADU_CURRENT = NRLADU_CURRENT +
     &          int(NELIM,8) * int(NFR-NELIM,8)
              ELSE
                SIZEFRNOCBLU = 0_8
                NRLADU_CURRENT = int(NELIM,8)*int(NFR,8)
              ENDIF
              NRLADU = NRLADU + NRLADU_CURRENT
              MAX_SIZE_FACTOR=max(MAX_SIZE_FACTOR,NRLADU_CURRENT)
      ENDIF
      IF (FORCE_BFR.OR.COMPRESS_PANEL) THEN
               NRLADU_if_LR_LU = NRLADU_if_LR_LU + 0_8
               NRLADU_CURRENT_MISSING = NRLADU_CURRENT
               NRLADULR_UD     = NRLADULR_UD + 
     &                           NRLADU_CURRENT*K464_8/1000_8
               NRLADULR_WC     = NRLADULR_WC + 
     &                           NRLADU_CURRENT
      ELSE
               NRLADU_CURRENT_MISSING = 0_8
               NRLADU_if_LR_LU = NRLADU_if_LR_LU + NRLADU_CURRENT
      ENDIF
      SIZECBI        = 2* NCB  + SIZEHEADER
      NIRNEC = max(NIRNEC,
     &             NIRADU+ISTKI+SIZECBI+MAXITEMPCB)
      NIRNEC_OOC = max(NIRNEC_OOC,
     &             NIRADU_OOC+ISTKI_OOC+SIZECBI+MAXITEMPCB +
     &             (XSIZE_OOC-XSIZE_IC) ) 
      CURRENT_ACTIVE_MEM = ISTKR+SIZECBINFR
      IF ( KEEP(50).NE.0 .AND. LEVEL.EQ.1 ) THEN
          SIZEFRNOCBLU = SIZEFRNOCBLU + int(NELIM,8)*int(NCB,8)
      ENDIF
      CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM + SIZEFRNOCBLU
      NRLNEC_if_LR_CB   = 
     &             max(NRLNEC_if_LR_CB  ,NRLADU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB)
      NRLNEC_if_LR_LUCB = 
     &             max(NRLNEC_if_LR_LUCB,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB+
     &             NRLADU_CURRENT_MISSING)
      NRLNECOOC_if_LR_LUCB = 
     &             max(NRLNECOOC_if_LR_LUCB,
     &             NRLADU_CURRENT+NRLADU_ROOT_3+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB)
      IF (NSTK .NE. 0 .AND. INSSARBR .AND.
     &     KEEP(234).NE.0 .AND. KEEP(55).EQ.0) THEN
          CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM - LSTKR(ITOP)
      ENDIF
      NRLNEC = max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM)
      NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_CURRENT+
     &             NRLADU_ROOT_3+CURRENT_ACTIVE_MEM)
      NRLNEC_if_LR_LU = 
     &             max(NRLNEC_if_LR_LU,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM+NRLADU_CURRENT_MISSING)
      NRLNECLR_LUCB_UD =  max(NRLNECLR_LUCB_UD,
     &                           NRLADULR_UD+ISTKRLR_CB_UD)
      NRLNECLR_LUCB_WC =  max(NRLNECLR_LUCB_WC,
     &                           NRLADULR_WC+ISTKRLR_CB_WC)
      PEAK_LRLU_UD = 
     &             max(PEAK_LRLU_UD,
     &             NRLNEC_if_LR_LU + NRLADULR_UD)
      PEAK_OOC_LRLU_UD = 
     &             max(PEAK_OOC_LRLU_UD,
     &             NRLNEC_ACTIVE + NRLADULR_UD)
      PEAK_OOC_LRLU_WC = 
     &             max(PEAK_OOC_LRLU_WC,
     &             NRLNEC_ACTIVE + NRLADULR_WC)
      PEAK_LRLUCB_UD = 
     &             max(PEAK_LRLUCB_UD,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_UD)
      PEAK_LRLUCB_WC = 
     &             max(PEAK_LRLUCB_WC,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_WC)
      PEAK_OOC_LRLUCB_UD = 
     &             max(PEAK_OOC_LRLUCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_UD)
      PEAK_OOC_LRLUCB_WC = 
     &             max(PEAK_OOC_LRLUCB_WC,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_WC)
      PEAK_LRCB_UD =
     &             max(PEAK_LRCB_UD,
     &             NRLNEC_if_LR_CB + NRLNECLR_CB_UD)
      PEAK_OOC_LRCB_UD =
     &             max(PEAK_OOC_LRCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_CB_UD)
      IF (NFR.GT.MAXFR) MAXFR = NFR
      IF (NSTK.GT.0) THEN
            DO 70 K=1,NSTK
               LSTK = LSTKR(ITOP)
               ISTKR = ISTKR - LSTK
               IF (K==1 .AND. INSSARBR.AND.KEEP(234).NE.0
     &            .AND.KEEP(55).EQ.0) THEN
               ELSE
                 CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM - LSTK
               ENDIF
               LSTK          = LSTKR_if_LRCB(ITOP)
               ISTKR_if_LRCB = ISTKR_if_LRCB - LSTK
               LSTK          = LSTKRLR_CB_UD(ITOP) 
               ISTKRLR_CB_UD = ISTKRLR_CB_UD - LSTK
               LSTK          = LSTKRLR_CB_WC(ITOP) 
               ISTKRLR_CB_WC = ISTKRLR_CB_WC - LSTK
               STKI  = LSTKI( ITOP )
               ISTKI = ISTKI - STKI
               ISTKI_OOC = ISTKI_OOC - STKI - (XSIZE_OOC-XSIZE_IC)
               ITOP = ITOP - 1
               IF (ITOP.LT.0) THEN
                  write(*,*) MYID,
     &            ': ERROR 2 in DMUMPS_ANA_DISTM. ITOP = ',ITOP
                  CALL MUMPS_ABORT()
               ENDIF
 70         CONTINUE
      ENDIF 
      IF (
     &        ( (INODE.NE.KEEP(20)).OR.(KEEP(60).EQ.0) ) 
     &       .AND.
     &        ( (INODE.NE.KEEP(38)).OR.(KEEP(60).LE.1) ) 
     &      )
     &  THEN
            ENTRIES_NODE_LOWER_PART = int(NFR-NELIM,8) * int(NELIM,8)
            IF ( KEEP(50).EQ.0 ) THEN
              ENTRIES_NODE_UPPER_PART = int(NFR,8) * int(NELIM,8)
            ELSE
              ENTRIES_NODE_UPPER_PART =
     &        (int(NELIM,8)*int(NELIM+1,8))/2_8
            ENDIF
            IF (KEEP(50).EQ.2 .AND. LEVEL.EQ.3) THEN
              CALL MUMPS_GET_FLOPS_COST(NFR, 
     &           NELIM, NELIM, 0,
     &           1,OPS_NODE)
            ELSE
              CALL MUMPS_GET_FLOPS_COST(NFR, 
     &           NELIM, NELIM,KEEP(50),
     &           1,OPS_NODE)
            ENDIF
      ELSE
           OPS_NODE = 0.0D0
           ENTRIES_NODE_UPPER_PART = 0_8
           ENTRIES_NODE_LOWER_PART = 0_8
      ENDIF
      ENTRIES_IN_FACTORS_LOC_MASTERS = 
     &                     ENTRIES_IN_FACTORS_LOC_MASTERS +
     &                            ENTRIES_NODE_UPPER_PART +
     &                            ENTRIES_NODE_LOWER_PART
      ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC_MASTERS
      OPSA_LOC = OPSA_LOC + dble(OPS_NODE)
      IF (MUMPS_IN_OR_ROOT_SSARBR(PROCNODE(STEP(INODE)),
     &        KEEP(199)) .OR. NE(STEP(INODE))==0) THEN
             OPS_SBTR_LOC = OPS_SBTR_LOC + OPS_NODE
      ENDIF
      IF (IFATH .EQ. 0) THEN
       RETURN
      ELSE
         NFRF = ND(STEP(IFATH))+KEEP(253)
         IF (DAD(STEP(IFATH)).EQ.0) THEN
           NELIMF = NFRF
         ELSE
           NELIMF = 0
           IN = IFATH
           DO WHILE (IN.GT.0)
              IN = FILS(IN)
              NELIMF = NELIMF+1
           ENDDO
         ENDIF
         LEVELF = MUMPS_TYPENODE(PROCNODE(STEP(IFATH)),KEEP(199))
         MASTERF= MUMPS_PROCNODE(PROCNODE(STEP(IFATH)),
     &                           KEEP(199)).EQ.MYID
         IF ( INSSARBR .AND. KEEP(234).NE.0) THEN
               NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_CURRENT+
     &           NRLADU_ROOT_3+CURRENT_ACTIVE_MEM)
               NRLNEC = max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM)
               NRLNEC_if_LR_LU = 
     &             max(NRLNEC_if_LR_LU,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM+
     &             NRLADU_CURRENT_MISSING)
               NRLNEC_if_LR_CB   = 
     &             max(NRLNEC_if_LR_CB  ,NRLADU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB)
               NRLNEC_if_LR_LUCB = 
     &             max(NRLNEC_if_LR_LUCB,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB+
     &             NRLADU_CURRENT_MISSING)
               NRLNECOOC_if_LR_LUCB = max(NRLNECOOC_if_LR_LUCB,
     &             NRLADU_CURRENT+NRLADU_ROOT_3+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB)
         ELSE
               NRLNEC = max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM+SIZECB)
               NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_CURRENT+
     &           NRLADU_ROOT_3+CURRENT_ACTIVE_MEM+SIZECB)
               NRLNEC_if_LR_LU = 
     &             max(NRLNEC_if_LR_LU,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM+SIZECB+
     &             NRLADU_CURRENT_MISSING)
               NRLNEC_if_LR_CB   = 
     &             max(NRLNEC_if_LR_CB  ,NRLADU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ 
     &             ISTKR_if_LRCB+ SIZECB)
               NRLNEC_if_LR_LUCB = 
     &             max(NRLNEC_if_LR_LUCB,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ 
     &             ISTKR_if_LRCB+ SIZECB+
     &             NRLADU_CURRENT_MISSING)
               NRLNECOOC_if_LR_LUCB = max(NRLNECOOC_if_LR_LUCB,
     &             NRLADU_CURRENT+NRLADU_ROOT_3+
     &             CURRENT_ACTIVE_MEM-ISTKR+ 
     &             ISTKR_if_LRCB+ SIZECB)
         ENDIF
         PEAK_LRLU_UD = 
     &             max(PEAK_LRLU_UD,
     &             NRLNEC_if_LR_LU + NRLADULR_UD)
         PEAK_OOC_LRLU_UD = 
     &             max(PEAK_OOC_LRLU_UD,
     &             NRLNEC_ACTIVE + NRLADULR_UD)
         PEAK_OOC_LRLU_WC = 
     &             max(PEAK_OOC_LRLU_WC,
     &             NRLNEC_ACTIVE + NRLADULR_WC)
         PEAK_LRLUCB_UD = 
     &             max(PEAK_LRLUCB_UD,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_UD)
         PEAK_LRLUCB_WC = 
     &             max(PEAK_LRLUCB_WC,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_WC)
         PEAK_OOC_LRLUCB_UD = 
     &             max(PEAK_OOC_LRLUCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_UD)
         PEAK_OOC_LRLUCB_WC = 
     &             max(PEAK_OOC_LRLUCB_WC,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_WC)
         PEAK_LRCB_UD =
     &             max(PEAK_LRCB_UD,
     &             NRLNEC_if_LR_CB + NRLNECLR_CB_UD)
         PEAK_OOC_LRCB_UD =
     &             max(PEAK_OOC_LRCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_CB_UD)
         IF (KEEP(50).NE.0.AND.LEVELF.NE.3
     &           .AND.PACKED_CB) THEN
                 SIZECB = (NCB8*(NCB8+1_8))/2_8
         ELSE
                 SIZECB = NCB8*NCB8
         ENDIF
         SIZECBI     = 2+ XSIZE_IC
         IF (STACKCB) THEN
            IF (FRERE(STEP(INODE)).EQ.0) THEN
                  write(*,*) ' ERROR 3 in DMUMPS_ANA_DISTM'
                  CALL MUMPS_ABORT()
           ENDIF
           ITOP = ITOP + 1
           IF ( ITOP .GT. NSTEPS ) THEN
             WRITE(*,*) 'ERROR 4 in DMUMPS_ANA_DISTM '
             CALL MUMPS_ABORT()
           ENDIF
           LSTKI(ITOP) = SIZECBI
           ISTKI=ISTKI + SIZECBI
           ISTKI_OOC = ISTKI_OOC + SIZECBI + (XSIZE_OOC-XSIZE_IC)
           LSTKR(ITOP) = SIZECB
           ISTKR = ISTKR + LSTKR(ITOP)
           NRLNEC = max(NRLNEC,NRLADU+ISTKR)
           NIRNEC = max(NIRNEC,NIRADU+ISTKI)
           NIRNEC_OOC = max(NIRNEC_OOC,NIRADU_OOC+ISTKI_OOC)
           NRLNEC_if_LR_LU  =  max(NRLNEC_if_LR_LU,
     &                             NRLADU_if_LR_LU+ISTKR)
           LSTKR_if_LRCB(ITOP) = SIZECB_if_LRCB
           ISTKR_if_LRCB       = ISTKR_if_LRCB +  LSTKR_if_LRCB(ITOP)
           LSTKRLR_CB_UD(ITOP) = SIZECBLR_UD
           ISTKRLR_CB_UD       = ISTKRLR_CB_UD + LSTKRLR_CB_UD(ITOP)
           LSTKRLR_CB_WC(ITOP) = SIZECBLR_WC
           ISTKRLR_CB_WC       = ISTKRLR_CB_WC + LSTKRLR_CB_WC(ITOP)
           NRLNECLR_CB_UD    =  max(NRLNECLR_CB_UD, ISTKRLR_CB_UD)
           NRLNECLR_LUCB_UD  =  max(NRLNECLR_LUCB_UD, 
     &                              NRLADULR_UD+ISTKRLR_CB_UD)
           NRLNECLR_LUCB_WC  =  max(NRLNECLR_LUCB_WC, 
     &                              NRLADULR_WC+ISTKRLR_CB_WC)
         ENDIF 
      ENDIF 
      END SUBROUTINE DMUMPS_PROCESS_NODE_UNDERL0
      END SUBROUTINE DMUMPS_ANA_DISTM_UNDERL0_1THR
      SUBROUTINE DMUMPS_PREP_ANA_DISTM_ABOVEL0 (
     &    N, SLAVEF, COMM, MYID, 
     &    STEP, DAD, ICNTL, LP, LPOK, INFO, 
     &    PHYS_L0_OMP, L_PHYS_L0_OMP, 
     &    IPOOL_A_L0_OMP, LPOOL_A_L0_OMP,
     &    KEEP, TNSTK_afterL0, 
     &    FLAGGED_LEAVES
     &  )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER, INTENT(IN) :: N, SLAVEF, COMM, MYID, ICNTL(60), 
     &                       LP 
      LOGICAL, INTENT(IN) :: LPOK
      INTEGER, INTENT(IN) :: L_PHYS_L0_OMP, LPOOL_A_L0_OMP
      INTEGER, INTENT(IN) :: PHYS_L0_OMP(max(1,L_PHYS_L0_OMP)),
     &                   IPOOL_A_L0_OMP(max(1,LPOOL_A_L0_OMP))
      INTEGER, INTENT(IN) :: KEEP(500)
      INTEGER, INTENT(IN) :: STEP(N), DAD(KEEP(28))
      INTEGER, INTENT(OUT)   :: FLAGGED_LEAVES(KEEP(28))
      INTEGER, INTENT(INOUT) :: TNSTK_afterL0(KEEP(28)), INFO(80)
      INTEGER :: ISLAVE, IERR, INODE, I, NSTEPS, allocok
      INTEGER :: SIZE_BUFREC,  Itemp, SIZE_RECEIVED
      INTEGER, ALLOCATABLE, DIMENSION(:) :: BUFREC
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IREQ
      INTEGER :: MPI_STATUS(MPI_STATUS_SIZE)
      SIZE_BUFREC = 0
      CALL MPI_ALLREDUCE(L_PHYS_L0_OMP, Itemp, 1, 
     &        MPI_INTEGER, MPI_MAX, COMM, IERR)
      SIZE_BUFREC = Itemp
      CALL MPI_ALLREDUCE(LPOOL_A_L0_OMP, Itemp, 1, 
     &        MPI_INTEGER, MPI_MAX, COMM, IERR)
      SIZE_BUFREC= max(SIZE_BUFREC, Itemp)
      ALLOCATE(IREQ(SLAVEF), BUFREC(SIZE_BUFREC), stat=allocok)
      IF (allocok.GT.0) THEN
             IF ( LPOK ) THEN
               WRITE(LP, '(A)') 
     &    ' Allocation failed in DMUMPS_PREP_ANA_DISTM_ABOVEL0'
             END IF
             INFO(1)= -7
             INFO(2)= SLAVEF+SIZE_BUFREC
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) RETURN
      NSTEPS = KEEP(28)
      DO I=1, NSTEPS
        FLAGGED_LEAVES(I) = 0
      ENDDO
      DO ISLAVE=1, SLAVEF
        IF (MYID.EQ.ISLAVE-1) CYCLE
        CALL MPI_ISEND( IPOOL_A_L0_OMP, LPOOL_A_L0_OMP,
     &          MPI_INTEGER, ISLAVE - 1, F_IPOOLAFTER, COMM,
     &          IREQ( ISLAVE ), IERR )
      ENDDO
      DO ISLAVE=1, SLAVEF
        IF (MYID.EQ.ISLAVE-1) CYCLE
        CALL MPI_RECV( BUFREC(1), SIZE_BUFREC,
     &           MPI_INTEGER, ISLAVE-1,
     &           F_IPOOLAFTER, COMM, MPI_STATUS, IERR )
        CALL MPI_GET_COUNT(MPI_STATUS, MPI_INTEGER,
     &                       SIZE_RECEIVED, IERR)
        DO I=1,SIZE_RECEIVED
          INODE = BUFREC(I)
          FLAGGED_LEAVES(STEP(INODE))=INODE
        ENDDO
      ENDDO
      IF (LPOOL_A_L0_OMP.GT.0) THEN
        DO I=1, LPOOL_A_L0_OMP
          INODE = IPOOL_A_L0_OMP(I)
          FLAGGED_LEAVES(STEP(INODE))=INODE
        ENDDO
      ENDIF
      DO ISLAVE=1, SLAVEF
         IF (MYID.EQ.ISLAVE-1) CYCLE
         CALL MPI_WAIT( IREQ( ISLAVE ), MPI_STATUS, IERR )
      ENDDO
      DO ISLAVE=1, SLAVEF
        IF (MYID.EQ.ISLAVE-1) CYCLE
        CALL MPI_ISEND(PHYS_L0_OMP, L_PHYS_L0_OMP,
     &          MPI_INTEGER, ISLAVE - 1, F_PHYS_L0, COMM,
     &          IREQ( ISLAVE ), IERR )
      ENDDO
      DO ISLAVE=1, SLAVEF
        IF (MYID.EQ.ISLAVE-1) CYCLE
        CALL MPI_RECV( BUFREC(1), SIZE_BUFREC,
     &           MPI_INTEGER, ISLAVE-1,
     &           F_PHYS_L0, COMM, MPI_STATUS, IERR )
        CALL MPI_GET_COUNT(MPI_STATUS, MPI_INTEGER,
     &                       SIZE_RECEIVED, IERR)
        DO I=1,SIZE_RECEIVED
          INODE = BUFREC(I)
          IF (DAD(STEP(INODE)).NE.0) THEN
           TNSTK_afterL0(STEP(DAD(STEP(INODE)))) 
     &       =  TNSTK_afterL0(STEP(DAD(STEP(INODE)))) - 1
          ENDIF
        ENDDO
      ENDDO
      DO ISLAVE=1, SLAVEF
         IF (MYID.EQ.ISLAVE-1) CYCLE
         CALL MPI_WAIT( IREQ( ISLAVE ), MPI_STATUS, IERR )
      ENDDO
      IF (allocated(IREQ))   DEALLOCATE(IREQ)
      IF (allocated(BUFREC)) DEALLOCATE(BUFREC)
      RETURN
      END SUBROUTINE DMUMPS_PREP_ANA_DISTM_ABOVEL0 
