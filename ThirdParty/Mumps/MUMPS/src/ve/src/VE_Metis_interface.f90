!
! -----------------------------------------------------------------------
! 
! NEC Vector Host MUMPS Library Interface
! 
! -----------------------------------------------------------------------
!

!
! -----------------------------------------------------------------------
! Interface to C metis_setdefaultoptions
! -----------------------------------------------------------------------
!

INTEGER FUNCTION VE_FMETIS_SETDEFAULTOPTIONS( OPTIONS )
  USE ISO_C_BINDING, ONLY : C_INT
  ! INTERFACE TO C ROUTINE
  INTEGER(KIND=C_INT) :: OPTIONS( * ), IERR(1)
  INTERFACE
    INTEGER(KIND=C_INT) FUNCTION VE_METIS_SETDEFAULTOPTIONS( OPTIONS ) &
                                 BIND(C, NAME='VE_Metis_setdefaultoptions')
      USE ISO_C_BINDING, ONLY : C_INT
      INTEGER(KIND=C_INT) :: OPTIONS( * )
    END FUNCTION VE_METIS_SETDEFAULTOPTIONS
  END INTERFACE

  WRITE(6,'(A)')'-----------------> Performing VH call of METIS_SETDEFAULTOPTIONS'
  IERR(1) =  VE_METIS_SETDEFAULTOPTIONS( OPTIONS )

  VE_FMETIS_SETDEFAULTOPTIONS = IERR(1)

END FUNCTION VE_FMETIS_SETDEFAULTOPTIONS

!
! -----------------------------------------------------------------------
! Interface to C metis_nodend
! -----------------------------------------------------------------------
!

FUNCTION VE_FMETIS_NODEND( NVTXS , XADJ , ADJNCY , VWGT , OPTIONS , PERM , IPERM )
  USE ISO_C_BINDING, ONLY : C_INT
  ! INTERFACE TO C ROUTINE
  INTEGER(KIND=C_INT) :: OPTIONS( * ), IERR(1)
  INTEGER(KIND=C_INT) :: NVTXS( 1 ) , XADJ( * )
  INTEGER(KIND=C_INT) :: ADJNCY( * ), VWGT( * )
  INTEGER(KIND=C_INT) :: PERM( * )  , IPERM( * )
  INTERFACE
    INTEGER(KIND=C_INT) FUNCTION VE_METIS_NODEND( NVTXS , XADJ , ADJNCY , VWGT , OPTIONS , PERM , IPERM ) &
                                 BIND(C, NAME='VE_Metis_nodend')
      USE ISO_C_BINDING, ONLY : C_INT
      INTEGER(KIND=C_INT) :: OPTIONS( * )
      INTEGER(KIND=C_INT) :: NVTXS( 1 ), XADJ( * )
      INTEGER(KIND=C_INT) :: ADJNCY( * ), VWGT( * )
      INTEGER(KIND=C_INT) :: PERM( * ), IPERM( * )
    END FUNCTION VE_METIS_NODEND
  END INTERFACE

  WRITE(6,'(A)')'-----------------> Performing VH call of METIS_NODEND'
  IERR(1) = VE_METIS_NODEND( NVTXS , XADJ , ADJNCY , VWGT , OPTIONS , PERM , IPERM )

  VE_FMETIS_NODEND = IERR(1)

END FUNCTION VE_FMETIS_NODEND

