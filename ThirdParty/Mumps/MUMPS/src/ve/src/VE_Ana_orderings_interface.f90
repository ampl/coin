!
! -----------------------------------------------------------------------
! 
! NEC Vector Host MUMPS Library Interface
! 
! -----------------------------------------------------------------------
!

!
! -----------------------------------------------------------------------
! Interface to MUMPS_ANA_H
! -----------------------------------------------------------------------
!

 SUBROUTINE VE_FMUMPS_ANA_H(                                          &
                  TOTEL  , COMPUTE_PERM  , N       , IWLEN  , PE    , &
                  PFREE  , LEN           , IW      , NV     , ELEN  , &
                  LAST   , NCMPA         , DEGREE  , HEAD   , NEXT  , &
                  W      , PARENT        )
  USE ISO_C_BINDING, ONLY : C_INT,C_BOOL,C_INT64_T
  IMPLICIT NONE
  ! INTERFACE TO C ROUTINE
  !     Input not modified
        INTEGER(KIND=C_INT), INTENT(IN)      :: TOTEL, N
        INTEGER(KIND=C_INT64_T), INTENT(IN)  :: IWLEN
        LOGICAL(KIND=C_BOOL), INTENT(IN)     :: COMPUTE_PERM
  !     Input undefined on output 
        INTEGER(KIND=C_INT), INTENT(INOUT) :: LEN(N), IW(IWLEN)
  ! 
  !     Output only 
        INTEGER(KIND=C_INT), INTENT(OUT) :: NCMPA(1)
        INTEGER(KIND=C_INT), INTENT(OUT) :: ELEN(N), LAST(N), PARENT(N)
  ! 
  !     Input/output
        INTEGER(KIND=C_INT64_T), INTENT(INOUT) :: PFREE(1)
        INTEGER(KIND=C_INT64_T), INTENT(INOUT) :: PE(N)
  !     NV also meaningful as input to encode compressed graphs
        INTEGER(KIND=C_INT), INTENT(INOUT) :: NV(N)
  ! 
  !     Internal Workspace only
        INTEGER(KIND=C_INT) :: NEXT(N), DEGREE(N), HEAD(TOTEL), W(N)
  INTERFACE
    SUBROUTINE VE_MUMPS_ANA_H(                                        &
                  TOTEL  , COMPUTE_PERM  , N       , IWLEN  , PE    , &
                  PFREE  , LEN           , IW      , NV     , ELEN  , &
                  LAST   , NCMPA         , DEGREE  , HEAD   , NEXT  , &
                  W      , PARENT        )                            &
                  BIND(C, NAME='VE_Mumps_ana_h')
      USE ISO_C_BINDING, ONLY : C_INT,C_BOOL,C_INT64_T
        IMPLICIT NONE
  !     Input not modified
        INTEGER(KIND=C_INT)    , VALUE, INTENT(IN) :: TOTEL, N
        INTEGER(KIND=C_INT64_T), VALUE, INTENT(IN) :: IWLEN
        LOGICAL(KIND=C_BOOL)   , VALUE, INTENT(IN) :: COMPUTE_PERM
  !     Input undefined on output 
        INTEGER(KIND=C_INT), INTENT(INOUT) :: LEN(N), IW(IWLEN)
  ! 
  !     Output only 
        INTEGER(KIND=C_INT), INTENT(OUT) :: NCMPA(1)
        INTEGER(KIND=C_INT), INTENT(OUT) :: ELEN(N), LAST(N), PARENT(N)
  ! 
  !     Input/output
        INTEGER(KIND=C_INT64_T), INTENT(INOUT) :: PFREE(1)
        INTEGER(KIND=C_INT64_T), INTENT(INOUT) :: PE(N)
  !     NV also meaningful as input to encode compressed graphs
        INTEGER(KIND=C_INT), INTENT(INOUT) :: NV(N)
  ! 
  !     Internal Workspace only
        INTEGER(KIND=C_INT) :: NEXT(N), DEGREE(N), HEAD(TOTEL), W(N)
    END SUBROUTINE VE_MUMPS_ANA_H
  END INTERFACE

  WRITE(6,'(A)')'-----------------> Performing VH call of MUMPS_ANA_H'
  CALL VE_MUMPS_ANA_H(                                                &
                  TOTEL  , COMPUTE_PERM  , N       , IWLEN  , PE    , &
                  PFREE  , LEN           , IW      , NV     , ELEN  , &
                  LAST   , NCMPA         , DEGREE  , HEAD   , NEXT  , &
                  W      , PARENT        )
END SUBROUTINE VE_FMUMPS_ANA_H

!
! -----------------------------------------------------------------------
! Interface to MUMPS_SYMQAMD
! -----------------------------------------------------------------------
!

 SUBROUTINE VE_FMUMPS_SYMQAMD(                                        &
          THRESH  , NDENSE  , N      , TOTEL          , IWLEN       , &
          PE      , PFREE   , LEN    , IW             , NV          , &
          ELEN    , LAST    , NCMPA  , DEGREE         , HEAD        , &
          NEXT    , W       , PERM   , LISTVAR_SCHUR  , SIZE_SCHUR  , &
          AGG6    , PARENT  )
  USE ISO_C_BINDING, ONLY : C_INT,C_BOOL,C_INT64_T
  IMPLICIT NONE
  ! INTERFACE TO C ROUTINE
  !     Input not modified
        INTEGER(KIND=C_INT)    , INTENT(IN) :: N, TOTEL, SIZE_SCHUR
        LOGICAL(KIND=C_BOOL)   , INTENT(IN) :: AGG6
        INTEGER(KIND=C_INT)    , INTENT(IN) :: THRESH
        INTEGER(KIND=C_INT64_T), INTENT(IN) :: IWLEN
        INTEGER(KIND=C_INT)    , INTENT(IN) :: LISTVAR_SCHUR(max(1,SIZE_SCHUR))
  !     Input undefined on output 
        INTEGER, INTENT(INOUT)  :: LEN(N), IW(IWLEN)
  ! 
  !     Output only 
        INTEGER(KIND=C_INT), INTENT(OUT)   :: NCMPA(1)
        INTEGER(KIND=C_INT), INTENT(OUT)   :: ELEN(N), LAST(TOTEL), PARENT(N)
  ! 
  !     Input/output
        INTEGER(KIND=C_INT), INTENT(INOUT)     :: NV(N)
        INTEGER(KIND=C_INT64_T), INTENT(INOUT) :: PFREE(1)
        INTEGER(KIND=C_INT64_T), INTENT(INOUT) :: PE(N)
        INTEGER(KIND=C_INT), INTENT(INOUT)     :: PERM(N)
  ! 
  !     Internal Workspace only
        INTEGER(KIND=C_INT), INTENT(OUT) :: NDENSE(N), DEGREE(N)
        INTEGER(KIND=C_INT), INTENT(OUT) :: HEAD(TOTEL), NEXT(N), W(N)
  INTERFACE
    SUBROUTINE VE_MUMPS_SYMQAMD(                                              &
                  THRESH  , NDENSE  , N      , TOTEL          , IWLEN       , &
                  PE      , PFREE   , LEN    , IW             , NV          , &
                  ELEN    , LAST    , NCMPA  , DEGREE         , HEAD        , &
                  NEXT    , W       , PERM   , LISTVAR_SCHUR  , SIZE_SCHUR  , &
                  AGG6    , PARENT  )                                         &
                  BIND(C, NAME='VE_Mumps_symqamd')
      USE ISO_C_BINDING, ONLY : C_INT,C_BOOL,C_INT64_T
        IMPLICIT NONE
  !     Input not modified
        INTEGER(KIND=C_INT)    , VALUE , INTENT(IN) :: N, TOTEL, SIZE_SCHUR
        LOGICAL(KIND=C_BOOL)   , VALUE , INTENT(IN) :: AGG6
        INTEGER(KIND=C_INT)    , VALUE , INTENT(IN) :: THRESH
        INTEGER(KIND=C_INT64_T), VALUE , INTENT(IN) :: IWLEN
        INTEGER(KIND=C_INT)    , INTENT(IN) :: LISTVAR_SCHUR(max(1,SIZE_SCHUR))
  !     Input undefined on output 
        INTEGER, INTENT(INOUT)  :: LEN(N), IW(IWLEN)
  ! 
  !     Output only 
        INTEGER(KIND=C_INT), INTENT(OUT)   :: NCMPA(1)
        INTEGER(KIND=C_INT), INTENT(OUT)   :: ELEN(N), LAST(TOTEL), PARENT(N)
  ! 
  !     Input/output
        INTEGER(KIND=C_INT), INTENT(INOUT)     :: NV(N)
        INTEGER(KIND=C_INT64_T), INTENT(INOUT) :: PFREE(1)
        INTEGER(KIND=C_INT64_T), INTENT(INOUT) :: PE(N)
        INTEGER(KIND=C_INT), INTENT(INOUT)     :: PERM(N)
  ! 
  !     Internal Workspace only
        INTEGER(KIND=C_INT), INTENT(OUT) :: NDENSE(N), DEGREE(N)
        INTEGER(KIND=C_INT), INTENT(OUT) :: HEAD(TOTEL), NEXT(N), W(N)
    END SUBROUTINE VE_MUMPS_SYMQAMD
  END INTERFACE

  WRITE(6,'(A)')'-----------------> Performing VH call of MUMPS_SYMQAMD'
  CALL VE_MUMPS_SYMQAMD(                                                      &
                  THRESH  , NDENSE  , N      , TOTEL          , IWLEN       , &
                  PE      , PFREE   , LEN    , IW             , NV          , &
                  ELEN    , LAST    , NCMPA  , DEGREE         , HEAD        , &
                  NEXT    , W       , PERM   , LISTVAR_SCHUR  , SIZE_SCHUR  , &
                  AGG6    , PARENT  )
END SUBROUTINE VE_FMUMPS_SYMQAMD

!
! -----------------------------------------------------------------------
! Interface to MUMPS_WRAP_GINP94
! -----------------------------------------------------------------------
!

 SUBROUTINE VE_FMUMPS_WRAP_GINP94(                                         &
      N             , IPE     , IW             , LIW8        , PERM      , &
      SizeOfBlocks  , KEEP60  , LISTVAR_SCHUR  , SIZE_SCHUR  , KEEP378   , &
      COLCOUNT      , PARENT  , PORDER         , IWTMP1      , IWTMP2    , &
      IWTMP3        , IWTMP4  , IWTMP5         , INFO        )
  USE ISO_C_BINDING, ONLY : C_INT,C_INT64_T
  IMPLICIT NONE
  ! INTERFACE TO C ROUTINE
  !     Input not modified
        INTEGER(KIND=C_INT)    , INTENT(IN) :: N, KEEP60, SIZE_SCHUR, KEEP378
        INTEGER(KIND=C_INT)    , INTENT(IN) :: SizeOfBlocks(N)
        INTEGER(KIND=C_INT)    , INTENT(IN) :: LISTVAR_SCHUR(SIZE_SCHUR)
        INTEGER(KIND=C_INT64_T), INTENT(IN) :: LIW8,IPE(N+1)
        INTEGER(KIND=C_INT)    , INTENT(IN) :: IW(LIW8)
  ! 
  !     Output only 
        INTEGER(KIND=C_INT), INTENT(OUT)   :: COLCOUNT(N),PARENT(N)
        INTEGER(KIND=C_INT), INTENT(OUT)   :: PORDER(N), IWTMP1(N), IWTMP2(N)
        INTEGER(KIND=C_INT), INTENT(OUT)   :: IWTMP3(N), IWTMP4(N), IWTMP5(N)
  ! 
  !     Input/output
        INTEGER(KIND=C_INT), INTENT(INOUT)     :: PERM(N)
        INTEGER(KIND=C_INT), INTENT(INOUT)     :: INFO(2)
  INTERFACE
    SUBROUTINE VE_MUMPS_WRAP_GINP94(                                          &
         N             , IPE     , IW             , LIW8        , PERM      , &
         SizeOfBlocks  , KEEP60  , LISTVAR_SCHUR  , SIZE_SCHUR  , KEEP378   , &
         COLCOUNT      , PARENT  , PORDER         , IWTMP1      , IWTMP2    , &
         IWTMP3        , IWTMP4  , IWTMP5         , INFO        )             &
         BIND(C, NAME='VE_Mumps_wrap_ginp94')
      USE ISO_C_BINDING, ONLY : C_INT,C_INT64_T
        IMPLICIT NONE
  !     Input not modified
        INTEGER(KIND=C_INT)    , VALUE , INTENT(IN) :: N, KEEP60, SIZE_SCHUR, KEEP378
        INTEGER(KIND=C_INT64_T), VALUE , INTENT(IN) :: LIW8
        INTEGER(KIND=C_INT)            , INTENT(IN) :: SizeOfBlocks(N)
        INTEGER(KIND=C_INT)            , INTENT(IN) :: LISTVAR_SCHUR(SIZE_SCHUR)
        INTEGER(KIND=C_INT64_T)        , INTENT(IN) :: IPE(N+1)
        INTEGER(KIND=C_INT)            , INTENT(IN) :: IW(LIW8)
  ! 
  !     Output only 
        INTEGER(KIND=C_INT), INTENT(OUT)   :: COLCOUNT(N),PARENT(N)
        INTEGER(KIND=C_INT), INTENT(OUT)   :: PORDER(N), IWTMP1(N), IWTMP2(N)
        INTEGER(KIND=C_INT), INTENT(OUT)   :: IWTMP3(N), IWTMP4(N), IWTMP5(N)
  ! 
  !     Input/output
        INTEGER(KIND=C_INT), INTENT(INOUT)     :: PERM(N)
        INTEGER(KIND=C_INT), INTENT(INOUT)     :: INFO(2)
    END SUBROUTINE VE_MUMPS_WRAP_GINP94
  END INTERFACE

  WRITE(6,'(A)')'-----------------> Performing VH call of MUMPS_WRAP_GINP94'
  CALL VE_MUMPS_WRAP_GINP94(                                                &
             N             , IPE     , IW             , LIW8        , PERM      , &
             SizeOfBlocks  , KEEP60  , LISTVAR_SCHUR  , SIZE_SCHUR  , KEEP378   , &
             COLCOUNT      , PARENT  , PORDER         , IWTMP1      , IWTMP2    , &
             IWTMP3        , IWTMP4  , IWTMP5         , INFO        )
END SUBROUTINE VE_FMUMPS_WRAP_GINP94
