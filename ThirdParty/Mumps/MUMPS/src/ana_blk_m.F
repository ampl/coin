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
      MODULE MUMPS_ANA_BLK_M
      TYPE COL_LMATRIX_T
       INTEGER :: NBINCOL   
       INTEGER, POINTER :: IRN(:) => null()  
      END TYPE COL_LMATRIX_T
      TYPE LMATRIX_T
       INTEGER    :: NBCOL             = 0
       LOGICAL    :: CONTIGUOUS_COL    = .FALSE.
       INTEGER    :: NBCOL_LOC         = 0
       INTEGER    :: FIRST             = 0
       INTEGER(8) :: NZL               = 0_8
       TYPE(COL_LMATRIX_T), POINTER :: COL(:) => null()
      END TYPE LMATRIX_T
      TYPE COMPACT_GRAPH_T
       INTEGER(8) :: NZG, SIZEADJALLOCATED
       INTEGER    :: NG      =0_8
       INTEGER    :: FIRST   =0
       INTEGER    :: LAST    =0
       INTEGER(8), POINTER :: IPE(:) => null()  
       INTEGER, POINTER :: ADJ(:) => null()  
      END TYPE COMPACT_GRAPH_T
      END MODULE MUMPS_ANA_BLK_M
