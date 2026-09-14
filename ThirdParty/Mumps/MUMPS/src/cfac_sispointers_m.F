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
       MODULE CMUMPS_FAC_S_IS_POINTERS_M
C      ----------------------------------
C      This module defines a type used in 
C      CMUMPS_FAC_DRIVER and CMUMPS_FAC_B
C      ----------------------------------
       TYPE CMUMPS_S_IS_POINTERS_T
          COMPLEX, POINTER, DIMENSION(:) :: A
          INTEGER, POINTER, DIMENSION(:) :: IW
       END TYPE CMUMPS_S_IS_POINTERS_T
       END MODULE CMUMPS_FAC_S_IS_POINTERS_M
