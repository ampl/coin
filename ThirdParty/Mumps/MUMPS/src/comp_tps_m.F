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
      MODULE CMUMPS_TPS_M
      TYPE CMUMPS_TPS_T
          COMPLEX, DIMENSION(:), POINTER :: A
      END TYPE CMUMPS_TPS_T
      END MODULE CMUMPS_TPS_M
      SUBROUTINE CMUMPS_TPS_M_RETURN()
      RETURN
      END SUBROUTINE CMUMPS_TPS_M_RETURN
