# Copyright (C) 2011 Stefan Vigerske and others
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: configure.ac 937 2010-10-02 19:05:12Z stefan $

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([ThirdPartyScip],[1.8.1],[http://projects.coin-or.org/BuildTools/newticket])

AC_COPYRIGHT([
Copyright 2011-2012 Stefan Vigerske and others.
All Rights Reserved.
This file is part of the open source package COIN-OR which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(scip/src/scip/scip.h)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

AC_COIN_PROJECTDIR_INIT(SCIP,9:1:8)

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(SCIP)

# Get the name of the C compiler and appropriate compiler options
AC_COIN_PROG_CC

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Initialize automake
AC_COIN_INIT_AUTO_TOOLS

#############################################################################
#                         Check SCIP version number                         #
#############################################################################

AC_MSG_CHECKING([SCIP version])
scipver=`awk '/define SCIP_VERSION/ {print $3}' ${srcdir}/scip/src/scip/def.h`
scipsubver=`awk '/define SCIP_SUBVERSION/ {print $3}' ${srcdir}/scip/src/scip/def.h`
AC_MSG_RESULT([${scipver}.${scipsubver}])
if test $scipver -lt 600 -o $scipver -ge 610 ; then
  AC_MSG_ERROR(["Incorrect SCIP version. This build system is made for SCIP 5.0.x"])
fi

#############################################################################
#                       Check for external libraries                        #
#############################################################################

# check for some LP solvers
AC_COIN_CHECK_PACKAGE(SoPlex, coinsoplex, SCIPLib)
AC_COIN_CHECK_PACKAGE(Clp, clp, SCIPLib)
AC_COIN_CHECK_USER_LIBRARY([Cplex], [CPX],[cplex.h],   [CPXgetstat],,   [SCIPLib])
AC_COIN_CHECK_USER_LIBRARY([Mosek], [MSK],[mosek.h],   [MSK_makeenv],,  [SCIPLib])
AC_COIN_CHECK_USER_LIBRARY([Xpress],[XPR],[xprs.h],    [XPRSinit],,     [SCIPLib])
AC_COIN_CHECK_USER_LIBRARY([Gurobi],[GRB],[gurobi_c.h],[GRBloadenv],,   [SCIPLib])
AC_COIN_CHECK_USER_LIBRARY([QSopt], [QSO],[qsopt.h],   [QScreate_prob],,[SCIPLib])

AC_ARG_WITH([scip-lpsolver],
  AS_HELP_STRING([--with-scip-lpsolver],
    [specify default LP solver to use in SCIP (soplex, soplex2, clp, cplex, mosek, xpress, gurobi, qsopt, none)]),
  [if test "$withval" = no ; then
     scip_lpsolver=none
   else
     scip_lpsolver="$withval"
   fi],
  [if test "$coin_has_soplex" = yes ; then
     scip_lpsolver=soplex2
   elif test "$coin_has_clp" = yes ; then
     scip_lpsolver=clp
   elif test "$coin_has_cpx" = true ; then
     scip_lpsolver=cplex
   elif test "$coin_has_msk" = true ; then
     scip_lpsolver=mosek
   elif test "$coin_has_xpr" = true ; then
     scip_lpsolver=xpress
   elif test "$coin_has_grb" = true ; then
     scip_lpsolver=gurobi
   elif test "$coin_has_qso" = true ; then
     scip_lpsolver=qsopt
   else
     scip_lpsolver=none
   fi
  ])

case "$scip_lpsolver" in
  soplex)
    if test "$coin_has_soplex" = no ; then
      AC_MSG_ERROR(["LP solver SoPlex not available.])
    fi
    AC_DEFINE(SCIP_LPISW_LPDEFAULT, SCIP_LPISW_SOPLEX2, [Default LP solver in SCIP])
  ;;
  soplex2)
    if test "$coin_has_soplex" = no ; then
      AC_MSG_ERROR(["LP solver SoPlex not available.])
    fi
    AC_DEFINE(SCIP_LPISW_LPDEFAULT, SCIP_LPISW_SOPLEX2, [Default LP solver in SCIP])
  ;;
  clp)
    if test "$coin_has_clp" = no ; then
      AC_MSG_ERROR(["LP solver CLP not available.])
    fi
    AC_DEFINE(SCIP_LPISW_LPDEFAULT, SCIP_LPISW_CLP, [Default LP solver in SCIP])
  ;;
  cplex)
    if test "$coin_has_cpx" = false ; then
      AC_MSG_ERROR(["LP solver CPLEX not available.])
    fi
    AC_DEFINE(SCIP_LPISW_LPDEFAULT, SCIP_LPISW_CPLEX, [Default LP solver in SCIP])
  ;;
  mosek)
    if test "$coin_has_msk" = false ; then
      AC_MSG_ERROR(["LP solver Mosek not available.])
    fi
    AC_DEFINE(SCIP_LPISW_LPDEFAULT, SCIP_LPISW_MOSEK, [Default LP solver in SCIP])
  ;;
  xpress)
    if test "$coin_has_xpr" = false ; then
      AC_MSG_ERROR(["LP solver Xpress not available.])
    fi
    AC_DEFINE(SCIP_LPISW_LPDEFAULT, SCIP_LPISW_XPRESS, [Default LP solver in SCIP])
  ;;
  gurobi)
    if test "$coin_has_grb" = false ; then
      AC_MSG_ERROR(["LP solver Gurobi not available.])
    fi
    AC_DEFINE(SCIP_LPISW_LPDEFAULT, SCIP_LPISW_GUROBI, [Default LP solver in SCIP])
  ;;
  qsopt)
    if test "$coin_has_qso" = false ; then
      AC_MSG_ERROR(["LP solver QSopt not available.])
    fi
    AC_DEFINE(SCIP_LPISW_LPDEFAULT, SCIP_LPISW_QSOPT, [Default LP solver in SCIP])
  ;;
  none)
    AC_DEFINE(SCIP_LPISW_LPDEFAULT, SCIP_LPISW_NONE, [Default LP solver in SCIP])
  ;;
  *)
    AC_MSG_ERROR([LP solver $scip_lpsolver unknown.])
  ;;
esac
AC_MSG_NOTICE([default LP solver is $scip_lpsolver])

# check for NLP solvers
AC_COIN_CHECK_PACKAGE(Ipopt, ipopt, SCIPLib)

#############################################################################
#                       Check for system header and libraries               #
#############################################################################

EXTRAFLAGS="-DWITH_SCIPDEF -DNO_CONFIG_HEADER"

AC_ARG_ENABLE([scip-threadsafe],
  AS_HELP_STRING([--disable-scip-threadsafe],[do not build SCIP threadsafe (default=auto)]),
  [scip_threadsafe=$enableval],
  [scip_threadsafe=auto])  

if test $scip_threadsafe = yes || test $scip_threadsafe = auto ; then 
  AC_CHECK_LIB([pthread],[pthread_create],
               [SCIPLIB_LIBS="-lpthread $SCIPLIB_LIBS"
                SCIPLIB_PCLIBS="-lpthread $SCIPLIB_PCLIBS"
               ],
               [if test $scip_threadsafe = yes ; then
                  AC_MSG_ERROR([--scip-threadsafe selected, but -lpthreads unavailable])
                else
                  AC_MSG_NOTICE([no pthreads, SCIP will not be threadsafe])
                  EXTRAFLAGS="$EXTRAFLAGS -DNPARASCIP"
                fi
               ])
else
  EXTRAFLAGS="$EXTRAFLAGS -DNPARASCIP"
fi

AC_COIN_CHECK_GNU_ZLIB(SCIPlib)
if test $coin_has_zlib = "yes"; then
  EXTRAFLAGS="$EXTRAFLAGS -DWITH_ZLIB -DSCIP_WITH_ZLIB"  
fi

AC_COIN_CHECK_GMP(SCIPLib)
if test $coin_has_gmp = yes; then
  EXTRAFLAGS="$EXTRAFLAGS -DWITH_GMP -DSCIP_WITH_GMP"  
fi

AC_COIN_CHECK_GNU_READLINE(SCIPLib)
if test $coin_has_readline = yes; then
  EXTRAFLAGS="$EXTRAFLAGS -DWITH_READLINE -DSCIP_WITH_READLINE"
fi

AC_CHECK_MEMBER(struct sigaction.sa_handler,,
  [EXTRAFLAGS="$EXTRAFLAGS -DNO_SIGACTION -DSCIP_NO_SIGACTION"],
  [#include <signal.h>])
    
AC_CHECK_FUNC(rand_r,,
  [EXTRAFLAGS="$EXTRAFLAGS -DNO_RAND_R"])
    
AC_CHECK_FUNC(strtok_r,,
  [EXTRAFLAGS="$EXTRAFLAGS -DNO_STRTOK_R -DSCIP_NO_STRTOK_R"])

AC_CHECK_FUNC(strerror_r,,
  [EXTRAFLAGS="$EXTRAFLAGS -DNO_STRERROR_R"])
    
AC_CHECK_HEADER([fenv.h],[EXTRAFLAGS="$EXTRAFLAGS -DSCIP_ROUNDING_FE"],
  [AC_CHECK_HEADER([float.h],[EXTRAFLAGS="$EXTRAFLAGS -DSCIP_ROUNDING_MS"])])

# check if user specifies bliss
AC_LANG_PUSH(C++)
AC_COIN_CHECK_USER_LIBRARY([Bliss],[BLISS],[bliss/graph.hh],[bliss_new],,[SCIPLib])
AC_LANG_POP

# check for bliss system installation
if test $coin_has_bliss = false ; then
  AC_LANG_PUSH(C++)
  AC_CHECK_LIB(bliss, bliss_new,
    [SCIPLIB_LIBS="$SCIPLIB_LIBS -lbliss"
     SCIPLIB_PCLIBS="$SCIPLIB_PCLIBS -lbliss"
     coin_has_bliss=true
    ])
  AC_CHECK_HEADER([bliss/orbit.hh],[have_bliss=yes],[have_bliss=no])
  AC_LANG_POP
fi
AM_CONDITIONAL(COIN_HAS_BLISS, [test $coin_has_bliss = true])
AM_CONDITIONAL(COIN_HAS_BLISSINCDIR, [test -n "$BLISSINCDIR"])
AC_ARG_ENABLE([bliss-gmp],
  AS_HELP_STRING([--enable-bliss-gmp],[indicate whether to assume that Bliss was build with GMP (default=yes, if GMP available, otherwise no)]),
  [bliss_has_gmp=$enableval],
  [bliss_has_gmp=$coin_has_gmp])
if test $bliss_has_gmp = yes ;
then
  EXTRAFLAGS="$EXTRAFLAGS -DBLISS_USE_GMP"
fi

case "$CC" in
  *gcc* )
    EXTRAFLAGS="$EXTRAFLAGS -Wno-overlength-strings"
    ;;
esac

# recent autotools versions can recognize OpenMP with the following macro
# this would be useful to set TPI
# AC_OPENMP
EXTRAFLAGS="$EXTRAFLAGS -DTPI_NONE"

CFLAGS="$CFLAGS $EXTRAFLAGS"
CXXFLAGS="$CXXFLAGS $EXTRAFLAGS"

#############################################################################
#                                Finish up                                  #
#############################################################################

AC_CONFIG_FILES([Makefile coinscip.pc coinscip-uninstalled.pc])

AC_CONFIG_HEADER([config_scip.h])
 
AC_COIN_FINALIZE
