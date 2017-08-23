# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: configure.ac 3527 2015-08-09 14:32:55Z stefan $

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([ThirdPartyFilterSQP],[1.3.4],[http://projects.coin-or.org/BuildTools/newticket])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(source/src/filter.f)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

AC_COIN_PROJECTDIR_INIT(FilterSQP,4:4:3)

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(FilterSQP)

# Get the name of the Fortran compiler and appropriate compiler options
AC_COIN_PROG_F77

# Get the linker flags (FLIBS) that are required for linking against this blas library
AC_COIN_F77_WRAPPERS

# Add FLIBS to FILTERSQP_PCLIBS, so that they get into the .pc files
FILTERSQP_PCLIBS="$FLIBS"
AC_SUBST(FILTERSQP_PCLIBS)

# Check if pkg-config is available (to setup addlibs file)
AC_COIN_HAS_PKGCONFIG

# Initialize automake
AC_COIN_INIT_AUTOMAKE

# Initialize automake
AC_COIN_INIT_AUTO_TOOLS

AC_CONFIG_FILES([Makefile coinfiltersqp.pc coinfiltersqp-uninstalled.pc])

AC_COIN_FINALIZE

##############################################################################
##                         LIBRARY FOR MINLP/FILTER                          #
##############################################################################
#
#if test -z "$FSQPMAKEFLAGS" ; then
#  FSQPMAKEFLAGS=
#
#  AC_MSG_CHECKING([whether $am_make is the GNU make])
#  rm -rf conftest.mak
#  ($am_make -v -f conftest.mak > conftest.grp) 2>&5
#  if $EGREP GNU conftest.grp >&5 ; then
#    AC_MSG_RESULT([yes])
#    FSQPMAKEFLAGS="-j 1 $FSQPMAKEFLAGS"
#  else
#    AC_MSG_RESULT([no])
#  fi
#  rm -f conftest.mak conftest.grp
#fi
#
#if test x"$AR" = x; then
#  AR=ar
#fi
#if test x"$ARFLAGS" = x; then
#  ARFLAGS=r
#fi
#
#FSQPMAKEFILECMDS="sed -e 's%ar r%$AR $ARFLAGS%' makefile.u "
#FSQPMAKEFILECMDS="$FSQPMAKEFILECMDS > makefile"
#
#AC_ARG_VAR(FSQPMAKEFLAGS,[Makefile flags necessary for compiling the filterSQP library])
#AC_SUBST(FSQPMAKEFILECMDS)
#
#AC_CONFIG_FILES([Makefile])
#
#AC_COIN_FINALIZE
