# $Id: Makefile.am 765 2011-08-14 11:12:19Z stefan $
#
# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libCouenneTwoImplied.la

# List all source files for this library, including headers
libCouenneTwoImplied_la_SOURCES = \
	TwoImpliedConstructors.cpp \
	TwoImpliedIneqCombine.cpp \
	TwoImpliedGenCuts.cpp 

# This is for libtool
libCouenneTwoImplied_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/../..` \
	-I`$(CYGPATH_W) $(srcdir)/../../expression` \
	-I`$(CYGPATH_W) $(srcdir)/../../branch` \
	-I`$(CYGPATH_W) $(srcdir)/../../problem` \
	-I`$(CYGPATH_W) $(srcdir)/../../bound_tightening` \
	$(COUENNELIB_CFLAGS)

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/src
