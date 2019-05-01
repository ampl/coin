# Copyright (C) 2010,2011 Stefan Vigerske
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: Makefile.am 928 2010-09-14 12:18:54Z stefan $

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST =

lib_LTLIBRARIES = libcoinsoplex.la

#ls -1 soplex/src/soplex/*cpp | sort | awk '{printf("  %s \\\n", $0)}'
libcoinsoplex_la_SOURCES = \
  soplex/src/soplex.cpp \
  soplex/src/soplex/changesoplex.cpp \
  soplex/src/soplex/clufactor.cpp \
  soplex/src/soplex/clufactor_rational.cpp \
  soplex/src/soplex/didxset.cpp \
  soplex/src/soplex/enter.cpp \
  soplex/src/soplex/gzstream.cpp \
  soplex/src/soplex/idxset.cpp \
  soplex/src/soplex/leave.cpp \
  soplex/src/soplex/mpsinput.cpp \
  soplex/src/soplex/nameset.cpp \
  soplex/src/soplex/rational.cpp \
  soplex/src/soplex/ratrecon.cpp \
  soplex/src/soplex/slufactor.cpp \
  soplex/src/soplex/slufactor_rational.cpp \
  soplex/src/soplex/solvedbds.cpp \
  soplex/src/soplex/solverational.cpp \
  soplex/src/soplex/solvereal.cpp \
  soplex/src/soplex/spxautopr.cpp \
  soplex/src/soplex/spxbasis.cpp \
  soplex/src/soplex/spxboundflippingrt.cpp \
  soplex/src/soplex/spxbounds.cpp \
  soplex/src/soplex/spxchangebasis.cpp \
  soplex/src/soplex/spxdantzigpr.cpp \
  soplex/src/soplex/spxdefaultrt.cpp \
  soplex/src/soplex/spxdefines.cpp \
  soplex/src/soplex/spxdesc.cpp \
  soplex/src/soplex/spxdevexpr.cpp \
  soplex/src/soplex/spxequilisc.cpp \
  soplex/src/soplex/spxfastrt.cpp \
  soplex/src/soplex/spxfileio.cpp \
  soplex/src/soplex/spxgeometsc.cpp \
  soplex/src/soplex/spxgithash.cpp \
  soplex/src/soplex/spxharrisrt.cpp \
  soplex/src/soplex/spxhybridpr.cpp \
  soplex/src/soplex/spxid.cpp \
  soplex/src/soplex/spxleastsqsc.cpp \
  soplex/src/soplex/spxlpbase_rational.cpp \
  soplex/src/soplex/spxlpbase_real.cpp \
  soplex/src/soplex/spxmainsm.cpp \
  soplex/src/soplex/spxout.cpp \
  soplex/src/soplex/spxparmultpr.cpp \
  soplex/src/soplex/spxquality.cpp \
  soplex/src/soplex/spxscaler.cpp \
  soplex/src/soplex/spxshift.cpp \
  soplex/src/soplex/spxsolve.cpp \
  soplex/src/soplex/spxsolver.cpp \
  soplex/src/soplex/spxstarter.cpp \
  soplex/src/soplex/spxsteeppr.cpp \
  soplex/src/soplex/spxsumst.cpp \
  soplex/src/soplex/spxvecs.cpp \
  soplex/src/soplex/spxvectorst.cpp \
  soplex/src/soplex/spxweightpr.cpp \
  soplex/src/soplex/spxweightst.cpp \
  soplex/src/soplex/spxwritestate.cpp \
  soplex/src/soplex/statistics.cpp \
  soplex/src/soplex/updatevector.cpp \
  soplex/src/soplex/usertimer.cpp \
  soplex/src/soplex/wallclocktimer.cpp

AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/soplex/src`

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)`

libcoinsoplex_la_LDFLAGS = $(LT_LDFLAGS)

# List all additionally required libraries
if DEPENDENCY_LINKING
libcoinsoplex_la_LIBADD = $(SOPLEX_LIBS)
endif

thirdpartyincludedir = $(includedir)/coin/ThirdParty
thirdpartyinclude_HEADERS = \
  soplex/src/soplex.h

thirdpartyspxincludedir = $(includedir)/coin/ThirdParty/soplex
thirdpartyspxinclude_HEADERS = \
  soplex/src/soplex/array.h \
  soplex/src/soplex/basevectors.h \
  soplex/src/soplex/classarray.h \
  soplex/src/soplex/clufactor.h \
  soplex/src/soplex/clufactor_rational.h \
  soplex/src/soplex/cring.h \
  soplex/src/soplex/dataarray.h \
  soplex/src/soplex/datahashtable.h \
  soplex/src/soplex/datakey.h \
  soplex/src/soplex/dataset.h \
  soplex/src/soplex/didxset.h \
  soplex/src/soplex/dsvectorbase.h \
  soplex/src/soplex/dsvector.h \
  soplex/src/soplex/dvectorbase.h \
  soplex/src/soplex/dvector.h \
  soplex/src/soplex/exceptions.h \
  soplex/src/soplex/gzstream.h \
  soplex/src/soplex/idlist.h \
  soplex/src/soplex/idxset.h \
  soplex/src/soplex/islist.h \
  soplex/src/soplex/lpcolbase.h \
  soplex/src/soplex/lpcol.h \
  soplex/src/soplex/lpcolsetbase.h \
  soplex/src/soplex/lpcolset.h \
  soplex/src/soplex/lprowbase.h \
  soplex/src/soplex/lprow.h \
  soplex/src/soplex/lprowsetbase.h \
  soplex/src/soplex/lprowset.h \
  soplex/src/soplex/mpsinput.h \
  soplex/src/soplex/nameset.h \
  soplex/src/soplex/notimer.h \
  soplex/src/soplex/random.h \
  soplex/src/soplex/rational.h \
  soplex/src/soplex/ratrecon.h \
  soplex/src/soplex/slinsolver.h \
  soplex/src/soplex/slinsolver_rational.h \
  soplex/src/soplex/slufactor.h \
  soplex/src/soplex/slufactor_rational.h \
  soplex/src/soplex/solbase.h \
  soplex/src/soplex/sol.h \
  soplex/src/soplex/sorter.h \
  soplex/src/soplex/spxalloc.h \
  soplex/src/soplex/spxautopr.h \
  soplex/src/soplex/spxbasis.h \
  soplex/src/soplex/spxboundflippingrt.h \
  soplex/src/soplex/spxdantzigpr.h \
  soplex/src/soplex/spxdefaultrt.h \
  soplex/src/soplex/spxdefines.h \
  soplex/src/soplex/spxdevexpr.h \
  soplex/src/soplex/spxequilisc.h \
  soplex/src/soplex/spxfastrt.h \
  soplex/src/soplex/spxfileio.h \
  soplex/src/soplex/spxgeometsc.h \
  soplex/src/soplex/spxgithash.h \
  soplex/src/soplex/spxharrisrt.h \
  soplex/src/soplex/spxhybridpr.h \
  soplex/src/soplex/spxid.h \
  soplex/src/soplex/spxleastsqsc.h \
  soplex/src/soplex/spxlpbase.h \
  soplex/src/soplex/spxlp.h \
  soplex/src/soplex/spxmainsm.h \
  soplex/src/soplex/spxout.h \
  soplex/src/soplex/spxparmultpr.h \
  soplex/src/soplex/spxpricer.h \
  soplex/src/soplex/spxratiotester.h \
  soplex/src/soplex/spxscaler.h \
  soplex/src/soplex/spxsimplifier.h \
  soplex/src/soplex/spxsolver.h \
  soplex/src/soplex/spxstarter.h \
  soplex/src/soplex/spxsteepexpr.h \
  soplex/src/soplex/spxsteeppr.h \
  soplex/src/soplex/spxsumst.h \
  soplex/src/soplex/spxvectorst.h \
  soplex/src/soplex/spxweightpr.h \
  soplex/src/soplex/spxweightst.h \
  soplex/src/soplex/ssvectorbase.h \
  soplex/src/soplex/ssvector.h \
  soplex/src/soplex/statistics.h \
  soplex/src/soplex/svectorbase.h \
  soplex/src/soplex/svector.h \
  soplex/src/soplex/svsetbase.h \
  soplex/src/soplex/svset.h \
  soplex/src/soplex/timerfactory.h \
  soplex/src/soplex/timer.h \
  soplex/src/soplex/unitvectorbase.h \
  soplex/src/soplex/unitvector.h \
  soplex/src/soplex/updatevector.h \
  soplex/src/soplex/usertimer.h \
  soplex/src/soplex/vectorbase.h \
  soplex/src/soplex/vector.h \
  soplex/src/soplex/wallclocktimer.h \
  soplex/src/soplex/git_hash.cpp

# soplex/stablesum.h was added with SoPlex 4.0.1
if HAVE_STABLESUM_H
  thirdpartyspxinclude_HEADERS += soplex/src/soplex/stablesum.h
endif

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = coinsoplex.pc

test:
	@echo "No test available for SoPlex."

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES = $(MAKEOKFILE)

include BuildTools/Makemain.inc
