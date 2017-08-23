# Copyright (C) 2010,2011 Stefan Vigerske
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: Makefile.am 928 2010-09-14 12:18:54Z stefan $

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST =

lib_LTLIBRARIES = libcoinsoplex.la

libcoinsoplex_la_SOURCES = \
  soplex/src/changesoplex.cpp \
  soplex/src/clufactor.cpp \
  soplex/src/didxset.cpp \
  soplex/src/enter.cpp \
  soplex/src/gzstream.cpp \
  soplex/src/idxset.cpp \
  soplex/src/leave.cpp \
  soplex/src/mpsinput.cpp \
  soplex/src/nameset.cpp \
  soplex/src/rational.cpp \
  soplex/src/slufactor.cpp \
  soplex/src/solverational.cpp \
  soplex/src/solvereal.cpp \
  soplex/src/soplex.cpp \
  soplex/src/soplexlegacy.cpp \
  soplex/src/spxautopr.cpp \
  soplex/src/spxbasis.cpp \
  soplex/src/spxboundflippingrt.cpp \
  soplex/src/spxbounds.cpp \
  soplex/src/spxchangebasis.cpp \
  soplex/src/spxdantzigpr.cpp \
  soplex/src/spxdefaultrt.cpp \
  soplex/src/spxdefines.cpp \
  soplex/src/spxdesc.cpp \
  soplex/src/spxdevexpr.cpp \
  soplex/src/spxequilisc.cpp \
  soplex/src/spxfastrt.cpp \
  soplex/src/spxfileio.cpp \
  soplex/src/spxgeometsc.cpp \
  soplex/src/spxgithash.cpp \
  soplex/src/spxharrisrt.cpp \
  soplex/src/spxhybridpr.cpp \
  soplex/src/spxid.cpp \
  soplex/src/spxlpbase_rational.cpp \
  soplex/src/spxlpbase_real.cpp \
  soplex/src/spxmainsm.cpp \
  soplex/src/spxout.cpp \
  soplex/src/spxparmultpr.cpp \
  soplex/src/spxquality.cpp \
  soplex/src/spxscaler.cpp \
  soplex/src/spxshift.cpp \
  soplex/src/spxsolver.cpp \
  soplex/src/spxsolve.cpp \
  soplex/src/spxstarter.cpp \
  soplex/src/spxsteeppr.cpp \
  soplex/src/spxsumst.cpp \
  soplex/src/spxvecs.cpp \
  soplex/src/spxvectorst.cpp \
  soplex/src/spxweightpr.cpp \
  soplex/src/spxweightst.cpp \
  soplex/src/spxwritestate.cpp \
  soplex/src/statistics.cpp \
  soplex/src/unitvector.cpp \
  soplex/src/updatevector.cpp

if TIMER_IS_USERTIMER
libcoinsoplex_la_SOURCES += soplex/src/usertimer.cpp soplex/src/wallclocktimer.cpp
else
libcoinsoplex_la_SOURCES += soplex/src/timer.cpp
endif

AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/soplex/src`

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)`

libcoinsoplex_la_LDFLAGS = $(LT_LDFLAGS)

thirdpartyincludedir = $(includedir)/coin/ThirdParty
thirdpartyinclude_HEADERS = \
  $(srcdir)/soplex/src/array.h \
  $(srcdir)/soplex/src/basevectors.h \
  $(srcdir)/soplex/src/classarray.h \
  $(srcdir)/soplex/src/clufactor.h \
  $(srcdir)/soplex/src/cring.h \
  $(srcdir)/soplex/src/dataarray.h \
  $(srcdir)/soplex/src/datahashtable.h \
  $(srcdir)/soplex/src/datakey.h \
  $(srcdir)/soplex/src/dataset.h \
  $(srcdir)/soplex/src/didxset.h \
  $(srcdir)/soplex/src/dsvectorbase.h \
  $(srcdir)/soplex/src/dsvector.h \
  $(srcdir)/soplex/src/dvectorbase.h \
  $(srcdir)/soplex/src/dvector.h \
  $(srcdir)/soplex/src/exceptions.h \
  $(srcdir)/soplex/src/gzstream.h \
  $(srcdir)/soplex/src/idlist.h \
  $(srcdir)/soplex/src/idxset.h \
  $(srcdir)/soplex/src/islist.h \
  $(srcdir)/soplex/src/lpcolbase.h \
  $(srcdir)/soplex/src/lpcol.h \
  $(srcdir)/soplex/src/lpcolsetbase.h \
  $(srcdir)/soplex/src/lpcolset.h \
  $(srcdir)/soplex/src/lprowbase.h \
  $(srcdir)/soplex/src/lprow.h \
  $(srcdir)/soplex/src/lprowsetbase.h \
  $(srcdir)/soplex/src/lprowset.h \
  $(srcdir)/soplex/src/mpsinput.h \
  $(srcdir)/soplex/src/nameset.h \
  $(srcdir)/soplex/src/random.h \
  $(srcdir)/soplex/src/rational.h \
  $(srcdir)/soplex/src/slinsolver.h \
  $(srcdir)/soplex/src/slufactor.h \
  $(srcdir)/soplex/src/solbase.h \
  $(srcdir)/soplex/src/sol.h \
  $(srcdir)/soplex/src/soplex.h \
  $(srcdir)/soplex/src/soplexlegacy.h \
  $(srcdir)/soplex/src/sorter.h \
  $(srcdir)/soplex/src/spxalloc.h \
  $(srcdir)/soplex/src/spxautopr.h \
  $(srcdir)/soplex/src/spxbasis.h \
  $(srcdir)/soplex/src/spxboundflippingrt.h \
  $(srcdir)/soplex/src/spxdantzigpr.h \
  $(srcdir)/soplex/src/spxdefaultrt.h \
  $(srcdir)/soplex/src/spxdefines.h \
  $(srcdir)/soplex/src/spxdevexpr.h \
  $(srcdir)/soplex/src/spxequilisc.h \
  $(srcdir)/soplex/src/spxfastrt.h \
  $(srcdir)/soplex/src/spxfileio.h \
  $(srcdir)/soplex/src/spxgeometsc.h \
  $(srcdir)/soplex/src/spxgithash.h \
  $(srcdir)/soplex/src/spxharrisrt.h \
  $(srcdir)/soplex/src/spxhybridpr.h \
  $(srcdir)/soplex/src/spxid.h \
  $(srcdir)/soplex/src/spxlpbase.h \
  $(srcdir)/soplex/src/spxlp.h \
  $(srcdir)/soplex/src/spxmainsm.h \
  $(srcdir)/soplex/src/spxout.h \
  $(srcdir)/soplex/src/spxparmultpr.h \
  $(srcdir)/soplex/src/spxpricer.h \
  $(srcdir)/soplex/src/spxratiotester.h \
  $(srcdir)/soplex/src/spxscaler.h \
  $(srcdir)/soplex/src/spxsimplifier.h \
  $(srcdir)/soplex/src/spxsolver.h \
  $(srcdir)/soplex/src/spxstarter.h \
  $(srcdir)/soplex/src/spxsteepexpr.h \
  $(srcdir)/soplex/src/spxsteeppr.h \
  $(srcdir)/soplex/src/spxsumst.h \
  $(srcdir)/soplex/src/spxvectorst.h \
  $(srcdir)/soplex/src/spxweightpr.h \
  $(srcdir)/soplex/src/spxweightst.h \
  $(srcdir)/soplex/src/ssvectorbase.h \
  $(srcdir)/soplex/src/ssvector.h \
  $(srcdir)/soplex/src/statistics.h \
  $(srcdir)/soplex/src/svectorbase.h \
  $(srcdir)/soplex/src/svector.h \
  $(srcdir)/soplex/src/svsetbase.h \
  $(srcdir)/soplex/src/svset.h \
  $(srcdir)/soplex/src/timer.h \
  $(srcdir)/soplex/src/unitvector.h \
  $(srcdir)/soplex/src/updatevector.h \
  $(srcdir)/soplex/src/vectorbase.h \
  $(srcdir)/soplex/src/vector.h \
  $(srcdir)/soplex/src/git_hash.cpp

if TIMER_IS_USERTIMER
thirdpartyinclude_HEADERS += $(srcdir)/soplex/src/timerfactory.h \
  $(srcdir)/soplex/src/notimer.h \
  $(srcdir)/soplex/src/wallclocktimer.h \
  $(srcdir)/soplex/src/usertimer.h
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
