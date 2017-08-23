# Copyright (C) 2003, 2006, 2009 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: Makefile.am 2080 2011-01-08 17:58:57Z stefan $

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST =

lib_LTLIBRARIES = libfiltersqp.la

libfiltersqp_la_SOURCES = \
	source/src/auxil.f \
	source/src/bqpd.f \
	source/src/filteraux.f \
	source/src/filter.f \
	source/src/QPsolves.f \
	source/src/scaling.f \
	source/src/sparseA.f \
	source/src/sparseL.f \
	source/src/util.f \
	source/src/wdotd.f

libfiltersqp_la_LDFLAGS = $(LT_LDFLAGS)

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = coinfiltersqp.pc

test:
	@echo "No test available for FilterSQP."

########################################################################
#                         Maintainer Stuff                             #
########################################################################

CLEANFILES =

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include BuildTools/Makemain.inc
