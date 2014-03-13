# Copyright (C) 2007, 2011 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author: Andreas Waechter     2007-10-08     IBM

# $Id: Makefile.in 2276 2013-05-05 12:33:44Z stefan $

## Installation of JDK
#JAVA_HOME = /usr/lib/jvm/java-7-openjdk-amd64

# NOTE: We assume here that there are no spaces in $JAVA_HOME.
#   If you are on Windows and have Java in a path like C:\Program Files\Java,
#   try setting JAVA_HOME to the DOS equivalent C:\Progra~1 or similar.

# Java compiler
JAVAC = $(JAVA_HOME)/bin/javac

# JNI header file generator
JAVAH = $(JAVA_HOME)/bin/javah -jni

# Jar utility
JAR = $(JAVA_HOME)/bin/jar

# Java runtime
JAVA = $(JAVA_HOME)/bin/java

# Java documentation utility
JAVADOC = $(JAVA_HOME)/bin/javadoc

# libtool
top_builddir = ../..
SHELL = /bin/bash
LIBTOOL = $(SHELL) $(top_builddir)/../libtool

# C++ Compiler command
CXX = g++

# C++ Compiler options
CXXFLAGS = `echo -O3 -pipe -DNDEBUG -pedantic-errors -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -Wno-unknown-pragmas -Wno-long-long   -DIPOPT_BUILD | sed -e s'|-pedantic-errors||'`

# additional C++ Compiler options for linking
LINKFLAGS = $(CXXFLAGS)  -Wl,--rpath -Wl,/home/viz/Downloads/Ipopt-3.11.7/lib -shared

# Include directories (we use the CYGPATH_W variables to allow compilation with Windows compilers)
INCL = `PKG_CONFIG_PATH=/home/viz/Downloads/Ipopt-3.11.7/lib64/pkgconfig:/home/viz/Downloads/Ipopt-3.11.7/lib/pkgconfig:/home/viz/Downloads/Ipopt-3.11.7/share/pkgconfig: pkg-config --cflags ipopt` -I`$(CYGPATH_W) $(JAVA_HOME)/include` -I`$(CYGPATH_W) $(JAVA_HOME)/include/win32` -I.
#INCL = -I`$(CYGPATH_W) /home/viz/Downloads/Ipopt-3.11.7/include/coin`  -I`$(CYGPATH_W) $(JAVA_HOME)/include` -I`$(CYGPATH_W) $(JAVA_HOME)/include/win32` -I.

# Linker flags
LIBS = `PKG_CONFIG_PATH=/home/viz/Downloads/Ipopt-3.11.7/lib64/pkgconfig:/home/viz/Downloads/Ipopt-3.11.7/lib/pkgconfig:/home/viz/Downloads/Ipopt-3.11.7/share/pkgconfig: pkg-config --libs ipopt`
##LIBS = -link -libpath:`$(CYGPATH_W) /home/viz/Downloads/Ipopt-3.11.7/lib` libipopt.lib -lblas -lm  -ldl
#LIBS = -L/home/viz/Downloads/Ipopt-3.11.7/lib -lipopt -lblas -lm  -ldl

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

SRCDIR = .


##########################################################################

ORGDIR = org/coinor
JNIHEADER = org_coinor_Ipopt.h

JIPOPTOBJ = jipopt.lo
JIPOPTDLL = libjipopt.so
LIBJIPOPT = lib/$(JIPOPTDLL)

IPOPTCLASS = $(ORGDIR)/Ipopt.class
HS071 = $(ORGDIR)/examples/HS071
SCALABLE = $(ORGDIR)/examples/scalable/SolveProblem

HS071CLASSES = $(ORGDIR)/examples/HS071.class

SCALCLASSES = \
	$(ORGDIR)/examples/scalable/LuksanVlcek1.class \
	$(ORGDIR)/examples/scalable/Scalable.class \
	$(ORGDIR)/examples/scalable/SolveProblem.class

SCALJAR = ScalableProblems.jar

$(LIBJIPOPT): $(JIPOPTOBJ) lib
	$(LIBTOOL) --tag=CXX --mode=link $(CXX) $(LINKFLAGS) -o $(LIBJIPOPT) $(JIPOPTOBJ) $(LIBS)

$(JNIHEADER): $(IPOPTCLASS)
	$(JAVAH) `echo $(IPOPTCLASS) | sed -e 's|/|.|g' -e 's|.class||'`

$(JIPOPTOBJ): $(JNIHEADER) 

lib:
	mkdir lib

.java.class:
	$(JAVAC) $<

.SUFFIXES: .cpp .lo .obj .java .class

$(SCALJAR): $(SCALCLASSES)
	$(JAR) cfm $(SCALJAR) $(SRCDIR)/manifest.txt $(IPOPTCLASS) $(SCALCLASSES)

test: $(LIBJIPOPT) $(HS071CLASSES) $(SCALJAR)
	$(JAVA) $(HS071)
	$(JAVA) -jar $(SCALJAR) LukVlI1 100
	$(JAVA) -jar $(SCALJAR) LukVlE1 101

clean:
	rm -rf $(LIBJIPOPT) $(JIPOPTOBJ) $(JNIHEADER) $(IPOPTCLASS) $(HS071CLASSES) $(SCALCLASSES) $(SCALJAR)

doc:
	$(JAVADOC) -d javadoc -windowtitle "JIpopt API documentation" org/coinor/Ipopt.java

.cpp.lo:
	$(LIBTOOL) --tag=CXX --mode=compile $(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `test -f '$<' || echo '$(SRCDIR)/'`$<

.cpp.obj:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `if test -f '$<'; then $(CYGPATH_W) '$<'; else $(CYGPATH_W) '$(SRCDIR)/$<'; fi`
