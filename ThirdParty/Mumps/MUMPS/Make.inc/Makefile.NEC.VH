#
# ------------------------------------------------------------------------------------
# This Makefile part is dedicated to the use of VH offload of scalar routines of MUMPS
# ------------------------------------------------------------------------------------
#

#
# Some symbols of METIS are still used on the VE. METIS on VE is thus needed.
#

LMETISDIR     = -L/path/to/ve/metis/lib
IMETIS        = -I/path/to/ve/metis/include

#
# We need the metis on the Vector Host which will run when offloading is triggered.
# Make sure VH metis is compiled with -fPIC
#

VHMETISLIBDIR = /path/to/vh/metis/lib

#
# The next section is related to the compilation of the interface of the code to be
# offloaded (VE to VH) and to the generation of the Vector Host library (libvh.so).
# The GNU compiler on your system will be used to generate the VH library.
# Note that in order to work properly, you must add the VHOUTPUTDIR path to
# the LD_LIBRARY_PATH.
#

VHOUTPUTDIR = $(topdir)/lib
VEINCDIR    = $(topdir)/src/ve/include
VESRCDIR    = $(topdir)/src/ve/src

VHMUMPSSRCDIR = $(topdir)/src
VHMUMPSOBJDIR = $(topdir)/src
VHMUMPSOBJS   = $(VHMUMPSOBJDIR)/ana_orderings.vh.o $(VHMUMPSOBJDIR)/ana_AMDMF.vh.o

VHCC    := gcc
VHFC    := gfortran
VHCCOPT := -O2 -fPIC
VHFCOPT := -O2 -fPIC
VEOBJS  := $(VESRCDIR)/VE_Metis_interface.o         \
           $(VESRCDIR)/VE_Metis_nodend.o            \
           $(VESRCDIR)/VE_Metis_setdefaultoptions.o \
           $(VESRCDIR)/VE_Ana_orderings_interface.o \
           $(VESRCDIR)/VE_Mumps_ana_h.o             \
           $(VESRCDIR)/VE_Mumps_symqamd.o           \
           $(VESRCDIR)/VE_Mumps_wrap_ginp94.o

$(VHMUMPSOBJDIR)/%.vh.o : $(VHMUMPSSRCDIR)/%.F
	$(VHFC) $(VHFCOPT) -c $< -o $@

$(VESRCDIR)/%.o : $(VESRCDIR)/%.c
	$(CC) $(OPTC) $(IMETIS) -I$(VEINCDIR) $(CDEFS) -c $< -o $@

$(VESRCDIR)/%.o : $(VESRCDIR)/%.f90
	$(FC) $(OPTF) -c $< -o $@

d c s z : vhlib
clean   : cleanvhlib

vhlib: $(VEOBJS) $(VHMUMPSOBJS)
	@mkdir -p $(VHOUTPUTDIR)
	@$(VHFC) -shared -Wl,-soname,libvh.so -o $(VHOUTPUTDIR)/libvh.so $(VHMUMPSOBJS) -Wl,--whole-archive $(VHMETISLIBDIR)/libmetis.a -Wl,--no-whole-archive

OBJS_COMMON_OTHER += $(VEOBJS)

cleanvhlib:
	@$(RM) $(VEOBJS) $(VHMUMPSOBJS)

ORDERINGSF  += -Dmetis
OPTF        += -DVHOFFLOAD -I$(VEINCDIR)
LMETIS       = $(LMETISDIR) -lmetis
