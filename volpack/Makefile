# (c) 1996-2000 The Regents of
# the University of California (through E.O. Lawrence Berkeley National
# Laboratory), subject to approval by the U.S. Department of Energy.

# Your use of this software is under license -- the license agreement is
# attached and included in the directory as license.txt or you may
# contact Berkeley Lab's Technology Transfer Department at TTD@lbl.gov.

# NOTICE OF U.S. GOVERNMENT RIGHTS.  The Software was developed under
# funding from the U.S. Government which consequently retains certain
# rights as follows: the U.S. Government has been granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable,
# worldwide license in the Software to reproduce, prepare derivative
# works, and perform publicly and display publicly.  Beginning five (5)
# years after the date permission to assert copyright is obtained from
# the U.S. Department of Energy, and subject to any subsequent five (5)
# year renewals, the U.S. Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide
# license in the Software to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
  
COMP  := KCC
COMP  := xlc
COMP  := g++
#COMP  := pgi
#COMP  := Intel
DEBUG := TRUE
DEBUG := FALSE
#PROFILE := TRUE

include Include.mk

M4 	:= m4
M4FLAGS := -s
M4FLAGS :=
RANLIB 	:= ranlib

# loop unrolling optimization: if the following line is not commented
# out then the library will be built with a special version of the
# rendering code that has been hand-optimized; use this unless your
# compiler's optimizer has trouble with procedures that have a large
# number of basic blocks
#CPPFLAGS += -DUNROLL_RUN_LOOP
CPPFLAGS += -DUSE_TIMER
ifeq ($(DEBUG),TRUE)
  CPPFLAGS += -DDEBUG
endif
ifeq ($(MPI),ON)
  CPPFLAGS += -DVP_MPI
endif

# list of optional C source files to build with code specialized
# for particular shaders and resampling filters
optional_comp_sources = vp_compAC11B.c vp_compAC31B.c vp_compAC32B.c \
			vp_compAR11B.c vp_compAR31B.c vp_compAR32B.c

optional_comp_objects =	$(optional_comp_sources:.c=.o)

# C define flags corresponding to COMP_SRCS; these tell the dispatch
# routines in the library which object files are available
CPPFLAGS +=	-DCOMP_AC11B   -DCOMP_AC31B   -DCOMP_AC32B \
		-DCOMP_AR11B   -DCOMP_AR31B   -DCOMP_AR32B

CPPFLAGS += -I.

CPPFLAGS += -DBL_BYTEORDER_LITTLE_ENDIAN

c_sources = vp_check.c vp_context.c vp_extract.c vp_file.c vp_linalg.c vp_octree.c \
              vp_renderA.c vp_renderB.c vp_renderC.c vp_renderR.c vp_resample.c \
              vp_rle.c vp_shade.c vp_transpose.c vp_util.c vp_view.c vp_warp.c
c_headers = volpack.h vp_global.h
m4_sources = vp_CompA.m4 vp_WarpA.m4

comp_sources = vp_compAC1NB.c vp_compAC3NB.c vp_compAR1NB.c vp_compAR3NB.c \
       vp_compAC1PB.c vp_compAC3PB.c vp_compAR1PB.c vp_compAR3PB.c \
       vp_compAC00G.c vp_compAR00G.c vp_compAC1NS.c vp_compAC3NS.c \
       vp_compAR1NS.c vp_compAR3NS.c $(optional_comp_sources)
comp_objects = $(comp_sources:.c=.o) $(optional_comp_objects)

warp_sources = vp_warpA101N.c vp_warpA301N.c vp_warpA110N.c vp_warpA111N.c \
       vp_warpA330N.c vp_warpA331N.c vp_warpA330R.c vp_warpA331R.c
warp_objects  = $(warp_sources:.c=.o)

generated = vp_opts.c $(comp_sources) $(warp_sources)

objects = $(c_sources:.c=.o) vp_opts.o $(comp_objects) $(warp_objects)
comp_input = no_edit_header vp_CompA.m4
warp_input = no_edit_header vp_WarpA.m4

sources = $(c_sources) $(c_headers) $(m4_sources)

all: libvolpack.a

TAGS: $(sources)
	etags $(sources)

libvolpack.a: $(objects)
	$(RM) libvolpack.a
	$(AR) $(ARFLAGS) libvolpack.a $(objects)
	$(RANLIB) libvolpack.a

$(objects): volpack.h vp_global.h

makeopts: makeopts.c
	$(CC) -o makeopts makeopts.c

vp_opts.c: makeopts Makefile Include.mk
	./makeopts vp_opts.c $(CXXFLAGS) $(CPPFLAGS)

vp_compA%.c: $(comp_input)
	@$(RM) $@
	@echo "$(M4) $(M4FLAGS) -DSourceFile=$@ $(comp_input) > $@"
	@$(M4) $(M4FLAGS) -DSourceFile=$@ $(comp_input) > $@
#	@chmod -w $@

vp_warpA%.c: $(warp_input)
	@$(RM) $@
	@echo "$(M4) $(M4FLAGS) -DSourceFile=$@ $(warp_input) > $@"
	@$(M4) $(M4FLAGS) -DSourceFile=$@ $(warp_input) > $@
#	@chmod -w $@

examples: libvolpack.a
	make -C examples COMP=$(COMP) MPI=$(MPI) DEBUG=$(DEBUG) PROFILE=$(PROFILE)

clean:
	$(RM) $(objects) $(generated) makeopts
	rm -f libvolpack.a
	$(RM) TAGS
	$(RM) *~
	make -C examples clean

realclean: clean

.PHONY: clean examples

.PRECIOUS: $(generated)

ifndef COMPILE.cpp
  LINK.cpp = $(LINK.cc)
  COMPILE.cpp=$(COMPILE.cc)
  %.cpp:
  %: %.cpp
	$(LINK.cpp) $^ $(LOADLIBES) $(LDLIBS) -o $@
  %.o: %.cpp
	$(COMPILE.cpp) $< $(OUTPUT_OPTION)
endif
