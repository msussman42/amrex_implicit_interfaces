#
# Generic setup for using Intel
#
CXX = icpc
CC  = icc
FC  = ifort
F90 = ifort

CXXFLAGS =
CFLAGS   =
FFLAGS   =
F90FLAGS =

########################################################################

intel_version = $(shell $(CXX) -dumpversion)

COMP_VERSION = $(intel_version)

########################################################################

ifeq ($(DEBUG),TRUE)

  CXXFLAGS += -g -O0 -std=c++14
  CFLAGS   += -g -O0 -std=c++14
  FFLAGS   += -g -O0 -check all
  F90FLAGS += -g -O0 -check all

else

  CXXFLAGS += -g -O2 -ip -qopt-report=5 -qopt-report-phase=vec -std=c++14
  CFLAGS   += -g -O2 -ip -qopt-report=5 -qopt-report-phase=vec -std=c++14
  FFLAGS   += -g -O2 -ip -qopt-report=5 -qopt-report-phase=vec
  F90FLAGS += -g -O2 -ip -qopt-report=5 -qopt-report-phase=vec

endif

########################################################################

#CFLAGS   += -std=c99

F90FLAGS += -implicitnone

FMODULES = -module $(fmoddir) -I$(fmoddir)

########################################################################

GENERIC_COMP_FLAGS =

ifeq ($(USE_OMP),TRUE)
  ifeq ($(firstword $(sort 16.0 $(intel_version))), 16.0) 
    GENERIC_COMP_FLAGS += -qopenmp
  else
    GENERIC_COMP_FLAGS += -openmp
  endif
endif

CXXFLAGS += $(GENERIC_COMP_FLAGS) -pthread
CFLAGS   += $(GENERIC_COMP_FLAGS)
FFLAGS   += $(GENERIC_COMP_FLAGS)
F90FLAGS += $(GENERIC_COMP_FLAGS)

########################################################################

override XTRALIBS += -lifcore -lmpi

ifeq ($(USE_OMP),TRUE)
  override XTRALIBS += -lifcoremt
endif

LINK_WITH_FORTRAN_COMPILER ?= $(USE_F_INTERFACES)

ifeq ($(LINK_WITH_FORTRAN_COMPILER),TRUE)
  override XTRALIBS += -lstdc++
endif