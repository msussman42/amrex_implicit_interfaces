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

ifndef COMP
  COMP = gcc
endif

ifndef DEBUG
  DEBUG = TRUE
endif

ifndef PROFILE
  PROFILE = FALSE
endif

ifeq ($(COMP),KCC)
  CXX = KCC
  CC  = KCC --c
endif
ifeq ($(COMP),g++)
  CXX = g++
  CC  = gcc
endif
ifeq ($(COMP),CC)
  CXX = CC
  CC  = cc
endif
ifeq ($(COMP),xlc)
  CXX = xlC
  CC  = xlc
endif
ifeq ($(COMP),Intel)
  CXX = icc
  CC  = icc
endif
ifeq ($(COMP),pgi)
  CXX = pgCC
  CC  = pgcc
endif

ifeq ($(PROFILE),TRUE)
  CXXFLAGS += -pg
endif

ifneq ($(DEBUG),TRUE)
  CPPFLAGS += -DNDEBUG
endif

ifeq ($(COMP),KCC)
  CXXFLAGS += --strict
  CXXFLAGS += --display_error_number
  CXXFLAGS += --diag_suppress 450
  CFLAGS += --strict
  CFLAGS   += --display_error_number
  # CFLAGS   += --diag_suppress 267
  ifneq ($(DEBUG),FALSE)
    CXXFLAGS += +K0 -g
    CFLAGS += +K0 -g
  else
    CXXFLAGS += +K3 -O3 -g
    CFLAGS += +K3 -O3 -g
  endif
else
  ifeq ($(COMP),g++)
    #  CXXFLAGS += -ansi
    #  CXXFLAGS += -Wall
  endif
  CXXFLAGS += -g -O
  CFLAGS   += -g -O
#  CFLAGS   += -ansi
#  CFLAGS   += -Wall
  ifeq ($(DEBUG),FALSE)
    CXXFLAGS += -O3
    CFLAGS   += -O3
  endif
  ifeq ($(COMP),Intel)
#    CXXFLAGS += -axK
#    CFLAGS   += -axK
  endif
endif

