# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME = ../../amrex-master

DEBUG     = TRUE
USE_MPI   = FALSE
USE_OMP   = FALSE
COMP      = gnu
DIM       = 2

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
VPATH_LOCATIONS  += ./
INCLUDE_LOCATIONS += ./ /usr/local/include ../adept-2.0.5/include

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules


