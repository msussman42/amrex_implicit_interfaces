TINY_PROFILE = TRUE

PRECISION=FLOAT

DEBUG = FALSE
DIM=3
USE_MPI = FALSE
USE_OMP = FALSE

USE_CUDA = TRUE
COMP = gnu

# ask gfortran the name of the library to link in.  First check for the
# static version.  If it returns only the name w/o a path, then it
# was not found.  In that case, ask for the shared-object version.
gfortran_liba  = $(shell $(F90) -print-file-name=libgfortran.a)
gfortran_libso = $(shell $(F90) -print-file-name=libgfortran.so)

ifneq ($(gfortran_liba),libgfortran.a)  # if found the full path is printed, thus `neq`.
  LIBRARY_LOCATIONS += $(dir $(gfortran_liba))
else
  LIBRARY_LOCATIONS += $(dir $(gfortran_libso))
endif

override XTRALIBS += -lgfortran

quadmath_liba  = $(shell $(F90) -print-file-name=libquadmath.a)
quadmath_libso = $(shell $(F90) -print-file-name=libquadmath.so)

ifneq ($(quadmath_liba),libquadmath.a)
  override XTRALIBS += -lquadmath
else ifneq ($(quadmath_libso),libquadmath.so)
  override XTRALIBS += -lquadmath
endif

include ./Make.package
include ./Make.CNS
