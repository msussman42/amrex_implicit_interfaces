# in bash:
# git add -u  --> deletes all removed files, updates modifications, and adds 
#     new files.
# git rm $(git ls-files --deleted)
# git clone https://github.com/AMReX-Codes/amrex.git
# For post processing:
# git clone https://github.com/AMReX-Codes/amrex.git
# git clone https://github.com/AMReX-Combustion/PeleAnalysis.git 
AMREX_HOME = ../amrex-master
MOF_HOME = ..

# in Make.machines:
# which_site := unknown
# which_computer := unknown
#
# ifdef HOSTNAME
#  host_name := $(strip $(HOSTNAME))
# else ifdef HOST
#  host_name := $(strip $(HOST))
# else
#  host_name := $(shell hostname -f)
# endif
#

PRECISION = DOUBLE
DEBUG	  = FALSE
USE_ASSERTION = TRUE
PROFILE = FALSE
MEM_PROFILE = FALSE
DIM       = 2

BL_NO_FORT = TRUE

MEHDI_VAHAB_FSI = TRUE
CTML_FIB = TRUE
CTML_ESH = FALSE
CTML_FSH = FALSE
CTML_FBC = FALSE

# LAMMPS:
# November 9, 8:00am, 2021
# lammps-29Sep2021
# cd src
# cd MAKE
# cd ..
# make yes-molecule
#  "Installing package molecule"
# make yes-extra-dump
#  "Installing package extra-dump"
# make serial
# files created: liblammps.a, liblammps_serial.a, lmp_serial
# mkdir ~sussman/lammps_run
# cp lmp_serial ~sussman/lammps_run/.
# put rod.LAMMPS_script and rod.LAMMPS_read_data in ~sussman/lammps_run
# cd ~sussman/lammps_run
# ./lmp_serial -in rod.LAMMPS_script
# VMD to visualize "rod_LAMMPS_atom.dcd"
#
# heaviside: scl enable devtoolset-2 python27 bash
# bromwich or compute1: scl enable devtoolset-8 bash
# Bromwich PGI:
#  http://www.pgroup.com/resources/docs/17.10/x86/pgi-install-guide/index.htm
#  bash
#  PGI=/opt/pgi; export PGI 
#  MANPATH=$MANPATH:$PGI/linux86-64/17.10/man; export MANPATH 
#  LM_LICENSE_FILE=$PGI/license.dat; export LM_LICENSE_FILE 
#  PATH=$PGI/linux86-64/17.10/bin:$PATH; export PATH
#  open mpi:
#  export PATH=$PGI/linux86-64/17.10/mpi/openmpi/bin:$PATH
#  export MANPATH=$MANPATH:$PGI/linux86-64/17.10/mpi/openmpi/man
#  hpc-rcc using Mehdi's script: CMOF/run_example/gaslens/create_sims.sh
#
#  Parallel debugger: ARM DDT or ARM Forge DDT (Allinea), 
#    Totalview, Intel Parallel Studio.
# -fno-range-check
#    Disable range checking on results of simplification of constant
#    expressions during compilation.

MAC_LAPTOP = FALSE
LINUX_DESKTOP = TRUE
HPC_GNU = FALSE

ifeq ($(HPC_GNU),TRUE)

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

FFLAGS   += $(GENERIC_GNU_FLAGS)
F90FLAGS += $(GENERIC_GNU_FLAGS)
endif
ifeq ($(LINUX_DESKTOP),TRUE)

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

FFLAGS   += $(GENERIC_GNU_FLAGS)
F90FLAGS += $(GENERIC_GNU_FLAGS)

endif

ifeq ($(MAC_LAPTOP),TRUE)
 COMP = g++-8
 FCOMP = gfortran-8
 AMREX_CCOMP=g++-8
 AMREX_FCOMP=gfortran-8
 CXX = g++-8
 CC = g++-8
 FC = gfortran-8
 F90 = gfortran-8
endif

# total view:
#  1. ssh -YX sussman@hpc-login.rcc.fsu.edu (same password as FSUID)
#      https://rcc.fsu.edu/manage
#  2. module load gnu openmpi  (version 8.3.1)
#  2.5 cd multimaterialTV
#  2.6 get latest code, USE_MPI=TRUE
#  3. make clean
#  4. srun  --pty -t 30:00 -n 4  -p genacc_q /bin/bash
#  5. module load tv
#  6. module load gnu openmpi (version 8.3.1)
#  7. mpirun -tv -np 4 ./amr2<tab> inputs<tab>
#     (same directory where program compiled)
#  8. A few windows appear, press "GO", then program output will appear
#     in the xterminal.  "HALT" to pause, "GO" to continue.
#  For memory debugging: "DEBUG" "OPEN MemoryScape"
#  Memory Debugging Options=Extreme
# module load gnu openmpi (version 8.3.1)
# module load intel openmpi (for main program only; libraries use gnu openmpi)
# mecfd_q  (dedicated nodes)
# mecfd18_q (dedicated nodes)
# backfill
# engineering_q
# engineering_long
# hpc-login.rcc.fsu.edu  "module load intel openmpi" (main program)  
#                        "module load gnu openmpi (version 8.3.1)"
# sussman + same password as FSUID
#   https://rcc.fsu.edu/manage
# bromwich or compute1 GNU MPI:
#  "module load mpi/openmpi-4-gcc8"
# csh: command >& out_and_err_file &
# spear-login.rcc.fsu.edu
# ifort -g -O0 -traceback -check bounds,uninit,pointers
# mpirun -np 2 valgrind -v --leak-check=full --show-reachable=yes 
#   --log-file=run.out.%p ../../amr3d.gnu.DEBUG.MPI.ex inputs3d.skirted  
#   1> run.out 2> run.err &
# ``Memory error checking in C and C++: Comparing Sanitizers and Valgrind.''
#
#FCOMP = intel
#COMP = intel
#COMP = PGI 
#FCOMP = PGI
USE_MPI = TRUE
# set OMP_NUM_THREADS=<nthreads>
# export OMP_NUM_THREADS=<nthreads>
USE_OMP = FALSE
#EBASE = amr
EBASE = amr

include $(AMREX_HOME)/Tools/GNUMake/Make.defs 

$(info MODDEP=$(MODDEP))

ifeq ($(MEHDI_VAHAB_FSI),TRUE)
 CPPFLAGS+=-DMVAHABFSI
 CXXFLAGS+=-DMVAHABFSI
 CFLAGS+=-DMVAHABFSI
 FFLAGS+=-DMVAHABFSI
 F90FLAGS+=-DMVAHABFSI
endif

ifeq ($(CTML_FIB),TRUE)
 CPPFLAGS+=-DINCLUDE_FIB
 FFLAGS+=-DINCLUDE_FIB
 F90FLAGS+=-DINCLUDE_FIB
endif
ifeq ($(CTML_ESH),TRUE)
 CPPFLAGS+=-DINCLUDE_ESH
 FFLAGS+=-DINCLUDE_ESH
 F90FLAGS+=-DINCLUDE_ESH
endif
ifeq ($(CTML_FSH),TRUE)
 CPPFLAGS+=-DINCLUDE_FSH
 FFLAGS+=-DINCLUDE_FSH
 F90FLAGS+=-DINCLUDE_FSH
endif
ifeq ($(CTML_FBC),TRUE)
 CPPFLAGS+=-DINCLUDE_FBC
 FFLAGS+=-DINCLUDE_FBC
 F90FLAGS+=-DINCLUDE_FBC
endif

ifeq ($(MEHDI_VAHAB_FSI),TRUE)
include $(MOF_HOME)/StructureCodeShoele/Make.package
endif

include $(AMREX_HOME)/Src/Base/Make.package

include ./Make.package
include $(MOF_HOME)/amrlib/Make.package

include $(MOF_HOME)/mglib/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

