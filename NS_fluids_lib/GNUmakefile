# in bash:
# git add -u  --> deletes all removed files, updates modifications, and adds 
#     new files.
# git rm $(git ls-files --deleted)
# git clone https://github.com/AMReX-Codes/amrex.git
AMREX_HOME = ../amrex-master
MOF_HOME = ..

PRECISION = DOUBLE
DEBUG	  = TRUE
USE_ASSERTION = TRUE
PROFILE = FALSE
MEM_PROFILE = FALSE
DIM       = 3
CYGWIN    = FALSE
CYGWIN_ZENBOOK = FALSE
MEHDI_VAHAB_FSI = FALSE
CTML_FIB = TRUE
CTML_ESH = FALSE
CTML_FSH = FALSE
CTML_FBC = FALSE

# CYGWIN:
# CYGWIN=TRUE
#  (causes: -Wl,--image-base -Wl,0x10000000)
#  CYGWIN packages needed:
#  1. default installation
#  2. cmake
#  3. editors
#  4. gcc (search)
#  5. mingw (search)
#  6. X11 
#  7. libc++
#  8. cpp (search)
#  9. gfortran
# 10. g++
# 11. make
# 12. vim
# 13. emacs
# 14. procps (search)  (corresponds to "top")
# heaviside: scl enable devtoolset-2 python27 bash
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
COMP	  = g++
FCOMP = gfortran
# total view:
#  1. ssh -YX sussman@hpc-login.rcc.fsu.edu
#  2. module load gnu-openmpi  (can I leave out the /2.1.0 ?)
#  2.5 cd multimaterialTV
#  2.6 get latest code, USE_MPI=TRUE
#  3. make clean
#  4. srun  --pty -t 30:00 -n 4  -p genacc_q /bin/bash
#  5. module load tv
#  6. module load gnu-openmpi
#  7. mpirun -tv -np 4 ./amr2<tab> inputs<tab>
#     (same directory where program compiled)
#  8. A few windows appear, press "GO", then program output will appear
#     in the xterminal.  "HALT" to pause, "GO" to continue.
#  For memory debugging: "DEBUG" "OPEN MemoryScape"
#  Memory Debugging Options=Extreme
# module load gnu-openmpi/2.1.0
# module load intel-openmpi/2.1.0
# mecfd_q  (dedicated nodes)
# mecfd18_q (dedicated nodes)
# backfill
# engineering_q
# engineering_long
# hpc-login.rcc.fsu.edu  "module load intel"
# sussman + old password
# bromwich or compute1 GNU MPI:
#  "module load mpi/openmpi-x86_64"
# csh: command >& out_and_err_file &
# if "ImportError: no module named site" then give command:
#  "export PYTHONHOME="
# spear-login.rcc.fsu.edu
# ifort -g -O0 -traceback -check bounds,uninit,pointers
#valgrind -v --leak-check=full --show-reachable=yes ./a.out 1> run.out 2> run.err  -Mpreprocess -> pgf95 
# bromwich intel:
# source /usr/local/intel/parallel_studio_xe_2018/psxevars.csh (csh) or
# source /usr/local/intel/parallel_studio_xe_2018/psxevars.sh (bash)
#FCOMP = intel
#COMP = intel
#COMP = PGI 
#FCOMP = PGI
USE_MPI   = TRUE
# set OMP_NUM_THREADS=<nthreads>
# export OMP_NUM_THREADS=<nthreads>
USE_OMP = FALSE
EBASE = amr
#EBASE = amrcompute

AMREX_PARTICLES = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs 

ifeq ($(CYGWIN),TRUE)
 CXXFLAGS+=-Wl,--image-base -Wl,0x10000000
 CFLAGS+=-Wl,--image-base -Wl,0x10000000
 FFLAGS+=-Wl,--image-base -Wl,0x10000000
 F90FLAGS+=-Wl,--image-base -Wl,0x10000000
endif
ifeq ($(CYGWIN_ZENBOOK),TRUE)
 CXXFLAGS+=-DCYGWINZENBOOK
 CFLAGS+=-DCYGWINZENBOOK
 FFLAGS+=-DCYGWINZENBOOK
 F90FLAGS+=-DCYGWINZENBOOK
endif

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
include $(MOF_HOME)/Vicar3D/Make.package
endif

include ./Make.package

include $(MOF_HOME)/mglib/Make.package
include $(MOF_HOME)/amrlib/Make.package

include $(AMREX_HOME)/Src/Base/Make.package
#include $(AMREX_HOME)/Src/C_ParticleLib/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package
include $(AMREX_HOME)/Tools/GNUMake/Make.rules

