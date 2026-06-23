
grover:
	mpifort -o GN.out -O3 mpi_vars.f90 qcmpisubs.f90 Guniv.f90

shor:
	mpifort -o SH.out -O3 mpi_vars.f90 qcmpisubs.f90 Suniv.f90

groverpople:
	mpifort -o GN.out Guniv.f90  qcmpisubs.f90 mpi_vars.f90 -L${MKL_PATH} -lmkl_scalapack -lmkl_blacs_sgimpt_lp64 -lmkl -lguide -lmpi -lpthread

shorpople:
	mpifort -o SH.out Suniv.f90 qcmpisubs.f90 mpi_vars.f90 -L${MKL_PATH} -lmkl_scalapack -lmkl_blacs_sgimpt_lp64 -lmkl -lguide -lmpi -lpthread

clean:
	rm -f *.o *.out 
cleanfort:
	rm -f fort.*
