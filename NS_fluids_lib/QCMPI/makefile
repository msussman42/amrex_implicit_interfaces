
grover:
	ftn -o GN.out  -Msave  -O3  Guniv.f90  qcmpisubs.f90  mpi_vars.f90

shor:
	ftn -o SH.out  -Msave  -O3  Suniv.f90   qcmpisubs.f90 mpi_vars.f90   

groverpople:
	ifort -o GN.out Guniv.f90  qcmpisubs.f90 mpi_vars.f90 -L${MKL_PATH} -lmkl_scalapack -lmkl_blacs_sgimpt_lp64 -lmkl -lguide -lmpi -lpthread

shorpople:
	ifort -o SH.out Suniv.f90 qcmpisubs.f90 mpi_vars.f90 -L${MKL_PATH} -lmkl_scalapack -lmkl_blacs_sgimpt_lp64 -lmkl -lguide -lmpi -lpthread

clean:
	rm -f *.o *.out 
cleanfort:
	rm -f fort.*
