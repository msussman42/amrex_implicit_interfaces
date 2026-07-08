grover:
	gfortran -o GN.out -cpp -O3 mpi_vars.f90 qcmpisubs.f90 Guniv.f90

shor:
	gfortran -o SH.out -cpp -O3 mpi_vars.f90 qcmpisubs.f90 Suniv.f90

truth:
	gfortran -o TR.out -cpp -O3 mpi_vars.f90 qcmpisubs.f90 truth_table.f90

clean:
	rm -f *.o *.out 
cleanfort:
	rm -f fort.*
