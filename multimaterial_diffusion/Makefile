test_opt:
	gfortran  -g -O3 -ffree-line-length-none -fno-range-check -fno-second-underscore -fimplicit-none  -DINCLUDE_FIB -DBL_LANG_FORT -DNDEBUG -DBL_GCC_VERSION='4.8.5' -DBL_GCC_MAJOR_VERSION=4 -DBL_GCC_MINOR_VERSION=8 -DBL_SPACEDIM=2 -DAMREX_SPACEDIM=2 -I. -DBL_FORT_USE_UNDERSCORE -DBL_Linux  multimat_FVM.F90 -c -o multimat_FVM.o
	gfortran  -g -O3 -ffree-line-length-none -fno-range-check -fno-second-underscore -fimplicit-none  -DINCLUDE_FIB -DBL_LANG_FORT -DNDEBUG -DBL_GCC_VERSION='4.8.5' -DBL_GCC_MAJOR_VERSION=4 -DBL_GCC_MINOR_VERSION=8 -DBL_SPACEDIM=2 -DAMREX_SPACEDIM=2 -I. -DBL_FORT_USE_UNDERSCORE -DBL_Linux  main.F90 -c -o main.o 
	gfortran  -o test_opt heightfunctioncurvature.o tetrahedron_keast_rule.o triangle_fekete_rule.o PROBCOMMON.o PROBMAIN.o GLOBALUTIL.o TECPLOTUTIL.o supercooled_exact_sol.o variable_temperature_drop.o mof2d_analytic_centroid.o mof3d_analytic_centroid.o MOF.o STACKVOLUME.o BICGSTAB_Yang_MULTI.o vof_cisl.o TSAT_UTIL.o vfrac_pair.o multimat_FVM.o MOF_REDIST_3D.o PLIC_3D.o sci_clsvof.o SOLIDFLUID.o MASS_TRANSFER_3D.o MOVE_INTERFACE_3D.o integrate.o main.o

run:
	./test_opt
clean:
	rm *.mod core.* test_opt *.jpg *.o *.tec nodedata* cendata*
clear:	
	rm *.dat
