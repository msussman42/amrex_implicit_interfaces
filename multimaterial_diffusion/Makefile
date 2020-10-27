test_debug:
	gfortran -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv -fimplicit-none -fno-second-underscore -DBL_SPACEDIM=2 -DAMREX_SPACEDIM=2 -I. -DBL_FORT_USE_UNDERSCORE  -g multimat_FVM.F90 -c -o multimat_FVM.o
	gfortran -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv -fimplicit-none -fno-second-underscore -DBL_SPACEDIM=2 -DAMREX_SPACEDIM=2 -I. -DBL_FORT_USE_UNDERSCORE  -g BICGSTAB_Yang_MULTI.F90 -c -o BICGSTAB_Yang_MULTI.o 
	gfortran -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv -fimplicit-none -fno-second-underscore -DBL_SPACEDIM=2 -DAMREX_SPACEDIM=2 -I. -DBL_FORT_USE_UNDERSCORE  -g MOF_REDIST_3D.F90 -c -o MOF_REDIST_3D.o 
	gfortran -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv -fimplicit-none -fno-second-underscore -DBL_SPACEDIM=2 -DAMREX_SPACEDIM=2 -I. -DBL_FORT_USE_UNDERSCORE  -g PLIC_3D.F90 -c -o PLIC_3D.o 
	gfortran -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv -fimplicit-none -fno-second-underscore -DBL_SPACEDIM=2 -DAMREX_SPACEDIM=2 -I. -DBL_FORT_USE_UNDERSCORE  -g sci_clsvof.F90 -c -o sci_clsvof.o 
	gfortran -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv -fimplicit-none -fno-second-underscore -DBL_SPACEDIM=2 -DAMREX_SPACEDIM=2 -I. -DBL_FORT_USE_UNDERSCORE  -g SOLIDFLUID.F90 -c -o SOLIDFLUID.o 
	gfortran -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv -fimplicit-none -fno-second-underscore -DBL_SPACEDIM=2 -DAMREX_SPACEDIM=2 -I. -DBL_FORT_USE_UNDERSCORE  -g MOVE_INTERFACE_3D.F90 -c -o MOVE_INTERFACE_3D.o 
	gfortran -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv -fimplicit-none -fno-second-underscore -DBL_SPACEDIM=2 -DAMREX_SPACEDIM=2 -I. -DBL_FORT_USE_UNDERSCORE  -g integrate.F90 -c -o integrate.o 
	gfortran -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv -fimplicit-none -fno-second-underscore -DBL_SPACEDIM=2 -DAMREX_SPACEDIM=2 -I. -DBL_FORT_USE_UNDERSCORE  -g main.F90 -c -o main.o 
	gfortran  -o test_debug heightfunctioncurvature.o tetrahedron_keast_rule.o triangle_fekete_rule.o PROBCOMMON.o PROBMAIN.o GLOBALUTIL.o TECPLOTUTIL.o supercooled_exact_sol.o variable_temperature_drop.o mof2d_analytic_centroid.o mof3d_analytic_centroid.o MOF.o STACKVOLUME.o BICGSTAB_Yang_MULTI.o vof_cisl.o TSAT_UTIL.o vfrac_pair.o multimat_FVM.o MOF_REDIST_3D.o PLIC_3D.o sci_clsvof.o SOLIDFLUID.o MASS_TRANSFER_3D.o MOVE_INTERFACE_3D.o integrate.o main.o

run:
	./test_debug
clean:
	rm *.mod core.* test_debug *.jpg *.o *.tec nodedata* cendata*
clear:	
	rm *.dat
