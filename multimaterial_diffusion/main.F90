#undef BL_LANG_CC
#define BL_LANG_FORT

#include "AMReX_BC_TYPES.H"

PROGRAM test 
USE GeneralClass ! vof_cisl.F90
USE probmain_module 
USE probcommon_module 
USE LegendreNodes
USE global_utility_module 
USE MOF_pair_module ! vfrac_pair.F90
USE mmat_FVM  ! multimat_FVM.F90
USE bicgstab_module
USE supercooled_exact_sol
USE variable_temperature_drop

IMPLICIT NONE

! Dai and Scannapieco: no graph or table 
!  discussing convergence characteristics
!  of the gradient for multimaterial problems.
! Garimella and Lipnikov: 2nd order, but rate of convergence not 
!  investigated for filament problem?
! Dawes: surrogate supermesh - filament problem?
! Kinkinzon: filament over 1 cell thick
! Zhiliakov et al: ???
! Yang Liu
! problem type
! 13 = star with thin filament
!      (vof_cisl.F90: pentaeps, dirichlet_pentafoil)
! 14 = star for two material sanity check
! 15 = hypocycloid with 2 materials
! 16 = nucleate boiling diffusion with thin 
!      filament between vapor bubble and substrate
!      (vof_cisl.F90: thermal_delta, declared in vof_cisl but defined in main)
! 17 = hypocycloid with 5 materials
! 19 = annulus cvg test
! 20 = hypocycloid with 6 materials
!

! July 5, 2019 Infinity norm for LSexact( x_face_centroid ), probtype=4
!  512  9.3E-7  (areaface: 1.1E-7, dir,side=1,2 i,j=346,439)
!  256  3.83E-6  (areaface: 1.0E-8, dir,side=1,2 i,j=114,26)
!  128  1.46E-5  (areaface: 1.0E-8, dir,side=1,2 i,j=114,26)
!  64   5.96E-5  (areaface: 2.9E-6, dir,side=1,2 i,j=21,8)
!
! July 10, 2019 Infinity norm for LSexact( x_face_centroid ), probtype=16
! 256   3.0E-5 (areaface: 3.0E-6)
! 128   1.1E-4 (areaface: 1.4E-5)
!  64   4.4E-4 (areaface: 7.5E-5)
! YANG LIU MODIFY TABLE
! probtype==13 (pentafoil), dirichlet bc, ERRTOL=0.999999D0, thick (0.2d0),
! TSTOP=1.25E-2, operator_ext=1 operator_int=3 linear exact=1,
!       L1      L2      LINF    L1_grd  L2_grd   LINF_grd Linf_Igrd L1_Igrd
!32(1)  8.7E-5  9.7E-5  3.5E-4   6.5E-4  1.1E-3   5.6E-3    0.48    8.0E-2
!64(4)  2.3E-5  2.4E-5  4.4E-5   1.5E-4  2.1E-4   8.4E-4    0.57    6.7E-2
!128(16) 5.9E-6 6.0E-6  2.2E-5   3.8E-5  6.3E-5   1.1E-3    0.54    4.0E-2
!256(64) 1.4E-6 1.5E-6  4.5E-6   1.0E-5  1.9E-5   4.3E-4    0.62    4.0E-2
!
! YANG LIU MODIFY TABLE
! probtype==13 (pentafoil), dirichlet bc, ERRTOL=0.01D0, thin (0.05d0),
! TSTOP=1.25E-2, operator_ext=1 operator_int=3 linear exact=1,
! note: no stencil available for getting interface flux.
!       L1      L2      LINF    L1_grd  L2_grd   LINF_grd 
! 32(1) 3.6E-4 5.9E-4  2.4E-3   2.4E-2 3.6E-2    1.5E-1   
! 64(2) 1.6E-4 2.2E-4  6.0E-4   1.5E-2 2.3E-2    8.7E-2   
!128(4) 1.5E-5 2.2E-5  1.1E-4   2.3E-3 4.5E-3    3.1E-2
!256(8) 1.6E-6 1.7E-6  9.8E-6   1.7E-4 3.0E-4    3.5E-3  
!
! YANG LIU MODIFY TABLE
! probtype==13 (pentafoil), neumann bc, ERRTOL=0.999999D0, thick (0.2d0),
! TSTOP=1.25E-2, operator_ext=1 operator_int=3 linear exact=1
! timesteps refined by factor of 2
! (OTHERWISE SIMULATIONS WILL TAKE A LONG TIME)
! 32(1), 64(2), 128(4), 256(8), 512(16)
!        L1         L2   Linf   L1 grad    L2 grad Linf grad  
!32/64   2.7E-3  3.0E-3  5.2E-3 3.7E-2     3.7E-2  3.9E-2              
!64/128  1.7E-3  1.9E-3  3.0E-3 1.8E-2     1.8E-2  2.6E-2              
!128/256 9.5E-4  1.1E-3  1.7E-3 9.3E-3     1.0E-2  1.6E-2              
!256/512 5.1E-4  5.6E-4  8.8E-4 4.9E-3     5.4E-3  9.2E-3              
!
! YANG LIU MODIFY TABLE
! probtype==16 (multiscale), error for material 2, kratio=1000:1, 
! material 2 is smooth.  filament thickness=0.02
! TSTOP=2.0,  
! timesteps refined by factor of 2.  16(4) 32(8), 64(16), 128(32), 256(64)
! SIMPLE (smooth):
!            dx    L1         L2      Linf   L1 grad    L2 grad Linf grad  
!16/32  0.03125    1.2E-1   1.2E-1   1.2E-1  1.3E-2     1.5E-2  2.5E-2
!32/64  0.015625   3.1E-1   3.1E-1   3.1E-1  2.2E-2     2.2E-2  2.7E-2
!64/128 0.0078125  1.9E-1   1.9E-1   1.9E-1  1.3E-2     1.3E-2  1.8E-2
!128/256 3.90625E-3 8.0E-2  8.0E-2   8.1E-2  5.2E-3     5.3E-3  9.3E-3
!average flux:(16) -0.131 (32) -0.155 (64) -0.177 (128) -0.192 (256) -0.197
!
! YANG LIU MODIFY TABLE
! SIMPLE (material 2 has corners):
!            dx    L1         L2      Linf   L1 grad    L2 grad Linf grad  
!16/32  0.03125      STENCIL NOT AVAILABLE FOR GRADIENT ERROR
!32/64  0.015625   3.7E-1   3.7E-1   3.8E-1  1.3E-2     1.3E-2  1.3E-2
!64/128 0.0078125  2.4E-1   2.4E-1   2.4E-1  6.4E-3     6.8E-3  1.3E-2
!128/256 3.90625E-3 1.3E-1  1.3E-1   1.3E-1  3.7E-3     4.0E-3  9.2E-3
!average flux:(16) -0.121 (32) -0.156 (64) -0.177 (128) -0.187 (256) -0.191
!
! YANG LIU MODIFY TABLE
! DS:
!            dx    L1         L2      Linf   L1 grad    L2 grad Linf grad  
!16/32  0.03125    2.1E-1   2.1E-1   2.2E-1  2.1E-2     2.3E-2  3.4E-2
!32/64  0.015625   3.3E-1   3.3E-1   3.3E-1  2.5E-2     2.5E-2  3.0E-2
!64/128 0.0078125  2.0E-1   2.0E-1   2.0E-1  1.4E-2     1.4E-2  2.0E-2
!128/256 3.90625E-3 1.2E-1  1.2E-1   1.2E-1  7.7E-3     7.8E-3  2.1E-2
!average flux:(16) -0.131 (32) -0.148 (64) -0.172 (128) -0.186 (256) -0.195
!
! YANG LIU MODIFY TABLE
! OP:
!            dx    L1         L2      Linf   L1 grad    L2 grad Linf grad  
!16/32  0.03125    4.6E-3   4.7E-3   5.5E-3  4.3E-3     4.7E-3  7.8E-3
!32/64  0.015625   1.3E-1   1.3E-1   1.3E-1  9.6E-3     9.6E-3  1.3E-2
!64/128 0.0078125  1.2E-1   1.2E-1   1.2E-1  7.9E-3     8.1E-3  1.3E-2
!128/256 3.90625E-3 6.2E-2  6.2E-2   6.2E-2  3.9E-3     4.0E-3  7.1E-3
!average flux:(16) -0.137 (32) -0.173 (64) -0.185 (128) -0.194 (256)-0.199
!
!
! YANG LIU MODIFY TABLE
! FUTURE WORK: LINEAR EXACT METHOD WHEN K=constant
! AMR free test: does the error decrease without refining the Eulerian grid?
!  i.e. does adding extra materials cause the error to go down?
! kratio=1
! SIMPLE (smooth):
!            dx    L1         L2      Linf   L1 grad    L2 grad Linf grad  
!16     0.0625     4.7E-2  5.6E-2   1.0E-1  2.5E-1     3.3E-1  5.5E-1
!32     0.03125    3.2E-2  4.0E-2   8.6E-2  1.7E-1     2.3E-1  4.2E-1
!64     0.015625   1.7E-2  2.1E-2   4.7E-2  8.7E-2     1.1E-1  2.7E-1
!128    0.0078125  7.9E-3  9.8E-3   2.3E-2  4.1E-2     5.4E-2  2.5E-1
!256    0.00390625 4.0E-3  4.9E-3   1.1E-2  2.1E-2     2.8E-2  2.2E-1
!average flux:(16) -9.06 (32) -9.31 (64) -9.46 (128) -9.54 (256) 9.58
!
! kratio=1
! SIMPLE (corner):
!            dx    L1         L2      Linf   L1 grad    L2 grad Linf grad  
!16     0.0625     6.9E-2  7.5E-2   1.0E-1  4.8E-1     7.7E-1  3.3
!32     0.03125    5.2E-2  5.5E-2   9.4E-2  3.5E-1     6.1E-1  3.4
!64     0.015625   2.3E-2  2.5E-2   5.0E-2  1.8E-1     3.8E-1  3.5
!128    0.0078125  1.2E-2  1.2E-2   2.5E-2  1.0E-1     2.9E-1  3.4
!256    0.00390625 5.8E-3  6.2E-3   1.3E-2  5.6E-2     2.0E-1  3.7
!average flux:(16) -8.95 (32) -9.73 (64) -9.90 (128) -9.99 (256) -10.03

!
! YANG: DO NOT ADD THIS DATA
! DS:
!            dx    L1         L2      Linf   L1 grad    L2 grad Linf grad  
!16     0.0625     5.0E-2  6.0E-2   1.1E-1  2.8E-1     3.8E-1  6.4E-1
!32     0.03125    4.7E-2  5.9E-2   1.5E-1  3.2E-1     4.1E-1  1.1
!64     0.015625   2.1E-2  2.6E-2   6.6E-2  1.2E-1     1.5E-1  4.9E-1
!128    0.0078125  1.1E-2  1.4E-2   4.8E-2  7.0E-3     1.2E-1  2.3
!256    0.00390625 5.7E-3  6.9E-3   2.0E-2  3.3E-2     5.5E-2  9.3E-1
!average flux:(16) -9.00 (32) -9.10 (64) -9.42 (128)-9.51  (256) -9.56
!
! YANG: DO NOT ADD THIS DATA
! OP:
!            dx    L1         L2      Linf   L1 grad    L2 grad Linf grad  
!16     0.0625     4.0E-2  4.9E-2   1.1E-1  2.5E-1     3.5E-1  5.5E-1
!32     0.03125    2.5E-2  3.1E-2   7.0E-2  1.5E-1     2.0E-1  3.7E-1
!64     0.015625   1.4E-2  1.7E-2   4.3E-2  7.8E-2     1.0E-1  2.1E-1
!128    0.0078125  6.9E-3  8.7E-3   2.2E-2  3.9E-2     5.1E-2  1.6E-1
!256    0.00390625 3.5E-3  4.4E-3   1.1E-2  2.0E-2     2.6E-2  1.5E-1
!average flux:(16) -9.01 (32) -9.36 (64) -9.48 (128) -9.55 (256) -9.58
!
! YANG: DO NOT ADD THIS DATA
! "LE":
!            dx    L1         L2      Linf   L1 grad    L2 grad Linf grad  
!16     0.0625     1.9E-1  1.9E-1   2.4E-1  2.4E-1     3.1E-1  5.6E-1
!32     0.03125    7.7E-2  8.2E-2   1.6E-1  2.1E-1     2.7E-1  9.8E-1
!64     0.015625   3.3E-3  3.9E-3   9.5E-3  2.3E-2     4.2E-2  4.5E-1
!128    0.0078125  3.2E-4  4.5E-4   3.6E-3  6.3E-3     1.6E-2  2.7E-1
!256    0.00390625 1.9E-4  2.6E-4   1.6E-3  3.6E-3     8.9E-3  2.1E-1
!average flux:(16) -9.06 (32) -9.34 (64) -9.59 (128) -9.61 (256) -9.61
!
! YANG LIU CHECK TABLE IF AGREEMENT
! probtype==15 (hypocycloid 2 mat), 2000 markers
! ERRTOL=0.999999D0, 
! TSTOP=1.25E-2, operator_ext=1 operator_int=1 linear exact=0,
! timesteps refined by factor of 2.  32(1), 64(2), 128(4), 256(8), 512(16)
!        L1         L2   Linf   L1 grad    L2 grad Linf grad  
!32/64   1.4E-1  2.1E-1  5.5E-1 2.7        3.7     8.8              
!64/128  1.1E-1  1.6E-1  4.1E-1 3.0        3.9     7.7              
!128/256 7.2E-2  9.7E-2  2.1E-1 1.8        2.2     4.3              
!256/512 4.0E-2  5.3E-2  1.1E-1 1.0        1.2     2.5
!
! YANG LIU CHECK TABLE IF AGREEMENT
! probtype==15 (hypocycloid 2 mat), 4000 markers
! ERRTOL=0.999999D0, 
! TSTOP=1.25E-2, operator_ext=1 operator_int=1 linear exact=1,
! timesteps refined by factor of 2.  32(1), 64(2), 128(4), 256(8), 512(16)
!        L1         L2   Linf   L1 grad    L2 grad Linf grad  
!32/64   1.4E-1  2.1E-1  5.5E-1 2.7        3.7     8.8              
!64/128  1.1E-1  1.6E-1  4.1E-1 3.0        3.9     7.8              
!128/256 7.2E-2  9.7E-2  2.1E-1 1.8        2.2     4.3              
!256/512 4.0E-2  5.3E-2  1.1E-1 1.0        1.2     2.6
!
! YANG LIU CHECK TABLE IF AGREEMENT
! probtype==20 (hypocycloid 6 mat), 2000 markers,type 3
! ERRTOL=0.999999D0, 
! TSTOP=1.25E-2, operator_ext=1 operator_int=1 linear exact=0,
! timesteps refined by factor of 2.  32(1), 64(2), 128(4), 256(8), 512(16)
!        L1         L2   Linf   L1 grad    L2 grad Linf grad  
!32/64   4.3E-1  4.3E-1  5.5E-1 3.1        3.2      4.2             
!64/128  2.7E-1  2.7E-1  3.1E-1 1.2        1.4      2.3
!128/256 1.5E-1  1.5E-1  1.8E-1 0.79       0.89     1.8              
!256/512 7.9E-2  8.1E-2  1.0E-1 5.2E-1     5.9E-1   1.4              
!
! YANG LIU CHECK TABLE IF AGREEMENT
! probtype==20 (hypocycloid 6 mat), 2000 markers,type 2
! ERRTOL=0.999999D0, 
! TSTOP=1.25E-2, operator_ext=1 operator_int=1 linear exact=0,
! timesteps refined by factor of 2.  32(1), 64(2), 128(4), 256(8), 512(16)
!        L1         L2   Linf   L1 grad    L2 grad Linf grad  
!32/64   1.9E-1  2.5E-1  5.5E-1 5.7        6.0      8.8             
!64/128  1.3E-1  1.8E-1  4.1E-1 4.2        4.9      7.8
!128/256 8.1E-2  1.0E-1  2.1E-1 2.2        2.5      4.3              
!256/512 4.5E-2  5.6E-2  1.1E-1 1.1        1.3      2.5              
!
! probtype==1, dirichlet bc, ERRTOL=0.999999D0, thick, TSTOP=1.25E-2
! operator_ext=1 operator_int=3 linear exact=1
!        L1         L2   Linf   L1 grad    L2 grad Linf grad  Linf flux
! 32(1)  2.6E-3 2.9E-3 5.5E-3   1.6E-2     2.1E-2   6.6E-2    1.2E-1
! 64(4)  7.1E-4 8.2E-4 1.7E-3   5.7E-3     7.3E-3   2.8E-2    4.8E-2 
!128(16) 1.7E-4 2.0E-4 4.2E-4   1.6E-3     2.1E-3   8.7E-3    2.4E-2
!256(64) 4.3E-5 5.1E-5 1.1E-4   4.2E-4     5.5E-4   2.7E-3    1.3E-2 
!
! probtype==1, dirichlet bc, ERRTOL=0.01D0, thin, TSTOP=1.25E-2
! operator_ext=1 operator_int=3 linear exact=1
!        L1         L2   Linf   L1 grad    L2 grad Linf grad 
! 32(1) 4.2E-4  6.9E-4   1.8E-3  3.8E-2    1.0E-1  5.4E-1
! 64(2) 5.4E-4  7.1E-4   1.9E-3  6.5E-2    1.0E-1  4.5E-1
!128(4) 1.8E-4  2.5E-4   6.6E-4  3.6E-2    5.4E-2  1.9E-1
!256(8) 3.0E-6  3.6E-6   7.5E-6  4.7E-4    8.0E-4  3.6E-3
!
! probtype==1, dirichlet bc, ERRTOL=0.01D0, thin, TSTOP=1.25E-2
! operator_ext=1 operator_int=1 linear exact=0
!        L1         L2   Linf   L1 grad    L2 grad Linf grad 
! 32(1) 1.1E-3  1.4E-3   3.3E-3  7.6E-2    1.3E-1  5.6E-1
! 64(2) 4.4E-4  6.3E-4   1.7E-3  6.2E-2    9.5E-2  4.0E-1
!128(4) 2.0E-4  3.1E-4   1.2E-3  4.2E-2    6.0E-2  2.1E-1
!256(8) 8.6E-5  1.3E-4   6.2E-4  3.3E-2    5.4E-2  3.3E-1
!
!
! probtype==1, thin
! ERRTOL=0.01D0, 2*radeps=0.01, Neumann, T=1.25E-2, linear exact
!        dx          L1      L2   Linf   L1 grad  L2 grad Linf grad 
! 32   1 0.03125    1.1E-2  1.3E-2 2.5E-2  2.5E-1  3.9E-1   1.2
! 64   2 0.015625   1.0E-2  1.3E-2 2.7E-2  2.7E-1  3.4E-1   9.7E-1
! 128  4 0.0078125  8.4E-3  1.0E-2 2.4E-2  1.4E-1  1.9E-1   5.9E-1
! 256  8 0.00390625 1.3E-3  1.5E-3 2.1E-3  4.3E-3  5.8E-3   5.5E-2
!
! probtype==1
! ERRTOL=0.999999, thick, Neumann, T=1.25E-2
! T=1.25E-2    L1      L2     LINF    L1_grd  L2_grd   LINF_grd
! 32   1 step  9.6E-3  1.1E-2 1.6E-2  2.5E-2  3.2E-2   8.4E-2
! 64   4 step  2.6E-3  2.9E-3 4.4E-3  6.8E-3  8.8E-3   2.7E-2
! 128 16 step  6.7E-4  7.5E-4 1.1E-3  1.8E-3  2.3E-3   1.1E-2
! 256 64 step  1.7E-4  1.9E-4 2.8E-4  4.5E-4  5.8E-4   3.1E-3
!
! YANG LIU MODIFY TABLE (MAKE GRAPH TOO!)
! probtype=19, ERRTOL=0.999999D0, polar solver: 128x256, tstop=0.004
! T=0.004        L1      L2     LINF    L1_grd  L2_grd   LINF_grd Linf_int_grd
! 32   1 step  6.9E-2  8.0E-2  1.5E-1   4.7E-1   6.2E-1 2.1       2.7
! 64   4 step  2.1E-2  2.4E-2  4.9E-2   1.8E-1   2.4E-1 7.5E-1    8.3E-1
! 128 16 step  5.3E-3  6.4E-3  1.4E-2   5.2E-2   6.8E-2 2.2E-1    2.5E-1
! 256 64 step  1.3E-3  1.6E-3  3.5E-3   1.4E-2   1.8E-2 6.1E-2    7.4E-2
!
! YANG LIU MODIFY TABLE 
! probtype=19, ERRTOL=0.999999D0, polar solver: 256x512, tstop=0.004
! T=0.004        L1      L2     LINF    L1_grd  L2_grd   LINF_grd Linf_int_grd
! 32   1 step  6.9E-2  8.0E-2  1.5E-1   4.7E-1   6.2E-1 2.1       2.7
! 64   4 step  2.1E-2  2.4E-2  5.0E-2   1.8E-1   2.4E-1 7.5E-1    8.4E-1
! 128 16 step  5.3E-3  6.4E-3  1.3E-2   5.2E-2   6.8E-2 2.2E-1    2.6E-1
! 256 64 step  1.3E-3  1.6E-3  3.4E-3   1.4E-2   1.8E-2 6.0E-2    8.4E-2
!
! 0=flat interface  
! 1=annulus  
! 2=vertical interface
! 3=expanding circle (material 1 inside, T=TSAT inside)
! 4=expanding or shrinking circle (material 1 inside, T=TSAT outside)
!    (TSTOP=1.25D-3)
! 5=phase change for vertical planar interface (initial location xblob_in)
!    (TSTOP=0.5d0)
! 13=pentafoil (pentaeps, dirichlet_pentafoil init below) 
! 16=multiscale (thin filament effect)
!    (filament thickness is thermal_delta init below)
! 19=polar solver
! 15=hypocycloid with 2 materials
! 20=hypocycloid with 6 materials
! 400=melting gingerbread (material 1 inside, T=TSAT initially)
! 401=ice melt (material 1 liquid, material 2 gas, material 3 ice)
INTEGER,PARAMETER          :: probtype_in=401
INTEGER,PARAMETER          :: stefan_flag=1
! 0.1 if probtype_in=3  0.4 if probtype_in=4
real(kind=8),PARAMETER     :: radblob_in = 0.4d0
! buffer for probtype_in=3
real(kind=8),PARAMETER     :: radblob2_in = 0.05d0  
real(kind=8),PARAMETER     :: xblob_in = 0.2d0
real(kind=8),PARAMETER     :: yblob_in = 0.5d0
! for probtype=16 , top and bot temperature profile
real(kind=8),parameter     :: NB_top=0.0d0, NB_bot=10.0d0  
! 1.0d0 for probtype==3
! -4.0d0 for probtype==4 (TDIFF=T_DISK_CENTER - TSAT)
! 10.0d0 for probtype==16
! for probtype==5, T(x=0)=273.0  T(x=1)=272.0+e^{-V(1-Vt)}
! material 1 on the left, material 2 on the right.
! 1.0d0 for probtype==400 (gingerbread)
!  (T=TSAT interior domain initially, T=TSAT+TDIFF on walls)
! 10.0d0 for probtype==401 (ice melt)
real(kind=8),PARAMETER     :: TDIFF_in = 10.0d0
! 10.0d0 for probtype==3
! 1.0d0 for probtype==4 (stationary benchmark)
! 1.0d0 for probtype==4 (shrinking material 1)
! 1.0d0 for probtype==400 (melting gingerbread)
real(kind=8),PARAMETER     :: latent_heat_in = 1.0d0
!0=low,1=simple,2=Dai and Scannapieco,3=orthogonal probe
INTEGER,PARAMETER          :: local_operator_internal = 3
INTEGER,PARAMETER          :: local_operator_external = 1
INTEGER,PARAMETER          :: local_linear_exact = 1
INTEGER                    :: ilev,max_ncell
INTEGER                    :: N_START,N_FINISH,N_CURRENT
! M=1 non-deforming boundary tests
! M=40 probtype_in=3 test with N=64
INTEGER                    :: M_START,M_FACTOR,M_CURRENT
INTEGER,PARAMETER          :: M_MAX_TIME_STEP = 2000
INTEGER,PARAMETER          :: plot_int = 1
! TSTOP=1.25d-2 for probtype_in=1 (annulus)
! TSTOP=1.25d-2 for probtype_in=13,15,20 (pentafoil, Hypocycloid)
! explicit time step for N=512 grid: 4 dt/dx^2 < 1
! dt<dx^2/4=1/(4 * 512^2)=.95E-6
! for N=256, explicit time step=3.8E-6
!
! multiscale: (probtype_in.eq.16)
! TSTOP=2.0d0  if filament thickness is 0.02.
!
! non-axisymmetric, polar solver for validation (probtype_in.eq.19):
! TSTOP=0.004d0
real(kind=8),parameter     :: TSTOP = 0.5d0
! fixed_dt=0.0d0 => use CFL condition
! fixed_dt=-1.0d0 => use TSTOP/M
real(kind=8)               :: fixed_dt_main,fixed_dt_current
real(kind=8),parameter     :: CFL = 0.5d0
real(kind=8),parameter     :: problo= 0.0d0, probhi= 1.0d0
integer,parameter          :: sdim_in = 2

INTEGER :: nmax
INTEGER :: nmat_in
INTEGER :: precond_type_in
INTEGER :: dir
INTEGER :: side
REAL(kind=8) :: xcen,ycen
REAL(kind=8) :: xcen_vec(2)
REAL(kind=8) :: time_init,xgrid,ygrid
REAL(kind=8) :: deltat_in
REAL(kind=8) :: deltat_polar
INTEGER      :: subcycling_step
REAL(kind=8) :: bicgstab_tol_in
REAL(kind=8) :: current_time_in
REAL(kind=8) :: alpha_in(100)

INTEGER                    :: i,j,tm
REAL(KIND=8)               :: h_in
REAL(KIND=8)               :: time_n,time_np1
!REAL(KIND=8),dimension(-1:N+1) :: XLINE,YLINE 
REAL(KIND=8),dimension(:), allocatable :: XLINE,YLINE ! nodes
!real(kind=8),dimension(-1:N) :: xCC,yCC       
real(kind=8),dimension(:), allocatable :: xCC,yCC       ! cell centers
!REAL(KIND=8)               :: Ts(M+1)
REAL(KIND=8),dimension(:), allocatable :: Ts
real(kind=8)               :: dx_in(sdim_in)
real(kind=8)               :: dx_local(sdim_in)
real(kind=8)               :: dx_coarse
!TYPE(POLYGON),dimension(-1:N,-1:N):: CELL_FAB
TYPE(POLYGON),dimension(:,:), allocatable :: CELL_FAB
real(kind=8),external      :: exact_temperature
real(kind=8)               :: max_front_vel
real(kind=8)               :: test_vel
real(kind=8)               :: lmSt
real(kind=8)               :: rstefan
real(kind=8)               :: T_FIELD
real(kind=8)               :: stefan_time
real(kind=8)               :: local_vof

real(kind=8)                :: xsten_cache(-1:1)
integer                     :: nhalf
integer                     :: imof
integer                     :: im
integer                     :: im1
integer                     :: im_opp
real(kind=8)                :: sumT,sumvf,sumvf2,voltotal,local_Pi
real(kind=8)                :: eff_radius
real(kind=8)                :: expect_radius
real(kind=8)                :: test_radblob

!---------------------------------------------------
integer                     :: local_nten
integer                     :: iten

!----------------------------------------
INTEGER order_algorithm(1000)
INTEGER MOFITERMAX
INTEGER ngeom_recon_in
INTEGER bfactmax
INTEGER domlo_in(2)
INTEGER domhi_in(2)
real(kind=8) :: problo_arr(2)

integer local_state_ncomp
integer nx_in,ny_in,lox_in,loy_in,hix_in,hiy_in
integer hflag
integer vofcomp
integer vofcomp2
integer scomp
integer nsteps
integer total_nsteps_parm
integer inode
real(kind=8) :: cc(2)
real(kind=8) :: dtemp1,dtemp2
real(kind=8) :: sum_alpha
real(kind=8) :: flxavg1,flxavg2
real(kind=8) :: flxtot1,flxtot2
real(kind=8) :: xlo_fluxtest,xhi_fluxtest
real(kind=8) :: y_fluxtest1
real(kind=8) :: y_fluxtest2
real(kind=8) :: LL
real(kind=8) :: TSAT

integer j_fluxtest,ilo_fluxtest,ihi_fluxtest,isum
integer icen,jcen

integer ireverse,isink

real(kind=8), dimension(:,:,:), allocatable :: UNEW_in
real(kind=8), dimension(:,:,:), allocatable :: UOLD_in
real(kind=8), dimension(:,:,:), allocatable :: beta_in
real(kind=8), dimension(:,:,:), allocatable :: VFRAC_MOF_in

! -1:N,-1:N,nmat
real(kind=8),dimension(:,:,:),allocatable :: vf
! -1:N,-1:N,nmat*ngeom_recon
real(kind=8),dimension(:,:,:),allocatable :: mofdata_FAB_in
! -1:N,-1:N,nmat
TYPE(POINTS),DIMENSION(:,:,:),allocatable :: CENTROID_FAB
! -1:N,-1:N,nmat,sdim
real(kind=8),dimension(:,:,:,:),allocatable :: centroid_mult   
! -1:N,-1:N,nmat
real(kind=8),dimension(:,:,:),allocatable :: T
real(kind=8),dimension(:,:,:),allocatable :: T_new

integer N_COARSE
real(kind=8) :: err1T,err2T,err3T
real(kind=8) :: err1T_gradient,err2T_gradient,err3T_gradient
real(kind=8) :: local_interp
real(kind=8) :: fine_gradientx,fine_gradienty
real(kind=8) :: coarse_gradientx,coarse_gradienty
real(kind=8) :: gradient_err
integer count1,count1_gradient
integer cni,cnj
integer vfcheck,gradient_check
real(kind=8),dimension(:,:,:),allocatable :: coarse_data
real(kind=8),dimension(:,:,:),allocatable :: fine_data
integer :: im_measure
integer :: constant_K_test
integer :: iter
integer :: finished_flag
real(kind=8) :: iter_average

integer :: sci_max_level

print *,"PROTOTYPE CODE DATE= Feb 29, 2020, 15:40pm"


im_measure=2
constant_K_test=0

print *,"im_measure= ",im_measure
print *,"constant_K_test= ",constant_K_test

! N space
! M time
N_START=64
N_FINISH=64
M_START=64
M_FACTOR=2

if (probtype_in.eq.4) then
 fixed_dt_main=-1.0d0 ! dt=1.25D-4 N=64  M=10  TSTOP=1.25D-3
else if ((probtype_in.eq.13).or. & ! hypocycloid
         (probtype_in.eq.15).or. &
         (probtype_in.eq.20)) then
 fixed_dt_main=-1.0d0 ! dt=1.25D-2 N=32  M=1  TSTOP=1.25D-2
else if (probtype_in.eq.16) then ! multiscale
 fixed_dt_main=-1.0d0 ! TSTOP=2.0
else if (probtype_in.eq.1) then ! annulus
 fixed_dt_main=-1.0d0 ! TSTOP=1.25D-2 
else if (probtype_in.eq.19) then ! polar solver
 fixed_dt_main=-1.0d0 ! TSTOP=0.004
else if (probtype_in.eq.0) then ! flat interface
 fixed_dt_main=-1.0d0 
else if (probtype_in.eq.2) then ! vertical
 fixed_dt_main=-1.0d0 
else if (probtype_in.eq.3) then ! expanding circle
 fixed_dt_main=-1.0d0  ! TSTOP=1.25D-3
else if (probtype_in.eq.4) then ! expanding or shrinking circle
 fixed_dt_main=-1.0d0  ! TSTOP=1.25D-3
else if (probtype_in.eq.5) then ! phase change vertical planar interface
 fixed_dt_main=-1.0d0  ! TSTOP=0.5d0
else if (probtype_in.eq.400) then ! gingerbread man
 fixed_dt_main=0.0d0
 print *,"gingeroutline should be in run directory"
else if (probtype_in.eq.401) then ! melting block of ice
 fixed_dt_main=0.0d0
else
 print *,"probtype_in invalid"
 stop
endif

if (fixed_dt_main.eq.-1.0d0) then
 ! do nothing fixed_dt_current=TSTOP/M_CURRENT
else if (fixed_dt_main.eq.0.0d0) then
 ! do nothing fixed_dt_current=0.0d0
else if (fixed_dt_main.gt.0.0d0) then
 ! do nothing fixed_dt_current=fixed_dt_main
else
 print *,"fixed_dt_main invalid"
 stop
endif

! INITIALIZE VARIABLES DECLARED IN vof_cisl.F90 (Module GeneralClass)

! r1=radcen-radeps
! r2=radcen+radeps
radcen=0.25d0
radeps=0.005d0  ! ! thick:0.1d0  thin:0.005d0
rlo=radcen-radeps
rhi=radcen+radeps
! radial_variation=0 for thin annulus Dirichlet test problem.
radial_variation=0
dirichlet_pentafoil=1
dirichlet_annulus=1

pentaeps=0.05d0 ! thick: 0.2d0  thin: 0.05d0 penta_foil
! for probtype=16   thickness of the thermal layer
! 0.001d0 is thin.  0.02d0 is thicker.
thermal_delta= 0.02d0  
! filament_test_type==0 for irregular material 2
! filament_test_type==1 for circular material 2
filament_test_type = 0

! INITIALIZE VARIABLES DECLARED IN BICGSTAB_Yang_MULTI.F90 
! (Module bicgstab_module)
! 0.01D0, 0.5D0, 0.99D0, 0.99999999D0 are some options
!ERRTOL=0.99999999D0
ERRTOL=0.01D0



VOFTOL_local=VOFTOL

local_Pi=4.0d0*atan(1.0d0)

! pcurve_ls defined in vof_cisl.F90
pcurve_ls = 0.0d0
call asteroidshape(pcurve_ls)
do i=1,pcurve_num+1
   pcurve_ls(1,i)=(pcurve_ls(1,i)+1.0d0)/2.0d0
   pcurve_ls(2,i)=(pcurve_ls(2,i)+1.0d0)/2.0d0
enddo

problox=problo
probloy=problo
probhix=probhi
probhiy=probhi
probloz=0.0d0
probhiz=0.0d0
problenx=probhix-problox
probleny=probhiy-probloy
problenz=probhiz-probloz

sci_max_level=0
fort_max_level=0
fort_finest_level=0
do im=1,100
  FSI_flag(im)=0
enddo

if (local_linear_exact.eq.1) then
 if ((local_operator_internal.eq.3).and. &
     (local_operator_external.eq.1)) then
  ! do nothing
 else
  print *,"local_operator_internal or local_operator_external bad"
  stop
 endif
else if (local_linear_exact.eq.0) then
 ! check nothing
else
 print *,"local_linear_exact invalid"
 stop
endif

print *,"N_START,N_FINISH ",N_START,N_FINISH
print *,"M_START,M_FACTOR ",M_START,M_FACTOR
print *,"fixed_dt_main= ",fixed_dt_main

N_CURRENT=N_START
M_CURRENT=M_START

DO WHILE (N_CURRENT.le.N_FINISH)

 if (fixed_dt_main.eq.-1.0d0) then
  fixed_dt_current=TSTOP/M_CURRENT
 else if (fixed_dt_main.ge.0.0d0) then
  fixed_dt_current=fixed_dt_main
 else
  print *,"fixed_dt_main invalid"
  stop
 endif

 nmax=POLYGON_LIST_MAX
 levelrz=0
 probtype=probtype_in  ! defined in probdataf95.H (probcommon)
 radblob=radblob_in
 radblob2=radblob2_in
 xblob=xblob_in
 yblob=yblob_in
 order_algorithm = 0

 max_front_vel=0.0

 isink=0

 do dir=1,sdim_in
   do side=1,2
    physbc(dir,side)=REFLECT_EVEN
    physbc_value(dir,side)=0.0d0
   enddo
 enddo

 if (probtype_in.eq.0) then

   nmat_in=2
   fort_heatviscconst(1)=1.0
   fort_heatviscconst(2)=0.1
   physbc(2,1)=EXT_DIR
   physbc_value(2,1)=3.0d0
   physbc(2,2)=EXT_DIR
   physbc_value(2,2)=2.0d0

 else if (probtype_in.eq.1) then

   order_algorithm(1)=1
   order_algorithm(2)=3
   order_algorithm(3)=2
   nmat_in=3
   fort_heatviscconst(1)=0.0 
   fort_heatviscconst(2)=1.0 
   fort_heatviscconst(3)=0.0 
   do dir=1,sdim_in
   do side=1,2
    physbc(dir,side)=EXT_DIR
    physbc_value(dir,side)=0.0d0
   enddo
   enddo

 else if (probtype_in.eq.2) then
   nmat_in=2
   fort_heatviscconst(1)=1.0
   fort_heatviscconst(2)=0.1
   physbc(1,1)=EXT_DIR
   physbc_value(1,1)=3.0d0
   physbc(1,2)=EXT_DIR
   physbc_value(1,2)=2.0d0
 else if (probtype_in.eq.3) then
   nmat_in=2
   fort_heatviscconst(1)=1.0
   fort_heatviscconst(2)=1.0
   do dir=1,sdim_in
   do side=1,2
    physbc(dir,side)=REFLECT_EVEN
    physbc_value(dir,side)=0.0
   enddo
   enddo

   ! distance function: see subroutine dist_concentric in multimat_FVM.F90
 else if (probtype_in.eq.4) then
   nmat_in=2
   fort_heatviscconst(1)=1.0  ! inside material (supercooled liquid)
   fort_heatviscconst(2)=0.0  ! outside (ice)
   do dir=1,sdim_in
   do side=1,2
    physbc(dir,side)=REFLECT_EVEN
    physbc_value(dir,side)=0.0
   enddo
   enddo
   if ((stefan_flag.eq.1).and. &
       (local_operator_internal.eq.3).and. &
       (local_operator_external.eq.1).and. &
       (local_linear_exact.eq.1)) then
    ! do nothing
   else
    print *,"stefan_flag,op int,op ext,or local_linear_exact invalid prob==4"
    stop
   endif

 else if (probtype_in.eq.400) then

   sci_max_level=2
   nmat_in=2
   fort_heatviscconst(1)=1.0  ! inside gingerbread (T=TSAT init.)
   fort_heatviscconst(2)=1.0  ! outside gingerbread (T=TSAT init.)
   do dir=1,sdim_in
   do side=1,2
    physbc(dir,side)=REFLECT_EVEN
    physbc_value(dir,side)=0.0
   enddo
   enddo
   physbc(2,1)=EXT_DIR
   physbc_value(2,1)=273.0d0+TDIFF_in
   physbc(2,2)=EXT_DIR
   physbc_value(2,2)=273.0d0+TDIFF_in
   if ((stefan_flag.eq.1).and. &
       (local_operator_internal.eq.3).and. &
       (local_operator_external.eq.1).and. &
       (local_linear_exact.eq.1)) then
    ! do nothing
   else
    print *,"stefan_flag,op int,op ext,or local_linear_exact invalid prob==400"
    stop
   endif

   FSI_flag(1)=7 ! gingerbread (in the man)

 else if (probtype_in.eq.401) then ! melting block of ice

   nmat_in=3
   fort_heatviscconst(1)=10.0  ! liquid
   fort_heatviscconst(2)=1.0  ! gas
   fort_heatviscconst(3)=100.0  ! ice
   do dir=1,sdim_in
   do side=1,2
    physbc(dir,side)=REFLECT_EVEN
    physbc_value(dir,side)=0.0
   enddo
   enddo
   physbc(2,1)=EXT_DIR
   physbc_value(2,1)=273.0d0+TDIFF_in
   if ((stefan_flag.eq.1).and. &
       (local_operator_internal.eq.3).and. &
       (local_operator_external.eq.1).and. &
       (local_linear_exact.eq.1)) then
    ! do nothing
   else
    print *,"stefan_flag,op int,op ext,or local_linear_exact invalid prob==400"
    stop
   endif

 else if (probtype_in.eq.5) then

   nmat_in=2
   fort_heatviscconst(1)=0.0  ! inside material (left material)
   fort_heatviscconst(2)=1.0  ! outside (right material)
   do dir=1,sdim_in
   do side=1,2
    physbc(dir,side)=REFLECT_EVEN
    physbc_value(dir,side)=0.0
   enddo
   enddo
   dir=1
   side=1
   physbc(dir,side)=EXT_DIR
   physbc_value(dir,side)=273.0d0

   if ((stefan_flag.eq.1).and.(local_linear_exact.eq.1)) then
    ! do nothing
   else
    print *,"stefan_flag or local_linear_exact invalid for probtype==4"
    stop
   endif

 else if(probtype_in.eq.19)then

   do dir=1,sdim_in
   do side=1,2
    physbc(dir,side)=REFLECT_EVEN
    physbc_value(dir,side)=0.0
   enddo
   enddo
   physbc(1,1)=EXT_DIR
   physbc_value(1,1)=0.0d0
   physbc(1,2)=EXT_DIR
   physbc_value(1,2)=0.0d0
   physbc(2,1)=EXT_DIR
   physbc_value(2,1)=0.0d0
   physbc(2,2)=EXT_DIR
   physbc_value(2,2)=0.0d0

   order_algorithm(1)=1
   order_algorithm(2)=3
   order_algorithm(3)=2
   nmat_in=3
   fort_heatviscconst(1)=0.0d0
   fort_heatviscconst(2)=1.0d0
   fort_heatviscconst(3)=0.0d0

 else if(probtype_in.eq.13)then
   physbc(1,1)=EXT_DIR
   physbc_value(1,1)=0.0d0
   physbc(1,2)=EXT_DIR
   physbc_value(1,2)=0.0d0
   physbc(2,1)=EXT_DIR
   physbc_value(2,1)=0.0d0
   physbc(2,2)=EXT_DIR
   physbc_value(2,2)=0.0d0
   order_algorithm(1)=1
   order_algorithm(2)=3
   order_algorithm(3)=2
   nmat_in=3
   fort_heatviscconst(1) = 0.0d0
   fort_heatviscconst(2) = 1.0d0
   fort_heatviscconst(3) = 0.0d0

 elseif(probtype_in.eq.14)then
   physbc(1,1)=EXT_DIR
   physbc_value(1,1)=0.0d0
   physbc(1,2)=EXT_DIR
   physbc_value(1,2)=0.0d0
   physbc(2,1)=EXT_DIR
   physbc_value(2,1)=0.0d0
   physbc(2,2)=EXT_DIR
   physbc_value(2,2)=0.0d0
   nmat_in=2 
   fort_heatviscconst(1) = 1.0d0           ! interior region
   fort_heatviscconst(2) = 2.0d0          ! exterior region
 elseif(probtype_in.eq.15)then
   physbc(1,1)=EXT_DIR
   physbc_value(1,1)=10.0d0
   physbc(1,2)=EXT_DIR
   physbc_value(1,2)=10.0d0
   physbc(2,1)=EXT_DIR
   physbc_value(2,1)=10.0d0
   physbc(2,2)=EXT_DIR
   physbc_value(2,2)=10.0d0
   nmat_in=2
   fort_heatviscconst(1) = 0.1d0           ! interior region
   fort_heatviscconst(2) = 10.0d0          ! exterior region
 else if (probtype_in.eq.16) then
   do dir=1,sdim_in
   do side=1,2
    physbc(dir,side)=REFLECT_EVEN
    physbc_value(dir,side)=0.0
   enddo
   enddo
   physbc(2,1)=EXT_DIR
   physbc_value(2,1)=NB_bot
   physbc(2,2)=EXT_DIR
   physbc_value(2,2)=NB_top
   nmat_in=3
   order_algorithm(1)=1
   order_algorithm(2)=2
   order_algorithm(3)=3
   if (constant_K_test.eq.1) then
    fort_heatviscconst(1) = 1.0d0
    fort_heatviscconst(2) = 1.0d0
    fort_heatviscconst(3) = 1.0d0 ! 0.001d0 for standard case
   else if (constant_K_test.eq.0) then
    fort_heatviscconst(1) = 1.0d0
    fort_heatviscconst(2) = 1.0d0
    fort_heatviscconst(3) = 0.001d0 ! 0.001d0 for standard case
   else
    print *,"constant_K_test invalid"
    stop
   endif

   print *,"filament_test_type= ",filament_test_type
 elseif(probtype_in.eq.17)then
   nmat_in=5
   physbc(1,1)=EXT_DIR
   physbc_value(1,1)=0.0d0
   physbc(1,2)=EXT_DIR
   physbc_value(1,2)=0.0d0
   physbc(2,1)=EXT_DIR
   physbc_value(2,1)=0.0d0
   physbc(2,2)=EXT_DIR
   physbc_value(2,2)=0.0d0
   order_algorithm(1)=1
   order_algorithm(2)=2
   order_algorithm(3)=3
   order_algorithm(4)=4
   order_algorithm(5)=5
   fort_heatviscconst(1) = 1.0d0
   fort_heatviscconst(2) = 0.1d0
   fort_heatviscconst(3) = 1.0d0
   fort_heatviscconst(4) = 0.1d0
   fort_heatviscconst(5) = 1.0d0
 elseif(probtype_in.eq.20)then
   physbc(1,1)=EXT_DIR
   physbc_value(1,1)=10.0d0
   physbc(1,2)=EXT_DIR
   physbc_value(1,2)=10.0d0
   physbc(2,1)=EXT_DIR
   physbc_value(2,1)=10.0d0
   physbc(2,2)=EXT_DIR
   physbc_value(2,2)=10.0d0
   nmat_in=6
   order_algorithm(1)=1
   order_algorithm(6)=2
   fort_heatviscconst(1) = 10.0d0
   fort_heatviscconst(6) = 0.1d0

   if (1.eq.0) then  ! type 3
    fort_heatviscconst(2) = 0.5d0  ! 0.1=type ii   0.5=type iii
    fort_heatviscconst(3) = 1.0d0  ! 0.1=type ii   1.0=type iii
    fort_heatviscconst(4) = 0.5d0  ! 0.1=type ii   0.5=type iii
    fort_heatviscconst(5) = 1.0d0  ! 0.1=type ii   1.0=type iii
   else if (1.eq.1) then ! type 2
    fort_heatviscconst(2) = 0.1d0  ! 0.1=type ii   0.5=type iii
    fort_heatviscconst(3) = 0.1d0  ! 0.1=type ii   1.0=type iii
    fort_heatviscconst(4) = 0.1d0  ! 0.1=type ii   0.5=type iii
    fort_heatviscconst(5) = 0.1d0  ! 0.1=type ii   1.0=type iii
   else
    print *,"must be type 3 or type 2"
    stop
   endif

 else
   print *,"probtype_in invalid1 ",probtype_in
   stop
 endif

 num_state_material=2
 normal_probe_size=1
 num_materials=nmat_in
 local_nten=( (nmat_in-1)*(nmat_in-1)+nmat_in-1 )/2
 global_nten=local_nten

 ngrow_expansion=2
 num_species_var=0

 bfact_time_order=1
 bfact_space_order(0)=1
 bfact_space_order(1)=1

 bfactmax=8

 call sanity_check(bfactmax+2)

 call init_cache(bfactmax+2)

 cache_max_level=fort_max_level
 if (cache_max_level.lt.sci_max_level) then
  cache_max_level=sci_max_level
 endif

 max_ncell=N_CURRENT
 do ilev=1,cache_max_level
  max_ncell=max_ncell*2
 enddo

 cache_index_low=-4*bfactmax
 cache_index_high=2*max_ncell+4*bfactmax

 allocate(grid_cache(0:cache_max_level, &
    cache_index_low:cache_index_high,sdim_in))

 h_in = (probhi-problo)/N_CURRENT
 do dir=1,sdim_in
  dx_in(dir) = h_in
 enddo

 do dir=1,sdim_in
  domlo_in(dir)=0
  domhi_in(dir)=N_CURRENT-1
 enddo

 allocate(dxlevel(0:cache_max_level,sdim_in))
 allocate(domlo_level(0:cache_max_level,sdim_in))
 allocate(domhi_level(0:cache_max_level,sdim_in))

 do dir=1,sdim_in
  dxlevel(0,dir)=dx_in(dir)
  domlo_level(0,dir)=domlo_in(dir)
  domhi_level(0,dir)=domhi_in(dir)
 enddo

 do ilev=0,cache_max_level
  do dir=1,sdim_in
   dx_local(dir)=dxlevel(ilev,dir)
  enddo
  do dir=1,sdim_in
   if (domlo_level(ilev,dir).ne.0) then
    print *,"domlo_level invalid"
    stop
   endif
   do i=domlo_level(ilev,dir)-2*bfactmax, &
        domhi_level(ilev,dir)+2*bfactmax
    inode=2*i
    if ((inode.lt.cache_index_low).or. &
        (inode+1.gt.cache_index_high)) then
     print *,"icell outside of cache range"
     stop
    endif
    nhalf=1
    problo_arr(1)=0.0d0
    problo_arr(2)=0.0d0
    call gridsten1D(xsten_cache,problo_arr,i, &
      domlo_in, &
      bfact_space_order(0),dx_local,dir,nhalf) 
    grid_cache(ilev,inode,dir)=xsten_cache(0)
    grid_cache(ilev,inode+1,dir)=xsten_cache(1)
   enddo ! i
  enddo ! dir=1..sdim_in

  if (ilev.lt.cache_max_level) then
   do dir=1,sdim_in
    dxlevel(ilev+1,dir)=0.5d0*dxlevel(ilev,dir)
    domlo_level(ilev+1,dir)=2*domlo_level(ilev,dir)
    domhi_level(ilev+1,dir)=2*(domhi_level(ilev,dir)+1)-1
   enddo
  else if (ilev.eq.cache_max_level) then
   ! do nothing
  else
   print *,"ilev invalid"
   stop
  endif

 enddo !ilev=0...cache_max_level

 grid_cache_allocated=1

 adv_vel=0.0d0
 gravity=0.0d0
 visc_coef=0.0d0
 num_state_base=2
 num_materials_vel=1
 num_materials_scalar_solve=nmat_in
 do im=1,nmat_in
   fort_drhodt(im)=0.0d0
   fort_drhodz(im)=0.0d0
   fort_tempconst(im)=1.0d0
   fort_initial_temperature(im)=1.0d0
   fort_tempcutoff(im)=1.0E+20
   fort_tempcutoffmax(im)=1.0E+20
   fort_denconst(im)=1.0d0
   fort_density_floor(im)=1.0d0
   fort_density_ceiling(im)=1.0d0
   fort_viscconst(im)=0.0d0
   fort_stiffCP(im)=1.0d0
   fort_material_type(im)=0
   override_density(im)=0
 enddo ! im=1..nmat_in

 do iten=1,local_nten
   use_exact_temperature(iten)=0
   use_exact_temperature(iten+local_nten)=0
   saturation_temp(iten)=0.0d0
   saturation_temp(iten+local_nten)=0.0d0
   latent_heat(iten)=0.0d0
   latent_heat(iten+local_nten)=0.0d0
   reaction_rate(iten)=0.0d0
   reaction_rate(iten+local_nten)=0.0d0
   freezing_model(iten)=0
   freezing_model(iten+local_nten)=0
 enddo

 if (probtype_in.eq.3) then

   saturation_temp(1)=0.0d0
   saturation_temp(2)=273.0d0
   fort_tempconst(1)=273.0
   fort_tempconst(2)=273.0-TDIFF_in
   fort_initial_temperature(1)=fort_tempconst(1)
   fort_initial_temperature(2)=fort_tempconst(2)
   latent_heat(1)=0.0d0 ! material 1 converted to material 2
   latent_heat(2)=-abs(latent_heat_in) ! material 2 converted to material 1
   fort_alpha(1)=1.0d0
   fort_alpha(2)=1.0d0
   fort_stefan_number(1)=TDIFF_in/abs(latent_heat_in)
   fort_stefan_number(2)=TDIFF_in/abs(latent_heat_in)
   fort_jacob_number(1)=fort_stefan_number(1)
   fort_jacob_number(2)=fort_stefan_number(2)

   call find_lambda(lmSt,fort_stefan_number(2))

   print *,"lmSt= ",lmSt

   fort_beta(1)=lmSt
   fort_beta(2)=lmSt
    ! sqrt(alpha time_radblob)*two*lmSt=radblob
   call solidification_front_time(lmSt, &
    fort_alpha(2),fort_time_radblob(2),radblob)
   fort_time_radblob(1)=fort_time_radblob(2)

   call solidification_front_radius_driver(fort_time_radblob(2),test_radblob)
   print *,"radblob= ",radblob
   print *,"test_radblob=",test_radblob

   call solidification_front_speed_driver(fort_time_radblob(2),max_front_vel) 

   fort_time_radblob(1)=fort_time_radblob(2)
   print *,"probtype_in=",probtype_in
   print *,"Stefan_number= ",fort_stefan_number(2)
   print *,"Jacob_number= ",fort_jacob_number(2)
   print *,"lmSt= ",lmSt
   print *,"alpha= ",fort_alpha(2)
   print *,"beta= ",fort_beta(2)
   print *,"time_radblob is the time for the front to grow from"
   print *,"r=0 to r=radblob"
   print *,"time_radblob=",fort_time_radblob(2)
   print *,"radius doubling time=4 * time_radblob - time_radblob=", &
      3.0d0*fort_time_radblob(2)
   print *,"max_front_vel=",max_front_vel
   print *,"front location: 2 beta sqrt(alpha t) "

 else if (probtype_in.eq.4) then

   saturation_temp(1)=273.0d0
   saturation_temp(2)=273.0d0
   fort_tempconst(1)=273.0
   fort_tempconst(2)=273.0
   fort_initial_temperature(1)=fort_tempconst(1)
   fort_initial_temperature(2)=fort_tempconst(2)
   latent_heat(1)=-abs(latent_heat_in) ! material 1 converted to material 2
   latent_heat(2)=0.0d0 ! material 2 converted to material 1
   ireverse=0
   isink=0
   fort_alpha(1)=1.0d0
   fort_alpha(2)=1.0d0
   fort_stefan_number(1)=TDIFF_in/abs(latent_heat_in)
   fort_stefan_number(2)=TDIFF_in/abs(latent_heat_in)
   fort_jacob_number(1)=fort_stefan_number(1)
   fort_jacob_number(2)=fort_stefan_number(2)

     ! variable_temperature_drop.F90
     ! TDIFF=T_DISK_CENTER - TSAT
   call axisymmetric_disk_init(latent_heat(1),saturation_temp(1), &
     TDIFF_in,fort_heatviscconst(1),radblob_in,stefan_flag, &
     ireverse,isink,probtype_in,1)  ! polar_flag=1

   fort_beta(1)=0.0d0
   fort_beta(2)=0.0d0
   fort_time_radblob(1)=0.0d0
   fort_time_radblob(2)=0.0d0

   call disk_get_speed(1,max_front_vel)
   max_front_vel=abs(max_front_vel) 

   print *,"probtype_in=",probtype_in
   print *,"max_front_vel=",max_front_vel

 else if (probtype_in.eq.400) then

    ! max_front_vel
   if ((abs(latent_heat_in).gt.0.0d0).and. &
       (fort_tempconst(1).ge.0.0d0).and. &
       (fort_tempconst(2).ge.0.0d0)) then
    max_front_vel=abs(TDIFF_in)* &
      (fort_tempconst(1)+fort_tempconst(2))/abs(latent_heat_in)
    if (max_front_vel.gt.0.0d0) then
     ! do nothing
    else
     print *,"max_front_vel invalid probtype_in=",probtype_in
     stop
    endif
   else
    print *,"latent_heat_in or fort_tempconst invalid"
    stop
   endif

   saturation_temp(1)=273.0d0
   saturation_temp(2)=273.0d0
   fort_tempconst(1)=273.0
   fort_tempconst(2)=273.0
   fort_initial_temperature(1)=fort_tempconst(1)
   fort_initial_temperature(2)=fort_tempconst(2)
   latent_heat(1)=abs(latent_heat_in) ! material 1 converted to material 2
   latent_heat(2)=0.0d0 ! material 2 converted to material 1
   ireverse=0
   isink=0
   fort_alpha(1)=1.0d0
   fort_alpha(2)=1.0d0
   fort_stefan_number(1)=TDIFF_in/abs(latent_heat_in)
   fort_stefan_number(2)=TDIFF_in/abs(latent_heat_in)
   fort_jacob_number(1)=fort_stefan_number(1)
   fort_jacob_number(2)=fort_stefan_number(2)

   fort_beta(1)=0.0d0
   fort_beta(2)=0.0d0
   fort_time_radblob(1)=0.0d0
   fort_time_radblob(2)=0.0d0

   print *,"probtype_in=",probtype_in
   print *,"max_front_vel=",max_front_vel


 else if (probtype_in.eq.401) then

    ! 1=liquid  2=gas  3=ice
    ! 12 13 23 21 31 32
    !  1  2  3  4  5  6
   saturation_temp(1)=275.0d0
   saturation_temp(2)=273.0d0
   saturation_temp(3)=273.0d0
   saturation_temp(4)=273.0d0
   saturation_temp(5)=273.0d0
   saturation_temp(6)=273.0d0
   fort_tempconst(1)=273.0
   fort_tempconst(2)=273.0
   fort_tempconst(3)=273.0
   fort_initial_temperature(1)=fort_tempconst(1)
   fort_initial_temperature(2)=fort_tempconst(2)
   fort_initial_temperature(3)=fort_tempconst(3)
   latent_heat(1)=1.0D0  ! liquid to gas
   latent_heat(2)=0.0D0
   latent_heat(3)=0.0D0
   latent_heat(4)=0.0D0
   latent_heat(5)=1.0D0  ! ice to liquid
   latent_heat(6)=0.0D0
  
   ireverse=0
   isink=0


   fort_alpha(1)=1.0d0
   fort_alpha(2)=1.0d0
   fort_alpha(3)=1.0d0
   fort_stefan_number(1)=1.0D0
   fort_stefan_number(2)=1.0D0
   fort_stefan_number(3)=1.0D0
   fort_jacob_number(1)=fort_stefan_number(1)
   fort_jacob_number(2)=fort_stefan_number(2)
   fort_jacob_number(3)=fort_stefan_number(3)

   fort_beta(1)=0.0d0
   fort_beta(2)=0.0d0
   fort_beta(3)=0.0d0

   fort_time_radblob(1)=0.0d0
   fort_time_radblob(2)=0.0d0
   fort_time_radblob(3)=0.0d0

    ! max_front_vel
   if (abs(latent_heat(5)).gt.0.0d0) then
     ! FOR YANG: bound on initial speed is 
     ! abs(TDIFF_in)*(k_water + k_ice)/(L_{ice,water} * h)
    max_front_vel=abs(TDIFF_in)* &
      (fort_heatviscconst(1)+fort_heatviscconst(3))/abs(latent_heat(5))
    if (max_front_vel.gt.0.0d0) then
     ! do nothing
    else
     print *,"max_front_vel invalid probtype_in=",probtype_in
     stop
    endif
   else
    print *,"latent_heat_in or fort_heatviscconst invalid"
    stop
   endif
   print *,"probtype_in=",probtype_in
   print *,"max_front_vel=",max_front_vel

 else if (probtype_in.eq.5) then

   saturation_temp(1)=273.0d0
   saturation_temp(2)=273.0d0
   fort_tempconst(1)=273.0d0  ! initial temperature on the left
   fort_tempconst(2)=273.0d0  ! should not be used.
   fort_initial_temperature(1)=fort_tempconst(1)
   fort_initial_temperature(2)=fort_tempconst(2)
   latent_heat(1)=0.0d0 ! material 1 converted to material 2
   latent_heat(2)=-abs(latent_heat_in) ! material 2 converted to material 1
   ireverse=0
   isink=0
   fort_alpha(1)=1.0d0
   fort_alpha(2)=1.0d0
   fort_stefan_number(1)=TDIFF_in/abs(latent_heat_in)
   fort_stefan_number(2)=TDIFF_in/abs(latent_heat_in)
   fort_jacob_number(1)=fort_stefan_number(1)
   fort_jacob_number(2)=fort_stefan_number(2)

   fort_beta(1)=0.0d0
   fort_beta(2)=0.0d0
   fort_time_radblob(1)=0.0d0
   fort_time_radblob(2)=0.0d0

   max_front_vel=1.0d0

   print *,"probtype_in=",probtype_in
   print *,"max_front_vel=",max_front_vel

 else if (probtype_in.eq.0) then
   ! do nothing
 else if (probtype_in.eq.1) then

   if (dirichlet_annulus.eq.1) then

     ! 12,13,23,21,31,32
     ! material 2 is in the middle
    if (local_nten.eq.3) then
     latent_heat(1)=1.0d0  ! material 1 converted to material 2
     use_exact_temperature(1)=2  ! pass material id=2 to exact_temperature
     latent_heat(2)=0.0d0  ! no 1-3 interface
     latent_heat(3)=1.0d0  ! material 2 converted to material 3
     use_exact_temperature(3)=2 ! pass material id=2 to exact_temperature
     latent_heat(4)=0.0d0  ! no 2 -> 1
     latent_heat(5)=0.0d0  ! no 3 -> 1
     latent_heat(6)=0.0d0  ! no 3 -> 2
    else
     print *,"local_nten invalid"
     stop
    endif

   else if (dirichlet_annulus.eq.0) then
    ! do nothing
   else
    print *,"dirichlet_annulus invalid"
    stop
   endif

 else if (probtype_in.eq.2) then
   ! do nothing
 else if(probtype_in.eq.19)then

   if (local_nten.eq.3) then
    ! material 2 is in the middle
    latent_heat(1)=1.0d0  ! material 1 converted to material 2
    use_exact_temperature(1)=2  ! pass material id=2 to exact_temperature
    latent_heat(2)=0.0d0  ! no 1-3 interface
    latent_heat(3)=1.0d0  ! material 2 converted to material 3
    use_exact_temperature(3)=2 ! pass material id=2 to exact_temperature
    latent_heat(4)=0.0d0  ! no 2 -> 1
    latent_heat(5)=0.0d0  ! no 3 -> 1
    latent_heat(6)=0.0d0  ! no 3 -> 2
   else
    print *,"local_nten invalid"
    stop
   endif

 else if(probtype_in.eq.13)then

   if (dirichlet_pentafoil.eq.1) then

     ! material 2 is in the middle
    latent_heat(1)=1.0d0  ! material 1 converted to material 2
    use_exact_temperature(1)=2  ! pass material id=2 to exact_temperature
    latent_heat(2)=0.0d0  ! no 1-3 interface
    latent_heat(3)=1.0d0  ! material 2 converted to material 3
    use_exact_temperature(3)=2 ! pass material id=2 to exact_temperature
    latent_heat(4)=0.0d0  ! no 2 -> 1
    latent_heat(5)=0.0d0  ! no 3 -> 1
    latent_heat(6)=0.0d0  ! no 3 -> 2

   else if (dirichlet_pentafoil.eq.0) then
    ! do nothing
   else
    print *,"dirichlet_pentafoil invalid"
    stop
   endif

 else if(probtype_in.eq.14)then
   ! do nothing
 else if(probtype_in.eq.15)then
   ! do nothing
 else if (probtype_in.eq.16) then ! multiscale

  ! variable_temperature_drop.F90
  ! TDIFF=T_DISK_CENTER - TSAT
  ! L=1.0d0
  ! TSAT=NB_top
  ! TDIFF=NB_bot
  ! RR=probleny
  ! ireverse=0
  ! isink=0
  call axisymmetric_disk_init(1.0d0,NB_top, &
    NB_bot,fort_heatviscconst(1),probleny,stefan_flag, &
    0,0,probtype_in,0)  ! polar_flag=0

 else if(probtype_in.eq.17)then
   ! do nothing
 else if(probtype_in.eq.20)then
   ! do nothing
 else
   print *,"probtype_in invalid2 ",probtype_in
   stop
 endif

 MOFITERMAX=15
 ngeom_raw=1+AMREX_SPACEDIM
 ngeom_recon=3+2*AMREX_SPACEDIM
 ngeom_recon_in=3+2*AMREX_SPACEDIM

 call initmof(order_algorithm, &
          nmat_in,MOFITERMAX, &
          0, &  ! MOF_DEBUG_RECON_in=0
          1, &  ! MOF_TURN_OFF_LS_in=1
          1, &  ! nthreads=1
          nmax) 

! initmof(..., mof_debug_recon_in, mof_turn_off_ls_in,nthreads,nmax)  MOF.F90
!  mof_debug_recon_in = 1 , output a lot ,    
!                     = 0 , nothing
! mof_turn_off_ls_in = 1 ,  not use levelset as input,  use centroid.. ,    
!                    = 0 ,  use levelset as input
! 

 print *,"in main.F90: probtype_in= ",probtype_in
 print *,"BEFORE TIME LOOP, N_CURRENT= ",N_CURRENT
 print *,"BEFORE TIME LOOP, M_CURRENT= ",M_CURRENT
 print *,"TSTOP= ",TSTOP
 print *,"fixed_dt_main= ",fixed_dt_main
 print *,"fixed_dt_current= ",fixed_dt_current
 print *,"these vars declared in vof_cisl.F90 are init. in main.F90"
 print *,"radcen, radeps, thermal_delta, pentaeps declared: vof_cisl.F90"
 print *,"dirichlet_pentafoil declared: vof_cisl.F90"
 print *,"dirichlet_annulus declared: vof_cisl.F90"
 print *,"radial_variation declared: vof_cisl.F90"
 print *,"radcen= ",radcen
 print *,"radeps= ",radeps
 print *,"thermal_delta= ",thermal_delta
 print *,"pentaeps= ",pentaeps
 print *,"dirichlet_pentafoil= ",dirichlet_pentafoil
 print *,"dirichlet_annulus= ",dirichlet_annulus
 print *,"radial_variation= ",radial_variation
 print *,"ERRTOL defined in BICGSTAB_Yang_MULTI.F90 and init in main.F90"
 print *,"ERRTOL= ",ERRTOL
 print *,"VOFTOL_local= ",VOFTOL_local
 print *,"pcurve_num= ",pcurve_num

 if (dirichlet_annulus.eq.1) then
   if ((radial_variation.eq.0).or.(radial_variation.eq.1)) then
    ! do nothing
   else
    print *,"radial_variation invalid"
    stop
   endif
 else if (dirichlet_annulus.eq.0) then
   if (radial_variation.eq.0) then
    ! do nothing
   else
    print *,"radial_variation invalid"
    stop
   endif
 else
   print *,"dirichlet_annulus invalid"
   stop
 endif

 allocate(XLINE(-1:N_CURRENT+1))
 allocate(YLINE(-1:N_CURRENT+1))
 ! init nodes
 do i = -1 , N_CURRENT+1
     XLINE(i)= problo + (i)*h_in
     YLINE(i)= problo + (i)*h_in
 enddo
 allocate(xCC(-1:N_CURRENT))
 allocate(yCC(-1:N_CURRENT))

 do i= -1, N_CURRENT
     xCC(i)= problo + (i+0.5d0)*h_in
     yCC(i)= problo + (i+0.5d0)*h_in
 enddo

 !--intial cell
 allocate(CELL_FAB(-1:N_CURRENT,-1:N_CURRENT))
 do i = -1 , N_CURRENT
    do j= -1 , N_CURRENT
       call init_cell(N_CURRENT,h_in,xCC,yCC,i,j,CELL_FAB(i,j))
    enddo
 enddo

 ! temperature, velocity, interface reconstruction, level set
 local_state_ncomp=nmat_in+local_nten*sdim_in+ &
    ngeom_recon_in*nmat_in+nmat_in*(sdim_in+1)

 allocate(vf(-1:N_CURRENT,-1:N_CURRENT,nmat_in)) 
 allocate(mofdata_FAB_in(-1:N_CURRENT,-1:N_CURRENT,ngeom_recon_in*nmat_in)) 
 allocate(CENTROID_FAB(-1:N_CURRENT,-1:N_CURRENT,nmat_in)) 
 allocate(centroid_mult(-1:N_CURRENT,-1:N_CURRENT,nmat_in,sdim_in)) 
 allocate(T(-1:N_CURRENT,-1:N_CURRENT,local_state_ncomp)) 
 allocate(T_new(-1:N_CURRENT,-1:N_CURRENT,local_state_ncomp)) 

 call convert_lag_to_eul(cache_max_level,sdim_in)

 ! init velocity in: vof_cisl.F90
 do iten=1,local_nten
   scomp=nmat_in+(iten-1)*sdim_in+1
   CALL INIT_V(N_CURRENT,xCC,yCC,probtype_in,iten,scomp,sdim_in,T)
 enddo

  ! STEP 0 FOR YANG: initial time step
 if (fixed_dt_current.eq.0.0) then
    if (max_front_vel.gt.0.0d0) then
     deltat_in = h_in*0.25d0/max_front_vel
    else if (max_front_vel.eq.0.0d0) then
     deltat_in = 0.5d0*h_in
    else
     print *,"max_front_vel invalid"
     stop
    endif
 else if (fixed_dt_current.gt.0.0) then
    deltat_in=fixed_dt_current
    if (abs(TSTOP-fixed_dt_current*M_current).gt.VOFTOL*fixed_dt_current) then
     print *,"TSTOP and fixed_dt_current are inconsistent"
     stop
    endif
 else
    print *,"fixed_dt_current invalid"
    stop
 endif

 print *,"deltat_in=",deltat_in

 if (probtype_in.eq.3) then
    print *,"approx number time steps to double radius: ", &
     NINT(3.0d0*fort_time_radblob(2)/deltat_in)
 else if (probtype_in.eq.4) then
     ! max_front_vel * N * dt \approx radblob
     ! N=2*radblob/(max_front_vel * dt)
    print *,"approx number time steps to double radius: ", &
     NINT(radblob/(max_front_vel*deltat_in))
 else if (probtype_in.eq.400) then
  ! do nothing
 else if (probtype_in.eq.401) then
  ! do nothing
 else if (probtype_in.eq.5) then
  print *,"Velocity is 1"
  print *,"number of steps to move 1 unit: ", &
    NINT(1.0d0/(deltat_in*max_front_vel))
 else if (probtype_in.eq.0) then
    ! do nothing
 else if (probtype_in.eq.1) then
    ! do nothing
 else if (probtype_in.eq.2) then
    ! do nothing
 else if(probtype_in.eq.19)then
    ! do nothing
 else if(probtype_in.eq.13)then
    ! do nothing
 else if(probtype_in.eq.14)then
    ! do nothing
 else if(probtype_in.eq.15)then
    ! do nothing
 else if (probtype_in.eq.16) then
    ! do nothing
 else if(probtype_in.eq.17)then
    ! do nothing
 else if(probtype_in.eq.20)then
    ! do nothing
 else
    print *,"probtype_in invalid3 ",probtype_in
    stop
 endif

 if (M_MAX_TIME_STEP.ge.M_CURRENT) then
  allocate(Ts(M_MAX_TIME_STEP+1))
  do i = 1,M_MAX_TIME_STEP+1
    Ts(i) = (i-1)*deltat_in
  enddo
 else
  print *,"M_MAX_TIME_STEP or M_CURRENT invalid"
  stop
 endif
 
    ! STEP 1 FOR YANG: initialize zeroeth and first order moments
    ! STEP 2 FOR YANG: initialize levelset functions
    ! in: multimat_FVM.F90
    ! init_vfncen calls:
    !  AdaptQuad_2d  (in multimat_FVM.F90)
    !  AdaptQuad_2d calls
    !   dist_concentric (in multimat_FVM.F90)
    ! TYPE(POINTS),DIMENSION(:,:,:),allocatable :: CENTROID_FAB 
    ! The centroid here is in an absolute coordinate system.
 call init_vfncen(N_CURRENT,CELL_FAB,nmat_in,dx_in,CENTROID_FAB,vf,probtype_in)

    ! real(kind=8),dimension(:,:,:,:),allocatable :: centroid_mult
    ! convert_cen is declared in: vfrac_pair.F90
    ! centroid_mult is still in absolute coordinate system.
 call convert_cen(nmat_in,sdim_in,N_CURRENT,CENTROID_FAB,centroid_mult)
    
    ! INIT_MOFdata (in multimat_FVM.F90)
    ! mofdata_FAB_in centroids are relative to a given cell's centroid.
 call init_mofdata(N_CURRENT,sdim_in,dx_in,nmat_in,local_nten,CELL_FAB, &
    vf,CENTROID_FAB,mofdata_FAB_in)
   
 do i = 0,N_CURRENT-1
    do imof=1,ngeom_recon_in*nmat_in
     mofdata_FAB_in(i,-1,imof) = mofdata_FAB_in(i,0,imof)
     mofdata_FAB_in(i,N_CURRENT,imof) = mofdata_FAB_in(i,N_CURRENT-1,imof)
    enddo
 enddo

 do i  = -1,N_CURRENT
    do imof=1,ngeom_recon_in*nmat_in
     mofdata_FAB_in(-1,i,imof) = mofdata_FAB_in(0,i,imof)
     mofdata_FAB_in(N_CURRENT,i,imof) = mofdata_FAB_in(N_CURRENT-1,i,imof)   
    enddo
 enddo

 if (probtype_in.eq.19) then
    ! Np,Mp,r_polar,z_polar,dr_polar,dz_polar,upolar
    ! declared in vof_cisl.F90
  call set_polar_2d(sdim_in,Np,Mp,fort_heatviscconst(2),deltat_in, &
      r_polar,z_polar,dr_polar,dz_polar,upolar,deltat_polar, &
      subcycling_step)
 endif

  ! STEP 3: Initialize Temperature
 do i= -1,N_CURRENT
 do j= -1,N_CURRENT

    do im = 1,nmat_in
     T(i,j,im)=1.0d0
    enddo

    do im=1,nmat_in
     vofcomp=(im-1)*ngeom_recon_in+1
     local_vof=mofdata_FAB_in(i,j,vofcomp)

     if (local_vof.gt.0.0d0) then
      xcen=centroid_mult(i,j,im,1)
      ycen=centroid_mult(i,j,im,2)
     else if (local_vof.eq.0.0d0) then
      xcen=xCC(i)
      ycen=yCC(j)
     else
      print *,"local_vof invalid"
      stop
     endif

     xcen_vec(1)=xcen
     xcen_vec(2)=ycen

     time_init=0.0

     if ((probtype_in.eq.0).or. &
         (probtype_in.eq.2)) then
      T(i,j,im)=2.0
      if (1.eq.0) then
        ! in: multimat_FVM.F90
       T(i,j,im)=exact_temperature(xcen_vec,time_init,im,probtype_in, &
        nmat_in,fort_heatviscconst)
      endif
     else if (probtype_in.eq.1) then
        ! in: multimat_FVM.F90
      T(i,j,im)=exact_temperature(xcen_vec,time_init,im,probtype_in, &
       nmat_in,fort_heatviscconst)
     else if (probtype_in.eq.3) then
      rstefan=sqrt((xcen-xblob)**2+(ycen-yblob)**2)
      if (im.eq.1) then
       T_FIELD=saturation_temp(2)
      else if (im.eq.2) then
       if (rstefan.le.radblob) then
        T_FIELD=saturation_temp(2)
       else if (rstefan.ge.radblob) then
        call liquid_temperature( &
         fort_beta(2), &
         fort_tempconst(2), &
         abs(latent_heat_in), &
         fort_stiffCP(2), &
         fort_stefan_number(2), &
         rstefan, &
         fort_time_radblob(2), &
         fort_heatviscconst(2), &
         T_FIELD)
       else
        print *,"rstefan invalid"
        stop
       endif
      else
       print *,"im invalid"
       stop
      endif
      T(i,j,im)=T_FIELD
      if (1.eq.0) then
       print *,"i,j,r,im,T_FIELD ",i,j,rstefan,im,T_FIELD
      endif
      if (T_FIELD.lt.fort_tempconst(2)-1.0D-7) then
       print *,"bust: fort_time_radblob(2)= ",fort_time_radblob(2)
       stop
      endif

     else if (probtype_in.eq.4) then

      rstefan=sqrt((xcen-xblob)**2+(ycen-yblob)**2)
      if (im.eq.2) then
       T_FIELD=saturation_temp(2)
      else if (im.eq.1) then
       if (rstefan.ge.radblob) then
        T_FIELD=saturation_temp(2)
       else if (rstefan.le.radblob) then
         ! TDIFF=T_DISK_CENTER - TSAT 
        call disk_eval_initial_temp(rstefan,T_FIELD)
       else
        print *,"rstefan invalid"
        stop
       endif
      else
       print *,"im invalid"
       stop
      endif
      T(i,j,im)=T_FIELD
      if (1.eq.0) then
       print *,"i,j,r,im,T_FIELD ",i,j,rstefan,im,T_FIELD
      endif
      if (T_FIELD.lt.fort_tempconst(2)-abs(TDIFF_in)) then
       print *,"bust: fort_tempconst(2)= ",fort_tempconst(2)
       print *,"bust: TDIFF_in= ",TDIFF_in
       stop
      endif

     else if (probtype_in.eq.400) then

      T_FIELD=saturation_temp(1)
      T(i,j,im)=T_FIELD

     else if (probtype_in.eq.401) then

      T_FIELD=273.0d0
      T(i,j,im)=T_FIELD

     else if (probtype_in.eq.5) then

      if (im.eq.1) then
       T_FIELD=saturation_temp(1)
      else if (im.eq.2) then
       ! (xcen,ycen)
       T_FIELD=272.0d0+exp(-(xcen-0.2d0))
      else
       print *,"im invalid"
       stop
      endif
      T(i,j,im)=T_FIELD
      if (1.eq.0) then
       print *,"i,j,r,im,T_FIELD ",i,j,rstefan,im,T_FIELD
      endif

     else if (probtype_in.eq.19) then   ! annulus cvg test

       ! pcenter declared in vof_cisl.F90
      pcenter(1)=0.5d0
      pcenter(2)=0.5d0

      T(i,j,im)=0.0d0

      if (im.eq.2) then   
       if (local_vof.gt.0.0d0) then
        call polar_cart_interpolate(Np,Mp,upolar,pcenter,rlo,rhi, &
                 xcen_vec,T(i,j,im))
       else if (local_vof.eq.0.0d0) then
        ! do nothing
       else
        print *,"local_vof invalid"
        stop
       endif
      else if ((im.eq.1).or.(im.eq.3)) then
       ! do nothing
      else
       print *,"im invalid"
       stop
      endif

     elseif (probtype_in.eq.13)then
       T(i,j,im)=exact_temperature(xcen_vec,time_init,im,probtype_in, &
        nmat_in,fort_heatviscconst)

     elseif(probtype_in.eq.14)then
      T(i,j,im)=2.0

     elseif(probtype_in.eq.15)then

      cc=0.5d0
      if((xcen .eq. 0.5d0).and. &
         (ycen .eq. 0.5d0))then
        T(i,j,im)=1.0d0
      else
        call dist_to_boundary(xcen_vec,dtemp1)
        call l2normd(2,xcen_vec,cc, dtemp2)
        T(i,j,im)= 1.0d0+dtemp2/dtemp1*(10.0d0-1.0d0)
      endif

     else if (probtype_in.eq.16) then

        ! (xcen,ycen) is centroid of material im region.
       call disk_eval_initial_temp(ycen,T_FIELD)
       T(i,j,im)=T_FIELD

     elseif(probtype_in.eq.17)then
      T(i,j,im)=2.0
     elseif(probtype_in.eq.20)then
      cc =0.5d0
      if((xcen .eq. 0.5d0).and. &
         (ycen .eq. 0.5d0))then
        T(i,j,im)=1.0d0
      else
        call dist_to_boundary(xcen_vec,dtemp1)
        call l2normd(2,xcen_vec,cc, dtemp2)
        T(i,j,im)= 1.0d0+dtemp2/dtemp1*(10.0d0-1.0d0)
      endif

     else
      print *,"probtype_in invalid4 ",probtype_in
      stop
     endif
    enddo ! im=1..nmat_in

 enddo
 enddo

 do i= -1,N_CURRENT
 do j= -1,N_CURRENT
     do im = 1,nmat_in

      if(vf(i,j,im).le.VOFTOL) then
       sumT  = 0.0d0
       sumvf = 0.0d0
       do im1 = 1,nmat_in
        sumT  = sumT  + vf(i,j,im1)*T(i,j,im1) 
        sumvf = sumvf + vf(i,j,im1)            
       enddo
       T(i,j,im) = sumT/sumvf
      endif

     enddo ! im
 enddo
 enddo

 do i= -1,N_CURRENT
 do j= -1,N_CURRENT
     do imof=1,ngeom_recon_in*nmat_in
      scomp=nmat_in+local_nten*sdim_in+imof
      T(i,j,scomp)=mofdata_FAB_in(i,j,imof)
     enddo
 enddo
 enddo

    ! INIT_LS declared in vof_cisl.F90
 scomp=nmat_in+local_nten*sdim_in+ngeom_recon*nmat_in+1
 CALL INIT_LS(N_CURRENT,xCC,yCC,dx_in, &
     probtype_in,nmat_in,scomp,sdim_in,T,local_state_ncomp)

 flxavg1=0.0d0
 flxavg2=0.0d0

 print *,"nmat_in=",nmat_in
 do im=1,nmat_in
    print *,"im,fort_heatviscconst(im) ",im,fort_heatviscconst(im)
 enddo

 iter_average=0.0d0

! STEP 4: TIME LOOP
! BEGIN TIME LOOP - ABOVE INITIALIZATION
!                   BELOW INTEGRATION IN TIME

 tm=1
 finished_flag=0
 do while (finished_flag.eq.0)

    current_time_in=Ts(tm) ! t^{n} (Ts(i)=(i-1) * deltat)

    print *,"STEP (>=1), TIME, DT ",tm,current_time_in,deltat_in

    ngeom_recon_in=2*sdim_in+3
    nx_in=N_CURRENT
    ny_in=N_CURRENT
    bicgstab_tol_in=1.0D-10
    precond_type_in=1 ! 0 M=I  1=Jacobi precond.
    hflag=0
    do im=1,nmat_in
     alpha_in(im)=fort_heatviscconst(im)
    enddo
    lox_in=0
    loy_in=0
    hix_in=nx_in-1
    hiy_in=ny_in-1
    allocate(UNEW_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1,local_state_ncomp)) 
    allocate(UOLD_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1,local_state_ncomp)) 
    allocate(beta_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1,nmat_in)) 
    allocate(VFRAC_MOF_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1, &
       nmat_in)) 

    do i=lox_in-1,hix_in+1
    do j=loy_in-1,hiy_in+1
     do im=1,nmat_in
      xgrid=(i+0.5)*h_in
      ygrid=(j+0.5)*h_in

      sumvf=0.0
      sum_alpha=0.0
      do im1=1,nmat_in
       vofcomp=ngeom_recon_in*(im1-1)+1
       sumvf=sumvf+mofdata_FAB_in(i,j,vofcomp)
       sum_alpha=sum_alpha+mofdata_FAB_in(i,j,vofcomp)/ &
               (alpha_in(im1)+1.0E-10)
      enddo
      sum_alpha=sum_alpha/sumvf
      sum_alpha=1.0/sum_alpha 
      beta_in(i,j,im)=sum_alpha

      vofcomp=ngeom_recon_in*(im-1)+1
      VFRAC_MOF_in(i,j,im)=mofdata_FAB_in(i,j,vofcomp)
     enddo ! im=1..nmat_in
     do im=1,local_state_ncomp
      UNEW_in(i,j,im)=T(i,j,im)
      UOLD_in(i,j,im)=T(i,j,im)
     enddo
    enddo
    enddo 

     ! in: BICGSTAB_Yang_MULTI.F90
    call INIT_GLOBALS( &
     local_state_ncomp, &
     local_operator_internal, &
     local_operator_external, &
     local_linear_exact, &
     probtype_in, &
     sdim_in,ngeom_recon_in, &
     nx_in,ny_in,lox_in,loy_in,hix_in,hiy_in, &
     UNEW_in,UOLD_in, &
     beta_in,h_in,precond_type_in,bicgstab_tol_in, &
     VFRAC_MOF_in,nmat_in,alpha_in,deltat_in, &
     mofdata_FAB_in,current_time_in)

    time_n=current_time_in
    time_np1=current_time_in+deltat_in

    nsteps=tm-1
    if (tm.eq.1) then
        ! in: BICGSTAB_Yang_MULTI.F90
     if (fixed_dt_main.eq.0.0d0) then
      total_nsteps_parm=M_MAX_TIME_STEP
     else
      total_nsteps_parm=M_CURRENT
     endif
     call output_solution(UNEW_in,time_n,nsteps,plot_int, &
             total_nsteps_parm, &
             fixed_dt_main)
    endif

     ! Dirichlet BC use t^n+1 data
     ! polar solver called before bicgstab is called.
    if (probtype_in.eq.19) then
     print *,"deltat_polar",deltat_polar                            
     print *,"subcycling_step", subcycling_step    
     print *,"Np,Mp= ",Np,Mp
     do i=1,subcycling_step                               
      call polar_2d_heat(sdim_in,Np,Mp,fort_heatviscconst(2),deltat_polar, &
       r_polar,z_polar,dr_polar,dz_polar,upolar)
      if (i.eq.(i/1000)*1000) then
       print *,"subcycling_step number: i=",i
      endif
     enddo
    endif

    call bicgstab(UNEW_in,hflag,iter)

    iter_average=iter_average+iter

     ! hflag==0
    call set_boundary(UNEW_in,0,local_state_ncomp)
    do i=lox_in-1,hix_in+1
    do j=loy_in-1,hiy_in+1
     do im=1,local_state_ncomp
      UNEW(i,j,im)=UNEW_in(i,j,im)
      UOLD_in(i,j,im)=UNEW_in(i,j,im)
      UOLD(i,j,im)=UOLD_in(i,j,im)
     enddo
    enddo
    enddo

    if (probtype_in.eq.0) then
     ! do nothing
    else if (probtype_in.eq.1) then
     ! do nothing
    else if (probtype_in.eq.2) then
     ! do nothing
    else if (probtype_in.eq.3) then
     ! do nothing
    else if (probtype_in.eq.4) then
     call axisymmetric_disk_advance(deltat_in)
    else if (probtype_in.eq.400) then
     ! do nothing
    else if (probtype_in.eq.401) then
     ! do nothing
    else if (probtype_in.eq.5) then
     ! do nothing
    else if (probtype_in.eq.19) then   ! annulus cvg test
     ! do nothing (polar solver called before bicgstab is called)
    elseif (probtype_in.eq.13)then
     ! do nothing
    elseif(probtype_in.eq.14)then
     ! do nothing
    elseif(probtype_in.eq.15)then
     ! do nothing
    else if (probtype_in.eq.16) then
     call axisymmetric_disk_advance(deltat_in)
    elseif(probtype_in.eq.17)then
     ! do nothing
    elseif(probtype_in.eq.20)then
     ! do nothing
    else
     print *,"probtype_in invalid"
     stop
    endif 

     ! interface updated here: 
     ! input: UOLD
     ! output: UNEW
    call update_interface(UOLD,UNEW,N_CURRENT,local_state_ncomp, &
      dx_in,time_n,deltat_in,nsteps,local_nten,stefan_flag)

    do i= -1,N_CURRENT
    do j= -1,N_CURRENT
     do im=1,nmat_in
      xcen=xCC(i)
      ycen=yCC(j)
      if (probtype_in.eq.0) then
       ! do nothing
      else if (probtype_in.eq.2) then
       ! do nothing
      else if (probtype_in.eq.1) then
       ! do nothing
      else if (probtype_in.eq.3) then
       rstefan=sqrt((xcen-xblob)**2+(ycen-yblob)**2)
       if (rstefan.ge.0.5d0-radblob2) then
        stefan_time=fort_time_radblob(2)+Ts(tm+1)
        call liquid_temperature_driver( &
         rstefan, &
         stefan_time, &
         T_FIELD)
        UNEW(i,j,im)=T_FIELD
       endif
      else if (probtype_in.eq.4) then
       ! do nothing - this is a shrinking disk, outside temperature
       ! is uniform, grad T dot n=0 on the outer walls.
      else if (probtype_in.eq.400) then
       ! do nothing
      else if (probtype_in.eq.401) then
       ! do nothing
      else if (probtype_in.eq.5) then
       if (xcen.ge.1.0-2.0d0*h_in) then
        T_FIELD=272.0d0+exp(-(xcen-0.2d0-Ts(tm+1)))
        UNEW(i,j,im)=T_FIELD
       endif 
      else if (probtype_in.eq.19) then   ! annulus cvg test
       ! do nothing
      elseif (probtype_in.eq.13)then
       ! do nothing
      elseif(probtype_in.eq.14)then
       ! do nothing
      elseif(probtype_in.eq.15)then
       ! do nothing
      elseif(probtype_in.eq.16)then
       ! do nothing
      elseif(probtype_in.eq.17)then
       ! do nothing
      elseif(probtype_in.eq.20)then
       ! do nothing
      else
       print *,"probtype_in invalid5 ",probtype_in
       stop
      endif
     enddo !im=1..nmat_in
    enddo ! j
    enddo ! i

     ! hflag==0
    call set_boundary(UNEW,0,local_state_ncomp)

    do i=lox_in-1,hix_in+1
    do j=loy_in-1,hiy_in+1
     do im=1,local_state_ncomp
      UNEW_in(i,j,im)=UNEW(i,j,im)
      UOLD_in(i,j,im)=UNEW(i,j,im)
      UOLD(i,j,im)=UOLD_in(i,j,im)
     enddo
    enddo
    enddo

    do i= -1,N_CURRENT
     do j= -1,N_CURRENT
      do imof=1,ngeom_recon_in*nmat_in
       scomp=nmat_in+local_nten*sdim_in+imof
       mofdata_FAB_in(i,j,imof)=UNEW_in(i,j,scomp)
       mofdata_FAB(i,j,imof)=UNEW_in(i,j,scomp)
      enddo
      do im=1,nmat_in
       vofcomp=ngeom_recon_in*(im-1)+1
       VFRAC_MOF_in(i,j,im)=mofdata_FAB_in(i,j,vofcomp)
       VFRAC_MOF(i,j,im)=mofdata_FAB_in(i,j,vofcomp)
      enddo
     enddo
    enddo

    nsteps=tm
        ! in: BICGSTAB_Yang_MULTI.F90
    if (fixed_dt_main.eq.0.0d0) then
     if (Ts(tm+1).ge.TSTOP-1.0D-14) then
      total_nsteps_parm=nsteps
     else
      total_nsteps_parm=M_MAX_TIME_STEP
     endif
    else
     total_nsteps_parm=M_CURRENT
    endif
    call output_solution(UNEW_in,time_np1,nsteps,plot_int, &
            total_nsteps_parm, &
            fixed_dt_main)

    call DEALLOCATE_GLOBALS()
   
    deallocate(UOLD_in) 
    deallocate(beta_in) 
    deallocate(VFRAC_MOF_in) 

    do i=lox_in-1,hix_in+1
    do j=loy_in-1,hiy_in+1
     do im=1,local_state_ncomp
      T_new(i,j,im)=UNEW_in(i,j,im)
     enddo
    enddo
    enddo
    icen=(hix_in+lox_in)/2
    jcen=(hiy_in+loy_in)/2
    if (isink.eq.1) then
     do im=1,nmat
      T_new(icen,jcen,im)=TDIFF_in+saturation_temp(1)
     enddo
    else if (isink.eq.0) then
     ! do nothing
    else
     print *,"isink invalid"
     stop
    endif

    deallocate(UNEW_in) 

   T = T_new

   do im=1,nmat_in
    voltotal=0.0d0
    do i= 0,N_CURRENT-1
     do j= 0,N_CURRENT-1
      vofcomp=nmat_in+local_nten*sdim_in+(im-1)*ngeom_recon_in+1
      sumvf=T(i,j,vofcomp)
      voltotal=voltotal+sumvf*h_in*h_in
     enddo
    enddo
    print *,"TIME= ",Ts(tm+1)," MAT= ",im," VOLUME= ",voltotal
    eff_radius=sqrt(voltotal/local_Pi)   ! pi r^2 = V  r=(V/pi)^(1/2)
    if (probtype_in.eq.5) then
     if (im.eq.1) then
      eff_radius=voltotal
     else if (im.eq.2) then
      eff_radius=1.0d0-voltotal
     else
      print *,"im invalid"
      stop
     endif
    endif
    print *,"TIME= ",Ts(tm+1)," MAT= ",im," EFF RADIUS= ",eff_radius
   enddo ! im=1..nmat_in

   if (probtype_in.eq.0) then
    ! do nothing
   else if (probtype_in.eq.2) then
    ! do nothing
   else if (probtype_in.eq.1) then
    ! do nothing
   else if (probtype_in.eq.3) then
    im=1
    stefan_time=fort_time_radblob(2)+Ts(tm+1)
    call solidification_front_radius_driver(stefan_time,expect_radius) 
    print *,"TIME= ",Ts(tm+1)," MAT= ",im," EXACT RADIUS= ",expect_radius
   else if (probtype_in.eq.4) then
    im=1
    expect_radius=axisymmetric_disk_radblob(2)
    print *,"TIME= ",Ts(tm+1)," MAT= ",im," EXACT RADIUS= ",expect_radius
   else if (probtype_in.eq.400) then
    ! do nothing
   else if (probtype_in.eq.401) then
    ! do nothing
   else if (probtype_in.eq.5) then
    expect_radius=0.2d0+Ts(tm+1)
    print *,"TIME= ",Ts(tm+1)," MAT= ",im," EXACT RADIUS= ",expect_radius
   else if (probtype_in.eq.19) then   ! annulus cvg test
    ! do nothing
   elseif (probtype_in.eq.13)then
    ! do nothing
   elseif(probtype_in.eq.14)then
    ! do nothing
   elseif(probtype_in.eq.15)then
    ! do nothing
   elseif(probtype_in.eq.16)then

    flxtot1=0.0d0 
    flxtot2=0.0d0 

    if ((N_CURRENT.eq.16).or. &
        (N_CURRENT.eq.32).or. &
        (N_CURRENT.eq.64).or. &
        (N_CURRENT.eq.128).or. &
        (N_CURRENT.eq.256).or. &
        (N_CURRENT.eq.512).or. &
        (N_CURRENT.eq.1024)) then

      ! filament_test_type declared in vof_cisl.F90
     if (filament_test_type.eq.0) then ! irregular material 2
      y_fluxtest1=0.40625d0
     else if (filament_test_type.eq.1) then ! circular material 2
      y_fluxtest1=0.59375d0
     else
      print *,"filament_test_type invalid"
      stop
     endif

     j_fluxtest=NINT(y_fluxtest1/dx_in(2))
     xlo_fluxtest=0.375d0
     xhi_fluxtest=0.625d0
     ilo_fluxtest=NINT(xlo_fluxtest/dx_in(1))
     ihi_fluxtest=NINT(xhi_fluxtest/dx_in(1))

     isum=0
     do i=ilo_fluxtest,ihi_fluxtest-1
      flxtot1=flxtot1+(T(i,j_fluxtest,2)-T(i,j_fluxtest-1,2))/dx_in(2)
      isum=isum+1
     enddo
     flxtot1=flxtot1/real(isum,8)
     flxavg1=flxavg1+flxtot1
     print *,"TIME,flxtot1 ",Ts(tm+1),flxtot1

     y_fluxtest2=0.84375d0
     j_fluxtest=NINT(y_fluxtest2/dx_in(2))
     xlo_fluxtest=0.0d0
     xhi_fluxtest=1.0d0
     ilo_fluxtest=NINT(xlo_fluxtest/dx_in(1))
     ihi_fluxtest=NINT(xhi_fluxtest/dx_in(1))

     isum=0
     do i=ilo_fluxtest,ihi_fluxtest-1
      flxtot2=flxtot2+(T(i,j_fluxtest,3)-T(i,j_fluxtest-1,3))/dx_in(2)
      isum=isum+1
     enddo
     flxtot2=flxtot2/real(isum,8)
     flxavg2=flxavg2+flxtot2
     print *,"TIME,flxtot2 ",Ts(tm+1),flxtot2

    else
     print *,"N_CURRENT out of range N_CURRENT=",N_CURRENT
     stop
    endif
   
   elseif(probtype_in.eq.17)then
    ! do nothing
   elseif(probtype_in.eq.20)then
    ! do nothing
   else
    print *,"probtype_in invalid20 ",probtype_in
    stop
   endif

   tm=tm+1

   finished_flag=0
   if (Ts(tm).ge.TSTOP-1.0D-14) then
    finished_flag=1
   endif
   if (tm-1.ge.M_MAX_TIME_STEP) then
    finished_flag=1
   endif
   if (tm-1.ge.M_CURRENT) then
    if ((fixed_dt_main.eq.-1.0d0).or. &
        (fixed_dt_main.gt.0.0d0)) then
     if ((finished_flag.ne.1).or. &
         (abs(Ts(tm)-TSTOP).gt.1.0D-14)) then
      print *,"expecting Ts(M_CURRENT+1) == TSTOP"
      stop
     endif
    else if (fixed_dt_main.eq.0.0d0) then
     ! check nothing
    else
     print *,"fixed_dt_main invalid"
     stop
    endif
   else if (tm.ge.2) then
    !check nothing
   else
    print *,"tm invalid"
    stop
   endif

   if ((probtype_in.eq.400).or. &
       (probtype_in.eq.401)) then ! gingerbread man or ice melt

    max_front_vel=0.0
    do i= 0,N_CURRENT-1
    do j= 0,N_CURRENT-1
     do im=1,nmat_in
      vofcomp=nmat_in+local_nten*sdim_in+(im-1)*ngeom_recon_in+1
      sumvf=T(i,j,vofcomp)
      if ((sumvf.ge.VOFTOL_REDIST).and. &
          (sumvf.le.1.0d0-VOFTOL_REDIST)) then
       do im_opp=im+1,nmat_in
        vofcomp2=nmat_in+local_nten*sdim_in+(im_opp-1)*ngeom_recon_in+1
        sumvf2=T(i,j,vofcomp2)
        if ((sumvf2.ge.VOFTOL_REDIST).and. &
            (sumvf2.le.1.0d0-VOFTOL_REDIST)) then
         call get_iten(im,im_opp,iten,nmat_in)
         do ireverse=0,1
          LL=abs(latent_heat(iten+ireverse*local_nten))
          TSAT=saturation_temp(iten+ireverse*local_nten)
          if (LL.gt.0.0d0) then
           test_vel=fort_heatviscconst(im)*abs(T(i,j,im)-TSAT)/LL+ &
                    fort_heatviscconst(im_opp)*abs(T(i,j,im_opp)-TSAT)/LL
           test_vel=test_vel*4.0d0/dx_in(1)
           if (test_vel.gt.max_front_vel) then
            max_front_vel=test_vel
           endif
          else if (LL.eq.0.0d0) then
           ! do nothing
          else
           print *,"LL invalid"
           stop
          endif
         enddo ! ireverse=0..1
        else if ((sumvf2.ge.-VOFTOL_REDIST).and. &
                 (sumvf2.le.VOFTOL_REDIST)) then
         ! do nothing
        else if ((sumvf2.ge.1.0d0-VOFTOL_REDIST).and. &
                 (sumvf2.le.1.0d0+VOFTOL_REDIST)) then
          ! do nothing
        else
         print *,"sumvf2 invalid"
         stop
        endif
       enddo ! im_opp=1..nmat_in
      else if ((sumvf.ge.-VOFTOL_REDIST).and. &
               (sumvf.le.VOFTOL_REDIST)) then
       ! do nothing
      else if ((sumvf.ge.1.0d0-VOFTOL_REDIST).and. &
               (sumvf.le.1.0d0+VOFTOL_REDIST)) then
       ! do nothing
      else
       print *,"sumvf invalid"
       stop
      endif
     enddo ! im=1..nmat_in
    enddo
    enddo

    if (max_front_vel.gt.0.0d0) then
     deltat_in=h_in*0.25d0/max_front_vel
     if (finished_flag.eq.0) then
      Ts(tm+1)=Ts(tm)+deltat_in
     else if (finished_flag.eq.1) then
      ! do nothing
     else
      print *,"finished_flag invalid"
      stop
     endif
    else
     print *,"max_front_vel invalid"
     stop
    endif

   else if (probtype_in.ne.400) then
    ! do not alter dt
   else
    print *,"probtype_in invalid"
    stop
   endif

 enddo ! do while (finished_flag.eq.0)

 iter_average=iter_average/real(M_CURRENT,8)

 print *,"N_CURRENT,M_CURRENT, iter_average=",iter_average

 flxavg1=flxavg1/real(M_CURRENT,8)
 flxavg2=flxavg2/real(M_CURRENT,8)

 if (probtype_in.eq.16) then
    ! time should be 0.0125
   print *,"M_CURRENT,time,dx,y_fluxtest1,flxavg1 ", &
           M_CURRENT,Ts(M_CURRENT+1),dx_in(1),y_fluxtest1,flxavg1
   print *,"M_CURRENT,time,dx,y_fluxtest2,flxavg2 ", &
           M_CURRENT,Ts(M_CURRENT+1),dx_in(1),y_fluxtest2,flxavg2
    ! 32,1  -1.74E-2  VPerr=0.112
    ! 64,2  -6.28E-2  VPerr=0.067
    ! 128,4 -12.6E-2  VPerr=0.009
 endif

 print *,"PROCESSING FOR RELATIVE ERRORS N_CURRENT, im_measure= ", &
   N_CURRENT,im_measure

 allocate(fine_data(lox_in-1:hix_in+1,loy_in-1:hiy_in+1,2))
 do i=lox_in-1,hix_in+1
 do j=loy_in-1,hiy_in+1
  scomp=im_measure
  fine_data(i,j,1)=T(i,j,scomp)
  scomp=nmat_in+local_nten*sdim_in+(im_measure-1)*ngeom_recon+1
  fine_data(i,j,2)=T(i,j,scomp)
 enddo
 enddo

 if (N_CURRENT.eq.N_START) then
  ! do nothing
 else if (N_CURRENT.gt.N_START) then

  N_COARSE=N_CURRENT/2
  print *,"PROCESSING FOR RELATIVE ERRORS N_COARSE, im_measure= ", &
   N_COARSE,im_measure

  err1T=0.0d0
  err2T=0.0d0
  err3T=0.0d0

  err1T_gradient=0.0d0
  err2T_gradient=0.0d0
  err3T_gradient=0.0d0

  local_interp=0.0d0

  count1=0
  count1_gradient=0

  do i=0,N_COARSE-1
  do j=0,N_COARSE-1

   cni=i*2
   cnj=j*2

   if (coarse_data(i,j,2).gt.0.99999d0) then
   
    vfcheck=0

    if ((fine_data(cni,cnj+1,2).lt.0.99999d0).or. &
        (fine_data(cni,cnj,2).lt.0.99999d0).or. &
        (fine_data(cni+1,cnj,2).lt.0.99999d0).or. &
        (fine_data(cni+1,cnj+1,2).lt.0.99999d0)) then
     vfcheck=1
    endif    

    gradient_check=0
    if ((i.gt.0).and.(i.lt.N_COARSE-1).and. &
        (j.gt.0).and.(j.lt.N_COARSE-1)) then

     if ((coarse_data(i,j-1,2).lt.0.99999d0).or. &
         (coarse_data(i,j+1,2).lt.0.99999d0).or. &
         (coarse_data(i+1,j,2).lt.0.99999d0).or. &
         (coarse_data(i-1,j,2).lt.0.99999d0)) then
      gradient_check=1
     endif    
    else
     gradient_check=1
    endif

    if (vfcheck.eq.0) then
     local_interp= &
           (fine_data(cni,cnj+1,1)+ &
            fine_data(cni,cnj,1)+ &
            fine_data(cni+1,cnj,1)+ &
            fine_data(cni+1,cnj+1,1))/4.0d0
     err2T= err2T + (coarse_data(i,j,1)-local_interp)**2
     err1T= err1T + abs(coarse_data(i,j,1)-local_interp)
     if(abs(coarse_data(i,j,1)-local_interp).gt.err3T) then
      err3T=abs(coarse_data(i,j,1)-local_interp)
     endif
     count1=count1+1

     if (gradient_check.eq.0) then
      fine_gradientx= &
           (fine_data(cni+1,cnj+1,1)- &
            fine_data(cni,cnj+1,1)+ &
            fine_data(cni+1,cnj,1)- &
            fine_data(cni,cnj,1))/(2.0d0*dx_in(1))
      fine_gradienty= &
           (fine_data(cni+1,cnj+1,1)- &
            fine_data(cni+1,cnj,1)+ &
            fine_data(cni,cnj+1,1)- &
            fine_data(cni,cnj,1))/(2.0d0*dx_in(2))

      coarse_gradientx=(coarse_data(i+1,j,1)- &
                        coarse_data(i-1,j,1))/(4.0d0*dx_in(1))
      coarse_gradienty=(coarse_data(i,j+1,1)- &
                        coarse_data(i,j-1,1))/(4.0d0*dx_in(2))

      gradient_err=sqrt((coarse_gradientx-fine_gradientx)**2+ &
                        (coarse_gradienty-fine_gradienty)**2)

      err2T_gradient= err2T_gradient + gradient_err**2
      err1T_gradient= err1T_gradient + gradient_err
      if(gradient_err.gt.err3T_gradient) then
       err3T_gradient=gradient_err
      endif
      count1_gradient=count1_gradient+1
     else if (gradient_check.eq.1) then
      ! do nothing
     else
      print *,"gradient_check invalid"
      stop
     endif

    else if (vfcheck.eq.1) then
     ! do nothing
    else
     print *,"vfcheck invalid"
     stop
    endif

   else if ((coarse_data(i,j,2).le.0.99999d0).and. &
            (coarse_data(i,j,2).ge.-VOFTOL)) then
    ! do nothing
   else
    print *,"coarse_data bust"
    stop
   endif

  enddo
  enddo

  dx_coarse=dx_in(1)*2.0d0

  print *,"RELATIVE ERRORS: count1,count1_gradient,dx_coarse ", &
      count1,count1_gradient,dx_coarse
   
  if ((count1.gt.0).and.(count1_gradient.gt.0)) then
   
!   err1T=err1T*(dx_coarse**2)
   err1T=err1T/count1
!   err2T=sqrt(err2T*(dx_coarse**2))
   err2T=sqrt(err2T/count1)
   print *,"RELATIVE ERRORS L1,L2,LINF ",err1T,err2T,err3T

!   err1T_gradient=err1T_gradient*(dx_coarse**2)
   err1T_gradient=err1T_gradient/count1_gradient
!   err2T_gradient=sqrt(err2T_gradient*(dx_coarse**2))
   err2T_gradient=sqrt(err2T_gradient/count1_gradient)
   print *,"RELATIVE ERRORS GRAD L1,L2,LINF ", &
      err1T_gradient,err2T_gradient,err3T_gradient
  else if ((count1.eq.0).or.(count1_gradient.eq.0)) then
   print *,"count1 or count1_gradient is 0"
  else
   print *,"count1 or count1_gradient invalid"
   stop
  endif

  deallocate(coarse_data)
 else
  print *,"N_CURRENT bust"
  stop
 endif
  
 allocate(coarse_data(lox_in-1:hix_in+1,loy_in-1:hiy_in+1,2))
 do i=lox_in-1,hix_in+1
 do j=loy_in-1,hiy_in+1
  coarse_data(i,j,1)=fine_data(i,j,1) 
  coarse_data(i,j,2)=fine_data(i,j,2) 
 enddo
 enddo
 deallocate(fine_data)

 deallocate(vf)
 deallocate(mofdata_FAB_in)
 deallocate(CENTROID_FAB)
 deallocate(centroid_mult)
 deallocate(T)
 deallocate(T_new)

 deallocate(Ts)

 deallocate(CELL_FAB)

 deallocate(xCC)
 deallocate(yCC)

 deallocate(XLINE)
 deallocate(YLINE)

 call delete_mof()

 call deallocate_FSI()

 deallocate(grid_cache)

 call delete_cache()

 deallocate(dxlevel)
 deallocate(domlo_level)
 deallocate(domhi_level)


 if ((probtype_in.eq.4).or. &
     (probtype_in.eq.16)) then
  call axisymmetric_disk_close()
 else
  ! do nothing
 endif

 print *,"AFTER TIME LOOP, N_CURRENT= ",N_CURRENT
 print *,"AFTER TIME LOOP, M_CURRENT= ",M_CURRENT

 N_CURRENT=N_CURRENT*2
 M_CURRENT=M_CURRENT*M_FACTOR

ENDDO ! N_CURRENT.le.N_FINISH

END PROGRAM
