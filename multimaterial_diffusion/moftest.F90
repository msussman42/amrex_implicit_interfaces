
PROGRAM mof_test 
USE probcommon_module
use MOF_routines_module
use geometry_intersect_module
IMPLICIT NONE

real(kind=8),parameter     :: lo= 0.0d0, hi=100.0d0
INTEGER                    :: i,j
real(kind=8)               :: h

real(kind=8)               :: CENTROID_D(2),CENTROID_L(2)

! MOF value transfer
real(kind=8)                          :: xtetlist(3,2,500)
real(kind=8)                          :: mofdata(14)
real(kind=8)                          :: dx(2)
real(kind=8)                          :: multi_centroidA(2)
real(kind=8)                          :: cell_center(2)
real(kind=8)                          :: vf



INTEGER nmat
INTEGER order_algorithm(2)
INTEGER MOFITERMAX
INTEGER imaterial_override

nmat=2
order_algorithm(1)=0
order_algorithm(2)=0
MOFITERMAX=10
imaterial_override=0
ngeom_raw=1+BL_SPACEDIM
ngeom_recon=3+2*BL_SPACEDIM

call initmof(order_algorithm,nmat,MOFITERMAX,imaterial_override)



vf = 0.5d0

h = hi-lo
 cell_center(1) = (lo+hi)/2.0d0
 cell_center(2) = (lo+hi)/2.0d0

 centroid_d(1) = 1.0d0/3.0d0
 centroid_d(2) = 1.0d0/3.0d0
 centroid_l(1) = 1.0d0 - 1.0d0/3.0d0
 centroid_l(2) = 1.0d0 - 1.0d0/3.0d0


xtetlist(:,:,:) = 0.0d0 
dx(1) = h
dx(2) = h


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      mofdata(1) = vf
      mofdata(2) = centroid_d(1)
      mofdata(3) = centroid_d(2)
      mofdata(4) = 0
      mofdata(5) = 0.0d0
      mofdata(6) = 0.0d0
      mofdata(7) = 0.0d0
      mofdata(8) = 1-vf
      mofdata(9) = centroid_l(1)
      mofdata(10) = centroid_l(2)
      mofdata(11) = 0
      mofdata(12) = 0.0d0
      mofdata(13) = 0.0d0
      mofdata(14) = 0.0d0

      call multimaterial_MOF(xtetlist, &
                              500, &
                              cell_center, & 
                              dx, &
                              mofdata,&
                              multi_centroidA, &  
                              0.0d0, &
                              0, &
                              cell_center,&
                              dx, &
                              0, &
                              2, &
                              2)
 
 write(*,*) "mx",mofdata(5),"my",mofdata(6),"inter",mofdata(7)
 write(*,*)  "centroid_dark", mofdata(2),mofdata(3)
 write(*,*)  "centroid_light", mofdata(9),mofdata(10) 



end program
