! in .vmrc file in home dir: set tabstop=1, set shiftwidth=1, set expandtab
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_ArrayLim.H"

      module dsdt_module
      use openacc

      contains

      subroutine fort_dsdtfab( &
       dsdtfab,DIMS(dsdtfab), &
       sdatafab,DIMS(sdatafab), &
       dx, &
       dt, &
       c, &
       lo,hi) &
      bind(c,name='fort_dsdtfab')

      IMPLICIT NONE

      integer, INTENT(in) :: DIMDEC(dsdtfab)
      integer, INTENT(in) :: DIMDEC(sdatafab)
      REAL_T, INTENT(out) :: dsdtfab(DIMV(dsdtfab),AMREX_SPACEDIM+1) 
      REAL_T, INTENT(in) :: sdatafab(DIMV(sdatafab),AMREX_SPACEDIM+1) 
      integer, INTENT(in) :: lo(AMREX_SPACEDIM),hi(AMREX_SPACEDIM)
      REAL_T, INTENT(in) :: dx(AMREX_SPACEDIM)
      REAL_T, INTENT(in) :: c
      REAL_T, INTENT(in) :: dt

      integer i,j,k
      REAL_T c2,h,h2,local_dsdt

      c2=c*c
      h=dx(1)
      h2=dx(1)*dx(1)
 
      print *,"size of real: ",SIZEOF(c2)

#if (AMREX_SPACEDIM==3)

!DISABLE!$acc data copyin(sdatafab,c2,h,h2,dt) copyout(dsdtfab) 
!DISABLE!$acc kernels
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)

       local_dsdt= &
         -c2*(sdatafab(i+1,j,k,2)-sdatafab(i-1,j,k,2))/(two*h)+ &
           half*dt*(sdatafab(i+1,j,k,1)- &
           two*sdatafab(i,j,k,1)+ &
           sdatafab(i-1,j,k,1))/h2
       local_dsdt=local_dsdt- &
         c2*(sdatafab(i,j+1,k,3)-sdatafab(i,j-1,k,3))/(two*h)+ &
          half*dt*(sdatafab(i,j+1,k,1)- &
         two*sdatafab(i,j,k,1)+ &
         sdatafab(i,j-1,k,1))/h2
       local_dsdt=local_dsdt- &
         c2*(sdatafab(i,j,k+1,4)-sdatafab(i,j,k-1,4))/(two*h)+ &
           half*dt*(sdatafab(i,j,k+1,1)- &
           two*sdatafab(i,j,k,1)+ &
        sdatafab(i,j,k-1,1))/h2
       dsdtfab(i,j,k,1)=local_dsdt

       dsdtfab(i,j,k,2)= &
        -(sdatafab(i+1,j,k,1)-sdatafab(i-1,j,k,1))/(two*h)+ &
         half*dt*(sdatafab(i+1,j,k,2)- &
        two*sdatafab(i,j,k,2)+ &
        sdatafab(i-1,j,k,2))/h2

       dsdtfab(i,j,k,3)= &
        -(sdatafab(i,j+1,k,1)-sdatafab(i,j-1,k,1))/(two*h)+ &
         half*dt*(sdatafab(i,j+1,k,3)- &
        two*sdatafab(i,j,k,3)+ &
        sdatafab(i,j-1,k,3))/h2

       dsdtfab(i,j,k,4)= &
        -(sdatafab(i,j,k+1,1)-sdatafab(i,j,k-1,1))/(two*h)+ &
         half*dt*(sdatafab(i,j,k+1,4)- &
        two*sdatafab(i,j,k,4)+ &
        sdatafab(i,j,k-1,4))/h2

      enddo
      enddo
      enddo  !i,j,k 
!DISABLE!$acc end kernels
!DISABLE!$acc end data
#endif

      return
      end subroutine fort_dsdtfab

      end module dsdt_module


