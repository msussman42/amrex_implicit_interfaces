#undef BL_LANG_CC
#define BL_LANG_FORT

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else  
print *,"dimension bust"
stop
#endif

module bubbleControl_module
use amrex_fort_module, only : amrex_real

implicit none

      integer, PARAMETER :: maxDrop_pack=400
      real(amrex_real), PARAMETER :: period_pack=0.0004
      real(amrex_real), PARAMETER :: rDrop_pack=0.1
      real(amrex_real), PARAMETER :: uDrop_pack=-600.0
      real(amrex_real) xRange_pack,xShift_pack,xMax_pack,xMin_pack 

      real(amrex_real) xSphere(500)
      real(amrex_real) ySphere(500)
      real(amrex_real) zSphere(500)
      real(amrex_real) rSphere(500)
      integer nSphere
      real(amrex_real) startPosi_pack(maxDrop_pack)
      real(amrex_real) xrand_pack(maxDrop_pack)
      real(amrex_real) LSDrop_pack(maxDrop_pack)
 
contains

      subroutine Pack_velbc(dir,side,veldir,velcell)
      use probcommon_module
      IMPLICIT NONE

      integer dir,side,veldir,dir2
      real(amrex_real) velcell(SDIM)

      do dir2=1,SDIM
       velcell(dir2)=zero
      enddo
      if ((dir.ge.1).and.(dir.le.SDIM).and.(adv_dir.ge.1).and. &
          (adv_dir.le.SDIM).and.(veldir.ge.1).and.(veldir.le.SDIM).and. &
          (side.ge.1).and.(side.le.2)) then
       if ((dir.eq.adv_dir).and.(veldir.eq.adv_dir)) then
        velcell(veldir)=adv_vel
       endif
      else
       print *,"parameters invalid Pack_velbc"
       stop
      endif

      return
      end subroutine Pack_velbc

      subroutine get_pack_vfrac(x,y,z,dx,vfrac,cenbc,time,dir)
      use probcommon_module
      IMPLICIT NONE

      integer dir
      real(amrex_real) x,y,z,time
      integer im
      real(amrex_real) dx(SDIM)
      real(amrex_real) vfrac(num_materials)
      real(amrex_real) cenbc(num_materials,SDIM)

      real(amrex_real) LS
      real(amrex_real) xp(5)
      integer iFound,ip

      if (dir.eq.2) then   ! y-direction
          !LS=0.05    ! liquid jet diameter   ysl is this correct?
          LS=0.02/2.0
          xp(1)=0.05 !0.55 !0.3535
          !xp(2)=0.7070
          !xp(3)=0.85
          !xp(5)=4.4
          ip=1
          iFound=0
          do while (ip.le.1.and.iFound.eq.0)
            if (abs(x-xp(ip)).le.LS) then ! GAS JET (material 1 is liquid)
              vfrac(1)=zero
              iFound=1
            else
              vfrac(1)=one
            endif
            ip=ip+1
          enddo
      else
         print*,"wrong type vfrac, stop"
         stop
      endif

      vfrac(2)=1.0-vfrac(1)
      do im=3,num_materials
         vfrac(im)=zero
      enddo

      return
      end subroutine get_pack_vfrac


      subroutine xloLS_pack(x,y,z,VOF,time,bigdist)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real) x,y,z,time
      real(amrex_real) VOF(num_materials)
      real(amrex_real) LS

      integer nDrop,iLoc,nLoc,insiderDrop
      real(amrex_real) bigdist,LS0
      real(amrex_real) nowTime,rainDuTime
      integer im

      rainDuTime=2.0*rDrop_pack/abs(adv_vel)
      nLoc=ceiling(rainDuTime/period_pack)
      nowTime=MOD(time,rainDuTime)
      !nowTime=MOD(time,period_pack)
      Call random_seed
      call random_number(startPosi_pack)
      nDrop=mod(ceiling(time/rainDuTime),maxDrop_pack)
      xrand_pack(1)=startPosi_pack(nDrop)*(xMax_pack-xMin_pack)+xMin_pack
      do iLoc=2,nLoc
        xrand_pack(iLoc)=xrand_pack(1)+ &
          (xMax_pack-xMin_pack)*(iLoc-1.0)/(nLoc-1.0)
        if (xrand_pack(iLoc).gt.xMax_pack) then
         xrand_pack(iLoc)=xrand_pack(iLoc)-(xMax_pack-xMin_pack)
        endif
      enddo
      if (nowTime.le.rainDuTime) then
       LS0=sqrt(rDrop_pack*rDrop_pack- &
           (rDrop_pack-abs(adv_vel)*nowTime)**2)
       ! there are nLoc drops here, need to find the level set distance
       do iLoc=1,nLoc
        LSDrop_pack(iLoc)=LS0-abs(y-(xrand_pack(iLoc)+ &
         uDrop_pack*nowTime))
       enddo
      
       insiderDrop=0
       do iLoc=1,nLoc
         if (LSDrop_pack(iLoc).ge.0.0) then  ! assume drops are not close
           insiderDrop=1
           VOF(1)=LSDrop_pack(iLoc)
           VOF(2)=-LSDrop_pack(iLoc)
           !print*, "iLoc2=", iLOC, "xra=",xrand_pack(iLoc),"adv",adv_vel
         endif
       enddo
     
       if (insiderDrop.eq.0) then ! point not in any drops find the smallest
         LS=bigdist
         do iLoc=1,nLoc
           if (abs(LSDrop_pack(iLoc)).le.abs(LS)) then
            LS=LSDrop_pack(iLoc)
           endif
         enddo
         VOF(1)=LS
         VOF(2)=-LS
       endif

       do im=3,num_materials
        VOF(im)=-bigdist
       enddo       
      endif 

      end subroutine xloLS_pack


      subroutine yloLS_pack(x,y,z,VOF,time,bigdist)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real) x,y,z,time
      real(amrex_real) VOF(num_materials)
      real(amrex_real) LS

      real(amrex_real) bigdist
      real(amrex_real) xp(5)
      integer im


      !LS=0.3-sqrt(x*x)  ! 0.3 is the gas jet diameter
      !VOF(1)=LS
      !VOF(2)=-LS
      xp(1)=.05
      !xp(2)=-0.5
      !xp(3)=0.0
      !xp(4)=0.5
      !xp(5)=1.0

      !LS=-1.0e9
      !do i=1,5
      !   LS=max(1.9/10.0-abs(x-xp(i)),LS)
      !enddo 
       
      LS=0.01-abs(x-xp(1))   ! 
      VOF(1)=-LS  ! GAS JET
      VOF(2)=LS 
      do im=3,num_materials
        VOF(im)=-bigdist
      enddo

      end subroutine yloLS_pack



end module bubbleControl_module
