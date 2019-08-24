#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#define STANDALONE 1

#include "REAL.H"
#include "CONSTANTS.H"
#include "SPACE.H"
#include "BC_TYPES.H"
#include "TECPLOTUTIL_F.H"
#include "ArrayLim.H"

#if (BL_SPACEDIM==3)
#define SDIM 3
#elif (BL_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

#if (STANDALONE==1)
      module tecplotutil_cpp_module
      contains
#endif

      subroutine FORT_TECPLOTFAB( &
       time, &
       fabdata,DIMS(fabdata), &
       growlo,growhi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx, &
       dir,ncomp,interior_only,nsteps)

      use global_utility_module

      IMPLICIT NONE

      INTEGER_T dir,ncomp,interior_only,nsteps,bfact
      INTEGER_T growlo(SDIM),growhi(SDIM) 
      INTEGER_T fablo(SDIM),fabhi(SDIM) 
      INTEGER_T plotlo(SDIM),plothi(SDIM) 
      INTEGER_T lo(SDIM),hi(SDIM) 
      INTEGER_T DIMDEC(fabdata)
      REAL_T fabdata(DIMV(fabdata),ncomp)
      REAL_T xlo(SDIM)
      REAL_T xsten(-1:1,SDIM)
      INTEGER_T nhalf
      REAL_T dx(SDIM)
      REAL_T time
      INTEGER_T klo,khi
    

      INTEGER_T ih

      character*17 newfilename
      character*2 matstr
      character*6 stepstr

      INTEGER_T i,j,k,im,nwrite,dir2

! Guibo

      character*80 Title,Varname,Zonename
      REAL*4 ZONEMARKER,EOHMARKER
      integer*4 :: iz_gb,ivar_gb
      integer*4, dimension(:,:), allocatable :: lo_gb,hi_gb

      ! define zone structure
      type zone_t
         real*8, pointer :: var(:,:,:,:)
      end type zone_t
      type(zone_t), dimension(:), allocatable :: zone_gb

! Guibo

      nhalf=1

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif

      if (ncomp.lt.1) then
       print *,"ncomp invalid 30"
       stop
      endif
      if (nsteps.lt.0) then
       print *,"nsteps invalid"
       stop
      endif
      if (time.lt.zero) then
       print *,"time invalid"
       stop
      endif
      write(stepstr,'(I6)') nsteps
      do i=1,6
       if (stepstr(i:i).eq.' ') then
        stepstr(i:i)='0'
       endif
      enddo

      write(newfilename,'(A7,A6,A4)') 'fabdata',stepstr,'.plt'

      if (interior_only.eq.0) then
       do dir2=1,SDIM
        plotlo(dir2)=growlo(dir2)
        plothi(dir2)=growhi(dir2)
       enddo
      else if (interior_only.eq.1) then
       do dir2=1,SDIM
        plotlo(dir2)=fablo(dir2)
        plothi(dir2)=fabhi(dir2)
       enddo
      else
       print *,"interior_only invalid"
       stop
      endif

      nwrite=SDIM+ncomp

      print *,"fabdata: ",newfilename

      !--------------------------------------------------
      ! Determine nzones_gb and allocate zone_gb, lo_gb, hi_gb

      allocate(zone_gb(1))
      allocate(lo_gb(1,SDIM))
      allocate(hi_gb(1,SDIM))

      ! Determine lo_gb, hi_gb
      ! Allocate zone_gb%var later
      iz_gb=1
      do dir2=1,SDIM
       lo_gb(iz_gb,dir2)=plotlo(dir2)
       hi_gb(iz_gb,dir2)=plothi(dir2)
      enddo

      !-----------------------------------------------------------
      ZONEMARKER = 299.0
      EOHMARKER  = 357.0
      open(unit=11,file=newfilename,form="unformatted",access="stream")

      ! +++++++ HEADER SECTION ++++++

      ! i.  Magic number, Version number
      write(11) "#!TDV112"
 
      ! ii. Integer value of 1.
      write(11) 1

      ! iii. Title and variable names.
      ! File type 0 = FULL,1 = GRID,2 = SOLUTION
      write(11) 0
      ! The TITLE
      Title = "CLSVOF data"
      call dumpstring(Title)
      ! Number of variables 
      write(11) nwrite

      ! Variable names.
      Varname='X'
      call dumpstring(Varname)
      Varname='Y'
      call dumpstring(Varname)

      if (SDIM.eq.3) then
       Varname='Z'
       call dumpstring(Varname)
      else if (SDIM.eq.2) then
       ! do nothing
      else
       print *,"dimension bust"
       stop
      endif

      do im=1,ncomp
       write(matstr,'(I2)') im
       do i=1,2
        if (matstr(i:i).eq.' ') then
         matstr(i:i)='0'
        endif
       enddo

       ih=1
       Varname='F'
       ih=ih+1
       Varname(ih:ih)='D'
       ih=ih+1
       do i=1,2
        Varname(ih:ih)=matstr(i:i)
        ih=ih+1
       enddo
       call dumpstring(Varname)
      enddo ! im


      ! Zones
      iz_gb=1
       ! Zone marker. Value = 299.0
      write(11) ZONEMARKER
       ! Zone name
      Zonename = "ZONE"
      call dumpstring(Zonename)
 
      write(11) -1   ! Parent Zone
      write(11) 0    ! StrandID
      write(11) time ! Solution time
      write(11) -1   ! Not used. Set to -1
      write(11) 0    ! Zone Type
      write(11) 0    ! Specify Var Location. 0 = Don't specify, 
                     ! all data is located at the nodes.
      write(11) 0    ! Are raw local 1-to-1 face neighbors supplied?
      write(11) 0    ! Number of miscellaneous user-defined  
                      !face neighbor connections

        ! ----- IMax,JMax,KMax
      write(11) hi_gb(iz_gb,1)-lo_gb(iz_gb,1)+1
      write(11) hi_gb(iz_gb,2)-lo_gb(iz_gb,2)+1
      if (SDIM.eq.3) then
       write(11) hi_gb(iz_gb,SDIM)-lo_gb(iz_gb,SDIM)+1
      else if (SDIM.eq.2) then
       write(11) 1
      else
       print *,"dimension bust"
       stop
      endif
 
      write(11) 0

      write(11) EOHMARKER
 
      ! +++++++ DATA SECTION ++++++

      iz_gb=1

      do dir2=1,SDIM
       lo(dir2)=lo_gb(iz_gb,dir2)
       hi(dir2)=hi_gb(iz_gb,dir2)
      enddo

      if (SDIM.eq.3) then
       klo=lo(SDIM)
       khi=hi(SDIM)
      else if (SDIM.eq.2) then
       klo=0
       khi=0
      else
       print *,"dimension bust"
       stop
      endif

      allocate(zone_gb(iz_gb)% &
        var(nwrite,lo(1):hi(1),lo(2):hi(2),klo:khi))

       ! order is IMPORTANT.
      do k=klo,khi
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)

       if (dir.eq.-1) then
        call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)
       else if ((dir.ge.0).and.(dir.lt.SDIM)) then
        call gridstenMAC(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf,dir+1)
       else
        print *,"dir invalid tecplotfab 2"
        stop
       endif

       zone_gb(iz_gb)%var(1,i,j,k)=xsten(0,1)
       zone_gb(iz_gb)%var(2,i,j,k)=xsten(0,2)
       if (SDIM.eq.3) then
        zone_gb(iz_gb)%var(SDIM,i,j,k)=xsten(0,SDIM)
       else if (SDIM.eq.2) then
        ! do nothing
       else
        print *,"dimension bust"
        stop
       endif

       do im=1,ncomp
        zone_gb(iz_gb)%var(SDIM+im,i,j,k)=fabdata(D_DECL(i,j,k),im)
       enddo
      enddo
      enddo
      enddo


       ! Zone marker  Value = 299.0
      write(11) ZONEMARKER
       ! Data format
      do i=1,nwrite
        write(11) 2
      enddo
      write(11) 0  ! Has passive variables: 0 = no, 1 = yes.
      write(11) 0  ! Has variable sharing 0 = no, 1 = yes.
      write(11) -1 ! Share connectivity list (-1 = no sharing). 
  
      do ivar_gb=1,nwrite
       write(11) minval(zone_gb(1)%var(ivar_gb,:,:,:))
       write(11) maxval(zone_gb(1)%var(ivar_gb,:,:,:))
      enddo

       ! order is IMPORTANT
      do ivar_gb=1,nwrite
       do k=klo,khi
       do j=lo_gb(iz_gb,2),hi_gb(iz_gb,2)
       do i=lo_gb(iz_gb,1),hi_gb(iz_gb,1)
        write(11) zone_gb(iz_gb)%var(ivar_gb,i,j,k)
       enddo
       enddo
       enddo
      enddo

      deallocate(zone_gb(iz_gb)%var)

      deallocate(zone_gb)
      deallocate(lo_gb)
      deallocate(hi_gb)

      close(11)
      return
      end subroutine FORT_TECPLOTFAB

#if (STANDALONE==1)
      end module tecplotutil_cpp_module
#endif

#undef STANDALONE
