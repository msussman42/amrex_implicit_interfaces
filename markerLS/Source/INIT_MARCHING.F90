#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

! get rid of autoindent   :setl noai nocin nosi inde=

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif



      subroutine fort_init_marching( &
       xlo,dx, &
       n_cell) &
      bind(c,name='fort_init_marching')

      use probcommon_module
      use probmain_module
      use global_utility_module
      use LegendreNodes

      IMPLICIT NONE

      INTEGER_T, intent(in) :: n_cell(SDIM)
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T bfactmax
      INTEGER_T dir
      INTEGER_T i
      INTEGER_T ilev
      INTEGER_T inode
      INTEGER_T max_ncell
      INTEGER_T domlo_in(3)
      REAL_T problo_arr(3)
      INTEGER_T nhalf
      INTEGER_T sci_max_level
      REAL_T dx_local(3)
      REAL_T xsten_cache(-1:1)

      problox=xlo(1)
      probloy=xlo(2)
      probhix=problox+n_cell(1)*dx(1)
      probhiy=probloy+n_cell(2)*dx(2)

      probloz=probloy
      probhiz=probhiy
      if (SDIM.eq.3) then
       probloz=xlo(SDIM)
       probhiz=probloz+n_cell(SDIM)*dx(SDIM)
      endif
      problenx=probhix-problox
      probleny=probhiy-probloy
      problenz=probhiz-probloz

      sci_max_level=0
      fort_max_level=0
      fort_finest_level=0

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

      max_ncell=n_cell(1)
      if (max_ncell.lt.n_cell(2)) then
       max_ncell=n_cell(2)
      endif
      if (max_ncell.lt.n_cell(SDIM)) then
       max_ncell=n_cell(SDIM)
      endif
              
      do ilev=1,cache_max_level
       max_ncell=max_ncell*2
      enddo

      cache_index_low=-4*bfactmax
      cache_index_high=2*max_ncell+4*bfactmax

      allocate(grid_cache(0:cache_max_level, &
       cache_index_low:cache_index_high,SDIM))

      allocate(dxlevel(0:cache_max_level,SDIM))
      allocate(domlo_level(0:cache_max_level,SDIM))
      allocate(domhi_level(0:cache_max_level,SDIM))

      domlo_in(1)=0
      domlo_in(2)=0
      domlo_in(3)=0

      do dir=1,SDIM
       dxlevel(0,dir)=dx(dir)
       domlo_level(0,dir)=0
       domhi_level(0,dir)=n_cell(dir)-1
      enddo

      do ilev=0,cache_max_level
       do dir=1,SDIM
        dx_local(dir)=dxlevel(ilev,dir)
       enddo
       do dir=1,SDIM
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
         problo_arr(3)=0.0d0
         call gridsten1D(xsten_cache,problo_arr,i, &
           domlo_in, &
           bfact_space_order(0),dx_local,dir,nhalf) 
         grid_cache(ilev,inode,dir)=xsten_cache(0)
         grid_cache(ilev,inode+1,dir)=xsten_cache(1)
        enddo ! i
       enddo ! dir=1..SDIM

       if (ilev.lt.cache_max_level) then
        do dir=1,SDIM
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


      return
      end subroutine fort_init_marching
