
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"
#include "INTERP_F.H"
#include "EXTRAP_COMP.H"

#define POLYGON_LIST_MAX (1000)


#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else  
print *,"dimension bust"
stop
#endif

      module interp_module
      use amrex_fort_module, only : amrex_real

      contains

      subroutine fort_override_finest_level(cc_finest_level) &
      bind(c,name='fort_override_finest_level')

      use probcommon_module

      IMPLICIT NONE
      
      integer, INTENT(in) :: cc_finest_level

      if ((cc_finest_level.lt.0).or.(cc_finest_level.gt.1000)) then
       print *,"cc_finest_level invalid"
       stop
      endif

      fort_finest_level=cc_finest_level

      return
      end subroutine fort_override_finest_level

      subroutine fort_gl_slab( &
       time_array, &
       slab_dt_type, &
       cc_time_order, &
       slablow,slabhigh) &
      bind(c,name='fort_gl_slab')

      use LegendreNodes
      use probcommon_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: slab_dt_type
      integer, INTENT(in) :: cc_time_order
      integer i1
      real(amrex_real), INTENT(out) :: time_array(0:cc_time_order)
      real(amrex_real), INTENT(in) :: slablow,slabhigh
      real(amrex_real) yGL(0:cc_time_order)
      real(amrex_real) slab_dt

      if (cc_time_order.ne.bfact_time_order) then
       print *,"cc_time_order invalid"
       stop
      endif
      if (cc_time_order.lt.1) then
       print *,"cc_time_order out of range"
       stop
      endif

      do i1=0,bfact_time_order
       if (slab_dt_type.eq.0) then
        yGL(i1)=cache_gauss_lobatto(bfact_time_order,i1,TMTYPE)
       else if (slab_dt_type.eq.1) then
        slab_dt=two/bfact_time_order
        yGL(i1)=-one+i1*slab_dt
       else
        print *,"slab_dt_type invalid"
        stop
       endif
       if ((yGL(i1).lt.-one).or.(yGL(i1).gt.one)) then
        print *,"yGL invalid in gl_slab"
        stop
       endif
       if (i1.gt.0) then
        if (yGL(i1).le.yGL(i1-1)) then
         print *,"yGL invalid in gl_slab2"
         stop
        endif
       endif
      enddo  ! i1

      if ((yGL(0)+one.ne.zero).or. &
          (yGL(bfact_time_order)-one.ne.zero)) then
       print *,"yGL invalid in gl_slab3"
       stop
      endif

      do i1=0,bfact_time_order
       time_array(i1)=slablow+(yGL(i1)+one)*(slabhigh-slablow)*half
      enddo

      return
      end subroutine fort_gl_slab

      subroutine fort_multimofinterp( &
       tid_in, &
       time, &
       datamof,DIMS(dmof), &
       clo,chi, &
       fdatamof,DIMS(fdmof), &
       flo,fhi, &
       datarecon,DIMS(drecon), &
       problo,dxf,dxc, &
       ngeom_recon_test,ngeom_raw_test, &
       levelc,levelf, &
       bfact_coarse,bfact_fine) &
      bind(c,name='fort_multimofinterp')

      use geometry_intersect_module
      use MOF_routines_module
      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: tid_in
      integer, INTENT(in) :: levelc,levelf
      integer, INTENT(in) :: bfact_coarse,bfact_fine
      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: ngeom_recon_test,ngeom_raw_test
      integer, INTENT(in) :: DIMDEC(dmof)
      integer, INTENT(in) :: DIMDEC(drecon)
      integer, INTENT(in) :: DIMDEC(fdmof)
      integer, INTENT(in) :: clo(SDIM),chi(SDIM)
      integer, INTENT(in) :: flo(SDIM),fhi(SDIM)
      integer domlo(SDIM)
      real(amrex_real), INTENT(in) :: &
              datamof(DIMV(dmof),num_materials*ngeom_raw)
      real(amrex_real), INTENT(out) :: &
              datarecon(DIMV(drecon),num_materials*ngeom_recon)
      real(amrex_real), INTENT(out) :: &
              fdatamof(DIMV(fdmof),num_materials*ngeom_raw)
      real(amrex_real), INTENT(in) :: problo(SDIM)
      real(amrex_real), INTENT(in) :: dxf(SDIM),dxc(SDIM)
      integer crse_growlo(3),crse_growhi(3)
      integer growlo(3),growhi(3)
      integer stenlo(3),stenhi(3)
      real(amrex_real) wt(SDIM)

      integer :: grid_index(SDIM)
      integer :: data_needed

      integer, parameter :: grid_level=-1

      integer i,j,k
      integer ifine,jfine,kfine
      integer ic,jc,kc
      real(amrex_real) testwt

      integer dir
      integer nmax,im,vofcomp_old,vofcomp_new
      real(amrex_real) mofdata(num_materials*ngeom_recon)
      real(amrex_real) mofdatafine(num_materials*ngeom_recon)
      real(amrex_real) multi_centroidA(num_materials,SDIM)
      real(amrex_real) volcell
      real(amrex_real) cencell(SDIM)
      real(amrex_real) volfine
      real(amrex_real) cenfine(SDIM)
      real(amrex_real) voltemp
      real(amrex_real) vof_super(num_materials)
      real(amrex_real) multi_volume(num_materials)
      real(amrex_real) multi_cen(SDIM,num_materials)
      integer, PARAMETER :: continuous_mof=STANDARD_MOF
      integer cmofsten(D_DECL(-1:1,-1:1,-1:1))
      integer, PARAMETER :: use_ls_data=0
      integer, PARAMETER ::  mof_verbose=0
      real(amrex_real) LS_stencil(D_DECL(-1:1,-1:1,-1:1),num_materials)
      integer, parameter :: nhalf=3
      integer, parameter :: nhalfgrid=1
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xstenfine(-nhalf:nhalf,SDIM)
      real(amrex_real) xstengrid(-nhalfgrid:nhalfgrid,SDIM)
      integer n_overlap
      integer tessellate

      if ((tid_in.ge.geom_nthreads).or.(tid_in.lt.0)) then
       print *,"tid_in invalid: ",tid_in
       stop
      endif 

      tessellate=0

      if (ngeom_recon.ne.2*SDIM+3) then
       print *,"ngeom_recon invalid"
       stop
      endif
      if (ngeom_recon_test.ne.2*SDIM+3) then
       print *,"ngeom_recon_test invalid"
       stop
      endif
      if (ngeom_raw.ne.SDIM+1) then
       print *,"ngeom_raw invalid"
       stop
      endif
      if (ngeom_raw_test.ne.SDIM+1) then
       print *,"ngeom_raw_test invalid"
       stop
      endif

      if ((levelc.ne.levelf-1).or.(levelc.lt.0)) then
       print *,"levelc or levelf invalid"
       stop
      endif
      if (levelf.gt.fort_finest_level) then
       print *,"levelf invalid"
       stop
      endif
      if (bfact_coarse.lt.1) then
       print *,"bfact_coarse invalid"
       stop
      endif
      if (bfact_fine.lt.1) then
       print *,"bfact_fine invalid"
       stop
      endif
      if (bfact_fine.gt.bfact_coarse) then
       print *,"bfact_fine invalid"
       stop
      endif

      do dir=1,SDIM
       domlo(dir)=0
      enddo

      nmax=POLYGON_LIST_MAX

      call growntilebox(clo,chi,clo,chi,crse_growlo,crse_growhi,0)
      call growntilebox(flo,fhi,flo,fhi,growlo,growhi,0) 

      do dir=1,SDIM
       if ((2*clo(dir).le.flo(dir)).and. &
           (2*(chi(dir)+1).ge.fhi(dir)+1)) then
        !do nothing
       else
        print *,"clo,chi,flo,fhi mismatch: ",clo,chi,flo,fhi 
        stop
       endif
      enddo !dir=1,sdim

      do k=crse_growlo(3),crse_growhi(3)
      do j=crse_growlo(2),crse_growhi(2)
      do i=crse_growlo(1),crse_growhi(1)

       grid_index(1)=i
       grid_index(2)=j
       if (SDIM.eq.3) then
        grid_index(SDIM)=k
       endif

        ! coarse cell ic covers 2*ic and 2*ic+1
       data_needed=1
       do dir=1,SDIM
        if ((2*grid_index(dir)+1.ge.growlo(dir)).and. &
            (2*grid_index(dir).le.growhi(dir))) then
         !do nothing
        else
         data_needed=0
        endif
       enddo !dir=1,SDIM

       if (data_needed.eq.1) then

        do im=1,num_materials
         vofcomp_old=(im-1)*ngeom_raw+1
         vofcomp_new=(im-1)*ngeom_recon+1
         mofdata(vofcomp_new)=datamof(D_DECL(i,j,k),vofcomp_old)
         do dir=1,SDIM
          mofdata(vofcomp_new+dir)=datamof(D_DECL(i,j,k),vofcomp_old+dir)
           !slope=0
          mofdata(vofcomp_new+SDIM+1+dir)=zero
         enddo

           ! order=0
         mofdata(vofcomp_new+SDIM+1)=zero

        enddo  ! im

!       call gridsten(xsten,problo,i,j,k,domlo,bfact_coarse,dxc,nhalf)
        call gridsten_level(xsten,i,j,k,levelc,nhalf)

         ! sum F_fluid=1  sum F_solid <= 1
        call make_vfrac_sum_ok_base( &
          cmofsten, &
          xsten,nhalf, &
          continuous_mof, &
          bfact_coarse,dxc, &
          tessellate, & !=0
          mofdata,SDIM)

        do im=1,num_materials
         vofcomp_new=(im-1)*ngeom_recon+1
         vof_super(im)=mofdata(vofcomp_new)
        enddo

        call multimaterial_MOF( &
          tid_in, &
          bfact_coarse,dxc,xsten,nhalf, &
          mof_verbose, &
          use_ls_data, & ! use_ls_data=0
          LS_stencil, &
          geom_xtetlist(1,1,1,tid_in+1), &
          geom_xtetlist(1,1,1,tid_in+1), &
          nmax, &
          nmax, &
          mofdata, & !intent(inout)
          vof_super, &
          multi_centroidA, &
          continuous_mof, & ! continuous_mof=STANDARD_MOF
          cmofsten, &
          grid_index, &
          grid_level, & !grid_level=-1
          SDIM)

        do dir=1,num_materials*ngeom_recon
         datarecon(D_DECL(i,j,k),dir)=mofdata(dir)
        enddo

       else if (data_needed.eq.0) then

        do dir=1,num_materials*ngeom_recon
         datarecon(D_DECL(i,j,k),dir)=zero
        enddo

       else
        print *,"data_needed invalid: ",data_needed
        stop
       endif

      enddo
      enddo
      enddo ! i,j,k


      do kfine=growlo(3),growhi(3)
      do jfine=growlo(2),growhi(2)
      do ifine=growlo(1),growhi(1)

       call gridsten(xstenfine,problo,ifine,jfine,kfine, &
         domlo,bfact_fine,dxf,nhalf)
       call Box_volumeFAST(bfact_fine,dxf,xstenfine,nhalf, &
         volfine,cenfine,SDIM)

       do dir=1,num_materials*ngeom_recon
        mofdatafine(dir)=zero
       enddo

       volcell=zero
       do dir=1,SDIM
        cencell(dir)=zero
       enddo

       n_overlap=0

       call coarse_subelement_stencil(ifine,jfine,kfine,stenlo,stenhi, &
         bfact_coarse,bfact_fine)

       do ic=stenlo(1),stenhi(1)
        call intersect_weight_interp(ic,ifine, &
          bfact_coarse,bfact_fine,wt(1))
        if (wt(1).gt.zero) then
         do jc=stenlo(2),stenhi(2)
          call intersect_weight_interp(jc,jfine, &
           bfact_coarse,bfact_fine,wt(2))
          if (wt(2).gt.zero) then
           do kc=stenlo(3),stenhi(3)
            if (SDIM.eq.3) then
             call intersect_weight_interp(kc,kfine, &
              bfact_coarse,bfact_fine,wt(SDIM))
            endif
            if (wt(SDIM).gt.zero) then
             if (SDIM.eq.2) then
              testwt=wt(1)*wt(2)
             else if (SDIM.eq.3) then
              testwt=wt(1)*wt(2)*wt(SDIM)
             else
              print *,"dimension bust"
              stop
             endif
             if (testwt.gt.zero) then

              grid_index(1)=ic
              grid_index(2)=jc
              if (SDIM.eq.3) then
               grid_index(SDIM)=kc
              endif

              data_needed=1
              do dir=1,SDIM
               if ((2*grid_index(dir)+1.ge.growlo(dir)).and. &
                   (2*grid_index(dir).le.growhi(dir))) then
                !do nothing
               else
                data_needed=0
               endif
              enddo !dir=1,SDIM

              if (data_needed.eq.1) then
               !do nothing
              else
               print *,"data_needed bust"
               print *,"grid_index=",grid_index
               print *,"ifine,jfine,kfine ",ifine,jfine,kfine
               print *,"growlo ",growlo
               print *,"growhi ",growhi
               print *,"crse_growlo ",crse_growlo
               print *,"crse_growhi ",crse_growhi
               print *,"wt=",wt
               print *,"testwt=",testwt
               stop
              endif
              n_overlap=n_overlap+1

              do dir=1,num_materials*ngeom_recon
               mofdata(dir)=datarecon(D_DECL(ic,jc,kc),dir)
              enddo

              call gridsten(xsten,problo,ic,jc,kc, &
               domlo,bfact_coarse,dxc,nhalf)

              do dir=1,SDIM
               xstengrid(-1,dir)=max(xsten(-1,dir),xstenfine(-1,dir))
               xstengrid(1,dir)=min(xsten(1,dir),xstenfine(1,dir))
               xstengrid(0,dir)=half*(xstengrid(-1,dir)+xstengrid(1,dir))
               if ((bfact_fine.eq.1).and.(bfact_coarse.eq.1)) then
                if ((abs(xstengrid(-1,dir)-xstenfine(-1,dir)).gt. &
                     VOFTOL*dxf(dir)).or. &
                    (abs(xstengrid(1,dir)-xstenfine(1,dir)).gt. &
                     VOFTOL*dxf(dir))) then
                 print *,"fine cell should be completely covered by coarse"
                 stop
                endif
               endif
              enddo ! dir

              call multi_get_volume_grid_simple( &
               tid_in, &
               EPS3, &
               tessellate, &  !=0
               bfact_coarse,dxc,xsten,nhalf, &
               mofdata, &
               xstengrid,nhalfgrid, &
               multi_volume,multi_cen, &
               geom_xtetlist(1,1,1,tid_in+1), &
               nmax, &
               nmax, &
               SDIM)

              do im=1,num_materials
               vofcomp_new=(im-1)*ngeom_recon+1
               mofdatafine(vofcomp_new)=mofdatafine(vofcomp_new)+ &
                multi_volume(im)
               do dir=1,SDIM
                mofdatafine(vofcomp_new+dir)= &
                 mofdatafine(vofcomp_new+dir)+ &
                 multi_cen(dir,im)*multi_volume(im)
               enddo

               if (is_rigid(im).eq.0) then
                volcell=volcell+multi_volume(im)
                do dir=1,SDIM
                 cencell(dir)=cencell(dir)+ &
                  multi_cen(dir,im)*multi_volume(im)
                enddo
               else if (is_rigid(im).eq.1) then
                ! do nothing
               else
                print *,"is_rigid(im) invalid"
                stop
               endif
              enddo ! im=1..num_materials
             else if (testwt.eq.zero) then
              print *,"testwt should not be 0"
              stop
             else
              print *,"testwt invalid"
              stop
             endif
            else if (wt(SDIM).eq.zero) then
             ! do nothing
            else 
             print *,"wt(SDIM) invalid"
             stop
            endif
           enddo ! kc
          endif
         enddo ! jc
        endif
       enddo ! ic

       if ((bfact_fine.eq.1).and.(bfact_coarse.eq.1)) then
        if (n_overlap.ne.1) then
         print *,"n_overlap invalid"
         stop
        endif
       endif

       if (volcell.gt.zero) then
        !do nothing
       else
        print *,"volcell must be positive: ",volcell
        stop
       endif
       do dir=1,SDIM
        cencell(dir)=cencell(dir)/volcell
       enddo
       if (abs(volcell-volfine).le.EPS2*volfine) then
        ! do nothing
       else
        print *,"volcell, volfine bad (multimofinterp):",volcell,volfine
        stop
       endif

       do im=1,num_materials
        vofcomp_old=(im-1)*ngeom_raw+1
        vofcomp_new=(im-1)*ngeom_recon+1
        voltemp=mofdatafine(vofcomp_new)/volcell
        fdatamof(D_DECL(ifine,jfine,kfine),vofcomp_old)=voltemp
        if (voltemp.gt.zero) then
         do dir=1,SDIM
          fdatamof(D_DECL(ifine,jfine,kfine),vofcomp_old+dir)= &
            mofdatafine(vofcomp_new+dir)/mofdatafine(vofcomp_new)- &
            cencell(dir)
         enddo
        else
         do dir=1,SDIM
          fdatamof(D_DECL(ifine,jfine,kfine),vofcomp_old+dir)=zero
         enddo 
        endif
       enddo ! im

      enddo
      enddo
      enddo ! ifine,jfine,kfine
 
      return 
      end subroutine fort_multimofinterp


      subroutine fort_lsinterp( &
       tid_in, &
       clsdata,DIMS(clsdata), &
       clo,chi, &
       flsdata,DIMS(flsdata), &
       flo,fhi, &
       problo, &
       dxf,dxc, &
       ncomp, &
       levelc,levelf, &
       bfact_coarse,bfact_fine) &
      bind(c,name='fort_lsinterp')

      use geometry_intersect_module
      use MOF_routines_module
      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: tid_in
      integer, INTENT(in) :: levelc,levelf
      integer, INTENT(in) :: bfact_coarse,bfact_fine
      integer, INTENT(in) :: ncomp
      integer, INTENT(in) :: DIMDEC(clsdata)
      integer, INTENT(in) :: DIMDEC(flsdata)
      integer, INTENT(in) :: clo(SDIM),chi(SDIM)
      integer, INTENT(in) :: flo(SDIM),fhi(SDIM)
      integer domlo(SDIM)
      real(amrex_real), INTENT(in) :: clsdata(DIMV(clsdata),ncomp)
      real(amrex_real), INTENT(out) :: flsdata(DIMV(flsdata),ncomp)
      real(amrex_real), INTENT(in) :: problo(SDIM)
      real(amrex_real), INTENT(in) :: dxf(SDIM),dxc(SDIM)
      integer growlo(3),growhi(3)
      integer stenlo(3),stenhi(3)
      real(amrex_real) wt(SDIM)

      integer ifine,jfine,kfine
      integer ic,jc,kc
      real(amrex_real) testwt

      integer dir
      real(amrex_real) volcell
      real(amrex_real) cencell(SDIM)
      real(amrex_real) volfine
      real(amrex_real) cenfine(SDIM)
      real(amrex_real) volcoarse
      real(amrex_real) cencoarse(SDIM)
      real(amrex_real) voltemp
      real(amrex_real) centemp(SDIM)
      integer, parameter :: nhalf=3
      integer, parameter :: nhalfgrid=1
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xstenfine(-nhalf:nhalf,SDIM)
      real(amrex_real) xstengrid(-nhalfgrid:nhalfgrid,SDIM)
      integer n_overlap
      real(amrex_real) LS_FINE(ncomp)
      real(amrex_real) LS_COARSE(ncomp)

      if ((tid_in.ge.geom_nthreads).or.(tid_in.lt.0)) then
       print *,"tid_in invalid: ",tid_in
       stop
      endif

      if (ncomp.ne.(SDIM+1)*num_materials) then
       print *,"ncomp invalid"
       stop
      endif
      if ((levelc.ne.levelf-1).or.(levelc.lt.0)) then
       print *,"levelc or levelf invalid"
       stop
      endif
      if (levelf.gt.fort_finest_level) then
       print *,"levelf invalid"
       stop
      endif
      if (bfact_coarse.lt.1) then
       print *,"bfact_coarse invalid"
       stop
      endif
      if (bfact_fine.lt.1) then
       print *,"bfact_fine invalid"
       stop
      endif
      if (bfact_fine.gt.bfact_coarse) then
       print *,"bfact_fine invalid"
       stop
      endif

      do dir=1,SDIM
       domlo(dir)=0
      enddo

      call growntilebox(flo,fhi,flo,fhi,growlo,growhi,0) 

      do kfine=growlo(3),growhi(3)
      do jfine=growlo(2),growhi(2)
      do ifine=growlo(1),growhi(1)

       call gridsten(xstenfine,problo,ifine,jfine,kfine, &
         domlo,bfact_fine,dxf,nhalf)
       call Box_volumeFAST(bfact_fine,dxf,xstenfine,nhalf, &
         volfine,cenfine,SDIM)

       do dir=1,ncomp
        LS_FINE(dir)=zero
       enddo

       volcell=zero
       do dir=1,SDIM
        cencell(dir)=zero
       enddo

       n_overlap=0

       call coarse_subelement_stencil(ifine,jfine,kfine,stenlo,stenhi, &
         bfact_coarse,bfact_fine)

       do ic=stenlo(1),stenhi(1)
        call intersect_weight_interp(ic,ifine, &
          bfact_coarse,bfact_fine,wt(1))
        if (wt(1).gt.zero) then
         do jc=stenlo(2),stenhi(2)
          call intersect_weight_interp(jc,jfine, &
           bfact_coarse,bfact_fine,wt(2))
          if (wt(2).gt.zero) then
           do kc=stenlo(3),stenhi(3)
            if (SDIM.eq.3) then
             call intersect_weight_interp(kc,kfine, &
              bfact_coarse,bfact_fine,wt(SDIM))
            endif
            if (wt(SDIM).gt.zero) then
             if (SDIM.eq.2) then
              testwt=wt(1)*wt(2)
             else if (SDIM.eq.3) then
              testwt=wt(1)*wt(2)*wt(SDIM)
             else
              print *,"dimension bust"
              stop
             endif
             if (testwt.gt.zero) then
 
              n_overlap=n_overlap+1

              do dir=1,ncomp
               LS_COARSE(dir)=clsdata(D_DECL(ic,jc,kc),dir)
              enddo

              call gridsten(xsten,problo,ic,jc,kc, &
               domlo,bfact_coarse,dxc,nhalf)

              call Box_volumeFAST(bfact_coarse,dxc,xsten,nhalf, &
                 volcoarse,cencoarse,SDIM)

              do dir=1,SDIM
               xstengrid(-1,dir)=max(xsten(-1,dir),xstenfine(-1,dir))
               xstengrid(1,dir)=min(xsten(1,dir),xstenfine(1,dir))
               xstengrid(0,dir)=half*(xstengrid(-1,dir)+xstengrid(1,dir))
               if ((bfact_fine.eq.1).and.(bfact_coarse.eq.1)) then
                if ((abs(xstengrid(-1,dir)-xstenfine(-1,dir)).gt. &
                     VOFTOL*dxf(dir)).or. &
                    (abs(xstengrid(1,dir)-xstenfine(1,dir)).gt. &
                     VOFTOL*dxf(dir))) then
                 print *,"fine cell should be completely covered by coarse"
                 stop
                endif
               endif
              enddo ! dir

              call Box_volumeFAST(bfact_fine,dxf,xstengrid,nhalfgrid, &
                 voltemp,centemp,SDIM)
              volcell=volcell+voltemp
              do dir=1,SDIM
               cencell(dir)=cencell(dir)+voltemp*centemp(dir)
              enddo

               ! normals
              do dir=num_materials+1,ncomp
               LS_FINE(dir)=LS_FINE(dir)+voltemp*LS_COARSE(dir) 
              enddo

               ! levelsets
              do dir=1,num_materials
               LS_FINE(dir)=LS_FINE(dir)+voltemp*LS_COARSE(dir)
              enddo ! dir=1..num_materials
 
             else if (testwt.eq.zero) then
              print *,"testwt cannot be 0"
              stop
             else
              print *,"testwt invalid"
              stop
             endif
            else if (wt(SDIM).eq.zero) then
             ! do nothing
            else 
             print *,"wt(SDIM) invalid"
             stop
            endif
           enddo ! kc
          endif
         enddo ! jc
        endif
       enddo ! ic

       if ((bfact_fine.eq.1).and.(bfact_coarse.eq.1)) then
        if (n_overlap.ne.1) then
         print *,"n_overlap invalid"
         stop
        endif
       endif

       if (volcell.gt.zero) then
        !do nothing
       else
        print *,"volcell must be positive: ",volcell
        stop
       endif
       do dir=1,SDIM
        cencell(dir)=cencell(dir)/volcell
       enddo 
       if (abs(volcell-volfine).le.EPS2*volfine) then
        !do nothing
       else
        print *,"volcell, volfine bad (lshointerp):",volcell,volfine
        stop
       endif

       do dir=1,ncomp
        flsdata(D_DECL(ifine,jfine,kfine),dir)=LS_FINE(dir)/volcell
       enddo 

      enddo
      enddo
      enddo ! ifine,jfine,kfine
 
      return 
      end subroutine fort_lsinterp



      ! in NavierStokes::VOF_Recon
      ! 1. get MOF data with 1 ghost cell (so that CMOF can be chosen)
      ! 2. reconstruct interior cells only.
      ! 3. do extended filpatch; MOF used for coarse/fine and ext_dir cells.
      subroutine fort_multiextmofinterp( &
       tid_in, &
       time, &
       datamof,DIMS(dmof), &
       fdatamof,DIMS(fdmof), &
       flo,fhi, &
       problo, &
       dxf,dxc, &
       ngeom_recon_test,ngeom_raw_test, &
       levelc,levelf, &
       bfact_coarse,bfact_fine) &
      bind(c,name='fort_multiextmofinterp')

      use geometry_intersect_module
      use MOF_routines_module
      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: tid_in
      integer, INTENT(in) :: levelc,levelf
      integer, INTENT(in) :: bfact_coarse,bfact_fine
      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: ngeom_recon_test,ngeom_raw_test
      integer, INTENT(in) :: DIMDEC(dmof)
      integer, INTENT(in) :: DIMDEC(fdmof)
      integer, INTENT(in) :: flo(SDIM),fhi(SDIM)
      integer domlo(SDIM)
      real(amrex_real), INTENT(in) :: &
        datamof(DIMV(dmof),num_materials*ngeom_recon)
      real(amrex_real), INTENT(out) :: &
        fdatamof(DIMV(fdmof),num_materials*ngeom_recon)
      real(amrex_real), INTENT(in) :: problo(SDIM),dxf(SDIM),dxc(SDIM)
      integer growlo(3),growhi(3)
      integer stenlo(3),stenhi(3)
      real(amrex_real) wt(SDIM)

      integer ifine,jfine,kfine
      integer ic,jc,kc
      real(amrex_real) testwt

      integer :: grid_index(SDIM)
      integer, parameter :: grid_level=-1

      integer dir
      integer nmax,im,vofcomp_old,vofcomp_new
      real(amrex_real) mofdata(num_materials*ngeom_recon)
      real(amrex_real) mofdatafine(num_materials*ngeom_recon)
      real(amrex_real) volcell
      real(amrex_real) cencell(SDIM)
      real(amrex_real) volfine
      real(amrex_real) cenfine(SDIM)
      real(amrex_real) voltemp
      real(amrex_real) vof_super(num_materials)
      real(amrex_real) multi_volume(num_materials)
      real(amrex_real) multi_cen(SDIM,num_materials)
      integer, parameter :: nhalf=3
      integer, parameter :: nhalfgrid=1
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xstenfine(-nhalf:nhalf,SDIM)
      real(amrex_real) xstengrid(-nhalfgrid:nhalfgrid,SDIM)
      integer n_overlap
      integer tessellate

      integer, PARAMETER :: use_ls_data=0
      integer, PARAMETER ::  mof_verbose=0
      integer, PARAMETER :: continuous_mof=STANDARD_MOF
      integer cmofsten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) multi_centroidA(num_materials,SDIM)
      real(amrex_real) LS_stencil(D_DECL(-1:1,-1:1,-1:1),num_materials)

      if ((tid_in.ge.geom_nthreads).or.(tid_in.lt.0)) then
       print *,"tid_in invalid: ",tid_in
       stop
      endif 

      tessellate=0

      if (ngeom_recon.ne.2*SDIM+3) then
       print *,"ngeom_recon invalid"
       stop
      endif
      if (ngeom_raw.ne.SDIM+1) then
       print *,"ngeom_raw invalid"
       stop
      endif
      if (ngeom_recon_test.ne.2*SDIM+3) then
       print *,"ngeom_recon_test invalid"
       stop
      endif
      if (ngeom_raw_test.ne.SDIM+1) then
       print *,"ngeom_raw_test invalid"
       stop
      endif

      if ((levelc.ne.levelf-1).or.(levelc.lt.0)) then
       print *,"levelc or levelf invalid"
       stop
      endif
      if (levelf.gt.fort_finest_level) then
       print *,"levelf invalid"
       stop
      endif
      if (bfact_coarse.lt.1) then
       print *,"bfact_coarse invalid"
       stop
      endif
      if (bfact_fine.lt.1) then
       print *,"bfact_fine invalid"
       stop
      endif
      if (bfact_fine.gt.bfact_coarse) then
       print *,"bfact_fine invalid"
       stop
      endif

      do dir=1,SDIM
       domlo(dir)=0
      enddo

      nmax=POLYGON_LIST_MAX

      call growntilebox(flo,fhi,flo,fhi,growlo,growhi,0) 

      do kfine=growlo(3),growhi(3)
      do jfine=growlo(2),growhi(2)
      do ifine=growlo(1),growhi(1)

       call gridsten(xstenfine,problo,ifine,jfine,kfine, &
         domlo,bfact_fine,dxf,nhalf)
       call Box_volumeFAST(bfact_fine,dxf,xstenfine,nhalf, &
         volfine,cenfine,SDIM)

       do dir=1,num_materials*ngeom_recon
        mofdatafine(dir)=zero
       enddo
       volcell=zero
       do dir=1,SDIM
        cencell(dir)=zero
       enddo

       n_overlap=0

       call coarse_subelement_stencil(ifine,jfine,kfine,stenlo,stenhi, &
         bfact_coarse,bfact_fine)

       do ic=stenlo(1),stenhi(1)
        call intersect_weight_interp(ic,ifine, &
          bfact_coarse,bfact_fine,wt(1))
        if (wt(1).gt.zero) then
         do jc=stenlo(2),stenhi(2)
          call intersect_weight_interp(jc,jfine, &
           bfact_coarse,bfact_fine,wt(2))
          if (wt(2).gt.zero) then
           do kc=stenlo(3),stenhi(3)
            if (SDIM.eq.3) then
             call intersect_weight_interp(kc,kfine, &
              bfact_coarse,bfact_fine,wt(SDIM))
            endif
            if (wt(SDIM).gt.zero) then
             if (SDIM.eq.2) then
              testwt=wt(1)*wt(2)
             else if (SDIM.eq.3) then
              testwt=wt(1)*wt(2)*wt(SDIM)
             else
              print *,"dimension bust"
              stop
             endif

             if (testwt.gt.zero) then

              n_overlap=n_overlap+1

              do dir=1,num_materials*ngeom_recon
               mofdata(dir)=datamof(D_DECL(ic,jc,kc),dir)
              enddo
              call gridsten(xsten,problo,ic,jc,kc, &
               domlo,bfact_coarse,dxc,nhalf)
              do dir=1,SDIM
               xstengrid(-1,dir)=max(xsten(-1,dir),xstenfine(-1,dir))
               xstengrid(1,dir)=min(xsten(1,dir),xstenfine(1,dir))
               xstengrid(0,dir)=half*(xstengrid(-1,dir)+xstengrid(1,dir))
               if ((bfact_fine.eq.1).and.(bfact_coarse.eq.1)) then
                if ((abs(xstengrid(-1,dir)-xstenfine(-1,dir)).gt. &
                     VOFTOL*dxf(dir)).or. &
                    (abs(xstengrid(1,dir)-xstenfine(1,dir)).gt. &
                     VOFTOL*dxf(dir))) then
                 print *,"fine cell should be completely covered by coarse"
                 stop
                endif
               endif
              enddo ! dir

              call multi_get_volume_grid_simple( &
               tid_in, &
               EPS3, &
               tessellate, &  !=0
               bfact_coarse,dxc,xsten,nhalf, &
               mofdata, &
               xstengrid,nhalfgrid, &
               multi_volume,multi_cen, &
               geom_xtetlist(1,1,1,tid_in+1), &
               nmax, &
               nmax, &
               SDIM)

              do im=1,num_materials
               vofcomp_new=(im-1)*ngeom_recon+1
               mofdatafine(vofcomp_new)=mofdatafine(vofcomp_new)+ &
                multi_volume(im)
               do dir=1,SDIM
                mofdatafine(vofcomp_new+dir)= &
                 mofdatafine(vofcomp_new+dir)+ &
                 multi_cen(dir,im)*multi_volume(im)
               enddo

               if (is_rigid(im).eq.0) then
                volcell=volcell+multi_volume(im)
                do dir=1,SDIM
                 cencell(dir)=cencell(dir)+ &
                  multi_cen(dir,im)*multi_volume(im)
                enddo
               else if (is_rigid(im).eq.1) then
                ! do nothing
               else
                print *,"is_rigid(im) invalid"
                stop
               endif
              enddo ! im=1..num_materials
             else if (testwt.eq.zero) then
              ! do nothing
             else
              print *,"testwt invalid"
              stop
             endif
            else if (wt(SDIM).eq.zero) then
             ! do nothing
            else 
             print *,"wt(SDIM) invalid"
             stop
            endif

           enddo ! kc
          endif
         enddo ! jc
        endif
       enddo ! ic

       if ((bfact_fine.eq.1).and.(bfact_coarse.eq.1)) then
        if (n_overlap.ne.1) then
         print *,"n_overlap invalid"
         stop
        endif
       endif

       if (volcell.gt.zero) then
        !do nothing
       else
        print *,"volcell must be positive: ",volcell
        stop
       endif
       do dir=1,SDIM
        cencell(dir)=cencell(dir)/volcell
       enddo
       if (abs(volcell-volfine).le.EPS2*volfine) then
        !do nothing
       else
        print *,"volcell, volfine bad (multiextmofinterp):",volcell,volfine
        stop
       endif

       do im=1,num_materials
        vofcomp_old=(im-1)*ngeom_recon+1
        vofcomp_new=(im-1)*ngeom_recon+1
        voltemp=mofdatafine(vofcomp_new)/volcell
        fdatamof(D_DECL(ifine,jfine,kfine),vofcomp_old)=voltemp
        if (voltemp.gt.zero) then
         do dir=1,SDIM
          fdatamof(D_DECL(ifine,jfine,kfine),vofcomp_old+dir)= &
            mofdatafine(vofcomp_new+dir)/mofdatafine(vofcomp_new)- &
            cencell(dir)
         enddo
        else
         do dir=1,SDIM
          fdatamof(D_DECL(ifine,jfine,kfine),vofcomp_old+dir)=zero
         enddo 
        endif

        do dir=SDIM+1,2*SDIM+2
         fdatamof(D_DECL(ifine,jfine,kfine),vofcomp_old+dir)=zero
        enddo

       enddo ! im

       do im=1,num_materials
        vofcomp_old=(im-1)*ngeom_recon+1
        vofcomp_new=(im-1)*ngeom_recon+1
        mofdata(vofcomp_new)= &
          fdatamof(D_DECL(ifine,jfine,kfine),vofcomp_old)
        do dir=1,SDIM
         mofdata(vofcomp_new+dir)= &
           fdatamof(D_DECL(ifine,jfine,kfine),vofcomp_old+dir)
          !slope=0
         mofdata(vofcomp_new+SDIM+1+dir)=zero
        enddo

          ! order=0
        mofdata(vofcomp_new+SDIM+1)=zero

       enddo  ! im=1...num_materials

!      call gridsten(xstenfine,problo,ifine,jfine,kfine, &
!        domlo,bfact_fine,dxf,nhalf)
       call gridsten_level(xstenfine,ifine,jfine,kfine,levelf,nhalf)

       grid_index(1)=ifine
       grid_index(2)=jfine
       if (SDIM.eq.3) then
        grid_index(SDIM)=kfine
       endif

        ! sum F_fluid=1  sum F_solid<=1
       call make_vfrac_sum_ok_base( &
         cmofsten, &
         xstenfine,nhalf, &
         continuous_mof, &
         bfact_fine,dxf, &
         tessellate, & !=0
         mofdata,SDIM)

       do im=1,num_materials
        vofcomp_new=(im-1)*ngeom_recon+1
        vof_super(im)=mofdata(vofcomp_new)
       enddo

       call multimaterial_MOF( &
         tid_in, &
         bfact_fine,dxf,xstenfine,nhalf, &
         mof_verbose, &
         use_ls_data, &
         LS_stencil, &
         geom_xtetlist(1,1,1,tid_in+1), &
         geom_xtetlist(1,1,1,tid_in+1), &
         nmax, &
         nmax, &
         mofdata, & !intent(inout)
         vof_super, &
         multi_centroidA, &
         continuous_mof, & ! continuous_mof=STANDARD_MOF
         cmofsten, &
         grid_index, &
         grid_level, & !grid_level=-1
         SDIM)

       do dir=1,num_materials*ngeom_recon
        fdatamof(D_DECL(ifine,jfine,kfine),dir)=mofdata(dir)
       enddo

      enddo
      enddo
      enddo ! ifine,jfine,kfine
 
      return 
      end subroutine fort_multiextmofinterp


      subroutine fort_ext_burnvel_interp( &
       velflag, &
       time, &
       cburn,DIMS(cburn), &
       fburn,DIMS(fburn), &
       flo,fhi, &
       problo, &
       dxf,dxc, &
       nburning, &
       levelc,levelf, &
       bfact_coarse,bfact_fine) &
      bind(c,name='fort_ext_burnvel_interp')

      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: velflag
      integer, INTENT(in) :: levelc,levelf
      integer, INTENT(in) :: bfact_coarse,bfact_fine
      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: nburning
      integer, INTENT(in) :: DIMDEC(cburn)
      integer, INTENT(in) :: DIMDEC(fburn)
      integer, INTENT(in) :: flo(SDIM),fhi(SDIM)
      integer domlo(SDIM)
      real(amrex_real), INTENT(in) :: cburn(DIMV(cburn),nburning)
      real(amrex_real), INTENT(out) :: fburn(DIMV(fburn),nburning)
      real(amrex_real), INTENT(in) :: problo(SDIM)
      real(amrex_real), INTENT(in) :: dxf(SDIM),dxc(SDIM)
      integer growlo(3),growhi(3)
      integer stenlo(3),stenhi(3)
      real(amrex_real) wt(SDIM)

      integer ifine,jfine,kfine
      integer ic,jc,kc
      real(amrex_real) testwt

      integer dir
      integer iten
      integer im
      integer, parameter :: nhalf=3
      integer, parameter :: nhalfgrid=1
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xstenfine(-nhalf:nhalf,SDIM)
      real(amrex_real) xstengrid(-nhalfgrid:nhalfgrid,SDIM)
      integer n_overlap

      integer iflag,iflag_sign,hitflag
      real(amrex_real) burn_fine(nburning,-2:2)
      real(amrex_real) burn_coarse(nburning)
      integer n_burn_interface(num_interfaces,-2:2)
      integer n_burn_material(num_materials,0:2)
      integer burnstat
      integer bcomp
      real(amrex_real) rburnstat
      integer ncomp_per
      integer ncomp_expect
      integer drag_type,drag_im

      if (velflag.eq.0) then
       ncomp_per=EXTRAP_PER_TSAT ! interface temperature and mass fraction
       ncomp_expect=num_interfaces*(ncomp_per+1)
      else if (velflag.eq.1) then
       ncomp_per=EXTRAP_PER_BURNING
       ncomp_expect=num_interfaces*(ncomp_per+1)
      else if (velflag.eq.2) then
       ncomp_per=0
       ncomp_expect=N_DRAG
      else
       print *,"velflag invalid"
       stop
      endif

      if ((levelc.ne.levelf-1).or.(levelc.lt.0)) then
       print *,"levelc or levelf invalid"
       stop
      endif
      if (levelf.gt.fort_finest_level) then
       print *,"levelf invalid"
       stop
      endif
      if (bfact_coarse.lt.1) then
       print *,"bfact_coarse invalid"
       stop
      endif
      if (bfact_fine.lt.1) then
       print *,"bfact_fine invalid"
       stop
      endif
      if (bfact_fine.gt.bfact_coarse) then
       print *,"bfact_fine invalid"
       stop
      endif
      if (nburning.eq.ncomp_expect) then
       ! do nothing
      else
       print *,"nburning invalid"
       stop
      endif


      do dir=1,SDIM
       domlo(dir)=0
      enddo

      call growntilebox(flo,fhi,flo,fhi,growlo,growhi,0) 

      do kfine=growlo(3),growhi(3)
      do jfine=growlo(2),growhi(2)
      do ifine=growlo(1),growhi(1)

       do dir=1,nburning
        do iflag=-2,2
         burn_fine(dir,iflag)=zero
        enddo
       enddo

       n_overlap=0

       do iten=1,num_interfaces
        do iflag=-2,2
         n_burn_interface(iten,iflag)=0
        enddo
       enddo
       do im=1,num_materials
        do iflag=0,2
         n_burn_material(im,iflag)=0
        enddo
       enddo

       call coarse_subelement_stencil(ifine,jfine,kfine,stenlo,stenhi, &
         bfact_coarse,bfact_fine)

       do ic=stenlo(1),stenhi(1)
        call intersect_weight_interp(ic,ifine, &
          bfact_coarse,bfact_fine,wt(1))
        if (wt(1).gt.zero) then
         do jc=stenlo(2),stenhi(2)
          call intersect_weight_interp(jc,jfine, &
           bfact_coarse,bfact_fine,wt(2))
          if (wt(2).gt.zero) then
           do kc=stenlo(3),stenhi(3)
            if (SDIM.eq.3) then
             call intersect_weight_interp(kc,kfine, &
              bfact_coarse,bfact_fine,wt(SDIM))
            endif
            if (wt(SDIM).gt.zero) then
             if (SDIM.eq.2) then
              testwt=wt(1)*wt(2)
             else if (SDIM.eq.3) then
              testwt=wt(1)*wt(2)*wt(SDIM)
             else
              print *,"dimension bust"
              stop
             endif

             if (testwt.gt.zero) then

              n_overlap=n_overlap+1

              do dir=1,nburning
               burn_coarse(dir)=cburn(D_DECL(ic,jc,kc),dir)
              enddo

              call gridsten(xsten,problo,ic,jc,kc, &
               domlo,bfact_coarse,dxc,nhalf)
              call gridsten(xstenfine,problo,ifine,jfine,kfine, &
               domlo,bfact_fine,dxf,nhalf)
              do dir=1,SDIM
               xstengrid(-1,dir)=max(xsten(-1,dir),xstenfine(-1,dir))
               xstengrid(1,dir)=min(xsten(1,dir),xstenfine(1,dir))
               xstengrid(0,dir)=half*(xstengrid(-1,dir)+xstengrid(1,dir))
               if ((bfact_fine.eq.1).and.(bfact_coarse.eq.1)) then
                if ((abs(xstengrid(-1,dir)-xstenfine(-1,dir)).gt. &
                     VOFTOL*dxf(dir)).or. &
                    (abs(xstengrid(1,dir)-xstenfine(1,dir)).gt. &
                     VOFTOL*dxf(dir))) then
                 print *,"fine cell should be completely covered by coarse"
                 stop
                endif
               endif
              enddo ! dir=1..sdim

              if ((velflag.eq.0).or. &
                  (velflag.eq.1)) then

               do iten=1,num_interfaces
                rburnstat=burn_coarse(iten)
                burnstat=NINT(rburnstat)
                if ((burnstat.eq.0).and.(rburnstat.eq.zero)) then
                 ! do nothing (burn_fine and n_burn_interface already init to 0)
                else if (((burnstat.eq.1).and.(rburnstat.eq.one)).or.  &
                         ((burnstat.eq.2).and.(rburnstat.eq.two)).or.  &
                         ((burnstat.eq.-1).and.(rburnstat.eq.-one)).or. &
                         ((burnstat.eq.-2).and.(rburnstat.eq.-two))) then
                 n_burn_interface(iten,burnstat)= &
                   n_burn_interface(iten,burnstat)+1
                 do dir=1,ncomp_per
                  bcomp=num_interfaces+(iten-1)*ncomp_per+dir
                  burn_fine(bcomp,burnstat)= &
                   burn_fine(bcomp,burnstat)+burn_coarse(bcomp)
                 enddo ! dir
                else
                 print *,"burnstat or rburnstatus invalid"
                 print *, &
                  "num_materials,num_interfaces,iten,burnstat,rburnstat ",  &
                  num_materials,num_interfaces,iten,burnstat,rburnstat
                 do iflag=-2,2
                  print *,"iflag ",iflag
                  print *,"n_burn_interface(iten,iflag) ", &
                    n_burn_interface(iten,iflag)
                  print *,"n_overlap ",n_overlap
                  print *,"ifine,jfine,kfine ",ifine,jfine,kfine
                  print *,"ic,jc,kc ",ic,jc,kc
                  stop
                 enddo ! iflag=-2,2
                endif
               enddo ! iten=1..num_interfaces

              else if (velflag.eq.2) then

               do im=1,num_materials
                rburnstat=burn_coarse(DRAGCOMP_FLAG+im)
                burnstat=NINT(rburnstat)
                if ((burnstat.eq.0).and.(rburnstat.eq.zero)) then
                 ! do nothing (burn_fine and n_burn_material already init to 0)
                else if (((burnstat.eq.1).and.(rburnstat.eq.one)).or.  &
                         ((burnstat.eq.2).and.(rburnstat.eq.two))) then
                 n_burn_material(im,burnstat)= &
                    n_burn_material(im,burnstat)+1

                 do bcomp=0,nburning-1
                  drag_type=fort_drag_type(bcomp,drag_im)
                  if (drag_im+1.eq.im) then 
                   if ((drag_type.ge.0).and.(drag_type.lt.DRAG_TYPE_NEXT)) then
                    burn_fine(bcomp+1,burnstat)= &
                     burn_fine(bcomp+1,burnstat)+burn_coarse(bcomp+1)
                   else
                    print *,"drag_type invalid"
                    stop
                   endif
                  else if ((drag_im.ge.0).and.(drag_im.lt.num_materials)) then
                   if ((drag_type.ge.0).and.(drag_type.lt.DRAG_TYPE_NEXT)) then
                    ! do nothing
                   else
                    print *,"drag_type invalid"
                    stop
                   endif
                  else
                   print *,"drag_im invalid"
                   stop
                  endif
                 enddo ! bcomp=1..nburning
                else
                 print *,"burnstat or rburnstatus invalid"
                 print *, &
                   "num_materials,num_interfaces,im,burnstat,rburnstat",  &
                   num_materials,num_interfaces,im,burnstat,rburnstat
                 do iflag=0,2
                  print *,"iflag ",iflag
                  print *,"n_burn_material(im,iflag) ", &
                    n_burn_material(im,iflag)
                  print *,"n_overlap ",n_overlap
                  print *,"ifine,jfine,kfine ",ifine,jfine,kfine
                  print *,"ic,jc,kc ",ic,jc,kc
                  stop
                 enddo ! iflag=0,2
                endif
               enddo ! im=1..num_materials
              else
               print *,"velflag invalid"
               stop
              endif

             else if (testwt.eq.zero) then
              print *,"testwt cannot be 0"
              stop
             else
              print *,"testwt invalid"
              stop
             endif
            else if (wt(SDIM).eq.zero) then
             ! do nothing
            else 
             print *,"wt(SDIM) invalid"
             stop
            endif

           enddo ! kc
          endif
         enddo ! jc
        endif
       enddo ! ic

       if ((bfact_fine.eq.1).and.(bfact_coarse.eq.1)) then
        if (n_overlap.ne.1) then
         print *,"n_overlap invalid"
         stop
        endif
       endif

       if (n_overlap.ge.1) then

        if ((velflag.eq.0).or. &
            (velflag.eq.1)) then

         do iten=1,num_interfaces
          hitflag=0
          do iflag=1,2
           do iflag_sign=-1,1,2
            if (hitflag.eq.0) then
             if ((n_burn_interface(iten,iflag*iflag_sign).gt.0).and. &
                 (n_burn_interface(iten,iflag*iflag_sign).le.n_overlap)) then
              hitflag=1
              fburn(D_DECL(ifine,jfine,kfine),iten)=iflag*iflag_sign
              do dir=1,ncomp_per
               bcomp=num_interfaces+(iten-1)*ncomp_per+dir
               fburn(D_DECL(ifine,jfine,kfine),bcomp)= &
                burn_fine(bcomp,iflag*iflag_sign)/ &
                n_burn_interface(iten,iflag*iflag_sign)
              enddo
             else if (n_burn_interface(iten,iflag*iflag_sign).eq.0) then
              ! do nothing
             else
              print *,"n_burn_interface(iten,iflag*iflag_sign) invalid"
              stop
             endif
            else if (hitflag.eq.1) then
             ! do nothing
            else
             print *,"hitflag invalid"
             stop
            endif
           enddo ! iflag_sign=-1,1
          enddo ! iflag=1,2

          if (hitflag.eq.0) then
           fburn(D_DECL(ifine,jfine,kfine),iten)=zero
           do dir=1,ncomp_per
            bcomp=num_interfaces+(iten-1)*ncomp_per+dir
            fburn(D_DECL(ifine,jfine,kfine),bcomp)=zero
           enddo
          else if (hitflag.eq.1) then
           ! do nothing
          else     
           print *,"hitflag invalid"
           stop
          endif

         enddo ! iten=1..num_interfaces

        else if (velflag.eq.2) then

         do im=1,num_materials
          hitflag=0
          do iflag=1,2
           if (hitflag.eq.0) then
            if ((n_burn_material(im,iflag).gt.0).and. &
                (n_burn_material(im,iflag).le.n_overlap)) then
             hitflag=1
             fburn(D_DECL(ifine,jfine,kfine),DRAGCOMP_FLAG+im)=iflag
             do bcomp=0,nburning-1
              drag_type=fort_drag_type(bcomp,drag_im)
              if (drag_im+1.eq.im) then 
               if ((drag_type.ge.0).and. &
                   (drag_type.lt.DRAG_TYPE_NEXT).and. &
                   (drag_type.ne.DRAG_TYPE_FLAG)) then
                fburn(D_DECL(ifine,jfine,kfine),bcomp+1)= &
                  burn_fine(bcomp+1,iflag)/n_burn_material(im,iflag)
               else if (drag_type.eq.DRAG_TYPE_FLAG) then
                if (bcomp+1.eq.DRAGCOMP_FLAG+im) then
                 ! do nothing
                else
                 print *,"bcomp or get_drag_type invalid"
                 stop
                endif
               else
                print *,"drag_type invalid"
                stop
               endif
              else if ((drag_im.ge.0).and.(drag_im.lt.num_materials)) then
               if ((drag_type.ge.0).and.(drag_type.lt.DRAG_TYPE_NEXT)) then
                ! do nothing
               else
                print *,"drag_type invalid"
                stop
               endif
              else
               print *,"drag_im invalid"
               stop
              endif
             enddo ! bcomp=0..nburning-1
            else if (n_burn_material(im,iflag).eq.0) then
             ! do nothing
            else
             print *,"n_burn_material(im,iflag) invalid"
             stop
            endif
           else if (hitflag.eq.1) then
            ! do nothing
           else
            print *,"hitflag invalid"
            stop
           endif
          enddo ! iflag=1,2

          if (hitflag.eq.0) then

           fburn(D_DECL(ifine,jfine,kfine),DRAGCOMP_FLAG+im)=zero
           do bcomp=0,nburning-1
            drag_type=fort_drag_type(bcomp,drag_im)
            if (drag_im+1.eq.im) then 
             if ((drag_type.ge.0).and. &
                 (drag_type.lt.DRAG_TYPE_NEXT).and. &
                 (drag_type.ne.DRAG_TYPE_FLAG)) then
              fburn(D_DECL(ifine,jfine,kfine),bcomp+1)=zero
             else if (drag_type.eq.DRAG_TYPE_FLAG) then
              if (bcomp+1.eq.DRAGCOMP_FLAG+im) then
               ! do nothing
              else
               print *,"bcomp or get_drag_type invalid"
               stop
              endif
             else
              print *,"drag_type invalid"
              stop
             endif
            else if ((drag_im.ge.0).and.(drag_im.lt.num_materials)) then
             if ((drag_type.ge.0).and.(drag_type.lt.DRAG_TYPE_NEXT)) then
              ! do nothing
             else
              print *,"drag_type invalid"
              stop
             endif
            else
             print *,"drag_im invalid"
             stop
            endif
           enddo ! bcomp=0..nburning-1

          else if (hitflag.eq.1) then
           ! do nothing
          else     
           print *,"hitflag invalid"
           stop
          endif

         enddo ! iten=1..num_interfaces

        else
         print *,"velflag invalid"
         stop
        endif

       else
        print *,"n_overlap invalid"
        stop
       endif

      enddo
      enddo
      enddo ! ifine,jfine,kfine
 
      return 
      end subroutine fort_ext_burnvel_interp


      subroutine fort_pcinterp ( &
       grid_type, &
       zapflag, &
       crse_data, &
       DIMS(crse_data), &
       crse_bx_lo, & ! crse_bx=CoarseBox(fine_bx)
       crse_bx_hi, &
       fine_data, &
       DIMS(fine_data), &
       fblo,fbhi, & ! fine_bx=fine_region & fine.box()
       problo, &
       dxf,dxc, &
       nvar, &
       levelc,levelf, &
       bfact_coarse,bfact_fine) &
      bind(c,name='fort_pcinterp')

      use global_utility_module
      use probcommon_module

      implicit none

      integer, INTENT(in) :: grid_type  ! -1..5
      integer, INTENT(in) :: levelc,levelf
      integer, INTENT(in) :: bfact_coarse,bfact_fine
      integer, INTENT(in) :: zapflag
      integer, INTENT(in) :: crse_bx_lo(SDIM)
      integer, INTENT(in) :: crse_bx_hi(SDIM)
      integer clo(SDIM),chi(SDIM)
      integer, INTENT(in) :: DIMDEC(crse_data)
      integer, INTENT(in) :: DIMDEC(fine_data)
      integer, INTENT(in) :: fblo(SDIM), fbhi(SDIM)
      integer flo(SDIM),fhi(SDIM)
      integer, INTENT(in) :: nvar
      real(amrex_real), INTENT(in) :: crse_data(DIMV(crse_data),nvar)
      real(amrex_real), INTENT(out) :: fine_data(DIMV(fine_data),nvar)
      real(amrex_real), INTENT(in) :: problo(SDIM)
      real(amrex_real), INTENT(in) :: dxf(SDIM)
      real(amrex_real), INTENT(in) :: dxc(SDIM)
      integer stenlo(3),stenhi(3)
      integer stenlen(3)
      integer growlo(3),growhi(3)
      integer :: box_type(SDIM)

      integer ifine,jfine,kfine
      integer ic,jc,kc
      integer dir2
      integer n

      real(amrex_real) wt(SDIM)

      real(amrex_real) voltotal,volall
      real(amrex_real) fine_value(nvar)

      integer, parameter :: nhalf=1
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xstenND(-nhalf:nhalf,SDIM)
      real(amrex_real) xfine(SDIM)
      integer chi_loc(SDIM)

      if (bfact_coarse.lt.1) then
       print *,"bfact_coarse invalid"
       stop
      endif
      if (bfact_fine.lt.1) then
       print *,"bfact_fine invalid"
       stop
      endif
      if (bfact_fine.gt.bfact_coarse) then
       print *,"bfact_fine invalid"
       stop
      endif
      if ((levelc.ne.levelf-1).or.(levelc.lt.0)) then
       print *,"levelc or levelf invalid"
       stop
      endif
      if (levelf.gt.fort_finest_level) then
       print *,"levelf invalid"
       stop
      endif
      if ((nvar.lt.1).or.(nvar.gt.9999)) then
       print *,"nvar invalid in pcinterp"
       stop
      endif

      if ((grid_type.lt.-1).or.(grid_type.gt.5)) then
       print *,"grid_type invalid pcinterp"
       stop
      endif
      call grid_type_to_box_type(grid_type,box_type)

      do dir2=1,SDIM
       chi_loc(dir2)=bfact_coarse-1+box_type(dir2)
      enddo

      do dir2=1,SDIM
       clo(dir2)=crse_bx_lo(dir2)
       chi(dir2)=crse_bx_hi(dir2)
       flo(dir2)=fblo(dir2)
       fhi(dir2)=fbhi(dir2)
       if (box_type(dir2).eq.1) then
        chi(dir2)=chi(dir2)-1 
        fhi(dir2)=fhi(dir2)-1 
       else if (box_type(dir2).eq.0) then
        ! do nothing
       else
        print *,"box_type invalid"
        stop
       endif
      enddo ! dir2
      
      call growntileboxMAC(flo,fhi,flo,fhi,growlo,growhi,0,grid_type) 

      do kfine=growlo(3),growhi(3)
      do jfine=growlo(2),growhi(2)
      do ifine=growlo(1),growhi(1)

       call coarse_subelement_stencilMAC(ifine,jfine,kfine,stenlo,stenhi, &
         bfact_coarse,bfact_fine,grid_type)
       do dir2=1,SDIM
        stenlen(dir2)=stenhi(dir2)-stenlo(dir2)+1
        if (box_type(dir2).eq.0) then
         if (stenlen(dir2).ne.bfact_coarse) then
          print *,"stenlen invalid"
          stop
         endif
        else if (box_type(dir2).eq.1) then
         if ((stenlen(dir2).ne.bfact_coarse+1).and. &
             (stenlen(dir2).ne.1)) then
          print *,"stenlen invalid"
          stop
         endif
        else
         print *,"box_type invalid"
         stop
        endif
       enddo ! dir2=1..sdim

       call gridstenMAC_level(xsten,ifine,jfine,kfine,levelf,nhalf,grid_type)

       ic=stenlo(1)
       jc=stenlo(2)
       kc=stenlo(SDIM)
       call gridstenND_level(xstenND,ic,jc,kc,levelc,nhalf)
       do dir2=1,SDIM
        xfine(dir2)=xsten(0,dir2)-xstenND(0,dir2)
        if ((xfine(dir2).ge.-EPS3*dxc(dir2)).and. &
            (xfine(dir2).le.(EPS3+bfact_coarse)*dxc(dir2))) then
         !do nothing
        else
         print *,"xfine out of bounds"
         stop
        endif
       enddo ! dir2=1..sdim

       do n=1,nvar
        fine_value(n)=zero
       enddo
       voltotal=zero

       do ic=stenlo(1),stenhi(1)
        if (box_type(1).eq.1) then
         call intersect_weightMAC_interp(ic,ifine, &
           bfact_coarse,bfact_fine,wt(1))
        else if (box_type(1).eq.0) then
         call intersect_weight_interp(ic,ifine, &
           bfact_coarse,bfact_fine,wt(1))
        else
         print *,"box_type(1) invalid"
         stop
        endif
        if (wt(1).gt.zero) then
         do jc=stenlo(2),stenhi(2)
          if (box_type(2).eq.1) then
           call intersect_weightMAC_interp(jc,jfine, &
            bfact_coarse,bfact_fine,wt(2))
          else if (box_type(2).eq.0) then
           call intersect_weight_interp(jc,jfine, &
            bfact_coarse,bfact_fine,wt(2))
          else
           print *,"box_type(2) invalid"
           stop
          endif
          if (wt(2).gt.zero) then
           do kc=stenlo(3),stenhi(3)
            if (SDIM.eq.3) then
             if (box_type(SDIM).eq.1) then
              call intersect_weightMAC_interp(kc,kfine, &
               bfact_coarse,bfact_fine,wt(SDIM))
             else if (box_type(SDIM).eq.0) then
              call intersect_weight_interp(kc,kfine, &
               bfact_coarse,bfact_fine,wt(SDIM))
             else
              print *,"box_type(SDIM) invalid"
              stop
             endif
            endif
            if (wt(SDIM).gt.zero) then
             volall=wt(1)
             do dir2=2,SDIM
              volall=volall*wt(dir2)
             enddo
             do n=1,nvar
              if (zapflag.eq.0) then
               fine_value(n)=fine_value(n)+ &
                 volall*crse_data(D_DECL(ic,jc,kc),n)
              else if (zapflag.eq.1) then
               ! do nothing
              else
               print *,"zapflag invalid"
               stop
              endif
             enddo
             voltotal=voltotal+volall
            endif
           enddo ! kc
          endif
         enddo ! jc
        endif
       enddo ! ic

       if (voltotal.gt.zero) then
        do n=1,nvar
         fine_value(n)=fine_value(n)/voltotal
         fine_data(D_DECL(ifine,jfine,kfine),n)=fine_value(n) 
        enddo
       else
        print *,"voltotal invalid: ",voltotal
        stop
       endif

      enddo
      enddo
      enddo ! looping ifine,jfine,kfine

      end subroutine fort_pcinterp

      subroutine fort_refine_density_interp ( &
       crse_data, &
       DIMS(crse_data), &
       crse_bx_lo, & ! crse_bx=CoarseBox(fine_bx)
       crse_bx_hi, &
       fine_data, &
       DIMS(fine_data), &
       fblo,fbhi, & ! fine_bx=fine_region & fine.box()
       problo, &
       dxf,dxc, &
       nvar, &
       levelc,levelf, &
       bfact_coarse,bfact_fine) &
      bind(c,name='fort_refine_density_interp')

      use global_utility_module
      use probcommon_module

      implicit none

      integer, INTENT(in) :: levelc,levelf
      integer, INTENT(in) :: bfact_coarse,bfact_fine
      integer, INTENT(in) :: crse_bx_lo(SDIM)
      integer, INTENT(in) :: crse_bx_hi(SDIM)
      integer clo(SDIM),chi(SDIM)
      integer, INTENT(in) :: DIMDEC(crse_data)
      integer, INTENT(in) :: DIMDEC(fine_data)
      integer, INTENT(in) :: fblo(SDIM), fbhi(SDIM)
      integer flo(SDIM),fhi(SDIM)
      integer, INTENT(in) :: nvar
      real(amrex_real), INTENT(in) :: crse_data(DIMV(crse_data),nvar)
      real(amrex_real), INTENT(out) :: fine_data(DIMV(fine_data),nvar)
      real(amrex_real), INTENT(in) :: problo(SDIM)
      real(amrex_real), INTENT(in) :: dxf(SDIM)
      real(amrex_real), INTENT(in) :: dxc(SDIM)
      integer stenlo(3),stenhi(3)
      integer stenlen(3)
      integer growlo(3),growhi(3)

      integer ifine,jfine,kfine
      integer ifine2,jfine2,kfine2
      integer ic,jc,kc
      integer ic2,jc2,kc2
      integer dir2
      integer nfine,ncrse

      real(amrex_real) wt(SDIM)

      real(amrex_real) voltotal,volall
      real(amrex_real) fine_value

      integer, parameter :: nhalf=1
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xstenND(-nhalf:nhalf,SDIM)
      real(amrex_real) xfine(SDIM)
      integer chi_loc(SDIM)

      if (bfact_coarse.lt.1) then
       print *,"bfact_coarse invalid: ",bfact_coarse
       stop
      endif
      if (bfact_fine.lt.1) then
       print *,"bfact_fine invalid: ",bfact_fine
       stop
      endif
      if (bfact_fine.gt.bfact_coarse) then
       print *,"bfact_fine invalid: ",bfact_fine
       stop
      endif
      if ((levelc.ne.levelf-1).or.(levelc.lt.0)) then
       print *,"levelc or levelf invalid"
       stop
      endif
      if (levelf.gt.fort_finest_level) then
       print *,"levelf invalid"
       stop
      endif
      if (nvar.ne.4*(AMREX_SPACEDIM-1)) then
       print *,"nvar invalid in fort_refine_density_interp"
       stop
      endif

      do dir2=1,SDIM
       chi_loc(dir2)=bfact_coarse-1
      enddo

      do dir2=1,SDIM
       clo(dir2)=crse_bx_lo(dir2)
       chi(dir2)=crse_bx_hi(dir2)
       flo(dir2)=fblo(dir2)
       fhi(dir2)=fbhi(dir2)
      enddo ! dir2
      
      call growntilebox(flo,fhi,flo,fhi,growlo,growhi,0) 

      do kfine=growlo(3),growhi(3)
      do jfine=growlo(2),growhi(2)
      do ifine=growlo(1),growhi(1)

         !coarse_subelement_stencil is declared in GLOBALUTIL.F90
       call coarse_subelement_stencil(ifine,jfine,kfine,stenlo,stenhi, &
         bfact_coarse,bfact_fine)
       do dir2=1,SDIM
        stenlen(dir2)=stenhi(dir2)-stenlo(dir2)+1
        if (stenlen(dir2).ne.bfact_coarse) then
         print *,"stenlen invalid"
         stop
        endif
       enddo ! dir2=1..sdim

       call gridsten_level(xsten,ifine,jfine,kfine,levelf,nhalf)

       ic=stenlo(1)
       jc=stenlo(2)
       kc=stenlo(SDIM)
       call gridstenND_level(xstenND,ic,jc,kc,levelc,nhalf)
       do dir2=1,SDIM
        xfine(dir2)=xsten(0,dir2)-xstenND(0,dir2)
        if ((xfine(dir2).ge.-EPS3*dxc(dir2)).and. &
            (xfine(dir2).le.(EPS3+bfact_coarse)*dxc(dir2))) then
         !do nothing
        else
         print *,"xfine out of bounds: ",dir2,xfine(dir2)
         stop
        endif
       enddo ! dir2=1..sdim

       nfine=0
       kfine2=0
#if (AMREX_SPACEDIM==3)
       do kfine2=0,1
#endif
       do jfine2=0,1
       do ifine2=0,1
        nfine=nfine+1

        if (nfine.eq.4*kfine2+2*jfine2+ifine2+1) then
         !do nothing
        else
         print *,"nfine invalid: ",nfine
         stop
        endif

        fine_value=zero
        voltotal=zero

        do ic=stenlo(1),stenhi(1)
         do ic2=0,1
          call intersect_weight_interp_refine( &
           ic,ic2,ifine,ifine2,  &
           bfact_coarse,bfact_fine,wt(1))

          if (wt(1).gt.zero) then
           do jc=stenlo(2),stenhi(2)
            do jc2=0,1
             call intersect_weight_interp_refine( &
              jc,jc2,jfine,jfine2, &
              bfact_coarse,bfact_fine,wt(2))

             if (wt(2).gt.zero) then
              do kc=stenlo(3),stenhi(3)
               kc2=0
#if (AMREX_SPACEDIM==3)
               do kc2=0,1
#endif
                if (SDIM.eq.3) then
                 call intersect_weight_interp_refine( &
                  kc,kc2,kfine,kfine2, &
                  bfact_coarse,bfact_fine,wt(SDIM))
                endif

                if (wt(SDIM).gt.zero) then
                 volall=wt(1)
                 do dir2=2,SDIM
                  volall=volall*wt(dir2)
                 enddo
                 ncrse=kc2*4+jc2*2+ic2+1
                 fine_value=fine_value+ &
                  volall*crse_data(D_DECL(ic,jc,kc),ncrse)
                 voltotal=voltotal+volall
                endif
#if (AMREX_SPACEDIM==3)
               enddo ! kc2
#endif
              enddo ! kc
             endif
            enddo ! jc2
           enddo ! jc
          endif
         enddo ! ic2
        enddo ! ic

        if (voltotal.gt.zero) then
         fine_value=fine_value/voltotal
         fine_data(D_DECL(ifine,jfine,kfine),nfine)=fine_value
        else
         print *,"voltotal invalid: ",voltotal
         stop
        endif
       enddo !ifine2
       enddo !jfine2
#if (AMREX_SPACEDIM==3)
       enddo !kfine2
#endif

      enddo
      enddo
      enddo ! looping ifine,jfine,kfine

      end subroutine fort_refine_density_interp

      ! enable_spectral:
      ! 0 - low order
      ! 1 - space/time spectral
      subroutine fort_seminterp ( &
       enable_spectral, &
       dxc,dxf, &
       crse,DIMS(crse), &
       fine,DIMS(fine), &
       fblo,fbhi, &
       nvar, &
       levelc,levelf, &
       bfact_coarse,bfact_fine) &
      bind(c,name='fort_seminterp')

      use global_utility_module
      use probcommon_module

      implicit none

      integer, INTENT(in) :: enable_spectral
      integer, INTENT(in) :: levelc,levelf
      integer, INTENT(in) :: bfact_coarse,bfact_fine
      integer, INTENT(in) :: DIMDEC(crse)
      integer, INTENT(in) :: DIMDEC(fine)
      integer, INTENT(in) :: fblo(SDIM), fbhi(SDIM)
      integer, INTENT(in) :: nvar
      real(amrex_real), INTENT(in) :: dxc(SDIM)
      real(amrex_real), INTENT(in) :: dxf(SDIM)
      real(amrex_real), INTENT(in) :: crse(DIMV(crse), nvar)
      real(amrex_real), INTENT(out) :: fine(DIMV(fine), nvar)
      integer growlo(3),growhi(3)
      integer stenlo(3),stenhi(3),stenlen(3)

      real(amrex_real) wt(SDIM)

      integer ic,jc,kc
      integer ifine,jfine,kfine
      integer ilocal,jlocal,klocal
      integer n,dir

      real(amrex_real) fine_value(nvar)
      real(amrex_real) voltotal,volall
      integer grid_type

      integer, parameter :: nhalf=1
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xstenND(-nhalf:nhalf,SDIM)
      real(amrex_real) xfine(SDIM)
      real(amrex_real), dimension(D_DECL(:,:,:),:),allocatable :: fcoarse
      integer chi_loc(SDIM)
      integer do_spectral_interp

      if (bfact_coarse.lt.1) then
       print *,"bfact_coarse invalid"
       stop
      endif
      if (bfact_fine.lt.1) then
       print *,"bfact_fine invalid"
       stop
      endif
      if (bfact_fine.gt.bfact_coarse) then
       print *,"bfact_fine invalid"
       stop
      endif
      if ((levelc.ne.levelf-1).or.(levelc.lt.0)) then
       print *,"levelc or levelf invalid"
       stop
      endif
      if ((nvar.lt.1).or.(nvar.gt.9999)) then
       print *,"nvar invalid"
       stop
      endif
      if (levelf.gt.fort_finest_level) then
       print *,"levelf invalid"
       stop
      endif

      grid_type=-1  ! ggg  (Gauss in all directions)

      do dir=1,SDIM
       chi_loc(dir)=bfact_coarse-1
      enddo
      allocate(fcoarse(D_DECL(0:chi_loc(1),0:chi_loc(2),0:chi_loc(3)),nvar))
 
      call growntilebox(fblo,fbhi,fblo,fbhi,growlo,growhi,0) 

      do kfine=growlo(3),growhi(3)
      do jfine=growlo(2),growhi(2)
      do ifine=growlo(1),growhi(1)

       call coarse_subelement_stencil(ifine,jfine,kfine,stenlo,stenhi, &
        bfact_coarse,bfact_fine)
       do dir=1,SDIM
        stenlen(dir)=stenhi(dir)-stenlo(dir)+1
        if (stenlen(dir).ne.bfact_coarse) then
         print *,"stenlen invalid"
         stop
        endif
       enddo ! dir

       call gridsten_level(xsten,ifine,jfine,kfine,levelf,nhalf)

       ic=stenlo(1) 
       jc=stenlo(2) 
       kc=stenlo(SDIM) 
       call gridstenND_level(xstenND,ic,jc,kc,levelc,nhalf)

       do dir=1,SDIM
        xfine(dir)=xsten(0,dir)-xstenND(0,dir)
        if ((xfine(dir).ge.-EPS3*dxc(dir)).and. &
            (xfine(dir).le.(EPS3+bfact_coarse)*dxc(dir))) then
         ! do nothing
        else
         print *,"xfine out of bounds"
         stop
        endif
       enddo ! dir

       do_spectral_interp=1

       if (bfact_coarse.eq.1) then
        do_spectral_interp=0
       endif

       if (enable_spectral.eq.0) then
        do_spectral_interp=0
       else if (enable_spectral.eq.1) then
        ! do nothing
       else
        print *,"enable_spectral invalid sem interp"
        stop
       endif

       if (do_spectral_interp.eq.1) then

        do kc=stenlo(3),stenhi(3)
        do jc=stenlo(2),stenhi(2)
        do ic=stenlo(1),stenhi(1)
         ilocal=ic-stenlo(1)
         jlocal=jc-stenlo(2)
         klocal=kc-stenlo(3)
         do n=1,nvar
          fcoarse(D_DECL(ilocal,jlocal,klocal),n)=crse(D_DECL(ic,jc,kc),n)
         enddo
        enddo
        enddo
        enddo

        call SEM_INTERP_ELEMENT( &
         nvar,bfact_coarse,grid_type, &
         chi_loc,dxc,xfine,fcoarse,fine_value)

        voltotal=one

       else if (do_spectral_interp.eq.0) then

        do n=1,nvar
         fine_value(n)=zero
        enddo
        voltotal=zero

        do ic=stenlo(1),stenhi(1)
         call intersect_weight_interp(ic,ifine, &
          bfact_coarse,bfact_fine,wt(1))
         if (wt(1).gt.zero) then
          do jc=stenlo(2),stenhi(2)
           call intersect_weight_interp(jc,jfine, &
            bfact_coarse,bfact_fine,wt(2))
           if (wt(2).gt.zero) then
            do kc=stenlo(3),stenhi(3)
             if (SDIM.eq.3) then
              call intersect_weight_interp(kc,kfine, &
               bfact_coarse,bfact_fine,wt(SDIM))
             endif
             if (wt(SDIM).gt.zero) then
              volall=wt(1)
              do dir=2,SDIM
               volall=volall*wt(dir)
              enddo
              do n=1,nvar
               fine_value(n)=fine_value(n)+volall*crse(D_DECL(ic,jc,kc),n)
              enddo
              voltotal=voltotal+volall
             endif
            enddo ! kc
           endif
          enddo ! jc
         endif
        enddo ! ic

       else
        print *,"do_spectral_interp invalid"
        stop
       endif

       if (voltotal.gt.zero) then
        do n=1,nvar
         fine_value(n)=fine_value(n)/voltotal
         fine(D_DECL(ifine,jfine,kfine),n)=fine_value(n)
        enddo
       else
        print *,"voltotal invalid: ",voltotal
        stop
       endif

      enddo
      enddo
      enddo ! looping ifine,jfine,kfine

      deallocate(fcoarse)

      end subroutine fort_seminterp

      subroutine fort_maskinterppc ( &
       crse,DIMS(crse), &
       fine,DIMS(fine), &
       fblo,fbhi, &
       nvar, &
       levelc,levelf, &
       bfact_coarse,bfact_fine) &
      bind(c,name='fort_maskinterppc')

      use global_utility_module
      use probcommon_module

      implicit none

      integer, INTENT(in) :: levelc,levelf
      integer, INTENT(in) :: bfact_coarse,bfact_fine
      integer dir
      integer, INTENT(in) :: DIMDEC(crse)
      integer cdlo(SDIM), cdhi(SDIM)
      integer, INTENT(in) :: DIMDEC(fine)
      integer, INTENT(in) :: fblo(SDIM), fbhi(SDIM)
      integer, INTENT(in) :: nvar
      real(amrex_real), INTENT(in) :: crse(DIMV(crse), nvar)
      real(amrex_real), INTENT(out) :: fine(DIMV(fine), nvar)
      integer stenlo(3),stenhi(3)
      integer growlo(3),growhi(3)

      integer ifine,jfine,kfine,ic,jc,kc,n
      real(amrex_real) wt(SDIM)

      real(amrex_real) voltotal,volall
      real(amrex_real) fine_value
      integer first_hit,m1,m2

      if (bfact_coarse.lt.1) then
       print *,"bfact_coarse invalid"
       stop
      endif
      if (bfact_fine.lt.1) then
       print *,"bfact_fine invalid"
       stop
      endif
      if (bfact_fine.gt.bfact_coarse) then
       print *,"bfact_fine invalid"
       stop
      endif
      if ((levelc.ne.levelf-1).or.(levelc.lt.0)) then
       print *,"levelc or levelf invalid"
       stop
      endif
      if (levelf.gt.fort_finest_level) then
       print *,"levelf invalid"
       stop
      endif

      cdlo(1)=ARG_L1(crse)
      cdlo(2)=ARG_L2(crse)
      cdhi(1)=ARG_H1(crse)
      cdhi(2)=ARG_H2(crse)
#if (AMREX_SPACEDIM==3)
      cdlo(SDIM)=ARG_L3(crse)
      cdhi(SDIM)=ARG_H3(crse)
#endif

      call growntilebox(fblo,fbhi,fblo,fbhi,growlo,growhi,0) 

      if (nvar.ne.1) then
       print *,"nvar invalid"
       stop
      endif
      n=1

      do kfine=growlo(3),growhi(3)
      do jfine=growlo(2),growhi(2)
      do ifine=growlo(1),growhi(1)

       fine_value=zero
       first_hit=0
       voltotal=zero

       call coarse_subelement_stencil(ifine,jfine,kfine,stenlo,stenhi, &
         bfact_coarse,bfact_fine)
       do ic=stenlo(1),stenhi(1)
        call intersect_weight_interp(ic,ifine, &
         bfact_coarse,bfact_fine,wt(1))
        if (wt(1).gt.zero) then
         do jc=stenlo(2),stenhi(2)
          call intersect_weight_interp(jc,jfine, &
           bfact_coarse,bfact_fine,wt(2))
          if (wt(2).gt.zero) then
           do kc=stenlo(3),stenhi(3)
            if (SDIM.eq.3) then
             call intersect_weight_interp(kc,kfine, &
              bfact_coarse,bfact_fine,wt(SDIM))
            endif
            if (wt(SDIM).gt.zero) then
             volall=wt(1)
             do dir=2,SDIM
              volall=volall*wt(dir)
             enddo

             dir=1
             if ((ic.lt.cdlo(dir)).or.(ic.gt.cdhi(dir))) then
              print *,"ic out of range ic=",ic
              print *,"cdlo,cdhi ",cdlo(dir),cdhi(dir)
              print *,"bfact_fine,bfact_coarse ",bfact_fine,bfact_coarse
              stop
             endif
             dir=2
             if ((jc.lt.cdlo(dir)).or.(jc.gt.cdhi(dir))) then
              print *,"jc out of range jc=",jc
              print *,"cdlo,cdhi ",cdlo(dir),cdhi(dir)
              print *,"bfact_fine,bfact_coarse ",bfact_fine,bfact_coarse
              stop
             endif
             if (SDIM.eq.3) then
              dir=SDIM
              if ((kc.lt.cdlo(dir)).or.(kc.gt.cdhi(dir))) then
               print *,"kc out of range"
               print *,"cdlo,cdhi ",cdlo(dir),cdhi(dir)
               print *,"bfact_fine,bfact_coarse ",bfact_fine,bfact_coarse
               stop
              endif
             endif

             if (volall.gt.zero) then
              if (first_hit.eq.0) then
               fine_value=crse(D_DECL(ic,jc,kc),n)
               first_hit=1
              else if (first_hit.eq.1) then
               if (fine_value.eq.zero) then
                ! do nothing
               else
                m1=NINT(fine_value)
                m2=NINT(crse(D_DECL(ic,jc,kc),n))
                if ((m1.le.0).or.(m2.lt.0)) then
                 print *,"m1 or m2 are invalid"
                 stop
                endif
                if (m1.ne.m2) then
                 fine_value=zero
                endif 
               endif
              else
               print *,"first_hit invalid"
               stop
              endif
              voltotal=voltotal+volall
             else if (volall.eq.zero) then
              ! do nothing
             else if (volall.lt.zero) then
              print *,"volall invalid"
              stop
             endif
              
            endif ! wt(sdim)>0
           enddo ! kc
          endif
         enddo ! jc
        endif
       enddo ! ic

       if (voltotal.gt.zero) then
        m1=NINT(fine_value)
        if (m1.ge.0) then
         fine(D_DECL(ifine,jfine,kfine),n)=fine_value
        else
         print *,"m1 invalid"
         stop
        endif
       else
        print *,"voltotal invalid: ",voltotal 
        stop
       endif

      enddo
      enddo
      enddo ! looping ifine,jfine,kfine

      end subroutine fort_maskinterppc


      ! cloMAC,chiMAC and floMAC,fhiMAC are face centered boxes
      ! data and finedata are face data
      ! enable_spectral:
      ! 0 - low order
      ! 1 - space/time spectral
      subroutine fort_edgeinterp( &
       enable_spectral, &
       grid_type, & ! -1...5
       cdata, &
       DIMS(cdata), &
       cloMAC,chiMAC, &
       finedata, &
       DIMS(fdata), &
       floMAC,fhiMAC, &
       problo, &
       dxf,dxc, &
       nvar, &
       levelc,levelf, &
       bfact_coarse,bfact_fine) &
      bind(c,name='fort_edgeinterp')

      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: enable_spectral
      integer, INTENT(in) :: levelc,levelf
      integer, INTENT(in) :: bfact_coarse,bfact_fine
      integer, INTENT(in) :: nvar
      integer, INTENT(in) :: grid_type ! -1..5
      integer, INTENT(in) :: DIMDEC(cdata)
      integer, INTENT(in) :: DIMDEC(fdata)
      integer, INTENT(in) :: cloMAC(SDIM),chiMAC(SDIM)
      integer clo(SDIM),chi(SDIM)
      integer, INTENT(in) :: floMAC(SDIM),fhiMAC(SDIM)
      integer flo(SDIM),fhi(SDIM)
      real(amrex_real), INTENT(in) :: cdata(DIMV(cdata),nvar)
      real(amrex_real), INTENT(out) :: finedata(DIMV(fdata),nvar)
      real(amrex_real), INTENT(in) :: problo(SDIM)
      real(amrex_real), INTENT(in) :: dxf(SDIM)
      real(amrex_real), INTENT(in) :: dxc(SDIM)
      integer growlo(3),growhi(3)
      integer stenlo(3),stenhi(3)
      integer stenlen(3)
      integer :: box_type(SDIM)

      real(amrex_real) wt(SDIM)

      integer dir2
      integer ic,jc,kc
      integer ifine,jfine,kfine
      integer ilocal,jlocal,klocal
      integer n

      real(amrex_real) fine_value(nvar)
      real(amrex_real) voltotal,volall

      integer, parameter :: nhalf=1
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xstenND(-nhalf:nhalf,SDIM)
      real(amrex_real) xfine(SDIM)
      real(amrex_real), dimension(D_DECL(:,:,:),:),allocatable :: fcoarse
      integer chi_loc(SDIM)
      integer khi
      integer do_spectral_interp

      if (bfact_coarse.lt.1) then
       print *,"bfact_coarse invalid"
       stop
      endif
      if (bfact_fine.lt.1) then
       print *,"bfact_fine invalid"
       stop
      endif
      if (bfact_fine.gt.bfact_coarse) then
       print *,"bfact_fine invalid"
       stop
      endif
      if ((levelc.ne.levelf-1).or.(levelc.lt.0)) then
       print *,"levelc or levelf invalid"
       stop
      endif
      if (levelf.gt.fort_finest_level) then
       print *,"levelf invalid"
       stop
      endif
      if ((nvar.lt.1).or.(nvar.gt.9999)) then
       print *,"nvar invalid in edge interp"
       stop
      endif

      if ((grid_type.lt.0).or.(grid_type.ge.SDIM)) then
       print *,"grid_type invalid edgeinterp"
       stop
      endif

      call grid_type_to_box_type(grid_type,box_type)

      do dir2=1,SDIM
       chi_loc(dir2)=bfact_coarse-1+box_type(dir2)
      enddo

      allocate(fcoarse(D_DECL(0:chi_loc(1),0:chi_loc(2),0:chi_loc(3)),nvar))

      do dir2=1,SDIM
       clo(dir2)=cloMAC(dir2)
       chi(dir2)=chiMAC(dir2)
       flo(dir2)=floMAC(dir2)
       fhi(dir2)=fhiMAC(dir2)
       if (box_type(dir2).eq.1) then
        chi(dir2)=chi(dir2)-1 
        fhi(dir2)=fhi(dir2)-1 
       else if (box_type(dir2).eq.0) then
        ! do nothing
       else
        print *,"box_type invalid"
        stop
       endif
      enddo ! dir2
      
      call growntileboxMAC(flo,fhi,flo,fhi,growlo,growhi,0,grid_type) 

      do kfine=growlo(3),growhi(3)
      do jfine=growlo(2),growhi(2)
      do ifine=growlo(1),growhi(1)

       call coarse_subelement_stencilMAC(ifine,jfine,kfine,stenlo,stenhi, &
         bfact_coarse,bfact_fine,grid_type)
       do dir2=1,SDIM
        stenlen(dir2)=stenhi(dir2)-stenlo(dir2)+1
        if (box_type(dir2).eq.0) then
         if (stenlen(dir2).ne.bfact_coarse) then
          print *,"stenlen invalid"
          stop
         endif
        else if (box_type(dir2).eq.1) then
         if ((stenlen(dir2).ne.bfact_coarse+1).and. &
             (stenlen(dir2).ne.1)) then
          print *,"stenlen invalid"
          stop
         endif
        else
         print *,"box_type invalid"
         stop
        endif
       enddo ! dir2=1..sdim

       call gridstenMAC_level(xsten,ifine,jfine,kfine,levelf,nhalf,grid_type)

       ic=stenlo(1)
       jc=stenlo(2)
       kc=stenlo(SDIM)
       call gridstenND_level(xstenND,ic,jc,kc,levelc,nhalf)
       do dir2=1,SDIM
        xfine(dir2)=xsten(0,dir2)-xstenND(0,dir2)
        if ((xfine(dir2).ge.-EPS3*dxc(dir2)).and. &
            (xfine(dir2).le.(EPS3+bfact_coarse)*dxc(dir2))) then
         !do nothing
        else
         print *,"xfine out of bounds"
         stop
        endif
       enddo ! dir2=1..sdim

       do_spectral_interp=1

       if (bfact_coarse.eq.1) then
        do_spectral_interp=0
       endif

       if (enable_spectral.eq.0) then
        do_spectral_interp=0
       else if (enable_spectral.eq.1) then
        ! do nothing
       else
        print *,"enable_spectral invalid edge interp"
        stop
       endif

       if (do_spectral_interp.eq.1) then

        if (SDIM.eq.2) then
         khi=0
        else if (SDIM.eq.3) then
         khi=chi_loc(SDIM)
        else
         print *,"dimension bust"
         stop
        endif
        do klocal=0,khi
        do jlocal=0,chi_loc(2)
        do ilocal=0,chi_loc(1)
         do n=1,nvar
          fcoarse(D_DECL(ilocal,jlocal,klocal),n)=zero
         enddo
        enddo
        enddo
        enddo

        do kc=stenlo(3),stenhi(3)
        do jc=stenlo(2),stenhi(2)
        do ic=stenlo(1),stenhi(1)
         ilocal=ic-stenlo(1)
         jlocal=jc-stenlo(2)
         klocal=kc-stenlo(3)
         do n=1,nvar
          fcoarse(D_DECL(ilocal,jlocal,klocal),n)=cdata(D_DECL(ic,jc,kc),n)
         enddo
        enddo
        enddo
        enddo

        call SEM_INTERP_ELEMENT( &
         nvar,bfact_coarse,grid_type, &
         chi_loc,dxc,xfine,fcoarse,fine_value)

        voltotal=one

       else if (do_spectral_interp.eq.0) then

        do n=1,nvar
         fine_value(n)=zero
        enddo
        voltotal=zero

        do ic=stenlo(1),stenhi(1)
         if (box_type(1).eq.1) then
          call intersect_weightMAC_interp(ic,ifine, &
           bfact_coarse,bfact_fine,wt(1))
         else if (box_type(1).eq.0) then
          call intersect_weight_interp(ic,ifine, &
           bfact_coarse,bfact_fine,wt(1))
         else
          print *,"box_type(1) invalid"
          stop
         endif
         if (wt(1).gt.zero) then
          do jc=stenlo(2),stenhi(2)
           if (box_type(2).eq.1) then
            call intersect_weightMAC_interp(jc,jfine, &
             bfact_coarse,bfact_fine,wt(2))
           else if (box_type(2).eq.0) then
            call intersect_weight_interp(jc,jfine, &
             bfact_coarse,bfact_fine,wt(2))
           else
            print *,"box_type(2) invalid"
            stop
           endif
           if (wt(2).gt.zero) then
            do kc=stenlo(3),stenhi(3)
             if (SDIM.eq.3) then
              if (box_type(SDIM).eq.1) then
               call intersect_weightMAC_interp(kc,kfine, &
                bfact_coarse,bfact_fine,wt(SDIM))
              else if (box_type(SDIM).eq.0) then
               call intersect_weight_interp(kc,kfine, &
                bfact_coarse,bfact_fine,wt(SDIM))
              else
               print *,"box_type(SDIM) invalid"
               stop
              endif
             endif
             if (wt(SDIM).gt.zero) then
              volall=wt(1)
              do dir2=2,SDIM
               volall=volall*wt(dir2)
              enddo
              do n=1,nvar
               fine_value(n)=fine_value(n)+volall*cdata(D_DECL(ic,jc,kc),n)
              enddo
              voltotal=voltotal+volall
             endif
            enddo ! kc
           endif
          enddo ! jc
         endif
        enddo ! ic

       else
        print *,"do_spectral_interp invalid"
        stop
       endif

       if (voltotal.gt.zero) then
        do n=1,nvar
         fine_value(n)=fine_value(n)/voltotal
         finedata(D_DECL(ifine,jfine,kfine),n)=fine_value(n) 
        enddo
       else
        print *,"voltotal invalid: ",voltotal
        stop
       endif

      enddo
      enddo
      enddo ! looping ifine,jfine,kfine

      deallocate(fcoarse)

      return 
      end subroutine fort_edgeinterp

      end module interp_module

