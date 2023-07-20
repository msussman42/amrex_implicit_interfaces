#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else  
print *,"dimension bust"
stop
#endif

#include "AMReX_FORT_INTEGER.H"
#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"
#include "EXTRAP_COMP.H"

#include "CTMLFSI_F.H"

      module CTML_module
      use probf90_module

      contains

      subroutine CTML_INIT_SOLID(&
       dx_max_level,&
       prob_lo,&
       prob_hi, &
       io_proc, &
       n_fib_bodies,&
       max_n_fib_nodes)

      use probcommon_module
      use dummy_module
      use probf90_module
      use flow_parameters 
      use grid_arrays

      IMPLICIT NONE
      REAL_T, INTENT(in) :: dx_max_level(SDIM)
      REAL_T, INTENT(in) :: prob_lo(SDIM)
      REAL_T, INTENT(in) :: prob_hi(SDIM)
      INTEGER_T, INTENT(in) :: io_proc
      INTEGER_T, INTENT(out) :: n_fib_bodies
      INTEGER_T, INTENT(out) :: max_n_fib_nodes
      
      INTEGER_T n_Read_in,i,dir
      logical the_boss

      if (1.eq.1) then
       print *,"calling CTML_INIT_SOLID"
      endif

      do dir=1,SDIM
       if (dx_max_level(dir).gt.zero) then
        !do nothing
       else
        print *,"dx_max_level(dir) invalid: ",dir,dx_max_level(dir)
        stop
       endif
      enddo

      !! call CTML_InitVicarVariables(dx_max_level,prob_lo,prob_hi)
      !! ******** x direction **********
      nxc_GLBL = (prob_hi(1) - prob_lo(1)) / dx_max_level(1)
      nx_GLBL = nxc_GLBL+1

      allocate( x(0:nx_GLBL+1))
      allocate(dx(0:nx_GLBL+1))

      do i=0,nxc_GLBL+1
       dx(i)=dx_max_level(1)
      enddo
       
      do i=0,nx_GLBL+1
        x(i)=prob_lo(1) + ((i-1)*dx_max_level(1))
      enddo

!! ******** y direction **********
      nyc_GLBL = (prob_hi(2) - prob_lo(2)) / dx_max_level(2)
      ny_GLBL = nyc_GLBL+1

      allocate( y(0:ny_GLBL+1))
      allocate(dy(0:ny_GLBL+1))

      do i=0,nyc_GLBL+1
       dy(i)=dx_max_level(2)
      enddo

      do i=0,ny_GLBL+1
       y(i)=prob_lo(2) + ((i-1)*dx_max_level(2))
      enddo

!! ******** z direction **********
!! if it is a 2D problem, fake it as a 3D problem!

#if (AMREX_SPACEDIM==2)
      nzc = 2
      nz = nzc+1
      allocate( z(0:nz+1))
      allocate(dz(0:nz+1))

      do i=0,nzc+1
       dz(i)=sqrt(dx_max_level(1)*dx_max_level(2)) 
      enddo

      do i=0,nz+1
       z(i)= (i-1)*sqrt(dx_max_level(1)*dx_max_level(2))
      enddo


#elif (AMREX_SPACEDIM==3)
      nzc = (prob_hi(3) - prob_lo(3)) / dx_max_level(3)
      nz = nzc+1

      allocate( z(0:nz+1))
      allocate(dz(0:nz+1))

      do i=0,nzc+1
       dz(i)=dx_max_level(3)
      enddo

      do i=0,nz+1
       z(i)=problo(3) + ((i-1)*dx_max_level(3))
      enddo

#endif


      n_Read_in = 0
      the_boss = .true.

       ! :set ignorecase  (ignore case for vi searches)
       ! :set ic  (ignore case for vi searches)
       ! grep -i  (ignore case)      
       ! amrex_implicit_interfaces/Vicar3D/UTIL_BOUNDARY_FORCE_FSI.F90
       ! amrex_implicit_interfaces/Vicar3D/distFSI/grid_def
       ! amrex_implicit_interfaces/Vicar3D/distFSI/initialize_ibm.F
       ! amrex_implicit_interfaces/Vicar3D/distFSI/create_grid2.F 
       ! open (9997,file='./inputdistibm.dat',status="unknown")
      call init_membrane_solver(SDIM,n_Read_in,the_boss)

      n_fib_bodies = nrIBM_fib
      max_n_fib_nodes = nIBM_fib

      vel_fib = zero
      force_fib = zero

      if (1.eq.1) then
       print *,"done with CTML_INIT_SOLID"
      endif

      return
      end subroutine CTML_INIT_SOLID

      subroutine CTML_GET_FIB_NODE_COUNT(&
       n_fib_bodies,&
       n_fib_nodes)

      use dummy_module
  
      IMPLICIT NONE
      INTEGER_T, INTENT(in) :: n_fib_bodies
      INTEGER_T, INTENT(out) :: n_fib_nodes(n_fib_bodies)

      integer i
      do i=1,nrIBM_fib
       n_fib_nodes(i) = nIBM_r_fib(i)
      end do
 
      return
      end subroutine CTML_GET_FIB_NODE_COUNT

       ! dummy_module declared in: ../Vicar3D/UTIL_BOUNDARY_FORCE_FSI.F90
      subroutine CTML_GET_POS_VEL_WT(&
       fib_pst,&
       fib_vel,&
       fib_wt,&
       n_fib_bodies,&
       max_n_fib_nodes,&
       ifib)

      use dummy_module
      use probcommon_module

      IMPLICIT NONE
      INTEGER_T, INTENT(in) :: n_fib_bodies
      INTEGER_T, INTENT(in) :: max_n_fib_nodes
      INTEGER_T, INTENT(in) :: ifib
      REAL_T, INTENT(inout) :: fib_pst(n_fib_bodies,max_n_fib_nodes,SDIM)
      REAL_T, INTENT(inout) :: fib_vel(n_fib_bodies,max_n_fib_nodes,SDIM)
      REAL_T, INTENT(inout) :: fib_wt(n_fib_bodies,max_n_fib_nodes)
      INTEGER_T inode,idir,inode_cutoff

      if (1.eq.1) then
       print *,"calling CTML_GET_POS_VEL_WT"
       print *,"ifib=",ifib
      endif

      if (n_fib_bodies.ne.nrIBM_fib) then
       print *,"n_fib_bodies.ne.nrIBM_fib"
       stop
      endif
      if (max_n_fib_nodes.ne.nIBM_fib) then
       print *,"max_n_fib_nodes.ne.nIBM_fib"
       stop
      endif
      if ((ifib.ge.1).and.(ifib.le.nrIBM_fib)) then
       inode_cutoff=0
       do inode=1,nIBM_fib
        do idir=1,SDIM
         fib_pst(ifib,inode,idir)=coord_fib(ifib,inode,idir)
           ! section=1
         fib_vel(ifib,inode,idir)=vel_fib(ifib,1,inode,idir)
        end do
        fib_wt(ifib,inode)=ds_fib(ifib,inode)
        if (fib_wt(ifib,inode).gt.zero) then
         ! do nothing
        else if ((fib_wt(ifib,inode).eq.zero).and.(inode.gt.1)) then
         if (inode_cutoff.eq.0) then
          inode_cutoff=inode
         endif
        else 
         print *,"fib_wt(ifib,inode) invalid"
         print *,"ifib= ",ifib
         print *,"inode= ",inode
         print *,"nIBM_fib= ",nIBM_fib
         print *,"nrIBM_fib= ",nrIBM_fib
         print *,"ds_fib(ifib,inode)=",ds_fib(ifib,inode)
         stop
        endif
       end do
       if (inode_cutoff.ne.0) then
        print *,"WARNING, ds_fib==0 for some nodes"
        print *,"inode_cutoff=",inode_cutoff 
        print *,"ifib= ",ifib
        print *,"nIBM_fib= ",nIBM_fib
        print *,"nrIBM_fib= ",nrIBM_fib
       endif
      else
       print *,"ifib invalid"
       stop
      endif

      if (1.eq.1) then
       print *,"done with CTML_GET_POS_VEL_FORCE_WT"
       print *,"ifib=",ifib
      endif

      return
      end subroutine CTML_GET_POS_VEL_WT


       ! dummy_module declared in: ../Vicar3D/UTIL_BOUNDARY_FORCE_FSI.F90
      subroutine CTML_PUT_FORCE(&
        !force_fib(_prev)=fib_frc
       fib_frc,& !caller: ctml_fib_frc
       n_fib_bodies,&
       max_n_fib_nodes,&
       ifib)

      use dummy_module
      use probcommon_module

      IMPLICIT NONE
      INTEGER_T, INTENT(in) :: n_fib_bodies
      INTEGER_T, INTENT(in) :: max_n_fib_nodes
      INTEGER_T, INTENT(in) :: ifib
      REAL_T, INTENT(inout) :: fib_frc(n_fib_bodies,max_n_fib_nodes,SDIM)
      INTEGER_T inode,idir,inode_cutoff

      if (1.eq.1) then
       print *,"calling CTML_PUT_FORCE"
       print *,"ifib=",ifib
      endif

      if (n_fib_bodies.ne.nrIBM_fib) then
       print *,"n_fib_bodies.ne.nrIBM_fib"
       stop
      endif
      if (max_n_fib_nodes.ne.nIBM_fib) then
       print *,"max_n_fib_nodes.ne.nIBM_fib"
       stop
      endif
      if ((ifib.ge.1).and.(ifib.le.nrIBM_fib)) then
       inode_cutoff=0
       do inode=1,nIBM_fib
        do idir=1,SDIM
          ! section=1
         force_fib(ifib,1,inode,idir)= &
           fib_frc(ifib,inode,idir)
        enddo
       enddo
      else
       print *,"ifib invalid"
       stop
      endif

      if (1.eq.1) then
       print *,"done with CTML_PUT_FORCE"
       print *,"ifib=",ifib
      endif

      return
      end subroutine CTML_PUT_FORCE



      subroutine CTML_RESET_ARRAYS()

      use probcommon_module
      use dummy_module
      use probf90_module
      use flow_parameters 
      use grid_arrays

      IMPLICIT NONE

      if (1.eq.1) then
       print *,"calling CTML_RESET_ARRAYS"
      endif

      vel_fib = zero
      force_fib = zero

      if (1.eq.1) then
       print *,"done with CTML_RESET_ARRAYS"
      endif

      return
      end subroutine CTML_RESET_ARRAYS

       ! CTML_SOLVE_SOLID is called from CLSVOF_ReadNodes
      subroutine CTML_SOLVE_SOLID(&
       cur_time,& ! t^{n+1}
       dt,&
       step, &
       verbose, &
       plot_int,&
       io_proc)
      use dummy_module

      IMPLICIT NONE

      REAL_T, INTENT(in) :: cur_time ! t^{n+1}
      REAL_T, INTENT(in) :: dt
      INTEGER_T, INTENT(in) :: step
      INTEGER_T, INTENT(in) :: verbose
      INTEGER_T, INTENT(in) :: plot_int
      INTEGER_T, INTENT(in) :: io_proc
      INTEGER_T debug_tick
      INTEGER_T ifib,isec,inode,idir

      logical monitorON,theboss

      debug_tick=2

      if (debug_tick.ge.1) then
       print *,"calling CTML_SOLVE_SOLID"
       print *,"cur_time,dt,step,io_proc ",cur_time,dt,step,io_proc
      endif

      if(verbose.gt.0) then
       monitorON=.true.
      else
       monitorON=.false.
      end if

      if(io_proc.eq.1) then
       theboss=.true.
      else
       theboss=.false.
      end if

      if (debug_tick.ge.2) then
       print *,"tick cur_time,dt,step,io_proc = ",cur_time,dt,step,io_proc
       print *,"nrIBM_fib,nsecIBMmax,nIBM_fib ", &
        nrIBM_fib,nsecIBMmax,nIBM_fib
       do ifib=1,nrIBM_fib
        do isec=1,nsecIBMmax
         do inode=1,nIBM_fib
          print *,"ifib,inode,ds_fib ",ifib,inode,ds_fib(ifib,inode)
          if (ds_fib(ifib,inode).gt.zero) then
           do idir=1,3
            print *,"ifib,isec,inode,idir ", &
             ifib,isec,inode,idir
            print *,"vel_fib ",vel_fib(ifib,isec,inode,idir)
            print *,"coord_fib ",coord_fib(ifib,inode,idir)
           enddo  ! idir
          else if (ds_fib(ifib,inode).eq.zero) then
           ! do nothing
          else
           print *,"ds_fib invalid"
           stop
          endif
         enddo
        enddo
       enddo
      endif

       ! case insensitive search in vi:  ":set ignorecase" or ":set ic"
       ! case sensitive search in vi:  ":set noignorecase" or ":set noic"
       ! dummy_module declared in: ../Vicar3D/UTIL_BOUNDARY_FORCE_FSI.F90
       ! tick declared in: ../Vicar3D/distFSI/tick.F
       ! keyword: "INCLUDE_FIB"
      if (1.eq.1) then
#ifdef INCLUDE_FIB
       call tick_fib(cur_time, &  ! t^{n+1}
        dt,step,monitorON,plot_int, &
        !ctml_fib_vel_halftime_prev
        vel_fib_halftime_prev(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,1), & 
        vel_fib_halftime_prev(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,2), &
        vel_fib_halftime_prev(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,3), & 
        !ctml_fib_vel_prev
        vel_fib_prev(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,1), & 
        vel_fib_prev(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,2), &
        vel_fib_prev(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,3), & 
        !ctml_fib_vel_prev
        vel_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,1), & 
        vel_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,2), &
        vel_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,3), &
        force_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,1), & 
        force_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,2), &
        force_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,3), & 
        !ctml_fib_pst_prev 
        coord_fib_prev(1:nrIBM_fib,1:nIBM_fib,1), & 
        coord_fib_prev(1:nrIBM_fib,1:nIBM_fib,2), & 
        coord_fib_prev(1:nrIBM_fib,1:nIBM_fib,3), &
        coord_fib(1:nrIBM_fib,1:nIBM_fib,1), & 
        coord_fib(1:nrIBM_fib,1:nIBM_fib,2), & 
        coord_fib(1:nrIBM_fib,1:nIBM_fib,3), &
        theboss)
#else 
       print *,"include_fib not defined"
       stop
#endif
      else if (1.eq.0) then
       call tick(cur_time,  & ! t^{n+1}
        dt,step,monitorON,plot_int, &
        vel_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,1), & !ctml_fib_vel 
        vel_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,2), & !ctml_fib_vel
        vel_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,3), & !ctml_fib_vel
        vel_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,1), &
        vel_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,2), &
        vel_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,3),  &
        vel_esh(1:nrIBM_esh,1:nIBM_esh,1), &
        vel_esh(1:nrIBM_esh,1:nIBM_esh,2), &
        vel_esh(1:nrIBM_esh,1:nIBM_esh,3),  &
        vel_fbc(1:nrIBM_fbc,1:nIBM_fbc,1), &
        vel_fbc(1:nrIBM_fbc,1:nIBM_fbc,2), &
        vel_fbc(1:nrIBM_fbc,1:nIBM_fbc,3),  &
        force_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,1), & !ctml_fib_frc
        force_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,2), & !ctml_fib_frc
        force_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,3), & !ctml_fib_frc
        force_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,1), &
        force_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,2), &
        force_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,3), &
        force_esh(1:nrIBM_esh,1:nIBM_esh,1), &
        force_esh(1:nrIBM_esh,1:nIBM_esh,2), &
        force_esh(1:nrIBM_esh,1:nIBM_esh,3), &
        force_fbc(1:nrIBM_fbc,1:nIBM_fbc,1), &
        force_fbc(1:nrIBM_fbc,1:nIBM_fbc,2), &
        force_fbc(1:nrIBM_fbc,1:nIBM_fbc,3), &
        coord_fib(1:nrIBM_fib,1:nIBM_fib,1), & !ctml_fib_pst  
        coord_fib(1:nrIBM_fib,1:nIBM_fib,2), & !ctml_fib_pst
        coord_fib(1:nrIBM_fib,1:nIBM_fib,3), & !ctml_fib_pst
        coord_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,1),   &
        coord_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,2),   &
        coord_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,3),   &
        coord_esh(1:nrIBM_esh,1:nIBM_esh,1),   &
        coord_esh(1:nrIBM_esh,1:nIBM_esh,2),   &
        coord_esh(1:nrIBM_esh,1:nIBM_esh,3),   &
        coord_fbc(1:nrIBM_fbc,1:nIBM_fbc,1),   &
        coord_fbc(1:nrIBM_fbc,1:nIBM_fbc,2),   &
        coord_fbc(1:nrIBM_fbc,1:nIBM_fbc,3),   &
        theboss)
      else
       print *,"tick not called; if statement corrupt"
       stop
      endif

      if (debug_tick.ge.2) then
       print *,"tick cur_time,dt,step = ",cur_time,dt,step
       print *,"nrIBM_fib,nsecIBMmax,nIBM_fib ", &
        nrIBM_fib,nsecIBMmax,nIBM_fib
       do ifib=1,nrIBM_fib
        do isec=1,nsecIBMmax
         do inode=1,nIBM_fib
          if (ds_fib(ifib,inode).gt.zero) then
           do idir=1,3
            print *,"ifib,isec,inode,idir ", &
             ifib,isec,inode,idir
            print *,"force_fib ",force_fib(ifib,isec,inode,idir)
            print *,"coord_fib (after) ",coord_fib(ifib,inode,idir)
           enddo ! idir
          else if (ds_fib(ifib,inode).eq.zero) then
           ! do nothing
          else
           print *,"ds_fib invalid"
           stop
          endif
         enddo ! inode
        enddo ! isec
       enddo ! ifib
      endif

      if (debug_tick.ge.1) then
       print *,"done with CTML_SOLVE_SOLID"
       print *,"cur_time,dt,step ",cur_time,dt,step
      endif

      return
      end subroutine CTML_SOLVE_SOLID


      end module CTML_module


subroutine fort_ctml_max_nodes(&
 nmat_in,&
 FSI_flag_in,&
 FSI_CTML_max_num_nodes_list,&
 FSI_CTML_max_num_elements_list) &
bind(c,name='fort_ctml_max_nodes')

 use probcommon_module
 use global_utility_module
 use dummy_module

 IMPLICIT NONE
 INTEGER_T, INTENT(in) :: nmat_in
 INTEGER_T, INTENT(in) :: FSI_flag_in(nmat_in)
 INTEGER_T, INTENT(inout) :: FSI_CTML_max_num_nodes_list(nmat_in)
 INTEGER_T, INTENT(inout) :: FSI_CTML_max_num_elements_list(nmat_in)
 INTEGER_T :: im

 do im=1,nmat_in
  if (FSI_flag_in(im).eq.FSI_SHOELE_CTML) then
#ifdef MVAHABFSI
   if (nIBM_fib.ge.2) then
    FSI_CTML_max_num_nodes_list(im)=nIBM_fib
    FSI_CTML_max_num_elements_list(im)=nIBM_fib-1
   else
    print *,"nIBM_fib invalid"
    stop
   endif
#else
   print *,"CTML(F): define MEHDI_VAHAB_FSI in GNUmakefile"
   stop
#endif
  else
   FSI_CTML_max_num_nodes_list(im)=0 
   FSI_CTML_max_num_elements_list(im)=0 
  endif
 enddo !im=1,nmat_in

 return
end subroutine fort_ctml_max_nodes


