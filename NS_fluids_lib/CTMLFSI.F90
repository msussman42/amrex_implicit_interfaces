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

      subroutine CTML_INTERNAL_MAX_NODES(&
       nmat_in,&
       FSI_flag_in,&
       CTML_max_num_nodes_list,&
       CTML_max_num_elements_list)

      IMPLICIT NONE
      INTEGER_T, INTENT(in) :: nmat_in
      INTEGER_T, INTENT(in) :: FSI_flag_in(nmat_in)
      INTEGER_T, INTENT(inout) :: CTML_max_num_nodes_list
      INTEGER_T, INTENT(inout) :: CTML_max_num_elements_list
      INTEGER_T :: im
      INTEGER_T :: Ns_IBM_fib_out
      INTEGER_T :: Ns_IBM_fsh_out
      INTEGER_T :: Nq_IBM_fsh_out
      INTEGER_T :: Ns_IBM_esh_out
      INTEGER_T :: Ns_IBM_fbc_out

      INTEGER_T :: Nr_IBM_out
      INTEGER_T :: Nr_IBM_fib_out
      INTEGER_T :: Nr_IBM_fsh_out
      INTEGER_T :: Nr_IBM_esh_out
      INTEGER_T :: Nr_IBM_fbc_out

      INTEGER_T :: CTML_num_solids_local

      CTML_max_num_nodes_list=0
      CTML_max_num_elements_list=0

      CTML_num_solids_local=0
      do im=1,nmat_in
       if (FSI_flag_in(im).eq.FSI_SHOELE_CTML) then
        CTML_num_solids_local=CTML_num_solids_local+1
       endif
      endif

      if (CTML_num_solids_local.eq.0) then
       !do nothing
      else if ((CTML_num_solids_local.ge.1).and. &
               (CTML_num_solids_local.le.nmat_in-1)) then

#ifdef MVAHABFSI
        call copy_nmaxIBM( &
          nr_IBM_out, &
          nr_IBM_fib_out, &
          nr_IBM_fsh_out, &
          nr_IBM_esh_out, &
          nr_IBM_fbc_out, &
          ns_IBM_fib_out, &
          ns_IBM_fsh_out, &
          nq_IBM_fsh_out, &
          ns_IBM_esh_out, &
          ns_IBM_fbc_out)

        if (Nr_IBM_out.eq.CTML_num_solids_local) then
         ! do nothing
        else
         print *,"Nr_IBM_out invalid"
         stop
        endif
        if (Nr_IBM_fib_out.eq.Nr_IBM_out) then
         ! do nothing
        else
         print *,"expecting Nr_IBM_fib_out.eq.Nr_IBM_out"
         stop
        endif
        if (Nr_IBM_fsh_out.eq.0) then
         ! do nothing
        else
         print *,"Nr_IBM_fsh_out invalid"
         stop
        endif
        if (Nr_IBM_esh_out.eq.0) then
         ! do nothing
        else
         print *,"Nr_IBM_esh_out invalid"
         stop
        endif
        if (Nr_IBM_fbc_out.eq.0) then
         ! do nothing
        else
         print *,"Nr_IBM_fbc_out invalid"
         stop
        endif

        if (Ns_IBM_fib_out.ge.2) then
         CTML_max_num_nodes_list=Ns_IBM_fib_out
         CTML_max_num_elements_list=Ns_IBM_fib_out-1
        else
         print *,"Ns_IBM_fib_out invalid"
         stop
        endif
#else
        print *,"CTML(F): define MEHDI_VAHAB_FSI in GNUmakefile"
        stop
#endif
      else
       print *,"CTML_num_solids_local invalid"
       stop
      endif

      return
      end subroutine CTML_INTERNAL_MAX_NODES

      end module CTML_module

subroutine fort_ctml_max_nodes(&
 nmat_in,&
 FSI_flag_in,&
 CTML_max_num_nodes_list,&
 CTML_max_num_elements_list) &
bind(c,name='fort_ctml_max_nodes')

use CTML_module

IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat_in
INTEGER_T, INTENT(in) :: FSI_flag_in(nmat_in)
INTEGER_T, INTENT(inout) :: CTML_max_num_nodes_list
INTEGER_T, INTENT(inout) :: CTML_max_num_elements_list

CTML_max_num_nodes_list=0
CTML_max_num_elements_list=0

call CTML_INTERNAL_MAX_NODES(&
 nmat_in,&
 FSI_flag_in,&
 CTML_max_num_nodes_list,&
 CTML_max_num_elements_list) 

return
end subroutine fort_ctml_max_nodes


