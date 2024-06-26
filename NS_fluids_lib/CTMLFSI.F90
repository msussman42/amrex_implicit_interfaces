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

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"
#include "EXTRAP_COMP.H"

#include "CTMLFSI_F.H"

      module CTML_module
      use amrex_fort_module, only : amrex_real
      use probf90_module

      contains

      subroutine CTML_INTERNAL_MAX_NODES(&
       nmat_in,&
       FSI_flag_in,&
       CTML_num_solids_out, &
       CTML_num_scalars_out, &
       CTML_max_num_nodes_list,&
       CTML_max_num_elements_list)

      IMPLICIT NONE
      integer, INTENT(in) :: nmat_in
      integer, INTENT(in) :: FSI_flag_in(nmat_in)
      integer, INTENT(inout) :: CTML_num_solids_out
      integer, INTENT(inout) :: CTML_num_scalars_out
      integer, INTENT(inout) :: CTML_max_num_nodes_list(3)
      integer, INTENT(inout) :: CTML_max_num_elements_list
      integer :: im
      integer :: Ns_IBM_fib_out
      integer :: Ns_IBM_fsh_out
      integer :: Nq_IBM_fsh_out
      integer :: Ns_IBM_esh_out
      integer :: Ns_IBM_fbc_out

      integer :: Nr_IBM_fib_out
      integer :: Nr_IBM_fsh_out
      integer :: Nr_IBM_esh_out
      integer :: Nr_IBM_fbc_out

      integer :: CTML_num_solids_local
      integer :: dir

      CTML_num_solids_out=0
      CTML_num_scalars_out=0

      do dir=1,3
       CTML_max_num_nodes_list(dir)=0
      enddo

      CTML_max_num_elements_list=0

      CTML_num_solids_local=0
      do im=1,nmat_in
       if (FSI_flag_in(im).eq.FSI_SHOELE_CTML) then
        CTML_num_solids_local=CTML_num_solids_local+1
       endif
      enddo

      if (CTML_num_solids_local.eq.0) then
       !do nothing
      else if ((CTML_num_solids_local.ge.1).and. &
               (CTML_num_solids_local.le.nmat_in-1)) then

#ifdef MVAHABFSI
        call copy_nmaxIBM( &
          CTML_num_solids_out, &
          CTML_num_scalars_out, &
          nr_IBM_fib_out, &
          nr_IBM_fsh_out, &
          nr_IBM_esh_out, &
          nr_IBM_fbc_out, &
          ns_IBM_fib_out, &
          ns_IBM_fsh_out, &
          nq_IBM_fsh_out, &
          ns_IBM_esh_out, &
          ns_IBM_fbc_out)

        if (CTML_num_solids_out.eq.CTML_num_solids_local) then
         ! do nothing
        else
         print *,"CTML_num_solids_out invalid"
         print *,"CTML_num_solids_out: ",CTML_num_solids_out
         print *,"CTML_num_solids_local: ",CTML_num_solids_local
         stop
        endif

        if (AMREX_SPACEDIM.eq.2) then
         if (Nr_IBM_fib_out.eq.CTML_num_solids_out) then
          ! do nothing
         else
          print *,"expecting Nr_IBM_fib_out.eq.CTML_num_solids_out"
          stop
         endif
         if (Nr_IBM_fsh_out.eq.0) then
          ! do nothing
         else
          print *,"Nr_IBM_fsh_out invalid"
          stop
         endif
        else if (AMREX_SPACEDIM.eq.3) then
         if (Nr_IBM_fib_out.eq.0) then
          ! do nothing
         else
          print *,"expecting Nr_IBM_fib_out.eq.0"
          stop
         endif
         if (Nr_IBM_fsh_out.eq.CTML_num_solids_out) then
          ! do nothing
         else
          print *,"Nr_IBM_fsh_out invalid"
          stop
         endif
        else
         print *,"AMREX_SPACEDIM invalid"
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

        if (AMREX_SPACEDIM.eq.2) then
         if (Ns_IBM_fib_out.ge.2) then
          CTML_max_num_nodes_list(1)=Ns_IBM_fib_out
          CTML_max_num_elements_list=Ns_IBM_fib_out-1
         else
          print *,"Ns_IBM_fib_out invalid"
          stop
         endif
        else if (AMREX_SPACEDIM.eq.3) then
         if ((Nq_IBM_fsh_out.ge.2).and. &
             (Ns_IBM_fsh_out.ge.2)) then
          CTML_max_num_nodes_list(1)=Nq_IBM_fsh_out
          CTML_max_num_nodes_list(2)=Ns_IBM_fsh_out
          CTML_max_num_elements_list= &
              (Ns_IBM_fsh_out-1)*(Nq_IBM_fsh_out-1)
         else
          print *,"N[sq]_IBM_fsh_out invalid"
          stop
         endif
        else
         print *,"AMREX_SPACEDIM invalid"
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
 CTML_num_solids_out, &
 CTML_num_scalars_out, &
 CTML_max_num_nodes_list,&
 CTML_max_num_elements_list) &
bind(c,name='fort_ctml_max_nodes')

use CTML_module

IMPLICIT NONE

integer, INTENT(in) :: nmat_in
integer, INTENT(in) :: FSI_flag_in(nmat_in)
integer, INTENT(inout) :: CTML_num_solids_out
integer, INTENT(inout) :: CTML_num_scalars_out
integer, INTENT(inout) :: CTML_max_num_nodes_list(3)
integer, INTENT(inout) :: CTML_max_num_elements_list
integer :: dir

CTML_num_solids_out=0
CTML_num_scalars_out=0

do dir=1,3
 CTML_max_num_nodes_list(dir)=0
enddo

CTML_max_num_elements_list=0

call CTML_INTERNAL_MAX_NODES(&
 nmat_in,&
 FSI_flag_in,&
 CTML_num_solids_out,&
 CTML_num_scalars_out,&
 CTML_max_num_nodes_list,&
 CTML_max_num_elements_list) 

return
end subroutine fort_ctml_max_nodes

