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

      subroutine CTML_INTERNAL_MAX_NODES(&
       nmat_in,&
       FSI_flag_in,&
       CTML_num_solids_out, &
       CTML_max_num_nodes_list,&
       CTML_max_num_elements_list)

      IMPLICIT NONE
      INTEGER_T, INTENT(in) :: nmat_in
      INTEGER_T, INTENT(in) :: FSI_flag_in(nmat_in)
      INTEGER_T, INTENT(inout) :: CTML_num_solids_out
      INTEGER_T, INTENT(inout) :: CTML_max_num_nodes_list(3)
      INTEGER_T, INTENT(inout) :: CTML_max_num_elements_list
      INTEGER_T :: im
      INTEGER_T :: Ns_IBM_fib_out
      INTEGER_T :: Ns_IBM_fsh_out
      INTEGER_T :: Nq_IBM_fsh_out
      INTEGER_T :: Ns_IBM_esh_out
      INTEGER_T :: Ns_IBM_fbc_out

      INTEGER_T :: Nr_IBM_fib_out
      INTEGER_T :: Nr_IBM_fsh_out
      INTEGER_T :: Nr_IBM_esh_out
      INTEGER_T :: Nr_IBM_fbc_out

      INTEGER_T :: CTML_num_solids_local
      INTEGER_T :: dir

      CTML_num_solids_out=0

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
         CTML_max_num_nodes_list(1)=Ns_IBM_fib_out
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
 CTML_num_solids_out, &
 CTML_max_num_nodes_list,&
 CTML_max_num_elements_list) &
bind(c,name='fort_ctml_max_nodes')

use CTML_module

IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat_in
INTEGER_T, INTENT(in) :: FSI_flag_in(nmat_in)
INTEGER_T, INTENT(inout) :: CTML_num_solids_out
INTEGER_T, INTENT(inout) :: CTML_max_num_nodes_list(3)
INTEGER_T, INTENT(inout) :: CTML_max_num_elements_list
INTEGER_T :: dir

CTML_num_solids_out=0

do dir=1,3
 CTML_max_num_nodes_list(dir)=0
enddo

CTML_max_num_elements_list=0

call CTML_INTERNAL_MAX_NODES(&
 nmat_in,&
 FSI_flag_in,&
 CTML_num_solids_out,&
 CTML_max_num_nodes_list,&
 CTML_max_num_elements_list) 

return
end subroutine fort_ctml_max_nodes

