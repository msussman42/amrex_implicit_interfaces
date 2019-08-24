
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "ArrayLim.H"

#define SDIM 2

! ::: -----------------------------------------------------------
! ::: This routine is intended to be a generic fill function
! ::: for cell-centered data.  It knows how to extrapolate
! ::: and reflect data and is used to supplement the problem-specific
! ::: fill functions which call it.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: q           <=  array to fill
! ::: lo,hi        => index extent of q array
! ::: domlo,domhi  => index extent of problem domain
! ::: bc	   => array of boundary flags bc(SPACEDIM,lo:hi)
! ::: 
! ::: NOTE: all corner as well as edge data is filled if not EXT_DIR
! ::: -----------------------------------------------------------

      subroutine filcc( &
       q,DIMS(q), &
       domlo,domhi,bc)

      INTEGER_T    DIMDEC(q)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      INTEGER_T    bc(SDIM,2)
      REAL_T     q(DIMV(q))

      INTEGER_T    nlft, nrgt, nbot, ntop
      INTEGER_T    ilo, ihi, jlo, jhi
      INTEGER_T    i, j
      INTEGER_T    is, ie, js, je

      nlft = max(0,domlo(1)-ARG_L1(q))
      nrgt = max(0,ARG_H1(q)-domhi(1))
      nbot = max(0,domlo(2)-ARG_L2(q))
      ntop = max(0,ARG_H2(q)-domhi(2))

      is = max(ARG_L1(q),domlo(1))
      ie = min(ARG_H1(q),domhi(1))
      js = max(ARG_L2(q),domlo(2))
      je = min(ARG_H2(q),domhi(2))

!     ::::: first fill sides
      if (nlft .gt. 0) then
         ilo = domlo(1)

	 if (bc(1,1) .eq. FOEXTRAP) then
	    do i = 1, nlft
	    do j = ARG_L2(q), ARG_H2(q)
	       q(ilo-i,j) = q(ilo,j)
	    end do
	    end do
	 else if (bc(1,1) .eq. HOEXTRAP) then
	    do i = 1, nlft
	    do j = ARG_L2(q), ARG_H2(q)
	       q(ilo-i,j) = two*q(ilo-i+1,j)-q(ilo-i+2,j)
	    end do
	    end do
	 else if (bc(1,1) .eq. REFLECT_EVEN) then
	    do i = 1, nlft
	     do j = ARG_L2(q), ARG_H2(q)
	       q(ilo-i,j) = q(ilo+i-1,j)
	    end do
	    end do
	 else if (bc(1,1) .eq. REFLECT_ODD) then
	    do i = 1, nlft
	    do j = ARG_L2(q), ARG_H2(q)
	       q(ilo-i,j) = -q(ilo+i-1,j)
	    end do
	    end do
         else if ((bc(1,1).ne.INT_DIR).and.(bc(1,1).ne.EXT_DIR)) then
          print *,"bc invalid"
          stop
	 end if
      end if

      if (nrgt .gt. 0) then
         ihi = domhi(1)

	 if (bc(1,2) .eq. FOEXTRAP) then
	    do i = 1, nrgt
	    do j = ARG_L2(q), ARG_H2(q)
	       q(ihi+i,j) = q(ihi,j)
	    end do
	    end do
         else if (bc(1,2) .eq. HOEXTRAP) then
	    do i = 1, nrgt
	    do j = ARG_L2(q), ARG_H2(q)
	       q(ihi+i,j) = two*q(ihi+i-1,j)-q(ihi+i-2,j)
	    end do
	    end do
	 else if (bc(1,2) .eq. REFLECT_EVEN) then
	    do i = 1, nrgt
            do j = ARG_L2(q), ARG_H2(q)
	       q(ihi+i,j) = q(ihi-i+1,j)
	    end do
	    end do
	 else if (bc(1,2) .eq. REFLECT_ODD) then
	    do i = 1, nrgt
            do j = ARG_L2(q), ARG_H2(q)
	       q(ihi+i,j) = -q(ihi-i+1,j)
	    end do
	    end do
         else if ((bc(1,2).ne.INT_DIR).and.(bc(1,2).ne.EXT_DIR)) then
          print *,"bc invalid"
          stop
	 end if
      end if

      if (nbot .gt. 0) then
         jlo = domlo(2)

	 if (bc(2,1) .eq. FOEXTRAP) then
	    do j = 1, nbot
	    do i = ARG_L1(q), ARG_H1(q)
	       q(i,jlo-j) = q(i,jlo)
	    end do
	    end do
         else if (bc(2,1) .eq. HOEXTRAP) then
	    do j = 1, nbot
	    do i = ARG_L1(q), ARG_H1(q)
	       q(i,jlo-j) = two*q(i,jlo-j+1)-q(i,jlo-j+2)
	    end do
	    end do
	 else if (bc(2,1) .eq. REFLECT_EVEN) then
	    do j = 1, nbot
 	    do i = ARG_L1(q), ARG_H1(q)
	       q(i,jlo-j) = q(i,jlo+j-1)
	    end do
	    end do
	 else if (bc(2,1) .eq. REFLECT_ODD) then
	    do j = 1, nbot
 	    do i = ARG_L1(q), ARG_H1(q)
	       q(i,jlo-j) = -q(i,jlo+j-1)
	    end do
	    end do
         else if ((bc(2,1).ne.INT_DIR).and.(bc(2,1).ne.EXT_DIR)) then
          print *,"bc invalid"
          stop
	 end if
      end if

      if (ntop .gt. 0) then
         jhi = domhi(2)

	 if (bc(2,2) .eq. FOEXTRAP) then
	    do j = 1, ntop
 	    do i = ARG_L1(q), ARG_H1(q)
	       q(i,jhi+j) = q(i,jhi)
	    end do
	    end do
         else if (bc(2,2) .eq. HOEXTRAP) then
	    do j = 1, ntop
 	    do i = ARG_L1(q), ARG_H1(q)
	       q(i,jhi+j) = two*q(i,jhi+j-1)-q(i,jhi+j-2)
	    end do
	    end do
	 else if (bc(2,2) .eq. REFLECT_EVEN) then
	    do j = 1, ntop
 	    do i = ARG_L1(q), ARG_H1(q)
	       q(i,jhi+j) = q(i,jhi-j+1)
	    end do
	    end do
	 else if (bc(2,2) .eq. REFLECT_ODD) then
	    do j = 1, ntop
 	    do i = ARG_L1(q), ARG_H1(q)
	       q(i,jhi+j) = -q(i,jhi-j+1)
	    end do
	    end do
         else if ((bc(2,2).ne.INT_DIR).and.(bc(2,2).ne.EXT_DIR)) then
          print *,"bc invalid"
          stop
	 end if
      end if

      return
      end

