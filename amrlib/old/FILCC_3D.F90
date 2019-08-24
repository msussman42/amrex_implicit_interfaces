#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "ArrayLim.H"

#define SDIM 3

! ::: -----------------------------------------------------------
! ::: This routine is intended to be a generic fill function
! ::: for cell centered data.  It knows how to exrapolate,
! ::: and reflect data and can be used to suppliment problem
! ::: specific fill functions (ie. EXT_DIR).
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: q        <=  array to fill
! ::: DIMS(q)   => index extent of q array
! ::: domlo,hi  => index extent of problem domain
! ::: bc	=> array of boundary flags bc(SPACEDIM,lo:hi)
! ::: 
! ::: NOTE: corner data not used in computing soln but must have
! :::       reasonable values for arithmetic to live
! ::: -----------------------------------------------------------

      subroutine filcc( &
       q,DIMS(q), &
       domlo,domhi,bc)

      INTEGER_T    DIMDEC(q)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      REAL_T     q(DIMV(q))
      INTEGER_T    bc(SDIM,2)

      INTEGER_T    nlft, nrgt, nbot, ntop, nup, ndwn
      INTEGER_T    ilo, ihi, jlo, jhi, klo, khi
      INTEGER_T    is,  ie,  js,  je,  ks,  ke
      INTEGER_T    i, j, k

      is = max(ARG_L1(q),domlo(1))
      ie = min(ARG_H1(q),domhi(1))
      js = max(ARG_L2(q),domlo(2))
      je = min(ARG_H2(q),domhi(2))
      ks = max(ARG_L3(q),domlo(3))
      ke = min(ARG_H3(q),domhi(3))

      nlft = max(0,domlo(1)-ARG_L1(q))
      nrgt = max(0,ARG_H1(q)-domhi(1))
      nbot = max(0,domlo(2)-ARG_L2(q))
      ntop = max(0,ARG_H2(q)-domhi(2))
      ndwn = max(0,domlo(3)-ARG_L3(q))
      nup  = max(0,ARG_H3(q)-domhi(3))
!
!     ::::: first fill sides
!
      if (nlft .gt. 0) then
         ilo = domlo(1)

	 if (bc(1,1) .eq. FOEXTRAP) then
	    do i = 1, nlft
               do k = ARG_L3(q),ARG_H3(q)
                  do j = ARG_L2(q),ARG_H2(q)
                     q(ilo-i,j,k) = q(ilo,j,k)
                  end do
               end do
	    end do
	 else if (bc(1,1) .eq. HOEXTRAP) then
	    do i = 1, nlft
               do k = ARG_L3(q),ARG_H3(q)
                  do j = ARG_L2(q),ARG_H2(q)
                     q(ilo-i,j,k) = two*q(ilo-i+1,j,k)-q(ilo-i+2,j,k)
                  end do
               end do
	    end do
	 else if (bc(1,1) .eq. REFLECT_EVEN) then
	    do i = 1, nlft
               do k = ARG_L3(q),ARG_H3(q)
                  do j = ARG_L2(q),ARG_H2(q)
                     q(ilo-i,j,k) = q(ilo+i-1,j,k)
                  end do
               end do
	    end do
	 else if (bc(1,1) .eq. REFLECT_ODD) then
	    do i = 1, nlft
               do k = ARG_L3(q),ARG_H3(q)
                  do j = ARG_L2(q),ARG_H2(q)
                     q(ilo-i,j,k) = -q(ilo+i-1,j,k)
                  end do
               end do
	    end do
	 end if
      end if

      if (nrgt .gt. 0) then
         ihi = domhi(1)

	 if (bc(1,2) .eq. FOEXTRAP) then
	    do i = 1, nrgt
               do k = ARG_L3(q),ARG_H3(q)
                  do j = ARG_L2(q),ARG_H2(q)
                     q(ihi+i,j,k) = q(ihi,j,k)
                  end do
               end do
	    end do
         else if (bc(1,2) .eq. HOEXTRAP) then
            do i = 1, nrgt
               do k = ARG_L3(q),ARG_H3(q)
                  do j = ARG_L2(q),ARG_H2(q)
                     q(ihi+i,j,k) = two*q(ihi+i-1,j,k)-q(ihi+i-2,j,k)
                  end do
               end do
            end do
	 else if (bc(1,2) .eq. REFLECT_EVEN) then
	    do i = 1, nrgt
               do k = ARG_L3(q),ARG_H3(q)
                  do j = ARG_L2(q),ARG_H2(q)
                     q(ihi+i,j,k) = q(ihi-i+1,j,k)
                  end do
               end do
	    end do
	 else if (bc(1,2) .eq. REFLECT_ODD) then
	    do i = 1, nrgt
               do k = ARG_L3(q),ARG_H3(q)
                  do j = ARG_L2(q),ARG_H2(q)
                     q(ihi+i,j,k) = -q(ihi-i+1,j,k)
                  end do
               end do
	    end do
	 end if
      end if

      if (nbot .gt. 0) then
         jlo = domlo(2)
         
	 if (bc(2,1) .eq. FOEXTRAP) then
	    do j = 1, nbot
               do k = ARG_L3(q),ARG_H3(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,jlo-j,k) = q(i,jlo,k)
                  end do
               end do
	    end do
         else if (bc(2,1) .eq. HOEXTRAP) then
            do j = 1, nbot
               do k = ARG_L3(q),ARG_H3(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,jlo-j,k) = two*q(i,jlo-j+1,k)-q(i,jlo-j+2,k)
                  end do
               end do
            end do
	 else if (bc(2,1) .eq. REFLECT_EVEN) then
	    do j = 1, nbot 
               do k = ARG_L3(q),ARG_H3(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,jlo-j,k) = q(i,jlo+j-1,k)
                  end do
               end do
	    end do
	 else if (bc(2,1) .eq. REFLECT_ODD) then
	    do j = 1, nbot
               do k = ARG_L3(q),ARG_H3(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,jlo-j,k) = -q(i,jlo+j-1,k)
                  end do
               end do
	    end do
	 end if
      end if

      if (ntop .gt. 0) then
         jhi = domhi(2)

	 if (bc(2,2) .eq. FOEXTRAP) then
	    do j = 1, ntop
               do k = ARG_L3(q),ARG_H3(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,jhi+j,k) = q(i,jhi,k)
                  end do
               end do
	    end do
         else if (bc(2,2) .eq. HOEXTRAP) then
            do j = 1, ntop
               do k = ARG_L3(q),ARG_H3(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,jhi+j,k) = two*q(i,jhi+j-1,k)-q(i,jhi+j-2,k)
                  end do
               end do
            end do
	 else if (bc(2,2) .eq. REFLECT_EVEN) then
	    do j = 1, ntop
               do k = ARG_L3(q),ARG_H3(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,jhi+j,k) = q(i,jhi-j+1,k)
                  end do
               end do
	    end do
	 else if (bc(2,2) .eq. REFLECT_ODD) then
	    do j = 1, ntop
               do k = ARG_L3(q),ARG_H3(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,jhi+j,k) = -q(i,jhi-j+1,k)
                  end do
               end do
	    end do
	 end if
      end if

      if (ndwn .gt. 0) then
         klo = domlo(3)

	 if (bc(3,1) .eq. FOEXTRAP) then
	    do k = 1, ndwn
               do j = ARG_L2(q),ARG_H2(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,j,klo-k) = q(i,j,klo)
                  end do
               end do
	    end do
         else if (bc(3,1) .eq. HOEXTRAP) then
            do k = 1, ndwn
               do j = ARG_L2(q),ARG_H2(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,j,klo-k) = two*q(i,j,klo-k+1)-q(i,j,klo-k+2)
                  end do
               end do
            end do
	 else if (bc(3,1) .eq. REFLECT_EVEN) then
	    do k = 1, ndwn
               do j = ARG_L2(q),ARG_H2(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,j,klo-k) = q(i,j,klo+k-1)
                  end do
               end do
	    end do
	 else if (bc(3,1) .eq. REFLECT_ODD) then
	    do k = 1, ndwn
               do j = ARG_L2(q),ARG_H2(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,j,klo-k) = -q(i,j,klo+k-1)
                  end do
               end do
	    end do
	 end if
      end if

      if (nup .gt. 0) then
         khi = domhi(3)

	 if (bc(3,2) .eq. FOEXTRAP) then
	    do k = 1, nup
               do j = ARG_L2(q),ARG_H2(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,j,khi+k) = q(i,j,khi)
                  end do
               end do
	    end do
         else if (bc(3,2) .eq. HOEXTRAP) then
            do k = 1, nup
               do j = ARG_L2(q),ARG_H2(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,j,khi+k) = two*q(i,j,khi+k-1)-q(i,j,khi+k-2)
                  end do
               end do
            end do
	 else if (bc(3,2) .eq. REFLECT_EVEN) then
	    do k = 1, nup
               do j = ARG_L2(q),ARG_H2(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,j,khi+k) = q(i,j,khi-k+1)
                  end do
               end do
	    end do
	 else if (bc(3,2) .eq. REFLECT_ODD) then
	    do k = 1, nup
               do j = ARG_L2(q),ARG_H2(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,j,khi+k) = -q(i,j,khi-k+1)
                  end do
               end do
	    end do
	 end if
      end if

      end

