
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "SPACE.H"
#include "BC_TYPES.H"
#include "ArrayLim.H"

#define SDIM 2


      subroutine filcc(bfact, &
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
      INTEGER_T    bfact

      nlft = max(0,domlo(1)-ARG_L1(q))
      nrgt = max(0,ARG_H1(q)-domhi(1))
      nbot = max(0,domlo(2)-ARG_L2(q))
      ntop = max(0,ARG_H2(q)-domhi(2))

      is = max(ARG_L1(q),domlo(1))
      ie = min(ARG_H1(q),domhi(1))
      js = max(ARG_L2(q),domlo(2))
      je = min(ARG_H2(q),domhi(2))

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif

!     ::::: first fill sides
      if (nlft .gt. 0) then
         ilo = domlo(1)

         if ((bc(1,1).eq.FOEXTRAP).or. &
             (bc(1,1).eq.EXT_DIR)) then
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
         else if (bc(1,1).eq.INT_DIR) then
          ! do nothing
         else
          print *,"bc invalid"
          stop
         end if
      end if

      if (nrgt .gt. 0) then
         ihi = domhi(1)

         if ((bc(1,2).eq.FOEXTRAP).or. &
             (bc(1,2).eq.EXT_DIR)) then
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
         else if (bc(1,2).eq.INT_DIR) then
          ! do nothing
         else
          print *,"bc invalid"
          stop
	 end if
      end if

      if (nbot .gt. 0) then
         jlo = domlo(2)

         if ((bc(2,1).eq.FOEXTRAP).or. &
             (bc(2,1).eq.EXT_DIR)) then
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
         else if (bc(2,1).eq.INT_DIR) then
          ! do nothing
         else
          print *,"bc invalid"
          stop
	 end if
      end if

      if (ntop .gt. 0) then
         jhi = domhi(2)

         if ((bc(2,2).eq.FOEXTRAP).or. &
             (bc(2,2).eq.EXT_DIR)) then
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
         else if (bc(2,2).eq.INT_DIR) then
          ! do nothing
         else
          print *,"bc invalid"
          stop
	 end if
      end if

      return
      end


! domlo,domhi are dimensions for face quantity (not cell) 
      subroutine efilcc(bfact, &
       q,DIMS(q), &
       domlo,domhi,bc,dir)
      IMPLICIT NONE

      INTEGER_T    DIMDEC(q)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      INTEGER_T    bc(SDIM,2),dir,bfact
      REAL_T     q(DIMV(q))

      INTEGER_T    nlft, nrgt, nbot, ntop
      INTEGER_T    ilo, ihi, jlo, jhi
      INTEGER_T    i, j
      INTEGER_T    is, ie, js, je

      if ((dir.lt.0).or.(dir.gt.1)) then
       print *,"dir out of range in EFILCC"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      nlft = max(0,domlo(1)-ARG_L1(q))
      nrgt = max(0,ARG_H1(q)-domhi(1))
      nbot = max(0,domlo(2)-ARG_L2(q))
      ntop = max(0,ARG_H2(q)-domhi(2))

      is = max(ARG_L1(q),domlo(1))
      ie = min(ARG_H1(q),domhi(1))
      js = max(ARG_L2(q),domlo(2))
      je = min(ARG_H2(q),domhi(2))

!     ::::: first fill sides
      if (nlft .ge. 0) then
         ilo = domlo(1)

         if ((bc(1,1).eq.FOEXTRAP).or. &
             (bc(1,1).eq.EXT_DIR).or. &
             (bc(1,1).eq.HOEXTRAP)) then
	    do i = 1, nlft
	    do j = ARG_L2(q), ARG_H2(q)
	       q(ilo-i,j) = q(ilo,j)
	    enddo
	    enddo
	 elseif (bc(1,1) .eq. REFLECT_EVEN) then

	    do j = ARG_L2(q), ARG_H2(q)
	    do i = 1, nlft
              if (dir.eq.0) then
	       q(ilo-i,j) = q(ilo+i,j)
              else
	       q(ilo-i,j) = q(ilo+i-1,j)
              endif
	    enddo
	    enddo
	 elseif (bc(1,1) .eq. REFLECT_ODD) then
	    do j = ARG_L2(q), ARG_H2(q)
              if (dir.eq.0) then
               q(ilo,j)=zero
              endif
	    do i = 1, nlft
              if (dir.eq.0) then
	       q(ilo-i,j) = -q(ilo+i,j)
              else
	       q(ilo-i,j) = -q(ilo+i-1,j)
              endif
	    enddo
	    enddo
         else if (bc(1,1).ne.INT_DIR) then
          print *,"bc invalid"
          stop
	 endif
      endif

      if (nrgt .ge. 0) then
         ihi = domhi(1)

         if ((bc(1,2).eq.FOEXTRAP).or. &
             (bc(1,2).eq.EXT_DIR).or. &
             (bc(1,2).eq.HOEXTRAP)) then
	    do i = 1, nrgt
	    do j = ARG_L2(q), ARG_H2(q)
	       q(ihi+i,j) = q(ihi,j)
	    enddo
	    enddo
	 elseif (bc(1,2) .eq. REFLECT_EVEN) then

	    do i = 1, nrgt
            do j = ARG_L2(q), ARG_H2(q)
             if (dir.eq.0) then
	       q(ihi+i,j) = q(ihi-i,j)
             else
	       q(ihi+i,j) = q(ihi-i+1,j)
             endif
	    enddo
	    enddo
	 elseif (bc(1,2) .eq. REFLECT_ODD) then
            do j = ARG_L2(q), ARG_H2(q)
              if (dir.eq.0) then
               q(ihi,j)=zero
              endif
	    do i = 1, nrgt
             if (dir.eq.0) then
	       q(ihi+i,j) = -q(ihi-i,j)
             else
	       q(ihi+i,j) = -q(ihi-i+1,j)
             endif
	    enddo
	    enddo
         else if (bc(1,2).ne.INT_DIR) then
          print *,"bc invalid"
          stop
	 endif
      endif

      if (nbot .ge. 0) then
         jlo = domlo(2)

         if ((bc(2,1).eq.FOEXTRAP).or. &
             (bc(2,1).eq.EXT_DIR).or. &
             (bc(2,1).eq.HOEXTRAP)) then
	    do j = 1, nbot
	    do i = ARG_L1(q), ARG_H1(q)
	       q(i,jlo-j) = q(i,jlo)
	    enddo
	    enddo
	 elseif (bc(2,1) .eq. REFLECT_EVEN) then
	    do j = 1, nbot
 	    do i = ARG_L1(q), ARG_H1(q)
              if (dir.eq.1) then
	       q(i,jlo-j) = q(i,jlo+j)
              else
	       q(i,jlo-j) = q(i,jlo+j-1)
              endif
	    enddo
	    enddo
	 elseif (bc(2,1) .eq. REFLECT_ODD) then
 	    do i = ARG_L1(q), ARG_H1(q)
              if (dir.eq.1) then
               q(i,jlo)=zero
              endif
	    do j = 1, nbot
              if (dir.eq.1) then
	       q(i,jlo-j) = -q(i,jlo+j)
              else
	       q(i,jlo-j) = -q(i,jlo+j-1)
              endif
	    enddo
	    enddo
         else if (bc(2,1).ne.INT_DIR) then
          print *,"bc invalid"
          stop
	 endif
      endif

      if (ntop .ge. 0) then
         jhi = domhi(2)

         if ((bc(2,2).eq.FOEXTRAP).or. &
             (bc(2,2).eq.EXT_DIR).or. &
             (bc(2,2).eq.HOEXTRAP)) then
	    do j = 1, ntop
 	    do i = ARG_L1(q), ARG_H1(q)
	       q(i,jhi+j) = q(i,jhi)
	    enddo
	    enddo
	 elseif (bc(2,2) .eq. REFLECT_EVEN) then
	    do j = 1, ntop
 	    do i = ARG_L1(q), ARG_H1(q)
              if (dir.eq.1) then
	       q(i,jhi+j) = q(i,jhi-j)
              else
	       q(i,jhi+j) = q(i,jhi-j+1)
              endif
	    enddo
	    enddo
	 elseif (bc(2,2) .eq. REFLECT_ODD) then
 	    do i = ARG_L1(q), ARG_H1(q)
              if (dir.eq.1) then
               q(i,jhi)=zero
              endif
	    do j = 1, ntop
              if (dir.eq.1) then
	       q(i,jhi+j) = -q(i,jhi-j)
              else
	       q(i,jhi+j) = -q(i,jhi-j+1)
              endif
	    enddo
	    enddo
         else if (bc(2,2).ne.INT_DIR) then
          print *,"bc invalid"
          stop
	 endif
      endif

      return
      end


! domlo,domhi are dimensions for node quantity (not cell) 
      subroutine ndfilcc(bfact, &
       q,DIMS(q), &
       domlo,domhi,bc)
      IMPLICIT NONE

      INTEGER_T bfact
      INTEGER_T    DIMDEC(q)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      INTEGER_T    bc(SDIM,2)
      REAL_T       q(DIMV(q))

      INTEGER_T    nlft, nrgt, nbot, ntop, nup, ndwn
      INTEGER_T    ilo, ihi, jlo, jhi, klo, khi
      INTEGER_T    is,  ie,  js,  je,  ks,  ke
      INTEGER_T    i, j, k
      INTEGER_T    kfirst,klast

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      is = max(ARG_L1(q),domlo(1))
      ie = min(ARG_H1(q),domhi(1))
      js = max(ARG_L2(q),domlo(2))
      je = min(ARG_H2(q),domhi(2))

#if (BL_SPACEDIM==3)
      ks = max(ARG_L3(q),domlo(SDIM))
      ke = min(ARG_H3(q),domhi(SDIM))
      ndwn = max(0,domlo(SDIM)-ARG_L3(q))
      nup  = max(0,ARG_H3(q)-domhi(SDIM))
      kfirst=ARG_L3(q)
      klast=ARG_H3(q)
#elif (BL_SPACEDIM==2)
      ks=0
      ke=0
      kfirst=0
      klast=0
#else
      print *,"dimension bust"
      stop
#endif

      nlft = max(0,domlo(1)-ARG_L1(q))
      nrgt = max(0,ARG_H1(q)-domhi(1))
      nbot = max(0,domlo(2)-ARG_L2(q))
      ntop = max(0,ARG_H2(q)-domhi(2))

!     ::::: first fill sides
      if (nlft .ge. 0) then
         ilo = domlo(1)

         if ((bc(1,1).eq.FOEXTRAP).or. &
             (bc(1,1).eq.EXT_DIR).or. &
             (bc(1,1).eq.HOEXTRAP)) then
	    do i = 1, nlft
	    do k = kfirst,klast
	    do j = ARG_L2(q),ARG_H2(q)
	       q(D_DECL(ilo-i,j,k)) = q(D_DECL(ilo,j,k))
	    enddo
	    enddo
	    enddo
	 elseif (bc(1,1) .eq. REFLECT_EVEN) then
	    do i = 1, nlft
	    do k = kfirst,klast
	    do j = ARG_L2(q),ARG_H2(q)
	     q(D_DECL(ilo-i,j,k)) = q(D_DECL(ilo+i,j,k))
	    enddo
	    enddo
	    enddo
	 elseif (bc(1,1) .eq. REFLECT_ODD) then
	    do k = kfirst,klast
	    do j = ARG_L2(q),ARG_H2(q)
               q(D_DECL(ilo,j,k))=zero
	    do i = 1, nlft
	     q(D_DECL(ilo-i,j,k)) = -q(D_DECL(ilo+i,j,k))
	    enddo
	    enddo
	    enddo
         else if (bc(1,1).ne.INT_DIR) then
          print *,"bc invalid"
          stop
	 endif
      endif

      if (nrgt .ge. 0) then
         ihi = domhi(1)

         if ((bc(1,2).eq.FOEXTRAP).or. &
             (bc(1,2).eq.EXT_DIR).or. &
             (bc(1,2).eq.HOEXTRAP)) then
	    do i = 1, nrgt
	    do k = kfirst,klast
	    do j = ARG_L2(q),ARG_H2(q)
	       q(D_DECL(ihi+i,j,k)) = q(D_DECL(ihi,j,k))
	    enddo
	    enddo
	    enddo
	 elseif (bc(1,2) .eq. REFLECT_EVEN) then
	    do i = 1, nrgt
	    do k = kfirst,klast
	    do j = ARG_L2(q),ARG_H2(q)
	     q(D_DECL(ihi+i,j,k)) = q(D_DECL(ihi-i,j,k))
	    enddo
	    enddo
	    enddo
	 elseif (bc(1,2) .eq. REFLECT_ODD) then
	    do k = kfirst,klast
	    do j = ARG_L2(q),ARG_H2(q)
              q(D_DECL(ihi,j,k))=zero
	    do i = 1, nrgt
	     q(D_DECL(ihi+i,j,k)) = -q(D_DECL(ihi-i,j,k))
	    enddo
	    enddo
	    enddo
         else if (bc(1,2).ne.INT_DIR) then
          print *,"bc invalid"
          stop
	 endif
      endif

      if (nbot .ge. 0) then
         jlo = domlo(2)
	
         if ((bc(2,1).eq.FOEXTRAP).or. &
             (bc(2,1).eq.EXT_DIR).or. &
             (bc(2,1).eq.HOEXTRAP)) then
	    do j = 1, nbot
	    do k = kfirst,klast
	    do i = ARG_L1(q),ARG_H1(q)
	       q(D_DECL(i,jlo-j,k)) = q(D_DECL(i,jlo,k))
	    enddo
	    enddo
	    enddo
	 elseif (bc(2,1) .eq. REFLECT_EVEN) then
	    do j = 1, nbot 
	    do k = kfirst,klast
	    do i = ARG_L1(q),ARG_H1(q)
	     q(D_DECL(i,jlo-j,k)) = q(D_DECL(i,jlo+j,k))
	    enddo
	    enddo
	    enddo
	 elseif (bc(2,1) .eq. REFLECT_ODD) then
	    do k = kfirst,klast
	    do i = ARG_L1(q),ARG_H1(q)
              q(D_DECL(i,jlo,k))=zero
	    do j = 1, nbot
	     q(D_DECL(i,jlo-j,k)) = -q(D_DECL(i,jlo+j,k))
	    enddo
	    enddo
	    enddo
         else if (bc(2,1).ne.INT_DIR) then
          print *,"bc invalid"
          stop
	 endif
      endif

      if (ntop .ge. 0) then
         jhi = domhi(2)

         if ((bc(2,2).eq.FOEXTRAP).or. &
             (bc(2,2).eq.EXT_DIR).or. &
             (bc(2,2).eq.HOEXTRAP)) then
	    do j = 1, ntop
	    do k = kfirst,klast
	    do i = ARG_L1(q),ARG_H1(q)
	       q(D_DECL(i,jhi+j,k)) = q(D_DECL(i,jhi,k))
	    enddo
	    enddo
	    enddo
	 elseif (bc(2,2) .eq. REFLECT_EVEN) then
	    do j = 1, ntop
	    do k = kfirst,klast
	    do i = ARG_L1(q),ARG_H1(q)
	     q(D_DECL(i,jhi+j,k)) = q(D_DECL(i,jhi-j,k))
	    enddo
	    enddo
	    enddo
	 elseif (bc(2,2) .eq. REFLECT_ODD) then
	    do k = kfirst,klast
	    do i = ARG_L1(q),ARG_H1(q)
              q(D_DECL(i,jhi,k))=zero
	    do j = 1, ntop
	     q(D_DECL(i,jhi+j,k)) = -q(D_DECL(i,jhi-j,k))
	    enddo
	    enddo
	    enddo
         else if (bc(2,2).ne.INT_DIR) then
          print *,"bc invalid"
          stop
	 endif
      endif

      if (SDIM.eq.3) then

       if (ndwn .ge. 0) then
         klo = domlo(SDIM)

         if ((bc(SDIM,1).eq.FOEXTRAP).or. &
             (bc(SDIM,1).eq.EXT_DIR).or. &
             (bc(SDIM,1).eq.HOEXTRAP)) then
	    do k = 1, ndwn
	    do j = ARG_L2(q),ARG_H2(q)
	    do i = ARG_L1(q),ARG_H1(q)
	       q(D_DECL(i,j,klo-k)) = q(D_DECL(i,j,klo))
	    enddo
	    enddo
	    enddo
	 elseif (bc(SDIM,1) .eq. REFLECT_EVEN) then
	    do k = 1, ndwn
	    do j = ARG_L2(q),ARG_H2(q)
	    do i = ARG_L1(q),ARG_H1(q)
	     q(D_DECL(i,j,klo-k)) = q(D_DECL(i,j,klo+k))
	    enddo
	    enddo
	    enddo
	 elseif (bc(SDIM,1) .eq. REFLECT_ODD) then
	    do j = ARG_L2(q),ARG_H2(q)
	    do i = ARG_L1(q),ARG_H1(q)
              q(D_DECL(i,j,klo))=zero
	    do k = 1, ndwn
	     q(D_DECL(i,j,klo-k)) = -q(D_DECL(i,j,klo+k))
	    enddo
	    enddo
	    enddo
         else if (bc(SDIM,1).ne.INT_DIR) then
          print *,"bc invalid"
          stop
	 endif
       endif

       if (nup .ge. 0) then
         khi = domhi(SDIM)

         if ((bc(SDIM,2).eq.FOEXTRAP).or. &
             (bc(SDIM,2).eq.EXT_DIR).or. &
             (bc(SDIM,2).eq.HOEXTRAP)) then
	    do k = 1, nup
	    do j = ARG_L2(q),ARG_H2(q)
	    do i = ARG_L1(q),ARG_H1(q)
	       q(D_DECL(i,j,khi+k)) = q(D_DECL(i,j,khi))
	    enddo
	    enddo
	    enddo
	 elseif (bc(SDIM,2) .eq. REFLECT_EVEN) then
	    do j = ARG_L2(q),ARG_H2(q)
	    do i = ARG_L1(q),ARG_H1(q)
	    do k = 1, nup
	     q(D_DECL(i,j,khi+k)) = q(D_DECL(i,j,khi-k))
	    enddo
	    enddo
	    enddo
	 elseif (bc(SDIM,2) .eq. REFLECT_ODD) then
	    do j = ARG_L2(q),ARG_H2(q)
	    do i = ARG_L1(q),ARG_H1(q)
              q(D_DECL(i,j,khi))=zero
	    do k = 1, nup
	     q(D_DECL(i,j,khi+k)) = -q(D_DECL(i,j,khi-k))
	    enddo
	    enddo
	    enddo
         else if (bc(SDIM,2).ne.INT_DIR) then
          print *,"bc invalid"
          stop
	 endif
       endif

      else if (SDIM.eq.2) then
       ! do nothing
      else
       print *,"dimension bust"
       stop
      endif


      return
      end



