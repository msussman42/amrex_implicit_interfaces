#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_FORT_INTEGER.H"
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

      module filcc_module
      contains

      subroutine local_filcc(bfact, &
       q, &
       domlo,domhi,bc)
      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: domlo(SDIM), domhi(SDIM)
       ! q inherits attributes from the target.
      REAL_T, INTENT(in), pointer :: q(D_DECL(:,:,:))
      INTEGER_T, INTENT(in) :: bc(SDIM,2)

      INTEGER_T, INTENT(in) :: bfact

      INTEGER_T    nlft, nrgt, nbot, ntop, nup, ndwn
      INTEGER_T    ilo, ihi, jlo, jhi
      INTEGER_T    i, j
#if (AMREX_SPACEDIM==3)
      INTEGER_T    k,klo,khi
#endif

      if (bfact.lt.1) then
       print *,"bfact invalid710"
       stop
      endif

      nlft = max(0,domlo(1)-LBOUND(q,1))
      nrgt = max(0,UBOUND(q,1)-domhi(1))
      nbot = max(0,domlo(2)-LBOUND(q,2))
      ntop = max(0,UBOUND(q,2)-domhi(2))
      ndwn=0
      nup=0
#if (AMREX_SPACEDIM==3)
      ndwn = max(0,domlo(SDIM)-LBOUND(q,SDIM))
      nup  = max(0,UBOUND(q,SDIM)-domhi(SDIM))
#elif (AMREX_SPACEDIM==2)
      ! do nothing
#else  
print *,"dimension bust"
stop
#endif

!
!     ::::: first fill sides
!
      if (nlft .gt. 0) then
       ilo = domlo(1)

       if ((bc(1,1).eq.FOEXTRAP).or. &
           (bc(1,1).eq.EXT_DIR)) then
        do i = 1, nlft
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do j = LBOUND(q,2),UBOUND(q,2)
           q(D_DECL(ilo-i,j,k)) = q(D_DECL(ilo,j,k))
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(1,1) .eq. HOEXTRAP) then
        do i = 1, nlft
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do j = LBOUND(q,2),UBOUND(q,2)
           q(D_DECL(ilo-i,j,k)) = &
             two*q(D_DECL(ilo-i+1,j,k))-q(D_DECL(ilo-i+2,j,k))
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(1,1) .eq. REFLECT_EVEN) then
        do i = 1, nlft
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do j = LBOUND(q,2),UBOUND(q,2)
           q(D_DECL(ilo-i,j,k)) = q(D_DECL(ilo+i-1,j,k))
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(1,1) .eq. REFLECT_ODD) then
        do i = 1, nlft
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do j = LBOUND(q,2),UBOUND(q,2)
           q(D_DECL(ilo-i,j,k)) = -q(D_DECL(ilo+i-1,j,k))
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
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
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do j = LBOUND(q,2),UBOUND(q,2)
           q(D_DECL(ihi+i,j,k)) = q(D_DECL(ihi,j,k))
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(1,2) .eq. HOEXTRAP) then
        do i = 1, nrgt
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do j = LBOUND(q,2),UBOUND(q,2)
           q(D_DECL(ihi+i,j,k)) = &
            two*q(D_DECL(ihi+i-1,j,k))-q(D_DECL(ihi+i-2,j,k))
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(1,2) .eq. REFLECT_EVEN) then
        do i = 1, nrgt
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do j = LBOUND(q,2),UBOUND(q,2)
           q(D_DECL(ihi+i,j,k)) = q(D_DECL(ihi-i+1,j,k))
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(1,2) .eq. REFLECT_ODD) then
        do i = 1, nrgt
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do j = LBOUND(q,2),UBOUND(q,2)
           q(D_DECL(ihi+i,j,k)) = -q(D_DECL(ihi-i+1,j,k))
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
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
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do i = LBOUND(q,1),UBOUND(q,1)
           q(D_DECL(i,jlo-j,k)) = q(D_DECL(i,jlo,k))
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(2,1) .eq. HOEXTRAP) then
        do j = 1, nbot
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do i = LBOUND(q,1),UBOUND(q,1)
           q(D_DECL(i,jlo-j,k)) = &
             two*q(D_DECL(i,jlo-j+1,k))-q(D_DECL(i,jlo-j+2,k))
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(2,1) .eq. REFLECT_EVEN) then
        do j = 1, nbot 
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do i = LBOUND(q,1),UBOUND(q,1)
           q(D_DECL(i,jlo-j,k)) = q(D_DECL(i,jlo+j-1,k))
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(2,1) .eq. REFLECT_ODD) then
        do j = 1, nbot
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do i = LBOUND(q,1),UBOUND(q,1)
           q(D_DECL(i,jlo-j,k)) = -q(D_DECL(i,jlo+j-1,k))
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
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
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do i = LBOUND(q,1),UBOUND(q,1)
           q(D_DECL(i,jhi+j,k)) = q(D_DECL(i,jhi,k))
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(2,2) .eq. HOEXTRAP) then
        do j = 1, ntop
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do i = LBOUND(q,1),UBOUND(q,1)
           q(D_DECL(i,jhi+j,k)) = &
             two*q(D_DECL(i,jhi+j-1,k))-q(D_DECL(i,jhi+j-2,k))
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(2,2) .eq. REFLECT_EVEN) then
        do j = 1, ntop
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do i = LBOUND(q,1),UBOUND(q,1)
           q(D_DECL(i,jhi+j,k)) = q(D_DECL(i,jhi-j+1,k))
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(2,2) .eq. REFLECT_ODD) then
        do j = 1, ntop
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do i = LBOUND(q,1),UBOUND(q,1)
           q(D_DECL(i,jhi+j,k)) = -q(D_DECL(i,jhi-j+1,k))
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(2,2).eq.INT_DIR) then
        ! do nothing
       else
        print *,"bc invalid"
        stop
       end if
      end if

#if (AMREX_SPACEDIM==3)
      if (ndwn .gt. 0) then
       klo = domlo(SDIM)

       if ((bc(SDIM,1).eq.FOEXTRAP).or. &
           (bc(SDIM,1).eq.EXT_DIR)) then
        do k = 1, ndwn
        do j = LBOUND(q,2),UBOUND(q,2)
        do i = LBOUND(q,1),UBOUND(q,1)
         q(i,j,klo-k) = q(i,j,klo)
        end do
        end do
        end do
       else if (bc(SDIM,1) .eq. HOEXTRAP) then
        do k = 1, ndwn
        do j = LBOUND(q,2),UBOUND(q,2)
        do i = LBOUND(q,1),UBOUND(q,1)
         q(i,j,klo-k) = two*q(i,j,klo-k+1)-q(i,j,klo-k+2)
        end do
        end do
        end do
       else if (bc(SDIM,1) .eq. REFLECT_EVEN) then
        do k = 1, ndwn
        do j = LBOUND(q,2),UBOUND(q,2)
        do i = LBOUND(q,1),UBOUND(q,1)
         q(i,j,klo-k) = q(i,j,klo+k-1)
        end do
        end do
        end do
       else if (bc(SDIM,1) .eq. REFLECT_ODD) then
        do k = 1, ndwn
        do j = LBOUND(q,2),UBOUND(q,2)
        do i = LBOUND(q,1),UBOUND(q,1)
         q(i,j,klo-k) = -q(i,j,klo+k-1)
        end do
        end do
        end do
       else if (bc(SDIM,1).eq.INT_DIR) then
        ! do nothing
       else
        print *,"bc invalid"
        stop
       end if
      end if

      if (nup .gt. 0) then
       khi = domhi(SDIM)

       if ((bc(SDIM,2).eq.FOEXTRAP).or. &
           (bc(SDIM,2).eq.EXT_DIR)) then
        do k = 1, nup
        do j = LBOUND(q,2),UBOUND(q,2)
        do i = LBOUND(q,1),UBOUND(q,1)
         q(i,j,khi+k) = q(i,j,khi)
        end do
        end do
        end do
       else if (bc(SDIM,2) .eq. HOEXTRAP) then
        do k = 1, nup
        do j = LBOUND(q,2),UBOUND(q,2)
        do i = LBOUND(q,1),UBOUND(q,1)
         q(i,j,khi+k) = two*q(i,j,khi+k-1)-q(i,j,khi+k-2)
        end do
        end do
        end do
       else if (bc(SDIM,2) .eq. REFLECT_EVEN) then
        do k = 1, nup
        do j = LBOUND(q,2),UBOUND(q,2)
        do i = LBOUND(q,1),UBOUND(q,1)
         q(i,j,khi+k) = q(i,j,khi-k+1)
        end do
        end do
        end do
       else if (bc(SDIM,2) .eq. REFLECT_ODD) then
        do k = 1, nup
        do j = LBOUND(q,2),UBOUND(q,2)
        do i = LBOUND(q,1),UBOUND(q,1)
         q(i,j,khi+k) = -q(i,j,khi-k+1)
        end do
        end do
        end do
       else if (bc(SDIM,2).eq.INT_DIR) then
        ! do nothing
       else
        print *,"bc invalid"
        stop
       end if
      end if
#endif

      end subroutine local_filcc

      subroutine local_filcc4D( &
       bfact, &
       q,scomp, &
       domlo,domhi,bc)
      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: scomp
      INTEGER_T, INTENT(in) :: domlo(SDIM), domhi(SDIM)
       ! q inherits attributes from the target.
      REAL_T, INTENT(in), pointer :: q(D_DECL(:,:,:),:)
      INTEGER_T, INTENT(in) :: bc(SDIM,2,scomp)

      INTEGER_T, INTENT(in) :: bfact

      INTEGER_T    nlft, nrgt, nbot, ntop, nup, ndwn
      INTEGER_T    ilo, ihi, jlo, jhi
      INTEGER_T    i, j
#if (AMREX_SPACEDIM==3)
      INTEGER_T    k,klo,khi
#endif

      if (bfact.lt.1) then
       print *,"bfact invalid710"
       stop
      endif

      nlft = max(0,domlo(1)-LBOUND(q,1))
      nrgt = max(0,UBOUND(q,1)-domhi(1))
      nbot = max(0,domlo(2)-LBOUND(q,2))
      ntop = max(0,UBOUND(q,2)-domhi(2))
      ndwn=0
      nup=0
#if (AMREX_SPACEDIM==3)
      ndwn = max(0,domlo(SDIM)-LBOUND(q,SDIM))
      nup  = max(0,UBOUND(q,SDIM)-domhi(SDIM))
#elif (AMREX_SPACEDIM==2)
      ! do nothing
#else  
print *,"dimension bust"
stop
#endif

!
!     ::::: first fill sides
!
      if (nlft .gt. 0) then
       ilo = domlo(1)

       if ((bc(1,1,scomp).eq.FOEXTRAP).or. &
           (bc(1,1,scomp).eq.EXT_DIR)) then
        do i = 1, nlft
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do j = LBOUND(q,2),UBOUND(q,2)
           q(D_DECL(ilo-i,j,k),scomp) = q(D_DECL(ilo,j,k),scomp)
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(1,1,scomp) .eq. HOEXTRAP) then
        do i = 1, nlft
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do j = LBOUND(q,2),UBOUND(q,2)
           q(D_DECL(ilo-i,j,k),scomp) = &
             two*q(D_DECL(ilo-i+1,j,k),scomp)-q(D_DECL(ilo-i+2,j,k),scomp)
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(1,1,scomp) .eq. REFLECT_EVEN) then
        do i = 1, nlft
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do j = LBOUND(q,2),UBOUND(q,2)
           q(D_DECL(ilo-i,j,k),scomp) = q(D_DECL(ilo+i-1,j,k),scomp)
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(1,1,scomp) .eq. REFLECT_ODD) then
        do i = 1, nlft
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do j = LBOUND(q,2),UBOUND(q,2)
           q(D_DECL(ilo-i,j,k),scomp) = -q(D_DECL(ilo+i-1,j,k),scomp)
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(1,1,scomp).eq.INT_DIR) then
        ! do nothing
       else
        print *,"bc invalid"
        stop
       end if
      end if

      if (nrgt .gt. 0) then
       ihi = domhi(1)

       if ((bc(1,2,scomp).eq.FOEXTRAP).or. &
           (bc(1,2,scomp).eq.EXT_DIR)) then
        do i = 1, nrgt
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do j = LBOUND(q,2),UBOUND(q,2)
           q(D_DECL(ihi+i,j,k),scomp) = q(D_DECL(ihi,j,k),scomp)
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(1,2,scomp) .eq. HOEXTRAP) then
        do i = 1, nrgt
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do j = LBOUND(q,2),UBOUND(q,2)
           q(D_DECL(ihi+i,j,k),scomp) = &
            two*q(D_DECL(ihi+i-1,j,k),scomp)-q(D_DECL(ihi+i-2,j,k),scomp)
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(1,2,scomp) .eq. REFLECT_EVEN) then
        do i = 1, nrgt
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do j = LBOUND(q,2),UBOUND(q,2)
           q(D_DECL(ihi+i,j,k),scomp) = q(D_DECL(ihi-i+1,j,k),scomp)
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(1,2,scomp) .eq. REFLECT_ODD) then
        do i = 1, nrgt
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do j = LBOUND(q,2),UBOUND(q,2)
           q(D_DECL(ihi+i,j,k),scomp) = -q(D_DECL(ihi-i+1,j,k),scomp)
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(1,2,scomp).eq.INT_DIR) then
        ! do nothing
       else
        print *,"bc invalid"
        stop
       end if
      end if

      if (nbot .gt. 0) then
       jlo = domlo(2)
         
       if ((bc(2,1,scomp).eq.FOEXTRAP).or. &
           (bc(2,1,scomp).eq.EXT_DIR)) then
        do j = 1, nbot
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do i = LBOUND(q,1),UBOUND(q,1)
           q(D_DECL(i,jlo-j,k),scomp) = q(D_DECL(i,jlo,k),scomp)
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(2,1,scomp) .eq. HOEXTRAP) then
        do j = 1, nbot
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do i = LBOUND(q,1),UBOUND(q,1)
           q(D_DECL(i,jlo-j,k),scomp) = &
             two*q(D_DECL(i,jlo-j+1,k),scomp)-q(D_DECL(i,jlo-j+2,k),scomp)
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(2,1,scomp) .eq. REFLECT_EVEN) then
        do j = 1, nbot 
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do i = LBOUND(q,1),UBOUND(q,1)
           q(D_DECL(i,jlo-j,k),scomp) = q(D_DECL(i,jlo+j-1,k),scomp)
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(2,1,scomp) .eq. REFLECT_ODD) then
        do j = 1, nbot
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do i = LBOUND(q,1),UBOUND(q,1)
           q(D_DECL(i,jlo-j,k),scomp) = -q(D_DECL(i,jlo+j-1,k),scomp)
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(2,1,scomp).eq.INT_DIR) then
        ! do nothing
       else
        print *,"bc invalid"
        stop
       end if
      end if

      if (ntop .gt. 0) then
       jhi = domhi(2)

       if ((bc(2,2,scomp).eq.FOEXTRAP).or. &
           (bc(2,2,scomp).eq.EXT_DIR)) then
        do j = 1, ntop
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do i = LBOUND(q,1),UBOUND(q,1)
           q(D_DECL(i,jhi+j,k),scomp) = q(D_DECL(i,jhi,k),scomp)
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(2,2,scomp) .eq. HOEXTRAP) then
        do j = 1, ntop
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do i = LBOUND(q,1),UBOUND(q,1)
           q(D_DECL(i,jhi+j,k),scomp) = &
             two*q(D_DECL(i,jhi+j-1,k),scomp)-q(D_DECL(i,jhi+j-2,k),scomp)
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(2,2,scomp) .eq. REFLECT_EVEN) then
        do j = 1, ntop
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do i = LBOUND(q,1),UBOUND(q,1)
           q(D_DECL(i,jhi+j,k),scomp) = q(D_DECL(i,jhi-j+1,k),scomp)
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(2,2,scomp) .eq. REFLECT_ODD) then
        do j = 1, ntop
#if (AMREX_SPACEDIM==3)
         do k = LBOUND(q,SDIM),UBOUND(q,SDIM)
#endif
          do i = LBOUND(q,1),UBOUND(q,1)
           q(D_DECL(i,jhi+j,k),scomp) = -q(D_DECL(i,jhi-j+1,k),scomp)
          end do
#if (AMREX_SPACEDIM==3)
         end do
#endif
        end do
       else if (bc(2,2,scomp).eq.INT_DIR) then
        ! do nothing
       else
        print *,"bc invalid"
        stop
       end if
      end if

#if (AMREX_SPACEDIM==3)
      if (ndwn .gt. 0) then
       klo = domlo(SDIM)

       if ((bc(SDIM,1,scomp).eq.FOEXTRAP).or. &
           (bc(SDIM,1,scomp).eq.EXT_DIR)) then
        do k = 1, ndwn
        do j = LBOUND(q,2),UBOUND(q,2)
        do i = LBOUND(q,1),UBOUND(q,1)
         q(i,j,klo-k,scomp) = q(i,j,klo,scomp)
        end do
        end do
        end do
       else if (bc(SDIM,1,scomp) .eq. HOEXTRAP) then
        do k = 1, ndwn
        do j = LBOUND(q,2),UBOUND(q,2)
        do i = LBOUND(q,1),UBOUND(q,1)
         q(i,j,klo-k,scomp) = two*q(i,j,klo-k+1,scomp)-q(i,j,klo-k+2,scomp)
        end do
        end do
        end do
       else if (bc(SDIM,1,scomp) .eq. REFLECT_EVEN) then
        do k = 1, ndwn
        do j = LBOUND(q,2),UBOUND(q,2)
        do i = LBOUND(q,1),UBOUND(q,1)
         q(i,j,klo-k,scomp) = q(i,j,klo+k-1,scomp)
        end do
        end do
        end do
       else if (bc(SDIM,1,scomp) .eq. REFLECT_ODD) then
        do k = 1, ndwn
        do j = LBOUND(q,2),UBOUND(q,2)
        do i = LBOUND(q,1),UBOUND(q,1)
         q(i,j,klo-k,scomp) = -q(i,j,klo+k-1,scomp)
        end do
        end do
        end do
       else if (bc(SDIM,1,scomp).eq.INT_DIR) then
        ! do nothing
       else
        print *,"bc invalid"
        stop
       end if
      end if

      if (nup .gt. 0) then
       khi = domhi(SDIM)

       if ((bc(SDIM,2,scomp).eq.FOEXTRAP).or. &
           (bc(SDIM,2,scomp).eq.EXT_DIR)) then
        do k = 1, nup
        do j = LBOUND(q,2),UBOUND(q,2)
        do i = LBOUND(q,1),UBOUND(q,1)
         q(i,j,khi+k,scomp) = q(i,j,khi,scomp)
        end do
        end do
        end do
       else if (bc(SDIM,2,scomp) .eq. HOEXTRAP) then
        do k = 1, nup
        do j = LBOUND(q,2),UBOUND(q,2)
        do i = LBOUND(q,1),UBOUND(q,1)
         q(i,j,khi+k,scomp) = two*q(i,j,khi+k-1,scomp)-q(i,j,khi+k-2,scomp)
        end do
        end do
        end do
       else if (bc(SDIM,2,scomp) .eq. REFLECT_EVEN) then
        do k = 1, nup
        do j = LBOUND(q,2),UBOUND(q,2)
        do i = LBOUND(q,1),UBOUND(q,1)
         q(i,j,khi+k,scomp) = q(i,j,khi-k+1,scomp)
        end do
        end do
        end do
       else if (bc(SDIM,2,scomp) .eq. REFLECT_ODD) then
        do k = 1, nup
        do j = LBOUND(q,2),UBOUND(q,2)
        do i = LBOUND(q,1),UBOUND(q,1)
         q(i,j,khi+k,scomp) = -q(i,j,khi-k+1,scomp)
        end do
        end do
        end do
       else if (bc(SDIM,2,scomp).eq.INT_DIR) then
        ! do nothing
       else
        print *,"bc invalid"
        stop
       end if
      end if
#endif

      end subroutine local_filcc4D


       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      subroutine local_grid_type_to_box_type(grid_type,box_type)
      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: grid_type
      INTEGER_T, INTENT(out) :: box_type(SDIM)
      INTEGER_T dir

      do dir=1,SDIM
       box_type(dir)=0  ! default to CELL
      enddo
      if (grid_type.eq.-1) then
       ! do nothing
      else if ((grid_type.ge.0).and. &
               (grid_type.lt.SDIM)) then
       box_type(grid_type+1)=1  ! NODE
      else if (grid_type.eq.3) then
       box_type(1)=1 ! NODE
       box_type(2)=1 ! NODE
      else if ((grid_type.eq.4).and.(SDIM.eq.3)) then
       box_type(1)=1 ! NODE
       box_type(SDIM)=1 ! NODE
      else if ((grid_type.eq.5).and.(SDIM.eq.3)) then
       box_type(2)=1 ! NODE
       box_type(SDIM)=1 ! NODE
      else
       print *,"grid_type invalid"
       stop
      endif
     
      return 
      end subroutine local_grid_type_to_box_type

      subroutine check_arr_idx(i,j,k,fablo,fabhi)
      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: i,j,k
      INTEGER_T, INTENT(in) :: fablo(3),fabhi(3)
      INTEGER_T :: ii(3)
      INTEGER_T :: dir
     
      ii(1)=i 
      ii(2)=j 
      ii(3)=k
      do dir=1,SDIM
       if ((ii(dir).ge.fablo(dir)).and. &
           (ii(dir).le.fabhi(dir))) then
        ! do nothing
       else
        print *,"bust:dir,ii,fablo,fabhi ",dir,ii(dir),fablo(dir),fabhi(dir)
        stop
       endif
      enddo ! dir=1..sdim
        
      return
      end subroutine check_arr_idx
 
! domlo,domhi are dimensions for face quantity (not cell) 
      subroutine efilcc(bfact, &
       q, &
       domlo,domhi,bc,grid_type)
      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: domlo(SDIM), domhi(SDIM)
      REAL_T, INTENT(in), pointer :: q(D_DECL(:,:,:))
      INTEGER_T, INTENT(in) :: bc(SDIM,2)

      INTEGER_T, INTENT(in) :: grid_type
      INTEGER_T, INTENT(in) :: bfact

      INTEGER_T :: box_type(SDIM)
      INTEGER_T :: ntofill(3,2)
      INTEGER_T :: fablo_declare(3)
      INTEGER_T :: fabhi_declare(3)
      INTEGER_T :: int_lo(3)
      INTEGER_T :: int_hi(3)
      INTEGER_T :: local_dir
      INTEGER_T :: dir_side
      INTEGER_T :: side
      INTEGER_T :: sidelo(3),sidehi(3),side_inc(3)
      INTEGER_T :: side_first(3)
      INTEGER_T :: side_last(3)
      INTEGER_T :: imult(3)
      INTEGER_T :: idist(3)

      INTEGER_T :: i,j,k
      INTEGER_T :: isrc,jsrc,ksrc
      INTEGER_T :: idx_norm


      call local_grid_type_to_box_type(grid_type,box_type);

      if (bfact.lt.1) then
       print *,"bfact invalid in EFILCC"
       stop
      endif

      fablo_declare(1)=LBOUND(q,1)
      fablo_declare(2)=LBOUND(q,2)
      fablo_declare(3)=0
#if (AMREX_SPACEDIM==3)
      fablo_declare(3)=LBOUND(q,SDIM)
#endif
      fabhi_declare(1)=UBOUND(q,1)
      fabhi_declare(2)=UBOUND(q,2)
      fabhi_declare(3)=0
#if (AMREX_SPACEDIM==3)
      fabhi_declare(3)=UBOUND(q,SDIM)
#endif
     
      int_lo(3)=0
      int_hi(3)=0
      do local_dir=1,SDIM
       int_lo(local_dir)=max(fablo_declare(local_dir),domlo(local_dir))
       int_hi(local_dir)=min(fabhi_declare(local_dir),domhi(local_dir))
      enddo

      ntofill(3,1)=0
      ntofill(3,2)=0
      do local_dir=1,SDIM
       side=1
       ntofill(local_dir,side)= &
         max(0,domlo(local_dir)-fablo_declare(local_dir))
       side=2
       ntofill(local_dir,side)= &
         max(0,fabhi_declare(local_dir)-domhi(local_dir))
      enddo ! local_dir=1..sdim
 
      do local_dir=1,SDIM

       do side=1,2

        if (ntofill(local_dir,side).ge.0) then

         sidelo(3)=0
         sidehi(3)=0
         side_inc(3)=1
         do dir_side=1,SDIM
          sidelo(dir_side)=fablo_declare(dir_side)
          sidehi(dir_side)=fabhi_declare(dir_side)
          side_inc(dir_side)=1
         enddo

         if (side.eq.1) then
          sidelo(local_dir)=domlo(local_dir)-ntofill(local_dir,side)
          sidehi(local_dir)=domlo(local_dir)
          side_inc(local_dir)=-1
         else if (side.eq.2) then
          sidehi(local_dir)=domhi(local_dir)+ntofill(local_dir,side)
          sidelo(local_dir)=domhi(local_dir)
          side_inc(local_dir)=1
         else
          print *,"side invalid"
          stop
         endif
          
         do dir_side=1,3
          imult(dir_side)=0
          idist(dir_side)=0

          if (sidelo(dir_side).le.sidehi(dir_side)) then
           if (side_inc(dir_side).eq.1) then
            side_first(dir_side)=sidelo(dir_side)
            side_last(dir_side)=sidehi(dir_side)
           else if (side_inc(dir_side).eq.-1) then
            side_first(dir_side)=sidehi(dir_side)
            side_last(dir_side)=sidelo(dir_side)
           else
            print *,"side_inc invalid"
            stop
           endif
          else
           print *,"sidelo or sidehi invalid"
           stop
          endif
         enddo ! dir_side=1..3
           
         imult(local_dir)=-side_inc(local_dir)  

         do i=side_first(1),side_last(1),side_inc(1) 
         do j=side_first(2),side_last(2),side_inc(2) 
         do k=side_first(3),side_last(3),side_inc(3) 
          if (local_dir.eq.1) then
           idx_norm=i
          else if (local_dir.eq.2) then
           idx_norm=j
          else if ((local_dir.eq.SDIM).and.(SDIM.eq.3)) then
           idx_norm=k
          else
           print *,"local_dir invalid"
           stop
          endif
          idist(local_dir)=abs(idx_norm-side_first(local_dir))

          if ((bc(local_dir,side).eq.FOEXTRAP).or. &
              (bc(local_dir,side).eq.EXT_DIR).or. &
              (bc(local_dir,side).eq.HOEXTRAP)) then
           if (idx_norm.ne.side_first(local_dir)) then
            if (idist(local_dir).gt.0) then
             isrc=i+imult(1)
             jsrc=j+imult(2)
             ksrc=k+imult(3)
             call check_arr_idx(i,j,k,fablo_declare,fabhi_declare)
             call check_arr_idx(isrc,jsrc,ksrc,fablo_declare,fabhi_declare)
             q(D_DECL(i,j,k)) = q(D_DECL(isrc,jsrc,ksrc))
            else 
             print *,"idist(local_dir) invalid"
             stop
            endif
           else if (idx_norm.eq.side_first(local_dir)) then
            if (idist(local_dir).eq.0) then
             ! do nothing
            else 
             print *,"idist(local_dir) invalid"
             stop
            endif
           else
            print *,"idx_norm invalid"
            stop
           endif 
          elseif (bc(local_dir,side).eq.REFLECT_EVEN) then
           if (idx_norm.ne.side_first(local_dir)) then
            if (idist(local_dir).gt.0) then
             if (box_type(local_dir).eq.1) then !NODE
              isrc=i+2*idist(1)*imult(1)
              jsrc=j+2*idist(2)*imult(2)
              ksrc=k+2*idist(3)*imult(3)
             else if (box_type(local_dir).eq.0) then !CELL
              isrc=i+(2*idist(1)-1)*imult(1)
              jsrc=j+(2*idist(2)-1)*imult(2)
              ksrc=k+(2*idist(3)-1)*imult(3)
             else
              print *,"box_type(local_dir) invalid"
              stop
             endif
             call check_arr_idx(i,j,k,fablo_declare,fabhi_declare)
             call check_arr_idx(isrc,jsrc,ksrc,fablo_declare,fabhi_declare)
             q(D_DECL(i,j,k)) = q(D_DECL(isrc,jsrc,ksrc))
            else 
             print *,"idist(local_dir) invalid"
             stop
            endif
           else if (idx_norm.eq.side_first(local_dir)) then
            if (idist(local_dir).eq.0) then
             ! do nothing
            else 
             print *,"idist(local_dir) invalid"
             stop
            endif
           else
            print *,"idx_norm invalid"
            stop
           endif 
          elseif (bc(local_dir,side).eq.REFLECT_ODD) then
           if (idx_norm.ne.side_first(local_dir)) then
            if (idist(local_dir).gt.0) then
             if (box_type(local_dir).eq.1) then !NODE
              isrc=i+2*idist(1)*imult(1)
              jsrc=j+2*idist(2)*imult(2)
              ksrc=k+2*idist(3)*imult(3)
             else if (box_type(local_dir).eq.0) then !CELL
              isrc=i+(2*idist(1)-1)*imult(1)
              jsrc=j+(2*idist(2)-1)*imult(2)
              ksrc=k+(2*idist(3)-1)*imult(3)
             else
              print *,"box_type(local_dir) invalid"
              stop
             endif
             call check_arr_idx(i,j,k,fablo_declare,fabhi_declare)
             call check_arr_idx(isrc,jsrc,ksrc,fablo_declare,fabhi_declare)
             q(D_DECL(i,j,k)) = -q(D_DECL(isrc,jsrc,ksrc))
            else 
             print *,"idist(local_dir) invalid"
             stop
            endif
           else if (idx_norm.eq.side_first(local_dir)) then
            if (idist(local_dir).eq.0) then
             if (box_type(local_dir).eq.1) then !NODE
              call check_arr_idx(i,j,k,fablo_declare,fabhi_declare)
              q(D_DECL(i,j,k))=zero
             else if (box_type(local_dir).eq.0) then !CELL
              ! do nothing
             else
              print *,"box_type(local_dir) invalid"
              stop
             endif
            else 
             print *,"idist(local_dir) invalid"
             stop
            endif
           else
            print *,"idx_norm invalid"
            stop
           endif 
          else if (bc(local_dir,side).eq.INT_DIR) then
           ! do nothing
          else
           print *,"bc invalid"
           stop
          endif
         enddo ! k
         enddo ! j
         enddo ! i
        else 
         print *,"expecting ntofill>=0"
         stop
        endif

       enddo ! side=1,2
 
      enddo ! local_dir=1,SDIM

      return
      end subroutine efilcc

      end module filcc_module
