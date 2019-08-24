!
! $Id: LO_UTIL.F,v 1.1 1998/03/24 07:07:23 almgren Exp $
!
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif
      
#undef TESTING_POLYNOMINTERPCOEFFS
#ifdef TESTING_POLYNOMINTERPCOEFFS
#define NORDER 3
#define REAL_T real*8
#define zero 0.d0
#define one 1.d0
#else
#include <CONSTANTS.H>
#include <REAL.H>
#endif
    
      

!     polyInterpCoeff:
!  
!     This routine returns the Lagrange interpolating coefficients for a
!     polynomial through N points, evaluated at xInt (see Numerical Recipes,
!     v2, p102, e.g.):
!
!            (x-x2)(x-x3)...(x-xN)              (x-x1)(x-x2)...(x-x(N-1))
!    P(x) = ----------------------- y1  + ... + ------------------------  yN
!           (x1-x2)(x1-x3)...(x1-xN)            (x1-x2)(x1-x3)...(x1-xN)
!
!     P(xInt) = sum_(i=1)^(N) y[i]*c[i]
!
      subroutine polyInterpCoeff(xInt, x, N, c)
      implicit none
      INTEGER_T N, i, j
      REAL_T xInt, x(N), c(N), num, den
      do j=1,N
         num = one
         den = one
         do i = 1,j-1
            num = num*(xInt - x(i))
            den = den*(x(j) - x(i))
         end do
         do i = j+1,N
            num = num*(xInt - x(i))
            den = den*(x(j) - x(i))
         end do
         if (den .eq. zero) then
          print *,"polyInterpCoeff::invalid data"
          stop
         endif
         c(j) = num/den
      end do
      return
      end

      subroutine polyInterpCoeffDeriv(xInt, x, N, c)
      implicit none
      INTEGER_T N, i, j, k
      REAL_T xInt, x(N), c(N), num, num2, den
      do j=1,N
         den = one
         do i = 1,j-1
            den = den*(x(j) - x(i))
         end do
         do i = j+1,N
            den = den*(x(j) - x(i))
         end do
         if (den .eq. zero) then
          print *,"polyInterpCoeffDeriv::invalid data"
          stop
         endif

         num=zero
         do k=1,N
          if (k.ne.j) then
           num2=one
           do i=1,N
            if ((i.ne.k).and.(i.ne.j)) then
             num2=num2*(xInt-x(i))
            endif
           enddo
           num=num+num2
          endif
         enddo
         c(j) = num/den
      end do
      return
      end


      

#ifdef TESTING_POLYNOMINTERPCOEFFS

!     
!     This is a test driver for the routine polyInterpCoeff.  Sample data
!     is created from the statement function, and the location of the 
!     boundary node and internal nodes are set, as apporpriate.  The
!     number of points created is equal to the test NORDER set at the
!     top of this file through a define.  The coefficients are computed,
!     and then the ghost cell value is constructed from the resulting
!     coefficients and written out.
!

      program polyInterpCoeffTest
      INTEGER_T i, j
      REAL_T c(NORDER), y(NORDER), yInt
      REAL_T x(NORDER), const, yfunc, xx, xInt
!
!     Make data, locations
      yfunc(xx) = (2*xx-one)**(NORDER-1) + const
      const = -one
      x(1) = zero
      y(1) = yfunc(x(1))
      do j=2,NORDER
         x(j) = j - 1.5*one
         y(j) = yfunc(x(j))
      end do
!      
!     Set interpolation point to ghost cell location
      xInt = (x(2) - x(3))*0.5
!      
!     Call routine
      call polyInterpCoeff(xInt,x,NORDER,c)
!
!     Evaluate result
      yInt = zero
      do j=1,NORDER
         yInt = yInt + c(j)*y(j)
      end do
!
!     Dump output
      write(6,*) 'x = ',x
      write(6,*) 'y = ',y
      write(6,*) 'c = ',c
      write(6,*) 'Interpolated y = ',yInt
      end

#endif
