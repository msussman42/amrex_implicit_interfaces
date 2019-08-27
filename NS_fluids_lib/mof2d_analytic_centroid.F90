!This file is part of Notus 0.3.0

!Copyright Bordeaux-INP, Université de Bordeaux, CNRS
!
!Contributors:
!Antoine Lemoine, 14-06-2017, antoine.lemoine@bordeaux-inp.fr

!This software is a computer program whose purpose is to simulate fluid flows.

!This software is governed by the CeCILL license under French law and
!abiding by the rules of distribution of free software.  You can  use,
!modify and/ or redistribute the software under the terms of the CeCILL
!license as circulated by CEA, CNRS and INRIA at the following URL
!"http://www.cecill.info".

!As a counterpart to the access to the source code and  rights to copy,
!modify and redistribute granted by the license, users are provided only
!with a limited warranty  and the software's author,  the holder of the
!economic rights,  and the successive licensors  have only  limited
!liability.

!In this respect, the user's attention is drawn to the risks associated
!with loading,  using,  modifying and/or developing or reproducing the
!software by the user in light of its specific status of free software,
!that may mean  that it is complicated to manipulate,  and  that  also
!therefore means  that it is reserved for developers  and  experienced
!professionals having in-depth computer knowledge. Users are therefore
!encouraged to load and test the software's suitability as regards their
!requirements in conditions enabling the security of their systems and/or
!data to be ensured and,  more generally, to use and operate it in the
!same conditions as regards security.

!The fact that you are presently reading this means that you have had
!knowledge of the CeCILL license and that you accept its terms.

#define DEBUG_ANALYTICAL_GRAD 0

module mod_mof2d_analytic_centroid
   implicit none

   private

   double precision, parameter :: PI_2 = 2d0*atan(1.0d0)
   double precision, parameter :: PI = 4d0*atan(1.0d0)

   public :: mof2d_compute_analytic_gradient

contains

   !> Compute the centroid and the gradient of the objective function in rectangular hexahedral cell.
   !!
   !! @param[in]  angles:        Spherical angles.
   !! @param[in]  ref_centroid:  Coordinates of the reference centroid.
   !! @param[in]  ref_volume:    Reference volume.
   !! @param[in]  c:             Dimensions of the cell.
   !! @param[out] centroids:     Coordinates of the centroid.
   !! @param[out] gradient:      Gradient of the objective function.
   !! @ingroup moment_of_fluid
      subroutine mof2d_compute_analytic_gradient(angles, &
        volume, c, centroid)
      double precision, dimension(1), intent(in) :: angles
      double precision, intent(in) :: volume
      double precision, dimension(2), intent(in) :: c
      double precision, dimension(2), intent(out) :: centroid

      double precision :: n1,n2
      integer :: trap_flag,tri_flag
      integer :: swap_x,swap_y
      double precision :: volcell,volfrac
      double precision :: vol_rect,vol_tri,volume_test
      double precision :: a,b,d

#if (DEBUG_ANALYTICAL_GRAD==1)
      print *,"angles ",angles(1)
      print *,"ref_centroid ",ref_centroid(1), &
              ref_centroid(2)
      print *,"volume ",volume
      print *,"c ",c(1),c(2)
#endif 

      n1=cos(angles(1))
      n2=sin(angles(1))
      trap_flag=0
      tri_flag=0
      volcell=c(1)*c(2)
      volfrac=volume/volcell

      if ((volfrac.gt.0.0d0).and. &
          (volfrac.le.0.5d0)) then

       swap_x=0
       swap_y=0
       if (n2.lt.0.0d0) then
        swap_y=1
        n2=-n2
       endif
       if (n1.lt.0.0d0) then
        swap_x=1
        n1=-n1
       endif
       
       if (n1.gt.0.0d0) then
        d=(volume+0.5d0*(c(2)**2)*n2/n1)*n1/c(2)
        a=d/n1
        b=(d-n2*c(2))/n1
        if ((a.ge.0.0d0).and. &
            (a.le.c(1)).and. &
            (b.ge.0.0d0).and. &
            (b.le.c(1)).and. &
            (a.ge.b)) then
         trap_flag=1
         vol_rect=b*c(2)
         vol_tri=0.5d0*(a-b)*c(2)
         volume_test=vol_rect+vol_tri
         if (abs(volume_test-volume).le.1.0D-10*volcell) then
          centroid(1)=(vol_rect*0.5d0*b+ &
                       vol_tri*(2.0d0*b+a)/3.0d0)/volume
          centroid(2)=(vol_rect*0.5d0*c(2)+ &
                       vol_tri*c(2)/3.0d0)/volume
         else
          print *,"volume_test invalid"
          stop
         endif
        else if ((a.ge.c(1)).or. &
                 (b.ge.c(1)).or. &
                 (b.le.0.0d0)) then
         ! do nothing
        else
         print *,"a,b invalid"
         stop
        endif
       else if (n1.eq.0.0d0) then
        ! do nothing
       else
        print *,"n1 invalid"
        stop
       endif

       if (trap_flag.eq.0) then
        if (n2.gt.0.0d0) then
         d=(volume+0.5d0*(c(1)**2)*n1/n2)*n2/c(1)
         b=d/n2
         a=(d-n1*c(1))/n2
         if ((a.ge.0.0d0).and. &
             (a.le.c(2)).and. &
             (b.ge.0.0d0).and. &
             (b.le.c(2)).and. &
             (a.le.b)) then
          trap_flag=2
          vol_rect=a*c(1)
          vol_tri=0.5d0*(b-a)*c(1)
          volume_test=vol_rect+vol_tri
          if (abs(volume_test-volume).le.1.0D-10*volcell) then
           centroid(2)=(vol_rect*0.5d0*a+ &
                        vol_tri*(2.0d0*a+b)/3.0d0)/volume
           centroid(1)=(vol_rect*0.5d0*c(1)+ &
                        vol_tri*c(1)/3.0d0)/volume
          else
           print *,"volume_test invalid"
           stop
          endif
         else if ((a.ge.c(2)).or. &
                  (b.ge.c(2)).or. &
                  (a.le.0.0d0)) then
          ! do nothing
         else
          print *,"a,b invalid"
          stop
         endif
        else if (n2.eq.0.0d0) then
         ! do nothing
        else
         print *,"n2 invalid"
         stop
        endif

        if (trap_flag.eq.0) then
         if ((n1.gt.0.0d0).and. &
             (n2.gt.0.0d0)) then
          d=sqrt(2.0d0*n1*n2*volume)
          a=d/n1
          b=d/n2
          if ((a.ge.0.0d0).and. &
              (a.le.c(1)).and. &
              (b.ge.0.0d0).and. &
              (b.le.c(2))) then

           tri_flag=1
           vol_tri=0.5d0*a*b
           volume_test=vol_tri
           if (abs(volume_test-volume).le.1.0D-10*volcell) then
            centroid(1)=a/3.0d0
            centroid(2)=b/3.0d0
           else
            print *,"volume_test invalid"
            stop
           endif
          else
           print *,"no orientation found"
           stop
          endif
         else
          print *,"n1 or n2 invalid"
          stop
         endif
        else if ((trap_flag.eq.1).or.(trap_flag.eq.2)) then
         ! do nothing
        else
         print *,"trap_flag invalid"
         stop
        endif
       else if (trap_flag.eq.1) then
        ! do nothing
       else
        print *,"trap_flag invalid"
        stop
       endif

       if ((trap_flag.eq.1).or. &
           (trap_flag.eq.2).or. &
           (tri_flag.eq.1)) then
        if (swap_x.eq.0) then
         ! do nothing
        else if (swap_x.eq.1) then
         centroid(1)=c(1)-centroid(1)
        else
         print *,"swap_x invalid"
         stop
        endif
        if (swap_y.eq.0) then
         ! do nothing
        else if (swap_y.eq.1) then
         centroid(2)=c(2)-centroid(2)
        else
         print *,"swap_y invalid"
         stop
        endif
       else
        print *,"trap_flag or tri_flag invalid"
        stop
       endif
      else
       print *,"volfrac invalid"
       stop
      endif

#if (DEBUG_ANALYTICAL_GRAD==1)
      print *,"AFTER"
      print *,"angles ",angles(1)
      print *,"ref_centroid ",ref_centroid(1), &
              ref_centroid(2)
      print *,"volume ",volume
      print *,"c ",c(1),c(2)
      print *,"centroid ",centroid(1), &
              centroid(2)
#endif

   end subroutine mof2d_compute_analytic_gradient

end module mod_mof2d_analytic_centroid

#undef DEBUG_ANALYTICAL_GRAD
