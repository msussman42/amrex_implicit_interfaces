!This file is part of Notus 0.5.0

!Copyright Bordeaux-INP, Université de Bordeaux, CNRS
!
!Contributors:
!Antoine Lemoine, 16-02-2016, antoine.lemoine@bordeaux-inp.fr

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

!> @defgroup point_3d Points formula in 3D
!! @brief Geometric tools relative to points
!! @ingroup computational_geometry_3d

module mod_cg3_points
use amrex_fort_module, only : amrex_real
   implicit none

contains

   !> Compute the cross product of two points
   !!
   !! @param[in] p1, p2: coordinates of two points of the space
   !! @ingroup point_3d
   pure function cg3_cross_product(p1, p2) result(r)
      real(amrex_real), dimension(3), intent(in) :: p1, p2
      real(amrex_real), dimension(3) :: r

      r(1) = p1(2)*p2(3) - p1(3)*p2(2)
      r(2) = p1(3)*p2(1) - p1(1)*p2(3)
      r(3) = p1(1)*p2(2) - p1(2)*p2(1)
   end function cg3_cross_product

   !> Compute the spherical angles from a direction in Cartesian coordinates
   !!
   !! @param[in]  direction: unit vector
   !! @param[out] angles: spherical angles (θ,φ)
   pure subroutine cg3_direction_to_spherical_angles(direction, angles)
      real(amrex_real), dimension(3), intent(in) :: direction
      real(amrex_real), dimension(2), intent(out) :: angles

      real(amrex_real), parameter :: PI = 2d0*acos(0d0)

      if ((abs(direction(3)) - 1d0) < epsilon(1d0)) then
         angles = [atan2(direction(2), direction(1)), acos(direction(3))]
      else
         if (direction(3) > 0d0) then
                 angles(1)=0.0
                 angles(2)=0.0
         else
                 angles(1)=0.0
                 angles(2)=PI
         end if
      end if
   end subroutine cg3_direction_to_spherical_angles

   !> Compute the direction in Cartesian coordinates from spherical angles.
   !!
   !! @param[in]   angles: spherical angles (θ,φ)
   !! @param[out]  direction: resulting unit vector
   pure subroutine cg3_spherical_angles_to_direction(angles, direction)
      real(amrex_real), dimension(2), intent(in) :: angles
      real(amrex_real), dimension(3), intent(out) :: direction

      direction =                            &
         & [                                 &
         &    sin(angles(2))*cos(angles(1)), &
         &    sin(angles(2))*sin(angles(1)), &
         &    cos(angles(2))                 &
         & ]
   end subroutine cg3_spherical_angles_to_direction

end module mod_cg3_points
