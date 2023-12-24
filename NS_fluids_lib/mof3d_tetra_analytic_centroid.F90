!This file is part of Notus 0.5.0

!Copyright Bordeaux-INP, Université de Bordeaux, CNRS
!
!Contributors:
!Antoine Lemoine, 23-06-2021, antoine.lemoine@bordeaux-inp.fr
!Thomas Milcent, 23-06-2021, thomas.milcent@u-bordeaux.fr

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

module mod_mof3d_tetra_analytic_centroid
use amrex_fort_module, only : amrex_real
   implicit none

   private

   real(amrex_real), parameter :: PI_2 = acos(0d0)
   real(amrex_real), parameter :: PI = 2d0*acos(0d0)
   real(amrex_real), parameter :: TAU = 4d0*acos(0d0)

   integer, parameter :: C_THETA = 1
   integer, parameter :: S_THETA = 2
   integer, parameter :: C_PHI = 3
   integer, parameter :: S_PHI = 4

   public :: mof3d_tetra_compute_analytic_gradient
   public :: mof3d_tetra_compute_analytic_gradient_reference
   public :: mof3d_tetra_compute_residual
   public :: mof3d_tetra_compute_transformation
   public :: mof3d_tetra_transform_angles_to_reference
   public :: mof3d_tetra_transform_angles_to_original
   public :: mof3d_tetra_find_best_transformation

contains

   !> Compute the centroid and the gradient of the objective function in a tetrahedral cell.
   !!
   !! @param[in]  p0, p1, p2, p3: Vertices of the tetrahedron.
   !! @param[in]  angles:         Spherical angles.
   !! @param[in]  ref_centroid1:  Coordinates of the reference centroid of the material 1.
   !! @param[in]  ref_centroid2:  Coordinates of the reference centroid of the material 2.
   !! @param[in]  ref_volume:     Reference volume.
   !! @param[out] objective:      Objective function.
   !! @param[out] gradient:       Gradient of the objective function.
   !! @ingroup moment_of_fluid
   pure subroutine mof3d_tetra_compute_analytic_gradient( &
      &               p0, p1, p2, p3,                     &
      &               angles,                             &
      &               ref_centroid1, ref_centroid2,       &
      &               ref_volume,                         &
      &               objective, gradient, t_centroid     &
      &            )
      use mod_cg3_points, only: cg3_cross_product
      logical, parameter :: mof_use_symmetric_reconstruction=.false.
      real(amrex_real), dimension(3), intent(in) :: p0, p1, p2, p3
      real(amrex_real), dimension(2), intent(in) :: angles
      real(amrex_real), dimension(3), intent(in) :: ref_centroid1, ref_centroid2
      real(amrex_real), intent(in) :: ref_volume
      real(amrex_real), intent(out) :: objective
      real(amrex_real), dimension(3), intent(out) :: t_centroid
      real(amrex_real), dimension(2), intent(out) :: gradient

      real(amrex_real), dimension(3,2) :: partial_derivative
      real(amrex_real), dimension(3)   :: centroid, diff1, diff2, sum_diff, special_centroid
      real(amrex_real), dimension(3,3) :: transformation
      real(amrex_real), dimension(2) :: t_angles
      real(amrex_real) :: cell_volume, volume, determinant, dual_coef

      determinant = dot_product(cg3_cross_product(p1 - p0, p2 - p0), p3 - p0)

      ! Ensure that the transformation preserves the orientation of the tetrahedron.
      transformation(:,1) = p1 - p0
      if (determinant < 0d0) then
         transformation(:,2) = p3 - p0
         transformation(:,3) = p2 - p0
      else
         transformation(:,2) = p2 - p0
         transformation(:,3) = p3 - p0
      end if

      volume = ref_volume/abs(determinant)

      ! Ensure that the angles are in their domains of definition:
      ! → θ ∈ [-π,π]
      ! → φ ∈ [0,π]
      !-----------------------------------------------------------

      t_angles = angles

      if (t_angles(2) < 0d0 .or. t_angles(2) > PI) then
         ! First, ensure that φ belongs to [0,2π[.
         t_angles(2) = modulo_tau(t_angles(2))

         ! Ιf φ belongs to ]π,2π[, transform to ]0,π[.
         ! Otherwise, φ belongs to [0,π].
         if (t_angles(2) > PI) then
            t_angles(2) = TAU - t_angles(2)
            t_angles(1) = PI + t_angles(1)
         end if
      end if

      ! Check whether θ belongs to [-π,π].
      if (t_angles(1) < -PI .or. t_angles(1) > PI) then
         t_angles(1) = modulo_tau(t_angles(1) + PI) - PI
      end if

      call compute_derivatives(t_angles, volume, centroid, partial_derivative)

      ! Transform the centroid back to the original configuration.
      t_centroid = matmul(transformation, centroid) + p0

      diff1 = t_centroid - ref_centroid1

      if (mof_use_symmetric_reconstruction) then
         cell_volume = abs(determinant)/6d0
         special_centroid = (cell_volume*(p0 + p1 + p2 + p3)/4d0 - (cell_volume - ref_volume)*ref_centroid2)/ref_volume
         diff2 = t_centroid - special_centroid
         dual_coef = (ref_volume/(cell_volume - ref_volume))**2
         sum_diff = diff1 + dual_coef*diff2

         objective = dot_product(diff1, diff1) + dual_coef*dot_product(diff2, diff2)
      else
         sum_diff = diff1

         objective = dot_product(diff1, diff1)
      end if

      gradient = [2d0*dot_product(sum_diff, matmul(transformation, partial_derivative(:,1))), &
         &        2d0*dot_product(sum_diff, matmul(transformation, partial_derivative(:,2)))]
   end subroutine mof3d_tetra_compute_analytic_gradient

   !> Compute the centroid and the derivatives in a tetrahedral cell.
   !!
   !! @param[in]  p0, p1, p2, p3: vertices of the tetrahedron.
   !! @param[in]  angles:         Spherical angles.
   !! @param[in]  ref_centroid1:  Coordinates of the reference centroid of the material 1.
   !! @param[in]  ref_centroid2:  Coordinates of the reference centroid of the material 2.
   !! @param[in]  ref_volume:     Reference volume.
   !! @param[out] residual:       Residual function.
   !! @param[out] jacobian:       Partial derivatives of the residual.
   !! @ingroup moment_of_fluid
   pure subroutine mof3d_tetra_compute_residual( &
      &               p0, p1, p2, p3,            &
      &               angles,                    &
      &               ref_centroid1,             &
      &               ref_centroid2,             &
      &               ref_volume,                &
      &               residual,                  &
      &               jacobian                   &
      &            )
      use mod_cg3_points, only: cg3_cross_product
      logical, parameter :: mof_use_symmetric_reconstruction=.false.
      real(amrex_real), dimension(3), intent(in) :: p0, p1, p2, p3
      real(amrex_real), dimension(2), intent(in) :: angles
      real(amrex_real), dimension(3), intent(in) :: ref_centroid1, ref_centroid2
      real(amrex_real), intent(in) :: ref_volume
      real(amrex_real), dimension(:), intent(out) :: residual
      real(amrex_real), dimension(:,:), intent(out) :: jacobian

      real(amrex_real), dimension(3,2) :: partial_derivative
      real(amrex_real), dimension(3)   :: centroid, dual_centroid
      real(amrex_real), dimension(3,3) :: transformation
      real(amrex_real), dimension(2) :: t_angles
      real(amrex_real) :: cell_volume, volume, determinant

      determinant = dot_product(cg3_cross_product(p1 - p0, p2 - p0), p3 - p0)

      ! Ensure that the transformation preserves the orientation of the tetrahedron.
      transformation(:,1) = p1 - p0
      if (determinant < 0d0) then
         transformation(:,2) = p3 - p0
         transformation(:,3) = p2 - p0
      else
         transformation(:,2) = p2 - p0
         transformation(:,3) = p3 - p0
      end if

      volume = ref_volume/abs(determinant)

      ! Ensure that the angles are in their domains of definition:
      ! → θ ∈ [-π,π]
      ! → φ ∈ [0,π]
      !-----------------------------------------------------------

      t_angles = angles

      if (t_angles(2) < 0d0 .or. t_angles(2) > PI) then
         ! First, ensure that φ belongs to [0,2π[.
         t_angles(2) = modulo_tau(t_angles(2))

         ! Ιf φ belongs to ]π,2π[, transform to ]0,π[.
         ! Otherwise, φ belongs to [0,π].
         if (t_angles(2) > PI) then
            t_angles(2) = TAU - t_angles(2)
            t_angles(1) = PI + t_angles(1)
         end if
      end if

      ! Check whether θ belongs to [-π,π].
      if (t_angles(1) < -PI .or. t_angles(1) > PI) then
         t_angles(1) = modulo_tau(t_angles(1) + PI) - PI
      end if

      call compute_derivatives(t_angles, volume, centroid, partial_derivative)

      ! Transform the centroid back to the original configuration.
      centroid = matmul(transformation, centroid) + p0

      ! Compute the residual.
      residual(1:3) = centroid - ref_centroid1

      ! Transform the partial derivatives back to the original configuration.
      jacobian(1:3,1) = matmul(transformation, partial_derivative(:,1))
      jacobian(1:3,2) = matmul(transformation, partial_derivative(:,2))

      ! Symmetric contribution.
      if (mof_use_symmetric_reconstruction) then
         ! Compute the volume of the tetrahedron.
         cell_volume = abs(determinant)/6d0
         ! Compute the dual centroid.
         dual_centroid = (cell_volume*(p0 + p1 + p2 + p3)/4d0 - ref_volume*centroid)/(cell_volume - ref_volume)
         ! Compute the residual.
         residual(4:6) = dual_centroid - ref_centroid2
         ! Compute the jacobian.
         jacobian(4:6,:) = ref_volume*jacobian(1:3,:)/(ref_volume - cell_volume)
      end if
   end subroutine mof3d_tetra_compute_residual

   pure subroutine mof3d_tetra_find_best_transformation(tetra, cell_centroid, ref_centroid, best_tetra)
      use mod_cg3_polyhedron, only: t_polyhedron, cg3_create_tetrahedron
      type(t_polyhedron), intent(in) :: tetra
      type(t_polyhedron), intent(out) :: best_tetra
      real(amrex_real), dimension(3), intent(in) :: cell_centroid
      real(amrex_real), dimension(3), intent(in) :: ref_centroid

      ! p1 and p2 as a function of p0 and p3.
      integer, dimension(2,4,4), parameter :: T = &
         & reshape(     &
         &    [         & ! p0, p3
         &       0, 0,  & !  1,  1
         &       4, 3,  & !  2,  1
         &       2, 4,  & !  3,  1
         &       3, 2,  & !  4,  1
         &       3, 4,  & !  1,  2
         &       0, 0,  & !  2,  2
         &       4, 1,  & !  3,  2
         &       1, 3,  & !  4,  2
         &       4, 2,  & !  1,  3
         &       1, 4,  & !  2,  3
         &       0, 0,  & !  3,  3
         &       2, 1,  & !  4,  3
         &       2, 3,  & !  1,  4
         &       3, 1,  & !  2,  4
         &       1, 2,  & !  3,  4
         &       0, 0   & !  4,  4
         &    ],        &
         &    [2, 4, 4] &
         & )

      real(amrex_real), dimension(3,6) :: edge
      real(amrex_real), dimension(3) :: direction
      real(amrex_real) :: angle, norm, max_val, min_val
      integer :: i, max_edge, min_edge

      direction = cell_centroid - ref_centroid
      norm = norm2(direction)

      if (norm <= tiny(norm)) then
         best_tetra = tetra
         return
      end if

      direction = direction/norm

      edge(:,1) = tetra%point(:,1) - cell_centroid
      edge(:,2) = tetra%point(:,2) - cell_centroid
      edge(:,3) = tetra%point(:,3) - cell_centroid
      edge(:,4) = tetra%point(:,4) - cell_centroid

      min_val = huge(min_val)
      max_val = 0d0
      min_edge = 1
      max_edge = 2
      do i = 1, 4
         edge(:,i) = edge(:,i)/norm2(edge(:,i))
         angle = abs(dot_product(direction, edge(:,i)))
         if (angle >= max_val) then
            max_val = angle
            max_edge = i
         end if
         if (angle <= min_val) then
            min_val = angle
            min_edge = i
         end if
      end do

      if (min_edge == max_edge) min_edge = modulo(max_edge, 4) + 1

      call cg3_create_tetrahedron(                     &
         &    tetra%point(:,max_edge),                 &
         &    tetra%point(:,T(1, max_edge, min_edge)), &
         &    tetra%point(:,T(2, max_edge, min_edge)), &
         &    tetra%point(:,min_edge),                 &
         &    best_tetra                               &
         & )
   end subroutine mof3d_tetra_find_best_transformation

   !> Compute the matrix to transform the reference tetrahedron to the original (transformed) tetrahedron.
   !!
   !! The transformation is given by X(ξ) = Aξ + b. This routine returns the matrix A providing the
   !! vertices of the original (transformed) tetrahedron.
   !!
   !! @param[in]  p0, p1, p2, p3: Vertices of the tetrahedron.
   !! @param[out] transformation: Tranformation matrix.
   !! @ingroup moment_of_fluid
   pure subroutine mof3d_tetra_compute_transformation(p0, p1, p2, p3, transformation)
      use mod_cg3_points, only: cg3_cross_product
      real(amrex_real), dimension(3), intent(in) :: p0, p1, p2, p3
      real(amrex_real), dimension(3,3), intent(out) :: transformation

      real(amrex_real) :: determinant

      determinant = dot_product(cg3_cross_product(p1 - p0, p2 - p0), p3 - p0)

      ! Ensure that the transformation preserves the orientation of the tetrahedron.
      transformation(:,1) = p1 - p0
      if (determinant > 0d0) then
         transformation(:,2) = p2 - p0
         transformation(:,3) = p3 - p0
      else
         transformation(:,2) = p3 - p0
         transformation(:,3) = p2 - p0
      end if
   end subroutine mof3d_tetra_compute_transformation

   !> Transform the spherical angles in the original configuration to the reference configuration.
   !!
   !! @param[in]  transformation: Tranformation matrix.
   !! @param[in]  orig_angles:    Angles in the original configuration.
   !! @param[out] ref_angles:     Angles in the reference configuration.
   !! @ingroup moment_of_fluid
   pure subroutine mof3d_tetra_transform_angles_to_reference(transformation, orig_angles, ref_angles)
      use mod_cg3_points, only: cg3_direction_to_spherical_angles, cg3_spherical_angles_to_direction
      real(amrex_real), dimension(3,3), intent(in) :: transformation
      real(amrex_real), dimension(2), intent(in) :: orig_angles
      real(amrex_real), dimension(2), intent(out) :: ref_angles

      real(amrex_real), dimension(3) :: direction

      call cg3_spherical_angles_to_direction(orig_angles, direction)

      direction = matmul(transpose(transformation), direction)
      direction = direction/norm2(direction)

      call cg3_direction_to_spherical_angles(direction, ref_angles)
   end subroutine mof3d_tetra_transform_angles_to_reference

   !> Transform the spherical angles in the reference configuration to the original configuration.
   !!
   !! @param[in]  transformation: Tranformation matrix.
   !! @param[in]  ref_angles:     Angles in the reference configuration.
   !! @param[out] orig_angles:    Angles in the original configuration.
   !! @ingroup moment_of_fluid
   pure subroutine mof3d_tetra_transform_angles_to_original(transformation, ref_angles, orig_angles)
      use mod_cg3_points, only: cg3_direction_to_spherical_angles, cg3_spherical_angles_to_direction
      real(amrex_real), dimension(3,3), intent(in) :: transformation
      real(amrex_real), dimension(2), intent(in) :: ref_angles
      real(amrex_real), dimension(2), intent(out) :: orig_angles

      real(amrex_real), dimension(3) :: direction
      real(amrex_real), dimension(3,3) :: cotransformation

      call cg3_spherical_angles_to_direction(ref_angles, direction)

      ! Compute the co-factor matrix.
      associate(t => transformation)
         cotransformation(:,1) = [t(2,2)*t(3,3) - t(2,3)*t(3,2), t(1,3)*t(3,2) - t(1,2)*t(3,3), t(1,2)*t(2,3) - t(1,3)*t(2,2)]
         cotransformation(:,2) = [t(2,3)*t(3,1) - t(2,1)*t(3,3), t(1,1)*t(3,3) - t(1,3)*t(3,1), t(1,3)*t(2,1) - t(1,1)*t(2,3)]
         cotransformation(:,3) = [t(2,1)*t(3,2) - t(2,2)*t(3,1), t(1,2)*t(3,1) - t(1,1)*t(3,2), t(1,1)*t(2,2) - t(1,2)*t(2,1)]
      end associate

      direction = matmul(cotransformation, direction)
      direction = direction/norm2(direction)

      call cg3_direction_to_spherical_angles(direction, orig_angles)
   end subroutine mof3d_tetra_transform_angles_to_original

   !> Compute the centroid and the gradient of the objective function in the reference tetrahedron.
   !!
   !! @param[in]  angles:         Spherical angles.
   !! @param[in]  ref_centroid:   Coordinates of the reference centroid.
   !! @param[in]  ref_volume:     Reference volume.
   !! @param[out] objective:      Objective function.
   !! @param[out] gradient:       Gradient of the objective function.
   !! @ingroup moment_of_fluid
   pure subroutine mof3d_tetra_compute_analytic_gradient_reference(angles, ref_centroid, volume, objective, gradient)
      real(amrex_real), dimension(2), intent(in) :: angles
      real(amrex_real), dimension(3), intent(in) :: ref_centroid
      real(amrex_real), intent(in) :: volume
      real(amrex_real), intent(out) :: objective
      real(amrex_real), dimension(2), intent(out) :: gradient

      real(amrex_real), dimension(3,2) :: partial_derivative
      real(amrex_real), dimension(3) :: centroid, diff

      call compute_derivatives(angles, volume, centroid, partial_derivative)

      ! Compute the difference between the centroid and the reference centroid in the local chart.
      diff = centroid - ref_centroid
      objective = dot_product(diff, diff)

      gradient = [2d0*dot_product(diff, partial_derivative(:,1)), &
         &        2d0*dot_product(diff, partial_derivative(:,2))]
   end subroutine mof3d_tetra_compute_analytic_gradient_reference

   pure subroutine compute_derivatives(angles, volume, centroid, derivative)
      real(amrex_real), dimension(2), intent(in) :: angles
      real(amrex_real), intent(in) :: volume
      real(amrex_real), dimension(3), intent(out) :: centroid
      real(amrex_real), dimension(3,2), intent(out) :: derivative

      real(amrex_real), dimension(4) :: trigo
      real(amrex_real) :: theta, phi, dual_volume, atan_1m6v, atan6v
      real(amrex_real) :: tan_theta, cos_theta, sin_theta, cot_theta
      real(amrex_real) :: lim1_t1, lim2_t1, lim3_t1
      real(amrex_real) :: lim1_t2, lim2_t2, lim3_t2
      real(amrex_real) :: lim1_t3, lim2_t3, lim3_t3
      real(amrex_real) :: lim1_t4, lim2_t4, lim3_t4
      real(amrex_real) :: lim1_t1dual, lim2_t1dual, lim3_t1dual
      real(amrex_real) :: lim1_t2dual, lim2_t2dual, lim3_t2dual
      real(amrex_real) :: lim1_t3dual, lim2_t3dual, lim3_t3dual
      real(amrex_real) :: lim1_t4dual, lim2_t4dual, lim3_t4dual

      ! Aliases for θ and φ.
      theta = angles(1)
      phi = angles(2)

      ! Compute trigonometric functions.
      trigo(C_THETA) = cos(theta)
      trigo(S_THETA) = sin(theta)
      trigo(C_PHI) = cos(phi)
      trigo(S_PHI) = sin(phi)

      ! Aliases of trigonometric functions.
      cos_theta = trigo(C_THETA)
      sin_theta = trigo(S_THETA)
      tan_theta = sin_theta/cos_theta
      cot_theta = cos_theta/sin_theta
      dual_volume = 1d0/6d0 - volume

      if (theta <= 0d0) then ! θ < θ6
         atan_1m6v = atan(1d0 - 6d0*volume)
         if (theta < -atan_1m6v - PI_2) then ! θ < θ3
            if (theta < atan_1m6v - PI) then ! θ0 < θ < θ1
               lim3_t2 = PI_2 - atan((1d0 - sqrt(6d0*volume*(1d0 - tan_theta)))*cos_theta)

               if (phi <= lim3_t2) then ! φ < φlim3_t2
                  lim1_t4dual = PI_2                                               &
                     &        - atan((-3d0*dual_volume*(sin_theta + cos_theta)     &
                     &        + sqrt(3d0*dual_volume*(3d0*dual_volume*cos_theta**2 &
                     &        + (2d0 - 6d0*dual_volume)*sin_theta*cos_theta        &
                     &        + 3d0*dual_volume*sin_theta**2)))                    &
                     &        / (1d0 - 6d0*dual_volume))

                  if (phi <= lim1_t4dual) then ! 0 < φ < φlim1_t4*
                     call derivatives_triangle4dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else
                     lim1_t2 = PI_2 - atan((1d0 - (1d0 - tan_theta)**2/(6d0*volume))*cos_theta)

                     if (phi < lim1_t2) then ! φlim1_t4* < φ < φlim1_t2
                        call derivatives_quad_edge1dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! φlim1_t2 < φ < φlim3_t2
                        call derivatives_triangle2(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if
                  end if

               else ! φ > φlim3_t2
                  lim3_t4 = PI_2 - atan((-3d0*volume*sin_theta + cos_theta                                           &
                     &    - sqrt(3d0*volume*(2d0*cos_theta**2 - 2d0*sin_theta*cos_theta + 3d0*volume*sin_theta**2))) &
                     &    / (1d0 - 6d0*volume))

                  if (phi < lim3_t4) then ! φlim3_t2 < φ < φlim3_t4
                     call derivatives_quad_edge3dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! φlim3_t4 < φ < π
                     call derivatives_triangle4(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if
               end if
            else if (theta <= -0.75d0*PI) then ! θ1 < θ < θ2
               lim3_t1dual = PI_2 - atan(tan_theta**2/(6d0*dual_volume)*cos_theta)

               if (phi <= lim3_t1dual) then ! φ < φlim3_t1*
                  lim1_t4dual = PI_2 - atan((-3d0*dual_volume*(sin_theta + cos_theta)                         &
                     &        + sqrt(3d0*dual_volume*(3d0*dual_volume*cos_theta**2                            &
                     &        + (2d0 - 6d0*dual_volume)*sin_theta*cos_theta + 3d0*dual_volume*sin_theta**2))) &
                     &        / (1d0 - 6d0*dual_volume))

                  if (phi <= lim1_t4dual) then ! 0 < φ < φlim1_t4*
                     call derivatives_triangle4dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else
                     lim2_t1dual = PI_2 - atan(sqrt(6d0*dual_volume*tan_theta)*cos_theta)

                     if (phi < lim2_t1dual) then ! φlim1_t4* < φ < φlim2_t1*
                        call derivatives_quad_edge1dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! φlim2_t1* < φ < φlim3_t1*
                        call derivatives_triangle1dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if
                  end if

               else ! φ > φlim3_t1*
                  lim3_t4 = PI_2 - atan((-3d0*volume*sin_theta + cos_theta                                           &
                     &    - sqrt(3d0*volume*(2d0*cos_theta**2 - 2d0*sin_theta*cos_theta + 3d0*volume*sin_theta**2))) &
                     &    / (1d0 - 6d0*volume))

                  if (phi < lim3_t4) then ! φlim3_t1* < φ < φlim3_t4
                     call derivatives_quad_edge3dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! φlim3_t4 < φ < π
                     call derivatives_triangle4(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if
               end if
            else ! θ2 < θ < θ3
               lim1_t1dual = PI_2 - atan(1d0/(6d0*dual_volume*tan_theta)*cos_theta)

               if (phi <= lim1_t1dual) then ! φ < φlim1_t1*
                  lim1_t4dual = PI_2 - atan((-3d0*dual_volume*(sin_theta + cos_theta)                         &
                     &        + sqrt(3d0*dual_volume*(3d0*dual_volume*cos_theta**2                            &
                     &        + (2d0 - 6d0*dual_volume)*sin_theta*cos_theta + 3d0*dual_volume*sin_theta**2))) &
                     &        / (1d0 - 6d0*dual_volume))

                  if (phi <= lim1_t4dual) then ! 0 < φ < φlim1_t4*
                     call derivatives_triangle4dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else
                     lim2_t1dual = PI_2 - atan(sqrt(6d0*dual_volume*tan_theta)*cos_theta)

                     if (phi < lim2_t1dual) then ! φlim1_t4* < φ < φlim2_t1*
                        call derivatives_quad_edge1dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! φlim2_t1* < φ < φlim1_t1*
                        call derivatives_triangle1dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if
                  end if
               else ! φ > φlim1_t1*
                  lim2_t4 = PI_2 - atan((sin_theta - 3d0*volume*cos_theta                                            &
                     &    - sqrt(3d0*volume*(3d0*volume*cos_theta**2 - 2d0*sin_theta*cos_theta + 2d0*sin_theta**2))) &
                     &    / (1d0 - 6d0*volume))

                  if (phi < lim2_t4) then ! φlim1_t1* < φ < φlim2_t4
                     call derivatives_quad_edge2dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! φlim2_t4 < φ < π
                     call derivatives_triangle4(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if
               end if
            end if
         else ! θ > θ3
            if (theta < -PI_2) then ! θ3 < θ < θ4
               lim1_t3 = PI_2 - atan((1d0 - sqrt(6d0*volume*(1d0 - cot_theta)))*sin_theta)

               if (phi <= lim1_t3) then ! φ < φlim1_t3
                  lim1_t4dual = PI_2 - atan((-3d0*dual_volume*(sin_theta + cos_theta) &
                     &        + sqrt(3d0*dual_volume*(3d0*dual_volume*cos_theta**2    &
                     &        + (2d0 - 6d0*dual_volume)*sin_theta*cos_theta           &
                     &        + 3d0*dual_volume*sin_theta**2)))                       &
                     &        / (1d0 - 6d0*dual_volume))

                  if (phi <= lim1_t4dual) then ! 0 < φ < φlim1_t4*
                     call derivatives_triangle4dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else
                     lim2_t3 = PI_2 - atan((1d0 - (1d0 - cot_theta)**2/(6d0*volume))*sin_theta)

                     if (phi < lim2_t3) then ! φlim1_t4* < φ < φlim2_t3
                        call derivatives_quad_edge1dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! φlim2_t3 < φ < φlim1_t3
                        call derivatives_triangle3(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if
                  end if
               else ! φ > φlim1_t3
                  lim2_t4 = PI_2 - atan((sin_theta - 3d0*volume*cos_theta                                            &
                     &    - sqrt(3d0*volume*(3d0*volume*cos_theta**2 - 2d0*sin_theta*cos_theta + 2d0*sin_theta**2))) &
                     &    / (1d0 - 6d0*volume))

                  if (phi < lim2_t4) then ! φlim1_t3 < φ < φlim2_t4
                     call derivatives_quad_edge2dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! φlim2_t4 < φ < π
                     call derivatives_triangle4(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if
               end if
            else if (theta < -atan(1d0-1d0/(6d0*volume)) - PI_2) then ! θ4 < θ < θ5
               lim1_t3 = PI_2 - atan((1d0 - sqrt(6d0*volume*(1d0 - cot_theta)))*sin_theta)

               if (phi <= lim1_t3) then ! φ < φlim1_t3
                  lim3_t4dual = PI_2 - atan((cos_theta - 3d0*dual_volume*sin_theta &
                     &        + sqrt(3d0*dual_volume*(2d0*cos_theta**2             &
                     &        - 2d0*sin_theta*cos_theta                            &
                     &        + 3d0*dual_volume*sin_theta**2)))                    &
                     &        / (1d0 - 6d0*dual_volume))

                  if (phi <= lim3_t4dual) then ! 0 < φ < φlim3_t4*
                     call derivatives_triangle4dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else
                     lim3_t3 = PI_2 - atan((1d0 - 1d0/(6d0*volume*(1d0 - cot_theta)))*sin_theta)

                     if (phi < lim3_t3) then ! φlim3_t4* < φ < φlim3_t3
                        call derivatives_quad_edge3(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! φlim3_t3 < φ < φlim1_t3
                        call derivatives_triangle3(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if
                  end if
               else ! φ > φlim1_t3
                  lim2_t4 = PI_2 - atan((sin_theta - 3d0*volume*cos_theta &
                     &    - sqrt(3d0*volume*(3d0*volume*cos_theta**2      &
                     &    - 2d0*sin_theta*cos_theta + 2d0*sin_theta**2))) &
                     &    / (1d0 - 6d0*volume))

                  if (phi < lim2_t4) then ! φlim1_t3 < φ < φlim2_t4
                     call derivatives_quad_edge2dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! φlim2_t4 < φ < π
                     call derivatives_triangle4(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if
               end if
            else ! θ5 < θ < θ6
               lim2_t2dual = PI_2 - atan((1d0 - 1d0/(6d0*dual_volume*(1d0 - tan_theta)))*cos_theta)

               if (phi <= lim2_t2dual) then ! φ < φlim2_t2*
                  lim3_t4dual = PI_2 - atan((cos_theta - 3d0*dual_volume*sin_theta        &
                     &        + sqrt(3d0*dual_volume*(2d0*cos_theta**2                    &
                     &        - 2d0*sin_theta*cos_theta + 3d0*dual_volume*sin_theta**2))) &
                     &        / (1d0 - 6d0*dual_volume))

                  if (phi <= lim3_t4dual) then ! 0 < φ < φlim3_t4*
                     call derivatives_triangle4dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else
                     lim3_t2dual = PI_2 - atan((1d0 - sqrt(6d0*dual_volume*(1d0 - tan_theta)))*cos_theta)

                     if (phi < lim3_t2dual) then ! φlim3_t4* < φ < φlim3_t2*
                        call derivatives_quad_edge3(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! φlim3_t2* < φ < φlim2_t2*
                        call derivatives_triangle2dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if
                  end if
               else ! φ > φlim2_t2*
                  lim2_t4 = PI_2 - atan((sin_theta - 3d0*volume*cos_theta &
                     &    - sqrt(3d0*volume*(3d0*volume*cos_theta**2      &
                     &    - 2d0*sin_theta*cos_theta + 2d0*sin_theta**2))) &
                     &    / (1d0 - 6d0*volume))

                  if (phi < lim2_t4) then ! φlim2_t2* < φ < φlim2_t4
                     call derivatives_quad_edge2dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! φlim2_t4 < φ < π
                     call derivatives_triangle4(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if
               end if
            end if
         end if
      else ! θ > θ6
         atan6v = atan(6d0*volume)

         if (theta <= -atan6v + PI_2) then ! θ < θ9
            if (theta <= atan6v) then ! θ6 < θ < θ7
               lim1_t2dual = PI_2 - atan((1d0 - (1d0 - tan_theta)**2/(6d0*dual_volume))*cos_theta)

               if (phi <= lim1_t2dual) then ! φ < φlim1_t2*
                  lim3_t4dual = PI_2 - atan((cos_theta - 3d0*dual_volume*sin_theta        &
                     &        + sqrt(3d0*dual_volume*(2d0*cos_theta**2                    &
                     &        - 2d0*sin_theta*cos_theta + 3d0*dual_volume*sin_theta**2))) &
                     &        / (1d0 - 6d0*dual_volume))

                  if (phi <= lim3_t4dual) then ! 0 < φ < φlim3_t4*
                     call derivatives_triangle4dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else
                     lim3_t2dual = PI_2 - atan((1d0 - sqrt(6d0*dual_volume*(1d0 - tan_theta)))*cos_theta)

                     if (phi < lim3_t2dual) then ! φlim3_t4* < φ < φlim3_t2*
                        call derivatives_quad_edge3(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! φlim3_t2* < φ < φlim1_t2*
                        call derivatives_triangle2dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if
                  end if
               else ! φ > φlim1_t2*
                  lim1_t4 = PI_2 - atan((-3d0*volume*(sin_theta + cos_theta)                    &
                     &    - sqrt(3d0*volume*(3d0*volume*cos_theta**2                            &
                     &    + (2d0 - 6d0*volume)*cos_theta*sin_theta + 3d0*volume*sin_theta**2))) &
                     &    / (1d0 - 6d0*volume))

                  if (phi < lim1_t4) then ! φlim1_t2* < φ < φlim1_t4
                     call derivatives_quad_edge1(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! φlim1_t4 < φ < π
                     call derivatives_triangle4(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if
               end if
            else if (theta <= 0.25d0*PI) then ! θ7 < θ < θ8
               lim2_t1 = PI_2 - atan(sqrt(6d0*volume*tan_theta)*cos_theta)

               if (phi <= lim2_t1) then ! φ < φlim2_t1
                  lim3_t4dual = PI_2 - atan((cos_theta - 3d0*dual_volume*sin_theta        &
                     &        + sqrt(3d0*dual_volume*(2d0*cos_theta**2                    &
                     &        - 2d0*sin_theta*cos_theta + 3d0*dual_volume*sin_theta**2))) &
                     &        / (1d0 - 6d0*dual_volume))

                  if (phi <= lim3_t4dual) then ! 0 < φ < φlim3_t4*
                     call derivatives_triangle4dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else
                     lim3_t1 = PI_2 - atan(tan_theta**2/(6d0*volume)*cos_theta)

                     if (phi < lim3_t1) then ! φlim3_t4* < φ < φlim3_t1
                        call derivatives_quad_edge3(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! φlim3_t1 < φ < φlim2_t1
                        call derivatives_triangle1(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if
                  end if
               else ! φ > φlim2_t1
                  lim1_t4 = PI_2 - atan((-3d0*volume*(sin_theta + cos_theta)                    &
                     &    - sqrt(3d0*volume*(3d0*volume*cos_theta**2                            &
                     &    + (2d0 - 6d0*volume)*cos_theta*sin_theta + 3d0*volume*sin_theta**2))) &
                     &    / (1d0 - 6d0*volume))

                  if (phi < lim1_t4) then ! φlim2_t1 < φ < φlim1_t4
                     call derivatives_quad_edge1(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! φlim1_t4 < φ < π
                     call derivatives_triangle4(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if
               end if
            else ! θ8 < θ < θ9
               lim2_t1 = PI_2 - atan(sqrt(6d0*volume*tan_theta)*cos_theta)

               if (phi <= lim2_t1) then ! φ < φlim2_t1
                  lim2_t4dual = PI_2 - atan((sin_theta - 3d0*dual_volume*cos_theta &
                     &        + sqrt(3d0*dual_volume*(3d0*dual_volume*cos_theta**2 &
                     &        - 2d0*sin_theta*cos_theta + 2d0*sin_theta**2)))      &
                     &        / (1d0 - 6d0*dual_volume))

                  if (phi <= lim2_t4dual) then ! 0 < φ < φlim2_t4*
                     call derivatives_triangle4dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else
                     lim1_t1 = PI_2 - atan(1d0/(6d0*volume*tan_theta)*cos_theta)

                     if (phi < lim1_t1) then ! φlim2_t4* < φ < φlim1_t1
                        call derivatives_quad_edge2(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! φlim1_t1 < φ < φlim2_t1
                        call derivatives_triangle1(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if
                  end if
               else ! φ > φlim2_t1
                  lim1_t4 = PI_2 - atan((-3d0*volume*(sin_theta + cos_theta)                    &
                     &    - sqrt(3d0*volume*(3d0*volume*cos_theta**2                            &
                     &    + (2d0 - 6d0*volume)*cos_theta*sin_theta + 3d0*volume*sin_theta**2))) &
                     &    / (1d0 - 6d0*volume))

                  if (phi < lim1_t4) then ! φlim2_t1 < φ < φlim1_t4
                     call derivatives_quad_edge1(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! φlim1_t4 < φ < π
                     call derivatives_triangle4(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if
               end if
            end if
         else ! θ > θ9
            if (theta <= PI_2) then ! θ9 < θ < θ10
               lim2_t3dual = PI_2 - atan((1d0 - (1d0 - cot_theta)**2/(6d0*dual_volume))*sin_theta)

               if (phi <= lim2_t3dual) then ! φ < φlim2_t3*
                  lim2_t4dual = PI_2 - atan((sin_theta - 3d0*dual_volume*cos_theta &
                     &        + sqrt(3d0*dual_volume*(3d0*dual_volume*cos_theta**2 &
                     &        - 2d0*sin_theta*cos_theta + 2d0*sin_theta**2)))      &
                     &        / (1d0 - 6d0*dual_volume))

                  if (phi <= lim2_t4dual) then ! 0 < φ < φlim2_t4*
                     call derivatives_triangle4dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else
                     lim1_t3dual = PI_2 - atan((1d0 - sqrt(6d0*dual_volume*(1d0 - cot_theta)))*sin_theta)

                     if (phi < lim1_t3dual) then ! φlim2_t4* < φ < φlim1_t3*
                        call derivatives_quad_edge2(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! φlim1_t3* < φ < φlim2_t3*
                        call derivatives_triangle3dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if
                  end if
               else ! φ > φlim2_t3*
                  lim1_t4 = PI_2 - atan((-3d0*volume*(sin_theta + cos_theta) &
                     &    - sqrt(3d0*volume*(3d0*volume*cos_theta**2         &
                     &    + (2d0 - 6d0*volume)*cos_theta*sin_theta           &
                     &    + 3d0*volume*sin_theta**2)))                       &
                     &    / (1d0 - 6d0*volume))

                  if (phi < lim1_t4) then ! φlim2_t3* < φ < φlim1_t4
                     call derivatives_quad_edge1(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! φlim1_t4 < φ < π
                     call derivatives_triangle4(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if
               end if
            else if (theta <= atan(1d0-1d0/(6d0*volume)) + PI) then ! θ10 < θ < θ11
               lim3_t3dual = PI_2 - atan((1d0 - 1d0/(6d0*dual_volume*(1d0 - cot_theta)))*sin_theta)

               if (phi <= lim3_t3dual) then ! φ < φlim3_t3*
                  lim2_t4dual = PI_2 - atan((sin_theta - 3d0*dual_volume*cos_theta &
                     &        + sqrt(3d0*dual_volume*(3d0*dual_volume*cos_theta**2 &
                     &        - 2d0*sin_theta*cos_theta + 2d0*sin_theta**2)))      &
                     &        / (1d0 - 6d0*dual_volume))

                  if (phi <= lim2_t4dual) then ! 0 < φ < φlim2_t4*
                     call derivatives_triangle4dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else
                     lim1_t3dual = PI_2 - atan((1d0 - sqrt(6d0*dual_volume*(1d0 - cot_theta)))*sin_theta)

                     if (phi < lim1_t3dual) then ! φlim2_t4* < φ < φlim1_t3*
                        call derivatives_quad_edge2(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! φlim1_t3* < φ < φlim3_t3*
                        call derivatives_triangle3dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if
                  end if
               else ! φ > φlim3_t3*
                  lim3_t4 = PI_2 - atan((-3d0*volume*sin_theta + cos_theta              &
                     &    - sqrt(3d0*volume*(2d0*cos_theta**2 - 2d0*sin_theta*cos_theta &
                     &    + 3d0*volume*sin_theta**2)))                                  &
                     &    / (1d0 - 6d0*volume))

                  if (phi < lim3_t4) then ! φlim3_t3* < φ < φlim3_t4
                     call derivatives_quad_edge3dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! φlim3_t4 < φ < π
                     call derivatives_triangle4(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if
               end if
            else ! θ11 < θ < θ12
               lim3_t2 = PI_2 - atan((1d0 - sqrt(6d0*volume*(1d0 - tan_theta)))*cos_theta)

               if (phi <= lim3_t2) then ! φ < φlim3_t2
                  lim2_t4dual = PI_2 - atan((sin_theta - 3d0*dual_volume*cos_theta &
                     &        + sqrt(3d0*dual_volume*(3d0*dual_volume*cos_theta**2 &
                     &        - 2d0*sin_theta*cos_theta + 2d0*sin_theta**2)))      &
                     &        / (1d0 - 6d0*dual_volume))

                  if (phi <= lim2_t4dual) then ! 0 < φ < φlim2_t4*
                     call derivatives_triangle4dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else
                     lim2_t2 = PI_2 - atan((1d0 - 1d0/(6d0*volume*(1d0 - tan_theta)))*cos_theta)

                     if (phi < lim2_t2) then ! φlim2_t4* < φ < φlim2_t2
                        call derivatives_quad_edge2(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! φlim2_t2 < φ < φlim3_t2
                        call derivatives_triangle2(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if
                  end if
               else ! φ > φlim3_t2
                  lim3_t4 = PI_2 - atan((-3d0*volume*sin_theta + cos_theta       &
                     &    - sqrt(3d0*volume*(2d0*cos_theta**2                    &
                     &    - 2d0*sin_theta*cos_theta + 3d0*volume*sin_theta**2))) &
                     &    / (1d0 - 6d0*volume))

                  if (phi < lim3_t4) then ! φlim3_t2 < φ < φlim3_t4
                     call derivatives_quad_edge3dual(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! φlim3_t4 < φ < π
                     call derivatives_triangle4(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if
               end if
            end if
         end if
      end if
   end subroutine compute_derivatives

   pure subroutine derivatives_triangle1(trigo, volume, derivative_theta, derivative_phi, centroid)
      real(amrex_real), dimension(4), intent(in) :: trigo
      real(amrex_real), intent(in) :: volume
      real(amrex_real), dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      real(amrex_real) :: alpha, beta, gamma, t2l, t3l, t2sq, coef
      real(amrex_real) :: tan_theta, sec_theta, cot_phi, csc_phi, cos_theta, sin_theta

      tan_theta = trigo(S_THETA)/trigo(C_THETA)
      cos_theta = trigo(C_THETA)
      sin_theta = trigo(S_THETA)
      sec_theta = 1d0/trigo(C_THETA)
      cot_phi = trigo(C_PHI)/trigo(S_PHI)
      csc_phi = 1d0/trigo(S_PHI)

      t2l = tan_theta
      t3l = cot_phi*sec_theta

      alpha = (6d0*volume*sin_theta*cot_phi/cos_theta**2)**(1d0/3d0)
      beta = alpha*cos_theta/sin_theta
      gamma = alpha*cos_theta/cot_phi
      ! Centroid
      centroid = [alpha, beta, gamma]/4d0

      t2sq = t2l**2
      coef = volume/2d0/(alpha**2)
      derivative_theta = coef*[t3l*(1d0 + 2d0*t2sq), -t3l/t2l*(2d0 + t2sq), 1d0 - t2sq]
      derivative_phi = coef*csc_phi**2*sec_theta*[-t2l, -1d0, 2d0*t2l/t3l]
   end subroutine derivatives_triangle1

   pure subroutine derivatives_triangle1dual(trigo, volume, derivative_theta, derivative_phi, centroid)
      real(amrex_real), dimension(4), intent(in) :: trigo
      real(amrex_real), intent(in) :: volume
      real(amrex_real), dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      real(amrex_real) :: alpha, beta, gamma, dual_volume, t2l, t3l, t2sq, coef
      real(amrex_real) :: tan_theta, sec_theta, cot_phi, csc_phi, cos_theta, sin_theta

      tan_theta = trigo(S_THETA)/trigo(C_THETA)
      sec_theta = 1d0/trigo(C_THETA)
      cot_phi = trigo(C_PHI)/trigo(S_PHI)
      csc_phi = 1d0/trigo(S_PHI)
      cos_theta = trigo(C_THETA)
      sin_theta = trigo(S_THETA)

      t2l = tan_theta
      t3l = cot_phi*sec_theta
      dual_volume = 1d0/6d0 - volume

      alpha = (6d0*dual_volume*sin_theta*cot_phi/cos_theta**2)**(1d0/3d0)
      beta = alpha*cos_theta/sin_theta
      gamma = alpha*cos_theta/cot_phi

      ! Centroid
      centroid = ([1d0, 1d0, 1d0]/6d0 - dual_volume*[alpha, beta, gamma])/(4d0*volume)

      t2sq = t2l**2
      coef = -dual_volume/volume*dual_volume/2d0/(alpha**2)

      derivative_theta = coef*[t3l*(1d0 + 2d0*t2sq), -t3l/t2l*(2d0 + t2sq), 1d0 - t2sq]
      derivative_phi = coef*csc_phi**2*sec_theta*[-t2l, -1d0, 2d0*t2l/t3l]
   end subroutine derivatives_triangle1dual

   pure subroutine derivatives_triangle2(trigo, volume, derivative_theta, derivative_phi, centroid)
      real(amrex_real), dimension(4), intent(in) :: trigo
      real(amrex_real), intent(in) :: volume
      real(amrex_real), dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      real(amrex_real) :: alpha, beta, gamma, t2l, t3l, coef, t23l
      real(amrex_real) :: tan_theta, sec_theta, cot_phi, csc_phi, cos_theta, sin_theta
      real(amrex_real) :: dthetat2l, dthetat3l
      real(amrex_real) :: dtheta_alpha, dtheta_beta, dtheta_gamma

      tan_theta = trigo(S_THETA)/trigo(C_THETA)
      sec_theta = 1d0/trigo(C_THETA)
      cot_phi = trigo(C_PHI)/trigo(S_PHI)
      csc_phi = 1d0/trigo(S_PHI)
      cos_theta = trigo(C_THETA)
      sin_theta = trigo(S_THETA)

      t2l = tan_theta
      t3l = cot_phi*sec_theta

      alpha = (6d0*volume*(cos_theta - sin_theta)*(cos_theta - cot_phi)/cos_theta**2)**(1d0/3d0)
      beta = alpha*cos_theta/(cos_theta - cot_phi)
      gamma = alpha*cos_theta/(cos_theta - sin_theta)

      ! Centroid
      centroid = [4d0 - alpha - beta - gamma, gamma, beta]/4d0

      dthetat2l = 1d0 + t2l**2
      dthetat3l = t2l*t3l
      t23l = (1d0 - t2l)/(1d0 - t3l)
      coef = -volume/(2d0*alpha**2)
      dtheta_alpha = (dthetat3l - 2d0/t23l*dthetat2l)
      dtheta_beta = (dthetat2l - 2d0*t23l*dthetat3l)
      dtheta_gamma = (1d0 - t3l)*dthetat2l + (1d0 - t2l)*dthetat3l

      derivative_theta = coef*[-dtheta_alpha - dtheta_beta - dtheta_gamma, dtheta_alpha, dtheta_beta]
      derivative_phi = -coef*csc_phi**2*sec_theta*[-1d0 + 2d0*t23l - (1d0 - t2l), 1d0, -2d0*t23l]
   end subroutine derivatives_triangle2

   pure subroutine derivatives_triangle2dual(trigo, volume, derivative_theta, derivative_phi, centroid)
      real(amrex_real), dimension(4), intent(in) :: trigo
      real(amrex_real), intent(in) :: volume
      real(amrex_real), dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      real(amrex_real) :: alpha, beta, gamma, dual_volume, t2l, t3l, t23l, coef, dthetat2l, dthetat3l
      real(amrex_real) :: tan_theta, sec_theta, cot_phi, csc_phi, cos_theta, sin_theta
      real(amrex_real) :: dtheta_alpha, dtheta_beta, dtheta_gamma

      tan_theta = trigo(S_THETA)/trigo(C_THETA)
      sec_theta = 1d0/trigo(C_THETA)
      cot_phi = trigo(C_PHI)/trigo(S_PHI)
      csc_phi = 1d0/trigo(S_PHI)
      cos_theta = trigo(C_THETA)
      sin_theta = trigo(S_THETA)

      t2l = tan_theta
      t3l = cot_phi*sec_theta
      dual_volume = 1d0/6d0 - volume

      alpha = (6d0*dual_volume*(cos_theta - sin_theta)*(cos_theta - cot_phi)/cos_theta**2)**(1d0/3d0)
      beta = alpha*cos_theta/(cos_theta - cot_phi)
      gamma = alpha*cos_theta/(cos_theta - sin_theta)

      ! Centroid
      centroid = ([1d0, 1d0, 1d0]/6d0 - dual_volume*[4d0 - alpha - beta - gamma, gamma, beta])/(4d0*volume)

      dthetat2l = 1d0 + t2l**2
      dthetat3l = t2l*t3l
      t23l = (1d0 - t2l)/(1d0 - t3l)
      coef = dual_volume/volume*dual_volume/(2d0*alpha**2)
      dtheta_alpha = (dthetat3l - 2d0/t23l*dthetat2l)
      dtheta_beta = (dthetat2l - 2d0*t23l*dthetat3l)
      dtheta_gamma = ((1d0 - t3l)*dthetat2l + (1d0 - t2l)*dthetat3l)

      derivative_theta = coef*[-dtheta_alpha - dtheta_beta - dtheta_gamma, dtheta_alpha, dtheta_beta]
      derivative_phi = -coef*csc_phi**2*sec_theta*[-1d0 + 2d0*t23l - (1d0 - t2l), 1d0, - 2d0*t23l]
   end subroutine derivatives_triangle2dual

   pure subroutine derivatives_triangle3(trigo, volume, derivative_theta, derivative_phi, centroid)
      real(amrex_real), dimension(4), intent(in) :: trigo
      real(amrex_real), intent(in) :: volume
      real(amrex_real), dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      real(amrex_real) :: alpha, beta, gamma
      real(amrex_real) :: tan_theta, sec_theta, cot_phi, cos_theta, sin_theta, csc_theta, cot_theta, csc_phi
      real(amrex_real) :: dtheta_alpha, dtheta_beta, dtheta_gamma, dphi_alpha, dphi_beta, dphi_gamma

      tan_theta = trigo(S_THETA)/trigo(C_THETA)
      cot_theta = trigo(C_THETA)/trigo(S_THETA)
      sec_theta = 1d0/trigo(C_THETA)
      csc_theta = 1d0/trigo(S_THETA)
      cos_theta = trigo(C_THETA)
      sin_theta = trigo(S_THETA)
      cot_phi = trigo(C_PHI)/trigo(S_PHI)
      csc_phi = 1d0/trigo(S_PHI)

      alpha = (6d0*volume*(sin_theta - cos_theta)*sin_theta/(sin_theta - cot_phi)**2)**(1d0/3d0)
      beta = alpha*(sin_theta - cot_phi)/sin_theta
      gamma = alpha*(sin_theta - cot_phi)/(sin_theta - cos_theta)

      ! Centroid
      centroid = [gamma, 4d0 - alpha - beta - gamma, alpha]/4d0

      dtheta_alpha = (-1d0 + cot_phi*(2d0*cos_theta + 2d0*sin_theta - csc_theta)) &
         &         / (3d0*(cot_phi - sin_theta)*(-cos_theta + sin_theta))
      dtheta_beta = (csc_theta + cot_phi*(1d0 + cot_theta - 2d0*csc_theta**2))/(3d0*(-cos_theta + sin_theta))
      dtheta_gamma = (-2d0 + cot_phi*(cos_theta + csc_theta + sin_theta))/(3d0*(-cos_theta + sin_theta)**2)

      derivative_theta = alpha/4d0*[dtheta_gamma, - dtheta_alpha - dtheta_beta - dtheta_gamma, dtheta_alpha]

      dphi_alpha = 2d0/(3d0*(cot_phi - sin_theta))
      dphi_beta = csc_theta/3d0
      dphi_gamma = 1d0/(3d0*(-cos_theta + sin_theta))

      derivative_phi = alpha/4d0*csc_phi**2*[dphi_gamma, - dphi_alpha - dphi_beta - dphi_gamma, dphi_alpha]
   end subroutine derivatives_triangle3

   pure subroutine derivatives_triangle3dual(trigo, volume, derivative_theta, derivative_phi, centroid)
      real(amrex_real), dimension(4), intent(in) :: trigo
      real(amrex_real), intent(in) :: volume
      real(amrex_real), dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      real(amrex_real) :: alpha, beta, gamma, dual_volume, coef
      real(amrex_real) :: tan_theta, sec_theta, cot_phi, cos_theta, sin_theta, csc_theta, cot_theta, csc_phi
      real(amrex_real) :: dtheta_alpha, dtheta_beta, dtheta_gamma, dphi_alpha, dphi_beta, dphi_gamma

      tan_theta = trigo(S_THETA)/trigo(C_THETA)
      cot_theta = trigo(C_THETA)/trigo(S_THETA)
      sec_theta = 1d0/trigo(C_THETA)
      csc_theta = 1d0/trigo(S_THETA)
      cos_theta = trigo(C_THETA)
      sin_theta = trigo(S_THETA)
      cot_phi = trigo(C_PHI)/trigo(S_PHI)
      csc_phi = 1d0/trigo(S_PHI)

      dual_volume = 1d0/6d0 - volume

      alpha = (6d0*dual_volume*(sin_theta - cos_theta)*sin_theta/(sin_theta - cot_phi)**2)**(1d0/3d0)
      beta = alpha*(sin_theta - cot_phi)/sin_theta
      gamma = alpha*(sin_theta - cot_phi)/(sin_theta - cos_theta)

      ! Centroid
      centroid = ([1d0, 1d0, 1d0]/6d0 - dual_volume*[gamma, 4d0 - alpha - beta - gamma, alpha])/(4d0*volume)

      coef = -dual_volume/volume/4d0*alpha

      dtheta_alpha = (-1d0 + cot_phi*(2d0*cos_theta + 2d0*sin_theta - csc_theta)) &
         &         / (3d0*(cot_phi - sin_theta)*(-cos_theta + sin_theta))
      dtheta_beta = (csc_theta + cot_phi*(1d0 + cot_theta - 2d0*csc_theta**2))/(3d0*(-cos_theta + sin_theta))
      dtheta_gamma = (-2d0 + cot_phi*(cos_theta + csc_theta + sin_theta))/(3d0*(-cos_theta + sin_theta)**2)

      derivative_theta = coef*[dtheta_gamma, - dtheta_alpha - dtheta_beta - dtheta_gamma, dtheta_alpha]

      dphi_alpha = 2d0/(3d0*(cot_phi - sin_theta))
      dphi_beta = csc_theta/3d0
      dphi_gamma = 1d0/(3d0*(-cos_theta + sin_theta))

      derivative_phi = coef*csc_phi**2*[dphi_gamma, - dphi_alpha - dphi_beta - dphi_gamma, dphi_alpha]
   end subroutine derivatives_triangle3dual

   pure subroutine derivatives_triangle4(trigo, volume, derivative_theta, derivative_phi, centroid)
      real(amrex_real), dimension(4), intent(in) :: trigo
      real(amrex_real), intent(in) :: volume
      real(amrex_real), dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      real(amrex_real) :: alpha, beta, gamma, tc2, cst, isct, isst, ct, st, coef
      real(amrex_real) :: cos_theta, sin_theta, tan_phi, sec_phi
      real(amrex_real) :: dtheta_alpha, dtheta_beta, dtheta_gamma, dphi_alpha, dphi_beta, dphi_gamma

      cos_theta = trigo(C_THETA)
      sin_theta = trigo(S_THETA)
      tan_phi = trigo(S_PHI)/trigo(C_PHI)
      sec_phi = 1d0/trigo(C_PHI)

      ct = cos_theta*tan_phi
      st = sin_theta*tan_phi

      alpha = (6d0*volume*(1d0 - ct)/(1d0 - st)**2)**(1d0/3d0)
      beta = alpha*(1 - st)/(1 - ct)
      gamma = alpha*(1 - st)
      isst = 1d0/(1d0 - st)

      ! Centroid
      centroid = [beta, alpha, 4d0 - alpha - beta - gamma]/4d0

      isct = 1d0/(ct - 1d0)
      tc2 = ct*cos_theta
      coef = volume/(alpha*(1d0 - st))**2/2d0

      dtheta_alpha = sin_theta - cos_theta + 2d0*tc2 - tan_phi
      dtheta_beta =  (cos_theta + 2d0*sin_theta -2d0*tan_phi + tc2)*isct
      dtheta_gamma = (sin_theta + 2d0*cos_theta -    tan_phi - tc2)*isst

      derivative_theta = coef*tan_phi*[dtheta_beta, dtheta_gamma, - dtheta_alpha - dtheta_beta - dtheta_gamma]

      cst = ct*sin_theta
      dphi_alpha = 2d0*cst - cos_theta - sin_theta
      dphi_beta = (cst + sin_theta - 2d0*cos_theta)*isct
      dphi_gamma = -(cst + cos_theta - 2d0*sin_theta)*isst

      derivative_phi = coef*sec_phi**2*[dphi_beta, dphi_gamma, - dphi_alpha - dphi_beta - dphi_gamma]
   end subroutine derivatives_triangle4

   pure subroutine derivatives_triangle4dual(trigo, volume, derivative_theta, derivative_phi, centroid)
      real(amrex_real), dimension(4), intent(in) :: trigo
      real(amrex_real), intent(in) :: volume
      real(amrex_real), dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      real(amrex_real) :: alpha, beta, gamma, dual_volume, tc2, cst, isct, isst, ct, st, coef
      real(amrex_real) :: cos_theta, sin_theta, tan_phi, sec_phi
      real(amrex_real) :: dtheta_alpha, dtheta_beta, dtheta_gamma, dphi_alpha, dphi_beta, dphi_gamma

      cos_theta = trigo(C_THETA)
      sin_theta = trigo(S_THETA)
      tan_phi = trigo(S_PHI)/trigo(C_PHI)
      sec_phi = 1d0/trigo(C_PHI)

      dual_volume = 1d0/6d0 - volume
      ct = cos_theta*tan_phi
      st = sin_theta*tan_phi

      alpha = (6d0*dual_volume*(1d0 - ct)/(1d0 - st)**2)**(1d0/3d0)
      beta = alpha*(1 - st)/(1 - ct)
      gamma = alpha*(1 - st)
      isst = 1d0/(1d0 - st)

      ! Centroid
      centroid = ([1d0, 1d0, 1d0]/6d0 - dual_volume*[beta, alpha, 4d0 - alpha - beta - gamma])/(4d0*volume)

      isct = 1d0/(ct - 1d0)
      tc2 = ct*cos_theta
      coef = -dual_volume/(2d0*volume)*dual_volume/(alpha*(1d0 - st))**2

      dtheta_alpha = sin_theta - cos_theta + 2d0*tc2 - tan_phi
      dtheta_beta =  (cos_theta + 2d0*sin_theta - 2d0*tan_phi + tc2)*isct
      dtheta_gamma = (sin_theta + 2d0*cos_theta -     tan_phi - tc2)*isst

      derivative_theta = coef*tan_phi*[dtheta_beta, dtheta_gamma, - dtheta_alpha - dtheta_beta - dtheta_gamma]

      cst = ct*sin_theta
      dphi_alpha = 2d0*cst - cos_theta - sin_theta
      dphi_beta = (cst + sin_theta - 2d0*cos_theta)*isct
      dphi_gamma = -(cst + cos_theta - 2d0*sin_theta)*isst

      derivative_phi = coef*sec_phi**2*[dphi_beta, dphi_gamma, - dphi_alpha - dphi_beta - dphi_gamma]
   end subroutine derivatives_triangle4dual

   pure subroutine derivatives_quad_edge1(trigo, volume, derivative_theta, derivative_phi, centroid)
      real(amrex_real), dimension(4), intent(in) :: trigo
      real(amrex_real), intent(in) :: volume
      real(amrex_real), dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      real(amrex_real) :: alpha, alphacos, cossin, beta, gamma, delta
      real(amrex_real) :: sin_theta, cos_theta, cot_phi
      real(amrex_real) :: sqrt_t23, xqe, t232t

      sin_theta = trigo(S_THETA)
      cos_theta = trigo(C_THETA)
      cot_phi = trigo(C_PHI)/trigo(S_PHI)

      t232t = (cos_theta + sin_theta - cot_phi)**2
      cossin = cos_theta*sin_theta
      sqrt_t23 = sqrt(cossin*(sin_theta - cot_phi)*(cos_theta- cot_phi))
      xqe = cos((acos(1d0/(2d0*sqrt_t23)*(2d0*cossin - cot_phi*(cos_theta + sin_theta - cot_phi) - 6d0*volume*t232t)) + 4d0*PI)/3d0)

      alpha = (2d0*sqrt_t23*xqe/cos_theta + sin_theta)/(cos_theta + sin_theta - cot_phi)
      alphacos = alpha*cos_theta
      beta = alphacos/sin_theta
      gamma = (cot_phi - alphacos)/(cot_phi - cos_theta)
      delta = (sin_theta - alphacos)/(sin_theta - cot_phi)

      centroid(1) = f0qex(sin_theta, cot_phi, cos_theta, volume)                   &
         &        + xqe*(f1qexbis(sin_theta, cot_phi, cos_theta, volume, sqrt_t23) &
         &        + xqe*f2qex(sin_theta, cot_phi, cos_theta))
      centroid(2) = f0qex(cos_theta, cot_phi, sin_theta, volume)                   &
         &        + xqe*(f1qexbis(cos_theta, cot_phi, sin_theta, volume, sqrt_t23) &
         &        + xqe*f2qex(cos_theta, cot_phi, sin_theta))
      centroid(3) = f0qez(sin_theta, cot_phi, cos_theta, volume)                   &
         &        + xqe*(f1qezbis(sin_theta, cot_phi, cos_theta, volume, sqrt_t23) &
         &        + xqe*f2qez(sin_theta, cot_phi, cos_theta))

      centroid = centroid/(24d0*volume*t232t**2)

      call quadedge_analytic_derivatives( &
         &    trigo,                      &
         &    [alpha, 0d0, 0d0],          &
         &    [0d0, beta, 0d0],           &
         &    [gamma, 0d0, 1d0 -  gamma], &
         &    [0d0, 1d0 - delta, delta],  &
         &    volume,                     &
         &    derivative_theta,           &
         &    derivative_phi              &
         & )
   end subroutine derivatives_quad_edge1

   pure subroutine derivatives_quad_edge1dual(trigo, volume, derivative_theta, derivative_phi, centroid)
      real(amrex_real), dimension(4), intent(in) :: trigo
      real(amrex_real), intent(in) :: volume
      real(amrex_real), dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      real(amrex_real) :: alpha, alphacos, cossin, beta, gamma, delta, dual_volume
      real(amrex_real) :: cos_theta, sin_theta, cot_phi
      real(amrex_real) :: sqrt_t23, xqe, t232t

      sin_theta = trigo(S_THETA)
      cos_theta = trigo(C_THETA)
      cot_phi = trigo(C_PHI)/trigo(S_PHI)

      dual_volume = 1d0/6d0 - volume
      t232t = (cos_theta + sin_theta - cot_phi)**2
      cossin = cos_theta*sin_theta
      sqrt_t23 = sqrt(cossin*(sin_theta - cot_phi)*(cos_theta - cot_phi))
      xqe = cos((acos((2d0*cossin - cot_phi*(cos_theta + sin_theta - cot_phi) - 6d0*dual_volume*t232t)/2d0/sqrt_t23) + 4d0*PI)/3d0)

      alpha = (2d0*sqrt_t23*xqe/cos_theta + sin_theta)/(cos_theta + sin_theta - cot_phi)
      alphacos = alpha*cos_theta
      beta = alphacos/sin_theta
      gamma = (cot_phi - alphacos)/(cot_phi - cos_theta)
      delta = (sin_theta - alphacos)/(sin_theta - cot_phi)

      centroid(1) = f0qex(sin_theta, cot_phi, cos_theta, dual_volume)                   &
         &        + xqe*(f1qexbis(sin_theta, cot_phi, cos_theta, dual_volume, sqrt_t23) &
         &        + xqe*f2qex(sin_theta, cot_phi, cos_theta))
      centroid(2) = f0qex(cos_theta, cot_phi, sin_theta, dual_volume)                   &
         &        + xqe*(f1qexbis(cos_theta, cot_phi, sin_theta, dual_volume, sqrt_t23) &
         &        + xqe*f2qex(cos_theta, cot_phi, sin_theta))
      centroid(3) = f0qez(sin_theta, cot_phi, cos_theta, dual_volume) &
         &        + xqe*(f1qezbis(sin_theta, cot_phi, cos_theta, dual_volume, sqrt_t23) &
         &        + xqe*f2qez(sin_theta, cot_phi, cos_theta))

      centroid = centroid/(24d0*dual_volume*t232t**2)
      centroid = ([1d0,1d0,1d0]/24d0 - dual_volume*centroid)/volume

      call quadedge_analytic_derivatives( &
         &    trigo,                      &
         &    [alpha, 0d0, 0d0],          &
         &    [0d0, beta, 0d0],           &
         &    [gamma, 0d0, 1d0 - gamma],  &
         &    [0d0, 1d0 - delta, delta],  &
         &    volume,                     &
         &    derivative_theta,           &
         &    derivative_phi              &
         & )
   end subroutine derivatives_quad_edge1dual

   pure subroutine derivatives_quad_edge2(trigo, volume, derivative_theta, derivative_phi, centroid)
      real(amrex_real), dimension(4), intent(in) :: trigo
      real(amrex_real), intent(in) :: volume
      real(amrex_real), dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      real(amrex_real) :: alpha, alphasin, sincot, beta, gamma, delta
      real(amrex_real) :: cot_phi, cos_theta, sin_theta
      real(amrex_real) :: sqrt_t23, xqe, t232t

      cos_theta = trigo(C_THETA)
      sin_theta = trigo(S_THETA)
      cot_phi = trigo(C_PHI)/trigo(S_PHI)

      t232t = (cos_theta - sin_theta - cot_phi)**2
      sincot = sin_theta*cot_phi
      sqrt_t23 = sqrt(sincot*(cos_theta-sin_theta)*(cos_theta-cot_phi))
      xqe = cos((acos(1d0/(2d0*sqrt_t23)*(2d0*sincot + cos_theta*(cos_theta - sin_theta - cot_phi) &
         &- 6d0*volume*t232t)) + 4d0*PI)/3d0)

      alpha = (2d0*sqrt_t23*xqe/sin_theta + cot_phi)/(cot_phi - cos_theta + sin_theta)
      alphasin = sin_theta*alpha
      beta = alphasin/cot_phi
      gamma = (cos_theta - alphasin)/(cos_theta - sin_theta)
      delta = (cot_phi - alphasin)/(cot_phi - cos_theta)

      centroid(1) = f0qez(sin_theta, cos_theta, cot_phi, volume)                   &
         &        + xqe*(f1qezbis(sin_theta, cos_theta, cot_phi, volume, sqrt_t23) &
         &        + xqe*f2qez(sin_theta, cos_theta, cot_phi))
      centroid(2) = f0qex(cot_phi, cos_theta, sin_theta, volume)                   &
         &        + xqe*(f1qexbis(cot_phi, cos_theta, sin_theta, volume, sqrt_t23) &
         &        + xqe*f2qex(cot_phi, cos_theta, sin_theta))
      centroid(3) = f0qex(sin_theta, cos_theta, cot_phi, volume)                   &
         &        + xqe*(f1qexbis(sin_theta, cos_theta, cot_phi, volume, sqrt_t23) &
         &        + xqe*f2qex(sin_theta, cos_theta, cot_phi))

      centroid = centroid/(24d0*volume*t232t**2)

      call quadedge_analytic_derivatives( &
         &    trigo,                      &
         &    [0d0, alpha, 0d0],          &
         &    [0d0, 0d0, beta],           &
         &    [1d0-gamma, gamma, 0d0],    &
         &    [delta, 0d0, 1d0 - delta],  &
         &    volume,                     &
         &    derivative_theta,           &
         &    derivative_phi              &
         & )
   end subroutine derivatives_quad_edge2

   pure subroutine derivatives_quad_edge2dual(trigo, volume, derivative_theta, derivative_phi, centroid)
      real(amrex_real), dimension(4), intent(in) :: trigo
      real(amrex_real), intent(in) :: volume
      real(amrex_real), dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      real(amrex_real) :: alpha, alphasin, sincot, beta, gamma, delta, dual_volume
      real(amrex_real) :: cos_theta, sin_theta, cot_phi
      real(amrex_real) :: sqrt_t23, xqe, t232t

      cos_theta = trigo(C_THETA)
      sin_theta = trigo(S_THETA)
      cot_phi = trigo(C_PHI)/trigo(S_PHI)

      dual_volume = 1d0/6d0 - volume
      t232t = (cos_theta - sin_theta - cot_phi)**2
      sincot = sin_theta*cot_phi
      sqrt_t23 = sqrt(sincot*(cos_theta-sin_theta)*(cos_theta-cot_phi))
      xqe = cos((acos(1d0/(2d0*sqrt_t23)*(2d0*sincot + cos_theta*(cos_theta - sin_theta - cot_phi) &
         &- 6d0*dual_volume*t232t)) + 4d0*PI)/3d0)

      alpha = (2d0*sqrt_t23*xqe/sin_theta + cot_phi)/(cot_phi - cos_theta + sin_theta)
      alphasin = sin_theta*alpha
      beta = alphasin/cot_phi
      gamma = (cos_theta - alphasin)/(cos_theta - sin_theta)
      delta = (cot_phi - alphasin)/(cot_phi - cos_theta)

      centroid(1) = f0qez(sin_theta, cos_theta, cot_phi, dual_volume)                   &
         &        + xqe*(f1qezbis(sin_theta, cos_theta, cot_phi, dual_volume, sqrt_t23) &
         &        + xqe*f2qez(sin_theta, cos_theta, cot_phi))
      centroid(2) = f0qex(cot_phi, cos_theta, sin_theta, dual_volume)                   &
         &        + xqe*(f1qexbis(cot_phi, cos_theta, sin_theta, dual_volume, sqrt_t23) &
         &        + xqe*f2qex(cot_phi, cos_theta, sin_theta))
      centroid(3) = f0qex(sin_theta, cos_theta, cot_phi, dual_volume)                   &
         &        + xqe*(f1qexbis(sin_theta, cos_theta, cot_phi, dual_volume, sqrt_t23) &
         &        + xqe*f2qex(sin_theta, cos_theta, cot_phi))

      centroid = centroid/(24d0*dual_volume*t232t**2)

      centroid = ([1d0,1d0,1d0]/24d0 - dual_volume*centroid)/volume

      call quadedge_analytic_derivatives( &
         &    trigo,                      &
         &    [0d0, alpha, 0d0],          &
         &    [0d0, 0d0, beta],           &
         &    [1d0 - gamma, gamma, 0d0],  &
         &    [delta, 0d0, 1d0 - delta],  &
         &    volume,                     &
         &    derivative_theta,           &
         &    derivative_phi              &
         & )
   end subroutine derivatives_quad_edge2dual

   pure subroutine derivatives_quad_edge3(trigo, volume, derivative_theta, derivative_phi, centroid)
      real(amrex_real), dimension(4), intent(in) :: trigo
      real(amrex_real), intent(in) :: volume
      real(amrex_real), dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      real(amrex_real) :: alpha, alphacot, coscot, beta, gamma, delta
      real(amrex_real) :: cos_theta, sin_theta, cot_phi
      real(amrex_real) :: sqrt_t23, xqe, t232t

      cos_theta = trigo(C_THETA)
      sin_theta = trigo(S_THETA)
      cot_phi = trigo(C_PHI)/trigo(S_PHI)

      t232t = (cos_theta - sin_theta + cot_phi)**2
      coscot = cos_theta*cot_phi
      sqrt_t23 = sqrt(coscot*(cos_theta-sin_theta)*(cot_phi-sin_theta))
      xqe = cos((acos(1d0/(2d0*sqrt_t23)*(2d0*coscot - sin_theta*(cos_theta - sin_theta + cot_phi) &
         &- 6d0*volume*t232t)) + 4d0*PI)/3d0)

      alpha = (2d0*sqrt_t23*xqe/cot_phi + cos_theta)/(cos_theta - sin_theta + cot_phi)
      alphacot = cot_phi*alpha
      beta = alphacot/cos_theta
      gamma = (sin_theta - alphacot)/(sin_theta - cot_phi)
      delta = (cos_theta - alphacot)/(cos_theta - sin_theta)

      centroid(1) = f0qex(cot_phi, sin_theta, cos_theta, volume)                   &
         &        + xqe*(f1qexbis(cot_phi, sin_theta, cos_theta, volume, sqrt_t23) &
         &        + xqe*f2qex(cot_phi, sin_theta, cos_theta))
      centroid(2) = f0qez(cot_phi, sin_theta, cos_theta, volume)                   &
         &        + xqe*(f1qezbis(cot_phi, sin_theta, cos_theta, volume, sqrt_t23) &
         &        + xqe*f2qez(cot_phi, sin_theta, cos_theta))
      centroid(3) = f0qex(cos_theta, sin_theta, cot_phi, volume)                   &
         &        + xqe*(f1qexbis(cos_theta, sin_theta, cot_phi, volume, sqrt_t23) &
         &        + xqe*f2qex(cos_theta, sin_theta, cot_phi))

      centroid = centroid/(24d0*volume*t232t**2)

      call quadedge_analytic_derivatives( &
         &    trigo,                      &
         &    [0d0, 0d0, alpha],          &
         &    [beta, 0d0, 0d0],           &
         &    [0d0, 1d0 - gamma, gamma],  &
         &    [1d0 - delta, delta, 0d0],  &
         &    volume,                     &
         &    derivative_theta,           &
         &    derivative_phi              &
         & )
   end subroutine derivatives_quad_edge3

   pure subroutine derivatives_quad_edge3dual(trigo, volume, derivative_theta, derivative_phi, centroid)
      real(amrex_real), dimension(4), intent(in) :: trigo
      real(amrex_real), intent(in) :: volume
      real(amrex_real), dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      real(amrex_real) :: alpha, coscot, alphacot, beta, gamma, delta, dual_volume
      real(amrex_real) :: cos_theta, sin_theta, cot_phi
      real(amrex_real) :: xqe, sqrt_t23, t232t

      cos_theta = trigo(C_THETA)
      sin_theta = trigo(S_THETA)
      cot_phi = trigo(C_PHI)/trigo(S_PHI)

      dual_volume = 1d0/6d0 - volume
      t232t = (cos_theta - sin_theta + cot_phi)**2
      coscot = cos_theta*cot_phi
      sqrt_t23 = sqrt(coscot*(cos_theta - sin_theta)*(cot_phi - sin_theta))
      xqe = cos((acos(1d0/(2d0*sqrt_t23)*(2d0*coscot - sin_theta*(cos_theta - sin_theta + cot_phi) &
         &- 6d0*dual_volume*t232t)) + 4d0*PI)/3d0)

      alpha = (2d0*sqrt_t23*xqe/cot_phi + cos_theta)/(cos_theta - sin_theta + cot_phi)
      alphacot = cot_phi*alpha
      beta = alphacot/cos_theta
      gamma = (sin_theta - alphacot)/(sin_theta - cot_phi)
      delta = (cos_theta - alphacot)/(cos_theta - sin_theta)

      centroid(1) = f0qex(cot_phi, sin_theta, cos_theta, dual_volume)                   &
         &        + xqe*(f1qexbis(cot_phi, sin_theta, cos_theta, dual_volume, sqrt_t23) &
         &        + xqe*f2qex(cot_phi, sin_theta, cos_theta))
      centroid(2) = f0qez(cot_phi, sin_theta, cos_theta, dual_volume)                   &
         &        + xqe*(f1qezbis(cot_phi, sin_theta, cos_theta, dual_volume, sqrt_t23) &
         &        + xqe*f2qez(cot_phi, sin_theta, cos_theta))
      centroid(3) = f0qex(cos_theta, sin_theta, cot_phi, dual_volume)                   &
         &        + xqe*(f1qexbis(cos_theta, sin_theta, cot_phi, dual_volume, sqrt_t23) &
         &        + xqe*f2qex(cos_theta, sin_theta, cot_phi))

      centroid = centroid/(24d0*dual_volume*t232t**2)
      centroid = ([1d0,1d0,1d0]/24d0 - dual_volume*centroid)/volume

      call quadedge_analytic_derivatives( &
         &    trigo,                      &
         &    [0d0, 0d0, alpha],          &
         &    [beta, 0d0, 0d0],           &
         &    [0d0, 1d0 - gamma, gamma],  &
         &    [1d0 - delta, delta, 0d0],  &
         &    volume,                     &
         &    derivative_theta,           &
         &    derivative_phi              &
         & )
   end subroutine derivatives_quad_edge3dual

   real(amrex_real) pure function g0(x, y, z)
      real(amrex_real), intent(in) :: x, y, z
      g0 = (2d0*x + y + z)**2 - z*(4d0*y + 9d0*x)
   end function g0

   real(amrex_real) pure function g1(x, y, z)
      real(amrex_real), intent(in) :: x, y, z
      g1 = 12d0*x*(x + 2d0*y - z)
   end function g1

   real(amrex_real) pure function f0qex(x, y, z, v)
      real(amrex_real), intent(in) :: x, y, z, v
      real(amrex_real) :: x1, y1, z0
      x1 = x*(x-y)
      y1 = z*(z-y)
      z0 = (z+x-y)**2
      f0qex = g0(x1,y1,z0) + 24d0*v*x1*z0
   end function f0qex

   real(amrex_real) pure function f0qez(x, y, z, v)
      real(amrex_real), intent(in) :: x, y, z, v
      real(amrex_real) :: x2, y2, z0
      x2 = x*z
      y2 = (x-y)*(z-y)
      z0 = (z+x-y)**2
      f0qez = -g0(x2,y2,z0) + z0*(z0 - 4d0*x2) + 24d0*v*x2*z0
   end function f0qez

   real(amrex_real) pure function f1qex(x, y, z, v)
      real(amrex_real), intent(in) :: x, y, z, v
      real(amrex_real) :: x1, y1, z0
      x1 = x*(x-y)
      y1 = z*(z-y)
      z0 = (z+x-y)**2
      f1qex = -2d0*sqrt(x1*y1)/y1*(g0(y1,x1,z0) + 9d0*x1*y1 + 6d0*v*z0*(x1 - z0))
   end function f1qex

   real(amrex_real) pure function f1qez(x, y, z, v)
      real(amrex_real), intent(in) :: x, y, z, v
      real(amrex_real) :: x2, y2, z0
      x2 = x*z
      y2 = (x-y)*(z-y)
      z0 = (z+x-y)**2
      f1qez = -2d0*sqrt(x2*y2)/y2*(g0(y2,x2,z0) + 9d0*x2*y2 + (1d0 - 6d0*v)*z0*(x2 - z0))
   end function f1qez

   real(amrex_real) pure function f1qexbis(x, y, z, v, sqrt_t23)
      real(amrex_real), intent(in) :: x, y, z, v, sqrt_t23
      real(amrex_real) :: x1, y1, z0
      x1 = x*(x-y)
      y1 = z*(z-y)
      z0 = (z+x-y)**2
      f1qexbis = -2d0*sqrt_t23/y1*( g0(y1,x1,z0) + 9d0*x1*y1 + 6d0*v*z0*(x1 - z0))
   end function f1qexbis

   real(amrex_real) pure function f1qezbis(x, y, z, v, sqrt_t23)
      real(amrex_real), intent(in) :: x, y, z, v, sqrt_t23
      real(amrex_real) :: x2, y2, z0
      x2 = x*z
      y2 = (x-y)*(z-y)
      z0 = (z+x-y)**2
      f1qezbis = -2d0*sqrt_t23/y2*(g0(y2,x2,z0) + 9d0*x2*y2 + (1d0 - 6d0*v)*z0*(x2 - z0))
   end function f1qezbis

   real(amrex_real) pure function f2qex(x, y, z)
      real(amrex_real), intent(in) :: x, y, z
      real(amrex_real) :: x1, y1, z0
      x1 = x*(x-y)
      y1 = z*(z-y)
      z0 = (z+x-y)**2
      f2qex = g1(x1,y1,z0)
   end function f2qex

   real(amrex_real) pure function f2qez(x, y, z)
      real(amrex_real), intent(in) :: x, y, z
      real(amrex_real) :: x2, y2, z0
      x2 = x*z
      y2 = (x-y)*(z-y)
      z0 = (z+x-y)**2
      f2qez = -g1(x2,y2,z0)
   end function f2qez

   pure subroutine quadedge_analytic_derivatives(trigo, a, b, c, d, volume, derivative_theta, derivative_phi)
      use mod_cg3_points, only: cg3_cross_product
      real(amrex_real), dimension(4), intent(in) :: trigo
      real(amrex_real), dimension(3), intent(in) :: a,b, c, d
      real(amrex_real), intent(in) :: volume
      real(amrex_real), dimension(3), intent(out) :: derivative_theta, derivative_phi

      real(amrex_real), dimension(3) ::  p1, p2, p3, p4, p5, xg, dtheta_n, dphi_n
      real(amrex_real), dimension(2) :: p1t, p2t, p3t, p4t, p5t, xgt
      real(amrex_real) :: cos_theta, sin_theta, cos_phi, sin_phi
      real(amrex_real) :: integralt_xy, integralt_xx, integralt_yy
      real(amrex_real) ::  s1, s2, s, i11t, i12t, i22t

      s1 = norm2(cg3_cross_product(b - a, d - a))/2d0
      s2 = norm2(cg3_cross_product(d - a, c - a))/2d0
      s = s1 + s2

      xg = (a + d + (s1*b + s2*c)/s)/3d0

      cos_theta = trigo(c_theta)
      sin_theta = trigo(s_theta)
      cos_phi = trigo(c_phi)
      sin_phi = trigo(s_phi)

      dtheta_n = [-sin_theta, cos_theta, 0d0]
      dphi_n = [cos_theta*cos_phi, sin_theta*cos_phi, -sin_phi]

      p1 = (a + b)/2d0
      p2 = (b + d)/2d0
      p3 = (d + a)/2d0
      p4 = (d + c)/2d0
      p5 = (c + a)/2d0

      xgt = [dot_product(xg, dtheta_n), dot_product(xg, dphi_n)]

      p1t = [dot_product(p1, dtheta_n), dot_product(p1, dphi_n)]
      p2t = [dot_product(p2, dtheta_n), dot_product(p2, dphi_n)]
      p3t = [dot_product(p3, dtheta_n), dot_product(p3, dphi_n)]
      p4t = [dot_product(p4, dtheta_n), dot_product(p4, dphi_n)]
      p5t = [dot_product(p5, dtheta_n), dot_product(p5, dphi_n)]

      integralt_xx = (s1*(p1t(1)**2 + p2t(1)**2) + s*p3t(1)**2 + s2*(p4t(1)**2 + p5t(1)**2))/3d0
      integralt_yy = (s1*(p1t(2)**2 + p2t(2)**2) + s*p3t(2)**2 + s2*(p4t(2)**2 + p5t(2)**2))/3d0
      integralt_xy = (s1*(p1t(1)*p1t(2) + p2t(1)*p2t(2)) + s*p3t(1)*p3t(2) + s2*(p4t(1)*p4t(2) + p5t(1)*p5t(2)))/3d0

      i11t = s*xgt(1)*xgt(1) - integralt_xx
      i12t = s*xgt(1)*xgt(2) - integralt_xy
      i22t = s*xgt(2)*xgt(2) - integralt_yy

      derivative_theta = sin_phi/volume*(i11t*dtheta_n + i12t*dphi_n)
      derivative_phi = (i12t*dtheta_n + i22t*dphi_n)/volume
   end subroutine quadedge_analytic_derivatives

   ! map x in [0, 2π[
   real(amrex_real) pure function modulo_tau(x) result(r)
      real(amrex_real), intent(in) :: x

      ! 10 is an arbitrary value
      if (x < 10d0*TAU .or. x > 10d0*TAU) then
         r = modulo(x, TAU)
      else
         r = x

         do while (r >= TAU)
            r = r - TAU
         end do

         do while (r < 0)
            r = r + TAU
         end do
      end if
   end function modulo_tau

end module mod_mof3d_tetra_analytic_centroid
