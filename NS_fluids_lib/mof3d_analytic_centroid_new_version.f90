!This file is part of Notus 0.6.0

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

module mod_mof3d_analytic_centroid
   implicit none

   private

   double precision, parameter :: PI_4 = acos(0d0)/2d0
   double precision, parameter :: PI_2 = acos(0d0)
   double precision, parameter :: PI_3_4 = 3d0*acos(0d0)/2d0
   double precision, parameter :: PI = 2d0*acos(0d0)

   integer, parameter :: C_THETA = 1
   integer, parameter :: S_THETA = 2
   integer, parameter :: C_PHI   = 3
   integer, parameter :: S_PHI   = 4

   public :: mof3d_compute_analytic_gradient, mof3d_compute_residual_analytic
   public :: mof3d_transform_angles_to_reference, mof3d_transform_angles_to_original

   ! Routines hierarchy
   ! ------------------
   !
   ! mof3d_compute_analytic_gradient                  → Compute gradient from input angles [public]
   ! mof3d_compute_residual_analytic                  → Compute the residual and its jacobian [public]
   ! mof3d_transform_angles_to_reference              → Transform the angles to the reference configuration [public]
   ! mof3d_transform_angles_to_original               → Transform the angles to the original configuration [public]
   ! ├── mof3d_transform_to_reference_map             → Transform input angles to local chart
   ! │   └── mof3d_is_point_inside_map                → Check if the couple of angles belongs to the local chart
   ! ├── mof3d_compute_analytic_derivatives_reference → Compute the partial derivatives on the local chart
   ! │   ├── mof3d_derivatives_quad_face_*            → Compute partial derivatives on the QuadFace sub-regions of the local chart
   ! │   ├── mof3d_derivatives_triangle               → Compute partial derivatives on the Triangle sub-region of the local chart
   ! │   ├── mof3d_derivatives_quad_edge_*            → Compute partial derivatives on the QuadEdge sub-region of the local chart
   ! │   ├── mof3d_derivatives_penta_*                → Compute partial derivatives on the Penta sub-regions of the local chart
   ! │   └── mof3d_derivatives_hexa                   → Compute partial derivatives on the Hexa sub-region of the local chart
   ! └── mof3d_transform_gradient_to_global           → Transform the gradient from the local chart to the global map

contains

   !> Compute the centroid and the gradient of the objective function in rectangular hexahedral cell.
   !!
   !! @param[in]  s_angles:      Spherical angles in the unit cube.
   !! @param[in]  ref_centroid1: Coordinates of the reference centroid of the material 1.
   !! @param[in]  ref_centroid2: Coordinates of the reference centroid of the material 2.
   !! @param[in]  ref_volume:    Reference volume.
   !! @param[in]  c:             Dimensions of the cell.
   !! @param[out] objective:     Objective function.
   !! @param[out] gradient:      Gradient of the objective function.
   !! @ingroup moment_of_fluid
   pure subroutine mof3d_compute_analytic_gradient(s_angles, ref_centroid1, ref_centroid2, ref_volume, c, objective, gradient)
      use variables_mof, only: mof_use_symmetric_reconstruction
      double precision, dimension(2), intent(in) :: s_angles
      double precision, dimension(3), intent(in) :: ref_centroid1, ref_centroid2
      double precision, intent(in) :: ref_volume
      double precision, dimension(3), intent(in) :: c
      double precision, intent(out) :: objective
      double precision, dimension(2), intent(out) :: gradient

      double precision, dimension(6) :: residual
      double precision, dimension(6,2) :: jacobian
      integer :: rows

      rows = merge(6, 3, mof_use_symmetric_reconstruction)

      call mof3d_compute_residual_analytic( &
         &    s_angles,                     &
         &    ref_centroid1,                &
         &    ref_centroid2,                &
         &    ref_volume,                   &
         &    c,                            &
         &    residual(1:rows),             &
         &    jacobian(1:rows,:)            &
         & )

      objective = dot_product(residual(1:rows), residual(1:rows))

      gradient = [                                      &
         &    2d0*dot_product(residual(1:rows), jacobian(1:rows,1)), &
         &    2d0*dot_product(residual(1:rows), jacobian(1:rows,2))  &
         & ]
   end subroutine mof3d_compute_analytic_gradient

   !> Transform the spherical angles in the original configuration to the reference configuration.
   !!
   !! @param[in]  transformation: Tranformation matrix.
   !! @param[in]  orig_angles:    Angles in the original configuration.
   !! @param[out] ref_angles:     Angles in the reference configuration.
   !! @ingroup moment_of_fluid
   pure subroutine mof3d_transform_angles_to_reference(transformation, orig_angles, ref_angles)
      use mod_cg3_points, only: cg3_direction_to_spherical_angles, cg3_spherical_angles_to_direction
      double precision, dimension(3), intent(in) :: transformation
      double precision, dimension(2), intent(in) :: orig_angles
      double precision, dimension(2), intent(out) :: ref_angles

      double precision, dimension(3) :: direction

      call cg3_spherical_angles_to_direction(orig_angles, direction)

      direction = transformation*direction
      direction = direction/norm2(direction)

      call cg3_direction_to_spherical_angles(direction, ref_angles)
   end subroutine mof3d_transform_angles_to_reference

   !> Transform the spherical angles in the reference configuration to the original configuration.
   !!
   !! @param[in]  transformation: Tranformation matrix.
   !! @param[in]  ref_angles:     Angles in the reference configuration.
   !! @param[out] orig_angles:    Angles in the original configuration.
   !! @ingroup moment_of_fluid
   pure subroutine mof3d_transform_angles_to_original(transformation, ref_angles, orig_angles)
      use mod_cg3_points, only: cg3_direction_to_spherical_angles, cg3_spherical_angles_to_direction
      double precision, dimension(3), intent(in) :: transformation
      double precision, dimension(2), intent(in) :: ref_angles
      double precision, dimension(2), intent(out) :: orig_angles

      double precision, dimension(3) :: direction

      call cg3_spherical_angles_to_direction(ref_angles, direction)

      direction = direction/transformation
      direction = direction/norm2(direction)

      call cg3_direction_to_spherical_angles(direction, orig_angles)
   end subroutine mof3d_transform_angles_to_original

   !> Compute the residual and its jacobian in a rectangular hexahedral cell.
   !!
   !! @param[in]  s_angles:      Spherical angles in the unit cube.
   !! @param[in]  ref_centroid1: Coordinates of the reference centroid of the material 1.
   !! @param[in]  ref_centroid2: Coordinates of the reference centroid of the material 2.
   !! @param[in]  ref_volume:    Reference volume of the material 1.
   !! @param[in]  c:             Dimensions of the cell.
   !! @param[out] residual:      Residual function.
   !! @param[out] jacobian:      Jacobian of the residual function.
   !! @ingroup moment_of_fluid
   pure subroutine mof3d_compute_residual_analytic(s_angles, ref_centroid1, ref_centroid2, ref_volume, c, residual, jacobian)
      use variables_mof, only: mof_use_symmetric_reconstruction
      double precision, dimension(2), intent(in) :: s_angles
      double precision, dimension(3), intent(in) :: ref_centroid1, ref_centroid2
      double precision, intent(in) :: ref_volume
      double precision, dimension(3), intent(in) :: c
      double precision, dimension(:), intent(out) :: residual
      double precision, dimension(:,:), intent(out) :: jacobian

      ! Prefixes:
      ! s_: scaled in unit cube.
      ! p_: permutated.

      double precision, dimension(3,2) :: p_s_partial_derivative, s_partial_derivative
      double precision, dimension(3) :: p_s_centroid, s_centroid, centroid
      double precision, dimension(3) :: dual_centroid
      double precision, dimension(2) :: p_s_angles
      double precision :: cell_volume, s_volume
      integer, dimension(3) :: sign_permutation

      ! Compute transformed angle and permutation.
      call mof3d_transform_to_local_chart(s_angles, p_s_angles, sign_permutation)

      ! Compute the volume of the hexahedron.
      cell_volume = product(c)

      ! Compute the reference volume in the reference configuration.
      s_volume = ref_volume/cell_volume

      ! Compute the gradient and the centroid in the local chart.
      call mof3d_compute_analytic_derivatives_local( &
         &    p_s_angles,                            &
         &    s_volume,                              &
         &    p_s_centroid,                          &
         &    p_s_partial_derivative                 &
         & )

      ! Transform the centroid in the orginal configuration.
      call mof3d_transform_point(sign_permutation, p_s_centroid, s_centroid)

      ! Transform the centroid in the permutated original configuration.
      centroid = s_centroid*c

      ! Compute the residual.
      residual(1:3) = centroid - ref_centroid1

      ! Transform the partial derivatives in the permutated orginal configuration.
      call mof3d_transform_vector(sign_permutation, p_s_partial_derivative(:,1), s_partial_derivative(:,1))
      call mof3d_transform_vector(sign_permutation, p_s_partial_derivative(:,2), s_partial_derivative(:,2))

      ! Partial derivatives in the original configuration.
      jacobian(1:3,1) = s_partial_derivative(:,1)*c
      jacobian(1:3,2) = s_partial_derivative(:,2)*c*sign_permutation(3)

      ! Symmetric contribution.
      if (mof_use_symmetric_reconstruction) then
         ! Compute the dual centroid.
         dual_centroid = (cell_volume*c/2d0 - ref_volume*centroid)/(cell_volume - ref_volume)
         ! Compute the residual.
         residual(4:6) = dual_centroid - ref_centroid2
         ! Compute the jacobian.
         jacobian(4:6,:) = ref_volume*jacobian(1:3,:)/(ref_volume - cell_volume)
      end if
   end subroutine mof3d_compute_residual_analytic

   pure subroutine mof3d_transform_to_local_chart(angles, t_angles, sign_permutation)
      double precision, dimension(2), intent(in) :: angles
      double precision, dimension(2), intent(out) :: t_angles
      integer, dimension(3), intent(out) :: sign_permutation

      ! Initialize the sign of the permutation
      sign_permutation = 1

      ! Compute permutation
      t_angles = angles

      ! Ensure that the angles are in their domains of definition:
      ! → θ ∈ [-π,π]
      ! → φ ∈ [0,π]
      !-----------------------------------------------------------

      ! Check if φ belongs to [0,π]
      if (t_angles(2) < 0d0 .or. t_angles(2) > PI) then
         ! First, ensure that φ belongs to [0,2π[
         t_angles(2) = modulo_tau(t_angles(2))

         ! Ιf φ belongs to ]π,2π[, transform to ]0,π[
         ! Otherwise, φ belongs to [0,π]
         if (t_angles(2) > PI) then
            t_angles(2) = 2d0*PI - t_angles(2)
            sign_permutation(3) = -1
         end if
      end if

      ! Check if θ belongs to [-π,π]
      if (t_angles(1) < -PI .or. t_angles(1) > PI) then
         t_angles(1) = modulo_tau(t_angles(1) + PI) - PI
      end if

      ! Crop to [0,π/2] × [0,π/2]
      !--------------------------

      ! Here θ belongs to [-π,π] and φ belongs to [0,π]

      ! Check if φ belongs to [0,π/2]
      if (t_angles(2) > PI_2) then
         ! Transform ]π/2,π] → [0,π/2[
         ! Reflection on third axis
         t_angles(2) = PI - t_angles(2)
         sign_permutation(3) = -sign_permutation(3)
      end if

      ! Check if θ belongs to [0,π/2]
      if (t_angles(1) > PI_2) then
         ! Transform ]π/2,π] → ]0,π/2]
         ! Rotate by -π/2 around third axis
         t_angles(1) = t_angles(1) - PI_2
         sign_permutation(1:2) = [1, -1]
      else if (t_angles(1) < -PI_2) then
         ! Transform [-π,-π/2[ → [0,π/2[
         ! Rotate by π around third axis
         t_angles(1) = t_angles(1) + PI
         sign_permutation(1:2) = [-1, -1]
      else if (t_angles(1) < 0d0) then
         ! Transform [-π/2,0[ → [0,π/2[
         ! Rotate by π/2 around third axis
         t_angles(1) = t_angles(1) + PI_2
         sign_permutation(1:2) = [-1, 1]
      end if
   end subroutine mof3d_transform_to_local_chart

   pure subroutine mof3d_transform_vector(sign_permutation, v, t_v)
      integer, dimension(3), intent(in) :: sign_permutation
      double precision, dimension(3), intent(in) :: v
      double precision, dimension(3), intent(out) :: t_v

      if (sign_permutation(1)*sign_permutation(2) < 0) then
         t_v(1) = sign_permutation(2)*v(2)
         t_v(2) = sign_permutation(1)*v(1)
         t_v(3) = sign_permutation(3)*v(3)
      else
         t_v(1) = sign_permutation(1)*v(1)
         t_v(2) = sign_permutation(2)*v(2)
         t_v(3) = sign_permutation(3)*v(3)
      end if
   end subroutine mof3d_transform_vector

   pure subroutine mof3d_transform_point(sign_permutation, p, t_p)
      integer, dimension(3), intent(in) :: sign_permutation
      double precision, dimension(3), intent(in) :: p
      double precision, dimension(3), intent(out) :: t_p

      if (sign_permutation(1)*sign_permutation(2) < 0) then
         t_p(1) = p(2)
         t_p(2) = p(1)
         t_p(3) = p(3)
      else
         t_p = p
      end if

      if (sign_permutation(1) == -1) t_p(2) = 1d0 - t_p(2)
      if (sign_permutation(2) == -1) t_p(1) = 1d0 - t_p(1)
      if (sign_permutation(3) == -1) t_p(3) = 1d0 - t_p(3)
   end subroutine mof3d_transform_point

   !                Volume fraction < 1/6                                    Volume fraction > 1/6
   !       ↑ φ                                                       ↑ φ
   !
   !       |                                                         |
   !
   !   A → +--------+-----------------------------+  -  → θ      A → +------------------+---------------+  -  → θ
   !       |Quad- .' \                           /                   |                .' \             /
   !       |Face.'    \                         /                    |              .'    \ QuadEdge  /
   !       |  .'       \                       /                     |  QuadFace  .'       \         /
   !       |.'  Penta   \      QuadEdge       /                      |          .'          \       /
   !   B → +-._          \                   /                       |        .'             \     /
   !           `-._       \                 /                        |      .'                \   /
   !       |       `-._    \               /                         |    .'                   \ /
   !                   `-._ \             /                      B → |  .'        Penta         X
   !   C → |               `-._         _.                           |.'                       / \
   !   D →                    \`-._ _.-'/                        C → +-._                     /   \
   !       |                   \   '   /                                 `-._                /     \
   !                            \     / Triangle                     |       `-._           /       \
   !       |                     \   /                                           `-._      /  Hexa   \
   !                              \ /                            D → |               `-.../           \
   !   E → |                       +                                                      `-._     _.-'
   !                                                                 |                        `---'
   !       |
   !       ↑        ↑        ↑     ↑     ↑        ↑                  ↑                  ↑ ↑     ↑     ↑ ↑
   !       0        1        2     3     4        5                  0                  1 2     3     4 5
   !
   !

   pure subroutine mof3d_compute_analytic_derivatives_local(angles, volume, centroid, derivative)
      double precision, dimension(2), intent(in) :: angles
      double precision, intent(in) :: volume
      double precision, dimension(3), intent(out) :: centroid
      double precision, dimension(3,2), intent(out) :: derivative

      if (6d0*volume < 1d0) then
         call mof3d_compute_analytic_derivatives_below_one_sixth(angles, volume, centroid, derivative)
      else if (volume <= 0.5d0 - 1d1*epsilon(1d0)) then
         call mof3d_compute_analytic_derivatives_above_one_sixth(angles, volume, centroid, derivative)
      else ! volume = 1/2
         call mof3d_compute_analytic_derivatives_half(angles, volume, centroid, derivative)
      end if
   end subroutine mof3d_compute_analytic_derivatives_local

   pure subroutine mof3d_compute_analytic_derivatives_below_one_sixth(angles, volume, centroid, derivative)
      double precision, dimension(2), intent(in) :: angles
      double precision, intent(in) :: volume
      double precision, dimension(3), intent(out) :: centroid
      double precision, dimension(3,2), intent(out) :: derivative

      double precision, dimension(4) :: trigo
      double precision :: l13, l23, l12, t1, t2
      double precision :: theta, phi, theta3
      double precision :: tan_theta, cos_theta, sin_theta
      double precision :: limit_penta_left1, limit_penta_left2, limit_penta_left3
      double precision :: limit_penta_right1, limit_penta_right2, limit_penta_right3
      double precision :: limit_penta_below1, limit_penta_below2, limit_penta_below3
      double precision :: limit_triangle1, limit_triangle2, limit_triangle3

      ! ℓij = 2V
      l23 = 2d0*volume
      l13 = 2d0*volume
      l12 = 2d0*volume

      ! Aliases for θ and φ
      theta = angles(1)
      phi = angles(2)

      ! Compute trigonometric functions
      trigo(C_THETA) = cos(theta)
      trigo(S_THETA) = sin(theta)
      trigo(C_PHI)   = cos(phi)
      trigo(S_PHI)   = sin(phi)

      ! Aliases of trigonometric functions
      cos_theta = trigo(C_THETA)
      sin_theta = trigo(S_THETA)
      tan_theta = sin_theta/cos_theta

      ! T1 and T2
      t1 = 1d0/tan_theta
      t2 = tan_theta

      ! θ3
      theta3 = PI_4

      ! Left of θ3
      if (theta <= theta3) then
         ! Left of θ1
         if (theta <= atan(l23)) then
            ! φlim_pl3
            limit_penta_left3 = PI_2 - atan((l23 + t2 + sqrt((l23 + t2)**2 - 4d0/3d0*t2**2))/2d0*cos_theta)

            ! Below φlim_pl3
            if (phi <= limit_penta_left3) then
               ! φlim_pb1
               limit_penta_below1 = PI_2 - atan((sin_theta + cos_theta)/l12)

               ! Below φlim_pb1
               if (phi <= limit_penta_below1) then
                  call mof3d_derivatives_quad_face_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
               else ! Above φlim_pb1
                  ! φlim_pb3
                  limit_penta_below3 = PI_2 - atan((cos_theta - sin_theta + (sin_theta)**2/(3d0*cos_theta))/l12)

                  ! Below φlim_pb3
                  if (phi <= limit_penta_below3) then
                     call mof3d_derivatives_penta_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pb3
                     call mof3d_derivatives_quad_edge_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if ! φlim_pb3
               end if ! φlim_pb1
            else ! Above φlim_pl3
               ! φlim_pl1
               limit_penta_left1 = PI_2 - atan((l23 - t2)*cos_theta)

               ! Below φlim_pl1
               ! Note: the inequality must be strict here
               if (phi < limit_penta_left1) then
                  call mof3d_derivatives_penta_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
               else ! Above φlim_pl1
                  call mof3d_derivatives_quad_face_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
               end if ! φlim_pl1
            end if ! φlim_pl3
         else if (theta <= atan(3d0*l23)) then ! Left of θ2t and right of θ1
            ! φlim_pl3
            limit_penta_left3 = PI_2 - atan(0.5d0*(l23 + t2 + sqrt((l23 + t2)**2 - 4d0/3d0*t2**2))*cos_theta)

            ! Below φlim_pl3
            if (phi <= limit_penta_left3) then
               ! φlim_pb1
               limit_penta_below1 = PI_2 - atan((sin_theta + cos_theta)/l12)

               ! Below φlim_pb1
               if (phi <= limit_penta_below1) then
                  call mof3d_derivatives_quad_face_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
               else ! Above φlim_pb1
                  ! φlim_pb3
                  limit_penta_below3 = PI_2 - atan((cos_theta - sin_theta + sin_theta**2/(3d0*cos_theta))/l12)

                  ! Below φlim_pb3
                  if (phi <= limit_penta_below3) then
                     call mof3d_derivatives_penta_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pb3
                     call mof3d_derivatives_quad_edge_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if ! φlim_pb3
               end if ! φlim_pb1
            else ! Above φlim_pl3
               ! φlim_pl2
               limit_penta_left2 = PI_2 - atan((t2 - (-t2 + sqrt(12d0*l23*t2 - 3d0*t2**2))/2d0)*cos_theta)

               ! Below φlim_pl2
               ! Note: the inequality must be strict here
               if (phi < limit_penta_left2) then
                  call mof3d_derivatives_penta_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
               else ! Above φlim_pl2
                  call mof3d_derivatives_quad_edge_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
               end if ! φlim_pl2
            end if
         else ! Right of θ2t and left of θ3
            ! φlim_t1
            limit_triangle1 = PI_2 - atan(t2**2/(3d0*l23)*cos_theta)

            ! Below φlim_t1
            if (phi <= limit_triangle1) then
               ! φlim_pb1
               limit_penta_below1 = PI_2 - atan((sin_theta + cos_theta)/l12 )

               ! Below φlim_pb1
               if (phi <= limit_penta_below1) then
                  call mof3d_derivatives_quad_face_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
               else ! Above φlim_pb1
                  ! φlim_pb3
                  limit_penta_below3 = PI_2 - atan((cos_theta - sin_theta + sin_theta**2/(3d0*cos_theta))/l12)

                  ! Below φlim_pb3
                  if (phi <= limit_penta_below3) then
                     call mof3d_derivatives_penta_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pb3
                     call mof3d_derivatives_quad_edge_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if ! φlim_pb3
               end if ! φlim_pb
            else ! Above φlim_t1
               ! φlim_t3
               limit_triangle3 = PI_2 - atan(sqrt(6d0*volume*tan_theta)*cos_theta)

               ! Below φlim_t3
               if (phi <= limit_triangle3) then
                  call mof3d_derivatives_triangle(trigo, volume, derivative(:,1), derivative(:,2), centroid)
               else ! Above φlim_t3
                  call mof3d_derivatives_quad_edge_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
               end if ! φlim_t3
            end if ! φlim_t1
         end if ! 0 < θ ≤ θ3
      else ! Right of θ3
         ! Left of θ4t
         if (theta <= atan(1d0/(3d0*l13))) then
            ! φlim_t2
            limit_triangle2 = PI_2 - atan(t1**2/(3d0*l13)*sin_theta)

            ! Below φlim_t2
            if (phi <= limit_triangle2) then
               ! φlim_pb1
               limit_penta_below1 = PI_2 - atan((sin_theta + cos_theta)/l12)

               ! Below φlim_pb1
               if (phi <= limit_penta_below1) then
                  call mof3d_derivatives_quad_face_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
               else ! Above φlim_pb1
                  ! φlim_pb2
                  limit_penta_below2 = PI_2 - atan((sin_theta - cos_theta + cos_theta**2/(3d0*sin_theta))/l12)

                  ! Below φlim_pb2
                  if (phi <= limit_penta_below2) then
                     call mof3d_derivatives_penta_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pb2
                     call mof3d_derivatives_quad_edge_right(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if ! φlim_pb2
               end if ! φlim_pb1
            else ! Above φlim_t2
               ! φlim_t3
               limit_triangle3 = PI_2 - atan(sqrt(6d0*volume*tan_theta)*cos_theta)

               ! Below φlim_t3
               if (phi <= limit_triangle3) then
                  call mof3d_derivatives_triangle(trigo, volume, derivative(:,1), derivative(:,2), centroid)
               else ! Above φlim_t3
                  call mof3d_derivatives_quad_edge_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
               end if ! φlim_t3
            end if ! φlim_t2
         else if (theta <= atan(1d0/l13)) then ! θ4t < θ ≤ θ5
            ! φlim_pr3
            limit_penta_right3 = PI_2 - atan((l13 + t1 + sqrt((l13 + t1)**2 - 4d0/3d0*t1**2))/2d0*sin_theta)

            ! Below φlim_pr3
            if (phi <= limit_penta_right3) then
               ! φlim_pb1
               limit_penta_below1 = PI_2 - atan((sin_theta + cos_theta)/l12)

               ! Below φlim_pb1
               if (phi <= limit_penta_below1) then
                  call mof3d_derivatives_quad_face_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
               else ! Above φlim_pb1
                  ! φlim_pb2
                  limit_penta_below2 = PI_2 - atan((sin_theta - cos_theta + cos_theta**2/(3d0*sin_theta))/l12)

                  ! Below φlim_pb2
                  if (phi <= limit_penta_below2) then
                     call mof3d_derivatives_penta_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pb2
                     call mof3d_derivatives_quad_edge_right(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if ! φlim_pb2
               end if ! φlim_pb1
            else ! Above φlim_pr3
               ! φlim_pr2
               limit_penta_right2 = PI_2 - atan((t1 - (-t1 + sqrt(12d0*l13*t1 - 3d0*t1**2))/2d0)*sin_theta)

               ! Below φlim_pr2
               ! Note: the inequality must be strict here
               if (phi < limit_penta_right2) then
                  call mof3d_derivatives_penta_right(trigo, volume, derivative(:,1), derivative(:,2), centroid)
               else ! Above φlim_pr2
                  call mof3d_derivatives_quad_edge_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
               end if ! φlim_pr2
            end if ! φlim_pr3
         else ! θ5 < θ ≤ π/2
            ! Below φlim_pr3
            limit_penta_right3 = PI_2 - atan((l13 + t1 + sqrt((l13 + t1)**2 - 4d0/3d0*t1**2))/2d0*sin_theta)

            ! Below φlim_pr3
            if (phi <= limit_penta_right3) then
               ! φlim_pb1
               limit_penta_below1 = PI_2 - atan((sin_theta + cos_theta)/l12)

               ! Below φlim_pb1
               if (phi <= limit_penta_below1) then
                  call mof3d_derivatives_quad_face_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
               else ! Above φlim_pb1
                  ! φlim_pb2
                  limit_penta_below2 = PI_2 - atan((sin_theta - cos_theta + cos_theta**2/(3d0*sin_theta))/l12)

                  ! Below φlim_pb2
                  if (phi <= limit_penta_below2) then
                     call mof3d_derivatives_penta_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pb2
                     call mof3d_derivatives_quad_edge_right(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if ! φlim_pb2
               end if ! φlim_pb1
            else ! φlim_pr3
               ! φlim_pr1
               limit_penta_right1 = PI_2 - atan((l13 - t1)*sin_theta)

               ! Below φlim_pr1
               ! Note: the inequality must be strict here
               if (phi < limit_penta_right1) then
                  call mof3d_derivatives_penta_right(trigo, volume, derivative(:,1), derivative(:,2), centroid)
               else ! Above φlim_pr1
                  call mof3d_derivatives_quad_face_right(trigo, volume, derivative(:,1), derivative(:,2), centroid)
               end if ! φlim_pr1
            end if ! φlim_pr3
         end if ! θ3 < θ ≤ π/2
      end if ! θ ≤ θ3 ?
   end subroutine mof3d_compute_analytic_derivatives_below_one_sixth

   pure subroutine mof3d_compute_analytic_derivatives_above_one_sixth(angles, volume, centroid, derivative)
      double precision, dimension(2), intent(in) :: angles
      double precision, intent(in) :: volume
      double precision, dimension(3), intent(out) :: centroid
      double precision, dimension(3,2), intent(out) :: derivative

      double precision, dimension(4) :: trigo
      double precision :: l13, l23, l12, t1, t2, sqrthexa
      double precision :: theta, phi, theta3, theta1, theta2h
      double precision :: tan_theta, cos_theta, sin_theta
      double precision :: limit_penta_left1, limit_penta_left2, limit_penta_left3
      double precision :: limit_penta_right1, limit_penta_right2, limit_penta_right3
      double precision :: limit_penta_below1, limit_penta_below2, limit_penta_below3
      double precision :: limit_hexa1, limit_hexa2, limit_hexa3

      ! ℓij = 2V
      l23 = 2d0*volume
      l13 = 2d0*volume
      l12 = 2d0*volume

      ! Aliases for θ and φ
      theta = angles(1)
      phi = angles(2)

      ! Compute trigonometric functions
      trigo(C_THETA) = cos(theta)
      trigo(S_THETA) = sin(theta)
      trigo(C_PHI)   = cos(phi)
      trigo(S_PHI)   = sin(phi)

      ! Aliases of trigonometric functions
      cos_theta = trigo(C_THETA)
      sin_theta = trigo(S_THETA)
      tan_theta = sin_theta/cos_theta

      ! T1 and T2
      t1 = 1d0/tan_theta
      t2 = tan_theta

      ! θ1
      theta1 = atan(l23)
      ! θ2h
      theta2h = atan((3d0 - sqrt(12d0*l23 - 3d0))/2d0)
      ! θ3
      theta3 = PI_4

      ! Case θ1 ≤ θ2h (implies θ4h ≤ θ5)
      if (theta1 <= theta2h) then
         ! θ ≤ θ3
         if (theta <= theta3) then
            ! θ ≤ θ1
            if (theta <= theta1) then
               ! φlim_pl3
               limit_penta_left3 = PI_2 - atan((l23 + t2 + sqrt((l23 + t2)**2 - 4d0/3d0*t2**2))/2d0*cos_theta)

               ! Below φlim_pl3
               if (phi <= limit_penta_left3) then
                  ! φlim_pb1
                  limit_penta_below1 = PI_2 - atan((sin_theta + cos_theta)/l12)

                  ! Below φlim_pb1
                  if (phi <= limit_penta_below1) then
                     call mof3d_derivatives_quad_face_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pl1
                     ! φlim_pb3
                     limit_penta_below3 = PI_2 - atan((cos_theta - sin_theta + sin_theta**2/(3d0*cos_theta))/l12)

                     ! Below φlim_pb3
                     if (phi < limit_penta_below3) then
                        call mof3d_derivatives_penta_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! Above φlim_pb3
                        call mof3d_derivatives_quad_edge_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if ! φlim_pb3
                  end if ! φlim_pl1
               else ! Above φlim_pl3
                  ! φlim_pl1
                  limit_penta_left1 = PI_2 - atan((l23 - t2)*cos_theta)

                  ! Below φlim_pl1
                  ! Note: the inequality must be strict here
                  if (phi < limit_penta_left1) then
                     call mof3d_derivatives_penta_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pl1
                     call mof3d_derivatives_quad_face_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if ! φlim_pl1
               end if ! φlim_pl3
            else if (theta <= theta2h) then ! θ1 < θ ≤ θ2h
               ! φlim_pl3
               limit_penta_left3 = PI_2 - atan((l23 + t2 + sqrt((l23 + t2)**2 - 4d0/3d0*t2**2))/2d0*cos_theta)

               ! Below φlim_pl3
               if (phi <= limit_penta_left3) then
                  ! φlim_pb1
                  limit_penta_below1 = PI_2 - atan((sin_theta + cos_theta)/l12)

                  ! Below φlim_pb1
                  if (phi <= limit_penta_below1) then
                     call mof3d_derivatives_quad_face_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pb1
                     ! φlim_pb3
                     limit_penta_below3 = PI_2 - atan((cos_theta - sin_theta + sin_theta**2/(3d0*cos_theta))/l12)

                     ! Below φlim_pb3
                     if (phi < limit_penta_below3) then
                        call mof3d_derivatives_penta_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! Above φlim_pb3
                        call mof3d_derivatives_quad_edge_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if ! φlim_pb3
                  end if ! φlim_pb1
               else ! Above φlim_pl3
                  ! φlim_pl2
                  limit_penta_left2 = PI_2 - atan((t2 - (-t2 + sqrt(12d0*l23*t2 - 3d0*t2**2))/2d0)*cos_theta)

                  ! Below φlim_pl2
                  ! Note: the inequality must be strict here
                  if (phi < limit_penta_left2) then
                     call mof3d_derivatives_penta_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pl2
                     call mof3d_derivatives_quad_edge_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if ! φlim_pl2
               end if ! φlim_pl3
            else ! θ2h < θ ≤ θ3
               ! φlim_h1
               sqrthexa = sqrt(l23*t2)
               limit_hexa1 = PI_2 - atan((1d0 - 2d0*sqrthexa*cos((acos((3d0*(1d0 - l23 - t2) &
                  &        + t2*t2)/(2d0*l23*sqrthexa)) + 4d0*PI)/3d0))*cos_theta)

               ! Below φlim_h1
               if (phi < limit_hexa1) then
                  ! φlim_pb1
                  limit_penta_below1 = PI_2 - atan((sin_theta + cos_theta)/l12)

                  ! Below φlim_pb1
                  if (phi <= limit_penta_below1) then
                     call mof3d_derivatives_quad_face_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pb1
                     ! φlim_h3
                     sqrthexa = sqrt((2d0 - l23)*t2)
                     limit_hexa3 = PI_2 - atan((1d0 + t2 + 2d0*sqrthexa*cos((acos(3d0*(1d0 - l23)*(1d0 + t2) &
                        &        / (2d0*(2d0 - l23)*sqrthexa)) + 4d0*PI)/3d0))*cos_theta)

                     ! Below φlim_h3
                     if (phi < limit_hexa3) then
                        call mof3d_derivatives_penta_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! Above φlim_h3
                        call mof3d_derivatives_hexa(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if ! φlim_h3
                  end if ! φlim_pb1
               else ! Above φlim_h1
                  ! φlim_pl2
                  limit_penta_left2 = PI_2 - atan((t2 - (-t2 + sqrt(12d0*l23*t2 - 3d0*t2**2))/2d0)*cos_theta)

                  ! Below φlim_pl2
                  ! Note: the inequality must be strict here
                  if (phi < limit_penta_left2) then
                     call mof3d_derivatives_penta_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pl2
                     call mof3d_derivatives_quad_edge_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if ! φlim_pl2
               end if ! φlim_h1
            end if ! Left of θ3 cases
         else ! Right of θ3
            ! θ ≤ θ4h
            if (theta < PI_2 - atan((3d0 - sqrt(12d0*l13 - 3d0))/2d0)) then
               ! φlim_h2
               sqrthexa = sqrt(l13*t1)
               limit_hexa2 = PI_2 - atan((1d0 - 2d0*sqrthexa*cos((acos((3d0*(1d0 - l13 - t1) + t1*t1) &
                  &        / (2d0*l13*sqrthexa)) + 4d0*PI)/3d0))*sin_theta)

               ! Below φlim_h2
               if (phi < limit_hexa2) then
                  ! φlim_pb1
                  limit_penta_below1 = PI_2 - atan((sin_theta + cos_theta)/l12)

                  ! Below φlim_pb1
                  if (phi <= limit_penta_below1) then
                     call mof3d_derivatives_quad_face_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pb1
                     ! φlim_h3
                     sqrthexa = sqrt((2d0 - l23)*t2)
                     limit_hexa3 = PI_2 - atan((1d0 + t2 + 2d0*sqrthexa*cos((acos(3d0*(1d0 - l23)*(1d0 + t2) &
                        &        / (2d0*(2d0 - l23)*sqrthexa)) + 4d0*PI)/3d0))*cos_theta)

                     ! Below φlim_h3
                     if (phi <= limit_hexa3) then
                        call mof3d_derivatives_penta_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! Above φlim_h3
                        call mof3d_derivatives_hexa(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if ! φlim_h3
                  end if ! φlim_pb1
               else ! Above φlim_h2
                  ! φlim_pr2
                  limit_penta_right2 = PI_2 - atan((t1 - (-t1 + sqrt(12d0*l13*t1 - 3d0*t1**2))/2d0)*sin_theta)

                  ! Below φlim_pr2
                  ! Note: the inequality must be strict here
                  if (phi < limit_penta_right2) then
                     call mof3d_derivatives_penta_right(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pr2
                     call mof3d_derivatives_quad_edge_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if ! φlim_pr2
               end if ! φlim_h2
            else if (theta <= atan(1d0/l13)) then ! θ4h < θ ≤ θ5
               ! φlim_pr3
               limit_penta_right3 = PI_2 - atan((l13 + t1 + sqrt((l13 + t1)**2 - 4d0/3d0*t1**2))/2d0*sin_theta)

               ! Below φlim_pr3
               if (phi <= limit_penta_right3) then
                  ! φlim_pb1
                  limit_penta_below1 = PI_2 - atan((sin_theta + cos_theta)/l12)

                  ! Below φlim_pb1
                  if (phi <= limit_penta_below1) then
                     call mof3d_derivatives_quad_face_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pb1
                     ! φlim_pb2
                     limit_penta_below2 = PI_2 - atan((sin_theta - cos_theta + cos_theta**2/(3d0*sin_theta))/l12)

                     ! Above φlim_pb2
                     if (phi < limit_penta_below2) then
                        call mof3d_derivatives_penta_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! Below φlim_pb2
                        call mof3d_derivatives_quad_edge_right(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if ! φlim_pb2
                  end if ! φlim_pb1
               else ! Above φlim_pr3
                  ! φlim_pr2
                  limit_penta_right2 = PI_2 - atan((t1 - (-t1 + sqrt(12d0*l13*t1 - 3d0*t1**2))/2d0)*sin_theta)

                  ! Below φlim_pr2
                  ! Note: the inequality must be strict here
                  if (phi < limit_penta_right2) then
                     call mof3d_derivatives_penta_right(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pr2
                     call mof3d_derivatives_quad_edge_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if ! φlim_pr2
               end if ! φlim_pr3
            else ! θ5 < θ ≤ π/2
               ! φlim_pr3
               limit_penta_right3 = PI_2 - atan((l13 + t1 + sqrt((l13 + t1)**2 - 4d0/3d0*t1**2))/2d0*sin_theta)

               ! Below φlim_pr3
               if (phi <= limit_penta_right3) then
                  ! φlim_pb1
                  limit_penta_below1 = PI_2 - atan((sin_theta + cos_theta)/l12)

                  ! Below φlim_pb1
                  if (phi <= limit_penta_below1) then
                     call mof3d_derivatives_quad_face_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pb1
                     ! φlim_pb2
                     limit_penta_below2 = PI_2 - atan((sin_theta - cos_theta + cos_theta**2/(3d0*sin_theta))/l12)

                     ! Below φlim_pb2
                     if (phi < limit_penta_below2) then
                        call mof3d_derivatives_penta_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! Above φlim_pb2
                        call mof3d_derivatives_quad_edge_right(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if ! φlim_pb2
                  end if ! φlim_pb1
               else ! Above φlim_pr3
                  ! φlim_pr1
                  limit_penta_right1 = PI_2 - atan((l13 - t1)*sin_theta)

                  ! Below φlim_pr1
                  ! Note: the inequality must be strict here
                  if (phi < limit_penta_right1) then
                     call mof3d_derivatives_penta_right(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pr1
                     call mof3d_derivatives_quad_face_right(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if ! φlim_pr1
               end if ! φlim_pr3
            end if ! Right of θ3
         end if ! θ3
      else ! Case θ1 > θ2h (implies θ4h > θ5)
         ! Left of θ3
         if (theta <= theta3) then
            ! Left of θ2h
            if (theta <= theta2h) then
               ! φlim_pl3
               limit_penta_left3 = PI_2 - atan((l23 + t2 + sqrt((l23 + t2)**2 - 4d0/3d0*t2**2))/2d0*cos_theta)

               ! Below φlim_pl3
               if (phi <= limit_penta_left3) then
                  ! φlim_pb1
                  limit_penta_below1 = PI_2 - atan((sin_theta + cos_theta)/l12)

                  ! Below φlim_pb1
                  if (phi <= limit_penta_below1) then
                     call mof3d_derivatives_quad_face_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pb1
                     ! φlim_pb3
                     limit_penta_below3 = PI_2 - atan((cos_theta - sin_theta + sin_theta**2/(3d0*cos_theta))/l12)

                     ! Below φlim_pb3
                     if (phi < limit_penta_below3) then
                        call mof3d_derivatives_penta_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! Above φlim_pb3
                        call mof3d_derivatives_quad_edge_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if ! φlim_pb3
                  end if ! φlim_pb1
               else ! φlim_pl3
                  ! φlim_pl1
                  limit_penta_left1 = PI_2 - atan((l23 - t2)*cos_theta)

                  ! Below φlim_pl1
                  ! Note: the inequality must be strict here
                  if (phi < limit_penta_left1) then
                     call mof3d_derivatives_penta_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pl1
                     call mof3d_derivatives_quad_face_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if ! φlim_pl1
               end if ! φlim_pl3
            else if (theta <= theta1) then ! θ2h < θ ≤ θ1
               ! φlim_h1
               sqrthexa = sqrt(l23*t2)
               limit_hexa1 = PI_2 - atan((1d0 - 2d0*sqrthexa*cos((acos((3d0*(1d0 - l23 - t2) + t2*t2) &
                  &        / (2d0*l23*sqrthexa)) + 4d0*PI)/3d0))*cos_theta)

               ! Below φlim_h1
               if (phi < limit_hexa1) then
                  ! φlim_pb1
                  limit_penta_below1 = PI_2 - atan((sin_theta + cos_theta)/l12)

                  ! Below φlim_pb1
                  if (phi <= limit_penta_below1) then
                     call mof3d_derivatives_quad_face_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pb1
                     ! φlim_h3
                     sqrthexa = sqrt((2d0 - l23)*t2)
                     limit_hexa3 = PI_2 - atan((1d0 + t2 + 2d0*sqrthexa*cos((acos(3d0*(1d0 - l23)*(1d0 + t2) &
                        &        / (2d0*(2d0 - l23)*sqrthexa)) + 4d0*PI)/3d0))*cos_theta)

                     ! Below φlim_h3
                     if (phi <= limit_hexa3) then
                        call mof3d_derivatives_penta_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! Above φlim_h3
                        call mof3d_derivatives_hexa(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if ! φlim_h3
                  end if ! φlim_pb1
               else ! φlim_h1
                  ! φlim_pl1
                  limit_penta_left1 = PI_2 - atan((l23 - t2)*cos_theta)

                  ! Below φlim_pl1
                  ! Note: the inequality must be strict here
                  if (phi < limit_penta_left1) then
                     call mof3d_derivatives_penta_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pl1
                     call mof3d_derivatives_quad_face_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if ! φlim_pl1
               end if ! φlim_h1
            else ! θ1 < θ ≤ θ3
               ! φlim_h1
               sqrthexa = sqrt(l23*t2)
               limit_hexa1 = PI_2 - atan((1d0 - 2d0*sqrthexa*cos((acos((3d0*(1d0 - l23 - t2) + t2*t2) &
                  &        / (2d0*l23*sqrthexa)) + 4d0*PI)/3d0))*cos_theta)

               ! Below φlim_h1
               if (phi < limit_hexa1) then
                  ! φlim_pb1
                  limit_penta_below1 = PI_2 - atan((sin_theta + cos_theta)/l12)

                  ! Below φlim_pb1
                  if (phi <= limit_penta_below1) then
                     call mof3d_derivatives_quad_face_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pb1
                     ! φlim_h3
                     sqrthexa = sqrt((2d0 - l23)*t2)
                     limit_hexa3 = PI_2 - atan((1d0 + t2 + 2d0*sqrthexa*cos((acos(3d0*(1d0 - l23)*(1d0 + t2) &
                        &        / (2d0*(2d0 - l23)*sqrthexa)) + 4d0*PI)/3d0))*cos_theta)

                     ! Below φlim_h3
                     if (phi <= limit_hexa3) then
                        call mof3d_derivatives_penta_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! Above φlim_h3
                        call mof3d_derivatives_hexa(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if ! φlim_h3
                  end if ! φlim_pb1
               else ! Above φlim_h1
                  ! φlim_pl2
                  limit_penta_left2 = PI_2 - atan((t2 - (-t2 + sqrt(12d0*l23*t2 - 3d0*t2**2))/2d0)*cos_theta)

                  ! Below φlim_pl2
                  ! Note: the inequality must be strict here
                  if (phi < limit_penta_left2) then
                     call mof3d_derivatives_penta_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pl2
                     call mof3d_derivatives_quad_edge_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if ! φlim_pl2
               end if ! φlim_h1
            end if ! Left of θ3
         else ! Right of θ3
            ! θ ≤ θ5
            if (theta <= atan(1d0/l13)) then
               ! φlim_h2
               sqrthexa = sqrt(l13*t1)
               limit_hexa2 = PI_2 - atan((1d0 - 2d0*sqrthexa*cos((acos((3d0*(1d0 - l13 - t1) + t1*t1) &
                  &        / (2d0*l13*sqrthexa)) + 4d0*PI)/3d0))*sin_theta)

               ! Below φlim_h2
               if (phi < limit_hexa2) then
                  ! φlim_pb1
                  limit_penta_below1 = PI_2 - atan((sin_theta + cos_theta)/l12)

                  ! Below φlim_pb1
                  if (phi <= limit_penta_below1) then
                     call mof3d_derivatives_quad_face_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! φlim_pb1
                     ! φlim_h3
                     sqrthexa = sqrt((2d0 - l23)*t2)
                     limit_hexa3 = PI_2 - atan((1d0 + t2 + 2d0*sqrthexa*cos((acos(3d0*(1d0 - l23)*(1d0 + t2) &
                        &        / (2d0*(2d0 - l23)*sqrthexa)) + 4d0*PI)/3d0))*cos_theta)

                     ! Below φlim_h3
                     if (phi <= limit_hexa3) then
                        call mof3d_derivatives_penta_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! Above φlim_h3
                        call mof3d_derivatives_hexa(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if ! φlim_h3
                  end if ! φlim_pb1
               else ! Above φlim_h2
                  ! φlim_pr2
                  limit_penta_right2 = PI_2 - atan((t1 - (-t1 + sqrt(12d0*l13*t1 - 3d0*t1**2))/2d0)*sin_theta)

                  ! Below φlim_pr2
                  ! Note: the inequality must be strict here
                  if (phi < limit_penta_right2) then
                     call mof3d_derivatives_penta_right(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pr2
                     call mof3d_derivatives_quad_edge_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if ! φlim_pr2
               end if ! φlim_h2
            else if (theta < PI_2 - atan((3d0 - sqrt(12d0*l13 - 3d0))/2d0)) then ! θ5 < θ ≤ θ4h
               ! φlim_h2
               sqrthexa = sqrt(l13*t1)
               limit_hexa2 = PI_2 - atan((1d0 - 2d0*sqrthexa*cos((acos((3d0*(1d0 - l13 - t1) + t1*t1) &
                  &        / (2d0*l13*sqrthexa)) + 4d0*PI)/3d0))*sin_theta)

               ! Below φlim_h2
               if (phi < limit_hexa2) then
                  ! φlim_pb1
                  limit_penta_below1 = PI_2 - atan((sin_theta + cos_theta)/l12)

                  ! Below φlim_pb1
                  if (phi <= limit_penta_below1) then
                     call mof3d_derivatives_quad_face_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pb1
                     ! φlim_h3
                     sqrthexa = sqrt((2d0 - l23)*t2)
                     limit_hexa3 = PI_2 - atan((1d0 + t2 + 2d0*sqrthexa*cos((acos(3d0*(1d0 - l23)*(1d0 + t2) &
                        &        / (2d0*(2d0 - l23)*sqrthexa)) + 4d0*PI)/3d0))*cos_theta)

                     ! Below φlim_h3
                     if (phi <= limit_hexa3) then
                        call mof3d_derivatives_penta_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! Above φlim_h3
                        call mof3d_derivatives_hexa(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if ! φlim_h3
                  end if ! φlim_pb1
               else ! Above φlim_h2
                  ! φlim_pr1
                  limit_penta_right1 = PI_2 - atan((l13 - t1)*sin_theta)

                  ! Below φlim_pr1
                  ! Note: the inequality must be strict here
                  if (phi < limit_penta_right1) then
                     call mof3d_derivatives_penta_right(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pr1
                     call mof3d_derivatives_quad_face_right(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if ! φlim_pr1
               end if ! φlim_h2
            else ! θ4h < θ ≤ π/2
               ! φlim_pr3
               limit_penta_right3 = PI_2 - atan((l13 + t1 + sqrt((l13 + t1)**2 - 4d0/3d0*t1**2))/2d0*sin_theta)

               ! Below φlim_pr3
               if (phi <= limit_penta_right3) then
                  ! φlim_pb1
                  limit_penta_below1 = PI_2 - atan((sin_theta + cos_theta)/l12)

                  ! Below φlim_pb1
                  if (phi <= limit_penta_below1) then
                     call mof3d_derivatives_quad_face_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pb1
                     ! φlim_pb2
                     limit_penta_below2 = PI_2 - atan((sin_theta - cos_theta &
                        &               + cos_theta**2/(3d0*sin_theta))/l12)

                     ! Below φlim_pb2
                     if (phi < limit_penta_below2) then
                        call mof3d_derivatives_penta_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     else ! Above φlim_pb2
                        call mof3d_derivatives_quad_edge_right(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                     end if ! φlim_pb2
                  end if ! φlim_pb1
               else ! Above φlim_pr3
                  ! φlim_pr1
                  limit_penta_right1 = PI_2 - atan((l13 - t1)*sin_theta)

                  ! Below φlim_pr1
                  ! Note: the inequality must be strict here
                  if (phi < limit_penta_right1) then
                     call mof3d_derivatives_penta_right(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  else ! Above φlim_pr1
                     call mof3d_derivatives_quad_face_right(trigo, volume, derivative(:,1), derivative(:,2), centroid)
                  end if ! φlim_pr1
               end if ! φlim_pr3
            end if ! Right of θ3
         end if ! θ ≤ θ3 ?
      end if ! θ1 ≤ θ2h and θ4h ≤ θ5 ?
   end subroutine mof3d_compute_analytic_derivatives_above_one_sixth

   pure subroutine mof3d_compute_analytic_derivatives_half(angles, volume, centroid, derivative)
      double precision, dimension(2), intent(in) :: angles
      double precision, intent(in) :: volume
      double precision, dimension(3), intent(out) :: centroid
      double precision, dimension(3,2), intent(out) :: derivative

      double precision, dimension(4) :: trigo
      double precision :: theta, phi, theta3
      double precision :: tan_theta, cos_theta, sin_theta
      double precision :: limit_penta_left1, limit_penta_right1, limit_penta_below1

      ! Aliases for θ and φ
      theta = angles(1)
      phi = angles(2)

      ! Compute trigonometric functions
      trigo(C_THETA) = cos(theta)
      trigo(S_THETA) = sin(theta)
      trigo(C_PHI)   = cos(phi)
      trigo(S_PHI)   = sin(phi)

      ! Aliases of trigonometric functions
      cos_theta = trigo(C_THETA)
      sin_theta = trigo(S_THETA)
      tan_theta = sin_theta/cos_theta

      ! θ3
      theta3 = PI_4

      ! φlim_pb1
      limit_penta_below1 = PI_2 - atan(cos_theta + sin_theta)

      ! Left of θ3
      if (theta <= theta3) then
         ! Below φlim_pb1
         if (phi <= limit_penta_below1) then
            call mof3d_derivatives_quad_face_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
         else ! Above φlim_pb1
            ! φlim_pl1
            limit_penta_left1 = PI_2 - atan(cos_theta - sin_theta)

            ! Below φlim_pl1
            ! Note: the inequality must be strict here
            if (phi < limit_penta_left1) then
               call mof3d_derivatives_hexa(trigo, volume, derivative(:,1), derivative(:,2), centroid)
            else ! Above φlim_pl1
               call mof3d_derivatives_quad_face_left(trigo, volume, derivative(:,1), derivative(:,2), centroid)
            end if ! φlim_pl1
         end if ! φlim_pb1
      else ! Right of θ3
         ! Below φlim_pb1
         if (phi <= limit_penta_below1) then
            call mof3d_derivatives_quad_face_below(trigo, volume, derivative(:,1), derivative(:,2), centroid)
         else ! Above φlim_pb1
            ! φlim_pr1
            limit_penta_right1 = PI_2 - atan(sin_theta - cos_theta)

            ! Below φlim_pr1
            ! Note: the inequality must be strict here
            if (phi < limit_penta_right1) then
               call mof3d_derivatives_hexa(trigo, volume, derivative(:,1), derivative(:,2), centroid)
            else ! Above φlim_pr1
               call mof3d_derivatives_quad_face_right(trigo, volume, derivative(:,1), derivative(:,2), centroid)
            end if ! φlim_pr1
         end if ! φlim_pb1
      end if ! θ ≤ θ3 ?
   end subroutine mof3d_compute_analytic_derivatives_half

   pure subroutine mof3d_derivatives_triangle(trigo, volume, derivative_theta, derivative_phi, centroid)
      double precision, dimension(4), intent(in) :: trigo
      double precision, intent(in) :: volume
      double precision, dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      double precision :: coef, t2l, t3l, l23, tan_theta, sec_theta, csc_phi, cot_phi, dthetat2l, dthetat3l, dphit3l

      tan_theta = trigo(S_THETA)/trigo(C_THETA)
      sec_theta = 1d0/trigo(C_THETA)
      csc_phi   = 1d0/trigo(S_PHI)
      cot_phi   = trigo(C_PHI)/trigo(S_PHI)

      t2l = tan_theta
      t3l = cot_phi*sec_theta
      l23 = 2d0*volume

      coef = (3d0*l23*t2l*t3l)**(1d0/3d0)

      ! Centroid
      centroid = coef/4d0*[1d0, 1d0/t2l, 1d0/t3l]

      ! d/dθ
      dthetat2l = (1d0 + t2l**2)
      dthetat3l = t2l*t3l

      derivative_theta = l23/4d0/(coef**2)*[t3l*dthetat2l + t2l*dthetat3l    , &
         &                                  dthetat3l - 2d0*t3l/t2l*dthetat2l, &
         &                                  dthetat2l - 2d0*t2l/t3l*dthetat3l]

      ! d/dφ
      dphit3l = -csc_phi**2*sec_theta

      derivative_phi = l23/4d0/(coef**2)*[ &
         &    t2l*dphit3l        , &
         &    dphit3l            , &
         &    -2d0*t2l/t3l*dphit3l &
         & ]
   end subroutine mof3d_derivatives_triangle

   pure subroutine mof3d_derivatives_quad_face_left(trigo, volume, derivative_theta, derivative_phi, centroid)
      double precision, dimension(4), intent(in) :: trigo
      double precision, intent(in) :: volume
      double precision, dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      double precision :: t2l, t3l, l23, tan_theta, sec_theta, csc_phi, cot_phi, dthetat2l, dthetat3l, dphit3l

      tan_theta = trigo(S_THETA)/trigo(C_THETA)
      sec_theta = 1d0/trigo(C_THETA)
      csc_phi   = 1d0/trigo(S_PHI)
      cot_phi   = trigo(C_PHI)/trigo(S_PHI)

      t2l = tan_theta
      t3l = cot_phi*sec_theta

      l23 = 2d0*volume

      ! Centroid
      centroid = 1d0/(12d0*l23)*[           &
         &    3d0*l23**2 + t2l**2 + t3l**2, &
         &    2d0*(3d0*l23 - t2l)         , &
         &    2d0*(3d0*l23 - t3l)           &
         & ]

      ! d/dθ
      dthetat2l = 1d0 + tan_theta**2
      dthetat3l = tan_theta*t3l

      derivative_theta = 1d0/(6d0*l23)*[t2l*dthetat2l + t3l*dthetat3l, -dthetat2l, -dthetat3l]

      ! d/dφ
      dphit3l = -csc_phi**2*sec_theta

      derivative_phi = 1d0/(6d0*l23)*[t3l*dphit3l, 0d0, -dphit3l]
   end subroutine mof3d_derivatives_quad_face_left

   pure subroutine mof3d_derivatives_quad_face_right(trigo, volume, derivative_theta, derivative_phi, centroid)
      double precision, dimension(4), intent(in) :: trigo
      double precision, intent(in) :: volume
      double precision, dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      double precision :: t2r, t3r, l13, csc_theta, cot_theta, csc_phi, cot_phi, dthetat2r, dthetat3r, dphit2r

      csc_theta = 1d0/trigo(S_THETA)
      cot_theta = trigo(C_THETA)/trigo(S_THETA)

      csc_phi   = 1d0/trigo(S_PHI)
      cot_phi   = trigo(C_PHI)/trigo(S_PHI)

      t2r = csc_theta*cot_phi
      t3r = cot_theta

      l13 = 2d0*volume

      ! Centroid
      centroid = 1d0/(12d0*l13)*[           &
         &    2d0*(3d0*l13 - t3r)         , &
         &    3d0*l13**2 + t2r**2 + t3r**2, &
         &    2d0*(3d0*l13 - t2r)           &
         & ]

      ! d/dθ
      dthetat2r = -cot_theta*csc_theta*cot_phi
      dthetat3r = -csc_theta**2

      derivative_theta = 1d0/(6d0*l13)*[-dthetat3r, t2r*dthetat2r + t3r*dthetat3r, -dthetat2r]

      ! d/dφ
      dphit2r = -csc_theta*csc_phi**2

      derivative_phi = 1d0/(6d0*l13)*[0d0, t2r*dphit2r, -dphit2r]
   end subroutine mof3d_derivatives_quad_face_right

   pure subroutine mof3d_derivatives_quad_face_below(trigo, volume, derivative_theta, derivative_phi, centroid)
      double precision, dimension(4), intent(in) :: trigo
      double precision, intent(in) :: volume
      double precision, dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      double precision :: t2b, t3b, l12, cos_theta, sin_theta, sec_phi, tan_phi, dthetat2b, dthetat3b, dphit2b, dphit3b

      cos_theta = trigo(C_THETA)
      sin_theta = trigo(S_THETA)

      sec_phi   = 1d0/trigo(C_PHI)
      tan_phi   = trigo(S_PHI)/trigo(C_PHI)

      t2b = cos_theta*tan_phi
      t3b = sin_theta*tan_phi

      l12 = 2d0*volume

      ! Centroid
      centroid = 1d0/(12d0*l12)*[2d0*(3d0*l12 - t2b), 2d0*(3d0*l12 - t3b), 3d0*l12**2 + t2b**2 + t3b**2]

      ! d/dθ
      dthetat2b = -t3b
      dthetat3b = t2b

      derivative_theta = 1d0/(6d0*l12)*[-dthetat2b, -dthetat3b, t2b*dthetat2b + t3b*dthetat3b]

      ! d/dφ
      dphit2b = cos_theta*sec_phi**2
      dphit3b = sin_theta*sec_phi**2

      derivative_phi = 1d0/(6d0*l12)*[-dphit2b, -dphit3b, t2b*dphit2b + t3b*dphit3b]
   end subroutine mof3d_derivatives_quad_face_below

   pure subroutine mof3d_derivatives_quad_edge_left(trigo, volume, derivative_theta, derivative_phi, centroid)
      double precision, dimension(4), intent(in) :: trigo
      double precision, intent(in) :: volume
      double precision, dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      double precision :: cot_phi, csc_phi, tan_theta, sec_theta
      double precision :: t2l, t3l, fqe, xqel, dtt2l, dtt3l, dtxqel, dft3l, dfxqel, l23

      tan_theta = trigo(S_THETA)/trigo(C_THETA)
      sec_theta = 1d0/trigo(C_THETA)
      csc_phi   = 1d0/trigo(S_PHI)
      cot_phi   = trigo(C_PHI)/trigo(S_PHI)

      l23 = 2d0*volume
      t2l = tan_theta
      t3l = cot_phi*sec_theta
      fqe = t3l**2 + 6d0*l23*t2l
      xqel = sqrt(36d0*l23*t2l - 3d0*t3l**2)

      ! Centroid
      centroid = 1d0/(108d0*l23*t2l)* [      &
         &    xqel*fqe                     , &
         &    xqel*fqe/t2l                 , &
         &    3d0*(18d0*l23*t2l - t3l*xqel)  &
         & ]

      ! d/dθ
      dtt2l = 1d0 + t2l**2
      dtt3l = t2l*t3l
      dtxqel = (18d0*l23*dtt2l - 3d0*t3l*dtt3l)/xqel

      derivative_theta = 1d0/(108d0*l23*t2l**2)                                               &
         &             * [(2d0*t2l*dtt3l - t3l*dtt2l)*t3l*xqel + fqe*t2l*dtxqel             , &
         &                (2d0*(t3l*dtt3l - (t3l**2/t2l + 3d0*l23)*dtt2l)*xqel + fqe*dtxqel), &
         &                3d0*((t3l*dtt2l - t2l*dtt3l)*xqel - t2l*t3l*dtxqel)                 &
         &               ]

      ! d/dφ
      dft3l = -csc_phi**2*sec_theta
      dfxqel = -3d0*t3l*dft3l/xqel

      derivative_phi = 1d0/(108d0*l23*t2l**2)                        &
         &           * [(2d0*t2l*dft3l)*t3l*xqel + fqe*t2l*dfxqel  , &
         &              2d0*t3l*dft3l*xqel + fqe*dfxqel            , &
         &              3d0*(-t2l*dft3l*xqel - t2l*t3l*dfxqel)]
   end subroutine mof3d_derivatives_quad_edge_left

   pure subroutine mof3d_derivatives_quad_edge_right(trigo, volume, derivative_theta, derivative_phi, centroid)
      double precision, dimension(4), intent(in) :: trigo
      double precision, intent(in) :: volume
      double precision, dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      double precision :: csc_phi, cot_phi, cot_theta, csc_theta
      double precision :: t2r, t3r, fqe, xqer, dtt2r, dtt3r, dtxqer, dft2r, dfxqer, l13

      csc_theta = 1d0/trigo(S_THETA)
      cot_theta = trigo(C_THETA)/trigo(S_THETA)
      csc_phi   = 1d0/trigo(S_PHI)
      cot_phi   = trigo(C_PHI)/trigo(S_PHI)

      t2r = csc_theta*cot_phi
      t3r = cot_theta

      l13 = 2d0*volume

      fqe = t3r**2 + 6d0*l13*t2r
      xqer = sqrt(36d0*l13*t2r - 3d0*t3r**2)

      ! Centroid
      centroid = 1d0/(108d0*l13*t2r)             &
         &     * [3d0*(18d0*l13*t2r - t3r*xqer), &
         &        xqer*fqe                     , &
         &        xqer*fqe/t2r                   &
         &       ]

      ! d/dθ
      dtt2r = -cot_theta*csc_theta*cot_phi
      dtt3r = -csc_theta**2
      dtxqer = (18d0*l13*dtt2r - 3d0*t3r*dtt3r)/xqer

      derivative_theta = 1d0/(108d0*l13*t2r**2)                                             &
         &             * [3d0*((t3r*dtt2r - t2r*dtt3r)*xqer - t2r*t3r*dtxqer)             , &
         &                (2d0*t2r*dtt3r - t3r*dtt2r)*t3r*xqer + fqe*t2r*dtxqer           , &
         &                2d0*(t3r*dtt3r - (t3r**2/t2r + 3d0*l13)*dtt2r)*xqer + fqe*dtxqer]

      ! d/dφ
      dft2r = -csc_theta*csc_phi**2
      dfxqer = 18d0*l13*dft2r/xqer

      derivative_phi = 1d0/(108d0*l13*t2r**2)                                &
         &           * [3d0*(t3r*dft2r*xqer - t2r*t3r*dfxqer)              , &
         &              -t3r*dft2r*t3r*xqer + fqe*t2r*dfxqer               , &
         &              -2d0*(t3r**2/t2r + 3d0*l13)*dft2r*xqer + fqe*dfxqer]
   end subroutine mof3d_derivatives_quad_edge_right

   pure subroutine mof3d_derivatives_quad_edge_below(trigo, volume, derivative_theta, derivative_phi, centroid)
      double precision, dimension(4), intent(in) :: trigo
      double precision, intent(in) :: volume
      double precision, dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      double precision :: sec_phi, tan_phi, cos_theta, sin_theta
      double precision :: t2b, t3b, fqe, xqeb, dtt2b, dtt3b, dtxqeb, dft2b, dft3b, dfxqeb, l12

      cos_theta = trigo(C_THETA)
      sin_theta = trigo(S_THETA)
      sec_phi   = 1d0/trigo(C_PHI)
      tan_phi   = trigo(S_PHI)/trigo(C_PHI)

      l12 = 2d0*volume

      t2b = cos_theta*tan_phi
      t3b = sin_theta*tan_phi
      fqe = t3b**2 + 6d0*l12*t2b
      xqeb = sqrt(36d0*l12*t2b - 3d0*t3b**2)

      ! Centroid
      centroid = 1d0/(108d0*l12*t2b)             &
         &     * [xqeb*fqe/t2b                 , &
         &        3d0*(18d0*l12*t2b - t3b*xqeb), &
         &        xqeb*fqe                     ]

      ! d/dθ
      dtt2b = -t3b
      dtt3b = t2b
      dtxqeb = (18d0*l12*dtt2b - 3d0*t3b*dtt3b)/xqeb

      derivative_theta = 1d0/(108d0*l12*t2b**2)                                             &
         &             * [2d0*(t3b*dtt3b - (t3b**2/t2b + 3d0*l12)*dtt2b)*xqeb + fqe*dtxqeb, &
         &                3d0*((t3b*dtt2b - t2b*dtt3b)*xqeb - t2b*t3b*dtxqeb)             , &
         &                (2d0*t2b*dtt3b - t3b*dtt2b)*t3b*xqeb + fqe*t2b*dtxqeb           ]

      ! d/dφ
      dft2b = cos_theta*sec_phi**2
      dft3b = sin_theta*sec_phi**2
      dfxqeb = (18d0*l12*dft2b - 3d0*t3b*dft3b)/xqeb

      derivative_phi = 1d0/(108d0*l12*t2b**2)                                              &
         &           * [2d0*(t3b*dft3b - (t3b**2/t2b + 3d0*l12)*dft2b)*xqeb + fqe*dfxqeb , &
         &              3d0*((t3b*dft2b - t2b*dft3b)*xqeb - t2b*t3b*dfxqeb)              , &
         &              (2d0*t2b*dft3b - t3b*dft2b)*t3b*xqeb + fqe*t2b*dfxqeb            ]
   end subroutine mof3d_derivatives_quad_edge_below

   ! Chen & Zhang formulas
   pure subroutine mof3d_derivatives_with_points(trigo, point, volume, derivative_theta, derivative_phi)
      use mod_cg3_points
      double precision, dimension(4), intent(in) :: trigo
      double precision, dimension(:,:), intent(in) :: point
      double precision, intent(in) :: volume
      double precision, dimension(3), intent(out) :: derivative_theta, derivative_phi

      double precision, dimension(3) :: i_phi, j_phi, k_phi
      double precision, dimension(3) :: origin, pref, translated_point
      double precision, dimension(2) :: d_theta, d_phi
      double precision :: cos_theta, sin_theta, cos_phi, sin_phi
      double precision :: integral_xy, integral_xx, integral_yy
      double precision :: area, triangle_area, x1, x2, x3, y1, y2, y3
      integer :: i
      double precision, dimension(3,size(point,2)) :: p

      ! Compute the origin of the coordinate system (the centroid of the polygon)
      area = 0d0
      origin = 0d0
      pref = point(:,1)

      do i = 2, size(point,2) - 1
         triangle_area = norm2(cg3_cross_product(point(:,i) - point(:,1), point(:,i+1) - point(:,1)))
         area = area + triangle_area
         origin = origin + (point(:,i) + point(:,i+1) - 2d0*point(:,1))*triangle_area
      end do

      area = area/2d0

      if (area > 0d0) then
         origin = pref + origin/(6d0*area)
      end if

      cos_theta = trigo(C_THETA)
      sin_theta = trigo(S_THETA)
      cos_phi   = trigo(C_PHI)
      sin_phi   = trigo(S_PHI)

      i_phi = [-sin_theta, cos_theta, 0d0]
      j_phi = [cos_theta*cos_phi, sin_theta*cos_phi, -sin_phi]
      k_phi = [cos_theta*sin_phi, sin_theta*sin_phi, cos_phi]

      ! Project all the points
      do i = 1, size(point,2)
         translated_point = point(:,i) - origin
         p(1,i) = dot_product(translated_point, i_phi)
         p(2,i) = dot_product(translated_point, j_phi)
         p(3,i) = dot_product(translated_point, k_phi)
      end do

      ! Compute integrals
      integral_xx = 0d0
      integral_xy = 0d0
      integral_yy = 0d0

      x3 = (p(1,1) + p(1,2))/2d0
      y3 = (p(2,1) + p(2,2))/2d0
      do i = 2, size(point,2) - 1
         x1 = x3
         x2 = (p(1,i  ) + p(1,i+1))/2d0
         x3 = (p(1,i+1) + p(1,1  ))/2d0
         y1 = y3
         y2 = (p(2,i  ) + p(2,i+1))/2d0
         y3 = (p(2,i+1) + p(2,1  ))/2d0
         triangle_area = norm2(cg3_cross_product(p(:,i) - p(:,1), p(:,i+1) - p(:,1)))/2d0
         integral_xx = integral_xx + triangle_area*(x1**2 + x2**2 + x3**2)/3d0
         integral_xy = integral_xy + triangle_area*(x1*y1 + x2*y2 + x3*y3)/3d0
         integral_yy = integral_yy + triangle_area*(y1**2 + y2**2 + y3**2)/3d0
      end do

      ! Compute the partial derivatives
      d_theta = -sin_phi/volume*[integral_xx, integral_xy]
      d_phi   = -1d0/volume*[integral_xy, integral_yy]

      derivative_theta = d_theta(1)*i_phi + d_theta(2)*j_phi
      derivative_phi = d_phi(1)*i_phi + d_phi(2)*j_phi
   end subroutine mof3d_derivatives_with_points

   pure subroutine mof3d_derivatives_penta_left(trigo, volume, derivative_theta, derivative_phi, centroid)
      double precision, dimension(4), intent(in) :: trigo
      double precision, intent(in) :: volume
      double precision, dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      double precision :: t2, t3, sqrt2t2t3, xp, tan_theta, sec_theta, csc_phi, cot_phi
      double precision :: l23, alp, bet, gam, del, eps
      double precision, dimension(3,5) :: point

      tan_theta = trigo(S_THETA)/trigo(C_THETA)
      sec_theta = 1d0/trigo(C_THETA)
      csc_phi   = 1d0/trigo(S_PHI)
      cot_phi   = trigo(C_PHI)/trigo(S_PHI)

      t2 = tan_theta
      t3 = cot_phi*sec_theta
      sqrt2t2t3 = sqrt(2d0*t2*t3)

      l23 = 2d0*volume
      xp = cos((acos(3d0*(t2 + t3 - l23)/(4d0*sqrt2t2t3)) + 4d0*PI)/3d0)

      ! Centroid
      centroid = 1d0/(6d0*l23)*[f0(t2,t3) + xp*(f1(t2,t3) + 24d0*t2*t3*xp), &
         &                      f2(t2,t3) + xp*(f3(t2,t3) - 24d0*t3*xp)   , &
         &                      f2(t3,t2) + xp*(f3(t3,t2) - 24d0*t2*xp)   ]

      ! Partial derivatives
      bet = t3 + 2d0*sqrt2t2t3*xp
      alp = bet + t2
      gam = alp - t3
      del = gam/t2
      eps = bet/t3

      ! Points of the top polygon
      point(:,1) = [alp, 0d0, 0d0]
      point(:,2) = [bet, 1d0, 0d0]
      point(:,3) = [0d0, 1d0, eps]
      point(:,4) = [0d0, del, 1d0]
      point(:,5) = [gam, 0d0, 1d0]

      call mof3d_derivatives_with_points(trigo, point, volume, derivative_theta, derivative_phi)

   contains

      double precision pure function f0(x, y)
         double precision, intent(in) :: x, y
         f0 = 2d0*(x + y)**2 - x*y
      end function f0

      double precision pure function f1(x, y)
         double precision, intent(in) :: x, y
         f1 = 3d0*sqrt2t2t3*(3d0*(x + y) + l23)
      end function f1

      double precision pure function f2(x, y)
         double precision, intent(in) :: x, y
         f2 = 6d0*l23 - 4d0*x - 3d0*y
      end function f2

      double precision pure function f3(x, y)
         double precision, intent(in) :: x, y
         f3 = 6d0*y/sqrt2t2t3*(l23 - 5d0*x - y)
      end function f3

   end subroutine mof3d_derivatives_penta_left

   pure subroutine mof3d_derivatives_penta_right(trigo, volume, derivative_theta, derivative_phi, centroid)
      double precision, dimension(4), intent(in) :: trigo
      double precision, intent(in) :: volume
      double precision, dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      double precision :: p1, p2, t2, t3, sqrt2p1p2, xp, csc_theta, cot_theta, csc_phi, cot_phi
      double precision :: l13, alp, bet, gam, del, eps
      double precision, dimension(3,5) :: point

      csc_theta = 1d0/trigo(S_THETA)
      cot_theta = trigo(C_THETA)/trigo(S_THETA)

      csc_phi   = 1d0/trigo(S_PHI)
      cot_phi   = trigo(C_PHI)/trigo(S_PHI)

      p1 = csc_theta*cot_phi
      p2 = cot_theta
      t2 = p1
      t3 = p2

      sqrt2p1p2 = sqrt(2d0*p1*p2)

      l13 = 2d0*volume
      xp = cos((acos(3d0*(p1 + p2 - l13)/(4d0*sqrt2p1p2)) + 4d0*PI)/3d0)

      ! Centroid
      centroid = 1d0/(6d0*l13)*[f2(p2,p1) + xp*(f3(p2,p1) - 24d0*p1*xp)   , &
         &                      f0(p1,p2) + xp*(f1(p1,p2) + 24d0*p1*p2*xp), &
         &                      f2(p1,p2) + xp*(f3(p1,p2) - 24d0*p2*xp)   ]

      ! Partial derivatives
      bet = t3 + 2d0*sqrt2p1p2*xp
      alp = bet + t2
      gam = alp - t3
      del = gam/t2
      eps = bet/t3

      ! Points of the top polygon
      point(:,1) = [0d0, alp, 0d0]
      point(:,2) = [0d0, bet, 1d0]
      point(:,3) = [eps, 0d0, 1d0]
      point(:,4) = [1d0, 0d0, del]
      point(:,5) = [1d0, gam, 0d0]

      call mof3d_derivatives_with_points(trigo, point, volume, derivative_theta, derivative_phi)

   contains

      double precision pure function f0(x, y)
         double precision, intent(in) :: x, y
         f0 = 2d0*(x + y)**2 - x*y
      end function f0

      double precision pure function f1(x, y)
         double precision, intent(in) :: x, y
         f1 = 3d0*sqrt2p1p2*(3d0*(x + y) + l13)
      end function f1

      double precision pure function f2(x, y)
         double precision, intent(in) :: x, y
         f2 = 6d0*l13 - 4d0*x - 3d0*y
      end function f2

      double precision pure function f3(x, y)
         double precision, intent(in) :: x, y
         f3 = 6d0*y/sqrt2p1p2*(l13 - 5d0*x - y)
      end function f3

   end subroutine mof3d_derivatives_penta_right

   pure subroutine mof3d_derivatives_penta_below(trigo, volume, derivative_theta, derivative_phi, centroid)
      double precision, dimension(4), intent(in) :: trigo
      double precision, intent(in) :: volume
      double precision, dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      double precision :: t2, t3, p1, p2, sqrt2p1p2, xp, cos_theta, sin_theta, tan_phi, sec_phi
      double precision :: l12, alp, bet, gam, del, eps
      double precision, dimension(3,5) :: point

      cos_theta = trigo(C_THETA)
      sin_theta = trigo(S_THETA)

      sec_phi   = 1d0/trigo(C_PHI)
      tan_phi   = trigo(S_PHI)/trigo(C_PHI)

      p1 = cos_theta*tan_phi
      p2 = sin_theta*tan_phi
      t2 = p1
      t3 = p2

      sqrt2p1p2 = sqrt(2d0*p1*p2)

      l12 = 2d0*volume
      xp = cos((acos(3d0*(p1 + p2 - l12)/(4d0*sqrt2p1p2)) + 4d0*PI)/3d0)

      ! Centroid
      centroid = 1d0/(6d0*l12)*[f2(p1,p2) + xp*(f3(p1,p2) - 24d0*p2*xp)   , &
         &                      f2(p2,p1) + xp*(f3(p2,p1) - 24d0*p1*xp)   , &
         &                      f0(p1,p2) + xp*(f1(p1,p2) + 24d0*p1*p2*xp)]

      ! Partial derivatives
      bet = t3 + 2d0*sqrt2p1p2*xp
      alp = bet + t2
      gam = alp - t3
      del = gam/t2
      eps = bet/t3

      ! Points of the top polygon
      point(:,1) = [0d0, 0d0, alp]
      point(:,2) = [1d0, 0d0, bet]
      point(:,3) = [1d0, eps, 0d0]
      point(:,4) = [del, 1d0, 0d0]
      point(:,5) = [0d0, 1d0, gam]

      call mof3d_derivatives_with_points(trigo, point, volume, derivative_theta, derivative_phi)

   contains

      double precision pure function f0(x, y)
         double precision, intent(in) :: x, y
         f0 = 2d0*(x + y)**2 - x*y
      end function f0

      double precision pure function f1(x, y)
         double precision, intent(in) :: x, y
         f1 = 3d0*sqrt2p1p2*(3d0*(x + y) + l12)
      end function f1

      double precision pure function f2(x, y)
         double precision, intent(in) :: x, y
         f2 = 6d0*l12 - 4d0*x - 3d0*y
      end function f2

      double precision pure function f3(x, y)
         double precision, intent(in) :: x, y
         f3 = 6d0*y/sqrt2p1p2*(l12 - 5d0*x - y)
      end function f3

   end subroutine mof3d_derivatives_penta_below

   pure subroutine mof3d_derivatives_hexa(trigo, volume, derivative_theta, derivative_phi, centroid)
      double precision, dimension(4), intent(in) :: trigo
      double precision, intent(in) :: volume
      double precision, dimension(3), intent(out) :: derivative_theta, derivative_phi, centroid

      double precision :: t2, t3, t4, sqrtt4, l23, xh, tan_theta, sec_theta, csc_phi, cot_phi
      double precision :: alp, bet, gam, del, eps, lam
      double precision, dimension(3,6) :: point

      tan_theta = trigo(S_THETA)/trigo(C_THETA)
      sec_theta = 1d0/trigo(C_THETA)
      csc_phi   = 1d0/trigo(S_PHI)
      cot_phi   = trigo(C_PHI)/trigo(S_PHI)

      t2  = tan_theta
      t3  = cot_phi*sec_theta
      t4  = 4d0*t3 - (1d0 - t2 + t3)**2
      l23 = 2d0*volume
      xh  = cos((acos(6d0*t2*t3*(1d0 - l23)/t4**1.5d0) + 4d0*PI)/3d0)
      sqrtt4 = sqrt(t4)

      centroid(1) = f0(1d0, t2, t3)/t2 - 24d0*(1d0 - l23)*(sqrtt4*xh + 2d0) &
         &        + 12d0*xh*xh*f1(1d0, t2, t3, t4)/(t2*t3)
      centroid(2) = f0(t2, 1d0, t3)/(t2*t2) - 24d0*(1d0 - l23)*(sqrtt4*xh/t2 + 2d0) &
         &        + 12d0*xh*xh*f1(t2, 1d0, t3, t4)/(t2*t2*t3)
      centroid(3) = f2(1d0, t2, t3)/t2 - 24d0*(1d0 - l23)*(sqrtt4*xh/t3 + 2d0) &
         &        + 12d0*xh*xh*f1(t3, 1d0, t2, t4)/(t2*t3*t3)

      centroid = centroid/(96d0*l23)

      bet = (1d0 - t2 + t3)/2d0 + sqrt(t4)*xh
      gam = bet + t2 - t3
      alp = (bet - 1d0 + t2)/t3
      del = gam/t2
      eps = bet/t3
      lam = 1d0/t2*(bet + t2 - 1d0)

      ! Points of the top polygon
      point(:,1) = [1d0, 0d0, alp]
      point(:,2) = [1d0, lam, 0d0]
      point(:,3) = [bet, 1d0, 0d0]
      point(:,4) = [0d0, 1d0, eps]
      point(:,5) = [0d0, del, 1d0]
      point(:,6) = [gam, 0d0, 1d0]

      call mof3d_derivatives_with_points(trigo, point, volume, derivative_theta, derivative_phi)

   contains

      double precision pure function f0(x, y, z)
         double precision, intent(in) :: x, y, z

         f0 = (x - y)**3*(3d0*x + y)/z - (8d0*x**3 + z**3 - 4d0*y*(y*y + z*z + 9d0*x*x) + 6d0*z*(y*y - x*x))
      end function f0

      double precision pure function f1(x, y, z, t)
         double precision, intent(in) :: x, y, z, t

         f1 = (2d0*(x*x - (y - z)**2) - t)*t
      end function f1

      double precision pure function f2(x, y, z)
         double precision, intent(in) :: x, y, z

         f2 = -((x - y)**2/z)**2 + 6d0*((x + y)**2 + 4d0*x*y) - 8d0*z*(x + y) + 3d0*z*z
      end function f2

   end subroutine mof3d_derivatives_hexa

   ! map x in [0, 2π[
   double precision pure function modulo_tau(x) result(r)
      double precision, intent(in) :: x

      double precision, parameter :: TAU = 4d0*acos(0d0)

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

end module mod_mof3d_analytic_centroid
