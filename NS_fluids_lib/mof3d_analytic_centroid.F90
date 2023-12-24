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

module mod_mof3d_analytic_centroid
use amrex_fort_module, only : amrex_real
   implicit none

   private

   real(amrex_real), parameter :: PI_2 = acos(0d0)
   real(amrex_real), parameter :: PI = 2d0*acos(0d0)

   public :: mof3d_compute_analytic_gradient, mof3d_compute_analytic_gradient_symmetric

   ! Routines hierarchy
   ! ------------------
   !
   ! mof3d_compute_analytic_gradient                  → Compute gradient from input angles [public]
   ! mof3d_compute_analytic_gradient_symmetric        → Compute symmetric gradient from input angles [public]
   ! ├── mof3d_transform_to_reference_map             → Transform input angles to reference map
   ! │   └── mof3d_is_point_inside_map                → Check if the couple of angles belongs to the reference map
   ! ├── mof3d_compute_analytic_derivatives_reference → Compute the partial derivatives on the reference map
   ! │   ├── mof3d_derivatives_quad_face              → Compute partial derivatives on the QuadFace sub-region of the reference map
   ! │   ├── mof3d_derivatives_triangle               → Compute partial derivatives on the Triangle sub-region of the reference map
   ! │   ├── mof3d_derivatives_quad_edge              → Compute partial derivatives on the QuadEdge sub-region of the reference map
   ! │   ├── mof3d_derivatives_penta                  → Compute partial derivatives on the Penta sub-region of the reference map
   ! │   └── mof3d_derivatives_hexa                   → Compute partial derivatives on the Hexa sub-region of the reference map
   ! └── mof3d_transform_gradient_to_global           → Transform the gradient from the reference map to the global map

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
#if (DEBUG_ANALYTICAL_GRAD==0)
   pure &
#endif
      subroutine mof3d_compute_analytic_gradient(angles,ref_centroid, &
        volume, c, centroid, gradient)
      real(amrex_real), dimension(2), INTENT(in) :: angles
      real(amrex_real), dimension(3), INTENT(in) :: ref_centroid
      real(amrex_real), INTENT(in) :: volume
      real(amrex_real), dimension(3), INTENT(in) :: c
      real(amrex_real), dimension(3), INTENT(out) :: centroid
      real(amrex_real), dimension(2), INTENT(out) :: gradient

      real(amrex_real), dimension(3,2) :: t_partial_derivative
      real(amrex_real), dimension(3) :: t_c, t_centroid, t_ref_centroid
      real(amrex_real), dimension(2) :: t_angles, t_gradient
      integer, dimension(3) :: permutation, sign_permutation, inverse_permutation, sign_inverse_permutation

#if (DEBUG_ANALYTICAL_GRAD==1)
      print *,"angles ",angles(1),angles(2)
      print *,"ref_centroid ",ref_centroid(1), &
              ref_centroid(2),ref_centroid(3)
      print *,"volume ",volume
      print *,"c ",c(1),c(2),c(3)
#endif 

      ! Compute transformed angle and permutation
      call mof3d_transform_to_reference_map(angles, c, volume, t_angles, permutation, sign_permutation)

      ! Compute the inverse of the permutation
      inverse_permutation(permutation([1,2,3])) = [1,2,3]
      sign_inverse_permutation = sign_permutation(inverse_permutation)

      ! Apply the transformations to the cell dimensions
      t_c = c(permutation)

      ! Compute the gradient and the centroid in the reference map
      call mof3d_compute_analytic_derivatives_reference(t_angles, volume, t_c, t_centroid, t_partial_derivative)

      ! Transform the centroid from the reference map to the global map
      centroid = t_centroid(inverse_permutation)
      where (sign_inverse_permutation == -1) centroid = c - centroid

      ! Apply the transformations to the reference centroid
      t_ref_centroid = ref_centroid(permutation)
      where (sign_permutation == -1) t_ref_centroid = t_c - t_ref_centroid

      ! Compute the gradient in the reference map
      t_gradient = [2d0*dot_product(t_centroid - t_ref_centroid, t_partial_derivative(:,1)), &
         &          2d0*dot_product(t_centroid - t_ref_centroid, t_partial_derivative(:,2))]

      ! Transform the gradient from the reference map to the global map
      call mof3d_transform_gradient_to_global(angles, permutation, sign_permutation, t_gradient, gradient)


#if (DEBUG_ANALYTICAL_GRAD==1)
      print *,"AFTER"
      print *,"angles ",angles(1),angles(2)
      print *,"ref_centroid ",ref_centroid(1), &
              ref_centroid(2),ref_centroid(3)
      print *,"volume ",volume
      print *,"c ",c(1),c(2),c(3)
      print *,"centroid ",centroid(1), &
              centroid(2),centroid(3)
      print *,"gradient ",gradient(1),gradient(2)
#endif

   end subroutine mof3d_compute_analytic_gradient

   !> Compute the centroid and the gradient of the objective function in rectangular hexahedral cell using
   !! the symmetric objective function.
   !!
   !! @param[in]  angles:        Spherical angles.
   !! @param[in]  ref_centroid1: Coordinates of the reference centroid of the material 1.
   !! @param[in]  ref_centroid2: Coordinates of the reference centroid of the material 2.
   !! @param[in]  ref_volume:    Reference volume of the material 1.
   !! @param[in]  c:             Dimensions of the cell.
   !! @param[out] centroids:     Coordinates of the centroid.
   !! @param[out] gradient:      Gradient of the objective function.
   !! @ingroup moment_of_fluid
   pure subroutine mof3d_compute_analytic_gradient_symmetric(angles, ref_centroid1, ref_centroid2, volume, c, centroid, gradient)
      real(amrex_real), dimension(2), INTENT(in) :: angles
      real(amrex_real), dimension(3), INTENT(in) :: ref_centroid1, ref_centroid2
      real(amrex_real), INTENT(in) :: volume
      real(amrex_real), dimension(3), INTENT(in) :: c
      real(amrex_real), dimension(3), INTENT(out) :: centroid
      real(amrex_real), dimension(2), INTENT(out) :: gradient

      real(amrex_real), dimension(3,2) :: t_partial_derivative
      real(amrex_real), dimension(3) :: t_c, t_centroid, t_ref_centroid, h_ref_centroid
      real(amrex_real), dimension(2) :: t_angles, t_gradient, gradient1, gradient2
      integer, dimension(3) :: permutation, sign_permutation, inverse_permutation, sign_inverse_permutation
      real(amrex_real) :: cell_volume

      ! Compute the volume of the cell
      cell_volume = c(1)*c(2)*c(3)

      ! Compute transformed angle and permutation
      call mof3d_transform_to_reference_map(angles, c, volume, t_angles, permutation, sign_permutation)

      ! Compute the inverse of the permutation
      inverse_permutation(permutation([1,2,3])) = [1,2,3]
      sign_inverse_permutation = sign_permutation(inverse_permutation)

      ! Apply the transformations to the cell dimensions
      t_c = c(permutation)

      ! Compute the gradient and the centroid in the reference map
      call mof3d_compute_analytic_derivatives_reference(t_angles, volume, t_c, t_centroid, t_partial_derivative)

      ! Transform the centroid from the reference map to the global map
      centroid = t_centroid(inverse_permutation)
      where (sign_inverse_permutation == -1) centroid = c - centroid

      ! Contribution of the first reference centroid

      ! Apply the transformations to the reference centroid
      t_ref_centroid = ref_centroid1(permutation)
      where (sign_permutation == -1) t_ref_centroid = t_c - t_ref_centroid

      ! Compute the gradient in the reference map
      t_gradient = [2d0*dot_product(t_centroid - t_ref_centroid, t_partial_derivative(:,1)), &
         &          2d0*dot_product(t_centroid - t_ref_centroid, t_partial_derivative(:,2))]

      ! Transform the gradient from the reference map to the global map
      call mof3d_transform_gradient_to_global(angles, permutation, sign_permutation, t_gradient, gradient1)

      ! Contribution of the second reference centroid

      ! Apply the transformations to the reference centroid
      h_ref_centroid = (cell_volume*[c(1)/2d0, c(2)/2d0, c(3)/2d0] - (cell_volume - volume)*ref_centroid2)/volume
      t_ref_centroid = h_ref_centroid(permutation)
      where (sign_permutation == -1) t_ref_centroid = t_c - t_ref_centroid

      ! Compute the gradient in the reference map
      t_gradient = [2d0*dot_product(t_centroid - t_ref_centroid, t_partial_derivative(:,1)), &
         &          2d0*dot_product(t_centroid - t_ref_centroid, t_partial_derivative(:,2))]

      ! Transform the gradient from the reference map to the global map
      call mof3d_transform_gradient_to_global(angles, permutation, sign_permutation, t_gradient, gradient2)

      ! Compute the gradient
      gradient = gradient1 + 2d0*(volume/(cell_volume - volume))**2*gradient2
   end subroutine mof3d_compute_analytic_gradient_symmetric

   pure subroutine mof3d_transform_gradient_to_global(angles, permutation, sign_permutation, t_gradient, gradient)
      real(amrex_real), dimension(2), INTENT(in) :: angles
      integer, dimension(3), INTENT(in) :: permutation
      integer, dimension(3), INTENT(in) :: sign_permutation
      real(amrex_real), dimension(2), INTENT(in) :: t_gradient
      real(amrex_real), dimension(2), INTENT(out) :: gradient

      real(amrex_real) :: sin_theta, cos_theta, tan_theta, sin_phi, cos_phi, tan_phi
      real(amrex_real) :: qsqrt
      real(amrex_real), dimension(2,2) :: matrix

      ! Transform the gradient back to the original domain
      select case (permutation(1))
      case(1)
         if (permutation(2) == 2) then
            ! σ = (1, 2, 3)
            if (modulo_tau(angles(2)) >= PI) then
               gradient = [sign_permutation(1)*sign_permutation(2), -sign_permutation(3)]*t_gradient
            else
               gradient = [sign_permutation(1)*sign_permutation(2),  sign_permutation(3)]*t_gradient
            end if
         else
            ! σ = (1, 3, 2)
            cos_theta = cos(angles(1))
            sin_theta = sin(angles(1))
            tan_theta = sin_theta/cos_theta
            cos_phi   = cos(angles(2))
            sin_phi   = sin(angles(2))
            tan_phi   = sin_phi/cos_phi

            qsqrt = 1d0/sqrt(1d0 - (sin_theta*sin_phi)**2)

            matrix(1,1) =  sign_permutation(1)*sign_permutation(2)*sin_theta*sin_phi*cos_phi/((cos_theta*sin_phi)**2 + cos_phi**2)
            matrix(2,1) = -sign_permutation(1)*sign_permutation(2)*cos_theta/((cos_theta*sin_phi)**2 + cos_phi**2)
            matrix(1,2) = -sign_permutation(3)*cos_theta*sin_phi*qsqrt
            matrix(2,2) = -sign_permutation(3)*sin_theta*cos_phi*qsqrt

            gradient = matmul(matrix, t_gradient)
         end if
      case(2)
         if (permutation(2) == 1) then
            ! σ = (2, 1, 3)
            if (modulo_tau(angles(2)) >= PI) then
               gradient = [-sign_permutation(1)*sign_permutation(2), -sign_permutation(3)]*t_gradient
            else
               gradient = [-sign_permutation(1)*sign_permutation(2),  sign_permutation(3)]*t_gradient
            end if
         else
            ! σ = (2, 3, 1)
            cos_theta = cos(angles(1))
            sin_theta = sin(angles(1))
            tan_theta = sin_theta/cos_theta
            cos_phi   = cos(angles(2))
            sin_phi   = sin(angles(2))
            tan_phi   = sin_phi/cos_phi

            qsqrt = 1d0/sqrt(1d0 - (cos_theta*sin_phi)**2)

            matrix(1,1) = -sign_permutation(1)*sign_permutation(2)*cos_theta*sin_phi*cos_phi/((sin_theta*sin_phi)**2 + cos_phi**2)
            matrix(2,1) = -sign_permutation(1)*sign_permutation(2)*sin_theta/((sin_theta*sin_phi)**2 + cos_phi**2)
            matrix(1,2) =  sign_permutation(3)*sin_theta*sin_phi*qsqrt
            matrix(2,2) = -sign_permutation(3)*cos_theta*cos_phi*qsqrt

            gradient = matmul(matrix, t_gradient)
         end if
      case default ! 3
         cos_theta = cos(angles(1))
         sin_theta = sin(angles(1))
         cos_phi   = cos(angles(2))
         sin_phi   = sin(angles(2))
         tan_phi   = sin_phi/cos_phi

         if (permutation(2) == 1) then
            ! σ = (3, 1, 2)
            qsqrt = 1d0/sqrt(1d0 - (sin_theta*sin_phi)**2)

            matrix(1,1) = -sign_permutation(1)*sign_permutation(2)*sin_theta*tan_phi/(1d0 + (cos_theta*tan_phi)**2)
            matrix(2,1) =  sign_permutation(1)*sign_permutation(2)*cos_theta/(cos_phi**2 + (cos_theta*sin_phi)**2)
            matrix(1,2) = -sign_permutation(3)*cos_theta*sin_phi*qsqrt
            matrix(2,2) = -sign_permutation(3)*sin_theta*cos_phi*qsqrt

            gradient = matmul(matrix, t_gradient)
         else
            ! σ = (3, 2, 1)
            qsqrt = 1d0/sqrt(1d0 - (cos_theta*sin_phi)**2)

            matrix(1,1) =  sign_permutation(1)*sign_permutation(2)*cos_theta*sin_phi*cos_phi/((sin_theta*sin_phi)**2 + cos_phi**2)
            matrix(2,1) =  sign_permutation(1)*sign_permutation(2)*sin_theta/((sin_theta*sin_phi)**2 + cos_phi**2)
            matrix(1,2) =  sign_permutation(3)*sin_theta*sin_phi*qsqrt
            matrix(2,2) = -sign_permutation(3)*cos_theta*cos_phi*qsqrt

            gradient = matmul(matrix, t_gradient)
         end if
      end select
   end subroutine mof3d_transform_gradient_to_global

   pure subroutine mof3d_transform_to_reference_map(angles, c, volume, t_angles, permutation, sign_permutation)
      real(amrex_real), dimension(2), INTENT(in) :: angles
      real(amrex_real), dimension(3), INTENT(in) :: c
      real(amrex_real), INTENT(in) :: volume
      real(amrex_real), dimension(2), INTENT(out) :: t_angles
      integer, dimension(3), INTENT(out) :: permutation, sign_permutation

      real(amrex_real), dimension(3) :: t_c
      real(amrex_real) :: theta_limit, phi_limit, sin_theta, cos_theta, cot_theta
      real(amrex_real) :: q, l13, t1, sqrtp3
      logical :: rotate_yzx

      ! Initialize to identity
      permutation = [1, 2, 3]
      sign_permutation = 1

      ! Compute permutation
      t_angles = angles

      ! Ensure that the angles are in their domains of definition:
      ! → θ ∈ ]-π,π]
      ! → φ ∈ [0,π]
      !-----------------------------------------------------------

      ! Check if φ belongs to [0,π]
      ! We start by φ since it transforms θ too
      if (t_angles(2) < 0d0 .or. t_angles(2) >= PI) then
         ! First, ensure that φ belongs to [0,2π[
         t_angles(2) = modulo_tau(t_angles(2))

         ! Ιf φ belongs to [π,2π[, transform to ]0,π]
         ! Otherwise, φ belongs to [0,π[
         if (t_angles(2) >= PI) then
            t_angles(2) = 2d0*PI - t_angles(2)
            t_angles(1) = t_angles(1) + PI
         end if
      end if

      ! Check if θ belongs to [-π,π]
      if (t_angles(1) < -PI .or. t_angles(1) > PI) then
         t_angles(1) = modulo_tau(t_angles(1) + PI) - PI
      end if

      ! Check if θ belongs to ]-π,π] (-π excluded)
      if (t_angles(1) == -PI) t_angles(1) = PI

      ! Crop to [0,π/2[ × [0,π/2]
      !--------------------------

      ! Here θ belongs to ]-π,π] and φ belongs to [0,π]

      ! Check if φ belongs to [0,π/2]
      if (t_angles(2) >= PI_2) then
         ! Transform [π/2,π] → [0,π/2]
         ! Reflection on third axis
         t_angles(2) = PI - t_angles(2)
         sign_permutation(3) = -sign_permutation(3)
      end if

      ! Warning: special case when φ = 0
      ! Note: we are sure that φ belongs to ]0,π/2] after this point
      if (t_angles(2) == 0d0) then
         ! Rotate by π/2 around the second axis
         t_angles(1) = 0d0
         t_angles(2) = PI_2
         sign_permutation = [sign_permutation(3), sign_permutation(2), -sign_permutation(1)]
         permutation = permutation([3,2,1])
         ! Return here
         return
      end if

      ! Check if θ belongs to [0,π/2[
      if (t_angles(1) == PI) then
         ! Rotate by π around third axis
         t_angles(1) = 0d0
         sign_permutation(1) = -1
         sign_permutation(2) = -1
      else if (t_angles(1) >= PI_2) then
         ! Transform [π/2,π[ → [0,π/2[
         ! Rotate by -π/2 around third axis
         t_angles(1) = t_angles(1) - PI_2
         permutation(1) = 2
         permutation(2) = 1
         sign_permutation(2) = -1
      else if (t_angles(1) < -PI_2) then
         ! Transform ]-π,-π/2[ → ]0,π/2[
         ! Rotate by π around third axis
         t_angles(1) = t_angles(1) + PI
         sign_permutation(1) = -1
         sign_permutation(2) = -1
      else if (t_angles(1) < 0d0) then
         ! Transform [-π/2,0[ → [0,π/2[
         ! Rotate by π/2 around third axis
         t_angles(1) = t_angles(1) + PI_2
         permutation(1) = 2
         permutation(2) = 1
         sign_permutation(1) = -1
      end if

      ! Here θ belongs to [0,π/2[ and φ belongs to ]0,π/2]

      ! Permute the dimensions of the cell
      t_c = c(permutation)

      ! Check if the angles belongs to the reference map
      !-------------------------------------------------

      if (mof3d_is_point_inside_map(t_angles, volume, t_c)) return

      ! Rotate around the diagonal of the cuboid
      !-----------------------------------------
      ! We have to decide which rotation to apply to be in the reference map
      !
      !           Z                          Z
      !           ×                          ×
      !           |                          |
      !       ↗   |   ↘                  ↙   |   ↖
      !           ×                          ×
      !        .-' `-.                    .-' `-.
      !      ×'   ←   `×                ×'   →   `×
      !    X             Y            X             Y

      ! Determinate in which direction to turn around the diagonal (c1,c2,c3) of the cuboid
      ! We suppose that the rotation (or permutation) is (2,3,1) (right picture)
      rotate_yzx = .false.

      ! Compute θ3
      theta_limit = atan(t_c(1)/t_c(2))

      ! Pre-compute trigonometric functions
      sin_theta = sin(t_angles(1))
      cos_theta = cos(t_angles(1))
      cot_theta = cos_theta/sin_theta

      ! If θ < θ3, the rotation is zxy. Otherwise, we need more information to decide which rotation to apply
      if (t_angles(1) >= theta_limit) then
         l13 = 2d0*volume/(t_c(1)*t_c(3))
         t1 = t_c(1)*cot_theta

         ! If the volume fraction < 1/6, the limit is between the right QuadEdge and the bottom Penta.
         ! If the volume fraction > 1/6, there are two limits to consider.
         if (3d0*l13 <= t_c(2)) then
            ! Limit between the right quad edge and the bottom Penta (φ_rot)
            phi_limit = PI_2 - atan((3d0*t_c(2)**2 - 3d0*t_c(2)*t1 + t1*t1)*sin_theta/(3d0*t_c(3)*l13))
         else
            ! Θ4^h simplified with arctan(1/x) = π/2 - arctan(x)
            theta_limit = PI_2 - atan((3d0*t_c(2) - sqrt(12d0*l13*t_c(2) - 3d0*t_c(2)**2))/(2d0*t_c(1)))

            ! If θ ≥ θ4, the limit to consider is the limit between the right QuadEdge and the bottom Penta
            if (t_angles(1) >= theta_limit) then
               ! Limit between the right QuadEdge and the bottom Penta (φ_rot)
               phi_limit = PI_2 - atan((3d0*t_c(2)**2 - 3d0*t_c(2)*t1 + t1*t1)*sin_theta/(3d0*t_c(3)*l13))
            else
               ! Limit between the right Penta and the Hexa (φlim_h2)
               sqrtp3 = sqrt(l13*t1)
               q = (3d0*t_c(2)*(t_c(2) - l13 - t1) + t1*t1)/l13
               phi_limit = PI_2 - atan((t_c(2) - 2d0*sqrtp3*cos((acos(q/(2d0*sqrtp3)) + 4d0*PI)/3d0))*sin_theta/t_c(3))
            end if
         end if

         if (t_angles(2) > phi_limit) rotate_yzx = .true.
      end if

      if (rotate_yzx) then
         !           Z
         !           ×
         !           |
         !       ↗   |   ↘          σ = (3,1,2)
         !           ×            σ⁻¹ = (2,3,1)
         !        .-' `-.
         !      ×'   ←   `×
         !    X             Y
         !
         t_angles = [PI_2 - atan(sin_theta*tan(t_angles(2))), acos(cos_theta*sin(t_angles(2)))]
         sign_permutation = sign_permutation([2,3,1])
         permutation = permutation([2,3,1])
      else
         !           Z
         !           ×
         !           |
         !       ↙   |   ↖          σ = (2,3,1)
         !           ×            σ⁻¹ = (3,1,2)
         !        .-' `-.
         !      ×'   →   `×
         !    X             Y
         !
         t_angles = [atan(cos_theta*tan(t_angles(2))), acos(sin_theta*sin(t_angles(2)))]
         t_angles(2) = min(t_angles(2), PI_2) ! Avoid round-off problems
         sign_permutation = sign_permutation([3,1,2])
         permutation = permutation([3,1,2])
      end if
   end subroutine mof3d_transform_to_reference_map

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

   ! Hyp: θ in -π, π
   ! Hyp: φ in 0, π
   ! Hyp: volume fraction <= 0.5
   logical pure function mof3d_is_point_inside_map(angles, volume, c) result(is_inside)
      real(amrex_real), dimension(2), INTENT(in) :: angles
      real(amrex_real), INTENT(in) :: volume
      real(amrex_real), dimension(3), INTENT(in) :: c

      real(amrex_real) :: l13, l23, t1, t2
      real(amrex_real) :: l_theta2, l_theta3, l_theta4, l_theta5
      real(amrex_real) :: limit, limit_b, sin_theta, cos_theta, tan_theta
      real(amrex_real) :: q, sqrtp3, c123

      is_inside = .false.

      if (angles(1) < 0d0) return
      if (angles(2) > PI_2) return

      c123 = c(1)*c(2)*c(3)

      tan_theta = tan(angles(1))

      l23 = 2d0*volume/(c(2)*c(3))
      l13 = 2d0*volume/(c(1)*c(3))
      t1 = c(1)/tan_theta
      t2 = c(2)*tan_theta

      if (6d0*volume < c123) then
         l_theta5 = atan(c(1)/l13)

         if (angles(1) > l_theta5) return

         l_theta2 = atan(3d0*l23/c(2))
         l_theta3 = atan(c(1)/c(2))
         l_theta4 = atan(c(1)/(3d0*l13))

         if (angles(1) <= l_theta2) then
            ! Quad-face, QuadEdge, Penta
            ! Limit Penta 3 (φlim_p3)
            limit = PI_2 - atan((l23 + t2 + sqrt((l23 + t2)**2 - 4d0/3d0*t2*t2))*cos(angles(1))/(2d0*c(3)))

            if (angles(2) < limit) return
         else if (angles(1) <= l_theta3) then
            ! Quad-edge, Triangle
            ! Limit Triangle 1 (φlim_t1)
            limit = PI_2 - atan(c(2)*t2*sin(angles(1))/(3d0*l23*c(3)))

            if (angles(2) < limit) return
         else if (angles(1) <= l_theta4) then
            ! Quad-edge, Triangle
            ! Limit Triangle 2 (φlim_t2)
            limit = PI_2 - atan(c(1)*t1*cos(angles(1))/(3d0*l13*c(3)))

            if (angles(2) < limit) return
         else
            ! Limit QuadEdge (φlim_qe)
            limit = PI_2 - atan((3d0*t1 - sqrt(12d0*l13*t1 - 3d0*t1**2))*sin(angles(1))/(2d0*c(3)))

            if (angles(2) < limit) return
         end if
      else if (volume/c123 <= 0.5d0 - 1d1*epsilon(1d0)) then
         l_theta4 = PI_2 - atan((3d0*c(2) - sqrt(12d0*c(2)*l13 - 3d0*c(2)*c(2)))/(2d0*c(1)))
         l_theta5 = atan(c(1)/l13)

         if (angles(1) > max(l_theta4, l_theta5)) return

         l_theta2 = atan((3d0*c(1) - sqrt(12d0*c(1)*l23 - 3d0*c(1)*c(1)))/(2d0*c(2)))
         cos_theta = cos(angles(1))

         if (angles(1) <= l_theta2) then
            ! Limit Penta 3 (φlim_p3)
            limit = PI_2 - atan((l23 + t2 + sqrt((l23 + t2)**2 - 4d0/3d0*t2*t2))*cos_theta/(2d0*c(3)))

            if (angles(2) < limit) return
         else if (angles(1) <= l_theta4) then
            ! Limit Hexa 3 (φlim_h3)
            sqrtp3 = sqrt((2d0*c(1) - l23)*t2)
            q = (3d0*(c(1) - l23)*(c(1) + t2))/(2d0*c(1) - l23)
            limit = PI_2 - atan((c(1) + t2 + 2d0*sqrtp3*cos((acos(q/(2d0*sqrtp3)) + 4d0*PI)/3d0))*cos_theta/c(3))

            if (angles(2) < limit) return
         end if

         l_theta3 = atan(c(1)/c(2))

         if (angles(1) > l_theta3) then ! θ > 3 and θ < 5
            sin_theta = sin(angles(1))

            ! Limit QuadEdge (φlim_qe)
            limit = PI_2 - atan((3d0*t1 - sqrt(12d0*l13*t1 - 3d0*t1**2))*sin_theta/(2d0*c(3)))

            ! Limit of Hexa 2 (φlim_h2)
            sqrtp3 = sqrt(l13*t1)
            q = (3d0*c(2)*(c(2) - l13 - t1) + t1*t1)/l13
            limit_b = PI_2 - atan((c(2) - 2d0*sqrtp3*cos((acos(q/(2d0*sqrtp3)) + 4d0*PI)/3d0))*sin_theta/c(3))

            if (angles(1) <= l_theta4) then
               if (angles(1) <= l_theta5) then
                  ! θ < 4 and Θ < 5
                  if (angles(2) < limit .and. angles(2) > limit_b) return
               else
                  ! θ < 4 and θ > 5
                  if (angles(2) > limit_b) return
               end if
            else
               ! θ > 4 and θ < 5
               if (angles(2) < limit) return
            end if
         end if
      else ! Volume fraction = 1/2
         l_theta4 = PI_2

         if (angles(1) > l_theta4) return

         cos_theta = cos(angles(1))

         ! Limit of Hexa 3 (φlim_h3)
         if (angles(1) < PI_2) then
            limit = atan(c(3)/(cos_theta*(c(1) + t2)))
         else
            limit = atan(c(3)/c(2))
         end if

         if (angles(2) < limit) return

         l_theta3 = atan(c(1)/c(2))

         if (angles(1) > l_theta3) then
            ! Limit of Hexa 2 (φlim_h2)
            limit_b = atan(c(3)/(c(2)*sin(angles(1)) - c(1)*cos_theta))

            if (angles(2) > limit_b) return
         end if
      end if

      is_inside = .true.
   end function mof3d_is_point_inside_map

   pure subroutine mof3d_compute_analytic_derivatives_reference(angles, volume, c, centroid, partial_derivative)
      real(amrex_real), dimension(2), INTENT(in) :: angles
      real(amrex_real), INTENT(in) :: volume
      real(amrex_real), dimension(3), INTENT(in) :: c
      real(amrex_real), dimension(3), INTENT(out) :: centroid
      real(amrex_real), dimension(3,2), INTENT(out) :: partial_derivative

      real(amrex_real) :: l13, l23, t1, t2, q, sqrtp3, c123
      real(amrex_real) :: l_theta1, l_theta2, l_theta3, l_theta4
      real(amrex_real) :: limit, tan_theta, cos_theta

      c123 = c(1)*c(2)*c(3)

      tan_theta = tan(angles(1))

      l23 = 2d0*volume/(c(2)*c(3))
      l13 = 2d0*volume/(c(1)*c(3))
      t1 = c(1)/tan_theta
      t2 = c(2)*tan_theta

      if (6d0*volume < c123) then
         l_theta1 = atan(l23/c(2))
         l_theta2 = atan(3d0*l23/c(2))
         l_theta3 = atan(c(1)/c(2))
         l_theta4 = atan(c(1)/(3d0*l13))

         if (angles(1) <= l_theta1) then
            ! Limit Penta 1 (φlim_p1)
            limit = PI_2 - atan((l23 - t2)*cos(angles(1))/c(3))

            if (angles(2) >= limit) then
               call mof3d_derivatives_quad_face(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
            else
               call mof3d_derivatives_penta(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
            end if
         else if (angles(1) <= l_theta2) then
            ! Limit Penta 2 (φlim_p2)
            limit = PI_2 - atan((3d0*t2 - sqrt(12d0*l23*t2 - 3d0*t2*t2))*cos(angles(1))/(2d0*c(3)))

            if (angles(2) >= limit) then
               call mof3d_derivatives_quad_edge(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
            else
               call mof3d_derivatives_penta(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
            end if
         else if (angles(1) <= l_theta4) then
            ! Limit Triangle 3 (φlim_t3)
            limit = PI_2 - atan(sqrt(3d0*l23*t2)*cos(angles(1))/c(3))

            if (angles(2) >= limit) then
               call mof3d_derivatives_quad_edge(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
            else
               call mof3d_derivatives_triangle(angles, volume, partial_derivative(:,1), partial_derivative(:,2), centroid)
            end if
         else
            call mof3d_derivatives_quad_edge(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
         end if
      else if (volume/c123 <= 0.5d0 - 1d1*epsilon(1d0)) then
         l_theta1 = atan(l23/c(2))
         l_theta2 = atan((3d0*c(1) - sqrt(12d0*c(1)*l23 - 3d0*c(1)*c(1)))/(2d0*c(2)))

         if (angles(1) <= max(l_theta1, l_theta2)) then
            if (l_theta1 <= l_theta2) then
               if (angles(1) <= l_theta1) then
                  ! Limit Penta 1 (φlim_p1)
                  limit = PI_2 - atan((l23 - t2)*cos(angles(1))/c(3))

                  if (angles(2) >= limit) then
                     call mof3d_derivatives_quad_face(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
                  else
                     call mof3d_derivatives_penta(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
                  end if
               else if (angles(1) <= l_theta2) then
                  ! Limit Penta 2 (φlim_p2)
                  limit = PI_2 - atan((3d0*t2 - sqrt(12d0*l23*t2 - 3d0*t2*t2))*cos(angles(1))/(2d0*c(3)))

                  if (angles(2) >= limit) then
                     call mof3d_derivatives_quad_edge(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
                  else
                     call mof3d_derivatives_penta(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
                  end if
               end if
            else ! l_theta1 > l_theta2
               cos_theta = cos(angles(1))

               ! Limit Penta 1 (φlim_p1)
               limit = PI_2 - atan((l23 - t2)*cos_theta/c(3))

               if (angles(1) <= l_theta2) then
                  if (angles(2) >= limit) then
                     call mof3d_derivatives_quad_face(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
                  else
                     call mof3d_derivatives_penta(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
                  end if
               else if (angles(2) >= limit) then
                  call mof3d_derivatives_quad_face(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
               else
                  ! Limit Hexa 1 (φlim_h1)
                  sqrtp3 = sqrt(l23*t2)
                  q = (3d0*c(1)*(c(1) - l23 - t2) + t2*t2)/l23
                  limit = PI_2 - atan((c(1) - 2d0*sqrtp3*cos((acos(q/(2d0*sqrtp3)) + 4d0*PI)/3d0))*cos_theta/c(3))

                  if (angles(2) >= limit) then
                     call mof3d_derivatives_penta(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
                  else
                     call mof3d_derivatives_hexa(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
                  end if
               end if
            end if
         else ! θ > l_theta1 and θ > l_theta2
            l_theta3 = atan(c(1)/c(2))

            if (angles(1) <= l_theta3) then
               cos_theta = cos(angles(1))

               ! Limit Penta 2 (φlim_p2)
               limit = PI_2 - atan((3d0*t2 - sqrt(12d0*l23*t2 - 3d0*t2*t2))*cos_theta/(2d0*c(3)))

               if (angles(2) >= limit) then
                  call mof3d_derivatives_quad_edge(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
               else
                  ! Limit Hexa 1 (φlim_h1)
                  sqrtp3 = sqrt(l23*t2)
                  q = (3d0*c(1)*(c(1) - l23 - t2) + t2*t2)/l23
                  limit = PI_2 - atan((c(1) - 2d0*sqrtp3*cos((acos(q/(2d0*sqrtp3)) + 4d0*PI)/3d0))*cos_theta/c(3))

                  if (angles(2) >= limit) then
                     call mof3d_derivatives_penta(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
                  else
                     call mof3d_derivatives_hexa(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
                  end if
               end if
            else ! θ > l_theta3
               ! Limit QuadEdge (φlim_qe)
               limit = PI_2 - atan((3d0*t1 - sqrt(12d0*l13*t1 - 3d0*t1**2))*sin(angles(1))/(2d0*c(3)))

               if (angles(2) >= limit) then
                  call mof3d_derivatives_quad_edge(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
               else
                  call mof3d_derivatives_hexa(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
               end if
            end if
         end if
      else ! volume = c123/2
         l_theta3 = atan(c(1)/c(2))

         if (angles(1) <= l_theta3) then
            ! Limit Penta 1 (φlim_p1)
            limit = PI_2 - atan((l23 - t2)*cos(angles(1))/c(3))

            if (angles(2) >= limit) then
               call mof3d_derivatives_quad_face(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
            else
               call mof3d_derivatives_hexa(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
            end if
         else
            call mof3d_derivatives_hexa(angles, volume, c, partial_derivative(:,1), partial_derivative(:,2), centroid)
         end if
      end if
   end subroutine mof3d_compute_analytic_derivatives_reference

   pure subroutine mof3d_derivatives_triangle(angles, volume, derivative_theta, derivative_phi, centroid)
      real(amrex_real), dimension(2), INTENT(in) :: angles
      real(amrex_real), INTENT(in) :: volume
      real(amrex_real), dimension(3), INTENT(out) :: derivative_theta, derivative_phi, centroid

      real(amrex_real) :: coef, t1, t2, tan_theta, sec_theta, csc_phi, cot_phi, dthetat1, dthetat2, dphit1, dphit2

      tan_theta = tan(angles(1))
      sec_theta = 1d0/cos(angles(1))
      csc_phi   = 1d0/sin(angles(2))
      cot_phi   = 1d0/tan(angles(2))

      t1 = tan_theta
      t2 = cot_phi*sec_theta

      coef = (6d0*volume*t1*t2)**(1d0/3d0)

      ! Centroid
      centroid = 0.25d0*coef*[1d0, 1d0/t1, 1d0/t2]

      ! d/dθ
      dthetat1 = 1d0 + t1*t1
      dthetat2 = t1*t2

      derivative_theta = 0.5d0*volume/(coef**2)*[t2*dthetat1 + t1*dthetat2         , &
         &                                       (t1*dthetat2 - 2d0*t2*dthetat1)/t1, &
         &                                       (t2*dthetat1 - 2d0*t1*dthetat2)/t2]

      ! d/dφ
      dphit1 = 0d0
      dphit2 = -csc_phi**2*sec_theta

      derivative_phi = 0.5d0*volume/(coef**2)*[t2*dphit1 + t1*dphit2         , &
         &                                     (t1*dphit2 - 2d0*t2*dphit1)/t1, &
         &                                     (t2*dphit1 - 2d0*t1*dphit2)/t2]
   end subroutine mof3d_derivatives_triangle

   pure subroutine mof3d_derivatives_quad_face(angles, volume, c, derivative_theta, derivative_phi, centroid)
      real(amrex_real), dimension(2), INTENT(in) :: angles
      real(amrex_real), INTENT(in) :: volume
      real(amrex_real), dimension(3), INTENT(in) :: c
      real(amrex_real), dimension(3), INTENT(out) :: derivative_theta, derivative_phi, centroid

      real(amrex_real) :: qf1, qf2, cv, tan_theta, sec_theta, csc_phi, cot_phi, dthetaqf1, dthetaqf2, dphiqf1, dphiqf2

      tan_theta = tan(angles(1))
      sec_theta = 1d0/cos(angles(1))
      csc_phi   = 1d0/sin(angles(2))
      cot_phi   = 1d0/tan(angles(2))

      qf1 = c(2)*tan_theta
      qf2 = c(3)*cot_phi*sec_theta

      cv = 2d0*volume/(c(2)*c(3))

      ! Centroid
      centroid = 1d0/(12d0*cv)*[3d0*cv**2 + qf1**2 + qf2**2, 2d0*c(2)*(3d0*cv - qf1), 2d0*c(3)*(3d0*cv - qf2)]

      ! d/dθ
      dthetaqf1 = c(2)*(1d0 + tan_theta**2)
      dthetaqf2 = tan_theta*qf2

      derivative_theta = 1d0/(6d0*cv)*[qf1*dthetaqf1 + qf2*dthetaqf2, -c(2)*dthetaqf1, -c(3)*dthetaqf2]

      ! d/dφ
      dphiqf1 = 0d0
      dphiqf2 = -c(3)*csc_phi**2*sec_theta

      derivative_phi = 1d0/(6d0*cv)*[qf1*dphiqf1 + qf2*dphiqf2, -c(2)*dphiqf1, -c(3)*dphiqf2]
   end subroutine mof3d_derivatives_quad_face

   pure subroutine mof3d_derivatives_quad_edge(angles, volume, c, derivative_theta, derivative_phi, centroid)
      real(amrex_real), dimension(2), INTENT(in) :: angles
      real(amrex_real), INTENT(in) :: volume
      real(amrex_real), dimension(3), INTENT(in) :: c
      real(amrex_real), dimension(3), INTENT(out) :: derivative_theta, derivative_phi, centroid

      real(amrex_real) :: cot_phi, csc_phi, tan_theta, sec_theta
      real(amrex_real) :: qe1, qe2, qe3, xqe, dtqe1, dtqe2, dtxqe, dfqe1, dfqe2, dfxqe

      tan_theta = tan(angles(1))
      sec_theta = 1d0/cos(angles(1))
      csc_phi   = 1d0/sin(angles(2))
      cot_phi   = 1d0/tan(angles(2))

      qe1 = tan_theta/c(3)
      qe2 = c(3)*cot_phi*sec_theta
      qe3 = qe2**2 + 12d0*volume*qe1
      xqe = sqrt(72d0*volume*qe1 - 3d0*qe2**2)

      ! Centroid
      centroid = 1d0/(216d0*volume*qe1)                     &
         &     * [xqe*qe3                                 , &
         &        xqe/(c(3)*qe1)*qe3                      , &
         &        3d0*c(3)*(36d0*volume*qe1 - qe2*xqe)]

      ! d/dθ
      dtqe1 = (1d0 + (c(3)*qe1)**2)/c(3)
      dtqe2 = c(3)*qe1*qe2
      dtxqe = (36d0*volume*dtqe1 - 3d0*qe2*dtqe2)/xqe

      derivative_theta = 1d0/(216d0*volume*qe1**2)                                                                  &
         &             * [(2d0*qe1*dtqe2 - qe2*dtqe1)*qe2*xqe + qe3*qe1*dtxqe                                     , &
         &                (2d0*(qe1*qe2*dtqe2 - (qe2**2 + 6d0*volume*qe1)*dtqe1)*xqe + qe3*qe1*dtxqe)/(c(3)*qe1)  , &
         &                3d0*c(3)*((qe2*dtqe1 - qe1*dtqe2)*xqe - qe1*qe2*dtxqe)                                  ]

      ! d/dφ
      dfqe1 = 0d0
      dfqe2 = -c(3)*csc_phi**2*sec_theta
      dfxqe = (36d0*volume*dfqe1 - 3d0*qe2*dfqe2)/xqe

      derivative_phi = 1d0/(216d0*volume*qe1**2)                                                                  &
         &           * [(2d0*qe1*dfqe2 - qe2*dfqe1)*qe2*xqe + qe3*qe1*dfxqe                                     , &
         &              (2d0*(qe1*qe2*dfqe2 - (qe2**2 + 6d0*volume*qe1)*dfqe1)*xqe + qe3*qe1*dfxqe)/(c(3)*qe1)  , &
         &              3d0*c(3)*((qe2*dfqe1 - qe1*dfqe2)*xqe - qe1*qe2*dfxqe)                                  ]
   end subroutine mof3d_derivatives_quad_edge

   pure subroutine mof3d_derivatives_penta(angles, volume, c, derivative_theta, derivative_phi, centroid)
      real(amrex_real), dimension(2), INTENT(in) :: angles
      real(amrex_real), INTENT(in) :: volume
      real(amrex_real), dimension(3), INTENT(in) :: c
      real(amrex_real), dimension(3), INTENT(out) :: derivative_theta, derivative_phi, centroid

      real(amrex_real) :: p1, p2, sqrt2p1p2, xp, cv, tan_theta, sec_theta, csc_phi, cot_phi
      real(amrex_real) :: dtp1, dtp2, dtxp, dtf0, dtf1, dtf21, dtf22, dtf31, dtf32
      real(amrex_real) :: dpp1, dpp2, dpxp, dpf0, dpf1, dpf21, dpf22, dpf31, dpf32

      tan_theta = tan(angles(1))
      sec_theta = 1d0/cos(angles(1))
      csc_phi   = 1d0/sin(angles(2))
      cot_phi   = 1d0/tan(angles(2))

      p1 = c(2)*tan_theta
      p2 = c(3)*cot_phi*sec_theta
      sqrt2p1p2 = sqrt(2d0*p1*p2)

      cv = 2d0*volume/(c(2)*c(3))
      xp = cos((acos(3d0*(p1 + p2 - cv)/(4d0*sqrt2p1p2)) + 4d0*PI)/3d0)

      ! Centroid
      centroid = 1d0/(6d0*cv)*[f0(p1,p2) + xp*(f1(p1,p2) + 24d0*p1*p2*xp)    , &
         &                     c(2)*(f2(p1,p2) + xp*(f3(p1,p2) - 24d0*p2*xp)), &
         &                     c(3)*(f2(p2,p1) + xp*(f3(p2,p1) - 24d0*p1*xp))]

      ! d/dθ
      dtp1  = c(2)*(1d0 + tan_theta**2)
      dtp2  = tan_theta*p2
      dtxp  = (dtp1/p1*(cv + p1 - p2) + dtp2/p2*(cv - p1 + p2))/(8d0*sqrt2p1p2*(4d0*xp**2 - 1d0))
      dtf0  = df0(p1,p2)*dtp1 + df0(p2,p1)*dtp2
      dtf1  = df1(p1,p2)*dtp1 + df1(p2,p1)*dtp2
      dtf21 = -4d0*dtp1 - 3d0*dtp2
      dtf22 = -4d0*dtp2 - 3d0*dtp1
      dtf31 = dxf3(p1,p2)*dtp1 + dyf3(p1,p2)*dtp2
      dtf32 = dxf3(p2,p1)*dtp2 + dyf3(p2,p1)*dtp1

      derivative_theta = 1d0/(6d0*cv)*[dtf0 + dtf1*xp + f1(p1,p2)*dtxp + 24d0*xp*((p1*dtp2 + p2*dtp1)*xp + 2d0*p1*p2*dtxp), &
         &                             c(2)*(dtf21 + dtf31*xp + f3(p1,p2)*dtxp - 24d0*xp*(dtp2*xp + 2d0*p2*dtxp))         , &
         &                             c(3)*(dtf22 + dtf32*xp + f3(p2,p1)*dtxp - 24d0*xp*(dtp1*xp + 2d0*p1*dtxp))         ]

      ! d/dφ
      dpp1  = 0d0
      dpp2  = -c(3)*csc_phi**2*sec_theta
      dpxp  = (dpp1/p1*(cv + p1 - p2) + dpp2/p2*(cv - p1 + p2))/(8d0*sqrt2p1p2*(4d0*xp**2 - 1d0))
      dpf0  = df0(p1,p2)*dpp1 + df0(p2,p1)*dpp2
      dpf1  = df1(p1,p2)*dpp1 + df1(p2,p1)*dpp2
      dpf21 = -4d0*dpp1 - 3d0*dpp2
      dpf22 = -4d0*dpp2 - 3d0*dpp1
      dpf31 = dxf3(p1,p2)*dpp1 + dyf3(p1,p2)*dpp2
      dpf32 = dxf3(p2,p1)*dpp2 + dyf3(p2,p1)*dpp1

      derivative_phi = 1d0/(6d0*cv)*[dpf0 + dpf1*xp + f1(p1,p2)*dpxp + 24d0*xp*((p1*dpp2 + p2*dpp1)*xp + 2d0*p1*p2*dpxp), &
         &                           c(2)*(dpf21 + dpf31*xp + f3(p1,p2)*dpxp - 24d0*xp*(dpp2*xp + 2d0*p2*dpxp))         , &
         &                           c(3)*(dpf22 + dpf32*xp + f3(p2,p1)*dpxp - 24d0*xp*(dpp1*xp + 2d0*p1*dpxp))         ]

   contains

      real(amrex_real) pure function f0(x, y)
         real(amrex_real), INTENT(in) :: x, y
         f0 = 2d0*(x + y)**2 - x*y
      end function f0

      real(amrex_real) pure function f1(x, y)
         real(amrex_real), INTENT(in) :: x, y
         f1 = 3d0*sqrt2p1p2*(3d0*(x + y) + cv)
      end function f1

      real(amrex_real) pure function f2(x, y)
         real(amrex_real), INTENT(in) :: x, y
         f2 = 6d0*cv - 4d0*x - 3d0*y
      end function f2

      real(amrex_real) pure function f3(x, y)
         real(amrex_real), INTENT(in) :: x, y
         f3 = 6d0*y/sqrt2p1p2*(cv - 5d0*x - y)
      end function f3

      real(amrex_real) pure function df0(x, y)
         real(amrex_real), INTENT(in) :: x, y
         df0 = 4d0*x + 3d0*y
      end function df0

      real(amrex_real) pure function df1(x, y)
         real(amrex_real), INTENT(in) :: x, y
         df1 = 3d0*y/sqrt2p1p2*(cv + 9d0*x + 3d0*y)
      end function df1

      real(amrex_real) pure function dxf3(x, y)
         real(amrex_real), INTENT(in) :: x, y
         dxf3 = 3d0*y/(x*sqrt2p1p2)*(y - 5d0*x - cv)
      end function dxf3

      real(amrex_real) pure function dyf3(x, y)
         real(amrex_real), INTENT(in) :: x, y
         dyf3 = 3d0/sqrt2p1p2*(cv - 5d0*x - 3d0*y)
      end function dyf3

   end subroutine mof3d_derivatives_penta

   pure subroutine mof3d_derivatives_hexa(angles, volume, c, derivative_theta, derivative_phi, centroid)
      real(amrex_real), dimension(2), INTENT(in) :: angles
      real(amrex_real), INTENT(in) :: volume
      real(amrex_real), dimension(3), INTENT(in) :: c
      real(amrex_real), dimension(3), INTENT(out) :: derivative_theta, derivative_phi, centroid

      real(amrex_real) :: h1, h2, h3, sqh3, xh, cv, tan_theta, sec_theta, csc_phi, cot_phi
      real(amrex_real) :: g01, g02, g03, g11, g12, g13
      real(amrex_real) :: dxg02, dxg03, dyg01, dzg01, dzg02, dzg03
      real(amrex_real) :: dxg12, dxg13, dyg11, dzg11, dzg12, dzg13, dtg11, dtg12, dtg13
      real(amrex_real) :: dth1, dth2, dth3, dtxh, dtsqh3xh, dtsqh3xh1, dtsqh3xh2, dtg01q, dtg02q, dtg03q, dtg11q, dtg12q, dtg13q
      real(amrex_real) :: dfh1, dfh2, dfh3, dfxh, dfsqh3xh, dfsqh3xh1, dfsqh3xh2, dfg01q, dfg02q, dfg03q, dfg11q, dfg12q, dfg13q

      tan_theta = tan(angles(1))
      sec_theta = 1d0/cos(angles(1))
      csc_phi   = 1d0/sin(angles(2))
      cot_phi   = 1d0/tan(angles(2))

      h1   = tan_theta
      h2   = cot_phi*sec_theta
      h3   = 4d0*c(1)*c(3)*h2 - (c(1) - c(2)*h1 + c(3)*h2)**2
      sqh3 = sqrt(h3)
      cv   = c(1)*c(2)*c(3) - 2d0*volume
      xh   = cos((acos((6d0*h1*h2*cv)/sqh3**3) + PI)/3d0)

      g01 = g0(c(1)   , c(2)*h1 ,c(3)*h2)
      g02 = g0(c(2)*h1, c(1)    ,c(3)*h2)
      g03 = g0(c(3)*h2, c(1)    ,c(2)*h1)
      g11 = g1(c(1)   , c(2)*h1 ,c(3)*h2, h3)
      g12 = g1(c(2)*h1, c(1)    ,c(3)*h2, h3)
      g13 = g1(c(3)*h2, c(1)    ,c(2)*h1, h3)

      ! Centroid
      centroid(1) = ((g01 + 12d0*g11*xh*xh)/(h1*h2)    + 24d0*cv*(sqh3*xh    - 2d0*c(1)))/(192d0*volume)
      centroid(2) = ((g02 + 12d0*g12*xh*xh)/(h1*h1*h2) + 24d0*cv*(sqh3*xh/h1 - 2d0*c(2)))/(192d0*volume)
      centroid(3) = ((g03 + 12d0*g13*xh*xh)/(h1*h2*h2) + 24d0*cv*(sqh3*xh/h2 - 2d0*c(3)))/(192d0*volume)

      dxg02 = dxg0(c(2)*h1, c(1)   , c(3)*h2)
      dxg03 = dxg0(c(3)*h2, c(1)   , c(2)*h1)
      dyg01 = dyg0(c(1)   , c(2)*h1, c(3)*h2)
      dzg01 = dzg0(c(1)   , c(2)*h1, c(3)*h2)
      dzg02 = dzg0(c(2)*h1, c(1)   , c(3)*h2)
      dzg03 = dzg0(c(3)*h2, c(1)   , c(2)*h1)

      dxg12 = dxg1(c(2)*h1                  , h3)
      dxg13 = dxg1(c(3)*h2                  , h3)
      dyg11 = dyg1(         c(2)*h1, c(3)*h2, h3)
      dzg11 = dzg1(         c(2)*h1, c(3)*h2, h3)
      dzg12 = dzg1(         c(1)   , c(3)*h2, h3)
      dzg13 = dzg1(         c(1)   , c(2)*h1, h3)
      dtg11 = dtg1(c(1)   , c(2)*h1, c(3)*h2, h3)
      dtg12 = dtg1(c(2)*h1, c(1)   , c(3)*h2, h3)
      dtg13 = dtg1(c(3)*h2, c(1)   , c(2)*h1, h3)

      ! d/dΘ
      dth1 = 1d0 + h1*h1
      dth2 = h1*h2
      dth3 = 4d0*c(1)*c(3)*dth2 - 2d0*(c(1) - c(2)*h1 + c(3)*h2)*(c(3)*dth2 - c(2)*dth1)
      dtxh = cv*(3d0*h1*h2*dth3 - 2d0*h3*(h1*dth2 + h2*dth1))/(sqh3*h3*h3*(4d0*xh*xh - 1d0))
      dtsqh3xh  = (xh*dth3 + 2d0*h3*dtxh)/(2d0*sqh3)
      dtsqh3xh1 = (h1*dtsqh3xh - sqh3*xh*dth1)/(h1*h1)
      dtsqh3xh2 = (h2*dtsqh3xh - sqh3*xh*dth2)/(h2*h2)

      dtg01q = (dyg01*c(2)*dth1 + dzg01*c(3)*dth2              - (h1*dth2    + h2*dth1       )*g01/(h1*h2)   )/(h1*h2)
      dtg02q = (dxg02*c(2)*dth1 + dzg02*c(3)*dth2              - (h1*h1*dth2 + 2d0*h1*h2*dth1)*g02/(h1*h1*h2))/(h1*h1*h2)
      dtg03q = (dxg03*c(3)*dth2 + dzg03*c(2)*dth1              - (h2*h2*dth1 + 2d0*h1*h2*dth2)*g03/(h1*h2*h2))/(h1*h2*h2)
      dtg11q = (dyg11*c(2)*dth1 + dzg11*c(3)*dth2 + dtg11*dth3 - (h1*dth2    + h2*dth1       )*g11/(h1*h2)   )/(h1*h2)
      dtg12q = (dxg12*c(2)*dth1 + dzg12*c(3)*dth2 + dtg12*dth3 - (h1*h1*dth2 + 2d0*h1*h2*dth1)*g12/(h1*h1*h2))/(h1*h1*h2)
      dtg13q = (dxg13*c(3)*dth2 + dzg13*c(2)*dth1 + dtg13*dth3 - (h2*h2*dth1 + 2d0*h1*h2*dth2)*g13/(h1*h2*h2))/(h1*h2*h2)

      derivative_theta(1) = dtg01q + 24d0*cv*dtsqh3xh  + 12d0*dtg11q*xh*xh + 24d0*g11*dtxh*xh/(h1*h2)
      derivative_theta(2) = dtg02q + 24d0*cv*dtsqh3xh1 + 12d0*dtg12q*xh*xh + 24d0*g12*dtxh*xh/(h1*h1*h2)
      derivative_theta(3) = dtg03q + 24d0*cv*dtsqh3xh2 + 12d0*dtg13q*xh*xh + 24d0*g13*dtxh*xh/(h1*h2*h2)
      derivative_theta = derivative_theta/(192d0*volume)

      ! d/dφ
      dfh1 = 0d0
      dfh2 = -csc_phi**2*sec_theta
      dfh3 = 4d0*c(1)*c(3)*dfh2 - 2d0*(c(1) - c(2)*h1 + c(3)*h2)*(c(3)*dfh2 - c(2)*dfh1)
      dfxh = cv*(3d0*h1*h2*dfh3 - 2d0*h3*(h1*dfh2 + h2*dfh1))/(sqh3*h3*h3*(4d0*xh*xh - 1d0))
      dfsqh3xh  = (xh*dfh3 + 2d0*h3*dfxh)/(2d0*sqh3)
      dfsqh3xh1 = (h1*dfsqh3xh - sqh3*xh*dfh1)/(h1*h1)
      dfsqh3xh2 = (h2*dfsqh3xh - sqh3*xh*dfh2)/(h2*h2)

      dfg01q = (dyg01*c(2)*dfh1 + dzg01*c(3)*dfh2              - (h1*dfh2    + h2*dfh1       )*g01/(h1*h2)   )/(h1*h2)
      dfg02q = (dxg02*c(2)*dfh1 + dzg02*c(3)*dfh2              - (h1*h1*dfh2 + 2d0*h1*h2*dfh1)*g02/(h1*h1*h2))/(h1*h1*h2)
      dfg03q = (dxg03*c(3)*dfh2 + dzg03*c(2)*dfh1              - (h2*h2*dfh1 + 2d0*h1*h2*dfh2)*g03/(h1*h2*h2))/(h1*h2*h2)
      dfg11q = (dyg11*c(2)*dfh1 + dzg11*c(3)*dfh2 + dtg11*dfh3 - (h1*dfh2    + h2*dfh1       )*g11/(h1*h2)   )/(h1*h2)
      dfg12q = (dxg12*c(2)*dfh1 + dzg12*c(3)*dfh2 + dtg12*dfh3 - (h1*h1*dfh2 + 2d0*h1*h2*dfh1)*g12/(h1*h1*h2))/(h1*h1*h2)
      dfg13q = (dxg13*c(3)*dfh2 + dzg13*c(2)*dfh1 + dtg13*dfh3 - (h2*h2*dfh1 + 2d0*h1*h2*dfh2)*g13/(h1*h2*h2))/(h1*h2*h2)

      derivative_phi(1) = dfg01q + 24d0*cv*dfsqh3xh  + 12d0*dfg11q*xh*xh + 24d0*g11*dfxh*xh/(h1*h2)
      derivative_phi(2) = dfg02q + 24d0*cv*dfsqh3xh1 + 12d0*dfg12q*xh*xh + 24d0*g12*dfxh*xh/(h1*h1*h2)
      derivative_phi(3) = dfg03q + 24d0*cv*dfsqh3xh2 + 12d0*dfg13q*xh*xh + 24d0*g13*dfxh*xh/(h1*h2*h2)
      derivative_phi = derivative_phi/(192d0*volume)

   contains

      real(amrex_real) pure function g0(x, y, z)
         real(amrex_real), INTENT(in) :: x, y, z
         g0 = ((3d0*x - 8d0*(y + z))*x + 6d0*(8d0*y*z + (y - z)**2))*x*x - (y - z)**4
      end function g0

      real(amrex_real) pure function g1(x, y, z, t)
         real(amrex_real), INTENT(in) :: x, y, z, t
         g1 = (2d0*(x*x - (y - z)**2) - t)*t
      end function g1

      real(amrex_real) pure function dxg0(x, y, z)
         real(amrex_real), INTENT(in) :: x, y, z
         dxg0 = ((12d0*x - 24d0*(y + z))*x + 12d0*(8d0*y*z + (y - z)**2))*x
      end function dxg0

      real(amrex_real) pure function dyg0(x, y, z)
         real(amrex_real), INTENT(in) :: x, y, z
         dyg0 = (36d0*z + 12d0*y - 8d0*x)*x*x - 4d0*(y - z)**3
      end function dyg0

      real(amrex_real) pure function dzg0(x, y, z)
         real(amrex_real), INTENT(in) :: x, y, z
         dzg0 = (36d0*y + 12d0*z - 8d0*x)*x*x + 4d0*(y - z)**3
      end function dzg0

      real(amrex_real) pure function dxg1(x, t)
         real(amrex_real), INTENT(in) :: x, t
         dxg1 = 4d0*t*x
      end function dxg1

      real(amrex_real) pure function dyg1(y, z, t)
         real(amrex_real), INTENT(in) :: y, z, t
         dyg1 = 4d0*t*(z - y)
      end function dyg1

      real(amrex_real) pure function dzg1(y, z, t)
         real(amrex_real), INTENT(in) :: y, z, t
         dzg1 = 4d0*t*(y - z)
      end function dzg1

      real(amrex_real) pure function dtg1(x, y, z, t)
         real(amrex_real), INTENT(in) :: x, y, z, t
         dtg1 = 2d0*(x*x - (y - z)**2 - t)
      end function dtg1
   end subroutine mof3d_derivatives_hexa

   ! map x in [0, 2π[
   real(amrex_real) pure function modulo_tau(x) result(r)
      real(amrex_real), INTENT(in) :: x

      real(amrex_real), parameter :: TAU = 4d0*acos(0d0)

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

#undef DEBUG_ANALYTICAL_GRAD
