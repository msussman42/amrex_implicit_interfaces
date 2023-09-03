!This file is part of Notus 0.5.0

!Copyright Bordeaux-INP, Université de Bordeaux, CNRS
!
!Contributors:
!Antoine Lemoine, 02-02-2016, antoine.lemoine@bordeaux-inp.fr

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

!> @defgroup polyhedron_group Polyhedron
!! @brief Geometric tools relative to polyhedron
!! @ingroup computational_geometry_3d

module mod_cg3_polyhedron
   use mod_cg3_points
   implicit none

   !> Definition of an incidence matrix
   !! @ingroup polyhedron_group
   type t_incidence_matrix
      !> List of objects
      integer, dimension(:), allocatable :: id
      !> Number of objects
      integer :: size
   end type t_incidence_matrix

   !> Definition of a 3-dimensional polyhedron
   !! @ingroup polyhedron_group
   type t_polyhedron
      !> Coordinates of the polyhedron vertices
      double precision, dimension(:,:), allocatable :: point
      !> List of the edges defined by the indices of their two nodes (such that edge(1,i) < edge(2,i))
      integer, dimension(:,:), allocatable :: edge
      !> List of the faces defined by the indices of their points defined in counterclockwise order
      type(t_incidence_matrix), dimension(:), allocatable :: face

      !> List of the faces defined by the indices of their edges defined in counterclockwise order
      type(t_incidence_matrix), dimension(:), allocatable :: face_to_edge
      !> List of neighbor faces of an edge
      integer, dimension(:,:), allocatable :: edge_to_face
      !> List of neighbor edges of a point
      type(t_incidence_matrix), dimension(:), allocatable :: point_to_edge

      !> List of edge tangents
      double precision, dimension(:,:), allocatable :: tangent
      !> List of face normals
      double precision, dimension(:,:), allocatable :: normal

      !> Number of points
      integer :: nb_points = 0
      !> Number of edges
      integer :: nb_edges = 0
      !> Number of faces
      integer :: nb_faces = 0
   end type t_polyhedron

   interface unalloc
      module procedure :: cg3_finalize_polyhedron
   end interface unalloc

contains

   !> Finalize a polyhedron
   !!
   !! @param[inout] polyhedron: any polyhedron
   !! @ingroup polyhedron_group
   subroutine cg3_finalize_polyhedron(polyhedron)
      type(t_polyhedron), intent(inout) :: polyhedron

      integer :: i

      if (allocated(polyhedron%point)) deallocate(polyhedron%point)
      if (allocated(polyhedron%edge)) deallocate(polyhedron%edge)
      if (allocated(polyhedron%face)) then
         do i = 1, size(polyhedron%face)
            if (allocated(polyhedron%face(i)%id)) deallocate(polyhedron%face(i)%id)
         end do
         deallocate(polyhedron%face)
      end if
      if (allocated(polyhedron%face_to_edge)) then
         do i = 1, size(polyhedron%face_to_edge)
            if (allocated(polyhedron%face_to_edge(i)%id)) deallocate(polyhedron%face_to_edge(i)%id)
         end do
         deallocate(polyhedron%face_to_edge)
      end if
      if (allocated(polyhedron%edge_to_face)) deallocate(polyhedron%edge_to_face)
      if (allocated(polyhedron%point_to_edge)) then
         do i = 1, size(polyhedron%point_to_edge)
            if (allocated(polyhedron%point_to_edge(i)%id)) deallocate(polyhedron%point_to_edge(i)%id)
         end do
         deallocate(polyhedron%point_to_edge)
      end if
      polyhedron%nb_points = 0
      polyhedron%nb_edges = 0
      polyhedron%nb_faces = 0
   end subroutine cg3_finalize_polyhedron

   !> Move allocation from source to target
   !!
   !! @param[inout] polyhedron_source: source polyhedron
   !! @param[inout] polyhedron_target: target polyhedron
   !! @ingroup polyhedron_group
   pure subroutine cg3_polyhedron_move_alloc(polyhedron_source, polyhedron_target)
      type(t_polyhedron), intent(inout) :: polyhedron_source
      type(t_polyhedron), intent(inout) :: polyhedron_target

      call move_alloc(polyhedron_source%point, polyhedron_target%point)
      call move_alloc(polyhedron_source%edge, polyhedron_target%edge)
      call move_alloc(polyhedron_source%face, polyhedron_target%face)

      call move_alloc(polyhedron_source%face_to_edge, polyhedron_target%face_to_edge)
      call move_alloc(polyhedron_source%edge_to_face, polyhedron_target%edge_to_face)
      call move_alloc(polyhedron_source%point_to_edge, polyhedron_target%point_to_edge)

      call move_alloc(polyhedron_source%tangent, polyhedron_target%tangent)
      call move_alloc(polyhedron_source%normal, polyhedron_target%normal)

      polyhedron_target%nb_points = polyhedron_source%nb_points
      polyhedron_target%nb_edges = polyhedron_source%nb_edges
      polyhedron_target%nb_faces = polyhedron_source%nb_faces
   end subroutine cg3_polyhedron_move_alloc

   !> Write a polyhedron to a VTK file
   !!
   !! @param[in] polyhedron: any initialized polyhedron
   !! @param[in] filename: VTK file name (with extension)
   !! @ingroup polyhedron_group
   subroutine cg3_polyhedron_write_vtk_file(polyhedron, filename)
      type(t_polyhedron), intent(in) :: polyhedron
      character(len=*), intent(in) :: filename

      integer :: gridunit, i, j, nb_face_points

      open(newunit=gridunit, file=trim(filename), status="replace")

      write(gridunit,'("# vtk DataFile Version 3.0")')
      write(gridunit,'(a)') trim(filename)
      write(gridunit,'("ASCII")')
      write(gridunit,'("DATASET POLYDATA")')
      write(gridunit,*)
      write(gridunit,'("POINTS ",i8," float")') polyhedron%nb_points

      ! Write the list of points
      do i = 1, polyhedron%nb_points
         write(gridunit,*) polyhedron%point(:,i)
      end do

      ! Count the number of integer for the faces
      nb_face_points = 0
      do i = 1, polyhedron%nb_faces
         nb_face_points = nb_face_points + polyhedron%face(i)%size + 1
      end do

      write(gridunit,'("POLYGONS ",i8," ",i8)') polyhedron%nb_faces, nb_face_points

      do j = 1, polyhedron%nb_faces
         write(gridunit,'(i8)',advance="no") polyhedron%face(j)%size
         do i = 1, polyhedron%face(j)%size
            write(gridunit,'(i8)',advance="no") polyhedron%face(j)%id(i) - 1
         end do
         write(gridunit,*)
      end do

      close(gridunit)
   end subroutine cg3_polyhedron_write_vtk_file

   !> Write a polyhedron to a OBJ Wavefront file
   !!
   !! @param[in] polyhedron: any initialized polyhedron
   !! @param[in] filename: OBJ file name (with extension)
   !! @ingroup polyhedron_group
   subroutine cg3_polyhedron_write_obj_file(polyhedron, filename)
      type(t_polyhedron), intent(in) :: polyhedron
      character(len=*), intent(in) :: filename

      integer :: gridunit, i

      open(newunit=gridunit, file=trim(filename), status="unknown")

      write(gridunit,'("####")')
      write(gridunit,'("#")')
      write(gridunit,'("# OBJ file generated by Notus")')
      write(gridunit,'("#")')
      write(gridunit,'("####")')
      write(gridunit,'("# Object ",a)') trim(filename)
      write(gridunit,'("#")')
      write(gridunit,'("# Vertices: ",g0)') polyhedron%nb_points
      write(gridunit,'("# Faces: ",g0)') polyhedron%nb_faces
      write(gridunit,'("####")')

      ! Write the list of points
      do i = 1, polyhedron%nb_points
         !write(gridunit,'("v ",f19.16," ",f19.16," ",f19.16)') polyhedron%point(:,i)
         write(gridunit,'("v ",g0," ",g0," ",g0)') polyhedron%point(:,i)
      end do

      write(gridunit,'("#")')

      ! Count the number of integer for the faces
      do i = 1, polyhedron%nb_faces
         write(gridunit,'("f ",*(g0,:," "))') polyhedron%face(i)%id
      end do

      write(gridunit,'("#")')
      write(gridunit,*)
      write(gridunit,'("# End of file")')

      close(gridunit)
   end subroutine cg3_polyhedron_write_obj_file

   !> Compute the volume of a polyhedron
   !!
   !! @param[in]  polyhedron: any initialized polyhedron
   !! @param[out] volume: volume of the polyhedron
   !! @ingroup polyhedron_group
   pure subroutine cg3_polyhedron_compute_volume(polyhedron, volume)
      type(t_polyhedron), intent(in) :: polyhedron
      double precision, intent(out) :: volume

      integer :: i, j
      double precision, dimension(3) :: area, origin

      volume = 0d0

      do j = 1, polyhedron%nb_faces
         area = 0d0
         origin = polyhedron%point(:,polyhedron%face(j)%id(1))
         do i = 2, polyhedron%face(j)%size-1
            area = area + cg3_cross_product(                                &
               &   polyhedron%point(:,polyhedron%face(j)%id(i))   - origin, &
               &   polyhedron%point(:,polyhedron%face(j)%id(i+1)) - origin)
         end do
         volume = volume + dot_product(area, origin)
      end do

      volume = volume/6d0
   end subroutine cg3_polyhedron_compute_volume

   !> Compute the volume and the centroid of a polyhedron
   !!
   !! @param[in]  polyhedron: any initialized polyhedron
   !! @param[out] volume: volume of the polyhedron
   !! @param[out] centroid: centroid of the polyhedron
   !! @ingroup polyhedron_group
   pure subroutine cg3_polyhedron_compute_centroid(polyhedron, volume, centroid)
      type(t_polyhedron), intent(in) :: polyhedron
      double precision, intent(out) :: volume
      double precision, dimension(3), intent(out) :: centroid

      integer :: i, j
      double precision, dimension(3) :: local_area, area, origin, left, right, a, b, c

      volume = 0d0
      centroid = 0d0

      do j = 1, polyhedron%nb_faces
         area = 0d0
         origin = polyhedron%point(:,polyhedron%face(j)%id(1))
         do i = 2, polyhedron%face(j)%size-1
            left = polyhedron%point(:,polyhedron%face(j)%id(i))
            right = polyhedron%point(:,polyhedron%face(j)%id(i+1))
            local_area = cg3_cross_product(left  - origin, right - origin)
            a = origin + left
            b = left + right
            c = origin + right
            area = area + local_area
            centroid = centroid + local_area*(a**2 + b**2 + c**2)
         end do
         volume = volume + dot_product(area, origin)
      end do

      centroid = centroid/(8d0*volume)
      volume = volume/6d0
   end subroutine cg3_polyhedron_compute_centroid

   !> Compute the faces normals of the polyhedron
   !!
   !! @param[inout] polyhedron: any polyhedron
   !! @param[out]   zero_area_face: contains the face id with zero area. 0 otherwise.
   !! @ingroup polyhedron_group
   pure subroutine cg3_polyhedron_compute_normals(polyhedron, zero_area_face)
      type(t_polyhedron), intent(inout) :: polyhedron
      integer, intent(out) :: zero_area_face

      integer :: i, j
      double precision, dimension(3) :: center
      double precision :: norm

      zero_area_face = 0

      if (allocated(polyhedron%normal)) then
         if (size(polyhedron%normal) /= polyhedron%nb_faces) then
            deallocate(polyhedron%normal)
            allocate(polyhedron%normal(3,polyhedron%nb_faces))
         end if
      else
         allocate(polyhedron%normal(3,polyhedron%nb_faces))
      end if

      polyhedron%normal = 0d0

      do i = 1, polyhedron%nb_faces
         ! Compute the center of the face vertices (not the centroid of the face)
         center = 0d0
         do j = 1, polyhedron%face(i)%size
            center = center + polyhedron%point(:,polyhedron%face(i)%id(j))
         end do
         center = center/polyhedron%face(i)%size

         ! Compute the normal of each sub-triangle of the cell
         do j = 1, polyhedron%face(i)%size - 1
            polyhedron%normal(:,i) = polyhedron%normal(:,i) + cg3_cross_product(              &
               &                     polyhedron%point(:,polyhedron%face(i)%id(j  )) - center, &
               &                     polyhedron%point(:,polyhedron%face(i)%id(j+1)) - center)
         end do

         ! Normalize
         norm = norm2(polyhedron%normal(:,i))

         if (norm < tiny(1d0)) then
            zero_area_face = i
            polyhedron%normal(:,i) = [1d0, 0d0, 0d0]
         end if

         polyhedron%normal(:,i) = polyhedron%normal(:,i)/norm
      end do
   end subroutine cg3_polyhedron_compute_normals

   !> Compute the edges tangents of the polyhedron
   !!
   !! @param[inout] polyhedron: any polyhedron
   !! @param[out]   zero_length_edge: contains the edge id with zero length. 0 otherwise.
   !! @ingroup polyhedron_group
   pure subroutine cg3_polyhedron_compute_tangents(polyhedron, zero_length_edge)
      type(t_polyhedron), intent(inout) :: polyhedron
      integer, intent(out) :: zero_length_edge

      integer :: i
      double precision :: norm

      zero_length_edge = 0

      if (allocated(polyhedron%tangent)) then
         if (size(polyhedron%tangent) /= polyhedron%nb_edges) then
            deallocate(polyhedron%tangent)
            allocate(polyhedron%tangent(3,polyhedron%nb_edges))
         end if
      else
         allocate(polyhedron%tangent(3,polyhedron%nb_edges))
      end if

      do i = 1, polyhedron%nb_edges
         polyhedron%tangent(:,i) = polyhedron%point(:,polyhedron%edge(2,i)) - polyhedron%point(:,polyhedron%edge(1,i))

         ! Normalize
         norm = norm2(polyhedron%tangent(:,i))

         if (norm < tiny(1d0)) then
            zero_length_edge = i
            return
         end if

         polyhedron%tangent(:,i) = polyhedron%tangent(:,i)/norm
      end do
   end subroutine cg3_polyhedron_compute_tangents

   pure subroutine cg3_create_tetrahedron(p1, p2, p3, p4, tetra)
      double precision, dimension(3), intent(in) :: p1, p2, p3, p4
      type(t_polyhedron), intent(out) :: tetra

      integer :: dum

      tetra%nb_points = 4
      tetra%nb_edges = 6
      tetra%nb_faces = 4

      allocate(tetra%point(3, tetra%nb_points))

      tetra%point(:,1) = p1
      tetra%point(:,2) = p2
      tetra%point(:,3) = p3
      tetra%point(:,4) = p4

      allocate(tetra%face(tetra%nb_faces))

      tetra%face(1)%size = 3
      allocate(tetra%face(1)%id(tetra%face(1)%size))
      tetra%face(1)%id = [2, 3, 4]
      tetra%face(2)%size = 3
      allocate(tetra%face(2)%id(tetra%face(2)%size))
      tetra%face(2)%id = [3, 2, 1]
      tetra%face(3)%size = 3
      allocate(tetra%face(3)%id(tetra%face(3)%size))
      tetra%face(3)%id = [4, 1, 2]
      tetra%face(4)%size = 3
      allocate(tetra%face(4)%id(tetra%face(4)%size))
      tetra%face(4)%id = [1, 4, 3]

      allocate(tetra%edge(2,tetra%nb_edges))

      tetra%edge(:,1) = [1, 2]
      tetra%edge(:,2) = [1, 3]
      tetra%edge(:,3) = [1, 4]
      tetra%edge(:,4) = [2, 3]
      tetra%edge(:,5) = [2, 4]
      tetra%edge(:,6) = [3, 4]

      allocate(tetra%edge_to_face(2,tetra%nb_edges))

      tetra%edge_to_face(:,1) = [3, 2]
      tetra%edge_to_face(:,2) = [2, 4]
      tetra%edge_to_face(:,3) = [3, 4]
      tetra%edge_to_face(:,4) = [1, 2]
      tetra%edge_to_face(:,5) = [1, 3]
      tetra%edge_to_face(:,6) = [4, 1]

      allocate(tetra%face_to_edge(tetra%nb_faces))

      tetra%face_to_edge(1)%size = 3
      allocate(tetra%face_to_edge(1)%id(3))
      tetra%face_to_edge(1)%id = [4, 6, 5]
      tetra%face_to_edge(2)%size = 3
      allocate(tetra%face_to_edge(2)%id(3))
      tetra%face_to_edge(2)%id = [2, 4, 1]
      tetra%face_to_edge(3)%size = 3
      allocate(tetra%face_to_edge(3)%id(3))
      tetra%face_to_edge(3)%id = [1, 5, 3]
      tetra%face_to_edge(4)%size = 3
      allocate(tetra%face_to_edge(4)%id(3))
      tetra%face_to_edge(4)%id = [3, 6, 2]

      allocate(tetra%point_to_edge(tetra%nb_points))

      tetra%point_to_edge(1)%size = 3
      allocate(tetra%point_to_edge(1)%id(3))
      tetra%point_to_edge(1)%id = [1, 2, 3]
      tetra%point_to_edge(2)%size = 3
      allocate(tetra%point_to_edge(2)%id(3))
      tetra%point_to_edge(2)%id = [1, 4, 5]
      tetra%point_to_edge(3)%size = 3
      allocate(tetra%point_to_edge(3)%id(3))
      tetra%point_to_edge(3)%id = [2, 4, 6]
      tetra%point_to_edge(4)%size = 3
      allocate(tetra%point_to_edge(4)%id(3))
      tetra%point_to_edge(4)%id = [3, 5, 6]

      call cg3_polyhedron_compute_normals(tetra, dum)
      call cg3_polyhedron_compute_tangents(tetra, dum)
   end subroutine cg3_create_tetrahedron

   ! Create a cuboid from cell dimensions
   subroutine cg3_create_cuboid(c, cuboid)
      double precision, dimension(3), intent(in) :: c
      type(t_polyhedron), intent(out) :: cuboid

      cuboid%nb_points = 8
      cuboid%nb_edges = 12
      cuboid%nb_faces = 6

      ! Generate the points
      allocate(cuboid%point(3, cuboid%nb_points))

      cuboid%point(:,1) = [0d0 , 0d0 , 0d0 ]
      cuboid%point(:,2) = [0d0 , 0d0 , c(3)]
      cuboid%point(:,3) = [0d0 , c(2), 0d0 ]
      cuboid%point(:,4) = [0d0 , c(2), c(3)]
      cuboid%point(:,5) = [c(1), 0d0 , 0d0 ]
      cuboid%point(:,6) = [c(1), 0d0 , c(3)]
      cuboid%point(:,7) = [c(1), c(2), 0d0 ]
      cuboid%point(:,8) = [c(1), c(2), c(3)]

      ! Generate the faces
      allocate(cuboid%face(cuboid%nb_faces))

      cuboid%face(1)%size = 4
      allocate(cuboid%face(1)%id(cuboid%face(1)%size))
      cuboid%face(1)%id = [8, 4, 2, 6]
      cuboid%face(2)%size = 4
      allocate(cuboid%face(2)%id(cuboid%face(2)%size))
      cuboid%face(2)%id= [8, 6, 5, 7]
      cuboid%face(3)%size = 4
      allocate(cuboid%face(3)%id(cuboid%face(3)%size))
      cuboid%face(3)%id = [8, 7, 3, 4]
      cuboid%face(4)%size = 4
      allocate(cuboid%face(4)%id(cuboid%face(4)%size))
      cuboid%face(4)%id = [4, 3, 1, 2]
      cuboid%face(5)%size = 4
      allocate(cuboid%face(5)%id(cuboid%face(5)%size))
      cuboid%face(5)%id = [1, 3, 7, 5]
      cuboid%face(6)%size = 4
      allocate(cuboid%face(6)%id(cuboid%face(6)%size))
      cuboid%face(6)%id = [2, 1, 5, 6]

      ! Generate the edges
      allocate(cuboid%edge(2,cuboid%nb_edges))

      cuboid%edge(:,1 ) = [1, 2]
      cuboid%edge(:,2 ) = [1, 3]
      cuboid%edge(:,3 ) = [1, 5]
      cuboid%edge(:,4 ) = [2, 4]
      cuboid%edge(:,5 ) = [2, 6]
      cuboid%edge(:,6 ) = [3, 4]
      cuboid%edge(:,7 ) = [3, 7]
      cuboid%edge(:,8 ) = [4, 8]
      cuboid%edge(:,9 ) = [5, 6]
      cuboid%edge(:,10) = [5, 7]
      cuboid%edge(:,11) = [6, 8]
      cuboid%edge(:,12) = [7, 8]

      ! Edges -> faces incidence matrix
      allocate(cuboid%edge_to_face(2,cuboid%nb_edges))

      cuboid%edge_to_face(:, 1) = [6, 4]
      cuboid%edge_to_face(:, 2) = [5, 4]
      cuboid%edge_to_face(:, 3) = [6, 5]
      cuboid%edge_to_face(:, 4) = [4, 1]
      cuboid%edge_to_face(:, 5) = [6, 1]
      cuboid%edge_to_face(:, 6) = [3, 4]
      cuboid%edge_to_face(:, 7) = [5, 3]
      cuboid%edge_to_face(:, 8) = [1, 3]
      cuboid%edge_to_face(:, 9) = [2, 6]
      cuboid%edge_to_face(:,10) = [5, 2]
      cuboid%edge_to_face(:,11) = [2, 1]
      cuboid%edge_to_face(:,12) = [2, 3]

      ! Faces -> edges incidence matrix
      allocate(cuboid%face_to_edge(cuboid%nb_faces))

      cuboid%face_to_edge(1)%size = 4
      allocate(cuboid%face_to_edge(1)%id(4))
      cuboid%face_to_edge(1)%id(:) = [5,11,8,4]
      cuboid%face_to_edge(2)%size = 4
      allocate(cuboid%face_to_edge(2)%id(4))
      cuboid%face_to_edge(2)%id(:) = [10,12,11,9]
      cuboid%face_to_edge(3)%size = 4
      allocate(cuboid%face_to_edge(3)%id(4))
      cuboid%face_to_edge(3)%id(:) = [6,8,12,7]
      cuboid%face_to_edge(4)%size = 4
      allocate(cuboid%face_to_edge(4)%id(4))
      cuboid%face_to_edge(4)%id(:) = [1,4,6,2]
      cuboid%face_to_edge(5)%size = 4
      allocate(cuboid%face_to_edge(5)%id(4))
      cuboid%face_to_edge(5)%id(:) = [2,7,10,3]
      cuboid%face_to_edge(6)%size = 4
      allocate(cuboid%face_to_edge(6)%id(4))
      cuboid%face_to_edge(6)%id(:) = [3,9,5,1]

      ! Vertices -> edges incidence matrix
      allocate(cuboid%point_to_edge(cuboid%nb_points))

      cuboid%point_to_edge(1)%size = 3
      allocate(cuboid%point_to_edge(1)%id(3))
      cuboid%point_to_edge(1)%id(:) = [1,2,3]
      cuboid%point_to_edge(2)%size = 3
      allocate(cuboid%point_to_edge(2)%id(3))
      cuboid%point_to_edge(2)%id(:) = [1,4,5]
      cuboid%point_to_edge(3)%size = 3
      allocate(cuboid%point_to_edge(3)%id(3))
      cuboid%point_to_edge(3)%id(:) = [2,6,7]
      cuboid%point_to_edge(4)%size = 3
      allocate(cuboid%point_to_edge(4)%id(3))
      cuboid%point_to_edge(4)%id(:) = [4,6,8]
      cuboid%point_to_edge(5)%size = 3
      allocate(cuboid%point_to_edge(5)%id(3))
      cuboid%point_to_edge(5)%id(:) = [3,9,10]
      cuboid%point_to_edge(6)%size = 3
      allocate(cuboid%point_to_edge(6)%id(3))
      cuboid%point_to_edge(6)%id(:) = [5,9,11]
      cuboid%point_to_edge(7)%size = 3
      allocate(cuboid%point_to_edge(7)%id(3))
      cuboid%point_to_edge(7)%id(:) = [7,10,12]
      cuboid%point_to_edge(8)%size = 3
      allocate(cuboid%point_to_edge(8)%id(3))
      cuboid%point_to_edge(8)%id(:) = [8,11,12]

      ! Tangents
      allocate(cuboid%tangent(3,cuboid%nb_edges))

      cuboid%tangent(:, 1) = [0d0, 0d0, 1d0]
      cuboid%tangent(:, 2) = [0d0, 1d0, 0d0]
      cuboid%tangent(:, 3) = [1d0, 0d0, 0d0]
      cuboid%tangent(:, 4) = [0d0, 1d0, 0d0]
      cuboid%tangent(:, 5) = [1d0, 0d0, 0d0]
      cuboid%tangent(:, 6) = [0d0, 0d0, 1d0]
      cuboid%tangent(:, 7) = [1d0, 0d0, 0d0]
      cuboid%tangent(:, 8) = [1d0, 0d0, 0d0]
      cuboid%tangent(:, 9) = [0d0, 0d0, 1d0]
      cuboid%tangent(:,10) = [0d0, 1d0, 0d0]
      cuboid%tangent(:,11) = [0d0, 1d0, 0d0]
      cuboid%tangent(:,12) = [0d0, 0d0, 1d0]

      ! Normals
      allocate(cuboid%normal(3,cuboid%nb_faces))

      cuboid%normal(:, 1) = [ 0d0, 0d0, 1d0]
      cuboid%normal(:, 2) = [ 1d0, 0d0, 0d0]
      cuboid%normal(:, 3) = [ 0d0, 1d0, 0d0]
      cuboid%normal(:, 4) = [-1d0, 0d0, 0d0]
      cuboid%normal(:, 5) = [ 0d0, 0d0,-1d0]
      cuboid%normal(:, 6) = [ 0d0,-1d0, 0d0]
   end subroutine cg3_create_cuboid

end module mod_cg3_polyhedron
