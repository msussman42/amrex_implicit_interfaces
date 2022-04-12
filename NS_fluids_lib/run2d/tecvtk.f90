      PROGRAM tecvtk

      IMPLICIT NONE

      integer :: Nodes,Cells,i,j
      integer :: n1,n2,n3,n4
      integer :: stat
      real, allocatable :: node_list(:,:)
      character(80) :: discard
      character(80) :: points_line
      real :: node_center(3)
      real :: node_min(3)
      real :: node_max(3)
      real :: max_side_len
      real :: min_side_len
      real :: buffer_factor
      real :: interior_radius
      integer :: dir

      print *,"This program reads input.vtk"
      print *,"and outputs an ascii tecplot file: output.plt"
      open(unit=2, file= 'input.vtk',status='old',iostat=stat)
      if (stat.ne.0) then
       print *,"input.vtk can not be opened"
       stop
      endif

      open(unit=4, file= 'output.plt',status='unknown',iostat=stat)
      if (stat.ne.0) then
       print *,"output.plt can not be opened"
       stop
      endif

      do i=1,4
       read(2,*) discard
       print *,"line ",i
       print *,discard
      enddo
      read(2,'(a6)',advance='no') points_line
      read(2,*) Nodes

      print *,"number of nodes=",Nodes

      allocate(node_list(Nodes,3))
      do i=1,Nodes
       read(2,*) node_list(i,1),node_list(i,2),node_list(i,3)
       print *,"Node i=",i,"x,y,z=", &
         node_list(i,1),node_list(i,2),node_list(i,3)
      enddo
      read(2,*) discard,Cells
      print *,"number of elements=",Cells

!     Nodes is number of nodes, Cells is number of cells
!     Shape: 3=Triangle, 4=Quadrilateral

      write(4,*) 'TITLE = "3D surface" '
      write(4,*) 'VARIABLES = "X", "Y", "Z" '
      if (1.eq.0) then
       write(4,*) 'ZONE T="Triangles", NODES= ', Nodes, ', ELEMENTS= ', &
        Cells, ', DATAPACKING=POINT, '
      else
       write(4,*) 'ZONE T="TRIANGLES", N= ', Nodes, ', E= ', &
        Cells, ', DATAPACKING=POINT, '
      endif

      write(4,*) 'ZONETYPE=FETRIANGLE' 

      print *,"outputting nodes..." 
      do i=1,Nodes
       write(4,*)  node_list(i,1), node_list(i,2), node_list(i,3)
      enddo
      
      print *,"inputting and outputting elements (cells)..." 

      do i=1,Cells
      
        read(2,*) n1,n2,n3,n4
        if (n1.eq.3) then
         ! do nothing
        else
         print *,"n1 invalid"
         stop
        endif
        write(4,*) n2+1, n3+1, n4+1
        if (1.eq.0) then
         print *,"Element i=",i,"n2+1,n3+1,n4+1=",n2+1,n3+1,n4+1
        endif

      end do  

      do dir=1,3
       node_center(dir)=0.0d0
       node_min(dir)=node_list(1,dir)
       node_max(dir)=node_list(1,dir)
      enddo

      do i=1,Nodes
       do dir=1,3
        node_center(dir)=node_center(dir)+node_list(i,dir)
        node_min(dir)=min(node_min(dir),node_list(i,dir))
        node_max(dir)=max(node_max(dir),node_list(i,dir))
       enddo
      enddo
      max_side_len=node_max(1)-node_min(1)
      min_side_len=node_max(1)-node_min(1)
      buffer_factor=0.2d0
      do dir=1,3
       node_center(dir)=node_center(dir)/Nodes
       max_side_len=max(max_side_len,node_max(dir)-node_min(dir))
       min_side_len=min(min_side_len,node_max(dir)-node_min(dir))
      enddo
      print *,"min_side_len,max_side_len ",min_side_len,max_side_len
      do dir=1,3
       print *,"dir,node_min ",dir,node_min(dir)
       print *,"dir,node_max ",dir,node_max(dir)
       print *,"dir,node_center ",dir,node_center(dir)
       print *,"recommended sign box  dir,lo,hi ",dir, &
               node_min(dir)-0.25d0*buffer_factor*max_side_len, &
               node_max(dir)+0.25d0*buffer_factor*max_side_len
       print *,"recommended interior sign box dir,lo,hi ",dir, &
               node_min(dir)+0.2d0*min_side_len, &
               node_max(dir)-0.2d0*min_side_len
      enddo
      print *,"recommended sign box LS: ", &
              -0.25d0*buffer_factor*max_side_len

      deallocate(node_list)
      close(2)
      close(4)

      END PROGRAM
