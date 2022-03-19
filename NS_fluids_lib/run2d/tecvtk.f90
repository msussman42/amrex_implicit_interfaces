      PROGRAM tecvtk

      IMPLICIT NONE

      integer :: Nodes,Cells,i,j
      integer :: n1,n2,n3,n4
      integer :: stat
      real, allocatable :: node_list(:,:)
      character(80) :: discard
      character(80) :: points_line

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

      end do
      deallocate(node_list)
      
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
        print *,"Element i=",i,"n2+1,n3+1,n4+1=",n2+1,n3+1,n4+1

      end do  

      close(2)
      close(4)

      END PROGRAM
