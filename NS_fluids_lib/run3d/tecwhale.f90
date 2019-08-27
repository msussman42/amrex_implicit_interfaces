      PROGRAM tecwhale


      IMPLICIT NONE

      Integer :: Nodes, Cells, Shape, N, i, garbage
      Integer :: n1, n2, n3
      Real :: X, Y, Z

      print *,"This program reads input.txt"
      print *,"and outputs a tecplot file: output.dat"
      open(unit=2, file= 'input.txt')
      open(unit=4, file= 'output.dat')

      print *,"reading number of nodes..."
      read(2,*) Nodes
      print *,"number of nodes=",Nodes
      print *,"reading number of elements (cells)..."
      read(2,*) Cells
      print *,"number of cells=",Cells
      print *,"reading number of nodes per element..."
      read(2,*) Shape
      print *,"number of nodes per element=",Shape
      if (Shape.ne.3) then
       print *,"this code is for triangles, use tecquad.f90 instead"
       stop
      endif

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

      print *,"inputting and outputting nodes..." 
      do i=1,Nodes

       read(2,*) N, x, y, z
       write(4,*)  x, y, z

      end do
      

      do i=1,3
        read(2,*) garbage
      end do

      print *,"inputting and outputting elements (cells)..." 

      do i=1,Cells-1
      
        read(2,*) n1
        read(2,*) n2
        read(2,*) n3
        write(4,*) n1, n2, n3

        read(2,*) garbage
        read(2,*) garbage

      end do  

      

      read(2,*) n1
      read(2,*) n2
      read(2,*) n3
      write(4,*) n1, n2, n3



      END PROGRAM
