      PROGRAM tecwhale


      IMPLICIT NONE

      Integer :: Nodes, Cells, Shape, N, i, garbage
      Integer :: n1, n2, n3, n4
      Real :: X, Y, Z

      print *,"This program reads naca0012.sci"
      print *,"and outputs a tecplot file: naca0012.tec"
      open(unit=2, file= 'naca0012.sci')
      open(unit=4, file= 'naca0012.tec')

      print *,"reading number of nodes..."
      read(2,*) Nodes
      print *,"number of nodes=",Nodes
      print *,"reading number of elements (cells)..."
      read(2,*) Cells
      print *,"number of cells=",Cells
      Shape=4
      print *,"number of nodes per element=",Shape
      if (Shape.ne.4) then
       print *,"this code is for quads"
       stop
      endif

!     Nodes is number of nodes, Cells is number of cells
!     Shape: 3=Triangle, 4=Quadrilateral

      write(4,*) 'TITLE = "3D surface" '
      write(4,*) 'VARIABLES = "X", "Y", "Z" '
      if (1.eq.0) then
       write(4,*) 'ZONE T="Quadrilaterals", NODES= ', Nodes, &
        ', ELEMENTS= ', &
        Cells, ', DATAPACKING=POINT, '
      else
       write(4,*) 'ZONE T="Quadrilaterals", N= ', Nodes, ', E= ', &
        Cells, ', DATAPACKING=POINT, '
      endif

      write(4,*) 'ZONETYPE=FEQUADRILATERAL' 

      print *,"inputting and outputting nodes..." 
      do i=1,Nodes

       read(2,*)  x, y, z
       write(4,*)  x, y, z

      end do
      

      print *,"inputting and outputting elements (cells)..." 

      do i=1,Cells-1
      
        read(2,*) n1,n2,n3,n4
        write(4,*) n1, n2, n3, n4

      end do  

      read(2,*) n1,n2,n3,n4
      write(4,*) n1, n2, n3, n4


      END PROGRAM
