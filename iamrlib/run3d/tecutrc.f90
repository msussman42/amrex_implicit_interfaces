      PROGRAM tecwhale


      IMPLICIT NONE

      Integer :: Nodes, Cells, Shape, N, i, garbage
      Integer :: n1, n2, n3, n4
      Real :: X, Y, Z

      print *,"This program reads injectorgeom.dat"
      print *,"and outputs a tecplot file: output.dat"
      open(unit=2, file= 'injectorgeom.dat')
      open(unit=4, file= 'output.dat')

      print *,"reading number of nodes and cells ..."
      read(2,*) Nodes,Cells
      print *,"number of nodes=",Nodes
      print *,"number of cells=",Cells

!     Nodes is number of nodes, Cells is number of cells

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

       read(2,*) x, y, z
       write(4,*)  x, y, z

      end do
      
      print *,"inputting and outputting elements (cells)..." 

      do i=1,Cells-1
      
        read(2,*) n1, n2, n3, n4
        write(4,*) n1, n2, n3

      end do  

      

      read(2,*) n1,n2,n3,n4
      write(4,*) n1, n2, n3

      close(2)
      close(4)

      END PROGRAM
