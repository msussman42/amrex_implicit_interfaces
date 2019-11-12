      PROGRAM tecwhale


      IMPLICIT NONE

      Integer :: Nodes, Cells, Shape, N, i, garbage
      Integer :: n1, n2, n3,n4
      Real :: X, Y, Z

      print *,"This program reads 2D input.cas"
      print *,"and outputs a 2D tecplot file: output.plt"
      open(unit=2, file= 'input.cas')
      open(unit=4, file= 'output.plt')

      print *,"reading number of nodes and number of lines ..."
      read(2,*) Nodes,Cells
      print *,"number of nodes=",Nodes
      print *,"number of lines=",Cells

!     Nodes is number of nodes, Cells is number of lines

      write(4,*) 'TITLE = "2D surface" '
      write(4,*) 'VARIABLES = "X", "Y" '
      write(4,*) 'ZONE T="TRIANGLES", N= ', Nodes, ', E= ', &
        Cells, ', DATAPACKING=POINT, '

      write(4,*) 'ZONETYPE=FELINESEG' 

      print *,"inputting and outputting nodes..." 
      do i=1,Nodes

        ! Yang had a z coordinate for the ginger bread file
       read(2,*) x, y, z
       write(4,*)  x, y

      end do
      
      print *,"inputting and outputting elements (cells)..." 

      do i=1,Cells
      
        read(2,*) n1,n2
        write(4,*) n1, n2

      end do  

      END PROGRAM
