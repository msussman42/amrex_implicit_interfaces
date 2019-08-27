      PROGRAM tecwhale


      IMPLICIT NONE

      Integer :: Nodes, Cells, Shape, N, i, garbage
      Integer :: n1, n2, n3, n4,numquads,dir
      Real :: X, Y, Z
      Real :: minnode(3),maxnode(3)

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
      do dir=1,3
       minnode(dir)=1.0E+20
       maxnode(dir)=-1.0E+20
      enddo

      print *,"first the program will figure out how many quads to split"
      numquads=0
      print *,"leafing through nodes..." 
      do i=1,Nodes
       read(2,*) N, x, y, z
       if (x.lt.minnode(1)) then
        minnode(1)=x
       endif
       if (x.gt.maxnode(1)) then
        maxnode(1)=x
       endif
       if (y.lt.minnode(2)) then
        minnode(2)=y
       endif
       if (y.gt.maxnode(2)) then
        maxnode(2)=y
       endif
       if (z.lt.minnode(3)) then
        minnode(3)=z
       endif
       if (z.gt.maxnode(3)) then
        maxnode(3)=z
       endif

      end do
      do dir=1,3
       print *,"dir,min,max ",dir,minnode(dir),maxnode(dir)
      enddo
      do i=1,2
        read(2,*) garbage
      end do
      print *,"leafing elements (cells)..." 
      do i=1,Cells-1
        read(2,*) Shape
        if (Shape.eq.3) then 
         read(2,*) n1
         read(2,*) n2
         read(2,*) n3
        else if (Shape.eq.4) then
         read(2,*) n1
         read(2,*) n2
         read(2,*) n3
         read(2,*) n4
         numquads=numquads+1
        else
         print *,"Shape invalid"
         stop
        endif
        read(2,*) garbage
      end do  
      read(2,*) Shape
      if (Shape.eq.3) then
         read(2,*) n1
         read(2,*) n2
         read(2,*) n3
      else if (Shape.eq.4) then
         read(2,*) n1
         read(2,*) n2
         read(2,*) n3
         read(2,*) n4
         numquads=numquads+1
      else
         print *,"Shape invalid"
         stop
      endif
      print *,"number of cells ",Cells
      print *,"number of quads ",numquads

      close(2)
      open(unit=2, file= 'input.txt')
      read(2,*) Nodes
      read(2,*) Cells
      read(2,*) Shape

!     Nodes is number of nodes, Cells is number of cells
!     Shape: 3=Triangle, 4=Quadrilateral  (file could contain both)

      write(4,*) 'TITLE = "3D surface" '
      write(4,*) 'VARIABLES = "X", "Y", "Z" '
      if (1.eq.0) then
       write(4,*) 'ZONE T="Triangles", NODES= ', Nodes, ', ELEMENTS= ', &
        Cells+numquads, ', DATAPACKING=POINT, '
      else
       write(4,*) 'ZONE T="TRIANGLES", N= ', Nodes, ', E= ', &
        Cells+numquads, ', DATAPACKING=POINT, '
      endif

      write(4,*) 'ZONETYPE=FETRIANGLE' 

      print *,"inputting and outputting nodes..." 
      do i=1,Nodes
       read(2,*) N, x, y, z
       write(4,*)  x, y, z
      end do
      

      do i=1,2
        read(2,*) garbage
      end do

      print *,"inputting and outputting elements (cells)..." 
      do i=1,Cells-1
        read(2,*) Shape
        if (Shape.eq.3) then 
         read(2,*) n1
         read(2,*) n2
         read(2,*) n3
         write(4,*) n1, n2, n3
        else if (Shape.eq.4) then
         read(2,*) n1
         read(2,*) n2
         read(2,*) n3
         read(2,*) n4
         write(4,*) n1,n2,n3
         write(4,*) n3,n4,n1
        else
         print *,"Shape invalid"
         stop
        endif
        read(2,*) garbage
      end do  
      read(2,*) Shape
      if (Shape.eq.3) then
         read(2,*) n1
         read(2,*) n2
         read(2,*) n3
         write(4,*) n1, n2, n3
      else if (Shape.eq.4) then
         read(2,*) n1
         read(2,*) n2
         read(2,*) n3
         read(2,*) n4
         write(4,*) n1,n2,n3
         write(4,*) n3,n4,n1
      else
         print *,"Shape invalid"
         stop
      endif

      END PROGRAM
