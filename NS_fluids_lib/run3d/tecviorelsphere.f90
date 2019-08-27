      PROGRAM tecwhale


      IMPLICIT NONE

      Integer :: Nodes, Cells, Shape, N, i, garbage
      Integer :: n1, n2, n3,n4
      Real :: X, Y, Z
      Integer :: num_delete,num_delete_end,num_delete_elements

      num_delete=48
      num_delete_end=4
      num_delete_elements=0

      print *,"This program reads input.txt"
      print *,"and outputs a tecplot file: output.plt"
      open(unit=2, file= 'input.txt')
      open(unit=4, file= 'output.plt')

      print *,"reading number of nodes and number of faces ..."
      read(2,*) Nodes,Cells
      print *,"number of nodes=",Nodes
      print *,"number of cells=",Cells
      read(2,*) garbage
      if (garbage.ne.3) then
       print *," this code only converts surfaces with triangular faces"
       stop
      endif
      Shape=3

!     Nodes is number of nodes, Cells is number of cells
!     Shape: 3=Triangle, 4=Quadrilateral

      write(4,*) 'TITLE = "3D surface" '
      write(4,*) 'VARIABLES = "X", "Y", "Z" '
      if (1.eq.0) then
       write(4,*) 'ZONE T="Triangles", NODES= ', &
         Nodes-num_delete-num_delete_end, &
          ', ELEMENTS= ', &
        Cells, ', DATAPACKING=POINT, '
      else
       write(4,*) 'ZONE T="TRIANGLES", N= ',  &
         Nodes-num_delete-num_delete_end, ', E= ', &
        Cells, ', DATAPACKING=POINT, '
      endif

      write(4,*) 'ZONETYPE=FETRIANGLE' 

      print *,"inputting and outputting nodes..." 
      do i=1,Nodes

       read(2,*) garbage,x, y, z
       if (garbage.ne.i) then
        print *,"node index invalid"
        stop
       endif
       if ((garbage.gt.num_delete).and. &
           (garbage.le.Nodes-num_delete_end))  then
        write(4,*)  x, y, z
       endif

      end do
      
      print *,"inputting and outputting elements (cells)..." 

      read(2,*) garbage
      if (garbage.ne.Cells) then
       print *,"Cells count invalid"
       stop
      endif
      do i=1,Cells
     
        read(2,*) garbage
        if (garbage.ne.i) then
         print *,"cell index invalid"
         stop
        endif
        read(2,*) garbage
        if (garbage.ne.3) then
         print *,"all elements must be triangles"
         stop
        endif 
        read(2,*) n1
        read(2,*) n2
        read(2,*) n3
        if ((n1.gt.num_delete).and.(n2.gt.num_delete).and. &
            (n3.gt.num_delete).and. &
            (n1.le.Nodes-num_delete_end).and. &
            (n2.le.Nodes-num_delete_end).and. &
            (n3.le.Nodes-num_delete_end)) then
         write(4,*) n1-num_delete, n2-num_delete, n3-num_delete
        else
         num_delete_elements=num_delete_elements+1
        endif
      end do  
      print *,"num_delete_elements= ",num_delete_elements

      END PROGRAM
