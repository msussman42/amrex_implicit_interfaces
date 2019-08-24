      PROGRAM tecpath


      IMPLICIT NONE

      Integer :: Node_start,Node_end, N, i, garbage
      Integer :: Node_len
      Integer :: n1, n2, n3
      Real :: X, Y, Z,time,base,offset,vel,zhi,endtime

      print *,"This program reads cenx,ceny,cenz"
      print *,"and outputs a tecplot file: output.dat"
      open(unit=2, file= 'cenx')
      open(unit=3, file= 'ceny')
      open(unit=4, file= 'cenz')
      open(unit=5, file= 'output.dat')

      Node_start=1  
      Node_end=45217  ! coarse
      Node_end=37568  ! fine
      Node_len=Node_end-Node_start+1
      print *,"Node_start=",Node_start
      print *,"Node_end=",Node_end
      print *,"Node_len=",Node_len
     

      write(5,*) 'TITLE = "3D moments" '
      write(5,*) 'VARIABLES = "X", "Y", "Z" '
      write(5,*) 'ZONE F="POINT", I= ', Node_len, ', J=1, K=1, ', &
        'SOLUTIONTIME= 0.0 STRANDID= 1'


      print *,"inputting and outputting the moments..." 
      vel=0.03
      zhi=0.5
      endtime=111.8  ! coarse
!      endtime=93.9  ! fine

      do i=1,Node_len
       read(2,*) time,X
       read(3,*) time,Y
       read(4,*) time,Z
       Z=time*zhi/endtime

       write(5,*)  X, Y, Z

      end do
      

      close(2)
      close(3)
      close(4)
      close(5)

      END PROGRAM
