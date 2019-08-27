      PROGRAM atantest

      IMPLICIT NONE

      Real*8 :: X, Y, angle
      integer :: i

      X=1.0E-20
      Y=1.0
      angle=atan(Y/X)
      print *,"X,Y,angle ",X,Y,angle

      do i=0,1000
       X=(i-500.0)/10.0
       print *,"X,angle ",X,atan(X)
      enddo

      END PROGRAM
