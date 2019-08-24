      PROGRAM least_squares

      IMPLICIT NONE

      integer, PARAMETER :: firstline=137
      integer, PARAMETER :: lastline=659
      Real*8 :: A(2,2)
      Real*8 :: B(2)
      Real*8 :: X(2)
      Real*8 :: tdata,cendata,det
      integer :: i,j,k

      print *,"looking for file: cen"
      print *,"linear least squares best fit:"
      print *,"firstline: ",firstline
      print *,"lastline: ",lastline

      open(unit=17, file= 'cen')

      do i=1,firstline-1
       read(17,*) tdata,cendata
      enddo
      do j=1,2
      do k=1,2
       A(j,k)=0.0
      enddo
      enddo
      do j=1,2
       B(j)=0.0
       X(j)=0.0
      enddo
      do i=1,lastline-firstline+1
       read(17,*) tdata,cendata
       if (i.eq.1) then
        print *,"first point: ",tdata,cendata
       endif
       A(1,1)=A(1,1)+1.0
       A(1,2)=A(1,2)+tdata
       A(2,1)=A(2,1)+tdata
       A(2,2)=A(2,2)+tdata*tdata
       B(1)=B(1)+cendata
       B(2)=B(2)+cendata*tdata
      enddo
      print *,"last point: ",tdata,cendata
      det=A(1,1)*A(2,2)-A(2,1)*A(1,2)
      X(1)=(A(2,2)*B(1)-A(1,2)*B(2))/det
      X(2)=(-A(2,1)*B(1)+A(1,1)*B(2))/det
      print *,"X1,X2: ",X(1),X(2)

      close(17)

      END PROGRAM
