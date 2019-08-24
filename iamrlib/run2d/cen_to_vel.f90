      PROGRAM cen_to_vel

      IMPLICIT NONE

      integer, PARAMETER :: firstline=2
      integer, PARAMETER :: lastline=4136
      Real*8 :: T1,T2,C1,C2
      Real*8 :: tdata,cendata,det,veldata
      integer :: i,j,k

      print *,"looking for file: cen"
      print *,"converting cen to vel:"
      print *,"firstline: ",firstline
      print *,"lastline: ",lastline

      open(unit=17, file= 'cen')
      open(unit=18, file= 'vel')

      do i=1,firstline-1
       read(17,*) T1,C1
      enddo
      do i=1,lastline-firstline+1
       read(17,*) T2,C2
       if (i.eq.1) then
        print *,"first point: ",T2,C2
       endif
       veldata=(C2-C1)/(T2-T1)
       write(18,*) T2,veldata
       T1=T2
       C1=C2
      enddo
      print *,"last point: ",T2,C2

      close(17)
      close(18)

      END PROGRAM
