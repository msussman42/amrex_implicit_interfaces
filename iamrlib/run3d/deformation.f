
! grep "MAT=1 mindistcen=" run.out > mindist
! grep "MAT=1 maxdistcen=" run.out > maxdist

      program main
      IMPLICIT NONE

      integer ntime,i
      real*8 time,time2
      real*8 mindistcen,maxdistcen
      real*8 deformation

      print *,"this program reads 'mindist','maxdist'"
      print *,"and outputs 'deformation'"

      ntime=5090
      print *,"ntime is the number of entries in the file"
      print *,"ntime=",ntime
 
      open(unit=12, file='mindist')
      open(unit=13, file='maxdist')
      open(unit=14, file='deformation')

      do i=1,ntime
       read(12,*) time,mindistcen
       read(13,*) time2,maxdistcen
       if (abs(time-time2).gt.1.0E-10) then
        print *,"time discrepancy"
        stop
       endif
       deformation=(maxdistcen-mindistcen)/(maxdistcen+mindistcen)
       write(14,*) time,deformation
      enddo

      close(12)
      close(13)
      close(14)

      end

