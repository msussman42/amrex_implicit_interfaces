PROGRAM deleterad
IMPLICIT NONE

character(20) :: dwave,dwavenew
integer :: numlines,i,j
real :: time,rad,timelast

 dwave="rad.txt"
 dwavenew="radnew.txt"
 print *,"opening ",dwave
 OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')
 print *,"opening ",dwavenew
 OPEN(unit=15,file=dwavenew,access='sequential',form="formatted",status='new')
 READ(14,*) numlines
 print *,"numlines=",numlines
 j=0
 timelast=0.0
 do i=1,numlines
  READ(14,*) time,rad
  if ((j.eq.0).or.(timelast.ne.time)) then
   WRITE(15,*) time,rad
   j=j+1
  endif
  timelast=time
 enddo
 close(14)
 close(15)

end PROGRAM deleterad

