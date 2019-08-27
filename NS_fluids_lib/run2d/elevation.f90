
      subroutine get_left_elevation(t,elevation_left)

      IMPLICIT NONE
      INTEGER,PARAMETER :: DP=SELECTED_REAL_KIND(15)
      REAL(KIND=DP),INTENT(IN) :: t
      REAL(KIND=DP),DIMENSION(2297) :: x,y 
      REAL(KIND=DP),DIMENSION(2297*3+1) :: m
      INTEGER :: n,i,k1,k2,k   
      REAL(KIND=DP),INTENT(OUT) :: elevation_left
      
      open(unit=2,file='InflowBC.dat',action='read')
      read(2,*,end=33) m
      33 close(2)
        x(1)=m(2)
        y(1)=m(3)

      do i=2,2297
        x(i)=m(i*3-1)
        y(i)=m(i*3)
      end do

     n=2297
   !first locate the interpolation interval, then linear interpolate 
      k1=1
      k2=n
     do while (k2-k1>1)
       k=floor((real(k2)+real(k1))/real(2))
        if (t>x(k)) then
          k1=k
        else
          k2=k
        end if
     end do
     elevation_left=y(k1)+((y(k2)-y(k1))/(x(k2)-x(k1)))*(t-x(k1))

      end subroutine get_left_elevation


      subroutine get_right_elevation(t,elevation_right)

      IMPLICIT NONE
      INTEGER,PARAMETER :: DP=SELECTED_REAL_KIND(15)
      INTEGER :: n,i,k1,k2,k   
      REAL(KIND=DP),INTENT(IN) :: t
      REAL(KIND=DP),DIMENSION(2297) :: x,y 
      REAL(KIND=DP),DIMENSION(2297*3+1) :: m
      REAL(KIND=DP),INTENT(OUT) :: elevation_right
      
      open(unit=2,file='OutflowBC.dat',action='read')
      read(2,*,end=33) m
      33 close(2)
        x(1)=m(2)
        y(1)=m(3)

      do i=2,2297
        x(i)=m(i*3-1)
        y(i)=m(i*3)
      end do

     n=2297
   !first locate the interpolation interval, then linear interpolate 
      k1=1
      k2=n
     do while (k2-k1>1)
       k=floor((real(k2)+real(k1))/real(2))
        if (t>x(k)) then
          k1=k
        else
          k2=k
        end if
     end do
     elevation_right=y(k1)+((y(k2)-y(k1))/(x(k2)-x(k1)))*(t-x(k1))

      end subroutine get_right_elevation

