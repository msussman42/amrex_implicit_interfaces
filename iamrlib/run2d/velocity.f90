      subroutine get_left_velocity(t,u_left,v_left)

      IMPLICIT NONE
      INTEGER,PARAMETER :: DP=SELECTED_REAL_KIND(15)
      INTEGER :: n,i,k1,k2,k  
      REAL(KIND=DP),INTENT(IN) :: t
      REAL(KIND=DP),DIMENSION(2297) :: x,y,z
      REAL(KIND=DP),DIMENSION(2297*3+1) :: m
      REAL(KIND=DP),INTENT(OUT) :: u_left,v_left
      
      open(unit=2,file='InflowBC.dat',action='read')
      read(2,*,end=33) m
      33 close(2)
      n=2297
      x(1)=m(2) 
      y(1)=m(3)
      z(1)=m(4) 
      do i=2,n
        x(i)=m(i*3-1)
        y(i)=m(i*3)
        z(i)=m(i*3+1)
      end do

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
      u_left=z(k1)+((z(k2)-z(k1))/(x(k2)-x(k1)))*(t-x(k1))
      v_left=y(k1)+((y(k2)-y(k1))/(x(k2)-x(k1)))*(t-x(k1))
      end subroutine get_left_velocity

      subroutine get_right_velocity(t,u_right,v_right)

      IMPLICIT NONE
      INTEGER,PARAMETER :: DP=SELECTED_REAL_KIND(15)
      INTEGER :: n,i,k1,k2,k   
      REAL(KIND=DP),INTENT(IN) :: t
      REAL(KIND=DP),DIMENSION(2297) :: x,y,z
      REAL(KIND=DP),DIMENSION(2297*3+1) :: m
      REAL(KIND=DP),INTENT(OUT) :: u_right,v_right
      
      open(unit=2,file='OutflowBC.dat',action='read')
      read(2,*,end=33) m
      33 close(2)
      n=2297
      x(1)=m(2) 
      y(1)=m(3)
      z(1)=m(4) 
      do i=2,n
        x(i)=m(i*3-1)
        y(i)=m(i*3)
        z(i)=m(i*3+1)
      end do

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
      u_right=z(k1)+((z(k2)-z(k1))/(x(k2)-x(k1)))*(t-x(k1))
      v_right=y(k1)+((y(k2)-y(k1))/(x(k2)-x(k1)))*(t-x(k1))
     end subroutine get_right_velocity
