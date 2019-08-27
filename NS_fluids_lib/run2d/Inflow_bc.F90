      subroutine inflow_bc(system,x,phi,delx)
      IMPLICIT NONE
!
      INTEGER_T   system
      REAL_T    aveQ,  aveV,  radius,  diameter,  x,  delx,  phi
      REAL_T    density,  viscosity,  sigma
      REAL_T Weber,   Reynolds,   Froude

#include "probdataf95.H" 

      if (probtype.ne.25) then
       print *,"probtype invalid in inflow_bc"
       stop
      endif
 
        radius=0.05
        diameter=2.0*radius

      go to (1,2,3,4,5,6,7,8,9,10), system
!
    1   aveV = 10
        go to 20
    2   aveV = 20
        go to 20
    3   aveV = 30
        go to 20
    4  aveV = 40
        go to 20
    5   aveV = 50
        go to 20
    6  aveV = 60
        go to 20
    7   aveV = 70
        go to 20
    8  aveV = 80
        go to 20
    9  aveV = 90
        go to 20
   10  aveV = 100
        go to 20

   20 continue

         density=0.9
         viscosity=10.0
         sigma=25.0

         Weber=(aveV**2)*radius*density/sigma
         Reynolds=density*radius*aveV/viscosity
         Froude=(aveV**2)/(radius*980.0)


         if (1.eq.0) then
          print *,"Weber,Reynolds,Froude ",Weber,Reynolds,Froude
          print *,"1/Weber,1/Reynolds,1/Froude ",one/Weber, &
             one/Reynolds,one/Froude
          stop
         endif
         if (radblob.ne.one) then
          print *,"dimensionless radius=1"
          stop
         endif
         if (fort_denconst(1).ne.one) then
          print *,"dimensionless density must be 1"
          stop
         endif
         phi=0.0
         aveV=one
         radius=one
         if (abs(x).le.radius) then
          phi=2.0*aveV*(1.0-(abs(x)/radius)**2-((delx/radius)**2)/4.0)       
         endif
!
      return
      end
