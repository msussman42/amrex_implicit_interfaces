      PROGRAM wavespeed


      IMPLICIT NONE

      Real :: pi,g,h,lambda,k
      Real :: omega,c
      Real :: mu,endtime,growth

      pi=4.0*atan(1.0)
      g=980.0
      h=50.0
      lambda=25.0
      k=2.0*pi/lambda
      omega=sqrt(g*k*tanh(k*h))
      c=omega/k
      
      mu=0.01
      endtime=20.0
      growth=exp(0.5*mu*omega*endtime)
      
      print *,"g,h,lambda,k ",g,h,lambda,k
      print *,"omega, c ",omega,c
      print *,"distance travelled per second ",c
      print *,"distance travelled per 1/10 second ",c/10.0
      print *,"mu,endtime ",mu,endtime
      print *,"growth ",growth

      growth=2.5
      mu=log(growth)*2.0/(endtime*omega)
      print *,"based on growth=",growth
      print *,"mu is ",mu

      END PROGRAM
