      program main
      IMPLICIT NONE

       ! fit: d=a t^p + b
       ! d1=a f(t1) + b
       ! d2=a f(t2) + b
       ! a=(d2-d1)/(f(t2)-f(t1)) 
      real*8 a,b
      real*8 t1,t2,d1,d2,p

      t1=11.0
      d1=0.026
      t2=40.0
      d2=0.046
      p=1.0/3.0
      a=(d2-d1)/(t2**p - t1**p)
      b=d1-a*(t1**p)
      print *,"t1,d1 ",t1,d1
      print *,"t2,d2 ",t2,d2
      print *,"d=a t^p + b"
      print *,"a,b,p ",a,b,p
      return
      end

