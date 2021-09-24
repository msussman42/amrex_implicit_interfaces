      PROGRAM least_squares

      IMPLICIT NONE

      Integer :: N,i
      Real :: sumx,sumy,sumx2,sumxy,x,y,xraw,yraw,A,p

      print *,"This program reads input_data"
      print *,"and fits data with form y=A x^p"
      open(unit=2, file= 'input_data')

      print *,"reading number of points..."
      read(2,*) N
      print *,"N=",N
      print *,"reading the data..."
      sumx=0.0d0
      sumy=0.0d0
      sumx2=0.0d0
      sumxy=0.0d0
      do i=1,N
       read(2,*)  xraw, yraw
       x=log(xraw)
       y=log(yraw)
       print *,"x,y,xraw,yraw ",x,y,xraw,yraw
       sumx=sumx+x
       sumy=sumy+y
       sumx2=sumx2+x**2
       sumxy=sumxy+x*y
      enddo
       !y=A x^p
       !log y=logA + p log x
       ! alpha=logA
       ! A=e^alpha
      p=(N*sumxy-sumx*sumy)/(N*sumx2-sumx**2)
      A=(sumy*sumx2-sumx*sumxy)/(N*sumx2-sumx**2)
      A=exp(A)
      print *,"A,p ",A,p

      END PROGRAM
