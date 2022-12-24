      PROGRAM least_squares

      IMPLICIT NONE

      integer, PARAMETER :: firstline=1
      integer, PARAMETER :: lastline=4
      Real*8 :: AA
      Real*8 :: A(2,2)
      Real*8 :: B(2)
      Real*8 :: X(2)
      Real*8 :: tdata,cendata,det
      Real*8 :: tdata_arr(lastline-firstline+1)
      Real*8 :: cendata_arr(lastline-firstline+1)
      REAL*8 :: LS_ERR,predict
      REAL*8 :: LS_ERR_raw,predict_raw
      REAL*8 :: sigma_1,sigma_2
      REAL*8 :: sd_1,sd_2
      integer :: i,j,k
      integer :: fit_type
      integer :: ndata

      ndata=lastline-firstline+1

      print *,"looking for file: cen"
      print *,"firstline: ",firstline
      print *,"lastline: ",lastline
      print *,"ndata: ",ndata

      open(unit=17, file= 'cen')

      do i=1,firstline-1
       read(17,*) tdata,cendata
      enddo
      do i=1,ndata
       read(17,*) tdata,cendata
       if (i.eq.1) then
        print *,"first point: ",tdata,cendata
       endif
       tdata_arr(i)=tdata
       cendata_arr(i)=cendata
      enddo
      print *,"last point: ",tdata,cendata
      close(17)

      do fit_type=0,2

       do j=1,2
       do k=1,2
        A(j,k)=0.0
       enddo
       enddo
       do j=1,2
        B(j)=0.0
        X(j)=0.0
       enddo
       do i=1,ndata
        tdata=tdata_arr(i)
        cendata=cendata_arr(i)
        if (fit_type.eq.0) then
         ! do nothing y=a+bx
        else if (fit_type.eq.1) then 
         ! log y = log a + b log x  (y=a x^b)
         tdata=log(tdata)
         cendata=log(cendata)
        else if (fit_type.eq.2) then
         ! log y= log a + b x    y=a e^(bx)
         cendata=log(cendata)
        else
         print *,"fit_type invalid"
         stop
        endif

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
       print *,"fit_type= ",fit_type
       LS_ERR=0.0
       LS_ERR_raw=0.0
       do i=1,ndata
        predict_raw=X(1)+X(2)*tdata
        if (fit_type.eq.0) then
         AA=X(1)
         predict=AA+X(2)*tdata_arr(i)
        else if (fit_type.eq.1) then
         AA=exp(X(1))
         predict=AA*(tdata_arr(i)**X(2))
        else if (fit_type.eq.2) then
         AA=exp(X(1))
         predict=AA*exp(X(2)*tdata_arr(i))
        endif
        LS_ERR=LS_ERR+(predict-cendata_arr(i))**2
        LS_ERR_raw=LS_ERR_raw+(predict_raw-cendata)**2
       enddo
       LS_ERR=sqrt(LS_ERR)

       print *,"fit_type",fit_type
       print *,"AA,X2: ",AA,X(2)
       print *,"LS_ERR= ",LS_ERR
       print *,"LS_ERR_raw= ",LS_ERR_raw
       print *,"X(1),X(2) ",X(1),X(2)
       sigma_1=(LS_ERR_raw/(ndata-2.0d0))*A(2,2)/det
       sd_1=sqrt(sigma_1)
       sigma_2=(LS_ERR_raw/(ndata-2.0d0))*A(1,1)/det
       sd_2=sqrt(sigma_2)
       print *,"variance,standard dev of X(1) ",sigma_1,sd_1
       print *,"variance,standard dev of X(2) ",sigma_2,sd_2
       if ((fit_type.eq.1).or.(fit_type.eq.2)) then
        print *,"estimated variance of e^X(1): ", &
          max((exp(X(1)+sd_1)-AA)**2,(exp(X(1)-sd_1)-AA)**2)
        print *,"estimated standard dev of e^X(1): ", &
          max(abs(exp(X(1)+sd_1)-AA),abs(exp(X(1)-sd_1)-AA))
       endif
      enddo

      END PROGRAM
