      subroutine polyinterpdata(tdata,cendata,ndata,x,cen_interp)
      IMPLICIT NONE

      integer, INTENT(in) :: ndata
      Real*8, INTENT(in) :: tdata(ndata)
      Real*8, INTENT(in) :: cendata(ndata)
      Real*8, INTENT(in) :: x
      Real*8, INTENT(out) :: cen_interp
      integer :: i,j
      Real*8 prod

      cen_interp=0.0d0
      do i=1,ndata
       prod=1.0d0
       do j=1,ndata
        if (i.ne.j) then
         prod=prod*(x-tdata(j))/(tdata(i)-tdata(j))
        endif
       enddo
       cen_interp=cen_interp+cendata(i)*prod
      enddo

      end subroutine polyinterpdata

      PROGRAM poly_interp

      IMPLICIT NONE

      integer, PARAMETER :: firstline=1
      integer, PARAMETER :: lastline=4
      integer, PARAMETER :: nplot=100
      Real*8 :: tdata_arr(lastline-firstline+1)
      Real*8 :: cendata_arr(lastline-firstline+1)
      Real*8 :: tdata,cendata,xi,cen_interp
      integer :: i,ndata

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

      do i=0,nplot
       xi=tdata_arr(1)+i*(tdata_arr(ndata)-tdata_arr(1))/nplot
       call polyinterpdata(tdata_arr,cendata_arr,ndata,xi,cen_interp)
       print *,"xi,cen_interp ",xi,cen_interp
      enddo

      END PROGRAM
