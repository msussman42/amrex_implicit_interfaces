      subroutine overwrite(data_in)
      IMPLICIT NONE

      real*8, dimension(-2:2,-2:3) :: data_in

      data_in(-2,3)=5.0

      return
      end
 
      program main
      IMPLICIT NONE

      real*8, dimension(:,:), allocatable :: datatest
      integer sizearr

      sizearr=2

      allocate(datatest(-sizearr:sizearr,-sizearr:sizearr))
      call overwrite(datatest)
      print *,"datatest= ",datatest(-2,-2)

      
      return
      end

