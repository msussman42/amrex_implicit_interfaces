c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine adjustsub(das)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c this is temp subroutine and is not called by other codes (FOR TESTS)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      USE HeaderFSI
      

      integer i,j,i2
      real*8 das

      das=time2
      goto 121
      do i=1,nr_ibm_esh
       do j=1,ns_ibm_r_esh(i)
          das=0.d0
             do i2=1,nposele_esh(i,j,1)
                 das=das+area(i,nposele_esh(i,j,1+i2))
             end do
             ds_ibm_esh(i,j)=das/dble(nposele_esh(i,j,1))
       end do
      enddo
      call  BucketSorting
      timebucket=0.0
 121  continue

      return
      end subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
