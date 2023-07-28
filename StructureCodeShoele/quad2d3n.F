c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine quad2d3n(ngau,gau,wei,iflag)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c calculate gauss quadrature points 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'grid_def'
      integer  ngau,iigau,iflag

      real(8)  gau(2,ngaumax), wei(ngaumax)
      if(iflag .eq. 0) then !main gauss quad
           iigau=igau
      elseif(iflag .eq. -1) then !main gauss quad
           iigau=10
      else  !for contact calculations
           iigau=igauCont
      endif

      if(iigau .le. 9) then      
      go to (1,2),iigau
 1    continue
      ngau = 1
      gau(1,1)=0.333333333333333
      gau(2,1)=0.333333333333333
      wei(  1)=0.5
      return
 2    continue
      ngau = 3
      gau(1,1)= 0.166666666666667
      gau(2,1)= 0.166666666666667
      gau(1,2)= 0.166666666666667
      gau(2,2)= 0.666666666666667
      gau(1,3)= 0.666666666666667
      gau(2,3)= 0.166666666666667
      wei(1)  = 0.166666666666667
      wei(2)  = 0.166666666666667
      wei(3)  = 0.166666666666667
      return
      else
      ngau = 3

      gau(1,1)= 0.00000000000000
      gau(2,1)= 0.00000000000000
      gau(1,2)= 1.00000000000000
      gau(2,2)= 0.00000000000000
      gau(1,3)= 0.00000000000000
      gau(2,3)= 1.00000000000000
      wei(1)  = 1.00000000000000
      wei(2)  = 1.00000000000000
      wei(3)  = 1.00000000000000
      return
      endif
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
