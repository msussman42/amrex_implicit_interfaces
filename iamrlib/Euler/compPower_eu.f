      function compPower(comp,power)
      implicit none     

      include 'parameters.h'

c     Input      
      real*8 comp,power
c     Output
      real*8 compPower    

      if (comp .gt. small) then
         compPower = comp**power
      elseif (comp .ge. -small) then
         compPower = 0.d0
      else
         stop 'compPower: comp < -small'
      endif

      end
