c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
        subroutine matreadContact
     &                  (jtype,itype,l,iv,phi,jtypemax,maxordern)

        integer iorder, i2, i, nordermax
        integer itype, jtype
        real  l(jtypemax,maxordern)
        real   iv(jtypemax,maxordern,maxordern)
        real  phi(jtypemax,0:2,maxordern,12)

       write(*,*) 'CONTACT SHAPE FUNCTION', itype          
       read(8813,*) nordermax
       read(8813,*) 
! 2*n +29 lines
       do iorder=3,nordermax
         write(*,*) itype, iorder
         if (itype .eq. iorder) then
           read(8813,*) norder 
           if(norder .ne. itype) then
              write(*,*)'worng input shape'
              stop
           end if
           read(8813,*) 
           read(8813,*) 
           read(8813,*) 
           read(8813,*) l(jtype,1:itype+6)
           read(8813,*) 
           read(8813,*) 
           read(8813,*) 
           do i=1,itype+6
              read(8813,*) iv(jtype,i,1:itype+6)
           end do
           read(8813,*) 
           read(8813,*) 
           read(8813,*) 
           read(8813,*) 
           do i2=0,2
              read(8813,*) 
              read(8813,*) 
              do i=1,itype+6
                 read(8813,*) phi(jtype,i2,i,1:12)
              end do
              read(8813,*) 
              read(8813,*) 
           end do
           read(8813,*) 
           read(8813,*) 
         else
           do i=1,4*(iorder+6)+26
              read(8813,*) 
           end do
         end if
       end do

       return
       end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
