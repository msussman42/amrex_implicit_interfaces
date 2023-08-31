c  August 31, 2023 grid_def.F changed
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
        subroutine matread(jtype,itype,l,iv,phi)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c read required data 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
       include 'grid_def'

       integer iorder, i2, i, nordermax
       integer itype, jtype
       real  l(jtypemax,maxordern)
       real   iv(jtypemax,maxordern,maxordern)
       real  phi(jtypemax,0:2,maxordern,12)

       read(13,*) nordermax
       read(13,*) 
! 2*n +29 lines
       do iorder=3,nordermax
         if (itype .eq. iorder) then
           read(13,*) norder 
           if(norder .ne. itype) then
              write(*,*)'worng input shape'
              stop
           end if
           read(13,*) 
           read(13,*) 
           read(13,*) 
           read(13,*) l(jtype,1:itype+6)
           read(13,*) 
           read(13,*) 
           read(13,*) 
           do i=1,itype+6
              read(13,*) iv(jtype,i,1:itype+6)
           end do
           read(13,*) 
           read(13,*) 
           read(13,*) 
           read(13,*) 
           do i2=0,2
              read(13,*) 
              read(13,*) 
              do i=1,itype+6
                 read(13,*) phi(jtype,i2,i,1:12)
              end do
              read(13,*) 
              read(13,*) 
           end do
           read(13,*) 
           read(13,*) 
         else
           do i=1,4*(iorder+6)+26
              read(13,*) 
           end do
         end if
       end do

       return
       end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
