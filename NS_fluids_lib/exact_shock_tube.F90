#define NCELL 400
#define probtype (2)
#define Pi (3.14159265358979d0)
#define gamma (1.4d0)
#define N_STATES 4
#define N_WAVES 3
#define N_EDGES 6

      subroutine xt_riemann( &
        rhol,vell,prel,soundl, &
        rhor,velr,prer,soundr, &
        rhoexact,velexact,preexact,x,t, &
        problo,probhi,probcen, &
        rho_0,vel_0,pre_0, &
        snd_0,w_0,signif_0,wave_type_0,IsVacuum)
      IMPLICIT NONE

      real*8 x,t,problo,probhi,probcen
      real*8 rhol,vell,prel
      real*8 rhor,velr,prer
      real*8 rhoexact,velexact,preexact
      real*8 soundl,soundr,soundexact
      real*8 rho_0(N_STATES),vel_0(N_STATES),pre_0(N_STATES)
      real*8 snd_0(N_STATES)
      real*8 w_0(N_EDGES)
      real*8 xfront(N_EDGES)
      real*8 wexpand,xmid
      integer signif_0(N_WAVES),wave_type_0(N_WAVES)
      integer nfronts
      logical IsVacuum

      if (rhol.le.0.0) then
       print *,"rhol cannot be negative riemann_flux: ",rhol
       stop
      endif
      if (prel.le.0.0) then
       print *,"prel cannot be negative riemann_flux: ",prel
       stop
      endif
      if (rhor.le.0.0) then
       print *,"rhor cannot be negative riemann_flux: ",rhor
       stop
      endif
      if (prer.le.0.0) then
       print *,"prer cannot be negative riemann_flux: ",prer
       stop
      endif
      
      xmid=probcen
      do nfronts=1,N_EDGES
       xfront(nfronts)=xmid+w_0(nfronts)*t
      enddo

      if (x.le.xfront(1)) then
       rhoexact=rhol
       velexact=vell
       preexact=prel
      else if (x.le.xfront(2)) then
       wexpand=(x-xmid)/t
       soundexact=(vell+2.0*soundl/(gamma-1.0)-wexpand)/ &
                  (1.0+2.0/(gamma-1.0))
       velexact=soundexact+wexpand
       rhoexact=(soundexact**2)*(rhol**gamma)/(gamma*prel)
       rhoexact=rhoexact**(1.0/(gamma-1.0))
       preexact=prel*(rhoexact**gamma)/(rhol**gamma)
      else if ((x.le.xfront(3)).or.(x.le.xfront(4))) then
       rhoexact=rho_0(2) 
       velexact=vel_0(2) 
       preexact=pre_0(2)
      else if (x.le.xfront(5)) then
       rhoexact=rho_0(3) 
       velexact=vel_0(3) 
       preexact=pre_0(3)
      else if (x.le.xfront(6)) then
       soundexact=(-velr+2.0*soundr/(gamma-1.0)+wexpand)/ &
                  (1.0+2.0/(gamma-1.0))
       velexact=wexpand-soundexact
       rhoexact=(soundexact**2)*(rhor**gamma)/(gamma*prer)
       rhoexact=rhoexact**(1.0/(gamma-1.0))
       preexact=prer*(rhoexact**gamma)/(rhor**gamma)
      else
       rhoexact=rhor
       velexact=velr
       preexact=prer
      endif

      return
      end

      subroutine density_bc(xlo,dx,lo,hi,den)
      IMPLICIT NONE

      integer lo,hi
      real*8 xlo,dx           
      real*8 den(lo-1:hi+1,3)
      integer i,j

      i=lo-1
      do j=1,3
       den(i,j)=den(lo,j)
      enddo
      i=hi+1
      do j=1,3
       den(i,j)=den(hi,j)
      enddo

      return
      end

      subroutine output_density(xlo,dx,lo,hi,state,nframe)
      IMPLICIT NONE

      integer lo,hi,nframe
      real*8 xlo,dx
      real*8 state(lo-1:hi+1,3)
      real*8 denfab(lo-1:hi+1,3)
      integer i,j
      real*8 x_output
      character*3 namestr
      character*5 framestr
      character*8 filename
      real*8 den,vel,pres,mom,energy,internal_energy,KE

      namestr='EXA'
      write(framestr,300) nframe
      do i=1,5
       if (framestr(i:i).eq.' ') then
        framestr(i:i)='0'
       endif
      enddo
      write(filename,400) namestr,framestr
      print *,"filename ",filename
      open(unit=11,file=filename)

      do i=lo,hi
       do j=1,3
        denfab(i,j)=state(i,j)
       enddo
      enddo
      call density_bc(xlo,dx,lo,hi,denfab)
      do i=lo,hi
       x_output=xlo+dx*(i-lo+0.5)
       den=denfab(i,1)
       mom=denfab(i,2)
       energy=denfab(i,3)
       vel=mom/den
       internal_energy=energy/den-0.5*vel*vel
       pres=den*internal_energy*(gamma-1.0) 
       write(11,*) x_output,den,vel,pres,energy
      enddo

      close(11)
300   FORMAT(I5)
400   FORMAT(A3,A5)

      return
      end
     

      program main
      IMPLICIT NONE
      real*8 dx,xlo,xhi
      integer lo,hi
      real*8 den(-1:NCELL,3)
      real*8 time,x
      integer i,j,nstep
      real*8 denl,vell,presl,denr,velr,presr
      real*8 soundl,soundr
      real*8 problo,probhi,probcen
      real*8 den_exact,vel_exact,p_exact

      real*8 rho_0(N_STATES),vel_0(N_STATES),pre_0(N_STATES)
      real*8 snd_0(N_STATES)
      real*8 w_0(N_EDGES)
      integer signif_0(N_WAVES),wave_type_0(N_WAVES)
      logical IsVacuum

      if (probtype.eq.0) then
       denl=1.0
       vell=0.0
       presl=1.0
       denr=1.0/8.0
       velr=0.0
       presr=0.1
       time=0.15
       problo=0.0
       probhi=1.0
       probcen=0.5*(problo+probhi)
      else if (probtype.eq.1) then  ! mach 4 shock tube
       denl=10.0
       vell=5.0
       presl=10.0*(gamma-1.0)
       denr=1.0
       velr=0.0
       presr=1.0*(gamma-1.0)
       time=0.1
       problo=0.0
       probhi=2.0
       probcen=1.0
      else if (probtype.eq.2) then ! strong shock tube
       denl=1.0
       vell=0.0
       presl=1.0D+10
       denr=0.125
       velr=0.0
       presr=0.1
       time=2.5E-6
       problo=0.0
       probhi=1.0
       probcen=0.5
      else
       print *,"probtype invalid"
       stop
      endif
      
      soundl=sqrt(gamma*presl/denl)
      soundr=sqrt(gamma*presr/denr)

      print *,"LEFT: den,u,p,sound ",denl,vell,presl,soundl
      print *,"RIGHT: den,u,p,sound ",denr,velr,presr,soundr
      print *,"problo,probhi,time ",problo,probhi,time

      call riemann( &
       denl,vell,presl,soundl, &
       denr,velr,presr,soundr, &
       rho_0,vel_0,pre_0,snd_0, w_0,signif_0, &
       wave_type_0,IsVacuum)
 
      lo=0
      hi=NCELL-1
      xlo=problo
      dx=(probhi-problo)/NCELL

      do i=lo,hi
       x=xlo+(i-lo+0.5d0)*dx
       call xt_riemann( &
        denl,vell,presl,soundl, &
        denr,velr,presr,soundr, &
        den_exact,vel_exact,p_exact,x,time, &
        problo,probhi,probcen, &
        rho_0,vel_0,pre_0, &
        snd_0,w_0,signif_0,wave_type_0,IsVacuum)
       den(i,1)=den_exact
       den(i,2)=den_exact*vel_exact   
       den(i,3)=p_exact/(gamma-1.0)+0.5*den_exact*(vel_exact**2)
      enddo

      nstep=0
      print *,"output of x,den,u,p,KE"
      call output_density(xlo,dx,lo,hi,den,nstep) 

      return
      end

