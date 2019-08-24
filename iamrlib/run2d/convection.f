      program main
      IMPLICIT NONE
      real*8 k,DeltaT,mu,den,grav,beta,alpha,cp
      real*8 cv,latentheat,tension,prandtl,reynolds,uscale
      real*8 tscale,lmicro,twait,vaporcavitysize
      real*8 denboil,densupersix,denbubble
      real*8 thickness,lscale,delpower,tempvar

      lmicro=1.0E-7
      vaporcavitysize=0.03
      k=0.67888E+5
      cp=4.21850E+7
      cv=3.77472E+7
      latentheat=2257.0E+7
      DeltaT=6.2
      mu=0.0089
      denbubble=0.00063
      den=0.95881
      alpha=k/(cp*den)
      grav=980.0
      denboil=0.9583665
      densupersix=0.953989
      beta=(1.0-densupersix/denboil)/6.0
      thickness=7.14*( (mu*alpha/(grav*beta*DeltaT*den))**(1.0/3.0) )
      tension=58.0
      lscale=sqrt(tension/(grav*den))
      uscale=sqrt(grav*lscale)
      tscale=lscale/uscale
      print *,"k=",k
      print *,"alpha=",alpha
      print *,"DeltaT=",DeltaT
      print *,"mu=",mu
      print *,"den=",den
      print *,"grav=",grav
      print *,"beta=",beta
      print *,"lscale=",lscale
      print *,"uscale=",uscale
      print *,"tscale=",tscale
      print *,"following Dhir"
      print *,"thickness=",thickness
      print *,"thickness/lscale=",thickness/lscale

      thickness=(grav*beta*DeltaT*den/(mu*alpha))**0.24
      thickness=thickness*0.0014
      delpower=1.0-3.0*0.24
      thickness=thickness**(1.0/delpower)
      print *,"following Kozanoglu and Lopez"
      print *,"thickness=",thickness
      print *,"thickness/lscale=",thickness/lscale

      print *,"units of erg: dyne * cm, 1J=10^7 erg"
      print *,"units of cp: erg/(g K)"
      print *,"units of latent heat: erg/g"
      print *,"units of k: erg/(cm s K)"
      print *,"units of alpha: cm^2/s"

      reynolds=den*lscale*uscale/mu
      prandtl=cp*mu/k
      print *,"Reynolds number, 1/Reynolds : ",reynolds,1.0/reynolds
      print *,"Prandtl number, 1/(Pr Re)  : ",prandtl,
     &   1.0/(prandtl*reynolds)
      print *,"scaled latent heat: ",latentheat/(cp*DeltaT)
      print *,"scaled thermal expansion coefficient: ",
     &   beta*DeltaT 
      print *,"scaled microscale length: ",lmicro/lscale
      tempvar=2.0*tension/(vaporcavitysize*den*latentheat)
      twait=(vaporcavitysize/(1.0-tempvar))**2
      twait=twait*2.0/(4.0*3.14*alpha)
      print *,"waiting time, scaled waiting time:",twait,
     &  twait/tscale
      
      end

