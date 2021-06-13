!get rid of autoindent   :setl noai nocin nosi inde=
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

#include "N_EXTRA_REAL.H"
#include "GODUNOV_F.H"


#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif



      function wallFunc(u_tau,u,y,K,B,rho_w,mu_w) result(f)
      implicit none
        !incompressible wall function, to be called in 
        !Newton's method, solving for friction 
        !velocity u_tau at the image point
        !u_tau: friction velocity, u: velocity parallel to wall, 
        !y: coord direction normal to wall
       REAL_T, intent(in) :: u_tau, u, y 
        !rho_w: wall density, mu_w: wall molecular viscosity 
       REAL_T, intent(in) :: K, B, rho_w, mu_w 
       REAL_T :: f
       REAL_T :: u_plus, y_plus

       if (u_tau.eq.zero) then
        print *,"u_tau.eq.zero"
        stop
       endif
       if (mu_w.eq.zero) then
        print *,"mu_w.eq.zero"
        stop
       endif
       u_plus = u/u_tau
       y_plus = rho_w*u_tau*y/mu_w

       f = -y_plus + u_plus + exp(-K*B)* &
        ( exp(K*u_plus)-one-K*u_plus-half*(K*u_plus)**2-(K*u_plus)**3/six )
      end function wallFunc
  
      function wallFuncDeriv(u_tau,u,y,K,B,rho_w,mu_w) result(fprime)
      implicit none
        !derivative of incompressible wall function w.r.t. u_tau, to be 
        !called in Newton's method
       REAL_T, intent(in) :: u_tau, u, y 
       REAL_T, intent(in) :: K, B, rho_w, mu_w
       REAL_T :: fprime
       REAL_T :: u_plus
   
       if (u_tau.eq.zero) then
        print *,"u_tau cannot be 0"
        stop
       endif 
       u_plus = u/u_tau
    
       fprime = -rho_w*y/mu_w - u_plus/u_tau + &
         exp(-K*B)*( -(K*u_plus/u_tau)*exp(K*u_plus)+ &
                     K*u_plus/u_tau+(K*u_plus)**2/u_tau+ &
                     (K*u_plus)**3/(two*u_tau) )
      end function wallFuncDeriv
      
      subroutine wallFunc_NewtonsMethod(u,y,tau_w,im_fluid)
       use probf90_module
       use global_utility_module
       use MOF_routines_module
       implicit none
       INTEGER_T, intent(in) :: im_fluid
       REAL_T, intent(in) :: u, y !uimage_tngt, delta_r
       REAL_T, intent(out) :: tau_w
       REAL_T :: u_tau,x_n,x_np1 !x_n, x_(n+1) --> u_tau
       REAL_T :: f, fprime
       REAL_T :: wallFunc, wallFuncDeriv
       INTEGER_T :: iter, iter_max=1000
        !initialize wall function parameters here
       REAL_T :: K=0.41, B=5.5, &
       rho_w, mu_w !rho_w: wall density, mu_w: wall 
                   !molecular viscosity, tau_w: wall shear stress

       if ((im_fluid.lt.1).or.(im_fluid.gt.num_materials)) then
        print *,"im_fluid invalid in wallFunc_NewtonsMethod"
        stop
       endif

       mu_w=fort_viscconst(im_fluid) 

       x_n = u !initial guess for u_tau
       do while ((abs(x_np1-x_n).gt.VOFTOL).and.(iter.lt.iter_max))
        f = wallFunc(x_n,u,y,K,B,rho_w,mu_w)
        fprime = wallFuncDeriv(x_n,u,y,K,B,rho_w,mu_w)
   
        if (abs(fprime).le.1.0e-15) then 
         print *, "divide by zero error: wallFunc_NewtonsMethod" !avoid /0
         stop
        else
         x_np1 = x_n-f/fprime
        endif
    
        iter = iter+1
        if(iter.ge.iter_max)then
          print *, "wallFunc_NewtonsMethod: no convergence"
          stop
        endif
       enddo ! while (abs(x_np1-x_n)>VOFTOL .and. iter<iter_max)
    
       u_tau = x_np1
       tau_w = rho_w*u_tau**2
      end subroutine wallFunc_NewtonsMethod


      module godunov_module
      use probf90_module

      type law_of_wall_parm_type
      INTEGER_T :: level
      INTEGER_T :: finest_level
      INTEGER_T :: bfact
      INTEGER_T :: nmat
      INTEGER_T :: nten
      REAL_T :: visc_coef
      REAL_T :: time
      REAL_T :: dt
      REAL_T, pointer :: usolid_raster(:)
      REAL_T, pointer :: n_raster(:)
      REAL_T, pointer :: x_image_raster(:)
      REAL_T, pointer :: x_projection_raster(:)
      REAL_T, pointer :: dx(:)
      REAL_T :: dxmin
      REAL_T, pointer :: xlo(:)
      INTEGER_T, pointer :: fablo(:)
      INTEGER_T, pointer :: fabhi(:)
      INTEGER_T :: ngrowfab ! should be 4
      REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: LSCP
      REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: LSFD
      REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: state ! nden comp.
      REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: ufluid
      REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: usolid
      end type law_of_wall_parm_type

      contains


      subroutine interp_face(x,xnode_array, &
        nodedata,iten,nten,interpdata,dirnormal)
      IMPLICIT NONE

      INTEGER_T dirnormal
      INTEGER_T iten,nten
      REAL_T x(SDIM)
      REAL_T xnode_array(0:1,0:1,SDIM)
      REAL_T nodedata(0:1,0:1,nten)
      REAL_T interpdata
      INTEGER_T dir,dirtan,dirloc,side
      REAL_T dx(SDIM)
      REAL_T xfrac(SDIM)
      REAL_T lenx,delta
 
      if ((dirnormal.lt.0).or.(dirnormal.gt.SDIM-1)) then
       print *,"dirnormal invalid"
       stop
      endif
 
      do dir=1,2
       dx(dir)=zero
       do side=0,1
        lenx=zero
        do dirloc=1,SDIM
         if (dir.eq.1) then
          delta=xnode_array(1,side,dirloc)-xnode_array(0,side,dirloc)
         else if (dir.eq.2) then 
          delta=xnode_array(side,1,dirloc)-xnode_array(side,0,dirloc)
         else
          print *,"dir invalid"
          stop
         endif
         if (delta.lt.zero) then
          print *,"delta invalid"
          stop
         endif
         lenx=lenx+delta**2
        enddo ! dirloc=1..sdim
        lenx=sqrt(lenx)
        dx(dir)=dx(dir)+lenx
       enddo ! side=0..1
       dx(dir)=half*dx(dir)
      enddo ! dir=1,2
        
      if (dx(1).le.zero) then
       print *,"dx(1) invalid"
       stop
      endif
      if (SDIM.eq.2) then
       if (dx(2).ne.zero) then
        print *,"dx(2) invalid"
        stop
       endif
      else if (SDIM.eq.3) then
       if (dx(2).le.zero) then
        print *,"dx(2) invalid"
        stop
       endif
      else
       print *,"dimension bust"
       stop
      endif

      do dir=1,2
       if (dirnormal.eq.0) then
        if (dir.eq.1) then
         dirtan=2
        else if ((dir.eq.2).and.(SDIM.eq.3)) then
         dirtan=3
        else if ((dir.eq.2).and.(SDIM.eq.2)) then
         dirtan=0
        else
         print *,"dir invalid"
         stop
        endif
       else if (dirnormal.eq.1) then
        if (dir.eq.1) then
         dirtan=1
        else if ((dir.eq.2).and.(SDIM.eq.3)) then
         dirtan=3
        else if ((dir.eq.2).and.(SDIM.eq.2)) then
         dirtan=0
        else
         print *,"dir invalid"
         stop
        endif
       else if ((dirnormal.eq.2).and.(SDIM.eq.3)) then
        if (dir.eq.1) then
         dirtan=1
        else if (dir.eq.2) then
         dirtan=2
        else
         print *,"dir invalid"
         stop
        endif
       else
        print *,"dirnormal invalid"
        stop
       endif

       if ((dirtan.eq.0).and.(SDIM.eq.2).and.(dir.eq.2)) then
        xfrac(dir)=xfrac(1)
       else if ((dirtan.ge.1).and.(dirtan.le.SDIM)) then
        xfrac(dir)=(x(dirtan)-xnode_array(0,0,dirtan))/dx(dir)
       else
        print *,"dirtan invalid"
        stop
       endif
       if (xfrac(dir).lt.zero) then
        xfrac(dir)=zero
       endif
       if (xfrac(dir).gt.one) then
        xfrac(dir)=one
       endif
      enddo ! dir=1,2

      if ((iten.ge.1).and.(iten.le.nten)) then
       interpdata= &
        (one-xfrac(1))*(one-xfrac(2))*nodedata(0,0,iten)+ &
        (one-xfrac(1))*xfrac(2)*nodedata(0,1,iten)+ &
        xfrac(1)*(one-xfrac(2))*nodedata(1,0,iten)+ &
        xfrac(1)*xfrac(2)*nodedata(1,1,iten)
      else
       print *,"iten invalid"
       stop
      endif

      return
      end subroutine interp_face


      subroutine interp_face_normal(x,xnode_array, &
        normaldata,iten,nten,interpnormal,dirnormal)
      IMPLICIT NONE

      INTEGER_T dirnormal
      INTEGER_T iten,nten
      REAL_T x(SDIM)
      REAL_T xnode_array(0:1,0:1,SDIM)
      REAL_T normaldata(0:1,0:1,SDIM,nten)
      REAL_T interpnormal(SDIM)
      INTEGER_T dir,dirtan,dirloc,side
      REAL_T dx(SDIM)
      REAL_T xfrac(SDIM)
      REAL_T lenx,delta
  
      if ((dirnormal.lt.0).or.(dirnormal.gt.SDIM-1)) then
       print *,"dirnormal invalid"
       stop
      endif
 
      do dir=1,2
       dx(dir)=zero
       do side=0,1
        lenx=zero
        do dirloc=1,SDIM
         if (dir.eq.1) then
          delta=xnode_array(1,side,dirloc)-xnode_array(0,side,dirloc)
         else if (dir.eq.2) then 
          delta=xnode_array(side,1,dirloc)-xnode_array(side,0,dirloc)
         else
          print *,"dir invalid"
          stop
         endif
         if (delta.lt.zero) then
          print *,"delta invalid"
          stop
         endif
         lenx=lenx+delta**2
        enddo ! dirloc=1..sdim
        lenx=sqrt(lenx)
        dx(dir)=dx(dir)+lenx
       enddo ! side=0..1
       dx(dir)=half*dx(dir)
      enddo ! dir=1,2
        
      if (dx(1).le.zero) then
       print *,"dx(1) invalid"
       stop
      endif
      if (SDIM.eq.2) then
       if (dx(2).ne.zero) then
        print *,"dx(2) invalid"
        stop
       endif
      else if (SDIM.eq.3) then
       if (dx(2).le.zero) then
        print *,"dx(2) invalid"
        stop
       endif
      else
       print *,"dimension bust"
       stop
      endif

      do dir=1,2
       if (dirnormal.eq.0) then
        if (dir.eq.1) then
         dirtan=2
        else if ((dir.eq.2).and.(SDIM.eq.3)) then
         dirtan=3
        else if ((dir.eq.2).and.(SDIM.eq.2)) then
         dirtan=0
        else
         print *,"dir invalid"
         stop
        endif
       else if (dirnormal.eq.1) then
        if (dir.eq.1) then
         dirtan=1
        else if ((dir.eq.2).and.(SDIM.eq.3)) then
         dirtan=3
        else if ((dir.eq.2).and.(SDIM.eq.2)) then
         dirtan=0
        else
         print *,"dir invalid"
         stop
        endif
       else if ((dirnormal.eq.2).and.(SDIM.eq.3)) then
        if (dir.eq.1) then
         dirtan=1
        else if (dir.eq.2) then
         dirtan=2
        else
         print *,"dir invalid"
         stop
        endif
       else
        print *,"dirnormal invalid"
        stop
       endif

       if ((dirtan.eq.0).and.(SDIM.eq.2).and.(dir.eq.2)) then
        xfrac(dir)=xfrac(1)
       else if ((dirtan.ge.1).and.(dirtan.le.SDIM)) then
        xfrac(dir)=(x(dirtan)-xnode_array(0,0,dirtan))/dx(dir)
       else
        print *,"dirtan invalid"
        stop
       endif
       if (xfrac(dir).lt.zero) then
        xfrac(dir)=zero
       endif
       if (xfrac(dir).gt.one) then
        xfrac(dir)=one
       endif
      enddo ! dir=1,2
 
      if ((iten.ge.1).and.(iten.le.nten)) then
       do dirloc=1,SDIM
        interpnormal(dirloc)= &
         (one-xfrac(1))*(one-xfrac(2))*normaldata(0,0,dirloc,iten)+ &
         (one-xfrac(1))*xfrac(2)*normaldata(0,1,dirloc,iten)+ &
         xfrac(1)*(one-xfrac(2))*normaldata(1,0,dirloc,iten)+ &
         xfrac(1)*xfrac(2)*normaldata(1,1,dirloc,iten)
       enddo ! dirloc=1..sdim
      else
       print *,"iten invalid"
       stop
      endif

      return
      end subroutine interp_face_normal



      subroutine add_crossing(ncrossing,xcrossing,iten, &
       LSint_node,xnode_array, &
       inode1,jnode1,inode2,jnode2,nten)
      IMPLICIT NONE

      INTEGER_T inode1,jnode1,inode2,jnode2
      INTEGER_T ncrossing
      REAL_T xcrossing(4,SDIM)
      INTEGER_T iten,nten
      REAL_T LSint_node(0:1,0:1,nten)
      REAL_T xnode_array(0:1,0:1,SDIM)
      REAL_T LS1,LS2,frac
      INTEGER_T dir

      if ((iten.ge.1).and.(iten.le.nten)) then
       LS1=LSint_node(inode1,jnode1,iten)
       LS2=LSint_node(inode2,jnode2,iten)
       if (((LS1.ge.zero).and.(LS2.lt.zero)).or. &
           ((LS1.lt.zero).and.(LS2.ge.zero))) then
        if ((ncrossing.ge.0).and.(ncrossing.lt.4)) then
         ncrossing=ncrossing+1
         frac=LS1/(LS1-LS2)
         do dir=1,SDIM
          xcrossing(ncrossing,dir)= &
           (one-frac)*xnode_array(inode1,jnode1,dir)+ &
           frac*xnode_array(inode2,jnode2,dir)
         enddo
        else
         print *,"ncrossing invalid"
         stop
        endif
       else if (((LS1.lt.zero).and.(LS2.lt.zero)).or. &
                ((LS1.ge.zero).and.(LS2.ge.zero))) then
        ! do nothing
       else
        print *,"LS1 or LS2 invalid"
        stop
       endif
      else
       print *,"iten invalid"
       stop
      endif

      return
      end subroutine add_crossing
 

      subroutine velfunc(vel,x,time,dx)
      IMPLICIT NONE

      REAL_T vel(SDIM)
      REAL_T x(SDIM),dx(SDIM)
      REAL_T time


      if (probtype.eq.29) then
       if (SDIM.eq.3) then
        call deform3duu(vel(1),x(1),x(2),x(SDIM),time,dx)
        call deform3dvv(vel(2),x(1),x(2),x(SDIM),time,dx)
        call deform3dww(vel(SDIM),x(1),x(2),x(SDIM),time,dx)
       else if (SDIM.eq.2) then
        call deformuu(vel(1),x(1),x(2),time,dx)
        call deformvv(vel(2),x(1),x(2),time,dx)
       else
        print *,"dimension bust"
        stop
       endif
      else if (probtype.eq.28) then
       call zalesakuu(vel(1),x(1),x(2),x(SDIM),time,dx)
       call zalesakvv(vel(2),x(1),x(2),x(SDIM),time,dx)
       if (SDIM.eq.3) then
        call zalesakww(vel(SDIM),x(1),x(2),x(SDIM),time,dx)
       endif
      else if (probtype.eq.31) then
       call circleuu(vel(1),x(1),x(2),x(SDIM))
       call circlevv(vel(2),x(1),x(2),x(SDIM))
       if (SDIM.eq.3) then
        call circleww(vel(SDIM),x(1),x(2),x(SDIM))
       endif
      else
       print *,"probtype invalid velfunc"
       stop
      endif

      return
      end subroutine velfunc

      subroutine init_passive_advect_flag(passive_advect_flag)
      IMPLICIT NONE

      INTEGER_T passive_advect_flag


      passive_advect_flag=0

      if ((probtype.eq.29).or.(probtype.eq.28).or.(probtype.eq.31)) then 
       passive_advect_flag=1
      endif

      return 
      end subroutine init_passive_advect_flag
 

      subroutine departure_node_split( &
        xstenMAC,nhalf,dx,bfact, &
        delta,passive_veltime, &
        normdir,dt,map_forward)
      use mass_transfer_module
      IMPLICIT NONE

      INTEGER_T normdir,nhalf,bfact,dir
      REAL_T xstenMAC(-nhalf:nhalf,SDIM)
      REAL_T delta
      REAL_T dx(SDIM) 
      REAL_T vel0(SDIM)
      REAL_T x0(SDIM)
      INTEGER_T map_forward
      REAL_T dt,passive_veltime,RR
      INTEGER_T passive_advect_flag

      if (nhalf.lt.1) then
       print *,"nhalf invalid departure node split"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid41"
       stop
      endif
      if ((normdir.lt.0).or.(normdir.ge.SDIM)) then
       print *,"normdir invalid"
       stop
      endif
      if ((map_forward.ne.0).and.(map_forward.ne.1)) then
       print *,"map_forward invalid"
       stop
      endif

      RR=one
      if (levelrz.eq.0) then
       ! do nothing
      else if (levelrz.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (levelrz.eq.3) then
       if (normdir.eq.1) then
        RR=xstenMAC(0,1)
       endif
      else
       print *,"levelrz invalid departure node split"
       stop
      endif

      do dir=1,SDIM
       x0(dir)=xstenMAC(0,dir)
      enddo

      call init_passive_advect_flag(passive_advect_flag)
       
      if (passive_advect_flag.eq.1) then 
       call velfunc(vel0,x0,passive_veltime,dx)
       delta=dt*vel0(normdir+1)/RR
      else if (passive_advect_flag.eq.0) then
       ! do nothing
      else
       print *,"passive_advect_flag invalid"
       stop
      endif

      if ((levelrz.eq.1).or. &
          (levelrz.eq.3)) then
       if (normdir.eq.0) then
        if (x0(normdir+1).le.VOFTOL*dx(normdir+1)) then
         delta=zero
        else
         call adjust_du(delta,normdir,x0(normdir+1),map_forward)
        endif
       else if ((normdir.eq.1).or.(normdir.eq.SDIM-1)) then
        if (x0(1).le.zero) then
         delta=zero
        endif
       else
        print *,"normdir invalid"
        stop
       endif
      endif  ! levelrz=1 or 3

      return
      end subroutine departure_node_split

      subroutine derive_density( &
       voldepart,voltarget,voltotal_depart, &
       override_density, &
       constant_density_all_time, &
       massdepart, &
       im,nmat, &
       density)
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: im,nmat
      INTEGER_T, intent(in) :: override_density(nmat)
      INTEGER_T, intent(in) :: constant_density_all_time(nmat)
      REAL_T, intent(in) :: voldepart,voltarget,voltotal_depart
      REAL_T, intent(in) :: massdepart
      REAL_T, intent(out) :: density
     
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im invalid22"
       stop
      endif 

       ! sanity check
      if (is_rigid(nmat,im).eq.1) then
       if (fort_material_type(im).ne.999) then
        print *,"fort_material_type(im) invalid"
        stop
       endif
      else if (is_rigid(nmat,im).eq.0) then
       ! do nothing
      else
       print *,"is_rigid(nmat,im) invalid"
       stop
      endif

      if ((is_rigid(nmat,im).eq.1).or. &
          (fort_material_type(im).eq.999)) then
       if (constant_density_all_time(im).eq.1) then
        density=fort_denconst(im)
       else
        print *,"constant_density_all_time invalid"
        stop
       endif
      else if ((voldepart.le.VOFTOL*voltotal_depart).or. &
               (voltarget.le.VOFTOL*voltotal_depart)) then
       density=fort_denconst(im)
      else if ((fort_material_type(im).ge.0).and. &
               (fort_material_type(im).le.MAX_NUM_EOS)) then

       if (fort_material_type(im).eq.0) then

         if (constant_density_all_time(im).eq.1) then
          density=fort_denconst(im)
         else if (constant_density_all_time(im).eq.0) then 
          density=massdepart/voldepart
         else
          print *,"constant_density_all_time invalid"
          stop
         endif

       else if (fort_material_type(im).gt.0) then
        if (constant_density_all_time(im).eq.0) then 
         density=massdepart/voltarget
        else
         print *,"constant_density_all_time invalid"
         stop
        endif
       else 
        print *,"fort_material_type is corrupt"
        stop
       endif

      else
       print *,"fort_material_type(im) invalid"
       stop
      endif

      return
      end subroutine derive_density



        ! i,j,k is upper adjoining cell
        ! i-ii,j-jj,k-kk is lower adjoining cell
        ! is,ie,js,je,ks,ke are bounds for face stencil in the
        ! TANGENTIAL (dirtan) direction.
        ! levelpc(im).
      subroutine slopecrossterm( &
        ntensor, &
        nmat,  & 
        massfrac, &
        total_mass, &
        levelpc,DIMS(levelpc), &
        faceLS,DIMS(faceLS), & ! =0 if imL<>imR or coarse/fine bc or ext. bc
        mdata,DIMS(mdata), &
        tdata,DIMS(tdata), &
        ii,jj,kk, &
        i,j,k,dir, &
        dirtan, &
        tcomp, &
        is,ie,js,je,ks,ke, &
        slopeterm)
      use global_utility_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: ntensor
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: DIMDEC(levelpc)
      INTEGER_T, intent(in) :: DIMDEC(faceLS)
      INTEGER_T, intent(in) :: DIMDEC(mdata)
      INTEGER_T, intent(in) :: DIMDEC(tdata)

      REAL_T, intent(in) :: faceLS(DIMV(faceLS),SDIM)
      REAL_T, intent(in) :: mdata(DIMV(mdata),SDIM)
      REAL_T, intent(in) :: tdata(DIMV(tdata),ntensor)
      REAL_T, intent(in) :: levelpc(DIMV(levelpc),nmat)
      REAL_T, intent(in) :: massfrac(nmat)
      REAL_T, intent(in) :: total_mass

      INTEGER_T, intent(in) :: ii,jj,kk
      INTEGER_T, intent(in) :: i,j,k,dir,dirtan
      INTEGER_T, intent(in) :: tcomp
      INTEGER_T, intent(in) :: is,ie,js,je,ks,ke
      REAL_T, intent(out) :: slopeterm

      INTEGER_T im1,jm1,km1,i1,j1,k1
      INTEGER_T im,im_primary,im_face,imL,imR
      REAL_T weight_total,sumslope
      REAL_T testslope
      INTEGER_T try_stencil
      REAL_T LSLEFT(nmat)
      REAL_T LSRIGHT(nmat)

      if (dir.eq.dirtan) then
       print *,"dir or dirtan invalid"
       stop
      endif
      if ((dir.lt.1).or.(dir.gt.SDIM)) then
       print *,"dir invalid slopecrossterm"
       stop
      endif
      if ((dirtan.lt.1).or.(dirtan.gt.SDIM)) then
       print *,"dirtan invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif

      if (ntensor.ne.SDIM*SDIM) then
       print *,"ntensor invalid"
       stop
      endif

      if (SDIM.eq.3) then
       if ((i.lt.is).or.(i.gt.ie).or. &
           (j.lt.js).or.(j.gt.je).or. &
           (k.lt.ks).or.(k.gt.ke)) then
        print *,"i,j,k invalid"
        stop
       endif
      else if (SDIM.eq.2) then
       if ((i.lt.is).or.(i.gt.ie).or. &
           (j.lt.js).or.(j.gt.je)) then
        print *,"i,j invalid"
        stop
       endif
      else
       print *,"dimension bust"
       stop
      endif

      im1=i-ii
      jm1=j-jj
      km1=k-kk
   
      if (total_mass.gt.zero) then ! not in RZ or r>0.

        ! massfrac derived from tessellating interfaces.
       im_primary=0
       do im=1,nmat
        if (im_primary.eq.0) then
         im_primary=im
        else if ((im_primary.ge.1).and.(im_primary.le.nmat)) then
         if (massfrac(im).gt.massfrac(im_primary)) then
          im_primary=im
         endif
        else
         print *,"im_primary invalid"
         stop
        endif
       enddo ! im
    
       if ((im_primary.ge.1).and.(im_primary.le.nmat).and. &
           (massfrac(im_primary).gt.zero)) then

        do im=1,nmat
         LSLEFT(im)=levelpc(D_DECL(im1,jm1,km1),im)
         LSRIGHT(im)=levelpc(D_DECL(i,j,k),im)
        enddo
        call get_primary_material(LSLEFT,nmat,imL)
        call get_primary_material(LSRIGHT,nmat,imR)
 
        if ((imL.ne.im_primary).and. &
            (imR.ne.im_primary)) then 
         ! do not increment slopeterm
        else if ((imL.eq.im_primary).or. &
                 (imR.eq.im_primary)) then

         if (mdata(D_DECL(i,j,k),dir).eq.zero) then
          ! do not increment slopeterm if both adjoining cells are 
          ! solid or the cell pair is outside the grid.
         else if (mdata(D_DECL(i,j,k),dir).eq.one) then 
          ! one of adjoining cells is fluid
          weight_total=zero
          sumslope=zero

          do i1=is,ie 
          do j1=js,je 
          do k1=ks,ke 

            !im_face=0 if, wrt i1,j1,k1, imL_1<>imR_1 or 
            !coarse/fine or exterior BC or is_clamped_face>=1
           im_face=NINT(faceLS(D_DECL(i1,j1,k1),dirtan))
           if ((im_face.ge.0).and.(im_face.le.nmat)) then
            ! do nothing
           else
            print *,"im_face invalid"
            stop
           endif

            ! the solid is represented as a rasterized boundary.
            ! if n=(1 0 0) then 2 mu D dot n=
            ! ( 2 u_x   u_y + v_x  u_z + w_x
            !   v_x+u_y   2 v_y    v_z + w_y
            !   u_z + w_x  v_z+w_y   2 w_z   ) dot (1 0 0) =
            ! (2 u_x  v_x+u_y  u_z+w_x)
           if (im_face.eq.0) then
            try_stencil=0
           else if (is_prescribed(nmat,im_face).eq.1) then 
            try_stencil=0  
           else if (is_prescribed(nmat,im_face).eq.0) then ! im_face=fluid
            if ((im_face.ge.1).and.(im_face.le.nmat)) then
             if (im_face.eq.im_primary) then
              try_stencil=1
             else if (im_face.ne.im_primary) then
              if (is_prescribed(nmat,im_primary).eq.1) then ! im_primary=solid
               try_stencil=1
              else if (is_prescribed(nmat,im_primary).eq.0) then
               try_stencil=0
              else
               print *,"is_prescribed(nmat,im_primary) invalid"
               stop
              endif
             else
              print *,"im_face invalid"
              stop
             endif
            else
             print *,"im_face invalid"
             stop
            endif
           else
            print *,"is_prescribed(nmat,im_face) invalid"
            stop
           endif
    
           if (try_stencil.eq.0) then
            ! do nothing - solid face or face not in stencil.
           else if (try_stencil.eq.1) then
    
            weight_total=weight_total+one
            testslope=tdata(D_DECL(i1,j1,k1),tcomp)
            sumslope=sumslope+testslope

           else
            print *,"try_stencil invalid"
            stop
           endif
          enddo
          enddo
          enddo ! i1,j1,k1

          if (weight_total.gt.zero) then

           sumslope=sumslope/weight_total
           slopeterm=slopeterm+sumslope

          else if (weight_total.eq.zero) then
           ! do nothing
          else
           print *,"weight_total invalid"
           stop
          endif 
         else
          print *,"mdata invalid"
          stop
         endif

        else
         print *,"imL or imR bust"
         stop 
        endif 
 
       else
        print *,"im_primary invalid"
        stop
       endif

      else if (total_mass.eq.zero) then
       ! do nothing
      else
       print *,"total_mass invalid"
       stop
      endif

      return
      end subroutine slopecrossterm



      !xsten_accept=stencil for ipart,jpart,kpart,iside_part
      !xsten_recon=stencil for idonate,jdonate,kdonate
      !xsten_donate=stencil for idonate,jdonate,kdonate,isidedonate
      !xsten_depart default value: xsten_accept
      !xsten_target default value: xsten_accept
      subroutine derive_mappings( &
       xsten_accept, &
       xsten_donate, &
       xsten_target, &
       xsten_depart, &
       usten_accept, &
       usten_donate, &
       xdepartsize, &
       xtargetsize, &
       xloint, &
       xhiint, &
       volint, &
       coeff, &
       bfact,dx,map_forward,normdir)
      IMPLICIT NONE

      INTEGER_T map_forward,bfact,normdir
      REAL_T dx(SDIM)
      REAL_T xsten_accept(-1:1,SDIM)
      REAL_T xsten_donate(-1:1,SDIM)
      REAL_T xsten_target(-1:1,SDIM)
      REAL_T xsten_depart(-1:1,SDIM)
      REAL_T usten_accept(-1:1)
      REAL_T usten_donate(-1:1)
      REAL_T xdepartsize,xtargetsize,xloint,xhiint
      REAL_T volint
      REAL_T coeff(2)
      REAL_T coeffINV(2)
      INTEGER_T ihalf

      if ((normdir.lt.0).or.(normdir.ge.SDIM)) then
       print *,"normdir invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid42"
       stop
      endif
      do ihalf=-1,1
       xsten_target(ihalf,normdir+1)=xsten_accept(ihalf,normdir+1)
       xsten_depart(ihalf,normdir+1)=xsten_donate(ihalf,normdir+1)
      enddo

      if (map_forward.eq.0) then ! backwards tracing

       do ihalf=-1,1
        xsten_depart(ihalf,normdir+1)=xsten_accept(ihalf,normdir+1)- &
         usten_accept(ihalf)
       enddo

      else if (map_forward.eq.1) then ! forwards tracing

       do ihalf=-1,1
        xsten_target(ihalf,normdir+1)=xsten_donate(ihalf,normdir+1)+ &
         usten_donate(ihalf)
       enddo

      else
       print *,"map_forward invalid"
       stop
      endif

      xdepartsize=xsten_depart(1,normdir+1)-xsten_depart(-1,normdir+1)
      xtargetsize=xsten_target(1,normdir+1)-xsten_target(-1,normdir+1)

      if ((xdepartsize.lt.VOFTOL*dx(normdir+1)).or. &
          (xdepartsize.gt.two*bfact*dx(normdir+1)).or. &
          (xtargetsize.lt.VOFTOL*dx(normdir+1)).or. &
          (xtargetsize.gt.two*bfact*dx(normdir+1)).or. &
          (xdepartsize.le.VOFTOL*xtargetsize).or. &
          (xtargetsize.le.VOFTOL*xdepartsize)) then
       print *,"WARNING xtarget or xdepart invalid size"
       print *,"xdepartsize ",xdepartsize
       print *,"xtargetsize ",xtargetsize
       print *,"dx= ",dx(normdir+1)
       print *,"map_forward= ",map_forward
       print *,"try lowering cfl<1/2 "
      endif

      if (map_forward.eq.0) then 
        ! xdonate is box that is intersected with the departure region.
       xloint=max(xsten_donate(-1,normdir+1),xsten_depart(-1,normdir+1))
       xhiint=min(xsten_donate(1,normdir+1),xsten_depart(1,normdir+1))
      else if (map_forward.eq.1) then
       xloint=max(xsten_accept(-1,normdir+1),xsten_target(-1,normdir+1))
       xhiint=min(xsten_accept(1,normdir+1),xsten_target(1,normdir+1))
      else
       print *,"map_forward invalid"
       stop
      endif
      volint=xhiint-xloint

      if (volint.gt.VOFTOL*dx(normdir+1)) then
       coeff(1)=xtargetsize/xdepartsize
       coeff(2)=xsten_target(-1,normdir+1)- &
         coeff(1)*xsten_depart(-1,normdir+1)
       coeffINV(1)=xdepartsize/xtargetsize
       coeffINV(2)=xsten_depart(-1,normdir+1)- &
         coeffINV(1)*xsten_target(-1,normdir+1)

       if (map_forward.eq.0) then ! xloint,xhiint in depart region
        xsten_target(-1,normdir+1)=coeff(1)*xloint+coeff(2)
        xsten_target(1,normdir+1)=coeff(1)*xhiint+coeff(2)
        xsten_depart(-1,normdir+1)=xloint
        xsten_depart(1,normdir+1)=xhiint
       else if (map_forward.eq.1) then ! xloint,xhiint in target region
        xsten_depart(-1,normdir+1)=coeffINV(1)*xloint+coeffINV(2)
        xsten_depart(1,normdir+1)=coeffINV(1)*xhiint+coeffINV(2)
        xsten_target(-1,normdir+1)=xloint
        xsten_target(1,normdir+1)=xhiint
       else
        print *,"map_forward invalid"
        stop
       endif

      else
       volint=zero
      endif

      return
      end subroutine derive_mappings

      subroutine get_default_scalar_diffusion(project_option, &
                     LS1,im_source,im_dest, &
                     den, &
                     heatcoeff)
      use global_utility_module
      IMPLICIT NONE
      INTEGER_T, intent(in) :: project_option
      INTEGER_T, intent(in) :: im_source
      INTEGER_T, intent(in) :: im_dest
      REAL_T, intent(in) :: den
      REAL_T, intent(in) :: LS1
      REAL_T, intent(out) :: heatcoeff
      INTEGER_T :: ispec
      INTEGER_T :: nmat

      nmat=num_materials

      if (den.gt.zero) then
       ! do nothing
      else
       print *,"den invalid"
       stop
      endif

      if (project_option.eq.2) then ! thermal diffusion
       if (LS1.ge.zero) then ! center cell owned by im_source
        heatcoeff=get_user_heatviscconst(im_source)
       else  ! center cell owned by im_dest
        heatcoeff=get_user_heatviscconst(im_dest)
       endif
      else if ((project_option.ge.100).and. & ! species diffusion
               (project_option.lt.100+num_species_var)) then
       ispec=project_option-100
       if (LS1.ge.zero) then ! center cell owned by im_source
        heatcoeff=fort_speciesviscconst(ispec*nmat+im_source)*den
       else ! center cell owned by im_dest
        heatcoeff=fort_speciesviscconst(ispec*nmat+im_dest)*den
       endif
      else
       print *,"project_option invalid"
       stop
      endif

      end subroutine get_default_scalar_diffusion


      subroutine interp_from_fluid( &
       LOW, &
       x_fluid, &
       im_secondary_image, &
       thermal_interp, &
       im_fluid, &
       im_solid, &
       LS_interp)
       
      use probf90_module
      use global_utility_module
      use MOF_routines_module
       
      implicit none
       
      type(law_of_wall_parm_type), intent(in) :: LOW
      INTEGER_T, intent(in) :: im_fluid
      INTEGER_T, intent(in) :: im_solid
      REAL_T, intent(in) :: x_fluid(SDIM)
      INTEGER_T, intent(inout) :: im_secondary_image
      REAL_T, intent(out) :: thermal_interp(num_materials)
      REAL_T, intent(out) :: LS_interp(num_materials*(1+SDIM))

      REAL_T :: xsten(-3:3,SDIM)
      REAL_T :: xsten_center(-3:3,SDIM)
      INTEGER_T nhalf
      INTEGER_T im
      INTEGER_T dir
      INTEGER_T cell_index(SDIM)
      INTEGER_T stencil_offset(SDIM)
      INTEGER_T istenlo(3)
      INTEGER_T istenhi(3)
      REAL_T WT,total_WT
      INTEGER_T isten,jsten,ksten
      REAL_T LS_sten(num_materials*(SDIM+1))
      INTEGER_T im_primary_sten
      REAL_T local_temperature

      nhalf=3

      call containing_cell( &
        LOW%bfact, &
        LOW%dx, &
        LOW%xlo, &
        LOW%fablo, &
        x_fluid, &
        cell_index)

      istenlo(3)=0
      istenhi(3)=0
      do dir=1,SDIM
       if (cell_index(dir)-1.lt.LOW%fablo(dir)-LOW%ngrowfab) then
        cell_index(dir)=LOW%fablo(dir)-LOW%ngrowfab+1
       endif
       if (cell_index(dir)+1.gt.LOW%fabhi(dir)+LOW%ngrowfab) then
        cell_index(dir)=LOW%fabhi(dir)+LOW%ngrowfab-1
       endif
       istenlo(dir)=cell_index(dir)-1
       istenhi(dir)=cell_index(dir)+1
      enddo ! dir=1..sdim

      total_WT=zero
      do im=1,LOW%nmat*(1+SDIM)
       LS_interp(im)=zero
      enddo 
      do im=1,LOW%nmat
       thermal_interp(im)=zero
      enddo

      isten=cell_index(1)
      jsten=cell_index(2)
      ksten=cell_index(SDIM)

      call gridsten_level(xsten_center,isten,jsten,ksten,LOW%level,nhalf)

      do isten=istenlo(1),istenhi(1)
      do jsten=istenlo(2),istenhi(2)
      do ksten=istenlo(3),istenhi(3)

       call gridsten_level(xsten,isten,jsten,ksten,LOW%level,nhalf)
       stencil_offset(1)=isten-cell_index(1)
       stencil_offset(2)=jsten-cell_index(2)
       if (SDIM.eq.3) then
        stencil_offset(SDIM)=ksten-cell_index(SDIM)
       endif
       call bilinear_interp_WT(xsten_center,nhalf,stencil_offset, &
        x_fluid,WT)
       if ((WT.ge.zero).and.(WT.le.one)) then
        ! do nothing
       else
        print *,"WT invalid"
        stop
       endif
 
       do im=1,LOW%nmat*(1+SDIM)
        LS_sten(im)=LOW%LSCP(D_DECL(isten,jsten,ksten),im)
       enddo

       call get_primary_material(LS_sten,LOW%nmat,im_primary_sten)
       if (im_primary_sten.eq.im_fluid) then
        ! do nothing
       else if (im_primary_sten.eq.im_solid) then
        ! do nothing
       else if ((im_primary_sten.ge.1).and. &
                (im_primary_sten.le.LOW%nmat)) then
        if (is_rigid(LOW%nmat,im_primary_sten).eq.1) then
         ! do nothing
        else if (is_rigid(LOW%nmat,im_primary_sten).eq.0) then
         if (abs(LS_sten(im_fluid)).le.LOW%dxmin*GNBC_RADIUS) then
          if (im_secondary_image.eq.0) then
           im_secondary_image=im_primary_sten
          endif
         else if (abs(LS_sten(im_fluid)).ge.LOW%dxmin*GNBC_RADIUS) then
          ! do nothing
         else
          print *,"LS_sten became corrupt"
          stop
         endif
        else
         print *,"is_rigid(LOW%nmat,im_primary_sten) invalid"
         stop
        endif
       else
        print *,"im_primary_sten invalid"
        stop
       endif
       do im=1,LOW%nmat*(1+SDIM)
        LS_interp(im)=LS_interp(im)+WT*LS_sten(im)
       enddo
       do im=1,LOW%nmat
        local_temperature= &
          LOW%state(D_DECL(isten,jsten,ksten), &
             (im-1)*num_state_material+2) 
        if (local_temperature.ge.zero) then
         thermal_interp(im)=thermal_interp(im)+WT*local_temperature
        else
         print *,"local_temperature invalid"
         stop
        endif
       enddo ! im=1..nmat
       total_WT=total_WT+WT

      enddo ! ksten
      enddo ! jsten
      enddo ! isten

      if (total_WT.gt.zero) then

       do im=1,LOW%nmat*(1+SDIM)
        LS_interp(im)=LS_interp(im)/total_WT
       enddo
       do im=1,LOW%nmat
        thermal_interp(im)=thermal_interp(im)/total_WT
       enddo

      else
       print *,"total_WT invalid"
       stop
      endif

      return
      end subroutine interp_from_fluid


       ! This routine transforms all velocities to a frame of reference with
       ! respect to the solid.
      subroutine getGhostVel( &
       LOW, &
       law_of_the_wall, &
       iSD,jSD,kSD, &
       iFD,jFD,kFD, &
       side_solid, &
       side_image, &
       data_dir, &
       uimage_raster, &
       usolid_law_of_wall, & ! out
       angle_ACT, & ! aka angle_ACT_cell "out"
       im_fluid, &
       im_solid)
       
       use probf90_module
       use global_utility_module
       use MOF_routines_module
       
       implicit none
       
       type(law_of_wall_parm_type), intent(in) :: LOW
       INTEGER_T, intent(in) :: law_of_the_wall
       INTEGER_T, intent(in) :: data_dir
       INTEGER_T, intent(in) :: im_fluid
       INTEGER_T, intent(in) :: im_solid
       INTEGER_T, intent(in) :: side_solid
       INTEGER_T, intent(in) :: side_image
       INTEGER_T, intent(in) :: iSD,jSD,kSD
       INTEGER_T, intent(in) :: iFD,jFD,kFD
       REAL_T, dimension(SDIM), intent(inout) :: uimage_raster
       REAL_T, dimension(SDIM), intent(out) :: usolid_law_of_wall
       REAL_T, intent(out) :: angle_ACT
      
       !D_DECL is defined in SPACE.H in the BoxLib/Src/C_BaseLib


       REAL_T, dimension(SDIM) :: u_tngt
       REAL_T :: uimage_nrml, ughost_nrml
       REAL_T :: ughost_tngt
       REAL_T :: uimage_tngt_mag
       REAL_T :: tau_w
       REAL_T :: viscosity_molecular, viscosity_eddy
       REAL_T :: density_fluid
       INTEGER_T :: dir
       REAL_T :: nrm_sanity
       REAL_T :: nrm_sanity_crossing
       REAL_T :: predict_deriv_utan
       REAL_T :: max_deriv_utan
       REAL_T :: critical_length
       REAL_T :: uimage_mag
       INTEGER_T :: im
       INTEGER_T :: im_fluid_crossing
       INTEGER_T :: im_primary_image
       INTEGER_T :: im_secondary_image
       INTEGER_T :: im_fluid1,im_fluid2
       INTEGER_T :: iten
       INTEGER_T :: iten_13,iten_23
       INTEGER_T :: nten_test
       REAL_T :: cos_angle,sin_angle
       INTEGER_T :: near_contact_line
       REAL_T :: mag
       REAL_T :: sinthetaACT
       REAL_T :: costhetaACT
       REAL_T :: dist_to_CL
       REAL_T :: nf_dot_ns
       REAL_T :: nf_crossing_dot_ns
       REAL_T :: nf_dot_nCL_perp
       REAL_T :: nf_crossing_dot_nCL_perp
       REAL_T, dimension(3) :: nCL
       REAL_T, dimension(3) :: nCL_crossing
       REAL_T, dimension(3) :: nCL_raster
       REAL_T, dimension(3) :: nCL_perp
       REAL_T, dimension(3) :: nCL_perp_crossing
       REAL_T, dimension(3) :: nCL_perp2
       REAL_T, dimension(3) :: nCL_perp2_crossing
       REAL_T, dimension(3) :: nf_prj
       REAL_T, dimension(3) :: nf_prj_crossing
       REAL_T :: ZEYU_mu_l, ZEYU_mu_g, ZEYU_sigma
       REAL_T :: ZEYU_thet_s,ZEYU_lambda,ZEYU_l_macro, ZEYU_l_micro
       REAL_T :: ZEYU_dgrid, ZEYU_d_closest, ZEYU_thet_d_apparent
       REAL_T :: ZEYU_u_cl, ZEYU_u_slip, ZEYU_thet_d
       REAL_T :: angle_im1
       INTEGER_T :: ZEYU_imodel
       INTEGER_T :: ZEYU_ifgnbc
       INTEGER_T :: im_vapor,im_liquid
       REAL_T :: delta_r_raster
       REAL_T :: nCL_dot_n_raster
       INTEGER_T :: nrad
       INTEGER_T nhalf
       REAL_T :: xstenFD(-3:3,SDIM)
       REAL_T :: xstenSD(-3:3,SDIM)
       REAL_T :: thermal_interp(num_materials)
       REAL_T :: LSPLUS_interp(num_materials*(1+SDIM))
       REAL_T :: LSMINUS_interp(num_materials*(1+SDIM))
       REAL_T :: LSTRIPLE_interp(num_materials*(1+SDIM))
       REAL_T :: LS_fluid(num_materials*(1+SDIM))
       REAL_T :: LS_solid(num_materials*(1+SDIM))
       REAL_T :: LS_crossing(num_materials*(1+SDIM))
       REAL_T :: LS_triple(num_materials*(1+SDIM))
       REAL_T :: nrm_solid(3)
       REAL_T :: nrm_fluid(3)
       REAL_T :: nrm_fluid_crossing(3)
       REAL_T, allocatable, dimension(:) :: user_tension
       REAL_T :: cross_denom
       REAL_T :: cross_factor
       REAL_T :: cross_factorMINUS
       REAL_T :: cross_factorPLUS
       INTEGER_T :: cross_factor_flag
       REAL_T :: xcrossing(SDIM)
       REAL_T :: xprobeMINUS_crossing(SDIM)
       REAL_T :: xprobePLUS_crossing(SDIM)
       REAL_T :: xprobe_triple(SDIM)
       REAL_T :: xtriple(SDIM)
       INTEGER_T :: debug_slip_velocity_enforcement

       debug_slip_velocity_enforcement=0
    
       nhalf=3 
       nten_test=( (LOW%nmat-1)*(LOW%nmat-1)+LOW%nmat-1 )/2
       allocate(user_tension(nten_test))
       if (LOW%nten.eq.nten_test) then
        ! do nothing
       else
        print *,"nten invalid"
        stop
       endif
        
       if ((data_dir.ge.0).and.(data_dir.lt.SDIM)) then
        ! do nothing
       else
        print *,"data_dir invalid"
        stop
       endif
       if ((im_fluid.lt.1).or.(im_fluid.gt.num_materials)) then
        print *,"im_fluid invalid in getGhostVel"
        stop
       endif
       if ((im_solid.lt.1).or.(im_solid.gt.num_materials)) then
        print *,"im_solid invalid in getGhostVel"
        stop
       endif
       if (LOW%ngrowfab.eq.4) then
        ! do nothing
       else
        print *,"LOW%ngrowfab invalid"
        stop
       endif
       if (is_rigid(LOW%nmat,im_solid).eq.1) then
        ! do nothing
       else
        print *,"is_rigid(nmat,im_solid) invalid"
        stop
       endif
       if (LOW%nmat.ne.num_materials) then
        print *,"nmat invalid"
        stop
       endif
       if (LOW%dt.gt.zero) then
        ! do nothing
       else
        print *,"dt invalid"
        stop
       endif 
       if (LOW%time.ge.zero) then
        ! do nothing
       else
        print *,"time invalid"
        stop
       endif 
       if (LOW%visc_coef.ge.zero) then
        ! do nothing
       else
        print *,"visc_coef invalid"
        stop
       endif 
       if ((law_of_the_wall.eq.1).or. &
           (law_of_the_wall.eq.2)) then
        ! do nothing
       else
        print *,"law_of_the_wall invalid"
        stop
       endif

       if (LOW%dxmin.gt.zero) then
        ! do nothing
       else
        print *,"LOW%dxmin invalid"
        stop
       endif

       nrad=3

        ! uimage_raster and LOW%usolid_raster are already initialized.
       do dir=1,SDIM
        usolid_law_of_wall(dir)=zero
        u_tngt(dir)=zero
       enddo
       ughost_nrml=zero
       ughost_tngt=zero
       angle_ACT=zero

       call gridsten_level(xstenFD,iFD,jFD,kFD,LOW%level,nhalf)
       do im=1,LOW%nmat*(1+SDIM)
        LS_fluid(im)=LOW%LSCP(D_DECL(iFD,jFD,kFD),im)
       enddo
       call gridsten_level(xstenSD,iSD,jSD,kSD,LOW%level,nhalf)
       do im=1,LOW%nmat*(1+SDIM)
        LS_solid(im)=LOW%LSCP(D_DECL(iSD,jSD,kSD),im)
       enddo

       near_contact_line=0
       im_primary_image=im_fluid

       if (LS_solid(im_solid).ge.zero) then
        if (LS_fluid(im_fluid).ge.zero) then ! im_fluid dominates the fluids
         if (LS_fluid(im_solid).le.zero) then
          cross_denom=LS_solid(im_solid)-LS_fluid(im_solid)
          if (cross_denom.gt.zero) then
           cross_factor=LS_solid(im_solid)/cross_denom
          else if (cross_denom.eq.zero) then
           cross_factor=half
          else
           print *,"cross_factor invalid"
           stop
          endif
          if ((cross_factor.ge.zero).and. &
              (cross_factor.le.one)) then
              ! xcrossing is where the solid interface passes inbetween
              ! the ghost (solid) point and the fluid (image) point.
           do dir=1,SDIM
            xcrossing(dir)=cross_factor*xstenFD(0,dir)+ &
                    (one-cross_factor)*xstenSD(0,dir)
           enddo

           do im=1,LOW%nmat*(1+SDIM)
            LS_crossing(im)=cross_factor*LS_fluid(im)+ &
                     (one-cross_factor)*LS_solid(im)
           enddo ! im=1..nmat*(1+SDIM)

           im_fluid_crossing=-1
           do im=1,LOW%nmat
            if (is_rigid(LOW%nmat,im).eq.1) then
             ! do nothing
            else if (is_rigid(LOW%nmat,im).eq.0) then
             if (im_fluid_crossing.eq.-1) then
              im_fluid_crossing=im
             else if ((im_fluid_crossing.ge.1).and. &
                      (im_fluid_crossing.le.LOW%nmat)) then
              if (LS_crossing(im_fluid_crossing).le.LS_crossing(im)) then
               im_fluid_crossing=im
              endif
             else
              print *,"im_fluid_crossing invalid"
              stop
             endif
            else 
             print *,"is_rigid invalid"
             stop
            endif
           enddo ! im=1..nmat

           call normalize_LS_normals(LOW%nmat,LS_crossing)

           if ((im_fluid_crossing.ge.1).and. &
               (im_fluid_crossing.le.LOW%nmat)) then
            ! do nothing
           else
            print *,"im_fluid_crossing invalid"
            stop
           endif

           nf_dot_ns=zero
           nf_crossing_dot_ns=zero
           nrm_fluid(3)=zero
           nrm_fluid_crossing(3)=zero
           nrm_solid(3)=zero
           do dir=1,SDIM
             ! points into im_fluid material
            nrm_fluid(dir)= &
                LS_crossing(LOW%nmat+(im_fluid-1)*SDIM+dir)
            nrm_fluid_crossing(dir)= &
                LS_crossing(LOW%nmat+(im_fluid_crossing-1)*SDIM+dir)
             ! points into im_solid material
            nrm_solid(dir)=LS_crossing(LOW%nmat+(im_solid-1)*SDIM+dir)
            nf_dot_ns=nf_dot_ns+ &
                  nrm_fluid(dir)*nrm_solid(dir)
            nf_crossing_dot_ns=nf_crossing_dot_ns+ &
                  nrm_fluid_crossing(dir)*nrm_solid(dir)
           enddo ! dir=1..sdim

              !    \
              !nCL  \
              !  <-- \
              ! ------O----
              ! 
              ! nCL is tangential to the solid/fluid interface
              ! since nrm_solid has been projected away from nrm_fluid
              ! nCL_perp points out of the screen.
           nrm_sanity=zero
           nrm_sanity_crossing=zero
           do dir=1,3
            nCL(dir)=nrm_fluid(dir)- &
                    nf_dot_ns*nrm_solid(dir)
            nrm_sanity=nrm_sanity+nCL(dir)**2
            nCL_crossing(dir)=nrm_fluid_crossing(dir)- &
                    nf_crossing_dot_ns*nrm_solid(dir)
            nrm_sanity_crossing=nrm_sanity_crossing+nCL_crossing(dir)**2
           enddo 
           nrm_sanity=sqrt(nrm_sanity)
           nrm_sanity_crossing=sqrt(nrm_sanity_crossing)

           if (nrm_sanity.gt.zero) then
            do dir=1,3
             nCL(dir)=nCL(dir)/nrm_sanity
            enddo
           else if (nrm_sanity.eq.zero) then
            ! do nothing
           else
            print *,"nrm_sanity invalid"
            stop
           endif

           if (nrm_sanity_crossing.gt.zero) then
            do dir=1,3
             nCL_crossing(dir)=nCL_crossing(dir)/nrm_sanity_crossing
            enddo
           else if (nrm_sanity_crossing.eq.zero) then
            ! do nothing
           else
            print *,"nrm_sanity_crossing invalid"
            stop
           endif

           ! nCL_perp is tangent to the contact line in the substrate plane.
           ! nCL is normal to the contact line in the substrate plane
           ! nrm_solid is normal to the substrate
           call crossprod(nCL,nrm_solid,nCL_perp)
           call crossprod(nCL_crossing,nrm_solid,nCL_perp_crossing)

           nrm_sanity=zero
           nrm_sanity_crossing=zero
           do dir=1,3
            nrm_sanity=nrm_sanity+nCL_perp(dir)**2
            nrm_sanity_crossing=nrm_sanity_crossing+nCL_perp_crossing(dir)**2
           enddo 
           nrm_sanity=sqrt(nrm_sanity)
           nrm_sanity_crossing=sqrt(nrm_sanity_crossing)

           if (nrm_sanity.gt.zero) then
            do dir=1,3
             nCL_perp(dir)=nCL_perp(dir)/nrm_sanity
            enddo
           else if (nrm_sanity.eq.zero) then
            ! do nothing
           else
            print *,"nrm_sanity invalid"
            stop
           endif

           if (nrm_sanity_crossing.gt.zero) then
            do dir=1,3
             nCL_perp_crossing(dir)=nCL_perp_crossing(dir)/nrm_sanity_crossing
            enddo
           else if (nrm_sanity_crossing.eq.zero) then
            ! do nothing
           else
            print *,"nrm_sanity_crossing invalid"
            stop
           endif

           nf_dot_nCL_perp=zero
           nf_crossing_dot_nCL_perp=zero
           do dir=1,3

            nf_dot_nCL_perp= &
               nf_dot_nCL_perp+ &
               nrm_fluid(dir)*nCL_perp(dir)

            nf_crossing_dot_nCL_perp= &
               nf_crossing_dot_nCL_perp+ &
               nrm_fluid_crossing(dir)*nCL_perp_crossing(dir)

           enddo
           ! nCL_perp is tangent to the contact line in the substrate plane.
           ! nrm_fluid points into im_fluid material.
           ! nf_prj is the fluid normal with the tangent contact line
           ! vector projected away; nf_prj, in a plane perpendicular
           ! to the substrate and perpendicular to the contact line,
           ! is the normal point into im_primary
           ! "im_primary_image"
           nrm_sanity=zero
           nrm_sanity_crossing=zero
           do dir=1,3
            nf_prj(dir)=nrm_fluid(dir)- &
                 nf_dot_nCL_perp*nCL_perp(dir)
            nrm_sanity=nrm_sanity+nf_prj(dir)**2

            nf_prj_crossing(dir)=nrm_fluid_crossing(dir)- &
                 nf_crossing_dot_nCL_perp*nCL_perp_crossing(dir)
            nrm_sanity_crossing=nrm_sanity_crossing+nf_prj_crossing(dir)**2
           enddo

           nrm_sanity=sqrt(nrm_sanity)
           nrm_sanity_crossing=sqrt(nrm_sanity_crossing)

           if (nrm_sanity.gt.zero) then
            do dir=1,3
             nf_prj(dir)=nf_prj(dir)/nrm_sanity
            enddo
           else if (nrm_sanity.eq.zero) then
            ! do nothing
           else
            print *,"nrm_sanity invalid"
            stop
           endif

           if (nrm_sanity_crossing.gt.zero) then
            do dir=1,3
             nf_prj_crossing(dir)=nf_prj_crossing(dir)/nrm_sanity_crossing
            enddo
           else if (nrm_sanity_crossing.eq.zero) then
            ! do nothing
           else
            print *,"nrm_sanity_crossing invalid"
            stop
           endif


           ! nCL_perp2 will also be tangent to the contact line.
           call crossprod(nrm_solid,nf_prj,nCL_perp2)
           call crossprod(nrm_solid,nf_prj_crossing,nCL_perp2_crossing)

            ! now we update nrm_fluid by finding the triple point
            ! and measuring nrm_fluid just above the solid.
            ! LS_crossing is the levelset function at the point in 
            ! between the ghost point and the image point at which
            ! LS_crossing(im_solid)=0
            ! NOTE: LS_IMAGE(im_fluid) > 0
            !      x  image point in the fluid , LS(im_fluid)>0 here
            !   --------   substrate interface
            !      x  ghost point in the substrate, LS(im_fluid) expected to
            !      be positive here too since the level set function is
            !      extrapolated into the substrate normal to the substrate.
            ! This routine is only called if LS_IMAGE(im_solid) <0 and
            ! LS_GHOST(im_solid)>0.
            ! xprobe_crossing is another point on the substrate:
            !
            !                       n_im_fluid_project
            !    xprobe_crossing     <----         x image  im_fluid
            !   ----x------------------------------x---------- crossing
            !                                      x ghost  im_fluid
            !  if LS_{im_fluid}(xprobe_crossing) *
            !     LS_{im_fluid}(xcrossing) < 0 =>
            !  a triple points exists in between.
            !  if LS_{im_fluid}(xcrossing) * LS_{im_fluid}(xghost) < 0 =>
            !  the triple point is already in the immediate vicinity, no
            !  need to search for it.
           if ((GNBC_RADIUS.ge.one).and. &
               (GNBC_RADIUS.le.three)) then
            do dir=1,SDIM
             xprobeMINUS_crossing(dir)=xcrossing(dir)- &
               GNBC_RADIUS*LOW%dxmin*nCL_crossing(dir)
             xprobePLUS_crossing(dir)=xcrossing(dir)+ &
               GNBC_RADIUS*LOW%dxmin*nCL_crossing(dir)
            enddo ! dir=1..sdim

            im_secondary_image=0

            call interp_from_fluid( &
             LOW, &
             xprobeMINUS_crossing, &
             im_secondary_image, &
             thermal_interp, &
             im_fluid, &
             im_solid, &
             LSMINUS_interp)

            call normalize_LS_normals(LOW%nmat,LSMINUS_interp)

            call interp_from_fluid( &
             LOW, &
             xprobePLUS_crossing, &
             im_secondary_image, &
             thermal_interp, &
             im_fluid, &
             im_solid, &
             LSPLUS_interp)

            call normalize_LS_normals(LOW%nmat,LSPLUS_interp)

            cross_factorMINUS=-one
            cross_factorPLUS=-one

            if (LSMINUS_interp(im_fluid_crossing)* &
                LS_crossing(im_fluid_crossing).le.zero) then
             cross_denom=LS_crossing(im_fluid_crossing)- &
                         LSMINUS_interp(im_fluid_crossing)
             if (cross_denom.ne.zero) then
              cross_factorMINUS=LS_crossing(im_fluid_crossing)/cross_denom
             else if (cross_denom.eq.zero) then
              cross_factorMINUS=half
             else
              print *,"cross_denom invalid"
              stop
             endif
            endif

            if (LSPLUS_interp(im_fluid_crossing)* &
                LS_crossing(im_fluid_crossing).le.zero) then
             cross_denom=LS_crossing(im_fluid_crossing)- &
                         LSPLUS_interp(im_fluid_crossing)
             if (cross_denom.ne.zero) then
              cross_factorPLUS=LS_crossing(im_fluid_crossing)/cross_denom
             else if (cross_denom.eq.zero) then
              cross_factorPLUS=half
             else
              print *,"cross_denom invalid"
              stop
             endif
            endif

            if ((cross_factorPLUS.eq.-one).and. &
                (cross_factorMINUS.eq.-one)) then
             cross_factor_flag=0
            else if ((cross_factorPLUS.eq.-one).and. &
                     (cross_factorMINUS.ge.zero)) then
             cross_factor_flag=-1
            else if ((cross_factorMINUS.eq.-one).and. &
                     (cross_factorPLUS.ge.zero)) then
             cross_factor_flag=1
            else if ((cross_factorMINUS.ge.zero).and. &
                     (cross_factorPLUS.ge.zero)) then
             if (cross_factorPLUS.le.cross_factorMINUS) then
              cross_factor_flag=1
             else if (cross_factorPLUS.ge.cross_factorMINUS) then
              cross_factor_flag=-1
             else
              print *,"cross_factor bust"
              stop
             endif
            else
             print *,"cross_factor bust"
             stop
            endif
           
            if (cross_factor_flag.eq.1) then 

             if ((cross_factorPLUS.ge.zero).and. &
                 (cross_factorPLUS.le.one)) then
              do dir=1,SDIM
               xtriple(dir)=cross_factorPLUS*xprobePLUS_crossing(dir)+ &
                  (one-cross_factorPLUS)*xcrossing(dir)
              enddo
              do im=1,LOW%nmat*(1+SDIM)
               LS_triple(im)=cross_factorPLUS*LSPLUS_interp(im)+ &
                   (one-cross_factorPLUS)*LS_crossing(im)
              enddo
             else
              print *,"cross_factorPLUS bust"
              stop
             endif
             
            else if (cross_factor_flag.eq.-1) then 

             if ((cross_factorMINUS.ge.zero).and. &
                 (cross_factorMINUS.le.one)) then
              do dir=1,SDIM
               xtriple(dir)=cross_factorMINUS*xprobeMINUS_crossing(dir)+ &
                  (one-cross_factorMINUS)*xcrossing(dir)
              enddo
              do im=1,LOW%nmat*(1+SDIM)
               LS_triple(im)=cross_factorMINUS*LSMINUS_interp(im)+ &
                   (one-cross_factorMINUS)*LS_crossing(im)
              enddo
             else
              print *,"cross_factorMINUS bust"
              stop
             endif

            else if (cross_factor_flag.eq.0) then
             ! do nothing
            else
             print *,"cross_factor_flag invalid"
             stop
            endif

            if ((cross_factor_flag.eq.1).or. &
                (cross_factor_flag.eq.-1)) then
              call normalize_LS_normals(LOW%nmat,LS_triple)

              nrm_solid(3)=zero
              do dir=1,SDIM
               nrm_solid(dir)=LS_triple(LOW%nmat+(im_solid-1)*SDIM+dir)
               xprobe_triple(dir)=xtriple(dir)-LOW%dxmin*nrm_solid(dir)
              enddo
              
              call interp_from_fluid( &
               LOW, &
               xprobe_triple, &
               im_secondary_image, &
               thermal_interp, &
               im_fluid, &
               im_solid, &
               LSTRIPLE_interp)

              call normalize_LS_normals(LOW%nmat,LSTRIPLE_interp)
              nrm_fluid(3)=zero
              do dir=1,SDIM
               nrm_fluid(dir)=LSTRIPLE_interp(LOW%nmat+(im_fluid-1)*SDIM+dir)
              enddo

              nf_dot_ns=zero
              do dir=1,SDIM
               nf_dot_ns=nf_dot_ns+nrm_fluid(dir)*nrm_solid(dir)
              enddo ! dir=1..sdim
 
              nrm_sanity=zero
              do dir=1,3
               nCL(dir)=nrm_fluid(dir)-nf_dot_ns*nrm_solid(dir)
               nrm_sanity=nrm_sanity+nCL(dir)**2
              enddo 
              nrm_sanity=sqrt(nrm_sanity)
              if (nrm_sanity.gt.zero) then
               do dir=1,3
                nCL(dir)=nCL(dir)/nrm_sanity
               enddo
              else if (nrm_sanity.eq.zero) then
               ! do nothing
              else
               print *,"nrm_sanity invalid"
               stop
              endif

              call crossprod(nCL,nrm_solid,nCL_perp)

              nrm_sanity=zero
              do dir=1,3
               nrm_sanity=nrm_sanity+nCL_perp(dir)**2
              enddo 
              nrm_sanity=sqrt(nrm_sanity)
              if (nrm_sanity.gt.zero) then
               do dir=1,3
                nCL_perp(dir)=nCL_perp(dir)/nrm_sanity
               enddo
              else if (nrm_sanity.eq.zero) then
               ! do nothing
              else
               print *,"nrm_sanity invalid"
               stop
              endif

              nf_dot_nCL_perp=zero
              do dir=1,3
               nf_dot_nCL_perp=nf_dot_nCL_perp+nrm_fluid(dir)*nCL_perp(dir)
              enddo
              nrm_sanity=zero
              do dir=1,3
               nf_prj(dir)=nrm_fluid(dir)-nf_dot_nCL_perp*nCL_perp(dir)
               nrm_sanity=nrm_sanity+nf_prj(dir)**2
              enddo
              nrm_sanity=sqrt(nrm_sanity)
              if (nrm_sanity.gt.zero) then
               do dir=1,3
                nf_prj(dir)=nf_prj(dir)/nrm_sanity
               enddo
              else if (nrm_sanity.eq.zero) then
               ! do nothing
              else
               print *,"nrm_sanity invalid"
               stop
              endif
              call crossprod(nrm_solid,nf_prj,nCL_perp2)

              near_contact_line=1
              if (im_secondary_image.eq.0) then
               near_contact_line=0
              else if ((im_secondary_image.ge.1).and. &
                       (im_secondary_image.le.LOW%nmat)) then
               if (is_rigid(LOW%nmat,im_secondary_image).eq.1) then
                near_contact_line=0
               else if (is_rigid(LOW%nmat,im_secondary_image).eq.0) then
                ! do nothing
               else
                print *,"is_rigid(nmat,im_secondary_image) invalid"
                stop
               endif
              else
               print *,"im_secondary_image invalid"
               stop
              endif

              if (near_contact_line.eq.1) then
               if (im_primary_image.lt.im_secondary_image) then
                im_fluid1=im_primary_image
                im_fluid2=im_secondary_image
               else if (im_primary_image.gt.im_secondary_image) then
                im_fluid2=im_primary_image
                im_fluid1=im_secondary_image
               else
                print *,"im_primary_image or im_secondary_image invalid"
                stop
               endif
               call get_iten(im_fluid1,im_fluid2,iten,LOW%nmat)
                ! in: subroutine getGhostVel
               call get_user_tension(xprobe_triple,LOW%time, &
                 fort_tension,user_tension, &
                 thermal_interp,LOW%nmat,LOW%nten,6)
               ! cos_angle and sin_angle correspond to the angle in im_fluid1
               call get_CL_iten(im_fluid1,im_fluid2,im_solid, &
                 iten_13,iten_23, &
                 user_tension,LOW%nten,cos_angle,sin_angle)
               sinthetaACT=zero
               costhetaACT=zero
                ! because nrm_solid and nf_prj have unit magnitude,
                ! sinthetaACT is the sine of the angle between nrm_solid
                ! and nf_prj; fluid normal in plane perpendicular to 
                ! contact line and solid.
               do dir=1,3
                sinthetaACT=sinthetaACT+nCL_perp2(dir)**2
                costhetaACT=costhetaACT+nrm_solid(dir)*nf_prj(dir)
               enddo
               sinthetaACT=sqrt(sinthetaACT)

               dist_to_CL=zero
               do dir=1,SDIM
                dist_to_CL=dist_to_CL+(xtriple(dir)-xcrossing(dir))**2
               enddo
               dist_to_CL=sqrt(dist_to_CL)
                
               if ((sinthetaACT.ge.zero).and.(costhetaACT.ge.zero)) then
                angle_ACT=asin(sinthetaACT)
               else if ((sinthetaACT.ge.zero).and.(costhetaACT.le.zero)) then
                angle_ACT=Pi-asin(sinthetaACT)
               else
                print *,"sinthetaACT or costhetaACT invalid"
                stop
               endif
               if (DEBUG_DYNAMIC_CONTACT_ANGLE.eq.1) then
                print *,"xcrossing ",xcrossing(1),xcrossing(2),xcrossing(SDIM)
                print *,"xtriple ",xtriple(1),xtriple(2),xtriple(SDIM)
                print *,"nrm_solid ",nrm_solid(1),nrm_solid(2),nrm_solid(SDIM)
                print *,"nrm_fluid ",nrm_fluid(1),nrm_fluid(2),nrm_fluid(SDIM)
                print *,"im_fluid,angle_ACT(rad,deg) ",im_fluid,angle_ACT, &
                        angle_ACT*180.0d0/Pi
                print *,"dx(1),dist_to_CL ",LOW%dx(1),dist_to_CL
                print *,"im_primary_image,im_secondary_image ", &
                        im_primary_image,im_secondary_image
               endif
              else if (near_contact_line.eq.0) then
               ! do nothing
              else
               print *,"near_contact_line invalid"
               stop
              endif
            else if (cross_factor_flag.eq.0) then
              ! do nothing
            else
             print *,"cross_factor_flag invalid"
             stop
            endif
           else
            print *,"GNBC_RADIUS invalid"
            stop
           endif
          else
           print *,"cross_factor invalid"
           stop
          endif
         else
          print *,"expecting LS_fluid(im_solid)<=0"
          stop
         endif
        else
         print *,"expecting LS_fluid(im_fluid)>=0"
         stop
        endif
       else
        print *,"expecting LS_solid(im_solid)>=0"
        stop
       endif

       viscosity_molecular=fort_viscconst(im_fluid)
       viscosity_eddy=fort_viscconst_eddy(im_fluid)
       density_fluid=fort_denconst(im_fluid)

       if (density_fluid.gt.zero) then
        ! do nothing
       else
        print *,"density_fluid invalid"
        stop
       endif

       !x_projection is closest point on the fluid/solid interface. 
       delta_r_raster=abs(LOW%x_image_raster(data_dir+1)- &
                          LOW%x_projection_raster(data_dir+1))
       if (delta_r_raster.gt.zero) then
        ! do nothing
       else
        print *,"delta_r_raster invalid"
        stop
       endif
      
       ! convert to solid velocity frame of reference.
       uimage_mag=zero
       do dir=1,SDIM
        uimage_raster(dir)=uimage_raster(dir)-LOW%usolid_raster(dir)
        uimage_mag=uimage_mag+uimage_raster(dir)**2
       enddo
       uimage_mag=sqrt(uimage_mag)
       if (uimage_mag.ge.zero) then
        ! do nothing
       else
        print *,"uimage_mag invalid"
        stop
       endif

       !normal and tangential velocity components of image point
       !magnitude, normal velocity component
       uimage_nrml = DOT_PRODUCT(uimage_raster,LOW%n_raster) 
       do dir = 1,SDIM
        u_tngt(dir) = uimage_raster(dir)-uimage_nrml*LOW%n_raster(dir)
       enddo

       uimage_tngt_mag = DOT_PRODUCT(u_tngt,u_tngt)
       uimage_tngt_mag = sqrt(uimage_tngt_mag) 

       if (uimage_tngt_mag.lt.zero) then
        print *,"uimage_tngt_mag.lt.zero"
        stop
       else if (uimage_tngt_mag.eq.zero) then
        ! do nothing
       else if (uimage_tngt_mag.gt.zero) then 
        !normalize tangential velocity vector
        do dir = 1,SDIM
         u_tngt(dir) = u_tngt(dir)/uimage_tngt_mag 
        enddo
       else
        print *,"uimage_tngt_mag invalid"
        print *,"uimage_tngt_mag=",uimage_tngt_mag
        print *,"uimage_nrml=",uimage_nrml
        do dir=1,SDIM
         print *,"dir,uimage_raster(dir) ",dir,uimage_raster(dir)
         print *,"dir,usolid_raster(dir) ",dir,LOW%usolid_raster(dir)
         print *,"dir,u_tngt(dir) ",dir,u_tngt(dir)
        enddo
        stop
       endif

       ! laminar boundary layer thickness:
       ! delta=4.91 (mu x / (U rho))^(1/2)
       ! dx=4.91 (mu dx/ (U rho))^(1/2)
       ! dx=4.91^2 mu/(U rho)
       !   =(g/(cm s))/(cm/s  g/cm^3)=(g/(cm s))/(g/(cm^2 s))=cm
       if (uimage_tngt_mag.gt.zero) then
        critical_length=((4.91D0)**2)*viscosity_molecular/ &
              (density_fluid*uimage_tngt_mag)
       else if (uimage_tngt_mag.eq.zero) then
        critical_length=LOW%dxmin*1.0D+10
       else
        print *,"uimage_tngt_mag invalid"
        stop
       endif

       if ((viscosity_molecular.gt.zero).and. &
           (viscosity_eddy.ge.zero)) then
        ! do nothing
       else
        print *,"viscosity_molecular.le.zero or viscosity_eddy.lt.zero"
        stop
       endif

        ! ghost velocity lives *on* the rasterized interface.
       ughost_nrml = zero

       ! From Spaldings' paper, a representative size for the linear 
       ! region is:
       ! y+ = 10.0
       ! y+ = y sqrt(tau rho)/mu_molecular
       !  or y+ is the intersection point of the linear (viscous) 
       !  and log layer profiles.

       if (law_of_the_wall.eq.1) then ! turbulence modeling here.

        if (critical_length.lt.LOW%dxmin) then

         if (viscosity_eddy.gt.zero) then
          !obtain wall shear stress tau_w
          !out tau_w

          if (delta_r_raster.gt.zero) then
           call wallFunc_NewtonsMethod(uimage_tngt_mag, &
                 delta_r_raster,tau_w,im_fluid) 
           ughost_tngt = uimage_tngt_mag - &
            tau_w*delta_r_raster/ &
            (viscosity_molecular+viscosity_eddy)

           predict_deriv_utan=abs(ughost_tngt-uimage_tngt_mag)/ &
              delta_r_raster
           max_deriv_utan=two*uimage_tngt_mag/critical_length
           if (predict_deriv_utan.lt.max_deriv_utan) then
            ! do nothing
           else
            print *,"predict_deriv_utan or max_deriv_utan invalid"
            print *,"predict_deriv_utan= ",predict_deriv_utan
            print *,"max_deriv_utan= ",max_deriv_utan
            stop
           endif
          else
           print *,"delta_r_raster invalid"
           stop
          endif

         else if (viscosity_eddy.eq.zero) then
           ! ghost velocity lives *on* the rasterized interface.
          ughost_tngt = zero
         else
          print *,"viscosity_eddy invalid"
          stop
         endif

        else if (critical_length.ge.LOW%dxmin) then
          ! ghost velocity lives *on* the rasterized interface.
         ughost_tngt = zero
        else
         print *,"critical_length invalid"
         stop
        endif

       else if (law_of_the_wall.eq.2) then ! GNBC model

        if (near_contact_line.eq.1) then

         if ((fort_denconst(im_fluid1).gt.zero).and. &
             (fort_denconst(im_fluid2).gt.zero)) then
 
          if ((sin_angle.ge.zero).and.(cos_angle.ge.zero)) then
           angle_im1=asin(sin_angle)
           ! sin_angle=sin(a)  cos_angle=cos(a)
           ! a=pi-asin(sin_angle)
          else if ((sin_angle.ge.zero).and.(cos_angle.le.zero)) then
           angle_im1=Pi-asin(sin_angle)
          else
           print *,"sin_angle or cos_angle invalid"
           stop
          endif

          ZEYU_imodel=1 ! GNBC
          ZEYU_ifgnbc=1 ! GNBC
          ZEYU_lambda=8.0D-7  ! slip length
          ZEYU_lambda=LOW%dxmin  ! slip length
          ZEYU_l_macro=LOW%dxmin
          ZEYU_l_micro=1.0D-9
          ZEYU_dgrid=LOW%dxmin 
          ZEYU_d_closest=abs(dist_to_CL)

          if (fort_denconst(im_fluid1).ge. &
              fort_denconst(im_fluid2)) then
           im_liquid=im_fluid1
           im_vapor=im_fluid2
           ZEYU_thet_s=angle_im1  ! thet_s in the liquid.
          else if (fort_denconst(im_fluid2).ge. &
                   fort_denconst(im_fluid1)) then
           im_liquid=im_fluid2
           im_vapor=im_fluid1
           ZEYU_thet_s=Pi-angle_im1
          else
           print *,"fort_denconst bust"
           stop
          endif
          ZEYU_mu_l=fort_viscconst(im_liquid)
          ZEYU_mu_g=fort_viscconst(im_vapor)
          ZEYU_sigma=user_tension(iten)
          if (im_primary_image.eq.im_liquid) then
           ZEYU_thet_d_apparent=angle_ACT
          else if (im_primary_image.eq.im_vapor) then
           ZEYU_thet_d_apparent=Pi-angle_ACT
          else
           print *,"im_primary_image or im_vapor invalid"
           stop
          endif
          ZEYU_u_cl=zero
          ZEYU_u_slip=zero
          ZEYU_thet_d=ZEYU_thet_d_apparent
           ! in: PROB.F90
          call dynamic_contact_angle(ZEYU_mu_l, ZEYU_mu_g, ZEYU_sigma, &
           ZEYU_thet_s, &
           ZEYU_imodel, ZEYU_ifgnbc, ZEYU_lambda, &
           ZEYU_l_macro, ZEYU_l_micro, &
           ZEYU_dgrid, ZEYU_d_closest, ZEYU_thet_d_apparent, &
           ZEYU_u_cl, ZEYU_u_slip, ZEYU_thet_d)

          ! nCL is normal to the contact line in the substrate plane.
          ! nCL points to the im_primary_image material.
          ! ZEYU_u_slip is positive if the contact line is advancing into
          ! the gas.
          ! NOTE: if the velocity, ZEYU_u_slip=0.0, the interface might
          ! still move slowly since the curvature will not be numerically
          ! a constant on the interface.  This is what people call 
          ! "parasitic currents" when the interface moves due to surface
          ! tension, even though the curvature = constant.
          ughost_tngt=ZEYU_u_slip     ! debugging: ZEYU_u_slip  * 1000.0 ??

          nCL_dot_n_raster=zero
          do dir=1,SDIM
           nCL_dot_n_raster=nCL_dot_n_raster+nCL(dir)*LOW%n_raster(dir)
          enddo
          mag=zero
          do dir=1,SDIM
           nCL_raster(dir)=nCL(dir)-nCL_dot_n_raster*LOW%n_raster(dir)
           mag=mag+nCL_raster(dir)**2
          enddo
          mag=sqrt(mag)
          if (mag.gt.zero) then
           do dir=1,SDIM
            nCL_raster(dir)=nCL_raster(dir)/mag
           enddo
          endif
          do dir=1,SDIM
           if (im_primary_image.eq.im_liquid) then
            u_tngt(dir)=-nCL_raster(dir)
           else if (im_primary_image.eq.im_vapor) then
            u_tngt(dir)=nCL_raster(dir)
           else
            print *,"im_primary_image or im_vapor invalid"
            stop
           endif
          enddo ! dir=1..sdim
         else
          print *,"fort_denconst invalid"
          stop
         endif
         if (DEBUG_DYNAMIC_CONTACT_ANGLE.eq.1) then
          print *,"xcrossing ",xcrossing(1),xcrossing(2),xcrossing(SDIM)
          print *,"xtriple ",xtriple(1),xtriple(2),xtriple(SDIM)
          print *,"nrm_solid ",nrm_solid(1),nrm_solid(2),nrm_solid(SDIM)
          print *,"nrm_fluid ",nrm_fluid(1),nrm_fluid(2),nrm_fluid(SDIM)
          print *,"im_fluid,angle_ACT(rad,deg) ",im_fluid,angle_ACT, &
                    angle_ACT*180.0d0/Pi
          print *,"im_liquid,ZEYU_thet_d_apparent(rad,deg) ",im_liquid, &
               ZEYU_thet_d_apparent,ZEYU_thet_d_apparent*180.0d0/Pi
          print *,"dx(1),dist_to_CL ",LOW%dx(1),dist_to_CL
          print *,"im_primary_image,im_secondary_image ", &
               im_primary_image,im_secondary_image
          print *," ZEYU_thet_s(rad,deg) ", &
           ZEYU_thet_s,ZEYU_thet_s*180.0d0/Pi
          print *,"u_tngt ",u_tngt(1),u_tngt(2),u_tngt(SDIM)
          print *,"ughost_tngt=",ughost_tngt
         endif

        else if (near_contact_line.eq.0) then
          ! ghost velocity lives *on* the rasterized interface.
         ughost_tngt = zero
        else
         print *,"near_contact_line invalid"
         stop
        endif

        if (debug_slip_velocity_enforcement.eq.1) then
         ughost_tngt = ten
         u_tngt(1)=one
         u_tngt(2)=zero
        endif

       else
        print *,"law_of_the_wall invalid"
        stop
       endif

        !get ghost velocity using normal and tangential components
        !default:
        !ughost dot n = ughost_nrml + 
        !               ughost_tngt * u_tngt dot n +
        !               usolid dot n=
        ! -(uimage-usolid) dot n +
        ! -(I-P)(uimage-usolid)dot n+
        ! usolid dot n = (2 usolid - uimage) dot n
        ! (I-P)ughost=(I-P)usolid+
        ! -(I-P)(uimage-usolid)=(I-P)(2 usolid - uimage)
       do dir=1,SDIM
        usolid_law_of_wall(dir) = &
                ughost_nrml*LOW%n_raster(dir)+ &
                ughost_tngt*u_tngt(dir)+ &
                LOW%usolid_raster(dir)
       enddo
      
       deallocate(user_tension)
 
      end subroutine getGhostVel

        ! FUTURE: Storing both the displacement and tensor at the
        !  cell centers (colocated) results in a large stencil for 
        !  the elastic force.  THIS MUST BE CHANGED:
        !    a) store displacement at cell centers.
        !    b) store tensor on the MAC grid
        !    c) store force at cell centers
        !    
        ! CURRENT STENCIL WIDTH: 5x5x5 stencil
        ! FUTURE STENCIL WIDTH: 3x3x3 stencil
        !
        ! u_t + u dot grad u = -grad p/density + div tau
        ! div u = 0
        ! phi_t + u dot grad phi=0   
        ! phi(t,x) > 0 in the deformable solid
        ! phi(t,x) < 0 in the fluid
        ! density=density_fluid * (1-H(phi)) + 
        !         density_solid * H(phi)
        ! H(phi)=1  if phi>0
        !       =0  if phi<0
        ! tau=tau_fluid * (1-H) + tau_solid * H
        ! tau_fluid=2 * mu (grad u + (grad u)^T)/2  mu=dynamic viscosity
        ! tau_solid=
        ! F=deformation gradient=I + grad XD 
        ! note: Michael Lai and David Rubin "intro to continuum mech"
        ! "u equiv XD"  "v equiv u"
        ! C=F^T F right Cauchy Green Tensor
        ! B=F F^T left Cauchy Green Tensor
        ! I1=tr(C),I2=(1/2)Tr(C)^2−tr(C^2),and I3=det(C)
      subroutine local_tensor_from_xdisplace( &
        LS_or_VOF_flag, & ! =0 => LS   =1 => VOF
        im_elastic, & ! 1..nmat
        tilelo,tilehi, &  ! tile box dimensions
        fablo,fabhi, &    ! fortran array box dimensions containing the tile
        bfact, &          ! space order
        level, &          ! 0<=level<=finest_level
        finest_level, &
        xlo,dx, &         ! xlo is lower left hand corner coordinate of fab
        ncomp_tensor, &  ! ncomp_tensor=4 in 2D (11,12,22,33) and 6 in 3D 
        nmat, &
        LS, &
        DIMS(LS), &
        VOF, &
        DIMS(VOF), &
        TNEWfab, &       ! FAB that holds elastic tensor, Q, when complete
        DIMS(TNEWfab), &
        XDISP_fab, &      
        DIMS(XDISP_fab)) 

      use global_utility_module
      use probcommon_module
      implicit none

      INTEGER_T, intent(in) :: LS_or_VOF_flag ! =0 LS   =1 VOF
      INTEGER_T, intent(in) :: im_elastic ! 1..nmat
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: ncomp_tensor
      INTEGER_T, intent(in), target :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in), target :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      REAL_T, intent(in), target :: xlo(SDIM)
      REAL_T, intent(in), target :: dx(SDIM)
      INTEGER_T, intent(in) :: DIMDEC(LS) 
      INTEGER_T, intent(in) :: DIMDEC(VOF) 
      INTEGER_T, intent(in) :: DIMDEC(TNEWfab) 
      INTEGER_T, intent(in) :: DIMDEC(XDISP_fab) 
      REAL_T, intent(in) :: LS( &  
        DIMV(LS), &
        nmat*(1+SDIM))
      REAL_T, intent(in) :: VOF( &  
        DIMV(VOF), &
        nmat*ngeom_recon)
      REAL_T, intent(inout) :: TNEWfab( &  ! Q assimilated from particles/cells
        DIMV(TNEWfab), &
        ncomp_tensor)
      REAL_T, intent(in) :: XDISP_fab( &
        DIMV(XDISP_fab), &
        SDIM)

      INTEGER_T growlo(3)
      INTEGER_T growhi(3)
      INTEGER_T i,j,k
      INTEGER_T ibase
      INTEGER_T dir_x,dir_space
      INTEGER_T iii,jjj,kkk
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf

      REAL_T dxminus,dxplus
      REAL_T grad_cen,grad_plus,grad_minus
      REAL_T LS_plus,LS_minus
      REAL_T gradu(SDIM,SDIM)  ! dir_x (displace), dir_space
      REAL_T DISP_TEN(SDIM,SDIM) ! dir_x (displace), dir_space
      REAL_T xdisplace_local,ydisplace_local
      REAL_T hoop_22
      INTEGER_T vofcomp
      REAL_T x_stress(SDIM)

      nhalf=3

      if (nmat.eq.num_materials) then
       ! do nothing
      else
       print *,"nmat invalid"
       stop
      endif

      if ((im_elastic.ge.1).and.(im_elastic.le.nmat)) then
       ! do nothing
      else
       print *,"im_elastic invalid"
       stop
      endif

       ! 6 in 3D, 4 in 2D
      if (ncomp_tensor.eq.2*SDIM) then
       ! do nothing
      else
       print *,"ncomp_tensor invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(VOF),1,-1,1271)
      call checkbound(fablo,fabhi,DIMS(LS),1,-1,1271)
      call checkbound(fablo,fabhi,DIMS(TNEWfab),1,-1,1271)
      call checkbound(fablo,fabhi,DIMS(XDISP_fab),1,-1,1271)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

! grad u=| u_r  u_t/r-v/r  u_z  |
!        | v_r  v_t/r+u/r  v_z  |
!        | w_r  w_t/r      w_z  |
! in RZ:  T33 gets u/r=x_displace/r
! in RTZ: T12=u_t/r - v/r
!         T22=v_t/r + u/r
! later:
! div S = | (r S_11)_r/r + (S_12)_t/r - S_22/r  + (S_13)_z |
!         | (r S_21)_r/r + (S_22)_t/r + S_12/r  + (S_23)_z |
!         | (r S_31)_r/r + (S_32)_t/r +           (S_33)_z |

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridsten_level(xsten,i,j,k,level,nhalf)
       do dir_space=1,SDIM
        x_stress(dir_space)=xsten(0,dir_space)
       enddo

       do dir_x=1,SDIM ! velocity (displacement) component u,v,w
       do dir_space=1,SDIM ! direction x,y,z

        iii=0
        jjj=0
        kkk=0
        if (dir_space.eq.1) then
         iii=1
        else if (dir_space.eq.2) then
         jjj=1
        else if ((dir_space.eq.3).and.(SDIM.eq.3)) then
         kkk=1
        else
         print *,"dir_space invalid"
         stop
        endif 
        dxplus=xsten(2,dir_space)-xsten(0,dir_space)
        dxminus=xsten(0,dir_space)-xsten(-2,dir_space)
        if ((dxplus.gt.zero).and.(dxminus.gt.zero)) then
         grad_plus=(XDISP_fab(D_DECL(i+iii,j+jjj,k+kkk),dir_x)- &
                    XDISP_fab(D_DECL(i,j,k),dir_x))/dxplus
         grad_minus=(XDISP_fab(D_DECL(i,j,k),dir_x)- &
                     XDISP_fab(D_DECL(i-iii,j-jjj,k-kkk),dir_x))/dxminus

         if (LS_or_VOF_flag.eq.0) then
          LS_plus=LS(D_DECL(i+iii,j+jjj,k+kkk),im_elastic)
          LS_minus=LS(D_DECL(i-iii,j-jjj,k-kkk),im_elastic)
         else if (LS_or_VOF_flag.eq.1) then
          vofcomp=(im_elastic-1)*ngeom_recon+1
          LS_plus=VOF(D_DECL(i+iii,j+jjj,k+kkk),vofcomp)-half
          LS_minus=VOF(D_DECL(i-iii,j-jjj,k-kkk),vofcomp)-half
         else
          print *,"LS_or_VOF_flag invalid"
          stop
         endif

         if ((LS_plus.ge.zero).and.(LS_minus.ge.zero)) then
          grad_cen=half*(grad_plus+grad_minus)
         else if ((LS_plus.lt.zero).and.(LS_minus.lt.zero)) then
          grad_cen=half*(grad_plus+grad_minus)
         else if ((LS_plus.ge.zero).and.(LS_minus.lt.zero)) then
          grad_cen=grad_plus
         else if ((LS_plus.lt.zero).and.(LS_minus.ge.zero)) then
          grad_cen=grad_minus
         else
          print *,"LS_plus or LS_minus invalid"
          stop
         endif

         gradu(dir_x,dir_space)=grad_cen  ! dir_x (displace), dir_space
        else
         print *,"dxplus or dxminus invalid"
         stop
        endif
       enddo ! dir_space
       enddo ! dir_x (displace)

       xdisplace_local=XDISP_fab(D_DECL(i,j,k),1)
       ydisplace_local=XDISP_fab(D_DECL(i,j,k),2)

        ! declared in GLOBALUTIL.F90
       call stress_from_strain( &
         im_elastic, &
         x_stress, &
         dx, &
         gradu, &  ! dir_x (displace),dir_space
         xdisplace_local, &
         ydisplace_local, &
         DISP_TEN, &  ! dir_x (displace),dir_space
         hoop_22) ! output: "theta-theta" component (xdisp/r)

       ibase=1
       TNEWfab(D_DECL(i,j,k),ibase)=DISP_TEN(1,1)
       ibase=ibase+1
       TNEWfab(D_DECL(i,j,k),ibase)=DISP_TEN(1,2)
       ibase=ibase+1
       TNEWfab(D_DECL(i,j,k),ibase)=DISP_TEN(2,2)

       ibase=ibase+1
       if (SDIM.eq.3) then
        TNEWfab(D_DECL(i,j,k),ibase)=DISP_TEN(SDIM,SDIM)
       else if (SDIM.eq.2) then
        if (levelrz.eq.0) then
          ! T33 (theta coordinate)
         TNEWfab(D_DECL(i,j,k),ibase)=zero
        else if (levelrz.eq.1) then
          ! T33 (theta coordinate)
          ! dX/dx + dX/dx
         TNEWfab(D_DECL(i,j,k),ibase)=two*hoop_22 ! 2 * (xdisp/r)
        else if (levelrz.eq.3) then
          ! T33 (z coordinate)
         TNEWfab(D_DECL(i,j,k),ibase)=zero
        else
         print *,"levelrz invalid"
         stop
        endif
       else
        print *,"dimension bust"
        stop
       endif
                    
       if (SDIM.eq.3) then                
        ibase=ibase+1
        TNEWfab(D_DECL(i,j,k),ibase)=DISP_TEN(1,SDIM)
        ibase=ibase+1
        TNEWfab(D_DECL(i,j,k),ibase)=DISP_TEN(2,SDIM)
       endif
       if (ibase.eq.2*SDIM) then
        ! do nothing
       else
        print *,"ibase invalid (7) ibase=",ibase
        stop
       endif

      enddo
      enddo
      enddo

      end subroutine local_tensor_from_xdisplace


      end module godunov_module


      ! enable_spectral:
      ! 0 - low order
      ! 1 - space/time spectral
      ! 2 - space spectral only
      ! 3 - time spectral only
      subroutine FORT_BUILD_MASKSEM( &
       spectral_cells_level, &
       mask_sweep, &
       level, &
       finest_level, &
       cur_time, &
       enable_spectral, &
       domlo,domhi, &
       vofbc, &
       maskcov,DIMS(maskcov), &
       masknbr,DIMS(masknbr), &
       mask,DIMS(mask), &
       oldmask,DIMS(oldmask), &
       vfrac,DIMS(vfrac), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       bfact_fine, &
       nmat)
      use probcommon_module
      use global_utility_module
      use MOF_routines_module
      use probf90_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: nmat
      REAL_T, intent(inout) :: spectral_cells_level(nmat)
      INTEGER_T, intent(in) :: mask_sweep
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      REAL_T, intent(in) :: cur_time
      INTEGER_T, intent(in) :: enable_spectral
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(in) :: domlo(SDIM),domhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: bfact_fine
      INTEGER_T, intent(in) :: vofbc(SDIM,2)
      INTEGER_T, intent(in) :: DIMDEC(maskcov) 
      INTEGER_T, intent(in) :: DIMDEC(masknbr) 
      INTEGER_T, intent(in) :: DIMDEC(mask) 
      INTEGER_T, intent(in) :: DIMDEC(oldmask) 
      INTEGER_T, intent(in) :: DIMDEC(vfrac) 
      REAL_T, intent(in) :: maskcov(DIMV(maskcov))
      REAL_T, intent(in) :: masknbr(DIMV(masknbr))
      REAL_T, intent(out) :: mask(DIMV(mask))
      REAL_T, intent(in) :: oldmask(DIMV(oldmask))
      REAL_T, intent(in) :: vfrac(DIMV(vfrac),nmat)

      INTEGER_T im,imcrit,im_max
      INTEGER_T sumtag
      INTEGER_T tag(nmat)
      INTEGER_T local_maskSEM
      INTEGER_T old_maskSEM
      INTEGER_T i,j,k
      INTEGER_T inbr,jnbr,knbr
      INTEGER_T inormal
      INTEGER_T dir,side
      INTEGER_T iofs,jofs,kofs,kofs_hi
      INTEGER_T stripstat
      INTEGER_T test_maskcov
      INTEGER_T test_masknbr
      INTEGER_T test_maskSEM
      INTEGER_T covered_count
      INTEGER_T uncovered_count
      REAL_T vfracmax
      REAL_T vfractest(nmat)
      REAL_T vfrac_fluid
      REAL_T vfrac_sum_fluid,vfrac_sum_solid
      INTEGER_T touch_coarsefine
      INTEGER_T localbc
      INTEGER_T clamped_cell_in_element
      REAL_T xclamped(SDIM)
      REAL_T LS_clamped
      REAL_T vel_clamped(SDIM)
      REAL_T temperature_clamped
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf

      nhalf=3

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif

      if (cur_time.ge.zero) then
       ! do nothing
      else
       print *,"cur_time invalid"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid build masksem"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid43"
       stop
      endif
      if ((bfact_fine.gt.bfact).or.(bfact_fine.lt.1)) then
       print *,"bfact_fine invalid"
       stop
      endif

      if ((enable_spectral.lt.0).or. &
          (enable_spectral.gt.3)) then
       print *,"enable_spectral invalid build masksem"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(maskcov),1,-1,1271)
      call checkbound(fablo,fabhi,DIMS(masknbr),1,-1,1271)
      call checkbound(fablo,fabhi,DIMS(mask),1,-1,1271)
      call checkbound(fablo,fabhi,DIMS(oldmask),1,-1,1271)
      call checkbound(fablo,fabhi,DIMS(vfrac),0,-1,1272)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call strip_status(i,j,k,bfact,stripstat)

       if (stripstat.eq.1) then

        do im=1,nmat
         tag(im)=0
        enddo
        imcrit=0

        covered_count=0
        uncovered_count=0

        im_max=0
        vfracmax=zero

        if (SDIM.eq.2) then
         kofs_hi=0
        else if (SDIM.eq.3) then
         kofs_hi=bfact-1
        else
         print *,"dimension bust"
         stop
        endif

        clamped_cell_in_element=0

        do iofs=0,bfact-1
        do jofs=0,bfact-1
        do kofs=0,kofs_hi

         call gridsten_level(xsten,i+iofs,j+jofs,k+kofs,level,nhalf)
         do dir=1,SDIM
          xclamped(dir)=xsten(0,dir)
         enddo
          ! LS>0 if clamped
         call SUB_clamped_LS(xclamped,cur_time,LS_clamped, &
             vel_clamped,temperature_clamped)

         if (LS_clamped.ge.zero) then
          clamped_cell_in_element=1
         else if (LS_clamped.lt.zero) then
          ! do nothing
         else
          print *,"LS_clamped is NaN"
          stop
         endif

         vfrac_sum_fluid=zero
         vfrac_sum_solid=zero

         do im=1,nmat
          vfractest(im)=vfrac(D_DECL(i+iofs,j+jofs,k+kofs),im)
          if (is_rigid(nmat,im).eq.1) then
           vfrac_sum_solid=vfrac_sum_solid+vfractest(im)

           if (vfractest(im).gt.vfracmax) then
            im_max=im
            vfracmax=vfractest(im)
           endif
           if (vfractest(im).gt.VOFTOL) then
            imcrit=im
            tag(im)=1
           endif

          else if (is_rigid(nmat,im).eq.0) then

           vfrac_sum_fluid=vfrac_sum_fluid+vfractest(im)

          else
           print *,"is_rigid(nmat,im) invalid"
           stop
          endif
         enddo ! im=1..nmat

         if ((vfrac_sum_solid.lt.-VOFTOL).or. &
             (vfrac_sum_solid.gt.one+VOFTOL)) then
          print *,"vfrac_sum_solid out of range"
          stop
         else if ((vfrac_sum_solid.ge.-VOFTOL).and. &
                  (vfrac_sum_solid.le.one+VOFTOL)) then
          ! do nothing
         else
          print *,"vfrac_sum_solid is NaN"
          stop
         endif

         if (abs(one-vfrac_sum_fluid).gt.0.01) then
          print *,"vfrac_sum_fluid out of range in build mask sem"
          print *,"vfrac_sum_solid=",vfrac_sum_solid
          print *,"vfrac_sum_fluid=",vfrac_sum_fluid
          print *,"i,j,k ",i,j,k
          print *,"iofs,jofs,kofs ",iofs,jofs,kofs
          print *,"level,finest_level ",level,finest_level
          stop
         endif
 
         do im=1,nmat
          if (is_rigid(nmat,im).eq.1) then
           ! do nothing
          else if (is_rigid(nmat,im).eq.0) then
           vfrac_fluid=(one-vfrac_sum_solid)*vfractest(im)
           if (vfrac_fluid.gt.vfracmax) then
            im_max=im
            vfracmax=vfrac_fluid
           endif  
           if (vfrac_fluid.gt.VOFTOL) then
            imcrit=im
            tag(im)=1
           endif
          else
           print *,"is_rigid(nmat,im) invalid"
           stop
          endif
         enddo ! im=1..nmat

         test_maskcov=NINT(maskcov(D_DECL(i+iofs,j+jofs,k+kofs))) 
         if (test_maskcov.eq.1) then
          uncovered_count=uncovered_count+1
         else if (test_maskcov.eq.0) then
          covered_count=covered_count+1
         else
          print *,"test_maskcov invalid"
          stop
         endif
           
        enddo
        enddo
        enddo ! iofs,jofs,kofs

        if ((imcrit.lt.1).or.(imcrit.gt.nmat)) then
         print *,"imcrit invalid"
         stop
        endif
        if ((im_max.lt.1).or.(im_max.gt.nmat)) then
         print *,"im_max invalid"
         stop
        endif
        if (vfracmax.le.VOFTOL) then
         print *,"vfracmax invalid"
         stop
        endif

        sumtag=0
        do im=1,nmat
         sumtag=sumtag+tag(im)
        enddo

        if ((sumtag.le.0).or.(sumtag.gt.nmat)) then
         print *,"sumtag invalid"
         stop
        else if (sumtag.eq.1) then
         if (clamped_cell_in_element.eq.1) then
          local_maskSEM=0
         else if (clamped_cell_in_element.eq.0) then
          if (is_rigid(nmat,imcrit).eq.1) then
           local_maskSEM=0
          else if (is_ice(nmat,imcrit).eq.1) then
           local_maskSEM=0
          else if (is_FSI_rigid(nmat,imcrit).eq.1) then
           local_maskSEM=0
          else if ((FSI_flag(imcrit).eq.0).or. &
                   (FSI_flag(imcrit).eq.7)) then
           local_maskSEM=imcrit
          else
           print *,"FSI_flag invalid"
           stop
          endif
         else
          print *,"clamped_cell_in_element invalid"
          stop
         endif

        else if (sumtag.gt.1) then
         local_maskSEM=0
        else
         print *,"sumtag bust"
         stop
        endif

        if (covered_count.gt.0) then
         if (uncovered_count.ne.0) then
          print *,"cannot have an element partially covered"
          stop
         endif
        else if (uncovered_count.gt.0) then
         if (covered_count.ne.0) then
          print *,"cannot have an element partially covered"
          stop
         endif
        else
         print *,"covered_count or uncovered_count invalid"
         stop
        endif 

        test_maskcov=NINT(maskcov(D_DECL(i,j,k))) 
        if ((test_maskcov.eq.0).or.(test_maskcov.eq.1)) then
         ! do nothing
        else
         print *,"test_maskcov invalid"
         stop
        endif

        if (mask_sweep.eq.0) then
         ! do nothing
        else if (mask_sweep.eq.1) then ! check neighboring elements

         old_maskSEM=NINT(oldmask(D_DECL(i,j,k))) 

         if (old_maskSEM.ne.local_maskSEM) then
          print *,"old_maskSEM.ne.local_maskSEM"
          stop
         endif

         if (test_maskcov.eq.1) then
          ! do nothing

          ! high order stencil never includes covered values.
         else if (test_maskcov.eq.0) then 
          local_maskSEM=0
         else
          print *,"test_maskcov invalid"
          stop
         endif

         if ((local_maskSEM.ge.1).and.(local_maskSEM.le.nmat)) then

          if (SDIM.eq.2) then
           kofs_hi=-1
          else if (SDIM.eq.3) then
           kofs_hi=bfact
          else
           print *,"dimension bust"
           stop
          endif

          do iofs=-1,bfact
          do jofs=-1,bfact
          do kofs=-1,kofs_hi
           test_maskSEM=NINT(oldmask(D_DECL(i+iofs,j+jofs,k+kofs))) 
           if (test_maskSEM.ne.local_maskSEM) then
            local_maskSEM=0
           endif
          enddo
          enddo
          enddo ! iofs,jofs,kofs

         else if (local_maskSEM.eq.0) then
          ! do nothing
         else
          print *,"local_maskSEM invalid"
          stop
         endif

        else
         print *,"mask_sweep invalid"
         stop
        endif


        if (SDIM.eq.2) then
         kofs_hi=0
        else if (SDIM.eq.3) then
         kofs_hi=bfact-1
        else
         print *,"dimension bust"
         stop
        endif

        do iofs=0,bfact-1
        do jofs=0,bfact-1
        do kofs=0,kofs_hi

         mask(D_DECL(i+iofs,j+jofs,k+kofs))=local_maskSEM

         if ((local_maskSEM.ge.1).and.(local_maskSEM.le.nmat)) then
          spectral_cells_level(local_maskSEM)= &
            spectral_cells_level(local_maskSEM)+one
         else if (local_maskSEM.eq.0) then
          ! do nothing
         else
          print *,"local maskSEM invalid"
          stop
         endif
  
        enddo
        enddo
        enddo ! iofs,jofs,kofs

         ! if order (bfact) > 1,
         ! 1. make sure no finest grid elements that neighbor uncovered
         !    coarse grid elements are tagged as CISL-MOF.
         ! 2. make sure no uncovered coarse grid elements are tagged as
         !    CISL-MOF.
        if (bfact.eq.1) then
         ! do nothing
        else if ((bfact.ge.2).and.(bfact.le.16)) then
         if (mask_sweep.eq.0) then
          ! do nothing
         else if (mask_sweep.eq.1) then
          if ((local_maskSEM.ge.1).and.(local_maskSEM.le.nmat)) then
           ! do nothing
          else if (local_maskSEM.eq.0) then
           if (level.eq.finest_level) then
            if (test_maskcov.eq.1) then

             touch_coarsefine=0

             do dir=1,SDIM
              do side=1,2

               inbr=i
               jnbr=j
               knbr=k

               localbc=vofbc(dir,side)

               if (side.eq.1) then
                if (dir.eq.1) then
                 inormal=i
                 inbr=inormal-1
                else if (dir.eq.2) then 
                 inormal=j
                 jnbr=inormal-1
                else if ((dir.eq.3).and.(SDIM.eq.3)) then 
                 inormal=k
                 knbr=inormal-1
                else
                 print *,"dir invalid"
                 stop
                endif

                 ! =1 fine-fine ghost =0 otherwise
                test_masknbr=NINT(masknbr(D_DECL(inbr,jnbr,knbr)))

                if (inormal.eq.fablo(dir)) then
                 if (localbc.eq.INT_DIR) then
                  if (test_masknbr.eq.1) then
                   ! do nothing (fine-fine)
                  else if (test_masknbr.eq.0) then
                   touch_coarsefine=1
                  else
                   print *,"test_masknbr invalid"
                   stop
                  endif
                 else if ((localbc.eq.EXT_DIR).or. &
                          (localbc.eq.FOEXTRAP).or. &
                          (localbc.eq.REFLECT_EVEN).or. &
                          (localbc.eq.REFLECT_ODD)) then
                  if (inormal.eq.domlo(dir)) then
                   ! do nothing
                  else
                   print *,"inormal invalid"
                   stop
                  endif
                 else
                  print *,"localbc invalid"
                  stop
                 endif
                else if ((inormal.gt.fablo(dir)).and. &
                         (inormal.lt.fabhi(dir))) then
                 ! do nothing
                else
                 print *,"inormal invalid"
                 stop
                endif
               else if (side.eq.2) then
                if (dir.eq.1) then
                 inormal=i+bfact-1
                 inbr=inormal+1
                else if (dir.eq.2) then 
                 inormal=j+bfact-1
                 jnbr=inormal+1
                else if ((dir.eq.3).and.(SDIM.eq.3)) then 
                 inormal=k+bfact-1
                 knbr=inormal+1
                else
                 print *,"dir invalid"
                 stop
                endif

                 ! =1 fine-fine ghost =0 otherwise
                test_masknbr=NINT(masknbr(D_DECL(inbr,jnbr,knbr)))

                if (inormal.eq.fabhi(dir)) then
                 if (localbc.eq.INT_DIR) then
                  if (test_masknbr.eq.1) then
                   ! do nothing (fine-fine)
                  else if (test_masknbr.eq.0) then
                   touch_coarsefine=1
                  else
                   print *,"test_masknbr invalid"
                   stop
                  endif
                 else if ((localbc.eq.EXT_DIR).or. &
                          (localbc.eq.FOEXTRAP).or. &
                          (localbc.eq.REFLECT_EVEN).or. &
                          (localbc.eq.REFLECT_ODD)) then
                  if (inormal.eq.domhi(dir)) then
                   ! do nothing
                  else
                   print *,"inormal invalid"
                   stop
                  endif
                 else
                  print *,"localbc invalid"
                  stop
                 endif
                else if ((inormal.gt.fablo(dir)).and. &
                         (inormal.lt.fabhi(dir))) then
                 ! do nothing
                else
                 print *,"inormal invalid"
                 stop
                endif
               else
                print *,"side invalid"
                stop
               endif
              enddo ! side
             enddo ! dir

             if (touch_coarsefine.eq.1) then
              print *,"cant mask finest lev element next to coarse element"
              print *,"i,j,k = ",i,j,k
              print *,"level,finest_level ",level,finest_level
              print *,"increase n_error_buf"
              stop
             else if (touch_coarsefine.eq.0) then
              ! do nothing
             else
              print *,"touch_coarsefine invalid"
              stop
             endif
            else 
             print *,"test_maskcov invalid"
             stop
            endif
           else if ((level.ge.0).and.(level.lt.finest_level)) then
            if (test_maskcov.eq.1) then
             print *,"cannot mask a coarse element that has bfact>1."
             print *,"i,j,k ",i,j,k
             print *,"level=",level
             print *,"finest_level=",finest_level
             print *,"bfact=",bfact
             print *,"increase n_error_buf"
             stop
            else if (test_maskcov.eq.0) then
             ! do nothing
            else
             print *,"test_maskcov invalid"
             stop
            endif
           else
            print *,"level invalid"
            stop
           endif
          else
           print *,"local_maskSEM invalid"
           stop
          endif
         else
          print *,"mask_sweep invalid"
          stop
         endif

        else
         print *,"bfact out of range"
         stop
        endif

       else if (stripstat.eq.0) then
        ! do nothing
       else
        print *,"stripstat invalid"
        stop
       endif

      enddo
      enddo
      enddo ! i,j,k (only interior cells)

      return
      end subroutine FORT_BUILD_MASKSEM


      subroutine FORT_BUILD_CONSERVE( &
       iden_base, &
       override_density, &
       constant_density_all_time, &
       temperature_primitive_variable, &
       conserve,DIMS(conserve), &
       den, &
       DIMS(den), &
       mom_den, &
       DIMS(mom_den), &
       vel,DIMS(vel), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       nmat,ngrow, &
       normdir, &
       nc_conserve, &
       nc_den)
      use probf90_module
      use global_utility_module
      use MOF_routines_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: nmat,ngrow
      INTEGER_T, intent(in) :: normdir
      INTEGER_T, intent(in) :: nc_conserve
      INTEGER_T, intent(in) :: nc_den
      INTEGER_T, intent(in) :: temperature_primitive_variable(nmat) 
      INTEGER_T, intent(in) :: override_density(nmat)
      INTEGER_T, intent(in) :: constant_density_all_time(nmat)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: DIMDEC(conserve) 
      INTEGER_T, intent(in) :: DIMDEC(den) 
      INTEGER_T, intent(in) :: DIMDEC(mom_den) 
      INTEGER_T, intent(in) :: DIMDEC(vel) 
      REAL_T, intent(out) :: conserve(DIMV(conserve),nc_conserve)
      REAL_T, intent(in) :: den(DIMV(den),nc_den)
      REAL_T, intent(in) :: mom_den(DIMV(mom_den),nmat)
      REAL_T, intent(in) :: vel(DIMV(den),SDIM)

      INTEGER_T i,j,k,im
      INTEGER_T istate,ispecies
      INTEGER_T dencomp,tempcomp,speccomp
      INTEGER_T veldir
      INTEGER_T iden_base
      INTEGER_T igridlo(3),igridhi(3)
      REAL_T dencore(nmat)
      REAL_T mom_dencore(nmat)
      REAL_T KE,vel1D,local_temperature,local_internal
      REAL_T :: massfrac_parm(num_species_var+1)

      if (nc_den.ne.num_state_material*nmat) then
       print *,"nc_den invalid"
       stop
      endif
      if (nc_conserve.ne.SDIM+nmat*num_state_material) then
       print *,"nc_conserve invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid44"
       stop
      endif
      if (ngrow.lt.1) then
       print *,"ngrow out of range in BUILD_CONSERVE ngrow=",ngrow
       stop
      endif
      if ((normdir.ge.0).and.(normdir.lt.SDIM)) then
       ! do nothing 
      else
       print *,"normdir invalid"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      do im=1,nmat
       if ((fort_material_type(im).eq.0).or. &
           (is_rigid(nmat,im).eq.1).or. &
           (fort_material_type(im).eq.999)) then
        if (temperature_primitive_variable(im).ne.1) then
         print *,"temperature_primitive_variable(im) invalid"
         stop
        endif
       else if ((fort_material_type(im).gt.0).and. &
                (is_rigid(nmat,im).eq.0).and. &
                (fort_material_type(im).ne.999)) then
        if ((temperature_primitive_variable(im).eq.0).or. &
            (temperature_primitive_variable(im).eq.1)) then
         ! do nothing
        else
         print *,"temperature_primitive_variable(im) invalid"
         stop
        endif
       else
        print *,"fort_material_type(im) or is_rigid invalid"
        stop
       endif
      enddo ! im=1..nmat

      call checkbound(fablo,fabhi,DIMS(conserve),ngrow,-1,1271)
      call checkbound(fablo,fabhi,DIMS(den),ngrow,-1,1272)
      call checkbound(fablo,fabhi,DIMS(mom_den),ngrow,-1,1272)
      call checkbound(fablo,fabhi,DIMS(vel),ngrow,-1,1272)

      call growntilebox(tilelo,tilehi,fablo,fabhi, &
        igridlo,igridhi,ngrow)

      if (iden_base.ne.SDIM) then
       print *,"iden_base invalid"
       stop
      endif

      do i=igridlo(1),igridhi(1)
      do j=igridlo(2),igridhi(2)
      do k=igridlo(3),igridhi(3)

        ! KE=u dot u/2
       KE=zero
       do veldir=1,SDIM
        vel1D=vel(D_DECL(i,j,k),veldir)
        conserve(D_DECL(i,j,k),veldir)=vel1D
        KE=KE+vel1D**2
       enddo
       KE=half*KE

       do im=1,nmat
        istate=1
        dencomp=(im-1)*num_state_material+istate
        dencore(im)=den(D_DECL(i,j,k),dencomp)
        mom_dencore(im)=mom_den(D_DECL(i,j,k),im)
          ! sanity check
        if (constant_density_all_time(im).eq.1) then
         if (abs(dencore(im)-fort_denconst(im)).le.VOFTOL) then
          ! do nothing
         else
          print *,"dencore(im) invalid"
          print *,"im,i,j,k,den ",im,i,j,k,dencore(im)
          print *,"dencomp=",dencomp
          print *,"normdir=",normdir
          stop
         endif
        else if (constant_density_all_time(im).eq.0) then 
         ! do nothing
        else
         print *,"constant_density_all_time invalid"
         stop
        endif

        if (dencore(im).gt.zero) then
         ! do nothing
        else
         print *,"density must be positive build_conserve"
         print *,"im,dencore(im) ",im,dencore(im)
         print *,"im,fort_denconst(im) ",im,fort_denconst(im)
         stop
        endif  

        if (mom_dencore(im).gt.zero) then
         ! do nothing
        else
         print *,"mom_density must be positive build_conserve"
         print *,"im,mom_dencore(im) ",im,mom_dencore(im)
         print *,"im,fort_denconst(im) ",im,fort_denconst(im)
         stop
        endif  

       enddo ! im=1..nmat

       do im=1,nmat

         ! in: FORT_BUILD_CONSERVE
        istate=1
        do while (istate.le.num_state_material)

         if (istate.eq.1) then ! Density
          dencomp=(im-1)*num_state_material+istate
          conserve(D_DECL(i,j,k),iden_base+dencomp)=dencore(im)
          istate=istate+1
         else if (istate.eq.2) then ! Temperature
          tempcomp=(im-1)*num_state_material+istate
          local_temperature=den(D_DECL(i,j,k),tempcomp)
          if (temperature_primitive_variable(im).eq.1) then ! non-conservative
            ! den * T
           conserve(D_DECL(i,j,k),iden_base+tempcomp)= &
            dencore(im)*local_temperature
          else if (temperature_primitive_variable(im).eq.0) then ! conservative

           call init_massfrac_parm(dencore(im),massfrac_parm,im)
           do ispecies=1,num_species_var
            massfrac_parm(ispecies)=den(D_DECL(i,j,k),tempcomp+ispecies)
            if (massfrac_parm(ispecies).ge.zero) then
             ! do nothing
            else
             print *,"massfrac_parm(ispecies) invalid"
             stop
            endif
           enddo

            ! den * (u dot u/2 + cv T)
           call INTERNAL_material(dencore(im),massfrac_parm, &
            local_temperature,local_internal, &
            fort_material_type(im),im)
           conserve(D_DECL(i,j,k),iden_base+tempcomp)= &
             dencore(im)*(KE+local_internal) 
          else
           print *,"temperature_primitive_variable invalid"
           stop
          endif
          istate=istate+1
         else if ((istate.eq.num_state_base+1).and. &
                  (num_species_var.gt.0)) then 
           ! den * Y
          do ispecies=1,num_species_var
           speccomp=(im-1)*num_state_material+num_state_base+ispecies
           conserve(D_DECL(i,j,k),iden_base+speccomp)= &
             dencore(im)*den(D_DECL(i,j,k),speccomp)
           istate=istate+1
          enddo ! ispecies=1..num_species_var
         else 
          print *,"istate invalid"
          stop
         endif

        enddo ! do while (istate.le.num_state_material)

        if (dencore(im).gt.zero) then 
         ! do nothing
        else
         print *,"dencore must be positive"
         stop
        endif

       enddo ! im=1..nmat

      enddo         
      enddo         
      enddo ! i,j,k (cell center "conserved" variables) 

      return
      end subroutine FORT_BUILD_CONSERVE



      subroutine FORT_BUILD_NEWMAC( &
       num_MAC_vectors, & ! num_MAC_vectors=1 or 2
       normdir, & ! 0..sdim-1
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       unode,DIMS(unode), &
       xmomside,DIMS(xmomside), &
       ymomside,DIMS(ymomside), &
       zmomside,DIMS(zmomside), &
       xmassside,DIMS(xmassside), &
       ymassside,DIMS(ymassside), &
       zmassside,DIMS(zmassside), &
       xvmac,DIMS(xvmac), &
       yvmac,DIMS(yvmac), &
       zvmac,DIMS(zvmac), &
       xdmac,DIMS(xdmac), &
       ydmac,DIMS(ydmac), &
       zdmac,DIMS(zdmac), &
       mask,DIMS(mask), &
       xlo,dx, &
       cur_time, &
       nmat, &
       level, &
       finest_level)
      use probcommon_module
      use probf90_module
      use global_utility_module
      use MOF_routines_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: num_MAC_vectors !num_MAC_vectors=1 or 2
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: normdir !normdir=0..sdim-1
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T :: veldir
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(in) :: bfact

      INTEGER_T, intent(in) :: DIMDEC(unode) 

      INTEGER_T, intent(in) :: DIMDEC(xmomside) 
      INTEGER_T, intent(in) :: DIMDEC(ymomside) 
      INTEGER_T, intent(in) :: DIMDEC(zmomside) 

      INTEGER_T, intent(in) :: DIMDEC(xmassside) 
      INTEGER_T, intent(in) :: DIMDEC(ymassside) 
      INTEGER_T, intent(in) :: DIMDEC(zmassside) 

      INTEGER_T, intent(in) :: DIMDEC(xvmac) 
      INTEGER_T, intent(in) :: DIMDEC(yvmac) 
      INTEGER_T, intent(in) :: DIMDEC(zvmac) 

      INTEGER_T, intent(in) :: DIMDEC(xdmac) 
      INTEGER_T, intent(in) :: DIMDEC(ydmac) 
      INTEGER_T, intent(in) :: DIMDEC(zdmac) 

      INTEGER_T, intent(in) :: DIMDEC(mask) 

      REAL_T, intent(in) :: unode(DIMV(unode))

      REAL_T, intent(in) :: xmomside(DIMV(xmomside), &
        2*num_MAC_vectors)
      REAL_T, intent(in) :: ymomside(DIMV(xmomside), &
        2*num_MAC_vectors)
      REAL_T, intent(in) :: zmomside(DIMV(xmomside), &
        2*num_MAC_vectors)

      REAL_T, intent(in) :: xmassside(DIMV(xmassside), &
        2*num_MAC_vectors)
      REAL_T, intent(in) :: ymassside(DIMV(ymassside), &
        2*num_MAC_vectors)
      REAL_T, intent(in) :: zmassside(DIMV(zmassside), &
        2*num_MAC_vectors)

      REAL_T, intent(out) :: xvmac(DIMV(xvmac))
      REAL_T, intent(out) :: yvmac(DIMV(yvmac))
      REAL_T, intent(out) :: zvmac(DIMV(zvmac))

      REAL_T, intent(out) :: xdmac(DIMV(xdmac))
      REAL_T, intent(out) :: ydmac(DIMV(ydmac))
      REAL_T, intent(out) :: zdmac(DIMV(zdmac))

      REAL_T, intent(in) :: mask(DIMV(mask))

      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: cur_time

      INTEGER_T igridlo(3),igridhi(3)
      INTEGER_T iii,jjj,kkk
      INTEGER_T i,j,k
      INTEGER_T iside
      INTEGER_T icell,jcell,kcell
      INTEGER_T ibucket
      INTEGER_T ibucket_map
      INTEGER_T ileft,jleft,kleft
      INTEGER_T iright,jright,kright
      INTEGER_T ivec
      INTEGER_T zapvel
      REAL_T maskleft
      REAL_T maskright
      REAL_T momface_total(num_MAC_vectors)
      REAL_T massface_total(num_MAC_vectors)
      REAL_T massquarter,momquarter
      REAL_T xsten(-1:1,SDIM)

      REAL_T xclamped(SDIM)
      REAL_T LS_clamped
      REAL_T vel_clamped(SDIM)
      REAL_T temperature_clamped

      INTEGER_T dir_local
      INTEGER_T nhalf

      nhalf=1

      if (bfact.lt.1) then
       print *,"bfact invalid45"
       stop
      endif
      if ((normdir.ge.0).and.(normdir.lt.SDIM)) then
       ! do nothing
      else
       print *,"normdir invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid build newmac"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if ((num_MAC_vectors.eq.1).or. &
          (num_MAC_vectors.eq.2)) then
       ! do nothing
      else
       print *,"num_MAC_vectors invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(unode),0,normdir,1271)

      call checkbound(fablo,fabhi,DIMS(xmomside),1,-1,1271)
      call checkbound(fablo,fabhi,DIMS(ymomside),1,-1,1271)
      call checkbound(fablo,fabhi,DIMS(zmomside),1,-1,1271)

      call checkbound(fablo,fabhi,DIMS(xmassside),1,-1,1271)
      call checkbound(fablo,fabhi,DIMS(ymassside),1,-1,1271)
      call checkbound(fablo,fabhi,DIMS(zmassside),1,-1,1271)

      call checkbound(fablo,fabhi,DIMS(mask),1,-1,1271)

      call checkbound(fablo,fabhi,DIMS(xvmac),0,0,1271)
      call checkbound(fablo,fabhi,DIMS(yvmac),0,1,1271)
      call checkbound(fablo,fabhi,DIMS(zvmac),0,SDIM-1,1271)

      call checkbound(fablo,fabhi,DIMS(xdmac),0,0,1271)
      call checkbound(fablo,fabhi,DIMS(ydmac),0,1,1271)
      call checkbound(fablo,fabhi,DIMS(zdmac),0,SDIM-1,1271)

      do veldir=1,SDIM

        iii=0
        jjj=0
        kkk=0 

        call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
          igridlo,igridhi,0,veldir-1,19)

        if (veldir.eq.1) then
         iii=1
        else if (veldir.eq.2) then
         jjj=1
        else if ((veldir.eq.3).and.(SDIM.eq.3)) then
         kkk=1
        else
         print *,"veldir invalid"
         stop
        endif

        do i=igridlo(1),igridhi(1)
        do j=igridlo(2),igridhi(2)
        do k=igridlo(3),igridhi(3)

          ! veldir=1..sdim
         call gridstenMAC_level(xsten,i,j,k,level,nhalf,veldir-1,25)
         do dir_local=1,SDIM
          xclamped(dir_local)=xsten(0,dir_local)
         enddo

         zapvel=0
         if (levelrz.eq.0) then
          ! do nothing
         else if (levelrz.eq.1) then
          if (SDIM.ne.2) then
           print *,"dimension bust"
           stop
          endif
          if ((xsten(0,1).le.VOFTOL*dx(1)).and. &
              (veldir.eq.1)) then
           zapvel=1
          endif
         else if (levelrz.eq.3) then
          if ((xsten(0,1).le.VOFTOL*dx(1)).and. &
              (veldir.eq.1)) then
           zapvel=1
          endif
         else
          print *,"levelrz invalid"
          stop
         endif

         if (zapvel.eq.1) then
          if (veldir.eq.1) then
           xvmac(D_DECL(i,j,k))=zero
           if (num_MAC_vectors.eq.2) then
            xdmac(D_DECL(i,j,k))=zero
           else if (num_MAC_vectors.eq.1) then
            ! do nothing
           else
            print *,"num_MAC_vectors invalid"
            stop
           endif
          else
           print *,"veldir invalid"
           stop
          endif
         else if (zapvel.eq.0) then
          iright=i
          jright=j
          kright=k
          ileft=iright-iii
          jleft=jright-jjj
          kleft=kright-kkk
           ! mask=1 if cell is not covered by level+1 or cell is 
           ! outside the domain.
          maskleft=mask(D_DECL(ileft,jleft,kleft))
          maskright=mask(D_DECL(iright,jright,kright))

          if ((maskleft.eq.one).and.(maskright.eq.one)) then

            do ivec=1,num_MAC_vectors
             momface_total(ivec)=zero
             massface_total(ivec)=zero
            enddo

             ! iside=-1 (left of face)
             ! iside=1  (right of face)
            do iside=-1,1,2

              ! iside=-1 (left of face)
             if (iside.eq.-1) then
              icell=ileft
              jcell=jleft
              kcell=kleft
              ibucket=2  ! right side of cell ileft,jleft,kleft

               ! iside=1  (right of face)
             else if (iside.eq.1) then
              icell=iright
              jcell=jright
              kcell=kright
              ibucket=1  ! left side of cell iright,jright,kright
             else 
              print *,"iside invalid"
              stop
             endif

             do ivec=1,num_MAC_vectors

              ibucket_map=ibucket+2*(ivec-1)
           
              if (veldir.eq.1) then
               massquarter=xmassside(D_DECL(icell,jcell,kcell),ibucket_map)
               momquarter=xmomside(D_DECL(icell,jcell,kcell),ibucket_map)
              else if (veldir.eq.2) then
               massquarter=ymassside(D_DECL(icell,jcell,kcell),ibucket_map)
               momquarter=ymomside(D_DECL(icell,jcell,kcell),ibucket_map)
              else if ((veldir.eq.3).and.(SDIM.eq.3)) then
               massquarter=zmassside(D_DECL(icell,jcell,kcell),ibucket_map)
               momquarter=zmomside(D_DECL(icell,jcell,kcell),ibucket_map)
              else
               print *,"veldir invalid"
               stop
              endif
              if (massquarter.ge.zero) then
               ! do nothing
              else
               print *,"massquarter cannot be negative, ivec=",ivec
               stop
              endif

              massface_total(ivec)=massface_total(ivec)+massquarter
              momface_total(ivec)=momface_total(ivec)+momquarter

             enddo !ivec=1,num_MAC_vectors

            enddo ! iside: do iside=-1,1,2

             ! mac velocity
            if (massface_total(1).gt.zero) then

             momface_total(1)=momface_total(1)/massface_total(1)
             ! LS>0 if clamped
             call SUB_clamped_LS(xclamped,cur_time,LS_clamped, &
                vel_clamped,temperature_clamped)
             if (LS_clamped.ge.zero) then
              momface_total(1)=vel_clamped(veldir)
             else if (LS_clamped.lt.zero) then
              ! do nothing
             else
              print *,"LS_clamped is NaN"
              stop
             endif

            else
             print *,"massface_total(1) invalid"
             stop
            endif
          
            if (num_MAC_vectors.eq.1) then
             ! do nothing
            else if (num_MAC_vectors.eq.2) then 
             if (massface_total(num_MAC_vectors).gt.zero) then
              momface_total(num_MAC_vectors)= &
                momface_total(num_MAC_vectors)/massface_total(num_MAC_vectors)
              if (veldir.eq.normdir+1) then
               momface_total(num_MAC_vectors)= &
                momface_total(num_MAC_vectors)+unode(D_DECL(i,j,k))
              else if ((veldir.ge.1).and.(veldir.le.SDIM)) then
               ! do nothing
              else
               print *,"veldir invalid"
               stop
              endif
             else if (massface_total(num_MAC_vectors).eq.zero) then
              momface_total(num_MAC_vectors)=zero
             else
              print *,"massface_total(num_MAC_vectors) invalid"
              stop
             endif
            else 
             print *,"num_MAC_vectors invalid"
             stop
            endif
             
            if (veldir.eq.1) then
             xvmac(D_DECL(i,j,k))=momface_total(1)

             if (num_MAC_vectors.eq.1) then
              ! do nothing
             else if (num_MAC_vectors.eq.2) then 
              xdmac(D_DECL(i,j,k))=momface_total(num_MAC_vectors)
             else
              print *,"num_MAC_vectors invalid"
              stop
             endif

            else if (veldir.eq.2) then
             yvmac(D_DECL(i,j,k))=momface_total(1)

             if (num_MAC_vectors.eq.1) then
              ! do nothing
             else if (num_MAC_vectors.eq.2) then 
              ydmac(D_DECL(i,j,k))=momface_total(num_MAC_vectors)
             else
              print *,"num_MAC_vectors invalid"
              stop
             endif

            else if ((veldir.eq.3).and.(SDIM.eq.3)) then
             zvmac(D_DECL(i,j,k))=momface_total(1)

             if (num_MAC_vectors.eq.1) then
              ! do nothing
             else if (num_MAC_vectors.eq.2) then 
              zdmac(D_DECL(i,j,k))=momface_total(num_MAC_vectors)
             else
              print *,"num_MAC_vectors invalid"
              stop
             endif

            else
             print *,"veldir invalid"
             stop
            endif

          else if ((maskleft.eq.zero).or.(maskright.eq.zero)) then
            ! do nothing
          else 
           print *,"maskleft or maskright invalid"
           stop
          endif

         else
          print *,"zapvel invalid build newmac"
          stop
         endif

        enddo
        enddo
        enddo  ! i,j,k

      enddo ! veldir=1..sdim

      return
      end subroutine FORT_BUILD_NEWMAC

      ! called from split_scalar_advection after 
      !  BUILD_SEMIREFINEVOF(tessellate==0)
      subroutine FORT_BUILD_MACVOF( &
       level, &
       finest_level, &
       normdir, &
       nrefine_vof, &
       nrefine_cen, &
       vofF,DIMS(vofF), &
       cenF,DIMS(cenF), &
       x_mac_old, &
       DIMS(x_mac_old), &
       xd_mac_old, & 
       DIMS(xd_mac_old), &
       xvof,DIMS(xvof), &
       xvel,DIMS(xvel), &  ! 1..num_MAC_vectors
       xvelslp,DIMS(xvelslp), &  ! xvelslope,xcen
       xlo,dx, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       nmat, &
       ngrow, &
       num_MAC_vectors, & ! num_MAC_vectors=1 or 2
       ngrowmac, &
       veldir)
      use probcommon_module
      use godunov_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: num_MAC_vectors !num_MAC_vectors=1 or 2.
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: normdir
      INTEGER_T, intent(in) :: nrefine_vof
      INTEGER_T, intent(in) :: nrefine_cen
      INTEGER_T, intent(in) :: nmat,ngrow
      INTEGER_T, intent(in) :: ngrowmac,veldir
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: DIMDEC(vofF) 
      INTEGER_T, intent(in) :: DIMDEC(cenF) 
      INTEGER_T, intent(in) :: DIMDEC(x_mac_old) 
      INTEGER_T, intent(in) :: DIMDEC(xd_mac_old) 
      INTEGER_T, intent(in) :: DIMDEC(xvof) 
      INTEGER_T, intent(in) :: DIMDEC(xvel) !1..num_MAC_vectors
      INTEGER_T, intent(in) :: DIMDEC(xvelslp) ! xvelslope,xcen
      REAL_T, intent(in) :: vofF(DIMV(vofF),nrefine_vof)
      REAL_T, intent(in) :: cenF(DIMV(cenF),nrefine_cen)
      REAL_T, intent(in) :: x_mac_old(DIMV(x_mac_old))
      REAL_T, intent(in) :: xd_mac_old(DIMV(xd_mac_old))
      REAL_T, intent(out) :: xvof(DIMV(xvof),nmat)
      REAL_T, intent(out) :: xvel(DIMV(xvel),num_MAC_vectors) 
      REAL_T, intent(out) :: xvelslp(DIMV(xvelslp),1+nmat) ! xvelslope,xcen
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)

      INTEGER_T icell,jcell,kcell,i,j,k,ii,jj,kk
      INTEGER_T igridlo(3),igridhi(3)
      INTEGER_T im
      REAL_T velmac(num_MAC_vectors)
      REAL_T volmatCV(nmat)
      REAL_T cenmatCV(nmat)
      REAL_T volCV_fluid
      REAL_T volCV_solid
      REAL_T volquarter
      INTEGER_T irefine,irefinecen,iside
      REAL_T xsten(-1:1,SDIM)
      INTEGER_T nhalf
      INTEGER_T ivec

      nhalf=1

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid build macvof"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid46"
       stop
      endif

      if ((num_materials_viscoelastic.ge.1).and. &
          (num_materials_viscoelastic.le.nmat)) then
       ! do nothing
      else
       print *,"num_materials_viscoelastic invalid"
       stop
      endif

      if ((num_MAC_vectors.eq.1).or. &
          (num_MAC_vectors.eq.2)) then
       ! do nothing
      else
       print *,"num_MAC_vectors invalid"
       stop
      endif

      if ((normdir.ge.0).and.(normdir.lt.SDIM)) then
       ! do nothing
      else
       print *,"normdir invalid BUILD_MACVOF"
       stop
      endif
      if (ngrow.lt.0) then
       print *,"ngrow invalid"
       stop
      endif
      if ((veldir.lt.1).or.(veldir.gt.SDIM)) then
       print *,"veldir invalid"
       stop
      endif
      if ((ngrowmac.lt.0).or.(ngrowmac.ge.ngrow)) then
       print *,"ngrowmac invalid"
       stop
      endif
      if (nrefine_vof.ne.2*nmat*SDIM) then
       print *,"nrefine_vof invalid"
       stop
      endif
      if (nrefine_cen.ne.2*nmat*SDIM*SDIM) then
       print *,"nrefine_cen invalid in build_macvof"
       stop
      endif
      call checkbound(fablo,fabhi,DIMS(x_mac_old),ngrowmac,veldir-1,1271)
      call checkbound(fablo,fabhi,DIMS(xd_mac_old),ngrowmac,veldir-1,1271)
      call checkbound(fablo,fabhi,DIMS(xvof),ngrowmac,veldir-1,1271)
      call checkbound(fablo,fabhi,DIMS(xvel),ngrowmac,veldir-1,1271)
      call checkbound(fablo,fabhi,DIMS(xvelslp),ngrowmac,veldir-1,1271)
      call checkbound(fablo,fabhi,DIMS(vofF),ngrow,-1,1272)
      call checkbound(fablo,fabhi,DIMS(cenF),ngrow,-1,1272)

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
        igridlo,igridhi,ngrowmac,veldir-1,20)

      ii=0
      jj=0
      kk=0
      if (veldir.eq.1) then
       ii=1
      else if (veldir.eq.2) then
       jj=1
      else if ((veldir.eq.3).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"veldir invalid in build_macvof"
       stop
      endif
       
      do i=igridlo(1),igridhi(1)
      do j=igridlo(2),igridhi(2)
      do k=igridlo(3),igridhi(3)

        ! veldir=1..sdim
       call gridstenMAC_level(xsten,i,j,k,level,nhalf,veldir-1,26)

       velmac(1)=x_mac_old(D_DECL(i,j,k))
       if (num_MAC_vectors.eq.2) then
        velmac(num_MAC_vectors)=xd_mac_old(D_DECL(i,j,k))
       else if (num_MAC_vectors.eq.1) then
        ! do nothing
       else
        print *,"num_MAC_vectors invalid"
        stop
       endif

       if (veldir.eq.1) then
        if (levelrz.eq.0) then
         ! do nothing
        else if (levelrz.eq.1) then
         if (SDIM.ne.2) then
          print *,"dimension bust"
          stop
         endif
         if (xsten(0,1).le.VOFTOL*dx(1)) then
          do ivec=1,num_MAC_vectors
           velmac(ivec)=zero
          enddo
         endif
        else if (levelrz.eq.3) then
         if (xsten(0,1).le.VOFTOL*dx(1)) then
          do ivec=1,num_MAC_vectors
           velmac(ivec)=zero
          enddo
         endif
        else
         print *,"levelrz invalid build macvof"
         stop
        endif
       else if (veldir.eq.2) then
        if (levelrz.eq.0) then
         ! do nothing
        else if (levelrz.eq.1) then
         if (SDIM.ne.2) then
          print *,"dimension bust"
          stop
         endif
         if (xsten(0,1).le.VOFTOL*dx(1)) then
          do ivec=1,num_MAC_vectors
           velmac(ivec)=zero
          enddo
         endif
        else if (levelrz.eq.3) then
         if (xsten(0,1).le.VOFTOL*dx(1)) then
          do ivec=1,num_MAC_vectors
           velmac(ivec)=zero
          enddo
         endif
        else
         print *,"levelrz invalid build macvof 2"
         stop
        endif
       else if ((veldir.eq.3).and.(SDIM.eq.3)) then
        if (levelrz.eq.0) then
         ! do nothing
        else if (levelrz.eq.1) then
         print *,"dimension bust"
         stop
        else if (levelrz.eq.3) then
         if (xsten(0,1).le.VOFTOL*dx(1)) then
          do ivec=1,num_MAC_vectors
           velmac(ivec)=zero
          enddo
         endif
        else
         print *,"levelrz invalid build macvof 3"
         stop
        endif
       else
        print *,"veldir invalid"
        stop
       endif

       do im=1,nmat
        volmatCV(im)=zero
        cenmatCV(im)=zero
       enddo
       volCV_fluid=zero
       volCV_solid=zero

       do iside=-1,1,2
        if (iside.eq.-1) then
         icell=i-ii
         jcell=j-jj
         kcell=k-kk
        else if (iside.eq.1) then
         icell=i
         jcell=j
         kcell=k
        else
         print *,"iside invalid"
         stop
        endif

        do im=1,nmat
         if (iside.eq.-1) then
          irefine=(veldir-1)*2*nmat+nmat+im
         else if (iside.eq.1) then
          irefine=(veldir-1)*2*nmat+im
         else
          print *,"side invalid"
          stop
         endif

           ! tessellate==0  (see BUILD_SEMIREFINEVOF)
         volquarter=vofF(D_DECL(icell,jcell,kcell),irefine)
         volmatCV(im)=volmatCV(im)+volquarter

         ! centroid in absolute coordinate system
         if ((normdir.ge.0).and.(normdir.lt.SDIM)) then

          if (iside.eq.-1) then !right side of left cell.
           irefinecen=(veldir-1)*2*nmat*SDIM+ &
            nmat*SDIM+(im-1)*SDIM+normdir+1
          else if (iside.eq.1) then !left side of right cell.
           irefinecen=(veldir-1)*2*nmat*SDIM+ &
            (im-1)*SDIM+normdir+1
          else
           print *,"iside invalid"
           stop
          endif

          cenmatCV(im)=cenmatCV(im)+ &
           volquarter*cenF(D_DECL(icell,jcell,kcell),irefinecen)
         else
          print *,"normdir invalid BUILD_MACVOF (2)"
          stop
         endif

         
         if (is_rigid(nmat,im).eq.0) then 
          volCV_fluid=volCV_fluid+volquarter
         else if (is_rigid(nmat,im).eq.1) then
          volCV_solid=volCV_solid+volquarter
         else
          print *,"is_rigid(nmat,im) invalid"
          stop
         endif

        enddo ! im=1..nmat

       enddo ! iside=-1,1,2

       if (volCV_solid.lt.zero) then
        print *,"volCV_solid invalid"
        stop
       endif

       if (volCV_fluid.lt.zero) then
        print *,"volCV_fluid invalid"
        stop
       else if (volCV_fluid.ge.zero) then

        do ivec=1,num_MAC_vectors
         xvel(D_DECL(i,j,k),ivec)=velmac(ivec)
        enddo

        do im=1,nmat

         if (volmatCV(im).lt.zero) then
          print *,"volmatCV invalid"
          print *,"im= ",im
          print *,"volmatCV= ",volmatCV(im)
          stop
         else if (volmatCV(im).eq.zero) then
          cenmatCV(im)=zero  ! placeholder
         else if (volmatCV(im).gt.zero) then
          if (volCV_fluid.le.zero) then
           print *,"volCV_fluid bust"
           stop
          endif
          cenmatCV(im)=cenmatCV(im)/volmatCV(im)
          volmatCV(im)=volmatCV(im)/volCV_fluid
          if ((volmatCV(im).lt.zero).or. &
              (volmatCV(im).gt.one+VOFTOL_SLOPES)) then
           print *,"volmatCV(im) invalid"
           stop
          endif
         else
          print *,"volmatCV(im) bust"
          stop
         endif
      
         if ((normdir.ge.0).and.(normdir.lt.SDIM)) then
           xvelslp(D_DECL(i,j,k),1+im)=cenmatCV(im) ! xcen
         else
           print *,"normdir invalid BUILD_MACVOF (3)"
           stop
         endif

         xvof(D_DECL(i,j,k),im)=volmatCV(im)
        enddo  ! im=1..nmat
       else
        print *,"volCV_fluid bust"
        stop
       endif
      enddo         
      enddo         
      enddo ! i,j,k (face center "conserved" variables) 

      return
      end subroutine FORT_BUILD_MACVOF


      ! masknbr:
      ! (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
      ! (2) =1 interior  =0 otherwise
      ! (3) =1 interior+ngrow-1  =0 otherwise
      ! (4) =1 interior+ngrow    =0 otherwise
      subroutine FORT_BUILD_SLOPES( &
       masknbr,DIMS(masknbr), &
       recon,DIMS(recon), &
       slsrc,DIMS(slsrc), &
       sldst,DIMS(sldst), &
       nc_conserve, &
       nmat, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       level, &
       finest_level, &
       velbc, &
       xlo,dx, &
       normdir, &
       ngrow,  &
       advection_order, &
       density_advection_order, &
       slope_limiter_option) 
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T nc_conserve
      INTEGER_T nmat
      INTEGER_T level
      INTEGER_T finest_level
      INTEGER_T normdir,ngrow
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T bfact
      INTEGER_T velbc(SDIM,2)
      REAL_T dx(SDIM)
      REAL_T xlo(SDIM)
      INTEGER_T advection_order(nmat)
      INTEGER_T density_advection_order(nmat)
      INTEGER_T slope_limiter_option
      INTEGER_T DIMDEC(masknbr) 
      INTEGER_T DIMDEC(recon) 
      INTEGER_T DIMDEC(slsrc)
      INTEGER_T DIMDEC(sldst)
      REAL_T masknbr(DIMV(masknbr))
      REAL_T recon(DIMV(recon),nmat*ngeom_recon)
      REAL_T slsrc(DIMV(slsrc),nc_conserve)
      REAL_T sldst(DIMV(sldst),nc_conserve)

      INTEGER_T igridlo(3),igridhi(3)
      INTEGER_T i,j,k
      INTEGER_T n
      INTEGER_T nden
      INTEGER_T ii,jj,kk
      INTEGER_T local_order
      INTEGER_T local_masknbr
      INTEGER_T dircheck
      INTEGER_T icheck
      INTEGER_T icrit
      REAL_T fcen(nmat)
      REAL_T fmixed
      REAL_T fsolid
      REAL_T fsolid_mixed
      REAL_T splus,sminus,scen,sleft,sright
      REAL_T dencen,denplus,denminus
      REAL_T slope_temp
      INTEGER_T im
      INTEGER_T im_primary
      INTEGER_T im_local
      INTEGER_T istate_local
      INTEGER_T nhalf
      REAL_T xsub(-3:3)
      REAL_T xsten(-1:1,SDIM)
      REAL_T dxleft,dxright
      INTEGER_T vofcomp

      if (bfact.lt.1) then
       print *,"bfact invalid47"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(masknbr),ngrow,-1,1247)
      call checkbound(fablo,fabhi,DIMS(slsrc),ngrow,-1,1248)
      call checkbound(fablo,fabhi,DIMS(sldst),1,-1,1249)
      call checkbound(fablo,fabhi,DIMS(recon),ngrow,-1,1251)

      if (nc_conserve.ne.SDIM+nmat*num_state_material) then
       print *,"nc_conserve invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if ((normdir.lt.0).or.(normdir.ge.SDIM)) then
       print *,"normdir invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid build slopes"
       stop
      endif

      if (num_state_material.ne.num_state_base+num_species_var) then
       print *,"num_state_material invalid"
       stop
      endif

      ii=0
      jj=0
      kk=0
      if (normdir.eq.0) then
       ii=1
      else if (normdir.eq.1) then
       jj=1
      else if ((normdir.eq.SDIM-1).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"normdir invalid"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi, &
       igridlo,igridhi,1)

      do i=igridlo(1),igridhi(1)
      do j=igridlo(2),igridhi(2)
      do k=igridlo(3),igridhi(3)
       do n=1,nc_conserve
        sldst(D_DECL(i,j,k),n)=zero
       enddo
      enddo
      enddo
      enddo

      do i=igridlo(1),igridhi(1)
      do j=igridlo(2),igridhi(2)
      do k=igridlo(3),igridhi(3)

       nhalf=1
       call gridsten_level(xsten,i,j,k,level,nhalf)

       local_order=2 
       
       local_masknbr=NINT(masknbr(D_DECL(i,j,k)))

       if (levelrz.eq.0) then
        ! do nothing
       else if ((levelrz.eq.1).or.(levelrz.eq.3)) then
        if (xsten(0,1).le.VOFTOL*dx(1)) then
         local_order=1
        endif
       else
        print *,"levelrz invalid build slopes"
        stop
       endif
       do dircheck=1,SDIM
        if (dircheck.eq.1) then
         icheck=i
        else if (dircheck.eq.2) then
         icheck=j
        else if ((dircheck.eq.3).and.(SDIM.eq.3)) then
         icheck=k
        else
         print *,"dircheck invalid"
         stop
        endif
        if (icheck.le.fablo(dircheck)) then
         if (velbc(dircheck,1).eq.INT_DIR) then
          if (local_masknbr.eq.0) then ! coarse/fine bdry
           local_order=1
          else if (local_masknbr.eq.1) then ! fine/fine bdry
           ! do nothing
          else
           print *,"masknbr invalid"
           stop
          endif
         else
          local_order=1 ! 1st order slopes touching or outside the domain
         endif
        else if (icheck.ge.fabhi(dircheck)) then
         if (velbc(dircheck,2).eq.INT_DIR) then
          if (local_masknbr.eq.0) then ! coarse/fine bdry
           local_order=1
          else if (local_masknbr.eq.1) then ! fine/fine bdry
           ! do nothing
          else
           print *,"masknbr invalid"
           stop
          endif
         else
          local_order=1 ! 1st order touching or outside the domain
         endif
        else
         ! do nothing
        endif
       enddo ! dircheck=1..sdim

       if (local_order.eq.1) then
        ! do nothing
       else if (local_order.eq.2) then

        if (normdir.eq.0) then
         icrit=i
        else if (normdir.eq.1) then
         icrit=j
        else if ((normdir.eq.2).and.(SDIM.eq.3)) then
         icrit=k
        else
         print *,"normdir invalid"
         stop
        endif

        nhalf=3 
        call gridsten1D_level(xsub,icrit,level,normdir+1,nhalf)
        dxleft=xsub(0)-xsub(-2)
        dxright=xsub(2)-xsub(0)
        if ((dxleft.le.zero).or.(dxright.le.zero)) then
         print *,"(dxleft.le.zero).or.(dxright.le.zero)"
         stop
        endif

        fsolid=zero 
        fsolid_mixed=zero 
        im_primary=0
        do im=1,nmat
         vofcomp=(im-1)*ngeom_recon+1
         fcen(im)=recon(D_DECL(i,j,k),vofcomp)
         if (is_rigid(nmat,im).eq.1) then
          fsolid=fsolid+fcen(im)
          fsolid_mixed=fsolid_mixed+fcen(im)+ &
            recon(D_DECL(i+ii,j+jj,k+kk),vofcomp)+ &
            recon(D_DECL(i-ii,j-jj,k-kk),vofcomp)
         else if (is_rigid(nmat,im).eq.0) then 
          if (im_primary.eq.0) then
           im_primary=im
          else if ((im_primary.ge.1).and.(im_primary.le.nmat)) then
           if (fcen(im).gt.fcen(im_primary)) then
            im_primary=im
           endif
          else
           print *,"im_primary invalid"
           stop
          endif
         else
          print *,"is_rigid(nmat,im) invalid"
          stop
         endif
        enddo ! im=1..nmat
        if ((im_primary.lt.1).or.(im_primary.gt.nmat)) then
         print *,"im_primary invalid"
         stop
        endif

        fsolid_mixed=fsolid_mixed/three

        vofcomp=(im_primary-1)*ngeom_recon+1

        fmixed=(recon(D_DECL(i+ii,j+jj,k+kk),vofcomp)+ &
                recon(D_DECL(i-ii,j-jj,k-kk),vofcomp)+ &
                fcen(im_primary))/three

        if ((fsolid_mixed.lt.zero).or. &
            (fsolid_mixed.gt.one+VOFTOL_SLOPES).or. &
            (fmixed.lt.zero).or. &
            (fmixed.gt.one+VOFTOL_SLOPES)) then
         print *,"fsolid_mixed or fmixed invalid"
         stop
        endif

        if ((fsolid.ge.half).or.(fsolid_mixed.ge.half)) then
         local_order=1
        else if ((fsolid.le.half).and.(fsolid_mixed.le.half)) then
         if ((fcen(im_primary).lt.half).or. &
             (fmixed.lt.half)) then
          local_order=1
         endif
        else
         print *,"fsolid or fsolid_mixed invalid"
         stop
        endif

        if (local_order.eq.2) then 
 
         do n=1,nc_conserve

          if ((n.ge.1).and.(n.le.SDIM)) then ! velocity
           local_order=advection_order(im_primary)
           if (local_order.eq.2) then
            if ((slope_limiter_option.eq.2).or. &
                (slope_limiter_option.eq.3)) then
             if ((abs(fmixed-one).ge.VOFTOL_SLOPES).or. &
                 (abs(fsolid_mixed).ge.VOFTOL_SLOPES)) then
              local_order=1
             endif
            else if ((slope_limiter_option.eq.1).or. &
                     (slope_limiter_option.eq.0)) then
             ! do nothing
            else
             print *,"slope_limiter_option invalid"
             stop
            endif 
           else if (local_order.eq.1) then
            ! do nothing
           else 
            print *,"local_order invalid"
            stop
           endif

           if (local_order.eq.1) then
            ! do nothing, slopes already 0
           else if (local_order.eq.2) then
            splus=slsrc(D_DECL(i+ii,j+jj,k+kk),n)
            sminus=slsrc(D_DECL(i-ii,j-jj,k-kk),n)
            scen=slsrc(D_DECL(i,j,k),n)

            nden=SDIM+(im_primary-1)*num_state_material+1
            dencen=slsrc(D_DECL(i,j,k),nden)
          
            if ((abs(fmixed-one).ge.VOFTOL_SLOPES).or. &
                (abs(fsolid_mixed).ge.VOFTOL_SLOPES)) then
             denplus=dencen
             denminus=dencen
            else if ((abs(fmixed-one).le.VOFTOL_SLOPES).and. &
                     (abs(fsolid_mixed).le.VOFTOL_SLOPES)) then
             denplus=slsrc(D_DECL(i+ii,j+jj,k+kk),nden)
             denminus=slsrc(D_DECL(i-ii,j-jj,k-kk),nden)
            else
             print *,"fmixed or fsolid_mixed invalid"
             stop
            endif 
            if ((dencen.le.zero).or. &
                (denplus.le.zero).or. &
                (denminus.le.zero)) then
             print *,"dencen,denplus, or denminus invalid"
             stop
            endif
            
            sleft=(scen*dencen-sminus*denminus)/dxleft
            sright=(splus*denplus-scen*dencen)/dxright

            if ((slope_limiter_option.eq.1).or. &
                (slope_limiter_option.eq.2)) then
             call minmod(sleft,sright,slope_temp)
            else if ((slope_limiter_option.eq.0).or. &
                     (slope_limiter_option.eq.3)) then
             slope_temp=half*(sleft+sright)
            else
             print *,"slope_limiter_option invalid"
             stop
            endif
            sldst(D_DECL(i,j,k),n)=slope_temp/dencen
           else
            print *,"local_order invalid"
            stop
           endif
          else if ((n.gt.SDIM).and.(n.le.nc_conserve)) then
           im_local=(n-SDIM-1)/num_state_material+1
           istate_local=n-SDIM-(im_local-1)*num_state_material 
           if ((im_local.lt.1).or.(im_local.gt.nmat)) then
            print *,"im_local invalid"
            stop
           endif

           if ((istate_local.lt.1).or. &
               (istate_local.gt.num_state_material)) then
            print *,"istate_local invalid"
            stop
           endif

           local_order=density_advection_order(im_local)

           if (local_order.eq.2) then

            if (im_local.eq.im_primary) then

             if ((istate_local.eq.1).or. &
                 (slope_limiter_option.eq.2).or. &
                 (slope_limiter_option.eq.3)) then ! density
              if ((abs(fmixed-one).ge.VOFTOL_SLOPES).or. &
                  (abs(fsolid_mixed).ge.VOFTOL_SLOPES)) then
               local_order=1
              endif
             else if ((istate_local.gt.1).and. &
                      ((slope_limiter_option.eq.0).or. &
                       (slope_limiter_option.eq.1))) then
              ! do nothing
             else
              print *,"istate_local or slope_limiter_option invalid"
              stop
             endif

             if ((istate_local.eq.1).or. &
                 (istate_local.eq.2)) then
              ! use density_advection_order
             else if ((istate_local.ge.num_state_base+1).and. &
                      (istate_local.le.num_state_base+num_species_var)) then
              local_order=1 ! must have 0<=massfrac<=1
             else
              print *,"istate_local invalid"
              stop
             endif  

            else if (im_local.ne.im_primary) then
             local_order=1
            else
             print *,"im_local or im_primary invalid"
             stop
            endif
              
           else if (local_order.eq.1) then
            ! do nothing
           else
            print *,"local_order invalid"
            stop
           endif

           if (local_order.eq.1) then
            ! do nothing, slopes already 0
           else if (local_order.eq.2) then
            splus=slsrc(D_DECL(i+ii,j+jj,k+kk),n)
            sminus=slsrc(D_DECL(i-ii,j-jj,k-kk),n)
            scen=slsrc(D_DECL(i,j,k),n)

            if ((istate_local.eq.1).or.  & ! density
                (istate_local.eq.2)) then ! energy
             if ((scen.le.zero).or. &
                 (splus.le.zero).or. &
                 (sminus.le.zero)) then
              print *,"scen,splus, or sminus.le.zero"
              stop
             endif
            else if ((istate_local.ge.3).and. &
                     (istate_local.le.num_state_material)) then
             ! check nothing
            else
             print *,"istate_local invalid"
             stop
            endif

            sleft=(scen-sminus)/dxleft
            sright=(splus-scen)/dxright

            if ((slope_limiter_option.eq.1).or. &
                (slope_limiter_option.eq.2)) then
             call minmod(sleft,sright,slope_temp)
            else if ((slope_limiter_option.eq.0).or. &
                     (slope_limiter_option.eq.3)) then
             slope_temp=half*(sleft+sright)
            else
             print *,"slope_limiter_option invalid"
             stop
            endif
            sldst(D_DECL(i,j,k),n)=slope_temp
           else
            print *,"local_order invalid"
            stop
           endif

          else
           print *,"n invalid"
           stop
          endif

         enddo ! n=1..nc_conserve

        else if (local_order.eq.1) then
         ! do nothing
        else
         print *,"local order invalid"
         stop
        endif

       else
        print *,"local order invalid"
        stop
       endif

      enddo
      enddo
      enddo ! i,j,k

      return
      end subroutine FORT_BUILD_SLOPES



      subroutine FORT_BUILD_SLOPES_FACE( &
       masknbr,DIMS(masknbr), &
       vfrac,DIMS(vfrac), &
       slsrc,DIMS(slsrc), &
       sldst,DIMS(sldst), &
       nmat, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       level, &
       finest_level, &
       velbc, &
       xlo,dx, &
       normdir, &
       slopedir, &
       ngrow,  &
       advection_order, &
       slope_limiter_option) 
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE


      INTEGER_T nmat
      INTEGER_T level
      INTEGER_T finest_level
      INTEGER_T normdir
      INTEGER_T slopedir
      INTEGER_T ngrow
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T bfact
      INTEGER_T velbc(SDIM,2)
      REAL_T dx(SDIM)
      REAL_T xlo(SDIM)
      INTEGER_T advection_order(nmat)
      INTEGER_T slope_limiter_option
      INTEGER_T DIMDEC(masknbr) 
      INTEGER_T DIMDEC(vfrac) 
      INTEGER_T DIMDEC(slsrc)
      INTEGER_T DIMDEC(sldst)
      REAL_T masknbr(DIMV(masknbr))
      REAL_T vfrac(DIMV(vfrac),nmat)
      REAL_T slsrc(DIMV(slsrc)) !slsrc(1)=xvel  slsrc(2)=xdisp
      REAL_T sldst(DIMV(sldst)) !sldst(1)=xvelslope,slpdst(2..nmat+1)=xcen

      INTEGER_T igridlo(3),igridhi(3)
      INTEGER_T i,j,k
      INTEGER_T icell,jcell,kcell
      INTEGER_T ii,jj,kk
      INTEGER_T iii,jjj,kkk
      INTEGER_T local_order
      INTEGER_T local_masknbr
      INTEGER_T side
      INTEGER_T dircheck
      INTEGER_T icheck
      INTEGER_T icrit
      REAL_T fcen(nmat)
      REAL_T fmixed
      REAL_T fsolid
      REAL_T fsolid_mixed
      REAL_T splus,sminus,scen,sleft,sright
      REAL_T slope_temp
      INTEGER_T im
      INTEGER_T im_primary
      INTEGER_T nhalf
      REAL_T xsub(-3:3)
      REAL_T xsten(-1:1,SDIM)
      REAL_T dxleft,dxright

      if (bfact.lt.1) then
       print *,"bfact invalid48"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(masknbr),ngrow,-1,1247)
      call checkbound(fablo,fabhi,DIMS(slsrc),ngrow,slopedir,1248)
      call checkbound(fablo,fabhi,DIMS(sldst),ngrow-1,slopedir,1249)
      call checkbound(fablo,fabhi,DIMS(vfrac),ngrow,slopedir,1251)

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if ((normdir.lt.0).or.(normdir.ge.SDIM)) then
       print *,"normdir invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid build slopes_face"
       stop
      endif

      iii=0
      jjj=0
      kkk=0
       ! mask is at cells, need to interpolate to faces if slopedir>=0
      if (slopedir.eq.0) then
       iii=1
      else if (slopedir.eq.1) then
       jjj=1
      else if ((slopedir.eq.SDIM-1).and.(SDIM.eq.3)) then
       kkk=1
      else
       print *,"slopedir invalid"
       stop
      endif

      ii=0
      jj=0
      kk=0
      if (normdir.eq.0) then
       ii=1
      else if (normdir.eq.1) then
       jj=1
      else if ((normdir.eq.SDIM-1).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"normdir invalid"
       stop
      endif

      if ((slopedir.ge.0).and.(slopedir.lt.SDIM)) then
       call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
        igridlo,igridhi,ngrow-1,slopedir,21)
      else
       print *,"slopedir invalid"
       stop
      endif

      do i=igridlo(1),igridhi(1)
      do j=igridlo(2),igridhi(2)
      do k=igridlo(3),igridhi(3)
       sldst(D_DECL(i,j,k))=zero
      enddo
      enddo
      enddo

      do i=igridlo(1),igridhi(1)
      do j=igridlo(2),igridhi(2)
      do k=igridlo(3),igridhi(3)

       nhalf=1
       if ((slopedir.ge.0).and.(slopedir.le.SDIM-1)) then
        call gridstenMAC_level(xsten,i,j,k,level,nhalf,slopedir,27)
       else
        print *,"slopedir invalid"
        stop
       endif

       local_order=2 
       do side=-1,1,2
       
        if (side.eq.-1) then
         icell=i-iii
         jcell=j-jjj
         kcell=k-kkk
        else if (side.eq.1) then
         icell=i
         jcell=j
         kcell=k
        else
         print *,"side invalid"
         stop
        endif

        local_masknbr=NINT(masknbr(D_DECL(icell,jcell,kcell)))

        if (levelrz.eq.0) then
         ! do nothing
        else if ((levelrz.eq.1).or.(levelrz.eq.3)) then
         if (xsten(0,1).le.VOFTOL*dx(1)) then
          local_order=1
         endif
        else
         print *,"levelrz invalid build slopes_face"
         stop
        endif
        do dircheck=1,SDIM
         if (dircheck.eq.1) then
          icheck=icell
         else if (dircheck.eq.2) then
          icheck=jcell
         else if ((dircheck.eq.3).and.(SDIM.eq.3)) then
          icheck=kcell
         else
          print *,"dircheck invalid"
          stop
         endif
         if (icheck.le.fablo(dircheck)) then
          if (velbc(dircheck,1).eq.INT_DIR) then
           if (local_masknbr.eq.0) then ! coarse/fine bdry
            local_order=1
           else if (local_masknbr.eq.1) then
            ! do nothing
           else
            print *,"masknbr invalid"
            stop
           endif
          else
           local_order=1 ! 1st order touching or outside the domain
          endif
         else if (icheck.ge.fabhi(dircheck)) then
          if (velbc(dircheck,2).eq.INT_DIR) then
           if (local_masknbr.eq.0) then ! coarse/fine bdry
            local_order=1
           else if (local_masknbr.eq.1) then
            ! do nothing
           else
            print *,"masknbr invalid"
            stop
           endif
          else
           local_order=1 ! 1st order touching or outside the domain
          endif
         else
          ! do nothing
         endif
        enddo ! dircheck=1..sdim
       enddo ! side

       if (local_order.eq.1) then
        ! do nothing
       else if (local_order.eq.2) then

        if (normdir.eq.0) then
         icrit=i
        else if (normdir.eq.1) then
         icrit=j
        else if ((normdir.eq.2).and.(SDIM.eq.3)) then
         icrit=k
        else
         print *,"normdir invalid"
         stop
        endif

        nhalf=3 
        if (normdir.eq.slopedir) then
         call gridsten1DMAC_level(xsub,icrit,level,normdir+1,nhalf)
        else if ((slopedir.ge.0).and.(slopedir.lt.SDIM)) then
         call gridsten1D_level(xsub,icrit,level,normdir+1,nhalf)
        else
         print *,"slopedir invalid"
         stop
        endif
        dxleft=xsub(0)-xsub(-2)
        dxright=xsub(2)-xsub(0)
        if ((dxleft.le.zero).or.(dxright.le.zero)) then
         print *,"(dxleft.le.zero).or.(dxright.le.zero)"
         stop
        endif

        fsolid=zero 
        fsolid_mixed=zero 
        im_primary=0
        do im=1,nmat
         fcen(im)=vfrac(D_DECL(i,j,k),im)
         if (is_rigid(nmat,im).eq.1) then
          fsolid=fsolid+fcen(im)
          fsolid_mixed=fsolid_mixed+fcen(im)+ &
            vfrac(D_DECL(i+ii,j+jj,k+kk),im)+ &
            vfrac(D_DECL(i-ii,j-jj,k-kk),im)
         else if (is_rigid(nmat,im).eq.0) then 
          if (im_primary.eq.0) then
           im_primary=im
          else if ((im_primary.ge.1).and.(im_primary.le.nmat)) then
           if (fcen(im).gt.fcen(im_primary)) then
            im_primary=im
           endif
          else
           print *,"im_primary invalid"
           stop
          endif
         else
          print *,"is_rigid(nmat,im) invalid"
          stop
         endif
        enddo ! im=1..nmat
        if ((im_primary.lt.1).or.(im_primary.gt.nmat)) then
         print *,"im_primary invalid"
         stop
        endif

        fsolid_mixed=fsolid_mixed/three

        fmixed=(vfrac(D_DECL(i+ii,j+jj,k+kk),im_primary)+ &
                vfrac(D_DECL(i-ii,j-jj,k-kk),im_primary)+ &
                fcen(im_primary))/three

        if ((fsolid_mixed.lt.zero).or. &
            (fsolid_mixed.gt.one+VOFTOL_SLOPES).or. &
            (fmixed.lt.zero).or. &
            (fmixed.gt.one+VOFTOL_SLOPES)) then
         print *,"fsolid_mixed or fmixed invalid"
         stop
        endif

        if ((fsolid.ge.half).or.(fsolid_mixed.ge.half)) then
         local_order=1
        else if ((fsolid.le.half).and.(fsolid_mixed.le.half)) then
         if ((fcen(im_primary).lt.half).or. &
             (fmixed.lt.half)) then
          local_order=1
         endif
        else
         print *,"fsolid or fsolid_mixed invalid"
         stop
        endif

        if (local_order.eq.2) then 

         local_order=advection_order(im_primary)
         if (local_order.eq.2) then
          if ((slope_limiter_option.eq.2).or. &
              (slope_limiter_option.eq.3)) then
           if ((abs(fmixed-one).ge.VOFTOL_SLOPES).or. &
               (abs(fsolid_mixed).ge.VOFTOL_SLOPES)) then
            local_order=1
           endif
          else if ((slope_limiter_option.eq.1).or. &
                   (slope_limiter_option.eq.0)) then
           ! do nothing
          else
           print *,"slope_limiter_option invalid"
           stop
          endif 
         else if (local_order.eq.1) then
          ! do nothing
         else 
          print *,"local_order invalid"
          stop
         endif

         if (local_order.eq.1) then
          ! do nothing, slopes already 0
         else if (local_order.eq.2) then

          splus=slsrc(D_DECL(i+ii,j+jj,k+kk))
          sminus=slsrc(D_DECL(i-ii,j-jj,k-kk))
          scen=slsrc(D_DECL(i,j,k))

          sleft=(scen-sminus)/dxleft
          sright=(splus-scen)/dxright

          if ((slope_limiter_option.eq.1).or. &
              (slope_limiter_option.eq.2)) then
           call minmod(sleft,sright,slope_temp)
          else if ((slope_limiter_option.eq.0).or. &
                   (slope_limiter_option.eq.3)) then
           slope_temp=half*(sleft+sright)
          else
           print *,"slope_limiter_option invalid"
           stop
          endif
          sldst(D_DECL(i,j,k))=slope_temp
         else
          print *,"local_order invalid"
          stop
         endif
        else if (local_order.eq.1) then
         ! do nothing
        else
         print *,"local order invalid"
         stop
        endif
       else
        print *,"local order invalid"
        stop
       endif

      enddo
      enddo
      enddo ! i,j,k

      return
      end subroutine FORT_BUILD_SLOPES_FACE




         ! 1=T11 2=T12 3=T22 4=T33 5=T13 6=T23
         ! rhoinverse is 1/den
      subroutine FORT_TENSORFORCE( &
       massface_index, &
       vofface_index, &
       ncphys, &
       nstate, &
       xlo,dx,  &
       xface,DIMS(xface), &
       yface,DIMS(yface), &
       zface,DIMS(zface), &
       xflux,DIMS(xflux), &
       yflux,DIMS(yflux), &
       zflux,DIMS(zflux), &
       lsfab,DIMS(lsfab), &
       rhoinverse, &
       DIMS(rhoinverse), &
       velnew,DIMS(velnew), &
       tensor,DIMS(tensor), &
       tilelo,tilehi, &
       fablo,fabhi,bfact,level, &
       dt,irz, &
       im_parm, & ! 0..nmat-1
       viscoelastic_model, &
       nmat,nden)
      use probcommon_module
      use global_utility_module
      use MOF_routines_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: massface_index
      INTEGER_T, intent(in) :: vofface_index
      INTEGER_T, intent(in) :: ncphys
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: im_parm
      INTEGER_T, intent(in) :: viscoelastic_model
      INTEGER_T, intent(in) :: nden,nstate
      INTEGER_T, intent(in) :: level
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)

      INTEGER_T, intent(in) :: DIMDEC(xface)
      INTEGER_T, intent(in) :: DIMDEC(yface)
      INTEGER_T, intent(in) :: DIMDEC(zface)
      INTEGER_T, intent(in) :: DIMDEC(xflux)
      INTEGER_T, intent(in) :: DIMDEC(yflux)
      INTEGER_T, intent(in) :: DIMDEC(zflux)
      INTEGER_T, intent(in) :: DIMDEC(lsfab)
      INTEGER_T, intent(in) :: DIMDEC(rhoinverse)
      INTEGER_T, intent(in) :: DIMDEC(velnew)
      INTEGER_T, intent(in) :: DIMDEC(tensor)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: xface(DIMV(xface),ncphys)
      REAL_T, intent(in) :: yface(DIMV(yface),ncphys)
      REAL_T, intent(in) :: zface(DIMV(zface),ncphys)
      REAL_T, intent(in) :: xflux(DIMV(xflux),SDIM*SDIM)
      REAL_T, intent(in) :: yflux(DIMV(yflux),SDIM*SDIM)
      REAL_T, intent(in) :: zflux(DIMV(zflux),SDIM*SDIM)
      REAL_T, intent(in) :: lsfab(DIMV(lsfab),nmat*(1+SDIM))
      REAL_T, intent(in) :: rhoinverse(DIMV(rhoinverse),nmat+1)
      REAL_T, intent(inout) :: velnew(DIMV(velnew),SDIM)
      REAL_T, intent(in) :: tensor(DIMV(tensor),FORT_NUM_TENSOR_TYPE)
      REAL_T, intent(in) :: dt
      INTEGER_T, intent(in) :: irz
      INTEGER_T :: i,j,k
      INTEGER_T :: xflux_comp
      INTEGER_T :: veldir
      REAL_T    :: deninv
      INTEGER_T dir
      REAL_T rplus,rminus,rval,RRX,RRY
      INTEGER_T mask_left,mask_right
      INTEGER_T mask_center(SDIM)
      INTEGER_T local_mask
      INTEGER_T mask_array(D_DECL(-1:1,-1:1,-1:1))
      INTEGER_T i1,j1,k1
      REAL_T hx,hy,hz
      REAL_T force(SDIM)
      REAL_T bodyforce
      INTEGER_T klo_stencil,khi_stencil
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf
      INTEGER_T imlocal
      REAL_T LScen(nmat)
      REAL_T Q(D_DECL(-1:1,-1:1,-1:1),3,3)

      REAL_T xflux_local(-1:1,SDIM,SDIM)
      REAL_T yflux_local(-1:1,SDIM,SDIM)
      REAL_T zflux_local(-1:1,SDIM,SDIM)

      REAL_T n_elastic(SDIM)
      INTEGER_T ii,jj
      INTEGER_T iQ_minus,iQ_plus
      INTEGER_T jQ_minus,jQ_plus
      INTEGER_T kQ_minus,kQ_plus
      INTEGER_T dir_outer

      nhalf=3

      if (bfact.lt.1) then
       print *,"bfact invalid49"
       stop
      endif
      if (ncphys.ne.vofface_index+2*nmat) then
       print *,"ncphys invalid"
       stop
      endif
      if (vofface_index.ne.massface_index+2*nmat) then
       print *,"vofface_index or massface_index invalid"
       stop
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      if (irz.ne.levelrz) then
       print *,"irz invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (nstate.ne.(SDIM+1)+ &
          nmat*(num_state_material+ngeom_raw)+1) then
       print *,"nstate invalid"
       stop
      endif
      if ((im_parm.lt.0).or.(im_parm.ge.nmat)) then
       print *,"im_parm invalid23"
       stop
      endif

      if (nden.ne.nmat*num_state_material) then
       print *,"nden invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(xface),0,0,263)
      call checkbound(fablo,fabhi,DIMS(yface),0,1,263)
      call checkbound(fablo,fabhi,DIMS(zface),0,SDIM-1,263)

      call checkbound(fablo,fabhi,DIMS(xflux),0,0,263)
      call checkbound(fablo,fabhi,DIMS(yflux),0,1,263)
      call checkbound(fablo,fabhi,DIMS(zflux),0,SDIM-1,263)

      call checkbound(fablo,fabhi,DIMS(lsfab),1,-1,7)
      call checkbound(fablo,fabhi, &
       DIMS(rhoinverse), &
       1,-1,7)
      call checkbound(fablo,fabhi,DIMS(velnew),1,-1,7)
      call checkbound(fablo,fabhi,DIMS(tensor),1,-1,7)

      if (SDIM.eq.3) then
       klo_stencil=-1
       khi_stencil=1
      else if (SDIM.eq.2) then
       klo_stencil=0
       khi_stencil=0
      else
       print *,"dimension bust"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
       call gridsten_level(xsten,i,j,k,level,nhalf)

       if (irz.eq.0) then
        ! do nothing
       else if (irz.eq.1) then
        if (SDIM.ne.2) then
         print *,"dimension bust"
         stop
        endif
        if (xsten(0,1).le.zero) then
         print *,"no neg domain in r-z"
         stop
        endif
       else if (irz.eq.3) then
        if (xsten(-2,1).le.zero) then
         print *,"no neg domain in r-T"
         stop
        endif
       else
        print *,"irz invalid"
        stop
       endif

       do i1=-1,1
       do j1=-1,1
       do k1=klo_stencil,khi_stencil
        do imlocal=1,nmat
         LScen(imlocal)=lsfab(D_DECL(i+i1,j+j1,k+k1),imlocal)
        enddo
        call get_primary_material(LScen,nmat,local_mask)
        if ((local_mask.eq.im_parm+1).and. &
            (LScen(im_parm+1).gt.zero)) then
         local_mask=1
        else if ((local_mask.ge.1).and.(local_mask.le.nmat)) then
         local_mask=0
        else
         print *,"local_mask invalid"
         stop
        endif
        mask_array(D_DECL(i1,j1,k1))=local_mask

        do ii=1,3
        do jj=1,3
         Q(D_DECL(i1,j1,k1),ii,jj)=zero
        enddo
        enddo
        Q(D_DECL(i1,j1,k1),1,1)=tensor(D_DECL(i+i1,j+j1,k+k1),1)
        Q(D_DECL(i1,j1,k1),1,2)=tensor(D_DECL(i+i1,j+j1,k+k1),2)
        Q(D_DECL(i1,j1,k1),2,2)=tensor(D_DECL(i+i1,j+j1,k+k1),3)
        Q(D_DECL(i1,j1,k1),3,3)=tensor(D_DECL(i+i1,j+j1,k+k1),4)
#if (AMREX_SPACEDIM==3)
        Q(D_DECL(i1,j1,k1),1,3)=tensor(D_DECL(i+i1,j+j1,k+k1),5)
        Q(D_DECL(i1,j1,k1),2,3)=tensor(D_DECL(i+i1,j+j1,k+k1),6)
#endif
        Q(D_DECL(i1,j1,k1),2,1)=Q(D_DECL(i1,j1,k1),1,2)
        Q(D_DECL(i1,j1,k1),3,1)=Q(D_DECL(i1,j1,k1),1,3)
        Q(D_DECL(i1,j1,k1),3,2)=Q(D_DECL(i1,j1,k1),2,3)

       enddo !k1=klo_stencil..khi_stencil
       enddo !j1=-1..1
       enddo !i1=-1..1

       if (irz.eq.0) then
        rval=one
        rplus=one
        rminus=one
        RRX=one
        RRY=one
       else if ((irz.eq.1).or.(irz.eq.3)) then
        RRX=one
        RRY=one
        rval=xsten(0,1)
        rplus=xsten(1,1)
        rminus=xsten(-1,1)
        if ((rval.lt.zero).or. &
            (rplus.lt.-VOFTOL*dx(1)).or. &
            (rminus.lt.-VOFTOL*dx(1))) then
         print *,"r values cannot be negative"
         stop
        endif
        RRX=rval
        if (irz.eq.3) then
         RRY=rval
        endif 
       else
        print *,"irz invalid"
        stop
       endif

       local_mask=mask_array(D_DECL(0,0,0))

        ! We do not have to check to see if the current cell is
        ! clamped since:
        !  (a) interp to MAC will override the MAC velocity if
        !      either one of the adjoining cells is clamped.
        !  (b) The new cell velocity after project is the interpolation
        !      of the adjoining MAC velocities if one of the adjoining
        !      MAC faces is a clamped face.  i.e., if an adjoining
        !      cell of a cell is clamped, then that cell gets the
        !      average of the adjoining MAC velocities.
        ! im_parm dominates the center cell.
       if (local_mask.eq.1) then
 
         ! im_parm=0..nmat-1
        do dir=1,SDIM
         n_elastic(dir)=lsfab(D_DECL(i,j,k),nmat+im_parm*SDIM+dir)
        enddo
        call normalize_vector(n_elastic)

        hx=(xsten(1,1)-xsten(-1,1))*RRX
        hy=(xsten(1,2)-xsten(-1,2))*RRY
        hz=xsten(1,SDIM)-xsten(-1,SDIM)

        if ((hx.gt.zero).and. &
            (hy.gt.zero).and. &
            (hz.gt.zero)) then
         ! do nothing
        else
         print *,"hx,hy, or hz invalid"
         stop
        endif

        if ((viscoelastic_model.eq.0).or. &
            (viscoelastic_model.eq.1).or. &
            (viscoelastic_model.eq.3)) then !incremental

         do dir_outer=1,SDIM

          iQ_minus=0
          jQ_minus=0
          kQ_minus=0

          iQ_plus=0
          jQ_plus=0
          kQ_plus=0

          if (dir_outer.eq.1) then
           iQ_minus=-1
           iQ_plus=1
          else if (dir_outer.eq.2) then
           jQ_minus=-1
           jQ_plus=1
          else if ((dir_outer.eq.3).and.(SDIM.eq.3)) then
           kQ_minus=-1
           kQ_plus=1
          else
           print *,"dir_outer invalid"
           stop
          endif

          do veldir=1,SDIM

           do dir=1,SDIM

            if (dir_outer.eq.1) then
             xflux_local(1,veldir,dir)= &
               (Q(D_DECL(iQ_plus,jQ_plus,kQ_plus),veldir,dir)+ &
                Q(D_DECL(0,0,0),veldir,dir))/two
             xflux_local(-1,veldir,dir)= &
               (Q(D_DECL(iQ_minus,jQ_minus,kQ_minus),veldir,dir)+ &
                Q(D_DECL(0,0,0),veldir,dir))/two
             xflux_local(0,veldir,dir)=Q(D_DECL(0,0,0),veldir,dir)
            else if (dir_outer.eq.2) then
             yflux_local(1,veldir,dir)= &
               (Q(D_DECL(iQ_plus,jQ_plus,kQ_plus),veldir,dir)+ &
                Q(D_DECL(0,0,0),veldir,dir))/two
             yflux_local(-1,veldir,dir)= &
               (Q(D_DECL(iQ_minus,jQ_minus,kQ_minus),veldir,dir)+ &
                Q(D_DECL(0,0,0),veldir,dir))/two
             yflux_local(0,veldir,dir)=Q(D_DECL(0,0,0),veldir,dir)
            else if ((dir_outer.eq.3).and.(SDIM.eq.3)) then
             zflux_local(1,veldir,dir)= &
               (Q(D_DECL(iQ_plus,jQ_plus,kQ_plus),veldir,dir)+ &
                Q(D_DECL(0,0,0),veldir,dir))/two
             zflux_local(-1,veldir,dir)= &
               (Q(D_DECL(iQ_minus,jQ_minus,kQ_minus),veldir,dir)+ &
                Q(D_DECL(0,0,0),veldir,dir))/two
             zflux_local(0,veldir,dir)=Q(D_DECL(0,0,0),veldir,dir)
            else
             print *,"dimension bust"
             stop
            endif
           enddo ! dir=1..sdim
          enddo ! veldir=1..sdim
         enddo ! dir_outer=1..sdim

        else if (viscoelastic_model.eq.2) then ! displacement gradient

          ! xflux,yflux,zflux come from "FORT_CROSSTERM_ELASTIC" in which
          ! "stress_from_strain" is applied to the displacement gradient
          ! matrix, then the resulting stress is multiplied by
          ! "elastic_viscosity * visc_coef" (see DERVISC)
         xflux_comp=1
         do veldir=1,SDIM
         do dir=1,SDIM
          xflux_local(1,veldir,dir)=xflux(D_DECL(i+1,j,k),xflux_comp)
          xflux_local(-1,veldir,dir)=xflux(D_DECL(i,j,k),xflux_comp)
          xflux_local(0,veldir,dir)=half*(xflux_local(1,veldir,dir)+ &
                  xflux_local(-1,veldir,dir))
          yflux_local(1,veldir,dir)=yflux(D_DECL(i,j+1,k),xflux_comp)
          yflux_local(-1,veldir,dir)=yflux(D_DECL(i,j,k),xflux_comp)
          yflux_local(0,veldir,dir)=half*(yflux_local(1,veldir,dir)+ &
                  yflux_local(-1,veldir,dir))

          if (SDIM.eq.3) then
           zflux_local(1,veldir,dir)=zflux(D_DECL(i,j,k+1),xflux_comp)
           zflux_local(-1,veldir,dir)=zflux(D_DECL(i,j,k),xflux_comp)
           zflux_local(0,veldir,dir)=half*(zflux_local(1,veldir,dir)+ &
                  zflux_local(-1,veldir,dir))
          else if (SDIM.eq.2) then
           ! do nothing
          else
           print *,"dimension bust"
           stop
          endif
          xflux_comp=xflux_comp+1
         enddo ! dir=1..sdim
         enddo ! veldir=1..sdim
         if (xflux_comp-1.eq.SDIM*SDIM) then
          ! do nothing
         else
          print *,"xflux_comp invalid"
          stop
         endif

        else
         print *,"viscoelastic_model invalid"
         stop
        endif

        dir=1
        mask_left=mask_array(D_DECL(-1,0,0))
        mask_right=mask_array(D_DECL(1,0,0))
         ! declared in: GLOBALUTIL.F90
         ! mask_center=1 if mask_left==1 or mask_right==1
         ! mask_center=0 if mask_left==0 and mask_right==0
        call project_tensor(mask_center(dir),n_elastic, &
         mask_left,mask_right,xflux_local,dir)

        dir=2
        mask_left=mask_array(D_DECL(0,-1,0))
        mask_right=mask_array(D_DECL(0,1,0))
        call project_tensor(mask_center(dir),n_elastic, &
         mask_left,mask_right,yflux_local,dir)

        if (SDIM.eq.3) then
         dir=SDIM
         mask_left=mask_array(D_DECL(0,0,-1))
         mask_right=mask_array(D_DECL(0,0,1))
         call project_tensor(mask_center(dir),n_elastic, &
           mask_left,mask_right,zflux_local,dir)
        endif

        do veldir=1,SDIM

         force(veldir)=zero

         dir=1
         force(veldir)=force(veldir)+ &
          mask_center(dir)* &
            (rplus*xflux_local(1,veldir,dir)- &
             rminus*xflux_local(-1,veldir,dir))/hx

         dir=2
         force(veldir)=force(veldir)+ &
          mask_center(dir)* &
            (yflux_local(1,veldir,dir)- &
             yflux_local(-1,veldir,dir))/hy

         if (SDIM.eq.3) then
          dir=SDIM
          force(veldir)=force(veldir)+ &
           mask_center(dir)* &
             (zflux_local(1,veldir,dir)- &
              zflux_local(-1,veldir,dir))/hz
         endif

        enddo ! veldir=1..sdim
                 
        if (irz.eq.0) then
         ! do nothing
        else if (irz.eq.1) then
         if (SDIM.ne.2) then
          print *,"dimension bust"
          stop
         endif
          ! in RZ,
          !  Q(3,3) gets 2 * hoop_22 (2 * xdisp/r) in 
          !  subroutine local_tensor_from_xdisplace
          !  -T33/r
          ! in XY,
          !  Q(3,3)=0.0
          ! in R-Theta,
          !  Q(3,3)=0.0
          ! local_tensor_from_xdisplace called from FORT_UPDATETENSOR
          ! In FORT_MAKETENSOR, Q is multiplied by the elastic bulk modulus.
         veldir=1
         bodyforce=-Q(D_DECL(0,0,0),3,3)/rval
         if (abs(bodyforce).lt.OVERFLOW_CUTOFF) then
          ! do nothing
         else
          print *,"bodyforce overflow bodyforce,rval:",bodyforce,rval
          stop
         endif
         force(veldir)=force(veldir)+bodyforce
        else if (irz.eq.3) then
          ! -T22/r
         veldir=1
         bodyforce=-Q(D_DECL(0,0,0),2,2)/rval
         force(veldir)=force(veldir)+bodyforce
        else
         print *,"irz invalid"
         stop
        endif 

        if (irz.eq.0) then
         ! do nothing
        else if (irz.eq.1) then
         if (SDIM.ne.2) then
          print *,"dimension bust"
          stop
         endif
         ! do nothing
        else if (irz.eq.3) then ! T12/r
         veldir=2
         bodyforce=Q(D_DECL(0,0,0),1,2)/rval
         force(veldir)=force(veldir)+bodyforce
        else
         print *,"irz invalid"
         stop
        endif 

        if (is_prescribed(nmat,im_parm+1).eq.1) then
         print *,"im_parm should not be an is_prescribed material"
         stop
        else if (is_prescribed(nmat,im_parm+1).eq.0) then
         do veldir=1,SDIM
          force(veldir)=force(veldir)*dt
         enddo
        else
         print *,"is_prescribed invalid"
         stop
        endif

        do veldir=1,SDIM
         deninv=rhoinverse(D_DECL(i,j,k),1)

         if (deninv.ge.zero) then 
          if (abs(force(veldir)).lt.OVERFLOW_CUTOFF) then
           ! do nothing
          else
           print *,"viscoelastic overflow veldir,force ",veldir,force(veldir)
           print *,"i,j,k,deninv ",i,j,k,deninv
           stop
          endif

          velnew(D_DECL(i,j,k),veldir)= &
           velnew(D_DECL(i,j,k),veldir)+force(veldir)*deninv
         else
          print *,"deninv invalid"
          stop
         endif
        enddo ! veldir

        ! im_parm does not dominate the center cell.
       else if (local_mask.eq.0) then
        ! do nothing (velnew is not incremented)
       else
        print *,"local_mask invalid"
        stop
       endif

      enddo
      enddo
      enddo ! i,j,k
 
      return
      end subroutine FORT_TENSORFORCE



         ! 1=T11 2=T12 3=T22 4=T33 5=T13 6=T23
         ! rhoinverse is 1/den
      subroutine FORT_TENSORHEAT( &
       massface_index, &
       vofface_index, &
       ncphys, &
       ntensor, &
       nstate, &
       xlo,dx,  &
       xface,DIMS(xface), &
       yface,DIMS(yface), &
       zface,DIMS(zface), &
       lsfab,DIMS(lsfab), &
       DeDTinverse, &
       DIMS(DeDTinverse), &
       vischeat,DIMS(vischeat), &
       tensor,DIMS(tensor), &
       gradu,DIMS(gradu), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       level, &
       dt, &
       irz, &
       im_parm, &
       nmat,nden)
      use probcommon_module
      use global_utility_module
      use MOF_routines_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: massface_index
      INTEGER_T, intent(in) :: vofface_index
      INTEGER_T, intent(in) :: ncphys
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: im_parm
      INTEGER_T, intent(in) :: nden,nstate,level
      INTEGER_T, intent(in) :: ntensor
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T, intent(in) :: DIMDEC(xface)
      INTEGER_T, intent(in) :: DIMDEC(yface)
      INTEGER_T, intent(in) :: DIMDEC(zface)
      INTEGER_T, intent(in) :: DIMDEC(lsfab)
      INTEGER_T, intent(in) :: DIMDEC(DeDTinverse)
      INTEGER_T, intent(in) :: DIMDEC(vischeat)
      INTEGER_T, intent(in) :: DIMDEC(tensor)
      INTEGER_T, intent(in) :: DIMDEC(gradu)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: xface(DIMV(xface),ncphys)
      REAL_T, intent(in) :: yface(DIMV(yface),ncphys)
      REAL_T, intent(in) :: zface(DIMV(zface),ncphys)
      REAL_T, intent(in) :: lsfab(DIMV(lsfab),nmat*(1+SDIM))
      REAL_T, intent(in) :: DeDTinverse(DIMV(DeDTinverse),nmat+1)
      REAL_T, intent(inout) :: vischeat(DIMV(vischeat))
      REAL_T, intent(in) :: tensor(DIMV(tensor),FORT_NUM_TENSOR_TYPE)
      REAL_T, intent(in) :: gradu(DIMV(gradu),ntensor)
      REAL_T, intent(in) :: dt
      INTEGER_T, intent(in) :: irz
      INTEGER_T :: i,j,k
      INTEGER_T veldir,dir
      INTEGER_T local_mask
      REAL_T IEforce,Tforce
      REAL_T one_over_DeDT
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf
      REAL_T local_gradu(3,3)
      INTEGER_T ux,vx,wx,uy,vy,wy,uz,vz,wz
      INTEGER_T nbase
      INTEGER_T grdcomp
      INTEGER_T imlocal
      REAL_T LScen(nmat)
      REAL_T Q(3,3)
      INTEGER_T ii,jj

      nhalf=3

      if (bfact.lt.1) then
       print *,"bfact invalid50"
       stop
      endif

      if (ncphys.ne.vofface_index+2*nmat) then
       print *,"ncphys invalid"
       stop
      endif
      if (vofface_index.ne.massface_index+2*nmat) then
       print *,"vofface_index or massface_index invalid"
       stop
      endif
     
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (irz.ne.levelrz) then
       print *,"irz invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (nstate.ne.(SDIM+1)+ &
          nmat*(num_state_material+ngeom_raw)+1) then
       print *,"nstate invalid"
       stop
      endif
      if ((im_parm.lt.0).or.(im_parm.ge.nmat)) then
       print *,"im_parm invalid24"
       stop
      endif

      if (nden.ne.nmat*num_state_material) then
       print *,"nden invalid"
       stop
      endif
      if (ntensor.ne.SDIM*SDIM) then
       print *,"ntensor invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(xface),0,0,263)
      call checkbound(fablo,fabhi,DIMS(yface),0,1,263)
      call checkbound(fablo,fabhi,DIMS(zface),0,SDIM-1,263)

      call checkbound(fablo,fabhi,DIMS(lsfab),1,-1,7)
      call checkbound(fablo,fabhi, &
       DIMS(DeDTinverse), &
       0,-1,7)
      call checkbound(fablo,fabhi,DIMS(vischeat),0,-1,7)
      call checkbound(fablo,fabhi,DIMS(tensor),0,-1,7)
      call checkbound(fablo,fabhi,DIMS(gradu),0,-1,7)

! u_x,v_x,w_x, u_y,v_y,w_y, u_z,v_z,w_z;  
      call tensorcomp_matrix(ux,uy,uz,vx,vy,vz,wx,wy,wz)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
       call gridsten_level(xsten,i,j,k,level,nhalf)

       do veldir=1,3
       do dir=1,3
        local_gradu(veldir,dir)=zero
       enddo 
       enddo 
       do dir=1,SDIM
        if (dir.eq.1) then
         nbase=ux-1
        else if (dir.eq.2) then
         nbase=uy-1
        else if ((dir.eq.SDIM).and.(SDIM.eq.3)) then
         nbase=uz-1
        else
         print *,"dir invalid viscoelastic heating"
         stop
        endif
        do veldir=1,SDIM
         grdcomp=nbase+veldir
         local_gradu(veldir,dir)=gradu(D_DECL(i,j,k),grdcomp)
        enddo ! veldir
       enddo ! dir

       if (irz.eq.0) then
        ! do nothing
       else if (irz.eq.1) then
        if (SDIM.ne.2) then
         print *,"dimension bust"
         stop
        endif
        if (xsten(0,1).le.zero) then
         print *,"no neg domain in r-z"
         stop
        endif
       else if (irz.eq.3) then
        if (xsten(-2,1).le.zero) then
         print *,"no neg domain in r-T"
         stop
        endif
       else
        print *,"irz invalid"
        stop
       endif

       do imlocal=1,nmat
        LScen(imlocal)=lsfab(D_DECL(i,j,k),imlocal)
       enddo
       call get_primary_material(LScen,nmat,local_mask)

       if ((local_mask.eq.im_parm+1).and. &
           (LScen(im_parm+1).gt.zero)) then
        local_mask=1
       else if ((local_mask.ge.1).and.(local_mask.le.nmat)) then
        local_mask=0
       else
        print *,"local_mask invalid"
        stop
       endif

       do ii=1,3
       do jj=1,3
        Q(ii,jj)=zero
       enddo
       enddo
       Q(1,1)=tensor(D_DECL(i,j,k),1)
       Q(1,2)=tensor(D_DECL(i,j,k),2)
       Q(2,2)=tensor(D_DECL(i,j,k),3)
       Q(3,3)=tensor(D_DECL(i,j,k),4)
#if (AMREX_SPACEDIM==3)
       Q(1,3)=tensor(D_DECL(i,j,k),5)
       Q(2,3)=tensor(D_DECL(i,j,k),6)
#endif
       Q(2,1)=Q(1,2)
       Q(3,1)=Q(1,3)
       Q(3,2)=Q(2,3)

        ! E: div( u dot tau )=(u tau_11)_x+(u tau_12)_y+(u tau_13)_z+...
        ! e: grad u : tau

       if (local_mask.eq.0) then
        IEforce=zero
       else if (local_mask.eq.1) then

        IEforce=zero
        do veldir=1,SDIM
        do dir=1,SDIM
         IEforce=IEforce+local_gradu(veldir,dir)*Q(veldir,dir)
        enddo
        enddo

        one_over_DeDT=DeDTinverse(D_DECL(i,j,k),1)  ! 1/(rho cv)

        if (one_over_DeDT.le.zero) then
         print *,"one_over_DeDT invalid"
         stop
        endif

        Tforce=IEforce*dt*one_over_DeDT
    
        vischeat(D_DECL(i,j,k))=vischeat(D_DECL(i,j,k))+Tforce

       else
        print *,"local_mask invalid"
        stop
       endif

      enddo
      enddo
      enddo ! i,j,k
 
      return
      end subroutine FORT_TENSORHEAT


         ! rhoinverse is 1/den
      subroutine FORT_VISCTENSORHEAT( &
       ntensor, &
       nsolve, &
       nstate, &
       xlo,dx,  &
       lsfab,DIMS(lsfab), &
       DeDTinverse, &
       DIMS(DeDTinverse), &
       vischeat,DIMS(vischeat), &
       xstress,DIMS(xstress), &
       ystress,DIMS(ystress), &
       zstress,DIMS(zstress), &
       gradu,DIMS(gradu), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       dt,irz, &
       nmat,nden)
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T ntensor
      INTEGER_T nsolve
      INTEGER_T nmat
      INTEGER_T nden,nstate,level
      REAL_T xlo(SDIM),dx(SDIM)
      INTEGER_T i,j,k
      INTEGER_T DIMDEC(lsfab)
      INTEGER_T DIMDEC(DeDTinverse)
      INTEGER_T DIMDEC(vischeat)
      INTEGER_T DIMDEC(xstress)
      INTEGER_T DIMDEC(ystress)
      INTEGER_T DIMDEC(zstress)
      INTEGER_T DIMDEC(gradu)
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T bfact
      REAL_T lsfab(DIMV(lsfab),nmat)
      REAL_T DeDTinverse(DIMV(DeDTinverse),nmat+1)
      REAL_T vischeat(DIMV(vischeat))
      REAL_T xstress(DIMV(xstress),nsolve)
      REAL_T ystress(DIMV(ystress),nsolve)
      REAL_T zstress(DIMV(zstress),nsolve)
      REAL_T gradu(DIMV(gradu),ntensor)
      REAL_T dt
      INTEGER_T irz
      INTEGER_T veldir,dir
      REAL_T IEforce,Tforce
      REAL_T one_over_DeDT
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf
      REAL_T local_gradu(3,3)
      REAL_T tensor(SDIM,SDIM)
      INTEGER_T ux,vx,wx,uy,vy,wy,uz,vz,wz,nbase
      INTEGER_T grdcomp

      nhalf=3

      if (bfact.lt.1) then
       print *,"bfact invalid51"
       stop
      endif
      if (nsolve.ne.SDIM) then
       print *,"nsolve invalid"
       stop
      endif
      if (ntensor.ne.SDIM*SDIM) then
       print *,"ntensor invalid"
       stop
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (irz.ne.levelrz) then
       print *,"irz invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (nstate.ne.(SDIM+1)+ &
          nmat*(num_state_material+ngeom_raw)+1) then
       print *,"nstate invalid"
       stop
      endif

      if (nden.ne.nmat*num_state_material) then
       print *,"nden invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(lsfab),1,-1,7)
      call checkbound(fablo,fabhi, &
       DIMS(DeDTinverse), &
       0,-1,7)
      call checkbound(fablo,fabhi,DIMS(vischeat),0,-1,7)
      call checkbound(fablo,fabhi,DIMS(xstress),0,0,7)
      call checkbound(fablo,fabhi,DIMS(ystress),0,1,7)
      call checkbound(fablo,fabhi,DIMS(zstress),0,SDIM-1,7)
      call checkbound(fablo,fabhi,DIMS(gradu),0,-1,7)

! u_x,v_x,w_x, u_y,v_y,w_y, u_z,v_z,w_z;  
      call tensorcomp_matrix(ux,uy,uz,vx,vy,vz,wx,wy,wz)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
       call gridsten_level(xsten,i,j,k,level,nhalf)

       if (irz.eq.0) then
        ! do nothing
       else if (irz.eq.1) then
        if (SDIM.ne.2) then
         print *,"dimension bust"
         stop
        endif
        if (xsten(0,1).le.zero) then
         print *,"no neg domain in r-z"
         stop
        endif
       else if (irz.eq.3) then
        if (xsten(-2,1).le.zero) then
         print *,"no neg domain in r-T"
         stop
        endif
       else
        print *,"irz invalid"
        stop
       endif

        ! E: div( u dot tau )=(u tau_11)_x+(u tau_12)_y+(u tau_13)_z+...
        ! e: grad u : tau

       do veldir=1,3
       do dir=1,3
        local_gradu(veldir,dir)=zero
       enddo 
       enddo 
       do dir=1,SDIM
        if (dir.eq.1) then
         nbase=ux-1
        else if (dir.eq.2) then
         nbase=uy-1
        else if ((dir.eq.SDIM).and.(SDIM.eq.3)) then
         nbase=uz-1
        else
         print *,"dir invalid viscoelastic heating"
         stop
        endif
        do veldir=1,SDIM
         grdcomp=nbase+veldir
         local_gradu(veldir,dir)=gradu(D_DECL(i,j,k),grdcomp)
        enddo ! veldir
       enddo ! dir=1..sdim

       do veldir=1,SDIM
        tensor(veldir,1)=half*(xstress(D_DECL(i,j,k),veldir)+ &
               xstress(D_DECL(i+1,j,k),veldir))
        tensor(veldir,2)=half*(ystress(D_DECL(i,j,k),veldir)+ &
               ystress(D_DECL(i,j+1,k),veldir))
        if (SDIM.eq.3) then
         tensor(veldir,SDIM)=half*(zstress(D_DECL(i,j,k),veldir)+ &
                zstress(D_DECL(i,j,k+1),veldir))
        endif
       enddo ! veldir

       IEforce=zero
       do veldir=1,SDIM
       do dir=1,SDIM
        IEforce=IEforce+local_gradu(veldir,dir)*tensor(veldir,dir)
       enddo
       enddo

       one_over_DeDT=DeDTinverse(D_DECL(i,j,k),1)  ! 1/(rho cv)

       if (one_over_DeDT.le.zero) then
        print *,"one_over_DeDT invalid"
        stop
       endif
       Tforce=-IEforce*dt*one_over_DeDT
       vischeat(D_DECL(i,j,k))= &
        vischeat(D_DECL(i,j,k))+Tforce

      enddo
      enddo
      enddo ! i,j,k
 
      return
      end subroutine FORT_VISCTENSORHEAT

      ! rhoinverse is 1/den
      ! curv(nten*(SDIM+5)+3+dir)=mgoni_force(dir)=
      !  (I-nn^T)(grad sigma)delta 
      ! masknbr:
      ! (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
      ! (2) =1 interior  =0 otherwise
      ! (3) =1 interior+ngrow-1  =0 otherwise
      ! (4) =1 interior+ngrow    =0 otherwise
      subroutine FORT_MARANGONIFORCE( &
       conservative_tension_force, &
       isweep, &
       nstate, &
       nten, &
       num_curv, &
       xlo,dx,  &
       facecut_index, &
       icefacecut_index, &
       curv_index, &
       pforce_index, &
       faceden_index, &  ! 1/rho
       icemask_index, &
       massface_index, &
       vofface_index, &
       ncphys, &
       xface,DIMS(xface), &
       yface,DIMS(yface), &
       zface,DIMS(zface), &
       xp,DIMS(xp), &
       yp,DIMS(yp), &
       zp,DIMS(zp), &
       maskcov,DIMS(maskcov), &
       masknbr,DIMS(masknbr), &
       vol,DIMS(vol), &
       areax,DIMS(areax), &
       areay,DIMS(areay), &
       areaz,DIMS(areaz), &
       xflux,DIMS(xflux), &
       yflux,DIMS(yflux), &
       zflux,DIMS(zflux), &
       vel,DIMS(vel), &
       den,DIMS(den), &
       ls,DIMS(ls), &
       lsho,DIMS(lsho), &
       rhoinverse, &
       DIMS(rhoinverse), &
       vof,DIMS(vof), &
       curv,DIMS(curv), &
       velnew,DIMS(velnew), &
       umac,DIMS(umac), &
       vmac,DIMS(vmac), &
       wmac,DIMS(wmac), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       bfact_grid, &
       level, &
       finest_level, &
       dt, &
       cur_time, &
       visc_coef, &
       presbc_in, &  
       velbc, &
       vofbc, &
       nmat,nden)
      use probcommon_module
      use probf90_module
      use global_utility_module
      use godunov_module
      use geometry_intersect_module
      use MOF_routines_module
      IMPLICIT NONE

      INTEGER_T conservative_tension_force
      INTEGER_T isweep
      INTEGER_T nstate
      INTEGER_T nten
      INTEGER_T num_curv
      INTEGER_T nmat
      INTEGER_T nden
      INTEGER_T level
      INTEGER_T finest_level
      INTEGER_T facecut_index
      INTEGER_T icefacecut_index
      INTEGER_T curv_index
      INTEGER_T pforce_index
      INTEGER_T faceden_index ! 1/rho
      INTEGER_T icemask_index
      INTEGER_T massface_index
      INTEGER_T vofface_index
      INTEGER_T ncphys
      REAL_T xlo(SDIM),dx(SDIM)
      INTEGER_T DIMDEC(xface)
      INTEGER_T DIMDEC(yface)
      INTEGER_T DIMDEC(zface)
      INTEGER_T DIMDEC(xp)
      INTEGER_T DIMDEC(yp)
      INTEGER_T DIMDEC(zp)
      INTEGER_T DIMDEC(maskcov)
      INTEGER_T DIMDEC(masknbr)
      INTEGER_T DIMDEC(vol)
      INTEGER_T DIMDEC(areax)
      INTEGER_T DIMDEC(areay)
      INTEGER_T DIMDEC(areaz)
      INTEGER_T DIMDEC(xflux)
      INTEGER_T DIMDEC(yflux)
      INTEGER_T DIMDEC(zflux)
      INTEGER_T DIMDEC(vel)
      INTEGER_T DIMDEC(den)
      INTEGER_T DIMDEC(ls)
      INTEGER_T DIMDEC(lsho)
      INTEGER_T DIMDEC(rhoinverse)
      INTEGER_T DIMDEC(vof)
      INTEGER_T DIMDEC(curv)
      INTEGER_T DIMDEC(velnew)
      INTEGER_T DIMDEC(umac)
      INTEGER_T DIMDEC(vmac)
      INTEGER_T DIMDEC(wmac)
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T bfact
      INTEGER_T bfact_grid
      REAL_T xface(DIMV(xface),ncphys)
      REAL_T yface(DIMV(yface),ncphys)
      REAL_T zface(DIMV(zface),ncphys)
      REAL_T xp(DIMV(xp),2*SDIM)
      REAL_T yp(DIMV(yp),2*SDIM)
      REAL_T zp(DIMV(zp),2*SDIM)
      REAL_T maskcov(DIMV(maskcov))
      REAL_T masknbr(DIMV(masknbr),4)
      REAL_T vol(DIMV(vol))
      REAL_T areax(DIMV(areax))
      REAL_T areay(DIMV(areay))
      REAL_T areaz(DIMV(areaz))
      REAL_T xflux(DIMV(xflux),SDIM)
      REAL_T yflux(DIMV(yflux),SDIM)
      REAL_T zflux(DIMV(zflux),SDIM)
      REAL_T ls(DIMV(ls),nmat*(SDIM+1))
      REAL_T lsho(DIMV(lsho),nmat*(SDIM+1))
      REAL_T rhoinverse(DIMV(rhoinverse),nmat+1)
      REAL_T vof(DIMV(vof),nmat*ngeom_recon)
      REAL_T curv(DIMV(curv),num_curv)
      REAL_T vel(DIMV(vel),SDIM)
      REAL_T den(DIMV(den),nden)
      REAL_T velnew(DIMV(velnew),SDIM)
      REAL_T umac(DIMV(umac))
      REAL_T vmac(DIMV(vmac))
      REAL_T wmac(DIMV(wmac))
      REAL_T dt,cur_time
      REAL_T visc_coef
      INTEGER_T presbc_in(SDIM,2) 
      INTEGER_T velbc(SDIM,2,SDIM)
      INTEGER_T vofbc(SDIM,2)
      REAL_T xsten(-3:3,SDIM)
      REAL_T xstenMAC(-3:3,SDIM)
      INTEGER_T nten_test,nhalf
      REAL_T LScen(nmat)
      REAL_T LS_star(nmat)
      INTEGER_T im,im_opp
      INTEGER_T im_trial,im_opp_trial
      INTEGER_T im_primary_left,im_primary_right
      INTEGER_T im_primary,im_secondary,im_tertiary
      INTEGER_T im_solid
      INTEGER_T im_prescribed
      INTEGER_T im_solid_valid
      INTEGER_T im_prescribed_valid
      INTEGER_T im_local
      INTEGER_T iforce
      INTEGER_T iten
      INTEGER_T iten_local
      INTEGER_T dirloc,dirtan
      INTEGER_T ix,jx
      INTEGER_T i,j,k
      INTEGER_T i_norm
      INTEGER_T i_norm_flux
      INTEGER_T local_i_norm
      INTEGER_T ic,jc,kc
      INTEGER_T ii,jj,kk
      INTEGER_T iii,jjj,kkk
      INTEGER_T iface,jface,kface
      INTEGER_T tcomp
      INTEGER_T itan,jtan
      INTEGER_T side_flag
      INTEGER_T sideface
      INTEGER_T sideflux
      INTEGER_T sidebc
      INTEGER_T sidecomp
      INTEGER_T inode,jnode
      INTEGER_T maskleft,maskright
      REAL_T x_array(-1:1,-1:1,0:1,SDIM)
      REAL_T xnode_array(0:1,0:1,SDIM)
      REAL_T xtension(SDIM)
      REAL_T xnode(SDIM)
      INTEGER_T im_primary_array(-1:1,-1:1,0:1)
      REAL_T LSint_array(-1:1,-1:1,0:1,nten)
      REAL_T LSmat_array(-1:1,-1:1,0:1,nmat)
      REAL_T thermal_array(-1:1,-1:1,0:1,nmat)
      REAL_T tension_array(-1:1,-1:1,0:1,nten)
      REAL_T LSmat_node(0:1,0:1,nmat)
      REAL_T thermal_node(0:1,0:1,nmat)
      REAL_T tension_node(0:1,0:1,nten)
      REAL_T LSint_node(0:1,0:1,nten)
      REAL_T LSpair_node(0:1,0:1,nten)
      REAL_T Nint_node(0:1,0:1,SDIM,nten)
      REAL_T Nmat_node(0:1,0:1,SDIM,nmat)
      REAL_T LSleftPC(nmat)
      REAL_T LSleft_face(nmat)
      REAL_T LSrightPC(nmat)
      REAL_T LSright_face(nmat)
      REAL_T LSmat(nmat)
      REAL_T LSmat_tess(nmat)
      REAL_T thermalmat(nmat)
      REAL_T LSint(nten)
      REAL_T wtsum,wt,xtest,signterm
      REAL_T wtsum_deriv(SDIM)
      REAL_T xlo_sten(SDIM)
      REAL_T xhi_sten(SDIM)
      REAL_T dx_sten(SDIM)
      REAL_T RR
      REAL_T nhold(SDIM)
      REAL_T mag,mag12
      INTEGER_T Nmat_node_valid(nmat)
      INTEGER_T Nint_node_valid(nten)
      INTEGER_T LSint_plus(nten)
      INTEGER_T LSint_minus(nten)
      INTEGER_T LSpair_plus(nten)
      INTEGER_T LSpair_minus(nten)
      REAL_T LStest
      INTEGER_T use_DCA
      REAL_T dotprod,udotn
      REAL_T nsolid(SDIM)
      REAL_T nproject(SDIM)
      REAL_T nghost(SDIM)
      REAL_T nperp(SDIM)
      REAL_T nface(SDIM)
      REAL_T liquid_viscosity
      REAL_T cos_angle,sin_angle
      INTEGER_T iten_13,iten_23
      REAL_T user_tension(nten)
      INTEGER_T nside
      INTEGER_T ncrossing
      REAL_T xcrossing(4,SDIM)
      INTEGER_T inode1,inode2,jnode1,jnode2
      REAL_T maxdist
      INTEGER_T cross1max,cross2max
      INTEGER_T icross1,icross2
      REAL_T x1(SDIM)
      REAL_T x2(SDIM)
      REAL_T abs_x1,abs_x2
      REAL_T xmid(SDIM)
      REAL_T LS1,LS2,frac
      REAL_T ten12
      REAL_T normal12(SDIM)
      REAL_T t_interface(SDIM)
      REAL_T areaface_MAC
      REAL_T delta_face,delta_mgoni
      REAL_T areaface_slice
      REAL_T areaface_slice_check
      REAL_T dx_node(2)
      REAL_T dxmax
      REAL_T dxmaxLS
      INTEGER_T curv_sanity_check,debug_curv
      REAL_T velloc
      REAL_T nexpect(SDIM)
      REAL_T texpect(SDIM)
      INTEGER_T veldir,dir_flux
      INTEGER_T local_masknbr
      INTEGER_T at_rz_center
      INTEGER_T at_reflect_wall
      INTEGER_T at_wall
      INTEGER_T at_ext_wall
      INTEGER_T is_solid_face
      INTEGER_T is_prescribed_face
      INTEGER_T partid_solid
      INTEGER_T partid_prescribed
      INTEGER_T klo_sten,khi_sten
      INTEGER_T FVM_surface_tension
      REAL_T AFACE
      REAL_T AFACE_ICE
      REAL_T massface,masscell,volface,denface,dencell
      REAL_T hx
      REAL_T hx_hoop
      REAL_T local_xcen,local_xleft,local_xright,local_dx
      REAL_T theta,theta_crit
      REAL_T vol_sten,volpos,facearea
      REAL_T cen_sten(SDIM)
      REAL_T cenpos(SDIM)
      REAL_T areacentroid(SDIM)
      REAL_T gradh
      REAL_T tension_scaled
      REAL_T nfluid_clsvof(SDIM)
      REAL_T nfluid(SDIM)
      INTEGER_T lsnormal_valid(2)
      REAL_T ls_intercept(2)
      REAL_T mass_local(2)
      REAL_T vol_local(2)
      REAL_T den_local(2)
      REAL_T fside
      REAL_T ffaceten(2)
      REAL_T surface_tension_force(SDIM)
      REAL_T mgoni_force(SDIM)
      REAL_T mgoni_force_local(SDIM)
      REAL_T hoop_force(SDIM)
      REAL_T force_hold
      INTEGER_T hoop_force_iten
      INTEGER_T hoop_force_ok
      REAL_T mgoni_temp(nmat)
      REAL_T gradT(SDIM)
      REAL_T forcelocal
      REAL_T force_mag
      REAL_T local_face(ncphys)
      REAL_T LS_save(D_DECL(-1:1,-1:1,-1:1))
      REAL_T localMAC_LS(D_DECL(-1:1,-1:1,-1:1),nmat)
      REAL_T localMAC_thermal(D_DECL(-1:1,-1:1,-1:1),nmat)
      REAL_T localMAC_velocity(SDIM)
     
      INTEGER_T im_solid_node
      INTEGER_T im_primary_cell
      REAL_T test_dist,im_solid_dist
      REAL_T vel_for_CL_model(SDIM)
      INTEGER_T imCL
      REAL_T denCL,T_CL

       ! sanity check is sensitive to FACETOL_DVOL used in LEVELSTRIP
       ! see PROBCOMMON.F90 for other possible sensitivities.
      curv_sanity_check=0
      debug_curv=0
 
      nhalf=3

      if (bfact.lt.1) then
       print *,"bfact invalid52"
       stop
      endif
      if (bfact_grid.lt.4) then
       print *,"bfact_grid invalid"
       stop
      endif

      do dirloc=1,SDIM
       if ((fablo(dirloc)/bfact_grid)*bfact_grid.ne.fablo(dirloc)) then
        print *,"fablo mod bfact_grid not 0"
        stop
       endif
       if (((fabhi(dirloc)+1)/bfact_grid)*bfact_grid.ne.fabhi(dirloc)+1) then
        print *,"fabhi+1 mod bfact_grid not 0"
        stop
       endif
      enddo ! dirloc=1..sdim

      ! height function curvature
      ! finite difference curvature
      ! pforce
      ! marangoni force
      ! dir/side flag
      ! im3
      ! x nten
      if (num_curv.ne.nten*(5+SDIM)) then
       print *,"num_curv invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid turbulenceforce nten, nten_test ",nten,nten_test
       stop
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (nstate.ne.(SDIM+1)+ &
          nmat*(num_state_material+ngeom_raw)+1) then
       print *,"nstate invalid"
       stop
      endif

      ! indexes start at 0
      if ((curv_index.ne.0).or. &
          (pforce_index.ne.1).or. &
          (facecut_index.ne.3).or. &
          (icefacecut_index.ne.4).or. &
          (icemask_index.ne.5).or. &
          (faceden_index.ne.2).or. &
          (vofface_index.ne.massface_index+2*nmat)) then
       print *,"face_index bust 2"
       stop
      endif
      if (ncphys.ne.vofface_index+2*nmat) then
       print *,"ncphys invalid"
       stop
      endif

       ! -1 if use static angle
      call get_use_DCA(use_DCA)

      if (dt.le.zero) then
       print *,"dt invalid"
       stop
      endif
      if (cur_time.lt.zero) then
       print *,"cur_time invalid"
       stop
      endif
      if (visc_coef.ge.zero) then
       ! do nothing
      else
       print *,"visc_coef invalid"
       stop
      endif
      if (nden.ne.nmat*num_state_material) then
       print *,"nden invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid marangoni force"
       stop
      endif
      if (finest_level.ne.fort_finest_level) then
       print *,"finest_level invalid marangoni force"
       stop
      endif
      if ((conservative_tension_force.ne.0).and. &
          (conservative_tension_force.ne.1)) then
       print *,"conservative_tension_force invalid"
       stop
      endif
      if ((isweep.ne.0).and.(isweep.ne.1)) then
       print *,"isweep invalid"
       stop
      endif

      if (SDIM.eq.2) then
       klo_sten=0
       khi_sten=0
      else if (SDIM.eq.3) then
       klo_sten=-1
       khi_sten=1
      else
       print *,"dimension bust"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(maskcov),2,-1,7)
      call checkbound(fablo,fabhi,DIMS(masknbr),2,-1,7)
      call checkbound(fablo,fabhi,DIMS(vol),2,-1,7)

      call checkbound(fablo,fabhi,DIMS(areax),2,0,7)
      call checkbound(fablo,fabhi,DIMS(areay),2,1,7)
      call checkbound(fablo,fabhi,DIMS(areaz),2,SDIM-1,7)

      call checkbound(fablo,fabhi,DIMS(xp),0,0,33)
      call checkbound(fablo,fabhi,DIMS(yp),0,1,33)
      call checkbound(fablo,fabhi,DIMS(zp),0,SDIM-1,33)

      call checkbound(fablo,fabhi,DIMS(umac),0,0,33)
      call checkbound(fablo,fabhi,DIMS(vmac),0,1,33)
      call checkbound(fablo,fabhi,DIMS(wmac),0,SDIM-1,33)

      call checkbound(fablo,fabhi,DIMS(xflux),0,0,7)
      call checkbound(fablo,fabhi,DIMS(yflux),0,1,7)
      call checkbound(fablo,fabhi,DIMS(zflux),0,SDIM-1,7)

      call checkbound(fablo,fabhi,DIMS(xface),0,0,33)
      call checkbound(fablo,fabhi,DIMS(yface),0,1,33)
      call checkbound(fablo,fabhi,DIMS(zface),0,SDIM-1,33)

      call checkbound(fablo,fabhi,DIMS(ls),2,-1,7)
      call checkbound(fablo,fabhi,DIMS(lsho),2,-1,7)
      call checkbound(fablo,fabhi, &
       DIMS(rhoinverse), &
       1,-1,7)
      call checkbound(fablo,fabhi,DIMS(vof),1,-1,7)
      call checkbound(fablo,fabhi,DIMS(curv),1,-1,7)
      call checkbound(fablo,fabhi,DIMS(vel),2,-1,7)
      call checkbound(fablo,fabhi,DIMS(den),2,-1,7)
      call checkbound(fablo,fabhi,DIMS(velnew),1,-1,7)

      call get_dxmax(dx,bfact,dxmax)
      call get_dxmaxLS(dx,bfact,dxmaxLS)

       ! isweep==0: surface tension force on MAC grid
       ! isweep==1: surface tension force on CELL grid
      if (isweep.eq.0) then

       if (conservative_tension_force.eq.1) then

        do veldir=0,SDIM-1

         ii=0
         jj=0
         kk=0
         if (veldir.eq.0) then
          ii=1
         else if (veldir.eq.1) then
          jj=1
         else if ((veldir.eq.2).and.(SDIM.eq.3)) then
          kk=1
         else
          print *,"veldir invalid"
          stop
         endif

         call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
          growlo,growhi,0,veldir,22)

         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)
          ! mask=1 if not covered by level+1 or outside the domain.
          maskleft=NINT(maskcov(D_DECL(i-ii,j-jj,k-kk)))
          maskright=NINT(maskcov(D_DECL(i,j,k)))
          if ((maskleft.eq.1).and.(maskright.eq.1)) then

           call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,veldir,28)
           hx=xstenMAC(1,veldir+1)-xstenMAC(-1,veldir+1)

           sidebc=0  ! non zero if EXT_DIR BC
           side_flag=0 ! non zero if at the grid boundary.
           if (veldir.eq.0) then
            i_norm=i
           else if (veldir.eq.1) then
            i_norm=j
           else if ((veldir.eq.2).and.(SDIM.eq.3)) then
            i_norm=k
           else
            print *,"veldir invalid"
            stop
           endif
 
           if (i_norm.eq.fablo(veldir+1)) then
            side_flag=1
            local_masknbr=NINT(masknbr(D_DECL(i-ii,j-jj,k-kk),1))
            if (vofbc(veldir+1,side_flag).eq.EXT_DIR) then
             sidebc=side_flag
             if (local_masknbr.ne.0) then
              print *,"local_masknbr.ne.0"
              stop
             endif
            endif
           else if (i_norm.eq.fabhi(veldir+1)+1) then
            side_flag=2
            local_masknbr=NINT(masknbr(D_DECL(i,j,k),1))
            if (vofbc(veldir+1,side_flag).eq.EXT_DIR) then
             sidebc=side_flag
             if (local_masknbr.ne.0) then
              print *,"local_masknbr.ne.0"
              stop
             endif
            endif
           else if ((i_norm.gt.fablo(veldir+1)).and. &
                    (i_norm.lt.fabhi(veldir+1)+1)) then
            ! do nothing
           else
            print *,"i_norm invalid"
            stop
           endif

           do im=1,nmat
            LSleftPC(im)=ls(D_DECL(i-ii,j-jj,k-kk),im)
            LSrightPC(im)=ls(D_DECL(i,j,k),im)
           enddo
           call get_primary_material(LSleftPC,nmat,im_primary_left)
           call get_primary_material(LSrightPC,nmat,im_primary_right)

           do dirloc=1,SDIM
            surface_tension_force(dirloc)=zero
            mgoni_force(dirloc)=zero
            mgoni_force_local(dirloc)=zero
            hoop_force(dirloc)=zero ! -S22/r, S12/r, 0
           enddo
           hoop_force_iten=0
           hoop_force_ok=0

           areaface_slice_check=zero

           if ((is_rigid(nmat,im_primary_left).eq.1).or. &
               (is_rigid(nmat,im_primary_right).eq.1)) then
            ! do nothing
           else if ((is_rigid(nmat,im_primary_left).eq.0).and. &
                    (is_rigid(nmat,im_primary_right).eq.0)) then

            at_rz_center=0
            if (levelrz.eq.0) then
             ! do nothing
            else if (levelrz.eq.1) then
             if (SDIM.eq.2) then
              if (veldir.eq.0) then
               if (i.eq.0) then
                if (side_flag.eq.1) then
                 if (vofbc(veldir+1,side_flag).eq.REFLECT_EVEN) then
                  at_rz_center=1
                 else
                  print *,"vofbc invalid"
                  stop
                 endif
                else
                 print *,"side_flag invalid"
                 stop
                endif
               else if ((i.ge.1).and.(i.le.fabhi(veldir+1)+1)) then
                ! do nothing
               else
                print *,"i invalid"
                stop
               endif
              else if (veldir.eq.SDIM-1) then
               ! do nothing
              else
               print *,"veldir invalid"
               stop
              endif
             else
              print *,"dimension bust"
              stop
             endif
            else if (levelrz.eq.3) then
             if (xstenMAC(0,1).gt.zero) then
              ! do nothing
             else
              print *,"xstenMAC(0,1) invalid"
              stop
             endif
            else
             print *,"levelrz invalid"
             stop
            endif
              
            if (at_rz_center.eq.1) then
             ! do nothing
            else if (at_rz_center.eq.0) then 

             call Box_volumeFAST(bfact,dx,xstenMAC,nhalf,vol_sten, &
               cen_sten,SDIM)

             if (veldir.eq.0) then
              areaface_slice=areax(D_DECL(i,j,k))
             else if (veldir.eq.1) then
              areaface_slice=areay(D_DECL(i,j,k))
             else if ((veldir.eq.2).and.(SDIM.eq.3)) then
              areaface_slice=areaz(D_DECL(i,j,k))
             else
              print *,"veldir invalid"
              stop
             endif

             do im=1,ncphys
              if (veldir.eq.0) then
               local_face(im)=xface(D_DECL(i,j,k),im)
              else if (veldir.eq.1) then
               local_face(im)=yface(D_DECL(i,j,k),im)
              else if ((veldir.eq.2).and.(SDIM.eq.3)) then
               local_face(im)=zface(D_DECL(i,j,k),im)
              else
               print *,"veldir invalid"
               stop
              endif
             enddo ! im=1..ncphys

             AFACE=local_face(facecut_index+1)
             if ((AFACE.ge.zero).and.(AFACE.le.half)) then
              AFACE=zero
             else if ((AFACE.ge.half).and.(AFACE.le.one)) then
              ! do nothing
             else
              print *,"AFACE invalid"
              stop
             endif

             AFACE_ICE=local_face(icefacecut_index+1)

             if ((AFACE.ge.zero).and. &
                 (AFACE.le.one).and. &
                 (AFACE_ICE.ge.zero).and. &
                 (AFACE_ICE.le.one)) then

              ! sideface=1 right half of cell that is left of the face.
              ! sideface=2 left half of cell that is right of the face.
              massface=zero
              volface=zero
              denface=zero

              do sideface=1,2
               mass_local(sideface)=zero
               vol_local(sideface)=zero
               do im=1,nmat
                mass_local(sideface)=mass_local(sideface)+ &
                 local_face(massface_index+2*(im-1)+sideface)
                vol_local(sideface)=vol_local(sideface)+ &
                 local_face(vofface_index+2*(im-1)+sideface)
               enddo 
               if (mass_local(sideface).ge.zero) then
                if (vol_local(sideface).ge.zero) then
                 if (vol_local(sideface).eq.zero) then
                  den_local(sideface)=zero
                 else
                  den_local(sideface)=mass_local(sideface)/vol_local(sideface)
                 endif
                 massface=massface+mass_local(sideface)
                 volface=volface+vol_local(sideface)
                 denface=denface+den_local(sideface)
                else
                 print *,"vol_local(sideface) invalid"
                 stop
                endif
               else
                print *,"mass(sideface) invalid"
                stop
               endif
              enddo  ! sideface=1,2

              call fixed_face( &
               nmat, &
               AFACE, &
               AFACE, &
               LSleftPC,LSrightPC, &
               is_solid_face, &
               is_prescribed_face, &
               im_solid, &
               im_prescribed, &
               im_solid_valid, &
               im_prescribed_valid, &
               partid_solid, &
               partid_prescribed) 

              at_reflect_wall=0
              at_wall=0
              at_ext_wall=0
              if ((side_flag.eq.1).or.(side_flag.eq.2)) then
               if ((presbc_in(veldir+1,side_flag).eq.REFLECT_EVEN).or. &
                   (presbc_in(veldir+1,side_flag).eq.FOEXTRAP)) then
                at_wall=1
               else if ((presbc_in(veldir+1,side_flag).eq.INT_DIR).or. &
                        (presbc_in(veldir+1,side_flag).eq.EXT_DIR)) then
                ! do nothing
               else
                print *,"presbc_in invalid"
                stop
               endif
        
               if (velbc(veldir+1,side_flag,veldir+1).eq.REFLECT_ODD) then
                at_reflect_wall=side_flag
               else if (velbc(veldir+1,side_flag,veldir+1).eq.EXT_DIR) then
                at_ext_wall=side_flag
               else if (velbc(veldir+1,side_flag,veldir+1).eq.INT_DIR) then
                ! do nothing
               else if (velbc(veldir+1,side_flag,veldir+1).eq.FOEXTRAP) then
                ! do nothing
               else if (velbc(veldir+1,side_flag,veldir+1).eq.REFLECT_EVEN) then
                ! do nothing
               else
                print *,"velbc invalid"
                stop
               endif
              else if (side_flag.eq.0) then
               ! do nothing
              else
               print *,"side_flag invalid"
               stop
              endif

               ! note: expecting AFACE = 0 at wall.
               ! at_wall==1 if presbc==REFLECT_EVEN or FOEXTRAP
              if (at_wall.eq.1) then
               is_solid_face=1
              else if (at_reflect_wall.eq.1) then
               is_solid_face=1
              else if (at_reflect_wall.eq.2) then
               is_solid_face=1
              else if ((at_reflect_wall.eq.0).and. &
                       (at_wall.eq.0)) then
               ! do nothing
              else
               print *,"bad:at_wall or at_reflect_wall"
               stop
              endif

              if (is_solid_face.eq.0) then

               do dirloc=1,SDIM
                vel_for_CL_model(dirloc)=zero
               enddo
               wtsum=zero

               do ic=-1,1
               do jc=-1,1
               do kc=klo_sten,khi_sten
                if (veldir.eq.0) then
                 local_i_norm=ic
                else if (veldir.eq.1) then
                 local_i_norm=jc
                else if ((veldir.eq.2).and.(SDIM.eq.3)) then
                 local_i_norm=kc
                else
                 print *,"veldir invalid"
                 stop
                endif
               
                local_xcen=xstenMAC(2*local_i_norm,veldir+1)
                local_xleft=xstenMAC(2*local_i_norm-1,veldir+1)
                local_xright=xstenMAC(2*local_i_norm+1,veldir+1)
                local_dx=local_xright-local_xleft
                if (local_dx.gt.zero) then
                 theta=local_xcen-local_xleft
                 if ((theta.gt.zero).and.(theta.lt.local_dx)) then
                  theta=theta/local_dx
                  do im=1,nmat
                   LSmat(im)= &
                    (one-theta)*lsho(D_DECL(i+ic-ii,j+jc-jj,k+kc-kk),im)+ &
                    theta*lsho(D_DECL(i+ic,j+jc,k+kc),im)
                   tcomp=(im-1)*num_state_material+2
                   localMAC_thermal(D_DECL(ic,jc,kc),im)= &
                    (one-theta)*den(D_DECL(i+ic-ii,j+jc-jj,k+kc-kk),tcomp)+ &
                    theta*den(D_DECL(i+ic,j+jc,k+kc),tcomp)
                  enddo ! im=1..nmat
                  call FIX_LS_tessellate(LSmat,LSmat_tess,nmat)
                  do im=1,nmat
                   localMAC_LS(D_DECL(ic,jc,kc),im)=LSmat_tess(im)
                  enddo
                  call get_primary_material(LSmat_tess,nmat,im_primary_cell)
                  do dirloc=1,SDIM
                   localMAC_velocity(dirloc)= &
                    (one-theta)*vel(D_DECL(i+ic-ii,j+jc-jj,k+kc-kk),dirloc)+ &
                    theta*vel(D_DECL(i+ic,j+jc,k+kc),dirloc)
                  enddo
                  if (is_rigid(nmat,im_primary_cell).eq.1) then
                   wt=zero
                  else if (is_rigid(nmat,im_primary_cell).eq.0) then
                   wt=one
                  else
                   print *,"is_rigid(nmat,im_primary_cell) invalid"
                   stop
                  endif
                  wtsum=wtsum+wt
                  do dirloc=1,SDIM
                   vel_for_CL_model(dirloc)=vel_for_CL_model(dirloc)+ &
                     wt*localMAC_velocity(dirloc)
                  enddo
                 else
                  print *,"theta invalid"
                  stop
                 endif
                else
                 print *,"local_dx invalid"
                 stop
                endif
               enddo ! kc=klo_sten,khi_sten
               enddo ! jc=-1,1
               enddo ! ic=-1,1

               if (wtsum.gt.zero) then
                do dirloc=1,SDIM
                 vel_for_CL_model(dirloc)=vel_for_CL_model(dirloc)/wtsum
                enddo
               else if (wtsum.eq.zero) then
                ! do nothing
               else
                print *,"wtsum invalid"
                stop
               endif

               FVM_surface_tension=1

               if ((level.ge.0).and.(level.lt.finest_level)) then
                FVM_surface_tension=0
               else if (level.eq.finest_level) then
                ic=i
                jc=j
                kc=k
                if ((side_flag.eq.1).or.(side_flag.eq.2)) then

                 if (side_flag.eq.1) then
                  ! do nothing, cell (ic,jc,kc) in the domain
                 else if (side_flag.eq.2) then
                  ic=ic-ii 
                  jc=jc-jj 
                  kc=kc-kk
                 else
                  print *,"side_flag invalid"
                  stop
                 endif 

                 if ((vofbc(veldir+1,side_flag).eq.EXT_DIR).or. &
                     (vofbc(veldir+1,side_flag).eq.REFLECT_EVEN).or. &
                     (vofbc(veldir+1,side_flag).eq.FOEXTRAP)) then
                  FVM_surface_tension=0
                 else if (vofbc(veldir+1,side_flag).eq.INT_DIR) then
                  if (local_masknbr.eq.0) then ! coarse/fine
                   FVM_surface_tension=0
                  else if (local_masknbr.eq.1) then ! fine/fine
                   ! do nothing
                  else
                   print *,"local_masknbr invalid"
                   stop
                  endif
                 else
                  print *,"vofbc(veldir+1,side_flag) invalid"
                  stop
                 endif
                else if (side_flag.eq.0) then
                 ! do nothing
                else
                 print *,"side_flag invalid"
                 stop
                endif

                do dir_flux=0,SDIM-1

                 if (dir_flux.eq.veldir) then
                  ! do nothing (already checked the boundary conditions in the
                  ! veldir direction up above)
                 else if (dir_flux.ne.veldir) then

                  iii=0
                  jjj=0
                  kkk=0
                  if (dir_flux.eq.0) then
                   iii=1
                   i_norm_flux=i
                  else if (dir_flux.eq.1) then
                   jjj=1
                   i_norm_flux=j
                  else if ((dir_flux.eq.2).and.(SDIM.eq.3)) then
                   kkk=1
                   i_norm_flux=k
                  else
                   print *,"dir_flux invalid"
                   stop
                  endif

                  if (i_norm_flux.eq.fablo(dir_flux+1)) then
                   if ((vofbc(dir_flux+1,1).eq.EXT_DIR).or. &
                       (vofbc(dir_flux+1,1).eq.REFLECT_EVEN).or. &
                       (vofbc(dir_flux+1,1).eq.FOEXTRAP)) then
                    FVM_surface_tension=0
                   else if (vofbc(dir_flux+1,1).eq.INT_DIR) then
                    local_masknbr=NINT(masknbr(D_DECL(ic-iii,jc-jjj,kc-kkk),1))
                    if (local_masknbr.eq.0) then ! coarse/fine
                     FVM_surface_tension=0
                    else if (local_masknbr.eq.1) then ! fine/fine
                     ! do nothing
                    else
                     print *,"local_masknbr invalid"
                     stop
                    endif
                   else
                    print *,"vofbc(dir_flux+1,1) invalid"
                    print *,"vofbc(dir_flux+1,1)=",vofbc(dir_flux+1,1)
                    print *,"dir_flux=",dir_flux
                    print *,"i_norm_flux=",i_norm_flux
                    print *,"fablo(dir_flux+1)=",fablo(dir_flux+1)
                    stop
                   endif
                  else if ((i_norm_flux.gt.fablo(dir_flux+1)).and. &
                           (i_norm_flux.le.fabhi(dir_flux+1))) then
                   ! do nothing
                  else
                   print *,"i_norm_flux invalid"
                   stop
                  endif
                  if (i_norm_flux.eq.fabhi(dir_flux+1)) then
                   if ((vofbc(dir_flux+1,2).eq.EXT_DIR).or. &
                       (vofbc(dir_flux+1,2).eq.REFLECT_EVEN).or. &
                       (vofbc(dir_flux+1,2).eq.FOEXTRAP)) then
                    FVM_surface_tension=0
                   else if (vofbc(dir_flux+1,2).eq.INT_DIR) then
                    local_masknbr=NINT(masknbr(D_DECL(ic+iii,jc+jjj,kc+kkk),1))
                    if (local_masknbr.eq.0) then ! coarse/fine
                     FVM_surface_tension=0
                    else if (local_masknbr.eq.1) then ! fine/fine
                     ! do nothing
                    else
                     print *,"local_masknbr invalid"
                     stop
                    endif
                   else
                    print *,"vofbc(dir_flux+1,2) invalid"
                    stop
                   endif
                  else if ((i_norm_flux.ge.fablo(dir_flux+1)).and. &
                           (i_norm_flux.lt.fabhi(dir_flux+1))) then
                   ! do nothing
                  else
                   print *,"i_norm_flux invalid"
                   stop
                  endif

                 else
                  print *,"dir_flux invalid"
                  stop
                 endif

                enddo ! dir_flux=0..sdim-1

               else
                print *,"level invalid 32"
                stop
               endif

               do im=1,nmat
                mgoni_temp(im)=localMAC_thermal(D_DECL(0,0,0),im)
               enddo 
               do dirloc=1,SDIM
                xtension(dirloc)=xstenMAC(0,dirloc)
               enddo
               call get_user_tension(xtension,cur_time, &
                fort_tension,user_tension, &
                mgoni_temp,nmat,nten,4)

               do im=1,nmat
                LSmat(im)=localMAC_LS(D_DECL(0,0,0),im)
               enddo

               theta_crit=-one
               call get_primary_material(LSmat,nmat,im)
               if (is_rigid(nmat,im).eq.0) then
                do ic=-1,1
                do jc=-1,1
                do kc=klo_sten,khi_sten
                 if ((abs(ic).eq.1).or. &
                     (abs(jc).eq.1).or. &
                     (abs(kc).eq.1)) then
                  do im_trial=1,nmat
                   LS_star(im_trial)= &
                    localMAC_LS(D_DECL(ic,jc,kc),im_trial) 
                  enddo
                  call get_primary_material(LS_star,nmat,im_opp_trial)
                  if (is_rigid(nmat,im_opp_trial).eq.0) then
                   if (im.ne.im_opp_trial) then
                    if ((LSmat(im).eq.zero).and. &
                        (LS_star(im_opp_trial).eq.zero)) then
                     theta_crit=zero
                     im_opp=im_opp_trial
                    else if ((LSmat(im).gt.zero).or. &
                             (LS_star(im_opp_trial).gt.zero)) then
                     if ((LSmat(im).ge.zero).and. &
                         (LS_star(im_opp_trial).ge.zero)) then
                      theta=LSmat(im)/(LSmat(im)+LS_star(im_opp_trial))
                      mag=zero
                      do dirloc=1,SDIM
                       if (dirloc.eq.1) then
                        x1(dirloc)=xstenMAC(2*ic,dirloc)
                       else if (dirloc.eq.2) then
                        x1(dirloc)=xstenMAC(2*jc,dirloc)
                       else if ((dirloc.eq.3).and.(SDIM.eq.3)) then
                        x1(dirloc)=xstenMAC(2*kc,dirloc)
                       else
                        print *,"dirloc invalid"
                        stop
                       endif
                       mag=mag+(xtension(dirloc)-x1(dirloc))**2
                      enddo ! dirloc=1..sdim
                      mag=sqrt(mag)
                      theta=theta*mag
                      if (theta_crit.eq.-one) then
                       theta_crit=theta
                       im_opp=im_opp_trial
                      else if (theta_crit.ge.zero) then
                       if (theta.lt.theta_crit) then
                        theta_crit=theta
                        im_opp=im_opp_trial
                       endif
                      else
                       print *,"theta_crit invalid"
                       stop
                      endif
                     else
                      print *,"LSmat or LS_star invalid"
                      stop
                     endif
                    else
                     print *,"LSmat or LS_star invalid"
                     stop
                    endif
                   else if (im.eq.im_opp_trial) then
                    ! do nothing
                   else
                    print *,"im or im_opp_trial invalid"
                    stop
                   endif
                  else if (is_rigid(nmat,im_opp_trial).eq.1) then
                   ! do nothing
                  else
                   print *,"is_rigid(nmat,im_opp_trial) invalid"
                   stop
                  endif
                 else if ((ic.eq.0).and.(jc.eq.0).and.(kc.eq.0)) then
                  ! do nothing
                 else
                  print *,"ic,jc, or kc invalid"
                  stop
                 endif
                enddo ! kc=klo_sten,khi_sten
                enddo ! jc=-1,1
                enddo ! ic=-1,1
               else if (is_rigid(nmat,im).eq.1) then
                ! do nothing
               else
                print *,"is_rigid(nmat,im) invalid"
                stop
               endif

               if (theta_crit.ge.zero) then
                if ((im.ge.1).and.(im.le.nmat).and. &
                    (im_opp.ge.1).and.(im_opp.le.nmat).and. &
                    (im.ne.im_opp)) then

                 call get_iten(im,im_opp,hoop_force_iten,nmat)

                 do ic=-1,1
                 do jc=-1,1
                 do kc=klo_sten,khi_sten
                  do im_trial=1,nmat
                   LS_star(im_trial)=localMAC_LS(D_DECL(ic,jc,kc),im_trial)
                  enddo
                  call get_LS_extend(LS_star,nmat,hoop_force_iten, &
                   LS_save(D_DECL(ic,jc,kc)))
                 enddo ! kc=klo_sten,khi_sten
                 enddo ! jc=-1,1
                 enddo ! ic=-1,1

                 ! centroid in absolute coordinate system
                 ! returns a volume fraction
                 call getvolume( &
                  bfact,dx,xstenMAC,nhalf, &
                  LS_save,volpos,facearea, &
                  cenpos,areacentroid,VOFTOL,SDIM)

                 if (facearea.ge.zero) then
                  if (vol_sten.gt.zero) then
                   delta_mgoni=facearea/vol_sten
                  else
                   print *,"vol_sten invalid"
                   stop
                  endif
                 else
                  print *,"facearea invalid"
                  stop
                 endif

                 ! normal points from (-) to (+)
                 ! finds grad phi/|grad phi| where grad=(d/dx,d/dy,d/dz) or
                 ! grad=(d/dr,d/dz) or
                 ! grad=(d/dr,d/dtheta,d/dz)
                 im_local=1
                 call find_cut_geom_slope_CLSVOF( &
                   LS_save,nfluid_clsvof, &
                   lsnormal_valid, &
                   ls_intercept, & 
                   bfact,dx,xstenMAC,nhalf, &
                   im_local, &
                   dxmaxLS, &
                   im_local, &
                   SDIM) 
                 if (lsnormal_valid(1).eq.1) then
                  RR=xstenMAC(0,1)
                  call prepare_normal(nfluid_clsvof,RR,mag)
                  if (mag.gt.zero) then
                   ! Marangoni force:
                   ! (I-nn^T)(grad sigma) delta=
                   ! (grad sigma - (grad sigma dot n)n ) delta
                   if (abs(fort_tension_slope(hoop_force_iten)).le.1.0D+20) then
                    dotprod=zero

                    if (fort_tension_slope(hoop_force_iten).ne.zero) then

                     do dirloc=1,SDIM
                      ic=0
                      jc=0
                      kc=0
                      RR=one
                      if (dirloc.eq.1) then
                       ic=1
                       RR=one
                      else if (dirloc.eq.2) then 
                       jc=1
                       if ((levelrz.eq.0).or. &
                           (levelrz.eq.1)) then
                        RR=one
                       else if (levelrz.eq.3) then
                        RR=xstenMAC(0,1)
                        if (RR.le.zero) then
                         print *,"RR invalid"
                         stop
                        endif
                       else
                        print *,"gradT: levelrz invalid"
                        stop
                       endif
                      else if ((dirloc.eq.3).and.(SDIM.eq.3)) then
                       kc=1
                       RR=one
                      else
                       print *,"dirloc invalid"
                       stop
                      endif 

                      hx_hoop=xstenMAC(2,dirloc)-xstenMAC(-2,dirloc)
                      if (hx_hoop.gt.zero) then
                       gradT(dirloc)=( &
                        localMAC_thermal(D_DECL(ic,jc,kc),im)+ &
                        localMAC_thermal(D_DECL(ic,jc,kc),im_opp)- &
                        localMAC_thermal(D_DECL(-ic,-jc,-kc),im)- &
                        localMAC_thermal(D_DECL(-ic,-jc,-kc),im_opp))/ &
                        (two*RR*hx_hoop)
                      else
                       print *,"hx_hoop invalid"
                       stop 
                      endif

                      dotprod=dotprod+ &
                       gradT(dirloc)* &
                       fort_tension_slope(hoop_force_iten)*nfluid_clsvof(dirloc)
                     enddo ! dirloc=1..sdim

                     ! tension=sigma_0 + slope*(T-T0)
                     do dirloc=1,SDIM
                      mgoni_force_local(dirloc)= &
                       (fort_tension_slope(hoop_force_iten)*gradT(dirloc)- &
                        nfluid_clsvof(dirloc)*dotprod)*delta_mgoni
                     enddo
                    else if (fort_tension_slope(hoop_force_iten).eq.zero) then
                     ! do nothing
                    else
                     print *,"fort_tension_slope(hoop_force_iten) invalid"
                     stop
                    endif

                    RR=xstenMAC(0,1)
                    call get_scaled_tension(user_tension(hoop_force_iten), &
                      tension_scaled)
                    if (SDIM.eq.2) then
                     if (levelrz.eq.0) then
                      ! do nothing
                     else if (levelrz.eq.1) then
                      if (RR.gt.zero) then
                       hoop_force(1)=-tension_scaled*delta_mgoni/RR
                      else
                       print *,"RR invalid"
                       stop
                      endif
                     else if (levelrz.eq.3) then
                      if (RR.gt.zero) then
                       hoop_force(1)=-tension_scaled*delta_mgoni* &
                        (one-nfluid_clsvof(2)**2)/RR
                       hoop_force(2)=-tension_scaled*delta_mgoni* &
                        nfluid_clsvof(1)*nfluid_clsvof(2)/RR
                      else
                       print *,"RR invalid"
                       stop
                      endif
                     else
                      print *,"levelrz invalid"
                      stop
                     endif
                    else if (SDIM.eq.3) then
                     if (levelrz.eq.0) then
                      ! do nothing
                     else if (levelrz.eq.3) then
                      if (RR.gt.zero) then
                       hoop_force(1)=-tension_scaled*delta_mgoni* &
                        (one-nfluid_clsvof(2)**2)/RR
                       hoop_force(2)=-tension_scaled*delta_mgoni* &
                        nfluid_clsvof(1)*nfluid_clsvof(2)/RR
                      else
                       print *,"RR invalid"
                       stop
                      endif
                     else
                      print *,"levelrz invalid"
                      stop
                     endif
                    else
                     print *,"dimension bust"
                     stop
                    endif

                   else
                    print *,"fort_tension_slopes bust"
                    stop
                   endif 
                  else if (mag.eq.zero) then
                   ! do nothing
                  else
                   print *,"mag invalid 2; mag=",mag
                   stop
                  endif
                 else if (lsnormal_valid(1).eq.0) then
                  ! do nothing
                 else
                  print *,"lsnormal_valid(1) invalid"
                  stop
                 endif
                else if ((im.eq.0).and.(im_opp.eq.0)) then
                 ! do nothing
                else
                 print *,"im or im_opp invalid 1"
                 print *,"im= ",im
                 print *,"im_opp= ",im_opp
                 print *,"nmat= ",nmat
                 stop
                endif
               else if (theta_crit.eq.-one) then
                ! do nothing
               else
                print *,"theta_crit invalid"
                stop
               endif
   
               if (FVM_surface_tension.eq.0) then

                call fluid_interface_tension(LSleftPC,LSrightPC,gradh, &
                 im_opp,im,nmat,nten)

                if (gradh.ne.zero) then

                 call get_iten(im,im_opp,iten,nmat)
                 call get_scaled_tension(user_tension(iten),tension_scaled)

                  ! approximation to -sigma kappa grad h
                 surface_tension_force(veldir+1)= &
                   -tension_scaled*local_face(curv_index+1)*gradh/hx
         
                else if (gradh.eq.zero) then
                 ! do nothing
                else
                 print *,"gradh invalid"
                 stop
                endif

                do dirloc=1,SDIM
                 mgoni_force(dirloc)=mgoni_force_local(dirloc)
                enddo

               else if (FVM_surface_tension.eq.1) then

                do dir_flux=0,SDIM-1
                do sideflux=0,1

                 iii=0
                 jjj=0
                 kkk=0
                 if (dir_flux.eq.0) then
                  iii=1
                 else if (dir_flux.eq.1) then
                  jjj=1
                 else if ((dir_flux.eq.2).and.(SDIM.eq.3)) then
                  kkk=1
                 else
                  print *,"dir_flux invalid"
                  stop
                 endif

                  ! sideflux=0..1 
                 do itan=-1,1
                 do jtan=-1,1
                 do sideface=0,1
                  if (dir_flux.eq.0) then
                   ic=sideflux-iii+sideface
                   jc=itan
                   kc=jtan
                  else if (dir_flux.eq.1) then
                   ic=itan
                   jc=sideflux-jjj+sideface
                   kc=jtan
                  else if ((dir_flux.eq.2).and.(SDIM.eq.3)) then
                   ic=itan
                   jc=jtan
                   kc=sideflux-kkk+sideface
                  else
                   print *,"dir_flux invalid"
                   stop
                  endif
                  if (SDIM.eq.2) then
                   kc=0
                  else if (SDIM.eq.3) then
                   ! do nothing
                  else
                   print *,"dimension bust"
                   stop
                  endif
                  do dirloc=1,SDIM
                   if (dirloc.eq.1) then
                    ix=ic
                   else if (dirloc.eq.2) then
                    ix=jc
                   else if ((dirloc.eq.SDIM).and.(SDIM.eq.3)) then
                    ix=kc
                   else
                    print *,"dirloc invalid"
                    stop
                   endif
                   xtension(dirloc)=xstenMAC(2*ix,dirloc)
                   x_array(itan,jtan,sideface,dirloc)=xtension(dirloc)
                  enddo ! dirloc=1..sdim

                  do im=1,nmat
                   LSmat(im)=localMAC_LS(D_DECL(ic,jc,kc),im)
                   if (abs(itan)+abs(jtan).eq.0) then
                    if (sideface.eq.0) then
                     LSleft_face(im)=LSmat(im)
                    else if (sideface.eq.1) then
                     LSright_face(im)=LSmat(im)
                    else
                     print *,"sideface invalid"
                     stop
                    endif
                   else if (abs(itan)+abs(jtan).eq.1) then
                    ! do nothing
                   else if (abs(itan)+abs(jtan).eq.2) then
                    ! do nothing
                   else
                    print *,"itan or jtan invalid"
                    stop
                   endif
                   thermalmat(im)=localMAC_thermal(D_DECL(ic,jc,kc),im)
                   if (thermalmat(im).le.zero) then
                    print *,"temperature must be positive"
                    stop
                   else if (thermalmat(im).gt.zero) then
                    ! do nothing
                   else
                    print *,"thermalmat(im) invalid"
                    stop
                   endif 
                  enddo !im=1,nmat

                  call get_user_tension(xtension,cur_time, &
                   fort_tension,user_tension, &
                   thermalmat, &
                   nmat,nten,5)

                  do im=1,nmat
                   do im_opp=im+1,nmat
                    call get_iten(im,im_opp,iten,nmat)
                    if (is_rigid(nmat,im).eq.1) then
                     LSint(iten)=LSmat(im)
                    else if (is_rigid(nmat,im_opp).eq.1) then
                     LSint(iten)=-LSmat(im_opp)
                    else if ((is_rigid(nmat,im).eq.0).and. &
                             (is_rigid(nmat,im_opp).eq.0)) then
                     call get_LS_extend(LSmat,nmat,iten,LSint(iten)) 
                    else
                     print *,"is_rigid invalid"
                     stop
                    endif
                    LSint_array(itan,jtan,sideface,iten)=LSint(iten)
                    tension_array(itan,jtan,sideface,iten)=user_tension(iten)
                   enddo ! im_opp=im+1..nmat

                   LSmat_array(itan,jtan,sideface,im)=LSmat(im)
                   thermal_array(itan,jtan,sideface,im)=thermalmat(im)
                  enddo ! im=1..nmat
                  call get_primary_material(LSmat,nmat, &
                    im_primary_array(itan,jtan,sideface))
                 enddo ! sideface=0..1
                 enddo ! jtan=-1..1
                 enddo ! itan=-1..1

                 if (sideflux.eq.0) then
                  ix=-1
                 else if (sideflux.eq.1) then
                  ix=1
                 else
                  print *,"sideflux invalid"
                  stop
                 endif

                 xnode(dir_flux+1)=xstenMAC(ix,dir_flux+1)

                 do iten=1,nten
                  LSint_plus(iten)=0
                  LSint_minus(iten)=0
                  LSpair_plus(iten)=0
                  LSpair_minus(iten)=0
                 enddo ! iten=1..nten

                 do inode=0,1
                 do jnode=0,1
                  ix=2*inode-1 ! -1 or 1
                  jx=2*jnode-1 ! -1 or 1
                  if (dir_flux.eq.0) then
                   xnode(2)=xstenMAC(ix,2)
                   if (SDIM.eq.3) then
                    xnode(SDIM)=xstenMAC(jx,SDIM)
                   endif
                  else if (dir_flux.eq.1) then
                   xnode(1)=xstenMAC(ix,1)
                   if (SDIM.eq.3) then
                    xnode(SDIM)=xstenMAC(jx,SDIM)
                   endif
                  else if ((dir_flux.eq.2).and.(SDIM.eq.3)) then
                   xnode(1)=xstenMAC(ix,1)
                   xnode(2)=xstenMAC(jx,2)
                  else
                   print *,"dir_flux invalid"
                   stop
                  endif
                  ! in 2d if dir_flux==0,
                  ! 0,0 -> xface,ylo
                  ! 1,0 -> xface,yhi
                  ! 0,1 -> xface,ylo
                  ! 1,1 -> xface,yhi
                  ! in 2d if dir_flux==1,
                  ! 0,0 -> xlo,yface
                  ! 1,0 -> xhi,yface
                  ! 0,1 -> xlo,yface
                  ! 1,1 -> xhi,yface
                  do dirloc=1,SDIM
                   xnode_array(inode,jnode,dirloc)=xnode(dirloc)
                  enddo

                  do im=1,nmat
                   LSmat_node(inode,jnode,im)=zero
                   thermal_node(inode,jnode,im)=zero
                   do dirloc=1,SDIM
                    Nmat_node(inode,jnode,dirloc,im)=zero
                   enddo
                   do im_opp=im+1,nmat
                    call get_iten(im,im_opp,iten,nmat)
                    LSint_node(inode,jnode,iten)=zero
                    tension_node(inode,jnode,iten)=zero
                    do dirloc=1,SDIM
                     Nint_node(inode,jnode,dirloc,iten)=zero
                    enddo
                   enddo ! im_opp=im+1..nmat
                  enddo ! im=1..nmat

                  im_solid_node=0
                  im_solid_dist=zero
 
                  wtsum=zero
                  do dirloc=1,SDIM
                   wtsum_deriv(dirloc)=zero
                  enddo
                
                  do itan=-1+inode,inode
                  do jtan=-1+jnode,jnode
                  do sideface=0,1 
                   wt=zero
                   do dirloc=1,SDIM
                    wt=wt+(xnode(dirloc)-x_array(itan,jtan,sideface,dirloc))**2
                   enddo
                   if (wt.gt.zero) then
                    wt=one/wt
                   else
                    print *,"wt invalid"
                    stop
                   endif
                   wtsum=wtsum+wt
                   do im=1,nmat
                    LSmat_node(inode,jnode,im)=LSmat_node(inode,jnode,im)+ &
                     wt*LSmat_array(itan,jtan,sideface,im)
                    thermal_node(inode,jnode,im)=thermal_node(inode,jnode,im)+ &
                     wt*thermal_array(itan,jtan,sideface,im)
                    do im_opp=im+1,nmat
                     call get_iten(im,im_opp,iten,nmat)
                     LSint_node(inode,jnode,iten)= &
                      LSint_node(inode,jnode,iten)+ &
                      wt*LSint_array(itan,jtan,sideface,iten)
                     tension_node(inode,jnode,iten)= &
                      tension_node(inode,jnode,iten)+ &
                      wt*tension_array(itan,jtan,sideface,iten)
                    enddo ! im_opp=im+1..nmat
                   enddo ! im=1..nmat

                   do dirloc=1,SDIM !dirloc is normal component to be filled.

                    wt=zero
                    do dirtan=1,SDIM
                     if (dirtan.ne.dirloc) then
                      wt=wt+ &
                       (xnode(dirtan)-x_array(itan,jtan,sideface,dirtan))**2
                     endif
                    enddo
                    if (wt.gt.zero) then
                     wt=one/wt
                    else
                     print *,"wt invalid"
                     stop
                    endif
                    wtsum_deriv(dirloc)=wtsum_deriv(dirloc)+wt
         
                    xtest=x_array(itan,jtan,sideface,dirloc) 

                    if (dirloc.eq.dir_flux+1) then
                     if (sideface.eq.0) then
                      signterm=-one
                      xlo_sten(dirloc)=xtest
                     else if (sideface.eq.1) then
                      signterm=one
                      xhi_sten(dirloc)=xtest
                     else
                      print *,"sideface invalid"
                      stop
                     endif
                    else if (dirloc.ne.dir_flux+1) then
                     if (dirloc.eq.1) then ! x direction
                      if ((dir_flux.ne.1).and. &
                          (dir_flux.ne.SDIM-1)) then
                       print *,"dir_flux invalid"
                       stop
                      endif
                      if (itan.eq.-1+inode) then
                       signterm=-one
                       xlo_sten(dirloc)=xtest
                      else if (itan.eq.inode) then
                       signterm=one
                       xhi_sten(dirloc)=xtest
                      else
                       print *,"itan invalid"
                       stop
                      endif
                     else if (dirloc.eq.2) then ! y direction
                      if (dir_flux.eq.0) then
                       if (itan.eq.-1+inode) then
                        signterm=-one
                        xlo_sten(dirloc)=xtest
                       else if (itan.eq.inode) then
                        signterm=one
                        xhi_sten(dirloc)=xtest
                       else
                        print *,"itan invalid"
                        stop
                       endif
                      else if ((dir_flux.eq.SDIM-1).and.(SDIM.eq.3)) then
                       if (jtan.eq.-1+jnode) then
                        signterm=-one
                        xlo_sten(dirloc)=xtest
                       else if (jtan.eq.jnode) then
                        signterm=one
                        xhi_sten(dirloc)=xtest
                       else
                        print *,"jtan invalid"
                        stop
                       endif
                      else
                       print *,"dir_flux invalid"
                       stop
                      endif
                     else if ((dirloc.eq.3).and.(SDIM.eq.3)) then
                      if ((dir_flux.ne.0).and. &
                          (dir_flux.ne.1)) then
                       print *,"dir_flux invalid"
                       stop
                      endif
                      if (jtan.eq.-1+jnode) then
                       signterm=-one
                       xlo_sten(dirloc)=xtest
                      else if (jtan.eq.jnode) then
                       signterm=one
                       xhi_sten(dirloc)=xtest
                      else
                       print *,"jtan invalid"
                       stop
                      endif
                     else
                      print *,"dirloc invalid"
                      stop
                     endif
                    else
                     print *,"dirloc or dir_flux invalid"
                     stop
                    endif

                    do im=1,nmat
                     Nmat_node(inode,jnode,dirloc,im)= &
                      Nmat_node(inode,jnode,dirloc,im)+ &
                      signterm*wt*LSmat_array(itan,jtan,sideface,im)
                     do im_opp=im+1,nmat
                      call get_iten(im,im_opp,iten,nmat)
                      Nint_node(inode,jnode,dirloc,iten)= &
                       Nint_node(inode,jnode,dirloc,iten)+ &
                       signterm*wt*LSint_array(itan,jtan,sideface,iten)
                     enddo ! im_opp=im+1..nmat
                    enddo ! im=1..nmat

                   enddo ! dirloc=1..sdim

                   im_primary_cell=im_primary_array(itan,jtan,sideface)
                   if ((im_primary_cell.ge.1).and. &
                       (im_primary_cell.le.nmat))  then
                    if (is_rigid(nmat,im_primary_cell).eq.1) then
                     test_dist=LSmat_array(itan,jtan,sideface,im_primary_cell)
                     if (im_solid_node.eq.0) then
                      im_solid_node=im_primary_cell
                      im_solid_dist=test_dist
                     else if ((im_solid_node.ge.1).and. &
                              (im_solid_node.le.nmat)) then
                      if (test_dist.gt.im_solid_dist) then
                       im_solid_node=im_primary_cell
                       im_solid_dist=test_dist
                      else if (test_dist.le.im_solid_dist) then
                       ! do nothing
                      else
                       print *,"test_dist invalid"
                       stop
                      endif
                     else
                      print *,"im_solid_node invalid"
                      stop
                     endif
                    else if (is_rigid(nmat,im_primary_cell).eq.0) then
                     ! do nothing
                    else
                     print *,"(is_rigid(nmat,im_primary_cell)) invalid"
                     stop
                    endif
                   else
                    print *,"im_primary_cell invalid"
                    stop
                   endif

                  enddo ! sideface
                  enddo ! jtan
                  enddo ! itan

                  if (wtsum.le.zero) then
                   print *,"wtsum invalid"
                   stop
                  endif

                  do im=1,nmat
                   LSmat_node(inode,jnode,im)= &
                    LSmat_node(inode,jnode,im)/wtsum
                   thermal_node(inode,jnode,im)= &
                    thermal_node(inode,jnode,im)/wtsum
                   LSmat(im)=LSmat_node(inode,jnode,im)
                   do im_opp=im+1,nmat
                    call get_iten(im,im_opp,iten,nmat)
                    LSint_node(inode,jnode,iten)= &
                     LSint_node(inode,jnode,iten)/wtsum
                    tension_node(inode,jnode,iten)= &
                      tension_node(inode,jnode,iten)/wtsum
                    LStest=LSint_node(inode,jnode,iten)
                    if (LStest.ge.zero) then
                     LSint_plus(iten)=1
                    else if (LStest.lt.zero) then
                     LSint_minus(iten)=1
                    else
                     print *,"LStest invalid"
                     stop
                    endif
                   enddo ! im_opp=im+1..nmat
                  enddo ! im=1..nmat
                  ! fluid materials should tessellate the domain.
                  ! for fluid materials:
                  ! LS_i^new = (LS_i^old - max_j,i<>j  LS_j^old)/2 
                  call FIX_LS_tessellate(LSmat,LSmat_tess,nmat)
                  do im=1,nmat
                   LSmat_node(inode,jnode,im)=LSmat_tess(im)
                   LSmat(im)=LSmat_tess(im)
                  enddo

                  ! init LSpair
                  im_primary=0
                  do im=1,nmat
                   if (is_rigid(nmat,im).eq.0) then
                    if (im_primary.eq.0) then
                     im_primary=im
                    else if ((im_primary.ge.1).and. &
                             (im_primary.le.nmat)) then
                     if (LSmat(im).gt.LSmat(im_primary)) then
                      im_primary=im
                     endif
                    else
                     print *,"im_primary invalid"
                     stop
                    endif
                   else if (is_rigid(nmat,im).eq.1) then
                    ! do nothing
                   else
                    print *,"is_rigid(nmat,im) invalid"
                    stop
                   endif
                  enddo ! im=1..nmat
                  if ((im_primary.ge.1).and. &
                      (im_primary.le.nmat)) then
        
                   if (LSmat(im_primary).lt.zero) then
                    print *,"LSmat(im_primary).lt.zero"
                    stop
                   endif

                   im_secondary=0
                   do im=1,nmat
                    if (is_rigid(nmat,im).eq.0) then
                     if (im.ne.im_primary) then
                      if (im_secondary.eq.0) then
                       im_secondary=im
                      else if ((im_secondary.ge.1).and. &
                               (im_secondary.le.nmat)) then
                       if (LSmat(im).gt.LSmat(im_secondary)) then
                        im_secondary=im
                       endif
                      else
                       print *,"im_secondary invalid"
                       stop
                      endif 
                     else if (im.eq.im_primary) then
                      ! do nothing
                     else
                      print *,"im or im_primary invalid"
                      stop
                     endif
                    else if (is_rigid(nmat,im).eq.1) then
                     ! do nothing
                    else
                     print *,"is_rigid(nmat,im) invalid"
                     stop
                    endif
                   enddo ! im=1..nmat

                   if ((im_secondary.ge.1).and. &
                       (im_secondary.le.nmat)) then
                    call get_tertiary_material(LSmat,nmat, &
                     im_primary,im_secondary,im_tertiary)
                   else if (im_secondary.eq.0) then
                    im_tertiary=0
                   else
                    print *,"im_secondary invalid"
                    stop
                   endif

                   if ((im_tertiary.ge.0).and.(im_tertiary.le.nmat)) then
                    ! do nothing
                   else
                    print *,"im_tertiary invalid"
                    stop
                   endif

                   do im=1,nmat
                    do im_opp=im+1,nmat
                     call get_iten(im,im_opp,iten,nmat)
                     if (is_rigid(nmat,im).eq.0) then
                      if (is_rigid(nmat,im_opp).eq.0) then

                       if ((im_secondary.ge.1).and. &
                           (im_secondary.le.nmat).and. &
                           (im_primary.ge.1).and. &
                           (im_primary.le.nmat).and. &
                           (im_primary.ne.im_secondary)) then
                      
                        if (im_tertiary.eq.0) then 

                         if (((im.eq.im_primary).and. &
                              (im_opp.eq.im_secondary)).or. &
                             ((im.eq.im_secondary).and. &
                              (im_opp.eq.im_primary))) then
                          LSpair_node(inode,jnode,iten)=99999.0
                         else
                          print *,"im_prim, im_sec, or im_tert invalid"
                          stop
                         endif

                        else if ((im_tertiary.ge.1).and. &
                                 (im_tertiary.le.nmat)) then

                         if (LSmat(im_tertiary).gt.zero) then
                          print *,"LSmat(im_tertiary).gt.zero"
                          stop
                         endif
                         if (LSmat(im_secondary).gt.zero) then
                          print *,"LSmat(im_secondary).gt.zero"
                          stop
                         endif
                         if (LSmat(im_primary).lt.zero) then
                          print *,"LSmat(im_primary).lt.zero"
                          stop
                         endif

                         if (((im.eq.im_primary).and. &
                              (im_opp.eq.im_secondary)).or. &
                             ((im.eq.im_secondary).and. &
                              (im_opp.eq.im_primary))) then
                         
                          LSpair_node(inode,jnode,iten)=-LSmat(im_tertiary)
                         
                         else if ((im.eq.im_secondary).and. &
                                  (im_opp.eq.im_tertiary)) then
                          LSpair_node(inode,jnode,iten)=LSmat(im_secondary)
                         else if ((im_opp.eq.im_secondary).and. &
                                  (im.eq.im_tertiary)) then
                          LSpair_node(inode,jnode,iten)=LSmat(im_secondary)
                         else if ((im.eq.im_primary).and. &
                                  (im_opp.eq.im_tertiary)) then
                          LSpair_node(inode,jnode,iten)=-LSmat(im_tertiary)
                         else if ((im_opp.eq.im_primary).and. &
                                  (im.eq.im_tertiary)) then
                          LSpair_node(inode,jnode,iten)=-LSmat(im_tertiary)
                         else
                          print *,"im or im_opp invalid"
                          stop
                         endif
         
                        else
                         print *,"im_tertiary invalid"
                         stop
                        endif
 
                       else 
                        print *,"im_primary or im_secondary invalid"
                        stop
                       endif

                      else if (is_rigid(nmat,im_opp).eq.1) then
                       LSpair_node(inode,jnode,iten)=-LSmat(im_opp)
                      else
                       print *,"is_rigid(nmat,im_opp) invalid"
                       stop
                      endif
                     else if (is_rigid(nmat,im).eq.1) then
                      LSpair_node(inode,jnode,iten)=LSmat(im)
                     else
                      print *,"is_rigid(nmat,im) invalid"
                      stop
                     endif
                     LStest=LSpair_node(inode,jnode,iten)
                     if (LStest.ge.zero) then
                      LSpair_plus(iten)=1
                     else if (LStest.lt.zero) then
                      LSpair_minus(iten)=1
                     else
                      print *,"LStest invalid"
                      stop
                     endif
                    
                    enddo ! im_opp=im+1 ... nmat
                   enddo ! im=1..nmat
                  else
                   print *,"im_primary invalid"
                   stop
                  endif

                  do dirloc=1,SDIM
                   if (wtsum_deriv(dirloc).le.zero) then
                    print *,"wtsum_deriv invalid"
                    stop
                   endif
                   dx_sten(dirloc)=xhi_sten(dirloc)-xlo_sten(dirloc) 
                   if (dx_sten(dirloc).gt.zero) then
                    do im=1,nmat
                     Nmat_node(inode,jnode,dirloc,im)= &
                      Nmat_node(inode,jnode,dirloc,im)/ &
                      (dx_sten(dirloc)*wtsum_deriv(dirloc))
                     do im_opp=im+1,nmat
                      call get_iten(im,im_opp,iten,nmat)
                      Nint_node(inode,jnode,dirloc,iten)= &
                       Nint_node(inode,jnode,dirloc,iten)/ &
                       (dx_sten(dirloc)*wtsum_deriv(dirloc))
                     enddo ! im_opp=im+1..nmat
                    enddo ! im=1..nmat
                   else
                    print *,"dx_sten(dirloc) invalid"
                    print *,"veldir=",veldir
                    print *,"dir_flux=",dir_flux
                    print *,"dirloc=",dirloc
                    print *,"dx_sten(dirloc)=",dx_sten(dirloc)
                    print *,"xlo_sten(dirloc)=",xlo_sten(dirloc)
                    print *,"xhi_sten(dirloc)=",xhi_sten(dirloc)
                    stop
                   endif
                  enddo ! dirloc=1..sdim
                  RR=xnode(1)
                  do im=1,nmat
                   do dirloc=1,SDIM
                    nhold(dirloc)=Nmat_node(inode,jnode,dirloc,im)
                   enddo
                   call prepare_normal(nhold,RR,mag)
                   do dirloc=1,SDIM
                    Nmat_node(inode,jnode,dirloc,im)=nhold(dirloc)
                   enddo
                   if (mag.eq.zero) then
                    Nmat_node_valid(im)=0
                   else if (mag.gt.zero) then
                    Nmat_node_valid(im)=1
                   else
                    print *,"mag invalid 3"
                    stop
                   endif
 
                   do im_opp=im+1,nmat
                    call get_iten(im,im_opp,iten,nmat)
                    do dirloc=1,SDIM
                     nhold(dirloc)=Nint_node(inode,jnode,dirloc,iten)
                    enddo
                    call prepare_normal(nhold,RR,mag)
                    do dirloc=1,SDIM
                     Nint_node(inode,jnode,dirloc,iten)=nhold(dirloc)
                    enddo
                    if (mag.eq.zero) then
                     Nint_node_valid(iten)=0
                    else if (mag.gt.zero) then
                     Nint_node_valid(iten)=1
                    else
                     print *,"mag invalid 4"
                     stop
                    endif
                   enddo ! im_opp=im+1..nmat
                  enddo ! im=1..nmat

                  if ((im_solid_node.ge.1).and. &
                      (im_solid_node.le.nmat)) then
                   if (is_rigid(nmat,im_solid_node).eq.1) then
                    im=0
                    im_opp=0
                    do im_trial=1,nmat
                     if (is_rigid(nmat,im_trial).eq.0) then
                      if (im.eq.0) then
                       im=im_trial
                      else if ((im.ge.1).and.(im.le.nmat)) then
                       if (LSmat(im_trial).gt.LSmat(im)) then
                        im=im_trial
                       endif
                      else
                       print *,"im invalid25"
                       stop
                      endif
                     else if (is_rigid(nmat,im_trial).eq.1) then
                      ! do nothing
                     else
                      print *,"is_rigid(nmat,im_trial) invalid"
                      stop
                     endif
                    enddo ! im_trial=1..nmat

                    if ((im.ge.1).and.(im.le.nmat)) then
                     do im_opp_trial=1,nmat
                      if ((is_rigid(nmat,im_opp_trial).eq.0).and. &
                          (im_opp_trial.ne.im)) then
                       if (im_opp.eq.0) then
                        im_opp=im_opp_trial
                       else if ((im_opp.ge.1).and.(im_opp.le.nmat)) then
                        if (LSmat(im_opp_trial).gt.LSmat(im_opp)) then
                         im_opp=im_opp_trial
                        endif
                       else
                        print *,"im_opp invalid"
                        stop
                       endif
                      else if ((is_rigid(nmat,im_opp_trial).eq.1).or. &
                               (im_opp_trial.eq.im)) then
                       ! do nothing
                      else
                       print *,"im_opp_trial bust"
                       stop
                      endif
                     enddo ! im_opp_trial=1..nmat

                     if ((im_opp.ge.1).and.(im_opp.le.nmat)) then
                      if (im.lt.im_opp) then
                       ! do nothing
                      else if (im.gt.im_opp) then
                       im_opp_trial=im_opp
                       im_opp=im
                       im=im_opp_trial
                      else
                       print *,"im or im_opp invalid"
                       stop
                      endif
                      call get_iten(im,im_opp,iten,nmat)
                      if ((Nint_node_valid(iten).eq.1).and. &
                          (Nmat_node_valid(im_solid_node).eq.1)) then
                       do iten_local=1,nten
                        user_tension(iten_local)= &
                         tension_node(inode,jnode,iten_local)
                       enddo
                       do im_local=1,nmat
                        thermalmat(im_local)=thermal_node(inode,jnode,im_local)
                       enddo
                       do dirloc=1,SDIM
                        xnode(dirloc)=xnode_array(inode,jnode,dirloc)
                       enddo
                       cos_angle=zero
                       sin_angle=zero
                       iten_13=0
                       iten_23=0
                       if (is_rigid(nmat,im_solid_node).eq.1) then
                        if ((im_solid_node.eq.im).or. &
                            (im_solid_node.eq.im_opp)) then
                         print *,"im_solid_node invalid"
                         stop
                        endif
                        call get_CL_iten(im,im_opp,im_solid_node, &
                         iten_13,iten_23, &
                         user_tension,nten,cos_angle,sin_angle)

                        if (user_tension(iten).eq.zero) then
                         ! do nothing
                        else if (user_tension(iten).gt.zero) then
 
                         ! implement dynamic contact angle algorithm here.
                         ! 1. project nfluid onto the 
                         !    solid (im_solid_node) material

                         ! use_DCA=0 static angle
                         ! use_DCA=1 Jiang
                         ! use_DCA=2 Kistler
                         if ((use_DCA.eq.0).or. &
                             (use_DCA.eq.1).or. &
                             (use_DCA.eq.2)) then

                          dotprod=zero
                          do dirloc=1,SDIM
                           nfluid(dirloc)=Nint_node(inode,jnode,dirloc,iten)
                           nsolid(dirloc)= &
                             Nmat_node(inode,jnode,dirloc,im_solid_node)
                           dotprod=dotprod+nfluid(dirloc)*nsolid(dirloc) 
                          enddo
                          ! nproject=(I-nsolid nsolid^T)nfluid
                          do dirloc=1,SDIM
                           nproject(dirloc)=nfluid(dirloc)- &
                             dotprod*nsolid(dirloc)
                          enddo
                          RR=one
                          call prepare_normal(nproject,RR,mag)
          
                          if (mag.gt.zero) then
                           ! find u dot nproject
                           udotn=zero
                           do dirloc=1,SDIM
                            udotn=udotn+ &
                             vel_for_CL_model(dirloc)*nproject(dirloc)
                           enddo
                           if (fort_denconst(im).ge.fort_denconst(im_opp)) then
                            imCL=im
                            denCL=fort_denconst(im)
                            T_CL=thermalmat(im)
                           else
                            imCL=im_opp
                            denCL=fort_denconst(im_opp)
                            T_CL=thermalmat(im_opp)
                           endif 
                           liquid_viscosity=visc_coef* &
                            get_user_viscconst(imCL,denCL,T_CL) 
        
                           ! modify cos_angle 
                           !  (initialized above as static angle)
                           ! use_DCA=0 static angle
                           ! use_DCA=1 Jiang
                           ! use_DCA=2 Kistler
                           call DCA_select_model(nproject,udotn,cos_angle, &
                            liquid_viscosity,user_tension(iten),cos_angle, &
                            use_DCA)
       
                           if (cos_angle.gt.one) then 
                            cos_angle=one
                           else if (cos_angle.lt.-one) then
                            cos_angle=-one
                           endif
       
                          else if (mag.eq.zero) then
                           ! do nothing (nproject has mag=0)
                          else
                           print *,"mag cannot be negative"
                           stop
                          endif 

                          ! fort_ZEYU_DCA_SELECT>=1
                         else if ((use_DCA.ge.101).and. & 
                                  (use_DCA.le.108)) then
                          if (use_DCA.eq.101) then
                           ! do nothing
                          else if ((use_DCA.ge.102).and.(use_DCA.le.108)) then
                           print *,"use_DCA>=102 and <=107 not supported"
                           print *,"for conservative surface tension alg."
                           stop
                          else
                           print *,"use_DCA bust"
                           stop
                          endif
                         else if (use_DCA.eq.-1) then
                          ! do nothing
                         else
                          print *,"use_DCA invalid"
                          stop
                         endif 

                         do dirloc=1,SDIM
                          nfluid(dirloc)=Nint_node(inode,jnode,dirloc,iten)
                          nsolid(dirloc)= &
                           -Nmat_node(inode,jnode,dirloc,im_solid_node)
                         enddo
                         call ghostnormal(nfluid,nsolid,cos_angle,nghost,nperp)
                         do dirloc=1,SDIM
                          Nint_node(inode,jnode,dirloc,iten)=nghost(dirloc)
                         enddo
                        else
                         print *,"user_tension(iten) cannot be negative"
                         stop
                        endif
                       else 
                        print *,"is_rigid(nmat,im_solid_node) invalid"
                        stop
                       endif
                      else if ((Nint_node_valid(iten).eq.0).or. &
                               (Nmat_node_valid(im_solid_node).eq.0)) then
                       ! do nothing
                      else
                       print *,"Nint_node_valid or Nmat_node_valid bust"
                       stop
                      endif
                     else if (im_opp.eq.0) then
                      ! do nothing
                     else
                      print *,"im_opp invalid"
                      stop
                     endif
                    else
                     print *,"im invalid (must have >=1 fluids)"
                     stop
                    endif 
                   else
                    print *,"is_rigid(nmat,im_solid_node) invalid"
                    stop
                   endif
                  else if (im_solid_node.eq.0) then
                   ! do nothing
                  else
                   print *,"im_solid_node invalid"
                   stop
                  endif

                 enddo ! jnode=0..1
                 enddo ! inode=0..1

                 if (dir_flux.eq.0) then
                  dx_node(1)=xnode_array(1,0,2)- &
                             xnode_array(0,0,2) ! dy
                  dx_node(2)=xnode_array(0,1,SDIM)- &
                             xnode_array(0,0,SDIM) ! dz
                  RR=xnode_array(0,0,1)
                  if (SDIM.eq.2) then
                   areaface_MAC=dx_node(1) ! dy
                   if (levelrz.eq.0) then
                    ! do nothing
                   else if (levelrz.eq.1) then
                    if (RR.ge.-VOFTOL*dx(1)) then
                     areaface_MAC=areaface_MAC*two*Pi*abs(RR) !2 pi r dz
                    else
                     print *,"RR invalid RR=",RR
                     stop
                    endif
                   else if (levelrz.eq.3) then
                    if (RR.le.zero) then
                     print *,"RR invalid"
                     stop
                    endif
                    areaface_MAC=areaface_MAC*RR ! r dtheta
                   else
                    print *,"levelrz invalid"
                    stop
                   endif
                  else if (SDIM.eq.3) then
                   areaface_MAC=dx_node(1)*dx_node(2) ! dy dz
                   if (levelrz.eq.0) then
                    ! do nothing
                   else if (levelrz.eq.3) then
                    if (RR.le.zero) then
                     print *,"RR invalid"
                     stop
                    endif
                    areaface_MAC=areaface_MAC*RR ! r dtheta dz
                   else
                    print *,"levelrz invalid"
                    stop
                   endif
                  else
                   print *,"dimension bust"
                   stop
                  endif
                 else if (dir_flux.eq.1) then
                  dx_node(1)=xnode_array(1,0,1)- &
                    xnode_array(0,0,1) ! dx
                  dx_node(2)=xnode_array(0,1,SDIM)- &
                    xnode_array(0,0,SDIM) ! dz
                  RR=half*(xnode_array(0,0,1)+xnode_array(1,0,1))
                  
                  if (SDIM.eq.2) then
                   areaface_MAC=dx_node(1) ! dx
                   if (levelrz.eq.0) then
                    ! do nothing
                   else if (levelrz.eq.1) then
                    if (RR.gt.zero) then
                     areaface_MAC=areaface_MAC*two*Pi*RR !2 pi r dr
                    else
                     print *,"RR invalid"
                     stop
                    endif
                   else if (levelrz.eq.3) then
                    ! do nothing (dr)
                   else
                    print *,"levelrz invalid"
                    stop
                   endif
                  else if (SDIM.eq.3) then
                   areaface_MAC=dx_node(1)*dx_node(2) ! dx dz
                   if (levelrz.eq.0) then
                    ! do nothing
                   else if (levelrz.eq.3) then
                    ! do nothing (dr dz)
                   else
                    print *,"levelrz invalid"
                    stop
                   endif
                  else
                   print *,"dimension bust"
                   stop
                  endif
                 else if ((dir_flux.eq.2).and.(SDIM.eq.3)) then
                  dx_node(1)=xnode_array(1,0,1)-xnode_array(0,0,1) ! dx
                  dx_node(2)=xnode_array(0,1,2)-xnode_array(0,0,2) ! dy
                  RR=half*(xnode_array(0,0,1)+xnode_array(1,0,1))
                  
                  if (SDIM.eq.3) then
                   areaface_MAC=dx_node(1)*dx_node(2) ! dx dy
                   if (levelrz.eq.0) then
                    ! do nothing
                   else if (levelrz.eq.3) then
                    if (RR.gt.zero) then
                     areaface_MAC=areaface_MAC*RR ! r dr dtheta
                    else
                     print *,"RR invalid"
                     stop
                    endif
                   else
                    print *,"levelrz invalid"
                    stop
                   endif
                  else
                   print *,"dimension bust"
                   stop
                  endif
                 else
                  print *,"dir_flux invalid"
                  stop
                 endif

                 if (dir_flux.eq.veldir) then
                  areaface_slice_check=areaface_slice_check+ &
                    half*areaface_MAC
                 endif
                
                 do im=1,nmat
                  do im_opp=im+1,nmat
                   if ((is_rigid(nmat,im).eq.0).and. &
                       (is_rigid(nmat,im_opp).eq.0)) then
                    call get_iten(im,im_opp,iten,nmat)
                    if ((LSpair_plus(iten).eq.1).and. &
                        (LSint_plus(iten).eq.1).and. &
                        (LSint_minus(iten).eq.1).and. &
                        (Nint_node_valid(iten).eq.1)) then
                     ncrossing=0
                     do nside=1,4
                      if (nside.eq.1) then
                       inode1=0
                       jnode1=0
                       inode2=1
                       jnode2=0
                      else if (nside.eq.2) then
                       inode1=1
                       jnode1=0
                       inode2=1
                       jnode2=1
                      else if (nside.eq.3) then
                       inode1=1
                       jnode1=1
                       inode2=0
                       jnode2=1
                      else if (nside.eq.4) then
                       inode1=0
                       jnode1=1
                       inode2=0
                       jnode2=0
                      else
                       print *,"nside invalid"
                       stop
                      endif
                      call add_crossing(ncrossing,xcrossing,iten, &
                       LSint_node,xnode_array, &
                       inode1,jnode1,inode2,jnode2,nten)
                     enddo ! nside=1..4

                     if (ncrossing.eq.0) then
                      maxdist=zero
                     else if (SDIM.eq.2) then
                      if (ncrossing.ne.2) then
                       print *,"(ncrossing.ne.2)"
                       stop
                      endif
                      maxdist=one
                      cross1max=1
                      cross2max=2
                     else if (SDIM.eq.3) then
                      maxdist=zero
                      cross1max=0
                      cross2max=0
                      do icross1=1,ncrossing
                       do icross2=icross1+1,ncrossing
                        mag=zero
                        do dirloc=1,SDIM 
                         mag=mag+ &
                          (xcrossing(icross1,dirloc)- &
                           xcrossing(icross2,dirloc))**2
                        enddo
                        mag=sqrt(mag)
                        if (maxdist.eq.zero) then
                         maxdist=mag
                         cross1max=icross1
                         cross2max=icross2
                        else if (maxdist.gt.zero) then
                         if (mag.gt.maxdist) then
                          maxdist=mag
                          cross1max=icross1
                          cross2max=icross2
                         else if ((mag.le.maxdist).and.(mag.ge.zero)) then
                          ! do nothing
                         else
                          print *,"mag invalid 5"
                          stop
                         endif
                        else
                         print *,"maxdist invalid"
                         stop
                        endif
                       enddo ! icross2
                      enddo ! icross1
                     else
                      print *,"dimension bust"
                      stop
                     endif

                     if (maxdist.gt.LS_CURV_TOL*dxmax) then
                      do dirloc=1,SDIM
                       x1(dirloc)=xcrossing(cross1max,dirloc)
                       x2(dirloc)=xcrossing(cross2max,dirloc)
                      enddo
                      call interp_face(x1,xnode_array,LSpair_node, &
                       iten,nten,LS1,dir_flux)
                      call interp_face(x2,xnode_array,LSpair_node, &
                       iten,nten,LS2,dir_flux)
                      if ((LS1.lt.zero).and.(LS2.lt.zero)) then
                       ! do nothing
                      else if ((LS1.ge.zero).or.(LS2.ge.zero)) then
                       if ((LS1.ge.zero).and.(LS2.ge.zero)) then
                        ! do nothing
                       else if ((LS1.ge.zero).and.(LS2.lt.zero)) then
                        frac=LS1/(LS1-LS2)
                        do dirloc=1,SDIM
                         x2(dirloc)=(one-frac)*x1(dirloc)+frac*x2(dirloc)
                        enddo
                       else if ((LS1.lt.zero).and.(LS2.ge.zero)) then
                        frac=LS1/(LS1-LS2)
                        do dirloc=1,SDIM
                         x1(dirloc)=(one-frac)*x1(dirloc)+frac*x2(dirloc)
                        enddo
                       else
                        print *,"LS1 or LS2 invalid"
                        stop
                       endif

                       mag=zero
                       do dirloc=1,SDIM
                        xmid(dirloc)=half*(x1(dirloc)+x2(dirloc))
                        abs_x1=x1(dirloc)
                        abs_x2=x2(dirloc)
                        if (SDIM.eq.2) then
                         ! do nothing
                        else if (SDIM.eq.3) then
                         if (levelrz.eq.0) then
                          ! do nothing
                         else if (levelrz.eq.3) then
                          if ((dirloc.eq.1).or.(dirloc.eq.3)) then
                           ! do nothing
                          else if (dirloc.eq.2) then 
                           if ((x1(1).gt.zero).and.(x2(1).gt.zero)) then
                            abs_x1=abs_x1*x1(1)
                            abs_x2=abs_x2*x2(1)
                           else
                            print *,"x1(1),x2(1) must be positive"
                            stop
                           endif
                          else
                           print *,"dirloc invalid"
                           stop
                          endif
                         else
                          print *,"levelrz invalid"
                          stop
                         endif
                        else
                         print *,"dimension bust"
                         stop
                        endif
 
                        mag=mag+(abs_x2-abs_x1)**2
                       enddo
                       mag=sqrt(mag)
                       call interp_face(xmid,xnode_array, &
                        tension_node,iten,nten,ten12,dir_flux)
                       call interp_face_normal(xmid,xnode_array, &
                        Nint_node,iten,nten,normal12,dir_flux)
                       RR=one
                       call prepare_normal(normal12,RR,mag12)
     
                       ! non-conservative surface tension force (2 mat):
                       !   (sigma constant)
                       !   sigma kappa grad H(phi)
                       ! units: sigma * (1/L) * (1/L) = sigma/L^2
                       !
                       ! conservative surface tension force (2 mat):
                       !   div[ (I-nn^T)sigma delta(phi) ] = div S
                       ! units: (1/L) * sigma * (1/L) = sigma/L^2
                       !
                       ! cylindrical coordinates: (r,theta,z)
                       ! (div S)rhat=(r S11)_r/r + S12,theta/r + S13,z - S22/r
                       ! (div S)that=(r S21)_r/r + S22,theta/r + S23,z + S12/r
                       ! (div S)zhat=(r S31)_r/r + S32,theta/r + S33,z 
                       !
                       ! div S = (1/|Omega|)sum_face int_Gamma_face S dot nface
                       !
                       ! delta=H'(phi)=div(Hn)  (phi is a distance function)
                       ! Let m be n projected to a face:
                       ! m=(I-nf nf^T)n   ( m dot nf=0 )
                       ! delta=div_face (Hm)/(n dot m)
                       ! int_Gamma_face(I-nn^T)delta dot nface=
                       ! int_Gamma_face(I-nn^T)dot nface div_face(Hm)/
                       !               (n dot m)=
                       ! |Gamma intersect Gamma_face|(I-nn^T)dot nface/
                       !                             (n dot m)
                       ! n dot m=n dot (n-(n dot nf)nf)=1-(n dot nf)^2
                       ! (I-nn^T)dot nf dot (I-nn^T)dot nf =
                       ! (nf-(n dot nf)n) dot (nf-(n dot nf)n)=
                       ! 1-2(n dot nf)^2+(n dot nf)^2)=1-(n dot nf)^2=n dot m
                       ! int_Gamma_face S dot nface=
                       !  |Gamma intersect Gamma_face| *
                       !  normalize( (I-nn^T)dot nface )
                       ! 

                        ! dotprod=n dot nface
                       dotprod=zero
                       do dirloc=1,SDIM
                        if (dirloc.eq.dir_flux+1) then
                         nface(dirloc)=one
                        else
                         nface(dirloc)=zero
                        endif
                        dotprod=dotprod+nface(dirloc)*normal12(dirloc)
                       enddo
                        
                       if (SDIM.eq.2) then
                        if (levelrz.eq.0) then
                         if (areaface_MAC.le.zero) then
                          print *,"areaface_MAC invalid"
                          stop
                         endif
                         delta_face=one/areaface_MAC
                        else if (levelrz.eq.1) then
                         if ((i.eq.0).and. &
                             (dir_flux.eq.0).and. &
                             (sideflux.eq.0).and. &
                             (veldir.ge.1).and. &
                             (veldir.lt.SDIM)) then
                          if (abs(areaface_MAC).lt.VOFTOL*dx(1)) then
                           delta_face=zero
                          else
                           print *,"areaface_MAC failed sanity check"
                           stop
                          endif
                         else if ((i.gt.0).or. &
                                  (dir_flux.eq.1).or. &
                                  (sideflux.eq.1).or. &
                                  (veldir.eq.0)) then 
                          if (areaface_MAC.le.zero) then
                           print *,"areaface_MAC invalid"
                           stop
                          endif
                          delta_face=two*Pi*abs(xmid(1))/areaface_MAC
                         else
                          print *,"i, dir_flux or sideflux invalid"
                          stop
                         endif
                        else if (levelrz.eq.3) then
                         if (areaface_MAC.gt.zero) then
                          delta_face=one/areaface_MAC
                         else
                          print *,"areaface_MAC invalid"
                          stop
                         endif
                        else
                         print *,"levelrz invalid"
                         stop
                        endif
                       else if (SDIM.eq.3) then
                         ! Gamma=interface
                         ! Gamma_face=face of cell
                         ! mag=|Gamma intersect Gamma_face|
                         ! areaface_MAC=|Gamma_face|
                        if (areaface_MAC.gt.zero) then
                         if (levelrz.eq.0) then
                          delta_face=mag/areaface_MAC
                         else if (levelrz.eq.3) then
                          delta_face=mag/areaface_MAC
                         else
                          print *,"levelrz invalid"
                          stop
                         endif
                        else
                         print *,"areaface_MAC invalid"
                         stop
                        endif
                       else 
                        print *,"dimension bust"
                        stop
                       endif
                       ! t_interface is tangent to the interface and points
                       ! out of the face.
                       ! t_interface=(I-nn^T)dot nface
                       ! (projection of nface onto the plane whose 
                       !  normal is n)
                       do dirloc=1,SDIM
                        t_interface(dirloc)=nface(dirloc)- &
                         dotprod*normal12(dirloc)
                       enddo
                        ! normalize t_interface
                       RR=one
                       call prepare_normal(t_interface,RR,mag12)
                    
                       if ((delta_face.gt.zero).and. &
                           (ten12.gt.zero).and. &
                           (hoop_force_iten.gt.0).and. &
                           (iten.eq.hoop_force_iten)) then
                        hoop_force_ok=1
                       else if ((delta_face.eq.zero).or. &
                                (ten12.eq.zero).or. &
                                (hoop_force_iten.eq.0).or. &
                                (iten.ne.hoop_force_iten)) then
                        ! do nothing
                       else
                        print *,"delta_face,ten12,hoop_force_iten or iten bad"
                        stop
                       endif
                        ! areaface_MAC=|Gamma_face|
                        ! delta_face=|Gamma intersect Gamma_face|/|Gamma_face| 
                       do dirloc=1,SDIM
                        force_mag=areaface_MAC*delta_face*ten12* &
                          t_interface(dirloc)/vol_sten

                         ! approximation to div (sigma(I-nn)delta)=
                         ! sigma grad delta - sigma kappa grad h 
                         ! if sigma is constant. 
                        if (sideflux.eq.0) then
                         surface_tension_force(dirloc)= &
                          surface_tension_force(dirloc)-force_mag
                        else if (sideflux.eq.1) then
                         surface_tension_force(dirloc)= &
                          surface_tension_force(dirloc)+force_mag
                        else
                         print *,"sideflux invalid"
                         stop
                        endif
                       enddo ! dirloc=1..sdim

                       if (debug_curv.eq.1) then
                        if ((i.eq.3).and. &
                            (j.eq.6).and. &
                            (k.eq.10)) then
                          ! phi_expect=y+radblob(x-blob)-yblob  or
                          ! phi_expect=z+radblob(x-xblob)+ 
                          !              radblob2(y-yblob)-zblob
                         nexpect(1)=radblob
                         nexpect(SDIM)=one
                         if (SDIM.eq.2) then
                          ! do nothing
                         else if (SDIM.eq.3) then
                          nexpect(2)=radblob2
                         else
                          print *,"dimension bust"
                          stop
                         endif
                         RR=one
                         call prepare_normal(nexpect,RR,mag12)
                         dotprod=zero
                         do dirloc=1,SDIM
                          dotprod=dotprod+nface(dirloc)*nexpect(dirloc)
                         enddo
                         do dirloc=1,SDIM
                          texpect(dirloc)=nface(dirloc)-dotprod*nexpect(dirloc)
                         enddo
                         RR=one
                         call prepare_normal(texpect,RR,mag12)
                         print *,"debug_curv: i,j,k,dir_flux ",i,j,k,dir_flux
                         print *,"sideflux=",sideflux
                         print *,"veldir: ",veldir
                         print *,"delta_face,areaface_MAC,ten12 ", &
                          delta_face,areaface_MAC,ten12
                         do dirloc=1,SDIM
                          print *,"dirloc,nexpect,texpect,t_interface ", &
                           dirloc,nexpect(dirloc),texpect(dirloc), &
                           t_interface(dirloc) 
                         enddo
                         do dirloc=1,SDIM
                          print *,"dirloc,fablo,fabhi ",dirloc, &
                            fablo(dirloc),fabhi(dirloc)
                         enddo
                         do dirloc=1,SDIM
                          print *,"dirloc,xmid ",dirloc,xmid(dirloc)
                         enddo
                         do inode=0,1
                         do jnode=0,1
                         do dirloc=1,SDIM
                          print *,"inode,jnode,dirloc,xnode_array ", &
                           inode,jnode,dirloc,xnode_array(inode,jnode,dirloc)
                         enddo
                         enddo
                         enddo
                        endif
                       else if (debug_curv.eq.0) then
                        ! do nothing
                       else
                        print *,"debug_curv invalid"
                        stop
                       endif
                      else
                       print *,"LS1 or LS2 invalid"
                       stop
                      endif
                     else if (maxdist.ge.zero) then
                      ! do nothing
                     else
                      print *,"maxdist invalid"
                      stop
                     endif
                    else if ((LSpair_plus(iten).eq.0).or. &
                             (LSint_plus(iten).eq.0).or. &
                             (LSint_minus(iten).eq.0).or. &
                             (Nint_node_valid(iten).eq.0)) then
                     ! do nothing
                    else
                     print *,"LSpair, LSint, or Nint_node_valid invalid"
                     stop
                    endif
                   else if ((is_rigid(nmat,im).eq.1).or. &
                            (is_rigid(nmat,im_opp).eq.1)) then
                    ! do nothing
                   else
                    print *,"is_rigid invalid"
                    stop
                   endif
                  enddo ! im_opp=im+1..nmat
                 enddo ! im=1..nmat
                enddo ! sideflux=0..1
                enddo ! dir_flux=0..sdim-1

                if (abs(areaface_slice_check-areaface_slice).gt. &
                    VOFTOL*dx(veldir+1)) then
                 print *,"areaface_slice_check or areaface_slice invalid"
                 print *,"areaface_slice_check=",areaface_slice_check
                 print *,"areaface_slice=",areaface_slice
                 stop
                endif

               else
                print *,"FVM_surface_tension invalid"
                stop
               endif

               ! surface_tension_force is an approximation to:
               !   -sigma kappa grad H (if sigma=constant)
               ! sideface=1 is right half of cell that is to left of face.
               ! sideface=2 is left half of cell that is to right of face.
               do dirloc=1,SDIM
                force_hold=mgoni_force(dirloc)+surface_tension_force(dirloc)
                if (hoop_force_ok.eq.1) then
                 if (hoop_force_iten.gt.0) then
                  force_hold=force_hold+hoop_force(dirloc)
                 else
                  print *,"hoop_force_iten invalid"
                  stop
                 endif
                else if (hoop_force_ok.eq.0) then
                 ! do nothing
                else
                 print *,"hoop_force_ok invalid"
                 stop
                endif

                if ((local_face(facecut_index+1).ge.zero).and. &
                    (local_face(facecut_index+1).le.half)) then
                 forcelocal=zero
                else if ((local_face(facecut_index+1).ge.half).and. &
                         (local_face(facecut_index+1).le.one)) then
                 forcelocal= &
                  dt*force_hold* &
                  local_face(facecut_index+1)* &
                  local_face(faceden_index+1)
                else
                 print *,"local_face(facecut_index+1) invalid"
                 stop
                endif

                if (veldir.eq.0) then
                 xflux(D_DECL(i,j,k),dirloc)=forcelocal
                 if (dirloc.eq.veldir+1) then
                  umac(D_DECL(i,j,k))=umac(D_DECL(i,j,k))+forcelocal
                 endif
                else if (veldir.eq.1) then
                 yflux(D_DECL(i,j,k),dirloc)=forcelocal
                 if (dirloc.eq.veldir+1) then
                  vmac(D_DECL(i,j,k))=vmac(D_DECL(i,j,k))+forcelocal
                 endif
                else if ((veldir.eq.2).and.(SDIM.eq.3)) then
                 zflux(D_DECL(i,j,k),dirloc)=forcelocal
                 if (dirloc.eq.veldir+1) then
                  wmac(D_DECL(i,j,k))=wmac(D_DECL(i,j,k))+forcelocal
                 endif
                else
                 print *,"veldir invalid"
                 stop
                endif
                 ! denface=den_local(1)+den_local(2)
                do sideface=1,2
                 forcelocal= &
                  dt*force_hold* &
                  hx*den_local(sideface)/denface

                 if (veldir.eq.0) then
                  xp(D_DECL(i,j,k),SDIM*(sideface-1)+dirloc)=forcelocal
                 else if (veldir.eq.1) then
                  yp(D_DECL(i,j,k),SDIM*(sideface-1)+dirloc)=forcelocal
                 else if ((veldir.eq.2).and.(SDIM.eq.3)) then
                  zp(D_DECL(i,j,k),SDIM*(sideface-1)+dirloc)=forcelocal
                 else
                  print *,"veldir invalid"
                  stop
                 endif
                enddo ! sideface=1,2
               enddo ! dirloc=1..sdim

              else if (is_solid_face.eq.1) then
               ! do nothing
              else
               print *,"is_solid_face invalid"
               stop
              endif

             else
              print *,"AFACE or AFACE_ICE invalid"
              stop
             endif

            else
             print *,"at_rz_center invalid"
             stop
            endif

           else
            print *,"is_rigid invalid"
            stop
           endif
       
          else if ((maskleft.eq.0).or.(maskright.eq.0)) then
           ! do nothing
          else
           print *,"maskleft or maskright invalid"
           stop
          endif

         enddo ! k
         enddo ! j
         enddo ! i

        enddo ! veldir=0..sdim-1

       else if (conservative_tension_force.eq.0) then
        ! do nothing
       else
        print *,"conservative_tension_force invalid"
        stop
       endif

      else if (isweep.eq.1) then
 
       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)

        call gridsten_level(xsten,i,j,k,level,nhalf)

        do dirloc=1,SDIM
         surface_tension_force(dirloc)=zero
        enddo

        do im=1,nmat
         LScen(im)=ls(D_DECL(i,j,k),im)
        enddo
        call get_primary_material(LScen,nmat,im)

        if (is_rigid(nmat,im).eq.0) then ! fluid region

         if (conservative_tension_force.eq.1) then

          do dirloc=0,SDIM-1 
           ii=0
           jj=0
           kk=0
           if (dirloc.eq.0) then
            ii=1
           else if (dirloc.eq.1) then
            jj=1
           else if ((dirloc.eq.2).and.(SDIM.eq.3)) then
            kk=1
           else
            print *,"dirloc out of range in MARANGONIFORCE"
            stop
           endif

           hx=xsten(1,dirloc+1)-xsten(-1,dirloc+1)
           RR=one
           if ((levelrz.eq.0).or. &
               (levelrz.eq.1)) then
            RR=one
           else if (levelrz.eq.3) then
            if (dirloc.eq.1) then ! theta direction
             RR=xsten(0,1)
            else if ((dirloc.eq.0).or.(dirloc.eq.SDIM-1)) then
             RR=one
            else
             print *,"dirloc invalid MARANGONIFORCE"
             stop
            endif
           else
            print *,"levelrz invalid MARANGONIFORCE 3"
            stop
           endif
           hx=hx*RR
           if (hx.le.zero) then
            print *,"hx invalid"
            stop
           endif

           ! sideface=1 left half of cell, sideface=2 right half of cell
           do sideface=1,2

            if (sideface.eq.1) then
             iface=i
             jface=j
             kface=k
            else if (sideface.eq.2) then
             iface=i+ii
             jface=j+jj
             kface=k+kk
            else
             print *,"sideface invalid"
             stop
            endif

            if (sideface.eq.1) then ! left half of cell
             sidecomp=SDIM+dirloc+1
            else if (sideface.eq.2) then ! right half of cell
             sidecomp=dirloc+1
            else
             print *,"sidecomp invalid"
             stop
            endif
 
            if (dirloc.eq.0) then 
             do im=1,ncphys
              local_face(im)=xface(D_DECL(iface,jface,kface),im)
             enddo
             fside=xp(D_DECL(iface,jface,kface),sidecomp)
            else if (dirloc.eq.1) then
             do im=1,ncphys
              local_face(im)=yface(D_DECL(iface,jface,kface),im)
             enddo
             fside=yp(D_DECL(iface,jface,kface),sidecomp)
            else if ((dirloc.eq.2).and.(SDIM.eq.3)) then
             do im=1,ncphys
              local_face(im)=zface(D_DECL(iface,jface,kface),im)
             enddo
             fside=zp(D_DECL(iface,jface,kface),sidecomp)
            else
             print *,"dirloc invalid MARANGONIFORCE 5"
             stop
            endif

            ffaceten(sideface)=fside

            mass_local(sideface)=zero
            do im=1,nmat 
             if (sideface.eq.1) then  ! left half of cell
              sidecomp=massface_index+2*(im-1)+2
             else if (sideface.eq.2) then ! right half of cell
              sidecomp=massface_index+2*(im-1)+1
             else
              print *,"sideface invalid"
              stop
             endif
             mass_local(sideface)=mass_local(sideface)+local_face(sidecomp) 
            enddo ! im=1..nmat

           enddo ! sideface=1..2

           masscell=mass_local(1)+mass_local(2)

           if ((mass_local(1).gt.zero).and. &
               (mass_local(2).gt.zero).and. &
               (masscell.gt.zero)) then

            volpos=vol(D_DECL(i,j,k))
            if (volpos.gt.zero) then
             dencell=masscell/volpos
             if (dencell.gt.zero) then
              surface_tension_force(dirloc+1)= &
                (ffaceten(2)+ffaceten(1))/(dencell*hx)
             else
              print *,"dencell invalid"
              stop
             endif
            else
             print *,"volpos invalid"
             stop
            endif
           else
            print *,"mass_local or masscell invalid"
            stop
           endif

          enddo ! dirloc=0..sdim-1
 
         else if (conservative_tension_force.eq.0) then

          ! curv: nten * (5+SDIM)
          !  curv_cellHT,curv_cellFD, pforce_cell,
          !  marangoni force(sdim),
          !  dir * side, im3
          do iten=1,nten

           do dirloc=1,SDIM
            iforce=(iten-1)*(5+SDIM)+3+dirloc
            surface_tension_force(dirloc)= &
             surface_tension_force(dirloc)+ &
             curv(D_DECL(i,j,k),iforce)*dt*rhoinverse(D_DECL(i,j,k),1)
           enddo  ! dirloc

          enddo !iten
  
         else
          print *,"conservative_tension_force invalid"
          stop
         endif

        else if (is_rigid(nmat,im).eq.1) then 
         ! do nothing
        else
         print *,"is_rigid(nmat,im) invalid"
         stop
        endif 
      
        do dirloc=1,SDIM
         velnew(D_DECL(i,j,k),dirloc)= &
          velnew(D_DECL(i,j,k),dirloc)+surface_tension_force(dirloc)

         ! sanity check is sensitive to FACETOL_DVOL used in LEVELSTRIP
         ! see PROBCOMMON.F90 for other possible sensitivities.
         if (curv_sanity_check.eq.1) then
          velloc=velnew(D_DECL(i,j,k),dirloc)
          if (abs(velloc).gt.1.0e-10) then
           print *,"curv_sanity_check: i,j,k,dirloc,velloc ", &
            i,j,k,dirloc,velloc
           do im=1,nmat
            print *,"im,F= ",im,vof(D_DECL(i,j,k),(im-1)*ngeom_recon+1)
           enddo
          endif
         else if (curv_sanity_check.eq.0) then
          ! do nothing
         else
          print *,"curv_sanity_check invalid"
          stop
         endif
        enddo ! dirloc

       enddo
       enddo
       enddo ! i,j,k

      else
       print *,"isweep invalid"
       stop
      endif

      return
      end subroutine FORT_MARANGONIFORCE

       ! MEHDI VAHAB HEAT SOURCE
       ! T^new=T^* + dt * Q/(rho cv)
       ! Q units: J/(m^3 s)
       ! called from: make_heat_source
       ! make_heat_source is called from veldiffuseALL
      subroutine FORT_HEATSOURCE( &
       nstate, &
       nmat, &
       nden, &
       xlo,dx,  &
       temperature_source, &
       temperature_source_cen, &
       temperature_source_rad, &
       rhoinverse, &
       DIMS(rhoinverse), &
       DeDTinverse, &
       DIMS(DeDTinverse), &
       Tnew,DIMS(Tnew), &
       lsfab,DIMS(lsfab), &
       recon,DIMS(recon), &
       vol,DIMS(vol), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       finest_level, &
       dt,time)
      use probf90_module
      use global_utility_module
      use MOF_routines_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: nstate
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nden
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T, intent(in) :: DIMDEC(rhoinverse)
      INTEGER_T, intent(in) :: DIMDEC(DeDTinverse)
      INTEGER_T, intent(in) :: DIMDEC(Tnew)
      INTEGER_T, intent(in) :: DIMDEC(lsfab)
      INTEGER_T, intent(in) :: DIMDEC(recon)
      INTEGER_T, intent(in) :: DIMDEC(vol)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: temperature_source
      REAL_T, intent(in) :: temperature_source_cen(SDIM)
      REAL_T, intent(in) :: temperature_source_rad(SDIM)
      REAL_T, intent(in) :: rhoinverse(DIMV(rhoinverse),nmat+1)
      REAL_T, intent(in) :: DeDTinverse(DIMV(DeDTinverse),nmat+1) ! 1/(rho cv)
      REAL_T, intent(inout) :: Tnew(DIMV(Tnew),nden)
      REAL_T, intent(in) :: lsfab(DIMV(lsfab),nmat)
      REAL_T, intent(in) :: recon(DIMV(recon),nmat*ngeom_recon)
      REAL_T, intent(in) :: vol(DIMV(vol))
      REAL_T, intent(in) :: dt,time

      REAL_T xsten(-3:3,SDIM)
      REAL_T xsten_cell(SDIM)
      INTEGER_T nhalf
      REAL_T LS(nmat)
      REAL_T VFRAC(nmat)
      REAL_T heat_source_local(nmat)
      REAL_T T_local(nmat)
      REAL_T den_local(nmat)
      INTEGER_T i,j,k
      INTEGER_T im
      INTEGER_T vofcomp,dencomp
      INTEGER_T dirloc
      INTEGER_T imattype
      REAL_T heat_source_total,vfrac_total
      REAL_T DeDT_local(nmat)
      INTEGER_T ispec
      REAL_T massfrac_parm(num_species_var+1)

      nhalf=3

      if (bfact.lt.1) then
       print *,"bfact invalid53"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (nstate.ne.(SDIM+1)+ &
          nmat*(num_state_material+ngeom_raw)+1) then
       print *,"nstate invalid"
       stop
      endif
      if (nden.ne.nmat*num_state_material) then
       print *,"nden invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid heat source"
       stop
      endif
      if (dt.le.zero) then
       print *,"dt invalid"
       stop
      endif
      if (time.lt.zero) then
       print *,"time invalid"
       stop
      endif
      if (temperature_source.lt.zero) then
       print *,"temperature_source invalid"
       stop
      endif

      call checkbound(fablo,fabhi, &
       DIMS(rhoinverse), &
       1,-1,7)
      call checkbound(fablo,fabhi, &
       DIMS(DeDTinverse), &
       1,-1,7)
      call checkbound(fablo,fabhi,DIMS(Tnew),1,-1,7)
      call checkbound(fablo,fabhi,DIMS(lsfab),1,-1,7)
      call checkbound(fablo,fabhi,DIMS(recon),1,-1,7)
      call checkbound(fablo,fabhi,DIMS(vol),1,-1,7)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
       call gridsten_level(xsten,i,j,k,level,nhalf)

       do im=1,nmat
        LS(im)=lsfab(D_DECL(i,j,k),im)
        vofcomp=(im-1)*ngeom_recon+1
        VFRAC(im)=recon(D_DECL(i,j,k),vofcomp) 
        dencomp=(im-1)*num_state_material+1
        den_local(im)=Tnew(D_DECL(i,j,k),dencomp)
        T_local(im)=Tnew(D_DECL(i,j,k),dencomp+1)
        if ((VFRAC(im).ge.-VOFTOL).and. &
            (VFRAC(im).le.one+VOFTOL)) then
         ! do nothing
        else
         print *,"VFRAC invalid"
         stop
        endif 
        if (den_local(im).gt.zero) then
         ! do nothing
        else
         print *,"den_local must be positive"
         stop
        endif
        if (T_local(im).gt.zero) then
         ! do nothing
        else
         print *,"T_local must be positive"
         stop
        endif
        call init_massfrac_parm(den_local(im),massfrac_parm,im)
        do ispec=1,num_species_var
         massfrac_parm(ispec)=Tnew(D_DECL(i,j,k),dencomp+1+ispec)
         if (massfrac_parm(ispec).ge.zero) then
          ! do nothing
         else
          print *,"massfrac_parm(ispec) invalid"
          stop
         endif
        enddo

        imattype=fort_material_type(im) 
          ! in otherwords, find c_v for material im
        call DeDT_material(den_local(im),massfrac_parm, &
         T_local(im), &
         DeDT_local(im),imattype,im)
       enddo ! im=1..nmat

       do dirloc=1,SDIM
        xsten_cell(dirloc)=xsten(0,dirloc)
       enddo

       ! MEHDI VAHAB HEAT SOURCE
       ! for right flux boundary condition:  
       !   T_new = T_old + dt (-q_right + qleft)/dx
       ! T^new=T^* + dt * Q/(rho cv)
       ! Q units: J/(m^3 s)
       ! get_local_heat_source in PROB.F90
       call get_local_heat_source( &
         time,dt, &
         nmat, &
         xsten_cell, &
         xsten, &  ! xsten(-nhalf:nhalf,SDIM) xsten(0,dir)=cell center
                   ! xsten(1,dir)=xcell + dx/2   xsten(-1,dir)=xcell-dx/2
                   ! xsten(2,dir)=xcell + dx     xsten(-2,dir)=xcell-dx
         nhalf, &
         temperature_source, &
         temperature_source_cen, &
         temperature_source_rad, &
         LS,VFRAC,T_local,den_local, &
         DeDT_local, &
         heat_source_local)

       heat_source_total=zero
       vfrac_total=zero 
       do im=1,nmat
        heat_source_total=heat_source_total+VFRAC(im)*heat_source_local(im)
        vfrac_total=vfrac_total+VFRAC(im)
       enddo

       if (vfrac_total.gt.zero) then
        ! do nothing
       else
        print *,"vfrac_total invalid"
        stop
       endif

       heat_source_total=heat_source_total/vfrac_total
 
         ! DeDTinverse = 1/(rho cv)
       do im=1,nmat
        T_local(im)=T_local(im)+ &
          dt*DeDTinverse(D_DECL(i,j,k),1)*heat_source_total
        if (1.eq.0) then
         if (heat_source_local(im).ne.zero) then
          print *,"x,im,heat_source_local ",xsten(0,1),xsten(0,2), &
           im,heat_source_local(im)
         endif
        endif

        dencomp=(im-1)*num_state_material+1
        Tnew(D_DECL(i,j,k),dencomp+1)=T_local(im)
       enddo ! im=1..nmat

      enddo
      enddo
      enddo ! i,j,k
 
      return
      end subroutine FORT_HEATSOURCE


         ! rhoinverse is 1/den
      subroutine FORT_SEMDELTAFORCE( &
       nstate, &
       nfluxSEM, &
       nstate_SDC, &
       nmat, &
       project_option, &
       xlo,dx,  &
       deltafab,DIMS(deltafab), &
       maskSEM,DIMS(maskSEM), &
       rhoinverse, &
       DIMS(rhoinverse), &
       DeDTinverse, &
       DIMS(DeDTinverse), &
       velnew,DIMS(velnew), &
       tilelo,tilehi, &
       fablo,fabhi,bfact,level, &
       dt)
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE


      INTEGER_T, intent(in) :: nmat,nstate
      INTEGER_T, intent(in) :: nfluxSEM
      INTEGER_T, intent(in) :: nstate_SDC
      INTEGER_T, intent(in) :: project_option
      INTEGER_T, intent(in) :: level
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T :: i,j,k
      INTEGER_T :: veldir
      INTEGER_T, intent(in) :: DIMDEC(deltafab)
      INTEGER_T, intent(in) :: DIMDEC(maskSEM)
      INTEGER_T, intent(in) :: DIMDEC(rhoinverse)
      INTEGER_T, intent(in) :: DIMDEC(DeDTinverse)
      INTEGER_T, intent(in) :: DIMDEC(velnew)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: deltafab(DIMV(deltafab),nstate_SDC)
      REAL_T, intent(in) :: maskSEM(DIMV(maskSEM))
      REAL_T, intent(in) :: rhoinverse(DIMV(rhoinverse),nmat+1)
      REAL_T, intent(in) :: DeDTinverse(DIMV(DeDTinverse),nmat+1)
      REAL_T, intent(inout) :: velnew(DIMV(velnew),nstate)
      REAL_T, intent(in) :: dt
      INTEGER_T idst,isrc,im
      INTEGER_T local_maskSEM


      if (bfact.lt.1) then
       print *,"bfact invalid54"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (nstate.ne.(SDIM+1)+ &
          nmat*(num_state_material+ngeom_raw)+1) then
       print *,"nstate invalid"
       stop
      endif
      if (dt.le.zero) then
       print *,"dt invalid"
       stop
      endif
      if (nfluxSEM.ne.SDIM+num_state_base) then
       print *,"nfluxSEM invalid sem delta force"
       stop
      endif
       ! (1) I scheme
       ! (2) temperature
       ! (3) viscosity
       ! (4) div(up)
       ! (5) pressure gradient
       ! (6) momentum force
      if (nstate_SDC.ne.nfluxSEM+1+SDIM+1+SDIM+SDIM) then
       print *,"nstate_SDC invalid"
       stop
      endif
      if ((project_option.eq.0).or. &
          (project_option.eq.2).or. &
          (project_option.eq.3).or. &
          (project_option.eq.4)) then
       ! do nothing
      else
       print *,"project_option invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(deltafab),0,-1,7)
      call checkbound(fablo,fabhi, &
       DIMS(rhoinverse), &
       1,-1,7)
      call checkbound(fablo,fabhi, &
       DIMS(DeDTinverse), &
       1,-1,7)
      call checkbound(fablo,fabhi,DIMS(velnew),1,-1,7)
      call checkbound(fablo,fabhi,DIMS(maskSEM),1,-1,1264)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       local_maskSEM=NINT(maskSEM(D_DECL(i,j,k)))
       if ((local_maskSEM.ge.1).and. &
           (local_maskSEM.le.nmat)) then

        if ((project_option.eq.3).or. &
            (project_option.eq.4)) then ! viscosity or -momentum force

         do veldir=1,SDIM
          velnew(D_DECL(i,j,k),veldir)= &
            velnew(D_DECL(i,j,k),veldir)- &
            rhoinverse(D_DECL(i,j,k),1)* &
            deltafab(D_DECL(i,j,k),veldir)
         enddo ! veldir=1..sdim

        else if ((project_option.eq.2).or. & ! thermal conduction
                 (project_option.eq.0)) then ! div(up), grad p

         isrc=1

         do im=1,nmat 
          idst=(SDIM+1)+ &
           (im-1)*num_state_material+2
          velnew(D_DECL(i,j,k),idst)= &
           velnew(D_DECL(i,j,k),idst)- &
           DeDTinverse(D_DECL(i,j,k),1)*deltafab(D_DECL(i,j,k),isrc)
         enddo ! im

         if (project_option.eq.0) then ! grad p

          do veldir=1,SDIM
           velnew(D_DECL(i,j,k),veldir)= &
             velnew(D_DECL(i,j,k),veldir)- &
             rhoinverse(D_DECL(i,j,k),1)* &
             deltafab(D_DECL(i,j,k),isrc+veldir)
          enddo ! veldir
         endif
        else
         print *,"project_option invalid"
         stop
        endif
       else if (local_maskSEM.eq.0) then
        ! do nothing
       else
        print *,"local_maskSEM invalid"
        stop
       endif

      enddo
      enddo
      enddo ! i,j,k

      return
      end subroutine FORT_SEMDELTAFORCE



      subroutine FORT_SEMDELTAFORCE_FACE( &
       dir, &
       faceden_index, &
       ncphys, &
       xlo,dx,  &
       deltafab,DIMS(deltafab), &
       maskSEM,DIMS(maskSEM), &
       xface,DIMS(xface), &
       xmac,DIMS(xmac), &
       tilelo,tilehi, &
       fablo,fabhi,bfact,level, &
       dt)
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: dir
      INTEGER_T, intent(in) :: faceden_index
      INTEGER_T, intent(in) :: ncphys
      INTEGER_T, intent(in) :: level
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T :: i,j,k
      INTEGER_T :: ii,jj,kk
      INTEGER_T, intent(in) :: DIMDEC(deltafab)
      INTEGER_T, intent(in) :: DIMDEC(maskSEM)
      INTEGER_T, intent(in) :: DIMDEC(xface)
      INTEGER_T, intent(in) :: DIMDEC(xmac)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: deltafab(DIMV(deltafab))
      REAL_T, intent(in) :: maskSEM(DIMV(maskSEM))
      REAL_T, intent(in) :: xface(DIMV(xface),ncphys)
      REAL_T, intent(inout) :: xmac(DIMV(xmac))
      REAL_T, intent(in) :: dt
      INTEGER_T maskleft,maskright
      INTEGER_T nmat

      nmat=num_materials

      if (bfact.lt.1) then
       print *,"bfact invalid55"
       stop
      endif
      if (faceden_index.ne.2) then
       print *,"faceden_index invalid"
       stop
      endif
      if (ncphys.lt.faceden_index+1) then
       print *,"ncphys invalid"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (dt.le.zero) then
       print *,"dt invalid"
       stop
      endif
      if ((dir.lt.0).or.(dir.ge.SDIM)) then
       print *,"dir invalid sem delta force face"
       stop
      endif

      ii=0
      jj=0
      kk=0
      if (dir.eq.0) then
       ii=1
      else if (dir.eq.1) then
       jj=1
      else if ((dir.eq.2).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"dir invalid sem delta force face 2"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(deltafab),0,dir,7)
      call checkbound(fablo,fabhi,DIMS(xface),0,dir,7)
      call checkbound(fablo,fabhi,DIMS(xmac),0,dir,7)
      call checkbound(fablo,fabhi,DIMS(maskSEM),1,-1,1264)

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
       growlo,growhi,0,dir,23)

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
       maskleft=NINT(maskSEM(D_DECL(i-ii,j-jj,k-kk)))
       maskright=NINT(maskSEM(D_DECL(i,j,k)))
       if ((maskleft.lt.0).or.(maskleft.gt.nmat)) then
        print *,"maskleft invalid"
        stop
       endif
       if ((maskright.lt.0).or.(maskright.gt.nmat)) then
        print *,"maskright invalid"
        stop
       endif
       if (((maskleft.ge.1).and.(maskleft.le.nmat)).or. &
           ((maskright.ge.1).and.(maskright.le.nmat))) then

        if ((maskleft.ne.0).and.(maskright.ne.0)) then
         if (maskleft.ne.maskright) then
          print *,"cannot have maskleft.ne.maskright"
          stop
         endif
        endif

        xmac(D_DECL(i,j,k))= &
          xmac(D_DECL(i,j,k))- &
          xface(D_DECL(i,j,k),faceden_index+1)* &
          deltafab(D_DECL(i,j,k))

       else if ((maskleft.eq.0).and.(maskright.eq.0)) then
        ! do nothing
       else
        print *,"maskleft or right invalid"
       endif

      enddo
      enddo
      enddo ! i,j,k

      return
      end subroutine FORT_SEMDELTAFORCE_FACE

       ! called from: NavierStokes::update_SEM_delta_force (NavierStokes.cpp)
       !   which is
       ! called from: NavierStokes::update_SEM_forces (MacProj.cpp)
       !   which is 
       ! called from: NavierStokes::update_SEM_forcesALL (MacProj.cpp)
       !   which is 
       ! called from: NavierStokes::veldiffuseALL (NavierStokes3.cpp)
       !   or 
       ! called from: NavierStokes::do_the_advance (NavierStokes3.cpp)
      subroutine FORT_UPDATESEMFORCE( &
       ns_time_order, &
       slab_step, &
       nsolve, &
       update_spectral, &
       update_stable, &
       nstate, &
       nfluxSEM, &
       nstate_SDC, &
       nmat, &
       project_option, &
       xlo,dx,  &
       gpfab,DIMS(gpfab), &
       divfab,DIMS(divfab), &
       hoopfab,DIMS(hoopfab), &
       HOfab,DIMS(HOfab), &
       LOfab,DIMS(LOfab), &
       maskSEM,DIMS(maskSEM), &
       tilelo,tilehi, &
       fablo,fabhi,bfact,level, &
       dt)
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T ns_time_order
      INTEGER_T slab_step
      INTEGER_T nsolve
      INTEGER_T update_spectral
      INTEGER_T update_stable
      INTEGER_T nmat,nstate
      INTEGER_T nfluxSEM
      INTEGER_T nstate_SDC
      INTEGER_T project_option,level
      REAL_T xlo(SDIM),dx(SDIM)
      INTEGER_T i,j,k
      INTEGER_T veldir
      INTEGER_T DIMDEC(gpfab)
      INTEGER_T DIMDEC(divfab)
      INTEGER_T DIMDEC(hoopfab)
      INTEGER_T DIMDEC(HOfab)
      INTEGER_T DIMDEC(LOfab)
      INTEGER_T DIMDEC(maskSEM)
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T bfact
      REAL_T gpfab(DIMV(gpfab),SDIM)
      REAL_T divfab(DIMV(divfab),nsolve)
      REAL_T hoopfab(DIMV(hoopfab),nsolve)
      REAL_T HOfab(DIMV(HOfab),nstate_SDC)
      REAL_T LOfab(DIMV(LOfab),nstate_SDC)
      REAL_T maskSEM(DIMV(maskSEM))
      REAL_T dt
      INTEGER_T local_maskSEM
      INTEGER_T velcomp


       ! in: UPDATESEMFORCE
      if ((ns_time_order.ge.2).and.(ns_time_order.le.32)) then
       ! do nothing
      else
       print *,"ns_time_order invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid56"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (nstate.ne.(SDIM+1)+ &
          nmat*(num_state_material+ngeom_raw)+1) then
       print *,"nstate invalid"
       stop
      endif
      if (dt.ne.one) then
       print *,"dt invalid"
       stop
      endif
      if (nfluxSEM.ne.SDIM+num_state_base) then
       print *,"nfluxSEM invalid update sem force"
       stop
      endif
       ! I-scheme,thermal conduction,viscosity,div(up),gp,-momforce
      if (nstate_SDC.ne.nfluxSEM+1+SDIM+1+SDIM+SDIM) then
       print *,"nstate_SDC invalid"
       stop
      endif
      if ((project_option.eq.0).or. &
          (project_option.eq.2).or. &
          (project_option.eq.3).or. &
          (project_option.eq.4)) then
       ! do nothing
      else
       print *,"project_option invalid"
       stop
      endif
      if ((nsolve.ne.1).and. &
          (nsolve.ne.SDIM)) then
       print *,"nsolve invalid 1"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(gpfab),0,-1,7)
      call checkbound(fablo,fabhi,DIMS(divfab),0,-1,7)
      call checkbound(fablo,fabhi,DIMS(hoopfab),0,-1,7)
      call checkbound(fablo,fabhi,DIMS(HOfab),0,-1,7)
      call checkbound(fablo,fabhi,DIMS(LOfab),0,-1,7)
      call checkbound(fablo,fabhi,DIMS(maskSEM),1,-1,1264)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
       local_maskSEM=NINT(maskSEM(D_DECL(i,j,k)))
       if ((local_maskSEM.ge.1).and. &
           (local_maskSEM.le.nmat)) then

        if (project_option.eq.3) then !viscosity

         do veldir=1,SDIM

          if (update_spectral.eq.1) then
           if ((slab_step.lt.-1).or.(slab_step.ge.bfact_time_order)) then
            print *,"slab_step invalid"
            stop
           endif
           ! I-scheme,thermal conduction,viscosity,div(up),gp,-force
           ! HOfab=-div 2 mu D-HOOP_FORCE_MARK_MF
           HOfab(D_DECL(i,j,k),nfluxSEM+1+veldir)= &
            divfab(D_DECL(i,j,k),veldir)- &
            hoopfab(D_DECL(i,j,k),veldir)
          else if (update_spectral.eq.0) then
           ! do nothing
          else
           print *,"update_spectral invalid"
           stop
          endif

          if (update_stable.eq.1) then
           if ((slab_step.lt.0).or.(slab_step.ge.bfact_time_order)) then
            print *,"slab_step invalid"
            stop
           endif
           ! I-scheme,thermal conduction,viscosity,div(up),gp,-force
           ! LOfab=-div 2 mu D-HOOP_FORCE_MARK_MF
           LOfab(D_DECL(i,j,k),nfluxSEM+1+veldir)= &
            divfab(D_DECL(i,j,k),veldir)- &
            hoopfab(D_DECL(i,j,k),veldir)
          else if (update_stable.eq.0) then
           ! do nothing
          else
           print *,"update_stable invalid"
           stop
          endif
         enddo ! veldir=1..sdim

        else if (project_option.eq.4) then ! -momforce

         do veldir=1,SDIM

          if (update_spectral.eq.1) then
           if ((slab_step.lt.-1).or.(slab_step.ge.bfact_time_order)) then
            print *,"slab_step invalid"
            stop
           endif
           ! I-scheme,thermal conduction,viscosity,div(up),gp,-momforce
           ! HOfab=-momforce
           HOfab(D_DECL(i,j,k),nfluxSEM+1+SDIM+1+SDIM+veldir)= &
             divfab(D_DECL(i,j,k),veldir)
          else if (update_spectral.eq.0) then
           ! do nothing
          else
           print *,"update_spectral invalid"
           stop
          endif

          if (update_stable.eq.1) then
           if ((slab_step.lt.0).or.(slab_step.ge.bfact_time_order)) then
            print *,"slab_step invalid"
            stop
           endif
           ! I-scheme,thermal conduction,viscosity,div(up),gp,-momforce
           ! LOfab=-momforce
           LOfab(D_DECL(i,j,k),nfluxSEM+1+SDIM+1+SDIM+veldir)= &
            divfab(D_DECL(i,j,k),veldir)
          else if (update_stable.eq.0) then
           ! do nothing
          else
           print *,"update_stable invalid"
           stop
          endif
         enddo ! veldir=1..sdim
 
        else if (project_option.eq.2) then  ! temperature diffusion

         velcomp=1

         if (update_spectral.eq.1) then
          if ((slab_step.lt.-1).or.(slab_step.ge.bfact_time_order)) then
           print *,"slab_step invalid"
           stop
          endif

          ! I-scheme,thermal conduction,viscosity,div(up),gp,-momforce
          ! HOfab=-div k grad T-THERMAL_FORCE_MF
          HOfab(D_DECL(i,j,k),nfluxSEM+1)= &
           divfab(D_DECL(i,j,k),velcomp)- &
           hoopfab(D_DECL(i,j,k),velcomp)
         else if (update_spectral.eq.0) then
          ! do nothing
         else
          print *,"update_spectral invalid"
          stop
         endif

         if (update_stable.eq.1) then
          if ((slab_step.lt.0).or.(slab_step.ge.bfact_time_order)) then
           print *,"slab_step invalid"
           stop
          endif

          ! I-scheme,thermal conduction,viscosity,div(up),gp,-momforce
          ! LOfab=-div k grad T-THERMAL_FORCE_MF
          LOfab(D_DECL(i,j,k),nfluxSEM+1)= &
           divfab(D_DECL(i,j,k),velcomp)- &
           hoopfab(D_DECL(i,j,k),velcomp)
         else if (update_stable.eq.0) then
          ! do nothing
         else
          print *,"update_stable invalid"
          stop
         endif

        else if (project_option.eq.0) then

         do veldir=1,SDIM

          if (update_spectral.eq.1) then
           if ((slab_step.lt.-1).or.(slab_step.ge.bfact_time_order)) then
            print *,"slab_step invalid"
            stop
           endif
           ! I-scheme,thermal conduction,viscosity,div(up),gp,-momforce
           ! HOfab=grad p
           HOfab(D_DECL(i,j,k),nfluxSEM+1+SDIM+1+veldir)= &
             gpfab(D_DECL(i,j,k),veldir)
          else if (update_spectral.eq.0) then
           ! do nothing
          else
           print *,"update_spectral invalid"
           stop
          endif

          if (update_stable.eq.1) then
           if ((slab_step.lt.0).or.(slab_step.ge.bfact_time_order)) then
            print *,"slab_step invalid"
            stop
           endif
           ! I-scheme,thermal conduction,viscosity,div(up),gp,-momforce
           ! LOfab=grad p
           LOfab(D_DECL(i,j,k),nfluxSEM+1+SDIM+1+veldir)= &
             gpfab(D_DECL(i,j,k),veldir)
          else if (update_stable.eq.0) then
           ! do nothing
          else
           print *,"update_stable invalid"
           stop
          endif

         enddo ! veldir=1..sdim

         velcomp=1

         if (update_spectral.eq.1) then
          if ((slab_step.lt.-1).or.(slab_step.ge.bfact_time_order)) then
           print *,"slab_step invalid"
           stop
          endif
          ! I-scheme,thermal conduction,viscosity,div(up),gp,-momforce
          ! HOfab=div (up)
          HOfab(D_DECL(i,j,k),nfluxSEM+1+SDIM+1)=divfab(D_DECL(i,j,k),1)
         else if (update_spectral.eq.0) then
          ! do nothing
         else
          print *,"update_spectral invalid"
          stop
         endif

         if (update_stable.eq.1) then
          if ((slab_step.lt.0).or.(slab_step.ge.bfact_time_order)) then
           print *,"slab_step invalid"
           stop
          endif
          ! I-scheme,thermal conduction,viscosity,div(up),gp,-momforce
          ! LOfab=div (up)
          LOfab(D_DECL(i,j,k),nfluxSEM+1+SDIM+1)=divfab(D_DECL(i,j,k),1)
         else if (update_stable.eq.0) then
          ! do nothing
         else
          print *,"update_stable invalid"
          stop
         endif

        else
         print *,"project_option invalid"
         stop
        endif

       else if (local_maskSEM.eq.0) then
        ! do nothing
       else
        print *,"local_maskSEM invalid"
        stop
       endif
      enddo
      enddo
      enddo

      return
      end subroutine FORT_UPDATESEMFORCE

       ! delta_MF, stableF_MF, spectralF_MF are initialized
       ! to zero in NavierStokes::init_delta_SDC()
      subroutine FORT_UPDATESEMFORCE_FACE( &
       project_option, &
       ns_time_order, &
       dir, &
       slab_step, &
       update_spectral, &
       update_stable, &
       nmat, &
       xlo,dx,  &
       gpfab,DIMS(gpfab), &
       HOfab,DIMS(HOfab), &
       LOfab,DIMS(LOfab), &
       maskSEM,DIMS(maskSEM), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       level, &
       dt)
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: project_option
      INTEGER_T, intent(in) :: ns_time_order
      INTEGER_T, intent(in) :: dir
      INTEGER_T, intent(in) :: slab_step
      INTEGER_T, intent(in) :: update_spectral
      INTEGER_T, intent(in) :: update_stable
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: level
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T, intent(in) :: DIMDEC(gpfab)
      INTEGER_T, intent(in) :: DIMDEC(HOfab)
      INTEGER_T, intent(in) :: DIMDEC(LOfab)
      INTEGER_T, intent(in) :: DIMDEC(maskSEM)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: gpfab(DIMV(gpfab))
      REAL_T, intent(out) :: HOfab(DIMV(HOfab))
      REAL_T, intent(out) :: LOfab(DIMV(LOfab))
      REAL_T, intent(in) :: maskSEM(DIMV(maskSEM))
      REAL_T, intent(in) :: dt
      INTEGER_T maskleft,maskright
      INTEGER_T :: i,j,k
      INTEGER_T ii,jj,kk
      INTEGER_T im_crit

      if (project_option.eq.0) then
       ! do nothing
      else if (project_option.eq.2) then ! thermal diffusion
       print *,"thermal diffusion force is only cell centered"
       stop
      else if (project_option.eq.3) then ! viscosity
       print *,"viscosity force is only cell centered"
       stop
      else
       print *,"project_option invalid"
       stop
      endif

      if ((ns_time_order.ge.2).and.(ns_time_order.le.32)) then
       ! do nothing
      else
       print *,"ns_time_order invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid57"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (dt.ne.one) then
       print *,"dt invalid"
       stop
      endif
      if ((dir.lt.0).or.(dir.ge.SDIM)) then
       print *,"dir invalid update sem force face"
       stop
      endif

      ii=0
      jj=0
      kk=0
      if (dir.eq.0) then
       ii=1
      else if (dir.eq.1) then
       jj=1
      else if ((dir.eq.2).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"dir invalid update sem force face 2"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(gpfab),0,dir,7)
      call checkbound(fablo,fabhi,DIMS(HOfab),0,dir,7)
      call checkbound(fablo,fabhi,DIMS(LOfab),0,dir,7)
      call checkbound(fablo,fabhi,DIMS(maskSEM),1,-1,1264)

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
       growlo,growhi,0,dir,24)
 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       maskleft=NINT(maskSEM(D_DECL(i-ii,j-jj,k-kk)))
       maskright=NINT(maskSEM(D_DECL(i,j,k)))
       if ((maskleft.lt.0).or.(maskleft.gt.nmat)) then
        print *,"maskleft invalid"
        stop
       endif
       if ((maskright.lt.0).or.(maskright.gt.nmat)) then
        print *,"maskright invalid"
        stop
       endif
       im_crit=0
       if ((maskleft.eq.0).and.(maskright.eq.0)) then
        ! do nothing
       else if (maskleft.eq.0) then
        im_crit=maskright
       else if (maskright.eq.0) then
        im_crit=maskleft
       else if (maskleft.eq.maskright) then
        im_crit=maskleft
       else
        print *,"maskleft or maskright invalid"
        stop
       endif

       if ((im_crit.ge.1).and.(im_crit.le.nmat)) then

        if (update_spectral.eq.1) then
         if ((slab_step.lt.-1).or.(slab_step.ge.bfact_time_order)) then
          print *,"slab_step invalid"
          stop
         endif 
         HOfab(D_DECL(i,j,k))=gpfab(D_DECL(i,j,k))
        else if (update_spectral.eq.0) then
         ! do nothing
        else
         print *,"update_spectral invalid"
         stop
        endif

        if (update_stable.eq.1) then
         if ((slab_step.lt.0).or.(slab_step.ge.bfact_time_order)) then
          print *,"slab_step invalid"
          stop
         endif

         LOfab(D_DECL(i,j,k))=gpfab(D_DECL(i,j,k))
        else if (update_stable.eq.0) then
         ! do nothing
        else
         print *,"update_stable invalid"
         stop
        endif

       else if (im_crit.eq.0) then
        ! do nothing
       else
        print *,"im_crit invalid"
        stop
       endif

      enddo
      enddo
      enddo ! i,j,k

      return
      end subroutine FORT_UPDATESEMFORCE_FACE


      subroutine FORT_SDC_TIME_QUAD( &
       HOncomp, &
       LOncomp, &
       delta_ncomp, &
       nstate, &
       nfluxSEM, &
       nstate_SDC, &
       nmat, &
       xlo,dx,  &
       deltafab,DIMS(deltafab), &
       HOfab,DIMS(HOfab), &
       LOfab,DIMS(LOfab), &
       maskSEM,DIMS(maskSEM), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       level, &
       finest_level, &
       dt)
      use LegendreNodes
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T HOncomp
      INTEGER_T LOncomp
      INTEGER_T delta_ncomp
      INTEGER_T nmat,nstate
      INTEGER_T nfluxSEM
      INTEGER_T nstate_SDC
      INTEGER_T level
      INTEGER_T finest_level
      REAL_T xlo(SDIM),dx(SDIM)
      INTEGER_T i,j,k
      INTEGER_T DIMDEC(deltafab)
      INTEGER_T DIMDEC(HOfab)
      INTEGER_T DIMDEC(LOfab)
      INTEGER_T DIMDEC(maskSEM)
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T bfact
      REAL_T deltafab(DIMV(deltafab),delta_ncomp)
      REAL_T HOfab(DIMV(HOfab),HOncomp)
      REAL_T LOfab(DIMV(LOfab),LOncomp)
      REAL_T maskSEM(DIMV(maskSEM))
      REAL_T dt
      INTEGER_T local_maskSEM
      INTEGER_T slab_step
      INTEGER_T ibase,ibase2,icomp
      INTEGER_T kinterval,jstencil,iquad,i1
      REAL_T force_integral,dt_sub
      REAL_T sanity_sum
      REAL_T GQws(0:bfact_time_order,0:bfact_time_order-1, &
                  1:bfact_time_order)
      REAL_T GQwsQUAD(0:bfact_time_order,1:bfact_time_order)
      REAL_T yGL(0:bfact_time_order)
      REAL_T ydiff(1:bfact_time_order)
      INTEGER_T ok_to_update(nstate_SDC)
      INTEGER_T imattype

      if (bfact.lt.1) then
       print *,"bfact invalid58"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (nstate.ne.(SDIM+1)+ &
          nmat*(num_state_material+ngeom_raw)+1) then
       print *,"nstate invalid"
       stop
      endif
      if (dt.le.zero) then
       print *,"dt invalid"
       stop
      endif
      if (nfluxSEM.ne.SDIM+num_state_base) then
       print *,"nfluxSEM invalid sdc time quad"
       stop
      endif 
       ! I-scheme, thermal conduction, viscosity, -div(up) ,-grad p,-momforce
      if (nstate_SDC.ne.nfluxSEM+1+SDIM+1+SDIM+SDIM) then
       print *,"nstate_SDC invalid"
       stop
      endif
      if (HOncomp.ne.nstate_SDC*(bfact_time_order+1)) then
       print *,"HOncomp invalid"
       stop
      endif
      if (LOncomp.ne.nstate_SDC*bfact_time_order) then
       print *,"LOncomp invalid"
       stop
      endif
      if (delta_ncomp.ne.nstate_SDC*bfact_time_order) then
       print *,"delta_ncomp invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid sdc time quad"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(HOfab),0,-1,7)
      call checkbound(fablo,fabhi,DIMS(LOfab),0,-1,7)
      call checkbound(fablo,fabhi,DIMS(deltafab),0,-1,7)
      call checkbound(fablo,fabhi,DIMS(maskSEM),1,-1,1264)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)

      call SDC_GQweights(bfact_time_order,bfact_time_order,GQws)

      do kinterval=1,bfact_time_order
       sanity_sum=zero
       do jstencil=0,bfact_time_order
        GQwsQUAD(jstencil,kinterval)=zero
        do iquad=0,bfact_time_order-1
         GQwsQUAD(jstencil,kinterval)=GQwsQUAD(jstencil,kinterval)+ &
          cache_gauss_w(bfact_time_order,iquad,TMTYPE)* &
          GQws(jstencil,iquad,kinterval)
        enddo ! iquad
        sanity_sum=sanity_sum+GQwsQUAD(jstencil,kinterval)
       enddo  ! jstencil
       if (abs(sanity_sum-two)>1.0E-12) then
        print *,"SDC sanity check failed"
        stop
       endif
      enddo  ! kinterval

      do i1=0,bfact_time_order
       yGL(i1)=cache_gauss_lobatto(bfact_time_order,i1,TMTYPE)
      enddo
      do i1=1,bfact_time_order
       ydiff(i1)=yGL(i1)-yGL(i1-1)
       if (ydiff(i1).le.zero) then
        print *,"yGL bust"
        stop
       endif
      enddo ! i1

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       local_maskSEM=NINT(maskSEM(D_DECL(i,j,k)))

       if (local_maskSEM.eq.0) then
        ! do nothing
       else if ((local_maskSEM.ge.1).and.(local_maskSEM.le.nmat)) then
        ! do nothing
       else
        print *,"local_maskSEM invalid"
        stop
       endif

       do slab_step=1,bfact_time_order
        ibase=(slab_step-1)*nstate_SDC
        do icomp=1,nstate_SDC
         deltafab(D_DECL(i,j,k),ibase+icomp)=zero
        enddo
       enddo

       if (local_maskSEM.eq.0) then
        ! do nothing
       else if ((local_maskSEM.ge.1).and. &
                (local_maskSEM.le.nmat)) then

        do slab_step=1,bfact_time_order

         ibase=(slab_step-1)*nstate_SDC

         do icomp=1,nstate_SDC
          ok_to_update(icomp)=0
         enddo

         do icomp=1,SDIM
          ok_to_update(icomp)=1 ! velocity (I-scheme)
         enddo

         icomp=SDIM+1

         imattype=fort_material_type(local_maskSEM)
         if (imattype.eq.999) then
          ! do nothing
         else if (imattype.eq.0) then
          ! do nothing
         else if ((imattype.gt.0).and. &
                  (imattype.le.MAX_NUM_EOS)) then
          ok_to_update(icomp)=1 ! density (I-scheme)
         else
          print *,"imattype invalid fort_sdc_time_quad"
          stop
         endif

         icomp=SDIM+2
         ok_to_update(icomp)=1 ! temperature (I-scheme)

         if (icomp.ne.nfluxSEM) then
          print *,"icomp invalid"
          stop
         endif

         ! thermal conduction, viscosity, div(up), grad p,-momforce
         do icomp=nfluxSEM+1,nstate_SDC
          ok_to_update(icomp)=1
         enddo

         do icomp=1,nstate_SDC
          if (ok_to_update(icomp).eq.1) then
           force_integral=zero
           do jstencil=0,bfact_time_order
            ibase2=jstencil*nstate_SDC+icomp
            force_integral=force_integral+ &
             GQwsQUAD(jstencil,slab_step)* &
             HOfab(D_DECL(i,j,k),ibase2)
           enddo ! jstencil
           dt_sub=dt*ydiff(slab_step)/two

           deltafab(D_DECL(i,j,k),ibase+icomp)= &
            force_integral*dt_sub/two- &
            dt_sub*LOfab(D_DECL(i,j,k),ibase+icomp)
          else if (ok_to_update(icomp).eq.0) then
           ! do nothing
          else
           print *,"ok_to_update invalid"
           stop
          endif
         enddo ! icomp

        enddo ! slab_step

       else
        print *,"local_maskSEM invalid"
        stop
       endif
 
      enddo
      enddo
      enddo

      return
      end subroutine FORT_SDC_TIME_QUAD



      subroutine FORT_SDC_TIME_QUAD_FACE( &
       dir, &
       HOncomp, &
       LOncomp, &
       delta_ncomp, &
       nstate, &
       nfluxSEM, &
       nstate_SDC, &
       nmat, &
       xlo,dx,  &
       deltafab,DIMS(deltafab), &
       HOfab,DIMS(HOfab), &
       LOfab,DIMS(LOfab), &
       maskSEM,DIMS(maskSEM), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       level, &
       finest_level, &
       dt)
      use LegendreNodes
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T dir
      INTEGER_T HOncomp
      INTEGER_T LOncomp
      INTEGER_T delta_ncomp
      INTEGER_T nmat,nstate
      INTEGER_T nfluxSEM
      INTEGER_T nstate_SDC
      INTEGER_T level
      INTEGER_T finest_level
      REAL_T xlo(SDIM),dx(SDIM)
      INTEGER_T DIMDEC(deltafab)
      INTEGER_T DIMDEC(HOfab)
      INTEGER_T DIMDEC(LOfab)
      INTEGER_T DIMDEC(maskSEM)
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T bfact
      REAL_T deltafab(DIMV(deltafab),delta_ncomp)
      REAL_T HOfab(DIMV(HOfab),HOncomp)
      REAL_T LOfab(DIMV(LOfab),LOncomp)
      REAL_T maskSEM(DIMV(maskSEM))
      REAL_T dt
      INTEGER_T i,j,k
      INTEGER_T ii,jj,kk
      INTEGER_T slab_step
      INTEGER_T kinterval,jstencil,iquad,i1
      REAL_T force_integral,dt_sub
      REAL_T sanity_sum
      REAL_T GQws(0:bfact_time_order,0:bfact_time_order-1, &
                  1:bfact_time_order)
      REAL_T GQwsQUAD(0:bfact_time_order,1:bfact_time_order)
      REAL_T yGL(0:bfact_time_order)
      REAL_T ydiff(1:bfact_time_order)
      INTEGER_T maskleft,maskright

      if (bfact.lt.1) then
       print *,"bfact invalid59"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (nstate.ne.(SDIM+1)+ &
          nmat*(num_state_material+ngeom_raw)+1) then
       print *,"nstate invalid"
       stop
      endif
      if (dt.le.zero) then
       print *,"dt invalid"
       stop
      endif
      if (nfluxSEM.ne.SDIM+num_state_base) then
       print *,"nfluxSEM invalid sdc time quad face"
       stop
      endif
       ! I-scheme, thermal conductivity, viscosity, div(up), gp, -momforce
      if (nstate_SDC.ne.nfluxSEM+1+SDIM+1+SDIM+SDIM) then
       print *,"nstate_SDC invalid"
       stop
      endif
      if (HOncomp.ne.(bfact_time_order+1)) then
       print *,"HOncomp invalid"
       stop
      endif
      if (LOncomp.ne.bfact_time_order) then
       print *,"LOncomp invalid"
       stop
      endif
      if (delta_ncomp.ne.bfact_time_order) then
       print *,"delta_ncomp invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid sdc time quad face"
       stop
      endif

      if ((dir.lt.0).or.(dir.ge.SDIM)) then
       print *,"dir invalid sdc time quad face"
       stop
      endif

      ii=0
      jj=0
      kk=0
      if (dir.eq.0) then
       ii=1
      else if (dir.eq.1) then
       jj=1
      else if ((dir.eq.2).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"dir invalid sdc time quad face 2"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(HOfab),0,dir,7)
      call checkbound(fablo,fabhi,DIMS(LOfab),0,dir,7)
      call checkbound(fablo,fabhi,DIMS(deltafab),0,dir,7)
      call checkbound(fablo,fabhi,DIMS(maskSEM),1,-1,1264)

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,0,dir,25)

      call SDC_GQweights(bfact_time_order,bfact_time_order,GQws)

      do kinterval=1,bfact_time_order
       sanity_sum=zero
       do jstencil=0,bfact_time_order
        GQwsQUAD(jstencil,kinterval)=zero
        do iquad=0,bfact_time_order-1
         GQwsQUAD(jstencil,kinterval)=GQwsQUAD(jstencil,kinterval)+ &
          cache_gauss_w(bfact_time_order,iquad,TMTYPE)* &
          GQws(jstencil,iquad,kinterval)
        enddo ! iquad
        sanity_sum=sanity_sum+GQwsQUAD(jstencil,kinterval)
       enddo  ! jstencil
       if (abs(sanity_sum-two)>1.0E-12) then
        print *,"SDC sanity check failed"
        stop
       endif
      enddo  ! kinterval

      do i1=0,bfact_time_order
       yGL(i1)=cache_gauss_lobatto(bfact_time_order,i1,TMTYPE)
      enddo
      do i1=1,bfact_time_order
       ydiff(i1)=yGL(i1)-yGL(i1-1)
       if (ydiff(i1).le.zero) then
        print *,"yGL bust"
        stop
       endif
      enddo ! i1

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       maskleft=NINT(maskSEM(D_DECL(i-ii,j-jj,k-kk)))
       maskright=NINT(maskSEM(D_DECL(i,j,k)))
       if ((maskleft.lt.0).or.(maskleft.gt.nmat)) then
        print *,"maskleft invalid"
        stop
       endif
       if ((maskright.lt.0).or.(maskright.gt.nmat)) then
        print *,"maskright invalid"
        stop
       endif

       do slab_step=1,bfact_time_order
        deltafab(D_DECL(i,j,k),slab_step)=zero
       enddo
       if (((maskleft.ge.1).and.(maskleft.le.nmat)).or. &
           ((maskright.ge.1).and.(maskright.le.nmat))) then

        do slab_step=1,bfact_time_order

         force_integral=zero
         do jstencil=0,bfact_time_order
          force_integral=force_integral+ &
           GQwsQUAD(jstencil,slab_step)* &
           HOfab(D_DECL(i,j,k),jstencil+1)
         enddo ! jstencil
         dt_sub=dt*ydiff(slab_step)/two

         deltafab(D_DECL(i,j,k),slab_step)= &
          force_integral*dt_sub/two- &
          dt_sub*LOfab(D_DECL(i,j,k),slab_step)

        enddo ! slab_step

       endif
 
      enddo
      enddo
      enddo

      return
      end subroutine FORT_SDC_TIME_QUAD_FACE

      subroutine FORT_MAKETENSOR( &
       partid, & ! 0..num_materials_viscoelastic-1
       level, &
       finest_level, &
       MAC_grid_displacement, &
       ncomp_visc, &
       im_parm, & ! 0..nmat-1
       xlo,dx, &
       recon,DIMS(recon), &  
       xdfab,DIMS(xdfab), &
       ydfab,DIMS(ydfab), &
       zdfab,DIMS(zdfab), &
       visc,DIMS(visc), &
       tensor,DIMS(tensor), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       elastic_viscosity,etaS, &
       elastic_time, &
       viscoelastic_model, &
       polymer_factor, &
       irz,ngrow,nmat)
      use probcommon_module
      use global_utility_module
      use mass_transfer_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: partid !0..num_materials_viscoelastic-1
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: MAC_grid_displacement
      INTEGER_T, intent(in) :: ncomp_visc
      INTEGER_T, intent(in) :: im_parm
      INTEGER_T, intent(in) :: nmat
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T, intent(in) :: ngrow
      INTEGER_T, intent(in) :: DIMDEC(recon)
      INTEGER_T, intent(in) :: DIMDEC(xdfab)
      INTEGER_T, intent(in) :: DIMDEC(ydfab)
      INTEGER_T, intent(in) :: DIMDEC(zdfab)
      INTEGER_T, intent(in) :: DIMDEC(visc)
      INTEGER_T, intent(in) :: DIMDEC(tensor)
      INTEGER_T, intent(in) :: tilelo(SDIM), tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM), fabhi(SDIM)
      INTEGER_T :: growlo(3), growhi(3)
      INTEGER_T, intent(in) :: bfact

      REAL_T, intent(in) :: recon(DIMV(recon),nmat*ngeom_recon)
      REAL_T, intent(in) :: xdfab(DIMV(xdfab))
      REAL_T, intent(in) :: ydfab(DIMV(ydfab))
      REAL_T, intent(in) :: zdfab(DIMV(zdfab))
      REAL_T, intent(in) :: visc(DIMV(visc),ncomp_visc)
      REAL_T, intent(inout) :: tensor(DIMV(tensor),FORT_NUM_TENSOR_TYPE)

      REAL_T, intent(in) :: elastic_viscosity,etaS
      REAL_T, intent(in) :: elastic_time,polymer_factor
      INTEGER_T, intent(in) :: viscoelastic_model
      INTEGER_T, intent(in) :: irz

      INTEGER_T ii,jj
      REAL_T Q(3,3),TQ(3,3)
      INTEGER_T i,j,k
      INTEGER_T dir_flux,side_flux,dir_local,dir_XD
      INTEGER_T im_elastic_p1
      REAL_T xcenter(SDIM)
      REAL_T XDcenter(SDIM)
      REAL_T delta_flux
      REAL_T xflux(SDIM)
      REAL_T XDside(SDIM)
      REAL_T XDside_stencil(0:1,0:SDIM-1,SDIM)
      REAL_T xflux_stencil(0:1,0:SDIM-1,SDIM)
      REAL_T gradXDtensor(SDIM,SDIM) ! dir_xdisp,dir_space
      REAL_T avgXDtensor(SDIM,SDIM)
      REAL_T DISP_TEN(SDIM,SDIM) ! dir_x (displace), dir_space
      REAL_T hoop_22
       ! if use_A==0 then force is div(mu H Q)/rho
       ! if use_A==1 then force is div(mu H A)/rho
      INTEGER_T use_A  
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf

      nhalf=3

      if ((partid.ge.0).and. &
          (partid.lt.num_materials_viscoelastic)) then
       ! do nothing
      else
       print *,"partid invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid MAKETENSOR"
       stop
      endif
      if (finest_level.ne.fort_finest_level) then
       print *,"finest_level invalid MAKETENSOR"
       stop
      endif

      if ((viscoelastic_model.eq.0).or. & ! (visc-etaS)/lambda
          (viscoelastic_model.eq.1)) then ! (visc-etaS)
       ! For incompressible flow, the equations using Q or A are equivalent.
       ! For compressible flow, one should probably set: use_A=1.
       use_A=0
      else if (viscoelastic_model.eq.2) then ! elastic model
       use_A=0
      else if (viscoelastic_model.eq.3) then ! incremental model
       use_A=0
      else
       print *,"viscoelastic_model invalid"
       stop
      endif

      if (ngrow.ne.1) then
       print *,"ngrow invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if ((im_parm.lt.0).or. &
          (im_parm.ge.nmat).or. &
          (is_rigid(nmat,im_parm+1).eq.1)) then
       print *,"im_parm invalid26"
       stop
      endif

      im_elastic_p1=im_parm+1

      if (ncomp_visc.ne.3*nmat) then
       print *,"ncomp_visc invalid"
       stop
      endif

      if (MAC_grid_displacement.eq.0) then
       call checkbound(fablo,fabhi,DIMS(xdfab),2,-1,11)
       call checkbound(fablo,fabhi,DIMS(ydfab),2,-1,11)
       call checkbound(fablo,fabhi,DIMS(zdfab),2,-1,11)
      else if (MAC_grid_displacement.eq.1) then
       call checkbound(fablo,fabhi,DIMS(xdfab),1,0,11)
       call checkbound(fablo,fabhi,DIMS(ydfab),1,1,11)
       call checkbound(fablo,fabhi,DIMS(zdfab),1,SDIM-1,11)
      else
       print *,"MAC_grid_displacement invalid"
       stop
      endif
      call checkbound(fablo,fabhi,DIMS(recon),2,-1,1277)

      call checkbound(fablo,fabhi,DIMS(visc),ngrow,-1,11)
      call checkbound(fablo,fabhi,DIMS(tensor),ngrow,-1,8)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,ngrow) 

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridsten_level(xsten,i,j,k,level,nhalf)
       do dir_local=1,SDIM
        xcenter(dir_local)=xsten(0,dir_local)
       enddo

       do ii=1,3
       do jj=1,3
        Q(ii,jj)=zero
       enddo
       enddo
       Q(1,1)=tensor(D_DECL(i,j,k),1)
       Q(1,2)=tensor(D_DECL(i,j,k),2)
       Q(2,2)=tensor(D_DECL(i,j,k),3)
       Q(3,3)=tensor(D_DECL(i,j,k),4)
#if (AMREX_SPACEDIM==3)
       Q(1,3)=tensor(D_DECL(i,j,k),5)
       Q(2,3)=tensor(D_DECL(i,j,k),6)
#endif
       Q(2,1)=Q(1,2)
       Q(3,1)=Q(1,3)
       Q(3,2)=Q(2,3)

       if ((viscoelastic_model.eq.0).or. & ! (visc-etaS)/lambda
           (viscoelastic_model.eq.1)) then ! (visc-etaS)
        ! do nothing
       else if (viscoelastic_model.eq.3) then ! incremental model
        ! do nothing
       else if (viscoelastic_model.eq.2) then ! elastic model
        if (MAC_grid_displacement.eq.0) then
         ! do nothing
        else if (MAC_grid_displacement.eq.1) then
         do dir_flux=0,SDIM-1
         do side_flux=0,1
          do dir_local=1,SDIM
           xflux(dir_local)=xsten(0,dir_local)
          enddo
          if (side_flux.eq.0) then
           xflux(dir_flux+1)=xsten(-1,dir_flux+1)
          else if (side_flux.eq.1) then
           xflux(dir_flux+1)=xsten(1,dir_flux+1)
          else
           print *,"side_flux invalid"
           stop
          endif

           ! interpfab_XDISP declared in MASS_TRANSFER_3D.F90
          call interpfab_XDISP( &
            bfact, & ! determines positioning of Gauss Legendre nodes
            level, &
            finest_level, &
            dx, &
            xlo, &
            xflux, &
            im_elastic_p1, &!1..nmat(prescribed as a fluid in the inputs file)
            nmat, &
            partid, & ! 0..num_materials_viscoelastic-1
            fablo,fabhi, &
            xdfab,DIMS(xdfab), &
            ydfab,DIMS(ydfab), &
            zdfab,DIMS(zdfab), &
            recon,DIMS(recon), &
            XDside) ! XD(xflux),YD(xflux),ZD(xflux): XDside(SDIM)

          do dir_XD=1,SDIM
           XDside_stencil(side_flux,dir_flux,dir_XD)= &
             XDside(dir_XD) 
          enddo
          do dir_local=1,SDIM
           xflux_stencil(side_flux,dir_flux,dir_local)=xflux(dir_local)
          enddo ! dir_local=1..sdim
                    
         enddo ! side_flux=0,1
         enddo ! dir_flux=0..sdim-1

         do dir_XD=1,SDIM
          XDcenter(dir_XD)=zero
         enddo
         do dir_flux=0,SDIM-1
          delta_flux=xflux_stencil(1,dir_flux,dir_flux+1)- &
                     xflux_stencil(0,dir_flux,dir_flux+1)
          if (delta_flux.gt.zero) then
           do dir_XD=1,SDIM
            gradXDtensor(dir_XD,dir_flux+1)= &
             (XDside_stencil(1,dir_flux,dir_XD)- &
              XDside_stencil(0,dir_flux,dir_XD))/delta_flux
            avgXDtensor(dir_XD,dir_flux+1)= &
             half*(XDside_stencil(1,dir_flux,dir_XD)+ &
                   XDside_stencil(0,dir_flux,dir_XD)) 
            XDcenter(dir_XD)=XDcenter(dir_XD)+avgXDtensor(dir_XD,dir_flux+1)
           enddo ! dir_XD=1..sdim
          else
           print *,"delta_flux invalid"
           stop
          endif
         enddo ! dir_flux=0..sdim-1
         do dir_XD=1,SDIM
          XDcenter(dir_XD)=XDcenter(dir_XD)/SDIM
         enddo

          ! declared in: GLOBALUTIL.F90
         call stress_from_strain( &
          im_elastic_p1, & ! =1..nmat
          xcenter, &
          dx, &
          gradXDtensor, &
          XDcenter(1), &
          XDcenter(2), &
          DISP_TEN, &  ! dir_x (displace),dir_space
          hoop_22)  ! output: "theta-theta" component xdisp/r if RZ

         do ii=1,3
         do jj=1,3
          Q(ii,jj)=zero
         enddo
         enddo

         Q(1,1)=DISP_TEN(1,1)
         Q(1,2)=DISP_TEN(1,2)
         Q(2,2)=DISP_TEN(2,2)
         if (SDIM.eq.3) then
          Q(3,3)=DISP_TEN(SDIM,SDIM)
         else if (SDIM.eq.2) then
          if (levelrz.eq.0) then
           ! T33 (theta coordinate)
           Q(3,3)=zero
          else if (levelrz.eq.1) then
           ! T33 (theta coordinate)
           ! dX/dx + dX/dx
           Q(3,3)=two*hoop_22 ! 2 * (xdisp/r)
          else if (levelrz.eq.3) then
           ! T33 (z coordinate)
           Q(3,3)=zero
          else
           print *,"levelrz invalid"
           stop
          endif
         else
          print *,"dimension bust"
          stop
         endif
                    
         if (SDIM.eq.3) then 
          Q(1,SDIM)=DISP_TEN(1,SDIM)
          Q(2,SDIM)=DISP_TEN(2,SDIM)
         endif
         Q(2,1)=Q(1,2)
         Q(3,1)=Q(1,3)
         Q(3,2)=Q(2,3)

        else
         print *,"MAC_grid_displacement invalid"
         stop
        endif
       else
        print *,"viscoelastic_model invalid"
        stop
       endif

       if (use_A.eq.0) then
        ! do nothing
       else if (use_A.eq.1) then
        do ii=1,3
         Q(ii,ii)=Q(ii,ii)+one
        enddo
       else
        print *,"use_A invalid"
        stop
       endif

        ! viscoelastic_model==0: 
        !   visc(nmat+im_parm+1)=(eta/lambda_mod)*visc_coef
        ! For a purely elastic model, it is assumed that the Deborah number
        ! is \infty.
        ! For a purely viscous model, it is assumed that the Deborah number is
        ! zero.
        ! viscoelastic_model==2:  (displacement gradient model)
        !   visc(nmat+im_parm+1)=eta*visc_coef
        ! viscoelastic_model==3:  (incremental model)
        !   visc(nmat+im_parm+1)=eta*visc_coef
       do ii=1,3
       do jj=1,3
         TQ(ii,jj)=Q(ii,jj)*visc(D_DECL(i,j,k),nmat+im_parm+1)
       enddo
       enddo

       tensor(D_DECL(i,j,k),1)=TQ(1,1)
       tensor(D_DECL(i,j,k),2)=TQ(1,2)
       tensor(D_DECL(i,j,k),3)=TQ(2,2)
       tensor(D_DECL(i,j,k),4)=TQ(3,3)
#if (AMREX_SPACEDIM==3)
       tensor(D_DECL(i,j,k),5)=TQ(1,3)
       tensor(D_DECL(i,j,k),6)=TQ(2,3)
#endif
      enddo
      enddo
      enddo

      return
      end subroutine MAKETENSOR

      subroutine FORT_COPY_VEL_ON_SIGN( &
       im_part, &
       nparts, &
       partid, &
       ngrowFSI, &
       nFSI, &
       nFSI_sub, &
       xlo,dx, &
       snew,DIMS(snew), &
       fsi,DIMS(fsi), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       nmat,nstate)
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T im_part
      INTEGER_T nparts
      INTEGER_T partid
      INTEGER_T ngrowFSI
      INTEGER_T nFSI
      INTEGER_T nFSI_sub
      INTEGER_T nmat
      INTEGER_T nstate
      REAL_T xlo(SDIM),dx(SDIM)
      INTEGER_T DIMDEC(snew)
      INTEGER_T DIMDEC(fsi)
      INTEGER_T tilelo(SDIM), tilehi(SDIM)
      INTEGER_T fablo(SDIM), fabhi(SDIM)
      INTEGER_T growlo(3), growhi(3)
      INTEGER_T bfact

      REAL_T snew(DIMV(snew),nstate)
      REAL_T fsi(DIMV(fsi),nFSI)

      INTEGER_T i,j,k,dir,ibase
      REAL_T LS

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if (nstate.ne.(SDIM+1)+ &
          nmat*(num_state_material+ngeom_raw)+1) then
       print *,"nstate invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (nFSI.ne.nparts*nFSI_sub) then
       print *,"nFSI invalid"
       stop
      endif
      if (ngrowFSI.ne.3) then
       print *,"ngrowFSI.ne.3"
       stop
      endif
      if ((nparts.lt.1).or.(nparts.gt.nmat)) then
       print *,"nparts invalid FORT_COPY_VEL_ON_SIGN"
       stop
      endif
      if ((partid.lt.0).or.(partid.ge.nparts)) then
       print *,"partid invalid"
       stop
      endif
      if (nFSI_sub.ne.12) then
       print *,"nFSI_sub invalid"
       stop
      endif
      if ((im_part.lt.0).or. &
          (im_part.ge.nmat).or. &
          (is_lag_part(nmat,im_part+1).ne.1)) then
       print *,"im_part invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(snew),1,-1,8)
      call checkbound(fablo,fabhi,DIMS(fsi),ngrowFSI,-1,8)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      ibase=partid*nFSI_sub

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

        ! nmat x (velocity + LS + temperature + flag + stress)
       LS=fsi(D_DECL(i,j,k),ibase+4)
       if (LS.ge.zero) then
        do dir=1,SDIM
         snew(D_DECL(i,j,k),dir)=fsi(D_DECL(i,j,k),ibase+dir)
        enddo
       else if (LS.le.zero) then
        ! do nothing
       else
        print *,"LS bust"
        stop
       endif

      enddo
      enddo
      enddo

      return
      end subroutine FORT_COPY_VEL_ON_SIGN


      subroutine FORT_BUILD_MOMENT( &
       level, &
       finest_level, &
       nFSI, &
       nFSI_sub, &
       nparts, &
       ngrowFSI, &
       im_solid_map, &
       xlo,dx, &
       snew,DIMS(snew), &
       lsnew,DIMS(lsnew), &
       fsi,DIMS(fsi), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       nmat,nstate)
      use probcommon_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      IMPLICIT NONE

      INTEGER_T level 
      INTEGER_T finest_level 
      INTEGER_T nFSI
      INTEGER_T nFSI_sub
      INTEGER_T nparts
      INTEGER_T ngrowFSI
      INTEGER_T im_solid_map(nparts)
      INTEGER_T nmat
      INTEGER_T nstate
      REAL_T xlo(SDIM),dx(SDIM)
      INTEGER_T DIMDEC(snew)
      INTEGER_T DIMDEC(lsnew)
      INTEGER_T DIMDEC(fsi)
      INTEGER_T tilelo(SDIM), tilehi(SDIM)
      INTEGER_T fablo(SDIM), fabhi(SDIM)
      INTEGER_T growlo(3), growhi(3)
      INTEGER_T bfact

      REAL_T snew(DIMV(snew),nstate)
      REAL_T lsnew(DIMV(lsnew),nmat*(1+SDIM))
      REAL_T fsi(DIMV(fsi),nFSI)

      INTEGER_T i,j,k
      INTEGER_T i1,j1,k1
      INTEGER_T im_part
      INTEGER_T partid
      INTEGER_T k1lo,k1hi
      INTEGER_T dir,ibase
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf
      REAL_T ldata(D_DECL(3,3,3))
      REAL_T volume_frac,facearea
      REAL_T areacentroid(SDIM)
      REAL_T centroid(SDIM)
      REAL_T cencell(SDIM)
      REAL_T volcell
      REAL_T LS_center,mag,delta
      REAL_T nrm(SDIM)
      INTEGER_T vofcomp

      nhalf=3

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if (nstate.ne.(SDIM+1)+ &
          nmat*(num_state_material+ngeom_raw)+1) then
       print *,"nstate invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (nFSI_sub.ne.12) then
       print *,"nFSI_sub invalid"
       stop
      endif
      if (nFSI.ne.nparts*nFSI_sub) then
       print *,"nFSI invalid"
       stop
      endif
      if ((nparts.lt.1).or.(nparts.gt.nmat)) then
       print *,"nparts invalid FORT_BUILD_MOMENT"
       stop
      endif
      if (ngrowFSI.ne.3) then
       print *,"ngrowFSI invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid 33"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(lsnew),1,-1,8)
      call checkbound(fablo,fabhi,DIMS(snew),1,-1,8)
      call checkbound(fablo,fabhi,DIMS(fsi),ngrowFSI,-1,8)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      if (SDIM.eq.2) then
       k1lo=0
       k1hi=0
      else if (SDIM.eq.3) then
       k1lo=-1
       k1hi=1
      else
       print *,"dimension bust"
       stop
      endif

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)

       do partid=1,nparts
        im_part=im_solid_map(partid)+1
        if ((im_part.lt.1).or.(im_part.gt.nmat)) then
         print *,"im_part invalid FORT_BUILD_MOMENT"
         stop
        endif
        if ((FSI_flag(im_part).eq.2).or. & ! prescribed solid CAD
            (FSI_flag(im_part).eq.4).or. & ! CTML FSI
            (FSI_flag(im_part).eq.6).or. & ! ice CAD
            (FSI_flag(im_part).eq.7)) then ! fluid CAD

         if (is_lag_part(nmat,im_part).ne.1) then
          print *,"is_lag_part(nmat,im_part).ne.1"
          stop
         endif

         ! nparts x (vel+LS+temperature+flag+stress)
         ibase=(partid-1)*nFSI_sub
         do i1=-1,1
         do j1=-1,1
         do k1=k1lo,k1hi
          ldata(D_DECL(i1+2,j1+2,k1+2))= &
            fsi(D_DECL(i+i1,j+j1,k+k1),ibase+4)
         enddo 
         enddo 
         enddo  ! i1,j1,k1=-1..1

         mag=zero
         do dir=1,SDIM
          i1=0
          j1=0
          k1=0
          if (dir.eq.1) then
           i1=1
          else if (dir.eq.2) then
           j1=1
          else if ((dir.eq.3).and.(SDIM.eq.3)) then
           k1=1
          else
           print *,"dir invalid"
           stop
          endif
          delta=xsten(2,dir)-xsten(-2,dir)
          if (delta.le.zero) then
           print *,"delta invalid delta=",delta
           print *,"dir=",dir
           print *,"xsten(2,dir)=",xsten(2,dir)
           print *,"xsten(-2,dir)=",xsten(-2,dir)
           print *,"i,j,k ",i,j,k
           print *,"xlo=",xlo(1),xlo(2),xlo(SDIM)
           print *,"fablo=",fablo(1),fablo(2),fablo(SDIM)
           print *,"bfact=",bfact
           print *,"nhalf=",nhalf
           print *,"dx=",dx(1),dx(2),dx(SDIM)
           stop
          endif
          nrm(dir)=(ldata(D_DECL(2+i1,2+j1,2+k1))- &
                    ldata(D_DECL(2-i1,2-j1,2-k1)))/delta
          mag=mag+nrm(dir)**2
         enddo ! dir=1..sdim
         mag=sqrt(mag)
         if (mag.eq.zero) then
          nrm(1)=one
         else if (mag.gt.zero) then
          do dir=1,SDIM
           nrm(dir)=nrm(dir)/mag
          enddo
         else
          print *,"mag invalid 6"
          stop
         endif

         do dir=1,SDIM
          lsnew(D_DECL(i,j,k),nmat+(im_part-1)*SDIM+dir)=nrm(dir)
         enddo
 
         LS_center=fsi(D_DECL(i,j,k),ibase+4)

         call getvolume( &
          bfact,dx,xsten,nhalf, &
          ldata,volume_frac,facearea, &
          centroid,areacentroid,VOFTOL,SDIM)
         call CISBOX(xsten,nhalf, &
          xlo,dx,i,j,k, &
          bfact,level, &
          volcell,cencell,SDIM)
         vofcomp=(SDIM+1)+nmat*num_state_material+ &
            (im_part-1)*ngeom_raw

         if (volume_frac.lt.zero) then
          print *,"volume_frac.lt.zero"
          stop
         else if (volume_frac.le.VOFTOL) then
          if (LS_center.ge.-VOFTOL*dx(1)) then
           volume_frac=VOFTOL_SLOPES
          endif
         else if (volume_frac.gt.one) then
          print *,"volume_frac.gt.one"
          stop
         else if (volume_frac.ge.one-VOFTOL) then
          if (LS_center.le.VOFTOL*dx(1)) then
           volume_frac=one-VOFTOL_SLOPES
          endif
         else if ((volume_frac.ge.VOFTOL).and. &
                  (volume_frac.le.one-VOFTOL)) then
          ! do nothing
         else
          print *,"volume_frac invalid"
          stop
         endif
 
         snew(D_DECL(i,j,k),vofcomp+1)=volume_frac 
         do dir=1,SDIM
          snew(D_DECL(i,j,k),vofcomp+1+dir)=centroid(dir)-cencell(dir)
         enddo 
        else if (FSI_flag(im_part).eq.1) then ! prescribed solid EUL
         ! do nothing
        else
         print *,"FSI_flag invalid"
         stop
        endif
       enddo ! partid=1..nparts

      enddo
      enddo
      enddo

      return
      end subroutine FORT_BUILD_MOMENT


        ! called form tensor_advecton_update() in NavierStokes.cpp
        ! vel is the advective velocity
      subroutine FORT_UPDATETENSOR( &
       level, &
       finest_level, &
       nmat, &
       im_critical, &  ! 0<=im_critical<=nmat-1
       ncomp_visc, & 
       vof,DIMS(vof), &
       visc,DIMS(visc), &
       tendata,DIMS(tendata), & ! tendata: FORT_GETSHEAR,iproject=onlyscalar=0
       dx,xlo, &
       vel,DIMS(vel), &
       tnew,DIMS(tnew), &
       told,DIMS(told), &
       xdisplace,DIMS(xdisplace), &
       tilelo, tilehi,  &
       fablo, fabhi, &
       bfact,  &
       dt, &
       elastic_time, &
       viscoelastic_model, &
       polymer_factor,irz,bc, &
       transposegradu)
      use probcommon_module
      use global_utility_module
      use godunov_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: nmat,im_critical
      INTEGER_T, intent(in) :: ncomp_visc
      INTEGER_T, intent(in) :: DIMDEC(vof)
      INTEGER_T, intent(in) :: DIMDEC(visc)
      INTEGER_T, intent(in) :: DIMDEC(tendata)
      INTEGER_T, intent(in) :: DIMDEC(vel)
      INTEGER_T, intent(in) :: DIMDEC(tnew)
      INTEGER_T, intent(in) :: DIMDEC(told)
      INTEGER_T, intent(in) :: DIMDEC(xdisplace)
      INTEGER_T, intent(in) :: tilelo(SDIM), tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM), fabhi(SDIM)
      INTEGER_T :: growlo(3), growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: dx(SDIM),xlo(SDIM)

      REAL_T, intent(in) :: vof(DIMV(vof),nmat*ngeom_recon)
      REAL_T, intent(in) :: visc(DIMV(visc),ncomp_visc)
       ! 1: sqrt(2 * D : D)
       ! 2..2+9-1: D11,D12,D13,D21,D22,D23,D31,D32,D33
       ! 11..11+9-1: ux,uy,uz,vx,vy,vz,wx,wy,wz
      REAL_T, intent(in) :: tendata(DIMV(tendata),20)
      REAL_T, intent(in) :: vel(DIMV(vel),SDIM)
      REAL_T, intent(out) :: tnew(DIMV(tnew),FORT_NUM_TENSOR_TYPE)
      REAL_T, intent(in) :: told(DIMV(told),FORT_NUM_TENSOR_TYPE)
      REAL_T, intent(in) :: xdisplace(DIMV(xdisplace),SDIM)

      INTEGER_T :: i,j,k,n
      REAL_T dt,elastic_time
      INTEGER_T viscoelastic_model
      REAL_T polymer_factor
      INTEGER_T transposegradu
      INTEGER_T bc(SDIM,2,SDIM)
      INTEGER_T irz
      INTEGER_T ii,jj,kk
      INTEGER_T iii,jjj
      REAL_T rsign 
      REAL_T visctensor(3,3)
      REAL_T gradu_FENECR(3,3)
      REAL_T gradV(3,3)
      REAL_T Q(3,3)
      REAL_T W_Jaumann(3,3)  ! W=(1/2)(grad V - (grad V)^T)
      REAL_T Aadvect(3,3)
      REAL_T Smult(3,3)
      REAL_T SA(3,3)
      REAL_T SAS(3,3)
      REAL_T shear
      REAL_T modtime,traceA
      INTEGER_T vofcomp
      REAL_T vfrac
      REAL_T growthrate,rr,uu

      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf

      INTEGER_T LS_or_VOF_flag
      INTEGER_T im_elastic

      nhalf=3

      if (irz.ne.levelrz) then
       print *,"irz invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid60"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid 34"
       stop
      endif

      if (viscoelastic_model.lt.0) then
       print *,"viscoelastic_model invalid"
       stop
      endif

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if ((im_critical.lt.0).or.(im_critical.ge.nmat)) then
       print *,"im_critical invalid27"
       stop
      endif
      if (ncomp_visc.ne.3*nmat) then
       print *,"ncomp visc invalid"
       stop
      endif
      if (dt.le.zero) then
       print *,"dt invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(vof),1,-1,9)
      call checkbound(fablo,fabhi,DIMS(visc),0,-1,9)
      call checkbound(fablo,fabhi,DIMS(tendata),0,-1,9)
      call checkbound(fablo,fabhi,DIMS(vel),1,-1,61)
      call checkbound(fablo,fabhi,DIMS(tnew),0,-1,62)
      call checkbound(fablo,fabhi,DIMS(told),0,-1,63)
      call checkbound(fablo,fabhi,DIMS(xdisplace),1,-1,63)

      if ((transposegradu.ne.0).and.(transposegradu.ne.1)) then
       print *,"transposegradu invalid"
       stop
      endif

      if (viscoelastic_model.eq.2) then ! elastic material

       LS_or_VOF_flag=1 ! =1 => use VOF for upwinding near interfaces
       im_elastic=im_critical+1
        ! elastic bulk modulus not included.
       call local_tensor_from_xdisplace( &
        LS_or_VOF_flag, &
        im_elastic, &
        tilelo,tilehi, &  ! tile box dimensions
        fablo,fabhi, &    ! fortran array box dimensions containing the tile
        bfact, &          ! space order
        level, &          ! 0<=level<=finest_level
        finest_level, &
        xlo,dx, &         ! xlo is lower left hand corner coordinate of fab
        FORT_NUM_TENSOR_TYPE, & !ncomp_tensor=4 in 2D (11,12,22,33) and 6 in 3D 
        nmat, &
        vof, & ! LS placeholder
        DIMS(vof), & ! LS placeholder
        vof, &  
        DIMS(vof), &
        tnew, &       ! FAB that holds elastic tensor, Q, when complete
        DIMS(tnew), &
        xdisplace, &      
        DIMS(xdisplace)) 

      else if ((viscoelastic_model.eq.0).or. &
               (viscoelastic_model.eq.1)) then
       ! do nothing
      else if (viscoelastic_model.eq.3) then ! incremental model
       ! do nothing
      else
       print *,"viscoelastic_model invalid"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridsten_level(xsten,i,j,k,level,nhalf)

       ! tendata has: |D|, D, grad U
       ! 1: sqrt(2 * D : D)
       ! 2..2+9-1: D11,D12,D13,D21,D22,D23,D31,D32,D33
       ! 11..11+9-1: ux,uy,uz,vx,vy,vz,wx,wy,wz
       shear=tendata(D_DECL(i,j,k),1)
       n=2
       do ii=1,3
       do jj=1,3
         ! (1/2) (grad U + (grad U)^T)
        visctensor(ii,jj)=tendata(D_DECL(i,j,k),n)
        n=n+1
       enddo 
       enddo 
       do ii=1,3
       do jj=1,3
        gradV(ii,jj)=tendata(D_DECL(i,j,k),n) !(veldir,dir)
        if (transposegradu.eq.0) then
          !gradu(veldir,dir)
         gradu_FENECR(ii,jj)=tendata(D_DECL(i,j,k),n)
        else if (transposegradu.eq.1) then
          !gradu(dir,veldir)
         gradu_FENECR(jj,ii)=tendata(D_DECL(i,j,k),n)
        else
         print *,"transposegradu invalid"
         stop
        endif

        n=n+1
       enddo 
       enddo 
       do ii=1,3
       do jj=1,3
        W_Jaumann(ii,jj)=half*(gradV(ii,jj)-gradV(jj,ii))
       enddo 
       enddo 
       
       do ii=1,3
       do jj=1,3
         Q(ii,jj)=zero
       enddo
       enddo

       Q(1,1)=told(D_DECL(i,j,k),1)
       Q(1,2)=told(D_DECL(i,j,k),2)
       Q(2,2)=told(D_DECL(i,j,k),3)
       Q(3,3)=told(D_DECL(i,j,k),4)
#if (AMREX_SPACEDIM==3)
       Q(1,3)=told(D_DECL(i,j,k),5)
       Q(2,3)=told(D_DECL(i,j,k),6)
#endif
       Q(2,1)=Q(1,2)
       Q(3,1)=Q(1,3)
       Q(3,2)=Q(2,3)

        ! viscoelastic_model==0 => modtime<=elastic_time
        ! viscoelastic_model==2 => modtime==elastic_time
        ! viscoelastic_model==3 => modtime==elastic_time
       modtime=visc(D_DECL(i,j,k),2*nmat+im_critical+1)
       vofcomp=im_critical*ngeom_recon+1
       vfrac=vof(D_DECL(i,j,k),vofcomp)
       if (modtime.lt.zero) then
        print *,"modtime invalid"
        stop
       endif

       if (abs(vfrac).le.VOFTOL) then

        do ii=1,3
        do jj=1,3
         Q(ii,jj)=zero
        enddo
        enddo

       else if ((vfrac.gt.zero).and.(vfrac.le.one+VOFTOL)) then

        if (viscoelastic_model.eq.2) then ! elastic material

         do ii=1,3 
         do jj=1,3 
          Q(ii,jj)=zero
         enddo
         enddo
         Q(1,1)=tnew(D_DECL(i,j,k),1)
         Q(1,2)=tnew(D_DECL(i,j,k),2)
         Q(2,2)=tnew(D_DECL(i,j,k),3)
         Q(3,3)=tnew(D_DECL(i,j,k),4)
#if (AMREX_SPACEDIM==3)
         Q(1,3)=tnew(D_DECL(i,j,k),5)
         Q(2,3)=tnew(D_DECL(i,j,k),6)
#endif
         Q(2,1)=Q(1,2)
         Q(3,1)=Q(1,3)
         Q(3,2)=Q(2,3)

        else if ((viscoelastic_model.eq.0).or. &
                 (viscoelastic_model.eq.1).or. &
                 (viscoelastic_model.eq.3)) then

         do ii=1,3 
          do jj=1,3 
           Aadvect(ii,jj)=Q(ii,jj)
            !cfl cond: |u|dt<dx and dt|gradu|<1
           if ((viscoelastic_model.eq.0).or. &
               (viscoelastic_model.eq.1)) then
            Smult(ii,jj)=dt*gradu_FENECR(ii,jj) 
           else if (viscoelastic_model.eq.3) then
            Smult(ii,jj)=dt*W_Jaumann(ii,jj) 
           else
            print *,"viscoelastic_model invalid"
            stop
           endif

           if (Smult(ii,jj).le.-one+VOFTOL) then
            Smult(ii,jj)=-one+VOFTOL
           else if (Smult(ii,jj).ge.one-VOFTOL) then
            Smult(ii,jj)=one-VOFTOL
           else if (abs(Smult(ii,jj)).le.one) then
            ! do nothing
           else
            print *,"Smult(ii,jj) became corrupt"
            stop
           endif
    
          enddo ! jj=1,3
          Smult(ii,ii)=Smult(ii,ii)+one
          Aadvect(ii,ii)=Aadvect(ii,ii)+one
         enddo  ! ii=1,3

         if ((viscoelastic_model.eq.0).or. &
             (viscoelastic_model.eq.1)) then

          do ii=1,3 
           if (Aadvect(ii,ii).lt.zero) then
            Aadvect(ii,ii)=zero
            print *,"WARNING Q^advect+I no longer positive definite"
            print *,"viscoelastic_model=",viscoelastic_model
           endif
          enddo  ! ii=1,3

          if (SDIM.eq.3) then

           rsign=Aadvect(1,3)
           if (Aadvect(1,3)**2.gt.Aadvect(1,1)*Aadvect(3,3)) then
            Aadvect(1,3)=sqrt(Aadvect(1,1)*Aadvect(3,3))
            if (rsign.lt.zero) then
             Aadvect(1,3)=-Aadvect(1,3)
            endif
            Aadvect(3,1)=Aadvect(1,3)
           endif
           rsign=Aadvect(2,3)
           if (Aadvect(2,3)**2.gt.Aadvect(2,2)*Aadvect(3,3)) then
            Aadvect(2,3)=sqrt(Aadvect(2,2)*Aadvect(3,3))
            if (rsign.lt.zero) then
             Aadvect(2,3)=-Aadvect(2,3)
            endif
            Aadvect(3,2)=Aadvect(2,3)
           endif

          else if (SDIM.eq.2) then
           ! do nothing
          else
           print *,"dimension bust"
           stop
          endif

          rsign=Aadvect(1,2)
          if (Aadvect(1,2)**2.gt.Aadvect(1,1)*Aadvect(2,2)) then
           Aadvect(1,2)=sqrt(Aadvect(1,1)*Aadvect(2,2))
           if (rsign.lt.zero) then
            Aadvect(1,2)=-Aadvect(1,2)
           endif
           Aadvect(2,1)=Aadvect(1,2)
          endif
 
         else if (viscoelastic_model.eq.3) then
          ! do nothing
         else
          print *,"viscoelastic_model invalid"
          stop
         endif

         if ((viscoelastic_model.eq.0).or. &
             (viscoelastic_model.eq.1)) then
          ! do nothing
         else if (viscoelastic_model.eq.3) then
          do ii=1,3
          do jj=1,3
           Aadvect(ii,jj)=Aadvect(ii,jj)+dt*two*visctensor(ii,jj) 
          enddo
          enddo
         else
          print *,"viscoelastic_model invalid"
          stop
         endif

         do ii=1,3
         do jj=1,3
          SA(ii,jj)=zero
          do kk=1,3
           SA(ii,jj)=SA(ii,jj)+Smult(ii,kk)*Aadvect(kk,jj)
          enddo
         enddo
         enddo
        
         do ii=1,3
         do jj=1,3
          SAS(ii,jj)=zero
          do kk=1,3
           SAS(ii,jj)=SAS(ii,jj)+SA(ii,kk)*Smult(jj,kk)
          enddo
         enddo
         enddo

         do ii=1,3
         do jj=1,3
          Q(ii,jj)=SAS(ii,jj)

          if ((ii.eq.3).and.(jj.eq.3)) then

           if (SDIM.eq.3) then
            ! do nothing
           else if (SDIM.eq.2) then
            if (levelrz.eq.0) then
             ! do nothing
            else if (levelrz.eq.1) then
             ! GETSHEAR put u/r in gradu(3,3)
             ! for hoop stress term: 
             ! dA/dt=2 u A/r
             ! if u<0:
             ! A^n+1 - A^n = dt 2uA^{n+1}/r
             ! (1-dt 2u/r)A^n+1=A^n
             ! A^n+1=A^n/(1-dt 2u/r)
             ! if u>0:
             ! A^n+1=(1+2u dt/r)A^n
             rr=xsten(0,1)
             if (rr.le.zero) then
              print *,"rr invalid"
              stop
             endif
             uu=gradu_FENECR(3,3)*rr
             growthrate=two*uu*dt/rr
             if (uu.gt.zero) then
              Q(ii,jj)=(one+growthrate)*Aadvect(ii,jj)
             else if (uu.lt.zero) then   
              Q(ii,jj)=Aadvect(ii,jj)/(one-growthrate)
             else if (uu.eq.zero) then
              Q(ii,jj)=Aadvect(ii,jj)
             else
              print *,"uu bust"
              stop
             endif
              ! Q=S A S^T at this stage
             if ((viscoelastic_model.eq.0).or. &
                 (viscoelastic_model.eq.1)) then
              if (Q(ii,jj).gt.zero) then
               ! do nothing
              else if (Q(ii,jj).le.zero) then
               print *,"Q(ii,jj)<=0"
               print *,"viscoelastic_model=",viscoelastic_model
               stop
              else
               print *,"Q(ii,jj) bust"
               stop
              endif
             else if (viscoelastic_model.eq.3) then
              ! do nothing
             else
              print *,"viscoelastic_model invalid"
              stop
             endif
            else if (levelrz.eq.3) then
             ! do nothing
            else
             print *,"levelrz invalid"
             stop
            endif 
           else
            print *,"dimension bust"
            stop
           endif
          else if ((ii.ne.3).or.(jj.ne.3)) then
           ! do nothing
          else
           print *,"ii or jj invalid"
           stop
          endif

         enddo  ! jj=1..3
         enddo  ! ii=1..3

         if ((viscoelastic_model.eq.0).or. &
             (viscoelastic_model.eq.1)) then
          do ii=1,3
           if (Q(ii,ii).lt.zero) then
            Q(ii,ii)=zero
            print *,"WARNING Q+I no longer positive definite"
            print *,"viscoelastic_model=",viscoelastic_model
            do iii=1,3
            do jjj=1,3
             print *,"iii,jjj,Q,Aadvect,Smult ",iii,jjj,Q(iii,jjj), &
              Aadvect(iii,jjj),Smult(iii,jjj)
            enddo
            enddo
            print *,"i,j,k,lo,hi ",i,j,k, &
             growlo(1),growlo(2),growlo(SDIM), &
             growhi(1),growhi(2),growhi(SDIM)
           else if (Q(ii,ii).ge.zero) then
            ! do nothing
           else
            print *,"Q invalid"
            stop
           endif

          enddo ! ii=1..3

          if (SDIM.eq.3) then

           rsign=Q(1,3)
           if (Q(1,3)**2.gt.Q(1,1)*Q(3,3)) then
            Q(1,3)=sqrt(Q(1,1)*Q(3,3))
            if (rsign.lt.zero) then
             Q(1,3)=-Q(1,3)
            endif
            Q(3,1)=Q(1,3)
           endif
           rsign=Q(2,3)
           if (Q(2,3)**2.gt.Q(2,2)*Q(3,3)) then
            Q(2,3)=sqrt(Q(2,2)*Q(3,3))
            if (rsign.lt.zero) then
             Q(2,3)=-Q(2,3)
            endif
            Q(3,2)=Q(2,3)
           endif

          else if (SDIM.eq.2) then
           ! do nothing
          else
           print *,"dimension bust"
           stop
          endif

          rsign=Q(1,2)
          if (Q(1,2)**2.gt.Q(1,1)*Q(2,2)) then
           Q(1,2)=sqrt(Q(1,1)*Q(2,2))
           if (rsign.lt.zero) then
            Q(1,2)=-Q(1,2)
           endif
           Q(2,1)=Q(1,2)
          endif
         else if (viscoelastic_model.eq.3) then
          ! do nothing
         else
          print *,"viscoelastic_model"
          stop
         endif 
 
         do ii=1,3
          Q(ii,ii)=Q(ii,ii)-one
         enddo

         ! note: for viscoelastic_model==2 or
         !           viscoelastic_model==3,
         !  modtime=lambda=elastic_time >> 1
         !
         ! Q^n+1=lambda Q^n/(lambda+dt)
         ! lambda (Q^n+1-Q^n)=-dt Q^n+1
         ! (Q^n+1-Q^n)/dt = -Q^n+1/lambda
         ! Q^n+1/lambda=Q^n/(lambda+dt) 
         do ii=1,3
         do jj=1,3
          Q(ii,jj)=modtime*Q(ii,jj)/(modtime+dt)
         enddo
         enddo

         traceA=zero
         do ii=1,3
          traceA=traceA+Q(ii,ii)+one
         enddo
         if ((viscoelastic_model.eq.0).or. &
             (viscoelastic_model.eq.1)) then
          if (traceA.lt.zero) then
           print *,"traceA cannot be negative!"
           print *,"viscoelastic_model=",viscoelastic_model
           stop
          endif
         else if (viscoelastic_model.eq.3) then
          !check nothing
         else
          print *,"viscoelastic_model invalid"
          stop
         endif

        else 
         print *,"viscoelastic_model invalid"
         stop
        endif

       else
        print *,"vfrac invalid"
        stop
       endif
 
       tnew(D_DECL(i,j,k),1)=Q(1,1)
       tnew(D_DECL(i,j,k),2)=Q(1,2)
       tnew(D_DECL(i,j,k),3)=Q(2,2)
       tnew(D_DECL(i,j,k),4)=Q(3,3)
#if (AMREX_SPACEDIM==3)
       tnew(D_DECL(i,j,k),5)=Q(1,3)
       tnew(D_DECL(i,j,k),6)=Q(2,3)
#endif

      enddo
      enddo
      enddo

      return
      end subroutine FORT_UPDATETENSOR

      subroutine FORT_FIX_HOOP_TENSOR( &
       level, &
       finest_level, &
       nmat,im, & 
       vof,DIMS(vof), &
       dx,xlo, &
       tnew,DIMS(tnew), &
       tilelo, tilehi,  &
       fablo, fabhi, &
       bfact,  &
       irz)
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T level
      INTEGER_T finest_level
      INTEGER_T nmat,im
      INTEGER_T i,j,k
      INTEGER_T DIMDEC(vof)
      INTEGER_T DIMDEC(tnew)
      INTEGER_T tilelo(SDIM), tilehi(SDIM)
      INTEGER_T fablo(SDIM), fabhi(SDIM)
      INTEGER_T growlo(3), growhi(3)
      INTEGER_T bfact
      REAL_T dx(SDIM),xlo(SDIM)

      REAL_T vof(DIMV(vof),nmat*ngeom_recon)
      REAL_T tnew(DIMV(tnew),FORT_NUM_TENSOR_TYPE)

      INTEGER_T irz
      INTEGER_T vofcomp
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf
      REAL_T Q,Qp1,Qnew
      REAL_T vfrac,vfracp1
      REAL_T rr,rrp1

      nhalf=3

      if (irz.ne.levelrz) then
       print *,"irz invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid61"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid 35"
       stop
      endif

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if ((im.lt.0).or.(im.ge.nmat)) then
       print *,"im invalid28"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(vof),0,-1,9)
      call checkbound(fablo,fabhi,DIMS(tnew),0,-1,62)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)

      if ((SDIM.eq.2).and.(levelrz.eq.1)) then

       if (growlo(1).eq.0) then
        i=0
        do j=growlo(2),growhi(2)
        do k=growlo(3),growhi(3)

         call gridsten_level(xsten,i,j,k,level,nhalf)

         if (abs(xsten(-1,1)).le.VOFTOL) then

          if (fabhi(1).ge.1) then

           Q=tnew(D_DECL(i,j,k),4)
           Qp1=tnew(D_DECL(i+1,j,k),4)
            
           vofcomp=im*ngeom_recon+1
           vfrac=vof(D_DECL(i,j,k),vofcomp)
           vfracp1=vof(D_DECL(i+1,j,k),vofcomp)

           if ((vfrac.ge.-VOFTOL).and. &
               (vfrac.le.one+VOFTOL).and. &
               (vfracp1.ge.-VOFTOL).and. &
               (vfracp1.le.one+VOFTOL)) then

            if ((vfrac.lt.half).or.(vfracp1.lt.half)) then
             Qnew=zero
            else if ((vfrac.ge.half).and.(vfracp1.ge.half)) then
             rr=xsten(0,1)
             rrp1=xsten(2,1)
             if ((rr.gt.zero).and.(rrp1.gt.rr)) then
              ! (Qnew-0)/rr = (Qp1-Qnew)/rrp1
              ! Qnew(1/rr+1/rrp1)=Qp1/rrp1
              ! Qnew(rrp1/rr + 1)=Qp1
              ! Qnew=Qp1/(rrp1/rr + 1)
              Qnew=Qp1/(rrp1/rr+one)
              if (abs(Qnew).gt.abs(Q)) then
               Qnew=Q
              endif
             else
              print *,"rr or rrp1 invalid"
              stop
             endif
            else
             print *,"vfrac invalid"
             stop
            endif
             ! we do not want the negative reflected
             ! ghost value of Q to be less than or equal
             ! to -1.0, otherwise the code will complain about loss of 
             ! the positive definite property.  (the code might still
             ! complain, because "average down" and advection might bring
             ! large values of Q33 to the i=0 column prior to this routine)
            if (Qnew.ge.half) then
             Qnew=half
            endif

            tnew(D_DECL(i,j,k),4)=Qnew
           else
            print *,"vfrac or vfracp1 invalid"
            stop
           endif 

          else
           print *,"fabhi invalid"
           stop
          endif

         else
          print *,"abs(xsten(-1,1)) invalid"
          stop
         endif

        enddo
        enddo
 
       else if (growlo(1).gt.0) then
        ! do nothing
       else
        print *,"growlo(1) invalid"
        stop
       endif
      else
       print *,"sdim or levelrz invalid"
       stop
      endif

      return
      end subroutine FORT_FIX_HOOP_TENSOR

        ! u_max(1..sdim) is max vel in dir.
        ! u_max(sdim+1) is max c^2
      subroutine FORT_ESTDT ( &
        enable_spectral, &
        AMR_min_phase_change_rate, &
        AMR_max_phase_change_rate, &
        elastic_time, &
        microlayer_substrate, &
        microlayer_angle, &
        microlayer_size, &
        macrolayer_size, &
        latent_heat, &
        reaction_rate, &
        freezing_model, &
        Tanasawa_or_Schrage_or_Kassemi, &
        distribute_from_target, &
        saturation_temp, &
        mass_fraction_id, &
        molar_mass, &
        species_molar_mass, &
        velmac,DIMS(velmac), &
        velcell,DIMS(velcell), &
        solidfab,DIMS(solidfab), &
        den,DIMS(den), &
        vof,DIMS(vof), &
        dist,DIMS(dist), &
        xlo,dx, &
        tilelo,tilehi, &
        fablo,fabhi, &
        bfact, &
        min_stefan_velocity_for_dt, &
        cap_wave_speed, &
        u_max, &
        u_max_estdt, &
        u_max_cap_wave, &
        dt_min, &
        rzflag, &
        Uref,Lref, &
        nten, &
        use_lsa, &
        denconst, &
        denconst_gravity, &
        visc_coef, &
        ns_gravity, &
        terminal_velocity_dt, &
        dirnormal, &
        nmat, &
        nparts, &
        nparts_def, &
        im_solid_map, &
        material_type, &
        time, &
        shock_timestep, &
        cfl, &
        EILE_flag, &
        level,finest_level)
      use probf90_module
      use global_utility_module
      use MOF_routines_module
      use hydrateReactor_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: nparts
      INTEGER_T, intent(in) :: nparts_def
      INTEGER_T, intent(in) :: im_solid_map(nparts_def)
      INTEGER_T, intent(in) :: enable_spectral
      INTEGER_T, intent(in) :: use_lsa
      INTEGER_T, intent(in) :: level,finest_level
      REAL_T, intent(in) :: cfl
      INTEGER_T, intent(in) :: EILE_flag
      INTEGER_T, intent(in) :: nmat,nten
      REAL_T, intent(in) :: AMR_min_phase_change_rate(SDIM)
      REAL_T, intent(in) :: AMR_max_phase_change_rate(SDIM)
      REAL_T, intent(in) :: elastic_time(nmat)
      INTEGER_T, intent(in) :: shock_timestep(nmat)
      INTEGER_T, intent(in) :: material_type(nmat)
      INTEGER_T, intent(in) :: microlayer_substrate(nmat)
      REAL_T, intent(in) :: microlayer_angle(nmat)
      REAL_T, intent(in) :: microlayer_size(nmat)
      REAL_T, intent(in) :: macrolayer_size(nmat)
      REAL_T, intent(in) :: latent_heat(2*nten)
      REAL_T, intent(in) :: reaction_rate(2*nten)
      REAL_T :: K_f
      INTEGER_T, intent(in) :: freezing_model(2*nten)
      INTEGER_T, intent(in) :: Tanasawa_or_Schrage_or_Kassemi(2*nten)
      INTEGER_T, intent(in) :: distribute_from_target(2*nten)
      REAL_T, intent(in) :: saturation_temp(2*nten)
      INTEGER_T, intent(in) :: mass_fraction_id(2*nten)
      REAL_T, intent(in) :: molar_mass(nmat)
      REAL_T, intent(in) :: species_molar_mass(num_species_var+1)
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      REAL_T, intent(in) :: time
      REAL_T u_core,u_core_estdt,uu,uu_estdt,c_core
      REAL_T cc,cleft,cright
      REAL_T cc_diag,cleft_diag,cright_diag
      INTEGER_T i,j,k
      INTEGER_T icell,jcell,kcell
      INTEGER_T ialt,jalt,kalt
      INTEGER_T, intent(in) :: rzflag
      INTEGER_T, intent(in) :: dirnormal
      INTEGER_T side,dir2
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(inout) :: u_max(SDIM+1)
      REAL_T, intent(inout) :: u_max_estdt(SDIM+1)
      REAL_T, intent(inout) :: u_max_cap_wave
      REAL_T, intent(inout) :: dt_min
      REAL_T user_tension(nten)
      REAL_T Uref,Lref
      REAL_T, intent(in) :: denconst(nmat)
      REAL_T, intent(in) :: denconst_gravity(nmat)
      REAL_T, intent(in) :: visc_coef
      REAL_T, intent(in) :: ns_gravity
      INTEGER_T, intent(inout) :: terminal_velocity_dt
      INTEGER_T, intent(in) :: DIMDEC(velmac)
      INTEGER_T, intent(in) :: DIMDEC(velcell)
      INTEGER_T, intent(in) :: DIMDEC(vof)
      INTEGER_T, intent(in) :: DIMDEC(dist)
      INTEGER_T, intent(in) :: DIMDEC(solidfab)
      INTEGER_T, intent(in) :: DIMDEC(den)
      REAL_T, intent(in) :: velmac(DIMV(velmac))
      REAL_T, intent(in) :: velcell(DIMV(velcell),SDIM)
      REAL_T, intent(in) :: solidfab(DIMV(solidfab),nparts_def*SDIM) 
       ! den,denA,E,temp
      REAL_T, intent(in) :: den(DIMV(den),num_state_material*nmat)  
      REAL_T, intent(in) :: vof(DIMV(vof),nmat*ngeom_raw)
      REAL_T, intent(in) :: dist(DIMV(dist),nmat)
      REAL_T, intent(in) :: min_stefan_velocity_for_dt
      REAL_T, intent(inout) :: cap_wave_speed(nten)
      REAL_T hx,hxmac
      REAL_T dthold
      INTEGER_T ii,jj,kk
      INTEGER_T im,im_primaryL,im_primaryR
      INTEGER_T ibase
      REAL_T xstenMAC(-3:3,SDIM)
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf
      REAL_T LSleft(nmat)
      REAL_T LSright(nmat)
      INTEGER_T im_opp
      INTEGER_T iten
      INTEGER_T im_source,im_dest
      REAL_T temperature_left,temperature_right
      REAL_T density_left,density_right
      REAL_T internal_energy_left,internal_energy_right
      REAL_T massfrac_parm_left(num_species_var+1)
      REAL_T massfrac_parm_right(num_species_var+1)
      REAL_T gradh
      INTEGER_T nten_test
      REAL_T weymouth_factor,weymouth_cfl
      REAL_T dxmin,dxmax,dxmaxLS,den1,den2,visc1,visc2
      INTEGER_T recompute_wave_speed
      REAL_T uulocal
      REAL_T smallestL
      REAL_T denjump
      REAL_T denjump_gravity
      REAL_T denjump_terminal
      REAL_T denjump_temp
      REAL_T denmax
      REAL_T denmax_gravity
      REAL_T USTEFAN,USTEFAN_hold
      REAL_T LS1,LS2,Tsrc,Tdst,Dsrc,Ddst,Csrc,Cdst,delta
      REAL_T VOFsrc,VOFdst
      REAL_T LL
      INTEGER_T velcomp
      INTEGER_T dcompsrc,dcompdst
      INTEGER_T tcompsrc,tcompdst
      INTEGER_T ireverse
      INTEGER_T ifaceR,jfaceR,kfaceR
      REAL_T uleft,uright
      REAL_T C_w0,PHYDWATER,Cmethane_in_hydrate
      INTEGER_T local_freezing_model
      INTEGER_T local_Tanasawa_or_Schrage_or_Kassemi
      INTEGER_T distribute_from_targ
      INTEGER_T vofcompsrc,vofcompdst
      REAL_T TSAT,Tsrcalt,Tdstalt
      REAL_T uleftcell,urightcell,udiffcell,umaxcell
      REAL_T velsum
      REAL_T RR
      REAL_T level_cap_wave_speed(nten)
      REAL_T ksource,kdest,alpha,beta,dt_heat
      INTEGER_T for_estdt
      REAL_T xI(SDIM)
      REAL_T mu
      INTEGER_T partid
      INTEGER_T ispec
      REAL_T vapor_den
      REAL_T elastic_wave_speed
      REAL_T source_perim_factor
      REAL_T dest_perim_factor
      REAL_T v_terminal
      REAL_T effective_velocity
      REAL_T local_elastic_time
      REAL_T ugrav

      nhalf=3

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if ((nparts.lt.0).or.(nparts.gt.nmat)) then
       print *,"nparts invalid FORT_ESTDT"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.nmat)) then
       print *,"nparts_def invalid FORT_ESTDT"
       stop
      endif
      if ((enable_spectral.lt.0).or. &
          (enable_spectral.gt.3)) then
       print *,"enable_spectral invalid estdt"
       stop
      endif

      if ((terminal_velocity_dt.eq.0).or. &
          (terminal_velocity_dt.eq.1)) then
       ! do nothing
      else
       print *,"terminal_velocity_dt invalid"
       stop
      endif

      vapor_den=one

      denjump=zero
      denjump_gravity=zero
      denmax=zero
      denmax_gravity=zero

       ! dxmin=min_d min_i dxsub_{gridtype,d,i} d=1..sdim  i=0..bfact-1
       ! gridtype=MAC or CELL
       ! if cylindrical coordinates, then dx_{\theta}*=problox
      call get_dxmin(dx,bfact,dxmin)
      if (dxmin.gt.zero) then
       ! do nothing
      else
       print *,"dxmin must be positive"
       stop
      endif

      call get_dxmax(dx,bfact,dxmax)
      call get_dxmaxLS(dx,bfact,dxmaxLS)

      if (Lref.eq.zero) then
       smallestL=dxmin
      else if (Lref.lt.dxmin) then
       smallestL=Lref
      else
       smallestL=dxmin
      endif

      if (visc_coef.ge.zero) then
       ! do nothing
      else
       print *,"visc_coef invalid"
       stop
      endif

      if ((finest_level.lt.0).or. &
          (level.lt.0).or. &
          (level.gt.finest_level)) then
       print *,"level or finest level invalid estdt"
       stop
      endif

      if ((use_lsa.ne.0).and.(use_lsa.ne.1)) then
       print *,"use_lsa invalid"
       stop
      endif

      if (Uref.lt.zero) then
       print *,"Uref invalid"
       stop
      endif
      if (Lref.lt.zero) then
       print *,"Lref invalid"
       stop
      endif

      if ((EILE_flag.eq.-1).or. & ! Weymouth and Yue
          (EILE_flag.eq.1).or.  & ! EILE
          (EILE_flag.eq.2).or.  & ! always EI
          (EILE_flag.eq.3)) then  ! always LE
       ! do nothing
      else 
       print *,"EILE flag invalid"
       stop
      endif

      if (cfl.le.zero) then
       print *,"cfl invalid"
       stop
      endif

      if ((EILE_flag.eq.1).or. & ! EI-LE
          (EILE_flag.eq.2).or. & ! always EI
          (EILE_flag.eq.3)) then ! always LE
       weymouth_cfl=half  ! we advect half cells.
      else if (EILE_flag.eq.-1) then ! Weymouth and Yue
       weymouth_cfl=one/(two*SDIM)
      else
       print *,"EILE_flag invalid"
       stop
      endif
      weymouth_factor=max(weymouth_cfl,weymouth_cfl/cfl)

      do im=1,nmat 

       if ((shock_timestep(im).ne.0).and. &
           (shock_timestep(im).ne.1).and. &
           (shock_timestep(im).ne.2)) then
        print *,"shock_timestep invalid"
        stop
       endif
       if (denconst(im).le.zero) then
        print *,"denconst invalid"
        stop
       endif
       if (denconst_gravity(im).lt.zero) then
        print *,"denconst_gravity invalid"
        stop
       endif
       mu=get_user_viscconst(im,fort_denconst(im),fort_tempconst(im))
       if (mu.lt.zero) then
        print *,"viscconst invalid"
        stop
       endif

       if (is_rigid(nmat,im).eq.0) then
        if (denconst(im).gt.denmax) then
         denmax=denconst(im)
        endif
        if (denconst_gravity(im).gt.denmax_gravity) then
         denmax_gravity=denconst_gravity(im)
        endif
       else if (is_rigid(nmat,im).eq.1) then
        ! do nothing
       else
        print *,"is_rigid invalid"
        stop
       endif

      enddo  ! im=1..nmat

      call get_max_denjump(denjump,nmat)

      denjump_gravity=zero
      do im=1,nmat 
       if (is_rigid(nmat,im).eq.0) then
        do im_opp=im+1,nmat
         if (is_rigid(nmat,im_opp).eq.0) then
          denjump_temp=abs(denconst_gravity(im)-denconst_gravity(im_opp))
          if (denjump_temp.gt.denjump_gravity) then
           denjump_gravity=denjump_temp
          endif
         else if (is_rigid(nmat,im_opp).eq.1) then
          ! do nothing
         else
          print *,"is_rigid invalid"
          stop
         endif
        enddo ! im_opp
       else if (is_rigid(nmat,im).eq.1) then
        ! do nothing
       else
        print *,"is_rigid invalid"
        stop
       endif
      enddo ! im=1..nmat

      call get_max_user_tension(fort_tension,user_tension,nmat,nten)

         ! finest_level is first level tried.
      recompute_wave_speed=0
      if (level.eq.finest_level) then
       recompute_wave_speed=1
      endif
      
      if (recompute_wave_speed.eq.1) then

       do im=1,nmat-1
        do im_opp=im+1,nmat
         if ((im.gt.nmat).or.(im_opp.gt.nmat)) then
          print *,"im or im_opp bust 6"
          stop
         endif
         call get_iten(im,im_opp,iten,nmat)
         if ((is_rigid(nmat,im).eq.1).or. &
             (is_rigid(nmat,im_opp).eq.1)) then
          cap_wave_speed(iten)=zero
         else if (user_tension(iten).eq.zero) then
          cap_wave_speed(iten)=zero
         else if (user_tension(iten).gt.zero) then
          den1=denconst(im)
          den2=denconst(im_opp)
          mu=get_user_viscconst(im,den1,fort_tempconst(im))
          visc1=visc_coef*mu+1.0D-10
          mu=get_user_viscconst(im_opp,den2,fort_tempconst(im_opp))
          visc2=visc_coef*mu+1.0D-10
            ! typically smallestL=dxmin
          call capillary_wave_speed(smallestL,den1,den2,visc1,visc2, &
            user_tension(iten),cap_wave_speed(iten),use_lsa)
         else
          print *,"user_tension invalid"
          stop
         endif
         level_cap_wave_speed(iten)=cap_wave_speed(iten)
        enddo ! im_opp=im+1..nmat
       enddo ! im=1..nmat-1
      else if (recompute_wave_speed.eq.0) then
       do im=1,nmat-1
        do im_opp=im+1,nmat
         call get_iten(im,im_opp,iten,nmat)
         level_cap_wave_speed(iten)=zero
        enddo
       enddo
      else
       print *,"recompute wave speed invalid"
       print *,"use_lsa=",use_lsa
       print *,"Uref=",Uref
       print *,"Lref=",Lref
       print *,"nmat=",nmat
       print *,"level=",level
       print *,"finest_level=",finest_level
       print *,"smallestL= ",smallestL
       print *,"dxmin= ",dxmin
       print *,"cfl= ",cfl
       print *,"denconst(1)= ",denconst(1)
       print *,"get_user_viscconst(1)= ", &
        get_user_viscconst(1,fort_denconst(1),fort_tempconst(1))
       print *,"visc_coef= ",visc_coef
       print *,"user_tension(1)= ",user_tension(1)
       print *,"EILE_flag= ",EILE_flag
       print *,"nten=",nten
       print *,"dirnormal=",dirnormal
       print *,"rzflag=",rzflag
       stop
      endif

      if (rzflag.eq.0) then
       ! do nothing
      else if (rzflag.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (rzflag.eq.3) then
       ! do nothing
      else
       print *,"rzflag invalid"
       stop
      endif

      if ((dirnormal.lt.0).or.(dirnormal.ge.SDIM)) then
       print *,"dirnormal invalid estdt"
       stop
      endif
      if ((adv_dir.lt.1).or.(adv_dir.gt.2*SDIM+1)) then
       print *,"adv_dir invalid ESTDT"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid estdt nten, nten_test ",nten,nten_test
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(velmac),0,dirnormal,4)
      call checkbound(fablo,fabhi,DIMS(velcell),1,-1,4)
      call checkbound(fablo,fabhi,DIMS(solidfab),0,dirnormal,4)
      call checkbound(fablo,fabhi,DIMS(den),1,-1,4)
      call checkbound(fablo,fabhi,DIMS(vof),1,-1,4)
       ! need enough ghost cells for the calls to derive_dist.
      call checkbound(fablo,fabhi,DIMS(dist),2,-1,4)

      if (rzflag.ne.levelrz) then
       print *,"rzflag invalid"
       stop
      endif

      ii=0
      jj=0
      kk=0
      if (dirnormal.eq.0) then
        ii=1
      else if (dirnormal.eq.1) then
        jj=1
      else if ((dirnormal.eq.2).and.(SDIM.eq.3)) then
        kk=1
      else
       print *,"dirnormal invalid estdt 2"
       stop
      endif

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi, &
              0,dirnormal,27)

      if (1.eq.0) then
       print *,"dt_min before estdt loop: ",dt_min
       print *,"weymouth_factor= ",weymouth_factor
      endif

      u_core = zero
      u_core_estdt = zero
      c_core = zero  ! max of c^2

! if rz and dirnormal=0 and u>0, need u dt r/(r-u dt) < dx
! u dt r < dx(r-u dt)
! u dt dx + udt r < dx r
! u dt (dx/r+1) < dx
!
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,dirnormal,31)
       hx=xstenMAC(1,dirnormal+1)-xstenMAC(-1,dirnormal+1)

       RR=one
       if (levelrz.eq.0) then
        ! do nothing
       else if (levelrz.eq.1) then
        ! do nothing
       else if (levelrz.eq.3) then
        if (dirnormal.eq.1) then ! theta direction
         RR=xstenMAC(0,1)
        endif
       else
        print *,"levelrz invalid"
        stop
       endif
       if (RR.le.zero) then
        print *,"RR invalid"
        stop
       endif
 
       hx=hx*RR
       hxmac=hx

       if (hx.le.(one-VOFTOL)*dxmin) then
        print *,"xstenMAC invalid estdt"
        print *,"hx= ",hx
        print *,"hx= ",dxmin
        print *,"i,j,k ",i,j,k
        print *,"dirnormal=",dirnormal
        print *,"xright ",xstenMAC(1,dirnormal+1)
        print *,"xleft ",xstenMAC(-1,dirnormal+1)
        stop
       endif

       if ((enable_spectral.eq.1).or. &
           (enable_spectral.eq.2)) then
        if (bfact.ge.2) then
         hx=dxmin
        else if (bfact.eq.1) then
         ! do nothing
        else
         print *,"bfact invalid63"
         stop
        endif
       else if (enable_spectral.eq.0) then
        ! do nothing
       else if (enable_spectral.eq.3) then
        ! do nothing
       else
        print *,"enable_spectral invalid"
        stop
       endif

       uu=zero
       uu_estdt=zero

        ! first check that characteristics do not collide or that
        ! a cell does not become a vacuum.
       ifaceR=i+ii
       jfaceR=j+jj
       kfaceR=k+kk
       if ((ifaceR.le.growhi(1)).and. &
           (jfaceR.le.growhi(2)).and. &
           (kfaceR.le.growhi(3))) then
        uleft=velmac(D_DECL(i,j,k))
        uright=velmac(D_DECL(ifaceR,jfaceR,kfaceR))
        if (uleft*uright.le.zero) then
         uu=max(uu,abs(uleft-uright))
         uu_estdt=max(uu_estdt,abs(uleft-uright))
        endif
       endif
      
       uleftcell=velcell(D_DECL(i-ii,j-jj,k-kk),dirnormal+1)
       urightcell=velcell(D_DECL(i,j,k),dirnormal+1)
       uu=max(uu,abs(uleftcell))
       uu=max(uu,abs(urightcell))
       uu_estdt=max(uu_estdt,abs(uleftcell))
       uu_estdt=max(uu_estdt,abs(urightcell))

       if ((enable_spectral.eq.1).or. &
           (enable_spectral.eq.2)) then
        if (bfact.ge.2) then
         velsum=zero
         do dir2=1,SDIM
          uleftcell=velcell(D_DECL(i-ii,j-jj,k-kk),dir2)
          urightcell=velcell(D_DECL(i,j,k),dir2)
          udiffcell=abs(uleftcell-urightcell)
          umaxcell=abs(uleftcell)
          if (abs(urightcell).gt.umaxcell) then
           umaxcell=abs(urightcell)
          endif
          if (abs(udiffcell).gt.umaxcell) then
           umaxcell=abs(udiffcell)
          endif
          velsum=velsum+umaxcell
         enddo ! dir2
         uu_estdt=max(uu_estdt,velsum)
        else if (bfact.eq.1) then
         ! do nothing
        else
         print *,"bfact invalid64"
         stop
        endif
       else if ((enable_spectral.eq.0).or. &
                (enable_spectral.eq.3)) then
        ! do nothing
       else
        print *,"enable_spectral invalid"
        stop
       endif

       do partid=0,nparts-1
        velcomp=partid*SDIM+dirnormal+1
        uu=max(uu,abs(solidfab(D_DECL(i,j,k),velcomp)))
        uu_estdt=max(uu_estdt,abs(solidfab(D_DECL(i,j,k),velcomp)))
       enddo ! partid=0..nparts-1

       uulocal=abs(velmac(D_DECL(i,j,k)))
       if (levelrz.eq.0) then
        ! do nothing
       else if ((levelrz.eq.1).or. &
                (levelrz.eq.3)) then
        if ((dirnormal.eq.0).and. &
            (xstenMAC(0,1).gt.VOFTOL*dx(1))) then
         if (xstenMAC(0,1).ge.hxmac) then
          uulocal=uulocal/( one-three*hxmac/(four*xstenMAC(0,1)) )  
         else if ((xstenMAC(0,1).gt.zero).and. &
                  (xstenMAC(-1,1).gt.zero)) then
          uulocal=four*uulocal
         else
          print *,"xstenMAC invalid estdt 2"
          print *,"i,j,k,dirnormal ",i,j,k,dirnormal
          print *,"hx=",hx
          print *,"hxmac=",hxmac
          print *,"xstenMAC(0,1)= ",xstenMAC(0,1)
          stop
         endif
        endif
       else
        print *,"levelrz invalid estdt"
        stop
       endif
       uu=max(uu,uulocal)
       uu_estdt=max(uu_estdt,uulocal)

       cleft=zero
       cright=zero
       cc=zero

       cleft_diag=zero
       cright_diag=zero
       cc_diag=zero

       do im=1,nmat
        if (fort_denconst(im).gt.zero) then
         if (fort_viscosity_state_model(im).ge.0) then
          if (visc_coef.ge.zero) then
           if (elastic_time(im).ge.zero) then
            if (elastic_time(im).eq.zero) then
             local_elastic_time=one
            else if (elastic_time(im).ge.one) then
             local_elastic_time=one
            else
             local_elastic_time=elastic_time(im)
            endif

             ! rho u_t = div beta (grad X + grad X^T)/2
             ! kg/m^3  m/s^2  = (1/m^2) beta m
             ! kg/(m^2 s^2) = (1/m) beta
             ! beta = kg/(m s^2)
             ! beta/rho = kg/(m s^2)   / (kg/m^3) = m^2/s^2
            elastic_wave_speed=visc_coef*fort_elastic_viscosity(im)/ &
                (local_elastic_time*fort_denconst(im))
            if (elastic_wave_speed.gt.zero) then
             elastic_wave_speed=sqrt(elastic_wave_speed)
             dthold=hx/elastic_wave_speed
             dt_min=min(dt_min,dthold)
            else if (elastic_wave_speed.eq.zero) then
             ! do nothing
            else
             print *,"elastic_wave_speed invalid"
             stop
            endif
           else
            print *,"elastic_time invalid"
            print *,"im= ",im
            print *,"elastic_time(im)=",elastic_time(im)
            stop
           endif
          else
           print *,"visc_coef invalid"
           stop
          endif
         else
          print *,"fort_viscosity_state_model invalid"
          stop
         endif
        else
         print *,"fort_denconst(im) invalid"
         stop
        endif
       enddo ! im=1..nmat

       do im=1,nmat
        LSleft(im)=dist(D_DECL(i-ii,j-jj,k-kk),im)
        LSright(im)=dist(D_DECL(i,j,k),im)
       enddo
       call get_primary_material(LSleft,nmat,im_primaryL)
       call get_primary_material(LSright,nmat,im_primaryR)

       USTEFAN=zero

        ! fluid region
       if ((is_rigid(nmat,im_primaryL).eq.0).and. &
           (is_rigid(nmat,im_primaryR).eq.0)) then 

        do ireverse=0,1
        do im=1,nmat-1
        do im_opp=im+1,nmat

         if ((im.gt.nmat).or.(im_opp.gt.nmat)) then
          print *,"im or im_opp bust 7"
          stop
         endif
         call get_iten(im,im_opp,iten,nmat)

         LL=latent_heat(iten+ireverse*nten)
         K_f=reaction_rate(iten+ireverse*nten)
         local_freezing_model=freezing_model(iten+ireverse*nten)
         local_Tanasawa_or_Schrage_or_Kassemi= &
           Tanasawa_or_Schrage_or_Kassemi(iten+ireverse*nten)
         distribute_from_targ=distribute_from_target(iten+ireverse*nten)
         TSAT=saturation_temp(iten+ireverse*nten)

         if ((local_freezing_model.eq.2).and.(num_species_var.ne.1)) then
          print *,"must define species var if hydrate model"
          stop
         endif

         if (ireverse.eq.0) then
          im_source=im
          im_dest=im_opp
         else if (ireverse.eq.1) then
          im_source=im_opp
          im_dest=im
         else
          print *,"ireverse invalid"
          stop
         endif
         dcompsrc=(im_source-1)*num_state_material+1
         tcompsrc=(im_source-1)*num_state_material+2
         vofcompsrc=(im_source-1)*ngeom_raw+1 
         dcompdst=(im_dest-1)*num_state_material+1
         tcompdst=(im_dest-1)*num_state_material+2
         vofcompdst=(im_dest-1)*ngeom_raw+1 

         ispec=mass_fraction_id(iten+ireverse*nten)
         if ((ispec.ge.1).and.(ispec.le.num_species_var)) then
          ! do nothing
         else if (ispec.eq.0) then
          if ((local_freezing_model.eq.4).or. & ! Tanasawa or Schrage
              (local_freezing_model.eq.5).or. & ! Stefan evap model
              (local_freezing_model.eq.2)) then ! hydrate
           print *,"ispec invalid"
           stop
          endif
         else
          print *,"ispec invalid"
          stop
         endif

         if ((is_rigid(nmat,im).eq.1).or. &
             (is_rigid(nmat,im_opp).eq.1)) then
          ! do nothing
         else if (LL.ne.zero) then

          do side=1,2
           if (side.eq.1) then
            icell=i-ii
            jcell=j-jj
            kcell=k-kk
            ialt=i
            jalt=j
            kalt=k
           else if (side.eq.2) then
            ialt=i-ii
            jalt=j-jj
            kalt=k-kk
            icell=i
            jcell=j
            kcell=k
           else
            print *,"side invalid"
            stop
           endif

           VOFsrc=vof(D_DECL(icell,jcell,kcell),vofcompsrc)
           VOFdst=vof(D_DECL(icell,jcell,kcell),vofcompdst)

           call gridsten_level(xsten,icell,jcell,kcell,level,nhalf)

            ! adjust LS1 if R-theta.
           call derive_dist(xsten,nhalf, &
            dist,DIMS(dist),icell,jcell,kcell,im_source,LS1)
            ! adjust LS2 if R-theta.
           call derive_dist(xsten,nhalf, &
            dist,DIMS(dist),icell,jcell,kcell,im_dest,LS2)
           
           Tsrc=den(D_DECL(icell,jcell,kcell),tcompsrc)
           Tdst=den(D_DECL(icell,jcell,kcell),tcompdst)
           Tsrcalt=den(D_DECL(ialt,jalt,kalt),tcompsrc)
           Tdstalt=den(D_DECL(ialt,jalt,kalt),tcompdst)
           Dsrc=den(D_DECL(icell,jcell,kcell),dcompsrc)
           Ddst=den(D_DECL(icell,jcell,kcell),dcompdst)

           if (LL.gt.zero) then ! evaporation
            vapor_den=Ddst
           else if (LL.lt.zero) then ! condensation
            vapor_den=Dsrc
           else
            print *,"LL invalid"
            stop
           endif

           if (local_freezing_model.eq.5) then ! stefan evap model
            if ((ispec.ge.1).and.(ispec.le.num_species_var)) then
             if (vapor_den.gt.zero) then
              ! do nothing
             else
              print *,"vapor_den invalid"
              stop
             endif  
            else
             print *,"ispec invalid"
             stop
            endif
           endif ! local_freezing_model==5

           Csrc=zero
           Cdst=zero

            ! hydrates
           if (local_freezing_model.eq.2) then
            if (distribute_from_targ.ne.0) then
             print *,"distribute_from_targ invalid"
             stop
            endif
            Csrc=den(D_DECL(icell,jcell,kcell),tcompsrc+1)
            Cdst=den(D_DECL(icell,jcell,kcell),tcompdst+1)
           endif

           if ((Dsrc.le.zero).or.(Ddst.le.zero)) then
            print *,"density must be positive estdt "
            print *,"im,im_opp ",im,im_opp
            print *,"im_source,im_dest ",im_source,im_dest
            print *,"Dsrc(source),Ddst(dest) ",Dsrc,Ddst
            stop
           endif

           if ((abs(LS1).le.two*dxmaxLS).and. &
               (abs(LS2).le.two*dxmaxLS)) then 

            delta=dxmin

            if ((LS1.ge.zero).or. &
                (LS2.ge.zero)) then

              ! coordinate of center of cell adjacent to face (i,j,k)
             do dir2=1,SDIM
              xI(dir2)=xsten(0,dir2)
             enddo
              ! either: 1-den_dst/den_src
              !     or: 1-den_src/den_dst
             if (fort_expansion_factor(iten+ireverse*nten).ge.one) then
              print *,"fort_expansion_factor(iten+ireverse*nten) bad"
              stop
             endif

             ksource=get_user_heatviscconst(im_source)
             kdest=get_user_heatviscconst(im_dest)
             alpha=fort_alpha(iten+ireverse*nten)
             beta=fort_beta(iten+ireverse*nten)
             if ((alpha.le.zero).or.(beta.le.zero)) then
              print *,"alpha or beta are invalid"
              stop
             endif
             dt_heat=one

             C_w0=fort_denconst(1)  ! density of water
             PHYDWATER=2.0D+19
             Cmethane_in_hydrate=Cdst ! hydrate is destination material.
       
             for_estdt=1

             source_perim_factor=one
             dest_perim_factor=one

             call get_vel_phasechange( &
              for_estdt, &
              xI, &
              ispec, &
              molar_mass, &
              species_molar_mass, &
              local_freezing_model, &
              local_Tanasawa_or_Schrage_or_Kassemi, &
              distribute_from_targ, &
              USTEFAN_hold, &
              Dsrc,Ddst, &
              Dsrc,Ddst, &
              ksource,kdest, &
              Tsrc,Tdst, &
              TSAT, &
              Tsrcalt,Tdstalt, &
              LL, &
              source_perim_factor, &
              dest_perim_factor, &
              microlayer_substrate(im_source), &
              microlayer_angle(im_source), &
              microlayer_size(im_source), &
              macrolayer_size(im_source), &
              microlayer_substrate(im_dest), &
              microlayer_angle(im_dest), &
              microlayer_size(im_dest), &
              macrolayer_size(im_dest), &
              delta, &  ! dxprobe_source
              delta, &  ! dxprobe_dest
              im_source,im_dest, &
              time,dt_heat, &
              alpha, &
              beta, &
              fort_expansion_factor(iten+ireverse*nten), &
              K_f, &
              Cmethane_in_hydrate, &
              C_w0, &
              PHYDWATER, &
              VOFsrc,VOFdst)

             if (USTEFAN_hold.lt.zero) then
              print *,"USTEFAN_hold.lt.zero"
              stop
             endif

             USTEFAN=USTEFAN+USTEFAN_hold

            else if ((LS1.lt.zero).and.(LS2.lt.zero)) then
             ! do nothing
            else
             print *,"LS1 or LS2 invalid"
             stop
            endif

           endif ! LS1,LS2 both close to 0
 
          enddo ! side 
         endif ! latent_heat<>0
        enddo ! im_opp
        enddo ! im
        enddo ! ireverse=0,1

       else if ((is_rigid(nmat,im_primaryL).eq.1).or. &
                (is_rigid(nmat,im_primaryR).eq.1)) then
        ! do nothing
       else
        print *,"im_primaryL, or im_primaryR invalid"
        stop 
       endif

       if (time.eq.zero) then
        ! do nothing
       else if (time.gt.zero) then
        USTEFAN=zero
        do dir2=1,SDIM
         if (USTEFAN.lt.abs(AMR_min_phase_change_rate(dir2))) then
          USTEFAN=abs(AMR_min_phase_change_rate(dir2))
         endif
         if (USTEFAN.lt.abs(AMR_max_phase_change_rate(dir2))) then
          USTEFAN=abs(AMR_max_phase_change_rate(dir2))
         endif
        enddo ! dir2=1..sdim
       else
        print *,"time invalid in ESTDT"
        stop
       endif

        ! factor of 4 in order to guarantee that characteristics do not
        ! collide.
        ! also, the factor of 4 should guarantee that a swept cell is not
        ! full at the end of CONVERTMATERIAL.
       if (min_stefan_velocity_for_dt.ge.zero) then
        if (USTEFAN.lt.min_stefan_velocity_for_dt) then
         USTEFAN=min_stefan_velocity_for_dt
        else if (USTEFAN.ge.min_stefan_velocity_for_dt) then
         ! do nothing
        else
         print *,"USTEFAN bust"
         stop
        endif
       else
        print *,"min_stefan_velocity_for_dt invalid"
        stop
       endif

       uu=abs(uu)+two*USTEFAN
       uu_estdt=abs(uu_estdt)+two*USTEFAN

       if (is_rigid(nmat,im_primaryL).eq.0) then
        ibase=(im_primaryL-1)*num_state_material
        density_left= &
          den(D_DECL(i-ii,j-jj,k-kk),ibase+1)

        if (material_type(im_primaryL).gt.0) then

         call init_massfrac_parm(density_left,massfrac_parm_left,im_primaryL)
         do ispec=1,num_species_var
          massfrac_parm_left(ispec)= &
            den(D_DECL(i-ii,j-jj,k-kk),ibase+2+ispec)
         enddo

         temperature_left=den(D_DECL(i-ii,j-jj,k-kk),ibase+2)
         call INTERNAL_material(density_left,massfrac_parm_left, &
          temperature_left, &
          internal_energy_left, &
          material_type(im_primaryL),im_primaryL)
         call SOUNDSQR_material(density_left,massfrac_parm_left, &
          internal_energy_left, &
          cleft_diag, &
          material_type(im_primaryL),im_primaryL)

         if ((shock_timestep(im_primaryL).eq.1).or. &
             ((shock_timestep(im_primaryL).eq.0).and.(time.eq.zero))) then
          cleft=cleft_diag
         endif
        else if (material_type(im_primaryL).eq.0) then
         ! do nothing
        else
         print *,"material_type(im_primaryL) invalid"
         stop
        endif 
       else if (is_rigid(nmat,im_primaryL).eq.1) then
        ! do nothing
       else
        print *,"is_rigid invalid"
        stop
       endif  

       if (is_rigid(nmat,im_primaryR).eq.0) then
        ibase=(im_primaryR-1)*num_state_material
        density_right=den(D_DECL(i,j,k),ibase+1)

        if (material_type(im_primaryR).gt.0) then

         call init_massfrac_parm(density_right,massfrac_parm_right,im_primaryR)
         do ispec=1,num_species_var
          massfrac_parm_right(ispec)= &
            den(D_DECL(i,j,k),ibase+2+ispec)
         enddo

         temperature_right=den(D_DECL(i,j,k),ibase+2)
         call INTERNAL_material(density_right,massfrac_parm_right, &
          temperature_right, &
          internal_energy_right, &
          material_type(im_primaryR),im_primaryR)
         call SOUNDSQR_material(density_right,massfrac_parm_right, &
          internal_energy_right, &
          cright_diag, &
          material_type(im_primaryR),im_primaryR)

         if ((shock_timestep(im_primaryR).eq.1).or. &
             ((shock_timestep(im_primaryR).eq.0).and.(time.eq.zero))) then
          cright=cright_diag
         endif

         if (im_primaryR.eq.im_primaryL) then
          if (density_left.gt.denmax) then
           denmax=density_left
          endif
          if (density_right.gt.denmax) then
           denmax=density_right
          endif
          denjump_temp=abs(density_left-density_right)
          if (denjump_temp.gt.denjump) then
           denjump=denjump_temp
          endif
         endif
        else if (material_type(im_primaryR).eq.0) then
         ! do nothing
        else
         print *,"material_type(im_primaryR) invalid"
         stop
        endif 
       else if (is_rigid(nmat,im_primaryR).eq.1) then
        ! do nothing
       else
        print *,"is_rigid invalid"
        stop
       endif

       cc_diag=max(cleft_diag,cright_diag)  ! c^2
       cc=max(cleft,cright)  ! c^2

       call check_user_defined_velbc(time,dirnormal,uu_estdt,dx)
       call check_user_defined_velbc(time,dirnormal,uu,dx)

       u_core = max(u_core,abs(uu))
       u_core_estdt = max(u_core_estdt,abs(uu_estdt))
       c_core = max(c_core,abs(cc_diag))  ! c^2

       if (uu_estdt.lt.uu) then
        print *,"uu_estdt invalid"
        stop
       endif
       if (u_core_estdt.lt.u_core) then
        print *,"u_core_estdt invalid"
        stop
       endif
       effective_velocity=abs(uu_estdt/weymouth_factor)+sqrt(cc)
       if (effective_velocity.gt.zero) then
        dt_min=min(dt_min,hx/effective_velocity)
       else if (effective_velocity.eq.zero) then
        ! do nothing
       else
        print *,"effective_velocity invalid"
        stop
       endif

        ! fluid region
       if ((is_rigid(nmat,im_primaryL).eq.0).and. &
           (is_rigid(nmat,im_primaryR).eq.0)) then 
        call fluid_interface(LSleft,LSright,gradh,im_opp,im,nmat)
       else if ((is_rigid(nmat,im_primaryL).eq.1).or. &
                (is_rigid(nmat,im_primaryR).eq.1)) then
        gradh=zero
       else
        print *,"im_primaryL, or im_primaryR invalid"
        stop 
       endif

       if (gradh.ne.zero) then

        if ((im.gt.nmat).or.(im_opp.gt.nmat)) then
         print *,"im or im_opp bust 8"
         stop
        endif
        call get_iten(im,im_opp,iten,nmat)

        if (level_cap_wave_speed(iten).lt.zero) then
         print *,"level_cap wave speed not initialized"
         stop
        else if (level_cap_wave_speed(iten).eq.zero) then
         ! do nothing
        else
         dthold=hx/level_cap_wave_speed(iten)
         dt_min=min(dt_min,dthold)
        endif

       else if (gradh.eq.zero) then
        ! do nothing
       else
        print *,"gradh invalid"
        stop
       endif 

      enddo
      enddo
      enddo  ! i,j,k

      if (abs(ns_gravity).gt.zero) then

       if (denmax.gt.zero) then

        if (terminal_velocity_dt.eq.0) then

         if (denjump.gt.zero) then
          if (denmax.gt.zero) then
           denjump=denjump/denmax
          else
           print *,"denmax became corrupt"
           stop
          endif
!   "delta rho" accounts for buoancy effects.
!   u dt < dx  u=(ubase + g (delta rho) dt)  
!   (ubase*dt+g (delta rho) dt^2) = dx   
!   g (delta rho) dt^2 + ubase dt - dx =0 
!   dt=(-ubase + sqrt(ubase^2 + 4 g (delta rho) dx))/(2 g drho) =
!      4 g (delta rho) dx/( 2 g (delta rho) 
!             (ubase+sqrt(ubase^2 + 4 g (delta rho) dx)))=
!      2 dx / (ubase+sqrt(ubase^2 + 4 g (delta rho) dx))=
!      sqrt(dx/(g (delta rho))) if ubase=0
!      dx/ubase          if g=0
          ugrav=half*(uu_estdt+sqrt(uu_estdt**2+  &
                  four*abs(denjump*ns_gravity)*dxmin))
          if (ugrav.gt.zero) then
           dthold=dxmin/ugrav 
           dt_min=min(dt_min,dthold)
          else
           print *,"ugrav invalid 1"
           print *,"uu_estdt ",uu_estdt
           print *,"denjump ",denjump
           print *,"ns_gravity ",ns_gravity
           print *,"dxmin ",dxmin
           print *,"denmax ",denmax
           stop
          endif
         else if (denjump.eq.zero) then
          ! do nothing
         else
          print *,"denjump cannot be negative, denjump = ",denjump
          stop
         endif

         if (denjump_gravity.gt.zero) then
          if (denmax_gravity.gt.zero) then
           denjump_gravity=denjump_gravity/denmax_gravity
           ugrav=half*(uu_estdt+sqrt(uu_estdt**2+  &
                  four*abs(denjump_gravity*ns_gravity)*dxmin))
           if (ugrav.gt.zero) then
            dthold=dxmin/ugrav 
            dt_min=min(dt_min,dthold)
           else
            print *,"ugrav invalid 2"
            print *,"uu_estdt ",uu_estdt
            print *,"denjump_gravity ",denjump_gravity
            print *,"ns_gravity ",ns_gravity
            print *,"dxmin ",dxmin
            print *,"denmax ",denmax
            stop
           endif
          else if (denmax_gravity.eq.zero) then
           ! do nothing
          else
           print *,"denmax_gravity invalid"
           stop
          endif
         else if (denjump_gravity.eq.zero) then
          ! do nothing
         else
          print *,"denjump_gravity cant be neg, denjump_gravity = ", &
            denjump_gravity
          stop
         endif

        else if (terminal_velocity_dt.eq.1) then

         denjump_terminal=fort_denconst(2)-fort_denconst(1)
         if (denjump_terminal.eq.zero) then
          denjump_terminal=denconst_gravity(2)-denconst_gravity(1)
          if (denjump_terminal.eq.zero) then
           print *,"denjump_terminal invalid"
           stop
          else if (denjump_terminal.gt.zero) then
           ! do nothing
          else
           print *,"denjump_terminal invalid"
           stop
          endif
         else if (denjump_terminal.gt.zero) then
          ! do nothing
         else
          print *,"denjump_terminal invalid"
          stop
         endif
          
         if (denjump_terminal.gt.zero) then
          if (fort_viscconst(1).gt.zero) then
           if (abs(ns_gravity).gt.zero) then
            if (radblob.gt.zero) then
             ! units: kg/m^3  m s/kg   m/s^2  m^2=(s/m^2) m^3/s^2=m/s
             v_terminal=(two/nine)*(denjump_terminal)* &
              (one/fort_viscconst(1))*abs(ns_gravity)*(radblob**2) 
             v_terminal=two*v_terminal
             if (v_terminal.gt.zero) then
              dthold=dxmin/v_terminal
              dt_min=min(dt_min,dthold)

              if (1.eq.0) then
               print *,"v_terminal= ",v_terminal
               print *,"dthold= ",dthold
               print *,"dt_min= ",dt_min
               print *,"dirnormal= ",dirnormal
               print *,"u_core= ",u_core
               print *,"u_core_estdt= ",u_core_estdt
               print *,"uu ",uu
               print *,"uu_estdt ",uu_estdt
               print *,"u_max(dirnormal+1)= ",u_max(dirnormal+1)
               print *,"u_max_estdt(dirnormal+1)= ",u_max_estdt(dirnormal+1)
              endif
             else
              print *,"v_terminal invalid"
              stop
             endif
            else
             print *,"radblob invalid"
             stop
            endif
           else
            print *,"ns_gravity invalid"
            stop
           endif
          else
           print *,"fort_viscconst(1) invalid"
           stop
          endif
         else
          print *,"denjump_terminal invalid"
          stop
         endif

        else
         print *,"terminal_velocity_dt invalid"
         stop
        endif

       else
        print *,"denmax invalid"
        stop
       endif

      else if (abs(ns_gravity).eq.zero) then
       ! do nothing
      else
       print *,"ns_gravity bust"
       stop
      endif 

      if (u_max(dirnormal+1).lt.u_core) then
       u_max(dirnormal+1)=u_core
      endif
      if (u_max_estdt(dirnormal+1).lt.u_core_estdt) then
       u_max_estdt(dirnormal+1)=u_core_estdt
      endif
      if (u_max(SDIM+1).lt.c_core) then
       u_max(SDIM+1)=c_core  ! c^2
      endif
      if (u_max_estdt(SDIM+1).lt.c_core) then
       u_max_estdt(SDIM+1)=c_core  ! c^2
      endif

      return
      end subroutine FORT_ESTDT


       ! (vel/RR) if R-THETA
       ! vel=vel*dt , call adjust_du if RZ, override if passive advect.
      subroutine FORT_VELMAC_OVERRIDE( &
       nmat, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       velbc, &
       dt,time, &
       passive_veltime, &
       vel_time, &
       dir_absolute_direct_split, &
       normdir, &
       utemp,DIMS(utemp), &
       unode,DIMS(unode), &
       ucell,DIMS(ucell), &
       xlo,dx, &
       mac_grow, &
       map_forward, &
       level, &
       finest_level, &
       SDC_outer_sweeps, &
       ns_time_order, &
       divu_outer_sweeps, &
       num_divu_outer_sweeps)
      use godunov_module
      use probf90_module
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: SDC_outer_sweeps
      INTEGER_T, intent(in) :: ns_time_order
      INTEGER_T, intent(in) :: divu_outer_sweeps
      INTEGER_T, intent(in) :: num_divu_outer_sweeps
      INTEGER_T, intent(in) :: level,finest_level
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: dir_absolute_direct_split
      INTEGER_T, intent(in) :: normdir
      INTEGER_T, intent(in) :: mac_grow,map_forward
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: dt,time,vel_time,passive_veltime
      INTEGER_T, intent(in) :: DIMDEC(utemp)
      INTEGER_T, intent(in) :: DIMDEC(unode)
      INTEGER_T, intent(in) :: DIMDEC(ucell)
     
      REAL_T, intent(in) :: utemp(DIMV(utemp)) 
      REAL_T, intent(inout) :: unode(DIMV(unode)) 
      REAL_T, intent(inout) :: ucell(DIMV(ucell),SDIM) 
      INTEGER_T, intent(in) :: velbc(SDIM,2)

      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
     
      INTEGER_T i,j,k
      INTEGER_T ii,jj,kk
      INTEGER_T idx,side
      REAL_T delta
      REAL_T hx
      REAL_T RR
      REAL_T xstenMAC(-3:3,SDIM)
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf
      INTEGER_T localbc
      INTEGER_T local_mac_grow

      nhalf=3

      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact invalid65"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
 
      if ((mac_grow.ne.1).and.(mac_grow.ne.2)) then
       print *,"mac_grow invalid mac_grow=",mac_grow
       stop
      endif
      if ((map_forward.ne.0).and.(map_forward.ne.1)) then
       print *,"map_forward invalid"
       stop
      endif
      if ((normdir.ge.0).and.(normdir.lt.SDIM)) then
       ! do nothing
      else
       print *,"normdir invalid"
       stop
      endif
      if ((dir_absolute_direct_split.ge.0).and. &
          (dir_absolute_direct_split.lt.SDIM)) then
       ! do nothing
      else
       print *,"dir_absolute_direct_split invalid"
       stop
      endif
      if (fabhi(normdir+1)-fablo(normdir+1)+1.lt.4) then
       print *,"blocking factor should be at least 4"
       stop
      endif
      if (level.gt.finest_level) then
       print *,"finest_level invalid velmac override"
       stop
      else if (level.lt.0) then
       print *,"level invalid velmac override"
       stop
      endif
      if ((SDC_outer_sweeps.ge.0).and. &
          (SDC_outer_sweeps.lt.ns_time_order)) then
       ! do nothing
      else
       print *,"SDC_outer_sweeps invalid"
       stop
      endif
      if ((divu_outer_sweeps.lt.0).or. &
          (divu_outer_sweeps.ge.num_divu_outer_sweeps)) then
       print *,"divu_outer_sweeps invalid FORT_VELMAC_OVERRIDE"
       stop
      endif
  
      if (levelrz.eq.0) then
       ! do nothing
      else if (levelrz.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (levelrz.eq.3) then
       ! do nothing
      else
       print *,"levelrz invalid velmac override"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(utemp),mac_grow+1,normdir,12)
      call checkbound(fablo,fabhi,DIMS(unode),mac_grow+1,normdir,12)
      call checkbound(fablo,fabhi,DIMS(ucell),mac_grow,-1,12)

      if (dt.le.zero) then
        print *,"dt invalid"
        stop
      endif

      ii=0
      jj=0
      kk=0
      if (normdir.eq.0) then
       ii=1
      else if (normdir.eq.1) then
       jj=1
      else if ((normdir.eq.2).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"normdir invalid"
       stop
      endif

        ! 1. multiply velocity by dt.
        ! 2. adjust velocity if RZ.
        ! 3. override velocity if it is a passive advection problem.
        ! 4. copy into mac_velocity
        ! 5. repeat for cell_velocity

      local_mac_grow=mac_grow+1

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,local_mac_grow,normdir,28)

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

        if (normdir.eq.0) then
         idx=i
        else if (normdir.eq.1) then
         idx=j
        else if ((normdir.eq.2).and.(SDIM.eq.3)) then
         idx=k
        else
         print *,"normdir invalid"
         stop
        endif

        call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,normdir,32)
        hx=xstenMAC(1,normdir+1)-xstenMAC(-1,normdir+1)
        if (hx.gt.zero) then
         ! do nothing
        else
         print *,"xstenMAC bust"
         stop
        endif

        RR=one
        if (levelrz.eq.0) then
         ! do nothing
        else if (levelrz.eq.1) then
         if (SDIM.ne.2) then
          print *,"dimension bust"
          stop
         endif
        else if (levelrz.eq.3) then
         if (normdir.eq.1) then
          RR=xstenMAC(0,1)
         endif
        else
         print *,"levelrz invalid velmac override"
         stop
        endif

        delta=utemp(D_DECL(i,j,k))

        side=0

        if (idx.le.fablo(normdir+1)) then
         side=1
        else if (idx.ge.fabhi(normdir+1)+1) then
         side=2
        else if ((idx.gt.fablo(normdir+1)).and. &
                 (idx.lt.fabhi(normdir+1)+1)) then
         ! do nothing
        else
         print *,"idx invalid"
         stop
        endif

        if ((side.eq.1).or.(side.eq.2)) then
         localbc=velbc(normdir+1,side)
         if (localbc.eq.REFLECT_ODD) then
          delta=zero
         else if (localbc.eq.EXT_DIR) then
          call velbc_override(vel_time,normdir+1,side,normdir+1, &
           delta, &
           xstenMAC,nhalf,dx,bfact)
         else if (localbc.eq.INT_DIR) then
          ! do nothing
         else if (localbc.eq.REFLECT_EVEN) then
          ! do nothing
         else if (localbc.eq.FOEXTRAP) then
          ! do nothing
         else
          print *,"localbc invalid"
          stop
         endif  ! cases for localbc 
        else if (side.eq.0) then
         ! do nothing
        else
         print *,"side invalid"
         stop
        endif  

        delta=dt*delta/RR

        ! modifies "delta" if "passive_advect_flag=1" or RZ. 
        call departure_node_split( &
          xstenMAC,nhalf,dx,bfact, &
          delta,passive_veltime, &
          normdir,dt,map_forward)

        if (abs(delta).ge.(one-0.001)*hx) then
         print *,"in: velmac_override"
         print *,"MAC: displacement exceeds grid cell"
         print *,"reduce cfl"
         print *,"utemp ",utemp(D_DECL(i,j,k))
         print *,"delta (u dt) = ",delta
         print *,"hx=    ",hx
         print *,"dt=    ",dt
         print *,"dir_absolute_direct_split= ",dir_absolute_direct_split
         print *,"normdir= ",normdir
         print *,"i,j,k ",i,j,k 
         print *,"level,finest_level ",level,finest_level
         print *,"SDC_outer_sweeps,ns_time_order ", &
            SDC_outer_sweeps,ns_time_order
         print *,"levelrz=",levelrz
         stop
        endif
      
         ! find displacements 
        unode(D_DECL(i,j,k))=delta

      enddo
      enddo
      enddo  ! i,j,k

      call growntilebox(tilelo,tilehi,fablo,fabhi, &
         growlo,growhi,mac_grow) 

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

        if (normdir.eq.0) then
         idx=i
        else if (normdir.eq.1) then
         idx=j
        else if ((normdir.eq.2).and.(SDIM.eq.3)) then
         idx=k
        else
         print *,"normdir invalid"
         stop
        endif

        call gridsten_level(xsten,i,j,k,level,nhalf)
        hx=xsten(1,normdir+1)-xsten(-1,normdir+1)
        if (hx.le.zero) then
         print *,"xsten bust"
         stop
        endif

        RR=one
        if (levelrz.eq.0) then
         ! do nothing
        else if (levelrz.eq.1) then
         if (SDIM.ne.2) then
          print *,"dimension bust"
          stop
         endif
        else if (levelrz.eq.3) then
         if (normdir.eq.1) then
          RR=xsten(0,1)
         endif
        else
         print *,"levelrz invalid velmac override"
         stop
        endif

        delta=ucell(D_DECL(i,j,k),normdir+1)

        side=0

        if (idx.lt.fablo(normdir+1)) then
         side=1
        else if (idx.gt.fabhi(normdir+1)) then
         side=2
        else if ((idx.ge.fablo(normdir+1)).and. &
                 (idx.le.fabhi(normdir+1))) then
         ! do nothing
        else
         print *,"idx invalid"
         stop
        endif

        if ((side.eq.1).or.(side.eq.2)) then
         localbc=velbc(normdir+1,side)

         if (localbc.eq.REFLECT_ODD) then
          delta=zero
         else if (localbc.eq.EXT_DIR) then
          call velbc_override(vel_time,normdir+1,side,normdir+1, &
           delta, &
           xsten,nhalf,dx,bfact)
         else if (localbc.eq.INT_DIR) then
          ! do nothing
         else if (localbc.eq.REFLECT_EVEN) then
          ! do nothing
         else if (localbc.eq.FOEXTRAP) then
          ! do nothing
         else
          print *,"localbc invalid"
          stop
         endif  ! cases for localbc 
        else if (side.eq.0) then
         ! do nothing
        else
         print *,"side invalid"
         stop
        endif  

        delta=delta*dt/RR

          ! modifies "delta" if "passive_advect_flag=1" or RZ. 
        call departure_node_split( &
          xsten,nhalf,dx,bfact, &
          delta,passive_veltime, &
          normdir,dt,map_forward)

        if (abs(delta).ge.(one-0.001)*hx) then
         print *,"in: velmac_override"
         print *,"CELL: displacement exceeds grid cell"
         print *,"reduce cfl"
         print *,"ucell ",ucell(D_DECL(i,j,k),normdir+1)
         print *,"delta= ",delta
         print *,"hx=    ",hx
         print *,"dt=    ",dt
         print *,"dir_absolute_direct_split= ",dir_absolute_direct_split
         print *,"normdir= ",normdir
         print *,"i,j,k ",i,j,k 
         print *,"level,finest_level ",level,finest_level
         print *,"SDC_outer_sweeps,ns_time_order ", &
            SDC_outer_sweeps,ns_time_order
         stop
        endif
        
        ucell(D_DECL(i,j,k),normdir+1)=delta
      enddo
      enddo
      enddo  ! i,j,k

      return
      end subroutine FORT_VELMAC_OVERRIDE


       ! mask=1 if cell not covered.
       ! masknbr=1 fine-fine border cells and interior cells.
       ! masknbr=0 coarse-fine cells and cells outside domain.
       ! called from getStateMOM_DEN
      subroutine FORT_DERIVE_MOM_DEN( &
       im_parm, &
       ngrow, &
       constant_density_all_time, & ! 1..nmat
       spec_material_id_AMBIENT, &  ! 1..num_species_var
       presbc_arr, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       dt, &
       mask,DIMS(mask), &
       masknbr,DIMS(masknbr), &
       vol,DIMS(vol), &
       eosdata,DIMS(eosdata), &
       momden,DIMS(momden), &
       recon,DIMS(recon), &
       xlo,dx, &
       gravity_normalized, &
       DrhoDT, &
       override_density, &
       nmat, &
       level,finest_level)
      use probf90_module
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: im_parm
      INTEGER_T, intent(in) :: ngrow
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: constant_density_all_time(nmat)
      INTEGER_T, intent(in) :: spec_material_id_AMBIENT(num_species_var+1)
      INTEGER_T, intent(in) :: level,finest_level
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact

      REAL_T, intent(in) :: dt
      INTEGER_T, intent(in) :: DIMDEC(vol)
      INTEGER_T, intent(in) :: DIMDEC(eosdata)
      INTEGER_T, intent(in) :: DIMDEC(momden)
      INTEGER_T, intent(in) :: DIMDEC(recon)
      INTEGER_T, intent(in) :: DIMDEC(mask)
      INTEGER_T, intent(in) :: DIMDEC(masknbr)
     
      REAL_T, intent(in) ::  mask(DIMV(mask)) 
      REAL_T, intent(in) ::  masknbr(DIMV(masknbr)) 
      REAL_T, intent(in) ::  vol(DIMV(vol)) 
      REAL_T, intent(in) :: eosdata(DIMV(eosdata),num_state_material*nmat)
      REAL_T, intent(out) :: momden(DIMV(momden),nmat)
      REAL_T, intent(in) :: recon(DIMV(recon),nmat*ngeom_recon)

      INTEGER_T, intent(in) :: presbc_arr(SDIM,2)

      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T, intent(in) :: override_density(nmat)
      REAL_T, intent(in) :: DrhoDT(nmat)
      REAL_T, intent(in) :: gravity_normalized
     
      INTEGER_T i,j,k
      INTEGER_T dir

      INTEGER_T dencomp
      REAL_T xpos(SDIM)
      REAL_T xsten(-3:3,SDIM)
      REAL_T rhohydro,preshydro,temperature
      INTEGER_T nhalf
      REAL_T density_of_TZ
      INTEGER_T caller_id
      REAL_T rho_base
      INTEGER_T vofcomp
      REAL_T local_vfrac

      nhalf=3

      if (bfact.lt.1) then
       print *,"bfact invalid66"
       stop
      endif
      if (ngrow.ge.1) then
       ! do nothing
      else
       print *,"ngrow>=1 required"
       stop
      endif

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid dencor"
       stop
      endif
      if ((im_parm.ge.1).and.(im_parm.le.nmat)) then
       ! do nothing
      else
       print *,"FORT_DERIVE_MOM_DEN: im_parm invalid, im_parm=",im_parm
       stop
      endif
      vofcomp=(im_parm-1)*ngeom_recon+1

      if (levelrz.eq.0) then
       ! do nothing
      else if (levelrz.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (levelrz.eq.3) then
       ! do nothing
      else
       print *,"levelrz invalid dencor"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(vol),ngrow,-1,12)
      call checkbound(fablo,fabhi,DIMS(eosdata),ngrow,-1,12)
      call checkbound(fablo,fabhi,DIMS(momden),ngrow,-1,12)
      call checkbound(fablo,fabhi,DIMS(recon),ngrow,-1,12)
      call checkbound(fablo,fabhi,DIMS(mask),ngrow,-1,12)
      call checkbound(fablo,fabhi,DIMS(masknbr),ngrow,-1,12)

      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,ngrow) 

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridsten_level(xsten,i,j,k,level,nhalf)
       do dir=1,SDIM
        xpos(dir)=xsten(0,dir)
       enddo

       if (mask(D_DECL(i,j,k)).eq.zero) then
          ! a default for coarse grid cells covered by finer levels.
        momden(D_DECL(i,j,k),im_parm)=fort_denconst(im_parm)
       else if (mask(D_DECL(i,j,k)).eq.one) then
     
        dencomp=(im_parm-1)*num_state_material+1

        if (constant_density_all_time(im_parm).eq.1) then
         rho_base=fort_denconst(im_parm)
        else if (constant_density_all_time(im_parm).eq.0) then
         rho_base=eosdata(D_DECL(i,j,k),dencomp)
        else
         print *,"constant_density_all_time(im_parm) invalid"
         stop
        endif

        ! rho=rho(T,z)
        if (override_density(im_parm).eq.1) then

         if (fort_material_type(im_parm).eq.0) then
          ! do nothing
         else
          print *,"override_density==1 for incomp material only"
          stop
         endif

         local_vfrac=recon(D_DECL(i,j,k),vofcomp)

         if ((local_vfrac.ge.VOFTOL).and.(local_vfrac.le.one+VOFTOL)) then

           ! den,T
          temperature=eosdata(D_DECL(i,j,k),dencomp+1)

           ! defined in: GLOBALUTIL.F90
           ! only takes into account fort_drhodz
          caller_id=0
          call default_hydrostatic_pressure_density( &
            xpos, &
            rho_base, &
            rhohydro, &
            preshydro, &
            temperature, &
            gravity_normalized, &
            im_parm, &
            override_density(im_parm), &
            caller_id)

          if (DrhoDT(im_parm).le.zero) then
           ! do nothing
          else
           print *,"DrhoDT invalid"
           stop
          endif

          density_of_TZ=rhohydro+ &
            rho_base*DrhoDT(im_parm)* &
            (temperature-fort_tempconst(im_parm))

          if ((temperature.ge.zero).and. &
              (rhohydro.gt.zero).and. &
              (fort_tempconst(im_parm).ge.zero).and. &
              (fort_denconst(im_parm).gt.zero).and. &
              (rho_base.gt.zero)) then 
           ! do nothing
          else
           print *,"invalid parameters to get the density"
           print *,"im_parm=",im_parm
           print *,"temperature=",temperature
           print *,"density_of_TZ=",density_of_TZ
           print *,"rho_base=",rho_base
           print *,"rhohydro=",rhohydro
           print *,"fort_tempconst(im_parm)=",fort_tempconst(im_parm)
           stop
          endif

          if (density_of_TZ.gt.zero) then
           ! do nothing
          else if (density_of_TZ.le.zero) then
           print *,"WARNING density_of_TZ.le.zero"
           print *,"im_parm=",im_parm
           print *,"temperature=",temperature
           print *,"density_of_TZ=",density_of_TZ
           print *,"rho_base=",rho_base
           print *,"rhohydro=",rhohydro
           print *,"fort_tempconst(im_parm)=",fort_tempconst(im_parm)
           print *,"fort_tempcutoffmax(im_parm)=",fort_tempcutoffmax(im_parm)
          
           temperature=fort_tempcutoffmax(im_parm)
  
           density_of_TZ=rhohydro+ &
            rho_base*DrhoDT(im_parm)* &
            (temperature-fort_tempconst(im_parm))

           if (density_of_TZ.gt.zero) then
            ! do nothing
           else
            print *,"density_of_TZ.le.zero (STILL)"
            stop
           endif

          else
           print *,"density_of_TZ bust"
           stop
          endif

          momden(D_DECL(i,j,k),im_parm)=density_of_TZ

         else if (abs(local_vfrac).le.VOFTOL) then

          momden(D_DECL(i,j,k),im_parm)=rho_base

         else
          print *,"local_vfrac invalid"
          stop
         endif

        else if ((override_density(im_parm).eq.0).or. &
                 (override_density(im_parm).eq.2)) then
         momden(D_DECL(i,j,k),im_parm)=rho_base
        else
         print *,"override_density invalid"
         stop
        endif

       else
        print *,"mask invalid"
        stop
       endif

      enddo
      enddo
      enddo ! i,j,k

      return
      end subroutine FORT_DERIVE_MOM_DEN


       ! adjust_temperature==1  modify temperature (Snew and coeff)
       ! adjust_temperature==0  modify coefficient (coeff)
       ! adjust_temperature==-1 modify heatx,heaty,heatz
       ! if project_option==2:
       !  heatxyz correspond to thermal diffusivity
       ! else if project_option>=100:
       !  heatxyz correspond to rho D
       !
      subroutine FORT_STEFANSOLVER( &
       project_option, & ! 2=thermal diffusion or 100...100+num_species_var-1
       solidheat_flag, & ! 0=diffuse in solid 1=dirichlet 2=Neumann
       microlayer_size, & ! 1..nmat
       microlayer_substrate, & ! 1..nmat
       microlayer_temperature_substrate, & ! 1..nmat
       adjust_temperature, &
       nmat, &
       nten, &
       nstate, &
       ntsat, &
       nden, &
       latent_heat, &
       freezing_model, &
       distribute_from_target, &
       saturation_temp, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       finest_level, &
       xlo,dx, &
       dt, &
       STATEFAB,DIMS(STATEFAB), &
       TgammaFAB,DIMS(TgammaFAB), &
       swept,DIMS(swept), &
       LS,DIMS(LS),  &
       T_fab,DIMS(T_fab),  &
       TorY_fab,DIMS(TorY_fab),  &
       Snew,DIMS(Snew), & 
       DeDT,DIMS(DeDT), &
       den,DIMS(den), &
       coeff,DIMS(coeff), &
       vol,DIMS(vol), &
       heatx,DIMS(heatx), &
       heaty,DIMS(heaty), &
       heatz,DIMS(heatz), &
       areax,DIMS(areax), &
       areay,DIMS(areay), &
       areaz,DIMS(areaz) )

      use probf90_module
      use global_utility_module
      use MOF_routines_module
      use mass_transfer_module
      use godunov_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: project_option
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: nstate
      INTEGER_T, intent(in) :: ntsat
      INTEGER_T, intent(in) :: nden

      INTEGER_T, intent(in) :: solidheat_flag
      REAL_T, intent(in) :: microlayer_size(nmat)
      INTEGER_T, intent(in) :: microlayer_substrate(nmat)
      REAL_T, intent(in) :: microlayer_temperature_substrate(nmat)
      
      INTEGER_T, intent(in) :: adjust_temperature
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      REAL_T, intent(in) :: latent_heat(2*nten)
      INTEGER_T, intent(in) :: freezing_model(2*nten)
      INTEGER_T, intent(in) :: distribute_from_target(2*nten)
      REAL_T, intent(in) :: saturation_temp(2*nten)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: dt
      INTEGER_T, intent(in) :: DIMDEC(STATEFAB)
      INTEGER_T, intent(in) :: DIMDEC(TgammaFAB)
      INTEGER_T, intent(in) :: DIMDEC(swept)
      INTEGER_T, intent(in) :: DIMDEC(LS)
      INTEGER_T, intent(in) :: DIMDEC(T_fab)
      INTEGER_T, intent(in) :: DIMDEC(TorY_fab)
      INTEGER_T, intent(in) :: DIMDEC(Snew)
      INTEGER_T, intent(in) :: DIMDEC(DeDT)
      INTEGER_T, intent(in) :: DIMDEC(den)
      INTEGER_T, intent(in) :: DIMDEC(coeff)
      INTEGER_T, intent(in) :: DIMDEC(vol)
      INTEGER_T, intent(in) :: DIMDEC(heatx)
      INTEGER_T, intent(in) :: DIMDEC(heaty)
      INTEGER_T, intent(in) :: DIMDEC(heatz)
      INTEGER_T, intent(in) :: DIMDEC(areax)
      INTEGER_T, intent(in) :: DIMDEC(areay)
      INTEGER_T, intent(in) :: DIMDEC(areaz)
      REAL_T, intent(in) :: STATEFAB(DIMV(STATEFAB),nden) 
      REAL_T, intent(in) :: TgammaFAB(DIMV(TgammaFAB),ntsat) 
      REAL_T, intent(in) :: swept(DIMV(swept),nmat)
      REAL_T, intent(in) :: LS(DIMV(LS),nmat*(SDIM+1))
      REAL_T, intent(in) :: T_fab(DIMV(T_fab),nmat)
      REAL_T, intent(in) :: TorY_fab(DIMV(TorY_fab),nmat)
      REAL_T, intent(out) :: Snew(DIMV(Snew),nstate)
      REAL_T, intent(in) :: DeDT(DIMV(DeDT),nmat+1)  ! 1/(rho cv) (cv=DeDT)
      ! 1/den (i.e. den actually stores 1/den)
      REAL_T, intent(in) :: den(DIMV(den),nmat+1)  
       ! alphanovolume or outer_iter_pressure
      REAL_T, intent(out) :: coeff(DIMV(coeff))  
      REAL_T, intent(in) :: vol(DIMV(vol))
       ! thermal conductivity
      REAL_T, intent(out) :: heatx(DIMV(heatx))
      REAL_T, intent(out) :: heaty(DIMV(heaty))
      REAL_T, intent(out) :: heatz(DIMV(heatz))
      REAL_T, intent(in) :: areax(DIMV(areax))
      REAL_T, intent(in) :: areay(DIMV(areay))
      REAL_T, intent(in) :: areaz(DIMV(areaz))

      INTEGER_T i,j,k
      INTEGER_T i1,j1,k1
      INTEGER_T k1lo,k1hi
      INTEGER_T dir,side
      INTEGER_T ii,jj,kk
      INTEGER_T ic,jc,kc
      INTEGER_T iface,jface,kface
      INTEGER_T at_interface
      INTEGER_T im_loop
      INTEGER_T im_primary  ! LS_center(im_primary)=max_im LS_center(im)
      INTEGER_T im_side_primary ! LS_side(im_side_primary)=max_im LS_side(im)
      INTEGER_T im_crit
      INTEGER_T im_crit_save
      INTEGER_T im_adjust
      INTEGER_T im
      INTEGER_T im_opp
      INTEGER_T ireverse
      INTEGER_T iten
      INTEGER_T im_source
      INTEGER_T im_dest
      INTEGER_T im_source_substrate
      INTEGER_T im_dest_substrate
      INTEGER_T ireverse_crit
      INTEGER_T iten_crit
      INTEGER_T im_source_crit
      INTEGER_T im_dest_crit
      INTEGER_T im_source_substrate_crit
      INTEGER_T im_dest_substrate_crit
      INTEGER_T local_freezing_model
      INTEGER_T distribute_from_targ
      REAL_T TGRAD_MAX,TGRAD_test,TorY_test
      REAL_T LL
      REAL_T Tgamma
      REAL_T TorYgamma_BC
      INTEGER_T Tgamma_STATUS
      INTEGER_T tsat_comp
      INTEGER_T ngrow_tsat
      REAL_T T_MIN(nmat)
      REAL_T T_MAX(nmat)
      REAL_T TorY_MIN(nmat)
      REAL_T TorY_MAX(nmat)
      INTEGER_T T_STATUS(nmat)
      INTEGER_T TorY_STATUS(nmat)

      INTEGER_T nten_test
      REAL_T over_den,over_cv
      REAL_T single_material_den
      REAL_T local_vol
      REAL_T original_coeff,delta_coeff,coeff_Tgamma
      REAL_T aface
      REAL_T LS_center(nmat)
      REAL_T LS_side(nmat)
      REAL_T LS_no_tess(nmat)
      REAL_T LS1,LS2
      REAL_T theta
      REAL_T theta_cutoff
      REAL_T heatcoeff
      REAL_T side_coeff,T_adjust
      INTEGER_T tcomp
      INTEGER_T start_freezing
      REAL_T SWEPTFACTOR,hx
      INTEGER_T ncomp_per_tsat
      REAL_T xsten(-3:3,SDIM)
      REAL_T xsten_side(-3:3,SDIM)
      REAL_T x_interface(SDIM)
      INTEGER_T nhalf
      INTEGER_T dir_inner
      REAL_T T_or_Y_min_sanity

      nhalf=3

      theta_cutoff=0.001

      ncomp_per_tsat=2
      if (ntsat.eq.nten*(ncomp_per_tsat+1)) then
       ! do nothing
      else
       print *,"nstat invalid"
       stop
      endif
      if (nden.eq.nmat*num_state_material) then
       ! do nothing
      else
       print *,"nden invalid"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid67"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base must be 2"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid ratemass nten, nten_test ",nten,nten_test
       stop
      endif
      if (nstate.ne. &
          (SDIM+1)+ &
          nmat*(num_state_material+ngeom_raw)+1) then
       print *,"nstate invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level or finest_level invalid stefan solver"
       stop
      endif
      if (dt.le.zero) then
       print *,"dt invalid"
       stop
      endif
       ! solidheat_flag==0 diffuse in solid
       ! solidheat_flag==1 dirichlet solid/fluid
       ! solidheat_flag==2 Neumann solid/fluid
      if ((solidheat_flag.lt.0).or. &
          (solidheat_flag.gt.2)) then
       print *,"solidheat_flag invalid"
       stop
      endif
      do im=1,nmat
       im_crit=microlayer_substrate(im)
       if (im_crit.eq.0) then
        ! do nothing
       else if ((im_crit.ge.1).and.(im_crit.le.nmat)) then
        if (is_rigid(nmat,im_crit).ne.1) then
         print *,"is_rigid(nmat,im_crit) invalid"
         stop
        endif
       else
        print *,"microlayer_substrate invalid"
        stop
       endif
      enddo ! im=1..nmat

      if (project_option.eq.2) then ! thermal diffusion
       T_or_Y_min_sanity=zero
      else if ((project_option.ge.100).and. & ! species diffusion
               (project_option.le.100+num_species_var-1)) then
       T_or_Y_min_sanity=zero
      else
       print *,"project_option invalid"
       stop
      endif

      if (SDIM.eq.2) then 
       k1lo=0
       k1hi=0
      else if (SDIM.eq.3) then
       k1lo=-1
       k1hi=1
      else
       print *,"dimension bust"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(STATEFAB),1,-1,234)
      call checkbound(fablo,fabhi,DIMS(TgammaFAB),1,-1,234)

      call checkbound(fablo,fabhi,DIMS(swept),0,-1,234)

      call checkbound(fablo,fabhi,DIMS(LS),1,-1,1226)
      call checkbound(fablo,fabhi,DIMS(T_fab),1,-1,1226)
      call checkbound(fablo,fabhi,DIMS(TorY_fab),1,-1,1226)
      call checkbound(fablo,fabhi,DIMS(Snew),1,-1,1227)
      call checkbound(fablo,fabhi,DIMS(DeDT),1,-1,1228) ! 1/(density * cv)
      call checkbound(fablo,fabhi,DIMS(den),1,-1,1229)  ! 1/(density)
      call checkbound(fablo,fabhi,DIMS(coeff),0,-1,1230)
      call checkbound(fablo,fabhi,DIMS(vol),0,-1,1231)
      call checkbound(fablo,fabhi,DIMS(heatx),0,0,1232)
      call checkbound(fablo,fabhi,DIMS(heaty),0,1,1233)
      call checkbound(fablo,fabhi,DIMS(heatz),0,SDIM-1,1234)
      call checkbound(fablo,fabhi,DIMS(areax),0,0,1235)
      call checkbound(fablo,fabhi,DIMS(areay),0,1,1236)
      call checkbound(fablo,fabhi,DIMS(areaz),0,SDIM-1,1237)
 
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridsten_level(xsten,i,j,k,level,nhalf)

       im_crit_save=-1 ! the dominant material in the center cell
       im_dest_crit=-1
       im_source_crit=-1
       im_dest_substrate_crit=-1
       im_source_substrate_crit=-1
       iten_crit=-1
       ireverse_crit=-1
       TGRAD_MAX=-1.0e+10

       do im=1,nmat
        LS_no_tess(im)=LS(D_DECL(i,j,k),im)
       enddo
       call LS_tessellate(LS_no_tess,LS_center,nmat)
       im_primary=1
       do im=2,nmat
        if (LS_center(im).gt.LS_center(im_primary)) then
         im_primary=im
        endif
       enddo

       if (LS_center(im_primary).ge.zero) then

        do im=1,nmat
         T_MIN(im)=zero
         T_MAX(im)=zero
         TorY_MIN(im)=zero
         TorY_MAX(im)=zero
         T_STATUS(im)=0
         TorY_STATUS(im)=0
        enddo ! im=1..nmat

        do i1=-1,1 
        do j1=-1,1 
        do k1=k1lo,k1hi
         do im=1,nmat
          LS_no_tess(im)=LS(D_DECL(i+i1,j+j1,k+k1),im)
         enddo
         call LS_tessellate(LS_no_tess,LS_side,nmat)
         im_side_primary=1
         do im=2,nmat
          if (LS_side(im).gt.LS_side(im_side_primary)) then
           im_side_primary=im
          endif
         enddo
         TorY_test=T_fab(D_DECL(i+i1,j+j1,k+k1),im_side_primary)
         if (TorY_test.ge.T_or_Y_min_sanity) then
          ! do nothing
         else
          print *,"TorY_test= ",TorY_test
          print *,"adjust_temperature=",adjust_temperature
          print *,"i,j,k ",i,j,k
          print *,"i1,j1,k1 ",i1,j1,k1
          print *,"im_side_primary ",im_side_primary
          print *,"im_primary ",im_primary
          do im=1,nmat
           print *,"im,T_fab ",im, &
            T_fab(D_DECL(i+i1,j+j1,k+k1),im) 
          enddo
          print *,"TorY_test.le.zero STEFANSOLVER"
          stop
         endif
         if (T_STATUS(im_side_primary).eq.0) then
          T_MIN(im_side_primary)=TorY_test
          T_MAX(im_side_primary)=TorY_test
         else if (T_STATUS(im_side_primary).eq.1) then
          if (TorY_test.gt.T_MAX(im_side_primary)) then
           T_MAX(im_side_primary)=TorY_test
          endif
          if (TorY_test.lt.T_MIN(im_side_primary)) then
           T_MIN(im_side_primary)=TorY_test
          endif
         else
          print *,"T_STATUS(im_side_primary) invalid"
          stop
         endif
         T_STATUS(im_side_primary)=1

         TorY_test=TorY_fab(D_DECL(i+i1,j+j1,k+k1),im_side_primary)
         if (TorY_test.ge.T_or_Y_min_sanity) then
          ! do nothing
         else
          print *,"TorY_test= ",TorY_test
          print *,"adjust_temperature=",adjust_temperature
          print *,"i,j,k ",i,j,k
          print *,"i1,j1,k1 ",i1,j1,k1
          print *,"im_side_primary ",im_side_primary
          print *,"im_primary ",im_primary
          do im=1,nmat
           print *,"im,TorY_fab ",im, &
            TorY_fab(D_DECL(i+i1,j+j1,k+k1),im) 
          enddo
          print *,"TorY_test.le.zero STEFANSOLVER"
          stop
         endif
         if (TorY_STATUS(im_side_primary).eq.0) then
          TorY_MIN(im_side_primary)=TorY_test
          TorY_MAX(im_side_primary)=TorY_test
         else if (TorY_STATUS(im_side_primary).eq.1) then
          if (TorY_test.gt.TorY_MAX(im_side_primary)) then
           TorY_MAX(im_side_primary)=TorY_test
          endif
          if (TorY_test.lt.TorY_MIN(im_side_primary)) then
           TorY_MIN(im_side_primary)=TorY_test
          endif
         else
          print *,"TorY_STATUS(im_side_primary) invalid"
          stop
         endif
         TorY_STATUS(im_side_primary)=1

        enddo
        enddo
        enddo ! i1,j1,k1

        do im=1,nmat-1
         do im_opp=im+1,nmat
          do ireverse=0,1
           if ((im.gt.nmat).or.(im_opp.gt.nmat)) then
            print *,"im or im_opp bust 8"
            stop
           endif
           call get_iten(im,im_opp,iten,nmat)
           LL=latent_heat(iten+ireverse*nten)
           Tgamma_STATUS=NINT(TgammaFAB(D_DECL(i,j,k),iten))
           if (ireverse.eq.0) then
            ! do nothing
           else if (ireverse.eq.1) then
            Tgamma_STATUS=-Tgamma_STATUS
           else
            print *,"ireverse invalid"
            stop
           endif

           if ((Tgamma_STATUS.eq.1).or.(Tgamma_STATUS.eq.2)) then

            local_freezing_model=freezing_model(iten+ireverse*nten)
            distribute_from_targ=distribute_from_target(iten+ireverse*nten)

            if ((distribute_from_targ.ne.0).and. &
                (distribute_from_targ.ne.1)) then
             print *,"distribute_from_targ invalid"
             stop
            endif

            if ((is_rigid(nmat,im).eq.0).and. &
                (is_rigid(nmat,im_opp).eq.0)) then

             if (LL.ne.zero) then

              if (ireverse.eq.0) then
               im_source=im
               im_dest=im_opp
              else if (ireverse.eq.1) then
               im_source=im_opp
               im_dest=im
              else
               print *,"ireverse invalid"
               stop
              endif
              call check_recalesce_status(im_source,start_freezing)

              if (start_freezing.eq.1) then

               ! local_freezing_model=0 (sharp interface stefan model)
               ! local_freezing_model=1 (source term model)
               ! local_freezing_model=2 (hydrate model)
               ! local_freezing_model=3 (wildfire)
               ! local_freezing_model=4 (Tanasawa/Schrage)
               ! local_freezing_model=5 (evaporation/condensation)
               ! local_freezing_model=6 (Palmore Desjardins)
               ! local_freezing_model=7 (Cavitation)
               if ((local_freezing_model.eq.0).or. &
                   (local_freezing_model.eq.5).or. &
                   (local_freezing_model.eq.6)) then

                if (project_option.eq.2) then
                   ! default Tgamma
                 Tgamma=saturation_temp(iten+ireverse*nten)
                 TorYgamma_BC=Tgamma
                 if (Tgamma.gt.zero) then
                  tsat_comp=nten+(iten-1)*ncomp_per_tsat+1
                  Tgamma=TgammaFAB(D_DECL(i,j,k),tsat_comp)
                  TorYgamma_BC=Tgamma
                  if (Tgamma.gt.zero) then
                   ! do nothing
                  else
                   print *,"Tgamma must be positive1"
                   stop
                  endif
                 else
                  print *,"saturation temperature must be positive2"
                  stop
                 endif
                else if ((project_option.ge.100).and. &
                         (project_option.lt.100+num_species_var)) then
                 Tgamma=saturation_temp(iten+ireverse*nten)
                 TorYgamma_BC=one
                 if (Tgamma.gt.zero) then
                  tsat_comp=nten+(iten-1)*ncomp_per_tsat+1
                  Tgamma=TgammaFAB(D_DECL(i,j,k),tsat_comp)
                  tsat_comp=nten+(iten-1)*ncomp_per_tsat+2
                  TorYgamma_BC=TgammaFAB(D_DECL(i,j,k),tsat_comp)
                  if (Tgamma.gt.zero) then
                   ! do nothing
                  else
                   print *,"Tgamma must be positive22"
                   stop
                  endif
                  if ((TorYgamma_BC.ge.zero).and.(TorYgamma_BC.le.one)) then
                   ! do nothing
                  else
                   print *,"TorYgamma_BC (aka Y) must be >= 0 and <=1"
                   stop
                  endif
                 else
                  print *,"saturation temperature must be positive33"
                  stop
                 endif
                else
                 print *,"project_option invalid"
                 stop
                endif

                im_source_substrate=im_source
                im_dest_substrate=im_dest

                if (local_freezing_model.eq.0) then ! stefan

                  ! Tgamma BC at thin filament interface.
                  if (solidheat_flag.eq.0) then ! diffuse in solid
                   if (microlayer_substrate(im_source).ne.0) then
                    im_source_substrate=microlayer_substrate(im_source)
                    if (is_rigid(nmat,im_source_substrate).ne.1) then
                     print *,"is_rigid(nmat,im_source_substrate) invalid"
                     stop
                    endif
                   endif
                   if (microlayer_substrate(im_dest).ne.0) then
                    im_dest_substrate=microlayer_substrate(im_dest)
                    if (is_rigid(nmat,im_dest_substrate).ne.1) then
                     print *,"is_rigid(nmat,im_dest_substrate) invalid"
                     stop
                    endif
                   endif
                  else if ((solidheat_flag.eq.1).or. & ! dirichlet 
                           (solidheat_flag.eq.2)) then ! neumann 
                   ! do nothing
                  else
                   print *,"solidheat_flag invalid"
                   stop
                  endif

                else if (local_freezing_model.eq.5) then ! stefan evap/diff
                  ! do nothing
                else if (local_freezing_model.eq.6) then !Palmore/Desjardins
                  ! do nothing
                else
                  print *,"local_freezing_model invalid 15"
                  stop
                endif

                if ((im_source.eq.im_primary).or. &
                    (im_source_substrate.eq.im_primary).or. &
                    (im_dest.eq.im_primary).or. &
                    (im_dest_substrate.eq.im_primary)) then

                  im_crit=im_primary

                   ! im_primary is found in the stencil
                  if (TorY_STATUS(im_crit).eq.1) then

                   if (LL.lt.zero) then ! freezing or condensation
                    TGRAD_test=Tgamma-T_MIN(im_crit)
                   else if (LL.gt.zero) then  ! melting or boiling
                    TGRAD_test=T_MAX(im_crit)-Tgamma
                   else
                    print *,"LL invalid"
                    stop
                   endif

                   if (DEBUG_EVAPORATION.eq.1) then
                    print *,"DEBUG_EVAPORATION STATEMENT 1"
                    print *,"i,j,k,x,y,z ",i,j,k, &
                        xsten(0,1),xsten(0,2),xsten(0,SDIM)
                    print *,"im,im_opp,ireverse ",im,im_opp,ireverse
                    print *,"im_source,im_dest ",im_source,im_dest
                    print *,"LL ",LL
                    print *,"project_option=",project_option
                    print *,"Tgamma,TorYgamma_BC ",Tgamma,TorYgamma_BC
                    print *,"im_crit=",im_crit
                    print *,"T_MIN(im_crit) ",T_MIN(im_crit)
                    print *,"T_MAX(im_crit) ",T_MAX(im_crit)
                    print *,"TorY_MIN(im_crit) ",TorY_MIN(im_crit)
                    print *,"TorY_MAX(im_crit) ",TorY_MAX(im_crit)
                    print *,"TGRAD_test=",TGRAD_TEST
                   endif

                    ! local_freezing_model=0 (sharp interface stefan model)
                    ! local_freezing_model=1 (source term model)
                    ! local_freezing_model=2 (hydrate model)
                    ! local_freezing_model=3 (wildfire)
                    ! local_freezing_model=4 (Tanasawa or Schrage)
                    ! local_freezing_model=5 (stefan evaporation/condensation)
                    ! local_freezing_model=6 (Palmore/Desjardins)
                    ! local_freezing_model=7 (Cavitation)
                   if ((local_freezing_model.eq.0).or. & !stefan model
                       ((local_freezing_model.eq.5).and. & !stefan:evap or cond
                        (TGRAD_test.gt.zero)).or. &
                       ((local_freezing_model.eq.6).and. & !Palmore/Desjardins
                        (TGRAD_test.gt.zero))) then 

                    if (im_dest_crit.eq.-1) then
                     im_crit_save=im_crit
                     im_dest_crit=im_dest
                     im_source_crit=im_source
                     im_dest_substrate_crit=im_dest_substrate
                     im_source_substrate_crit=im_source_substrate
                     iten_crit=iten
                     ireverse_crit=ireverse
                     TGRAD_MAX=TGRAD_test
                    else if ((im_dest_crit.ge.1).and. &
                             (im_dest_crit.le.nmat)) then
                     if (TGRAD_test.gt.TGRAD_MAX) then
                      im_crit_save=im_crit
                      im_dest_crit=im_dest
                      im_source_crit=im_source
                      im_dest_substrate_crit=im_dest_substrate
                      im_source_substrate_crit=im_source_substrate
                      iten_crit=iten
                      ireverse_crit=ireverse
                      TGRAD_MAX=TGRAD_test
                     endif
                    else
                     print *,"im_dest_crit invalid"
                     stop
                    endif

                   else if ((local_freezing_model.eq.5).and. &
                            (TGRAD_test.le.zero)) then
                    ! do nothing
                   else if ((local_freezing_model.eq.6).and. &
                            (TGRAD_test.le.zero)) then
                    ! do nothing
                   else
                    print *,"local_freezing_model invalid 16"
                    stop
                   endif

                  else if (TorY_STATUS(im_crit).eq.0) then
                   ! do nothing
                  else
                   print *,"TorY_STATUS(im_crit) invalid"
                   stop
                  endif

                else if ((im_source.ne.im_primary).and. &
                         (im_source_substrate.ne.im_primary).and. &
                         (im_dest.ne.im_primary).and. &
                         (im_dest_substrate.ne.im_primary)) then
                  ! do nothing
                else
                 print *,"LS_center invalid"
                 stop
                endif

               else if ((local_freezing_model.eq.1).or. &
                        (local_freezing_model.eq.2).or. &
                        (local_freezing_model.eq.4)) then !Tanasawa or Schrage
                ! do nothing
               else
                print *,"freezing_model invalid in stefansolver"
                print *,"local_freezing_model= ",local_freezing_model
                print *,"iten,ireverse,nten ",iten,ireverse,nten
                stop
               endif

              else if (start_freezing.eq.0) then
               ! do nothing
              else
               print *,"start_freezing invalid"
               stop
              endif

             else if (LL.eq.zero) then
              ! do nothing
             else
              print *,"LL invalid"
              stop
             endif

            else if ((is_rigid(nmat,im).eq.1).or. &
                     (is_rigid(nmat,im_opp).eq.1)) then
             ! do nothing
            else
             print *,"is_rigid(nmat,im) or is_rigid(nmat,im_opp) invalid"
             stop
            endif

           else if ((Tgamma_STATUS.eq.-1).or.(Tgamma_STATUS.eq.-2)) then
            ! do nothing
           else if (Tgamma_STATUS.eq.0) then
            ! do nothing
           else
            print *,"Tgamma_STATUS invalid"
            stop
           endif

          enddo !ireverse=0,1
         enddo ! im_opp=im+1..nmat
        enddo ! im=1..nmat-1

        if (im_dest_crit.eq.-1) then
         ! do nothing
        else if ((im_dest_crit.ge.1).and. &
                 (im_dest_crit.le.nmat)) then

         SWEPTFACTOR=swept(D_DECL(i,j,k),im_dest_crit) !default:SWEPTFACTOR==1
         if ((SWEPTFACTOR.lt.LSTOL).or. &
             (SWEPTFACTOR.gt.one)) then
          print *,"SWEPTFACTOR INVALID"
          stop
         endif
         over_den=den(D_DECL(i,j,k),1)  ! 1/(rho)
         over_cv=DeDT(D_DECL(i,j,k),1)  ! 1/(rho cv)
         local_vol=vol(D_DECL(i,j,k))
         single_material_den=STATEFAB(D_DECL(i,j,k), &
           (im_primary-1)*num_state_material+1)

         if ((over_den.gt.zero).and. &
             (over_cv.gt.zero).and. &
             (local_vol.gt.zero).and. &
             (single_material_den.gt.zero)) then
          ! do nothing
         else
          print *,"over_den, over_cv, local_vol, or single_mat_den invalid"
          stop
         endif

         original_coeff=one/(dt*SWEPTFACTOR)
         if (project_option.eq.2) then
          original_coeff=original_coeff/over_cv
         else if ((project_option.ge.100).and. &
                  (project_option.lt.100+num_species_var)) then
          original_coeff=original_coeff/over_den
         else
          print *,"project_option invalid"
          stop
         endif

         delta_coeff=zero
         coeff_Tgamma=zero

         im_crit=im_crit_save

         iten=iten_crit
         ireverse=ireverse_crit 

         im_source=im_source_crit
         im_source_substrate=im_source_substrate_crit
         im_dest=im_dest_crit
         im_dest_substrate=im_dest_substrate_crit

         LL=latent_heat(iten+ireverse*nten)
         local_freezing_model=freezing_model(iten+ireverse*nten)
         distribute_from_targ=distribute_from_target(iten+ireverse*nten)

         if (project_option.eq.2) then
          TorYgamma_BC=saturation_temp(iten+ireverse*nten)
         else if ((project_option.ge.100).and. &
                  (project_option.lt.100+num_species_var)) then
          TorYgamma_BC=one
         else
          print *,"project_option invalid"
          stop
         endif


         do dir=1,SDIM
           ii=0
           jj=0
           kk=0
           if (dir.eq.1) then
            ii=1
           else if (dir.eq.2) then
            jj=1
           else if ((dir.eq.SDIM).and.(SDIM.eq.3)) then
            kk=1
           else
            print *,"dir invalid stefansolver"
            stop
           endif 
           do side=-1,1,2
            if (side.eq.-1) then
             ic=i-ii
             jc=j-jj
             kc=k-kk
             iface=i
             jface=j
             kface=k
            else if (side.eq.1) then
             ic=i+ii
             jc=j+jj
             kc=k+kk
             iface=i+ii
             jface=j+jj
             kface=k+kk
            else
             print *,"side invalid"
             stop
            endif

            call gridsten_level(xsten_side,ic,jc,kc,level,nhalf)

            if (dir.eq.1) then
             aface=areax(D_DECL(iface,jface,kface))
            else if (dir.eq.2) then
             aface=areay(D_DECL(iface,jface,kface))
            else if ((dir.eq.3).and.(SDIM.eq.3)) then
             aface=areaz(D_DECL(iface,jface,kface))
            else
             print *,"dir invalid stefansolver 2"
             stop
            endif

            do im_loop=1,nmat
             LS_no_tess(im_loop)=LS(D_DECL(ic,jc,kc),im_loop)
            enddo
            call LS_tessellate(LS_no_tess,LS_side,nmat)

            ! thermal diffusivity==0 where LS changes sign
            ! and latent_heat <> 0.
            at_interface=0
            if ((LS_center(im_source)*LS_side(im_source).le.zero).and. &
                (LS_center(im_dest)*LS_side(im_dest).le.zero)) then
             at_interface=1
             LS1=LS_center(im_source)-LS_center(im_dest)
             LS2=LS_side(im_source)-LS_side(im_dest)
             call get_default_scalar_diffusion(project_option, &
                     LS1,im_source,im_dest, &
                     single_material_den, &
                     heatcoeff)
            else if ((LS_center(im_source_substrate)* &
                      LS_side(im_source_substrate).le.zero).and. &
                     (LS_center(im_dest_substrate)* &
                      LS_side(im_dest_substrate).le.zero)) then
             at_interface=1
             LS1=LS_center(im_source_substrate)-LS_center(im_dest_substrate)
             LS2=LS_side(im_source_substrate)-LS_side(im_dest_substrate)
             call get_default_scalar_diffusion(project_option, &
                     LS1,im_source_substrate,im_dest_substrate, &
                     single_material_den, &
                     heatcoeff)
            endif

            if (at_interface.eq.1) then

              ! cannot do tiling here.
             if (adjust_temperature.eq.-1) then ! modify heatx,heaty,heatz

              if (dir.eq.1) then
               heatx(D_DECL(iface,jface,kface))=zero
              else if (dir.eq.2) then
               heaty(D_DECL(iface,jface,kface))=zero
              else if ((dir.eq.3).and.(SDIM.eq.3)) then
               heatz(D_DECL(iface,jface,kface))=zero
              else
               print *,"dir invalid stefansolver 2B"
               stop
              endif
             
             else if ((adjust_temperature.eq.0).or. & ! modify coeff
                      (adjust_temperature.eq.1)) then ! modify Snew and coeff
              ! do nothing
             else
              print *,"adjust_temperature invalid"
              stop
             endif
        
             if (heatcoeff.lt.zero) then
              print *,"heatcoeff invalid"
              stop
             endif

             if (LS1*LS2.le.zero) then

              if ((LS1.eq.zero).and.(LS2.eq.zero)) then
               theta=half
              else if ((LS1.ne.zero).or.(LS2.ne.zero)) then
               theta=LS1/(LS1-LS2)
              else
               print *,"LS1 or LS2 invalid"
               stop
              endif
              if ((theta.ge.zero).and.(theta.le.one+VOFTOL)) then
               ! do nothing
              else
               print *,"theta invalid"
               stop
              endif
              if (theta.le.theta_cutoff) then
               theta=theta_cutoff
              endif

              do dir_inner=1,SDIM
               x_interface(dir_inner)=theta*xsten_side(0,dir_inner)+ &
                       (one-theta)*xsten(0,dir_inner)
              enddo

              if (project_option.eq.2) then
               tsat_comp=nten+(iten-1)*ncomp_per_tsat+1
              else if ((project_option.ge.100).and. &
                       (project_option.lt.100+num_species_var)) then
               tsat_comp=nten+(iten-1)*ncomp_per_tsat+2
              else
               print *,"project_option invalid"
               stop
              endif

              ngrow_tsat=1
              call interpfab_tsat( &
               i,j,k, &
               ireverse, &
               iten, &
               nten, &
               ntsat, &
               bfact, &
               level, &
               finest_level, &
               dx,xlo, &
               x_interface, &
               tsat_comp, &
               ngrow_tsat, &
               fablo,fabhi, &
               TgammaFAB,DIMS(TgammaFAB), &
               TorYgamma_BC)  ! TorYgamma_BC here is an output

              hx=abs(xsten(0,dir)-xsten(2*side,dir))
              if ((levelrz.eq.3).and.(dir.eq.2)) then
               hx=hx*xsten(side,1)
              endif
              side_coeff=aface*heatcoeff/(theta*hx)
              delta_coeff=delta_coeff+side_coeff
              coeff_Tgamma=coeff_Tgamma+TorYgamma_BC*side_coeff

              if (DEBUG_EVAPORATION.eq.1) then
               print *,"DEBUG_EVAPORATION STATEMENT 2"
               print *,"project_option,i,j,k,dir,side ", &
                       project_option,i,j,k,dir,side
               print *,"im_source,im_dest,TorYgamma_BC ", &
                       im_source,im_dest,TorYgamma_BC
              endif
             endif ! LS1 * LS2 <=0

            else if (at_interface.eq.0) then
             ! do nothing
            else
             print *,"at_interface invalid in FORT_STEFANSOLVER"
             print *,"project_option=",project_option
             print *,"solidheat_flag=",solidheat_flag
             stop
            endif

           enddo ! side=-1,1,2
         enddo ! dir=1..sdim
       
         if (adjust_temperature.eq.-1) then ! modify heatxyz
           ! do nothing
         else if ((adjust_temperature.eq.0).or. & ! modify coeff
                  (adjust_temperature.eq.1)) then ! modify Snew and coeff
 
           if (delta_coeff.gt.zero) then

             ! im_crit is material that dominates center cell.
            if ((im_crit.lt.1).or. &
                (im_crit.gt.nmat)) then
             print *,"im_crit invalid"
             stop
            endif

            delta_coeff=delta_coeff/local_vol
            coeff_Tgamma=coeff_Tgamma/local_vol
   
            if (adjust_temperature.eq.1) then

             TorY_test=TorY_fab(D_DECL(i,j,k),im_crit)
             if (TorY_test.ge.T_or_Y_min_sanity) then
              ! do nothing
             else
              print *,"TorY_test<T_or_Y_min_sanity"
              stop
             endif

             T_adjust=(original_coeff*TorY_test+coeff_Tgamma)/ &
                      (original_coeff+delta_coeff) 

             do im_adjust=1,nmat
              if (project_option.eq.2) then
               tcomp=(SDIM+1)+ &
                (im_adjust-1)*num_state_material+2
              else if ((project_option.ge.100).and. &
                       (project_option.le.100+num_species_var-1)) then
               tcomp=(SDIM+1)+ &
                (im_adjust-1)*num_state_material+3+project_option-100
              else
               print *,"project_option invalid"
               stop
              endif
              Snew(D_DECL(i,j,k),tcomp)=T_adjust
             enddo  ! im_adjust=1..nmat

             coeff(D_DECL(i,j,k))=T_adjust

            else if (adjust_temperature.eq.0) then

             coeff(D_DECL(i,j,k))=original_coeff+delta_coeff

            else
             print *,"adjust_temperature invalid"
             stop
            endif  

           else if (delta_coeff.eq.zero) then
            ! do nothing
           else
            print *,"delta_coeff invalid"
            stop
           endif

         else
           print *,"adjust_temperature invalid"
           stop
         endif

        else 
         print *,"im_dest_crit invalid"
         stop
        endif

       else if (LS_center(im_primary).lt.zero) then
        ! do nothing (vacuum should not be found)
       else
        print *,"LS_center(im_primary) invalid"
        stop
       endif

      enddo ! k=growlo(3),growhi(3)
      enddo ! j=growlo(2),growhi(2)
      enddo ! i=growlo(1),growhi(1)

      return
      end subroutine FORT_STEFANSOLVER

! MEHDI VAHAB HEAT SOURCE
! T^new=T^* + dt A Q/(rho cv V) 
! Q units: J/(m^2 s)
      subroutine FORT_HEATSOURCE_FACE( &
       nmat,nten,nstate, &
       latent_heat, &
       saturation_temp, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx, &
       dt, &
       time, &
       level, &
       finest_level, &
       LS,DIMS(LS),  &
       Snew,DIMS(Snew), & 
       DeDT,DIMS(DeDT), &
       den,DIMS(den), &
       vol,DIMS(vol), &
       heatx,DIMS(heatx), &
       heaty,DIMS(heaty), &
       heatz,DIMS(heatz), &
       areax,DIMS(areax), &
       areay,DIMS(areay), &
       areaz,DIMS(areaz) )

      use probf90_module
      use global_utility_module
      use MOF_routines_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: nmat,nten,nstate
      REAL_T, intent(in) :: latent_heat(2*nten)
      REAL_T, intent(in) :: saturation_temp(2*nten)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: dt
      REAL_T, intent(in) :: time
      INTEGER_T, intent(in) :: DIMDEC(LS)
      INTEGER_T, intent(in) :: DIMDEC(Snew)
      INTEGER_T, intent(in) :: DIMDEC(DeDT)
      INTEGER_T, intent(in) :: DIMDEC(den)
      INTEGER_T, intent(in) :: DIMDEC(vol)
      INTEGER_T, intent(in) :: DIMDEC(heatx)
      INTEGER_T, intent(in) :: DIMDEC(heaty)
      INTEGER_T, intent(in) :: DIMDEC(heatz)
      INTEGER_T, intent(in) :: DIMDEC(areax)
      INTEGER_T, intent(in) :: DIMDEC(areay)
      INTEGER_T, intent(in) :: DIMDEC(areaz)
      REAL_T, intent(in) :: LS(DIMV(LS),nmat*(1+SDIM))
      REAL_T, intent(inout) :: Snew(DIMV(Snew),nstate)
      REAL_T, intent(in) :: DeDT(DIMV(DeDT),nmat+1)  ! 1/(rho cv) (cv=DeDT)
       ! 1/den (i.e. den actually stores 1/den)
      REAL_T, intent(in) :: den(DIMV(den),nmat+1) 
      REAL_T, intent(in) :: vol(DIMV(vol))
      REAL_T, intent(in) :: heatx(DIMV(heatx))
      REAL_T, intent(in) :: heaty(DIMV(heaty))
      REAL_T, intent(in) :: heatz(DIMV(heatz))
      REAL_T, intent(in) :: areax(DIMV(areax))
      REAL_T, intent(in) :: areay(DIMV(areay))
      REAL_T, intent(in) :: areaz(DIMV(areaz))

      INTEGER_T i,j,k
      INTEGER_T heat_dir
      INTEGER_T heat_side
      INTEGER_T ii,jj,kk
      INTEGER_T iface,jface,kface
      INTEGER_T icell,jcell,kcell
      INTEGER_T im
      INTEGER_T im_primary
      INTEGER_T im_primary_cell
      INTEGER_T nten_test
      REAL_T over_den,over_cv,local_vol
      REAL_T aface,hface
      INTEGER_T tcomp
      REAL_T heat_source_term
      REAL_T heat_flux
      REAL_T flux_sign
      REAL_T ls_cell_or_face(nmat)
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf

      if (bfact.lt.1) then
       print *,"bfact invalid68"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base must be 2"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid ratemass nten, nten_test ",nten,nten_test
       stop
      endif
      if (nstate.ne. &
          (SDIM+1)+ &
          nmat*(num_state_material+ngeom_raw)+1) then
       print *,"nstate invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid heat source face"
       stop
      endif
      if (dt.le.zero) then
       print *,"dt invalid"
       stop
      endif
      if (time.lt.zero) then
       print *,"time invalid"
       stop
      endif

      nhalf=3

      call checkbound(fablo,fabhi,DIMS(LS),1,-1,1239)
      call checkbound(fablo,fabhi,DIMS(Snew),1,-1,1240)
      call checkbound(fablo,fabhi,DIMS(DeDT),0,-1,1241)
      call checkbound(fablo,fabhi,DIMS(den),0,-1,1242)
      call checkbound(fablo,fabhi,DIMS(vol),0,-1,1243)

       ! thermal conductivity
      call checkbound(fablo,fabhi,DIMS(heatx),0,0,1244)
      call checkbound(fablo,fabhi,DIMS(heaty),0,1,1245)
      call checkbound(fablo,fabhi,DIMS(heatz),0,SDIM-1,1246)

      call checkbound(fablo,fabhi,DIMS(areax),0,0,1247)
      call checkbound(fablo,fabhi,DIMS(areay),0,1,1248)
      call checkbound(fablo,fabhi,DIMS(areaz),0,SDIM-1,1249)
 
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
       call gridsten_level(xsten,i,j,k,level,nhalf)

       ! MEHDI VAHAB HEAT SOURCE
       ! heat_dir=1,2,3
       ! heat_side=1,2
       if (is_in_probtype_list().eq.1) then
        call SUB_EB_heat_source(time,dt,xsten,nhalf, &
               heat_flux,heat_dir,heat_side)
       else
        heat_flux=zero
       endif

       if (heat_flux.gt.zero) then
 
        do im=1,nmat
         ls_cell_or_face(im)=LS(D_DECL(i,j,k),im)
        enddo
        call get_primary_material(ls_cell_or_face,nmat,im_primary)

        if (is_rigid(nmat,im_primary).eq.1) then

         ii=0
         jj=0
         kk=0
         if (heat_dir.eq.1) then
          ii=1
         else if (heat_dir.eq.2) then
          jj=1
         else if ((heat_dir.eq.3).and.(SDIM.eq.3)) then
          kk=1
         else
          print *,"heat_dir invalid heatsource"
          stop
         endif
         if (heat_side.eq.2) then
          iface=i+ii
          jface=j+jj
          kface=k+kk
          icell=i+ii
          jcell=j+jj
          kcell=k+kk
          flux_sign=one
         else if (heat_side.eq.1) then
          iface=i
          jface=j
          kface=k
          icell=i-ii
          jcell=j-jj
          kcell=k-kk
          flux_sign=-one
         else
          print *,"heat_side invalid"
          stop
         endif

         do im=1,nmat
          ls_cell_or_face(im)=LS(D_DECL(icell,jcell,kcell),im)
         enddo
         call get_primary_material(ls_cell_or_face,nmat,im_primary_cell)

         if (is_rigid(nmat,im_primary_cell).eq.0) then

           ! in: subroutine FORT_HEATSOURCE_FACE
          if (heat_dir.eq.1) then
           aface=areax(D_DECL(iface,jface,kface))
           hface=heatx(D_DECL(iface,jface,kface))
          else if (heat_dir.eq.2) then
           aface=areay(D_DECL(iface,jface,kface))
           hface=heaty(D_DECL(iface,jface,kface))
          else if ((heat_dir.eq.3).and.(SDIM.eq.3)) then
           aface=areaz(D_DECL(iface,jface,kface))
           hface=heatz(D_DECL(iface,jface,kface))
          else
           print *,"heat_dir invalid heatsource 2"
           stop
          endif

          if ((aface.lt.zero).or.(hface.lt.zero)) then
           print *,"aface or hface (thermal conductivity) invalid"
           stop
          endif

          over_den=den(D_DECL(i,j,k),1)
          over_cv=DeDT(D_DECL(i,j,k),1)  ! 1/(rho cv)
          local_vol=vol(D_DECL(i,j,k))

          if ((over_den.le.zero).or.(over_cv.le.zero).or. &
               (local_vol.le.zero)) then
            print *,"over_den, over_cv, or local_vol invalid"
            stop
          endif

          heat_source_term=flux_sign*dt*over_cv*aface*heat_flux/local_vol

          tcomp=(SDIM+1)+2

          Snew(D_DECL(i,j,k),tcomp)= &
           Snew(D_DECL(i,j,k),tcomp)+heat_source_term

         else if (is_rigid(nmat,im_primary_cell).eq.1) then
          ! do nothing
         else
          print *,"is_rigid(nmat,im_primary_cell) invalid"
          stop
         endif 
        else if (is_rigid(nmat,im_primary).eq.0) then
         ! do nothing
        else
         print *,"is_rigid(nmat,im_primary) invalid"
         stop
        endif 
       else if (heat_flux.eq.zero) then
        ! do nothing
       else
        print *,"heat_flux invalid"
        stop
       endif 

      enddo 
      enddo 
      enddo 

      return
      end subroutine FORT_HEATSOURCE_FACE



       ! called from:NavierStokes::init_FSI_GHOST_MAC_MF(int ngrow) 
       ! (in NavierStokes.cpp)
       ! called when "law_of_the_wall=0,1,2"
       ! if nparts==0 => interpolate state cell velocity to MAC grid.
       ! if nparts>0 and law_of_the_wall==0 => interpolate solid cell velocity
       ! to MAC grid.
      subroutine FORT_WALLFUNCTION( &
       data_dir, &
       law_of_the_wall, &
       im_solid_map, &
       level, &
       finest_level, &
       ngrow_law_of_wall, &
       ngrow_distance, &
       nmat, &
       nparts, &
       nparts_ghost, &
       nden, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       xlo,dx, &
       dt, &
       time, &
       LSCP,DIMS(LSCP),  &
       LSFD,DIMS(LSFD),  &
       state,DIMS(state), &
       ufluid,DIMS(ufluid), &
       usolid,DIMS(usolid), &
       ughost,DIMS(ughost), &
       history_dat, &
       DIMS(history_dat), &
       nhistory, &
       visc_coef)
      use probf90_module
      use global_utility_module
      use MOF_routines_module
      use godunov_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: data_dir
      INTEGER_T, intent(in) :: nhistory
      INTEGER_T, intent(in) :: law_of_the_wall
      INTEGER_T, intent(in) :: level,finest_level
      INTEGER_T, intent(in) :: ngrow_law_of_wall
      INTEGER_T, intent(in) :: ngrow_distance
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nparts
      INTEGER_T, intent(in) :: nparts_ghost
      INTEGER_T, intent(in) :: nden
      INTEGER_T, intent(in) :: im_solid_map(nparts_ghost)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in), target :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in), target :: xlo(SDIM)
      REAL_T, intent(in), target :: dx(SDIM)
      REAL_T, intent(in) :: dt
      REAL_T, intent(in) :: time
      REAL_T, intent(in) :: visc_coef
       ! DIMDEC is defined in ArrayLim.H in the BoxLib/Src/C_BaseLib
      INTEGER_T, intent(in) :: DIMDEC(LSCP)
      INTEGER_T, intent(in) :: DIMDEC(LSFD)
      INTEGER_T, intent(in) :: DIMDEC(state)
      INTEGER_T, intent(in) :: DIMDEC(ufluid) ! declare x,y,z dimensions of LS
      INTEGER_T, intent(in) :: DIMDEC(usolid)
      INTEGER_T, intent(in) :: DIMDEC(ughost)
      INTEGER_T, intent(in) :: DIMDEC(history_dat)

        ! LS1,LS2,..,LSn,normal1,normal2,...normal_n 
        ! normal points from negative to positive
        !DIMV(LS)=x,y,z  nmat=num. materials
      !CP=Closest Point
      REAL_T, intent(in), target :: LSCP(DIMV(LSCP),nmat*(SDIM+1)) 
      ! FD=Finite Difference
      REAL_T, intent(in), target :: LSFD(DIMV(LSFD),nmat*SDIM)  
      REAL_T, intent(in), target :: state(DIMV(state),nden)
      REAL_T, intent(in), target :: ufluid(DIMV(ufluid),SDIM+1) ! u,v,w,p
      REAL_T, intent(in), target :: usolid(DIMV(usolid),nparts_ghost*SDIM) 
      REAL_T, intent(out) :: ughost(DIMV(ughost),nparts_ghost*SDIM) 
       ! nhistory=nparts_ghost * (usolid_law_of_wall,uimage,usolid,angle)
      REAL_T, intent(out) :: history_dat(DIMV(history_dat),nhistory) 
      INTEGER_T i,j,k
      INTEGER_T ii,jj,kk
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf
      REAL_T LS_left(nmat)
      REAL_T LS_right(nmat)
      INTEGER_T side_solid,side_image
      INTEGER_T partid
      INTEGER_T im_solid
      INTEGER_T im_fluid
      INTEGER_T im_primary_left
      INTEGER_T im_primary_right
      INTEGER_T im
      INTEGER_T dir
      INTEGER_T isideSOL,jsideSOL,ksideSOL
      INTEGER_T isideFD,jsideFD,ksideFD
      REAL_T, target :: n_raster(SDIM)
      REAL_T, target :: x_projection_raster(SDIM)
      REAL_T, target :: x_image_raster(SDIM)
      REAL_T usolid_law_of_wall(SDIM)
      REAL_T uimage_raster(SDIM)
      REAL_T, target :: usolid_raster(SDIM)
      REAL_T angle_ACT_cell
      INTEGER_T nten
      INTEGER_T nhistory_sub
      type(law_of_wall_parm_type) :: law_of_wall_parm

      nten=( (nmat-1)*(nmat-1)+nmat-1 )/2

      nhalf=3

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in ratemasschange"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (ngrow_distance.eq.4) then
       ! do nothing
      else
       print *,"ngrow_distance invalid"
       stop
      endif
      if (ngrow_law_of_wall.ne.4) then
       print *,"ngrow_law_of_wall invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (nden.ne.nmat*num_state_material) then
       print *,"nden invalid"
       stop
      endif
      if ((nparts.ge.0).and.(nparts.le.nmat)) then 
       ! do nothing
      else
       print *,"nparts invalid FORT_WALLFUNCTION"
       stop
      endif
      if ((nparts_ghost.ge.1).and. &
          (nparts_ghost.le.nmat).and. &
          (nparts_ghost.ge.nparts)) then 
       ! do nothing
      else
       print *,"nparts_ghost invalid FORT_WALLFUNCTION"
       stop
      endif

      if ((nparts_ghost.eq.nparts).or.(nparts_ghost.eq.1)) then
       ! do nothing
      else
       print *,"nparts_ghost invalid"
       stop
      endif

       ! ughost,imgVR,solVR,angle
      nhistory_sub=3*SDIM+1

      if (nhistory.eq.nparts_ghost*nhistory_sub) then
       ! do nothing
      else
       print *,"nhistory invalid"
       stop
      endif

      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif 
      if (time.ge.zero) then
       ! do nothing
      else
       print *,"time invalid"
       stop
      endif 
      if (visc_coef.ge.zero) then
       ! do nothing
      else
       print *,"visc_coef invalid"
       stop
      endif 
      if ((law_of_the_wall.eq.0).or. &
          (law_of_the_wall.eq.1).or. &
          (law_of_the_wall.eq.2)) then
       ! do nothing
      else
       print *,"law_of_the_wall invalid"
       stop
      endif
      if ((data_dir.ge.0).and.(data_dir.le.SDIM-1)) then
       ! do nothing
      else
       print *,"data_dir invalid"
       stop
      endif
      ii=0
      jj=0
      kk=0
      if (data_dir.eq.0) then
       ii=1
      else if (data_dir.eq.1) then
       jj=1
      else if ((data_dir.eq.SDIM-1).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"data_dir invalid"
       stop
      endif

       ! check that LS has ngrow_distance border cells.
       ! the "-1" means that LS is a cell centered variable
       ! instead of a face centered (staggared) variable.
       ! valid values for position type for staggared variables
       ! are 0,1,..,SDIM-1. The last parameter is a "unique" 
       ! caller id that is printed to the screen if the sanity
       ! check fails.
      call checkbound(fablo,fabhi,DIMS(LSCP),ngrow_distance,-1,1252)
      call checkbound(fablo,fabhi,DIMS(LSFD),ngrow_distance,-1,1252)
      call checkbound(fablo,fabhi,DIMS(state),ngrow_law_of_wall,-1,1253)
      call checkbound(fablo,fabhi,DIMS(ufluid),ngrow_law_of_wall,-1,1254)
      call checkbound(fablo,fabhi,DIMS(usolid),ngrow_law_of_wall,-1,1255)
      call checkbound(fablo,fabhi,DIMS(ughost),0,data_dir,1255)
      call checkbound(fablo,fabhi,DIMS(history_dat),0,data_dir,1255)

      law_of_wall_parm%visc_coef=visc_coef
      law_of_wall_parm%time=time
      law_of_wall_parm%dt=dt
      law_of_wall_parm%nmat=nmat
      law_of_wall_parm%nten=nten
      law_of_wall_parm%level=level
      law_of_wall_parm%finest_level=finest_level
      law_of_wall_parm%bfact=bfact
      law_of_wall_parm%dx=>dx
      law_of_wall_parm%xlo=>xlo
      law_of_wall_parm%fablo=>fablo
      law_of_wall_parm%fabhi=>fabhi
      law_of_wall_parm%ngrowfab=ngrow_distance
      law_of_wall_parm%LSCP=>LSCP
      law_of_wall_parm%LSFD=>LSFD
      law_of_wall_parm%state=>state
      law_of_wall_parm%ufluid=>ufluid
      law_of_wall_parm%usolid=>usolid

      law_of_wall_parm%dxmin=dx(1)
      if (dx(2).lt.law_of_wall_parm%dxmin) then
       law_of_wall_parm%dxmin=dx(2)
      endif
      if (dx(SDIM).lt.law_of_wall_parm%dxmin) then
       law_of_wall_parm%dxmin=dx(SDIM)
      endif
       
      call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi, &
              0,data_dir,29) 

       ! A FAB (fortran array box) is tessellated into tiles.
       ! i.e. a single FAB can contain multiple tiles.
       ! BOXLIB stores the FAB.   FAB's store data on a rectangular grid
       ! traverse faces of a given tile.
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       if (nparts.eq.0) then
        if (nparts_ghost.eq.1) then
         do dir=1,SDIM
          ughost(D_DECL(i,j,k),dir)= &
            half*(ufluid(D_DECL(i,j,k),dir)+ &
                  ufluid(D_DECL(i-ii,j-jj,k-kk),dir))
         enddo
        else
         print *,"nparts_ghost invalid"
         stop
        endif
       else if ((nparts.ge.1).and.(nparts.le.nmat)) then

        do partid=1,nparts

         do dir=1,SDIM
          usolid_raster(dir)= &
            half*(usolid(D_DECL(i,j,k),(partid-1)*SDIM+dir)+ &
                  usolid(D_DECL(i-ii,j-jj,k-kk),(partid-1)*SDIM+dir))
          ughost(D_DECL(i,j,k),(partid-1)*SDIM+dir)=usolid_raster(dir)
         enddo  

         if (law_of_the_wall.eq.0) then
          ! do nothing
         else if ((law_of_the_wall.eq.1).or. & !turbulent boundary layer CODY
                  (law_of_the_wall.eq.2)) then !GNBC for contact line ZEYU

          do im=1,nmat
           LS_right(im)=LSCP(D_DECL(i,j,k),im)
           LS_left(im)=LSCP(D_DECL(i-ii,j-jj,k-kk),im)
          enddo
          call get_primary_material(LS_right,nmat,im_primary_right)
          call get_primary_material(LS_left,nmat,im_primary_left)
          if ((im_primary_right.ge.1).and.(im_primary_right.le.nmat).and. &
              (im_primary_left.ge.1).and.(im_primary_left.le.nmat)) then
           ! do nothing
          else
           print *,"im_primary_left or im_primary_right invalid"
           stop
          endif

          im_solid=im_solid_map(partid)+1  ! type integer: material id
          if ((im_solid.ge.1).and.(im_solid.le.nmat)) then
           ! do nothing
          else
           print *,"im_solid invalid FORT_WALLFUNCTION"
           stop
          endif

          side_solid=-1
          side_image=-1

          if (is_lag_part(nmat,im_solid).eq.1) then

           ! law of wall or dynamic contact angle treatment
           ! only for rigid substrates; not flexible substrates.
           if (is_rigid(nmat,im_solid).eq.1) then
            ! im_solid=material id of a rigid solid.
            ! Here, we test if cell center is in the solid.
            ! zero is defined in CONSTANTS.H
            ! CONSTANTS.H is defined in: ./BoxLib/Src/C_BaseLib/CONSTANTS.H

            do dir=1,SDIM
             n_raster(dir)=zero ! points to solid
            enddo      
            if ((LS_right(im_solid).ge.zero).and. &
                (im_primary_right.eq.im_solid)) then
             side_solid=1 ! right side
             isideSOL=i
             jsideSOL=j
             ksideSOL=k
             isideFD=i-ii
             jsideFD=j-jj
             ksideFD=k-kk

             if (is_rigid(nmat,im_primary_left).eq.0) then
              side_image=0  ! left side
              im_fluid=im_primary_left
              n_raster(data_dir+1)=one
              do dir=1,SDIM
               uimage_raster(dir)=ufluid(D_DECL(i-ii,j-jj,k-kk),dir)
              enddo  
             else if (is_rigid(nmat,im_primary_left).eq.1) then
              ! do nothing
             else
              print *,"is_rigid(nmat,im_primary_left) invalid"
              stop
             endif

            else if ((LS_left(im_solid).ge.zero).and. &
                     (im_primary_left.eq.im_solid)) then 
             side_solid=0  ! left side
             isideSOL=i-ii
             jsideSOL=j-jj
             ksideSOL=k-kk
             isideFD=i
             jsideFD=j
             ksideFD=k

             if (is_rigid(nmat,im_primary_right).eq.0) then
              side_image=1 ! right side
              im_fluid=im_primary_right
              n_raster(data_dir+1)=-one
              do dir=1,SDIM
               uimage_raster(dir)=ufluid(D_DECL(i,j,k),dir)
              enddo  
             else if (is_rigid(nmat,im_primary_right).eq.1) then
              ! do nothing
             else
              print *,"is_rigid(nmat,im_primary_right) invalid"
              stop
             endif

            else
             side_solid=-1
             side_image=-1
            endif

            if (((side_solid.eq.0).and.(side_image.eq.1)).or. &
                ((side_solid.eq.1).and.(side_image.eq.0))) then

             ! xsten(0,dir) gives dir'th component of coordinate of the storage
             ! location of cell (i,j,k)
             ! e.g. 1D:
             !      xsten(-2,1)  xsten(0,1)  xsten(2,1)
             !     |     .     |      .     |     .    |
             !          i-1           i          i+1
             ! xsten(-3,1)   xsten(-1,1) xsten(1,1)  xsten(3,1)
             call gridsten_level(xsten,isideSOL,jsideSOL,ksideSOL,level,nhalf)

             do dir=1,SDIM
              x_projection_raster(dir)=xsten(0,dir)
              x_image_raster(dir)=xsten(0,dir)
             enddo
             if ((side_solid.eq.0).and.(side_image.eq.1)) then
              x_projection_raster(data_dir+1)=xsten(1,data_dir+1)
              x_image_raster(data_dir+1)=xsten(2,data_dir+1)
             else if ((side_solid.eq.1).and.(side_image.eq.0)) then
              x_projection_raster(data_dir+1)=xsten(-1,data_dir+1)
              x_image_raster(data_dir+1)=xsten(-2,data_dir+1)
             else
              print *,"side_solid or side_image invalid"
              stop
             endif

             law_of_wall_parm%x_image_raster=>x_image_raster
             law_of_wall_parm%x_projection_raster=>x_projection_raster
             law_of_wall_parm%usolid_raster=>usolid_raster
             law_of_wall_parm%n_raster=>n_raster

              ! call CODY ESTEBEs routine here 
              ! (defined in this file: GODUNOV_3D.F90)
             call getGhostVel( &
               law_of_wall_parm, &
               law_of_the_wall, &
               isideSOL, &
               jsideSOL, &
               ksideSOL, &
               isideFD, &
               jsideFD, &
               ksideFD, &
               side_solid, &
               side_image, &
               data_dir, &
               uimage_raster, &
               usolid_law_of_wall, &
               angle_ACT_cell, &   ! actual contact angle at image point
               im_fluid, &
               im_solid)

               ! solid "ghost" velocity in the solid regions.
             do dir=1,SDIM
              ughost(D_DECL(i,j,k),(partid-1)*SDIM+dir)= &
               usolid_law_of_wall(dir)

              history_dat(D_DECL(i,j,k),(partid-1)*nhistory_sub+dir)= &
                usolid_law_of_wall(dir)
              history_dat(D_DECL(i,j,k),(partid-1)*nhistory_sub+SDIM+dir)= &
                uimage_raster(dir)
              history_dat(D_DECL(i,j,k),(partid-1)*nhistory_sub+2*SDIM+dir)= &
                usolid_raster(dir)
             enddo  ! dir=1..sdim

             history_dat(D_DECL(i,j,k), &
                (partid-1)*nhistory_sub+nhistory_sub)=angle_ACT_cell
           
            else if ((side_solid.eq.-1).or.(side_image.eq.-1)) then
             ! do nothing
            else
             print *,"side_solid or side_image invalid"
             stop
            endif

           else if (is_rigid(nmat,im_solid).eq.0) then
            ! do nothing
           else
            print *,"is_rigid(nmat,im_solid) invalid"
            stop
           endif
          else 
           print *,"is_lag_part(nmat,im_solid) invalid"
           stop
          endif
         else
          print *,"law_of_the_wall invalid"
          stop
         endif
        enddo ! partid=1..nparts
       else
        print *,"nparts invalid"
        stop
       endif
      enddo ! k
      enddo ! j
      enddo ! i

      return
      end subroutine FORT_WALLFUNCTION

        ! called from NavierStokes.cpp
        ! put ns.wall_slip_weight=0.5 for example in the inputs file.
        ! ns.wall_slip_weight=0.0 => do not strengthen the slip BC
        ! ns.wall_slip_weight=1.0 => strongest imposition of slip BC
      subroutine FORT_ASSIMILATE_STATEDATA( &
       isweep, &
       law_of_the_wall, &
       wall_slip_weight, &
       damping_coefficient, &
       im_solid_map, &
       level, &
       finest_level, &
       nstate, &
       nmat, &
       nparts, &
       nparts_ghost, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       xlo,dx, &
       dt, &
       time, & ! cur_time
       LS_state,DIMS(LS_state), &
       state,DIMS(state), &
       macx,DIMS(macx), &
       macy,DIMS(macy), &
       macz,DIMS(macz), &
       ughostx,DIMS(ughostx), &  ! stores the slip velocity
       ughosty,DIMS(ughosty), &  
       ughostz,DIMS(ughostz))
      use probf90_module
      use global_utility_module
      use MOF_routines_module
      use godunov_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: isweep
      INTEGER_T, intent(in) :: law_of_the_wall
      REAL_T, intent(in) :: wall_slip_weight
      INTEGER_T, intent(in) :: level,finest_level
      INTEGER_T, intent(in) :: nstate
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nparts
      INTEGER_T, intent(in) :: nparts_ghost
      REAL_T, intent(in) :: damping_coefficient(nmat)
      INTEGER_T, intent(in), target :: im_solid_map(nparts_ghost)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in), target :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in), target :: xlo(SDIM)
      REAL_T, intent(in), target :: dx(SDIM)
      REAL_T, intent(in) :: dt
      REAL_T, intent(in) :: time
      INTEGER_T, intent(in) :: DIMDEC(LS_state)
      INTEGER_T, intent(in) :: DIMDEC(state)
      INTEGER_T, intent(in) :: DIMDEC(macx)
      INTEGER_T, intent(in) :: DIMDEC(macy)
      INTEGER_T, intent(in) :: DIMDEC(macz)
      INTEGER_T, intent(in) :: DIMDEC(ughostx)
      INTEGER_T, intent(in) :: DIMDEC(ughosty)
      INTEGER_T, intent(in) :: DIMDEC(ughostz)

      REAL_T, intent(inout), target :: LS_state(DIMV(LS_state), &
           nmat*(1+SDIM))
      REAL_T, intent(inout), target :: state(DIMV(state),nstate)
      REAL_T, intent(inout), target :: macx(DIMV(macx))
      REAL_T, intent(inout), target :: macy(DIMV(macy))
      REAL_T, intent(inout), target :: macz(DIMV(macz))
      REAL_T, intent(in), target :: ughostx(DIMV(ughostx),nparts_ghost*SDIM) 
      REAL_T, intent(in), target :: ughosty(DIMV(ughosty),nparts_ghost*SDIM) 
      REAL_T, intent(in), target :: ughostz(DIMV(ughostz),nparts_ghost*SDIM) 
      INTEGER_T i,j,k
      INTEGER_T iface,jface,kface
      INTEGER_T icell,jcell,kcell
      INTEGER_T ileft,jleft,kleft
      INTEGER_T iright,jright,kright
      INTEGER_T i_nbr,j_nbr,k_nbr
      INTEGER_T ii,jj,kk
      INTEGER_T iii,jjj,kkk
      INTEGER_T im
      INTEGER_T im_primary
      INTEGER_T im_stencil
      INTEGER_T im_stencil_left
      INTEGER_T im_stencil_right
      REAL_T :: local_damping
      REAL_T :: LS_local(nmat)
      REAL_T :: LS_stencil(nmat)
      REAL_T :: LS_LEFT(nmat)
      REAL_T :: LS_RIGHT(nmat)
      REAL_T, target :: xsten(-3:3,SDIM)
      INTEGER_T nhalf
      INTEGER_T nstate_test
      type(assimilate_parm_type) :: assimilate_parm
      type(assimilate_out_parm_type) :: assimilate_out_parm
      INTEGER_T cell_flag
      INTEGER_T veldir
      INTEGER_T dir
      INTEGER_T dirtan
      INTEGER_T side
      INTEGER_T side_nbr
      INTEGER_T partid
      REAL_T velsum(SDIM)
      REAL_T wtsum
      REAL_T velface

      nhalf=3

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in ratemasschange"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      nstate_test=(SDIM+1)+ &
        nmat*(num_state_material+ngeom_raw)+1
      if (nstate.ne.nstate_test) then
       print *,"nstate invalid in GODUNOV_3D.F90 "
       print *,"nstate=",nstate
       print *,"nstate_test=",nstate_test
       stop
      endif
      if ((nparts.ge.0).and.(nparts.le.nmat)) then 
       ! do nothing
      else
       print *,"nparts invalid FORT_ASSIMILATE_STATEDATA"
       stop
      endif
      if ((nparts_ghost.ge.1).and. &
          (nparts_ghost.le.nmat).and. &
          (nparts_ghost.ge.nparts)) then 
       ! do nothing
      else
       print *,"nparts_ghost invalid FORT_WALLFUNCTION"
       stop
      endif

      if ((nparts_ghost.eq.nparts).or.(nparts_ghost.eq.1)) then
       ! do nothing
      else
       print *,"nparts_ghost invalid"
       stop
      endif

      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif 
      if (time.ge.zero) then
       ! do nothing
      else
       print *,"time invalid"
       stop
      endif 
      if ((law_of_the_wall.eq.0).or. &
          (law_of_the_wall.eq.1).or. &
          (law_of_the_wall.eq.2)) then
       ! do nothing
      else
       print *,"law_of_the_wall invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(LS_state),1,-1,1253)
      call checkbound(fablo,fabhi,DIMS(state),1,-1,1253)
      call checkbound(fablo,fabhi,DIMS(macx),0,0,1253)
      call checkbound(fablo,fabhi,DIMS(macy),0,1,1253)
      call checkbound(fablo,fabhi,DIMS(macz),0,SDIM-1,1253)
      call checkbound(fablo,fabhi,DIMS(ughostx),0,0,1255)
      call checkbound(fablo,fabhi,DIMS(ughosty),0,1,1255)
      call checkbound(fablo,fabhi,DIMS(ughostz),0,SDIM-1,1255)

      assimilate_parm%time=time
      assimilate_parm%dt=dt
      assimilate_parm%nhalf=nhalf
      assimilate_parm%nstate=nstate
      assimilate_parm%nmat=nmat
      assimilate_parm%nparts=nparts
      assimilate_parm%nparts_ghost=nparts_ghost
      assimilate_parm%im_solid_map=>im_solid_map
      assimilate_parm%level=level
      assimilate_parm%finest_level=finest_level
      assimilate_parm%bfact=bfact
      assimilate_parm%dx=>dx
      assimilate_parm%xlo=>xlo
      assimilate_parm%fablo=>fablo
      assimilate_parm%fabhi=>fabhi
      assimilate_parm%ughostx=>ughostx
      assimilate_parm%ughosty=>ughosty
      assimilate_parm%ughostz=>ughostz

      assimilate_parm%dxmin=dx(1)
      if (dx(2).lt.assimilate_parm%dxmin) then
       assimilate_parm%dxmin=dx(2)
      endif
      if (dx(SDIM).lt.assimilate_parm%dxmin) then
       assimilate_parm%dxmin=dx(SDIM)
      endif
       
      assimilate_out_parm%state=>state
      assimilate_out_parm%macx=>macx
      assimilate_out_parm%macy=>macy
      assimilate_out_parm%macz=>macz

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      cell_flag=-1

      if ((wall_slip_weight.ge.zero).and. &
          (wall_slip_weight.le.one)) then
       ! do nothing
      else
       print *,"wall_slip_weight invalid"
       stop
      endif

      if (isweep.eq.0) then

       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)

        call gridsten_level(xsten,i,j,k,level,nhalf)

        assimilate_parm%xsten=>xsten
        if (is_in_probtype_list().eq.1) then
         call SUB_ASSIMILATE(assimilate_parm,assimilate_out_parm, &
          i,j,k,cell_flag)
        else
         ! do nothing
        endif

         ! check if a fluid cell neighbors a solid cell
        do im=1,nmat
         LS_local(im)=LS_state(D_DECL(i,j,k),im)
        enddo
        ! first checks the rigid materials for a positive LS; if none
        ! exist, then the subroutine checks the fluid materials.
        call get_primary_material(LS_local,nmat,im_primary) 

        if (is_rigid(nmat,im_primary).eq.0) then
         ! check all neighbors in "star stencil" for solid cells.
         ! for each solid cell, update the present cell center velocity
         ! with a weighted average of the slip velocity and the present
         ! velocity.
         do veldir=1,SDIM
          velsum(veldir)=zero
         enddo
         wtsum=zero

         do dir=1,SDIM
          ii=0
          jj=0
          kk=0
          if (dir.eq.1) then
           ii=1
          else if (dir.eq.2) then
           jj=1
          else if ((dir.eq.3).and.(SDIM.eq.3)) then
           kk=1
          else
           print *,"dir invalid"
           stop
          endif
          do side=1,2
           if (side.eq.1) then
            iface=i
            jface=j
            kface=k
            icell=i-ii
            jcell=j-jj
            kcell=k-kk
           else if (side.eq.2) then
            iface=i+ii
            jface=j+jj
            kface=k+kk
            icell=i+ii
            jcell=j+jj
            kcell=k+kk
           else
            print *,"side invalid"
            stop
           endif

           do im=1,nmat
            LS_stencil(im)=LS_state(D_DECL(icell,jcell,kcell),im)
           enddo
           call get_primary_material(LS_stencil,nmat,im_stencil) 
           if (is_rigid(nmat,im_stencil).eq.0) then
            ! do nothing
           else if (is_rigid(nmat,im_stencil).eq.1) then
            partid=1
            do while ((im_solid_map(partid)+1.ne.im_stencil).and. &
                      (partid.lt.nparts_ghost))
             partid=partid+1
            enddo
            if (im_solid_map(partid)+1.eq.im_stencil) then
             do veldir=1,SDIM
              if (dir.eq.1) then
               velface=ughostx(D_DECL(iface,jface,kface), &
                       (partid-1)*SDIM+veldir)
              else if (dir.eq.2) then
               velface=ughosty(D_DECL(iface,jface,kface), &
                       (partid-1)*SDIM+veldir)
              else if ((dir.eq.3).and.(SDIM.eq.3)) then
               velface=ughostz(D_DECL(iface,jface,kface), &
                       (partid-1)*SDIM+veldir)
              else
               print *,"dir invalid"
               stop
              endif

              velsum(veldir)=velsum(veldir)+velface
             enddo ! veldir=1..sdim
             wtsum=wtsum+one
            else
             print *,"im_solid_map(partid) invalid"
             stop
            endif
           else
            print *,"is_rigid(nmat,im_stencil) invalid"
            stop
           endif 
          enddo ! side=1,2
         enddo ! dir=1..SDIM
         if (wtsum.eq.zero) then
          ! do nothing
         else if (wtsum.gt.zero) then
          do veldir=1,SDIM
           state(D_DECL(i,j,k),veldir)= &
                  (one-wall_slip_weight)*state(D_DECL(i,j,k),veldir)+ &
                  wall_slip_weight*velsum(veldir)/wtsum
          enddo
         else
          print *,"wtsum invalid"
          stop
         endif 

        else if (is_rigid(nmat,im_primary).eq.1) then
         ! do nothing
        else
         print *,"is_rigid(nmat,im_primary) invalid"
         stop
        endif

         ! cell centered velocity
        if (damping_coefficient(im_primary).gt.zero) then
         ! v_t = -mu v =>  v^{n+1} - v^{n} = -mu dt v^{n+1}
         ! v^{n+1}=v^{n}/(1+mu dt)
         do veldir=1,SDIM
          state(D_DECL(i,j,k),veldir)=state(D_DECL(i,j,k),veldir)/ &
                  (one+damping_coefficient(im_primary)*dt)
         enddo
        else if (damping_coefficient(im_primary).eq.zero) then
         ! do nothing
        else
         print *,"damping_coefficient(im_primary) invalid"
         stop
        endif
       enddo ! k
       enddo ! j
       enddo ! i

      else if (isweep.eq.1) then

       do cell_flag=0,SDIM-1

        ii=0
        jj=0
        kk=0
        if (cell_flag.eq.0) then
         ii=1
        else if (cell_flag.eq.1) then
         jj=1
        else if ((cell_flag.eq.2).and.(SDIM.eq.3)) then
         kk=1
        else
         print *,"cell_flag invalid"
         stop
        endif
        call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
             growlo,growhi,0,cell_flag,30) 

        do i=growlo(1),growhi(1)
        do j=growlo(2),growhi(2)
        do k=growlo(3),growhi(3)

          !cell_flag=0..sdim-1
         call gridstenMAC_level(xsten,i,j,k,level,nhalf,cell_flag,33)

         assimilate_parm%xsten=>xsten
         if (is_in_probtype_list().eq.1) then
          call SUB_ASSIMILATE(assimilate_parm,assimilate_out_parm, &
           i,j,k,cell_flag)
         else
          ! do nothing
         endif

          ! check that both adjoining cells are fluid cells
         ileft=i-ii
         jleft=j-jj
         kleft=k-kk
         iright=i
         jright=j
         kright=k
         do im=1,nmat
          LS_LEFT(im)=LS_state(D_DECL(ileft,jleft,kleft),im)
         enddo
         call get_primary_material(LS_LEFT,nmat,im_stencil_left)
         do im=1,nmat
          LS_RIGHT(im)=LS_state(D_DECL(iright,jright,kright),im)
         enddo
         call get_primary_material(LS_RIGHT,nmat,im_stencil_right)

         if (is_rigid(nmat,im_stencil_left).eq.0) then
          if (is_rigid(nmat,im_stencil_right).eq.0) then
           ! check if any neighbor solid faces oriented perpendicular to the
           ! given face.
           !   side_nbr=1   side_nbr=2
           !     --------------------
           !     |    s?   |  s?    |
           !     |         |        |
           !     --------------------
           !     |    f   f|f  f    |
           !     |        f|f       |
           !     --------------------
           !     |    s?   |  s?    |
           !     |         |        |
           !     --------------------
           ! loop through the "s?" cells, wtsum will contain the count
           wtsum=zero
           do dirtan=1,SDIM
            if (dirtan.ne.cell_flag+1) then
             do side_nbr=1,2
              if (side_nbr.eq.1) then
               i_nbr=ileft
               j_nbr=jleft
               k_nbr=kleft
              else if (side_nbr.eq.2) then
               i_nbr=iright
               j_nbr=jright
               k_nbr=kright
              else
               print *,"side_nbr invalid"
               stop
              endif
              iii=0
              jjj=0
              kkk=0
              if (dirtan.eq.1) then
               iii=1
              else if (dirtan.eq.2) then
               jjj=1
              else if ((dirtan.eq.3).and.(SDIM.eq.3)) then
               kkk=1
              else
               print *,"dirtan invalid"
               stop
              endif
            
              do side=1,2
               if (side.eq.1) then
                icell=i_nbr-iii
                jcell=j_nbr-jjj
                kcell=k_nbr-kkk
               else if (side.eq.2) then
                icell=i_nbr+iii
                jcell=j_nbr+jjj
                kcell=k_nbr+kkk
               else
                print *,"side invalid"
                stop
               endif

               do im=1,nmat
                LS_stencil(im)=LS_state(D_DECL(icell,jcell,kcell),im)
               enddo
               call get_primary_material(LS_stencil,nmat,im_stencil) 
               if (is_rigid(nmat,im_stencil).eq.0) then
                ! do nothing
               else if (is_rigid(nmat,im_stencil).eq.1) then
                partid=1
                do while ((im_solid_map(partid)+1.ne.im_stencil).and. &
                          (partid.lt.nparts_ghost))
                 partid=partid+1
                enddo
                if (im_solid_map(partid)+1.eq.im_stencil) then
                 wtsum=wtsum+one
                else
                 print *,"im_solid_map(partid) invalid"
                 stop
                endif
               else
                print *,"is_rigid(nmat,im_stencil) invalid"
                stop
               endif 
              enddo ! side=1,2
             enddo ! side_nbr=1..2
            else if (dirtan.eq.cell_flag+1) then
             ! do nothing
            else
             print *,"dirtan invalid"
             stop
            endif
           enddo ! dirtan=1..SDIM
           if (wtsum.eq.zero) then
            ! do nothing
           else if (wtsum.gt.zero) then
            do veldir=1,SDIM

             velsum(veldir)= &
                half*(state(D_DECL(ileft,jleft,kleft),veldir)+ &
                      state(D_DECL(iright,jright,kright),veldir))

             if (cell_flag+1.eq.veldir) then
              if (veldir.eq.1) then
               macx(D_DECL(i,j,k))= &
                   (one-wall_slip_weight)*macx(D_DECL(i,j,k))+ &
                   wall_slip_weight*velsum(veldir)
              else if (veldir.eq.2) then
               macy(D_DECL(i,j,k))= &
                   (one-wall_slip_weight)*macy(D_DECL(i,j,k))+ &
                   wall_slip_weight*velsum(veldir)
              else if ((veldir.eq.3).and.(SDIM.eq.3)) then
               macz(D_DECL(i,j,k))= &
                   (one-wall_slip_weight)*macz(D_DECL(i,j,k))+ &
                   wall_slip_weight*velsum(veldir)
              else
               print *,"veldir invalid"
               stop
              endif
             else if ((cell_flag+1.ge.1).and. &
                      (cell_flag+1.le.SDIM)) then
              ! do nothing
             else
              print *,"cell_flag invalid"
              stop
             endif
            enddo ! veldir=1..sdim
           else
            print *,"wtsum invalid"
            stop
           endif 
          else if (is_rigid(nmat,im_stencil_right).eq.1) then
           ! do nothing
          else
           print *,"is_rigid(nmat,im_stencil_right) (right) invalid"
           stop
          endif
         else if (is_rigid(nmat,im_stencil_left).eq.1) then
          ! do nothing
         else
          print *,"is_rigid(nmat,im_stencil_left) (left) invalid"
          stop
         endif
         do im=1,nmat
          LS_stencil(im)=half*(LS_LEFT(im)+LS_RIGHT(im))
         enddo
         call get_primary_material(LS_stencil,nmat,im_stencil)

         local_damping=damping_coefficient(im_stencil)
         if (local_damping.gt.zero) then
          if (cell_flag.eq.0) then
           macx(D_DECL(i,j,k))=macx(D_DECL(i,j,k))/(one+dt*local_damping)
          else if (cell_flag.eq.1) then
           macy(D_DECL(i,j,k))=macy(D_DECL(i,j,k))/(one+dt*local_damping)
          else if ((cell_flag.eq.2).and.(SDIM.eq.3)) then
           macz(D_DECL(i,j,k))=macz(D_DECL(i,j,k))/(one+dt*local_damping)
          else
           print *,"cell_flag invalid"
           stop
          endif
         else if (local_damping.eq.zero) then
          ! do nothing
         else
          print *,"local_damping invalid"
          stop
         endif

        enddo ! k
        enddo ! j
        enddo ! i

       enddo ! cell_flag=0...sdim-1

      else
       print *,"isweep invalid"
       stop
      endif


      return
      end subroutine FORT_ASSIMILATE_STATEDATA


        ! recon:
        ! vof,ref centroid,order,slope,intercept  x nmat
        ! icemask_index component initialized to 1 in "init_physics_vars"
        ! if nmat=2, nten=1
        ! if nmat=3, nten=3    12 13 23
        ! if nmat=4, nten=6    12 13 14 23 24 34
      subroutine FORT_INITJUMPTERM( &
       mdotplus, &
       mdotminus, &
       mdotcount, &
       ngrow_expansion, &
       time, &
       level,finest_level, &
       nmat,nten, &
       latent_heat, &
       saturation_temp, &
       freezing_model, &
       distribute_from_target, &
       constant_volume_mdot, &
       constant_density_all_time, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx,dt, &
       maskcov,DIMS(maskcov), &
       JUMPFAB,DIMS(JUMPFAB), &
       mdot,DIMS(mdot), &
       LSnew,DIMS(LSnew) )
      use probf90_module
      use global_utility_module
      use MOF_routines_module

      IMPLICIT NONE

      REAL_T, intent(inout) :: mdotplus
      REAL_T, intent(inout) :: mdotminus
      REAL_T, intent(inout) :: mdotcount
      INTEGER_T, intent(in) :: ngrow_expansion
      REAL_T, intent(in) :: time
      INTEGER_T, intent(in) :: level,finest_level
      INTEGER_T, intent(in) :: nmat,nten
      REAL_T, intent(in) :: latent_heat(2*nten)
      REAL_T, intent(in) :: saturation_temp(2*nten)
      INTEGER_T, intent(in) :: freezing_model(2*nten)
      INTEGER_T, intent(in) :: distribute_from_target(2*nten)
      INTEGER_T, intent(in) :: constant_volume_mdot(2*nten)
      INTEGER_T, intent(in) :: constant_density_all_time(nmat)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: dt
      INTEGER_T, intent(in) :: DIMDEC(maskcov)
      INTEGER_T, intent(in) :: DIMDEC(JUMPFAB)
      INTEGER_T, intent(in) :: DIMDEC(mdot)
      INTEGER_T, intent(in) :: DIMDEC(LSnew)
      REAL_T, intent(in) :: maskcov(DIMV(maskcov))
      REAL_T, intent(in) :: JUMPFAB(DIMV(JUMPFAB),2*nten)
      REAL_T, intent(inout) :: mdot(DIMV(mdot))
      REAL_T, intent(in) :: LSnew(DIMV(LSnew),nmat)

      INTEGER_T i,j,k
      INTEGER_T im,im_opp,ireverse,iten
      INTEGER_T iten_shift
      INTEGER_T im_source,im_dest
      INTEGER_T nten_test
      REAL_T jump_strength

      REAL_T LL
      REAL_T divu_material
      INTEGER_T local_mask


      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif

      if (ngrow_expansion.ne.2) then
       print *,"ngrow_expansion invalid"
       stop
      endif
      if (time.lt.zero) then
       print *,"time invalid"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in convert_material"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid ratemass nten, nten_test ",nten,nten_test
       stop
      endif
      if (dt.le.zero) then
       print *,"dt invalid"
       stop
      endif

      do im=1,nmat-1
       do im_opp=im+1,nmat
        do ireverse=0,1
         call get_iten(im,im_opp,iten,nmat)
         iten_shift=iten+ireverse*nten
         if ((freezing_model(iten_shift).lt.0).or. &
             (freezing_model(iten_shift).gt.7)) then
          print *,"freezing_model invalid init jump term"
          print *,"iten,ireverse,nten ",iten,ireverse,nten
          stop
         endif
         if ((distribute_from_target(iten_shift).lt.0).or. &
             (distribute_from_target(iten_shift).gt.1)) then
          print *,"distribute_from_target invalid init jump term"
          print *,"iten,ireverse,nten ",iten,ireverse,nten
          stop
         endif
         if (constant_volume_mdot(iten_shift).eq.0) then 
          ! do nothing
         else if (constant_volume_mdot(iten_shift).eq.1) then
          ! distribute -sum mdot to the source:
         else if (constant_volume_mdot(iten_shift).eq.-1) then
          ! distribute -sum mdot to the dest:
         else
          print *,"constant_volume_mdot(iten_shift) invalid"
          stop
         endif
        enddo ! ireverse
       enddo ! im_opp
      enddo ! im

      call checkbound(fablo,fabhi,DIMS(maskcov),1,-1,1264)
      call checkbound(fablo,fabhi,DIMS(JUMPFAB),ngrow_expansion,-1,1268)
      call checkbound(fablo,fabhi,DIMS(mdot),0,-1,1269)
      call checkbound(fablo,fabhi,DIMS(LSnew),1,-1,1270)
 
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       local_mask=NINT(maskcov(D_DECL(i,j,k)))

       if (local_mask.eq.1) then

        do im=1,nmat-1
         do im_opp=im+1,nmat
          do ireverse=0,1
           if ((im.gt.nmat).or.(im_opp.gt.nmat)) then
            print *,"im or im_opp bust 10"
            stop
           endif

           call get_iten(im,im_opp,iten,nmat)
           iten_shift=iten+ireverse*nten

           LL=latent_heat(iten_shift)

           jump_strength=JUMPFAB(D_DECL(i,j,k),iten_shift)
  
           if ((is_rigid(nmat,im).eq.1).or. &
               (is_rigid(nmat,im_opp).eq.1)) then 
            ! do nothing
           else if (LL.ne.zero) then
            if (ireverse.eq.0) then
             im_source=im
             im_dest=im_opp
            else if (ireverse.eq.1) then
             im_source=im_opp
             im_dest=im
            else
             print *,"ireverse invalid"
             stop
            endif

              ! jump_strength units: cm^3/s^2
            divu_material=jump_strength

            if (divu_material.gt.zero) then
             mdotplus=mdotplus+divu_material
             mdotcount=mdotcount+one
            else if (divu_material.lt.zero) then
             mdotminus=mdotminus+divu_material
             mdotcount=mdotcount+one
            else if (divu_material.eq.zero) then
             ! do nothing
            else
             print *,"divu_material bust"
             stop
            endif
            
            mdot(D_DECL(i,j,k))=mdot(D_DECL(i,j,k))+divu_material
           else if (LL.eq.zero) then
            ! do nothing
           else
            print *,"LL bust"
            stop
           endif  
          enddo ! ireverse=0...1
         enddo ! im_opp=im+1...nmat
        enddo ! im=1...nmat-1

       else if (local_mask.eq.0) then
        ! do nothing
       else
        print *,"local_mask invalid"
        stop
       endif

      enddo ! k
      enddo ! j
      enddo ! i

      return
      end subroutine FORT_INITJUMPTERM

        ! recon:
        ! vof,ref centroid,order,slope,intercept  x nmat
        ! icemask_index component initialized to 1 in "init_physics_vars"
        ! if nmat=2, nten=1
        ! if nmat=3, nten=3    12 13 23
        ! if nmat=4, nten=6    12 13 14 23 24 34
        ! called from NavierStokes.cpp: NavierStokes::level_init_icemask()
        !   which is called from
        !     NavierStokes::make_physics_varsALL
      subroutine FORT_INIT_ICEMASK( &
       time, &
       facecut_index, &
       icefacecut_index, &
       icemask_index, &
       massface_index, &
       vofface_index, &
       ncphys, &
       level,finest_level, &
       nmat,nten, &
       latent_heat, &
       saturation_temp, &
       freezing_model, &
       distribute_from_target, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx,dt, &
       maskcov,DIMS(maskcov), &
       xface,DIMS(xface), &
       yface,DIMS(yface), &
       zface,DIMS(zface), &
       LSnew,DIMS(LSnew), &
       recon,DIMS(recon) )
      use probf90_module
      use global_utility_module
      use MOF_routines_module

      IMPLICIT NONE

      REAL_T time
      INTEGER_T facecut_index
      INTEGER_T icefacecut_index
      INTEGER_T icemask_index
      INTEGER_T massface_index
      INTEGER_T vofface_index
      INTEGER_T ncphys
      INTEGER_T level,finest_level
      INTEGER_T nmat,nten
      REAL_T latent_heat(2*nten)
      REAL_T saturation_temp(2*nten)
      INTEGER_T freezing_model(2*nten)
      INTEGER_T distribute_from_target(2*nten)
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growloMAC(3),growhiMAC(3)
      INTEGER_T bfact
      REAL_T xlo(SDIM)
      REAL_T dx(SDIM)
      REAL_T dt
      INTEGER_T DIMDEC(maskcov)
      INTEGER_T DIMDEC(xface)
      INTEGER_T DIMDEC(yface)
      INTEGER_T DIMDEC(zface)
      INTEGER_T DIMDEC(LSnew)
      INTEGER_T DIMDEC(recon)
      REAL_T maskcov(DIMV(maskcov))
      REAL_T xface(DIMV(xface),ncphys)
      REAL_T yface(DIMV(yface),ncphys)
      REAL_T zface(DIMV(zface),ncphys)
      REAL_T LSnew(DIMV(LSnew),nmat)
      REAL_T recon(DIMV(recon),nmat*ngeom_recon)
      INTEGER_T i,j,k
      INTEGER_T dir
      INTEGER_T dir2
      INTEGER_T ii,jj,kk
      INTEGER_T ireverse
      INTEGER_T iten
      INTEGER_T im,im_opp
      INTEGER_T im_left,im_opp_left
      INTEGER_T im_right,im_opp_right
      INTEGER_T ireverse_left,ireverse_right
      INTEGER_T nten_test
      REAL_T LSleft(nmat)
      REAL_T LSright(nmat)

      REAL_T ice_test,cut_test
      REAL_T icemask_left
      REAL_T icemask_right
      REAL_T icefacecut_left
      REAL_T icefacecut_right
      REAL_T icemask
      REAL_T icefacecut
      REAL_T xstenMAC(-3:3,SDIM)
      REAL_T xmac(SDIM)
      INTEGER_T nhalf
      INTEGER_T local_mask_right
      INTEGER_T local_mask_left

      nhalf=3

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif

      if (time.lt.zero) then
       print *,"time invalid"
       stop
      endif

      if (FORT_MUSHY_THICK.lt.one) then
       print *,"FORT_MUSHY_THICK too small"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in init_icemask"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      if (vofface_index.ne.massface_index+2*nmat) then
       print *,"face_index bust 1"
       stop
      endif
      if (icemask_index.ne.5) then
       print *,"icemask_index bust"
       stop
      endif
      if (facecut_index.ne.3) then
       print *,"facecut_index bust"
       stop
      endif
      if (icefacecut_index.ne.4) then
       print *,"icefacecut_index bust"
       stop
      endif
      if (ncphys.ne.vofface_index+2*nmat) then
       print *,"ncphys invalid"
       stop
      endif

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid ratemass nten, nten_test ",nten,nten_test
       stop
      endif
      if (dt.le.zero) then
       print *,"dt invalid"
       stop
      endif

      do im=1,nmat-1
       do im_opp=im+1,nmat
        do ireverse=0,1
         call get_iten(im,im_opp,iten,nmat)
         if ((freezing_model(iten+ireverse*nten).lt.0).or. &
             (freezing_model(iten+ireverse*nten).gt.7)) then
          print *,"freezing_model invalid init ice mask"
          print *,"iten,ireverse,nten ",iten,ireverse,nten
          stop
         endif
         if ((distribute_from_target(iten+ireverse*nten).lt.0).or. &
             (distribute_from_target(iten+ireverse*nten).gt.1)) then
          print *,"distribute_from_target invalid init ice mask"
          print *,"iten,ireverse,nten ",iten,ireverse,nten
          stop
         endif
        enddo ! ireverse
       enddo ! im_opp
      enddo ! im

      call checkbound(fablo,fabhi,DIMS(maskcov),1,-1,1264)
      call checkbound(fablo,fabhi,DIMS(xface),0,0,213)
      call checkbound(fablo,fabhi,DIMS(yface),0,1,214)
      call checkbound(fablo,fabhi,DIMS(zface),0,SDIM-1,215)
      call checkbound(fablo,fabhi,DIMS(LSnew),1,-1,1270)
      call checkbound(fablo,fabhi,DIMS(recon),1,-1,1271)
 
      do dir=0,SDIM-1
       ii=0
       jj=0
       kk=0
       if (dir.eq.0) then
        ii=1
       else if (dir.eq.1) then
        jj=1
       else if ((dir.eq.2).and.(SDIM.eq.3)) then
        kk=1
       else
        print *,"dir invalid init ice mask"
        stop
       endif

       call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
        growloMAC,growhiMAC,0,dir,31) 
       do i=growloMAC(1),growhiMAC(1)
       do j=growloMAC(2),growhiMAC(2)
       do k=growloMAC(3),growhiMAC(3)

        local_mask_right=NINT(maskcov(D_DECL(i,j,k))) 
        local_mask_left=NINT(maskcov(D_DECL(i-ii,j-jj,k-kk))) 

         ! check if this is an uncovered face.
        if ((local_mask_right.eq.1).or. &
            (local_mask_left.eq.1)) then
         
         call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,dir,34)
         do dir2=1,SDIM
          xmac(dir2)=xstenMAC(0,dir2)
         enddo 
         do im=1,nmat
          LSleft(im)=LSnew(D_DECL(i-ii,j-jj,k-kk),im)
          LSright(im)=LSnew(D_DECL(i,j,k),im)
         enddo

          ! get_icemask defined in PROB.F90
          ! in: FORT_INIT_ICEMASK
         call get_icemask( &
          xmac, &
          time, &
          dx,bfact, &
          icemask_left, &  ! 0 or 1
          icefacecut_left, & ! 0<=f<=1
          im_left, &
          im_opp_left, &
          ireverse_left, &
          LSleft, &
          latent_heat, &
          distribute_from_target, &
          nmat,nten)

          ! get_icemask defined in PROB.F90
          ! in: FORT_INIT_ICEMASK
         call get_icemask( &
          xmac, &
          time, &
          dx,bfact, &
          icemask_right, &  ! 0 or 1
          icefacecut_right, & ! 0<=f<=1
          im_right, &
          im_opp_right, &
          ireverse_right, &
          LSright, &
          latent_heat, &
          distribute_from_target, &
          nmat,nten)

         icemask=min(icemask_left,icemask_right)

         if ((icefacecut_left.ge.zero).and. &
             (icefacecut_right.ge.zero).and. &
             (icefacecut_left.le.one).and. &
             (icefacecut_right.le.one)) then

          if ((icefacecut_left.le.two*ICEFACECUT_EPS).and. &
              (icefacecut_right.le.two*ICEFACECUT_EPS)) then
           icefacecut=zero
          else
           icefacecut=min(icefacecut_left,icefacecut_right)
          endif

         else
          print *,"icefacecut_left or icefacecut_right invalid"
          print *,"icefacecut_left ",icefacecut_left
          print *,"icefacecut_right ",icefacecut_right
          stop
         endif

         if (icemask.eq.zero) then
          ! do nothing
         else if (icemask.eq.one) then
          ! do nothing
         else
          print *,"icemask invalid icemask=",icemask
          print *,"icefacecut=",icefacecut
          stop
         endif
         if ((icefacecut.ge.zero).and.(icefacecut.le.one)) then
          ! do nothing
         else
          print *,"icefacecut invalid icefacecut=",icefacecut 
          stop
         endif

         if (dir.eq.0) then
          ice_test=xface(D_DECL(i,j,k),icemask_index+1)
          cut_test=xface(D_DECL(i,j,k),icefacecut_index+1)
         else if (dir.eq.1) then
          ice_test=yface(D_DECL(i,j,k),icemask_index+1)
          cut_test=yface(D_DECL(i,j,k),icefacecut_index+1)
         else if ((dir.eq.2).and.(SDIM.eq.3)) then
          ice_test=zface(D_DECL(i,j,k),icemask_index+1)
          cut_test=zface(D_DECL(i,j,k),icefacecut_index+1)
         else
          print *,"dir invalid init_icemask 2"
          stop
         endif

         if (ice_test.eq.zero) then
          ! do nothing
         else if (ice_test.eq.one) then
          ! do nothing
         else
          print *,"ice_test invalid ice_test=",ice_test
          print *,"cut_test=",cut_test
          stop
         endif
         if ((cut_test.ge.zero).and.(cut_test.le.one)) then
          ! do nothing
         else
          print *,"cut_test invalid" 
          print *,"cut_test ",cut_test
          stop
         endif
  
         if (dir.eq.0) then
          xface(D_DECL(i,j,k),icefacecut_index+1)=icefacecut
          xface(D_DECL(i,j,k),icemask_index+1)=icemask
         else if (dir.eq.1) then
          yface(D_DECL(i,j,k),icefacecut_index+1)=icefacecut
          yface(D_DECL(i,j,k),icemask_index+1)=icemask
         else if ((dir.eq.2).and.(SDIM.eq.3)) then
          zface(D_DECL(i,j,k),icefacecut_index+1)=icefacecut
          zface(D_DECL(i,j,k),icemask_index+1)=icemask
         else
          print *,"dir invalid init_icemask 3"
          stop
         endif

        else if ((local_mask_right.eq.0).and. &
                 (local_mask_left.eq.0)) then
         ! do nothing
        else
         print *,"local_mask_right or local_mask_left invalid"
         stop
        endif

       enddo ! k
       enddo ! j
       enddo ! i
      enddo ! dir=0..sdim-1

      return
      end subroutine FORT_INIT_ICEMASK


      ! tag = 1 -> donor cell
      ! tag = 2 -> receving cell
      ! tag = 0 -> none of above
      subroutine FORT_TAGEXPANSION(&
      latent_heat, &
      freezing_model, &
      distribute_from_target, &
      ngrow_expansion, &
      time, &
      vofbc, &
      expect_mdot_sign, &
      mdot_sum, &
      mdot_sum_comp, &
      im_source, &
      im_dest, &
      indexEXP, &
      level,finest_level, &
      nmat,nten, &
      tilelo,tilehi, &
      fablo,fabhi, &
      bfact, &
      xlo,dx,dt, &
      maskcov,DIMS(maskcov), &
      tag, &
      DIMS(tag), &
      tag_comp, &
      DIMS(tag_comp), &
      expan,DIMS(expan), &
      expan_comp,DIMS(expan_comp), &
      LS,DIMS(LS), &  ! newdistfab=(*localMF[LSNEW_MF])[mfi]
      recon,DIMS(recon))
      use probf90_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: ngrow_expansion
      REAL_T, intent(in) :: time
      REAL_T, intent(inout) :: mdot_sum
      REAL_T, intent(inout) :: mdot_sum_comp
      REAL_T, intent(in) :: expect_mdot_sign
      INTEGER_T, intent(in) :: im_source,im_dest
      INTEGER_T :: im_ice
      INTEGER_T, intent(in) :: indexEXP
      INTEGER_T, intent(in) :: level,finest_level
      INTEGER_T, intent(in) :: nmat,nten
      INTEGER_T :: nten_test
      REAL_T, intent(in) :: latent_heat(2*nten)
      INTEGER_T, intent(in) :: freezing_model(2*nten)
      INTEGER_T, intent(in) :: distribute_from_target(2*nten)
      INTEGER_T :: local_distribute_from_target(2*nten)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: dt
      INTEGER_T, intent(in) :: DIMDEC(maskcov)
      INTEGER_T, intent(in) :: DIMDEC(tag)
      INTEGER_T, intent(in) :: DIMDEC(tag_comp)
      INTEGER_T, intent(in) :: DIMDEC(expan)
      INTEGER_T, intent(in) :: DIMDEC(expan_comp)
      INTEGER_T, intent(in) :: DIMDEC(LS)
      INTEGER_T, intent(in) :: DIMDEC(recon)
      REAL_T, intent(in) :: maskcov(DIMV(maskcov))
      REAL_T, intent(out) :: tag(DIMV(tag))
      REAL_T, intent(out) :: tag_comp(DIMV(tag_comp))
      REAL_T, intent(in) :: expan(DIMV(expan),2*nten)
      REAL_T, intent(in) :: expan_comp(DIMV(expan_comp),2*nten)
      REAL_T, intent(in) :: LS(DIMV(LS),nmat*(1+SDIM))
      REAL_T, intent(in) :: recon(DIMV(recon),nmat*ngeom_recon)

      INTEGER_T local_freezing_model
      INTEGER_T vofbc(SDIM,2)
      INTEGER_T i,j,k
      REAL_T VFRAC(nmat)
      REAL_T VDOT
      REAL_T LSCELL(nmat)
      REAL_T ICEMASK
      REAL_T icefacecut
      INTEGER_T im,im_opp
      INTEGER_T ireverse
      INTEGER_T iten
      INTEGER_T im_primary
      INTEGER_T vofcomp
      INTEGER_T nhalf
      REAL_T xsten(-3:3,SDIM)
      REAL_T xsten_center(SDIM)
      INTEGER_T local_mask
      INTEGER_T dir
      INTEGER_T index_compare
      INTEGER_T complement_flag

      nhalf=3

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      !! Sanity checks

      if (time.lt.zero) then
       print *,"time invalid"
       stop
      endif
      if (ngrow_expansion.ne.2) then
       print *,"ngrow_expansion invalid"
       stop
      endif

      if ((im_source.ge.1).and.(im_source.le.nmat)) then
       ! do nothing
      else
       print *,"im_source invalid"
       stop
      endif
      if ((im_dest.ge.1).and.(im_dest.le.nmat)) then
       ! do nothing
      else
       print *,"im_dest invalid"
       stop
      endif
      if (im_dest.eq.im_source) then
       print *,"im_dest or im_source invalid"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in distribute_expansion"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid tagexpansion nten, nten_test ",nten,nten_test
       stop
      endif
      if ((indexEXP.lt.0).or.(indexEXP.ge.2*nten)) then
       print *,"indexEXP invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(maskcov),1,-1,1272)
      call checkbound(fablo,fabhi,DIMS(LS),1,-1,1272)
      call checkbound(fablo,fabhi,DIMS(recon),1,-1,1273)
      call checkbound(fablo,fabhi,DIMS(expan),ngrow_expansion,-1,1274)
      call checkbound(fablo,fabhi,DIMS(expan_comp),ngrow_expansion,-1,1274)
      call checkbound(fablo,fabhi,DIMS(tag),2*ngrow_expansion,-1,1275)
      call checkbound(fablo,fabhi,DIMS(tag_comp),2*ngrow_expansion,-1,1275)

      local_freezing_model=freezing_model(indexEXP+1)
      if ((local_freezing_model.lt.0).or. &
          (local_freezing_model.gt.7)) then
       print *,"local_freezing_model invalid 17"
       stop
      endif
      if ((distribute_from_target(indexEXP+1).lt.0).or. &
          (distribute_from_target(indexEXP+1).gt.1)) then
       print *,"distribute_from_target invalid"
       stop
      endif

      ! Iterate over the box
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       local_mask=NINT(maskcov(D_DECL(i,j,k)))

       if (local_mask.eq.1) then

        call gridsten_level(xsten,i,j,k,level,nhalf)
        do dir=1,SDIM
         xsten_center(dir)=xsten(0,dir)
        enddo

        !! water freezing to ice
        !! im_source -> index for water material
        !! im_dest   -> index for ice material
            
         ! check for being a donor cell:
         ! 1. non-zero expansion term
         ! 2. (F_ice < 0.5) or
         !    (icemask=0.0)
         ! 
         ! check for being a receiving cell:
         ! (F_ice > 0.5)&&(icemask>0.0)
         ! 
         ! weight for redistribution: inversely proportional to distance
         ! between donor and receiver.
              
        tag(D_DECL(i,j,k)) = zero
        tag_comp(D_DECL(i,j,k)) = zero

        do im=1,nmat
         vofcomp=(im-1)*ngeom_recon+1
         VFRAC(im)=recon(D_DECL(i,j,k),vofcomp)
         LSCELL(im)=LS(D_DECL(i,j,k),im)
        enddo

         ! first checks the rigid materials for a positive LS; if none
         ! exist, then the subroutine checks the fluid materials.
        call get_primary_material(LSCELL,nmat,im_primary)

        VDOT=expan(D_DECL(i,j,k),indexEXP+1)
        if (expect_mdot_sign.eq.one) then
         if (VDOT.lt.zero) then
          print *,"expecting VDOT>=0"
          print *,"sign+ VDOT invalid i,j,k,vdot ",i,j,k,VDOT
          stop
         else if (VDOT.ge.zero) then
          ! do nothing
         else
          print *,"VDOT bust"
          stop
         endif
        else if (expect_mdot_sign.eq.-one) then
         if (VDOT.gt.zero) then
          print *,"expecting VDOT<=0"
          print *,"sign- VDOT invalid i,j,k,vdot ",i,j,k,VDOT
          stop
         else if (VDOT.le.zero) then
          ! do nothing
         else
          print *,"VDOT bust"
          stop
         endif
        else
         print *,"expect_mdot_sign"
         stop
        endif

        do complement_flag=0,1
         do iten=1,2*nten
          local_distribute_from_target(iten)=distribute_from_target(iten)
          if (complement_flag.eq.0) then
           ! do nothing
          else if (complement_flag.eq.1) then
           local_distribute_from_target(iten)= &
              1-local_distribute_from_target(iten)
          else
           print *,"complement_flag invalid"
           stop
          endif
         enddo ! iten=1..2*nten

         im_ice=0
         if (is_ice(nmat,im_dest).eq.1) then
          im_ice=im_dest
         else if (is_ice(nmat,im_source).eq.1) then
          im_ice=im_source
         else if ((is_ice(nmat,im_dest).eq.0).and. &
                  (is_ice(nmat,im_source).eq.0)) then
          im_ice=0
         else
          print *,"is_ice invalid"
          stop
         endif

         if (im_ice.eq.0) then ! neither source nor dest are ice materials.
          ICEMASK=one
         else if ((im_ice.ge.1).and. &
                  (im_ice.le.nmat)) then

          ! in: FORT_TAGEXPANSION
          ! ICEMASK=0 => mask off this cell.
          ! ICEMASK=1 => do nothing
          ! get_icemask declared in PROB.F90
          call get_icemask( &
           xsten_center, &
           time, &
           dx,bfact, &
           ICEMASK, &
           icefacecut, &
           im,im_opp, &
           ireverse, &
           LSCELL, &
           latent_heat, &
           local_distribute_from_target, &
           nmat,nten)

          if (ireverse.eq.-1) then
           ICEMASK=one
          else if ((ireverse.eq.0).or.(ireverse.eq.1)) then
           call get_iten(im,im_opp,iten,nmat)
           index_compare=iten+ireverse*nten-1
           if ((index_compare.ge.0).and.(index_compare.lt.2*nten)) then
            if (index_compare.eq.indexEXP) then
             ! do nothing
            else
             ICEMASK=one
            endif
           else
            print *,"index_compare invalid"
            stop
           endif
          else
           print *,"ireverse invalid"
           stop
          endif
 
         else 
          print *,"im_ice invalid"
          stop
         endif

         if (complement_flag.eq.0) then
          mdot_sum=mdot_sum+VDOT
         else if (complement_flag.eq.1) then
          mdot_sum_comp=mdot_sum_comp+VDOT
         else
          print *,"complement_flag invalid"
          stop
         endif

         ! FSI_flag=0,7,3,6 ok
         if ((is_rigid(nmat,im_primary).eq.0).and. &
             (is_FSI_rigid(nmat,im_primary).eq.0)) then 

          if (VDOT.ne.zero) then ! nonzero source

           if (local_distribute_from_target(indexEXP+1).eq.0) then

            if ((VFRAC(im_dest).lt.half).or. &
                (ICEMASK.le.zero)) then
             if (complement_flag.eq.0) then
              tag(D_DECL(i,j,k)) = one ! donor cell
             else if (complement_flag.eq.1) then
              tag_comp(D_DECL(i,j,k)) = one ! donor cell
             else
              print *,"complement_flag invalid"
              stop
             endif
            else if ((VFRAC(im_dest).ge.half).and. &
                     (ICEMASK.gt.zero)) then
             ! do nothing - acceptor cell
            else
             print *,"VFRAC or ICEMASK bust"
             stop      
            endif

           else if (local_distribute_from_target(indexEXP+1).eq.1) then

            if ((VFRAC(im_source).lt.half).or. &
                (ICEMASK.le.zero)) then
             if (complement_flag.eq.0) then
              tag(D_DECL(i,j,k)) = one ! donor cell
             else if (complement_flag.eq.1) then
              tag_comp(D_DECL(i,j,k)) = one ! donor cell
             else
              print *,"complement_flag invalid"
              stop
             endif
            else if ((VFRAC(im_source).ge.half).and. &
                     (ICEMASK.gt.zero)) then
             ! do nothing - acceptor cell
            else
             print *,"VFRAC or ICEMASK bust"     
             stop
            endif

           else
            print *,"local_distribute_from_target(indexEXP+1) invalid"
            stop
           endif

          else if (VDOT.eq.zero) then
           ! do nothing
          else
           print *,"VDOT became corrupt"
           stop
          endif 

          if (local_distribute_from_target(indexEXP+1).eq.0) then

           if ((VFRAC(im_dest).ge.half).and. &
               (ICEMASK.gt.zero)) then
            if (complement_flag.eq.0) then
             tag(D_DECL(i,j,k)) = two ! receiver
            else if (complement_flag.eq.1) then
             tag_comp(D_DECL(i,j,k)) = two ! receiver
            else
             print *,"complement_flag invalid"
             stop
            endif
           else if ((VFRAC(im_dest).lt.half).or. &
                    (ICEMASK.le.zero)) then
            ! do nothing - donor cell
           else
            print *,"VFRAC or ICEMASK bust"
            stop      
           endif
 
          else if (local_distribute_from_target(indexEXP+1).eq.1) then

           if ((VFRAC(im_source).ge.half).and. &
               (ICEMASK.gt.zero)) then
            if (complement_flag.eq.0) then
             tag(D_DECL(i,j,k)) = two ! receiver
            else if (complement_flag.eq.1) then
             tag_comp(D_DECL(i,j,k)) = two ! receiver
            else
             print *,"complement_flag invalid"
             stop
            endif
           else if ((VFRAC(im_source).lt.half).or. &
                    (ICEMASK.le.zero)) then
            ! do nothing - donor cell
           else
            print *,"VFRAC or ICEMASK bust"     
            stop
           endif
 
          else
           print *,"local_distribute_from_target(indexEXP+1) invalid"
           stop
          endif

         !in the prescribed solid.
         else if ((is_rigid(nmat,im_primary).eq.1).or. &
                  (is_FSI_rigid(nmat,im_primary).eq.1)) then 
          ! do nothing (tag initialized to 0, neither donor nor receiver)
         else
          print *,"is_rigid(nmat,im_primary) invalid or"
          print *,"is_FSI_rigid(nmat,im_primary) invalid"
          stop
         endif 

        enddo ! complement_flag=0..1

       else if (local_mask.eq.0) then
        ! do nothing
       else
        print *,"local_mask invalid"
        stop
       endif

      enddo ! k
      enddo ! j
      enddo ! i

      return
      end subroutine FORT_TAGEXPANSION


      ! recon( nmat * ngeom_recon )
      ! ngeom_recon=2 * SDIM + 3
      ! volume fraction, reference centroid, order, slope, intercept
      ! nten =  (nmat * (nmat-1))/2  -> number of possible surface 
      !                                 tension coefficients
      ! if nmat=2, nten=1
      ! if nmat=3, nten=3    12 13 23
      ! if nmat=4, nten=6    12 13 14 23 24 34
      
      ! This will be called before FORT_INITJUMPTERM and after
      ! FORT_CONVERTMATERIAL
      ! tag = 1 -> donor cell
      ! tag = 2 -> receving cell
      ! tag = 0 -> non of above
      subroutine FORT_DISTRIBUTEEXPANSION(&
       ngrow_expansion, &
       im_source, &
       im_dest, &
       indexEXP, &
       level,finest_level, &
       nmat,nten, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx,dt, &
       maskcov,DIMS(maskcov),&
       LS,DIMS(LS),&
       tag, &
       DIMS(tag),&
       tag_comp, &
       DIMS(tag_comp),&
       expan, &
       DIMS(expan), &
       expan_comp, &
       DIMS(expan_comp) )
       use probf90_module
       use global_utility_module
       use geometry_intersect_module


       IMPLICIT NONE

       INTEGER_T, intent(in) :: im_source,im_dest,indexEXP,ngrow_expansion
       INTEGER_T, intent(in) :: level,finest_level
       INTEGER_T, intent(in) :: nmat,nten
       INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
       INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
       INTEGER_T :: growlo(3),growhi(3)
       INTEGER_T :: stenlo(3),stenhi(3)
       INTEGER_T :: stenlo2(3),stenhi2(3)
       INTEGER_T, intent(in) :: bfact
       REAL_T, intent(in) :: xlo(SDIM)
       REAL_T, intent(in) :: dx(SDIM)
       REAL_T, intent(in) :: dt
       INTEGER_T, intent(in) :: DIMDEC(maskcov)
       INTEGER_T, intent(in) :: DIMDEC(LS)
       INTEGER_T, intent(in) :: DIMDEC(tag)
       INTEGER_T, intent(in) :: DIMDEC(tag_comp)
       INTEGER_T, intent(in) :: DIMDEC(expan)
       INTEGER_T, intent(in) :: DIMDEC(expan_comp)
       REAL_T, intent(in) :: maskcov(DIMV(maskcov))
       REAL_T, intent(in) :: LS(DIMV(LS),nmat)
       REAL_T, intent(in) :: tag(DIMV(tag))
       REAL_T, intent(in) :: tag_comp(DIMV(tag_comp))
       REAL_T, intent(inout) :: expan(DIMV(expan),2*nten)
       REAL_T, intent(inout) :: expan_comp(DIMV(expan_comp),2*nten)

       INTEGER_T i,j,k,isweep
       REAL_T DLS
       REAL_T maxgrad,maxgrad2,curgrad
       INTEGER_T dir
       INTEGER_T i_n,j_n,k_n
       INTEGER_T i_nn,j_nn,k_nn
       REAL_T weight,total_weight,crit_weight
       REAL_T total_weight2,crit_weight2
       REAL_T TAGLOC,TAGSIDE
       REAL_T xsten_n(-1:1,SDIM)
       REAL_T xsten_nn(-1:1,SDIM)
       REAL_T crit_ratio,factor
       INTEGER_T is_inner,is_inner_main,nhalf
       REAL_T LS_receiver,LS_donor
       INTEGER_T local_mask

       nhalf=1

         ! if normal is close to inbetween two receiving cells, then include
         ! both receiving cells in stencil for the donor.
       crit_ratio=sqrt(4.0/5.0)

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

       !! Sanity checks

       if (ngrow_expansion.ne.2) then
        print *,"ngrow_expansion invalid"
        stop
       endif
       if ((im_source.ge.1).and.(im_source.le.nmat)) then
        ! do nothing
       else
        print *,"im_source invalid"
        stop
       endif
       if ((im_dest.ge.1).and.(im_dest.le.nmat)) then
        ! do nothing
       else
        print *,"im_dest invalid"
        stop
       endif
       if (im_dest.eq.im_source) then
        print *,"im_dest or im_source invalid"
        stop
       endif

       if ((level.lt.0).or.(level.gt.finest_level)) then
        print *,"level invalid in distribute_expansion"
        stop
       end if
       if ((indexEXP.lt.0).or.(indexEXP.ge.2*nten)) then
        print *,"indexEXP invalid"
        stop
       endif
       if (bfact.lt.1) then
        print *,"bfact too small"
        stop
       endif

       call checkbound(fablo,fabhi,DIMS(maskcov),1,-1,1272)
       call checkbound(fablo,fabhi,DIMS(LS),2*ngrow_expansion,-1,1272)
       call checkbound(fablo,fabhi,DIMS(expan),ngrow_expansion,-1,122)
       call checkbound(fablo,fabhi,DIMS(expan_comp),ngrow_expansion,-1,122)
       call checkbound(fablo,fabhi,DIMS(tag),2*ngrow_expansion,-1,122)
       call checkbound(fablo,fabhi,DIMS(tag_comp),2*ngrow_expansion,-1,122)

       ! Iterate over the box
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)

        local_mask=NINT(maskcov(D_DECL(i,j,k)))

        if (local_mask.eq.1) then

         ! if a receiving cell
         TAGLOC=tag(D_DECL(i,j,k))
         if(TAGLOC.eq.two) then
          call stencilbox(i,j,k,fablo,fabhi,stenlo,stenhi,ngrow_expansion)

          do i_n=stenlo(1),stenhi(1)
          do j_n=stenlo(2),stenhi(2)
          do k_n=stenlo(3),stenhi(3)

           is_inner_main=1
           if ((abs(i_n-i).gt.1).or. &
               (abs(j_n-j).gt.1).or. &
               (abs(k_n-k).gt.1)) then
            is_inner_main=0
           endif

           ! if there is donor neighbor cell
           TAGSIDE=tag(D_DECL(i_n,j_n,k_n))
           if(TAGSIDE.eq.one) then

            call gridsten_level(xsten_n,i_n,j_n,k_n,level,nhalf)
            call stencilbox(i_n,j_n,k_n,fablo,fabhi,stenlo2,stenhi2, &
              ngrow_expansion)

              ! first: find the maximum of |d phi/d n_grid|
              ! second: use the maximum to find the weights.
            crit_weight=zero
            total_weight=zero
            maxgrad=-one

            crit_weight2=zero
            total_weight2=zero
            maxgrad2=-one

            do isweep=0,1

             do i_nn=stenlo2(1),stenhi2(1)
             do j_nn=stenlo2(2),stenhi2(2)
             do k_nn=stenlo2(3),stenhi2(3)

              is_inner=1
              if ((abs(i_nn-i_n).gt.1).or. &
                  (abs(j_nn-j_n).gt.1).or. &
                  (abs(k_nn-k_n).gt.1)) then
               is_inner=0
              endif

              if (tag(D_DECL(i_nn,j_nn,k_nn)).eq.two) then ! receiver
               call gridsten_level(xsten_nn,i_nn,j_nn,k_nn,level,nhalf)
               weight=zero
               do dir=1,SDIM
                weight=weight+(xsten_n(0,dir)-xsten_nn(0,dir))**2
               enddo
               if (weight.le.zero) then
                print *,"weight invalid"
                stop
               endif
               LS_receiver=LS(D_DECL(i_nn,j_nn,k_nn),im_dest)
               LS_donor=LS(D_DECL(i_n,j_n,k_n),im_dest)
               DLS=(LS_receiver-LS_donor)**2
               curgrad=sqrt(DLS/weight)
               if (isweep.eq.0) then

                if (is_inner.eq.1) then
                 if ((curgrad.gt.maxgrad).or.(maxgrad.eq.-one)) then
                  maxgrad=curgrad
                 endif
                endif 

                if ((curgrad.gt.maxgrad2).or.(maxgrad2.eq.-one)) then
                 maxgrad2=curgrad
                endif
               
               else if (isweep.eq.1) then

                if (is_inner.eq.1) then

                 if (maxgrad.lt.zero) then
                  weight=zero
                  if (is_inner_main.eq.1) then
                   print *,"maxgrad invalid"
                   stop
                  endif
                 else if (maxgrad.eq.zero) then
                  weight=one
                 else if (curgrad/maxgrad.lt.crit_ratio) then
                  weight=zero
                 else
                  weight=(curgrad/maxgrad)**2
                 endif
                 total_weight=total_weight+weight
                 if ((i.eq.i_nn).and.(j.eq.j_nn).and.(k.eq.k_nn)) then
                  crit_weight=weight
                  if (is_inner_main.eq.0) then
                   print *,"is_inner_main invalid"
                   stop
                  endif
                 endif

                endif  ! is_inner=1

                if ((maxgrad2.lt.maxgrad).or.(maxgrad2.lt.zero)) then
                 print *,"maxgrad2 invalid"
                 stop
                else if (maxgrad2.eq.zero) then
                 weight=one
                else if (curgrad/maxgrad2.lt.crit_ratio) then
                 weight=zero
                else
                 weight=(curgrad/maxgrad2)**2
                endif
                total_weight2=total_weight2+weight
                if ((i.eq.i_nn).and.(j.eq.j_nn).and.(k.eq.k_nn)) then
                 crit_weight2=weight
                endif

               else
                print *,"isweep invalid"
                stop
               endif
              end if  ! receiver
             end do ! k_nn
             end do ! j_nn
             end do ! i_nn

             if (isweep.eq.0) then
              if (maxgrad.lt.zero) then
               if (is_inner_main.eq.1) then
                print *,"maxgrad: a donor cell has no receiving cell"
                stop
               endif
              end if
              if (maxgrad2.lt.zero) then
               print *,"maxgrad2: a donor cell has no receiving cell"
               stop
              end if
             else if (isweep.eq.1) then
              if (crit_weight.lt.zero) then
               print *,"crit_weight invalid"
               stop
              endif
              if (crit_weight2.lt.zero) then
               print *,"crit_weight2 invalid"
               stop
              endif
              if (total_weight.gt.zero) then
               factor=crit_weight/total_weight
              else if ((is_inner_main.eq.0).and. &
                       (total_weight2.gt.zero)) then
               factor=crit_weight2/total_weight2
              else
               print *,"a donor cell has no receiving cell"
               stop
              end if

               ! transfer from donor to receiving cell
               ! with weight crit_weight/total_weight
              expan(D_DECL(i,j,k),indexEXP+1) = &
               expan(D_DECL(i,j,k),indexEXP+1) + &
               expan(D_DECL(i_n,j_n,k_n),indexEXP+1)*factor
             else
              print *,"isweep invalid"
              stop
             endif
            enddo ! isweep=0..1

           else if ((TAGSIDE.eq.two).or.(TAGSIDE.eq.zero)) then
             ! do nothing
           else
            print *,"TAGSIDE invalid"
            stop
           end if ! donor cell
          end do ! k_n
          end do ! j_n
          end do ! i_n

         else if ((TAGLOC.eq.one).or.(TAGLOC.eq.zero)) then
          ! do nothing
         else
          print *,"TAGLOC invalid"
          stop
         end if ! receiving cell

          ! ---------------- DISTRIBUTE FOR COMPLEMENT ----------------

         ! if a receiving cell
         TAGLOC=tag_comp(D_DECL(i,j,k))
         if(TAGLOC.eq.two) then
          call stencilbox(i,j,k,fablo,fabhi,stenlo,stenhi,ngrow_expansion)

          do i_n=stenlo(1),stenhi(1)
          do j_n=stenlo(2),stenhi(2)
          do k_n=stenlo(3),stenhi(3)

           is_inner_main=1
           if ((abs(i_n-i).gt.1).or. &
               (abs(j_n-j).gt.1).or. &
               (abs(k_n-k).gt.1)) then
            is_inner_main=0
           endif

           ! if there is donor neighbor cell
           TAGSIDE=tag_comp(D_DECL(i_n,j_n,k_n))
           if(TAGSIDE.eq.one) then

            call gridsten_level(xsten_n,i_n,j_n,k_n,level,nhalf)
            call stencilbox(i_n,j_n,k_n,fablo,fabhi,stenlo2,stenhi2, &
              ngrow_expansion)

              ! first: find the maximum of |d phi/d n_grid|
              ! second: use the maximum to find the weights.
            crit_weight=zero
            total_weight=zero
            maxgrad=-one

            crit_weight2=zero
            total_weight2=zero
            maxgrad2=-one

            do isweep=0,1

             do i_nn=stenlo2(1),stenhi2(1)
             do j_nn=stenlo2(2),stenhi2(2)
             do k_nn=stenlo2(3),stenhi2(3)

              is_inner=1
              if ((abs(i_nn-i_n).gt.1).or. &
                  (abs(j_nn-j_n).gt.1).or. &
                  (abs(k_nn-k_n).gt.1)) then
               is_inner=0
              endif

              if (tag_comp(D_DECL(i_nn,j_nn,k_nn)).eq.two) then ! receiver
               call gridsten_level(xsten_nn,i_nn,j_nn,k_nn,level,nhalf)
               weight=zero
               do dir=1,SDIM
                weight=weight+(xsten_n(0,dir)-xsten_nn(0,dir))**2
               enddo
               if (weight.le.zero) then
                print *,"weight invalid"
                stop
               endif
               LS_receiver=LS(D_DECL(i_nn,j_nn,k_nn),im_dest)
               LS_donor=LS(D_DECL(i_n,j_n,k_n),im_dest)
               DLS=(LS_receiver-LS_donor)**2
               curgrad=sqrt(DLS/weight)
               if (isweep.eq.0) then

                if (is_inner.eq.1) then
                 if ((curgrad.gt.maxgrad).or.(maxgrad.eq.-one)) then
                  maxgrad=curgrad
                 endif
                endif 

                if ((curgrad.gt.maxgrad2).or.(maxgrad2.eq.-one)) then
                 maxgrad2=curgrad
                endif
               
               else if (isweep.eq.1) then

                if (is_inner.eq.1) then

                 if (maxgrad.lt.zero) then
                  weight=zero
                  if (is_inner_main.eq.1) then
                   print *,"maxgrad invalid"
                   stop
                  endif
                 else if (maxgrad.eq.zero) then
                  weight=one
                 else if (curgrad/maxgrad.lt.crit_ratio) then
                  weight=zero
                 else
                  weight=(curgrad/maxgrad)**2
                 endif
                 total_weight=total_weight+weight
                 if ((i.eq.i_nn).and.(j.eq.j_nn).and.(k.eq.k_nn)) then
                  crit_weight=weight
                  if (is_inner_main.eq.0) then
                   print *,"is_inner_main invalid"
                   stop
                  endif
                 endif

                endif  ! is_inner=1

                if ((maxgrad2.lt.maxgrad).or.(maxgrad2.lt.zero)) then
                 print *,"maxgrad2 invalid"
                 stop
                else if (maxgrad2.eq.zero) then
                 weight=one
                else if (curgrad/maxgrad2.lt.crit_ratio) then
                 weight=zero
                else
                 weight=(curgrad/maxgrad2)**2
                endif
                total_weight2=total_weight2+weight
                if ((i.eq.i_nn).and.(j.eq.j_nn).and.(k.eq.k_nn)) then
                 crit_weight2=weight
                endif

               else
                print *,"isweep invalid"
                stop
               endif
              end if  ! receiver
             end do ! k_nn
             end do ! j_nn
             end do ! i_nn

             if (isweep.eq.0) then
              if (maxgrad.lt.zero) then
               if (is_inner_main.eq.1) then
                print *,"maxgrad: a donor cell has no receiving cell"
                stop
               endif
              end if
              if (maxgrad2.lt.zero) then
               print *,"maxgrad2: a donor cell has no receiving cell"
               stop
              end if
             else if (isweep.eq.1) then
              if (crit_weight.lt.zero) then
               print *,"crit_weight invalid"
               stop
              endif
              if (crit_weight2.lt.zero) then
               print *,"crit_weight2 invalid"
               stop
              endif
              if (total_weight.gt.zero) then
               factor=crit_weight/total_weight
              else if ((is_inner_main.eq.0).and. &
                       (total_weight2.gt.zero)) then
               factor=crit_weight2/total_weight2
              else
               print *,"a donor cell has no receiving cell"
               stop
              end if

               ! transfer from donor to receiving cell
               ! with weight crit_weight/total_weight
              expan_comp(D_DECL(i,j,k),indexEXP+1) = &
               expan_comp(D_DECL(i,j,k),indexEXP+1) + &
               expan_comp(D_DECL(i_n,j_n,k_n),indexEXP+1)*factor
             else
              print *,"isweep invalid"
              stop
             endif
            enddo ! isweep=0..1

           else if ((TAGSIDE.eq.two).or.(TAGSIDE.eq.zero)) then
             ! do nothing
           else
            print *,"TAGSIDE invalid"
            stop
           end if ! donor cell
          end do ! k_n
          end do ! j_n
          end do ! i_n

         else if ((TAGLOC.eq.one).or.(TAGLOC.eq.zero)) then
          ! do nothing
         else
          print *,"TAGLOC invalid"
          stop
         end if ! receiving cell


        else if (local_mask.eq.0) then
         ! do nothing
        else
         print *,"local_mask invalid"
         stop
        endif

       end do ! k
       end do ! j
       end do ! i
      end subroutine FORT_DISTRIBUTEEXPANSION


      ! tag = 1 -> donor cell
      ! tag = 2 -> receving cell
      ! tag = 0 -> non of above
      subroutine FORT_CLEAREXPANSION(&
       ngrow_expansion, &
       mdot_sum, &
       mdot_lost, &
       mdot_sum_comp, &
       mdot_lost_comp, &
       im_source, &
       im_dest, &
       indexEXP, &
       level,finest_level, &
       nmat,nten, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx,dt, &
       maskcov,DIMS(maskcov),&
       tag, &
       DIMS(tag),&
       tag_comp, &
       DIMS(tag_comp),&
       expan, &
       DIMS(expan), &
       expan_comp, &
       DIMS(expan_comp) )
       use probf90_module
       use global_utility_module
       use geometry_intersect_module

       IMPLICIT NONE

       INTEGER_T, intent(in) :: ngrow_expansion
       REAL_T, intent(inout) :: mdot_sum,mdot_lost
       REAL_T, intent(inout) :: mdot_sum_comp,mdot_lost_comp
       INTEGER_T, intent(in) :: im_source,im_dest,indexEXP
       INTEGER_T, intent(in) :: level,finest_level
       INTEGER_T, intent(in) :: nmat,nten
       INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
       INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
       INTEGER_T :: growlo(3),growhi(3)
       INTEGER_T :: stenlo(3),stenhi(3)
       INTEGER_T, intent(in) :: bfact
       REAL_T, intent(in) :: xlo(SDIM)
       REAL_T, intent(in) :: dx(SDIM)
       REAL_T, intent(in) :: dt
       INTEGER_T, intent(in) :: DIMDEC(maskcov)
       INTEGER_T, intent(in) :: DIMDEC(tag)
       INTEGER_T, intent(in) :: DIMDEC(tag_comp)
       INTEGER_T, intent(in) :: DIMDEC(expan)
       INTEGER_T, intent(in) :: DIMDEC(expan_comp)
       REAL_T, intent(in) :: maskcov(DIMV(maskcov))
       REAL_T, intent(in) :: tag(DIMV(tag))
       REAL_T, intent(in) :: tag_comp(DIMV(tag_comp))
       REAL_T, intent(inout) :: expan(DIMV(expan),2*nten)
       REAL_T, intent(inout) :: expan_comp(DIMV(expan_comp),2*nten)

       INTEGER_T i,j,k
       INTEGER_T i_n,j_n,k_n,receive_flag,nhalf
       INTEGER_T TAGLOC,TAGSIDE
       REAL_T xsten(-1:1,SDIM)

       INTEGER_T local_mask

       nhalf=1

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

       !! Sanity checks
       if (ngrow_expansion.ne.2) then
        print *,"ngrow_expansion invalid"
        stop
       endif
       if ((im_source.ge.1).and.(im_source.le.nmat)) then
        ! do nothing
       else
        print *,"im_source invalid"
        stop
       endif
       if ((im_dest.ge.1).and.(im_dest.le.nmat)) then
        ! do nothing
       else
        print *,"im_dest invalid"
        stop
       endif
       if (im_dest.eq.im_source) then
        print *,"im_dest or im_source invalid"
        stop
       endif

       if ((level.lt.0).or.(level.gt.finest_level)) then
        print *,"level invalid in clear_expansion"
        stop
       end if
       if ((indexEXP.lt.0).or.(indexEXP.ge.2*nten)) then
        print *,"indexEXP invalid"
        stop
       endif
       if (bfact.lt.1) then
        print *,"bfact too small"
        stop
       endif

       call checkbound(fablo,fabhi,DIMS(maskcov),1,-1,122)
       call checkbound(fablo,fabhi,DIMS(expan),ngrow_expansion,-1,122)
       call checkbound(fablo,fabhi,DIMS(expan_comp),ngrow_expansion,-1,122)
       call checkbound(fablo,fabhi,DIMS(tag),2*ngrow_expansion,-1,122)
       call checkbound(fablo,fabhi,DIMS(tag_comp),2*ngrow_expansion,-1,122)

       ! Iterate over the box
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)

        local_mask=NINT(maskcov(D_DECL(i,j,k)))

        if (local_mask.eq.1) then

         ! if a donor cell
         TAGLOC=NINT(tag(D_DECL(i,j,k)))
         if(TAGLOC.eq.1) then
          if (1.eq.1) then  ! SANITY CHECK
           call stencilbox(i,j,k,fablo,fabhi,stenlo,stenhi,ngrow_expansion)
           receive_flag=0
           do i_n=stenlo(1),stenhi(1)
           do j_n=stenlo(2),stenhi(2)
           do k_n=stenlo(3),stenhi(3)
            TAGSIDE=NINT(tag(D_DECL(i_n,j_n,k_n)))
            if (TAGSIDE.eq.2) then
             receive_flag=1
            endif
           enddo
           enddo
           enddo
           if (receive_flag.eq.0) then
            mdot_lost=mdot_lost+expan(D_DECL(i,j,k),indexEXP+1)
            call gridsten_level(xsten,i,j,k,level,nhalf)
            if (1.eq.0) then
             print *,"donor has no receiver"
             print *,"i,j,k ",i,j,k
             print *,"x,y,z ",xsten(0,1),xsten(0,2),xsten(0,SDIM)
             print *,"LOST: ",expan(D_DECL(i,j,k),indexEXP+1)
            endif
           endif
          endif
          expan(D_DECL(i,j,k),indexEXP+1)=zero
         else if ((TAGLOC.eq.2).or.(TAGLOC.eq.0)) then
          ! do nothing
         else
          print *,"TAGLOC invalid"
          stop
         endif
         mdot_sum=mdot_sum+expan(D_DECL(i,j,k),indexEXP+1)

          ! --------------- COMPLEMENT -----------------

          ! if a donor cell
         TAGLOC=NINT(tag_comp(D_DECL(i,j,k)))
         if(TAGLOC.eq.1) then
          if (1.eq.1) then  ! SANITY CHECK
           call stencilbox(i,j,k,fablo,fabhi,stenlo,stenhi,ngrow_expansion)
           receive_flag=0
           do i_n=stenlo(1),stenhi(1)
           do j_n=stenlo(2),stenhi(2)
           do k_n=stenlo(3),stenhi(3)
            TAGSIDE=NINT(tag_comp(D_DECL(i_n,j_n,k_n)))
            if (TAGSIDE.eq.2) then
             receive_flag=1
            endif
           enddo
           enddo
           enddo
           if (receive_flag.eq.0) then
            mdot_lost_comp=mdot_lost_comp+ &
                  expan_comp(D_DECL(i,j,k),indexEXP+1)
            call gridsten_level(xsten,i,j,k,level,nhalf)
            if (1.eq.0) then
             print *,"donor has no receiver (complement) "
             print *,"i,j,k ",i,j,k
             print *,"x,y,z ",xsten(0,1),xsten(0,2),xsten(0,SDIM)
             print *,"LOST: ",expan_comp(D_DECL(i,j,k),indexEXP+1)
            endif
           endif
          endif
          expan_comp(D_DECL(i,j,k),indexEXP+1)=zero
         else if ((TAGLOC.eq.2).or.(TAGLOC.eq.0)) then
          ! do nothing
         else
          print *,"TAGLOC invalid"
          stop
         endif
         mdot_sum_comp=mdot_sum_comp+ &
            expan_comp(D_DECL(i,j,k),indexEXP+1)

        else if (local_mask.eq.0) then
         ! do nothing
        else
         print *,"local_mask invalid"
         stop
        endif

       enddo
       enddo
       enddo ! i,j,k

      end subroutine FORT_CLEAREXPANSION


      subroutine FORT_SOD_SANITY( &
        id,nc,lo,hi, &
        snew,DIMS(snew))
      use probf90_module
      use CISL_SANITY_MODULE

      IMPLICIT NONE
      INTEGER_T id,nc
      INTEGER_T lo(SDIM),hi(SDIM)
      INTEGER_T DIMDEC(snew)
      REAL_T snew(DIMV(snew),nc)
      INTEGER_T i,j,k,nmat
      REAL_T, dimension(:,:), allocatable :: comparestate

      nmat=num_materials
      if (nc.ne.(SDIM+1)+ &
          nmat*(num_state_material+ngeom_raw)+1) then
       print *,"nc invalid in sod sanity"
       stop
      endif

      if (DO_SANITY_CHECK.eq.1) then
       print *,"SANITY CHECK FROM C++ id=",id
       allocate(comparestate(lo(1)-1:hi(1)+1,5))
       j=1
       k=0
       do i=lo(1),hi(1)
        comparestate(i,1)=snew(D_DECL(i,j,k),1)
        comparestate(i,2)=snew(D_DECL(i,j,k),SDIM+2)
        comparestate(i,3)=snew(D_DECL(i,j,k),SDIM+3)
       enddo 
       call compare_sanity(comparestate,1,3,8)
       deallocate(comparestate)
       print *,"AFTER SANITY CHECK FROM C++ id=",id
      else
       print *,"sod sanity check should not be called"
       stop
      endif

      return
      end

      subroutine FORT_AGGRESSIVE( &
       datatype, &
       warning_cutoff, &
       tilelo,tilehi, &
       fablo,fabhi, &
       growlo,growhi, &
       bfact, &
       dx, &
       scomp, &
       ncomp, &
       ndefined, &
       ngrow,dir,id, &
       verbose, &
       force_check, &
       gridno, &
       ngrid,level,finest_level, &
       mf,DIMS(mf))
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T datatype
      REAL_T warning_cutoff
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(SDIM),growhi(SDIM)
      INTEGER_T bfact
      REAL_T dx(SDIM)
      INTEGER_T scomp,ncomp,ndefined
      INTEGER_T ngrow,dir,id
      INTEGER_T verbose
      INTEGER_T force_check
      INTEGER_T gridno,ngrid,level,finest_level
      INTEGER_T DIMDEC(mf)
      REAL_T mf(DIMV(mf),ndefined)

      if (bfact.lt.1) then
       print *,"bfact invalid69"
       stop
      endif

      if (((verbose.eq.0).or.(verbose.eq.1)).and.(force_check.eq.0)) then
       ! do nothing
      else if ((verbose.eq.2).or.(force_check.eq.1)) then
       call aggressive_worker( &
        datatype, &
        warning_cutoff, &
        tilelo,tilehi, &
        fablo,fabhi, &
        growlo,growhi, &
        bfact, &
        dx, &
        scomp, &
        ncomp, &
        ndefined, &
        ngrow,dir,id, &
        verbose, &
        force_check, &
        gridno,ngrid,level,finest_level, &
        mf,DIMS(mf))
      else
       print *,"verbose invalid"
       stop
      endif

      return
      end subroutine FORT_AGGRESSIVE

         ! "coarray fortran"  (MPI functionality built in)
         ! masknbr=1.0 in the interior
         !        =1.0 fine-fine ghost cells
         !        =0.0 coarse-fine ghost cells and outside the domain.
         ! mask=tag if not covered by level+1 or outside the domain.
      subroutine FORT_VFRAC_SPLIT( &
       nprocessed, &
       tid, &
       added_weight, &
       density_floor, &
       density_ceiling, &
       solidheat_flag, &
       temperature_primitive_variable, &
       dencomp,mofcomp,errcomp, & 
       latent_heat, &
       freezing_model, &
       distribute_from_target, &
       nten, &
       face_flag, &
       override_density, &
       constant_density_all_time, &
       velbc, &
       EILE_flag, &
       dir_counter, &
       normdir, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       bfact_f, &
       dt, &
       time, &
       passive_veltime, &
       LS,DIMS(LS), &  ! original data
       den, &
       DIMS(den), &
       mom_den, &
       DIMS(mom_den), &
       tensor,DIMS(tensor), &
       velfab,DIMS(velfab), &
       PLICSLP,DIMS(PLICSLP), &  ! slope data
       snew,DIMS(snew), &  ! this is the result
       tennew,DIMS(tennew), & 
       LSnew,DIMS(LSnew), &
       ucell,DIMS(ucell), &  ! other vars
       vofls0,DIMS(vofls0), &  
       mask,DIMS(mask), & !mask=1 if not covered by level+1 or outside domain
       masknbr,DIMS(masknbr), &
       unode,DIMS(unode), & ! vel*dt
       xlo,dx, &
       conserve,DIMS(conserve), & ! local variables
       xvel,DIMS(xvel), & ! 1..num_MAC_vectors
       yvel,DIMS(yvel), &
       zvel,DIMS(zvel), &
       xvelslp,DIMS(xvelslp), & ! xvelslope,xcen
       yvelslp,DIMS(yvelslp), &
       zvelslp,DIMS(zvelslp), &
       momslope,DIMS(momslope), &
       xmomside,DIMS(xmomside), & ! 1..2*num_MAC_vectors
       ymomside,DIMS(ymomside), &
       zmomside,DIMS(zmomside), &
       xmassside,DIMS(xmassside), & ! 1..2*num_MAC_vectors
       ymassside,DIMS(ymassside), &
       zmassside,DIMS(zmassside), &
       ngrow, &
       ngrow_mac_old, &
       nc_conserve, &
       iden_base, &
       nmat, &
       map_forward, &
       recon_ncomp, &
       den_recon_ncomp, &
       ncomp_state, &
       ntensor, &
       nc_bucket, &
       nrefine_vof, &
       num_MAC_vectors, & !=1 or 2 (VFRAC_SPLIT)
       NUM_CELL_ELASTIC, &
       verbose, &
       gridno,ngrid, &
       level, &
       finest_level, &
       dombc, &
       domlo,domhi)
      use probf90_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      use godunov_module
      use CISL_SANITY_MODULE


      IMPLICIT NONE

      INTEGER_T, intent(in) :: num_MAC_vectors !=1 or 2
      INTEGER_T, intent(in) :: NUM_CELL_ELASTIC
      INTEGER_T, intent(inout) :: nprocessed
      INTEGER_T, intent(in) :: tid

      INTEGER_T, intent(in) :: iden_base
      INTEGER_T, intent(in) :: nc_conserve
      INTEGER_T, intent(in) :: nrefine_vof
      INTEGER_T, intent(in) :: ngrow,ngrow_mac_old
      INTEGER_T, intent(in) :: solidheat_flag
      INTEGER_T, intent(in) :: nten
      REAL_T, intent(in) :: latent_heat(2*nten)
      INTEGER_T, intent(in) :: freezing_model(2*nten)
      INTEGER_T, intent(in) :: distribute_from_target(2*nten)

      INTEGER_T, intent(in) :: dencomp,mofcomp,errcomp
      INTEGER_T, intent(in) :: face_flag
      INTEGER_T, intent(in) :: domlo(SDIM),domhi(SDIM)
      INTEGER_T, intent(in) :: dombc(SDIM,2)
      INTEGER_T, intent(in) :: EILE_flag
      REAL_T, intent(in) :: passive_veltime
      INTEGER_T, intent(in) :: dir_counter
      INTEGER_T, intent(in) :: normdir
      INTEGER_T, intent(in) :: verbose
      INTEGER_T :: force_check
      INTEGER_T, intent(in) :: gridno,ngrid
      INTEGER_T, intent(in) :: level,finest_level
      INTEGER_T, intent(in) :: recon_ncomp
      INTEGER_T, intent(in) :: den_recon_ncomp
      INTEGER_T, intent(in) :: ncomp_state
      INTEGER_T, intent(in) :: ntensor
      INTEGER_T, intent(in) :: nc_bucket
      INTEGER_T :: nc_bucket_test
      INTEGER_T, intent(in) :: map_forward
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: bfact_f
      INTEGER_T, intent(in) :: override_density(nmat)
      INTEGER_T, intent(in) :: constant_density_all_time(nmat)
      REAL_T, intent(in) :: dt,time
       ! original data
      INTEGER_T, intent(in) :: DIMDEC(LS)
      INTEGER_T, intent(in) :: DIMDEC(den)
      INTEGER_T, intent(in) :: DIMDEC(mom_den)
      INTEGER_T, intent(in) :: DIMDEC(tensor)
      INTEGER_T, intent(in) :: DIMDEC(velfab)
       ! slope data
      INTEGER_T, intent(in) :: DIMDEC(PLICSLP)
       ! new data
      INTEGER_T, intent(in) :: DIMDEC(snew)
      INTEGER_T, intent(in) :: DIMDEC(tennew)
      INTEGER_T, intent(in) :: DIMDEC(LSnew)
       ! other vars
      INTEGER_T, intent(in) :: DIMDEC(ucell)
      INTEGER_T, intent(in) :: DIMDEC(vofls0)
      INTEGER_T, intent(in) :: DIMDEC(mask)
      INTEGER_T, intent(in) :: DIMDEC(masknbr)
      INTEGER_T, intent(in) :: DIMDEC(unode)
       ! local variables
      INTEGER_T, intent(in) :: DIMDEC(conserve)
      INTEGER_T, intent(in) :: DIMDEC(xvel)
      INTEGER_T, intent(in) :: DIMDEC(yvel)
      INTEGER_T, intent(in) :: DIMDEC(zvel)
      INTEGER_T, intent(in) :: DIMDEC(xvelslp)
      INTEGER_T, intent(in) :: DIMDEC(yvelslp)
      INTEGER_T, intent(in) :: DIMDEC(zvelslp)
      INTEGER_T, intent(in) :: DIMDEC(momslope)

      INTEGER_T, intent(in) :: DIMDEC(xmomside) 
      INTEGER_T, intent(in) :: DIMDEC(ymomside) 
      INTEGER_T, intent(in) :: DIMDEC(zmomside) 
      INTEGER_T, intent(in) :: DIMDEC(xmassside) 
      INTEGER_T, intent(in) :: DIMDEC(ymassside) 
      INTEGER_T, intent(in) :: DIMDEC(zmassside) 

       ! FABS
       ! original data
      REAL_T, intent(in) :: LS(DIMV(LS),nmat)
      REAL_T, intent(in) :: den(DIMV(den),den_recon_ncomp)
      REAL_T, intent(in) :: mom_den(DIMV(mom_den),nmat)
      REAL_T, intent(in) :: tensor(DIMV(tensor),ntensor)
      REAL_T, intent(in) :: velfab(DIMV(velfab),SDIM+1)
       ! slope data
      REAL_T, intent(in) :: PLICSLP(DIMV(PLICSLP),recon_ncomp)
       ! new data
      REAL_T, intent(inout) :: snew(DIMV(snew),ncomp_state)
      REAL_T, intent(inout) :: tennew(DIMV(tennew),ntensor)
      REAL_T, intent(inout) :: LSnew(DIMV(LSnew),nmat)
       ! other vars
       ! displacement
      REAL_T, intent(in) :: ucell(DIMV(ucell),SDIM)
      REAL_T, intent(in) :: vofls0(DIMV(vofls0),2*nmat)
      REAL_T, intent(in) :: mask(DIMV(mask))
      ! =1 int. =1 fine-fine in domain =0 o.t.
      REAL_T, intent(in) :: masknbr(DIMV(masknbr)) 
      REAL_T, intent(in) :: unode(DIMV(unode))
       ! local variables
      REAL_T, intent(in) :: conserve(DIMV(conserve),nc_conserve)
      REAL_T, intent(in) :: xvel(DIMV(xvel), &
         num_MAC_vectors) 
      REAL_T, intent(in) :: yvel(DIMV(yvel), &
         num_MAC_vectors)  
      REAL_T, intent(in) :: zvel(DIMV(zvel), &
         num_MAC_vectors) 
      REAL_T, intent(in) :: xvelslp(DIMV(xvelslp),1+nmat)  ! xvelslope,xcen
      REAL_T, intent(in) :: yvelslp(DIMV(yvelslp),1+nmat)  
      REAL_T, intent(in) :: zvelslp(DIMV(zvelslp),1+nmat)  
      REAL_T, intent(in) :: momslope(DIMV(momslope),nc_conserve)

      REAL_T, intent(inout) :: xmomside(DIMV(xmomside), &
        2*num_MAC_vectors)
      REAL_T, intent(inout) :: ymomside(DIMV(ymomside), &
        2*num_MAC_vectors)
      REAL_T, intent(inout) :: zmomside(DIMV(zmomside), &
        2*num_MAC_vectors)

      REAL_T, intent(inout) :: xmassside(DIMV(xmassside), &
        2*num_MAC_vectors)
      REAL_T, intent(inout) :: ymassside(DIMV(ymassside), &
        2*num_MAC_vectors)
      REAL_T, intent(inout) :: zmassside(DIMV(zmassside), &
        2*num_MAC_vectors)
    
      INTEGER_T, intent(in) :: temperature_primitive_variable(nmat) 
      REAL_T, intent(in) :: density_floor(nmat)
      REAL_T, intent(in) :: density_ceiling(nmat)
      REAL_T, intent(in) :: added_weight(nmat)

      INTEGER_T, intent(in) :: velbc(SDIM,2)

      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)

      INTEGER_T ibucket
      INTEGER_T ibucket_map
      REAL_T xsten_crse(-1:1,SDIM)
      INTEGER_T dir2
      INTEGER_T iside
      INTEGER_T iside_part
      INTEGER_T isidedonate
      INTEGER_T iside_low,iside_high
      INTEGER_T vofcomp
      INTEGER_T im
      REAL_T mom2(SDIM)
      REAL_T xsten_MAC(-1:1,SDIM)
      REAL_T xsten_accept(-1:1,SDIM)
      REAL_T xsten_donate(-1:1,SDIM)
      REAL_T xsten_target(-1:1,SDIM)
      REAL_T xsten_depart(-1:1,SDIM)
      REAL_T usten_accept(-1:1)
      REAL_T usten_donate(-1:1)

      INTEGER_T null_velocity_flag

      REAL_T xdepartsize,xtargetsize,xloint,xhiint
      REAL_T volint
      REAL_T coeff(2)
      INTEGER_T tessellate
      INTEGER_T nmax
      INTEGER_T ii,jj,kk
     
      INTEGER_T veldir
      INTEGER_T veldir_comp

      REAL_T totalmass_depart
     
      INTEGER_T istate,ispecies,igeom
    
      REAL_T KE,vel1D,local_internal
      INTEGER_T no_material_flag

      INTEGER_T dencomp_data,statecomp_data,tempcomp_data,speccomp_data
      INTEGER_T imap

      REAL_T volcell_recon
      REAL_T cencell_recon(SDIM)
      REAL_T volcell_accept
      REAL_T cencell_accept(SDIM)
      REAL_T volcell_donate
      REAL_T cencell_donate(SDIM)
      INTEGER_T iii,jjj,kkk

      REAL_T massdepart
      REAL_T massdepart_mom

      INTEGER_T idonate,jdonate,kdonate
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T datatype

      INTEGER_T istencil
      INTEGER_T maskcell
      REAL_T donate_data
      REAL_T donate_data_MAC(SDIM,num_MAC_vectors)
      REAL_T donate_density
      REAL_T donate_mom_density
      REAL_T donate_slope,donate_cen
      REAL_T ETcore

      INTEGER_T idonate_MAC,jdonate_MAC,kdonate_MAC
      INTEGER_T icrse,jcrse,kcrse
      INTEGER_T ipart,jpart,kpart

      INTEGER_T idonatelow
      INTEGER_T idonatehigh

      INTEGER_T itensor_base,imof_base,iLS_base,iFtarget_base
      INTEGER_T iden_mom_base

      REAL_T voltotal_target
      REAL_T voltotal_depart
      REAL_T LS_voltotal_depart
      REAL_T moment_grid_diff_face

      REAL_T mofdata_grid(recon_ncomp)
      REAL_T snew_hold(ncomp_state)
      REAL_T tennew_hold(ntensor)
      REAL_T dencore(nmat)
      REAL_T newLS(nmat)
      REAL_T newvfrac_weymouth(nmat)
      REAL_T newvfrac_cor(nmat)
      REAL_T newvfrac(nmat)
      REAL_T moment_grid_diff(nmat)
      REAL_T volmat_depart(nmat)
      REAL_T volmat_target(nmat)
      REAL_T volmat_depart_cor(nmat)
      REAL_T volmat_target_cor(nmat)
      REAL_T multi_volume(nmat)
      REAL_T multi_volume_grid(nmat)
      REAL_T multi_cen(SDIM,nmat)
      REAL_T multi_cen_grid(SDIM,nmat)
      REAL_T newcen(SDIM,nmat)
      REAL_T veldata(nc_bucket)
      REAL_T veldata_MAC(SDIM,num_MAC_vectors)
      REAL_T veldata_MAC_mass(SDIM,num_MAC_vectors)

      REAL_T, dimension(:,:), allocatable :: compareconserve
      REAL_T, dimension(:,:), allocatable :: comparestate

      INTEGER_T nten_test

      INTEGER_T ihalf
      INTEGER_T nhalf
      INTEGER_T check_intersection
      INTEGER_T check_accept
      REAL_T xsten_recon(-1:1,SDIM)

      INTEGER_T ivec

      REAL_T warning_cutoff
      INTEGER_T momcomp
      REAL_T u_minimum,u_maximum

      REAL_T cutoff,DXMAXLS
      INTEGER_T all_incomp
      REAL_T divuterm,vol_target_local
      INTEGER_T k1lo,k1hi

      REAL_T massfrac_parm(num_species_var+1)

      INTEGER_T caller_id
    
! VFRAC_SPLIT code starts here


      if ((tid.lt.0).or. &
          (tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      nmax=POLYGON_LIST_MAX ! in: VFRAC_SPLIT

      k1lo=0
      k1hi=0
      if (SDIM.eq.2) then
       ! do nothing
      else if (SDIM.eq.3) then
       k1lo=-1
       k1hi=1
      else
       print *,"dimension bust"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid70"
       stop
      endif
      if (bfact_f.lt.1) then
       print *,"bfact_f invalid"
       stop
      endif
      if ((bfact.ne.bfact_f).and.(bfact.ne.2*bfact_f)) then
       print *,"bfact invalid71"
       stop
      endif
      if (nrefine_vof.ne.nmat*2*SDIM) then
       print *,"nrefine_vof invalid"
       stop
      endif
      if ((solidheat_flag.lt.0).or. &
          (solidheat_flag.gt.2)) then
       print *,"solidheat_flag invalid"
       stop
      endif

      if ((level.lt.0).or. &
          (level.gt.finest_level)) then
       print *,"level invalid vfrac split"
       stop
      endif
      if ((verbose.lt.0).or.(verbose.gt.2)) then
       print *,"verbose invalid"
       stop
      endif
      if ((gridno.lt.0).or.(gridno.ge.ngrid)) then
       print *,"gridno invalid in vfrac split"
       stop
      endif

      if ((face_flag.lt.0).or.(face_flag.gt.1)) then
       print *,"face_flag invalid F4"
       stop
      endif

      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid vfrac_split nten, nten_test ",nten,nten_test
       stop
      endif

      if ((ngrow.ne.2).and.(ngrow.ne.3)) then
       print *,"ngrow invalid"
       stop
      endif

      if ((ngrow_mac_old.ne.0).and. &
          (ngrow_mac_old.ne.2)) then
       print *,"ngrow_mac_old invalid"
       stop
      endif
 
      if (ncomp_state.ne.(SDIM+1)+ &
          nmat*(num_state_material+ngeom_raw)+1) then
       print *,"ncomp_state invalid"
       stop
      endif

      all_incomp=1

      do im=1,nmat

       if (fort_material_type(im).eq.0) then
        ! do nothing
       else if (fort_material_type(im).eq.999) then
        ! do nothing
       else if ((fort_material_type(im).ge.1).and. &
                (fort_material_type(im).le.MAX_NUM_EOS)) then
        all_incomp=0
       else
        print *,"fort_material_type invalid"
        stop
       endif

       if ((density_floor(im).lt.zero).or. &
           (density_floor(im).ge.fort_denconst(im))) then
        print *,"density_floor invalid"
        stop
       endif
       if ((density_ceiling(im).le.zero).or. &
           (density_ceiling(im).lt.fort_denconst(im))) then
        print *,"density_ceiling invalid"
        stop
       endif
       if ((temperature_primitive_variable(im).ne.0).and. &
           (temperature_primitive_variable(im).ne.1)) then
        print *,"temperature_primitive_variable invalid"
        stop
       endif
       if (added_weight(im).le.zero) then
        print *,"added_weight invalid"
        stop
       endif

       if ((fort_material_type(im).eq.0).or. &
           (is_rigid(nmat,im).eq.1).or. &
           (fort_material_type(im).eq.999)) then
        if (temperature_primitive_variable(im).ne.1) then
         print *,"temperature_primitive_variable(im) invalid"
         stop
        endif
       else if ((fort_material_type(im).gt.0).and. &
                (is_rigid(nmat,im).eq.0).and. &
                (fort_material_type(im).ne.999)) then
        if ((temperature_primitive_variable(im).eq.0).or. &
            (temperature_primitive_variable(im).eq.1)) then
         ! do nothing
        else
         print *,"temperature_primitive_variable(im) invalid"
         stop
        endif
       else
        print *,"fort_material_type(im) or is_rigid invalid"
        stop
       endif

       if ((override_density(im).ne.0).and. &
           (override_density(im).ne.1).and. &
           (override_density(im).ne.2)) then
        print *,"override_density invalid"
        stop
       endif

      enddo  ! im=1..nmat

      if (num_state_material.ne. &
          num_state_base+num_species_var) then
       print *,"num_state_material invalid"
       stop
      endif

      if ((num_materials_viscoelastic.ge.1).and. &
          (num_materials_viscoelastic.le.nmat)) then
       if ((ntensor.eq. &
            num_materials_viscoelastic*FORT_NUM_TENSOR_TYPE+SDIM).or. &
           (ntensor.eq. &
            num_materials_viscoelastic*FORT_NUM_TENSOR_TYPE)) then
        ! do nothing
       else
        print *,"ntensor invalid"
        stop
       endif
       if (ntensor.eq.NUM_CELL_ELASTIC) then
        ! do nothing
       else
        print *,"ntensor invalid"
        stop
       endif
      else
       print *,"num_materials_viscoelastic invalid"
       stop
      endif

      if ((EILE_flag.eq.-1).or. & ! Weymouth and Yue
          (EILE_flag.eq.1).or.  & ! EILE
          (EILE_flag.eq.2).or.  & ! always EI
          (EILE_flag.eq.3)) then  ! always LE
       ! do nothing
      else 
       print *,"EILE flag invalid"
       stop
      endif

      if ((dir_counter.lt.0).or.(dir_counter.ge.SDIM)) then
       print *,"dir_counter invalid"
       stop
      endif

      if ((normdir.lt.0).or.(normdir.ge.SDIM)) then
       print *,"normdir invalid"
       stop
      endif

      if (den_recon_ncomp.ne.nmat*num_state_material) then
       print *,"den_recon_ncomp invalid"
       stop
      endif

      if (recon_ncomp.ne.nmat*ngeom_recon) then
       print *,"recon_ncomp invalid"
       stop
      endif

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif

      if ((map_forward.ne.0).and.(map_forward.ne.1)) then
       print *,"map_forward invalid"
       stop
      endif

      if (levelrz.eq.0) then
       ! do nothing
      else if (levelrz.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension crash"
        stop
       endif
      else if (levelrz.eq.3) then
       ! do nothing
      else
       print *,"levelrz invalid vfrac split"
       stop
      endif

       ! SANITY CHECKS TO MAKE SURE THAT input FABs have the expected number of
       ! ghost cells.

       ! original data
      call checkbound(fablo,fabhi,DIMS(LS),1,-1,136)
      call checkbound(fablo,fabhi,DIMS(den),ngrow,-1,1231)
      call checkbound(fablo,fabhi,DIMS(mom_den),ngrow,-1,1231)
      call checkbound(fablo,fabhi,DIMS(tensor),1,-1,1231)
      call checkbound(fablo,fabhi,DIMS(velfab),ngrow,-1,125)
       ! slope data
      call checkbound(fablo,fabhi,DIMS(PLICSLP),ngrow,-1,132)
       ! new data
      call checkbound(fablo,fabhi,DIMS(snew),1,-1,130)
      call checkbound(fablo,fabhi,DIMS(tennew),1,-1,130)
      call checkbound(fablo,fabhi,DIMS(LSnew),1,-1,138)
       ! other vars
      call checkbound(fablo,fabhi,DIMS(ucell),ngrow-1,-1,135)
      call checkbound(fablo,fabhi,DIMS(vofls0),1,-1,135)
      call checkbound(fablo,fabhi,DIMS(mask),ngrow,-1,133)
      call checkbound(fablo,fabhi,DIMS(masknbr),ngrow,-1,134)
      call checkbound(fablo,fabhi,DIMS(unode),ngrow-1,normdir,121)

      if (dt.le.zero) then
        print *,"dt invalid"
        stop
      endif

      if (dencomp.ne.SDIM+1) then
       print *,"dencomp invalid"
       stop
      endif
      if (mofcomp.ne.dencomp+nmat*num_state_material) then
       print *,"mofcomp invalid"
       stop
      endif
      if (errcomp+1.ne.ncomp_state) then
       print *,"errcomp invalid"
       stop
      endif

      ii=0
      jj=0
      kk=0
      if (normdir.eq.0) then
       ii=1
      else if (normdir.eq.1) then
       jj=1
      else if ((normdir.eq.2).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"normdir invalid"
       stop
      endif

      growlo(3)=0
      growhi(3)=0

      if (face_flag.eq.1) then

       call checkbound(fablo,fabhi,DIMS(xvel),ngrow_mac_old,0,1243)
       call checkbound(fablo,fabhi,DIMS(yvel),ngrow_mac_old,1,125)
       call checkbound(fablo,fabhi,DIMS(zvel),ngrow_mac_old,SDIM-1,126)

       call checkbound(fablo,fabhi,DIMS(xvelslp),ngrow_mac_old,0,1243)
       call checkbound(fablo,fabhi,DIMS(yvelslp),ngrow_mac_old,1,125)
       call checkbound(fablo,fabhi,DIMS(zvelslp),ngrow_mac_old,SDIM-1,126)

       call checkbound(fablo,fabhi,DIMS(xmomside),1,-1,1271)
       call checkbound(fablo,fabhi,DIMS(ymomside),1,-1,1271)
       call checkbound(fablo,fabhi,DIMS(zmomside),1,-1,1271)

       call checkbound(fablo,fabhi,DIMS(xmassside),1,-1,1271)
       call checkbound(fablo,fabhi,DIMS(ymassside),1,-1,1271)
       call checkbound(fablo,fabhi,DIMS(zmassside),1,-1,1271)
 
       if ((NUM_CELL_ELASTIC.eq. &
            2*SDIM*num_materials_viscoelastic+SDIM).and. &
           (NUM_CELL_ELASTIC.eq. &
            FORT_NUM_TENSOR_TYPE* &
            num_materials_viscoelastic+SDIM).and. &
           (ntensor.eq.NUM_CELL_ELASTIC).and. &
           (num_MAC_vectors.eq.1)) then
        ! do nothing
       else if ((NUM_CELL_ELASTIC.eq. &
                 2*SDIM*num_materials_viscoelastic).and. &
                (NUM_CELL_ELASTIC.eq. &
                 FORT_NUM_TENSOR_TYPE* &
                 num_materials_viscoelastic).and. &
                (ntensor.eq.NUM_CELL_ELASTIC).and. &
                (num_MAC_vectors.eq.2)) then
        ! do nothing
       else
        print *,"expecting displacement at cell centers or face centers"
        stop
       endif
      else if (face_flag.eq.0) then
       if ((NUM_CELL_ELASTIC.eq. &
            2*SDIM*num_materials_viscoelastic+SDIM).and. &
           (NUM_CELL_ELASTIC.eq. &
            FORT_NUM_TENSOR_TYPE* &
            num_materials_viscoelastic+SDIM).and. &
           (ntensor.eq.NUM_CELL_ELASTIC).and. &
           (num_MAC_vectors.eq.1)) then
        ! do nothing
       else
        print *,"expecting displacement at cell centers"
        stop
       endif
      else
       print *,"face_flag invalid F6"
       stop
      endif

      if (nc_conserve.ne.SDIM+nmat*num_state_material) then
       print *,"nc_conserve invalid"
       stop
      endif
      if (iden_base.ne.SDIM) then
       print *,"iden_base invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(momslope),1,-1,1238)
      call checkbound(fablo,fabhi,DIMS(conserve),ngrow,-1,1238)
     
      if (DO_SANITY_CHECK.eq.1) then
       print *,"SANITY CHECK AFTER CONSERVE dir_counter= ",dir_counter
       print *,"SANITY CHECK AFTER CONSERVE normdir= ",normdir
       call CISL_sanity(dt,normdir,dir_counter,1)
       allocate(compareconserve(fablo(1)-1:fabhi(1)+1,5))
       allocate(comparestate(fablo(1)-1:fabhi(1)+1,5))
       jcrse=1
       kcrse=0
       do icrse=fablo(1),fabhi(1)
        compareconserve(icrse,1)= &
         conserve(D_DECL(icrse,jcrse,kcrse),1)* &
         conserve(D_DECL(icrse,jcrse,kcrse),iden_base+1)
        compareconserve(icrse,2)= &
          conserve(D_DECL(icrse,jcrse,kcrse),iden_base+1)
        compareconserve(icrse,3)= &
          conserve(D_DECL(icrse,jcrse,kcrse),iden_base+2)
        comparestate(icrse,1)=snew(D_DECL(icrse,jcrse,kcrse),1)
        comparestate(icrse,2)=snew(D_DECL(icrse,jcrse,kcrse),dencomp+1)
        comparestate(icrse,3)=snew(D_DECL(icrse,jcrse,kcrse),dencomp+2)
       enddo
       call compare_sanity(compareconserve,1,3,2)
       call compare_sanity(comparestate,1,3,1)
       deallocate(compareconserve)
       deallocate(comparestate)
       print *,"AFTER SANITY CHECK AFTER CONSERVE dir_counter= ",dir_counter
       print *,"AFTER SANITY CHECK AFTER CONSERVE normdir= ",normdir
      endif ! if do_sanity_check
    
      force_check=0
      datatype=0 
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,ngrow)
 
      warning_cutoff=two
      call aggressive_worker( &
       datatype, &
       warning_cutoff, &
       tilelo,tilehi, &
       fablo,fabhi, &
       growlo,growhi, &
       bfact, &
       dx, &
       0, &
       nmat*ngeom_recon, &
       nmat*ngeom_recon, &
       ngrow,-1,19, &
       verbose, &
       force_check, &
       gridno,ngrid, &
       level,finest_level, &
       PLICSLP,DIMS(PLICSLP))

      warning_cutoff=1.0e+15
      call aggressive_worker( &
       datatype, & 
       warning_cutoff, &
       tilelo,tilehi, &
       fablo,fabhi, &
       growlo,growhi, &
       bfact, &
       dx, &
       0, &
       nc_conserve, &
       nc_conserve, &
       ngrow,-1,20, &
       verbose, &
       force_check, &
       gridno,ngrid, &
       level, &
       finest_level, &
       conserve,DIMS(den))

      if (iden_base.ne.SDIM) then
       print *,"iden_base invalid"
       stop
      endif

       ! the band thickness in FORT_ADVECTIVE_PRESSURE is 2 * DXMAXLS
      call get_dxmaxLS(dx,bfact,DXMAXLS)
      cutoff=DXMAXLS

      itensor_base=iden_base+nmat*num_state_material
      imof_base=itensor_base+NUM_CELL_ELASTIC
      iLS_base=imof_base+nmat*ngeom_raw
      iFtarget_base=iLS_base+nmat
      iden_mom_base=iFtarget_base+nmat
      nc_bucket_test=iden_mom_base+nmat
      if (nc_bucket_test.ne.nc_bucket) then
       print *,"nc_bucket invalid"
       stop
      endif

      if (face_flag.eq.1) then

        ! xmomside,ymomside,zmomside already init to 0 
        ! xmassside,ymassside,zmassside already init to 0 
       do veldir=1,SDIM

        iii=0
        jjj=0
        kkk=0
        if (veldir.eq.1) then
         iii=1
        else if (veldir.eq.2) then
         jjj=1
        else if ((veldir.eq.3).and.(SDIM.eq.3)) then
         kkk=1
        else
         print *,"veldir invalid"
         stop
        endif

        call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
         growlo,growhi,0,veldir-1,32)

         ! directionally_split_module vars:
         ! icrse,jcrse,kcrse,iside is a MAC index.
         ! ipart,jpart,kpart,iside_part is a CELL index.
        do icrse=growlo(1),growhi(1)
        do jcrse=growlo(2),growhi(2)
        do kcrse=growlo(3),growhi(3)
         nprocessed=nprocessed+1

         nhalf=1
         call gridstenMAC_level(xsten_MAC,icrse,jcrse,kcrse, &
           level,nhalf,veldir-1,35)

         check_accept=1

         if (levelrz.eq.0) then
          ! do nothing
         else if (levelrz.eq.1) then
          if (SDIM.ne.2) then
           print *,"dimension bust"
           stop
          endif
          if ((xsten_MAC(0,1).le.VOFTOL*dx(1)).and. &
              (veldir.eq.1)) then
           check_accept=0
          endif
         else if (levelrz.eq.3) then
          if ((xsten_MAC(0,1).le.VOFTOL*dx(1)).and. &
              (veldir.eq.1)) then
           check_accept=0
          endif
         else
          print *,"levelrz invalid vfrac split 2"
          stop
         endif

         if (check_accept.eq.1) then

           do iside=-1,1,2

            iside_part=-iside

             ! control volume left of face
            if (iside.eq.-1) then
             ipart=icrse-iii 
             jpart=jcrse-jjj 
             kpart=kcrse-kkk
 
             ! control volume right of face
            else if (iside.eq.1) then
             ipart=icrse
             jpart=jcrse
             kpart=kcrse
            else
             print *,"iside invalid"
             stop
            endif

            nhalf=1
            call CISBOXHALF(xsten_accept,nhalf, &
             xlo,dx,ipart,jpart,kpart,iside_part,veldir, &
             bfact,level, &
             volcell_accept,cencell_accept,SDIM)
 
            if (volcell_accept.le.zero) then
             print *,"volcell_accept invalid"
             stop
            endif

            do ihalf=-1,1
            do dir2=1,SDIM
             xsten_target(ihalf,dir2)=xsten_accept(ihalf,dir2)
             xsten_depart(ihalf,dir2)=xsten_accept(ihalf,dir2)
            enddo
            enddo

            ! veldata(momcomp)  momcomp=(im-1)*sdim+veldir
            ! veldata(iden_mom_base+im)

            usten_accept(-1)=unode(D_DECL(ipart,jpart,kpart))
            nhalf=1
             ! normdir=0..sdim-1
            call gridstenMAC_level(xsten_MAC,ipart,jpart,kpart, &
              level,nhalf,normdir,36)

            if (levelrz.eq.0) then
             ! do nothing
            else if ((levelrz.eq.1).or. &
                     (levelrz.eq.3)) then
             if ((xsten_MAC(0,1).le.VOFTOL*dx(1)).and. &
                 (normdir.eq.0)) then
              usten_accept(-1)=zero
             endif
            else
             print *,"levelrz invalid add to bucket mac"
             stop
            endif
            usten_accept(1)=unode(D_DECL(ipart+ii,jpart+jj,kpart+kk))

            if (usten_accept(-1).gt.zero) then
             idonatelow=-1
            else
             idonatelow=0
            endif
            if (usten_accept(1).lt.zero) then
             idonatehigh=1
            else
             idonatehigh=0
            endif

            if (veldir.eq.normdir+1) then
             veldir_comp=normdir+1
             usten_accept(0)=ucell(D_DECL(ipart,jpart,kpart),veldir_comp)
             u_minimum=min(usten_accept(-1),usten_accept(1))
             u_maximum=max(usten_accept(-1),usten_accept(1))
             if (usten_accept(0).lt.u_minimum) then
              usten_accept(0)=u_minimum
             endif
             if (usten_accept(0).gt.u_maximum) then
              usten_accept(0)=u_maximum
             endif
      
             if (iside_part.eq.-1) then
              usten_accept(1)=usten_accept(0)
              usten_accept(0)=half*(usten_accept(-1)+usten_accept(1))
             else if (iside_part.eq.1) then
              usten_accept(-1)=usten_accept(0)
              usten_accept(0)=half*(usten_accept(-1)+usten_accept(1))
             else
              print *,"iside_part invalid"
              stop
             endif
            else
             usten_accept(0)=half*(usten_accept(-1)+usten_accept(1))
            endif

            idonate=ipart
            jdonate=jpart
            kdonate=kpart

            do istate=1,nc_bucket
             veldata(istate)=zero
            enddo
            do ivec=1,num_MAC_vectors
             veldata_MAC(veldir,ivec)=zero
            enddo

            do istencil=idonatelow,idonatehigh

             if (normdir.eq.0) then
              idonate=ipart+istencil
             else if (normdir.eq.1) then
              jdonate=jpart+istencil
             else if ((normdir.eq.2).and.(SDIM.eq.3)) then
              kdonate=kpart+istencil
             else
              print *,"normdir invalid"
              stop
             endif

             iside_low=iside_part
             iside_high=iside_part
             if (veldir.eq.normdir+1) then
              iside_low=-1
              iside_high=1
             endif

             do isidedonate=iside_low,iside_high,2

              call CISBOX(xsten_recon,1, &
               xlo,dx,idonate,jdonate,kdonate, &
               bfact,level, &
               volcell_recon,cencell_recon,SDIM)

              call CISBOXHALF(xsten_donate,1, &
               xlo,dx,idonate,jdonate,kdonate,isidedonate,veldir, &
               bfact,level, &
               volcell_donate,cencell_donate,SDIM)

              check_intersection=1

              if (levelrz.eq.0) then
               ! do nothing
              else if (levelrz.eq.1) then
               if (SDIM.ne.2) then
                print *,"dimension bust"
                stop
               endif
               if (xsten_recon(0,1).le.VOFTOL*dx(1)) then
                check_intersection=0
               endif
              else if (levelrz.eq.3) then
               if (xsten_recon(0,1).le.VOFTOL*dx(1)) then
                check_intersection=0
               endif
              else
               print *,"levelrz invalid add to bucket MAC 2"
               stop
              endif

              if (check_intersection.eq.1) then 

               usten_donate(-1)=unode(D_DECL(idonate,jdonate,kdonate))

               nhalf=1
               call gridstenMAC_level(xsten_MAC,idonate,jdonate,kdonate, &
                level,nhalf,normdir,37)

               if (levelrz.eq.0) then
                ! do nothing
               else if ((levelrz.eq.1).or. &
                        (levelrz.eq.3)) then
                if ((xsten_MAC(0,1).le.VOFTOL*dx(1)).and. &
                    (normdir.eq.0)) then
                 usten_donate(-1)=zero
                endif
               else
                print *,"levelrz invalid add to bucket mac 3"
                stop
               endif

               usten_donate(1)=unode(D_DECL(idonate+ii,jdonate+jj,kdonate+kk))

               if (veldir.eq.normdir+1) then

                veldir_comp=normdir+1
                usten_donate(0)= &
                 ucell(D_DECL(idonate,jdonate,kdonate),veldir_comp)
                u_minimum=min(usten_donate(-1),usten_donate(1))
                u_maximum=max(usten_donate(-1),usten_donate(1))
                if (usten_donate(0).lt.u_minimum) then
                 usten_donate(0)=u_minimum
                endif
                if (usten_donate(0).gt.u_maximum) then
                 usten_donate(0)=u_maximum
                endif

                if (isidedonate.eq.-1) then
                 usten_donate(1)=usten_donate(0)
                 usten_donate(0)= &
                  half*(usten_donate(-1)+usten_donate(1))
                else if (isidedonate.eq.1) then
                 usten_donate(-1)=usten_donate(0)
                 usten_donate(0)= &
                  half*(usten_donate(-1)+usten_donate(1))
                else
                 print *,"isidedonate invalid"
                 stop
                endif
               else
                usten_donate(0)= &
                 half*(usten_donate(-1)+usten_donate(1))
               endif

               call derive_mappings( &
                xsten_accept, &
                xsten_donate, &
                xsten_target, &
                xsten_depart, &
                usten_accept, &
                usten_donate, &
                xdepartsize, &
                xtargetsize, &
                xloint, &
                xhiint, &
                volint, &
                coeff, &
                bfact,dx,map_forward,normdir)

               if (volint.gt.zero) then  

                do dir2=1,nmat*ngeom_recon
                 mofdata_grid(dir2)= &
                  PLICSLP(D_DECL(idonate,jdonate,kdonate),dir2)
                enddo
 
                ! find departure volume within xsten_xrecon
                tessellate=0
                if (nmax.lt.10) then
                 print *,"nmax bust 3"
                 stop
                endif
                call multi_get_volume_grid_simple( &
                  tessellate, &  !=0
                  bfact,dx, &
                  xsten_recon,1, &
                  mofdata_grid, &
                  xsten_depart,1, &
                  multi_volume_grid,multi_cen_grid, &
                  geom_xtetlist_uncapt(1,1,1,tid+1), &
                  nmax, &
                  nmax, &
                  nmat,SDIM,3)

                 ! we are inside the istencil + isidedonate loops.
                LS_voltotal_depart=zero
                do im=1,nmat
                 if (is_rigid(nmat,im).eq.0) then 
                  LS_voltotal_depart=LS_voltotal_depart+ &
                   multi_volume_grid(im)
                 else if (is_rigid(nmat,im).eq.1) then
                  ! do nothing (fluids tessellate)
                 else
                  print *,"is_rigid invalid"
                  stop
                 endif
                enddo !im=1..nmat

                if (LS_voltotal_depart.le.zero) then
                 print *,"LS_voltotal_depart bust "
                 print *,"LS_voltotal_depart ",LS_voltotal_depart
                 stop
                endif

                 ! initialize momentum for each material.
                do im=1,nmat

                 ! density
                 dencomp_data=(im-1)*num_state_material+1
                 donate_data= &
                  conserve(D_DECL(idonate,jdonate,kdonate), &
                           iden_base+dencomp_data) 
                 massdepart_mom=donate_data*added_weight(im)

                 if (massdepart_mom.le.zero) then
                  print *,"density invalid in vfrac split 1"
                  print *,"idonate,jdonate,kdonate ",idonate,jdonate,kdonate
                  print *,"im ",im
                  print *,"donate_data ",donate_data
                  print *,"multi_volume_grid(im) ",multi_volume_grid(im)
                  print *,"massdepart_mom ",massdepart_mom
                  stop
                 endif

                 mom2(veldir)=zero 

                  !left half of cell (idonate,jdonate,kdonate)
                 if (isidedonate.eq.-1) then
                  idonate_MAC=idonate 
                  jdonate_MAC=jdonate 
                  kdonate_MAC=kdonate 

                  !right half of cell (idonate,jdonate,kdonate)
                 else if (isidedonate.eq.1) then
                  idonate_MAC=idonate+iii 
                  jdonate_MAC=jdonate+jjj 
                  kdonate_MAC=kdonate+kkk 
                 else
                  print *,"isidedonate invalid"
                  stop
                 endif

                 ! xvel   : xvel  1..num_MAC_vectors
                 ! xvelslp: xvelslope,xcen   1..1+nmat
                 
                 if (veldir.eq.1) then
                  do ivec=1,num_MAC_vectors
                   donate_data_MAC(veldir,ivec)= &
                    xvel(D_DECL(idonate_MAC,jdonate_MAC,kdonate_MAC),ivec) 
                  enddo
                  donate_cen= &
                    xvelslp(D_DECL(idonate_MAC,jdonate_MAC,kdonate_MAC),1+im) 
                  donate_slope= &
                    xvelslp(D_DECL(idonate_MAC,jdonate_MAC,kdonate_MAC),1) 
                 else if (veldir.eq.2) then
                  do ivec=1,num_MAC_vectors
                   donate_data_MAC(veldir,ivec)= &
                    yvel(D_DECL(idonate_MAC,jdonate_MAC,kdonate_MAC),ivec) 
                  enddo
                  donate_cen= &
                    yvelslp(D_DECL(idonate_MAC,jdonate_MAC,kdonate_MAC),1+im) 
                  donate_slope= &
                    yvelslp(D_DECL(idonate_MAC,jdonate_MAC,kdonate_MAC),1) 
                 else if ((veldir.eq.3).and.(SDIM.eq.3)) then
                  do ivec=1,num_MAC_vectors
                   donate_data_MAC(veldir,ivec)= &
                    zvel(D_DECL(idonate_MAC,jdonate_MAC,kdonate_MAC),ivec) 
                  enddo
                  donate_cen= &
                    zvelslp(D_DECL(idonate_MAC,jdonate_MAC,kdonate_MAC),1+im) 
                  donate_slope= &
                    zvelslp(D_DECL(idonate_MAC,jdonate_MAC,kdonate_MAC),1) 
                 else
                  print *,"veldir invalid"
                  stop
                 endif
                 moment_grid_diff_face= &
                  multi_cen_grid(normdir+1,im)-donate_cen 

                 ! a few lines before:
                 !   massdepart_mom=donate_data*added_weight(im)
                 massdepart_mom=massdepart_mom*multi_volume_grid(im)

                 mom2(veldir)=massdepart_mom* &
                  (donate_data_MAC(veldir,1)+ &
                   donate_slope*moment_grid_diff_face)

                 veldata(iden_mom_base+im)= &
                  veldata(iden_mom_base+im)+massdepart_mom

                 momcomp=(im-1)*SDIM+veldir
                 veldata(momcomp)=veldata(momcomp)+mom2(veldir) 

                  ! this gets incremented for im=1..nmat
                 veldata_MAC_mass(veldir,1)= &
                  veldata_MAC_mass(veldir,1)+massdepart_mom

                  ! this gets incremented for im=1..nmat
                 veldata_MAC(veldir,1)= &
                  veldata_MAC(veldir,1)+mom2(veldir)

                 if ((num_materials_viscoelastic.ge.1).and. &
                     (num_materials_viscoelastic.le.nmat)) then

                  if (fort_store_elastic_data(im).eq.1) then
                   imap=1
                   do while ((fort_im_elastic_map(imap)+1.ne.im).and. &
                             (imap.le.num_materials_viscoelastic))
                    imap=imap+1
                   enddo
                   if ((imap.ge.1).and. &
                       (imap.le.num_materials_viscoelastic)) then
                    if (num_MAC_vectors.eq.1) then
                     ! do nothing
                    else if (num_MAC_vectors.eq.2) then
                      ! only increment displacement ONCE.
                     if (imap.eq.1) then
                      veldata_MAC_mass(veldir,2)= &
                       veldata_MAC_mass(veldir,2)+LS_voltotal_depart
                      veldata_MAC(veldir,2)= &
                       veldata_MAC(veldir,2)+ &
                       donate_data_MAC(veldir,2)*LS_voltotal_depart
                     else if ((imap.ge.2).and. &
                              (imap.le.num_materials_viscoelastic)) then
                      ! do nothing
                     else
                      print *,"imap invalid"
                      stop
                     endif
                    else
                     print *,"num_MAC_vectors invalid"
                     stop
                    endif
                   else 
                    print *,"imap invalid"
                    stop
                   endif
                  else if (fort_store_elastic_data(im).eq.0) then
                   ! do nothing
                  else
                   print *,"fort_store_elastic_data(im) invalid"
                   stop
                  endif

                 else
                  print *,"num_materials_viscoelastic invalid"
                  stop
                 endif

                enddo ! im=1,..,nmat 

               else if (volint.eq.zero) then
                ! do nothing
               else
                print *,"volint bust"
                stop
               endif 

              else if (check_intersection.eq.0) then
               ! do nothing
              else
               print *,"check_intersection invalid"
               stop
              endif

             enddo  ! isidedonate=iside_low,iside_high,2

            enddo  ! istencil=idonatelow,idonatehigh

            do im=1,nmat

              ! left of cell
             if (iside_part.eq.-1) then
              ibucket=1

              ! right of cell
             else if (iside_part.eq.1) then
              ibucket=2
             else
              print *,"iside_part invalid"
              stop
             endif

             momcomp=(im-1)*SDIM+veldir

             if (is_prescribed(nmat,im).eq.0) then

              if ((num_materials_viscoelastic.ge.1).and. &
                  (num_materials_viscoelastic.le.nmat)) then

               if (fort_store_elastic_data(im).eq.1) then
                imap=1
                do while ((fort_im_elastic_map(imap)+1.ne.im).and. &
                          (imap.le.num_materials_viscoelastic))
                 imap=imap+1
                enddo

                if ((imap.ge.2).and. & 
                    (imap.le.num_materials_viscoelastic)) then
                 ! do nothing

                 ! increment just ONCE. (since using veldata_MAC)
                else if (imap.eq.1) then
                 if (num_MAC_vectors.eq.1) then
                  ! do nothing
                 else if (num_MAC_vectors.eq.2) then
                  ibucket_map=ibucket+2
                  if (veldir.eq.1) then
                   xmomside(D_DECL(ipart,jpart,kpart),ibucket_map)= &
                    xmomside(D_DECL(ipart,jpart,kpart),ibucket_map)+ &
                    veldata_MAC(veldir,2)
                   xmassside(D_DECL(ipart,jpart,kpart),ibucket_map)= &
                    xmassside(D_DECL(ipart,jpart,kpart),ibucket_map)+ &
                    veldata_MAC_mass(veldir,2)
                  else if (veldir.eq.2) then
                   ymomside(D_DECL(ipart,jpart,kpart),ibucket_map)= &
                    ymomside(D_DECL(ipart,jpart,kpart),ibucket_map)+ &
                    veldata_MAC(veldir,2)
                   ymassside(D_DECL(ipart,jpart,kpart),ibucket_map)= &
                    ymassside(D_DECL(ipart,jpart,kpart),ibucket_map)+ &
                    veldata_MAC_mass(veldir,2)
                  else if ((veldir.eq.3).and.(SDIM.eq.3)) then
                   zmomside(D_DECL(ipart,jpart,kpart),ibucket_map)= &
                    zmomside(D_DECL(ipart,jpart,kpart),ibucket_map)+ &
                    veldata_MAC(veldir,2)
                   zmassside(D_DECL(ipart,jpart,kpart),ibucket_map)= &
                    zmassside(D_DECL(ipart,jpart,kpart),ibucket_map)+ &
                    veldata_MAC_mass(veldir,2)
                  else
                   print *,"veldir invalid"
                   stop
                  endif

                 else
                  print *,"num_MAC_vectors invalid"
                  stop
                 endif
                else 
                 print *,"imap invalid"
                 stop
                endif
               else if (fort_store_elastic_data(im).eq.0) then
                ! do nothing
               else
                print *,"fort_store_elastic_data(im) invalid"
                stop
               endif

              else
               print *,"num_materials_viscoelastic invalid"
               stop
              endif
         
               ! increment for each material:im=1..nmat (since using veldata)
              if (veldir.eq.1) then
               xmomside(D_DECL(ipart,jpart,kpart),ibucket)= &
                xmomside(D_DECL(ipart,jpart,kpart),ibucket)+ &
                veldata(momcomp)
               xmassside(D_DECL(ipart,jpart,kpart),ibucket)= &
                xmassside(D_DECL(ipart,jpart,kpart),ibucket)+ &
                veldata(iden_mom_base+im)
              else if (veldir.eq.2) then
               ymomside(D_DECL(ipart,jpart,kpart),ibucket)= &
                ymomside(D_DECL(ipart,jpart,kpart),ibucket)+ &
                veldata(momcomp)
               ymassside(D_DECL(ipart,jpart,kpart),ibucket)= &
                ymassside(D_DECL(ipart,jpart,kpart),ibucket)+ &
                veldata(iden_mom_base+im)
              else if ((veldir.eq.3).and.(SDIM.eq.3)) then
               zmomside(D_DECL(ipart,jpart,kpart),ibucket)= &
                zmomside(D_DECL(ipart,jpart,kpart),ibucket)+ &
                veldata(momcomp)
               zmassside(D_DECL(ipart,jpart,kpart),ibucket)= &
                zmassside(D_DECL(ipart,jpart,kpart),ibucket)+ &
                veldata(iden_mom_base+im)
              else
               print *,"veldir invalid"
               stop
              endif
     
             else if (is_prescribed(nmat,im).eq.1) then
              ! do nothing
             else
              print *,"is_prescribed(nmat,im) invalid"
              stop
             endif

            enddo ! im=1..nmat

           enddo ! iside=-1,1,2

         else if (check_accept.eq.0) then
          ! do nothing (buckets already init to 0)
         else
          print *,"check_accept bust"
          stop
         endif
        enddo
        enddo
        enddo  ! icrse,jcrse,kcrse -> growntileboxMAC(0 ghost)

       enddo ! veldir=1..sdim

      else if (face_flag.eq.0) then
       ! do nothing
      else
       print *,"face_flag invalid"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi, &
       growlo,growhi,0)

      do icrse=growlo(1),growhi(1)
      do jcrse=growlo(2),growhi(2)
      do kcrse=growlo(3),growhi(3)
       nprocessed=nprocessed+1

       maskcell=NINT(mask(D_DECL(icrse,jcrse,kcrse)))

       if (maskcell.eq.1) then

        nhalf=1
        call gridsten_level(xsten_crse,icrse,jcrse,kcrse,level,nhalf)

        check_accept=1
        if (levelrz.eq.0) then
         ! do nothing
        else if (levelrz.eq.1) then
         if (SDIM.ne.2) then
          print *,"dimension bust"
          stop
         endif
         if (xsten_crse(0,1).le.VOFTOL*dx(1)) then
          check_accept=0
          print *,"icrse invalid"
          stop
         endif
        else if (levelrz.eq.3) then
         if (xsten_crse(0,1).le.VOFTOL*dx(1)) then
          check_accept=0
          print *,"icrse invalid"
          stop
         endif
        else
         print *,"levelrz invalid vfrac split 3"
         stop
        endif
 
        if (check_accept.eq.1) then 

         voltotal_depart=zero
         do istate=1,nc_bucket
          veldata(istate)=zero
         enddo

         nhalf=1
         call CISBOX(xsten_accept,nhalf, &
          xlo,dx,icrse,jcrse,kcrse, &
          bfact,level, &
          volcell_accept,cencell_accept,SDIM)

         if (volcell_accept.le.zero) then
          print *,"volcell_accept invalid"
          stop
         endif

         do ihalf=-1,1
         do dir2=1,SDIM
          xsten_target(ihalf,dir2)=xsten_accept(ihalf,dir2)
          xsten_depart(ihalf,dir2)=xsten_accept(ihalf,dir2)
         enddo
         enddo

         usten_accept(-1)=unode(D_DECL(icrse,jcrse,kcrse))

         nhalf=1
         call gridstenMAC_level(xsten_MAC,icrse,jcrse,kcrse, &
           level,nhalf,normdir,38)

         if (levelrz.eq.0) then
          ! do nothing
         else if ((levelrz.eq.1).or. &
                  (levelrz.eq.3)) then
          if ((xsten_MAC(0,1).le.VOFTOL*dx(1)).and. &
              (normdir.eq.0)) then
           usten_accept(-1)=zero
          endif
         else
          print *,"levelrz invalid add to bucket"
          stop
         endif

         usten_accept(1)=unode(D_DECL(icrse+ii,jcrse+jj,kcrse+kk))

         if (usten_accept(-1).gt.zero) then
          idonatelow=-1
         else if (usten_accept(-1).le.zero) then
          idonatelow=0
         else
          print *,"usten_accept(-1) invalid"
          stop
         endif

         if (usten_accept(1).lt.zero) then
          idonatehigh=1
         else if (usten_accept(1).ge.zero) then
          idonatehigh=0
         else
          print *,"usten_accept(1) invalid"
          stop
         endif

         usten_accept(0)=half*(usten_accept(-1)+usten_accept(1))

         idonate=icrse
         jdonate=jcrse
         kdonate=kcrse

         null_velocity_flag=0

         if ((usten_accept(1).eq.zero).and. &
             (usten_accept(-1).eq.zero)) then
          null_velocity_flag=1
          if ((idonatelow.eq.0).and.(idonatehigh.eq.0)) then
           ! do nothing
          else
           print *,"idonatelow or idonatehigh invalid"
           stop
          endif
         else if ((usten_accept(1).ne.zero).or. &
                  (usten_accept(-1).ne.zero)) then
          ! do nothing
         else
          print *,"usten_accept is corrupt"
          stop
         endif

         do istencil=idonatelow,idonatehigh
 
          if (normdir.eq.0) then
           idonate=icrse+istencil
          else if (normdir.eq.1) then 
           jdonate=jcrse+istencil
          else if ((normdir.eq.2).and.(SDIM.eq.3)) then
           kdonate=kcrse+istencil
          else
           print *,"normdir invalid"
           stop
          endif
      
          call CISBOX(xsten_recon,1, &
           xlo,dx,idonate,jdonate,kdonate, &
           bfact,level, &
           volcell_recon,cencell_recon,SDIM)
 
          call CISBOX(xsten_donate,1, &
           xlo,dx,idonate,jdonate,kdonate, &
           bfact,level, &
           volcell_donate,cencell_donate,SDIM)

          check_intersection=1

          if (levelrz.eq.0) then
           ! do nothing
          else if (levelrz.eq.1) then
           if (SDIM.ne.2) then
            print *,"dimension bust"
            stop
           endif
           if (xsten_recon(0,1).le.VOFTOL*dx(1)) then
            check_intersection=0
            print *,"idonate invalid"
            stop
           endif
          else if (levelrz.eq.3) then
           if (xsten_recon(0,1).le.VOFTOL*dx(1)) then
            check_intersection=0
            print *,"idonate invalid"
            stop
           endif
          else
           print *,"levelrz invalid add to bucket 2"
           stop
          endif

          if (check_intersection.eq.1) then 

           usten_donate(-1)=unode(D_DECL(idonate,jdonate,kdonate))

           nhalf=1
           call gridstenMAC_level(xsten_MAC,idonate,jdonate,kdonate, &
             level,nhalf,normdir,39)

           if (levelrz.eq.0) then
            ! do nothing
           else if ((levelrz.eq.1).or. &
                    (levelrz.eq.3)) then

            if ((xsten_MAC(0,1).le.VOFTOL*dx(1)).and. &
                (normdir.eq.0)) then
             usten_donate(-1)=zero
            endif

           else
            print *,"levelrz invalid add to bucket 3"
            stop
           endif

           usten_donate(1)=unode(D_DECL(idonate+ii,jdonate+jj,kdonate+kk))

           usten_donate(0)=half*(usten_donate(-1)+usten_donate(1))

             ! normdir=0..sdim-1
           call derive_mappings( &
            xsten_accept, &
            xsten_donate, &
            xsten_target, &
            xsten_depart, &
            usten_accept, &
            usten_donate, &
            xdepartsize, &
            xtargetsize, &
            xloint, &
            xhiint, &
            volint, &
            coeff, &
            bfact,dx,map_forward,normdir)

           if (volint.gt.zero) then  

             ! we are inside the istencil loop.
            LS_voltotal_depart=zero

            do dir2=1,nmat*ngeom_recon
             mofdata_grid(dir2)= &
              PLICSLP(D_DECL(idonate,jdonate,kdonate),dir2)
            enddo

              ! the volumes and centroids are tessellating for the fluid
              ! materials, but not the solid materials.  Solid materials are
              ! immersed into the domain.

            if ((istencil.eq.-1).and. &
                (icrse.eq.11).and.(jcrse.eq.16).and.(kcrse.eq.12).and. &
                (level.eq.1).and. &
                (1.eq.0)) then
             caller_id=-1
            else
             caller_id=3
            endif
             
            call multi_get_volume_grid_and_map( &
              normdir, & ! normdir=0..sdim-1
              coeff, &
              bfact,dx, &
              xsten_recon,1, &
              mofdata_grid, &
              xsten_depart,1, &
              multi_volume_grid, & ! intersection of departure with grid.
              multi_cen_grid, &
              multi_volume, & ! intersection of target with grid.
              multi_cen, &
              geom_xtetlist_uncapt(1,1,1,tid+1), &
              nmax, &
              nmax, &
              nmat,SDIM,caller_id)

            if (null_velocity_flag.eq.1) then
             if (istencil.eq.0) then
              do im=1,nmat
               vofcomp=(im-1)*ngeom_recon+1
               multi_volume_grid(im)=mofdata_grid(vofcomp)*volcell_recon
               multi_volume(im)=multi_volume_grid(im)
               do dir2=1,SDIM
                multi_cen_grid(dir2,im)=cencell_recon(dir2)+ &
                    mofdata_grid(vofcomp+dir2)
                multi_cen(dir2,im)=multi_cen_grid(dir2,im)
               enddo ! dir2=1..sdim
              enddo !im=1..nmat
             else
              print *,"istencil invalid"
              stop
             endif
            else if (null_velocity_flag.eq.0) then
             ! do nothing
            else
             print *,"null_velocity_flag invalid"
             stop
            endif 

             ! normdir=0..sdim-1
            do im=1,nmat
             vofcomp=(im-1)*ngeom_recon+1
             moment_grid_diff(im)= &
              multi_cen_grid(normdir+1,im)- &
              mofdata_grid(vofcomp+normdir+1)- &
              cencell_recon(normdir+1) 
  
              ! fluid materials tessellate the domain. 
             if (is_rigid(nmat,im).eq.0) then 
              LS_voltotal_depart=LS_voltotal_depart+ &
               multi_volume_grid(im)
             else if (is_rigid(nmat,im).eq.1) then
              ! do nothing
             else
              print *,"is_rigid invalid"
              stop
             endif

            enddo  ! im=1,..,nmat

            if (LS_voltotal_depart.le.zero) then
             print *,"LS_voltotal_depart bust "
             print *,"LS_voltotal_depart ",LS_voltotal_depart
             print *,"map_forward,volint ",map_forward,volint
             print *,"istencil ",istencil
             print *,"icrse,jcrse,kcrse ",icrse,jcrse,kcrse
             print *,"level,finest_level ",level,finest_level
             do im=1,nmat 
              vofcomp=(im-1)*ngeom_recon+1
              print *,"im,multi_volume_grid ",im,multi_volume_grid(im)
              print *,"im,multi_volume ",im,multi_volume(im)
              print *,"im,vfrac ",im,mofdata_grid(vofcomp)
              print *,"im,flag ",im,mofdata_grid(vofcomp+SDIM+1)
             enddo
             stop
            endif
            voltotal_depart=voltotal_depart+ &
             LS_voltotal_depart

            do im=1,nmat
             ! density
             dencomp_data=(im-1)*num_state_material+1
              ! conserve is initialized in FORT_BUILD_CONSERVE.
              ! donate_density is equal to the density that is stored in the
              ! old state variable.
             donate_density= &
              conserve(D_DECL(idonate,jdonate,kdonate),iden_base+dencomp_data) 
             donate_mom_density= &
              mom_den(D_DECL(idonate,jdonate,kdonate),im) 
             if (donate_density.gt.zero) then
              ! do nothing
             else
              print *,"donate_density must be positive"
              stop
             endif
             if (donate_mom_density.gt.zero) then
              ! do nothing
             else
              print *,"donate_mom_density must be positive"
              stop
             endif
             donate_slope= &
              momslope(D_DECL(idonate,jdonate,kdonate),iden_base+dencomp_data)
             massdepart=donate_density+donate_slope*moment_grid_diff(im)
             massdepart_mom=donate_density*added_weight(im)

             if (massdepart.le.zero) then
              print *,"density invalid in vfrac split 2"
              print *,"icrse,jcrse,kcrse ",icrse,jcrse,kcrse
              print *,"idonate,jdonate,kdonate ",idonate,jdonate,kdonate
              print *,"im ",im
              print *,"donate_density ",donate_density
              print *,"donate_slope ",donate_slope
              print *,"multi_volume_grid(im) ",multi_volume_grid(im)
              print *,"massdepart ",massdepart
              stop
             endif

             if (massdepart_mom.le.zero) then
              print *,"density (for mom) invalid in vfrac split 2"
              print *,"icrse,jcrse,kcrse ",icrse,jcrse,kcrse
              print *,"idonate,jdonate,kdonate ",idonate,jdonate,kdonate
              print *,"im ",im
              print *,"donate_density ",donate_density
              print *,"multi_volume_grid(im) ",multi_volume_grid(im)
              print *,"massdepart_mom ",massdepart_mom
              stop
             endif


              ! momslope is zero in partial volume fraction cells.
              ! In pure cells, momslope is (rho u)_x/rho_i
             do veldir=1,SDIM
              mom2(veldir)=zero 
              donate_data= &
               conserve(D_DECL(idonate,jdonate,kdonate),veldir)
               ! density has a piecewise constant reconstruction in mixed
               ! cells.   In pure cells, velocity slope is (rho u)_x/rho_i
               ! massdepart_mom is a density in the departure region.
              donate_slope= &
               momslope(D_DECL(idonate,jdonate,kdonate),veldir)
              mom2(veldir)=multi_volume_grid(im)*massdepart_mom* &
               (donate_data+donate_slope*moment_grid_diff(im))
             enddo  ! veldir=1..sdim (velocity)

              ! slope might not be zero.
             massdepart=massdepart*multi_volume_grid(im)
              ! slope is zero.
             massdepart_mom=massdepart_mom*multi_volume_grid(im)

             veldata(iden_mom_base+im)= &
              veldata(iden_mom_base+im)+massdepart_mom
             veldata(iden_base+dencomp_data)= &
              veldata(iden_base+dencomp_data)+massdepart

             ! skip density,then do energy,scalars,Q, ...
             ! for temperature:
             ! if incompressible: conserve=den * temp
             ! if compressible  : conserve=0.5 den |u|^2 + den * temp
             ! this is energy in the departure region.
             istate=2
             do while (istate.le.num_state_material)
              statecomp_data=(im-1)*num_state_material+istate

               ! conserve initialized in FORT_BUILD_CONSERVE.
               ! Temperature and species variables are multiplied by 
               ! dencore(im) in BUILD_CONSERVE.  (dencore(im) is the
               ! value of density stored in the state variable)
              donate_data= &
               conserve(D_DECL(idonate,jdonate,kdonate), &
                        iden_base+statecomp_data) 
              donate_slope= &
               momslope(D_DECL(idonate,jdonate,kdonate), &
                        iden_base+statecomp_data)
              
              veldata(iden_base+statecomp_data)= &
               veldata(iden_base+statecomp_data)+ &
               multi_volume_grid(im)*(donate_data+ & 
               donate_slope*moment_grid_diff(im))

              if (istate.eq.2) then
               if (veldata(iden_base+statecomp_data).lt.zero) then
                print *,"energy became negative "
                print *,"im,comp2 ",im,iden_base+statecomp_data
                print *,"current donated value ", &
                 veldata(iden_base+statecomp_data)
                print *,"icrse,jcrse,kcrse ",icrse,jcrse,kcrse 
                print *,"idonate,jdonate,kdonate ", &
                 idonate,jdonate,kdonate
                print *,"normdir,dir_counter ",normdir,dir_counter
                print *,"nmat ",nmat
                stop
               endif
              endif

              istate=istate+1
             enddo  !do while (istate.le.num_state_material)

             if ((num_materials_viscoelastic.ge.1).and. &
                 (num_materials_viscoelastic.le.nmat)) then

              if (fort_store_elastic_data(im).eq.1) then
               imap=1
               do while ((fort_im_elastic_map(imap)+1.ne.im).and. &
                         (imap.le.num_materials_viscoelastic))
                imap=imap+1
               enddo

               if (imap.le.num_materials_viscoelastic) then

                do istate=1,FORT_NUM_TENSOR_TYPE
                 statecomp_data=(imap-1)*FORT_NUM_TENSOR_TYPE+istate
                 ! configuration tensor is stored at the cell centers, not the
                 ! corresponding material centroid.
                 donate_data= &
                  tensor(D_DECL(idonate,jdonate,kdonate),statecomp_data)
                 veldata(itensor_base+statecomp_data)= &
                  veldata(itensor_base+statecomp_data)+ &
                  LS_voltotal_depart*donate_data 
                enddo !istate=1..FORT_NUM_TENSOR_TYPE

               else 
                print *,"imap invalid"
                stop
               endif
              else if (fort_store_elastic_data(im).eq.0) then
               ! do nothing
              else
               print *,"fort_store_elastic_data(im) invalid"
               stop
              endif

             else
              print *,"num_materials_viscoelastic invalid"
              stop
             endif

             ! level set function for im material.
             ! level set function is stored at the cell centers, not the
             ! corresponding material centroid.
             donate_data=LS(D_DECL(idonate,jdonate,kdonate),im) 
             veldata(iLS_base+im)=veldata(iLS_base+im)+ &
              LS_voltotal_depart*donate_data

             vofcomp=(im-1)*ngeom_raw+1
             ! material volume from departure (donating) region
             veldata(imof_base+vofcomp)= &
              veldata(imof_base+vofcomp)+multi_volume_grid(im)
             ! material volume from target (accepting) region
             veldata(iFtarget_base+im)= &
              veldata(iFtarget_base+im)+multi_volume(im)

             ! material centroid from target (accepting) region
             do dir2=1,SDIM
              veldata(imof_base+vofcomp+dir2)= &
               veldata(imof_base+vofcomp+dir2)+ &
               multi_volume(im)*multi_cen(dir2,im)
             enddo 

             do veldir=1,SDIM
               ! fluid materials tessellate the domain.
              if (is_rigid(nmat,im).eq.0) then
               veldata(veldir)=veldata(veldir)+mom2(veldir) 
              else if (is_rigid(nmat,im).eq.1) then
               ! do nothing
              else
               print *,"is_rigid invalid"
               stop
              endif
             enddo ! veldir=1..sdim
    
            enddo ! im=1,..,nmat (state variables, geometry, velocity)

             ! displacement is stored at the element centers, not the
             ! corresponding material centroids.
            if (FORT_NUM_TENSOR_TYPE*num_materials_viscoelastic+SDIM.eq. &
                NUM_CELL_ELASTIC) then
             if (num_MAC_vectors.eq.1) then
              do istate=1,SDIM
               statecomp_data= &
                 num_materials_viscoelastic*FORT_NUM_TENSOR_TYPE+istate
               donate_data= &
                tensor(D_DECL(idonate,jdonate,kdonate),statecomp_data)
               veldata(itensor_base+statecomp_data)= &
                veldata(itensor_base+statecomp_data)+ &
                LS_voltotal_depart*donate_data 
              enddo !istate=1..sdim
             else
              print *,"num_MAC_vectors invalid"
              stop
             endif
            else if (FORT_NUM_TENSOR_TYPE*num_materials_viscoelastic.eq. &
                     NUM_CELL_ELASTIC) then
             if (num_MAC_vectors.eq.2) then
              ! do nothing
             else
              print *,"num_MAC_vectors invalid"
              stop
             endif
            else
             print *,"NUM_CELL_ELASTIC invalid"
             stop
            endif

           endif ! volint>0

          else if (check_intersection.eq.0) then
           ! do nothing
          else
           print *,"check_intersection invalid"
           stop
          endif

         enddo  ! istencil=idonatelow..idonatehigh

         voltotal_depart=zero
         voltotal_target=zero
         do im=1,nmat
          vofcomp=(im-1)*ngeom_raw+1 
          volmat_target(im)=veldata(iFtarget_base+im)
          volmat_depart(im)=veldata(imof_base+vofcomp)

          volmat_target_cor(im)=volmat_target(im)
          volmat_depart_cor(im)=volmat_depart(im)

           ! fluid materials tessellate the domain.
          if (is_rigid(nmat,im).eq.0) then
           voltotal_target=voltotal_target+volmat_target(im)
           voltotal_depart=voltotal_depart+volmat_depart(im)
          else if (is_rigid(nmat,im).eq.1) then
           ! do nothing
          else
           print *,"is_rigid invalid"
           stop
          endif
         enddo ! im=1..nmat
        
         if (voltotal_depart.le.zero) then
          print *,"voltotal_depart bust "
          stop
         endif
         if (voltotal_target.le.zero) then
          print *,"voltotal_target bust "
          stop
         endif

          ! volDm/volD+(volDm/volT-volDm/volD)=volDm/volD+
          ! (volDm/volD)*(volD/volT-1) 
          ! dt div u = 1 - vol_depart/vol_target 
          !  or in 1D,
          ! vol_depart=xright - uright dt  - (xleft - uleft dt)=
          !            vol_target-div u dt dx
          ! 1-vol_depart/vol_target=1-(vol_target-div u dt dx)/vol_target=
          !                         dt div u

         divuterm=one-voltotal_depart/voltotal_target

         do im=1,nmat
          vofcomp=(im-1)*ngeom_raw+1 

          newvfrac(im)=volmat_target(im)/voltotal_target
          newvfrac_cor(im)=volmat_target_cor(im)/voltotal_target

          if (vofls0(D_DECL(icrse,jcrse,kcrse),im).le.half) then
           newvfrac_weymouth(im)=volmat_depart_cor(im)/voltotal_target
           if (newvfrac_weymouth(im).gt.one) then
            newvfrac_weymouth(im)=one
           endif
          else if (vofls0(D_DECL(icrse,jcrse,kcrse),im).ge.half) then
           newvfrac_weymouth(im)=one- &
            (voltotal_depart-volmat_depart_cor(im))/voltotal_target
           if (newvfrac_weymouth(im).lt.zero) then
            newvfrac_weymouth(im)=zero
           endif
          else
           print *,"vofls0 bust"
           stop
          endif

          if (newvfrac(im).le.VOFTOL) then
           newvfrac_weymouth(im)=newvfrac(im)
           newvfrac_cor(im)=newvfrac(im)
          endif
     
          nhalf=1 
          call CISBOX(xsten_accept,nhalf, &
           xlo,dx,icrse,jcrse,kcrse, &
           bfact,level, &
           volcell_accept,cencell_accept,SDIM)
  
          do dir2=1,SDIM
           if (newvfrac(im).gt.VOFTOL) then
            newcen(dir2,im)= &
             veldata(imof_base+vofcomp+dir2)/ &
             volmat_target(im)- &
             cencell_accept(dir2)
           else
            newcen(dir2,im)=zero
           endif
          enddo ! dir2

         enddo  ! im (geometry)

         call consistent_materials(newvfrac_cor,newcen,nmat)

         if ((EILE_flag.eq.1).or. & ! EILE
             (EILE_flag.eq.2).or. & ! EI
             (EILE_flag.eq.3)) then ! LE
          ! do nothing
         else if (EILE_flag.eq.-1) then ! weymouth and Yue
          do im=1,nmat
           newvfrac_cor(im)=newvfrac_weymouth(im)
          enddo
          call consistent_materials(newvfrac_cor,newcen,nmat)
         else
          print *,"EILE_flag invalid"
          stop
         endif
   
         ! pressure
         statecomp_data=SDIM+1
         snew_hold(statecomp_data)= &
           velfab(D_DECL(icrse,jcrse,kcrse),statecomp_data)

         ! density
         do im=1,nmat

          dencomp_data=(im-1)*num_state_material+1
          massdepart=veldata(iden_base+dencomp_data)
          if (massdepart.ge.zero) then
           ! do nothing
          else
           print *,"new mass cannot be negative"
           print *,"im= ",im
           print *,"new mass= ",massdepart
           stop
          endif
           ! if is_rigid==1 or voldepart<eps or voltarget<eps then
           !  den=fort_denconst(im)
           ! else if mat_type==0 then
           !  if override==0 or 2 then
           !   den=fort_denconst(im)
           !  else if override==1 then
           !   den=massdepart/voldepart
           !  endif
           ! else if mat_type>0 then
           !  den=massdepart/voltarget
           ! endif
          if (all_incomp.eq.1) then
           vol_target_local=volmat_depart_cor(im)
          else if (all_incomp.eq.0) then
           vol_target_local=volmat_target_cor(im)
          else
           print *,"all_incomp invalid"
           stop
          endif

          if (temperature_primitive_variable(im).eq.1) then
           vol_target_local=volmat_depart_cor(im)
          else if (temperature_primitive_variable(im).eq.0) then
           ! do nothing
          else
           print *,"temperature_primitive_variable(im) invalid"
           stop
          endif

          ! if is_rigid(im), density=fort_denconst(im)
          ! if incompressible,
          !   if constant_density_all_time==1 then density=fort_denconst(im)
          !   if constant_density_all_time==0 then 
          !                                  density=mass_depart/vol_depart
          ! if compressible,
          !  if constant_density_all_time==0 then
          !   density=massdepart/voltarget
          !  else
          !   return error.
          ! subroutine derive_density declared in GODUNOV_3D.F90 (this file)
          call derive_density(volmat_depart_cor(im), &
           vol_target_local,voltotal_depart, &
           override_density, &
           constant_density_all_time, &
           massdepart,im,nmat, &
           dencore(im))
          istate=dencomp+(im-1)*num_state_material+1
          if (dencore(im).gt.zero) then
           ! do nothing
          else
           print *,"density must be positive vfrac_split 2"
           print *,"im,dencore(im) ",im,dencore(im)
           stop
          endif
          if (dencore(im).lt.density_floor(im)) then
           dencore(im)=density_floor(im)
          endif
          if (density_ceiling(im).gt.zero) then
           if (dencore(im).gt.density_ceiling(im)) then
            dencore(im)=density_ceiling(im)
           endif
          else
           print *,"density_ceiling(im) invalid"
           stop
          endif
          snew_hold(istate)=dencore(im)

         enddo ! im, updating density

         ! levelset function
         ! voltotal_depart=sum_{fluid mat} volmat_depart(im)
         ! iLS_base=nmat*sdim+nmat*num_state_material+nmat*ngeom_raw
         do im=1,nmat
          if (voltotal_depart.gt.zero) then
           newLS(im)=veldata(iLS_base+im)/voltotal_depart
          else
           print *,"voltotal_depart invalid"
           stop
          endif 
         enddo  ! im=1..nmat (updating levelset vars)

         do im=1,nmat

          if ((num_materials_viscoelastic.ge.1).and. &
              (num_materials_viscoelastic.le.nmat)) then

           if (fort_store_elastic_data(im).eq.1) then
            imap=1
            do while ((fort_im_elastic_map(imap)+1.ne.im).and. &
                      (imap.le.num_materials_viscoelastic))
             imap=imap+1
            enddo

            if (imap.le.num_materials_viscoelastic) then

             do istate=1,FORT_NUM_TENSOR_TYPE
              statecomp_data=(imap-1)*FORT_NUM_TENSOR_TYPE+istate
              if (voltotal_depart.gt.zero) then
               tennew_hold(statecomp_data)= &
                 veldata(itensor_base+statecomp_data)/voltotal_depart
              else
               print *,"voltotal_depart invalid"
               stop
              endif 
 
             enddo !istate=1..FORT_NUM_TENSOR_TYPE

            else 
             print *,"imap invalid"
             stop
            endif
           else if (fort_store_elastic_data(im).eq.0) then
            ! do nothing
           else
            print *,"fort_store_elastic_data(im) invalid"
            stop
           endif

          else
           print *,"num_materials_viscoelastic invalid"
           stop
          endif
 
         enddo ! im=1..nmat (updating viscoelastic vars)

          ! displacement: (F_m X_m)_t + div(F_m u X_m)= u F_m
          !  (X_m)_t + u dot grad X_m = u
         if (FORT_NUM_TENSOR_TYPE*num_materials_viscoelastic+SDIM.eq. &
             NUM_CELL_ELASTIC) then
          if (num_MAC_vectors.eq.1) then
           do istate=1,SDIM
            statecomp_data= &
             num_materials_viscoelastic*FORT_NUM_TENSOR_TYPE+istate
            if (voltotal_depart.gt.zero) then
             tennew_hold(statecomp_data)= &
              veldata(itensor_base+statecomp_data)/voltotal_depart
            else
             print *,"voltotal_depart invalid"
             stop
            endif 
            if (istate.eq.normdir+1) then
             ! ucell is already the local displacement
             tennew_hold(statecomp_data)= &
              tennew_hold(statecomp_data)+ &
              ucell(D_DECL(icrse,jcrse,kcrse),istate)
            endif
           enddo !istate=1..SDIM
          else
           print *,"num_MAC_vectors invalid"
           stop
          endif

         else if (FORT_NUM_TENSOR_TYPE*num_materials_viscoelastic.eq. &
                  NUM_CELL_ELASTIC) then
          if (num_MAC_vectors.eq.2) then
           ! do nothing
          else
           print *,"num_MAC_vectors invalid"
           stop
          endif
         else
          print *,"NUM_CELL_ELASTIC invalid"
          stop
         endif

         ! velocity=mom/mass
         ! fluid materials tessellate the domain.
         totalmass_depart=zero
         do im=1,nmat
          if (is_rigid(nmat,im).eq.0) then
           massdepart_mom=veldata(iden_mom_base+im)
           totalmass_depart=totalmass_depart+massdepart_mom
          else if (is_rigid(nmat,im).eq.1) then
           ! do nothing
          else
           print *,"is_rigid invalid"
           stop
          endif
         enddo ! im=1..nmat

         if (totalmass_depart.le.zero) then
          print *,"totalmass_depart bust totalmass_depart=",totalmass_depart
          do dir2=1,SDIM
           print *,"dir,fablo,fabhi ",dir2,fablo(dir2),fabhi(dir2)
          enddo
          print *,"icrse,jcrse,kcrse ",icrse,jcrse,kcrse
          print *,"num_state_material=",num_state_material
          print *,"normdir=",normdir
          print *,"dir_counter=",dir_counter
          print *,"nmat,map_forward,level,finest_level ", &
           nmat,map_forward,level,finest_level
          stop
         endif

         do veldir=1,SDIM
          snew_hold(veldir)=veldata(veldir)/totalmass_depart
         enddo

          ! make sure 0<=F<=1 and sum F_i = 1.
          ! also truncation 1.0e-8 to 0 and 1-1.0e-8 to 1.
         call consistent_materials(newvfrac_cor,newcen,nmat)

         do im=1,nmat

          vofcomp=(im-1)*ngeom_raw+1

          KE=zero
          do veldir=1,SDIM
           momcomp=veldir
           vel1D=snew_hold(momcomp)
           KE=KE+vel1D**2
          enddo ! veldir
          KE=half*KE

          if (ngeom_raw.eq.SDIM+1) then
           snew_hold(mofcomp+vofcomp)=newvfrac_cor(im)
           do dir2=1,SDIM
            snew_hold(mofcomp+vofcomp+dir2)=newcen(dir2,im)
           enddo
          else
           print *,"ngeom_raw invalid in vfrac split"
           print *,"ngeom_raw= ",ngeom_raw
           stop
          endif

          no_material_flag=0

          if ( (volmat_depart(im).le. &
                VOFTOL*voltotal_depart).or. &
               (volmat_depart_cor(im).le. &
                VOFTOL*voltotal_depart).or. &
               (newvfrac_cor(im).le.VOFTOL).or. &
               (volmat_target(im).le. &
                VOFTOL*voltotal_depart).or. &
               (volmat_target_cor(im).le. &
                VOFTOL*voltotal_depart) ) then

           no_material_flag=1

          endif

           ! in: FORT_VFRAC_SPLIT
          dencomp_data=(im-1)*num_state_material+1

          istate=1
          do while (istate.le.num_state_material)

           if (istate.eq.1) then
            ! do nothing, density updated above
            istate=istate+1
           else if (istate.eq.2) then 

            do ispecies=1,num_species_var
             speccomp_data=(im-1)*num_state_material+num_state_base+ &
               ispecies
             if (no_material_flag.eq.1) then ! no material (im)
              snew_hold(dencomp+speccomp_data)=zero
             else if (no_material_flag.eq.0) then
              if (is_rigid(nmat,im).eq.1) then ! mass fraction=0 in solids.
               snew_hold(dencomp+speccomp_data)=zero
              else if (is_rigid(nmat,im).eq.0) then
               massdepart=veldata(iden_base+dencomp_data)
               if (massdepart.gt.zero) then
                snew_hold(dencomp+speccomp_data)= &
                 veldata(iden_base+speccomp_data)/massdepart
               else
                print *,"massdepart invalid"
                stop
               endif 
              else
               print *,"is_rigid invalid"
               stop
              endif
             else 
              print *,"no_material_flag invalid"
              stop
             endif

            enddo ! ispecies=1..num_species_var

            tempcomp_data=(im-1)*num_state_material+istate

            if (no_material_flag.eq.1) then
             snew_hold(dencomp+tempcomp_data)=fort_tempconst(im)
            else if (no_material_flag.eq.0) then
             if (is_rigid(nmat,im).eq.1) then
              if (fort_material_type(im).ne.999) then
               print *,"fort_material_type(im).ne.999"
               stop
              endif

              ! solidheat_flag==0 diffuse in solid
              ! solidheat_flag==1 dirichlet bc at solid-fluid
              ! solidheat_flag==2 insulating bc at solid-fluid
              if (solidheat_flag.eq.0) then ! diffuse in solid

               massdepart=veldata(iden_base+dencomp_data)
               if (massdepart.le.zero) then
                print *,"massdepart invalid"
                stop
               endif 
               ETcore=veldata(iden_base+tempcomp_data)/massdepart

              else if (solidheat_flag.eq.2) then ! neumann

               ! placeholder
               ETcore=fort_tempconst(im)

              else if (solidheat_flag.eq.1) then ! dirichlet

               ! placeholder
               ETcore=fort_tempconst(im)

              else
               print *,"solidheat_flag invalid"
               stop
              endif

              if (ETcore.lt.fort_tempcutoff(im)) then
               ETcore=fort_tempcutoff(im)
              endif
              if (ETcore.gt.fort_tempcutoffmax(im)) then
               ETcore=fort_tempcutoffmax(im)
              endif

              if (ETcore.le.zero) then
               print *,"Energy went negative"
               stop
              endif

              snew_hold(dencomp+tempcomp_data)=ETcore
             else if (is_rigid(nmat,im).eq.0) then
              if ((fort_material_type(im).ge.0).and. &
                  (fort_material_type(im).le.MAX_NUM_EOS)) then
               ! do nothing
              else
               print *,"fort_material_type invalid"
               stop
              endif
              massdepart=veldata(iden_base+dencomp_data)
              if (massdepart.le.zero) then
               print *,"massdepart invalid"
               stop
              endif 
              ! integral_omega_depart rho T F_m /
              ! integral_omega_depart rho F_m
              if (temperature_primitive_variable(im).eq.1) then
               ETcore=veldata(iden_base+tempcomp_data)/massdepart
              else if (temperature_primitive_variable(im).eq.0) then
               ! integral_omega_depart rho (u dot u/2 + c_v T) F_m /
               ! integral_omega_depart rho F_m
               ETcore=veldata(iden_base+tempcomp_data)/massdepart
               local_internal=ETcore-KE
               if (local_internal.gt.zero) then

                call init_massfrac_parm(dencore(im),massfrac_parm,im)
                do ispecies=1,num_species_var
                 speccomp_data=(im-1)*num_state_material+num_state_base+ &
                    ispecies
                  !dencomp=(SDIM+1)
                 massfrac_parm(ispecies)= &
                    snew_hold(dencomp+speccomp_data)
                 if (massfrac_parm(ispecies).ge.zero) then
                  ! do nothing
                 else
                  print *,"massfrac_parm(ispecies) invalid"
                  stop
                 endif
                enddo ! ispecies=1..num_species_var

                call TEMPERATURE_material(dencore(im),massfrac_parm, &
                 ETcore,local_internal,fort_material_type(im),im) 
               else
                ETcore=fort_tempcutoff(im)
               endif
              else
               print *,"temperature_primitive_variable invalid"
               stop
              endif
              if (ETcore.lt.fort_tempcutoff(im)) then
               ETcore=fort_tempcutoff(im)
              endif
              if (ETcore.gt.fort_tempcutoffmax(im)) then
               ETcore=fort_tempcutoffmax(im)
              endif
              if (ETcore.le.zero) then
               print *,"Energy went negative"
               stop
              endif
              snew_hold(dencomp+tempcomp_data)=ETcore
             else
              print *,"is_rigid invalid"
              stop
             endif
            else 
             print *,"no_material_flag invalid"
             stop
            endif

            istate=istate+1+num_species_var
           else
            print *,"istate invalid"
            stop
           endif

          enddo ! do while (istate.le.num_state_material)

         enddo  ! im=1..nmat


         do istate=1,(SDIM+1)
          snew(D_DECL(icrse,jcrse,kcrse),istate)=snew_hold(istate)
         enddo

          ! density, temperature, other scalars
          ! volume fractions, centroids
         do im=1,nmat
          if (is_rigid(nmat,im).eq.0) then
           do istate=1,num_state_material
            statecomp_data=dencomp+(im-1)*num_state_material+istate
            snew(D_DECL(icrse,jcrse,kcrse),statecomp_data)= &
              snew_hold(statecomp_data)
           enddo ! istate=1..num_state_material

           do igeom=1,ngeom_raw
            statecomp_data=dencomp+ &
              nmat*num_state_material+ &
              (im-1)*ngeom_raw+igeom
            snew(D_DECL(icrse,jcrse,kcrse),statecomp_data)= &
              snew_hold(statecomp_data)
           enddo
          else if (is_rigid(nmat,im).eq.1) then

           if (solidheat_flag.eq.0) then ! diffuse in solid
            tempcomp_data=dencomp+(im-1)*num_state_material+2
            snew(D_DECL(icrse,jcrse,kcrse),tempcomp_data)= &
              snew_hold(tempcomp_data)
           else if (solidheat_flag.eq.2) then ! neumann
            ! do nothing
           else if (solidheat_flag.eq.1) then ! dirichlet
            ! do nothing
           else
            print *,"solidheat_flag invalid"
            stop
           endif

          else
           print *,"is_rigid invalid"
           stop
          endif

          if ((num_materials_viscoelastic.ge.1).and. &
              (num_materials_viscoelastic.le.nmat)) then

           if (fort_store_elastic_data(im).eq.1) then
            imap=1
            do while ((fort_im_elastic_map(imap)+1.ne.im).and. &
                      (imap.le.num_materials_viscoelastic))
             imap=imap+1
            enddo
            if (imap.le.num_materials_viscoelastic) then

             do istate=1,FORT_NUM_TENSOR_TYPE
              statecomp_data=(imap-1)*FORT_NUM_TENSOR_TYPE+istate
              tennew(D_DECL(icrse,jcrse,kcrse),statecomp_data)= &
               tennew_hold(statecomp_data)
             enddo !istate=1..FORT_NUM_TENSOR_TYPE

            else 
             print *,"imap invalid"
             stop
            endif
           else if (fort_store_elastic_data(im).eq.0) then
            ! do nothing
           else
            print *,"fort_store_elastic_data(im) invalid"
            stop
           endif

          else
           print *,"num_materials_viscoelastic invalid"
           stop
          endif

         enddo ! im=1..nmat

         if (FORT_NUM_TENSOR_TYPE*num_materials_viscoelastic+SDIM.eq. &
             NUM_CELL_ELASTIC) then

          if (num_MAC_vectors.eq.1) then

           do istate=1,SDIM
            statecomp_data= &
             num_materials_viscoelastic*FORT_NUM_TENSOR_TYPE+istate
            tennew(D_DECL(icrse,jcrse,kcrse),statecomp_data)= &
             tennew_hold(statecomp_data)
           enddo !istate=1..SDIM

          else
           print *,"num_MAC_vectors invalid"
           stop
          endif

         else if (FORT_NUM_TENSOR_TYPE*num_materials_viscoelastic.eq. &
                  NUM_CELL_ELASTIC) then

          if (num_MAC_vectors.eq.2) then
           ! do nothing
          else
           print *,"num_MAC_vectors invalid"
           stop
          endif

         else
          print *,"NUM_CELL_ELASTIC invalid"
          stop
         endif

         do im=1,nmat
          if (is_rigid(nmat,im).eq.0) then
           LSnew(D_DECL(icrse,jcrse,kcrse),im)=newLS(im)

           ! level set comes from Lagrangian representation.
          else if (is_rigid(nmat,im).eq.1) then
           ! do nothing
          else
           print *,"is_rigid(nmat,im) invalid"
           stop
          endif
         enddo ! im=1..nmat

        else if (check_accept.eq.0) then
         ! do nothing
        else
         print *,"check_accept invalid"
         stop
        endif

       else if (maskcell.eq.0) then
        ! do nothing
       else
        print *,"maskcell invalid"
        stop
       endif

      enddo
      enddo
      enddo ! icrse,jcrse,kcrse -> growntilebox(0 ghost cells)

      if (DO_SANITY_CHECK.eq.1) then
       print *,"SANITY CHECK AFTER ADVECTION dir_counter= ",dir_counter
       print *,"SANITY CHECK AFTER ADVECTION normdir= ",normdir
       call CISL_sanity(dt,normdir,dir_counter,2)
       allocate(compareconserve(fablo(1)-1:fabhi(1)+1,5))
       allocate(comparestate(fablo(1)-1:fabhi(1)+1,5))
       jcrse=1
       kcrse=0
       do icrse=fablo(1),fabhi(1)
        compareconserve(icrse,1)= &
          conserve(D_DECL(icrse,jcrse,kcrse),1)* &
          conserve(D_DECL(icrse,jcrse,kcrse),iden_base+1)
        compareconserve(icrse,2)= &
          conserve(D_DECL(icrse,jcrse,kcrse),iden_base+1)
        compareconserve(icrse,3)= &
          conserve(D_DECL(icrse,jcrse,kcrse),iden_base+2)
        comparestate(icrse,1)=snew(D_DECL(icrse,jcrse,kcrse),1)
        comparestate(icrse,2)=snew(D_DECL(icrse,jcrse,kcrse),dencomp+1)
        comparestate(icrse,3)=snew(D_DECL(icrse,jcrse,kcrse),dencomp+2)
       enddo
       call compare_sanity(compareconserve,1,3,2)
       call compare_sanity(comparestate,1,3,1)
       deallocate(compareconserve)
       deallocate(comparestate)
       print *,"AFTER SANITY CHECK AFTER ADVECTION dir_counter= ",dir_counter
       print *,"AFTER SANITY CHECK AFTER ADVECTION normdir= ",normdir
      endif

      return
      end subroutine FORT_VFRAC_SPLIT


      ! combine_flag==0 (FVM -> GFM)
      ! combine_flag==1 (GFM -> FVM)
      ! combine_flag==2 (combine if vfrac<VOFTOL)
      ! project_option==3 (cell centered velocity)
      ! project_option==0 (MAC velocity - COMBINEVELFACE is called)
      subroutine FORT_COMBINEVEL( &
       tid, &
       hflag, &
       num_materials_combine, &
       mass_fraction_id, &
       latent_heat, &
       freezing_model, &
       distribute_from_target, &
       saturation_temp, &
       hydrate_flag, & ! scalar
       nmat, &
       nparts, &
       nparts_def, &
       im_solid_map, &
       nten, &
       nsolve, &
       project_option, &
       combine_idx, &
       combine_flag, &
       interface_cond_avail, &
       nstate_main, &
       ncomp_cell, &
       scomp, &
       ncomp, &
       scomp_size, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       finest_level, &
       ntsat, &
       TgammaFAB,DIMS(TgammaFAB), &
       maskcov,DIMS(maskcov), &
       solxfab,DIMS(solxfab), &
       solyfab,DIMS(solyfab), &
       solzfab,DIMS(solzfab), &
       LSNEW,DIMS(LSNEW), &
       LS,DIMS(LS), &
       vof,DIMS(vof), &
       cellfab,DIMS(cellfab), &
       newcell,DIMS(newcell), &
       state,DIMS(state), &  !Snew
       velbc, &
       listbc, &
       xlo,dx, &
       cur_time)
      use probf90_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      use godunov_module
      use mass_transfer_module
 
      IMPLICIT NONE

      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in) :: num_materials_combine
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nparts
      INTEGER_T, intent(in) :: nparts_def
      INTEGER_T, intent(in) :: im_solid_map(nparts_def)
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: hflag
      INTEGER_T, intent(in) :: mass_fraction_id(2*nten)
      REAL_T, intent(in) :: latent_heat(2*nten)
      INTEGER_T, intent(in) :: freezing_model(2*nten)
      INTEGER_T, intent(in) :: distribute_from_target(2*nten)
      REAL_T, intent(in) :: saturation_temp(2*nten)
      INTEGER_T, intent(in) :: hydrate_flag
      INTEGER_T, intent(in) :: nsolve
      INTEGER_T, intent(in) :: project_option
      INTEGER_T, intent(in) :: combine_idx
      INTEGER_T, intent(in) :: combine_flag
      INTEGER_T, intent(in) :: interface_cond_avail
      INTEGER_T, intent(in) :: nstate_main
      INTEGER_T, intent(in) :: ncomp_cell
      INTEGER_T, intent(in) :: scomp_size
      INTEGER_T, intent(in) :: scomp(scomp_size)
      INTEGER_T, intent(in) :: ncomp(scomp_size)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: ntsat

      INTEGER_T, intent(in) :: DIMDEC(TgammaFAB)
      INTEGER_T, intent(in) :: DIMDEC(maskcov)
      INTEGER_T, intent(in) :: DIMDEC(solxfab)
      INTEGER_T, intent(in) :: DIMDEC(solyfab)
      INTEGER_T, intent(in) :: DIMDEC(solzfab)
      INTEGER_T, intent(in) :: DIMDEC(LSNEW)
      INTEGER_T, intent(in) :: DIMDEC(LS)
      INTEGER_T, intent(in) :: DIMDEC(vof)
      INTEGER_T, intent(in) :: DIMDEC(cellfab)
      INTEGER_T, intent(in) :: DIMDEC(newcell)
      INTEGER_T, intent(in) :: DIMDEC(state)

      INTEGER_T, intent(in) :: velbc(SDIM,2,SDIM)
      INTEGER_T, intent(in) :: listbc(SDIM,2, &
             nsolve*num_materials_combine)

      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM) 
      REAL_T :: dxmaxLS

      REAL_T, intent(in) :: cur_time

      REAL_T, intent(in) :: TgammaFAB(DIMV(TgammaFAB),ntsat)
      REAL_T, intent(in) :: maskcov(DIMV(maskcov))
      REAL_T, intent(in) :: solxfab(DIMV(solxfab),nparts_def*SDIM)
      REAL_T, intent(in) :: solyfab(DIMV(solyfab),nparts_def*SDIM)
      REAL_T, intent(in) :: solzfab(DIMV(solzfab),nparts_def*SDIM)
      REAL_T, intent(in) :: LSNEW(DIMV(LSNEW),nmat*(1+SDIM))
      REAL_T, intent(in) :: LS(DIMV(LS),nmat*(1+SDIM))
      REAL_T, intent(in) :: vof(DIMV(vof),nmat*ngeom_recon)
      REAL_T, intent(inout) :: cellfab(DIMV(cellfab),ncomp_cell)
      REAL_T, intent(out) :: newcell(DIMV(newcell), &
              nsolve*num_materials_combine)
      REAL_T, intent(in) :: state(DIMV(state),nstate_main) !Snew

      REAL_T DATA_FLOOR
 
      INTEGER_T nten_test 

      INTEGER_T i,j,k
      INTEGER_T dir
      INTEGER_T side
      INTEGER_T i1,j1,k1
      INTEGER_T k1lo,k1hi
      INTEGER_T im,im_opp
      INTEGER_T im_crit
      INTEGER_T im_primary
      INTEGER_T im_secondary
      INTEGER_T ireverse,iten

      INTEGER_T im_source
      INTEGER_T im_source_master
      INTEGER_T im_dest
      INTEGER_T im_dest_master
      INTEGER_T tsat_flag

      INTEGER_T vofcomp
      INTEGER_T start_freezing
      INTEGER_T nhalf

      REAL_T xsten(-3:3,SDIM)
      REAL_T xsten_ofs(-3:3,SDIM)

      REAL_T total_vol_cell
      REAL_T fluid_vfrac_cell
      REAL_T mass_sum
      REAL_T weight_sum

      REAL_T local_volume
      REAL_T local_mass

      REAL_T volcell
      REAL_T cencell(SDIM)

      REAL_T solid_dist

      REAL_T cell_LS(nmat)
      REAL_T cell_vfrac(nmat)
      REAL_T cell_mfrac(nmat)

      INTEGER_T is_solid_cell ! =0 no solid  1<=is_solid_cell<=nmat+1 o.t.

      INTEGER_T dencomp
      INTEGER_T cellcomp

      REAL_T vsol(nsolve)
      REAL_T velsum(nsolve)
      REAL_T ucombine(nsolve)

      REAL_T test_density
      REAL_T test_temp
      REAL_T LS_source,LS_dest

      REAL_T LL
      INTEGER_T local_freezing_model
      INTEGER_T distribute_from_targ
      REAL_T TSAT_master
      REAL_T TDIFF
      REAL_T TDIFF_master
      REAL_T T_out(1)
      REAL_T Tcenter(nmat)
      REAL_T thermal_state(nmat)

      REAL_T xtarget(SDIM)
      REAL_T xI(SDIM)
      REAL_T nrm(SDIM)

      REAL_T mofdata(nmat*ngeom_recon)
      INTEGER_T nmax

      INTEGER_T mask_test

      REAL_T XC_sten(D_DECL(-1:1,-1:1,-1:1),SDIM)
      REAL_T VF_sten(D_DECL(-1:1,-1:1,-1:1))
      REAL_T LS_sten(D_DECL(-1:1,-1:1,-1:1))
      REAL_T T_sten(D_DECL(-1:1,-1:1,-1:1))

      INTEGER_T partid
      INTEGER_T partid_vel_plus
      INTEGER_T partid_vel
      INTEGER_T im_solid_vel_plus
      INTEGER_T im_solid_vel
      INTEGER_T im_solid_vel_max
      REAL_T LSCRIT_solid_plus
      REAL_T LSCRIT_solid
      REAL_T LSTEST
      INTEGER_T ncomp_per_tsat
      INTEGER_T Tgamma_STATUS
      INTEGER_T ispec
      REAL_T Tgamma
      REAL_T TorYgamma_BC
      INTEGER_T tsat_comp
      INTEGER_T local_tessellate
      REAL_T xclamped(SDIM)
      REAL_T LS_clamped
      REAL_T vel_clamped(SDIM)
      REAL_T temperature_clamped

      DATA_FLOOR=zero

      nhalf=3
      nmax=POLYGON_LIST_MAX ! in: COMBINEVEL

      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      if ((hflag.ne.0).and.(hflag.ne.1)) then
       print *,"hflag invalid1  hflag=",hflag
       stop
      endif 
      if (combine_idx.lt.-1) then
       print *,"combine_idx invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid75"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid combinevel nten, nten_test ",nten,nten_test
       stop
      endif

      ncomp_per_tsat=2
      if (ntsat.eq.nten*(ncomp_per_tsat+1)) then
       ! do nothing
      else
       print *,"nstat invalid"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid combine vel"
       stop
      endif
      if ((hydrate_flag.ne.0).and.(hydrate_flag.ne.1)) then
       print *,"hydrate_flag invalid"
       stop
      endif
      if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
       print *,"nsolve invalid"
       stop
      endif

      if (nstate_main.ne.(SDIM+1)+ &
          nmat*(num_state_material+ngeom_raw)+1) then
       print *,"nstate_main invalid"
       stop
      endif

      if ((ncomp_cell.ne.num_materials_combine).and. &
          (ncomp_cell.ne.nstate_main)) then
       print *,"ncomp_cell invalid"
       stop
      endif
   
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      if ((nparts.lt.0).or.(nparts.gt.nmat)) then
       print *,"nparts invalid FORT_COMBINEVEL"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.nmat)) then
       print *,"nparts_def invalid FORT_COMBINEVEL"
       stop
      endif

      if (combine_idx.eq.-1) then

       if (project_option.eq.0) then
        print *,"project_option==0 not allowed here"
        stop
       else if (project_option.eq.2) then ! thermal conduction
        if (scomp_size.ne.nmat) then
         print *,"scomp_size invalid"
         stop
        endif
        do im=1,nmat
         if (ncomp(im).ne.1) then
          print *,"ncomp invalid37"
          stop
         endif
         if (scomp(im).ne.(SDIM+1)+ &
             (im-1)*num_state_material+1) then
          print *,"scomp invalid"
          stop
         endif
        enddo ! im=1..nmat
       else if (project_option.eq.3) then  ! viscosity
        if (scomp(1).ne.0) then
         print *,"scomp invalid"
         stop
        endif
        if (ncomp(1).ne.SDIM) then
         print *,"ncomp invalid38"
         stop
        endif
        if (scomp_size.ne.1) then
         print *,"scomp_size invalid"
         stop
        endif
       else if ((project_option.ge.100).and. &
                (project_option.le.100+num_species_var-1)) then
        if (scomp_size.ne.nmat) then
         print *,"scomp_size invalid"
         stop
        endif
        do im=1,nmat
         if (ncomp(im).ne.1) then
          print *,"ncomp invalid39"
          stop
         endif
         if (scomp(im).ne.(SDIM+1)+ &
             (im-1)*num_state_material+2+project_option-100) then
          print *,"scomp invalid"
          stop
         endif
        enddo ! im=1..nmat
       else
        print *,"project_option invalid"
        stop
       endif

      else if (combine_idx.ge.0) then
       ! do nothing
      else
       print *,"combine_idx invalid"
       stop
      endif

      if ((interface_cond_avail.eq.0).or. &
          (interface_cond_avail.eq.1)) then
       ! do nothing
      else
       print *,"interface_cond_avail invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(maskcov),1,-1,1276)
      call checkbound(fablo,fabhi,DIMS(solxfab),0,0,1276)
      call checkbound(fablo,fabhi,DIMS(solyfab),0,1,1276)
      call checkbound(fablo,fabhi,DIMS(solzfab),0,SDIM-1,1276)
      call checkbound(fablo,fabhi,DIMS(LSNEW),1,-1,1276)

      call checkbound(fablo,fabhi,DIMS(LS),1,-1,1276)

      call checkbound(fablo,fabhi,DIMS(vof),1,-1,1276)
      call checkbound(fablo,fabhi,DIMS(cellfab),1,-1,1273)
      call checkbound(fablo,fabhi,DIMS(newcell),0,-1,1273)
      call checkbound(fablo,fabhi,DIMS(state),1,-1,1273)
      call checkbound(fablo,fabhi,DIMS(TgammaFAB),1,-1,234)

      call get_dxmaxLS(dx,bfact,dxmaxLS)

      if (SDIM.eq.2) then
       k1lo=0
       k1hi=0
      else if (SDIM.eq.3) then
       k1lo=-1
       k1hi=1
      else
       print *,"dimension bust"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi, &
       growlo,growhi,0)

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       mask_test=NINT(maskcov(D_DECL(i,j,k)))

       !mask=tag if not covered by level+1 or outside the domain.
       if (mask_test.eq.1) then

        call gridsten_level(xsten,i,j,k,level,nhalf)
        do dir=1,SDIM
         xclamped(dir)=xsten(0,dir)
        enddo

        do im=1,nmat*ngeom_recon
         mofdata(im)=vof(D_DECL(i,j,k),im)
        enddo
        local_tessellate=3
        call multi_get_volume_tessellate( &
         local_tessellate, &  ! =3
         bfact, &
         dx,xsten,nhalf, &
         mofdata, &
         geom_xtetlist(1,1,1,tid+1), &
         nmax, &
         nmax, &
         nmat,SDIM,3)

        mass_sum=zero
        fluid_vfrac_cell=zero
        total_vol_cell=zero

        do im=1,nmat

         vofcomp=(im-1)*ngeom_recon+1
         local_volume=mofdata(vofcomp)
         
         if ((local_volume.ge.-VOFTOL).and. &
             (local_volume.le.one+VOFTOL)) then
          if (local_volume.lt.zero) then
           local_volume=zero
          endif
          if (local_volume.gt.one) then
           local_volume=one
          endif
         else
          print *,"local_volume invalid"
          stop
         endif
         cell_vfrac(im)=local_volume

         dencomp=(SDIM+1)+ &
             (im-1)*num_state_material+1
         test_density=state(D_DECL(i,j,k),dencomp)
         if (test_density.le.zero) then
          print *,"test_density invalid"
          stop
         endif
         local_mass=test_density*local_volume ! local_volume is a volume frac.
         total_vol_cell=total_vol_cell+local_volume
         if (is_rigid(nmat,im).eq.0) then
          fluid_vfrac_cell=fluid_vfrac_cell+local_volume
         else if (is_rigid(nmat,im).eq.1) then
          ! do nothing
         else
          print *,"is_rigid invalid"
          stop
         endif
         mass_sum=mass_sum+local_mass

         cell_mfrac(im)=local_mass

        enddo ! im=1,nmat

        if (total_vol_cell.le.zero) then
         print *,"total_vol_cell invalid: ",total_vol_cell
         stop
        endif

        if (mass_sum.le.zero) then
         print *,"mass_sum invalid"
         stop
        endif

        fluid_vfrac_cell=fluid_vfrac_cell/total_vol_cell

        if ((fluid_vfrac_cell.ge.zero).and. &
            (fluid_vfrac_cell.le.one+VOFTOL)) then
         ! do nothing
        else
         print *,"fluid_vfrac_cell invalid"
         stop
        endif

        do im=1,nmat
         cell_vfrac(im)=cell_vfrac(im)/total_vol_cell
         cell_LS(im)=cell_vfrac(im)-half
        enddo ! im

        partid=0
        partid_vel_plus=0
        partid_vel=0
        im_solid_vel_plus=0
        im_solid_vel=0
        LSCRIT_solid_plus=-1.0D+99
        LSCRIT_solid=-1.0D+99
        
        do im=1,nmat
         if (is_lag_part(nmat,im).eq.1) then
          if (is_rigid(nmat,im).eq.1) then
           LSTEST=LSNEW(D_DECL(i,j,k),im)
           if (LSTEST.ge.zero) then
            if (im_solid_vel_plus.eq.0) then
             partid_vel_plus=partid
             im_solid_vel_plus=im
             LSCRIT_solid_plus=LSTEST
            else if ((im_solid_vel_plus.ge.1).and. &
                     (im_solid_vel_plus.le.nmat)) then
             if (LSTEST.ge.LSCRIT_solid_plus) then
              partid_vel_plus=partid
              im_solid_vel_plus=im
              LSCRIT_solid_plus=LSTEST
             endif
            else
             print *,"im_solid_vel_plus invalid combinevel:",im_solid_vel_plus
             stop
            endif
           else if (LSTEST.lt.zero) then
            ! do nothing
           else
            print *,"LSTEST invalid"
            stop
           endif

           if (im_solid_vel.eq.0) then
            partid_vel=partid
            im_solid_vel=im
            LSCRIT_solid=LSTEST
           else if ((im_solid_vel.ge.1).and. &
                    (im_solid_vel.le.nmat)) then
            if (LSTEST.ge.LSCRIT_solid) then
             partid_vel=partid
             im_solid_vel=im
             LSCRIT_solid=LSTEST
            endif
           else
            print *,"im_solid_vel invalid"
            stop
           endif

          else if (is_rigid(nmat,im).eq.0) then
           ! do nothing
          else
           print *,"is_rigid(nmat,im) invalid"
           stop
          endif
          partid=partid+1
         else if (is_lag_part(nmat,im).eq.0) then
          if (is_rigid(nmat,im).eq.0) then
           ! do nothing
          else
           print *,"is_rigid(nmat,im) invalid"
           stop
          endif
         else
          print *,"is_lag_part(nmat,im) invalid"
          stop
         endif
        enddo ! im=1..nmat

        if (partid.ne.nparts) then
         print *,"partid invalid"
         stop
        endif

        is_solid_cell=0

        if ((im_solid_vel_plus.ge.1).and. &
            (im_solid_vel_plus.le.nmat)) then
         if (is_prescribed(nmat,im_solid_vel_plus).eq.1) then
          is_solid_cell=im_solid_vel_plus
          if (im_solid_map(partid_vel_plus+1)+1.ne.im_solid_vel_plus) then
           print *,"im_solid_map(partid_vel_plus+1)+1.ne.im_solid_vel_plus"
           stop
          endif
         else if (is_prescribed(nmat,im_solid_vel_plus).eq.0) then
          ! do nothing
         else
          print *,"is_prescribed(nmat,im_solid_vel_plus) invalid"
          stop
         endif
        else if (im_solid_vel_plus.eq.0) then
         ! do nothing
        else
         print *,"im_solid_vel_plus invalid combinevel:2",im_solid_vel_plus
         stop
        endif

         ! solid_dist>0 in the solid
         ! returns: im_solid_vel_max=max_{is_rigid(im)==1} cell_LS(im)
        call combine_solid_LS(cell_LS,nmat,solid_dist,im_solid_vel_max)

        if (solid_dist.ge.zero) then
         if ((im_solid_vel_max.ge.1).and. &
             (im_solid_vel_max.le.nmat)) then
          if (is_prescribed(nmat,im_solid_vel_max).eq.1) then
           is_solid_cell=im_solid_vel_max
          else if (is_prescribed(nmat,im_solid_vel_max).eq.0) then
           ! do nothing
          else
           print *,"is_prescribed(im_solid_vel_max) invalid"
           stop
          endif
         else
          print *,"im_solid_vel_max invalid"
          stop
         endif
        else if (solid_dist.lt.zero) then
         ! do nothing
        else
         print *,"solid_dist invalid"
         stop
        endif

         ! first checks the rigid materials for a positive LS; if none
         ! exist, then the subroutine checks the fluid materials.
        call get_primary_material(cell_LS,nmat,im_primary)

        if (is_prescribed(nmat,im_primary).eq.1) then
         is_solid_cell=im_primary
        else if (is_prescribed(nmat,im_primary).eq.0) then
         ! do nothing
        else
         print *,"is_prescribed invalid"
         stop
        endif

        if (project_option.eq.3) then ! viscosity

         if (combine_flag.eq.2) then !combine if vfrac<VOFTOL

          if (nsolve.ne.SDIM) then
           print *,"nsolve invalid"
           stop
          endif
          if (scomp_size.ne.1) then
           print *,"scomp_size invalid"
           stop
          endif
          if (scomp(1).ne.0) then
           print *,"scomp invalid"
           stop
          endif
          if (ncomp(1).ne.nsolve) then
           print *,"ncomp(1) invalid"
           stop
          endif

          do cellcomp=1,SDIM
           if (hflag.eq.0) then
            if (im_solid_vel.eq.0) then
             vsol(cellcomp)=zero
            else if ((im_solid_vel.ge.1).and. &
                     (im_solid_vel.le.nmat)) then
             if (im_solid_map(partid_vel+1)+1.ne.im_solid_vel) then
              print *,"im_solid_map(partid_vel+1)+1.ne.im_solid_vel"
              stop
             endif
             if (cellcomp.eq.1) then
              vsol(cellcomp)=half* &
                (solxfab(D_DECL(i,j,k),partid_vel*SDIM+cellcomp)+ &
                 solxfab(D_DECL(i+1,j,k),partid_vel*SDIM+cellcomp))
             else if (cellcomp.eq.2) then
              vsol(cellcomp)=half* &
                (solyfab(D_DECL(i,j,k),partid_vel*SDIM+cellcomp)+ &
                 solyfab(D_DECL(i,j+1,k),partid_vel*SDIM+cellcomp))
             else if ((cellcomp.eq.SDIM).and.(SDIM.eq.3)) then
              vsol(cellcomp)=half* &
                (solzfab(D_DECL(i,j,k),partid_vel*SDIM+cellcomp)+ &
                 solzfab(D_DECL(i,j,k+1),partid_vel*SDIM+cellcomp))
             else
              print *,"cellcomp invalid"
              stop
             endif

            else
             print *,"im_solid_vel invalid"
             stop
            endif
           else if (hflag.eq.1) then
            vsol(cellcomp)=zero
           else
            print *,"hflag invalid2 hflag=",hflag
            stop
           endif
          enddo ! cellcomp=1..sdim
          
           ! LS>0 if clamped
          call SUB_clamped_LS(xclamped,cur_time,LS_clamped, &
                vel_clamped,temperature_clamped)

          if (LS_clamped.ge.zero) then
           is_solid_cell=nmat+1
           if (hflag.eq.0) then ! inhomogeneous
            do dir=1,SDIM
             vsol(dir)=vel_clamped(dir)
            enddo
           else if (hflag.eq.1) then
            do dir=1,SDIM
             vsol(dir)=zero
            enddo
           else
            print *,"hflag invalid3 hflag=",hflag
            stop
           endif
          else if (LS_clamped.lt.zero) then
           ! do nothing
          else
           print *,"LS_clamped is NaN"
           stop
          endif

          if ((is_solid_cell.ge.1).and. &
              (is_solid_cell.le.nmat+1)) then
   
           do cellcomp=1,SDIM
            ucombine(cellcomp)=vsol(cellcomp)
           enddo ! cellcomp

          else if (is_solid_cell.eq.0) then

           do cellcomp=1,SDIM
            ucombine(cellcomp)=cellfab(D_DECL(i,j,k),cellcomp)
           enddo ! cellcomp

          else
           print *,"is_solid_cell invalid: ",is_solid_cell
           stop
          endif

          do cellcomp=1,SDIM
           cellfab(D_DECL(i,j,k),cellcomp)=ucombine(cellcomp)
          enddo ! cellcomp

         else
          print *,"combine_flag invalid"
          stop
         endif
  
        else if ((project_option.eq.2).or. &     ! temperature
                 ((project_option.ge.100).and. & ! species
                  (project_option.le.100+num_species_var-1))) then

         do im=1,nmat

          if (cell_vfrac(im).le.one-VOFTOL) then

           vofcomp=(im-1)*ngeom_recon+1
          
           if (((im_primary.eq.im).and. &
               (combine_flag.eq.0)).or. &  ! FVM -> GFM
              ((cell_vfrac(im).ge.VOFTOL).and. &
               (combine_flag.eq.1))) then  ! GFM -> FVM

            if (combine_idx.ne.-1) then
             print *,"combine_idx invalid"
             stop
            endif
            if (num_materials_combine.ne.nmat) then
             print *,"num_materials_combine invalid"
             stop
            endif

            do i1=-1,1
            do j1=-1,1
            do k1=k1lo,k1hi
             call gridsten_level(xsten_ofs,i+i1,j+j1,k+k1,level,nhalf)
             call Box_volumeFAST(bfact,dx,xsten_ofs,nhalf, &
              volcell,cencell,SDIM)
             do dir=1,SDIM
              XC_sten(D_DECL(i1,j1,k1),dir)= &
               vof(D_DECL(i+i1,j+j1,k+k1),vofcomp+dir)+cencell(dir)
             enddo
             VF_sten(D_DECL(i1,j1,k1))= &
              vof(D_DECL(i+i1,j+j1,k+k1),vofcomp)
             LS_sten(D_DECL(i1,j1,k1))= &
              LS(D_DECL(i+i1,j+j1,k+k1),im)
            enddo
            enddo
            enddo ! i1,j1,k1

            if (combine_flag.eq.0) then ! centroid -> center
             do dir=1,SDIM
              xtarget(dir)=xsten(0,dir)
             enddo
            else if (combine_flag.eq.1) then ! center -> centroid
             do dir=1,SDIM
              xtarget(dir)=XC_sten(D_DECL(0,0,0),dir)
             enddo
            else
             print *,"combine_flag invalid"
             stop
            endif

            do dir=1,SDIM
             xI(dir)=xtarget(dir)
            enddo

            ! check for TSAT BC.
            tsat_flag=0
            im_source_master=0
            im_dest_master=0
            TSAT_master=273.0d0

             ! check for Tgamma or Ygamma boundary condition.

            do im_crit=1,nmat
             dencomp=(SDIM+1)+ &
              (im_crit-1)*num_state_material+1
             Tcenter(im_crit)=cellfab(D_DECL(i,j,k),scomp(im_crit)+1)
             if (Tcenter(im_crit).ge.zero) then
              ! do nothing
             else
              print *,"Tcenter(im_crit) invalid"
              stop
             endif
             if (project_option.eq.2) then
              thermal_state(im_crit)=Tcenter(im_crit)
             else if ((project_option.ge.100).and. & ! species
                      (project_option.le.100+num_species_var-1)) then
              thermal_state(im_crit)= &
                 state(D_DECL(i,j,k),dencomp+1)
             else
              print *,"project_option invalid"
              stop
             endif

            enddo ! im_crit=1..nmat

            do ireverse=0,1
             do im_opp=1,nmat
              if (im_opp.ne.im) then

               call get_iten(im,im_opp,iten,nmat)
               LL=latent_heat(iten+ireverse*nten)

               if (interface_cond_avail.eq.1) then

                if (LL.ne.zero) then

                 if (((ireverse.eq.0).and.(im.lt.im_opp)).or. &
                     ((ireverse.eq.1).and.(im.gt.im_opp))) then
                  im_source=im
                  im_dest=im_opp
                 else if (((ireverse.eq.0).and.(im.gt.im_opp)).or. &
                          ((ireverse.eq.1).and.(im.lt.im_opp))) then
                  im_source=im_opp
                  im_dest=im
                 else
                  print *,"ireverse invalid"
                  stop
                 endif

                 call check_recalesce_status(im_source,start_freezing)

                 if (start_freezing.eq.1) then

                  local_freezing_model=freezing_model(iten+ireverse*nten)
                  distribute_from_targ= &
                        distribute_from_target(iten+ireverse*nten)
                  if ((distribute_from_targ.lt.0).or. &
                      (distribute_from_targ.gt.1)) then
                   print *,"distribute_from_targ invalid"
                   stop
                  endif

                  if ((local_freezing_model.eq.0).or. &
                      (local_freezing_model.eq.5).or. &
                      (local_freezing_model.eq.6)) then ! Palmore/Desjardins

                   if ((im_primary.eq.im).or.(im_primary.eq.im_opp)) then

                    call get_secondary_material(cell_LS,nmat, &
                      im_primary,im_secondary)

                    if (im_primary.eq.im_secondary) then
                     print *,"cannot have im_primary.eq.im_secondary"
                     stop
                    endif

                    if ((im_secondary.eq.im).or. &
                        (im_secondary.eq.im_opp)) then

                     if ((cell_vfrac(im).ge.VOFTOL).and. &
                         (cell_vfrac(im_opp).ge.VOFTOL)) then

                      Tgamma_STATUS=NINT(TgammaFAB(D_DECL(i,j,k),iten))
                      if (ireverse.eq.0) then
                       ! do nothing
                      else if (ireverse.eq.1) then
                       Tgamma_STATUS=-Tgamma_STATUS
                      else
                       print *,"ireverse invalid"
                       stop
                      endif

                      if (project_option.eq.2) then
                       ! do nothing
                      else if ((project_option.ge.100).and. &
                               (project_option.lt.100+num_species_var)) then
                       if ((local_freezing_model.eq.0).or. & !saturated
                           (local_freezing_model.eq.5)) then !saturated
                        Tgamma_STATUS=0

                        ! Palmore/Desjardins, partial mass fraction
                       else if (local_freezing_model.eq.6) then 
                        ispec=mass_fraction_id(iten+ireverse*nten)
                        if (ispec.eq.project_option-100+1) then
                         ! do nothing
                        else if ((ispec.ge.1).and. &
                                 (ispec.le.num_species_var)) then
                         Tgamma_STATUS=0
                        else
                         print *,"ispec invalid"
                         stop
                        endif
                       else
                        print *,"local_freezing_model invalid"
                        stop
                       endif
                      else
                       print *,"project_option invalid"
                       stop
                      endif

                      if ((Tgamma_STATUS.eq.1).or.(Tgamma_STATUS.eq.2)) then

                       if (project_option.eq.2) then
                        ! default Tgamma
                        Tgamma=saturation_temp(iten+ireverse*nten)
                        TorYgamma_BC=Tgamma
                        if (Tgamma.gt.zero) then
                         tsat_comp=nten+(iten-1)*ncomp_per_tsat+1
                         Tgamma=TgammaFAB(D_DECL(i,j,k),tsat_comp)
                         TorYgamma_BC=Tgamma
                         if (Tgamma.gt.zero) then
                          ! do nothing
                         else
                          print *,"Tgamma must be positive1"
                          stop
                         endif
                        else
                         print *,"saturation temperature must be positive2"
                         stop
                        endif
                       else if ((project_option.ge.100).and. &
                                (project_option.lt.100+num_species_var)) then
                        Tgamma=saturation_temp(iten+ireverse*nten)
                        TorYgamma_BC=one
                        if (Tgamma.gt.zero) then
                         tsat_comp=nten+(iten-1)*ncomp_per_tsat+1
                         Tgamma=TgammaFAB(D_DECL(i,j,k),tsat_comp)
                         tsat_comp=nten+(iten-1)*ncomp_per_tsat+2
                         TorYgamma_BC=TgammaFAB(D_DECL(i,j,k),tsat_comp)
                         if (Tgamma.gt.zero) then
                          ! do nothing
                         else
                          print *,"Tgamma must be positive22"
                          stop
                         endif
                         if ((TorYgamma_BC.ge.zero).and. &
                             (TorYgamma_BC.le.one)) then
                          ! do nothing
                         else
                          print *,"TorYgamma_BC (aka Y) must be >= 0 and <=1"
                          stop
                         endif
                        else
                         print *,"saturation temperature must be positive33"
                         stop
                        endif
                       else
                        print *,"project_option invalid"
                        stop
                       endif

                       if (LL.lt.zero) then ! freezing
                        TDIFF=max(Tgamma-thermal_state(im), &
                                  Tgamma-thermal_state(im_opp))
                       else if (LL.gt.zero) then ! melting
                        TDIFF=max(thermal_state(im)-Tgamma, &
                                  thermal_state(im_opp)-Tgamma)
                       else
                        print *,"LL invalid"
                        stop
                       endif
                       if (tsat_flag.eq.0) then
                        tsat_flag=1
                        TSAT_master=TorYgamma_BC
                        TDIFF_master=TDIFF
                        im_source_master=im_source
                        im_dest_master=im_dest
                       else if (tsat_flag.eq.1) then
                        if (TDIFF.gt.TDIFF_master) then
                         TSAT_master=TorYgamma_BC
                         TDIFF_master=TDIFF
                         im_source_master=im_source
                         im_dest_master=im_dest
                        endif
                       else
                        print *,"tsat_flag invalid"
                        stop
                       endif 

                      else if (Tgamma_STATUS.eq.0) then
                       ! do nothing
                      else
                       print *,"Tgamma_STATUS invalid"
                       stop
                      endif

                     else if ((abs(cell_vfrac(im)).le.VOFTOL).or. &
                              (abs(cell_vfrac(im_opp)).le.VOFTOL)) then
                      ! do nothing
                     else
                      print *,"cell_vfrac invalid"
                      stop
                     endif
                    endif ! im_secondary==im or im_opp
                   endif ! im_primary=im or im_opp

                  else if ((local_freezing_model.eq.1).or. &
                           (local_freezing_model.eq.2).or. &
                           (local_freezing_model.eq.4).or. & !Tanasawa/Schrage
                           (local_freezing_model.eq.7)) then ! cavitation
                   ! do nothing
                  else 
                   print *,"local_freezing_model not supported"
                   stop
                  endif

                 else if (start_freezing.eq.0) then
                  ! do nothing
                 else
                  print *,"start_freezing invalid"
                  stop
                 endif

                else if (LL.eq.zero) then
                 ! do nothing
                else
                 print *,"LL invalid"
                 stop
                endif

               else if (interface_cond_avail.eq.0) then
                ! do nothing
               else
                print *,"interface_cond_avail invalid"
                stop
               endif

              else if (im_opp.eq.im) then
               ! do nothing
              else
               print *,"im_opp invalid"
               stop
              endif
         
             enddo ! im_opp
            enddo ! ireverse

             ! combine_flag==0 or 1.
            do i1=-1,1
            do j1=-1,1
            do k1=k1lo,k1hi
             test_temp=cellfab(D_DECL(i+i1,j+j1,k+k1),scomp(im)+1)

             if ((project_option.eq.2).or. & ! thermal combine
                 ((project_option.ge.100).and. &
                  (project_option.le.100+num_species_var-1))) then

              if (hflag.eq.0) then

               if (test_temp.ge.zero) then
                ! do nothing
               else
                print *,"project_option ",project_option
                print *,"test_temp invalid test_temp=",test_temp
                print *,"combine_flag=",combine_flag
                print *,"level,finest_level ",level,finest_level
                print *,"i,j,k ",i,j,k
                print *,"i1,j1,k1 ",i1,j1,k1
                do dir=1,SDIM
                 print *,"dir,growlo,growhi ",dir,growlo(dir),growhi(dir)
                enddo
                print *,"x= ",xsten(2*i1,1)
                print *,"y= ",xsten(2*j1,2)
                print *,"im= ",im
                print *,"INT_DIR= ",INT_DIR
                print *,"EXT_DIR= ",EXT_DIR
                print *,"REFLECT_EVEN= ",REFLECT_EVEN
                print *,"REFLECT_ODD= ",REFLECT_ODD
                print *,"FOEXTRAP= ",FOEXTRAP
                print *,"scomp_size=",scomp_size
                do im_opp=1,nmat
                 print *,"im_opp,F ",im_opp,cell_vfrac(im_opp)
                 print *,"im_opp,scomp,ncomp ", &
                  im_opp,scomp(im_opp),ncomp(im_opp)
                 do dir=1,SDIM
                  do side=1,2
                   print *,"dir,side,im_opp,listbc ",dir,side,im_opp, &
                     listbc(dir,side,im_opp)
                  enddo
                 enddo
                enddo ! im_opp
                 
                stop
               endif

              else if (hflag.eq.1) then
               ! do nothing
              else
               print *,"hflag invalid"
               stop
              endif

             else
              print *,"project_option invalid"
              stop
             endif
        
             T_sten(D_DECL(i1,j1,k1))=test_temp
            enddo
            enddo
            enddo ! i1,j1,k1

            if (nsolve.ne.1) then
             print *,"nsolve invalid"
             stop
            endif

            if (tsat_flag.eq.1) then

             if (is_rigid(nmat,im_primary).eq.0) then

              im_source=im_source_master
              im_dest=im_dest_master
              LS_source=LS(D_DECL(i,j,k),im_source)
              LS_dest=LS(D_DECL(i,j,k),im_dest)

              if ((abs(LS_source).le.two*dxmaxLS).and. &
                  (abs(LS_dest).le.two*dxmaxLS)) then
               if (LS_dest.ge.zero) then
                do dir=1,SDIM
                 nrm(dir)=LS(D_DECL(i,j,k),nmat+(im_source-1)*SDIM+dir)
                 xI(dir)=xsten(0,dir)-LS_source*nrm(dir)
                enddo
               else if (LS_source.ge.zero) then
                do dir=1,SDIM
                 nrm(dir)=LS(D_DECL(i,j,k),nmat+(im_dest-1)*SDIM+dir)
                 xI(dir)=xsten(0,dir)-LS_dest*nrm(dir)
                enddo
               else if ((LS_dest.le.zero).and. &
                        (LS_source.le.zero)) then
                tsat_flag=0
               else
                print *,"LS_dest or LS_source invalid"
                stop
               endif
              else if ((abs(LS_source).ge.two*dxmaxLS).or. &
                       (abs(LS_dest).ge.two*dxmaxLS)) then
               tsat_flag=0
              else
               print *,"LS_dest or LS_source invalid"
               stop
              endif

             else if (is_rigid(nmat,im_primary).eq.1) then
              tsat_flag=0
             else
              print *,"im_primary invalid"
              stop
             endif

            else if (tsat_flag.eq.0) then
             ! do nothing
            else
             print *,"tsat_flag invalid"
             stop
            endif

             ! in: FORT_COMBINEVEL
            call center_centroid_interchange( &
             DATA_FLOOR, &
             nsolve, &
             combine_flag,  & ! 0=>centroid -> center   1=>center->centroid
             tsat_flag, &
             bfact, &
             level, &
             finest_level, &
             dx,xlo, &
             xsten,nhalf, &
             T_sten, &  
             XC_sten, &  
             xI, &
             xtarget, &
             VF_sten, &
             LS_sten, &
             TSAT_master, &
             T_out)

            if (combine_flag.eq.0) then ! centroid -> center
             do im_opp=1,nmat
              newcell(D_DECL(i,j,k),im_opp)=T_out(1)
             enddo
            else if (combine_flag.eq.1) then  ! center -> centroid
             newcell(D_DECL(i,j,k),im)=T_out(1)
            else
             print *,"combine_flag invalid"
             stop
            endif

           else if ((im_primary.ne.im).and. &
                    (combine_flag.eq.0)) then  ! centroid -> center
            ! do nothing
           else if ((cell_vfrac(im).le.VOFTOL).and. &
                    ((combine_flag.eq.1).or. &   ! GFM->FVM
                     (combine_flag.eq.2))) then ! combine if F==0

            if (nsolve.ne.1) then
             print *,"nsolve invalid"
             stop
            endif
            if (num_materials_combine.ne.nmat) then
             print *,"num_materials_combine invalid"
             stop
            endif

            velsum(1)=zero
            weight_sum=zero

            do im_crit=1,nmat

             weight_sum=weight_sum+cell_mfrac(im_crit)

             if (combine_idx.eq.-1) then
              cellcomp=scomp(im_crit)+1
             else if (combine_idx.ge.0) then
              cellcomp=im_crit
             else
              print *,"combine_idx invalid"
              stop
             endif

             test_temp=cellfab(D_DECL(i,j,k),cellcomp)
             if (hflag.eq.0) then
              if (test_temp.ge.zero) then
               ! do nothing
              else
               print *,"test_temp must be positive: combinevel"
               print *,"test_temp=",test_temp
               print *,"im_crit=",im_crit
               print *,"cellcomp=",cellcomp
               stop
              endif
             else if (hflag.eq.1) then
              ! do nothing
             else
              print *,"hflag invalid3 hflag=",hflag
              stop
             endif

             velsum(1)=velsum(1)+cell_mfrac(im_crit)*test_temp

            enddo ! im_crit

            if (weight_sum.gt.zero) then
             velsum(1)=velsum(1)/weight_sum
            else
             print *,"weight_sum invalid 1"
             stop
            endif 

            if (hflag.eq.0) then
             if (velsum(1).ge.zero) then
              ! do nothing
             else
              print *,"velsum must be nonneg: combinevel"
              stop
             endif
            else if (hflag.eq.1) then
             ! do nothing
            else
             print *,"hflag invalid4 hflag=",hflag
             stop
            endif

            if (combine_idx.eq.-1) then
             cellcomp=scomp(im)+1
            else if (combine_idx.ge.0) then
             cellcomp=im
            else
             print *,"combine_idx invalid"
             stop
            endif

            if (combine_flag.eq.1) then
             newcell(D_DECL(i,j,k),im)=velsum(1)
            else if (combine_flag.eq.2) then
             cellfab(D_DECL(i,j,k),cellcomp)=velsum(1)
            else
             print *,"combine_flag invalid"
             stop
            endif

           else if ((cell_vfrac(im).ge.VOFTOL).and. &
                    (combine_flag.eq.2)) then
            ! do nothing
           else
            print *,"cell_vfrac or combine_flag bust"
            stop
           endif
   
          else if (cell_vfrac(im).ge.one-VOFTOL) then

           if (combine_flag.eq.0) then ! centroid -> center
            T_out(1)=cellfab(D_DECL(i,j,k),scomp(im)+1)
             ! since cell_vfrac(im)==1 => cell_vfrac(im_opp)==0.0 (im_opp!=im)
            do im_opp=1,nmat
             newcell(D_DECL(i,j,k),im_opp)=T_out(1)
            enddo
           else if (combine_flag.eq.1) then ! center->centroid
            ! do nothing
           else if (combine_flag.eq.2) then ! extrap where vfrac=0
            ! do nothing
           else
            print *,"combine_flag invalid"
            stop
           endif
          
          else
           print *,"cell_vfrac or combine_flag bust"
           stop
          endif

         enddo ! im=1,nmat

        else
         print *,"project_option invalid"
         stop
        endif

       else if (mask_test.eq.0) then
        ! do nothing
       else
        print *,"mask_test invalid"
        stop
       endif
 
      enddo
      enddo
      enddo ! i,j,k

      return
      end subroutine FORT_COMBINEVEL

       ! combine_flag==2 (only overwrite if F=0)
      subroutine FORT_COMBINEVELFACE( &
       tid, &
       hflag, &
       facecut_index, &
       icefacecut_index, &
       vofface_index, &
       massface_index, &
       ncphys, &
       nmat, &
       nparts, &
       nparts_def, &
       im_solid_map, &
       nten, &
       combine_idx, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       finest_level, &
       velbc, &
       vof,DIMS(vof), &
       mac,DIMS(mac), &
       xface,DIMS(xface), &
       LS,DIMS(LS), &
       solfab,DIMS(solfab), &
       xlo,dx, &
       dir, &
       cur_time)
      use probf90_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
 
      IMPLICIT NONE

      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in) :: hflag
      INTEGER_T, intent(in) :: facecut_index
      INTEGER_T, intent(in) :: icefacecut_index
      INTEGER_T, intent(in) :: massface_index
      INTEGER_T, intent(in) :: vofface_index
      INTEGER_T, intent(in) :: ncphys

      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nparts
      INTEGER_T, intent(in) :: nparts_def
      INTEGER_T, intent(in) :: im_solid_map(nparts_def)
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: combine_idx

      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level

      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: velbc(SDIM,2,SDIM)

      INTEGER_T, intent(in) :: DIMDEC(vof)
      INTEGER_T, intent(in) :: DIMDEC(mac)
      INTEGER_T, intent(in) :: DIMDEC(xface)
      INTEGER_T, intent(in) :: DIMDEC(LS)
      INTEGER_T, intent(in) :: DIMDEC(solfab)
  
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM) 
      INTEGER_T, intent(in) :: dir
      REAL_T, intent(in) :: cur_time

      REAL_T, intent(in) :: vof(DIMV(vof),nmat*ngeom_recon)
      REAL_T, intent(inout) :: mac(DIMV(mac))
      REAL_T, intent(in) :: xface(DIMV(xface),ncphys)
      REAL_T, intent(in) :: LS(DIMV(LS),nmat*(SDIM+1))
      REAL_T, intent(in) :: solfab(DIMV(solfab),nparts_def*SDIM)

      INTEGER_T nten_test 
      INTEGER_T ii,jj,kk 
      INTEGER_T i,j,k
      INTEGER_T icell,jcell,kcell
      INTEGER_T side
      INTEGER_T im
      INTEGER_T idx
      INTEGER_T nhalf

      REAL_T ucombine
      REAL_T xstenMAC(-3:3,SDIM)
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T is_solid_face
      REAL_T vsol
      INTEGER_T iboundary
      INTEGER_T at_RZ_face

      REAL_T face_vfrac_cell(nmat)
      REAL_T face_vfrac(2,nmat)
      REAL_T face_mfrac(2,nmat)
      REAL_T fluid_vfrac_face
      REAL_T total_vol_face
      REAL_T total_vol_face_cell
      REAL_T local_volume
      REAL_T local_mass
      INTEGER_T vofcomp
      INTEGER_T velcomp
      REAL_T mofdata(nmat*ngeom_recon)
      INTEGER_T nmax
      INTEGER_T im_solid_crit,imL,imR
      INTEGER_T partid
      INTEGER_T partid_crit
      INTEGER_T partidL,partidR
      REAL_T lsleft(nmat)
      REAL_T lsright(nmat)
      REAL_T xface_local
      INTEGER_T dir_local
      REAL_T xclamped_face(SDIM)
      REAL_T xclamped_minus(SDIM)
      REAL_T xclamped_plus(SDIM)
      REAL_T LS_clamped_minus
      REAL_T LS_clamped_plus
      REAL_T vel_clamped_minus(SDIM)
      REAL_T vel_clamped_plus(SDIM)
      REAL_T vel_clamped_face(SDIM)
      REAL_T temperature_clamped_minus
      REAL_T temperature_clamped_plus
      INTEGER_T is_clamped_face
      INTEGER_T local_tessellate

      nhalf=3
      nmax=POLYGON_LIST_MAX ! in: COMBINEVELFACE

      if ((tid.lt.0).or. &
          (tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      if (cur_time.ge.zero) then
       ! do nothing
      else
       print *,"cur_time invalid"
       stop
      endif

      if ((hflag.ne.0).and.(hflag.ne.1)) then
       print *,"hflag invalid5 hflag=",hflag
       stop
      endif 
      if (combine_idx.lt.-1) then
       print *,"combine_idx invalid"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid76"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid combinevelface nten, nten_test ",nten,nten_test
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid combinevel face"
       stop
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      if ((nparts.lt.0).or.(nparts.gt.nmat)) then
       print *,"nparts invalid FORT_COMBINEVELFACE"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.nmat)) then
       print *,"nparts_def invalid FORT_COMBINEVELFACE"
       stop
      endif

       ! indexes start at 0
      if ((facecut_index.ne.3).or. &
          (icefacecut_index.ne.4).or. &
          (vofface_index.ne.massface_index+2*nmat)) then
       print *,"face_index bust 1"
       stop
      endif
      if (ncphys.ne.vofface_index+2*nmat) then
       print *,"ncphys invalid"
       stop
      endif

      ii=0
      jj=0
      kk=0
      if (dir.eq.0) then
       ii=1
      else if (dir.eq.1) then
       jj=1
      else if ((dir.eq.2).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"dir out of range in COMBINEVELFACE"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(vof),1,-1,1273)
      call checkbound(fablo,fabhi,DIMS(mac),0,dir,1273)
      call checkbound(fablo,fabhi,DIMS(xface),0,dir,1273)
      call checkbound(fablo,fabhi,DIMS(solfab),0,dir,1276)
      call checkbound(fablo,fabhi,DIMS(LS),1,-1,1276)

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
       growlo,growhi,0,dir,33)

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,dir,50)
       do dir_local=1,SDIM
        xclamped_face(dir_local)=xstenMAC(0,dir_local)
       enddo

       do side=1,2
        if (side.eq.1) then
         icell=i-ii
         jcell=j-jj
         kcell=k-kk
        else if (side.eq.2) then
         icell=i
         jcell=j
         kcell=k
        else
         print *,"side invalid"
         stop
        endif

        call gridsten_level(xsten,icell,jcell,kcell,level,nhalf)
        do dir_local=1,SDIM
         if (side.eq.1) then
          xclamped_minus(dir_local)=xsten(0,dir_local)
         else if (side.eq.2) then
          xclamped_plus(dir_local)=xsten(0,dir_local)
         else
          print *,"side invalid"
          stop
         endif
        enddo ! dir_local=1..sdim
       enddo ! side=1,2

       if (dir.eq.0) then
        idx=i
       else if (dir.eq.1) then
        idx=j
       else if ((dir.eq.2).and.(SDIM.eq.3)) then
        idx=k
       else
        print *,"dir invalid combinevel"
        stop
       endif

       do im=1,nmat
        face_vfrac_cell(im)=zero
       enddo
       total_vol_face_cell=zero

       fluid_vfrac_face=zero
       total_vol_face=zero

       at_RZ_face=0
       if (levelrz.eq.0) then
        ! do nothing
       else if ((levelrz.eq.1).or. &
                (levelrz.eq.3)) then
        if ((dir.eq.0).and. &
            (xstenMAC(0,1).le.VOFTOL*dx(1))) then
         at_RZ_face=1
        endif
       else
        print *,"levelrz invalid combine vel"
        stop
       endif

       is_solid_face=0

       vsol=zero

       if (at_RZ_face.eq.1) then

        is_solid_face=1
        vsol=zero

       else if (at_RZ_face.eq.0) then

        do side=1,2
 
         if (side.eq.1) then
          icell=i-ii
          jcell=j-jj
          kcell=k-kk
         else if (side.eq.2) then
          icell=i
          jcell=j
          kcell=k
         else
          print *,"side invalid"
          stop
         endif

         call gridsten_level(xsten,icell,jcell,kcell,level,nhalf)

         do im=1,nmat*ngeom_recon
          mofdata(im)=vof(D_DECL(icell,jcell,kcell),im)
         enddo
         local_tessellate=3
         call multi_get_volume_tessellate( &
          local_tessellate, & !  =3
          bfact, &
          dx,xsten,nhalf, &
          mofdata, &
          geom_xtetlist(1,1,1,tid+1), &
          nmax, &
          nmax, &
          nmat,SDIM,3)

         do im=1,nmat

          vofcomp=(im-1)*ngeom_recon+1
          local_volume=mofdata(vofcomp)

          if ((local_volume.ge.-VOFTOL).and. &
              (local_volume.le.one+VOFTOL)) then
           if (local_volume.lt.zero) then
            local_volume=zero
           endif
           if (local_volume.gt.one) then
            local_volume=one
           endif
          else
           print *,"local_volume invalid"
           stop
          endif

          face_vfrac_cell(im)=face_vfrac_cell(im)+local_volume
          total_vol_face_cell=total_vol_face_cell+local_volume

          local_volume=xface(D_DECL(i,j,k),vofface_index+2*(im-1)+side)
          local_mass=xface(D_DECL(i,j,k),massface_index+2*(im-1)+side)

          total_vol_face=total_vol_face+local_volume
          if (is_prescribed(nmat,im).eq.0) then
           fluid_vfrac_face=fluid_vfrac_face+local_volume
          else if (is_prescribed(nmat,im).eq.1) then
           ! do nothing
          else
           print *,"is_prescribed invalid"
           stop
          endif

          face_vfrac(side,im)=local_volume
          face_mfrac(side,im)=local_mass

         enddo ! im=1,nmat

        enddo ! side=1,2

        if (total_vol_face.gt.zero) then
         ! do nothing
        else
         print *,"fluid_vfrac_face ",fluid_vfrac_face
         print *,"total_vol_face ",total_vol_face
         print *,"total_vol_face bust combinevel1" 
         stop
        endif
        if (total_vol_face_cell.gt.zero) then
         ! do nothing
        else
         print *,"total_vol_face_cell invalid"
         stop
        endif

        fluid_vfrac_face=fluid_vfrac_face/total_vol_face
        if ((fluid_vfrac_face.ge.zero).and. &
            (fluid_vfrac_face.le.one+VOFTOL)) then
         ! do nothing
        else
         print *,"fluid_vfrac_face invalid"
         stop
        endif

        im_solid_crit=0
        partid_crit=0
        partid=0

        do im=1,nmat
         face_vfrac_cell(im)=face_vfrac_cell(im)/total_vol_face_cell
         lsleft(im)=LS(D_DECL(i-ii,j-jj,k-kk),im)
         lsright(im)=LS(D_DECL(i,j,k),im)
         if (is_prescribed(nmat,im).eq.1) then
          if (im_solid_crit.eq.0) then
           im_solid_crit=im
           partid_crit=partid
          else if ((im_solid_crit.ge.1).and. &
                   (im_solid_crit.le.nmat)) then
           if (lsleft(im)+lsright(im).ge. &
               lsleft(im_solid_crit)+lsright(im_solid_crit)) then
            im_solid_crit=im
            partid_crit=partid
           endif
          else
           print *,"im_solid_crit invalid combinevelface:",im_solid_crit
           stop
          endif
         else if (is_prescribed(nmat,im).eq.0) then
          ! do nothing
         else
          print *,"is_prescribed(nmat,im) invalid"
          stop
         endif

         if (is_lag_part(nmat,im).eq.1) then
          partid=partid+1
         else if (is_lag_part(nmat,im).eq.0) then
          ! do nothing
         else
          print *,"is_lag_part invalid"
          stop
         endif

        enddo ! im=1..nmat

        if (partid.ne.nparts) then
         print *,"partid.ne.nparts"
         stop
        endif

        call get_primary_material(lsleft,nmat,imL)
        call get_primary_material(lsright,nmat,imR)

        partidL=0
        partidR=0

        if (is_prescribed(nmat,imL).eq.1) then ! is_rigid=1 CTML_FSI_mat=0
         do im=1,imL-1
          if (is_lag_part(nmat,im).eq.1) then
           partidL=partidL+1
          else if (is_lag_part(nmat,im).eq.0) then
           ! do nothing
          else
           print *,"is_lag_part(nmat,im) invalid"
           stop
          endif
         enddo ! im=1..imL-1
         if (im_solid_map(partidL+1)+1.ne.imL) then
          print *,"im_solid_map(partidL+1)+1.ne.imL"
          stop
         endif
        else if (is_prescribed(nmat,imL).eq.0) then
         ! do nothing 
        else
         print *,"is_prescribed(nmat,imL) invalid"
         stop
        endif

        if (is_prescribed(nmat,imR).eq.1) then
         do im=1,imR-1
          if (is_lag_part(nmat,im).eq.1) then
           partidR=partidR+1
          else if (is_lag_part(nmat,im).eq.0) then
           ! do nothing
          else
           print *,"is_lag_part(nmat,im) invalid"
           stop
          endif
         enddo ! im=1..imR-1
         if (im_solid_map(partidR+1)+1.ne.imR) then
          print *,"im_solid_map(partidR+1)+1.ne.imR"
          stop
         endif
        else if (is_prescribed(nmat,imR).eq.0) then
         ! do nothing 
        else
         print *,"is_rigid(nmat,imR) invalid"
         stop
        endif

        if (fluid_vfrac_face.le.VOFTOL_AREAFRAC) then
         is_solid_face=2
        else if ((fluid_vfrac_face.ge.VOFTOL_AREAFRAC).and. &
                 (fluid_vfrac_face.le.one+VOFTOL)) then
         xface_local=xface(D_DECL(i,j,k),facecut_index+1)
         if ((xface_local.ge.zero).and.(xface_local.le.half)) then
          xface_local=zero
          is_solid_face=3
         else if ((xface_local.ge.half).and.(xface_local.le.one)) then
          if ((is_prescribed(nmat,imR).eq.1).or. &
              (is_prescribed(nmat,imL).eq.1)) then
           is_solid_face=3
          else if ((is_prescribed(nmat,imR).eq.0).and. &
                   (is_prescribed(nmat,imL).eq.0)) then
           ! do nothing
          else
           print *,"imR or imL invalid"
           stop
          endif
         else
          print *,"xface_local invalid"
          stop
         endif
        else
         print *,"fluid_vfrac_face invalid"
         stop
        endif

        if (hflag.eq.0) then
         if ((im_solid_crit.ge.1).and. &
             (im_solid_crit.le.nmat)) then
          if (im_solid_map(partid_crit+1)+1.ne.im_solid_crit) then
           print *,"im_solid_map(partid_crit+1)+1.ne.im_solid_crit"
           stop
          endif
          velcomp=partid_crit*SDIM+dir+1
          vsol=solfab(D_DECL(i,j,k),velcomp)
         else if (im_solid_crit.eq.0) then
          vsol=zero
         else
          print *,"im_solid_crit invalid combinevelface2:",im_solid_crit
          stop
         endif
        else if (hflag.eq.1) then
         vsol=zero
        else
         print *,"hflag invalid6 hflag=",hflag
         stop
        endif

       else
        print *,"at_RZ_face invalid"
        stop
       endif

       iboundary=0
       if (idx.eq.fablo(dir+1)) then
        iboundary=1
        side=1
       else if (idx.eq.fabhi(dir+1)+1) then
        iboundary=1
        side=2
       else if ((idx.gt.fablo(dir+1)).and. &
                (idx.lt.fabhi(dir+1)+1)) then
        ! do nothing
       else
        print *,"idx invalid"
        stop
       endif

       if (iboundary.eq.1) then
        if (velbc(dir+1,side,dir+1).eq.REFLECT_ODD) then
         vsol=zero
         is_solid_face=4
        else if (velbc(dir+1,side,dir+1).eq.EXT_DIR) then
         if (hflag.eq.0) then
          call velbc_override(cur_time,dir+1,side,dir+1, &
           vsol, &
           xstenMAC,nhalf,dx,bfact)
         else if (hflag.eq.1) then
          vsol=zero
         else
          print *,"hflag invalid7 hflag=",hflag
          stop
         endif
         is_solid_face=5
        else if ((velbc(dir+1,side,dir+1).eq.INT_DIR).or. &
                 (velbc(dir+1,side,dir+1).eq.REFLECT_EVEN).or. &
                 (velbc(dir+1,side,dir+1).eq.FOEXTRAP)) then
         ! do nothing
        else
         print *,"velbc invalid"
         stop
        endif
       else if (iboundary.eq.0) then
        ! do nothing
       else
        print *,"iboundary invalid"
        stop
       endif

         ! LS>0 if clamped
       call SUB_clamped_LS(xclamped_minus,cur_time,LS_clamped_minus, &
           vel_clamped_minus,temperature_clamped_minus)
       call SUB_clamped_LS(xclamped_plus,cur_time,LS_clamped_plus, &
           vel_clamped_plus,temperature_clamped_plus)
       if ((LS_clamped_minus.ge.zero).or. &
           (LS_clamped_plus.ge.zero)) then
        is_clamped_face=1
        do dir_local=1,SDIM
         if (LS_clamped_minus.lt.zero) then
          vel_clamped_face(dir_local)=vel_clamped_plus(dir_local)
          is_clamped_face=2
         else if (LS_clamped_plus.lt.zero) then
          vel_clamped_face(dir_local)=vel_clamped_minus(dir_local)
          is_clamped_face=3
         else
          vel_clamped_face(dir_local)=half*(vel_clamped_plus(dir_local)+ &
            vel_clamped_minus(dir_local))
         endif
        enddo ! dir_local=1..sdim
       else if ((LS_clamped_minus.lt.zero).and. &
                (LS_clamped_plus.lt.zero)) then
        is_clamped_face=0
       else
        print *,"LS_clamped plus or minus is NaN"
        stop
       endif

       if ((is_clamped_face.eq.1).or. &
           (is_clamped_face.eq.2).or. &
           (is_clamped_face.eq.3)) then
        if (is_solid_face.eq.0) then
         is_solid_face=2
         if (hflag.eq.1) then
          vsol=zero
         else if (hflag.eq.0) then
          vsol=vel_clamped_face(dir+1)
         else
          print *,"hflag invalid"
          stop
         endif
        else if ((is_solid_face.ge.1).and. &
                 (is_solid_face.le.5)) then
         ! do nothing
        else
         print *,"is_solid_face invalid"
         stop
        endif
       else if (is_clamped_face.eq.0) then
        ! do nothing
       else
        print *,"is_clamped_face invalid"
        stop
       endif

       if (is_solid_face.eq.1) then ! RZ face
        ucombine=vsol
       else if (is_solid_face.eq.4) then ! REFLECT_ODD
        ucombine=vsol
       else if (is_solid_face.eq.5) then ! EXT_DIR
        ucombine=vsol
       else if ((is_solid_face.eq.2).or. &
                (is_solid_face.eq.3)) then
        ucombine=vsol
       else if (is_solid_face.eq.0) then
        ucombine=mac(D_DECL(i,j,k))
       else
        print *,"is_solid_face invalid 1 ",is_solid_face
        stop
       endif
       mac(D_DECL(i,j,k))=ucombine

      enddo
      enddo
      enddo

      return
      end subroutine FORT_COMBINEVELFACE

       ! PART I (low order):
       !
       ! itensor_iter==0 => face grad u
       !    tileloop==0  => low order
       !       spectral_loop==0 => low order
       !       spectral_loop==1 => do nothing
       !    tileloop==1  => do nothing
       ! itensor_iter==1 => cell grad u
       !    tileloop==0  => low order
       !
       ! PART II (high order)
       !
       ! tileloop==0 => do nothing
       !
       ! tileloop==1
       !   itensor_iter==0 => face grad u
       !     spectral_loop==0 initialize semflux
       !     spectral_loop==1 flux=average of values shared at face.
       !   itensor_iter==1 => cell grad u  
       !   (only called when spectral_loop==0)
       !     
      subroutine FORT_FACE_GRADIENTS( &
       im_tensor, &
       elastic_partid, &
       im_elastic_map, &
       ns_time_order, &
       divu_outer_sweeps, &
       num_divu_outer_sweeps, &
       SDC_outer_sweeps, &
       tileloop, &
       dir, &
       slab_step, &
       itensor_iter, & ! 0=>gradient face  1=>gradient at cell
       time, &
       temperature_primitive_variable, &
       face_flag, &
       enable_spectral, &
       velbc, &
       spectral_loop, &
       ncfluxreg, &
       semflux,DIMS(semflux), &
       amrsync,DIMS(amrsync), &
       mask0,DIMS(mask0), &  ! mask0=1 if not cov. by finer level or outside. 
       mask3,DIMS(mask3), &  ! 1 at fine-fine ghost cells 0 at other ghost.
       maskSEM,DIMS(maskSEM), &
       faceLS,DIMS(faceLS), & 
       mdata,DIMS(mdata), & 
       tdata,DIMS(tdata), & 
       c_tdata,DIMS(c_tdata), & 
       vel,DIMS(vel), &
       solidx,DIMS(solidx), &
       solidy,DIMS(solidy), &
       solidz,DIMS(solidz), &
       levelpc,DIMS(levelpc), &
       recon,DIMS(recon), &
       xlo,dx, &
       rzflag, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact,bfact_c,bfact_f, &
       level, &
       finest_level, &
       nmat, &
       nparts, &
       nparts_def, &
       im_solid_map, &
       homflag, &
       ntensor, &
       SEM_upwind, &
       SEM_advection_algorithm, &
       simple_AMR_BC_flag_viscosity)
      use probcommon_module
      use global_utility_module
      use probf90_module
      use MOF_routines_module
 
      IMPLICIT NONE

      REAL_T def_dt

      INTEGER_T, intent(in) :: im_tensor
      INTEGER_T, intent(in) :: elastic_partid
      INTEGER_T, intent(in) :: im_elastic_map(num_materials_viscoelastic)
      INTEGER_T, intent(in) :: ntensor
      INTEGER_T, intent(in) :: ns_time_order
      INTEGER_T, intent(in) :: divu_outer_sweeps
      INTEGER_T, intent(in) :: num_divu_outer_sweeps
      INTEGER_T, intent(in) :: SDC_outer_sweeps
      INTEGER_T, intent(in) :: tileloop
      INTEGER_T, intent(in) :: dir
      INTEGER_T, intent(in) :: itensor_iter
      INTEGER_T, intent(in) :: spectral_loop
      INTEGER_T, intent(in) :: slab_step
      INTEGER_T, intent(in) :: enable_spectral
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level

      INTEGER_T, intent(in) :: SEM_upwind
      INTEGER_T, intent(in) :: SEM_advection_algorithm
      INTEGER_T, intent(in) :: simple_AMR_BC_flag_viscosity

      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nparts
      INTEGER_T, intent(in) :: nparts_def
      INTEGER_T, intent(in) :: im_solid_map(nparts_def)
      INTEGER_T, intent(in) :: ncfluxreg
      INTEGER_T, intent(in) :: rzflag 
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact,bfact_c,bfact_f
      INTEGER_T, intent(in) :: DIMDEC(semflux)
      INTEGER_T, intent(in) :: DIMDEC(amrsync)
      INTEGER_T, intent(in) :: DIMDEC(maskSEM)
      INTEGER_T, intent(in) :: DIMDEC(mask0)
      !1 at fine-fine ghost cells 0 at other ghost.
      INTEGER_T, intent(in) :: DIMDEC(mask3) 
      INTEGER_T, intent(in) :: DIMDEC(faceLS)
      INTEGER_T, intent(in) :: DIMDEC(mdata)
      INTEGER_T, intent(in) :: DIMDEC(tdata)
      INTEGER_T, intent(in) :: DIMDEC(c_tdata)
      INTEGER_T, intent(in) :: DIMDEC(vel)
      INTEGER_T, intent(in) :: DIMDEC(solidx)
      INTEGER_T, intent(in) :: DIMDEC(solidy)
      INTEGER_T, intent(in) :: DIMDEC(solidz)
      INTEGER_T, intent(in) :: DIMDEC(levelpc)
      INTEGER_T, intent(in) :: DIMDEC(recon)
 
      INTEGER_T, intent(in) :: velbc(SDIM,2,SDIM)
      REAL_T, intent(in) :: time
      INTEGER_T, intent(in) :: temperature_primitive_variable(nmat) 
 
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM) 
      REAL_T, intent(inout) :: semflux(DIMV(semflux),ncfluxreg)
       ! intent(inout) instead of intent(in) since
       ! this parameter doubles as "xp" in SEM_CELL_TO_MAC
      REAL_T, intent(inout) :: amrsync(DIMV(amrsync),SDIM)

      REAL_T, intent(in) :: maskSEM(DIMV(maskSEM))
       ! mask0=tag if not covered by level+1 or outside the domain.
      REAL_T, intent(in) :: mask0(DIMV(mask0))
       ! mask3=tag at exterior fine/fine border.
       ! mask3=1-tag at other exterior boundaries.
      REAL_T, intent(in) :: mask3(DIMV(mask3))
      REAL_T, intent(inout) :: faceLS(DIMV(faceLS),SDIM)
      REAL_T, intent(out) :: mdata(DIMV(mdata),SDIM)
      REAL_T, intent(inout) :: tdata(DIMV(tdata),ntensor)
      REAL_T, intent(out) :: c_tdata(DIMV(c_tdata),ntensor)
      REAL_T, intent(in) :: vel(DIMV(vel),SDIM)
      REAL_T, intent(in) :: solidx(DIMV(solidx),nparts_def*SDIM)
      REAL_T, intent(in) :: solidy(DIMV(solidy),nparts_def*SDIM)
      REAL_T, intent(in) :: solidz(DIMV(solidz),nparts_def*SDIM)
      REAL_T, intent(in) :: levelpc(DIMV(levelpc),nmat)
      REAL_T, intent(in) :: recon(DIMV(recon),nmat*ngeom_recon)
  
      INTEGER_T i,j,k
      INTEGER_T dir2
      INTEGER_T dirtan
      INTEGER_T mask_boundary
      INTEGER_T local_maskSEM
      INTEGER_T maskcov
      INTEGER_T index_adjoin(2,3)
      INTEGER_T masktest,bctest,icrit
      INTEGER_T ii,jj,kk
      INTEGER_T ishift,jshift,kshift
      INTEGER_T shift_flag
      INTEGER_T im1,jm1,km1
      INTEGER_T ip1,jp1,kp1
      INTEGER_T nc
      REAL_T delta
      INTEGER_T side
      INTEGER_T ux,vx,wx,uy,vy,wy,uz,vz,wz,nbase
      REAL_T vplus(SDIM)
      REAL_T vminus(SDIM)

      REAL_T lsleft(nmat)
      REAL_T lsright(nmat)

      REAL_T solidvelleft(SDIM)
      REAL_T solidvelright(SDIM)
      INTEGER_T indexleft,indexright
      INTEGER_T im,imL,imR
      INTEGER_T homflag

      REAL_T hold_grad
      REAL_T RR
      REAL_T xstenMAC(-1:1,SDIM)
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf,nhalfcell
      INTEGER_T scomp,scomp_bc,dcomp
      INTEGER_T ncomp_source
      INTEGER_T ncomp_dest
      INTEGER_T ncomp_xvel
      INTEGER_T ncomp_denold
      INTEGER_T ncomp_veldest
      INTEGER_T ncomp_dendest
      INTEGER_T ncomp_cterm
      INTEGER_T operation_flag
      INTEGER_T energyflag 
      INTEGER_T face_flag 
      INTEGER_T project_option_vel
      REAL_T lspoint(nmat)
      INTEGER_T im_primary
      REAL_T gradu(SDIM)
      REAL_T int_xlo,int_xhi
      REAL_T leftwt,rightwt
      REAL_T slopeLT,slopeRT
      REAL_T theta_factor
      INTEGER_T stripstat
      INTEGER_T ielem,jelem,kelem
      INTEGER_T velcomp
      INTEGER_T tensorcomponent
      INTEGER_T ncomp_xgp
      INTEGER_T ncomp_xp
      INTEGER_T left_rigid,right_rigid
      INTEGER_T partidL,partidR
      INTEGER_T conservative_div_uu

      REAL_T xclamped_minus_sten(-1:1,SDIM)
      REAL_T xclamped_plus_sten(-1:1,SDIM)
      REAL_T xclamped_minus(SDIM)
      REAL_T xclamped_plus(SDIM)
      REAL_T xclamped_cen(SDIM)
      REAL_T LS_clamped_cen
      REAL_T LS_clamped_minus
      REAL_T LS_clamped_plus
      REAL_T vel_clamped_cen(SDIM)
      REAL_T vel_clamped_minus(SDIM)
      REAL_T vel_clamped_plus(SDIM)
      REAL_T vel_clamped_face(SDIM)
      REAL_T temperature_clamped_cen
      REAL_T temperature_clamped_minus
      REAL_T temperature_clamped_plus
      INTEGER_T is_clamped_face

      nhalf=1
      nhalfcell=3

      if ((dir.lt.1).or.(dir.gt.SDIM)) then
       print *,"dir invalid face gradients"
       stop
      endif

      if ((slab_step.lt.-1).or.(slab_step.gt.bfact_time_order)) then
       print *,"slab_step invalid face gradients "
       stop
      endif
      if (time.ge.zero) then
       ! do nothing
      else
       print *,"time invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid77"
       stop
      endif
      if ((bfact_c.ne.bfact).and.(bfact_c.ne.2*bfact)) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((bfact_f.ne.bfact).and.(bfact.ne.2*bfact_f)) then
       print *,"bfact_f invalid"
       stop
      endif
      if ((spectral_loop.ne.0).and. &
          (spectral_loop.ne.1)) then
       print *,"spectral_loop invalid"
       stop
      endif
      if (itensor_iter.eq.0) then ! face grad U
       if ((ncfluxreg/SDIM)*SDIM.ne.ncfluxreg) then
        print *,"ncfluxreg invalid15 ",ncfluxreg
        stop
       endif
       if (ncfluxreg.ne.ntensor) then
        print *,"ncfluxreg invalid16 ",ncfluxreg
        stop
       endif
      else if (itensor_iter.eq.1) then ! cell grad U
       if (ncfluxreg.lt.SDIM) then
        print *,"ncfluxreg invalid17 ",ncfluxreg
        stop
       endif
       if (spectral_loop.eq.0) then
        ! do nothing
       else
        print *,"spectral_loop invalid"
        print *,"unnecessary to average cell based gradients"
        stop
       endif
      else
       print *,"itensor_iter invalid"
       stop
      endif

      if (im_tensor.eq.-1) then
       ! do nothing
      else if ((im_tensor.ge.0).and.(im_tensor.lt.nmat)) then
       if ((elastic_partid.ge.0).and. &
           (elastic_partid.lt.num_materials_viscoelastic)) then
        if (im_tensor.eq.im_elastic_map(elastic_partid+1)) then
         if (im_elastic_map(elastic_partid+1).eq. &
             fort_im_elastic_map(elastic_partid+1)) then
          if (enable_spectral.eq.0) then
           if (homflag.eq.0) then
            ! do nothing
           else
            print *,"homflag invalid"
            stop
           endif
          else
           print *,"enable_spectral invalid"
           stop
          endif
         else
          print *,"im_elastic_map invalid"
          stop
         endif
        else
         print *,"im_tensor invalid"
         stop
        endif
       else
        print *,"elastic_partid invalid"
        stop
       endif
      else
       print *,"im_tensor invalid"
       stop
      endif

      if ((enable_spectral.lt.0).or. &
          (enable_spectral.gt.3)) then
       print *,"enable_spectral invalid face gradients"
       stop
      endif
 
      if ((simple_AMR_BC_flag_viscosity.eq.0).or. &
          (simple_AMR_BC_flag_viscosity.eq.1)) then
       ! do nothing
      else
       print *,"simple_AMR_BC_flag_viscosity invalid"
       stop
      endif

      if ((SEM_upwind.ne.0).and. &
          (SEM_upwind.ne.1)) then
       print *,"SEM_upwind invalid face gradients"
       stop
      endif
      if ((SEM_advection_algorithm.eq.0).or. &
          (SEM_advection_algorithm.eq.1)) then
       ! do nothing
      else
       print *,"SEM_advection_algorithm invalid"
       stop
      endif
      if ((face_flag.ne.0).and. &
          (face_flag.ne.1)) then
       print *,"face_flag invalid"
       stop
      endif

      do im=1,nmat

       if (fort_denconst(im).le.zero) then
        print *,"denconst invalid"
        stop
       endif

       if ((fort_material_type(im).eq.0).or. &
           (is_rigid(nmat,im).eq.1).or. &
           (fort_material_type(im).eq.999)) then
        if (temperature_primitive_variable(im).ne.1) then
         print *,"temperature_primitive_variable(im) invalid"
         stop
        endif
       else if ((fort_material_type(im).gt.0).and. &
                (is_rigid(nmat,im).eq.0).and. &
                (fort_material_type(im).ne.999)) then
        if ((temperature_primitive_variable(im).eq.0).or. &
            (temperature_primitive_variable(im).eq.1)) then
         ! do nothing
        else
         print *,"temperature_primitive_variable(im) invalid"
         stop
        endif
       else
        print *,"fort_material_type(im) or is_rigid invalid"
        stop
       endif

      enddo ! im=1..nmat

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
       ! in: FORT_FACE_GRADIENTS
      if ((ns_time_order.ge.1).and.(ns_time_order.le.32)) then
       ! do nothing
      else
       print *,"ns_time_order invalid"
       stop
      endif
      if ((divu_outer_sweeps.ge.0).and. &
          (divu_outer_sweeps.lt.num_divu_outer_sweeps)) then
       ! do nothing
      else
       print *,"divu_outer_sweeps invalid FORT_FACE_GRADIENTS"
       stop
      endif
       ! in: FORT_FACE_GRADIENTS
      if ((SDC_outer_sweeps.ge.0).and. &
          (SDC_outer_sweeps.lt.ns_time_order)) then
       ! do nothing
      else
       print *,"SDC_outer_sweeps invalid in FORT_FACE_GRADIENTS"
       stop
      endif

      if (rzflag.eq.0) then
       ! do nothing
      else if (rzflag.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (rzflag.eq.3) then
       ! do nothing
      else
       print *,"rzflag invalid"
       stop
      endif

      if (ntensor.ne.SDIM*SDIM) then
       print *,"ntensor invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(amrsync),0,dir-1,231)

      if ((dir.eq.1).and.(tileloop.eq.0)) then

       call checkbound(fablo,fabhi,DIMS(semflux),1,-1,231)
       call checkbound(fablo,fabhi,DIMS(faceLS),1,-1,1263)
       call checkbound(fablo,fabhi,DIMS(mask0),1,-1,1264)
       call checkbound(fablo,fabhi,DIMS(mask3),1,-1,1264)
       call checkbound(fablo,fabhi,DIMS(mdata),1,-1,1264)
       call checkbound(fablo,fabhi,DIMS(tdata),1,-1,1265)
       call checkbound(fablo,fabhi,DIMS(c_tdata),1,-1,1265)
       call checkbound(fablo,fabhi,DIMS(vel),1,-1,1266)
       call checkbound(fablo,fabhi,DIMS(solidx),0,0,1267)
       call checkbound(fablo,fabhi,DIMS(solidy),0,1,1267)
       call checkbound(fablo,fabhi,DIMS(solidz),0,SDIM-1,1267)
       call checkbound(fablo,fabhi,DIMS(levelpc),2,-1,1368)
       call checkbound(fablo,fabhi,DIMS(recon),2,-1,1368)
       call checkbound(fablo,fabhi,DIMS(maskSEM),1,-1,1264)

      endif

      call tensorcomp_matrix(ux,uy,uz,vx,vy,vz,wx,wy,wz)

! compute u_x,v_x,w_x, u_y,v_y,w_y, u_z,v_z,w_z;  

      ii=0
      jj=0
      kk=0
      if (dir.eq.1) then
       ii=1
       nbase=ux-1
      else if (dir.eq.2) then
       jj=1
       nbase=uy-1
      else if ((dir.eq.SDIM).and.(SDIM.eq.3)) then
       kk=1
       nbase=uz-1
      else
       print *,"dir invalid face gradients 2"
       stop
      endif

      if (itensor_iter.eq.0) then  ! face grad U

       if (tileloop.eq.0) then ! low order

        if (spectral_loop.eq.0) then

          ! same as growntileboxMAC, except includes one layer of ghost
          ! cells in the tangential directions.
         call growntileboxTENSOR(tilelo,tilehi,fablo,fabhi, &
           growlo,growhi,dir-1)

         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)

          !dir=1..sdim
          call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,dir-1,51)

          im1=i-ii
          jm1=j-jj
          km1=k-kk

           ! kluge to prevent reading out of the array
           ! bounds for coupling terms that have solid cells
           ! in the stencil.  
           ! The corner values will never be used if cell pair is
           ! owned by a different material (or not owned at all)
           ! from the material owned on the face.
           ! If the cell pair is owned by a prescribed solid, then
           ! it is not used in the coupling stencil.
          shift_flag=0
          ishift=i
          jshift=j
          kshift=k
          if ((dir.eq.1).or. &
              (dir.eq.3)) then
           dirtan=2
           if (j.eq.fablo(dirtan)-1) then
            jshift=j+1
            shift_flag=1
           else if (j.eq.fabhi(dirtan)+1) then
            jshift=j-1
            shift_flag=1
           else if ((j.ge.fablo(dirtan)).and. &
                    (j.le.fabhi(dirtan))) then
            ! do nothing
           else
            print *,"j invalid"
            stop
           endif
          endif

          if ((dir.eq.2).or. &
              (dir.eq.3)) then
           dirtan=1
           if (i.eq.fablo(dirtan)-1) then
            ishift=i+1
            shift_flag=1
           else if (i.eq.fabhi(dirtan)+1) then
            ishift=i-1
            shift_flag=1
           else if ((i.ge.fablo(dirtan)).and. &
                    (i.le.fabhi(dirtan))) then
            ! do nothing
           else
            print *,"i invalid"
            stop
           endif
          endif

          if (((dir.eq.1).or. &
               (dir.eq.2)).and.  &
              (SDIM.eq.3)) then
           dirtan=SDIM
           if (k.eq.fablo(dirtan)-1) then
            kshift=k+1
            shift_flag=1
           else if (k.eq.fabhi(dirtan)+1) then
            kshift=k-1
            shift_flag=1
           else if ((k.ge.fablo(dirtan)).and. &
                    (k.le.fabhi(dirtan))) then
            ! do nothing
           else
            print *,"k invalid FORT_FACE_GRADIENTS"
            stop
           endif
          endif

          side=1
          index_adjoin(side,1)=im1
          index_adjoin(side,2)=jm1
          index_adjoin(side,3)=km1
          side=2
          index_adjoin(side,1)=i
          index_adjoin(side,2)=j
          index_adjoin(side,3)=k

          mask_boundary=0

          do dir2=1,SDIM
           do side=1,2
            if (side.eq.1) then
             ! mask3=tag at exterior fine/fine border.
             ! mask3=1-tag at other exterior boundaries.
             ! mask0=tag if not covered by level+1 or outside the domain.
             masktest=NINT(mask3(D_DECL(im1,jm1,km1)))
             maskcov=NINT(mask0(D_DECL(im1,jm1,km1)))
            else if (side.eq.2) then
             masktest=NINT(mask3(D_DECL(i,j,k)))
             maskcov=NINT(mask0(D_DECL(i,j,k)))
            else
             print *,"side invalid"
             stop
            endif
            if (maskcov.eq.1) then
             ! do nothing
            else if (maskcov.eq.0) then
             mask_boundary=1 ! do not include covered cells in the stencil.
            else
             print *,"maskcov invalid"
             stop
            endif

            icrit=index_adjoin(side,dir2)

            if (icrit.lt.fablo(dir2)) then
             if (masktest.eq.0) then ! not fine/fine
              do nc=1,SDIM
               bctest=velbc(dir2,1,nc)
               if ((bctest.eq.INT_DIR).or. &!coarse/fine not fine/fine,periodic
                   (bctest.eq.FOEXTRAP).or. &
                   (bctest.eq.EXT_DIR).or. &
                   (bctest.eq.REFLECT_ODD).or. &
                   (bctest.eq.REFLECT_EVEN)) then
                mask_boundary=1
               else
                print *,"bctest invalid"
                stop
               endif
              enddo ! nc
             else if (masktest.eq.1) then
              ! do nothing fine-fine border
             else
              print *,"masktest invalid"
              stop
             endif
            endif
            if (icrit.gt.fabhi(dir2)) then
             if (masktest.eq.0) then
              do nc=1,SDIM
               bctest=velbc(dir2,2,nc)
               if ((bctest.eq.INT_DIR).or. &
                   (bctest.eq.FOEXTRAP).or. &
                   (bctest.eq.EXT_DIR).or. &
                   (bctest.eq.REFLECT_ODD).or. &
                   (bctest.eq.REFLECT_EVEN)) then
                mask_boundary=1
               else
                print *,"bctest invalid"
                stop
               endif
              enddo ! nc
             else if (masktest.eq.1) then
              ! do nothing fine-fine border
             else
              print *,"masktest invalid"
              stop
             endif
            endif
           enddo ! side
          enddo ! dir2
            
          mdata(D_DECL(i,j,k),dir)=one

          do im=1,nmat
           lsleft(im)=levelpc(D_DECL(im1,jm1,km1),im)
           lsright(im)=levelpc(D_DECL(i,j,k),im)
          enddo ! im=1..nmat

          call get_primary_material(lsleft,nmat,imL)
          call get_primary_material(lsright,nmat,imR)

          call gridsten_level(xclamped_minus_sten,im1,jm1,km1,level,nhalf)
          call gridsten_level(xclamped_plus_sten,i,j,k,level,nhalf)
          do dir2=1,SDIM
           xclamped_minus(dir2)=xclamped_minus_sten(0,dir2)
           xclamped_plus(dir2)=xclamped_plus_sten(0,dir2)
          enddo

           ! LS>0 if clamped
          call SUB_clamped_LS(xclamped_minus,time,LS_clamped_minus, &
                vel_clamped_minus,temperature_clamped_minus)
          call SUB_clamped_LS(xclamped_plus,time,LS_clamped_plus, &
                vel_clamped_plus,temperature_clamped_plus)
          if ((LS_clamped_minus.ge.zero).or. &
              (LS_clamped_plus.ge.zero)) then
           is_clamped_face=1
           do dir2=1,SDIM
            if (LS_clamped_minus.lt.zero) then
             vel_clamped_face(dir2)=vel_clamped_plus(dir2)
             is_clamped_face=2
            else if (LS_clamped_plus.lt.zero) then
             vel_clamped_face(dir2)=vel_clamped_minus(dir2)
             is_clamped_face=3
            else
             vel_clamped_face(dir2)=half*(vel_clamped_plus(dir2)+ &
               vel_clamped_minus(dir2))
            endif
           enddo
          else if ((LS_clamped_minus.lt.zero).and. &
                   (LS_clamped_plus.lt.zero)) then
           is_clamped_face=0
          else
           print *,"LS_clamped plus or minus is NaN"
           stop
          endif


          if ((is_clamped_face.eq.1).or. &
              (is_clamped_face.eq.2).or. &
              (is_clamped_face.eq.3)) then
           faceLS(D_DECL(i,j,k),dir)=zero ! mask off this gradient
          else if (is_clamped_face.eq.0) then
           if (mask_boundary.eq.0) then
            if (imL.eq.imR) then
             faceLS(D_DECL(i,j,k),dir)=imL
            else
             faceLS(D_DECL(i,j,k),dir)=zero ! mask off this gradient
            endif
           else if (mask_boundary.eq.1) then
            faceLS(D_DECL(i,j,k),dir)=zero ! mask off this gradient
           else
            print *,"mask_boundary invalid"
            stop
           endif
          else
           print *,"is_clamped_face invalid"
           stop
          endif

          indexleft=index_adjoin(1,dir)
          indexright=index_adjoin(2,dir)

          partidL=0
          partidR=0

          if (is_prescribed(nmat,imL).eq.1) then
           do im=1,imL-1
            if (is_lag_part(nmat,im).eq.1) then
             partidL=partidL+1
            else if (is_lag_part(nmat,im).eq.0) then
             ! do nothing
            else
             print *,"is_lag_part(nmat,im) invalid"
             stop
            endif
           enddo ! im=1..imL-1
           if (im_solid_map(partidL+1)+1.ne.imL) then
            print *,"im_solid_map(partidL+1)+1.ne.imL"
            stop
           endif
          else if (is_prescribed(nmat,imL).eq.0) then
           ! do nothing 
          else
           print *,"is_prescribed(nmat,imL) invalid"
           stop
          endif

          if (is_prescribed(nmat,imR).eq.1) then
           do im=1,imR-1
            if (is_lag_part(nmat,im).eq.1) then
             partidR=partidR+1
            else if (is_lag_part(nmat,im).eq.0) then
             ! do nothing
            else
             print *,"is_lag_part(nmat,im) invalid"
             stop
            endif
           enddo ! im=1..imR-1
           if (im_solid_map(partidR+1)+1.ne.imR) then
            print *,"im_solid_map(partidR+1)+1.ne.imR"
            stop
           endif
          else if (is_prescribed(nmat,imR).eq.0) then
           ! do nothing 
          else
           print *,"is_prescribed(nmat,imR) invalid"
           stop
          endif
 
          do nc=1,SDIM

           velcomp=nc
           left_rigid=0
           right_rigid=0

           if (is_prescribed(nmat,imL).eq.1) then
            left_rigid=1

            if (dir.eq.1) then
             solidvelleft(nc)= &
               solidx(D_DECL(ishift,jshift,kshift),partidL*SDIM+nc)
            else if (dir.eq.2) then
             solidvelleft(nc)= &
               solidy(D_DECL(ishift,jshift,kshift),partidL*SDIM+nc)
            else if ((dir.eq.3).and.(SDIM.eq.3)) then
             solidvelleft(nc)= &
               solidz(D_DECL(ishift,jshift,kshift),partidL*SDIM+nc)
            else
             print *,"dir invalid"
             stop
            endif
           else if (is_prescribed(nmat,imL).eq.0) then
            solidvelleft(nc)=zero
           else
            print *,"is_prescribed(nmat,imL) invalid"
            stop
           endif

           if (is_prescribed(nmat,imR).eq.1) then
            right_rigid=1

            if (dir.eq.1) then
             solidvelright(nc)= &
               solidx(D_DECL(ishift,jshift,kshift),partidR*SDIM+nc)
            else if (dir.eq.2) then
             solidvelright(nc)= &
               solidy(D_DECL(ishift,jshift,kshift),partidR*SDIM+nc)
            else if ((dir.eq.3).and.(SDIM.eq.3)) then
             solidvelright(nc)= &
               solidz(D_DECL(ishift,jshift,kshift),partidR*SDIM+nc)
            else
             print *,"dir invalid"
             stop
            endif

           else if (is_prescribed(nmat,imR).eq.0) then
            solidvelright(nc)=zero
           else
            print *,"is_prescribed(nmat,imR) invalid"
            stop
           endif

           if ((indexleft.lt.fablo(dir)).and. &
               (velbc(dir,1,nc).eq.EXT_DIR).and. &
               (is_prescribed(nmat,imL).eq.0)) then
            left_rigid=1
            solidvelleft(nc)=vel(D_DECL(im1,jm1,km1),velcomp)
           endif

           if ((indexright.gt.fabhi(dir)).and. &
               (velbc(dir,2,nc).eq.EXT_DIR).and. &
               (is_prescribed(nmat,imR).eq.0)) then
            right_rigid=1
            solidvelright(nc)=vel(D_DECL(i,j,k),velcomp)
           endif

           if (is_clamped_face.eq.1) then
            left_rigid=1
            right_rigid=1
            solidvelleft(nc)=vel_clamped_face(nc)
            solidvelright(nc)=vel_clamped_face(nc)
           else if (is_clamped_face.eq.2) then
            right_rigid=1
            solidvelright(nc)=vel_clamped_face(nc)
           else if (is_clamped_face.eq.3) then
            left_rigid=1
            solidvelleft(nc)=vel_clamped_face(nc)
           else if (is_clamped_face.eq.0) then
            ! do nothing
           else
            print *,"is_clamped_face invalid"
            stop
           endif

           if (homflag.eq.1) then
            solidvelleft(nc)=zero
            solidvelright(nc)=zero
           else if (homflag.eq.0) then
            ! do nothing
           else
            print *,"homflag invalid 1"
            stop
           endif

           ! lsleft_solid and lsright_solid < 0 in the solid.
           if ((left_rigid.eq.1).and.(right_rigid.eq.1)) then
            mdata(D_DECL(i,j,k),dir)=zero
           else if ((left_rigid.eq.1).or.(right_rigid.eq.1)) then
            if (shift_flag.eq.0) then
             ! do nothing

             !prevent reading solidvel outside the array bounds
            else if (shift_flag.eq.1) then 
             mdata(D_DECL(i,j,k),dir)=zero
            else
             print *,"shift_flag invalid"
             stop
            endif
           else if ((left_rigid.eq.0).and.(right_rigid.eq.0)) then
            ! do nothing
           else
            print *,"left_rigid or right_rigid invalid"
            stop
           endif

           ! extra factor of r for theta gradient in cylindrical coordinates.
           RR=one
           if (rzflag.eq.0) then
            ! do nothing
           else if (rzflag.eq.1) then
            if (SDIM.ne.2) then
             print *,"dimension bust"
             stop
            endif
            ! do nothing
           else if (rzflag.eq.3) then
            if (dir.eq.2) then ! theta direction  s_theta/r
             RR=xstenMAC(0,1)
            endif
           else
            print *,"rzflag invalid"
            stop
           endif
           if (RR.le.zero) then
            print *,"RR invalid"
            stop
           else if (RR.gt.zero) then
            delta=RR*(xstenMAC(1,dir)-xstenMAC(-1,dir))
           else
            print *,"RR bust"
            stop
           endif
           if (delta.le.zero) then
            print *,"delta invalid"
            stop
           endif

           theta_factor=one

            ! gradient of displacement vector.
           if ((im_tensor.ge.0).and.(im_tensor.lt.nmat)) then

            vplus(nc)=vel(D_DECL(i,j,k),velcomp)
            vminus(nc)=vel(D_DECL(im1,jm1,km1),velcomp)
            if ((left_rigid.eq.0).and.(right_rigid.eq.0)) then
             ! do nothing
            else if ((left_rigid.eq.1).and.(right_rigid.eq.1)) then
             theta_factor=zero ! zero out the gradient.
            else if ((left_rigid.eq.0).and.(right_rigid.eq.1)) then
             theta_factor=zero ! zero out the gradient.
            else if ((left_rigid.eq.1).and.(right_rigid.eq.0)) then
             theta_factor=zero ! zero out the gradient.
            else
             print *,"left_rigid or right_rigid invalid"
             stop
            endif

            ! gradient of velocity vector.
           else if (im_tensor.eq.-1) then

            if ((left_rigid.eq.0).and.(right_rigid.eq.0)) then
             vplus(nc)=vel(D_DECL(i,j,k),velcomp)
             vminus(nc)=vel(D_DECL(im1,jm1,km1),velcomp)
            else if ((left_rigid.eq.1).and.(right_rigid.eq.1)) then
             vplus(nc)=solidvelright(nc)
             vminus(nc)=solidvelleft(nc)
             theta_factor=zero ! zero out the gradient.
            else if ((left_rigid.eq.0).and.(right_rigid.eq.1)) then
             vplus(nc)=solidvelright(nc)
             vminus(nc)=vel(D_DECL(im1,jm1,km1),velcomp)
             theta_factor=two
            else if ((left_rigid.eq.1).and.(right_rigid.eq.0)) then
             vplus(nc)=vel(D_DECL(i,j,k),velcomp)
             vminus(nc)=solidvelleft(nc)
             theta_factor=two
            else
             print *,"left_rigid or right_rigid invalid"
             stop
            endif

           else
            print *,"im_tensor invalid"
            stop
           endif


 !ux,vx,wx,uy,vy,wy,uz,vz,wz
 ! grad u in cylindrical coordinates:
 !
 ! S= (grad u + grad u^T)/2 
 !
 ! grad u=| u_r  u_t/r-v/r  u_z  |
 !        | v_r  v_t/r+u/r  v_z  |
 !        | w_r  w_t/r      w_z  |
 !
 ! S=
 !   |u_r     (u_t/r+v_r-v/r)/2   (u_z+w_r)/2   |
 !   | .      v_t/r+u/r           (v_z+w_t/r)/2 |
 !   | .      .                   w_z           |
 !
 ! note: v_r-v/r=r(v/r)_r
 !
 ! 2S=
 !
 !   |2u_r     (u_t/r+v_r-v/r)   (u_z+w_r)   |
 !   | .       2v_t/r+2u/r       (v_z+w_t/r) |
 !   | .      .                    2w_z      |
 !
 ! 
 ! div S = | (r S_11)_r/r + (S_12)_t/r - S_22/r  + (S_13)_z |
 !         | (r S_21)_r/r + (S_22)_t/r + S_12/r  + (S_23)_z |
 !         | (r S_31)_r/r + (S_32)_t/r +           (S_33)_z |
 ! 
 ! ur =  costheta u + sintheta v
 ! ut = -sintheta u + costheta v
 ! 
 ! u = costheta ur - sintheta ut
 ! v = sintheta ur + costheta ut
 !
 ! e.g. theta=pi/2  ur=0 
 !   u=-utheta  v=0
 ! if constant viscosity:
 ! div u = (ru)_r/r + v_t/r + w_z= u_r + u/r +v_t/r + w_z=0
 ! (div u)_r=u_rr+u_r/r-u/r^2+v_tr/r-v_t/r^2+w_zr=0
 ! (div u)_t/r=u_rt/r+u_t/r^2+v_tt/r^2+w_zt/r=0
 ! (div u)_z=u_rz+u_z/r+v_tz/r+w_zz=0
 !
 ! div(2 S)=
 ! |2u_rr+2u_r/r+u_tt/r^2+v_tr/r-v_t/r^2-2v_t/r^2-2u/r^2+u_zz+w_rz |
 ! |u_tr/r+v_rr+v_r/r-v_r/r+2v_tt/r^2+2u_t/r^2+u_t/r^2+v_r/r-v/r^2+v_zz+w_tz/r|
 ! |u_zr+u_z/r+w_rr+w_r/r+v_zt/r+w_tt/r^2 + 2w_zz |=
 !
 ! |u_rr+u_r/r-u/r^2+u_tt/r^2-2v_t/r^2+u_zz |    
 ! |v_rr+v_r/r+v_tt/r^2+2u_t/r^2-v/r^2+v_zz |
 ! |w_rr+w_r/r+w_tt/r^2+w_zz                |
 !
 ! compromise: 
 !
 ! GU=| u_r       u_t/r  u_z  |
 !    | v_r       v_t/r  v_z  |
 !    | w_r       w_t/r  w_z  |
 !
 ! hoop term 1st component:  -3 v_t/r^2 - 2 u/r^2
 ! hoop term 2nd component:   3 u_t/r^2 - v/r^2
 ! 
 ! If constant_viscosity==true:
 ! hoop term 1st component:  -2 v_t/r^2 - u/r^2
 ! hoop term 2nd component:   2 u_t/r^2 - v/r^2
 ! No coupling terms.
 ! Diagonal terms not multiplied by 2.

           hold_grad=theta_factor*(vplus(nc)-vminus(nc))/delta

           tensorcomponent=nbase+nc

           tdata(D_DECL(i,j,k),tensorcomponent)=hold_grad

          enddo ! nc=1..sdim

         enddo
         enddo
         enddo  ! i,j,k faces

        else if (spectral_loop.eq.1) then
         ! do nothing
        else
         print *,"spectral_loop invalid"
         stop
        endif

       else if (tileloop.eq.1) then ! high order (see below)
        ! do nothing
       else
        print *,"tileloop invalid"
        stop
       endif

      else if (itensor_iter.eq.1) then  ! cell grad U

       if (spectral_loop.eq.0) then
        ! do nothing
       else
        print *,"spectral_loop invalid"
        stop
       endif

       if (tileloop.eq.0) then ! low order gradients at cells

         ! same as growntilebox, except includes one layer of ghost
         ! cells in the tangential directions.
        call growntileboxTENSOR_SEM( &
         tilelo,tilehi,fablo,fabhi,growlo,growhi,dir-1) 

        do i=growlo(1),growhi(1)
        do j=growlo(2),growhi(2)
        do k=growlo(3),growhi(3)
         call gridsten_level(xsten,i,j,k,level,nhalfcell)
         do dir2=1,SDIM
          xclamped_cen(dir2)=xsten(0,dir2)
         enddo

          ! LS>0 if clamped
         call SUB_clamped_LS(xclamped_cen,time,LS_clamped_cen, &
                vel_clamped_cen,temperature_clamped_cen)

         im1=i-ii
         jm1=j-jj
         km1=k-kk
         ip1=i+ii
         jp1=j+jj
         kp1=k+kk
         do im=1,nmat
          lspoint(im)=levelpc(D_DECL(i,j,k),im)
         enddo
         call get_primary_material(lspoint,nmat,im_primary)
         if ((im_primary.lt.1).or.(im_primary.gt.nmat)) then
          print *,"im_primary invalid"
          stop
         endif

         if (LS_clamped_cen.ge.zero) then
          do nc=1,SDIM
           gradu(nc)=zero
          enddo
         else if (LS_clamped_cen.lt.zero) then

          if (is_prescribed(nmat,im_primary).eq.1) then
           do nc=1,SDIM
            gradu(nc)=zero
           enddo
          else if (is_prescribed(nmat,im_primary).eq.0) then
           int_xlo=max(xsten(-2,dir),xsten(-1,dir))
           int_xhi=min(xsten(0,dir),xsten(1,dir))
           if (int_xhi.gt.int_xlo) then
            leftwt=int_xhi-int_xlo
           else
            leftwt=zero
           endif
           int_xlo=max(xsten(0,dir),xsten(-1,dir))
           int_xhi=min(xsten(2,dir),xsten(1,dir))
           if (int_xhi.gt.int_xlo) then
            rightwt=int_xhi-int_xlo
           else
            rightwt=zero
           endif
           if ((leftwt.le.zero).or.(rightwt.le.zero)) then
            print *,"weights invalid"
            stop
           endif
           lsleft(im_primary)=levelpc(D_DECL(im1,jm1,km1),im_primary)
           lsright(im_primary)=levelpc(D_DECL(ip1,jp1,kp1),im_primary)
            
           do nc=1,SDIM
            tensorcomponent=nbase+nc
            slopeLT=tdata(D_DECL(i,j,k),tensorcomponent)
            slopeRT=tdata(D_DECL(ip1,jp1,kp1),tensorcomponent)

            if (((lsleft(im_primary).ge.zero).and. &
                 (lsright(im_primary).ge.zero)).or. &
                ((lsleft(im_primary).lt.zero).and. &
                 (lsright(im_primary).lt.zero))) then
             gradu(nc)=(leftwt*slopeLT+rightwt*slopeRT)/(leftwt+rightwt)
            else if (lsleft(im_primary).ge.zero) then
             gradu(nc)=slopeLT
            else if (lsright(im_primary).ge.zero) then
             gradu(nc)=slopeRT
            else
             print *,"lsleft or lsright invalid"
             stop
            endif

           enddo ! nc=1..sdim
          else
           print *,"is_prescribed(im_primary) invalid"
           stop
          endif

         else
          print *,"LS_clamped_cen is NaN"
          stop
         endif
            
         do nc=1,SDIM
          tensorcomponent=nbase+nc
          c_tdata(D_DECL(i,j,k),tensorcomponent)=gradu(nc)
         enddo ! nc=1..sdim

        enddo ! k
        enddo ! j
        enddo ! i

       else if (tileloop.eq.1) then ! high order below
        ! do nothing
       else
        print *,"tileloop invalid"
        stop
       endif

      else
       print *,"itensor_iter invalid"
       stop
      endif

       ! in: FACE_GRADIENTS
      if ((enable_spectral.eq.1).or. &  ! SEM space and time
          (enable_spectral.eq.2)) then  ! SEM space only

       if (im_tensor.eq.-1) then
        ! do nothing
       else
        print *,"expecting im_tensor==-1 for spectral method"
        stop
       endif

       if (bfact.ge.2) then

        if (tileloop.eq.1) then ! tileloop==1 => high order

         if (itensor_iter.eq.0) then ! tensor on MAC grid.

          if ((dir.lt.1).or.(dir.gt.SDIM)) then
           print *,"dir invalid face gradients 3"
           stop
          endif

          ! same as growntilebox, except includes one layer of ghost
          ! cells in the tangential directions.
          call growntileboxTENSOR_SEM(tilelo,tilehi,fablo,fabhi, &
           growlo,growhi,dir-1)

          do i=growlo(1),growhi(1)
          do j=growlo(2),growhi(2)
          do k=growlo(3),growhi(3)
           local_maskSEM=NINT(maskSEM(D_DECL(i,j,k)))
           maskcov=NINT(mask0(D_DECL(i,j,k)))

           if ((local_maskSEM.ge.1).and. &
               (local_maskSEM.le.nmat).and. &
               (maskcov.eq.1)) then

            call strip_status_dir(i,j,k,bfact,dir-1,stripstat)

             ! stripstat==1 if (idx_dir/bfact)*bfact == idx_dir
             ! In otherwords, one can have stripstat==1 even if
             ! one of the tangential indexes is outside of the box.
            if (stripstat.eq.1) then

              ielem=i
              jelem=j
              kelem=k

              scomp=1  

              dcomp=nbase+1

              ncomp_source=SDIM
              ncomp_dest=SDIM
              ncomp_xgp=ntensor
              ncomp_xp=SDIM ! number of amrsync components
              scomp_bc=1
              operation_flag=6  ! evaluate tensor values
              energyflag=0
              project_option_vel=3
              def_dt=one
              conservative_div_uu=0

               ! in: FORT_FACE_GRADIENTS
               ! the boundary conditions for "vel" are already set in the
               ! ghost cell.  e.g. homogeneous versus inhomogeneous.
              call SEM_CELL_TO_MAC( &
               conservative_div_uu, &
               ncomp_xp, &  ! number of amrsync components
               simple_AMR_BC_flag_viscosity, &
               level, &
               finest_level, &
               nmat, &
               operation_flag, &  ! 6
               energyflag, &
               temperature_primitive_variable, &
               project_option_vel, &
               SEM_upwind, &
               SEM_advection_algorithm, &
               time, &  ! beta
               time, &  ! visc_coef
               time, &
               def_dt, &  ! dt
               ielem,jelem,kelem, &
               tilelo,tilehi, &
               fablo,fabhi, &
               xlo, &
               dx,dir, &
               bfact,bfact_c,bfact_f, &
               velbc, &  ! presbc
               velbc, &
               scomp, &
               scomp_bc, &
               dcomp, &
               ncomp_dest, &
               ncomp_source, &
               ncomp_xgp, &     ! =ntensor
               ncomp_dest, &    ! ncphys
               spectral_loop, &
               ncfluxreg, &
               semflux,DIMS(semflux), &
               mask3,DIMS(mask3), &
               mask0,DIMS(mask0), & !mask0=1 if not cov. by finer or outside.
               vel,DIMS(vel), &
               vel,DIMS(vel), &  ! pres
               vel,DIMS(vel), &  ! den
               tdata,DIMS(tdata), &  !xface
               tdata,DIMS(tdata), &  !xgp (destination)
               tdata,DIMS(tdata), &  !xcut
               amrsync,DIMS(amrsync), & !xp
               tdata,DIMS(tdata), &  !xvel
               maskSEM,DIMS(maskSEM))

            else if (stripstat.eq.0) then
             ! do nothing
            else
             print *,"stripstat invalid"
             stop
            endif
         
           else if ((local_maskSEM.eq.0).or. &
                    (maskcov.eq.0)) then

            ! do nothing

           else
            print *,"local_maskSEM invalid"
            stop
           endif 

          enddo
          enddo
          enddo ! i,j,k

         else if (itensor_iter.eq.1) then ! cell grad U

          if (spectral_loop.eq.0) then
           ! do nothing
          else
           print *,"spectral_loop invalid"
           stop
          endif

          if ((dir.lt.1).or.(dir.gt.SDIM)) then
           print *,"dir invalid face gradients 4"
           stop
          endif

           ! same as growntilebox, except includes one layer of ghost
           ! cells in the tangential directions.
          call growntileboxTENSOR_SEM(tilelo,tilehi,fablo,fabhi, &
           growlo,growhi,dir-1)

          do i=growlo(1),growhi(1)
          do j=growlo(2),growhi(2)
          do k=growlo(3),growhi(3)

           local_maskSEM=NINT(maskSEM(D_DECL(i,j,k)))
           maskcov=NINT(mask0(D_DECL(i,j,k)))

           if ((local_maskSEM.ge.1).and. &
               (local_maskSEM.le.nmat).and. &
               (maskcov.eq.1)) then

            call strip_status_dir(i,j,k,bfact,dir-1,stripstat)

             ! stripstat==1 if (idx_dir/bfact)*bfact == idx_dir
            if (stripstat.eq.1) then

              ielem=i
              jelem=j
              kelem=k

              ! u_x,v_x,w_x, u_y,v_y,w_y, u_z,v_z,w_z;  
              dcomp=nbase+1

              scomp=dcomp

              scomp_bc=1
              ncomp_dest=SDIM ! ncomp
              ncomp_xvel=ntensor
              ncomp_denold=ntensor
              ncomp_veldest=ntensor
              ncomp_dendest=ntensor
              ncomp_cterm=ntensor
              operation_flag=5  ! interpolate grad u from MAC to CELL
              project_option_vel=3
              energyflag=0
              def_dt=one
              conservative_div_uu=0

               ! 1<=dir<=sdim
              call SEM_MAC_TO_CELL( &
               ncomp_denold, &
               ncomp_veldest, &
               ncomp_dendest, &
               conservative_div_uu, &
               ns_time_order, &
               divu_outer_sweeps, &
               num_divu_outer_sweeps, &
               SDC_outer_sweeps, &
               SEM_advection_algorithm, &
               level, &
               finest_level, &
               nmat, &
               operation_flag, & ! operation_flag==5
               project_option_vel, &
               energyflag, &
               temperature_primitive_variable, &
               face_flag, &
               energyflag, & ! homflag
               local_maskSEM, &
               time, &
               slab_step, &
               def_dt, & ! dt
               ielem,jelem,kelem, &
               tilelo,tilehi, &
               fablo,fabhi, &
               xlo,dx,dir,bfact, &
               velbc, &
               velbc, & ! presbc
               scomp, &
               scomp_bc, &
               dcomp, &
               ncomp_dest, & ! ncomp
               ncomp_xvel, &
               ncomp_cterm, &
               tdata,DIMS(tdata), & ! vol
               tdata,DIMS(tdata), & ! xface
               tdata,DIMS(tdata), & ! xp
               tdata,DIMS(tdata), & ! xvel
               tdata,DIMS(tdata), & ! maskcoef
               tdata,DIMS(tdata), & ! cterm
               tdata,DIMS(tdata), & ! mdotcell
               tdata,DIMS(tdata), & ! pold
               tdata,DIMS(tdata), & ! denold
               tdata,DIMS(tdata), & ! ustar
               c_tdata,DIMS(c_tdata), & ! veldest
               c_tdata,DIMS(c_tdata), & ! dendest
               c_tdata,DIMS(c_tdata)) !divdest

            else if (stripstat.eq.0) then
             ! do nothing
            else
             print *,"stripststat invalid"
             stop
            endif
  
           else if ((local_maskSEM.eq.0).or. &
                    (maskcov.eq.0)) then
            ! do nothing
           else
            print *,"local_maskSEM invalid"
            stop
           endif

          enddo
          enddo
          enddo

         else
          print *,"itensor_iter invalid"
          stop
         endif 

        else if (tileloop.eq.0) then
         ! do nothing
        else
         print *,"tileloop invalid"
         stop
        endif

       else if (bfact.eq.1) then
        ! do nothing
       else
        print *,"bfact invalid78"
        stop
       endif

      else if ((enable_spectral.eq.0).or. &  ! low order
               (enable_spectral.eq.3)) then  ! SEM time only
       ! do nothing
      else
       print *,"enable_spectral invalid"
       stop
      endif

      return 
      end subroutine FORT_FACE_GRADIENTS

! Prior to calling this routine:
!  a) init_gradu_tensor(...,LOCAL_CELLTENSOR_MF,LOCAL_FACETENSOR_MF)
!     i)  doit_gradu_tensor  spectral_loop==0  itensor_iter==0
!     ii) doit_gradu_tensor  spectral_loop==1  itensor_iter==0
!     ii) doit_gradu_tensor  spectral_loop==0  itensor_iter==1
!     in doit_gradu_tensor:
!         FORT_FACE_GRADIENTS, tileloop==0 (low order), tileloop==1 (SEM)
!  b) spectral_loop=0,1
!     dir=1..sdim
!     tileloop=0...3
!       FORT_CROSSTERM
!
! in CROSSTERM:
!  if tileloop==0  spectral_loop==0 
!   low order fluxes (slopecrossterm)
!  if tileloop==0  spectral_loop==1
!   do nothing 
!  if tileloop==1  spectral_loop==0
!   high order fluxes
!  if tileloop==1  spectral_loop==1
!   do nothing 
!  if tileloop==2  spectral_loop==0
!   xflux=visc_coef*xface(facevisc_index+1)*(xflux+divterm)
!  if tileloop==2  spectral_loop==1
!   do nothing 
!  if tileloop==3  spectral_loop==0
!   semflux=xflux
!  if tileloop==3  spectral_loop==1
!   xflux=(semflux_in + semflux_out)/2
!
! radial velocity is negated if r<0
! -dt * visc_coef * viscface * (grad U + grad U^T)
! fluxes found on "dir" faces
      subroutine FORT_CROSSTERM( &
       nsolve, &
       tileloop, &
       dir, &  ! dir=1..sdim
       operation_flag, & ! 8
       enable_spectral, &
       spectral_loop, &
       ncfluxreg, &
       semflux,DIMS(semflux), & 
       mask,DIMS(mask), &  ! 1=fine/fine 0=coarse/fine
       maskcoef,DIMS(maskcoef), & ! 1=not cov by level+1 or outside.
       faceLS,DIMS(faceLS), & 
       mdata,DIMS(mdata), & 
       tdata,DIMS(tdata), & 
       c_tdata,DIMS(c_tdata), & 
       maskSEM,DIMS(maskSEM), &
       xlo,dx, &
       dt, &
       cur_time, &
       vel,DIMS(vel), &
       levelpc,DIMS(levelpc), &
       xflux,DIMS(xflux), &
       xface,DIMS(xface), &
       recon,DIMS(recon), &  
       facevisc_index, &
       vofface_index, &
       massface_index, &
       ncphys, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       rzflag, &
       velbc, &
       visc_coef, &
       nmat, &
       nden, &
       ntensor, &
       constant_viscosity, &
       homflag)
      use probcommon_module
      use global_utility_module
      use godunov_module
      use MOF_routines_module
 
      IMPLICIT NONE

      INTEGER_T, intent(in) :: nsolve
      INTEGER_T, intent(in) :: ntensor
      INTEGER_T, intent(in) :: tileloop
      INTEGER_T, intent(in) :: spectral_loop
      INTEGER_T, intent(in) :: ncfluxreg
      INTEGER_T, intent(in) :: operation_flag
      INTEGER_T, intent(in) :: enable_spectral
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: facevisc_index
      INTEGER_T, intent(in) :: massface_index
      INTEGER_T, intent(in) :: vofface_index
      INTEGER_T, intent(in) :: ncphys
      INTEGER_T, intent(in) :: homflag
      INTEGER_T :: nc
      INTEGER_T, intent(in) :: constant_viscosity
      INTEGER_T, intent(in) :: nmat,nden
      INTEGER_T, intent(in) :: velbc(SDIM,2,SDIM) 
      INTEGER_T, intent(in) :: rzflag 
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T :: indexlo(SDIM),indexhi(SDIM)
      INTEGER_T :: sideidx(SDIM)
      INTEGER_T :: indexmid(SDIM)
      INTEGER_T :: index_flux(SDIM)
      INTEGER_T :: index_edge(SDIM)
      INTEGER_T :: index_opp(SDIM)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: DIMDEC(semflux)
      INTEGER_T, intent(in) :: DIMDEC(mask)
      INTEGER_T, intent(in) :: DIMDEC(maskcoef)
      INTEGER_T, intent(in) :: DIMDEC(maskSEM)
      INTEGER_T, intent(in) :: DIMDEC(faceLS)
      INTEGER_T, intent(in) :: DIMDEC(mdata)
      INTEGER_T, intent(in) :: DIMDEC(tdata)
      INTEGER_T, intent(in) :: DIMDEC(c_tdata)
      INTEGER_T, intent(in) :: DIMDEC(vel)
      INTEGER_T, intent(in) :: DIMDEC(levelpc)
      INTEGER_T, intent(in) :: DIMDEC(xflux)
      INTEGER_T, intent(in) :: DIMDEC(xface)
      INTEGER_T, intent(in) :: DIMDEC(recon)
  
      REAL_T, intent(in) :: dt 
      REAL_T, intent(in) :: cur_time
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM) 
      REAL_T, intent(in) :: mask(DIMV(mask))
      REAL_T, intent(in) :: maskcoef(DIMV(maskcoef))
      REAL_T, intent(inout) :: semflux(DIMV(semflux),ncfluxreg)

      REAL_T, intent(in) :: faceLS(DIMV(faceLS),SDIM)
      REAL_T, intent(in) :: mdata(DIMV(mdata),SDIM)
      REAL_T, intent(in) :: tdata(DIMV(tdata),ntensor)
      REAL_T, intent(in) :: c_tdata(DIMV(c_tdata),ntensor)

      REAL_T, intent(in) :: maskSEM(DIMV(maskSEM))
      REAL_T, intent(in) :: vel(DIMV(vel),SDIM)
      REAL_T, intent(in) :: levelpc(DIMV(levelpc),nmat)

      REAL_T, intent(out) :: xflux(DIMV(xflux),nsolve)  ! u
      REAL_T, intent(in) :: xface(DIMV(xface),ncphys)

      REAL_T, intent(in) :: recon(DIMV(recon),nmat*ngeom_recon)

      REAL_T, intent(in) :: visc_coef

      INTEGER_T, intent(in) :: dir
 
      INTEGER_T ilo,ihi 
      INTEGER_T jlo,jhi 
      INTEGER_T klo,khi 

      INTEGER_T i,j,k
      INTEGER_T dir2
      INTEGER_T ic,jc,kc
      INTEGER_T dirtan(2)
      INTEGER_T coupling(2)
      INTEGER_T ii,jj,kk,im1,jm1,km1
      REAL_T gradterm,alpha
      INTEGER_T side
      INTEGER_T ux,vx,wx,uy,vy,wy,uz,vz,wz,nbase
      INTEGER_T uxMM
      INTEGER_T vyMM
      INTEGER_T wzMM

      INTEGER_T im
      REAL_T LSleft(nmat)
      REAL_T LSright(nmat)
      REAL_T divterm
      INTEGER_T compressible_face
      REAL_T uxterm,vyterm,wzterm
      REAL_T visc_constant
      REAL_T diff_flux(SDIM)
      INTEGER_T imL,imR
      REAL_T total_mass,DMface
      REAL_T massfrac(nmat)
      REAL_T massF(2*nmat)
      REAL_T xsten(-1:1,SDIM)
      REAL_T xstenMAC(-1:1,SDIM)
      REAL_T xclamped_minus_sten(-1:1,SDIM)
      REAL_T xclamped_plus_sten(-1:1,SDIM)
      REAL_T xclamped_minus(SDIM)
      REAL_T xclamped_plus(SDIM)

      REAL_T LS_clamped_minus
      REAL_T LS_clamped_plus
      REAL_T vel_clamped_minus(SDIM)
      REAL_T vel_clamped_plus(SDIM)
      REAL_T vel_clamped_face(SDIM)
      REAL_T temperature_clamped_minus
      REAL_T temperature_clamped_plus
      INTEGER_T is_clamped_face
      INTEGER_T nhalf
      INTEGER_T nbr_outside_domain_flag(2)
      INTEGER_T nbr_covered_flag ! 0=covered 1=not covered
      INTEGER_T isten

      INTEGER_T local_bctype(2)
      INTEGER_T local_maskSEM
      REAL_T x_sep(2)
      REAL_T local_bcval(2)
      REAL_T local_interp(0:bfact)
      REAL_T local_vel(0:bfact)
      REAL_T RRface(0:bfact)
      REAL_T lineflux(0:bfact,SDIM)
      REAL_T local_data(1:bfact)
      REAL_T local_data_side(2)
      REAL_T local_grad(0:bfact)

      INTEGER_T maskcov
      INTEGER_T mask_out
      INTEGER_T shared_face ! in: fort_crossterm
      INTEGER_T test_maskSEM
      INTEGER_T stripstat
      INTEGER_T elemlo(3),elemhi(3)
      INTEGER_T ielem,jelem,kelem
      REAL_T avgflux(SDIM)
      INTEGER_T i_in,j_in,k_in
      INTEGER_T i_out,j_out,k_out
      INTEGER_T iflux,jflux,kflux
      INTEGER_T velcomp
      INTEGER_T tcompMM
      REAL_T xflux_temp
      INTEGER_T constant_viscosity_override
      INTEGER_T side_cell,side_face
      INTEGER_T velcomp_alt
      INTEGER_T inorm
      INTEGER_T inorm_elem
      INTEGER_T local_bc
      INTEGER_T conservative_div_uu

      REAL_T local_flux_val
      REAL_T local_flux_val_in
      REAL_T local_flux_val_out
      INTEGER_T project_option

      nhalf=1

      project_option=3

      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact invalid79"
       stop
      endif

      if (nsolve.eq.SDIM) then
       ! do nothing
      else
       print *,"nsolve invalid in CROSSTERM"
       stop
      endif

      if ((enable_spectral.lt.0).or. &
          (enable_spectral.gt.3)) then
       print *,"enable_spectral invalid crossterm"
       stop
      endif
      if (operation_flag.ne.8) then
       print *,"operation_flag invalid5"
       stop
      endif
      if (ntensor.ne.SDIM*SDIM) then
       print *,"ntensor invalid"
       stop
      endif
      if (ncfluxreg.ne.ntensor) then
       print *,"ncfluxreg invalid18 ",ncfluxreg
       stop
      endif
      if ((spectral_loop.ne.0).and. &
          (spectral_loop.ne.1)) then
       print *,"spectral_loop invalid"
       stop
      endif

      if ((homflag.eq.0).or.(homflag.eq.1)) then
       ! do nothing
      else
       print *,"homflag invalid 2"
       stop
      endif
      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt must be positive in crossterm"
       stop
      endif
      if (cur_time.ge.zero) then
       ! do nothing
      else
       print *,"cur_time must be nonneg in crossterm"
       stop
      endif
      if ((constant_viscosity.ne.0).and. &
          (constant_viscosity.ne.1)) then
       print *,"constant_viscosity invalid"
       stop
      endif

      if (visc_coef.ge.zero) then
       ! do nothing
      else
       print *,"visc_coef invalid"
       stop
      endif

       ! indexes start at 0
      if ((facevisc_index.ne.6).or. &
          (vofface_index.ne.massface_index+2*nmat)) then
       print *,"face_index bust 8"
       stop
      endif
      if (ncphys.ne.vofface_index+2*nmat) then
       print *,"ncphys invalid"
       stop
      endif

      do im=1,nmat
       if (fort_denconst(im).le.zero) then
        print *,"denconst invalid"
        stop
       endif
      enddo
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (nden.ne.nmat*num_state_material) then
       print *,"nden invalid"
       stop
      endif

      if (rzflag.eq.0) then
       ! do nothing
      else if (rzflag.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (rzflag.eq.3) then
       ! do nothing
      else
       print *,"rzflag invalid"
       stop
      endif

      if ((tileloop.eq.0).and.(spectral_loop.eq.0)) then
       if (dir.eq.1) then
        call checkbound(fablo,fabhi,DIMS(faceLS),1,-1,1277)
        call checkbound(fablo,fabhi,DIMS(mdata),1,-1,1278)
        call checkbound(fablo,fabhi,DIMS(tdata),1,-1,1279)
        call checkbound(fablo,fabhi,DIMS(c_tdata),1,-1,1265)
        call checkbound(fablo,fabhi,DIMS(vel),1,-1,1281)
        call checkbound(fablo,fabhi,DIMS(levelpc),2,-1,1284)
        call checkbound(fablo,fabhi,DIMS(recon),1,-1,234)
        call checkbound(fablo,fabhi,DIMS(maskSEM),1,-1,1264)
        call checkbound(fablo,fabhi,DIMS(semflux),1,-1,231)
        call checkbound(fablo,fabhi,DIMS(mask),1,-1,234)
        call checkbound(fablo,fabhi,DIMS(maskcoef),1,-1,234)
       endif
       call checkbound(fablo,fabhi,DIMS(xflux),0,dir-1,1285)
       call checkbound(fablo,fabhi,DIMS(xface),0,dir-1,1288)
      endif

        ! mdata(i,j,k,dir)=1 if at least one adjoining cell is a fluid cell.
        ! order: ux,vx,wx,uy,vy,wy,uz,vz,wz
      call tensorcomp_matrix(ux,uy,uz,vx,vy,vz,wx,wy,wz)

      ii=0
      jj=0
      kk=0

      if (dir.eq.1) then
       ii=1
       nbase=ux-1
      else if (dir.eq.2) then
       jj=1
       nbase=uy-1
      else if ((dir.eq.3).and.(SDIM.eq.3)) then
       kk=1
       nbase=uz-1
      else
       print *,"dir invalid crossterm"
       stop
      endif

      if (dir.eq.1) then  ! fluxes on x-face
       coupling(1)=uy
       coupling(2)=uz
       dirtan(1)=2
       dirtan(2)=SDIM
      else if (dir.eq.2) then  ! fluxes on y-face
       coupling(1)=vx
       coupling(2)=vz
       dirtan(1)=1
       dirtan(2)=SDIM
      else if ((dir.eq.3).and.(SDIM.eq.3)) then ! fluxes on z-face
       coupling(1)=wx
       coupling(2)=wy
       dirtan(1)=1
       dirtan(2)=2
      else
       print *,"dir invalid crossterm 2"
       stop
      endif

      if (tileloop.eq.0) then ! low order (grad U + grad U^T)

       if (spectral_loop.eq.0) then

        call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
          growlo,growhi,0,dir-1,34)

        do i=growlo(1),growhi(1)
        do j=growlo(2),growhi(2)
        do k=growlo(3),growhi(3)

         call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,dir-1,52)

         if (dir.eq.1) then
          inorm=i
         else if (dir.eq.2) then
          inorm=j
         else if ((dir.eq.3).and.(SDIM.eq.3)) then
          inorm=k
         else
          print *,"dir invalid"
          stop
         endif

         constant_viscosity_override=0

         side_face=0
         if (inorm.eq.fablo(dir)) then
          side_face=1
         else if (inorm.eq.fabhi(dir)+1) then
          side_face=2
         else if ((inorm.gt.fablo(dir)).and. &
                  (inorm.lt.fabhi(dir)+1)) then
          ! do nothing
         else
          print *,"inorm invalid"
          stop
         endif

         do velcomp_alt=1,SDIM
          if (side_face.eq.0) then
           ! do nothing
          else if ((side_face.eq.1).or.(side_face.eq.2)) then
           local_bc=velbc(dir,side_face,velcomp_alt)
           if ((local_bc.eq.EXT_DIR).or. &
               (local_bc.eq.REFLECT_EVEN).or. &
               (local_bc.eq.REFLECT_ODD).or. &
               (local_bc.eq.FOEXTRAP)) then
            constant_viscosity_override=1
           else if (local_bc.eq.INT_DIR) then
            ! do nothing
           else
            print *,"local_bc invalid"
            stop
           endif
          else
           print *,"side_face invalid"
           stop
          endif
         enddo ! velcomp_alt=1..sdim

         im1=i-ii
         jm1=j-jj
         km1=k-kk

         do im=1,2*nmat
          massF(im)=xface(D_DECL(i,j,k),massface_index+im)
         enddo
         do im=1,nmat
          massfrac(im)=zero
         enddo
         total_mass=zero
         do side=1,2
          do im=1,nmat
           DMface=massF(2*(im-1)+side)
           if (DMface.lt.zero) then
            print *,"DMface bust"
            stop
           endif
           total_mass=total_mass+DMface
           massfrac(im)=massfrac(im)+DMface
          enddo ! im
         enddo ! side
         if (total_mass.gt.zero) then
          do im=1,nmat
           massfrac(im)=massfrac(im)/total_mass
          enddo
         else if (total_mass.eq.zero) then
          ! do nothing
         else
          print *,"total_mass invalid"
          stop
         endif

         do velcomp=1,SDIM
          diff_flux(velcomp)=zero
         enddo  ! velcomp

         if ((constant_viscosity.eq.0).and. &
             (constant_viscosity_override.eq.0)) then

          do nc=1,SDIM-1

            ! find face stencil
           ilo=i
           ihi=i
           jlo=j
           jhi=j
           klo=k
           khi=k

           if (dir.eq.1) then ! u component  x-face
            ilo=i-1
            if (nc.eq.1) then  ! du/dy
             jhi=j+1
            else if (nc.eq.2) then  ! du/dz
             khi=k+1
            else
             print *,"nc invalid"
             stop
            endif
           else if (dir.eq.2) then  ! v component y-face
            jlo=j-1
            if (nc.eq.1) then ! dv/dx
             ihi=i+1
            else if (nc.eq.2) then ! dv/dz
             khi=k+1
            else
             print *,"nc invalid"
             stop
            endif
           else if ((dir.eq.3).and.(SDIM.eq.3)) then ! w component z-face
            klo=k-1
            if (nc.eq.1) then ! dw/dx
             ihi=i+1
            else if (nc.eq.2) then  ! dw/dy
             jhi=j+1
            else
             print *,"nc invalid"
             stop
            endif
           else
            print *,"dir invalid crossterm 3"
            stop
           endif

           ! mdata=0 if both adjoining cells to a face are solid cells or
           ! a cell pair is outside the grid.
           call slopecrossterm( &
             ntensor, &
             nmat,  &
             massfrac, &
             total_mass, &
             levelpc,DIMS(levelpc), &
             faceLS,DIMS(faceLS), &
             mdata,DIMS(mdata), &
             tdata,DIMS(tdata), &
             ii,jj,kk, &
             i,j,k,dir, &
             dirtan(nc), &
             coupling(nc), &
             ilo,ihi, &
             jlo,jhi, &
             klo,khi, &
             diff_flux(dirtan(nc)))

          enddo ! nc=1..sdim-1

         else if ((constant_viscosity.eq.1).or. &
                  (constant_viscosity_override.eq.1)) then
          ! do nothing
         else
          print *,"constant_viscosity or constant_viscosity_override invalid"
          stop
         endif
    

! 2 u_x, v_x, w_x or
! u_y, 2 v_y, w_y or
! u_z, v_z, 2 w_z 

         do nc=1,SDIM

          if (constant_viscosity.eq.0) then
           if (nc.eq.dir) then
            alpha=two
           else
            alpha=one
           endif
          else if (constant_viscosity.eq.1) then
           alpha=one
          else
           print *,"constant_viscosity invalid"
           stop
          endif

           ! ux,vx, wx  or
           ! uy,vy, wy  or
           ! uz,vz, wz  
          tcompMM=nbase+nc

          local_flux_val=tdata(D_DECL(i,j,k),tcompMM)
          if (abs(local_flux_val).le.OVERFLOW_CUTOFF) then
           ! do nothing
          else
           print *,"tdata overflow: i,j,k,nbase,nc,tdata ", &
            i,j,k,nbase,nc,local_flux_val
           print *,"dt,cur_time,level ",dt,cur_time,level
           stop
          endif

           ! use_dt=0
           ! use_HO=0
          call SEM_VISC_SANITY(101,dt,xstenMAC,nhalf,local_flux_val, &
            dir,nc,0,0,project_option,bfact,enable_spectral, &
            constant_viscosity)

          gradterm=alpha*local_flux_val

          if (side_face.eq.0) then
           ! do nothing
          else if ((side_face.eq.1).or.(side_face.eq.2)) then
           local_bc=velbc(dir,side_face,nc)
           if ((local_bc.eq.INT_DIR).or. &
               (local_bc.eq.REFLECT_ODD).or. &
               (local_bc.eq.EXT_DIR)) then
            ! do nothing
           else if ((local_bc.eq.REFLECT_EVEN).or. &
                    (local_bc.eq.FOEXTRAP)) then
            gradterm=zero
           else
            print *,"local_bc invalid"
            stop
           endif
          else
           print *,"side_face invalid"
           stop
          endif
    
          if (constant_viscosity.eq.0) then 
           diff_flux(nc)=diff_flux(nc)+gradterm
          else if (constant_viscosity.eq.1) then
           if (diff_flux(nc).eq.zero) then
            diff_flux(nc)=gradterm
           else
            print *,"diff_flux should be zero"
            stop
           endif
          else
           print *,"constant_viscosity invalid"
           stop
          endif

         enddo  ! nc=1..sdim

         do velcomp=1,SDIM
          xflux(D_DECL(i,j,k),velcomp)=diff_flux(velcomp)
         enddo  ! velcomp

        enddo
        enddo
        enddo  ! i,j,k faces

       else if (spectral_loop.eq.1) then
        ! do nothing
       else
        print *,"spectral_loop invalid"
        stop
       endif

      else if (tileloop.eq.1) then ! high order (grad U + grad U^T)

       if (spectral_loop.eq.0) then

         ! overwrite fluxes in spectral elements 
         ! CROSSTERM -> TO_MAC  (interpolate CELL_TO_MAC)
        if ((enable_spectral.eq.1).or. &  ! SEM space and time
            (enable_spectral.eq.2)) then  ! SEM space

         if (bfact.ge.2) then

          call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)
          do i=growlo(1),growhi(1)
          do j=growlo(2),growhi(2)
          do k=growlo(3),growhi(3)

           if ((dir.ge.1).and.(dir.le.SDIM)) then
            ! do nothing
           else
            print *,"dir invalid crossterm 4"
            stop
           endif

           if (dir.eq.1) then
            inorm=i
           else if (dir.eq.2) then
            inorm=j
           else if ((dir.eq.3).and.(SDIM.eq.3)) then
            inorm=k
           else
            print *,"dir invalid"
            stop
           endif

           call strip_status(i,j,k,bfact,stripstat)

           if (stripstat.eq.1) then

            local_maskSEM=NINT(maskSEM(D_DECL(i,j,k)))
            maskcov=NINT(maskcoef(D_DECL(i,j,k)))

            if ((local_maskSEM.ge.1).and. &
                (local_maskSEM.le.nmat).and. &
                (maskcov.eq.1)) then

              ! elemhi(dir)=elemlo(dir)
             call elementbox(i,j,k,bfact,dir-1,elemlo,elemhi)
             do ielem=elemlo(1),elemhi(1)
             do jelem=elemlo(2),elemhi(2)
             do kelem=elemlo(3),elemhi(3)

              if (dir.eq.1) then ! x-fluxes
               inorm_elem=ielem
              else if (dir.eq.2) then ! y-fluxes
               inorm_elem=jelem
              else if ((dir.eq.SDIM).and.(SDIM.eq.3)) then ! z-fluxes
               inorm_elem=kelem
              else
               print *,"dir invalid crossterm 5"
               stop
              endif

              if ((inorm_elem/bfact)*bfact.ne.inorm_elem) then
               print *,"inorm_elem invalid"
               stop
              endif
              if (inorm_elem.lt.0) then
               print *,"inorm_elem invalid"
               stop
              endif

              do dir2=1,SDIM
               if (fablo(dir2).lt.0) then
                print *,"fablo invalid"
                stop
               endif
              enddo ! dir2

              indexlo(1)=ielem
              indexlo(2)=jelem
              if (SDIM.eq.3) then
               indexlo(SDIM)=kelem
              endif
              do dir2=1,SDIM
               indexhi(dir2)=indexlo(dir2)
              enddo ! dir2
              indexhi(dir)=indexlo(dir)+bfact-1

              do isten=0,bfact
               do velcomp=1,SDIM
                lineflux(isten,velcomp)=zero
               enddo
              enddo

              if (constant_viscosity.eq.0) then

               do nc=1,SDIM-1

                do side=1,2

                 nbr_outside_domain_flag(side)=0

                 do dir2=1,SDIM
                  sideidx(dir2)=indexlo(dir2)
                 enddo

                 if (side.eq.1) then

                  sideidx(dir)=indexlo(dir)-1

                  do dir2=1,SDIM
                   if (sideidx(dir2).lt.fablo(dir2)) then
                    if (velbc(dir2,side,dir2).ne.INT_DIR) then
                     nbr_outside_domain_flag(side)=1
                    endif
                   endif
                  enddo ! dir2

                  i_out=indexlo(1)-ii
                  j_out=indexlo(2)-jj
                  k_out=indexlo(SDIM)-kk

                 else if (side.eq.2) then

                  sideidx(dir)=indexhi(dir)+1
                  do dir2=1,SDIM
                   if (sideidx(dir2).gt.fabhi(dir2)) then
                    if (velbc(dir2,side,dir2).ne.INT_DIR) then
                     nbr_outside_domain_flag(side)=1
                    endif
                   endif
                  enddo ! dir2

                  i_out=indexhi(1)+ii
                  j_out=indexhi(2)+jj
                  k_out=indexhi(SDIM)+kk

                 else 
                  print *,"side invalid"
                  stop
                 endif

                 if (nbr_outside_domain_flag(side).eq.0) then

                  ic=i_out
                  jc=j_out
                  kc=k_out

                   ! 0=covered 1=not covered
                  nbr_covered_flag=NINT(maskcoef(D_DECL(ic,jc,kc)))

                  call gridsten(xsten,xlo, &
                   ic,jc,kc, &
                   fablo,bfact,dx,nhalf)

                  local_bctype(side)=0  ! interior

                  local_bcval(side)=zero

                  tcompMM=coupling(nc)
  
                  if (nbr_covered_flag.eq.1) then 
                   local_data_side(side)=c_tdata(D_DECL(ic,jc,kc),tcompMM)
                  else if (nbr_covered_flag.eq.0) then
                   local_data_side(side)=zero
                   local_bctype(side)=-1 ! extrap
                   local_bcval(side)=zero
                  else
                   print *,"nbr_covered_flag invalid"
                   stop
                  endif

                 else if (nbr_outside_domain_flag(side).eq.1) then 

                  call gridsten(xsten,xlo, &
                   i_out,j_out,k_out, &
                   fablo,bfact,dx,nhalf)

                  local_data_side(side)=zero

                  if (velbc(dir,side,dir).eq.INT_DIR) then
                   print *,"velbc bust "
                   print *,"side=",side
                   print *,"nbr_outside_domain_flag= ", &
                     nbr_outside_domain_flag(side)
                   stop
                  endif
                  if ((velbc(dir,side,dir).eq.REFLECT_EVEN).or. &
                      (velbc(dir,side,dir).eq.FOEXTRAP).or. &
                      (velbc(dir,side,dir).eq.REFLECT_ODD).or. &
                      (velbc(dir,side,dir).eq.EXT_DIR)) then
                   local_bctype(side)=-1 ! extrap
                   local_bcval(side)=zero
                  else
                   print *,"velbc is corrupt"
                   stop
                  endif
                 else
                  print *,"nbr_outside_domain_flag invalid 1"
                  stop
                 endif

                enddo ! side=1..2

                do isten=0,bfact-1
                 do dir2=1,SDIM
                  indexmid(dir2)=indexlo(dir2)
                 enddo
                 indexmid(dir)=indexlo(dir)+isten
               
                 ic=indexmid(1)
                 jc=indexmid(2)
                 kc=indexmid(SDIM)

                 call gridsten(xsten,xlo, &
                  ic,jc,kc, &
                  fablo,bfact,dx,nhalf)

                 tcompMM=coupling(nc)

                 local_data(isten+1)=c_tdata(D_DECL(ic,jc,kc),tcompMM)
                enddo !isten=0..bfact-1

                do isten=0,bfact
                 indexmid(dir)=indexlo(dir)+isten
                 ic=indexmid(1)
                 jc=indexmid(2)
                 kc=indexmid(SDIM)

                 call gridstenMAC(xstenMAC,xlo, &
                  ic,jc,kc, &
                  fablo,bfact,dx,nhalf,dir-1,53)

                 RRface(isten)=xstenMAC(0,1)

                 local_vel(isten)=zero
                enddo ! isten=0..bfact

                conservative_div_uu=0

                call lineGRAD( &
                 conservative_div_uu, &
                 levelrz, &
                 dir, &
                 nc, &
                 RRface, &
                 local_bctype, &
                 local_bcval, &
                 local_vel, &
                 local_data, &
                 local_data_side, &
                 local_grad, &
                 local_interp, &
                 bfact, &
                 dx(dir), &
                 x_sep, &
                 operation_flag)

                do isten=0,bfact

                 side_cell=0
                 if (inorm+isten.eq.fablo(dir)) then
                  side_cell=1
                 else if (inorm+isten.eq.fabhi(dir)+1) then
                  side_cell=2
                 else if ((inorm+isten.gt.fablo(dir)).and. &
                          (inorm+isten.lt.fabhi(dir)+1)) then
                  ! do nothing
                 else
                  print *,"inorm or isten invalid"
                  stop
                 endif

                 constant_viscosity_override=0
                
                 do velcomp_alt=1,SDIM
                  if (side_cell.eq.0) then
                   ! do nothing
                  else if (((side_cell.eq.1).and.(isten.eq.0)).or. &
                           ((side_cell.eq.2).and.(isten.eq.bfact))) then
                   local_bc=velbc(dir,side_cell,velcomp_alt)
                   if ((local_bc.eq.EXT_DIR).or. &
                       (local_bc.eq.REFLECT_EVEN).or. &
                       (local_bc.eq.REFLECT_ODD).or. &
                       (local_bc.eq.FOEXTRAP)) then
                    constant_viscosity_override=1
                   else if (local_bc.eq.INT_DIR) then
                    ! do nothing
                   else
                    print *,"local_bc invalid"
                    stop
                   endif
                  else if ((isten.gt.0).and.(isten.lt.bfact)) then
                   ! do nothing
                  else
                   print *,"side_cell invalid1"
                   stop
                  endif
                 enddo ! velcomp_alt=1..sdim

                 if (constant_viscosity_override.eq.0) then
                  lineflux(isten,dirtan(nc))= & 
                   lineflux(isten,dirtan(nc))+local_interp(isten)
                 else if (constant_viscosity_override.eq.1) then
                  ! do nothing
                 else
                  print *,"constant_viscosity_override invalid"
                  stop
                 endif
                enddo ! isten=0..bfact

               enddo ! nc=1..sdim-1

              else if (constant_viscosity.eq.1) then
               ! do nothing
              else
               print *,"constant_viscosity invalid"
               stop
              endif
                 
              do isten=0,bfact

               side_cell=0
               if (inorm+isten.eq.fablo(dir)) then
                side_cell=1
               else if (inorm+isten.eq.fabhi(dir)+1) then
                side_cell=2
               else if ((inorm+isten.gt.fablo(dir)).and. &
                        (inorm+isten.lt.fabhi(dir)+1)) then
                ! do nothing
               else
                print *,"inorm or isten invalid"
                stop
               endif

               ! in: crossterm; prevent race condition when doing tiling.
               ! shared_face==1 => two threads would compete for same
               ! point without intervention.
               shared_face=0 

               do dir2=1,SDIM
                indexmid(dir2)=indexlo(dir2)
               enddo
               indexmid(dir)=indexlo(dir)+isten

               ic=indexmid(1)
               jc=indexmid(2)
               kc=indexmid(SDIM)

               call gridstenMAC_level(xstenMAC, &
                ic,jc,kc,level,nhalf,dir-1,54)

               if (isten.eq.bfact) then ! right side of element

                test_maskSEM=NINT(maskSEM(D_DECL(ic,jc,kc)))
                maskcov=NINT(maskcoef(D_DECL(ic,jc,kc)))

                if ((indexmid(dir).gt.fablo(dir)).and. &
                    (indexmid(dir).le.fabhi(dir)).and. &
                    (test_maskSEM.eq.local_maskSEM).and. &
                    (maskcov.eq.1)) then
                 shared_face=1
                else if ((indexmid(dir).eq.fabhi(dir)+1).or. &
                         (test_maskSEM.ne.local_maskSEM).or. &
                         (maskcov.eq.0)) then
                 ! do nothing
                else
                 print *,"indexmid,test_maskSEM, or maskcov invalid"
                 stop
                endif

               else if ((isten.ge.0).and.(isten.lt.bfact)) then
                ! do nothing
               else
                print *,"isten invalid"
                stop
               endif

               do nc=1,SDIM

                if (constant_viscosity.eq.0) then
                 if (nc.eq.dir) then
                  alpha=two
                 else
                  alpha=one
                 endif
                else if (constant_viscosity.eq.1) then
                 alpha=one
                else 
                 print *,"constant_viscosity invalid"
                 stop
                endif

                ! ux,vx, wx  or
                ! uy,vy, wy  or
                ! uz,vz, wz  
                tcompMM=nbase+nc
        
                local_flux_val=tdata(D_DECL(ic,jc,kc),tcompMM)

                 ! use_dt=0
                 ! use_HO=1
                call SEM_VISC_SANITY(1,dt,xstenMAC,nhalf,local_flux_val, &
                  dir,nc,0,1,project_option,bfact,enable_spectral, &
                  constant_viscosity)

                gradterm=alpha*local_flux_val

                if (side_cell.eq.0) then
                 ! do nothing
                else if (((side_cell.eq.1).and.(isten.eq.0)).or. &
                         ((side_cell.eq.2).and.(isten.eq.bfact))) then
                 local_bc=velbc(dir,side_cell,nc)
                 if ((local_bc.eq.INT_DIR).or. &
                     (local_bc.eq.REFLECT_ODD).or. &
                     (local_bc.eq.EXT_DIR)) then
                  ! do nothing
                 else if ((local_bc.eq.REFLECT_EVEN).or. &
                          (local_bc.eq.FOEXTRAP)) then
                  gradterm=zero
                 else
                  print *,"local_bc invalid"
                  stop
                 endif
                else if ((isten.gt.0).and.(isten.lt.bfact)) then
                 ! do nothing
                else
                 print *,"side_cell invalid2"
                 print *,"side_cell=",side_cell
                 print *,"isten,bfact= ",isten,bfact
                 print *,"dir=",dir
                 print *,"nc=",nc
                 print *,"local_bc=",local_bc
                 print *,"int dir= ",INT_DIR
                 print *,"ref odd= ",REFLECT_ODD
                 print *,"ext dir= ",EXT_DIR
                 print *,"ref evn= ",REFLECT_EVEN
                 print *,"fo ext=",FOEXTRAP
                 print *,"i,j,k=",i,j,k
                 print *,"inorm=",inorm
                 print *,"growlo=",growlo(1),growlo(2),growlo(3)
                 print *,"growhi=",growhi(1),growhi(2),growhi(3)
                 print *,"fablo(dir),fabhi(dir) ",fablo(dir),fabhi(dir)
                 print *,"ielem,jelem,kelem ",ielem,jelem,kelem
                 print *,"inorm_elem ",inorm_elem
                 stop
                endif
               
                if (constant_viscosity.eq.0) then 
                 lineflux(isten,nc)=lineflux(isten,nc)+gradterm
                else if (constant_viscosity.eq.1) then
                 if (lineflux(isten,nc).eq.zero) then
                  lineflux(isten,nc)=gradterm
                 else
                  print *,"lineflux should be zero"
                  stop
                 endif
                else
                 print *,"constant_viscosity invalid"
                 stop
                endif

               enddo  ! nc=1..sdim
            
               ! shared_face=1 if right element edge and not the right
               ! edge of the grid or next to a low order or covered
               ! element.
               if (shared_face.eq.0) then

                do velcomp=1,SDIM
                 xflux(D_DECL(ic,jc,kc),velcomp)=lineflux(isten,velcomp)
                enddo  ! velcomp

               else if (shared_face.eq.1) then
                ! do nothing
               else
                print *,"shared_face invalid"
                stop
               endif

              enddo ! isten=0..bfact

             enddo
             enddo
             enddo !ielem,jelem,kelem

            else if ((local_maskSEM.eq.0).or. &
                     (maskcov.eq.0)) then
             ! do nothing
            else
             print *,"local_maskSEM or maskcov invalid"
             stop
            endif 

           else if (stripstat.eq.0) then
            ! do nothing
           else
            print *,"stripstat invalid"
            stop
           endif

          enddo
          enddo
          enddo ! i,j,k

         else if (bfact.eq.1) then
          ! do nothing
         else
          print *,"bfact invalid80"
          stop
         endif

        else if ((enable_spectral.eq.0).or. &  ! No SEM
                 (enable_spectral.eq.3)) then  ! SEM time only.
         ! do nothing
        else
         print *,"enable_spectral invalid"
         stop
        endif

       else if (spectral_loop.eq.1) then
        ! do nothing
       else
        print *,"spectral_loop invalid"
        stop
       endif

      ! (1) fluxes+=divu  (2) fluxes*=visccoef
      else if (tileloop.eq.2) then 

       if (spectral_loop.eq.0) then 

          ! does not include the right face of the tile unless
          ! tilehi==fabhi
        call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
         growlo,growhi,0,dir-1,35)

         ! u_x+v_y+w_z on flux face
         ! multiply by visc_constant
        do i=growlo(1),growhi(1)
        do j=growlo(2),growhi(2)
        do k=growlo(3),growhi(3)

          ! dir=1..sdim
         call gridstenMAC_level(xstenMAC, &
           i,j,k,level,nhalf,dir-1,55)

         im1=i-ii
         jm1=j-jj
         km1=k-kk

         if (dir.eq.1) then
          inorm=i
         else if (dir.eq.2) then
          inorm=j
         else if ((dir.eq.3).and.(SDIM.eq.3)) then
          inorm=k
         else
          print *,"dir invalid"
          stop
         endif

         side_face=0
         if (inorm.eq.fablo(dir)) then
          side_face=1
         else if (inorm.eq.fabhi(dir)+1) then
          side_face=2
         else if ((inorm.gt.fablo(dir)).and. &
                  (inorm.lt.fabhi(dir)+1)) then
          ! do nothing
         else
          print *,"inorm invalid"
          stop
         endif

         do im=1,nmat
          LSleft(im)=levelpc(D_DECL(im1,jm1,km1),im)
          LSright(im)=levelpc(D_DECL(i,j,k),im)
         enddo
         call get_primary_material(LSleft,nmat,imL)
         call get_primary_material(LSright,nmat,imR)

         if ((imL.lt.1).or.(imL.gt.nmat)) then
          print *,"imL invalid"
          stop
         endif
         if ((imR.lt.1).or.(imR.gt.nmat)) then
          print *,"imR invalid"
          stop
         endif

         call gridsten_level(xclamped_minus_sten,im1,jm1,km1,level,nhalf)
         call gridsten_level(xclamped_plus_sten,i,j,k,level,nhalf)
         do dir2=1,SDIM
          xclamped_minus(dir2)=xclamped_minus_sten(0,dir2)
          xclamped_plus(dir2)=xclamped_plus_sten(0,dir2)
         enddo

          ! LS>0 if clamped
         call SUB_clamped_LS(xclamped_minus,cur_time,LS_clamped_minus, &
              vel_clamped_minus,temperature_clamped_minus)
         call SUB_clamped_LS(xclamped_plus,cur_time,LS_clamped_plus, &
              vel_clamped_plus,temperature_clamped_plus)
         if ((LS_clamped_minus.ge.zero).or. &
             (LS_clamped_plus.ge.zero)) then
          is_clamped_face=1
          do dir2=1,SDIM
           if (LS_clamped_minus.lt.zero) then
            vel_clamped_face(dir2)=vel_clamped_plus(dir2)
            is_clamped_face=2
           else if (LS_clamped_plus.lt.zero) then
            vel_clamped_face(dir2)=vel_clamped_minus(dir2)
            is_clamped_face=3
           else
            vel_clamped_face(dir2)=half*(vel_clamped_plus(dir2)+ &
              vel_clamped_minus(dir2))
           endif
          enddo
         else if ((LS_clamped_minus.lt.zero).and. &
                  (LS_clamped_plus.lt.zero)) then
          is_clamped_face=0
         else
          print *,"LS_clamped plus or minus is NaN"
          stop
         endif

         compressible_face=1

         if ((is_clamped_face.eq.1).or. &
             (is_clamped_face.eq.2).or. &
             (is_clamped_face.eq.3)) then
          compressible_face=0
         else if (is_clamped_face.eq.0) then

          if ((is_rigid(nmat,imL).eq.1).or. &
              (is_rigid(nmat,imR).eq.1)) then
           compressible_face=0
          else if ((is_rigid(nmat,imL).eq.0).and. &
                   (is_rigid(nmat,imR).eq.0)) then
           ! do nothing
          else
           print *,"is_rigid invalid"
           stop
          endif
         else
          print *,"is_clamped_face invalid"
          stop
         endif

         if ((fort_material_type(imL).eq.0).or. &
             (fort_material_type(imR).eq.0)) then
          compressible_face=0
         endif

         do velcomp_alt=1,SDIM
          if (side_face.eq.0) then
           ! do nothing
          else if ((side_face.eq.1).or.(side_face.eq.2)) then
           local_bc=velbc(dir,side_face,velcomp_alt)
           if ((local_bc.eq.EXT_DIR).or. &
               (local_bc.eq.REFLECT_EVEN).or. &
               (local_bc.eq.REFLECT_ODD).or. &
               (local_bc.eq.FOEXTRAP)) then
            compressible_face=0
           else if (local_bc.eq.INT_DIR) then
            ! do nothing
           else
            print *,"local_bc invalid"
            stop
           endif
          else
           print *,"side_face invalid"
           stop
          endif
         enddo ! velcomp_alt=1..sdim

         uxMM=ux
         vyMM=vy
         wzMM=wz

         wzterm=zero

         if (dir.eq.1) then

          uxterm=tdata(D_DECL(i,j,k),uxMM)
          vyterm=(tdata(D_DECL(i,j,k),vyMM)+tdata(D_DECL(i,j+1,k),vyMM)+ &
           tdata(D_DECL(i-1,j,k),vyMM)+tdata(D_DECL(i-1,j+1,k),vyMM))/four
          if (SDIM.eq.3) then
           wzterm= &
            (tdata(D_DECL(i,j,k),wzMM)+tdata(D_DECL(i,j,k+1),wzMM)+ &
             tdata(D_DECL(i-1,j,k),wzMM)+tdata(D_DECL(i-1,j,k+1),wzMM))/four
          endif

         else if (dir.eq.2) then

          uxterm=(tdata(D_DECL(i,j,k),uxMM)+tdata(D_DECL(i+1,j,k),uxMM)+ &
           tdata(D_DECL(i,j-1,k),uxMM)+tdata(D_DECL(i+1,j-1,k),uxMM))/four
          vyterm=tdata(D_DECL(i,j,k),vyMM)
          if (SDIM.eq.3) then
           wzterm= &
            (tdata(D_DECL(i,j,k),wzMM)+tdata(D_DECL(i,j,k+1),wzMM)+ &
             tdata(D_DECL(i,j-1,k),wzMM)+tdata(D_DECL(i,j-1,k+1),wzMM))/four
          endif

         else if ((dir.eq.3).and.(SDIM.eq.3)) then

          uxterm=(tdata(D_DECL(i,j,k),uxMM)+tdata(D_DECL(i+1,j,k),uxMM)+ &
           tdata(D_DECL(i,j,k-1),uxMM)+tdata(D_DECL(i+1,j,k-1),uxMM))/four
          vyterm=(tdata(D_DECL(i,j,k),vyMM)+tdata(D_DECL(i,j+1,k),vyMM)+ &
           tdata(D_DECL(i,j,k-1),vyMM)+tdata(D_DECL(i,j+1,k-1),vyMM))/four
          wzterm=tdata(D_DECL(i,j,k),wzMM)

         else
          print *,"dir invalid crossterm 6"
          stop
         endif

         if (compressible_face.eq.0) then
          divterm=zero
         else if (mdata(D_DECL(i,j,k),dir).eq.zero) then
          divterm=zero
         else if (compressible_face.eq.1) then
          divterm=-(two/SDIM)*(uxterm+vyterm+wzterm)
         else
          print *,"compressible_face bust"
          stop
         endif

         visc_constant=visc_coef*xface(D_DECL(i,j,k),facevisc_index+1)
         if (visc_constant.lt.zero) then
          print *,"visc_constant cannot be negative"
          stop
         endif
         visc_constant=-dt*visc_constant

         do velcomp=1,SDIM

          local_flux_val=xflux(D_DECL(i,j,k),velcomp)+divterm

           ! use_dt=0
           ! use_HO=1
          call SEM_VISC_SANITY(2,dt,xstenMAC,nhalf,local_flux_val, &
            dir,velcomp,0,1,project_option,bfact,enable_spectral, &
            constant_viscosity)

          xflux(D_DECL(i,j,k),velcomp)=visc_constant*local_flux_val

         enddo  ! velcomp

        enddo
        enddo
        enddo  ! i,j,k add divterm to fluxes and multiply by visc_constant

       else if (spectral_loop.eq.1) then
        ! do nothing
       else
        print *,"spectral_loop invalid"
        stop
       endif

       ! average "left" and "right" fluxes at faces separating
       ! 2 spectral elements.
      else if (tileloop.eq.3) then 

       if ((enable_spectral.eq.1).or. & ! SEM space and time
           (enable_spectral.eq.2)) then ! SEM space only

        if (bfact.ge.2) then

         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)
         if ((dir.lt.1).or.(dir.gt.SDIM)) then
          print *,"dir invalid crossterm"
          stop
         endif
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)

          call strip_status(i,j,k,bfact,stripstat)

          if (stripstat.eq.1) then

           local_maskSEM=NINT(maskSEM(D_DECL(i,j,k)))
           maskcov=NINT(maskcoef(D_DECL(i,j,k)))

           if ((local_maskSEM.ge.1).and. &
               (local_maskSEM.le.nmat).and. &
               (maskcov.eq.1)) then

            call elementbox(i,j,k,bfact,dir-1,elemlo,elemhi)
            do ielem=elemlo(1),elemhi(1)
            do jelem=elemlo(2),elemhi(2)
            do kelem=elemlo(3),elemhi(3)

             if (dir.eq.1) then ! x-fluxes
              inorm_elem=ielem
             else if (dir.eq.2) then ! y-fluxes
              inorm_elem=jelem
             else if ((dir.eq.SDIM).and.(SDIM.eq.3)) then ! z-fluxes
              inorm_elem=kelem
             else
              print *,"dir invalid crossterm 7"
              stop
             endif

             if ((inorm_elem/bfact)*bfact.ne.inorm_elem) then
              print *,"inorm_elem invalid"
              stop
             endif
             if (inorm_elem.lt.0) then
              print *,"inorm_elem invalid"
              stop
             endif

             indexlo(1)=ielem
             indexlo(2)=jelem
             if (SDIM.eq.3) then
              indexlo(SDIM)=kelem
             endif
             do dir2=1,SDIM
              index_edge(dir2)=indexlo(dir2)
              index_opp(dir2)=indexlo(dir2)
              index_flux(dir2)=indexlo(dir2)
             enddo ! dir2

             do side=1,2

              ! in: crossterm  avoid race conditions when doing tiling.
              ! shared_face==1 => two threads would compete for same
              ! point without intervention.
              shared_face=0 

              if (side.eq.1) then
               index_flux(dir)=indexlo(dir)
               index_edge(dir)=indexlo(dir)
               index_opp(dir)=index_edge(dir)-1
              else if (side.eq.2) then
               index_flux(dir)=indexlo(dir)+bfact
               index_edge(dir)=indexlo(dir)+bfact-1
               index_opp(dir)=index_edge(dir)+1
              else
               print *,"side invalid"
               stop
              endif

              i_in=index_edge(1)
              j_in=index_edge(2)
              k_in=index_edge(SDIM)

              i_out=index_opp(1)
              j_out=index_opp(2)
              k_out=index_opp(SDIM)

              iflux=index_flux(1)
              jflux=index_flux(2)
              kflux=index_flux(SDIM)

              call gridstenMAC_level(xstenMAC, &
               iflux,jflux,kflux,level,nhalf,dir-1,56)

              test_maskSEM=NINT(maskSEM(D_DECL(i_out,j_out,k_out)))
              maskcov=NINT(maskcoef(D_DECL(i_out,j_out,k_out)))

              if (side.eq.2) then

               if ((index_edge(dir).ge.fablo(dir)).and. &
                   (index_edge(dir).lt.fabhi(dir)).and. &
                   (test_maskSEM.eq.local_maskSEM).and. &
                   (maskcov.eq.1)) then
                shared_face=1
               else if ((index_edge(dir).eq.fabhi(dir)).or. &
                        (test_maskSEM.ne.local_maskSEM).or. &
                        (maskcov.eq.0)) then
                ! do nothing
               else
                print *,"index_edge invalid"
                stop
               endif

              else if (side.eq.1) then
               ! do nothing
              else
               print *,"side invalid"
               stop
              endif

               ! set mask_out=0 at coarse-fine grid border.
               ! use the (already init) flux from the fine side only.
              mask_out=1
              do dir2=1,SDIM
               if ((index_opp(dir2).lt.fablo(dir2)).or. &
                   (index_opp(dir2).gt.fabhi(dir2))) then
                mask_out=NINT(mask(D_DECL(i_out,j_out,k_out)))
               endif
              enddo

               ! do not average with flux from 
               ! an element that is covered.
              if (maskcov.eq.1) then
               ! do nothing
              else if (maskcov.eq.0) then
               mask_out=0
              else
               print *,"maskcov invalid"
               stop
              endif

               ! set mask_out=0 at high/low order interface.
               ! use the (already init) flux from the high order side only.
              if (test_maskSEM.ne.local_maskSEM) then
               mask_out=0
              endif

               ! shared_face=1 for faces on right side of elements and not
               !  touching the right side of the grid, not touching
               !  a maskSEM==0 element, and not touching a covered element.
               ! shared_face=0 for faces on the left side of elements or for
               !  the face touching the right side of the grid, touching
               !  a maskSEM!=cen_maskSEM element, or touching a covered element.
               ! no need to average fluxes at a "shared_face"
               ! since iface=0 for following element corresponds to
               ! iface=bfact of the previous.
              if (shared_face.eq.1) then
               mask_out=0
              else if (shared_face.eq.0) then
               ! do nothing
              else
               print *,"shared_face invalid"
               stop
              endif

              if (spectral_loop.eq.0) then

               do nc=1,SDIM
                xflux_temp=xflux(D_DECL(iflux,jflux,kflux),nc)
              
                 ! use_dt=1
                 ! use_HO=1
                call SEM_VISC_SANITY(3,dt,xstenMAC,nhalf,xflux_temp, &
                  dir,nc,1,1,project_option,bfact,enable_spectral, &
                  constant_viscosity)

                semflux(D_DECL(i_in,j_in,k_in),nbase+nc)=xflux_temp
               enddo ! nc

              else if (spectral_loop.eq.1) then

               if (mask_out.eq.1) then

                do nc=1,SDIM

                 local_flux_val_in=semflux(D_DECL(i_in,j_in,k_in),nbase+nc)
                 local_flux_val_out=semflux(D_DECL(i_out,j_out,k_out),nbase+nc)
                  ! use_dt=1
                  ! use_HO=1
                 call SEM_VISC_SANITY(4,dt,xstenMAC,nhalf,local_flux_val_in, &
                         dir,nc,1,1,project_option,bfact,enable_spectral, &
                         constant_viscosity)
                 call SEM_VISC_SANITY(5,dt,xstenMAC,nhalf,local_flux_val_out, &
                         dir,nc,1,1,project_option,bfact,enable_spectral, &
                         constant_viscosity)

                 avgflux(nc)=half*(local_flux_val_in+local_flux_val_out)

                enddo ! nc=1..sdim

                do nc=1,SDIM
                 xflux(D_DECL(iflux,jflux,kflux),nc)=avgflux(nc)
                enddo ! nc=1..sdim
               else if (mask_out.eq.0) then
                ! do nothing
               else
                print *,"mask_out invalid"
                stop
               endif

              else
               print *,"spectral_loop invalid"
               stop
              endif

             enddo ! side=1,2

            enddo
            enddo
            enddo !elem,jelem,kelem

           else if ((local_maskSEM.eq.0).or. &
                    (maskcov.eq.0)) then
            ! do nothing
           else
            print *,"loca_maskSEM_cen or maskcov invalid"
            stop
           endif
          else if (stripstat.eq.0) then
           ! do nothing
          else
           print *,"stripstat invalid"
           stop
          endif

         enddo
         enddo
         enddo ! i,j,k

        else if (bfact.eq.1) then
         ! do nothing
        else
         print *,"bfact invalid81"
         stop
        endif

       else if ((enable_spectral.eq.0).or. &
                (enable_spectral.eq.3)) then
        ! do nothing
       else
        print *,"enable_spectral invalid"
        stop
       endif

      else 
       print *,"tileloop invalid"
       stop
      endif

      return 
      end subroutine FORT_CROSSTERM

       ! called from NavierStokes2.cpp:
       ! void NavierStokes::MAC_GRID_ELASTIC_FORCE
      subroutine FORT_MAC_ELASTIC_FORCE( &
       im_elastic, & ! 0..nmat-1
       partid, & ! 0..num_materials_viscoelastic-1
       dir, & ! 0..sdim-1
       ncomp_visc, &
       visc_coef, &
       facevisc_index, &
       faceden_index, &
       massface_index, &
       vofface_index, &
       ncphys, &
       velbc, &
       dt, &
       cur_time, &
       xlo,dx, &
       visc,DIMS(visc), &
       mask,DIMS(mask), &  ! 1=fine/fine 0=coarse/fine
       maskcoef,DIMS(maskcoef), & ! 1=not cov by level+1 or outside.
       levelpc,DIMS(levelpc), &
       XDfab,DIMS(XDfab), &
       YDfab,DIMS(YDfab), &
       ZDfab,DIMS(ZDfab), &
       xfacefab,DIMS(xfacefab), &
       UMACNEW, &
       DIMS(UMACNEW), &
       recon,DIMS(recon), &  
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       finest_level, &
       rzflag, &
       domlo,domhi, &
       nmat, &
       nten)
      use probcommon_module
      use global_utility_module
      use godunov_module
      use MOF_routines_module
      use mass_transfer_module
 
      IMPLICIT NONE

      INTEGER_T, intent(in) :: im_elastic !0..nmat-1
      INTEGER_T, intent(in) :: partid !0..num_materials_viscoelastic-1
      INTEGER_T, intent(in) :: dir  ! MAC force component, dir=0..sdim-1
      INTEGER_T, intent(in) :: ncomp_visc
      REAL_T, intent(in) :: visc_coef
      INTEGER_T, intent(in) :: facevisc_index
      INTEGER_T, intent(in) :: faceden_index
      INTEGER_T, intent(in) :: massface_index
      INTEGER_T, intent(in) :: vofface_index
      INTEGER_T, intent(in) :: ncphys
      INTEGER_T, intent(in) :: velbc(SDIM,2,SDIM) 
      REAL_T, intent(in) :: dt 
      REAL_T, intent(in) :: cur_time
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM) 
      INTEGER_T, intent(in) :: DIMDEC(visc)
      INTEGER_T, intent(in) :: DIMDEC(mask)
      INTEGER_T, intent(in) :: DIMDEC(maskcoef)
      INTEGER_T, intent(in) :: DIMDEC(levelpc)
      INTEGER_T, intent(in) :: DIMDEC(XDfab)
      INTEGER_T, intent(in) :: DIMDEC(YDfab)
      INTEGER_T, intent(in) :: DIMDEC(ZDfab)
      INTEGER_T, intent(in) :: DIMDEC(xfacefab)
      INTEGER_T, intent(in) :: DIMDEC(UMACNEW)
      INTEGER_T, intent(in) :: DIMDEC(recon)

      REAL_T, intent(in) :: visc(DIMV(visc),ncomp_visc)
      REAL_T, intent(in) :: mask(DIMV(mask))
      REAL_T, intent(in) :: maskcoef(DIMV(maskcoef))
      REAL_T, target, intent(in) :: levelpc(DIMV(levelpc),nmat*(1+SDIM))
      REAL_T, intent(in) :: XDfab(DIMV(XDfab))
      REAL_T, intent(in) :: YDfab(DIMV(YDfab))
      REAL_T, intent(in) :: ZDfab(DIMV(ZDfab))
      REAL_T, intent(in) :: xfacefab(DIMV(xfacefab),vofface_index+2*nmat)
      REAL_T, intent(inout) :: UMACNEW(DIMV(UMACNEW))
      REAL_T, intent(in) :: recon(DIMV(recon),nmat*ngeom_recon)

      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact

      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: rzflag 
      INTEGER_T, intent(in) :: domlo(SDIM),domhi(SDIM)
      INTEGER_T, intent(in) :: nmat,nten
      INTEGER_T :: i,j,k
      INTEGER_T :: dir_flux,side_flux
      INTEGER_T :: side_comp
      REAL_T :: xstenMAC(-3:3,SDIM)
      INTEGER_T nhalf
      INTEGER_T inormal
      INTEGER_T dircomp
      INTEGER_T dir_deriv,dir_pos,dir_XD
      REAL_T, target :: xflux(SDIM)
      REAL_T xplus(SDIM)
      REAL_T xminus(SDIM)
      REAL_T XDplus(SDIM)
      REAL_T XDminus(SDIM)
      REAL_T XDcenter(SDIM)
      REAL_T gradXDtensor(SDIM,SDIM)
      REAL_T DISP_TEN(SDIM,SDIM)
      REAL_T hoop_22 ! xdisp/r
      REAL_T eps_deriv,dxmin
      INTEGER_T im_elastic_p1
      INTEGER_T im_LS

      INTEGER_T dir_local
      REAL_T, target :: cell_data_interp(nmat*(1+SDIM))
      REAL_T, target :: dx_local(SDIM)
      REAL_T, target :: xlo_local(SDIM)
      INTEGER_T, target :: fablo_local(SDIM)
      INTEGER_T, target :: fabhi_local(SDIM)

      type(interp_from_grid_parm_type) :: data_in
      type(interp_from_grid_out_parm_type) :: data_out

      REAL_T xflux_local(-1:1,SDIM,SDIM)
      REAL_T yflux_local(-1:1,SDIM,SDIM)
      REAL_T zflux_local(-1:1,SDIM,SDIM)
      REAL_T center_flux(SDIM,SDIM)
      REAL_T center_hoop_22

      REAL_T LS_at_flux_point(2,SDIM,nmat*(1+SDIM))
      INTEGER_T mask_flux_point(2,SDIM)
      REAL_T x_at_flux_point(2,SDIM,SDIM)
      REAL_T n_elastic(SDIM)
      REAL_T hx,hy,hz,rplus,rminus,rval
      REAL_T LS_clamped
      REAL_T vel_clamped(SDIM)
      REAL_T temperature_clamped
      INTEGER_T local_mask
      INTEGER_T mask_left,mask_right
      INTEGER_T mask_center(SDIM)
      REAL_T force(SDIM)
      REAL_T bodyforce
      REAL_T deninv
      REAL_T XFORCE_local
      

      im_elastic_p1=im_elastic+1

      do dir_local=1,SDIM
       dx_local(dir_local)=dx(dir_local)
       xlo_local(dir_local)=xlo(dir_local)
       fablo_local(dir_local)=fablo(dir_local)
       fabhi_local(dir_local)=fabhi(dir_local)
      enddo

      call checkbound(fablo,fabhi,DIMS(visc),1,-1,11)
      call checkbound(fablo,fabhi,DIMS(mask),1,-1,1277)
      call checkbound(fablo,fabhi,DIMS(maskcoef),1,-1,1277)
      call checkbound(fablo,fabhi,DIMS(levelpc),2,-1,1277)
      call checkbound(fablo,fabhi,DIMS(recon),2,-1,1277)
      call checkbound(fablo,fabhi,DIMS(xfacefab),0,dir,1277)
      call checkbound(fablo,fabhi,DIMS(UMACNEW),0,dir,1277)
      call checkbound(fablo,fabhi,DIMS(XDfab),1,0,1277)
      call checkbound(fablo,fabhi,DIMS(YDfab),1,1,1277)
      call checkbound(fablo,fabhi,DIMS(ZDfab),1,SDIM-1,1277)

      nhalf=3
  
      call get_dxmin(dx,bfact,dxmin)
      if (dxmin.gt.zero) then
       eps_deriv=dxmin*(1.0D-2)
      else
       print *,"dxmin invalid"
       stop
      endif

      if ((im_elastic.ge.0).and. &
          (im_elastic.lt.nmat)) then
       ! do nothing
      else
       print *,"im_elastic invalid"
       stop
      endif
      if ((partid.ge.0).and. &
          (partid.lt.num_materials_viscoelastic)) then
       ! do nothing
      else
       print *,"partid invalid"
       stop
      endif

       ! dir=0..sdim-1 
      call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,0,dir,36)

       ! traverse the "dir" MAC grid.
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       if (dir.eq.0) then
        inormal=i
       else if (dir.eq.1) then
        inormal=j
       else if ((dir.eq.2).and.(SDIM.eq.3)) then
        inormal=k
       else
        print *,"dir invalid"
        stop
       endif

        ! dir=0..sdim-1 
       call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,dir,57)

       do dircomp=1,SDIM
        xflux(dircomp)=xstenMAC(0,dircomp)
       enddo
       rval=xflux(1)

       data_out%data_interp=>cell_data_interp

       data_in%level=level
       data_in%finest_level=finest_level
       data_in%bfact=bfact ! bfact=kind of spectral element grid 
       data_in%nmat=num_materials
       data_in%dx=>dx_local
       data_in%xlo=>xlo_local
       data_in%fablo=>fablo_local
       data_in%fabhi=>fabhi_local
       data_in%ngrowfab=2

       data_in%state=>levelpc
       data_in%LS=>levelpc  ! placeholder

       data_in%ncomp=nmat*(1+SDIM)
       data_in%scomp=1
       data_in%xtarget=>xflux
       data_in%interp_foot_flag=0
       call interp_from_grid_util(data_in,data_out)

       call get_primary_material(cell_data_interp,nmat,local_mask)
       if ((local_mask.eq.im_elastic_p1).and. &
           (cell_data_interp(im_elastic_p1).gt.zero)) then
        local_mask=1
       else if ((local_mask.ge.1).and.(local_mask.le.nmat)) then
        local_mask=0
       else
        print *,"local_mask invalid"
        stop
       endif
        ! LS>0 if clamped
       call SUB_clamped_LS(xflux,cur_time,LS_clamped, &
             vel_clamped,temperature_clamped)
       if (LS_clamped.ge.zero) then
        local_mask=0
       else if (LS_clamped.lt.zero) then
        ! do nothing
       else
        print *,"LS_clamped invalid"
        stop
       endif
       if ((rzflag.eq.1).or.(rzflag.eq.3)) then
        if (dir.eq.0) then
         if (abs(xflux(dir+1)).le.VOFTOL*dx(dir+1)) then
          local_mask=0
         else if (abs(xflux(dir+1)).ge.VOFTOL*dx(dir+1)) then
          ! do nothing
         else
          print *,"xflux bust"
          stop
         endif
        else if ((dir.eq.1).or.(dir.eq.SDIM-1)) then
         ! do nothing
        else
         print *,"dir invalid"
         stop
        endif
       else if (rzflag.eq.0) then
        ! do nothing
       else
        print *,"rzflag invalid"
        stop
       endif

       XFORCE_local=zero
               
        ! im_elastic_p1 dominates the center of the MAC control volume.
       if (local_mask.eq.1) then
 
         ! im_elastic_p1=1..nmat
         ! im_elastic=0..nmat-1
        do dir_local=1,SDIM
         n_elastic(dir_local)=cell_data_interp(nmat+im_elastic*SDIM+dir_local)
        enddo
         ! normalize_vector declared in GLOBALUTIL.F90
        call normalize_vector(n_elastic)

        center_hoop_22=zero

         ! traverse all the flux face centroids associated with the
         ! "dir" MAC grid control volume (i,j,k)
        do dir_flux=0,SDIM-1
        do side_flux=0,1

         ! for a uniform grid (i.e. not spectral element grid),
         ! if  dir=0
         !  xstenMAC(0,1)=i * dx        xstenMAC(1,1)=(i+1)*dx
         !  xstenMAC(0,2)=(j+1/2) * dy  xstenMAC(1,2)=(j+3/2)*dy
         !  xstenMAC(0,3)=(k+1/2) * dz  xstenMAC(1,3)=(k+3/2)*dz
         ! center of the current (dir) MAC grid control volume:
         do dircomp=1,SDIM
          xflux(dircomp)=xstenMAC(0,dircomp)
         enddo

         if (dir_flux.eq.dir) then
          if ((inormal.eq.domlo(dir+1)).and. &
              (side_flux.eq.0).and. &
              (velbc(dir+1,side_flux+1,dir+1).ne.INT_DIR)) then 
            ! 1/2 size control vol
           xflux(dir_flux+1)=xstenMAC(0,dir_flux+1)
          else if ((inormal.eq.domhi(dir+1)+1).and. &
                   (side_flux.eq.1).and. &
                   (velbc(dir+1,side_flux+1,dir+1).ne.INT_DIR)) then
           xflux(dir_flux+1)=xstenMAC(0,dir_flux+1)
          else if ((inormal.ge.domlo(dir+1)).and. &
                   (inormal.le.domhi(dir+1)+1)) then
           if (side_flux.eq.0) then
            xflux(dir_flux+1)=xstenMAC(-1,dir_flux+1)
           else if (side_flux.eq.1) then
            xflux(dir_flux+1)=xstenMAC(1,dir_flux+1)
           else
            print *,"side_flux invalid"
            stop
           endif
          else
           print *,"inormal invalid"
           stop
          endif
         else if (dir_flux.ne.dir) then
          if (side_flux.eq.0) then
           xflux(dir_flux+1)=xstenMAC(-1,dir_flux+1)
          else if (side_flux.eq.1) then
           xflux(dir_flux+1)=xstenMAC(1,dir_flux+1)
          else
           print *,"side_flux invalid"
           stop
          endif
         else
          print *,"dir_flux or dir invalid"
          stop
         endif

          ! dir=0..sdim-1 is the force component
          ! eps_deriv=1.0D-2 * dxmin
          ! (f(x+eps_deriv)-f(x-eps_deriv))/(2 * eps_deriv)
         do dir_deriv=1,SDIM ! d/dx, d/dy, d/dz
          do dir_pos=1,SDIM
           xplus(dir_pos)=xflux(dir_pos)
           xminus(dir_pos)=xflux(dir_pos)
          enddo
          xplus(dir_deriv)=xplus(dir_deriv)+eps_deriv
          xminus(dir_deriv)=xminus(dir_deriv)-eps_deriv
           ! interpfab_XDISP declared in MASS_TRANSFER_3D.F90
          call interpfab_XDISP( &
            bfact, & ! determines positioning of Gauss Legendre nodes
            level, &
            finest_level, &
            dx, &
            xlo, &
            xplus, &
            im_elastic_p1, &!1..nmat(prescribed as a fluid in the inputs file)
            nmat, &
            partid, & ! 0..num_materials_viscoelastic-1
            fablo,fabhi, &
            XDfab,DIMS(XDfab), &
            YDfab,DIMS(YDfab), &
            ZDfab,DIMS(ZDfab), &
            recon,DIMS(recon), &
            XDplus) ! XD(xplus),YD(xplus),ZD(xplus)
                    
          call interpfab_XDISP( &
            bfact, & ! determines positioning of Gauss Legendre nodes
            level, &
            finest_level, &
            dx, &
            xlo, &
            xminus, &
            im_elastic_p1, & ! 1..nmat
            nmat, &
            partid, & ! 0..num_materials_viscoelastic-1
            fablo,fabhi, &
            XDfab,DIMS(XDfab), &
            YDfab,DIMS(YDfab), &
            ZDfab,DIMS(ZDfab), &
            recon,DIMS(recon), &
            XDminus) ! XD(xminus),YD(xminus),ZD(xminus)

          do dir_XD=1,SDIM
           gradXDtensor(dir_XD,dir_deriv)= &
              (XDplus(dir_XD)-XDminus(dir_XD))/eps_deriv
          enddo
         enddo ! dir_deriv=1..sdim

         call interpfab_XDISP( &
           bfact, & ! determines positioning of Gauss Legendre nodes
           level, &
           finest_level, &
           dx, &
           xlo, &
           xflux, & ! MAC grid face center
           im_elastic_p1, & ! 1..nmat
           nmat, &
           partid, & ! 0..num_materials_viscoelastic-1
           fablo,fabhi, &
           XDfab,DIMS(XDfab), &
           YDfab,DIMS(YDfab), &
           ZDfab,DIMS(ZDfab), &
           recon,DIMS(recon), &
           XDcenter) ! XD(xflux),YD(xflux),ZD(xflux)

          ! declared in: GLOBALUTIL.F90
         call stress_from_strain( &
          im_elastic_p1, & ! =1..nmat
          xflux, &
          dx, &
          gradXDtensor, &
          XDcenter(1), &
          XDcenter(2), &
          DISP_TEN, &  ! dir_x (displace),dir_space
          hoop_22)  ! output: "theta-theta" component xdisp/r if RZ
         
         if (side_flux.eq.0) then
          side_comp=-1
         else if (side_flux.eq.1) then 
          side_comp=1
         else
          print *,"side_flux invalid"
          stop
         endif 
         do dir_XD=1,SDIM
          do dir_deriv=1,SDIM
           if (dir_flux.eq.0) then
            xflux_local(side_comp,dir_XD,dir_deriv)=DISP_TEN(dir_XD,dir_deriv)
           else if (dir_flux.eq.1) then
            yflux_local(side_comp,dir_XD,dir_deriv)=DISP_TEN(dir_XD,dir_deriv)
           else if ((dir_flux.eq.2).and.(SDIM.eq.3)) then
            zflux_local(side_comp,dir_XD,dir_deriv)=DISP_TEN(dir_XD,dir_deriv)
           else
            print *,"dir_flux invalid"
            stop
           endif
          enddo ! dir_deriv
         enddo ! dir_XD

         data_out%data_interp=>cell_data_interp

         data_in%level=level
         data_in%finest_level=finest_level
         data_in%bfact=bfact ! bfact=kind of spectral element grid 
         data_in%nmat=num_materials
         data_in%dx=>dx_local
         data_in%xlo=>xlo_local
         data_in%fablo=>fablo_local
         data_in%fabhi=>fabhi_local
         data_in%ngrowfab=2

         data_in%state=>levelpc
         data_in%LS=>levelpc  ! placeholder

         data_in%ncomp=nmat*(1+SDIM)
         data_in%scomp=1
         data_in%xtarget=>xflux
         data_in%interp_foot_flag=0
         call interp_from_grid_util(data_in,data_out)

         do im_LS=1,nmat*(1+SDIM)
          LS_at_flux_point(side_flux+1,dir_flux+1,im_LS)= &
                  cell_data_interp(im_LS)
         enddo
         call get_primary_material(cell_data_interp,nmat,local_mask)
         if ((local_mask.eq.im_elastic_p1).and. &
             (cell_data_interp(im_elastic_p1).gt.zero)) then
          local_mask=1
         else if ((local_mask.ge.1).and.(local_mask.le.nmat)) then
          local_mask=0
         else
          print *,"local_mask invalid"
          stop
         endif
         mask_flux_point(side_flux+1,dir_flux+1)=local_mask

         do dir_local=1,SDIM
          x_at_flux_point(side_flux+1,dir_flux+1,dir_local)=xflux(dir_local)
         enddo

          ! hoop_22=xdisp/r
         center_hoop_22=center_hoop_22+hoop_22

        enddo ! side_flux=0..1
        enddo ! dir_flux=0..sdim-1

        center_hoop_22=center_hoop_22/(2*SDIM)

         ! divergence of fluxes goes here
        do dir_XD=1,SDIM
         do dir_deriv=1,SDIM
          center_flux(dir_XD,dir_deriv)=zero

          xflux_local(0,dir_XD,dir_deriv)=half*( &
            xflux_local(-1,dir_XD,dir_deriv)+ &
            xflux_local(1,dir_XD,dir_deriv))
          center_flux(dir_XD,dir_deriv)= &
            center_flux(dir_XD,dir_deriv)+ &
            xflux_local(0,dir_XD,dir_deriv)

          yflux_local(0,dir_XD,dir_deriv)=half*( &
            yflux_local(-1,dir_XD,dir_deriv)+ &
            yflux_local(1,dir_XD,dir_deriv))
          center_flux(dir_XD,dir_deriv)= &
            center_flux(dir_XD,dir_deriv)+ &
            yflux_local(0,dir_XD,dir_deriv)

          if (SDIM.eq.3) then
           zflux_local(0,dir_XD,dir_deriv)=half*( &
            zflux_local(-1,dir_XD,dir_deriv)+ &
            zflux_local(1,dir_XD,dir_deriv))
           center_flux(dir_XD,dir_deriv)= &
            center_flux(dir_XD,dir_deriv)+ &
            zflux_local(0,dir_XD,dir_deriv)
          else if (SDIM.eq.2) then
           ! do nothing
          else
           print *,"dimension bust"
           stop
          endif
          center_flux(dir_XD,dir_deriv)= &
           center_flux(dir_XD,dir_deriv)/SDIM
         enddo ! dir_deriv
        enddo ! dir_XD

         ! [n dot tau dot n] = - sigma kappa
         ! [n dot tau dot tj] = 0
        
        dir_local=1
        mask_left=mask_flux_point(1,dir_local)
        mask_right=mask_flux_point(2,dir_local)
         ! declared in: GLOBALUTIL.F90
         ! mask_center=1 if mask_left==1 or mask_right==1
         ! mask_center=0 if mask_left==0 and mask_right==0
        call project_tensor(mask_center(dir_local),n_elastic, &
          mask_left,mask_right,xflux_local,dir_local)

        dir_local=2
        mask_left=mask_flux_point(1,dir_local)
        mask_right=mask_flux_point(2,dir_local)
         ! declared in: GLOBALUTIL.F90
         ! mask_center=1 if mask_left==1 or mask_right==1
         ! mask_center=0 if mask_left==0 and mask_right==0
        call project_tensor(mask_center(dir_local),n_elastic, &
          mask_left,mask_right,yflux_local,dir_local)

        if (SDIM.eq.3) then
         dir_local=SDIM
         mask_left=mask_flux_point(1,dir_local)
         mask_right=mask_flux_point(2,dir_local)
         call project_tensor(mask_center(dir_local),n_elastic, &
            mask_left,mask_right,zflux_local,dir_local)
        endif

        dir_local=1
        hx=x_at_flux_point(2,dir_local,dir_local)- &
           x_at_flux_point(1,dir_local,dir_local)
        dir_local=2
        hy=x_at_flux_point(2,dir_local,dir_local)- &
           x_at_flux_point(1,dir_local,dir_local)
        dir_local=SDIM
        hz=x_at_flux_point(2,dir_local,dir_local)- &
           x_at_flux_point(1,dir_local,dir_local)

        if ((hx.gt.zero).and.(hy.gt.zero).and.(hz.gt.zero)) then

         rplus=x_at_flux_point(2,1,1)
         rminus=x_at_flux_point(1,1,1)

         do dir_XD=1,SDIM

          force(dir_XD)=zero

          dir_local=1
          force(dir_XD)=force(dir_XD)+ &
           mask_center(dir_local)* &
              (rplus*xflux_local(1,dir_XD,dir_local)- &
               rminus*xflux_local(-1,dir_XD,dir_local))/hx

          dir_local=2
          force(dir_XD)=force(dir_XD)+ &
           mask_center(dir_local)* &
              (yflux_local(1,dir_XD,dir_local)- &
               yflux_local(-1,dir_XD,dir_local))/hy

          if (SDIM.eq.3) then
           dir_local=SDIM
           force(dir_XD)=force(dir_XD)+ &
            mask_center(dir_local)* &
               (zflux_local(1,dir_XD,dir_local)- &
                zflux_local(-1,dir_XD,dir_local))/hz
          endif

         enddo ! dir_XD=1..sdim
                   
         if (rzflag.eq.0) then
           ! do nothing
         else if (rzflag.eq.1) then

          if (SDIM.ne.2) then
           print *,"dimension bust"
           stop
          endif
           ! -T33/r
           ! center_hoop_22=xdisp/r
          dir_XD=1
          bodyforce=-(two*center_hoop_22)/rval
          if (abs(bodyforce).lt.OVERFLOW_CUTOFF) then
           ! do nothing
          else
           print *,"bodyforce overflow bodyforce,rval:",bodyforce,rval
           stop
          endif
          force(dir_XD)=force(dir_XD)+bodyforce

         else if (rzflag.eq.3) then
          ! -T22/r
          dir_XD=1
          bodyforce=-center_flux(2,2)/rval
          force(dir_XD)=force(dir_XD)+bodyforce
         else
          print *,"rzflag invalid"
          stop
         endif 

         if (rzflag.eq.0) then
          ! do nothing
         else if (rzflag.eq.1) then
          if (SDIM.ne.2) then
           print *,"dimension bust"
           stop
          endif
          ! do nothing
         else if (rzflag.eq.3) then ! T12/r
          dir_XD=2
          bodyforce=center_flux(1,2)/rval
          force(dir_XD)=force(dir_XD)+bodyforce
         else
          print *,"rzflag invalid"
          stop
         endif 

         if (is_prescribed(nmat,im_elastic_p1).eq.1) then
          print *,"im_elastic should not be an is_prescribed material"
          stop
         else if (is_prescribed(nmat,im_elastic_p1).eq.0) then
          do dir_XD=1,SDIM
           force(dir_XD)=force(dir_XD)*dt
          enddo
         else
          print *,"is_prescribed invalid"
          stop
         endif

         do dir_XD=1,SDIM
          deninv=xfacefab(D_DECL(i,j,k),faceden_index+1)

          if (deninv.ge.zero) then 
           if (abs(force(dir_XD)).lt.OVERFLOW_CUTOFF) then
            ! do nothing
           else
            print *,"elastic overflow dir_XD,force ",dir_XD,force(dir_XD)
            print *,"i,j,k,deninv ",i,j,k,deninv
            stop
           endif
 
           if (dir_XD.eq.dir+1) then
            XFORCE_local=force(dir_XD)*deninv
            UMACNEW(D_DECL(i,j,k))=UMACNEW(D_DECL(i,j,k))+XFORCE_local
           endif

          else
           print *,"deninv invalid"
           stop
          endif
         enddo ! dir_XD = 1 ..sdim

        else
         print *,"hx,hy, or hz invalid"
         stop
        endif

       else if (local_mask.eq.0) then
        ! do nothing
       else
        print *,"local_mask invalid"
        stop
       endif

      enddo ! k
      enddo ! j
      enddo ! i
      
      return 
      end subroutine FORT_MAC_ELASTIC_FORCE


      subroutine FORT_CROSSTERM_ELASTIC( &
       ncomp_visc, &
       im_tensor, & ! 0..nmat-1
       dir, &  ! dir=1..sdim
       visc,DIMS(visc), &
       mask,DIMS(mask), &  ! 1=fine/fine 0=coarse/fine
       maskcoef,DIMS(maskcoef), & ! 1=not cov by level+1 or outside.
       faceLS,DIMS(faceLS), & 
       mdata,DIMS(mdata), & 
       tdata,DIMS(tdata), & 
       c_tdata,DIMS(c_tdata), & 
       xlo,dx, &
       dt, &
       cur_time, &
       vel,DIMS(vel), &
       levelpc,DIMS(levelpc), &
       xflux,DIMS(xflux), &
       xface,DIMS(xface), &
       recon,DIMS(recon), &  
       facevisc_index, &
       vofface_index, &
       massface_index, &
       ncphys, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       rzflag, &
       velbc, &
       visc_coef, &
       nmat, &
       nden, &
       ntensor)
      use probcommon_module
      use global_utility_module
      use godunov_module
      use MOF_routines_module
 
      IMPLICIT NONE

      INTEGER_T, intent(in) :: ncomp_visc
      INTEGER_T, intent(in) :: im_tensor ! 0..nmat-1
      INTEGER_T, intent(in) :: ntensor
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: facevisc_index
      INTEGER_T, intent(in) :: massface_index
      INTEGER_T, intent(in) :: vofface_index
      INTEGER_T, intent(in) :: ncphys
      INTEGER_T :: nc
      INTEGER_T, intent(in) :: nmat,nden
      INTEGER_T, intent(in) :: velbc(SDIM,2,SDIM) 
      INTEGER_T, intent(in) :: rzflag 
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: DIMDEC(visc)
      INTEGER_T, intent(in) :: DIMDEC(mask)
      INTEGER_T, intent(in) :: DIMDEC(maskcoef)
      INTEGER_T, intent(in) :: DIMDEC(faceLS)
      INTEGER_T, intent(in) :: DIMDEC(mdata)
      INTEGER_T, intent(in) :: DIMDEC(tdata)
      INTEGER_T, intent(in) :: DIMDEC(c_tdata)
      INTEGER_T, intent(in) :: DIMDEC(vel)
      INTEGER_T, intent(in) :: DIMDEC(levelpc)
      INTEGER_T, intent(in) :: DIMDEC(xflux)
      INTEGER_T, intent(in) :: DIMDEC(xface)
      INTEGER_T, intent(in) :: DIMDEC(recon)
  
      REAL_T, intent(in) :: dt 
      REAL_T, intent(in) :: cur_time
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM) 

      REAL_T, intent(in) :: visc(DIMV(visc),ncomp_visc)
      REAL_T, intent(in) :: mask(DIMV(mask))
      REAL_T, intent(in) :: maskcoef(DIMV(maskcoef))

      REAL_T, intent(in) :: faceLS(DIMV(faceLS),SDIM)
      REAL_T, intent(in) :: mdata(DIMV(mdata),SDIM)
      REAL_T, intent(in) :: tdata(DIMV(tdata),ntensor)
      REAL_T, intent(in) :: c_tdata(DIMV(c_tdata),ntensor)

      REAL_T, intent(in) :: vel(DIMV(vel),SDIM)
      REAL_T, intent(in) :: levelpc(DIMV(levelpc),nmat)
      REAL_T, intent(out) :: xflux(DIMV(xflux),ntensor)
      REAL_T, intent(in) :: xface(DIMV(xface),ncphys)

      REAL_T, intent(in) :: recon(DIMV(recon),nmat*ngeom_recon)

      REAL_T, intent(in) :: visc_coef

      INTEGER_T, intent(in) :: dir  ! dir=1..sdim
 
      INTEGER_T ilo,ihi 
      INTEGER_T jlo,jhi 
      INTEGER_T klo,khi 

      INTEGER_T i,j,k
      INTEGER_T space_dir
      INTEGER_T xflux_comp
      INTEGER_T velcomp
      INTEGER_T dirtan(2)
      INTEGER_T coupling(SDIM,SDIM)
      INTEGER_T ii,jj,kk,im1,jm1,km1
      INTEGER_T side
      INTEGER_T ux,vx,wx,uy,vy,wy,uz,vz,wz

      INTEGER_T im
      REAL_T diff_flux(SDIM,SDIM) ! dir_x (displace),dir_space
      REAL_T DISP_TEN(SDIM,SDIM) ! dir_x (displace),dir_space
      REAL_T total_mass,DMface
      REAL_T massfrac(nmat)
      REAL_T massF(2*nmat)
      REAL_T xstenMAC(-1:1,SDIM)
      INTEGER_T nhalf

      INTEGER_T tcompMM
      INTEGER_T constant_viscosity_override
      INTEGER_T side_face
      INTEGER_T velcomp_alt
      INTEGER_T inorm
      INTEGER_T local_bc

      INTEGER_T project_option
      INTEGER_T im_elastic
      REAL_T hoop_22 ! xdisp/r
      REAL_T xdisplace_local
      REAL_T ydisplace_local
      REAL_T visc_local

      REAL_T x_stress(SDIM)

      nhalf=1

      project_option=3

      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact invalid79"
       stop
      endif

      if ((im_tensor.ge.0).and.(im_tensor.lt.num_materials)) then
       ! do nothing
      else
       print *,"im_tensor invalid"
       stop
      endif
      im_elastic=im_tensor+1

      if (ntensor.ne.SDIM*SDIM) then
       print *,"ntensor invalid"
       stop
      endif
      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt must be positive in crossterm elastic"
       stop
      endif
      if (cur_time.ge.zero) then
       ! do nothing
      else
       print *,"cur_time must be nonneg in crossterm"
       stop
      endif

      if (visc_coef.ge.zero) then
       ! do nothing
      else
       print *,"visc_coef invalid"
       stop
      endif

       ! indexes start at 0
      if ((facevisc_index.ne.6).or. &
          (vofface_index.ne.massface_index+2*nmat)) then
       print *,"face_index bust 8"
       stop
      endif
      if (ncphys.ne.vofface_index+2*nmat) then
       print *,"ncphys invalid"
       stop
      endif

      do im=1,nmat
       if (fort_denconst(im).le.zero) then
        print *,"denconst invalid"
        stop
       endif
      enddo
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (nden.ne.nmat*num_state_material) then
       print *,"nden invalid"
       stop
      endif

      if (rzflag.eq.0) then
       ! do nothing
      else if (rzflag.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (rzflag.eq.3) then
       ! do nothing
      else
       print *,"rzflag invalid"
       stop
      endif
      if (ncomp_visc.ne.3*nmat) then
       print *,"ncomp_visc invalid"
       stop
      endif

      if (dir.eq.1) then
        call checkbound(fablo,fabhi,DIMS(visc),1,-1,11)
        call checkbound(fablo,fabhi,DIMS(faceLS),1,-1,1277)
        call checkbound(fablo,fabhi,DIMS(mdata),1,-1,1278)
        call checkbound(fablo,fabhi,DIMS(tdata),1,-1,1279)
        call checkbound(fablo,fabhi,DIMS(c_tdata),1,-1,1265)
        call checkbound(fablo,fabhi,DIMS(vel),1,-1,1281)
        call checkbound(fablo,fabhi,DIMS(levelpc),2,-1,1284)
        call checkbound(fablo,fabhi,DIMS(recon),1,-1,234)
        call checkbound(fablo,fabhi,DIMS(mask),1,-1,234)
        call checkbound(fablo,fabhi,DIMS(maskcoef),1,-1,234)
      endif
      call checkbound(fablo,fabhi,DIMS(xflux),0,dir-1,1285)
      call checkbound(fablo,fabhi,DIMS(xface),0,dir-1,1288)

        ! mdata(i,j,k,dir)=1 if at least one adjoining cell is a fluid cell.
        ! order: ux,vx,wx,uy,vy,wy,uz,vz,wz
      call tensorcomp_matrix(ux,uy,uz,vx,vy,vz,wx,wy,wz)

      ii=0
      jj=0
      kk=0

       ! coupling(velcomp,space_dir)
      velcomp=1
      coupling(velcomp,1)=ux
      coupling(velcomp,2)=uy
      if (SDIM.eq.3) then
       coupling(velcomp,SDIM)=uz
      endif
      velcomp=2
      coupling(velcomp,1)=vx
      coupling(velcomp,2)=vy
      if (SDIM.eq.3) then
       coupling(velcomp,SDIM)=vz
      endif
      if (SDIM.eq.3) then
       velcomp=SDIM
       coupling(velcomp,1)=wx
       coupling(velcomp,2)=wy
       coupling(velcomp,SDIM)=wz
      endif

      if (dir.eq.1) then
       ii=1
      else if (dir.eq.2) then
       jj=1
      else if ((dir.eq.3).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"dir invalid crossterm"
       stop
      endif

      if (dir.eq.1) then  ! fluxes on x-face
       dirtan(1)=2
       dirtan(2)=SDIM
      else if (dir.eq.2) then  ! fluxes on y-face
       dirtan(1)=1
       dirtan(2)=SDIM
      else if ((dir.eq.3).and.(SDIM.eq.3)) then ! fluxes on z-face
       dirtan(1)=1
       dirtan(2)=2
      else
       print *,"dir invalid crossterm 2"
       stop
      endif

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,0,dir-1,37)

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,dir-1,58)
       do space_dir=1,SDIM
        x_stress(space_dir)=xstenMAC(0,space_dir)
       enddo

       if (dir.eq.1) then
        inorm=i
       else if (dir.eq.2) then
        inorm=j
       else if ((dir.eq.3).and.(SDIM.eq.3)) then
        inorm=k
       else
        print *,"dir invalid"
        stop
       endif

       constant_viscosity_override=0

       side_face=0
       if (inorm.eq.fablo(dir)) then
        side_face=1
       else if (inorm.eq.fabhi(dir)+1) then
        side_face=2
       else if ((inorm.gt.fablo(dir)).and. &
                (inorm.lt.fabhi(dir)+1)) then
        ! do nothing
       else
        print *,"inorm invalid"
        stop
       endif

       do velcomp_alt=1,SDIM
        if (side_face.eq.0) then
         ! do nothing
        else if ((side_face.eq.1).or.(side_face.eq.2)) then
         local_bc=velbc(dir,side_face,velcomp_alt)
         if ((local_bc.eq.EXT_DIR).or. &
             (local_bc.eq.REFLECT_EVEN).or. &
             (local_bc.eq.REFLECT_ODD).or. &
             (local_bc.eq.FOEXTRAP)) then
          constant_viscosity_override=1
         else if (local_bc.eq.INT_DIR) then
          ! do nothing
         else
          print *,"local_bc invalid"
          stop
         endif
        else
         print *,"side_face invalid"
         stop
        endif
       enddo ! velcomp_alt=1..sdim

       im1=i-ii
       jm1=j-jj
       km1=k-kk

       do im=1,2*nmat
        massF(im)=xface(D_DECL(i,j,k),massface_index+im)
       enddo
       do im=1,nmat
        massfrac(im)=zero
       enddo
       total_mass=zero
       do side=1,2
        do im=1,nmat
         DMface=massF(2*(im-1)+side)
         if (DMface.lt.zero) then
          print *,"DMface bust"
          stop
         endif
         total_mass=total_mass+DMface
         massfrac(im)=massfrac(im)+DMface
        enddo ! im
       enddo ! side
       if (total_mass.gt.zero) then
        do im=1,nmat
         massfrac(im)=massfrac(im)/total_mass
        enddo
       else if (total_mass.eq.zero) then
        ! do nothing
       else
        print *,"total_mass invalid"
        stop
       endif

       do velcomp=1,SDIM
       do space_dir=1,SDIM
        diff_flux(velcomp,space_dir)=zero
       enddo  ! space_dir
       enddo  ! velcomp

       if (constant_viscosity_override.eq.0) then

        do velcomp=1,SDIM
         
         do nc=1,SDIM-1

          ! find face stencil
          ilo=i
          ihi=i
          jlo=j
          jhi=j
          klo=k
          khi=k

          if (dir.eq.1) then ! x-face
           ilo=i-1
           if (dirtan(nc).eq.2) then  ! d/dy
            jhi=j+1
           else if ((dirtan(nc).eq.3).and.(SDIM.eq.3)) then  ! d/dz
            khi=k+1
           else
            print *,"dirtan(nc) invalid"
            stop
           endif
          else if (dir.eq.2) then  ! y-face
           jlo=j-1
           if (dirtan(nc).eq.1) then  ! d/dx
            ihi=i+1
           else if ((dirtan(nc).eq.3).and.(SDIM.eq.3)) then  ! d/dz
            khi=k+1
           else
            print *,"dirtan(nc) invalid"
            stop
           endif
          else if ((dir.eq.3).and.(SDIM.eq.3)) then ! z-face
           klo=k-1
           if (dirtan(nc).eq.1) then  ! d/dx
            ihi=i+1
           else if (dirtan(nc).eq.2) then  ! d/dy
            jhi=j+1
           else
            print *,"dirtan(nc) invalid"
            stop
           endif
          else
           print *,"dir invalid crossterm_elastic 3"
           stop
          endif

          ! mdata=0 if both adjoining cells to a face are solid cells or
          ! a cell pair is outside the grid.
          call slopecrossterm( &
           ntensor, &
           nmat,  &
           massfrac, &
           total_mass, &
           levelpc,DIMS(levelpc), &
           faceLS,DIMS(faceLS), &
           mdata,DIMS(mdata), &
           tdata,DIMS(tdata), &
           ii,jj,kk, &
           i,j,k,dir, &
           dirtan(nc), &
           coupling(velcomp,dirtan(nc)), &
           ilo,ihi, &
           jlo,jhi, &
           klo,khi, &
           diff_flux(velcomp,dirtan(nc)))

         enddo ! nc=1..sdim-1

        enddo ! velcomp=1..sdim

       else if (constant_viscosity_override.eq.1) then
        ! do nothing
       else
        print *,"constant_viscosity_override invalid"
        stop
       endif
    
       do velcomp=1,SDIM
        tcompMM=coupling(velcomp,dir)
        diff_flux(velcomp,dir)=tdata(D_DECL(i,j,k),tcompMM)

        if (side_face.eq.0) then
         ! do nothing
        else if ((side_face.eq.1).or.(side_face.eq.2)) then
         local_bc=velbc(dir,side_face,velcomp)
         if ((local_bc.eq.INT_DIR).or. &
             (local_bc.eq.REFLECT_ODD).or. &
             (local_bc.eq.EXT_DIR)) then
          ! do nothing
         else if ((local_bc.eq.REFLECT_EVEN).or. &
                  (local_bc.eq.FOEXTRAP)) then
          diff_flux(velcomp,dir)=zero
         else
          print *,"local_bc invalid"
          stop
         endif
        else
         print *,"side_face invalid"
         stop
        endif

       enddo ! velcomp=1..sdim

       xdisplace_local=half*(vel(D_DECL(i,j,k),1)+ &
         vel(D_DECL(im1,jm1,km1),1))
       ydisplace_local=half*(vel(D_DECL(i,j,k),2)+ &
         vel(D_DECL(im1,jm1,km1),2))

        ! visc(D_DECL(i,j,k),nmat+im_elastic)=
        !  elastic_viscosity * visc_coef (DERVISC)
       visc_local=half*(visc(D_DECL(i,j,k),nmat+im_elastic)+ &
         visc(D_DECL(im1,jm1,km1),nmat+im_elastic))

       if (visc_local.ge.zero) then
        ! do nothing
       else
        print *,"visc_local invalid"
        stop
       endif
        

        ! declared in GLOBALUTIL.F90
       call stress_from_strain( &
         im_elastic, & ! 1<=im_elastic<=nmat
         x_stress, &
         dx, &
         diff_flux, &
         xdisplace_local, &
         ydisplace_local, &
         DISP_TEN, & ! velcomp,space_dir
         hoop_22) ! in RZ, DISP_TEN does not have theta,theta component.
                  ! in RZ,
                  ! TNEWfab(3,3) gets 2 * hoop_22 in 
                  ! subroutine local_tensor_from_xdisplace
                  ! hoop_22=xdisp/r

       xflux_comp=1
       do velcomp=1,SDIM
       do space_dir=1,SDIM
        xflux(D_DECL(i,j,k),xflux_comp)= &
          visc_local*DISP_TEN(velcomp,space_dir)
        xflux_comp=xflux_comp+1
       enddo  ! space_dir
       enddo  ! velcomp
       if (xflux_comp-1.eq.SDIM*SDIM) then
        ! do nothing
       else
        print *,"xflux_comp invalid"
        stop
       endif

      enddo
      enddo
      enddo  ! i,j,k faces

      return 
      end subroutine FORT_CROSSTERM_ELASTIC




       ! called from: NavierStokes3.cpp
      subroutine FORT_HEATADVANCE( &
       level, &
       finest_level, &
       cur_time, &
       nmat, &
       nparts, &
       nparts_def, &
       im_solid_map, &
       nsolve, &
       nstate, &
       xlo,dx, &
       solxfab,DIMS(solxfab), &
       solyfab,DIMS(solyfab), &
       solzfab,DIMS(solzfab), &
       snew,DIMS(snew), &
       lsnew,DIMS(lsnew), &
       du,DIMS(du), &
       tilelo,tilehi, &
       fablo,fabhi,bfact)
      use probcommon_module
      use global_utility_module
      use probf90_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      REAL_T, intent(in) :: cur_time
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nparts
      INTEGER_T, intent(in) :: nparts_def
      INTEGER_T, intent(in) :: im_solid_map(nparts_def)
      INTEGER_T, intent(in) :: nsolve
      INTEGER_T, intent(in) :: nstate
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T :: growlo(3),growhi(3)

      INTEGER_T, intent(in) :: DIMDEC(solxfab)
      INTEGER_T, intent(in) :: DIMDEC(solyfab)
      INTEGER_T, intent(in) :: DIMDEC(solzfab)
      INTEGER_T, intent(in) :: DIMDEC(snew)
      INTEGER_T, intent(in) :: DIMDEC(lsnew)
      INTEGER_T, intent(in) :: DIMDEC(du)
      INTEGER_T :: i,j,k
      INTEGER_T :: dir
      INTEGER_T :: im
      INTEGER_T :: ibase
      INTEGER_T :: velcomp

      REAL_T, intent(in) :: solxfab(DIMV(solxfab),nparts_def*SDIM)
      REAL_T, intent(in) :: solyfab(DIMV(solyfab),nparts_def*SDIM)
      REAL_T, intent(in) :: solzfab(DIMV(solzfab),nparts_def*SDIM)
      REAL_T, intent(inout) :: snew(DIMV(snew),nstate)
      REAL_T, intent(in) :: lsnew(DIMV(lsnew),nmat*(SDIM+1))
      REAL_T, intent(in) :: du(DIMV(du),nsolve)

      REAL_T Tforce,new_temperature,TEMPERATURE
      REAL_T xclamped(SDIM)
      REAL_T LS_clamped
      REAL_T vel_clamped(SDIM)
      REAL_T temperature_clamped
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf

      nhalf=3

      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid veladvance"
       stop
      endif
      if (finest_level.ne.fort_finest_level) then
       print *,"finest_level invalid veladvance"
       stop
      endif
      if (nsolve.ne.1) then
       print *,"nsolve invalid"
       stop
      endif

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (nstate.ne.(SDIM+1)+ &
          nmat*(num_state_material+ngeom_raw)+1) then
       print *,"nstate invalid"
       stop
      endif

      if ((nparts.lt.0).or.(nparts.gt.nmat)) then
       print *,"nparts invalid FORT_HEATADVANCE"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.nmat)) then
       print *,"nparts_def invalid FORT_HEATADVANCE"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(solxfab),0,0,1251)
      call checkbound(fablo,fabhi,DIMS(solyfab),0,1,1251)
      call checkbound(fablo,fabhi,DIMS(solzfab),0,SDIM-1,1251)
      call checkbound(fablo,fabhi,DIMS(snew),1,-1,1251)
      call checkbound(fablo,fabhi,DIMS(lsnew),1,-1,1251)
      call checkbound(fablo,fabhi,DIMS(du),0,-1,1251)
     
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridsten_level(xsten,i,j,k,level,nhalf)
       do dir=1,SDIM
        xclamped(dir)=xsten(0,dir)
       enddo

        ! LS>0 if clamped
       call SUB_clamped_LS(xclamped,cur_time,LS_clamped, &
             vel_clamped,temperature_clamped)

       do im=1,nmat
        ibase=(SDIM+1)+(im-1)*num_state_material
        TEMPERATURE=snew(D_DECL(i,j,k),ibase+2)
        if (TEMPERATURE.le.zero) then
          print *,"HEATADVANCE: temperature must be positive"
          print *,"i,j,k,im ",i,j,k,im
          print *,"TEMPERATURE= ",TEMPERATURE
          stop
        endif
         ! viscous heating term.
        velcomp=1
        Tforce=du(D_DECL(i,j,k),velcomp)
        new_temperature=TEMPERATURE+Tforce
        if (new_temperature.le.zero) then
         new_temperature=TEMPERATURE
        endif
        snew(D_DECL(i,j,k),ibase+2)=new_temperature
       enddo ! im = 1..nmat

      enddo
      enddo
      enddo
  
      return 
      end subroutine FORT_HEATADVANCE

      module FSI_PC_module

       use iso_c_binding
       use amrex_fort_module, only : amrex_real,amrex_particle_real
       use iso_c_binding, only: c_int

       implicit none

       type, bind(C) :: particle_t
         real(amrex_particle_real) :: pos(SDIM)
           ! xfoot,dist,vel,den,T,insert time
         real(amrex_particle_real) :: extra_state(N_EXTRA_REAL)
         integer(c_int) :: id
         integer(c_int) :: cpu
       end type particle_t

       type accum_parm_type
        INTEGER_T, pointer :: fablo(:)
        INTEGER_T, pointer :: fabhi(:)
        INTEGER_T, pointer :: tilelo(:)
        INTEGER_T, pointer :: tilehi(:)
        INTEGER_T :: bfact
        INTEGER_T :: level
        INTEGER_T :: finest_level
        INTEGER_T :: matrix_points
        INTEGER_T :: RHS_points
        INTEGER_T :: ncomp_accumulate
        REAL_T, pointer :: dx(:)
        REAL_T, pointer :: xlo(:)
        INTEGER_T :: Npart
        type(particle_t), pointer, dimension(:) :: particles
        INTEGER_T :: nmat
       end type accum_parm_type

      contains

      subroutine traverse_particlesVEL( &
       accum_PARM, &
       matrixfab, &
       DIMS(matrixfab), &
       LS, &
       DIMS(LS), &
       VEL_fab, &
       DIMS(VEL_fab), &
       ncomp_accumulate)

      use probcommon_module
      use global_utility_module

      INTEGER_T, intent(in) :: ncomp_accumulate
      type(accum_parm_type), intent(in) :: accum_PARM
      INTEGER_T, intent(in) :: DIMDEC(LS) 
      INTEGER_T, intent(in) :: DIMDEC(VEL_fab) 
      INTEGER_T, intent(in) :: DIMDEC(matrixfab) 
      REAL_T, target, intent(in) :: LS( &
        DIMV(LS), &
        num_materials*(1+SDIM))
      REAL_T, target, intent(in) :: VEL_fab( &
        DIMV(VEL_fab), &
        SDIM)
      REAL_T, intent(inout) :: matrixfab( &
        DIMV(matrixfab), &
        ncomp_accumulate)

      INTEGER_T :: nhalf
      REAL_T :: eps
      INTEGER_T :: interior_ID
      INTEGER_T :: dir
      REAL_T, target :: xpart(SDIM)
      REAL_T xpartfoot(SDIM)
      REAL_T xdisp(SDIM)
      REAL_T velpart(SDIM)
      INTEGER_T cell_index(SDIM)
      INTEGER_T interior_ok
      INTEGER_T i,j,k
      REAL_T xsten(-3:3,SDIM)
      REAL_T tmp,w_p
      REAL_T xc(SDIM)
      INTEGER_T npart_local

      REAL_T, target :: cell_data_interp(SDIM)
      REAL_T, target :: dx_local(SDIM)
      REAL_T, target :: xlo_local(SDIM)
      INTEGER_T, target :: fablo_local(SDIM)
      INTEGER_T, target :: fabhi_local(SDIM)
      INTEGER_T, target :: tilelo_local(SDIM)
      INTEGER_T, target :: tilehi_local(SDIM)

      type(interp_from_grid_parm_type) :: data_in
      type(interp_from_grid_out_parm_type) :: data_out

      do dir=1,SDIM
       dx_local(dir)=accum_PARM%dx(dir)
       xlo_local(dir)=accum_PARM%xlo(dir)
       fablo_local(dir)=accum_PARM%fablo(dir)
       fabhi_local(dir)=accum_PARM%fabhi(dir)
       tilelo_local(dir)=accum_PARM%tilelo(dir)
       tilehi_local(dir)=accum_PARM%tilehi(dir)
      enddo

      call checkbound(fablo_local,fabhi_local,DIMS(LS),2,-1,1271)
      call checkbound(tilelo_local,tilehi_local,DIMS(matrixfab),1,-1,1271)
      call checkbound(fablo_local,fabhi_local,DIMS(VEL_fab),2,-1,1271)

      data_out%data_interp=>cell_data_interp

      data_in%level=accum_PARM%level
      data_in%finest_level=accum_PARM%finest_level
      data_in%bfact=accum_PARM%bfact
      data_in%nmat=accum_PARM%nmat

      data_in%dx=>dx_local
      data_in%xlo=>xlo_local
      data_in%fablo=>fablo_local
      data_in%fabhi=>fabhi_local
      data_in%ngrowfab=2

      data_in%state=>VEL_fab
      data_in%LS=>LS  ! not used

      data_in%ncomp=SDIM
      data_in%scomp=1

      nhalf=3

      eps=accum_PARM%dx(1)/10.0d0
      if (eps.gt.zero) then
       ! do nothing
      else
       print *,"eps invalid"
       stop
      endif

      if (accum_PARM%Npart.ge.0) then
       npart_local=accum_PARM%Npart
      else
       print *,"accum_PARM%Npart invalid"
       stop
      endif

      do interior_ID=1,npart_local

       if (accum_PARM%Npart.ge.0) then
        do dir=1,SDIM
         xpart(dir)=accum_PARM%particles(interior_ID)%pos(dir)
         xpartfoot(dir)=accum_PARM%particles(interior_ID)%extra_state(dir)
         xdisp(dir)=xpart(dir)-xpartfoot(dir)
         velpart(dir)= &
          accum_PARM%particles(interior_ID)%extra_state(SDIM+1+dir)
        enddo ! dir=1..sdim

        data_in%xtarget=>xpart
        data_in%interp_foot_flag=0
        call interp_from_grid_util(data_in,data_out)

        call containing_cell(accum_PARM%bfact, &
          accum_PARM%dx, &
          accum_PARM%xlo, &
          accum_PARM%fablo, &
          xpart, &
          cell_index)

        interior_ok=1
        do dir=1,SDIM
         if ((cell_index(dir).lt.accum_PARM%tilelo(dir)-1).or. &
             (cell_index(dir).gt.accum_PARM%tilehi(dir)+1)) then
          interior_ok=0
         endif
        enddo

        if (interior_ok.eq.1) then
         i=cell_index(1)
         j=cell_index(2)
         k=cell_index(SDIM)
         call gridsten_level(xsten,i,j,k,accum_PARM%level,nhalf)
         tmp=0.0d0
         do dir=1,SDIM
          xc(dir)=xsten(0,dir)
          tmp=tmp+(xpart(dir)-xc(dir))**2
         enddo
         tmp=sqrt(tmp)
         w_p=(1.0d0/(eps+tmp))

         if (w_p.gt.zero) then
          matrixfab(D_DECL(i,j,k),1)= &
           matrixfab(D_DECL(i,j,k),1)+w_p
          do dir=1,SDIM
           matrixfab(D_DECL(i,j,k),1+dir)= &
            matrixfab(D_DECL(i,j,k),1+dir)+ &
            w_p*(data_out%data_interp(dir)-velpart(dir))
          enddo
         else
          print *,"w_p invalid"
          stop
         endif
        else if (interior_ok.eq.0) then
         ! do nothing
        else
         print *,"interior_ok invalid"
         stop
        endif

       else
        print *,"accum_PARM%Npart invalid"
        stop
       endif

      enddo ! do interior_ID=1,accum_PARM%Npart

      return
      end subroutine traverse_particlesVEL

       ! called from NavierStokes.cpp:
       ! NavierStokes::assimilate_vel_from_particles()
      subroutine fort_assimilate_VEL_from_particles( &
        fluid_relaxation_time_to_particle, &
        dt, &
        tid, &  ! thread id
        tilelo,tilehi, &  ! tile box dimensions
        fablo,fabhi, &    ! fortran array box dimensions containing the tile
        bfact, &          ! space order
        level, &          ! 0<=level<=finest_level
        finest_level, &
        xlo,dx, &         ! xlo is lower left hand corner coordinate of fab
        particles_no_nbr, & 
        particles_nbr, & 
        particles_only_nbr, & 
        Np_no_nbr, & ! pass by value
        Np_nbr, & ! pass by value
        Nn, & ! pass by value
        matrix_points, & 
        RHS_points, &    
        ncomp_accumulate, & ! matrix_points+sdim * RHS_points
        nmat, &
        LS, &
        DIMS(LS), &
        SNEWfab, &  
        DIMS(SNEWfab), &
        UMACNEW, &  
        DIMS(UMACNEW), &
        VMACNEW, &  
        DIMS(VMACNEW), &
        WMACNEW, &  
        DIMS(WMACNEW), &
        vel_fab, &      
        DIMS(vel_fab), &
        matrixfab, &     ! accumulation FAB
        DIMS(matrixfab)) &
      bind(c,name='fort_assimilate_VEL_from_particles')

      use global_utility_module
      use probcommon_module
      use godunov_module
      implicit none

      INTEGER_T, intent(in) :: nmat
      REAL_T, intent(in) :: dt
      REAL_T, intent(in) :: fluid_relaxation_time_to_particle
      INTEGER_T, intent(in) :: matrix_points
      INTEGER_T, intent(in) :: RHS_points
      INTEGER_T, intent(in) :: ncomp_accumulate
      INTEGER_T, value, intent(in) :: Np_no_nbr,Np_nbr,Nn ! pass by value
      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in), target :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in), target :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      REAL_T, intent(in), target :: xlo(SDIM)
      REAL_T, intent(in), target :: dx(SDIM)
      INTEGER_T, intent(in) :: DIMDEC(LS) 
      INTEGER_T, intent(in) :: DIMDEC(matrixfab) 
      INTEGER_T, intent(in) :: DIMDEC(SNEWfab) 
      INTEGER_T, intent(in) :: DIMDEC(UMACNEW) 
      INTEGER_T, intent(in) :: DIMDEC(VMACNEW) 
      INTEGER_T, intent(in) :: DIMDEC(WMACNEW) 
      INTEGER_T, intent(in) :: DIMDEC(vel_fab) 
      REAL_T, intent(inout) :: matrixfab( &
        DIMV(matrixfab), &
        ncomp_accumulate)
      REAL_T, intent(in) :: LS( &  
        DIMV(LS), &
        nmat*(1+SDIM))
      REAL_T, intent(inout) :: SNEWfab( & 
        DIMV(SNEWfab), &
        SDIM)
      REAL_T, intent(inout) :: UMACNEW( & 
        DIMV(UMACNEW))
      REAL_T, intent(inout) :: VMACNEW( & 
        DIMV(VMACNEW))
      REAL_T, intent(inout) :: WMACNEW( & 
        DIMV(WMACNEW))
      REAL_T, intent(in) :: vel_fab( &
        DIMV(vel_fab), &
        SDIM)

      type(particle_t), intent(in), target :: particles_no_nbr(Np_no_nbr)
      type(particle_t), intent(in) :: particles_nbr(Np_nbr)
      type(particle_t), intent(in), target :: particles_only_nbr(Nn)

      type(accum_parm_type) :: accum_PARM

      INTEGER_T growlo(3)
      INTEGER_T growhi(3)
      INTEGER_T i,j,k
      INTEGER_T dir
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf

      REAL_T A_matrix,B_matrix
      REAL_T lambda
      REAL_T vel_local
      REAL_T local_wt
      INTEGER_T interior_ID
      REAL_T xpart1,xpart2
      INTEGER_T dirmac
      INTEGER_T ii,jj,kk

      nhalf=3

      if (nmat.eq.num_materials) then
       ! do nothing
      else
       print *,"nmat invalid"
       stop
      endif

      if (matrix_points.eq.1) then
       ! do nothing
      else
       print *,"matrix_points invalid"
       stop
      endif
      if (RHS_points.eq.1) then
       ! do nothing
      else
       print *,"RHS_points invalid"
       stop
      endif
      if (ncomp_accumulate.eq.matrix_points+SDIM*RHS_points) then
       ! do nothing
      else
       print *,"ncomp_accumulate invalid"
       stop
      endif
      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(LS),2,-1,1271)
      call checkbound(tilelo,tilehi,DIMS(matrixfab),1,-1,1271)
      call checkbound(fablo,fabhi,DIMS(SNEWfab),1,-1,1271)
      call checkbound(fablo,fabhi,DIMS(UMACNEW),0,0,1271)
      call checkbound(fablo,fabhi,DIMS(VMACNEW),0,1,1271)
      call checkbound(fablo,fabhi,DIMS(WMACNEW),0,SDIM-1,1271)
      call checkbound(fablo,fabhi,DIMS(vel_fab),2,-1,1271)

      accum_PARM%fablo=>fablo 
      accum_PARM%fabhi=>fabhi
      accum_PARM%tilelo=>tilelo 
      accum_PARM%tilehi=>tilehi
      accum_PARM%bfact=bfact
      accum_PARM%level=level
      accum_PARM%finest_level=finest_level
      accum_PARM%matrix_points=matrix_points
      accum_PARM%RHS_points=RHS_points
      accum_PARM%ncomp_accumulate=ncomp_accumulate
      accum_PARM%dx=>dx
      accum_PARM%xlo=>xlo
      accum_PARM%nmat=nmat

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      if (Np_no_nbr+Nn.eq.Np_nbr) then
       ! do nothing
      else
       print *,"Np_no_nbr+Nn.ne.Np_nbr"
       stop
      endif

       ! sanity check
      do interior_ID=1,Np_no_nbr
       do dir=1,SDIM
        xpart1=particles_nbr(interior_ID)%pos(dir)
        xpart2=particles_no_nbr(interior_ID)%pos(dir)
        if (abs(xpart1-xpart2).le.dx(dir)*VOFTOL) then
         ! do nothing
        else
         print *,"xpart1 or xpart2 invalid"
         stop
        endif
       enddo
      enddo ! interior_ID=1,Np_no_nbr

      accum_PARM%particles=>particles_no_nbr
      accum_PARM%Npart=Np_no_nbr

      call traverse_particlesVEL( &
        accum_PARM, &
        matrixfab, &
        DIMS(matrixfab), &
        LS, &
        DIMS(LS), &
        vel_fab, &
        DIMS(vel_fab), &
        ncomp_accumulate)

      accum_PARM%particles=>particles_only_nbr
      accum_PARM%Npart=Nn

      call traverse_particlesVEL( &
        accum_PARM, &
        matrixfab, &
        DIMS(matrixfab), &
        LS, &
        DIMS(LS), &
        vel_fab, &
        DIMS(vel_fab), &
        ncomp_accumulate)

      if (dt.gt.zero) then
       if (fluid_relaxation_time_to_particle.eq.zero) then
        local_wt=one
       else if (fluid_relaxation_time_to_particle.gt.zero) then
        local_wt=one-exp(-dt/fluid_relaxation_time_to_particle)
       else
        print *,"fluid_relaxation_time_to_particle invalid"
        stop
       endif
      else
       print *,"dt invalid"
       stop
      endif

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
       call gridsten_level(xsten,i,j,k,level,nhalf)
       A_matrix=matrixfab(D_DECL(i,j,k),1) ! sum w(xp)
       do dir=1,SDIM
        B_matrix=matrixfab(D_DECL(i,j,k),1+dir) !sum w*(vel(xp)-vel_p)
        vel_local=vel_fab(D_DECL(i,j,k),dir)
        if (A_matrix.eq.zero) then
         SNEWfab(D_DECL(i,j,k),dir)=vel_local
        else if (A_matrix.gt.zero) then
         ! lambda=sum (interp(vel)-vel_p)w_p/sum w_p
         lambda=B_matrix/A_matrix
         if ((local_wt.ge.zero).and.(local_wt.le.one)) then
          SNEWFAB(D_DECL(i,j,k),dir)= &
            vel_local-local_wt*lambda
         else
          print *,"local_wt invalid"
          stop
         endif
        else
         print *,"A_matrix invalid"
         stop
        endif
       enddo ! dir=1..SDIM
      enddo
      enddo
      enddo

      do dirmac=1,SDIM
       call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
          growlo,growhi,0,dirmac-1,38)

       ii=0
       jj=0
       kk=0
       if (dirmac.eq.1) then
        ii=1
       else if (dirmac.eq.2) then
        jj=1
       else if ((dirmac.eq.3).and.(SDIM.eq.3)) then
        kk=1
       else
        print *,"dirmac invalid"
        stop
       endif

       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
         !dirmac=1..sdim
        call gridstenMAC_level(xsten,i,j,k,level,nhalf,dirmac-1,59)

        A_matrix=matrixfab(D_DECL(i,j,k),1)+ &
                 matrixfab(D_DECL(i-ii,j-jj,k-kk),1) ! sum w(xp)
        B_matrix=matrixfab(D_DECL(i,j,k),1+dirmac)+ &
           matrixfab(D_DECL(i-ii,j-jj,k-kk),1+dirmac) !sum w*(vel(xp)-vel_p)
        if (A_matrix.eq.zero) then
         ! do nothing
        else if (A_matrix.gt.zero) then
         ! lambda=sum (interp(vel)-vel_p)w_p/sum w_p
         if ((local_wt.ge.zero).and.(local_wt.le.one)) then
          lambda=local_wt*B_matrix/A_matrix
         else
          print *,"local_wt invalid"
          stop
         endif
         if (dirmac.eq.1) then
          UMACNEW(D_DECL(i,j,k))=UMACNEW(D_DECL(i,j,k))-lambda
         else if (dirmac.eq.2) then
          VMACNEW(D_DECL(i,j,k))=VMACNEW(D_DECL(i,j,k))-lambda
         else if ((dirmac.eq.3).and.(SDIM.eq.3)) then
          WMACNEW(D_DECL(i,j,k))=WMACNEW(D_DECL(i,j,k))-lambda
         else
          print *,"dirmac invalid"
          stop
         endif
        else
         print *,"A_matrix invalid"
         stop
        endif
       enddo
       enddo
       enddo
      enddo ! dirmac=1..sdim

      end subroutine fort_assimilate_VEL_from_particles



      subroutine traverse_particles( &
       accum_PARM, &
       matrixfab, &
       DIMS(matrixfab), &
       LS, &
       DIMS(LS), &
       XDISP_fab, &
       DIMS(XDISP_fab), &
       ncomp_accumulate)

      use probcommon_module
      use global_utility_module

      INTEGER_T, intent(in) :: ncomp_accumulate
      type(accum_parm_type), intent(in) :: accum_PARM
      INTEGER_T, intent(in) :: DIMDEC(LS) 
      INTEGER_T, intent(in) :: DIMDEC(XDISP_fab) 
      INTEGER_T, intent(in) :: DIMDEC(matrixfab) 
      REAL_T, target, intent(in) :: LS( &
        DIMV(LS), &
        num_materials*(1+SDIM))
      REAL_T, target, intent(in) :: XDISP_fab( &
        DIMV(XDISP_fab), &
        SDIM)
      REAL_T, intent(inout) :: matrixfab( &
        DIMV(matrixfab), &
        ncomp_accumulate)

      INTEGER_T :: nhalf
      REAL_T :: eps
      INTEGER_T :: interior_ID
      INTEGER_T :: dir
      REAL_T, target :: xpart(SDIM)
      REAL_T xpartfoot(SDIM)
      REAL_T xdisp(SDIM)
      INTEGER_T cell_index(SDIM)
      INTEGER_T interior_ok
      INTEGER_T i,j,k
      REAL_T xsten(-3:3,SDIM)
      REAL_T tmp,w_p
      REAL_T xc(SDIM)
      INTEGER_T npart_local

      REAL_T, target :: cell_data_interp(SDIM)
      REAL_T, target :: dx_local(SDIM)
      REAL_T, target :: xlo_local(SDIM)
      INTEGER_T, target :: fablo_local(SDIM)
      INTEGER_T, target :: fabhi_local(SDIM)

      type(interp_from_grid_parm_type) :: data_in
      type(interp_from_grid_out_parm_type) :: data_out

      do dir=1,SDIM
       dx_local(dir)=accum_PARM%dx(dir)
       xlo_local(dir)=accum_PARM%xlo(dir)
       fablo_local(dir)=accum_PARM%fablo(dir)
       fabhi_local(dir)=accum_PARM%fabhi(dir)
      enddo

      call checkbound(fablo_local,fabhi_local,DIMS(LS),2,-1,1271)
      call checkbound(fablo_local,fabhi_local,DIMS(matrixfab),0,-1,1271)
      call checkbound(fablo_local,fabhi_local,DIMS(XDISP_fab),2,-1,1271)

      data_out%data_interp=>cell_data_interp

      data_in%level=accum_PARM%level
      data_in%finest_level=accum_PARM%finest_level
      data_in%bfact=accum_PARM%bfact
      data_in%nmat=accum_PARM%nmat

      data_in%dx=>dx_local
      data_in%xlo=>xlo_local
      data_in%fablo=>fablo_local
      data_in%fabhi=>fabhi_local
      data_in%ngrowfab=2

      data_in%state=>XDISP_fab
      data_in%LS=>LS ! not used

      data_in%ncomp=SDIM
      data_in%scomp=1

      nhalf=3

      eps=accum_PARM%dx(1)/10.0d0
      if (eps.gt.zero) then
       ! do nothing
      else
       print *,"eps invalid"
       stop
      endif

      if (accum_PARM%Npart.ge.0) then
       npart_local=accum_PARM%Npart
      else
       print *,"accum_PARM%Npart invalid"
       stop
      endif

      do interior_ID=1,npart_local

       if (accum_PARM%Npart.ge.0) then
        do dir=1,SDIM
         xpart(dir)=accum_PARM%particles(interior_ID)%pos(dir)
         xpartfoot(dir)=accum_PARM%particles(interior_ID)%extra_state(dir)
         xdisp(dir)=xpart(dir)-xpartfoot(dir)
        enddo

        data_in%xtarget=>xpart
        data_in%interp_foot_flag=0
        call interp_from_grid_util(data_in,data_out)

        call containing_cell(accum_PARM%bfact, &
          accum_PARM%dx, &
          accum_PARM%xlo, &
          accum_PARM%fablo, &
          xpart, &
          cell_index)

        interior_ok=1
        do dir=1,SDIM
         if ((cell_index(dir).lt.accum_PARM%tilelo(dir)).or. &
             (cell_index(dir).gt.accum_PARM%tilehi(dir))) then
          interior_ok=0
         endif
        enddo

        if (interior_ok.eq.1) then
         i=cell_index(1)
         j=cell_index(2)
         k=cell_index(SDIM)
         call gridsten_level(xsten,i,j,k,accum_PARM%level,nhalf)
         tmp=0.0d0
         do dir=1,SDIM
          xc(dir)=xsten(0,dir)
          tmp=tmp+(xpart(dir)-xc(dir))**2
         enddo
         tmp=sqrt(tmp)
         w_p=(1.0d0/(eps+tmp))

         if (w_p.gt.zero) then
          matrixfab(D_DECL(i,j,k),1)= &
           matrixfab(D_DECL(i,j,k),1)+w_p
          do dir=1,SDIM
           matrixfab(D_DECL(i,j,k),1+dir)= &
            matrixfab(D_DECL(i,j,k),1+dir)+ &
            w_p*(data_out%data_interp(dir)-xdisp(dir))
          enddo
         else
          print *,"w_p invalid"
          stop
         endif
        else if (interior_ok.eq.0) then
         ! do nothing
        else
         print *,"interior_ok invalid"
         stop
        endif

       else
        print *,"accum_PARM%Npart invalid"
        stop
       endif

      enddo ! do interior_ID=1,accum_PARM%Npart

      return
      end subroutine traverse_particles

       ! called from NavierStokes.cpp:
       !  NavierStokes::accumulate_PC_info(int im_elastic)
       ! 1. isweep==0: gets displacement data from particle data 
       !    and Eulerian data.
       ! 2. isweep==1: calculates the elastic stress tensor from 
       !    u=X(t,x0)-x0
      subroutine fort_assimilate_tensor_from_particles( &
        particles_weight_XD, &
        im_PLS_cpp, & ! 0..nmat-1
        isweep, &
        tid, &  ! thread id
        tilelo,tilehi, &  ! tile box dimensions
        fablo,fabhi, &    ! fortran array box dimensions containing the tile
        bfact, &          ! space order
        level, &          ! 0<=level<=finest_level
        finest_level, &
        xlo,dx, &         ! xlo is lower left hand corner coordinate of fab
        particles_no_nbr, & 
        particles_nbr, & 
        particles_only_nbr, & 
        Np_no_nbr, & ! pass by value
        Np_nbr, & ! pass by value
        Nn, & ! pass by value
        ncomp_tensor, &  ! ncomp_tensor=4 in 2D (11,12,22,33) and 6 in 3D 
        matrix_points, & 
        RHS_points, &    
        ncomp_accumulate, & ! matrix_points+sdim * RHS_points
        nmat, &
        LS, &
        DIMS(LS), &
        TNEWfab, &       ! FAB that holds elastic tensor, Q, when complete
        DIMS(TNEWfab), &
        XDNEWfab, &       
        DIMS(XDNEWfab), &
        XDISP_fab, &      
        DIMS(XDISP_fab), &
        matrixfab, &     ! accumulation FAB
        DIMS(matrixfab)) &
      bind(c,name='fort_assimilate_tensor_from_particles')

      use global_utility_module
      use probcommon_module
      use godunov_module
      implicit none

      INTEGER_T, intent(in) :: im_PLS_cpp
      INTEGER_T, intent(in) :: isweep
      INTEGER_T, intent(in) :: nmat
      REAL_T, intent(in) :: particles_weight_XD
      INTEGER_T, intent(in) :: ncomp_tensor
      INTEGER_T, intent(in) :: matrix_points
      INTEGER_T, intent(in) :: RHS_points
      INTEGER_T, intent(in) :: ncomp_accumulate
      INTEGER_T, value, intent(in) :: Np_no_nbr,Np_nbr,Nn ! pass by value
      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in), target :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in), target :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      REAL_T, intent(in), target :: xlo(SDIM)
      REAL_T, intent(in), target :: dx(SDIM)
      INTEGER_T, intent(in) :: DIMDEC(LS) 
      INTEGER_T, intent(in) :: DIMDEC(matrixfab) 
      INTEGER_T, intent(in) :: DIMDEC(TNEWfab) 
      INTEGER_T, intent(in) :: DIMDEC(XDNEWfab) 
      INTEGER_T, intent(in) :: DIMDEC(XDISP_fab) 
      REAL_T, intent(inout) :: matrixfab( &
        DIMV(matrixfab), &
        ncomp_accumulate)
      REAL_T, intent(in) :: LS( &  
        DIMV(LS), &
        nmat*(1+SDIM))
      REAL_T, intent(inout) :: TNEWfab( &  ! Q assimilated from particles/cells
        DIMV(TNEWfab), &
        ncomp_tensor)
      REAL_T, intent(inout) :: XDNEWfab( &  
        DIMV(XDNEWfab), &
        SDIM)
      REAL_T, intent(in) :: XDISP_fab( &
        DIMV(XDISP_fab), &
        SDIM)
      type(particle_t), intent(in), target :: particles_no_nbr(Np_no_nbr)
      type(particle_t), intent(in) :: particles_nbr(Np_nbr)
      type(particle_t), intent(in), target :: particles_only_nbr(Nn)

      type(accum_parm_type) :: accum_PARM

      INTEGER_T growlo(3)
      INTEGER_T growhi(3)
      INTEGER_T i,j,k
      INTEGER_T dir
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf

      REAL_T A_matrix,B_matrix
      REAL_T lambda
      REAL_T XDISP_local
      INTEGER_T LS_or_VOF_flag
      INTEGER_T im_elastic
      REAL_T local_wt
      INTEGER_T interior_ID
      REAL_T xpart1,xpart2

      nhalf=3

      if (nmat.eq.num_materials) then
       ! do nothing
      else
       print *,"nmat invalid"
       stop
      endif

      if ((im_PLS_cpp.ge.0).and.(im_PLS_cpp.lt.nmat)) then
       ! do nothing
      else
       print *,"im_PLS_cpp invalid"
       stop
      endif

       ! 6 in 3D, 4 in 2D
      if (ncomp_tensor.eq.2*SDIM) then
       ! do nothing
      else
       print *,"ncomp_tensor invalid"
       stop
      endif
      if (matrix_points.eq.1) then
       ! do nothing
      else
       print *,"matrix_points invalid"
       stop
      endif
      if (RHS_points.eq.1) then
       ! do nothing
      else
       print *,"RHS_points invalid"
       stop
      endif
      if (ncomp_accumulate.eq.matrix_points+SDIM*RHS_points) then
       ! do nothing
      else
       print *,"ncomp_accumulate invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(LS),2,-1,1271)
      call checkbound(fablo,fabhi,DIMS(matrixfab),0,-1,1271)
      call checkbound(fablo,fabhi,DIMS(TNEWfab),1,-1,1271)
      call checkbound(fablo,fabhi,DIMS(XDNEWfab),1,-1,1271)
      call checkbound(fablo,fabhi,DIMS(XDISP_fab),2,-1,1271)

      accum_PARM%fablo=>fablo 
      accum_PARM%fabhi=>fabhi
      accum_PARM%tilelo=>tilelo 
      accum_PARM%tilehi=>tilehi
      accum_PARM%bfact=bfact
      accum_PARM%level=level
      accum_PARM%finest_level=finest_level
      accum_PARM%matrix_points=matrix_points
      accum_PARM%RHS_points=RHS_points
      accum_PARM%ncomp_accumulate=ncomp_accumulate
      accum_PARM%dx=>dx
      accum_PARM%xlo=>xlo
      accum_PARM%nmat=nmat

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      if (Np_no_nbr+Nn.eq.Np_nbr) then
       ! do nothing
      else
       print *,"Np_no_nbr+Nn.ne.Np_nbr"
       stop
      endif
       ! sanity check
      do interior_ID=1,Np_no_nbr
       do dir=1,SDIM
        xpart1=particles_nbr(interior_ID)%pos(dir)
        xpart2=particles_no_nbr(interior_ID)%pos(dir)
        if (abs(xpart1-xpart2).le.dx(dir)*VOFTOL) then
         ! do nothing
        else
         print *,"xpart1 or xpart2 invalid"
         stop
        endif
       enddo
      enddo 

      if (isweep.eq.0) then

       accum_PARM%particles=>particles_no_nbr
       accum_PARM%Npart=Np_no_nbr

       call traverse_particles( &
         accum_PARM, &
         matrixfab, &
         DIMS(matrixfab), &
         LS, &
         DIMS(LS), &
         XDISP_fab, &
         DIMS(XDISP_fab), &
         ncomp_accumulate)

       accum_PARM%particles=>particles_only_nbr
       accum_PARM%Npart=Nn

       call traverse_particles( &
         accum_PARM, &
         matrixfab, &
         DIMS(matrixfab), &
         LS, &
         DIMS(LS), &
         XDISP_fab, &
         DIMS(XDISP_fab), &
         ncomp_accumulate)

       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        call gridsten_level(xsten,i,j,k,level,nhalf)
        A_matrix=matrixfab(D_DECL(i,j,k),1) ! sum w(xp)
        do dir=1,SDIM
         B_matrix=matrixfab(D_DECL(i,j,k),1+dir) ! sum w*(X_cell(xp)-X_cell_p)
         XDISP_local=XDISP_fab(D_DECL(i,j,k),dir)
         if (A_matrix.eq.zero) then
          XDNEWFAB(D_DECL(i,j,k),dir)=XDISP_local
         else if (A_matrix.gt.zero) then
           ! lambda=sum (interp(XD)-XD_p)w_p/sum w_p
          lambda=B_matrix/A_matrix
          local_wt=particles_weight_XD
          if ((local_wt.ge.zero).and.(local_wt.le.one)) then
           XDNEWFAB(D_DECL(i,j,k),dir)= &
            XDISP_local-local_wt*lambda
          else
           print *,"local_wt invalid"
           stop
          endif
         else
          print *,"A_matrix invalid"
          stop
         endif
        enddo ! dir=1..SDIM
       enddo
       enddo
       enddo

      else if (isweep.eq.1) then

! grad u=| u_r  u_t/r-v/r  u_z  |
!        | v_r  v_t/r+u/r  v_z  |
!        | w_r  w_t/r      w_z  |
! in RZ:  T33 gets u/r=x_displace/r
! in RTZ: T12=u_t/r - v/r
!         T22=v_t/r + u/r
! later:
! div S = | (r S_11)_r/r + (S_12)_t/r - S_22/r  + (S_13)_z |
!         | (r S_21)_r/r + (S_22)_t/r + S_12/r  + (S_23)_z |
!         | (r S_31)_r/r + (S_32)_t/r +           (S_33)_z |


       LS_or_VOF_flag=0 ! =0 => use LS for upwinding near interfaces
       im_elastic=im_PLS_cpp+1
       call local_tensor_from_xdisplace( &
        LS_or_VOF_flag, &
        im_elastic, &
        tilelo,tilehi, &  ! tile box dimensions
        fablo,fabhi, &    ! fortran array box dimensions containing the tile
        bfact, &          ! space order
        level, &          ! 0<=level<=finest_level
        finest_level, &
        xlo,dx, &         ! xlo is lower left hand corner coordinate of fab
        ncomp_tensor, &  ! ncomp_tensor=4 in 2D (11,12,22,33) and 6 in 3D 
        nmat, &
        LS, &
        DIMS(LS), &
        LS, &  ! VOF placeholder
        DIMS(LS), & ! VOF placeholder
        TNEWfab, &       ! FAB that holds elastic tensor, Q, when complete
        DIMS(TNEWfab), &
        XDISP_fab, &      
        DIMS(XDISP_fab)) 

      else 
       print *,"isweep invalid"
       stop
      endif

      end subroutine fort_assimilate_tensor_from_particles


       ! called from NavierStokes.cpp:
       !  NavierStokes::accumulate_PC_info(int im_elastic)
      subroutine fort_assimilate_tensor_from_xdisplace( &
        im_PLS_cpp, & ! 0..nmat-1
        tid, &  ! thread id
        tilelo,tilehi, &  ! tile box dimensions
        fablo,fabhi, &    ! fortran array box dimensions containing the tile
        bfact, &          ! space order
        level, &          ! 0<=level<=finest_level
        finest_level, &
        xlo,dx, &         ! xlo is lower left hand corner coordinate of fab
        ncomp_tensor, &  ! ncomp_tensor=4 in 2D (11,12,22,33) and 6 in 3D 
        nmat, &
        LS, &
        DIMS(LS), &
        TNEWfab, &       ! FAB that holds elastic tensor, Q, when complete
        DIMS(TNEWfab), &
        XDISP_fab, &      
        DIMS(XDISP_fab)) &
      bind(c,name='fort_assimilate_tensor_from_xdisplace')

      use global_utility_module
      use probcommon_module
      use godunov_module
      implicit none

      INTEGER_T, intent(in) :: im_PLS_cpp
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: ncomp_tensor
      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in), target :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in), target :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      REAL_T, intent(in), target :: xlo(SDIM)
      REAL_T, intent(in), target :: dx(SDIM)
      INTEGER_T, intent(in) :: DIMDEC(LS) 
      INTEGER_T, intent(in) :: DIMDEC(TNEWfab) 
      INTEGER_T, intent(in) :: DIMDEC(XDISP_fab) 
      REAL_T, intent(in) :: LS( &  
        DIMV(LS), &
        nmat*(1+SDIM))
      REAL_T, intent(inout) :: TNEWfab( &  ! Q assimilated from particles/cells
        DIMV(TNEWfab), &
        ncomp_tensor)
      REAL_T, intent(in) :: XDISP_fab( &
        DIMV(XDISP_fab), &
        SDIM)

      INTEGER_T LS_or_VOF_flag
      INTEGER_T im_elastic

      if (nmat.eq.num_materials) then
       ! do nothing
      else
       print *,"nmat invalid"
       stop
      endif

      if ((im_PLS_cpp.ge.0).and.(im_PLS_cpp.lt.nmat)) then
       ! do nothing
      else
       print *,"im_PLS_cpp invalid"
       stop
      endif

       ! 6 in 3D, 4 in 2D
      if (ncomp_tensor.eq.2*SDIM) then
       ! do nothing
      else
       print *,"ncomp_tensor invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(LS),2,-1,1271)
      call checkbound(fablo,fabhi,DIMS(TNEWfab),1,-1,1271)
      call checkbound(fablo,fabhi,DIMS(XDISP_fab),2,-1,1271)

      LS_or_VOF_flag=0  ! use LS for upwinding at interfaces
      im_elastic=im_PLS_cpp+1
       ! e.g. (grad XD + (grad XD)^{T})/2
       ! XD=X(x0,t) - x0  => displacement
       ! XD_{t} + u dot grad XD = u   u=velocity
       ! if XD approximated on an Eulerian mesh, then "reversibility property"
       ! is lost. 
       ! if XD approximated on a Lagrangian mesh, then "synchronous coupling
       ! property" is lost.
       ! Strategy: hybrid Eulerian, Lagrangian method for FSI
       ! prior hybrid approaches: "material point method"
       ! methods by J. Teran (UCLA), others ....
       ! note: Henshaw has already proven that classical asynchronous coupling
       ! results in a method that is unconditionally unstable.
       ! in: GODUNOV_3D.F90
      call local_tensor_from_xdisplace( &
        LS_or_VOF_flag, &
        im_elastic, &
        tilelo,tilehi, &  ! tile box dimensions
        fablo,fabhi, &    ! fortran array box dimensions containing the tile
        bfact, &          ! space order
        level, &          ! 0<=level<=finest_level
        finest_level, &
        xlo,dx, &         ! xlo is lower left hand corner coordinate of fab
        ncomp_tensor, &  ! ncomp_tensor=4 in 2D (11,12,22,33) and 6 in 3D 
        nmat, &
        LS, &
        DIMS(LS), &
        LS, &  ! VOF placeholder
        DIMS(LS), & ! VOF placeholder
        TNEWfab, &       ! FAB that holds elastic tensor, Q, when complete
        DIMS(TNEWfab), &
        XDISP_fab, &      
        DIMS(XDISP_fab)) 

      end subroutine fort_assimilate_tensor_from_xdisplace


      end module FSI_PC_module

