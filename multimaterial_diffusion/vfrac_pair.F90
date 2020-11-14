#undef BL_LANG_CC
#define BL_LANG_FORT

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"


! GeneralClass is in vof_cisl.F90
module MOF_pair_module
  USE GeneralClass
  use probcommon_module
  USE probmain_module 
  use mof_routines_module
  use geometry_intersect_module
  use tsat_module

  implicit none

  abstract interface
    subroutine sub_interface(imat,x,y,dist,probtype_in)
    integer,intent(in)               :: imat
    integer,intent(in)               :: probtype_in
    real(kind=8), intent(in)         :: x,y
    real(kind=8), intent(out)        :: dist
    end subroutine
  end interface


  REAL_T VP_max_LS_error
  INTEGER_T VP_i_max
  INTEGER_T VP_j_max

  INTEGER_T VP_i_debug
  INTEGER_T VP_j_debug

   ! VP_dir_max=VP_side_max=0 if internal flux
  INTEGER_T VP_dir_max
  INTEGER_T VP_side_max
  INTEGER_T VP_im_in_max
  INTEGER_T VP_im_out_max
  INTEGER_T VP_test_number_max

  INTEGER_T VP_i_current
  INTEGER_T VP_j_current
  INTEGER_T VP_dir_current
  INTEGER_T VP_side_current
  INTEGER_T VP_im_in_current
  INTEGER_T VP_im_out_current
  INTEGER_T VP_test_number_current
  REAL_T VP_areaface_max
  REAL_T VP_areaface_current
  REAL_T VP_centroid_face(2)

contains


  !===================================================================
  !   calc all the pairs along the faces of a cell
  !====================================================================
  subroutine vfrac_pair_cell( &
   nhalf, &
   nmat,sdim,dx, &
   ngeom_recon, &
   xsten_cell_sten, &
   mofdata_sten, &
   ext_face_sten, &
   thin_cen_sten, &  ! absolute coordinates
   frac_pair_cell, &
   x_pair_cell)

    implicit none

    integer, intent(in) :: nhalf
    integer, intent(in) :: sdim,nmat,ngeom_recon
    real(kind=8), intent(in) :: dx(sdim)
    REAL*8, intent(in) :: xsten_cell_sten(-1:1,-1:1,-nhalf:nhalf,sdim)
    REAL*8 :: xsten_left(-nhalf:nhalf,sdim)
    REAL*8 :: xsten_right(-nhalf:nhalf,sdim)
    REAL*8, intent(in) :: mofdata_sten(-1:1,-1:1,ngeom_recon*nmat)
    REAL*8 :: mofdata_left(ngeom_recon*nmat)
    REAL*8 :: mofdata_right(ngeom_recon*nmat)
    real(kind=8), intent(in) :: ext_face_sten(-1:1,sdim,nmat+1,sdim,2)
     ! absolute coordinates
    real(kind=8), intent(in)  :: thin_cen_sten(-1:1,sdim,sdim,nmat,sdim,2) 
    real(kind=8), intent(out) :: frac_pair_cell(nmat,nmat,sdim,2)
    real(kind=8), intent(out) :: x_pair_cell(nmat,nmat,sdim,sdim,2)

    real*8         :: vol_total
    integer        :: caller_id
    integer        :: outside_side,outside_side_nbr
    integer        :: left_side,right_side
    integer        :: inside_side_nbr
    integer        :: im
    integer        :: im_outside,im_inside
    integer        :: im_left,im_right
    integer        :: dir,side
    integer        :: dir2
    real(kind=8)   :: frac_outside(nmat),frac_inside(nmat)
    real(kind=8)   :: x_outside(sdim,nmat)
    real(kind=8)   :: x_inside(sdim,nmat)

    real(kind=8)   :: frac_pair(nmat,nmat)
    real(kind=8)   :: x_pair(nmat,nmat,sdim)
    real(kind=8)   :: L_face
    integer        :: tessellate
    integer        :: iii,jjj
    integer        :: ii
    integer        :: bfact,nmax
    REAL(kind=8)   :: xtrilistuncapt1(SDIM+1,SDIM,POLYGON_LIST_MAX)
    REAL(kind=8)   :: xtrilistuncapt2(SDIM+1,SDIM,POLYGON_LIST_MAX)

    bfact = 1
    nmax=POLYGON_LIST_MAX

    if (ngeom_recon.ne.2*sdim+3) then
     print *,"ngeom_recon invalid"
     stop
    endif
    if (nhalf.eq.3) then
     ! do nothing
    else
     print *,"nhalf invalid"
     stop
    endif

    tessellate=0

    do im_outside=1,nmat
    do im_inside=1,nmat
    do dir=1,sdim
    do side=1,2
     frac_pair_cell(im_outside,im_inside,dir,side) = 0.0d0 
     do dir2=1,sdim
      x_pair_cell(im_outside,im_inside,dir2,dir,side) = 0.0d0  
     enddo
    enddo
    enddo
    enddo
    enddo
    do dir = 1, sdim
     iii=0
     jjj=0
     if (dir.eq.1) then
      L_face=dx(2)
      iii=1
     else if (dir.eq.2) then
      L_face=dx(1)
      jjj=1
     else
      print *,"dir invalid"
      stop
     endif

     do side = 1,2

      if (side.eq.1) then
       outside_side=-1
       outside_side_nbr=2
       left_side=-1
       right_side=0
      else if (side.eq.2) then
       outside_side=1
       outside_side_nbr=1
       left_side=0
       right_side=1
      else
       print *,"side invalid"
       stop
      endif
      inside_side_nbr=3-outside_side_nbr

      do im=1,nmat*ngeom_recon
       mofdata_left(im)=mofdata_sten(iii*left_side,jjj*left_side,im)
       mofdata_right(im)=mofdata_sten(iii*right_side,jjj*right_side,im)
      enddo
      do ii=-nhalf,nhalf
      do dir2=1,sdim
       xsten_left(ii,dir2)= &
               xsten_cell_sten(iii*left_side,jjj*left_side,ii,dir2)
       xsten_right(ii,dir2)= &
               xsten_cell_sten(iii*right_side,jjj*right_side,ii,dir2)
      enddo
      enddo
      do im = 1, nmat
       frac_outside(im) = ext_face_sten(outside_side,dir,im,dir, &
               outside_side_nbr)
       frac_inside(im) = ext_face_sten(0,dir,im,dir, &
               inside_side_nbr)  ! inside
      enddo

      do dir2=1,sdim
      do im = 1, nmat
       x_outside(dir2,im)=thin_cen_sten(outside_side,dir,dir2,im,dir, &
               outside_side_nbr)
       x_inside(dir2,im)=thin_cen_sten(0,dir,dir2,im,dir,inside_side_nbr)
      enddo
      enddo
      if ((VP_i_debug.eq.VP_i_current).and. &
          (VP_j_debug.eq.VP_j_current)) then
       print *,"calling multi_get_area_pairs  dir,side= ",dir,side
      endif

       ! x_pair in absolute coordinate system.
      caller_id=12
      call multi_get_area_pairs( &
        bfact,dx, &
        xsten_right, &  ! aka xsten0_plus
        xsten_left, &   ! aka xsten0_minus
        nhalf, &
        mofdata_right, & ! aka mofdata_plus
        mofdata_left, &  ! aka mofdata_minus
        nmat, &
        dir, &  ! dir=1..sdim
        frac_pair, & ! left,right
        x_pair, & ! left,right
        sdim, &
        xtrilistuncapt1, &
        nmax, &
        xtrilistuncapt2, &
        nmax, &
        nmax, &
        caller_id)

      vol_total=zero
      do im_left=1,nmat
      do im_right=1,nmat
       vol_total=vol_total+frac_pair(im_left,im_right)
      enddo
      enddo
      if (vol_total.gt.zero) then

       do im_left=1,nmat
       do im_right=1,nmat
        frac_pair(im_left,im_right)= &
                frac_pair(im_left,im_right)/vol_total
       enddo
       enddo

      else
       print *,"vol_total invalid"
       stop
      endif


      do im_outside=1,nmat
      do im_inside=1,nmat
       if (side.eq.1) then
        im_left=im_outside
        im_right=im_inside
       else if (side.eq.2) then
        im_left=im_inside
        im_right=im_outside
       else
        print *,"side invalid"
        stop
       endif

       frac_pair_cell(im_outside,im_inside,dir,side) = &
               frac_pair(im_left,im_right)
       do dir2=1,sdim
        x_pair_cell(im_outside,im_inside,dir2,dir,side) = &
               x_pair(im_left,im_right,dir2)
       enddo
      enddo
      enddo

      if ((VP_i_debug.eq.VP_i_current).and. &
          (VP_j_debug.eq.VP_j_current)) then
       print *,"i,j,dir,side (vfrac_pair) ",VP_i_current,VP_j_current,dir,side
       do im_outside=1,nmat
       do im_inside=1,nmat
        if (frac_pair_cell(im_outside,im_inside,dir,side).gt.0.0d0) then
         print *,"im_outside,im_inside,dir,side,frac_pair_cell ", &
                im_outside,im_inside,dir,side, &
                frac_pair_cell(im_outside,im_inside,dir,side)
         do dir2=1,sdim
          print *,"im_outside,im_inside,dir2,dir,side,x_pair_cell ", &
              im_outside,im_inside,dir2,dir,side, &
              x_pair_cell(im_outside,im_inside,dir2,dir,side)
         enddo
        endif
       enddo 
       enddo 

       do im = 1, nmat
        if (frac_outside(im).gt.0.0d0) then
         print *,"im,frac_outside ",im,frac_outside(im)
         do dir2=1,sdim
          print *,"dir2,im,x_outside ",dir2,im,x_outside(dir2,im)
         enddo
        endif
        if (frac_inside(im).gt.0.0d0) then
         print *,"im,frac_inside ",im,frac_inside(im)
         do dir2=1,sdim
          print *,"dir2,im,x_inside ",dir2,im,x_inside(dir2,im)
         enddo
        endif
       enddo

      endif

     enddo ! side
    enddo ! dir

  end subroutine vfrac_pair_cell
  !/////////////////////////////////////////////////////////////////////////
  !/////////////////////////////////////////////////////////////////////////

  ! =========================================================================
  ! =========================================================================
  !   From legacy code, perturb interface and cell faces
  ! 
  ! =========================================================================
  !=========================================================================

  !------------------perturb faces of the cell-------------------
  subroutine ptb_ext(ngeom_recon, &
   nmat,sdim,dx,mofdata, &
   xsten, &
   ext_facefrac_cell, &
   multi_cen_cell)
    implicit none

    integer,         intent(in) :: ngeom_recon,nmat, sdim
    real(kind=8)                :: h       !  >__<
    real(kind=8),    intent(in) :: mofdata(ngeom_recon*nmat)
    real(kind=8),    intent(in) :: dx(sdim)
    real(kind=8)                :: mofdataproject(nmat*ngeom_recon)

    real(kind=8)                :: vcenter(nmat)
    INTEGER                     :: sorted_list(nmat)
    real(kind=8)                :: xsten(-3:3,SDIM)
    real(kind=8)                :: xsten_thin(-1:1,SDIM)
    REAL(kind=8)                :: dxthin
    integer                     :: vofcomp
    REAL(kind=8)                :: ext_facefrac_cell(nmat+1,sdim,2)

    integer                     :: im
    integer                     :: dir_side
    integer                     :: dir_local
    integer                     :: side,ivert,isten
    integer                     :: im_crit,im_exclude
    integer                     :: shapeflag
    integer                     :: tessellate
    integer                     :: nhalf, nhalf_thin
    integer                     :: nmax
    REAL(kind=8)                :: multi_volume(nmat)
    REAL(KIND=8)                :: multi_cen(SDIM,nmat)
    REAL(KIND=8)                :: multi_area(nmat)
    REAL(KIND=8)                :: total_vol
    REAL(kind=8)                :: xtrilistuncapt(SDIM+1,SDIM,POLYGON_LIST_MAX)
    REAL(kind=8)                :: dummy_tri(SDIM+1,SDIM)  
    real(kind=8)                :: multi_cen_cell(sdim,nmat,sdim,2) 

    integer                     :: bfact

    if (ngeom_recon.eq.2*sdim+3) then
     ! do nothing
    else
     print *,"ngeom_recon invalid"
     stop
    endif

    h = dx(1)

    bfact = 1
    nhalf=3 
    nhalf_thin=1 
    nmax=POLYGON_LIST_MAX

    do dir_local=1,sdim
    do im=1,nmat
    do dir_side=1,sdim
    do side=1,2
     multi_cen_cell(dir_local,im,dir_side,side) = 0.0d0
    enddo
    enddo
    enddo
    enddo
    do im=1,nmat+1
    do dir_side=1,sdim
    do side=1,2
     ext_facefrac_cell(im,dir_side,side) = 0.0d0
    enddo
    enddo
    enddo

    do ivert=1,SDIM+1
       do dir_local=1,SDIM
          dummy_tri(ivert,dir_local)=zero
       enddo
    enddo

    ! vcenter = volume fraction 
    do im=1,nmat
       vofcomp=(im-1)*ngeom_recon+1
       vcenter(im)=mofdata(vofcomp)    ! extract volume fraction
    enddo ! im
    ! sort the volume fractions.  
    im_exclude = 0
    call sort_vof(vcenter,im_exclude,sorted_list,nmat)
    im_crit=sorted_list(1)  ! vcenter(im_crit) is largest value.
    if ((im_crit.ge.1).and.(im_crit.le.nmat)) then
     ! do nothing
    else
     print *,"im_crit invalid"
     stop
    endif

    ! ext_facefrac_cell: nmat+1,sdim,2 components
    ! last 2*sdim components are areas of the cell face (dir,side)
    ! e.g. in 2d uniform grid: ext_facefrac_cell=h
    ! first nmat*2*sdim components are areafractions. 0<=A<=1
    ! *********cal  tot  area  of faces
    do dir_side=1,SDIM
     do side=1,2
      ext_facefrac_cell(nmat+1,dir_side,side) = h      ! 2d
     enddo
    enddo ! dir_side,side
    ! *************************

    if (vcenter(im_crit) .ge. one-FACETOL_DVOL) then
     do dir_side=1,SDIM
      dxthin=FACETOL_DVOL*(xsten(1,dir_side)-xsten(-1,dir_side))
      if (dxthin.gt.0.0d0) then
       do side=1,2
        ext_facefrac_cell(im_crit,dir_side,side)=one
        if(side .eq. 1)then
         do dir_local = 1,sdim
          if(dir_local .eq. dir_side)then
           multi_cen_cell(dir_local,im_crit,dir_side,side)=xsten(-1,dir_side) 
          else
           multi_cen_cell(dir_local,im_crit,dir_side,side)=xsten(0,dir_local)
          endif
         enddo
        elseif(side .eq. 2)then
         do dir_local = 1,sdim
          if(dir_local .eq. dir_side)then
           multi_cen_cell(dir_local,im_crit,dir_side,side)=xsten(1,dir_side) 
          else
           multi_cen_cell(dir_local,im_crit,dir_side,side)=xsten(0,dir_local)
          endif
         enddo
        else
         print *,"side invalid"
         stop
        endif
       enddo ! side=1..2
      else
       print *,"dxthin invalid"
       stop
      endif
     enddo  ! dir_side
       !----------------------------------------------
       !   multi_cen_cell(sdim,im_crit,sdim,2) 
       !---------------------------------------------
    else if (vcenter(im_crit).gt.zero) then
     shapeflag=0
     tessellate=1
     do dir_side=1,SDIM
     do side=1,2
      do isten=-1,1
       do dir_local=1,SDIM
        xsten_thin(isten,dir_local)=xsten(isten,dir_local)
       enddo
      enddo
      dxthin=FACETOL_DVOL*(xsten(1,dir_side)-xsten(-1,dir_side))
      if (dxthin.gt.0.0d0) then
       if (side.eq.1) then
        xsten_thin(1,dir_side)=xsten(-1,dir_side)+dxthin
       else if (side.eq.2) then
        xsten_thin(-1,dir_side)=xsten(1,dir_side)-dxthin
       else
        print *,"side invalid"
        stop
       endif
       xsten_thin(0,dir_side)= &
          half*(xsten_thin(-1,dir_side)+xsten_thin(1,dir_side))
       ! find volumes and areas (not scaled) of materials in
       ! xsten_thin box: xsten_thin(0,dir_side) = center of thin box
       ! xsten_thin(1,dir_side) right side in dir_side direction
       ! xsten_thin(-1,dir_side) left side 
       call project_slopes_to_face( &
        bfact,dx,xsten,nhalf, &
        mofdata,mofdataproject, &
        nmat,SDIM,dir_side,side)

       call multi_get_volume_grid( &
        tessellate, &
        bfact,dx,xsten,nhalf, &
        mofdataproject, &  ! was mofdata
        xsten_thin,nhalf_thin, &
        dummy_tri, &
        multi_volume, &
        multi_cen, &  ! absolute coordinates
        multi_area, &
        xtrilistuncapt, &
        nmax,nmax,nmat,SDIM, &
        shapeflag,3)

       do im=1,nmat
        do dir_local=1,SDIM
         if (dir_local.eq.dir_side) then
          if (side.eq.1) then
           multi_cen_cell(dir_local,im,dir_side,side)=xsten(-1,dir_side)
          else if (side.eq.2) then
           multi_cen_cell(dir_local,im,dir_side,side)=xsten(1,dir_side)
          else
           print *,"side invalid"
           stop
          endif
         else 
          multi_cen_cell(dir_local,im,dir_side,side)=multi_cen(dir_local,im)
         endif
        enddo  ! dir_local=1..sdim
       enddo ! im=1..nmat

       total_vol=zero
       do im=1,nmat
        total_vol=total_vol+multi_volume(im)
       enddo
       if (total_vol.gt.zero) then
        do im=1,nmat
         ext_facefrac_cell(im,dir_side,side)=multi_volume(im)/total_vol
        enddo
       else
        print *,"total_vol invalid"
        stop
       endif

       if ((VP_i_debug.eq.VP_i_current).and. &
           (VP_j_debug.eq.VP_j_current)) then
        print *,"i,j,dir_side,side (ptb_ext) ", &
               VP_i_current,VP_j_current,dir_side,side
        do im=1,nmat
         if (ext_facefrac_cell(im,dir_side,side).gt.0.0d0) then
          print *,"im,dir_side,side,ext_facefrac_cell ",im,dir_side,side, &
           ext_facefrac_cell(im,dir_side,side)
          do dir_local=1,sdim
           print *,"dir_local,im,dir_side,side,multi_cen_cell ", &
            dir_local,im,dir_side,side, &
            multi_cen_cell(dir_local,im,dir_side,side)
          enddo
         endif
        enddo ! im=1..nmat
       endif

      else
       print *,"dxthin invalid"
       stop
      endif

     enddo ! side
     enddo ! dir_side

    else
     print *,"vcenter(im_crit) out of range"
     stop
    endif

    return

  end subroutine ptb_ext
  !------------------------------------------------------------------
  subroutine sort_vof(vfrac_data,im_exclude,sorted_list,nmat)
    implicit none

    integer     ,intent(in)   :: nmat,im_exclude
    real(kind=8),intent(in)   :: vfrac_data(nmat)
    integer                   :: sorted_list(nmat)

    integer                   :: im,changed,nsweeps,swap,do_swap

    do im=1,nmat
       sorted_list(im)=im
    enddo

    changed=1
    nsweeps=0
    do while ((changed.eq.1).and.(nsweeps.lt.nmat-1))
       changed=0
       do im=1,nmat-nsweeps-1
          do_swap=0
          if (sorted_list(im).eq.im_exclude) then
             do_swap=1
          else if (sorted_list(im+1).eq.im_exclude) then
             do_swap=0
          else if (vfrac_data(sorted_list(im)).lt. &
               vfrac_data(sorted_list(im+1))) then

             do_swap=1
          endif
          if (do_swap.eq.1) then
             swap=sorted_list(im)
             sorted_list(im)=sorted_list(im+1)
             sorted_list(im+1)=swap
             changed=1
          endif
       enddo
       nsweeps=nsweeps+1
    enddo

    return

  end subroutine sort_vof


  subroutine get_kappa_simple( &
   sdim, &
   alpha, &
   nmat, &
   im1,im2, &
   kappa, &
   nf, & ! points out of cell
   xinside,xoutside, &
   xinside_face,xoutside_face)
  IMPLICIT NONE

  integer, intent(in)      :: sdim
  integer, intent(in)      :: nmat
  integer, intent(in)      :: im1,im2
  real(kind=8),intent(in)  :: alpha(nmat)
  real(kind=8),intent(out) :: kappa
  real(kind=8),intent(in)  :: nf(sdim)
  real(kind=8),intent(in)  :: xinside(sdim)
  real(kind=8),intent(in)  :: xoutside(sdim)
  real(kind=8),intent(in)  :: xinside_face(sdim)
  real(kind=8),intent(in)  :: xoutside_face(sdim)

  real(kind=8)      :: n1(sdim)
  real(kind=8)      :: n2(sdim)
  real(kind=8)      :: mag_nf,mag_inside,mag_outside
  real(kind=8)      :: dotprod1,dotprod2,k1,k2
  integer           :: dir

  if ((im1.ge.1).and.(im1.le.nmat).and. &
      (im2.ge.1).and.(im2.le.nmat)) then
   mag_nf=zero
   mag_inside=zero
   mag_outside=zero
   do dir=1,sdim
    mag_nf=mag_nf+nf(dir)**2
    n1(dir)=xinside(dir)-xinside_face(dir)
    mag_inside=mag_inside+n1(dir)**2
    n2(dir)=xoutside(dir)-xoutside_face(dir)
    mag_outside=mag_outside+n2(dir)**2
   enddo
   mag_nf=sqrt(mag_nf)
   mag_inside=sqrt(mag_inside)
   mag_outside=sqrt(mag_outside)
   if ((abs(1.0d0-mag_nf).le.VOFTOL).and. &
       (mag_inside.gt.0.0d0).and. &
       (mag_outside.gt.0.0d0)) then
    ! grad T1 approx n1(T1-T1I)/L1  
    ! grad T2 approx n2(T2-T2I)/L2
    ! T1I=T2I
    ! k1 grad T1 dot nf = k2 grad T2 dot nf  
    ! k1 <- k1 n1 dot nf/L1
    ! k2 <- k2 n2 dot nf/L2
    ! k1(T1-T1I)=-k2(T2-T2I)
    ! T1I=T2I
    ! T1I(k2+k1)=k1 T1 + k2 T2
    ! T1I=(k1 T1 + k2 T2)/(k1+k2)
    ! k1(T1-T1I)=k1 T1 - k1 (k1 T1 + k2 T2)/(k1+k2) =
    ! (k2 k1 T1 - k1 k2 T2)/(k1+k2)=
    ! (T1-T2) k1 k2/(k1+k2)
    ! k1 k2/(k1+k2)=1/(1/k1 + 1/k2)=(1/2) * harmonic average(k1,k2)
    ! e.g. k1=k2=1/(h/2)=2/h then 
    ! k1 k2/(k1+k2)= (4/h^2)/(4/h)=1/h
    dotprod1=0.0d0
    dotprod2=0.0d0
    do dir=1,sdim
     n1(dir)=n1(dir)/mag_inside
     n2(dir)=n2(dir)/mag_outside
     dotprod1=dotprod1+n1(dir)*nf(dir)/mag_nf
     dotprod2=dotprod2+n2(dir)*nf(dir)/mag_nf
    enddo
    k1=alpha(im1)*abs(dotprod1)/mag_inside
    k2=alpha(im2)*abs(dotprod2)/mag_outside
    if ((k1.ge.0.0d0).and. &
        (k2.ge.0.0d0).and. &
        (k1+k2.gt.0.0d0)) then
     kappa=k1*k2/(k1+k2)
    else if ((k1.eq.0.0d0).and.(k2.eq.0.0d0)) then
     kappa=0.0d0
    else
     print *,"k1 or k2 invalid"
     stop
    endif
   else
    print *,"mag_nf, mag_inside, or mag_outside invalid"
    stop
   endif
  else
   print *,"im1 or im2 invalid"
   stop
  endif

  return
  end subroutine get_kappa_simple

  subroutine get_kappa_dirichlet( &
   sdim, &
   alpha, &
   nmat, &
   im1, &
   kappa, &
   nf, & ! points out of cell
   xinside, &
   xinside_face)
  IMPLICIT NONE

  integer, intent(in)      :: sdim
  integer, intent(in)      :: nmat
  integer, intent(in)      :: im1
  real(kind=8),intent(in)  :: alpha(nmat)
  real(kind=8),intent(out) :: kappa
  real(kind=8),intent(in)  :: nf(sdim)
  real(kind=8),intent(in)  :: xinside(sdim)
  real(kind=8),intent(in)  :: xinside_face(sdim)

  real(kind=8)      :: n1(sdim)
  real(kind=8)      :: mag_nf,mag_inside
  real(kind=8)      :: dotprod1,k1
  integer           :: dir

  if ((im1.ge.1).and.(im1.le.nmat)) then
   mag_nf=zero
   mag_inside=zero
   do dir=1,sdim
    mag_nf=mag_nf+nf(dir)**2
    n1(dir)=xinside(dir)-xinside_face(dir)
    mag_inside=mag_inside+n1(dir)**2
   enddo
   mag_nf=sqrt(mag_nf)
   mag_inside=sqrt(mag_inside)
   if ((abs(1.0d0-mag_nf).le.VOFTOL).and. &
       (mag_inside.gt.0.0d0)) then
    ! k1(T1-T1I)=-k2(T2-T2I)
    ! T1I=T2I
    ! T1I(k2+k1)=k1 T1 + k2 T2
    ! T1I=(k1 T1 + k2 T2)/(k1+k2)
    ! k1(T1-T1I)=k1 T1 - k1 (k1 T1 + k2 T2)/(k1+k2) =
    ! (k2 k1 T1 - k1 k2 T2)/(k1+k2)=
    ! (T1-T2) k1 k2/(k1+k2)
    ! e.g. k1=k2=1/(h/2)=2/h then 
    ! k1 k2/(k1+k2)= (4/h^2)/(4/h)=1/h
    dotprod1=0.0d0
    do dir=1,sdim
     n1(dir)=n1(dir)/mag_inside
     dotprod1=dotprod1+n1(dir)*nf(dir)/mag_nf
    enddo
    k1=alpha(im1)*abs(dotprod1)/mag_inside
    if (k1.gt.0.0d0) then
     kappa=k1
    else if (k1.eq.0.0d0) then
     kappa=0.0d0
    else
     print *,"k1 invalid"
     stop
    endif
   else
    print *,"mag_nf or mag_inside invalid"
    stop
   endif
  else
   print *,"im1 invalid"
   stop
  endif

  return
  end subroutine get_kappa_dirichlet



  subroutine get_kappa_driver_complicated(alpha,theta1,theta2, &
                 im1,im2,kappa,nmat)
  IMPLICIT NONE

  integer im1,im2,nmat
  real(kind=8),intent(in)     :: alpha(nmat)
  real(kind=8),intent(out)    :: kappa
  real(kind=8),intent(in)     :: theta1,theta2
  real(kind=8)                :: t1,t2,theta,kappa_eps

  kappa_eps=1.0d-8

  t1=theta1
  t2=theta2 
  if ((t1.lt.0.0).or.(t2.lt.0.0)) then
   print *,"t1 or t2 invalid"
   stop
  endif
  theta=t1+t2
  if (theta.eq.0.0) then
   kappa=0.0
  else if (theta.gt.0.0) then
   t1=t1/theta 
   t2=t2/theta 
   if ((t1.le.kappa_eps).or.(t2.le.kappa_eps)) then
    kappa=0.0
   else if (alpha(im1).eq.0.0) then
    kappa=0.0
   else if (alpha(im2).eq.0.0) then
    kappa=0.0
   else if ((alpha(im1).gt.0.0).and.(alpha(im2).gt.0.0)) then
    kappa=(t1+t2)*alpha(im1)*alpha(im2)/ &
     (t1*alpha(im2)+t2*alpha(im1))
   else
    print *,"alpha invalid"
    stop
   endif
  else
   print *,"theta invalid"
   stop
  endif

  return
  end subroutine get_kappa_driver_complicated

    ! iten<0 if im1<>im2 and not a dirichlet boundary
    ! iten=0 if im1==im2
    ! iten>0 if im1<>im2 and a dirichlet boundary
  subroutine init_multimaterial_flag(im1,im2,nmat,internal_multimaterial, &
    im_source,im_dest)
  use global_utility_module, only : get_iten
  IMPLICIT NONE
  integer,intent(in)          :: im1,im2,nmat
  integer                     :: ireverse
  integer                     :: iten
  integer,intent(out)         :: internal_multimaterial
  integer,intent(out)         :: im_source,im_dest
  
  internal_multimaterial=0
  im_source=0
  im_dest=0

  ireverse_tsat=0
  iten_tsat=0

  if (im1.eq.im2) then
    ! do nothing
  else if (im1.ne.im2) then
    call get_iten(im1,im2,iten,nmat)

    iten_tsat=iten

    internal_multimaterial=-iten
    im_source=im1
    im_dest=im2

    if ((iten.ge.1).and.(iten.le.global_nten)) then
     do ireverse=0,1
      if (latent_heat(iten+ireverse*global_nten).eq.zero) then
       ! do nothing
      else if (latent_heat(iten+ireverse*global_nten).ne.zero) then
       ireverse_tsat=ireverse
       internal_multimaterial=iten+ireverse*global_nten
       if (ireverse.eq.0) then
        if (im1.lt.im2) then
         im_source=im1
         im_dest=im2
        else if (im1.gt.im2) then
         im_source=im2
         im_dest=im1
        else
         print *,"im1 or im2 bust"
         stop
        endif
       else if (ireverse.eq.1) then
        if (im1.lt.im2) then
         im_source=im2
         im_dest=im1
        else if (im1.gt.im2) then
         im_source=im1
         im_dest=im2
        else
         print *,"im1 or im2 bust"
         stop
        endif
       else
        print *,"ireverse invalid"
        stop
       endif
      else
       print *,"latent_heat NaN"
       stop
      endif
     enddo
    else
     print *,"iten invalid, init_multimaterial_flag"
     print *,"global_nten ",global_nten
     print *,"iten=",iten

     stop
    endif

    if (1.eq.0) then
     if ((probtype.eq.403).and. &
         (im_source.eq.1).and. &
         (im_dest.eq.2).and. &
         (internal_multimaterial.eq.1)) then
      ! do nothing
     else
      print *,"dendrite parameters invalid (1)"
      stop
     endif
    endif
          
  else
    print *,"im1 bust"
    stop
  endif

  end subroutine init_multimaterial_flag

   ! im1 is the target material
  subroutine get_kappa_complicated(alpha,xsten,dir,side, &
   cen1,cen2,im1,im2,kappa,nmat,sdim)
  IMPLICIT NONE

  integer im1,im2,nmat,dir,side,sdim
  real(kind=8),intent(in)     :: alpha(nmat)
  real(kind=8),intent(in)     :: xsten(-3:3,sdim)
  real(kind=8),intent(in)     :: cen1(sdim)
  real(kind=8),intent(in)     :: cen2(sdim)
  real(kind=8),intent(out)    :: kappa
  real(kind=8)                :: theta1,theta2
  integer                     :: internal_multimaterial
  integer                     :: im_source,im_dest

  if (nmat.eq.num_materials) then

   call init_multimaterial_flag(im1,im2,nmat,internal_multimaterial, &
    im_source,im_dest)

   if (internal_multimaterial.le.0) then 
    if (side.eq.1) then
     theta1=abs(xsten(-1,dir)-cen1(dir)) 
     theta2=abs(xsten(-1,dir)-cen2(dir)) 
    else if (side.eq.2) then
     theta1=abs(xsten(1,dir)-cen1(dir)) 
     theta2=abs(xsten(1,dir)-cen2(dir)) 
    else
     print *,"side invalid"
     stop
    endif

    call get_kappa_driver_complicated(alpha,theta1,theta2, &
            im1,im2,kappa,nmat)
   else if ((internal_multimaterial.ge.1).and. &
            (internal_multimaterial.le.2*global_nten)) then
    kappa=alpha(im1)
   else
    print *,"internal_multimaterial invalid"
    stop
   endif
  else
   print *,"nmat invalid"
   stop
  endif

  return
  end subroutine get_kappa_complicated

  subroutine get_kappa_int_complicated(alpha,dist,im1,im2,kappa,nmat)
  IMPLICIT NONE

  integer im1,im2,nmat
  real(kind=8),intent(in)     :: alpha(nmat)
  real(kind=8),intent(out)    :: kappa
  real(kind=8)                :: dist(nmat,nmat)           
  real(kind=8)                :: theta1,theta2
 
  theta1=dist(im1,im2)
  theta2=dist(im2,im1)
  call get_kappa_driver_complicated(alpha,theta1,theta2, &
          im1,im2,kappa,nmat)

  return
  end subroutine get_kappa_int_complicated

  !---------------------------------------------------------------
  ! dist(im_point1,im_point2)=distance from im_point1 centroid to im_line
  !   reconstructed interface.
  ! dist(im_point2,im_point1)=distance from im_point2 centroid to im_line
  !   reconstructed interface.
  subroutine dist_point_line(sdim,nmat,dx,mofdata, &
   im_line,im_point1,im_point2,dist, &
   xcentroid_cell, &
   xcenter_cell, &
   xclosest)  ! absolute coord
    implicit none

    integer,intent(in)          :: sdim,nmat
    integer,intent(in)          :: im_line
    integer,intent(in)          :: im_point1
    integer,intent(in)          :: im_point2
    real(kind=8),intent(in)     :: mofdata((sdim*2+3)*nmat)
    real(kind=8),intent(in)     :: dx(sdim)
    real(kind=8),intent(in)     :: xcentroid_cell(sdim)
    real(kind=8),intent(in)     :: xcenter_cell(sdim)
    real(kind=8)                :: xdiff(sdim)
    real(kind=8)                :: dist(nmat,nmat)           
    real(kind=8)                :: xclosest(nmat,nmat,sdim)           

    integer                     :: ngeom_recon
    integer                     :: dir,vofcomp
    real(kind=8)                :: m(sdim)
    real(kind=8)                :: cent(sdim)
    real(kind=8)                :: intercept
    real(kind=8)                :: modm
    real(kind=8)                :: LS

    ngeom_recon = sdim*2+3

    if ((im_line.ge.1).and.(im_line.le.nmat).and. &
        (im_point1.ge.1).and.(im_point1.le.nmat).and. &
        (im_point2.ge.1).and.(im_point2.le.nmat)) then

     ! vfrac,cen,order,slope,intercept
     vofcomp=ngeom_recon*(im_line-1)+1
     do dir = 1,sdim
      m(dir) = mofdata(vofcomp+sdim+1+dir) ! slope
     enddo
     intercept = mofdata(vofcomp+2*sdim+2)

     modm = 0.0d0
     do dir = 1,sdim
       modm = modm + m(dir)*m(dir)
     enddo
     modm = sqrt(modm)

     if (modm.gt.0.0d0) then

      vofcomp=ngeom_recon*(im_point1-1)+1
      cent = 0.0d0
      do dir = 1,sdim   
       cent(dir) =  mofdata(vofcomp + dir)
      enddo

      ! centroid coordinate is relative to the centroid of the cell.
      ! x0 is the center of the cell.
      ! LS=n dot (x-x0) + intercept
      ! distance=|n dot (xcentroid_absolute-x0)+intercept|
      if (dist(im_point1,im_point2).eq.zero) then

       do dir=1,sdim
        xdiff(dir)=cent(dir)+xcentroid_cell(dir)-xcenter_cell(dir)
       enddo
       LS=(dot_product(m,xdiff)+intercept)/modm
       dist(im_point1,im_point2) = abs(LS)
       do dir=1,sdim
        xclosest(im_point1,im_point2,dir) =  &
              cent(dir)+xcentroid_cell(dir)-LS*m(dir)/modm
       enddo

       if (dist(im_point1,im_point2).gt.zero) then

        vofcomp=ngeom_recon*(im_point2-1)+1
        cent = 0.0d0
        do dir = 1,sdim   
         cent(dir) =  mofdata(vofcomp + dir)
        enddo

        if (dist(im_point2,im_point1).eq.zero) then

         do dir=1,sdim
          xdiff(dir)=cent(dir)+xcentroid_cell(dir)-xcenter_cell(dir)
         enddo
         LS=(dot_product(m,xdiff)+intercept)/modm
         dist(im_point2,im_point1) = abs(LS)
         do dir=1,sdim
          xclosest(im_point2,im_point1,dir) = &
               cent(dir)+xcentroid_cell(dir)-LS*m(dir)/modm
         enddo

         if (dist(im_point2,im_point1).gt.zero) then

          ! do nothing

         else
          print *,"distance should be positive"
          stop
         endif

        else
         print *,"initial distance should be zero"
         stop
        endif
       else
        print *,"distance should be positive"
        stop
       endif
      else
       print *,"initial distance should be zero"
       stop
      endif

     else
      print *,"modm invalid"
      stop
     endif

    else
     print *,"im_line, im_point1, or im_point2 invalid"
     stop
    endif


  end subroutine dist_point_line
  !------------------------------------------------------------
  !   perturb the interface inside the cell
  !---------------------------------------------------------
  subroutine ptb_int(ngeom_recon, &
    nmat,sdim,mofdata,dx,xsten, &
    int_facefrac_cell, &
    int_centroid_cell, &
    int_face_normal_cell, &
    dist_to_int_cell, &
    xclosest)  ! absolute coordinates
    implicit none

    integer,         intent(in) :: ngeom_recon,nmat, sdim
    real(kind=8),    intent(in) :: mofdata(ngeom_recon*nmat)
    real(kind=8),    intent(in) :: dx(sdim)

    real(kind=8)                :: vcenter(nmat)
    INTEGER                     :: sorted_list(nmat)
    real(kind=8), intent(in)    :: xsten(-3:3,SDIM)
    real(kind=8)                :: xcentroid_cell(SDIM)
    real(kind=8)                :: xcenter_cell(SDIM)
    integer                     :: vofcomp
    REAL(kind=8)                :: int_facefrac_cell(nmat+1,nmat)
    REAL(kind=8)                :: int_centroid_cell(nmat,nmat,sdim)
    REAL(kind=8)                :: int_face_normal_cell(nmat,nmat,sdim)

    integer                     :: im,dir,im_opp
    integer                     :: im1,im2
    integer                     :: im_crit,im_exclude
    integer                     :: ivert
    integer                     :: shapeflag
    integer                     :: tessellate
    integer                     :: nhalf
    integer                     :: nmax
    REAL(kind=8)                :: multi_volume(nmat)
    real(kind=8)                :: multi_volume_offset(nmat)
    REAL(KIND=8)                :: multi_cen(SDIM,nmat)
    REAL(KIND=8)                :: multi_cen_offset(SDIM,nmat)
    REAL(KIND=8)                :: multi_area(nmat)
    REAL(kind=8)                :: xtrilistuncapt(SDIM+1,SDIM,POLYGON_LIST_MAX)
    REAL(kind=8)                :: dummy_tri(SDIM+1,SDIM) 


    real(kind=8)                :: mofdatavalid(ngeom_recon*nmat)
    REAL(kind=8)                :: multi_area_offset(nmat)
    REAL(kind=8)                :: total_face,uncaptured_volume
    INTEGER                     :: irank,testflag
    REAL(kind=8)                :: intercept
    REAL(kind=8)                :: local_facefrac(nmat)
    REAL(kind=8)                :: dist_to_int_cell(nmat,nmat)
    REAL(kind=8)                :: xclosest(nmat,nmat,SDIM) ! absolute coord
    integer                     :: bfact
    real(kind=8)                :: delta_volume
    real(kind=8)                :: LS

    integer                     :: is_processed(nmat)
    integer :: nhalf_box

    nhalf_box=1

    tessellate=1

    if (ngeom_recon.ne.2*sdim+3) then
     print *,"ngeom_recon invalid"
     stop
    endif

    bfact = 1
    nhalf=3
    nmax=POLYGON_LIST_MAX
    !print *,nmax
    !nmax = 400

    if (levelrz.eq.0) then
     do dir=1,SDIM
      xcentroid_cell(dir)=xsten(0,dir)
      xcenter_cell(dir)=xsten(0,dir)
     enddo
    else
     print *,"levelrz==0 only supported here"
     stop
    endif

    do im1=1,nmat
    do im2=1,nmat
     dist_to_int_cell(im1,im2) = 0.0d0
     do dir=1,SDIM
      xclosest(im1,im2,dir)=zero
     enddo
    enddo
    enddo

    do ivert=1,SDIM+1
     do dir=1,SDIM
      dummy_tri(ivert,dir)=zero
     enddo
    enddo
    ! xsten(0,dir)  dir=1,2  is center of cell
    ! xsten(1,dir)  is right coordinate in dir direction
    ! xsten(-1,dir) is left coordinate in dir direction.
    ! vcenter = volume fraction 
    do im=1,nmat
     vofcomp=(im-1)*ngeom_recon+1
     vcenter(im)=mofdata(vofcomp)
    enddo ! im

    ! sort the volume fractions.  
    im_exclude=0
    call sort_vof(vcenter,im_exclude,sorted_list,nmat)
    im_crit=sorted_list(1) ! vcenter(im_crit) is largest value.
    if ((im_crit.lt.1).or.(im_crit.gt.nmat)) then
       print *,"im_crit invalid"
       stop
    endif

    do im1=1,nmat+1
    do im2=1,nmat
     int_facefrac_cell(im1,im2)=zero
    enddo
    enddo
    do im1=1,nmat
    do im2=1,nmat
    do dir=1,sdim
     int_face_normal_cell(im1,im2,dir)=0.0
     int_centroid_cell(im1,im2,dir)=0.0
    enddo
    enddo
    enddo

    if (vcenter(im_crit).ge.one-FACETOL_DVOL) then
     ! do nothing, there are no internal faces
    else if (vcenter(im_crit).gt.zero) then
     ! normalize the volume fractions so that the sum is 1.
     call make_vfrac_sum_ok_copy( &
       xsten,nhalf,nhalf_box, &
       bfact,dx, &
       tessellate,mofdata,mofdatavalid,nmat,SDIM,3001)

     do im=1,nmat
      is_processed(im)=0
     enddo

     uncaptured_volume=one
     irank=1

     ! index: (im_inside-1)*(nmat+1)+im_outside
     do while ((irank.le.nmat).and.(uncaptured_volume.gt.zero))
      do im=1,nmat
       vofcomp=(im-1)*ngeom_recon+1
       testflag=NINT(mofdatavalid(vofcomp+SDIM+1))
       if (testflag.eq.irank) then

        is_processed(im)=1

        intercept=mofdatavalid(vofcomp+2*SDIM+2)
        uncaptured_volume=uncaptured_volume-mofdatavalid(vofcomp)
        if (uncaptured_volume.lt.FACETOL_DVOL) then
         uncaptured_volume=zero
        endif
        if (uncaptured_volume.gt.zero) then ! we have a valid interface.
         ! x0 is cell center
         ! dist=intercept+slopes dot (x-x0)
         ! perturb interface into the other materials
         mofdatavalid(vofcomp+2*SDIM+2)=intercept+half*FACETOL_DVOL*dx(1)
         shapeflag=0

         call multi_get_volume_grid( &
          tessellate, &
          bfact,dx,xsten,nhalf, &
          mofdatavalid, &
          xsten,nhalf, &
          dummy_tri, &
          multi_volume_offset, &
          multi_cen_offset, &  ! absolute coord.
          multi_area_offset, &
          xtrilistuncapt, &
          nmax,nmax,nmat,SDIM, &
          shapeflag,3)

         mofdatavalid(vofcomp+2*SDIM+2)=intercept

         call multi_get_volume_grid( &
          tessellate, &
          bfact,dx,xsten,nhalf, &
          mofdatavalid, &
          xsten,nhalf, &
          dummy_tri, &
          multi_volume, &
          multi_cen, &  ! absolute coord.
          multi_area, &
          xtrilistuncapt, &
          nmax,nmax,nmat,SDIM, &
          shapeflag,3)

         ! area centroid calculation for (im,im_opp) interface:
         ! multi_cen(im_opp)*multi_volume(im_opp)=
         !   multi_cen_offset(im_opp)*multi_volume_offset(im_opp)+
         !   area_centroid * 
         !   (multi_volume(im_opp)-multi_volume_offset(im_opp))
         if (multi_volume_offset(im).gt.multi_volume(im)) then

          if (multi_area(im).gt.zero) then

           do im_opp=1,nmat
            local_facefrac(im_opp)=zero
           enddo
           total_face=zero

           do im_opp=1,nmat
            if (im_opp.ne.im) then
             if (is_processed(im_opp).eq.0) then
              if (multi_volume(im_opp).ge.multi_volume_offset(im_opp)) then
               local_facefrac(im_opp)= &
                abs(multi_volume(im_opp)-multi_volume_offset(im_opp))
               total_face=total_face+local_facefrac(im_opp)
              else 
               local_facefrac(im_opp)=0.0d0
              endif
             else if (is_processed(im_opp).eq.1) then
              ! do nothing
             else
              print *,"is_processed invalid"
              stop
             endif
            endif
           enddo !im_opp

           if (total_face.gt.zero) then
            do im_opp=1,nmat
             if (im_opp.ne.im) then
              local_facefrac(im_opp)=local_facefrac(im_opp)/total_face
               ! note: int_face_adjust copy (im,im_opp) component into
               !       the (im_opp,im) component.
              int_facefrac_cell(im,im_opp)=local_facefrac(im_opp)

              if (local_facefrac(im_opp).gt.0.0d0) then
                ! vfrac,cen,order,slope,intercept
               vofcomp=(im-1)*ngeom_recon+1
               do dir=1,sdim
                 ! points from im to im_opp
                int_face_normal_cell(im,im_opp,dir)=-mofdata(vofcomp+sdim+1+dir)
                 ! points from im_opp to im
                int_face_normal_cell(im_opp,im,dir)=mofdata(vofcomp+sdim+1+dir)
               enddo
  ! note: if im and im_opp adjoin the im interface, then
  !  it should be impossible for im and im_opp to adjoin the im_opp interface.
  ! dist_to_int_cell(im,im_opp)=distance from im centroid to im interface.
  ! dist_to_int_cell(im_opp,im)=distance from im_opp centroid to im interface.
               call dist_point_line(sdim,nmat,dx,mofdata, &
                im,im,im_opp,dist_to_int_cell, &
                xcentroid_cell, &
                xcenter_cell, &
                xclosest)  ! absolute coord

               delta_volume=multi_volume(im_opp)-multi_volume_offset(im_opp)
               if (delta_volume.gt.0.0d0) then
                LS=intercept
                do dir=1,sdim
                 int_centroid_cell(im,im_opp,dir)= &
                  (multi_cen(dir,im_opp)*multi_volume(im_opp)- &
                   multi_cen_offset(dir,im_opp)* &
                   multi_volume_offset(im_opp))/delta_volume
                 LS=LS-(int_centroid_cell(im,im_opp,dir)- &
                        xsten(0,dir))* &
                       int_face_normal_cell(im,im_opp,dir) ! -n
                enddo ! dir=1..sdim
                do dir=1,sdim
                 int_centroid_cell(im,im_opp,dir)= &
                   int_centroid_cell(im,im_opp,dir)+ &
                   LS*int_face_normal_cell(im,im_opp,dir) ! -LS * n
                 int_centroid_cell(im_opp,im,dir)= &
                       int_centroid_cell(im,im_opp,dir)
                enddo ! dir=1..sdim
               else
                print *,"delta_volume invalid"
                stop
               endif
              else if (local_facefrac(im_opp).eq.0.0d0) then
               ! do nothing
              else
               print *,"local_facefrac invalid"
               stop
              endif 
             endif ! im_opp<>im
            enddo ! im_opp
            int_facefrac_cell(nmat+1,im)=multi_area(im)
           else
            print *,"opposite materials disappeared"
            stop
           endif
          else
           print *,"im boundary disappeared"
           stop
          endif
         else
          print *,"im region should grow"
          stop
         endif

        endif ! uncaptured_volume>0
       endif  ! testflag=irank
      enddo ! im
      irank=irank+1
     enddo  ! while irank<=nmat and uncaptured_volume>0 

     if ((VP_i_debug.eq.VP_i_current).and. &
         (VP_j_debug.eq.VP_j_current)) then
      print *,"i,j (ptb_int) ",VP_i_current,VP_j_current
      do im1=1,nmat+1
      do im2=1,nmat
       if (int_facefrac_cell(im1,im2).gt.0.0d0) then
        print *,"im1,im2,int_facefrac_cell ",im1,im2, &
         int_facefrac_cell(im1,im2)
        if (im1.le.nmat) then
         do dir=1,sdim
          print *,"im1,im2,dir,int_centroid_cell ",im1,im2,dir, &
           int_centroid_cell(im1,im2,dir)
         enddo
        endif
       endif
      enddo
      enddo
     endif
    
    else
     print *,"vcenter(im_crit) out of range"
     stop
    endif

    return
  end subroutine ptb_int
  !-----------------------------------------------------------
  !///////////////////////////////////////////////////////////////////////

  !==================================================================
  !===================================================================
  !  calculate the gradient and diffusion operator 
  ! 
  !==================================================================
  !===================================================================

  ! -div alpha grad u
  subroutine cell_div_cal_complicated( &
   ngeom_recon, &
   sdim,nmat,dx, &
   xsten, &
   mofdata_sten, &
   rho_box, &
   alpha, &
   mat_cen_sten,&
   im_in,  &
   frac_pair_cell, &
   int_face, &
   int_face_normal, &
   dist_to_int, &
   div_tot)

    implicit none 

    integer,intent(in)       :: sdim,nmat,ngeom_recon
    real(kind=8),intent(in)  :: rho_box(-1:1,-1:1,nmat)
    real(kind=8)             :: rho(-1:1,-1:1,nmat)
    real(kind=8),intent(in)  :: xsten(-3:3,sdim)

    real(kind=8),intent(in)  :: frac_pair_cell(nmat,nmat,sdim,2)
    real(kind=8),intent(in)  :: int_face(nmat,nmat)
    real(kind=8),intent(in)  :: int_face_normal(nmat,nmat,sdim)
    real(kind=8),intent(in)  :: mofdata_sten(-1:1,-1:1,ngeom_recon*nmat)
    real(kind=8),intent(in)  :: dx(sdim)
    real(kind=8),intent(in)  :: alpha(nmat)
    real(kind=8),intent(in)  :: dist_to_int(nmat,nmat)

    integer                  :: i,j,im1,im2
    integer                  :: dir,side
    integer                  :: dir2
    real(kind=8)             :: mat_cen_sten(-1:1,-1:1,nmat,sdim)

    real(kind=8)             :: kappa
    real(kind=8)             :: grad
    real(kind=8),intent(out) :: div_tot
    integer                  :: success_flag


    !integer,intent(in)       :: function_flag
    integer,intent(in)       :: im_in
    real(kind=8)             :: Ltemp
    real(kind=8)             :: coef
    real(kind=8)             :: AFRAC,vf1,vf2
    real(kind=8)             :: LOWTOL,grad_high
    real(kind=8)             :: nf(sdim),n1(sdim)

    real(kind=8),external    :: norm_2d
    type(polygon)            :: cell(-1:1,-1:1)
    integer                  :: vofcomp1,vofcomp2,debughigh
    integer                  :: faceid,ii,jj

    debughigh=0

    LOWTOL=0.01d0

    if (ngeom_recon.ne.2*sdim+3) then
     print *,"ngeom_recon invalid"
     stop
    endif

    div_tot = 0.0d0
    rho = rho_box
    success_flag=1 ! cell_div_cal_complicated

    do i=-1,1
    do j=-1,1
     cell(i,j)%center%val(1)=xsten(2*i,1) 
     cell(i,j)%center%val(2)=xsten(2*j,2) 
    enddo
    enddo

    ! cell boundary(external flux) first

    if(im_in .le. 0) then
      print *,"wrong flag"
      stop
    elseif(im_in .gt. 0) then
     im1 = im_in

     vofcomp1=(im1-1)*ngeom_recon+1
     vf1=mofdata_sten(0,0,vofcomp1)
     if ((vf1.lt.-LOWTOL).or.(vf1.gt.one+LOWTOL)) then
      print *,"vf1 invalid"
      stop
     endif
        
     do im2 = 1,nmat
      do dir = 1,sdim
       do side = 1,2

        AFRAC=frac_pair_cell(im2,im1,dir,side)

        if (abs(AFRAC).le.FACETOL_DVOL) then
         AFRAC=0.0d0
        else if (abs(AFRAC-1.0d0).le.FACETOL_DVOL) then
         AFRAC=1.0d0
        else if ((AFRAC.gt.0.0).and.(AFRAC.lt.1.0)) then
         ! do nothing
        else
         print *,"AFRAC invalid"
         stop
        endif

        if ((AFRAC.gt.0.0d0).and.(AFRAC.le.1.0d0)) then

         if((dir.eq.1).and.(side.eq.1)) then
          faceid=1
          ii=-1
          jj=0
         elseif((dir.eq.1).and.(side.eq.2)) then
          faceid=2
          ii=1
          jj=0
         elseif((dir.eq.2).and.(side.eq.1)) then
          faceid=3
          ii=0
          jj=-1
         elseif((dir.eq.2).and.(side.eq.2)) then
          faceid=4
          ii=0
          jj=1
         else
          print *,"dir or side invalid"
          stop
         endif

         nf(1)=ii
         nf(2)=jj

          ! success_flag=0 => failure
          ! success_flag=1 => success

         Ltemp=0.0d0
         do dir2=1,sdim
          n1(dir2) = mat_cen_sten(ii,jj,im2,dir2)-mat_cen_sten(0,0,im1,dir2)
          Ltemp=Ltemp+(mat_cen_sten(ii,jj,im2,dir2)- &
                       mat_cen_sten(0,0,im1,dir2))**2
         enddo
         Ltemp=sqrt(Ltemp)

         call get_kappa_complicated( &
          alpha,xsten,dir,side, &
          mat_cen_sten(0,0,im1,:), &
          mat_cen_sten(ii,jj,im2,:), &
          im1,im2,kappa,nmat,sdim)
         grad=rho(0,0,im1)-rho(ii,jj,im2)

         if (Ltemp.le.0.0) then
          print *,"Ltemp invalid"
          stop
         endif

         coef=kappa/Ltemp
         coef=coef*abs(dot_product(nf,n1))/Ltemp ! |n1/Ltemp|=1 

         grad=-grad*coef

         vofcomp2=(im2-1)*ngeom_recon+1
         vf2=mofdata_sten(ii,jj,vofcomp2)
         if ((vf2.lt.-LOWTOL).or.(vf2.gt.one+LOWTOL)) then
          print *,"vf2 invalid"
          stop
         endif
         if  ((vf1.le.1.0d0-LOWTOL).or.(vf2.le.1.0d0-LOWTOL)) then
          call ext_grad(nmat,sdim,faceid,alpha,im2,im1,rho, &
           cell,dx,mat_cen_sten,mofdata_sten, &
           grad_high,success_flag)
          if (success_flag.eq.1) then
           if (debughigh.eq.1) then
            print *,"high order dir,side,x,y,im1,im2 ",dir,side, &
             xsten(0,1),xsten(0,2),im1,im2
           endif
           grad=grad_high
          else if (success_flag.eq.0) then
           ! do nothing
          else
           print *,"success_flag invalid"
           stop
          endif
         else if ((abs(vf1-1.0d0).le.LOWTOL).and. &
                  (abs(vf2-1.0d0).le.LOWTOL)) then
          ! do nothing
         else
          print *,"vf1 or vf2 invalid"
          stop
         endif

          ! div_tot will be an approximation to 
          ! -integral_material_boundary alpha grad rho dot n dA
          ! where n points out.
         div_tot = div_tot-grad*AFRAC*dx(dir)

        elseif(AFRAC.gt.1.0d0)then
         print *,"frac_pair_cell bust in cell_div_cal_complicated"
         stop
        elseif(AFRAC.lt.0.0d0) then
         print *,"frac_pair_cell invalid in cell_div_cal_complicated"
         stop
        endif

       enddo ! side
      enddo ! dir
     enddo ! im2

      ! internal flux cal
     do im2 = 1,nmat
      AFRAC=int_face(im1,im2)
      if(AFRAC.gt.FACETOL_DVOL*dx(1)) then
       if (im1.eq.im2) then
        print *,"im1==im2 cannot happen"
        stop
       endif

       ! nf points from im1 to im2.
       Ltemp=0.0d0
       do dir2=1,sdim
        nf(dir2)=int_face_normal(im1,im2,dir2)
        Ltemp=Ltemp+(mat_cen_sten(0,0,im2,dir2)- &
                     mat_cen_sten(0,0,im1,dir2))**2
       enddo
       Ltemp=sqrt(Ltemp)

       do dir2=1,sdim
        n1(dir2) = mat_cen_sten(0,0,im2,dir2)-mat_cen_sten(0,0,im1,dir2)
       enddo
       call get_kappa_int_complicated(alpha,dist_to_int,im1,im2,kappa,nmat)

       coef=kappa/Ltemp
       coef=coef*abs(dot_product(nf,n1))/Ltemp ! |n1/Ltemp|=1

       grad=-coef*(rho(0,0,im1)-rho(0,0,im2))

       call int_grad(nmat,sdim,alpha,im2,im1,rho,&
        cell,dx,mat_cen_sten,mofdata_sten, &
        grad_high,success_flag)
       if (success_flag.eq.1) then
        if (debughigh.eq.1) then
         print *,"high order interior,x,y,im1,im2 ", &
           xsten(0,1),xsten(0,2),im1,im2
        endif
        grad=grad_high
       else if (success_flag.eq.0) then
        ! do nothing
       else
        print *,"success_flag invalid"
        stop
       endif
 
       div_tot = div_tot - grad*AFRAC
      endif
     enddo ! im2
       
    else
     print *,"material_in invalid" 
     stop
    endif

  end subroutine cell_div_cal_complicated

  subroutine check_for_least_squares( &
     sdim, &
     nmat, &
     alpha, &
     xclosest, & ! absolute coordinates
     x0, & 
     dist_to_int, &
     int_face, &
     im_main, &
     im_tilde, &
     facearea)
  IMPLICIT NONE

  integer, intent(in)  :: sdim
  integer, intent(in)  :: nmat
  integer, intent(in)  :: im_main
  integer, intent(out) :: im_tilde
  real(kind=8), intent(out) :: facearea
  real(kind=8), intent(in)  :: alpha(nmat)
  real(kind=8), intent(in)  :: int_face(nmat,nmat)
   ! dist_to_int(im1,im2) is distance from im1 centroid to im1,im2
   ! interface.
  real(kind=8), intent(in)  :: dist_to_int(nmat,nmat)
  real(kind=8), intent(in)  :: xclosest(nmat,nmat,sdim) ! absolute coord
  real(kind=8), intent(in)  :: x0(sdim)
  integer :: im_multimaterial
  integer :: im_opp
  integer :: internal_multimaterial
  integer :: im_source,im_dest

  im_tilde=im_main
  im_multimaterial=0
  facearea=zero

  do im_opp=1,nmat
   if (im_opp.eq.im_main) then
    ! do nothing
   else if (im_opp.ne.im_main) then
    if (dist_to_int(im_main,im_opp).eq.zero) then
     ! do nothing
    else if (dist_to_int(im_main,im_opp).gt.zero) then
     if (im_multimaterial.eq.0) then
      im_multimaterial=im_opp
     else if ((im_multimaterial.ge.1).and. &
              (im_multimaterial.le.nmat)) then
      if (dist_to_int(im_main,im_opp).lt. &
          dist_to_int(im_main,im_multimaterial)) then
       im_multimaterial=im_opp
      endif
     else
      print *,"im_multimaterial invalid"
      stop
     endif
    else 
     print *,"dist_to_int(im_main,im_opp) invalid" 
     stop
    endif

    if (im_multimaterial.eq.0) then
     ! do nothing
    else if ((im_multimaterial.ge.1).and. &
             (im_multimaterial.le.nmat).and. &
             (im_multimaterial.ne.im_main)) then
     call init_multimaterial_flag(im_main,im_multimaterial, &
             nmat,internal_multimaterial,im_source,im_dest)
     if (internal_multimaterial.eq.0) then 
      print *,"cannot have im_multimaterial==im_main"
      stop
     else if ((internal_multimaterial.ge.1).and. &
               (internal_multimaterial.le.2*global_nten)) then
      im_tilde=im_multimaterial
      facearea=int_face(im_main,im_multimaterial)
     else if ((internal_multimaterial.le.-1).and. &
              (internal_multimaterial.ge.-global_nten)) then
      im_tilde=im_multimaterial
      facearea=int_face(im_main,im_multimaterial)
     else
      print *,"internal_multimaterial invalid"
      stop
     endif
    else
     print *,"im_multimaterial invalid"
     stop
    endif
   else
    print *,"im_opp invalid"
    stop
   endif
  enddo ! im_opp=1..nmat

  return
  end subroutine check_for_least_squares

  subroutine VP_sanity_check(x_interface,LS_sub,xsten)
  implicit none 
  procedure(sub_interface) :: LS_sub

  REAL_T, intent(in) :: x_interface(2)
  REAL_T, intent(in) :: xsten(-3:3,2)
  REAL_T dist_sanity
  REAL_T dx_local
  INTEGER_T dir_local
  REAL_T VP_tol
  INTEGER_T in_box

  VP_tol=1.0D-5

  in_box=1
  do dir_local=1,2
   dx_local=xsten(1,dir_local)-xsten(-1,dir_local)
   if (dx_local.gt.0.0d0) then
    if ((x_interface(dir_local).le.xsten(-1,dir_local)-VP_tol*dx_local).or. &
        (x_interface(dir_local).ge.xsten(1,dir_local)+VP_tol*dx_local)) then
     if (VP_areaface_current.gt.0.0d0) then
      if (VP_areaface_current.le.0.1d0*dx_local) then
       in_box=0
      else if (VP_areaface_current.ge.0.1d0*dx_local) then
       ! do nothing
      else
       print *,"VP_areaface_current invalid"
       stop
      endif
     else
      print *,"VP_areaface_current invalid"
      stop
     endif
    else if ((x_interface(dir_local).ge. &
              xsten(-1,dir_local)-VP_tol*dx_local).and. &
             (x_interface(dir_local).le. &
              xsten(1,dir_local)+VP_tol*dx_local)) then
     ! do nothing
    else
     print *,"x_interface invalid"
     stop
    endif
   else
    print *,"dx_local invalid"
    stop
   endif
  enddo ! dir_local=1..2

  if (in_box.eq.1) then
   call LS_sub(VP_im_in_current,x_interface(1),x_interface(2), &
          dist_sanity,probtype)
   if (abs(dist_sanity).gt.VP_max_LS_error) then
    VP_i_max=VP_i_current
    VP_j_max=VP_j_current
    VP_dir_max=VP_dir_current
    VP_side_max=VP_side_current
    VP_im_in_max=VP_im_in_current
    VP_im_out_max=VP_im_out_current
    VP_max_LS_error=abs(dist_sanity)
    do dir_local=1,2
     VP_centroid_face(dir_local)=x_interface(dir_local)
    enddo
    VP_areaface_max=VP_areaface_current
    VP_test_number_max=VP_test_number_current
   endif
  else if (in_box.eq.0) then
   ! do nothing
  else
   print *,"in_box invalid"
   stop
  endif
  
  return
  end subroutine VP_sanity_check


  subroutine VP_inside_check(x_interface,xsten)
  implicit none 
  procedure(sub_interface) :: LS_sub

  REAL_T, intent(in) :: x_interface(2)
  REAL_T, intent(in) :: xsten(-3:3,2)
  REAL_T dx_local(2)
  INTEGER_T dir_local
  REAL_T VP_tol
  INTEGER_T in_box

  VP_tol=1.0D-5

  in_box=1
  do dir_local=1,2
   dx_local(dir_local)=xsten(1,dir_local)-xsten(-1,dir_local)
   if (dx_local(dir_local).gt.0.0d0) then
    if ((x_interface(dir_local).le. &
         xsten(-1,dir_local)-VP_tol*dx_local(dir_local)).or. &
        (x_interface(dir_local).ge. &
         xsten(1,dir_local)+VP_tol*dx_local(dir_local))) then
     if (VP_areaface_current.gt.0.0d0) then
      in_box=0
     else
      print *,"VP_areaface_current invalid"
      stop
     endif
    else if ((x_interface(dir_local).ge. &
              xsten(-1,dir_local)-VP_tol*dx_local(dir_local)).and. &
             (x_interface(dir_local).le. &
              xsten(1,dir_local)+VP_tol*dx_local(dir_local))) then
     ! do nothing
    else
     print *,"x_interface invalid"
     stop
    endif
   else
    print *,"dx_local invalid"
    stop
   endif
  enddo ! dir_local=1..2

  if (in_box.eq.1) then
   ! do nothing
  else if (in_box.eq.0) then
   print *,"critical face value not in the box"
   print *,"VP_i_current ",VP_i_current
   print *,"VP_j_current ",VP_j_current
   print *,"x=",xsten(0,1)
   print *,"y=",xsten(0,2)
   print *,"VP_xcentroid ",x_interface(1)
   print *,"VP_ycentroid ",x_interface(2)
   do dir_local=1,2
    print *,"dir_local= ",dir_local
    print *,"dx_local(dir_local) ",dx_local(dir_local)
    print *,"xlocut=",xsten(-1,dir_local)-VP_tol*dx_local(dir_local)
    print *,"xhicut=",xsten(1,dir_local)+VP_tol*dx_local(dir_local)
   enddo
   print *,"VP_areaface_current ",VP_areaface_current
   print *,"VP_dir_current ",VP_dir_current
   print *,"VP_side_current ",VP_side_current
   print *,"VP_im_in_current ",VP_im_in_current
   print *,"VP_im_out_current ",VP_im_out_current
   print *,"VP_test_number_current ",VP_test_number_current
   stop
  else
   print *,"in_box invalid"
   stop
  endif
  
  return
  end subroutine VP_inside_check


  subroutine get_multimaterial_temperature( &
          L1,L2,T1,T2,K1,K2,TI)
  real(kind=8), intent(in)  :: L1,L2,T1,T2,K1,K2
  real(kind=8), intent(out) :: TI

  real(kind=8) :: sigma1,sigma2

  if (L1.eq.0.0d0) then
   TI=T1
  else if (L2.eq.0.0d0) then
   TI=T2
  else if (K1.eq.0.0d0) then
   TI=T2
  else if (K2.eq.0.0d0) then
   TI=T1
  else if ((L1.gt.0.0d0).and.(L2.gt.0.0d0).and. &
           (K1.gt.0.0d0).and.(K2.gt.0.0d0)) then
   sigma1=K1/L1
   sigma2=K2/L2
   TI=(sigma1*T1+sigma2*T2)/(sigma1+sigma2)
  else
   print *,"parameters invalid in get_multimaterial_temperature"
   stop
  endif

  end subroutine get_multimaterial_temperature


   ! div_tot = -div alpha grad u
   ! diag_coeff_flag==1 => find diagonal coeff 
   !   (local_hflag=1, local_rho=0 except for target comp.)
  subroutine cell_div_cal_simple( &
   LS_sub, & 
   diag_coeff_flag, &  
   linear_exact, &
   operator_internal, &
   operator_external, &
   hflag, &
   ngeom_recon, &
   sdim,nmat, &
   time_div, &
   dx, &
   xsten, &
   mofdata_sten, &
   rho_box, &
   alpha, &
   mat_cen_sten,&  ! absolute coordinates
   im_in,  &
   frac_pair_cell, &
   x_pair_cell, & ! absolute coord.
   gap_alarm_sten, &
   int_face_sten, &
   int_centroid_sten, &
   int_face_normal, &
   dist_to_int_sten, &
   xclosest_sten, &  ! absolute coord
   div_tot)
   use global_utility_module, only : crossprod, &
           gradient_at_dirichlet_boundary
   use mass_transfer_module, only : get_interface_temperature

    implicit none 

    procedure(sub_interface) :: LS_sub

    integer,intent(in)       :: sdim
    integer,intent(in)       :: nmat,ngeom_recon
    integer,intent(in)       :: hflag
    integer                  :: local_hflag
    integer,intent(in)       :: linear_exact
    integer                  :: operator_internal
    integer                  :: operator_external
    integer,intent(in)       :: diag_coeff_flag
    real(kind=8),intent(in)  :: rho_box(-1:1,-1:1,nmat)
    real(kind=8)             :: local_rho(-1:1,-1:1,nmat)
    real(kind=8),intent(in)  :: xsten(-3:3,sdim)

      ! im_outside,im_inside
    real(kind=8),intent(in)  :: frac_pair_cell(nmat,nmat,sdim,2)
    real(kind=8),intent(in)  :: x_pair_cell(nmat,nmat,sdim,sdim,2)
    real(kind=8),intent(in)  :: gap_alarm_sten(-1:1,-1:1)
    real(kind=8),intent(in)  :: int_face_sten(-1:1,-1:1,nmat,nmat)
    real(kind=8),intent(in)  :: int_centroid_sten(-1:1,-1:1,nmat,nmat,sdim)
    real(kind=8),intent(in)  :: int_face_normal(nmat,nmat,sdim)
    real(kind=8),intent(in)  :: mofdata_sten(-1:1,-1:1,ngeom_recon*nmat)
    real(kind=8),intent(in)  :: time_div
    real(kind=8),intent(in)  :: dx(sdim)
    real(kind=8),intent(in)  :: alpha(nmat)
    real(kind=8),intent(in)  :: dist_to_int_sten(-1:1,-1:1,nmat,nmat)
     ! xclosest is in absolute coord.
    real(kind=8),intent(in)  :: xclosest_sten(-1:1,-1:1,nmat,nmat,sdim)

    real(kind=8)             :: int_face_inside(nmat,nmat)
    real(kind=8)             :: int_face_outside(nmat,nmat)
    real(kind=8)             :: dist_to_int_inside(nmat,nmat)
    real(kind=8)             :: dist_to_int_outside(nmat,nmat)
    real(kind=8)             :: xclosest_inside(nmat,nmat,sdim)
    real(kind=8)             :: xclosest_outside(nmat,nmat,sdim)
    real(kind=8)             :: x0_inside(sdim)
    real(kind=8)             :: x0_outside(sdim)

    integer                  :: i,j
    integer                  :: im
    integer                  :: im_inside,im_outside
    integer                  :: im1_face,im2_face
    integer                  :: im_tilde
    integer                  :: im_tilde_inside,im_tilde_outside
    integer                  :: dir_main,side
    integer                  :: dir2
    integer                  :: dircrit
    real(kind=8)             :: mat_cen_sten(-1:1,-1:1,nmat,sdim) ! absolute

    real(kind=8)             :: kappa,AFRAC
    real(kind=8)             :: vf1
    real(kind=8)             :: grad
    real(kind=8),intent(out) :: div_tot

    !integer,intent(in)       :: function_flag
    integer,intent(in)       :: im_in
    real(kind=8)             :: LOWTOL
    real(kind=8)             :: nf(sdim)
    real(kind=8)             :: nminus(3)
    real(kind=8)             :: nplus(3)
    integer                  :: vofcomp1
    integer                  :: faceid,ii,jj

    integer                  :: internal_multimaterial
    integer                  :: im_source
    integer                  :: im_dest

    real(kind=8),external    :: norm_2d
    type(polygon)            :: cell(-1:1,-1:1)

    real(kind=8)             :: Tdata(4)
    real(kind=8)             :: Tminus,Tplus
    integer                  :: idata
    real(kind=8)             :: xface(sdim)
    real(kind=8)             :: xinside(sdim)
    real(kind=8)             :: xoutside(sdim)
    real(kind=8)             :: xinside_face(sdim)
    real(kind=8)             :: xoutside_face(sdim)
    real(kind=8)             :: xminus(3)
    real(kind=8)             :: xminus_tilde(3)
    real(kind=8)             :: xplus(3)
    real(kind=8)             :: xplus_tilde(3)
    real(kind=8)             :: xminusI(3)
    real(kind=8)             :: xminusI_tilde(3)
    real(kind=8)             :: xplusI(3)
    real(kind=8)             :: xplusI_tilde(3)
    real(kind=8)             :: x2_I(3)
    real(kind=8)             :: x3_I(3)
    real(kind=8)             :: s_I(3)
    real(kind=8)             :: t1_vec(3)
    real(kind=8)             :: t2_vec(3)
    real(kind=8)             :: rho_inside,rho_outside,rho_I
    real(kind=8)             :: mag_ntilde
    real(kind=8)             :: mag_nminus
    real(kind=8)             :: mag_nplus
    real(kind=8)             :: mag1_vec
    real(kind=8)             :: mag2_vec
    real(kind=8)             :: len_inside,len_outside,len_total
    real(kind=8)             :: facearea_inside,facearea_outside
    real(kind=8)             :: facearea_minus,facearea_plus
    integer                  :: nplus_valid
    real(kind=8)             :: gradT_face(3)
    real(kind=8)             :: mag_debug
    real(kind=8)             :: xinside_face_debug(2)
    real(kind=8)             :: Tminus_tilde,Tplus_tilde
    real(kind=8)             :: len_minus,len_minus_tilde
    real(kind=8)             :: len_plus,len_plus_tilde


    LOWTOL=0.01d0

    if (ngeom_recon.ne.2*sdim+3) then
     print *,"ngeom_recon invalid"
     stop
    endif
    if (time_div.ge.0.0d0) then
     ! do nothing
    else
     print *,"time_div invalid"
     stop
    endif

    div_tot = 0.0d0

    local_hflag=hflag

    do i=-1,1
    do j=-1,1
    do im=1,nmat
     local_rho(i,j,im) = rho_box(i,j,im)
    enddo
    enddo
    enddo

    do i=-1,1
    do j=-1,1
     cell(i,j)%center%val(1)=xsten(2*i,1) 
     cell(i,j)%center%val(2)=xsten(2*j,2) 
    enddo
    enddo

    if ((im_in.ge.1).and.(im_in.le.nmat)) then

     if (diag_coeff_flag.eq.1) then
      do i=-1,1
      do j=-1,1
      do im=1,nmat
       local_rho(i,j,im) = zero
      enddo
      enddo
      enddo
      local_rho(0,0,im_in)=one
      local_hflag=1
     else if (diag_coeff_flag.eq.0) then
      ! do nothing
     else
      print *,"diag_coeff_flag invalid"
      stop
     endif

      ! cell boundary(external flux) first
     im_inside = im_in
    
     vofcomp1=(im_inside-1)*ngeom_recon+1
     vf1=mofdata_sten(0,0,vofcomp1)
     if ((vf1.ge.-LOWTOL).and.(vf1.le.one+LOWTOL)) then
      ! do nothing
     else
      print *,"vf1 invalid"
      stop
     endif

     rho_inside=local_rho(0,0,im_inside)
   
     do im1_face=1,nmat
     do im2_face=1,nmat
      int_face_inside(im1_face,im2_face)= &
        int_face_sten(0,0,im1_face,im2_face)
      dist_to_int_inside(im1_face,im2_face)= &
        dist_to_int_sten(0,0,im1_face,im2_face)
      do dir2=1,sdim
       xclosest_inside(im1_face,im2_face,dir2)= &
        xclosest_sten(0,0,im1_face,im2_face,dir2)
      enddo
     enddo ! im2_face
     enddo ! im1_face

     do dir2=1,sdim
      xinside(dir2)=mat_cen_sten(0,0,im_inside,dir2) ! absolute
     enddo ! dir2=1..sdim

     do dir2=1,3
      gradT_face(dir2)=0.0d0
      xminus(dir2)=0.0d0
      xminusI(dir2)=0.0d0
      xplus(dir2)=0.0d0
      xplusI(dir2)=0.0d0
      x2_I(dir2)=0.0d0
      x3_I(dir2)=0.0d0
      nminus(dir2)=0.0d0
      nplus(dir2)=0.0d0
     enddo

     do im_outside = 1,nmat
      do dir_main = 1,sdim
       do side = 1,2

        VP_dir_current=dir_main
        VP_side_current=side
        VP_im_in_current=im_inside
        VP_im_out_current=im_outside

        AFRAC=frac_pair_cell(im_outside,im_inside,dir_main,side)

        if (abs(AFRAC).le.FACETOL_DVOL) then
         AFRAC=0.0d0
        else if (abs(AFRAC-1.0d0).le.FACETOL_DVOL) then
         AFRAC=1.0d0
        else if ((AFRAC.gt.0.0).and.(AFRAC.lt.1.0)) then
         ! do nothing
        else
         print *,"AFRAC invalid"
         stop
        endif

        if ((AFRAC.gt.0.0d0).and.(AFRAC.le.1.0d0)) then

         do dir2=1,sdim
          xface(dir2)=xsten(0,dir2)
          x0_inside(dir2)=xsten(0,dir2)
          x0_outside(dir2)=xsten(0,dir2)
         enddo

         if((dir_main.eq.1).and.(side.eq.1)) then
          faceid=1
          ii=-1
          jj=0
          xface(dir_main)=xsten(ii,dir_main)
          x0_outside(dir_main)=xsten(2*ii,dir_main)
         elseif((dir_main.eq.1).and.(side.eq.2)) then
          faceid=2
          ii=1
          jj=0
          xface(dir_main)=xsten(ii,dir_main)
          x0_outside(dir_main)=xsten(2*ii,dir_main)
         elseif((dir_main.eq.2).and.(side.eq.1)) then
          faceid=3
          ii=0
          jj=-1
          xface(dir_main)=xsten(jj,dir_main)
          x0_outside(dir_main)=xsten(2*jj,dir_main)
         elseif((dir_main.eq.2).and.(side.eq.2)) then
          faceid=4
          ii=0
          jj=1
          xface(dir_main)=xsten(jj,dir_main)
          x0_outside(dir_main)=xsten(2*jj,dir_main)
         else
          print *,"dir_main or side invalid"
          stop
         endif

         nf(1)=ii
         nf(2)=jj

         rho_outside=local_rho(ii,jj,im_outside)

           ! SANITY CHECK
         do dir2=1,sdim
          xinside_face(dir2)=xface(dir2)
          if (dir2.ne.dir_main) then
           xinside_face(dir2)= &
              x_pair_cell(im_outside,im_inside,dir2,dir_main,side)
          endif
         enddo

         VP_areaface_current=AFRAC*dx(1)
         VP_test_number_current=1
         if (im_inside.ne.im_outside) then
          call VP_sanity_check(xinside_face,LS_sub,xsten)
         endif
         call VP_inside_check(xinside_face,xsten)
       
         do dir2=1,sdim
          xoutside(dir2)=mat_cen_sten(ii,jj,im_outside,dir2) ! absolute
          xinside_face(dir2)=xface(dir2)
          xoutside_face(dir2)=xface(dir2)
         enddo

         mag_ntilde=zero
         do dir2=1,sdim
          mag_ntilde=mag_ntilde+ &
           (xinside(dir2)-xoutside(dir2))**2
         enddo
         mag_ntilde=sqrt(mag_ntilde)

         if (mag_ntilde.gt.zero) then

          do im1_face=1,nmat
          do im2_face=1,nmat
           int_face_outside(im1_face,im2_face)= &
             int_face_sten(ii,jj,im1_face,im2_face)
           dist_to_int_outside(im1_face,im2_face)= &
             int_face_sten(ii,jj,im1_face,im2_face)
           do dir2=1,sdim
            xclosest_outside(im1_face,im2_face,dir2)= &
             xclosest_sten(ii,jj,im1_face,im2_face,dir2)
           enddo
          enddo ! im2_face
          enddo ! im1_face

           !   --------------
           !   |     |.     |
           !   |     |      |
           !   |   . |      |
           !   --------------
          if (operator_external.eq.1) then ! simple method
           len_total=abs(xoutside(dir_main)-xinside(dir_main))
           if (len_total.gt.zero) then
            len_inside=abs(xface(dir_main)-xinside(dir_main))
            len_outside=abs(xface(dir_main)-xoutside(dir_main))
            if ((len_inside.gt.zero).and.(len_outside.gt.zero)) then
             if (abs(len_inside+len_outside-len_total).le. &
                 VOFTOL*len_total) then
              len_total=len_inside+len_outside
              do dir2=1,sdim
               if (dir2.ne.dir_main) then
                xinside_face(dir2)= &
                 (len_inside*xoutside(dir2)+ &
                  len_outside*xinside(dir2))/len_total
                xoutside_face(dir2)=xinside_face(dir2)
               endif
              enddo
             else
              print *,"gap_alarm cutoff is 0.1"

              print *,"exterior flux"
              print *,"len_inside=",len_inside
              print *,"len_outside=",len_outside
              print *,"len_total=",len_total
              print *,"abs(len_inside+len_outside-len_total) ", &
                 abs(len_inside+len_outside-len_total)
              print *,"dir_main=",dir_main
              print *,"xinside(dir_main)= ",xinside(dir_main)
              print *,"xoutside(dir_main)= ",xoutside(dir_main)
              print *,"xface(dir_main)= ",xface(dir_main)
              print *,"VOFTOL=",VOFTOL
              print *,"len_inside,len_outside or len_total invalid"
              stop
             endif
            else
             print *,"len_inside or len_outside invalid"
             stop
            endif
           else
            print *,"len_total invalid"
            stop
           endif
          else if (operator_external.eq.2) then ! Dai and Scannapieco
           do dir2=1,sdim
            if (dir2.ne.dir_main) then
             xinside_face(dir2)= &
                 x_pair_cell(im_outside,im_inside,dir2,dir_main,side)
             xoutside_face(dir2)=xinside_face(dir2)
            endif
           enddo
           VP_areaface_current=AFRAC*dx(1)
           VP_test_number_current=2
           if (im_inside.ne.im_outside) then
            call VP_sanity_check(xinside_face,LS_sub,xsten)
           endif
           call VP_inside_check(xinside_face,xsten)
          else if (operator_external.eq.3) then ! orthogonal probe
           do dir2=1,sdim
            if (dir2.ne.dir_main) then
             xinside_face(dir2)=xinside(dir2)
             xoutside_face(dir2)=xoutside(dir2)
            endif
           enddo
           VP_areaface_current=AFRAC*dx(1)
           VP_test_number_current=22
           call VP_inside_check(xinside_face,xsten)
           VP_test_number_current=23
           call VP_inside_check(xoutside_face,xsten)
          else
           print *,"operator_external invalid"
           stop
          endif

          im_tilde=im_inside
          nplus_valid=0
          facearea_minus=0.0d0
          facearea_plus=0.0d0

          do idata=1,4
           Tdata(idata)=0.0d0
          enddo

          if (im_inside.ne.im_outside) then
           ! do nothing
          else if (linear_exact.eq.0) then
           ! do nothing
          else if (gap_alarm_sten(0,0).eq.1.0d0) then
           ! do nothing
          else if (gap_alarm_sten(ii,jj).eq.1.0d0) then
           ! do nothing
          else if ((im_inside.eq.im_outside).and. &
                   (linear_exact.eq.1).and. &
                   (gap_alarm_sten(0,0).eq.0.0d0).and. &
                   (gap_alarm_sten(ii,jj).eq.0.0d0)) then
             ! if a valid opposite material is found, then
             ! im_tilde!=im_inside and facearea>0
           call check_for_least_squares( &
             sdim, &
             nmat, &
             alpha, &
             xclosest_inside, &
             x0_inside, &
             dist_to_int_inside, &
             int_face_inside, &
             im_inside, &
             im_tilde_inside, &
             facearea_inside)
           call check_for_least_squares( &
             sdim, &
             nmat, &
             alpha, &
             xclosest_outside, &
             x0_outside, &
             dist_to_int_outside, &
             int_face_outside, &
             im_inside, &
             im_tilde_outside, &
             facearea_outside)

           if ((facearea_inside.eq.zero).and. &
               (facearea_outside.eq.zero)) then

            ! do nothing

           else if (facearea_inside.ge.facearea_outside) then

            if (facearea_inside.gt.zero) then
             ! do nothing
            else
             print *,"facearea_inside invalid"
             stop
            endif

            if ((im_tilde_inside.ge.1).and. &
                (im_tilde_inside.le.nmat).and. &
                (im_tilde_inside.ne.im_inside)) then
             ! do nothing
            else
             print *,"im_tilde invalid"
             stop
            endif
            if ((im_tilde_inside.eq.im_tilde_outside).or. &
                (im_tilde_outside.eq.im_inside)) then
             im_tilde=im_tilde_inside
             Tminus=rho_inside
             Tplus=rho_outside

             if (im_tilde_outside.eq.im_inside) then
              ! do nothing
             else if (im_tilde_outside.eq.im_tilde) then
              nplus_valid=1
             else
              print *,"im_tilde_outside invalid"
              stop
             endif

             Tminus_tilde=local_rho(0,0,im_tilde) ! inside
             Tplus_tilde=local_rho(ii,jj,im_tilde) ! outside

             len_minus=0.0d0
             len_minus_tilde=0.0d0

             len_plus=0.0d0
             len_plus_tilde=0.0d0

             do dir2=1,sdim
              xminus(dir2)=xinside(dir2)
              xplus(dir2)=xoutside(dir2)
              xminusI(dir2)=xclosest_inside(im_inside,im_tilde,dir2)
              xminusI_tilde(dir2)=xclosest_inside(im_tilde,im_inside,dir2)
              xminus_tilde(dir2)=mat_cen_sten(0,0,im_tilde,dir2) ! absolute
              len_minus=len_minus+(xminus(dir2)-xminusI(dir2))**2
              len_minus_tilde=len_minus_tilde+ &
                    (xminus_tilde(dir2)-xminusI_tilde(dir2))**2
              if (nplus_valid.eq.0) then
               xplusI(dir2)=zero
               xplusI_tilde(dir2)=zero
              else if (nplus_valid.eq.1) then
               xplusI(dir2)=xclosest_outside(im_inside,im_tilde,dir2)
               xplusI_tilde(dir2)=xclosest_outside(im_tilde,im_inside,dir2)
               xplus_tilde(dir2)=mat_cen_sten(ii,jj,im_tilde,dir2)
               len_plus=len_plus+(xplus(dir2)-xplusI(dir2))**2
               len_plus_tilde=len_plus_tilde+ &
                    (xplus_tilde(dir2)-xplusI_tilde(dir2))**2
              else
               print *,"nplus_valid invalid"
               stop
              endif
             enddo ! dir2=1..sdim

             len_minus=sqrt(len_minus)
             len_minus_tilde=sqrt(len_minus_tilde)

             len_plus=sqrt(len_plus)
             len_plus_tilde=sqrt(len_plus_tilde)

             facearea_minus=facearea_inside
             facearea_plus=facearea_outside

            else if ((im_tilde_inside.ne.im_tilde_outside).and. &
                     (im_tilde_outside.ne.im_inside)) then
             ! do nothing
            else
             print *,"im_tilde bust"
             stop
            endif

           else if (facearea_inside.le.facearea_outside) then

            if (facearea_outside.gt.zero) then
             ! do nothing
            else
             print *,"facearea_outside invalid"
             stop
            endif

            if ((im_tilde_outside.ge.1).and. &
                (im_tilde_outside.le.nmat).and. &
                (im_tilde_outside.ne.im_inside)) then
             ! do nothing
            else
             print *,"im_tilde_outside invalid"
             stop
            endif
            if ((im_tilde_inside.eq.im_tilde_outside).or. &
                (im_tilde_inside.eq.im_inside)) then
             im_tilde=im_tilde_outside
             Tplus=rho_inside
             Tminus=rho_outside

             if (im_tilde_inside.eq.im_inside) then
              ! do nothing
             else if (im_tilde_inside.eq.im_tilde) then
              nplus_valid=1
             else
              print *,"im_tilde_outside invalid"
              stop
             endif

             Tminus_tilde=local_rho(ii,jj,im_tilde) ! outside
             Tplus_tilde=local_rho(0,0,im_tilde) ! inside

             len_minus=0.0d0
             len_minus_tilde=0.0d0

             len_plus=0.0d0
             len_plus_tilde=0.0d0

             do dir2=1,sdim
              xplus(dir2)=xinside(dir2)
              xminus(dir2)=xoutside(dir2)
              xminusI(dir2)=xclosest_outside(im_inside,im_tilde,dir2)
              xminusI_tilde(dir2)=xclosest_outside(im_tilde,im_inside,dir2)
              xminus_tilde(dir2)=mat_cen_sten(ii,jj,im_tilde,dir2) ! absolute
              len_minus=len_minus+(xminus(dir2)-xminusI(dir2))**2
              len_minus_tilde=len_minus_tilde+ &
                    (xminus_tilde(dir2)-xminusI_tilde(dir2))**2
              if (nplus_valid.eq.0) then
               xplusI(dir2)=zero
               xplusI_tilde(dir2)=zero
              else if (nplus_valid.eq.1) then
               xplusI(dir2)=xclosest_inside(im_inside,im_tilde,dir2)
               xplusI_tilde(dir2)=xclosest_inside(im_tilde,im_inside,dir2)
               xplus_tilde(dir2)=mat_cen_sten(0,0,im_tilde,dir2)
               len_plus=len_plus+(xplus(dir2)-xplusI(dir2))**2
               len_plus_tilde=len_plus_tilde+ &
                    (xplus_tilde(dir2)-xplusI_tilde(dir2))**2
              else
               print *,"nplus_valid invalid"
               stop
              endif
             enddo ! dir2=1..sdim

             len_minus=sqrt(len_minus)
             len_minus_tilde=sqrt(len_minus_tilde)

             len_plus=sqrt(len_plus)
             len_plus_tilde=sqrt(len_plus_tilde)

             facearea_minus=facearea_outside
             facearea_plus=facearea_inside

            else if ((im_tilde_inside.ne.im_tilde_outside).and. &
                     (im_tilde_inside.ne.im_inside)) then
             ! do nothing
            else
             print *,"im_tilde bust"
             stop
            endif
           else
            print *,"facearea invalid"
            stop
           endif
          else
           print *,"im_inside or linear_exact invalid"
           stop
          endif

          if ((im_tilde.ne.im_inside).and. &
              (im_tilde.ge.1).and. &
              (im_tilde.le.nmat)) then
              ! input: im_inside, im_tilde, nmat
              ! output: internal_multimaterial,im_source,im_dest
           call init_multimaterial_flag(im_inside,im_tilde,nmat, &
            internal_multimaterial,im_source,im_dest)
           if (internal_multimaterial.eq.0) then 
            print *,"expecting im_inside,im_tilde a Multimaterial boundary"
            stop
           else if (((internal_multimaterial.ge.1).and. &
                     (internal_multimaterial.le.2*global_nten)).or. &
                    ((internal_multimaterial.le.-1).and. &
                     (internal_multimaterial.ge.-global_nten))) then
            mag_nminus=zero
            mag_nplus=zero
            dircrit=0
            do dir2=1,3
             nminus(dir2)=zero
             nplus(dir2)=zero
            enddo
            do dir2=1,sdim
             nminus(dir2)=xminus(dir2)-xminusI(dir2)
             mag_nminus=mag_nminus+nminus(dir2)**2
             if (nplus_valid.eq.0) then
              ! do nothing
             else if (nplus_valid.eq.1) then
              nplus(dir2)=xplus(dir2)-xplusI(dir2)
              mag_nplus=mag_nplus+nplus(dir2)**2
             else
              print *,"nplus_valid invalid"
              stop
             endif
             if (dircrit.eq.0) then
              dircrit=1
             else if (abs(nminus(dircrit)).le.abs(nminus(dir2))) then
              dircrit=dir2
             else if (abs(nminus(dircrit)).ge.abs(nminus(dir2))) then
              ! do nothing
             else
              print *,"nminus invalid"
              stop
             endif
            enddo ! dir2=1..sdim
            mag_nminus=sqrt(mag_nminus)
            mag_nplus=sqrt(mag_nplus)

            if ((mag_nminus.gt.zero).and. &
                (dircrit.ge.1).and. &
                (dircrit.le.sdim)) then
             do dir2=1,3
              nminus(dir2)=nminus(dir2)/mag_nminus
              if (nplus_valid.eq.0) then
               ! do nothing
              else if (nplus_valid.eq.1) then
               if (mag_nplus.gt.zero) then
                nplus(dir2)=nplus(dir2)/mag_nplus
               else
                print *,"mag_nplus invalid"
                stop
               endif
              else
               print *,"nplus_valid invalid"
               stop
              endif

              s_I(dir2)=one
             enddo ! dir2=1..3
             s_I(dircrit)=zero
             call crossprod(nminus,s_I,t1_vec)
             mag1_vec=zero
             do dir2=1,3
              mag1_vec=mag1_vec+t1_vec(dir2)**2
             enddo
             mag1_vec=sqrt(mag1_vec)
             if (mag1_vec.gt.zero) then
              do dir2=1,3
               t1_vec(dir2)=t1_vec(dir2)/mag1_vec
               x2_I(dir2)=xminusI(dir2)+TANGENT_EPS*mag_ntilde*t1_vec(dir2)
              enddo
              call crossprod(nminus,t1_vec,t2_vec)
              mag2_vec=zero
              do dir2=1,3
               mag2_vec=mag2_vec+t2_vec(dir2)**2
              enddo
              mag2_vec=sqrt(mag2_vec)
              ! the magnitude of the cross product of two orthogonal unit
              ! vectors should be one.
              if (abs(mag2_vec-one).le.VOFTOL) then
               do dir2=1,3
                t2_vec(dir2)=t2_vec(dir2)/mag2_vec
                x3_I(dir2)=xminusI(dir2)+TANGENT_EPS*mag_ntilde*t2_vec(dir2)
               enddo

               if ((local_hflag.eq.1).or. &
                   (internal_multimaterial.lt.0)) then
               
                if ((internal_multimaterial.le.-1).and. &
                    (internal_multimaterial.ge.-global_nten)) then
                 call get_multimaterial_temperature(len_minus, &
                        len_minus_tilde,Tminus,Tminus_tilde, &
                        alpha(im_inside),alpha(im_tilde),Tdata(3))
                 Tdata(1)=Tdata(3)
                 Tdata(2)=Tdata(3)
                 if (nplus_valid.eq.0) then
                  Tdata(4)=zero
                  facearea_plus=zero
                 else if (nplus_valid.eq.1) then 
                  call get_multimaterial_temperature(len_plus, &
                         len_plus_tilde,Tplus,Tplus_tilde, &
                         alpha(im_inside),alpha(im_tilde),Tdata(4))
                 else
                  print *,"nplus_valid invalid"
                  stop
                 endif

                else if ((internal_multimaterial.ge.1).and. &
                         (internal_multimaterial.le.2*global_nten)) then
                 if (local_hflag.eq.1) then
                  ! do nothing
                 else
                  print *,"local_hflag invalid"
                  stop
                 endif
                else
                 print *,"internal_multimaterial invalid"
                 stop
                endif 

               else if ((local_hflag.eq.0).and. &
                        (internal_multimaterial.ge.1).and. &
                        (internal_multimaterial.le.2*global_nten)) then
                  ! in: MASS_TRANSFER_3D.F90
                call get_interface_temperature( &
                  use_tsatfab, &
                  i_tsat,j_tsat,k_tsat, &
                  ireverse_tsat, &
                  iten_tsat, &
                  ntsat, &
                  bfact_tsat, &
                  level_tsat, &
                  finest_level_tsat, &
                  dx_tsat,xlo_tsat, &
                  ngrow_tsat, &
                  fablo_tsat,fabhi_tsat, &
                  TSATFAB,DIMS(TSATFAB), &
                  Tdata(1), &
                  internal_multimaterial, &
                  saturation_temp, &
                  use_exact_temperature, &
                  x2_I, &
                  time_div, &
                  nmat, &
                  global_nten,1)

                call get_interface_temperature( &
                  use_tsatfab, &
                  i_tsat,j_tsat,k_tsat, &
                  ireverse_tsat, &
                  iten_tsat, &
                  ntsat, &
                  bfact_tsat, &
                  level_tsat, &
                  finest_level_tsat, &
                  dx_tsat,xlo_tsat, &
                  ngrow_tsat, &
                  fablo_tsat,fabhi_tsat, &
                  TSATFAB,DIMS(TSATFAB), &
                  Tdata(2), &
                  internal_multimaterial, &
                  saturation_temp, &
                  use_exact_temperature, &
                  x3_I, &
                  time_div, &
                  nmat, &
                  global_nten,2)

                call get_interface_temperature( &
                  use_tsatfab, &
                  i_tsat,j_tsat,k_tsat, &
                  ireverse_tsat, &
                  iten_tsat, &
                  ntsat, &
                  bfact_tsat, &
                  level_tsat, &
                  finest_level_tsat, &
                  dx_tsat,xlo_tsat, &
                  ngrow_tsat, &
                  fablo_tsat,fabhi_tsat, &
                  TSATFAB,DIMS(TSATFAB), &
                  Tdata(3), &
                  internal_multimaterial, &
                  saturation_temp, &
                  use_exact_temperature, &
                  xminusI, &
                  time_div, &
                  nmat, &
                  global_nten,3)

                if (nplus_valid.eq.0) then
                 Tdata(4)=zero
                 facearea_plus=zero
                else if (nplus_valid.eq.1) then 
                 call get_interface_temperature( &
                  use_tsatfab, &
                  i_tsat,j_tsat,k_tsat, &
                  ireverse_tsat, &
                  iten_tsat, &
                  ntsat, &
                  bfact_tsat, &
                  level_tsat, &
                  finest_level_tsat, &
                  dx_tsat,xlo_tsat, &
                  ngrow_tsat, &
                  fablo_tsat,fabhi_tsat, &
                  TSATFAB,DIMS(TSATFAB), &
                  Tdata(4), &
                  internal_multimaterial, &
                  saturation_temp, &
                  use_exact_temperature, &
                  xplusI, &
                  time_div, &
                  nmat, &
                  global_nten,4)
                else
                 print *,"nplus_valid invalid"
                 stop
                endif

               else
                print *,"local_hflag or internal_multimaterial invalid"
                stop
               endif

               if (1.eq.0) then
                if (probtype.eq.403) then
                 if (local_hflag.eq.1) then
                  if ((Tdata(1).eq.zero).and. &
                      (Tdata(2).eq.zero).and. &
                      (Tdata(3).eq.zero).and. &
                      (Tdata(4).eq.zero)) then
                   ! do nothing
                  else
                   print *,"expecting homogeneous Tdata"
                   stop
                  endif
                 else if (local_hflag.eq.0) then
                  print *,"i_tsat,j_tsat,k_tsat ", &
                     i_tsat,j_tsat,k_tsat
                  print *,"im_inside,im_outside ",im_inside,im_outside
                  print *,"dir_main,side ",dir_main,side
                  print *,"Tdata1-4 ",Tdata(1),Tdata(2),Tdata(3),Tdata(4)
                  print *,"nplus_valid=",nplus_valid
                 else
                  print *,"local_hflag invalid"
                  stop
                 endif 
                else
                 print *,"probtype invalid not 403"
                 stop
                endif
               endif

               VP_areaface_current=facearea_minus
               VP_test_number_current=3
               call VP_sanity_check(xminusI,LS_sub,xsten)

               VP_areaface_current=facearea_minus
               VP_test_number_current=4
               call VP_sanity_check(x2_I,LS_sub,xsten)

               VP_areaface_current=facearea_minus
               VP_test_number_current=5
               call VP_sanity_check(x3_I,LS_sub,xsten)

               if (nplus_valid.eq.1) then 
                VP_areaface_current=facearea_plus
                VP_test_number_current=6
                call VP_sanity_check(xplusI,LS_sub,xsten)
               endif

               call gradient_at_dirichlet_boundary( &
                gradT_face, &
                facearea_minus, &
                facearea_plus, &
                xminus,xplus, &
                Tminus,Tplus, &
                nminus,xminusI, &
                nplus,xplusI, &
                x2_I, &
                x3_I, &
                Tdata)
              else
               print *,"mag2_vec invalid"
               stop
              endif
             else
              print *,"mag1_vec invalid"
              stop
             endif
            else
             print *,"mag_nminus invalid"
             stop
            endif
           else
            print *,"internal_multimaterial invalid"
            stop
           endif
          else if (im_tilde.eq.im_inside) then
           ! do nothing
          else
           print *,"im_tilde invalid"
           stop
          endif
         else
          print *,"mag_ntilde invalid: centroids of neighboring cells"
          print *,"cannot coincide at the same point."
          stop
         endif

         if (im_inside.ne.im_outside) then
          call init_multimaterial_flag(im_inside,im_outside,nmat, &
           internal_multimaterial, &
           im_source,im_dest)
          if ((internal_multimaterial.le.-1).and. &
              (internal_multimaterial.ge.-global_nten)) then 
           call get_kappa_simple( &
             sdim, &
             alpha, &
             nmat, &
             im_inside,im_outside, &
             kappa, &
             nf, &  ! points out of the cell
             xinside,xoutside, &
             xinside_face,xoutside_face)

            ! grad approximates: k grad T dot (-nf)
            ! later: div_tot=div_tot+grad*AFRAC*dx(dir_main)
           if (kappa.ge.0.0d0) then
            grad=kappa*(rho_inside-rho_outside)
           else
            print *,"kappa invalid"
            stop
           endif
          else if ((internal_multimaterial.ge.1).and. &
                   (internal_multimaterial.le.2*global_nten)) then
           call get_kappa_dirichlet( &
             sdim, &
             alpha, &
             nmat, &
             im_inside, &
             kappa, &
             nf, &  ! points out of the cell
             xinside, &
             xinside_face)

           if (kappa.ge.0.0d0) then

            grad=rho_inside
            if (local_hflag.eq.1) then
             ! do nothing
            else if (local_hflag.eq.0) then

             if (1.eq.0) then
              xinside_face_debug(1)=0.47656250000000000d0
              xinside_face_debug(2)=0.35594045622499831d0
              mag_debug=0.0d0
              do dir2=1,2
               mag_debug=mag_debug+ &
                 (xinside_face_debug(dir2)-xinside_face(dir2))**2
              enddo
              mag_debug=sqrt(mag_debug)
              if (mag_debug.le.1.0D-12) then
               mag_debug=0.0d0
               do dir2=1,2
                mag_debug=mag_debug+ &
                 (xinside_face_debug(dir2)-0.5d0)**2
               enddo
               mag_debug=sqrt(mag_debug)
               print *,"mag_debug=",mag_debug
              endif
             endif

             call get_interface_temperature( &
               use_tsatfab, &
               i_tsat,j_tsat,k_tsat, &
               ireverse_tsat, &
               iten_tsat, &
               ntsat, &
               bfact_tsat, &
               level_tsat, &
               finest_level_tsat, &
               dx_tsat,xlo_tsat, &
               ngrow_tsat, &
               fablo_tsat,fabhi_tsat, &
               TSATFAB,DIMS(TSATFAB), &
               rho_I, &
               internal_multimaterial, &
               saturation_temp, &
               use_exact_temperature, &
               xinside_face, &
               time_div, &
               nmat, &
               global_nten,5)

             if (1.eq.0) then
              if (probtype.eq.403) then
               if (local_hflag.eq.0) then
                print *,"NORMAL (ext_face) DIRICHLET: im_inside,im_outside ", &
                  im_inside,im_outside
                print *,"i_tsat,j_tsat,k_tsat ", &
                   i_tsat,j_tsat,k_tsat
                print *,"dir_main,side ",dir_main,side
                print *,"rho_I ",rho_I
               else
                print *,"local_hflag invalid"
                stop
               endif 
              else
               print *,"probtype invalid not 403"
               stop
              endif
             endif

             grad=grad-rho_I

            else
             print *,"local_hflag invalid"
             stop
            endif

            ! grad approximates: k grad T dot (-nf)
            ! later: div_tot=div_tot+grad*AFRAC*dx(dir_main)
            grad=kappa*grad
           else
            print *,"kappa invalid"
            stop
           endif
          else
           print *,"internal_multimaterial invalid"
           stop
          endif
         else if (im_inside.eq.im_outside) then
          if (im_tilde.eq.im_inside) then  
           call get_kappa_simple( &
             sdim, &
             alpha, &
             nmat, &
             im_inside,im_outside, &
             kappa, &
             nf, &  ! points out of the cell
             xinside,xoutside, &
             xinside_face,xoutside_face)

            ! grad approximates: k grad T dot (-nf)
            ! later: div_tot=div_tot+grad*AFRAC*dx(dir_main)
           if (kappa.ge.0.0d0) then
            grad=kappa*(rho_inside-rho_outside)
           else
            print *,"kappa invalid"
            stop
           endif
          else if ((im_tilde.ge.1).and. &
                   (im_tilde.le.nmat).and. &
                   (im_tilde.ne.im_inside)) then
           kappa=alpha(im_inside)
           if (kappa.ge.0.0d0) then
            grad=-kappa*(gradT_face(1)*nf(1)+gradT_face(2)*nf(2))
           else
            print *,"kappa invalid"
            stop
           endif
          else
           print *,"im_tilde invalid"
           stop
          endif
         else
          print *,"im_inside or im_outside invalid"
          stop
         endif

         div_tot = div_tot+grad*AFRAC*dx(dir_main)

        elseif(AFRAC.gt.1.0d0)then
         print *,"frac_pair_cell bust in cell_div_cal_simple"
         stop
        elseif(AFRAC.lt.0.0d0) then
         print *,"frac_pair_cell invalid in cell_div_cal_simple"
         stop
        endif

       enddo ! side
      enddo ! dir_main
     enddo ! im_outside=1..nmat

     VP_dir_current=0
     VP_side_current=0

     do im_outside = 1,nmat

      VP_im_in_current=im_inside
      VP_im_out_current=im_outside

      AFRAC=int_face_inside(im_inside,im_outside)

      if (AFRAC.ge.FACETOL_DVOL*dx(1)) then
 
       if (im_inside.ne.im_outside) then
        ! do nothing
       else
        print *,"im_inside==im_outside cannot happen"
        stop
       endif

       ! nf points from im_inside to im_outside.
       do dir2=1,sdim
        nf(dir2)=int_face_normal(im_inside,im_outside,dir2)
       enddo

       rho_outside=local_rho(0,0,im_outside)

       len_inside=0.0d0
       len_outside=0.0d0
       do dir2=1,sdim
        xoutside(dir2)=mat_cen_sten(0,0,im_outside,dir2)
        xinside_face(dir2)=xclosest_sten(0,0,im_inside,im_outside,dir2)
        xoutside_face(dir2)=xclosest_sten(0,0,im_outside,im_inside,dir2)
        len_inside=len_inside+(xinside(dir2)-xinside_face(dir2))**2
        len_outside=len_outside+(xoutside(dir2)-xoutside_face(dir2))**2
       enddo ! dir2=1..sdim

       VP_areaface_current=AFRAC
       VP_test_number_current=7
       call VP_sanity_check(xinside_face,LS_sub,xsten)
       VP_test_number_current=27
       call VP_sanity_check(xoutside_face,LS_sub,xsten)

       len_inside=sqrt(len_inside)
       len_outside=sqrt(len_outside)
       len_total=len_inside+len_outside
       if ((len_inside.gt.0.0d0).and. &
           (len_outside.gt.0.0d0)) then
        do dir2=1,sdim
         xface(dir2)=(len_outside*xinside_face(dir2)+ &
                      len_inside*xoutside_face(dir2))/len_total
        enddo
        do dir2=1,sdim
         if (operator_internal.eq.1) then ! simple method
          xinside_face(dir2)=xface(dir2) 
          xoutside_face(dir2)=xface(dir2) 
         else if (operator_internal.eq.2) then ! Dai and Scannapieco
          xinside_face(dir2)=int_centroid_sten(0,0,im_inside,im_outside,dir2)
          xoutside_face(dir2)=int_centroid_sten(0,0,im_outside,im_inside,dir2)
         else if (operator_internal.eq.3) then ! orthogonal probe
          ! do nothing
         else
          print *,"operator_internal invalid"
          stop
         endif
        enddo !dir2=1..sdim

        if ((operator_internal.eq.1).or. &  ! simple
            (operator_internal.eq.2)) then  ! DS
         mag_debug=0.0d0
         do dir2=1,sdim
          mag_debug=mag_debug+(xinside_face(dir2)-xoutside_face(dir2))**2
         enddo
         mag_debug=sqrt(mag_debug)
         if (mag_debug.ge.VOFTOL*dx(1)) then
          print *,"xinside_face and xoutside_face should be the same"
          print *,"AFRAC=",AFRAC
          print *,"xinside_face ",xinside_face(1),xinside_face(2)
          print *,"xoutside_face ",xoutside_face(1),xoutside_face(2)
          print *,"VP_i_current ",VP_i_current
          print *,"VP_j_current ",VP_j_current
          print *,"VP_im_in_current ",VP_im_in_current
          print *,"VP_im_out_current ",VP_im_out_current
          stop
         endif
         VP_areaface_current=AFRAC
         VP_test_number_current=37
         call VP_inside_check(xinside_face,xsten)
         VP_test_number_current=38
         call VP_inside_check(xoutside_face,xsten)
        else if (operator_internal.eq.3) then ! OP
         ! check nothing
        else
         print *,"operator_internal invalid"
         stop
        endif

        VP_areaface_current=AFRAC
        VP_test_number_current=8
        call VP_sanity_check(xinside_face,LS_sub,xsten)

        VP_areaface_current=AFRAC
        VP_test_number_current=9
        call VP_sanity_check(xoutside_face,LS_sub,xsten)

        call init_multimaterial_flag(im_inside,im_outside,nmat, &
          internal_multimaterial,im_source,im_dest)
        if ((internal_multimaterial.le.-1).and. &
            (internal_multimaterial.ge.-global_nten)) then 
         call get_kappa_simple( &
           sdim, &
           alpha, &
           nmat, &
           im_inside,im_outside, &
           kappa, &
           nf, &  ! points out of the cell
           xinside,xoutside, &
           xinside_face,xoutside_face)

           ! grad approximates: k grad T dot (-nf)
           ! later: div_tot=div_tot+grad*AFRAC*dx(dir_main)
         grad=kappa*(rho_inside-rho_outside)
        else if ((internal_multimaterial.ge.1).and. &
                 (internal_multimaterial.le.2*global_nten)) then
         call get_kappa_dirichlet( &
           sdim, &
           alpha, &
           nmat, &
           im_inside, &
           kappa, &
           nf, &  ! points out of the cell
           xinside, &
           xinside_face)

         grad=rho_inside
         if (local_hflag.eq.1) then
          ! do nothing
         else if (local_hflag.eq.0) then
          call get_interface_temperature( &
            use_tsatfab, &
            i_tsat,j_tsat,k_tsat, &
            ireverse_tsat, &
            iten_tsat, &
            ntsat, &
            bfact_tsat, &
            level_tsat, &
            finest_level_tsat, &
            dx_tsat,xlo_tsat, &
            ngrow_tsat, &
            fablo_tsat,fabhi_tsat, &
            TSATFAB,DIMS(TSATFAB), &
            rho_I, &
            internal_multimaterial, &
            saturation_temp, &
            use_exact_temperature, &
            xinside_face, &
            time_div, &
            nmat, &
            global_nten,6)

          if (1.eq.0) then
           if (probtype.eq.403) then
            if (local_hflag.eq.0) then
             print *,"NORMAL (int_face) DIRICHLET: im_inside,im_outside ", &
                  im_inside,im_outside
             print *,"i_tsat,j_tsat,k_tsat ", &
                i_tsat,j_tsat,k_tsat
             print *,"rho_I ",rho_I
            else
             print *,"local_hflag invalid"
             stop
            endif 
           else
            print *,"probtype invalid not 403"
            stop
           endif
          endif

          grad=grad-rho_I
         else
          print *,"local_hflag invalid"
          stop
         endif
          ! grad approximates: k grad T dot (-nf)
          ! later: div_tot=div_tot+grad*AFRAC
         grad=kappa*grad
        else
         print *,"internal_multimaterial invalid"
         stop
        endif

        div_tot = div_tot + grad*AFRAC
       else
        print *,"len_inside or len_outside invalid"
        stop
       endif
      else if ((AFRAC.ge.0.0d0).and. &
               (AFRAC.le.FACETOL_DVOL*dx(1))) then
       ! do nothing
      else
       print *,"AFRAC invalid"
       stop
      endif
     enddo ! im_outside=1..nmat

    else
     print *,"im_in invalid" 
     stop
    endif

  end subroutine cell_div_cal_simple



  !--------------------------------------------------------------
  !////////////////////////////////////////////////////////////////////
  !////////////////////////////////////////////////////////////////////
  subroutine ext_grad(nmat,sdim,flag,kappa,mat_osd,mat_isd, T,&
         cell, dx,mat_cen_sten,mofdata, &
         grad,success_flag)
    implicit none

    integer,intent(in)            :: nmat,sdim 
    real(kind=8),intent(in)       :: mat_cen_sten(-1:1,-1:1,nmat,sdim)
    real(kind=8),intent(in)       :: kappa(nmat)

    type(polygon)                 :: cell(-1:1, -1:1)
    real(kind=8),intent(in)       :: dx(sdim)
    integer,intent(in)            :: flag
    integer,intent(in)            :: mat_osd,mat_isd
    ! material number
    real(kind=8),intent(in)       :: mofdata(-1:1,-1:1,(sdim*2+3)*nmat)
    real(kind=8),dimension(-1:1,-1:1,nmat),intent(in) :: T
    integer :: success_flag

    integer                       :: i,j
    integer                       :: datalen
    integer                       :: tmat

    integer,dimension(-1:1,-1:1)           :: ttag
    integer,dimension(-1:1,-1:1)           :: osd_tag,isd_tag
    real(kind=8),dimension(-1:1,-1:1,sdim) :: lsdata_osd, lsdata_isd,lsdata
    real(kind=8),dimension(-1:1,-1:1)      :: w_osd, w_isd,w
  
    real(kind=8)                  :: normal(sdim)
    real(kind=8)                  :: local_eps
    type(points)                  :: fluxpoint
    type(points)                  :: x1,x2
    
    real(kind=8),allocatable    :: A(:,:) , b(:), x_solved(:)
    real(kind=8), intent(out)   :: grad
    real(kind=8) ,external      :: norm_2d

    local_eps=1.0d-8

    success_flag=1 ! ext_grad

    tmat = 0
    datalen = (sdim * 2 + 3)
    osd_tag = 0
    isd_tag = 0
    ttag = 0
    lsdata_isd = 0.0d0
    lsdata_osd = 0.0d0
    lsdata = 0.0d0
    normal = 0.0d0
    w_osd = 0.0d0
    w_isd = 0.0d0
    w = 0.0d0

!     print *,"mofdata", mofdata(0,0,1:7)
!     print *,"mofdata",mofdata(0,0,8:14)

    if (1.eq.0) then
     print *,"mat_cen_sten"

     do i = -1,1
     do j = -1,1
      print *,"i=",i,"j=",j
      print *, mat_cen_sten(i,j,1,:)
     enddo
     enddo
    endif

     if (flag .eq. 1) then

          normal(1) = -1
          normal(2) = 0

     elseif (flag .eq. 2) then

          normal(1) = 1
          normal(2) = 0

     elseif(flag .eq. 3)then

          normal(1) = 0
          normal(2) = -1

     elseif(flag .eq. 4)then

          normal(1) = 0
          normal(2) = 1
     else
          print *,"invalid flag in  ext_grad"
     endif


    if (mat_osd .eq. mat_isd) then
       tmat = mat_osd
       if (flag .eq. 1) then
                       
          do j = -1, 1     
             
           if (mofdata(-1, j, datalen * (tmat - 1) + 1) .gt. local_eps) then
                ttag(-1, j) = 1
                lsdata(-1, j, :) = mat_cen_sten(-1, j, tmat, :)
           endif

           if (mofdata(0, j, datalen * (tmat - 1) + 1) .gt. local_eps) then
                ttag(0, j) = 1
                lsdata(0, j, :) = mat_cen_sten(0, j, tmat, :)
           endif

           if (mofdata(-1, j, datalen * (mat_osd - 1) + 1) .gt. local_eps) then
                osd_tag(-1, j) = 1
                lsdata_osd(-1, j, :) = mat_cen_sten(-1, j, mat_osd, :)
           endif

           if (mofdata(0, j, datalen * (mat_isd - 1) + 1) .gt. local_eps) then
                isd_tag(0, j) = 1
                lsdata_isd(0, j, :) = mat_cen_sten(0, j, mat_isd, :)
           endif
          enddo


       elseif (flag .eq. 2) then

          do j = -1, 1     

           if (mofdata(1, j, datalen * (tmat-1) + 1) .gt. local_eps) then
                ttag(1, j) = 1
                lsdata(1, j, :) = mat_cen_sten(1, j, tmat, :)
           endif

           if(mofdata(0, j, datalen * (tmat - 1) + 1) .gt. local_eps) then
                ttag(0, j) = 1
                lsdata(0, j, :) = mat_cen_sten(0, j, tmat, :)
           endif
           if (mofdata(1, j, datalen * (mat_osd-1) + 1) .gt. local_eps) then
              osd_tag(1, j) = 1
              lsdata_osd(1, j, :) = mat_cen_sten(1, j, mat_osd, :)
           endif

           if(mofdata(0, j, datalen * (mat_isd - 1) + 1) .gt. local_eps) then
             isd_tag(0, j) = 1
             lsdata_isd(0, j, :) = mat_cen_sten(0, j, mat_isd, :)
           endif
          enddo

       elseif(flag .eq. 3)then

          do i = -1, 1     
           if (mofdata(i, -1, datalen * (tmat-1) + 1) .gt. local_eps) then
              ttag(i, -1) = 1
              lsdata(i, -1, :) = mat_cen_sten(i, -1, tmat, :)
           endif

           if(mofdata(i, 0, datalen * (tmat - 1) + 1) .gt. local_eps) then
            ttag(i, 0) = 1
            lsdata(i, 0, :) = mat_cen_sten(i, 0, tmat, :)
           endif
           if (mofdata(i, -1, datalen * (mat_osd-1) + 1) .gt. local_eps) then
            osd_tag(i, -1) = 1
            lsdata_osd(i, -1, :) = mat_cen_sten(i, -1, mat_osd, :)
           endif

           if(mofdata(i, 0, datalen * (mat_isd - 1) + 1) .gt. local_eps) then
              isd_tag(i, 0) = 1
              lsdata_isd(i, 0, :) = mat_cen_sten(i, 0, mat_isd, :)
           endif
        enddo


       elseif(flag .eq. 4)then

          do i = -1, 1     
             if (mofdata(i, 1, datalen * (tmat-1) + 1) .gt. local_eps) then
                ttag(i, 1) = 1
                lsdata(i, 1, :) = mat_cen_sten(i, 1, tmat, :)
             endif

             if(mofdata(i, 0, datalen * (tmat - 1) + 1) .gt. local_eps) then
                ttag(i, 0) = 1
                lsdata(i, 0, :) = mat_cen_sten(i, 0, tmat, :)
             endif

             if (mofdata(i, 1, datalen * (mat_osd-1) + 1) .gt. local_eps) then
                osd_tag(i, 1) = 1
                lsdata_osd(i, 1, :) = mat_cen_sten(i, 1, mat_osd, :)
             endif

             if(mofdata(i, 0, datalen * (mat_isd - 1) + 1) .gt. local_eps) then
                isd_tag(i, 0) = 1
                lsdata_isd(i, 0, :) = mat_cen_sten(i, 0, mat_isd, :)
             endif

          enddo

       else
          print *,"invalid flag in  ext_grad"
       endif

  else

       do i = -1 , 1
          do j = -1, 1     
           
           if (mofdata(i, j, datalen * (mat_osd - 1) + 1) .gt. local_eps) then
             osd_tag(i, j) = 1
             lsdata_osd(i, j, :) = mat_cen_sten(i, j, mat_osd, :)
           endif

           if (mofdata(i, j, datalen * (mat_isd - 1) + 1) .gt. local_eps) then
              isd_tag(i, j) = 1
              lsdata_isd(i, j, :) = mat_cen_sten(i, j, mat_isd, :)
           endif

          enddo
        enddo

   endif

! find flux point,  intersection of interface and straight line through two material

  if(flag .eq. 1) then
          x1%val(1) = lsdata_osd(-1,0,1) 
          x1%val(2) = lsdata_osd(-1,0,2) 
          x2%val(1) = lsdata_isd(0,0,1) 
          x2%val(2) = lsdata_isd(0,0,2) 
          call find_fluxpoint(-1.0d0,0.0d0,cell(0,0)%center%val(1)-0.5d0*dx(1),x1,x2,fluxpoint)
   
    
  elseif(flag .eq. 2) then
          x1%val(1) = lsdata_osd(1,0,1) 
          x1%val(2) = lsdata_osd(1,0,2) 
          x2%val(1) = lsdata_isd(0,0,1) 
          x2%val(2) = lsdata_isd(0,0,2) 
          call find_fluxpoint(-1.0d0,0.0d0,cell(0,0)%center%val(1)+0.5d0*dx(1),x1,x2,fluxpoint)

  elseif(flag .eq. 3) then

          x1%val(1) = lsdata_osd(0,-1,1) 
          x1%val(2) = lsdata_osd(0,-1,2) 
          x2%val(1) = lsdata_isd(0,0,1) 
          x2%val(2) = lsdata_isd(0,0,2) 
         call find_fluxpoint(0.0d0,-1.0d0,cell(0,0)%center%val(2)-0.5d0*dx(2),x1,x2,fluxpoint)

  elseif(flag .eq. 4) then
    
          x1%val(1) = lsdata_osd(0,1,1) 
          x1%val(2) = lsdata_osd(0,1,2) 
          x2%val(1) = lsdata_isd(0,0,1) 
          x2%val(2) = lsdata_isd(0,0,2) 
          call find_fluxpoint(0.0d0,-1.0d0,cell(0,0)%center%val(2)+0.5d0*dx(2),x1,x2,fluxpoint)
  else
    print *,"invalid flag in  ext_grad"
    stop
  endif
  
! calculate the weight
if(mat_osd .eq. mat_isd) then


  do i = -1 , 1
   do j = -1, 1
     if(ttag(i,j) .ne. 0)then
       CALL SDD(dx(1)/2, (NORM_2d(lsdata(i,j,1),lsdata(i,j,2), &
               & fluxpoint%val(1),fluxpoint%val(2))**2.0d0)/dx(1), w(i,j))
       if(w(i,j) .lt. local_eps) then
        print *,"weight invalid"
        print *,i,j
        print *,w(i,j)
        print *,lsdata(i,j,:)
        print *,fluxpoint%val
        stop
       endif
     endif
   enddo
  enddo
  
  allocate(A(3,3),b(3),x_solved(3))
  call Matrix_assemble_single(sdim, nmat, dx, tmat, ttag, &
                        & lsdata, w, T, &
                        & normal, fluxpoint, kappa(tmat),&
                        & A,b)
  call solve_determinant(3,A,b,x_solved,success_flag)
  grad = x_solved(2)
  deallocate(A,b,x_solved)
  

else

  do i = -1 , 1
   do j = -1 , 1
    if(osd_tag(i,j) .ne. 0) then
      CALL SDD(dx(1)/2, (NORM_2d(lsdata_osd(i,j,1),lsdata_osd(i,j,2), &
               & fluxpoint%val(1),fluxpoint%val(2)))**2.0d0/dx(1), w_osd(i,j))
    endif

    if(isd_tag(i,j) .ne. 0) then
      CALL SDD(dx(1)/2, (NORM_2d(lsdata_isd(i,j,1),lsdata_isd(i,j,2), &
               & fluxpoint%val(1),fluxpoint%val(2)))**2.0d0/dx(1), w_isd(i,j))
    endif   

   enddo
  enddo  
 !do i = -1,1
 !    print *, "osd",w_osd(i,:) 
 !enddo
 !do i = -1,1
 !    print *, "isd",w_isd(i,:) 
 !enddo

! assemble the matrix
  allocate(A(4,4),b(4),x_solved(4))
  call Matrix_assemble_double(sdim, nmat, dx, mat_osd,mat_isd,osd_tag, &
    isd_tag, &
    lsdata_osd, lsdata_isd, w_osd, w_isd, T, &
    normal, fluxpoint, kappa(mat_osd), kappa(mat_isd), &
    A,b)

! use cramer's rule to solve for the flux  
  call solve_determinant(4,A,b,x_solved,success_flag)
  grad = x_solved(3)
   deallocate(A,b,x_solved)
 endif

  end subroutine ext_grad

!-----------------------------------------------------------------------------
subroutine int_grad(nmat,sdim,kappa,mat_osd,mat_isd, T,&
                    & cell, dx,mat_cen_sten,mofdata, &
                    & grad,success_flag)
    implicit none

    integer                       :: success_flag
    integer,intent(in)            :: nmat,sdim
    real(kind=8),intent(in)       :: mat_cen_sten(-1:1,-1:1,nmat,sdim)
    real(kind=8),intent(in)       :: kappa(nmat)
    type(polygon)                 :: cell(-1:1,-1:1)
    real(kind=8)                  :: center(2)
    real(kind=8),intent(in)       :: dx(sdim)
    integer,intent(in)            :: mat_osd,mat_isd


    real(kind=8),intent(in)       :: mofdata(-1:1,-1:1,(sdim*2+3)*nmat)
    real(kind=8),dimension(-1:1,-1:1,nmat),intent(in) :: T

    integer                       :: i,j
    integer                       :: datalen

    integer,dimension(-1:1,-1:1)           :: osd_tag,isd_tag
    real(kind=8),dimension(-1:1,-1:1,sdim) :: lsdata_osd, lsdata_isd
    real(kind=8),dimension(-1:1,-1:1)      :: w_osd, w_isd

    real(kind=8)                  :: normal(sdim),centroid(sdim), centroid1(sdim)
    real(kind=8)                  :: c
    real(kind=8)                  :: local_eps
    type(points)                  :: fluxpoint
    type(points)                  :: x1,x2

    real(kind=8),external         :: norm_2d

    real(kind=8)                :: A(4,4) , b(4), x_solved(4)
    real(kind=8), intent(out)   :: grad

    success_flag=1 ! int_grad

    local_eps=1.0d-8

    datalen = (sdim * 2 + 3)
    osd_tag = 0
    isd_tag = 0
    lsdata_isd = 0.0d0
    lsdata_osd = 0.0d0
    normal = 0.0d0
    w_osd = 0.0d0
    w_isd = 0.0d0
    c = 0.0d0
    centroid = 0.0d0
    centroid1 = 0.0d0

    center(1) = cell(0,0)%center%val(1)
    center(2) = cell(0,0)%center%val(2)

  
    do i  = 1, sdim
      normal(i) = mofdata(0,0,datalen*(mat_isd - 1)+sdim+2+i)
      centroid(i) = mofdata(0,0,datalen*(mat_isd - 1) + 1+i)
      centroid(i) = centroid(i) + center(i)
      centroid1(i) = mofdata(0,0,datalen*(mat_osd - 1) + 1+i)
      centroid1(i) = centroid1(i) + center(i)
    enddo


    

    if(abs(normal(1)) .lt. local_eps) then
       normal(1)= 0.0d0
    endif
    if(abs(normal(2)) .lt. local_eps) then
       normal(2)= 0.0d0
    endif


    do  i = 1,sdim
      c =  c - normal(i)*center(i)  
    enddo
      c = c + mofdata(0,0,datalen*mat_isd)
!...........................................................................................
    do i = 1, sdim
     normal(i) = -1.0d0*normal(i)
    enddo                                                       ! ATTENSION
    C = -1.0D0*C                                           
!.......................................................................................
    x1%val = centroid
    x2%val = centroid1

    call find_fluxpoint(normal(1),normal(2),c,x1,x2,fluxpoint)

    do i  = -1 , 1
      do j = -1 , 1
         if(mofdata(i,j, datalen * (mat_osd-1) + 1) .gt. local_eps) then
            osd_tag(i, j) = 1
            lsdata_osd(i, j, :) = mat_cen_sten(i, j, mat_osd, :)                        
         endif
         if(mofdata(i,j, datalen * (mat_isd-1) + 1) .gt. local_eps) then
            isd_tag(i, j) = 1
            lsdata_isd(i, j, :) = mat_cen_sten(i, j, mat_isd, :)                        
         endif
      enddo
    enddo


! calculate the weight
  do i = -1 , 1
   do j = -1 , 1
    if(osd_tag(i,j) .ne. 0) then
      CALL SDD(dx(1)/2, (NORM_2d(lsdata_osd(i,j,1),lsdata_osd(i,j,2), &
               & fluxpoint%val(1),fluxpoint%val(2)))**2.0d0/dx(1), w_osd(i,j))
    endif

    if(isd_tag(i,j) .ne. 0) then
        CALL SDD(dx(1)/2, (NORM_2d(lsdata_isd(i,j,1),lsdata_isd(i,j,2), &
               & fluxpoint%val(1),fluxpoint%val(2)))**2.0d0/dx(1), w_isd(i,j))
    endif   

   enddo
  enddo



! assemble the matrix
  call Matrix_assemble_double(sdim, nmat, dx,mat_osd,mat_isd, osd_tag, isd_tag, &
                        & lsdata_osd, lsdata_isd, w_osd, w_isd, T, &
                        & normal, fluxpoint, kappa(mat_osd), kappa(mat_isd), &
                        & A,b)



! use cramer's rule to solve for the flux  
  call solve_determinant(4,A,b,x_solved,success_flag)
  grad = x_solved(3)



end subroutine int_grad
!---------------------------------------------------------------------------
subroutine Matrix_assemble_single(sdim, nmat, dx, imat,tag, &
                            & lsdata, w, T, &
                            & normal, fluxpoint, k,&
                            & A,b)
implicit none

    integer,intent(in)            :: nmat,sdim 
    real(kind=8),intent(in)       :: dx(sdim)
    integer,intent(in)                                :: imat
    integer,dimension(-1:1,-1:1),intent(in)           :: tag
    real(kind=8),dimension(-1:1,-1:1,sdim),intent(in) :: lsdata
    real(kind=8),dimension(-1:1,-1:1),intent(in)      :: w
    real(kind=8),dimension(-1:1,-1:1,nmat),intent(in) :: T
   
    real(kind=8),intent(in)                  :: normal(sdim)
    real(kind=8)                  :: m(sdim)
    type(points)                  :: fluxpoint
    real(kind=8),intent(in)       :: k

    

    integer                       :: i,j
    real(kind=8)                  :: A(3,3)
    real(kind=8)                  :: b(3)
    real(kind=8)                  :: x,y,x1,y1,local_eps

    local_eps=1.0d-8

    m = normal
    A = 0.0d0
    b = 0.0d0

    do i = -1 , 1
     do j = -1 , 1
       if(tag(i,j) .eq. 1) then
         
         x = lsdata(i,j,1) - fluxpoint%val(1)
         y = lsdata(i,j,2) - fluxpoint%val(2) 
         A(1,1) = A(1,1) + 2.0d0*w(i,j)*((m(1)*y - m(2)*x)**2.0d0)
         A(1,2) = A(1,2) + 2.0d0*w(i,j)*((m(1)*x + m(2)*y)/k)*(m(1)*y - m(2)*x)
         A(1,3) = A(1,3) + 2.0d0*w(i,j)*(m(1)*y - m(2)*x)
         b(1) = b(1) +  2.0d0*w(i,j)*(m(1)*y - m(2)*x)*T(i,j,imat)       ! ----------T(i,j)-------
       endif
      enddo
     enddo  


    A(2,1) = A(1,2)
    do i = -1 , 1
      do j = -1 , 1
        if(tag(i,j) .eq. 1) then
         x1 = lsdata(i,j,1) - fluxpoint%val(1)
         y1 = lsdata(i,j,2) - fluxpoint%val(2) 
         A(2,2) = A(2,2) + 2.0d0*w(i,j)*((m(1)*x1 + m(2)*y1)**2.0d0)/(k*k)
         A(2,3) = A(2,3) + 2.0d0*w(i,j)*(m(1)*x1 + m(2)*y1)/k
         b(2) = b(2) + 2.0d0*w(i,j)*(m(1)*x1 + m(2)*y1)/k*T(i,j,imat)
        endif       
      enddo
    enddo


    A(3,1) = A(1,3) 
    A(3,2) = A(2,3)
    do i = -1 ,1
     do j = -1 ,1
        if(tag(i,j) .eq. 1) then          
         A(3,3) = A(3,3) + 2.0d0*w(i,j)
         b(3) = b(3) + 2.0d0*w(i,j)*T(i,j,imat)
        endif             
     enddo
    enddo

  do i = 1 , 3
   do j = 1 , 3
     if(abs(A(i,j)) .lt. local_eps)then
       A(i,j) = 0.0d0
     endif
   enddo
   if(abs(b(i)) .lt. local_eps) then
     b(i) = 0.0d0
   endif
  enddo
 

end subroutine Matrix_assemble_single


!------------------------------------------------------------------------------
 
 subroutine Matrix_assemble_double(sdim, nmat, dx, mat_osd,mat_isd,osd_tag, isd_tag, &
                            & lsdata_osd, lsdata_isd, w_osd, w_isd, T, &
                            & normal, fluxpoint, k_osd, k_isd,&
                            & A,b)
 implicit none
 
    integer,intent(in)            :: nmat,sdim 
    real(kind=8),intent(in)       :: dx(sdim)
    integer,intent(in)                                :: mat_osd,mat_isd
    integer,dimension(-1:1,-1:1),intent(in)           :: osd_tag,isd_tag
    real(kind=8),dimension(-1:1,-1:1,sdim),intent(in) :: lsdata_osd, lsdata_isd
    real(kind=8),dimension(-1:1,-1:1),intent(in)      :: w_osd, w_isd
    real(kind=8),dimension(-1:1,-1:1,nmat),intent(in) :: T
   
    real(kind=8),intent(in)                  :: normal(sdim)
    real(kind=8)                  :: m(sdim)
    type(points)                  :: fluxpoint
    real(kind=8),intent(in)       :: k_osd, k_isd
    

    integer                       :: i,j
    real(kind=8)                  :: A(4,4)
    real(kind=8)                  :: b(4)
    real(kind=8)                  :: x,y,x1,x2,y1,y2
    real(kind=8)                  :: local_eps

    local_eps=1.0d-8

    m = normal
    A = 0.0d0
    b = 0.0d0


    do i = -1 , 1
     do j = -1 , 1
       if(osd_tag(i,j) .eq. 1) then
         x = lsdata_osd(i,j,1) - fluxpoint%val(1)
         y = lsdata_osd(i,j,2) - fluxpoint%val(2) 
         A(1,1) = A(1,1) + 2.0d0*w_osd(i,j)*((m(1)*y - m(2)*x)**2.0d0)
         A(1,3) = A(1,3) + 2.0d0*w_osd(i,j)*((m(1)*x + m(2)*y)/k_osd)*(m(1)*y - m(2)*x)
         A(1,4) = A(1,4) + 2.0d0*w_osd(i,j)*(m(1)*y - m(2)*x)
         b(1) = b(1) +  2.0d0*w_osd(i,j)*(m(1)*y - m(2)*x)*T(i,j,mat_osd)       ! ----------T(i,j)-------
       endif
      enddo
     enddo
  
    do i = -1 , 1
     do j = -1 , 1       
       if(isd_tag(i,j) .eq. 1) then
         x = lsdata_isd(i,j,1) - fluxpoint%val(1)
         y = lsdata_isd(i,j,2) - fluxpoint%val(2)
         A(2,2) = A(2,2) + 2.0d0*w_isd(i,j)*((m(1)*y - m(2)*x)**2.0d0)
         A(2,3) = A(2,3) + 2.0d0*w_isd(i,j)*((m(1)*x + m(2)*y)/k_isd)*(m(1)*y - m(2)*x)
         A(2,4) = A(2,4) + 2.0d0*w_isd(i,j)*(m(1)*y - m(2)*x)
         b(2) = b(2) +  2.0d0*w_isd(i,j)*(m(1)*y - m(2)*x)*T(i,j,mat_isd)          
       endif
     enddo
    enddo


    A(3,1) = A(1,3)
    A(3,2) = A(2,3)
    do i = -1 , 1
      do j = -1 , 1
        if(osd_tag(i,j) .eq. 1) then
         x1 = lsdata_osd(i,j,1) - fluxpoint%val(1)
         y1 = lsdata_osd(i,j,2) - fluxpoint%val(2)           
         A(3,3) = A(3,3) + 2.0d0*w_osd(i,j)*((m(1)*x1 + m(2)*y1)**2.0d0)/((k_osd)**2.0d0)
         A(3,4) = A(3,4) + 2.0d0*w_osd(i,j)*(m(1)*x1 + m(2)*y1)/k_osd
         b(3) = b(3) + 2.0d0*w_osd(i,j)*(m(1)*x1 + m(2)*y1)/k_osd*T(i,j,mat_osd)
        endif
 
        if(isd_tag(i,j) .eq. 1) then
         x2 = lsdata_isd(i,j,1) - fluxpoint%val(1)
         y2 = lsdata_isd(i,j,2) - fluxpoint%val(2)           
         A(3,3) = A(3,3) + 2.0d0*w_isd(i,j)*((m(1)*x2 + m(2)*y2)**2.0d0)/((k_isd)**2.0d0)
         A(3,4) = A(3,4) + 2.0d0*w_isd(i,j)*(m(1)*x2 + m(2)*y2)/k_isd
         b(3) = b(3) + 2.0d0*w_isd(i,j)*(m(1)*x2 + m(2)*y2)/k_isd*T(i,j,mat_isd)
        endif        
      enddo
    enddo


    A(4,1) = A(1,4) 
    A(4,2) = A(2,4)
    A(4,3) = A(3,4)
    do i = -1 ,1
     do j = -1 ,1
        if(osd_tag(i,j) .eq. 1) then          
         A(4,4) = A(4,4) + 2.0d0*w_osd(i,j)
         b(4) = b(4) + 2.0d0*w_osd(i,j)*T(i,j,mat_osd)
        endif

        if(isd_tag(i,j) .eq. 1) then           
         A(4,4) = A(4,4) + 2.0d0*w_isd(i,j)
         b(4) = b(4) + 2.0d0*w_isd(i,j)*T(i,j,mat_isd)
        endif              
     enddo
    enddo

  do i = 1 , 4
   do j = 1 , 4
     if(abs(A(i,j)) .lt. local_eps)then
       A(i,j) = 0.0d0
     endif
   enddo
   if(abs(b(i)) .lt. local_eps) then
     b(i) = 0.0d0
   endif
  enddo




  end subroutine Matrix_assemble_double



  subroutine int_face_adjust(nmat,int_facefrac, int_face )
   implicit none
   
   integer,intent(in)            :: nmat 
   real(kind=8),intent(in)       :: int_facefrac(nmat+1,nmat)
   real(kind=8)                  :: int_face(nmat,nmat)
   real(kind=8)                  :: total_area,areafrac

   integer                       :: i , j

   do i = 1, nmat
    do j = 1, nmat
     int_face(i,j)=0.0
    enddo
   enddo

   do i = 1, nmat

     total_area=int_facefrac(nmat+1,i)

     do j = 1, nmat

      if (i.ne.j) then

       areafrac=int_facefrac(i,j)
       if (areafrac.lt.0.0) then
        print *,"areafrac invalid"
        stop
       else if (areafrac.eq.0.0) then
        ! do nothing
       else if ((areafrac.gt.0.0).and.(areafrac.le.1.0)) then
        if (total_area.le.0.0) then
         print *,"total_area invalid"
         stop
        endif
        int_face(i,j)=areafrac*total_area
        int_face(j,i)=int_face(i,j)
       else
        print *,"areafrac invalid"
        stop
       endif

      endif ! i<>j

     enddo ! j=1..nmat

   enddo ! i=1..nmat


  end subroutine int_face_adjust





  !------------------------------------------------------------
  subroutine convert_cen(nmat,sdim,N,centroid,centroid_mult)
    implicit none

    integer,intent(in)       :: nmat, sdim, N
    type(points),intent(in),dimension(-1:N,-1:N,nmat)  :: centroid
    real(kind=8),dimension(-1:N,-1:N,nmat,sdim)        :: centroid_mult

    integer                   :: i,j,im,dir

    do i = -1,N
       do j= -1,N
          do im = 1,nmat
           do dir=1,sdim
             centroid_mult(i,j,im,dir) = centroid(i,j,im)%val(dir)
           enddo
          enddo
       enddo
    enddo


  end subroutine convert_cen

  !-----------------------------------------------------------
  SUBROUTINE INIT_KAPPA(NMAT,thermal_cond,kappa)
    IMPLICIT NONE

    INTEGER,INTENT(IN)          :: NMAT
    real(kind=8),intent(in)     :: thermal_cond(nmat)
    REAL(KIND=8)                :: KAPPA(nmat,nmat)

    integer                     :: i,j

    !write(16,*) thermal_cond

    do i = 1,nmat
       do j = 1, nmat
          kappa(i,j) = 2.0d0/( 1.0d0/thermal_cond(i) + 1.0d0/thermal_cond(j) )
       enddo
    enddo



  END SUBROUTINE INIT_KAPPA







  !///////////////////////////////////////////////////////////////////
  !-------------------------------------------------------------------
  !              flux calculation
  !-------------------------------------------------------------------
  !///////////////////////////////////////////////////////////////////


  !-----------------------------------------------------------------

  subroutine find_fluxpoint(a,b,c,x1,x2,thetaI)
    ! ax + by + c = 0

    implicit none

    real(kind=8) ,intent(in)    :: a,b,c
    type(points) ,intent(in)    :: x1,x2

    type(points) ,intent(out)   :: thetaI
    real(kind=8)                :: local_eps

    local_eps=1.0d-8

    if(abs(a) .lt. local_eps  .and. abs(b) .lt. local_eps)then
       print *,"wrong slope in find_fluxpoint"
    elseif(abs(a) .lt. local_eps) then
       if(abs(x1%val(2)-x2%val(2)) .lt. local_eps) then
          print *,"Error, the two lines are parallel"
       else
          thetaI%val(2) = c/(-b)
          thetaI%val(1) = ((x1%val(1)-x2%val(1))*c/(-b) + x1%val(1)*(x1%val(2)-x2%val(2)) &
               & - (x1%val(1)-x2%val(1))*x1%val(2))/(x1%val(2)-x2%val(2))
       endif
    elseif(abs(b) .lt. local_eps) then
       if(abs(x1%val(1)-x2%val(1)) .lt. local_eps) then
          print *,"Error, the two lines are parallel"
       else
          thetaI%val(1) = c/(-a)
          thetaI%val(2) = ((x1%val(2)-x2%val(2))*c/(-a) + (x1%val(1)-x2%val(1))*x1%val(2) &
               & -x1%val(1)*(x1%val(2)-x2%val(2)))/(x1%val(1)-x2%val(1))
       endif
    else
       thetaI%val(1) = (-(x1%val(1)-x2%val(1))*x1%val(2) + x1%val(1)*(x1%val(2)-x2%val(2)) &
            & - c/b*(x1%val(1)-x2%val(1)))   &
            &/(x1%val(2)-x2%val(2)+(x1%val(1)-x2%val(1))*a/b)
       thetaI%val(2) = -a/b*thetaI%val(1) - c/b
    endif




  end subroutine find_fluxpoint





  !////////////////////////////////////////////////////////////
  ! use cramer's rule to solve linear system
  !////////////////////////////////////////////////////////////
  !-------------------------------------------------------
  subroutine solve_determinant(N,A_in,b_in,x,success_flag)
    implicit none

    integer, intent(in)          :: N
    real(kind=8),intent(in)      :: A_in(N,N)
    real(kind=8),intent(in)      :: b_in(N)
    real(kind=8),intent(out)     :: x(N)
    integer :: success_flag
    integer                      :: i
    real(kind=8)                :: A(N,N)
    real(kind=8)                :: b(N)


    integer                      :: flag

    real(kind=8)                 :: deter(N)
    real(kind=8)                 :: den

    flag = 0

    A = A_in
    b = b_in
    
    do i = 1,N
      if(b(i) .ne. 0.0d0 ) then
       flag = 1
      endif
    enddo

    if (1.eq.0) then
     print *, "A = "
     do i = 1, N
      print *, A(:,i)
     enddo
    endif

     !  print *,"b"
     !  print *,b
    call deter_cal(N,A,den)
    if(abs(den) .le. 0.0) then
      if (1.eq.0) then
       print *,"warning, denominator determinant is zero"
       print *, den
      endif
      success_flag=0
    endif

    if(flag .eq. 0) then 
     if (1.eq.0) then
      print *,"Only zero solution exists"
     endif
     x = 0.0d0
    elseif(flag .eq. 1)then
     do i = 1,N
      A(:,i) = b(:)
      call deter_cal(N,A,deter(i))
      x(i) = deter(i)/den
      A = A_in  
     enddo
    else
     print *,"flag error in solve_determinant"
    endif

  end subroutine solve_determinant
  !--------------------------------------------------
  subroutine deter_cal(N,A_in,deter)
    implicit none

    integer, intent(in)       :: N
    real(kind=8),intent(in)   :: A_in(N,N)

    real(kind=8)              :: A(N,N)
    integer                   :: i,j,k
    integer                   :: scl
    !real(kind=8)              :: 
    integer                   :: sgn
    real(kind=8)                 :: coef

    real(kind=8),intent(out)   :: deter
    integer                    :: mark
    real(kind=8) :: local_eps

    local_eps=1.0d-8

    A = A_in
    sgn = 1
    mark = 1
  !     print *, "A_in"
  !     do i = 1, N
  !      print *, A_in(:,i) 
  !     enddo

    do j = 1,N-1
       scl = j
       if(abs(A(scl,scl)) .lt. local_eps) then
          sgn = sgn*(-1)
          call zero_pivot(N,A,scl,mark)
    !    print *, "mark",mark
       endif
      if(mark .ne. 0) then
       do i = scl+1,N
          coef = A(i,j)/A(scl,j)
          do k = scl,N
             A(i,k) = A(i,k) - A(scl,k)*coef
          enddo
       enddo
      else
        exit
      endif
    enddo

    if(mark .ne. 0) then
     deter = sgn
     do i = 1,N
       deter = deter*A(i,i)
     enddo
    else
     deter = 0.0d0
    endif

  end subroutine deter_cal

  !--------------------------------------------------
  subroutine zero_pivot(N,A,iin,mark)
    implicit none

    integer,intent(in)           :: N
    !real(kind=8),intent(in)      :: A_in(N,N)
    real(kind=8)                 :: A(N,N)
    real(kind=8)                 :: B(N),C(N)
    integer,intent(in)           :: iin
    integer                      :: i,mark
    real(kind=8)                 :: local_eps

    local_eps=1.0d-8
    mark = 0

    do i = iin+1,N       
       if(abs(A(i,iin)) .gt. local_eps) then
          mark = i
          exit
       endif
    enddo
    if(mark .ne. 0) then
!       print *,"warning, the determinant is 0"
!       print *, "where:"
!       do i = 1, N
!        print *, A(:,i) 
!       enddo
!      stop
    B(:) = A(iin,:)
    C(:) = A(mark,:)
    A(iin,:) = C(:)
    A(mark,:) = B(:)
    endif


  end subroutine zero_pivot




end module MOF_pair_module






