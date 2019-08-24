#undef BL_LANG_CC
#define BL_LANG_FORT

#include "BC_TYPES.H"

      module bicgstab_module
      USE probmain_module

      REAL*8 :: VOFTOL_local

       ! set ERRTOL in main.F90
       ! 0.01D0, 0.5D0, 0.99D0, 0.999999D0 are some options
      REAL*8 :: ERRTOL

      integer state_ncomp
      !0=low,1=simple,2=Dai and Scannapieco,3=orthogonal probe
      integer operator_internal
      integer operator_external
      integer linear_exact
      integer probtypeCG ! 0=flat interface  1=annulus 2=vertical interface
      integer sdim
      integer ngeom_reconCG
      integer nx,ny,nz,lox,loy,loz,hix,hiy,hiz
      integer nmat
      integer precond_type,nsmooth
      REAL*8 deltat
      REAL*8 bicgstab_tol
      REAL*8 h,meshvol
      REAL*8 alpha(100)
      REAL*8, dimension(:,:,:), allocatable :: beta
      REAL*8, dimension(:,:,:), allocatable :: UNEW
      REAL*8, dimension(:,:,:), allocatable :: UOLD
      REAL*8, dimension(:,:,:), allocatable :: VFRAC_MOF
      REAL*8, dimension(:,:,:), allocatable :: G
      REAL*8, dimension(:,:,:), allocatable :: DIAG_FIELD

      REAL*8, dimension(:,:,:), allocatable :: mofdata_FAB
      REAL*8                                :: current_time 
      REAL*8                                :: tnp1
      REAL*8                                :: tn_mid
      REAL*8                                :: time_source_global

      REAL*8, dimension(:,:,:,:),   allocatable :: xsten_FAB
      REAL*8, dimension(:,:,:,:),   allocatable :: int_face_FAB
       ! absolute coord
      REAL*8, dimension(:,:,:,:,:), allocatable :: int_centroid_FAB
      REAL*8, dimension(:,:,:,:,:), allocatable :: int_face_normal_FAB
      REAL*8, dimension(:,:,:,:),   allocatable :: dist_to_int_FAB
       ! absolute coord
      REAL*8, dimension(:,:,:,:,:), allocatable :: xclosest_FAB
       ! nmat+1,sdim,2
      REAL*8, dimension(:,:,:,:,:), allocatable :: ext_face_FAB
      REAL*8, dimension(:,:,:,:,:,:), allocatable :: multi_cen_cell_FAB
       ! i,j,im_outside,im_inside,dir,side
      REAL*8, dimension(:,:,:,:,:,:), allocatable :: frac_pair_cell_FAB
       ! i,j,im_outside,im_inside,xdir,dir,side
      REAL*8, dimension(:,:,:,:,:,:,:), allocatable :: x_pair_cell_FAB
      REAL*8, dimension(:,:,:,:), allocatable :: centroid_mult_FAB
      REAL*8 dx(3)
 
      contains


      subroutine WALLBC_DIAG(xpoint,ypoint,coeff2,bctype,dir,sidesten)
      IMPLICIT NONE

      integer dir,sidesten
      REAL*8 xpoint,ypoint
      integer bctype
      REAL*8 coeff2

      if (bctype.eq.REFLECT_EVEN) then
       coeff2=1.0
      else if (bctype.eq.EXT_DIR) then
       coeff2=-1.0
      else
       print *,"bctype invalid bctype=",bctype
       stop
      endif
      
      return
      end subroutine WALLBC_DIAG


   
! bctype=REFLECT_EVEN Neumann 
! bctype=EXT_DIR Dirichlet

      subroutine DEALLOCATE_GLOBALS()
      IMPLICIT NONE

      DEALLOCATE(beta)
      DEALLOCATE(UNEW)
      DEALLOCATE(UOLD)
      DEALLOCATE(VFRAC_MOF)
      DEALLOCATE(G)
      DEALLOCATE(DIAG_FIELD)

      DEALLOCATE(mofdata_FAB)

      DEALLOCATE(xsten_FAB)
      DEALLOCATE(int_face_FAB)
      DEALLOCATE(int_centroid_FAB)
      DEALLOCATE(int_face_normal_FAB)
      DEALLOCATE(dist_to_int_FAB)
      DEALLOCATE(xclosest_FAB)
      DEALLOCATE(ext_face_FAB)
      DEALLOCATE(multi_cen_cell_FAB)
      DEALLOCATE(frac_pair_cell_FAB)
      DEALLOCATE(x_pair_cell_FAB)
      DEALLOCATE(centroid_mult_FAB)
      
      return
      end subroutine DEALLOCATE_GLOBALS

      subroutine build_MM_diag()
      use MOF_pair_module
      use mmat_FVM, only: dist_concentric

      IMPLICIT NONE

      REAL*8 mat_cen_sten(-1:1,-1:1,nmat,sdim)

      REAL*8 frac_pair_cell(nmat,nmat,sdim,2)
      REAL*8 x_pair_cell(nmat,nmat,sdim,sdim,2)

      REAL*8 int_face_cell(nmat,nmat)
      REAL*8 int_face_sten(-1:1,-1:1,nmat,nmat)
      REAL*8 int_centroid_sten(-1:1,-1:1,nmat,nmat,sdim)
      REAL*8 int_face_normal_cell(nmat,nmat,sdim)
      REAL*8 int_facefrac_cell(nmat+1,nmat)
      REAL*8 int_centroid_cell(nmat,nmat,sdim)
      REAL*8 ext_facefrac_cell(nmat+1,sdim,2)
      REAL*8 dist_to_int_cell(nmat,nmat)
      REAL*8 dist_to_int_sten(-1:1,-1:1,nmat,nmat)
        ! absolute coord
      REAL*8 xclosest(nmat,nmat,sdim)
      REAL*8 xclosest_sten(-1:1,-1:1,nmat,nmat,sdim)
      REAL*8 xsten_cell(-3:3,sdim)
      REAL*8 multi_cen_cell(sdim,nmat,sdim,2)
      integer i,j,imof
      integer im1,im2,im
      integer im_outside,im_inside
      integer dir,sidesten,isten
      integer dircen,dir2
      integer ii,jj,vofcomp
      REAL*8 vf,diag_local
      REAL*8 mofdata_cell(ngeom_reconCG*nmat)
      REAL*8 mofdata_sten(-1:1,-1:1,ngeom_reconCG*nmat)
      real(kind=8) :: rho_box(-1:1,-1:1,nmat)

      REAL*8 thin_cen_sten(-1:1,sdim,sdim,nmat,sdim,2)

      REAL*8 ext_face_sten(-1:1,sdim,nmat+1,sdim,2)
      integer im_in
      integer diag_coeff_flag
      integer local_hflag
      integer local_linear_exact
      REAL*8 time_placeholder

      diag_coeff_flag=1
      local_hflag=1
      local_linear_exact=0
      time_placeholder=0.0d0

      VP_i_debug=-9999
      VP_j_debug=-9999

      do i=lox-1,hix+1
      do j=loy-1,hiy+1
       do im=1,nmat+1
       do dir=1,sdim
       do sidesten=1,2
        ext_face_FAB(i,j,im,dir,sidesten)=0.0
       enddo
       enddo
       enddo
       do im1=1,nmat
       do im2=1,nmat
        int_face_FAB(i,j,im1,im2)=0.0
        do dir=1,sdim
         int_centroid_FAB(i,j,im1,im2,dir)=0.0
        enddo
       enddo
       enddo
       do im1=1,nmat
       do im2=1,nmat
        dist_to_int_FAB(i,j,im1,im2)=0.0
        do dir=1,sdim
         xclosest_FAB(i,j,im1,im2,dir)=0.0
        enddo
       enddo
       enddo
       do im1=1,nmat
       do im2=1,nmat
       do dir=1,sdim
        int_face_normal_FAB(i,j,im1,im2,dir)=0.0
       enddo
       enddo
       enddo
       do dircen=1,sdim
       do im=1,nmat
       do dir2=1,sdim
       do sidesten=1,2
        multi_cen_cell_FAB(i,j,dircen,im,dir2,sidesten)=0.0
       enddo
       enddo
       enddo
       enddo
       do im_outside=1,nmat
       do im_inside=1,nmat
       do dir=1,sdim
       do sidesten=1,2
        frac_pair_cell_FAB(i,j,im_outside,im_inside,dir,sidesten)=0.0
        do dir2=1,sdim
         x_pair_cell_FAB(i,j,im_outside,im_inside,dir2,dir,sidesten)=0.0
        enddo
       enddo
       enddo
       enddo
       enddo
      enddo
      enddo ! i,j

      do i=lox-1,hix+1
      do j=loy-1,hiy+1
     
       VP_i_current=i 
       VP_j_current=j 

       do imof=1,ngeom_reconCG*nmat
        mofdata_cell(imof)=mofdata_FAB(i,j,imof)
       enddo
       do isten=-3,3
       do dir=1,sdim
        xsten_cell(isten,dir)=xsten_FAB(i,j,isten,dir)
       enddo
       enddo

       call ptb_int(ngeom_reconCG, &
         nmat,sdim,mofdata_cell,dx,xsten_cell, &
         int_facefrac_cell, &
         int_centroid_cell, &
         int_face_normal_cell, &
         dist_to_int_cell, &
         xclosest) ! absolute coord

       do im1=1,nmat
       do im2=1,nmat
        if (dist_to_int_cell(im1,im2).lt.0.0) then
         print *,"dist_to_int_cell(im1,im2) invalid"
         stop
        endif
        dist_to_int_FAB(i,j,im1,im2)=dist_to_int_cell(im1,im2)
        do dir=1,sdim
         xclosest_FAB(i,j,im1,im2,dir)=xclosest(im1,im2,dir)
        enddo
       enddo
       enddo

       do im1=1,nmat
       do im2=1,nmat
       do dir=1,sdim
        int_face_normal_FAB(i,j,im1,im2,dir)= &
          int_face_normal_cell(im1,im2,dir)
       enddo
       enddo
       enddo

       call ptb_ext(ngeom_reconCG, &
        nmat,sdim,dx,mofdata_cell, &
        xsten_cell, &
        ext_facefrac_cell, &
        multi_cen_cell)  ! absolute coordinates

       do im=1,nmat+1
       do dir=1,sdim
       do sidesten=1,2
        if (ext_facefrac_cell(im,dir,sidesten).lt.0.0) then
         print *,"ext_facefrac_cell(p,dir,sidesten) invalid"
         stop
        endif
        ext_face_FAB(i,j,im,dir,sidesten)= &
         ext_facefrac_cell(im,dir,sidesten)
       enddo
       enddo
       enddo

       do dircen=1,sdim
       do im=1,nmat
       do dir2=1,sdim
       do sidesten=1,2
        multi_cen_cell_FAB(i,j,dircen,im,dir2,sidesten)= &
         multi_cen_cell(dircen,im,dir2,sidesten)
       enddo
       enddo
       enddo
       enddo
      
        ! in: vfrac_pair.F90 
       call int_face_adjust(nmat,int_facefrac_cell,int_face_cell)

       do im1=1,nmat
       do im2=1,nmat
        if (int_face_cell(im1,im2).lt.0.0) then
         print *,"int_face_cell(im1,im2) invalid"
         stop
        endif
        int_face_FAB(i,j,im1,im2)=int_face_cell(im1,im2)
        do dircen=1,sdim
         int_centroid_FAB(i,j,im1,im2,dircen)= &
           int_centroid_cell(im1,im2,dircen)
        enddo
       enddo
       enddo

      enddo ! j=loy-1,hiy+1
      enddo ! i=lox-1,hix+1

      do i=lox,hix
       do im=1,nmat+1
       do dir=1,sdim
       do sidesten=1,2
        ext_face_FAB(i,loy-1,im,dir,sidesten)= &
          ext_face_FAB(i,loy,im,dir,sidesten)
        ext_face_FAB(i,hiy+1,im,dir,sidesten)= &
          ext_face_FAB(i,hiy,im,dir,sidesten)
       enddo
       enddo
       enddo
       do im1=1,nmat
       do im2=1,nmat
        int_face_FAB(i,loy-1,im1,im2)=int_face_FAB(i,loy,im1,im2)
        int_face_FAB(i,hiy+1,im1,im2)=int_face_FAB(i,hiy,im1,im2)
        do dir=1,sdim
         int_centroid_FAB(i,loy-1,im1,im2,dir)= &
           int_centroid_FAB(i,loy,im1,im2,dir)
         int_centroid_FAB(i,hiy+1,im1,im2,dir)= &
           int_centroid_FAB(i,hiy,im1,im2,dir)
        enddo
       enddo
       enddo
      enddo ! i

      do j=loy-1,hiy+1
       do im=1,nmat+1
       do dir=1,sdim
       do sidesten=1,2
        ext_face_FAB(lox-1,j,im,dir,sidesten)= &
         ext_face_FAB(lox,j,im,dir,sidesten)
        ext_face_FAB(hix+1,j,im,dir,sidesten)= &
         ext_face_FAB(hix,j,im,dir,sidesten)
       enddo
       enddo
       enddo
       do im1=1,nmat
       do im2=1,nmat
        int_face_FAB(lox-1,j,im1,im2)=int_face_FAB(lox,j,im1,im2)
        int_face_FAB(hix+1,j,im1,im2)=int_face_FAB(hix,j,im1,im2)
        do dir=1,sdim
         int_centroid_FAB(lox-1,j,im1,im2,dir)= &
           int_centroid_FAB(lox,j,im1,im2,dir)
         int_centroid_FAB(hix+1,j,im1,im2,dir)= &
           int_centroid_FAB(hix,j,im1,im2,dir)
        enddo
       enddo
       enddo
      enddo

      do i=lox,hix
      do j=loy,hiy

       VP_i_current=i 
       VP_j_current=j 

       do dir = 1,sdim

        do ii = -1,1
         if (dir .eq. 1) then   
          do dircen=1,sdim
          do im=1,nmat
          do dir2=1,sdim
          do sidesten=1,2
           thin_cen_sten(ii,dir,dircen,im,dir2,sidesten)= &
             multi_cen_cell_FAB(i+ii,j,dircen,im,dir2,sidesten)
          enddo
          enddo
          enddo
          enddo

          do im=1,nmat+1
          do dir2=1,sdim
          do sidesten=1,2
           ext_face_sten(ii,dir,im,dir2,sidesten) =  &
             ext_face_FAB(i+ii,j,im,dir2,sidesten)
          enddo
          enddo
          enddo

         elseif(dir .eq. 2)then

          do dircen=1,sdim
          do im=1,nmat
          do dir2=1,sdim
          do sidesten=1,2
           thin_cen_sten(ii,dir,dircen,im,dir2,sidesten)= &
            multi_cen_cell_FAB(i,j+ii,dircen,im,dir2,sidesten)
          enddo
          enddo
          enddo
          enddo
          do im=1,nmat+1
          do dir2=1,sdim
          do sidesten=1,2
           ext_face_sten(ii,dir,im,dir2,sidesten) = &
             ext_face_FAB(i,j+ii,im,dir2,sidesten)
          enddo
          enddo
          enddo
         endif
        enddo   ! ii
       enddo   ! dir

       call vfrac_pair_cell( &
        nmat,sdim,dx, &
        ext_face_sten, &
        thin_cen_sten,  &  ! absolute coordinates
        frac_pair_cell, &
        x_pair_cell)

       do im_outside=1,nmat
       do im_inside=1,nmat
       do dir=1,sdim
       do sidesten=1,2
        frac_pair_cell_FAB(i,j,im_outside,im_inside,dir,sidesten)= &
          frac_pair_cell(im_outside,im_inside,dir,sidesten)
        do dir2=1,sdim
         x_pair_cell_FAB(i,j,im_outside,im_inside,dir2,dir,sidesten)= &
          x_pair_cell(im_outside,im_inside,dir2,dir,sidesten)
        enddo
       enddo
       enddo
       enddo
       enddo
      enddo
      enddo ! i,j

      do i=lox,hix
       do im_outside=1,nmat
       do im_inside=1,nmat
       do dir=1,sdim
       do sidesten=1,2
        frac_pair_cell_FAB(i,loy-1,im_outside,im_inside,dir,sidesten)= &
         frac_pair_cell_FAB(i,loy,im_outside,im_inside,dir,sidesten)
        frac_pair_cell_FAB(i,hiy+1,im_outside,im_inside,dir,sidesten)= &
         frac_pair_cell_FAB(i,hiy,im_outside,im_inside,dir,sidesten)  
        do dir2=1,sdim
         x_pair_cell_FAB(i,loy-1,im_outside,im_inside,dir2,dir,sidesten)= &
          x_pair_cell_FAB(i,loy,im_outside,im_inside,dir2,dir,sidesten)
         x_pair_cell_FAB(i,hiy+1,im_outside,im_inside,dir2,dir,sidesten)= &
          x_pair_cell_FAB(i,hiy,im_outside,im_inside,dir2,dir,sidesten)
        enddo
       enddo
       enddo
       enddo
       enddo
      enddo

      do j=loy-1,hiy+1
       do im_outside=1,nmat
       do im_inside=1,nmat
       do dir=1,sdim
       do sidesten=1,2
        frac_pair_cell_FAB(lox-1,j,im_outside,im_inside,dir,sidesten)= &
         frac_pair_cell_FAB(lox,j,im_outside,im_inside,dir,sidesten)
        frac_pair_cell_FAB(hix+1,j,im_outside,im_inside,dir,sidesten)= &
         frac_pair_cell_FAB(hix,j,im_outside,im_inside,dir,sidesten) 
        do dir2=1,sdim
         x_pair_cell_FAB(lox-1,j,im_outside,im_inside,dir2,dir,sidesten)= &
          x_pair_cell_FAB(lox,j,im_outside,im_inside,dir2,dir,sidesten)
         x_pair_cell_FAB(hix+1,j,im_outside,im_inside,dir2,dir,sidesten)= &
          x_pair_cell_FAB(hix,j,im_outside,im_inside,dir2,dir,sidesten) 
        enddo
       enddo
       enddo
       enddo
       enddo
      enddo

      VP_max_LS_error=0.0d0

      do i=lox,hix
      do j=loy,hiy

       VP_i_current=i
       VP_j_current=j

       do imof=1,ngeom_reconCG*nmat
        mofdata_cell(imof)=mofdata_FAB(i,j,imof)
       enddo
       do ii=-1,1
       do jj=-1,1
        do imof=1,ngeom_reconCG*nmat
         mofdata_sten(ii,jj,imof)=mofdata_FAB(i+ii,j+jj,imof)
        enddo
        do im=1,nmat
        do dir=1,sdim
         mat_cen_sten(ii,jj,im,dir) =  &
           centroid_mult_FAB(i+ii,j+jj,im,dir)
        enddo
        enddo
       enddo
       enddo

       do im_outside=1,nmat
       do im_inside=1,nmat
       do dir=1,sdim
       do sidesten=1,2
        frac_pair_cell(im_outside,im_inside,dir,sidesten)= &
         frac_pair_cell_FAB(i,j,im_outside,im_inside,dir,sidesten)
        do dir2=1,sdim
         x_pair_cell(im_outside,im_inside,dir2,dir,sidesten)= &
          x_pair_cell_FAB(i,j,im_outside,im_inside,dir2,dir,sidesten)
        enddo
       enddo
       enddo
       enddo
       enddo

       do ii=-1,1
       do jj=-1,1

        do im1=1,nmat
        do im2=1,nmat
         dist_to_int_sten(ii,jj,im1,im2)= &
           dist_to_int_FAB(i+ii,j+jj,im1,im2)
         do dir=1,sdim
          xclosest_sten(ii,jj,im1,im2,dir)= &
           xclosest_FAB(i+ii,j+jj,im1,im2,dir)
          int_centroid_sten(ii,jj,im1,im2,dir)= &
           int_centroid_FAB(i+ii,j+jj,im1,im2,dir)
         enddo
         int_face_sten(ii,jj,im1,im2)=int_face_FAB(i+ii,j+jj,im1,im2)
        enddo ! im2
        enddo ! im1

       enddo ! jj
       enddo ! ii

       do im1=1,nmat
       do im2=1,nmat
       do dir=1,sdim
        int_face_normal_cell(im1,im2,dir)= &
          int_face_normal_FAB(i,j,im1,im2,dir)
       enddo
       enddo
       enddo
       do isten=-3,3
       do dir=1,sdim
        xsten_cell(isten,dir)=xsten_FAB(i,j,isten,dir)
       enddo
       enddo

       do im_in=1,nmat
        call cell_div_cal_simple( &
         dist_concentric, &
         diag_coeff_flag, &
         local_linear_exact, &
         operator_internal, &
         operator_external, &
         local_hflag, &
         ngeom_reconCG, &
         sdim,nmat, &
         time_placeholder, &
         dx, &
         xsten_cell, &
         mofdata_sten, &
         rho_box, &
         alpha,  &
         mat_cen_sten, &
         im_in, &
         frac_pair_cell, &
         x_pair_cell, &
         int_face_sten, &
         int_centroid_sten, &
         int_face_normal_cell, &
         dist_to_int_sten, &
         xclosest_sten, & ! absolute coord
         diag_local)
        vofcomp=(im_in-1)*ngeom_reconCG+1
        vf=mofdata_FAB(i,j,vofcomp)
        if (vf.le.VOFTOL_local) then
         diag_local=meshvol/deltat
        else
         diag_local=(vf*meshvol/deltat)+diag_local
        endif
        if ((operator_internal.ge.1).and. &
            (operator_internal.le.3).and. &
            (operator_external.ge.1).and. &
            (operator_external.le.3)) then
         DIAG_FIELD(i,j,im_in)=diag_local
        else if ((operator_internal.eq.0).and. &
                 (operator_external.eq.0)) then
         ! do nothing
        else
         print *,"operator_internal or operator_external invalid"
         stop
        endif
       enddo ! im_in

      enddo
      enddo

      if (1.eq.1) then
       print *,"nmat=",nmat
       print *,"dx(1),dx(2) ",dx(1),dx(2)
       print *,"VP_max_LS_error ",VP_max_LS_error
       print *,"VP_i_max ",VP_i_max
       print *,"VP_j_max ",VP_j_max
       print *,"x=",(VP_i_max+0.5d0)*dx(1)
       print *,"y=",(VP_j_max+0.5d0)*dx(2)
       print *,"VP_xcentroid ",VP_centroid_face(1)
       print *,"VP_ycentroid ",VP_centroid_face(2)
       print *,"VP_areaface_max ",VP_areaface_max
       print *,"VP_dir_max ",VP_dir_max
       print *,"VP_side_max ",VP_side_max
       print *,"VP_im_in_max ",VP_im_in_max
       print *,"VP_im_out_max ",VP_im_out_max
       print *,"VP_test_number_max ",VP_test_number_max
!       stop
      endif

      return
      end subroutine build_MM_diag

      subroutine get_filament_source_local(xpoint,time_source,im,G_in)
      USE mmat_FVM

      IMPLICIT NONE

      integer, intent(in):: im
      REAL*8, intent(in) :: time_source
      REAL*8, intent(out) :: G_in 
      REAL*8, intent(in) :: xpoint(sdim)
     
      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im invalid in get_filament_source_local im,nmat: ", &
        im,nmat  
      endif

      if (1.eq.0) then
       print *,"probtypeCG,im,sdim in get_filament_source: ", &
        probtypeCG,im,sdim
      endif
      call get_filament_source(xpoint,time_source,probtypeCG, &
       im,sdim,G_in)

      return
      end subroutine get_filament_source_local
 
      subroutine INIT_GLOBALS( &
        local_state_ncomp, &
        local_operator_internal, &
        local_operator_external, &
        local_linear_exact, &
        probtype_in, &
        sdim_in,ngeom_recon_in, &
        nx_in,ny_in,lox_in,loy_in,hix_in,hiy_in, &
        UNEW_in,UOLD_in, &
        beta_in,h_in,precond_type_in,bicgstab_tol_in, &
        VFRAC_MOF_in,nmat_in,alpha_in,deltat_in,&
        mofdata_FAB_in, current_time_in)
      IMPLICIT NONE

      integer, intent(in) :: local_operator_internal
      integer, intent(in) :: local_operator_external
      integer, intent(in) :: local_linear_exact
      integer, intent(in) :: local_state_ncomp
      integer, intent(in) :: probtype_in
      integer isten
      integer sdim_in,ngeom_recon_in
      integer nx_in,ny_in,lox_in,loy_in,hix_in,hiy_in
      integer nmat_in,precond_type_in
      integer im,i,j
      REAL*8 deltat_in,bicgstab_tol_in
      REAL*8 h_in
      REAL*8 alpha_in(nmat_in)
      REAL*8 beta_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1,nmat_in)
      REAL*8 UNEW_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1, &
        local_state_ncomp)
      REAL*8 UOLD_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1, &
        local_state_ncomp)
      REAL*8 VFRAC_MOF_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1,nmat_in)

      REAL*8 DIAGCOEFF,betahalf,UWALL,coeff2
      integer dir
      integer sidesten
      integer bctype
      integer imof
      integer basecomp

      REAL*8 vf
      REAL*8 xsrc(sdim_in)
      REAL*8 xgrid,ygrid,xbc,ybc

       ! ngeom_recon_in=vfrac,centroid,order,slope,intercept=2*sdim+3
      REAL*8 mofdata_FAB_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1, &
       ngeom_recon_in*nmat_in) 
      REAL*8 current_time_in,tn,GSRC

      operator_internal=local_operator_internal
      operator_external=local_operator_external
      linear_exact=local_linear_exact

      print *,"in INIT_GLOBALS: probtype_in= ",probtype_in
      probtypeCG=probtype_in
      print *,"in INIT_GLOBALS: probtype_in, probtypeCG= ", &
        probtype_in,probtypeCG

      state_ncomp=local_state_ncomp
      sdim=sdim_in
      ngeom_reconCG=2*sdim+3
      if (ngeom_reconCG.ne.ngeom_recon_in) then
       print *,"ngeom_recon_in invalid"
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif
      nsmooth=8

      nx=nx_in
      ny=ny_in
      lox=lox_in
      loy=loy_in
      hix=hix_in
      hiy=hiy_in
      nz=0
      loz=0
      hiz=0

      precond_type=precond_type_in
      deltat=deltat_in
      bicgstab_tol=bicgstab_tol_in
      h=h_in      
      meshvol=h*h

      nmat=nmat_in

      if (state_ncomp.eq.nmat+global_nten*sdim+ &
                         ngeom_reconCG*nmat+nmat*(sdim+1)) then
       ! do nothing
      else
       print *,"state_ncomp invalid"
       stop
      endif

      do dir=1,sdim
       dx(dir)=h
      enddo
 
      current_time = current_time_in 

      if (deltat.le.0.0) then
       print *,"deltat invalid"
       stop
      endif
      if (current_time.lt.0.0) then
       print *,"current_time invalid"
       stop
      endif
      tn=current_time
      tnp1=current_time+deltat
      tn_mid=0.5d0*(tn+tnp1)

      time_source_global=tnp1

      do im=1,nmat
       alpha(im)=alpha_in(im)
      enddo

      allocate(beta(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(UNEW(lox-1:hix+1,loy-1:hiy+1,state_ncomp)) 
      allocate(UOLD(lox-1:hix+1,loy-1:hiy+1,state_ncomp)) 
      allocate(VFRAC_MOF(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(G(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(DIAG_FIELD(lox-1:hix+1,loy-1:hiy+1,nmat)) 

      allocate(mofdata_FAB(lox-1:hix+1,loy-1:hiy+1, &
        ngeom_reconCG*nmat))

      allocate(xsten_FAB(lox-1:hix+1,loy-1:hiy+1,-3:3,sdim))
      allocate(int_face_FAB(lox-1:hix+1,loy-1:hiy+1,nmat,nmat))
      allocate(int_centroid_FAB(lox-1:hix+1,loy-1:hiy+1,nmat,nmat,sdim))
      allocate(int_face_normal_FAB(lox-1:hix+1,loy-1:hiy+1, &
        nmat,nmat,sdim))
      allocate(dist_to_int_FAB(lox-1:hix+1,loy-1:hiy+1,nmat,nmat))
      allocate(xclosest_FAB(lox-1:hix+1,loy-1:hiy+1,nmat,nmat,sdim))
      allocate(ext_face_FAB(lox-1:hix+1,loy-1:hiy+1,nmat+1,sdim,2))
      allocate(multi_cen_cell_FAB(lox-1:hix+1,loy-1:hiy+1, &
        sdim,nmat,sdim,2))
      allocate(frac_pair_cell_FAB(lox-1:hix+1,loy-1:hiy+1, &
        nmat,nmat,sdim,2))
      allocate(x_pair_cell_FAB(lox-1:hix+1,loy-1:hiy+1, &
        nmat,nmat,sdim,sdim,2))
      allocate(centroid_mult_FAB(lox-1:hix+1,loy-1:hiy+1,nmat,sdim))

      do i=lox-1,hix+1
      do j=loy-1,hiy+1
      do im=1,nmat
       DIAG_FIELD(i,j,im)=1.0E+10
      enddo 
      enddo 
      enddo 

      do i=lox-1,hix+1
      do j=loy-1,hiy+1

       do isten=-3,3
        xsten_FAB(i,j,isten,1)=(i+0.5)*h+isten*h*0.5
        xsten_FAB(i,j,isten,2)=(j+0.5)*h+isten*h*0.5
       enddo

       do im=1,state_ncomp
        UNEW(i,j,im)=UNEW_in(i,j,im)
        UOLD(i,j,im)=UOLD_in(i,j,im)
       enddo

       do im=1,nmat
        xgrid=xsten_FAB(i,j,0,1)
        ygrid=xsten_FAB(i,j,0,2)
      
        beta(i,j,im)=beta_in(i,j,im)
        VFRAC_MOF(i,j,im)=VFRAC_MOF_in(i,j,im)
        vf=VFRAC_MOF(i,j,im)

        basecomp=ngeom_reconCG*(im-1)
        do imof=1,ngeom_reconCG 
         mofdata_FAB(i,j,basecomp+imof) = &
          mofdata_FAB_in(i,j,basecomp+imof)
        enddo
        do dir=1,sdim
         centroid_mult_FAB(i,j,im,dir)= &
          mofdata_FAB(i,j,basecomp+1+dir)+xsten_FAB(i,j,0,dir)
        enddo

        do dir=1,sdim
         xsrc(dir)=centroid_mult_FAB(i,j,im,dir)
        enddo

        if (abs(vf).le.VOFTOL_local) then
         G(i,j,im)=(meshvol/deltat)*UOLD(i,j,im)
        else if ((vf.ge.VOFTOL_local).and.(vf.le.1.0+VOFTOL_local)) then
         if (1.eq.0) then
          print *,"probtypeCG,im,sdim before source: ", &
            probtypeCG,im,sdim
         endif 
         if ((vf.ge.VOFTOL_local).and. &
             (vf.le.1.0d0+VOFTOL_local)) then
          call get_filament_source_local(xsrc,time_source_global, &
            im,GSRC)
         else if ((vf.ge.0.0d0).and.(vf.le.VOFTOL_local)) then
          GSRC=0.0d0
         else
          print *,"vf became corrupt"
          stop
         endif
          
         G(i,j,im)=(meshvol/deltat)*UOLD(i,j,im)+meshvol*GSRC
         if ((operator_internal.ge.1).and. &
             (operator_internal.le.3).and. &
             (operator_external.ge.1).and. &
             (operator_external.le.3)) then
          G(i,j,im)=G(i,j,im)*vf
         else if ((operator_internal.eq.0).and. &
                  (operator_external.eq.0)) then
          ! do nothing
         else
          print *,"operator_internal or operator_external invalid"
          stop
         endif

        else
         print *,"vf invalid"
         stop
        endif
 
        DIAGCOEFF=meshvol/deltat

        if ((i.ge.lox).and.(i.le.hix).and. &
            (j.ge.loy).and.(j.le.hiy)) then

         xbc=xgrid
         ybc=ygrid
         betahalf=2.0*beta_in(i-1,j,im)*beta_in(i,j,im)/ &
          (beta_in(i-1,j,im)+beta_in(i,j,im))
         DIAGCOEFF=DIAGCOEFF+betahalf
         if (i.gt.lox) then
          ! do nothing
         else
          dir=1
          sidesten=1
          bctype=physbc(dir,sidesten)
          UWALL=physbc_value(dir,sidesten)
          xbc=0.0
          call WALLBC_DIAG(xbc,ybc,coeff2,bctype,dir,sidesten)
          DIAGCOEFF=DIAGCOEFF-betahalf*coeff2
         endif

         betahalf=2.0*beta_in(i+1,j,im)*beta_in(i,j,im)/ &
          (beta_in(i+1,j,im)+beta_in(i,j,im))
         DIAGCOEFF=DIAGCOEFF+betahalf
         if (i.lt.hix) then
          ! do nothing
         else
          dir=1
          sidesten=2
          bctype=physbc(dir,sidesten)
          UWALL=physbc_value(dir,sidesten)
          xbc=nx*h
          call WALLBC_DIAG(xbc,ybc,coeff2,bctype,dir,sidesten)
          DIAGCOEFF=DIAGCOEFF-betahalf*coeff2
         endif

         xbc=xgrid
         ybc=ygrid

         betahalf=2.0*beta_in(i,j-1,im)*beta_in(i,j,im)/ &
          (beta_in(i,j-1,im)+beta_in(i,j,im))
 
         DIAGCOEFF=DIAGCOEFF+betahalf
         if (j.gt.loy) then
          ! do nothing
         else
          dir=2
          sidesten=1
          bctype=physbc(dir,sidesten)
          UWALL=physbc_value(dir,sidesten)
          ybc=0.0
          call WALLBC_DIAG(xbc,ybc,coeff2,bctype,dir,sidesten)
          DIAGCOEFF=DIAGCOEFF-betahalf*coeff2
         endif

         betahalf=2.0*beta_in(i,j+1,im)*beta_in(i,j,im)/ &
          (beta_in(i,j+1,im)+beta_in(i,j,im))
         DIAGCOEFF=DIAGCOEFF+betahalf
         if (j.lt.hiy) then
          ! do nothing
         else
          dir=2
          sidesten=2
          bctype=physbc(dir,sidesten)
          UWALL=physbc_value(dir,sidesten)
          ybc=ny*h
          call WALLBC_DIAG(xbc,ybc,coeff2,bctype,dir,sidesten)
          DIAGCOEFF=DIAGCOEFF-betahalf*coeff2
         endif
        endif ! i,j in the interior

        if (DIAGCOEFF.le.0.0) then
         print *,"matrix should be positive definite"
         print *,"DIAGCOEFF=",DIAGCOEFF
         print *,"i,j,im,xgrid,ygrid,vf ",i,j,im,xgrid,ygrid,vf
         print *,"beta= ",beta(i,j,im)
         print *,"h,meshvol ",h,meshvol
         print *,"dx= ",dx(1),dx(2)
         stop
        endif
        DIAG_FIELD(i,j,im)=DIAGCOEFF
       enddo ! im
      enddo
      enddo

      call build_MM_diag()

      return
      end subroutine INIT_GLOBALS

      subroutine WALLBC(xpoint,ypoint,UGHOST,UIN,UWALL,coeff1,coeff2, &
        hflag,bctype,dir,sidesten)
      IMPLICIT NONE

      integer dir,sidesten
      REAL*8 xpoint,ypoint
      integer hflag,bctype
      REAL*8 UGHOST,UWALL,UIN,coeff1,coeff2

      if (bctype.eq.REFLECT_EVEN) then
       UGHOST=UIN
       coeff1=0.0
       coeff2=1.0
      else if (bctype.eq.EXT_DIR) then
       if (hflag.eq.1) then
        UGHOST=-UIN
        coeff1=0.0
        coeff2=-1.0
       else if (hflag.eq.0) then
        UGHOST=2.0*UWALL-UIN
        coeff1=2.0*UWALL
        coeff2=-1.0
       else
        print *,"hflag invalid"
        stop
       endif
      else
       print *,"bctype invalid bctype=",bctype
       stop
      endif
      
      return
      end subroutine WALLBC

      subroutine set_boundary(U,hflag,data_ncomp)
      IMPLICIT NONE

      integer i,j,dir,sidesten,im,data_ncomp
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,data_ncomp)
      REAL*8 xbc,ybc
      REAL*8 coeff1,coeff2
      integer hflag
      integer ilo,ihi,jlo,jhi,ii,jj
      integer bctype
      REAL*8 UWALL

      do im=1,nmat

      do dir=1,2
      do sidesten=1,2
       bctype=physbc(dir,sidesten)
       UWALL=physbc_value(dir,sidesten)
        ! xlo
       if ((dir.eq.1).and.(sidesten.eq.1)) then
        ilo=lox-1
        ihi=lox-1
        jlo=loy
        jhi=hiy
        ii=1
        jj=0
        ! xhi
       else if ((dir.eq.1).and.(sidesten.eq.2)) then
        ilo=hix+1
        ihi=hix+1
        jlo=loy
        jhi=hiy
        ii=-1
        jj=0
        ! jlo
       else if ((dir.eq.2).and.(sidesten.eq.1)) then
        jlo=loy-1
        jhi=loy-1
        ilo=lox
        ihi=hix
        ii=0
        jj=1
        ! jhi
       else if ((dir.eq.2).and.(sidesten.eq.2)) then
        jlo=hiy+1
        jhi=hiy+1
        ilo=lox
        ihi=hix
        ii=0
        jj=-1
       else
        print *,"dir or sidesten invalid"
        stop
       endif

       do i=ilo,ihi
       do j=jlo,jhi
        xbc=(i+0.5)*h+0.5*ii*h
        ybc=(j+0.5)*h+0.5*jj*h 
        call WALLBC(xbc,ybc,U(i,j,im),U(i+ii,j+jj,im),UWALL, &
         coeff1,coeff2, &
         hflag,bctype,dir,sidesten)
       enddo 
       enddo 
      enddo 
      enddo 

      enddo ! im

      return
      end subroutine set_boundary

      subroutine DOTPROD(U,V,dsum)
      IMPLICIT NONE

      integer i,j,im
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 V(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 dsum

      dsum=0.0
      do im=1,nmat
      do i=lox,hix
      do j=loy,hiy
       if (VFRAC_MOF(i,j,im).ge.VOFTOL_local) then
        dsum=dsum+U(i,j,im)*V(i,j,im)
       endif
      enddo
      enddo 
      enddo 
  
      return
      end subroutine DOTPROD


      subroutine DOTPROD_wt(U,V,dsum)
      IMPLICIT NONE

      integer i,j,im
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 V(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 dsum,vf

      dsum=0.0
      do im=1,nmat
      do i=lox,hix
      do j=loy,hiy
       vf=VFRAC_MOF(i,j,im)
       if (vf.ge.VOFTOL_local) then
        dsum=dsum+U(i,j,im)*V(i,j,im)*vf
       endif
      enddo
      enddo 
      enddo 
  
      return
      end subroutine DOTPROD_wt

      subroutine NORMPROD_wt(U,dnorm)
      IMPLICIT NONE

      REAL*8 dnorm
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 dsum

      call DOTPROD_wt(U,U,dsum) 
      dnorm=sqrt(dsum)

      return
      end subroutine NORMPROD_wt




      subroutine NORMPROD(U,dnorm)
      IMPLICIT NONE

      REAL*8 dnorm
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 dsum

      call DOTPROD(U,U,dsum) 
      dnorm=sqrt(dsum)

      return
      end subroutine NORMPROD

       ! W=AA*U + BB*V
      subroutine LINCOMB(U,V,W,AA,BB) 
      IMPLICIT NONE

      integer i,j,im
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 V(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 W(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 AA,BB

      do im=1,nmat
      do i=lox,hix
      do j=loy,hiy
       if (VFRAC_MOF(i,j,im).ge.VOFTOL_local) then
        W(i,j,im)=AA*U(i,j,im)+BB*V(i,j,im)
       endif
      enddo
      enddo
      enddo

      return
      end subroutine LINCOMB

        ! V=U
      subroutine COPYVEC(U,V) 
      IMPLICIT NONE

      integer i,j,im
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 V(lox-1:hix+1,loy-1:hiy+1,nmat)
      do im=1,nmat
      do i=lox-1,hix+1
      do j=loy-1,hiy+1
       V(i,j,im)=U(i,j,im)
      enddo
      enddo
      enddo

      return
      end subroutine COPYVEC

      subroutine ZAPVEC(U) 
      IMPLICIT NONE

      integer i,j,im
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      do im=1,nmat
      do i=lox-1,hix+1
      do j=loy-1,hiy+1
       U(i,j,im)=0.0
      enddo
      enddo
      enddo

      return
      end subroutine ZAPVEC

       ! for k=1,2,3 ....
       ! Z^{k}=Z^{k-1}+D^{-1}(R-AZ^{k-1})
       ! e.g. A=D-L-U
       ! D Z^k = D Z^k-1 + R - (D-L-U)Z^k-1
       ! D Z^k - (L+U)Z^k-1=R 
      subroutine JACPRECOND(Z,R,hflag)
      IMPLICIT NONE

      integer i,j,im,iter
      REAL*8 R(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 Z(lox-1:hix+1,loy-1:hiy+1,nmat)
      integer hflag
      integer save_linear_exact

      REAL*8, dimension(:,:,:), allocatable :: ZN
      REAL*8, dimension(:,:,:), allocatable :: ZNP1
      REAL*8, dimension(:,:,:), allocatable :: RSTAR

      allocate(ZN(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(ZNP1(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(RSTAR(lox-1:hix+1,loy-1:hiy+1,nmat)) 

       ! ZN=Z
      call COPYVEC(Z,ZN)

      do iter=1,nsmooth
       call set_boundary(ZN,hflag,nmat)

       save_linear_exact=linear_exact
       linear_exact=0
       call RESID(RSTAR,R,ZN,hflag)
       linear_exact=save_linear_exact

       do im=1,nmat
       do i=lox,hix
       do j=loy,hiy
        ZNP1(i,j,im)=ZN(i,j,im)+(1.0/DIAG_FIELD(i,j,im))*RSTAR(i,j,im)
       enddo
       enddo
       enddo
        ! ZN=ZNP1
       call COPYVEC(ZNP1,ZN)
      enddo  ! iter
        ! Z=ZN
      call COPYVEC(ZN,Z)

      deallocate(ZN) 
      deallocate(ZNP1) 
      deallocate(RSTAR) 

      return
      end subroutine JACPRECOND

      subroutine MAKE_CONSISTENT(U,hflag)
      IMPLICIT NONE

      integer i,j,im,hflag
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 UAVG,FAVG

      call set_boundary(U,hflag,nmat)
      do i=lox-1,hix+1
      do j=loy-1,hiy+1
       UAVG=0.0
       FAVG=0.0
       do im=1,nmat
        if (VFRAC_MOF(i,j,im).gt.VOFTOL_local) then
         UAVG=UAVG+U(i,j,im)*VFRAC_MOF(i,j,im)
         FAVG=FAVG+VFRAC_MOF(i,j,im)
        endif
       enddo
       if (abs(FAVG-1.0).ge.1.0E-4) then
        print *,"FAVG invalid"
        print *,"i,j ",i,j
        print *,"FAVG=",FAVG
        print *,"VFRAC(1): ",VFRAC_MOF(i,j,1)
        print *,"VFRAC(2): ",VFRAC_MOF(i,j,2)
        stop
       endif
       UAVG=UAVG/FAVG
       do im=1,nmat
        if (VFRAC_MOF(i,j,im).le.VOFTOL_local) then
         U(i,j,im)=UAVG
        endif 
       enddo

      enddo 
      enddo 

      return
      end subroutine MAKE_CONSISTENT


      subroutine ATIMESU(AU,U,hflag)
      use MOF_pair_module
      use mmat_FVM, only: dist_concentric

      IMPLICIT NONE

      integer i,j,im,ii,jj,imof
      integer dir,dir_pair
      integer isten
      integer im1,im2
      integer im_outside,im_inside
      integer sidesten
      REAL*8 betahalfxlo
      REAL*8 betahalfxhi
      REAL*8 betahalfylo
      REAL*8 betahalfyhi
      REAL*8 DIAGCOEFF
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 AU(lox-1:hix+1,loy-1:hiy+1,nmat)
      integer hflag
      REAL*8 xsten_cell(-3:3,sdim)
      REAL*8 vf
      REAL*8 rho_box(-1:1,-1:1,nmat)
      REAL*8 mofdata_sten(-1:1,-1:1,ngeom_reconCG*nmat)
      REAL*8 mat_cen_sten(-1:1,-1:1,nmat,sdim)
      REAL*8 frac_pair_cell(nmat,nmat,sdim,2)
      REAL*8 x_pair_cell(nmat,nmat,sdim,sdim,2)
      REAL*8 int_face_sten(-1:1,-1:1,nmat,nmat)
      REAL*8 int_centroid_sten(-1:1,-1:1,nmat,nmat,sdim)
      REAL*8 int_face_normal_cell(nmat,nmat,sdim)
      REAL*8 dist_to_int_sten(-1:1,-1:1,nmat,nmat)
      REAL*8 xclosest_sten(-1:1,-1:1,nmat,nmat,sdim)
      REAL*8 div_tot
      integer diag_coeff_flag

      diag_coeff_flag=0

      call MAKE_CONSISTENT(U,hflag)

      call set_boundary(U,hflag,nmat)

      if ((operator_internal.eq.0).and. &
          (operator_external.eq.0)) then

       do im=1,nmat
       do i=lox,hix
       do j=loy,hiy
 
        betahalfxlo=2.0*beta(i-1,j,im)*beta(i,j,im)/ &
         (beta(i-1,j,im)+beta(i,j,im))
        betahalfxhi=2.0*beta(i+1,j,im)*beta(i,j,im)/ &
         (beta(i+1,j,im)+beta(i,j,im))
        betahalfylo=2.0*beta(i,j-1,im)*beta(i,j,im)/ &
         (beta(i,j-1,im)+beta(i,j,im))
        betahalfyhi=2.0*beta(i,j+1,im)*beta(i,j,im)/ &
         (beta(i,j+1,im)+beta(i,j,im))
        DIAGCOEFF=meshvol/deltat+(betahalfxlo+betahalfxhi+ &
                 betahalfylo+betahalfyhi)
        AU(i,j,im)=DIAGCOEFF*U(i,j,im)-(betahalfxlo*U(i-1,j,im)+ &
         betahalfxhi*U(i+1,j,im)+betahalfylo*U(i,j-1,im)+ &
         betahalfyhi*U(i,j+1,im))
       enddo
       enddo
       enddo
      else if ((operator_internal.ge.1).and. &
               (operator_internal.le.3).and. &
               (operator_external.ge.1).and. &
               (operator_external.le.3)) then

       VP_max_LS_error=0.0d0

       do i=lox,hix
       do j=loy,hiy

        VP_i_current=i
        VP_j_current=j

        do isten=-3,3
        do dir=1,sdim
         xsten_cell(isten,dir)=xsten_FAB(i,j,isten,dir)
        enddo
        enddo
        do ii=-1,1
        do jj=-1,1
         do im=1,nmat
          rho_box(ii,jj,im)=U(i+ii,j+jj,im)
         enddo
         do imof=1,ngeom_reconCG*nmat
          mofdata_sten(ii,jj,imof)=mofdata_FAB(i+ii,j+jj,imof)
         enddo
        enddo
        enddo
        do ii=-1,1
        do jj=-1,1
         do im=1,nmat
         do dir=1,sdim
          mat_cen_sten(ii,jj,im,dir) = &
            centroid_mult_FAB(i+ii,j+jj,im,dir)
         enddo
         enddo
        enddo
        enddo
        do im_outside=1,nmat
        do im_inside=1,nmat
        do dir=1,sdim
        do sidesten=1,2
         frac_pair_cell(im_outside,im_inside,dir,sidesten)= &
          frac_pair_cell_FAB(i,j,im_outside,im_inside,dir,sidesten)
         do dir_pair=1,sdim
          x_pair_cell(im_outside,im_inside,dir_pair,dir,sidesten)= &
           x_pair_cell_FAB(i,j,im_outside,im_inside,dir_pair,dir,sidesten)
         enddo
        enddo
        enddo
        enddo
        enddo
        do im1=1,nmat
        do im2=1,nmat
        do dir=1,sdim
         int_face_normal_cell(im1,im2,dir)= &
           int_face_normal_FAB(i,j,im1,im2,dir)
        enddo
        enddo
        enddo
        do ii=-1,1
        do jj=-1,1
         do im1=1,nmat
         do im2=1,nmat
          dist_to_int_sten(ii,jj,im1,im2)= &
            dist_to_int_FAB(i+ii,j+jj,im1,im2)
          do dir=1,sdim
           xclosest_sten(ii,jj,im1,im2,dir)= &
            xclosest_FAB(i+ii,j+jj,im1,im2,dir)
           int_centroid_sten(ii,jj,im1,im2,dir)= &
            int_centroid_FAB(i+ii,j+jj,im1,im2,dir)
          enddo
          int_face_sten(ii,jj,im1,im2)=int_face_FAB(i+ii,j+jj,im1,im2)
         enddo
         enddo
        enddo ! jj
        enddo ! ii

        do im=1,nmat

         call cell_div_cal_simple( &
           dist_concentric, &
           diag_coeff_flag, &
           linear_exact, &
           operator_internal, &
           operator_external, &
           hflag, &
           ngeom_reconCG, &
           sdim,nmat, &
           tnp1, &
           dx, &
           xsten_cell, &
           mofdata_sten, &
           rho_box, &
           alpha, &
           mat_cen_sten, &
           im, &
           frac_pair_cell, &
           x_pair_cell, &
           int_face_sten, &
           int_centroid_sten, & ! absolute coord.
           int_face_normal_cell, &
           dist_to_int_sten, &
           xclosest_sten, &  ! absolute coord.
           div_tot)

         vf=VFRAC_MOF(i,j,im)
         if (vf.le.VOFTOL_local) then
          AU(i,j,im)=(meshvol/deltat)*U(i,j,im)
         else        
          AU(i,j,im)=(meshvol/deltat)*vf*U(i,j,im)+div_tot 
         endif

        enddo ! im
       enddo
       enddo ! i,j

      else
       print *,"operator_internal or operator_external invalid"
       stop
      endif
        
      return
      end subroutine ATIMESU



        ! R=RHS-A U
      subroutine RESID(R,RHS,U,hflag)
      IMPLICIT NONE

      REAL*8 R(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 RHS(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 AA,BB
      integer hflag
      REAL*8, dimension(:,:,:), allocatable :: AU

      allocate(AU(lox-1:hix+1,loy-1:hiy+1,nmat)) 

      call ATIMESU(AU,U,hflag)
      AA=1.0
      BB=-1.0
      call LINCOMB(RHS,AU,R,AA,BB)

      deallocate(AU)

      return
      end subroutine RESID


      subroutine preconditioner(Z,R,hflag)
      IMPLICIT NONE

      REAL*8 R(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 Z(lox-1:hix+1,loy-1:hiy+1,nmat)
      integer hflag

      call ZAPVEC(Z)
      if (precond_type.eq.0) then
       call COPYVEC(R,Z)
      else if (precond_type.eq.1) then
       call JACPRECOND(Z,R,hflag)
      else
       print *,"precond_type invalid" 
       stop
      endif

      return
      end subroutine preconditioner

      subroutine check_fab(U,id)
      IMPLICIT NONE

      integer id
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 biggest,smallest,biginv
      integer i,j,im

      biggest=-1.0E+10
      smallest=1.0E+10
      biginv=0.0
      do im=1,nmat
      do i=lox-1,hix+1
      do j=loy-1,hiy+1
       if (U(i,j,im).gt.biggest) then
        biggest=U(i,j,im)
       endif
       if (U(i,j,im).lt.smallest) then
        smallest=U(i,j,im)
       endif
       if (U(i,j,im).ne.0.0) then
        if (1.0/abs(U(i,j,im)).gt.biginv) then
         biginv=1.0/abs(U(i,j,im))
        endif
       endif
      enddo
      enddo
      enddo

      print *,"id,biggest,smallest ",id,biggest,smallest
      print *,"id,biginv ",id,biginv

      return
      end subroutine check_fab

        ! precond_type=0 M=I, =1 Jacobi
      subroutine bicgstab(U,hflag,iter)
      IMPLICIT NONE

      integer, intent(out) :: iter
      integer, intent(in) :: hflag
      integer i,j,im
      integer maxiter
      integer hflagcg,restart_flag
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 dnorm,dnorm0,dnorm_wt

      REAL*8, dimension(:,:,:), allocatable :: U0
      REAL*8, dimension(:,:,:), allocatable :: V0
      REAL*8, dimension(:,:,:), allocatable :: P0
      REAL*8, dimension(:,:,:), allocatable :: R0
      REAL*8, dimension(:,:,:), allocatable :: U1
      REAL*8, dimension(:,:,:), allocatable :: V1
      REAL*8, dimension(:,:,:), allocatable :: P1
      REAL*8, dimension(:,:,:), allocatable :: R1
      REAL*8, dimension(:,:,:), allocatable :: R0hat
      REAL*8, dimension(:,:,:), allocatable :: UINIT
      REAL*8, dimension(:,:,:), allocatable :: RHS
      REAL*8, dimension(:,:,:), allocatable :: Y
      REAL*8, dimension(:,:,:), allocatable :: Hvec
      REAL*8, dimension(:,:,:), allocatable :: S
      REAL*8, dimension(:,:,:), allocatable :: T
      REAL*8, dimension(:,:,:), allocatable :: Z

      REAL*8 rho0,w0,rho1,AA,w1,BB,a1,a2

      allocate(U0(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(V0(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(P0(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(R0(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(U1(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(V1(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(P1(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(R1(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(R0hat(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(UINIT(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(RHS(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(Y(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(Hvec(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(S(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(T(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(Z(lox-1:hix+1,loy-1:hiy+1,nmat)) 

      print *,"bicgstab: operator_internal= ",operator_internal
      print *,"bicgstab: operator_external= ",operator_external
      print *,"bicgstab: linear_exact= ",linear_exact

        ! U0=V0=P0=0 
      AA=0.0 
      do im=1,nmat
      do i=lox-1,hix+1
      do j=loy-1,hiy+1
       U0(i,j,im)=AA
       V0(i,j,im)=AA
       P0(i,j,im)=AA
      enddo
      enddo
      enddo

      if (1.eq.0) then
       call check_fab(U0,1)
       call check_fab(G,2)
       call check_fab(beta,3)
      endif

       ! R0=G-A U0
      call RESID(R0,G,U0,hflag)

      if (1.eq.0) then
       call check_fab(R0,4)
      endif
       ! R0hat=R0
      call COPYVEC(R0,R0hat)
       ! UINIT=U0
      call COPYVEC(U0,UINIT)
      call ZAPVEC(U0)
       ! RHS=R0
      call COPYVEC(R0,RHS)

       ! rho0=AA=w0=1
      rho0=1.0
      AA=1.0
      w0=1.0

      call NORMPROD(R0,dnorm0)
      print *,"initial,dnorm0 ",dnorm0
      dnorm=1.0
      dnorm_wt=1.0
      maxiter=2000
      iter=0

      hflagcg=1
      do while ((dnorm_wt.gt.bicgstab_tol).and.(iter.lt.maxiter))
       print *,"iter,dnorm ",iter,dnorm
       print *,"iter,dnorm_wt ",iter,dnorm_wt

         ! rho1= R0hat^H R0
       call DOTPROD(R0hat,R0,rho1)
       if (1.eq.0) then
        print *,"rho1=",rho1
       endif
       
       restart_flag=0
       if ((sqrt(abs(rho0)).lt.bicgstab_tol*0.01).or. &
           (sqrt(abs(w0)).lt.bicgstab_tol*0.01)) then
        restart_flag=1
       endif 

       if (restart_flag.eq.0) then
          ! (R0hat^H R0)/(R0hat^H dot R0_before)  *   (AA/w0)
        BB=(rho1/rho0)*(AA/w0)
        a1=1.0
        a2=-w0
 
         ! P1=P0-w0 V0
        call LINCOMB(P0,V0,P1,a1,a2)
         ! P1=R0+BB P1
        call LINCOMB(R0,P1,P1,a1,BB)
        ! Y=M^{-1}P1
        call preconditioner(Y,P1,hflagcg)
         
         ! V1=A Y
        call ATIMESU(V1,Y,hflagcg)

         ! AA=rho1/R0hat dot V1
        call DOTPROD(R0hat,V1,AA)

        if (sqrt(abs(AA)).lt.bicgstab_tol*0.01) then
         restart_flag=1
        endif

        if (restart_flag.eq.0) then
         AA=rho1/AA

         ! Hvec=U0+AA Y
         a1=1.0
         a2=AA
         call LINCOMB(U0,Y,Hvec,a1,a2) 
         ! U1=Hvec
         call COPYVEC(Hvec,U1)
         ! R1=RHS-A U1
         call RESID(R1,RHS,U1,hflagcg)
         call NORMPROD(R1,dnorm)
         call NORMPROD_wt(R1,dnorm_wt)
         dnorm=dnorm/dnorm0
         dnorm_wt=dnorm_wt/dnorm0

         if (dnorm_wt.gt.bicgstab_tol) then
          ! S=R0-AA V1
          a1=1.0
          a2=-AA
          call LINCOMB(R0,V1,S,a1,a2) 

           ! Z=M^{-1}S
          call preconditioner(Z,S,hflagcg)

           ! T=A Z
          call ATIMESU(T,Z,hflagcg)

           ! simple case is: (T,S)/(T,T)=(AZ,S)/(AZ,AZ)   (MZ=S)
          call DOTPROD(T,S,a1)
          call DOTPROD(T,T,a2)
          if (sqrt(abs(a2)).lt.bicgstab_tol*0.01) then
           restart_flag=1
          endif

          if (restart_flag.eq.0) then
           w1=a1/a2
           ! U1=Hvec+w1 Z
           a1=1.0
           a2=w1
           call LINCOMB(Hvec,Z,U1,a1,a2) 
          endif
         endif ! dnorm_wt>bicgstab_tol
          ! R1=RHS-A U1
         call RESID(R1,RHS,U1,hflagcg)
         call NORMPROD(R1,dnorm)
         call NORMPROD_wt(R1,dnorm_wt)
         dnorm=dnorm/dnorm0
         dnorm_wt=dnorm_wt/dnorm0
         rho0=rho1
         w0=w1
          ! R0=R1
         call COPYVEC(R1,R0) 
         call COPYVEC(P1,P0) 
         call COPYVEC(V1,V0) 
         call COPYVEC(U1,U0) 
        endif  ! restart_flag=0
       endif ! restart_flag=0

       if (restart_flag.eq.0) then
        ! do nothing
       else if (restart_flag.eq.1) then
        call RESID(R0,RHS,U0,hflagcg) ! R0=RHS-A U0
         ! R0hat=R0
        call COPYVEC(R0,R0hat)
        call COPYVEC(U0,U1) 
         ! rho0=AA=w0=1
        rho0=1.0
        AA=1.0
        w0=1.0
        call NORMPROD(R0,dnorm)
        call NORMPROD_wt(R0,dnorm_wt)
        dnorm=dnorm/dnorm0
        dnorm_wt=dnorm_wt/dnorm0
        call ZAPVEC(V0)
        call ZAPVEC(P0)
       else
        print *,"restart_flag invalid"
        stop
       endif

       iter=iter+1
      enddo

      print *,"at the end: iter,dnorm ",iter,dnorm
      print *,"at the end: iter,dnorm_wt ",iter,dnorm_wt
       ! U=UINIT+U1
      a1=1.0
      a2=a1
      call LINCOMB(UINIT,U1,U,a1,a2)

      call MAKE_CONSISTENT(U,hflag)
 
      deallocate(U0) 
      deallocate(V0) 
      deallocate(P0) 
      deallocate(R0) 
      deallocate(U1) 
      deallocate(V1) 
      deallocate(P1) 
      deallocate(R1) 
      deallocate(R0hat) 
      deallocate(UINIT) 
      deallocate(RHS) 
      deallocate(Y) 
      deallocate(Hvec) 
      deallocate(S) 
      deallocate(T) 
      deallocate(Z) 

      return
      end subroutine bicgstab


      subroutine vinterp(valu,gridval,gridx,gridy,gridz, &
        IV0,IV1,XX,YY,ZZ,icomp)
      IMPLICIT NONE

      REAL*8 valu
      REAL*8 gridx(8),gridy(8),gridz(8),gridval(8)
      integer IV0,IV1
      REAL*8 XX(sdim),YY(sdim),ZZ(sdim)
      integer icomp
      REAL*8 tt

      if (abs(gridval(IV0)-valu).le.1.0D-10) then
       XX(icomp)=gridx(IV0)
       YY(icomp)=gridy(IV0)
#if (BL_SPACEDIM==3)
       ZZ(icomp)=gridz(IV0)
#endif
      else if (abs(gridval(IV1)-valu).le.1.0D-10) then
       XX(icomp)=gridx(IV1)
       YY(icomp)=gridy(IV1)
#if (BL_SPACEDIM==3)
       ZZ(icomp)=gridz(IV1)
#endif
      else
       tt=(gridval(IV0)-valu)/(gridval(IV0)-gridval(IV1))
       if ((tt.lt.0.0).or.(tt.gt.1.0)) then
        print *,"tt invalid"
        stop
       endif
 
       XX(icomp)=tt*gridx(IV1)+(1.0-tt)*gridx(IV0)
       YY(icomp)=tt*gridy(IV1)+(1.0-tt)*gridy(IV0)
#if (BL_SPACEDIM==3)
       ZZ(icomp)=tt*gridz(IV1)+(1.0-tt)*gridz(IV0)
#endif
      endif

      return
      end subroutine vinterp

! in 2d, NORMAL(3)=0
! fictitious node at x1,y1,1
      subroutine addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
      IMPLICIT NONE

      integer itri
      integer imaxtri
      REAL*8 geom(sdim,imaxtri)
      REAL*8 XX(sdim),YY(sdim),ZZ(sdim),NORMAL(sdim)
      REAL*8 CPROD(sdim),VEC1(3),VEC2(3),DOTPROD

#if (BL_SPACEDIM==3)
      VEC1(1)=XX(sdim)-XX(1)
      VEC1(2)=YY(sdim)-YY(1)
      VEC1(3)=ZZ(sdim)-ZZ(1)
#elif (BL_SPACEDIM==2)
      VEC1(1)=0.0
      VEC1(2)=0.0
      VEC1(3)=1.0
#else
      print *,"dimension bust"
      stop
#endif

      VEC2(1)=XX(2)-XX(1)
      VEC2(2)=YY(2)-YY(1)
#if (BL_SPACEDIM==3)
      VEC2(3)=ZZ(2)-ZZ(1)
#elif (BL_SPACEDIM==2)
      VEC2(3)=0.0
#else
      print *,"dimension bust"
      stop
#endif
      CPROD(1)=VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
      CPROD(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
#if (BL_SPACEDIM==3)
      CPROD(sdim)=VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
#endif
      DOTPROD=CPROD(1)*NORMAL(1)+CPROD(2)*NORMAL(2)
#if (BL_SPACEDIM==3)
      DOTPROD=DOTPROD+CPROD(sdim)*NORMAL(sdim)
#endif

      if (itri+sdim.gt.imaxtri-1) then
       print *,"itri too big"
       stop
      endif

      geom(1,itri+1)=XX(1)
      geom(2,itri+1)=YY(1)
#if (BL_SPACEDIM==3)
      geom(sdim,itri+1)=ZZ(1)
#endif
#if (BL_SPACEDIM==3)
      if (DOTPROD.gt.0.0) then
       geom(1,itri+2)=XX(2)
       geom(2,itri+2)=YY(2)
       geom(sdim,itri+2)=ZZ(2)
       geom(1,itri+sdim)=XX(sdim)
       geom(2,itri+sdim)=YY(sdim)
       geom(sdim,itri+sdim)=ZZ(sdim)
      else
       geom(1,itri+2)=XX(sdim)
       geom(2,itri+2)=YY(sdim)
       geom(3,itri+2)=ZZ(sdim)
       geom(1,itri+sdim)=XX(2)
       geom(2,itri+sdim)=YY(2)
       geom(sdim,itri+sdim)=ZZ(2)
      endif
#elif (BL_SPACEDIM==2)
      geom(1,itri+2)=XX(2)
      geom(2,itri+2)=YY(2)
#else
      print *,"dimension bust"
      stop
#endif


      itri=itri+sdim

      return
      end subroutine addgeom


      subroutine polysegment(gridx,gridy,gridz,gridval,valu,geom, &
        itri,IV0,IV1,IV2,IV3,imaxtri)
      use global_utility_module
      IMPLICIT NONE

      integer itri
      integer IV0,IV1,IV2,IV3
      integer imaxtri
      REAL*8 geom(sdim,imaxtri)
      REAL*8 valu
      REAL*8 gridx(8),gridy(8),gridz(8),gridval(8)
      REAL*8 AA(4,4),sourcex(4),bb(4),NORMAL(sdim)
      REAL*8 XX(sdim),YY(sdim),ZZ(sdim)
      integer istat,idxtri

      idxtri=0
      if (gridval(IV0).lt.valu) then
       idxtri=idxtri+1
      endif
      if (gridval(IV1).lt.valu) then
       idxtri=idxtri+2
      endif
      if (gridval(IV2).lt.valu) then
       idxtri=idxtri+4
      endif

      if ((idxtri.ne.7).and.(idxtri.ne.0)) then
       AA(1,1)=gridx(IV0) 
       AA(1,2)=gridy(IV0) 
       AA(1,3)=0.0
       AA(1,4)=1.0
       AA(2,1)=gridx(IV1) 
       AA(2,2)=gridy(IV1) 
       AA(2,3)=0.0
       AA(2,4)=1.0
       AA(3,1)=gridx(IV2) 
       AA(3,2)=gridy(IV2) 
       AA(3,3)=0.0
       AA(3,4)=1.0
       AA(4,1)=gridx(IV0) 
       AA(4,2)=gridy(IV0) 
       AA(4,3)=1.0
       AA(4,4)=1.0
       bb(1)=gridval(IV0)
       bb(2)=gridval(IV1)
       bb(3)=gridval(IV2)
       bb(4)=gridval(IV0)
       call matrix_solve(AA,sourcex,bb,istat,4)
       if (istat.ne.0) then
        NORMAL(1)=sourcex(1)
        NORMAL(2)=sourcex(2)
        if ((idxtri.eq.1).or.(idxtri.eq.6)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV1,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV2,XX,YY,ZZ,2)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.2).or.(idxtri.eq.5)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV0,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV2,XX,YY,ZZ,2)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.3).or.(idxtri.eq.4)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV2,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV2,XX,YY,ZZ,2)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
        endif
       endif ! istat.ne.0.0
      endif ! idxtri.ne.7 or 0

      return
      end subroutine polysegment


      subroutine polytri(gridx,gridy,gridz,gridval,valu,geom, &
        itri,IV0,IV1,IV2,IV3,imaxtri)
      use global_utility_module
      IMPLICIT NONE

      integer itri
      integer IV0,IV1,IV2,IV3
      integer imaxtri
      REAL*8 geom(sdim,imaxtri)
      REAL*8 valu
      REAL*8 gridx(8),gridy(8),gridz(8),gridval(8)
      REAL*8 AA(4,4),sourcex(4),bb(4),NORMAL(sdim)
      REAL*8 XX(sdim),YY(sdim),ZZ(sdim)
      integer istat,idxtri

      idxtri=0
      if (gridval(IV0).lt.valu) then
       idxtri=idxtri+1
      endif
      if (gridval(IV1).lt.valu) then
       idxtri=idxtri+2
      endif
      if (gridval(IV2).lt.valu) then
       idxtri=idxtri+4
      endif
      if (gridval(IV3).lt.valu) then
       idxtri=idxtri+8
      endif

      if ((idxtri.ne.15).and.(idxtri.ne.0)) then
       AA(1,1)=gridx(IV0) 
       AA(1,2)=gridy(IV0) 
       AA(1,3)=gridz(IV0) 
       AA(1,4)=1.0
       AA(2,1)=gridx(IV1) 
       AA(2,2)=gridy(IV1) 
       AA(2,3)=gridz(IV1) 
       AA(2,4)=1.0
       AA(3,1)=gridx(IV2) 
       AA(3,2)=gridy(IV2) 
       AA(3,3)=gridz(IV2) 
       AA(3,4)=1.0
       AA(4,1)=gridx(IV3) 
       AA(4,2)=gridy(IV3) 
       AA(4,3)=gridz(IV3) 
       AA(4,4)=1.0
       bb(1)=gridval(IV0)
       bb(2)=gridval(IV1)
       bb(3)=gridval(IV2)
       bb(4)=gridval(IV3)
       call matrix_solve(AA,sourcex,bb,istat,4)
       if (istat.ne.0) then
        NORMAL(1)=sourcex(1)
        NORMAL(2)=sourcex(2)
        NORMAL(sdim)=sourcex(sdim)
        if ((idxtri.eq.1).or.(idxtri.eq.14)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV1,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV2,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV3,XX,YY,ZZ,3)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.2).or.(idxtri.eq.13)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV0,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV3,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV2,XX,YY,ZZ,3)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.3).or.(idxtri.eq.12)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV3,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV2,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV3,XX,YY,ZZ,3)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
 
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV3,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV2,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV2,XX,YY,ZZ,3)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.4).or.(idxtri.eq.11)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV0,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV1,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV3,XX,YY,ZZ,3)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.5).or.(idxtri.eq.10)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV1,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV3,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV3,XX,YY,ZZ,3)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
  
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV1,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV2,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV3,XX,YY,ZZ,3)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.6).or.(idxtri.eq.9)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV1,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV3,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV3,XX,YY,ZZ,3)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
 
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV1,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV2,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV3,XX,YY,ZZ,3)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.7).or.(idxtri.eq.8)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV3,IV0,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV3,IV2,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV3,IV1,XX,YY,ZZ,3)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
        endif
       endif ! istat.ne.0.0
      endif ! idxtri.ne.15 or 0

      return
      end subroutine polytri



      subroutine isogrid(im,nsteps,out_time,plot_int,total_nsteps)
      IMPLICIT NONE

      integer, intent(in) :: im
      integer, intent(in) :: nsteps
      integer, intent(in) :: total_nsteps
      integer, intent(in) :: plot_int
      integer strandid 
      REAL*8 valu
      REAL*8, intent(in) :: out_time

      character*2 matstr
      character*6 stepstr
      character*3 Mstr

      character*7 newnamestr7
      character*11 newcennamestr11
      character*20 newfilename20
      character*24 newcenfilename24

      real(8), dimension(:,:), allocatable :: Node  ! dir,index
      integer(4), dimension(:,:), allocatable :: IntElem  ! 1 or 2, index
      integer(4) :: NumNodes
      integer(4) :: NumIntElems
      integer(4) :: CurNodes
      integer(4) :: CurIntElems
      integer imaxtri,itri,curtri,ipass
      REAL*8 geom(sdim,200)
      REAL*8 lnode(0:1,0:1,0:1)
      REAL*8 xnode(0:1,0:1,0:1,sdim)
      integer i,j,k,ii,jj,kk,dir,ISUM,itemp,jtemp
      integer kklo,kkhi,nodehi
      integer N8(8)
      REAL*8 gridval(8)
      REAL*8 gridx(8)
      REAL*8 gridy(8)
      REAL*8 gridz(8)

      REAL*8 xtarget(sdim)
      REAL*8 nn(sdim)
      REAL*8 intercept
      REAL*8 xref(sdim)
      integer nparticles
      integer interface_found

      REAL*8 vfrac,vfrac_side
      integer iside,jside,kside
      REAL*8 xoutput(sdim)
      integer ibase,isten
      REAL*8 xsten(-3:3,sdim)

      if (im.lt.1) then
       print *,"im out of range"
       stop
      endif
      if ((nsteps.ge.0).and.(nsteps.le.total_nsteps)) then
       ! do nothing
      else
       print *,"nsteps invalid"
       stop
      endif

      write(Mstr,'(I3)') total_nsteps
      do i=1,3
       if (Mstr(i:i).eq.' ') then
        Mstr(i:i)='0'
       endif
      enddo

      write(matstr,'(I2)') im
      do i=1,2
       if (matstr(i:i).eq.' ') then
        matstr(i:i)='0'
       endif
      enddo

      imaxtri=200
      valu=0.0

      write(newnamestr7,'(A3,A2,A2)') 'mat',matstr,'ls'
      write(newcennamestr11,'(A6,A2,A3)') 'refcen',matstr,'pos'

      write(stepstr,'(I6)') nsteps
      do i=1,6
       if (stepstr(i:i).eq.' ') then
        stepstr(i:i)='0'
       endif
      enddo
      write(newfilename20,'(A7,A6,A3,A4)') &
        newnamestr7,stepstr,Mstr,'.tec'
      print *,"newfilename20 ",newfilename20
      open(unit=11,file=newfilename20)

      if (sdim.eq.3) then
       write(11,*) 'TITLE = "3D surface" '
       write(11,*) 'VARIABLES = "X", "Y", "Z" '
      else if (sdim.eq.2) then
       write(11,*) 'TITLE = "2D surface" '
       write(11,*) 'VARIABLES = "X", "Y" '
      else
       print *,"dimension bust"
       stop
      endif

      NumNodes=0
      NumIntElems=0
      nparticles=0

      do ipass=0,1
       CurNodes=0
       CurIntElems=0

       if (ipass.eq.1) then
        allocate(Node(sdim,NumNodes))
        allocate(IntElem(sdim,NumIntElems))
        write(11,'(A23,I14,A5,I14,A21)')  &
          'ZONE T="TRIANGLES", N= ',NumNodes,  &
          ', E= ',NumIntElems,', DATAPACKING=POINT, '
        if (plot_int.le.0) then
         strandid=1
        else
         strandid=(nsteps/plot_int)+1
        endif
        if (sdim.eq.3) then
         write(11,'(A33,G25.16e3,A10,I10)')  &
          'ZONETYPE=FETRIANGLE SOLUTIONTIME=',out_time, &
          " STRANDID=",strandid
        else if (sdim.eq.2) then
         write(11,'(A32,G25.16e3,A10,I10)')  &
          'ZONETYPE=FELINESEG SOLUTIONTIME=',out_time, &
          " STRANDID=",strandid
        else
         print *,"dimension bust"
         stop
        endif
        write(newcenfilename24,'(A11,A6,A3,A4)') newcennamestr11, &
          stepstr,Mstr,'.tec'
        print *,"newcenfilename24 ",newcenfilename24
        open(unit=12,file=newcenfilename24)
        if (sdim.eq.3) then
         write(12,*) 'TITLE = "3D moments" '
         write(12,*) 'VARIABLES = "X", "Y", "Z" '
        else if (sdim.eq.2) then
         write(12,*) 'TITLE = "2D moments" '
         write(12,*) 'VARIABLES = "X", "Y" '
        else
         print *,"dimension bust"
         stop
        endif
        if (plot_int.le.0) then
         strandid=1
        else
         strandid=(nsteps/plot_int)+1
        endif
        write(12,'(A19,I14,A26,G25.16e3,A10,I10)') & 
          'ZONE F="POINT", I= ', nparticles,  &
          ', J=1, K=1, SOLUTIONTIME= ',out_time,' STRANDID=',strandid

       endif ! ipass=1

       if (sdim.eq.3) then
        kklo=0
        kkhi=1
        nodehi=8
       else if (sdim.eq.2) then
        kklo=0
        kkhi=0
        nodehi=4
       else
        print *,"dimension bust"
        stop
       endif

       do i=lox,hix
       do j=loy,hiy
       do k=loz,hiz

        do dir=1,sdim
        do isten=-3,3
         xsten(isten,dir)=xsten_FAB(i,j,isten,dir)
        enddo 
        enddo 
        ibase=(im-1)*ngeom_reconCG

         ! vfrac,centroid,order,slope,intercept
        do dir=1,sdim
         nn(dir)=mofdata_FAB(i,j,ibase+sdim+2+dir)
        enddo
        intercept=mofdata_FAB(i,j,ibase+2*sdim+3)

        vfrac=mofdata_FAB(i,j,ibase+1)

        itri=0
        do ii=0,1
        do jj=0,1
        do kk=kklo,kkhi

         do dir=1,sdim
          xtarget(dir)=xsten(-1,dir)
         enddo
         if (ii.eq.1) then
          dir=1 
          xtarget(dir)=xsten(1,dir)
         endif
         if (jj.eq.1) then
          dir=2 
          xtarget(dir)=xsten(1,dir)
         endif
         if (kk.eq.1) then
          dir=sdim
          xtarget(dir)=xsten(1,dir)
         endif
         do dir=1,sdim
          xnode(ii,jj,kk,dir)=xtarget(dir)
         enddo

         if (ii.eq.0) then
          iside=i-1
         else if (ii.eq.1) then
          iside=i+1
         else
          print *,"ii invalid"
          stop
         endif
         if (jj.eq.0) then
          jside=j-1
         else if (jj.eq.1) then
          jside=j+1
         else
          print *,"jj invalid"
          stop
         endif
         if (kk.eq.0) then
          kside=k-1
         else if (kk.eq.1) then
          kside=k+1
         else
          print *,"kk invalid"
          stop
         endif
         
         if ((iside.ge.lox).and.(iside.le.hix).and. &
             (jside.ge.loy).and.(jside.le.hiy)) then 
          vfrac_side=mofdata_FAB(iside,jside,ibase+1)
         else
          vfrac_side=vfrac
         endif

         interface_found=0

         if ((vfrac.ge.1.0-VOFTOL_local).and. &
             (vfrac_side.le.VOFTOL_local)) then
          lnode(ii,jj,kk)=-VOFTOL_local
          interface_found=1
         else if ((vfrac.le.VOFTOL_local).and. &
                  (vfrac_side.ge.1.0-VOFTOL_local)) then
          lnode(ii,jj,kk)=VOFTOL_local
          interface_found=1
         else if (vfrac.ge.1.0-VOFTOL_local) then
          lnode(ii,jj,kk)=1.0
         else if (vfrac.le.VOFTOL_local) then
          lnode(ii,jj,kk)=-1.0
         else
          lnode(ii,jj,kk)=intercept
          do dir=1,sdim
           lnode(ii,jj,kk)=lnode(ii,jj,kk)+ &
             nn(dir)*(xtarget(dir)-xsten(0,dir))
          enddo 
          interface_found=1
         endif
        enddo
        enddo
        enddo ! ii,jj,kk

        gridval(1)=lnode(0,0,0)
        gridval(2)=lnode(1,0,0)
        gridval(3)=lnode(1,1,0)
        gridval(4)=lnode(0,1,0)
        gridval(5)=lnode(0,0,kkhi)
        gridval(6)=lnode(1,0,kkhi)
        gridval(7)=lnode(1,1,kkhi)
        gridval(8)=lnode(0,1,kkhi)

        ISUM=0
        do itemp=1,nodehi
         if (gridval(itemp).ge.valu) then
          N8(itemp)=1
         else 
          N8(itemp)=0
         endif
         ISUM=ISUM+N8(itemp)
        enddo

        if ((ISUM.eq.0).or.(ISUM.eq.nodehi)) then
         goto 999
        endif

        if (ipass.eq.0) then
         nparticles=nparticles+1
        else if (ipass.eq.1) then
         do dir=1,sdim
          xref(dir)=xsten(0,dir)+mofdata_FAB(i,j,ibase+1+dir)
         enddo
         if (sdim.eq.3) then
          write(12,*) xref(1),xref(2),xref(sdim)
         else if (sdim.eq.2) then
          write(12,*) xref(1),xref(2)
         else
          print *,"dimension bust"
          stop
         endif
        else
         print *,"ipass invalid"
         stop
        endif
       
        if (interface_found.eq.1) then
         if (1.eq.0) then
          print *,"im,x,y,vf,xcen,ycen ",im,xsten(0,1),xsten(0,2), &
           mofdata_FAB(i,j,ibase+1), &
           mofdata_FAB(i,j,ibase+2),mofdata_FAB(i,j,ibase+3)
         endif
        else if (interface_found.eq.0) then
         ! do nothing
        else
         print *,"interface_found invalid"
         stop
        endif
 
        dir=1
        gridx(1)=xnode(0,0,0,dir)
        gridx(2)=xnode(1,0,0,dir)
        gridx(3)=xnode(1,1,0,dir)
        gridx(4)=xnode(0,1,0,dir)
        gridx(5)=xnode(0,0,kkhi,dir)
        gridx(6)=xnode(1,0,kkhi,dir)
        gridx(7)=xnode(1,1,kkhi,dir)
        gridx(8)=xnode(0,1,kkhi,dir)
        dir=2
        gridy(1)=xnode(0,0,0,dir)
        gridy(2)=xnode(1,0,0,dir)
        gridy(3)=xnode(1,1,0,dir)
        gridy(4)=xnode(0,1,0,dir)
        gridy(5)=xnode(0,0,kkhi,dir)
        gridy(6)=xnode(1,0,kkhi,dir)
        gridy(7)=xnode(1,1,kkhi,dir)
        gridy(8)=xnode(0,1,kkhi,dir)
        dir=sdim
        gridz(1)=xnode(0,0,0,dir)
        gridz(2)=xnode(1,0,0,dir)
        gridz(3)=xnode(1,1,0,dir)
        gridz(4)=xnode(0,1,0,dir)
        gridz(5)=xnode(0,0,kkhi,dir)
        gridz(6)=xnode(1,0,kkhi,dir)
        gridz(7)=xnode(1,1,kkhi,dir)
        gridz(8)=xnode(0,1,kkhi,dir)

        if (sdim.eq.2) then
         call polysegment(gridx,gridy,gridz,gridval,valu,geom, &
          itri,1,2,3,3,imaxtri)
         call polysegment(gridx,gridy,gridz,gridval,valu,geom, &
          itri,1,4,3,3,imaxtri)
        else if (sdim.eq.3) then
         call polytri(gridx,gridy,gridz,gridval,valu,geom, &
          itri,1,3,4,8,imaxtri)
         call polytri(gridx,gridy,gridz,gridval,valu,geom, &
          itri,1,3,7,8,imaxtri)
         call polytri(gridx,gridy,gridz,gridval,valu,geom, &
          itri,1,5,7,8,imaxtri)
         call polytri(gridx,gridy,gridz,gridval,valu,geom, &
          itri,1,7,2,3,imaxtri)
         call polytri(gridx,gridy,gridz,gridval,valu,geom, &
          itri,1,7,2,5,imaxtri)
         call polytri(gridx,gridy,gridz,gridval,valu,geom, &
          itri,6,7,2,5,imaxtri)
        else
         print *,"dimension bust"
         stop
        endif
!
999     continue
 
        if (ipass.eq.0) then
         NumNodes=NumNodes+itri
         NumIntElems=NumIntElems+itri/sdim
        else if (ipass.eq.1) then
         curtri=0
         do itemp=1,itri/sdim
          CurIntElems=CurIntElems+1
          if (CurIntElems.gt.NumIntElems) then
           print *,"CurIntElems invalid"
           stop
          endif
          do jtemp=1,sdim
           CurNodes=CurNodes+1
           if (CurNodes.gt.NumNodes) then
            print *,"CurNodes invalid"
            stop
           endif
           curtri=curtri+1
           if (curtri.gt.itri) then
            print *,"curtri invalid"
            stop
           endif
           do dir=1,sdim
            Node(dir,CurNodes)=geom(dir,curtri)
           enddo
           IntElem(jtemp,CurIntElems)=CurNodes
          enddo
         enddo  
        else
         print *,"ipass invalid"
         stop
        endif
  
       enddo
       enddo
       enddo
      enddo ! ipass

      close(12)

      do i=1,NumNodes 
       do dir=1,sdim
        xoutput(dir)=Node(dir,i)
       enddo
       if (sdim.eq.3) then
        write(11,*) xoutput(1),xoutput(2),xoutput(sdim)
       else if (sdim.eq.2) then
        write(11,*) xoutput(1),xoutput(2)
       else
        print *,"dimension bust"
        stop
       endif
      enddo
      do i=1,NumIntElems
       if (sdim.eq.3) then
        write(11,*) IntElem(1,i),IntElem(2,i),IntElem(sdim,i)
       else if (sdim.eq.2) then
        write(11,*) IntElem(1,i),IntElem(2,i)
       else
        print *,"dimension bust"
        stop
       endif
      enddo
      close(11)

      deallocate(Node)
      deallocate(IntElem)
 
      return
      end subroutine isogrid

      subroutine report_error(local_UNEW,out_time)

        ! declared in multimat_FVM.F90
      REAL(kind=8), external :: exact_temperature
      REAL*8 out_time
      REAL*8 local_UNEW(lox-1:hix+1,loy-1:hiy+1,state_ncomp)
      integer vofcomp,i,j,im,dir,sidesten,isten,ii,jj,dir2
      integer icrit(nmat)
      integer jcrit(nmat)
      integer icrit_gradient(nmat)
      integer jcrit_gradient(nmat)
      REAL*8 vf,vf_side,UEXACT,UEXACT_side,err,mag
      REAL*8 GRD,GRD_EXACT
      REAL*8 linf_error(nmat)
      REAL*8 l1_error(nmat)
      REAL*8 l2_error(nmat)
      REAL*8 error_denom(nmat)
      REAL*8 linf_gradient_error(nmat)
      REAL*8 l1_gradient_error(nmat)
      REAL*8 l2_gradient_error(nmat)
      REAL*8 gradient_error_denom(nmat)
      REAL*8 xref(sdim)
      REAL*8 xref_side(sdim)
      REAL*8 xsten(-3:3,sdim)
      REAL*8 xsten_side(-3:3,sdim)
     
      if (ngeom_reconCG.ne.2*sdim+3) then
       print *,"ngeom_reconCG invalid"
       stop
      endif
      if (h.le.0.0) then
       print *,"h invalid"
       stop
      endif
 
      do im=1,nmat
       icrit(im)=-1
       jcrit(im)=-1
       icrit_gradient(im)=-1
       jcrit_gradient(im)=-1
       linf_error(im)=-1.0
       l1_error(im)=0.0 
       l2_error(im)=0.0 
       error_denom(im)=0.0 
       linf_gradient_error(im)=-1.0
       l1_gradient_error(im)=0.0 
       l2_gradient_error(im)=0.0 
       gradient_error_denom(im)=0.0 
      enddo ! im

      do i=lox,hix
      do j=loy,hiy
      do im=1,nmat
       do dir=1,sdim
       do isten=-3,3
        xsten(isten,dir)=xsten_FAB(i,j,isten,dir)
       enddo 
       enddo 
       vofcomp=(im-1)*ngeom_reconCG+1
       do dir=1,sdim
        xref(dir)=xsten(0,dir)+mofdata_FAB(i,j,vofcomp+dir)
       enddo
       vf=mofdata_FAB(i,j,vofcomp) 
       if (vf.gt.ERRTOL) then
        UEXACT=exact_temperature(xref, &
          out_time,im,probtypeCG,nmat,alpha)
        err=abs(local_UNEW(i,j,im)-UEXACT)
        if (err.gt.linf_error(im)) then
         icrit(im)=i  
         jcrit(im)=j
         linf_error(im)=err
        endif
        error_denom(im)=error_denom(im)+meshvol*vf
        l1_error(im)=l1_error(im)+err*meshvol*vf  
        l2_error(im)=l2_error(im)+(err**2)*meshvol*vf  

        do dir=1,sdim
        do sidesten=-1,1,2
         ii=0
         jj=0
         if (dir.eq.1) then
          ii=sidesten
         else if (dir.eq.2) then
          jj=sidesten
         else
          print *,"dir invalid"
          stop
         endif
         do dir2=1,sdim
         do isten=-3,3
          xsten_side(isten,dir2)=xsten_FAB(i+ii,j+jj,isten,dir2)
         enddo 
         enddo 
         do dir2=1,sdim
          xref_side(dir2)= &
            xsten_side(0,dir2)+mofdata_FAB(i+ii,j+jj,vofcomp+dir2)
         enddo
         vf_side=mofdata_FAB(i+ii,j+jj,vofcomp) 
         if (vf_side.gt.ERRTOL) then
          UEXACT_side=exact_temperature(xref_side, &
           out_time,im,probtypeCG,nmat,alpha)
          mag=0.0
          do dir2=1,sdim
           mag=mag+(xref(dir2)-xref_side(dir2))**2
          enddo
          mag=sqrt(mag)
          GRD_EXACT=(UEXACT_side-UEXACT)/mag
          GRD=(local_UNEW(i+ii,j+jj,im)-local_UNEW(i,j,im))/mag
          err=abs(GRD-GRD_EXACT)
          if (err.gt.linf_gradient_error(im)) then
           icrit_gradient(im)=i
           jcrit_gradient(im)=j
           linf_gradient_error(im)=err
          endif
          gradient_error_denom(im)=gradient_error_denom(im)+meshvol*vf
          l1_gradient_error(im)=l1_gradient_error(im)+err*meshvol*vf  
          l2_gradient_error(im)=l2_gradient_error(im)+ &
             (err**2)*meshvol*vf
         endif
        enddo
        enddo

       endif ! vf.gt.errtol
      enddo ! im
      enddo ! j
      enddo ! i

      do im=1,nmat
       print *,"----------------------------------------------"
       i=icrit(im)
       j=jcrit(im)
       do dir=1,sdim
       do isten=-3,3
        xsten(isten,dir)=xsten_FAB(i,j,isten,dir)
       enddo 
       enddo 
       do dir=1,sdim
        xref(dir)=xsten(0,dir)+mofdata_FAB(i,j,vofcomp+dir)
       enddo
       vofcomp=(im-1)*ngeom_reconCG+1
       vf=mofdata_FAB(i,j,vofcomp) 
       print *,"time,im,x,y,vf,linf ",out_time, &
        im,xref(1),xref(2),vf,linf_error(im)
       print *,"time,im,error_denom ",out_time,im,error_denom(im)
       if (error_denom(im).gt.0.0d0) then
        print *,"time,im,l1 ",out_time,im,l1_error(im)/error_denom(im)
        print *,"time,im,l2 ",out_time,im, &
         sqrt(l2_error(im)/error_denom(im))
       endif

       i=icrit_gradient(im)
       j=jcrit_gradient(im)
       do dir=1,sdim
       do isten=-3,3
        xsten(isten,dir)=xsten_FAB(i,j,isten,dir)
       enddo
       enddo
       do dir=1,sdim
        xref(dir)=xsten(0,dir)+mofdata_FAB(i,j,vofcomp+dir)
       enddo
       vofcomp=(im-1)*ngeom_reconCG+1
       vf=mofdata_FAB(i,j,vofcomp)
       print *,"time,im,x,y,vf,linf_gradient ",out_time, &
        im,xref(1),xref(2),vf,linf_gradient_error(im)
       print *,"time,im,gradient_error_denom ", &
        out_time,im,gradient_error_denom(im)
       if (gradient_error_denom(im).gt.0.0d0) then
        print *,"time,im,l1_gradient ", &
          out_time,im,l1_gradient_error(im)/gradient_error_denom(im)
        print *,"time,im,l2_gradient ", &
         out_time,im, &
         sqrt(l2_gradient_error(im)/gradient_error_denom(im))
       endif

       print *,"----------------------------------------------"
      enddo ! im

      return
      end subroutine report_error


      subroutine output_solution(local_UNEW,out_time,nsteps,plot_int, &
        total_nsteps)
      IMPLICIT NONE

      character*8 nodedatastr
      character*14 nodedatafile1
      character*21 nodedatafile2
      character*7 cendatastr
      character*13 cendatafile1
      character*20 cendatafile2
      REAL*8, intent(in) ::  &
        local_UNEW(lox-1:hix+1,loy-1:hiy+1,state_ncomp)
      integer i,j
      integer, intent(in) :: nsteps
      integer, intent(in) :: total_nsteps
      integer, intent(in) :: plot_int
      integer strandid,im,lscomp
      REAL*8 xpoint,ypoint
      REAL*8, intent(in) :: out_time
      character*6 stepstr
      character*3 Mstr
      character*2 matstr
      integer iten
      integer ii,jj
      REAL*8 sum_wt
      REAL*8 wt_tecplot(0:1,0:1)
      REAL*8 node_data
      REAL*8 cen_data

      if (global_nten.eq.((nmat-1)*(nmat-1)+(nmat-1))/2) then
       ! do nothing
      else
       print *,"global_nten invalid"
       stop
      endif
      if (sdim.eq.2) then
       ! do nothing
      else
       print *,"only 2d supported in the prototype code"
       stop
      endif
      if (state_ncomp.eq. &
          nmat+global_nten*sdim+ngeom_reconCG*nmat+nmat*(sdim+1)) then
       ! do nothing
      else
       print *,"state_ncomp invalid"
       stop
      endif
      if ((nsteps.ge.0).and.(nsteps.le.total_nsteps)) then
       ! do nothing
      else
       print *,"nsteps invalid"
       stop
      endif

      call report_error(local_UNEW,out_time)

       ! hflag=0
      call set_boundary(local_UNEW,0,state_ncomp)

      if (plot_int.gt.0) then
       if (((nsteps/plot_int)*plot_int.eq.nsteps).or. &
           (nsteps.eq.total_nsteps)) then

        do im=1,nmat
         call isogrid(im,nsteps,out_time,plot_int,total_nsteps)
        enddo

        write(stepstr,'(I6)') nsteps
        do i=1,6
         if (stepstr(i:i).eq.' ') then
          stepstr(i:i)='0'
         endif
        enddo

        write(Mstr,'(I3)') total_nsteps
        do i=1,3
         if (Mstr(i:i).eq.' ') then
          Mstr(i:i)='0'
         endif
        enddo

        if (plot_int.le.0) then
         strandid=1
        else
         strandid=(nsteps/plot_int)+1
        endif

        write(nodedatastr,'(A8)') 'nodedata'
        write(nodedatafile1,'(A8,A6)') nodedatastr,stepstr
        write(nodedatafile2,'(A14,A3,A4)') nodedatafile1,Mstr,'.tec'
        print *,"nodedatafile2 ",nodedatafile2
        open(unit=11,file=nodedatafile2)

        write(11,'(A17)',ADVANCE="NO") 'VARIABLES="X","Y"'
        do im=1,nmat
         write(11,'(A3)',ADVANCE="NO") ',"T'
         write(matstr,'(I2)') im
         do i=1,2
          if (matstr(i:i).eq.' ') then
           matstr(i:i)='0'
          endif
         enddo
         write(11,'(A2)',ADVANCE="NO") matstr
         write(11,'(A1)',ADVANCE="NO") '"'
        enddo ! im=1..nmat
        do iten=1,global_nten
         write(11,'(A3)',ADVANCE="NO") ',"U'
         write(matstr,'(I2)') iten
         do i=1,2
          if (matstr(i:i).eq.' ') then
           matstr(i:i)='0'
          endif
         enddo
         write(11,'(A2)',ADVANCE="NO") matstr
         write(11,'(A1)',ADVANCE="NO") '"'

         write(11,'(A3)',ADVANCE="NO") ',"V'
         write(matstr,'(I2)') iten
         do i=1,2
          if (matstr(i:i).eq.' ') then
           matstr(i:i)='0'
          endif
         enddo
         write(11,'(A2)',ADVANCE="NO") matstr
         write(11,'(A1)',ADVANCE="NO") '"'
        enddo ! iten=1..global_nten

        do im=1,nmat
         write(11,'(A3)',ADVANCE="NO") ',"F'
         write(matstr,'(I2)') im
         do i=1,2
          if (matstr(i:i).eq.' ') then
           matstr(i:i)='0'
          endif
         enddo
         write(11,'(A2)',ADVANCE="NO") matstr
         write(11,'(A1)',ADVANCE="NO") '"'
        enddo ! im=1..nmat

        do im=1,nmat
         write(11,'(A4)',ADVANCE="NO") ',"LS'
         write(matstr,'(I2)') im
         do i=1,2
          if (matstr(i:i).eq.' ') then
           matstr(i:i)='0'
          endif
         enddo
         write(11,'(A2)',ADVANCE="NO") matstr
         write(11,'(A1)',ADVANCE="NO") '"'
        enddo ! im=1..nmat

        do im=1,nmat
         write(11,'(A4)',ADVANCE="NO") ',"NX'
         write(matstr,'(I2)') im
         do i=1,2
          if (matstr(i:i).eq.' ') then
           matstr(i:i)='0'
          endif
         enddo
         write(11,'(A2)',ADVANCE="NO") matstr
         write(11,'(A1)',ADVANCE="NO") '"'

         write(11,'(A4)',ADVANCE="NO") ',"NY'
         write(matstr,'(I2)') im
         do i=1,2
          if (matstr(i:i).eq.' ') then
           matstr(i:i)='0'
          endif
         enddo
         write(11,'(A2)',ADVANCE="NO") matstr
         if (im.lt.nmat) then
          write(11,'(A1)',ADVANCE="NO") '"'
         else
          write(11,'(A1)') '"'
         endif
        enddo ! im=1..nmat

        write(11,*) 'zone i=',hix-lox+2,' j=', &
         hiy-loy+2,' f=point'
        write(11,*) 'SOLUTIONTIME=',out_time," STRANDID=",strandid

        do i=lox,hix+1
        do j=loy,hiy+1
         xpoint=i*h
         ypoint=j*h
         lscomp=nmat+global_nten*sdim+ngeom_reconCG*nmat
         write(11,'(D25.16)',ADVANCE="NO") xpoint
         write(11,'(D25.16)',ADVANCE="NO") ypoint
         do im=1,nmat+global_nten*sdim
          sum_wt=0.0d0
          do ii=0,1
          do jj=0,1
           wt_tecplot(ii,jj)=1.0d0
           if ((im.ge.1).and.(im.le.nmat)) then
            wt_tecplot(ii,jj)=VFRAC_MOF(i+ii-1,j+jj-1,im)
           endif
           sum_wt=sum_wt+wt_tecplot(ii,jj)
          enddo
          enddo
          if (sum_wt.eq.0.0d0) then
           sum_wt=0.0d0
           do ii=0,1
           do jj=0,1
            wt_tecplot(ii,jj)=1.0d0
            sum_wt=sum_wt+wt_tecplot(ii,jj)
           enddo
           enddo
          endif
          node_data=0.0d0
          do ii=0,1
          do jj=0,1
           node_data=node_data+ &
             local_UNEW(i+ii-1,j+jj-1,im)*wt_tecplot(ii,jj)
          enddo
          enddo
          node_data=node_data/sum_wt
          write(11,'(D25.16)',ADVANCE="NO") node_data
         enddo
         do im=1,nmat
          sum_wt=0.0d0
          node_data=0.0d0
          do ii=0,1
          do jj=0,1
           wt_tecplot(ii,jj)=1.0d0
           sum_wt=sum_wt+wt_tecplot(ii,jj)
           node_data=node_data+ &
             VFRAC_MOF(i+ii-1,j+jj-1,im)*wt_tecplot(ii,jj)
          enddo
          enddo
          node_data=node_data/sum_wt
          write(11,'(D25.16)',ADVANCE="NO") node_data
         enddo
         do im=1,nmat*(sdim+1)
          sum_wt=0.0d0
          node_data=0.0d0
          do ii=0,1
          do jj=0,1
           wt_tecplot(ii,jj)=1.0d0
           sum_wt=sum_wt+wt_tecplot(ii,jj)
           node_data=node_data+ &
            local_UNEW(i+ii-1,j+jj-1,lscomp+im)*wt_tecplot(ii,jj)
          enddo
          enddo
          node_data=node_data/sum_wt
          if (im.lt.nmat*(sdim+1)) then
           write(11,'(D25.16)',ADVANCE="NO") node_data
          else
           write(11,'(D25.16)') node_data
          endif
         enddo
        enddo
        enddo

        close(11)

        write(cendatastr,'(A7)') 'cendata'
        write(cendatafile1,'(A7,A6)') cendatastr,stepstr
        write(cendatafile2,'(A13,A3,A4)') cendatafile1,Mstr,'.tec'
        print *,"cendatafile2 ",cendatafile2
        open(unit=11,file=cendatafile2)

        write(11,'(A17)',ADVANCE="NO") 'VARIABLES="X","Y"'
        do im=1,nmat
         write(11,'(A3)',ADVANCE="NO") ',"T'
         write(matstr,'(I2)') im
         do i=1,2
          if (matstr(i:i).eq.' ') then
           matstr(i:i)='0'
          endif
         enddo
         write(11,'(A2)',ADVANCE="NO") matstr
         write(11,'(A1)',ADVANCE="NO") '"'
        enddo ! im=1..nmat
        do iten=1,global_nten
         write(11,'(A3)',ADVANCE="NO") ',"U'
         write(matstr,'(I2)') iten
         do i=1,2
          if (matstr(i:i).eq.' ') then
           matstr(i:i)='0'
          endif
         enddo
         write(11,'(A2)',ADVANCE="NO") matstr
         write(11,'(A1)',ADVANCE="NO") '"'

         write(11,'(A3)',ADVANCE="NO") ',"V'
         write(matstr,'(I2)') iten
         do i=1,2
          if (matstr(i:i).eq.' ') then
           matstr(i:i)='0'
          endif
         enddo
         write(11,'(A2)',ADVANCE="NO") matstr
         write(11,'(A1)',ADVANCE="NO") '"'
        enddo ! iten=1..global_nten

        do im=1,nmat
         write(11,'(A3)',ADVANCE="NO") ',"F'
         write(matstr,'(I2)') im
         do i=1,2
          if (matstr(i:i).eq.' ') then
           matstr(i:i)='0'
          endif
         enddo
         write(11,'(A2)',ADVANCE="NO") matstr
         write(11,'(A1)',ADVANCE="NO") '"'
        enddo ! im=1..nmat

        do im=1,nmat
         write(11,'(A4)',ADVANCE="NO") ',"LS'
         write(matstr,'(I2)') im
         do i=1,2
          if (matstr(i:i).eq.' ') then
           matstr(i:i)='0'
          endif
         enddo
         write(11,'(A2)',ADVANCE="NO") matstr
         write(11,'(A1)',ADVANCE="NO") '"'
        enddo ! im=1..nmat

        do im=1,nmat
         write(11,'(A4)',ADVANCE="NO") ',"NX'
         write(matstr,'(I2)') im
         do i=1,2
          if (matstr(i:i).eq.' ') then
           matstr(i:i)='0'
          endif
         enddo
         write(11,'(A2)',ADVANCE="NO") matstr
         write(11,'(A1)',ADVANCE="NO") '"'

         write(11,'(A4)',ADVANCE="NO") ',"NY'
         write(matstr,'(I2)') im
         do i=1,2
          if (matstr(i:i).eq.' ') then
           matstr(i:i)='0'
          endif
         enddo
         write(11,'(A2)',ADVANCE="NO") matstr
         if (im.lt.nmat) then
          write(11,'(A1)',ADVANCE="NO") '"'
         else
          write(11,'(A1)') '"'
         endif
        enddo ! im=1..nmat

        write(11,*) 'zone i=',hix-lox+3,' j=', &
         hiy-loy+3,' f=point'
        write(11,*) 'SOLUTIONTIME=',out_time," STRANDID=",strandid

        do i=lox-1,hix+1
        do j=loy-1,hiy+1
         xpoint=(i+0.5d0)*h
         ypoint=(j+0.5d0)*h
         lscomp=nmat+global_nten*sdim+ngeom_reconCG*nmat
           ! D25.16=" -0.1234..(16)D+123"
         write(11,'(D25.16)',ADVANCE="NO") xpoint
         write(11,'(D25.16)',ADVANCE="NO") ypoint
          ! temperature and velocity
         do im=1,nmat+global_nten*sdim
          cen_data=local_UNEW(i,j,im)
          write(11,'(D25.16)',ADVANCE="NO") cen_data
         enddo
          ! volume fractions
         do im=1,nmat
          cen_data=VFRAC_MOF(i,j,im)
          write(11,'(D25.16)',ADVANCE="NO") cen_data
         enddo
          ! level set data
         do im=1,nmat*(sdim+1)
          cen_data=local_UNEW(i,j,lscomp+im)
          if (im.lt.nmat*(sdim+1)) then
           write(11,'(D25.16)',ADVANCE="NO") cen_data
          else
           write(11,'(D25.16)') cen_data
          endif
         enddo
        enddo
        enddo

        close(11)
     
       endif
 
      else if (plot_int.eq.0) then
       ! do nothing
      else
       print *,"plot_int invalid"
       stop
      endif

      return
      end subroutine output_solution

      end module bicgstab_module

