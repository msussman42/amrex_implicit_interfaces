#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "SPACE.H"

! Author: Mark Sussman sussman@math.fsu.edu
! Department of Mathematics
! Florida State University
! Tallahassee, FL 32306
!

module fully_connected_module

implicit none

type weight_matrix_type
 INTEGER_T :: size_l,size_lm1
 REAL_T, pointer :: weights(:,:) ! 1..size_l,1..size_lm1
 REAL_T, pointer :: bias(:) ! 1..size_l
end type weight_matrix_type

type layer_type
 INTEGER_T :: n_nodes
 REAL_T, pointer :: state_data_predict(:,:) ! i=1..n_nodes,r=1..n_predict
 REAL_T, pointer :: hidden_data_predict(:,:) ! i=1..n_nodes,r=1..n_predict
 REAL_T, pointer :: state_data(:,:)  ! i=1..n_nodes,r=1..n_samples
 REAL_T, pointer :: hidden_data(:,:)  ! i=1..n_nodes,r=1..n_samples
 REAL_T, pointer :: DCDS(:,:)  ! i=1..n_nodes,r=1..n_samples
 type(weight_matrix_type) :: wt_data
 type(weight_matrix_type) :: DCDW
 type(weight_matrix_type) :: moment1
 type(weight_matrix_type) :: moment2
end type layer_type

type NN_type
! parameter for activation function, 
! e.g. for sigmoid function:
! sigma(x)=1/(1+exp(-x))
! use sigma(x/eps_sigmoid) for activation.
 REAL_T eps_activation
 INTEGER_T activation_id ! 0=sigmoid 1=tanh 2=parametric RELU
 INTEGER_T n_input
 INTEGER_T n_output
 INTEGER_T n_middle
 ! l=1,n_layers l=1 input layer l=n_layers output layer  
 INTEGER_T n_layers 
 INTEGER_T n_predict_samples
 INTEGER_T n_total_samples
 INTEGER_T n_batch_samples
 INTEGER_T predict_flag
 INTEGER_T verbose
 INTEGER_T, pointer :: nodes_per_layer(:) ! l=1..n_layers
 type(layer_type), pointer :: layer_array(:) ! l=1..n_layers
 REAL_T, pointer :: Xtrain(:,:)  ! i=1..n_input,r=1..n_total_samples
 REAL_T, pointer :: Ytrain(:,:)  ! i=1..n_output,r=1..n_total_samples
 REAL_T, pointer :: Xpredict(:,:)  ! i=1..n_input,r=1..n_predict_samples
 REAL_T, pointer :: Ypredict(:,:)  ! i=1..n_output,r=1..n_predict_samples
 REAL_T, pointer :: min_Xtrain(:) ! i=1..n_input
 REAL_T, pointer :: max_Xtrain(:) ! i=1..n_input
 REAL_T, pointer :: X_random_basis(:,:) ! i=1..n_input,j=1..nodes_per_layer(2)

 REAL_T, pointer :: A(:,:) !i,j=1..nodes_per_layer(2)+1
 REAL_T, pointer :: B(:) !i=1..nodes_per_layer(2)+1

 REAL_T cost
 REAL_T cost_scale

 REAL_T adam_alpha  ! step size
 REAL_T adam_beta1
 REAL_T adam_beta2
 REAL_T adam_eps
 INTEGER_T adam_timestep
 INTEGER_T steepest_descent_flag
 INTEGER_T full_cycle_min_sweeps
end type NN_type

contains

      REAL_T function radial_basis_fn(x,c,n,eps)
      IMPLICIT NONE

      INTEGER_T n
      REAL_T eps
      REAL_T x(n)
      REAL_T c(n)
      INTEGER_T i
      REAL_T e_power

      if (n.ge.1) then
       e_power=zero
       do i=1,n
        e_power=e_power+(x(i)-c(i))**2
       enddo
       radial_basis_fn=exp(-e_power/eps)
      else
       print *,"n invalid"
       stop
      endif

      return
      end function radial_basis_fn

      REAL_T function NN_test_func(x,probtype)
      IMPLICIT NONE

      REAL_T x
      INTEGER_T probtype

      if (probtype.eq.1) then
       NN_test_func=sin(x)
      else if (probtype.eq.2) then
       NN_test_func=x**3/two-two*(x**2)/three+x/four+two
      else
       print *,"probtype invalid"
       stop
      endif

      return
      end function NN_test_func


      REAL_T function activation_func(x,r,level,node,NN)
      IMPLICIT NONE

      type(NN_type) :: NN
      INTEGER_T i,r
      REAL_T x,y
      INTEGER_T level,node
      REAL_T eps
      INTEGER_T id
      REAL_T, allocatable :: x_data(:)
      REAL_T, allocatable :: x_center(:)
      INTEGER_T n_input

      eps=NN%eps_activation
      id=NN%activation_id
      n_input=NN%n_input

      if ((level.ge.2).and.(n_input.ge.1).and. &
          (r.ge.1).and. &
          (r.le.NN%n_total_samples).and. &
          (node.ge.1).and. &
          (node.le.NN%nodes_per_layer(level))) then

       if (eps.gt.zero) then

        if (id.eq.0) then
         y=x/eps
         activation_func=exp(y)/(exp(y)+one) ! sigmoid
        else if (id.eq.1) then
         y=x/eps
         activation_func=tanh(y)
        else if (id.eq.2) then ! parametric RELU
         if (x.lt.zero) then
          activation_func=eps*x
         else if (x.ge.zero) then
          activation_func=x
         else
          print *,"x bust"
          stop
         endif
        else if (id.eq.3) then ! polynomial least sqr emulator
         if (level.eq.2) then
          allocate(x_data(n_input))
          do i=1,n_input
           if (NN%predict_flag.eq.0) then
            x_data(i)=NN%Xtrain(i,r)
           else if (NN%predict_flag.eq.1) then
            x_data(i)=NN%Xpredict(i,r)
           else
            print *,"predict_flag invalid"
            stop
           endif
          enddo
          if (n_input.eq.1) then
           activation_func=x_data(1)**node
          else
           print *,"n_input out of range"
           stop
          endif
          deallocate(x_data)
         else if (level.eq.3) then
          activation_func=x
         else
          print *,"level invalid"
          stop
         endif
         ! radial basis functions, randomly distributed
        else if (id.eq.4) then 
         if (level.eq.2) then
          allocate(x_data(n_input))
          allocate(x_center(n_input))
          do i=1,n_input
           if (NN%predict_flag.eq.0) then
            x_data(i)=NN%Xtrain(i,r)
           else if (NN%predict_flag.eq.1) then
            x_data(i)=NN%Xpredict(i,r)
           else
            print *,"predict_flag invalid"
            stop
           endif
           x_center(i)=NN%X_random_basis(i,node)
          enddo

          activation_func=radial_basis_fn(x_data,x_center,n_input,eps)

          deallocate(x_center)
          deallocate(x_data)
         else if (level.eq.3) then
          activation_func=x
         else
          print *,"level invalid"
          stop
         endif

        else
         print *,"id invalid"
         stop
        endif
       else
        print *,"eps invalid"
        stop
       endif
      else
       print *,"level, n_input, node or r invalid"
       stop
      endif

      return
      end function activation_func

       ! for sigmoid:
       ! (ex/(1+ex))'=(ex(1+ex)-ex(ex))/(1+ex)^2=
       ! ex/(1+ex)^2=(1-s)s
      REAL_T function activation_prime(x,r,level,node,NN)
      IMPLICIT NONE

      type(NN_type) :: NN
      INTEGER_T i,r
      REAL_T x
      INTEGER_T level,node
      REAL_T eps
      INTEGER_T id
      REAL_T, allocatable :: x_data(:)
      REAL_T, allocatable :: x_center(:)
      INTEGER_T n_input
      REAL_T s

      eps=NN%eps_activation
      id=NN%activation_id
      n_input=NN%n_input

      if ((level.ge.2).and.(n_input.ge.1).and. &
          (r.ge.1).and. &
          (r.le.NN%n_total_samples).and. &
          (node.ge.1).and. &
          (node.le.NN%nodes_per_layer(level))) then

       if (eps.gt.zero) then

        if (id.eq.0) then ! sigmoid
         s=activation_func(x,r,level,node,NN)
         activation_prime=s*(one-s)/eps
        else if (id.eq.1) then ! tanh
         s=activation_func(x,r,level,node,NN)
         activation_prime=(one-s*s)/eps
        else if (id.eq.2) then ! parametric RELU
         if (x.lt.zero) then
          activation_prime=eps
         else if (x.ge.zero) then
          activation_prime=one
         else
          print *,"x invalid"
          stop
         endif
        else if (id.eq.3) then ! polynomial least sqr emulator
         if (level.eq.2) then
          allocate(x_data(n_input))
          do i=1,n_input
           if (NN%predict_flag.eq.0) then
            x_data(i)=NN%Xtrain(i,r)
           else if (NN%predict_flag.eq.1) then
            x_data(i)=NN%Xpredict(i,r)
           else
            print *,"predict_flag invalid"
            stop
           endif
          enddo
          if (n_input.eq.1) then
           activation_prime=node*(x_data(1)**(node-1))
          else
           print *,"n_input out of range"
           stop
          endif
          deallocate(x_data)
         else if (level.eq.3) then
          activation_prime=one
         else
          print *,"level invalid"
          stop
         endif

         ! radial basis functions, randomly distributed
        else if (id.eq.4) then 
         if (level.eq.2) then
          allocate(x_data(n_input))
          allocate(x_center(n_input))
          do i=1,n_input
           if (NN%predict_flag.eq.0) then
            x_data(i)=NN%Xtrain(i,r)
           else if (NN%predict_flag.eq.1) then
            x_data(i)=NN%Xpredict(i,r)
           else
            print *,"predict_flag invalid"
            stop
           endif
           x_center(i)=NN%X_random_basis(i,node)
          enddo

          activation_prime=zero

          deallocate(x_center)
          deallocate(x_data)
         else if (level.eq.3) then
          activation_prime=one
         else
          print *,"level invalid"
          stop
         endif

        else
         print *,"id invalid"
         stop
        endif
       else
        print *,"eps invalid"
        stop
       endif
      else
       print *,"level, n_input, node or r invalid"
       stop
      endif

      return
      end function activation_prime

      subroutine init_RBF_centroids(NN)
      IMPLICIT NONE

      type(NN_type) :: NN
      real (kind=8), allocatable :: r_normal_vec(:)
      integer ( kind = 4 ) seed
      INTEGER_T i,j,counter
      REAL_T xrange

      seed=1
      ALLOCATE(r_normal_vec(NN%n_input*NN%nodes_per_layer(2)))
      call r8vec_uniform_01(NN%n_input*NN%nodes_per_layer(2),seed,r_normal_vec)
      counter=1
      do i=1,NN%n_input
       xrange=NN%max_Xtrain(i)-NN%min_Xtrain(i)
       if (xrange.gt.zero) then
        do j=1,NN%nodes_per_layer(2)
         NN%X_random_basis(i,j)=NN%min_Xtrain(i)+ &
          r_normal_vec(counter)*xrange
         counter=counter+1
        enddo
       else
        print *,"xrange invalid"
        stop
       endif
      enddo  ! i=1..n_input
      DEALLOCATE(r_normal_vec)

      return
      end subroutine init_RBF_centroids

       ! after clusters are found:
       ! y=sum wi phi_i(x)
       ! C=sum_j ||y_j-sum_i wi phi_i(x_j)||^2
       ! dC/dwk=0: sum_j y_j phi_k(x_j)=sum_i wi sum_j phi_i(x_j)phi_k(x_j)
       ! Aki=sum_j phi_i(x_j)phi_k(x_j)=sum_j phi_j phi_j^T
       ! z^T A z=sum_j z^T phi_j phi_j^T z=sum_j (phi_j^T z)^T (phi_j^T z)>0
       ! so A is SPD => we can use CG or PCG
       ! (but to be most general, we will use PBiCGSTAB)
       ! P. Yin, M. Pham, A. Oberman, S. Osher, JSC 2018
      subroutine stochastic_BE_for_K_means(NN)
      IMPLICIT NONE

       ! Xtrain(:,:) i=1..n_input, r=1..n_total_samples
       ! n_total_samples is the mini-batch size
       ! X_random_basis(:,:) i=1..n_input,j=1..nodes_per_layer(2) are 
       !  the centroids which are found.
      type(NN_type) :: NN
      INTEGER_T d,K,M,imaxit,omaxit
      REAL_T gammaBE,alpha,beta
      REAL_T, allocatable :: X_random_basis_new(:,:)
      REAL_T, allocatable :: Y_random_basis(:,:)
      REAL_T, allocatable :: gradphi(:,:)
      INTEGER_T, allocatable :: Xtrain_map(:)
      INTEGER_T outer_iter,dir,j,l,r
      INTEGER_T j_crit
      REAL_T dist,dist_crit,phi,mag_gradphi
      REAL_T min_phi,min_mag_gradphi
      character*9 filename9
      character*8 filename8

       ! randomly populate NN%X_random_basis
      call init_RBF_centroids(NN)
      K=NN%nodes_per_layer(2)  ! number of clusters
      d=NN%n_input  ! dimension of the centroids
      gammaBE=K
      M=NN%n_total_samples ! mini-batch size
      alpha=0.75 ! averaging parameter
      beta=one/1.01  ! decay parameter
      imaxit=10
      omaxit=150

      ALLOCATE(X_random_basis_new(d,K))
      ALLOCATE(Y_random_basis(d,K))
      ALLOCATE(gradphi(d,K))
      ALLOCATE(Xtrain_map(M))

      do outer_iter=1,omaxit

       do dir=1,d
       do j=1,K
        Y_random_basis(dir,j)=NN%X_random_basis(dir,j)
        X_random_basis_new(dir,j)=Y_random_basis(dir,j)
       enddo
       enddo

       do l=1,imaxit

        ! find gradphi(Y_random_basis) 
        ! first we associate to each Xtrain point a corresponding closest 
        ! Y_random_basis point.
        do r=1,M
         dist_crit=zero
         j_crit=0
         do j=1,K
          dist=zero
          do dir=1,d
           dist=dist+(NN%Xtrain(dir,r)-Y_random_basis(dir,j))**2
          enddo
          dist=sqrt(dist)
          if (j_crit.eq.0) then
           j_crit=j
           dist_crit=dist
          else if ((j_crit.ge.1).and.(j_crit.le.K).and. &
                   (dist_crit.ge.zero)) then
           if (dist.lt.dist_crit) then
            j_crit=j
            dist_crit=dist
           endif
          else
           print *,"j_crit or dist_crit invalid"
           stop
          endif
         enddo ! j=1..K
         if ((j_crit.ge.1).and.(j_crit.le.K).and.(dist_crit.ge.zero)) then
          Xtrain_map(r)=j_crit
         else
          print *,"j_crit or dist_crit invalid"
          stop
         endif 
        enddo !r=1..M

        do dir=1,d
        do j=1,K
         gradphi(dir,j)=zero
        enddo
        enddo
        phi=zero
        do dir=1,d
        do r=1,M
         j_crit=Xtrain_map(r)
         gradphi(dir,j_crit)=gradphi(dir,j_crit)+ &
            (Y_random_basis(dir,j_crit)-NN%Xtrain(dir,r))
         phi=phi+ &
          (Y_random_basis(dir,j_crit)-NN%Xtrain(dir,r))**2
        enddo    
        enddo 
        phi=phi/(two*M)
        mag_gradphi=zero   
        do dir=1,d
        do j=1,K
         gradphi(dir,j)=gradphi(dir,j)/M
         mag_gradphi=mag_gradphi+gradphi(dir,j)**2
        enddo
        enddo
        mag_gradphi=sqrt(mag_gradphi)
        do dir=1,d
        do j=1,K
         Y_random_basis(dir,j)=NN%X_random_basis(dir,j)- &
            gammaBE*gradphi(dir,j)
         X_random_basis_new(dir,j)=alpha*X_random_basis_new(dir,j)+ &
          (one-alpha)*Y_random_basis(dir,j)
        enddo
        enddo

        print *,"outer_iter,l,mag_gradphi,phi ",outer_iter,l,mag_gradphi,phi

        if ((l.eq.1).and.(outer_iter.eq.1)) then
         min_phi=phi
         min_mag_gradphi=mag_gradphi
        else
         if (phi.lt.min_phi) then
          min_phi=phi
         endif
         if (mag_gradphi.lt.min_mag_gradphi) then
          min_mag_gradphi=mag_gradphi
         endif
        endif 

       enddo ! l=1..imaxit

       gammaBE=gammaBE*beta

       do dir=1,d
       do j=1,K
        NN%X_random_basis(dir,j)=X_random_basis_new(dir,j)
       enddo
       enddo

      enddo ! outer_iter=1..omaxit
      print *,"min_mag_gradphi,min_phi ",min_mag_gradphi,min_phi

      print *,"K means clustering points in the file: NNCLUSTER"
      write(filename9,'(A9)') 'NNCLUSTER'
      open(unit=11,file=filename9)
      do j=1,K
       do dir=1,d
        if (dir.lt.d) then
         write(11,'(G25.16e3)',ADVANCE="NO") NN%X_random_basis(dir,j)
        else if (dir.eq.d) then
         write(11,'(G25.16e3)',ADVANCE="NO") NN%X_random_basis(dir,j)
         write(11,'(G25.16e3)') zero
        endif
       enddo
      enddo
      close(11)

      print *,"original Xtrain points in the file: NNXTRAIN"
      write(filename8,'(A8)') 'NNXTRAIN'
      open(unit=11,file=filename8)
      do j=1,M
       do dir=1,d
        if (dir.lt.d) then
         write(11,'(G25.16e3)',ADVANCE="NO") NN%Xtrain(dir,j)
        else if (dir.eq.d) then
         write(11,'(G25.16e3)',ADVANCE="NO") NN%Xtrain(dir,j)
         write(11,'(G25.16e3)') zero
        endif
       enddo
      enddo
      close(11)

      DEALLOCATE(Xtrain_map)
      DEALLOCATE(gradphi)
      DEALLOCATE(Y_random_basis)
      DEALLOCATE(X_random_basis_new)

      end subroutine stochastic_BE_for_K_means

      subroutine ATIMESU(AU,U,NN,K)

      IMPLICIT NONE

      type(NN_type) :: NN
      INTEGER_T K
      REAL_T AU(K+1)
      REAL_T U(K+1)
      INTEGER_T i,j

      if (K.ne.NN%nodes_per_layer(2)) then
       print *,"K.ne.NN%nodes_per_layer(2)"
       stop
      endif

      do i=1,K+1
       AU(i)=zero
       do j=1,K+1
        AU(i)=AU(i)+NN%A(i,j)*U(j)
       enddo
      enddo
      
      return
      end subroutine ATIMESU


      subroutine LINCOMB(U,V,W,AA,BB,NN,K)

      IMPLICIT NONE

      type(NN_type) :: NN
      INTEGER_T K
      REAL_T U(K+1)
      REAL_T V(K+1)
      REAL_T W(K+1)
      REAL_T AA,BB
      INTEGER_T i

      if (K.ne.NN%nodes_per_layer(2)) then
       print *,"K.ne.NN%nodes_per_layer(2)"
       stop
      endif

      do i=1,K+1
       W(i)=AA*U(i)+BB*V(i)
      enddo
      
      return
      end subroutine LINCOMB


      subroutine DOTPROD(U,V,dsum,NN,K)

      IMPLICIT NONE

      type(NN_type) :: NN
      INTEGER_T K
      REAL_T U(K+1)
      REAL_T V(K+1)
      REAL_T dsum
      INTEGER_T i

      if (K.ne.NN%nodes_per_layer(2)) then
       print *,"K.ne.NN%nodes_per_layer(2)"
       stop
      endif

      dsum=zero
      do i=1,K+1
       dsum=dsum+U(i)*V(i)
      enddo
      
      return
      end subroutine DOTPROD


      subroutine NORMPROD(U,dnorm,NN,K)

      IMPLICIT NONE

      type(NN_type) :: NN
      INTEGER_T K
      REAL_T U(K+1)
      REAL_T dnorm,dsum

      if (K.ne.NN%nodes_per_layer(2)) then
       print *,"K.ne.NN%nodes_per_layer(2)"
       stop
      endif
 
      call DOTPROD(U,U,dsum,NN,K)
      dnorm=sqrt(dsum)
      
      return
      end subroutine NORMPROD

        ! R=RHS-A U
      subroutine RESID(R,RHS,U,NN,K)
      IMPLICIT NONE

      type(NN_type) :: NN
      INTEGER_T K
      REAL_T R(K+1)
      REAL_T RHS(K+1)
      REAL_T U(K+1)
      REAL_T, allocatable :: AU(:)
      REAL_T AA,BB

      if (K.ne.NN%nodes_per_layer(2)) then
       print *,"K.ne.NN%nodes_per_layer(2)"
       stop
      endif

      allocate(AU(K+1))

      call ATIMESU(AU,U,NN,K)
      AA=one
      BB=-one
      call LINCOMB(RHS,AU,R,AA,BB,NN,K)

      deallocate(AU)

      return
      end subroutine RESID

      subroutine ZAPVEC(Z,NN,K)
      IMPLICIT NONE

      type(NN_type) :: NN
      INTEGER_T K
      REAL_T Z(K+1)
      INTEGER_T i

      if (K.ne.NN%nodes_per_layer(2)) then
       print *,"K.ne.NN%nodes_per_layer(2)"
       stop
      endif

      do i=1,K+1
       Z(i)=zero
      enddo

      return
      end subroutine ZAPVEC


      subroutine COPYVEC(U,V,NN,K)
      IMPLICIT NONE

      type(NN_type) :: NN
      INTEGER_T K
      REAL_T U(K+1)
      REAL_T V(K+1)
      INTEGER_T i

      if (K.ne.NN%nodes_per_layer(2)) then
       print *,"K.ne.NN%nodes_per_layer(2)"
       stop
      endif

      do i=1,K+1
       V(i)=U(i)
      enddo

      return
      end subroutine COPYVEC


       ! for k=1,2,3 ....
       ! Z^{k}=Z^{k-1}+D^{-1}(R-AZ^{k-1})
       ! e.g. A=D-L-U
       ! D Z^k = D Z^k-1 + R - (D-L-U)Z^k-1
       ! D Z^k - (L+U)Z^k-1=R 
      subroutine JACPRECOND(Z,R,NN,K)
      IMPLICIT NONE

      type(NN_type) :: NN
      INTEGER_T K
      REAL_T Z(K+1)
      REAL_T R(K+1)

      REAL_T, allocatable :: ZN(:)
      REAL_T, allocatable :: ZNP1(:)
      REAL_T, allocatable :: RSTAR(:)
      INTEGER_T iter,nsmooth,i

      if (K.ne.NN%nodes_per_layer(2)) then
       print *,"K.ne.NN%nodes_per_layer(2)"
       stop
      endif
      allocate(ZN(K+1))
      allocate(ZNP1(K+1))
      allocate(RSTAR(K+1))

      nsmooth=8

       ! ZN=Z
      call COPYVEC(Z,ZN,NN,K)

      do iter=1,nsmooth
       call RESID(RSTAR,R,ZN,NN,K)

       do i=1,K+1
        if (NN%A(i,i).gt.zero) then
         ZNP1(i)=ZN(i)+(one/NN%A(i,i))*RSTAR(i)
        else
         print *,"NN%A(i,i) invalid"
         stop
        endif
       enddo
        ! ZN=ZNP1
       call COPYVEC(ZNP1,ZN,NN,K)
      enddo  ! iter
        ! Z=ZN
      call COPYVEC(ZN,Z,NN,K)

      deallocate(ZN) 
      deallocate(ZNP1) 
      deallocate(RSTAR) 

      return
      end subroutine JACPRECOND


      subroutine preconditioner(Z,R,NN,K)
      IMPLICIT NONE

      type(NN_type) :: NN
      INTEGER_T K
      REAL_T Z(K+1)
      REAL_T R(K+1)
      INTEGER_T precond_type

      precond_type=1

      call ZAPVEC(Z,NN,K)
      if (precond_type.eq.0) then
       call COPYVEC(R,Z,NN,K)
      else if (precond_type.eq.1) then
       call JACPRECOND(Z,R,NN,K)
      else
       print *,"precond_type invalid" 
       stop
      endif

      return
      end subroutine preconditioner


        ! precond_type=0 M=I, =1 Jacobi
      subroutine BICGSTAB(cost_tol,grad_tol,max_steps,NN)
      IMPLICIT NONE

      type(NN_type) :: NN

      REAL_T cost_tol,grad_tol
      INTEGER_T max_steps
      INTEGER_T K

      REAL_T, allocatable :: G(:)
      REAL_T, allocatable :: U(:)
      REAL_T, allocatable :: U0(:)
      REAL_T, allocatable :: V0(:)
      REAL_T, allocatable :: P0(:)
      REAL_T, allocatable :: R0(:)
      REAL_T, allocatable :: U1(:)
      REAL_T, allocatable :: V1(:)
      REAL_T, allocatable :: P1(:)
      REAL_T, allocatable :: R1(:)
      REAL_T, allocatable :: R0hat(:)
      REAL_T, allocatable :: UINIT(:)
      REAL_T, allocatable :: RHS(:)
      REAL_T, allocatable :: Y(:)
      REAL_T, allocatable :: Hvec(:)
      REAL_T, allocatable :: S(:)
      REAL_T, allocatable :: T(:)
      REAL_T, allocatable :: Z(:)

      REAL_T rho0,w0,rho1,AA,w1,BB,a1,a2
      REAL_T dnorm,dnorm0
      INTEGER_T i
      INTEGER_T iter,restart_flag

      K=NN%nodes_per_layer(2)

      allocate(G(K+1))
      allocate(U(K+1))
      allocate(U0(K+1))
      allocate(V0(K+1))
      allocate(P0(K+1))
      allocate(R0(K+1))
      allocate(U1(K+1))
      allocate(V1(K+1))
      allocate(P1(K+1))
      allocate(R1(K+1))
      allocate(R0hat(K+1))
      allocate(UINIT(K+1))
      allocate(RHS(K+1))
      allocate(Y(K+1))
      allocate(Hvec(K+1))
      allocate(S(K+1))
      allocate(T(K+1))
      allocate(Z(K+1))

        ! U0=V0=P0=0 
      AA=zero
      do i=1,K+1
       U0(i)=AA
       V0(i)=AA
       P0(i)=AA
       G(i)=NN%B(i)
      enddo

       ! R0=G-A U0
      call RESID(R0,G,U0,NN,K)

       ! R0hat=R0
      call COPYVEC(R0,R0hat,NN,K)
       ! UINIT=U0
      call COPYVEC(U0,UINIT,NN,K)
      call ZAPVEC(U0,NN,K)
       ! RHS=R0
      call COPYVEC(R0,RHS,NN,K)

       ! rho0=AA=w0=1
      rho0=one
      AA=one
      w0=one

      call NORMPROD(R0,dnorm0,NN,K)
      print *,"initial,dnorm0 ",dnorm0
      dnorm=one
      iter=0

      do while ((dnorm.gt.grad_tol).and.(iter.lt.max_steps))
       print *,"iter,dnorm ",iter,dnorm

         ! rho1= R0hat^H R0
       call DOTPROD(R0hat,R0,rho1,NN,K)
       
       restart_flag=0
       if ((sqrt(abs(rho0)).lt.grad_tol*0.01).or. &
           (sqrt(abs(w0)).lt.grad_tol*0.01)) then
        restart_flag=1
       endif 

       if (restart_flag.eq.0) then
          ! (R0hat^H R0)/(R0hat^H dot R0_before)  *   (AA/w0)
        BB=(rho1/rho0)*(AA/w0)
        a1=one
        a2=-w0
 
         ! P1=P0-w0 V0
        call LINCOMB(P0,V0,P1,a1,a2,NN,K)
         ! P1=R0+BB P1
        call LINCOMB(R0,P1,P1,a1,BB,NN,K)
        ! Y=M^{-1}P1
        call preconditioner(Y,P1,NN,K)
         
         ! V1=A Y
        call ATIMESU(V1,Y,NN,K)

         ! AA=rho1/R0hat dot V1
        call DOTPROD(R0hat,V1,AA,NN,K)

        if (sqrt(abs(AA)).lt.grad_tol*0.01) then
         restart_flag=1
        endif

        if (restart_flag.eq.0) then
         AA=rho1/AA

         ! Hvec=U0+AA Y
         a1=one
         a2=AA
         call LINCOMB(U0,Y,Hvec,a1,a2,NN,K) 
         ! U1=Hvec
         call COPYVEC(Hvec,U1,NN,K)
         ! R1=RHS-A U1
         call RESID(R1,RHS,U1,NN,K)
         call NORMPROD(R1,dnorm,NN,K)
         dnorm=dnorm/dnorm0

         if (dnorm.gt.grad_tol) then
          ! S=R0-AA V1
          a1=one
          a2=-AA
          call LINCOMB(R0,V1,S,a1,a2,NN,K) 

           ! Z=M^{-1}S
          call preconditioner(Z,S,NN,K)

           ! T=A Z
          call ATIMESU(T,Z,NN,K)

           ! simple case is: (T,S)/(T,T)=(AZ,S)/(AZ,AZ)   (MZ=S)
          call DOTPROD(T,S,a1,NN,K)
          call DOTPROD(T,T,a2,NN,K)
          if (sqrt(abs(a2)).lt.grad_tol*0.01) then
           restart_flag=1
          endif

          if (restart_flag.eq.0) then
           w1=a1/a2
           ! U1=Hvec+w1 Z
           a1=one
           a2=w1
           call LINCOMB(Hvec,Z,U1,a1,a2,NN,K) 
          endif
         endif ! dnorm>bicgstab_tol
          ! R1=RHS-A U1
         call RESID(R1,RHS,U1,NN,K)
         call NORMPROD(R1,dnorm,NN,K)
         dnorm=dnorm/dnorm0
         rho0=rho1
         w0=w1
          ! R0=R1
         call COPYVEC(R1,R0,NN,K) 
         call COPYVEC(P1,P0,NN,K) 
         call COPYVEC(V1,V0,NN,K) 
         call COPYVEC(U1,U0,NN,K) 
        endif  ! restart_flag=0
       endif ! restart_flag=0

       if (restart_flag.eq.0) then
        ! do nothing
       else if (restart_flag.eq.1) then
        call RESID(R0,RHS,U0,NN,K) ! R0=RHS-A U0
         ! R0hat=R0
        call COPYVEC(R0,R0hat,NN,K)
        call COPYVEC(U0,U1,NN,K) 
         ! rho0=AA=w0=1
        rho0=one
        AA=one
        w0=one
        call NORMPROD(R0,dnorm,NN,K)
        dnorm=dnorm/dnorm0
        call ZAPVEC(V0,NN,K)
        call ZAPVEC(P0,NN,K)
       else
        print *,"restart_flag invalid"
        stop
       endif

       iter=iter+1
      enddo
      print *,"at the end: iter,dnorm ",iter,dnorm
       ! U=UINIT+U1
      a1=one
      a2=a1
      call LINCOMB(UINIT,U1,U,a1,a2,NN,K)
      do i=1,K+1
       if (i.eq.K+1) then
        NN%layer_array(3)%wt_data%bias(1)=U(i)
       else
        NN%layer_array(3)%wt_data%weights(1,i)=U(i)
       endif
      enddo

      deallocate(G) 
      deallocate(U) 
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
      end subroutine BICGSTAB


       ! Aij=sum_k phi_i(x_k)phi_j(x_k)
       ! Bi=sum_k y_k phi_i(x_k)
      subroutine BICGSTAB_init(NN)
      IMPLICIT NONE

      type(NN_type) :: NN

      REAL_T eps,phi_i,phi_j
      INTEGER_T d,M,K,dir,i,j,l
      REAL_T, allocatable :: x_data(:)
      REAL_T, allocatable :: x_center(:)

      if (NN%n_output.eq.1) then
       if (NN%n_layers.eq.3) then
        if (NN%nodes_per_layer(3).ne.1) then
         print *,"NN%nodes_per_layer(3).ne.1"
         stop
        endif
        eps=NN%eps_activation
        K=NN%nodes_per_layer(2)
        M=NN%n_total_samples
        d=NN%n_input
        allocate(x_data(d))
        allocate(x_center(d))
        ALLOCATE(NN%A(K+1,K+1)) 
        ALLOCATE(NN%B(K+1)) 
        do i=1,K+1
         do j=1,K+1
          NN%A(i,j)=zero
          do l=1,M
           if (i.eq.K+1) then
            phi_i=one
           else
            do dir=1,d
             x_data(dir)=NN%Xtrain(dir,l)
             x_center(dir)=NN%X_random_basis(dir,i)
            enddo
            phi_i=radial_basis_fn(x_data,x_center,d,eps)
           endif
           if (j.eq.K+1) then
            phi_j=one
           else
            do dir=1,d
             x_data(dir)=NN%Xtrain(dir,l)
             x_center(dir)=NN%X_random_basis(dir,j)
            enddo
            phi_j=radial_basis_fn(x_data,x_center,d,eps)
           endif
           NN%A(i,j)=NN%A(i,j)+phi_i*phi_j
          enddo ! l=1..M
         enddo ! j=1..K+1
         NN%B(i)=zero
         do l=1,M
          if (i.eq.K+1) then
           phi_i=one
          else
           do dir=1,d
            x_data(dir)=NN%Xtrain(dir,l)
            x_center(dir)=NN%X_random_basis(dir,i)
           enddo
           phi_i=radial_basis_fn(x_data,x_center,d,eps)
          endif
          NN%B(i)=NN%B(i)+phi_i*NN%Ytrain(1,l)
         enddo ! l=1..M

        enddo !i=1..K+1
           
        deallocate(x_center)
        deallocate(x_data)
       else 
        print *,"NN%n_layers invalid"
        stop
       endif       
      else
       print *,"if dimension(y)>1, then call RBF dimension(y) times"
       stop
      endif

      return
      end subroutine BICGSTAB_init

       ! ADAM: A method for stochastic optimization 
       ! authors: Diederik Kingma and Jimmy Lei Ba
      subroutine ADAM_step(sample_start,sample_end,mag_change,NN)
      IMPLICIT NONE

      INTEGER_T sample_start,sample_end
      type(NN_type) :: NN
      REAL_T mag_change
      INTEGER_T i,j,l
      REAL_T moment1_hat,moment2_hat,delta_theta

      NN%adam_timestep=NN%adam_timestep+1
      NN%predict_flag=0
      call forward_propagate(sample_start,sample_end,NN)
      if (NN%verbose.ge.1) then
       print *,"after forward_prop: timestep,cost,rel cost ", &
         NN%adam_timestep,NN%cost,sqrt(NN%cost/NN%cost_scale)
      endif
      call back_propagate(sample_start,sample_end,NN)
      mag_change=zero
      do l=2,NN%n_layers

       do j=1,NN%nodes_per_layer(l) 

        do i=1,NN%nodes_per_layer(l-1)
         if (NN%steepest_descent_flag.eq.1) then
          delta_theta=NN%layer_array(l)%DCDW%weights(j,i)
         else if (NN%steepest_descent_flag.eq.0) then
          NN%layer_array(l)%moment1%weights(j,i)= &  
           NN%adam_beta1*NN%layer_array(l)%moment1%weights(j,i)+ &  
           (one-NN%adam_beta1)*NN%layer_array(l)%DCDW%weights(j,i)
          NN%layer_array(l)%moment2%weights(j,i)= &  
           NN%adam_beta2*NN%layer_array(l)%moment2%weights(j,i)+ &  
           (one-NN%adam_beta2)*(NN%layer_array(l)%DCDW%weights(j,i)**2)
          moment1_hat= &
           NN%layer_array(l)%moment1%weights(j,i)/ &
           (one-NN%adam_beta1**NN%adam_timestep)
          moment2_hat= &
           NN%layer_array(l)%moment2%weights(j,i)/ &
           (one-NN%adam_beta2**NN%adam_timestep)
          delta_theta=moment1_hat/(sqrt(moment2_hat)+NN%adam_eps)
         else
          print *,"NN%steepest_descent_flag invalid"
          stop
         endif

         NN%layer_array(l)%wt_data%weights(j,i)= &  
           NN%layer_array(l)%wt_data%weights(j,i)- &
           NN%adam_alpha*delta_theta

         mag_change=mag_change+delta_theta**2
        enddo ! i=1..size_lm1

        if (NN%steepest_descent_flag.eq.1) then
         delta_theta=NN%layer_array(l)%DCDW%bias(j)
        else if (NN%steepest_descent_flag.eq.0) then
         NN%layer_array(l)%moment1%bias(j)= &  
           NN%adam_beta1*NN%layer_array(l)%moment1%bias(j)+ &  
           (one-NN%adam_beta1)*NN%layer_array(l)%DCDW%bias(j)
         NN%layer_array(l)%moment2%bias(j)= &  
           NN%adam_beta2*NN%layer_array(l)%moment2%bias(j)+ &  
           (one-NN%adam_beta2)*(NN%layer_array(l)%DCDW%bias(j)**2)
         moment1_hat= &
           NN%layer_array(l)%moment1%bias(j)/ &
           (one-NN%adam_beta1**NN%adam_timestep)
         moment2_hat= &
           NN%layer_array(l)%moment2%bias(j)/ &
           (one-NN%adam_beta2**NN%adam_timestep)
         delta_theta=moment1_hat/(sqrt(moment2_hat)+NN%adam_eps)
        else
         print *,"NN%steepest_descent_flag invalid"
         stop
        endif

        NN%layer_array(l)%wt_data%bias(j)= &  
           NN%layer_array(l)%wt_data%bias(j)- &
           NN%adam_alpha*delta_theta
        mag_change=mag_change+delta_theta**2
       enddo ! j=1..size_l

      enddo ! l=2..n_layers
      mag_change=sqrt(mag_change)

      return
      end subroutine ADAM_step

       ! stop when min number of sweeps done and 
       ! ||ytrain-yNN||/||ytrain|| < cost_tol or
       ! ||grad w||/||grad w init|| <grad_tol
      subroutine ADAM_loop(cost_tol,grad_tol,max_steps,NN)
      IMPLICIT NONE

      REAL_T cost_tol
      REAL_T grad_tol
      type(NN_type) :: NN
      INTEGER_T tolerance_met
      INTEGER_T sample_start,sample_end
      REAL_T mag_change
      REAL_T mag_change_init
      character*5 step_str
      character*11 filename11
      INTEGER_T l,i,j,r,rstate,rwrap
      INTEGER_T counter
      INTEGER_T size_l,size_lm1
      REAL_T rel_cost
      REAL_T rel_grad
      INTEGER_T max_steps
      INTEGER_T first_iter

      sample_start=1
      sample_end=NN%n_batch_samples

      tolerance_met=0
      first_iter=1

      if (NN%activation_id.eq.4) then

       call BICGSTAB(cost_tol,grad_tol,max_steps,NN)

      else

       do while (tolerance_met.eq.0)

        if (NN%verbose.ge.2) then
         write(step_str,'(I5)') NN%adam_timestep
         do i=1,5
          if (step_str(i:i).eq.' ') then
           step_str(i:i)='0'
          endif
         enddo
         write(filename11,'(A6,A5)') 'WTDATA',step_str
         open(unit=11,file=filename11)
         counter=0
         do l=2,NN%n_layers
          size_l=NN%layer_array(l)%wt_data%size_l
          size_lm1=NN%layer_array(l)%wt_data%size_lm1
          if ((size_l.eq.NN%nodes_per_layer(l)).and. &
              (size_lm1.eq.NN%nodes_per_layer(l-1))) then
           do i=1,size_l
            do j=1,size_lm1
             counter=counter+1
             write(11,*)counter, &
               NN%layer_array(l)%wt_data%weights(i,j)
            enddo ! j
           enddo ! i
           do i=1,size_l
            counter=counter+1
            write(11,*)counter, &
               NN%layer_array(l)%wt_data%bias(i)
           enddo ! i
          else
           print *,"size_l or size_lm1 invalid"
           stop
          endif
         enddo ! l
         close(11)
        endif

        call ADAM_step(sample_start,sample_end,mag_change,NN)

        if ((NN%cost.ge.zero).and.(NN%cost_scale.gt.zero)) then
         rel_cost=sqrt(NN%cost/NN%cost_scale)
        else
         print *,"cost or cost_scale invalid"
         stop
        endif
        if (first_iter.eq.1) then
         mag_change_init=mag_change
        else if (first_iter.eq.0) then
         ! do nothing
        else
         print *,"first_iter invalid"
         stop
        endif
        first_iter=0

        if (mag_change.eq.zero) then
         rel_grad=zero
        else if (mag_change.gt.zero) then
         if (mag_change_init.gt.zero) then
          rel_grad=mag_change/mag_change_init
         else
          print *,"mag_change_init invalid"
          stop
         endif
        else
         print *,"mag_change invalid"
         stop
        endif

        if ((NN%adam_timestep/10)*10.eq.NN%adam_timestep) then
         print *,"timestep,cost,rel cost,rel grad ",NN%adam_timestep, &
           NN%cost,rel_cost,rel_grad
        endif

        if (NN%verbose.ge.2) then
         write(step_str,'(I5)') NN%adam_timestep
         do i=1,5
          if (step_str(i:i).eq.' ') then
           step_str(i:i)='0'
          endif
         enddo
         write(filename11,'(A6,A5)') 'NNDATA',step_str
         open(unit=11,file=filename11)
         l=NN%n_layers
         do i=1,NN%n_output
          do r=sample_start,sample_end
           rstate=r-sample_start+1
           rwrap=r
           if (rwrap.gt.NN%n_total_samples) then
            rwrap=rwrap-NN%n_total_samples
           endif
           write(11,*)NN%Xtrain(i,rwrap), &
               NN%layer_array(l)%state_data(i,rstate), &
               NN%Ytrain(i,rwrap)
          enddo ! r
         enddo ! i
         close(11)
        endif

        if (mag_change.ge.zero) then
         sample_start=sample_start+NN%n_batch_samples
         if (sample_start.gt.NN%n_total_samples) then
          sample_start=sample_start-NN%n_total_samples
         endif
         sample_end=sample_start+NN%n_batch_samples-1
       
         if (NN%n_batch_samples*NN%adam_timestep.gt. &
             NN%n_total_samples*max_steps) then
          tolerance_met=1
         endif
         if (NN%n_batch_samples*NN%adam_timestep.gt. &
             NN%n_total_samples*NN%full_cycle_min_sweeps) then
          if (mag_change.eq.zero) then
           tolerance_met=1
          else if (mag_change_init.gt.zero) then
           if (rel_grad.le.grad_tol) then
            tolerance_met=1
           endif
          else
           print *,"mag_change_init invalid"
           stop
          endif
          if (rel_cost.le.cost_tol) then
           tolerance_met=1
          endif
         endif 
        else
         print *,"mag_change invalid"
         stop
        endif

       enddo  ! while (tolerance_met==0)

      endif
     
      return
      end subroutine ADAM_loop
 
       ! Ypredict_parm is both an input (for testing accuracy)
       ! and an output (the predicted value)
      subroutine NN_predict(n_predict_samples_in, &
       Xpredict_in,Ypredict_parm,NN)
      IMPLICIT NONE

      type(NN_type) :: NN
      INTEGER_T n_predict_samples_in
      REAL_T Xpredict_in(NN%n_input,n_predict_samples_in)
      REAL_T Ypredict_parm(NN%n_output,n_predict_samples_in)
      INTEGER_T i,r
      INTEGER_T sample_start,sample_end
      
      NN%n_predict_samples=n_predict_samples_in
      ALLOCATE(NN%Xpredict(NN%n_input,NN%n_predict_samples))
      ALLOCATE(NN%Ypredict(NN%n_output,NN%n_predict_samples))
      do i=1,NN%n_input
      do r=1,NN%n_predict_samples
       NN%Xpredict(i,r)=Xpredict_in(i,r)
      enddo
      enddo
      do i=1,NN%n_output
      do r=1,NN%n_predict_samples
       NN%Ypredict(i,r)=Ypredict_parm(i,r)
      enddo
      enddo

      NN%predict_flag=1
      sample_start=1
      sample_end=NN%n_predict_samples
      call forward_propagate(sample_start,sample_end,NN)
      do i=1,NN%n_output
      do r=1,NN%n_predict_samples
       Ypredict_parm(i,r)=NN%Ypredict(i,r)
      enddo
      enddo
      NN%predict_flag=0

      deallocate(NN%Xpredict)
      deallocate(NN%Ypredict)

      return
      end subroutine NN_predict

      subroutine init_weight_matrix(wt_matrix,size_l,size_lm1) 
      IMPLICIT NONE

      type(weight_matrix_type) :: wt_matrix
      INTEGER_T size_l
      INTEGER_T size_lm1
      INTEGER_T i,j

      if ((size_l.ge.1).and.(size_lm1.ge.1)) then
       wt_matrix%size_l=size_l
       wt_matrix%size_lm1=size_lm1
       ALLOCATE(wt_matrix%weights(size_l,size_lm1))
       do i=1,size_l
       do j=1,size_lm1
        wt_matrix%weights(i,j)=zero
       enddo 
       enddo 
       ALLOCATE(wt_matrix%bias(size_l))
       do i=1,size_l
        wt_matrix%bias(i)=zero
       enddo 
      else
       print *,"size_l or size_lm1 invalid"
       stop
      endif

      end subroutine init_weight_matrix

      subroutine scale_training_min_max( &
        Xtrain,Ytrain,Ytrain_scale, &
        n_input,n_output,n_total_samples, &
        min_Xtrain,max_Xtrain, &
        min_Ytrain,max_Ytrain)
      IMPLICIT NONE

      INTEGER_T n_input,n_output,n_total_samples
      REAL_T Xtrain(n_input,n_total_samples)
      REAL_T Ytrain(n_output,n_total_samples)
      REAL_T Ytrain_scale(n_output,n_total_samples)
      REAL_T min_Xtrain(n_input)
      REAL_T max_Xtrain(n_input)
      REAL_T min_Ytrain(n_output)
      REAL_T max_Ytrain(n_output)
      REAL_T Xrange
      REAL_T Yrange
      INTEGER_T i,r

      if (n_output.ge.1) then
       do i=1,n_output
        if (n_total_samples.ge.1) then
         min_Ytrain(i)=Ytrain(i,1)
         max_Ytrain(i)=Ytrain(i,1)
         do r=1,n_total_samples
          if (min_Ytrain(i).gt.Ytrain(i,r)) then
           min_Ytrain(i)=Ytrain(i,r)
          endif
          if (max_Ytrain(i).lt.Ytrain(i,r)) then
           max_Ytrain(i)=Ytrain(i,r)
          endif
         enddo ! r
         Yrange=max_Ytrain(i)-min_Ytrain(i)
         if (Yrange.gt.zero) then
          do r=1,n_total_samples
           Ytrain_scale(i,r)=(Ytrain(i,r)-min_Ytrain(i))/Yrange
          enddo 
         else
          print *,"Yrange invalid"
          stop
         endif
        else
         print *,"n_total_samples invalid"
         stop
        endif
       enddo ! i=1..n_output
      else
       print *,"n_output invalid"
       stop
      endif

      if (n_input.ge.1) then
       do i=1,n_input
        if (n_total_samples.ge.1) then
         min_Xtrain(i)=Xtrain(i,1)
         max_Xtrain(i)=Xtrain(i,1)
         do r=1,n_total_samples
          if (min_Xtrain(i).gt.Xtrain(i,r)) then
           min_Xtrain(i)=Xtrain(i,r)
          endif
          if (max_Xtrain(i).lt.Xtrain(i,r)) then
           max_Xtrain(i)=Xtrain(i,r)
          endif
         enddo ! r
         Xrange=max_Xtrain(i)-min_Xtrain(i)
         if (Xrange.gt.zero) then
          ! do nothing
         else
          print *,"Xrange invalid"
          stop
         endif
        else
         print *,"n_total_samples invalid"
         stop
        endif
       enddo ! i=1..n_input
      else
       print *,"n_input invalid"
       stop
      endif

      return
      end subroutine scale_training_min_max

      subroutine init_NN( &
       probtype, &
       verbose_in, &
       eps_activation_in,activation_id_in, &
       n_input_in,n_output_in,n_middle_in, &
       n_layers_in,n_total_samples_in,n_batch_samples_in, &
       Xtrain_in,Ytrain_in, &
       min_Xtrain_in,max_Xtrain_in, &
       NN)
      IMPLICIT NONE

      INTEGER_T probtype
      INTEGER_T verbose_in
      type(NN_type) :: NN
      REAL_T eps_activation_in
      INTEGER_T activation_id_in
      INTEGER_T n_input_in,n_output_in,n_middle_in
      INTEGER_T n_layers_in
      INTEGER_T n_total_samples_in,n_batch_samples_in
      REAL_T Xtrain_in(n_input_in,n_total_samples_in)
      REAL_T Ytrain_in(n_output_in,n_total_samples_in)
      REAL_T min_Xtrain_in(n_input_in)
      REAL_T max_Xtrain_in(n_input_in)
      INTEGER_T i,j,r,l,size_l,size_lm1
      INTEGER_T counter
      INTEGER_T steepest_descent_flag
      integer ( kind = 4 ) seed
      real (kind=8), allocatable :: r_normal_array(:,:)
      real (kind=8), allocatable :: r_normal_vec(:)
      REAL_T RELU_factor

      if (probtype.lt.1) then
       print *,"probtype invalid"
       stop
      endif

      NN%full_cycle_min_sweeps=1000

      NN%steepest_descent_flag=1

      if (NN%steepest_descent_flag.eq.0) then  ! ADAM
       NN%adam_alpha=0.001
      else if (NN%steepest_descent_flag.eq.1) then ! steepest descent
       NN%adam_alpha=0.0001
      else
       print *,"steepest_descent_flag invalid"
       stop
      endif

      NN%adam_beta1=0.9
      NN%adam_beta2=0.999
      NN%adam_eps=1.0E-8
      NN%adam_timestep=0

      NN%verbose=verbose_in

      NN%n_input=n_input_in
      NN%n_output=n_output_in
      NN%n_middle=n_middle_in
      NN%n_layers=n_layers_in
      NN%n_total_samples=n_total_samples_in
      NN%n_batch_samples=n_batch_samples_in
      NN%eps_activation=eps_activation_in
      NN%activation_id=activation_id_in

      NN%predict_flag=0

      ALLOCATE(NN%Xtrain(NN%n_input,NN%n_total_samples))
      ALLOCATE(NN%Ytrain(NN%n_output,NN%n_total_samples))
      ALLOCATE(NN%min_Xtrain(NN%n_input))
      ALLOCATE(NN%max_Xtrain(NN%n_input))

      do i=1,NN%n_input
       NN%min_Xtrain(i)=min_Xtrain_in(i)
       NN%max_Xtrain(i)=max_Xtrain_in(i)
      enddo

      do i=1,NN%n_input
      do r=1,NN%n_total_samples
       NN%Xtrain(i,r)=Xtrain_in(i,r)
      enddo
      enddo
      do i=1,NN%n_output
      do r=1,NN%n_total_samples
       NN%Ytrain(i,r)=Ytrain_in(i,r)
      enddo
      enddo

      ALLOCATE(NN%layer_array(NN%n_layers))
      ALLOCATE(NN%nodes_per_layer(NN%n_layers))
      NN%nodes_per_layer(1)=NN%n_input
      NN%nodes_per_layer(NN%n_layers)=NN%n_output
      do l=2,NN%n_layers-1
       NN%nodes_per_layer(l)=NN%n_middle
      enddo

      ALLOCATE(NN%X_random_basis(NN%n_input,NN%nodes_per_layer(2)))

      if (1.eq.0) then
       call init_RBF_centroids(NN)
      else
       call stochastic_BE_for_K_means(NN)
      endif
 
      do l=1,NN%n_layers
       ALLOCATE(NN%layer_array(l)%state_data( &
                NN%nodes_per_layer(l),NN%n_batch_samples))
       ALLOCATE(NN%layer_array(l)%hidden_data( &
                NN%nodes_per_layer(l),NN%n_batch_samples))
       ALLOCATE(NN%layer_array(l)%DCDS( &
                NN%nodes_per_layer(l),NN%n_batch_samples))

       size_l=NN%nodes_per_layer(l)
       if ((l.ge.2).and.(l.le.NN%n_layers)) then
        size_lm1=NN%nodes_per_layer(l-1)
       else if (l.eq.1) then
        size_lm1=1
       else
        print *,"l invalid l=",l
        stop
       endif

        ! wt_data bias can be init to 0.
        ! other wt_data values must be initialized with values pulled
        ! from standard normal distribution, mean=0, variance=1,
        ! then multiplied by sqrt(2/size_lm1) (for RELU)
       call init_weight_matrix(NN%layer_array(l)%wt_data,size_l,size_lm1)
       call init_weight_matrix(NN%layer_array(l)%DCDW,size_l,size_lm1)
       call init_weight_matrix(NN%layer_array(l)%moment1,size_l,size_lm1)
       call init_weight_matrix(NN%layer_array(l)%moment2,size_l,size_lm1)

       if (NN%activation_id.eq.2) then ! parametric RELU
        ALLOCATE(r_normal_array(size_l,size_lm1))
        seed=l
         ! mean 0 variance 1
        call r8mat_normal_01(size_l,size_lm1,seed,r_normal_array)
        RELU_factor=sqrt(two/size_lm1)
        do i=1,size_l
        do j=1,size_lm1
         NN%layer_array(l)%wt_data%weights(i,j)=r_normal_array(i,j)*RELU_factor
        enddo 
        enddo 
        DEALLOCATE(r_normal_array)
       else if (NN%activation_id.eq.0) then ! sigmoid
        ALLOCATE(r_normal_vec(size_l*size_lm1))
        seed=l+23
        call r8vec_uniform_01(size_l*size_lm1,seed,r_normal_vec)
        counter=1
        do i=1,size_l
        do j=1,size_lm1
         NN%layer_array(l)%wt_data%weights(i,j)=r_normal_vec(counter)/size_lm1
         counter=counter+1
        enddo 
        enddo 
        DEALLOCATE(r_normal_vec)
       else if (NN%activation_id.eq.4) then ! RBF
        if (l.eq.2) then
         do i=1,size_l
         do j=1,size_lm1
          NN%layer_array(l)%wt_data%weights(i,j)=one
         enddo
         enddo
         do i=1,size_l
          NN%layer_array(l)%wt_data%bias(i)=zero
         enddo
        else if (l.eq.3) then
         ALLOCATE(r_normal_vec(size_l*size_lm1))
         seed=23
         call r8vec_uniform_01(size_l*size_lm1,seed,r_normal_vec)
         counter=1
         do i=1,size_l
         do j=1,size_lm1
          NN%layer_array(l)%wt_data%weights(i,j)=r_normal_vec(counter)/size_lm1
          counter=counter+1
         enddo 
         enddo 
         DEALLOCATE(r_normal_vec)
         do i=1,size_l
          NN%layer_array(l)%wt_data%bias(i)=zero
         enddo
        else if (l.eq.1) then
         ! do nothing
        else
         print *,"l invalid"
         stop
        endif

       else
        print *,"do not know how to make an initial guess"
        stop
       endif

      enddo ! l=1..n_layers

      if (NN%activation_id.eq.4) then ! RBF
       call BICGSTAB_init(NN)
      endif

       ! sanity check
      if (1.eq.0) then
       if (probtype.eq.2) then
        do l=2,NN%n_layers
         size_l=NN%nodes_per_layer(l)
         size_lm1=NN%nodes_per_layer(l-1)
         if (l.eq.2) then
          do i=1,size_l
          do j=1,size_lm1
           NN%layer_array(l)%wt_data%weights(i,j)=one
          enddo
          enddo
          do i=1,size_l
           NN%layer_array(l)%wt_data%bias(i)=zero
          enddo
         else if (l.eq.3) then
          NN%layer_array(l)%wt_data%weights(1,1)=one/four
          NN%layer_array(l)%wt_data%weights(1,2)=-two/three
          NN%layer_array(l)%wt_data%weights(1,3)=one/two
          NN%layer_array(l)%wt_data%bias(1)=two
         else
          print *,"l invalid in sanity check"
          stop
         endif
        enddo
       else
        print *,"probtype invalid in sanity check"
        stop
       endif
      endif
        
      return
      end subroutine init_NN

      subroutine copy_weight_matrix(WT_dest,WT_source)
      IMPLICIT NONE

      type(weight_matrix_type) :: WT_dest
      type(weight_matrix_type) :: WT_source
      INTEGER_T size_l,size_lm1,i,j

      size_l=WT_source%size_l
      size_lm1=WT_source%size_lm1
      if ((size_l.ge.1).and.(size_lm1.ge.1)) then
       call init_weight_matrix(WT_dest,size_l,size_lm1)
       do i=1,size_l
       do j=1,size_lm1
        WT_dest%weights(i,j)=WT_source%weights(i,j)
       enddo
       enddo
       do i=1,size_l
        WT_dest%bias(i)=WT_source%bias(i)
       enddo
      else
       print *,"size_l or size_lm1 invalid"
       stop
      endif

      return
      end subroutine copy_weight_matrix

      subroutine copy_NN(NN_source,NN_dest)
      IMPLICIT NONE

      type(NN_type) :: NN_source
      type(NN_type) :: NN_dest
      INTEGER_T i,j,r,l

      NN_dest%full_cycle_min_sweeps=NN_source%full_cycle_min_sweeps
      NN_dest%steepest_descent_flag=NN_source%steepest_descent_flag
      NN_dest%adam_alpha=NN_source%adam_alpha
      NN_dest%adam_beta1=NN_source%adam_beta1
      NN_dest%adam_beta2=NN_source%adam_beta2
      NN_dest%adam_eps=NN_source%adam_eps
      NN_dest%adam_timestep=NN_source%adam_timestep
      NN_dest%n_input=NN_source%n_input
      NN_dest%n_output=NN_source%n_output
      NN_dest%n_middle=NN_source%n_middle
      NN_dest%n_layers=NN_source%n_layers
      NN_dest%n_total_samples=NN_source%n_total_samples
      NN_dest%n_batch_samples=NN_source%n_batch_samples
      NN_dest%eps_activation=NN_source%eps_activation
      NN_dest%activation_id=NN_source%activation_id

      NN_dest%verbose=NN_source%verbose

      NN_dest%predict_flag=NN_source%predict_flag

      ALLOCATE(NN_dest%Xtrain(NN_dest%n_input,NN_dest%n_total_samples))
      ALLOCATE(NN_dest%Ytrain(NN_dest%n_output,NN_dest%n_total_samples))

      ALLOCATE(NN_dest%min_Xtrain(NN_dest%n_input))
      ALLOCATE(NN_dest%max_Xtrain(NN_dest%n_input))

      do i=1,NN_dest%n_input
       NN_dest%min_Xtrain(i)=NN_source%min_Xtrain(i)
       NN_dest%max_Xtrain(i)=NN_source%max_Xtrain(i)
      enddo

      do i=1,NN_dest%n_input
      do r=1,NN_dest%n_total_samples
       NN_dest%Xtrain(i,r)=NN_source%Xtrain(i,r)
      enddo
      enddo
      do i=1,NN_dest%n_output
      do r=1,NN_dest%n_total_samples
       NN_dest%Ytrain(i,r)=NN_source%Ytrain(i,r)
      enddo
      enddo

      ALLOCATE(NN_dest%layer_array(NN_dest%n_layers))
      ALLOCATE(NN_dest%nodes_per_layer(NN_dest%n_layers))
      do l=1,NN_dest%n_layers
       NN_dest%nodes_per_layer(l)=NN_source%nodes_per_layer(l)
      enddo

      ALLOCATE(NN_dest%X_random_basis(NN_dest%n_input, &
                                      NN_dest%nodes_per_layer(2)))
      do i=1,NN_dest%n_input
      do j=1,NN_dest%nodes_per_layer(2)
       NN_dest%X_random_basis(i,j)=NN_source%X_random_basis(i,j)
      enddo
      enddo

      do l=1,NN_dest%n_layers
       ALLOCATE(NN_dest%layer_array(l)%state_data( &
                NN_dest%nodes_per_layer(l),NN_dest%n_batch_samples))
       ALLOCATE(NN_dest%layer_array(l)%hidden_data( &
                NN_dest%nodes_per_layer(l),NN_dest%n_batch_samples))
       ALLOCATE(NN_dest%layer_array(l)%DCDS( &
                NN_dest%nodes_per_layer(l),NN_dest%n_batch_samples))

       do i=1,NN_dest%nodes_per_layer(l)
        do r=1,NN_dest%n_batch_samples
         NN_dest%layer_array(l)%state_data(i,r)= &
            NN_source%layer_array(l)%state_data(i,r)
         NN_dest%layer_array(l)%hidden_data(i,r)= &
            NN_source%layer_array(l)%hidden_data(i,r)
         NN_dest%layer_array(l)%DCDS(i,r)= &
            NN_source%layer_array(l)%DCDS(i,r)
        enddo ! r
       enddo ! i

       call copy_weight_matrix(NN_dest%layer_array(l)%wt_data, &
         NN_source%layer_array(l)%wt_data)
       call copy_weight_matrix(NN_dest%layer_array(l)%DCDW, &
         NN_source%layer_array(l)%DCDW)
       call copy_weight_matrix(NN_dest%layer_array(l)%moment1, &
         NN_source%layer_array(l)%moment1)
       call copy_weight_matrix(NN_dest%layer_array(l)%moment2, &
         NN_source%layer_array(l)%moment2)

      enddo ! l=1..n_layers
        
      return
      end subroutine copy_NN

      subroutine delete_weight_matrix(wt_matrix)
      IMPLICIT NONE

      type(weight_matrix_type) :: wt_matrix

      DEALLOCATE(wt_matrix%weights)
      DEALLOCATE(wt_matrix%bias)

      end subroutine delete_weight_matrix

      subroutine delete_NN(NN)
      IMPLICIT NONE

      type(NN_type) :: NN
      INTEGER_T l

      DEALLOCATE(NN%X_random_basis)

      DEALLOCATE(NN%Xtrain)
      DEALLOCATE(NN%Ytrain)
      DEALLOCATE(NN%min_Xtrain)
      DEALLOCATE(NN%max_Xtrain)

      do l=1,NN%n_layers
       DEALLOCATE(NN%layer_array(l)%state_data)
       DEALLOCATE(NN%layer_array(l)%hidden_data)
       DEALLOCATE(NN%layer_array(l)%DCDS)
       call delete_weight_matrix(NN%layer_array(l)%wt_data)
       call delete_weight_matrix(NN%layer_array(l)%DCDW)
       call delete_weight_matrix(NN%layer_array(l)%moment1)
       call delete_weight_matrix(NN%layer_array(l)%moment2)
      enddo ! l=1,NN%n_layers

      DEALLOCATE(NN%layer_array)
      DEALLOCATE(NN%nodes_per_layer)

      if (NN%activation_id.eq.4) then ! RBF
       DEALLOCATE(NN%A)
       DEALLOCATE(NN%B)
      endif
        
      return
      end subroutine delete_NN

      subroutine forward_propagate(sample_start,sample_end,NN)
      IMPLICIT NONE

      INTEGER_T sample_start,sample_end

      type(NN_type) :: NN
      INTEGER_T l,i,j,r,rwrap,rstate

      if (NN%predict_flag.eq.0) then
       if ((sample_end.gt.NN%n_total_samples+NN%n_batch_samples-1).or. &
           (sample_start.gt.sample_end).or. &
           (sample_start.lt.1).or. &
           (sample_end-sample_start+1.gt.NN%n_batch_samples)) then
        print *,"sample_end or sample_start invalid"
        stop
       endif
      else if (NN%predict_flag.eq.1) then
       if ((sample_start.ne.1).or. &
           (sample_end.ne.NN%n_predict_samples)) then
        print *,"sample_end or sample_start invalid"
        stop
       endif
       do l=1,NN%n_layers
        ALLOCATE(NN%layer_array(l)%state_data_predict( &
                 NN%nodes_per_layer(l),NN%n_predict_samples))
        ALLOCATE(NN%layer_array(l)%hidden_data_predict( &
                 NN%nodes_per_layer(l),NN%n_predict_samples))
       enddo
      else
       print *,"predict_flag invalid"
       stop
      endif

      l=1
      do i=1,NN%n_input
       do r=sample_start,sample_end
        rstate=r-sample_start+1
        rwrap=r
        if (rwrap.gt.NN%n_total_samples) then
         rwrap=rwrap-NN%n_total_samples
        endif
        if (NN%predict_flag.eq.0) then
         NN%layer_array(l)%state_data(i,rstate)=NN%Xtrain(i,rwrap)
         NN%layer_array(l)%hidden_data(i,rstate)=NN%Xtrain(i,rwrap)
        else if (NN%predict_flag.eq.1) then
         NN%layer_array(l)%state_data_predict(i,r)=NN%Xpredict(i,r)
         NN%layer_array(l)%hidden_data_predict(i,r)=NN%Xpredict(i,r)
        else
         print *,"predict_flag invalid"
         stop
        endif
       enddo ! r=sample_start...sample_end
      enddo ! i=1..n_input
      do l=2,NN%n_layers
       do j=1,NN%nodes_per_layer(l)
        do r=sample_start,sample_end
         rstate=r-sample_start+1
         rwrap=r
         if (rwrap.gt.NN%n_total_samples) then
          rwrap=rwrap-NN%n_total_samples
         endif

         if (NN%predict_flag.eq.0) then
          NN%layer_array(l)%hidden_data(j,rstate)=zero
          do i=1,NN%nodes_per_layer(l-1)
           NN%layer_array(l)%hidden_data(j,rstate)= &
            NN%layer_array(l)%hidden_data(j,rstate)+ &
            NN%layer_array(l-1)%state_data(i,rstate)* &
            NN%layer_array(l)%wt_data%weights(j,i)
          enddo  ! i=1..size_lm1
          ! bias
          NN%layer_array(l)%hidden_data(j,rstate)= &
           NN%layer_array(l)%hidden_data(j,rstate)+ &
           NN%layer_array(l)%wt_data%bias(j)

          NN%layer_array(l)%state_data(j,rstate)= &
           activation_func(NN%layer_array(l)%hidden_data(j,rstate),r,l,j,NN) 
         else if (NN%predict_flag.eq.1) then
          NN%layer_array(l)%hidden_data_predict(j,r)=zero
          do i=1,NN%nodes_per_layer(l-1)
           NN%layer_array(l)%hidden_data_predict(j,r)= &
            NN%layer_array(l)%hidden_data_predict(j,r)+ &
            NN%layer_array(l-1)%state_data_predict(i,r)* &
            NN%layer_array(l)%wt_data%weights(j,i)
          enddo  ! i=1..size_lm1
          ! bias
          NN%layer_array(l)%hidden_data_predict(j,r)= &
           NN%layer_array(l)%hidden_data_predict(j,r)+ &
           NN%layer_array(l)%wt_data%bias(j)

          NN%layer_array(l)%state_data_predict(j,r)= &
           activation_func(NN%layer_array(l)%hidden_data_predict(j,r),r,l,j,NN) 
         else
          print *,"predict_flag invalid"
          stop
         endif

        enddo ! r=sample_start ... sample_end
       enddo ! j=1..nodes_per_layer(l)
      enddo ! l=2..NN%n_layers
      NN%cost=zero
      NN%cost_scale=zero
      l=NN%n_layers

      if (NN%n_output.ne.NN%nodes_per_layer(l)) then
       print *,"NN%n_output.ne.NN%nodes_per_layer(l)"
       stop
      endif

      do i=1,NN%n_output
       do r=sample_start,sample_end
        rstate=r-sample_start+1
        rwrap=r
        if (rwrap.gt.NN%n_total_samples) then
         rwrap=rwrap-NN%n_total_samples
        endif
        if (NN%predict_flag.eq.0) then
         NN%cost=NN%cost+ &
          (NN%layer_array(l)%state_data(i,rstate)-NN%Ytrain(i,rwrap))**2
         NN%cost_scale=NN%cost_scale+NN%Ytrain(i,rwrap)**2
        else if (NN%predict_flag.eq.1) then
         NN%cost=NN%cost+ &
          (NN%layer_array(l)%state_data_predict(i,r)-NN%Ypredict(i,r))**2
         NN%cost_scale=NN%cost_scale+NN%Ypredict(i,r)**2
         NN%Ypredict(i,r)=NN%layer_array(l)%state_data_predict(i,r)
        else
         print *,"predict_flag invalid"
         stop
        endif
       enddo ! r=sample_start ... sample_end
      enddo ! i=1..n_output

      if (NN%predict_flag.eq.1) then
       do l=1,NN%n_layers
        DEALLOCATE(NN%layer_array(l)%state_data_predict)
        DEALLOCATE(NN%layer_array(l)%hidden_data_predict)
       enddo
      else if (NN%predict_flag.eq.0) then
       ! do nothing
      else
       print *,"predict_flag invalid"
       stop
      endif

      return
      end subroutine forward_propagate

      subroutine back_propagate(sample_start,sample_end,NN)
      IMPLICIT NONE

      INTEGER_T sample_start,sample_end
      type(NN_type) :: NN
      type(NN_type) :: NN_plus
      type(NN_type) :: NN_minus
      INTEGER_T l,i,j,r,rwrap,rstate
      INTEGER_T brute_force_flag
      REAL_T dw

      brute_force_flag=0

      if (brute_force_flag.eq.0) then

       l=NN%n_layers
       do i=1,NN%n_output
        do r=sample_start,sample_end
         rstate=r-sample_start+1
         rwrap=r
         if (rwrap.gt.NN%n_total_samples) then
          rwrap=rwrap-NN%n_total_samples
         endif
         NN%layer_array(l)%DCDS(i,rstate)= &
          two*(NN%layer_array(l)%state_data(i,rstate)-NN%Ytrain(i,rwrap))
        enddo
       enddo ! i=1..n_output

        ! first calculate DCDS on the older generations.
       do l=NN%n_layers-1,1,-1

         ! p=parent  c=child
         ! if two layers,
         ! C=sum_r (S^c_r - Y_r)^2
         ! dC/dw^(c,p)=sum_r dC/dS^(c,r) dS^(c,r)/dw^(c,p)
         ! S^(c,r)=f(H^c,r)  H^c_r=sum w_c,p S^(p,r) + w^c_bias
         ! dS^(c,r)/dw^(c,p)=f'(H^c,r)S^(p,r)
         ! dC/dS^p,r=sum_c dC/dS^c,r dS^c,r/dS^p,r
         ! dS^c,r/dS^p,r=f'(H^c,r)w^(c,p)
         ! dC/dwbias^c=sum_r dC/dS^(c,r) f'(H^c_r)  
        do j=1,NN%nodes_per_layer(l)  ! traverse parent nodes
         do r=sample_start,sample_end
          rstate=r-sample_start+1
          rwrap=r
          if (rwrap.gt.NN%n_total_samples) then
           rwrap=rwrap-NN%n_total_samples
          endif
          NN%layer_array(l)%DCDS(j,rstate)=zero
          do i=1,NN%nodes_per_layer(l+1)  ! traverse children nodes
           NN%layer_array(l)%DCDS(j,rstate)= &  ! parent
            NN%layer_array(l)%DCDS(j,rstate)+ & ! parent
            NN%layer_array(l+1)%DCDS(i,rstate)* & ! child
            activation_prime(NN%layer_array(l+1)%hidden_data(i,rstate), &
             r,l+1,i,NN)* &
            NN%layer_array(l+1)%wt_data%weights(i,j)
          enddo ! i=1..size_lp1
         enddo ! r
        enddo ! j=1..size_l
       enddo ! l=NN%n_layers-1 ... 1

        ! second: calculate DCDW

        ! p=parent  c=child
        ! if two layers,
        ! C=sum_r (S^c_r - Y_r)^2
        ! dC/dw^(c,p)=sum_r dC/dS^(c,r) dS^(c,r)/dw^(c,p)
        ! S^(c,r)=f(H^c,r)  H^c_r=sum w_c,p S^(p,r) + w^c_bias
        ! dS^(c,r)/dw^(c,p)=f'(H^c,r)S^(p,r)
        ! dC/dS^p,r=sum_c dC/dS^c,r dS^c,r/dS^p,r
        ! dS^c,r/dS^p,r=f'(H^c,r)w^(c,p)
        ! dC/dwbias^c=sum_r dC/dS^(c,r) f'(H^c_r)  

       do l=NN%n_layers-1,1,-1

        do j=1,NN%nodes_per_layer(l)  ! traverse parent nodes
         do i=1,NN%nodes_per_layer(l+1)  ! traverse children nodes
          NN%layer_array(l+1)%DCDW%weights(i,j)=zero
         enddo
        enddo
        do i=1,NN%nodes_per_layer(l+1)  ! traverse children nodes
         NN%layer_array(l+1)%DCDW%bias(i)=zero
        enddo

         ! DC/Dw(c,p)
        do j=1,NN%nodes_per_layer(l)  ! traverse parent nodes
         do i=1,NN%nodes_per_layer(l+1)  ! traverse children nodes
          do r=sample_start,sample_end
           rstate=r-sample_start+1
           rwrap=r
           if (rwrap.gt.NN%n_total_samples) then
            rwrap=rwrap-NN%n_total_samples
           endif
           NN%layer_array(l+1)%DCDW%weights(i,j)= &
            NN%layer_array(l+1)%DCDW%weights(i,j)+ &
            NN%layer_array(l+1)%DCDS(i,rstate)* &  ! child
            activation_prime(NN%layer_array(l+1)%hidden_data(i,rstate), &
             r,l+1,i,NN)* &
            NN%layer_array(l)%state_data(j,rstate)
          enddo ! r=sample_start ... sample_end
         enddo ! i=1..size_lp1
        enddo ! j=1..size_l

         ! DC/Dwbias
        do i=1,NN%nodes_per_layer(l+1)  ! traverse children nodes
         do r=sample_start,sample_end
          rstate=r-sample_start+1
          rwrap=r
          if (rwrap.gt.NN%n_total_samples) then
           rwrap=rwrap-NN%n_total_samples
          endif
          NN%layer_array(l+1)%DCDW%bias(i)= &
            NN%layer_array(l+1)%DCDW%bias(i)+ &
            NN%layer_array(l+1)%DCDS(i,rstate)* &  ! child
            activation_prime(NN%layer_array(l+1)%hidden_data(i,rstate), &
             r,l+1,i,NN)
         enddo ! r=sample_start ... sample_end
        enddo ! i=1..size_lp1

       enddo ! l=NN%n_layers-1 ... 1

      else if (brute_force_flag.eq.1) then

       do l=2,NN%n_layers
        do j=1,NN%nodes_per_layer(l) ! children

         do i=1,NN%nodes_per_layer(l-1) ! parents
           ! copy_NN(source,dest)
          call copy_NN(NN,NN_plus)
          call copy_NN(NN,NN_minus)
          dw=1.0E-6
          NN_plus%layer_array(l)%wt_data%weights(j,i)= &
            NN%layer_array(l)%wt_data%weights(j,i)+dw
          NN_minus%layer_array(l)%wt_data%weights(j,i)= &
             NN%layer_array(l)%wt_data%weights(j,i)-dw
          NN%predict_flag=0
          NN_plus%predict_flag=0
          NN_minus%predict_flag=0
          call forward_propagate(sample_start,sample_end,NN_plus) 
          call forward_propagate(sample_start,sample_end,NN_minus) 
          NN%layer_array(l)%DCDW%weights(j,i)= &
              (NN_plus%cost-NN_minus%cost)/(two*dw)

          call delete_NN(NN_plus)
          call delete_NN(NN_minus)
         enddo ! i=1..size_lm1

          ! copy_NN(source,dest)
         call copy_NN(NN,NN_plus)
         call copy_NN(NN,NN_minus)
         dw=1.0E-6
         NN_plus%layer_array(l)%wt_data%bias(j)= &
            NN%layer_array(l)%wt_data%bias(j)+dw
         NN_minus%layer_array(l)%wt_data%bias(j)= &
            NN%layer_array(l)%wt_data%bias(j)-dw
         NN_plus%predict_flag=0
         NN_minus%predict_flag=0
         call forward_propagate(sample_start,sample_end,NN_plus) 
         call forward_propagate(sample_start,sample_end,NN_minus) 
         NN%layer_array(l)%DCDW%bias(j)= &
           (NN_plus%cost-NN_minus%cost)/(two*dw)

         call delete_NN(NN_plus)
         call delete_NN(NN_minus)

        enddo ! j=1..size_l
       enddo ! l=2..n_layer
      else
       print *,"brute_force_flag invalid"
       stop
      endif

      return
      end subroutine back_propagate

end module fully_connected_module

