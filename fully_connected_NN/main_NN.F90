#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "SPACE.H"


      program main
      use fully_connected_module

      IMPLICIT NONE

      type(NN_type) :: NN
      REAL_T, allocatable :: Xtrain(:,:)
      REAL_T, allocatable :: Ytrain(:,:)
      REAL_T, allocatable :: Ytrain_scale(:,:)
      REAL_T, allocatable :: Xpredict(:,:)
      REAL_T, allocatable :: Ypredict_scale(:,:)
      REAL_T, allocatable :: min_Ytrain(:)
      REAL_T, allocatable :: max_Ytrain(:)
      REAL_T, allocatable :: min_Xtrain(:)
      REAL_T, allocatable :: max_Xtrain(:)
      REAL_T Ypredict_unscale
      REAL_T Yexpect
      REAL_T Yrange
      REAL_T eps_activation
       ! 0=sigmoid  1=tanh 2=parametric RELU 
       ! 3=polynomial least squares emulator
       ! 4=radial basis function interp
      INTEGER_T activation_id  
      INTEGER_T n_output,n_input,n_layers,n_middle
      INTEGER_T n_total_samples,n_batch_samples
      INTEGER_T n_predict
      REAL_T xlo,xhi,h
      REAL_T rand_num
      REAL_T grad_tol,cost_tol
      INTEGER_T max_steps
      INTEGER_T seed
      INTEGER_T i
       ! probtype==1 sin(x)
       ! probtype==2 cubic polynomial
      INTEGER_T probtype
      INTEGER_T verbose
      character*9 filename9

      verbose=0
      n_output=1
      n_input=1
      n_middle=20
      n_layers=3
      n_total_samples=1000
      n_batch_samples=1000
      n_predict=100

      ALLOCATE(Xtrain(n_input,n_total_samples))
      ALLOCATE(Ytrain(n_output,n_total_samples))
      ALLOCATE(Ytrain_scale(n_output,n_total_samples))
      ALLOCATE(Xpredict(n_input,n_predict))
      ALLOCATE(Ypredict_scale(n_output,n_predict))
      ALLOCATE(min_Ytrain(n_output))
      ALLOCATE(max_Ytrain(n_output))
      ALLOCATE(min_Xtrain(n_input))
      ALLOCATE(max_Xtrain(n_input))

      probtype=1
      if (probtype.eq.1) then
       xlo=zero
       xhi=two*Pi
      else if (probtype.eq.2) then
       xlo=-one
       xhi=one
      else
       print *,"probtype invalid"
       stop
      endif

      h=(xhi-xlo)/n_predict

!      eps_activation=(xhi-xlo)/10.0
      eps_activation=0.5 !parametric RELU: f(eps,x)=eps x if x<0 =x if x>0
       ! 0=sigmoid 1=tanh 2=parametric RELU 
       ! 3=polynomial least squares emulator
       ! 4=radial basis function
      activation_id=4

      print *,"activation_id=",activation_id

      seed=86456
      call srand(seed)
      do i=1,n_total_samples
       ! 0<=rand()<=1
       rand_num=rand()
       Xtrain(1,i)=xlo+(xhi-xlo)*rand_num
       Ytrain(1,i)=NN_test_func(Xtrain(1,i),probtype)
       if (verbose.ge.1) then
        print *,"Xtrain,Ytrain ",Xtrain(1,i),Ytrain(1,i)
       endif
      enddo

       ! 0<=Ytrain_scale<=1 
      call scale_training_min_max( &
        Xtrain,Ytrain,Ytrain_scale, &
        n_input,n_output,n_total_samples, &
        min_Xtrain,max_Xtrain, &
        min_Ytrain,max_Ytrain)

      if (verbose.ge.1) then
       do i=1,n_total_samples
        print *,"Xtrain,Ytrain_scale ",Xtrain(1,i),Ytrain_scale(1,i)
       enddo
      endif

      call init_NN( &
       probtype, &
       verbose, &
       eps_activation,activation_id, &
       n_input,n_output,n_middle, &
       n_layers, &
       n_total_samples,n_batch_samples, &
       Xtrain,Ytrain_scale, &
       min_Xtrain,max_Xtrain,NN)

      cost_tol=0.01
      grad_tol=1.0E-8
      max_steps=1500
      call ADAM_loop(cost_tol,grad_tol,max_steps,NN)

      do i=1,n_predict
       Xpredict(1,i)=xlo+(i-half)*h
       Yrange=max_Ytrain(1)-min_Ytrain(1)
       if (Yrange.gt.zero) then
        Ypredict_scale(1,i)= &
         (NN_test_func(Xpredict(1,i),probtype)-min_Ytrain(1))/Yrange
       else
        print *,"Yrange invalid"
        stop
       endif
      enddo ! i=1..n_predict

      call NN_predict(n_predict,Xpredict,Ypredict_scale,NN)

      print *,"NN_predict rel cost: ",sqrt(NN%cost/NN%cost_scale)
      print *,"prediction results in the file: NNPREDICT"

      write(filename9,'(A9)') 'NNPREDICT'
      open(unit=11,file=filename9)

      do i=1,n_predict
       Yrange=max_Ytrain(1)-min_Ytrain(1)
       if (Yrange.gt.zero) then
        Ypredict_unscale=min_Ytrain(1)+Yrange*Ypredict_scale(1,i)
        Yexpect=NN_test_func(Xpredict(1,i),probtype)
        write(11,*) Xpredict(1,i),Ypredict_unscale,Yexpect
       else
        print *,"Yrange invalid"
        stop
       endif
      enddo ! i=1..n_predict

      close(11)

      call delete_NN(NN)

      DEALLOCATE(Xtrain)
      DEALLOCATE(Ytrain)
      DEALLOCATE(Ytrain_scale)
      DEALLOCATE(Xpredict)
      DEALLOCATE(Ypredict_scale)
      DEALLOCATE(min_Xtrain)
      DEALLOCATE(max_Xtrain)
      DEALLOCATE(min_Ytrain)
      DEALLOCATE(max_Ytrain)

      return 
      end
