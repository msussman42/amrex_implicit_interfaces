#ifdef SINGLE_PRECISION
#  define sp           4
#else
#  define sp           8
#endif

program main
  use NeuralNetwork
  use DecisionTree
  use RandomForest
  implicit none

  Integer :: i
  Real(sp) :: vof
  Real(sp) :: exact_angle(2)
  Real(sp) :: initial_angle(2)
  Real(sp) :: inputs(3)

  ! Derived data type for NN, DT, RF
  Type(Neural_Network) :: NN
  Type(Decision_Tree) :: DT
  type(Random_Forest) :: RF

  ! Initialzation
  Call NN%Initialization()
  Call DT%Initialization()
  Call RF%Initialization()

  ! Read traning data
  open(11,file='exact_f.dat',status='old')
  open(12,file='exact_angle.dat',status='old')
  open(13,file='initial_angle.dat',status='old')
  Read(11,*)vof
  Read(12,*)exact_angle(1:2)
  Read(13,*)initial_angle(1:2)
  close(11)
  close(12)
  close(13)

  ! See prediction result, compare with training.py
  ! given a unit cell, [0,1]^3
  ! 1. centroid of the unit cell is: x_{cell}=(1/2,1/2,1/2)
  ! 2. "relative centroid" is x_{relative}=x_{act} - x_{cell} 
  !    where x_{act} is the 
  !    centroid of material located in the cell [0,1]^3.
  ! 3. n_{predict}=x_{relative}/||x_{relative}||
  ! 4. call slope_to_angle(n_{predict},angle_{predict})
  ! 5. angle_output=DT%predict(angle_{predict},F_{ref})
  ! 6. call angle_to_slope(angle_{output},n_{MachineLearning}
  inputs(1:2) = initial_angle(1:2)
  inputs(3) = vof
  Print *, 'Exact angle', exact_angle(1:2) 
  Print *, 'NN Fortran predict', NN%predict(inputs)
  Print *, 'DT Fortran predict', DT%predict(inputs)
  Print *, 'RF Fortran predict', RF%predict(inputs)

end program main
