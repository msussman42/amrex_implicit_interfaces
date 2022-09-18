Module NeuralNetwork

  Implicit None
  Private
  Public :: Neural_Network

  !↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓Neural Networks variables↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
  Type Ragged_Vector
    real(8), Allocatable :: Vec(:)
  End Type Ragged_vector

  Type Ragged_Matrix
    real(8), Allocatable :: Mat(:,:)
  End Type Ragged_Matrix

  Type :: NN_activation
    Procedure(Sub_Interface), Pointer, NoPass :: activate => NULL()
  End Type NN_activation

  ! interface for choose activation type
  Interface
    Function Sub_Interface(n, X)
      Integer,  Intent(in) :: n
      real(8), Intent(in), Dimension(n) :: X
      real(8), Dimension(n) :: Sub_Interface
    End Function Sub_Interface
  End Interface

  ! neural network coefficients
  Type :: Neural_Network
    Integer :: n_inputs
    Integer :: n_outputs
    real(8), Allocatable :: Inputs(:)
    real(8), Allocatable :: Outputs(:)
    Integer :: layers
    Integer, Allocatable :: Layer_Size(:)
    Type(Ragged_Vector), Allocatable :: Activations(:)
    Type(Ragged_Vector), Allocatable :: Intercepts(:)
    Type(Ragged_Matrix), Allocatable :: Coefs(:)
    Character(10) :: Activation_type
    Character(10) :: out_Activation_type
    Type(NN_Activation) :: Activation
    Type(NN_Activation) :: Out_Activation
    Contains
      Procedure :: Initialization
      Procedure :: Predict
  End Type Neural_Network
  !↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑End Neural Network Variables↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑

  Contains

  Subroutine Initialization(self)
    Implicit None

    Class(Neural_Network) :: self
    Character(100) :: tmp
    Character(100) :: string
    real(8),allocatable :: line(:)
    Integer :: error
    Integer :: i,j

    tmp = 'nn_coef.dat'

    open(81,file=tmp,status='unknown')

    Read(81,*,iostat=error) string
    Read(81,*) self%layers
    Allocate(self%layer_size(self%layers))

    Read(81,*,iostat=error) string
    Read(81,*) self%layer_size
    self%n_inputs = self%layer_size(1)
    self%n_outputs = self%layer_size(self%layers)


    Allocate(self%activations(self%layers))
    Do i = 1,self%layers
      Allocate(self%activations(i)%vec(self%layer_size(i)))
    End do

    Read(81,*,iostat=error) string
    Allocate(self%intercepts(self%layers-1))
    Do i = 1,self%layers-1
      Allocate(self%intercepts(i)%vec(self%layer_size(i+1)))
      Read(81,*) self%intercepts(i)%vec
    End do

    Allocate(self%coefs(self%layers-1))
    Read(81,*,iostat=error) string
    Do i = 1,self%layers-1
      Allocate(self%coefs(i)%mat(self%layer_size(i+1),self%layer_size(i)))
      Read(81,*,iostat=error) string
      Allocate(line(self%layer_size(i)))
      Do j = 1,self%layer_size(i+1)
        Read(81,*) line
        self%coefs(i)%mat(j,:) = line
      End do
      DeAllocate(line)
    End do

    Read(81,*,iostat=error) string
    Read(81,*) self%Activation_type
    Read(81,*,iostat=error) string
    Read(81,*) self%out_Activation_type

    Close(81)

    If (Trim(self%Activation_type).eq.'logistic') Then
      self%activation%Activate => Activation_logistic
    Else If (Trim(self%Activation_type).eq.'tanh') Then
      self%activation%Activate => Activation_tanh
    Else If (Trim(self%Activation_type).eq.'softmax') Then
      self%activation%Activate => Activation_softmax
    Else If (Trim(self%Activation_type).eq.'relu') Then
      self%activation%Activate => Activation_ReLU
    Else If (Trim(self%Activation_type).eq.'identity') Then
      self%activation%Activate => Activation_identity
    Else
      Write(*,*) 'invalid activation type'
    End If

    If (Trim(self%out_Activation_type).eq.'logistic') Then
      self%out_activation%Activate => Activation_logistic
    Else If (Trim(self%out_Activation_type).eq.'tanh') Then
      self%out_activation%Activate => Activation_tanh
    Else If (Trim(self%out_Activation_type).eq.'softmax') Then
      self%out_activation%Activate => Activation_softmax
    Else If (Trim(self%out_Activation_type).eq.'relu') Then
      self%out_activation%Activate => Activation_ReLU
    Else If (Trim(self%out_Activation_type).eq.'identity') Then
      self%out_activation%Activate => Activation_identity
    Else
      Write(*,*) 'invalid output activation type'
    End If

  End Subroutine Initialization


  !↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓ Predict Subroutines↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
  function Predict(self, input)
    Implicit None
    Class(Neural_Network) :: self
    real(8) :: input(self%n_inputs)
    real(8) :: Predict(self%n_outputs)
    integer :: i

    self%activations(1)%vec = input

    Do i = 1, self%layers-2
      self%activations(i+1)%vec = matmul(self%coefs(i)%mat, &
          self%activations(i)%vec) + self%intercepts(i)%vec
      self%activations(i+1)%vec = self%activation%activate( &
          self%layer_size(i+1), self%activations(i+1)%vec)
    End do
    self%activations(Self%layers)%vec = &
        matmul(Self%coefs(Self%layers-1)%mat, &
        Self%activations(Self%layers-1)%vec) + &
        Self%intercepts(Self%layers-1)%vec
    Self%activations(Self%layers)%vec = &
        Self%out_activation%activate(Self%n_outputs, &
        Self%activations(Self%layers)%vec)

    Predict = Self%activations(Self%layers)%vec

  End function Predict
  !↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑End Prediction Subroutines↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑


  !↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓ Activation functions↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
  !-----------------------------------------------------
  ! In this part, the activation functions are defined.
  !   all the activation functions is called by a
  !   uniform procedure "Activation".
  !
  !-----------------------------------------------------
  ! I/O
  !-----------------------------------------------------
  !
  ! inputs:
  !   X:
  !     input vector, should be 1-D
  !   n:
  !     length of the input vector
  ! output:
  !    function result:
  !     output vector, should be 1-D with n length
  !
  !-----------------------------------------------------
  ! Activation functions
  !-----------------------------------------------------
  !
  ! Logistic function:
  !    f(x) = 1/(1+e^-x)
  !
  ! Tanh function:
  !    f(x) = tanh(x)
  !
  ! ReLU function:
  !    f(x) = max(0,x)
  !
  ! Identity function:
  !    f(x) = x
  !
  ! Softmax function:
  !    f(x)_i = exp(x_i)/sigma(x_i)
  !    Only used for classification.
  !    (Which means useless currently)
  !-----------------------------------------------------
  function Activation_logistic(n,X)
    Implicit None
    integer, Intent(in) :: n
    real(8), Intent(in), dimension(n) :: X
    real(8), dimension(n) :: Activation_logistic
    Activation_logistic = 1.0 / (1.0+exp(-X))
  End function Activation_logistic

  function Activation_tanh(n,X)
    Implicit None
    integer, Intent(in) :: n
    real(8), Intent(in), dimension(n) :: X
    real(8), dimension(n) :: Activation_tanh
    Activation_tanh = tanh(X)
  End function Activation_tanh

  function Activation_ReLU(n,X)
    Implicit None
    integer, Intent(in) :: n
    real(8), Intent(in), dimension(n) :: X
    real(8), dimension(n) :: Activation_ReLU
      Activation_ReLU = max(X,0.d0)
    End function Activation_ReLU

  function Activation_identity(n,X)
    Implicit None
    integer, Intent(in) :: n
    real(8), Intent(in), dimension(n) :: X
    real(8), dimension(n) :: Activation_identity
    Activation_identity = X
  End function Activation_identity

  function Activation_softmax(n,X)
    Implicit None
    Integer, Intent(in) :: n
    real(8), Intent(in), dimension(n) :: X
    real(8), dimension(n) :: tmp
    real(8), dimension(n) :: Activation_softmax
    tmp = exp(X - maxval(X))/sum(tmp)
    Activation_softmax = 1.0 / (1.0+exp(-X))
  End function Activation_softmax
  !↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑End Activation function↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑

End Module NeuralNetwork

Module DecisionTree

  Implicit None
  Private
  Public :: Decision_Tree

  !↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓Decision Tree Variables↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
  Type Nodes
    Integer  :: children_left
    Integer  :: children_right
    Integer  :: feature
    real(8) :: threshold
    real(8), Allocatable :: Values(:)
    ! Contains
  End Type Nodes

  type Trees
    Integer :: node_count
    Integer :: max_depth
    Type(Nodes), Allocatable :: Node(:)
  End type Trees

  Type :: Decision_Tree
    Integer :: n_inputs
    Integer :: n_outputs
    real(8), Allocatable :: Inputs(:)
    real(8), Allocatable :: Outputs(:)
    Type(Trees) :: Tree
  Contains
    Procedure :: Initialization
    Procedure :: Predict
  End Type Decision_Tree
  !↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑End Decision Tree Variables↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑


  Contains

  Subroutine Initialization(self)
    Implicit None

    Class(Decision_Tree) :: self
    Character(100) :: tmp
    Character(100) :: string
    Integer :: error
    Integer :: i

    tmp = 'dt_coef.dat'
    Open(82, file=tmp, status='unknown')

    Read(82,*,iostat=error) string
    Read(82,*) self%Tree%node_count
    Read(82,*,iostat=error) string
    Read(82,*) self%n_inputs
    Read(82,*,iostat=error) string
    Read(82,*) self%n_outputs
    Read(82,*,iostat=error) string
    Read(82,*) self%Tree%max_depth

    Allocate(self%Tree%node(self%Tree%node_count))

    Do i = 1,self%Tree%node_count
      Allocate(self%Tree%node(i)%values(self%n_outputs))
      Read(82,*,iostat=error) string
      Read(82,*) self%Tree%node(i)%children_left
      Read(82,*) self%Tree%node(i)%children_right
      Read(82,*) self%Tree%node(i)%feature
      Read(82,*) self%Tree%node(i)%threshold
      Read(82,*) self%Tree%node(i)%values
    End do

  End Subroutine Initialization

  !↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓ Predict Subroutines↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
  function Predict(self, input)
    Implicit None
    Class(Decision_Tree) :: self
    real(8) :: input(self%n_inputs)
    real(8) :: Predict(self%n_outputs)

    integer :: i,n

    n = 1
    Do i = 1, Self%Tree%max_depth
      if (Self%Tree%node(n)%feature .eq. -1) Exit
      if (input(Self%Tree%node(n)%feature) .le. SELF%Tree%node(n)%threshold) Then
        n = Self%Tree%node(n)%children_left
      else
        n = Self%Tree%node(n)%children_right
      End if
    End do

    Predict = Self%Tree%node(n)%values

  End function Predict
  !↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑End Prediction Subroutines↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑

End Module DecisionTree


Module RandomForest

  Implicit None
  Private
  Public :: Random_Forest 


  !↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓Decision Tree Variables↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
  Type Nodes
    Integer  :: children_left
    Integer  :: children_right
    Integer  :: feature
    real(8) :: threshold
    real(8), Allocatable :: Values(:)
    ! Contains
  End Type Nodes

  type Trees
    Integer :: node_count
    Integer :: max_depth
    Type(Nodes), Allocatable :: Node(:)
  End type Trees

  Type :: Random_Forest
    Integer :: n_inputs
    Integer :: n_outputs
    Integer :: tree_count
    real(8), Allocatable :: Inputs(:)
    real(8), Allocatable :: Outputs(:)
    Type(Trees), Allocatable :: Tree(:)
  Contains
    Procedure :: Initialization
    Procedure :: Predict
  End Type Random_Forest
  !↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑End Decision Tree Variables↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑


  Contains

  Subroutine Initialization(self)
    Implicit None

    Class(Random_Forest) :: self
    Character(100) :: tmp
    Character(100) :: string
    Integer :: error
    Integer :: i,j

    tmp = 'rf_coef.dat'
    Open(83, file=tmp, status='unknown')

    Read(83, *,iostat=error) string
    Read(83, *) self%tree_count

    Allocate(self%tree(self%tree_count))

    Do j = 1, self%Tree_count
      Read(83,*,iostat=error) string
      Read(83,*) self%Tree(j)%node_count
      Read(83,*,iostat=error) string
      Read(83,*) self%n_inputs
      Read(83,*,iostat=error) string
      Read(83,*) self%n_outputs
      Read(83,*,iostat=error) string
      Read(83,*) self%tree(j)%max_depth

      Allocate(self%Tree(j)%node(self%tree(j)%node_count))

      Do i = 1,self%Tree(j)%node_count
        Allocate(self%Tree(j)%node(i)%values(self%n_outputs))
        Read(83,*,iostat=error) string
        Read(83,*) self%Tree(j)%node(i)%children_left
        Read(83,*) self%Tree(j)%node(i)%children_right
        Read(83,*) self%Tree(j)%node(i)%feature
        Read(83,*) self%Tree(j)%node(i)%threshold
        Read(83,*) self%Tree(j)%node(i)%values
      End do
    End do

  End Subroutine Initialization

  !↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓ Predict Subroutines↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
  function Predict(self, input)
    Implicit None
    Class(Random_Forest) :: self
    real(8) :: input(self%n_inputs)
    real(8) :: Predict(self%n_outputs)

    integer :: i, j, n

    Predict = 0.0d0
    Do j = 1, Self%tree_count
      n=1
      Do i = 1, Self%tree(j)%max_depth
        if (Self%tree(j)%node(n)%feature .eq. -1) Exit
        if (input(Self%tree(j)%node(n)%feature) .le. Self%tree(j)%node(n)%threshold) Then
          n = Self%tree(j)%node(n)%children_left
        else
          n = Self%tree(j)%node(n)%children_right
        End if
      End do
      Predict = Predict + Self%tree(j)%node(n)%values
    End do

    Predict = Predict / Self%tree_count

  End function Predict
  !↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑End Prediction Subroutines↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑

End Module RandomForest
