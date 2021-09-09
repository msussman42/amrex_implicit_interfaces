module constants  
 use amrex_fort_module, only : amrex_real
 implicit none 

 real(amrex_real), parameter :: M_PI = 3.14159265358979323846
 real(amrex_real), parameter :: M_EPSILON = 10E-6
   
end module constants 

module IntfPtc_module
 use amrex_fort_module, only: amrex_real, amrex_particle_real
 use iso_c_binding ,    only: c_int
   
 implicit none
    !!public particle_t?
    
 type, bind(C)  :: particle_t
  real(amrex_particle_real) :: pos(3)     !< Position
  
  real(amrex_particle_real) :: phi(1)     !< LS value
  real(amrex_particle_real) :: G_Ptc(3)   !< for RK-iterations
  
  integer(c_int)            :: id         !< Particle id
  integer(c_int)            :: cpu        !< Particle cpu
  
  integer(c_int)            :: sorted     !< Particle is in the right cell
  integer(c_int)            :: i_cell     !< Particle cell x
  integer(c_int)            :: j_cell     !< Particle cell y
  integer(c_int)            :: k_cell     !< Particle cell z
  
 end type particle_t
    
end module IntfPtc_module