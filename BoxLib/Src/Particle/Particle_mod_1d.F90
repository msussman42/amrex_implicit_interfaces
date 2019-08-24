module BoxLib_particle_module

  use iso_c_binding
  use bl_fort_module, only : c_real

  implicit none

  private

  public :: BoxLib_particle_set_position, BoxLib_particle_get_position

contains

  subroutine BoxLib_particle_set_position (particles, ns, np, x) &
       bind(c,name='BoxLib_particle_set_position')
    integer(c_int)  , intent(in   ), value :: ns, np
    real(c_real), intent(inout)        :: particles(ns,np)
    real(c_real), intent(in   )        :: x(np)

    integer :: i

    do i = 1, np
       particles(1,i) = x(i)
    end do
  end subroutine BoxLib_particle_set_position

  subroutine BoxLib_particle_get_position (particles, ns, np, x) &
       bind(c,name='BoxLib_particle_get_position')
    integer(c_int)  , intent(in   ), value :: ns, np
    real(c_real), intent(in   )        :: particles(ns,np)
    real(c_real), intent(  out)        :: x(np)

    integer :: i

    do i = 1, np
       x(i) = particles(1,i)
    end do
  end subroutine BoxLib_particle_get_position

end module BoxLib_particle_module
