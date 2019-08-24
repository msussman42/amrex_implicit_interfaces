c-----------------------------------------------------------------------
c
c     General nondimensional small and large numbers:
c
      real*8     small,       large
      parameter( small=1.d-9, large=1.d+9 )

c-----------------------------------------------------------------------
c
c     NONDIMENSIONAL small/large numbers:
c
      
      real*8     rho_min,          vel_min,       pre_min
      parameter( rho_min=small,    vel_min=small, pre_min=small )
      
      real*8     x_min,               t_min
      parameter( x_min = small,       t_min = small )
      real*8     x_infty,             t_infty
      parameter( x_infty = large,     t_infty = large )
