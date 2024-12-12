      program main
      IMPLICIT NONE

      real*8 mypi,vert
      real*8 wetting_radius,vol,effective_hemi_r,static_angle,angle_term

      mypi=4.d0*atan(1.d0)
       !volume is given
      vol=1.0d0
       !wetting radius if angle=90 degrees
      effective_hemi_r=(3.0*vol/(2.0*mypi))**0.33333d0
      static_angle=mypi/2.0d0
      angle_term=2.0d0/3.0d0-cos(static_angle)+(cos(static_angle)**3)/3.0d0
      wetting_radius=(vol/(mypi*angle_term))**0.333333d0
       ! center of virtual sphere is a distance "vert" from the
       ! substrate.
      vert=-wetting_radius*cos(static_angle)

      print *,"vol= ",vol
      print *,"effective_hemi_r= ",effective_hemi_r
      print *,"static_angle= ",static_angle
      print *,"angle_term= ",angle_term
      print *,"wetting_radius= ",wetting_radius
      print *,"vert= ",vert
  
      end

