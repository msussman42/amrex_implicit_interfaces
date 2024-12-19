      program main
      IMPLICIT NONE

      real*8 mypi,vert
      real*8 wetting_radius,vol,effective_hemi_r,static_angle,angle_term
      real*8 height,virtual_sphere_radius
      real*8 sanity_vol
      real*8 sanity_vol2

      mypi=4.d0*atan(1.d0)
       !volume is given
      vol=8.5e-3
       !wetting radius if angle=90 degrees
      effective_hemi_r=(3.0*vol/(2.0*mypi))**0.33333d0
      static_angle=mypi/2.0d0
      static_angle=mypi*50.0/180.0
      angle_term=(2.0+cos(static_angle))*((1.0-cos(static_angle))**2)/
     &  (3.0*sin(static_angle)**3)

      wetting_radius=(vol/(mypi*angle_term))**0.333333d0
 
       ! a/r=sin(theta)
      virtual_sphere_radius=wetting_radius/sin(static_angle)

       ! center of virtual sphere is a distance "vert" from the
       ! substrate.
      vert=-virtual_sphere_radius*cos(static_angle)
      height=virtual_sphere_radius-abs(vert)

      print *,"vol= ",vol
      print *,"effective_hemi_r= ",effective_hemi_r
      print *,"static_angle= ",static_angle
      print *,"angle_term= ",angle_term
      print *,"wetting_radius= ",wetting_radius
      print *,"vert= ",vert
      print *,"height= ",height
      print *,"virtual_sphere_radius= ",virtual_sphere_radius
      sanity_vol=mypi*height*(3.0*(wetting_radius**2)+height**2)/6.0
      print *,"sanity_vol=",sanity_vol
      sanity_vol2=mypi*(height**2)*(3.0*virtual_sphere_radius-height)/3.0
      print *,"sanity_vol2=",sanity_vol2
 
      print *,"surface tension: ",72.8*(1.0-cos(static_angle)) 
      end

