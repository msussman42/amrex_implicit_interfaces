#ifndef _DGMOF_H_
#define _DGMOF_H_ 

#include <fstream>

#include <Box.H>
#include <Geometry.H>
#include <BoxArray.H>
#include <Array.H>

class Polyhedra {

public:

 // constructor
 Polyhedra(Array<Polyhedron> polyhedron_array);
 // destructor
 ~Polyhedra();

  // polyhedron_array_plus is intersection of this polyhedron_array
  // with the region where slope * (x-xcenter)+intercept>0
 void Cut_with_plane(Real slope[3],Real intercept,
  Real xcenter[3],
  Array<Polyhedron>& polyhedron_array_plus,
  Array<Polyhedron>& polyhedron_array_minus);

  // npoly is the number of polynomials to integrate.
 void Integrate_polynomial(Polynomial poly[],int npoly,
  Real& definite_integral[]);
