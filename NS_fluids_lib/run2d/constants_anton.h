#ifndef constants_anton_h
#define constants_anton_h

# define npoints 10
/* NOTE: need to make sure 'constants.h:npoints' has the same value as global.f90:npoints */
double ignoreimpact;
double hardwirecharge;
int recastdist;
int skipfar;
/*  c++ equivalents to "global.f90" values. Copied over in routine "setconstants". */

const double pi=3.14159265359;
const double perm=8.8542e-12;	/* epsilon_0 */
const double elec=1.602178e-19;
const double invmass=1.0/2.18e-25;
/* const double chargeconstant=elec/(2.0*pi*pi*perm); */
const double chargeconstant=1.602178e-19/(2.0*3.14159265359*3.14159265359*8.8542e-12);
const double bigeps=1e-6; 
const double eps=1e-14; /* Epsilon to account for floating point error */
const double kgauss=8.99e9;		/* electrostatic Coulomb constant, [Nm^2/C^2] */

#endif /* !defined(constants_anton_h) */
