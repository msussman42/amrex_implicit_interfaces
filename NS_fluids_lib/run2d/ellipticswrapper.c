#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "axisym.h"
#include "derivs.h"
#include "elliptics.h"

extern void adoublesetupwrap_(int*, PANEL[], double[], double[], double[], int*);
extern void panelconstantswrap_(int*, PANEL[], double[], double[], int*);
extern void panelsolverealwrap_(int*, double[], PANEL[]);
extern void calcewrap_(double*, double*, double*, double*, double*, double*, double*);
extern void calcpotwrap_(double*, double*, double*, int*, double[], PANEL[]);
extern void calcpotdnwrap_(double*, double*, double*, int*, double[], PANEL[]);
extern void setconstantswrap_(const int*, const double*, const int*, const double*);
extern void quicksortwrap_(double[], double[], int*, int*);

/* pan[*numpan]  panelarray[(*numpan)*(*numpan)]  */
void panelconstantswrap_(int* numpan, PANEL pan[], 
	double panelarray[], double wts[npoints], int* invertyes)  
{	panelconstants(numpan, pan, panelarray, wts, invertyes);
}	

void adoublesetupwrap_(int* numpan, PANEL pan[], double bside[], 
	double panelarray[], double wts[npoints], int* invertyes)
{	adoublesetup(numpan, pan, bside, panelarray, wts, invertyes);
}	

void panelsolverealwrap_(int* numpan2, double panelarray[], 
	PANEL pan[])
{	panelsolvereal(numpan2, panelarray, pan);
}	

void calcewrap_(double* Fnorm2, double* Elocal, double* zp0, double* zp1, double* rp0, 
	double* rp1, double* qZpt) 
{	calcE(Fnorm2, Elocal, zp0, zp1, rp0, rp1, qZpt); 
}	

void calcpotwrap_(double* Pot, double* zp, double* rp, int* numpan, double wts[npoints], 
	PANEL pan[])
{	calcpot(Pot, zp, rp, numpan, wts, pan); 
}	

void calcpotdnwrap_(double* Pot2, double* zp2, double* rp2, int* numpan2, 
	double wts2[npoints], PANEL pan2[])
{	calcpotdn(Pot2, zp2, rp2, numpan2, wts2, pan2); 
}	

void setconstantswrap_(const int* SkipFar, const double* IgnoreImpact, 
	const int* RecastDist, const double* HardwireCharge)
{	setconstants(SkipFar, IgnoreImpact, RecastDist, HardwireCharge);
}	
