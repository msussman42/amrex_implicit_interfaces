#ifndef axisym_h
#define axisym_h

#include <math.h>               /*  Math library */
#include "constants_anton.h"

typedef struct PANEL {
	double type;		/* D/N */
	double z0;		/* panel z start */
	double z1;		/* panel z end */
	double r0;		/* panel r start */
	double r1;	
	double str; 	/* amount of charge. (SOLVE for) */
	double znrm;	/* z normal */
	double rnrm;
	double midz;
	double midr;
	double length;
	double midpot; 	/* panel potential or flux */
	double zpoints[npoints]; /* z locations of sub panel integration points */
	double rpoints[npoints];
	double rscaled[npoints]; /* array of rpoints/r_panel_center */
} PANEL;

#endif /* defined(axisym_h) */
