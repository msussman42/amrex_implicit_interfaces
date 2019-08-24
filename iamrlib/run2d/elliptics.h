#include "axisym.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* function prototypes */
extern void adoublesetup(int*, PANEL[], double[], double[], double[], int*);
extern void calcE(double*, double*, double*, double*, double*, double*, double*);
extern void calcpot(double*, double*, double*, int*, double[], PANEL []);
extern void calcpotdn(double*, double*, double*, int*, double[], PANEL []);
extern double ellipke(double, double*, double*);
extern void forcetometers(PANEL [], int*, int*);
extern void invertmatrix(double[], int);
extern void panelconstants(int*, PANEL[], double[], double[], int*);
extern void panelsolvereal(int*, double[], PANEL[]);
extern void quicksort(double[], double[], int*, int*);
extern void setconstants(const int*, const double*, const int*, const double*);

/* internal function declarations */
int Pivot( double[], double[], int*, int* );
void Swap( double, double );

double ellipke(double m, double* k, double* e)  
{
  /*  Complete elliptic integral calculation
  L. Shure 1-9-88
  Modified to include the second kind by Bjorn Bonnevier
  from the Alfven Laboratory, KTH, Stockholm, Sweden
  Copyright 1984-2001 The MathWorks, Inc.
  $Revision: 5.17 $  $Date: 2001/04/15 12:01:44 $

  ELLIPKE uses the method of the arithmetic-geometric mean
  described in [1].

  References:
  [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
	Functions" Dover Publications", 1965, 17.6.

	Code adapted from MATLAB by Emhoff, with the addition of handling for
	values less than zero, as given by [1].*/

	double tol=1e-5;
	double fact=1.0;
	if (m < 0) { 
		m=-m;
		fact=sqrt(1.0+m); 
		m=m/(1.0+m);
	}

	double a0 = 1.0;
	double b0 = sqrt(1-m);
	double s0 = m;
	double mm = 1.0;
	double a1;
	double b1;
	double c1;
	double w1;
	double pwr2=1.0;

	while (mm > tol)  {
		a1 = (a0+b0)/2;
		b1 = sqrt(a0*b0);
		c1 = (a0-b0)/2;
		pwr2*=2.0;
		w1 = pwr2*c1*c1;
		mm = w1;
		s0 +=w1;
		a0 = a1;
		b0 = b1;
	}

	*k = pi/(2*a1*fact);
	*e = *k *fact*fact*(1-s0/2);

	return(*k);
};	/* ellipke */

/* dimensions of pan are *numpan */
/* dimensions of panelarray are (*numpan)**2 */
void panelconstants(int* numpan, PANEL pan[], double panelarray[], 
	double wts[npoints], int* invertyes)  
{
	int i;
	int j;
	int k;
	int ind;

	double z;
	double r;
	double r1;
	double r2;
	double R;
	double m;
	double km;
	double em;
	double temp;
	double denom;

	for (i=0; i<*numpan; i++)  {
		r1=pan[i].midr;
		if (abs(1.0-pan[i].type)<bigeps)  {  /* ith panel is Neumann */
			for (j=0; j<*numpan; j++)  {
				ind=*numpan *i+j;
				panelarray[ind]=0.0;
				temp=0.0;

				if (i==j)  {
					for (k=1; k<npoints; k++)  {
						r2=pan[j].rpoints[k];
						z=pan[i].midz-pan[j].zpoints[k];
						r=r1+r2;
						R=z*z+r*r;

						m=4.0*r1*r2/R;

						ellipke(m, &km, &em);
						em=em/(1.0-m);
						denom=R*sqrt(R);

            			/* nri*dphi/dr1 + nzi*dphi/dz1 */
						temp+=pan[j].rscaled[k]*wts[k]*(
								pan[i].rnrm*dphidr1(em,km,m,r1,r2,denom)+
								pan[i].znrm*dphidz1(em,z,denom));
					}
					panelarray[ind]=temp*pan[j].length*0.5-0.5*pi/r1;
				}

				else if (abs(1.0-pan[j].type)<bigeps)  {  /* jth panel is Neumann */
					for (k=0; k<npoints; k++)  {
						r2=pan[j].rpoints[k];
						z=pan[i].midz-pan[j].zpoints[k];
						r=r1+r2;
						R=z*z+r*r;
						m=4.0*r1*r2/R;
						ellipke(m, &km, &em);
						em=em/(1.0-m);
						denom=R*sqrt(R);

            			/* nri*dphi/dr1 + nzi*dphi/dz1 */
						temp+=pan[j].rscaled[k]*wts[k]*(
								pan[i].rnrm*dphidr1(em,km,m,r1,r2,denom)
								+pan[i].znrm*dphidz1(em,z,denom));
					}

					panelarray[ind]=temp*pan[j].length*0.5;
				}
				else  {  /* jth panel is Dirichlet */
					for (k=0; k<npoints; k++)  {
						r2=pan[j].rpoints[k];
						z=pan[i].midz-pan[j].zpoints[k];
						r=r1+r2;
						R=z*z+r*r;
						m=4.0*r1*r2/R;
						ellipke(m, &km, &em);
						em=em/(1.0-m);
						denom=R*R*sqrt(R)*(1-m);

	/* nzi*nzj*dphi/dz1dz2 + nri*nzj*dphi/dr1dz2 + */
        /* nzi*nrj*dphi/dz1dr2 + nri*nrj*dphi/dr1dr2 */
						temp+=pan[j].rscaled[k]*wts[k]*(
								pan[i].znrm*pan[j].znrm*dphidz2dz1(em,km,m,R,z,denom)+
								pan[i].znrm*pan[j].rnrm*dphidr2dz1(em,km,m,R,z,r1,r2,denom)+
								pan[i].rnrm*pan[j].znrm*dphidz2dr1(em,km,m,R,z,r1,r2,denom)+
								pan[i].rnrm*pan[j].rnrm*dphidr2dr1(em,km,m,R,r1,r2,denom));
					}

					panelarray[ind]=temp*pan[j].length*0.5;
				}
			}
		}

		else  {  /* ith panel is Dirichlet */
			for (j=0; j<*numpan; j++)  {
				ind=*numpan*i+j;
				panelarray[ind]=0.0;
				temp=0.0;

				if (i==j)  {
					for (k=0; k<npoints; k++)  {

						r2=pan[j].rpoints[k];
						z=pan[i].midz-pan[j].zpoints[k];
						r=r1+r2;
						R=z*z+r*r;

						m=4.0*r1*r2/R;
						ellipke(m, &km, &em);
						em=em/(1-m);
						denom=R*sqrt(R);

           		 		/* nrj*dphi/dr2+nzj*dphi/dz2 */
						temp+=pan[j].rscaled[k]*wts[k]*(
								pan[j].rnrm*dphidr2(em,km,m,r1,r2,denom)-
								pan[j].znrm*dphidz1(em,z,denom));
					}
					panelarray[ind]=temp*pan[j].length*0.5+0.5*pi/r1;
				}

				else if (abs(1.0-pan[j].type)<bigeps)  {  /* jth panel is Neumann */
					for (k=0; k<npoints; k++)  {
						r2=pan[j].rpoints[k];
						z=pan[i].midz-pan[j].zpoints[k];
						r=r1+r2;
						R=z*z+r*r;
						m=4.0*r1*r2/R;
						ellipke(m, &km, &em);
						temp+=pan[j].rscaled[k]*wts[k]*km/sqrt(R);
					}
					panelarray[ind]=temp*pan[j].length*0.5;
				}

				else  {  /* jth panel is Dirichlet */
					for (k=0; k<npoints; k++)  {
						r2=pan[j].rpoints[k];
						z=pan[i].midz-pan[j].zpoints[k];
						r=r1+r2;
						R=z*z+r*r;
						m=4.0*r1*r2/R;
						ellipke(m, &km, &em);
						em=em/(1-m);
						denom=R*sqrt(R);

            			/* nrj*dphi/dr2+nzj*dphi/dz2 */
						temp+=pan[j].rscaled[k]*wts[k]*(
								pan[j].rnrm*dphidr2(em,km,m,r1,r2,denom)-
								pan[j].znrm*dphidz1(em,z,denom));
					}
					panelarray[ind]=temp*pan[j].length*0.5;
				}
			}
		}
	}
	/* printf("A(i,j)= \n");
	for (i=0; i<*numpan; i++)  {
		for (j=0; j<*numpan; j++)  {
			ind=*numpan*i+j;
			printf("%11.4g", panelarray[ind]);
		}
	printf(" ... \n");	
	} */

	if (*invertyes) {
		invertmatrix(panelarray, *numpan);
	}
}	/* END panelconstants */

/* dimensions of pan,bside: *numpan
    "   "         panelarray: (*numpan)^2  */
              
void adoublesetup(int* numpan, PANEL pan[], double bside[], 
	double panelarray[], double wts[npoints], int* invertyes)
{
	int i, j, k;
	int panfromto, countpan;

	double z, r, r1, r2, R, m, km, em, dm;
	double denom, baseChunk;
	double ai, bi, ci, di, dGdn_in;
        double* RHS_gamma;
        double* RHS_alpha;
        int freeflag;

        RHS_gamma=(double *) calloc(*numpan,sizeof(double));
        RHS_alpha=(double *) calloc(*numpan,sizeof(double));

	for (i=0; i<*numpan; i++)  {		
		r1=pan[i].midr;
		RHS_gamma[i]=0.; RHS_alpha[i]=0;

		for (j=0; j<*numpan; j++)  {	
			countpan=*numpan*i +j;
			panelarray[countpan]=0.;
			panfromto=0;
			if ( (abs(1.0-pan[i].type)<bigeps) && (abs(1.0-pan[j].type)<bigeps) ) 
				panfromto=11;
			else if (abs(1.0-pan[i].type)<bigeps)
				panfromto=10;
			else if (abs(1.0-pan[j].type)<bigeps)
				panfromto=1;
			else if ( (abs(pan[i].type))<bigeps && (abs(pan[i].type)<bigeps) )
				panfromto=0;
			else
				printf("DANGERDANGERDANGER\n");

			ai=0.; bi=0.; ci=0.; di=0.;

			for (k=0; k<npoints; k++)  {
				r2=pan[j].rpoints[k];
				z=pan[i].midz-pan[j].zpoints[k];
				r=r1+r2;
				R=z*z+r*r;
				m=4.0*r1*r2/R;

				ellipke(m, &km, &em);
				dm=em/(1.0-m);	/* = E(m)/W = D(m) */

				baseChunk=pan[j].rscaled[k]*wts[k]*pan[j].length*0.5;

				denom=R*sqrt(R);
				dGdn_in=pan[j].rnrm*dphidr1(dm,km,m,r1,r2,denom)-
					pan[j].znrm*dphidz1(dm,z,denom);

				if (i != j) {
					switch(panfromto)
					{
						case(0):		/* dirchlet(D)-dirichlet. =G */
							ai += baseChunk*km/sqrt(R);
							RHS_gamma[i] += baseChunk*dGdn_in;

							break;
						case(1): 		/* D-neuman(N). =gradG.nj */
							bi += baseChunk*dGdn_in;

							break;
						case(10):		/* N-D. =gradG.ni */
							ci += baseChunk*km/sqrt(R);

							break;
						case(11):		/* N-N. =grad.gradG.ni.nj */
							denom=R*R*sqrt(R)*(1-m);

							di +=baseChunk*(
									pan[i].znrm*pan[i].znrm*dphidz2dz1(dm,km,m,R,z,denom)+
									pan[i].znrm*pan[i].rnrm*dphidr2dz1(dm,km,m,R,z,r1,r2,denom)+
									pan[i].rnrm*pan[i].znrm*dphidz2dr1(dm,km,m,R,z,r1,r2,denom)+
									pan[i].rnrm*pan[i].rnrm*dphidr2dr1(dm,km,m,R,r1,r2,denom));

							RHS_alpha[i] += baseChunk*dGdn_in;

							break;
						default:
							break;
					} 					/* end switch(panfromto) */
				} else {

				}
			} 						/* end Kdo */
			panelarray[countpan]= (ai+bi+ci+di);
		}							/* end GreenJ */
	}								/* end GreenI */

	for (i=0; i<*numpan; i++)  {		
		if (abs(1.0-pan[i].type)<bigeps) {		/* Neumann */
			bside[i]= pan[i].midpot*(0.5 - RHS_alpha[i]);
		} else {					/* Dirichlet */
			bside[i]= pan[i].midpot*(0.5 - RHS_gamma[i]);
		}
		/* printf("%9.2f %9.2f %9.2f %9.2f \n",
			RHS_alpha[i], RHS_gamma[i], pan[i].midpot, bside[i]); */
	}

	/* invert matrix */
	if (*invertyes) {
		invertmatrix(panelarray, *numpan);
	}
        free(RHS_gamma);
        free(RHS_alpha);
}	/* END adoublesetup */

void invertmatrix(double a[], int n)  
{	/*  Inverts matrix A with dimension n */
	int i;
	int j;
	int k;
	int p=0;
	int temp;
	int t1;
	int t2;

	double mult=0.0;

        int* ind;
        double* invert;
        double* diag;

        ind=(int *) calloc(n,sizeof(int));
        invert=(double *) calloc(n*n,sizeof(double));
        diag=(double *) calloc(n,sizeof(double));

	/*  Initialize to the identity */
	for (i=0; i<n; i++)  {
		for (j=0; j<n; j++)
			invert[i*n+j]=0.0;
		invert[i*n+i]=1.0;
	}
	for (i=0; i<n; i++)
		ind[i]=i;

	/*  Do Gaussian elimination to make A upper-diagonal */
	for (i=0; i<n; i++)  {
		p=i;
		for (j=i+1; j<n; j++)  {
			if (fabs(a[ind[j]*n+i]) > fabs(a[ind[p]*n+i]))
				p=j;
		}

		if (p != i)  {
			temp=ind[i];
			ind[i]=ind[p];
			ind[p]=temp;
		}

		if (fabs(a[ind[i]*n+i]) < eps)  {
			printf("elliptics:invertmatrix small err. Element %4i = A(%4i,%4i) \n", ind[i]*n+i+1,ind[i]+1,i+1);
			/*exit; */
		}
		t1=ind[i]*n;
		for (j=i+1; j<n; j++)  {
			t2=ind[j]*n;
			mult=a[t2+i]/a[t1+i];
			for (k=i; k<n; k++)  {
				a[t2+k]-=mult*a[t1+k];
			}
			for (k=0; k<j; k++)
				invert[t2+ind[k]]-=mult*invert[t1+ind[k]];
		}
	}
	/*  Do Gaussian elimination to make A diagonal */
	for (i=n-1; i>-1; i--)  {
		t1=ind[i]*n;
		for (j=i-1; j>-1; j--)  {
			t2=ind[j]*n;
			mult=a[t2+i]/a[t1+i];
			a[t2+i]=0.0;
			for (k=0; k<n; k++)  {
				invert[t2+k]-=mult*invert[t1+k];
			}
		}
	}
	/*  Divide inverted rows by the diagonal element in A */
	/*    to recover the identity in A, and then store in A */
	for(i=0; i<n; i++)
		diag[i]=a[ind[i]*n+i];

	for (i=0; i<n; i++)  {
		t1=i*n;
		t2=ind[i]*n;
		for (j=0; j<n; j++)  {
			a[t1+j]=invert[t2+j]/=diag[i];
		}
	}
        free(ind);
        free(invert);
        free(diag);
}	/* invertmatrix */

/* panelarray[(*numpan2)*(*numpan2)]   pan[*numpan2]  */
void panelsolvereal(int* numpan2, double panelarray[], 
	PANEL pan[])
{
	int i=0;
	int j=0;
	int k=0;

        double* gammas;
        double* phis;
        gammas=(double *) calloc(*numpan2,sizeof(double));
        phis=(double *) calloc(*numpan2,sizeof(double));

	for (i=0; i<*numpan2; i++)  {
		phis[i]=pan[i].midpot/chargeconstant;
		gammas[i]=0.0;
	}

  	/* Solve system of equations */
    /* Multiply inverted array by solution vector */
	for (i=0; i<*numpan2; i++)  {
		k=*numpan2*i;
		for (j=0; j<*numpan2; j++)
			gammas[i]+=panelarray[k+j]*phis[j];
	}

	for (i=0; i<*numpan2; i++)
		pan[i].str=gammas[i];	
        free(gammas);
        free(phis);
}	/* panelsolvereal */

/* calculates normal force on panel due to electric field 'E2' at 
	point (rpp,zpp) with 'pan' being the panel location. Called 
	when using hardwired charge at height r=0, z=qZpt */

void calcE( double* Fnorm2, double* E2, double* zp0, double* zp1, 
	double* rp0, double* rp1, double* qZpt ) 
{
	double delR, delZ, QtoPt, midR, midZ, dot1;
	double scnorm[2], spannorm[2];

	midR=0.5*(*rp0+*rp1);
	midZ=0.5*(*zp0+*zp1);
	delR=midR-0.;
	delZ=midZ-*qZpt;
	QtoPt =  sqrt(delR*delR + delZ*delZ);

	*E2 = kgauss*hardwirecharge/(QtoPt*QtoPt);
	/* *E2 = kgauss*hardwirecharge/(QtoPt); */

	/* calculating how much of this electric field is normal to the surface */
	scnorm[0] = delR/QtoPt;
	scnorm[1] = delZ/QtoPt;

	/* panel norm itself */
	delR=*rp0-*rp1;
	delZ=*zp0-*zp1;
	QtoPt =  sqrt(delR*delR + delZ*delZ);

	spannorm[0] = -delZ/QtoPt;
	spannorm[1] = delR/QtoPt;

	/*dot1=DOT_PRODUCT(scnorm,spannorm); */
	dot1=scnorm[0]*spannorm[0] + scnorm[1]*spannorm[1];

	*Fnorm2=-perm*fabs(dot1)*(*E2)*(*E2);
	/* *Fnorm2=-perm*(*E2)*(*E2); */
	return;
} 

/* calculates potential 'Pot2' at point (rpp,zpp) when 'pan' is sigma/mu values, 
     NOT  Dirichlet and Neumann*/
/*  pan[*numpan3] */
void calcpot(double* Pot2, double* zpp, double* rpp, int* numpan3, 
	double wts[npoints], PANEL pan[]) 
{
	int i;
	int k;

	double zp;
	double rp;

	double z;
	double r;
	double r2;
	double pot=0.0;
	double temp=0.0;

	double m;
	double R;
	double km;
	double em;
	double emm;
	double denom, rdist;

	rp=*rpp;
	zp=*zpp;

	for (i=0; i<*numpan3; i++)  {
		temp=0.0;

		if (1==skipfar)	/* skips far-away points */
			rdist=sqrt((pan[i].midr-rp)*(pan[i].midr-rp) + (pan[i].midz-zp)*(pan[i].midz-zp));
		else
			rdist=0;

		if (rdist<=ignoreimpact) {
			if (abs(1.0-pan[i].type)<bigeps)  {  /* panel is Neumann */
				for (k=0; k<npoints; k++)  {
					r2=pan[i].rpoints[k];
					z=zp-pan[i].zpoints[k];
					r=rp+r2;
					R=z*z+r*r;
					m=4.0*rp*r2/R;
					ellipke(m, &km,&em);
					temp+=pan[i].rscaled[k]*wts[k]*km/sqrt(R);
				}
				pot+=pan[i].str*temp*pan[i].length*0.5;
			}
			else  {  /* panel is Dirichlet */
				for (k=0; k<npoints; k++)  {
					r2=pan[i].rpoints[k];
					z=zp-pan[i].zpoints[k];
					r=rp+r2;
					R=z*z+r*r;
					m=4.0*rp*r2/R;

					ellipke(m, &km,&em);
					emm=em/(1-m);
					denom=sqrt(R)*R;

					temp+=pan[i].rscaled[k]*wts[k]*(
							pan[i].rnrm*dphidr1(emm,km,m,r2,rp,denom)-
							pan[i].znrm*dphidz1(emm,z,denom));
				}
				pot+=pan[i].str*temp*pan[i].length*0.5;
			}
		}
	}
	*Pot2=pot*chargeconstant;
} 

/* calculates potential 'Pot2' at point (rpp,zpp) when 'pan' is Dirichlet and Neumann
     values, NOT sigma/mu */
/* pan[*numpan4] */
void calcpotdn(double* Pot3, double* zp2, double* rp2, int* numpan4, double wts[npoints], 
	PANEL pan[]) 
{
	int i, k;

	double zp, rp;
	double z, r, r2, pot=0.0;
	double m, R, km, em, denom;
	double phiod, phion, dphiod, dphion;

	rp=*rp2; zp=*zp2;

	for (i=0; i<*numpan4; i++)  {
		phiod=0.; phion=0.;
		dphiod=0.; dphion=0.;

		for (k=0; k<npoints; k++)  
		{
			r2=pan[i].rpoints[k];
			z=zp-pan[i].zpoints[k];
			r=rp+r2;
			R=z*z+r*r;
			denom=R*sqrt(R);

			m=4.0*rp*r2/R;
			ellipke(m, &km,&em);

			if (abs(1.0-pan[i].type)<bigeps) {	/* Neumann */
				phion  +=pan[i].rscaled[k]*wts[k]*(
						pan[i].rnrm*dphidr1(em,km,m,rp,r2,denom)-
						pan[i].znrm*dphidz1(em,z,denom));

				dphion +=pan[i].rscaled[k]*wts[k]*km/sqrt(R);
			} else {
				phiod  +=pan[i].rscaled[k]*wts[k]*(
						pan[i].rnrm*dphidr1(em,km,m,rp,r2,denom)-
						pan[i].znrm*dphidz1(em,z,denom));

				dphiod +=pan[i].rscaled[k]*wts[k]*km/sqrt(R);
			}
		}
		if (abs(1.0-pan[i].type)<bigeps)  {	/* panel is Neumann */
			pot += pan[i].length*0.5*(pan[i].str*phion - pan[i].midpot*dphion);
		} else {			/* panel is Dirichlet */
			pot += pan[i].length*0.5*(pan[i].midpot*phiod - pan[i].str*dphiod);
		}
	}
	*Pot3=pot*chargeconstant;
} 

/* Automatically copies Fortran global.f90 constants to c++ equivalent. */
void setconstants(const int* skipfar2, const double* ignoreimpact2, 
	const int* recastdist2, const double* hardwirecharge2) 
{
	skipfar=*skipfar2;
	ignoreimpact=*ignoreimpact2;
	recastdist=*recastdist2;
	hardwirecharge=*hardwirecharge2;
} 

/* Quicksort:  Sort an array a, using the quicksort algorithm. a2 is a second column
     that gets copied along when sorting a. */
void quicksort( double a[], double a2[], int* first, int* last ) 
{
	int pivot2;

	if( first < last ) {
		pivot2 = Pivot( a, a2, first, last );
		quicksort( a, a2, first, &pivot2-1 );
		quicksort( a, a2, &pivot2+1, last );
	}
}


/*  Pivot:  Find and return the index of pivot element. */
int Pivot( double a[], double a2[], int* first, int* last ) 
{
	int i;
	int  p = *first;
	double pivot = a[*first];

	for( i = *first+1 ; i <= *last ; i++ ) {
		if( a[i] <= pivot ) {
			p++;
			Swap( a[i], a[p] );
			Swap( a2[i], a2[p] ); 
                             /* swaps a2 based on the ordering of a */
		}
	}

	Swap( a[p], a[*first] );
	Swap( a2[p], a2[*first] );

	return p;
}


/*  Swap:  Swap two item (by reference). */
void Swap( double v1, double v2 )
{
	double  tmpVal;

	tmpVal = v1;
	v1 = v2;
	v2 = tmpVal;
}
