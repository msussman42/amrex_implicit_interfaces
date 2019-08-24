#ifndef derivs_h
#define derivs_h

#include "constants_anton.h"

/*****************************************************************************/
/* First order derivatives */

/*  Derivative with respect to z1 - flip the sign of z for d/dz2 */
inline double dphidz1(double dm, double z, double denom)  {
/* Denom = R*sqrt(R) */
	return(-dm*z/denom);
}

/*  Derivative with respect to r1, reverse r1 and r2 to get d/dr2  */
inline double dphidr1(double dm, double km, double m,
double r1, double r2, double denom)  {
/* Denom = R*sqrt(R) */
	if (r1 < eps)
		return(0);
	else if (r2 < eps)
		return(-r1*dm/denom);

	return((2*r2*(dm-km)/m-(r1+r2)*dm)/denom);
}

/*  Derivative with respect to r2, reverse r1 and r2 to get d/dr1 */
inline double dphidr2(double dm, double km, double m, 
double r1, double r2, double denom)  {
/* Denom = R*sqrt(R) */
	if (r2 < eps)
		return(0);
	else if (r1 < eps)
		return(-r2*dm/denom);

	return((2*r1*(dm-km)/m-(r1+r2)*dm)/denom);
}

/*****************************************************************************/
/* Second order derivatives 

  Second Derivative with respect to z2 and z1 */
inline double dphidz2dz1(double dm, double km, double m, double R, double z, 
double denom)  {
/* Denom = R*R*sqrt(R)*(1-m) */
	return((z*z*km+(R*(1-m)-2*z*z*(2-m))*dm)/denom);
}

/*  Second derivative with respect to z2 and r1 */
inline double dphidz2dr1(double dm, double km, double m, double R, double z, 
double r1, double r2, double denom)  {
/* Denom = R*R*sqrt(R)*(1-m) */
	if (r1 < eps)
		return(0);

	double r=r1+r2;
	return(-z*(0.5*R*(km-dm)/r1-r*km-(2*r2-2*r*(2-m))*dm)/denom);
}

/*  Second derivative with respect to r2 and z1 */
inline double dphidr2dz1(double dm, double km, double m, double R, double z, 
double r1, double r2, double denom)  {
/* Denom = R*R*sqrt(R)*(1-m) */
	if (r2 < eps)
		return(0);

	double r=r1+r2;
	return(z*(0.5*R*(km-dm)/r2-r*km-(2*r1-2*r*(2-m))*dm)/denom);
}

/*  Second derivative with respect to r2 and r1 */
inline double dphidr2dr1(double dm, double km, double m, double R, 
double r1, double r2, double denom)  {
/* Denom = R*R*sqrt(R)*(1-m) */
	if (r1 < eps || r2 < eps)
		return(0);

	double r=r1+r2;
	double rs=r*r;
	return(-(2*rs*(dm-km)/m+(R+rs)*km+(2*rs*(m-1)-R*(1+m))*dm)/denom);
}

/* Second derivative with respect to r2 and z2 */
inline double dphidz2dr2(double dm, double km, double m, double R, double z, 
double r1, double r2, double denom)  {
/* Denom = R*R*sqrt(R)*(1-m) */
	if (r2 < eps)
		return(0);

	double r=r1+r2;
	return(-z*(0.5*R*(km-dm)/r2-r*km-(2*r1-2*r*(2-m))*dm)/denom);
}

/* Second derivative with respect to r2 and r2 */
inline double dphidr2dr2(double dm, double km, double m, double R, double z, 
double r1, double r2, double denom)  {
/* Denom = R*R*sqrt(R)*(1-m) */
	double zs=z*z;
	double r2s=r2*r2;
	return(((km-dm)*R*R*(1.0-m)/(2.0*r2s)+(R*r1/r2+2.0*(zs*m-R+2.0*r2s))*dm+zs*km)/denom);
}

/*****************************************************************************/
/* Third order derivatives

 Third derivative with respect to z2, z2, and z1 */
inline double dphidz2dz2dz1(double dm, double km, double m, double R,
double z, double denom)  {
  /* denom=R*R*R*sqrt(R)*(1-m)*(1-m) */
	double zs=z*z;
	return(((R*(12.0-18.0*m)-23*zs*(1.0-m)+m*m*(6*R-8*zs))*dm
			+(4.0*zs*(2.0-m)+3*R*(m-1.0))*km)*z/denom);
}

/* Third derivative with respect to z2, z2, and r1 */
inline double dphidz2dz2dr1(double dm, double km, double m, double z, 
double r1, double r2, double denom)  {
  /* denom=R*R*R*sqrt(R)*(1-m)*(1-m) */
	double r=r1+r2;
	double rs=r*r;
	double zs=z*z;

	if (r1 < eps)
		return(0.0);

	return((
			2*(km - dm)*r2*(rs - 2*zs)/m +
			dm*(
				r*(2*rs*(m*m - 3*m + 2) + zs*(-6*m*m + 17*m - 19)) +
				2*r2*(zs*(7 - m) + rs*m)
			) +
			km*(
				rs*r*(m-1) - 2*r2*(rs + 2*zs) + r*zs*(7 - 3*m)
			)
		)/denom);
}

/* Third derivative with respect to z2, r2, and z1 */
inline double dphidz2dr2dz1(double dm, double km, double m, double z, 
double r1, double r2, double denom)  {
  /* denom=R*R*R*sqrt(R)*(1-m)*(1-m) */
	double r=r1+r2;
	double rs=r*r;
	double zs=z*z;

	if (r2 < eps)
		return(0.0);

	if (m < eps)
    /*return((
      dm*(
        r*(2*rs*(m*m - 3*m + 2) + zs*(-6*m*m + 17*m - 19)) +
        2*r1*(zs*(7 - m) + rs*m)
      ) +
      km*(
        rs*r*(m-1) - 2*r1*(rs + 2*zs) + r*zs*(7 - 3*m)
      )
    )/denom);*/
		return(-0.5*(3*3.14*r2*(rs-4*zs))/denom);

	return(-(
			2*(km - dm)*r1*(rs - 2*zs)/m +
			dm*(
				r*(2*rs*(m*m - 3*m + 2) + zs*(-6*m*m + 17*m - 19)) +
				2*r1*(zs*(7 - m) + rs*m)
			) +
			km*(
				rs*r*(m-1) - 2*r1*(rs + 2*zs) + r*zs*(7 - 3*m)
			)
		)/denom);
}

/* Third derivative with respect to z2, r2, and r1 */
inline double dphidz2dr2dr1(double dm, double km, double m, double R,
double z,double r1, double r2, double denom)  {
  /* denom=R*R*R*sqrt(R)*(1-m)*(1-m) */
	double r=r1+r2;
	double rs=r*r;

	if (r1 < eps)
		return(0.0);

	return((
			6*(km - dm)*rs/m -
			dm*(R*(2*m*m - 7*m - 3) + ((19 - 8*m)*m - 9)*rs) -
			km*(R*(3 + m) + 2*(3 - 2*m)*rs)
		)*z/denom);
}

/* Third derivative with respect to r2, r2, and z1 */
inline double dphidr2dr2dz1(double dm, double km, double m, double R,
double z, double r1, double r2, double denom)  {
  /* denom=R*R*R*sqrt(R)*(1-m)*(1-m) */
	double r=r1+r2;
	double r2s=r2*r2;
	double rs=r*r;
	double zs=z*z;
	double mm1=1-m;

	return((
			(
				2*R*R*(m*(2 + m) - 1) +
				16*R*mm1*mm1*r2s -
				4*r2s*(4*(3 - m)*(rs - 2*r2s) + mm1*(11 - 8*m)*zs)
			)*dm +
			(
				R*R*(2 - 3*m) -
				8*R*mm1*r2s +
				8*r2s*(rs - 2*r2s + (3 - 2*m)*zs)
			)*km
		)*(-z/(4*r2s*denom))); 
}

/* Third derivative with respect to r2, r2, and r1 */
inline double dphidr2dr2dr1(double dm, double km, double m, double R,
double r1, double r2, double denom)  {
	if (r1 < eps)
		return(0.0);

	double r1p2=r1*r1;
	double r2p2=r2*r2;
	double r2p4=r2p2*r2p2;

	double r=r1+r2;
	double rp3=r*r*r;
	double mp2=m*m;

  /* denom=R*R*R*sqrt(R)*(1-m)*(1-m) */
	return((
			(4*(dm - km)*(R*r*(r1 - r2) + 3*r2p4))/m - 
			(
				R*(R*m*(7*m - 4) + 2*mp2*r*r1 - (m - 8)*r1p2) - 
				(10*R*mp2*r - 2*(9 + m*(8*m - 19))*rp3 + R*(9 + 4*m)*r1)*r2 + 
				R*(29*m - 8)*r2p2
			)*dm -
			(
				m*(5*R*R + 8*rp3*r2 - R*r*(2*r1 + 11*r2)) - 
				12*r*r2*r2p2 - R*(6*r1p2 + 7*r1*r2 - 6*r2p2)
			)*km
		)/(2*r2*denom));
}

/* Third derivative with respect to r2, r2, and r2 */
inline double dphidr2dr2dr2(double dm, double km, double m, double R,
double r1, double r2, double denom)  {
  /* denom=R*R*R*sqrt(R)*(1-m)*(1-m) */
	double Rp2=R*R;
	double Rp3=Rp2*R;
	double r=r1+r2;
	double rp3=r*r*r;

	double r1p2=r1*r1;
	double r2p2=r2*r2;
	double r2p3=r2p2*r2;
	double r2p4=r2p2*r2p2;
	double r2p5=r2p4*r2;

  /* denom=R*R*R*sqrt(R)*(1-m)*(1-m) */
	return((
			(
				Rp3*(2 - m*(5 - 4*m)) -
				Rp2*m*(1 + m)*r1*r2 -
				R*(R*(3 - m)*m + (1 + m*(5 + 2*m))*r1p2)*r2p2 +
				2*((1 + (15 - 8*m)*m)*rp3 + 2*R*(3 - m)*(1 - 3*m)*r1)*r2p3 +
				R*(33 - m*(39 - 14*m))*r2p4 - 16*(3 - m)*r*r2p5
			)*dm -
			(
				Rp3*(2 - 4*m) +
				Rp2*m*(10*r1 - r2)*r2 -
				8*(1 - m)*rp3*r2p3 - 8*r*r2p5 -
				R*r2p2*((1 - m)*r*r - 4*r2p2 + 2*r*(3*r1 - 6*r2 + 4*m*r2))
			)*km
		)/(2*r2p3*denom)); 

}

/*****************************************************************************
Fourth order derivatives

	Fourth derivative with respect to z2, z2, z2, and z1 */
inline double dphidz2dz2dz2dz1(double dm, double km, double m, double z,
double r, double denom)  { 
  /* denom=R*R*R*R*sqrt(R)*(1-m)*(1-m)*(1-m) */
	double zp2=z*z;
	double zp4=zp2*zp2;

	double rp2=r*r;
	double rp4=rp2*rp2;
	double mm1=m-1;

	return((
			(
				6*(m-2)*mm1*mm1*rp4-
				6*mm1*(19+m*(6*m-17))*rp2*zp2 +
				(m*(18+m*(6*m-22))-50)*zp4
			)*dm + 
			(
				3*mm1*mm1*rp4-
				6*mm1*(3*m-7)*rp2*zp2+
				(26+m*(3*m-5))*zp4
			)*km
		)/denom);
}

/* Fourth derivative with respect to z2, z2, z2, and r1 */
inline double dphidz2dz2dz2dr1(double dm, double km, double m, double z,
double r1, double r2, double denom)  {
	if (r1 < eps)
		return(0.0);

	double r=r1+r2;
	double rp2=r*r;
	double zp2=z*z;

  /* denom=R*R*R*R*sqrt(R)*(1-m)*(1-m)*(1-m) */
	return((
			(6*(km - dm)*r2*(3*rp2 - 2*zp2))/m -
			(
				3*(m - 1)*(23 + m*(8*m - 23))*rp2*r +
				6*(m - 4)*(2*m - 1)*rp2*r2 +
				(107 + m*(-126 + (91 - 24*m)*m))*r*zp2 -
				4*(23 + (m - 3)*m)*r2*zp2
			)*dm +
			(-6*rp2*(2*(m - 2)*(m - 1)*r + (2 + m)*r2) +
				((47 + m*(12*m - 35))*r + 2*(m - 19)*r2)*zp2
			)*km
		)*z/denom);
}

/* Fourth derivative with respect to z2, z2, r2, and z1 */
inline double dphidz2dz2dr2dz1(double dm, double km, double m, double R,
double z, double r1, double r2, double denom)  {

	if (r2 < eps)
		return(0.0);

	double r=r1+r2;
	double rp2=r*r;
	double zp2=z*z;

  /* denom=R*R*R*R*sqrt(R)*(1-m)*(1-m)*(1-m) */
	if (r1 < eps)  {
		return((
				(6*R*(dm - km)*(3*rp2 - 2*zp2))/(4*r2) +
				(
					3*(m - 1)*(23 + m*(8*m - 23))*rp2*r +
					6*(m - 4)*(2*m - 1)*rp2*r1 +
					(107 + m*(-126 + (91 - 24*m)*m))*r*zp2 -
					4*(23 + (m - 3)*m)*r1*zp2
				)*dm -
				(-6*rp2*(2*(m - 2)*(m - 1)*r + (2 + m)*r1) +
					((47 + m*(12*m - 35))*r + 2*(m - 19)*r1)*zp2
				)*km
			)*z/denom);
	}
	return((
			(6*(dm - km)*r1*(3*rp2 - 2*zp2))/m +
			(
				3*(m - 1)*(23 + m*(8*m - 23))*rp2*r +
				6*(m - 4)*(2*m - 1)*rp2*r1 +
				(107 + m*(-126 + (91 - 24*m)*m))*r*zp2 -
				4*(23 + (m - 3)*m)*r1*zp2
			)*dm -
			(-6*rp2*(2*(m - 2)*(m - 1)*r + (2 + m)*r1) +
				((47 + m*(12*m - 35))*r + 2*(m - 19)*r1)*zp2
			)*km
		)*z/denom);
}

/* Fourth derivative with respect to z2, z2, r2, and r1 */
inline double dphidz2dz2dr2dr1(double dm, double km, double m, double R,
double z, double r1, double r2, double denom)  {

	if (r1 < eps)
		return(0.0);

	double r=r1+r2;
	double rp2=r*r;
	double zp2=z*z;

  /* denom=R*R*R*R*sqrt(R)*(1-m)*(1-m)*(1-m) */
	return((
			(6*(dm - km)*rp2*(rp2 - 4*zp2))/m -
			(
				rp2*(R*(15-m*(28-m*(27-8*m))) - 6*(10-m*(33 - 4*(7 - 2*m)*m))*zp2) +
				R*(R*(3 + (4 - m)*m*(1 - 2*m)) - (15 + m*(58 - m*(33 - 8*m)))*zp2)
			)*dm +
			(R*(R*(1 - m)*(3 + m) - (15 + (13 - 4*m)*m)*zp2) +
				rp2*(2*R*(6 - m*(5 - 2*m)) - 3*(15 - m*(21 - 8*m))*zp2)
			)*km
		)/denom); 
}

/* Fourth derivative with respect to z2, r2, r2, and z1 */
inline double dphidz2dr2dr2dz1(double dm, double km, double m, double R,
double z, double r1, double r2, double denom)  {
  /* denom=R*R*R*R*sqrt(R)*(1-m)*(1-m)*(1-m) */
	double Rp2=R*R;
	double Rp3=Rp2*R;
	double r=r1+r2;
	double rp2=r*r;
	double r2p2=r2*r2;
	double zp2=z*z;
	double mm1=1.0-m;

	return((
			(
				-2*Rp3*(1 + m*(m*m + m - 3)) +
				4*R*mm1*r2p2*(4*(m - 3)*(rp2 - 2*r2p2) + 5*mm1*(8*m - 11)*zp2) +
				8*r2p2*zp2*((47+m*(8*m-31))*(rp2-2*r2p2) + mm1*(41-12*m*(5-2*m))*zp2) +
				Rp2*(16*mm1*mm1*mm1*r2p2 + (6 - m*(15 + (19 - 4*m)*m))*zp2)
			)*dm + 
			(
				R*mm1*(Rp2*(2 - 3*m) - 8*(R*mm1 - rp2)*r2p2 - 16*r2p2*r2p2) +
				2*(Rp2*(-3 + m*(6 + m)) + 20*R*mm1*(3 - 2*m)*r2p2 + 16*(m - 3)*r2p2*(rp2 - 2*r2p2))*zp2 - 
				4*(47 + 3*m*(8*m - 21))*r2p2*zp2*zp2
			)*km
		)/(4*r2p2*denom)); 
}

/* Fourth derivative with respect to z2, r2, r2, and r1 */
inline double dphidz2dr2dr2dr1(double dm, double km, double m, double R,
double z, double r1, double r2, double denom)  {
	if (r1 < eps)
		return(0.0);

  /* denom=R*R*R*R*sqrt(R)*(1-m)*(1-m)*(1-m) */
	double r=r1+r2;
	double r1p2=r1*r1;
	double r2p2=r2*r2;
	double r2p4=r2p2*r2p2;
	double mm1=1.0-m;
	double temp=1.0+m;
	temp=temp*temp;


	return((
			(12*(dm - km)*(R*r*(r1 - r2) + 5*r2p4))/m +
			(
				R*(R*m*(34 + m*(m*(105 - m*(71 - 18*m)) - 69)) +
					2*(m*(m*(39 - 2*m*(19 - 6*m)) - 16) - 15)*r1p2) +
				R*temp*r1*r2 +
				2*R*(15 + m*(m*(225 - 2*m*(73 - 18*m)) - 160))*r2p2 +
				12*(m*(33 - 4*m*(7 - 2*m)) - 10)*r2p4
			)*dm +
			(
				R*(R*m*(m*(29 - 9*(3 - m)*m) - 19) + 3*(8 + (1 - m)*m*(5 - 4*m))*r1p2) -
				r1*r2*(R*(1 + m + 2*m*m) - 2*mm1*r1p2) +
				2*(3*R*(m*(21 - m*(19 - 6*m)) - 4) + mm1*r1p2)*r2p2 +
				2*mm1*r1*r2p2*r2 + 6*(15 - m*(21 - 8*m))*r2p4
			)*km 
		)*z/(2*r2*denom));
}

/* Fourth derivative with respect to r2, r2, r2, and z1 */
inline double dphidr2dr2dr2dz1(double dm, double km, double m, double R,
double z, double r1, double r2, double denom)  {
  /* denom=R*R*R*R*sqrt(R)*(1-m)*(1-m)*(1-m) */
	double r=r1+r2;
	double rp2=r*r;
	double rp3=rp2*r;
	double Rp2=R*R;
	double Rp3=Rp2*R;
	double mm1=1.0-m;
	double mp2=m*m;

	double r1p2=r1*r1;
	double r2p2=r2*r2;
	double r2p3=r2p2*r2;
	double r2p5=r2p3*r2p2;

	return((
			(
				2*Rp3*mm1*mm1*mm1 +
				8*rp3*r2p3*(3 + m*(35 - 2*m*(19 - 6*m))) -
				8*r*r2p5*(47 - m*(31 - 8*m)) +
				Rp2*r2*(2*mp2*m*(1 + m)*r1 - m*(m*(73 - 12*m*(4 - m)) - 1)*r2) +
				R*r2p2*(mm1*mm1*rp2 - 2*r1p2*(3 + 5*mp2) -
					2*(m*(154 - m*(109 - 28*m)) - 97)*r2p2 + 2*mm1*r*(3*r1 -6*r2+4*m*r2))
			)*dm +
			(
				Rp3*((5 - 6*m)*m - 2) +
				2*rp3*r2p3*(23 - m*(55 - 24*m)) +
				32*r*r2p5*(3 - m) +
				Rp2*m*r2*(6*r + r1 + 2*m*r1 + (4 + 3*(3 - m)*m)*r2) +
				R*r2p2*(2*(2 - m)*(1 - m)*rp2 + (1 + m + 2*mp2)*r1p2 -
					16*(3 - (3 - m)*m)*r*r2 + ((39 - 14*m)*m - 49)*r2p2)
			)*km
		)*z/(-2*r2p3*denom));
}

/* Fourth derivative with respect to r2, r2, r2, and r1 */
inline double dphidr2dr2dr2dr1(double dm, double km, double m, double R,
double r1, double r2, double denom)  {
	if (r1 < eps)
		return(0.0);

  /* denom=R*R*R*R*sqrt(R)*(1-m)*(1-m)*(1-m) */
	double Rp2=R*R;
	double mp2=m*m;
	double mp3=mp2*m;
	double temp=R*r1*r2;

	double r1p2=r1*r1;
	double r1p4=r1p2*r1p2;
	double r2p2=r2*r2;
	double r2p4=r2p2*r2p2;
	double r2p6=r2p4*r2p2;

	return((
			(km - dm)*(4*Rp2*r1p2 + 30*r2p6 - 2*R*r2p2*(r1p2 + 9*r2p2))/m -
			(2*(dm + km)*r1p4*r2p2)/m -
			(
				Rp2*(R*mp2*(1 - m*(13 - m*(17 - 3*(4 - m)*m))) -
					(3 - (3 - m)*m)*(4 - 3*mp3)*r1p2) -
				temp*(R*(1 - 2*mp3*(1 + m)) + mp2*(1 - 2*m)*r1p2) + 
				r2p2*(9*Rp2*m*(1 - (1 - m)*m)*(7 - m*(7 - 2*m)) -
					R*(1 - m*(4 + m*(11 - 8*m)))*r1p2 - 2*(2 + m)*r1p4) +
				R*(45 - m*(260 - m*(345 - 16*m*(13 - 3*m))))*r2p4 -
				6*r2p6*(10 - m*(33 - 4*m*(7 - 2*m)))
			)*dm +
			(
				Rp2*(R*mp2*(1 - (3 - m)*(2 - m)*m) -r1p2*(10 - m*(8 - m + 3*mp2 - mp3)))-
				temp*(R*(1 + m*(1 - 2*m*(2 - m - mp2))) + m*(3 + 2*mp2)*r1p2) +
				r2p2*(Rp2*m*(30 - m*(44 - 33*m + 9*mp2)) - R*m*(1 - 3*m)*r1p2 - (1 + m)*r1p4) -
				3*temp*r2p2 +
				r2p4*(R*(36 - m*(101 - m*(83 - 24*m))) + 2*(1 - 3*m)*r1p2) -
				3*r2p6*(15 - m*(21 - 8*m))
			)*km
		)/(r2p2*denom));
}

#endif /*!defined(derivs_h) */
