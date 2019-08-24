#include <math.h>
#include <iostream>
#include <assert.h>
#include "grid2d.H"

#define SWAP(x0,x) {float ** tmp=x0;x0=x;x=tmp;}
#define FOR_EACH_CELL for ( i=1 ; i<m-1 ; i++ ) { for ( j=1 ; j<n-1 ; j++ ) {
#define FOR_EACH_CELL1 for ( i=0 ; i<m-1 ; i++ ) { for ( j=0 ; j<n-1 ; j++ ) {
#define END_FOR }}
#define MAX(x,y) ( (x > y) ? (x) : (y) )
#define MIN(x,y) ( (x < y) ? (x) : (y) )
#define SIGN(x) ( (x > 0) ? (1) : ( (x < 0) ? (-1) : (0)) )
#define PI 3.14159f
#define ABS(x) ( (x > 0) ? (x) : (-x) )


int a = 0;
int MAX_CELL_PARTICLES = 16;

float cubic(float x){
	float ax = fabs(x);

	if (ax <= 1)
		return 1.5*ax*ax*ax-2.5*ax*ax+1;
	else if (ax <= 2)
		return -0.5*ax*ax*ax+2.5*ax*ax-4*ax+2;
	else
		return 0;
}

float quartic(float x){
	float ax = fabs(x);

	if (ax <= 1)
		return 4/3.0*ax*ax*ax-7/3.0*ax*ax+1;
	else if (ax <= 2)
		return -7/12.0*ax*ax*ax+3.0*ax*ax-59/12.0*ax+2.5;
	else if (ax <= 3)
		return 1/12.0*ax*ax*ax-2/3.0*ax*ax+7/4.0*ax-1.5;
	else
		return 0;
}

float linear(float x){
	float ax = fabs(x);

	if (ax <= 1)
		return 1-ax;
	else
		return 0;
}

Point& Point::operator = (const Point& pt)
{ x = pt.x; y = pt.y; return *this; }

//Point::~Point();
Cell::Cell(){}

Cell::Cell(int ii, int jj, float dx, float dy){
	i = ii;
	j = jj;
	phi = 1;
	phi_recon = 3.0f*dx;

	particles = new Point[MAX_CELL_PARTICLES];
	initial_particles = new bool[MAX_CELL_PARTICLES];
	new_particles = new Point[MAX_CELL_PARTICLES];
	initial_particles_new = new bool[MAX_CELL_PARTICLES];
	segments = new Point[4];
	segments1 = new Point[4];
	color = new float *[MAX_CELL_PARTICLES];
	new_color = new float *[MAX_CELL_PARTICLES];
	for (int k = 0; k < MAX_CELL_PARTICLES; k++){
		color[k] = new float [3];
		new_color[k] = new float [3];
	}
	num_segmentpts = 0;
	num_segmentpts1 = 0;
	num_particles = 0;
	num_newparticles = 0;
}

Cell::~Cell()
{
//	delete [] particles;
//	delete [] initial_particles;
//	delete [] new_particles;
//	delete [] segments;
//	delete [] segments1;
//	for (int i = 0; i < MAX_CELL_PARTICLES; i++){
//		delete [] color[i];
//		delete [] new_color[i];
//	}
	
}

//constructors
Grid::Grid() {};
Grid::Grid(int mm, int nn, float xx_lo, float xx_hi, float yy_lo, float yy_hi, int bbcxlo, int bbcxhi,int bbcylo, int bbcyhi){
	int i,j;
	m = mm;
	n = nn;
	trsize = 0;
	num_markers = 0;
	num_free_markers = 0;
	ssize = 0;
	ssize1 = 0;
	ssize2 = 0;
	mass = 0;
	x_lo = xx_lo;
	x_hi = xx_hi;
	y_lo = yy_lo;
	y_hi = yy_hi;
	dx = (x_hi-x_lo)/(m-1);
	dy = (y_hi-y_lo)/(n-1);
	bcxlo = bbcxlo;
	bcxhi = bbcxhi;
	bcylo = bbcylo;
	bcyhi = bbcyhi;
	iter = 0;

	phi = new float *[m];
	solidphi = new float *[m];
	den = new float *[m];	
	phi_recon = new float *[m];
	color = new float **[m];
	cell_list = new Cell *[m-1];
	for (i = 0; i < m-1; i++){
		cell_list[i] = new Cell [n-1];
		for (j = 0; j < n-1; j++)
			cell_list[i][j] = Cell(i,j,dx,dy);
	}

	temp_phi1 = new float *[m];
	temp_phi2 = new float *[m];
	u = new float *[m];
	v = new float *[m];
	F1 = new float *[m];
	F2 = new float *[m];	
	s = new Point [4*m*n];
	s1 = new Point [4*m*n];
	s2 = new Point [4*m*n];
	pointcolor = new float *[4*m*n];
	plist = new Point[9*MAX_CELL_PARTICLES];
	psize = 0;
	for (i = 0; i < m; i++){
	  phi[i] = new float [n];
	  solidphi[i] = new float [n];
	  den[i] = new float [n];
	  phi_recon[i] = new float [n];
	  color[i] = new float *[n];
	  temp_phi1[i] = new float [n];
	  temp_phi2[i] = new float [n];
	  u[i] = new float [n];
	  v[i] = new float [n];
	  F1[i] = new float [n];
	  F2[i] = new float [n];
	  for (j = 0; j < n; j++){
		  phi[i][j] = 3*dx;
		  solidphi[i][j] = 3*dx;
		  den[i][j] = 3*dx;
		  phi_recon[i][j] =3*dx;
		  color[i][j] = new float [3];
		  temp_phi1[i][j] = phi[i][j];
		  temp_phi2[i][j] = phi[i][j];
		  u[i][j] = 0.0f;
		  v[i][j] = 0.0f;
		  F1[i][j] = 0.0f;
		  F2[i][j] = 0.0f;
	  }
	}
	for (i = 0; i < 4*m*n; i++){
		s[i] = Point(0,0);
		s1[i] = Point(0,0);
		s2[i] = Point(0,0);
		pointcolor[i] = new float [3];
	}

	//the derivatives...
	dpx = new float *[m];
	dmx = new float *[m];
	dpy = new float *[m];
	dmy = new float *[m];
	for (i = 0; i < m; i++){
		dmx[i] = new float [n];
		dmy[i] = new float [n];
		dpx[i] = new float [n];
		dpy[i] = new float [n];
		for (j = 0; j < n; j++){
			dmx[i][j] = 0;
			dmy[i][j] = 0;
			dpx[i][j] = 0;
			dpy[i][j] = 0;
		}
	}
}

void Grid::InitializeProblem(char * problem_name){
	int i,j;
	double d1;
	float radius = 0.08f*MIN(m,n), dlt = 0.25f, left = 47.499f, right = 52.501f, notch = 85.001f, maxval = 0;
	float **temp_phi = new float *[m];
	for (i = 0; i < m; i++)
		temp_phi[i] = new float [n];


	if (problem_name == "zalesak"){	
		float xc = 50.00000f, yc = 75.00000f, r = 15.000000f, angle;
		//initialize notched circle
		for (i = 0; i < m; i++){
			for ( j = 0; j < n; j++){

				d1=sqrt(pow(x_lo+i*dx-xc,2)+pow(y_lo+j*dy-yc,2))-r;
				
				if ((x_lo+i*dx>=left)&&(x_lo+i*dx<=right)){
				  if (y_lo+j*dy<=60.0)
				   d1=( (x_lo+i*dx<xc) ? sqrt(pow(y_lo+j*dy-60.0,2)+pow(x_lo+i*dx-left,2)) : sqrt(pow(y_lo+j*dy-60.0,2)+pow(x_lo+i*dx-right,2)) );
				  else if (y_lo+j*dy<=notch)
				   d1=MIN( ( (x_lo+i*dx<50.0) ? (x_lo+i*dx-left) : (right-(x_lo+i*dx)) ) , notch-(y_lo+j*dy) );
				  else if ((y_lo+j*dy<=90.0)&&(d1<=0.0))
				   d1=MAX(d1,notch-(y_lo+j*dy));
				 }
				else if ((d1<0.0)&&(x_lo+i*dx<left)){
				  if (y_lo+j*dy<=notch)
				   d1=MAX(d1,(x_lo+i*dx-left));
				  else
				   d1=MAX(d1,-sqrt(pow(x_lo+i*dx-left,2)+pow(y_lo+j*dy-notch,2)) );
				 }
				else if ((d1<0.0)&&(x_lo+i*dx>right)){
				  if (y_lo+j*dy<=notch)
				   d1=MAX(d1,(right-(x_lo+i*dx)));
				  else
				   d1=MAX(d1,-sqrt(pow(x_lo+i*dx-right,2)+pow(y_lo+j*dy-notch,2)) );
				 }
				
				phi[i][j] = -float(d1);
				angle = atan2(x_lo+i*dx-xc, y_lo+j*dy-yc);
				if ((-PI < angle)&&(angle<=-PI/2)){
					color[i][j][0] = 1.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;}
				else if ((-PI/2 < angle) && (angle <= 0)){
					color[i][j][0] = 1.0f; color[i][j][1] = 0.0f; color[i][j][2] = 0.0f;}
				else if ((0 < angle) && (angle <= PI/2)){
					color[i][j][0] = 0.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;}
				else{ 
					color[i][j][0] = 0.0f; color[i][j][1] = 0.0f; color[i][j][2] = 1.0f;}
				color[i][j][0] = 0.0f; color[i][j][1] =0.0f; color[i][j][2] = 0.0f;

			}
		}

		//initialize velocity
		for (i = 0; i < m; i++)
		  for (j = 0; j < n; j++){
			  u[i][j] = PI/314.0f*(50-(y_lo+j*dy));
			  v[i][j] = PI/314.0f*(x_lo+i*dx-50);
		  }
		//initialize cells
		for ( i = 0; i < m-1; i++)
			for ( j = 0; j < n-1; j++){
				Cell* c = &cell_list[i][j];
				c->phi = (phi[i][j]+phi[i][j+1]+phi[i+1][j]+phi[i+1][j+1])/4;
			}
	}
	else if (problem_name == "vortex in a box"){
		float angle, xc = 0.5, yc = 0.75, r = 0.15;
		for (i = 0; i < m; i++)
			for ( j = 0; j < n; j++){
				phi[i][j] = -float((sqrt(pow(x_lo+i*dx-0.5,2)+pow(y_lo+j*dy-0.75,2))-0.15));
				angle = atan2(x_lo+i*dx-xc, y_lo+j*dy-yc);
				if ((-PI < angle)&&(angle<=-PI/2)){
					color[i][j][0] = 1.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;}
				else if ((-PI/2 < angle) && (angle <= 0)){
					color[i][j][0] = 1.0f; color[i][j][1] = 0.0f; color[i][j][2] = 0.0f;}
				else if ((0 < angle) && (angle <= PI/2)){
					color[i][j][0] = 0.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;}
				else{ 
					color[i][j][0] = 0.0f; color[i][j][1] = 0.0f; color[i][j][2] = 1.0f;}
				color[i][j][0] = 0.0f; color[i][j][1] = 0.0f; color[i][j][2] = 0.0f;
			}
		//initialize velocity
		for (i = 0; i < m; i++)
		  for (j = 0; j < n; j++){
			  v[i][j] = float(2*sin(PI*(x_lo+i*dx))*cos(PI*(x_lo+i*dx))*pow(sin(PI*(y_lo+j*dy)),2));
			  u[i][j] = float(-2*sin(PI*(y_lo+j*dy))*cos(PI*(y_lo+j*dy))*pow(sin(PI*(x_lo+i*dx)),2));
		  }
		//initialize cells
		for ( i = 0; i < m-1; i++)
			for ( j = 0; j < n-1; j++){
				Cell* c = &cell_list[i][j];
				c->phi = (phi[i][j]+phi[i][j+1]+phi[i+1][j]+phi[i+1][j+1])/4;
			}	
	}
	else if (problem_name == "contact"){
		float xc1 = 0.5, yc1 = 0.89, xc2 = 0.5, yc2 = 0.65, r1 = 0.1, r2 = 0.11, angle;
		for (j = 0; j < n; j++){
			if (j < int((yc2+r2-x_lo)/dy)+3){
				for ( i = 0; i < m; i++){
					phi[i][j] = -float((sqrt(pow(x_lo+i*dx-xc1,2)+pow(y_lo+j*dy-yc1,2))-r1));
					phi[i][j] = MAX(phi[i][j],-float((sqrt(pow(x_lo+i*dx-xc2,2)+pow(y_lo+j*dy-yc2,2))-r2)));
					angle = atan2(x_lo+i*dx-xc2, y_lo+j*dy-yc2);
					if ((-PI < angle)&&(angle<=-PI/2)){
						color[i][j][0] = 1.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;}
					else if ((-PI/2 < angle) && (angle <= 0)){
						color[i][j][0] = 1.0f; color[i][j][1] = 0.0f; color[i][j][2] = 0.0f;}
					else if ((0 < angle) && (angle <= PI/2)){
						color[i][j][0] = 0.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;}
					else{ 
						color[i][j][0] = 0.0f; color[i][j][1] = 0.0f; color[i][j][2] = 1.0f;}
					color[i][j][0] = 0.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;
				}
			}
			else{
				for ( i = 0; i < m; i++){
					phi[i][j] = -float((sqrt(pow(x_lo+i*dx-xc1,2)+pow(y_lo+j*dy-yc1,2))-r1));
					phi[i][j] = MAX(phi[i][j],-float((sqrt(pow(x_lo+i*dx-xc2,2)+pow(y_lo+j*dy-yc2,2))-r2)));
					angle = atan2(x_lo+i*dx-xc1, y_lo+j*dy-yc1);
					if ((-PI < angle)&&(angle<=-PI/2)){
						color[i][j][0] = 1.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;}
					else if ((-PI/2 < angle) && (angle <= 0)){
						color[i][j][0] = 1.0f; color[i][j][1] = 0.0f; color[i][j][2] = 0.0f;}
					else if ((0 < angle) && (angle <= PI/2)){
						color[i][j][0] = 0.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;}
					else{ 
						color[i][j][0] = 0.0f; color[i][j][1] = 0.0f; color[i][j][2] = 1.0f;}
					color[i][j][0] = 0.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;
				}
			
			}
			
		}
		//initialize velocity
		for (i = 0; i < m; i++)
		  for (j = 0; j < n; j++){
			  v[i][j] = 0;
			  u[i][j] = 0;
		  }
		//initialize cells
		for ( i = 0; i < m-1; i++)
			for ( j = 0; j < n-1; j++){
				Cell* c = &cell_list[i][j];
				c->phi = (phi[i][j]+phi[i][j+1]+phi[i+1][j]+phi[i+1][j+1])/4;
			}	
	}
	else if (problem_name == "contact1"){
		float xc1 = 0.5, yc1 = 0.8, xc2 = 0.5, yc2 = 0.2, r1 = 0.15, r2 = 0.15, angle;
		for (j = 0; j < n; j++){
			if (j < int((yc2+r2-x_lo)/dy)+3){
				for ( i = 0; i < m; i++){
					phi[i][j] = -float((sqrt(pow(x_lo+i*dx-xc1,2)+pow(y_lo+j*dy-yc1,2))-r1));
					phi[i][j] = MAX(phi[i][j],-float((sqrt(pow(x_lo+i*dx-xc2,2)+pow(y_lo+j*dy-yc2,2))-r2)));
					angle = atan2(x_lo+i*dx-xc2, y_lo+j*dy-yc2);
					if ((-PI < angle)&&(angle<=-PI/2)){
						color[i][j][0] = 1.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;}
					else if ((-PI/2 < angle) && (angle <= 0)){
						color[i][j][0] = 1.0f; color[i][j][1] = 0.0f; color[i][j][2] = 0.0f;}
					else if ((0 < angle) && (angle <= PI/2)){
						color[i][j][0] = 0.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;}
					else{ 
						color[i][j][0] = 0.0f; color[i][j][1] = 0.0f; color[i][j][2] = 1.0f;}
					color[i][j][0] = 0.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;
				}
			}
			else{
				for ( i = 0; i < m; i++){
					phi[i][j] = -float((sqrt(pow(x_lo+i*dx-xc1,2)+pow(y_lo+j*dy-yc1,2))-r1));
					phi[i][j] = MAX(phi[i][j],-float((sqrt(pow(x_lo+i*dx-xc2,2)+pow(y_lo+j*dy-yc2,2))-r2)));
					angle = atan2(x_lo+i*dx-xc1, y_lo+j*dy-yc1);
					if ((-PI < angle)&&(angle<=-PI/2)){
						color[i][j][0] = 1.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;}
					else if ((-PI/2 < angle) && (angle <= 0)){
						color[i][j][0] = 1.0f; color[i][j][1] = 0.0f; color[i][j][2] = 0.0f;}
					else if ((0 < angle) && (angle <= PI/2)){
						color[i][j][0] = 0.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;}
					else{ 
						color[i][j][0] = 0.0f; color[i][j][1] = 0.0f; color[i][j][2] = 1.0f;}
				//	color[i][j][0] = 0.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;
				}
			
			}
			
		}
		//initialize velocity
		for (i = 0; i < m; i++)
		  for (j = 0; j < n; j++){
			  u[i][j] = 0;
			  if (j < n/2)
				  v[i][j] = 1;
			  else
				  v[i][j] = -1;
			  
		  }
		//initialize cells
		for ( i = 0; i < m-1; i++)
			for ( j = 0; j < n-1; j++){
				Cell* c = &cell_list[i][j];
				c->phi = (phi[i][j]+phi[i][j+1]+phi[i+1][j]+phi[i+1][j+1])/4;
			}	
	}
	else if (problem_name == "deformation field"){
		float xc = 0.50000, yc = 0.5, r = 0.15, angle;
		for (i = 0; i < m; i++)
			for ( j = 0; j < n; j++){
				phi[i][j] = -float((sqrt(pow(x_lo+i*dx-xc,2)+pow(y_lo+j*dy-yc,2))-r));
				angle = atan2(x_lo+i*dx-xc, y_lo+j*dy-yc);
				if ((-PI < angle)&&(angle<=-PI/2)){
					color[i][j][0] = 1.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;}
				else if ((-PI/2 < angle) && (angle <= 0)){
					color[i][j][0] = 1.0f; color[i][j][1] = 0.0f; color[i][j][2] = 0.0f;}
				else if ((0 < angle) && (angle <= PI/2)){
					color[i][j][0] = 0.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;}
				else{ 
					color[i][j][0] = 0.0f; color[i][j][1] = 0.0f; color[i][j][2] = 1.0f;}
				color[i][j][0] = 0.0f; color[i][j][1] = 0.0f; color[i][j][2] = 0.0f;
			}
		//initialize velocity
		for (i = 0; i < m; i++)
		  for (j = 0; j < n; j++){
			  v[i][j] = float(cos(4*PI*(x_lo+i*dx+0.5))*cos(4*PI*(y_lo+j*dy+0.5)));
			  u[i][j] = float(sin(4*PI*(x_lo+i*dx+0.5))*sin(4*PI*(y_lo+j*dy+0.5)));
			    maxval = MAX(maxval,fabs(v[i][j]));
			  maxval = MAX(maxval,fabs(u[i][j]));
		  }
		//initialize cells
		for ( i = 0; i < m-1; i++)
			for ( j = 0; j < n-1; j++){
				Cell* c = &cell_list[i][j];
				c->phi = (phi[i][j]+phi[i][j+1]+phi[i+1][j]+phi[i+1][j+1])/4;
			}	
	}
	else if (problem_name == "rotation"){
		float xc = 0.5, yc = 0.5, r = 0.3, angle;
		for (i = 0; i < m; i++)
			for ( j = 0; j < n; j++){
				phi[i][j] = -float((sqrt(pow(x_lo+i*dx-xc,2)+pow(y_lo+j*dy-yc,2))-r));
				angle = atan2(x_lo+i*dx-xc, y_lo+j*dy-yc);
				if ((-PI < angle)&&(angle<=-PI/2)){
					color[i][j][0] = 1.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;}
				else if ((-PI/2 < angle) && (angle <= 0)){
					color[i][j][0] = 1.0f; color[i][j][1] = 0.0f; color[i][j][2] = 0.0f;}
				else if ((0 < angle) && (angle <= PI/2)){
					color[i][j][0] = 0.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;}
				else{ 
					color[i][j][0] = 0.0f; color[i][j][1] = 0.0f; color[i][j][2] = 1.0f;}
				//color[i][j][0] = 0.0f; color[i][j][1] = 0.0f; color[i][j][2] = 0.0f;
			}
		//initialize velocity
		for (i = 0; i < m; i++)
		  for (j = 0; j < n; j++){
			  v[i][j] = float(cos(PI*(x_lo+i*dx+0.0))*cos(PI*(y_lo+j*dy+0.5)));
			  u[i][j] = float(sin(PI*(x_lo+i*dx+0.0))*sin(PI*(y_lo+j*dy+0.5)));
			  //u[i][j] = 0.0f;
			  //v[i][j] = 0.05f;
		  }
		for (i = 0; i < m; i++)
		  for (j = 0; j < n; j++){
			  u[i][j] = PI*(0.5-(y_lo+j*dy));
			  v[i][j] = PI*(x_lo+i*dx-0.5);
		  }
		//initialize cells
		for ( i = 0; i < m-1; i++)
			for ( j = 0; j < n-1; j++){
				Cell* c = &cell_list[i][j];
				c->phi = (phi[i][j]+phi[i][j+1]+phi[i+1][j]+phi[i+1][j+1])/4;
			}	
	}
	else if (problem_name == "ellipse"){
		float xc = 3.5000001, yc = 2.00001, A = 4.0, B = 2.0;
		for (i = 0; i < m; i++)
			for ( j = 0; j < n; j++){
				phi[i][j] = -float((sqrt(pow(x_lo+i*dx,2)/A/A+pow(y_lo+j*dy,2)/B/B)-1.0)*
					(0.1+pow(x_lo+i*dx-xc,2)+pow(y_lo+j*dy-yc,2)));
				color[i][j][0] = 0.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;
				color[i][j][0] = 0.0f; color[i][j][1] = 0.0f; color[i][j][2] = 0.0f;
			}
		//initialize velocity
		for (i = 0; i < m; i++)
		  for (j = 0; j < n; j++){
			  v[i][j] = float(cos(PI*(x_lo+i*dx+0.0))*cos(PI*(y_lo+j*dy+0.5)));
			  u[i][j] = float(sin(PI*(x_lo+i*dx+0.0))*sin(PI*(y_lo+j*dy+0.5)));
		  }
		//initialize cells
		for ( i = 0; i < m-1; i++)
			for ( j = 0; j < n-1; j++){
				Cell* c = &cell_list[i][j];
				c->phi = (phi[i][j]+phi[i][j+1]+phi[i+1][j]+phi[i+1][j+1])/4;
			}	
	}
	else if (problem_name == "line"){
		float xc = 0.000001;
		for (i = 0; i < m; i++)
			for ( j = 0; j < n; j++){
				phi[i][j] = -4.0*float(fabs(x_lo+i*dx-xc)-3);
				color[i][j][0] = 0.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;
			}
		//initialize velocity
		for (i = 0; i < m; i++)
		  for (j = 0; j < n; j++){
			  v[i][j] = 0;
			  u[i][j] = 1;
		  }
		//initialize cells
		for ( i = 0; i < m-1; i++)
			for ( j = 0; j < n-1; j++){
				Cell* c = &cell_list[i][j];
				c->phi = (phi[i][j]+phi[i][j+1]+phi[i+1][j]+phi[i+1][j+1])/4;
			}	
	}
	else if (problem_name == "square"){
		float side = 1.5f, r = 2, angle;
		float xc = 5.55f, yc = 2.55f;
		//initialize square
		for (i = 0; i < m; i++)
			for ( j = 0; j < n; j++){
				phi[i][j] = -float((sqrt(pow(x_lo+i*dx-xc,2)+pow(y_lo+j*dy-yc,2))-r));
				angle = atan2(x_lo+i*dx-xc, y_lo+j*dy-yc);
				if ((-PI < angle)&&(angle<=-PI/2)){
					color[i][j][0] = 1.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;}
				else if ((-PI/2 < angle) && (angle <= 0)){
					color[i][j][0] = 1.0f; color[i][j][1] = 0.0f; color[i][j][2] = 0.0f;}
				else if ((0 < angle) && (angle <= PI/2)){
					color[i][j][0] = 0.0f; color[i][j][1] = 1.0f; color[i][j][2] = 0.0f;}
				else{ 
					color[i][j][0] = 0.0f; color[i][j][1] = 0.0f; color[i][j][2] = 1.0f;}
			}
		//initialize velocity
		for (i = 0; i < m; i++)
		  for (j = 0; j < n; j++){
			  u[i][j] = 0.0f;
			  v[i][j] = 0.05f;
		  }
		//initialize cells
		for ( i = 0; i < m-1; i++)
			for ( j = 0; j < n-1; j++){
				Cell* c = &cell_list[i][j];
				c->phi = (phi[i][j]+phi[i][j+1]+phi[i+1][j]+phi[i+1][j+1])/4;
			}
	}
	else if (problem_name == "square1"){
		float side = 5.999f;
		float xc = 5.00001f, yc = 5.4999f;
		//initialize square
		for (i = 0; i < m; i++){
			for ( j = 0; j < n; j++){
				if ((fabs(x_lo+i*dx-xc) <= side)&&(fabs(j-yc) <= side))
					phi[i][j] = float(MIN(fabs(side-fabs(i-xc)),fabs(side-fabs(j-yc))));
				else if ((fabs(i-xc) <= side)&&(fabs(j-yc) > side))
					phi[i][j] = -float(fabs(side-fabs(j-yc)));
				else if ((fabs(i-xc) > side)&&(fabs(j-yc) <= side))
					phi[i][j] = -float(fabs(side-fabs(i-xc)));
				else{
					float a = float(sqrt(pow(i-xc+side,2)+pow(j-yc+side,2)));
					float b = float(sqrt(pow(i-xc-side,2)+pow(j-yc-side,2)));
					float c = float(sqrt(pow(i-xc+side,2)+pow(j-yc-side,2)));
					float d = float(sqrt(pow(i-xc-side,2)+pow(j-yc+side,2)));
					phi[i][j] = -MIN(a,MIN(b,MIN(c,d)));
				}
				phi[i][j] *= -1;
			}
		}
		//initialize velocity
		for (i = 0; i < m; i++)
		  for (j = 0; j < n; j++){
			  u[i][j] = PI/314.0f*(n/2-j);
			  v[i][j] = PI/314.0f*(i-m/2);
			  u[i][j] = 0;
			  v[i][j] = 0.05;
		  }
		//initialize cells
		for ( i = 0; i < m-1; i++)
			for ( j = 0; j < n-1; j++){
				Cell* c = &cell_list[i][j];
				c->phi = (phi[i][j]+phi[i][j+1]+phi[i+1][j]+phi[i+1][j+1])/4;
			}
	}
	else if (problem_name == "squarec"){
		float side = 0.999f;
		float xc = m/float(2)-0.4999f, yc = 4.4999f;
		//initialize square
		for (i = 0; i < m-1; i++){
			for ( j = 0; j < n-1; j++){
				Cell* cl = &cell_list[i][j];
				if ((fabs(i+0.5-xc) <= side)&&(fabs(j+0.5-yc) <= side))
					cl->phi = float(MIN(fabs(side-fabs(i+0.5-xc)),fabs(side-fabs(j+0.5-yc))));
				else if ((fabs(i+0.5-xc) <= side)&&(fabs(j+0.5-yc) > side))
					cl->phi = -float(fabs(side-fabs(j+0.5-yc)));
				else if ((fabs(i+0.5-xc) > side)&&(fabs(j+0.5-yc) <= side))
					cl->phi = -float(fabs(side-fabs(i+0.5-xc)));
				else{
					float a = float(sqrt(pow(i+0.5-xc+side,2)+pow(j+0.5-yc+side,2)));
					float b = float(sqrt(pow(i+0.5-xc-side,2)+pow(j+0.5-yc-side,2)));
					float c = float(sqrt(pow(i+0.5-xc+side,2)+pow(j+0.5-yc-side,2)));
					float d = float(sqrt(pow(i+0.5-xc-side,2)+pow(j+0.5-yc+side,2)));
					cl->phi = -MIN(a,MIN(b,MIN(c,d)));
				}
				cl->phi *= -1;
			}
		}
		//initialize velocity
		for (i = 0; i < m; i++)
		  for (j = 0; j < n; j++){
			  u[i][j] = PI/314.0f*(n/2-j);
			  v[i][j] = PI/314.0f*(i-m/2);
		  }
	}
		else if (problem_name == "diamond"){
		float side = 0.999f;
		float xc = m/float(2)-0.4999f, yc = 4.4999f;
		//initialize square
		for (i = 0; i < m-1; i++){
			for ( j = 0; j < n-1; j++){
				Cell* cl = &cell_list[i][j];
				if ((fabs(i+0.5-xc) <= side)&&(fabs(j+0.5-yc) <= side))
					cl->phi = float(MIN(fabs(side-fabs(i+0.5-xc)),fabs(side-fabs(j+0.5-yc))));
				else if ((fabs(i+0.5-xc) <= side)&&(fabs(j+0.5-yc) > side))
					cl->phi = -float(fabs(side-fabs(j+0.5-yc)));
				else if ((fabs(i+0.5-xc) > side)&&(fabs(j+0.5-yc) <= side))
					cl->phi = -float(fabs(side-fabs(i+0.5-xc)));
				else{
					float a = float(sqrt(pow(i+0.5-xc+side,2)+pow(j+0.5-yc+side,2)));
					float b = float(sqrt(pow(i+0.5-xc-side,2)+pow(j+0.5-yc-side,2)));
					float c = float(sqrt(pow(i+0.5-xc+side,2)+pow(j+0.5-yc-side,2)));
					float d = float(sqrt(pow(i+0.5-xc-side,2)+pow(j+0.5-yc+side,2)));
					cl->phi = -MIN(a,MIN(b,MIN(c,d)));
				}
				cl->phi *= -1;
			}
		}
		//initialize velocity
		for (i = 0; i < m; i++)
		  for (j = 0; j < n; j++){
			  u[i][j] = PI/314.0f*(n/2-j);
			  v[i][j] = PI/314.0f*(i-m/2);
		  }
	}
	else{		
		for (i=0 ; i<m ; i++) 
			for (j = 0; j < n; j++)
				phi[i][j] = -float((sqrt(pow(i-(m-1)/2.0f,2)/16+pow(j-(n-1)/2.0f,2))-radius));

		for (i=0 ; i<m-1 ; i++) 
			for (j = 0; j < n-1; j++){
				Cell* cl = &cell_list[i][j];
				cl->phi = (phi[i][j]+phi[i][j+1]+phi[i+1][j]+phi[i+1][j+1])/4;
			}
	}

	for (i = 0; i < m; i++)
		delete [] temp_phi[i];
	delete [] temp_phi;
}

void Grid::InitializeParticles(void){
	int i,j,k,l,np,index,s;
	Cell* c;
	Point a,b,p;
	double d;
	float t,dist,dist0;

	for (i = 0; i < m-1; i++)
		for (j = 0; j < n-1; j++){
			c = &cell_list[i][j];
			for (s = 0; s < MAX_CELL_PARTICLES; s++)
				c->initial_particles_new[s] = 0;
			if (CellChangesSignInside(c, 0)){
				index = 0;
				for (k = 0; k < c->num_segmentpts/2; k++){//for each segment
					a = c->segments[2*k];
					b = c->segments[2*k+1];
					d = sqrt(pow(a.x-b.x,2)+pow(a.y-b.y,2));
					np = int(floor(d/0.1f));//10 particles per unit length
					if (np < 2)
						np = 2;
					for (l = 0; l < np; l++){
						t = l/float(np);
						p = Point(x_lo+dx*(t*a.x+(1-t)*b.x),y_lo+dy*(t*a.y+(1-t)*b.y));
						c->particles[c->num_particles++] = p;
						c->initial_particles[index++] = true;
						dist0 = 3*dx;
						for (int ix = i; ix <= i+1; ix++)
							for (int iy = j; iy <= j+1; iy++){
								dist = float(sqrt((p.x-ix*dx-x_lo)*(p.x-ix*dx-x_lo)+(p.y-iy*dy-y_lo)*(p.y-iy*dy-y_lo)));
								if (dist < dist0){
									dist0 = dist;
									c->color[c->num_particles-1] = color[ix][iy];
								}
							}
						//c->color[c->num_particles-1] = color[i][j];
					}
				}
			}	
		}
	psize = CountParticles();
}

bool Grid::CellChangesSignInside(Cell* c, float threshold){
	int i = c->i, j = c->j;
	int t2,t3,t4,p=SIGN(phi[i][j]-threshold);
		
	t2 = SIGN(phi[i][j+1]-threshold);
	t3 = SIGN(phi[i+1][j]-threshold);
	t4 = SIGN(phi[i+1][j+1]-threshold);

	if ((p*t2<=0)||(p*t3<=0)||(p*t4<=0))
		return true;
	else
		return false;
}

void Grid::DPlusx (float **dplus, float **u1){
	int i,j;
	for (i = 1; i < m-1; i++)
		for (j = 1; j < n-1; j++)
			dplus[i][j] = (u1[i+1][j]-u1[i][j])/dx;
	set_bnd ( 0, dplus );
}

void Grid::DMinusx (float** dminus, float **u1){
	int i,j;
	for (i = 1; i < m-1; i++)
		for (j = 1; j < n-1; j++)
			dminus[i][j] = (u1[i][j]-u1[i-1][j])/dx;
	set_bnd ( 0, dminus );
}

void Grid::DPlusy (float** dplus, float **v1){
	int i,j;
	for (i = 1; i < m-1; i++)
		for (j = 1; j < n-1; j++)
			dplus[i][j] = (v1[i][j+1]-v1[i][j])/dy;
	set_bnd ( 0, dplus );
}

void Grid::DMinusy (float** dminus, float **v1){
	int i,j;
	for (i = 1; i < m-1; i++)
		for (j = 1; j < n-1; j++)
			dminus[i][j] = (v1[i][j]-v1[i][j-1])/dy;	
	set_bnd ( 0, dminus );
}

void Grid::Reinitialize(int iter, float dtau){
	int i,j,ii;
	float fact, epsilon = dx, grad;
	
	for (ii = 0; ii < iter; ii++){
		DPlusx(dpx,phi);
		DMinusx(dmx,phi);
		DPlusy(dpy,phi);
		DMinusy(dmy,phi);
		for (i=1 ; i<m-1 ; i++) 
			for (j = 1; j < n-1; j++){
				if (phi[i][j] > 0)
					fact = float(sqrt(MAX(pow(MAX(dmx[i][j],0),2),pow(MIN(dpx[i][j],0),2))+
						MAX(pow(MAX(dmy[i][j],0),2),pow(MIN(dpy[i][j],0),2)))-1);
				else 
					fact = float(sqrt(MAX(pow(MAX(dpx[i][j],0),2),pow(MIN(dmx[i][j],0),2))+
						MAX(pow(MAX(dpy[i][j],0),2),pow(MIN(dmy[i][j],0),2)))-1);

				phi[i][j] -= dtau*SIGN(phi[i][j])*fact;
				grad = float(sqrt(pow((dpx[i][j]+dmx[i][j])/2,2)+pow((dpy[i][j]+dmy[i][j])/2,2)));
		//		phi[i][j] -= float(dtau*phi[i][j]/sqrt(pow(phi[i][j],2)+
		//			pow(grad,2)*pow(epsilon,2))*fact);
		//		phi[i][j] -= float(dtau*phi[i][j]/sqrt(pow(phi[i][j],2)+
		//			pow(epsilon,2))*fact);
				//the pow(1-grad,2) or pow(grad,2) stuff is the peng/merriman... fix
				//doesn't seem to do much for the zalesak's problem though.... well, for
				//dt = 0.09 it fixes something though 1-grad actually blows things up....???
			}
	}
}

void Grid::ReinitializeCIR(int iter, float dtau){
	int i,j,ii,ix,iy;
	float epsilon = dx/6, grad, gradx, grady, ldx, ldy;
	Point p;
	for (ii = 0; ii < iter; ii++){
		DPlusx(dpx,phi);
		DMinusx(dmx,phi);
		DPlusy(dpy,phi);
		DMinusy(dmy,phi);
		for (i=1 ; i<m-1 ; i++) 
			for (j = 1; j < n-1; j++){
				if (cell_list[i][j].num_particles > 0)
					int aa = 0;
			//	gradx = (dpx[i][j]+dmx[i][j])/2;
			//	grady = (dpy[i][j]+dmy[i][j])/2;
				if (phi[i][j] > 0){
					gradx = dmx[i][j];
					grady =	dmy[i][j];
				}
				else{ 
					gradx = dpx[i][j];
					grady = dpy[i][j];
				}
				grad = float(sqrt(pow(gradx,2)+pow(grady,2)));
				p.x = x_lo+i*dx-dtau*gradx/grad*phi[i][j]/sqrt(pow(phi[i][j],2)+pow(epsilon,2));
				p.y = y_lo+j*dy-dtau*grady/grad*phi[i][j]/sqrt(pow(phi[i][j],2)+pow(epsilon,2));
				ix = (int)((p.x-x_lo)/dx); if (ix < 1) ix = 1; if (ix > m-3) ix = m-3;
				iy = (int)((p.y-y_lo)/dy); if (iy < 1) iy = 1; if (iy > n-3) iy = n-3;
				ldx = p.x-x_lo-ix*dx; ldx /= dx;
				ldy = p.y-y_lo-iy*dy; ldy /= dy;
				phi_recon[i][j] = Interpolate(phi,1,ix,iy,ldx,ldy) + dtau*phi[i][j]/sqrt(pow(phi[i][j],2)+pow(epsilon,2));
			}
		for (i=1 ; i<m-1 ; i++) 
			for (j = 1; j < n-1; j++){
				phi[i][j] = phi_recon[i][j];
				phi_recon[i][j] = 3*dx;
			}
	}
}

void Grid::ReinitializeRS(int iter, float dtau){
	int i,j,ii;
	float fact=0, grad, den, den1, eps = 1e-2f, epsilon = dx, d;
	Cell *c;
	Point p;

	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++){
			phi_recon[i][j] = phi[i][j];
			temp_phi1[i][j] = phi[i][j];
			temp_phi2[i][j] = phi[i][j];
		}
	for (i = 1; i < m-1; i++)
		for (j = 1; j < n-1; j++){
			p = Point(i,j);
			
			den = sqrt(pow(phi[i+1][j]-phi[i-1][j],2)+pow(phi[i][j+1]-phi[i][j-1],2))/dx/2;
			if (phi[i][j] > 0)
				den = float(sqrt(MAX(pow(MAX((phi[i][j]-phi[i-1][j])/dx,0),2),pow(MIN((phi[i+1][j]-phi[i][j])/dx,0),2))+
					MAX(pow(MAX((phi[i][j]-phi[i][j-1])/dy,0),2),pow(MIN((phi[i][j+1]-phi[i][j])/dy,0),2))));
			else 
				den = float(sqrt(MAX(pow(MAX((phi[i+1][j]-phi[i][j])/dx,0),2),pow(MIN((phi[i][j]-phi[i-1][j])/dx,0),2))+
					MAX(pow(MAX((phi[i][j+1]-phi[i][j])/dy,0),2),pow(MIN((phi[i][j]-phi[i][j-1])/dy,0),2))));
			if (den < eps){//for robustness....
				den = MAX(den,MAX((phi[i+1][j]-phi[i][j])/dx,MAX((phi[i][j]-phi[i-1][j])/dx,MAX((phi[i][j+1]-phi[i][j])/dy,MAX((phi[i][j]-phi[i][j-1])/dy,eps)))));
			}
			if ((phi[i][j]*phi[i-1][j] <= 0)||(phi[i][j]*phi[i+1][j] <= 0)||
					(phi[i][j]*phi[i][j-1] <= 0)||(phi[i][j]*phi[i][j+1] <= 0)){
				den1 = 3*dx;
				for (int h = -1; h <= 0; h++)
					for (int l = -1; l <= 0; l++){
						c = &cell_list[i+h][j+l];
						den1 = MIN(den1,FindLineDist(c,p,0));//exact distance to interface...
					}
				temp_phi2[i][j] = den1*SIGN(phi[i][j]);//phi[i][j]/sqrt(pow(phi[i][j],2)+pow(epsilon,2));//SIGN(phi[i][j]);
			//	temp_phi2[i][j] = SIGN(phi[i][j]);
			//	temp_phi2[i][j] = den1*phi[i][j]/sqrt(pow(phi[i][j],2)+pow(epsilon,2));//SIGN(phi[i][j]);
			}
			temp_phi1[i][j] = phi[i][j]/den;			
		}
	
		

	for (ii = 0; ii < iter; ii++){
		for (i=1 ; i<m-1 ; i++){ 
			for (j = 1; j < n-1; j++){
				if ((phi[i][j]*phi[i-1][j] <= 0)||(phi[i][j]*phi[i+1][j] <= 0)||
					(phi[i][j]*phi[i][j-1] <= 0)||(phi[i][j]*phi[i][j+1] <= 0)){
					//we are inside the narrow band				
					fact = dtau/dx*(SIGN((phi[i][j]))*fabs(phi_recon[i][j])-temp_phi1[i][j]);
				//	fact = dtau/dx*(phi[i][j]/sqrt(pow(phi[i][j],2)+pow(epsilon,2))*fabs(phi_recon[i][j])-temp_phi2[i][j]);
				//	fact = dtau/dx*(phi[i][j]/sqrt(pow(phi[i][j],2)+den*den*pow(epsilon,2))*fabs(phi_recon[i][j])-temp_phi1[i][j]);
				
				}
				else{
					//outside the narrow band
					if (phi[i][j] > 0)
						grad = float(sqrt(MAX(pow(MAX((phi_recon[i][j]-phi_recon[i-1][j])/dx,0),2),pow(MIN((phi_recon[i+1][j]-phi_recon[i][j])/dx,0),2))+
						MAX(pow(MAX((phi_recon[i][j]-phi_recon[i][j-1])/dy,0),2),pow(MIN((phi_recon[i][j+1]-phi_recon[i][j])/dy,0),2))))-1;
					else 
						grad = float(sqrt(MAX(pow(MAX((phi_recon[i+1][j]-phi_recon[i][j])/dx,0),2),pow(MIN((phi_recon[i][j]-phi_recon[i-1][j])/dx,0),2))+
						MAX(pow(MAX((phi_recon[i][j+1]-phi_recon[i][j])/dy,0),2),pow(MIN((phi_recon[i][j]-phi_recon[i][j-1])/dy,0),2))))-1;
					fact = dtau*grad*phi[i][j]/sqrt(pow(phi[i][j],2)+pow(epsilon,2));
				}
				
				phi_recon[i][j] -= fact;
			}
		}
	}
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++){
			phi[i][j] = phi_recon[i][j];
			phi_recon[i][j] = 3*dx;
			temp_phi1[i][j] = 3*dx;
			temp_phi2[i][j] = 3*dx;
		}
}

void Grid::ReinitializeV(int iter, float dtau){
	int i,j,ii,h,l,k;
	float fact, grad, den, den1, eps = 1e-1f, epsilon = dx, ldx, ldy, ldxi, ldyj, dist,q,c=0.5,rho=2.0;
	Point p;
	Cell *d;

	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++){
			phi_recon[i][j] = phi[i][j];
			temp_phi2[i][j] = 0;
		}
	//work inside the band	
	for (i = 1; i < m-1; i++)
		for (j = 1; j < n-1; j++){
			float lambda = 0, wght = 0, wght1 = 0,twght = 0;  
			int lim = 1;
			//use the cells about (i,j) to find correction for phi
			if ((phi[i][j]*phi[i-1][j] <= 0)||(phi[i][j]*phi[i+1][j] <= 0)||
				(phi[i][j]*phi[i][j-1] <= 0)||(phi[i][j]*phi[i][j+1] <= 0)){
					for (h = -lim; h <=lim-1; h++){
					if ((0<i+h)&&(i+h<m-1)){
						for (l = -lim; l <=lim-1; l++){
							if ((0<j+l)&&(j+l<n-1)){
								d = &cell_list[i+h][j+l];
								for (k=0; k < d->num_particles; k++){
									 p = d->particles[k];
									 ldx = p.x-x_lo-(i+h)*dx; ldx /= dx;
									 ldy = p.y-y_lo-(j+l)*dy; ldy /= dy;
									 ldxi = p.x-x_lo-i*dx; ldxi /= dx;
									 ldyj = p.y-y_lo-j*dy; ldyj /= dy;
									 wght = 1/(1+pow(ldxi*ldxi+ldyj*ldyj,3));
									 q = sqrt(ldxi*ldxi+ldyj*ldyj)/rho;
									 if (q<1)
										 wght1= (exp(-q*q/(c*c))-exp(-1/(c*c)))/(1-exp(-1/(c*c)));
									 else wght1 = 0;
									 twght += wght1;
									 lambda += wght1*Interpolate(phi,1,i+h,j+l,ldx,ldy);
								}
							}
						}
					}
				}
				if (twght > 0){								
					phi_recon[i][j] = phi[i][j]-lambda/twght;
					temp_phi2[i][j] = 1;
				}
				if (phi[i][j] > 0)
					den = float(sqrt(MAX(pow(MAX((phi[i][j]-phi[i-1][j])/dx,0),2),pow(MIN((phi[i+1][j]-phi[i][j])/dx,0),2))+
						MAX(pow(MAX((phi[i][j]-phi[i][j-1])/dy,0),2),pow(MIN((phi[i][j+1]-phi[i][j])/dy,0),2))));
				else 
					den = float(sqrt(MAX(pow(MAX((phi[i+1][j]-phi[i][j])/dx,0),2),pow(MIN((phi[i][j]-phi[i-1][j])/dx,0),2))+
						MAX(pow(MAX((phi[i][j+1]-phi[i][j])/dy,0),2),pow(MIN((phi[i][j]-phi[i][j-1])/dy,0),2))));								
				phi_recon[i][j] /= den;
	
			}						
		}
	for (i=0; i < m; i++)
		for (j=0; j<n; j++)
				phi[i][j] = phi_recon[i][j];
	
	//work only outside the band
	for (ii = 0; ii < iter; ii++){
		for (i=1 ; i<m-1 ; i++){ 
			for (j = 1; j < n-1; j++){
				
				if ((phi[i][j]*phi[i-1][j] <= 0)||(phi[i][j]*phi[i+1][j] <= 0)||
					(phi[i][j]*phi[i][j-1] <= 0)||(phi[i][j]*phi[i][j+1] <= 0)){
					fact = 0;
				}
				else{
					//outside the narrow band
					if (phi[i][j] > 0)
						grad = float(sqrt(MAX(pow(MAX((phi_recon[i][j]-phi_recon[i-1][j])/dx,0),2),pow(MIN((phi_recon[i+1][j]-phi_recon[i][j])/dx,0),2))+
						MAX(pow(MAX((phi_recon[i][j]-phi_recon[i][j-1])/dy,0),2),pow(MIN((phi_recon[i][j+1]-phi_recon[i][j])/dy,0),2))))-1;
					else 
						grad = float(sqrt(MAX(pow(MAX((phi_recon[i+1][j]-phi_recon[i][j])/dx,0),2),pow(MIN((phi_recon[i][j]-phi_recon[i-1][j])/dx,0),2))+
						MAX(pow(MAX((phi_recon[i][j+1]-phi_recon[i][j])/dy,0),2),pow(MIN((phi_recon[i][j]-phi_recon[i][j-1])/dy,0),2))))-1;
					fact = dtau*grad*phi[i][j]/sqrt(pow(phi[i][j],2)+pow(epsilon,2));
				}
				
				phi_recon[i][j] -= fact;
			}
		}
	}

	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			phi[i][j] = phi_recon[i][j];
}

void Grid::CheckReinit(int iter, float dtau){
	int i,j,ii;
	float fact, grad, den, den1, eps = 1e-1f, epsilon = dx/6.0, d;
	Cell *c;
	Point p;

	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			phi_recon[i][j] = phi[i][j];
	//work inside the band	
	for (i = 1; i < m-1; i++)
		for (j = 1; j < n-1; j++){
			p = Point(i,j);
		//	den = sqrt(pow(phi[i+1][j]-phi[i-1][j],2)+pow(phi[i][j+1]-phi[i][j-1],2))/dx/2;
			if (phi[i][j] > 0)
					den = float(sqrt(MAX(pow(MAX((phi[i][j]-phi[i-1][j])/dx,0),2),pow(MIN((phi[i+1][j]-phi[i][j])/dx,0),2))+
						MAX(pow(MAX((phi[i][j]-phi[i][j-1])/dy,0),2),pow(MIN((phi[i][j+1]-phi[i][j])/dy,0),2))));
				else 
					den = float(sqrt(MAX(pow(MAX((phi[i+1][j]-phi[i][j])/dx,0),2),pow(MIN((phi[i][j]-phi[i-1][j])/dx,0),2))+
						MAX(pow(MAX((phi[i][j+1]-phi[i][j])/dy,0),2),pow(MIN((phi[i][j]-phi[i][j-1])/dy,0),2))));
			if ((phi[i][j]*phi[i-1][j] <= 0)||(phi[i][j]*phi[i+1][j] <= 0)||
				(phi[i][j]*phi[i][j-1] <= 0)||(phi[i][j]*phi[i][j+1] <= 0)){								
				phi_recon[i][j] /= den;//or /1
			}	
			else{
					//outside the narrow band
					if (phi[i][j] > 0)
						grad = float(sqrt(MAX(pow(MAX((phi[i][j]-phi[i-1][j])/dx,0),2),pow(MIN((phi[i+1][j]-phi[i][j])/dx,0),2))+
						MAX(pow(MAX((phi[i][j]-phi[i][j-1])/dy,0),2),pow(MIN((phi[i][j+1]-phi[i][j])/dy,0),2))))-1;
					else 
						grad = float(sqrt(MAX(pow(MAX((phi[i+1][j]-phi[i][j])/dx,0),2),pow(MIN((phi[i][j]-phi[i-1][j])/dx,0),2))+
						MAX(pow(MAX((phi[i][j+1]-phi[i][j])/dy,0),2),pow(MIN((phi[i][j]-phi[i][j-1])/dy,0),2))))-1;
					fact = dtau*grad*phi[i][j]/sqrt(pow(phi[i][j],2)+pow(epsilon,2));
					phi_recon[i][j] -= fact;
			}
		}
	//work only outside the band
	for (ii = 0; ii < iter; ii++){
		for (i=1 ; i<m-1 ; i++){ 
			for (j = 1; j < n-1; j++){
				
				if ((phi[i][j]*phi[i-1][j] <= 0)||(phi[i][j]*phi[i+1][j] <= 0)||
					(phi[i][j]*phi[i][j-1] <= 0)||(phi[i][j]*phi[i][j+1] <= 0)){
					fact = 0;
				}
				else{
					//outside the narrow band
					if (phi[i][j] > 0)
						grad = float(sqrt(MAX(pow(MAX((phi_recon[i][j]-phi_recon[i-1][j])/dx,0),2),pow(MIN((phi_recon[i+1][j]-phi_recon[i][j])/dx,0),2))+
						MAX(pow(MAX((phi_recon[i][j]-phi_recon[i][j-1])/dy,0),2),pow(MIN((phi_recon[i][j+1]-phi_recon[i][j])/dy,0),2))))-1;
					else 
						grad = float(sqrt(MAX(pow(MAX((phi_recon[i+1][j]-phi_recon[i][j])/dx,0),2),pow(MIN((phi_recon[i][j]-phi_recon[i-1][j])/dx,0),2))+
						MAX(pow(MAX((phi_recon[i][j+1]-phi_recon[i][j])/dy,0),2),pow(MIN((phi_recon[i][j]-phi_recon[i][j-1])/dy,0),2))))-1;
					fact = dtau*grad*phi[i][j]/sqrt(pow(phi[i][j],2)+pow(epsilon,2));
				}
				
				phi_recon[i][j] -= fact;
			}
		}
	}
	for (i = 1; i < m-1; i++)
		for (j = 1; j < n-1; j++)
			phi[i][j] = phi_recon[i][j];
		
}

void Grid::Redistance(void){
	int i,j,k, ix, iy, ixmin, ixmax, iymin, iymax, size = 2;
	Cell *c;
	Point p;
	float maxdist = 3.0f*dx, phi_x, phi_y,dist,err,den1,den2;

	// find the minimum distance between a node and the neighbouring particles

	for (i=0; i < m; i++)
		for (j=0; j<n; j++)
			phi_recon[i][j] = maxdist;
	for (i=0; i<m-1; i++)
		for (j=0; j< n-1; j++){
			c = &cell_list[i][j];
		//	if ((fabs(phi[i][j])==3*dx)&&(fabs(phi[i+1][j])==3*dx)&&(fabs(phi[i][j+1])==3*dx)&&(fabs(phi[i+1][j+1])==3*dx))
		//		continue;
			for (k = 0; k < c->num_particles; k++){
				p = c->particles[k];
				ixmin = MAX(0,i-size);//essentially i-2...
				ixmax = MIN(m-2,i+size);//essentially i+2...
				iymin = MAX(0,j-size);
				iymax = MIN(n-2,j+size);
				for (ix = ixmin; ix <= ixmax; ix++)
					for (iy = iymin; iy <= iymax; iy++){
						dist = float(sqrt((p.x-ix*dx-x_lo)*(p.x-ix*dx-x_lo)+(p.y-iy*dy-y_lo)*(p.y-iy*dy-y_lo)));
						//phi_recon[ix][iy] = MIN(dist,phi_recon[ix][iy]);
						if (dist < phi_recon[ix][iy]){
							phi_recon[ix][iy] = dist;
							color[ix][iy] = c->color[k];
						}
					}
			}
		}



	for (i = 1; i < m-1; i++)
		for ( j = 1; j < n-1; j++){
			if ((fabs(phi[i][j]) < 1.0*dx)){
				phi_x = (phi[i+1][j]-phi[i][j]+phi[i+1][j+1]-phi[i][j+1])/(2*dx);
				phi_y = (phi[i][j+1]-phi[i][j]+phi[i+1][j+1]-phi[i+1][j])/(2*dy);
				float norm = float(sqrt(phi_x*phi_x+phi_y*phi_y));
				if (norm > 0){
					phi_x /= norm;
					phi_y /= norm;
				}
				p = Point(i,j);
				dist = 3*dx;
				err = phi[i-1][j]+dx;
				if (err <= 0) dist = MIN(dist,err);
				err = phi[i+1][j]+dx;
				if (err <= 0) dist = MIN(dist,err);
				err = phi[i][j-1]+dx;
				if (err <= 0) dist = MIN(dist,err);
				err = phi[i][j+1]+dx;
				if (err <= 0) dist = MIN(dist,err);
				den1 = 3*dx;
				den2 = 3*dx;
				for (int h = -1; h <= 0; h++)
					for (int l = -1; l <= 0; l++){
						c = &cell_list[i+h][j+l];
						den1 = MIN(den1,FindLineDist(c,p,1));//exact distance to interface...
						den2 = MIN(den2,FindLineDist(c,p,0));
					}
			//	temp_phi1[i][j] = dist;
				temp_phi1[i][j] = MIN(den1-dx,den2);
			}
		}

	//update phi, using the new distance but keeping the old sign
	//note that this will update correctly outside nodes farther than dx from the interface
	//doing only phi[i][j] = phi_recon[i][j]*SIGN(phi[i][j]) is not robust
	for (i = 0; i < m; i++)
		for ( j = 0; j < n; j++){
		//	if ((fabs(phi[i][j]) < 1.0*dx))
		//		phi[i][j] = phi_recon[i][j]*SIGN(temp_phi1[i][j]);
		//	else		
				phi[i][j] = phi_recon[i][j]*SIGN(phi[i][j]);
			phi_recon[i][j] = maxdist;	
		}
}

void Grid::Redistance0(void){
	//first idea from Mark - doesn't work well
	int i,j,k, ix, iy, ixmin, ixmax, iymin, iymax, size = 2;
	Cell *c;
	Point p;
	float maxdist = 3.0f*dx, a, phi_x, phi_y,dist;
		
	
	//find the sign for interface cells
	for (i=1; i<m-1; i++)
		for (j=1; j<n-1; j++){
			//find the gradient at the node
			if ((fabs(phi[i][j])==3*dx)&&(fabs(phi[i+1][j])==3*dx)&&(fabs(phi[i][j+1])==3*dx)&&(fabs(phi[i+1][j+1])==3*dx))
				continue;
			ix = i; iy = j;
			if (i > m-3) ix = m-3; if (j > n-3) iy = n-3; if (i < 1) ix = 1; if (j < 1) iy = 1;
			phi_x = (2*phi[ix+1][iy]-2*phi[ix-1][iy]+phi[ix+1][iy+1]-phi[ix-1][iy+1]+phi[ix+1][iy-1]-phi[ix-1][iy-1])/(4*dx);
			phi_y = (2*phi[ix][iy+1]-2*phi[ix][iy-1]+phi[ix+1][iy+1]-phi[ix+1][iy-1]+phi[ix-1][iy+1]-phi[ix-1][iy-1])/(4*dy);
		//	phi_x = (2*phi[ix+1][iy]-2*phi[ix-1][iy])/(4*dx);
		//	phi_y = (2*phi[ix][iy+1]-2*phi[ix][iy-1])/(4*dy);
			float norm = float(sqrt(pow(phi_x,2)+pow(phi_y,2)));
			if (norm > 0){
				phi_x /= norm;
				phi_y /= norm;
			}
			//populate the particle list
			psize = 0;
			for (int h = -1; h <=0; h++)
				for (int l = -1; l <=0; l++){
					c = &cell_list[i+h][j+l];
					for (k=0; k < c->num_particles; k++){
						 p = c->particles[k];
						 p.x = p.x-x_lo;
						 p.y = p.y-y_lo;
						 plist[psize++] = p;
					}
				}
			temp_phi1[i][j] = 0;
			if (psize > 0){
				a = 0;
				for ( k = 0; k < psize; k++)
					a += -plist[k].x*phi_x-plist[k].y*phi_y;
							
				a = a/psize;
				phi_recon[i][j] = float(fabs(phi[i][j])*SIGN((a + phi_x*i*dx + phi_y*j*dy)));
				temp_phi1[i][j] = 1;
			}
			
		}
	


	//reset phi_recon
	for (i=0; i < m; i++)
		for (j=0; j<n; j++){
			if (temp_phi1[i][j] > 0)
				phi[i][j] = phi_recon[i][j];
		//	temp_phi1[i][j] = 1;
			phi_recon[i][j] = maxdist;
		}

	for (i=0; i<m-1; i++)
		for (j=0; j< n-1; j++){
			c = &cell_list[i][j];
		//	if ((fabs(phi[i][j])==3*dx)&&(fabs(phi[i+1][j])==3*dx)&&(fabs(phi[i][j+1])==3*dx)&&(fabs(phi[i+1][j+1])==3*dx))
		//		continue;
			for (k = 0; k < c->num_particles; k++){
				p = c->particles[k];
				ixmin = MAX(0,i-size);//essentially i-2...
				ixmax = MIN(m-2,i+size);//essentially i+2...
				iymin = MAX(0,j-size);
				iymax = MIN(n-2,j+size);
				for (ix = ixmin; ix <= ixmax; ix++)
					for (iy = iymin; iy <= iymax; iy++){
						dist = float(sqrt((p.x-ix*dx-x_lo)*(p.x-ix*dx-x_lo)+(p.y-iy*dy-y_lo)*(p.y-iy*dy-y_lo)));
						//phi_recon[ix][iy] = MIN(dist,phi_recon[ix][iy]);
						if (dist < phi_recon[ix][iy]){
							phi_recon[ix][iy] = dist;
							color[ix][iy] = c->color[k];
						}
					}
			}
		}
	for (i=0; i < m; i++)
		for (j=0; j<n; j++){
			if (phi[i][j] < dx)
				phi[i][j] = phi_recon[i][j]*SIGN((phi[i][j]));
			phi_recon[i][j] = maxdist;
		}
}

int Grid::CountParticles(void){
	int i,j,num;
	Cell *c;
	num = 0;
	for (i=0; i<m-1; i++)
		for (j=0; j< n-1; j++){
			c = &cell_list[i][j];
			if (c->num_particles > 0)
				num += c->num_particles;
		}
	return num; 
}


void Grid::Redistance1(void){
	int i,j,k, ix, iy, ixmin, ixmax, iymin, iymax, size = 2, h, l;
	Cell *d;
	Point p;
	float maxdist = 3.0f*dx, ldx, ldy, ldxi, ldyj, dist,q,c=0.3,rho=1.5;
		
	for (i=0; i<m; i++)
		for (j=0; j< n; j++)
			temp_phi2[i][j] = 0;
	for (i=1; i<m-1; i++)
		for (j=1; j< n-1; j++){
			float lambda = 0, wght = 0, wght1 = 0, twght = 0;  
			int lim = 1;
			//use the cells about (i,j) to find correction for phi
			for (h = -lim; h <=lim-1; h++){
				if ((0<(i+h))&&((i+h)<(m-1))){
					for (l = -lim; l <=lim-1; l++){
						if ((0<(j+l))&&((j+l)<(n-1))){
							d = &cell_list[i+h][j+l];
							for (k=0; k < d->num_particles; k++){
								 p = d->particles[k];
								 ldx = p.x-x_lo-(i+h)*dx; ldx /= dx;
								 ldy = p.y-y_lo-(j+l)*dy; ldy /= dy;
								 ldxi = p.x-x_lo-i*dx; ldxi /= dx;
								 ldyj = p.y-y_lo-j*dy; ldyj /= dy;
								 wght = 1/(1+pow(ldxi*ldxi+ldyj*ldyj,3));
								 q = sqrt(ldxi*ldxi+ldyj*ldyj)/rho;
								 if (q<1)
									 wght1= (exp(-q*q/(c*c))-exp(-1/(c*c)))/(1-exp(-1/(c*c)));
								 else wght1 = 0;
								 twght += wght1;
								 if ((i+h<m-2)&&(j+l<n-2))
									 lambda += wght1*Interpolate(phi,3,i+h,j+l,ldx,ldy);//1st or 3rd order interp
								 else
									lambda += wght1*Interpolate(phi,1,i+h,j+l,ldx,ldy);
							}
						}
					}
				}
			}
			if (twght > 0){
				phi_recon[i][j] = phi[i][j]-lambda/twght;
			//	phi[i][j] -= lambda/twght;
				temp_phi2[i][j] = 1;
			}
		}

	//reset phi_recon
	for (i=0; i < m; i++)
		for (j=0; j<n; j++){
			if (temp_phi2[i][j])
				phi[i][j] = phi_recon[i][j];
			phi_recon[i][j] = maxdist;
		}
		//new fix
/*	for (i=0; i<m-1; i++)
		for (j=0; j< n-1; j++){
			d = &cell_list[i][j];
			for (k = 0; k < d->num_particles; k++){
				p = d->particles[k];
				ixmin = MAX(0,i-size);//essentially i-2...
				ixmax = MIN(m-2,i+size);//essentially i+2...
				iymin = MAX(0,j-size);
				iymax = MIN(n-2,j+size);
				for (ix = ixmin; ix <= ixmax; ix++)
					for (iy = iymin; iy <= iymax; iy++){
						dist = float(sqrt((p.x-ix*dx-x_lo)*(p.x-ix*dx-x_lo)+(p.y-iy*dy-y_lo)*(p.y-iy*dy-y_lo)));
						if (dist < phi_recon[ix][iy]){
							phi_recon[ix][iy] = dist;
						}
					}
			}
		}
	for (i=0; i < m; i++)
		for (j=0; j<n; j++){
			if (phi_recon[i][j] < maxdist){
				phi_recon[i][j] = MAX((fabs(phi[i][j])),(phi_recon[i][j]-dx));
			}
			else
				phi_recon[i][j] = MAX(dx,(fabs(phi[i][j])));
		}
	for (i=0; i < m; i++)
		for (j=0; j<n; j++){
			phi[i][j] = phi_recon[i][j]*SIGN((phi[i][j]));
			phi_recon[i][j] = maxdist;
		}*/
}

void Grid::CorrectLevelSet(void){
	float dist,err;
	int i, j;
	Point p;

	for (i = 0; i < m; i++){
		for (j = 0; j < n; j++){
			if (fabs(phi[i][j]) < dx){
				dist = 3*dx;
			 	err = phi[i-1][j]+dx;
				if (err <= 0) dist = MIN(dist,err);
				err = phi[i+1][j]+dx;
				if (err <= 0) dist = MIN(dist,err);
				err = phi[i][j-1]+dx;
				if (err <= 0) dist = MIN(dist,err);
				err = phi[i][j+1]+dx;
				if (err <= 0) dist = MIN(dist,err);
				//if (dist < 3*dx)
				//	phi[i][j] = dist;
				//phi[i][j] = SIGN((phi[i-1][j]+phi[i+1][j]+phi[i][j-1]+phi[i][j+1])/4)*fabs(phi[i][j]);
			}
		}
	}
}

void Grid::FixLevelSetInside(float threshold, float eps){
	Cell *d;
	float dist, dist0;
	int i, j, h, l, lim;
	Point p;

	for (i = 1; i < m-1; i++){
		for (j = 1; j < n-1; j++){
			p = Point(i,j);//a node
			dist = 3*dx;
			dist0 = 3*dx;
			if (fabs(phi[i][j]) < 1.0000f*threshold){//or 1.4142f
			    lim = 2;
				//use the cells inside a 4x4 nbhd of p to find its minimal distance to the threshold2 level set
				for (h = -lim; h <=lim-1; h++){
					if ((0<i+h)&&(i+h<m-1)){
						for (l = -lim; l <=lim-1; l++){
							if ((0<j+l)&&(j+l<n-1)){
								d = &cell_list[i+h][j+l];
								if ((phi[i+h][j+l] <= 0)||(phi[i+h+1][j+l] <= 0)||(phi[i+h][j+l+1] <= 0)||(phi[i+h+1][j+l+1] <= 0))
									dist = MIN(dist,FindLineDist(d,p,threshold));
							}
						}
					}
				}

			//	if ((phi[i][j]<eps*threshold))
			//		continue;
			//	else if ((phi[i][j] >= eps*threshold) && (dist < 2*threshold))
				//	phi[i][j] = dist-threshold;//this is the "true" signed distance, based on the "good" outside level -1 set
				//	phi[i][j] = MIN(phi[i][j],dist0);
					temp_phi1[i][j] = dist-threshold;
				//		temp_phi1[i][j] = dist;
			}
		}
	}
}

void Grid::FixLevelSetOutside(float threshold, float eps){
	Cell *d;
	float dist, dist0,epsilon=dx;//eps=0.01 for zalesak
	int i,j;
	Point p;

	for (i = 1; i < m-1; i++){
		for (j = 1; j < n-1; j++){
			p = Point(i,j);//a node
			dist = 3*dx;
			dist0 = 3*dx;
			if (fabs(phi[i][j]) < threshold){
				int lim = 2;
				//use the cells inside a 4x4 nbhd of p to find its minimal distance to the threshold2 level set
				for (int h = -lim; h <=lim-1; h++){
					if ((0<i+h)&&(i+h<m-1)){
						for (int l = -lim; l <=lim-1; l++){
							if ((0<j+l)&&(j+l<n-1)){
								d = &cell_list[i+h][j+l];
								if ((phi[i+h][j+l] >= 0)||(phi[i+h+1][j+l] >= 0)||(phi[i+h][j+l+1] >= 0)||(phi[i+h+1][j+l+1] >= 0))
									dist = MIN(dist,FindLineDist(d,p,threshold));
							}
						}
					}
				}
				/*if (dist > fabs(temp_phi1[i][j])&&(fabs(temp_phi1[i][j])<epsilon))
					phi[i][j] = temp_phi1[i][j];
				else if (dist < epsilon)
					phi[i][j] = threshold-dist;*/
			/*	if (dist > fabs(temp_phi1[i][j]))//{
					if (fabs(temp_phi1[i][j]) < fabs(phi[i][j]))
						phi[i][j] = temp_phi1[i][j];
				//}
				else{ 
					if (dist < epsilon)
						if (dist < fabs(phi[i][j]))
							phi[i][j] = threshold-dist;
				}*/
				if (fabs(threshold-dist-phi[i][j]) > fabs(temp_phi1[i][j]-phi[i][j]))
						phi[i][j] = temp_phi1[i][j];
				else
						phi[i][j] = threshold-dist;
			
			}
		}
	}
}


bool Grid::RedistributeParticlesGrid(float threshold){
	//this is done after solving the reconnection problem and after we get the new cell segments from marching squares
	//it simply uses the reconstructed segments to insert new particles along them
	
	int i,j,k,l,index;
	Cell* c;
	Point a,b;
	bool done = 0;
	float dist,ldx,ldy,lim1 = 1.5;
	Point p;

	for (i = 0; i < m-1; i++)
		for (j = 0; j < n-1; j++){
			c = &cell_list[i][j];
			index = 0;
			if ((fabs(phi[i][j]) <= 1.5*dx)&&(fabs(phi[i+1][j]) <= 1.5*dx)&&(fabs(phi[i][j+1]) <= 1.5*dx)&&(fabs(phi[i+1][j+1]) <= 1.5*dx)){
				if (CellChangesSignInside(c, 0)){
					if (c->num_particles <= 1){
						done = 1;
						AddParticlesToCell(c, threshold);
					}
				}
				else{
					if (c->num_particles > 0){
						bool ip_true = 0;
						for (l = 0; l < c->num_particles; l++){
							if (c->initial_particles[l])
								ip_true = 1;
						}
						
						//if there are initial particles in this cell
						if (ip_true){//keep just the initial particles
							for (l = 0; l < c->num_particles; l++){
								if (c->initial_particles[l]){
									c->particles[index] = c->particles[l];
									c->initial_particles[index] = true;
									index++;
								}
							}
							//make sure all the rest of the particles are marked "added" 
							c->num_particles = index;
							for (l = index; l < MAX_CELL_PARTICLES; l++)
								c->initial_particles[l] = 0;
							//if inside particles (phi>0) delete everything, else keep
							if (phi[i][j] > 0)
								c->num_particles = 0;
						}
						else{//delete all particles
							for (l = 0; l < MAX_CELL_PARTICLES; l++)
								c->initial_particles[l] = 0;
							c->num_particles = 0;
						}

						//delete everything
					//	for (l = 0; l < MAX_CELL_PARTICLES; l++)
					//		c->initial_particles[l] = 0;
					//	c->num_particles = 0;
					/*
						for (k=0; k < c->num_particles; k++)
							c->initial_particles_new[k] = 1;
						for (k=0; k < c->num_particles; k++){
							 p = c->particles[k];
							 ldx = p.x-x_lo-i*dx; ldx /= dx;
							 ldy = p.y-y_lo-j*dy; ldy /= dy;
							 p = Point(i+ldx,j+ldy);
							 dist = FindLineDist(c,p,0)/dx;
							 if (dist > lim1)//mark the particle to be deleted
								 c->initial_particles_new[k] = 0;
						}

						index = 0;
						for (k=0; k < c->num_particles; k++){
							if (c->initial_particles_new[k])
								c->particles[index++] = c->particles[k];
						}
						c->num_particles = index;
					*/
					}
				}
			}
			else{
				float lim = 1.5;
				//delete all particles if too far from zero level set
			/*	if ((fabs(phi[i][j])> lim*dx)&&(fabs(phi[i][j+1])> lim*dx)&&(fabs(phi[i+1][j])> lim*dx)&&(fabs(phi[i+1][j+1])> lim*dx)){
					for (l = 0; l < MAX_CELL_PARTICLES; l++)
						c->initial_particles[l] = 0;
					c->num_particles = 0;
				}*/
/*
				for (k=0; k < c->num_particles; k++)
					c->initial_particles_new[k] = 1;
				for (k=0; k < c->num_particles; k++){
					 p = c->particles[k];
					 ldx = p.x-x_lo-i*dx; ldx /= dx;
					 ldy = p.y-y_lo-j*dy; ldy /= dy;
					 p = Point(i+ldx,j+ldy);
					 dist = FindLineDist(c,p,0)/dx;
					 if (dist > lim1)//mark the particle to be deleted
						 c->initial_particles_new[k] = 0;
				}

				index = 0;
				for (k=0; k < c->num_particles; k++){
					if (c->initial_particles_new[k])
						c->particles[index++] = c->particles[k];
				}
				c->num_particles = index;
*/
			}
		}
	psize = CountParticles();
	return done;
}

bool Grid::RedistributeParticles(float threshold){
	//this is done after solving the reconnection problem and after we get the new cell segments from marching squares
	//it simply uses the reconstructed segments to insert new particles along them
	
	int i,j,index;
	Cell* c;
	Point a,b;
	bool done = 0;
	for (i = 0; i < m-1; i++)
		for (j = 0; j < n-1; j++){
			c = &cell_list[i][j];
			if ((fabs(phi[i][j]) <= 1.5*dx)&&(fabs(phi[i+1][j]) <= 1.5*dx)&&(fabs(phi[i][j+1]) <= 1.5*dx)&&(fabs(phi[i+1][j+1]) <= 1.5*dx)){
				if (CellChangesSignInside(c, 0)){
					if (c->num_particles <= 1){
						done = 1;
						AddParticlesToCell(c, threshold);
					}
				}
			}
		}
	return done;
}

void Grid::AddParticlesToCell(Cell* c, float threshold){
	//version1, based on the level set linear reconstruction
	Point a,b,p;
	float d,t,dist,eps = 0.1*dx;
	int k,l,np,lim = 1;
	Cell *dd;
	for (k = 0; k < c->num_segmentpts/2; k++){
		a = c->segments[2*k];
		b = c->segments[2*k+1];
		d = sqrt(pow(a.x-b.x,2)+pow(a.y-b.y,2));
		if (d > 0.5){
			np = int(floor(d/0.5f));//i.e. 2 particles per UNIT segment on the average....
			if (np < 2)
				np = 2;
		//	if (np > 4)
		//		np = 4;
			//c->num_particles = 0;
			for (l = 0; l < np; l++){
				t = l/float(np-1);	
			//	if (c->num_particles <= MAX_CELL_PARTICLES){
				//first check if the point is not bogus, i.e. too far from the 1 level set
					dist = 3*dx;
					//p = Point(x_lo+dx*(t*a.x+(1-t)*b.x),y_lo+dy*(t*a.y+(1-t)*b.y));
					p = Point(t*a.x+(1-t)*b.x,t*a.y+(1-t)*b.y);
					//use the cells inside a 3x3 nbhd of p to find its minimal distance to the threshold (-dx) level set
					for (int h = -lim; h <=lim; h++)
						for (int l = -lim; l <=lim; l++){
							if ((0<=c->i+h)&&(c->i+h<m-1)&&(0<=c->j+l)&&(c->j+l<n-1)){
								dd = &cell_list[c->i+h][c->j+l];
								if (dd->num_segmentpts1 > 0)
									dist = MIN(dist,FindLineDist(dd,p,1));
							}
						}
					//	if ((dist < threshold + eps)&&(dist > threshold - eps)&&(c->num_particles < MAX_CELL_PARTICLES)){
						if (c->num_particles < MAX_CELL_PARTICLES){
							c->particles[c->num_particles] = Point(x_lo+dx*(t*a.x+(1-t)*b.x),y_lo+dy*(t*a.y+(1-t)*b.y));
							c->initial_particles[c->num_particles] = false;
							c->color[c->num_particles] = color[c->i][c->j];
							c->num_particles++;
						}
			//	}
			}
		}
	}
}

void Grid::AdvectParticles(float dt){
	int i,j,k,ctr;
	Cell* c;
	bool* stays_inside = new bool [MAX_CELL_PARTICLES];
	for ( k = 0; k < MAX_CELL_PARTICLES; k++)
		stays_inside[k] = 0;

	for (i = 0; i < m-1; i++)
		for ( j = 0; j < n-1; j++){

			c = &cell_list[i][j];
			for ( k = 0; k < MIN(c->num_particles,MAX_CELL_PARTICLES); k++){
				if (AdvectParticle(c, k ,dt))
					stays_inside[k] = 1;
				else
					stays_inside[k] = 0;
			}
			//update the cell particles (some stayed, some flew away)
			ctr = 0;
			for ( k = 0; k < MIN(c->num_particles,MAX_CELL_PARTICLES); k++)
				if (stays_inside[k]){
					c->particles[ctr] = c->particles[k];
					c->color[ctr] = c->color[k];
					if (c->initial_particles[k] == 1)
						c->initial_particles[ctr] = 1;
					else
						c->initial_particles[ctr] = 0;
					ctr++;
				}
			c->num_particles = ctr;
		}
	//add the new particles to each cell
	for (i = 0; i < m-1; i++)
		for ( j = 0; j < n-1; j++){
			c = &cell_list[i][j];
			ctr = MIN(c->num_particles,MAX_CELL_PARTICLES);
			for ( k = 0; k < MIN(c->num_newparticles,MAX_CELL_PARTICLES-ctr); k++){
				c->particles[ctr+k] = c->new_particles[k];
				c->color[ctr+k] = c->new_color[k];
				if (c->initial_particles_new[k] == 1)
					c->initial_particles[ctr+k] = 1;
				else
					c->initial_particles[ctr+k] = 0;
			}
			c->num_particles += c->num_newparticles;
			c->num_particles = MIN(c->num_particles,MAX_CELL_PARTICLES); 
			c->num_newparticles = 0;
			for ( k = 0; k < MAX_CELL_PARTICLES; k++){
				c->initial_particles_new[k] = 0;
			}
		}

	delete [] stays_inside;
}

bool Grid::AdvectParticle(Cell* c, int k, float dt){
 
	Point p = c->particles[k];
	int ix = int((p.x-x_lo)/dx), iy = int((p.y-y_lo)/dy);
	if (ix < 0 ) ix = 0; if (ix > m-2) ix = m-2;	
	if (iy < 0 ) iy = 0; if (iy > n-2) iy = n-2;
	//cout << ix-c->i << " " << iy - c->j << " " << endl;
	float ldx = (p.x-x_lo)/dx-ix, ldy = (p.y-y_lo)/dy-iy, x, y;
	float u_vel = (1-ldy)*((1-ldx)*u[ix][iy]+ldx*u[ix+1][iy])+ldy*((1-ldx)*u[ix][iy+1]+ldx*u[ix+1][iy+1]);
	float v_vel = (1-ldy)*((1-ldx)*v[ix][iy]+ldx*v[ix+1][iy])+ldy*((1-ldx)*v[ix][iy+1]+ldx*v[ix+1][iy+1]);
	//add the new point
	float px = p.x+dt*u_vel, py = p.y+dt*v_vel;
	int ix1 = int((px-x_lo)/dx), iy1 = int((py-y_lo)/dy);
	if (ix1 < 0 ) ix1 = 0; if (ix1 > m-2) ix1 = m-2;	
	if (iy1 < 0 ) iy1 = 0; if (iy1 > n-2) iy1 = n-2;
	float ldx1 = (px-x_lo)/dx-ix1, ldy1 = (py-y_lo)/dy-iy1;
	float u_vel1 = (1-ldy1)*((1-ldx1)*u[ix1][iy1]+ldx1*u[ix1+1][iy1])+ldy1*((1-ldx1)*u[ix1][iy1+1]+ldx1*u[ix1+1][iy1+1]);
	float v_vel1 = (1-ldy1)*((1-ldx1)*v[ix1][iy1]+ldx1*v[ix1+1][iy1])+ldy1*((1-ldx1)*v[ix1][iy1+1]+ldx1*v[ix1+1][iy1+1]);
	px = p.x + 0.5f*(u_vel+u_vel1)*dt;
	py = p.y + 0.5f*(v_vel+v_vel1)*dt;
	//px = p.x + u_vel1*dt;
	//py = p.y + v_vel1*dt;
	//this was Crank-Nicholson (not Runge Kutta 2)...
	x = (px-x_lo)/dx; y = (py-y_lo)/dy;
	//if (x<0.5f) x=0.5f; if (x>m-1.5f) x=m-1.5f;
	//if (y<0.5f) y=0.5f; if (y>n-1.5f) y=n-1.5f; 
	if (x<0.0f){ 
		if (bcxlo == 0)
			x = x+(m-1);
		else
			x=0.0f;
	}
	if (x>=m-1){ 
		if (bcxhi == 0)
			x = x-(m-1);
		else
			x=m-1.00001f;
	}
	if (y<0.0f){ 
		if (bcylo == 0)
			y = y+(n-1);
		else
			y=0.0f;
	}
	if (y>=n-1){
		if (bcyhi == 0)
			y = y-(n-1);
		else
			y=n-1.00001f;
	}
	ix = int(x); iy = int(y);
	bool periodicx = ((bcxlo == 0)||(bcxhi == 0)),periodicy = ((bcylo == 0)||(bcyhi == 0));
	if ((ix == c->i) && (iy == c->j)){//no deletions, the particle stays inside the cell
		c->particles[k].x = px;
		c->particles[k].y = py;
		return true;
	}
	else{//point p will be deleted from the cell list at the end of the upper function loop
		int ind = cell_list[ix][iy].num_newparticles;
		if (ind < MAX_CELL_PARTICLES){
			cell_list[ix][iy].new_particles[ind].x = (periodicx) ? x_lo+x*dx : px;
			cell_list[ix][iy].new_particles[ind].y = (periodicy) ? y_lo+y*dy : py;
			cell_list[ix][iy].new_color[ind] = c->color[k];
			if (c->initial_particles[k] == 1){
				cell_list[ix][iy].initial_particles_new[ind] = 1;
				c->initial_particles[k] = 0;
			}
			else
				cell_list[ix][iy].initial_particles_new[ind] = 0;
			cell_list[ix][iy].num_newparticles++;
		}
		return false;
	}
}

int Grid::CountInterfaceCells(void){
	int totalcount = 0;
	Cell* c;
	for (int i = 0; i<m-1; i++)
		for (int j = 0; j < n-1; j++){
			c = &cell_list[i][j];
			if (CellChangesSignInside(c,0))
				totalcount += 1;
		}
	return totalcount;
}

void Grid::FindZeroLevelCurve(float scalex, float scaley, float threshold, bool cellbased){
	
	//all the segments are oriented ccw for a level set negative inside, positive outside
	int u, r, M, N, i;
	float cc, cr, uc, ur, xnew, ynew, x, y, area = 0, l1, l2, l3, l4;
	float sx = scalex, sy = scaley;//why integer??
	Cell *c;
	Point * ss;
	int size;
	if (threshold == 0)
		ss = s;
	else if (threshold > 0)
		ss = s1;
	else
		ss = s2;

	size = 0;
	if (cellbased){
		N = n-2; M = m-2; }
	else{
		N = n-1; M = m-1; }

	for (int h=0; h<N; h++) {
		for (int w=0; w<M; w++) {
			u = h + 1;
			r = w + 1;
			if (cellbased){
				cc = cell_list[w][h].phi-threshold;
				cr = cell_list[r][h].phi-threshold;
				uc = cell_list[w][u].phi-threshold;
				ur = cell_list[r][u].phi-threshold;	
				x = (w+0.5f)*sx;
				y = (h+0.5f)*sy;
			}
			else{
				cc = phi[w][h]-threshold;
				cr = phi[r][h]-threshold;
				uc = phi[w][u]-threshold;
				ur = phi[r][u]-threshold;
				x = w*sx;
				y = h*sy;
			}
			c = &cell_list[w][h];
			if (threshold == 0)
				c->num_segmentpts = 0;
			else
				c->num_segmentpts1 = 0;
			
			if (cc > 0 && cr > 0 && uc > 0 && ur > 0){
				area += sx*sy;
				continue;
			}
			else if (cc > 0 && cr > 0 && uc <= 0 && ur > 0){
				xnew = x;
				ynew = y + (cc / (cc - uc) * sy );
				ss[size++] = Point(xnew, ynew);
				l1 = y+sy-ynew;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cc/(cc-uc)*color[w][u][i]-uc/(cc-uc)*color[w][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);					
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x + (uc / (uc - ur) * sx );
				ynew = y + sy;
				ss[size++] = Point(xnew,ynew);
				l2 = xnew-x;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = uc/(uc-ur)*color[r][u][i]-ur/(uc-ur)*color[w][u][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				area += sx*sy-l1*l2/2;
				continue;
			}
			else if (cc > 0 && cr > 0 && uc > 0 && ur <= 0){
				xnew = x + (uc / (uc - ur) * sx );
				ynew = y + sy;
				ss[size++] = Point(xnew,ynew);
				l1 = x+sx-xnew;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = uc/(uc-ur)*color[r][u][i]-ur/(uc-ur)*color[w][u][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x + sx;
				ynew = y + (cr / (cr - ur) * sy );
				ss[size++] = Point(xnew,ynew);
				l2 = y+sy-ynew;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cr/(cr-ur)*color[r][u][i]-ur/(cr-ur)*color[r][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				area += sx*sy-l1*l2/2;
				continue;
			}
			else if (cc > 0 && cr > 0 && uc <= 0 && ur <= 0){
				xnew = x;
				ynew = y + (cc / (cc - uc) * sy );
				ss[size++] = Point(xnew,ynew);
				l1 = ynew-y;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cc/(cc-uc)*color[w][u][i]-uc/(cc-uc)*color[w][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x + sx;
				ynew = y + (cr / (cr - ur) * sy );
				ss[size++] = Point(xnew,ynew);
				l2 = ynew-y;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cr/(cr-ur)*color[r][u][i]-ur/(cr-ur)*color[r][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				area += sx*(l1+l2)/2;
				continue;
			}
			else if (cc > 0 && cr <= 0 && uc > 0 && ur > 0){
				xnew = x + sx;
				ynew = y + (cr / (cr - ur) * sy );
				ss[size++] = Point(xnew,ynew);
				l1 = ynew-y;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cr/(cr-ur)*color[r][u][i]-ur/(cr-ur)*color[r][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x + (cc / (cc - cr) * sx );
				ynew = y;
				ss[size++] = Point(xnew,ynew);
				l2 = x+sx-xnew;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cc/(cc-cr)*color[r][h][i]-cr/(cc-cr)*color[w][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				area += sx*sy-l1*l2/2;
				continue;
			}
		/*	else if (cc > 0 && cr <= 0 && uc <= 0 && ur > 0){
				xnew = x;
				ynew = y + (cc / (cc - uc) * sy );
				ss[size++] = Point(xnew,ynew);
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x + (uc / (uc - ur) * sx );
				ynew = y + sy;
				ss[size++] = Point(xnew,ynew);
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x + sx;
				ynew = y + (cr / (cr - ur) * sy );
				ss[size++] = Point(xnew,ynew);
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x + (cc / (cc - cr) * sx );
				ynew = y;
				ss[size++] = Point(xnew,ynew);
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				continue;
			}*/
			else if (cc > 0 && cr <= 0 && uc <= 0 && ur > 0){
				xnew = x + sx;
				ynew = y + (cr / (cr - ur) * sy );
				ss[size++] = Point(xnew,ynew);
				l1 = y+sy-ynew;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cr/(cr-ur)*color[r][u][i]-ur/(cr-ur)*color[r][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x + (uc / (uc - ur) * sx );
				ynew = y + sy;
				ss[size++] = Point(xnew,ynew);
				l2 = x+sx-xnew;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = uc/(uc-ur)*color[r][u][i]-ur/(uc-ur)*color[w][u][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x;
				ynew = y + (cc / (cc - uc) * sy );
				ss[size++] = Point(xnew,ynew);
				l3 = ynew-y;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cc/(cc-uc)*color[w][u][i]-uc/(cc-uc)*color[w][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x + (cc / (cc - cr) * sx );
				ynew = y;
				ss[size++] = Point(xnew,ynew);
				l4 = xnew-x;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cc/(cc-cr)*color[r][h][i]-cr/(cc-cr)*color[w][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				area += l1*l2/2+l3*l4/2;
				continue;
			}
			else if (cc > 0 && cr <= 0 && uc > 0 && ur <= 0){
				xnew = x + (uc / (uc - ur) * sx );
				ynew = y + sy;
				ss[size++] = Point(xnew,ynew);
				l1 = xnew-x;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = uc/(uc-ur)*color[r][u][i]-ur/(uc-ur)*color[w][u][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x + (cc / (cc - cr) * sx );
				ynew = y;
				ss[size++] = Point(xnew,ynew);
				l2 = xnew-x;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cc/(cc-cr)*color[r][h][i]-cr/(cc-cr)*color[w][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				area += sy*(l1+l2)/2;
				continue;
			}
			else if (cc > 0 && cr <= 0 && uc <= 0 && ur <= 0){
				xnew = x;
				ynew = y + (cc / (cc - uc) * sy );
				l1 = ynew-y;
				ss[size++] = Point(xnew,ynew);
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cc/(cc-uc)*color[w][u][i]-uc/(cc-uc)*color[w][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x + (cc / (cc - cr) * sx );
				ynew = y;
				ss[size++] = Point(xnew,ynew);
				l2 = xnew-x;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cc/(cc-cr)*color[r][h][i]-cr/(cc-cr)*color[w][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				area += l1*l2/2;
				continue;
			}
			else if (cc <= 0 && cr > 0 && uc > 0 && ur > 0){
				xnew = x + (cc / (cc - cr) * sx );
				ynew = y;
				ss[size++] = Point(xnew,ynew);
				l1 = xnew-x;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cc/(cc-cr)*color[r][h][i]-cr/(cc-cr)*color[w][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x;
				ynew = y + (cc / (cc - uc) * sy );
				ss[size++] = Point(xnew,ynew);	
				l2 = ynew-y;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cc/(cc-uc)*color[w][u][i]-uc/(cc-uc)*color[w][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				area += sx*sy-l1*l2/2;
				continue;
			}
			else if (cc <= 0 && cr > 0 && uc <= 0 && ur > 0){
				xnew = x + (cc / (cc - cr) * sx );
				ynew = y;
				ss[size++] = Point(xnew,ynew);
				l1 = xnew-x;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cc/(cc-cr)*color[r][h][i]-cr/(cc-cr)*color[w][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x + (uc / (uc - ur) * sx );
				ynew = y + sy;
				ss[size++] = Point(xnew,ynew);
				l2 = xnew-x;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = uc/(uc-ur)*color[r][u][i]-ur/(uc-ur)*color[w][u][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				area += sx*sy-sy*(l1+l2)/2;
				continue;
			}
			else if (cc <= 0 && cr > 0 && uc > 0 && ur <= 0){
				xnew = x + (uc / (uc - ur) * sx );
				ynew = y + sy;
				ss[size++] = Point(xnew,ynew);
				l1 = xnew-x;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = uc/(uc-ur)*color[r][u][i]-ur/(uc-ur)*color[w][u][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x;
				ynew = y + (cc / (cc - uc) * sy );
				ss[size++] = Point(xnew,ynew);
				l2 = y+sy-ynew;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cc/(cc-uc)*color[w][u][i]-uc/(cc-uc)*color[w][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x + (cc / (cc - cr) * sx );
				ynew = y;
				ss[size++] = Point(xnew,ynew);
				l3 = x+sx-xnew;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cc/(cc-cr)*color[r][h][i]-cr/(cc-cr)*color[w][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x + sx;
				ynew = y + (cr / (cr - ur) * sy );
				ss[size++] = Point(xnew,ynew);
				l4 = ynew-y;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cr/(cr-ur)*color[r][u][i]-ur/(cr-ur)*color[r][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				area += l1*l2/2+l3*l4/2;
				continue;
			}
			else if (cc <= 0 && cr > 0 && uc <= 0 && ur <= 0){
				xnew = x + (cc / (cc - cr) * sx );
				ynew = y;
				ss[size++] = Point(xnew,ynew);
				l1 = x+sx-xnew;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cc/(cc-cr)*color[r][h][i]-cr/(cc-cr)*color[w][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x + sx;
				ynew = y + (cr / (cr - ur) * sy );
				l2 = ynew-y;
				ss[size++] = Point(xnew,ynew);
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cr/(cr-ur)*color[r][u][i]-ur/(cr-ur)*color[r][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				area += l1*l2/2;
				continue;
			}
			else if (cc <= 0 && cr <= 0 && uc > 0 && ur > 0){
				xnew = x + sx;
				ynew = y + (cr / (cr - ur) * sy );
				ss[size++] = Point(xnew,ynew);
				l1 = ynew-y;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cr/(cr-ur)*color[r][u][i]-ur/(cr-ur)*color[r][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x;
				ynew = y + (cc / (cc - uc) * sy );
				ss[size++] = Point(xnew,ynew);	
				l2 = ynew-y;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cc/(cc-uc)*color[w][u][i]-uc/(cc-uc)*color[w][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				area += sx*sy-sx*(l1+l2)/2;
				continue;
			}
			else if (cc <= 0 && cr <= 0 && uc <= 0 && ur > 0){
				xnew = x + sx;
				ynew = y + (cr / (cr - ur) * sy );
				ss[size++] = Point(xnew,ynew);
				l1 = y+sy-ynew;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cr/(cr-ur)*color[r][u][i]-ur/(cr-ur)*color[r][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x + (uc / (uc - ur) * sx );
				ynew = y + sy;
				ss[size++] = Point(xnew,ynew);
				l2 = x+sx-xnew;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = uc/(uc-ur)*color[r][u][i]-ur/(uc-ur)*color[w][u][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				area += l1*l2/2;
				continue;
			}
			else if (cc <= 0 && cr <= 0 && uc > 0 && ur <= 0){
				xnew = x + (uc / (uc - ur) * sx );
				ynew = y + sy;
				ss[size++] = Point(xnew,ynew);
				l1 = xnew-x;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = uc/(uc-ur)*color[r][u][i]-ur/(uc-ur)*color[w][u][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				xnew = x;
				ynew = y + (cc / (cc - uc) * sy );
				ss[size++] = Point(xnew,ynew);
				l2 = y+sy-ynew;
				for (i = 0; i < 3; i++)
					pointcolor[size-1][i] = cc/(cc-uc)*color[w][u][i]-uc/(cc-uc)*color[w][h][i];
				if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
				area += l1*l2/2;
				continue;
			}
			else if (cc <= 0 && cr <= 0 && uc <= 0 && ur <= 0){
				if (cc == 0 && cr == 0 && uc < 0 && ur < 0){
					xnew = x;
					ynew = y;
					ss[size++] = Point(xnew,ynew);
					for (i = 0; i < 3; i++)
						pointcolor[size-1][i] = color[w][h][i];
					if (threshold == 0)
						c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
					else 
						c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
					xnew = x+sx;
					ynew = y;
					ss[size++] = Point(xnew,ynew);	
					for (i = 0; i < 3; i++)
						pointcolor[size-1][i] = color[r][h][i];
					if (threshold == 0)
						c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
					else 
						c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
					continue;
				}
				else if (cc == 0 && cr < 0 && uc == 0 && ur < 0){
					xnew = x;
					ynew = y+sy;
					ss[size++] = Point(xnew,ynew);
					for (i = 0; i < 3; i++)
						pointcolor[size-1][i] = color[w][u][i];
					if (threshold == 0)
						c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
					else 
						c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
					xnew = x;
					ynew = y;
					ss[size++] = Point(xnew,ynew);	
					for (i = 0; i < 3; i++)
						pointcolor[size-1][i] = color[w][h][i];
					if (threshold == 0)
						c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
					else 
						c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
					continue;
				}
				else if (cc < 0 && cr == 0 && uc < 0 && ur == 0){
					xnew = x+sx;
					ynew = y;
					ss[size++] = Point(xnew,ynew);
					for (i = 0; i < 3; i++)
						pointcolor[size-1][i] = color[r][h][i];
					if (threshold == 0)
						c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
					else 
						c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
					xnew = x+sx;
					ynew = y+sy;
					ss[size++] = Point(xnew,ynew);	
					for (i = 0; i < 3; i++)
						pointcolor[size-1][i] = color[r][u][i];
					if (threshold == 0)
						c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
					else 
						c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
					continue;
				}
				else if (cc < 0 && cr < 0 && uc == 0 && ur == 0){
					xnew = x+sx;
					ynew = y+sy;
					ss[size++] = Point(xnew,ynew);
					for (i = 0; i < 3; i++)
						pointcolor[size-1][i] = color[r][u][i];
					if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
					xnew = x;
					ynew = y+sy;
					for (i = 0; i < 3; i++)
						pointcolor[size-1][i] = color[w][u][i];
					ss[size++] = Point(xnew,ynew);	
					if (threshold == 0)
					c->segments[c->num_segmentpts++] = Point(xnew/sx,ynew/sy);
				else 
					c->segments1[c->num_segmentpts1++] = Point(xnew/sx,ynew/sy);
					continue;
				}
				else continue;
			}
		}	
	}
	area /= sx*sy;
	mass = area*dx*dy;
	if (threshold == 0)
		ssize = size;
	else if (threshold > 0)
		ssize1 = size;
	else
		ssize2 = size;
}

float Grid::FindLineDist(Cell* d, Point point, float threshold){
	//find the minimum distance between point and the segments in the cell d
	//point is given in its grid coordinates but the result is in physical coords

	float mindist = 3.0f, dist, t;
	Point a,b,p,pp;
	float ix = point.x,iy = point.y;
	p = Point(ix,iy);
	int num;
	if (threshold == 0)
		num = d->num_segmentpts/2;
	else
		num = d->num_segmentpts1/2;
	for (int k = 0; k < num; k++){
		//work on each segment
		if (threshold == 0){
			a = d->segments[2*k]; 
			b = d->segments[2*k+1];
		}
		else{
			a = d->segments1[2*k]; 
			b = d->segments1[2*k+1];
		}

		//find the projection p=ta+(1-t)b on [ab]; t = (b-a)(b-p)/(b-a)^2 
		t = ((b.x-a.x)*(b.x-p.x)+(b.y-a.y)*(b.y-p.y))/((b.x-a.x)*(b.x-a.x)+(b.y-a.y)*(b.y-a.y));
		pp = Point(t*a.x+(1-t)*b.x,t*a.y+(1-t)*b.y);
		if ((t < 0) || (t > 1)){//if the projection falls outside [ab]
			dist = float(MIN(sqrt((a.x-ix)*(a.x-ix)+(a.y-iy)*(a.y-iy)),sqrt((b.x-ix)*(b.x-ix)+(b.y-iy)*(b.y-iy))));
		}
		else
			dist = float(sqrt((pp.x-ix)*(pp.x-ix)+(pp.y-iy)*(pp.y-iy)));
		mindist = MIN(mindist,dist);
	}
	mindist = mindist*dx;
	return mindist;
}

void Grid::set_bnd (int b, float ** x ){
	int i;

	for ( i=1 ; i<n-1 ; i++ ) {
		x[0][i] = b==1 ? -x[1][i] : x[1][i];
		x[m-1][i] = b==1 ? -x[m-2][i] : x[m-2][i];
	}
	for ( i = 1; i < m-1; i++){
		x[i][0] = b==2 ? -x[i][1] : x[i][1];
		x[i][n-1] = b==2 ? -x[i][n-2] : x[i][n-2];
	}
	//set the corner values
	x[0][0] = 0.5f*(x[1][0]+x[0][1]);
	x[0][n-1] = 0.5f*(x[1][n-1]+x[0][n-2]);
	x[m-1][0] = 0.5f*(x[m-2][0]+x[m-1][1]);
	x[m-1][n-1] = 0.5f*(x[m-2][n-1]+x[m-1][n-2]);
}

void Grid::advect ( int b, float ** d, float ** d0, float ** u, float ** v, float dt ){
	int i, j, i0, j0, i1, j1;
	float x, y, s0, t0, s1, t1, thr = 0.000001f;
	
	if (bcxlo != 0){
	for (i = 1; i < m-1; i++)
		for (j = 1; j < n-1; j++){
			x = i-dt*u[i][j]/dx; y = j-dt*v[i][j]/dy;
			if (x<thr){ 
				if (bcxlo == 0)
					x = x+(m-1);
				else
					x=thr;
				}
			if (x>=m-1-thr){ 
				if (bcxhi == 0)
					x = x-(m-1);
				else
					x=m-1-thr;
			}
			if (y<thr){ 
				if (bcylo == 0)
					y = y+(n-1);
				else
					y=thr;
			}
			if (y>=n-1-thr){
				if (bcyhi == 0)
					y = y-(n-1);
				else
					y=n-1-thr;
			}
			i0=(int)x; i1=i0+1;
			j0=(int)y; j1=j0+1;
			s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
			d[i][j] = s0*(t0*d0[i0][j0]+t1*d0[i0][j1])+
						 s1*(t0*d0[i1][j0]+t1*d0[i1][j1]);
		//	if ((i0>1)&&(j0>1)&&(i0 < m-3)&&(j0<n-3))
		//		d[i][j] = Interpolate(d0,3,i0,j0,s1,t1);
		//	else
		//		d[i][j] = Interpolate(d0,1,i0,j0,s1,t1);
				
		}
	set_bnd ( b, d );
	}
	else{
		for (i = 0; i < m; i++)
		for (j = 1; j < n-1; j++){
			x = i-dt*u[i][j]/dx; y = j-dt*v[i][j]/dy;
			if (x<thr){ 
				if (bcxlo == 0)
					x = x+(m-1);
				else
					x=thr;
				}
			if (x>=m-1-thr){ 
				if (bcxhi == 0)
					x = x-(m-1);
				else
					x=m-1-thr;
			}
			if (y<thr){ 
				if (bcylo == 0)
					y = y+(n-1);
				else
					y=thr;
			}
			if (y>=n-1-thr){
				if (bcyhi == 0)
					y = y-(n-1);
				else
					y=n-1-thr;
			}
			i0=(int)x; i1=i0+1;
			j0=(int)y; j1=j0+1;
			s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
			d[i][j] = s0*(t0*d0[i0][j0]+t1*d0[i0][j1])+
						 s1*(t0*d0[i1][j0]+t1*d0[i1][j1]);	
		}
		//upper and lower boundaries
		for ( i = 0; i < m; i++){
			d[i][0] = d0[i][1];
			d[i][n-1] = d0[i][n-2];
		}
	}
}

void Grid::update_levelsetSL(float dt){
	int i,j;
	for (i = 0; i < m; i++){
		for (j = 0; j < n; j++)
			temp_phi1[i][j] = phi[i][j];
	}

	//do the advection
	advect(0,phi,temp_phi1,u,v,dt);
	iter++;
}

void Grid::update_levelsetDL(float dt){//DL = Dupont & Liu
	int i,j;
	float extraphi = 0.0f;

	for (i = 0; i < m; i++){
		for (j = 0; j < n; j++){
			temp_phi1[i][j] = phi[i][j];
			temp_phi2[i][j] = phi[i][j];
		}
	}

	//do the advection
	advect(0,temp_phi2,temp_phi1,u,v,dt);
	advect(0,temp_phi1,temp_phi2,u,v,-dt);
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++){
			extraphi = (phi[i][j]-temp_phi1[i][j])/2.0f;
			if (fabs(extraphi) < 0.001)
				extraphi = 0;
			temp_phi1[i][j] = phi[i][j]+extraphi;
			temp_phi2[i][j] = extraphi;
		}
	advect(0,phi,temp_phi1,u,v,dt);
	iter++;
}

void Grid::InterpolateGrid2Cells(void){
	for (int i = 0; i < m-1; i++)
		for (int j = 0; j < n-1; j++)
			cell_list[i][j].phi = (phi[i][j]+phi[i][j+1]+phi[i+1][j]+phi[i+1][j+1])/4;
}

float Grid::Interpolate(float** field, int order, int i, int j, float ldx, float ldy){
		float interp = 0;
		float (*funct)(float x);
		int st;
		
		if (order ==4){
			funct = quartic;
			st = 2;
		}
		if (order ==3){
			funct = cubic;
			st = 1;
		}
		else{
			funct = linear;
			st = 0;
		}
		for (int ii = -st; ii <=1+st; ii++)
			for (int jj = -st; jj <=1+st; jj++)
				interp += field[i+ii][j+jj]*(funct(ldx-ii))*(funct(ldy-jj));
		return interp;
}
