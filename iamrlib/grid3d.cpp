#include "grid3d.H"
#include <math.h>
#include <iostream>
#include <assert.h>
#include <stdlib.h>

using namespace std;

#define SWAP(x0,x) {float *** tmp=x0;x0=x;x=tmp;}
#define FOR_EACH_CELL for ( i=1 ; i<m-1 ; i++ ) { for ( j=1 ; j<n-1 ; j++ ) { for ( k=1 ; k<p-1 ; k++ ) {
#define END_FOR }}}
#define MAX(x,y) ( (x > y) ? (x) : (y) )
#define MIN(x,y) ( (x < y) ? (x) : (y) )
#define SIGN(x) ( (x > 0) ? (1) : ( (x < 0) ? (-1) : (0)) )
#define PI 3.14159f
#define MAX_CELL_PARTICLES 50

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

//These tables are used so that everything can be done in little loops that you can look at all at once
// rather than in pages and pages of unrolled code.

//a2fVertexOffset lists the positions, relative to vertex0, of each of the 8 vertices of a cube
static const float VertexOffset[8][3] =
{
        {0.0, 0.0, 0.0},{1.0, 0.0, 0.0},{1.0, 1.0, 0.0},{0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},{1.0, 0.0, 1.0},{1.0, 1.0, 1.0},{0.0, 1.0, 1.0}
};

//a2iEdgeConnection lists the index of the endpoint vertices for each of the 12 edges of the cube
static const int EdgeConnection[12][2] = 
{
        {0,1}, {1,2}, {2,3}, {3,0},
        {4,5}, {5,6}, {6,7}, {7,4},
        {0,4}, {1,5}, {2,6}, {3,7}
};

//a2fEdgeDirection lists the direction vector (vertex1-vertex0) for each edge in the cube
static const float EdgeDirection[12][3] =
{
        {1.0, 0.0, 0.0},{0.0, 1.0, 0.0},{-1.0, 0.0, 0.0},{0.0, -1.0, 0.0},
        {1.0, 0.0, 0.0},{0.0, 1.0, 0.0},{-1.0, 0.0, 0.0},{0.0, -1.0, 0.0},
        {0.0, 0.0, 1.0},{0.0, 0.0, 1.0},{ 0.0, 0.0, 1.0},{0.0,  0.0, 1.0}
};
int a = 0;

Point& Point::operator = (const Point& pt)
{ x = pt.x; y = pt.y; z = pt.z; return *this; }

Triangle& Triangle::operator = (const Triangle& tr)
{ a = tr.a; b = tr.b; c = tr.c; na = tr.na; nb = tr.nb; nc = tr.nc; return *this; } 

Marker& Marker::operator = (const Marker& mk)
{ p = mk.p; color = mk.color; transp = mk.transp; return *this; } 

Cell::Cell(){}

Cell::Cell(int ii, int jj, int kk, float dx, float dy, float dz){
	i = ii;
	j = jj;
	k = kk;
	phi = 1;
	phi_recon = 3.0f*dx;
	markers = new int[MAX_CELL_PARTICLES];
	triangles = new int[5];
	num_triangles = 0;
	num_particles = 0;
}

Cell::~Cell()
{
}

//constructors
Grid::Grid() {};
Grid::Grid(int mm, int nn, int pp, float xx_lo, float xx_hi, float yy_lo, float yy_hi, float zz_lo, float zz_hi, int bbcxlo, int bbcxhi, int bbcylo, int bbcyhi,int bbczlo, int bbczhi){
	int i,j,k;
	m = mm;
	n = nn;
	p = pp;
	trsize = 0;
	x_lo = xx_lo;
	x_hi = xx_hi;
	y_lo = yy_lo;
	y_hi = yy_hi;
	z_lo = zz_lo;
	z_hi = zz_hi;
	dx = (x_hi-x_lo)/(m-1);
	dy = (y_hi-y_lo)/(n-1);
	dz = (z_hi-z_lo)/(p-1);
	bcxlo = bbcxlo;
	bcxhi = bbcxhi;
	bcylo = bbcylo;
	bcyhi = bbcyhi;
	bczlo = bbczlo;
	bczhi = bbczhi;
	iter = 0;
	
	phi = new float **[m];
	solidphi = new float **[m];
	phi_recon = new float **[m];
	temp_phi = new float **[m];
	den = new float **[m];
	u = new float **[m];
	v = new float **[m];
	w = new float **[m];
	F1 = new float **[m];
	F2 = new float **[m];
	F3 = new float **[m];
	color = new float ***[m];
	tr = new Triangle [5*m*n*3];
	marker_list = new Marker [MAX_CELL_PARTICLES*m*n*3];
	free_markers = new Marker [m*n*4*40];
	num_markers = 0;
	num_free_markers = 0;
	cell_list = new Cell **[m-1];
	for (i = 0; i < m-1; i++){
		cell_list[i] = new Cell *[n-1];
		for (j = 0; j < n-1; j++){
			cell_list[i][j] = new Cell [p-1];
			for (k = 0; k < p-1; k++)	
					cell_list[i][j][k] = Cell(i,j,k,dx,dy,dz);
		}
	}
	for (i = 0; i < m; i++){
	  phi[i] = new float *[n];
	  solidphi[i] = new float *[n];
	  phi_recon[i] = new float *[n];
	  temp_phi[i] = new float *[n];
	  den[i] = new float *[n];
	  u[i] = new float *[n];
	  v[i] = new float *[n];
	  w[i] = new float *[n];
	  F1[i] = new float *[n];
	  F2[i] = new float *[n];
	  F3[i] = new float *[n];
	  color[i] = new float **[n];
	  for (j = 0; j < n; j++){
		  phi[i][j] = new float [p];
		  solidphi[i][j] = new float [p];
		  phi_recon[i][j] = new float [p];
		  temp_phi[i][j] = new float [p];
    		  den[i][j] = new float [p];
		  u[i][j] = new float [p];
		  v[i][j] = new float [p];
		  w[i][j] = new float [p];
		  F1[i][j] = new float [p];
		  F2[i][j] = new float [p];
		  F3[i][j] = new float [p];
		  color[i][j] = new float *[p];
		  for (k = 0; k < p; k++){
			phi[i][j][k] = 1.0f;
			solidphi[i][j][k] = 1.0f;
			phi_recon[i][j][k] = 1.0f;
			temp_phi[i][j][k] = 1.0f;
			den[i][j][k] = 0.001f;
			u[i][j][k] = 0.0f;
			v[i][j][k] = 0.0f;
			w[i][j][k] = 0.0f;
			F1[i][j][k] = 0.0f;
			F2[i][j][k] = 0.0f;
			F3[i][j][k] = 0.0f;
			color[i][j][k] = new float[3];
		  }
	  }
	}
}

float pow2float(float rr) {

 return rr*rr;

}

void Grid::initialize_problem(char * problem_name){
	int i,j,k;
	float radius = 6, dlt = 0.25f, left = 47.4999f, right = 52.5001f, notch = 72.5001f, maxval = 0;
	double d1;
	float ***temp_phi = new float **[m];

	if (problem_name == "zalesak"){	
		float xc = 50.0f, yc = 75.000001f, zc = 50.0f, r = 15.0f;
		//initialize notched circle
		for (i = 0; i < m; i++){
			for ( j = 0; j < n; j++)
				for (k = 0; k < p; k++){

				d1=sqrt(pow2float(x_lo+i*dx-xc)+
                                        pow2float(y_lo+j*dy-yc)+
                                        pow2float(z_lo+k*dz-zc))-r;
				
				if ((x_lo+i*dx>=left)&&(x_lo+i*dx<=right)){
				  if (y_lo+j*dy<=60.0)
				   d1=( (x_lo+i*dx<xc) ? sqrt(pow2float(y_lo+j*dy-60.0)+pow2float(x_lo+i*dx-left)) : sqrt(pow2float(y_lo+j*dy-60.0)+pow2float(x_lo+i*dx-right)) );
				  else if (y_lo+j*dy<=notch)
				   d1=MIN( ( (x_lo+i*dx<50.0) ? (x_lo+i*dx-left) : (right-(x_lo+i*dx)) ) , notch-(y_lo+j*dy) );
				  else if ((y_lo+j*dy<=90.0)&&(d1<=0.0))
				   d1=MAX(d1,notch-(y_lo+j*dy));
				 }
				else if ((d1<0.0)&&(x_lo+i*dx<left)){
				  if (y_lo+j*dy<=notch)
				   d1=MAX(d1,(x_lo+i*dx-left));
				  else
				   d1=MAX(d1,-sqrt(pow2float(x_lo+i*dx-left)+pow2float(y_lo+j*dy-notch)) );
				 }
				else if ((d1<0.0)&&(x_lo+i*dx>right)){
				  if (y_lo+j*dy<=notch)
				   d1=MAX(d1,(right-(x_lo+i*dx)));
				  else
				   d1=MAX(d1,-sqrt(pow2float(x_lo+i*dx-right)+pow2float(y_lo+j*dy-notch)) );
				 }
				
				phi[i][j][k] = -float(d1);
				color[i][j][k][0] = 1.0f; color[i][j][k][1] =1.0f; color[i][j][k][2] = 0.0f;

			}
		}

		//initialize velocity
		for (i = 0; i < m; i++)
		  for (j = 0; j < n; j++)
			for (k = 0; k < p; k++){
			  u[i][j][k] = PI/314.0f*(50-(y_lo+j*dy));
			  v[i][j][k] = PI/314.0f*(x_lo+i*dx-50);
			  w[i][j][k] = 0.0f;
		  }
	}
	else if (problem_name == "vortex in a box"){
		for (i = 0; i < m; i++)
			for ( j = 0; j < n; j++)
				for (k = 0; k < p; k++){
				phi[i][j][k] = -float((sqrt(pow2float(x_lo+i*dx-0.35)+pow2float(y_lo+j*dy-0.35)+pow2float(z_lo+k*dz-0.35))-0.15));
				color[i][j][k][0] = 1.0f; color[i][j][k][1] =1.0f; color[i][j][k][2] = 0.0f;
				}
		//initialize velocity
		for (i = 0; i < m; i++)
		  for (j = 0; j < n; j++)
			for (k = 0; k < p; k++){
				u[i][j][k] = float(2*sin(2*PI*(y_lo+j*dy))*sin(2*PI*(z_lo+k*dz))*pow2float(sin(PI*(x_lo+i*dx))));
				v[i][j][k] = -float(sin(2*PI*(x_lo+i*dx))*sin(2*PI*(z_lo+k*dz))*pow2float(sin(PI*(y_lo+j*dy))));
			    w[i][j][k] = -float(sin(2*PI*(x_lo+i*dx))*sin(2*PI*(y_lo+j*dy))*pow2float(sin(PI*(z_lo+k*dz))));
		  }	
	}
	else if (problem_name == ""){
		for (i=0 ; i<m ; i++) 
				for (j = 0; j < n; j++)
					for (k = 0; k < p; k++)
						phi[i][j][k] = float(-(sqrt(pow2float(i-8)+pow2float(j-(n-1)/2.0f)+pow2float(k-(p-1)/2.0f))-radius));
	}
}


void Grid::CreateSphere(float xc, float yc, float zc, float radius){
	int i,j,k;
	float dist;
	for (i = 0; i < m; i++)
		for ( j = 0; j < n; j++)
			for (k = 0; k < p; k++){
			dist = float((sqrt(pow2float(x_lo+i*dx-xc)+pow2float(y_lo+j*dy-yc)+pow2float(z_lo+k*dz-zc))-radius));
			if (phi[i][j][k] > dist)
				phi[i][j][k] = dist;
			}
}

bool Grid::InsideDomain(Point* pp){
	bool inside_domain;
	float dx_2=dx/2;

	if (((pp->x-x_lo)>=dx_2)&&((x_hi-pp->x)>=dx_2)&&((pp->y-y_lo)>=dx_2)&&
		((y_hi-pp->y)>=dx_2)&&((pp->z-z_lo)>=dx_2)&&((z_hi-pp->z)>=dx_2))
		inside_domain = true;
	else
		inside_domain = false;
	
	return inside_domain;
}

void Grid::InitializeParticles(void){
	int i,j,kk,k,li,lj,np,index;
	Cell* cl;
	Point pp,ctr;
	vect a,b,c;
	double area,la,lb,lc,sp;
	float t,dist,dist0,eps=dx/16.0,dx_2=dx/2;

	for (i = 0; i < m-1; i++)
		for (j = 0; j < n-1; j++)
			for (k = 0; k < p-1; k++){
				cl = &cell_list[i][j][k];
				if (CellChangesSignInside(cl, 0)){
					index = 0;
					for (kk = 0; kk < cl->num_triangles; kk++){//for each triangle
						a = tr[cl->triangles[kk]].a;
						b = tr[cl->triangles[kk]].b;
						c = tr[cl->triangles[kk]].c;
						la = sqrt(pow2float(a.fX-b.fX)+pow2float(a.fY-b.fY)+pow2float(a.fZ-b.fZ));
						lb = sqrt(pow2float(a.fX-c.fX)+pow2float(a.fY-c.fY)+pow2float(a.fZ-c.fZ));
						lc = sqrt(pow2float(c.fX-b.fX)+pow2float(c.fY-b.fY)+pow2float(c.fZ-b.fZ));
						//compute the semiperimeter
						sp = (la+lb+lc)/2;
						//and area
						area = sqrt(sp*(sp-la)*(sp-lb)*(sp-lc));
						//estimate number of particles (on any triangle edge) 
						np = int(floor(area/0.11f));
						if (np < 2)
							np = 2;
						//we place np*(np+1)/2 particles
						for (li = 0; li <= np; li++)
							for (lj = 0; lj <= li; lj++){
								t = li/float(np);
								pp = Point(x_lo+dx*(a.fX+li*(b.fX-a.fX)/np+lj*(c.fX-b.fX)/(li+1)),
									y_lo+dy*(a.fY+li*(b.fY-a.fY)/np+lj*(c.fY-b.fY)/(li+1)),
									z_lo+dz*(a.fZ+li*(b.fZ-a.fZ)/np+lj*(c.fZ-b.fZ)/(li+1)));
								ctr = Point(x_lo+dx*(a.fX+b.fX+c.fX)/3.0,y_lo+dx*(a.fY+b.fY+c.fY)/3.0,z_lo+dx*(a.fZ+b.fZ+c.fZ)/3.0);
								pp.x += eps*(ctr.x-pp.x);
								pp.y += eps*(ctr.y-pp.y);
								pp.z += eps*(ctr.z-pp.z);
								if ((cl->num_particles < MAX_CELL_PARTICLES-16)&&(InsideDomain(&pp))){
									marker_list[num_markers].p = pp;
									cl->markers[cl->num_particles++] = num_markers;									
									
									dist0 = 3*dx;
									for (int ix = i; ix <= i+1; ix++)
										for (int iy = j; iy <= j+1; iy++)
											for (int iz = k; iz <= k+1; iz++){
												dist = float(sqrt((pp.x-ix*dx-x_lo)*(pp.x-ix*dx-x_lo)+
													(pp.y-iy*dy-y_lo)*(pp.y-iy*dy-y_lo)+(pp.z-iz*dz-z_lo)*(pp.z-iz*dz-z_lo)));
												if (dist < dist0){
													dist0 = dist;
													marker_list[num_markers].color = vect(color[ix][iy][iz][0],color[ix][iy][iz][1],color[ix][iy][iz][2]);

												}
											}
									num_markers++;
								}
							}
					}
				}	
			}
	int aa = 0;
}

void Grid::AssignParticleTexture(int type, char * pxls, int usize, int vsize){
	int i, ui, vi, Rc, Gc, Bc; 
	float uc, vc, cu;
	Marker* mk;
	//for spheres I use (x,y,z) = (sinu*cosv,sinu*sinv,cosu) with u in (-pi/2,pi/2) v in (0,2pi)
	//so that u = 1/R*acos(z-zc) and v = atan(y-yc/x-xc)
	if (type == 1){
		for (i = 0; i < num_markers; i++){
			mk = &marker_list[i];
			mk->transp = 0;
		//	if (mk->p.z < 2.6){//this is usual
			if (mk->p.z > 2.0){//or < 2.6
				mk->color.fX = 0.85;
				mk->color.fY = 0.95;
				mk->color.fZ = 1.0;
			}
			else {
				//cu = (mk->p.z-7.6)/0.3;//the usual
				//cu = (mk->p.z-4.6)/1.0;
				cu = (mk->p.z-5.1)/2.0;
				if (cu > 0)
					cu = MIN(cu,1);
				else if (cu < 0)
					cu = MAX(cu,-1);
				uc = acos(cu);
				vc = atan2(mk->p.y-1.55,mk->p.x-1.45);//milkdrop
				//vc = atan2(mk->p.y-5.0,mk->p.x-5.0);//coffee
				//vc = atan2(mk->p.y-2.5,mk->p.x-5.0);//dam
				uc = PI-uc;
				vc = vc+PI;
				if (vc < 0)
					vc = 0;
				if (vc > 2*PI)
					vc = 2*PI;
				ui = int((usize-1)*uc/PI);
				vi = int((vsize-1)*vc/2/PI);
				Rc = int(pxls[3*(ui*vsize+vi)+0]);
				Gc = int(pxls[3*(ui*vsize+vi)+1]); 
				Bc = int(pxls[3*(ui*vsize+vi)+2]);
				if (Rc<0)
					Rc +=255;
				if (Gc<0)
					Gc += 255;
				if (Bc<0)
					Bc += 255;
				mk->color.fX = Rc/255.0;
				mk->color.fY = Gc/255.0;
				mk->color.fZ = Bc/255.0;
			}
		}
	}
	else if (type == 3){
		for (i = 0; i < num_markers; i++){
			mk = &marker_list[i];
			if (mk->p.z > 1.5){
				mk->color.fX = 0.85;
				mk->color.fY = 0.95;
				mk->color.fZ = 1.0;
			}
			else{
				ui = int((usize-1)*(mk->p.x-x_lo)/(x_hi-x_lo));
				if (ui>=(usize-1)) ui = usize -2;	
				vi = int((vsize-1)*(mk->p.y-y_lo)/(y_hi-y_lo));
				if (vi>=(vsize-1)) vi = vsize -2;
				//std::cout << ui << " " << vi << '\n';
				Rc = int(pxls[3*(ui*vsize+vi)+0]);
				Gc = int(pxls[3*(ui*vsize+vi)+1]); 
				Bc = int(pxls[3*(ui*vsize+vi)+2]);
				if (Rc<0)
					Rc +=255;
				if (Gc<0)
					Gc += 255;
				if (Bc<0)
					Bc += 255;
				mk->color.fX = Rc/255.0;
				mk->color.fY = Gc/255.0;
				mk->color.fZ = Bc/255.0;
			}
		}
	}
}

void Grid::OutputMarkers(char *filename){

	Cell *c;
	FILE *fid = fopen (filename, "w+t");
	for (int i = 0; i < m-1; i++)
		for (int j = 0; j < n-1; j++)
			for (int k = 0; k < p-1; k++){
				c = &cell_list[i][j][k];
				if ( c->num_particles >= 0){
					fprintf(fid, "%d ", c->num_particles);
					for (int kk=0; kk < c->num_particles; kk++)	
						fprintf(fid, "%f %f %f ", marker_list[c->markers[kk]].p.x,marker_list[c->markers[kk]].p.y,marker_list[c->markers[kk]].p.z);		
					for (int kk=0; kk < c->num_particles; kk++)	
						fprintf(fid, "%f %f %f ", marker_list[c->markers[kk]].color.fX,marker_list[c->markers[kk]].color.fY,marker_list[c->markers[kk]].color.fZ);
					fprintf(fid, "\n");
				}
			}
	fclose(fid);
}

void Grid::OutputFreeMarkers(char *filename){

	FILE *fid = fopen (filename, "w+t");
	fprintf(fid, "%d ", num_free_markers);
	fprintf(fid, "\n");
	for (int i = 0; i < num_free_markers; i++){
		fprintf(fid, "%f %f %f %f", free_markers[i].p.x, free_markers[i].p.y, free_markers[i].p.z, free_markers[i].transp);		
		fprintf(fid, "\n");
	}
	fclose(fid);
}

void Grid::OutputLevelSet(char *filename){

	FILE *fid = fopen (filename, "w+t");
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < p; k++)
				fprintf(fid, "%f\n", phi[i][j][k]);	
	fclose(fid);
}

void Grid::OutputSolidLevelSet(char *filename){

 FILE *fid = fopen (filename, "w+t");
 for (int i = 0; i < m; i++)
 for (int j = 0; j < n; j++)
 for (int k = 0; k < p; k++)
  fprintf(fid, "%f\n", solidphi[i][j][k]);	
 fclose(fid);
}

void Grid::OutputDensity(char *filename){

 FILE *fid = fopen (filename, "w+t");
 for (int i = 0; i < m; i++)
 for (int j = 0; j < n; j++)
 for (int k = 0; k < p; k++)
  fprintf(fid, "%f\n", den[i][j][k]);	
 fclose(fid);
}

void Grid::OutputWeber(char *filename){
 Cell* c;
 float vel_var, curv;
 FILE *fid = fopen (filename, "w+t");
 for (int i = 0; i < m-1; i++)
 for (int j = 0; j < n-1; j++)
 for (int k = 0; k < p-1; k++){
  c = &cell_list[i][j][k];
  if (NotBoundaryCell(c)){
   vel_var = VelocityVariation(c);
   curv = AbsCurvature(c);
   fprintf(fid, "%f %f\n",vel_var,curv);	
  }
 }
 fclose(fid);
}

void Grid::OutputVelocity(char *filename){

 FILE *fid = fopen (filename, "w+t");
 for (int i = 0; i < m; i++)
 for (int j = 0; j < n; j++)
 for (int k = 0; k < p; k++)
  fprintf(fid, "%f %f %f\n", u[i][j][k], v[i][j][k], w[i][j][k]);	
 fclose(fid);
}

void Grid::OutputVelocityMag(char *filename){
	
 float mag;
 FILE *fid = fopen (filename, "w+t");
 for (int i = 0; i < m; i++)
 for (int j = 0; j < n; j++)
 for (int k = 0; k < p; k++){
  mag = float(sqrt(pow2float(u[i][j][k])+ pow2float(v[i][j][k])+ 
        pow2float(w[i][j][k])));
  fprintf(fid, "%f\n", mag);	
 }
 fclose(fid);
}

void Grid::OutputFVelocity(char *filename){

 FILE *fid = fopen (filename, "w+t");
 for (int i = 0; i < m; i++)
 for (int j = 0; j < n; j++)
 for (int k = 0; k < p; k++)
  fprintf(fid, "%f %f %f\n", u[i][j][k], v[i][j][k], w[i][j][k]);	
 fclose(fid);
}

void Grid::OutputColor(char *filename){

 FILE *fid = fopen (filename, "w+t");
 for (int i = 0; i < m; i++)
 for (int j = 0; j < n; j++)
 for (int k = 0; k < p; k++)
  fprintf(fid, "%f %f %f\n", color[i][j][k][0],color[i][j][k][1],
     color[i][j][k][2]);	
 fclose(fid);
}

void Grid::StopTopVel(void){
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			for (int k = p-2; k < p; k++){
				u[i][j][k] = 0;
				v[i][j][k] = 0;
			}
}
vect Grid::InterpVelocity(Point pp){
	int ix = int((pp.x-x_lo)/dx), iy = int((pp.y-y_lo)/dy), iz = int((pp.z-z_lo)/dz);
	float u_vel=0,v_vel=0,w_vel=0;
	if ((ix < 0 )||(ix >= m-1)||(iy < 0 )||(iy >= n-1)||(iz < 0 )||(iz >= p-1))
		return vect(0,0,0);	
	float ldx = (pp.x-x_lo)/dx-ix, ldy = (pp.y-y_lo)/dy-iy, ldz = (pp.z-z_lo)/dz-iz, x, y, z;
	//cout << " ix iy iz " << ix << " " << iy << " " << iz << '/n';  
	u_vel = (1-ldz)*((1-ldy)*((1-ldx)*u[ix][iy][iz]+ldx*u[ix+1][iy][iz])+ldy*((1-ldx)*u[ix][iy+1][iz]+ldx*u[ix+1][iy+1][iz]))+
		ldz*((1-ldy)*((1-ldx)*u[ix][iy][iz+1]+ldx*u[ix+1][iy][iz+1])+ldy*((1-ldx)*u[ix][iy+1][iz+1]+ldx*u[ix+1][iy+1][iz+1]));
	v_vel = (1-ldz)*((1-ldy)*((1-ldx)*v[ix][iy][iz]+ldx*v[ix+1][iy][iz])+ldy*((1-ldx)*v[ix][iy+1][iz]+ldx*v[ix+1][iy+1][iz]))+
		(ldz)*((1-ldy)*((1-ldx)*v[ix][iy][iz+1]+ldx*v[ix+1][iy][iz+1])+ldy*((1-ldx)*v[ix][iy+1][iz+1]+ldx*v[ix+1][iy+1][iz+1]));
	w_vel = (1-ldz)*((1-ldy)*((1-ldx)*w[ix][iy][iz]+ldx*w[ix+1][iy][iz])+ldy*((1-ldx)*w[ix][iy+1][iz]+ldx*w[ix+1][iy+1][iz]))+
		(ldz)*((1-ldy)*((1-ldx)*w[ix][iy][iz+1]+ldx*w[ix+1][iy][iz+1])+ldy*((1-ldx)*w[ix][iy+1][iz+1]+ldx*w[ix+1][iy+1][iz+1]));
	return vect(u_vel,v_vel,w_vel);	
}
void Grid::AdvectMarkers(float dt){

	int dropped_markers = 0,num_newmarkers = num_markers;

	for (int i = 0; i < num_markers; i++){
		if (!AdvectMarker(i,dropped_markers,dt)){
			num_newmarkers--;
			dropped_markers++;
		}
	}
	num_markers = num_newmarkers;
}

void Grid::AdvectFreeMarkers(float dt){

	int dropped_markers = 0,num_newmarkers = num_free_markers;

	for (int i = 0; i < num_free_markers; i++){
		if (!AdvectFreeMarker(i,dropped_markers,dt)){
			num_newmarkers--;
			dropped_markers++;
		}
	}
	num_free_markers = num_newmarkers;
	//num_free_markers = 0;
}

bool Grid::AdvectMarker(int i, int dropped, float dt){
	
	Point pp = marker_list[i].p;
	int ix = int((pp.x-x_lo)/dx), iy = int((pp.y-y_lo)/dy), iz = int((pp.z-z_lo)/dz);
	if ((ix < 0 )||(ix >= m-1)||(iy < 0 )||(iy >= n-1)||(iz < 0 )||(iz >= p-1))
		return false;	
	float ldx = (pp.x-x_lo)/dx-ix, ldy = (pp.y-y_lo)/dy-iy, ldz = (pp.z-z_lo)/dz-iz, x, y, z, dx_2=dx/2;

	//interpolate the velocity at the marker position (linear interp)
	float u_vel = (1-ldz)*((1-ldy)*((1-ldx)*u[ix][iy][iz]+ldx*u[ix+1][iy][iz])+ldy*((1-ldx)*u[ix][iy+1][iz]+ldx*u[ix+1][iy+1][iz]))+
		ldz*((1-ldy)*((1-ldx)*u[ix][iy][iz+1]+ldx*u[ix+1][iy][iz+1])+ldy*((1-ldx)*u[ix][iy+1][iz+1]+ldx*u[ix+1][iy+1][iz+1]));
	float v_vel = (1-ldz)*((1-ldy)*((1-ldx)*v[ix][iy][iz]+ldx*v[ix+1][iy][iz])+ldy*((1-ldx)*v[ix][iy+1][iz]+ldx*v[ix+1][iy+1][iz]))+
		(ldz)*((1-ldy)*((1-ldx)*v[ix][iy][iz+1]+ldx*v[ix+1][iy][iz+1])+ldy*((1-ldx)*v[ix][iy+1][iz+1]+ldx*v[ix+1][iy+1][iz+1]));
	float w_vel = (1-ldz)*((1-ldy)*((1-ldx)*w[ix][iy][iz]+ldx*w[ix+1][iy][iz])+ldy*((1-ldx)*w[ix][iy+1][iz]+ldx*w[ix+1][iy+1][iz]))+
		(ldz)*((1-ldy)*((1-ldx)*w[ix][iy][iz+1]+ldx*w[ix+1][iy][iz+1])+ldy*((1-ldx)*w[ix][iy+1][iz+1]+ldx*w[ix+1][iy+1][iz+1]));
	
	//compute the new marker position
	float px = pp.x+dt*u_vel, py = pp.y+dt*v_vel, pz = pp.z+dt*w_vel;	
	float px1 = (px-x_lo)/dx, py1 = (py-y_lo)/dy, pz1 = (pz-z_lo)/dz;
	//if ((ix1 < 0 )||(ix1 >= m-1)||(iy1 < 0 )||(iy1 >= n-1)||(iz1 < 0 )||(iz1 >= p-1))//out of the domain
	//	return false;
	//if (ix1 < 0 ) ix1 = 0; if (ix1 > m-2) ix1 = m-2;	
	//if (iy1 < 0 ) iy1 = 0; if (iy1 > n-2) iy1 = n-2;
	//if (iz1 < 0 ) iz1 = 0; if (iz1 > p-2) iz1 = p-2;
	if (px1<0.5f){ 
		if (bcxlo == 0)
			px1 = px1+(m-1);
		else
			return false;
	}
	if (px1>=m-1-0.5f){ 
		if (bcxhi == 0)
			px1 = px1-(m-1);
		else
			return false;
	}
	if (py1<0.5f){ 
		if (bcylo == 0)
			py1 = py1+(n-1);	
		else
			return false; 
	}
	if (py1>=n-1-0.5f){
		if (bcyhi == 0)
			py1 = py1-(n-1);
		else
			return false;
	}
	if (pz1<0.5f){ 
		if (bczlo == 0)
			pz1 = pz1+(p-1);	
		else
			return false; 
	}
	if (pz1>=p-1-0.5f){
		if (bczhi == 0)
			pz1 = pz1-(p-1);
		else
			return false;
	}
	int ix1 = int(px1), iy1 = int(py1), iz1 = int(pz1);
	//int ix1 = int((px-x_lo)/dx), iy1 = int((py-y_lo)/dy), iz1 = int((pz-z_lo)/dz);
	if ((ix1 < 0 )||(ix1 >= m-1)||(iy1 < 0 )||(iy1 >= n-1)||(iz1 < 0 )||(iz1 >= p-1))//out of the domain
		return false;
	float ldx1 = (px-x_lo)/dx-ix1, ldy1 = (py-y_lo)/dy-iy1, ldz1 = (pz-z_lo)/dz-iz1;
	float u_vel1 = (1-ldz1)*((1-ldy1)*((1-ldx1)*u[ix1][iy1][iz1]+ldx1*u[ix1+1][iy1][iz1])+ldy1*((1-ldx1)*u[ix1][iy1+1][iz1]+ldx1*u[ix1+1][iy1+1][iz1]))+
		(ldz1)*((1-ldy1)*((1-ldx1)*u[ix1][iy1][iz1+1]+ldx1*u[ix1+1][iy1][iz1+1])+ldy1*((1-ldx1)*u[ix1][iy1+1][iz1+1]+ldx1*u[ix1+1][iy1+1][iz1+1]));
	float v_vel1 = (1-ldz1)*((1-ldy1)*((1-ldx1)*v[ix1][iy1][iz1]+ldx1*v[ix1+1][iy1][iz1])+ldy1*((1-ldx1)*v[ix1][iy1+1][iz1]+ldx1*v[ix1+1][iy1+1][iz1]))+
		(ldz1)*((1-ldy1)*((1-ldx1)*v[ix1][iy1][iz1+1]+ldx1*v[ix1+1][iy1][iz1+1])+ldy1*((1-ldx1)*v[ix1][iy1+1][iz1+1]+ldx1*v[ix1+1][iy1+1][iz1+1]));
	float w_vel1 = (1-ldz1)*((1-ldy1)*((1-ldx1)*w[ix1][iy1][iz1]+ldx1*w[ix1+1][iy1][iz1])+ldy1*((1-ldx1)*w[ix1][iy1+1][iz1]+ldx1*w[ix1+1][iy1+1][iz1]))+
		(ldz1)*((1-ldy1)*((1-ldx1)*w[ix1][iy1][iz1+1]+ldx1*w[ix1+1][iy1][iz1+1])+ldy1*((1-ldx1)*w[ix1][iy1+1][iz1+1]+ldx1*w[ix1+1][iy1+1][iz1+1]));
	px = pp.x + 0.5f*(u_vel+u_vel1)*dt;
	py = pp.y + 0.5f*(v_vel+v_vel1)*dt;
	pz = pp.z + 0.5f*(w_vel+w_vel1)*dt;
	//this was Crank-Nicholson (not Runge Kutta 2)...?
	x = (px-x_lo)/dx; y = (py-y_lo)/dy; z = (pz-z_lo)/dz;
	if (x<0.5f){ 
		if (bcxlo == 0)
			x = x+(m-1);
		else
			return false;
	}
	if (x>=m-1-0.5f){ 
		if (bcxhi == 0)
			x = x-(m-1);
		else
			return false;
	}
	if (y<0.5f){ 
		if (bcylo == 0)
			y = y+(n-1);	
		else
			return false; 
	}
	if (y>=n-1-0.5f){
		if (bcyhi == 0)
			y = y-(n-1);
		else
			return false;
	}
	if (z<0.5f){ 
		if (bczlo == 0)
			z = z+(p-1);	
		else
			return false; 
	}
	if (z>=p-1-0.5f){
		if (bczhi == 0)
			z = z-(p-1);
		else
			return false;
	}
	marker_list[i-dropped].p = Point(px,py,pz);
	marker_list[i-dropped].color = marker_list[i].color;
	return true;
}

bool Grid::AdvectFreeMarker(int i, int dropped, float dt){
	//physical parameters are in cm/g/s
	float gravity = -980, density = 1, radiusbub = 2.5e-4, alphadt = 1e-4*dt, oneovervisc=1e6;
	float phi_marker, newvelx, newvely, newvelz, px, py, pz, uinit = 0, vinit=0, winit=0;
	vect newvel;
	Point pp = free_markers[i].p;
	int ix = int((pp.x-x_lo)/dx), iy = int((pp.y-y_lo)/dy), iz = int((pp.z-z_lo)/dz);
	if ((ix < 0 )||(ix >= m-1)||(iy < 0 )||(iy >= n-1)||(iz < 0 )||(iz >= p-1))
		return false;	
	
	float ldx = (pp.x-x_lo)/dx-ix, ldy = (pp.y-y_lo)/dy-iy, ldz = (pp.z-z_lo)/dz-iz, x, y, z;

	//scale radius (and mass) with random size
	radiusbub *= fabs(free_markers[i].transp);
	float mass = density*radiusbub;

	//interpolate the velocity at the marker position (linear interp)
	float u_vel = (1-ldz)*((1-ldy)*((1-ldx)*F1[ix][iy][iz]+ldx*F1[ix+1][iy][iz])+ldy*((1-ldx)*F1[ix][iy+1][iz]+ldx*F1[ix+1][iy+1][iz]))+
		ldz*((1-ldy)*((1-ldx)*F1[ix][iy][iz+1]+ldx*F1[ix+1][iy][iz+1])+ldy*((1-ldx)*F1[ix][iy+1][iz+1]+ldx*F1[ix+1][iy+1][iz+1]));
	float v_vel = (1-ldz)*((1-ldy)*((1-ldx)*F2[ix][iy][iz]+ldx*F2[ix+1][iy][iz])+ldy*((1-ldx)*F2[ix][iy+1][iz]+ldx*F2[ix+1][iy+1][iz]))+
		(ldz)*((1-ldy)*((1-ldx)*F2[ix][iy][iz+1]+ldx*F2[ix+1][iy][iz+1])+ldy*((1-ldx)*F2[ix][iy+1][iz+1]+ldx*F2[ix+1][iy+1][iz+1]));
	float w_vel = (1-ldz)*((1-ldy)*((1-ldx)*F3[ix][iy][iz]+ldx*F3[ix+1][iy][iz])+ldy*((1-ldx)*F3[ix][iy+1][iz]+ldx*F3[ix+1][iy+1][iz]))+
		(ldz)*((1-ldy)*((1-ldx)*F3[ix][iy][iz+1]+ldx*F3[ix+1][iy][iz+1])+ldy*((1-ldx)*F3[ix][iy+1][iz+1]+ldx*F3[ix+1][iy+1][iz+1]));
	
	//initial velocity 
	if (free_markers[i].transp > 0){//droplet, use own velocity
		uinit = free_markers[i].color.fX;
		vinit = free_markers[i].color.fY;
		winit = free_markers[i].color.fZ;
	}
	else{//bubble, use liquid velocity
		uinit = u_vel;
		vinit = v_vel;
		winit = w_vel;
	}
	
	//compute (intermediate) new marker  position 
	px = pp.x+dt*uinit;
	py = pp.y+dt*vinit;
	pz = pp.z+dt*winit;

	//check for disappearance
	int ix1 = int((px-x_lo)/dx), iy1 = int((py-y_lo)/dy), iz1 = int((pz-z_lo)/dz);
	if ((ix1 < 0 )||(ix1 > m-2)||(iy1 < 0 )||(iy1 > n-2)||(iz1 < 0 )||(iz1 > p-2))//out of the domain
		return false;
	
	float ldx1 = (px-x_lo)/dx-ix1, ldy1 = (py-y_lo)/dy-iy1, ldz1 = (pz-z_lo)/dz-iz1;
	float u_vel1 = (1-ldz1)*((1-ldy1)*((1-ldx1)*F1[ix1][iy1][iz1]+ldx1*F1[ix1+1][iy1][iz1])+ldy1*((1-ldx1)*F1[ix1][iy1+1][iz1]+ldx1*F1[ix1+1][iy1+1][iz1]))+
		(ldz1)*((1-ldy1)*((1-ldx1)*F1[ix1][iy1][iz1+1]+ldx1*F1[ix1+1][iy1][iz1+1])+ldy1*((1-ldx1)*F1[ix1][iy1+1][iz1+1]+ldx1*F1[ix1+1][iy1+1][iz1+1]));
	float v_vel1 = (1-ldz1)*((1-ldy1)*((1-ldx1)*F2[ix1][iy1][iz1]+ldx1*F2[ix1+1][iy1][iz1])+ldy1*((1-ldx1)*F2[ix1][iy1+1][iz1]+ldx1*F2[ix1+1][iy1+1][iz1]))+
		(ldz1)*((1-ldy1)*((1-ldx1)*F2[ix1][iy1][iz1+1]+ldx1*F2[ix1+1][iy1][iz1+1])+ldy1*((1-ldx1)*F2[ix1][iy1+1][iz1+1]+ldx1*F2[ix1+1][iy1+1][iz1+1]));
	float w_vel1 = (1-ldz1)*((1-ldy1)*((1-ldx1)*F3[ix1][iy1][iz1]+ldx1*F3[ix1+1][iy1][iz1])+ldy1*((1-ldx1)*F3[ix1][iy1+1][iz1]+ldx1*w[ix1+1][iy1+1][iz1]))+
		(ldz1)*((1-ldy1)*((1-ldx1)*F3[ix1][iy1][iz1+1]+ldx1*F3[ix1+1][iy1][iz1+1])+ldy1*((1-ldx1)*F3[ix1][iy1+1][iz1+1]+ldx1*F3[ix1+1][iy1+1][iz1+1]));
	
	//update gravity and drag terms
	if (free_markers[i].transp > 0){//droplets
		newvelx = (alphadt*0.5f*(u_vel+u_vel1)+free_markers[i].color.fX*(mass-0.5f*alphadt))/(mass+0.5f*alphadt);
		newvely = (alphadt*0.5f*(v_vel+v_vel1)+free_markers[i].color.fY*(mass-0.5f*alphadt))/(mass+0.5f*alphadt);
		//newvelz = (alphadt*0.5f*(w_vel+w_vel1)+free_markers[i].color.fZ*(mass-0.5f*alphadt)+mass*gravity*dt)/(mass+0.5f*alphadt);
		float xcoeff = 1.0f;//was 0.5 for dog, etc...? 0.25 for bubbles?
		newvelz = (alphadt*0.5f*(w_vel+w_vel1)+free_markers[i].color.fZ*(mass-0.5f*alphadt)+mass*gravity*dt*xcoeff)/(mass+0.5f*alphadt);
	}
	else{//bubbles
		newvelx = u_vel1;
		newvely = v_vel1;
		newvelz = w_vel1-4.0/9.0f*gravity*radiusbub*radiusbub*oneovervisc;
	}
	
	//update final position
	px = pp.x + 0.5f*(newvelx+uinit)*dt;
	py = pp.y + 0.5f*(newvely+vinit)*dt;
	pz = pp.z + 0.5f*(newvelz+winit)*dt;

	x = (px-x_lo)/dx; y = (py-y_lo)/dy; z = (pz-z_lo)/dz;
	if (x<0.5f){ 
		if (bcxlo == 0)
			x = x+(m-1);
		else
			return false;
	}
	if (x>=m-1-0.5f){ 
		if (bcxhi == 0)
			x = x-(m-1);
		else
			return false;
	}
	if (y<0.5f){ 
		if (bcylo == 0)
			y = y+(n-1);	
		else
			return false; 
	}
	if (y>=n-1-0.5f){
		if (bcyhi == 0)
			y = y-(n-1);
		else
			return false;
	}
	if (z<0.5f){ 
		if (bczlo == 0)
			z = z+(p-1);	
		else
			return false; 
	}
	if (z>=p-1-0.5f){
		if (bczhi == 0)
			z = z-(p-1);
		else
			return false;
	}

	//check if inside solid
	//if (InsideSolid(&cell_list[(int)x][(int)y][(int)z],0))
	//	BounceMarker(px,py,pz,newvelx,newvely,newvelz);
	
	x = (px-x_lo)/dx; y = (py-y_lo)/dy; z = (pz-z_lo)/dz;
	if (InsideSolid(&cell_list[(int)x][(int)y][(int)z],0))
		return false;

	//check for phase change
	//interpolate the level set value at the marker position (linear interp)
	ix = (int)x; iy = (int)y; iz = (int)z;
	ldx = x-ix; ldy = y-iy; ldz = z-iz;
	phi_marker= (1-ldz)*((1-ldy)*((1-ldx)*phi[ix][iy][iz]+ldx*phi[ix+1][iy][iz])+ldy*((1-ldx)*phi[ix][iy+1][iz]+ldx*phi[ix+1][iy+1][iz]))+
		ldz*((1-ldy)*((1-ldx)*phi[ix][iy][iz+1]+ldx*phi[ix+1][iy][iz+1])+ldy*((1-ldx)*phi[ix][iy+1][iz+1]+ldx*phi[ix+1][iy+1][iz+1]));
	if (phi_marker*(free_markers[i].transp)>0)
		return false;
	free_markers[i-dropped].p = Point(px,py,pz);
	free_markers[i-dropped].color = vect(newvelx,newvely,newvelz);
	free_markers[i-dropped].transp = free_markers[i].transp;
	return true;
}

void Grid::BounceMarker(float& px, float& py, float& pz, float& newvelx, float& newvely, float& newvelz){
	float x = (px-x_lo)/dx, y = (py-y_lo)/dy, z = (pz-z_lo)/dz;
	int i = x, j = y, k = z;
	float gradx=(solidphi[i+1][j][k]-solidphi[i-1][j][k])/(2*dx), grady = (solidphi[i][j+1][k]-solidphi[i][j-1][k])/(2*dy), 
	gradz = (solidphi[i][j][k+1]-solidphi[i][j][k-1])/(2*dz);
	float grad = sqrt(gradx*gradx+grady*grady+gradz*gradz)+1e-8;
	float scalprodovergrad = (newvelx*gradx+newvely*grad+newvelz*gradz)/(grad*grad);
	px = px + newvelx-2*scalprodovergrad*gradx;
	py = py + newvely-2*scalprodovergrad*grady;
	pz = pz + newvelz-2*scalprodovergrad*gradz;
}

void Grid::RebuildCellMarkerListPointers(void){
	Point pp;
	int i,j,k,dropped_markers = 0,num_newmarkers = num_markers;
	Cell* c;
	float half=0.5f,ix,iy,iz;

	for (i = 0; i < m-1; i++)
		for ( j = 0; j < n-1; j++)
			for ( k = 0; k < p-1; k++){
				c = &cell_list[i][j][k];
				c->num_particles = 0;
			}
	for (i = 0; i < num_markers; i++){
		pp = marker_list[i].p;
		ix = (pp.x-x_lo)/dx, iy = (pp.y-y_lo)/dy, iz = (pp.z-z_lo)/dz;
		//cout << ix << " " << iy << " " << iz << endl;
		if ((ix >= half)&&(ix < m-1-half)&&(iy >= half)&&(iy < n-1-half)&&(iz >= half)&&(iz < p-1-half)){
			c = &cell_list[int(ix)][int(iy)][int(iz)];
			if (c->num_particles < MAX_CELL_PARTICLES){
				c->markers[c->num_particles++] = i-dropped_markers;	
				marker_list[i-dropped_markers+1] = marker_list[i+1];
			}
			else{//delete marker
				num_newmarkers--;
				if (i<num_markers-1)
					marker_list[i-dropped_markers] = marker_list[i+1];
				dropped_markers++;
			}
		}
		else{//delete marker
				num_newmarkers--;
				if (i<num_markers-1)
					marker_list[i-dropped_markers] = marker_list[i+1];
				dropped_markers++;
			}
	}
	num_markers = num_newmarkers;
}

int Grid::CountFullCells(void){
	int i,j,k,fullcells=0;
	Cell* c;

	for (i = 0; i < m-1; i++)
		for ( j = 0; j < n-1; j++)
			for ( k = 0; k < p-1; k++){
				c = &cell_list[i][j][k];
				if (c->num_particles == MAX_CELL_PARTICLES)
					fullcells++;
			}
	return fullcells;
}

int Grid::CountNonemptyCells(void){
	int i,j,k,necells=0;
	Cell* c;

	for (i = 0; i < m-1; i++)
		for ( j = 0; j < n-1; j++)
			for ( k = 0; k < p-1; k++){
				c = &cell_list[i][j][k];
				if (c->num_particles > 0)
					necells++;
			}
	return necells;
}

bool Grid::CellChangesSignInside(Cell* c, float threshold){
	int i = c->i, j = c->j, k = c->k;
	int t2,t3,t4,t5,t6,t7,t8,p=SIGN(phi[i][j][k]-threshold);
		
	t2 = SIGN(phi[i][j+1][k]-threshold);
	t3 = SIGN(phi[i+1][j][k]-threshold);
	t4 = SIGN(phi[i+1][j+1][k]-threshold);
	t5 = SIGN(phi[i][j][k+1]-threshold);
	t6 = SIGN(phi[i][j+1][k+1]-threshold);
	t7 = SIGN(phi[i+1][j][k+1]-threshold);
	t8 = SIGN(phi[i+1][j+1][k+1]-threshold);

	if ((p*t2<=0)||(p*t3<=0)||(p*t4<=0)||(p*t5<=0)||(p*t6<=0)||(p*t7<=0)||(p*t8<=0))
		return true;
	else
		return false;
}

bool Grid::InsideSolid(Cell* c, float threshold){
	int i = c->i, j = c->j, k = c->k;
	int t2,t3,t4,t5,t6,t7,t8,p=SIGN(solidphi[i][j][k]-threshold);
		
	t2 = SIGN(solidphi[i][j+1][k]-threshold);
	t3 = SIGN(solidphi[i+1][j][k]-threshold);
	t4 = SIGN(solidphi[i+1][j+1][k]-threshold);
	t5 = SIGN(solidphi[i][j][k+1]-threshold);
	t6 = SIGN(solidphi[i][j+1][k+1]-threshold);
	t7 = SIGN(solidphi[i+1][j][k+1]-threshold);
	t8 = SIGN(solidphi[i+1][j+1][k+1]-threshold);

	if ((p*t2<=0)||(p*t3<=0)||(p*t4<=0)||(p*t5<=0)||(p*t6<=0)||(p*t7<=0)||(p*t8<=0))
		return false;
	else if (p>0)//inside of solid has negative sign 
		return false;
	else
		return true;
}

void Grid::WriteMeshFile(char* filename){
	int i, j, vindex=0, flag;
	float M = MAX(m,MAX(n,p));
	vect** vertex = new vect *[100000], *temp = new vect(0,0,0);
	for (i = 0; i < 100000; i++)
		vertex[i] = new vect(0,0,0);
	int** trlist = new int *[300000];
	for (i = 0; i < 300000; i++)
		trlist[i] = new int[3];
	int* trlistbool = new int [300000];
	for (i = 0; i < 300000; i++)
		trlistbool[i] = 1;

	vertex[0] = &tr[0].a;
	vertex[1] = &tr[0].b;
	vertex[2] = &tr[0].c;
	trlist[0][0] = 0;
	trlist[0][1] = 1;
	trlist[0][2] = 2;
	vindex = 4;
	for (i = 1; i < trsize; i++) {
		//work on a
		temp->fX = tr[i].a.fX;
		temp->fY = tr[i].a.fY;
		temp->fZ = tr[i].a.fZ;
		flag = 0;
		for ( j = 0; j < vindex; j++){
			if ((temp->fX == vertex[j]->fX)&&(temp->fY == vertex[j]->fY)&&(temp->fZ == vertex[j]->fZ)){
				trlist[i][0] = j;
				flag = 1;
				break;
			}
		}
		if (flag == 0){
			trlist[i][0] = vindex;
			vertex[vindex] = &tr[i].a;
			vindex++;
		}
		//work on b
		temp->fX = tr[i].b.fX;
		temp->fY = tr[i].b.fY;
		temp->fZ = tr[i].b.fZ;
		flag = 0;
		for ( j = 0; j < vindex; j++){
			if ((temp->fX == vertex[j]->fX)&&(temp->fY == vertex[j]->fY)&&(temp->fZ == vertex[j]->fZ)){
				trlist[i][1] = j;
				flag = 1;
				break;
			}
		}
		if (flag == 0){
			trlist[i][1] = vindex;
			vertex[vindex] = &tr[i].b;
			vindex++;
		}
		//work on c
		temp->fX = tr[i].c.fX;
		temp->fY = tr[i].c.fY;
		temp->fZ = tr[i].c.fZ;
		flag = 0;
		for ( j = 0; j < vindex; j++){
			if ((temp->fX == vertex[j]->fX)&&(temp->fY == vertex[j]->fY)&&(temp->fZ == vertex[j]->fZ)){
				trlist[i][2] = j;
				flag = 1;
				break;
			}
		}
		if (flag == 0){
			trlist[i][2] = vindex;
			vertex[vindex] = &tr[i].c;
			vindex++;
		}
		//now (a,b,c) are indexed as new or old vertices	
    }

	//fix the mesh - delete triangles with repeated vertices
	int bad_trngl = 0;
	for (i = 0; i < trsize; i++){if ((trlist[i][0] == trlist[i][1])||(trlist[i][1] == trlist[i][2])||(trlist[i][0] == trlist[i][2]))
			trlistbool[i] = 0;
	}
	FILE *fid = fopen (filename, "w+t");
	for (i = 0; i < vindex; i++)
		fprintf(fid, "%s %f %f %f \n","v", vertex[i]->fX, vertex[i]->fZ, - vertex[i]->fY );
	for (i = 0; i < trsize; i++){
		if (trlistbool[i])
			fprintf(fid, "%s %d %d %d \n","f", trlist[i][0]+1, trlist[i][1]+1,trlist[i][2]+1);
	}
	fclose ( fid );
}

void Grid::update_levelsetSL(float dt){
	int i,j,k;
	float ***temp_phi = new float **[m];
	for (i = 0; i < m; i++){
		temp_phi[i] = new float *[n];
		for (j = 0; j < n; j++){
			temp_phi[i][j] = new float [p];
			for (k = 0; k < p; k++)
				temp_phi[i][j][k] = phi[i][j][k];
		}
	}

	//do the advection
	advect(0,phi,temp_phi,u,v,w,dt);
	iter++;

	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			delete [] temp_phi[i][j];
    for (i = 0; i < m; i++)
		delete [] temp_phi[i];
	delete [] temp_phi;
}

void Grid::update_levelsetDL(float dt){//DL = Dupont & Liu
	int i,j,k;
	float ***temp_phi1 = new float **[m], extraphi = 0.0f;
	float ***temp_phi2 = new float **[m];
	for (i = 0; i < m; i++){
		temp_phi1[i] = new float *[n];
		temp_phi2[i] = new float *[n];
		for (j = 0; j < n; j++){
			temp_phi1[i][j] = new float [p];
			temp_phi2[i][j] = new float [p];
			for (k = 0; k < p; k++){
				temp_phi1[i][j][k] = phi[i][j][k];
				temp_phi2[i][j][k] = phi[i][j][k];
			}
		}
	}
	//do the advection
	advect(0,temp_phi2,temp_phi1,u,v,w,dt);
	advect(0,temp_phi1,temp_phi2,u,v,w,-dt);
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			for (k = 0; k < p; k++){
				extraphi = (phi[i][j][k]-temp_phi1[i][j][k])/2.0;
				if (fabs(extraphi) < 0.001)
					extraphi = 0;
				temp_phi1[i][j][k] = phi[i][j][k]+extraphi;
				temp_phi2[i][j][k] = extraphi;
		}
	advect(0,phi,temp_phi1,u,v,w,dt);
	iter++;

	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++){
			delete [] temp_phi1[i][j];
			delete [] temp_phi2[i][j];
		}
	for (i = 0; i < m; i++){
		delete [] temp_phi1[i];
		delete [] temp_phi2[i];
	}
	delete [] temp_phi1;delete [] temp_phi2;

}

float Grid::Interpolate(float*** field, int order, int i, int j, int k, float ldx, float ldy, float ldz){
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
				for (int kk = -st; kk <=1+st; kk++){
					if ((i+ii>=0)&&(i+ii<m)&&(j+jj>=0)&&(j+jj<n)&&(k+kk>0)&&(k+kk<p))
						interp += field[i+ii][j+jj][k+kk]*funct(ldx-ii)*funct(ldy-jj)*funct(ldz-kk);
				}
		return interp;
}

void Grid::Redistance1(void){
	int i,j,k, h, l, lk, kk;
	Cell *d;
	Point pp;
	float maxdist = 3.0f*dx, ldx, ldy, ldz, ldxi, ldyj, ldzk, q, c=0.2f, rho=1.8f;
		
	for (i=0; i<m; i++)
		for (j=0; j< n; j++)
			for (k=0; k<p; k++)
				temp_phi[i][j][k] = 0;
	for (i=1; i<m-1; i++)
		for (j=1; j< n-1; j++)
			for (k=1; k<p-1; k++){
				float lambda = 0, wght1 = 0, twght = 0;  
				int lim = 1;
				//use the cells about (i,j) to find correction for phi
				for (h = -lim; h <=lim-1; h++){
					if ((0<i+h)&&(i+h<m-2)){
						for (l = -lim; l <=lim-1; l++){
							if ((0<j+l)&&(j+l<n-2)){
								for (lk = -lim; lk <=lim-1; lk++){
									if ((0<k+lk)&&(k+lk<p-2)){
										d = &cell_list[i+h][j+l][k+lk];
										for (kk=0; kk < d->num_particles; kk++){
											 pp = marker_list[d->markers[kk]].p;
											 ldx = pp.x-x_lo-(i+h)*dx; ldx /= dx;
											 ldy = pp.y-y_lo-(j+l)*dy; ldy /= dy;
											 ldz = pp.z-z_lo-(k+lk)*dz; ldz /= dz;
											 ldxi = pp.x-x_lo-i*dx; ldxi /= dx;
											 ldyj = pp.y-y_lo-j*dy; ldyj /= dy;
											 ldzk = pp.z-z_lo-k*dz; ldzk /= dz;
											 q = sqrt(ldxi*ldxi+ldyj*ldyj+ldzk*ldzk)/rho;
											 if (q<1) wght1= (exp(-q*q/(c*c))-exp(-1/(c*c)))/(1-exp(-1/(c*c)));
											 else wght1 = 0;
											 twght += wght1;
											 lambda += wght1*Interpolate(phi,3,i+h,j+l,k+lk,ldx,ldy,ldz);//1st or 3rd order interp
										}
									}
								}
							}
						}
					}
				}
				if (twght > 0){
					phi_recon[i][j][k] = phi[i][j][k]-lambda/twght;
					temp_phi[i][j][k] = 1;
				}
			}

	//reset phi_recon
	for (i=0; i < m; i++)
		for (j=0; j<n; j++)
			for (k=0; k<p; k++){
				if (temp_phi[i][j][k])
					phi[i][j][k] = phi_recon[i][j][k];
				phi_recon[i][j][k] = maxdist;
			}
}

bool Grid::ChangeDensity(float level){
	int ix, iy, iz, i, j, k;
	Point pp;
	float ldx, ldy, ldz, q, c=0.2f, rho=1.74f, oneoverdx=1.0f/dx, wght, den_liquid=1.0;
	float expc = exp(-1/(c*c));
	float ooexpc= 1/(1-expc);
	for (int ifm = 0; ifm < num_free_markers; ifm++){
		pp = free_markers[ifm].p;
		ix = int((pp.x-x_lo)*oneoverdx), iy = int((pp.y-y_lo)*oneoverdx), iz = int((pp.z-z_lo)*oneoverdx);
		if ((ix < 0 )||(ix >= m-1)||(iy < 0 )||(iy >= n-1)||(iz < 0 )||(iz >= p-1))
		//if ((ix < 1 )||(ix >= m-2)||(iy < 1 )||(iy >= n-2)||(iz < 1 )||(iz >= p-2))
			return false;	
		for (i = 0; i <= 1; i++)
			for (j = 0; j <= 1; j++)
				for (k = 0; k <= 1; k++){
					ldx = (pp.x-x_lo)*oneoverdx-(ix+i);
					ldy = (pp.y-y_lo)*oneoverdx-(iy+j);
					ldz = (pp.z-z_lo)*oneoverdx-(iz+k);

					q = sqrt(ldx*ldx+ldy*ldy+ldz*ldz)/rho;
					if (q<1) wght= (exp(-q*q/(c*c))-expc)*ooexpc;
					else wght = 0;
										
					den[ix+i][iy+j][iz+k] += wght;
					
					if (den[ix+i][iy+j][iz+k]>=den_liquid)
						den[ix+i][iy+j][iz+k]=den_liquid;
				}
	}
	//check density sign
	for (i=0; i<m; i++)
		for (j=0; j< n; j++)
			for (k=0; k<p; k++){
				if (den[i][j][k] <= 0.00125)//was <=0
					den[i][j][k]=0.00125;
			}
	UpdateGasLevel(level);
	UpdateGeneralLevel(level);
	return true;
}

void Grid::UpdateGasLevel(float level){
	int i,j,k;
	Cell* c;
	float center_zval = 0;

	cout << "level = " << level << '\n';
	for (i=0; i<m; i++)
		for (j=0; j< n; j++)
			for (k=0; k<p; k++){
				if (solidphi[i][j][k]<0){//node is inside solid
					center_zval = z_lo+dz*k;
					if (center_zval > level)
						phi[i][j][k] = -fabs(phi[i][j][k]);
					else
						phi[i][j][k] = fabs(phi[i][j][k]);
				}
			}
}

void Grid::UpdateGeneralLevel(float level){
	int i,j,k;
	Cell* c;
	float center_zval = 0;

	for (i=0; i<m; i++)
		for (j=0; j< n; j++)
			for (k=0; k<p; k++){
				center_zval = z_lo+dz*k;
				if ((phi[i][j][k]<0)&&(center_zval>level+2*dx)&&(den[i][j][k]>0.95))//node is inside air
						phi[i][j][k] = 0.95-den[i][j][k];
			}
}

void Grid::Reinit(int itern, float dtau){
	int i,j,k,ii;
	float fact, grad, den, den1, eps = 1e-1f, epsilon = dx/6.0, d;

	for (i=0; i < m; i++)
		for (j=0; j<n; j++)
			for (k=0; k<p; k++)
				phi_recon[i][j][k] = phi[i][j][k];
	//work inside the band	
	/*for (i = 1; i < m-1; i++)
		for (j = 1; j < n-1; j++)
			for (k = 1; k < p-1; k++){
				if (phi[i][j][k] > 0)
					den = float(sqrt(MAX(pow(MAX((phi[i][j][k]-phi[i-1][j][k])/dx,0),2),pow(MIN((phi[i+1][j][k]-phi[i][j][k])/dx,0),2))+
						MAX(pow(MAX((phi[i][j][k]-phi[i][j-1][k])/dy,0),2),pow(MIN((phi[i][j+1][k]-phi[i][j][k])/dy,0),2))+
						MAX(pow(MAX((phi[i][j][k]-phi[i][j][k-1])/dz,0),2),pow(MIN((phi[i][j][k+1]-phi[i][j][k])/dz,0),2))));
				else 
					den = float(sqrt(MAX(pow(MAX((phi[i+1][j][k]-phi[i][j][k])/dx,0),2),pow(MIN((phi[i][j][k]-phi[i-1][j][k])/dx,0),2))+
						MAX(pow(MAX((phi[i][j+1][k]-phi[i][j][k])/dy,0),2),pow(MIN((phi[i][j][k]-phi[i][j-1][k])/dy,0),2))+
						MAX(pow(MAX((phi[i][j][k+1]-phi[i][j][k])/dz,0),2),pow(MIN((phi[i][j][k]-phi[i][j][k-1])/dz,0),2))));
				if ((phi[i][j][k]*phi[i-1][j][k] <= 0)||(phi[i][j][k]*phi[i+1][j][k] <= 0)||
					(phi[i][j][k]*phi[i][j-1][k] <= 0)||(phi[i][j][k]*phi[i][j+1][k] <= 0)||
					(phi[i][j][k]*phi[i][j][k-1] <= 0)||(phi[i][j][k]*phi[i][j][k+1] <= 0)){								
					phi_recon[i][j][k] /= den;//or /1
				}	
				else{
						//outside the narrow band
						if (phi[i][j][k] > 0)
							grad = float(sqrt(MAX(pow(MAX((phi[i][j][k]-phi[i-1][j][k])/dx,0),2),pow(MIN((phi[i+1][j][k]-phi[i][j][k])/dx,0),2))+
							MAX(pow(MAX((phi[i][j][k]-phi[i][j-1][k])/dy,0),2),pow(MIN((phi[i][j+1][k]-phi[i][j][k])/dy,0),2))+
							MAX(pow(MAX((phi[i][j][k]-phi[i][j][k-1])/dz,0),2),pow(MIN((phi[i][j][k+1]-phi[i][j][k])/dz,0),2))))-1;
						else 
							grad = float(sqrt(MAX(pow(MAX((phi[i+1][j][k]-phi[i][j][k])/dx,0),2),pow(MIN((phi[i][j][k]-phi[i-1][j][k])/dx,0),2))+
							MAX(pow(MAX((phi[i][j+1][k]-phi[i][j][k])/dy,0),2),pow(MIN((phi[i][j][k]-phi[i][j-1][k])/dy,0),2))+
							MAX(pow(MAX((phi[i][j][k+1]-phi[i][j][k])/dz,0),2),pow(MIN((phi[i][j][k]-phi[i][j][k-1])/dz,0),2))))-1;
						fact = dtau*grad*phi[i][j][k]/sqrt(pow(phi[i][j][k],2)+pow(epsilon,2));
						phi_recon[i][j][k] -= fact;
				}
			}*/
	//work only outside the band
	for (ii = 0; ii < itern; ii++){
		for (i=1 ; i<m-1 ; i++){ 
			for (j = 1; j < n-1; j++)
				for (k = 1; k < p-1; k++){
					if ((phi[i][j][k]*phi[i-1][j][k] <= 0)||(phi[i][j][k]*phi[i+1][j][k] <= 0)||
						(phi[i][j][k]*phi[i][j-1][k] <= 0)||(phi[i][j][k]*phi[i][j+1][k] <= 0)||
						(phi[i][j][k]*phi[i][j][k-1] <= 0)||(phi[i][j][k]*phi[i][j][k+1] <= 0)){
						fact = 0;
					}
					else{
						//outside the narrow band
						if (phi[i][j][k] > 0)
							grad = float(sqrt(MAX(pow(MAX((phi_recon[i][j][k]-phi_recon[i-1][j][k])/dx,0),2),pow(MIN((phi_recon[i+1][j][k]-phi_recon[i][j][k])/dx,0),2))+
							MAX(pow(MAX((phi_recon[i][j][k]-phi_recon[i][j-1][k])/dy,0),2),pow(MIN((phi_recon[i][j+1][k]-phi_recon[i][j][k])/dy,0),2))+
							MAX(pow(MAX((phi_recon[i][j][k]-phi_recon[i][j][k-1])/dz,0),2),pow(MIN((phi_recon[i][j][k+1]-phi_recon[i][j][k])/dz,0),2))))-1;
						else 
							grad = float(sqrt(MAX(pow(MAX((phi_recon[i+1][j][k]-phi_recon[i][j][k])/dx,0),2),pow(MIN((phi_recon[i][j][k]-phi_recon[i-1][j][k])/dx,0),2))+
							MAX(pow(MAX((phi_recon[i][j+1][k]-phi_recon[i][j][k])/dy,0),2),pow(MIN((phi_recon[i][j][k]-phi_recon[i][j-1][k])/dy,0),2))+
							MAX(pow(MAX((phi_recon[i][j][k+1]-phi_recon[i][j][k])/dz,0),2),pow(MIN((phi_recon[i][j][k]-phi_recon[i][j][k-1])/dz,0),2))))-1;
						fact = dtau*grad*phi[i][j][k]/sqrt(pow(phi[i][j][k],2)+pow(epsilon,2));
					}
					
					phi_recon[i][j][k] -= fact;
				}
		}
	}
	for (i = 1; i < m-1; i++)
		for (j = 1; j < n-1; j++)
			for (k = 1; k < p-1; k++){//update only outside the narrow band
				if ((phi[i][j][k]*phi[i-1][j][k] <= 0)||(phi[i][j][k]*phi[i+1][j][k] <= 0)||
						(phi[i][j][k]*phi[i][j-1][k] <= 0)||(phi[i][j][k]*phi[i][j+1][k] <= 0)||
						(phi[i][j][k]*phi[i][j][k-1] <= 0)||(phi[i][j][k]*phi[i][j][k+1] <= 0)){}
				else{
					phi[i][j][k] = phi_recon[i][j][k];
				}
			}	
}


void Grid::DeleteExtraMarkers(void){
	int i,j,k,ix,iy,iz,side = 1, dropped_markers=0,new_nummarkers = num_markers; 
	Point pp;
	Cell* c;
	float rand_radius, web = 100.0f;//4.0 for ballnew
	/**************
	//num_free_markers = 0;//(UN)COMMENT THIS IF NEEDED, DEPENDING ON POST-PROCESSING OPTION
	**************/
	for (i = 0; i < num_markers; i++){
		pp = marker_list[i].p;	
		ix = int((pp.x-x_lo)/dx), iy = int((pp.y-y_lo)/dx), iz = int((pp.z-z_lo)/dx);
		c = &cell_list[ix][iy][iz];
		if ((!CellChangesSignInside(c, 0))&&(c->num_particles>0)&&(phi[ix][iy][iz]!=0)){
			//update free marker list
				
			if ((ix > side)&&(ix < m-2-side)&&(iy > side)&&(iy < n-2-side)&&(iz > side)&&(iz < p-2-side)){//avoid boundaries
				//if ((num_free_markers<m*n*4*40)&&(!InsideSolid(c,0))){//we don't delete markers inside the solid
				float vel_var = VelocityVariation(c), curv = AbsCurvature(c);			
				if ((num_free_markers<m*n*4*40)&&(!InsideSolid(c,0))&&(vel_var>1000)){
				//	std::cout << "vel_var " << vel_var << " curv " << curv << '\n';
				//if ((num_free_markers<m*n*4*40)&&(!InsideSolid(c,0))&&(WeberThreshold(c,web))){//we also delete markers inside the solid
					if (phi[ix][iy][iz]>0)
						marker_list[i].transp = -1;
					else
						marker_list[i].transp = 1;//-1 means gas, 1 means liquid - assuming phi(liquid)>0
					rand_radius = float(rand())/RAND_MAX+0.5;
					marker_list[i].transp *= rand_radius;
					free_markers[num_free_markers] = marker_list[i];
					free_markers[num_free_markers].color = InterpVelocity(pp);//we use the color vector for velocity
					num_free_markers++;
				}
			}
			//inside particle, delete it
			c->num_particles--; 
			new_nummarkers--;
			if (i<num_markers-1)
				marker_list[i-dropped_markers] = marker_list[i+1];
			dropped_markers++;
		
		}
		else if (i<num_markers-1)
			marker_list[i-dropped_markers+1] = marker_list[i+1];
	}
	cout << "DROPPED MARKERS " << dropped_markers << endl;
	num_markers = new_nummarkers;

}

bool Grid::WeberThreshold(Cell* c, float web){
	float vel_var, curv;
	vel_var = VelocityVariation(c);
	curv = AbsCurvature(c);
	std::cout << "vel_var " << vel_var << " curv " << curv << '\n';
	if (NotBoundaryCell(c))
		if (vel_var*curv > web)
			return true;
		else{
			std::cout << " rejected " << '\n';
			return false;
		}
	else
		return false;
}

bool Grid::NotBoundaryCell(Cell* c){
	int ix = c->i, iy = c->j, iz = c->k, side = 0;
	if ((ix > side)&&(ix < m-2-side)&&(iy > side)&&(iy < n-2-side)&&(iz > side)&&(iz < p-2-side))
		return true;
	else
		return false;
}

float Grid::AbsCurvature(Cell* c){
	//this routine assumes c is not a boundary cell
	int i = c->i, j = c->j, k = c->k, ii, jj;
	float oodx = 1/dx;
	float phix = ((phi[i+1][j][k]-phi[i][j][k])+(phi[i+1][j][k+1]-phi[i][j][k+1])+
				(phi[i+1][j+1][k]-phi[i][j+1][k])+(phi[i+1][j+1][k+1]-phi[i][j+1][k+1]))*0.25f*oodx;
	float phiy = ((phi[i][j+1][k]-phi[i][j][k])+(phi[i][j+1][k+1]-phi[i][j][k+1])+
				(phi[i+1][j+1][k]-phi[i+1][j][k])+(phi[i+1][j+1][k+1]-phi[i+1][j][k+1]))*0.25f*oodx;
	float phiz = ((phi[i][j][k+1]-phi[i][j][k])+(phi[i+1][j][k+1]-phi[i+1][j][k])+
				(phi[i][j+1][k+1]-phi[i][j+1][k])+(phi[i+1][j+1][k+1]-phi[i+1][j+1][k]))*0.25f*oodx;
	float phixx = ((phi[i+2][j][k]+phi[i-1][j][k]-phi[i][j][k]-phi[i+1][j][k])+
				(phi[i+2][j][k+1]+phi[i-1][j][k+1]-phi[i][j][k+1]-phi[i+1][j][k+1])+
				(phi[i+2][j+1][k]+phi[i-1][j+1][k]-phi[i][j+1][k]-phi[i+1][j+1][k])+
				(phi[i+2][j+1][k+1]+phi[i-1][j+1][k+1]-phi[i][j+1][k+1]-phi[i+1][j+1][k+1]))*0.125f*oodx*oodx;
	float phiyy = ((phi[i][j+2][k]+phi[i][j-1][k]-phi[i][j][k]-phi[i][j+1][k])+
				(phi[i][j+2][k+1]+phi[i][j-1][k+1]-phi[i][j][k+1]-phi[i][j+1][k+1])+
				(phi[i+1][j+2][k]+phi[i+1][j-1][k]-phi[i+1][j][k]-phi[i+1][j+1][k])+
				(phi[i+1][j+2][k+1]+phi[i+1][j-1][k+1]-phi[i+1][j][k+1]-phi[i+1][j+1][k+1]))*0.125f*oodx*oodx;
	float phizz = ((phi[i][j][k+2]+phi[i][j][k-1]-phi[i][j][k]-phi[i][j][k+1])+
				(phi[i+1][j][k+2]+phi[i+1][j][k-1]-phi[i+1][j][k]-phi[i+1][j][k+1])+
				(phi[i][j+1][k+2]+phi[i][j+1][k-1]-phi[i][j+1][k]-phi[i][j+1][k+1])+
				(phi[i+1][j+1][k+2]+phi[i+1][j+1][k-1]-phi[i+1][j+1][k]-phi[i+1][j+1][k+1]))*0.125f*oodx*oodx;
	float phixy = ((phi[i+1][j+1][k]+phi[i][j][k]-phi[i+1][j][k]-phi[i][j+1][k])+
				(phi[i+1][j+1][k+1]+phi[i][j][k+1]-phi[i+1][j][k+1]-phi[i][j+1][k+1]))*0.25f*oodx*oodx;
	float phixz = ((phi[i+1][j][k+1]+phi[i][j][k]-phi[i+1][j][k]-phi[i][j][k+1])+
				(phi[i+1][j+1][k+1]+phi[i][j+1][k]-phi[i+1][j+1][k]-phi[i][j+1][k+1]))*0.25f*oodx*oodx;
	float phiyz = ((phi[i][j+1][k+1]+phi[i][j][k]-phi[i][j+1][k]-phi[i][j][k+1])+
				(phi[i+1][j+1][k+1]+phi[i+1][j][k]-phi[i+1][j+1][k]-phi[i+1][j][k+1]))*0.25f*oodx*oodx;
	float gradphi = sqrt(phix*phix+phiy*phiy+phiz*phiz);
	float curv = (phix*phix*phiyy-2*phix*phiy*phixy+phiy*phiy*phixx+phix*phix*phizz-2*phix*phiz*phixz+
				phiz*phiz*phixx+phiy*phiy*phizz-2*phiy*phiz*phiyz+phiz*phiz*phiyy)/(gradphi*gradphi*gradphi);
	float abscurv = (float)fabs(curv);
	return abscurv;
}

float Grid::VelocityVariation(Cell* c){
	int i = c->i, j = c->j, k = c->k, ii, jj;
	vect v[8]={vect(F1[i][j][k],F2[i][j][k],F3[i][j][k]),vect(F1[i+1][j][k],F2[i+1][j][k],F3[i+1][j][k]), 
		vect(F1[i][j+1][k],F2[i][j+1][k],F3[i][j+1][k]),vect(F1[i+1][j+1][k],F2[i+1][j+1][k],F3[i+1][j+1][k]),
		vect(F1[i][j][k+1],F2[i][j][k+1],F3[i][j][k+1]),vect(F1[i+1][j][k+1],F2[i+1][j][k+1],F3[i+1][j][k+1]),
		vect(F1[i][j+1][k+1],F2[i][j+1][k+1],F3[i][j+1][k+1]),vect(F1[i+1][j+1][k+1],F2[i+1][j+1][k+1],F3[i+1][j+1][k+1])};
	float len, min, max;
	/*min = 1e10;
	for (ii= 0; ii<8; ii++)
		for (jj=0; jj < 8; jj++){
			len = (v[ii].fX-v[jj].fX)*(v[ii].fX-v[jj].fX)+(v[ii].fY-v[jj].fY)*(v[ii].fY-v[jj].fY)+
				(v[ii].fZ-v[jj].fZ)*(v[ii].fZ-v[jj].fZ);	
			if ((len > 0) && (min > len))
					min = len;
		}
	min = sqrt(min);
	return min;
	*/	
	/*max = 0;
	for (ii= 0; ii<8; ii++)
		for (jj=0; jj < 8; jj++){
			len = (v[ii].fX-v[jj].fX)*(v[ii].fX-v[jj].fX)+(v[ii].fY-v[jj].fY)*(v[ii].fY-v[jj].fY)+
				(v[ii].fZ-v[jj].fZ)*(v[ii].fZ-v[jj].fZ);	
			if (max < len)
					max = len;
		}
	return max;*/
	/**/
	float oodx = 1/dx;
	float phix = ((phi[i+1][j][k]-phi[i][j][k])+(phi[i+1][j][k+1]-phi[i][j][k+1])+
				(phi[i+1][j+1][k]-phi[i][j+1][k])+(phi[i+1][j+1][k+1]-phi[i][j+1][k+1]))*0.25f*oodx;
	float phiy = ((phi[i][j+1][k]-phi[i][j][k])+(phi[i][j+1][k+1]-phi[i][j][k+1])+
				(phi[i+1][j+1][k]-phi[i+1][j][k])+(phi[i+1][j+1][k+1]-phi[i+1][j][k+1]))*0.25f*oodx;
	float phiz = ((phi[i][j][k+1]-phi[i][j][k])+(phi[i+1][j][k+1]-phi[i+1][j][k])+
				(phi[i][j+1][k+1]-phi[i][j+1][k])+(phi[i+1][j+1][k+1]-phi[i+1][j+1][k]))*0.25f*oodx;
	float gradphi = sqrt(phix*phix+phiy*phiy+phiz*phiz);
	if (gradphi == 0)
		gradphi = 1;
	max = 0;
	vect diff, normal = vect(phix/gradphi,phiy/gradphi,phiz/gradphi), proj;
	float coeff;
	for (ii= 0; ii<8; ii++)
		for (jj=0; jj < 8; jj++){
			diff = vect(v[ii].fX-v[jj].fX, v[ii].fY-v[jj].fY, v[ii].fZ-v[jj].fZ);
			coeff = (normal.fX*diff.fX+normal.fY*diff.fY+normal.fZ*diff.fZ)/gradphi/gradphi;//v*n/n^2
			proj = vect(diff.fX-coeff*normal.fX,diff.fY-coeff*normal.fY,diff.fZ-coeff*normal.fZ);
			len = proj.fX*proj.fX+proj.fY*proj.fY+proj.fZ*proj.fZ;
			if (max < len){
					max = len;
			}
		}
	return max;	
}

bool Grid::AddParticles(int boiling){
	//this is done after solving the reconnection problem and after we get the new cell segments from marching squares
	//it simply uses the reconstructed segments to insert new particles along them
	
	int i,j,k,l,index, initial_nummarkers=num_markers;
	Cell* c;
	bool done = 0;
	float dist,ldx,ldy,ldz,lim1 = 1.75*dx;

	for (i = 0; i < m-1; i++)
		for (j = 0; j < n-1; j++)
			for (k = 0; k < p-1; k++){
			c = &cell_list[i][j][k];
			index = 0;
			if ((fabs(phi[i][j][k]) <= lim1)&&(fabs(phi[i+1][j][k]) <= lim1)&&
				(fabs(phi[i][j+1][k]) <= lim1)&&(fabs(phi[i+1][j+1][k]) <= lim1)&&
				(fabs(phi[i][j+1][k]) <= lim1)&&(fabs(phi[i+1][j+1][k]) <= lim1)){
				if (CellChangesSignInside(c, 0)&&(!InsideSolid(c,0))){//we don't add particles inside solid
					if (c->num_particles <= 1){
						done = 1;
						AddParticlesToCell(c,initial_nummarkers,boiling);
					}
				}
			}
		}
	return done;
}

void Grid::AddParticlesToCell(Cell* cl, int num_mark, int boiling){
	//version1, based on the level set linear reconstruction
	Point pp;
	vect a,b,c;
	double area,la,lb,lc,sp;
	float t,at,bt,ct;
	int kk,l,np,li,lj;

	for (kk = 0; kk < cl->num_triangles; kk++){	
		if (!InsideSolid(cl,0)){
			a = tr[cl->triangles[kk]].a;
			b = tr[cl->triangles[kk]].b;
			c = tr[cl->triangles[kk]].c;
			la = sqrt(pow2float(a.fX-b.fX)+pow2float(a.fY-b.fY)+pow2float(a.fZ-b.fZ));
			lb = sqrt(pow2float(a.fX-c.fX)+pow2float(a.fY-c.fY)+pow2float(a.fZ-c.fZ));
			lc = sqrt(pow2float(c.fX-b.fX)+pow2float(c.fY-b.fY)+pow2float(c.fZ-b.fZ));
			//compute the semiperimeter
			sp = (la+lb+lc)/2;
			//and area
			area = sqrt(sp*(sp-la)*(sp-lb)*(sp-lc));
			//estimate number of particles (on any triangle edge) 
			np = int(floor(area/0.1f));
			if (np < 2)
				np = 2;
			//we place np*(np+1)/2 particles
			for (li = 0; li <= np; li++)
				for (lj = 0; lj <= li; lj++){
					t = li/float(np);
					at = x_lo+dx*(a.fX+li*(b.fX-a.fX)/np+lj*(c.fX-b.fX)/(li+1));
					bt = y_lo+dy*(a.fY+li*(b.fY-a.fY)/np+lj*(c.fY-b.fY)/(li+1));
					ct = z_lo+dz*(a.fZ+li*(b.fZ-a.fZ)/np+lj*(c.fZ-b.fZ)/(li+1));
					if (at >= x_hi)
						at = x_hi-10e-6;
					if (bt >= y_hi)
						bt = y_hi-10e-6;
					if (ct >= z_hi)
						ct = z_hi-10e-6;
					pp = Point(at,bt,ct);
					if ((cl->num_particles < MAX_CELL_PARTICLES-8)&&(InsideDomain(&pp))){
						marker_list[num_markers].p = pp;
						cl->markers[cl->num_particles++] = num_markers;	
						if (boiling == 0)
							marker_list[num_markers].color = EstimateColor(pp,cl,num_mark);
						else
							marker_list[num_markers].color = BoilingColor(pp,cl,num_mark);
						num_markers++;
					}
				}
		}
	}
}

vect Grid::BoilingColor(Point pp, Cell *cl,int num_mark){
	float Rcol=1.0,Gcol=0.0,Bcol=0.0, xp = pp.x, yp = pp.y, zp = pp.z;
	if ((xp < 0)&&(yp < 0)){
		Rcol=1.0,Gcol=1.0,Bcol=1.0;}
	else if ((xp < 0)&&(yp > 0)){
		Rcol=1.0,Gcol=0.5,Bcol=0.5;}
	else if ((xp > 0)&&(yp < 0)){
		Rcol=0.5,Gcol=1.0,Bcol=0.5;}
	else{
		Rcol=0.5,Gcol=0.5,Bcol=1.0;}

	return vect(Rcol,Gcol,Bcol);
}

vect Grid::EstimateColor(Point pp, Cell *cl,int num_mark){
	//num_mark is the initial number of markers before RedistributeParticles is called
	//it is used to decide if a cetrain marker is newly added or not
	//in order to know if its color will be used or not
	float Rcol=0,Gcol=0,Bcol=0,rho=10.0,cc=0.08,q;
	float lambdaR=0.0,lambdaG=0.0,lambdaB=0.0,twght=0.0,wght1=0.0,ldxi,ldyj,ldzk;
	int i,j,k,kk,ix = cl->i, iy = cl->j, iz = cl->k, lim = 1,nm0; 
	Cell* d;

	for (i = ix-lim; i < ix+lim+1; i++)
		for (j = iy-lim; j < iy+lim+1; j++)
			for (k = iz-lim; k < iz+lim+1; k++){
				if ((i >= 0)&&(i < m-1)&&(j >= 0)&&(j < n-1)&&(k >= 0)&&(k < p-1)){
					d = &cell_list[i][j][k];
					nm0 = d->num_particles;
					if (nm0 > 0){
						for (kk = 0; kk < nm0; kk++){
							if (d->markers[kk] < num_mark){
								ldxi = (marker_list[d->markers[kk]].p.x-pp.x)/dx;
								ldyj = (marker_list[d->markers[kk]].p.y-pp.y)/dy;
								ldzk = (marker_list[d->markers[kk]].p.z-pp.z)/dz;
								q = sqrt(ldxi*ldxi+ldyj*ldyj+ldzk*ldzk)/rho;
								if (q<1)
									wght1= (exp(-q*q/(cc*cc))-exp(-1/(cc*cc)))/(1-exp(-1/(cc*cc)));
								else
									wght1 = 0;
								twght += wght1;
								lambdaR += wght1*marker_list[d->markers[kk]].color.fX;	
								lambdaG += wght1*marker_list[d->markers[kk]].color.fY;	
								lambdaB += wght1*marker_list[d->markers[kk]].color.fZ;	
							}
						}
					}
				}
			}

	if (twght > 0){
		Rcol = lambdaR/twght;
		Gcol = lambdaG/twght;
		Bcol = lambdaB/twght;
	}
	else{
		Rcol = 0.85; Gcol = 0.95; Bcol = 1.0;		
	}
	return vect(Rcol,Gcol,Bcol);
}

void Grid::set_bnd (int b, float *** x ){
	int i,j;
	//faces
	for ( i=1 ; i<n-1 ; i++ ) 
		for (j=1; j < p-1; j++){
			x[0][i][j] = b==1 ? -x[1][i][j] : x[1][i][j];
			x[m-1][i][j] = b==1 ? -x[m-2][i][j] : x[m-2][i][j];
		}
	for ( i = 1; i < m-1; i++)
		for (j = 1; j < p-1; j++){ 
		x[i][0][j] = b==2 ? -x[i][1][j] : x[i][1][j];
		x[i][n-1][j] = b==2 ? -x[i][n-2][j] : x[i][n-2][j];
		}
	for ( i = 1; i < m-1; i++)
		for (j = 1; j < n-1; j++){ 
		x[i][j][0] = b==3 ? -x[i][j][1] : x[i][j][1];
		x[i][j][p-1] = b==3 ? -x[i][j][p-2] : x[i][j][p-2];
		}
	//edges
	for ( i = 1 ; i<n-1 ; i++ ){ 
		x[0][i][0] = 0.5f*(x[0][i][1]+x[1][i][0]);
		x[m-1][i][0] = 0.5f*(x[m-1][i][1]+x[m-2][i][0]);
		x[0][i][p-1] = 0.5f*(x[0][i][p-1]+x[1][i][p-2]);
		x[m-1][i][p-1] = 0.5f*(x[m-1][i][p-2]+x[m-2][i][p-1]);
		}
	for ( i = 1; i < m-1; i++){
		x[i][0][0] = 0.5f*(x[i][1][0]+x[i][0][1]);
		x[i][n-1][0] = 0.5f*(x[i][n-2][0]+x[i][n-1][1]);
		x[i][0][p-1] = 0.5f*(x[i][1][p-1]+x[i][0][p-2]);
		x[i][n-1][p-1] = 0.5f*(x[i][n-1][p-2]+x[i][n-2][p-1]);
		}
	for (j = 1; j < p-1; j++){ 
		x[0][0][j] = 0.5f*(x[1][0][j]+x[0][1][j]);
		x[0][n-1][j] = 0.5f*(x[1][n-1][j]+x[0][n-2][j]);
		x[m-1][0][j] = 0.5f*(x[m-2][0][j]+x[m-1][1][j]);
		x[m-1][n-1][j] = 0.5f*(x[m-2][n-1][j]+x[m-1][n-2][j]);
		}
	//corners
	x[0][0][0] = (x[1][0][0]+x[0][1][0]+x[0][0][1])/3;
	x[0][n-1][0] = (x[1][n-1][0]+x[0][n-2][0]+x[0][n-1][1])/3;
	x[m-1][0][0] = (x[m-2][0][0]+x[m-1][1][0]+x[m-1][0][1])/3;
	x[m-1][n-1][0] = (x[m-2][n-1][0]+x[m-1][n-2][0]+x[m-1][n-1][1])/3;
	x[0][0][p-1] = (x[1][0][p-1]+x[0][1][p-1]+x[0][0][p-2])/3;
	x[0][n-1][p-1] = (x[1][n-1][p-1]+x[0][n-2][p-1]+x[0][n-1][p-2])/3;
	x[m-1][0][p-1] = (x[m-2][0][p-1]+x[m-1][1][p-1]+x[m-1][0][p-2])/3;
	x[m-1][n-1][p-1] = (x[m-2][n-1][p-1]+x[m-1][n-2][p-1]+x[m-1][n-1][p-2])/3;
}

void Grid::advect ( int b, float *** d, float *** d0, float *** u, float *** v, float *** w, float dt ){
	int i, j, k, i0, j0, k0, i1, j1, k1;
	float x, y, z, s0, t0, r0, s1, t1, r1, thr = 0.000001f;

	if ((bcxlo != 0)&&(bcylo != 0)){
		for (i = 1; i < m-1; i++)
			for (j = 1; j < n-1; j++)
				for (k = 1; k < p-1; k++){
					x = i-dt*u[i][j][k]/dx; y = j-dt*v[i][j][k]/dy; z = k-dt*w[i][j][k]/dz;
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
					if (z<thr){ 
						if (bczlo == 0)
							z = z+(p-1);
						else
							z=thr;
					}
					if (z>=p-1-thr){
						if (bczhi == 0)
							z = z-(p-1);
						else
							z=p-1-thr;
					}
					i0=(int)x; i1=i0+1;
					j0=(int)y; j1=j0+1;
					k0=(int)z; k1=k0+1;
					s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1; r1 = z-k0; r0 = 1-r1;
					d[i][j][k] = s0*(t0*(r0*d0[i0][j0][k0]+r1*d0[i0][j0][k1])
							 +t1*(r0*d0[i0][j1][k0]+r1*d0[i0][j1][k1]))+
						  s1*(t0*(r0*d0[i1][j0][k0]+r1*d0[i1][j0][k1])
							 +t1*(r0*d0[i1][j1][k0]+r1*d0[i1][j1][k1]));
				}
		set_bnd ( b, d );
	}
	else if ((bcxlo == 0)&&(bcylo != 0)){
		for (i = 0; i < m; i++)
			for (j = 1; j < n-1; j++)
				for (k = 1; k < p-1; k++){
					x = i-dt*u[i][j][k]/dx; y = j-dt*v[i][j][k]/dy; z = k-dt*w[i][j][k]/dz;
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
					if (z<thr){ 
						if (bczlo == 0)
							z = z+(p-1);
						else
							z=thr;
					}
					if (z>=p-1-thr){
						if (bczhi == 0)
							z = z-(p-1);
						else
							z=p-1-thr;
					}
					i0=(int)x; i1=i0+1;
					j0=(int)y; j1=j0+1;
					k0=(int)z; k1=k0+1;
					s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1; r1 = z-k0; r0 = 1-r1;
					d[i][j][k] = s0*(t0*(r0*d0[i0][j0][k0]+r1*d0[i0][j0][k1])
							 +t1*(r0*d0[i0][j1][k0]+r1*d0[i0][j1][k1]))+
						  s1*(t0*(r0*d0[i1][j0][k0]+r1*d0[i1][j0][k1])
							 +t1*(r0*d0[i1][j1][k0]+r1*d0[i1][j1][k1]));
				}
		//z boundaries
		for ( i = 0; i < m; i++)
			for (j = 0; j < n; j++){
				d[i][j][0] = d0[i][j][1];
				d[i][j][p-1] = d0[i][j][p-2];
			}
		//y boundaries
		for (i = 0; i < m; i++)
			for ( k = 0; k < p; k++){
				d[i][0][k] = d0[i][1][k];
				d[i][n-1][k] = d0[i][n-2][k];
			}
	}
	else if ((bcxlo != 0)&&(bcylo == 0)){
		for (i = 1; i < m-1; i++)
			for (j = 0; j < n; j++)
				for (k = 1; k < p-1; k++){
					x = i-dt*u[i][j][k]/dx; y = j-dt*v[i][j][k]/dy; z = k-dt*w[i][j][k]/dz;
					if (x<thr)
						x=thr;
					if (x>=m-1-thr)
						x=m-1-thr;
					if (y<thr)
						y = y+(n-1);
					if (y>=n-1-thr)
						y = y-(n-1);
					if (z<thr){ 
						if (bczlo == 0)
							z = z+(p-1);
						else
							z=thr;
					}
					if (z>=p-1-thr){
						if (bczhi == 0)
							z = z-(p-1);
						else
							z=p-1-thr;
					}
					i0=(int)x; i1=i0+1;
					j0=(int)y; j1=j0+1;
					k0=(int)z; k1=k0+1;
					s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1; r1 = z-k0; r0 = 1-r1;
					d[i][j][k] = s0*(t0*(r0*d0[i0][j0][k0]+r1*d0[i0][j0][k1])
							 +t1*(r0*d0[i0][j1][k0]+r1*d0[i0][j1][k1]))+
						  s1*(t0*(r0*d0[i1][j0][k0]+r1*d0[i1][j0][k1])
							 +t1*(r0*d0[i1][j1][k0]+r1*d0[i1][j1][k1]));
				}
		//z boundaries
		for ( i = 0; i < m; i++)
			for (j = 0; j < n; j++){
				d[i][j][0] = d0[i][j][1];
				d[i][j][p-1] = d0[i][j][p-2];
			}
		//x boundaries
		for ( k = 0; k < p; k++)
			for (j = 0; j < n; j++){
				d[0][j][k] = d0[1][j][k];
				d[m-1][j][k] = d0[m-2][j][k];
			}
	}
	else{
	}//the other periodic cases... y, and both x and y
	
}


//fGetOffset finds the approximate point of intersection of the surface
// between two points with the values fValue1 and fValue2
float Grid::GetOffset(float fValue1, float fValue2, float fValueDesired)
{
        float fDelta = fValue2 - fValue1;

        if(fDelta == 0.0)
        {
                return 0.5;
        }
        return (fValueDesired - fValue1)/fDelta;
}

void Grid::NormalizeVector(vect &rfVectorResult, vect &rfVectorSource)
{
        float fOldLength;
        float fScale;

        fOldLength = sqrtf( (rfVectorSource.fX * rfVectorSource.fX) +
                            (rfVectorSource.fY * rfVectorSource.fY) +
                            (rfVectorSource.fZ * rfVectorSource.fZ) );

        if(fOldLength == 0.0)
        {
                rfVectorResult.fX = rfVectorSource.fX;
                rfVectorResult.fY = rfVectorSource.fY;
                rfVectorResult.fZ = rfVectorSource.fZ;
        }
        else
        {
                fScale = 1.0/fOldLength;
                rfVectorResult.fX = rfVectorSource.fX*fScale;
                rfVectorResult.fY = rfVectorSource.fY*fScale;
                rfVectorResult.fZ = rfVectorSource.fZ*fScale;
        }
}


//MarchCube performs the Marching Cubes algorithm on a single cube
void Grid::MarchCube(int fX, int fY, int fZ, float fScale)
{
        extern int CubeEdgeFlags[256];
        extern int TriangleConnectionTable[256][16];

        int iCorner, iVertex, iVertexTest, iEdge, iTriangle, iFlagIndex, iEdgeFlags;
        float fOffset;
        vect sColor;
        float afCubeValue[8], afGradx[8], afGrady[8], afGradz[8];
        vect asEdgeVertex[12];
        vect asEdgeNorm[12];
		Cell *cell;

        //Make a local copy of the values at the cube's corners
        for(iVertex = 0; iVertex < 8; iVertex++)
        {
			int i = fX+(int)VertexOffset[iVertex][0];
			int j = fY+(int)VertexOffset[iVertex][1];
			int k = fZ+(int)VertexOffset[iVertex][2];
			afCubeValue[iVertex] = phi[i][j][k];
			
			if (i < 1 ) i = 1; if (i > m-2) i = m-2;
			if (j < 1 ) j = 1; if (j > n-2) j = n-2;
			if (k < 1 ) k = 1; if (k > p-2) k = p-2;
			afGradx[iVertex] = (phi[i+1][j][k]-phi[i-1][j][k])/2;
			afGrady[iVertex] = (phi[i][j+1][k]-phi[i][j-1][k])/2;
			afGradz[iVertex] = (phi[i][j][k+1]-phi[i][j][k-1])/2;
        }

        //Find which vertices are inside of the surface and which are outside
        iFlagIndex = 0;
        for(iVertexTest = 0; iVertexTest < 8; iVertexTest++)
        {
                if(afCubeValue[iVertexTest] <= 0.0f) 
                        iFlagIndex |= 1<<iVertexTest;
        }

        //Find which edges are intersected by the surface
        iEdgeFlags = CubeEdgeFlags[iFlagIndex];

        //If the cube is entirely inside or outside of the surface, then there will be no intersections
        if(iEdgeFlags == 0) 
        {
                return;
        }

        //Find the point of intersection of the surface with each edge
        //Then find the normal to the surface at those points
        for(iEdge = 0; iEdge < 12; iEdge++)
        {
                //if there is an intersection on this edge
                if(iEdgeFlags & (1<<iEdge))
                {
                        fOffset = GetOffset(afCubeValue[ EdgeConnection[iEdge][0] ], 
                                  afCubeValue[ EdgeConnection[iEdge][1] ], 0);

                        asEdgeVertex[iEdge].fX = (fX + VertexOffset[ EdgeConnection[iEdge][0] ][0]  +  fOffset * EdgeDirection[iEdge][0]) * fScale;
                        asEdgeVertex[iEdge].fY = (fY + VertexOffset[ EdgeConnection[iEdge][0] ][1]  +  fOffset * EdgeDirection[iEdge][1]) * fScale;
                        asEdgeVertex[iEdge].fZ = (fZ + VertexOffset[ EdgeConnection[iEdge][0] ][2]  +  fOffset * EdgeDirection[iEdge][2]) * fScale;
						
						fOffset = 1-fOffset;
						asEdgeNorm[iEdge].fX = fOffset*afGradx[EdgeConnection[iEdge][0]]+(1-fOffset)*afGradx[EdgeConnection[iEdge][1]];
						asEdgeNorm[iEdge].fY = fOffset*afGrady[EdgeConnection[iEdge][0]]+(1-fOffset)*afGrady[EdgeConnection[iEdge][1]];
						asEdgeNorm[iEdge].fZ = fOffset*afGradz[EdgeConnection[iEdge][0]]+(1-fOffset)*afGradz[EdgeConnection[iEdge][1]];
						NormalizeVector(asEdgeNorm[iEdge],asEdgeNorm[iEdge]);
						//  GetNormal(asEdgeNorm[iEdge], asEdgeVertex[iEdge].fX, asEdgeVertex[iEdge].fY, asEdgeVertex[iEdge].fZ);
                }
        }


		int ltrsize = 0;
        //Draw the triangles that were found.  There can be up to five per cube
        for(iTriangle = 0; iTriangle < 5; iTriangle++)
        {
                if(TriangleConnectionTable[iFlagIndex][3*iTriangle] < 0)
                        break;
				//get the triangle vertices and vertex normals
				iCorner = 0;
				iVertex = TriangleConnectionTable[iFlagIndex][3*iTriangle+iCorner];
				tr[trsize].na = asEdgeNorm[iVertex];
				tr[trsize].a = asEdgeVertex[iVertex];
				iCorner = 1;
				iVertex = TriangleConnectionTable[iFlagIndex][3*iTriangle+iCorner];
				tr[trsize].nb = asEdgeNorm[iVertex];
				tr[trsize].b = asEdgeVertex[iVertex];
				iCorner = 2;
				iVertex = TriangleConnectionTable[iFlagIndex][3*iTriangle+iCorner];
				tr[trsize].nc = asEdgeNorm[iVertex];
				tr[trsize].c = asEdgeVertex[iVertex];
				cell = &cell_list[fX][fY][fZ];
				cell->triangles[ltrsize++] = trsize; 
				trsize++;
        }
		cell->num_triangles = ltrsize;

}

//vMarchingCubes iterates over the entire dataset, calling vMarchCube on each cube
void Grid::MarchingCubes(void)
{
	trsize = 0;		
    int iX, iY, iZ;
    for(iX = 0; iX < m-1; iX++)
    for(iY = 0; iY < n-1; iY++)
    for(iZ = 0; iZ < p-1; iZ++)
    {
           // vMarchCube1(iX*fStepSize, iY*fStepSize, iZ*fStepSize, fStepSize);
			MarchCube(iX,iY,iZ,1);
    }

}

int CubeEdgeFlags[256]=
{
        0x000, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00, 
        0x190, 0x099, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c, 0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90, 
        0x230, 0x339, 0x033, 0x13a, 0x636, 0x73f, 0x435, 0x53c, 0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30, 
        0x3a0, 0x2a9, 0x1a3, 0x0aa, 0x7a6, 0x6af, 0x5a5, 0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0, 
        0x460, 0x569, 0x663, 0x76a, 0x066, 0x16f, 0x265, 0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60, 
        0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0x0ff, 0x3f5, 0x2fc, 0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0, 
        0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x055, 0x15c, 0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950, 
        0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0x0cc, 0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0, 
        0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0x0cc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0, 
        0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x055, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650, 
        0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc, 0x2fc, 0x3f5, 0x0ff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0, 
        0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c, 0x36c, 0x265, 0x16f, 0x066, 0x76a, 0x663, 0x569, 0x460, 
        0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac, 0x4ac, 0x5a5, 0x6af, 0x7a6, 0x0aa, 0x1a3, 0x2a9, 0x3a0, 
        0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c, 0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x033, 0x339, 0x230, 
        0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c, 0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x099, 0x190, 
        0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c, 0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x000
};

//  For each of the possible vertex states listed in aiCubeEdgeFlags there is a specific triangulation
//  of the edge intersection points.  a2iTriangleConnectionTable lists all of them in the form of
//  0-5 edge triples with the list terminated by the invalid value -1.
//  For example: a2iTriangleConnectionTable[3] list the 2 triangles formed when corner[0] 
//  and corner[1] are inside of the surface, but the rest of the cube is not.
//
//  I found this table in an example program someone wrote long ago.  It was probably generated by hand

int TriangleConnectionTable[256][16] =  
{
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
        {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
        {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
        {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
        {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
        {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
        {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
        {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
        {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
        {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
        {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
        {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
        {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
        {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
        {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
        {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
        {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
        {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
        {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
        {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
        {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
        {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
        {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
        {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
        {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
        {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
        {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
        {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
        {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
        {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
        {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
        {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
        {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
        {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
        {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
        {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
        {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
        {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
        {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
        {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
        {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
        {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
        {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
        {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
        {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
        {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
        {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
        {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
        {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
        {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
        {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
        {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
        {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
        {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
        {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
        {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
        {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
        {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
        {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
        {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
        {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
        {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
        {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
        {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
        {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
        {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
        {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
        {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
        {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
        {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
        {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
        {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
        {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
        {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
        {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
        {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
        {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
        {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
        {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
        {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
        {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
        {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
        {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
        {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
        {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
        {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
        {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
        {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
        {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
        {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
        {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
        {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
        {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
        {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
        {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
        {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
        {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
        {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
        {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
        {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
        {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
        {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
        {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
        {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
        {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
        {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
        {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
        {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
        {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
        {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
        {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
        {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
        {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
        {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
        {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
        {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
        {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
        {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
        {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
        {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
        {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
        {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
        {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
        {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
        {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
        {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
        {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
        {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
        {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
        {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
        {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
        {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
        {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
        {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
        {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
        {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
        {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
        {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
        {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
        {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
        {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
        {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
        {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
        {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
        {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
        {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
        {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
        {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
        {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
        {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
        {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
        {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
        {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
        {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
        {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
        {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
        {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
        {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
        {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
        {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
        {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
        {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
        {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
        {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
        {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
        {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
        {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
        {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
        {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
};
