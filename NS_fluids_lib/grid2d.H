#ifndef GRID2_H
#define GRID2_H

struct Point{
	Point(void) : x(0), y(0) {};
	~Point(){};
	float x,y;

	Point(float xx, float yy){
		x = xx;
		y = yy;
	}
	Point& operator	= ( const Point& pt );
//	Point Point::operator = (Point pt);
};

struct Cell{
	Cell();
	Cell(int ii,int jj, float dxx, float dyy);
	~Cell();
	int i,j;
	float phi;
	float phi_recon;

	Point *particles;
	bool *initial_particles;
	bool *initial_particles_new;
	Point *new_particles;
	int num_particles;
	int num_newparticles;

	int num_segmentpts,num_segmentpts1;
	Point *segments, *segments1;
	float **color;
	float **new_color;
};

class Grid{
public: 
	int m, n;
	float **phi;//the level set
	float **solidphi;
	float **den;	
	float **phi_recon;
	float **temp_phi1;
	float **temp_phi2;
	float ***color;
	float **u, **v;//velocity components
	float **F1, **F2;//force components	
	Point *s, *s1, *s2;//levelset points
	float **pointcolor;
	int ssize,ssize1,ssize2;
	float x_lo, x_hi, y_lo, y_hi;
	float dx,dy;
	float mass;
	int trsize, num_markers, num_free_markers;
	Cell **cell_list;
	Point *plist;
	int psize;
	float **dpx, **dpy, **dmx, **dmy;
	int bcxlo,bcxhi,bcylo,bcyhi;//boundary conditions
	int iter;

	Grid();
	Grid(int m,int n, float xx_lo, float xx_hi, float yy_lo, 
          float yy_hi, int bcxlo, int bcxhi,int bcylo, int bcyhi);
                       //, float xsize, float ysize);
	~Grid(){}
	void InitializeParticles(void);
	void FixLevelSetInside(float threshold, float eps);
	void FixLevelSetOutside(float threshold, float eps);
	bool CellChangesSignInside(Cell* c, float threshold);
	void Redistance0(void);
	void Redistance(void);
	void Redistance1(void);
	void Reinitialize(int iter, float dtau);
	void ReinitializeCIR(int iter, float dtau);
	bool RedistributeParticlesGrid(float threshold);
	void AddParticlesToCell(Cell* c, float threshold);
	float FindLineDist(Cell* d,  Point point, float threshold);
	void AdvectParticles(float dt);
	bool AdvectParticle(Cell* c, int k, float dt);
	void InterpolateCells2Grid(void);
	void InterpolateGrid2Cells(void);
	void InitializeProblem(char * problem_name);
	void update_levelsetDL(float dt);
	void update_levelsetSL(float dt);
	void FindZeroLevelCurve(float scalex, float scaley, float threshold, bool cellbased);
	void advect (int b, float ** d, float ** d0, float ** u, float ** v, float dt );
	void set_bnd (int b, float ** x );
	void DPlusx (float **dplus, float **u);
	void DPlusy (float **dplus, float **u);
	void DMinusx (float **dminus, float **u);
	void DMinusy (float **dminus, float **u);
	float Interpolate(float** field, int order, int i, int j, float ldx, float ldy);
	int CountParticles(void);
	void CorrectLevelSet(void);
	void ReinitializeRS(int iter, float dtau);
	void CheckReinit(int iter, float dtau);
	void ReinitializeV(int iter, float dtau);
	bool RedistributeParticles(float threshold);
	int CountInterfaceCells(void);
	
};

#endif // GRID2_H
