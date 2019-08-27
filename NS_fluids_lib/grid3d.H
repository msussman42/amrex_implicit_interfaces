#ifndef GRID3_H
#define GRID3_H

struct Point{
	Point(void) : x(0), y(0), z(0) {};
	float x,y,z;

	Point(float xx, float yy, float zz){
		x = xx;
		y = yy;
		z = zz;
	}
	Point& operator	= ( const Point& pt );
//	Point Point::operator = (Point pt);
};


struct vect {
    vect(void) : fX(0), fY(0), fZ(0) {};
    float fX,fY,fZ;     
    vect(float xx, float yy, float zz){
	fX = xx;
	fY = yy;
	fZ = zz;
    }
};

struct Triangle{
	vect a,b,c,na,nb,nc;
	Triangle(void ) {};
	Triangle(vect aa, vect bb, vect cc, 
          vect nna, vect nnb, vect nnc){
		a = aa;
		b = bb;
		c = cc;
		na = nna;
		nb = nnb;
		nc = nnc;
	}
	Triangle& operator	= ( const Triangle& tri );
};

struct Marker{
	Point p;
	vect color;
	float transp;
	Marker(void){};
	Marker(Point pp, vect ccolor, float ttransp){
		p = pp;
		color = ccolor;
		transp = ttransp;
	}
	Marker& operator	= ( const Marker& mk );
	~Marker(void){};
};

struct Cell{
	Cell();
	Cell(int ii,int jj, int kk, float dxx, float dyy, float dzz);
	~Cell();
	int i,j,k;
	float phi;
	float phi_recon;
	int num_particles;
	int num_triangles;
	int *triangles;//triangle list pointers
	int *markers;
};

class Grid{
public: 
	int m, n, p;
	float ***phi, ***solidphi, ***temp_phi, ***phi_recon, ***den;//the level set, etc
	float ***u, ***v, ***w, ***F1, ***F2, ***F3, ****color;//velocity and force components and color
	float x_lo, x_hi, y_lo, y_hi, z_lo, z_hi, dx, dy, dz;
	Triangle *tr;//triangle list
	int trsize, num_markers, num_free_markers;
	Cell ***cell_list;
	Marker *marker_list;
	Marker *free_markers;
	int bcxlo,bcxhi,bcylo,bcyhi,bczlo,bczhi;//boundary conditions
	int iter;

	Grid();
	Grid(int m,int n, int p, float x_lo, float x_hi, float y_lo, float y_hi, float z_lo, float z_hi, int bcxlo, int bcxhi,int bcylo, int bcyhi, int bczlo, int bczhi);
	~Grid(){}
	void initialize_problem(char * problem_name);
	void CreateSphere(float xc, float yc, float zc, float radius);
	vect BoilingColor(Point pp, Cell *cl,int num_mark);
	void InitializeParticles(void);
	void AssignParticleTexture(int type, char * pxls, int usize, int vsize);
	void OutputMarkers(char *filename);
	void OutputFreeMarkers(char *filename);
	void OutputColor(char *filename);
	void StopTopVel(void);
	vect InterpVelocity(Point pp);
	void AdvectMarkers(float dt);
	bool AdvectMarker(int i, int dropped, float dt);
	void AdvectFreeMarkers(float dt);
	bool AdvectFreeMarker(int i, int dropped, float dt);
	void RebuildCellMarkerListPointers(void);
	bool CellChangesSignInside(Cell* c, float threshold);
	bool InsideSolid(Cell* c, float threshold);
    int CountFullCells(void);
	int CountNonemptyCells(void);
	void OutputLevelSet(char *filename);
	void OutputSolidLevelSet(char *filename);
	void OutputWeber(char *filename);
	void OutputDensity(char *filename);
	void OutputVelocity(char *filename);
	void OutputVelocityMag(char *filename);
	void OutputFVelocity(char *filename);
	void WriteMeshFile(char* filename);
	void update_levelsetSL(float dt);
	void update_levelsetDL(float dt);
	float Interpolate(float*** field, int order, int i, int j, int k, float ldx, float ldy, float ldz);
	void Redistance1(void);
	void Reinit(int itern, float dtau);
	void DeleteExtraMarkers(void);
	bool AddParticles(int boiling);
	void AddParticlesToCell(Cell* c, int num_mark, int boiling);
	vect EstimateColor(Point pp, Cell *cl, int num_mark);
	void BounceMarker(float& px, float& py, float& pz, float& newvelx, float& newvely, float& newvelz);
	bool ChangeDensity(float level);
	void UpdateGasLevel(float level);
	void UpdateGeneralLevel(float level);
	bool WeberThreshold(Cell* c, float web);
	float VelocityVariation(Cell* c);
	float AbsCurvature(Cell* c);
	bool NotBoundaryCell(Cell* c);
	bool InsideDomain(Point* pp);

	void set_bnd (int b, float *** x );
	void advect (int b, float *** d, float *** d0, float *** u, float *** v, float *** w, float dt );
	float GetOffset(float fValue1, float fValue2, float fValueDesired);
	void MarchingCubes(void);
	void MarchCube(int fX, int fY, int fZ, float fScale);
	void NormalizeVector(vect &rfVectorResult, 
                             vect &rfVectorSource);
	void GetNormal(vect &rfNormal, float fX, float fY, float fZ);
};

#endif // GRID3_H
