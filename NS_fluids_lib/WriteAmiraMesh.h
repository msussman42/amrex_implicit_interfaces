//================================================================================
 
// $Author: erlebach $
// $Date: 2002/01/05 15:18:28 $
// $Header: /home/erlebach/D2/src/CVS/computeCP/WriteAmiraMesh.h,v 2.2 2002/01/05 15:18:28 erlebach Exp $
// Sticky tag, $Name:  $
// $RCSfile: WriteAmiraMesh.h,v $
// $Revision: 2.2 $
// $State: Exp $
 
//================================================================================

#ifndef _WRITEAMIRAMESH_H_
#define _WRITEAMIRAMESH_H_

#ifdef GORDON_FOURBYTEINT
#define MSLONG int
#else
#define MSLONG long
#endif

#include <string.h>

class WriteAmiraMesh {
private:
	int nx, ny, nz;
	MSLONG ntot;
	int nbComponents;

public:
	WriteAmiraMesh() {;}
	WriteAmiraMesh(int nx_, int ny_, int nz_, int nbComponents_);
	MSLONG getPts() {return ntot;}

        /// Rectilinear grid
        void writeRectilinearData(float* vx, char* name,const char* fileName_,
          float* x,float* y,float* z,int viorel_format=0);
        void writeRectilinearData(float* vx, float* vy, float* vz, 
          char* name,const char* fileName_,
          float* x,float* y,float* z,int viorel_format=0);
        void readRectilinearData(float* vx, char* name,const char* fileName_,
          float* x,float* y,float* z,int viorel_format);
        void readRectilinearData(float* vx, float* vy, float* vz, 
          char* name,const char* fileName_,
          float* x,float* y,float* z,int viorel_format);


};

#endif
