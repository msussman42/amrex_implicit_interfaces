//============================================================================

// $Author: erlebach $
// $Date: 2002/01/05 15:18:28 $
// $Header: /home/erlebach/D2/src/CVS/computeCP/WriteAmiraMesh.cpp,v 2.2 2002/01/05 15:18:28 erlebach Exp $
// Sticky tag, $Name:  $
// $RCSfile: WriteAmiraMesh.cpp,v $
// $Revision: 2.2 $
// $State: Exp $

//==========================================================================


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "WriteAmiraMesh.h"
#include "MyEndian.h"


WriteAmiraMesh::WriteAmiraMesh(int nx_, int ny_, int nz_,int nbComponents_)
   : nx(nx_), ny(ny_), nz(nz_), nbComponents(nbComponents_)
{
    if (sizeof(MSLONG)!=4) {
     printf("long does not have a size of 4 bytes");
     exit(1);
    }
    ntot = nx*ny*nz;
}

void WriteAmiraMesh::writeRectilinearData(float* vx, float* vy, float* vz,
  char* name,const char* fileName_,float* x,float* y,float* z,
  int viorel_format)
{

    if ((viorel_format!=0)&&(viorel_format!=1)) {
     printf("viorel_format invalid\n");
     exit(1);
    }

    if (sizeof(MSLONG)!=4) {
     printf("long does not have a size of 4 bytes");
     exit(1);
    }
    MyEndian myendian;
    
    if (viorel_format==1)
     myendian.disableConvert();

    printf("writeRectilinearData, filename= %s\n", fileName_);
    FILE* fd = fopen(fileName_, "w");
    if (!fd) {
        printf("file %s cannot be opened for writing\n", fileName_);
        exit(1);
    }
    if (viorel_format==0) {
     fprintf(fd, "# AmiraMesh BINARY\n");
     fprintf(fd, "# Dimensions in x, y, z, directions\n");
     fprintf(fd, "define Lattice %d %d %d\n", nx, ny, nz);
     fprintf(fd, "define Coordinates %d \n", nx+ny+nz);
     fprintf(fd, "Parameters {\n");
     fprintf(fd, " CoordType \"rectilinear\",\n");
     fprintf(fd, " }\n");
     fprintf(fd, "Coordinates { float xyz } = @1\n");
     fprintf(fd, "Lattice { float %s } = @2\n", name);
    }

    if (viorel_format==0) 
     fprintf(fd, "\n@1\n");

    myendian.fwrite(x, sizeof(float), nx, fd);
    myendian.fwrite(y, sizeof(float), ny, fd);
    myendian.fwrite(z, sizeof(float), nz, fd);

    ntot = nx*ny*nz;

    if (viorel_format==0) 
     fprintf(fd, "\n@2\n");

    printf("ntot= %d\n", ntot);
    printf("writeRectilinearData (vx,vy,vz) \n");

    {for (MSLONG i=0; i < ntot; i++) {
        myendian.fwrite(vx+i, sizeof(float), 1, fd);
        myendian.fwrite(vy+i, sizeof(float), 1, fd);
        myendian.fwrite(vz+i, sizeof(float), 1, fd);
    }}
    fclose(fd);
}


void WriteAmiraMesh::writeRectilinearData(float* scalar, char* name,
  const char* fileName_,float* x,float* y,float* z,
  int viorel_format)
{
    if ((viorel_format!=0)&&(viorel_format!=1)) {
     printf("viorel_format invalid\n");
     exit(1);
    }
    if (sizeof(MSLONG)!=4) {
     printf("long does not have a size of 4 bytes");
     exit(1);
    }
    MyEndian myendian;
    if (viorel_format==1)
     myendian.disableConvert();

    printf("writeRectilinearData, filename= %s\n", fileName_);
    FILE* fd = fopen(fileName_, "w");
    if (!fd) {
        printf("file %s cannot be opened for writing\n", fileName_);
        exit(1);
    }
    if (viorel_format==0) {
     fprintf(fd, "# AmiraMesh BINARY\n");
     fprintf(fd, "# Dimensions in x, y, z, directions\n");
     fprintf(fd, "define Lattice %d %d %d\n", nx, ny, nz);
     fprintf(fd, "define Coordinates %d \n", nx+ny+nz);
     fprintf(fd, "Parameters {\n");
     fprintf(fd, " CoordType \"rectilinear\",\n");
     fprintf(fd, " }\n");
     fprintf(fd, "Coordinates { float xyz } = @1\n");
     fprintf(fd, "Lattice { float %s } = @2\n", name);
    }

    if (viorel_format==0) 
     fprintf(fd, "\n@1\n");

    myendian.fwrite(x, sizeof(float), nx, fd);
    myendian.fwrite(y, sizeof(float), ny, fd);
    myendian.fwrite(z, sizeof(float), nz, fd);

    if (viorel_format==0) 
     fprintf(fd, "\n@2\n");

    ntot = nx*ny*nz;

    printf("ntot= %d\n", ntot);
    printf("writeRectilinearData (scalar) \n");

    myendian.fwrite(scalar, sizeof(float), ntot, fd);

    fclose(fd);
}



void WriteAmiraMesh::readRectilinearData(float* vx, float* vy, float* vz,
  char* name,const char* fileName_,float* x,float* y,float* z,
  int viorel_format)
{

    if ((viorel_format!=0)&&(viorel_format!=1)) {
     printf("viorel_format invalid\n");
     exit(1);
    }
    if (sizeof(MSLONG)!=4) {
     printf("long does not have a size of 4 bytes");
     exit(1);
    }

    MyEndian myendian;

    if (viorel_format==0) {
     printf("this is only for viorel's output format\n");
     exit(1);
    }

    if (viorel_format==1)
     myendian.disableConvert();

    printf("readRectilinearData, filename= %s\n", fileName_);
    FILE* fd = fopen(fileName_, "r");
    if (!fd) {
        printf("file %s cannot be opened for reading\n", fileName_);
        exit(1);
    }

    myendian.fread(x, sizeof(float), nx, fd);
    myendian.fread(y, sizeof(float), ny, fd);
    myendian.fread(z, sizeof(float), nz, fd);

    ntot = nx*ny*nz;

    printf("ntot= %d\n", ntot);

    printf("readRectilinearData (vx,vy,vz) \n");

    {for (MSLONG i=0; i < ntot; i++) {
        myendian.fread(vx+i, sizeof(float), 1, fd);
        myendian.fread(vy+i, sizeof(float), 1, fd);
        myendian.fread(vz+i, sizeof(float), 1, fd);
    }}
    fclose(fd);
}


void WriteAmiraMesh::readRectilinearData(float* scalar, char* name,
  const char* fileName_,float* x,float* y,float* z,
  int viorel_format)
{
    if ((viorel_format!=0)&&(viorel_format!=1)) {
     printf("viorel_format invalid\n");
     exit(1);
    }
    if (sizeof(MSLONG)!=4) {
     printf("long does not have a size of 4 bytes");
     exit(1);
    }

    MyEndian myendian;

    if (viorel_format==0) {
     printf("this is only for viorel's output format\n");
     exit(1);
    }

    if (viorel_format==1)
     myendian.disableConvert();

    printf("readRectilinearData, filename= %s\n", fileName_);
    FILE* fd = fopen(fileName_, "r");
    if (!fd) {
        printf("file %s cannot be opened for reading\n", fileName_);
        exit(1);
    }

    myendian.fread(x, sizeof(float), nx, fd);
    myendian.fread(y, sizeof(float), ny, fd);
    myendian.fread(z, sizeof(float), nz, fd);

    ntot = nx*ny*nz;

    printf("ntot= %d\n", ntot);

    printf("readRectilinearData (scalar) \n");

    myendian.fread(scalar, sizeof(float), ntot, fd);

    fclose(fd);
}

