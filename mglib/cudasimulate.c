#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef BL_USE_FLOAT
    typedef float Real;
#else
    typedef double Real;
#endif

typedef Real* pointertype;
typedef Real sizetype;


void gpu_jacobiCPU(
pointertype x, 
pointertype b,
pointertype alpha, 
pointertype alphasing, 
#if (BL_SPACEDIM==3)
pointertype Uz,pointertype Lz,
#endif
pointertype Ux,pointertype Uy,  
pointertype Lx, pointertype Ly, 
int zstride, int ystride, int blocksInY, pointertype Ax,
pointertype rhs,pointertype soln,pointertype red,pointertype black,
pointertype mask,long int N,int NX,int NY,int NZ) {

int idx;

int buffer=NX;
#if (BL_SPACEDIM==3)
buffer*=NY;
#endif

for (idx=0;idx<N;idx++) {
 rhs[idx]=0.0;
 soln[idx]=0.0;
 red[idx]=0.0;
 black[idx]=0.0;
 Ax[idx]=0.0;
}

for (idx=buffer;idx<N-buffer;idx++) {
 rhs[idx]=b[idx]-alphasing[idx]*x[idx]+
#if (BL_SPACEDIM==3)
  Uz[idx]*x[idx+zstride]+Lz[idx]*x[idx-zstride]+
#endif
  Ux[idx]*x[idx+1]+Lx[idx]*x[idx-1]+
  Uy[idx]*x[idx+ystride]+Ly[idx]*x[idx-ystride];
}

for (idx=buffer;idx<N-buffer;idx++) {
 red[idx]=rhs[idx]/alpha[idx]; 
}

for (idx=buffer;idx<N-buffer;idx++) {
 black[idx]=(rhs[idx]+
#if (BL_SPACEDIM==3)
  Uz[idx]*red[idx+zstride]+Lz[idx]*red[idx-zstride]+
#endif
  Ux[idx]*red[idx+1]+Lx[idx]*red[idx-1]+
  Uy[idx]*red[idx+ystride]+Ly[idx]*red[idx-ystride])/alpha[idx];
}

for (idx=buffer;idx<N-buffer;idx++) {
 red[idx]=(rhs[idx]+
#if (BL_SPACEDIM==3)
  Uz[idx]*black[idx+zstride]+Lz[idx]*black[idx-zstride]+
#endif 
  Ux[idx]*black[idx+1]+Lx[idx]*black[idx-1]+
  Uy[idx]*black[idx+ystride]+Ly[idx]*black[idx-ystride])/alpha[idx];
}

for (idx=buffer;idx<N-buffer;idx++) {
 soln[idx]=mask[idx]*red[idx]+(1.0-mask[idx])*black[idx];
}

for (idx=buffer;idx<N-buffer;idx++) {
 Ax[idx]=x[idx]+soln[idx];
}

}

void kernel_gpu_jacobiCPU(
 pointertype x, 
 pointertype b, 
 pointertype xx_d,
 pointertype bb_d,
 pointertype ax_d,
 pointertype rhs_d,
 pointertype soln_d,
 pointertype red_d,
 pointertype black_d,
 pointertype mask_d,
 pointertype alpha_d,
 pointertype alphasing_d,
#if (BL_SPACEDIM==3)
 pointertype Uz_d,pointertype Lz_d,
#endif
 pointertype Ux_d, pointertype Uy_d, 
 pointertype Lx_d, pointertype Ly_d, 
 long int Np, int NX, int NY, int NZ) {

    long int N=Np;
    int zstride=NX*NY;
    int ystride=NX;
    long int grid_bytes=sizeof(sizetype)*N;
    int blocksIny=0;

    memcpy(xx_d, x, grid_bytes);
    memcpy(bb_d, b, grid_bytes);

    gpu_jacobiCPU(xx_d,bb_d,
      alpha_d,alphasing_d, 
#if (BL_SPACEDIM==3)
      Uz_d,Lz_d,
#endif
      Ux_d, Uy_d, 
      Lx_d, Ly_d, 
      zstride,ystride, blocksIny, ax_d,
      rhs_d,soln_d,red_d,black_d,mask_d,
      N,NX,NY,NZ);

    memcpy(x, ax_d, grid_bytes);
}


void gpu_applyCPU(
pointertype x, 
pointertype alpha, 
#if (BL_SPACEDIM==3)
pointertype Uz,pointertype Lz,
#endif
pointertype Ux,pointertype Uy, 
pointertype Lx, pointertype Ly, 
int zstride, int ystride, int blocksInY, pointertype Ax,
long int N,int NX,int NY,int NZ)
{

int idx;
int buffer=NX;
#if (BL_SPACEDIM==3)
buffer*=NY;
#endif

for (idx=buffer;idx<N-buffer;idx++) {
 Ax[idx]=alpha[idx]*x[idx]-(
#if (BL_SPACEDIM==3)
         Uz[idx] * x[idx+zstride]+
         Lz[idx] * x[idx-zstride]+
#endif
         Ux[idx] * x[idx+1]+
         Lx[idx] * x[idx-1]+
         Uy[idx] * x[idx+ystride]+
         Ly[idx] * x[idx-ystride]);
}
 
}


void kernel_gpu_applyCPU(
 pointertype x, 
 pointertype xx_d,
 pointertype ax_d,
 pointertype alpha_d,
#if (BL_SPACEDIM==3)
 pointertype Uz_d,pointertype Lz_d,
#endif
 pointertype Ux_d, pointertype Uy_d, 
 pointertype Lx_d, pointertype Ly_d, 
 long int Np, pointertype Ax, int NX, int NY, int NZ) {

    long int N=Np;
    int zstride=NX*NY;
    int ystride=NX;
    long int grid_bytes=sizeof(sizetype)*N;
    int blocksIny=0;

    memcpy(xx_d, x, grid_bytes);

    gpu_applyCPU(xx_d,alpha_d, 
#if (BL_SPACEDIM==3)
      Uz_d,Lz_d,
#endif
      Ux_d, Uy_d, 
      Lx_d, Ly_d, 
      zstride,ystride, blocksIny, ax_d,
      N,NX,NY,NZ);

    memcpy(Ax, ax_d, grid_bytes);
}



// Function 1: Copy arrays to GPU and leave them there //////////////////////

void kernel_wrapperCPU(
 pointertype mask, 
 pointertype alpha, 
 pointertype alphasing, 
#if (BL_SPACEDIM==3)
 pointertype Uz,pointertype Lz,
#endif
 pointertype Ux, pointertype Uy, 
 pointertype Lx, pointertype Ly, 
 long int N, 
 pointertype* xx_d, 
 pointertype* bb_d, 
 pointertype* ax_d, 
 pointertype* rhs_d,
 pointertype* soln_d,
 pointertype* red_d,
 pointertype* black_d,
 pointertype* mask_d,
 pointertype* alpha_d, 
 pointertype* alphasing_d, 
#if (BL_SPACEDIM==3)
 pointertype* Uz_d,pointertype* Lz_d,
#endif
 pointertype* Ux_d, pointertype* Uy_d,  
 pointertype* Lx_d, pointertype* Ly_d)
{

   long int grid_bytes=sizeof(sizetype)*N;

   // Allocate memory on CPU
   *xx_d=(pointertype)malloc(grid_bytes);
   *bb_d=(pointertype)malloc(grid_bytes);
   *ax_d=(pointertype)malloc(grid_bytes);

   *rhs_d=(pointertype)malloc(grid_bytes);
   *soln_d=(pointertype)malloc(grid_bytes);
   *red_d=(pointertype)malloc(grid_bytes);
   *black_d=(pointertype)malloc(grid_bytes);
   *mask_d=(pointertype)malloc(grid_bytes);

   *alpha_d=(pointertype)malloc(grid_bytes);
   *alphasing_d=(pointertype)malloc(grid_bytes);
#if (BL_SPACEDIM==3)
   *Uz_d=(pointertype)malloc(grid_bytes);
   *Lz_d=(pointertype)malloc(grid_bytes);
#endif
   *Ux_d=(pointertype)malloc(grid_bytes);
   *Uy_d=(pointertype)malloc(grid_bytes);
   *Lx_d=(pointertype)malloc(grid_bytes);
   *Ly_d=(pointertype)malloc(grid_bytes);

   // copy arrays from CPU to CPU
   memcpy(*mask_d, mask, grid_bytes);
   memcpy(*alpha_d, alpha, grid_bytes);
   memcpy(*alphasing_d,alphasing,grid_bytes);
#if (BL_SPACEDIM==3)
   memcpy(*Uz_d, Uz, grid_bytes);
   memcpy(*Lz_d, Lz, grid_bytes);
#endif
   memcpy(*Ux_d, Ux, grid_bytes);
   memcpy(*Uy_d, Uy, grid_bytes);
   memcpy(*Lx_d, Lx, grid_bytes);
   memcpy(*Ly_d, Ly, grid_bytes);

  return;
}


// Function 3:  Free CPU /////////////////////////////////////////////////////

void free_deviceCPU(
pointertype* xx_d, 
pointertype* bb_d, 
pointertype* ax_d, 
pointertype* rhs_d,
pointertype* soln_d,
pointertype* red_d,
pointertype* black_d,
pointertype* mask_d,
pointertype* alpha_d, 
pointertype* alphasing_d, 
#if (BL_SPACEDIM==3)
pointertype* Uz_d,pointertype* Lz_d,
#endif
pointertype* Ux_d, pointertype* Uy_d, 
pointertype* Lx_d, pointertype* Ly_d)
{
 free(*xx_d);
 free(*bb_d);
 free(*ax_d);
 free(*rhs_d);
 free(*soln_d);
 free(*red_d);
 free(*black_d);
 free(*mask_d);
 free(*alpha_d);
 free(*alphasing_d);
#if (BL_SPACEDIM==3)
 free(*Uz_d);
 free(*Lz_d);
#endif
 free(*Ux_d);
 free(*Uy_d);
 free(*Lx_d);
 free(*Ly_d);
}

