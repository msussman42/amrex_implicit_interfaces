#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>

#ifdef BL_USE_FLOAT
    typedef float Real;
#else
    typedef double Real;
#endif

typedef Real* pointertype;
typedef Real sizetype;


__global__ void gpu_jacobi(
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

unsigned int blockIdxz = blockIdx.y / blocksInY;
unsigned int blockIdxy = blockIdx.y % blocksInY;
unsigned int k = blockIdxz *blockDim.z + threadIdx.z;
unsigned int j = blockIdxy *blockDim.y + threadIdx.y;
unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
int idx = k*zstride + j*ystride + i;

 rhs[idx]=0.0;
 soln[idx]=0.0;
 red[idx]=0.0;
 black[idx]=0.0;
 Ax[idx]=0.0;
 __syncthreads();
 
 rhs[idx]=b[idx]-alphasing[idx]*x[idx]+
#if (BL_SPACEDIM==3)
  Uz[idx]*x[idx+zstride]+Lz[idx]*x[idx-zstride]+
#endif
  Ux[idx]*x[idx+1]+Lx[idx]*x[idx-1]+
  Uy[idx]*x[idx+ystride]+Ly[idx]*x[idx-ystride];

 __syncthreads();

 red[idx]=rhs[idx]/alpha[idx]; 
 __syncthreads();

 black[idx]=(rhs[idx]+
#if (BL_SPACEDIM==3)
  Uz[idx]*red[idx+zstride]+Lz[idx]*red[idx-zstride]+
#endif
  Ux[idx]*red[idx+1]+Lx[idx]*red[idx-1]+
  Uy[idx]*red[idx+ystride]+Ly[idx]*red[idx-ystride])/alpha[idx];
 __syncthreads();

 red[idx]=(rhs[idx]+
#if (BL_SPACEDIM==3)
  Uz[idx]*black[idx+zstride]+Lz[idx]*black[idx-zstride]+
#endif 
  Ux[idx]*black[idx+1]+Lx[idx]*black[idx-1]+
  Uy[idx]*black[idx+ystride]+Ly[idx]*black[idx-ystride])/alpha[idx];
 __syncthreads();

 soln[idx]=mask[idx]*red[idx]+(1.0-mask[idx])*black[idx];
 __syncthreads();

 Ax[idx]=x[idx]+soln[idx];
 __syncthreads();
}

extern "C" void kernel_gpu_jacobi(
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

      /* grid dimensions divisible by 8 */
    int blocksIny = NY/4;
    dim3 dimBlock(512,1,1);
    dim3 dimGrid(NX/8,NZ*NY/64,1);

    cudaMemcpy(xx_d, x, grid_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(bb_d, b, grid_bytes, cudaMemcpyHostToDevice);

    gpu_jacobi <<< dimGrid, dimBlock >>> (xx_d,bb_d,
      alpha_d,alphasing_d, 
#if (BL_SPACEDIM==3)
      Uz_d,Lz_d,
#endif
      Ux_d, Uy_d, 
      Lx_d, Ly_d, 
      zstride,ystride, blocksIny, ax_d,
      rhs_d,soln_d,red_d,black_d,mask_d,
      N,NX,NY,NZ);

    cudaMemcpy(x, ax_d, grid_bytes, cudaMemcpyDeviceToHost);
}


__global__ void gpu_apply(
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
unsigned int blockIdxz = blockIdx.y / blocksInY;
unsigned int blockIdxy = blockIdx.y % blocksInY;
unsigned int k = blockIdxz *blockDim.z + threadIdx.z;
unsigned int j = blockIdxy *blockDim.y + threadIdx.y;
unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
int idx = k*zstride + j*ystride + i;

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


extern "C" void kernel_gpu_apply(
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

      /* grid dimensions divisible by 8 */
    int blocksIny = NY/4;
    dim3 dimBlock(512,1,1);
    dim3 dimGrid(NX/8,NZ*NY/64,1);

    cudaMemcpy(xx_d, x, grid_bytes, cudaMemcpyHostToDevice);

    gpu_apply <<< dimGrid, dimBlock >>> (xx_d,alpha_d, 
#if (BL_SPACEDIM==3)
      Uz_d,Lz_d,
#endif
      Ux_d, Uy_d, 
      Lx_d, Ly_d, 
      zstride,ystride, blocksIny, ax_d,
      N,NX,NY,NZ);

    cudaMemcpy(Ax, ax_d, grid_bytes, cudaMemcpyDeviceToHost);
}



// Function 1: Copy arrays to GPU and leave them there //////////////////////

extern "C" void kernel_wrapper(
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
 pointertype* Lx_d, pointertype* Ly_d) {

   long int grid_bytes=sizeof(sizetype)*N;

   // Allocate memory on GPU
   // if parameter passed as pointertype&, then (void **)&xx_d
   cudaMalloc( (void **) xx_d, grid_bytes );
   cudaMalloc( (void **) bb_d, grid_bytes );
   cudaMalloc( (void **) ax_d, grid_bytes );

   cudaMalloc( (void **) rhs_d, grid_bytes );
   cudaMalloc( (void **) soln_d, grid_bytes );
   cudaMalloc( (void **) red_d, grid_bytes );
   cudaMalloc( (void **) black_d, grid_bytes );
   cudaMalloc( (void **) mask_d, grid_bytes );

   cudaMalloc( (void **) alpha_d, grid_bytes );
   cudaMalloc( (void **) alphasing_d, grid_bytes );
#if (BL_SPACEDIM==3)
   cudaMalloc( (void **) Uz_d, grid_bytes );
   cudaMalloc( (void **) Lz_d, grid_bytes );
#endif
   cudaMalloc( (void **) Ux_d, grid_bytes );
   cudaMalloc( (void **) Uy_d, grid_bytes );
   cudaMalloc( (void **) Lx_d, grid_bytes );
   cudaMalloc( (void **) Ly_d, grid_bytes );

   // copy arrays from CPU to GPU
   cudaMemcpy(*mask_d, mask, grid_bytes, cudaMemcpyHostToDevice);
   cudaMemcpy(*alpha_d, alpha, grid_bytes, cudaMemcpyHostToDevice);
   cudaMemcpy(*alphasing_d,alphasing,grid_bytes,cudaMemcpyHostToDevice);
#if (BL_SPACEDIM==3)
   cudaMemcpy(*Uz_d, Uz, grid_bytes, cudaMemcpyHostToDevice);
   cudaMemcpy(*Lz_d, Lz, grid_bytes, cudaMemcpyHostToDevice);
#endif
   cudaMemcpy(*Ux_d, Ux, grid_bytes, cudaMemcpyHostToDevice);
   cudaMemcpy(*Uy_d, Uy, grid_bytes, cudaMemcpyHostToDevice);
   cudaMemcpy(*Lx_d, Lx, grid_bytes, cudaMemcpyHostToDevice);
   cudaMemcpy(*Ly_d, Ly, grid_bytes, cudaMemcpyHostToDevice);

  return;
}


// Function 3:  Free GPU /////////////////////////////////////////////////////

extern "C" void free_device(
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
 cudaFree(*xx_d);
 cudaFree(*bb_d);
 cudaFree(*ax_d);
 cudaFree(*rhs_d);
 cudaFree(*soln_d);
 cudaFree(*red_d);
 cudaFree(*black_d);
 cudaFree(*mask_d);
 cudaFree(*alpha_d);
 cudaFree(*alphasing_d);
#if (BL_SPACEDIM==3)
 cudaFree(*Uz_d);
 cudaFree(*Lz_d);
#endif
 cudaFree(*Ux_d);
 cudaFree(*Uy_d);
 cudaFree(*Lx_d);
 cudaFree(*Ly_d);
}

