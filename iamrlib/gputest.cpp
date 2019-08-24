#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>

typedef float* pointertype;
typedef float sizetype;


extern "C" void kernel_wrapper(pointertype& alpha, pointertype& Ux, 
 pointertype& Uy, pointertype& Uz, pointertype& Lx, pointertype& Ly,
 pointertype& Lz, long int& Np, pointertype& alpha_d, pointertype& Ux_d,
 pointertype& Uy_d, pointertype& Uz_d, pointertype& Lx_d, pointertype& Ly_d,
 pointertype& Lz_d);

extern "C" void kernel_wrapper_add(pointertype& x, pointertype& alpha_d, 
 pointertype& Ux_d, pointertype& Uy_d, pointertype& Uz_d,
 pointertype& Lx_d, pointertype& Ly_d,  pointertype& Lz_d, 
 long int& Np, pointertype& Ax, int& NX, int& NY, int& NZ);

extern "C" void free_device(pointertype& alpha_d, pointertype& Ux_d,
 pointertype& Uy_d, pointertype& Uz_d, pointertype& Lx_d, pointertype& Ly_d, 
 pointertype& Lz_d);


void cpu_kernel(pointertype x, pointertype alpha, pointertype Ux,
pointertype Uy, pointertype Uz, pointertype& Lx, pointertype Ly, pointertype Lz,
int zstride, int ystride, long int N, pointertype Ax)
{
  int idx;
  for (idx=zstride; idx<N-zstride; idx++){
    Ax[idx] = alpha[idx]*x[idx] + Ux[idx] * x[idx+1]
                             + Lx[idx] * x[idx-1]
                             + Uy[idx] * x[idx+ystride]
                             + Ly[idx] * x[idx-ystride]
                             + Uz[idx] * x[idx+zstride]
                             + Lz[idx] * x[idx-zstride];
  }
}


int main(void)
{
   int i, k, zstride, ystride;
   int NX=32, NY=32, NZ=64;
   long int Np=65536;
   pointertype x, alpha, alpha_d, Ax;
   pointertype Ux, Uy, Uz, Lx, Ly, Lz;
   pointertype Ux_d, Uy_d, Uz_d, Lx_d, Ly_d, Lz_d;

   alpha=(pointertype)malloc(sizeof(sizetype)*Np);
   x=(pointertype)malloc(sizeof(sizetype)*Np);
   Ux=(pointertype)malloc(sizeof(sizetype)*Np);
   Uy=(pointertype)malloc(sizeof(sizetype)*Np);
   Uz=(pointertype)malloc(sizeof(sizetype)*Np);
   Lx=(pointertype)malloc(sizeof(sizetype)*Np);
   Ly=(pointertype)malloc(sizeof(sizetype)*Np);
   Lz=(pointertype)malloc(sizeof(sizetype)*Np);
   Ax=(pointertype)malloc(sizeof(sizetype)*Np);

// indexing => A[i][j][k]=A[i+j*(NX-1)+k*(NX-1)*(NY-1)]

   zstride=NX*NY;
   ystride=NX;
//   Np=NX*NY*NZ;

   for (i=0; i<Np; i++){
    alpha[i]=8.0;
    x[i]=1.0*i;
    Ux[i]=-1.0;
    Uy[i]=-2.0;
    Uz[i]=-3.0;
    Lx[i]=4.0;
    Ly[i]=5.0;
    Lz[i]=6.0;
    Ax[i]=0.0;
   }

/*
   for (k=0; k<100; k++){
    cpu_kernel (x, alpha, Ux, Uy, Uz, Lx, Ly, Lz, zstride,
                 ystride, Np, Ax);
   }
*/


/*
   printf("Ax computed by cpu\n");
   for (i=0; i<Np; i++){
     printf("%d %f\n", i, Ax[i]);
   }
*/

   printf("Copying to GPU\n");


     kernel_wrapper(alpha, Ux, Uy, Uz, Lx, Ly, Lz, Np, alpha_d, 
                     Ux_d, Uy_d, Uz_d, Lx_d, Ly_d, Lz_d);


   printf("Alpha Copy Complete, Begin Vector Add\n");

   for (k=0; k<100; k++){
     kernel_wrapper_add(x, alpha_d, Ux_d, Uy_d, Uz_d, Lx_d, Ly_d, Lz_d,
                         Np, Ax, NX, NY, NZ);
   }



   printf("Kernel Complete\n");


     free_device(alpha_d, Ux_d, Uy_d, Uz_d, Lx_d, Ly_d, Lz_d);


/*   printf("Ax in main program finish\n");
   for (i=0; i<Np; i++){
     printf("%d %f\n", i, Ax[i]);
   }
*/

  return 0;
}

