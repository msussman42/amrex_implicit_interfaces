#ifndef CNS_INDEX_H_
#define CNS_INDEX_H_

#define N_CONSERVATIVE 2
#define N_CHARACTERISTIC 2

#if (N_CONSERVATIVE==2)
#define PVAR 1
#define UVAR 2
#define W1 1
#define W2 2

#define U0 1.0
#define c  1.0
#define rho0 1.0
#define lambda1 (U0+c)
#define lambda2 (U0-c)

#define X11 ...
#define X21 ...
#define X12 ...
#define X22 ...

#define XINV11 ...
#define XINV21 ...
#define XINV12 ...
#define XINV22 ...
#elif (N_CONSERVATIVE==1)
#define UVAR 1
#define W1 1
#define lambda1 1.0
#define X11 1.0
#define XINV11 1.0
#else
abort
#endif
