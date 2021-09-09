#ifndef MYFUNC_F_H_
#define MYFUNC_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

using std::vector;

    void init_phi(const int* lo, const int* hi,
                  amrex_real* data, const int* dlo, const int* dhi,
                  const amrex_real* dx, const amrex_real* prob_lo, const amrex_real* prob_hi);
				  
	void init_vel(const int* lo, const int* hi,
				  const int* domlo, const int* domhi,
				  amrex_real* vel_u, const int* ulo, const int* uhi,
				  amrex_real* vel_v, const int* vlo, const int* vhi,
#if (BL_SPACEDIM == 3)   
                  amrex_real* vel_w, const int* wlo, const int* whi,
#endif
				  const amrex_real* dx, const amrex_real* prob_lo, const amrex_real* prob_hi, const amrex_real* time);
    
	//void NarrowBand(const int* lo, const int* hi, 
	//				amrex_real* phi, const int* philo, const int* phihi,
	//				const amrex_real* dx, 
	//				int* mask, const int* masklo, const int* maskhi,
	//				int* mI1, int* mI2, int* mI3, 
	//				int* maskSize);
	
    void compute_flux(const int* lo, const int* hi,
                      const int* domlo, const int* domhi,
                      amrex_real* phi, const int* philo, const int* phihi,
                      amrex_real* fluxx, const int* fxlo, const int* fxhi,
                      amrex_real* fluxy, const int* fylo, const int* fyhi,
#if (BL_SPACEDIM == 3)   
                      amrex_real* fluxz, const int* fzlo, const int* fzhi,
#endif
					  amrex_real* vel_u, const int* ulo, const int* uhi,
                      amrex_real* vel_v, const int* vlo, const int* vhi,
#if (BL_SPACEDIM == 3)   
                      amrex_real* vel_w, const int* wlo, const int* whi,
#endif
                      const amrex_real* dx);

    void RK3_Williamson(const int* lo, const int* hi, //bx
                    amrex_real* phi_old, const int* phi_oldlo, const int* phi_oldhi,
                    amrex_real* phi, const int* philo, const int* phihi,
                    int* mask, const int* masklo, const int* maskhi,
                    amrex_real* vel_u, const int* ulo, const int* uhi,
                    amrex_real* vel_v, const int* vlo, const int* vhi,
                    amrex_real* vel_w, const int* wlo, const int* whi,
                    const amrex_real* dx, const amrex_real* dt, const amrex_real* time, 
                    int* mi_i, int* mi_j, int* mi_k, 
                    const int* maskSize,
                    int* iter, amrex_real* G, const int* Glo, const int* Ghi);
    
    void attractToIntf(void* particles_data, const int* np,
                  int** c_vectors, int* sizes, const int* clo, const int* chi,
                  amrex_real* phi, const int* philo, const int* phihi,
                  amrex_real* grid_normals, const int* nlo, const int* nhi,
                  amrex_real* vel_u, const int* ulo, const int* uhi,
                  amrex_real* vel_v, const int* vlo, const int* vhi,
                  amrex_real* vel_w, const int* wlo, const int* whi,
                  const int* maskSize, int* mi_i, int* mi_j, int* mi_k, 
                  const amrex_real* dx, const amrex_real* p_lo,
                  const int* polyOrder, int* numrows, int* numcols, amrex_real* w, amrex_real* P, amrex_real* A_inv);
    
    void RK3_Williamson_Ptc(void* particles_data,
                    amrex_real* vel_u, const int* ulo, const int* uhi,
                    amrex_real* vel_v, const int* vlo, const int* vhi,
                    amrex_real* vel_w, const int* wlo, const int* whi,
                    const amrex_real* dt, const amrex_real* time, 
                    const int* np,
                    int* pi_i, int* pi_j, int* pi_k, 
                    const amrex_real* dx, const amrex_real* p_lo, 
                    int* iter);//, amrex_real* G_Ptc, const int* Glo, const int* Ghi);
                    
    void Redistance(const int* lo, const int* hi, //bx
                    amrex_real* phi, const int* philo, const int* phihi,
                    int* mask, const int* masklo, const int* maskhi,
                    const amrex_real* dx, const amrex_real* dtau, amrex_real* tau,
                    int* mi_i, int* mi_j, int* mi_k, 
                    const int* maskSize);
                   
    void GridNormals(amrex_real* grid_normals, const int* nlo, const int* nhi, 
                      amrex_real* phi, const int* philo, const int* phihi, 
                      int* mi_i, int* mi_j, int* mi_k, const int* np, 
                      const amrex_real* dx);
    
    //finite diff
    void update_phi(const int* lo, const int* hi,
                    amrex_real* pold, const int* polo, const int* pohi,
                    amrex_real* pnew, const int* pnlo, const int* pnhi,
                    amrex_real* fluxx, const int* fxlo, const int* fxhi,
                    amrex_real* fluxy, const int* fylo, const int* fyhi,
#if (BL_SPACEDIM == 3)   
                    amrex_real* fluxz, const int* fzlo, const int* fzhi,
#endif
                    const amrex_real* dx, const amrex_real* dt);
                    
    void polyMatrixDim(int* numrows, int* numcols, const int* polyOrder);
    void init_Pmatrix(const int* polyOrder, int* numrows, int* numcols, 
                      const amrex_real* dx, amrex_real* w_poly, amrex_real* P);
    void init_LeastSquaresMatrix(int* numrows, int* numcols, amrex_real* w, amrex_real* P, amrex_real* A_inv);
    void C_getPolyInterpCoeffs(const int* polyOrder, int* numrows, int* numcols, 
                             const amrex_real* dx, amrex_real* w, amrex_real* P, amrex_real* A_inv,
                             amrex_real* phi, const int* philo, const int* phihi,
                             int* indexVec, amrex_real* a_coeffs);
                             
    void minimizeError(void* particles_data, const int* np,
                       const int* lo, const int* hi,
                       amrex_real* phi, const int* philo, const int* phihi, 
                       int** c_vectors, int* sizes, const int* clo, const int* chi,
                       int* mi_i, int* mi_j, int* mi_k, const int* maskSize,
                       const amrex_real* p_lo, const amrex_real* dx,
                       const int* polyOrder, int* numrows, int* numcols, amrex_real* w, amrex_real* P, amrex_real* A_inv); 
                       
    void crossingTime(const int* lo, const int* hi, const amrex_real* dx,
                        amrex_real* phi, const int* philo, const int* phihi, 
                        amrex_real* phi_old, const int* phi_oldlo, const int* phi_oldhi,
                        amrex_real* dist, const int* distlo, const int* disthi,
                        int* mask, const int* masklo, const int* maskhi,
                        amrex_real* grid_normals, const int* nlo, const int* nhi,
                        int* mi_i, int* mi_j, int* mi_k, const int* maskSize,
                        amrex_real* tau, const amrex_real* dtau,
                        int* iter_redist, const int* max_iter);


#ifdef __cplusplus
}
#endif

#endif
