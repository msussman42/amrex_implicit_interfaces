#ifndef _ABecLAPLACIAN_H_
#define _ABecLAPLACIAN_H_

#include <local_thread_class.H>
#include <AMReX_Geometry.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include <AMReX_BoxArray.H>
#include <AMReX_MultiFab.H>

namespace amrex{

class ABecLaplacian
{
public:
    int cfd_level;
    int cfd_project_option;
    Real cfd_mglib_min_coeff_factor;
    bool cfd_tiling;

      // we use cfd_max_grid_size[0] in order to determine dmap
      // on the coarser mg levels.
    Vector<int> cfd_max_grid_size;

    int laplacian_solvability; // 0=> do not project, 1=> project RHS

    ABecLaplacian (
     const Vector<int>& ns_max_grid_size,
     const BoxArray& grids,
     const Geometry& geom,
     const DistributionMapping& dmap,
     int bfact,
     int cfd_level_in,
     int cfd_project_option_in,
     Real cfd_mglib_min_coeff_factor_in,
     int nsolve,
     bool ns_tiling_in,
     int _use_mg_precond_at_top);

    virtual ~ABecLaplacian ();

    void aCoefficients (MultiFab& _a);
    void generateCoefficients();
    void bCoefficients (MultiFab& _b, int dir);

    void apply (MultiFab& out, MultiFab& in,
     int level,MultiFab& pbdry,Vector<int> bcpres_array);
    void applyBC (MultiFab& inout,int level,
       MultiFab& pbdry,Vector<int> bcpres_array);
    void residual (MultiFab& residL,MultiFab& rhsL,
       MultiFab& solnL,int level,
       MultiFab& pbdry,Vector<int> bcpres_array);
    void smooth (MultiFab& solnL,MultiFab& rhsL,
      int level,MultiFab& pbdry,Vector<int> bcpres_array,
      int smooth_type);

    void project_null_space(MultiFab& rhsL,int level);

    void Fsmooth (MultiFab& solnL,
                  MultiFab& rhsL,
                  int level,
                  int smooth_type);

    void pcg_solve(MultiFab* z,MultiFab* r,
     Real eps_abs,Real bot_atol,
     MultiFab* pbdryhom,
     Vector<int> bcpres_array,
     int usecg_at_bottom,
     int smooth_type,int bottom_smooth_type,
     int presmooth,int postsmooth,
     int use_PCG,
     int level, 
     int caller_id);

    void CG_solve(
      int& cg_cycles_out,
      int nsverbose,int is_bottom,
      MultiFab& solnL,
      MultiFab& rhsL,
      Real eps_abs,Real bot_atol,
      MultiFab& pbdry,Vector<int> bcpres_array,
      int usecg_at_bottom,
      int& meets_tol,int smooth_type,
      int bottom_smooth_type,int presmooth,
      int postsmooth,
      Real& error0,
      int level);

    void CG_dump_params(
      Real rnorm,
      Real rnorm0,
      Real eps_abs,Real relative_error,
      int is_bottom,Real bot_atol,
      int usecg_at_bottom,int smooth_type,
      int bottom_smooth_type,int presmooth,
      int postsmooth,MultiFab& mf1,
      MultiFab& mf2,
      int level);


    void MG_solve (int nsverbose,
     MultiFab& _sol, MultiFab& _rhs, 
     Real _eps_abs,Real _atol_b,
     int usecg_at_bottom,MultiFab& pbdry,
     Vector<int> bcpres_array,
     int smooth_type,
     int bottom_smooth_type,int presmooth,
     int postsmooth);

    Real LPnorm(MultiFab &in, int level = 0) const;

    const Geometry& getGeom (int level) const;
    const BoxArray& LPboxArray (int level) const;
    const DistributionMapping& dmapLevel (int level) const;
    int get_bfact_array (int level) const;

    void LP_dot(const MultiFab& w_in,
		const MultiFab& p_in,
		int level_in,
		Real& dot_result);

     // sol=y+alpha p
    void LP_update(MultiFab& sol,Real alpha,
     MultiFab& y,const MultiFab& p,int level);

    static int CG_def_maxiter;
    static int CG_def_restart_period;

protected:

   void Bottom_Krylov_checkpoint(
     int vcycle,
     Real krylov_error,
     Real& best_error,
     int& best_iter,
     MultiFab& sol,
     MultiFab& restart_sol,
     int& restart_flag);

    void CG_advance (
       MultiFab& p,
       Real beta,
       const MultiFab& z,
       MultiFab& y,
       int level);

    void CG_check_for_convergence(
      int coarsefine,
      int presmooth,int postsmooth,
      Real rnorm,
      Real rnorm0,
      Real eps_abs,
      Real relative_error,
      int nit,
      int& error_close_to_zero,
      int level);

    static int CG_MG_defaults_reported;

    static int CG_def_verbose;

    //
    // The data.
    //
    int        CG_maxiter;  // Current maximum number of allowed iterations.
    int        CG_restart_period;
    int        CG_verbose;  // Current verbosity level.
    int CG_use_mg_precond_at_top; // Use multigrid as a preconditioner.

    Vector<Real> pw_dotprod_var;

     //2*coarsefine   = rnorm
     //2*coarsefine+1 = eps_abs
    Vector< Array<Real,4> > CG_error_history;

    void MG_solve_ (int nsverbose,MultiFab& _sol,
     Real eps_abs,Real atol_b,MultiFab& pbdry,Vector<int> bcpres_array,
     int usecg_at_bottom,
     int smooth_type,int bottom_smooth_type,int presmooth,int postsmooth);

    void MG_residualCorrectionForm (MultiFab& newrhs,
      MultiFab& oldrhs,MultiFab& solnL,
      MultiFab& inisol,MultiFab& pbdry,Vector<int> bcpres_array,
      int level);

    int MG_numLevels (const Box& local_box) const;
    virtual Real MG_errorEstimate (int level,
        MultiFab& pbdry,Vector<int> bcpres_array);
    void MG_average (MultiFab& c,MultiFab& f,int clevel,int flevel);
    void MG_interpolate (MultiFab& f,MultiFab& c,int clevel,int flevel);
    void MG_coarsestSmooth(MultiFab& solL,MultiFab& rhsL,
      int level,Real eps_abs,Real atol_b,
      int usecg_at_bottom,
      MultiFab& pbdry,Vector<int> bcpres_array,
      int smooth_type,int bottom_smooth_type,
      int presmooth,int postsmooth);

    void MG_relax (MultiFab& solL,MultiFab& rhsL,
      int level,Real eps_abs,
      Real atol_b,int usecg_at_bottom,
      MultiFab& pbdry,Vector<int> bcpres_array,
      int smooth_type,
      int bottom_smooth_type,int presmooth,int postsmooth);

    static int MG_def_nu_0;
    static int MG_def_nu_f;
    static int MG_def_nu_b;
    static int MG_def_verbose;

    int MG_verbose;
    int MG_numlevels_var;
    int CG_numlevels_var;
    int MG_nu_0;
    int MG_nu_f;
    int MG_nu_b;

    MultiFab* MG_initialsolution;
    Vector< MultiFab* > MG_res;
    Vector< MultiFab* > MG_rhs;
    Vector< MultiFab* > MG_cor;
    Vector< MultiFab* > MG_pbdrycoarser;


 void makeCoefficients (
	MultiFab& crs,
        const MultiFab& fine,
        int             level,
        int             avg);
		 
 void buildMatrix();

 virtual void Fapply (MultiFab& out,
                      MultiFab& in,
                      int level);
 virtual void Fdiagsum(MultiFab&       out,
                       int             level);

 Vector< MultiFab* > acoefs;

// bxleft,byleft,bzleft,bxright,byright,bzright,icbx,icby,icbz,diag,diagsing
// mask,icdiag,icdiagrb,ax,solnsave,rhssave,redsoln,blacksoln
 Vector< MultiFab* > workcoefs;

 Vector< MultiFab* > laplacian_ones;

 Vector< Array< MultiFab*, AMREX_SPACEDIM> > bcoefs;

 Array< MultiFab*, 2> CG_delta_sol;
 Array< MultiFab*, 2> CG_r;
 Array< MultiFab*, 2> CG_z;
 Array< MultiFab*, 2> CG_Av_search;
 Array< MultiFab*, 2> CG_p_search;
 Array< MultiFab*, 2> CG_p_search_SOLN;
 Array< MultiFab*, 2> CG_v_search;
 Array< MultiFab*, 2> CG_rhs_resid_cor_form;
 Array< MultiFab*, 2> CG_pbdryhom;

  // MG_pdryhom only needed by multigrid on the finest level
 MultiFab* MG_pbdryhom;
 Vector< MultiFab* > MG_CG_diagsumL;
 Vector< MultiFab* > MG_CG_ones_mf_copy;

 Vector< BoxArray > gbox;
 Vector< DistributionMapping > dmap_array;
 Vector< int > bfact_array;
 Vector<MultiFab*> maskvals;
 Vector< Geometry > geomarray;

 int nsolve_ABec;

private:

    static int mglib_blocking_factor;

    static int nghostRHS;
    static int nghostSOLN;

    static Real a_def;
    static Real b_def;

};


inline
void
ABecLaplacian::aCoefficients (MultiFab& _a)
{
    if (! _a.ok())
     amrex::Error("_a not ok");

    if (_a.boxArray() != (acoefs[0])->boxArray())
     amrex::Error("_a and acoefs boxarrays do not match");
    if (_a.nComp()!=nsolve_ABec) 
     amrex::Error("_a has wrong ncomp");

     // src,scomp,dcomp,snghost,dnghost
    acoefs[0]->ParallelCopy(_a,0,0,nsolve_ABec,0,0);

} // subroutine aCoefficients

// 1. generate coefficients on the finest level.
// 2. average down coefficients and generate on the coarser levels.
inline
void
ABecLaplacian::generateCoefficients () {
 buildMatrix();
}


inline
void
ABecLaplacian::bCoefficients (MultiFab& _b,
                              int dir) {

    if (! _b.ok())
     amrex::Error("_b not ok");

    if (_b.boxArray() != (bcoefs[0][dir])->boxArray())
     amrex::Error("_b box arrays do not match");

    if (_b.nComp()!=nsolve_ABec) 
     amrex::Error("_b has wrong ncomp");

     // src,scomp,dcomp,snghost,dnghost
    bcoefs[0][dir]->ParallelCopy(_b,0,0,nsolve_ABec,0,0);

} // subroutine bCoefficients


inline
const BoxArray&
ABecLaplacian::LPboxArray (int level) const
{
 if ((level>=0)&&(level<MG_numlevels_var)) {
  // do nothing
 } else {
  amrex::Error("level invalid");
 }
 return gbox[level];
}

inline
const DistributionMapping&
ABecLaplacian::dmapLevel (int level) const
{
 if ((level>=0)&&(level<MG_numlevels_var)) {
  // do nothing
 } else
  amrex::Error("level invalid");
 return dmap_array[level];
}

inline int
ABecLaplacian::get_bfact_array (int level) const
{   
 if ((level>=0)&&(level<MG_numlevels_var)) {
  // do nothing 
 } else
  amrex::Error("level invalid");
 return bfact_array[level];
}   

inline
const Geometry&
ABecLaplacian::getGeom (int level) const
{
 if ((level>=0)&&(level<MG_numlevels_var)) {
  // do nothing
 } else
  amrex::Error("level invalid");
 return geomarray[level];
}


}/* namespace amrex */

#endif /*_ABecLAPLACIAN_H_*/
