/*
** (c) 1996-2000 The Regents of the University of California (through
** E.O. Lawrence Berkeley National Laboratory), subject to approval by
** the U.S. Department of Energy.  Your use of this software is under
** license -- the license agreement is attached and included in the
** directory as license.txt or you may contact Berkeley Lab's Technology
** Transfer Department at TTD@lbl.gov.  NOTICE OF U.S. GOVERNMENT RIGHTS.
** The Software was developed under funding from the U.S. Government
** which consequently retains certain rights as follows: the
** U.S. Government has been granted for itself and others acting on its
** behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
** Software to reproduce, prepare derivative works, and perform publicly
** and display publicly.  Beginning five (5) years after the date
** permission to assert copyright is obtained from the U.S. Department of
** Energy, and subject to any subsequent five (5) year renewals, the
** U.S. Government is granted for itself and others acting on its behalf
** a paid-up, nonexclusive, irrevocable, worldwide license in the
** Software to reproduce, prepare derivative works, distribute copies to
** the public, perform publicly and display publicly, and to permit
** others to do so.
*/


#ifndef _CGSOLVER_H_
#define _CGSOLVER_H_

//
// $Id: CGSolver.H,v 1.10 2001/08/01 21:51:03 lijewski Exp $
//

#if defined(BL_OLD_STL)
#include <math.h>
#else
#include <cmath>
#endif

#include <Array.H>
#include <Pointers.H>
#include <MultiFab.H>
#include <ABecLaplacian.H>

class MultiGrid;

class CGSolver
{
public:
  //
    //@ManDoc: constructor
    //
    CGSolver (ABecLaplacian& _Lp,
              int    _use_mg_precond,
              int    _lev);
    //
    //@ManDoc: destructor
    //
    virtual ~CGSolver ();

    void pcg_solve(MultiFab* z,MultiFab* r,
     Real eps_abs,Real bot_atol,
     MultiFab* pbdryhom,
     Array<int> bcpres_array,
     int usecg_at_bottom,
     int smooth_type,int bottom_smooth_type,
     int presmooth,int postsmooth,
     int use_PCG,
     int ncomp,int nghost,int nghostRHS);

    void solve(
      int bicgstab_flag,
      int nsverbose,int is_bottom,
      MultiFab& solnL,MultiFab& rhsL,
      Real eps_abs,Real bot_atol,
      MultiFab& pbdry,Array<int> bcpres_array,
      int usecg_at_bottom,
      int& meets_tol,int smooth_type,
      int bottom_smooth_type,int presmooth,
      int postsmooth,Real& error0);


    void dump_params(Real rnorm,Real rnorm0,
      Real eps_abs,Real relative_error,
      int is_bottom,Real bot_atol,
      int usecg_at_bottom,int smooth_type,
      int bottom_smooth_type,int presmooth,
      int postsmooth,MultiFab& mf1,
      MultiFab& mf2);

    //
    //@ManDoc: set maximum allowed number of CG iterations
    //
    void setMaxIter (int _maxiter);
    //
    //@ManDoc: get maximum allowed number of CG iterations
    //
    int getMaxIter () const;
    //
    //
    //@ManDoc: set flag determining whether MG preconditioning is used
    //
    void setUseMGPrecond (int _use_mg_precond);
    //
    //@ManDoc: get  flag determining whether MG preconditioning is used
    //
    int getUseMGPrecond () const;
    //
    //@ManDoc: set the verbosity value
    //
    void setVerbose (int _verbose);
    //
    //@ManDoc: return the verbosity value
    //
    int getVerbose ();

protected:

    //
    //@ManDoc: construct work space, initialize parameters
    //
    static void initialize ();
    //
    //@ManDoc: if  (use\_mg\_precond == 1) then define the MultiGrid * mg\_precond
    //
    virtual void set_mg_precond ();
    //
    //@ManDoc: compute p = z  +  beta p in the CG algorithm
    //
    void advance (MultiFab&       p,
                  Real            beta,
                  const MultiFab& z,
                  MultiFab& y);

    void check_for_convergence(Real rnorm,Real rnorm0,Real eps_abs,
           Real relative_error,int nit,int& error_close_to_zero);

    //
    //@ManDoc: Initialized yet?
    //
    static int initialized;
    //
    //@ManDoc: Default maximum number of allowed iterations.
    //
    static int def_maxiter;
    //
    //@ManDoc: default unstable_criterion used to test for loss of accuracy in
    //
    static double def_unstable_criterion;
    //
    //
    //@ManDoc: Default verbosity flag.
    //
    static int def_verbose;
    //
    // The data.
    //
    ABecLaplacian&     Lp;     // Operator for linear system to be solved.
    MultiGrid* mg_precond;     // MultiGrid solver to be used as preconditioner
    int        maxiter;        // Current maximum number of allowed iterations.
    int        verbose;        // Current verbosity level.
    int        lev;            // Level of the linear operator to use
    int        use_mg_precond; // Use multigrid as a preconditioner.
    
private:
    //
    // Disable copy constructor and assignment operator.
    //
    CGSolver (const CGSolver&);
    CGSolver& operator= (const CGSolver&);
};

inline
void
CGSolver::setMaxIter (int _maxiter)
{
    maxiter = _maxiter;
}

inline
int
CGSolver::getMaxIter () const
{
    return maxiter;
}

inline
void
CGSolver::setUseMGPrecond (int _use_mg_precond)
{
    use_mg_precond = _use_mg_precond;
    set_mg_precond();
}

inline
int
CGSolver::getUseMGPrecond () const
{
    return use_mg_precond;
}

inline
void
CGSolver::setVerbose (int _verbose)
{
    verbose = _verbose;
}

inline
int
CGSolver::getVerbose ()
{
    return verbose;
}

#endif /*_CGSOLVER_H_*/
