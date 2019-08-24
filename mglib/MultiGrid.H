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


#ifndef _MULTIGRID_H_
#define _MULTIGRID_H_

//
// $Id: MultiGrid.H,v 1.11 2001/08/01 21:51:04 lijewski Exp $
//

#include <Array.H>
#include <Pointers.H>
#include <MultiFab.H>
#include <ABecLaplacian.H>
#include <CGSolver.H>

#include <algorithm>

class MultiGrid
{
public:
    //
    //@ManDoc: constructor
    //
    MultiGrid (ABecLaplacian& _Lp);
    //
    //@ManDoc: destructor
    //
    virtual ~MultiGrid ();
    //
  void solve (int nsverbose,
    MultiFab& _sol, MultiFab& _rhs, 
    Real _eps_abs,Real _atol_b,
    int usecg_at_bottom,MultiFab& pbdry,
    Array<int> bcpres_array,
    int smooth_type,
    int bottom_smooth_type,int presmooth,
    int postsmooth);
    //
    //@ManDoc: return the linear operator
    //
    ABecLaplacian& linOp ();
    //
    //@ManDoc: set the maximum permitted multigrid iterations
    //
    void setMaxIter (int _maxiter);
    //
    //@ManDoc: return the maximum permitted multigrid iterations
    //
    int getMaxIter () const;
    //
    //@ManDoc: set the number of multigrid iterations to perform
    //
    void setNumIter (int _numiter);
    //
    //@ManDoc: return the number of multigrid iterations
    //
    int getNumIter () const;
    //
    //
    //@ManDoc: set/return the number of multigrid levels
    //
    int getNumLevels (int _numlevels);
    //
    //@ManDoc: return the number of multigrid levels
    //
    int getNumLevels () const;
    //
    //@ManDoc: set the verbosity value
    //
    void setVerbose (int _verbose);
    //
    //@ManDoc: return the verbosity value
    //
    int getVerbose ();
    //@ManDoc: Return the number of level relaxations (not implemented)
    //
    int cntRelax () const;
    //
    //@ManDoc: Return the number of smoothing passes at bottom of MG (if no cg)
    //
    int finalSmooth () const;
    //
    //
    //@ManDoc: set the number of post-cg relax passes
    //
    void set_nu_b (int _nu_b);
    //
    //@ManDoc: set the number of post-cg relax passes
    //
    int get_nu_b () const;

    void set_maxiter_b (int n);
    //
    //@ManDoc: get the maximum permitted relative tolerance
    //
    int  get_maxiter_b () const;
protected:
    //@ManMemo: Internal members
    
    //
    //@ManDoc: Solve the linear system to relative and absolute tolerance
    //
    int solve_ (int nsverbose,MultiFab& _sol,
     Real eps_abs,Real atol_b,MultiFab& pbdry,Array<int> bcpres_array,
     int usecg_at_bottom,
     int smooth_type,int bottom_smooth_type,int presmooth,int postsmooth);

    //
    //@ManDoc: Put the system in r-c form
    //
    void residualCorrectionForm (MultiFab& newrhs,
      MultiFab& oldrhs,MultiFab& solnL,
      MultiFab& inisol,MultiFab& pbdry,Array<int> bcpres_array,
      int level);
    //
    //@ManDoc: Make space, set switches for new solution level
    //
    void prepareForLevel (int level);
    //
    //@ManDoc: Compute the number of multigrid levels, assuming ratio=2
    //
    int numLevels () const;
    //
    //@ManDoc: Return scalar estimate of error
    //
    virtual Real errorEstimate (int level,
        MultiFab& pbdry,Array<int> bcpres_array);
    //
    //@ManDoc: Transfer MultiFab from fine to coarse level
    //
    void average (MultiFab& c,MultiFab& f,int clevel,int flevel);
    //
    //@ManDoc: Transfer MultiFab from coarse to fine level
    //
    void interpolate (MultiFab& f,MultiFab& c,int clevel,int flevel);

    void coarsestSmooth(MultiFab& solL,MultiFab& rhsL,
      int level,Real eps_abs,Real atol_b,
      int usecg_at_bottom,
      MultiFab& pbdry,Array<int> bcpres_array,
      int smooth_type,int bottom_smooth_type,
      int presmooth,int postsmooth);

    void relax (MultiFab& solL,MultiFab& rhsL,
      int level,Real eps_abs,
      Real atol_b,int usecg_at_bottom,
      MultiFab& pbdry,Array<int> bcpres_array,
      int smooth_type,
      int bottom_smooth_type,int presmooth,int postsmooth);

    //
    //@ManDoc: set flags, etc
    //
    static void initialize ();
    //
    //@ManDoc: default number of level, pre-, post- and bottom smooths
    //
    static int def_nu_0, def_nu_f;
    //
    //@ManDoc: default number of post-cg relax passes
    //
    static int def_nu_b;
    //
    //@ManDoc: default maximum number of complete MG cycle iterations
    //
    static int def_maxiter;
    //
    //@ManDoc: default number of complete MG cycle iterations to perform
    //
    static int def_numiter;
    //
    //@ManDoc: default verbosity
    //
    static int def_verbose;
    //
    //@ManDoc: default relative, absolute tolerance for cg solve
    //
    static int def_maxiter_b;
    //
    //@ManDoc: default maximum number of mg levels
    //
    static int def_numLevelsMAX;
    static int def_use_bicgstab_in_mglib;
    //
    //@ManDoc: verbosity
    //
    int verbose;
    //
    //@ManDoc: Number of MG levels
    //
    int numlevels;
    //
    //@ManDoc: current maximum number of allowed iterations
    //
    int maxiter;
    //
    //@ManDoc: current number of iterations to take
    //
    int numiter;
    //
    //@ManDoc: current number of level, pre-, post- and bottom smooths
    //
    int nu_0, nu_f;
    //
    //@ManDoc: current number of post-cg relax passes
    //
    int nu_b;
    //
    //
    //@ManDoc: relative, absolute tolerance, maxiter for cg bottom solver.
    //
    int maxiter_b;
    //
    //@ManDoc: maximum number of mg levels
    //
    int numLevelsMAX;
    //
    //@ManDoc: internal temp data to store initial guess of solution
    //
    MultiFab* initialsolution;
    //
    //@ManDoc: internal temp data
    //
    Array< MultiFab* > res;
    //
    //@ManDoc: internal temp data
    //
    Array< MultiFab* > rhs;
    //
    //@ManDoc: internal temp data
    //
    Array< MultiFab* > cor;
    //
    //@ManDoc: internal reference to linear operator
    //
    ABecLaplacian &Lp;

private:
    //
    // Flag, whether initialized.
    //
    static bool initialized;
    //
    // Disallow copy constructor, assignment operator
    //
    MultiGrid (const MultiGrid&);
    MultiGrid& operator= (const MultiGrid&);
};

inline
ABecLaplacian&
MultiGrid::linOp ()
{
    return Lp;
}

inline
void
MultiGrid::setMaxIter (int _maxiter)
{
    maxiter = _maxiter;
}

inline
int
MultiGrid::getMaxIter () const
{
    return maxiter;
}

inline
void
MultiGrid::setNumIter (int _numiter)
{
    numiter = _numiter;
}

inline
int
MultiGrid::getNumIter () const
{
    return numiter;
}

inline
int
MultiGrid::getNumLevels () const
{
    return numlevels;
}

inline
void
MultiGrid::setVerbose (int _verbose)
{
    verbose = _verbose;
}

inline
int
MultiGrid::getVerbose ()
{
    return verbose;
}

inline
int
MultiGrid::cntRelax () const
{
    return nu_0;
}

inline
int
MultiGrid::finalSmooth () const
{
    return nu_f;
}

inline
void
MultiGrid::set_maxiter_b (int n)
{
    maxiter_b = n;
}

inline
int
MultiGrid::get_maxiter_b () const
{
  return maxiter_b;
}

inline
void
MultiGrid::set_nu_b (int _nu_b)
{
    nu_b = _nu_b;
}

inline
int
MultiGrid::get_nu_b () const
{
    return nu_b;
}

inline
int
MultiGrid::getNumLevels (int _numlevels)
{
    BL_ASSERT(_numlevels >= 0);
    int oldnumlevels = numlevels;
    numlevels = std::min(_numlevels, numLevels());
    return oldnumlevels;
}

#endif /*_MULTIGRID_H_*/
