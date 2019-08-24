#ifndef _ABecLAPLACIAN_H_
#define _ABecLAPLACIAN_H_

#include <Geometry.H>
#include <Array.H>
#include <Tuple.H>
#include <Pointers.H>
#include <REAL.H>
#include <BoxArray.H>
#include <MultiFab.H>

class ABecLaplacian
{
public:
    int cfd_level;
    int cfd_project_option;
    bool cfd_tiling;

    int laplacian_solvability; // 0=> do not project, 1=> project RHS
    int check_for_singular;  // 0=> diag<>0  1=> diag can be 0
    Real diag_regularization;

    ABecLaplacian (const BoxArray& grids,const Geometry& geom,
     int bfact,int cfd_level_in,int cfd_project_option_in,
     int nsolve,
     bool ns_tiling_in);
    ABecLaplacian (const ABecLaplacian& _Lp,int level);
    virtual ~ABecLaplacian ();

    const MultiFab& aCoefficients (int level);
    const MultiFab& bCoefficients (int dir,int level);
    const MultiFab& matCoefficients (int level);
     // =1 if diagonal <> 0
     // =0 otherwise
    const MultiFab& onesCoefficients (int level);
    void aCoefficients (MultiFab& _a);
    void generateCoefficients();
    void bCoefficients (MultiFab& _b, int dir);
    void non_sing_coefficients(Real offdiag_coeff_);

    virtual void prepareForLevel (int level);
    virtual void clearToLevel (int level);

    void invalidate_a_to_level (int lev);
    void invalidate_b_to_level (int lev);
    void invalidate_non_sing_to_level (int lev);

    void apply (MultiFab& out, MultiFab& in,
     int level,MultiFab& pbdry,Array<int> bcpres_array);
    void applyBC (MultiFab& inout,int level,
       MultiFab& pbdry,Array<int> bcpres_array);
    void residual (MultiFab& residL,MultiFab& rhsL,
       MultiFab& solnL,int level,
       MultiFab& pbdry,Array<int> bcpres_array);
    void smooth (MultiFab& solnL,MultiFab& rhsL,
      int level,MultiFab& pbdry,Array<int> bcpres_array,
      int smooth_type);

    void project_null_space(MultiFab& rhsL,int level);

    void Fsmooth (MultiFab& solnL,
                  MultiFab& rhsL,
                  int level,
                  int smooth_type);
    Real norm(MultiFab &in, int level = 0) const;
    const Geometry& getGeom (int level) const;
    const BoxArray& boxArray (int level) const;
    int get_bfact_array (int level) const;
    int get_nsolve () const;
    int numGrids () const;
    int numLevels () const;

    void LP_dot(MultiFab& w,const MultiFab& p,int level,Real& result);
     // sol=y+alpha p
    void LP_update(MultiFab& sol,Real alpha,
     MultiFab& y,const MultiFab& p,int level);

protected:

 void makeCoefficients (MultiFab&       crs,
                        const MultiFab& fine,
                        int             level,
                        int             avg);
 void buildMatrix(MultiFab& ones_mf,MultiFab& work,
  Real offdiag_coeff_level,
  MultiFab& a,MultiFab& bx,MultiFab& by,MultiFab& bz,
  int define_flag,int level);

 void initCoefficients (const BoxArray &_ba);
 virtual void Fapply (MultiFab& out,
                      MultiFab& in,
                      int level);
 virtual void Fdiagsum(MultiFab&       out,
                       int             level);

 Array< MultiFab* > acoefs;

// bxleft,byleft,bzleft,bxright,byright,bzright,icbx,icby,icbz,diag,diagsing
// mask,icdiag,icdiagrb,ax,solnsave,rhssave,redsoln,blacksoln
 Array< MultiFab* > workcoefs;
 Array< MultiFab* > laplacian_ones;

 Array< Real > offdiag_coeff;

 Array< Tuple< MultiFab*, BL_SPACEDIM> > bcoefs;

 std::vector< BoxArray > gbox;
 std::vector< int > bfact_array;
 Array<MultiFab*> maskvals;
 std::vector< Geometry > geomarray;

 int nsolve_bicgstab;

private:
    std::vector<bool> a_valid;
    std::vector<bool> b_valid;
    std::vector<bool> non_sing_valid;

    static Real a_def;
    static Real b_def;

    ABecLaplacian (const ABecLaplacian& rhs);
    ABecLaplacian& operator= (const ABecLaplacian& rhs);
};

inline
void
ABecLaplacian::aCoefficients (MultiFab& _a)
{
    if (! _a.ok())
     BoxLib::Error("_a not ok");

    if (_a.boxArray() != (acoefs[0])->boxArray())
     BoxLib::Error("_a and acoefs boxarrays do not match");
    if (_a.nComp()!=nsolve_bicgstab) 
     BoxLib::Error("_a has wrong ncomp");

    invalidate_a_to_level(0);
    ParallelDescriptor::Barrier();
    MultiFab::Copy(*acoefs[0],_a,0,0,nsolve_bicgstab,0);
    ParallelDescriptor::Barrier();
}

// generate coefficients on the finest level.
inline
void
ABecLaplacian::generateCoefficients ()
{
 int define_flag=0;
 int lev=0; // lev==0 is the finest level

 buildMatrix(*laplacian_ones[lev],*workcoefs[lev],
   offdiag_coeff[lev],
   *acoefs[lev],*bcoefs[lev][0],*bcoefs[lev][1],
   *bcoefs[lev][BL_SPACEDIM-1],define_flag,lev);
}


inline
void
ABecLaplacian::bCoefficients (MultiFab& _b,
                              int dir) {

    
    if (! _b.ok())
     BoxLib::Error("_b not ok");

    if (_b.boxArray() != (bcoefs[0][dir])->boxArray())
     BoxLib::Error("_b box arrays do not match");

    if (_b.nComp()!=nsolve_bicgstab) 
     BoxLib::Error("_b has wrong ncomp");

    invalidate_b_to_level(0);
    ParallelDescriptor::Barrier();
    MultiFab::Copy(*bcoefs[0][dir],_b,0,0,nsolve_bicgstab,0);
    ParallelDescriptor::Barrier();
}


inline
void
ABecLaplacian::non_sing_coefficients(Real offdiag_coeff_) {

 if (offdiag_coeff_>0.0) {
  // do nothing
 } else
  BoxLib::Error("offdiag_coeff_ not ok");

 invalidate_non_sing_to_level(0);
 offdiag_coeff[0]=offdiag_coeff_;

}  // end subroutine non_sing_coefficients


inline
const MultiFab&
ABecLaplacian::aCoefficients (int level)
{
    prepareForLevel(level);
    ParallelDescriptor::Barrier();
    return *acoefs[level];
}

inline
const MultiFab&
ABecLaplacian::bCoefficients (int dir,int level)
{
    prepareForLevel(level);
    ParallelDescriptor::Barrier();
    return *bcoefs[level][dir];
}


inline
const MultiFab&
ABecLaplacian::onesCoefficients (int level) {

    prepareForLevel(level);
    ParallelDescriptor::Barrier();
    return *laplacian_ones[level];
}



inline
const MultiFab&
ABecLaplacian::matCoefficients (int level) {

    prepareForLevel(level);
    ParallelDescriptor::Barrier();
    return *workcoefs[level];
}

inline
int
ABecLaplacian::numLevels () const
{
    return gbox.size();
}

inline
const BoxArray&
ABecLaplacian::boxArray (int level) const
{
    if (level>=numLevels())
     BoxLib::Error("level invalid");
    return gbox[level];
}

inline int
ABecLaplacian::get_nsolve() const
{
 return nsolve_bicgstab;
}

inline int
ABecLaplacian::get_bfact_array (int level) const
{   
    if (level>=numLevels())
     BoxLib::Error("level invalid");
    return bfact_array[level];
}   

inline
const Geometry&
ABecLaplacian::getGeom (int level) const
{
    if (level>=numLevels())
     BoxLib::Error("level invalid");
    return geomarray[level];
}


inline
int
ABecLaplacian::numGrids () const
{
    return gbox[0].size();
}

#endif /*_ABecLAPLACIAN_H_*/
