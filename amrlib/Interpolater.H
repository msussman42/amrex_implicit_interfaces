#ifndef AMREX_INTERPOLATER_H_
#define AMREX_INTERPOLATER_H_ 

#include <AMReX_Box.H>
#include <AMReX_BCRec.H>
#include <AMReX_REAL.H>

namespace amrex {

class Geometry;
class FArrayBox;
class InterpolaterBoxCoarsener;

//
// Virtual base class for interpolaters.
//
// Specifies interpolater interface for coarse-to-fine interpolation in space.
//

class Interpolater
{
public:
 virtual ~Interpolater ();
 virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf,int grid_type)=0;
 virtual void interp (
   Real time,
   const FArrayBox& crse,
   int              crse_comp,
   FArrayBox&       fine,
   int              fine_comp,
   int              ncomp,
   const Box&       fine_region,
   const Geometry&  crse_geom,
   const Geometry&  fine_geom,
   Vector<BCRec>&    bcr,
   int levelc,int levelf,
   int bfactc,int bfactf,
   int grid_type) = 0;

 virtual InterpolaterBoxCoarsener BoxCoarsener (int bfactc,int bfactf,
   int grid_type);

};

class PCInterp
    :
    public Interpolater
{
public:
    virtual ~PCInterp ();
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf,int grid_type);
    virtual void interp (
     Real time,
     const FArrayBox& crse,
     int              crse_comp,
     FArrayBox&       fine,
     int              fine_comp,
     int              ncomp,
     const Box&       fine_region,
     const Geometry&  crse_geom,
     const Geometry&  fine_geom,
     Vector<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf,
     int grid_type); // -1,0,1,2,3,4,5
};


class LSInterp
    :
    public Interpolater
{
public:
    int LSInterp_nmat;

    virtual ~LSInterp ();
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf,int grid_type);
    virtual void interp (
     Real time,
     const FArrayBox& crse,
     int              crse_comp,
     FArrayBox&       fine,
     int              fine_comp,
     int              ncomp,
     const Box&       fine_region,
     const Geometry&  crse_geom,
     const Geometry&  fine_geom,
     Vector<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf,
     int grid_type);
};


class SEMInterp
    :
    public Interpolater
{
public:

      // 0 - low order
      // 1 - space/time spectral
      // 2 - space spectral only
      // 3 - time spectral only
    int interp_enable_spectral;

    virtual ~SEMInterp ();
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf,int grid_type);
    virtual void interp (
     Real time,
     const FArrayBox& crse,
     int              crse_comp,
     FArrayBox&       fine,
     int              fine_comp,
     int              ncomp,
     const Box&       fine_region,
     const Geometry&  crse_geom,
     const Geometry&  fine_geom,
     Vector<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf,
     int grid_type); // -1,0,1,2,3,4,5
};


class maskSEMInterp
    :
    public Interpolater
{
public:
    virtual ~maskSEMInterp ();
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf,int grid_type);
    virtual void interp (
     Real time,
     const FArrayBox& crse,
     int              crse_comp,
     FArrayBox&       fine,
     int              fine_comp,
     int              ncomp,
     const Box&       fine_region,
     const Geometry&  crse_geom,
     const Geometry&  fine_geom,
     Vector<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf,
     int grid_type);
};


class PCInterpNull
    :
    public Interpolater
{
public:
    virtual ~PCInterpNull ();
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf,int grid_type);
    virtual void interp (
     Real time,
     const FArrayBox& crse,
     int              crse_comp,
     FArrayBox&       fine,
     int              fine_comp,
     int              ncomp,
     const Box&       fine_region,
     const Geometry&  crse_geom,
     const Geometry&  fine_geom,
     Vector<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf,
     int grid_type);
};



class RefineDensityInterp
    :
    public Interpolater
{
public:

    virtual ~RefineDensityInterp ();
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf,
       int grid_type);
    virtual void interp (
     Real time,
     const FArrayBox& crse,
     int              crse_comp,
     FArrayBox&       fine,
     int              fine_comp,
     int              ncomp,
     const Box&       fine_region,
     const Geometry&  crse_geom,
     const Geometry&  fine_geom,
     Vector<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf,
     int grid_type);
};



class RefineElasticInterp
    :
    public Interpolater
{
public:

    virtual ~RefineElasticInterp ();
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf,
       int grid_type);
    virtual void interp (
     Real time,
     const FArrayBox& crse,
     int              crse_comp,
     FArrayBox&       fine,
     int              fine_comp,
     int              ncomp,
     const Box&       fine_region,
     const Geometry&  crse_geom,
     const Geometry&  fine_geom,
     Vector<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf,
     int grid_type);
};



class multiMOFInterp
    :
    public Interpolater
{
public:
    int multiMOFInterp_nmat;
    int multiMOFInterp_ngeom_raw;
    int multiMOFInterp_ngeom_recon;

    virtual ~multiMOFInterp ();
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf,
       int grid_type);
    virtual void interp (
     Real time,
     const FArrayBox& crse,
     int              crse_comp,
     FArrayBox&       fine,
     int              fine_comp,
     int              ncomp,
     const Box&       fine_region,
     const Geometry&  crse_geom,
     const Geometry&  fine_geom,
     Vector<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf,
     int grid_type);
};


class multiEXTMOFInterp
    :
    public Interpolater
{
public:
    int multiMOFInterp_nmat;
    int multiMOFInterp_ngeom_raw;
    int multiMOFInterp_ngeom_recon;

    virtual ~multiEXTMOFInterp ();
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf,int grid_type);
    virtual void interp (
     Real time,
     const FArrayBox& crse,
     int              crse_comp,
     FArrayBox&       fine,
     int              fine_comp,
     int              ncomp,
     const Box&       fine_region,
     const Geometry&  crse_geom,
     const Geometry&  fine_geom,
     Vector<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf,
     int grid_type);
};

class BurnVelInterp
    :
    public Interpolater
{
public:
    int burnvel_nmat;
    int burnvel_nten;
    int burnvel_ncomp_per;
    int burnvel_ncomp;

    virtual ~BurnVelInterp ();
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf,int grid_type);
    virtual void interp (
     Real time,
     const FArrayBox& crse,
     int              crse_comp,
     FArrayBox&       fine,
     int              fine_comp,
     int              ncomp,
     const Box&       fine_region,
     const Geometry&  crse_geom,
     const Geometry&  fine_geom,
     Vector<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf,
     int grid_type);
};


class UMACInterp
    :
    public Interpolater
{
public:

      // 0 - low order
      // 1 - space/time spectral
      // 2 - space spectral only
      // 3 - time spectral only
    int interp_enable_spectral;

    virtual ~UMACInterp ();
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf,int grid_type);
    virtual void interp (
     Real time,
     const FArrayBox& crse,
     int              crse_comp,
     FArrayBox&       fine,
     int              fine_comp,
     int              ncomp,
     const Box&       fine_region,
     const Geometry&  crse_geom,
     const Geometry&  fine_geom,
     Vector<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf,
     int grid_type);
};

//
// CONSTRUCT A GLOBAL OBJECT OF EACH VERSION.
//
extern LSInterp                  ls_interp;
extern PCInterp                  pc_interp;
extern SEMInterp                 sem_interp_DEFAULT;
extern SEMInterp                 sem_interp_LOW_PARM;
extern SEMInterp                 sem_interp_HIGH_PARM;
extern PCInterpNull              pc_interp_null;
extern RefineDensityInterp       refine_density_pc_interp;
extern RefineElasticInterp       refine_elastic_pc_interp;
extern multiMOFInterp            multi_mof_interp;
extern multiEXTMOFInterp         multi_extmof_interp;
extern BurnVelInterp             burnvel_interp;
extern BurnVelInterp             tsat_interp;
extern BurnVelInterp             drag_interp;
extern UMACInterp                umac_interp;
extern maskSEMInterp             mask_sem_interp;

//class BoxConverter is declared in: AMReX_Box.H
// the members of BoxConverter are: "doit", "clone", and "~BoxConverter"
class InterpolaterBoxCoarsener
    : public BoxConverter
{
public:
	// constructor
    InterpolaterBoxCoarsener (Interpolater* mapper_,
      int bfactc_,int bfactf_,int grid_type_);

    virtual Box doit (const Box& fine) const;

    virtual BoxConverter* clone () const;

private:
    Interpolater* mapper;
    int bfactf;
    int bfactc;
    int grid_type;
};

} // namespace amrex

#endif /*AMREX_INTERPOLATER_H_*/
