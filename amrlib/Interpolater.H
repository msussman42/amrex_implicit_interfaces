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
 virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf)=0;
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
   Array<BCRec>&    bcr,
   int levelc,int levelf,
   int bfactc,int bfactf) = 0;

 virtual InterpolaterBoxCoarsener BoxCoarsener (int bfactc,int bfactf);

};

class PCInterp
    :
    public Interpolater
{
public:
    virtual ~PCInterp ();
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf);
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
     Array<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf);
};


class LSHOInterp
    :
    public Interpolater
{
public:
    int LSHOInterp_nmat;
    int LSHOInterp_LO;  //0=>use normals  1=>do not use normals

    virtual ~LSHOInterp ();
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf);
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
     Array<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf);
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
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf);
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
     Array<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf);
};


class maskSEMInterp
    :
    public Interpolater
{
public:
    virtual ~maskSEMInterp ();
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf);
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
     Array<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf);
};


class PCInterpNull
    :
    public Interpolater
{
public:
    virtual ~PCInterpNull ();
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf);
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
     Array<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf);
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
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf);
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
     Array<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf);
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
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf);
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
     Array<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf);
};

class BurnVelInterp
    :
    public Interpolater
{
public:
    int burnvel_nmat;
    int burnvel_nten;

    virtual ~BurnVelInterp ();
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf);
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
     Array<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf);
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
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf);
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
     Array<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf);
};

class VMACInterp
    :
    public Interpolater
{
public:

      // 0 - low order
      // 1 - space/time spectral
      // 2 - space spectral only
      // 3 - time spectral only
    int interp_enable_spectral;

    virtual ~VMACInterp ();
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf);
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
     Array<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf);
};

class WMACInterp
    :
    public Interpolater
{
public:

      // 0 - low order
      // 1 - space/time spectral
      // 2 - space spectral only
      // 3 - time spectral only
    int interp_enable_spectral;

    virtual ~WMACInterp ();      
    virtual Box CoarseBox (const Box& fine,int bfactc,int bfactf);
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
     Array<BCRec>&    bcr,
     int levelc,int levelf,
     int bfactc,int bfactf);
};

//
// CONSTRUCT A GLOBAL OBJECT OF EACH VERSION.
//
extern LSHOInterp                ls_ho_interp_LOW_PARM;
extern LSHOInterp                ls_ho_interp_HIGH_PARM;
extern PCInterp                  pc_interp;
extern SEMInterp                 sem_interp_DEFAULT;
extern SEMInterp                 sem_interp_LOW_PARM;
extern SEMInterp                 sem_interp_HIGH_PARM;
extern PCInterpNull              pc_interp_null;
extern multiMOFInterp            multi_mof_interp;
extern multiEXTMOFInterp         multi_extmof_interp;
extern BurnVelInterp             burnvel_interp;
extern UMACInterp                umac_interp;
extern VMACInterp                vmac_interp;
extern WMACInterp                wmac_interp;
extern maskSEMInterp             mask_sem_interp;

class InterpolaterBoxCoarsener
    : public BoxConverter
{
public:
    InterpolaterBoxCoarsener (Interpolater* mapper_,int bfactc_,int bfactf_);
    virtual Box doit (const Box& fine) const;
    virtual BoxConverter* clone () const;
private:
    Interpolater* mapper;
    int bfactf;
    int bfactc;
};

} // namespace amrex

#endif /*AMREX_INTERPOLATER_H_*/
