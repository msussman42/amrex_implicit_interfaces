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

//
// $Id: AmrLevel.cpp,v 1.81 2002/02/22 20:53:54 lijewski Exp $
//
#include <winstd.H>

#if defined(BL_OLD_STL)
#include <stdio.h>
#include <string.h>
#else
#include <cstdio>
#include <cstring>
#endif

#include <strstream>

#include <AmrLevel.H>
#include <Derive.H>
#include <BoxDomain.H>
#include <ParallelDescriptor.H>
#include <Utility.H>
#include <ParmParse.H>

DescriptorList AmrLevel::desc_lst;
DeriveList     AmrLevel::derive_lst;

void
AmrLevel::postCoarseTimeStep (Real time)
{}

int
AmrLevel::Level () const
{
    return level;
}

const BoxArray&
AmrLevel::boxArray () const
{
    return grids;
}

int
AmrLevel::numGrids () const {
    return grids.size();
}

int
AmrLevel::get_fine_ratio (int dir) const {
    return fine_ratio[dir];
}

const Array<RealBox>&
AmrLevel::gridLocations () const
{
    return grid_loc;
}

const Box&
AmrLevel::Domain () const
{
    return geom.Domain();
}

long
AmrLevel::nStep () const
{
    return parent->levelSteps(level);
}

const Geometry&
AmrLevel::Geom () const
{
    return geom;
}

StateData&
AmrLevel::get_state_data (int state_indx)
{
    return state[state_indx];
}

MultiFab&
AmrLevel::get_old_data (int state_indx)
{
    return state[state_indx].oldData();
}

const MultiFab&
AmrLevel::get_old_data (int state_indx) const
{
    return state[state_indx].oldData();
}

MultiFab&
AmrLevel::get_new_data (int state_indx)
{
    return state[state_indx].newData();
}

const MultiFab&
AmrLevel::get_new_data (int state_indx) const
{
    return state[state_indx].newData();
}

const DescriptorList&
AmrLevel::get_desc_lst ()
{
    return desc_lst;
}

DeriveList&
AmrLevel::get_derive_lst ()
{
    return derive_lst;
}

AmrLevel::AmrLevel ()
{
   parent = 0;
   level = -1;
}

AmrLevel::AmrLevel (Amr&            papa,
                    int             lev,
                    const Geometry& level_geom,
                    const BoxArray& ba,
                    Real            time)
    :
    geom(level_geom),
    grids(ba)
{
    level  = lev;
    parent = &papa;

    fine_ratio = IntVect::TheUnitVector(); fine_ratio.scale(-1);
    crse_ratio = IntVect::TheUnitVector(); crse_ratio.scale(-1);

    if (level > 0)
    {
        crse_ratio = parent->refRatio(level-1);
    }
    if (level < parent->maxLevel())
    {
        fine_ratio = parent->refRatio(level);
    }

    state.resize(desc_lst.size());

    Real current_dt=parent->dtLevel(lev);

    for (int i = 0; i < state.size(); i++) {
        state[i].define(geom.Domain(),
                        grids,
                        desc_lst[i],
                        time,
                        current_dt);
    }

    finishConstructor();
}

void
AmrLevel::restart (Amr&          papa,
                   std::istream& is,
		   bool          bReadSpecial)
{
    parent = &papa;

    is >> level;
    is >> geom;

    fine_ratio = IntVect::TheUnitVector(); fine_ratio.scale(-1);
    crse_ratio = IntVect::TheUnitVector(); crse_ratio.scale(-1);

    if (level > 0)
    {
        crse_ratio = parent->refRatio(level-1);
    }
    if (level < parent->maxLevel())
    {
        fine_ratio = parent->refRatio(level);
    }

    if (bReadSpecial)
    {
        BoxLib::readBoxArray(grids, is, bReadSpecial);
    }
    else
    {
        grids.readFrom(is);
    }

    int nstate;
    is >> nstate;
    int ndesc = desc_lst.size();
    if (nstate!=ndesc) {
     std::cout << "checkpoint nstate<>current ndesc: nstate= " << nstate <<
        "ndesc= " << ndesc << '\n';
     BoxLib::Error("aborting read checkpoint nstate<>ndesc");
    }

    state.resize(ndesc);
    for (int i = 0; i < ndesc; i++)
    {
        state[i].restart(is, desc_lst[i], papa.theRestartFile(), bReadSpecial);
    }

    finishConstructor();
}

void
AmrLevel::finishConstructor ()
{
    //
    // Set physical locations of grids.
    //
    grid_loc.resize(grids.size());

    for (int i = 0; i < grid_loc.size(); i++)
    {
        grid_loc[i] = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
    }
}

void
AmrLevel::setTimeLevel (Real time,Real dt)
{
    for (int k = 0; k < desc_lst.size(); k++)
    {
        state[k].setTimeLevel(time,dt);
    }
}

int
AmrLevel::isStateVariable (const std::string& name,
                           int&           typ,
                           int&            n)
{
    for (typ = 0; typ < desc_lst.size(); typ++)
    {
        const StateDescriptor& desc = desc_lst[typ];

        for (n = 0; n < desc.nComp(); n++)
        {
            if (desc.name(n) == name)
                return true;
        }
    }
    return false;
}

long
AmrLevel::countCells () const
{
    long cnt = 0;
    for (int i = 0; i < grids.size(); i++)
    {
        cnt += grids[i].numPts();
    }
    return cnt;
}

void
AmrLevel::checkPoint (const std::string& dir,
                      std::ostream&  os,
                      VisMF::How     how)
{
    int ndesc = desc_lst.size(), i;
    //
    // Build directory to hold the MultiFabs in the StateData at this level.
    // The directory is relative the the directory containing the Header file.
    //
    char buf[64];
    sprintf(buf, "Level_%d", level);
    std::string Level = buf;
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
    {
        FullPath += '/';
    }
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(FullPath, 0755))
            BoxLib::CreateDirectoryFailed(FullPath);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << '\n' << geom  << '\n';
        grids.writeOn(os);
        os << ndesc << '\n';
    }
    //
    // Output state data.
    //
    for (i = 0; i < ndesc; i++)
    {
        //
        // Now build the full relative pathname of the StateData.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        // There is only one MultiFab written out at each level in HyperCLaw.
        //
        std::string PathNameInHeader = Level;
        sprintf(buf, "/SD_%d", i);
        PathNameInHeader += buf;
        std::string FullPathName = FullPath;
        FullPathName += buf;
        state[i].checkPoint(PathNameInHeader, FullPathName, os, how);
    }
}

AmrLevel::~AmrLevel ()
{
    parent = 0;
}


int AmrLevel::derive_ncomp(const std::string& name) {

  int return_ncomp=0;

  int index, scomp, ncomp;

  if (isStateVariable(name, index, scomp)) {
   return_ncomp=1;
  } else if (const DeriveRec* rec = derive_lst.get(name)) {
   rec->getRange(0, index, scomp, ncomp);
   return_ncomp=rec->numDerive();
  } else {
   std::string msg("AmrLevel::derive(MultiFab*): unknown variable: ");
   msg += name;
   BoxLib::Error(msg.c_str());
  }

  return return_ncomp;
}

MultiFab*
AmrLevel::derive (const std::string& name,int is_old,int ngrow)
{
    BL_ASSERT(ngrow >= 0);

    MultiFab* mf = 0;

    int index, scomp, ncomp;

    if (isStateVariable(name, index, scomp)) {
     MultiFab& S_temp=state[index].newData();
     mf = new MultiFab(state[index].boxArray(), 1, ngrow,
               S_temp.DistributionMap(),Fab_allocate);
     int piecewise_constant=0;
     MultiFab* oldfine_data=LinInterpData(is_old,scomp,1);
     StableFillPatch(oldfine_data,*mf,ngrow,is_old,
       scomp,0,1,mf->boxArray());
     delete oldfine_data;
    } else if (const DeriveRec* rec = derive_lst.get(name)) {
        rec->getRange(0, index, scomp, ncomp);
        BoxArray srcBA(state[index].boxArray());
        BoxArray dstBA(state[index].boxArray());
        MultiFab& S_temp=state[index].newData();
        Real time=0.0;
        if (is_old==1)
         time=state[index].prevTime();
        else if (is_old==0)
         time=state[index].curTime();
        else
         BoxLib::Error("is_old invalid");

        srcBA.convert(rec->boxMap());
        dstBA.convert(rec->deriveType());

        if (srcBA!=dstBA)
         BoxLib::Error("srcBA<>dstBA");
        if (rec->deriveType()!=IndexType::TheCellType())
         BoxLib::Error("derive utility only available for cell data"); 

        MultiFab srcMF(srcBA, rec->numState(), ngrow+1,
                   S_temp.DistributionMap(),Fab_allocate);

        for (int k = 0, dc = 0; k < rec->numRange(); k++, dc += ncomp)
        {
            rec->getRange(k, index, scomp, ncomp);
            int piecewise_constant=0;
  
            MultiFab* oldfine_data=LinInterpData(is_old,scomp,ncomp);
            StableFillPatch(oldfine_data,srcMF,ngrow+1,is_old,scomp,dc,
             ncomp,srcMF.boxArray());
            delete oldfine_data;
        }

        mf = new MultiFab(dstBA, rec->numDerive(), ngrow,
                  S_temp.DistributionMap(),Fab_allocate);

        for (MFIter mfi(srcMF); mfi.isValid(); ++mfi)
        {
            int         grid_no = mfi.index();
            Real*       ddat    = (*mf)[grid_no].dataPtr();
            const int*  dlo     = (*mf)[grid_no].loVect();
            const int*  dhi     = (*mf)[grid_no].hiVect();
            int         n_der   = rec->numDerive();
            Real*       cdat    = srcMF[mfi].dataPtr();
            const int*  clo     = srcMF[mfi].loVect();
            const int*  chi     = srcMF[mfi].hiVect();
            int         n_state = rec->numState();
            const int*  dom_lo  = state[index].getDomain().loVect();
            const int*  dom_hi  = state[index].getDomain().hiVect();
            const Real* dx      = geom.CellSize();
            const int*  bcr     = rec->getBC();
            const Real* xlo     = grid_loc[grid_no].lo();
            Real        dt      = parent->dtLevel(level);

            Box gridbox = dstBA[grid_no];
            const int* lo=gridbox.loVect();
            const int* hi=gridbox.hiVect();

            rec->derFunc()(
             ddat,ARLIM(dlo),ARLIM(dhi),&n_der,
             cdat,ARLIM(clo),ARLIM(chi),&n_state,
             lo,hi,&ngrow,dom_lo,dom_hi,dx,xlo,&time,&dt,bcr,
             &level,&grid_no);
        }
        ParallelDescriptor::Barrier();
    }
    else
    {
        //
        // If we got here, cannot derive given name.
        //
        std::string msg("AmrLevel::derive(MultiFab*): unknown variable: ");
        msg += name;
        BoxLib::Error(msg.c_str());
    }

    return mf;
}

void
AmrLevel::derive(const std::string& name,int is_old, MultiFab& mf, int dcomp)
{
    BL_ASSERT(dcomp < mf.nComp());

    const int ngrow = mf.nGrow();

    int index, scomp, ncomp;

    if (isStateVariable(name,index,scomp)) {
     int piecewise_constant=0;

     MultiFab* oldfine_data=LinInterpData(is_old,scomp,1);
     StableFillPatch(oldfine_data,mf,ngrow,is_old,scomp,dcomp,
             1,mf.boxArray());
     delete oldfine_data;
    } else if (const DeriveRec* rec = derive_lst.get(name)) {
     rec->getRange(0,index,scomp,ncomp);
     BoxArray srcBA(mf.boxArray());
     BoxArray dstBA(mf.boxArray());

     srcBA.convert(state[index].boxArray()[0].ixType());
     BL_ASSERT(rec->deriveType() == dstBA[0].ixType());
     if (srcBA!=dstBA)
      BoxLib::Error("srcBA<>dstBA");
     if (rec->deriveType()!=IndexType::TheCellType())
      BoxLib::Error("derive utility only available for cell data"); 

     MultiFab srcMF(srcBA,rec->numState(),ngrow+1,
                mf.DistributionMap(),Fab_allocate);

     Real time=0.0;
     if (is_old==0)
      time=state[index].curTime();
     else if (is_old==1)
      time=state[index].prevTime();
     else
      BoxLib::Error("is_old invalid");

     for (int k = 0, dc = 0; k < rec->numRange(); k++, dc += ncomp)
     {
        rec->getRange(k,index,scomp,ncomp);
        int piecewise_constant=0;
  
        MultiFab* oldfine_data=LinInterpData(is_old,scomp,ncomp);
        StableFillPatch(oldfine_data,srcMF,ngrow+1,is_old,scomp,dc,
           ncomp,srcMF.boxArray());
        delete oldfine_data;
     }

     for (MFIter mfi(srcMF); mfi.isValid(); ++mfi) {
            int         grid_no = mfi.index();
            Real*       ddat    = mf[grid_no].dataPtr(dcomp);
            const int*  dlo     = mf[grid_no].loVect();
            const int*  dhi     = mf[grid_no].hiVect();
            int         n_der   = rec->numDerive();
            Real*       cdat    = srcMF[mfi].dataPtr();
            const int*  clo     = srcMF[mfi].loVect();
            const int*  chi     = srcMF[mfi].hiVect();
            int         n_state = rec->numState();
            const int*  dom_lo  = state[index].getDomain().loVect();
            const int*  dom_hi  = state[index].getDomain().hiVect();
            const Real* dx      = geom.CellSize();
            const int*  bcr     = rec->getBC();
            const RealBox temp  = RealBox(mf[grid_no].box(),geom.CellSize(),
                        geom.ProbLo());
            const Real* xlo     = temp.lo();
            Real        dt      = parent->dtLevel(level);

            Box gridbox = dstBA[grid_no];
            const int* lo=gridbox.loVect();
            const int* hi=gridbox.hiVect();
         
            rec->derFunc()(
             ddat,ARLIM(dlo),ARLIM(dhi),&n_der,
             cdat,ARLIM(clo),ARLIM(chi),&n_state,
             lo,hi,&ngrow,dom_lo,dom_hi,dx,xlo,&time,&dt,bcr,
             &level,&grid_no);
     }
     ParallelDescriptor::Barrier();
    }
    else
    {
        //
        // If we got here, cannot derive given name.
        //
        std::string msg("AmrLevel::derive(MultiFab*): unknown variable: ");
        msg += name;
        BoxLib::Error(msg.c_str());
    }
}

Array<int>
AmrLevel::getBCArray (int state_idx,
                      int gridno,
                      int strt_comp,
                      int ncomp)
{
    Array<int> bc(2*BL_SPACEDIM*ncomp);

    for (int n = 0; n < ncomp; n++)
    {
        const int* b_rec = state[state_idx].getBC(strt_comp+n,gridno).vect();
        for (int m = 0; m < 2*BL_SPACEDIM; m++)
            bc[2*BL_SPACEDIM*n + m] = b_rec[m];
    }

    return bc;
}

int
AmrLevel::okToRegrid ()
{
    return true;
}

void
AmrLevel::setPlotVariables ()
{
    ParmParse pp("amr");

    if (pp.contains("plot_vars"))
    {
        std::string nm;
      
        int nPltVars = pp.countval("plot_vars");
      
        for (int i = 0; i < nPltVars; i++)
        {
            pp.get("plot_vars", nm, i);

            if (nm == "ALL") 
                parent->fillStatePlotVarList();
            else if (nm == "NONE")
                parent->clearStatePlotVarList();
            else
                parent->addStatePlotVar(nm);
        }
    }
    else 
    {
        //
        // The default is to add them all.
        //
        parent->fillStatePlotVarList();
    }
  
    if (pp.contains("derive_plot_vars"))
    {
        std::string nm;
      
        int nDrvPltVars = pp.countval("derive_plot_vars");
      
        for (int i = 0; i < nDrvPltVars; i++)
        {
            pp.get("derive_plot_vars", nm, i);

            if (nm == "ALL") 
                parent->fillDerivePlotVarList();
            else if (nm == "NONE")
                parent->clearDerivePlotVarList();
            else
                parent->addDerivePlotVar(nm);
        }
    }
    else 
    {
        //
        // The default is to add none of them.
        //
        parent->clearDerivePlotVarList();
    }
}

