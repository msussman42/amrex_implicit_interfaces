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
// $Id: StateData.cpp,v 1.33 2001/08/09 22:42:00 marc Exp $
//
#include <winstd.H>

#include <algorithm>

#include <RealBox.H>
#include <StateData.H>
#include <StateDescriptor.H>
#include <ParallelDescriptor.H>

const Real INVALID_TIME = -1.0e30;
#define BL_IGNORE_MAX 100000

StateData::StateData () 
{
   desc = 0;
   new_data = old_data = 0;
   new_time.start = INVALID_TIME;
   new_time.stop  = INVALID_TIME;
   old_time.start = INVALID_TIME;
   old_time.stop  = INVALID_TIME;
}

StateData::StateData (const Box&             p_domain,
                      const BoxArray&        grds,
                      const StateDescriptor* d,
                      Real                   cur_time,
                      Real                   dt)
{
    define(p_domain, grds, *d, cur_time, dt);
}

void
StateData::define (const Box&             p_domain,
                   const BoxArray&        grds,
                   const StateDescriptor& d,
                   Real                   time,
                   Real                   dt)
{
    domain = p_domain;
    desc = &d;
    grids.define(grds);
    //
    // Convert to proper type.
    //
    IndexType typ(desc->getType());
    StateDescriptor::TimeCenter t_typ(desc->timeType());
    if (!typ.cellCentered())
    {
        domain.convert(typ);
        grids.convert(typ);
    }
    if (t_typ == StateDescriptor::Point) {
        new_time.start = new_time.stop = time;
        old_time.start = old_time.stop = time;
    } else
     BoxLib::Error("only Point allowed");

    int ncomp = desc->nComp();

    new_data = new MultiFab(grids,ncomp,desc->nExtra(),Fab_allocate);
    old_data = new MultiFab(grids,ncomp,desc->nExtra(),
     new_data->DistributionMap(),Fab_allocate);

    buildBC();
}

void
StateData::restart (std::istream&          is,
                    const StateDescriptor& d,
                    const std::string&     chkfile,
                    bool                   bReadSpecial)
{
    if (bReadSpecial)
        ParallelDescriptor::Abort();  // not implemented

    desc = &d;

    is >> domain;

    if (bReadSpecial)
    {
        BoxLib::readBoxArray(grids, is, bReadSpecial);
    }
    else
    {
        grids.readFrom(is);
    }

    is >> old_time.start;
    is >> old_time.stop;
    is >> new_time.start;
    is >> new_time.stop;

    int nsets;
    is >> nsets;
    if (nsets!=2)
     BoxLib::Error("both old and new data should be checkpointed");

    old_data = 0;

    new_data = new MultiFab;
    std::string mf_name;
    is >> mf_name;
    //
    // Note that mf_name is relative to the Header file.
    // We need to prepend the name of the chkfile directory.
    //
    std::string FullPathName = chkfile;
    if (!chkfile.empty() && chkfile[chkfile.length()-1] != '/')
        FullPathName += '/';
    FullPathName += mf_name;
    VisMF::Read(*new_data, FullPathName);

    if (nsets == 2)
    {
        old_data = new MultiFab;
        is >> mf_name;
        //
        // Note that mf_name is relative to the Header file.
        // We need to prepend the name of the chkfile directory.
        //
        FullPathName = chkfile;
        if (!chkfile.empty() && chkfile[chkfile.length()-1] != '/')
            FullPathName += '/';
        FullPathName += mf_name;
        VisMF::Read(*old_data, FullPathName);
    }

    buildBC();
}

void
StateData::buildBC ()
{
    int ncomp = desc->nComp();
    bc.resize(ncomp);
    for (int i = 0; i < ncomp; i++)
    {
        bc[i].resize(grids.size());
        for (int j = 0; j < grids.size(); j++)
        {
            BCRec bcr;
            BoxLib::setBC(grids[j],domain,desc->getBC(i),bcr);
            bc[i].set(j,bcr);
        }
    }
}

StateData::~StateData()
{
   desc = 0;
   delete new_data;
   delete old_data;
}

const StateDescriptor*
StateData::descriptor () const
{
    return desc;
}

const Box&
StateData::getDomain () const
{
    return domain;
}

const BoxArray&
StateData::boxArray () const
{
    return grids;
}

Real
StateData::curTime () const
{
    return 0.5*(new_time.start + new_time.stop);
}

Real
StateData::prevTime () const
{
    return 0.5*(old_time.start + old_time.stop);
}

MultiFab&
StateData::newData ()
{
    if (new_data==0)
     BoxLib::Error("new_data==0");
    BL_ASSERT(new_data != 0);
    return *new_data;
}

const MultiFab&
StateData::newData () const
{
    if (new_data==0)
     BoxLib::Error("new_data==0");
    BL_ASSERT(new_data != 0);
    return *new_data;
}

MultiFab&
StateData::oldData ()
{
    if (old_data==0)
     BoxLib::Error("old_data==0");
    BL_ASSERT(old_data != 0);
    return *old_data;
}

const MultiFab&
StateData::oldData () const
{
    if (old_data==0)
     BoxLib::Error("old_data==0");
    BL_ASSERT(old_data != 0);
    return *old_data;
}

FArrayBox&
StateData::newGrid (int i)
{
    BL_ASSERT(new_data != 0);
    return (*new_data)[i];
}

FArrayBox&
StateData::oldGrid (int i)
{
    BL_ASSERT(old_data != 0);
    return (*old_data)[i];
}

Array<BCRec>&
StateData::getBCs (int comp)
{
    return bc[comp];
}

const BCRec&
StateData::getBC (int comp, int i) const
{
    return bc[comp][i];
}

bool
StateData::hasData () const
{
  if ((old_data==0)&&(new_data==0))
   return false;
  else if ((old_data!=0)&&(new_data!=0))
   return true;
  else
   BoxLib::Error("old/new data must both be allocated or not");
}

void
StateData::setOldTimeLevel (Real time)
{
    if (desc->timeType() == StateDescriptor::Point) {
        old_time.start = old_time.stop = time;
    } else {
        BoxLib::Error("StateData::setOldTimeLevel called with Interval");
    }
}

void
StateData::setNewTimeLevel (Real time)
{
    if (desc->timeType() == StateDescriptor::Point) {
        new_time.start = new_time.stop = time;
    } else {
        BoxLib::Error("StateData::setNewTimeLevel called with Interval");
    }
}

void
StateData::setTimeLevel (Real time,Real dt)
{
    if (desc->timeType() == StateDescriptor::Point)
    {
        new_time.start = new_time.stop = time;
        old_time.start = old_time.stop = time - dt;
    }
    else
    {
     BoxLib::Error("this option not available");
    }
}

void
StateData::CopyNewToOld(Real dt)
{
    old_time = new_time;
    if (desc->timeType() == StateDescriptor::Point)
    {
        new_time.start += dt;
        new_time.stop  += dt;
    }
    else
    {
        BoxLib::Error("invalid timeType");
    }
    MultiFab & oldmulti = *old_data;
    MultiFab & newmulti = *new_data;
    int ncomp=newmulti.nComp();
    int ngrow=newmulti.nGrow();
    MultiFab::Copy(oldmulti,newmulti,0,0,ncomp,ngrow);
}

void
StateData::checkPoint (const std::string& name,
                       const std::string& fullpathname,
                       std::ostream&  os,
                       VisMF::How     how)
{
    static const std::string NewSuffix("_New_MF");
    static const std::string OldSuffix("_Old_MF");

    if ((old_data==0)||(new_data==0))
     BoxLib::Error("both old/new data must be allocated");

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // The relative name gets written to the Header file.
        //
        std::string mf_name_old = name; mf_name_old += OldSuffix;
        std::string mf_name_new = name; mf_name_new += NewSuffix;

        os << domain << '\n';

        grids.writeOn(os);

        os << old_time.start << '\n'
           << old_time.stop  << '\n'
           << new_time.start << '\n'
           << new_time.stop  << '\n';

        os << 2 << '\n' << mf_name_new << '\n' << mf_name_old << '\n';
    }
    BL_ASSERT(new_data);
    std::string mf_fullpath_new = fullpathname; mf_fullpath_new += NewSuffix;
    VisMF::Write(*new_data,mf_fullpath_new,how);

    BL_ASSERT(old_data);
    std::string mf_fullpath_old = fullpathname; mf_fullpath_old += OldSuffix;
    VisMF::Write(*old_data,mf_fullpath_old,how);
}

void
StateData::printTimeInterval (std::ostream &os) const
{
    os << '['
       << old_time.start
       << ' '
       << old_time.stop
       << "] ["
       << new_time.start
       << ' '
       << new_time.stop
       << ']'
       << '\n';
}

//
// The following is from the asci version of StateData.C
//

void
BoxLib::readBoxArray (BoxArray&     ba,
                      std::istream& is,
                      bool          bReadSpecial)
{
    if (bReadSpecial == false)
    {
        ba.readFrom(is);
    }
    else
    {
        BL_ASSERT(ba.size() == 0);
        int maxbox;
        unsigned long in_hash; // will be ignored
        is.ignore(BL_IGNORE_MAX, '(') >> maxbox >> in_hash;
        ba.resize(maxbox);
        for (int i = 0; i < maxbox; i++)
        {
            Box b;
            is >> b;
            ba.set(i, b);
        }
        is.ignore(BL_IGNORE_MAX, ')');

        if (is.fail())
            BoxLib::Error("readBoxArray(BoxArray&,istream&,int) failed");
    }
}

