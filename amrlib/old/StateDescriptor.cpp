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
// $Id: StateDescriptor.cpp,v 1.16 2001/08/09 16:20:37 marc Exp $
//
#include <winstd.H>

#include <algorithm>
#include <string>
#include <iostream>

#include <StateDescriptor.H>
#include <BCRec.H>


DescriptorList::DescriptorList ()
    :
    desc(PArrayManage)
{}

void
DescriptorList::clear ()
{
    desc.clear();
}

int
DescriptorList::size () const
{
    return desc.size();
}

void
DescriptorList::setComponent (int                               indx,
                              int                               comp,
                              const std::string&                nm,
                              const BCRec&                      bc,
                              int interp_type)
{
    desc[indx].setComponent(comp,nm,bc,interp_type);
}

const StateDescriptor&
DescriptorList::operator[] (int k) const
{
    return desc[k];
}

void
DescriptorList::addDescriptor (int                         indx,
                               IndexType                   typ,
                               StateDescriptor::TimeCenter ttyp,
                               int                         nextra,
                               int                         num_comp)
{
    if (indx >= desc.size())
        desc.resize(indx+1);
    desc.set(indx, new StateDescriptor(typ,ttyp,indx,nextra,num_comp));
}  

StateDescriptor::StateDescriptor ()
    :
    id(-1),
    ncomp(0),
    ngrow(0),
    t_type(Point)
{}

StateDescriptor::StateDescriptor (IndexType                   btyp,
                                  StateDescriptor::TimeCenter ttyp,
                                  int                         ident,
                                  int                         nextra, 
                                  int                         num_comp)
    :
    type(btyp),
    t_type(ttyp),
    id(ident),
    ngrow(nextra),
    ncomp(num_comp)
{
    BL_ASSERT (num_comp > 0);
   
    names.resize(num_comp);
    bc.resize(num_comp);
    interp_type_array.resize(num_comp);
}

StateDescriptor::~StateDescriptor ()
{
}

IndexType
StateDescriptor::getType () const
{
    return type;
}

StateDescriptor::TimeCenter
StateDescriptor::timeType () const
{
    return t_type;
}

int
StateDescriptor::nComp () const
{
    return ncomp;
}

int
StateDescriptor::nExtra () const
{
    return ngrow;
}

const std::string&
StateDescriptor::name (int i) const
{
    return names[i];
}

const BCRec&
StateDescriptor::getBC (int i) const
{
    return bc[i];
}

int StateDescriptor::get_interp_type(int i) const {

 return interp_type_array[i];
}

const Array<BCRec>&
StateDescriptor::getBCs () const
{
    return bc;
}

int
StateDescriptor::inRange (int sc, int nc) const
{
    return sc>=0 && sc+nc<=ncomp;
}

void
StateDescriptor::define (IndexType                   btyp,
                         StateDescriptor::TimeCenter ttyp,
                         int                         ident,
                         int                         nextra,
                         int                         num_comp)
{
    type     = btyp;
    t_type   = ttyp;
    id       = ident;
    ngrow    = nextra;
    ncomp    = num_comp;

    BL_ASSERT (num_comp > 0);
   
    names.resize(num_comp);
    bc.resize(num_comp);
    interp_type_array.resize(num_comp);
}

void
StateDescriptor::setComponent (int                               comp,
                               const std::string&                nm,
                               const BCRec&                      bcr,
                               int interp_type)
{
    BL_ASSERT(comp >= 0 && comp < ncomp && names[comp].empty());
    names[comp]       = nm;
    bc[comp]          = bcr;
    interp_type_array[comp]=interp_type;
}

void
StateDescriptor::dumpNames (std::ostream& os,
                            int           start_comp,
                            int           num_comp) const
{
    BL_ASSERT(start_comp >= 0 && start_comp+num_comp <= ncomp);

    for (int k = 0; k < num_comp; k++)
    {
        os << names[start_comp+k] << ' ';
    }
}

