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
// $Id: Derive.cpp,v 1.16 2001/08/09 22:42:00 marc Exp $
//
#include <winstd.H>

#if defined(BL_OLD_STL)
#include <string.h>
#else
#include <cstring>
#endif

#include <Derive.H>
#include <StateDescriptor.H>

DeriveRec::DeriveRec (const std::string& name,
                      IndexType      result_type,
                      int            nvar_derive,
                      DeriveFunc     der_func,
                      DeriveBoxMap   box_map)
    :
    derive_name(name),
    der_type(result_type),
    n_derive(nvar_derive),
    variable_names(),
    func(der_func),
    bx_map(box_map),
    n_state(0),
    nsr(0),
    rng(0),
    bcr(0)
{}

DeriveRec::DeriveRec (const std::string& name,
                      IndexType      result_type,
                      int            nvar_derive,
		      Array<std::string>& var_names,
                      DeriveFunc     der_func,
                      DeriveBoxMap   box_map)
    :
    derive_name(name),
    der_type(result_type),
    n_derive(nvar_derive),
    variable_names(var_names),
    func(der_func),
    bx_map(box_map),
    n_state(0),
    nsr(0),
    rng(0),
    bcr(0)
{}

DeriveRec::~DeriveRec () 
{
   delete [] bcr;
   func     = 0;
   bx_map   = 0;
   while (rng != 0)
   {
       StateRange* r = rng;
       rng = rng->next;
       delete r;
   }
}

const std::string&
DeriveRec::name () const
{
    return derive_name;
}

IndexType
DeriveRec::deriveType () const
{
    return der_type;
}

DeriveFunc
DeriveRec::derFunc () const
{
    return func;
}

DeriveRec::DeriveBoxMap
DeriveRec::boxMap () const
{
    return bx_map;
}

int
DeriveRec::numDerive () const
{
    return n_derive;
}

int
DeriveRec::numRange () const
{
    return nsr;
}

int
DeriveRec::numState () const
{
    return n_state;
}

const int*
DeriveRec::getBC () const
{
    return bcr;
}

void
DeriveRec::addRange (const DescriptorList& d_list,
                     int                   state_indx,
                     int                   src_comp,
                     int                   num_comp) 
{
    const StateDescriptor& d = d_list[state_indx];

    StateRange* r = new StateRange;

    r->typ  = state_indx;
    r->sc   = src_comp;
    r->nc   = num_comp;
    r->next = 0;
    //
    // Add to end of list.
    //
    if (rng == 0)
    {
        rng = r;
    }
    else
    {
        StateRange* prev = rng;
        while (prev->next != 0)
            prev = prev->next;
        prev->next = r;
    }
    nsr++;
    n_state += num_comp;

    buildBC(d_list);
}

void
DeriveRec::getRange (int  k,
                     int& state_indx,
                     int& src_comp,
                     int& num_comp) const
{
    StateRange* r;

    for (r = rng; r != 0 && k > 0; k--, r = r->next)
        ;
    BL_ASSERT(r != 0);
    state_indx = r->typ;
    src_comp   = r->sc;
    num_comp   = r->nc;
}

void
DeriveRec::buildBC (const DescriptorList& d_list)
{
    BL_ASSERT(nsr > 0);
    delete [] bcr;
    bcr = new int[2*BL_SPACEDIM*n_state];
    int* bci = bcr;
    for (DeriveRec::StateRange* r = rng; r != 0; r = r->next)
    {
        const StateDescriptor& d = d_list[r->typ];

        for (int k = 0; k < r->nc; k++)
        {
            const int* bc = d.getBC(r->sc + k).vect();

            for (int j = 0; j < 2*BL_SPACEDIM; j++)
            {
                bci[j] = bc[j];
            }
            bci += 2*BL_SPACEDIM;
        }
    }
}

const
std::string&
DeriveRec::variableName(int comp) const
{
  if (comp < variable_names.size()) 
     return variable_names[comp];

  return derive_name;
}

DeriveList::DeriveList () {}

void
DeriveList::add (const std::string&      name,
                 IndexType               result_type,
                 int                     nvar_der,
                 DeriveFunc              der_func,
                 DeriveRec::DeriveBoxMap bx_map)
{
    lst.push_back(DeriveRec(name,result_type,nvar_der,der_func,bx_map));
}

void
DeriveList::add (const std::string&      name,
                 IndexType               res_typ,
                 int                     nvar_der,
                 Array<std::string>&     vars,
                 DeriveFunc              der_func,
                 DeriveRec::DeriveBoxMap bx_map)
{
    lst.push_back(DeriveRec(name,res_typ,nvar_der,vars,der_func,bx_map));
}

std::list<DeriveRec>&
DeriveList::dlist ()
{
    return lst;
}

bool
DeriveList::canDerive (const std::string& name) const 
{
    for (std::list<DeriveRec>::const_iterator li = lst.begin();
         li != lst.end();
         ++li)
    {
        if (li->derive_name == name)
            return true;
    }
    return false;
}

const DeriveRec*
DeriveList::get (const std::string& name) const
{
    for (std::list<DeriveRec>::const_iterator li = lst.begin();
         li != lst.end();
         ++li)
    {
        if (li->derive_name == name)
            return &(*li);
    }
    return 0;
}

void
DeriveList::addComponent (const std::string&    name,
                          const DescriptorList& d_list, 
                          int                   state_indx,
                          int                   s_comp,
                          int                   n_comp)
{
    std::list<DeriveRec>::iterator li = lst.begin();

    for ( ; li != lst.end(); ++li)
    {
        if (li->derive_name == name)
            break;
    }

    BL_ASSERT (li != lst.end());

    li->addRange(d_list, state_indx, s_comp, n_comp);
}

