
#ifndef _StateDescriptor_H_
#define _StateDescriptor_H_ 

#include <utility>
#include <vector>

#include <Array.H>
#include <REAL.H>
#include <Box.H>
#include <BCRec.H>
#include <ArrayLim.H>

class Interpolater;

// FILL or GROUP FILL
extern "C"
{
    typedef void (*BndryFuncDefault)(
      int* level,
      Real* data, ARLIM_P(lo), ARLIM_P(hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);
}

//
// Attributes of StateData.
//

class StateDescriptor
{
    friend class DescriptorList;

public:
    //
    // Type of function called by BCRec for user-supplied boundary data.
    //
    // Class wrapper around BndryFuncDefault.
    //
    class BndryFunc
    {
    public:
        //
        // Bogus constructor.
        //
        BndryFunc ();
        //
        // A Constructor.
        //
        BndryFunc (BndryFuncDefault inFunc);
        //
        // Another Constructor.
        //
        BndryFunc (BndryFuncDefault inFunc,BndryFuncDefault gFunc);
        //
        // Return a ptr to a clone of this object.
        // It is the responsibility of the caller to delete the result.
        //
        virtual BndryFunc* clone () const;
        //
        // Destructor.
        //
        virtual ~BndryFunc ();
        //
        // Fill boundary cells using "regular" function.
        //
        // FILL or GROUP FILL
        virtual void operator () (
           int* level,
           Real* data, const int* lo, const int* hi,
           const int* dom_lo, const int* dom_hi,
           const Real* dx, const Real* grd_lo,
           const Real* time, const int* bc,
           const int* scomp,int* ncomp,int* bfact) const;
        //
        // Fill boundary cells using "group" function.
        //
        // FILL or GROUP FILL
        virtual void operator () (
          int* level,
          Real* data, const int* lo, const int* hi,
          const int* dom_lo, const int* dom_hi,
          const Real* dx, const Real* grd_lo,
          const Real* time, const int* bc,
          const int* scomp,int* ncomp,int* bfact,
          bool) const;
    private:

        BndryFuncDefault  m_func;
        BndryFuncDefault  m_gfunc;
    };
    //
    // The default constructor.
    //
    StateDescriptor ();
    //
    // Constructor that sets all data members.
    //
    StateDescriptor (IndexType     btyp,
                     int           ident,
                     int           nextra,
                     int           num_comp,
                     Interpolater* interp,
                     bool          store_in_checkpoint); // default true
    //
    // The destructor.
    //
    ~StateDescriptor ();
    //
    // Define the data members if constructed with default construtor.
    //
    void define (IndexType     btyp,
                 int           ident,
                 int           nextra,
                 int           num_comp,
                 Interpolater* interp,
                 bool          store_in_checkpoint);  // default true
    //
    // Sets details of comp component.
    //
    void setComponent (int                comp,
                       const std::string& nm,
                       const BCRec&       bc,
                       const BndryFunc&   func,
                       Interpolater*      interp=0,  
                       int                max_map_start_comp=-1, 
                       int                min_map_end_comp=-1 
                      );  
    //
    // Sets details of comp component.
    //
    void setComponent (int                comp,
                       const std::string& nm,
                       const BCRec&       bc,
                       const BndryFunc&   func,
                       Interpolater*      interp,
                       bool               master_or_slave,
                       int                groupsize
                      );
    //
    // Resets boundary condition data for comp component.
    //
    void resetMapper(int comp,Interpolater* interp);
    void saveMapper(int comp,Interpolater*& interp);

    void reset_bcrecs(int comp,BCRec bcr);
    void save_bcrecs_statedesc(int comp,BCRec& bcr);

    void resetComponentBCs (int              comp,
                            const BCRec&     bcr,
                            const BndryFunc& func);
    //
    // Set interpolaters for a subset of the state vector components.
    //
    void setUpMaps (int&                use_default_map,
                    const Interpolater* default_map,
                    int                 start_comp,
                    int                 num_comp,
                    Interpolater**&     maps, 
                    int&                nmaps,
                    int*&               map_start_comp, 
                    int*&               map_num_comp,
                    int*&               max_start_comp,
                    int*&               min_end_comp) const;
    //
    // Cleanup interpolaters for a subset of the state vector components.
    //
    void cleanUpMaps (Interpolater**& maps, 
                      int*&           map_start_comp,
                      int*&           map_num_comp,
                      int*&           max_start_comp,
                      int*&           min_end_comp) const;
    //
    // Output names of components.
    //
    void dumpNames (std::ostream& os,
                    int           start_comp,
                    int           num_comp) const;
    //
    // Returns the IndexType.
    //
    IndexType getType () const;
    //


    int nComp () const;
    //
    // Returns the grow factor.
    //
    int nExtra () const;
    //
    // Returns the interpolater.
    //
    Interpolater* interp () const;
    //
    // Returns the interpolater of specified component.
    //
    Interpolater* interp (int i) const;
    //
    // Returns the name of specified component.
    //
    const std::string& name (int i) const;
    //
    // Returns the BCRec of specified component.
    //
    const BCRec& getBC (int i) const;
    //
    // Returns all BCRecs.
    //
    const Array<BCRec>& getBCs () const;
    //
    // Returns the BndryFunc of specified component.
    //
    const BndryFunc& bndryFill (int i) const;

    void check_inRange (Array<int> scompBC_map, int nc) const;

    std::vector< std::pair<int,int> > 
      sameInterps (Array<int> scompBC_map, int ncomp) const;

    //
    // Should store this StateData in the checkpoint file
    //
    bool store_in_checkpoint () const;

    bool master (int i) const { return m_master[i]; }

    int groupsize (int i) const { return m_groupsize[i]; }

private:

    IndexType          type;     // Cell centered, node centered ...
    int                id;       // Unique id
    int                ncomp;    // Number of components
    int                ngrow;    // Grow factor
    Interpolater*      mapper;   // Default interpolator
    bool               m_store_in_checkpoint; // store in checkpoint file?
    Array<std::string> names;    // Printable names of components
    Array<BCRec>       bc;       // Array of bndry types for entire level
     //!< Array of pointers to bndry fill functions
    Vector<std::unique_ptr<BndryFunc> >  bc_func; 
    std::vector<bool>  m_master; // Are we a master or slave?
    Array<int>         m_groupsize; // Groupsize if we're a master
    //
    // If mapper_comp[icomp] != 0, that map is used instead of mapper
    // when interpolating icomp; otherwise, mapper is used.
    //
    Array<Interpolater*> mapper_comp;
    //
    // For component icomp, max_map_start_comp[icomp] and
    // min_map_end_comp[icomp] represent the maximum starting index and
    // minimum ending index with which to use mapper_comp[icomp] when
    // interpolating component icomp.  These are convenient for using
    // with "coupled" multiple component maps like CellConservativeLinear.
    //
    // If there is no need to "couple" a component, max_map_start_comp[icomp]
    // and min_map_end_comp[icomp] simply equal icomp.
    //
    Array<int> max_map_start_comp;
    Array<int> min_map_end_comp;

};

//
//@Man:
//@Memo: A List of StateDescriptors.
/*@Doc:

  A container class for StateDescriptors.
*/

class DescriptorList
{
public:
    //
    // The constructor.
    //
    DescriptorList ();
    //
    // Set the list to its default state.
    //
    void clear ();
    //
    // Returns number of elements in the list.
    //
    int size () const;
    //
    // Adds new StateDescriptor at index indx to list.
    //
    void addDescriptor (int           indx,
          IndexType                   typ,
          int                         nextra,
          int                         num_comp,
          Interpolater*               interp,
          bool                        store_in_checkpoint); // default true
    //
    // Calls resetComponentBCs() on StateDescriptor at index indx.
    //
    void resetMapper(int indx,int comp,Interpolater* interp);
    void saveMapper(int indx,int comp,Interpolater*& interp);

    void reset_bcrecs(int indx,int comp,BCRec bcr);
    void save_bcrecs_statedesc(int indx,int comp,BCRec& bcr);

    void resetComponentBCs (int                               indx,
                            int                               comp,
                            const BCRec&                      bc,
                            const StateDescriptor::BndryFunc& func);
    //
    // Calls setComponent() on StateDescriptor at index indx.
    //
    void setComponent (int                  indx,
          int                               comp,
          const std::string&                nm,
          const BCRec&                      bc,
          const StateDescriptor::BndryFunc& func,
          Interpolater*                     interp=0,  
          int                               max_map_start_comp=-1, 
          int                               min_map_end_comp=-1 
         );
    //
    // Calls setComponent() on StateDescriptor at index indx on group.
    //
    void setComponent (int                               indx,
      int                               comp,
      const Array<std::string>&         nm,
      const Array<BCRec>&               bc,
      const StateDescriptor::BndryFunc& func,
      Interpolater*                     interp=0
     );
    //
    // Returns StateDescriptor at index k.
    //
    const StateDescriptor& operator[] (int k) const;

private:
    //
    // These are disallowed.
    //
    DescriptorList (const DescriptorList&);
    DescriptorList& operator= (const DescriptorList&);

    Vector<std::unique_ptr<StateDescriptor> > desc;
};

#endif /*_StateDescriptor_H_*/
