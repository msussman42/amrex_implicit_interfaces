
#ifndef AMREX_StateDescriptor_H_
#define AMREX_StateDescriptor_H_ 

#include <utility>
#include <memory>

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_BCRec.H>
#include <AMReX_ArrayLim.H>

namespace amrex {

class Interpolater;

// FILL or GROUP FILL
extern "C"
{
    typedef void (*BndryFuncDefaultSUSSMAN)(
      int* tid_in,
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(lo), AMREX_ARLIM_P(hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc_local,
      const int* scomp,int* ncomp_local,int* bfact);
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
    // Class wrapper around BndryFuncDefaultSUSSMAN.
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
        BndryFunc (BndryFuncDefaultSUSSMAN inFunc);
        //
        // Another Constructor.
        //
        BndryFunc (BndryFuncDefaultSUSSMAN inFunc,BndryFuncDefaultSUSSMAN gFunc);
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
           int* tid_in,
           int* grid_type,
           int* level,
           Real* data, const int* lo, const int* hi,
           const int* dom_lo, const int* dom_hi,
           const Real* dx, const Real* grd_lo,
           const Real* time, const int* bc_local,
           const int* scomp,int* ncomp_local,int* bfact) const;
        //
        // Fill boundary cells using "group" function.
        //
        // FILL or GROUP FILL
        virtual void operator () (
          int* tid_in,
          int* grid_type,
          int* level,
          Real* data, const int* lo, const int* hi,
          const int* dom_lo, const int* dom_hi,
          const Real* dx, const Real* grd_lo,
          const Real* time, const int* bc_local,
          const int* scomp,int* ncomp_local,int* bfact,
          bool) const;
    private:

        BndryFuncDefaultSUSSMAN  m_func;
        BndryFuncDefaultSUSSMAN  m_gfunc;
    };
    //
    // The default constructor.
    //
    StateDescriptor () noexcept;
    //
    // Constructor that sets all data members.
    //
    StateDescriptor (IndexType     btyp,
                     int           ident,
                     int           nextra,
                     int           num_comp,
                     Interpolater* interp,
		     int           state_holds_data_in,
		     int           blocking_in); 
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
		 int           state_holds_data_in,
		 int           blocking_in);  
    //
    // Sets details of comp component.
    //
    void setComponent (int                comp,
                       const std::string& nm,
                       const BCRec&       bc,
                       const BndryFunc&   func,
                       Interpolater*   interp_in=0,  
                       int             max_map_start_comp=-1, 
                       int             min_map_end_comp=-1 
                      );  
    //
    // Sets details of comp component.
    //
    void setComponent (int                comp,
                       const std::string& nm,
                       const BCRec&       bc,
                       const BndryFunc&   func,
                       Interpolater*      interp_in,
                       bool               master_or_slave,
                       int                groupsize_in
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
    IndexType getType () const noexcept;
    //


    int nComp () const noexcept;
    //
    // Returns the grow factor.
    //
    int nExtra () const noexcept;
    //
    // Returns the interpolater.
    //
    Interpolater* interp () const noexcept;
    //
    // Returns the interpolater of specified component.
    //
    Interpolater* interp (int i) const noexcept;
    //
    // Returns the name of specified component.
    //
    const std::string& name (int i) const noexcept;
    //
    // Returns the BCRec of specified component.
    //
    const BCRec& getBC (int i) const noexcept;
    //
    // Returns all BCRecs.
    //
    const Vector<BCRec>& getBCs () const noexcept;
    //
    // Returns the BndryFunc of specified component.
    //
    const BndryFunc& bndryFill (int i) const noexcept;

    void check_inRange (Vector<int> scompBC_map, int nc) const;

    std::vector< std::pair<int,int> > 
      sameInterps (Vector<int> scompBC_map, int ncomp_in) const;

    int get_state_holds_data () const noexcept;
    int get_state_blocking () const noexcept;

    bool master (int i) const { return m_master[i]; }

    int groupsize (int i) const { return m_groupsize[i]; }

private:

    IndexType          type;     // Cell centered, node centered ...
    int                id;       // Unique id
    int                ncomp;    // Number of components
    int                ngrow;    // Grow factor
    Interpolater*      mapper;   // Default interpolator
    int                state_holds_data; //used by "get_state_holds_data()"
    int                state_blocking;
    Vector<std::string> names;    // Printable names of components
    Vector<BCRec>       bc;       // Array of bndry types for entire level
     //!< Array of pointers to bndry fill functions
    Vector<std::unique_ptr<BndryFunc> >  bc_func; 
    std::vector<bool>  m_master; // Are we a master or slave?
    Vector<int>         m_groupsize; // Groupsize if we're a master
    //
    // If mapper_comp[icomp] != 0, that map is used instead of mapper
    // when interpolating icomp; otherwise, mapper is used.
    //
    Vector<Interpolater*> mapper_comp;
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
    Vector<int> max_map_start_comp;
    Vector<int> min_map_end_comp;

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
    DescriptorList () noexcept;
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
	  int                         state_holds_data_in,
	  int                         blocking_in); 
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
          Interpolater*                 interp_in=0,  
          int                           max_map_start_comp=-1, 
          int                           min_map_end_comp=-1 
         );
    //
    // Calls setComponent() on StateDescriptor at index indx on group.
    //
    void setComponent (int              indx,
      int                               comp,
      const Vector<std::string>&         nm,
      const Vector<BCRec>&               bc,
      const StateDescriptor::BndryFunc& func,
      Interpolater*                     interp_in=0
     );
    //
    // Returns StateDescriptor at index k.
    //
    const StateDescriptor& operator[] (int k) const noexcept;

private:
    //
    // These are disallowed.
    //
    DescriptorList (const DescriptorList&);
    DescriptorList& operator= (const DescriptorList&);

    Vector<std::unique_ptr<StateDescriptor> > desc;
};

} // namespace amrex

#endif /*AMREX_StateDescriptor_H_*/
