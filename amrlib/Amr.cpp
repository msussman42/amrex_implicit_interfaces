
#include <algorithm>
#include <cstdio>
#include <list>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <AMReX_Geometry.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_CoordSys.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BoxDomain.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_TagBox.H>

#include <AMReX_AmrParGDB.H>

#include <AMReX_Cluster.H>
#include <LevelBld.H>
#include <AmrLevel.H>
#include <Amr.H>
#include <StateData.H>

#include <INTERP_F.H>

#include <DRAG_COMP.H>
#include <EXTRAP_COMP.H>
#include <GLOBALUTIL_F.H>

#ifdef BL_USE_ARRAYVIEW
#include <DatasetClient.H>
#endif

namespace amrex {

//
// Static class members.  
// Set defaults in Initialize()!!!
//
Vector<BoxArray>       Amr::initial_ba;
Vector<BoxArray>       Amr::regrid_ba;


namespace
{
    const std::string CheckPointVersion("CheckPointVersion_1.0");

    bool initialized = false;
}

//namespace
//{
    //
    // These are all ParmParse'd in.  Set defaults in Initialize()!!!
    //
    int  plot_nfiles;
    int  mffile_nstreams;
    int  probinit_natonce;
    int  checkpoint_nfiles;
    int  regrid_on_restart;
    int  use_efficient_regrid;
    bool refine_grid_layout;
    int  plotfile_on_restart;
    int  insitu_on_restart;
    int  checkpoint_on_restart;
    bool checkpoint_files_output;
    int  compute_new_dt_on_regrid;
    bool precreateDirectories;
    bool prereadFAHeaders;
    VisMF::Header::Version plot_headerversion(VisMF::Header::Version_v1);
    VisMF::Header::Version checkpoint_headerversion(VisMF::Header::Version_v1);
    Vector<int> recalesce_flag;  //if set, then "recalesce_state" checkpointed.

//}

int
Amr::Old_blockingFactor(int lev) const
{
 int local_blocking=blocking_factor[lev][0];
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if (blocking_factor[lev][dir]==local_blocking) {
   // do nothing
  } else
   amrex::Error("expecting uniform blocking factor");
 }
 return local_blocking;
}


int
Amr::Old_maxGridSize(int lev) const
{
 int local_max_grid_size=max_grid_size[lev][0];
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if (max_grid_size[lev][dir]==local_max_grid_size) {
   // do nothing
  } else
   amrex::Error("expecting uniform dir=0,1,2 max_grid_size");
 }
 return local_max_grid_size;
}


int 
Amr::Time_blockingFactor () const 
{   
    return time_blocking_factor;
}
   
int 
Amr::get_MAX_NUM_SLAB () const 
{   
    return MAX_NUM_SLAB;
}
   
int 
Amr::get_slab_dt_type () const 
{   
    return slab_dt_type;
}   

void
Amr::Initialize ()
{
    if (initialized) return;
    //
    // Set all defaults here!!!
    //
    plot_nfiles              = 64;
    mffile_nstreams          = 1;
    probinit_natonce         = 512;
    checkpoint_nfiles        = 64;
    regrid_on_restart        = 0;
    use_efficient_regrid     = 0;
    plotfile_on_restart      = 0;
    checkpoint_on_restart    = 0;
    checkpoint_files_output  = true;
    compute_new_dt_on_regrid = 0;
    precreateDirectories     = true;
    prereadFAHeaders         = true;
    plot_headerversion       = VisMF::Header::Version_v1;
    checkpoint_headerversion = VisMF::Header::Version_v1;

    amrex::ExecOnFinalize(Amr::Finalize);

    initialized = true;
}

void
Amr::Finalize ()
{
    Amr::regrid_ba.clear();
    Amr::initial_ba.clear();

    initialized = false;
}


std::ostream&
Amr::DataLog (int i)
{
    return *datalog[i];
}

int 
Amr::AMR_recalesce_flag(int im) const {

 if ((im<1)||(im>global_AMR_num_materials))
  amrex::Error("im out of range");

 return recalesce_flag[im-1];

}


Vector<std::unique_ptr<AmrLevel> >&
Amr::getAmrLevels () noexcept
{
    return amr_level;
}

// 1. AmrMesh():
//  a. Geometry::Setup()  (rb==nullptr, coord=-1, is_per=nullptr)
//   Setup Geometry from ParmParse file.
//   (might be needed for variableSetup or getLevelBld)
//  b. InitAmrMesh
//     i. geom[i].define(index_domain)  i=0...max_level
// 2. Initialize()
// 3. InitAmr()
//  a. levelbld = getLevelBld();
Amr::Amr () 
  : AmrMesh() {

     // init default values for some parameters.
    Initialize();
     // levelbld = getLevelBld();
     // ...
    InitAmr();

} // end subroutine Amr::Amr () 

void
Amr::InitAmr () {

    AMR_max_phase_change_rate.resize(AMREX_SPACEDIM);
    AMR_min_phase_change_rate.resize(AMREX_SPACEDIM);
    for (int j=0;j<AMREX_SPACEDIM;j++) {
     AMR_max_phase_change_rate[j]=0.0;
     AMR_min_phase_change_rate[j]=0.0;
    }
    
    AMR_volume_history_recorded=0;
    AMR_volume_history.resize(1);
    AMR_volume_history[0]=0.0;

    //
    // Determine physics class.
    //
    std::fflush(NULL);
    if (1==1) {
     std::cout << "prior to levelbld = getLevelBld() on processor " <<
         ParallelDescriptor::MyProc() << "\n";
    }
    std::fflush(NULL);

    //LevelBld* levelbld
    //"getLevelBld" is declared in: NS_fluids_lib/NSBld.cpp
    levelbld = getLevelBld();
    //
    // Global function that define state variables.
    //
    std::fflush(NULL);
    if (1==1) {
     std::cout << "levelbld->variableSetUp() on processor " <<
         ParallelDescriptor::MyProc() << "\n";
    }
    std::fflush(NULL);
    levelbld->variableSetUp();
    //
    // Set default values.
    //
    grid_eff               = 0.7;
    plot_int               = -1;
    slice_int              = -1;
    n_proper               = 1;
    last_plotfile          = 0;
    last_checkpoint        = 0;
    record_run_info        = false;
    record_grid_info       = false;
    file_name_digits       = 5;
    record_run_info_terse  = false;

    int i;
    for (i = 0; i < AMREX_SPACEDIM; i++)
     isPeriodic[i] = false;

    ParmParse ppns("ns");
    ParmParse pp("amr");
    //
    // Check for command line flags.
    //
    verbose = 0;
    pp.query("v",verbose);

    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "Amr.verbose= " << verbose << '\n';
    }

    ppns.get("num_materials",global_AMR_num_materials);
    if ((global_AMR_num_materials<2)||(global_AMR_num_materials>999))
     amrex::Error("global_AMR_num_materials invalid");

    global_AMR_num_materials_viscoelastic=0;

    Vector<Real> elastic_viscosity_temp;
    Vector<Real> elastic_time_temp;
    Vector<int> viscoelastic_model_temp;
    Vector<int> store_elastic_data_temp;
    elastic_viscosity_temp.resize(global_AMR_num_materials);
    elastic_time_temp.resize(global_AMR_num_materials);
    viscoelastic_model_temp.resize(global_AMR_num_materials);
    store_elastic_data_temp.resize(global_AMR_num_materials);

    for (int im=0;im<global_AMR_num_materials;im++) {
     elastic_viscosity_temp[im]=0.0;
     elastic_time_temp[im]=0.0;
     viscoelastic_model_temp[im]=0;
     store_elastic_data_temp[im]=0;
    }
    ppns.queryarr("elastic_viscosity",elastic_viscosity_temp,0,
		  global_AMR_num_materials);
    ppns.queryarr("elastic_time",elastic_time_temp,0,
		  global_AMR_num_materials);
    ppns.queryarr("viscoelastic_model",viscoelastic_model_temp,0,
		  global_AMR_num_materials);

    for (int im=0;im<global_AMR_num_materials;im++) {
     if (elastic_viscosity_temp[im]>0.0) {
      if (fort_built_in_elastic_model(&elastic_viscosity_temp[im],
           &viscoelastic_model_temp[im])==1) {
       store_elastic_data_temp[im]=1;
      } else if (fort_built_in_elastic_model(&elastic_viscosity_temp[im],
           &viscoelastic_model_temp[im])==0) {
       // do nothing
      } else
       amrex::Error("fort_is_eulerian_elastic_model invalid");
     } else if (elastic_viscosity_temp[im]==0.0) {
      // do nothing
     } else
      amrex::Error("elastic_viscosity_temp[im] invalid");
    } // im=0..global_AMR_num_materials-1 

    for (int im=0;im<global_AMR_num_materials;im++) {
     if (store_elastic_data_temp[im]==1) {
      global_AMR_num_materials_viscoelastic++;
     } else if (store_elastic_data_temp[im]==0) {
      // do nothing
     } else
      amrex::Error("store_elastic_data_temp invalid");
    } // im=0..global_AMR_num_materials-1 

    recalesce_flag.resize(global_AMR_num_materials);
    for (int im=0;im<global_AMR_num_materials;im++) {
     recalesce_flag[im]=0;
    }
    ppns.queryarr("recalesce_flag",recalesce_flag,0,global_AMR_num_materials);
    for (int im=0;im<global_AMR_num_materials;im++) {
     if ((recalesce_flag[im]!=0)&&
         (recalesce_flag[im]!=1)&&
         (recalesce_flag[im]!=2))
      amrex::Error("recalesce_flag invalid");
    }

    pp.query("regrid_on_restart",regrid_on_restart);
    if ((regrid_on_restart!=0)&&(regrid_on_restart!=1))
     amrex::Error("regrid_on_restart invalid");

    pp.query("plotfile_on_restart",plotfile_on_restart);
    pp.query("checkpoint_on_restart",checkpoint_on_restart);

    pp.query("checkpoint_files_output", checkpoint_files_output);

    plot_nfiles = ParallelDescriptor::NProcs();
    checkpoint_nfiles = ParallelDescriptor::NProcs();

    pp.query("refine_grid_layout", refine_grid_layout);

    pp.query("mffile_nstreams", mffile_nstreams);
    pp.query("probinit_natonce", probinit_natonce);

    probinit_natonce = 
     std::max(1, std::min(ParallelDescriptor::NProcs(), probinit_natonce));

    pp.query("file_name_digits", file_name_digits);

    if (pp.contains("run_log")) {
     std::string log_file_name;
     pp.get("run_log",log_file_name);
     setRecordRunInfo(log_file_name);
    }
    if (pp.contains("run_log_terse")) {
     std::string log_file_name;
     pp.get("run_log_terse",log_file_name);
     setRecordRunInfoTerse(log_file_name);
    }
    if (pp.contains("grid_log")) {
     std::string grid_file_name;
     pp.get("grid_log",grid_file_name);
     setRecordGridInfo(grid_file_name);
    }

    if (pp.contains("data_log")) {
     int num_datalogs = pp.countval("data_log");
     datalog.resize(num_datalogs);
     Vector<std::string> data_file_names(num_datalogs);
     pp.queryarr("data_log",data_file_names,0,num_datalogs);
     for (int i_data = 0; i_data < num_datalogs; i_data++) 
      setRecordDataInfo(i_data,data_file_names[i_data]);
    }

    //
    // Restart or run from scratch?
    //
    pp.query("restart", restart_file);
    int nlev     = max_level+1;

    level_cells_advanced.resize(nlev);
    level_steps.resize(nlev);
    level_count.resize(nlev);
    space_blocking_factor.resize(nlev);
    amr_level.resize(nlev);
    //
    // Set bogus values.
    //

    dt_AMR=1.0;

    time_blocking_factor = 1;
    MAX_NUM_SLAB=33;
    slab_dt_type=0; // 0=SEM 1=evenly spaced

    for (i = 0; i < nlev; i++) {
     level_cells_advanced[i] = 0.0;
     level_steps[i] = 0;
     level_count[i] = 0;
     space_blocking_factor[i] = 1;
    }

    regrid_int=0;

    //
    // Read other amr specific values.
    //
    check_file_root = "chk";
    pp.query("check_file",check_file_root);

    check_int = -1;
    int got_check_int = pp.query("check_int",check_int);

    check_per = -1.0;
    int got_check_per = pp.query("check_per",check_per);

    if (got_check_int == 1 && got_check_per == 1) {
     amrex::Error("Must only specify amr.check_int OR amr.check_per");
    }

    plot_file_root = "plt";
    pp.query("plot_file",plot_file_root);

    plot_int = -1;
    int got_plot_int = pp.query("plot_int",plot_int);

    plot_per = -1.0;
    int got_plot_per = pp.query("plot_per",plot_per);

    if (got_plot_int == 1 && got_plot_per == 1) {
     amrex::Error("Must only specify amr.plot_int OR amr.plot_per");
    }

    slice_int=-1;

    if (got_plot_int==1) {
     slice_int=plot_int;

     int got_slice_int=pp.query("slice_int",slice_int);
     if ((got_slice_int!=0)&&(got_slice_int!=1))
      amrex::Error("got_slice_int invalid");

     if (slice_int>plot_int)
      amrex::Error("slice_int should be less than or equal to plot_int");
    }

    pp.queryarr("space_blocking_factor",space_blocking_factor);
    pp.query("time_blocking_factor",time_blocking_factor);
    pp.query("MAX_NUM_SLAB",MAX_NUM_SLAB);
    if (time_blocking_factor+1>MAX_NUM_SLAB)
     amrex::Error("MAX_NUM_SLAB too small");
    if (time_blocking_factor<1)
     amrex::Error("time_blocking_factor too small");

    pp.query("slab_dt_type",slab_dt_type);
    if ((slab_dt_type!=0)&&
        (slab_dt_type!=1))
     amrex::Error("slab_dt_type invalid");

    //
    // SUSSMAN: just one regrid_int value
    //
    if (max_level > 0) {

     int numvals = pp.countval("regrid_int");
     if (numvals == 1) {
      pp.query("regrid_int",regrid_int);
     } else {
      amrex::Error("specify just one regrid_int value");
     }

    }

    //
    // SUSSMAN:
    // Now check offset; CoordSys does not need it anymore though.
    //
    Vector<int> n_cell(AMREX_SPACEDIM);
    pp.getarr("n_cell",n_cell,0,AMREX_SPACEDIM);
    IntVect lo(IntVect::TheZeroVector()), hi(n_cell);
    hi -= IntVect::TheUnitVector();

    Real offset[AMREX_SPACEDIM];
    for (i = 0; i < AMREX_SPACEDIM; i++) {
     const Real delta = geom[0].ProbLength(i)/(Real)n_cell[i];
     offset[i]        = geom[0].ProbLo(i) + delta*lo[i];
     if ((lo[i]==0)&&(delta>0.0)) {
      // do nothing
     } else
      amrex::Error("expecting lo=0 and delta>0");
    }

    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "Amr.time_blocking_factor= " <<
        time_blocking_factor << '\n';
     std::cout << "Amr.MAX_NUM_SLAB= " <<
        MAX_NUM_SLAB << '\n';
     std::cout << "Amr.slab_dt_type= " <<
        slab_dt_type << '\n';
    }

    ppns.get("num_species_var",global_AMR_num_species_var);
    if ((global_AMR_num_species_var<0)||(global_AMR_num_species_var>999))
     amrex::Error("global_AMR_num_species_var invalid");

    global_AMR_particles_flag=0;
    ppns.query("particles_flag",global_AMR_particles_flag);
    if ((global_AMR_particles_flag==0)||
        (global_AMR_particles_flag==1)) {
     //do nothing
    } else
     amrex::Error("global_AMR_particles_flag invalid");

     //these variable definitions are needed by "SOA_NCOMP"
    int num_materials=global_AMR_num_materials;
    int num_species_var=global_AMR_num_species_var;
    int num_materials_viscoelastic=global_AMR_num_materials_viscoelastic;

    global_AMR_num_SoA_var=SOA_NCOMP;

    std::fflush(NULL);
    std::cout << "global_AMR_num_SoA_var= " << global_AMR_num_SoA_var <<
	 " on processor " << ParallelDescriptor::MyProc() << "\n";
    std::fflush(NULL);
    std::fflush(NULL);
    std::cout << "num_materials= " << num_materials <<
	 " on processor " << ParallelDescriptor::MyProc() << "\n";
    std::fflush(NULL);
    std::fflush(NULL);
    std::cout << "num_species_var= " << num_species_var <<
	 " on processor " << ParallelDescriptor::MyProc() << "\n";
    std::fflush(NULL);
    std::fflush(NULL);
    std::cout << "num_materials_viscoelastic= "<<num_materials_viscoelastic<<
	 " on processor " << ParallelDescriptor::MyProc() << "\n";
    std::fflush(NULL);

    m_gdb.reset(new AmrParGDB(this));

} // subroutine InitAmr

Amr::~Amr ()
{
    if (level_steps[0] > last_checkpoint)
        checkPoint();

    if (level_steps[0] > last_plotfile) {
     int do_plot=1;
     int do_slice=((slice_int>0) ? 1 : 0);
     int SDC_outer_sweeps=0;
     int slab_step=Time_blockingFactor()-1;
     writePlotFile(
      do_plot,do_slice,
      SDC_outer_sweeps,slab_step);
    }

    levelbld->variableCleanUp();

    Amr::Finalize();
}

void
Amr::setRecordGridInfo (const std::string& filename)
{
    record_grid_info = true;
    if (ParallelDescriptor::IOProcessor())
    {
        gridlog.open(filename.c_str(),std::ios::out|std::ios::app);
        if (!gridlog.good())
            amrex::FileOpenFailed(filename);
    }
    ParallelDescriptor::Barrier("Amr::setRecordGridInfo");
}

void
Amr::setRecordRunInfo (const std::string& filename)
{
    record_run_info = true;
    if (ParallelDescriptor::IOProcessor())
    {
        runlog.open(filename.c_str(),std::ios::out|std::ios::app);
        if (!runlog.good())
            amrex::FileOpenFailed(filename);
    }
    ParallelDescriptor::Barrier("Amr::setRecordRunInfo");
}

void
Amr::setRecordRunInfoTerse (const std::string& filename)
{
    record_run_info_terse = true;
    if (ParallelDescriptor::IOProcessor())
    {
        runlog_terse.open(filename.c_str(),std::ios::out|std::ios::app);
        if (!runlog_terse.good())
            amrex::FileOpenFailed(filename);
    }
    ParallelDescriptor::Barrier("Amr::setRecordRunInfoTerse");
}

void
Amr::setRecordDataInfo (int i, const std::string& filename)
{
    if (ParallelDescriptor::IOProcessor())
    {
        datalog[i].reset(new std::fstream);
        datalog[i]->open(filename.c_str(),std::ios::out|std::ios::app);
        if (!datalog[i]->good())
            amrex::FileOpenFailed(filename);
    }
    ParallelDescriptor::Barrier("Amr::setRecordDataInfo");
}

void
Amr::setDt(Real dt)
{
 dt_AMR=dt;
}

int
Amr::okToContinue () noexcept
{
    int ok = true;
    for (int i = 0; ok && (i <= finest_level); i++)
        ok = ok && amr_level[i]->okToContinue();
    return ok;
}

void
Amr::writeDEBUG_PlotFile(int num,int SDC_outer_sweeps,int slab_step) {

 int do_plot=1;
 int do_slice=((slice_int>0) ? 1 : 0);
 writePlotFile(
   do_plot,do_slice,
   SDC_outer_sweeps,slab_step);

}

void
Amr::writePlotFile (int do_plot,int do_slice,
                    int SDC_outer_sweeps,int slab_step) {

 if ((do_plot==1)||
     ((do_plot==0)&&(do_slice==1))) {

   for (int k(0); k <= finest_level; ++k) {
    amr_level[k]->writePlotFile(
	do_plot,do_slice,
        SDC_outer_sweeps,slab_step);
   }

 } else if ((do_plot==0)&&(do_slice==0)) {
  // do nothing
 } else
  amrex::Error("do_plot or do_slice invalid");
  
}  // subroutine writePlotFile


void
Amr::AMR_checkInput ()
{
    FabArrayBase::Initialize();

    if (max_level < 0)
        amrex::Error("checkInput: max_level not set");
    //
    // 1. Check that blocking_factor is a power of 2 and no smaller than 2.
    // 2. Check that blocking_factor[i+1]<=blocking_factor[i].
    // 3. Check that blocking_factor[i]>=8 if i<max_level.
    //    (this last check insures that there are at least 4 coarse 
    //     (level i) proper nesting cells next to a (level i+1) finer level)
    for (int i = 0; i <= max_level; i++) {
        int k = Old_blockingFactor(i);
	if (i<max_level) {
         if (k<8) {
  	  std::cout << "must have at least 4 proper nesting cells" << '\n';
	  amrex::Error("must have blocking_factor>=8 if lev<max_level");
	 }
	 if (k<Old_blockingFactor(i+1))
	  amrex::Error("blocking_factor[i]<blocking_factor[i+1]");
	}
        if (k<2)
         amrex::Error("blocking factor must be 2 or larger on max_level");

        while ( k > 0 && (k%2 == 0) )
            k /= 2;
        if (k != 1)
            amrex::Error("Amr::checkInputs: blocking_factor not power of 2");
    } // i=0 .. max_level

    for (int i = 0; i <= max_level; i++) {
        int k = space_blocking_factor[i];

        if (k>Old_blockingFactor(i))
         amrex::Error("space_blocking_factor[i]>Old_blockingFactor(i)");
	if (i<max_level) {
	 if (k<space_blocking_factor[i+1])
	  amrex::Error("space_blocking_factor[i]<space_blocking_factor[i+1]");
	}

         // the number of coarse grid proper nesting cells for level i+1
         // is blocking_factor[i]/2
        if ((i>=0)&&(i<=max_level)) {
         if (Old_blockingFactor(i)<2*k)
          amrex::Error("bfact_grid>=2*space_blocking_Factor required");
	} else
	 amrex::Error("i invalid");

         //cannot have an element that is partially covered by a finer
         //grid.
	if ((i>=0)&&(i<max_level)) {
         if (Old_blockingFactor(i+1)<2*k)
          amrex::Error("(Old_blockingFactor(i+1)<2*k)");
        } else if (i==max_level) {
         // do nothing
	} else
	 amrex::Error("i invalid");

        while ( k > 0 && (k%2 == 0) )
            k /= 2;
        if (k != 1)
            amrex::Error("space_blocking_factor not power of 2");
    } // i=0..max_level

    for (int i = 0; i < max_level; i++) {
        int k = time_blocking_factor;
        if (k>Old_blockingFactor(0))
         amrex::Error("time_blocking_factor too big");
        while ( k > 0 && (k%2 == 0) )
            k /= 2;
        if (k != 1)
            amrex::Error("time_blocking_factor not power of 2");
    }


    //
    // Check level dependent values.
    //
    int i;

    const Box& domain = geom[0].Domain();
    if (!domain.ok())
        amrex::Error("level 0 domain bad or not set");
    //
    // Check that domain size is a multiple of blocking_factor[0].
    //
    for (i = 0; i < AMREX_SPACEDIM; i++) {
        int len = domain.length(i);
        if (len%Old_blockingFactor(0) != 0)
            amrex::Error("domain size not divisible by blocking_factor");
    }
    for (i = 0; i < AMREX_SPACEDIM; i++) {
        int len = domain.length(i);
        if (len%space_blocking_factor[0] != 0)
         amrex::Error("domain size not divisible by space_blocking_factor");
    }
    for (i = 0; i < AMREX_SPACEDIM; i++) {
        int len = domain.length(i);
        if (len%time_blocking_factor != 0)
         amrex::Error("domain size not divisible by time_blocking_factor");
    }
    //
    // Check that max_grid_size is even.
    //
    for (i = 0; i < max_level; i++)
    {
        if (Old_maxGridSize(i)%2 != 0)
            amrex::Error("max_grid_size is not even");
    }

    //
    // Check that max_grid_size is a multiple of blocking_factor at every level.
    //
    for (i = 0; i < max_level; i++) {
     if (Old_maxGridSize(i)%Old_blockingFactor(i) != 0)
      amrex::Error("max_grid_size not divisible by blocking_factor");
    }
    for (i = 0; i < max_level; i++) {
     if (Old_maxGridSize(i)%space_blocking_factor[i] != 0)
      amrex::Error("max_grid_size not divisible by space_blocking_factor");
    }
    for (i = 0; i < max_level; i++) {
     if (Old_maxGridSize(i)%time_blocking_factor != 0)
      amrex::Error("max_grid_size not divisible by time_blocking_factor");
    }

    if (!geom[0].ProbDomain().ok())
        amrex::Error("checkInput: bad physical problem size");

    if (max_level > 0) 
     if (regrid_int <= 0)
      amrex::Error("checkinput: regrid_int must be positive if max_level>0");

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
       std::cout << "Successfully read inputs file ... " << '\n';
} // end subroutine AMR_checkInput ()

// called from main.cpp after having called:
// Amr* amrptr = new Amr(); which calls
// AmrMesh()
// InitAmr() 
//  they call,
// geom[i].define(index_domain)  i=0...max_level and
// m_gdb.reset(new AmrParGDB(this));
void
Amr::init (Real strt_time,
           Real stop_time)
{
    if (!restart_file.empty() && restart_file != "init")
    {
        restart(restart_file);
    }
    else
    {
        initialInit(strt_time,stop_time);
        checkPoint();
        if (plot_int > 0 || plot_per > 0) {
         int do_plot=1;
         int do_slice=((slice_int>0) ? 1 : 0);
         int SDC_outer_sweeps=0;
         int slab_step=Time_blockingFactor()-1;
         writePlotFile(
          do_plot,do_slice,
          SDC_outer_sweeps,slab_step);
        }
    }
}

void
Amr::initialInit (Real strt_time,
                  Real stop_time)
{
 AMR_checkInput();

 finest_level = 0;
 fort_override_finest_level(&finest_level);

 cumtime = strt_time;

 // define level=0 
 defBaseLevel(strt_time);
   
 // timestep for just level=0 (finest_level=0) 
 amr_level[0]->computeInitialDt(finest_level,
                               dt_AMR,
                               stop_time);

 if (max_level > 0)
  bldFineLevels(strt_time);

 for (int lev = 0; lev <= finest_level; lev++)
  amr_level[lev]->setTimeLevel(strt_time,dt_AMR);

 int initialInit_flag=1;

 for (int lev = 0; lev <= finest_level; lev++)
  amr_level[lev]->post_regrid(0,0,finest_level,initialInit_flag,strt_time);

 // recomputes the timestep on all levels and updates statedata and dt_AMR
 // ns.post_init calls: post_init_state, computeInitialDt, and
 // sum_integrated_quantities.
 for (int lev = 0; lev <= finest_level; lev++)
  amr_level[lev]->post_init(stop_time);

 for (int lev = 0; lev <= finest_level; lev++) {
  level_count[lev] = 0;
  level_steps[lev] = 0;
  level_cells_advanced[lev] = 0.0;
 }

 if (ParallelDescriptor::IOProcessor()) {
  if (verbose > 1) {
   std::cout << "INITIAL GRIDS \n";
   printGridInfo(std::cout,0,finest_level);
  } else if (verbose > 0) { 
   std::cout << "INITIAL GRIDS \n";
   printGridSummary(std::cout,0,finest_level);
  }
 }

 if (record_grid_info && ParallelDescriptor::IOProcessor()) {
  gridlog << "INITIAL GRIDS \n";
  printGridInfo(gridlog,0,finest_level);
 }

} // subroutine initialInit

void
Amr::restart (const std::string& filename)
{
    double dRestartTime0 = ParallelDescriptor::second();

    VisMF::SetMFFileInStreams(mffile_nstreams);

    int i;

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
     std::cout << "restarting calculation from file: " << 
       filename << std::endl;

    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "RESTART from file = " << filename << '\n';

    //
    // Start calculation from given restart file.
    //
    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "RESTART from file = " << filename << '\n';
    //
    // Open the checkpoint header file for reading.
    //
    std::string File = filename;
    std::string FullPathName=filename;

    File += '/';
    File += "Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);
    //
    // Read global data.
    //
    // Attempt to differentiate between old and new CheckPointFiles.
    //
    int         spdim;
    std::string first_line;

    std::getline(is,first_line);

    if (first_line == CheckPointVersion) {  // new_check_point_format==true
        is >> spdim;
    } else {  // new_check_point_format==false
        spdim = atoi(first_line.c_str());
    }

    if (spdim != AMREX_SPACEDIM) {
     std::cerr << "Amr::restart(): bad spacedim = " << spdim << '\n';
     amrex::Abort();
    }

    is >> cumtime;
    int mx_lev;
    is >> mx_lev;
    is >> finest_level;
    int old_finest_level=finest_level;

    fort_override_finest_level(&finest_level);

    Vector<Box> inputs_domain(max_level+1);
    for (int lev = 0; lev <= max_level; lev++) {
     Box bx(geom[lev].Domain().smallEnd(),geom[lev].Domain().bigEnd());
     inputs_domain[lev] = bx;
    }

// SUSSMAN KLUGE RESTART

    int local_num_materials=0;
    is >> local_num_materials;
    if (local_num_materials<=0)
     amrex::Error("local_num_materials invalid in restart");

    AMR_max_phase_change_rate.resize(AMREX_SPACEDIM);
    AMR_min_phase_change_rate.resize(AMREX_SPACEDIM);
    for (int j=0;j<AMREX_SPACEDIM;j++) {
     is >> AMR_max_phase_change_rate[j];
     is >> AMR_min_phase_change_rate[j];
    }

    is >> AMR_volume_history_recorded;
    AMR_volume_history.resize(local_num_materials);
    for (int j=0;j<local_num_materials;j++) {
     is >> AMR_volume_history[j];
     if (AMR_volume_history[j]<0.0)
      amrex::Error("cannot have negative volume_history");
    }

    if (local_num_materials==global_AMR_num_materials) {
     // do nothing
    } else
     amrex::Error("local_num_materials!=global_AMR_num_materials");

    int checkpoint_recalesce_data=0;
    for (int im=0;im<global_AMR_num_materials;im++) {
     if (recalesce_flag[im]!=0)
      checkpoint_recalesce_data=1;
    }

    if (checkpoint_recalesce_data==0) {
     // do nothing
    } else if (checkpoint_recalesce_data==1) {
     int array_size;

     is >> array_size;
     recalesce_state_old.resize(array_size);
     for (int j=0;j<array_size;j++)
      is >> recalesce_state_old[j];

     is >> array_size;
     recalesce_state_new.resize(array_size);
     for (int j=0;j<array_size;j++)
      is >> recalesce_state_new[j];

    } else
     amrex::Error("checkpoint_recalesce_data invalid");

// END SUSSMAN KLUGE RESTART

    is >> dt_AMR;

     //restarting with greater or equal to the previous number of levels.
    if (max_level >= mx_lev) {

       for (i = 0; i <= mx_lev; i++) is >> geom[i];
       for (i = 0; i <= mx_lev; i++) is >> level_steps[i];

        // level_cells_advanced not checkpointed
       for (i = 0; i <= mx_lev; i++) {
        level_cells_advanced[i]=0.0;
       }
       for (i = 0; i <= mx_lev; i++) is >> level_count[i];

       //
       // Set bndry conditions.
       //
       if (max_level > mx_lev) {

        for (i = mx_lev+1; i <= max_level; i++) {
         level_steps[i] = level_steps[i-1];
         level_count[i] = 0;
         level_cells_advanced[i] = 0.0;
        }

       }

       if (regrid_on_restart and max_level > 0)
           level_count[0] = regrid_int;

       AMR_checkInput();
       //
       // Read levels.
       //
       for (int lev(0); lev <= finest_level; ++lev)
       {
           amr_level[lev].reset((*levelbld)());
            // internal to amr_level -> restart are the commands:
            // parent->SetBoxArray(level, grids);
            // parent->SetDistributionMap(level, dmap);
	    // new_finest_level = finest_level
           amr_level[lev]->restart(*this, is,old_finest_level,finest_level);
           this->SetBoxArray(lev, amr_level[lev]->boxArray());
           this->SetDistributionMap(lev, amr_level[lev]->DistributionMap());
       }
       //
       // Build any additional data structures.
       //

       for (int lev(0); lev <= finest_level; lev++)
           amr_level[lev]->post_restart();

     // restarting with a coarser mesh (less levels)
    } else if ((max_level>=0)&&(max_level<mx_lev)) {

       if (ParallelDescriptor::IOProcessor())
          amrex::Warning("Amr::restart(): max_level is lower than before");

        //just in case max_level reduced to be below the restarted
        //max_level which is  "mx_lev"
       int new_finest_level = std::min(max_level,finest_level);

       finest_level = new_finest_level;
       fort_override_finest_level(&finest_level);
 
       // These are just used to hold the extra stuff we have to read in.
       Geometry   geom_dummy;
       int         int_dummy;
       IntVect intvect_dummy;

       for (i = 0          ; i <= max_level; i++) is >> geom[i];
         // discard the previous levels in which level>max_level.
       for (i = max_level+1; i <= mx_lev   ; i++) is >> geom_dummy;

       for (i = max_level; i <  mx_lev   ; i++) is >> intvect_dummy;

       for (i = 0          ; i <= max_level; i++) is >> level_steps[i];
       for (i = max_level+1; i <= mx_lev   ; i++) is >> int_dummy;

        // level_cells_advanced not checkpointed
       for (i = 0          ; i <= max_level; i++) {
        level_cells_advanced[i]=0.0;
       }
       for (i = max_level+1; i <= mx_lev   ; i++) {
        //level_cells_advanced[i]=0.0;
       }

       for (i = 0          ; i <= max_level; i++) is >> level_count[i];
       for (i = max_level+1; i <= mx_lev   ; i++) is >> int_dummy;

       if (regrid_on_restart and max_level > 0)
           level_count[0] = regrid_int;

       AMR_checkInput();

       //
       // Read levels.
       //
       int lev;
       for (lev = 0; lev <= new_finest_level; lev++)
       {
           amr_level[lev].reset((*levelbld)());
            // internal to amr_level -> restart are the commands:
            // parent->SetBoxArray(level, grids);
            // parent->SetDistributionMap(level, dmap);
	    // new_finest_level = finest_level
           amr_level[lev]->restart(*this, is,old_finest_level,finest_level);
           this->SetBoxArray(lev, amr_level[lev]->boxArray());
           this->SetDistributionMap(lev, amr_level[lev]->DistributionMap());
       }
       //
       // Build any additional data structures.
       //
       for (lev = 0; lev <= new_finest_level; lev++)
           amr_level[lev]->post_restart();

    } else
     amrex::Error("max_level or mx_lev invalid");

    for (int lev = 0; lev <= finest_level; lev++)
    {
       Box restart_domain(geom[lev].Domain());
       if (! (inputs_domain[lev] == restart_domain) )
       {
          if (ParallelDescriptor::IOProcessor())
          {
             std::cout << "Problem at level " << lev << '\n';
             std::cout << "Domain according to     inputs file is " <<  inputs_domain[lev] << '\n';
             std::cout << "Domain according to checkpoint file is " << restart_domain      << '\n';
             std::cout << "Amr::restart() failed -- box from inputs file does not equal box from restart file" << std::endl;
          }
          amrex::Abort();
       }
    }


    if (verbose > 0)
    {
        double dRestartTime = ParallelDescriptor::second() - dRestartTime0;

        ParallelDescriptor::ReduceRealMax(dRestartTime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "Restart time = " << dRestartTime << " seconds." << '\n';
        }
    }
}

void
Amr::checkPoint ()
{
 if (!checkpoint_files_output) return;

 VisMF::SetNOutFiles(checkpoint_nfiles);
 //
 // In checkpoint files always write out FABs in NATIVE format.
 //
 FABio::Format thePrevFormat = FArrayBox::getFormat();

 FArrayBox::setFormat(FABio::FAB_NATIVE);

 double dCheckPointTime0 = ParallelDescriptor::second();

 const std::string ckfile = amrex::Concatenate(check_file_root,level_steps[0],file_name_digits);

 std::string FullPathName=ckfile;

 if (verbose > 0 && ParallelDescriptor::IOProcessor())
     std::cout << "CHECKPOINT: file = " << ckfile << std::endl;

 if (record_run_info && ParallelDescriptor::IOProcessor())
     runlog << "CHECKPOINT: file = " << ckfile << '\n';
 //
 // Only the I/O processor makes the directory if it doesn't already exist.
 //
 if (ParallelDescriptor::IOProcessor())
     if (!amrex::UtilCreateDirectory(ckfile, 0755))
         amrex::CreateDirectoryFailed(ckfile);
 //
 // Force other processors to wait till directory is built.
 //
 ParallelDescriptor::Barrier();

 std::string HeaderFileName = ckfile + "/Header";

 VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

 std::ofstream HeaderFile;

 HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

 int old_prec = 0, i;

 if (ParallelDescriptor::IOProcessor()) {
  //
  // Only the IOProcessor() writes to the header file.
  //
  HeaderFile.open(HeaderFileName.c_str(), std::ios::out|std::ios::trunc|std::ios::binary);

  if (!HeaderFile.good())
      amrex::FileOpenFailed(HeaderFileName);

  old_prec = HeaderFile.precision(15);

  HeaderFile << CheckPointVersion << '\n'
             << AMREX_SPACEDIM       << '\n'
             << cumtime           << '\n'
             << max_level         << '\n'
             << finest_level      << '\n';

// SUSSMAN KLUGE CHECKPOINT

  int local_num_materials=AMR_volume_history.size();
  if (local_num_materials<=0)
   amrex::Error("local_num_materials invalid in checkpoint");
  HeaderFile << local_num_materials << '\n';

  if (local_num_materials==global_AMR_num_materials) {
   // do nothing
  } else
   amrex::Error("local_num_materials!=global_AMR_num_materials");

  for (int j=0;j<AMREX_SPACEDIM;j++) {
   HeaderFile << AMR_max_phase_change_rate[j] << '\n';
   HeaderFile << AMR_min_phase_change_rate[j] << '\n';
  }

  HeaderFile << AMR_volume_history_recorded << '\n';
  for (int j=0;j<local_num_materials;j++) {
   HeaderFile << AMR_volume_history[j] << '\n';
   if (AMR_volume_history[j]<0.0) 
    amrex::Error("AMR_volume_history cannot be negative");
  }

  int checkpoint_recalesce_data=0;
  for (int im=0;im<global_AMR_num_materials;im++) {
   if (recalesce_flag[im]!=0)
    checkpoint_recalesce_data=1;
  }

  if (checkpoint_recalesce_data==0) {
   // do nothing
  } else if (checkpoint_recalesce_data==1) {

   HeaderFile << recalesce_state_old.size() << '\n';
   for (int j=0;j<recalesce_state_old.size();j++)
    HeaderFile << recalesce_state_old[j] << '\n';
   HeaderFile << recalesce_state_new.size() << '\n';
   for (int j=0;j<recalesce_state_new.size();j++)
    HeaderFile << recalesce_state_new[j] << '\n';

  } else
   amrex::Error("checkpoint_recalesce_data invalid");

// END SUSSMAN KLUGE CHECKPOINT

  HeaderFile << dt_AMR << ' ' << '\n';

  //
  // Write out problem domain.
  //
  for (i = 0; i <= max_level; i++) HeaderFile << geom[i]        << ' ';
  HeaderFile << '\n';

  for (i = 0; i <= max_level; i++) HeaderFile << level_steps[i] << ' ';
  HeaderFile << '\n';

    // level_cells_advanced is not checkpointed.

  for (i = 0; i <= max_level; i++) HeaderFile << level_count[i] << ' ';
  HeaderFile << '\n';

 }  // if IOProcessor

  // dir="ckfile" (directory)
  // os=HeaderFile 
 for (i = 0; i <= finest_level; i++)
  amr_level[i]->checkPoint(ckfile, HeaderFile);

 if (ParallelDescriptor::IOProcessor()) {
  HeaderFile.precision(old_prec);

  if (!HeaderFile.good())
   amrex::Error("Amr::checkpoint() failed");
 }

 //
 // Don't forget to reset FAB format.
 //
 FArrayBox::setFormat(thePrevFormat);

 if (verbose > 0) {
  double dCheckPointTime = ParallelDescriptor::second() - dCheckPointTime0;

  ParallelDescriptor::ReduceRealMax(dCheckPointTime,
    ParallelDescriptor::IOProcessorNumber());

  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "checkPoint() time = " << dCheckPointTime << 
     " secs." << '\n';
  }
 }
 ParallelDescriptor::Barrier();
}


void
Amr::regrid_level_0_on_restart()
{


 if ((max_level==0)&&(regrid_on_restart==1)) {

  regrid_on_restart = 0;
  //
  // Coarsening before we split the grids ensures that each resulting
  // grid will have an even number of cells in each direction.
  //
  BoxArray lev0(amrex::coarsen(Geom(0).Domain(),2));
  //
  // Now split up into list of grids within max_grid_size[0] limit.
  //
  lev0.maxSize(Old_maxGridSize(0)/2);
  //
  // Now refine these boxes back to level 0.
  //
  lev0.refine(2);
  
   //
   // Construct skeleton of new level.
   //
  DistributionMapping dm(lev0);
  AmrLevel* a = (*levelbld)(*this,0,Geom(0),lev0,dm,cumtime);
      
  a->init(*amr_level[0],lev0,dm);
  amr_level[0].reset(a);
      
  this->SetBoxArray(0, amr_level[0]->boxArray());
  this->SetDistributionMap(0, amr_level[0]->DistributionMap());

  // calls CopyNewToOld 
  // calls setTimeLevel(cumtime,dt_AMR) 
  int initialInit_flag=0;
  amr_level[0]->post_regrid(0,0,0,initialInit_flag,cumtime);
      
  if (ParallelDescriptor::IOProcessor()) {
   if (verbose > 1) {
    printGridInfo(amrex::OutStream(),0,finest_level);
   } else if (verbose > 0) {
    printGridSummary(amrex::OutStream(),0,finest_level);
   }
  }
      
  if (record_grid_info && ParallelDescriptor::IOProcessor())
   printGridInfo(gridlog,0,finest_level);

 } else {
  std::cout << "max_level not 0 or regrid_on_restart not 1\n";
  amrex::Error("regrid_level_0_on_restart(): invalid environment");
 }

} // end subroutine regrid_level_0_on_restart()


void
Amr::timeStep (Real time,
               Real stop_time)
{

 if (std::abs(time-cumtime)>1.0e-13)
  amrex::Error("time<>cumtime");

 if ((max_level==0)&&(regrid_on_restart==1)) {

  regrid_level_0_on_restart();

  if (record_grid_info && ParallelDescriptor::IOProcessor())
   printGridInfo(gridlog,0,finest_level);
 
 } else if ((max_level>0)||(regrid_on_restart==0)) {
  // do nothing
 } else
  amrex::Error("finest_level or regrid_on_restart invalid");


 int max_coarsest = std::min(finest_level, max_level-1);

 for (int level = 0; level <= max_coarsest; level++) {

  const int old_finest = finest_level;

  if (level_count[level] >= regrid_int) {

   regrid(level,time); // new levels might be created here.
 
   if (finest_level>old_finest+1)
    amrex::Error("cannot create more than one new level at a time");

   for (int k = level; k <= finest_level; k++)
    level_count[k] = 0;

  }
  int max_coarsest_new = std::min(finest_level, max_level-1);
  if (max_coarsest_new<max_coarsest)
   max_coarsest=max_coarsest_new;

 } // level=0...max_coarsest

 if ((plotfile_on_restart)&&(!(restart_file.empty()))) {
  plotfile_on_restart = 0;
  int do_plot=1;
  int do_slice=((slice_int>0) ? 1 : 0);
  int SDC_outer_sweeps=0;
  int slab_step=Time_blockingFactor()-1;
  writePlotFile(
   do_plot,do_slice,
   SDC_outer_sweeps,slab_step);
 }

 for (int level=0;level<=finest_level;level++) {

  if (verbose > 0 && ParallelDescriptor::IOProcessor()) {
   std::cout << "ADVANCE grids at level "
             << level
             << " with dt = "
             << dt_AMR
             << std::endl;
  }
  if (level==0) {
    // calls CopyNewToOldALL, setTimeLevel(time+dt,dt)
   Real dt_AMR_new = amr_level[level]->advance(time,dt_AMR);
   dt_AMR = dt_AMR_new;
  } else if (level>0) {
   // do nothing
  } else
   amrex::Error("level invalid");

  level_steps[level]++;
  level_count[level]++;
  level_cells_advanced[level]+=amr_level[level]->countCells();

  if (verbose > 0 && ParallelDescriptor::IOProcessor()) {
   std::cout << "Advanced "
             << amr_level[level]->countCells()
             << " cells at level "
             << level
             << std::endl;
  }

 } // level=0..finest_level

  // synchronization can be done here (e.g. average down)
 for (int level=finest_level;level>=0;level--) {
  amr_level[level]->post_timestep(stop_time);
 }

}  // subroutine timeStep

void Amr::recalesce_copy_new_to_old() {

 int recalesce_num_state=6;

 if ((recalesce_state_old.size()!=
      global_AMR_num_materials*recalesce_num_state)||
     (recalesce_state_new.size()!=
      global_AMR_num_materials*recalesce_num_state)) {
   amrex::Error("recalesce sizes incorrect");
 } else {

  for (int im=0;im<global_AMR_num_materials*recalesce_num_state;im++) {

   recalesce_state_old[im]=recalesce_state_new[im];

  } // im
 } // everything correct size?
} // recalesce_copy_new_to_old



void Amr::recalesce_copy_old_to_new() {

 int recalesce_num_state=6;

 if ((recalesce_state_old.size()!=global_AMR_num_materials*recalesce_num_state)||
     (recalesce_state_new.size()!=global_AMR_num_materials*recalesce_num_state)) {
   amrex::Error("recalesce sizes incorrect");
 } else {
  for (int im=0;im<global_AMR_num_materials*recalesce_num_state;im++) {

   recalesce_state_new[im]=recalesce_state_old[im];

  } // im
 } // everything correct size?
} // recalesce_copy_old_to_new



void Amr::recalesce_init() {

 int recalesce_num_state=6;

 recalesce_state_old.resize(recalesce_num_state*global_AMR_num_materials);
 recalesce_state_new.resize(recalesce_num_state*global_AMR_num_materials);

 for (int im=0;im<recalesce_num_state*global_AMR_num_materials;im++) {
  recalesce_state_old[im]=-1.0;
  recalesce_state_new[im]=-1.0;
 }

} // recalesce_init


void Amr::recalesce_get_state(Vector<Real>& recalesce_state_out) { 

 int recalesce_num_state=6;

 if (recalesce_state_out.size()!=recalesce_num_state*global_AMR_num_materials)
  amrex::Error("recalesce_state_out has incorrect size");
 if (recalesce_state_old.size()!=recalesce_num_state*global_AMR_num_materials)
  amrex::Error("recalesce_state_old has incorrect size");

 for (int im=0;im<recalesce_num_state*global_AMR_num_materials;im++)
  recalesce_state_out[im]=recalesce_state_old[im];

} // recalesce_get_state


void Amr::recalesce_put_state(Vector<Real>& recalesce_state_in) {

 int recalesce_num_state=6;

 if (recalesce_state_new.size()!=recalesce_num_state*global_AMR_num_materials)
  amrex::Error("recalesce_state_new has incorrect size");
 if (recalesce_state_in.size()!=recalesce_num_state*global_AMR_num_materials)
  amrex::Error("recalesce_state_in has incorrect size");

 for (int im=0;im<recalesce_num_state*global_AMR_num_materials;im++)
  recalesce_state_new[im]=recalesce_state_in[im];

} // recalesce_put_state


void
Amr::coarseTimeStep (Real stop_time)
{
    const double run_strt = ParallelDescriptor::second() ;

    //SUSSMAN
    //in: AMReX_FabArrayBase.cpp
    if (1==0) {
     FabArrayBase::flushTileArrayCache();
     FabArrayBase::flushFBCache();
     FabArrayBase::flushCPCache();
    }
    if (1==0) {
     if (ParallelDescriptor::IOProcessor()) {
      FabArrayBase::m_FA_stats.print();
      FabArrayBase::m_TAC_stats.print();
      FabArrayBase::m_FBC_stats.print();
      FabArrayBase::m_CPC_stats.print();
      FabArrayBase::m_FPinfo_stats.print();
      FabArrayBase::m_CFinfo_stats.print();
     }
    }


     // check dt on all the levels.
    if (level_steps[0] > 0) {
         // in AmrLevel.H: virtual void computeNewDt
         // NavierStokes::computeNewDt
        amr_level[0]->computeNewDt(finest_level,dt_AMR,stop_time);
    } else if (level_steps[0]==0) {
     // do nothing since initial dt already calculated 
     // NavierStokes::computeInitialDt
    } else
      amrex::Error("level_steps invalid");

     // dt_AMR might be modified within "timeStep"
    timeStep(cumtime,stop_time);

    cumtime += dt_AMR;

    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        double    run_stop = ParallelDescriptor::second() - run_strt;

        ParallelDescriptor::ReduceRealMax(run_stop,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "\nCoarse TimeStep time: " << run_stop << '\n' ;

    }

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "\nSTEP = "
                 << level_steps[0]
                  << " TIME = "
                  << cumtime
                  << " DT = "
                  << dt_AMR
                  << '\n'
                  << std::endl;
    }
    if (record_run_info && ParallelDescriptor::IOProcessor())
    {
        runlog << "STEP = "
               << level_steps[0]
               << " TIME = "
               << cumtime
               << " DT = "
               << dt_AMR
               << '\n';
    }
    if (record_run_info_terse && ParallelDescriptor::IOProcessor())
        runlog_terse << level_steps[0] << " " << 
        cumtime << " " << dt_AMR << '\n';

    int check_test = 0;
    if (check_per > 0.0)
    {
      const int num_per_old = cumtime / check_per;
      const int num_per_new = (cumtime+dt_AMR) / check_per;

      if (num_per_old != num_per_new)
	{
	 check_test = 1;
	}
    }

    int to_stop       = 0;    
    int to_checkpoint = 0;
    if (ParallelDescriptor::IOProcessor())
    {
        FILE *fp;
        if ((fp=fopen("dump_and_continue","r")) != 0)
        {
            remove("dump_and_continue");
            to_checkpoint = 1;
            fclose(fp);
        }
        else if ((fp=fopen("stop_run","r")) != 0)
        {
            remove("stop_run");
            to_stop = 1;
            fclose(fp);
        }
        else if ((fp=fopen("dump_and_stop","r")) != 0)
        {
            remove("dump_and_stop");
            to_checkpoint = 1;
            to_stop = 1;
            fclose(fp);
        }
    }

    //SUSSMAN
    ParallelDescriptor::Bcast(&to_checkpoint, 1, 
	ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&to_stop,       1, 
        ParallelDescriptor::IOProcessorNumber());

    ParallelDescriptor::Barrier();

    //SUSSMAN
    int check_int_trigger=0;
    if (check_int>0) {
     if (level_steps[0] % check_int == 0)
      check_int_trigger=1;
    } else if (check_int==0) {
     // do nothing
    } else if (check_int==-1) {
     // do nothing
    } else
     amrex::Error("check_int invalid");

    if ((check_int_trigger==1) || 
        (check_test == 1) || 
        (to_checkpoint==1)) {
     last_checkpoint = level_steps[0];
     checkPoint();
    } else if ((check_int_trigger==0)&&
	       (check_test==0)&&
	       (to_checkpoint==0)) {
     // do nothing
    } else
     amrex::Error("check_int_trigger,check_test,or to_checkpoint bad");

    int plot_test = 0;
    if (plot_per > 0.0)
    {
#ifdef BL_USE_OLDPLOTPER
      const int num_per_old = (cumtime-dt_AMR) / plot_per;
      const int num_per_new = (cumtime       ) / plot_per;

      if (num_per_old != num_per_new)
#else
      Long lorig=(Long) cumtime/plot_per;
      Real rorig=cumtime/plot_per;
      Real rR=(Real) lorig;
      rR=rorig-rR;
      if (rR < (dt_AMR*0.001))
#endif
	{
	  plot_test = 1;
	}
    }

    int do_plot=0;
    if (((plot_int > 0)&&(level_steps[0] % plot_int == 0)) || 
        (plot_test == 1)|| 
        (to_checkpoint))
     do_plot=1;
    int do_slice=0;
    if ((slice_int > 0)&&(level_steps[0] % slice_int == 0)) 
     do_slice=1;
    
    if ((do_plot==1)||(do_slice==1)) { 
     if (do_plot==1)
      last_plotfile = level_steps[0];
     int SDC_outer_sweeps=0;
     int slab_step=Time_blockingFactor()-1;
     writePlotFile(
      do_plot,do_slice,
      SDC_outer_sweeps,slab_step);
    }

    if (to_stop)
    {
        ParallelDescriptor::Barrier();
        if (to_checkpoint)
        {
            amrex::Abort("Stopped by user w/ checkpoint");
        }
        else
        {
            amrex::Abort("Stopped by user w/o checkpoint");
        }
    }
} // end subroutine coarseTimeStep

void
Amr::defBaseLevel (Real strt_time)
{

    if (finest_level==0) {
     // do nothing
    } else
     amrex::Error("finest_level invalid");

    const Box& domain = Geom(0).Domain();
    IntVect d_length  = domain.size();
    const int* d_len  = d_length.getVect();

    for (int idir = 0; idir < AMREX_SPACEDIM; idir++)
     if (d_len[idir]%2 != 0)
      amrex::Error("defBaseLevel: must have even number of cells");

    BoxArray lev0=MakeBaseGrids();

    this->SetBoxArray(0, lev0);

     // SUSSMAN
    int nprocs=ParallelDescriptor::NProcs();
    DistributionMapping dm(lev0,nprocs);
    this->SetDistributionMap(0, dm);

    //
    // Now build level 0 grids.
    //
    amr_level[0].reset((*levelbld)(*this,0,Geom(0),grids[0],dmap[0],strt_time));

    amr_level[0]->initData();
} // subroutine defBaseLevel

void
Amr::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
    amr_level[lev]->errorEst(tags,TagBox::CLEAR,TagBox::SET,time, 
		    n_error_buf[lev][0], ngrow);
}

BoxArray
Amr::GetAreaNotToTag (int lev)
{
    return BoxArray(amr_level[lev]->getAreaNotToTag());
}

void
Amr::ManualTagsPlacement (int lev, TagBoxArray& tags, const Vector<IntVect>& bf_lev)
{
    amr_level[lev]->manual_tags_placement(tags, bf_lev);
}


// called from timeStep and bldFineLevels.
// bldFineLevels is called from initialInit.
//  (note: defBaseLevel is also called from initialInit)
// levelbld is called from defBaseLevel, timeStep, regrid, bldFineLevels
void
Amr::regrid (int  lbase,
             Real time,
             bool initial)
{

 BL_PROFILE("Amr::regrid()");

 if (std::abs(time-cumtime)>1.0e-13)
  amrex::Error("time<>cumtime in regrid");

 if (verbose > 0 && ParallelDescriptor::IOProcessor())
  std::cout << "REGRID: at level lbase = " << lbase << std::endl;

 if (record_run_info && ParallelDescriptor::IOProcessor())
  runlog << "REGRID: at level lbase = " << lbase << '\n';

 int max_coarsest=std::min(finest_level,max_level-1);
 if (lbase>max_coarsest)
  amrex::Error("cannot have lbase>max_coarsest");

 int new_finest;
 Vector<BoxArray> new_grid_places(max_level+1);
 Vector<DistributionMapping> new_dmap(max_level+1);

 grid_places(lbase,time,new_finest,new_grid_places);

 if (new_finest>finest_level+1)
  amrex::Error("cannot create more than one new level at a time");

 int regrid_level_zero=0;
 if (lbase==0) {
  if (new_grid_places[0] != amr_level[0]->boxArray())
   regrid_level_zero=1;
 }

  //
  // Reclaim all remaining storage for levels > new_finest.
  //
 for (int lev = new_finest+1; lev <= finest_level; lev++) {
  amr_level[lev].reset();
  this->ClearBoxArray(lev);
  this->ClearDistributionMap(lev);
 }

 const int start = ((regrid_level_zero==1) ? 0 : lbase+1);

 for (int lev = start, End = std::min(finest_level,new_finest); 
      lev <= End; lev++) {
  if (new_grid_places[lev] == amr_level[lev]->boxArray()) {
   new_grid_places[lev] = amr_level[lev]->boxArray();  // to avoid duplicates
   new_dmap[lev] = amr_level[lev]->DistributionMap();
  } else {
   // do nothing
  }
 }  // lev=start ... min(finest_level,new_finest)

 finest_level = new_finest;
 fort_override_finest_level(&finest_level);

 for (int lev = start; lev <= new_finest; lev++) {

  if (new_grid_places[lev].size()<1) {
   std::cout << "initial= " << initial << '\n';
   std::cout << "start= " << start << '\n';
   std::cout << "lbase= " << lbase << '\n';
   std::cout << "finest_level= " << finest_level << '\n';
   std::cout << "new_finest= " << new_finest << '\n';
   std::cout << "amr_level[0]->boxArray \n";
   std::cout << amr_level[0]->boxArray() << '\n';
   std::cout << "lev= " << lev << '\n';
  }

  if (new_dmap[lev].empty()) {
   new_dmap[lev].define(new_grid_places[lev]);
  }

  AmrLevel* a = (*levelbld)(*this,lev,geom[lev],
    new_grid_places[lev],new_dmap[lev],cumtime);

  if (initial) {
   //
   // We're being called on startup from bldFineLevels().
   // NOTE: The initData function may use a filPatch, and so needs to
   //       be officially inserted into the hierarchy prior to the call.
   //
   amr_level[lev].reset(a);

   this->SetBoxArray(lev, amr_level[lev]->boxArray());
   this->SetDistributionMap(lev, amr_level[lev]->DistributionMap());

   amr_level[lev]->initData();

  } else if (amr_level[lev]) {

    // amr_level[lev-1] should already be init.
   a->init(*amr_level[lev],new_grid_places[lev],new_dmap[lev]);

   amr_level[lev].reset(a);

   this->SetBoxArray(lev, amr_level[lev]->boxArray());
   this->SetDistributionMap(lev, amr_level[lev]->DistributionMap());

  } else {

   if (lev>0) {
    a->init(new_grid_places[lev],new_dmap[lev]);
    amr_level[lev].reset(a);
    this->SetBoxArray(lev, amr_level[lev]->boxArray());
    this->SetDistributionMap(lev, amr_level[lev]->DistributionMap());
   } else
    amrex::Error("expecting lev>0 if amr_level[lev]==FALSE");

  }

 } // lev=start..new_finest

 int initialInit_flag=0;

 for (int lev = start; lev <= new_finest; lev++) {
  amr_level[lev]->post_regrid(lbase,start,new_finest,
    initialInit_flag,time);
 }

 if (record_run_info && ParallelDescriptor::IOProcessor()) {
  runlog << "REGRID: at level lbase = " << lbase << '\n';
  printGridInfo(runlog,start,finest_level);
 }
 if (record_grid_info && ParallelDescriptor::IOProcessor()) {
  if (lbase == 0)
   gridlog << "STEP = " << level_steps[0] << ' ';

  gridlog << "TIME = "
          << time
          << " : REGRID  with lbase = "
          << lbase
          << '\n';

  printGridInfo(gridlog,start,finest_level);
 }
 if (verbose > 0 && ParallelDescriptor::IOProcessor()) {
  if (lbase == 0)
   std::cout << "STEP = " << level_steps[0] << ' ';

  std::cout << "TIME = "
            << time
            << " : REGRID  with lbase = "
            << lbase
            << std::endl;

  if (verbose > 1) {
   printGridInfo(std::cout,start,finest_level);
  } else {
   printGridSummary(std::cout,start,finest_level);
  }
 }
} // end subroutine regrid

void 
Amr::print_cells_advanced() {

 for (int lev = 0; lev <= max_level; lev++) {
  std::cout << "LEVEL= " << lev << " CELLS ADVANCED= " <<
   level_cells_advanced[lev] << '\n';
 } 
}

void
Amr::printGridInfo (std::ostream& os,
                    int           min_lev,
                    int           max_lev)
{
    for (int lev = min_lev; lev <= max_lev; lev++)
    {
        const BoxArray&           bs      = amr_level[lev]->boxArray();
        int                       numgrid = bs.size();
        Long                      ncells  = amr_level[lev]->countCells();
        double                    ntot    = geom[lev].Domain().d_numPts();
        Real                      frac    = 100.0*(Real(ncells) / ntot);
        int dummy_slab_step=0;
        const DistributionMapping& map = 
         amr_level[lev]->get_new_data(0,dummy_slab_step+1).DistributionMap();

        os << "  Level "
           << lev
           << "   "
           << numgrid
           << " grids  "
           << ncells
           << " cells  "
           << frac
           << " % of domain"
           << '\n';


        for (int k = 0; k < numgrid; k++)
        {
            const Box& b = bs[k];

            os << ' ' << lev << ": " << b << "   ";
                
            for (int i = 0; i < AMREX_SPACEDIM; i++)
                os << b.length(i) << ' ';

            os << ":: " << map[k] << '\n';
        }
    }

    os << std::endl; // Make sure we flush!
}

void
Amr::printGridSummary (std::ostream& os,
                       int           min_lev,
                       int           max_lev)
{
    for (int lev = min_lev; lev <= max_lev; lev++)
    {
        const BoxArray&           bs      = amr_level[lev]->boxArray();
        int                       numgrid = bs.size();
        Long                      ncells  = amr_level[lev]->countCells();
        double                    ntot    = geom[lev].Domain().d_numPts();
        Real                      frac    = 100.0*(Real(ncells) / ntot);

        os << "  Level "
           << lev
           << "   "
           << numgrid
           << " grids  "
           << ncells
           << " cells  "
           << frac
           << " % of domain"
           << '\n';
    }

    os << std::endl; // Make sure we flush!
}

// new_finest cannot be greater than finest_level+1
void
Amr::grid_places (int              lbase,
                  Real time,
                  int&             new_finest,
                  Vector<BoxArray>& new_grids)
{

    BL_PROFILE("Amr::grid_places()");

    const Real strttime = amrex::second();

    if (lbase == 0)
    {
	new_grids[0] = MakeBaseGrids();
    }

    if ( time == 0. && !initial_grids_file.empty() && !use_fixed_coarse_grids)
    {
        new_finest = std::min(max_level,(finest_level+1));
        new_finest = std::min<int>(new_finest,initial_ba.size());

        for (int lev = 1; lev <= new_finest; lev++)
        {
            BoxList bl;
            int ngrid = initial_ba[lev-1].size();
            for (int i = 0; i < ngrid; i++)
            {
                Box bx(initial_ba[lev-1][i]);
                if (lev > lbase)
                    bl.push_back(bx);
            }
            if (lev > lbase)
                new_grids[lev].define(bl);
        }
        return;
    }

    // Use grids in initial_grids_file as fixed coarse grids.
    if ( ! initial_grids_file.empty() && use_fixed_coarse_grids)
    {
        new_finest = std::min(max_level,(finest_level+1));
        new_finest = std::min<int>(new_finest,initial_ba.size());

        for (int lev = lbase+1; lev <= new_finest; lev++)
        {
            BoxList bl;
            int ngrid = initial_ba[lev-1].size();
            for (int i = 0; i < ngrid; i++)
            {
                Box bx(initial_ba[lev-1][i]);

                if (lev > lbase)
                    bl.push_back(bx);

            }
            if (lev > lbase)
                new_grids[lev].define(bl);
            new_grids[lev].maxSize(Old_maxGridSize(lev));
        }
    }
    else if ( !regrid_grids_file.empty() )  // Use grids in regrid_grids_file 
    {
        new_finest = std::min(max_level,(finest_level+1));
        new_finest = std::min<int>(new_finest,regrid_ba.size());
        for (int lev = 1; lev <= new_finest; lev++)
        {
            BoxList bl;
            int ngrid = regrid_ba[lev-1].size();
            for (int i = 0; i < ngrid; i++)
            {
                Box bx(regrid_ba[lev-1][i]);
                if (lev > lbase)
                    bl.push_back(bx);
            }
            if (lev > lbase)
                new_grids[lev].define(bl);
        }
        return;
    }

    MakeNewGrids(lbase, time, new_finest, new_grids);

    if (verbose > 0)
    {
        Real stoptime = amrex::second() - strttime;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
	amrex::Print() << "grid_places() time: " << stoptime << " new finest: " << new_finest<< '\n';
#ifdef BL_LAZY
	});
#endif
    }


}  // subroutine grid_places

// called from initialInit
void
Amr::bldFineLevels (Real strt_time)
{
 BL_PROFILE("Amr::bldFineLevels()");
 if (max_level<=0)
  amrex::Error("max_level invalid in bldFineLevels");

 finest_level = 0;
 fort_override_finest_level(&finest_level);

 Vector<BoxArray> new_grid_places(max_level+1);

 int new_finest=0;
 int previous_new_finest=0;
 int num_grid_places=0;
 int grid_places_done=0;

 while (grid_places_done==0) {

  grid_places(finest_level,strt_time,new_finest,new_grid_places);

  if (new_finest>finest_level) {

   if (new_finest>finest_level+1)
    amrex::Error("cannot create more than one new level at a time");

   finest_level = new_finest;
   fort_override_finest_level(&finest_level);

   if (new_grid_places[new_finest].size()<1) {
    std::cout << "IN bldFineLevels\n";
    std::cout << "finest_level= " << finest_level << '\n';
    std::cout << "new_finest= " << new_finest << '\n';
    std::cout << "new_grid_places[0] \n";
    std::cout << new_grid_places[0] << '\n';
    amrex::Error("new_grid_places[new_finest] invalid");
   }

    // SUSSMAN
   int nprocs=ParallelDescriptor::NProcs();
   DistributionMapping new_dm(new_grid_places[new_finest],nprocs);

    // see the constructor in AmrLevel.cpp:
    // AmrLevel::AmrLevel ( ....  )
   AmrLevel* a_level = (*levelbld)(*this,
                         new_finest,
                         geom[new_finest],
                         new_grid_places[new_finest],
                         new_dm,
                         strt_time);

   amr_level[new_finest].reset(a_level);
   this->SetBoxArray(new_finest, new_grid_places[new_finest]);
   this->SetDistributionMap(new_finest, new_dm);

   amr_level[new_finest]->initData();

  } else if ((new_finest>=0)&&(new_finest<=finest_level)) {
   // do nothing
  } else
   amrex::Error("new_finest invalid");

  if (finest_level>max_level)
   amrex::Error("finest_level is corrupt");
  
  if (finest_level==max_level)
   grid_places_done=1;

  if (num_grid_places>0) {
   if (previous_new_finest==new_finest)
    grid_places_done=1;
  }
  num_grid_places++;
  previous_new_finest=new_finest;

 }  // while (grid_places_done==0)

 if (verbose > 0 && ParallelDescriptor::IOProcessor())
  std::cout << "num_grid_places= " << num_grid_places << '\n';

 bool grids_the_same;

 const int MaxCnt = 4;

 int count = 0;

 if (max_level<=0)
  amrex::Error("max_level invalid in bldFineLevels");

 do {
   for (int i = 0; i <= finest_level; i++)
    new_grid_places[i] = amr_level[i]->boxArray();

   regrid(0,strt_time,true);

   grids_the_same = true;

   for (int i = 0; i <= finest_level && grids_the_same; i++)
    if (!(new_grid_places[i] == amr_level[i]->boxArray()))
     grids_the_same = false;

   count++;
 } while (!grids_the_same && count < MaxCnt);


}  // subroutine bldFineLevels

} // namespace amrex

