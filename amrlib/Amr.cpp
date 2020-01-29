
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
#include <Cluster.H>
#include <LevelBld.H>
#include <AmrLevel.H>
#include <Amr.H>
#include <StateData.H>

#include <INTERP_F.H>

#ifdef BL_USE_ARRAYVIEW
#include <DatasetClient.H>
#endif

namespace amrex {

//
// Static class members.  
// Set defaults in Initialize()!!!
//
std::list<std::string> Amr::state_plot_vars;
bool                   Amr::first_plotfile;

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
    bool plot_files_output;
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
    Vector<int> AMR_FSI_flag;
    int  AMR_num_materials;

//}


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
    Amr::first_plotfile      = true;
    plot_nfiles              = 64;
    mffile_nstreams          = 1;
    probinit_natonce         = 512;
    plot_files_output        = true;
    checkpoint_nfiles        = 64;
    regrid_on_restart        = 0;
    use_efficient_regrid     = 0;
    AMR_num_materials = 0;
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
    Amr::state_plot_vars.clear();

    initialized = false;
}

bool Amr::Plot_Files_Output () { return plot_files_output; }

std::ostream&
Amr::DataLog (int i)
{
    return *datalog[i];
}

int 
Amr::AMR_recalesce_flag(int im) const {

 if ((im<1)||(im>AMR_num_materials))
  amrex::Error("im out of range");

 return recalesce_flag[im-1];

}


Vector<std::unique_ptr<AmrLevel> >&
Amr::getAmrLevels () noexcept
{
    return amr_level;
}


#ifdef AMREX_PARTICLES
void 
Amr::RedistributeParticles () 
{
    amr_level[0]->particle_redistribute(0,true);
}
#endif


Amr::Amr () {

     // init default values for some parameters.
    Initialize();
     // Geometry::Setup()
     // levelbld = getLevelBld();
     // ...
    InitAmr();

}


void
Amr::InitAmr () {

    AMR_volume_history_recorded=0;
    AMR_volume_history.resize(1);
    AMR_volume_history[0]=0.0;

    //
    // Setup Geometry from ParmParse file.
    // May be needed for variableSetup or even getLevelBld.
    //
    Geometry::Setup();
    //
    // Determine physics class.
    //
    std::fflush(NULL);
    if (1==1) {
     std::cout << "levelbld = getLevelBld() on processor " <<
         ParallelDescriptor::MyProc() << "\n";
    }
    std::fflush(NULL);

    //LevelBld* levelbld
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
    max_level              = -1;
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

    ppns.get("num_materials",AMR_num_materials);
    if (AMR_num_materials<1)
     amrex::Error("AMR_num_materials invalid");

    recalesce_flag.resize(AMR_num_materials);
    for (int im=0;im<AMR_num_materials;im++) {
     recalesce_flag[im]=0;
    }
    ppns.queryarr("recalesce_flag",recalesce_flag,0,AMR_num_materials);
    for (int im=0;im<AMR_num_materials;im++) {
     if ((recalesce_flag[im]!=0)&&
         (recalesce_flag[im]!=1)&&
         (recalesce_flag[im]!=2))
      amrex::Error("recalesce_flag invalid");
    }

    AMR_FSI_flag.resize(AMR_num_materials);
    for (int im=0;im<AMR_num_materials;im++) {
     AMR_FSI_flag[im]=0;
    }
    ppns.queryarr("FSI_flag",AMR_FSI_flag,0,AMR_num_materials);
    for (int im=0;im<AMR_num_materials;im++) {
     if ((AMR_FSI_flag[im]!=0)&& 
         (AMR_FSI_flag[im]!=1)&&
         (AMR_FSI_flag[im]!=2)&&
         (AMR_FSI_flag[im]!=3)&&
         (AMR_FSI_flag[im]!=4)&&
	 (AMR_FSI_flag[im]!=5))
      amrex::Error("AMR_FSI_flag invalid in Amr.cpp");
    }

    pp.query("regrid_on_restart",regrid_on_restart);
    if ((regrid_on_restart!=0)&&(regrid_on_restart!=1))
     amrex::Error("regrid_on_restart invalid");

    pp.query("plotfile_on_restart",plotfile_on_restart);
    pp.query("checkpoint_on_restart",checkpoint_on_restart);

    pp.query("checkpoint_files_output", checkpoint_files_output);
    pp.query("plot_files_output", plot_files_output);

    plot_nfiles = ParallelDescriptor::NProcs();
    checkpoint_nfiles = ParallelDescriptor::NProcs();

    pp.query("refine_grid_layout", refine_grid_layout);

    pp.query("mffile_nstreams", mffile_nstreams);
    pp.query("probinit_natonce", probinit_natonce);

    probinit_natonce = std::max(1, std::min(ParallelDescriptor::NProcs(), probinit_natonce));

    pp.query("file_name_digits", file_name_digits);

    if (pp.contains("run_log"))
    {
        std::string log_file_name;
        pp.get("run_log",log_file_name);
        setRecordRunInfo(log_file_name);
    }
    if (pp.contains("run_log_terse"))
    {
        std::string log_file_name;
        pp.get("run_log_terse",log_file_name);
        setRecordRunInfoTerse(log_file_name);
    }
    if (pp.contains("grid_log"))
    {
        std::string grid_file_name;
        pp.get("grid_log",grid_file_name);
        setRecordGridInfo(grid_file_name);
    }

    if (pp.contains("data_log"))
    {
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
    //
    // Read max_level and alloc memory for container objects.
    //
    pp.get("max_level", max_level);
    int nlev     = max_level+1;

    geom.resize(nlev);
    dmap.resize(nlev);
    grids.resize(nlev);

    level_cells_advanced.resize(nlev);
    level_steps.resize(nlev);
    level_count.resize(nlev);
    blocking_factor.resize(nlev);
    space_blocking_factor.resize(nlev);
    max_grid_size.resize(nlev);
    n_error_buf.resize(nlev);
    amr_level.resize(nlev);
    //
    // Set bogus values.
    //

    dt_AMR=1.0;

    time_blocking_factor = 1;
    MAX_NUM_SLAB=33;
    slab_dt_type=0; // 0=SEM 1=evenly spaced

    for (i = 0; i < nlev; i++)
    {
        level_cells_advanced[i] = 0.0;
        level_steps[i] = 0;
        level_count[i] = 0;
        n_error_buf[i] = 1;
        blocking_factor[i] = 2;
        space_blocking_factor[i] = 1;
        max_grid_size[i] = (AMREX_SPACEDIM == 2) ? 128 : 32;
    }

    if (max_level > 0) 
    {
       regrid_int.resize(max_level);
       for (i = 0; i < max_level; i++)
           regrid_int[i]  = 0;
    }

    //
    // Read other amr specific values.
    //
    check_file_root = "chk";
    pp.query("check_file",check_file_root);

    check_int = -1;
    int got_check_int = pp.query("check_int",check_int);

    check_per = -1.0;
    int got_check_per = pp.query("check_per",check_per);

    if (got_check_int == 1 && got_check_per == 1)
    {
        amrex::Error("Must only specify amr.check_int OR amr.check_per");
    }

    plot_file_root = "plt";
    pp.query("plot_file",plot_file_root);

    plot_int = -1;
    int got_plot_int = pp.query("plot_int",plot_int);

    plot_per = -1.0;
    int got_plot_per = pp.query("plot_per",plot_per);

    if (got_plot_int == 1 && got_plot_per == 1)
    {
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

    pp.query("n_proper",n_proper);
    pp.query("grid_eff",grid_eff);
    pp.queryarr("n_error_buf",n_error_buf,0,max_level);
    //
    // Read in max_grid_size.  Use defaults if not explicitly defined.
    //
    int cnt = pp.countval("max_grid_size");

    if (cnt == 1)
    {
        //
        // Set all values to the single available value.
        //
        int the_max_grid_size = 0;

        pp.get("max_grid_size",the_max_grid_size);

        for (i = 0; i <= max_level; i++)
        {
            max_grid_size[i] = the_max_grid_size;
        }
    }
    else if (cnt > 1)
    {
        //
        // Otherwise we expect a vector of max_grid_size values.
        //
        pp.getarr("max_grid_size",max_grid_size,0,max_level+1);
    }
    //
    // Read in the blocking_factors.  Use defaults if not explicitly defined.
    //
    cnt = pp.countval("blocking_factor");

    if (cnt == 1)
    {
        //
        // Set all values to the single available value.
        //
        int the_blocking_factor = 0;

        pp.get("blocking_factor",the_blocking_factor);

        for (i = 0; i <= max_level; i++)
        {
            blocking_factor[i] = the_blocking_factor;
        }
    }
    else if (cnt > 1)
    {
        //
        // Otherwise we expect a vector of blocking factors.
        //
        pp.getarr("blocking_factor",blocking_factor,0,max_level+1);
    }

    pp.queryarr("space_blocking_factor",space_blocking_factor);
    pp.query("time_blocking_factor",time_blocking_factor);
    pp.query("MAX_NUM_SLAB",MAX_NUM_SLAB);
    if (time_blocking_factor+1>MAX_NUM_SLAB)
     amrex::Error("MAX_NUM_SLAB too small");

    pp.query("slab_dt_type",slab_dt_type);
    if ((slab_dt_type!=0)&&
        (slab_dt_type!=1))
     amrex::Error("slab_dt_type invalid");

    //
    // Read in the regrid interval if max_level > 0.
    //
    if (max_level > 0) 
    {
       int numvals = pp.countval("regrid_int");
       if (numvals == 1)
       {
           //
           // Set all values to the single available value.
           //
           int the_regrid_int = 0;
           pp.query("regrid_int",the_regrid_int);
           for (i = 0; i < max_level; i++)
           {
               regrid_int[i] = the_regrid_int;
           }
       }
       else if (numvals < max_level)
       {
           amrex::Error("You did not specify enough values of regrid_int");
       }
       else 
       {
           //
           // Otherwise we expect a vector of max_level values
           //
           pp.queryarr("regrid_int",regrid_int,0,max_level);
       }
    }
    //
    // Read computational domain and set geometry.
    //
    Vector<int> n_cell(AMREX_SPACEDIM);
    pp.getarr("n_cell",n_cell,0,AMREX_SPACEDIM);
    BL_ASSERT(n_cell.size() == AMREX_SPACEDIM);
    IntVect lo(IntVect::TheZeroVector()), hi(n_cell);
    hi -= IntVect::TheUnitVector();
    Box index_domain(lo,hi);
    for (i = 0; i <= max_level; i++)
    {
        geom[i].define(index_domain);
        if (i < max_level)
            index_domain.refine(2);
    }
    //
    // SUSSMAN:
    // Now check offset; CoordSys does not need it anymore though.
    //
    Real offset[AMREX_SPACEDIM];
    for (i = 0; i < AMREX_SPACEDIM; i++)
    {
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

} // subroutine InitAmr

void
Amr::SetDistributionMap (int lev, const DistributionMapping& dmap_in) noexcept
{
    if (dmap[lev] != dmap_in) dmap[lev] = dmap_in;
}

void
Amr::SetBoxArray (int lev, const BoxArray& ba_in) noexcept
{
    if (grids[lev] != ba_in) grids[lev] = ba_in;
}


void
Amr::ClearDistributionMap (int lev) noexcept
{
    dmap[lev] = DistributionMapping();
}

void
Amr::ClearBoxArray (int lev) noexcept
{
    grids[lev] = BoxArray();
}


bool
Amr::isStatePlotVar (const std::string& name)
{
    for (std::list<std::string>::const_iterator li = state_plot_vars.begin(), End = state_plot_vars.end();
         li != End;
         ++li)
    {
        if (*li == name)
            return true;
    }
    return false;
}

void
Amr::fillStatePlotVarList ()
{
    state_plot_vars.clear();
    const DescriptorList &desc_lst = AmrLevel::get_desc_lst();
    for (int typ(0); typ < desc_lst.size(); ++typ) {
        for (int comp(0); comp < desc_lst[typ].nComp(); ++comp) {
            if (desc_lst[typ].getType() == IndexType::TheCellType()) {
                state_plot_vars.push_back(desc_lst[typ].name(comp));
            }
        }
    }
}

void
Amr::clearStatePlotVarList ()
{
    state_plot_vars.clear();
}

void
Amr::addStatePlotVar (const std::string& name)
{
    if (!isStatePlotVar(name))
        state_plot_vars.push_back(name);
}

void
Amr::deleteStatePlotVar (const std::string& name)
{
    if (isStatePlotVar(name))
        state_plot_vars.remove(name);
}

Amr::~Amr ()
{
    if (level_steps[0] > last_checkpoint)
        checkPoint();

    if (level_steps[0] > last_plotfile) {
     int do_plot=1;
     int do_slice=((slice_int>0) ? 1 : 0);
     int SDC_outer_sweeps=0;
     int slab_step=Time_blockingFactor()-1;
     writePlotFile(plot_file_root,
      level_steps[0],do_plot,do_slice,
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

 const std::string debug_root("DEB");
 int do_plot=1;
 int do_slice=((slice_int>0) ? 1 : 0);
 writePlotFile(debug_root,
   num,do_plot,do_slice,
   SDC_outer_sweeps,slab_step);

}

void
Amr::writePlotFile (const std::string& root,
                    int                num,
                    int do_plot,int do_slice,
                    int SDC_outer_sweeps,int slab_step) {

 if (!Plot_Files_Output()) {
  return;
 }

 if (do_plot==1) {

  VisMF::SetNOutFiles(plot_nfiles);
  VisMF::Header::Version currentVersion(VisMF::GetHeaderVersion());
  VisMF::SetHeaderVersion(currentVersion);

  if (first_plotfile) {
     first_plotfile = false;
     amr_level[0]->setPlotVariables();
  }

  double dPlotFileTime0 = ParallelDescriptor::second();

  const std::string pltfile = amrex::Concatenate(root,num,file_name_digits);

  if (verbose>0) {
   amrex::Print() << "PLOTFILE: file = " << pltfile << '\n';
  }

  if (record_run_info && ParallelDescriptor::IOProcessor())
     runlog << "PLOTFILE: file = " << pltfile << '\n';

  int stream_max_tries=4;
  bool abort_on_stream_retry_failure=false;
  amrex::StreamRetry sretry(pltfile, abort_on_stream_retry_failure,
                             stream_max_tries);

  const std::string pltfileTemp(pltfile + ".temp");

  while(sretry.TryFileOutput()) {

   //
   //  if either the pltfile or pltfileTemp exists, rename them
   //  to move them out of the way.  then create pltfile
   //  with the temporary name, then rename it back when
   //  it is finished writing.  then stream retry can rename
   //  it to a bad suffix if there were stream errors.
   //

   if(precreateDirectories) {    // ---- make all directories at once
    amrex::UtilRenameDirectoryToOld(pltfile, false);      // dont call barrier
    if ((verbose > 1)||(1==1)) {
     amrex::Print() << "IOIOIOIO:  precreating directories for " << 
	     pltfileTemp << "\n";
    }
    amrex::PreBuildDirectorHierarchy(pltfileTemp, "Level_", 
	finest_level + 1, true);  // call barrier
   } else {
    amrex::UtilRenameDirectoryToOld(pltfile, false);     // dont call barrier
    amrex::UtilCreateCleanDirectory(pltfileTemp, true);  // call barrier
   }

   std::string HeaderFileName(pltfileTemp + "/Header");

   VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

   std::ofstream HeaderFile;

   HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

   int old_prec(0);

   if (ParallelDescriptor::IOProcessor()) {
     //
     // Only the IOProcessor() writes to the header file.
     //
    HeaderFile.open(HeaderFileName.c_str(), 
	std::ios::out | std::ios::trunc | std::ios::binary);
    if ( ! HeaderFile.good()) {
     amrex::FileOpenFailed(HeaderFileName);
    }
    old_prec = HeaderFile.precision(15);
   }

   for (int k(0); k <= finest_level; ++k) {
    amr_level[k]->writePlotFile(
	pltfileTemp, HeaderFile,
	do_plot,do_slice,
        SDC_outer_sweeps,slab_step);
   }

   if (ParallelDescriptor::IOProcessor()) {
    HeaderFile.precision(old_prec);
    if ( ! HeaderFile.good()) {
     amrex::Error("Amr::writePlotFile() failed");
    }
   }

    //last_plotfile = level_steps[0]; (set outside this routine)

   if (verbose > 0) {
    const int IOProc     = ParallelDescriptor::IOProcessorNumber();
    double dPlotFileTime = ParallelDescriptor::second() - dPlotFileTime0;

    ParallelDescriptor::ReduceRealMax(dPlotFileTime,IOProc);

    amrex::Print() << "Write plotfile time = " << dPlotFileTime << 
		"  seconds" << "\n\n";
   }
   ParallelDescriptor::Barrier("Amr::writePlotFile::end");

   if(ParallelDescriptor::IOProcessor()) {
     std::rename(pltfileTemp.c_str(), pltfile.c_str());
   }
   ParallelDescriptor::Barrier("Renaming temporary plotfile.");
   //
   // the plotfile file now has the regular name
   //
  }  // end while

  VisMF::SetHeaderVersion(currentVersion);

  BL_PROFILE_REGION_STOP("Amr::writePlotFile()");

 } else if ((do_plot==0)&&(do_slice==1)) {
  const std::string pltfile = amrex::Concatenate(root,num,file_name_digits);
  std::ofstream HeaderFile;
  for (int k = 0; k <= finest_level; k++)
   amr_level[k]->writePlotFile(pltfile,HeaderFile,
    do_plot,do_slice,
    SDC_outer_sweeps,slab_step);
 } else if ((do_plot==0)&&(do_slice==0)) {
  // do nothing
 } else
  amrex::Error("do_plot or do_slice invalid");
  
}  // subroutine writePlotFile


void
Amr::checkInput ()
{
    FabArrayBase::Initialize();

    if (max_level < 0)
        amrex::Error("checkInput: max_level not set");
    //
    // 1. Check that blocking_factor is a power of 2 and no smaller than 4.
    // 2. Check that blocking_factor[i+1]<=blocking_factor[i].
    // 3. Check that blocking_factor[i]>=8 if i<max_level.
    //    (this last check insures that there are at least 4 coarse 
    //     (level i) proper nesting cells next to a (level i+1) finer level)
    for (int i = 0; i <= max_level; i++) {
        int k = blocking_factor[i];
	if (i<max_level) {
         if (k<8) {
  	  std::cout << "must have at least 4 proper nesting cells" << '\n';
	  amrex::Error("must have blocking_factor>=8 if lev<max_level");
	 }
	 if (k<blocking_factor[i+1])
	  amrex::Error("blocking_factor[i]<blocking_factor[i+1]");
	}
        if (k<4)
         amrex::Error("blocking factor must be 4 or larger");

        while ( k > 0 && (k%2 == 0) )
            k /= 2;
        if (k != 1)
            amrex::Error("Amr::checkInputs: blocking_factor not power of 2");
    } // i=0 .. max_level

    for (int i = 0; i <= max_level; i++) {
        int k = space_blocking_factor[i];

        if (k>blocking_factor[i])
         amrex::Error("space_blocking_factor too big");
	if (i<max_level) {
	 if (k<space_blocking_factor[i+1])
	  amrex::Error("space_blocking_factor[i]<space_blocking_factor[i+1]");
	}

         // the number of coarse grid proper nesting cells for level i+1
         // is blocking_factor[i]/2
        if ((i>=0)&&(i<max_level)) {
         if (blocking_factor[i]<k)
          amrex::Error("bfact_grid>=space_blocking_Factor required");
         if (blocking_factor[i]<2*k)
          amrex::Error("bfact_grid>=2*space_blocking_Factor required");
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
        if (k>blocking_factor[0])
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
        if (len%blocking_factor[0] != 0)
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
        if (max_grid_size[i]%2 != 0)
            amrex::Error("max_grid_size is not even");
    }

    //
    // Check that max_grid_size is a multiple of blocking_factor at every level.
    //
    for (i = 0; i < max_level; i++) {
     if (max_grid_size[i]%blocking_factor[i] != 0)
      amrex::Error("max_grid_size not divisible by blocking_factor");
    }
    for (i = 0; i < max_level; i++) {
     if (max_grid_size[i]%space_blocking_factor[i] != 0)
      amrex::Error("max_grid_size not divisible by space_blocking_factor");
    }
    for (i = 0; i < max_level; i++) {
     if (max_grid_size[i]%time_blocking_factor != 0)
      amrex::Error("max_grid_size not divisible by time_blocking_factor");
    }

    if (!geom[0].ProbDomain().ok())
        amrex::Error("checkInput: bad physical problem size");

    if (max_level > 0) 
       if (regrid_int[0] <= 0)
          amrex::Error("checkinput: regrid_int not defined and max_level > 0");

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
       std::cout << "Successfully read inputs file ... " << '\n';
} // end subroutine checkInput ()

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
         writePlotFile(plot_file_root,level_steps[0],
          do_plot,do_slice,
          SDC_outer_sweeps,slab_step);
        }
    }
#ifdef HAS_XGRAPH
    if (first_plotfile)
    {
        first_plotfile = false;
        amr_level[0]->setPlotVariables();
    }
#endif
}

void
Amr::initialInit (Real strt_time,
                  Real stop_time)
{
 checkInput();

 finest_level = 0;
 FORT_OVERRIDE_FINEST_LEVEL(&finest_level);

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

 for (int lev = 0; lev <= finest_level; lev++)
  amr_level[lev]->post_regrid(0,finest_level,strt_time);

 // recomputes the timestep on all levels and updates statedata and dt_AMR
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
    FORT_OVERRIDE_FINEST_LEVEL(&finest_level);

    Vector<Box> inputs_domain(max_level+1);
    for (int lev = 0; lev <= max_level; lev++) {
     Box bx(geom[lev].Domain().smallEnd(),geom[lev].Domain().bigEnd());
     inputs_domain[lev] = bx;
    }

// SUSSMAN KLUGE RESTART

    int num_materials=0;
    is >> num_materials;
    if (num_materials<=0)
     amrex::Error("num_materials invalid in restart");
    is >> AMR_volume_history_recorded;
    AMR_volume_history.resize(num_materials);
    for (int j=0;j<num_materials;j++) {
     is >> AMR_volume_history[j];
     if (AMR_volume_history[j]<0.0)
      amrex::Error("cannot have negative volume_history");
    }

    int checkpoint_recalesce_data=0;
    for (int im=0;im<AMR_num_materials;im++) {
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
           level_count[0] = regrid_int[0];

       checkInput();
       //
       // Read levels.
       //
       for (int lev(0); lev <= finest_level; ++lev)
       {
           amr_level[lev].reset((*levelbld)());
           amr_level[lev]->restart(*this, is);
           this->SetBoxArray(lev, amr_level[lev]->boxArray());
           this->SetDistributionMap(lev, amr_level[lev]->DistributionMap());
       }
       //
       // Build any additional data structures.
       //

       for (int lev(0); lev <= finest_level; lev++)
           amr_level[lev]->post_restart();

    } else if ((max_level>=0)&&(max_level<mx_lev)) {

       if (ParallelDescriptor::IOProcessor())
          amrex::Warning("Amr::restart(): max_level is lower than before");

       int new_finest_level = std::min(max_level,finest_level);

       finest_level = new_finest_level;
       FORT_OVERRIDE_FINEST_LEVEL(&finest_level);
 
       // These are just used to hold the extra stuff we have to read in.
       Geometry   geom_dummy;
       int         int_dummy;
       IntVect intvect_dummy;

       for (i = 0          ; i <= max_level; i++) is >> geom[i];
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
           level_count[0] = regrid_int[0];

       checkInput();

       //
       // Read levels.
       //
       int lev;
       for (lev = 0; lev <= new_finest_level; lev++)
       {
           amr_level[lev].reset((*levelbld)());
           amr_level[lev]->restart(*this, is);
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

    if (ParallelDescriptor::IOProcessor())
    {
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

        int num_materials=AMR_volume_history.size();
        if (num_materials<=0)
         amrex::Error("num_materials invalid in checkpoint");
        HeaderFile << num_materials << '\n';
        HeaderFile << AMR_volume_history_recorded << '\n';
        for (int j=0;j<num_materials;j++) {
         HeaderFile << AMR_volume_history[j] << '\n';
         if (AMR_volume_history[j]<0.0) 
          amrex::Error("AMR_volume_history cannot be negative");
        }

        int checkpoint_recalesce_data=0;
        for (int im=0;im<AMR_num_materials;im++) {
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

    if (ParallelDescriptor::IOProcessor())
    {
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
Amr::timeStep (Real time,
               Real stop_time)
{

 if (fabs(time-cumtime)>1.0e-13)
  amrex::Error("time<>cumtime");

 if ((finest_level==0)&&(regrid_on_restart==1)) {

  regrid_on_restart = 0;

  const Box& domain = geom[0].Domain();
  IntVect d_length  = domain.size();
  const int* d_len  = d_length.getVect();

  for (int idir = 0; idir < AMREX_SPACEDIM; idir++)
   if (d_len[idir]%2 != 0)
    amrex::Error("timeStep: must have even number of cells");

  BoxArray lev0(1);
  lev0.set(0,amrex::coarsen(domain,2));
  lev0.maxSize(max_grid_size[0]/2);
  lev0.refine(2);

   // SUSSMAN
  int nprocs=ParallelDescriptor::NProcs();
  DistributionMapping dm(lev0,nprocs);

  AmrLevel* a = (*levelbld)(*this,0,geom[0],lev0,dm,cumtime);

   // calls setTimeLevel for level=0 using old level dt.
   // dm is the DistributionMapping on the new level 0.
   // lev0 is the BoxArray on the new level 0
  a->init(*amr_level[0],lev0,dm);
  amr_level[0].reset(a);
  this->SetBoxArray(0, amr_level[0]->boxArray());
  this->SetDistributionMap(0, amr_level[0]->DistributionMap());

   // calls CopyNewToOld 
   // calls setTimeLevel(cumtime,dt_AMR) 
  amr_level[0]->post_regrid(0,0,cumtime);

  if (ParallelDescriptor::IOProcessor()) {
   if (verbose > 1) {
    printGridInfo(std::cout,0,finest_level);
   } else if (verbose > 0) {
    printGridSummary(std::cout,0,finest_level);
   }
  }

  if (record_grid_info && ParallelDescriptor::IOProcessor())
   printGridInfo(gridlog,0,finest_level);
 
 } else if ((finest_level>0)||(regrid_on_restart==0)) {
  // do nothing
 } else
  amrex::Error("finest_level or regrid_on_restart invalid");


 int max_coarsest = std::min(finest_level, max_level-1);

 for (int level = 0; level <= max_coarsest; level++) {

  const int old_finest = finest_level;

  if (level_count[level] >= regrid_int[level]) {

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
  writePlotFile(plot_file_root,level_steps[0],
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

void Amr::recalesce_copy_new_to_old(int nmat) {

 if (nmat!=AMR_num_materials)
  amrex::Error("nmat invalid");

 int recalesce_num_state=6;

 if ((recalesce_state_old.size()!=nmat*recalesce_num_state)||
     (recalesce_state_new.size()!=nmat*recalesce_num_state)) {
   amrex::Error("recalesce sizes incorrect");
 } else {

  for (int im=0;im<nmat*recalesce_num_state;im++) {

   recalesce_state_old[im]=recalesce_state_new[im];

  } // im
 } // everything correct size?
} // recalesce_copy_new_to_old



void Amr::recalesce_copy_old_to_new(int nmat) {

 if (nmat!=AMR_num_materials)
  amrex::Error("nmat invalid");

 int recalesce_num_state=6;

 if ((recalesce_state_old.size()!=nmat*recalesce_num_state)||
     (recalesce_state_new.size()!=nmat*recalesce_num_state)) {
   amrex::Error("recalesce sizes incorrect");
 } else {
  for (int im=0;im<nmat*recalesce_num_state;im++) {

   recalesce_state_new[im]=recalesce_state_old[im];

  } // im
 } // everything correct size?
} // recalesce_copy_old_to_new



void Amr::recalesce_init(int nmat) {

 if (nmat!=AMR_num_materials)
  amrex::Error("nmat invalid");

 int recalesce_num_state=6;

 recalesce_state_old.resize(recalesce_num_state*nmat);
 recalesce_state_new.resize(recalesce_num_state*nmat);

 for (int im=0;im<recalesce_num_state*nmat;im++) {
  recalesce_state_old[im]=-1.0;
  recalesce_state_new[im]=-1.0;
 }

} // recalesce_init


void Amr::recalesce_get_state(Vector<Real>& recalesce_state_out,int nmat) { 

 if (nmat!=AMR_num_materials)
  amrex::Error("nmat invalid");

 int recalesce_num_state=6;

 if (recalesce_state_out.size()!=recalesce_num_state*nmat)
  amrex::Error("recalesce_state_out has incorrect size");
 if (recalesce_state_old.size()!=recalesce_num_state*nmat)
  amrex::Error("recalesce_state_old has incorrect size");

 for (int im=0;im<recalesce_num_state*nmat;im++)
  recalesce_state_out[im]=recalesce_state_old[im];

} // recalesce_get_state


void Amr::recalesce_put_state(Vector<Real>& recalesce_state_in,int nmat) {

 if (nmat!=AMR_num_materials)
  amrex::Error("nmat invalid");

 int recalesce_num_state=6;

 if (recalesce_state_new.size()!=recalesce_num_state*nmat)
  amrex::Error("recalesce_state_new has incorrect size");
 if (recalesce_state_in.size()!=recalesce_num_state*nmat)
  amrex::Error("recalesce_state_in has incorrect size");

 for (int im=0;im<recalesce_num_state*nmat;im++)
  recalesce_state_new[im]=recalesce_state_in[im];

} // recalesce_put_state


void
Amr::coarseTimeStep (Real stop_time)
{
    const double run_strt = ParallelDescriptor::second() ;

    //SUSSMAN
    //in: AMReX_FabArrayBase.cpp
    FabArrayBase::flushTileArrayCache();
    FabArrayBase::flushFBCache();
    FabArrayBase::flushCPCache();

     // check dt on all the levels.
    if (level_steps[0] > 0) {
        int post_regrid_flag = 0;
         // in AmrLevel.H: virtual void computeNewDt
         // NavierStokes::computeNewDt
        amr_level[0]->computeNewDt(finest_level,
         dt_AMR,stop_time,post_regrid_flag);
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
      long lorig=(long) cumtime/plot_per;
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
     writePlotFile(plot_file_root,level_steps[0],
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
    const Box& domain = geom[0].Domain();
    IntVect d_length  = domain.size();
    const int* d_len  = d_length.getVect();

    for (int idir = 0; idir < AMREX_SPACEDIM; idir++)
     if (d_len[idir]%2 != 0)
      amrex::Error("defBaseLevel: must have even number of cells");

    BoxArray lev0(1);
    lev0.set(0,amrex::coarsen(domain,2));
    lev0.maxSize(max_grid_size[0]/2);
    lev0.refine(2);

    this->SetBoxArray(0, lev0);

     // SUSSMAN
    int nprocs=ParallelDescriptor::NProcs();
    DistributionMapping dm(lev0,nprocs);
    this->SetDistributionMap(0, dm);

    //
    // Now build level 0 grids.
    //
    amr_level[0].reset((*levelbld)(*this,0,geom[0],grids[0],dmap[0],strt_time));

    amr_level[0]->initData();
} // subroutine defBaseLevel

// called from timeStep and bldFineLevels.
// bldFineLevels is called from initialInit.
//  (note: defBaseLevel is also called from initialInit)
// levelbld is called from defBaseLevel, timeStep, regrid, bldFineLevels
void
Amr::regrid (int  lbase,
             Real time,
             bool initial)
{

 if (fabs(time-cumtime)>1.0e-13)
  amrex::Error("time<>cumtime in regrid");

 if (verbose > 0 && ParallelDescriptor::IOProcessor())
  std::cout << "REGRID: at level lbase = " << lbase << std::endl;

 if (record_run_info && ParallelDescriptor::IOProcessor())
  runlog << "REGRID: at level lbase = " << lbase << '\n';

 int max_coarsest=std::min(finest_level,max_level-1);
 if (lbase>max_coarsest)
  amrex::Error("cannot have lbase>max_coarsest");

 int new_finest;
 Vector<BoxArray> new_grids(max_level+1);
 Vector<DistributionMapping> new_dmap(max_level+1);

 grid_places(lbase,new_finest,new_grids);
 if (new_finest>finest_level+1)
  amrex::Error("cannot create more than one new level at a time");

 int regrid_level_zero=0;
 if (lbase==0) {
  if (new_grids[0] != amr_level[0]->boxArray())
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
  if (new_grids[lev] == amr_level[lev]->boxArray()) {
   new_grids[lev] = amr_level[lev]->boxArray();  // to avoid duplicates
   new_dmap[lev] = amr_level[lev]->DistributionMap();
  } else {
   // do nothing
  }
 }  // lev=start ... min(finest_level,new_finest)

 finest_level = new_finest;
 FORT_OVERRIDE_FINEST_LEVEL(&finest_level);

 for (int lev = start; lev <= new_finest; lev++) {

  if (new_grids[lev].size()<1) {
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
   new_dmap[lev].define(new_grids[lev]);
  }

  AmrLevel* a = (*levelbld)(*this,lev,geom[lev],
    new_grids[lev],new_dmap[lev],cumtime);

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
   a->init(*amr_level[lev],new_grids[lev],new_dmap[lev]);
   amr_level[lev].reset(a);
   this->SetBoxArray(lev, amr_level[lev]->boxArray());
   this->SetDistributionMap(lev, amr_level[lev]->DistributionMap());

  } else {

   a->init(new_grids[lev],new_dmap[lev]);
   amr_level[lev].reset(a);
   this->SetBoxArray(lev, amr_level[lev]->boxArray());
   this->SetDistributionMap(lev, amr_level[lev]->DistributionMap());

  }

 } // lev=start..new_finest

 for (int lev = start; lev <= new_finest; lev++)
  amr_level[lev]->post_regrid(lbase,new_finest,time);

 if (record_run_info && ParallelDescriptor::IOProcessor()) {
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
} // subroutine regrid

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
        long                      ncells  = amr_level[lev]->countCells();
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
        long                      ncells  = amr_level[lev]->countCells();
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

void
Amr::ProjPeriodic (BoxList&        blout,
                   const Geometry& geom)
{
    //
    // Add periodic translates to blout.
    //
    Box domain = geom.Domain();

    BoxList blorig(blout);

    int nist,njst,nkst;
    int niend,njend,nkend;
    nist = njst = nkst = 0;
    niend = njend = nkend = 0;
    D_TERM( nist , =njst , =nkst ) = -1;
    D_TERM( niend , =njend , =nkend ) = +1;

    int ri,rj,rk;
    for (ri = nist; ri <= niend; ri++)
    {
        if (ri != 0 && !geom.isPeriodic(0))
            continue;
        if (ri != 0 && geom.isPeriodic(0))
            blorig.shift(0,ri*domain.length(0));
        for (rj = njst; rj <= njend; rj++)
        {
            if (rj != 0 && !geom.isPeriodic(1))
                continue;
            if (rj != 0 && geom.isPeriodic(1))
                blorig.shift(1,rj*domain.length(1));
            for (rk = nkst; rk <= nkend; rk++)
            {
                if (rk != 0 && !geom.isPeriodic(2))
                    continue;
                if (rk != 0 && geom.isPeriodic(2))
                    blorig.shift(2,rk*domain.length(2));

                BoxList tmp(blorig);
                tmp.intersect(domain);
                blout.catenate(tmp);
 
                if (rk != 0 && geom.isPeriodic(2))
                    blorig.shift(2,-rk*domain.length(2));
            }
            if (rj != 0 && geom.isPeriodic(1))
                blorig.shift(1,-rj*domain.length(1));
        }
        if (ri != 0 && geom.isPeriodic(0))
            blorig.shift(0,-ri*domain.length(0));
    }
}

// new_finest cannot be greater than finest_level+1
void
Amr::grid_places (int              lbase,
                  int&             new_finest,
                  Vector<BoxArray>& new_grids)
{
 int ilev=0;
 int max_crse = std::min(finest_level,max_level-1);

 const double strttime = ParallelDescriptor::second();

 if (lbase == 0) {
  const Box& domain = geom[0].Domain();
  IntVect d_length  = domain.size();
  const int* d_len  = d_length.getVect();

  for (int idir = 0; idir < AMREX_SPACEDIM; idir++)
   if (d_len[idir]%2 != 0)
    amrex::Error("grid_places: must have even number of cells");

  BoxArray lev0(1);
  lev0.set(0,amrex::coarsen(domain,2));
  lev0.maxSize(max_grid_size[0]/2);
  lev0.refine(2);

  new_grids[0] = lev0;

 } // lbase==0 

  // blocking_factor[i] is the blocking factor for level i.
  //  (why did I say level i+1 previously?)
  // The proper nesting buffer of level i cells next to level i+1
  // cells is bf_lev[i]*n_proper.
 Vector<int> bf_lev(max_level); 
 Vector<int> rr_lev(max_level);
 Vector<Box> pc_domain(max_level);  // Coarsened problem domain.

 // blocking_factor is a power of 2 and no smaller than 4.
 for (ilev = 0; ilev <= max_crse; ilev++) {
   bf_lev[ilev] = blocking_factor[ilev]/2;
   if (2*bf_lev[ilev]!=blocking_factor[ilev])
    amrex::Error("2*bf_lev[ilev]!=blocking_factor[ilev]");
 }

 for (ilev = lbase; ilev < max_crse; ilev++) {
   rr_lev[ilev] = (2*bf_lev[ilev])/bf_lev[ilev+1];
   if (rr_lev[ilev]<2)
    amrex::Error("rr_lev[ilev]<2");
 }

  //2*bf_lev[ilev]==blocking_factor[ilev]
 for (ilev = lbase; ilev <= max_crse; ilev++) {
  pc_domain[ilev] = amrex::coarsen(geom[ilev].Domain(),bf_lev[ilev]);
 }

 Vector<BoxList> p_n(max_level);      // Proper nesting domain.
 Vector<BoxList> p_n_comp(max_level); // Complement proper nesting domain.

 BoxList bl(amr_level[lbase]->boxArray());
 bl.simplify();
  //2*bf_lev[i]==blocking_factor[i]
 bl.coarsen(bf_lev[lbase]);
 p_n_comp[lbase].complementIn(pc_domain[lbase],bl);

 p_n_comp[lbase].simplify();
  // grow each box in p_n_comp[lbase] by n_proper
  // proper nesting size: n_proper*bf_lev
 p_n_comp[lbase].accrete(n_proper);
 Amr::ProjPeriodic(p_n_comp[lbase], Geometry(pc_domain[lbase]));
 p_n[lbase].complementIn(pc_domain[lbase],p_n_comp[lbase]);
 p_n[lbase].simplify();

 if (lbase==0) {
  if (p_n_comp[lbase].size()==0) {
   // do nothing
  } else
   amrex::Error("p_n_comp[lbase].size() should be 0");
 } // lbase==0

 bl.clear();

 for (ilev = lbase+1; ilev <= max_crse; ilev++) {

  p_n_comp[ilev] = p_n_comp[ilev-1];

  p_n_comp[ilev].simplify();

  p_n_comp[ilev].refine(rr_lev[ilev-1]);
   // grow each box in p_n_comp[ilev] by n_proper
  p_n_comp[ilev].accrete(n_proper);

  Amr::ProjPeriodic(p_n_comp[ilev], Geometry(pc_domain[ilev]));

  p_n[ilev].complementIn(pc_domain[ilev],p_n_comp[ilev]);
  p_n[ilev].simplify();

  if (lbase==0) {
   if (p_n_comp[ilev].size()==0) {
    // do nothing
   } else
    amrex::Error("p_n_comp[ilev].size() should be 0");
  } // lbase==0

 } // ilev=lbase+1 ... max_crse

 new_finest = lbase;

 for (int levc = max_crse; levc >= lbase; levc--) {

  int levf = levc+1;
  int ngrow = 0;

  if (levf < new_finest) {
   BoxArray ba_proj(new_grids[levf+1]);

   ba_proj.coarsen(2);
   ba_proj.grow(n_proper);
   ba_proj.coarsen(2);

   BoxArray levcBA = amr_level[levc]->boxArray();

   while (!levcBA.contains(ba_proj)) {
    BoxArray tmp = levcBA;
    tmp.grow(1);
    levcBA = tmp;
    ngrow++;
   }
  }  // levf<new_finest

   // TagBox.H: TagBoxArray (const BoxArray& bs,dm,int _ngrow=0)
    // SUSSMAN
  TagBoxArray tags(amr_level[levc]->boxArray(),
		   amr_level[levc]->DistributionMap(),
		   n_error_buf[levc]+ngrow);

  amr_level[levc]->errorEst(tags,
                           TagBox::CLEAR,TagBox::SET,
                           n_error_buf[levc],ngrow);

  if (levf < new_finest) {

   int nerr = n_error_buf[levf];

   BoxList bl_tagged(new_grids[levf+1]);
   bl_tagged.simplify();
   bl_tagged.coarsen(2);
   for (BoxList::iterator blt = bl_tagged.begin(), End = bl_tagged.end();
        blt != End;
        ++blt) {
    for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
     if (blt->smallEnd(idir) == geom[levf].Domain().smallEnd(idir))
      blt->growLo(idir,nerr);
     if (blt->bigEnd(idir) == geom[levf].Domain().bigEnd(idir))
      blt->growHi(idir,nerr);
    }
   } // blt

   Box mboxF = amrex::grow(bl_tagged.minimalBox(),1);
   BoxList blFcomp;
   blFcomp.complementIn(mboxF,bl_tagged);
   blFcomp.simplify();
   bl_tagged.clear();

   int iv=nerr/2;
    // Grow each Box in the BoxList by iv.
   blFcomp.accrete(iv);
   BoxList blF;
   blF.complementIn(mboxF,blFcomp);
   BoxArray baF(blF);
   blF.clear();
   baF.grow(n_proper);
   for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    if (nerr > n_error_buf[levc]*2) 
     baF.grow(idir,nerr-n_error_buf[levc]*2);
   }

   baF.coarsen(2);

   tags.setVal(baF,TagBox::SET);
  }  // levf < new_finest

  int nbuf_all=n_error_buf[levc]+ngrow;
  IntVect nbuf_vec(D_DECL(nbuf_all,nbuf_all,nbuf_all));

  tags.buffer(nbuf_vec);

  int bl_max = bf_lev[levc];
  if (bl_max>=2) {
   int bf_all=bf_lev[levc];
   IntVect bf_vec(D_DECL(bf_all,bf_all,bf_all));
   tags.coarsen(bf_vec); // guarantee proper nesting of n_proper*bf_lev
  } else
   amrex::Error("blocking_factor>=4 required => bf_lev>=2");

  amr_level[levc]->manual_tags_placement(tags, bf_lev);
  tags.mapPeriodic(Geometry(pc_domain[levc]));
  tags.setVal(p_n_comp[levc],TagBox::CLEAR);

  Vector<IntVect> tagvec;
  tags.collate(tagvec);
  tags.clear();

  if (tagvec.size() > 0) {
   new_finest = std::max(new_finest,levf);
   ClusterList clist(&tagvec[0],tagvec.size());
   clist.chop(grid_eff);
   BoxDomain bd;
   bd.add(p_n[levc]);
   clist.intersect(bd);
   bd.clear();

   BoxList new_bx;
   clist.boxList(new_bx);
   new_bx.refine(bf_lev[levc]);
   new_bx.simplify();
   BL_ASSERT(new_bx.isDisjoint());

   int largest_grid_size;
   largest_grid_size = max_grid_size[levf] / 2;
   new_bx.maxSize(largest_grid_size);

   new_bx.refine(2);
   BL_ASSERT(new_bx.isDisjoint());
   new_grids[levf].define(new_bx);
  }  // tagvec.size()>0
 } // levc=max_crse ... lbase; levc--

 const int NProcs = ParallelDescriptor::NProcs();

 if ((NProcs > 1)&&(refine_grid_layout)) {

  for (int cnt = 1; cnt <= 4; cnt *= 2) {

   for (ilev = lbase; ilev <= new_finest; ilev++) {

    const int ChunkSize = max_grid_size[ilev]/cnt;

    IntVect chunk(D_DECL(ChunkSize,ChunkSize,ChunkSize));

    for (int j = 0; j < AMREX_SPACEDIM; j++) {
     chunk[j] /= 2;

     if ((new_grids[ilev].size() < NProcs) && 
         (chunk[j]%blocking_factor[ilev] == 0)) {
      new_grids[ilev].maxSize(chunk);
     }
    }
   } // ilev
  } // cnt
 }

 if (verbose > 0) {
  double stoptime = ParallelDescriptor::second() - strttime;

  ParallelDescriptor::ReduceRealMax(stoptime,
        ParallelDescriptor::IOProcessorNumber());

  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "grid_places() time: " << stoptime << '\n';
  }
 }  // verbose>0

 if (new_finest>finest_level+1)
  amrex::Error("cannot create more than one new level at a time");

}  // subroutine grid_places

// called from initialInit
void
Amr::bldFineLevels (Real strt_time)
{
 if (max_level<=0)
  amrex::Error("max_level invalid in bldFineLevels");

 finest_level = 0;
 FORT_OVERRIDE_FINEST_LEVEL(&finest_level);

 Vector<BoxArray> new_grids(max_level+1);

 int new_finest=0;
 int previous_new_finest=0;
 int num_grid_places=0;
 int grid_places_done=0;

 while (grid_places_done==0) {

  grid_places(finest_level,new_finest,new_grids);

  if (new_finest>finest_level) {

   if (new_finest>finest_level+1)
    amrex::Error("cannot create more than one new level at a time");

   finest_level = new_finest;
   FORT_OVERRIDE_FINEST_LEVEL(&finest_level);

   if (new_grids[new_finest].size()<1) {
    std::cout << "IN bldFineLevels\n";
    std::cout << "finest_level= " << finest_level << '\n';
    std::cout << "new_finest= " << new_finest << '\n';
    std::cout << "new_grids[0] \n";
    std::cout << new_grids[0] << '\n';
    amrex::Error("new_grids[new_finest] invalid");
   }

    // SUSSMAN
   int nprocs=ParallelDescriptor::NProcs();
   DistributionMapping new_dm(new_grids[new_finest],nprocs);

    // see the constructor in AmrLevel.cpp:
    // AmrLevel::AmrLevel ( ....  )
   AmrLevel* a_level = (*levelbld)(*this,
                         new_finest,
                         geom[new_finest],
                         new_grids[new_finest],
                         new_dm,
                         strt_time);

   amr_level[new_finest].reset(a_level);
   this->SetBoxArray(new_finest, new_grids[new_finest]);
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
    new_grids[i] = amr_level[i]->boxArray();

   regrid(0,strt_time,true);

   grids_the_same = true;

   for (int i = 0; i <= finest_level && grids_the_same; i++)
    if (!(new_grids[i] == amr_level[i]->boxArray()))
     grids_the_same = false;

   count++;
 } while (!grids_the_same && count < MaxCnt);


}  // subroutine bldFineLevels

} // namespace amrex

