max_step  = 99999    # maximum timestep
#max_step  =  2    # maximum timestep
stop_time =  1000  # maximum problem time

# ------------------  INPUTS TO CLASS AMR ---------------------
# Figure 1E
# expected effective diameter: 0.499 cm
# nozzle radius 0.085cm
# multiply computed effective diameter by 0.085 
# We expect the run.out file to have 0.499/0.085=5.87
#
# The dimensionless radius is 1.0
# The dimensionless inflow average velocity is 1.0
#
geometry.coord_sys      = 1        # 0 => cart, 1 => RZ
geometry.prob_lo   =  0.0 0.0
geometry.prob_hi   =  16.0 48.0
geometry.is_periodic = 0 0 0     

ns.EILE_flag=-1
# =0 MOF =2 CMOF =4 CLSVOF/CMOF
ns.continuous_mof=2
ns.force_cmof_at_triple_junctions=1
ns.partial_cmof_stencil_at_walls=1
ns.MOFITERMAX=15
ns.vof_height_function=1
ns.FD_curv_interp=1

# 0MGPCG 1PCG 2 MINV=I
ns.project_solver_type=0
Lp.smooth_type=2
Lp.bottom_smooth_type=2
Lp.gpucount=1

# multigrid class
#mg.verbose = 2
#cg.verbose = 2
# set above to 2 for maximum verbosity
mg.nu_f = 40
mg.nu_0 = 1   # 1 - v-cycle 2 - w-cycle
cg.maxiter = 200
mg.bot_atol = 1.0e-12
mg.rtol_b = -0.01
#Lp.v = 1

ns.num_materials=3
ns.num_species_var=0

ns.post_init_pressure_solve=1 

ns.minimum_relative_error=1.0e-18
ns.initial_project_cycles=3  # number of initial jacobi iterations
ns.initial_cg_cycles=1  # number of MGPCG steps to do in order to get
                         # decent sized residual.

amr.n_cell    = 16 48 # min dx=16/(16 * 8)
amr.max_level = 3
ns.projection_pressure_scale=100.0

# 0- 1 level 1- 2 levels  2- 3 levels
amr.regrid_int      = 1       # how often to regrid
amr.n_error_buf     = 4 4 4 4 4    # number of buffer cells in error est
amr.grid_eff        = 0.8   # what constitutes an efficient grid
ns.enable_spectral=0
amr.time_blocking_factor=1
#amr.blocking_factor = 16 16 8 4       # block factor in grid generation
#amr.space_blocking_factor = 8 8 4 2  # block factor in grid generation
amr.blocking_factor = 8 8 8 8       # block factor in grid generation
amr.space_blocking_factor = 1 1 1 1  # block factor in grid generation
amr.check_int       = 400      # number of timesteps between checkpoints
amr.check_file      = chk     # root name of checkpoint file
amr.plot_int        = 200  # 200 for production run
amr.plot_file       = plt 
amr.grid_log        = grdlog  # name of grid logging file
amr.max_grid_size   = 1024
#amr.restart         = chk13600
#amr.trace   =1

# ------------------  INPUTS TO PHYSICS CLASS -------------------
ns.dt_cutoff      = 0.000001  # level 0 timestep below which we halt

mac.mac_abs_tol        = 1.0e-10   # tolerence for mac projections

amr.plotfile_on_restart=1
#ns.visual_revolve=32

ns.cfl            = 0.5      # cfl number for hyperbolic system
ns.init_shrink    = 1.0      # scale back initial timestep
ns.change_max     = 1.1      # scale back initial timestep
ns.visc_coef      = 1.0      # was 0.2749, now included in viscconst
                             # must have visc_coef=1 if some shear-thinning
                             # materials
mac.visc_abs_tol   = 1.0e-10
ns.gravity        = -0.0429    # body force  (gravity in MKS units)
# 12, 13, 23
# material 1 is liquid (ambient phase)
# material 2 will be gas
# material 3 is solid
# sigma_{i,j}cos(theta_{i,k})=sigma_{j,k}-sigma_{i,k}
# 90 degree contact angle:
ns.tension        = 0.3269 0.3269 0.3269 
#ns.fixed_dt	  = 0.0025     # hardwire dt
ns.sum_interval   = 5        # timesteps between computing mass 

ns.axis_dir=5  # case (e) from Helsby and Tuson
ns.vorterr=0.0 0.0 0.0
ns.rgasinlet=1.57
ns.vinletgas=0.0
ns.twall=0.0
ns.advbot=0.0  # inflow velocity prescribed in PROB_2D.F90 (inflow_bc)
               # "inflow_bc" sets up Poiseuille flow
               # quantity in cell is average velocity for Poiseuille flow

ns.adv_vel=0.0
ns.adv_dir=1

# these are shear thinning parameters
# viscconst*(1+(beta*gamma_dot)**alpha)**((n-1)/alpha)
# material 1 is the ambient material, so make
# beta>0 for material 1, and leave beta=0 for the other 2 materials.
ns.Carreau_alpha=0.0 0.0 0.0
ns.Carreau_beta=0.0 0.0 0.0
ns.Carreau_n=0.0 0.0 0.0

# these are viscoelastic parameters for future reference.
# Bingham model can be added too in the future.
ns.elastic_time=0.0 0.0 0.0
ns.elastic_viscosity=0.0 0.0 0.0
ns.polymer_factor=0.0 0.0 0.0

ns.heatviscconst=0.0 0.0 0.0 
ns.viscconst=0.2749 3.970655E-5 0.2749
ns.denconst=1.0 0.000985 1.0 
ns.tempconst=293.0 293.0 293.0 
ns.material_type=0 0 999  # 999 is reserved for prescribed solid. (nozzle)
ns.FSI_flag=0 0 1
ns.pressure_error_cutoff=0.0 0.0 0.0 
ns.xblob=0.0
ns.yblob=0.0
ns.zblob=1.0  # height of nozzle
ns.radblob=1.0  # radius of nozzle
ns.radblob2=0.25 # thickness of nozzle
ns.denfact=1.0
ns.velfact=0.0
ns.probtype=25

proj.bogus_value = 5.0e+5
proj.Pcode = 0
#proj.Pcode = 2

#ns.mem_debug = 1
ns.v = 1
#ns.d = 1
ns.output_drop_distribution=1

# ----------------  PROBLEM DEPENDENT INPUTS
ns.lo_bc          = 3 1 
ns.hi_bc          = 5 2

# >>>>>>>>>>>>>  BC FLAGS <<<<<<<<<<<<<<<<
# 0 = Interior           3 = Symmetry
# 1 = Inflow             4 = SlipWall
# 2 = Outflow            5 = NoSlipWall

# turn any of these on to generate run-time timing stats


# select single or double precision of FAB output data
#        default is whatever precision code is compiled with.
#fab.precision = FLOAT     # output in FLOAT or DOUBLE
fab.precision = DOUBLE    # output in FLOAT or DOUBLE

