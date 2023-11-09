module multidomain_design_mod
!!
!! Set the main variables controlling the multidomain design
!!
!! Most variables can be overridden at runtime via the namelist inputs, see 
!!   subroutine read_multidomain_design_variables
!! below.

use global_mod, only : dp, ip, charlen
use stop_mod, only: generic_stop
use logging_mod, only: log_output_unit

implicit none

public

!
! Variables in &MULTIDOMAIN_GLOBAL_PROPERTIES begin here
!

! Lower-left corner coordinate of multidomain in degrees lon,lat.
! Should match lower-left corner in ../multidomain_design/create_boxes.R
real(dp) ::  global_ll(2) = [16.0_dp, -75.0_dp]

! Length/width of multidomain in degrees lon,lat
real(dp) :: global_lw(2) = [156.0_dp , 33.0_dp] - [16.0_dp, -75.0_dp]

! Is the multidomain periodic in the east-west direction? The code will throw
! an error if set to .true. AND global_lw(1) is not very close to 360 degrees.
logical ::  use_periodic_EW_multidomain = .false.

! Cell-size of the coarsest domains. The cell size of nested domains is
! defined relative to this. If should exactly divide global_lw, but if 
! not is adjusted to a nearby value.
real(dp) ::  global_dx_arcmin = 1.0_dp !2.0_dp

! Global time-step in the multidomain. Must satisfy CFL condition
! everywhere, in combination with local timestepping (non-staggered-grid domains)
! as controlled by the nested domains' input files. Start with a guess, 
! then check the stationary_timestep information (in the logs) to refine it.
real(dp) ::  global_dt = 2.2_dp ! * global_dx_arcmin

! Final time (for non-test runs) in seconds. Start time is always 0.0
real(dp) ::  final_time_full_runs = 24.0_dp * 3600.0_dp
! Final time for test runs (e.g. used for load balancing or simple checks)
real(dp) ::  final_time_test_runs = 360.0_dp

! Elevation file names listed in priority order
character(len=charlen) :: swals_elevation_files_in_preference_order = &
    "../elevation/swals_elevation_files_in_preference_order.txt"

! Override the initial stage in user provided polygons. If non-empty, this
! should link to a csv file with 2 columns (and no header), being the polygon
! csv file, and the value to set. e.g.:
!     full_path_to_polygon_file_1.csv, -0.2 
!     full_path_to_polygon_file_2.csv, 1.2 
!     ....
! where the format of files like "full_path_to_polygon_file_1.csv" is the same
! as files specifying the breakwalls/inverts etc.
character(len=charlen) :: override_initial_stage_polygons_values_file = ""

! Workaround for raster missing data. 
! Many rasters denote NA with a large negative value, but it often becomes 
! slightly mangled (e.g. by precision transformations) and then may not be correctly 
! interpreted as NA. To workaround such cases, we can treat all raster values below a 
! large negative threshold as NA (so long as they are never "genuine" values).
real(dp) :: raster_na_below_limit = -1.0e+20_dp

! Optionally do local smoothing of elevation along coarse-to-fine nesting boundaries.
! This can improve stability
logical :: smooth_elevation_along_nesting_boundaries = .TRUE.

! Breakwall file names listed in priority order. These define 3D lines that will
! be burned into the elevation grid (unless the elevation is already higher)
character(len=charlen) :: breakwalls_file_list = "../breakwalls/swals_breakwall_files.txt"

! Invert file names listed in priority order. These define 3D lines that will be
! burned into the elevation grid (unless the elevation is already lower)
character(len=charlen) :: inverts_file_list = ""
 
! Point gauges (where we store time-series outputs)
character(len=charlen) :: swals_point_gauge_file = "../gauges/point_gauges_2022_12_14.csv"
 
! Approx timestep between any file outputs (seconds). The frequency of 
! particular kinds of outputs can be reduced relative to this (below).
real(dp) :: approximate_writeout_timestep = 30.0_dp

! Optionally write grids less often than the approximate_writeout_timestep, to 
! keep file-size down. 
integer(ip) :: write_grids_every_nth_writeout_timestep = 400_ip
 
! Optionally print outputs less often than the writeout timestep, to keep the
! log file-size down
integer(ip) :: print_every_nth_writeout_timestep = 10_ip
 
! Optionally write gauges less often than the writeout timestep, to keep the
! file-size down
integer(ip) :: write_gauges_every_nth_writeout_timestep = 1_ip

! Which gridded variables are stored every timestep? Possible values are
! [character(len=charlen) :: "stage", "uh", "vh", "elev"]
character(len=charlen) :: time_grids_to_store(4) = ""

! Which gridded variables are stored once? Possible values are
! [character(len=charlen):: 'max_stage', 'max_flux', 'max_speed', 'arrival_time', &
!  'manning_squared', 'elevation0', 'elevation_source_file_index', 'time_of_max_stage']
character(len=charlen) :: nontemporal_grids_to_store(8) = ""

! Number of nesting levels (including one for the global domain). For instance
! if you had a model with a coarse outer grid (dx = 1.0) and a set of nested 
! inner grids (all defined on the same nesting level) then then NNL=2. 
! More details are specified in MULTIDOMAIN_NESTING_PROPERTIES.
integer(ip) :: NNL = 4_ip

!
! Variables in &MULTIDOMAIN_NESTING_PROPERTIES begin here
!
!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Specify nesting information for the "NNL" nesting levels, and how the work is
! distributed among MPI processes.
!
! We make the following assumptions (some are more restrictive than what SWALS
! allows in general)
! 1. Nesting level 'i' contains one or more domains specified by the 
!    nesting_domain_extents_file(i) (except for nesting level 1).
! 2. Domains on nesting level 'i' are entirely contained by domains on nesting 
!    level 'i-1'
! 3. The default grid-size of nesting level 'i' domains is an integer divisor
!    of the default grid-size on nesting level 'i-1'. The integer is specified 
!    by dx_refinement_factor(i).
!
! It is possible to coarsen specific domains on nesting level 'i' by setting
! the "coarsen" column in nesting_domain_extents_file(i) to 1 (rather than 
! 0). The coarsened domains will have:
!       grid size = (coarsening_factors(i) * regular grid size on level 'i')
! The user must ensure the grid size still meets the constraint that when 
! domains are nesting, the finer domain has a cell size that is an integer 
! divisor of the coarser domain.
! 

! Input files defining the domains on each nesting level. For the global domain
! it should always be empty. (The meaning of columns in each file is partly
! specified further below).
character(len=charlen), allocatable :: nesting_domain_extents_file(:) 

! Timestepping method for domains on each nesting level. 
character(len=charlen), allocatable :: nesting_domain_timestepping_method(:) 
! In practice there may be stability problems if using staggered grid solvers 
! at any level other than the global domain. Also the staggered grid solvers 
! cannot do local timestepping (due to stability issues). But in global scale 
! tsunami problems, staggered solvers are usually most efficient for the global 
! domain.

! Grid size refinement of domains in each nesting level, relative to the next 
! coarser nesting level. The global domain should always have a value 1
integer(ip), allocatable :: dx_refinement_factors(:)

! If the input file tells us to coarsen some domains in the nesting layer, then 
! this value defines how much they coarsen. The coarsening_factors must be an integer
! divisor of the corresponding entry in dx_refinement_factors. Coarsening can be 
! useful to speed up some domains in deep water, which otherwise slow models down, if
! they don't require high resolution for accuracy. 
! Negative values imply no coarsening will be used. 
integer(ip), allocatable :: coarsening_factors(:)

! If 'rk2' timestepping is used, this controls the 'theta' limiter parameter.
! In many problems a non-TVD value (>2) is better to reduce dissipation. From
! experience a value of 4.0 behaves well. 
real(dp), allocatable :: theta_fv_limiter(:)
 
! Load balance file, used to distribute work in parallel with MPI. 
character(len=charlen) ::  load_balance_file = "" 
! The file contains one row per domain. The top row is the global domain, then 
! all domains on nesting level 2, nesting level 3, ... . The domain order matches 
! that in the nesting_domain_extents_file. 
! Each row has one or more integers that correspond the coarray image indices 
! (or equivalantly "mpi_rank - 1"). Those images do the calculations on that 
! domain. If there is more than one image on a domain, SWALS will partition 
! that domain into tiles with approximately equal size. Each image can appear 
! multiple times (both within a row and on different rows). 
!
! EXAMPLE load_balance_file with 5 domains, for 4 MPI processes
!     1, 2, 3, 4
!     1
!     2, 3
!     4,
!     3, 3 
! The above would split domain 1 into 6 pieces, and domains 3 and 5 into 2 
! pieces, while the other domains would have 1 piece. On domain 5, both pieces 
! would be run on image 3.
!
! If the file is empty then SWALS will split every domain with one piece per 
! image, which probably isn't what you want for complex models.
!
! If the file contains integers that exceed the number of MPI processes 
! (MPI_COMM_SIZE), then the integer MPI_COMM_SIZE+1 will be mapped to image 1,
! MPI_COMM_SIZE+2 will be mapped to image 2, etc.
!
! In practice a default load balance file is created, and a short model is run.
! The logfiles from that model run are the used to develop a new load balance 
! file with more equal work between processes, using 
! 

! To hold domain metadata we will make an array with one entry per nest level.
! Use a dedicated type, and make clear the mapping between columns and data
type real_table_type
    real(dp), allocatable :: metadata(:,:) ! For a single nesting level
end type

! Hold domain input metadata for each nest level (not included in namelist)
type(real_table_type), allocatable :: nesting_domain_extents(:)

! These help setup the nesting (not included in namelist)
integer(ip), allocatable :: num_domains_per_nesting_level(:), &
    domain_index_to_nest_with(:)

! Mapping between columns of files provided in nesting_domain_extents_file(:),
! and variables interest. This should only be changed if you alter the file format.
integer(ip), parameter :: &
    r_ll(2) = [1_ip, 2_ip], & ! index of Lower left lon,lat
    r_ur(2) = [3_ip, 4_ip], & ! index of Upper right lon,lat
    r_dpth = 5_ip, &          ! index of Max depth (estimate) in domain
    r_coarse = 7_ip           ! index of "Coarsen the domain?" (0.0 = No, 1.0 = Yes)

! Variables that can be optionally read from &MULTIDOMAIN_GLOBAL_PROPERTIES
namelist /MULTIDOMAIN_GLOBAL_PROPERTIES/ &
    global_ll, global_lw, use_periodic_EW_multidomain, &
    global_dx_arcmin, global_dt, &
    final_time_full_runs, final_time_test_runs, &
    swals_elevation_files_in_preference_order, &
    raster_na_below_limit, &
    smooth_elevation_along_nesting_boundaries, &
    breakwalls_file_list, inverts_file_list, swals_point_gauge_file, &
    approximate_writeout_timestep, &
    write_grids_every_nth_writeout_timestep, &
    print_every_nth_writeout_timestep, &
    write_gauges_every_nth_writeout_timestep, &
    time_grids_to_store, nontemporal_grids_to_store, &
    NNL

! Variables controlling details of the nesting levels. Their size depends on 
! NNL, so they have to be allocated after that is known
namelist /MULTIDOMAIN_NESTING_PROPERTIES/ &
    nesting_domain_extents_file, &
    nesting_domain_timestepping_method, &
    dx_refinement_factors, &
    coarsening_factors, &
    theta_fv_limiter, &
    load_balance_file

contains

subroutine read_multidomain_design_variables(multidomain_design_namelist_file)
    !! Enable multidomain variables to be set at run-time
    character(len=charlen), intent(in) :: multidomain_design_namelist_file
        !! File containing two namelists:
        !!     1. &MULTIDOMAIN_GLOBAL_PROPERTIES
        !!     2. &MULTIDOMAIN_NESTING_PROPERTIES
        !! Namelist 2 contains all information about nesting levels EXCEPT the 
        !! number of nesting levels (NNL), which must be specified in 
        !! namelist 1 (so we can allocate memory for arrays in namelist 2)
    integer :: fid, ierr

    ! Open the file
    open(newunit=fid, file=multidomain_design_namelist_file, action='read', &
        iostat=ierr)
    if(ierr /= 0) then
        write(log_output_unit, *) 'Invalid namelist input file ', &
            trim(multidomain_design_namelist_file)
        call generic_stop
    end if
   
    ! Read the global properties 
    read(nml=MULTIDOMAIN_GLOBAL_PROPERTIES, unit=fid, iostat=ierr)
    if(ierr /= 0) then
        write(log_output_unit, *) &
            'Problem reading &MULTIDOMAIN_GLOBAL_PROPERTIES in ', &
            trim(multidomain_design_namelist_file)
        call generic_stop
    end if

    write(unit=log_output_unit, nml=MULTIDOMAIN_GLOBAL_PROPERTIES)

    ! NNL is now known. 
    ! Allocate space for nesting variables
    allocate(nesting_domain_extents_file(NNL), &
        nesting_domain_timestepping_method(NNL), dx_refinement_factors(NNL), &
        coarsening_factors(NNL), theta_fv_limiter(NNL))

    ! Set defaults that will throw errors if later set incorrectly
    nesting_domain_extents_file = ""
    dx_refinement_factors = -1_ip
    coarsening_factors = -1_ip 
    ! Set good defaults for our application
    theta_fv_limiter = 4.0_dp
    nesting_domain_timestepping_method(1) = "leapfrog_linear_plus_nonlinear_friction"
    nesting_domain_timestepping_method(2:NNL) = "rk2"

    ! Read the nesting variables 
    read(nml=MULTIDOMAIN_NESTING_PROPERTIES, unit=fid, iostat=ierr)
    if(ierr /= 0) then
        write(log_output_unit, *) &
            'Problem reading &MULTIDOMAIN_NESTING_PROPERTIES in ', &
            trim(multidomain_design_namelist_file)
        call generic_stop
    end if

    write(unit=log_output_unit, nml=MULTIDOMAIN_NESTING_PROPERTIES)

    ! Setup some other useful variables (not set by the user)
    allocate(nesting_domain_extents(NNL), num_domains_per_nesting_level(NNL), &
        domain_index_to_nest_with(NNL))   
    ! Defaults will cause errors if not set properly later
    num_domains_per_nesting_level = -1_ip
    domain_index_to_nest_with = -1_ip

end subroutine

subroutine write_multidomain_design_variables_to_logfile
    !! Write all namelist variables to the logfile.
    
    write(unit=log_output_unit, nml=MULTIDOMAIN_GLOBAL_PROPERTIES)
    write(unit=log_output_unit, nml=MULTIDOMAIN_NESTING_PROPERTIES)

end subroutine

end module multidomain_design_mod
