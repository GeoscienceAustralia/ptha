module local_routines 
    !!
    !! NTHMP benchmark problem 9 -- Okushiri tsunami field test case
    !!

    use global_mod, only: dp, ip, charlen
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use read_raster_mod, only: multi_raster_type
    use logging_mod, only: log_output_unit
    implicit none

    contains 

    subroutine set_initial_conditions_BP09(domain, all_dx_md)
        type(domain_type), intent(inout):: domain
            !! The domain
        real(dp), intent(in), optional :: all_dx_md(:,:,:)
            !! Metadata on multidomain grid sizes, used to smooth along multidomain boundaries.

        integer(ip):: i, j
        character(len=charlen):: input_elevation(6), input_stage(1)
        real(dp), allocatable:: x(:), y(:)
        type(multi_raster_type):: elevation_data, stage_data
        real(dp), allocatable :: random_uniform(:)

        ! Use this to check the effect of very small perturbations of the initial condition
        ! Eventually  a tiny perturbation will eventually lead to differences in the results. 
        ! This seems comparable to the differences that eventually emerge with different 
        ! domain partitioning (with can be eliminated by providing a load_balance_partition.txt 
        ! file)
        real(dp), parameter :: random_perturbation_scale = 0.0e-10_dp

        ! Stage
        domain%U(:,:,STG) = 0.0e-0_dp

        ! Input rasters
        input_stage(1) = "../test_repository/BP09-FrankG-Okushiri_island/initial_condition_raster/HNO1993.tif"
        call stage_data%initialise(input_stage)

        ! Preference order for elevation: input_elevation(1) > input_elevation(2) > .... > input_elevation(6)
        input_elevation(6) = "../test_repository/BP09-FrankG-Okushiri_island/bathymetry_rasters_continuous/OK24.tif"
        input_elevation(5) = "../test_repository/BP09-FrankG-Okushiri_island/bathymetry_rasters_continuous/OK08.tif"
        input_elevation(4) = "../test_repository/BP09-FrankG-Okushiri_island/bathymetry_rasters_continuous/OK03.tif"
        input_elevation(3) = "../test_repository/BP09-FrankG-Okushiri_island/bathymetry_rasters_continuous/MO01.tif"
        input_elevation(2) = "../test_repository/BP09-FrankG-Okushiri_island/bathymetry_rasters_continuous/AO15.tif"
        input_elevation(1) = "../test_repository/BP09-FrankG-Okushiri_island/bathymetry_rasters_continuous/MB05.tif"
        call elevation_data%initialise(input_elevation)

        ! Make space for x/y coordinates, at which we will look-up the rasters
        allocate(x(domain%nx(1)), y(domain%nx(1)))
        x = domain%x

        allocate(random_uniform(domain%nx(1)))
        
        ! Set stage and elevation row-by-row.
        do j = 1, domain%nx(2)
            y = domain%y(j)

            ! Set stage - and clip 'NA' regions (since the stage raster does no cover the entire domain)
            call stage_data%get_xy(x,y, domain%U(:,j,STG), domain%nx(1), bilinear=1_ip, na_below_limit=-1.0e+20_dp)
            where(domain%U(:,j,STG) <= -1.0e+20_dp ) domain%U(:,j,STG) = 0.0_dp 

            ! Set elevation -- no need to clip NA.
            call elevation_data%get_xy(x,y, domain%U(:,j,ELV), domain%nx(1), bilinear=1_ip, &
                raster_index=domain%elevation_source_file_index(:,j))

            ! Here we can experiment with the random perturbation
            call random_number(random_uniform)
            domain%U(:,j,STG) = domain%U(:,j,STG) + random_perturbation_scale * (random_uniform-0.5_dp)
            call random_number(random_uniform)
            domain%U(:,j,ELV) = domain%U(:,j,ELV) + random_perturbation_scale * (random_uniform-0.5_dp)
        
        end do
        call elevation_data%finalise()
        call stage_data%finalise()

        if(domain%timestepping_method == 'cliffs') then
            call domain%smooth_elevation(smooth_method = 'cliffs')
        end if

        deallocate(x,y, random_uniform)

        ! Smooth near fine-to-coarse boundaries. Must do this BEFORE adjusting the stage to be >= elevation.
        if(present(all_dx_md)) call domain%smooth_elevation_near_nesting_fine2coarse_boundaries(all_dx_md)

        if(domain%timestepping_method /= 'linear') then
            domain%manning_squared = 0.02_dp * 0.02_dp
        else
            ! Avoid very shallow depths below MSL in linear solver
            where(domain%U(:,:,ELV) < 0.0_dp .and. domain%U(:,:,ELV) > -5.0_dp) domain%U(:,:,ELV) = -5.0_dp
        end if

        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV) + 1.0e-07_dp)

        write(log_output_unit,*) 'Stage range is: ', minval(domain%U(:,:,STG)), maxval(domain%U(:,:,STG))
        write(log_output_unit,*) 'Elev range is: ', minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))

    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program BP09
    !!
    !! NTHMP benchmark problem 9 -- Okushiri tsunami field test case
    !!

    use global_mod, only: ip, dp, charlen, default_nonlinear_timestepping_method,&
        default_linear_timestepping_method
    use multidomain_mod, only: multidomain_type
    use boundary_mod, only: flather_boundary
    use timer_mod, only: timer_type
    use logging_mod, only: log_output_unit, send_log_output_to_file
    use stop_mod, only: generic_stop
    use coarray_intrinsic_alternatives, only : swals_mpi_init, swals_mpi_finalize
    use iso_c_binding, only: C_DOUBLE
    use local_routines

    implicit none

    ! Type holding all domains 
    type(multidomain_type) :: md

    ! Local timing object
    type(timer_type) :: program_timer

    ! Increase mesh_refine to decrease the cell size (i.e. for convergence testing)
    real(dp), parameter :: mesh_refine = 0.4_dp ! Good for more efficient testing
    !real(dp), parameter :: mesh_refine = 1.0_dp ! Standard

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 7.50_dp
    real(dp), parameter :: final_time = 3600.0_dp * 1.0_dp

    ! Length/width
    real(dp), parameter :: global_lw(2) = [3.9_dp, 4.6_dp]
    ! Lower-left corner coordinate
    real(dp), parameter :: global_ll(2) = [137.6_dp, 39.6_dp]
    ! grid size (number of x/y cells). 
    integer(ip), parameter :: global_nx(2) = nint(global_lw*[77, 100]*1.12_dp*mesh_refine)

    ! Refinement factor for domains
    integer(ip), parameter :: nest_ratio = 5_ip

    ! The global (i.e. outer-domain) time-step in the multidomain 
    real(dp) ::  global_dt = 0.20_dp * (1.0_dp/mesh_refine)

    !! Optionally put a "very high res" domain around monai. Slows things down and the code requires care
    !! If this is true, then also use "mesh_refine = 1.0_dp"
    !logical, parameter:: very_high_res_monai = .false.
    !! If using very high res domain, reduce the entire model timestep when it becomes active
    !real(dp), parameter :: very_high_res_timestep_reduction = 10.0_dp
    !! If using a very high res domain, assume it is dry before this time. 
    !! If not set correctly, expect the wrong answer!
    !real(dp), parameter :: very_high_res_static_before_time = 235.0_dp

    ! Optionally smooth the elevation near fine 2 coarse nesting boundaries.
    logical, parameter :: smooth_elevation_near_fine2coarse_nesting_boundaries = .false.


    ! Useful misc variables
    integer(ip):: j, nd

    call swals_mpi_init()

    call program_timer%timer_start('setup')

#ifndef SPHERICAL
    write(log_output_unit,*) 'Code assumes spherical coordinates, but SPHERICAL is not defined'
    call generic_stop
#endif
    
    !! nd domains in this model
    !if(very_high_res_monai) then
    !    ! In this case we put a 30cm x 30cm cell domain around the monai inundation peak (with mesh_refine=1.0_dp)
    !    !
    !    ! It leads to a peak of > 30.0m with mesh_refine=1.0_dp (vs 31.7 obs, although the latter varies depending on which
    !    ! dataset is used).
    !    !
    !    ! The model needs various modifications to do this stably, and it takes more than twice as long.
    !    nd = 7
    !else
        ! This case has a domain res ~ 1.5x1.5m around the monai inundation peak (with mesh_refine = 1.0_dp)
        !
        ! It leads to a peak of >28.0m with mesh_refine=1.0_dp (vs 31.7 obs, although the latter varies depending on which
        ! dataset is used). So the very_high_res_monai case does better.
        nd = 6
    !end if

    allocate(md%domains(nd))

    ! Need to provide a load_balance_partition.txt file to get exact reproducibility 
    ! with varying number of openmp threads & MPI ranks
    call get_command_argument(1, md%load_balance_file)

    !
    ! Setup basic metadata
    !

    ! Linear domain
    md%domains(1)%lw = global_lw
    md%domains(1)%lower_left =global_ll
    md%domains(1)%nx = global_nx
    md%domains(1)%dx = md%domains(1)%lw/md%domains(1)%nx
    md%domains(1)%dx_refinement_factor = 1.0_dp
    md%domains(1)%timestepping_refinement_factor = 1_ip
    md%domains(1)%timestepping_method = default_linear_timestepping_method 
    md%domains(1)%cliffs_minimum_allowed_depth = 5.0_dp

    ! Higher res around region of interest
    call md%domains(2)%match_geometry_to_parent(&
        parent_domain=md%domains(1), &
        lower_left  = [139.2_dp, 41.5_dp], &
        upper_right = [140.5_dp, 43.5_dp], &
        dx_refinement_factor = nest_ratio, &
        timestepping_refinement_factor = 1_ip)
    md%domains(2)%timestepping_method = default_nonlinear_timestepping_method
    md%domains(2)%cliffs_minimum_allowed_depth = 2.0_dp

    ! Okushiri Island focus
    call md%domains(3)%match_geometry_to_parent(&
        parent_domain=md%domains(2), &
        lower_left  = [139.3_dp, 41.95_dp], &
        upper_right = [139.6_dp, 42.26_dp], &
        dx_refinement_factor = nest_ratio, &
        timestepping_refinement_factor = 2_ip)
    md%domains(3)%timestepping_method = default_nonlinear_timestepping_method
    md%domains(3)%cliffs_minimum_allowed_depth = 1.0_dp

    ! The Monai domain 
    call md%domains(4)%match_geometry_to_parent(&
        parent_domain=md%domains(3), &
        lower_left  = [139.38_dp, 42.0828_dp], &
        upper_right = [139.44_dp, 42.11592_dp], &
        dx_refinement_factor = nest_ratio, &
        timestepping_refinement_factor = 6_ip,&
        rounding_method='nearest')
    md%domains(4)%timestepping_method = default_nonlinear_timestepping_method  
    md%domains(4)%cliffs_minimum_allowed_depth = 1.0_dp
    
    ! A more detailed Monai domain 
    call md%domains(5)%match_geometry_to_parent(&
        parent_domain=md%domains(4), &
        lower_left  = [139.4150_dp, 42.0968_dp], &
        upper_right = [139.4264_dp, 42.1020_dp], &
        dx_refinement_factor = nest_ratio, &
        timestepping_refinement_factor = 6_ip,&
        rounding_method='nearest')
    md%domains(5)%timestepping_method = default_nonlinear_timestepping_method  
    md%domains(5)%cliffs_minimum_allowed_depth = 0.2_dp

    ! The Aonae domain
    call md%domains(6)%match_geometry_to_parent(&
        parent_domain=md%domains(3), &
        lower_left  = [139.440_dp, 42.030_dp], &
        upper_right = [139.5_dp, 42.070_dp], &
        dx_refinement_factor = nest_ratio, &
        timestepping_refinement_factor = 2_ip,&
        rounding_method = 'nearest')
    md%domains(6)%timestepping_method = default_nonlinear_timestepping_method  
    md%domains(6)%cliffs_minimum_allowed_depth = 0.2_dp

    !if(very_high_res_monai) then
    !    ! An even more detailed Monai domain
    !    call md%domains(7)%match_geometry_to_parent(&
    !        parent_domain=md%domains(5), &
    !        !lower_left  = [139.4176_dp, 42.0935_dp], &
    !        lower_left  = [139.4230_dp, 42.09896_dp], &
    !        !upper_right = [139.4299_dp, 42.1046_dp], &
    !        upper_right = [139.4254_dp, 42.09995_dp], &
    !        dx_refinement_factor = nest_ratio, &
    !        timestepping_refinement_factor = 6_ip,&
    !        rounding_method='nearest')
    !    ! Older versions of the code needed euler for stability.
    !    md%domains(7)%timestepping_method = 'euler'  
    !    ! Should be dry before this time -- note also we change the time-step 
    !    ! in the evolve loop on the basis of this value
    !    md%domains(7)%static_before_time = very_high_res_static_before_time
    !end if

    do j = 1, size(md%domains)
        md%domains(j)%nontemporal_grids_to_store = [character(len=charlen) :: &
            'max_stage', 'elevation0', 'manning_squared', 'elevation_source_file_index', &
            'time_of_max_stage', 'min_stage']
    end do

    ! Allocate domains and prepare comms
    call md%setup()

    call md%memory_summary()

    ! Give boundary condition to domans that have non-nesting boundaries
    do j = 1, size(md%domains)
        if(any(.not.md%domains(j)%is_nesting_boundary)) md%domains(j)%boundary_subroutine => flather_boundary
    end do

    ! Set initial conditions
    do j = 1, size(md%domains)
        if(smooth_elevation_near_fine2coarse_nesting_boundaries) then
            ! Option with local smoothing of elevation
            call set_initial_conditions_BP09(md%domains(j), md%all_dx_md)
        else
            ! No smoothing
            call set_initial_conditions_BP09(md%domains(j))
        end if
    end do
    call md%make_initial_conditions_consistent()

    ! For stability in 'null' regions, we set them to 'high land' that should be inactive. 
    call md%set_null_regions_to_dry()
   
    write(log_output_unit,*) 'End setup'

    ! Print the gravity-wave CFL limit, to guide timestepping
    do j = 1, size(md%domains)
        write(log_output_unit,*) 'domain: ', j, 'ts: ', &
            md%domains(j)%stationary_timestep_max()
    end do

    call program_timer%timer_end('setup')
    call program_timer%timer_start('evolve')

    flush(log_output_unit)

    !
    ! Evolve the code
    !

    do while (.true.)
       
        ! Print and write outputs each time-interval of "approximate_writeout_frequency"
        call md%write_outputs_and_print_statistics(&
            approximate_writeout_frequency=approximate_writeout_frequency, &
            write_grids_less_often = 1_ip, &
            write_gauges_less_often = 1_ip, &
            print_less_often = 1_ip,&
            timing_tol = 1.0e-06_dp)

        ! Finish looping at some point
        if (md%domains(1)%time > final_time) exit

        !if(very_high_res_monai) then
        !    ! Take a different time step once the high res domain starts evolving.
        !    ! This will introduce a (formal) first-order error into the linear leap-frog scheme,
        !    ! because that requires a fixed time step, but only at the change in time point.
        !    ! But practically no problem (?).
        !    if(md%domains(1)%time < md%domains(7)%static_before_time) then
        !        call md%evolve_one_step(global_dt)
        !    else
        !        call md%evolve_one_step(global_dt/very_high_res_timestep_reduction)
        !    end if
        !else
            ! Regular case
            call md%evolve_one_step(global_dt)
        !end if

    end do

    call program_timer%timer_end('evolve')
    call md%finalise_and_print_timers

    write(log_output_unit,*) ''
    write(log_output_unit, *) 'Program timer'
    write(log_output_unit, *) ''
    call program_timer%print(log_output_unit)

    call swals_mpi_finalize()
end program
