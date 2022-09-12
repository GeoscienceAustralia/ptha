module local_routines 
    !!
    !! Setup NTHMP currents test problem -- Tohoku tsunami at Hilo Harbour with an artificial forcing.
    !!
    use global_mod, only: dp, ip, charlen, wall_elevation
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use read_raster_mod, only: multi_raster_type
    use file_io_mod, only: count_file_lines
    use linear_interpolator_mod, only: linear_interpolator_type
    use points_in_poly_mod, only: points_in_poly

    implicit none

    ! Hold some data for the boundary condition
    type :: boundary_information_type
        character(charlen):: bc_file
        real(dp), allocatable:: boundary_data(:,:)
        type(linear_interpolator_type):: gauge_ts_function
        real(dp):: boundary_elev
        real(dp):: t0 = 0.0_dp
    end type

    ! Parameter controlling the northern boundary treatment
    character(len=charlen) :: boundary_type = 'boundary_stage_radiation_momentum'

    ! This will hold the information -- is seen by other parts of the module
    type(boundary_information_type):: boundary_information

    contains 

    ! Read files with boundary condition info, and make a BC function
    subroutine setup_boundary_information(bc_file, boundary_elev)
        character(charlen), intent(in):: bc_file
        real(dp), intent(in):: boundary_elev

        integer(ip):: bc_unit, nr, nc, skip, i, extra

        boundary_information%bc_file = bc_file
        boundary_information%boundary_elev = boundary_elev
        open(newunit=bc_unit, file=bc_file)
        nr = count_file_lines(bc_unit)
        nc = 2
        skip = 1
        extra = 1 ! One more data point to avoid exceeding time 
        allocate(boundary_information%boundary_data(nr - skip + extra, nc))
        do i = 1, nr
            if(i > skip) then
                read(bc_unit, *) boundary_information%boundary_data(i - skip,:)
            else
                read(bc_unit, *) 
            end if 
        end do
        close(bc_unit)
        ! Extend the time-series with a constant value, so that time does not
        ! exceed model run time
        boundary_information%boundary_data(nr - skip + extra,1) = 1.0e+06_dp + &
            boundary_information%boundary_data(nr - skip + extra - 1, 1)

        boundary_information%boundary_data(nr - skip + extra,2) = &
            boundary_information%boundary_data(nr - skip + extra - 1, 2)

        call boundary_information%gauge_ts_function%initialise(&
                boundary_information%boundary_data(:,1), &
                boundary_information%boundary_data(:,2))
        boundary_information%t0 = boundary_information%boundary_data(1,1)

    end subroutine
    
    ! Make a function to evaluate the boundary at the domain
    !
    function boundary_function(domain, t, i, j) result(stage_uh_vh_elev)
        type(domain_type), intent(in):: domain
        real(dp), intent(in):: t
        integer(ip), intent(in) :: i, j
        real(dp):: stage_uh_vh_elev(4)
        real(dp) :: local_elev, dhdt(1), dt, next_h(1)

        ! Set the stage
        call boundary_information%gauge_ts_function%eval(&
            [t + boundary_information%t0], stage_uh_vh_elev(1:1))
      
        ! Get the time-derivative of stage. Useful for some approaches
        dt = 1.0e-06
        call boundary_information%gauge_ts_function%eval(&
            [t + boundary_information%t0 + dt], next_h)
        dhdt = (next_h(1)/dt - stage_uh_vh_elev(1)/dt)
       
        ! Set the elevation 
        local_elev = domain%U(i,j,4)
        stage_uh_vh_elev(4) = local_elev

        if(local_elev >= stage_uh_vh_elev(1)) then
            ! Treat dry boundary case
            stage_uh_vh_elev(1) = local_elev
            stage_uh_vh_elev(2:3) = 0.0_dp
        else

            select case(boundary_type)

            case('boundary_stage_radiation_momentum') 
                ! These will never be used for this boundary
                stage_uh_vh_elev(2:3) = 0.0_dp

            case('flather_with_vh_from_continuity')
                ! Approach 4: Flat-free-surface continuity
                if(j == domain%nx(2)) then
                    stage_uh_vh_elev(2) = 0.0_dp
                    ! Assume flat free surface from offshore to the coast + 
                    ! estuary volume. That gives a rough estimate of -VH. 
                    ! However, in practice the factor needs tuning.
                    stage_uh_vh_elev(3) = 0.0_dp -dhdt(1) * 3800.0_dp
                else
                    stage_uh_vh_elev(2:3) = 0.0_dp
                end if

            case default
                print*, 'boundary_type not recognised'
                stop
            end select

        end if

    end function

    ! Main setup routine
    subroutine set_initial_conditions(domain)
        class(domain_type), target, intent(inout):: domain
        integer(ip):: j
        character(len=charlen):: input_elevation, input_stage
        real(dp), allocatable:: x(:), y(:)
        type(multi_raster_type):: elevation_data
        real(dp) :: wall

        ! Stage
        domain%U(:,:,STG) = 0.0_dp

        ! Set elevation with the raster
        input_elevation = './bathymetry/hilo_grid_1_3_arsec.tif' 
        call elevation_data%initialise([input_elevation])

        ! Make space for x/y coordinates, at which we will look-up the rasters
        allocate(x(domain%nx(1)), y(domain%nx(1)))
        x = domain%x

        do j = 1, domain%nx(2)
            y = domain%y(j)
            call elevation_data%get_xy(x, y, domain%U(:,j,ELV), domain%nx(1), &
                bilinear=1_ip)
        end do
        deallocate(x,y)

        if(domain%timestepping_method == 'cliffs') then
            domain%cliffs_minimum_allowed_depth = 0.1_dp
            call domain%smooth_elevation(smooth_method='cliffs')
        end if

        print*, 'Elevation range: ', minval(domain%U(:,:,ELV)), &
            maxval(domain%U(:,:,ELV))

        ! A bit unclear how this boundary should be treated -- a wall works ok
        wall = 20._dp
        domain%U((domain%nx(1)-1):domain%nx(1),:,ELV) = wall

        ! Friction 
        if(domain%timestepping_method /= 'linear') then
            domain%manning_squared = 0.025_dp**2
        end if

        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV) + 1.0e-07_dp)

    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program Hilo_harbour_Tohoku
    !!
    !! NTHMP currents test problem -- Tohoku tsunami at Hilo Harbour with an artificial forcing.
    !!

    use global_mod, only: ip, dp, minimum_allowed_depth, &
        default_nonlinear_timestepping_method
    use domain_mod, only: domain_type
    use multidomain_mod, only: multidomain_type, setup_multidomain
    use boundary_mod, only: flather_boundary, boundary_stage_radiation_momentum
    use local_routines
    use timer_mod
    use logging_mod, only: log_output_unit
    use coarray_intrinsic_alternatives, only: swals_mpi_init, swals_mpi_finalize
    implicit none

    ! Useful misc variables
    integer(ip):: j, nd

    ! Type holding all domains 
    type(multidomain_type) :: md

    type(timer_type) :: program_timer

    ! Increase resolution by this amount, e.g. 1.0 = 1/3 arcsec, 2.0 = 1/6 arcsec
    real(dp), parameter :: mesh_refine = 0.5_dp 
    
    real(dp) ::  global_dt = 0.27_dp / mesh_refine
    real(dp), parameter :: final_time = 23370.0_dp 

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 30.0_dp
    integer(ip), parameter :: only_write_grids_every_nth_output_step = 1_ip

    character(len=charlen) ::  bc_file = './boundary/se_dat_converted.csv'
    real(dp) :: bc_elev

    ! The elevation has resolution = 1/3 arc-sec. Make sure our model is inside this
    real(dp), parameter :: onethird_arcsec = 1.0_dp / (3.0_dp * (60.0_dp**2) )
    ! Length/width
    real(dp), parameter :: global_lw(2) = [0.065_dp, 0.05_dp] - (4 * onethird_arcsec) 
    ! Lower-left corner coordinate
    real(dp), parameter :: global_ll(2) = [204.9_dp, 19.71_dp] + 2*onethird_arcsec
    ! grid size (number of x/y cells)
    integer(ip), parameter :: global_nx(2) = nint(global_lw/onethird_arcsec * mesh_refine) 

    call swals_mpi_init
    call program_timer%timer_start('setup')

    ! nd domains in this model
    nd = 1
    allocate(md%domains(nd))
    md%load_balance_file = 'load_balance_partition.txt'

    !
    ! Setup basic metadata
    !

    md%domains(1)%lower_left = global_ll
    md%domains(1)%nx = global_nx
    md%domains(1)%lw = global_lw
    md%domains(1)%dx = md%domains(1)%lw/md%domains(1)%nx
    md%domains(1)%timestepping_method = default_nonlinear_timestepping_method
    md%domains(1)%static_before_time = 1920.0_dp ! No need to evolve prior to this time

    print*, 1, ' lw: ', md%domains(1)%lw, ' ll: ', md%domains(1)%lower_left, &
        ' dx: ', md%domains(1)%dx, ' nx: ', md%domains(1)%nx

    ! Allocate domains and prepare comms
    call md%setup()

    if(md%domains(1)%is_staggered_grid) then
        ! The 'boundary_stage_radiation_momentum' boundary is insufficiently 
        ! radiative with the leapfrog solvers. This bc is still not great, but better.
        boundary_type = 'flather_with_vh_from_continuity'
    else
        ! Good for rk2 and related type solvers -- but poor for 
        ! leapfrog_nonlinear, and not so good for cliffs.
        boundary_type = 'boundary_stage_radiation_momentum'
    end if

    ! Initial conditions
    do j = 1, size(md%domains)
        call set_initial_conditions(md%domains(j))
    end do
    call md%make_initial_conditions_consistent 
    call md%set_null_regions_to_dry()

    ! The elevation has a flat area near the boundary (see NTHMP problem blurb for discussion).
    ! Get min elevation there as we use it in the boundary forcing.
    bc_elev = HUGE(1.0_dp)
    do j = 1, size(md%domains)
        bc_elev = min(bc_elev, minval(md%domains(j)%U(:,:,4)))
    end do
    ! Build boundary conditions
    call setup_boundary_information(bc_file, bc_elev)

    do j = 1, size(md%domains) 
        select case(boundary_type)
        case('flather_with_vh_from_continuity')
            md%domains(j)%boundary_subroutine => flather_boundary
        case('boundary_stage_radiation_momentum')
            md%domains(j)%boundary_subroutine => boundary_stage_radiation_momentum
        case default
            stop "Invalid boundary_type value"
        end select
        md%domains(j)%boundary_function => boundary_function
    end do

    ! Setup hazard points
    call md%set_point_gauges_from_csv("point_gauges.csv", skip_header=1_ip)
   
    ! Print the gravity-wave CFL limit, to guide timestepping
    do j = 1, size(md%domains)
        write(log_output_unit, *) 'domain: ', j, 'ts: ', &
            md%domains(j)%stationary_timestep_max()
    end do

    print*, 'End setup'
    call program_timer%timer_end('setup')
    call program_timer%timer_start('evolve')

    ! Evolve the code
    do while (.true.)
       
        call program_timer%timer_start('IO')
        call md%write_outputs_and_print_statistics(&
            approximate_writeout_frequency=approximate_writeout_frequency, &
            write_grids_less_often = only_write_grids_every_nth_output_step, &
            write_gauges_less_often = 1_ip, &
            print_less_often = 1_ip, &
            timing_tol = 1.0e-06_dp)
        call program_timer%timer_end('IO')

        ! Finish at some point
        if (md%domains(1)%time > final_time) exit

        ! Main evolve
        call md%evolve_one_step(global_dt)

    end do

    call program_timer%timer_end('evolve')
    call md%finalise_and_print_timers

    write(log_output_unit,*) ''
    write(log_output_unit, *) 'Program timer'
    write(log_output_unit, *) ''
    call program_timer%print(log_output_unit)

    call swals_mpi_finalize
end program
