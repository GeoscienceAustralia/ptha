module local_routines 
    !! Test flow algorithms in a large area spherical box with complex topography
    use global_mod, only: dp, ip, wall_elevation, pi, charlen
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use read_raster_mod, only: multi_raster_type
    use logging_mod, only: log_output_unit
    implicit none

    contains 

    subroutine set_initial_conditions(domain)            
        class(domain_type), intent(inout):: domain
        integer(ip):: i, j
        real(dp), allocatable:: x(:), y(:)
        type(multi_raster_type) :: elevation_rast

        ! Make space for x/y coordinates, at which we will look-up the rasters
        allocate(x(domain%nx(1)), y(domain%nx(1)))
        x = domain%x

        call elevation_rast%initialise([character(len=charlen):: '../generic_example/japan_dem.tif'])
        
        ! Set stage and elevation row-by-row.
        do j = 1, domain%nx(2)
            y = domain%y(j)

            domain%U(:,j,STG) = exp(sin(2*pi*x/10.0_dp))*cos(2*pi*y/10.0_dp)

            call elevation_rast%get_xy(x, y, domain%U(:,j,ELV), domain%nx(1), bilinear=1_ip)
        end do

        deallocate(x,y)
       
        ! Wet everywhere and not too shallow
        domain%U(:,:,ELV) = min(domain%U(:,:,ELV), -100.0_dp)

        ! Ensure topography is smooth
        do j = 1,10
            call domain%smooth_elevation
        end do

        ! Wall boundaries
        domain%U(1:2, :, ELV) = 100.0_dp
        domain%U(:, 1:2, ELV) = 100.0_dp
        j = domain%nx(1)
        domain%U(j-1:j, :, ELV) = 100.0_dp
        j = domain%nx(2)
        domain%U(:, j-1:j, ELV) = 100.0_dp


        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV) + 1.0e-07_dp)

        write(log_output_unit,*) 'Stage range is: ', minval(domain%U(:,:,STG)), maxval(domain%U(:,:,STG))
        write(log_output_unit,*) 'Elev range is: ', minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))

    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program spherical_box
    !! Test flow algorithms in a large area spherical box with complex topography

    use global_mod, only: ip, dp, charlen
    use multidomain_mod, only: multidomain_type
    use timer_mod, only: timer_type
    use logging_mod, only: log_output_unit
    use stop_mod, only: generic_stop
    use local_routines
    implicit none

    ! Type holding all domains 
    type(multidomain_type) :: md

    ! Local timing object
    type(timer_type) :: program_timer

    character(len=charlen) :: timestepping_method

    ! Approx timestep between outputs
    real(dp) :: approximate_writeout_frequency = 1000
    real(dp) :: final_time = 3600.0_dp 
    real(dp) :: my_dt 

    ! Lower-left corner coordinate
    real(dp), parameter :: global_ll(2) = [122.0_dp, 14.0_dp] ![0.0_dp, 0.0_dp]
    ! Length/width
    real(dp), parameter :: global_lw(2) = [185.0_dp, 55.0_dp] - global_ll
    ! grid size (number of x/y cells)
    integer(ip), parameter :: global_nx(2) = [945_ip, 615_ip] 

    integer(ip):: j

    ! Set the model name
    md%output_basedir = './OUTPUTS/'

    call program_timer%timer_start('setup')

    call get_command_argument(1, timestepping_method)

#ifndef SPHERICAL
    write(log_output_unit,*) 'Code assumes spherical coordinates, must compile with -DSPHERICAL'
    call generic_stop
#endif

    ! Single domain model
    allocate(md%domains(1))

    ! Setup basic metadata
    md%domains(1)%lw = global_lw
    md%domains(1)%lower_left =global_ll
    md%domains(1)%nx = global_nx
    md%domains(1)%dx = md%domains(1)%lw/md%domains(1)%nx
    md%domains(1)%dx_refinement_factor = 1.0_dp
    md%domains(1)%timestepping_refinement_factor = 1_ip
    md%domains(1)%timestepping_method = timestepping_method

    ! Allocate domains and prepare comms
    call md%setup()

    ! Set initial conditions
    do j = 1, size(md%domains)
        call set_initial_conditions(md%domains(j))
    end do

    ! These steps are important in complex nested models, but not really needed here.
    call md%make_initial_conditions_consistent()
    call md%set_null_regions_to_dry()

    ! Fixed timestep
    my_dt = md%stationary_timestep_max()
   
    write(log_output_unit,*) 'End setup'

    call program_timer%timer_end('setup')
    call program_timer%timer_start('evolve')

    !
    ! Evolve the code
    !
    do while (.true.)

        call md%write_outputs_and_print_statistics(&
            approximate_writeout_frequency = approximate_writeout_frequency, &
            timing_tol=my_dt/3.0_dp)

        if (md%domains(1)%time > final_time) exit

        call md%evolve_one_step(my_dt)

    end do

    ! Write on the last time
    call md%write_outputs_and_print_statistics()

    call program_timer%timer_end('evolve')
    call md%finalise_and_print_timers

    write(log_output_unit,*) ''
    write(log_output_unit, *) 'Program timer'
    write(log_output_unit, *) ''
    call program_timer%print(log_output_unit)
end program
