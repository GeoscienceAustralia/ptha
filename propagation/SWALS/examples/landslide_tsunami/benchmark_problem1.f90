module local_routines 
    !!
    !! Setup routine for "Benchmark Problem 1" from the "Third International workshop
    !! on long-wave runup models, June 17-18 2004". It models the runup of an initial
    !! waveform on a linearly sloping beach. The analytical solution was produced
    !! using the techniques of Carrier, Wu and Yeh (2003). The problem descriptions
    !! and solution were sourced from http://isec.nacse.org/workshop/2004\_cornell/bmark1.html
    !!
    use global_mod, only: dp, ip, charlen, wall_elevation
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use linear_interpolator_mod, only: linear_interpolator_type

    implicit none

    contains 

    subroutine set_initial_conditions(domain, wall_width)
        class(domain_type), target, intent(inout):: domain
        real(dp), intent(in) :: wall_width

        integer(ip):: i, j, initial_stage_unit, file_header_lines
        character(len=charlen):: initial_stage
        real(dp), allocatable:: x(:)
        real(dp) :: wall
        real(dp) :: gauge_xy(3,3)
        real(dp), parameter :: beach_slope = -1.0_dp/10.0_dp
        type(linear_interpolator_type):: initial_stage_interpolator
        integer(ip), parameter :: nb = 1001 ! Number of points in boundary series
        real(dp) :: initial_x(nb + 1), initial_s(nb + 1)

        !
        ! Read the initial stage
        !
        initial_stage = './DATA/initial_condition.txt'
        open(newunit=initial_stage_unit, file=initial_stage)
        ! Skip the file header
        file_header_lines=13
        do i = 1, file_header_lines
            read(initial_stage_unit, *)
        end do
        ! Since the file starts at x=0, lets append large negative x, so the
        ! interpolator covers the whole domain
        initial_x(1) = -1.0e+06_dp
        initial_s(1) = 0.0_dp
        do i = 2, nb+1
            read(initial_stage_unit, *) initial_x(i), initial_s(i)
        end do
        close(initial_stage_unit)
        ! Make the interpolation function
        call initial_stage_interpolator%initialise(initial_x, initial_s)
        
        allocate(x(domain%nx(1)))
        x = domain%x
        wall = 50.0_dp
        do j = 1, domain%nx(2)
            ! Make the topography, with side-walls of fixed width (cannot just
            ! have a single cell, because of nesting)
            if( (domain%y(j) < wall_width) .or. (domain%lw(2) - domain%y(j) < wall_width)) then
                domain%U(:,j,ELV) = wall
            else
                domain%U(:,j,ELV) = x * beach_slope
            end if
            call initial_stage_interpolator%eval(x, domain%U(:,j,STG))
        end do
        call initial_stage_interpolator%finalise

        !if(maxval(domain%x) > 49950.0_dp) then
        !    domain%U(domain%nx(1), :, ELV) = wall
        !end if
        ! Use walls as lateral boundary conditions
        ! FIXME: Use periodic boundary conditions instead
        !domain%U(domain%nx(1), :, ELV) = wall
        
        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV) + 1.0e-07_dp)

        deallocate(x)

        print*, 'Elevation range: ', minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))

        if(allocated(domain%manning_squared)) then
            domain%manning_squared = 0.0_dp**2
        end if

        if(domain%timestepping_method == 'cliffs') then
            domain%cliffs_minimum_allowed_depth = 0.02_dp
        end if

        !! Gauges
        !gauge_xy(1:3, 1) = [4.521, 1.196, 5.0]
        !gauge_xy(1:3, 2) = [4.521, 1.696, 7.0]
        !gauge_xy(1:3, 3) = [4.521, 2.196, 9.0]
        !call domain%setup_point_gauges(xy_coords = gauge_xy(1:2,:), gauge_ids=gauge_xy(3,:))

    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program landslide_tsunami 
    !! "Benchmark Problem 1" from the "Third International workshop
    !! on long-wave runup models, June 17-18 2004". It models the runup of an initial
    !! waveform on a linearly sloping beach. The analytical solution was produced
    !! using the techniques of Carrier, Wu and Yeh (2003). The problem descriptions
    !! and solution were sourced from http://isec.nacse.org/workshop/2004\_cornell/bmark1.html

    use global_mod, only: ip, dp, minimum_allowed_depth, default_nonlinear_timestepping_method
    use domain_mod, only: domain_type
    use multidomain_mod, only: multidomain_type, setup_multidomain, test_multidomain_mod
    use boundary_mod, only: boundary_stage_transmissive_normal_momentum, flather_boundary
    use local_routines
    use timer_mod
    use logging_mod, only: log_output_unit
    implicit none

    ! Useful misc variables
    integer(ip):: j, i, i0, j0, centoff, nd, lg
    real(dp):: last_write_time, gx(4), gy(4)

    ! Type holding all domains 
    type(multidomain_type) :: md

    type(timer_type) :: program_timer

    real(dp), parameter :: mesh_refine = 1.0_dp ! Increase/decrease resolution by this amount
    
    real(dp) ::  global_dt != 0.040_dp / mesh_refine

    ! Approx timestep between outputs
    real(dp) :: approximate_writeout_frequency = 1.00_dp
    real(dp) :: final_time = 360.0_dp

    ! Length/width
    real(dp), parameter :: global_lw(2) = [400.0_dp + 50000.0_dp, 125.0_dp/mesh_refine]
    ! Lower-left corner coordinate
    real(dp), parameter :: global_ll(2) = [-400.0_dp, 0.0_dp]
    ! grid size (number of x/y cells)
    !integer(ip), dimension(2):: global_nx = [nint(global_lw(1)/25.0_dp)*nint(mesh_refine), 5_ip] !nint(global_lw/25.0_dp) * mesh_refine

    real(dp) :: res_d1 = 25.0_dp/mesh_refine, res_d2 = 5.0_dp/mesh_refine

    call program_timer%timer_start('setup')

    ! nd domains in this model
    nd = 2
    allocate(md%domains(nd))

    !
    ! Setup basic metadata
    !

    ! Main domain
    md%domains(1)%lower_left = global_ll  + [3000.0_dp, 0.0_dp] !global_ll
    md%domains(1)%lw = global_lw - [3000.0_dp, 0.0_dp]
    md%domains(1)%nx = nint(md%domains(1)%lw/res_d1)
    md%domains(1)%dx = md%domains(1)%lw/md%domains(1)%nx
    md%domains(1)%timestepping_refinement_factor = 1_ip
    md%domains(1)%dx_refinement_factor = 1.0_dp
    md%domains(1)%timestepping_method = default_nonlinear_timestepping_method ! Can set this to 'linear', but the difference with the analytical solution becomes obvious

    !print*, 1, ' lw: ', md%domains(1)%lw, ' ll: ', md%domains(1)%lower_left, ' dx: ', md%domains(1)%dx, &
    !    ' nx: ', md%domains(1)%nx

    ! Main domain
    md%domains(2)%lower_left = global_ll 
    md%domains(2)%lw = global_lw - [md%domains(1)%lw(1), 0.0_dp]
    md%domains(2)%nx = nint(md%domains(2)%lw/res_d2)
    md%domains(2)%dx = md%domains(2)%lw/md%domains(2)%nx
    md%domains(2)%timestepping_refinement_factor = 1_ip
    md%domains(2)%dx_refinement_factor = res_d1/res_d2
    md%domains(2)%timestepping_method = default_nonlinear_timestepping_method !'cliffs'
    
    ! Allocate domains and prepare comms
    call md%setup()

    ! Initial conditions
    do j = 1, size(md%domains)
        call set_initial_conditions(md%domains(j), wall_width=res_d1*2)
    end do
    call md%make_initial_conditions_consistent()

    md%domains(1)%boundary_subroutine => flather_boundary !boundary_subroutine

    ! NOTE: For stability in 'null' regions, we set them to 'high land' that
    ! should be inactive. 
    call md%set_null_regions_to_dry()

    print*, 'End setup'

    ! Print the gravity-wave CFL limit, to guide timestepping
    do j = 1, size(md%domains)
        print*, 'domain: ', j, 'ts: ', &
            md%domains(j)%stationary_timestep_max()
    end do

    global_dt = md%stationary_timestep_max()
    ! rk2n needs a shorter timestep
    if(md%domains(1)%timestepping_method == 'rk2n') global_dt = global_dt * 0.8_dp

    ! Trick to get the code to write out just after the first timestep
    last_write_time = -approximate_writeout_frequency

    call program_timer%timer_end('setup')
    call program_timer%timer_start('evolve')

    ! Trick to get the code to write out at the first timestep
    last_write_time = -approximate_writeout_frequency
    ! Evolve the code
    do while (.true.)
        
        ! IO 
        call program_timer%timer_start('IO')
        call md%write_outputs_and_print_statistics(approximate_writeout_frequency=approximate_writeout_frequency)
        call program_timer%timer_end('IO')

        if (md%domains(1)%time > final_time) exit
        call md%evolve_one_step(global_dt)

    end do

    call program_timer%timer_end('evolve')
    call md%finalise_and_print_timers

    print*, ''
    call program_timer%print(output_file_unit=log_output_unit)

end program
