module local_routines 
    !!
    !! Setup the undular bore problem
    !!

    use global_mod, only: dp, ip, charlen, pi
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use read_raster_mod, only: read_gdal_raster
    use which_mod, only: which
    use file_io_mod, only: count_file_lines
    use linear_interpolator_mod, only: linear_interpolator_type
    implicit none

    real(dp), parameter :: mean_elevation = -20.0_dp, length = 30000.0_dp
    real(dp), parameter :: amplitude = 2.0_dp, period = 780.0_dp

    contains

    function boundary_function(domain, t, i, j) result(stage_uh_vh_elev)
        ! Function to evaluate the boundary at the domain, passed to model boundary conditions
        type(domain_type), intent(in):: domain
        real(dp), intent(in):: t
        integer(ip), intent(in) :: i, j
        real(dp):: stage_uh_vh_elev(4)

        if(i == 1) then
            stage_uh_vh_elev (1) = amplitude * sin(2.0_dp*pi*t/period)
            stage_uh_vh_elev(2) = 0.0_dp ! Comes from boundary condition
            stage_uh_vh_elev(3) = 0.0_dp
            stage_uh_vh_elev(4) = domain%U(i,j,ELV)
        else
            stage_uh_vh_elev(STG) = domain%U(i,j,STG)
            stage_uh_vh_elev(ELV) = domain%U(i,j,ELV)
            stage_uh_vh_elev(UH:VH) = 0.0_dp
        end if
    end function

    subroutine set_initial_conditions(domain)
        type(domain_type), intent(inout):: domain

        real(dp):: x, y, elev
        integer(ip):: j, i, k
        real(dp):: gauge_xy(2,11), wall

        ! Stage
        domain%U(:,:,STG) = 0.0_dp
        domain%msl_linear = 0.0_dp

        domain%U(:,:,UH:VH) = 0.0_dp

        domain%U(:,:,ELV) = mean_elevation

        ! Reflective boundaries on 3 sides
        wall = 10.0_dp
        domain%U(:, 1, ELV) = wall
        domain%U(:, 2, ELV) = wall
        domain%U(:, domain%nx(2), ELV) = wall
        domain%U(:, domain%nx(2)-1, ELV) = wall
        domain%U(domain%nx(1), :, ELV) = wall
        domain%U(domain%nx(1)-1, :, ELV) = wall
   
        ! Stage >= bed 
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV))
        
        !! Define locations of gauge outputs
        !gauge_xy(2,:) = 0.0_dp ! Always y == 0
        !gauge_xy(1, 1:4) = 0.0_dp + 1.0e-06_dp ! Nudge x-coordinate inside the domain
        !gauge_xy(1,5) = tank_x(2)
        !gauge_xy(1,6) = 0.5_dp * (tank_x(3) + tank_x(2))
        !gauge_xy(1,7) = tank_x(3)
        !gauge_xy(1,8) = 0.5_dp * (tank_x(3) + tank_x(4))
        !gauge_xy(1,9) = tank_x(4)
        !gauge_xy(1,10) = 0.5_dp * (tank_x(4) + tank_x(5))
        !! Include just before the wall
        !gauge_xy(1,11) = domain%x(domain%nx(1) - 2)

        !call domain%setup_point_gauges(gauge_xy)

    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program undular_bore
    !!
    !! Undular bore development from an initial sinusoidal wave
    !!
    use global_mod, only: ip, dp, minimum_allowed_depth
    use multidomain_mod, only: multidomain_type
    use boundary_mod, only: boundary_stage_transmissive_momentum, flather_boundary, &
        transmissive_boundary, boundary_stage_transmissive_normal_momentum
    use linear_interpolator_mod, only: linear_interpolator_type
    use local_routines
    implicit none

    type(multidomain_type):: md

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 10.0_dp
    real(dp), parameter :: final_time = 2000.0_dp

    ! Domain info
    character(charlen) :: timestepping_method, buf

    ! Local variables 
    real(dp) :: dx_in, dx, timestep
    character(charlen):: test_case, bc_file
    real(dp):: tank_bases(4), tank_slopes(4) 


    ! Timestepping method
    call get_command_argument(1, timestepping_method)  
    ! Resolution
    call get_command_argument(2, buf)
    read(buf, *) dx_in
      
    
    ! Resolution
    dx = 1.0_dp*dx_in !5.0_dp !10.0_dp

    ! Setup model with 1 domain
    allocate(md%domains(1))
    md%domains(1)%lw = [length, 5*dx]
    md%domains(1)%lower_left =  [0.0_dp, -2.5_dp * dx]
    md%domains(1)%nx = nint(md%domains(1)%lw/dx)
    md%domains(1)%timestepping_method = timestepping_method
    md%domains(1)%use_dispersion = .true. !
    !md%domains(1)%nc_grid_output%flush_every_n_output_steps = 1_ip !

    ! Tapering off of dispersive terms (won't have any effect)
    md%domains(1)%ds%td1 = 1.0_dp
    md%domains(1)%ds%td2 = 0.5_dp

    ! Non-TVD limiting for finite volume schemes
    md%domains(1)%theta = 4.0_dp

    ! Output directory should record the solver and resolution
    md%output_basedir = 'OUTPUTS/' // trim(timestepping_method) // '_' // trim(buf)

    call md%setup

    call set_initial_conditions(md%domains(1))
    md%domains(1)%boundary_function => boundary_function
    md%domains(1)%boundary_subroutine => boundary_stage_transmissive_normal_momentum

    call md%make_initial_conditions_consistent() ! Get the initial volume right

    ! Fixed timestep  
    timestep = md%stationary_timestep_max() * 0.5_dp 
    !print*, trim(timestepping_method), ', dx = ', dx, ', timestep = ', timestep

    ! Evolve the code
    do while (.true.)

        ! Avoid storing grids often
        call md%write_outputs_and_print_statistics(&
            approximate_writeout_frequency=approximate_writeout_frequency, &
            write_grids_less_often = 1_ip)

        if (md%domains(1)%time > final_time) exit

        call md%evolve_one_step(timestep)

    end do

    call md%finalise_and_print_timers

end program
