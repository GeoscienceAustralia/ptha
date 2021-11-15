module local_routines 
    !!
    !! Setup overbank flow problem (aligned to grid). This can treat either NS-aligned or EW-aligned channels.
    !!
    use global_mod, only: dp, ip, wall_elevation, gravity
    use domain_mod, only: domain_type, STG, UH, VH, ELV 
    use which_mod, only: which
    implicit none
   
    ! Geometry and flow parameters 
    real(dp), parameter :: &
        fraction_floodplain = 0.4_dp, &
        fraction_banks = 0.4_dp, &
        fraction_flatbed = 0.2_dp, &
        floodplain_elev = 0.0_dp, &
        flatbed_elev = -5.0_dp, &
        bed_slope = 0.001_dp, &
        lambda = 10.0_dp, &
        manning_n = 0.02_dp, &
        initial_free_surface = 0.5_dp

    logical :: channel_is_NS_aligned

    contains 

    subroutine set_initial_flow_at_ij(domain, i, j)
        ! Useful worker subroutine -- set an approximate initial condition at [x,y]= [domain%x(i), domain%y(j)]
        ! We will also use this to set the boundary conditions
        type(domain_type), intent(inout) :: domain
        integer(ip), intent(in) :: i, j

        integer(ip):: UH_ind, VH_ind
        real(dp):: x, y, bank_gradient, depth, fricf, width
        real(dp), parameter :: wall_elevation = 20.0_dp
        logical :: is_a_side

        if(channel_is_NS_aligned) then
            ! NS aligned channel, flowing to the south
            bank_gradient = (floodplain_elev - flatbed_elev)/ &
                (domain%lw(1)*(1.0_dp - fraction_floodplain - fraction_flatbed)/2.0_dp)

            x = domain%x(i)
            y = domain%y(j)
            width = domain%lw(1)
            is_a_side = ( i == 1) .or. ( i == domain%nx(1) )
            UH_ind = UH
            VH_ind = VH
        else
            ! EW aligned channel, flowing to the east
            bank_gradient = (floodplain_elev - flatbed_elev)/ &
                (domain%lw(2)*(1.0_dp - fraction_floodplain - fraction_flatbed)/2.0_dp)

            y = domain%x(i)
            x = domain%y(j)
            width = domain%lw(2)
            is_a_side = ( j == 1) .or. ( j == domain%nx(2) )
            UH_ind = VH
            VH_ind = UH
        end if

        ! Trapezoidal geometry
        domain%U(i,j,ELV) = y * bed_slope + &
            max(flatbed_elev, min( &
               floodplain_elev, &
               flatbed_elev + bank_gradient * (abs(x) - width*fraction_flatbed/2)))

        ! Wall boundaries along the sides
        if( is_a_side ) domain%U(i,j,ELV) = domain%U(i,j,ELV) + wall_elevation

        ! Stage                       
        domain%U(i, j, STG) = max(y * bed_slope + initial_free_surface, domain%U(i,j,ELV))

        domain%U(i,j,UH_ind) = 0.0_dp

        ! Guess the steady-uniform flow solution
        depth = max(domain%U(i,j,STG) - domain%U(i,j,ELV), 0.0_dp)
        fricf = manning_n**2 * 8 * gravity * merge(depth**(-1.0_dp/3.0_dp), 999999.9_dp, depth > 0)
        domain%U(i,j,VH_ind) = -sqrt(gravity * depth * bed_slope * 8 / fricf) * depth

    end subroutine

    subroutine set_initial_conditions_overbank_flow(domain)            
        type(domain_type), intent(inout):: domain

        integer(ip):: j, i

        domain%manning_squared = manning_n**2

        do j = 1, domain%nx(2)
            do i = 1, domain%nx(1)
                call set_initial_flow_at_ij(domain, i, j)
            end do
        end do

    end subroutine

    subroutine boundary_forcing(domain)
        ! For the boundary condition, we use the approximate analytical solution
        type(domain_type), intent(inout) :: domain

        integer :: i

        if(channel_is_NS_aligned) then
            ! North-South flowing channel
            ! Set north and south boundary conditions
            do i = 1, domain%nx(1)
                call set_initial_flow_at_ij(domain, i, 1_ip)
            end do
            do i = 1, domain%nx(1)
                call set_initial_flow_at_ij(domain, i, domain%nx(2))
            end do
        else
            ! West-East flowing channel
            ! Set east and west boundary conditions
            do i = 1, domain%nx(2)
                call set_initial_flow_at_ij(domain, 1_ip, i)
            end do
            do i = 1, domain%nx(2)
                call set_initial_flow_at_ij(domain, domain%nx(1), i)
            end do
        end if
        
    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program uniform_channel
    !!
    !! Steady uniform overbank flow in a grid-aligned channel. 
    !! 
    !! We can compare with solution of Shiono and Knight (1991)
    !! Shiono, K. & Knight, D. W. Turbulent open-channel flows with variable depth across the channel.
    !! Journal of Fluid Mechanics, 1991, 222, 617-646 
    !!
    !! Both NS-aligned and EW-aligned channels are modelled, so we test both directions of the
    !! eddy-viscosity model.
    !!
    use global_mod, only: ip, dp, charlen, default_nonlinear_timestepping_method
    use multidomain_mod, only: multidomain_type
    use local_routines
    implicit none

    integer(ip):: i, j
    real(dp):: last_write_time
    type(multidomain_type):: md

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 60.0_dp
    real(dp), parameter :: final_time = 4000.0_dp

    ! Length/width
    real(dp) :: global_lw(2) 
    ! Lower-left corner coordinate
    real(dp) :: global_ll(2)
    ! grid size (number of x/y cells)
    integer(ip) :: global_nx(2)

    real(dp) :: global_dt
    character(len=charlen) :: channel_alignment

    call get_command_argument(1, channel_alignment)

    if(channel_alignment == 'NS_aligned') then
        channel_is_NS_aligned = .TRUE.

        global_lw = [1000.0_dp, 8000.0_dp] 
        global_ll = -global_lw/2.0_dp
        global_nx = [50*3, 200] ! Deliberate non-square dx spacing

    else if(channel_alignment == 'EW_aligned') then
        channel_is_NS_aligned = .FALSE.

        global_lw = [8000.0_dp, 1000.0_dp] 
        global_ll = -global_lw/2.0_dp
        global_nx = [200, 50*3] ! Deliberate non-square dx spacing
    else
        stop 'unknown value of channel_alignment'
    end if

    ! One domain for this problem
    allocate(md%domains(1))

    md%domains(1)%lw = global_lw
    md%domains(1)%nx = global_nx
    md%domains(1)%lower_left = global_ll
    md%domains(1)%timestepping_method = default_nonlinear_timestepping_method
    ! Set a higher eddy viscosity so we can more easily see the impact
    md%domains(1)%use_eddy_viscosity = .true.
    md%domains(1)%eddy_visc_constants(2) = lambda

    call md%setup

    ! Initial conditions
    do i = 1, size(md%domains)
        call set_initial_conditions_overbank_flow(md%domains(i))
        md%domains(i)%boundary_subroutine => boundary_forcing
    end do
    call md%make_initial_conditions_consistent()

    global_dt = 0.3_dp * md%stationary_timestep_max()

    ! Evolve the code
    do while (.true.)

        call md%write_outputs_and_print_statistics(&
            approximate_writeout_frequency = approximate_writeout_frequency)

        if (md%domains(1)%time > final_time) exit 

        call md%evolve_one_step(global_dt)

    end do

    call md%finalise_and_print_timers

END PROGRAM
