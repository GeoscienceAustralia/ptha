module local_routines 
    !!
    !! Setup NTHMP benchmark problem 4.
    !!
    use global_mod, only: dp, ip, charlen, wall_elevation, g => gravity
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    implicit none

    contains 

    subroutine set_initial_conditions_BP1(domain, d, beach_slope, land_length, sea_length, X0, X1, H, gamma0)
        class(domain_type), target, intent(inout):: domain
        real(dp), intent(in) :: d, beach_slope, land_length, sea_length, X0, X1, H, gamma0 

        integer(ip) :: i,j, ngauge
        real(dp):: x,y, elev, stage, vel
        real(dp), allocatable:: gauge_xy(:,:)

        do j = 1, domain%nx(2)
            do i = 1, domain%nx(1)
                x = domain%x(i)
                y = domain%y(j)
                elev = max(-x * beach_slope, -d)
                stage = H * (1.0_dp/cosh(gamma0*(x - X1)/d))**2
                vel = - sqrt(g/d) * stage
                domain%U(i,j,ELV) = elev
                domain%U(i,j,STG) = stage
                domain%U(i,j,UH) = vel*d
                domain%U(i,j,VH) = 0.0_dp
            end do
        end do
        !print*, 'Stage initial range: ', maxval(domain%U(:,:,1)), minval(domain%U(:,:,1))

        ! Add reflective walls, with zero velocity, and stage>=elev
        domain%U(:,1,ELV) = wall_elevation
        domain%U(:,domain%nx(2),ELV) = wall_elevation
        domain%U(domain%nx(1),:,ELV) = wall_elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV)+1.0e-07_dp)
        domain%U(:,1,UH) = 0.0_dp
        domain%U(domain%nx(1),:,UH) = 0.0_dp
        domain%U(:,domain%nx(2),UH) = 0.0_dp

        if(allocated(domain%manning_squared)) then
            ! The runup in this problem can be sensitive to friction
            ! This value is 'tuned' to work ok on both problems, but we could
            ! change runup for the large-amplitude case quite a bit by adjusting this.
            domain%manning_squared = 0.005_dp*0.005_dp
        end if

        ! Get gauge points -- every d/10
        ngauge = nint((maxval(domain%x) - minval(domain%x))/d * 10 )
        allocate(gauge_xy(2, ngauge))
        do i = 1, ngauge
            gauge_xy(1,i) = minval(domain%x) + d/20.0_dp + (i-1)*d/10.0_dp
            gauge_xy(2,i) = 0.0_dp
        end do

        call domain%setup_point_gauges(gauge_xy)

    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program BP04
    !!
    !! NTHMP benchmark problem 4 -- solitary wave on a simple beach. This has
    !! various experimental results, leading to a dimensionless relation between
    !! the onshore runup and offshore wave height.
    !!

    use global_mod, only: ip, dp
    use multidomain_mod, only: multidomain_type
    use local_routines
    implicit none

    type(multidomain_type) :: md

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 0.025_dp
    real(dp), parameter :: final_time = 50.0_dp

    ! Domain info
    character(charlen) :: timestepping_method, tempchar
    
    ! length/width
    real(dp) :: global_lw(2), global_ll(2)
    integer(ip) :: global_nx(2) 

    ! Local variables 
    real(dp) :: timestep, L, dx, tank_width, tank_length, initial_depth, &
        beach_slope, beta, X0, X1, gamma0, H, land_length, sea_length, h_on_d

    call get_command_argument(1, timestepping_method)  
    call get_command_argument(2, tempchar)
    read(tempchar, *) initial_depth
    call get_command_argument(3, tempchar)
    read(tempchar, *) h_on_d

    ! Geometric parameters as specified in the description
    beta = atan(1.0_dp/19.85_dp)
    beach_slope = tan(beta)
    X0 = initial_depth / beach_slope 
    H = h_on_d * initial_depth
    gamma0 = sqrt(3.0_dp * H / (4.0_dp * initial_depth))
    L = initial_depth * acosh(sqrt(20.0_dp))/gamma0
    X1 = X0 + L

    ! Tank geometry 
    tank_width = 0.14_dp
    land_length = 3.0_dp * X0 ! Set this so that inundation can occur
    sea_length = X0 + L + 10.0_dp * L ! Set this so the sea boundary is far enough away.
    tank_length = sea_length + land_length
 
    ! Resolution
    dx = 0.02_dp !(0.5_dp/12.0_dp)

    ! Model parameters. Note we apply the boundary on the right. x = 0 is the shoreline
    global_lw = [tank_length, tank_width]
    global_ll = [-land_length, -tank_width/2.0_dp]
    global_nx = nint(global_lw/dx)

    ! Setup md with 1 domain
    allocate(md%domains(1))
    md%domains(1)%lw = global_lw
    md%domains(1)%lower_left = global_ll
    md%domains(1)%nx = global_nx
    md%domains(1)%timestepping_method = timestepping_method

    call md%setup
    
    ! Call local routine to set initial conditions
    call set_initial_conditions_BP1(md%domains(1), initial_depth, beach_slope, &
        land_length, sea_length, X0, X1, H, gamma0)

    ! Fixed timestep
    timestep = md%stationary_timestep_max() * 0.5_dp

    ! Evolve the code
    do while (.true.)

        call md%write_outputs_and_print_statistics(&
            approximate_writeout_frequency=approximate_writeout_frequency, &
            write_grids_less_often = 99999_ip)

        if (md%domains(1)%time > final_time) exit

        call md%evolve_one_step(timestep)

    end do

    call md%finalise_and_print_timers

end program
