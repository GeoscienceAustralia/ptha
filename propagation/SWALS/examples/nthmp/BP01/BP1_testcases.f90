module local_routines 
    use global_mod, only: dp, ip, charlen, wall_elevation, g => gravity
    use domain_mod, only: domain_type
    use read_raster_mod, only: read_gdal_raster
    use which_mod, only: which
    use file_io_mod, only: count_file_lines
    use linear_interpolator_mod, only: linear_interpolator_type
    implicit none

    contains 

    subroutine set_initial_conditions_BP1(domain, d, beach_slope, land_length, sea_length, X0, X1, H, gamma0)
        class(domain_type), target, intent(inout):: domain
        real(dp), intent(in) :: d, beach_slope, land_length, sea_length, X0, X1, H, gamma0 

        integer(ip) :: i,j, ngauge
        real(dp):: x,y, elev, stg, vel
        real(dp), allocatable:: gauge_xy(:,:)

        do j = 1, domain%nx(2)
            do i = 1, domain%nx(1)
                x = domain%x(i)
                y = domain%y(j)
                elev = max(-x * beach_slope, -d)
                stg = H * (1.0_dp/cosh(gamma0*(x - X1)/d))**2
                vel = - sqrt(g/d) * stg
                domain%U(i,j, 4) = elev
                domain%U(i,j,1) = stg
                domain%U(i,j,2) = vel*d
                domain%U(i,j,3) = 0.0_dp
            end do
        end do
        print*, 'Stage initial range: ', maxval(domain%U(:,:,1)), minval(domain%U(:,:,1))

        ! Add reflective walls, with zero velocity, and stage>=elev
        domain%U(:,1,4) = wall_elevation
        domain%U(:,domain%nx(2),4) = wall_elevation
        domain%U(domain%nx(1),:,4) = wall_elevation
        domain%U(:,:,1) = max(domain%U(:,:,1), domain%U(:,:,4)+1.0e-07_dp)
        domain%U(:,1,2) = 0.0_dp
        domain%U(domain%nx(1),:,2) = 0.0_dp
        domain%U(:,domain%nx(2),2) = 0.0_dp

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program bp1
    use global_mod, only: ip, dp, minimum_allowed_depth, pi
    use domain_mod, only: domain_type
    use linear_interpolator_mod, only: linear_interpolator_type
    use local_routines
    implicit none

    integer(ip):: j
    real(dp):: last_write_time, rain_rate
    type(domain_type):: domain

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 0.025_dp
    real(dp), parameter :: final_time = 50.0_dp

    ! Domain info
    character(charlen) :: timestepping_method, tempchar != 'linear' !'rk2' !'linear'
    
    !! length/width
    real(dp), dimension(2) :: global_lw, global_ll
    integer(ip), dimension(2) :: global_nx 

    ! Local variables 
    real(dp) :: timestep, L, dx, tank_width, tank_length, initial_depth, &
        beach_slope, beta, X0, X1, gamma0, H, land_length, sea_length, h_on_d

    call get_command_argument(1, timestepping_method)  
    print*, 'timestepping method: ', TRIM(timestepping_method)
    call get_command_argument(2, tempchar)
    read(tempchar, *) initial_depth
    print*, 'initial depth: ', initial_depth
    call get_command_argument(3, tempchar)
    read(tempchar, *) h_on_d
    print*, 'h_on_d: ', h_on_d 

    ! Geometric parameters as specified in the description
    !initial_depth = 1.0_dp ! This is arbitrary
    beta = atan(1.0_dp/19.85_dp)
    beach_slope = tan(beta)
    X0 = initial_depth / beach_slope 
    !H = 0.019_dp * initial_depth
    H = h_on_d * initial_depth
    gamma0 = sqrt(3.0_dp * H / (4.0_dp * initial_depth))
    L = initial_depth * acosh(sqrt(20.0_dp))/gamma0
    X1 = X0 + L

    print*, 'd: ', initial_depth, ' H: ', H, ' beach slope: ', beach_slope, &
        ' X0 :', X0, ' X1: ', X1, ' L: ', L, ' gamma: ', gamma0
    !stop
 
    ! Tank geometry 
    tank_width = 2.0_dp
    land_length = 3.0_dp * X0 ! Set this so that inundation can occur
    sea_length = X0 + L + 10.0_dp * L ! Set this so the sea boundary is far enough away.
    tank_length = sea_length + land_length
 
    ! Resolution
    dx = 0.1_dp !(0.5_dp/12.0_dp)

    ! Model parameters. Note we apply the boundary on the right. x = 0 is the shoreline
    global_lw = [tank_length, tank_width]
    global_ll = [-land_length, -tank_width/2.0_dp]
    global_nx = nint(global_lw/dx)

    print*, 'll: ', global_ll

    domain%timestepping_method = timestepping_method
    !if(timestepping_method == 'euler') domain%theta = 0.0_dp
    !if(timestepping_method == 'rk2') domain%theta = 1.0_dp

    ! Allocate domain -- must have set timestepping method BEFORE this
    call domain%allocate_quantities(global_lw, global_nx, global_ll)

    ! Call local routine to set initial conditions
    call set_initial_conditions_BP1(domain, initial_depth, beach_slope, &
        land_length, sea_length, X0, X1, H, gamma0)

    ! Linear requires a fixed timestep 
    if (timestepping_method == 'linear') then
        timestep = domain%linear_timestep_max() 
    end if


    ! Trick to get the code to write out just after the first timestep
    last_write_time = domain%time - approximate_writeout_frequency

    ! Evolve the code
    do while (.true.)

        if(domain%time - last_write_time >= approximate_writeout_frequency) then

            last_write_time = last_write_time + approximate_writeout_frequency

            ! This avoids any artefacts in the numerical update of the model
            ! which should be overwritten by the boundary condition
            !call domain%update_boundary()

            call domain%print()
            call domain%write_to_output_files(time_only=.true.)
            call domain%write_gauge_time_series()
            print*, 'Mass balance: ', domain%mass_balance_interior()

        end if

        if (domain%time > final_time) exit

        ! Variable timestep
        if(timestepping_method == 'linear') then
            call domain%evolve_one_step(timestep = timestep)
        else
            call domain%evolve_one_step()
        end if

    end do

    call domain%write_max_quantities()

    ! Print timing info
    call domain%timer%print()

    call domain%finalise()

end program
