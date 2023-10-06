
module local_routines 
    !! Setup isolated building test.
    !! S. Soares-Frazao and Y. Zech, "Experimental study of dam-break flow against an isolated obstacle" Journal of Hydraulic research,
    !! 2007, VOL 45 Extra Issue, 27-36

    use global_mod, only: dp, ip, charlen, gravity
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use burn_into_grid_mod, only: xyz_lines_type
    use points_in_poly_mod, only: point_in_poly

    implicit none

    !
    ! Parameters defining the domain geometry.
    !

    ! Manning friction. This affects the placement of shocks, and so has
    ! significant on some gauge time-series (e.g. G2).
    real(dp), parameter :: manning_squared  = 0.01_dp**2 !0.015_dp**2

    ! Polygons defining the geometry
    character(len=charlen) :: domain_polygons(3) = [character(len=charlen):: &
        "poly/dam_1.csv", "poly/dam_2.csv", "poly/building.csv"]
   
    contains 

    ! Main setup routine
    subroutine set_initial_conditions(domain)
        type(domain_type), intent(inout):: domain

        integer(ip):: i, j, k, n
        character(len=charlen):: input_elevation, input_stage
        real(dp), allocatable:: x(:), y(:)
        real(dp) :: dd, gauges_1_to_4_x_coord
        real(dp) :: gauge_xy(3,8), leftmost_x(3)
        type(xyz_lines_type) :: polygons
        logical :: is_inside
        real(dp), parameter :: wall = 1.0_dp, initial_reservoir_stage = 0.4_dp, initial_downstream_stage = 0.02_dp !0.002_dp
        

        ! These define the dam-wall and the building
        call polygons%read_from_csv(domain_polygons, skip_header=1_ip)

        domain%U = 0.0_dp
        domain%manning_squared = manning_squared

        !
        ! Setup detailed elevation / manning 
        !
        allocate(x(domain%nx(1)), y(domain%nx(1)))
        x = domain%x
        do j = 1, domain%nx(2)
            y = domain%y(j)
            do i = 1, domain%nx(1)
                ! Side-slopes in Reservoir
                if(y(i) < 0.34_dp) domain%U(i,j, ELV) = 0.155_dp / 0.34_dp * (0.34_dp - y(i))
                if(y(i) > (3.6_dp - 0.34_dp)) domain%U(i,j, ELV) = 0.155_dp / 0.34_dp * (y(i) - (3.6_dp - 0.34_dp))

                ! If we are inside the polygons defining the dam/building, set the elevation to 'high'
                do k = 1, size(polygons%lines)
                    n = size(polygons%lines(k)%xyz,2)
                    call point_in_poly(n ,& 
                        polygons%lines(k)%xyz(1,:), polygons%lines(k)%xyz(2,:), &
                        x(i), y(i), is_inside)
                    if(is_inside) domain%U(i,j,ELV) = wall
                end do

                ! Initial stage
                if(x(i) < 6.75_dp) then
                    domain%U(i,j, STG) = initial_reservoir_stage
                else
                    domain%U(i,j, STG) = initial_downstream_stage
                end if

            end do
        end do

        deallocate(x,y)

        print*, 'Elevation range: ', minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))

        ! Wall boundaries (without boundary conditions)
        domain%U(:,1,ELV) = wall
        domain%U(:,domain%nx(2),ELV) = wall
        domain%U(domain%nx(1),:,ELV) = wall
        domain%U(1,:,ELV) = wall


        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV) + 1.0e-07_dp)

    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program isolated_building
    !! Run the isolated building test.
    !! S. Soares-Frazao and Y. Zech, "Experimental study of dam-break flow against an isolated obstacle" Journal of Hydraulic research,
    !! 2007, VOL 45 Extra Issue, 27-36

    use global_mod, only: ip, dp, minimum_allowed_depth, default_nonlinear_timestepping_method
    use domain_mod, only: domain_type
    use multidomain_mod, only: multidomain_type
    use timer_mod
    use logging_mod, only: log_output_unit
    use local_routines
    implicit none

    ! Useful misc variables
    integer(ip):: j, i, nd

    ! Type holding all domains 
    type(multidomain_type) :: md

    type(timer_type) :: program_timer

    real(dp), parameter :: mesh_refine = 1.0_dp ! 8.0_dp ! Increase resolution by this amount
   
    ! This can voilate the theoretical rk2 time-stepping limit, but still the solution is good. 
    ! That might be common for rk2, for instance see:
    ! Andrew Giuliani, Lilia Krivodonova, On the optimal CFL number of SSP methods for hyperbolic problems.
    real(dp) ::  global_dt = 0.02_dp / mesh_refine
    !real(dp) :: local_dt

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 0.01_dp
    real(dp), parameter :: final_time = 40._dp

    !
    ! Key geometric parameters (others defined in the 'local routines' module)
    !
    ! Length/Width
    real(dp), parameter :: global_lw(2) = [35.8_dp, 3.60_dp]
    ! Lower-left coordinate
    real(dp), parameter :: global_ll(2) = [0.0_dp, 0.0_dp] ! Matching sketch in supplementary material
    !real(dp), parameter :: global_ll(2) = [-0.15_dp, 0.0_dp] ! Matching Fig 1 in paper


    ! Grid size (number of x/y cells) in outer domain
    integer(ip), parameter:: global_nx(2) = nint(global_lw*10*mesh_refine, ip)
    integer(ip), parameter :: write_grids_and_print_every_nth_step = 10_ip


    call program_timer%timer_start('setup')

    ! nd domains in this model
    nd = 1
    allocate(md%domains(nd))

    !
    ! Setup basic metadata
    !

    ! Main domain
    md%domains(1)%lower_left =global_ll
    md%domains(1)%lw = global_lw
    md%domains(1)%nx = global_nx
    md%domains(1)%dx = md%domains(1)%lw/md%domains(1)%nx
    md%domains(1)%timestepping_refinement_factor = 1_ip
    md%domains(1)%dx_refinement_factor = 1.0_dp
    md%domains(1)%timestepping_method = default_nonlinear_timestepping_method !'cliffs' !'rk2'
    md%domains(1)%cliffs_minimum_allowed_depth = 0.01_dp
    !md%domains(1)%theta = 4.0_dp
    !md%domains(1)%use_eddy_viscosity = .true.
    !md%domains(1)%eddy_visc_coef(1:2) = [0.0_dp, 0.5_dp]

    print*, 1, ' lw: ', md%domains(1)%lw, ' ll: ', md%domains(1)%lower_left, ' dx: ', md%domains(1)%dx, &
        ' nx: ', md%domains(1)%nx

    if(md%domains(1)%timestepping_method == 'rk2n') global_dt = global_dt * 4.0_dp
     
    ! Allocate domains and prepare comms
    call md%setup()

    ! Initial conditions
    do j = 1, size(md%domains)
        call set_initial_conditions(md%domains(j))
    end do

    call md%make_initial_conditions_consistent()
    
    ! NOTE: For stability in 'null' regions, we set them to 'high land' that
    ! should be inactive. 
    call md%set_null_regions_to_dry()

    ! gauges
    call md%set_point_gauges_from_csv("point_gauges.csv", skip_header=1_ip)

    ! Print the gravity-wave CFL limit, to guide timestepping
    !local_dt = HUGE(1.0_dp)
    do j = 1, size(md%domains)
        !local_dt = min(local_dt, md%domains(j)%stationary_timestep_max() * 0.95_dp)
        print*, 'domain: ', j, 'ts: ', &
            md%domains(j)%stationary_timestep_max()
    end do

    print*, 'End setup'
    call program_timer%timer_end('setup')
    call program_timer%timer_start('evolve')

    ! Evolve the code
    do while (.true.)

        ! Write gauges every time-step, but print and write grids less often
        call program_timer%timer_start('IO')
        call md%write_outputs_and_print_statistics(approximate_writeout_frequency=approximate_writeout_frequency,&
            write_grids_less_often=write_grids_and_print_every_nth_step, &
            print_less_often = write_grids_and_print_every_nth_step,&
            write_gauges_less_often= 1_ip)
        call program_timer%timer_end('IO')

        if (md%domains(1)%time > final_time) exit

        call md%evolve_one_step(global_dt)
        !call md%evolve_one_step(local_dt)
        
    end do

    call program_timer%timer_end('evolve')
    call md%finalise_and_print_timers

    print*, ''
    call program_timer%print(output_file_unit=log_output_unit)

end program
