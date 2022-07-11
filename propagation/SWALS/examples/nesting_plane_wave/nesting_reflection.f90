module local_routines 
    !! Setup for plane wave propagation problem with nesting
    use global_mod, only: dp, ip, charlen, wall_elevation, gravity, pi
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use linear_interpolator_mod, only: linear_interpolator_type

    implicit none

    contains 

    ! Main setup routine
    subroutine set_initial_conditions(domain)
        class(domain_type), target, intent(inout):: domain

        integer(ip):: i, j
        character(len=charlen):: initial_stage
        real(dp), allocatable:: x(:)
        real(dp) :: wall, stage0, x0, d0, a0, k0
        real(dp) :: gauge_xy(3,3)

        wall = 50.0_dp
        allocate(x(domain%nx(1)))
        x = domain%x
      
        ! Define the initial wave 
        x0 = 0.0_dp
        d0 = -100.0_dp
        a0 = 0.001_dp
        k0 = 1.0_dp / (250.0_dp * 38.0_dp)

        do j = 1, domain%nx(2)
            do i = 1, domain%nx(1)
               if(abs(x(i) - x0)*k0 < 3.25_dp) then
                   stage0 = a0 * cos(2*pi*(x(i) - x0)*k0)
               else
                   stage0 = 0.0_dp
               end if
               domain%U(i,j,STG) = stage0
               domain%U(i,j,ELV) = d0 
               domain%U(i,j,UH) = (-d0) * sqrt(gravity/(-d0))*stage0
               domain%U(i,j,VH) = 0.0_dp
            end do    
        end do
        
        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV) + 1.0e-07_dp)

        deallocate(x)

        print*, 'Elevation range: ', minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))

        if(domain%timestepping_method /= 'linear') then
            domain%manning_squared = 0.0_dp
        end if

    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program nesting_reflection 
    !! Plane wave propagation problem with nesting. This problem can highlight weaknesses
    !! in structured grid finite-volume methods with TVD limiters, which tend to dissipate the wave too 
    !! quickly. That is unrelated to the use of nesting - rather it's due to limiter clipping.
    !! That issue can be solved using non-TVD limiting, or no limiting. Alternatively, the leapfrog schemes
    !! are better suited to this kind of problem.

    use global_mod, only: ip, dp, minimum_allowed_depth, charlen
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
    
    real(dp) ::  global_dt = 3.985_dp / mesh_refine

    ! Approx timestep between outputs
    real(dp) :: approximate_writeout_frequency = 5.00_dp
    real(dp) :: final_time = 3600.0_dp

    integer(ip), parameter :: inner_mesh_refine = 3_ip

    ! Length/width
    real(dp), parameter :: global_lw(2) = [100000.0_dp , 5000.0_dp/mesh_refine]
    ! Lower-left corner coordinate
    real(dp), parameter :: global_ll(2) = -global_lw/2.0_dp
    ! Resolution
    real(dp), parameter :: res_d1 = 250.0_dp/mesh_refine 
    real(dp), parameter :: res_d2 = 250.0_dp/mesh_refine * 1.0_dp/inner_mesh_refine

    real(dp) :: lower_left(2), upper_right(2)

    ! Standard test version
    ! Run the '2 domain' version
    integer(ip), parameter :: n_domains = 2_ip
    ! Allow overshoots to prevent excess dissipation
    real(dp), parameter :: nonlinear_theta = 4.0_dp
    character(len=charlen) :: compute_fluxes_inner_method = "DE1_low_fr_diffusion" ! "DE1" -- distinguishes 'rk2' and 'midpoint'

    !@ Alternate version -- shows the dissipation of the regular 'DE1' type approach.
    !@ Similar behaviour is seen in 1D "Basilisk", which employs a very similar algorithm.
    !@ Note the dissipation also occurs with >1 domain (but n_domains=1 for simplicity here)
    !integer(ip), parameter :: n_domains = 1_ip
    !real(dp), parameter :: nonlinear_theta = 1.3_dp
    !character(len=charlen), parameter :: compute_fluxes_inner_method = 'DE1' !'DE1_low_fr_diffusion'

    character(len=charlen) :: ts_method

    call program_timer%timer_start('setup')

    call get_command_argument(1, ts_method)

    ! Periodic BC
    md%periodic_xs = [global_ll(1) , global_ll(1) + global_lw(1)]
    md%periodic_ys = [global_ll(2) , global_ll(2) + global_lw(2)]

    !
    ! Setup basic metadata
    !
    if(n_domains == 4) then
        ! Version with 4 domains
        nd = 4
        allocate(md%domains(nd))

        ! Left domain
        md%domains(1)%lower_left = global_ll
        md%domains(1)%lw = [global_lw(1)/3.0_dp, global_lw(2)]
        md%domains(1)%nx = nint(md%domains(1)%lw/res_d1)
        md%domains(1)%dx = md%domains(1)%lw/md%domains(1)%nx
        md%domains(1)%timestepping_refinement_factor = 1_ip
        md%domains(1)%dx_refinement_factor = 1.0_dp
        md%domains(1)%timestepping_method = ts_method

        ! Main domain
        md%domains(2)%lower_left = [global_ll(1) + 1.0_dp/3.0_dp * global_lw(1), global_ll(2)]
        md%domains(2)%lw = [global_lw(1)/3.0_dp, global_lw(2)/2.0_dp] 
        md%domains(2)%nx = nint(md%domains(2)%lw/res_d2)
        md%domains(2)%dx = md%domains(2)%lw/md%domains(2)%nx
        md%domains(2)%timestepping_refinement_factor = 1_ip
        md%domains(2)%dx_refinement_factor = nint(res_d1/res_d2)
        md%domains(2)%timestepping_method = ts_method
        
        ! Right domain
        md%domains(3)%lower_left = [global_ll(1) + global_lw(1)*2.0/3.0, global_ll(2)]
        md%domains(3)%lw = [global_lw(1)/3.0_dp, global_lw(2)]
        md%domains(3)%nx = nint(md%domains(1)%lw/res_d1)
        md%domains(3)%dx = md%domains(1)%lw/md%domains(1)%nx
        md%domains(3)%timestepping_refinement_factor = 1_ip
        md%domains(3)%dx_refinement_factor = 1.0_dp
        md%domains(3)%timestepping_method = ts_method 

        ! Main domain
        md%domains(4)%lower_left = [global_ll(1) + 1.0_dp/3.0_dp * global_lw(1), 0.0_dp]
        md%domains(4)%lw = [global_lw(1)/3.0_dp, global_lw(2)/2.0_dp]
        md%domains(4)%nx = nint(md%domains(2)%lw/res_d2)
        md%domains(4)%dx = md%domains(2)%lw/md%domains(2)%nx
        md%domains(4)%timestepping_refinement_factor = 1_ip
        md%domains(4)%dx_refinement_factor = nint(res_d1/res_d2)
        md%domains(4)%timestepping_method = ts_method

        ! Allow overshoots to prevent excess dissipation
        md%domains(1)%theta = nonlinear_theta
        md%domains(2)%theta = nonlinear_theta
        md%domains(3)%theta = nonlinear_theta
        md%domains(4)%theta = nonlinear_theta

    else if(n_domains == 2) then
        ! Version with 2 domains. 
        nd = 2
        allocate(md%domains(nd))

        ! Large domain
        md%domains(1)%lower_left = global_ll
        md%domains(1)%lw = [global_lw(1), global_lw(2)]
        md%domains(1)%nx = nint(md%domains(1)%lw/res_d1)
        md%domains(1)%dx = md%domains(1)%lw/md%domains(1)%nx
        md%domains(1)%timestepping_refinement_factor = 1_ip
        md%domains(1)%dx_refinement_factor = 1.0_dp
        md%domains(1)%timestepping_method = ts_method

        ! Inner domain
        lower_left  = [global_ll(1) + 1.0_dp/3.0_dp * global_lw(1), global_ll(2)]
        upper_right = [global_ll(1) + 2.0_dp/3.0_dp * global_lw(1), global_ll(2) + global_lw(2)]
        call md%domains(2)%match_geometry_to_parent(&
            parent_domain = md%domains(1), &
            lower_left = lower_left, &
            upper_right = upper_right, &
            dx_refinement_factor = inner_mesh_refine, &
            timestepping_refinement_factor = inner_mesh_refine)
        md%domains(2)%timestepping_method = ts_method

        ! Allow overshoots to prevent excess dissipation
        md%domains(1)%theta = nonlinear_theta
        md%domains(2)%theta = nonlinear_theta

    else if(n_domains == 1) then
        ! Single domain
        nd = 1
        allocate(md%domains(nd))

        md%domains(1)%lower_left = global_ll
        md%domains(1)%lw = [global_lw(1), global_lw(2)]
        md%domains(1)%nx = nint(md%domains(1)%lw/res_d1)
        md%domains(1)%dx = md%domains(1)%lw/md%domains(1)%nx
        md%domains(1)%timestepping_refinement_factor = 1_ip
        md%domains(1)%dx_refinement_factor = 1.0_dp
        md%domains(1)%timestepping_method = ts_method 

        ! Allow overshoots to prevent excess dissipation
        md%domains(1)%theta = nonlinear_theta

    end if

    ! Experiment with non-default compute fluxes approaches
    do j = 1, size(md%domains)
        md%domains(j)%compute_fluxes_inner_method = compute_fluxes_inner_method
    end do

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

    print*, 'End setup'

    ! Print the gravity-wave CFL limit, to guide timestepping
    do j = 1, size(md%domains)
        print*, 'domain: ', j, 'ts: ', &
            md%domains(j)%stationary_timestep_max()
    end do

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
