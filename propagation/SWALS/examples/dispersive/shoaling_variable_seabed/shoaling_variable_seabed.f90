module local_routines 
    !!
    !! Coulaud et al (2025) http://dx.doi.org/10.1016/j.coastaleng.2024.104645
    !! Case 2: Unsteady shoaling of a wave train
    !!

    use global_mod, only: dp, ip, charlen, pi, gravity
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    implicit none

    ! Variables defined in paper
    integer(ip), parameter :: N = 10
    real(dp), parameter :: d1 = 1.0_dp, &
        a1 = 0.125_dp * d1, &
        dmin = 0.2_dp * d1, &
        L = 50.0_dp * d1, &
        xl = 0.0_dp, &
        xr = 74.0_dp * d1

    contains

    subroutine set_initial_conditions(domain, long_dimension)
        type(domain_type), intent(inout):: domain
        character(len=charlen), intent(in) :: long_dimension

        real(dp):: x, y, elev
        integer(ip):: j, i, k
        real(dp):: gauge_xy(2,11), wall

        domain%msl_linear = 0.0_dp
        domain%U(:,:,UH:VH) = 0.0_dp

        if(long_dimension == 'y') then
            ! Initial elevation
            do i = 1, domain%nx(1)
                where(domain%y <= L) domain%U(i,:,ELV) = -(dmin/d1 + (1.0_dp - dmin/d1)/cosh(tan((pi*domain%y/(2*L)))))
                where(domain%y > L) domain%U(i,:,ELV) = -(dmin/d1)
            end do

            ! Initial stage
            do i = 1, domain%nx(1)
                ! Plot from the paper makes clear this is only in the left part of the domain
                where(domain%y <= L) &
                    domain%U(i,:,STG) = d1 * (a1/d1)*(cos(2.0_dp*pi*N*domain%y/L)/cosh(tan(pi*domain%y/(2*L))))
            end do

        else if(long_dimension == 'x') then
            ! Initial elevation
            do j = 1, domain%nx(2)
                where(domain%x <= L) domain%U(:,j,ELV) = -(dmin/d1 + (1.0_dp - dmin/d1)/cosh(tan((pi*domain%x/(2*L)))))
                where(domain%x > L) domain%U(:,j,ELV) = -(dmin/d1)
            end do

            ! Initial stage
            do j = 1, domain%nx(2)
                ! Plot from the paper makes clear this is only in the left part of the domain
                where(domain%x <= L) &
                    domain%U(:,j,STG) = d1 * (a1/d1)*(cos(2.0_dp*pi*N*domain%x/L)/cosh(tan(pi*domain%x/(2*L))))
            end do
        else
            stop "unknown long dimension"
        end if

        ! Reflective boundaries on 4 sides
        wall = 1.0_dp
        domain%U(:, 1:2, ELV) = wall
        domain%U(1:2, :, ELV) = wall
        domain%U(:, domain%nx(2)-1:domain%nx(2), ELV) = wall
        domain%U(domain%nx(1)-1:domain%nx(1), :, ELV) = wall
   
        ! Stage >= bed 
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV))
        
    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program shoaling_variable_seabed
    !!
    !! Coulaud et al (2025) http://dx.doi.org/10.1016/j.coastaleng.2024.104645
    !! Case 2: Unsteady shoaling of a wave train
    !!
    use global_mod, only: ip, dp
    use multidomain_mod, only: multidomain_type
    use boundary_mod, only: boundary_stage_transmissive_normal_momentum
    use local_routines
    implicit none

    type(multidomain_type):: md

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 0.1_dp
    real(dp), parameter :: final_time = 45.0_dp*sqrt(d1/gravity)

    ! Domain info
    character(charlen) :: timestepping_method, buf, long_dimension

    ! Local variables 
    real(dp) :: dx_in, dx, timestep

    ! Timestepping method
    call get_command_argument(1, timestepping_method)  
    ! Resolution
    call get_command_argument(2, buf)
    read(buf, *) dx_in

    ! Set the long-dimension (either 'y' or 'x'). Use this run in 1D along x or y axes
    call get_command_argument(3, long_dimension) 
      
    
    ! Resolution
    dx = d1*dx_in 

    ! Setup model with 1 domain
    allocate(md%domains(1))

    if(long_dimension == 'y') then
        ! Uneven grid size == (dx along y), (3*dx along x)
        md%domains(1)%lw = [5*dx*3, xr + 2*dx]
        md%domains(1)%lower_left =  [-2.5_dp * dx*3, -2.5*dx] ! Real flow over x==0, y >= -(dx/2) (with space for reflective boundaries)
        md%domains(1)%nx = nint(md%domains(1)%lw/[3*dx, dx])
    else if(long_dimension == 'x') then
        ! Uneven grid size == (dx along x), (3*dx along y)
        md%domains(1)%lw = [xr + 2*dx, 5*dx*3]
        md%domains(1)%lower_left =  [-2.5_dp * dx, -2.5*dx*3] ! Real flow over y==0, x >= -(dx/2) (with space for reflective boundaries)
        md%domains(1)%nx = nint(md%domains(1)%lw/[dx, 3*dx])
    else
        stop "unknown long_dimension"
    end if

    md%domains(1)%timestepping_method = timestepping_method
    md%domains(1)%use_dispersion = .true. !
    !md%domains(1)%nc_grid_output%flush_every_n_output_steps = 1_ip !

    md%domains(1)%nontemporal_grids_to_store = [character(len=charlen):: 'max_stage', 'min_stage', 'elevation0']

    ! Tapering off of dispersive terms (won't have any effect)
    md%domains(1)%ds%td1 = 0.01_dp
    md%domains(1)%ds%td2 = 0.005_dp
    !md%domains(1)%ds%tridiagonal_inner_iter = 4_ip

    ! Non-TVD limiting for finite volume schemes
    md%domains(1)%theta = 4.0_dp

    ! Output directory should record the solver and resolution
    md%output_basedir = 'OUTPUTS/' // trim(timestepping_method) // &
        '_long_dimension_' // trim(long_dimension) // '_' // trim(buf)

    call md%setup

    call set_initial_conditions(md%domains(1), long_dimension)

    call md%make_initial_conditions_consistent() ! Get the initial volume right

    ! Fixed timestep  
    timestep = md%stationary_timestep_max() * 0.8_dp 
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
