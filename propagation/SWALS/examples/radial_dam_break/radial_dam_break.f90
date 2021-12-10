module local_routines 
    !!
    !! Setup the radial dam-break problem
    !!
    use global_mod, only: dp, ip, wall_elevation
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    implicit none

    contains 

    subroutine set_initial_conditions_radial_dam(domain)            
        !
        ! This function sets the initial conditions in a domain
        !
        class(domain_type), target, intent(inout):: domain
        integer(ip):: i,j
        real(dp):: x, y, cx, cy, initial_stage_1, initial_stage_2, radius

        initial_stage_1 = 1.0_dp 
        initial_stage_2 = 10.0_dp
        radius = 50.0_dp

        ! Stage
        domain%U(:,:,STG) = initial_stage_1 
        ! Depth integrated velocity
        domain%U(:,:, UH:VH) = 0.0_dp

        ! Elevation
        domain%U(:,:,ELV) = 0._dp
        ! Wall boundaries (without boundary conditions)
        domain%U(1,:,ELV) = 20.0_dp !wall_elevation 
        domain%U(domain%nx(1),:,4) = 20.0_dp !wall_elevation 
        domain%U(:,1,ELV) = 20.0_dp !wall_elevation
        domain%U(:,domain%nx(2),ELV) = 20.0_dp !wall_elevation

        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV))

        ! Add a stage perturbation
        cx = (domain%lw(1))*0.5
        cy = (domain%lw(2))*0.5
        do i = 1,domain%nx(1)
            do j = 1, domain%nx(2)
                ! Set perturbation based on distance to the centre
                x = (i-0.5_dp)*domain%dx(1) - cx
                y = (j-0.5_dp)*domain%dx(2) - cy 
                if( x*x + y*y < (radius)**2) then
                    domain%U(i,j,STG) = initial_stage_2 !10.0_dp
                end if
            end do
        end do

        ! Manning friction squared
        domain%manning_squared = 0.0_dp

    end subroutine


end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program radial_dam_break
    !!
    !! Radial dam-break problem
    !!
    use global_mod, only: ip, dp, charlen, default_nonlinear_timestepping_method
    use multidomain_mod, only: multidomain_type
    use domain_mod, only: UH, VH
    use local_routines
    implicit none

    integer(ip):: i, nsteps
    real(dp):: last_write_time
    type(multidomain_type) :: md

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 0.2_dp
    real(dp), parameter :: final_time = 2.0_dp * 1

    ! length/width
    real(dp), parameter, dimension(2):: global_lw = [400._dp, 400._dp] 
    ! lower-left corner coordinate
    real(dp), parameter, dimension(2):: global_ll = [0._dp, 0._dp]
    ! grid size (number of x/y cells)
    integer(ip), parameter, dimension(2):: global_nx = [400, 400] * 2 + 1

    ! Misc
    integer :: j, nd
    real(dp) :: global_dt

    nd = 1 ! Number of domains in model
    allocate(md%domains(nd))

    !
    ! Set the domain properties
    !
    md%domains(1)%timestepping_method = default_nonlinear_timestepping_method

    ! Domain Geometry
    md%domains(1)%lw = global_lw
    md%domains(1)%lower_left = global_ll
    md%domains(1)%nx = global_nx

    ! Output variables to store
    md%domains(1)%time_grids_to_store = [character(len=charlen):: 'stage', 'uh', 'vh']
    md%domains(1)%nontemporal_grids_to_store = [character(len=charlen):: 'max_stage', 'max_speed', 'max_flux', &
        'arrival_time', 'elevation0', 'manning_squared']

    ! Define how the "arrival time" statistic is calculated. It is the time when stage 
    ! exceeds the sum of these two variables.
    md%domains(1)%msl_linear = 0.0_dp
    md%domains(1)%arrival_stage_threshold_above_msl_linear = 1.01_dp

    ! Setup the multidomain -- note for some models this can change the number of domains (e.g. for parallel),
    ! although that is not done in this example
    call md%setup()

    do j = 1, size(md%domains)
        ! Set initial conditions on each domain
        call set_initial_conditions_radial_dam(md%domains(j))
    end do

    ! Time-step at which we evolve the solution
    global_dt = 0.75_dp * md%stationary_timestep_max()

    ! Evolve the solution
    do while (.TRUE.)

        call md%write_outputs_and_print_statistics(approximate_writeout_frequency=approximate_writeout_frequency)

        if (md%domains(1)%time > final_time) then
            exit 
        end if

        call md%evolve_one_step(global_dt)

    end do

    call md%finalise_and_print_timers

end program
