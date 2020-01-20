module local_routines 
    !!
    !! Simulate wave motion in a paraboloid bowl, which has an analytical solution 
    !!
    use global_mod, only: dp, ip, charlen, wall_elevation, gravity
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use read_raster_mod, only: multi_raster_type
    use logging_mod, only: log_output_unit
    implicit none

    contains 

    subroutine set_initial_conditions_paraboloid_bowl(domain)            
        class(domain_type), target, intent(inout):: domain
        integer(ip):: i, j
        real(dp) :: x, y, L, D, r0, r, A

        r0 = 2000.0_dp
        L = 2500.0_dp
        D = 1000.0_dp
        A = (L**4 - r0**4)/(r0**4 + L**4)

        ! Set stage and elevation 
        ! This saves memory compared to doing it all at once.
        do j = 1, domain%nx(2)
            y = domain%y(j)
            do i = 1, domain%nx(1)
                x = domain%x(i)

                r = sqrt(x*x + y*y)
                ! Bed elevation
                domain%U(i,j,ELV) = D*(r*r/(L*L) - 1.0_dp)
                ! Stage
                domain%U(i,j,STG) = D*( & 
                    (sqrt(1.0_dp-A*A))/(1.0_dp - A) & 
                    -1.0_dp - r*r/(L*L)*((1.0_dp - A*A)/((1.0_dp - A)**2) - 1.0_dp) &
                    )
            end do
        end do

        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV) + 1.0e-07_dp)

        write(log_output_unit,*) 'Stage range is: ', minval(domain%U(:,:,STG)), maxval(domain%U(:,:,STG))
        write(log_output_unit,*) 'Elev range is: ', minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))

        if(domain%timestepping_method == 'cliffs') then
            domain%cliffs_minimum_allowed_depth = 5_dp
        end if

    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program run_paraboloid_basin
    !!
    !! Simulate wave motion in a paraboloid bowl, which has an analytical solution 
    !!

    use global_mod, only: ip, dp, minimum_allowed_depth, charlen, default_nonlinear_timestepping_method
    use domain_mod, only: domain_type
    use multidomain_mod, only: multidomain_type, setup_multidomain, test_multidomain_mod
    use boundary_mod, only: flather_boundary, transmissive_boundary
    use local_routines
    use timer_mod
    use logging_mod, only: log_output_unit, send_log_output_to_file
    use stop_mod, only: generic_stop
    use iso_c_binding, only: C_DOUBLE !, C_INT, C_LONG
#ifdef COARRAY_PROVIDE_CO_ROUTINES
    use coarray_intrinsic_alternatives, only: co_min
#endif
    !use iso_fortran_env, only: int64, int32
    implicit none

    ! Type holding all domains 
    type(multidomain_type) :: md

    ! Local timing object
    type(timer_type) :: program_timer

    ! Change this to decrease the cell size by mesh_refine (i.e. for convergence testing)
    integer(ip), parameter :: mesh_refine = 4_ip

    ! The global (i.e. outer-domain) time-step in the multidomain 
    real(dp) ::  global_dt ! = (0.23_dp/mesh_refine) * 1.0_dp

    ! Approx timestep between outputs
    real(dp) :: approximate_writeout_frequency = 1.0_dp
    real(dp) :: final_time = 300.0_dp * 1.0_dp

    ! Length/width
    real(dp), parameter, dimension(2):: global_lw = 8000.0_dp * [1.0_dp, 1.0_dp]
    ! Lower-left corner coordinate
    real(dp), parameter, dimension(2):: global_ll = -global_lw/2.0_dp
    ! grid size (number of x/y cells)
    integer(ip), parameter, dimension(2):: global_nx = [100, 100]*mesh_refine + 1

    ! Useful misc variables
    integer(ip):: nd, i, j
    real(dp):: last_write_time

    call program_timer%timer_start('setup')

#ifdef SPHERICAL
    write(log_output_unit,*) 'Code assumes cartesian oordinates, but SPHERICAL is defined'
    call generic_stop
#endif
    
    ! nd domains in this model
    nd = 1
    allocate(md%domains(nd))

    !
    ! Setup basic metadata
    !

    ! Linear domain
    md%domains(1)%lw = global_lw
    md%domains(1)%lower_left =global_ll
    md%domains(1)%nx = global_nx
    md%domains(1)%dx = md%domains(1)%lw/md%domains(1)%nx
    md%domains(1)%dx_refinement_factor = 1.0_dp
    md%domains(1)%timestepping_refinement_factor = 1_ip
    md%domains(1)%timestepping_method = default_nonlinear_timestepping_method 

    !@ Splitting the domain the same way, irrespective of np, improves reproducibility
    ! md%load_balance_file = 'load_balance_partition.txt'
    call get_command_argument(1, md%load_balance_file)

    ! Allocate domains and prepare comms
    call md%setup(extra_halo_buffer=0_ip)

    ! Set initial conditions
    do j = 1, size(md%domains)
        call set_initial_conditions_paraboloid_bowl(md%domains(j))
    end do
    call md%make_initial_conditions_consistent()
    
    ! NOTE: For stability in 'null' regions, we set them to 'high land' that
    ! should be inactive. 
    call md%set_null_regions_to_dry()
   
    write(log_output_unit,*) 'End setup'

    ! Print the gravity-wave CFL limit, to guide timestepping
    do j = 1, size(md%domains)
        write(log_output_unit,*) 'domain: ', j, 'ts: ', &
            md%domains(j)%stationary_timestep_max()
    end do
    global_dt = md%stationary_timestep_max()

    call program_timer%timer_end('setup')
    call program_timer%timer_start('evolve')

#ifdef COARRAY
    sync all
    ! Get the minimum global_dt
    call co_min(global_dt)
#endif
    flush(log_output_unit)

    ! The water will rise (and fall) in the deepest area so we need a timestep 
    ! somewhat smaller than the initial stationary timestep. From experience 
    ! this timestep is stable
    global_dt = 0.71_dp * global_dt * md%domains(1)%cfl ! * 1.0_dp/3.0_dp

    !
    ! Evolve the code
    !

    do while (.true.)
        
        ! IO 
        call program_timer%timer_start('IO')
        call md%write_outputs_and_print_statistics(&
            approximate_writeout_frequency = approximate_writeout_frequency)
        call program_timer%timer_end('IO')

        call md%evolve_one_step(global_dt)

        if (md%domains(1)%time > final_time) exit
    end do

    call program_timer%timer_end('evolve')
    call md%finalise_and_print_timers

    write(log_output_unit,*) ''
    write(log_output_unit, *) 'Program timer'
    write(log_output_unit, *) ''
    call program_timer%print(log_output_unit)
end program
