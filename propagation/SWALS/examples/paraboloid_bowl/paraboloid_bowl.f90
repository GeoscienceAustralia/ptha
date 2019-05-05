module local_routines 
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

    end subroutine

end module 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program run_paraboloid_basin

    use global_mod, only: ip, dp, minimum_allowed_depth, charlen
    use domain_mod, only: domain_type
    use multidomain_mod, only: multidomain_type, setup_multidomain, test_multidomain_mod
    use boundary_mod, only: flather_boundary, transmissive_boundary
    use local_routines
    use timer_mod
    use logging_mod, only: log_output_unit, send_log_output_to_file
    use stop_mod, only: generic_stop
    use iso_c_binding, only: C_DOUBLE !, C_INT, C_LONG
#ifdef COARRAY_PROVIDE_CO_ROUTINES
    use coarray_intrinsic_alternatives, only: co_max, co_sum, co_broadcast
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
    real(dp) ::  global_dt = (0.23_dp/mesh_refine)

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
    md%domains(1)%timestepping_method = 'rk2' !'midpoint'

    ! Allocate domains and prepare comms
    call md%setup()

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
            md%domains(j)%linear_timestep_max()
    end do

    call program_timer%timer_end('setup')
    call program_timer%timer_start('evolve')

#ifdef COARRAY
    sync all
    flush(log_output_unit)
#endif

    !
    ! Evolve the code
    !

    ! Trick to get the code to write out just after the first timestep
    last_write_time = -approximate_writeout_frequency
    do while (.true.)
        
        ! IO 
        if(md%domains(1)%time - last_write_time >= approximate_writeout_frequency) then
            call program_timer%timer_start('IO')

            call md%print()
            !call md%domains(1)%print()

            do j = 1, size(md%domains)
                call md%domains(j)%write_to_output_files()
            end do
            last_write_time = last_write_time + approximate_writeout_frequency
            flush(log_output_unit)

            call program_timer%timer_end('IO')
#ifdef COARRAY
            ! This sync can be useful for debugging but is not a good idea in general
            !sync all
#endif
        end if

        call md%evolve_one_step(global_dt)

        if (md%domains(1)%time > final_time) exit
    end do

    call program_timer%timer_end('evolve')

    ! Print out timing info for each
    do i = 1, nd
        write(log_output_unit,*) ''
        write(log_output_unit,*) 'Timer ', i
        write(log_output_unit,*) ''
        call md%domains(i)%timer%print(log_output_unit)
        call md%domains(i)%write_max_quantities()
        call md%domains(i)%finalise()
    end do

    write(log_output_unit, *) ''
    write(log_output_unit, *) 'Multidomain timer'
    write(log_output_unit, *) ''
    call md%timer%print(log_output_unit)

    write(log_output_unit,*) ''
    write(log_output_unit, *) 'Program timer'
    write(log_output_unit, *) ''
    call program_timer%print(log_output_unit)
end program
