module local_routines 
    !! Test the convergence of the code in a periodic domain
    use global_mod, only: dp, ip, wall_elevation, pi
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use read_raster_mod, only: multi_raster_type
    use logging_mod, only: log_output_unit
    implicit none

    contains 

    subroutine set_initial_conditions(domain)            
        class(domain_type), intent(inout):: domain
        integer(ip):: i, j
        real(dp), allocatable:: x(:), y(:)

        ! Make space for x/y coordinates, at which we will look-up the rasters
        allocate(x(domain%nx(1)), y(domain%nx(1)))
        x = domain%x
        
        ! Set stage and elevation row-by-row.
        do j = 1, domain%nx(2)
            y = domain%y(j)
            domain%U(:,j,ELV) = 2.0_dp - sin(2*pi*x) - cos(2*pi*y)
            domain%U(:,j,STG) = domain%U(:,j,ELV) + 10.0_dp + exp(sin(2*pi*x))*cos(2*pi*y)
            domain%U(:,j,UH) = sin(cos(2*pi*x))*sin(2*pi*y)
            domain%U(:,j,VH) = cos(2*pi*x)*cos(sin(2*pi*y))
        end do

        deallocate(x,y)

        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV) + 1.0e-07_dp)

        write(log_output_unit,*) 'Stage range is: ', minval(domain%U(:,:,STG)), maxval(domain%U(:,:,STG))
        write(log_output_unit,*) 'Elev range is: ', minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))

    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program periodic_convergence
    !! Test the convergence of the code in a periodic domain

    use global_mod, only: ip, dp, default_nonlinear_timestepping_method
    use domain_mod, only: domain_type
    use multidomain_mod, only: multidomain_type
    use timer_mod, only: timer_type
    use logging_mod, only: log_output_unit
    use stop_mod, only: generic_stop
    use local_routines
    implicit none

    ! Type holding all domains 
    type(multidomain_type) :: md

    ! Local timing object
    type(timer_type) :: program_timer

    ! Change this to decrease the cell size by mesh_refine (i.e. for convergence testing)
    character(len=32) :: mesh_refine_input
    integer(ip) :: mesh_refine

    ! Approx timestep between outputs
    real(dp) :: approximate_writeout_frequency = 0.005_dp
    real(dp) :: final_time = 0.06_dp
    real(dp) :: my_dt 

    ! Length/width
    real(dp), parameter, dimension(2):: global_lw = [1.0_dp, 1.0_dp]
    ! Lower-left corner coordinate
    real(dp), parameter, dimension(2):: global_ll = [0.0_dp, 0.0_dp]
    ! grid size (number of x/y cells)
    integer(ip), dimension(2):: global_nx

    ! Number of domains in model; loop variable
    integer(ip):: nd, j

    call get_command_argument(1, mesh_refine_input)
    read(mesh_refine_input, '(I4)') mesh_refine

    my_dt = 4e-04_dp/mesh_refine
    global_nx = [100_ip, 100_ip]*mesh_refine + 1

    ! Set the model name
    md%output_basedir = './OUTPUTS/'

    call program_timer%timer_start('setup')

#ifdef SPHERICAL
    write(log_output_unit,*) 'Code assumes cartesian coordinates, but SPHERICAL is defined'
    call generic_stop
#endif

    ! Set periodic boundary condition
    md%periodic_xs = [0.0_dp, 1.0_dp]
    md%periodic_ys = [0.0_dp, 1.0_dp]
    
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
    md%domains(1)%timestepping_method = default_nonlinear_timestepping_method ! Do not use linear!

    ! Allocate domains and prepare comms
    call md%setup()
    call md%memory_summary()

    ! Set initial conditions
    do j = 1, size(md%domains)
        call set_initial_conditions(md%domains(j))
    end do
    call md%make_initial_conditions_consistent()
    
    ! NOTE: For stability in 'null' regions, we set them to 'high land' that
    ! should be inactive. 
    call md%set_null_regions_to_dry()
   
    write(log_output_unit,*) 'End setup'

    call program_timer%timer_end('setup')
    call program_timer%timer_start('evolve')

#ifdef COARRAY
    sync all
    flush(log_output_unit)
#endif

    !
    ! Evolve the code
    !
    do while (.true.)

        call md%write_outputs_and_print_statistics(approximate_writeout_frequency = 4.0e-04_dp, timing_tol=my_dt/3.0_dp)

        if (md%domains(1)%time > final_time) exit

        call md%evolve_one_step(my_dt)

    end do

    call program_timer%timer_end('evolve')
    call md%finalise_and_print_timers

    write(log_output_unit,*) ''
    write(log_output_unit, *) 'Program timer'
    write(log_output_unit, *) ''
    call program_timer%print(log_output_unit)
end program
