module local_routines 
    !!
    !! Setup a radially symmetric problem with potential flow
    !!
    use global_mod, only: dp, ip, charlen
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use read_raster_mod, only: multi_raster_type
    implicit none

    contains 

    subroutine set_initial_conditions_radial_potential(domain)
        !
        ! This function sets the initial conditions in a domain
        !
        class(domain_type), target, intent(inout):: domain

        real(dp), parameter :: bed_elev = -4000.0_dp
        type(multi_raster_type):: stage_data
        integer(ip) :: i, j
        real(dp), allocatable :: x(:), y(:)

        ! Elevation (flat bed)
        domain%U(:,:,ELV) = bed_elev 

        ! Initial fluxes are zero
        domain%U(:,:,UH:VH) = 0.0_dp

        ! Stage from a file
        call stage_data%initialise([character(len=charlen):: 'initial_condition_file.tif'])
        allocate(x(domain%nx(1)), y(domain%nx(1)))
        x = domain%x
        do j = 1, domain%nx(2)
            y = domain%y(j)
            call stage_data%get_xy(x,y, domain%U(:,j,STG), domain%nx(1), bilinear=1_ip)
        end do

    end subroutine


end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program radial_potential
    !!
    !! Radially symmetric problem for potential flow
    !!
    use global_mod, only: ip, dp, charlen, default_nonlinear_timestepping_method
    use multidomain_mod, only: multidomain_type
    use domain_mod, only: UH, VH
    use local_routines
    use coarray_intrinsic_alternatives, only: swals_mpi_init, swals_mpi_finalize
    implicit none

    integer(ip):: i, nsteps
    real(dp):: last_write_time
    type(multidomain_type) :: md

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 100.0_dp
    real(dp), parameter :: final_time = 100.0

    ! length/width
    real(dp), parameter, dimension(2):: global_lw = 200000.0_dp * [1.0_dp, 1.0_dp]
    ! lower-left corner coordinate
    real(dp), parameter, dimension(2):: global_ll = -global_lw/2.0_dp
    ! grid size (number of x/y cells)
    integer(ip), parameter, dimension(2):: global_nx = [100, 50]*2 + 1 ! Deliberately uneven

    real(dp) :: global_dt = 1.25_dp * 401.0_dp/global_nx(1)

    ! Misc
    integer :: j, nd

    call swals_mpi_init

    nd = 1 ! Number of domains in model
    allocate(md%domains(nd))

    !
    ! Set the domain properties
    !
    md%domains(1)%timestepping_method = 'linear' !default_nonlinear_timestepping_method

    ! Domain Geometry
    md%domains(1)%lw = global_lw
    md%domains(1)%lower_left = global_ll
    md%domains(1)%nx = global_nx
    md%domains(1)%msl_linear = 0.0_dp
    md%domains(1)%use_dispersion = .true.
    md%domains(1)%minimum_nesting_layer_thickness = 20_ip

    ! Output variables to store
    md%domains(1)%time_grids_to_store = [character(len=charlen):: 'stage', 'uh', 'vh']
    md%domains(1)%nontemporal_grids_to_store = [character(len=charlen):: 'max_stage', 'max_speed', 'max_flux', &
        'arrival_time', 'elevation0', 'manning_squared']

    !md%load_balance_file = 'load_balance_partition.txt'

    call md%setup()

    do j = 1, size(md%domains)
        ! Set initial conditions on each domain
        call set_initial_conditions_radial_potential(md%domains(j))
    end do

    call md%make_initial_conditions_consistent()

    ! Evolve the solution
    do while (.TRUE.)

        call md%write_outputs_and_print_statistics(approximate_writeout_frequency=approximate_writeout_frequency)

        if (md%domains(1)%time > final_time) exit

        call md%evolve_one_step(global_dt)

    end do

    call md%finalise_and_print_timers
    call swals_mpi_finalize

end program
