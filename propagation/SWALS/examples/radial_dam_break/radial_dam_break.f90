module local_routines 
    use global_mod, only: dp, ip, wall_elevation
    use domain_mod, only: domain_type
    implicit none

    contains 

    subroutine set_initial_conditions_radial_dam(domain)            
        class(domain_type), target, intent(inout):: domain
        integer(ip):: i,j
        real(dp):: x, y, cx, cy, initial_stage_1, initial_stage_2, radius

        initial_stage_1 = 1.0_dp 
        initial_stage_2 = 10.0_dp
        radius = 50.0_dp

        ! Stage
        domain%U(:,:,1) = initial_stage_1 
        domain%MSL_linear = initial_stage_1

        ! Elevation
        domain%U(:,:,4) = 0._dp
        ! Wall boundaries (without boundary conditions)
        domain%U(1,:,4) = 20.0_dp !wall_elevation 
        domain%U(domain%nx(1),:,4) = 20.0_dp !wall_elevation 
        domain%U(:,1,4) = 20.0_dp !wall_elevation
        domain%U(:,domain%nx(2),4) = 20.0_dp !wall_elevation

        ! Ensure stage >= elevation
        domain%U(:,:,1 ) = max(domain%U(:,:,1), domain%U(:,:,4))

        ! Add a stage perturbation
        cx = (domain%lw(1))*0.5
        cy = (domain%lw(2))*0.5
        do i = 1,domain%nx(1)
            do j = 1, domain%nx(2)
                ! Set perturbation based on distance to the centre
                x = (i-0.5_dp)*domain%dx(1) - cx
                y = (j-0.5_dp)*domain%dx(2) - cy 
                if( x*x + y*y < (radius)**2) then
                    domain%u(i,j,1) = initial_stage_2 !10.0_dp
                end if
            end do
        end do

        domain%manning_squared = 0.0_dp

    end subroutine


end module 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program radial_dam_break
    use global_mod, only: ip, dp, charlen
    use domain_mod, only: domain_type
    use file_io_mod, only: read_csv_into_array
    use local_routines
    implicit none

    integer(ip):: i, nsteps
    real(dp):: last_write_time
    type(domain_type):: domain

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 0.2_dp
    real(dp), parameter :: final_time = 2.0_dp * 1

    ! length/width
    real(dp), parameter, dimension(2):: global_lw = [400._dp, 400._dp] 
    ! lower-left corner coordinate
    real(dp), parameter, dimension(2):: global_ll = [0._dp, 0._dp]
    ! grid size (number of x/y cells)
    integer(ip), parameter, dimension(2):: global_nx = [400, 400] * 2 ! [400, 400] 

    ! analytical solution
    real(dp), allocatable :: analytical_solution(:,:)
    character(len=charlen):: analytical_solution_file


    !domain%theta = 1.0_dp
    domain%timestepping_method = 'rk2' !'euler' !'rk2n'

    !domain%compute_fluxes_inner_method='EEC'

    ! Allocate domain
    CALL domain%allocate_quantities(global_lw, global_nx, global_ll)

    ! Call local routine to set initial conditions
    CALL set_initial_conditions_radial_dam(domain)

    ! Trick to get the code to write out just after the first timestep
    last_write_time = -approximate_writeout_frequency

    ! Evolve the code
    do while (.TRUE.)

        if(domain%time - last_write_time >= approximate_writeout_frequency) then

            last_write_time = last_write_time + approximate_writeout_frequency

            call domain%print()
            call domain%write_to_output_files()

            if (domain%time > final_time) then
                exit 
            end if

        end if

        call domain%evolve_one_step()

    end do

    call domain%write_max_quantities()

    call domain%timer%print()

    !! Compare with analytical solution
    !analytical_solution_file = 'Ver_numerical_2.000000.csv'
    !call read_csv_into_array(analytical_solution, analytical_solution_file)
    !! FIXME: Need to add code to automate comparison

    ! Crude check on peak velocity. 7.75 is a numerical value, not a proper
    ! comparison with analytical solution
   
    call domain%compute_depth_and_velocity()
 
    if((abs(maxval(domain%velocity(:,:,2)) - 7.75_dp) < 0.05) .and. &
       (abs(maxval(domain%velocity(:,:,3)) - 7.75_dp) < 0.05) .and. &
       (abs(minval(domain%velocity(:,:,2)) + 7.75_dp) < 0.05) .and. &
       (abs(minval(domain%velocity(:,:,3)) + 7.75_dp) < 0.05) .and. &
       (abs(maxval(domain%velocity(:,:,2)) + minval(domain%velocity(:,:,2))) < 1.0e-03_dp) .and. &
       (abs(maxval(domain%velocity(:,:,3)) + minval(domain%velocity(:,:,3))) < 1.0e-03_dp) ) then
        print*, ' '
        print*, '##############'
        print*, 'PASS'
        print*, '##############'
        print*, ' '
    else
        print*, ' '
        print*, '##############'
        print*, 'FAIL'
        print*, '##############'
        print*, ' '
    end if

    CALL domain%finalise()

end program
