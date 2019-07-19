module local_routines 
    use global_mod, only: dp, ip, wall_elevation
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    implicit none

    contains 

    subroutine set_initial_conditions_dam(domain, h_upstream, h_downstream)
        class(domain_type), target, intent(inout):: domain
        real(dp), intent(in) :: h_upstream, h_downstream

        integer(ip):: i,j
        real(dp):: x, y, cx, cy, initial_stage_1, initial_stage_2, radius

        initial_stage_1 = h_downstream !0.0001_dp 
        initial_stage_2 = h_upstream !1.0_dp

        ! Stage
        domain%U(:,:,STG) = initial_stage_1 
        domain%MSL_linear = initial_stage_1
        do j = 1, domain%nx(2)
            where(domain%x > 0.0_dp) domain%U(:,j,STG) = initial_stage_2
        end do

        ! Elevation
        domain%U(:,:,ELV) = 0._dp
        !! Wall boundaries (without boundary conditions)
        !domain%U(1,:,4) = 20.0_dp !wall_elevation 
        !domain%U(domain%nx(1),:,4) = 20.0_dp !wall_elevation 
        !domain%U(:,1,4) = 20.0_dp !wall_elevation
        !domain%U(:,domain%nx(2),4) = 20.0_dp !wall_elevation

        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV))

    
        domain%manning_squared = 0.0_dp

    end subroutine


end module 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program dam_break
    use global_mod, only: ip, dp, charlen
    use domain_mod, only: domain_type
    use file_io_mod, only: read_csv_into_array
    use local_routines
    use boundary_mod, only: transmissive_boundary
    implicit none

    integer(ip):: i, nsteps
    real(dp):: last_write_time
    type(domain_type):: domain

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 0.2_dp
    real(dp), parameter :: final_time = 30.0_dp * 1

    ! length/width
    real(dp), parameter, dimension(2):: global_lw = [200._dp, 20._dp] 
    ! lower-left corner coordinate
    real(dp), parameter, dimension(2):: global_ll = -global_lw/2.0_dp
    ! grid size (number of x/y cells)
    integer(ip), parameter, dimension(2):: global_nx = [200, 20] * 4 ! [400, 400] 

    ! analytical solution
    real(dp), allocatable :: analytical_solution(:,:)
    character(len=charlen):: analytical_solution_file, input_char
    real(dp) :: h_upstream, h_downstream

    domain%timestepping_method = 'rk2' 

    ! Get the upstream/downstream initial depth from the command line
    call get_command_argument(1, input_char)
    read(input_char, *) h_upstream
    call get_command_argument(2, input_char)
    read(input_char, *) h_downstream

    ! Allocate domain
    CALL domain%allocate_quantities(global_lw, global_nx, global_ll)

    ! Call local routine to set initial conditions
    CALL set_initial_conditions_dam(domain, h_upstream, h_downstream)

    domain%boundary_subroutine => transmissive_boundary

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

    !call domain%compute_depth_and_velocity()
 

    CALL domain%finalise()

end program
