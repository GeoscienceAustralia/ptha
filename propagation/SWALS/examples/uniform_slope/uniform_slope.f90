module local_routines 
    !!
    !! Steady uniform flow down a planar slope
    !!
    use global_mod, only: dp, ip, wall_elevation
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    implicit none
    
    real(dp), parameter:: bed_slope = 0.01_dp


    contains 

    subroutine set_initial_conditions_uniform_slope(domain)            
        class(domain_type), target, intent(inout):: domain
        integer(ip):: i,j
        real(dp):: x, y, cx, cy, slope

        ! Set elevation
        slope = bed_slope 

        cx = (domain%lw(1))*0.5
        cy = (domain%lw(2))*0.5
        do i = 1,domain%nx(1)
            do j = 1, domain%nx(2)
                x = (i-0.5_dp)*domain%dx(1) - cx
                y = (j-0.5_dp)*domain%dx(2) - cy 
                domain%U(i,j,ELV) = y*slope
            end do
        end do

        domain%U(:,:,STG) = domain%U(:,:,ELV)
        if(domain%timestepping_method == 'leapfrog_linear_plus_nonlinear_friction') then
            ! linear_with_nonlinear_friction is fragile for this problem (however we want
            ! to test the friction correctness). So we "help" it by giving initial conditions
            ! close to the steady-state solution. Otherwise shocks develop and the solver
            ! performs very badly. We do not set it exactly
            domain%U(:,:,STG) = domain%U(:,:,STG) + 0.485_dp
            domain%U(:,:,VH) = -0.9999_dp
        end if

        ! Wall boundaries along the sides and top
        domain%U(1,:,ELV) = wall_elevation 
        domain%U(domain%nx(1),:,ELV) = wall_elevation 
        domain%U(:,domain%nx(2),ELV) = wall_elevation 
        domain%U(:,1,ELV) = wall_elevation 

        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV))

        if(domain%timestepping_method == 'leapfrog_linear_plus_nonlinear_friction') then
            ! More nursing of the leapfrog solver.
            ! Ensure the initial condition respects that wall boundaries have no discharge
            do j = 1, domain%nx(2)
                do i = 1, domain%nx(1)
                    if( (domain%U(i,j  ,STG) < domain%U(i,j  ,ELV) + 1.0e-03_dp) .or. &
                        (domain%U(i,j+1,STG) < domain%U(i,j+1,ELV) + 1.0e-03_dp) ) then
                        domain%U(i,j,VH) = 0.0_dp
                    end if
                end do
            end do
        end if

        domain%manning_squared = 0.03_dp**2

    end subroutine


end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program uniform_slope
    !!
    !! Steady uniform flow down a planar slope
    !!
    use global_mod, only: ip, dp, charlen
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use file_io_mod, only: read_csv_into_array
    use local_routines
    implicit none

    integer(ip):: i, j
    real(dp):: last_write_time
    type(domain_type):: domain

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 60.0_dp
    real(dp) :: final_time

    ! length/width
    real(dp), parameter, dimension(2):: global_lw = [50._dp, 1000._dp] 
    ! lower-left corner coordinate
    real(dp), parameter, dimension(2):: global_ll = [0._dp, 0._dp]
    ! grid size (number of x/y cells)
    integer(ip), parameter, dimension(2):: global_nx = [50, 1000] ! [400, 400] 

    ! Discharge of 1 cubic meter / second per metre width
    real(dp), parameter :: discharge_per_unit_width = 1.0_dp
    real(dp) :: Qin, model_vd, model_d, model_v, theoretical_d, &
        theoretical_v, theoretical_vol, model_vol, initial_vol, ts, error_thresh_vd
    integer(ip) :: nspread

    ! This gives the inflow discharge per cell, where it is applied
    Qin = discharge_per_unit_width * (global_lw(1)/global_nx(1))

    call get_command_argument(1, domain%timestepping_method)

    !domain%timestepping_method = 'rk2'
    !domain%timestepping_method = 'leapfrog_linear_plus_nonlinear_friction'

    if(domain%timestepping_method == 'leapfrog_linear_plus_nonlinear_friction') then
        ! This is not a good numerical method for this problem. But to facilitate testing,
        ! we "nurse it" when setting the initial conditions.
        domain%linear_solver_is_truely_linear = .false. ! Would not make sense to use purely-linear on this problem
        final_time = 600.0_dp
        error_thresh_vd = 1.0e-04_dp
    else
        ! The finite volume methods are good for this problem, no issues with dry states, etc.
        final_time = 600.0_dp
        error_thresh_vd = 1.0e-05_dp
    end if

    ! Prevent dry domain from causing enormous initial step
    domain%maximum_timestep = 5.0_dp

    ! Allocate domain
    call domain%allocate_quantities(global_lw, global_nx, global_ll)

    ! call local routine to set initial conditions
    call set_initial_conditions_uniform_slope(domain)

    initial_vol = sum(domain%U(:,:,STG) - domain%U(:,:,ELV))*product(domain%dx)

    ! Trick to get the code to write out just after the first timestep
    last_write_time = -approximate_writeout_frequency

    if(.not. domain%adaptive_timestepping) then
        ! The leapfrog solve needs a specified timestep
        ! Because this problem is "hard" for that solver, it turns out we need a fairly
        ! small timestep
        ts = 0.02_dp 
    end if

    ! Evolve the code
    do while (.true.)

        if(domain%time - last_write_time >= approximate_writeout_frequency) then

            last_write_time = last_write_time + approximate_writeout_frequency

            call domain%print()
            call domain%write_to_output_files()

            if (domain%time > final_time) then
                exit 
            end if

        end if

        ! Evolve 
        if(.not. domain%adaptive_timestepping) then
            call domain%evolve_one_step(timestep=ts)
            ! Add discharge inflow
            domain%U(2:(domain%nx(1)-1), domain%nx(2)-1, STG) = &
                domain%U(2:(domain%nx(1)-1), domain%nx(2)-1, STG) + Qin*ts / product(domain%dx)
        else
            call domain%evolve_one_step()
            ! Add discharge inflow
            domain%U(2:(domain%nx(1)-1), domain%nx(2)-1, STG) = &
                domain%U(2:(domain%nx(1)-1), domain%nx(2)-1, STG) + Qin*domain%evolve_step_dt / product(domain%dx)
        end if



    end do

    call domain%write_max_quantities()

    call domain%timer%print()

    theoretical_vol = Qin * domain%time * domain%lw(1) * (global_nx(1) - 2) * 1.0_dp / global_nx(1) + initial_vol
    model_vol = sum(domain%U(:,:,STG) - domain%U(:,:,ELV)) * product(domain%dx)

    ! Analytical solution
    ! vd = discharge_per_unit_width
    ! bed_slope = manning_squared * v * abs(v) / d^(4/3)
    ! bed_slope = manning_squared * discharge_per_unit_width**2 / (d^(10/3))
    ! d = (manning_squared * discharge_per_unit_width**2 / bed_slope )**(3/10)

    i = domain%nx(1)/2
    j = domain%nx(2)/2    
    model_vd = -1.0_dp * domain%U(i,j,VH)    
    model_d = domain%U(i,j,STG) - domain%U(i,j,ELV)
    model_v = model_vd/model_d

    theoretical_d = (domain%manning_squared(i,j) * discharge_per_unit_width**2 / bed_slope)**(3.0_dp/10.0_dp)

    print*, 'Initial vol: ', initial_vol
    print*, ''
    print*, ''

    print*, '## Testing ##'
    print*, ' '
    print*, '    Theoretical volume: ', theoretical_vol
    print*, '    Model volume: ', model_vol
    print*, '      error: ', theoretical_vol - model_vol
    if( abs(theoretical_vol - model_vol) < 1.0e-06_dp * theoretical_vol) then
        print*, 'PASS'
    else
        print*, 'FAIL: Mass conservation error', theoretical_vol, model_vol
    end if
    print*, ' '
    print*, '   Theoretical vd: ', discharge_per_unit_width
    print*, '   Model vd (near centre): ', model_vd
    if( abs(discharge_per_unit_width - model_vd) < error_thresh_vd*discharge_per_unit_width) then
        print*, 'PASS'
    else
        print*, 'FAIL: Flux error ', discharge_per_unit_width, model_vd
    end if

    print*, ' '
    print*, '   Theoretical d: ', theoretical_d
    print*, '   Model d: ', model_d
    if( abs(theoretical_d - model_d) < 1.0e-04_dp*theoretical_d) then
        print*, 'PASS'
    else
        print*, 'FAIL: Friction or steady flow error', theoretical_d, model_d
    end if

    call domain%finalise()    

end program
