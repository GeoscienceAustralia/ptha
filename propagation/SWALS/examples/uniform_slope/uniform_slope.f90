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
    use multidomain_mod, only: multidomain_type
    use file_io_mod, only: read_csv_into_array
    use local_routines
    implicit none

    integer(ip):: i, j
    real(dp):: last_write_time
    type(multidomain_type):: md

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 60.0_dp
    real(dp) :: final_time

    ! length/width
    real(dp), parameter :: global_lw(2) = [50._dp, 1000._dp] 
    ! lower-left corner coordinate
    real(dp), parameter :: global_ll(2) = [0._dp, 0._dp]
    ! grid size (number of x/y cells)
    integer(ip), parameter :: global_nx(2) = [50, 1000] ! [400, 400] 

    ! Discharge of 1 cubic meter / second per metre width
    real(dp), parameter :: discharge_per_unit_width = 1.0_dp
    real(dp) :: Qin, model_vd, model_d, model_v, theoretical_d, &
        theoretical_v, theoretical_vol, model_vol, initial_vol, ts, error_thresh_vd
    integer(ip) :: nspread, nd

    ! This gives the inflow discharge per cell, where it is applied
    Qin = discharge_per_unit_width * (global_lw(1)/global_nx(1))

    nd = 1 ! Number of domains in model
    allocate(md%domains(nd))

    ! Timestepping method from command line
    call get_command_argument(1, md%domains(1)%timestepping_method)

    ! Geometry
    md%domains(1)%lw = global_lw
    md%domains(1)%lower_left = global_ll
    md%domains(1)%nx = global_nx

    if(md%domains(1)%timestepping_method == 'leapfrog_linear_plus_nonlinear_friction') then
        ! This is not a good numerical method for this problem. But to facilitate testing,
        ! we "nurse it" when setting the initial conditions.
        md%domains(1)%linear_solver_is_truely_linear = .false. ! Would not make sense to use purely-linear on this problem
        final_time = 600.0_dp
        error_thresh_vd = 1.0e-04_dp
    else
        ! The finite volume methods are good for this problem, no issues with dry states, etc.
        final_time = 600.0_dp
        error_thresh_vd = 1.0e-05_dp
    end if

    call md%setup

    call set_initial_conditions_uniform_slope(md%domains(1))

    initial_vol = sum(md%domains(1)%U(:,:,STG) - md%domains(1)%U(:,:,ELV))*product(md%domains(1)%dx)

    if(md%domains(1)%timestepping_method == 'leapfrog_linear_plus_nonlinear_friction') then
        ! The leapfrog solve needs a specified timestep
        ! Because this problem is "hard" for that solver, it turns out we need a fairly
        ! small timestep
        ts = 0.02_dp
    else
        ts = 0.1_dp
    end if

    ! Evolve the solution
    do while (.TRUE.)

        call md%write_outputs_and_print_statistics(approximate_writeout_frequency=approximate_writeout_frequency)

        if (md%domains(1)%time > final_time) exit

        call md%evolve_one_step(ts)

        ! Add discharge inflow
        md%domains(1)%U(2:md%domains(1)%nx(1)-1, md%domains(1)%nx(2)-1, STG) = &
            md%domains(1)%U(2:md%domains(1)%nx(1)-1, md%domains(1)%nx(2)-1, STG) + &
            Qin*ts/product(md%domains(1)%dx)

    end do
    
    theoretical_vol = Qin * md%domains(1)%time * md%domains(1)%lw(1) * (global_nx(1) - 2) * 1.0_dp / global_nx(1) + initial_vol
    model_vol = sum(md%domains(1)%U(:,:,STG) - md%domains(1)%U(:,:,ELV)) * product(md%domains(1)%dx)

    ! Analytical solution
    ! vd = discharge_per_unit_width
    ! bed_slope = manning_squared * v * abs(v) / d^(4/3)
    ! bed_slope = manning_squared * discharge_per_unit_width**2 / (d^(10/3))
    ! d = (manning_squared * discharge_per_unit_width**2 / bed_slope )**(3/10)

    i = md%domains(1)%nx(1)/2
    j = md%domains(1)%nx(2)/2    
    model_vd = -1.0_dp * md%domains(1)%U(i,j,VH)    
    model_d = md%domains(1)%U(i,j,STG) - md%domains(1)%U(i,j,ELV)
    model_v = model_vd/model_d

    theoretical_d = (md%domains(1)%manning_squared(i,j) * discharge_per_unit_width**2 / bed_slope)**(3.0_dp/10.0_dp)

    call md%finalise_and_print_timers

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

end program
