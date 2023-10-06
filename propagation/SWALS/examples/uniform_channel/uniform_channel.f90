module local_routines 
    !!
    !! Setup uniform channel problem (aligned to grid)
    !!
    use global_mod, only: dp, ip, wall_elevation
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use which_mod, only: which
    implicit none
    
    real(dp), parameter:: bed_slope = 0.01_dp, channel_width=20.0_dp, channel_depth=2.0_dp

    contains 

    subroutine set_initial_conditions_uniform_channel(domain)            
        type(domain_type), intent(inout):: domain
        integer(ip):: i,j
        real(dp):: x, y, cx, cy, slope

        ! Set elevation
        slope = bed_slope 
        cx = (domain%lw(1))*0.5_dp
        cy = (domain%lw(2))*0.5_dp
        do j = 1, domain%nx(2)
            do i = 1,domain%nx(1)
                x = (i-0.5_dp)*domain%dx(1) - cx
                y = (j-0.5_dp)*domain%dx(2) - cy 
                domain%U(i,j,ELV) = y*slope
                if(x > -channel_width/2.0_dp .and. x < channel_width/2.0_dp) then
                    domain%U(i,j,ELV) = domain%U(i,j,ELV) - channel_depth 
                end if
            end do
        end do
        ! Wall boundaries along the sides and top
        domain%U(1,:,ELV) = wall_elevation 
        domain%U(domain%nx(1),:,ELV) = wall_elevation 
        domain%U(:,domain%nx(2),ELV) = wall_elevation 
        domain%U(:,1,ELV) = wall_elevation 

        ! Stage
        domain%U(:,:,STG) = domain%U(:,:,ELV)

        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV))

        ! Friction
        domain%manning_squared = 0.03_dp**2

    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program uniform_channel
    !!
    !! Steady uniform flow in a grid-aligned channel.
    !!
    use global_mod, only: ip, dp, charlen, default_nonlinear_timestepping_method
    use multidomain_mod, only: multidomain_type
    use file_io_mod, only: read_csv_into_array
    use local_routines
    implicit none

    integer(ip):: i, j, nd
    real(dp):: last_write_time, global_dt
    type(multidomain_type):: md

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 60.0_dp
    real(dp), parameter :: final_time = 900.0_dp

    ! Length/width
    real(dp), parameter:: global_lw(2) = [50._dp, 1000._dp] 
    ! lower-left corner coordinate
    real(dp), parameter:: global_ll(2) = [0._dp, 0._dp]
    ! grid size (number of x/y cells)
    integer(ip), parameter:: global_nx(2) = [50, 1000] ! [400, 400] 

    ! Discharge of 1 cubic meter / second per metre width
    real(dp) :: discharge_per_unit_width = 1.0_dp
    real(dp) :: Qin, model_vd, model_d, model_v, theoretical_d, theoretical_v, theoretical_vol, model_vol

    integer(ip), allocatable:: channel_x_indices(:)

    ! This gives the inflow discharge per cell, where it is applied
    Qin = discharge_per_unit_width * (global_lw(1)/global_nx(1))

    nd = 1 ! Number of domains in model
    allocate(md%domains(nd))

    ! Domain Geometry
    md%domains(1)%lw = global_lw
    md%domains(1)%lower_left = global_ll
    md%domains(1)%nx = global_nx
    md%domains(1)%timestepping_method = default_nonlinear_timestepping_method

    ! Allocate multidomain data structures
    call md%setup

    ! Set initial conditions on each domain
    call set_initial_conditions_uniform_channel(md%domains(1))

    ! Find channel indices 
    ! This corresponds to the lowest points on a horizontal cross-section
    call which(md%domains(1)%U(:,2,ELV) < (minval(md%domains(1)%U(:,2,ELV)) + 1.0e-06), channel_x_indices)
    print*, channel_x_indices

    ! Time-step at which we evolve the solution
    global_dt = 0.08_dp

    ! Evolve the solution
    do while (.TRUE.)

        call md%write_outputs_and_print_statistics(approximate_writeout_frequency=approximate_writeout_frequency)

        if (md%domains(1)%time > final_time) exit

        call md%evolve_one_step(global_dt)

        ! Add discharge inflow
        md%domains(1)%U(channel_x_indices, md%domains(1)%nx(2)-1, STG) = &
            md%domains(1)%U(channel_x_indices, md%domains(1)%nx(2)-1, STG) + &
            Qin*global_dt/product(md%domains(1)%dx)

    end do


    ! Mass conservation check
    theoretical_vol = Qin * md%domains(1)%time * channel_width
    model_vol = sum(md%domains(1)%U(:,:,STG) - md%domains(1)%U(:,:,ELV)) * md%domains(1)%dx(1) * md%domains(1)%dx(2)

    ! Analytical solution for uniform slope
    ! vd = discharge_per_unit_width
    ! bed_slope = manning_squared * v * abs(v) / d^(4/3)
    ! bed_slope = manning_squared * discharge_per_unit_width**2 / (d^(10/3))
    ! d = (manning_squared * discharge_per_unit_width**2 / bed_slope )**(3/10)
    i = md%domains(1)%nx(1)/2
    j = md%domains(1)%nx(2)/2    
    model_vd = -1.0_dp * md%domains(1)%U(i,j,3)    
    model_d = md%domains(1)%U(i,j,1) - md%domains(1)%U(i,j,4)
    model_v = model_vd/model_d
    theoretical_d = (md%domains(1)%manning_squared(i,j) * discharge_per_unit_width**2 / bed_slope)**(3.0_dp/10.0_dp)

    ! Cleanup domain
    call md%finalise_and_print_timers

    !
    ! Print test results
    !

    print*, '## Testing ##'
    print*, ''
    print*, '    Theoretical volume: ', theoretical_vol
    print*, '    Model volume: ', model_vol
    print*, '      error: ', theoretical_vol - model_vol
    if( abs(theoretical_vol - model_vol) < 1.0e-06_dp * theoretical_vol) then
        print*, 'PASS'
    else
        stop 'FAIL: Mass conservation error'
    end if
    

    print*, ' '
    print*, '   Theoretical vd: ', discharge_per_unit_width
    print*, '   Model vd (near centre): ', model_vd
    if( abs(discharge_per_unit_width - model_vd) < 1.0e-03_dp*discharge_per_unit_width) then
        print*, 'PASS'
    else
        stop 'FAIL: Flux error'
    end if

    print*, ' '
    print*, '   Theoretical d: ', theoretical_d
    print*, '   Model d: ', model_d
    if( abs(theoretical_d - model_d) < 1.0e-02_dp*theoretical_d) then
        print*, 'PASS'
    else
        stop 'FAIL: Friction or steady flow error'
    end if

end program
