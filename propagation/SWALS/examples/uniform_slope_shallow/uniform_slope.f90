module local_routines 
    !!
    !! Very shallow steady uniform flow down a uniform slope
    !!
    use global_mod, only: dp, ip, wall_elevation
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    implicit none
    

    contains 

    subroutine set_initial_conditions_uniform_slope(domain, bed_slope)
        class(domain_type), intent(inout):: domain
        real(dp), intent(in) :: bed_slope
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
        ! Wall boundaries along the sides and top
        domain%U(1,:,ELV) = wall_elevation 
        domain%U(domain%nx(1),:,ELV) = wall_elevation 
        domain%U(:,domain%nx(2),ELV) = wall_elevation 
        domain%U(:,1,ELV) = wall_elevation 

        ! Stage
        domain%U(:,:,STG) = domain%U(:,:,ELV)

        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV))

        domain%manning_squared = 0.03_dp**2

    end subroutine


end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program uniform_slope_shallow
    !!
    !! Very shallow steady uniform flow down a uniform slope
    !!

    use global_mod, only: ip, dp, charlen, gravity, default_nonlinear_timestepping_method
    use multidomain_mod, only: multidomain_type
    use file_io_mod, only: read_csv_into_array
    use local_routines
    implicit none

    integer(ip):: i, j, nd
    real(dp):: last_write_time
    type(multidomain_type):: md

    ! Slope of surface. Use higher value to make it supercritical
    real(dp), parameter:: bed_slope = 0.01_dp ! 0.05_dp

    ! Rate of input discharge
    ! 0.00001_dp fails with second order spatial (because limiter pushes it to
    ! first order). 2.0_dp fails with first order spatial
    real(dp), parameter :: discharge_per_unit_width = 0.0001_dp 
    ! At first order spatial accuracy with 'euler' timestepping, the test fails
    ! once the depth is around the elevation drop between cells (e.g. discharge_per_unit_width=2.0_dp)
    ! At second order accuracy, the test fails at much shallower depths.
    logical, parameter:: force_first_order_accurate = .false. 

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 60.0_dp
    real(dp), parameter :: final_time = 1000.0_dp/(discharge_per_unit_width)**(1.0_dp/3.0_dp)
    ! In final_time, the division is a 'rough way' to account for different timescales to equilibrium required. Might not work for
    ! every discharge.

    ! length/width
    real(dp), parameter :: global_lw(2) = [50._dp, 1000._dp] 
    ! lower-left corner coordinate
    real(dp), parameter :: global_ll(2) = [0._dp, 0._dp]
    ! grid size (number of x/y cells)
    integer(ip), parameter :: global_nx(2) = [50, 100] ! [400, 400] 

    ! Discharge of 1 cubic meter / second per metre width
    real(dp) :: Qin, model_vd, model_d, model_v, theoretical_d, theoretical_v, theoretical_vol, model_vol, ts


    ! This gives the inflow discharge per cell, where it is applied
    Qin = discharge_per_unit_width * (global_lw(1)/global_nx(1))

    nd = 1 ! Number of domains in model
    allocate(md%domains(nd))

    ! Geometry
    md%domains(1)%lw = global_lw
    md%domains(1)%lower_left = global_ll
    md%domains(1)%nx = global_nx
    ! Timestepping
    md%domains(1)%timestepping_method = default_nonlinear_timestepping_method !'euler' !'rk2n'

    ! If we force the code to be first-order accurate, we fail on this problem
    if(force_first_order_accurate) md%domains(1)%theta = 0.0_dp

    call md%setup

    ! call local routine to set initial conditions
    call set_initial_conditions_uniform_slope(md%domains(1), bed_slope)

    ts = 2.0_dp ! From experience this timestep is good for rk2

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

    theoretical_vol = Qin * md%domains(1)%time * (md%domains(1)%nx(1) - 2) 
    model_vol = sum(md%domains(1)%U(:,:,STG) - md%domains(1)%U(:,:,ELV)) * product(md%domains(1)%dx)

    print*, 'Cell area: ', md%domains(1)%dx, product(md%domains(1)%dx)
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

    print*, '## Testing ##'
    print*, '    Froude: ', (discharge_per_unit_width/theoretical_d)/sqrt(gravity * theoretical_d)
    print*, ' '
    print*, '    Theoretical volume: ', theoretical_vol
    print*, '    Model volume: ', model_vol
    print*, '      error: ', theoretical_vol - model_vol
    if( abs(theoretical_vol - model_vol) < 1.0e-06_dp * theoretical_vol) then
        print*, 'PASS'
    else
        stop 'FAIL: Mass conservation error'
    end if
    print*, ' '
    print*, '   Theoretical d: ', theoretical_d
    print*, '   Model d: ', model_d
    if( abs(theoretical_d - model_d) < 1.0e-04_dp*theoretical_d) then
        print*, 'PASS'
    else
        stop 'FAIL: Friction or steady flow error'
    end if

    print*, ' '
    print*, '   Theoretical vd: ', discharge_per_unit_width
    print*, '   Model vd (near centre): ', model_vd
    if( abs(discharge_per_unit_width - model_vd) < 1.0e-04_dp*discharge_per_unit_width) then
        print*, 'PASS'
    else
        stop 'FAIL: Flux error '
    end if

end program
