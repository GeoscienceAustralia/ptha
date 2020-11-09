!
! This file includes various timestepping routines. It is #included in the domain class to reduce the complexity there.
!
!



    subroutine one_euler_step(domain, timestep, update_nesting_boundary_flux_integral)
        !! 
        !! Standard forward euler 1st order timestepping scheme.
        !! This is also used as a component of more advanced timestepping schemes.
        !! Argument 'timestep' is optional, but if provided should satisfy the CFL condition
        !!
        type(domain_type), intent(inout):: domain
        real(dp), optional, intent(in):: timestep
            !! The timestep by which the solution is advanced. If not provided, advance by the CFL permitted timestep, computed by
            !! domain%compute_fluxes
        logical, optional, intent(in) :: update_nesting_boundary_flux_integral
            !! If TRUE (default), then call domain%nesting_boundary_flux_integral_tstep inside the routine. Sometimes it is more
            !! straightforward to switch this off, and do the computations separately (e.g in rk2).

        character(len=charlen):: timer_name
        real(dp):: ts
        logical :: nesting_bf

        ! By default, update the nesting_boundary_flux_integral tracker
        if(present(update_nesting_boundary_flux_integral)) then
            nesting_bf = update_nesting_boundary_flux_integral
        else
            nesting_bf = .TRUE.
        end if
            

        !TIMER_START('flux')
        call domain%compute_fluxes(ts)
        !TIMER_STOP('flux')

        !TIMER_START('update')

        ! Store the internally computed max timestep
        domain%max_dt = ts

        if(present(timestep)) then
            ts = timestep
        end if

        call domain%update_U(ts)

        ! Track flux through boundaries
        domain%boundary_flux_evolve_integral = domain%boundary_flux_evolve_integral + &
            ts * sum(domain%boundary_flux_store)
        ! Track flux through 'exterior' boundaries (i.e. not nesting)
        domain%boundary_flux_evolve_integral_exterior = domain%boundary_flux_evolve_integral_exterior + &
            ts * sum(domain%boundary_flux_store_exterior, mask=domain%boundary_exterior)

        !TIMER_STOP('update')

        if(nesting_bf) then 
            ! Update the nesting boundary flux
            call domain%nesting_boundary_flux_integral_tstep(&
                ts,&
                flux_NS=domain%flux_NS, flux_NS_lower_index=1_ip, &
                flux_EW=domain%flux_EW, flux_EW_lower_index=1_ip, &
                var_indices=[STG, VH],&
                flux_already_multiplied_by_dx=.TRUE.)
        end if

        ! Coarray communication, if required (this has been superceeded by the multidomain approach)
        !TIMER_START('partitioned_comms')
        if(domain%use_partitioned_comms) call domain%partitioned_comms%communicate(domain%U)
        !TIMER_STOP('partitioned_comms')

    end subroutine

    subroutine one_rk2_step(domain, timestep)
        !!
        !! Standard 2-step second order timestepping runge-kutta scheme.
        !! Argument timestep is optional, but if provided should satisfy the CFL condition
        !!
        type(domain_type), intent(inout):: domain
        real(dp), optional, intent(in):: timestep !! Advance this far in time

        real(dp):: backup_time, dt_first_step, max_dt_store
        integer(ip):: j, k
        character(len=charlen):: timer_name
        integer, parameter :: var_inds(2) = [STG, VH]

        ! Backup quantities
        backup_time = domain%time
        call domain%backup_quantities()

        !
        ! rk2 involves taking 2 euler steps, then halving the change.
        !

        ! First euler step -- 
        if(present(timestep)) then
            call one_euler_step(domain, timestep, &
                update_nesting_boundary_flux_integral = .FALSE.)
            dt_first_step = timestep
        else
            call one_euler_step(domain, update_nesting_boundary_flux_integral = .FALSE.)
            dt_first_step = domain%max_dt
        end if

        ! Update the nesting boundary flux, by only 1/2 tstep. By keeping it outside of
        ! the euler-step, we can accumulate the result over multiple time-steps without
        ! needing a temporary
        call domain%nesting_boundary_flux_integral_tstep(&
            dt=(dt_first_step*HALF_dp),&
            flux_NS=domain%flux_NS, flux_NS_lower_index=1_ip, &
            flux_EW=domain%flux_EW, flux_EW_lower_index=1_ip, &
            var_indices=var_inds,&
            flux_already_multiplied_by_dx=.TRUE.)


        max_dt_store = domain%max_dt
       
        ! Second euler step with the same timestep 
        call one_euler_step(domain, dt_first_step, &
            update_nesting_boundary_flux_integral=.FALSE.)

        ! Update the nesting boundary flux, by the other 1/2 tstep
        call domain%nesting_boundary_flux_integral_tstep(&
            dt=(dt_first_step*HALF_dp),&
            flux_NS=domain%flux_NS, flux_NS_lower_index=1_ip, &
            flux_EW=domain%flux_EW, flux_EW_lower_index=1_ip, &
            var_indices=var_inds,&
            flux_already_multiplied_by_dx=.TRUE.)


        !TIMER_START('average')

        ! Take average (but allow for openmp)
        !
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        do k = 1, size(domain%backup_U, 3)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, domain%nx(2)
                domain%U(:, j, k) = HALF_dp * (domain%U(:, j, k) +&
                    domain%backup_U(:, j, k))
            end do
            !$OMP END DO NOWAIT
        end do
        !$OMP END PARALLEL

        ! Fix time (since we updated twice) and boundary flux integral
        domain%time = backup_time + dt_first_step
        domain%boundary_flux_evolve_integral = HALF_dp * domain%boundary_flux_evolve_integral
        domain%boundary_flux_evolve_integral_exterior = HALF_dp * domain%boundary_flux_evolve_integral_exterior

        ! We want the CFL timestep that is reported to always be based on the first step -- so force that here
        domain%max_dt = max_dt_store

        !TIMER_STOP('average')

    end subroutine
    

    subroutine one_rk2n_step(domain, timestep)
        !!
        !! An n-step 2nd order runge-kutta timestepping scheme.
        !! This involves less flux calls per timestep advance than rk2.
        !! Argument timestep is optional, but if provided, (timestep/4.0) should satisfy the CFL condition
        !! (because this routine takes (n-1)=4 repeated time-steps of that size, where n=5).
        !! FIXME: Still need to implement nesting boundary flux integral timestepping, if we
        !! want to use this inside a multidomain
        type(domain_type), intent(inout):: domain
        real(dp), optional, intent(in):: timestep !! The timestep. Need to have (timestep/4.0) satisfying the CFL condition

        integer(ip), parameter:: n = rk2n_n_value ! number of substeps to take, must be > 2
        real(dp), parameter:: n_inverse = ONE_dp / (ONE_dp * n)
        real(dp):: backup_time, dt_first_step, reduced_dt, max_dt_store
        real(dp):: backup_flux_integral_exterior, backup_flux_integral
        integer(ip):: j,k
        character(len=charlen):: timer_name

        ! Backup quantities
        backup_time = domain%time
        call domain%backup_quantities()
     
        if(present(timestep)) then 
            ! store timestep
            dt_first_step = timestep/(n-1.0_dp)
            ! first step 
            call one_euler_step(domain, dt_first_step,&
                update_nesting_boundary_flux_integral=.FALSE.)
        else
            ! first step 
            call one_euler_step(domain, update_nesting_boundary_flux_integral=.FALSE.)
            ! store timestep
            dt_first_step = domain%max_dt
        end if

        ! Later in the timestepping we do a subtraction of (1/n) *(the change in flow over (n-1) steps)
        ! In terms of flux tracking, this is as though for those (n-1) steps we only timestepped (n-1)/n of what we did timestep,
        call domain%nesting_boundary_flux_integral_tstep(&
            dt=(dt_first_step*(n-1.0_dp)*n_inverse),& ! Rescaled
            flux_NS=domain%flux_NS, flux_NS_lower_index=1_ip, &
            flux_EW=domain%flux_EW, flux_EW_lower_index=1_ip, &
            var_indices=[STG, VH],&
            flux_already_multiplied_by_dx=.TRUE.)

        max_dt_store = domain%max_dt

        ! Steps 2, n-1
        do j = 2, n-1        
            call one_euler_step(domain, dt_first_step,&
                update_nesting_boundary_flux_integral=.FALSE.)

            call domain%nesting_boundary_flux_integral_tstep(&
                dt=(dt_first_step*(n-1.0_dp)*n_inverse),& !Rescaled
                flux_NS=domain%flux_NS, flux_NS_lower_index=1_ip, &
                flux_EW=domain%flux_EW, flux_EW_lower_index=1_ip, &
                var_indices=[STG, VH],&
                flux_already_multiplied_by_dx=.TRUE.)
        end do

        ! Store 1/n * (original_u - u), which will later be added to the solution [like subtracting 1/n of the
        ! evolved flow for the last (n-1) steps]

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        do k = 1, size(domain%backup_U, 3)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, domain%nx(2)
                domain%backup_U(:,j,k) = (domain%backup_U(:,j, k) - domain%U(:,j,k)) * n_inverse
            end do
            !$OMP END DO NOWAIT
        end do
        !$OMP END PARALLEL

        ! Later when we add backup_U, the effect on the boundary_flux_integral is like subtracting( 1/n*flux_at_this_point)
        domain%boundary_flux_evolve_integral = domain%boundary_flux_evolve_integral*(n-1.0_dp)*n_inverse
        domain%boundary_flux_evolve_integral_exterior = domain%boundary_flux_evolve_integral_exterior*(n-1.0_dp)*n_inverse

        ! Now take one step of duration (n-1)/n * dt.
        ! For this step we flux-track as normal
        reduced_dt = (ONE_dp * n - ONE_dp) * n_inverse * dt_first_step
        call one_euler_step(domain, reduced_dt, &
            update_nesting_boundary_flux_integral=.TRUE.)
        
        !TIMER_START('final_update')

        ! Final update
        ! domain%U = domain%U + domain%backup_U
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        do k = 1, size(domain%backup_U, 3)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, domain%nx(2)
                domain%U(:, j, k) = domain%U(:, j, k) + domain%backup_U(:, j, k)
            end do
            !$OMP END DO NOWAIT
        end do
        !$OMP END PARALLEL

        ! Fix time and timestep (since we updated (n-1)*dt regular timesteps)
        domain%time = backup_time + dt_first_step * (n-1) 

        ! We want the max_dt reported to always be based on the first-sub-step CFL (since that is
        ! what is required for stability, even though at later sub-steps the max_dt might reduce)
        domain%max_dt = max_dt_store

        !TIMER_STOP('final_update')

    end subroutine


    subroutine one_midpoint_step(domain, timestep)
        !!
        !! Another 2nd order timestepping scheme (like the trapezoidal rule).
        !! Argument timestep is optional, but if provided should satisfy the CFL condition.
        !!
        type(domain_type), intent(inout):: domain
        real(dp), optional, intent(in) :: timestep

        real(dp):: backup_time, dt_first_step
        integer(ip):: j, k

        ! Backup quantities
        backup_time = domain%time
        call domain%backup_quantities()
        
        !TIMER_START('flux')
        call domain%compute_fluxes(dt_first_step)
        !TIMER_STOP('flux')

        domain%max_dt = dt_first_step

        !TIMER_START('update')
        ! First euler sub-step
        if(present(timestep)) then

            dt_first_step = timestep
            call domain%update_U(dt_first_step*HALF_dp)
        else
            ! First euler sub-step
            call domain%update_U(dt_first_step*HALF_dp)
        end if
        !TIMER_STOP('update')

        !TIMER_START('partitioned_comms')
        if(domain%use_partitioned_comms) call domain%partitioned_comms%communicate(domain%U)
        !TIMER_STOP('partitioned_comms')

        ! Compute fluxes 
        !TIMER_START('flux')
        call domain%compute_fluxes()
        !TIMER_STOP('flux')

        ! Set U back to backup_U
        !
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        do k = 1, size(domain%backup_U, 3)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, domain%nx(2)
                domain%U(:, j, k) = domain%backup_U(:, j, k)
            end do
            !$OMP END DO NOWAIT
        end do
        !$OMP END PARALLEL

        ! Fix time
        domain%time = backup_time

        ! Update U
        !TIMER_START('update')
        call domain%update_U(dt_first_step)
        !TIMER_STOP('update')

        ! Update the nesting boundary flux
        call domain%nesting_boundary_flux_integral_tstep(&
            dt_first_step,&
            flux_NS=domain%flux_NS, flux_NS_lower_index=1_ip, &
            flux_EW=domain%flux_EW, flux_EW_lower_index=1_ip, &
            var_indices=[STG, VH],&
            flux_already_multiplied_by_dx=.TRUE.)


        !TIMER_START('partitioned_comms')
        if(domain%use_partitioned_comms) call domain%partitioned_comms%communicate(domain%U)
        !TIMER_STOP('partitioned_comms')

        domain%boundary_flux_evolve_integral = sum(domain%boundary_flux_store)*&
            dt_first_step
        domain%boundary_flux_evolve_integral_exterior = sum(domain%boundary_flux_store_exterior, mask=domain%boundary_exterior)*&
            dt_first_step


    end subroutine

    subroutine one_linear_leapfrog_step(domain, dt)
        !! 
        !! Linear shallow water equations leap-frog update.
        !! Update domain%U by timestep dt, using the linear shallow water equations.
        !! Note that unlike the other timestepping routines, this does not require a
        !! prior call to domain%compute_fluxes 
        !!
        type(domain_type), intent(inout):: domain
        real(dp), intent(in):: dt
            !! The timestep to advance. Should remain constant in between repeated calls
            !! to the function (since the numerical method assumes constant timestep)

        if(domain%linear_solver_is_truely_linear) then
            call one_truely_linear_leapfrog_step(domain, dt)
        else
            call one_not_truely_linear_leapfrog_step(domain, dt)
        end if

    end subroutine

    ! Truely-linear leap-frog solver
    subroutine one_truely_linear_leapfrog_step(domain, dt)
        type(domain_type), intent(inout):: domain
        real(dp), intent(in):: dt

        ! Do we represent pressure gradients with a 'truely' linear term g * depth0 * dStage/dx,
        ! or with a nonlinear term g * depth * dStage/dx (i.e. where the 'depth' varies)?
        logical, parameter:: truely_linear = .TRUE.

        ! The linear solver code has become complex [in order to reduce memory footprint, and
        ! include coriolis, while still having openmp work]. So it is moved here. 
#include "domain_mod_linear_solver_include.f90"

    end subroutine

    ! Not-truely-linear leap-frog solver
    subroutine one_not_truely_linear_leapfrog_step(domain, dt)
        type(domain_type), intent(inout):: domain
        real(dp), intent(in):: dt

        ! Do we represent pressure gradients with a 'truely' linear term g * depth0 * dStage/dx,
        ! or with a nonlinear term g * depth * dStage/dx (i.e. where the 'depth' varies)?
        logical, parameter:: truely_linear = .FALSE.

        ! The linear solver code has become complex [in order to reduce memory footprint, and
        ! include coriolis, while still having openmp work]. So it is moved here. 
#include "domain_mod_linear_solver_include.f90"

    end subroutine

    subroutine one_leapfrog_linear_plus_nonlinear_friction_step(domain, dt)
        !! Linear solver is combined with nonlinear friction. This might be useful to allow some dissipation
        !! in very large scale tsunami models.
        type(domain_type), intent(inout):: domain
        real(dp), intent(in):: dt

        if(domain%linear_solver_is_truely_linear) then
            call one_leapfrog_truely_linear_plus_nonlinear_friction_step(domain, dt)
        else
            call one_leapfrog_not_truely_linear_plus_nonlinear_friction_step(domain, dt)
        end if

    end subroutine

    ! Truely-linear leap-frog solver PLUS NONLINEAR FRICTION
    subroutine one_leapfrog_truely_linear_plus_nonlinear_friction_step(domain, dt)
        type(domain_type), intent(inout):: domain
        real(dp), intent(in):: dt
        ! Do we represent pressure gradients with a 'truely' linear term g * depth0 * dStage/dx,
        ! or with a nonlinear term g * depth * dStage/dx (i.e. where the 'depth' varies)?
        logical, parameter:: truely_linear = .true.
        ! The linear solver code has become complex [in order to reduce memory footprint, and
        ! include coriolis, while still having openmp work]. So it is moved here. 

#define LINEAR_PLUS_NONLINEAR_FRICTION
#include "domain_mod_linear_solver_include.f90"
#undef LINEAR_PLUS_NONLINEAR_FRICTION

    end subroutine

    ! Not-truely-linear leap-frog solver PLUS NONLINEAR FRICTION
    subroutine one_leapfrog_not_truely_linear_plus_nonlinear_friction_step(domain, dt)
        type(domain_type), intent(inout):: domain
        real(dp), intent(in):: dt
        ! Do we represent pressure gradients with a 'truely' linear term g * depth0 * dStage/dx,
        ! or with a nonlinear term g * depth * dStage/dx (i.e. where the 'depth' varies)?
        logical, parameter:: truely_linear = .false.
        ! The linear solver code has become complex [in order to reduce memory footprint, and
        ! include coriolis, while still having openmp work]. So it is moved here. 

#define LINEAR_PLUS_NONLINEAR_FRICTION
#include "domain_mod_linear_solver_include.f90"
#undef LINEAR_PLUS_NONLINEAR_FRICTION

    end subroutine

    !
    ! Nonlinear leapfrog
    !
    subroutine one_leapfrog_nonlinear_step(domain, dt)
        type(domain_type), intent(inout):: domain
        real(dp), intent(in):: dt

#define FROUDE_LIMITING
#include "domain_mod_leapfrog_nonlinear_include.f90"
#undef FROUDE_LIMITING

    end subroutine

    subroutine one_cliffs_step(domain, dt)
        !! Use the CLIFFS solver of Elena Tolkova (similar to MOST with a different wet/dry treatment).
        type(domain_type), intent(inout) :: domain
        real(dp), intent(in) :: dt !! Timestep to advance -- should satisfy the CFL condition

        integer(ip) :: i, j
        real(dp) :: dx, dy !, max_dt
        real(dp) :: zeta(domain%nx(2))
        ! Flag to remind us that some computations can be optimized later
        logical, parameter :: update_depth_velocity_celerity = .true.

#ifdef SPHERICAL              
        ! Spherical scaling factor in CLIFFS
        zeta = domain%tanlat_on_radius_earth * 0.125_dp
#else
        zeta = 0.0_dp
#endif

        if(update_depth_velocity_celerity) then
            ! In principle we can re-use these variables if they are already up-to-date. 
            ! FIXME: Do this in a more optimized implementation.

            ! Compute depth and velocity
            call domain%compute_depth_and_velocity(domain%cliffs_minimum_allowed_depth)
            ! Compute celerity -- stored in backup_U(:,:,1)
            !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(domain)
            do j = 1, domain%nx(2)
                domain%backup_U(:,j,1) = sqrt(gravity * domain%depth(:,j))
            end do
            !$OMP END PARALLEL DO
        end if

        ! Loop over the y coordinate
        !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(domain, zeta, dt)
        do j = 2, (domain%nx(2) - 1)
            ! 1D-slice along the x direction
            dx = 0.5_dp * (domain%distance_bottom_edge(j+1) + domain%distance_bottom_edge(j))
            dy = 0.0_dp ! Never used for x-sweep 
            ! Get that tan(lat) term
            call cliffs(dt, 1_ip, j, domain%nx(1), domain%cliffs_minimum_allowed_depth, &
                domain%U(:,:,ELV), &
                ! Celerity stored in backup_U
                domain%backup_U(:,:,1), &
                domain%velocity(:,:,UH), domain%velocity(:,:,VH), dx, dy, &
                zeta, domain%manning_squared, domain%linear_friction_coeff)
        end do
        !$OMP END PARALLEL DO
        

        ! Loop over the x-coordinate. 
        ! NOTE: This may have reduced performance because of all the strided lookups like domain%backup_U(i,:,1).
        !        Consider optimising later on
        !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(domain, zeta, dt)
        do i = 2, (domain%nx(1) - 1)
            ! 1D-slice along the y direction
            dx = 0.0_dp ! Never used for y-sweep
            dy = domain%distance_left_edge(i) 
            call cliffs(dt, 2_ip, i, domain%nx(2), domain%cliffs_minimum_allowed_depth, &
                domain%U(:,:,ELV), &
                ! Celerity stored in backup_U
                domain%backup_U(:,:,1), &
                domain%velocity(:,:,UH), domain%velocity(:,:,VH), dx, dy, &
                zeta, domain%manning_squared, domain%linear_friction_coeff)
        end do
        !$OMP END PARALLEL DO

        ! Update domain%U -- in principle part of this could be skipped if we 
        ! were to immediately were to take another step, except
        ! for some required boundary updates
        !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(domain)
        do j = 1, domain%nx(2)
            domain%depth(:,j) = (domain%backup_U(:,j,1) * domain%backup_U(:,j,1) / gravity)
            domain%U(:,j,STG) = domain%depth(:,j) + domain%U(:,j,ELV)
            domain%U(:,j,UH) = domain%velocity(:,j,UH) * domain%depth(:,j) * &
                merge(1, 0, domain%depth(:,j) >= domain%cliffs_minimum_allowed_depth)
            domain%U(:,j,VH) = domain%velocity(:,j,VH) * domain%depth(:,j) * &
                merge(1, 0, domain%depth(:,j) >= domain%cliffs_minimum_allowed_depth)
        end do
        !$OMP END PARALLEL DO

        domain%time = domain%time + dt

        call domain%update_boundary()
        call domain%apply_forcing(dt)

        domain%dt_last_update = dt

    end subroutine

