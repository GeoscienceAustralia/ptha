    !
    ! This file includes various alternatives to the 'domain%update_U' routine. 
    ! The latter has become complicated in trying to optimize the code.
    !


    subroutine update_U_manning(domain, dt)
        !!
        !! Update the values of domain%U (i.e. the main flow variables), based on the fluxes and sources in domain
        !!
        type(domain_type), intent(inout):: domain
        real(dp), intent(in):: dt !! Timestep to advance

        ! Friction power-law term as parameter makes the code more efficient.
        real(dp), parameter :: friction_power_depth = NEG_SEVEN_ON_THREE_dp

#include "domain_mod_update_U_DE1_inner_include.f90"
    end subroutine

    subroutine update_U_chezy(domain, dt)
        !!
        !! Update the values of domain%U (i.e. the main flow variables), based on the fluxes and sources in domain
        !!
        type(domain_type), intent(inout):: domain
        real(dp), intent(in):: dt !! Timestep to advance

        ! Friction power-law term as parameter makes the code more efficient.
        real(dp), parameter :: friction_power_depth = -2.0_dp

#include "domain_mod_update_U_DE1_inner_include.f90"
    end subroutine




!    !
!    ! Update the values of domain%U (i.e. the main flow variables), based on the fluxes and sources in domain
!    ! Same as update_U but slightly restructured (faster in some cases)
!    !
!    ! @param domain the domain to be updated
!    ! @param dt timestep to advance
!    !
!    subroutine update_U_restructured(domain, dt)
!        class(domain_type), intent(inout):: domain
!        real(dp), intent(in):: dt
!        real(dp):: inv_cell_area_dt, depth, implicit_factor(domain%nx(1)), dt_gravity, fs, depth_neg7on3
!        integer(ip):: j, i, kk
!
!        ! Friction power-law term as parameter makes the code more efficient.
!        real(dp), parameter :: friction_power_depth = NEG_SEVEN_ON_THREE_dp
!
!        if(domain%friction_type /= 'manning') stop "Alternative friction models not implemented here" 
!
!        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, dt)
!        dt_gravity = dt * gravity
!        !$OMP DO SCHEDULE(STATIC)
!        do j = 1, domain%nx(2)
!            ! For spherical coordiantes, cell area changes with y.
!            ! For cartesian coordinates this could be moved out of the loop
!            inv_cell_area_dt = dt / domain%area_cell_y(j)
!
!            !$OMP SIMD
!            do i = 1, domain%nx(1)
!        
!                !! Fluxes
!                do kk = 1, 3
!                    domain%U(i,j,kk) = domain%U(i,j,kk) - inv_cell_area_dt * ( & 
!                        (domain%flux_NS(i, j+1, kk) - domain%flux_NS(i, j, kk)) + &
!                        (domain%flux_EW(i+1, j, kk) - domain%flux_EW(i, j, kk) ))
!                end do
!            end do
!
!            !$OMP SIMD
!            do i = 1, domain%nx(1)
!
!#ifndef NOFRICTION
!                depth = domain%depth(i,j)
!                depth_neg7on3 = max(depth, minimum_allowed_depth)**friction_power_depth
!                !depth_neg7on3 = exp(log(max(depth, minimum_allowed_depth))*NEG_SEVEN_ON_THREE_dp)
!
!                ! Implicit friction slope update
!                ! U_new = U_last + U_explicit_update - dt*depth*friction_slope_multiplier*U_new
!
!                ! If we multiply this by UH or VH, we get the associated friction slope term
!                fs = domain%manning_squared(i,j) * &
!                    sqrt(domain%velocity(i,j,UH) * domain%velocity(i,j,UH) + &
!                         domain%velocity(i,j,VH) * domain%velocity(i,j,VH)) * &
!                    depth_neg7on3
!
!                ! Velocity clipping
!                implicit_factor(i) = merge( ONE_dp/(ONE_dp + dt_gravity*max(depth, ZERO_dp)*fs), ZERO_dp, &
!                    (domain%U(i,j,STG) > (domain%U(i,j,ELV) + minimum_allowed_depth)) )
!#else
!                implicit_factor(i) = ONE_dp
!#endif
!            end do
!
!            !$OMP SIMD
!            do i = 1, domain%nx(1)
!                    ! Pressure gradients 
!                    domain%U(i,j,UH) = domain%U(i,j,UH) + inv_cell_area_dt * domain%explicit_source(i, j, UH)
!                    ! Friction
!                    domain%U(i,j,UH) = domain%U(i,j,UH) * implicit_factor(i)
!
!                    ! Here we add an extra source to take care of the j-1 pressure gradient addition, which
!                    ! could not be updated in parallel (since j-1 might be affected by another OMP thread) 
!                    domain%U(i,j,VH) = domain%U(i,j,VH) + inv_cell_area_dt * &
!                        (domain%explicit_source(i, j, VH) + domain%explicit_source_VH_j_minus_1(i, j+1))
!                    ! Friction
!                    domain%U(i,j,VH) = domain%U(i,j,VH) * implicit_factor(i)
!                !end if
!            end do
!        end do
!        !$OMP END DO
!        !$OMP END PARALLEL
!        domain%time = domain%time + dt
!
!        call domain%update_boundary()
!        call domain%apply_forcing(dt)
!    
!    end subroutine
!    
!    !
!    ! Vectorized version of update_U (i.e. using arrays to avoid the 'i' loop)
!    ! Update the values of domain%U (i.e. the main flow variables), based on the fluxes and sources in domain
!    !
!    ! @param domain the domain to be updated
!    ! @param dt timestep to advance
!    !
!    subroutine update_U_vectorized(domain, dt)
!        class(domain_type), intent(inout):: domain
!        real(dp), intent(in):: dt
!
!        real(dp):: inv_cell_area_dt, implicit_factor(domain%nx(1)), dt_gravity, fs(domain%nx(1))
!        integer(ip):: j, kk
!        ! Friction power-law term as parameter makes the code more efficient.
!        real(dp), parameter :: friction_power_depth = NEG_SEVEN_ON_THREE_dp
!
!        if(domain%friction_type /= 'manning') stop "Alternative friction models not implemented here" 
!
!        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, dt)
!        dt_gravity = dt * gravity
!        !$OMP DO SCHEDULE(GUIDED)
!        do j = 1, domain%nx(2)
!            ! For spherical coordiantes, cell area changes with y.
!            ! For cartesian coordinates this could be moved out of the loop
!            inv_cell_area_dt = dt / domain%area_cell_y(j)
!
!            !! Fluxes
!            do kk = 1, 3
!                domain%U(:,j,kk) = domain%U(:,j,kk) - inv_cell_area_dt * ( & 
!                    (domain%flux_NS(:, j+1, kk) - domain%flux_NS(:, j, kk)) + &
!                    (domain%flux_EW(2:(domain%nx(1)+1), j, kk) - domain%flux_EW(1:domain%nx(1), j, kk) ))
!            end do
!
!#ifndef NOFRICTION
!            ! Implicit friction slope update
!            ! U_new = U_last + U_explicit_update - dt*depth*friction_slope_multiplier*U_new
!
!            !! If we multiply this by UH or VH, we get the associated friction slope term
!            fs = domain%manning_squared(:,j) * &
!                sqrt(domain%velocity(:,j,UH) * domain%velocity(:,j,UH) + &
!                     domain%velocity(:,j,VH) * domain%velocity(:,j,VH)) * &
!                !norm2(domain%velocity(i,j,UH:VH)) * &
!                max(domain%depth(:,j), minimum_allowed_depth)**friction_power_depth
!                !exp(log(max(domain%depth(:,j), minimum_allowed_depth))*NEG_SEVEN_ON_THREE_dp)
!
!            implicit_factor = ONE_dp/(ONE_dp + dt_gravity*max(domain%depth(:,j), ZERO_dp)*fs)
!#else
!            implicit_factor = ONE_dp
!#endif
!            ! Velocity clipping here
!            implicit_factor = merge(implicit_factor, ZERO_dp, domain%U(:,j,STG) > (domain%U(:,j,ELV) + minimum_allowed_depth))
!            ! Pressure gradients 
!            domain%U(:,j,UH) = domain%U(:,j,UH) + inv_cell_area_dt * domain%explicit_source(:, j, UH)
!            ! Friction
!            domain%U(:,j,UH) = domain%U(:,j,UH) * implicit_factor
!
!            ! Here we add an extra source to take care of the j-1 pressure gradient addition, which
!            ! could not be updated in parallel (since j-1 might be affected by another OMP thread) 
!            domain%U(:,j,VH) = domain%U(:,j,VH) + inv_cell_area_dt * &
!                (domain%explicit_source(:, j, VH) + domain%explicit_source_VH_j_minus_1(:, j+1))
!            ! Friction
!            domain%U(:,j,VH) = domain%U(:,j,VH) * implicit_factor
!        end do
!        !$OMP END DO
!        !$OMP END PARALLEL
!        domain%time = domain%time + dt
!
!        call domain%update_boundary()
!        call domain%apply_forcing(dt)
!    
!    end subroutine
!
!
