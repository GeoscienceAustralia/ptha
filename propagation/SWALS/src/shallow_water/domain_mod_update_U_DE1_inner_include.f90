! 
! This file contains a code-snippet that is #included in the subroutine "update_U" in domain_mod. It allows us to make the
! variable "friction_power_depth" a parameter. For efficiency it's good for the code to know the power-law power in the friction
! term.
!  

    !subroutine update_U(domain, dt)
    !    !!
    !    !! Update the values of domain%U (i.e. the main flow variables), based on the fluxes and sources in domain
    !    !!
    !    class(domain_type), intent(inout):: domain
    !    real(dp), intent(in):: dt !! Timestep to advance
        
        real(dp):: inv_cell_area_dt, depth, implicit_factor, dt_gravity, fs
        integer(ip):: j, i, kk

        domain%dt_last_update = dt


        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, dt)
        dt_gravity = dt * gravity
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, domain%nx(2)
            ! For spherical coordiantes, cell area changes with y.
            ! For cartesian coordinates this could be moved out of the loop
            inv_cell_area_dt = dt / domain%area_cell_y(j)

            !$OMP SIMD PRIVATE(kk)
            do i = 1, domain%nx(1)
                ! Fluxes
                do kk = 1, 3
                    domain%U(i,j,kk) = domain%U(i,j,kk) - inv_cell_area_dt * ( & 
                        (domain%flux_NS(i, j+1, kk) - domain%flux_NS(i, j, kk)) + &
                        (domain%flux_EW(i+1, j, kk) - domain%flux_EW(i, j, kk) ))
                end do
            end do

            !$OMP SIMD PRIVATE(depth, fs, implicit_factor)
            do i = 1, domain%nx(1)

                ! Velocity clipping
                depth = domain%depth(i,j)

                if (domain%U(i,j,STG) <= (domain%U(i,j,ELV) + minimum_allowed_depth)) then
                    domain%U(i,j,UH) = ZERO_dp 
                    domain%U(i,j,VH) = ZERO_dp 
                else
#ifndef NOFRICTION
                    ! Implicit friction slope update
                    ! U_new = U_last + U_explicit_update - dt_gravity*depth*friction_slope_multiplier*U_new - 
                    !         dt * linear_friction_coeff * U_new

                    ! If we multiply this by UH or VH, we get the associated friction slope term
                    fs = domain%manning_squared(i,j) * &
                        sqrt(domain%velocity(i,j,UH) * domain%velocity(i,j,UH) + &
                             domain%velocity(i,j,VH) * domain%velocity(i,j,VH)) * &
                        max(depth, minimum_allowed_depth)**friction_power_depth
                        !exp(log(max(depth, minimum_allowed_depth))*NEG_SEVEN_ON_THREE_dp)

                    implicit_factor = ONE_dp/(ONE_dp + dt_gravity*max(depth, ZERO_dp)*fs + dt*domain%linear_friction_coeff)
                    !if(domain%U(i,j,STG) <= (domain%U(i,j,ELV) + minimum_allowed_depth)) implicit_factor = ZERO_dp
#else
                    implicit_factor = ONE_dp
#endif
                    ! Pressure gradients 
                    domain%U(i,j,UH) = domain%U(i,j,UH) + inv_cell_area_dt * domain%explicit_source(i, j, UH)
                    ! Friction
                    domain%U(i,j,UH) = domain%U(i,j,UH) * implicit_factor

                    ! Here we add an extra source to take care of the j-1 pressure gradient addition, which
                    ! could not be updated in parallel (since j-1 might be affected by another OMP thread) 
                    domain%U(i,j,VH) = domain%U(i,j,VH) + inv_cell_area_dt * &
                        (domain%explicit_source(i, j, VH) + domain%explicit_source_VH_j_minus_1(i, j+1))
                    ! Friction
                    domain%U(i,j,VH) = domain%U(i,j,VH) * implicit_factor
                end if
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        domain%time = domain%time + dt

!    end subroutine


