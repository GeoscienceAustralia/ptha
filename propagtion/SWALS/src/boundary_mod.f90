
MODULE boundary_mod
    !
    ! Here we define subroutines that update the exterior boundaries
    ! of the domain.
    !
    ! In a typical application the user might provide their own, which
    !  is passed to the procedure pointer 'domain%boundary_subroutine'.
    ! The only requirement on the boundary_subroutine is that it 
    !  has 'domain' as its only input argument.
    ! If access to other data is required for boundary_subroutine, you 
    ! can create the boundary_subroutine in its own module, and put the
    ! required data in the same module.
    !
    USE global_mod, only: dp, ip, charlen, gravity, minimum_allowed_depth, wall_elevation
    USE domain_mod, only: domain_type, STG, UH, VH, ELV

    IMPLICIT NONE

    REAL(dp), PARAMETER, PRIVATE :: HALF_dp = 0.5_dp, ZERO_dp = 0.0_dp, ONE_dp=1.0_dp

    CONTAINS

    !
    ! Simple transmissive boundary condition.
    !
    ! Stage/velocity/elevation at the boundary cells
    ! is set to the values at cells just inside the boundary
    !
    ! This boundary condition is theoretically only valid for supercritical outflow
    ! If used for subcritical flows, it can have a tendency to drift or go unstable,
    ! but not always. Theoretically, the problem would seem to be that there is
    ! nothing to 'nudge' the boundary to the correct solution (e.g. consider open
    ! ocean boundaries in tsunami models -- we might know that the stage will be close
    ! to MSL -- but such information is not embedded in this boundary condition). Mathematically,
    ! with subcritical flows there will be incoming characteristics, so we should
    ! provide some information.
    !
    ! Considering the above, use with caution.
    !
    SUBROUTINE transmissive_boundary(domain)
        TYPE(domain_type), INTENT(INOUT):: domain
        INTEGER(ip):: i, j, k

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        !$OMP DO SCHEDULE(STATIC), COLLAPSE(2)
        DO k = 1, 4
            DO j = 1, domain%nx(2)
                ! Make transmissive boundary conditions
                ! For stability it seems important that the elevation is equal to
                ! the neighbour elevation [likely because this makes the current
                ! approach equivalent to velocity extrapolation]
            
                ! West boundary
                if(domain%boundary_exterior(4)) domain%U(1,j,k) = domain%U(2,j,k) !- domain%U(3,j,k)
                ! East boundary
                if(domain%boundary_exterior(2)) domain%U(domain%nx(1),j,k) = domain%U(domain%nx(1)-1,j,k) !- domain%U(domain%nx(1)-2,j,k)
            END DO
        END DO
        !$OMP END DO NOWAIT

        !$OMP DO SCHEDULE(STATIC), COLLAPSE(2)
        DO k = 1, 4
            DO i = 1, domain%nx(1)
                ! South boundary
                if(domain%boundary_exterior(3)) domain%U(i,1,k) =  domain%U(i,2,k) !- domain%U(i, 3, k)
                ! North boundary
                if(domain%boundary_exterior(1)) domain%U(i,domain%nx(2),k) = domain%U(i,domain%nx(2)-1,k) !- domain%U(i, domain%nx(2) - 2, k)
            END DO
        END DO
        !$OMP END DO
        !$OMP END PARALLEL

    END SUBROUTINE

    !
    ! Apply the boundary function domain%boundary_function at the boundaries.
    ! The latter function must take as input (domain, time, x, y) and 
    ! return a vector of length 4 giving the [stage,uh,vh,elev] values at
    ! boundary location x,y at time t.
    !
    ! If a nonlinear solver is used, then the normal velocity is set equal
    ! to the normal velocity inside the domain, and the transverse velocity
    ! is set to zero.
    !
    ! If the linear solver is used, ony stage and elevation are set, since
    ! with the staggered grid, that is enough to compute the interior UH/VH
    !
    SUBROUTINE boundary_stage_transmissive_normal_momentum(domain)
        TYPE(domain_type), INTENT(INOUT):: domain
        INTEGER(ip):: i, j, k
        REAL(dp):: bc_values(4)

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        !$OMP DO SCHEDULE(STATIC)
        DO j = 1, domain%nx(2)
            ! West boundary
            if(domain%boundary_exterior(4)) then
                i = 1
                bc_values = domain%boundary_function(domain, domain%time, domain%x(i), domain%y(j))
                domain%U(i,j,1) = bc_values(1)
                domain%U(i,j,4) = bc_values(4)
                !domain%U(i,j,2:3) = ZERO_dp
                if(domain%timestepping_method /= 'linear') then
                    domain%U(i,j,2) = domain%U(i+1,j,2) !2.0_dp * domain%U(i+1,j,2) - domain%U(i+2,j,2)
                    domain%U(i,j,3) = ZERO_dp
                end if
            end if

            ! East boundary
            if(domain%boundary_exterior(2)) then
                i = domain%nx(1)
                bc_values = domain%boundary_function(domain, domain%time, domain%x(i), domain%y(j))
                domain%U(i,j,1) = bc_values(1)
                domain%U(i,j,4) = bc_values(4)
                !domain%U(i,j,2:3) = ZERO_dp
                if(domain%timestepping_method /= 'linear') then
                    domain%U(i,j,2) = domain%U(i-1,j,2) !2.0_dp * domain%U(i-1,j,2) - domain%U(i-2, j, 2)
                    domain%U(i,j,3) = ZERO_dp
                end if
            end if

        END DO
        !$OMP END DO

        !$OMP DO SCHEDULE(STATIC)
        DO i = 1, domain%nx(1)
            ! South boundary
            if(domain%boundary_exterior(3)) then
                j = 1
                bc_values = domain%boundary_function(domain, domain%time, domain%x(i), domain%y(j))
                domain%U(i,j,1) = bc_values(1)
                domain%U(i,j,4) = bc_values(4)
                !domain%U(i,j,2:3) = ZERO_dp
                if(domain%timestepping_method /= 'linear') then
                    domain%U(i,j,3) = domain%U(i,j+1,3) !2.0_dp * domain%U(i,j+1,3) - domain%U(i,j+2,3)
                    domain%U(i,j,2) = ZERO_dp
                end if
            end if

            ! North boundary
            if(domain%boundary_exterior(1)) then
                j = domain%nx(2)
                bc_values = domain%boundary_function(domain, domain%time, domain%x(i), domain%y(j))
                domain%U(i,j,1) = bc_values(1)
                domain%U(i,j,4) = bc_values(4)
                !domain%U(i,j,2:3) = ZERO_dp
                if(domain%timestepping_method /= 'linear') then
                    domain%U(i,j,3) = domain%U(i,j-1,3) !2.0_dp * domain%U(i,j-1,3) - domain%U(i,j-2,3)
                    domain%U(i,j,2) = ZERO_dp
                end if
            end if
        END DO
        !$OMP END DO
        !$OMP END PARALLEL

    END SUBROUTINE
   
    !
    ! Flather type boundary condition.
    !
    ! Based on ideas in Blayo and Debreu (2005) Ocean Modelling 9:231-252, their Eqn 26
    !
    ! In practice it seems good when the free surface perturbation is small compared
    ! with the flow depth (e.g. domain of linear shallow water theory), but is partly reflective
    ! in the nonlinear regime. 
    ! This seems consistent with the theory outlined by Blayo and Debreu (2005)
    ! 
    SUBROUTINE flather_boundary(domain)
        TYPE(domain_type), INTENT(INOUT):: domain
        INTEGER(ip):: i, j, k
        REAL(dp) :: sqrt_g_on_di, w1, w2, w3, depth_inside, depth_inside_inv
        REAL(dp) :: stage_outside, u_outside, v_outside, depth_outside
        REAL(dp) :: boundary_temp(4)


        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        !$OMP DO SCHEDULE(STATIC)
        DO j = 1, domain%nx(2)

            ! West edge
            if(domain%boundary_exterior(4)) then
                i = 1
                depth_inside = domain%U(i+1,j,STG) - domain%U(i+1,j,ELV)
              
                ! Get the 'outside' values 
                if(associated(domain%boundary_function) .EQV. .FALSE.) then
                    ! default values
                    stage_outside = max(ZERO_dp, domain%U(i,j,ELV))!ZERO_dp 
                    u_outside = ZERO_dp
                    v_outside = ZERO_dp   
                    depth_outside = stage_outside - domain%U(i,j,ELV) !depth_inside
                else
                    ! use the boundary function
                    boundary_temp = domain%boundary_function(domain, domain%time, domain%x(i), domain%y(j))
                    stage_outside = boundary_temp(1)
                    depth_outside = boundary_temp(1) - boundary_temp(4)
                    if(depth_outside  <= minimum_allowed_depth) then
                        u_outside = ZERO_dp
                        v_outside = ZERO_dp
                    else
                        u_outside = boundary_temp(2)/depth_outside
                        v_outside = boundary_temp(3)/depth_outside
                    end if
                end if

                if((depth_inside > minimum_allowed_depth).and.(depth_outside > minimum_allowed_depth)) then
                    depth_inside_inv = ONE_dp / depth_inside
                    sqrt_g_on_di = sqrt(gravity * depth_inside_inv)
                    ! w1 = u + sqrt(g/depth) * stage_outside -- uses 'outside' info
                    w1 = u_outside + sqrt_g_on_di * stage_outside 
                    ! w2 = v (velocity parallel to boundary)
                    w2 = merge(domain%U(i+1,j,VH) * depth_inside_inv, ZERO_dp, domain%U(i+1,j,UH) < ZERO_dp)
                    ! w3 = u - sqrt(g/depth) * stage_inside -- uses inside info
                    w3 = domain%U(i+1,j,UH) * depth_inside_inv - sqrt_g_on_di * domain%U(i+1,j,STG)

                    !stage = -(w3-w1)/(2*sqrt_g_on_depth_inside)
                    domain%U(i,j,STG) = (-w3 + w1)/(2.0_dp * sqrt_g_on_di)
                    !uh = (w3+w1)/2.*depth_inside
                    domain%U(i,j,UH) = (w3 + w1) * HALF_dp * depth_inside
                    !vh =  w2*depth_inside
                    domain%U(i,j,VH) = w2 * depth_inside
                else
                    domain%U(i,j,STG) = stage_outside
                    domain%U(i,j,UH) = ZERO_dp
                    domain%U(i,j,VH) = ZERO_dp
                endif
            end if

            ! East edge
            if(domain%boundary_exterior(2)) then
                i = domain%nx(1)
                depth_inside = domain%U(i-1,j,STG) - domain%U(i-1,j,ELV)

                ! Get the 'outside' values 
                if(associated(domain%boundary_function) .EQV. .FALSE.) then
                    ! default values
                    stage_outside = max(ZERO_dp, domain%U(i,j,ELV))
                    u_outside = ZERO_dp
                    v_outside = ZERO_dp   
                    depth_outside = stage_outside - domain%U(i,j,ELV) !depth_inside
                else
                    ! use the boundary function
                    boundary_temp = domain%boundary_function(domain, domain%time, domain%x(i), domain%y(j))
                    stage_outside = boundary_temp(1)
                    depth_outside = boundary_temp(1) - boundary_temp(4)
                    if(depth_outside  <= minimum_allowed_depth) then
                        u_outside = ZERO_dp
                        v_outside = ZERO_dp
                    else
                        u_outside = boundary_temp(2)/depth_outside
                        v_outside = boundary_temp(3)/depth_outside
                    end if
                end if

                if((depth_inside > minimum_allowed_depth).AND.(depth_outside > minimum_allowed_depth)) then
                    depth_inside_inv = ONE_dp / depth_inside
                    sqrt_g_on_di = sqrt(gravity * depth_inside_inv)
                    ! w1 = u - sqrt(g/depth) * stage_outside -- uses 'outside' info
                    w1 = u_outside - sqrt_g_on_di * stage_outside 
                    ! w2 = v (velocity parallel to boundary)
                    w2 = merge(domain%U(i-1,j,VH) * depth_inside_inv, ZERO_dp, domain%U(i-1,j,UH) > ZERO_dp)
                    ! w3 = u + sqrt(g/depth) * stage_inside -- uses inside info
                    w3 = domain%U(i-1,j,UH) * depth_inside_inv + sqrt_g_on_di * domain%U(i-1,j,STG)

                    !stage = (w3-w1)/(2*sqrt_g_on_depth_inside)
                    domain%U(i,j,STG) = (w3 - w1)/(2.0_dp * sqrt_g_on_di)
                    !uh = (w3+w1)/2.*depth_inside
                    domain%U(i,j,UH) = (w3 + w1) * HALF_dp * depth_inside
                    !vh =  w2*depth_inside
                    domain%U(i,j,VH) = w2 * depth_inside
                else
                    domain%U(i,j,STG) = stage_outside
                    domain%U(i,j,UH) = ZERO_dp
                    domain%U(i,j,VH) = ZERO_dp
                endif
            end if
        END DO
        !$OMP END DO
        
        !$OMP DO SCHEDULE(STATIC)
        DO i = 1, domain%nx(1)
            ! South edge
            if(domain%boundary_exterior(3)) then
                j = 1
                depth_inside = domain%U(i,j+1,STG) - domain%U(i,j+1,ELV)

                ! Get the 'outside' values 
                if(associated(domain%boundary_function) .EQV. .FALSE.) then
                    ! default values
                    stage_outside = max(ZERO_dp, domain%U(i,j,ELV))
                    u_outside = ZERO_dp
                    v_outside = ZERO_dp   
                    depth_outside = stage_outside - domain%U(i,j,ELV) !depth_inside
                else
                    ! use the boundary function
                    boundary_temp = domain%boundary_function(domain, domain%time, domain%x(i), domain%y(j))
                    stage_outside = boundary_temp(1)
                    depth_outside = boundary_temp(1) - boundary_temp(4)
                    if(depth_outside  <= minimum_allowed_depth) then
                        u_outside = ZERO_dp
                        v_outside = ZERO_dp
                    else
                        u_outside = boundary_temp(2)/depth_outside
                        v_outside = boundary_temp(3)/depth_outside
                    end if
                end if

                if((depth_inside > minimum_allowed_depth).AND.(depth_outside > minimum_allowed_depth)) then
                    depth_inside_inv = ONE_dp / depth_inside
                    sqrt_g_on_di = sqrt(gravity * depth_inside_inv)
                    ! w1 = v + sqrt(g/depth) * stage_outside -- uses 'outside' info
                    w1 = v_outside + sqrt_g_on_di * stage_outside 
                    ! w2 = u (velocity parallel to boundary)
                    w2 = merge(domain%U(i,j+1,UH) * depth_inside_inv, ZERO_dp, domain%U(i,j+1,VH) < ZERO_dp)
                    ! w3 = v - sqrt(g/depth) * stage_inside -- uses inside info
                    w3 = domain%U(i,j+1,VH) * depth_inside_inv - sqrt_g_on_di * domain%U(i,j+1,STG)

                    !stage = -(w3-w1)/(2*sqrt_g_on_depth_inside)
                    domain%U(i,j,STG) = (- w3 + w1)/(2.0_dp * sqrt_g_on_di)
                    !vh = (w3+w1)/2.*depth_inside
                    domain%U(i,j,VH) = (w3 + w1) * HALF_dp * depth_inside
                    !uh =  w2*depth_inside
                    domain%U(i,j,UH) = w2 * depth_inside
                else
                    domain%U(i,j,STG) = stage_outside
                    domain%U(i,j,UH) = ZERO_dp
                    domain%U(i,j,VH) = ZERO_dp
                endif
            end if

            ! North edge
            if(domain%boundary_exterior(1)) then
                j = domain%nx(2)
                depth_inside = domain%U(i,j-1,STG) - domain%U(i,j-1,ELV)

                ! Get the 'outside' values 
                if(associated(domain%boundary_function) .EQV. .FALSE.) then
                    ! default values
                    stage_outside = max(ZERO_dp, domain%U(i,j, ELV))
                    u_outside = ZERO_dp
                    v_outside = ZERO_dp   
                    depth_outside = stage_outside - domain%U(i,j,ELV) !depth_inside
                else
                    ! use the boundary function
                    boundary_temp = domain%boundary_function(domain, domain%time, domain%x(i), domain%y(j))
                    stage_outside = boundary_temp(1)
                    depth_outside = boundary_temp(1) - boundary_temp(4)
                    if(depth_outside  <= minimum_allowed_depth) then
                        u_outside = ZERO_dp
                        v_outside = ZERO_dp
                    else
                        u_outside = boundary_temp(2)/depth_outside
                        v_outside = boundary_temp(3)/depth_outside
                    end if
                end if

                if((depth_inside > minimum_allowed_depth).AND.(depth_outside > minimum_allowed_depth)) then
                    depth_inside_inv = ONE_dp / depth_inside
                    sqrt_g_on_di = sqrt(gravity * depth_inside_inv)
                    ! w1 = v - sqrt(g/depth) * stage_outside -- uses 'outside' info
                    w1 = v_outside - sqrt_g_on_di * stage_outside 
                    ! w2 = u (velocity parallel to boundary)
                    w2 = merge(domain%U(i,j-1,UH) * depth_inside_inv, ZERO_dp, domain%U(i,j-1,VH) > ZERO_dp)
                    ! w3 = v + sqrt(g/depth) * stage_inside -- uses inside info
                    w3 = domain%U(i,j-1,VH) * depth_inside_inv + sqrt_g_on_di * domain%U(i,j-1,STG)

                    !stage = (w3-w1)/(2*sqrt_g_on_depth_inside)
                    domain%U(i,j,STG) = (w3 - w1)/(2.0_dp * sqrt_g_on_di)
                    !vh = (w3+w1)/2.*depth_inside
                    domain%U(i,j,VH) = (w3 + w1) * HALF_dp * depth_inside
                    !uh =  w2*depth_inside
                    domain%U(i,j,UH) = w2 * depth_inside
                else
                    domain%U(i,j,STG) = stage_outside
                    domain%U(i,j,UH) = ZERO_dp
                    domain%U(i,j,VH) = ZERO_dp
                endif
            endif

        END DO
        !$OMP END DO

        !$OMP END PARALLEL
    END SUBROUTINE

    !
    ! Impose Periodic EW boundaries, and reflective NS boundaries
    ! This can be useful in global scale tsunami models with poles
    ! cut off (although it would be better to use transmissive at north 
    !
    SUBROUTINE periodic_EW_reflective_NS(domain)
        TYPE(domain_type), INTENT(INOUT) :: domain
        INTEGER(ip):: j

        ! Set east boundary from west boundary, and vice-versa

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        !$OMP DO SCHEDULE(STATIC)
        DO j = 1, domain%nx(2)
            if(domain%boundary_exterior(4)) then
                domain%U(1, j, :) = domain%U(domain%nx(1)-3, j, :)
                domain%U(2, j, :) = domain%U(domain%nx(1)-2, j, :)
            end if
            if(domain%boundary_exterior(2)) then
                domain%U(domain%nx(1)-1, j, :) = domain%U(3, j, :)
                domain%U(domain%nx(1), j, :) = domain%U(4, j, :)
            end if

            ! Use a reflective NS condition
            if(domain%boundary_exterior(3)) then
                if((j == 1).OR.(j==2)) then
                    domain%U(:,j,ELV) = wall_elevation
                    domain%U(:,j,STG) = wall_elevation
                end if
            endif
            if(domain%boundary_exterior(1)) then
                if((j==domain%nx(2)-1).OR.(j==domain%nx(2))) then
                    domain%U(:,j,ELV) = wall_elevation 
                    domain%U(:,j,STG) = wall_elevation
                end if
            endif

        END DO    
        !$OMP END DO
        !$OMP END PARALLEL
    END SUBROUTINE
END MODULE
