! Compute fluxes for 2D shallow water equations on structured grid
!
! Use an Audusse type method, like ANUGA, but structured
!
! Updated variables:
!    domain%flux_NS, domain%flux_EW, domain%explicit_source, domain%explicit_source_VH_j_minus_1, max_dt_out
!
! @param domain the model domain type for which fluxes etc will be computed
! @param max_dt_out optional real scalar, if provided is set to the max_dt
!     allowed based on CFL limit.
!
!subroutine compute_fluxes_DE1(domain, max_dt_out)

    !type(domain_type), intent(inout) :: domain
    !real(dp), intent(inout) :: max_dt_out

    ! Option to use upwind flux for transverse momentum component, as in 
    ! e.g. "Chertock et al. (2017) Well-balanced schemes for
    ! the shallow water equations with Coriolis forces"
    !logical, parameter :: upwind_transverse_momentum = .false.
    logical, parameter :: upwind_normal_momentum = .false.
    !logical, parameter :: reduced_momentum_diffusion = .true.

    ! wavespeeds
    real(dp):: s_max, s_min, gs_pos, gs_neg, sminsmax
    ! stage/depth at + and - side of edge
    real(dp):: stage_pos, stage_neg, stage_pos_star, stage_neg_star, &
               depth_pos, depth_pos_star, depth_pos_c,&
               depth_neg, depth_neg_star, depth_neg_c,&
               z_half, z_neg, z_pos
    ! velocities and momenta at + and - sides of edge
    real(dp):: u_pos, v_pos, u_neg, v_neg, ud_neg, vd_neg, ud_pos, vd_pos, vel_beta_neg, vel_beta_pos
    ! convenience variables
    integer(ip):: i, j, nx, ny, jlast
    real(dp):: denom, inv_denom, max_speed, max_dt, dx_cfl_inv(2), z_pos_b, z_neg_b, stg_b, stg_a
    real(dp):: bed_slope_pressure_s, bed_slope_pressure_n, bed_slope_pressure_e, bed_slope_pressure_w
    real(dp) :: half_g_hh_edge, precision_factor, half_g_hh_edge_a, half_g_hh_edge_b
    real(dp):: half_cfl, max_dt_inv, dxb, common_multiple, fr2
    character(len=charlen):: timer_name
    real(dp), parameter :: diffusion_scale = ONE_dp
    real(dp), parameter :: EPS = 1.0e-06_dp 
    real(dp), parameter :: EPS_gs = sqrt(gravity * minimum_allowed_depth)  ! EPS relevant to sqrt(gH)
    real(dp), parameter:: half_gravity = HALF_dp * gravity
    integer(ip):: n_ext, loop_work_count, j_low, j_high, my_omp_id, n_omp_threads
    ! Option to use experimental non-conservative pressure gradient term (with .false.)
    logical, parameter :: flux_conservative_pressure = .true.
    
    ! Bottom edge values on 'positive' and 'negative' side (i.e. viewed from j and j-1 respectively)
    ! theta_wd controls the limiting
    real(dp) :: theta_wd_neg_B(domain%nx(1)), theta_wd_pos_B(domain%nx(1)), &
        stage_neg_B(domain%nx(1)), stage_pos_B(domain%nx(1)), &
        depth_neg_B(domain%nx(1)), depth_pos_B(domain%nx(1)), &
        u_neg_B(domain%nx(1)), u_pos_B(domain%nx(1)), &
        v_neg_B(domain%nx(1)), v_pos_B(domain%nx(1))
   
    ! Left edge values on 'positive' and 'negative' side (i.e. viewed from i and i-1 respectively) 
    ! theta_wd controls the limiting
    real(dp) :: theta_wd_neg_L(domain%nx(1)), theta_wd_pos_L(domain%nx(1)), &
        stage_neg_L(domain%nx(1)), stage_pos_L(domain%nx(1)), &
        depth_neg_L(domain%nx(1)), depth_pos_L(domain%nx(1)), &
        u_neg_L(domain%nx(1)), u_pos_L(domain%nx(1)), &
        v_neg_L(domain%nx(1)), v_pos_L(domain%nx(1))

    real(dp) :: bed_j_minus_1(domain%nx(1)), max_dt_inv_work(domain%nx(1)), explicit_source_im1_work(domain%nx(1))

    half_cfl = HALF_dp * domain%cfl
    ny = domain%nx(2)
    nx = domain%nx(1)
    loop_work_count = ny - 2 + 1 ! Number of indices between 'ny' and '2'


    ! Set dt to a high value (it will drop)
    max_dt = domain%maximum_timestep
    ! By computing the inverse we avoid division in the loop
    max_dt_inv = ONE_dp/max_dt


    ! NOTE: Here we manually determine the openmp loop indices, because
    ! some efficiency can be gained by moving through a contiguous chunk
    !
    !
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, half_cfl, nx, ny, loop_work_count, max_dt) REDUCTION(MAX:max_dt_inv)
    max_dt_inv_work = ONE_dp/max_dt
#ifdef NOOPENMP
    ! Serial loop from 2:(ny)
    j_low = 2 
    j_high = ny
#else
    ! Parallel loop from 2:ny
    !
    ! NOTE: In fortran, 
    !     DO i = start, end
    !       .... code here ...
    !     END DO
    ! will not do anything at all if start > end [that's the standard, if
    ! we don't specify a stride other than 1].
    ! This behaviour is important for ensuring the openmp loop sharing
    ! code here works, even if e.g. we have more omp threads than loop
    ! indices. In that case, we will have some threads with j_low >
    ! j_high -- which will not loop -- and that is the desired behaviour.
    !
    my_omp_id = omp_get_thread_num()
    n_omp_threads = omp_get_num_threads()
    ! j_low = lower loop index of the thread my_omp_id
    j_low = nint(loop_work_count * my_omp_id * 1.0_dp / n_omp_threads) + 2
    ! j_high = upper loop index of the thread my_omp_id
    j_high = nint(loop_work_count * (my_omp_id + 1) * 1.0_dp / n_omp_threads) + 2 - 1
#endif

    ! Use this in the loop to reduce computation when the previous j value was j-1.
    ! (which is not always true because of openmp)
    jlast = -HUGE(1_ip)

    ! Main loop
    do j = j_low, j_high

        ! Reset explicit source -- note we only refer to (:,j,:) in the loop, so there is no need
        ! to zero beforehand
        domain%explicit_source(:,j,:) = ZERO_dp
        domain%explicit_source_VH_j_minus_1(:,j) = ZERO_dp
        explicit_source_im1_work = ZERO_dp

        ! Assume distance_left_edge is constant but distance_bottom_edge
        ! might change with y (spherical coordinates). For cartesian coordinates
        ! we could move this outside the loop.
        dx_cfl_inv(1) = ONE_dp/(domain%distance_bottom_edge(j) * half_cfl)
        dx_cfl_inv(2) = ONE_dp/(domain%distance_left_edge(1) * half_cfl)

        ! Get bed at j-1 (might improve cache access) 
        bed_j_minus_1 = domain%U(:,j-1,ELV)

        ! Get bottom edge values
        call get_bottom_edge_values(domain, j, nx, ny, &
            theta_wd_neg_B, theta_wd_pos_B, &
            stage_neg_B, stage_pos_B, &
            depth_neg_B, depth_pos_B, &
            u_neg_B, u_pos_B, &
            v_neg_B, v_pos_B, &
            reuse_gradients = (jlast == j-1))
        ! By setting jlast, we allow the above subroutine call to reuse
        ! limited gradient values from the previous loop iterate
        jlast = j 

        ! Get left edge values
        call get_left_edge_values(domain, j, nx, &
            theta_wd_neg_L, theta_wd_pos_L, &
            stage_neg_L, stage_pos_L, &
            depth_neg_L, depth_pos_L, &
            u_neg_L, u_pos_L, &
            v_neg_L, v_pos_L)

        ! NOTE: much of the loop below can be vectorized, but care will be needed
        ! for the max_dt_inv term
        !$OMP SIMD
        do i = 2, nx
            !
            ! North-South flux computation
            ! Compute flux_NS(i,j) = flux(i,j-1/2)
            !
            ! Cell i,j has North-South flux(i,j+1/2) at flux_NS index i,j+1, 
            !          and North-South flux(i,j-1/2) at flux_NS index i,j
            !

            ! Negative/positive bottom edge values
            stage_neg = stage_neg_B(i)
            stage_pos = stage_pos_B(i)
            depth_neg = depth_neg_B(i)
            depth_pos = depth_pos_B(i)
            u_neg = u_neg_B(i)
            u_pos = u_pos_B(i)
            v_neg = v_neg_B(i)
            v_pos = v_pos_B(i)
            depth_neg_c = domain%depth(i, j-1)
            depth_pos_c = domain%depth(i, j)

            ! Bed elevation
            z_neg = stage_neg - depth_neg
            z_pos = stage_pos - depth_pos
            z_half = max(z_neg, z_pos)

            ! Audusse stage and depths, negative side
            stage_neg_star = max(stage_neg, z_half)
            depth_neg_star = stage_neg_star - z_half

            ! Velocity (in NS direction)
            if(depth_neg_star == ZERO_dp) then
                v_neg = ZERO_dp
                u_neg = ZERO_dp
            end if

            ! Audusse type depth_integrated_velocity correction
            vd_neg = v_neg * depth_neg_star
            ud_neg = u_neg * depth_neg_star

            ! Gravity wave celerity
            gs_neg = sqrt(gravity * depth_neg_star)

            ! Audusse stage and depths, positive side
            stage_pos_star = max(stage_pos, z_half)
            depth_pos_star = stage_pos_star - z_half
            
            ! Velocity (in NS direction)
            if(depth_pos_star == ZERO_dp) then
                v_pos = ZERO_dp
                u_pos = ZERO_dp
            end if
            ! Correct depth-integrated_velocity (Audusse type approach)
            vd_pos = v_pos * depth_pos_star
            ud_pos = u_pos * depth_pos_star

            ! Gravity wave celerity
            gs_pos = sqrt(gravity * depth_pos_star)

            vel_beta_neg = v_neg * advection_beta
            vel_beta_pos = v_pos * advection_beta

            ! Wave-celerities
            s_max = max(max(vel_beta_neg + gs_neg, vel_beta_pos + gs_pos), ZERO_dp)
            s_min = min(min(vel_beta_neg - gs_neg, vel_beta_pos - gs_pos), ZERO_dp)

            denom = s_max - s_min

            if (denom > EPS .and. (gs_pos > EPS_gs .or. gs_neg > EPS_gs)) then
                inv_denom = domain%distance_bottom_edge(j) / denom 
            else
                inv_denom = ZERO_dp
            end if

            sminsmax = s_min * s_max * diffusion_scale

            if(reduced_momentum_diffusion) then
                ! Use this to scale diffusion for uh/vh terms. Too much diffusion can cause artefacts
                fr2 = advection_beta * sqrt(max(0.001_dp, min(1.0_dp, &
                    1.0_dp*(v_neg**2 + v_pos**2 + u_neg**2 + u_pos**2)/(gs_neg**2 + gs_pos**2 + 1.0e-10_dp))))
            else
                fr2 = 1.0_dp
            end if

            ! Advection flux terms 
            domain%flux_NS(i, j, STG) = &
                (s_max * vd_neg - &
                 s_min * vd_pos + &
                 sminsmax * (stage_pos_star - stage_neg_star)) * inv_denom 

            if(.not. upwind_transverse_momentum) then
                domain%flux_NS(i, j, UH) = &
                    (s_max * ud_neg * vel_beta_neg - &
                     s_min * ud_pos * vel_beta_pos + &
                    fr2 * sminsmax* (ud_pos - ud_neg)) * inv_denom
            else
                !domain%flux_NS(i,j,UH) = merge(ud_neg*vel_beta_neg, ud_pos*vel_beta_pos, &
                !    vel_beta_neg+vel_beta_pos > ZERO_dp)
                domain%flux_NS(i,j,UH) = advection_beta * &
                    merge(u_neg, u_pos, domain%flux_NS(i,j,STG) > 0.0_dp) * domain%flux_NS(i,j,STG)
            end if

            if(.not. upwind_normal_momentum) then
                domain%flux_NS(i, j, VH) = &
                    (s_max * (vd_neg * vel_beta_neg) - & 
                     s_min * (vd_pos * vel_beta_pos) + &
                     fr2*sminsmax * (vd_pos - vd_neg))*inv_denom !- &
            else
                domain%flux_NS(i, j, VH) = advection_beta * &
                    merge(v_neg, v_pos, domain%flux_NS(i,j,STG) > 0.0_dp) * domain%flux_NS(i,j,STG)
            end if

            ! Here we put in the gravity/pressure terms. Can try a flux conservative treatment,
            ! or a good old fashioned { g x depth x grad(stage) } term
            if(flux_conservative_pressure) then
                ! Bring the NS gh^2/2 term out of the flux_NS. Note the term is already
                ! multiplied by edge-length. Later we put it in the 'explicit source' terms.
                half_g_hh_edge =  inv_denom * half_gravity * (&
                    s_max * depth_neg_star * depth_neg_star - &
                    s_min * depth_pos_star * depth_pos_star)

                !! 'bed slope' part of pressure gradient term at i,j-1, north side
                bed_slope_pressure_n = half_gravity * ( &
                    (depth_neg_star - depth_neg)*(depth_neg_star + depth_neg) - &
                    (depth_neg + depth_neg_c)*(z_neg - bed_j_minus_1(i)))

                ! To reduce floating point round-off, subtract this factor from half_g_hh_edge 
                ! prior to accumulating in explicit_source type terms. So long as we use the same factor
                ! in the north and south, this will not affect the result. But it can save us some decimal places!
                !! NOTE: This doesn't really offer improvements -- but it might if we combine it inside the half_g_hh_edge
                !! calculation. Consider for later
                !precision_factor =  0.5*(domain%distance_bottom_edge(j) + domain%distance_bottom_edge(j-1)) * &
                !    half_gravity * depth_neg_c * depth_neg_c
             
                !! NOTE regarding OPENMP!!
                !! We cannot naively do:
                ! domain%explicit_source(i, j - 1 ,VH) = domain%explicit_source(i, j-1, VH) + bed_slope_pressure_n + ...
                !! Since j-1 might also be updated on another OMP thread
                !! The solution is below
                !! Other pressure gradient cases refer to i,j, or i-1,j, so should be ok
                domain%explicit_source_VH_j_minus_1(i,j) = &
                    bed_slope_pressure_n * domain%distance_bottom_edge(j) - &
                    half_g_hh_edge 
                    !(half_g_hh_edge - precision_factor) ! Add to explicit_source later

                !! 'bed slope' part of pressure gradient term at i,j, south-side
                bed_slope_pressure_s = half_gravity * ( & 
                    (depth_pos_star - depth_pos) * (depth_pos_star + depth_pos) - &
                    (depth_pos + depth_pos_c)*(z_pos - domain%U(i, j, ELV)) ) 
                
                !precision_factor =  0.5*(domain%distance_bottom_edge(j+1) + domain%distance_bottom_edge(j)) * &
                !    half_gravity * depth_pos_c * depth_pos_c
                domain%explicit_source(i, j, VH) = domain%explicit_source(i, j, VH) - &
                    bed_slope_pressure_s * domain%distance_bottom_edge(j) + &
                    half_g_hh_edge 
                    !(half_g_hh_edge - precision_factor)
            else
                ! Non-conservative treatment
                ! area * g * depth * (stage_top - stage_bottom)/dy
                ! Note "stage_top = stage_neg"
                !   and
                !     "stage_bottom = stage_pos"
                !stg_b = stage_neg
                stg_b = (s_max * s_max * stage_neg + s_min * s_min * stage_pos)/max(s_max*s_max + s_min*s_min, 1.0e-06_dp)
                domain%explicit_source_VH_j_minus_1(i,j) = &
                    -domain%area_cell_y(j-1) * gravity * depth_neg_c * stg_b / domain%distance_left_edge(i)
                !stg_a = stage_pos 
                stg_a = stg_b
                domain%explicit_source(i,j,VH) = domain%explicit_source(i,j,VH) + &
                    domain%area_cell_y(j) * gravity * depth_pos_c * stg_a / domain%distance_left_edge(i)
            end if


            ! Timestep
            max_speed = max(s_max, -s_min)
            if (max_speed > EPS) then
                !max_dt_inv = max(max_dt_inv, max_speed * dx_cfl_inv(2))
                max_dt_inv_work(i) =  max(max_dt_inv_work(i), max_speed * dx_cfl_inv(2))
            end if

        !end do

    !
    !
    ! EW Flux computation
    !
    !
        !!!$OMP SIMD
        !do i = 2, nx
            ! left edge variables 
            stage_neg = stage_neg_L(i)
            stage_pos = stage_pos_L(i)
            depth_neg = depth_neg_L(i)
            depth_pos = depth_pos_L(i)
            u_neg = u_neg_L(i)
            u_pos = u_pos_L(i)
            v_neg = v_neg_L(i)
            v_pos = v_pos_L(i)
            depth_neg_c = domain%depth(i-1, j)
            depth_pos_c = domain%depth(i, j)

            ! Bed elevation
            z_neg = stage_neg - depth_neg
            z_pos = stage_pos - depth_pos
            z_half = max(z_neg, z_pos)

            ! Audusse stage and depths, negative side
            stage_neg_star = max(stage_neg, z_half)
            depth_neg_star = stage_neg_star - z_half

            ! Velocity (in EW direction)
            if(depth_neg_star == ZERO_dp) then
                u_neg = ZERO_dp
                v_neg = ZERO_dp
            end if
            ! Audusse type depth-integrated-velocity correction
            ud_neg = u_neg * depth_neg_star
            vd_neg = v_neg * depth_neg_star

            gs_neg = sqrt(gravity * depth_neg_star)

            ! Audusse stage and depths, positive side
            stage_pos_star = max(stage_pos, z_half)
            depth_pos_star = stage_pos_star - z_half

            ! Velocity (in NS direction)
            if(depth_pos_star == ZERO_dp) then
                u_pos = ZERO_dp
                v_pos = ZERO_dp
            end if
            ud_pos = u_pos * depth_pos_star
            vd_pos = v_pos * depth_pos_star

            gs_pos = sqrt(gravity * depth_pos_star)

            vel_beta_neg = u_neg * advection_beta
            vel_beta_pos = u_pos * advection_beta

            ! Wave-celerities
            s_max = max(max(vel_beta_neg + gs_neg, vel_beta_pos + gs_pos), ZERO_dp)
            s_min = min(min(vel_beta_neg - gs_neg, vel_beta_pos - gs_pos), ZERO_dp)

            denom = s_max - s_min

            if (denom > EPS .and. (gs_pos > EPS_gs .or. gs_neg > EPS_gs)) then
                inv_denom = domain%distance_left_edge(i) / denom 
            else
                inv_denom = ZERO_dp
            end if
            sminsmax = s_min * s_max * diffusion_scale

            if(reduced_momentum_diffusion) then
                ! Use this to scale diffusion for uh/vh terms. Too much diffusion can cause artefacts
                fr2 = advection_beta * sqrt(max(0.001_dp, min(1.0_dp, &
                    1.0_dp*(v_neg**2 + v_pos**2 + u_neg**2 + u_pos**2)/(gs_neg**2 + gs_pos**2 + 1.0e-10_dp))))
            else
                fr2 = 1.0_dp
            end if

            domain%flux_EW(i,j,STG) = &
                (s_max * ud_neg - &
                 s_min * ud_pos + &
                 sminsmax * (stage_pos_star - stage_neg_star)) * inv_denom

            if(.not. upwind_normal_momentum) then
                domain%flux_EW(i,j,UH) = &
                    (s_max * (ud_neg * vel_beta_neg ) - &
                     s_min * (ud_pos * vel_beta_pos ) + &
                     fr2 * sminsmax * (ud_pos - ud_neg) )*inv_denom
            else
                domain%flux_EW(i,j,UH) = advection_beta * &
                    merge(u_neg, u_pos, domain%flux_EW(i,j,STG) > 0.0_dp) * domain%flux_EW(i,j,STG)
            end if
            if(.not. upwind_transverse_momentum) then
                domain%flux_EW(i,j,VH) = &
                    (s_max * vd_neg * vel_beta_neg - &
                     s_min * vd_pos * vel_beta_pos + &
                    fr2 * sminsmax * (vd_pos - vd_neg)) * inv_denom
            else
                !domain%flux_EW(i,j,VH) = merge(vd_neg*vel_beta_neg, vd_pos*vel_beta_pos, &
                !    vel_beta_neg+vel_beta_pos > ZERO_dp)
                domain%flux_EW(i,j,VH) = advection_beta * &
                    merge(v_neg, v_pos, domain%flux_EW(i,j,STG) > 0.0_dp) * domain%flux_EW(i,j,STG)

            end if

            !if(i == 200 .and. j == 10) then
            !    print*, 'STG: ', domain%U(i-1:i+1, j, STG), stage_pos_star
            !    print*, 'VEL: ', domain%velocity(i-1:i+1, j, UH), vel_beta_pos
            !end if


            ! Here we put in the gravity/pressure terms. Can try a flux conservative treatment,
            ! or a good old fashioned { g x depth x grad(stage) } term
            if(flux_conservative_pressure) then
                ! Bring EW pressure gradient term out of flux_EW. Note it is already
                ! multipied by edge-length
                half_g_hh_edge = inv_denom * half_gravity * (&
                    s_max * depth_neg_star * depth_neg_star - &
                    s_min * depth_pos_star * depth_pos_star)

                !! Pressure gradient term at i-1,j -- east side of cell i-1,j
                bed_slope_pressure_e = half_gravity * (&
                    (depth_neg_star - depth_neg)*(depth_neg_star + depth_neg) - &
                    (depth_neg + depth_neg_c)*(z_neg - domain%U(i-1 , j, ELV)))

                ! Add this into explicit_source to reduce floating point round-off errors.  
                ! Mathematically its effect will cancel. 
                ! NOTE: This doesn't really offer improvements, but it might if we integrate it inside
                ! half_g_hh
                !precision_factor =  domain%distance_left_edge(i) * half_gravity * depth_neg_c * depth_neg_c

                !domain%explicit_source(i-1, j, UH) = domain%explicit_source(i-1, j ,UH) + &
                !    bed_slope_pressure_e * domain%distance_left_edge(i) - half_g_hh_edge
                !    !bed_slope_pressure_e * domain%distance_left_edge(i) - (half_g_hh_edge - precision_factor)
                explicit_source_im1_work(i-1) =  &
                    bed_slope_pressure_e * domain%distance_left_edge(i) - half_g_hh_edge
                !    !bed_slope_pressure_e * domain%distance_left_edge(i) - (half_g_hh_edge - precision_factor)

                !! Pressure gradient term at i,j -- west_side
                bed_slope_pressure_w = half_gravity * (&
                    (depth_pos_star - depth_pos) * (depth_pos_star + depth_pos) - &
                    (depth_pos + depth_pos_c)*(z_pos - domain%U(i, j, ELV)) )

                ! Add this into explicit_source to reduce floating point round-off errors.  
                ! Mathematically its effect will cancel
                !precision_factor =  domain%distance_left_edge(i) * half_gravity * depth_pos_c * depth_pos_c

                domain%explicit_source(i, j, UH) = domain%explicit_source(i, j, UH) - &
                    bed_slope_pressure_w * domain%distance_left_edge(i) + half_g_hh_edge
                    !bed_slope_pressure_w * domain%distance_left_edge(i) + (half_g_hh_edge - precision_factor)
            else
                ! Non conservative treatment
                ! area * g * depth * (stage_right - stage_left) / dx
                ! stage_right = stage_neg
                ! stage_pos = stage_left
                dxb = (0.5_dp * (domain%distance_bottom_edge(j+1) + domain%distance_bottom_edge(j)))
                !stg_b = stage_neg
                stg_b = (s_max*s_max*stage_neg + s_min*s_min*stage_pos)/max(s_max*s_max+s_min*s_min, 1.0e-6_dp)
                !domain%explicit_source(i-1,j,UH) = domain%explicit_source(i-1,j,UH) - &
                !    domain%area_cell_y(j) * gravity * depth_neg_c * stg_b / dxb 
                explicit_source_im1_work(i-1) = - &
                    domain%area_cell_y(j) * gravity * depth_neg_c * stg_b / dxb 
                !stg_a = stage_pos
                stg_a = stg_b !(s_max*stage_pos - s_min*stage_neg)/max(s_max-s_min, 1.0e-06_dp)
                domain%explicit_source(i,j,UH) = domain%explicit_source(i,j,UH) + &
                    domain%area_cell_y(j) * gravity * depth_pos_c * stg_a / dxb

            end if

#ifdef SPHERICAL

            !! Source term associated with integrating 'grad' in spherical coordinates
            !! This puts the equations in flux-conservative form as required for finite volumes.
            !! Say we need to integrate 1/R d(F)/dlat over the cell area.
            !! This gives an expression like { R cos(lat) dlon [F_{lat+} - F_{lat-}] }
            !! This is not in finite volume friendly form.
            !! Using the product rule for derivatives we can write
            ! R cos(lat) dlon * [(F)_{lat+} - (F)_{lat-}] = flux - source
            !! where
            ! (flux)  = [ (F R cos(lat) dlon)_{lat+} - (F R cos(lat) dlon)_{lat-} ] 
            ! (source) = F [(R cos(lat) dlon)_{lat+} - (R cos(lat) dlon)_{lat-}]
            if(flux_conservative_pressure) then
                ! 'Source' term associated with treating the pressure in
                ! flux conservative form in spherical coordinates
                domain%explicit_source(i, j, VH) = domain%explicit_source(i, j, VH) + &
                    half_gravity * domain%depth(i, j) * domain%depth(i, j) * &
                    (domain%distance_bottom_edge(j+1) - domain%distance_bottom_edge(j))
            end if

#ifdef CORIOLIS
            ! Coriolis terms
            common_multiple = domain%area_cell_y(j) * domain%coriolis(j)
            domain%explicit_source(i,j,UH) = domain%explicit_source(i,j,UH) + &
                common_multiple * domain%U(i,j,VH)
            domain%explicit_source(i,j,VH) = domain%explicit_source(i,j,VH) - &
                common_multiple * domain%U(i,j,UH)
#endif

            ! Extra Spherical terms from Williamson et al (1992) 
            common_multiple = domain%area_cell_y(j) * &
                domain%tanlat_on_radius_earth(j) * domain%velocity(i,j,UH)
            domain%explicit_source(i,j,UH) = domain%explicit_source(i,j,UH) + &
                common_multiple * domain%U(i,j,VH)
            domain%explicit_source(i,j,VH) = domain%explicit_source(i,j, VH) - &
                common_multiple * domain%U(i,j,UH)
#endif

            !! Timestep
            max_speed = max(s_max, -s_min)
            if (max_speed > EPS) then
                max_dt_inv_work(i) = max(max_dt_inv_work(i), max_speed * dx_cfl_inv(1))
            end if

        end do
        ! Reduction
        max_dt_inv = maxval(max_dt_inv_work)
        ! Apply the explicit source update that could not be done in the SIMD loop
        domain%explicit_source(:,j,UH) = domain%explicit_source(:,j,UH) + explicit_source_im1_work
    end do
    !$OMP END PARALLEL
    
    max_dt = ONE_dp/max_dt_inv
    max_dt_out = max_dt

!end subroutine
