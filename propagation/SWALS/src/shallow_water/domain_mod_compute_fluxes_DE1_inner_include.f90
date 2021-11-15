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
    !logical, parameter :: use_eddy_viscosity = .false.

    ! wavespeeds
    real(dp):: s_max, s_min, gs_pos, gs_neg, sminsmax
    ! stage/depth at + and - side of edge
    real(dp):: stage_pos, stage_neg, stage_pos_star, stage_neg_star, &
               depth_pos, depth_pos_star, depth_pos_c,&
               depth_neg, depth_neg_star, depth_neg_c,&
               z_half, z_neg, z_pos
    ! velocities and momenta at + and - sides of edge
    real(dp):: u_pos, v_pos, u_neg, v_neg, ud_neg, vd_neg, ud_pos, vd_pos, vel_beta_neg, vel_beta_pos, &
        uh_neg, uh_pos, vh_neg, vh_pos
    ! convenience variables
    integer(ip):: i, j, nx, ny, jlast
    real(dp):: denom, inv_denom, max_speed, max_dt, dx_cfl_half_inv(2), z_pos_b, z_neg_b, stg_b, stg_a
    real(dp):: bed_slope_pressure_s, bed_slope_pressure_n, bed_slope_pressure_e, bed_slope_pressure_w
    real(dp):: half_g_hh_edge, precision_factor, half_g_hh_edge_a, half_g_hh_edge_b
    real(dp):: half_cfl, max_dt_inv, dxb, common_multiple, fr2
    real(dp):: diffusion_hll, diffusion_turbulent, diffusion_total, diffusion_max
    character(len=charlen):: timer_name
    real(dp), parameter :: diffusion_scale = ONE_dp
    real(dp), parameter :: EPS = 1.0e-12_dp
    real(dp), parameter :: EPS_gs = sqrt(gravity * minimum_allowed_depth)  ! EPS relevant to sqrt(gH)
    real(dp), parameter:: half_gravity = HALF_dp * gravity
    integer(ip):: n_ext, loop_work_count, j_low, j_high, my_omp_id, n_omp_threads
    ! Option to use experimental non-conservative pressure gradient term (with .false.)
    logical, parameter :: flux_conservative_pressure = .true.
    logical, parameter :: extrapolate_uh_vh = .false.

    ! NS change over a cell of stage/depth/u/v, and the "theta" coefficient controlling limiting
    ! Store values for row 'j' and row 'j-1'
    real(dp) :: dstage_NS_lower(domain%nx(1)), dstage_NS(domain%nx(1)), &
        ddepth_NS_lower(domain%nx(1)), ddepth_NS(domain%nx(1)), &
        du_NS_lower(domain%nx(1)), du_NS(domain%nx(1)), &
        dv_NS_lower(domain%nx(1)), dv_NS(domain%nx(1)), &
        duh_NS_lower(domain%nx(1)), duh_NS(domain%nx(1)), &
        dvh_NS_lower(domain%nx(1)), dvh_NS(domain%nx(1)), &
        theta_wd_NS_lower(domain%nx(1)), theta_wd_NS(domain%nx(1))

    ! EW change over a cell of stage/depth/u/v, and the "theta" coefficient controlling limiting
    ! Only store values for row 'j'
    real(dp) :: dstage_EW(domain%nx(1)), ddepth_EW(domain%nx(1)), du_EW(domain%nx(1)), &
        dv_EW(domain%nx(1)), theta_wd_EW(domain%nx(1)), duh_EW(domain%nx(1)), dvh_EW(domain%nx(1))

    real(dp) :: bed_j_minus_1(domain%nx(1)), max_dt_inv_work(domain%nx(1)), explicit_source_im1_work(domain%nx(1))

    ! Should we reduce the hll numerical diffusion (if possible) when applying the physical diffusion? 
    ! If .true., then if numerical and turbulent have the same sign, only use the one with largest absolute value
    logical, parameter :: reduced_numerical_diffusion = .true. 
    real(dp) :: eddy_visc_j(domain%nx(1)), eddy_visc_jm1(domain%nx(1))

    half_cfl = HALF_dp * domain%cfl
    ny = domain%nx(2)
    nx = domain%nx(1)
    loop_work_count = ny - 2 + 1 ! Number of indices between 'ny' and '2'


    ! Set dt to a high value (it will drop)
    max_dt = domain%maximum_timestep
    ! By computing the inverse we avoid division in the loop
    max_dt_inv = ONE_dp/max_dt

    ! NOTE: Below we manually determine the openmp loop indices, because
    ! some efficiency can be gained by moving through a contiguous chunk

    ! This openmp REDUCTION was buggy in ifort2019 (was worked-around by introducing a scratch array)
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, half_cfl, nx, ny, loop_work_count, max_dt), &
    !$OMP REDUCTION(MAX:max_dt_inv)

    max_dt_inv_work = ONE_dp/max_dt
#ifdef NOOPENMP
    ! Serial loop from 2:(ny)
    j_low = 2
    j_high = ny
    my_omp_id = 0
    n_omp_threads = 1
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
    j_low = nint(loop_work_count * my_omp_id * ONE_dp / n_omp_threads) + 2
    ! j_high = upper loop index of the thread my_omp_id
    j_high = nint(loop_work_count * (my_omp_id + 1) * ONE_dp / n_omp_threads) + 2 - 1

#endif

    ! Pre-fetch the flow values at the NS edge from 'j-1'
    call get_NS_limited_gradient_dx(domain, j_low-1, nx, ny, &
        theta_wd_NS, dstage_NS, ddepth_NS, du_NS, dv_NS, duh_NS, dvh_NS, extrapolate_uh_vh)

    if(use_eddy_viscosity) then
        ! Get eddy viscosity for cells at j-1
        call get_eddy_viscosity(domain, j_low-1, eddy_visc_j)
    end if

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
        dx_cfl_half_inv(1) = ONE_dp/(domain%distance_bottom_edge(j) * half_cfl)
        dx_cfl_half_inv(2) = ONE_dp/(domain%distance_left_edge(1) * half_cfl)

        ! Get bed at j-1 (might improve cache access)
        bed_j_minus_1 = domain%U(:,j-1,ELV)

        if(use_eddy_viscosity) then
            ! Get the eddy viscosity, reusing value from previous j
            eddy_visc_jm1 = eddy_visc_j
            call get_eddy_viscosity(domain, j, eddy_visc_j)
        end if

        ! Get the NS change in stage, depth, u-vel, v-vel, at the row 'j-1'
        ! We can reuse the old values since the 'j' loop is ordered
        dstage_NS_lower = dstage_NS
        ddepth_NS_lower = ddepth_NS
        du_NS_lower = du_NS
        dv_NS_lower = dv_NS
        duh_NS_lower = duh_NS
        dvh_NS_lower = dvh_NS
        theta_wd_NS_lower = theta_wd_NS

        ! Get the NS change in stage, depth, u-vel, v-vel, at the row 'j'
        call get_NS_limited_gradient_dx(domain, j, nx, ny, &
            theta_wd_NS, dstage_NS, ddepth_NS, du_NS, dv_NS, duh_NS, dvh_NS, extrapolate_uh_vh)

        !$OMP SIMD PRIVATE(stage_neg, stage_pos, depth_neg, depth_pos, u_neg, u_pos, v_neg, v_pos, depth_neg_c, depth_pos_c, &
        !$OMP vd_neg, ud_neg, vd_pos, ud_pos, diffusion_turbulent, diffusion_hll, diffusion_total, diffusion_max, &
        !$OMP z_neg, z_pos, z_half, stage_neg_star, depth_neg_star, stage_pos_star, depth_pos_star, gs_neg, gs_pos, vel_beta_neg, &
        !$OMP vel_beta_pos, s_min, s_max, inv_denom, denom, sminsmax, fr2, half_g_hh_edge, bed_slope_pressure_e, &
        !$OMP bed_slope_pressure_w, dxb, stg_b, stg_a, common_multiple, max_speed)
        do i = 2, nx
            !
            ! North-South flux computation
            ! Compute flux_NS(i,j) = flux(i,j-1/2)
            !
            ! Cell i,j has North-South flux(i,j+1/2) at flux_NS index i,j+1,
            !          and North-South flux(i,j-1/2) at flux_NS index i,j
            !

            ! Negative/positive bottom edge values
            stage_neg = domain%U(i, j-1, STG) + HALF_dp * dstage_NS_lower(i)
            stage_pos = domain%U(i, j  , STG) - HALF_dp * dstage_NS(i)
            depth_neg = domain%depth(i, j-1) + HALF_dp * ddepth_NS_lower(i)
            depth_pos = domain%depth(i, j  ) - HALF_dp * ddepth_NS(i)
            u_neg = domain%velocity(i, j-1, UH) + HALF_dp * du_NS_lower(i)
            u_pos = domain%velocity(i, j  , UH) - HALF_dp * du_NS(i)
            v_neg = domain%velocity(i, j-1, VH) + HALF_dp * dv_NS_lower(i)
            v_pos = domain%velocity(i, j  , VH) - HALF_dp * dv_NS(i)
            !uh_neg = domain%U(i, j-1, UH) + HALF_dp * duh_NS_lower(i)
            !uh_pos = domain%U(i, j  , UH) - HALF_dp * duh_NS(i)
            !vh_neg = domain%U(i, j-1, VH) + HALF_dp * dvh_NS_lower(i)
            !vh_pos = domain%U(i, j  , VH) - HALF_dp * dvh_NS(i)
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
            !gs_neg = sqrt(gravity * max(depth_neg_c, depth_pos_c)) ! DEBUG

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
            !gs_pos = gs_neg !sqrt(gravity * depth_pos_c) ! DEBUG

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
                fr2 = advection_beta * sqrt(max(0.001_dp, min(ONE_dp, &
                    (v_neg**2 + v_pos**2 + u_neg**2 + u_pos**2)/(gs_neg**2 + gs_pos**2 + 1.0e-10_dp))))
            else
                fr2 = ONE_dp
            end if

            ! Advection flux terms
            domain%flux_NS(i, j, STG) = &
                (s_max * vd_neg - &
                 s_min * vd_pos + &
                 sminsmax * (stage_pos_star - stage_neg_star)) * inv_denom

            if(.not. use_eddy_viscosity) then

                ! Momentum fluxes without turbulent diffusion 

                ! UH term
                if(.not. upwind_transverse_momentum) then
                     domain%flux_NS(i, j, UH) = &
                         (s_max * ud_neg * vel_beta_neg - &
                          s_min * ud_pos * vel_beta_pos + &
                         fr2 * sminsmax * (ud_pos - ud_neg)) * inv_denom
                 else
                     domain%flux_NS(i,j,UH) = advection_beta * &
                         merge(u_neg, u_pos, domain%flux_NS(i,j,STG) > ZERO_dp) * domain%flux_NS(i,j,STG)
                 end if

                 ! VH term
                 if(.not. upwind_normal_momentum) then
                     domain%flux_NS(i, j, VH) = &
                         (s_max * vd_neg * vel_beta_neg - &
                          s_min * vd_pos * vel_beta_pos + &
                          fr2 * sminsmax * (vd_pos - vd_neg))*inv_denom
                 else
                     domain%flux_NS(i, j, VH) = advection_beta * &
                         merge(v_neg, v_pos, domain%flux_NS(i,j,STG) > ZERO_dp) * domain%flux_NS(i,j,STG)
                 end if

             else
                ! Momentum fluxes with turbulent diffusion 

                ! UH term
                if(.not. upwind_transverse_momentum) then
                    domain%flux_NS(i, j, UH) = &
                        (s_max * ud_neg * vel_beta_neg - &
                         s_min * ud_pos * vel_beta_pos ) * inv_denom
                    diffusion_hll =  fr2 * sminsmax * (ud_pos - ud_neg) * inv_denom
                else
                    domain%flux_NS(i,j,UH) = advection_beta * &
                        merge(u_neg, u_pos, domain%flux_NS(i,j,STG) > ZERO_dp) * domain%flux_NS(i,j,STG)
                    diffusion_hll = ZERO_dp
                end if
                !
                ! Turbulent diffusion
                ! Note the eddy-viscosity has been clipped based on explicit time-step constraints, as in e.g. SWASH and other
                ! well-known solvers that use explicit eddy-visc.
                diffusion_turbulent = - &
                    ! { eddy_visc * depth * dx * du/dy} = "turbulent_diffusion" * dx (geometric mean for depth)
                    HALF_dp * (eddy_visc_j(i) + eddy_visc_jm1(i)) * &
                    sqrt(depth_neg_star*depth_pos_star) * domain%distance_bottom_edge(j) * &
                        ((domain%velocity(i,j, UH) - domain%velocity(i, j-1, UH))/domain%distance_left_edge(i) )

                ! Compute diffusive terms.
                diffusion_max = max(abs(diffusion_turbulent), abs(diffusion_hll))
                ! Sum the diffusions -- unless "reduced_numerical_diffusion == .true." and both terms are of the same sign,
                ! in which case we take the max of their absolute value
                diffusion_total = merge(diffusion_turbulent + diffusion_hll, sign(diffusion_max, diffusion_turbulent), &
                    .not. ((diffusion_turbulent*diffusion_hll > 0) .and. reduced_numerical_diffusion))

                domain%flux_NS(i,j,UH) = domain%flux_NS(i,j,UH) + diffusion_total

                ! VH term
                if(.not. upwind_normal_momentum) then
                    domain%flux_NS(i, j, VH) = &
                        (s_max * vd_neg * vel_beta_neg - & 
                         s_min * vd_pos * vel_beta_pos ) * inv_denom
                    diffusion_hll = fr2 * sminsmax * (vd_pos - vd_neg)*inv_denom 
                else
                    domain%flux_NS(i, j, VH) = advection_beta * &
                        merge(v_neg, v_pos, domain%flux_NS(i,j,STG) > ZERO_dp) * domain%flux_NS(i,j,STG)
                    diffusion_hll = ZERO_dp
                end if
                !
                ! Turbulent diffusion
                diffusion_turbulent = - &
                    ! { eddy_visc * depth * dx * du/dy} = "turbulent_diffusion" * dx (geometric mean for depth)
                    HALF_dp * (eddy_visc_j(i) + eddy_visc_jm1(i)) * &
                    sqrt(depth_neg_star*depth_pos_star) * domain%distance_bottom_edge(j) * &
                        ((domain%velocity(i,j, VH) - domain%velocity(i, j-1, VH))/domain%distance_left_edge(i) )

                ! Compute diffusive terms.
                diffusion_max = max(abs(diffusion_turbulent), abs(diffusion_hll))
                ! Sum the diffusions -- unless "reduced_numerical_diffusion == .true. " and both terms are of the same sign,
                ! in which case we take the max of their absolute value
                diffusion_total = merge(diffusion_turbulent + diffusion_hll, sign(diffusion_max, diffusion_turbulent), &
                    .not. ((diffusion_turbulent*diffusion_hll > 0) .and. reduced_numerical_diffusion))

                domain%flux_NS(i,j,VH) = domain%flux_NS(i,j,VH) + diffusion_total
            end if

            ! Here we put in the gravity/pressure terms. Can try a flux conservative treatment,
            ! or a good old fashioned { g x depth x grad(stage) } term
            if(flux_conservative_pressure) then
                ! Bring the NS gh^2/2 term out of the flux_NS. Note the term is already
                ! multiplied by edge-length. Later we put it in the 'explicit source' terms.
                half_g_hh_edge =  inv_denom * half_gravity * (&
                    s_max * depth_neg_star * depth_neg_star - &
                    s_min * depth_pos_star * depth_pos_star)

                ! 'bed slope' part of pressure gradient term at i,j-1, north side
                bed_slope_pressure_n = half_gravity * ( &
                    (depth_neg_star - depth_neg)*(depth_neg_star + depth_neg) - &
                    (depth_neg + depth_neg_c)*(z_neg - bed_j_minus_1(i)))

                ! NOTE regarding OPENMP !
                ! We cannot naively do:
                ! domain%explicit_source(i, j - 1 ,VH) = domain%explicit_source(i, j-1, VH) + bed_slope_pressure_n + ...
                ! Since j-1 might also be updated on another OMP thread
                ! The solution is below
                ! Other pressure gradient cases refer to i,j, or i-1,j, so should be ok
                domain%explicit_source_VH_j_minus_1(i,j) = &
                    bed_slope_pressure_n * domain%distance_bottom_edge(j) - &
                    half_g_hh_edge

                ! 'bed slope' part of pressure gradient term at i,j, south-side
                bed_slope_pressure_s = half_gravity * ( &
                    (depth_pos_star - depth_pos) * (depth_pos_star + depth_pos) - &
                    (depth_pos + depth_pos_c)*(z_pos - domain%U(i, j, ELV)) )

                domain%explicit_source(i, j, VH) = domain%explicit_source(i, j, VH) - &
                    bed_slope_pressure_s * domain%distance_bottom_edge(j) + &
                    half_g_hh_edge
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
                stg_a = stg_b
                domain%explicit_source(i,j,VH) = domain%explicit_source(i,j,VH) + &
                    domain%area_cell_y(j) * gravity * depth_pos_c * stg_a / domain%distance_left_edge(i)
            end if


            ! Timestep
            max_speed = max(s_max, -s_min)
            if (max_speed > EPS) then
                max_dt_inv_work(i) =  max(max_dt_inv_work(i), max_speed * dx_cfl_half_inv(2))
            end if

        end do

        !
        !
        ! EW Flux computation
        !
        !

        ! Compute flux_EW(i,j) = flux(i-1/2,j)
        !
        ! Cell i,j has East-West flux(i+1/2,j) at flux_EW index i+1,j
        !          and East-West flux(i-1/2,j) at flux_EW index i,j

        ! Get EW change in stage, depth, u, v, as well as the "theta" limiter coefficient
        call get_EW_limited_gradient_dx(domain, j, nx, ny, &
            theta_wd_EW, dstage_EW, ddepth_EW, du_EW, dv_EW, duh_EW, dvh_EW, extrapolate_uh_vh)

        !$OMP SIMD PRIVATE(stage_neg, stage_pos, depth_neg, depth_pos, u_neg, u_pos, v_neg, v_pos, depth_neg_c, depth_pos_c, &
        !$OMP vd_neg, ud_neg, vd_pos, ud_pos, &
        !$OMP z_neg, z_pos, z_half, stage_neg_star, depth_neg_star, stage_pos_star, depth_pos_star, gs_neg, gs_pos, vel_beta_neg, &
        !$OMP vel_beta_pos, s_min, s_max, inv_denom, denom, sminsmax, fr2, half_g_hh_edge, bed_slope_pressure_e, &
        !$OMP bed_slope_pressure_w, dxb, stg_b, stg_a, common_multiple, max_speed)
        do i = 2, nx

            ! left edge variables
            stage_neg = domain%U(i-1,j,STG) + HALF_dp * dstage_EW(i-1)
            stage_pos = domain%U(i  ,j,STG) - HALF_dp * dstage_EW(i)
            depth_neg = domain%depth(i-1,j) + HALF_dp * ddepth_EW(i-1)
            depth_pos = domain%depth(i  ,j) - HALF_dp * ddepth_EW(i)
            u_neg = domain%velocity(i-1,j,UH) + HALF_dp * du_EW(i-1)
            u_pos = domain%velocity(i  ,j,UH) - HALF_dp * du_EW(i)
            v_neg = domain%velocity(i-1,j,VH) + HALF_dp * dv_EW(i-1)
            v_pos = domain%velocity(i  ,j,VH) - HALF_dp * dv_EW(i)
            !uh_neg = domain%U(i-1,j,UH) + HALF_dp * duh_EW(i-1)
            !uh_pos = domain%U(i  ,j,UH) - HALF_dp * duh_EW(i)
            !vh_neg = domain%U(i-1,j,VH) + HALF_dp * dvh_EW(i-1)
            !vh_pos = domain%U(i  ,j,VH) - HALF_dp * dvh_EW(i)
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
            !gs_neg = sqrt(gravity * max(depth_neg_c, depth_pos_c)) ! DEBUG

            ! Audusse stage and depths, positive side
            stage_pos_star = max(stage_pos, z_half)
            depth_pos_star = stage_pos_star - z_half

            ! Velocity (in EW direction)
            if(depth_pos_star == ZERO_dp) then
                u_pos = ZERO_dp
                v_pos = ZERO_dp
            end if
            ud_pos = u_pos * depth_pos_star
            vd_pos = v_pos * depth_pos_star

            gs_pos = sqrt(gravity * depth_pos_star)
            !gs_pos = gs_neg !sqrt(gravity * depth_pos_c) ! DEBUG

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
                fr2 = advection_beta * sqrt(max(0.001_dp, min(ONE_dp, &
                    (v_neg**2 + v_pos**2 + u_neg**2 + u_pos**2)/(gs_neg**2 + gs_pos**2 + 1.0e-10_dp))))
            else
                fr2 = ONE_dp
            end if

            domain%flux_EW(i,j,STG) = &
                (s_max * ud_neg - &
                 s_min * ud_pos + &
                 sminsmax * (stage_pos_star - stage_neg_star)) * inv_denom

            if(.not. use_eddy_viscosity) then

                ! Momentum fluxes without turbulent diffusion

                ! UH term
                if(.not. upwind_normal_momentum) then
                    domain%flux_EW(i,j,UH) = &
                        (s_max * ud_neg * vel_beta_neg  - &
                         s_min * ud_pos * vel_beta_pos + &
                         fr2 * sminsmax * (ud_pos - ud_neg) )*inv_denom
                else
                    domain%flux_EW(i,j,UH) = advection_beta * &
                        merge(u_neg, u_pos, domain%flux_EW(i,j,STG) > ZERO_dp) * domain%flux_EW(i,j,STG)
                end if

                ! VH term
                if(.not. upwind_transverse_momentum) then
                    domain%flux_EW(i,j,VH) = &
                        (s_max * vd_neg * vel_beta_neg - &
                         s_min * vd_pos * vel_beta_pos + &
                        fr2 * sminsmax * (vd_pos - vd_neg)) * inv_denom
                else
                    domain%flux_EW(i,j,VH) = advection_beta * &
                        merge(v_neg, v_pos, domain%flux_EW(i,j,STG) > ZERO_dp) * domain%flux_EW(i,j,STG)
                end if

            else

                ! Momentum fluxes with turbulent diffusion

                ! UH term
                if(.not. upwind_normal_momentum) then
                    domain%flux_EW(i,j,UH) = &
                        (s_max * ud_neg * vel_beta_neg  - &
                         s_min * ud_pos * vel_beta_pos )*inv_denom
                     diffusion_hll = fr2 * sminsmax * (ud_pos - ud_neg) * inv_denom
                else
                    domain%flux_EW(i,j,UH) = advection_beta * &
                        merge(u_neg, u_pos, domain%flux_EW(i,j,STG) > ZERO_dp) * domain%flux_EW(i,j,STG)
                    diffusion_hll = ZERO_dp
                end if
                !
                ! Turbulent diffusion
                diffusion_turbulent = - &
                    ! { eddy_visc * depth * dy * dU/dx } = "turbulent_diffusion" * dy (geometric mean for depth)
                    HALF_dp * (eddy_visc_j(i) + eddy_visc_j(i-1)) * &
                    sqrt(depth_neg_star*depth_pos_star) * domain%distance_left_edge(i) * &
                        ((domain%velocity(i,j, UH) - domain%velocity(i-1, j, UH))/&
                         (HALF_dp * (domain%distance_bottom_edge(j) + domain%distance_bottom_edge(j+1))))

                ! Compute diffusive terms.
                diffusion_max = max(abs(diffusion_turbulent), abs(diffusion_hll))
                ! Sum the diffusions -- unless "reduced_numerical_diffusion == .true." and both terms are of the same sign,
                ! in which case we take the max of their absolute value
                diffusion_total = merge(diffusion_turbulent + diffusion_hll, sign(diffusion_max, diffusion_turbulent), &
                    .not. ((diffusion_turbulent*diffusion_hll > 0) .and. reduced_numerical_diffusion))

                domain%flux_EW(i,j,UH) = domain%flux_EW(i,j,UH) + diffusion_total

                ! VH term
                if(.not. upwind_transverse_momentum) then
                    domain%flux_EW(i,j,VH) = &
                        (s_max * vd_neg * vel_beta_neg - &
                         s_min * vd_pos * vel_beta_pos ) * inv_denom
                    diffusion_hll =  fr2 * sminsmax * (vd_pos - vd_neg) * inv_denom
                else
                    domain%flux_EW(i,j,VH) = advection_beta * &
                        merge(v_neg, v_pos, domain%flux_EW(i,j,STG) > ZERO_dp) * domain%flux_EW(i,j,STG)
                    diffusion_hll = ZERO_dp
                end if
                !
                ! Turbulent diffusion, constant eddy viscosity
                diffusion_turbulent = - &
                    ! { eddy_visc * depth * dy * du/dx} = "turbulent_diffusion" * dy (geometric mean for depth)
                    HALF_dp * (eddy_visc_j(i) + eddy_visc_j(i-1)) * &
                    sqrt(depth_neg_star*depth_pos_star) * domain%distance_left_edge(i) * &
                        ((domain%velocity(i,j, VH) - domain%velocity(i-1, j, VH))/&
                         (HALF_dp * (domain%distance_bottom_edge(j) + domain%distance_bottom_edge(j+1))))

                ! Compute diffusive terms.
                diffusion_max = max(abs(diffusion_turbulent), abs(diffusion_hll))
                ! Sum the diffusions -- unless "reduced_numerical_diffusion == .true. " and both terms are of the same sign,
                ! in which case we take the max of their absolute value
                diffusion_total = merge(diffusion_turbulent + diffusion_hll, sign(diffusion_max, diffusion_turbulent), &
                    .not. ((diffusion_turbulent*diffusion_hll > 0) .and. reduced_numerical_diffusion))

                domain%flux_EW(i,j,VH) = domain%flux_EW(i,j,VH) + diffusion_total

            end if

            ! Here we put in the gravity/pressure terms. Can try a flux conservative treatment,
            ! or a good old fashioned { g x depth x grad(stage) } term
            if(flux_conservative_pressure) then
                ! Bring EW pressure gradient term out of flux_EW. Note it is already
                ! multipied by edge-length
                half_g_hh_edge = inv_denom * half_gravity * (&
                    s_max * depth_neg_star * depth_neg_star - &
                    s_min * depth_pos_star * depth_pos_star)

                ! Pressure gradient term at i-1,j -- east side of cell i-1,j
                bed_slope_pressure_e = half_gravity * (&
                    (depth_neg_star - depth_neg)*(depth_neg_star + depth_neg) - &
                    (depth_neg + depth_neg_c)*(z_neg - domain%U(i-1 , j, ELV)))

                ! To make this section OMP SIMD, we cannot update domain%explicit_source(i-1,j)
                ! Thus we introduce a term to store its value, which is added outside the loop
                explicit_source_im1_work(i-1) =  &
                    bed_slope_pressure_e * domain%distance_left_edge(i) - half_g_hh_edge

                ! Pressure gradient term at i,j -- west_side
                bed_slope_pressure_w = half_gravity * (&
                    (depth_pos_star - depth_pos) * (depth_pos_star + depth_pos) - &
                    (depth_pos + depth_pos_c)*(z_pos - domain%U(i, j, ELV)) )

                domain%explicit_source(i, j, UH) = domain%explicit_source(i, j, UH) - &
                    bed_slope_pressure_w * domain%distance_left_edge(i) + half_g_hh_edge
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

            ! Source term associated with integrating 'grad' in spherical coordinates
            ! This puts the equations in flux-conservative form as required for finite volumes.
            ! Say we need to integrate 1/R d(F)/dlat over the cell area.
            ! This gives an expression like { R cos(lat) dlon [F_{lat+} - F_{lat-}] }
            ! This is not in finite volume friendly form.
            ! Using the product rule for derivatives we can write
            ! R cos(lat) dlon * [(F)_{lat+} - (F)_{lat-}] = flux - source
            ! where
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

            ! Timestep
            max_speed = max(s_max, -s_min)
            if (max_speed > EPS) then
                max_dt_inv_work(i) = max(max_dt_inv_work(i), max_speed * dx_cfl_half_inv(1))
            end if

        end do

        ! Apply the explicit source update that would break vectorization potential of the above loop
        domain%explicit_source(:,j,UH) = domain%explicit_source(:,j,UH) + explicit_source_im1_work
    end do
    max_dt_inv = maxval(max_dt_inv_work)
    !$OMP END PARALLEL

    max_dt = ONE_dp/max_dt_inv
    max_dt_out = max_dt


    contains

    !
    ! Get the NS gradients for stage, depth, u-vel, v-vel, at row j
    !
    subroutine get_NS_limited_gradient_dx(domain, j, nx, ny, &
            theta_wd_NS, dstage_NS, ddepth_NS, du_NS, dv_NS, duh_NS, dvh_NS, extrapolate_uh_vh)

        type(domain_type), intent(in) :: domain
        integer(ip), intent(in) :: j, nx, ny
        real(dp), intent(inout) :: theta_wd_NS(nx), dstage_NS(nx), ddepth_NS(nx), du_NS(nx), dv_NS(nx), duh_NS(nx), dvh_NS(nx)
        logical, intent(in) :: extrapolate_uh_vh

        integer(ip) :: i
        real(dp) :: mindep, maxdep, theta_local

        if(j > 1 .and. j < ny) then
            ! Typical case

            ! limiter coefficient
            !$OMP SIMD PRIVATE(mindep, maxdep, theta_local)
            do i = 1, nx
                mindep = min(domain%depth(i, j-1), domain%depth(i,j), domain%depth(i,j+1)) - minimum_allowed_depth
                maxdep = max(domain%depth(i, j-1), domain%depth(i,j), domain%depth(i,j+1)) + &
                    limiter_coef3*minimum_allowed_depth
                theta_local = limiter_coef4 * (mindep/maxdep - limiter_coef1)
                theta_wd_NS(i) = max(domain%theta * min(ONE_dp, theta_local), ZERO_dp)
            end do

            ! stage
            call limited_gradient_dx_vectorized(domain%U(:,j,STG), domain%U(:,j-1,STG), domain%U(:,j+1,STG), &
                theta_wd_NS, dstage_NS, nx)
            ! depth
            call limited_gradient_dx_vectorized(domain%depth(:,j), domain%depth(:,j-1), domain%depth(:,j+1), &
                theta_wd_NS, ddepth_NS, nx)
            ! u velocity
            !theta_wd_NS = domain%theta ! DEBUG
            call limited_gradient_dx_vectorized(domain%velocity(:,j,UH), domain%velocity(:,j-1, UH), &
                domain%velocity(:,j+1, UH), theta_wd_NS, du_NS, nx)
            ! v velocity
            !theta_wd_NS = domain%theta ! DEBUG
            call limited_gradient_dx_vectorized(domain%velocity(:,j,VH), domain%velocity(:,j-1, VH), &
                domain%velocity(:,j+1, VH), theta_wd_NS, dv_NS, nx)
            if(extrapolate_uh_vh) then
                ! uh
                call limited_gradient_dx_vectorized(domain%U(:,j,UH), domain%U(:,j-1, UH), &
                    domain%U(:,j+1, UH), theta_wd_NS, duh_NS, nx)
                ! vh
                call limited_gradient_dx_vectorized(domain%U(:,j,VH), domain%U(:,j-1, VH), &
                    domain%U(:,j+1, VH), theta_wd_NS, dvh_NS, nx)
            end if

        else
            ! Border case -- all gradients = zero
            theta_wd_NS = ZERO_dp
            dstage_NS = ZERO_dp
            ddepth_NS = ZERO_dp
            du_NS = ZERO_dp
            dv_NS = ZERO_dp
            if(extrapolate_uh_vh) then
                duh_NS = ZERO_dp
                dvh_NS = ZERO_dp
            end if
        end if

    end subroutine

    !
    ! Get the EW gradients for stage, depth, u-vel, v-vel, at row j
    !
    subroutine get_EW_limited_gradient_dx(domain, j, nx, ny, &
            theta_wd_EW, dstage_EW, ddepth_EW, du_EW, dv_EW, duh_EW, dvh_EW, extrapolate_uh_vh)

        type(domain_type), intent(in) :: domain
        integer(ip), intent(in) :: j, nx, ny
        real(dp), intent(inout) :: theta_wd_EW(nx), dstage_EW(nx), ddepth_EW(nx), du_EW(nx), dv_EW(nx), duh_EW(nx), dvh_EW(nx)
        logical, intent(in) :: extrapolate_uh_vh

        real(dp) :: mindep, maxdep, theta_local
        integer(ip) :: i

        ! limiter coefficient
        !$OMP SIMD PRIVATE(mindep, maxdep, theta_local)
        do i = 2, nx-1
            mindep = min(domain%depth(i-1, j), domain%depth(i,j), domain%depth(i+1,j)) - minimum_allowed_depth
            maxdep = max(domain%depth(i-1, j), domain%depth(i,j), domain%depth(i+1,j)) + &
                limiter_coef3*minimum_allowed_depth
            theta_local = limiter_coef4 * (mindep/maxdep - limiter_coef1)
            theta_wd_EW(i) = max(domain%theta * min(ONE_dp, theta_local), ZERO_dp)
        end do

        ! stage
        call limited_gradient_dx_vectorized(domain%U(2:(nx-1),j,STG), domain%U(1:(nx-2),j,STG), domain%U(3:nx,j,STG), &
            theta_wd_EW(2:(nx-1)), dstage_EW(2:(nx-1)), nx)
        dstage_EW(1) = ZERO_dp
        dstage_EW(nx) = ZERO_dp

        ! depth
        call limited_gradient_dx_vectorized(domain%depth(2:(nx-1),j), domain%depth(1:(nx-2),j), domain%depth(3:nx,j), &
            theta_wd_EW(2:(nx-1)), ddepth_EW(2:(nx-1)), nx)
        ddepth_EW(1) = ZERO_dp
        ddepth_EW(nx) = ZERO_dp

        ! u velocity
        !theta_wd_EW = domain%theta ! DEBUG
        call limited_gradient_dx_vectorized(domain%velocity(2:(nx-1),j,UH), domain%velocity(1:(nx-2),j,UH), &
            domain%velocity(3:nx,j, UH), theta_wd_EW(2:(nx-1)), du_EW(2:(nx-1)), nx)
        du_EW(1) = ZERO_dp
        du_EW(nx) = ZERO_dp

        ! v velocity
        !theta_wd_EW = domain%theta  ! DEBUG
        call limited_gradient_dx_vectorized(domain%velocity(2:(nx-1),j,VH), domain%velocity(1:(nx-2),j,VH), &
            domain%velocity(3:nx,j,VH), theta_wd_EW(2:(nx-1)), dv_EW(2:(nx-1)), nx)
        dv_EW(1) = ZERO_dp
        dv_EW(nx) = ZERO_dp

        ! Optionally get uh/vh
        if(extrapolate_uh_vh) then
            ! uh
            call limited_gradient_dx_vectorized(domain%U(2:(nx-1),j,UH), domain%U(1:(nx-2),j,UH), &
                domain%U(3:nx,j, UH), theta_wd_EW(2:(nx-1)), duh_EW(2:(nx-1)), nx)
            duh_EW(1) = ZERO_dp
            duh_EW(nx) = ZERO_dp

            ! vh
            call limited_gradient_dx_vectorized(domain%U(2:(nx-1),j,VH), domain%U(1:(nx-2),j,VH), &
                domain%U(3:nx,j, VH), theta_wd_EW(2:(nx-1)), dvh_EW(2:(nx-1)), nx)
            dvh_EW(1) = ZERO_dp
            dvh_EW(nx) = ZERO_dp
        end if

    end subroutine

    subroutine get_eddy_viscosity(domain, j, eddy_visc)
        !! Compute the eddy viscosity for row j. 
        !!
        !! The supported models are:
        !! - Constant [ if (domain%eddy_visc_constants(2) == 0.0)]
        !! - Constant + Smagorinsky,  with coefficients determined by domain%eddy_visc_constants(1:2) (
        !!   where the first entry is the constant coefficient and the second is the smagorinsky coefficient).
        !! - Constant + "Shiono & Knight 1991" model (again domain%eddy_visc_constants(1) is constant model, and 
        !!   domain%eddy_visc_constants(2) is the constant in "constant * depth * shear_velocity")
        !!
        !! The viscosity is limited to min(visc, dx*dx/(2*dt)) so as to not interfere with the explicit timestepping.
        !! This is done in SWASH and other models.
        !!
        ! @param domain
        ! @param j integer -- get viscosity for domain%U(:,j,:)
        ! @param eddy_vis -- vector eddy viscosity
        !
        type(domain_type), intent(in) :: domain
            !! 
        integer(ip), intent(in) :: j
            !! Compute viscosity for the row associated with domain%U(:,j,:) and y-coordinate domain%y(j)
        real(dp), intent(out) :: eddy_visc(domain%nx(1))
            !! On output holds the eddy viscosity for row j.

        integer(ip) :: i
        real(dp) :: two_dx, two_dy
        real(dp), parameter :: eps_zerodiv = 1.0e-10_dp

        character(len=charlen), parameter :: eddy_viscosity_type = 'shionoknight' !'smagorinsky'

        if(domain%eddy_visc_constants(2) == 0.0_dp) then

            !
            ! Constant eddy-viscosity
            !

            eddy_visc = domain%eddy_visc_constants(1) 

            ! Clip for stability of explicit time-stepping (see SWASH, Zijelma et al (2011), and other codes)
            two_dx = domain%distance_bottom_edge(j) + domain%distance_bottom_edge(j+1)
            do i = 2, (domain%nx(1)-1)
                two_dy = domain%distance_left_edge(i) + domain%distance_left_edge(i+1)
                eddy_visc(i) = min(eddy_visc(i), &
                    0.25_dp * min(two_dx * two_dx, two_dy * two_dy) / (2.0_dp * domain%dt_last_update + eps_zerodiv))
            end do
            ! Naive extrapolation at the edges
            eddy_visc(1) = eddy_visc(2)
            eddy_visc(domain%nx(1)) = eddy_visc(domain%nx(1)-1)

        else if(eddy_viscosity_type == 'smagorinsky') then

            !
            ! Smagorinsky eddy viscosity
            !

            if(j /= 1 .and. j /= domain%nx(2)) then
                ! Typical case

                two_dx = domain%distance_bottom_edge(j) + domain%distance_bottom_edge(j+1)

                do i = 2, (domain%nx(1)-1)

                    two_dy = domain%distance_left_edge(i) + domain%distance_left_edge(i+1)

                    ! Smagorinsky type formula with extra constant
                    eddy_visc(i) = domain%eddy_visc_constants(1) + &
                        domain%eddy_visc_constants(2) * domain%area_cell_y(j) * sqrt( &
                        ( (domain%velocity(i+1,j,UH) - domain%velocity(i-1,j,UH))/two_dx)**2 + &
                        ( (domain%velocity(i,j+1,VH) - domain%velocity(i,j-1,VH))/two_dy)**2 + &
                        HALF_dp * ( &
                            abs(domain%velocity(i,j+1,UH) - domain%velocity(i,j-1,UH))/two_dy + &
                            abs(domain%velocity(i+1,j,VH) - domain%velocity(i-1,j,VH))/two_dx &
                            )**2)

                     ! Clip for stability of explicit time-stepping (see SWASH, Zijelma et al (2011), and other codes)
                     eddy_visc(i) = min(eddy_visc(i), &
                         0.25_dp * min(two_dx * two_dx, two_dy * two_dy) / (2.0_dp * domain%dt_last_update + eps_zerodiv))
                end do

            else if(j == 1) then
                ! South boundary case -- avoid j-1 indexing

                two_dx = domain%distance_bottom_edge(j) + domain%distance_bottom_edge(j+1)

                do i = 2, (domain%nx(1)-1)

                    two_dy = domain%distance_left_edge(i) + domain%distance_left_edge(i+1)

                    ! Smagorinsky type formula with extra constant
                    eddy_visc(i) = domain%eddy_visc_constants(1) + &
                        domain%eddy_visc_constants(2) * domain%area_cell_y(j) * sqrt( &
                        ( (domain%velocity(i+1,j,UH) - domain%velocity(i-1,j,UH))/two_dx)**2 + &
                        ( 2*(domain%velocity(i,j+1,VH) - domain%velocity(i,j,VH))/two_dy)**2 + &
                        HALF_dp * ( &
                            2*abs(domain%velocity(i,j+1,UH) - domain%velocity(i,j,UH))/two_dy + &
                            abs(domain%velocity(i+1,j,VH) - domain%velocity(i-1,j,VH))/two_dx )**2 )

                     ! Clip for stability of explicit time-stepping (see SWASH, Zijelma et al (2011), and other codes)
                     eddy_visc(i) = min(eddy_visc(i), &
                         0.25_dp * min(two_dx * two_dx, two_dy * two_dy) / (2.0_dp * domain%dt_last_update + eps_zerodiv))
                end do

            else if(j == domain%nx(2)) then
                ! North boundary case -- avoid j+1 indexing

                two_dx = domain%distance_bottom_edge(j-1) + domain%distance_bottom_edge(j)
                do i = 2, (domain%nx(1)-1)

                    two_dy = domain%distance_left_edge(i) + domain%distance_left_edge(i+1)

                    ! Smagorinsky type formula with extra constant
                    eddy_visc(i) = domain%eddy_visc_constants(1) + &
                        domain%eddy_visc_constants(2) * domain%area_cell_y(j) * sqrt( &
                        ( (domain%velocity(i+1,j,UH) - domain%velocity(i-1,j,UH))/two_dx)**2 + &
                        ( 2*(domain%velocity(i,j,VH) - domain%velocity(i,j-1,VH))/two_dy)**2 + &
                        HALF_dp * ( &
                            2*abs(domain%velocity(i,j,UH) - domain%velocity(i,j-1,UH))/two_dy + &
                            abs(domain%velocity(i+1,j,VH) - domain%velocity(i-1,j,VH))/two_dx )**2 )

                     ! Clip for stability of explicit time-stepping (see SWASH, Zijelma et al (2011), and other codes)
                     eddy_visc(i) = min(eddy_visc(i), &
                         0.25_dp * min(two_dx * two_dx, two_dy * two_dy) / (2.0_dp * domain%dt_last_update + eps_zerodiv))
                end do

            end if

            ! Naive extrapolation at the edges
            eddy_visc(1) = eddy_visc(2)
            eddy_visc(domain%nx(1)) = eddy_visc(domain%nx(1)-1)

        else if(eddy_viscosity_type == 'shionoknight') then
            
            !
            ! A common eddy-viscosity, used in the Compound-Channel-Flow model of Shiono and Knight (1991) and many other
            ! places.
            !     eddy_viscocity = constant * depth * shear_velocity
            ! where
            !     shear_velocity = sqrt(bed_shear_stress/rho) = sqrt( d ^(-1/3) * n^2 * g * speed^2 )
            !
            eddy_visc = domain%eddy_visc_constants(1) + &
                domain%eddy_visc_constants(2) * domain%depth(:,j)**(1.0_dp - 1.0_dp/6.0_dp) * &
                sqrt(domain%manning_squared(:,j) * gravity * (domain%velocity(:,j,UH)**2 + domain%velocity(:,j,VH)**2))

            ! Clip for stability of explicit time-stepping (see SWASH, Zijelma et al (2011), and other codes)
            two_dx = domain%distance_bottom_edge(j) + domain%distance_bottom_edge(j+1)
            do i = 2, (domain%nx(1)-1)
                two_dy = domain%distance_left_edge(i) + domain%distance_left_edge(i+1)
                eddy_visc(i) = min(eddy_visc(i), &
                    0.25_dp * min(two_dx * two_dx, two_dy * two_dy) / (2.0_dp * domain%dt_last_update + eps_zerodiv))
            end do

            ! Naive extrapolation at the edges
            eddy_visc(1) = eddy_visc(2)
            eddy_visc(domain%nx(1)) = eddy_visc(domain%nx(1)-1)

        end if

    end subroutine

!end subroutine
