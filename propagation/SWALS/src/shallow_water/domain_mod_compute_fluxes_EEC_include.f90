! Compute fluxes for 2D shallow water equations on structured grid
!
! Use an EEC method
!
! Updated variables:
!    domain%flux_NS, domain%flux_EW, domain%explicit_source, domain%explicit_source_VH_j_minus_1, max_dt_out
!
! @param domain the model domain type for which fluxes etc will be computed
! @param max_dt_out optional real scalar, if provided is set to the max_dt
!     allowed based on CFL limit.
!
subroutine compute_fluxes_EEC(domain, max_dt_out)

    type(domain_type), intent(inout) :: domain
    real(dp), intent(inout) :: max_dt_out

    ! wavespeeds
    real(dp):: s_max, s_min, gs_pos, gs_neg, sminsmax
    ! stage/depth at + and - side of edge
    real(dp):: stage_pos, stage_neg, stage_pos_star, stage_neg_star, &
               depth_pos, depth_pos_star, depth_pos_c,&
               depth_neg, depth_neg_star, depth_neg_c, depthsq_neg, &
               z_half, z_neg, z_pos, dz
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
    ! Option to use upwind flux for transverse momentum component, as in
    ! e.g. "Chertock et al. (2017) Well-balanced schemes for
    ! the shallow water equations with Coriolis forces"
    logical, parameter :: upwind_transverse_momentum = .false.
    logical, parameter :: upwind_normal_momentum = .false.
    logical, parameter :: reduced_momentum_diffusion = .true. !.false.

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

    real(dp) :: bed_j_minus_1(domain%nx(1))

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
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, half_cfl, nx, ny, loop_work_count) REDUCTION(MAX:max_dt_inv)
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

        ! Assume distance_left_edge is constant but distance_bottom_edge
        ! might change with y (spherical coordinates). For cartesian coordinates
        ! we could move this outside the loop.
        dx_cfl_inv(1) = ONE_dp/(domain%distance_bottom_edge(j) * half_cfl)
        dx_cfl_inv(2) = ONE_dp/(domain%distance_left_edge(1) * half_cfl)

        ! Get bed at j-1 (might improve cache access)

        ! Get bottom edge values
        !call get_bottom_edge_values(domain, j, nx, ny, &
        !    theta_wd_neg_B, theta_wd_pos_B, &
        !    stage_neg_B, stage_pos_B, &
        !    depth_neg_B, depth_pos_B, &
        !    u_neg_B, u_pos_B, &
        !    v_neg_B, v_pos_B, &
        !    reuse_gradients = (jlast == j-1))
        ! By setting jlast, we allow the above subroutine call to reuse
        ! limited gradient values from the previous loop iterate
        jlast = j

        ! Get left edge values
        !call get_left_edge_values(domain, j, nx, &
        !    theta_wd_neg_L, theta_wd_pos_L, &
        !    stage_neg_L, stage_pos_L, &
        !    depth_neg_L, depth_pos_L, &
        !    u_neg_L, u_pos_L, &
        !    v_neg_L, v_pos_L)

        ! Maybe better access in loop
        stage_neg_B = domain%U(:,j-1,STG)
        u_neg_B = domain%velocity(:,j-1,UH)
        v_neg_B = domain%velocity(:,j-1,VH)
        bed_j_minus_1 = domain%U(:,j-1,ELV)

        ! NOTE: much of the loop below can be vectorized, but care will be needed
        ! for the max_dt_inv term
        do i = 2, nx
            !
            ! North-South flux computation
            ! Compute flux_NS(i,j) = flux(i,j-1/2)
            !
            ! Cell i,j has North-South flux(i,j+1/2) at flux_NS index i,j+1,
            !          and North-South flux(i,j-1/2) at flux_NS index i,j
            !

            ! Negative/positive bottom edge values
            stage_neg = HALF_dp * (domain%U(i,j,STG) + stage_neg_B(i))
            z_neg = HALF_dp * (domain%U(i,j,ELV) + bed_j_minus_1(i))
            depth_neg = max(stage_neg - z_neg, ZERO_dp)
            depthsq_neg = HALF_dp * ( (domain%U(i,j,STG) - domain%U(i,j,ELV))**2 + (stage_neg_B(i) - bed_j_minus_1(i))**2)
            u_neg = HALF_dp * (domain%velocity(i,j,UH) + u_neg_B(i))
            v_neg = HALF_dp * (domain%velocity(i,j,VH) + v_neg_B(i))
            dz = domain%U(i,j,ELV) - bed_j_minus_1(i)

            if(depth_neg == ZERO_dp .or. (domain%U(i,j,STG) < bed_j_minus_1(i)) .or. (stage_neg_B(i) < domain%U(i,j,ELV))) then
                v_neg = ZERO_dp
                u_neg = ZERO_dp
                !depthsq_neg = ZERO_dp
                !depth_neg = ZERO_dp
                !z_neg = min(domain%U(i,j,STG), stage_neg_B(i))
            end if

            domain%flux_NS(i,j,STG) = v_neg * depth_neg * domain%distance_bottom_edge(j)
            domain%flux_NS(i,j,UH) = u_neg * v_neg * depth_neg * domain%distance_bottom_edge(j)
            domain%flux_NS(i,j,VH) = (v_neg * v_neg * depth_neg + half_gravity*depthsq_neg) * domain%distance_bottom_edge(j)

            ! Bed-slope
            domain%explicit_source_VH_j_minus_1(i,j) = -&
                domain%area_cell_y(j-1) * half_gravity * depth_neg * dz/domain%distance_left_edge(1)
            domain%explicit_source(i,j,VH) = domain%explicit_source(i,j,VH) - &
                domain%area_cell_y(j)   * half_gravity * depth_neg * dz/domain%distance_left_edge(1)


            ! Timestep
            max_speed = (max(abs(u_neg), abs(v_neg)) + sqrt(depth_neg * gravity))
            if (max_speed > EPS) then
                max_dt_inv = max(max_dt_inv, max_speed * dx_cfl_inv(2))
            end if

            !
            ! EW Flux computation
            !
            stage_neg = HALF_dp * (domain%U(i-1,j,STG) + domain%U(i,j,STG))
            z_neg = HALF_dp * (domain%U(i-1,j,ELV) + domain%U(i,j,ELV))
            depth_neg = max(stage_neg - z_neg, ZERO_dp)
            depthsq_neg = HALF_dp * ( (domain%U(i-1,j,STG) - domain%U(i-1,j,ELV))**2 + (domain%U(i,j,STG) - domain%U(i,j,ELV))**2)
            u_neg = HALF_dp * (domain%velocity(i-1,j,UH) + domain%velocity(i,j,UH))
            v_neg = HALF_dp * (domain%velocity(i-1,j,VH) + domain%velocity(i,j,VH))
            dz = domain%U(i,j,ELV) - domain%U(i-1,j,ELV)

            if(depth_neg == ZERO_dp .or. (domain%U(i,j,STG) < domain%U(i-1,j,ELV)) .or. &
               (domain%U(i-1,j,STG) < domain%U(i,j,ELV))) then
                v_neg = ZERO_dp
                u_neg = ZERO_dp
                !depthsq_neg = ZERO_dp
                !depth_neg = ZERO_dp
            end if

            domain%flux_EW(i,j,STG) = u_neg * depth_neg * domain%distance_left_edge(1)
            domain%flux_EW(i,j,UH) = (u_neg * u_neg * depth_neg + half_gravity*depthsq_neg) * domain%distance_left_edge(1)
            domain%flux_EW(i,j,VH) = (v_neg * u_neg * depth_neg ) * domain%distance_left_edge(1)

            ! Bed-slope
            dxb = (0.5_dp * (domain%distance_bottom_edge(j+1) + domain%distance_bottom_edge(j)))
            domain%explicit_source(i-1,j,UH) = domain%explicit_source(i-1,j,UH) - &
                domain%area_cell_y(j) * half_gravity * depth_neg * dz/dxb
            domain%explicit_source(i,j,UH) = domain%explicit_source(i,j,UH) - &
                domain%area_cell_y(j) * half_gravity * depth_neg * dz/dxb

            max_speed = max(abs(u_neg), abs(v_neg)) + sqrt(depth_neg * gravity)
            if (max_speed > EPS) then
                max_dt_inv = max(max_dt_inv, max_speed * dx_cfl_inv(2))
            end if

#ifdef SPHERICAL

            stop 'Not yet supporting spherical'
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

        end do
    end do
    !$OMP END PARALLEL

    ! Periodic BC
    !domain%flux_NS(2:nx,1, :) = domain%flux_NS(2:nx, ny, :)
    !domain%flux_NS(2:nx,ny+1, :) = domain%flux_NS(2:nx,2, :)
    !domain%flux_EW(1,:, :) = domain%flux_EW(nx, :, :)
    !domain%flux_EW(nx+1,:, :) = domain%flux_EW(2, :,:)

    max_dt = ONE_dp/max_dt_inv
    max_dt_out = max_dt

end subroutine
