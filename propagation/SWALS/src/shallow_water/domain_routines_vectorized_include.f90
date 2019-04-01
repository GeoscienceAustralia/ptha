! Include this in "domain_mod.f90"

! This file includes vectorized versions of some key routines
! for the domain. Note they are not always faster -- initial experience
! was that compute_fluxes_vectorized was slower, whereas 
! 'extrapolate_edge_second_order_vectorized' is much faster. However, it may
! depend on the machine and compiler.

    pure subroutine flux_NS_vectorize(stage_pos, stage_neg, depth_pos, depth_neg, &
        u_pos, u_neg, v_pos, v_neg, depth_pos_c, depth_neg_c, bed_j_minus_1_v, bed_j_v, &
        flux1, flux2, flux3, half_g_hh_edge, s_max, s_min, &
        bed_slope_pressure_n, bed_slope_pressure_s, &
        distance_bottom_edge, n)
        
        integer(ip), intent(in) :: n
        real(dp), intent(in) :: stage_pos(n), stage_neg(n), depth_pos(n), depth_neg(n), &
            depth_pos_c(n), depth_neg_c(n), distance_bottom_edge, &
            bed_j_minus_1_v(n), bed_j_v(n)
        real(dp), intent(inout) :: v_pos(n), v_neg(n), u_pos(n), u_neg(n)
        real(dp), intent(out) :: flux1(n), flux2(n), flux3(n), half_g_hh_edge(n), s_max(n), s_min(n)
        real(dp), intent(out) :: bed_slope_pressure_n(n), bed_slope_pressure_s(n)

        integer(ip), parameter :: n1 = 2 !vectorization_size
        real(dp) :: z_neg(n1), z_pos(n1), z_half(n1), stage_pos_star(n1), stage_neg_star(n1), &
            depth_pos_star(n1), depth_neg_star(n1), gs_pos(n1), gs_neg(n1), &
            ud_pos(n1), ud_neg(n1), vd_pos(n1), vd_neg(n1), sminsmax(n1), vel_beta_neg(n1), vel_beta_pos(n1), &
            denom(n1), inv_denom(n1)

        ! FIXME: Copied from compute_fluxes
        real(dp), parameter :: EPS = 1.0e-06_dp, diffusion_scale = ONE_dp
        real(dp), parameter:: half_gravity = HALF_dp * gravity

        ! Bed elevation
        z_neg = stage_neg - depth_neg
        z_pos = stage_pos - depth_pos
        z_half = max(z_neg, z_pos)

        stage_neg_star = max(stage_neg, z_half)
        depth_neg_star = stage_neg_star - z_half

        ! Velocity (in NS direction)
        where(depth_neg_star == ZERO_dp)
            v_neg = ZERO_dp
            u_neg = ZERO_dp
        end where
        ! Audusse type depth_integrated_velocity correction
        vd_neg = v_neg * depth_neg_star
        ud_neg = u_neg * depth_neg_star

        ! Gravity wave celerity
        gs_neg = sqrt(gravity * depth_neg_star)
        
        stage_pos_star = max(stage_pos, z_half)
        depth_pos_star = stage_pos_star - z_half

        ! Velocity (in NS direction)
        where(depth_pos_star == ZERO_dp)
            v_pos = ZERO_dp
            u_pos = ZERO_dp
        end where

        ! Correct depth-integrated_velocity (Audusse type approach)
        vd_pos = v_pos * depth_pos_star
        ud_pos = u_pos * depth_pos_star

        ! Gravity wave celerity
        gs_pos = sqrt(gravity * depth_pos_star)

        ! Wave-celerities
        s_max = max(max(v_neg + gs_neg, v_pos + gs_pos), ZERO_dp)
        s_min = min(min(v_neg - gs_neg, v_pos - gs_pos), ZERO_dp)

        denom = s_max - s_min

        where(denom > EPS)
            inv_denom = distance_bottom_edge / denom
        elsewhere
            inv_denom = ZERO_dp
        end where

        sminsmax = s_min * s_max * diffusion_scale
        vel_beta_neg = v_neg * advection_beta
        vel_beta_pos = v_pos * advection_beta

        ! Advection flux terms 
        flux1 = &
            (s_max * vd_neg - &
             s_min * vd_pos + &
            sminsmax * (stage_pos_star - stage_neg_star)) * inv_denom
        flux2 = &
            (s_max * ud_neg * vel_beta_neg - &
             s_min * ud_pos * vel_beta_pos + &
            sminsmax * (ud_pos - ud_neg)) * inv_denom
        flux3 = &
            (s_max * (vd_neg * vel_beta_neg) - & 
             s_min * (vd_pos * vel_beta_pos) + &
            sminsmax * (vd_pos - vd_neg)) * inv_denom

        ! Bring the NS gh^2/2 term out of the flux_NS. Note the term is already
        ! multiplied by edge-length. Later we put it in the 'explicit source' terms.
        half_g_hh_edge =  inv_denom * half_gravity * (&
            s_max * depth_neg_star * depth_neg_star - &
            s_min * depth_pos_star * depth_pos_star)
        
        !! Pressure gradient term at i,j-1, north side
        bed_slope_pressure_n = half_gravity * ( &
            (depth_neg_star - depth_neg)*(depth_neg_star + depth_neg) - &
            (depth_neg + depth_neg_c)*(z_neg- bed_j_minus_1_v))
       

        !! Pressure gradient term at i,j, south-side
        bed_slope_pressure_s = half_gravity * ( & 
            (depth_pos_star - depth_pos) * (depth_pos_star + depth_pos) - &
            (depth_pos + depth_pos_c)*(z_pos - bed_j_v) ) 

    end subroutine 

    pure subroutine flux_EW_vectorize(stage_pos, stage_neg, depth_pos, depth_neg, &
        u_pos, u_neg, v_pos, v_neg, depth_pos_c, depth_neg_c, bed_i_minus_1_v, bed_j_v, &
        flux1, flux2, flux3, half_g_hh_edge, s_max, s_min, &
        bed_slope_pressure_e, bed_slope_pressure_w, &
        dsl, n)
        
        integer(ip), intent(in) :: n
        real(dp), intent(in) :: stage_pos(n), stage_neg(n), depth_pos(n), depth_neg(n), &
            depth_pos_c(n), depth_neg_c(n), dsl(n), &
            bed_i_minus_1_v(n), bed_j_v(n)
        real(dp), intent(inout) :: v_pos(n), v_neg(n), u_pos(n), u_neg(n)
        real(dp), intent(out) :: flux1(n), flux2(n), flux3(n), half_g_hh_edge(n), s_max(n), s_min(n)
        real(dp), intent(out) :: bed_slope_pressure_e(n), bed_slope_pressure_w(n)

        integer(ip), parameter :: n1 = 2 !vectorization_size
        real(dp) :: z_neg(n1), z_pos(n1), z_half(n1), stage_pos_star(n1), stage_neg_star(n1), &
            depth_pos_star(n1), depth_neg_star(n1), gs_pos(n1), gs_neg(n1), &
            ud_pos(n1), ud_neg(n1), vd_pos(n1), vd_neg(n1), sminsmax(n1), vel_beta_neg(n1), vel_beta_pos(n1), &
            denom(n1), inv_denom(n1)

        ! FIXME: Copied from compute_fluxes
        real(dp), parameter :: EPS = 1.0e-06_dp, diffusion_scale = ONE_dp
        real(dp), parameter:: half_gravity = HALF_dp * gravity

        ! Bed elevation
        z_neg = stage_neg - depth_neg
        z_pos = stage_pos - depth_pos
        z_half = max(z_neg, z_pos)

        stage_neg_star = max(stage_neg, z_half)
        depth_neg_star = stage_neg_star - z_half

        ! Velocity (in EW direction)
        where(depth_neg_star == ZERO_dp)
            u_neg = ZERO_dp
            v_neg = ZERO_dp
        end where
        ! Audusse type depth-integrated-velocity correction
        ud_neg = u_neg * depth_neg_star
        vd_neg = v_neg * depth_neg_star

        gs_neg = sqrt(gravity * depth_neg_star)

        stage_pos_star = max(stage_pos, z_half)
        depth_pos_star = stage_pos_star - z_half
        ! Velocity (in NS direction)
        where(depth_pos_star == ZERO_dp)
            u_pos = ZERO_dp
            v_pos = ZERO_dp
        end where

        ud_pos = u_pos * depth_pos_star
        vd_pos = v_pos * depth_pos_star

        gs_pos = sqrt(gravity * depth_pos_star)

        ! Wave-celerities
        s_max = max(max(u_neg + gs_neg, u_pos + gs_pos), ZERO_dp)
        s_min = min(min(u_neg - gs_neg, u_pos - gs_pos), ZERO_dp)

        denom = s_max - s_min

        where(denom > EPS)
            inv_denom = dsl / denom
        elsewhere
            inv_denom = ZERO_dp
        end where

        sminsmax = s_min * s_max * diffusion_scale
        vel_beta_neg = u_neg * advection_beta
        vel_beta_pos = u_pos * advection_beta

        flux1 = &
            (s_max * ud_neg - &
             s_min * ud_pos + &
            sminsmax * (stage_pos_star - stage_neg_star)) * inv_denom
        flux2 = &
            (s_max * (ud_neg * vel_beta_neg ) - &
             s_min * (ud_pos * vel_beta_pos ) + &
             sminsmax * (ud_pos - ud_neg)) * inv_denom
        flux3 = &
            (s_max * vd_neg * vel_beta_neg - &
             s_min * vd_pos * vel_beta_pos + &
            sminsmax * (vd_pos - vd_neg)) * inv_denom

        ! Bring EW pressure gradient term out of flux_EW. Note it is already
        ! multipied by edge-length
        half_g_hh_edge = inv_denom * half_gravity * (&
            s_max * depth_neg_star * depth_neg_star - &
            s_min * depth_pos_star * depth_pos_star)
    
        !! Pressure gradient term at i-1,j -- east side of cell i-1,j
        bed_slope_pressure_e = half_gravity * (&
            (depth_neg_star - depth_neg)*(depth_neg_star + depth_neg) - &
            (depth_neg + depth_neg_c)*(z_neg - bed_i_minus_1_v))

        !! Pressure gradient term at i,j -- west_side

        bed_slope_pressure_w = half_gravity * (&
            (depth_pos_star- depth_pos) * (depth_pos_star + depth_pos) - &
            (depth_pos + depth_pos_c)*(z_pos - bed_j_v) )


    end subroutine

    !
    ! Compute the fluxes, and other things, in preparation for an update of U 
    !
    ! This is a vectorized version of compute_fluxes in "domain_mod.f90". Initial testing
    ! suggested it was slower than the original.
    !
    ! Update values of:
    ! domain%flux_NS, domain%flux_EW, domain%max_dt, domain%explicit_source, 
    ! domain%explicit_source_VH_j_minus_1, domain%boundary_flux_store
    !
    ! @param domain the model domain type for which fluxes etc will be computed
    ! @param max_dt_out optional real scalar, if provided is set to the max_dt
    !     allowed based on CFL limit.
    subroutine compute_fluxes_vectorized(domain, max_dt_out)
        ! Compute fluxes for 2D shallow water equations on structured grid
        ! Use an Audusse type method, like ANUGA, but structured

        class(domain_type), intent(inout):: domain
        real(dp), optional, intent(inout) :: max_dt_out

        integer(ip):: i, j, nx, ny, jlast, ii
        integer(ip), parameter :: v = 2 !vectorization_size ! Size used to try to vectorize the routine
        ! wavespeeds
        real(dp):: s_max(v), s_min(v), gs_pos(v), gs_neg(v), sminsmax(v)
        ! stage/depth at + and - side of edge
        real(dp):: stage_pos(v), stage_neg(v), &
                   depth_pos(v), depth_pos_c(v),&
                   depth_neg(v), depth_neg_c(v),&
                   z_half(v), z_neg(v), z_pos(v)
        ! velocities and momenta at + and - sides of edge
        real(dp):: u_pos(v), v_pos(v), u_neg(v), v_neg(v), &
            ud_neg(v), vd_neg(v), ud_pos(v), vd_pos(v), &
            vel_beta_neg(v), vel_beta_pos(v)
        ! convenience variables
        real(dp):: max_speed(v)
        real(dp):: bed_slope_pressure_s(v), bed_slope_pressure_n(v), bed_slope_pressure_e(v), &
            bed_slope_pressure_w(v), half_g_hh_edge(v), flux1(v), flux2(v), flux3(v), &
            bed_j_v(v), bed_j_minus_1_v(v), bed_i_minus_1_v(v), dsl(v)
        real(dp):: half_cfl, max_dt_inv, max_dt, dx_cfl_inv(2)
        character(len=charlen):: timer_name
        real(dp), parameter :: EPS = 1.0e-06_dp, diffusion_scale = ONE_dp
        real(dp), parameter:: half_gravity = HALF_dp * gravity
        integer(ip):: n_ext, loop_work_count, j_low, j_high, my_omp_id, n_omp_threads
        integer(ip) :: imn, imx, vsize
        
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

        real(dp) :: bed_j_minus_1(domain%nx(1)), bed_j(domain%nx(1))

        half_cfl = HALF_dp * domain%cfl
        ny = domain%nx(2)
        nx = domain%nx(1)


        ! Set dt to a high value (it will drop)
        max_dt = domain%maximum_timestep
        ! By computing the inverse we avoid division in the loop
        max_dt_inv = ONE_dp/max_dt

        ! Need to have depth/velocity up-to-date for flux computation
        call domain%compute_depth_and_velocity()


        ! NOTE: Here we manually determine the openmp loop indices, because
        ! some efficiency can be gained by moving through a contiguous chunk
        !
        loop_work_count = ny - 2 + 1 ! Number of indices between 'ny' and '2'
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
        !!!$OMP DO SCHEDULE(GUIDED)
        !!do j = 2, ny
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

            bed_j_minus_1 = domain%U(:,j-1,ELV)
            bed_j = domain%U(:,j, ELV)

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

            !!!!!!$OMP SIMD ! Might not speed things up?
            !do i = 2, nx
            ! Here we try to loop so that the vectors have length 32 (could
            ! also use 64 or something else). Alternatively
            ! one can remove this do loop, with imn = 2, imx = domain%nx(1)
            ! In either case, the non-vectorized loop seems faster
            do i = 2, domain%nx(1), v

                ! Indices for packing vectors (which are mostly of size v)
                imn = i
                imx = min(imn + v - 1, domain%nx(1))
                vsize = min(v, imx - imn + 1)

                bed_j_v(1:vsize) = bed_j(imn:imx)
                bed_i_minus_1_v(1:vsize) = bed_j((imn-1):(imx-1))
                bed_j_minus_1_v(1:vsize) = bed_j_minus_1(imn:imx)
                dsl(1:vsize) = domain%distance_left_edge(imn:imx)
                !
                ! North-South flux computation
                ! Compute flux_NS(i,j) = flux(i,j-1/2)
                !
                ! Cell i,j has North-South flux(i,j+1/2) at flux_NS index i,j+1, 
                !          and North-South flux(i,j-1/2) at flux_NS index i,j
                !

                ! Negative/positive bottom edge values
                stage_neg(1:vsize) = stage_neg_B(imn:imx)
                stage_pos(1:vsize) = stage_pos_B(imn:imx)
                depth_neg(1:vsize) = depth_neg_B(imn:imx)
                depth_pos(1:vsize) = depth_pos_B(imn:imx)
                u_neg(1:vsize) = u_neg_B(imn:imx)
                u_pos(1:vsize) = u_pos_B(imn:imx)
                v_neg(1:vsize) = v_neg_B(imn:imx)
                v_pos(1:vsize) = v_pos_B(imn:imx)
                depth_neg_c(1:vsize) = domain%depth(imn:imx, j-1)
                depth_pos_c(1:vsize) = domain%depth(imn:imx, j)

                ! Do the flux computation
                call flux_NS_vectorize(stage_pos, stage_neg, depth_pos, depth_neg, &
                    u_pos, u_neg, v_pos, v_neg, depth_pos_c, depth_neg_c, &
                    bed_j_minus_1_v, bed_j_v, &
                    flux1, flux2, flux3, half_g_hh_edge, s_max, s_min, &
                    bed_slope_pressure_n, bed_slope_pressure_s, &
                    domain%distance_bottom_edge(j), v)

                ! Pack fluxes into main flux array
                domain%flux_NS(imn:imx, j, STG) = flux1(1:vsize)
                domain%flux_NS(imn:imx, j, UH) = flux2(1:vsize)
                domain%flux_NS(imn:imx, j, VH) = flux3(1:vsize)

                !! NOTE regarding OPENMP!!
                !! We cannot naively do:
                ! domain%explicit_source(i, j - 1 ,VH) = domain%explicit_source(i, j-1, VH) + bed_slope_pressure_n + ...
                !! Since j-1 might also be updated on another OMP thread
                !! The solution is below
                !! Other pressure gradient cases refer to i,j, or i-1,j, so should be ok
                domain%explicit_source_VH_j_minus_1(imn:imx,j) = &
                    bed_slope_pressure_n(1:vsize) * domain%distance_bottom_edge(j) - half_g_hh_edge(1:vsize) ! Add to explicit_source later

                domain%explicit_source(imn:imx, j, VH) = domain%explicit_source(imn:imx, j, VH) - &
                    bed_slope_pressure_s(1:vsize) * domain%distance_bottom_edge(j) + half_g_hh_edge(1:vsize)

                ! Timestep
                max_speed = max(s_max, -s_min)
                do ii = 1, vsize
                    if (max_speed(ii) > EPS) then
                        max_dt_inv = max(max_dt_inv, max_speed(ii) * dx_cfl_inv(2))
                    end if
                end do

            end do

            do i = 2, domain%nx(1), v

                ! Indices for packing vectors (which are mostly of size v)
                imn = i
                imx = min(imn + v - 1, domain%nx(1))
                vsize = min(v, imx - imn + 1)
                !
                !
                ! EW Flux computation
                !
                !

                ! left edge variables 
                stage_neg(1:vsize) = stage_neg_L(imn:imx)
                stage_pos(1:vsize) = stage_pos_L(imn:imx)
                depth_neg(1:vsize) = depth_neg_L(imn:imx)
                depth_pos(1:vsize) = depth_pos_L(imn:imx)
                u_neg(1:vsize) = u_neg_L(imn:imx)
                u_pos(1:vsize) = u_pos_L(imn:imx)
                v_neg(1:vsize) = v_neg_L(imn:imx)
                v_pos(1:vsize) = v_pos_L(imn:imx)
                depth_neg_c(1:vsize) = domain%depth((imn-1):(imx-1), j)
                depth_pos_c(1:vsize) = domain%depth(imn:imx, j)

                bed_j_v(1:vsize) = bed_j(imn:imx)
                bed_i_minus_1_v(1:vsize) = bed_j((imn-1):(imx-1))
                bed_j_minus_1_v(1:vsize) = bed_j_minus_1(imn:imx)
                dsl(1:vsize) = domain%distance_left_edge(imn:imx)

                ! Do the flux computation
                call flux_EW_vectorize(stage_pos, stage_neg, depth_pos, depth_neg, &
                    u_pos, u_neg, v_pos, v_neg, depth_pos_c, depth_neg_c, &
                    bed_i_minus_1_v, bed_j_v, &
                    flux1, flux2, flux3, half_g_hh_edge, s_max, s_min, &
                    bed_slope_pressure_e, bed_slope_pressure_w, &
                    dsl, v)

                ! Pack fluxes into main array
                domain%flux_EW(imn:imx, j, STG) = flux1(1:vsize)
                domain%flux_EW(imn:imx, j, UH) = flux2(1:vsize)
                domain%flux_EW(imn:imx, j, VH) = flux3(1:vsize)

                domain%explicit_source((imn-1):(imx-1), j, UH) = domain%explicit_source((imn-1):(imx-1), j ,UH) + &
                    bed_slope_pressure_e(1:vsize) * dsl(1:vsize) - half_g_hh_edge(1:vsize)

                domain%explicit_source(imn:imx, j, UH) = domain%explicit_source(imn:imx, j, UH) - &
                    bed_slope_pressure_w(1:vsize) * dsl(1:vsize) + half_g_hh_edge(1:vsize)

#ifdef SPHERICAL

                print*, 'NEED TO ADD THE SPHERICAL TERMS IN THE VECTORIZED SOLVER'
                stop

                domain%explicit_source(imn:imx, j, VH) = domain%explicit_source(imn:imx, j, VH) + &
                    half_gravity * domain%depth(imn:imx, j) * domain%depth(imn:imx, j) * &
                    (domain%distance_bottom_edge(j+1) - domain%distance_bottom_edge(j))

#endif

                ! Timestep
                max_speed = max(s_max, -s_min)
                do ii = 1, vsize 
                    if (max_speed(ii) > EPS) then
                        max_dt_inv = max(max_dt_inv, max_speed(ii) * dx_cfl_inv(1))
                    end if
                end do

            end do
        end do
        !!!$OMP END DO
        !$OMP END PARALLEL
        
        max_dt = ONE_dp/max_dt_inv
        !domain%max_dt = max_dt
        if(present(max_dt_out)) max_dt_out = max_dt

        !
        ! Spatially integrate boundary fluxes. Order is N, E, S, W. 
        ! Note we integrate in the INTERIOR (i.e. cut the edge rows/columns, then
        ! integrate the boundary of what remains). 
        !
        n_ext = domain%exterior_cells_width
        ! Outward boundary flux over the north
        domain%boundary_flux_store(1) = sum(domain%flux_NS( (1+n_ext):(nx-n_ext), ny+1-n_ext, STG))
        ! Outward boundary flux over the east
        domain%boundary_flux_store(2) = sum(domain%flux_EW( nx+1-n_ext, (1+n_ext):(ny-n_ext), STG))
        ! Outward boundary flux over the south
        domain%boundary_flux_store(3) = -sum(domain%flux_NS( (1+n_ext):(nx-n_ext), 1+n_ext, STG))
        ! Outward boundary flux over the west
        domain%boundary_flux_store(4) = -sum(domain%flux_EW( 1+n_ext, (1+n_ext):(ny-n_ext), STG))

        !
        ! Compute fluxes relevant to the multidomain case. Note similar code occurs at the 
        ! end of subroutine compute_fluxes -- however, here we deal with the fact that 
        ! the fluxes are stored in domain%U
        !
        if(domain%nesting%my_index > 0) then
            ! We are in a multidomain -- carefully compute fluxes through
            ! exterior boundaries
            domain%boundary_flux_store_exterior = ZERO_dp
       
            ! Here we implement masked versions of the boundary flux sums above, only counting cells
            ! where the priority domain is receiving/sending the fluxes on actual physical boundaries 

            ! North boundary
            if(domain%boundary_exterior(1)) then
                domain%boundary_flux_store_exterior(1) = sum(&
                    domain%flux_NS( (1+n_ext):(nx-n_ext), ny+1-n_ext, STG),&
                    mask = (&
                        domain%nesting%priority_domain_index((1+n_ext):(nx-n_ext), ny-n_ext) == domain%nesting%my_index .and. &
                        domain%nesting%priority_domain_image((1+n_ext):(nx-n_ext), ny-n_ext) == domain%nesting%my_image ))
            end if

            ! East boundary
            if(domain%boundary_exterior(2)) then
                domain%boundary_flux_store_exterior(2) = sum(&
                    domain%flux_EW( nx+1-n_ext, (1+n_ext):(ny-n_ext), STG),&
                    mask = (&
                        domain%nesting%priority_domain_index(nx-n_ext, (1+n_ext):(ny-n_ext)) == domain%nesting%my_index .and. &
                        domain%nesting%priority_domain_image(nx-n_ext, (1+n_ext):(ny-n_ext)) == domain%nesting%my_image ))
            end if

            ! South boundary
            if(domain%boundary_exterior(3)) then
                domain%boundary_flux_store_exterior(3) = -sum(&
                    domain%flux_NS( (1+n_ext):(nx-n_ext), 1+n_ext, STG),&
                    mask = (&
                        domain%nesting%priority_domain_index((1+n_ext):(nx-n_ext), 1+n_ext) == domain%nesting%my_index .and. &
                        domain%nesting%priority_domain_image((1+n_ext):(nx-n_ext), 1+n_ext) == domain%nesting%my_image ))
            end if

            ! West boundary
            if(domain%boundary_exterior(4)) then
                domain%boundary_flux_store_exterior(4) = -sum(&
                    domain%flux_EW( 1+n_ext, (1+n_ext):(ny-n_ext), STG),&
                    mask = (&
                        domain%nesting%priority_domain_index(1+n_ext, (1+n_ext):(ny-n_ext)) == domain%nesting%my_index .and. &
                        domain%nesting%priority_domain_image(1+n_ext, (1+n_ext):(ny-n_ext)) == domain%nesting%my_image ))
            end if 
        else
            ! We are not doing nesting
            domain%boundary_flux_store_exterior = domain%boundary_flux_store
        end if

    end subroutine

!    !
!    ! Update the values of domain%U (i.e. the main flow variables), based on the fluxes and sources in domain
!    ! This is the same as domain%update_U, but the inner loop is instead written in vector form. Initial testing
!    ! suggested this is faster than the original version
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
!
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
!                (max(domain%depth(:,j), minimum_allowed_depth)**(NEG_SEVEN_ON_THREE_dp))
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
!    
!    END SUBROUTINE
!
!    !!!
!    !! This is an 'explicitly vectorized' version of extrapolate_edge_second_order
!    !!
!    !! Initial testing suggested this is faster than the original version
!    !!
!    !pure subroutine extrapolate_edge_second_order_vectorized(U_local, U_lower, U_upper, theta, &
!    !    extrapolation_sign, edge_value, n)
!    !    integer(ip), intent(in):: n
!    !    real(dp), intent(in):: U_local(n), U_lower(n), U_upper(n), theta(n) 
!    !    real(dp), intent(in):: extrapolation_sign(n)
!    !    real(dp), intent(out) :: edge_value(n)
!
!    !    integer(ip) :: i, imn, imx, vsize
!
!    !    ! Local 'small' vectors used to pack data and enhance vectorization
!    !    integer, parameter :: v = vectorization_size
!    !    real(dp):: c(v), d(v), a(v), b(v), e(v), th(v)
!
!
!    !    ! Strided loop by stride 'v', to help the compiler vectorize
!    !    do i = 1, n, v
!
!    !        ! Lower/upper indices
!    !        imn = i
!    !        imx = min(imn + v - 1, n)
!
!    !        ! If we have hit the end of the vector, then vsize < v, otherwise vsize=v
!    !        vsize = min(v, imx - imn + 1)
!
!    !        ! Pack data from input arrays into small arrays of size v. 
!    !        a(1:vsize) = U_upper(imn:imx) - U_local(imn:imx)
!    !        b(1:vsize) = U_local(imn:imx) - U_lower(imn:imx)
!    !        th(1:vsize) = theta(imn:imx)
!    !        !call minmod_sub(a, b, d) 
!    !        !d = minmod(a, b)
!    !        d = merge(min(abs(a), abs(b))*sign(ONE_dp,a), ZERO_dp, sign(ONE_dp,a) == sign(ONE_dp,b))
!
!    !        d = d * th ! Limit on the local gradient
!    !        e = HALF_dp * (a + b)
!    !        b = ZERO_dp
!    !        c = merge(b, e, d == ZERO_dp) 
!
!    !        ! NOTE: IF d /= 0, then clearly d, c have the same sign
!    !        ! We exploit this to avoid a further minmod call (which seems
!    !        ! expensive)
!    !        b = merge(min(c, d), max(c, d), d > ZERO_dp)
!
!    !        edge_value(imn:imx) = U_local(imn:imx) + HALF_dp * extrapolation_sign(imn:imx) * b(1:vsize)
!    !    end do
!    !        
!
!    !end subroutine
