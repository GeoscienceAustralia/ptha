!
! The code here solves the nonlinear shallow water equations in cartesian or spherical 
! coordinates, with or without coriolis. It is modified from the corresponding linear code.
!
! It has been moved out of domain_mod.f90 because: A) it became complex, and B) including this
! code in different subroutines with different parameters/preprocessor-variables facilitates
! efficient code generation
! 
!
! The subroutine header has been commented out


!    ! 
!    ! Nonlinear shallow water equations leap-frog update
!    !
!    ! Update domain%U by timestep dt, using the nonlinear shallow water equations.
!    !
!    ! @param domain the domain to advance
!    ! @param dt the timestep to advance. Should remain constant in between repeated calls
!    !     to the function (since the numerical method assumes constant timestep)
!    !
!    subroutine one_leapfrog_nonlinear_step(domain, dt)
!        class(domain_type), intent(inout):: domain
!        real(dp), intent(in):: dt

        ! Domain size
        integer(ip) :: nx, ny
        ! Useful constants for removing divisons from inner loop (e.g. dt/area)
        real(dp):: inv_cell_area_dt, inv_cell_area_dt_vh, inv_cell_area_dt_vh_g, inv_cell_area_dt_g
        real(dp):: fr_limit_uh, fr_limit_vh
        ! d(stage)/dy; depth_{i, j+1/2}, depth_{j, i+1/2, j}
        real(dp):: dw_j(domain%nx(1)), h_jph_vec(domain%nx(1)), h_iph_vec(domain%nx(1))
        real(dp) :: h_iph_wet_strict(domain%nx(1)), h_jph_wet_strict(domain%nx(1))
        real(dp) :: h_iph_wet(domain%nx(1)), h_jph_wet(domain%nx(1))
       
        integer(ip):: j, i, xl, xu, yl, yu, n_ext, my_omp_id, n_omp_threads, loop_work_count, jm1, im1, jp1, jp2
        integer(ip) :: yl_omp, yu_omp

        ! Vector to hold the (coriolis force x time-step x 0.5), which arises
        ! in an inner loop. Use a vector so we can zero it for dry cells
        real(dp):: dt_half_coriolis(domain%nx(1)), dt_half_coriolis_jph(domain%nx(1))

        ! Vector to hold the semi_implicit_friction_multiplier
        real(dp):: friction_multiplier_UH(domain%nx(1)), friction_multiplier_VH(domain%nx(1))

        ! For the coriolis term in the UH update, we store an interpolated version of
        ! the 'pre-update' VH at 'i+1/2', 'j-1/2'. We also need the same at 'i+1/2', 'j+1/2'
        real(dp) :: vh_iph_jmh(domain%nx(1)), vh_iph_jph(domain%nx(1))

        ! For the coriolis term in the VH update, we store and interpolated version of
        ! the 'pre-update' UH at 'i', 'j'. We also need that at 'i', 'j+1'
        real(dp) :: uh_i_j(domain%nx(1)), uh_i_jp1(domain%nx(1))
        real(dp) :: uh_i_jp1_max_loop_index(domain%nx(1))

        ! Work-array indices for domain%advection_work(:,:,index)
        ! FIXME: The code could be restructured to avoid the use of this memory.
        !        The current code was edited from the linear solver (which is much more memory efficient),
        !        and to simplify the implementation this large block of temporary storage was used.
        !        A re-write could follow the approach of the linear solver (where we store a few rows of data
        !        only, and take care with openmp to ensure the required data is available).
        integer(ip), parameter :: DEPTH_LAG = 1, UUH = 2, UVH = 3, VUH = 4, VVH = 5
        real(dp) :: depth_local, inv_depth_local

        domain%dt_last_update = dt


        ! For the nonlinear advection solver, we must store the previous stage
        domain%advection_work(:,:, DEPTH_LAG) = max(domain%U(:,:,STG) - domain%U(:,:,ELV), ZERO_dp)
        domain%advection_work(:,:, UUH) = ZERO_dp
        domain%advection_work(:,:, UVH) = ZERO_dp
        domain%advection_work(:,:, VUH) = ZERO_dp
        domain%advection_work(:,:, VVH) = ZERO_dp
        ! Clear uh/vh fluxes
        domain%flux_EW(:,:,UH:VH) = ZERO_dp
        domain%flux_NS(:,:,UH:VH) = ZERO_dp

        !
        ! idea: U(i, j, UH) = UH_{i+1/2, j}
        !     : U(i, j, VH) = VH_{i, j+1/2}
        ! We rely on boundary conditions for ALL boundary stage values, and
        ! apply the boundary condition after the stage update. We then update
        ! all uh/vh values which can be updated using those stage values [so avoid
        ! having to specify boundary conditions for them]

        ! If doing nesting, track fluxes through nesting regions (mass flux
        ! only, since momentum flux is zero by definition).
        ! For the leap-frog solver, the mass fluxes are given by the
        ! depth-integrated-velocity components.
        ! Note domain%U(:,1,VH) gives south-edge flux for cell j=2, and
        ! similarly domain%U(1,:,UH) gives the east-edge flux for cell i = 2.
        call domain%nesting_boundary_flux_integral_tstep(dt, &
           flux_NS=domain%U(:,:,VH:VH), flux_NS_lower_index=2_ip, &
           flux_EW=domain%U(:,:,UH:UH), flux_EW_lower_index=2_ip, &
           var_indices=[STG, STG], flux_already_multiplied_by_dx=.FALSE.)
        
        nx = domain%nx(1)
        ny = domain%nx(2)
     
        ! For now don't try to support evolving a subset of rows/columns
        ! (It may or may not work but needs to be tested).
        xL = 1 !domain%xL
        xU = nx !domain%xU
        yL = 1 !domain%yL
        yU = ny !domain%yU

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, dt, nx, ny, xL, xU, yL, yU)

        !
        ! Update stage
        !
        !$OMP DO SCHEDULE(STATIC)
        do j = (yL+1),(yU-1)
            ! For spherical coordinates, cell area changes with y.
            ! For cartesian coordinates this could be moved out of the loop
            inv_cell_area_dt = dt / domain%area_cell_y(j)

            do i = (xL+1) , (xU-1)
                ! dstage/dt = - 1/(R cos (lat)) [ d uh / dlon + dvhcos(lat)/dlat ]
                domain%U(i, j, STG) = domain%U(i, j, STG) - inv_cell_area_dt * ( &
                    (domain%U(i, j, UH) - domain%U(i-1, j, UH)) * &
                        domain%distance_left_edge(i) + &
                    (domain%U(i, j, VH)*domain%distance_bottom_edge(j+1) - &
                        domain%U(i, j-1, VH)*domain%distance_bottom_edge(j))&
                    )
        
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        domain%time = domain%time + dt * HALF_dp
        call domain%update_boundary()
        call domain%apply_forcing(dt)

        if(domain%use_partitioned_comms) call domain%partitioned_comms%communicate(domain%U)


        ! Boundary flux integration (relevant to single-domain case)
        ! note : U(i, j, UH) = UH_{i+1/2, j}
        !      : U(i, j, VH) = VH_{i, j+1/2}
        n_ext = domain%exterior_cells_width
        ! Outward boundary flux over the north -- integrate over the 'interior' cells
        domain%boundary_flux_store(1) = sum(domain%U((n_ext+1):(nx-n_ext),ny-n_ext,VH)) * &
            domain%distance_bottom_edge(ny-n_ext+1)
        ! Outward boundary flux over the east
        domain%boundary_flux_store(2) = sum(domain%U(nx-n_ext,(1+n_ext):(ny-n_ext),UH)) * &
            domain%distance_left_edge(nx-n_ext+1)
        ! Outward boundary flux over the south
        domain%boundary_flux_store(3) = -sum(domain%U((n_ext+1):(nx-n_ext),n_ext,VH)) * &
            domain%distance_bottom_edge(n_ext+1)
        ! Outward boundary flux over the west
        domain%boundary_flux_store(4) = -sum(domain%U(n_ext,(1+n_ext):(ny-n_ext),UH)) * &
            domain%distance_left_edge(n_ext+1)

        domain%boundary_flux_evolve_integral = sum(domain%boundary_flux_store) * dt

        !
        ! Boundary flux integration (relevant to multidomain case)
        !
        if(domain%nesting%my_index > 0) then
            ! We are in a multidomain -- carefully compute fluxes through
            ! exterior boundaries
            domain%boundary_flux_store_exterior = ZERO_dp
       
            ! Here we implement masked versions of the boundary flux sums above, only counting cells
            ! where the priority domain is receiving/sending the fluxes on actual physical boundaries 

            ! North boundary
            if(domain%boundary_exterior(1)) then
                domain%boundary_flux_store_exterior(1) = domain%distance_bottom_edge(ny-n_ext+1) * sum(&
                    domain%U( (1+n_ext):(nx-n_ext), ny-n_ext, VH),&
                    mask = (&
                        domain%nesting%priority_domain_index((1+n_ext):(nx-n_ext), ny-n_ext) == domain%nesting%my_index .and. &
                        domain%nesting%priority_domain_image((1+n_ext):(nx-n_ext), ny-n_ext) == domain%nesting%my_image ))
            end if

            ! East boundary
            if(domain%boundary_exterior(2)) then
                domain%boundary_flux_store_exterior(2) = domain%distance_left_edge(nx-n_ext+1) * sum(&
                    domain%U( nx-n_ext, (1+n_ext):(ny-n_ext), UH),&
                    mask = (&
                        domain%nesting%priority_domain_index(nx-n_ext, (1+n_ext):(ny-n_ext)) == domain%nesting%my_index .and. &
                        domain%nesting%priority_domain_image(nx-n_ext, (1+n_ext):(ny-n_ext)) == domain%nesting%my_image ))
            end if

            ! South boundary
            if(domain%boundary_exterior(3)) then
                domain%boundary_flux_store_exterior(3) = - domain%distance_bottom_edge(n_ext+1) * sum(&
                    domain%U( (1+n_ext):(nx-n_ext), n_ext, VH),&
                    mask = (&
                        domain%nesting%priority_domain_index((1+n_ext):(nx-n_ext), 1+n_ext) == domain%nesting%my_index .and. &
                        domain%nesting%priority_domain_image((1+n_ext):(nx-n_ext), 1+n_ext) == domain%nesting%my_image ))
            end if

            ! West boundary
            if(domain%boundary_exterior(4)) then
                domain%boundary_flux_store_exterior(4) = - domain%distance_left_edge(n_ext+1) * sum(&
                    domain%U(n_ext, (1+n_ext):(ny-n_ext), UH),&
                    mask = (&
                        domain%nesting%priority_domain_index(1+n_ext, (1+n_ext):(ny-n_ext)) == domain%nesting%my_index .and. &
                        domain%nesting%priority_domain_image(1+n_ext, (1+n_ext):(ny-n_ext)) == domain%nesting%my_image ))
            end if 
        else
            ! We are not doing nesting
            domain%boundary_flux_store_exterior = domain%boundary_flux_store
        end if

        domain%boundary_flux_evolve_integral_exterior = &
            sum(domain%boundary_flux_store_exterior, mask=domain%boundary_exterior) * dt

        !
        ! At this point, U(:,:,STG) updated everywhere (assuming xL/xU/yL/yU denote inactive parts of domain in an inclusive sense)
        ! Next update uh, vh from indices 1...(domain%nx - 1), which correspond to fluxes that have exterior stage values.
        !

        ! Predefine the friction terms. In the 'truely-linear' case we will only call this once,
        ! whereas in the not-truely-linear case it will be called every timestep (i.e. domain%friction_work_is_setup==FALSE)
        if(.not. domain%friction_work_is_setup) call precompute_friction_work(domain)

        ! NOTE: Here we manually determine the openmp loop indices. This is
        ! required to include the Coriolis terms without having to copy
        ! memory of size domain%U(:,:,UH:VH) -- which would increase memory use
        ! by 50% for the simple linear solver. 
        ! (For the nonlinear solver this is not useful because we store many things anyway)
        !
        loop_work_count = yU - yL ! Number of indices between 'yU - 1' and 'yL'
        !
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, dt, nx, ny, xL, xU, yL, yU, loop_work_count)
        !
#if defined(NOOPENMP)
        ! Serial loop from yL:(yU-1)
        yl_omp = yL ! 1
        yu_omp = yU - 1 ! ny - 1
#else
        ! Parallel loop from yL:(yU-1)
        !
        ! NOTE: In fortran, 
        !     DO i = start, end
        !       .... code here ...
        !     END DO
        ! will not do anything at all if start > end [that's the standard, if
        ! we don't specify a stride other than 1].
        ! This behaviour is important for ensuring the openmp loop sharing
        ! code here works, even if e.g. we have more omp threads than loop
        ! indices. In that case, we will have some threads with yl_omp >
        ! yu_omp -- which will not loop -- and that is the desired behaviour.
        !
        my_omp_id = omp_get_thread_num()
        n_omp_threads = omp_get_num_threads()
        ! yl_omp = lower loop index of the thread my_omp_id
        yl_omp = nint(loop_work_count * my_omp_id * 1.0_dp / n_omp_threads) + yL ! Like '1'
        ! yu_omp = upper loop index of the thread my_omp_id
        yu_omp = nint(loop_work_count * (my_omp_id + 1) * 1.0_dp / n_omp_threads) + yL - 1 ! Like 'ny - 1'
#endif
        !
        ! Now loop indices are determined
        !

        ! Store the depth (i,j) at the time of the previous momentum update
        do j = yl_omp, yu_omp
            domain%advection_work(:,j,DEPTH_LAG) = 0.5_dp * (domain%advection_work(:,j,DEPTH_LAG) + &
                max(domain%U(:,j,STG) - domain%U(:,j,ELV), ZERO_dp) )
        end do

        ! 
        ! Get uuh and uvh terms at {i+1/2, j} for i-inds 1...(domain%nx(1)) and j-inds 1...(domain%nx(2))
        ! Extrapolate as required to cover those indices
        !
        do j = yl_omp, yu_omp
            jm1 = max(j-1, 1) ! Avoid out-of-bounds for U*VH when j=1 {questionable value}
            do i = xL, xU-1
                ! depth at i+1/2, j
                depth_local = 0.5_dp * (domain%advection_work(i, j, DEPTH_LAG) + domain%advection_work(i+1, j, DEPTH_LAG))
                inv_depth_local = merge(1.0_dp/depth_local, ZERO_dp, depth_local > minimum_allowed_depth)

                ! U*UH at i+1/2, j
                domain%advection_work(i,j,UUH) = (domain%U(i,j,UH) * domain%U(i,j,UH)) * inv_depth_local

                ! U*VH at i+1/2, j. 
                if(i < nx-1) domain%advection_work(i,j,UVH) = domain%U(i,j,UH) * inv_depth_local * 0.25_dp * &
                    (domain%U(i,jm1,VH) + domain%U(i+1,jm1,VH) + domain%U(i,j,VH) + domain%U(i+1,j,VH))
             end do
             ! Extrapolation near boundary. If xU != nx, then anyway the values should be zero, so do nothing.
             if(xU == nx) then
                 domain%advection_work(nx  ,j,UUH) = domain%advection_work(nx-1, j, UUH)
                 domain%advection_work(nx  ,j,UVH) = domain%advection_work(nx-2, j, UVH)
                 domain%advection_work(nx-1,j,UVH) = domain%advection_work(nx-2, j, UVH)
             end if
        end do
        ! Barrier to ensure that ny-1 has been updated
        !$OMP BARRIER
        ! Extraoplation near boundary
        if(yu_omp == ny-1) then
            domain%advection_work(:, ny,UUH) = domain%advection_work(:, ny-1, UUH)
            domain%advection_work(:, ny,UVH) = domain%advection_work(:, ny-1, UVH)
        end if

        !
        ! Get vuh and vvh terms at {i+1/2, j+1/2} for i-inds 1...(domain%nx(1)) and j-inds 1...(domain%nx(2))
        ! Extrapolate as required
        !
        do j = yl_omp, yu_omp
             do i = xL, xU
                ! depth at i, j+1/2
                depth_local = 0.5_dp * (domain%advection_work(i, j, DEPTH_LAG) + domain%advection_work(i, j+1, DEPTH_LAG))
                inv_depth_local = merge(1.0_dp / depth_local, ZERO_dp, depth_local > minimum_allowed_depth)

                ! V*VH at i, j+1/2
                domain%advection_work(i,j,VVH) = domain%U(i,j,VH) * domain%U(i, j, VH) * inv_depth_local 

                ! V*UH at i, j+1/2. 
                ! Avoid out-of-bounds when xL=1
                if(i > 1 .and. j < (ny - 1)) domain%advection_work(i,j,VUH) = domain%U(i,j,VH) * inv_depth_local * 0.25_dp * &
                     (domain%U(i-1,j,UH) + domain%U(i,j,UH) + domain%U(i-1,j+1,UH) + domain%U(i,j+1,UH))
            end do
            ! Extrapolate near boundary for V*UH only. If xL != 1, then anyway the term should be zero
            if(xL == 1) domain%advection_work(1, j, VUH) = domain%advection_work(2, j, VUH)
        end do
        ! Barrier to ensure ny-1, ny-2 have been updated
        !$OMP BARRIER
        ! Extrapolate near boundary.
        if(yu_omp == ny-1) then
            domain%advection_work(xL:xU, ny  , VVH) = domain%advection_work(xL:xU,ny-1,VVH)
            ! Zero gradient extrapolation
            domain%advection_work(xL:xU, ny  , VUH) = domain%advection_work(xL:xU,ny-2,VUH)
            domain%advection_work(xL:xU, ny-1, VUH) = domain%advection_work(xL:xU,ny-2,VUH)
        end if

        ! Compute fluxes once advection_work terms are all computed
        !$OMP BARRIER

        !
        ! Simple upwind advection as in Tunawi and other schemes
        !

        ! EW flux of UH
        do j = yl_omp, yu_omp
            ! NB: distance_left_edge is constant, so do not worry about shifting for staggered grid
            domain%flux_EW(xL+1:xU-1, j, UH) = domain%distance_left_edge(xL+1:xU-1) * merge(&
                domain%advection_work(xL  :xU-2, j, UUH), & 
                domain%advection_work(xL+1:xU-1, j, UUH), &
                domain%U(xL+1:xU-1, j, UH) + domain%U(xL:xU-2, j, UH) > ZERO_dp)
            ! Extrapolate if we really are at the boundaries -- noting "distance_left_edge = constant"
            if(xL == 1 ) domain%flux_EW(1 ,j, UH) = domain%advection_work( 1, j, UUH) * domain%distance_left_edge( 1)
            if(xU == nx) domain%flux_EW(nx,j, UH) = domain%advection_work(nx, j, UUH) * domain%distance_left_edge(nx)
            !if(xL == 1 ) domain%flux_EW(1 ,j, UH) = domain%flux_EW(2,j,UH)
            !if(xU == nx ) domain%flux_EW(nx ,j, UH) = domain%flux_EW(nx-1,j,UH)


            ! Treat case where neighbour UH is of opposite sign and equal
            domain%flux_EW(xL+1:xU-1, j, UH) = merge(&
                domain%distance_left_edge(xL+1:xU-1) * 0.5_dp * &
                    (domain%advection_work(xL  :xU-2, j, UUH) + &
                     domain%advection_work(xL+1:xU-1, j, UUH)), &
                domain%flux_EW(xL+1:xU-1, j, UH), &
                (domain%U(xL+1:xU-1, j, UH) + domain%U(xL:xU-2, j, UH)) == ZERO_dp)

            ! Wetting and drying
            domain%flux_EW(xL:xU, j, UH) = merge(domain%flux_EW(xL:xU, j, UH), &
                ZERO_dp, &
                domain%U(xL:xU, j, STG) > domain%U(xL:xU, j, ELV) + minimum_allowed_depth)

        end do

        ! NS flux of uh -- here we work on the top-edge
        do j = yl_omp, yu_omp
            domain%flux_NS(xL:xU-1, j+1, UH) = domain%distance_bottom_edge(j+1) * merge(&
                domain%advection_work(xL:xU-1, j+0, UVH), &
                domain%advection_work(xL:xU-1, j+1, UVH), &
                !vh_iph_jph(xL:xU-1) > ZERO_dp)
                domain%U(xL:xU-1, j, VH) + domain%U(xL+1:xU, j, VH) > ZERO_dp)

            ! Treat case where neighbour VH is of opposite sign and equal
            domain%flux_NS(xL:xU-1, j+1, UH) = merge(&
                domain%distance_bottom_edge(j+1) * 0.5_dp * &
                    (domain%advection_work(xL:xU-1, j+0, UVH) + domain%advection_work(xL:xU-1, j+1, UVH)), &
                domain%flux_NS(xL:xU-1, j+1, UH), &
                (domain%U(xL:xU-1, j, VH) + domain%U(xL+1:xU, j, VH)) == ZERO_dp)

            ! Wetting and drying
            domain%flux_NS(xL:xU-1, j+1, UH) = merge(&
                domain%flux_NS(xL:xU-1, j+1, UH), &
                ZERO_dp, &
                (domain%U(xL:xU-1, j  , STG) > domain%U(xL:xU-1, j  , ELV) + minimum_allowed_depth) .and. &
                (domain%U(xL+1:xU, j  , STG) > domain%U(xL+1:xU, j  , ELV) + minimum_allowed_depth) .and. &
                (domain%U(xL:xU-1, j+1, STG) > domain%U(xL:xU-1, j+1, ELV) + minimum_allowed_depth) .and. &
                (domain%U(xL+1:xU, j+1, STG) > domain%U(xL+1:xU, j+1, ELV) + minimum_allowed_depth) )

        end do

        ! Take care of the bottom edge that we missed
        if(yl_omp == yL) then
            j = yL
            jm1 = max(j-1, 1) ! Safe j-1
            domain%flux_NS(xL:xU-1, j, UH) = domain%distance_bottom_edge(j) * merge(&
                domain%advection_work(xL:xU-1, jm1, UVH), &
                domain%advection_work(xL:xU-1, j  , UVH), &
                domain%U(xL:xU-1, jm1, VH) + domain%U(xL+1:xU, jm1, VH) > ZERO_dp)

            ! Treat case where neighbour VH is of opposite sign and equal
            domain%flux_NS(xL:xU-1, j, UH) = merge(&
                domain%distance_bottom_edge(j) * 0.5_dp * &
                    (domain%advection_work(xL:xU-1, jm1, UVH) + domain%advection_work(xL:xU-1, j  , UVH)), &
                domain%flux_NS(xL:xU-1, j, UH), &
                (domain%U(xL:xU-1, jm1, VH) + domain%U(xL+1:xU, jm1, VH)) == ZERO_dp)

            ! Wetting and drying
            domain%flux_NS(xL:xU-1, j, UH) = merge(&
                domain%flux_NS(xL:xU-1, j, UH), &
                ZERO_dp, &
                (domain%U(xL:xU-1, j  , STG) > domain%U(xL:xU-1, j  , ELV) + minimum_allowed_depth) .and. &
                (domain%U(xL+1:xU, j  , STG) > domain%U(xL+1:xU, j  , ELV) + minimum_allowed_depth) .and. &
                (domain%U(xL:xU-1, jm1, STG) > domain%U(xL:xU-1, jm1, ELV) + minimum_allowed_depth) .and. &
                (domain%U(xL+1:xU, jm1, STG) > domain%U(xL+1:xU, jm1, ELV) + minimum_allowed_depth) )
        end if

        ! NS flux of vh -- here we work with the bottom edge
        do j = yl_omp, yu_omp
            jm1 = max(j-1, 1)
            domain%flux_NS(xL:xU-1, j, VH) = &
                0.5_dp * (domain%distance_bottom_edge(j) + domain%distance_bottom_edge(j+1)) * merge(&
                    domain%advection_work(xL:xU-1, jm1, VVH), &
                    domain%advection_work(xL:xU-1, j  , VVH), &
                    domain%U(xL:xU-1, j, VH) + domain%U(xL:xU-1, jm1, VH) > ZERO_dp)

            ! Treat case where neighbour VH is of opposite sign and equal
            domain%flux_NS(xL:xU-1, j, VH) = merge(&
                    0.5_dp * (domain%distance_bottom_edge(j) + domain%distance_bottom_edge(j+1)) * 0.5_dp * &
                        (domain%advection_work(xL:xU-1, jm1, VVH) + domain%advection_work(xL:xU-1, j  , VVH)), &
                    domain%flux_NS(xL:xU-1, j, VH), &
                    (domain%U(xL:xU-1, j, VH) + domain%U(xL:xU-1, jm1, VH)) == ZERO_dp)

            ! Wetting and drying
            domain%flux_NS(xL:xU-1, j, VH) = merge(domain%flux_NS(xL:xU-1, j, VH), &
                ZERO_dp, &
                domain%U(xL:xU-1, j, STG) > domain%U(xL:xU-1, j, ELV) + minimum_allowed_depth)
        end do

        ! Take care of top-edge we missed
        if(yu_omp == (yU - 1)) then
            j = yu_omp + 1
            jp1 = min(j+1, ny+1) ! Safe j+1
            domain%flux_NS(xL:xU-1, j, VH) = &
                0.5_dp * (domain%distance_bottom_edge(j) + domain%distance_bottom_edge(jp1)) * &
                merge(domain%advection_work(xL:xU-1, j-1, VVH), &
                      domain%advection_work(xL:xU-1, j  , VVH), &
                      domain%U(xL:xU-1, j, VH) + domain%U(xL:xU-1, j-1, VH)  > ZERO_dp)

            ! Treat case where neighbour VH is of opposite sign and equal
            domain%flux_NS(xL:xU-1, j, VH) = merge(&
                0.5_dp * (domain%distance_bottom_edge(j) + domain%distance_bottom_edge(jp1)) * &
                    0.5_dp * (domain%advection_work(xL:xU-1, j-1  , VVH) + domain%advection_work(xL:xU-1, j, VVH)), &
                domain%flux_NS(xL:xU-1, j, VH), &
                (domain%U(xL:xU-1, j, VH) + domain%U(xL:xU-1, j-1, VH))  == ZERO_dp)

            ! Wetting and drying
            domain%flux_NS(xL:xU-1, j, VH) = merge(domain%flux_NS(xL:xU-1, j, VH), &
                ZERO_dp, &
                domain%U(xL:xU-1, j, STG) > domain%U(xL:xU-1, j, ELV) + minimum_allowed_depth)

        end if

        do j = yl_omp, yu_omp
            ! EW flux of vh
            ! distance_left_edge = constant for all i
            domain%flux_EW(xL+1:xU, j, VH) = domain%distance_left_edge(xL:xU-1) * merge(&
                domain%advection_work(xL  :xU-1, j, VUH), &
                domain%advection_work(xL+1:xU  , j, VUH), &
                domain%U(xL:xU-1,j,UH) + domain%U(xL:xU-1,j+1,UH) > ZERO_dp) 
            if(xL == 1) domain%flux_EW(1, j, VH) = domain%flux_EW(2, j, VH)

            ! Treat case where neighbour UH is of opposite sign and equal
            domain%flux_EW(xL+1:xU, j, VH) = merge(&
                 domain%distance_left_edge(xL:xU-1) * 0.5_dp * &
                     (domain%advection_work(xL  :xU-1, j, VUH) + domain%advection_work(xL+1:xU  , j, VUH)), &
                 domain%flux_EW(xL+1:xU, j, VH), &
                 (domain%U(xL:xU-1,j,UH) + domain%U(xL:xU-1,j+1,UH)) == ZERO_dp)

            ! Wetting and drying 
            domain%flux_EW(xL+1:xU, j, VH) = merge(&
                domain%flux_EW(xL+1:xU, j, VH), &
                ZERO_dp, &
                (domain%U(xL:xU-1, j  , STG) > domain%U(xL:xU-1, j  , ELV) + minimum_allowed_depth) .and. &
                (domain%U(xL+1:xU, j  , STG) > domain%U(xL+1:xU, j  , ELV) + minimum_allowed_depth) .and. &
                (domain%U(xL:xU-1, j+1, STG) > domain%U(xL:xU-1, j+1, ELV) + minimum_allowed_depth) .and. &
                (domain%U(xL+1:xU, j+1, STG) > domain%U(xL+1:xU, j+1, ELV) + minimum_allowed_depth) )

        end do

        ! Tricks to implement coriolis without increasing memory usage much
        !
        ! For coriolis, we need values of 'VH' at coorinates corresponding to UH,
        ! and also values of UH at coordinates corresponding to VH. 
        !
        ! This requires 4 point averaging of the 'OLD' values of UH and VH.
        ! 
        ! A simply way to do that is to store all the OLD values -- but that involves
        ! lots of memory. Or,
        ! We can do this in a loop with care, without storing those 'OLD' values,
        ! so long as we store values needed at omp thread loop boundaries, and
        ! control the order that the 'j' indices are iterated, and generally be
        ! careful about openmp.

        !
        ! Zero vectors which hold strategically stored information
        !
        vh_iph_jmh(xL:(xU-1)) = ZERO_dp
        vh_iph_jph(xL:(xU-1)) = ZERO_dp
        uh_i_j(xL:(xU-1)) = ZERO_dp
        uh_i_jp1(xL:(xU-1)) = ZERO_dp
        uh_i_jp1_max_loop_index(xL:(xU-1)) = ZERO_dp
        dt_half_coriolis(xL:(xU-1)) = ZERO_dp
        dt_half_coriolis_jph(xL:(xU-1)) = ZERO_dp

        !
        ! Before starting the loop, get the value of VH at 
        ! 'i+1/2', 'yl_omp-1/2', to prevent the possibility that it is
        ! updated by another openmp thread before we read it !
        if(yl_omp == 1) then
            ! In this case, there is no 'j-1/2' value, just use the boundary value
            vh_iph_jmh(xL:(xU-1)) = HALF_dp * &
                (domain%U(xL    :(xU-1), 1, VH) + &
                 domain%U((xL+1):xU    , 1, VH))
        else
            ! Average values at 'i' and 'i+1'
            vh_iph_jmh(xL:(xU-1)) = HALF_dp * &
                (domain%U(xL    :(xU-1), yl_omp - 1, VH) + &
                 domain%U((xL+1):xU    , yl_omp - 1, VH))
        end if

        !
        ! Before starting the loop, store the value of UH at 'i', 'yu_omp+1', to prevent the
        ! possibility that it is updated by another openmp thread before we read it
        !
        if(xL > 1) then
            ! Average values at 'i-1/2' and 'i+1/2'
            uh_i_jp1_max_loop_index(xL:(xU-1)) = HALF_dp * &
                (domain%U( (xL-1):(xU-2), yu_omp + 1, UH) + &
                 domain%U( xL    :(xU-1), yu_omp + 1, UH))
        else
            ! Avoid out-of bounds error
            uh_i_jp1_max_loop_index((xL+1):(xU-1)) = HALF_dp * &
                (domain%U( xL   :(xU-2), yu_omp + 1, UH) + &
                 domain%U( xL+1 :(xU-1), yu_omp + 1, UH))
            uh_i_jp1_max_loop_index(xL) = domain%U(1, yu_omp + 1, UH)
        end if

        !
        ! Get the initial value of UH at 'i, yl_omp'. This 'starting value' is required because 
        ! inside the loop we set 'uh_i_j = uh_i_jph' after the update.
        !
        uh_i_j((xL+1):(xU-1)) = HALF_dp * &
            ( domain%U((xL+1):(xU-1), yl_omp, UH) + &
              domain%U(    xL:(xU-2), yl_omp, UH) )
        !
        ! Special case of xL which is protective if xL == 1
        uh_i_j(xL) = HALF_dp * ( domain%U(xL, yl_omp, UH) + &
            domain%U(max(xL-1, 1), yl_omp, UH) )


        ! No thread can start updating the loop until all threads have their
        ! 'bounding' coriolis terms
        !$OMP BARRIER
        ! BEWARE -- After the above barrier, we cannot safetly work with domain%U(:,:, UH:VH) EXCEPT when the j-index = j,
        ! because it may have been updated on another openmp process. However we can work with domain%U(:,:,STG), and of
        ! course ELV

        ! Main update-flux loop
        !do j = 1, ny - 1
        do j = yl_omp, yu_omp

            ! For spherical coordiantes, cell area changes with y.
            ! For cartesian coordinates this could be moved out of the loop
            inv_cell_area_dt = dt / domain%area_cell_y(j)
            inv_cell_area_dt_g = gravity * inv_cell_area_dt
            !
            ! For leap-frog with spherical coordinates, the area associated
            ! with the NS momentum term is different to the area associated
            ! with the stage and EW momentum term
            inv_cell_area_dt_vh = dt / &
                (HALF_dp * (domain%area_cell_y(j) + domain%area_cell_y(j+1)))
            inv_cell_area_dt_vh_g = gravity * inv_cell_area_dt_vh
       
            ! 
            ! Try to keep control-flow and non-local memory jumps out of inner loop
            ! This improves speed on my machine with gfortran (11/08/2016)
            !
            dw_j(xL:(xU-1)) = domain%U(xL:(xU-1), j+1, STG) - domain%U(xL:(xU-1), j, STG)

            !
            ! In the g * d * dStage/dx type term, let d vary. This implies nonlinearity
            !

            ! Depth at j-plus-half
            h_jph_vec(xL:xU) = HALF_dp * ((domain%U(xL:xU,j+1,STG) + domain%U(xL:xU,j,STG)) - &
                           (domain%U(xL:xU,j+1,ELV) + domain%U(xL:xU,j,ELV)))
            where(h_jph_vec(xL:xU) < minimum_allowed_depth) h_jph_vec(xL:xU) = ZERO_dp

            ! Depth at i-plus-half
            h_iph_vec(xL:(xU-1)) = &
                HALF_dp * ((domain%U((xL+1):xU, j, STG) + domain%U(xL:(xU-1), j, STG)) -&
                           (domain%U((xL+1):xU, j, ELV) + domain%U(xL:(xU-1), j, ELV)))
            where(h_iph_vec(xL:xU-1) < minimum_allowed_depth) h_iph_vec(xL:xU-1) = ZERO_dp

            !
            ! Wet/dry criteria
            !

            h_jph_wet_strict(xL:xU    ) = merge(ONE_dp, ZERO_dp, &
                ! Both cells wet
                ((domain%U(xL:xU,j+1,STG) > (domain%U(xL:xU,j+1,ELV) + minimum_allowed_depth)) .and. &
                 (domain%U(xL:xU,j+0,STG) > (domain%U(xL:xU,j+0,ELV) + minimum_allowed_depth)) ) )

            h_jph_wet(xL:xU    ) = merge(ONE_dp, ZERO_dp, &
                ! Both cells wet, OR
                (h_jph_wet_strict(xL:xU) > ZERO_dp) .or. &
                ! Wet cell (@ j+1) has higher stage than dry cell (@ j), OR
                ((domain%U(xL:xU,j+1,STG) > (domain%U(xL:xU,j+1,ELV) + minimum_allowed_depth)) .and. &
                 (domain%U(xL:xU,j+1,STG) > (domain%U(xL:xU,j+0,ELV) + minimum_allowed_depth)) ) .or. &
                ! Wet cell (@ j) has higher stage than dry cell (@ j+1)
                ((domain%U(xL:xU,j+0,STG) > (domain%U(xL:xU,j+0,ELV) + minimum_allowed_depth)) .and. &
                 (domain%U(xL:xU,j+0,STG) > (domain%U(xL:xU,j+1,ELV) + minimum_allowed_depth)) ) )

            h_iph_wet_strict(xL:(xU-1)) = merge(ONE_dp, ZERO_dp, &
                ! Both cells wet
                ((domain%U(xL+1:xU,j,STG) > (domain%U(xL+1:xU,j,ELV) + minimum_allowed_depth)) .and. &
                 (domain%U(xL:xU-1,j,STG) > (domain%U(xL:xU-1,j,ELV) + minimum_allowed_depth)) ) )

            h_iph_wet(xL:(xU-1)) = merge(ONE_dp, ZERO_dp, &
                ! Both cells wet, OR
                (h_iph_wet_strict(xL:(xU-1)) > ZERO_dp) .or. &
                ! Wet cell (@ i+1) has higher stage than dry cell (@ i)
                ((domain%U(xL+1:xU,j,STG) > (domain%U(xL+1:xU,j,ELV) + minimum_allowed_depth)) .and. &
                 (domain%U(xL+1:xU,j,STG) > (domain%U(xL:xU-1,j,ELV) + minimum_allowed_depth)) ) .or. &
                ! Wet cell (@ i) has higher stage than dry cell (@ i+1)
                ((domain%U(xL:xU-1,j,STG) > (domain%U(xL:xU-1,j,ELV) + minimum_allowed_depth)) .and. &
                 (domain%U(xL:xU-1,j,STG) > (domain%U(xL+1:xU,j,ELV) + minimum_allowed_depth)) ) )


            ! 'Old' VH at (i+1/2, j+1/2) -- requires averaging values at 'i,j+1/2' and 'i+1, j+1/2'
            vh_iph_jph(xL:(xU-1)) = HALF_dp * &
                (domain%U(xL    :(xU-1), j, VH) + &
                 domain%U((xL+1):xU    , j, VH))

            ! 'Old' UH at (i, j+1). 
            ! First get it assuming we are not at the last loop index -- and then fix it 
            ! -- and account for the case that j+1 is on another openmp thread.
            ! Step1: Get everything except xL, which needs special treatment if xL == 1
            uh_i_jp1((xL+1):(xU-1)) = HALF_dp * &
                (domain%U((xL+1):(xU-1), j+1, UH) + &
                 domain%U(    xL:(xU-2), j+1, UH) )
            ! Special case of xL which is protective if xL == 1
            uh_i_jp1(xL) = HALF_dp * ( domain%U(xL          , j+1, UH) + &
                                       domain%U(max(xL-1, 1), j+1, UH) )
            ! Final step to ensure that if (j == yu_omp), so j+1 is on another openmp 
            ! thread, then we still get the non-updated value of UH.
            uh_i_jp1(xl:(xU-1)) = merge(&
                               uh_i_jp1(xL:(xU-1)), &
                uh_i_jp1_max_loop_index(xL:(xU-1)), &
                j /= yu_omp)

#if defined(CORIOLIS)
            ! Avoid recomputing coriolis coefficient in loop. Note we set the
            ! coriolis force to zero over dry cells, as is obviously desirable.
            dt_half_coriolis(xL:(xU-1)) = dt * domain%coriolis(j) * HALF_dp * &
                merge(ONE_dp, ZERO_dp, h_iph_wet_strict(xL:(xU-1)) > ZERO_dp)
            dt_half_coriolis_jph(xL:(xU-1)) = dt * &
                domain%coriolis_bottom_edge(j+1) * HALF_dp * &
                merge(ONE_dp, ZERO_dp, h_jph_wet_strict(xL:(xU-1)) > ZERO_dp)
#endif

            ! Compute a multiplier that applies a semi-implicit treatment of manning friction.
            !
            ! This is equivalent to adding a term "g * depth * friction_slope" to the equations.
            !   where friction_slope = {manning_sq * depth**(-4./3.) * speed * velocity_component}
            !
            ! The term is broken up into a fast, semi-implicit discretization. 
            ! For UH, we break it up as (terms in { })
            !    g * depth * friction_slope = {g * manning_sq * depth^(-7/3)} * { sqrt(uh^2 + vh^2) } * {uh}
            ! while for VH, the final term is vh
            !    g * depth * friction_slope = {g * manning_sq * depth^(-7/3)} * { sqrt(uh^2 + vh^2) } * {vh}
            ! The discretization is semi-implicit as follows:
            !    - The component 'g n^2 depth^(-7/3)' is a constant in time for the "truely-linear" solver, 
            !      because the effective depth never changes in this "linear" framework. Otherwise it is updated
            !      every time-step.
            !    - The term sqrt(uh^2+vh^2) {= speed*depth} is treated explicitly.
            !    - The remaining 'uh' or 'vh' term is treated implicitly. 
            !
            ! Thus, appending the term g * depth * friction slope to the equations can be reduced to
            ! a multiplication of the form { 1/(1 + explicit_part_of_friction_terms) }
            !
            !
            friction_multiplier_UH(xL:(xU-1)) = ONE_dp / ( ONE_dp + &
                ! Linear friction like in Fine et al (2012), Kulikov et al (2014)
                dt * domain%linear_friction_coeff + &
                ! Nonlinear friction
                dt * domain%friction_work(xL:(xU-1), j, UH) * &
                sqrt(domain%U(xL:(xU-1),j,UH)**2 + (0.5_dp * ( vh_iph_jmh(xL:(xU-1)) + vh_iph_jph(xL:(xU-1)) ) )**2 ) )

            friction_multiplier_VH(xL:(xU-1)) = ONE_dp / (ONE_dp + &
                ! Linear friction like in Fine et al (2012), Kulikov et al (2014)
                dt * domain%linear_friction_coeff + &
                ! Nonlinear friction
                dt * domain%friction_work(xL:(xU-1), j, VH) * &
                sqrt(domain%U(xL:(xU-1),j,VH)**2 + (0.5_dp * ( uh_i_j(xL:(xU-1)) + uh_i_jp1(xL:(xU-1)) ) )**2 ) )

            do i = xL , (xU-1)
                ! Update without coriolis or friction.

                ! duh/dt = - g * h0/(R cos (lat)) [ d stage / dlon ]
                domain%U(i, j, UH) = domain%U(i, j, UH) * h_iph_wet(i) - &
                    inv_cell_area_dt_g * h_iph_vec(i) * h_iph_wet(i) * &
                    (domain%U(i+1, j, STG) - domain%U(i, j, STG)) * &
                    domain%distance_left_edge(i+1)

                ! dvh/dt = - g * h0/(R) [ d stage / dlat ]
                domain%U(i, j, VH) = domain%U(i, j, VH) * h_jph_wet(i)  - &
                    inv_cell_area_dt_vh_g * h_jph_vec(i) * h_jph_wet(i) * &
                    dw_j(i) * domain%distance_bottom_edge(j+1)

                ! Append nonlinear advection terms
                domain%U(i, j, UH) = domain%U(i, j, UH) - &
                    inv_cell_area_dt * advection_beta * h_iph_wet(i) * (&
                        (domain%flux_EW(i+1, j  , UH) - domain%flux_EW(i, j, UH) ) + &
                        (domain%flux_NS(i  , j+1, UH) - domain%flux_NS(i, j, UH)))

                domain%U(i, j, VH) = domain%U(i, j, VH) - &
                    inv_cell_area_dt_vh * advection_beta * h_jph_wet(i) * (&
                        (domain%flux_EW(i+1,  j, VH) - domain%flux_EW(i, j, VH)) + &
                        (domain%flux_NS(i  ,j+1, VH) - domain%flux_NS(i, j, VH)))

#if defined(CORIOLIS)
                ! Append coriolis term
                ! duh/dt = f*vh
                domain%U(i, j, UH) = domain%U(i, j, UH)  &
                    + dt_half_coriolis(i) * (vh_iph_jmh(i) + vh_iph_jph(i))

                ! dvh/dt = -f*uh
                domain%U(i, j, VH) = domain%U(i, j, VH)  &
                    - dt_half_coriolis_jph(i) * (uh_i_j(i) + uh_i_jp1(i))

#endif

                ! Add the semi-implicit friction terms
                domain%U(i,j,UH) = friction_multiplier_UH(i) * domain%U(i,j,UH)
                domain%U(i,j,VH) = friction_multiplier_VH(i) * domain%U(i,j,VH)

            end do

#if defined(FROUDE_LIMITING)
            ! Froude-number limiter. Something like this is done in the JAGURS code.
            do i = xL , (xU-1)
                ! Limit UH
                depth_local = 0.5_dp * (domain%U(i  ,j,STG) - domain%U(i  ,j,ELV) + &
                                        domain%U(i+1,j,STG) - domain%U(i+1,j,ELV))
                if(depth_local < minimum_allowed_depth) depth_local = ZERO_dp
                fr_limit_uh = domain%leapfrog_froude_limit * sqrt(depth_local * gravity) * &
                    depth_local * sign(1.0_dp, domain%U(i,j,UH))
                if(domain%U(i,j,UH) * domain%U(i,j,UH) > fr_limit_uh * fr_limit_uh) domain%U(i,j,UH) = fr_limit_uh
            end do
            do i = xL , (xU-1)
                ! Limit VH 
                depth_local = 0.5_dp * (domain%U(i,j  ,STG) - domain%U(i,j  ,ELV) + &
                                        domain%U(i,j+1,STG) - domain%U(i,j+1,ELV))
                if(depth_local < minimum_allowed_depth) depth_local = ZERO_dp
                fr_limit_vh = domain%leapfrog_froude_limit * sqrt(depth_local * gravity) * &
                    depth_local * sign(1.0_dp, domain%U(i,j,VH))
                if(domain%U(i,j,VH) * domain%U(i,j,VH) > fr_limit_vh * fr_limit_vh) domain%U(i,j,VH) = fr_limit_vh
            end do
#endif                      

            ! On the next j iteration, the 'old' value of VH at i+1/2,
            ! j-1/2 can be derived using the current value of VH at i+1, j+1/2
            vh_iph_jmh(xL:(xU-1)) = vh_iph_jph(xL:(xU-1))
            uh_i_j(xL:(xU-1)) = uh_i_jp1(xL:(xU-1))

        end do

        ! For nonlinear solvers we apply a wetting-and-drying treatment
        do j = yl_omp, yu_omp                
            ! 
            ! Zero flux outflow from dry cells
            !
            do i = xL, xU !(xU-1)
                if(domain%U(i,j,STG) < domain%U(i,j,ELV) + minimum_allowed_depth) then
                    ! No UH outflow from dry cell to the west
                    if(domain%U(i,j,UH) > ZERO_dp) domain%U(i,j,UH) = ZERO_dp
                    if(i > 1) then
                        ! No UH outflow from dry cell to the east
                        if(domain%U(i-1, j, UH) < ZERO_dp) domain%U(i-1, j, UH) = ZERO_dp
                    end if

                    ! No VH outflow from dry cell to the south
                    if(domain%U(i,j,VH) > ZERO_dp) domain%U(i,j,VH) = ZERO_dp
                end if

                ! No VH outflow from dry cell to the north
                if(domain%U(i,j+1,STG) < domain%U(i,j+1,ELV) + minimum_allowed_depth) then
                    if(domain%U(i,j,VH) < ZERO_dp) domain%U(i,j,VH) = ZERO_dp
                end if
            end do
        end do
            

        !@!$OMP END DO
        !$OMP END PARALLEL
        domain%time = domain%time + HALF_dp*dt


contains

    ! Partial computation of friction term 
    !
    ! This assumes stage/UH/VH/elev have been set
    !
    ! The friction work term is of the form
    !     g * n^2 / constant_depth^(7/3)
    ! The key point is that when multiplied by ||UH|| * uh, it will be equal to the
    ! standard friction form: g*constant_depth*friction_slope
    !
    subroutine precompute_friction_work(domain)
        type(domain_type), intent(inout) :: domain

        integer(ip) :: i, j, jp1, ip1
        real(dp) :: depth_iph, depth_jph, nsq_iph, nsq_jph

        real(dp), parameter :: manning_depth_power = NEG_SEVEN_ON_THREE_dp, chezy_depth_power = -2.0_dp
        ! The "#include" code below can also work for a "linear with nonlinear friction" model, but here
        ! that is not done
        logical, parameter :: truely_linear = .FALSE.

        if(domain%friction_type == 'manning') then

#define _FRICTION_DEPTH_POWER_ manning_depth_power
#include "domain_mod_leapfrog_solver_friction_include.f90"
#undef _FRICTION_DEPTH_POWER_

        else if(domain%friction_type == 'chezy') then

#define _FRICTION_DEPTH_POWER_ chezy_depth_power
#include "domain_mod_leapfrog_solver_friction_include.f90"
#undef _FRICTION_DEPTH_POWER_

        end if

    end subroutine

!    end subroutine

