!
! The code here solves the linear shallow water equations in cartesian or spherical 
! coordinates, with or without coriolis. It can also include a nonlinear friction term,
! (which may be an efficient choice for modelling large scale tsunami propagation with dissipation).
!
! It has been moved out of domain_mod.f90, because it became complex
!
! The subroutine header has been commented out


!    ! 
!    ! Linear shallow water equations leap-frog update
!    !
!    ! Update domain%U by timestep dt, using the linear shallow water equations.
!    !
!    ! @param domain the domain to advance
!    ! @param dt the timestep to advance. Should remain constant in between repeated calls
!    !     to the function (since the numerical method assumes constant timestep)
!    !
!    subroutine one_linear_leapfrog_step(domain, dt)
!        class(domain_type), intent(inout):: domain
!        real(dp), intent(in):: dt
!
!        ! Do we represent pressure gradients with a 'truely' linear term g * depth0 * dStage/dx,
!        ! or with a nonlinear term g * depth * dStage/dx (i.e. where the 'depth' varies)?
!        logical, parameter:: truely_linear = .TRUE.

        ! Domain size
        integer(ip) :: nx, ny
        ! Useful constants for removing divisons from inner loop (e.g. dt/area)
        real(dp):: inv_cell_area_dt, inv_cell_area_dt_vh_g, inv_cell_area_dt_g
        ! d(stage)/dy; depth_{i, j+1/2}, depth_{j, i+1/2, j}
        real(dp):: dw_j(domain%nx(1)), h_jph_vec(domain%nx(1)), h_iph_vec(domain%nx(1))
        real(dp) :: h_iph_wet(domain%nx(1)), h_jph_wet(domain%nx(1))
       
        integer(ip):: j, i, xl, xu, yl, yu, n_ext, my_omp_id, n_omp_threads, loop_work_count
        integer(ip) :: yl_omp, yU_omp

#if defined(CORIOLIS) || defined(LINEAR_PLUS_NONLINEAR_FRICTION)
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

#endif

        domain%dt_last_update = dt

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
        


        !TIMER_START('LF_update')

        nx = domain%nx(1)
        ny = domain%nx(2)
      
        xL = domain%xL
        xU = domain%xU
        yL = domain%yL
        yU = domain%yU


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
        !TIMER_STOP('LF_update')

        !TIMER_START('partitioned_comms')
        if(domain%use_partitioned_comms) call domain%partitioned_comms%communicate(domain%U)
        !TIMER_STOP('partitioned_comms')

        !TIMER_START('LF_update')
        !! Boundary flux integration
        ! note : U(i, j, UH) = UH_{i+1/2, j}
        !      : U(i, j, VH) = VH_{i, j+1/2}
        n_ext = domain%exterior_cells_width
        ! Outward boundary flux over the north -- integrate over the 'interior'
        ! cells
        domain%boundary_flux_store(1) = &
            sum(domain%U((n_ext+1):(nx-n_ext),ny-n_ext,VH)) * &
            domain%distance_bottom_edge(ny-n_ext+1)
        ! Outward boundary flux over the east
        domain%boundary_flux_store(2) = &
            sum(domain%U(nx-n_ext,(1+n_ext):(ny-n_ext),UH)) * &
            domain%distance_left_edge(nx-n_ext+1)
        ! Outward boundary flux over the south
        domain%boundary_flux_store(3) = - &
            sum(domain%U((n_ext+1):(nx-n_ext),n_ext,VH)) * &
            domain%distance_bottom_edge(n_ext+1)
        ! Outward boundary flux over the west
        domain%boundary_flux_store(4) = - &
            sum(domain%U(n_ext,(1+n_ext):(ny-n_ext),UH)) * &
            domain%distance_left_edge(n_ext+1)

        domain%boundary_flux_evolve_integral = &
            sum(domain%boundary_flux_store) * dt

        !
        ! Compute fluxes relevant to the multidomain case
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
        ! Update uh, vh
        !

#if defined(LINEAR_PLUS_NONLINEAR_FRICTION)
        ! Compute some expensive parts of the friction term
        ! If domain%linear_solver_is_truely_linear, this will only happen
        ! once (on the first timestep). Otherwise it should happen every timestep
        if(.not. domain%friction_work_is_setup) call precompute_friction_work(domain)
#endif

        ! NOTE: Here we manually determine the openmp loop indices. This is
        ! required to include the Coriolis terms without having to copy
        ! memory of size domain%U(:,:,UH:VH) -- which would increase memory use
        ! by 50%
        !
        loop_work_count = yU - yL ! Number of indices between 'yU - 1' and 'yL'
        !
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, dt, nx, ny, xL, xU, yL, yU, loop_work_count)
        !
#ifdef NOOPENMP
        ! Serial loop from yL:(yU-1)
        yl_omp = yL
        yu_omp = yU - 1
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
        yl_omp = nint(loop_work_count * my_omp_id * 1.0_dp / n_omp_threads) + yL
        ! yu_omp = upper loop index of the thread my_omp_id
        yu_omp = nint(loop_work_count * (my_omp_id + 1) * 1.0_dp / n_omp_threads) + yL - 1
#endif
        !
        ! Now loop indices are determined
        !
#if defined(CORIOLIS) || defined(LINEAR_PLUS_NONLINEAR_FRICTION)
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


        !! No thread can start updating the loop until all threads have their
        !! 'bounding' coriolis terms
        !$OMP BARRIER

#endif
        !! Main update-flux loop
        !do j = 1, ny - 1
        do j = yl_omp, yu_omp
            ! For spherical coordiantes, cell area changes with y.
            ! For cartesian coordinates this could be moved out of the loop
            inv_cell_area_dt_g = gravity * dt / domain%area_cell_y(j)
            !
            ! For leap-frog with spherical coordinates, the area associated
            ! with the NS momentum term is different to the area associated
            ! with the stage and EW momentum term
            inv_cell_area_dt_vh_g = gravity * dt / &
                (HALF_dp * (domain%area_cell_y(j) + domain%area_cell_y(j+1)))
       
            ! 
            ! Try to keep control-flow and non-local memory jumps out of inner loop
            ! This improves speed on my machine with gfortran (11/08/2016)
            !
            dw_j(xL:(xU-1)) = domain%U(xL:(xU-1), j+1, STG) - domain%U(xL:(xU-1), j, STG)

            if(truely_linear) then
                !
                ! In the g * d * dStage/dx type term, let d be constant 
                !

                ! Depth at j-plus-half
                h_jph_vec(xL:xU) = merge( &
                    domain%msl_linear - HALF_dp * (domain%U(xL:xU,j+1,ELV) + domain%U(xL:xU,j,ELV)), &
                    ZERO_dp, &
                    (( domain%U(xL:xU,j+1,ELV) < -minimum_allowed_depth + domain%msl_linear).AND. &
                     ( domain%U(xL:xU,j,ELV)   < -minimum_allowed_depth + domain%msl_linear)))

                ! Depth at i-plus-half
                h_iph_vec(xL:(xU-1)) = merge( &
                    domain%msl_linear - HALF_dp * (domain%U((xL+1):xU, j, ELV) + domain%U((xL):(xU-1), j, ELV)), &
                    ZERO_dp, &
                    (( domain%U(xL:(xU-1),j,ELV) < -minimum_allowed_depth + domain%msl_linear).AND.&
                     ( domain%U((xL+1):xU,j,ELV) < -minimum_allowed_depth + domain%msl_linear)))  

                ! These variables can be used to zero UH/VH when stage < bed. However, this would introduce a nonlinearity into the
                ! equations, which seems undesirable for a 'truely-linear' approach.
                h_jph_wet(xL:xU    ) = ONE_dp
                h_iph_wet(xL:(xU-1)) = ONE_dp
            else
                !
                ! In the g * d * dStage/dx type term, let d vary. This means
                ! the equations are not actually linear!
                !

                ! Depth at j-plus-half
                h_jph_vec(xL:xU) = merge(&
                    HALF_dp * ((domain%U(xL:xU,j+1,STG) + domain%U(xL:xU,j,STG)) - &
                               (domain%U(xL:xU,j+1,ELV) + domain%U(xL:xU,j,ELV))), &
                    ZERO_dp, &
                    ((domain%U(xL:xU,j+1,STG) - domain%U(xL:xU,j+1,ELV) > minimum_allowed_depth).AND. &
                     (domain%U(xL:xU,j  ,STG) - domain%U(xL:xU,j  ,ELV) > minimum_allowed_depth)))

                ! Depth at i-plus-half
                h_iph_vec(xL:(xU-1)) = merge( &
                    HALF_dp * ((domain%U((xL+1):xU, j, STG) + domain%U(xL:(xU-1), j, STG)) -&
                               (domain%U((xL+1):xU, j, ELV) + domain%U(xL:(xU-1), j, ELV))), &
                    ZERO_dp, &
                    ((domain%U(    xL:(xU-1), j, STG) - domain%U(    xL:(xU-1),j,ELV) > minimum_allowed_depth).AND.&
                     (domain%U((xL+1):xU    , j, STG) - domain%U((xL+1):xU    ,j,ELV) > minimum_allowed_depth)))  

                ! Zero UH/VH when depths are < minimum_allowed_depth. This introduces an additional nonlinearity into the equations,
                ! but seems reasonable in the 'not-truely-linear case'
                h_jph_wet(xL:xU    ) = merge(ONE_dp, ZERO_dp, h_jph_vec(xL:xU)     > ZERO_dp)
                h_iph_wet(xL:(xU-1)) = merge(ONE_dp, ZERO_dp, h_iph_vec(xL:(xU-1)) > ZERO_dp)

            end if


#if defined(CORIOLIS) || defined(LINEAR_PLUS_NONLINEAR_FRICTION)
            ! 'Old' VH at (i+1/2, j+1/2) -- requires averaging values at 'i,j+1/2' and 'i+1, j+1/2'
            vh_iph_jph(xL:(xU-1)) = HALF_dp * &
                (domain%U(xL    :(xU-1), j, VH) + &
                 domain%U((xL+1):xU    , j, VH))

            ! 'Old' UH at (i, j+1). 
            ! First get it assuming we are not at the last loop index -- and then fix it
            ! Step1: Get everything except xL, which needs special treatment if xL == 1
            uh_i_jp1((xL+1):(xU-1)) = HALF_dp * &
                (domain%U((xL+1):(xU-1), j+1, UH) + &
                 domain%U(    xL:(xU-2), j+1, UH) )
            ! Special case of xL which is protective if xL == 1
            uh_i_jp1(xL) = HALF_dp * ( domain%U(xL          , j+1, UH) + &
                                       domain%U(max(xL-1, 1), j+1, UH) )
            ! Final step to ensure that if (j == yu_omp), then we get the 
            ! non-updated value of UH.
            uh_i_jp1(xl:(xU-1)) = merge(&
                               uh_i_jp1(xL:(xU-1)), &
                uh_i_jp1_max_loop_index(xL:(xU-1)), &
                j /= yu_omp)
#endif

#if defined(CORIOLIS)
            ! Avoid recomputing coriolis coefficient in loop. Note we set the
            ! coriolis force to zero over dry cells, as is obviously desirable.
            dt_half_coriolis(xL:(xU-1)) = dt * domain%coriolis(j) * HALF_dp * &
                merge(ONE_dp, ZERO_dp, h_iph_vec(xL:(xU-1)) > ZERO_dp)
            dt_half_coriolis_jph(xL:(xU-1)) = dt * &
                domain%coriolis_bottom_edge(j+1) * HALF_dp * &
                merge(ONE_dp, ZERO_dp, h_jph_vec(xL:(xU-1)) > ZERO_dp)
#endif

#if defined(LINEAR_PLUS_NONLINEAR_FRICTION)
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
            !    - The component 'g n^2 depth^(-7/3)' is a constant in time, because the effective depth
            !       never changes in this "linear" framework. 
            !    - The term sqrt(uh^2+vh^2) {= speed*depth} is treated explicitly.
            !    - The remaining 'uh' or 'vh' term is treated implicitly. 
            !
            ! Thus, appending the term g * depth * friction slope to the equations can be reduced to
            ! a multiplication of the form { 1/(1 + explicit_part_of_friction_terms) }
            !
            friction_multiplier_UH(xL:(xU-1)) = 1.0_dp / ( 1.0_dp + &
                dt * domain%friction_work(xL:(xU-1), j, UH) * &
                sqrt(domain%U(xL:(xU-1),j,UH)**2 + (0.5_dp * ( vh_iph_jmh(xL:(xU-1)) + vh_iph_jph(xL:(xU-1)) ) )**2 ) )

            friction_multiplier_VH(xL:(xU-1)) = 1.0_dp / (1.0_dp + &
                dt * domain%friction_work(xL:(xU-1), j, VH) * &
                sqrt(domain%U(xL:(xU-1),j,VH)**2 + (0.5_dp * ( uh_i_j(xL:(xU-1)) + uh_i_jp1(xL:(xU-1)) ) )**2 ) )

#endif


            do i = xL , (xU-1)
#ifndef CORIOLIS
                ! This update has no coriolis [other that that, it still 'works' in spherical coords]
                ! duh/dt = - g * h0/(R cos (lat)) [ d stage / dlon ]
                domain%U(i, j, UH) = domain%U(i, j, UH) * h_iph_wet(i) - &
                    inv_cell_area_dt_g * h_iph_vec(i) *&
                    (domain%U(i+1, j, STG) - domain%U(i, j, STG)) * &
                    domain%distance_left_edge(i+1)

                ! dvh/dt = - g * h0/(R) [ d stage / dlat ]
                domain%U(i, j, VH) = domain%U(i, j, VH) * h_jph_wet(i) - &
                    inv_cell_area_dt_vh_g * h_jph_vec(i) *&
                    dw_j(i) * domain%distance_bottom_edge(j+1)

#else        
                !
                ! This update has coriolis. 
                !

                ! duh/dt = - g * h0/(R cos (lat)) [ d stage / dlon ] + f*vh
                domain%U(i, j, UH) = domain%U(i, j, UH) * h_iph_wet(i) - &
                    inv_cell_area_dt_g * h_iph_vec(i) *&
                    (domain%U(i+1, j, STG) - domain%U(i, j, STG)) * &
                    domain%distance_left_edge(i+1) &
                    + dt_half_coriolis(i) * (vh_iph_jmh(i) + vh_iph_jph(i))

                ! dvh/dt = - g * h0/(R) [ d stage / dlat ] - f*uh
                domain%U(i, j, VH) = domain%U(i, j, VH) * h_jph_wet(i) - &
                    inv_cell_area_dt_vh_g * h_jph_vec(i) *&
                    dw_j(i) * domain%distance_bottom_edge(j+1) &
                    - dt_half_coriolis_jph(i) * (uh_i_j(i) + uh_i_jp1(i))

#endif

#ifdef LINEAR_PLUS_NONLINEAR_FRICTION
                ! Add the friction terms, if desired
                domain%U(i,j,UH) = friction_multiplier_UH(i) * domain%U(i,j,UH)
                domain%U(i,j,VH) = friction_multiplier_VH(i) * domain%U(i,j,VH)
#endif

            end do

#if defined(CORIOLIS) || defined(LINEAR_PLUS_NONLINEAR_FRICTION)
            ! On the next j iteration, the value of 'old' value of VH at i+1/2,
            ! j-1/2 can be derived using the current value of VH at i+1, j+1/2
            vh_iph_jmh(xL:(xU-1)) = vh_iph_jph(xL:(xU-1))
            uh_i_j(xL:(xU-1)) = uh_i_jp1(xL:(xU-1))
#endif

        end do

        !!!!!$OMP END DO
        !$OMP END PARALLEL
        domain%time = domain%time + HALF_dp*dt

        !TIMER_STOP('LF_update')
  

    contains


    ! Precompute the friction-work for the solver "leapfrog_linear_plus_nonlinear_friction"
    !
    ! This assumes stage/UH/VH/elev have been set
    !
    ! The friction work term is of the form
    !     g * n^2 / constant_depth^(7/3)
    ! The key point is that when multiplied by ||UH|| * uh, it will be equal to the
    ! standard friction form: g*constant_depth*friction_slope
    !
    ! The pre-computation removes an expensive power-law call from the inner loop
    !
    subroutine precompute_friction_work(domain)
        type(domain_type), intent(inout) :: domain

        integer(ip) :: i, j, jp1, ip1
        real(dp) :: depth_iph, depth_jph, nsq_iph, nsq_jph

        if(domain%timestepping_method /= 'leapfrog_linear_plus_nonlinear_friction') &
            stop 'precompute_friction_work can only be called with timestepping_method=leapfrog_linear_plus_nonlinear_friction'

        if(.not. allocated(domain%friction_work)) &
            stop 'friction_work is not allocated: ensure elevation is set before running this routine'

        if(truely_linear) then
            !
            ! Here we evaluate the depth assuming stage=domain%msl_linear
            ! This is the standard logic for the "truely linear" linear shallow water equations.
            !

            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, domain%nx(2)
                do i = 1, domain%nx(1)

                    ! UH component
                    ip1 = min(i+1, domain%nx(1))
                    depth_iph = 0.5_dp * (domain%msl_linear - domain%U(i,j,ELV) + domain%msl_linear - domain%U(ip1,j, ELV))
                    depth_iph = max(depth_iph, minimum_allowed_depth)
                    nsq_iph = (0.5_dp * (sqrt(domain%manning_squared(i,j)) + sqrt(domain%manning_squared(ip1,j))))**2
                    domain%friction_work(i,j,UH) = gravity * nsq_iph * (depth_iph**(-7.0_dp/3.0_dp))

                    ! VH component
                    jp1 = min(j+1, domain%nx(2))
                    depth_jph = 0.5_dp * (domain%msl_linear - domain%U(i,j,ELV) + domain%msl_linear - domain%U(i,jp1, ELV))
                    depth_jph = max(depth_jph, minimum_allowed_depth)
                    nsq_jph = (0.5_dp * (sqrt(domain%manning_squared(i,j)) + sqrt(domain%manning_squared(i,jp1))) )**2
                    domain%friction_work(i,j,VH) = gravity * nsq_jph * (depth_jph**(-7.0_dp/3.0_dp))

                end do
            end do
            !$OMP END DO
            !$OMP END PARALLEL

            ! For a truely-linear solver, this is only required once
            domain%friction_work_is_setup = .true.
        else
            !
            ! This differs from the above, in that we use domain%U(:,:,STG) to compute the depth.
            ! This means we generally need to recompute the friction-work terms at every time-step.
            !
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, domain%nx(2)
                do i = 1, domain%nx(1)

                    ! UH component
                    ip1 = min(i+1, domain%nx(1))
                    depth_iph = 0.5_dp * (domain%U(i,j,STG) - domain%U(i,j,ELV) + domain%U(ip1,j,STG) - domain%U(ip1,j, ELV))
                    depth_iph = max(depth_iph, minimum_allowed_depth)
                    nsq_iph = (0.5_dp * (sqrt(domain%manning_squared(i,j)) + sqrt(domain%manning_squared(ip1,j))))**2
                    domain%friction_work(i,j,UH) = gravity * nsq_iph * (depth_iph**(-7.0_dp/3.0_dp))

                    ! VH component
                    jp1 = min(j+1, domain%nx(2))
                    depth_jph = 0.5_dp * (domain%U(i,j,STG) - domain%U(i,j,ELV) + domain%U(i,jp1,STG) - domain%U(i,jp1, ELV))
                    depth_jph = max(depth_jph, minimum_allowed_depth)
                    nsq_jph = (0.5_dp * (sqrt(domain%manning_squared(i,j)) + sqrt(domain%manning_squared(i,jp1))) )**2
                    domain%friction_work(i,j,VH) = gravity * nsq_jph * (depth_jph**(-7.0_dp/3.0_dp))

                end do
            end do
            !$OMP END DO
            !$OMP END PARALLEL
        end if

    end subroutine

!end subroutine

