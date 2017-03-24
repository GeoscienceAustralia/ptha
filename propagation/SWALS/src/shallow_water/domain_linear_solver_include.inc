!
! The code here solves the linear shallow water equations in cartesian or spherical 
! coordinates, with or without coriolis. 
!
! It has been moved out of domain_mod.f90, because it became complex
!
! The subroutine header has been commented out


!    ! 
!    ! Linear shallow water equations leap-frog update
!    !
!    ! Update domain%U by timestep dt, using the linear shallow water equations.
!    ! Note that unlike the other timestepping routines, this does not require a
!    ! prior call to domain%compute_fluxes 
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

        integer(ip) :: nx, ny
        real(dp):: inv_cell_area_dt, inv_cell_area_dt_vh_g, inv_cell_area_dt_g
        real(dp):: dw_j(domain%nx(1)), h_jph_vec(domain%nx(1)), h_iph_vec(domain%nx(1))
        integer(ip):: j, i, xl, xu, yl, yu, n_ext, my_omp_id, n_omp_threads, loop_work_count
        integer(ip) :: yl_special, yU_special

#ifdef CORIOLIS
        ! Vector to hold the (coriolis force x time-step x 0.5), which arises
        ! in an inner loop. Use a vector so we can zero it for dry cells
        real(dp):: dt_half_coriolis(domain%nx(1)), dt_half_coriolis_jph(domain%nx(1))

        ! For the coriolis term in the UH update, we store an interpolated version of
        ! the 'pre-update' VH at 'i+1/2', 'j-1/2'. We also need the same at 'i+1/2', 'j+1/2'
        real(dp) :: vh_iph_jmh(domain%nx(1)), vh_iph_jph(domain%nx(1))

        ! For the coriolis term in the VH update, we store and interpolated version of
        ! the 'pre-update' UH at 'i', 'j'. We also need that at 'i', 'j+1'
        real(dp) :: uh_i_j(domain%nx(1)), uh_i_jp1(domain%nx(1))
        real(dp) :: uh_i_jp1_max_loop_index(domain%nx(1))
#endif

        !
        ! idea: U(i, j, UH) = UH_{i+1/2, j}
        !     : U(i, j, VH) = VH_{i, j+1/2}
        ! We rely on boundary conditions for ALL boundary stage values, and
        ! apply the boundary condition after the stage update. We then update
        ! all uh/vh values which can be updated using those stage values [so avoid
        ! having to specify boundary conditions for them]

        TIMER_START('LF_update')

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
        !do j = 2, ny-1
        do j = (yL+1),(yU-1)
        !do concurrent (j = 1:domain%nx(2))
            ! For spherical coordiantes, cell area changes with y.
            ! For cartesian coordinates this could be moved out of the loop
            inv_cell_area_dt = dt / domain%area_cell_y(j)

            !do concurrent (i = 2:(nx-1))
            do concurrent (i = (xL+1):(xU-1))
                ! dstage/dt = - 1/(R cos (lat)) [ d uh / dlon + dvhcos(lat)/dlat ]
                domain%U(i, j, STG) = domain%U(i, j, STG) - inv_cell_area_dt * &
                    ((domain%U(i, j, UH) - domain%U(i-1, j, UH))*domain%distance_left_edge(i) + &
                    (domain%U(i, j, VH)*domain%distance_bottom_edge(j+1) - &
                        domain%U(i, j-1, VH)*domain%distance_bottom_edge(j)))
        
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        
        domain%time = domain%time + dt * HALF_dp
        call domain%update_boundary()

        TIMER_START('partitioned_comms')
        call domain%partitioned_comms%communicate(domain%U)
        TIMER_STOP('partitioned_comms')

        !! Boundary flux integration
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
        ! Update uh, vh
        !
        ! NOTE: Here we manually determine the openmp loop indices. This is
        ! required to include the coriolis terms without having to copy
        ! memory of size domain%U(:,:,UH:VH) -- which would increase memory use by 50%
        !
        loop_work_count = yU - yL ! Number of indices between 'yU - 1' and 'yL'
        !
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, dt, nx, ny, xL, xU, yL, yU, loop_work_count)
        !
#ifdef NOOPENMP
        ! Serial loop from yL:(yU-1)
        yl_special = yL
        yu_special = yU - 1
#else
        ! Parallel loop from yL:(yU-1)
        !
        ! NOTE: In fortran, 
        !     DO i = start, end
        !       .... code here ...
        !     END DO
        ! will not do anything at all if start > end [that's the standard, if we don't
        !   specify a stride other than 1].
        ! This behaviour is important for ensuring the openmp loop sharing
        !   code here works, even if e.g. we have more omp threads than loop
        !   indices. In that case, we will have some threads with yl_special >
        !   yu_special -- which will not loop -- and that is the desired behaviour.
        !
        my_omp_id = omp_get_thread_num()
        n_omp_threads = omp_get_num_threads()
        ! yl_special = lower loop index of the thread my_omp_id
        yl_special = nint(loop_work_count * my_omp_id * 1.0_dp / n_omp_threads) + yL
        ! yu_special = upper loop index of the thread my_omp_id
        yu_special = nint(loop_work_count * (my_omp_id + 1) * 1.0_dp / n_omp_threads) + yL - 1
#endif
        !
        ! Now loop indices are determined
        !
#ifdef CORIOLIS
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
        ! Before starting the loop, get the value of VH at 'i+1/2', 'yl_special-1/2', to prevent the
        ! possibility that it is updated by another openmp thread before we read it
        !
        if(yl_special == 1) then
            ! In this case, there is no 'j-1/2' value, just use the boundary value
            vh_iph_jmh(xL:(xU-1)) = HALF_dp * (domain%U(xL     :(xU-1), 1, VH) + &
                                                domain%U((xL+1):xU    , 1, VH))
        else
            ! Average values at 'i' and 'i+1'
            vh_iph_jmh(xL:(xU-1)) = HALF_dp * (domain%U(xL    :(xU-1), yl_special - 1, VH) + &
                                               domain%U((xL+1):xU    , yl_special - 1, VH))
        end if

        !
        ! Before starting the loop, store the value of UH at 'i', 'yu_special+1', to prevent the
        ! possibility that it is updated by another openmp thread before we read it
        !
        if(xL > 1) then
            ! Average values at 'i-1/2' and 'i+1/2'
            uh_i_jp1_max_loop_index(xL:(xU-1)) = HALF_dp * (domain%U( (xL-1):(xU-2), yu_special + 1, UH) + &
                                                            domain%U( xL    :(xU-1), yu_special + 1, UH))
        else
            ! Avoid out-of bounds error
            uh_i_jp1_max_loop_index((xL+1):(xU-1)) = HALF_dp * (domain%U( xL   :(xU-2), yu_special + 1, UH) + &
                                                                domain%U( xL+1 :(xU-1), yu_special + 1, UH))
            uh_i_jp1_max_loop_index(xL) = domain%U(1, yu_special + 1, UH)
        end if

        !
        ! Get the initial value of UH at 'i, yl_special'. This 'starting value' is required because 
        ! inside the loop we set 'uh_i_j = uh_i_jph' after the update.
        !
        uh_i_j((xL+1):(xU-1)) = HALF_dp * ( domain%U((xL+1):(xU-1), yl_special, UH) + &
                                            domain%U(    xL:(xU-2), yl_special, UH) )
        !
        ! Special case of xL which is protective if xL == 1
        uh_i_j(xL) = HALF_dp * ( domain%U(xL, yl_special, UH) + domain%U(max(xL-1, 1), yl_special, UH) )


        !! No thread can start updating the loop until all threads have their
        !! 'bounding' coriolis terms
        !$OMP BARRIER

#endif
        !! Main update-flux loop
        !do j = 1, ny - 1
        do j = yl_special, yu_special
            ! For spherical coordiantes, cell area changes with y.
            ! For cartesian coordinates this could be moved out of the loop
            ! For leap-frog, the area associated with the NS momentum term is different
            inv_cell_area_dt_g = gravity * dt / domain%area_cell_y(j)
            inv_cell_area_dt_vh_g = gravity * dt / (HALF_dp * (domain%area_cell_y(j) + domain%area_cell_y(j+1)))
       
            ! 
            ! Try to keep control-flow and non-local memory jumps out of inner loop
            ! This improves speed on my machine with gfortran (11/08/2016)
            !
            dw_j = domain%U(:, j+1, STG) - domain%U(:, j, STG)

            if(truely_linear) then
                !
                ! In the g * d * dStage/dx type term, let d be constant 
                !

                ! Depth at j-plus-half
                h_jph_vec(xL:xU) = merge(domain%msl_linear - HALF_dp * (domain%U(xL:xU,j+1,ELV) + domain%U(xL:xU,j,ELV)), &
                    ZERO_dp, &
                    (( domain%U(xL:xU,j+1,ELV) < -minimum_allowed_depth + domain%msl_linear).AND. &
                     ( domain%U(xL:xU,j,ELV)   < -minimum_allowed_depth + domain%msl_linear)))

                ! Depth at i-plus-half
                h_iph_vec(xL:(xU-1)) = merge( &
                    domain%msl_linear - HALF_dp * (domain%U((xL+1):xU, j, ELV) + domain%U((xL):(xU-1), j, ELV)), &
                    ZERO_dp, &
                    (( domain%U(xL:(xU-1),j,ELV) < -minimum_allowed_depth + domain%msl_linear).AND.&
                     ( domain%U((xL+1):xU,j,ELV) < -minimum_allowed_depth + domain%msl_linear)))  
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
                     (domain%U(xL:xU,j,STG) - domain%U(xL:xU,j,ELV) > minimum_allowed_depth)))

                ! Depth at i-plus-half
                h_iph_vec(xL:(xU-1)) = merge( &
                    HALF_dp * ((domain%U((xL+1):xU, j, STG)  + domain%U(xL:(xU-1), j, STG)) -&
                               (domain%U((xL+1):xU, j, ELV) - domain%U(xL:(xU-1), j, ELV))), &
                    ZERO_dp, &
                    ((domain%U(xL:(xU-1), j, STG) - domain%U(xL:(xU-1),j,ELV) > minimum_allowed_depth).AND.&
                        (domain%U((xL+1):xU, j, STG) - domain%U((xL+1):xU,j,ELV) > minimum_allowed_depth)))  


            end if

#ifdef CORIOLIS
            ! 'Old' VH at (i+1/2, j+1/2) -- requires averaging values at 'i,j+1/2' and 'i+1, j+1/2'
            vh_iph_jph(xL:(xU-1)) = HALF_dp * (domain%U(xL    :(xU-1), j, VH) + &
                                               domain%U((xL+1):xU    , j, VH))

            ! 'Old' UH at (i, j+1). 
            ! First get it assuming we are not at the last loop index -- and then fix it
            ! Step1: Get everything except xL, which needs special treatment if xL == 1
            uh_i_jp1((xL+1):(xU-1)) = HALF_dp * (domain%U((xL+1):(xU-1), j+1, UH) + &
                                                 domain%U(    xL:(xU-2), j+1, UH) )
            ! Special case of xL which is protective if xL == 1
            uh_i_jp1(xL) = HALF_dp * ( domain%U(xL          , j+1, UH) + &
                                       domain%U(max(xL-1, 1), j+1, UH) )
            ! Final step to ensure that if (j == yu_special), then we get the 
            ! non-updated value of UH.
            uh_i_jp1(xl:(xU-1)) = merge(&
                               uh_i_jp1(xL:(xU-1)), &
                uh_i_jp1_max_loop_index(xL:(xU-1)), &
                j /= yu_special)

            ! Avoid recomputing coriolis coefficient in loop. Note we set the
            ! coriolis force to zero over dry cells, as is obviously desirable.
            dt_half_coriolis(xL:(xU-1)) = dt * domain%coriolis(j) * HALF_dp * &
                merge(ONE_dp, ZERO_dp, h_iph_vec(xL:(xU-1)) > ZERO_dp)
            dt_half_coriolis_jph(xL:(xU-1)) = dt * domain%coriolis_bottom_edge(j+1) * HALF_dp * &
                merge(ONE_dp, ZERO_dp, h_jph_vec(xL:(xU-1)) > ZERO_dp)

            !print*, ''
            !print*, 'j: ', j, xL, xU, yl_special, yu_special
            !print*, 'uh_i_j: ', minval(uh_i_j), maxval(uh_i_j)
            !print*, 'uh_i_jp1: ', minval(uh_i_jp1), maxval(uh_i_jp1), maxval(abs(domain%U(:,j,UH)))
            !print*, 'vh_iph_jmh: ', minval(vh_iph_jmh), maxval(vh_iph_jmh), maxval(abs(domain%U(:,j,VH)))
            !print*, 'vh_iph_jph: ', minval(vh_iph_jph), maxval(vh_iph_jph)
#endif



            !do concurrent (i = 1:(nx-1))
            do concurrent (i = xL:(xU-1))
#ifndef CORIOLIS
                ! This update has no coriolis [other that that, it still 'works' in spherical coords]
                ! duh/dt = - g * h0/(R cos (lat)) [ d stage / dlon ]
                domain%U(i, j, UH) = domain%U(i, j, UH) - inv_cell_area_dt_g * h_iph_vec(i) *&
                    (domain%U(i+1, j, STG) - domain%U(i, j, STG)) * domain%distance_left_edge(i+1)

                ! dvh/dt = - g * h0/(R) [ d stage / dlat ]
                domain%U(i, j, VH) = domain%U(i, j, VH) - inv_cell_area_dt_vh_g * h_jph_vec(i) *&
                    dw_j(i) * domain%distance_bottom_edge(j+1)
#else        
                !
                ! This update has coriolis. 
                !

                ! duh/dt = - g * h0/(R cos (lat)) [ d stage / dlon ] + f*vh
                domain%U(i, j, UH) = domain%U(i, j, UH) - inv_cell_area_dt_g * h_iph_vec(i) *&
                    (domain%U(i+1, j, STG) - domain%U(i, j, STG)) * domain%distance_left_edge(i+1) &
                    + dt_half_coriolis(i) * (vh_iph_jmh(i) + vh_iph_jph(i))

                ! dvh/dt = - g * h0/(R) [ d stage / dlat ] - f*uh
                domain%U(i, j, VH) = domain%U(i, j, VH) - inv_cell_area_dt_vh_g * h_jph_vec(i) *&
                    dw_j(i) * domain%distance_bottom_edge(j+1) &
                    - dt_half_coriolis_jph(i) * (uh_i_j(i) + uh_i_jp1(i))

#endif
            end do

#ifdef CORIOLIS
            ! On the next j iteration, the value of 'old' value of VH at i+1/2, j-1/2 can be derived
            ! using the current value of VH at i+1, j+1/2
            vh_iph_jmh(xL:(xU-1)) = vh_iph_jph(xL:(xU-1))
            uh_i_j(xL:(xU-1)) = uh_i_jp1(xL:(xU-1))
#endif

        end do

        !!!!!$OMP END DO
        !$OMP END PARALLEL
        domain%time = domain%time + HALF_dp*dt

        TIMER_STOP('LF_update')
  
!    end subroutine
