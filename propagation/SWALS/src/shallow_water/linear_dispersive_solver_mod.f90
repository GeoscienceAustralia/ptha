! Compile with -DEVOLVE_TIMER to add detailed timing of the evolve loop
#ifdef EVOLVE_TIMER
#   define EVOLVE_TIMER_START(tname) call ds%evolve_timer%timer_start(tname)
#   define EVOLVE_TIMER_STOP(tname)  call ds%evolve_timer%timer_end(tname)
#else
#   define EVOLVE_TIMER_START(tname)
#   define EVOLVE_TIMER_STOP(tname)
#endif

! Use definition below to use Peregrine dispersion which includes an additional topographic term
! (by default we ignore this term, similar to JAGURS and many other tsunami codes)
!#define DISPERSIVE_PEREGRINE_IN_FLUX_FORM

! FIXME: The tridiagonal solve method might have efficiency improved by precomputing
!     lower_uh, diag_uh, upper_uh, lower_vh, diag_vh, upper_vh
!     .... and related infor for RHS vh/uh terms....
! at the cost of more storage
!     Or instead of dgtsv we might use a mixture of dgttrf (factorize once, store the factors and pivots) and dgttrs (solve)

module linear_dispersive_solver_mod
    !!
    !! Type for solving the linear dispersive equation as per:
    !!   Baba, T.; Takahashi, N.; Kaneda, Y.; Ando, K.; Matsuoka, D. & Kato, T. Parallel Implementation of Dispersive Tsunami Wave
    !!       Modeling with a Nesting Algorithm for the 2011 Tohoku Tsunami Pure and Applied Geophysics, 2015, 172, 3455-3472
    !! Or
    !!   Baba, T.; Allgeyer, S.; Hossen, J.; Cummins, P. R.; Tsushima, H.; Imai, K.; Yamashita, K. & Kato, T. Accurate numerical
    !!   simulation of the far-field tsunami caused by the 2011 Tohoku earthquake, including the effects of Boussinesq dispersion,
    !!   seawater density stratification, elastic loading, and gravitational potential change Ocean Modelling, Elsevier BV, 2017, 111,
    !!   46–54
    !! The equations are:
    !!     d( uh ) / dt + {...nonlinear-shallow-water-terms...} = &
    !!        h0^2 / ( 3 R_coslat) * d/dlon [ 1/R_coslat ( d^2( uh ) / (dt dlon) + d^2(vh coslat) / (dt dlat) ) ]
    !!     d (vh ) / dt + {...nonlinear-shallow-water-terms...} = &
    !!        h0^2 / ( 3 R       ) * d/dlat [ 1/R_coslat ( d^2( uh ) / (dt dlon) + d^2(vh coslat) / (dt dlat) ) ]
    !!
    !! Where 'h0' is the depth below MSL_linear. Note that Baba et al use 'co-latitude' instead of latitude, so my
    !! cos(theta) is their sin(theta)
    !!
    !! In tensor notation (with time as a parameter), if writing the NLSW in contravariant form, this would be written as 
    !!     g^{lj} \nabla_j( \frac{\partial (\nabla_i q^{i})}{\partial t} )
    !! where 
    !!     q^i = hu^i = contravariant components of flux = [ uh/(Rcoslat), vh/R ]  and 
    !!     \nabla_{j} is the covariant derivative and 
    !!     g^{lj} is the metric tensor, used to raise an index (so the term is a contravariant tensor). 
    !! (Conversion to the conventional form in spherical coordinates would also involve multiplying the UH equation
    !! by Rcoslat and the VH equation by R). Notice the term inside the \nabla_j covariant derivative
    !! is a scalar (so this covariant derivative reduces to the partial derivative) and that scalar is the time derivative of the 
    !! mass flux terms (as appears in the mass conservation equation). So there are no other "hidden" terms caused by
    !! Christoffel symbols.
    !!
    !! ** Notice we can bring the time-derivative out to the front. This is key to solving the equation numerically. **
    !!     d( uh ) / dt + {...nonlinear-shallow-water-terms...} = &
    !!        d/dt * { h0^2 / ( 3 R_coslat) * d/dlon [ 1/R_coslat ( d( uh ) / (dlon) + d(vh coslat) / (dlat) ) ] }
    !!     d (vh ) / dt + {...nonlinear-shallow-water-terms...} = &
    !!        d/dt * { h0^2 / ( 3 R       ) * d/dlat [ 1/R_coslat ( d( uh ) / (dlon) + d(vh coslat) / (dlat) ) ] }
    !!
    !! Another dispersive model supported here is is Peregrine dispersion, which has the form
    !!     d( uh ) / dt + {...nonlinear-shallow-water-terms...} = &
    !!        h0^2 / ( 2 R_coslat) * d/dlon [ 1/R_coslat ( d^2( uh ) / (dt dlon) + d^2(vh coslat) / (dt dlat) ) ] &
    !!       -h0^3 / ( 6 R_coslat) * d/dlon [ 1/R_coslat ( d^2( u  ) / (dt dlon) + d^2(v  coslat) / (dt dlat) ) ]
    !!     d (vh ) / dt + {...nonlinear-shallow-water-terms...} = &
    !!        h0^2 / ( 2 R       ) * d/dlat [ 1/R_coslat ( d^2( uh ) / (dt dlon) + d^2(vh coslat) / (dt dlat) ) ] &
    !!       -h0^3 / ( 6 R       ) * d/dlat [ 1/R_coslat ( d^2( u  ) / (dt dlon) + d^2(v  coslat) / (dt dlat) ) ]
    !! and reduces to the other dispersive model with flat topography.
    !!
    !! To implement this in both cartesian and spherical coordinates, note that
    !!     radius_earth = distance_cell_left_edge / dy_in_radians ! This is independent of the latitude or longitude.
    !! and
    !!     R_coslat = distance_cell_bottom_edge / dx_in_radians ! Need to compute the cell-bottom distance at the right 'latitude'.
    !! So
    !!     coslat = (distance_cell_bottom_edge / dx_in_radians) / ( distance_cell_left_edge / dy_in_radians )
    !!     ! Note the factor to convert from degrees to radians will drop out in this equation.
    !! Also
    !!     R_coslat_dlon = distance_cell_bottom_edge ! Here we need to compute the cell-bottom distance at the right 'latitude'.
    !! and
    !!     R_dlat = distance_cell_left_edge ! This is independent of the latitude or longitude
    !! and
    !!     R_dlat_coslat = distance_cell_left_edge * coslat = ( distance_cell_bottom_edge / dx_in_radians) * dy_in_radians
    !!

    use global_mod, only: dp, ip
    use logging_mod, only: log_output_unit
    use quadratic_extrapolation_mod, only: quadratic_extrapolation_type
    use stop_mod, only: generic_stop
    use timer_mod, only: timer_type

    implicit none

    integer, parameter :: STG=1, UH=2, VH=3, ELV=4

    type dispersive_solver_type

        logical :: is_staggered_grid 
        real(dp), allocatable :: RHS(:,:,:), last_U(:,:,:)
            !! Right hand side, and value of U at the start of the (dispersive) timestep

        real(dp), allocatable :: Ax(:,:,:)
            !! For tridiagonal iteration, work array used to compute dispersive RHS terms

        real(dp) :: td1 = 0.0_dp, td2 = -1.0_dp 
            !! Depths (below MSL) used to linearly taper-off dispersion. 
            !! Set td1 > td2 to linearly reduce dispersive terms to zero between depths of td1 --> td2
            !! Default corresponds to no tapering. 

        integer(ip) :: tridiagonal_inner_iter = 2_ip
            !! The tridiagonal solver does this many iterations of 'x-sweep'+'ysweep' for each call to solve_tridiag

        type(quadratic_extrapolation_type) :: qet
            !! Allow quadratic extrapolation in time

        type(timer_type) :: evolve_timer
            !! For timing dispersive terms only

        !logical :: debug_me = .false. ! Hack to help debugging

        contains
        procedure, non_overridable :: store_last_U => store_last_U
        procedure, non_overridable :: setup => setup_dispersive_solver

        ! Tridiagonal solver
        procedure, non_overridable :: solve_tridiag => linear_dispersive_solve_TRIDIAG

    end type

    contains

    subroutine store_last_U(ds, U)
        ! For timestepping we need to retain U prior to the shallow water timestep.
        class(dispersive_solver_type), intent(inout) :: ds
        real(dp), intent(in) :: U(:,:,:)

        integer(ip) :: j

EVOLVE_TIMER_START('disp_store_last_U')
        !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(ds, U)
        do j = 1, size(U,2)
            ds%last_U(:,j,UH:ELV) = U(:,j,UH:ELV)
        end do
        !$OMP END PARALLEL DO
EVOLVE_TIMER_STOP('disp_store_last_U')

    end subroutine

    subroutine setup_dispersive_solver(ds, nx, ny, is_staggered_grid)
        !! Allocate workspace
        class(dispersive_solver_type), intent(inout) :: ds
        integer(ip), intent(in) :: nx, ny
        logical, intent(in) :: is_staggered_grid

        integer(ip) :: i, j

EVOLVE_TIMER_START('disp_setup')
        ds%is_staggered_grid = is_staggered_grid

        if(allocated(ds%RHS)) deallocate(ds%RHS)
        if(allocated(ds%last_U)) deallocate(ds%last_U)
        allocate(ds%RHS(nx, ny, UH:VH), ds%last_U(nx, ny, UH:ELV))
        ! OMP first-touch allocation
        !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(ds)
        do j = 1, ny
            ds%RHS(:,j,UH:VH) = 0.0_dp
            ds%last_U(:,j,UH:ELV) = 0.0_dp
        end do
        !$OMP END PARALLEL DO

        ! Workspace for tridiagonal iteration
        if(allocated(ds%Ax)) deallocate(ds%Ax)
        allocate(ds%Ax(nx,ny,UH:VH))
        !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(ds)
        do j = 1, ny
            ds%Ax(:,j,UH:VH) = 0.0_dp
        end do
        !$OMP END PARALLEL DO

        ! Prepare to extrapolate UH/VH in time
        call ds%qet%setup(nx, ny, 2)
EVOLVE_TIMER_STOP('disp_setup')

    end subroutine

    
    subroutine linear_dispersive_staggered_matmult(elev, uh, vh, msl_linear, distance_bottom_edge, distance_left_edge, &
            dlon, dlat, td1, td2, solution_uh, solution_vh, update_UH, update_VH)
        !
        ! Towards solving the linear dispersive equation as per:
        !   Baba, T.; Takahashi, N.; Kaneda, Y.; Ando, K.; Matsuoka, D. & Kato, T. Parallel Implementation of Dispersive Tsunami Wave
        !       Modeling with a Nesting Algorithm for the 2011 Tohoku Tsunami Pure and Applied Geophysics, 2015, 172, 3455-3472
        ! Or
        !   Baba, T.; Allgeyer, S.; Hossen, J.; Cummins, P. R.; Tsushima, H.; Imai, K.; Yamashita, K. & Kato, T. Accurate numerical
        !   simulation of the far-field tsunami caused by the 2011 Tohoku earthquake, including the effects of Boussinesq dispersion,
        !   seawater density stratification, elastic loading, and gravitational potential change Ocean Modelling, Elsevier BV, 2017,
        !   111, 46–54
        !
        ! BACKGROUND:
        !
        ! The shallow water equation for depth-integrated-x-velocity, with a dispersive term, can be organised like:
        !    d(uh)/dt + { other shallow water equations terms } = dispersive_term_uh
        ! where dispersive_term_uh contains a time-derivative of some second-order spatial derivatives, see further explanation
        ! in the comments at the start of this file.
        !
        ! This can be written as
        !     uh  = [ shallow_water_solution ] + dt*dispersive_term_uh
        ! Since dispersive_term_uh contains a single derivative wrt time, it turns out that (dt * dispersive_term_uh)
        ! does not contain the term 'dt' (although it depends on the flow variables at the last timestep)
        !
        ! Suppose that ( dt * dispersive_term_uh) is written as:
        !    (dt * dispersive_term_uh) = { A_uh%*%x - A_uh%*%x_lasttimestep }
        ! , where 'A_uh' is a sparse matrix that discretizes the spatial derivatives in the dispersive term, 'x' is the flow
        ! variables [stage, uh, vh, elevation] at the next time-step, and %*% is matrix multiplication.
        !
        ! Then this routine computes A_uh%*%x
        !
        ! The above is repeated for vh (with somewhat different equations in spherical coordinates).
        !
        real(dp), intent(in) :: elev(:,:), uh(:,:), vh(:,:), &
            dlon, dlat, msl_linear, td1, td2, &
            distance_bottom_edge(:), distance_left_edge(:)
        real(dp), intent(inout) :: solution_uh(:,:), solution_vh(:,:)
        logical, intent(in) :: update_UH, update_VH

        integer(ip) :: i, j, ip2, jp2
        real(dp) :: R_coslat_dlat, R_coslat_dlon, coslat_jph, coslat_jmh, d_iph_j, d_i_jph, dispersive_premult
        real(dp) :: r_dlat, R_coslat_dlon_jp1, R_coslat_dlon_j, R_coslat_dlat_jp1, R_coslat_dlat_j, coslat_jp1h
        real(dp) :: dinv_ip1h_j, dinv_iph_j, dinv_imh_j, &
                    dinv_ip1_jph, dinv_ip1_jmh, &
                    dinv_i_jph, dinv_i_jmh, &
                    dinv_iph_jp1, dinv_imh_jp1, &
                    dinv_i_jp1h

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(uh, vh, elev, msl_linear, solution_uh, solution_vh, &
        !$OMP                                  dlat, dlon, td1, td2, &
        !$OMP                                  distance_left_edge, distance_bottom_edge, update_UH, update_VH)

        if(update_UH) then
            !$OMP DO
            do j = 2, size(uh, 2)-1

                !if(j == 1 .or. j == size(uh, 2)) then
                !    solution_uh(:,j) = uh(:,j)
                !    cycle
                !end if

                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
                ! cartesian coordinates.
                R_coslat_dlon   = 0.5_dp * (distance_bottom_edge(j+1) + distance_bottom_edge(j  )) ! x-cell-distance at uh(i,j)
                R_coslat_dlat   = R_coslat_dlon * dlat / dlon
                coslat_jph = distance_bottom_edge(j+1) * dlat / (dlon * distance_left_edge(1))
                coslat_jmh = distance_bottom_edge(j+0) * dlat / (dlon * distance_left_edge(1))

                !solution_uh(1,j) = uh(1,j) ! Instead this is set in the caller routine
                !solution_uh(size(uh,1),j) = uh(size(uh,1),j)
                !$OMP SIMD
                do i = 2, size(uh, 1) - 1
                    ! UH dispersive term

                    ! Depth at location of uh(i,j) -- or zero if either neighbour is dry at MSL
                    d_iph_j = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i+1,j)), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i+1, j) < msl_linear)
#ifdef DISPERSIVE_PEREGRINE_IN_FLUX_FORM
                    !
                    ! Peregrine dispersion
                    !
                    dispersive_premult = d_iph_j*d_iph_j / (R_coslat_dlon) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_iph_j - td2)/(td1 - td2)))

                    ! Compute various inverse depth terms, which are zeroed at 'dry (at MSL)' cells.
                    ! FIXME: Could avoid recomputing some of these (noting the loop order)
                    ip2 = min(i+2, size(uh, 1)) ! Safe i+2 index
                    dinv_ip1h_j = merge(1.0_dp/(0.5_dp * (msl_linear - elev(i+1,j) + msl_linear - elev(ip2,j))), 0.0_dp, &
                        elev(i+1,j) < msl_linear .and. elev(ip2, j) < msl_linear)
                    dinv_iph_j = merge(1.0_dp/d_iph_j, 0.0_dp, d_iph_j > 0.0_dp)
                    dinv_imh_j = merge(1.0_dp/(0.5_dp * (msl_linear - elev(i-1,j) + msl_linear - elev(i,j))), 0.0_dp, &
                        elev(i-1,j) < msl_linear .and. elev(i, j) < msl_linear)
                    dinv_ip1_jph = merge(1.0_dp/(0.5_dp * (msl_linear - elev(i+1,j+1) + msl_linear - elev(i+1,j))), 0.0_dp, &
                        elev(i+1,j+1) < msl_linear .and. elev(i+1, j) < msl_linear)
                    dinv_ip1_jmh = merge(1.0_dp/(0.5_dp * (msl_linear - elev(i+1,j) + msl_linear - elev(i+1,j-1))), 0.0_dp, &
                        elev(i+1,j) < msl_linear .and. elev(i+1, j-1) < msl_linear)
                    dinv_i_jph = merge(1.0_dp/(0.5_dp * (msl_linear - elev(i,j+1) + msl_linear - elev(i,j))), 0.0_dp, &
                        elev(i,j+1) < msl_linear .and. elev(i, j) < msl_linear)
                    dinv_i_jmh = merge(1.0_dp/(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i,j-1))), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i, j-1) < msl_linear)

                    solution_uh(i,j) = &
                        ! Include all terms in A_uh%*%x here
                        dispersive_premult * 0.5_dp * ( &
                            ! lon-diff{ 1/(R_coslat) * d/dlon(uh) }
                            1.0_dp / R_coslat_dlon * ( uh(i+1,j) - 2.0_dp * uh(i,j) + uh(i-1,j) ) + &
                            ! lon-diff{ 1/(R_coslat) * d/dlat( coslat * vh ) }
                            1.0_dp / R_coslat_dlat * ( (coslat_jph * vh(i+1,j) - coslat_jmh * vh(i+1,j-1)) - &
                                                       (coslat_jph * vh(i  ,j) - coslat_jmh * vh(i  ,j-1)) ) ) &
                        - dispersive_premult * (1.0_dp/6.0_dp) * d_iph_j * ( &
                            ! lon-diff{ 1/(R_coslat) * d/dlon(u) }
                            1.0_dp / R_coslat_dlon*(uh(i+1,j)*dinv_ip1h_j - 2.0_dp*uh(i,j)*dinv_iph_j + uh(i-1,j)*dinv_imh_j ) + &
                            ! lon-diff{ 1/(R_coslat) * d/dlat( coslat * v ) }
                            1.0_dp / R_coslat_dlat*((coslat_jph*vh(i+1,j)*dinv_ip1_jph - coslat_jmh*vh(i+1,j-1)*dinv_ip1_jmh) - &
                                                    (coslat_jph*vh(i  ,j)*dinv_i_jph   - coslat_jmh*vh(i  ,j-1)*dinv_i_jmh )) )
#else
                    !
                    ! Simplest dispersive model, matching JAGURS
                    !
                    dispersive_premult = d_iph_j*d_iph_j / (3.0_dp * R_coslat_dlon) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_iph_j - td2)/(td1 - td2)))

                    solution_uh(i,j) = &
                        ! Include all terms in A_uh%*%x here
                        dispersive_premult * ( &
                            ! lon-diff{ 1/(R_coslat) * d/dlon(uh) }
                            1.0_dp / R_coslat_dlon * ( uh(i+1,j) - 2.0_dp * uh(i,j) + uh(i-1,j) ) + &
                            ! lon-diff{ 1/(R_coslat) * d/dlat( coslat * vh ) }
                            1.0_dp / R_coslat_dlat * ( (coslat_jph * vh(i+1,j) - coslat_jmh * vh(i+1,j-1)) - &
                                                       (coslat_jph * vh(i  ,j) - coslat_jmh * vh(i  ,j-1)) ) )
#endif
                end do
            end do
            !$OMP END DO NOWAIT
        end if

        if(update_VH) then
            !$OMP DO
            do j = 2, size(vh, 2)-1

                !if(j == 1 .or. j == size(vh, 2)) then
                !    solution_vh(:,j) = vh(:,j)
                !    cycle
                !end if

                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
                ! cartesian coordinates.
                r_dlat = distance_left_edge(1)
                R_coslat_dlon_jp1 = 0.5_dp * ( distance_bottom_edge(j+2) + distance_bottom_edge(j+1) )
                R_coslat_dlon_j   = 0.5_dp * ( distance_bottom_edge(j+1) + distance_bottom_edge(j+0) )
                R_coslat_dlat_jp1 = R_coslat_dlon_jp1 * dlat / dlon
                R_coslat_dlat_j   = R_coslat_dlon_j   * dlat / dlon

                coslat_jp1h = distance_bottom_edge(j+2) * dlat / ( distance_left_edge(1)  * dlon )
                coslat_jph  = distance_bottom_edge(j+1) * dlat / ( distance_left_edge(1)  * dlon )
                coslat_jmh  = distance_bottom_edge(j+0) * dlat / ( distance_left_edge(1)  * dlon )

                !solution_vh(1,j) = vh(1,j) ! Instead this is set in caller routine
                !solution_vh(size(vh,1), j) = vh(size(vh,1), j)
                !$OMP SIMD
                do i = 2, size(vh, 1) - 1

                    ! Depth at location of vh(i,j) -- or zero if either neighbour is dry (at MSL)
                    d_i_jph = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i,j+1)), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i, j+1) < msl_linear)
#ifdef DISPERSIVE_PEREGRINE_IN_FLUX_FORM
                    !
                    ! Peregrine dispersion
                    !
                    dispersive_premult = d_i_jph*d_i_jph / (r_dlat) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_jph - td2)/(td1 - td2)))

                    ! Compute various inverse depth terms, which are zeroed at 'dry (at MSL)' cells.
                    ! FIXME: Could avoid recomputing some of these (noting the loop order)
                    dinv_iph_jp1 = merge(1.0_dp / (0.5_dp * (msl_linear - elev(i,j+1) + msl_linear - elev(i+1,j+1))), 0.0_dp, &
                        elev(i,j+1) < msl_linear .and. elev(i+1, j+1) < msl_linear)
                    dinv_imh_jp1 = merge(1.0_dp / (0.5_dp * (msl_linear - elev(i,j+1) + msl_linear - elev(i-1,j+1))), 0.0_dp, &
                        elev(i,j+1) < msl_linear .and. elev(i-1, j+1) < msl_linear)
                    dinv_iph_j = merge(1.0_dp / (0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i+1,j))), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i+1, j) < msl_linear)
                    dinv_imh_j = merge(1.0_dp / (0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i-1,j))), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i-1, j) < msl_linear)
                    jp2 = min(j+2, size(vh, 2)) ! Safe j+2
                    dinv_i_jp1h = merge(1.0_dp / (0.5_dp * (msl_linear - elev(i,j+1) + msl_linear - elev(i,jp2))), 0.0_dp, &
                        elev(i,j+1) < msl_linear .and. elev(i, jp2) < msl_linear)
                    dinv_i_jph = merge(1.0_dp/d_i_jph, 0.0_dp, d_i_jph > 0.0_dp)
                    dinv_i_jmh = merge(1.0_dp / (0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i,j-1))), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i, j-1) < msl_linear)

                    solution_vh(i,j) = &
                        ! Include all terms in A_vh%*%x here
                        dispersive_premult * 0.5_dp * ( &
                            ! lat-diff{ 1/R_coslat duh/dlon }
                            ( 1.0_dp / R_coslat_dlon_jp1 * (uh(i  , j+1) - uh(i-1, j+1)) - &
                              1.0_dp / R_coslat_dlon_j   * (uh(i  , j  ) - uh(i-1, j  )) ) + &
                            ! lat-diff{ 1/R_coslat d(coslat vh)/dlat }
                            ( 1.0_dp / R_coslat_dlat_jp1 * ( coslat_jp1h * vh(i, j+1) - coslat_jph * vh(i, j  )) - &
                              1.0_dp / R_coslat_dlat_j   * ( coslat_jph  * vh(i, j  ) - coslat_jmh * vh(i, j-1)) ) ) & 
                        - dispersive_premult * (1.0_dp/6.0_dp) * d_i_jph * ( &
                            ! lat-diff{ 1/R_coslat du/dlon }
                            ( 1.0_dp / R_coslat_dlon_jp1 * (uh(i  , j+1)*dinv_iph_jp1 - uh(i-1, j+1)*dinv_imh_jp1) - &
                              1.0_dp / R_coslat_dlon_j   * (uh(i  , j  )*dinv_iph_j   - uh(i-1, j  )*dinv_imh_j  )  ) + &
                            ! lat-diff{ 1/R_coslat d(coslat v)/dlat }
                            ( 1.0_dp / R_coslat_dlat_jp1 *(coslat_jp1h*vh(i, j+1)*dinv_i_jp1h-coslat_jph*vh(i, j  )*dinv_i_jph) - &
                              1.0_dp / R_coslat_dlat_j   *(coslat_jph *vh(i, j  )*dinv_i_jph -coslat_jmh*vh(i, j-1)*dinv_i_jmh) ) )


#else
                    !
                    ! Simplest dispersive model, matching JAGURS
                    !
                    dispersive_premult = d_i_jph*d_i_jph / (3.0_dp * r_dlat) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_jph - td2)/(td1 - td2)))

                    solution_vh(i,j) = &
                        ! Include all terms in A_vh%*%x here
                        dispersive_premult * ( &
                            ! lat-diff{ 1/R_coslat duh/dlon }
                            ( 1.0_dp / R_coslat_dlon_jp1 * (uh(i  , j+1) - uh(i-1, j+1)) - &
                              1.0_dp / R_coslat_dlon_j   * (uh(i  , j  ) - uh(i-1, j  )) ) + &
                            ! lat-diff{ 1/R_coslat d(coslat vh)/dlat }
                            ( 1.0_dp / R_coslat_dlat_jp1 * ( coslat_jp1h * vh(i, j+1) - coslat_jph * vh(i, j  )) - &
                              1.0_dp / R_coslat_dlat_j   * ( coslat_jph  * vh(i, j  ) - coslat_jmh * vh(i, j-1)) ) )
#endif
                end do
            end do
            !$OMP END DO
        end if

        !$OMP END PARALLEL

    end subroutine

    subroutine linear_dispersive_cellcentred_matmult(elev, uh, vh, msl_linear, distance_bottom_edge, distance_left_edge, &
            dlon, dlat, td1, td2, solution_uh, solution_vh, update_UH, update_VH)
        !
        ! Towards solving the linear dispersive equation as per:
        !   Baba, T.; Takahashi, N.; Kaneda, Y.; Ando, K.; Matsuoka, D. & Kato, T. Parallel Implementation of Dispersive Tsunami Wave
        !       Modeling with a Nesting Algorithm for the 2011 Tohoku Tsunami Pure and Applied Geophysics, 2015, 172, 3455-3472
        ! Or
        !   Baba, T.; Allgeyer, S.; Hossen, J.; Cummins, P. R.; Tsushima, H.; Imai, K.; Yamashita, K. & Kato, T. Accurate numerical
        !   simulation of the far-field tsunami caused by the 2011 Tohoku earthquake, including the effects of Boussinesq dispersion,
        !   seawater density stratification, elastic loading, and gravitational potential change Ocean Modelling, Elsevier BV, 2017,
        !   111, 46–54
        !
        ! BACKGROUND:
        !   This routine is a cell-centred varient of subroutine "linear_dispersive_staggered_matmult". See documentation there
        !   for more information
        !
        real(dp), intent(in) :: elev(:,:), uh(:,:), vh(:,:), &
            dlon, dlat, msl_linear, td1, td2, &
            distance_bottom_edge(:), distance_left_edge(:)
        real(dp), intent(inout) :: solution_uh(:,:), solution_vh(:,:)
        logical, intent(in) :: update_UH, update_VH

        integer(ip) :: i, j
        real(dp) :: R_coslat_dlat, R_coslat_dlon, coslat_jp1, coslat_j, coslat_jm1, d_i_j, dispersive_premult
        real(dp) :: r_dlat, R_coslat_dlon_jph, R_coslat_dlon_jmh, R_coslat_dlat_jph, R_coslat_dlat_jmh
        real(dp) :: dinv_ip1_j, dinv_i_j, dinv_im1_j, &
                    dinv_ip1_jp1, dinv_i_jp1, dinv_im1_jp1, &
                    dinv_ip1_jm1, dinv_i_jm1, dinv_im1_jm1

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(uh, vh, elev, msl_linear, solution_uh, solution_vh, &
        !$OMP                                  dlat, dlon, td1, td2, &
        !$OMP                                  distance_left_edge, distance_bottom_edge, update_UH, update_VH)

        if(update_UH) then
            !$OMP DO
            do j = 2, size(uh, 2)-1

                !if(j == 1 .or. j == size(uh, 2)) then
                !    solution_uh(:,j) = uh(:,j)
                !    cycle
                !end if

                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
                ! cartesian coordinates.
                R_coslat_dlon   = 0.5_dp * (distance_bottom_edge(j+1) + distance_bottom_edge(j  ))
                R_coslat_dlat   = R_coslat_dlon * dlat / dlon

                ! NOTE regarding array bounds. 
                ! - We know j > 1 and j < size(uh,2). 
                ! - Also distance_bottom_edge has size == size(uh,2)+1
                coslat_jp1 = 0.5_dp * (distance_bottom_edge(j+2) + distance_bottom_edge(j+1)) * &
                    dlat / (dlon * distance_left_edge(1))
                !coslat_j   = 0.5_dp * (distance_bottom_edge(j+1) + distance_bottom_edge(j+0)) * &
                !    dlat / (dlon * distance_left_edge(1))
                coslat_jm1 = 0.5_dp * (distance_bottom_edge(j+0) + distance_bottom_edge(j-1)) * &
                    dlat / (dlon * distance_left_edge(1))

                !solution_uh(1,j) = uh(1,j) ! Instead this is set in the caller routine
                !solution_uh(size(uh,1),j) = uh(size(uh,1),j) ! Instead this is set in the caller routine
                !$OMP SIMD
                do i = 2, size(uh, 1) - 1
                    ! UH dispersive term

                    ! Depth at location of uh(i,j)
                    !d_i_j = merge(msl_linear - elev(i, j), 0.0_dp, elev(i,j) < msl_linear)
                    d_i_j = merge(msl_linear - elev(i, j), 0.0_dp, all(elev(i-1:i+1,j) < msl_linear))

#ifdef DISPERSIVE_PEREGRINE_IN_FLUX_FORM
                    !
                    ! Peregrine dispersion
                    !
                    dispersive_premult = d_i_j*d_i_j / (R_coslat_dlon) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_j - td2)/(td1 - td2)))

                    ! Compute various inverse depth terms, which are zeroed at 'dry (at MSL)' cells.
                    ! FIXME: Could avoid recomputing some of these (noting the loop order)
                    dinv_ip1_j = merge(1.0_dp/(msl_linear - elev(i+1,j)), 0.0_dp, elev(i+1,j) < msl_linear)
                    dinv_i_j = merge(1.0_dp/(msl_linear - elev(i,j)), 0.0_dp, elev(i,j) < msl_linear)
                    dinv_im1_j = merge(1.0_dp/(msl_linear - elev(i-1,j)), 0.0_dp, elev(i-1,j) < msl_linear)
                    dinv_ip1_jp1 = merge(1.0_dp/(msl_linear - elev(i+1,j+1)), 0.0_dp, elev(i+1,j+1) < msl_linear)
                    dinv_ip1_jm1 = merge(1.0_dp/(msl_linear - elev(i+1,j-1)), 0.0_dp, elev(i+1,j-1) < msl_linear)
                    dinv_im1_jp1 = merge(1.0_dp/(msl_linear - elev(i-1,j+1)), 0.0_dp, elev(i-1,j+1) < msl_linear)
                    dinv_im1_jm1 = merge(1.0_dp/(msl_linear - elev(i-1,j-1)), 0.0_dp, elev(i-1,j-1) < msl_linear)

                    solution_uh(i,j) = &
                        ! Include all terms in A_uh%*%x here
                        ! h0^2 / (2 R_coslat) * d/dlon
                        dispersive_premult * 0.5_dp * ( &
                            ! lon-diff { 1/(R_coslat) * d/dlon(uh) }
                            1.0_dp / R_coslat_dlon * ( uh(i+1,j) - 2.0_dp * uh(i,j) + uh(i-1,j) ) + &
                            ! lon-diff { 1/(R_coslat) * d/dlat( coslat * vh ) }
                            1.0_dp / R_coslat_dlat * ( &
                                (coslat_jp1 * vh(i+1,j+1) - coslat_jm1 * vh(i+1,j-1))*0.25_dp - &
                                (coslat_jp1 * vh(i-1,j+1) - coslat_jm1 * vh(i-1,j-1))*0.25_dp)  &
                                ) &
                        ! -h0^3 / (6 R_coslat) * d/dlon
                        - dispersive_premult * (1.0_dp/6.0_dp) * d_i_j * ( &
                            ! lon-diff { 1/(R_coslat) * d/dlon(u) }
                            1.0_dp / R_coslat_dlon * (uh(i+1,j)*dinv_ip1_j - 2.0_dp * uh(i,j)*dinv_i_j + uh(i-1,j)*dinv_im1_j ) + &
                            ! lon-diff { 1/(R_coslat) * d/dlat( coslat * v) }
                            1.0_dp / R_coslat_dlat * ( &
                                (coslat_jp1 * vh(i+1,j+1)*dinv_ip1_jp1 - coslat_jm1 * vh(i+1,j-1)*dinv_ip1_jm1)*0.25_dp - &
                                (coslat_jp1 * vh(i-1,j+1)*dinv_im1_jp1 - coslat_jm1 * vh(i-1,j-1)*dinv_im1_jm1)*0.25_dp)  &
                                )

#else
                    !
                    ! Simplest dispersive model, matching JAGURS
                    !

                    dispersive_premult = d_i_j*d_i_j / (3.0_dp * R_coslat_dlon) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_j - td2)/(td1 - td2)))

                    solution_uh(i,j) = &
                        ! Include all terms in A_uh%*%x here
                        ! h0^2 / (3 R_coslat) * d/dlon
                        dispersive_premult * ( &
                            ! lon-diff { 1/(R_coslat) * d/dlon(uh) }
                            1.0_dp / R_coslat_dlon * ( uh(i+1,j) - 2.0_dp * uh(i,j) + uh(i-1,j) ) + &
                            ! lon-diff { 1/(R_coslat) * d/dlat( coslat * vh ) }
                            1.0_dp / R_coslat_dlat * ( &
                                (coslat_jp1 * vh(i+1,j+1) - coslat_jm1 * vh(i+1,j-1))*0.25_dp - &
                                (coslat_jp1 * vh(i-1,j+1) - coslat_jm1 * vh(i-1,j-1))*0.25_dp)  &
                                )
#endif
                end do
            end do
            !$OMP END DO NOWAIT
        end if

        if(update_VH) then
            !$OMP DO
            do j = 2, size(vh, 2) - 1

                !if(j == 1 .or. j == size(vh, 2)) then
                !    solution_vh(:,j) = vh(:,j)
                !    cycle
                !end if

                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
                ! cartesian coordinates.
                r_dlat = distance_left_edge(1)
                R_coslat_dlon_jph = distance_bottom_edge(j+1)
                R_coslat_dlon_jmh = distance_bottom_edge(j+0)
                R_coslat_dlat_jph = R_coslat_dlon_jph * dlat / dlon
                R_coslat_dlat_jmh = R_coslat_dlon_jmh * dlat / dlon

                coslat_jp1 = 0.5_dp*(distance_bottom_edge(j+2)+distance_bottom_edge(j+1)) * dlat / (distance_left_edge(1) * dlon)
                coslat_j   = 0.5_dp*(distance_bottom_edge(j+1)+distance_bottom_edge(j+0)) * dlat / (distance_left_edge(1) * dlon)
                coslat_jm1 = 0.5_dp*(distance_bottom_edge(j+0)+distance_bottom_edge(j-1)) * dlat / (distance_left_edge(1) * dlon)

                !solution_vh(1,j) = vh(1,j) ! Instead this is set in the caller routine
                !solution_vh(size(vh,1), j) = vh(size(vh,1), j)
                !$OMP SIMD
                do i = 2, size(vh, 1) - 1

                    ! Depth at location of vh(i,j) -- or zero if either neighbour is dry (at MSL)
                    !d_i_jph = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i,j+1)), 0.0_dp, &
                    !    elev(i,j) < msl_linear .and. elev(i, j+1) < msl_linear)
                    !d_i_j = merge(msl_linear - elev(i,j), 0.0_dp, elev(i,j) < msl_linear)
                    d_i_j = merge(msl_linear - elev(i,j), 0.0_dp, all(elev(i,j-1:j+1) < msl_linear))

#ifdef DISPERSIVE_PEREGRINE_IN_FLUX_FORM
                    !
                    ! Peregrine dispersion
                    !
                    dispersive_premult = d_i_j*d_i_j / (r_dlat) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_j - td2)/(td1 - td2)))

                    ! Compute various inverse depth terms, which are zeroed at 'dry (at MSL)' cells.
                    ! FIXME: Could avoid recomputing some of these (noting the loop order)
                    dinv_ip1_jp1 = merge(1.0_dp / (msl_linear - elev(i+1,j+1)), 0.0_dp,  elev(i+1, j+1) < msl_linear)
                    dinv_i_jp1 = merge(1.0_dp / (msl_linear - elev(i,j+1)), 0.0_dp,  elev(i, j+1) < msl_linear)
                    dinv_im1_jp1 = merge(1.0_dp / (msl_linear - elev(i-1,j+1)), 0.0_dp,  elev(i-1, j+1) < msl_linear)
                    dinv_ip1_j = merge(1.0_dp / (msl_linear - elev(i+1,j)), 0.0_dp,  elev(i+1, j) < msl_linear)
                    dinv_im1_j = merge(1.0_dp / (msl_linear - elev(i-1,j)), 0.0_dp,  elev(i-1, j) < msl_linear)
                    dinv_i_j = merge(1.0_dp / (msl_linear - elev(i,j)), 0.0_dp,  elev(i, j) < msl_linear)
                    dinv_ip1_jm1 = merge(1.0_dp / (msl_linear - elev(i+1,j-1)), 0.0_dp,  elev(i+1, j-1) < msl_linear)
                    dinv_i_jm1 = merge(1.0_dp / (msl_linear - elev(i,j-1)), 0.0_dp,  elev(i, j-1) < msl_linear) 
                    dinv_im1_jm1 = merge(1.0_dp / (msl_linear - elev(i-1,j-1)), 0.0_dp,  elev(i-1, j-1) < msl_linear)

                    solution_vh(i,j) = &
                        ! Include all terms in A_vh%*%x here
                        ! h0**2 / (2 * R) * d/dlat * (
                        dispersive_premult * 0.5_dp * ( &
                            ! uh term = lat-diff{ 1/Rcoslat duh/dlon }
                            !   1/(R coslat) d(uh)/dlon @ i, j+1/2, by mean of centred differences at j+1 and j.
                            ( 1.0_dp / R_coslat_dlon_jph * (uh(i+1, j+1) - uh(i-1,j+1) + uh(i+1,j) - uh(i-1,j))*0.25_dp - &
                            !   1/(R coslat) d(uh)/dlon @ i, j-1/2, by mean of centred differences at j and j-1
                              1.0_dp / R_coslat_dlon_jmh * (uh(i+1, j-1) - uh(i-1,j-1) + uh(i+1,j) - uh(i-1,j))*0.25_dp ) + &
                            ! vh term = lat-diff{ 1/Rcoslat d(coslatvh)/dlat }
                            !   1/(R coslat ) d(coslat vh)/dlat @ i, j+1/2
                            ( 1.0_dp / R_coslat_dlat_jph * ( coslat_jp1 * vh(i, j+1) - coslat_j   * vh(i, j  )) - &
                            !   1/(R coslat ) d(coslat vh)/dlat @ i, j-1/2
                              1.0_dp / R_coslat_dlat_jmh * ( coslat_j   * vh(i, j  ) - coslat_jm1 * vh(i, j-1)) ) ) &
                        - dispersive_premult * (1.0_dp / 6.0_dp) * d_i_j * ( &
                            ! uh term = lat-diff{ 1/Rcoslat du/dlon }
                            !   1/(R coslat) d(u)/dlon @ i, j+1/2, by mean of centred differences at j+1 and j.
                            ( 1.0_dp / R_coslat_dlon_jph * (uh(i+1, j+1)*dinv_ip1_jp1 - uh(i-1,j+1)*dinv_im1_jp1 + &
                                                            uh(i+1,j)*dinv_ip1_j      - uh(i-1,j)*dinv_im1_j)*0.25_dp - &
                            !   1/(R coslat) d(u)/dlon @ i, j-1/2, by mean of centred differences at j and j-1
                              1.0_dp / R_coslat_dlon_jmh * (uh(i+1, j-1)*dinv_ip1_jm1 - uh(i-1,j-1)*dinv_im1_jm1 + &
                                                            uh(i+1,j)*dinv_ip1_j      - uh(i-1,j)*dinv_im1_j)*0.25_dp ) + &
                            ! vh term = lat-diff{ 1/Rcoslat d(coslatv)/dlat }
                            !   1/(R coslat ) d(coslat v)/dlat @ i, j+1/2
                            ( 1.0_dp / R_coslat_dlat_jph * ( coslat_jp1*vh(i, j+1)*dinv_i_jp1 - coslat_j*vh(i, j)*dinv_i_j) - &
                            !   1/(R coslat ) d(coslat v)/dlat @ i, j-1/2
                              1.0_dp / R_coslat_dlat_jmh * ( coslat_j*vh(i, j)*dinv_i_j - coslat_jm1*vh(i, j-1)*dinv_i_jm1) ) )
#else
                    !
                    ! Simplest dispersive model, matching JAGURS
                    !
                    dispersive_premult = d_i_j*d_i_j / (3.0_dp * r_dlat) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_j - td2)/(td1 - td2)))

                    solution_vh(i,j) = &
                        ! Include all terms in A_vh%*%x here
                        ! h0**2 / (3 * R) * d/dlat * (
                        dispersive_premult * ( &
                            ! uh term = lat-diff{ 1/Rcoslat duh/dlon }
                            !   1/(R coslat) d(uh)/dlon @ i, j+1/2, by mean of centred differences at j+1 and j.
                            ( 1.0_dp / R_coslat_dlon_jph * (uh(i+1, j+1) - uh(i-1,j+1) + uh(i+1,j) - uh(i-1,j))*0.25_dp - &
                            !   1/(R coslat) d(uh)/dlon @ i, j-1/2, by mean of centred differences at j and j-1
                              1.0_dp / R_coslat_dlon_jmh * (uh(i+1, j-1) - uh(i-1,j-1) + uh(i+1,j) - uh(i-1,j))*0.25_dp ) + &
                            ! vh term = lat-diff{ 1/Rcoslat d(coslatvh)/dlat }
                            !   1/(R coslat ) d(coslat vh)/dlat @ i, j+1/2
                            ( 1.0_dp / R_coslat_dlat_jph * ( coslat_jp1 * vh(i, j+1) - coslat_j   * vh(i, j  )) - &
                            !   1/(R coslat ) d(coslat vh)/dlat @ i, j-1/2
                              1.0_dp / R_coslat_dlat_jmh * ( coslat_j   * vh(i, j  ) - coslat_jm1 * vh(i, j-1)) ) )
#endif
                end do
            end do
            !$OMP END DO
        end if

        !$OMP END PARALLEL

    end subroutine


    subroutine linear_dispersive_sweep_staggered_TRIDIAG(elev, uh, vh, &
            RHS_uh, RHS_vh, &
            msl_linear, distance_bottom_edge, distance_left_edge, td1, td2, &
            dlon, dlat, update_UH, update_VH)

        real(dp), intent(in) :: elev(:,:), RHS_uh(:,:), RHS_vh(:,:), &
            dlon, dlat, msl_linear, td1, td2, &
            distance_bottom_edge(:), distance_left_edge(:)
        real(dp), intent(inout) :: uh(:,:), vh(:,:)
        logical, intent(in) :: update_UH, update_VH

        integer(ip) :: i, j, ip2, jp2
        integer :: N, INFO
        real(dp) :: R_coslat_dlat, R_coslat_dlon, coslat_jph, coslat_jmh, d_iph_j, d_i_jph, dispersive_premult
        real(dp) :: r_dlat, R_coslat_dlon_jp1, R_coslat_dlon_j, R_coslat_dlat_jp1, R_coslat_dlat_j, coslat_jp1h
        ! Arrays for tridagonal solves
        real(dp) :: diagonalA_uh(size(uh, 1)), lowerA_uh(size(uh, 1)), upperA_uh(size(uh, 1)), b_uh(size(uh,1))
        real(dp) :: diagonalA_vh(size(vh, 2)), lowerA_vh(size(vh, 2)), upperA_vh(size(vh, 2)), b_vh(size(vh,2))
        real(dp) :: dinv_ip1h_j, dinv_iph_j, dinv_imh_j, &
                    dinv_ip1_jph, dinv_ip1_jmh, &
                    dinv_i_jph, dinv_i_jmh, &
                    dinv_iph_jp1, dinv_imh_jp1, &
                    dinv_i_jp1h

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(uh, vh, elev, RHS_uh, RHS_vh, msl_linear, &
        !$OMP                                  dlat, dlon, td1, td2, &
        !$OMP                                  distance_left_edge, distance_bottom_edge, update_UH, update_VH)

        if(update_UH) then
            !$OMP DO
            do j = 2, size(uh, 2) - 1

                !if(j == 1 .or. j == size(uh,2)) cycle

                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
                ! cartesian coordinates.
                R_coslat_dlon   = 0.5_dp * (distance_bottom_edge(j+1) + distance_bottom_edge(j  )) ! x-cell-distance at uh(i,j)
                R_coslat_dlat   = R_coslat_dlon * dlat / dlon
                coslat_jph = distance_bottom_edge(j+1) * dlat / (dlon * distance_left_edge(1))
                coslat_jmh = distance_bottom_edge(j+0) * dlat / (dlon * distance_left_edge(1))

                ! We do not update the boundary values
                N = size(uh, 1)
                diagonalA_uh(1) = 1.0_dp
                diagonalA_uh(N) = 1.0_dp
                lowerA_uh(1) = 0.0_dp
                lowerA_uh(N) = 0.0_dp
                upperA_uh(1) = 0.0_dp
                upperA_uh(N) = 0.0_dp
                !b_uh(1) = uh(1,j)
                !b_uh(N) = uh(N,j)
                b_uh(1) = RHS_uh(1,j)
                b_uh(N) = RHS_uh(N,j)

                !$OMP SIMD
                do i = 2, size(uh, 1) - 1
                    ! UH dispersive term

                    ! Depth at location of uh(i,j) -- or zero if either neighbour is dry (at MSL)
                    d_iph_j = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i+1,j)), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i+1, j) < msl_linear)
#ifdef DISPERSIVE_PEREGRINE_IN_FLUX_FORM
                    !
                    ! Peregrine dispersion
                    !
                    dispersive_premult = d_iph_j*d_iph_j / (R_coslat_dlon) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_iph_j - td2)/(td1 - td2)))

                    ! Compute various inverse depth terms, which are zeroed at 'dry (at MSL)' cells.
                    ! FIXME: Could avoid recomputing some of these (noting the loop order)
                    ip2 = min(i+2, size(uh, 1)) ! Safe i+2 index
                    dinv_ip1h_j = merge(1.0_dp/(0.5_dp * (msl_linear - elev(i+1,j) + msl_linear - elev(ip2,j))), 0.0_dp, &
                        elev(i+1,j) < msl_linear .and. elev(ip2, j) < msl_linear)
                    dinv_iph_j = merge(1.0_dp/d_iph_j, 0.0_dp, d_iph_j > 0.0_dp)
                    dinv_imh_j = merge(1.0_dp/(0.5_dp * (msl_linear - elev(i-1,j) + msl_linear - elev(i,j))), 0.0_dp, &
                        elev(i-1,j) < msl_linear .and. elev(i, j) < msl_linear)
                    dinv_ip1_jph = merge(1.0_dp/(0.5_dp * (msl_linear - elev(i+1,j+1) + msl_linear - elev(i+1,j))), 0.0_dp, &
                        elev(i+1,j+1) < msl_linear .and. elev(i+1, j) < msl_linear)
                    dinv_ip1_jmh = merge(1.0_dp/(0.5_dp * (msl_linear - elev(i+1,j) + msl_linear - elev(i+1,j-1))), 0.0_dp, &
                        elev(i+1,j) < msl_linear .and. elev(i+1, j-1) < msl_linear)
                    dinv_i_jph = merge(1.0_dp/(0.5_dp * (msl_linear - elev(i,j+1) + msl_linear - elev(i,j))), 0.0_dp, &
                        elev(i,j+1) < msl_linear .and. elev(i, j) < msl_linear)
                    dinv_i_jmh = merge(1.0_dp/(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i,j-1))), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i, j-1) < msl_linear)

                    ! Terms related to the uh derivative (with 1.0 --> time update)
                    diagonalA_uh(i) = 1.0_dp + 2.0_dp * dispersive_premult * 0.5_dp/R_coslat_dlon &
                        -2.0_dp * dispersive_premult * (1.0_dp/6.0_dp) * d_iph_j / R_coslat_dlon * dinv_iph_j
                    lowerA_uh(i) = -dispersive_premult * 0.5_dp/R_coslat_dlon &
                        + dispersive_premult * (1.0_dp/6.0_dp) * d_iph_j / R_coslat_dlon * dinv_imh_j
                    upperA_uh(i) = -dispersive_premult * 0.5_dp/R_coslat_dlon &
                        + dispersive_premult * (1.0_dp/6.0_dp) * d_iph_j / R_coslat_dlon * dinv_ip1h_j

                    b_uh(i) = RHS_uh(i,j) + &
                        dispersive_premult * 0.5_dp * ( &
                            ! lon-diff{ 1/(R_coslat) * d/dlat( coslat * vh ) }
                            1.0_dp / R_coslat_dlat * ( (coslat_jph * vh(i+1,j) - coslat_jmh * vh(i+1,j-1)) - &
                                                       (coslat_jph * vh(i  ,j) - coslat_jmh * vh(i  ,j-1)) ) ) & 
                        -dispersive_premult * (1.0_dp/6.0_dp) * d_iph_j * ( &
                            ! lon-diff{ 1/(R_coslat) * d/dlat( coslat * v ) }
                            1.0_dp / R_coslat_dlat*((coslat_jph*vh(i+1,j)*dinv_ip1_jph - coslat_jmh*vh(i+1,j-1)*dinv_ip1_jmh) - &
                                                    (coslat_jph*vh(i  ,j)*dinv_i_jph   - coslat_jmh*vh(i  ,j-1)*dinv_i_jmh )) )

#else
                    !
                    ! Simplest dispersive model, matching JAGURS
                    !
                    dispersive_premult = d_iph_j*d_iph_j / (3.0_dp * R_coslat_dlon) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_iph_j - td2)/(td1 - td2)))

                    ! Terms related to the uh derivative (with 1.0 --> time update)
                    diagonalA_uh(i) = 1.0_dp + 2.0_dp * dispersive_premult * 1.0_dp/R_coslat_dlon
                    lowerA_uh(i) = -dispersive_premult * 1.0_dp/R_coslat_dlon
                    upperA_uh(i) = -dispersive_premult * 1.0_dp/R_coslat_dlon

                    b_uh(i) = RHS_uh(i,j) + &
                        dispersive_premult * ( &
                            ! lon-diff{ 1/(R_coslat) * d/dlat( coslat * vh ) }
                            1.0_dp / R_coslat_dlat * ( (coslat_jph * vh(i+1,j) - coslat_jmh * vh(i+1,j-1)) - &
                                                       (coslat_jph * vh(i  ,j) - coslat_jmh * vh(i  ,j-1)) ) )
#endif
                end do
                ! Compute the solution with lapack's tridiagonal solver
#ifdef REALFLOAT
                call sgtsv(N, 1, lowerA_uh(2:N), diagonalA_uh, upperA_uh(1:N-1), b_uh, N, INFO)
#else
                call dgtsv(N, 1, lowerA_uh(2:N), diagonalA_uh, upperA_uh(1:N-1), b_uh, N, INFO)
#endif
                uh(:,j) = b_uh

            end do
            !$OMP END DO
        end if

        if(update_VH) then
            !$OMP DO
            ! For the tridiagonal version, need to have 'i' as the outer loop, sweeping by j.
            do i = 2, size(vh, 1) - 1

                N = size(vh, 2)
                diagonalA_vh(1) = 1.0_dp
                diagonalA_vh(N) = 1.0_dp
                lowerA_vh(1) = 0.0_dp
                lowerA_vh(N) = 0.0_dp
                upperA_vh(1) = 0.0_dp
                upperA_vh(N) = 0.0_dp
                !b_vh(1) = vh(i,1)
                !b_vh(N) = vh(i,N)
                b_vh(1) = RHS_vh(i,1)
                b_vh(N) = RHS_vh(i,N)

                !$OMP SIMD
                do j = 2, size(vh, 2) - 1
                    !if(j == 1 .or. j == size(vh,2)) cycle

                    ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
                    ! cartesian coordinates.
                    r_dlat = distance_left_edge(1)
                    R_coslat_dlon_jp1 = 0.5_dp * ( distance_bottom_edge(j+2) + distance_bottom_edge(j+1) )
                    R_coslat_dlon_j   = 0.5_dp * ( distance_bottom_edge(j+1) + distance_bottom_edge(j+0) )
                    R_coslat_dlat_jp1 = R_coslat_dlon_jp1 * dlat / dlon
                    R_coslat_dlat_j   = R_coslat_dlon_j   * dlat / dlon

                    coslat_jp1h = distance_bottom_edge(j+2) * dlat / ( distance_left_edge(1)  * dlon )
                    coslat_jph  = distance_bottom_edge(j+1) * dlat / ( distance_left_edge(1)  * dlon )
                    coslat_jmh  = distance_bottom_edge(j+0) * dlat / ( distance_left_edge(1)  * dlon )

                    ! Depth at location of vh(i,j) -- or zero if either neighbour is dry (at MSL)
                    d_i_jph = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i,j+1)), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i, j+1) < msl_linear)
#ifdef DISPERSIVE_PEREGRINE_IN_FLUX_FORM
                    !
                    ! Peregrine dispersion
                    !
                    dispersive_premult = d_i_jph*d_i_jph / (r_dlat) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_jph - td2)/(td1 - td2)))

                    ! Compute various inverse depth terms, which are zeroed at 'dry (at MSL)' cells.
                    ! FIXME: Could avoid recomputing some of these (noting the loop order)
                    dinv_iph_jp1 = merge(1.0_dp / (0.5_dp * (msl_linear - elev(i,j+1) + msl_linear - elev(i+1,j+1))), 0.0_dp, &
                        elev(i,j+1) < msl_linear .and. elev(i+1, j+1) < msl_linear)
                    dinv_imh_jp1 = merge(1.0_dp / (0.5_dp * (msl_linear - elev(i,j+1) + msl_linear - elev(i-1,j+1))), 0.0_dp, &
                        elev(i,j+1) < msl_linear .and. elev(i-1, j+1) < msl_linear)
                    dinv_iph_j = merge(1.0_dp / (0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i+1,j))), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i+1, j) < msl_linear)
                    dinv_imh_j = merge(1.0_dp / (0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i-1,j))), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i-1, j) < msl_linear)
                    jp2 = min(j+2, size(vh, 2)) ! Safe j+2
                    dinv_i_jp1h = merge(1.0_dp / (0.5_dp * (msl_linear - elev(i,j+1) + msl_linear - elev(i,jp2))), 0.0_dp, &
                        elev(i,j+1) < msl_linear .and. elev(i, jp2) < msl_linear)
                    dinv_i_jph = merge(1.0_dp/d_i_jph, 0.0_dp, d_i_jph > 0.0_dp)
                    dinv_i_jmh = merge(1.0_dp / (0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i,j-1))), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i, j-1) < msl_linear)


                    ! Terms related to the vh derivative (with 1.0 --> time update)
                    diagonalA_vh(j) = 1.0_dp &
                        + dispersive_premult * 0.5_dp * &
                            (1.0_dp / R_coslat_dlat_jp1 * coslat_jph + &
                             1.0_dp / R_coslat_dlat_j   * coslat_jph ) &
                        -dispersive_premult * (1.0_dp/6.0_dp) * d_i_jph * (&
                            1.0_dp / R_coslat_dlat_jp1 * coslat_jph * dinv_i_jph + &
                            1.0_dp / R_coslat_dlat_j   * coslat_jph * dinv_i_jph ) 
                    lowerA_vh(j) = -dispersive_premult * 0.5_dp / R_coslat_dlat_j   * coslat_jmh &
                        +dispersive_premult * (1.0_dp/6.0_dp) * d_i_jph /R_coslat_dlat_j * coslat_jmh*dinv_i_jmh 
                    upperA_vh(j) = -dispersive_premult * 0.5_dp / R_coslat_dlat_jp1 * coslat_jp1h &
                        +dispersive_premult * (1.0_dp/6.0_dp) * d_i_jph /R_coslat_dlat_jp1 * coslat_jp1h*dinv_i_jp1h

                    b_vh(j) = RHS_vh(i,j) + &
                        dispersive_premult * 0.5_dp * ( &
                            ! uh term
                            ( 1.0_dp / R_coslat_dlon_jp1 * (uh(i  , j+1) - uh(i-1, j+1)) - &
                              1.0_dp / R_coslat_dlon_j   * (uh(i  , j  ) - uh(i-1, j  )) ) ) &
                        - dispersive_premult * (1.0_dp/6.0_dp) * d_i_jph * ( &
                            ! lat-diff{ 1/R_coslat du/dlon }
                            ( 1.0_dp / R_coslat_dlon_jp1 * (uh(i  , j+1)*dinv_iph_jp1 - uh(i-1, j+1)*dinv_imh_jp1) - &
                              1.0_dp / R_coslat_dlon_j   * (uh(i  , j  )*dinv_iph_j   - uh(i-1, j  )*dinv_imh_j  )  ))

#else
                    !
                    ! Simplest dispersive model, matching JAGURS
                    !

                    dispersive_premult = d_i_jph*d_i_jph / (3.0_dp * r_dlat) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_jph - td2)/(td1 - td2)))


                    ! Terms related to the vh derivative (with 1.0 --> time update)
                    diagonalA_vh(j) = 1.0_dp + dispersive_premult * (1.0_dp / R_coslat_dlat_jp1 * coslat_jph + &
                                                                     1.0_dp / R_coslat_dlat_j   * coslat_jph )
                    lowerA_vh(j) = -dispersive_premult * 1.0_dp / R_coslat_dlat_j   * coslat_jmh
                    upperA_vh(j) = -dispersive_premult * 1.0_dp / R_coslat_dlat_jp1 * coslat_jp1h

                    b_vh(j) = RHS_vh(i,j) + &
                        dispersive_premult * ( &
                            ! uh term
                            ( 1.0_dp / R_coslat_dlon_jp1 * (uh(i  , j+1) - uh(i-1, j+1)) - &
                              1.0_dp / R_coslat_dlon_j   * (uh(i  , j  ) - uh(i-1, j  )) ) )
#endif
                end do
                ! Compute the solution with lapack's tridiagonal solver
#ifdef REALFLOAT
                call sgtsv(N, 1, lowerA_vh(2:N), diagonalA_vh, upperA_vh(1:N-1), b_vh, N, INFO)
#else
                call dgtsv(N, 1, lowerA_vh(2:N), diagonalA_vh, upperA_vh(1:N-1), b_vh, N, INFO)
#endif
                vh(i,:) = b_vh

            end do
            !$OMP END DO
        end if

        !$OMP END PARALLEL

    end subroutine

    subroutine linear_dispersive_sweep_cellcentred_TRIDIAG(elev, uh, vh, &
            RHS_uh, RHS_vh, &
            msl_linear, distance_bottom_edge, distance_left_edge, td1, td2, &
            dlon, dlat, update_UH, update_VH)

        real(dp), intent(in) :: elev(:,:), RHS_uh(:,:), RHS_vh(:,:), &
            dlon, dlat, msl_linear, td1, td2, &
            distance_bottom_edge(:), distance_left_edge(:)
        real(dp), intent(inout) :: uh(:,:), vh(:,:)
        logical, intent(in) :: update_UH, update_VH

        integer(ip) :: i, j
        integer :: N, INFO
        real(dp) :: R_coslat_dlat, R_coslat_dlon, coslat_jp1, coslat_jm1, coslat_j, d_i_j,dispersive_premult
        real(dp) :: r_dlat, R_coslat_dlon_jp1, R_coslat_dlon_j, R_coslat_dlat_jp1, R_coslat_dlat_j
        real(dp) :: R_coslat_dlon_jph, R_coslat_dlon_jmh, R_coslat_dlat_jph, R_coslat_dlat_jmh
        real(dp) :: dinv_ip1_j, dinv_i_j, dinv_im1_j, &
                    dinv_ip1_jp1, dinv_i_jp1, dinv_im1_jp1, &
                    dinv_ip1_jm1, dinv_i_jm1, dinv_im1_jm1
        ! Arrays for tridagonal solves
        real(dp) :: diagonalA_uh(size(uh, 1)), lowerA_uh(size(uh, 1)), upperA_uh(size(uh, 1)), b_uh(size(uh,1))
        real(dp) :: diagonalA_vh(size(vh, 2)), lowerA_vh(size(vh, 2)), upperA_vh(size(vh, 2)), b_vh(size(vh,2))

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(uh, vh, elev, RHS_uh, RHS_vh, msl_linear, &
        !$OMP                                  dlat, dlon, td1, td2, &
        !$OMP                                  distance_left_edge, distance_bottom_edge, update_UH, update_VH)

        if(update_UH) then
            !$OMP DO
            do j = 2, size(uh, 2) - 1

                !if(j == 1 .or. j == size(uh, 2)) cycle

                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
                ! cartesian coordinates.
                R_coslat_dlon   = 0.5_dp * (distance_bottom_edge(j+1) + distance_bottom_edge(j  ))
                R_coslat_dlat   = R_coslat_dlon * dlat / dlon

                ! NOTE regarding array bounds. 
                ! - We know j > 1 and j < size(uh,2). 
                ! - Also distance_bottom_edge has size == size(uh,2)+1
                coslat_jp1 = 0.5_dp * (distance_bottom_edge(j+2) + distance_bottom_edge(j+1)) * &
                    dlat / (dlon * distance_left_edge(1))
                coslat_jm1 = 0.5_dp * (distance_bottom_edge(j+0) + distance_bottom_edge(j-1)) * &
                    dlat / (dlon * distance_left_edge(1))

                ! We do not update the boundary values
                N = size(uh, 1)
                diagonalA_uh(1) = 1.0_dp
                diagonalA_uh(N) = 1.0_dp
                lowerA_uh(1) = 0.0_dp
                lowerA_uh(N) = 0.0_dp
                upperA_uh(1) = 0.0_dp
                upperA_uh(N) = 0.0_dp
                !b_uh(1) = uh(1,j)
                !b_uh(N) = uh(N,j)
                b_uh(1) = RHS_uh(1,j)
                b_uh(N) = RHS_uh(N,j)

                !$OMP SIMD
                do i = 2, size(uh, 1) - 1
                    ! UH dispersive term

                    ! Depth at location of uh(i,j) -- or zero if either neighbour is dry (at MSL)
                    !d_i_j = merge((msl_linear - elev(i,j)), 0.0_dp, elev(i,j) < msl_linear)
                    d_i_j = merge((msl_linear - elev(i,j)), 0.0_dp, all(elev(i-1:i+1,j) < msl_linear))

#ifdef DISPERSIVE_PEREGRINE_IN_FLUX_FORM
                    !
                    ! Peregrine dispersion
                    !
                    dispersive_premult = d_i_j*d_i_j / (R_coslat_dlon) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_j - td2)/(td1 - td2)))

                    ! Compute various inverse depth terms, which are zeroed at 'dry (at MSL)' cells.
                    ! FIXME: Could avoid recomputing some of these (noting the loop order)
                    dinv_ip1_j = merge(1.0_dp/(msl_linear - elev(i+1,j)), 0.0_dp, elev(i+1,j) < msl_linear)
                    dinv_i_j = merge(1.0_dp/(msl_linear - elev(i,j)), 0.0_dp, elev(i,j) < msl_linear)
                    dinv_im1_j = merge(1.0_dp/(msl_linear - elev(i-1,j)), 0.0_dp, elev(i-1,j) < msl_linear)
                    dinv_ip1_jp1 = merge(1.0_dp/(msl_linear - elev(i+1,j+1)), 0.0_dp, elev(i+1,j+1) < msl_linear)
                    dinv_ip1_jm1 = merge(1.0_dp/(msl_linear - elev(i+1,j-1)), 0.0_dp, elev(i+1,j-1) < msl_linear)
                    dinv_im1_jp1 = merge(1.0_dp/(msl_linear - elev(i-1,j+1)), 0.0_dp, elev(i-1,j+1) < msl_linear)
                    dinv_im1_jm1 = merge(1.0_dp/(msl_linear - elev(i-1,j-1)), 0.0_dp, elev(i-1,j-1) < msl_linear)

                    ! Terms related to the uh derivative (with 1.0 --> time update)
                    diagonalA_uh(i) = 1.0_dp + 2.0_dp * dispersive_premult * 0.5_dp/R_coslat_dlon &
                         -2.0_dp * dispersive_premult * (1.0_dp/6.0_dp) * d_i_j /R_coslat_dlon * dinv_i_j
                    lowerA_uh(i) = -dispersive_premult * 0.5_dp/R_coslat_dlon &
                        + dispersive_premult * (1.0_dp/6.0_dp) * d_i_j /R_coslat_dlon * dinv_im1_j
                    upperA_uh(i) = -dispersive_premult * 0.5_dp/R_coslat_dlon &
                        + dispersive_premult * (1.0_dp/6.0_dp) * d_i_j /R_coslat_dlon * dinv_ip1_j

                    b_uh(i) = RHS_uh(i,j) + &
                        dispersive_premult * 0.5_dp * ( &
                            ! lon-diff{ 1/(R_coslat) * d/dlat( coslat * vh ) }
                            1.0_dp / R_coslat_dlat * ( &
                                (coslat_jp1 * vh(i+1,j+1) - coslat_jm1 * vh(i+1,j-1))*0.25_dp - &
                                (coslat_jp1 * vh(i-1,j+1) - coslat_jm1 * vh(i-1,j-1))*0.25_dp ) &
                            ) &
                        - dispersive_premult * (1.0_dp/6.0_dp) * d_i_j * ( &
                            ! lon-diff { 1/(R_coslat) * d/dlat( coslat * v) }
                            1.0_dp / R_coslat_dlat * ( &
                                (coslat_jp1 * vh(i+1,j+1)*dinv_ip1_jp1 - coslat_jm1 * vh(i+1,j-1)*dinv_ip1_jm1)*0.25_dp - &
                                (coslat_jp1 * vh(i-1,j+1)*dinv_im1_jp1 - coslat_jm1 * vh(i-1,j-1)*dinv_im1_jm1)*0.25_dp)  &
                                )
#else
                    !
                    ! Simplest dispersive model, matching JAGURS
                    !
                    dispersive_premult = d_i_j*d_i_j / (3.0_dp * R_coslat_dlon) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_j - td2)/(td1 - td2)))


                    ! Terms related to the uh derivative (with 1.0 --> time update)
                    diagonalA_uh(i) = 1.0_dp + 2.0_dp * dispersive_premult * 1.0_dp/R_coslat_dlon
                    lowerA_uh(i) = -dispersive_premult * 1.0_dp/R_coslat_dlon
                    upperA_uh(i) = -dispersive_premult * 1.0_dp/R_coslat_dlon

                    b_uh(i) = RHS_uh(i,j) + &
                        dispersive_premult * ( &
                            ! 1/(R_coslat) * d/dlat( coslat * vh )
                            1.0_dp / R_coslat_dlat * ( &
                                (coslat_jp1 * vh(i+1,j+1) - coslat_jm1 * vh(i+1,j-1))*0.25_dp - &
                                (coslat_jp1 * vh(i-1,j+1) - coslat_jm1 * vh(i-1,j-1))*0.25_dp ) &
                            )
#endif
                end do
                ! Compute the solution with lapack's tridiagonal solver
#ifdef REALFLOAT
                call sgtsv(N, 1, lowerA_uh(2:N), diagonalA_uh, upperA_uh(1:N-1), b_uh, N, INFO)
#else
                call dgtsv(N, 1, lowerA_uh(2:N), diagonalA_uh, upperA_uh(1:N-1), b_uh, N, INFO)
#endif
                uh(:,j) = b_uh

            end do
            !$OMP END DO
        end if

        if(update_VH) then
            !$OMP DO
            ! For the tridiagonal version, need to have 'i' as the outer loop, sweeping by j.
            do i = 2, size(vh, 1) - 1

                N = size(vh, 2)
                diagonalA_vh(1) = 1.0_dp
                diagonalA_vh(N) = 1.0_dp
                lowerA_vh(1) = 0.0_dp
                lowerA_vh(N) = 0.0_dp
                upperA_vh(1) = 0.0_dp
                upperA_vh(N) = 0.0_dp
                !b_vh(1) = vh(i,1)
                !b_vh(N) = vh(i,N)
                b_vh(1) = RHS_vh(i,1)
                b_vh(N) = RHS_vh(i,N)

                !$OMP SIMD
                do j = 2, size(vh, 2) - 1

                    !if(j == 1 .or. j == size(vh, 2)) cycle

                    ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
                    ! cartesian coordinates.
                    r_dlat = distance_left_edge(1)
                    R_coslat_dlon_jph = distance_bottom_edge(j+1)
                    R_coslat_dlon_jmh = distance_bottom_edge(j+0)
                    R_coslat_dlat_jph = R_coslat_dlon_jph * dlat / dlon
                    R_coslat_dlat_jmh = R_coslat_dlon_jmh * dlat / dlon

                    coslat_jp1 = 0.5_dp*(distance_bottom_edge(j+2)+distance_bottom_edge(j+1)) * &
                        dlat / (distance_left_edge(1) * dlon)
                    coslat_j   = 0.5_dp*(distance_bottom_edge(j+1)+distance_bottom_edge(j+0)) * &
                        dlat / (distance_left_edge(1) * dlon)
                    coslat_jm1 = 0.5_dp*(distance_bottom_edge(j+0)+distance_bottom_edge(j-1)) * &
                        dlat / (distance_left_edge(1) * dlon)

                    ! Depth at location of vh(i,j) -- or zero if either neighbour is dry (at MSL)
                    !d_i_j = merge((msl_linear - elev(i,j)), 0.0_dp, elev(i,j) < msl_linear)
                    d_i_j = merge((msl_linear - elev(i,j)), 0.0_dp, all(elev(i,j-1:j+1) < msl_linear))

#ifdef DISPERSIVE_PEREGRINE_IN_FLUX_FORM
                    !
                    ! Peregrine dispersion
                    !
                    dispersive_premult = d_i_j*d_i_j / (r_dlat) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_j - td2)/(td1 - td2)))

                    ! Compute various inverse depth terms, which are zeroed at 'dry (at MSL)' cells.
                    ! FIXME: Could avoid recomputing some of these (noting the loop order)
                    dinv_ip1_jp1 = merge(1.0_dp / (msl_linear - elev(i+1,j+1)), 0.0_dp,  elev(i+1, j+1) < msl_linear)
                    dinv_i_jp1 = merge(1.0_dp / (msl_linear - elev(i,j+1)), 0.0_dp,  elev(i, j+1) < msl_linear)
                    dinv_im1_jp1 = merge(1.0_dp / (msl_linear - elev(i-1,j+1)), 0.0_dp,  elev(i-1, j+1) < msl_linear)
                    dinv_ip1_j = merge(1.0_dp / (msl_linear - elev(i+1,j)), 0.0_dp,  elev(i+1, j) < msl_linear)
                    dinv_im1_j = merge(1.0_dp / (msl_linear - elev(i-1,j)), 0.0_dp,  elev(i-1, j) < msl_linear)
                    dinv_i_j = merge(1.0_dp / (msl_linear - elev(i,j)), 0.0_dp,  elev(i, j) < msl_linear)
                    dinv_ip1_jm1 = merge(1.0_dp / (msl_linear - elev(i+1,j-1)), 0.0_dp,  elev(i+1, j-1) < msl_linear)
                    dinv_i_jm1 = merge(1.0_dp / (msl_linear - elev(i,j-1)), 0.0_dp,  elev(i, j-1) < msl_linear) 
                    dinv_im1_jm1 = merge(1.0_dp / (msl_linear - elev(i-1,j-1)), 0.0_dp,  elev(i-1, j-1) < msl_linear)

                    ! Terms related to the vh derivative (with 1.0 --> time update)
                    diagonalA_vh(j) = 1.0_dp + dispersive_premult * 0.5_dp * &
                            (1.0_dp / R_coslat_dlat_jph * coslat_j + &
                             1.0_dp / R_coslat_dlat_jmh * coslat_j ) &
                        - dispersive_premult * (1.0_dp/6.0_dp) * d_i_j * &
                            (1.0_dp / R_coslat_dlat_jph * coslat_j * dinv_i_j + &
                             1.0_dp / R_coslat_dlat_jmh * coslat_j * dinv_i_j )
                    lowerA_vh(j) = -dispersive_premult * 0.5_dp / R_coslat_dlat_jmh * coslat_jm1 &
                        + dispersive_premult * (1.0_dp/6.0_dp) * d_i_j / R_coslat_dlat_jmh * coslat_jm1 * dinv_i_jm1
                    upperA_vh(j) = -dispersive_premult * 0.5_dp / R_coslat_dlat_jph * coslat_jp1 &
                        + dispersive_premult * (1.0_dp/6.0_dp) * d_i_j / R_coslat_dlat_jph * coslat_jp1 * dinv_i_jp1

                    b_vh(j) = RHS_vh(i,j) + &
                        dispersive_premult * 0.5_dp * ( &
                            ! uh term
                            ! 1/(R coslat) d(uh)/dlon @ i, j+1/2, by mean of central differences at j+1 and j.
                            ( 1.0_dp / R_coslat_dlon_jph * (uh(i+1, j+1) - uh(i-1,j+1) + uh(i+1,j) - uh(i-1,j))*0.25_dp - &
                            ! 1/(R coslat) d(uh)/dlon @ i, j-1/2, by mean of central differences at j and j-1
                              1.0_dp / R_coslat_dlon_jmh * (uh(i+1, j-1) - uh(i-1,j-1) + uh(i+1,j) - uh(i-1,j))*0.25_dp ) ) &
                        - dispersive_premult * (1.0_dp / 6.0_dp) * d_i_j * ( &
                            ! uh term = lat-diff{ 1/Rcoslat du/dlon }
                            !   1/(R coslat) d(u)/dlon @ i, j+1/2, by mean of centred differences at j+1 and j.
                            ( 1.0_dp / R_coslat_dlon_jph * (uh(i+1, j+1)*dinv_ip1_jp1 - uh(i-1,j+1)*dinv_im1_jp1 + &
                                                            uh(i+1,j)*dinv_ip1_j      - uh(i-1,j)*dinv_im1_j)*0.25_dp - &
                            !   1/(R coslat) d(u)/dlon @ i, j-1/2, by mean of centred differences at j and j-1
                              1.0_dp / R_coslat_dlon_jmh * (uh(i+1, j-1)*dinv_ip1_jm1 - uh(i-1,j-1)*dinv_im1_jm1 + &
                                                            uh(i+1,j)*dinv_ip1_j      - uh(i-1,j)*dinv_im1_j)*0.25_dp ) )

#else
                    !
                    ! Simplest dispersive model, matching JAGURS
                    !

                    dispersive_premult = d_i_j*d_i_j / (3.0_dp * r_dlat) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_j - td2)/(td1 - td2)))

                    ! Terms related to the vh derivative (with 1.0 --> time update)
                    diagonalA_vh(j) = 1.0_dp + dispersive_premult * (1.0_dp / R_coslat_dlat_jph * coslat_j + &
                                                                     1.0_dp / R_coslat_dlat_jmh * coslat_j )
                    lowerA_vh(j) = -dispersive_premult * 1.0_dp / R_coslat_dlat_jmh * coslat_jm1
                    upperA_vh(j) = -dispersive_premult * 1.0_dp / R_coslat_dlat_jph * coslat_jp1

                    b_vh(j) = RHS_vh(i,j) + &
                        dispersive_premult * ( &
                            ! uh term
                            ! 1/(R coslat) d(uh)/dlon @ i, j+1/2, by mean of central differences at j+1 and j.
                            ( 1.0_dp / R_coslat_dlon_jph * (uh(i+1, j+1) - uh(i-1,j+1) + uh(i+1,j) - uh(i-1,j))*0.25_dp - &
                            ! 1/(R coslat) d(uh)/dlon @ i, j-1/2, by mean of central differences at j and j-1
                              1.0_dp / R_coslat_dlon_jmh * (uh(i+1, j-1) - uh(i-1,j-1) + uh(i+1,j) - uh(i-1,j))*0.25_dp ) )
#endif
                end do
                ! Compute the solution with lapack's tridiagonal solver
#ifdef REALFLOAT
                call sgtsv(N, 1, lowerA_vh(2:N), diagonalA_vh, upperA_vh(1:N-1), b_vh, N, INFO)
#else
                call dgtsv(N, 1, lowerA_vh(2:N), diagonalA_vh, upperA_vh(1:N-1), b_vh, N, INFO)
#endif
                vh(i,:) = b_vh

            end do
            !$OMP END DO
        end if

        !$OMP END PARALLEL

    end subroutine

    subroutine linear_dispersive_solve_TRIDIAG(ds, is_staggered_grid, U, &
            dlon, dlat, distance_bottom_edge, distance_left_edge, msl_linear, rhs_is_up_to_date, &
            estimate_solution_forward_in_time, forward_time)
        !! Solve the dispersive equations with tridiagonal iteration, optionally
        !! using forward extrapolation in time to guess the solution
        class(dispersive_solver_type), intent(inout) :: ds
        logical, intent(in) :: is_staggered_grid
            ! If TRUE then assume U is on a staggered grid, otherwise cellcentred grid
        real(dp), intent(inout) :: U(:,:,:)
            ! domain%U after an explicit shallow water update
        real(dp), intent(in) :: dlon, dlat, distance_bottom_edge(:), distance_left_edge(:), msl_linear
            ! cell dx, dy, edge distances, and mean-sea-level for the linear solver
        logical, optional, intent(in) :: rhs_is_up_to_date
        logical, optional, intent(in) :: estimate_solution_forward_in_time
        real(dp), optional, intent(in) :: forward_time

        integer :: i, j, iter, nm1
        real(dp) :: max_err, last_U
        logical :: rhs_ready, quadratic_extrap, did_extrapolate
        integer(ip) :: niter

        rhs_ready = .FALSE.
        if(present(rhs_is_up_to_date)) rhs_ready = rhs_is_up_to_date

        quadratic_extrap = .FALSE.
        if(present(estimate_solution_forward_in_time)) then
            quadratic_extrap = estimate_solution_forward_in_time
            if(quadratic_extrap .and. (.not. present(forward_time))) call generic_stop
        end if


        if(.not. rhs_ready) then
EVOLVE_TIMER_START('disp_rhs')
            ! Get the explicit part of the dispersive terms
            if(is_staggered_grid) then
                call linear_dispersive_staggered_matmult(&
                    ds%last_U(:,:,ELV), ds%last_U(:,:,UH), ds%last_U(:,:,VH), &
                    msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, ds%td1, ds%td2, &
                    ds%Ax(:,:,UH), ds%Ax(:,:,VH), &
                    update_UH=.true., update_VH=.true.)
            else
                call linear_dispersive_cellcentred_matmult(&
                    ds%last_U(:,:,ELV), ds%last_U(:,:,UH), ds%last_U(:,:,VH), &
                    msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, ds%td1, ds%td2, &
                    ds%Ax(:,:,UH), ds%Ax(:,:,VH), &
                    update_UH=.true., update_VH=.true.)
            end if
            ! Setup the right-hand-side term, combining the shallow-water solution with the explicit part
            ! of the dispersive term
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(ds, U)
            nm1 = size(U,1) - 1
            !$OMP DO
            do j = 1, size(U,2)
                if((j == 1) .or. (j == size(U,2))) then
                    ! Do nothing boundary conditions north/south
                    ds%RHS(:,j,UH) = U(:,j,UH)
                else
                    ! Regular case
                    ds%RHS(2:nm1,j,UH) = U(2:nm1,j,UH) - ds%Ax(2:nm1,j,UH)
                    ! Do-nothing boundary conditions east/west
                    ds%RHS(1,j,UH) = U(1,j,UH)
                    ds%RHS(size(U,1),j,UH) = U(size(U,1),j,UH)
                end if
            end do
            !$OMP END DO NOWAIT

            !$OMP DO
            do j = 1, size(U,2)
                if((j == 1) .or. (j == size(U,2))) then
                    ! Do nothing boundary conditions north/south
                    ds%RHS(:,j,VH) = U(:,j,VH)
                else
                    ! Regular case
                    ds%RHS(2:nm1,j,VH) = U(2:nm1,j,VH) - ds%Ax(2:nm1,j,VH) 
                    ! Do nothing boundary conditions east/west
                    ds%RHS(1,j,VH) = U(1,j,VH)
                    ds%RHS(size(U,1),j,VH) = U(size(U,1),j,VH)
                end if
            end do
            !$OMP END DO

            !$OMP END PARALLEL
EVOLVE_TIMER_STOP('disp_rhs')
        end if


        if(quadratic_extrap) then 
EVOLVE_TIMER_START('disp_extrapolate_forward')
            ! Extrapolate forward in time to guess the solution
            call ds%qet%extrapolate_in_time(forward_time, U(:,:,UH:VH), &
                do_nothing_if_missing_times=.true., did_extrapolate=did_extrapolate)
                !Q: Is it possible that extrapolation leads to a sign change of UH or VH,
                !   that violates wet/dry logic? Such as:
                !    - nonzero UH if the cell has just gone dry
                !    - positive UH even if its neighbours are dry in a way that prevents flow
                ! ??
EVOLVE_TIMER_STOP('disp_extrapolate_forward')
            if(did_extrapolate) then
                ! If we could extrapolate forward in time then less iterations should be needed
                niter = ds%tridiagonal_inner_iter
            else
                ! If we could not extrapolate forward in time then more iterations may be needed
                niter = 10_ip*ds%tridiagonal_inner_iter
            end if
        else
            ! If we did not try to extrapolate forward in time, then assume the initial guess was good.
            niter = ds%tridiagonal_inner_iter
        end if

        iter_loop: do iter = 1, niter

#ifdef EVOLVE_TIMER
            !
            !UH and VH update, separated for evolve timing
            !
            if(is_staggered_grid) then
EVOLVE_TIMER_START('disp_sweep_uh')
                call linear_dispersive_sweep_staggered_TRIDIAG(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
                    ds%RHS(:,:,UH), ds%RHS(:,:,VH), &
                    msl_linear, distance_bottom_edge, distance_left_edge, ds%td1, ds%td2, dlon, dlat, &
                    update_UH=.true., update_VH=.false.)
EVOLVE_TIMER_STOP('disp_sweep_uh')
EVOLVE_TIMER_START('disp_sweep_vh')
                call linear_dispersive_sweep_staggered_TRIDIAG(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
                    ds%RHS(:,:,UH), ds%RHS(:,:,VH), &
                    msl_linear, distance_bottom_edge, distance_left_edge, ds%td1, ds%td2, dlon, dlat, &
                    update_UH=.false., update_VH=.true.)
EVOLVE_TIMER_STOP('disp_sweep_vh')
            else
EVOLVE_TIMER_START('disp_sweep_uh')
                call linear_dispersive_sweep_cellcentred_TRIDIAG(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
                    ds%RHS(:,:,UH), ds%RHS(:,:,VH), &
                    msl_linear, distance_bottom_edge, distance_left_edge, ds%td1, ds%td2, dlon, dlat, &
                    update_UH=.true., update_VH=.false.)
EVOLVE_TIMER_STOP('disp_sweep_uh')
EVOLVE_TIMER_START('disp_sweep_vh')
                call linear_dispersive_sweep_cellcentred_TRIDIAG(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
                    ds%RHS(:,:,UH), ds%RHS(:,:,VH), &
                    msl_linear, distance_bottom_edge, distance_left_edge, ds%td1, ds%td2, dlon, dlat, &
                    update_UH=.false., update_VH=.true.)
EVOLVE_TIMER_STOP('disp_sweep_vh')
            end if
#else
            !
            !UH and VH update
            !
            if(is_staggered_grid) then
                call linear_dispersive_sweep_staggered_TRIDIAG(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
                    ds%RHS(:,:,UH), ds%RHS(:,:,VH), &
                    msl_linear, distance_bottom_edge, distance_left_edge, ds%td1, ds%td2, dlon, dlat, &
                    update_UH=.true., update_VH=.true.)
            else
                call linear_dispersive_sweep_cellcentred_TRIDIAG(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
                    ds%RHS(:,:,UH), ds%RHS(:,:,VH), &
                    msl_linear, distance_bottom_edge, distance_left_edge, ds%td1, ds%td2, dlon, dlat, &
                    update_UH=.true., update_VH=.true.)
            end if
#endif

        end do iter_loop

    end subroutine

end module
